// drrateenv.cpp: implementation of the drrateenv class.
//
//////////////////////////////////////////////////////////////////////

#include "drrateenv.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const int CALIB_VOL_MAT = 10;

bool DRRateEnv::operator==(const DRRateEnv& a) const
{
	bool ans = (m_envDate == a.m_envDate) &&  
		(m_baseVols == a.m_baseVols) && (m_zeroCurves == a.m_zeroCurves) &&
		(m_calibParams == a.m_calibParams);

	return ans;
}

DRRateEnv::DRRateEnv (DRDate& envDate, DRZeroCurves& zeroCurves, DRBaseVols& baseVols, 
					  DRSwapVols& swapVols, double meanReversion) 
					  : m_envDate(envDate), m_swapVols(swapVols), m_baseVols(baseVols),
					  m_zeroCurves(zeroCurves)
{
	m_calibParams.resize (7);
	m_calibParams[0] = 6;
	m_calibParams[1] = 1;
	m_calibParams[2] = meanReversion;
	m_calibParams[3] = 5;
	m_calibParams[4] = 0;
	m_calibParams[5] = 1;
	m_calibParams[6] = 0;
}

DRRateEnv::DRRateEnv(MbsIO& mbsIO)
{
	MbsIOMap& mbsIOMap = mbsIO.get_map();
	DRString dt = mbsIOMap.get("envDate");
	m_envDate = toDate(dt);
	
	m_calibParams.resize(7);
	m_calibParams[0] = 6;
	m_calibParams[1] = 1;
	m_calibParams[2] = toNumber(mbsIOMap.get("meanReversion"));
	m_calibParams[3] = 5;
	m_calibParams[4] = 0;
	m_calibParams[5] = 1;
	m_calibParams[6] = 0;

	bool useBaseVols = toBool (mbsIOMap.get("useBaseVols"));
	
	MbsIOMap& envObjectMap = mbsIOMap["ENVOBJECT"].get_map();

	DRString envChoice = envObjectMap.get("EnvChoice");

	if (envChoice == "KAPITAL") {
		m_swapVols = DRSwapVols(KAPITAL);
		m_zeroCurves = DRZeroCurves(KAPITAL, m_envDate);
		
		if (useBaseVols) {
			m_baseVols = DRBaseVols(KAPITAL, m_envDate);
		}
	}
	else if (envChoice == "SWAP") {
		MbsIOMap& fileMap = envObjectMap["ARGS"].get_map();

		DRString path;
		if (fileMap.is_key("path")) path = fileMap.get("path");

		DRString swapVolFile = path + fileMap.get("swapvolfilename", "swo_vol.usd");
		m_swapVols = DRSwapVols(SWAP, swapVolFile);
		
		DRString yieldFile = path + fileMap.get("yieldfilename", "yldcrvlc.usd");
		DRString cmtFile = path + fileMap.get("humpfilename", "cmtspds.usd");
		m_zeroCurves = DRZeroCurves(SWAP, m_envDate, cmtFile, yieldFile);
		
		if (useBaseVols) {
			DRString baseVolFile = path + fileMap.get("basevolfilename", "base_vol.usd");
			m_baseVols = DRBaseVols (SWAP, m_envDate, baseVolFile);
		}		
	}
	else if (envChoice == "DIRECT") {
		MbsIO& args = envObjectMap["ARGS"];

		if (!useBaseVols)
			throw DRException ("Currently, must use Basevols for direct env");

		m_zeroCurves = DRZeroCurves (m_envDate, args);
		if (useBaseVols) {
			m_baseVols = DRBaseVols (m_envDate, args);
		}
	}
	else 
		throw DRException ("Illegal type of EnvChoice ") << envChoice;
	
	if (!useBaseVols) {
		m_baseVols = DRBaseVols(m_envDate, m_swapVols, CALIB_VOL_MAT);
	}
}

ostream& operator<<(ostream& s, DRRateEnv& a)
{
	s << "EnvDate\t" << a.m_envDate << endl;

	s << "SwapVols\n" << a.m_swapVols;
	s << "BaseVols\n" << a.m_baseVols;
	s << "ZeroCurves\n" << a.m_zeroCurves;

	s << "CalibParams:\t" << a.m_calibParams;


	return s;
}
