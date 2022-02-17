//	MlEqVolatilityStructure.cpp : Volatility structure & volatiltiy data
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqVolatilityStructure.h"
#include "MlEqInterpolator.h"
#include "MlEqMaths.h"
#include "MlEqNormal.h"
#include "MlEqObjects.h"
#include "MlEqStrike.h"


/****************************************************************
**	Class  : MlEqVolatilityStructure::Scenario + child classes
**	Routine: Total implementation
**	Returns: n / a
**	Action : n / a
**           
****************************************************************/

MlEqVolatilityStructure::scenario_type MlEqVolatilityStructure::Scenario_FsaLiquidityShift::GetScenarioType(void) const
{
	return MlEqVolatilityStructure::stPostVolShift;
}

MlEqVolatilityStructure::Scenario_FsaLiquidityShift::Scenario_FsaLiquidityShift(double fAverageVolume, double fCoreVolatility, double fNetGamma, double fPositionDirection, bool bIsStock)
{
	m_fAverageVolume = fAverageVolume;
	m_fCoreVolatility = fCoreVolatility;
	m_fNetGamma = fNetGamma;
	m_fPositionDirection = fPositionDirection;
	m_bIsStock = bIsStock;
}

double MlEqVolatilityStructure::Scenario_FsaLiquidityShift::Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const
{
	if (m_fAverageVolume < 1e-7) return f;

	double fMultiplier = !m_fCoreVolatility ? std::max(0.01, 0.07 * f) : m_fCoreVolatility;
	if (m_bIsStock) fMultiplier *= 2.0;
	fMultiplier += 0.01 * f * m_fNetGamma / m_fAverageVolume * sqrt(2.0 / PI) * 20.0;	// Why not, eh?
	f = std::max(0.0001, f - m_fPositionDirection * fMultiplier);
	return f;
}

MlEqVolatilityStructure::scenario_type MlEqVolatilityStructure::Scenario_GetVolType::GetScenarioType(void) const
{
	return MlEqVolatilityStructure::stVolShift;
}

MlEqVolatilityStructure::Scenario_ShiftForwardSkew::Scenario_ShiftForwardSkew(long nWindowStartDate, long nWindowEndDate, long nPeriodStartDays, long nPeriodEndDays, double fFwdSkewShiftAmount)
{			
	m_nWindowStartDate			=	nWindowStartDate;
	m_nWindowEndDate			=	nWindowEndDate;
	m_nPeriodStartDays			=	nPeriodStartDays;
	m_nPeriodEndDays			=	nPeriodEndDays;
	m_fFwdSkewShiftAmount		=	fFwdSkewShiftAmount;
}

double MlEqVolatilityStructure::Scenario_ShiftForwardSkew::Apply(double f, MlEqConstStrikeHandle hStrikeInvariant, DATE dateStart, DATE dateEnd, double fReferenceVolatility) const
{
	if (!((long)dateStart <= m_nWindowEndDate && (long)dateStart >= m_nWindowStartDate)){
		return f;  // dates are outside the window
	}
	long dif = (long)dateEnd - (long)dateStart;
	if (!( dif <= m_nPeriodEndDays && dif >= m_nPeriodStartDays ) ){
		return f;  // dates are outside the window
	}

	double shiftVal = -fReferenceVolatility * m_fFwdSkewShiftAmount * hStrikeInvariant->m_strike;
	f += shiftVal;
	return f;
}

MlEqVolatilityStructure::scenario_type MlEqVolatilityStructure::Scenario_ShiftForwardSkew::GetScenarioType(void) const
{	
	return MlEqVolatilityStructure::stForwardSkewShift;
}

MlEqVolatilityStructure::Scenario_ShiftBeta::Scenario_ShiftBeta(double bShift, estring szTenor)
{	
	CVector vTenor, vBump;

	vTenor.resize(4);
	vBump.resize(4);

	vTenor[0] = 1;	// ugly but simple...
	vTenor[1] = 3;
	vTenor[2] = 6;
	vTenor[3] = 12;

	szTenor.trim();
	szTenor.lc();

	if (szTenor.size()){
		int tenor = 0;

		if (szTenor == "1m"){
			tenor = 0;
		} else if (szTenor == "3m"){
			tenor = 1;
		} else if(szTenor == "6m"){
			tenor = 2;
		} else if (szTenor == "1y" || szTenor == "12m"){
			tenor = 3;
		} else {
			throw "Wrong tenor for beta shift : " + szTenor;
		}
		
		vBump[tenor] = bShift;				
	} else {
		vBump.resize(4, bShift);
	}
	m_hBetaShift = new MlEqInterpolator(vTenor, vBump);
}

double MlEqVolatilityStructure::Scenario_ShiftBeta::Apply(double f, MlEqConstStrikeHandle, DATE dateStart, DATE dateEnd, double) const
{
	double frequencyRU = MlEqStochBetaVolatilityStructure::getFrequencyRoundUp(dateStart, dateEnd);
	double df = m_hBetaShift->getValue(frequencyRU);
	f += df;
	return f;
}

MlEqVolatilityStructure::scenario_type MlEqVolatilityStructure::Scenario_ShiftBeta::GetScenarioType(void) const
{
	return stBetaShift;
}

/*explicit*/ MlEqVolatilityStructure::Scenario_RelativeShift::Scenario_RelativeShift(double fAmount) : m_fAmount(fAmount)
{
}

double MlEqVolatilityStructure::Scenario_RelativeShift::Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const
{
	f *= 1.0 + m_fAmount;
	if (f < 0.0) throw "Negative volatility encountered after shifting";
	return f;
}

/*explicit*/ MlEqVolatilityStructure::Scenario_AbsoluteShift::Scenario_AbsoluteShift(double fAmount) : m_fAmount(fAmount)
{
}

double MlEqVolatilityStructure::Scenario_AbsoluteShift::Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const
{
	f += m_fAmount;
	if (f < 0.0) throw "Negative volatility encountered after shifting";
	return f;
}

/*explicit*/ MlEqVolatilityStructure::Scenario_Set::Scenario_Set(double fAmount) : m_fAmount(fAmount)
{
}

double MlEqVolatilityStructure::Scenario_Set::Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const
{
	f = m_fAmount;
	if (f < 0.0) throw "Negative volatility encountered after shifting";
	return f;
}

bool MlEqVolatilityStructure::Scenario_Set::IsResetting(void) const
{
	return true;
}

/*explicit*/ MlEqVolatilityStructure::Scenario_ShiftSkew::Scenario_ShiftSkew(const std::vector<MlEqStrikeHandle>& Strikes, const CVector& Maturities, const CMatrix& Shift, MlEqConstAssetHandle hAsset)
{	
	m_bRecursionBlock = false;
	
	m_shiftStrikes	= Strikes;
	if (m_shiftStrikes.size() != Shift.cols() || Maturities.getsize() != Shift.rows()){
		throw "Invalid parameters supplied to the ShiftSkew scenario. The number of strikes, maturities and shifts differ";
	}	
	
	//
    // first we set up the one dimensional interpolator for maturities
	// i.e. interpolator across Time (=acrossY)
	//

	int Nmats = Maturities.getsize();
	int Nstrk = m_shiftStrikes.size();

	MlEqInterpolatorHandle MaturityInterpolator = new MlEqInterpolator(Maturities, Maturities);

	CVector xVals(Nstrk);
		
	for ( int i = 0 ; i < Nstrk; i++ )
	{
		xVals[i] = m_shiftStrikes[i]->m_strike;
	}

	CVector yVals(Nstrk);
	GVector<CVector> yValsArray(Nmats);
	
	for ( int i = 0 ; i < Nmats; i++ )
	{
		yValsArray[i] = Shift[i];
	}

    GVector<CVector> xValsASMATRIX(1);
	xValsASMATRIX[0] = xVals;

	m_hAsset = hAsset; 


	MlEqInterpolatorHandle ShiftMatrix = new MlEqInterpolator(xValsASMATRIX, yValsArray);
/*
	int			addTanhWings = 1;
	double		cL = 0.01;	
	double		cR = 3.;		
	int			yPower	= 1;			
	
	MlEqInterpolatorHandle ShiftMatrix = new MlEqCubicSplineInterpolator(xValsASMATRIX, yValsArray, addTanhWings, cL, cR, yPower );
*/
	// this is rather useless
//	GVector<CVector> MATURITIES(1);
//	MATURITIES[0] = Maturities ;
//	MaturityInterpolator = new MlEqCubicSplineInterpolator(MATURITIES, MATURITIES, addTanhWings, cL, cR, yPower );
  
	

	m_ShiftInterpolator = new MlEq2DInterpolator(MaturityInterpolator, ShiftMatrix);		
}


double MlEqVolatilityStructure::Scenario_ShiftSkew::Apply(double f, MlEqConstStrikeHandle hStrike, DATE date, DATE, double) const
{			
	if (m_bRecursionBlock){
		// This can happen if the m_shiftStrikes[0] is a normalised strike since
		// MlEqNormalizedStrike::reinitializeStrike calls GetVol which lands us
		// back here!		
		return f;
	}
		
	double shiftVal, Val;
	if (hStrike->m_strikeType != Normalised || hStrike->GetIsFloating()){
		flip_flop<bool> ff(m_bRecursionBlock, true, false);
		m_shiftStrikes[0]->reinitializeStrike(*m_hAsset);						
		MlEqStrikeHandle tmphStrike = m_shiftStrikes[0]->copy();
		MlEqStrike::convertStrikes(*tmphStrike, *hStrike);	// convert input stike into shift strike type
		Val = tmphStrike->m_strike;
	} else {
		Val = hStrike->m_strike;
	}

	shiftVal = m_ShiftInterpolator->getValue(Val, (long)date);
	f += shiftVal;

	if (f < 0.0) throw "Negative volatility encountered after shifting";
	return f;

}

/*explicit*/ MlEqVolatilityStructure::Scenario_ShiftStrikeSkew::Scenario_ShiftStrikeSkew(MlEqInterpolatorHandle hVolShiftInterpolator, MlEqInterpolatorHandle hVolSkewShiftInterpolator, MlEqConstAssetHandle hAsset, bool bIsRelative, bool bStrikeSkewInvariant)
{
	m_bRecursionBlock = false;	
	m_hShiftStrikeInvariantVolInterpolator = hVolShiftInterpolator;
	m_hShiftStrikeInvariantSkewInterpolator = hVolSkewShiftInterpolator;
	m_hAsset = hAsset;
	m_bStrikeSkewInvariant = bStrikeSkewInvariant;
	m_bRelativeShift = bIsRelative;
}

double MlEqVolatilityStructure::Scenario_ShiftStrikeSkew::Apply(double f, MlEqConstStrikeHandle hStrike, DATE date, DATE, double) const
{
	if (m_bRecursionBlock){
		// This happens when we create a new normalised strike object and call
		// the getNaturalATMVol function.		
		return f;
	}
	
	double shiftVal = 0.0, shiftSkew = 0.0;

	// shiftVal
	if (!!m_hShiftStrikeInvariantVolInterpolator && (!m_bStrikeSkewInvariant || !m_hShiftStrikeInvariantSkewInterpolator)){
		// We just do a vega report in this case (no skews are shifted)
		shiftVal = m_hShiftStrikeInvariantVolInterpolator->getValue(date);
	}
	
	// shiftSkew
	if (!!m_hShiftStrikeInvariantSkewInterpolator){		
		if (hStrike->m_strikeType != Normalised || hStrike->GetIsFloating()){
			flip_flop<bool> ff(m_bRecursionBlock, true, false);
			
			MlEqConstDateHandle	hDate = m_hAsset->GetDateHandle();
			MlEqStrikeHandle pStrike = new MlEqNormalizedStrike(*m_hAsset, 0.0, date, false, hDate);
			MlEqStrike::convertStrikes(*pStrike, *hStrike);	// convert input stike into shift strike type
			shiftSkew = pStrike->m_strike;
		} else {
			shiftSkew = hStrike->m_strike;
		}		
		
		flip_flop<bool> ff(m_bRecursionBlock, true, false);
		double volatm = m_hAsset->getNaturalATMVol(date);
		double fac = 1.0;
		if ( m_bStrikeSkewInvariant && !!m_hShiftStrikeInvariantVolInterpolator){
			fac = m_hShiftStrikeInvariantVolInterpolator->getValue(date); 
		}

		shiftSkew *= -fac * m_hShiftStrikeInvariantSkewInterpolator->getValue(date) * volatm;
	}

	if (!m_bRelativeShift){
		f += shiftVal + shiftSkew;
	} else {
		f *= 1.0 + shiftVal + shiftSkew;
	}	
	if (f < 0.0) throw "Negative volatility encountered after shifting";
	return f;
}


double MlEqVolatilityStructure::Scenarios::Apply(double f, MlEqConstStrikeHandle hStrike, DATE date1, DATE date2, double fReferenceVolatility, scenario_type st) const
{
	std::map<scenario_type, std::vector<RCPtr<Scenario> > >::const_iterator itMap = find(st);
	if (itMap != end()){
		for (std::vector<RCPtr<Scenario> >::const_iterator it = itMap->second.begin(); it != itMap->second.end(); ++it){
			f = (*it)->Apply(f, hStrike, date1, date2, fReferenceVolatility);			
		}
	}
	return f;
}

bool MlEqVolatilityStructure::Scenarios::IsFirstResetting(scenario_type st) const
{
	std::map<scenario_type, std::vector<RCPtr<Scenario> > >::const_iterator itMap = find(st);	
	return itMap != end() && itMap->second.size() && (*itMap->second.begin())->IsResetting();
}

void MlEqVolatilityStructure::Scenarios::push_back(const RCPtr<Scenario>& h)
{				
	scenario_type st = h->GetScenarioType();
	std::map<scenario_type, std::vector<RCPtr<Scenario> > >::iterator itMap = find(st);

	if (itMap == end()){
		std::pair<scenario_type, std::vector<RCPtr<Scenario> > > p(st, std::vector<RCPtr<Scenario> > ());		
		itMap = insert(itMap, p);
	} else if (h->IsResetting()){
		// There is no point in applying any scenarios before a resetting one.		
		itMap->second.clear();		
	}		
	itMap->second.push_back(h);	
}

long MlEqVolatilityStructure::Scenarios::size(MlEqVolatilityStructure::scenario_type st) const
{	
	std::map<scenario_type, std::vector<RCPtr<Scenario> > >::const_iterator itMap = find(st);
	if (itMap == end()) return 0L;
	return itMap->second.size();		
}

/****************************************************************
**	Class  :  MlEqVolatilityStructure 
**	Routines: All scenario perturbations.
**	Returns:  n / a
**	Action :  n / a
**           
****************************************************************/

void MlEqVolatilityStructure::FsaLiquidityShift(double fAverageVolume, double fCoreVolatility, double fNetGamma, double fPositionDirection, bool bIsStock)
{
	RCPtr<Scenario> h;
	h = new Scenario_FsaLiquidityShift(fAverageVolume, fCoreVolatility, fNetGamma, fPositionDirection, bIsStock);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::parallelShiftVol(double VolShift, bool shiftIsRelative)
{			
	RCPtr<Scenario> h;
	if (shiftIsRelative){
		h = new Scenario_RelativeShift(VolShift);
	} else {
		h = new Scenario_AbsoluteShift(VolShift);
	}	
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::Reset()
{
	m_scenarios = m_scenarios_original;
}

void MlEqVolatilityStructure::SetVol(double vol)
{
	RCPtr<Scenario> h = new Scenario_Set(vol);
	m_scenarios.push_back(h);	
}

void MlEqVolatilityStructure::setUpBetaShift(double bshift, const std::string& szTenor)	
{
	throw "Wrong volatility structure type";
}

void MlEqStochBetaVolatilityStructure::setUpBetaShift(double bshift, const std::string& szTenor)
{
	RCPtr<Scenario> h = new Scenario_ShiftBeta(bshift, szTenor);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::setUpFwdSkewShift(double shiftAmount, long WindowStartDate, long WindowEndDate, long PeriodStartDays, long PeriodEndDays)
{
	RCPtr<Scenario> h = new Scenario_ShiftForwardSkew(WindowStartDate, WindowEndDate, PeriodStartDays, PeriodEndDays, shiftAmount);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::shiftSkew(std::vector<MlEqStrikeHandle>& Strikes, CVector& Maturities, CMatrix& Shift, MlEqConstAssetHandle hAsset)
{
	RCPtr<Scenario> h = new Scenario_ShiftSkew(Strikes, Maturities, Shift, hAsset);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::shiftStrikeInvariantSkew(MlEqInterpolatorHandle VolShiftInterpolator, MlEqInterpolatorHandle VolSkewShiftInterpolator, MlEqAssetHandle hAsset, bool shiftIsRelative)
{
	RCPtr<Scenario> h = new Scenario_ShiftStrikeSkew(VolShiftInterpolator, VolSkewShiftInterpolator, hAsset, shiftIsRelative, true);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::shiftStrikeSkew(MlEqInterpolatorHandle VolShiftInterpolator, MlEqInterpolatorHandle VolSkewShiftInterpolator, MlEqAssetHandle hAsset)
{
	RCPtr<Scenario> h = new Scenario_ShiftStrikeSkew(VolShiftInterpolator, VolSkewShiftInterpolator, hAsset, false, false);
	m_scenarios.push_back(h);
}

void MlEqVolatilityStructure::Stick(void)
{	
	m_scenarios_original = m_scenarios;
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: MlEqVolatilityStructure
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

MlEqVolatilityStructure::MlEqVolatilityStructure() :  m_ds(NoDataSource)
{
	// Nothing is given. Bad constructor
	m_spot = 0.0;
	m_asymptodicVolIsInitialized = 0;							
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: MlEqVolatilityStructure
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

MlEqVolatilityStructure::MlEqVolatilityStructure(		
		const double					spot, 
		const CVector&					discountFactors, 
		const CVector&					forwards,
		const vector<DATE>&				terms, 
		MlEqDateHandle					hDate,
		const vector < MlEqVolDataHandle> &		volBidData,
		const vector < MlEqVolDataHandle> &		volMidData,
		const vector < MlEqVolDataHandle> &		volAskData,
		const MlEqInterpolatorHandle&	interpAcrossTime,
		const int interpAcrossTimeFlag,
		const VolatilityDataTypeEnum& voltype,
		MlEqVolatilityStructureHandle pAsymtodicVols,
		double decayFactor) : m_ds(NoDataSource)

{

	initialize(spot,discountFactors,forwards,terms,hDate,
			   volBidData,volMidData,volAskData,interpAcrossTime,
			   interpAcrossTimeFlag,voltype,pAsymtodicVols, decayFactor);

}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: MlEqVolatilityStructure
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

void MlEqVolatilityStructure::initialize( 
		const double					spot, 
		const CVector&					discountFactors, 
		const CVector&					forwards,
		const vector<DATE>&				terms, 
		MlEqDateHandle					hDate,
		const vector < MlEqVolDataHandle> &		volBidData,
		const vector < MlEqVolDataHandle> &		volMidData,
		const vector < MlEqVolDataHandle> &		volAskData,
		const MlEqInterpolatorHandle&	interpAcrossTime,
		const int interpAcrossTimeFlag,
		const VolatilityDataTypeEnum& voltype,
		MlEqVolatilityStructureHandle pAsymptotic,
		double decayFactor
)


{
	m_asymtodicVolStruct = pAsymptotic;
	if (!!pAsymptotic ){
		m_asymptodicVolIsInitialized=1;
	} else{
		m_asymptodicVolIsInitialized = 0;
	}

	
	m_decayFactor = decayFactor;

	m_spot					=  spot;
	m_discountFactors		=  discountFactors;
	m_forwards				=  forwards;
	m_Dates					=  terms;		
	m_hDate					=  hDate;
	m_InterpAcrossTime		=  interpAcrossTime;
	m_interpAcrossTimeFlag	=  interpAcrossTimeFlag;
	m_volType				=  voltype;
	m_BidVolData			=  volBidData;	
	m_MidVolData			=  volMidData;
	m_AskVolData			=  volAskData;
				
	m_scenarios_original = m_scenarios;
	

	if( discountFactors.getsize() != forwards.getsize() )
	{
		throw("What´s going on Paul:size of discount arra");
	}

	int ndates = terms.size();

	if( ndates == 0 )
		throw(  "No terms were entered." );

	int ndataHandles = volMidData.size();
	if ( ndataHandles != 0 )
	{
		if ( ndataHandles == 1)
		{
			if ( volMidData[0]->getNumberOfSlices() != ndates )
			{
				throw(  "number od Dates entered is inconsistent with number of lines encountered in VolData " );
			}
		}
		else if ( ndataHandles != ndates )
		{
			throw(  "number of date handles is inconsistent with number of dates entered" );
		}

	}


	ndataHandles = volBidData.size();
	if ( ndataHandles != 0 )
	{
		if ( ndataHandles == 1)
		{
			if ( volBidData[0]->getNumberOfSlices() != ndates )
			{
				throw(  "inconsistent VolData encountered" );
			}
		}
		else if ( ndataHandles != ndates )
		{
			throw(  "inconsistent VolData encountered" );
		}

	}

	ndataHandles = volAskData.size();
	if ( ndataHandles != 0 )
	{
		if ( ndataHandles == 1)
		{
			if ( volAskData[0]->getNumberOfSlices() != ndates )
			{
				throw(  "inconsistent VolData encountered" );
			}
		}
		else if ( ndataHandles != ndates )
		{
			throw(  "inconsistent VolData encountered" );
		}
	}

}

/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: ~MlEqVolatilityStructure
**	Returns: nothing
**	Action : destructor
**           
****************************************************************/

MlEqVolatilityStructure::~MlEqVolatilityStructure()
{
	m_asymptodicVolIsInitialized = 0;
	Reset();
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: ~MlEqVolatilityStructure
**	Returns: nothing
**	Action : copy operator
**           
****************************************************************/

void MlEqVolatilityStructure::setAsymtodicVolStruct(MlEqVolatilityStructureHandle volstruct)
{
	
	m_asymtodicVolStruct			= volstruct;
	m_asymptodicVolIsInitialized	= !!volstruct ? 1 : 0;
}

/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: ~MlEqVolatilityStructure
**	Returns: nothing
**	Action : copy operator
**           
****************************************************************/

MlEqVolatilityStructure& MlEqVolatilityStructure::operator=( const MlEqVolatilityStructure& rhs )
{
	if (&rhs == this) return *this;

	m_spot					= rhs.m_spot;
	m_discountFactors		= rhs.m_discountFactors;
	m_forwards				= rhs.m_forwards;
	m_Dates					= rhs.m_Dates;
	m_hDate					= rhs.m_hDate;
	m_InterpAcrossTime		= rhs.m_InterpAcrossTime;
	m_volType				=  rhs.m_volType;
	m_BidVolData			=  rhs.m_BidVolData;
	m_MidVolData			=  rhs.m_MidVolData;
	m_AskVolData			=  rhs.m_AskVolData;
    m_interpAcrossTimeFlag	=  rhs.m_interpAcrossTimeFlag;
	m_asymptodicVolIsInitialized = rhs.m_asymptodicVolIsInitialized;
		
	m_scenarios = rhs.m_scenarios;
	m_scenarios_original = rhs.m_scenarios_original;
	throw(" make sure this makes really a proper deep copy; copying all data");
	return *this;
}



/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: MlEqVolatilityStructure
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

MlEqVolatilityStructure::MlEqVolatilityStructure(const MlEqVolatilityStructure& rhs)
:
m_spot( rhs.m_spot ),
m_discountFactors( rhs.m_discountFactors ),
m_forwards( rhs.m_forwards ),
m_Dates( rhs.m_Dates ),
m_hDate( rhs.m_hDate ),
m_InterpAcrossTime ( rhs.m_InterpAcrossTime ),
m_volType( rhs.m_volType),
m_BidVolData(rhs.m_BidVolData),
m_MidVolData(rhs.m_MidVolData),
m_AskVolData(rhs.m_AskVolData),
m_interpAcrossTimeFlag(rhs.m_interpAcrossTimeFlag),
m_asymptodicVolIsInitialized(rhs.m_asymptodicVolIsInitialized),
m_scenarios(rhs.m_scenarios),
m_scenarios_original(rhs.m_scenarios_original)
{	
}

DataSourceEnum MlEqVolatilityStructure::GetDataSource(void) const
{
	return m_ds;
}

MlEqConstDateHandle MlEqVolatilityStructure::getDateToDouble(void) const
{
	return m_hDate;
}

long MlEqVolatilityStructure::GetDateNumber(void) const
{
	return m_hDate->GetDate();
}
	
/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/
	
double MlEqVolatilityStructure::getDensity(const MlEqStrike& strike, DATE date, double forward, BidAskMidEnum bidOrAsk, int cummulativeFlag) 
{	

	double res;	
	double eps = 0.005;
	double xstrike,xstrike0;

	MlEqConstDateHandle d2d			=	getDateToDouble();
	double mat				=	d2d->GetYearFraction(date);

	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);

	xstrike0				=	fixedStrike.m_strike;
	double vol				=	getVol(strike,date,bidOrAsk); 
	if ( vol < 0  ){
		throw("have you set up bid/offer vols properly? missing data encountered");
	}


	res						=   -2.0*Bs(forward,vol,mat,xstrike0,1.0,1);

	if ( cummulativeFlag ){
		res /= 2.0;
	}

	xstrike					=	xstrike0*(1.0+eps);
	fixedStrike.m_strike	=	xstrike;
	vol						=	getVol(fixedStrike,date,bidOrAsk); 
	res						+=  Bs(forward,vol,mat,xstrike,1.0,1);

	if ( !cummulativeFlag )
	{
		xstrike					=	xstrike0*(1.0-eps);
		fixedStrike.m_strike	=	xstrike;
		vol						=	getVol(fixedStrike,date,bidOrAsk); 
		res						+=  Bs(forward,vol,mat,xstrike,1.0,1);

		res /=  pow(fixedStrike.m_strike*eps,2.0);
	}
	else
	{
		res /=  fixedStrike.m_strike*eps;
		res  +=  1.0;
	}

	return res;
	
}	

/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/
	
double MlEqVolatilityStructure::StrikeInTermsOfStandartDeviations(const MlEqStrike& strike,const DATE& date,double forward,BidAskMidEnum bidOrAsk ) 
{	

	double res;	
	double eps = 0.005;
	double xstrike,xstrike0;

	MlEqConstDateHandle d2d			=	getDateToDouble();
	double mat				=	d2d->GetYearFraction(date);

	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);

	xstrike0				=	fixedStrike.m_strike;
	double vol				=	getVol(strike,date,bidOrAsk); 
	if ( vol < 0  ){
		throw("have you set up bid/offer vols properly? missing data encountered");
	}

	res						=   -Bs(forward,vol,mat,xstrike0,1.0,1);
	xstrike					=	xstrike0*(1.0+eps);
	fixedStrike.m_strike	=	xstrike;
	vol						=	getVol(fixedStrike,date,bidOrAsk); 
	res						+=  Bs(forward,vol,mat,xstrike,1.0,1);
	res						/=	fixedStrike.m_strike*eps;
	res						+=  1.0;
	res						=	normal_inv(res);

	return res;
	
}	

	
/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

double MlEqVolatilityStructure::getVol(const MlEqStrike& strike,int islice,BidAskMidEnum bidOrAsk) 
{
	switch(bidOrAsk){
	case Bid:	
		return getVol(strike,islice,m_BidVolData);	
	case Ask:	
		return getVol(strike,islice,m_AskVolData);	
	case Middle:
		return getVol(strike,islice,m_MidVolData);
	}	
	throw "unknown bid/offer type encountered";	
}


// searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
void dateLocate(const vector< DATE >& xx, const DATE& x, int& j)
{
	int ju,jm,jl;
	bool ascnd;

    int n=xx.size();
	jl=-1;						// initialise lower
	ju=n;						// and upper limits
	ascnd=(xx[n-1]>xx[0]);		// true if ascending order of table, false otherwise
	while (ju-jl > 1)			// if we are not yet done
	{
		jm=(ju+jl)>>1;			// compute a midpoint
		if(x>=xx[jm]==ascnd)	
			jl=jm;				// and replace either the lower limit
		else
			ju=jm;				// or upper limit, as appropriate
	}							// repeat until test condition satisfied
	
	if( x == xx[0] ) j = 0;		// set the output
	else if( x == xx[n-1] ) j = n-2;
	else j = jl;
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

void MlEqVolatilityStructure::getYDataAcrossTime(const MlEqStrike& strike,CVector& yData,CVector& xData,int startPoint,int numPoints,vector <MlEqVolDataHandle >	& volData)
{

	// Note that we delegate the handling of negative xData values (corresponding to a volatility surface slice before the surface reference date)
	// to this function.

	// ToDo - Alex to use InterpolateTypeEnum instead
	if ( m_interpAcrossTimeFlag == 0 )
	{
		yData.resize(numPoints);

		// this is just regular volinterpolation
		for(int i = 0; i < numPoints; i++ ){
			yData[i] =  getVol( strike, startPoint + i ,volData);
		}
	}
	else if ( m_interpAcrossTimeFlag == 1 )
	{
		// this corresponds to interpolation in variance
		yData.resize(numPoints);

		double volY,timeX;
		for(int i = 0; i < numPoints; i++ )
		{
			xData[i]	= timeX = m_hDate->GetYearFraction( m_Dates[ startPoint + i ] );

			volY		= getVol( strike, startPoint + i ,volData);

			if ( timeX < 1.0e-4 ){
				yData[i]	= volY;
			}
			else
			{
				yData[i]	= volY*volY*timeX;
			}
		}
	}
	else if ( m_interpAcrossTimeFlag == 2 ){
		// this corresponds to interpolation in the actual volData
		throw("this case is not implemented yet");
	}else
	{
		throw("incorrect m_interpAcrossTimeFlag encountered");
	}

}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getVolFromYData(const MlEqStrike& strike,double xtime,double yData,CVector& xValues)
{

	if ( m_interpAcrossTimeFlag == 0 )
	{
		return yData;
	}
	else if ( m_interpAcrossTimeFlag == 1 )
	{
		// this corresponds to interpolation in variance

		double res;
		if ( xtime < 1e-4 && xValues[0] < 1e-4)
		{
			res = yData;
		}
		else if ( xtime < xValues[0] ){
			res = sqrt(yData/xValues[0]);
		}
		else		{
			res = sqrt(yData/xtime);
		}

		return res;
	}
	else if ( m_interpAcrossTimeFlag == 2 )
	{
		throw("case not implemented yet");
//		return volData.getVol(strike,yParam,xtime) ;
	}
	else
	{
		throw("incorrect m_interpAcrossTimeFlag encountered");
	}

}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : This is the lowest level get volatility function
**           
****************************************************************/

double MlEqVolatilityStructure::getVol(const MlEqStrike& strike, DATE date, std::vector<MlEqVolDataHandle>& volData)
{	
	double								res;			
	MlEqStrikeHandle					hStrike = strike.copy();		// ToDo - remove the deep copy and replace with //MlEqConstStrikeHandle hStrike(&strike);
		
	// If the first element (if any) in m_scenarios is resetting then we don't need to get the
	// 'natural' volatility from the vol data collection (which is expensive).
	if (m_scenarios.IsFirstResetting(stVolShift)){
		res = 0.0;
	} else {
		// Get the 'natural' volatility value from the vol data objects.
		// We assume that the vectors of interpolators are many.
			
		if (!volData.size()) throw "no voldata has been set getVol function";
		
		CVector retVal(1);

		// find where this date is relative to m_date
		// get the m_dates that is of interest based on the lagPoints
		int endPoint = m_Dates.size() - 1;
		int whichPosition = -1;
		dateLocate( m_Dates, date, whichPosition );

		if ( m_interpAcrossTimeFlag == 1 && whichPosition == endPoint ){						
			if (m_Dates.size() > 1){		// To prevent intinite recursion
				res = getVol(strike,m_Dates[endPoint],volData) ;
				return res;
			}
		}

		int lapPoints = m_InterpAcrossTime->lagPoints();
		// do the m_InterpAcrossTime
		int startPoint = 0;

		startPoint = MlEqMaths::Max(whichPosition - lapPoints, startPoint);
		endPoint = MlEqMaths::Min(whichPosition + lapPoints, endPoint);	
		
		int numPoints = endPoint - startPoint + 1;

		CVector xData(numPoints);
		CVector yData(numPoints);

		double timeX;	
		for(int i = 0; i < numPoints; i++ ){
			xData[i] = timeX = m_hDate->GetYearFraction( m_Dates[ startPoint + i ] );			
		}
		getYDataAcrossTime(strike,yData,xData,startPoint,numPoints,volData);

		m_InterpAcrossTime->initialize(xData, yData);	
		double term = m_hDate->GetYearFraction(date);
		
		m_InterpAcrossTime->getValues(retVal, term);
		res = getVolFromYData(strike, term, retVal[0], xData);				
	}
	
	// Apply any scenarios.
	res = m_scenarios.Apply(res, hStrike, date, 0, 0.0, stVolShift);
	
	// Now applpy the post vol shift scenarios.
	res = m_scenarios.Apply(res, hStrike, date, 0, 0.0, stPostVolShift);
	
	return res;
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: 
**	Action : 
**           
****************************************************************/
double MlEqVolatilityStructure::getVol(const MlEqStrike& strike, int islice, std::vector<MlEqVolDataHandle>& VolData)
{
	if (!VolData.size()){
		throw "VolData not set";
	}

	double retVal = 0.0;

	// error checking for islice		
	if (VolData.size() == 1){
		retVal = VolData[0]->getVol( strike, islice);
	} else {
		retVal = VolData[islice]->getVol(strike, islice);
	}

	return retVal;
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/


double MlEqVolatilityStructure::getVol(const MlEqStrike& strike, DATE date, BidAskMidEnum bidOrAsk)
{
	switch(bidOrAsk){
	case Bid:	
		return getVol(strike, date, m_BidVolData);
	case Ask:	
		return getVol(strike, date, m_AskVolData);
	case Middle:
		return getVol(strike,date,m_MidVolData) ;		
	}	
	throw "Unknown bid/offer type encountered";	
}	





















/*
void MlEqVolatilityStructure::getVolDeriv(CVector& res, const MlEqStrike& strike, const DATE& date) 
{

	res.resize(4);
	MlEqVolDataHandle	volData = m_MidVolData[0];
	CVector volderivs(3);	

	int startPoint = 0;
	int endPoint = m_Dates.size() - 1;
	int whichPosition = -1;
	dateLocate( m_Dates, date, whichPosition );

	int nslices = volData->getNumberOfSlices();


	if(whichPosition < startPoint || whichPosition >= endPoint || nslices == 1)
	{
		whichPosition = std::max( startPoint, std::min( endPoint, whichPosition) );
		volData->getVolDeriv( volderivs, strike, whichPosition );

		for(int k=0; k<3; k++){
			res[k] = volderivs[k];
		}
		res[3] = 0.;// in that case we assume no time depedency because we are out of bounds...	
		return;
	}



	CVector volderivs_end(3);	

	startPoint	= MlEqMaths::Max(whichPosition, startPoint );
	endPoint	= MlEqMaths::Min(whichPosition + 1, endPoint );	

	volData->getVolDeriv( volderivs, strike, startPoint );
	volData->getVolDeriv( volderivs_end, strike, endPoint );

	double t  = m_hDate->GetYearFraction( date );
	double ts = m_hDate->GetYearFraction( m_Dates[startPoint] );
	double te = m_hDate->GetYearFraction( m_Dates[endPoint] );

	double wght = (t-ts)/(te-ts);

	if( m_interpAcrossTimeFlag == 0 )	// linear interpolation...
	{
		for(int k=0; k<3; k++)		{
			res[k] = volderivs[k] * (1.-wght) + volderivs_end[k] * wght;
		}

		res[3] = (volderivs_end[2] - volderivs[2]) / (te-ts) ;	// time derivative
	}
	else if( m_interpAcrossTimeFlag == 1 )	// variance interpolation
	{

		double vols = volderivs[2];
		double vole = volderivs_end[2];
		double vars = vols*vols*ts;
		double vare = vole*vole*te;

		double vol = wght * vare + (1.-wght) * vars;
		vol = sqrt( vol/t );

		double dvars = vols * ts * volderivs[0] ;
		double dvare = vole * te * volderivs_end[0] ;
		double dvar = wght * dvare + (1.-wght) * dvars ;

		double d2vars = ts * (volderivs[0]*volderivs[0] + vols * volderivs[1]) ;
		double d2vare = te * (volderivs_end[0]*volderivs_end[0] + vole * volderivs_end[1]) ;
		double d2var = wght * d2vare + (1.-wght) * d2vars ;


		res[0] = dvar / (vol*t) ;
		res[1] = ( d2var / t - res[0]*res[0] ) / vol ;
		res[2] = vol;

		res[3] = 0.5 * ( (vare - vars)/(vol*(te-ts)) - vol ) / t ;

	} else {
		throw "wrong interpolation type";
	}
	
}

*/





void MlEqVolatilityStructure::getVolDeriv(CVector& res, const MlEqStrike& strike, const DATE& date, const CVector& mktForwards) 
{

	res.resize(4);
	MlEqVolDataHandle	volData = m_MidVolData[0];
	CVector volderivs(3);	

	int startPoint = 0;
	int endPoint = m_Dates.size() - 1;
	int whichPosition = -1;
	dateLocate( m_Dates, date, whichPosition );

	int nslices = volData->getNumberOfSlices();


	if(whichPosition < startPoint || whichPosition >= endPoint || nslices == 1)
	{
		whichPosition = std::max( startPoint, std::min( endPoint, whichPosition) );
		volData->getVolDeriv( volderivs, strike, whichPosition, mktForwards[whichPosition] );
//		getVolDeriv( volderivs, strike, whichPosition, mktForwards[whichPosition]);

		for(int k=0; k<3; k++){
			res[k] = volderivs[k];
		}
		res[3] = 0.;// in that case we assume no time depedency because we are out of bounds...	

//		getVolBumpDeriv( res, strike, date);

		return;
	}



	CVector volderivs_end(3);	


	startPoint	= MlEqMaths::Max(whichPosition, startPoint );
	endPoint	= MlEqMaths::Min(whichPosition + 1, endPoint );	

	volData->getVolDeriv( volderivs, strike, startPoint, mktForwards[startPoint]  );
	volData->getVolDeriv( volderivs_end, strike, endPoint, mktForwards[endPoint]  );

//	getVolDeriv( volderivs, strike, startPoint, mktForwards[startPoint]  );
//	getVolDeriv( volderivs_end, strike, endPoint, mktForwards[endPoint]  );


	double t  = m_hDate->GetYearFraction( date );
	double ts = m_hDate->GetYearFraction( m_Dates[startPoint] );
	double te = m_hDate->GetYearFraction( m_Dates[endPoint] );

	double wght = (t-ts)/(te-ts);

	if( m_interpAcrossTimeFlag == 0 )	// linear interpolation...
	{
		for(int k=0; k<3; k++)		{
			res[k] = volderivs[k] * (1.-wght) + volderivs_end[k] * wght;
		}

		res[3] = (volderivs_end[2] - volderivs[2]) / (te-ts) ;	// time derivative
	}
	else if( m_interpAcrossTimeFlag == 1 )	// variance interpolation
	{

		double vols = volderivs[2];
		double vole = volderivs_end[2];
		double vars = vols*vols*ts;
		double vare = vole*vole*te;

		double vol = wght * vare + (1.-wght) * vars;
		vol = sqrt( vol/t );

		double dvars = vols * ts * volderivs[0] ;
		double dvare = vole * te * volderivs_end[0] ;
		double dvar = wght * dvare + (1.-wght) * dvars ;

		double d2vars = ts * (volderivs[0]*volderivs[0] + vols * volderivs[1]) ;
		double d2vare = te * (volderivs_end[0]*volderivs_end[0] + vole * volderivs_end[1]) ;
		double d2var = wght * d2vare + (1.-wght) * d2vars ;


		res[0] = dvar / (vol*t) ;
		res[1] = ( d2var / t - res[0]*res[0] ) / vol ;
		res[2] = vol;

		res[3] = 0.5 * ( (vare - vars)/(vol*(te-ts)) - vol ) / t ;

	} else {
		throw "wrong interpolation type";
	}

//	getVolBumpDeriv( res, strike, date);
	
}


void MlEqVolatilityStructure::getVolDeriv(CVector& res, const MlEqStrike& strike, int slice, double mktForward) 
{
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in local vol calculation";
	}

	const double eps = 1e-6;

	double xstrike = strike.m_strike;

	MlEqStrike fixedAtSlice( xstrike*mktForward );
	DATE sliceDate = m_Dates[slice];

	double vol  = getVol(fixedAtSlice, sliceDate );
	fixedAtSlice.m_strike = (xstrike + eps)*mktForward;	
	double vol_ = getVol(fixedAtSlice, sliceDate );
	fixedAtSlice.m_strike = (xstrike - eps)*mktForward ;	
	double _vol = getVol(fixedAtSlice, sliceDate );

	res[0] = (vol_-_vol) / (2.*eps);
	res[1] = (vol_+_vol-2.*vol) / (eps*eps);
	res[2] = vol;

}


void MlEqVolatilityStructure::getVolBumpDeriv( CVector& res, const MlEqStrike& strike, DATE date ) 
{
	
	MlEqStrikeHandle	hStrike = strike.copy();	// ToDo - remove the deep copy and replace with //MlEqConstStrikeHandle hStrike(&strike);
	
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in local vol calculation";
	}

	const double eps = 1e-6;

	double xstrike = strike.m_strike;

	double tmp = 0.0;;
	if (m_scenarios.IsFirstResetting(stVolShift)){
		tmp = 0.0;
	}

	double vol = res[2];
	double dvdk = res[0];
	double d2vdk2 = res[1];

	tmp = vol;
	
	double bump = m_scenarios.Apply(tmp, hStrike, date, 0, 0.0, stVolShift) - tmp ; // bump

	hStrike->m_strike = xstrike + eps ;
	tmp = vol + dvdk * eps + 0.5*d2vdk2 * eps*eps; 
	double bump_ = m_scenarios.Apply(tmp, hStrike, date, 0, 0.0, stVolShift) - tmp ; // bump

	hStrike->m_strike = xstrike - eps ;
	tmp = vol - dvdk * eps + 0.5*d2vdk2 * eps*eps;
	double _bump = m_scenarios.Apply(tmp, hStrike, date, 0, 0.0, stVolShift) - tmp ; // bump

	double dbumpdk		= (bump_-_bump) / (2.*eps);
	double d2bumpdk2	= (bump_+_bump-2.*bump) / (eps*eps);

	res[0] += dbumpdk ;
	res[1] += d2bumpdk2 ;
	res[2] += bump;

	double dt = 1.0; // in days...

	tmp = vol + res[3] * dt/ 365.0;
	hStrike->m_strike = xstrike;

	double dbumpdt = m_scenarios.Apply(tmp, hStrike, date + dt, 0, 0.0, stVolShift) - tmp  - bump ;

	res[3] += dbumpdt * 365.0/dt ; // change that !!!
}






/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: PutDataSource
**	Returns: nothing
**	Action : Sets the Sirius-style data source of this object.
**           This is, arguably, not a quant function.
**           
****************************************************************/
void MlEqVolatilityStructure::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}



/****************************************************************
**	Class  : MlEqVolData 
**	Routine: getNumberOfSlices
**	Returns: int
**	Action : number of slices
**           
****************************************************************/

int MlEqVolData::getNumberOfSlices() const
{
	return m_numberOfSlices;
}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: getNumberOfSlices
**	Returns: int
**	Action : number of slices
**           
****************************************************************/

int MlEqVolData::getNumberOfStikes(int idate) const
{
	if ( getNumberOfSlices() == 1 ){

		return m_strikes[0].size();
	}
	else if ( getNumberOfSlices() > idate ){

		return m_strikes[idate].size();
	}
	else{
		throw("incorrect vol data set up");
	}

}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
 
MlEqVolData::MlEqVolData(	vector<MlEqInterpolatorHandle>&	InterpAcrossStrike,
							vector< vector< MlEqStrikeHandle > >&	strikes	)

{ 
	initialize(	InterpAcrossStrike,strikes	);

}
/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
/*
void MlEqVolData::getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice) 
{
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in derivatives calculation (local vol)";
	}

	const int nsize = getNumberOfSlices();

	double volderiv,secondderiv;

	if ( nsize == 1){
		whichSlice = 0;	
	}
	
	if( ( nsize != 1 && whichSlice >= nsize )|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	

	
	MlEqStrikeHandle pStrike;
	
	double xVal;

	pStrike = getpStrike(whichSlice,0);
	xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,strike);

	double xstrike = pStrike->m_strike;
	getVolSpaceDerivs(res,*pStrike, whichSlice);

	pStrike->m_strike = xVal;
	
	CVector sres(2); 
	pStrike->deriv(sres);

	MlEqForwardBasedStrike* IPtr = (MlEqForwardBasedStrike*) (&strike);
	double fwd_strike = IPtr->m_fwdparameter;

	sres[0] /= fwd_strike;
	sres[1] /= fwd_strike * fwd_strike ;

	volderiv = res[0]*sres[0];

	secondderiv = res[1]*pow(sres[0],2.0) + res[0]*sres[1];

	res[0] = volderiv;
	res[1] = secondderiv;
	
}
*/

void MlEqVolData::getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice, double mktForward) 
{
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in derivatives calculation (local vol)";
	}

	const int nsize = getNumberOfSlices();

	double volderiv,secondderiv;

	if ( nsize == 1){
		whichSlice = 0;	
	}
	
	if( ( nsize != 1 && whichSlice >= nsize )|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	
	
	MlEqStrikeHandle pStrike;

	MlEqStrike fixedStrikeAtSlice(strike.m_strike * mktForward);
	
	pStrike = getpStrike(whichSlice,0);
	double xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,fixedStrikeAtSlice);

	double xstrike = pStrike->m_strike;
	getVolSpaceDerivs(res,*pStrike, whichSlice);

	CVector sres(2); 
	pStrike->deriv(sres);

	sres[0] *= mktForward;
	sres[1] *= mktForward * mktForward ;

	volderiv	= res[0]*sres[0];
	secondderiv = res[1]*pow(sres[0],2.0) + res[0]*sres[1];

	res[0] = volderiv;
	res[1] = secondderiv;

	pStrike->m_strike = xVal;
	
}

/*
void MlEqVolData::getVolSpaceDerivs(CVector& res, MlEqStrike& strike, int slice) 
{
	const double eps = 1e-6;

	double xstrike = strike.m_strike;

	double vol  = getVol(strike, slice);
	strike.m_strike = xstrike * (1.+eps);	
	double vol_ = getVol(strike, slice);
	strike.m_strike = xstrike * (1.-eps);	
	double _vol = getVol(strike, slice);

	res[0] = (vol_-_vol) / (2.*xstrike*eps);
	res[1] = (vol_+_vol-2.*vol) / (xstrike*xstrike*eps*eps);
	res[2] = vol;

	strike.m_strike = xstrike ;
}
*/


void MlEqVolData::getVolSpaceDerivs(CVector& res, MlEqStrike& strike, int slice) 
{
	const double eps = 1e-6;

	double xstrike = strike.m_strike;

	double vol  = getVol(strike, slice);
	strike.m_strike = xstrike + eps;	
	double vol_ = getVol(strike, slice);
	strike.m_strike = xstrike - eps ;	
	double _vol = getVol(strike, slice);

	res[0] = (vol_-_vol) / (2.*eps);
	res[1] = (vol_+_vol-2.*vol) / (eps*eps);
	res[2] = vol;

	strike.m_strike = xstrike ;
}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
 
void MlEqVolData::initialize(	vector<MlEqInterpolatorHandle>&	InterpAcrossStrike,
							 vector< vector< MlEqStrikeHandle > >&	strikes	)

{ 
	m_InterpAcrossStrike	=	InterpAcrossStrike;
	m_strikes				=	strikes;

	int ndim =  m_InterpAcrossStrike.size();

	if ( ndim == 1 )
		ndim = m_InterpAcrossStrike[0]->getNumberOfSlices();

	m_numberOfSlices	=	ndim;

}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/

void MlEqVolData::reinitializeStrikes(const MlEqAsset& asset)
{

//	vector< vector< MlEqStrikeHandle > >	m_strikes;				//[idate][istrike]

	for ( int idate = 0 ; idate < m_strikes.size(); idate++ )
	{
		for ( int istrike = 0 ; istrike < m_strikes[idate].size(); istrike++ )
		{
			m_strikes[idate][istrike]->reinitializeStrike(asset);
		}
	}
}
/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
 
MlEqVolData::MlEqVolData(	vector< MlEqStrikeHandle  >&	strikes	,int numberSlices)
{ 
	initialize(strikes,numberSlices	);
}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : initializes
**           
****************************************************************/

void MlEqVolData::initialize(vector< MlEqStrikeHandle  >&	strikes,int numberSlices )
{

// only one strike is entered per slice
// this constructor is meant to be used for "analytic" voldata

	m_strikes.resize(strikes.size());

	m_numberOfSlices	=	numberSlices;

	for ( int i = 0 ; i < strikes.size(); i++ )
	{
		m_strikes[i].resize(1);
		m_strikes[i] = strikes;//deep copy? alex
	}

}



/****************************************************************
**	Class  : MlEqVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : initializes
**           
****************************************************************/

void MlEqVolData::initialize(const vector < vector< MlEqStrikeHandle  > >&	strikes,int numberSlices )
{
	m_strikes.resize(strikes.size());
	m_numberOfSlices	=	numberSlices;
	m_strikes	= strikes;
	
}



/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
 
MlEqVolData::MlEqVolData( MlEqInterpolatorHandle&	InterpAcrossStrike,
							 vector< MlEqStrikeHandle > &	strikes	)
{ 
// only one timeslice entered

	m_InterpAcrossStrike.resize(1);
	m_InterpAcrossStrike[0] = InterpAcrossStrike;

	m_strikes.resize(1);

	m_strikes[0] = strikes;// deep copy? alex sos


	int ndim =  m_InterpAcrossStrike.size();

	if ( ndim == 1 )
		ndim = m_InterpAcrossStrike[0]->getNumberOfSlices();

	m_numberOfSlices = ndim;
}





/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : constructor
**           
****************************************************************/
 
MlEqVolData::MlEqVolData(const MlEqVolData& rhs)
:
m_InterpAcrossStrike( rhs.m_InterpAcrossStrike ),
m_strikes( rhs.m_strikes ),
m_numberOfSlices(rhs.m_numberOfSlices)
{ 

}




/****************************************************************
**	Class  : MlEqVolData 
**	Routine: MlEqVolData
**	Returns: nothing
**	Action : assignment operator
**           
****************************************************************/

MlEqVolData& MlEqVolData::operator=(const MlEqVolData& rhs)
{
	if ( this == &rhs)
		return *this;

	m_InterpAcrossStrike	= rhs.m_InterpAcrossStrike;
	m_strikes				= rhs.m_strikes;
	m_numberOfSlices		= rhs.m_numberOfSlices;

	return *this;
}

/****************************************************************
**	Class  : MlEqVolData 
**	Routine: ~MlEqVolData
**	Returns: nothing
**	Action : assignment operator
**           
****************************************************************/
	
	
MlEqVolData::~MlEqVolData()
{	
}	
	

	
/****************************************************************
**	Class  : MlEqVolData 
**	Routine: ~getVol
**	Returns: return volatility
**	Action : calculates volatility
**  
****************************************************************/
	
double MlEqVolData::getVol(const MlEqStrike& strike,int islice) 
{	
	//  what about the crap where it is 0th slice size only??? 

	if( islice < 0 || islice >= getNumberOfSlices() )
	{
		throw("function called for incorrect time slice");
	}
	
	int whichSlice = 0;
	if ( m_InterpAcrossStrike.size() > 1 )
	{
		whichSlice = islice;
	}
	
	
	MlEqStrikeHandle pStrike;	
	double xVal;
	pStrike = getpStrike(islice,0);	
	xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,strike);
	
	MlEqInterpolatorHandle pInterp = getObjectFromVector(m_InterpAcrossStrike,islice);

	double vol = pInterp->getValue(pStrike->m_strike,islice);

//	double vol = m_InterpAcrossStrike[whichSlice]->getValue(pStrike->m_strike,islice);


	pStrike->m_strike = xVal;
	
	return vol;
}	
	
	

/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: initNewSvl
**	Returns: returns nothing
**	Action : initializes new svl parameters from old svl parameters
**  
****************************************************************/ 

// this function initializes new svl parameters from old svl parameters
/*static*/ void MlSvlVolData::initNewSvl(SvlCoeffs& svl, bool bThrow/* =true*/)
{	
	svl.bNewValid = false;
	if (!svl.bOldValid){
		if (bThrow) throw "The internal svi parameterisation is not valid";
		return;
	}
	
	double a	= svl.a;
	double b	= svl.b;
	double rho	= svl.rho;
	double m	= svl.m;
	double sig	= svl.sig;
	double t	= svl.t;
	
	if ( sig < 1e-3 ){
		if (bThrow) throw "svl parameter sig parameter must be greater than zero";
		return;
	}

	svl.ATMVolF = a-b*rho*m+b*sqrt(sig*sig+m*m);
	if ( svl.ATMVolF < 0.0 ){
		if (bThrow) throw "inconsistent ram parameters entered";
		return;
	}
	svl.ATMVolF = sqrt(svl.ATMVolF);

	double skew;
	skew = sqrt(t)/(2.0*svl.ATMVolF)*b*(rho-m/sqrt(sig*sig+m*m));
	svl.Skew = skew;

	svl.CallWing = sqrt(t)/svl.ATMVolF*b*(rho+1.0);
	svl.PutWing  = -sqrt(t)/svl.ATMVolF*b*(rho-1.0);

	double volmin = a+b*sig*sqrt(1.0-rho*rho);

	if ( volmin < 0.0 ){
		if (bThrow) throw "Negative volmin encountered in svl parameter conversion";
		return;
	}
	svl.VolMin   = sqrt(volmin);

	svl.bNewValid = true;

}

/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: initOldSvl
**	Returns: returns nothing
**	Action : initializes old svl parameters from new svl parameters
**  
****************************************************************/ 


/*static*/ void MlSvlVolData::initOldSvl(SvlCoeffs& svl, bool bThrow/* = false*/)
//  this function initializes old svl parameters from new svl parameters
{
	svl.bOldValid = false;
	if (!svl.bNewValid){
		if (bThrow) throw "The svl parameterisation is invalid";
		return;
	}
	
	double cw		= svl.CallWing;
	double pw		= svl.PutWing;
	double volatm	= svl.ATMVolF;
	double skew		= svl.Skew;
	double volmin	= svl.VolMin;
	double t		= svl.t;

	if ( cw* pw < 0.0 ){
		if (bThrow) throw "callwing and putwing parameter must have same sign";
		return;
	}

	double eps = 1e-6;

	if ( pow(cw-pw-4.0*skew,2.0) > pow( cw+pw,2.0) ){
		if (bThrow) throw "inconsistency encountered in svl parameter set up";
		return;
	}

	if ( fabs(cw) < eps && fabs(pw) < eps )
	{
		// this is the flat case

		if ( fabs( svl.ATMVolF-svl.VolMin) > eps ){
			if (bThrow) throw "svl-volmin must be identical to svl-volatm in this case";
			return;
		}

		svl.b	= 0.0;
		svl.a	= pow(svl.ATMVolF,2.0);
		svl.rho = 0.0;
		svl.sig = 0.0;
		svl.m	= 0.0;

		svl.bOldValid = true;
		return;
	}



	if ( t < eps ){
		if (bThrow) throw "no svl conversion can be done for zero svl t-parameter";
		return;
	}

	svl.rho = (cw-pw)/(cw+pw);
	svl.b   = 1.0/2.0*volatm/(sqrt(t)*svl.rho)*(cw-pw);

	if ( fabs(svl.rho) > 1.0 ){
		if (bThrow) throw "svl rho-parameter must be between -1 and +1";
		return;
	}

	if ( fabs(svl.ATMVolF-svl.VolMin) < eps )
	{

		if ( fabs(svl.b) > eps )
		{
			double xrho = 1.0+2.0*svl.ATMVolF/sqrt(t)*svl.Skew/svl.b;
			if ( fabs(xrho-svl.rho) < eps )
			{
				svl.m	= svl.rho;// arbitrary
				svl.a	= pow(svl.ATMVolF,2.0);
				svl.sig = 0.0;
				svl.bOldValid = true;
				return;
			}
		}

		// everything is arbitrary in this case

		svl.b = 0.0;
		svl.rho = 0.0;
		svl.m = 0.0;
		svl.sig = 0.0;

		svl.bOldValid = true;
		return;

	}

	double x = svl.rho-2.0*volatm/sqrt(t)*skew/svl.b;

	if ( fabs(x*x-1.0) < 1e-3 ){
		if (bThrow) throw "inconsistent parameters encountered in svl conversion";
		return;
	}
	x = x/sqrt(1.0-x*x);


	if ( fabs(svl.b) < 1e-3 ){
		if (bThrow) throw "inconsistent parameters encountered in svl conversion";
	}

	svl.sig = (volatm*volatm-volmin*volmin)/svl.b;
	svl.sig /= -svl.rho*x+sqrt(1.0+x*x)-sqrt(1.0-svl.rho*svl.rho);
	svl.m = x*svl.sig;
	svl.a = volmin*volmin-svl.b*sqrt(1.0-svl.rho*svl.rho)*svl.sig;

	svl.bOldValid = true;
}


/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: MlSvlVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 
	
MlSvlVolData::MlSvlVolData(const MlSvlVolData& rhs) :	MlEqVolData( rhs )
{	
	m_numberOfSlices = rhs.m_numberOfSlices;

	m_numSvlParams = rhs.m_numSvlParams;

	MlSvlVolData::m_svl= rhs.MlSvlVolData::m_svl; 
	
};	

/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: MlSvlVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 
	
MlSvlVolData::MlSvlVolData( vector< vector< MlEqStrikeHandle > > &strikes, vector<SvlCoeffs> & svl)	
{	
	initialize(strikes,svl);
	
};	
	
/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: MlSvlVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 
	
void MlSvlVolData::initialize( const vector< vector< MlEqStrikeHandle > > &strikes, vector<SvlCoeffs> & svl)	
{	
	m_numberOfSlices		=	svl.size();
	m_svl.resize(m_numberOfSlices);

	for ( int i = 0 ; i < m_numberOfSlices; i++ )
	{
		if (!svl[i].bOldValid) initOldSvl(svl[i]);
		m_svl[i] = svl[i];
	}

	MlEqVolData::initialize(strikes	,m_numberOfSlices);
};	

	
double getSvlVar(double xstrike,const SvlCoeffs& svl)
{
	if (!svl.bOldValid) throw "The internal svi parameterisation is not valid";
	double vol = svl.a+svl.b*(svl.rho*(xstrike-svl.m)+
				 sqrt(pow(svl.sig,2.0)+pow(xstrike-svl.m,2.0)));

	return vol;
}


double getSvlVol(double xstrike,const SvlCoeffs& svl)
{
	if (!svl.bOldValid) throw "The internal svi parameterisation is not valid";
	double vol = svl.a+svl.b*(svl.rho*(xstrike-svl.m)+
				 sqrt(pow(svl.sig,2.0)+pow(xstrike-svl.m,2.0)));

	if ( fabs(svl.t) < 1e-12 )
	{
		throw("time zero encountered in getSvlVol");
	}

	vol			= sqrt(vol);
	return vol;
}

void getSvlVolDeriv(CVector& res,double xstrike,const SvlCoeffs& svl)
{
	if (!svl.bOldValid) throw "The internal svi parameterisation is not valid";
	double  a = svl.a,
			b = svl.b,
			r = svl.rho,
			ms = xstrike-svl.m,
			s = svl.sig*svl.sig;

	double quad = sqrt( s + ms*ms );

	double var = a + b * ( r * ms +	quad );

	double dvar  = b * ( r + ms/quad );
	double d2var = b * ( 1 - (ms*ms)/(quad*quad) ) / quad ;

	double vol		= sqrt(var);
	double dvol		= dvar/(2.*vol) ;
	double d2vol	= ( 0.5*d2var - 0.25*dvar*dvar/var ) / vol ;

	res[0] = dvol;
	res[1] = d2vol;
	res[2] = vol;
}


/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: getVol
**	Returns: vol
**	Action : 
**  
****************************************************************/

/*
void MlSvlVolData::getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice) 
{	
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in derivatives calculation (local vol)";
	}

	const int nsize = getNumberOfSlices();
	if ( nsize == 1){
		whichSlice = 0;	
	}	
	
	if( ( nsize != 1 && whichSlice >= nsize ) || whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}	
	

	double rstrike = strike.m_strike;
	MlEqForwardBasedStrike* IPtr = (MlEqForwardBasedStrike*) (&strike);
	double strike_fwd = IPtr->m_fwdparameter;

	MlEqStrikeHandle pStrike = getpStrike(whichSlice,0);
	IPtr = (MlEqForwardBasedStrike*)(&*pStrike);
	double vol_fwd = IPtr->m_fwdparameter;
	double xstrike = log( rstrike * strike_fwd/vol_fwd );	// be careful that this assumes log-based strikes 

	getSvlVolDeriv( res, xstrike, m_svl[whichSlice] );

	double volderiv	= res[0] / rstrike;
	double secondderiv = (res[1] - res[0]) / (rstrike*rstrike) ;

	res[0] = volderiv;
	res[1] = secondderiv;
}	
*/

void MlSvlVolData::getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice, double mktForward) 
{	
	if( strike.m_strikeType != ForwardBased ){
		throw "Wrong strike type in derivatives calculation (local vol)";
	}

	const int nsize = getNumberOfSlices();
	if ( nsize == 1){
		whichSlice = 0;	
	}	
	
	if( ( nsize != 1 && whichSlice >= nsize ) || whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}	

	double rstrike = strike.m_strike;
	MlEqStrikeHandle pStrike = getpStrike(whichSlice,0);
	MlEqForwardBasedStrike* IPtr = (MlEqForwardBasedStrike*)(&*pStrike);
	double vol_fwd = IPtr->m_fwdparameter;

	double xstrike = log( rstrike * mktForward/vol_fwd );	// be careful that this assumes log-based strikes 

	getSvlVolDeriv( res, xstrike, m_svl[whichSlice] );

	double volderiv	= res[0] / rstrike;
	double secondderiv = (res[1] - res[0]) / (rstrike*rstrike) ;

	res[0] = volderiv;
	res[1] = secondderiv;
}	

	
/****************************************************************
**	Class  : MlSvlVolData 
**	Routine: getVol
**	Returns: vol
**	Action : 
**  
****************************************************************/

double MlSvlVolData::getVol(const MlEqStrike& strike,int whichSlice) 
{	

	const int nsize = getNumberOfSlices();
	double vol;
	
	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	
	if ( nsize == 1)
		whichSlice = 0;	
	
	MlEqStrikeHandle pStrike;
	
	double xVal;

	pStrike = getpStrike(whichSlice,0);
	
	xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,strike);

	double xstrike = pStrike->m_strike;
	vol = getSvlVol(xstrike,m_svl[whichSlice]);	
	pStrike->m_strike = xVal;//sos do something conceptual here!

	return vol;
}



/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MlEqJumpWingVolData::MlEqJumpWingVolData(const MlEqJumpWingVolData& rhs)
:
MlEqVolData( rhs )
{
	m_jw				= rhs.m_jw;
	m_numberOfSlices	= rhs.m_numberOfSlices;
};



/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MLHullVolData::MLHullVolData(const MLHullVolData& rhs)
:
MlEqVolData( rhs )
{
	m_HullData			= rhs.m_HullData;
	m_numberOfSlices	= rhs.m_numberOfSlices;
}

/****************************************************************
**	Class  : MlRamVolData 
**	Routine: MlRamVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MLRamVolData::MLRamVolData(const MLRamVolData& rhs)
	: MlEqVolData( rhs )
{
	m_RamData			= rhs.m_RamData;
	m_numberOfSlices	= rhs.m_numberOfSlices;
}

/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MlEqJumpWingVolData::MlEqJumpWingVolData(vector< vector < MlEqStrikeHandle > > &strikes,vector<JumpWingCoeffs>& jmp)
{
		initialize(strikes,jmp);
}

/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MLHullVolData::MLHullVolData(vector< vector < MlEqStrikeHandle > > &strikes,vector<HullCoeffs>& hull)
{

	initialize(strikes,hull);

}


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MLRamVolData::MLRamVolData(vector< vector < MlEqStrikeHandle > > &strikes,vector<RamCoeffs>& ram)
{
	initialize(strikes,ram);
};



JumpWingCoeffs::JumpWingCoeffs()
{
	ATMVolF		= -1e99;
	a			= -1e99;
	CallWing	= -1e99;
	PutWing		= -1e99;
	ASym		= -1e99;
	t			= -1e99;	

	skew		= -1e99;
	volMin		= -1e99;
}




bool implied_jmpWing_fcn(double x,void* vp,double* f)
{			
	JumpWingCoeffs* pjw = (JumpWingCoeffs*)vp;

	double b = x;
	double c = pjw->CallWing;
	double p = pjw->PutWing;

	double F = pow(b*p/c,c/(c+p))+b*pow(b*p/c,-p/(p+c));
	F *= 1.0/(1.0+b);
	if ( F < 0.0 ){
		throw("Negative argument encountered in JumpWing conversion");
	}
	F = log(F);

	double r = pow(pjw->volMin/pjw->ATMVolF,2.0);
	if ( r >= 1.0 ){
		throw("inconsistent volmin encountered in jumpwing conversion");
	}
	double a = (r-1.0)/F;
    double skew = 1.0/2.0*a*(c-b*p)/(1.0+b);

//	skew = log(skew*skew);

	*f =  1.0/skew;
	return true;
}	


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::initialize(vector< vector < MlEqStrikeHandle > > &strikes,vector<JumpWingCoeffs>& jmp)
{

	m_numberOfSlices		=	jmp.size();

	if ( strikes.size() == jmp.size() )
	{
		m_strikes				=	strikes;
		}
	else
	{
		if ( strikes.size() != 1 )
		{
			throw("array size of strikeHandles should be either one or same size than input data");
		}

		m_strikes.resize(strikes.size());
		
		for (int  i = 0 ; i < strikes.size(); i++ )
		{
			m_strikes[i] = strikes[i];
		}

	}

	m_jw					=	jmp;

	for ( int i = 0 ; i < jmp.size(); i++ ){
		setSkewStrenght_and_Asymmetry(m_jw[i]);
	}
};


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::setSkewFactor(JumpWingCoeffs& jmp)
{
	double skew;

	if ( jmp.ASym <= 0.0 ){
		throw("asymmetry in jumpwing should be greater than zero");
	}

	skew = 1.0/2.0*jmp.a*(jmp.CallWing-jmp.ASym*jmp.PutWing)/(1.0+jmp.ASym);
	jmp.skew = skew;
}


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::setVolMin(JumpWingCoeffs& jmp)
{
	double zmin;

	zmin = 1.0/(jmp.CallWing+jmp.PutWing)*log(jmp.ASym*jmp.PutWing/jmp.CallWing);

	double volmin = getVol(jmp,zmin) ;
	jmp.volMin = volmin;
}



/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::setSkewStrenght_and_Asymmetry(JumpWingCoeffs& jmp)
{

	//  initilise jmp.a and jmp.Asym from jmp.skew and jmp.volMin if necessary

		if ( jmp.volMin < 0.0 ){
			return;
		}

		if ( fabs(jmp.skew) < 1e-6 )
		{
			if ( jmp.volMin != jmp.ATMVolF ){
				throw("if the skew is zero volMin must be identical to baseVol");
			}

			jmp.a = 0.0;
			jmp.ASym = 1.0; // b is arbitrary in this case

			return;
		}


		double lower_x	= 1e-4;
		double upper_x	= 8.0;
		double accuracy = 1e-4;
		int max_tries	= 50;
		double skew		= jmp.skew;
		double result_b;
		int rootfind_flag = eROOTFIND_BRENT_NOGROW;

		jmp.a		= 1;
		jmp.ASym		= 0.4;
		
//		double objVal = log(skew*skew);

		rootfind_solve( rootfind_flag, implied_jmpWing_fcn,  lower_x,  upper_x,
						accuracy,  max_tries,  1.0/skew, &jmp, &result_b,NULL);
	

		jmp.ASym		= result_b;
//		m_jw.ASym	= result_b;

		JumpWingCoeffs* pjw;
		pjw= &(jmp);

		double b = pjw->ASym;
		double c = pjw->CallWing;
		double p = pjw->PutWing;

		double F = pow(b*p/c,c/(c+p))+b*pow(b*p/c,-p/(p+c));
		F *= 1.0/(1.0+b);
		if ( F < 0.0 ){
			throw("Negative argument encountered in JumpWing conversion");
		}
		F = log(F);

		double r = pow(pjw->volMin/pjw->ATMVolF,2.0);
		if ( r >= 1.0 ){
			throw("inconsistent volmin encountered in jumpwing conversion");
		}

		jmp.a	= (r-1.0)/F;
//		m_jw.a	= jmp.a;

		// check min vol

		double zmin = 1.0/(jmp.CallWing+jmp.PutWing)*log(jmp.ASym*jmp.PutWing/jmp.CallWing);


		double checkVolmin = getVol(jmp,zmin);

		double checkSkew;

		checkSkew = 1.0/2.0*jmp.a*(jmp.CallWing-jmp.ASym*jmp.PutWing)/(1.0+jmp.ASym);


}



/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MLHullVolData::initialize(vector< vector < MlEqStrikeHandle > > &strikes,vector<HullCoeffs>& hull)
{

	m_numberOfSlices		=	hull.size();
	if ( strikes.size() == hull.size() )
	{
		m_strikes				=	strikes;
	}
	else
	{
		if ( strikes.size() != 1 )
		{
			throw("array size of strikeHandles should be either one or same size than input data");
		}

		m_strikes.resize(strikes.size());
		
		for (int  i = 0 ; i < strikes.size(); i++ )
		{
			m_strikes[i] = strikes[i];
		}

	}

	m_HullData					=	hull;

};
/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MLRamVolData::initialize(vector< vector < MlEqStrikeHandle > > &strikes,vector<RamCoeffs>& ram)
{

	m_numberOfSlices	=	ram.size();
	if ( strikes.size() == ram.size() )
	{
		m_strikes				=	strikes;
	}
	else
	{
		if ( strikes.size() != 1 ){
			throw("array size of strikeHandles should be either one or same size than input data");
		}

		m_strikes.resize(strikes.size());
		for (int  i = 0 ; i < strikes.size(); i++ ){
			m_strikes[i] = strikes[i];
		}
	}


	for (int i = 0 ; i < ram.size(); i++ ){
		ram[i].VolLowStrike		=	getRamVol(ram[i].DownCutoff->m_strike, ram[i]);
		ram[i].VolHighStrike	=	getRamVol(ram[i].UpCutoff->m_strike, ram[i]);
		ram[i].Vol90Strike		=	getRamVol(0.9, ram[i]);
		ram[i].Vol110Strike	    =	getRamVol(1.1, ram[i]);
	}

	m_RamData					=	ram;
};


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: MlEqJumpWingVolData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MlEqJumpWingVolData::MlEqJumpWingVolData( vector < MlEqStrikeHandle  > &strikes,JumpWingCoeffs& jmp)
{

//  constructor to use for one time slice only

	initialize(strikes,jmp);

};


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::initialize( vector < MlEqStrikeHandle  > &strikes,JumpWingCoeffs& jmp)
{

//  constructor to use for one time slice only

	m_numberOfSlices		=	1;
	m_strikes.resize(m_numberOfSlices);
	m_strikes[0]				=	strikes;
	m_jw.resize(1);

	m_jw[0]					=	jmp;

};

/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: initialize
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqJumpWingVolData::initialize( vector < MlEqStrikeHandle  > &strikes,vector < JumpWingCoeffs >& jmp)
{

// we need to have onle one strike per timeslice as Jumpwing is analytic
// the one strike is necessary to determine the striketype only

	m_numberOfSlices		=	jmp.size();

	m_strikes.resize(strikes.size());

	for ( int i = 0 ; i < m_numberOfSlices; i++ ){
		m_jw[i]	= jmp[i];
	}

	for (int i = 0 ; i < strikes.size(); i++ )
	{
		m_strikes[i][0] = strikes[i];
	}

}


/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: getVol
**	Returns: vol
**	Action : 
**  
****************************************************************/

double MlEqJumpWingVolData::getVol(const MlEqStrike& strike,int whichSlice) 
{	



	const int nsize = getNumberOfSlices();
	
	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	
	if ( nsize == 1)
		whichSlice = 0;	

	int whichStrikeSlice = 0;
	if ( m_strikes.size() > 1 )
	{
		whichStrikeSlice = whichSlice;
	}
	
	MlEqStrike::convertStrikes(*m_strikes[whichStrikeSlice][0],strike);
	
	double xstrike = m_strikes[whichStrikeSlice][0]->m_strike;
	double vol;

	vol = getVol(m_jw[whichSlice],xstrike) ;

	return vol;
}




/****************************************************************
**	Class  : MlEqJumpWingVolData 
**	Routine: getVol
**	Returns: vol
**	Action : 
**  
****************************************************************/

double MlEqJumpWingVolData::getVol(const JumpWingCoeffs& jmp,double xstrike) 
{	

// it's the users responsibility that xstrike is of the correct strike type: Franzl

  /****************************************************************************
  m_jw_0 = vol0 = volATMF;
  m_jw_1 = a
  m_jw_2 = c	call wing slope coefficient
  m_jw_3 = p      put  wing slope coefficient
  m_jw_4 = b	further asymmetry parameter, default to 1, perhaps not needed...
    
  vol^2 = f^2 = vol0^2 + a * log( (exp(c*z) + b*exp(-p*z)) / (1+b) )

  vol^2 = f^2 = vol0^2 * ( 1 + a * log( (exp(c*z) + b*exp(-p*z)) / (1+b) ) )
  
  where z = NS, so that vol0 = vol(z=0) = volATMF.

  ****************************************************************************/

	double vol;

	vol = 	pow(jmp.ATMVolF,2.0) * ( 1.0 + jmp.a * log( ( exp( jmp.CallWing * xstrike) + jmp.ASym * exp( - jmp.PutWing * xstrike ) ) / (1.0 + jmp.ASym)));
	vol = sqrt(vol);

	return vol;

}

/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

MlARPropVolData::MlARPropVolData(
								double	AtmVolLamda,
								double	AtmVolS0,
								double	AtmVolSInfty,
								double	SkewVolLambda,
								double	SkewVolS0,
								double  SkewVolInfty,	
								double	curvature,
								int		timeCutoff,
								CVector	fwds,
								MlEqDateHandle hDate, 
								vector<long> dates	

								)
								 
{

	initialize(	AtmVolLamda,
				AtmVolS0,
				AtmVolSInfty,
				SkewVolLambda,
				SkewVolS0,
				SkewVolInfty,	
				curvature,
				timeCutoff,
				fwds,
				hDate, 
				dates);	

}

/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

void MlARPropVolData::initialize(double	AtmVolLamda,
								double	AtmVolS0,
								double	AtmVolSInfty,
								double	SkewVolLambda,
								double	SkewVolS0,
								double  SkewVolInfty,	
								double	curvature,
								int		timeCutoff,
								CVector	fwds,
								MlEqDateHandle hDate, 
								vector<long> dates	
								)
{

		m_AtmVolLamda		=	AtmVolLamda;
		m_AtmVolS0			=	AtmVolS0;
		m_AtmVolSInfty		=	AtmVolSInfty;
		m_SkewVolLambda		=	SkewVolLambda;
		m_SkewVolS0			=	SkewVolS0;
		m_SkewVolInfty		=	SkewVolInfty;
		m_curvature			=	curvature;
		m_timeCutoff		=	timeCutoff/365.25;
		m_fwds				=	fwds;
		m_numberOfSlices	=	m_fwds.getsize();

		if ( m_fwds.getsize() != dates.size() ){
			throw("inconsitent input arraysizes encountered in MlARPropVoldata");
		}

		m_t.resize(dates.size());
		m_VolAtm.resize(m_fwds.getsize());

		for ( int i = 0 ; i < m_t.getsize(); i++ )
		{
			m_t[i] = hDate->GetYearFraction(dates[i]);
			m_VolAtm[i] = getAtmVol(m_t[i]);
		}

}

/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

MlARPropVolData::MlARPropVolData(const MlARPropVolData& rhs)
{

		m_AtmVolLamda		=	rhs.m_AtmVolLamda;
		m_AtmVolS0			=	rhs.m_AtmVolS0;
		m_AtmVolSInfty		=	rhs.m_AtmVolSInfty;
		m_SkewVolLambda		=	rhs.m_SkewVolLambda;
		m_SkewVolS0			=	rhs.m_SkewVolS0;
		m_SkewVolInfty		=	rhs.m_SkewVolInfty;
		m_curvature			=	rhs.m_curvature;
		m_numberOfSlices	=	rhs.m_numberOfSlices;
		m_timeCutoff		=	rhs.m_timeCutoff;
		m_fwds				=	rhs.m_fwds;	
		m_t					=	rhs.m_t	;
		m_VolAtm			=   rhs.m_VolAtm;
}

/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

double MlARPropVolData::S(double T)
{
	double x	= T/m_SkewVolLambda;
	double pT	= (1.0-exp(-x))/x;
	double s	= pT*pow(m_SkewVolS0,2.0) + (1.0-pT)*pow(m_SkewVolInfty,2.0);
		   s	= sqrt(s);
	return s;
}


/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/


double MlARPropVolData::getAtmVol(double T)
{

	double x	= (T+m_timeCutoff)/m_AtmVolLamda;
	double pT	= (1.0-exp(-x))/x;
	double s	= pT*pow(m_AtmVolS0,2.0) + (1.0-pT)*pow(m_AtmVolSInfty,2.0);
		   s	= sqrt(s);
	return s;
}

/****************************************************************
**	class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

double MlARPropVolData::getVol(const MlEqStrike& strike,int islice) 
{

	double fwd,vol,t;
	if ( m_fwds.getsize() > islice ){
		fwd = m_fwds[islice];
	}
	else{
		throw("incorrect arraysize encountered in MLARPropVolData");
	};

	if ( m_t.getsize() > islice ){
		t = m_t[islice];
	}
	else{
		throw("incorrect arraysize encountered in MLARPropVolData");
	};


	MlEqStrikeHandle pStrike;
	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);
	double xstrike;

	if ( fixedStrike.m_strike > fwd ){
		pStrike = new MlEqLogStrike(1.0,fwd,-123);
	}else{
		pStrike = new MlEqForwardBasedStrike(1.0,fwd,-123);
	}

	MlEqStrike::convertStrikes(*pStrike,strike);

	if ( fixedStrike.m_strike > fwd ){
		xstrike = pStrike->m_strike ;
	}else{
		xstrike = pStrike->m_strike-1.0;
	}

	double s  = S(t+m_timeCutoff);
	vol		  = m_VolAtm[islice] + 100.0*s*xstrike +
				  10.0*400.0/(365.25*pow(m_VolAtm[islice],2.0))*m_curvature/(t+m_timeCutoff)*pow(xstrike,2.0);

	if ( fixedStrike.m_strike > fwd ){
		vol     += 50.0*s*pow(xstrike,2.0);
	}

	return vol;

}
/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

MlARPropFit::MlARPropFit(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CVector& fittingFwds,MlEqDateHandle hDate,CMatrix& BoundsOnConstraints,int timeCutoff,double curvature)
{
	initialize(mktVols,FittingStrikes,fittingDates,fittingFwds,hDate,BoundsOnConstraints,timeCutoff,curvature);
};

/****************************************************************
**	Class  : MlARPropVolData 
**	Routine: constructor
**	Returns: 
**	Action : 
**  
****************************************************************/

void MlARPropFit::initialize(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CVector& fittingFwds,MlEqDateHandle hDate,CMatrix& BoundsOnConstraints,int timeCutoff,double curvature)
{

	m_mktVols	= mktVols;
	m_t.resize( fittingDates.size());
	for ( int i = 0 ; i < fittingDates.size(); i++ )
	{
		m_t[i] = hDate->GetYearFraction(fittingDates[i]);
	}
	
//	_asm int 3;
	m_fwds				=	fittingFwds;
	m_timeCutoff		=	timeCutoff/365.25;
	m_FittingStrikes	=   FittingStrikes;
	m_curvature			=	curvature;
	m_numberOfSlices	=	m_fwds.getsize();

	int ndates			=	m_t.getsize() ;
	
	m_VolAtmMkt.resize(ndates);
	for ( int i = 0 ; i < ndates; i++ )
	{		
		MlEqStrike strike(m_fwds[i]);
		m_VolAtmMkt[i] = mktVols->getVol(strike,(DATE)fittingDates[i]) ;	
	}	

	if ( FittingStrikes.size() != ndates ){
		throw("error in input format for fitting strikes");
	}

	m_VolSkewMkt.resize(ndates);
	double skew;
	for (int i = 0 ; i < ndates; i++ )
	{		
		if ( FittingStrikes[i].size() != 2 ){
			throw("error in format for fitting strikes");
		}
		skew = mktVols->getVol(*FittingStrikes[i][1],(DATE)fittingDates[i])-m_mktVols->getVol(*FittingStrikes[i][0],(DATE)fittingDates[i]) ;
		m_VolSkewMkt[i] = skew;
	}

	iVector linearVariableIndex(2);
	CMatrix NonZeroPartialDerivativesSpecification;
	int outputFlag = 0;
	double InitialTolerance		=	1e-6;
	double FinalTolerance		=	1e-4;	
	double StoppingTolerance	=	1e-16;

	int maximizeFlag = 0;
	int returnCode;

//  start fitting to atm vols first

	m_fitSkew = 0;

	int numVariables = 3;
	CVector initialGuess(numVariables);

	initialGuess[0] = 0.5;
	initialGuess[1] = m_VolAtmMkt[0];
	initialGuess[2] = m_VolAtmMkt[ndates-1];

	CMatrix xValBounds(numVariables,2);

	for (int i = 0 ; i < numVariables; i++ )
	{
		xValBounds[0][0]	=	 0.01;
		xValBounds[0][1]	=	 5;
		xValBounds[1][0]	=	 0.01;
		xValBounds[1][1]	=	 1.5;
		xValBounds[2][0]	=	 0.01;
		xValBounds[2][1]	=	 1.5;
	}

	LSGRGSolver::initialize(initialGuess,linearVariableIndex,
						   xValBounds,BoundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);


	solve( initialGuess,maximizeFlag,returnCode);

	m_VolAtm.resize(m_t.getsize());
	for (int i = 0 ; i < m_VolAtm.getsize(); i++ ){
		m_VolAtm[i] = getAtmVol(m_t[i]);
	}


	m_fitSkew = 1;
	InitialTolerance	=	1e-12;
	FinalTolerance		=	1e-12;	
	StoppingTolerance	=	1e-20;

	initialGuess[0] = 0.5;
	initialGuess[1] = m_VolSkewMkt[0];
	initialGuess[2] = m_VolSkewMkt[ndates-1];

	for (int i = 0 ; i < numVariables; i++ )
	{
		xValBounds[0][0]	=	 0.01;
		xValBounds[0][1]	=	 5;
		xValBounds[1][0]	=	 10.0*(initialGuess[1]-0.02);
		xValBounds[1][1]	=	 10.0*(initialGuess[1]+0.02);
		xValBounds[2][0]	=	 10.0*(initialGuess[2]-0.02);
		xValBounds[2][1]	=	 10.0*(initialGuess[2]+0.02);
	}
	
	LSGRGSolver::initialize(initialGuess,linearVariableIndex,
						   xValBounds,BoundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);


	solve( initialGuess,maximizeFlag,returnCode);

}

/****************************************************************
**	Class  : MlARPropFit 
**	Routine: Objective function
**	Returns: returns double
**	Action : 
**  
****************************************************************/
	
void MlARPropFit::ObjectiveFcn(double* gVals,double* xVals)
{	


	int idate;
	double objVal=0.0,diff,fitvol,fitskew,vfit=0.0;	

	if ( !m_fitSkew )
	{
		m_AtmVolLamda	=	xVals[0];
		m_AtmVolS0		=	xVals[1];
		m_AtmVolSInfty	=	xVals[2];

		for ( idate = 0 ;idate < m_t.getsize(); idate++ )
		{	
			fitvol  =	getAtmVol(m_t[idate]);
			diff	=	fitvol-m_VolAtmMkt[idate];
			vfit	+=	diff*diff;
		}
	}
	else
	{
			m_SkewVolLambda	=	xVals[0];
			m_SkewVolS0		=	xVals[1];
			m_SkewVolInfty	=	xVals[2];

			for ( idate = 0 ;idate < m_t.getsize(); idate++ )
			{	
				fitskew		=	getVol(*m_FittingStrikes[idate][1],idate)-getVol(*m_FittingStrikes[idate][0],idate) ;
				diff		=	fitskew-m_VolSkewMkt[idate];
				vfit		+=	diff*diff;
			}
	}
	
	*gVals	=	vfit;
}

/****************************************************************
**	Class  : MlEqVolMultiplierData 
**	Routine: MlEqVolMultiplierData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

MlEqVolMultiplierData::MlEqVolMultiplierData(const MlEqVolMultiplierData& rhs)
:
MlEqVolData( rhs )
{
	m_referenceVol		= rhs.m_referenceVol;
	m_numberOfSlices	= rhs.m_numberOfSlices;
};

/****************************************************************
**	Class  : MlEqVolMultiplierData 
**	Routine: MlEqVolMultiplierData
**	Returns: returns nothing
**	Action : constructor
**  
****************************************************************/ 

void MlEqVolMultiplierData::PutReferenceVol(CVector& refVols)
{
	int nsize = getNumberOfSlices();

	if ( nsize != refVols.getsize() ){
		throw("incorrect array size in PutReferenceVol encountered");
	}

	m_referenceVol = refVols;
}


/****************************************************************
**	Class  : MlEqVolMultiplierData 
**	Routine: getReferenceVol
**	Returns: returns getReferenceVol of a given slice
**	Action : 
**  
****************************************************************/

double MlEqVolMultiplierData::getReferenceVol(int whichSlice) const
{
	double vol;
	const int nsize = getNumberOfSlices();

	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}

	if ( nsize == 1) {
		vol = m_referenceVol[0];	
	} else {
		vol =  m_referenceVol[whichSlice];
	}

	return vol;
}


/****************************************************************
**	Class  : MlEqVolMultiplierData 
**	Routine: getVol
**	Returns: vol
**	Action : 
**  
****************************************************************/

double MlEqVolMultiplierData::getVol(const MlEqStrike& strike,int whichSlice) 
{
	const double refVol			= getReferenceVol(whichSlice); 
	const int nsize = getNumberOfSlices();

	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 )	{
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}


	if ( nsize == 1)
		whichSlice = 0;	


	MlEqStrike::convertStrikes(*m_strikes[whichSlice][0],strike);

	double vol = InterpAcrossStrikeGetValue(m_strikes[whichSlice][0]->m_strike,whichSlice);

	vol *= refVol;


	return vol;
}

/****************************************************************
**	Class  : MlEqVolData 
**	Routine: InterpAcrossStrikeGetValue
**	Returns: double
**	Action : 
**  
****************************************************************/

double MlEqVolData::InterpAcrossStrikeGetValue(double xval,int whichSlice) const
{
	const int nslices = getNumberOfSlices();

	if( ( nslices != 1 && whichSlice >= nslices )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}

	if ( m_InterpAcrossStrike.size() > 1 )
		return m_InterpAcrossStrike[whichSlice]->getValue(xval,0);
	else
		return m_InterpAcrossStrike[0]->getValue(xval,whichSlice);

}


/****************************************************************
**	Class  : MlEqVolData 
**	Routine: getpStrike
**	Returns: returns strike for given idate,istrike
**	Action : 
**  
****************************************************************/

const MlEqStrikeHandle& MlEqVolData::getpStrike(int idate,int istrike)
{
	return getObjectFromVectorVector(m_strikes,idate,istrike,getNumberOfSlices());

/*	int nsize = m_strikes.size();
	if ( nsize == 1 )
	{
		return m_strikes[0][istrike];
	}
	else if ( nsize == getNumberOfSlices() )
	{
		return m_strikes[idate][istrike];
	}
	else
	{
		throw("indexing error in strike array");
	}
*/
}




/****************************************************************
**	Class  : SvlFit 
**	Routine: getIdeal
**	Returns: returns ideal fitting function
**	Action : 
**  
****************************************************************/


double SvlFit::getIdeal(int ifitParam,const double time,const double parmconst,const double parmlong,const double parmlambda)
{
	double timefactor=1.0;
	if ( ifitParam == 1 || ifitParam == 3 ){
		timefactor = 0.0;
	}

	double idealparam = pow(time,timefactor)*(parmconst+parmlong*(1-((1-exp(-parmlambda*time))/(parmlambda*time))));
	return idealparam;
}

/****************************************************************
**	Class  : SvlFit 
**	Routine: Objective function
**	Returns: returns double
**	Action : 
**  
****************************************************************/
	
void SvlFit::ObjectiveFcn(double* gVals,double* xVals)
{	
	
	int idate,istrike;
	double time,d1,vega,sig,vfit,diff,vol,fitvol;
	
	int	numSlices = getNumberOfSlices();
	MlEqForwardBasedStrike strike;

	double objVal=0.0;

	for ( idate = 0 ; idate < numSlices; idate++ )
	{	
		time = m_doubleTime[idate];

		m_svl[idate].a		= xVals[0+m_numSvlParams*idate];
		m_svl[idate].b		= xVals[1+m_numSvlParams*idate];
		m_svl[idate].sig	= xVals[2+m_numSvlParams*idate];
		m_svl[idate].rho	= xVals[3+m_numSvlParams*idate];
		m_svl[idate].m		= xVals[4+m_numSvlParams*idate];
		m_svl[idate].t		= time;
	
		for ( istrike = 0 ;istrike < m_midVols[idate].getsize(); istrike++ )
		{	
			vfit = 0.0;
			MlEqStrike::convertStrikes(strike,*getpStrike(idate,istrike));

			if ( m_midVols.size() > idate ){
				vol = m_midVols[idate][istrike];
			}
			else if ( m_bidVols.size() > idate ){
				vol = m_bidVols[idate][istrike];
			}
			else if ( m_askVols.size() > idate ){
				vol = m_askVols[idate][istrike];
			}
			else{
				throw("inconsistent input data encountered in svl volfitter");
			}

			sig		=	vol*sqrt(time);
			d1		=	(log(1.0/strike.m_strike)+0.5*sig*sig)/sig;
			vega	=	exp(-0.5*d1*d1);
			fitvol  =	MlEqMaths::Max(getVol(*getpStrike(idate,istrike),idate),0.0);

			if ( m_midVols.size() > idate )
			{
				diff	=	fitvol-m_midVols[idate][istrike];
				vfit	+=	m_fitMidVols*diff*diff;
			}

			if ( m_bidVols.size() > idate )
			{
				diff	=	MlEqMaths::Max(m_bidVols[idate][istrike]-fitvol,0.0);//use these to force inside bid/ask -- usually makes no difference
				vfit	+=	m_fitSpread*diff*diff;
			}

			if ( m_askVols.size() > idate )
			{
				diff	=	MlEqMaths::Max(fitvol-m_askVols[idate][istrike],0.0);
				vfit	+=	m_fitSpread*diff*diff;
			}

			vfit *= vega*vega;
			objVal += vfit;
		}		
	}
	
	if ( m_fitIdeals )
	{					
			double idealparam,parmdiff,parmfactor,x,y,z;
			int offset = m_numSvlParams*numSlices;
			for ( idate = 0 ; idate < numSlices; idate++ )
			{	
				time = m_doubleTime[idate];
				for ( int iparam = 0; iparam < m_numSvlParams; iparam++ )
				{
					if ( iparam == 1 ){
						parmfactor = 10.0;
					}else{
						parmfactor = 1.0;
					}

					x = xVals[offset+0+3*iparam];
					y = xVals[offset+1+3*iparam];
					z = xVals[offset+2+3*iparam];

					idealparam	=	getIdeal(iparam,time,z,x,y);			
				    parmdiff	=	parmfactor*(xVals[idate*m_numSvlParams+iparam]-idealparam);
					parmdiff	*=	parmdiff;

					if ( iparam < 2 ){
						parmdiff /= time;
					}

					objVal		+=	parmdiff*m_fitIdeals;
				}
			}
	}	
			

	*gVals	=	objVal;

}


/****************************************************************
**	Class  : SvlFit 
**	Routine: Objective function
**	Returns: nothing
**	Action : relates "independent" variables to each other
**  
****************************************************************/

void SvlFit::ConstraintFcn(double* contraintVals, double* xVals)
{

	int N = 6;//number of constraints per slice

//  first timeslice has different constraints

	contraintVals[0] = xVals[2];// sigma bigger than zero
	contraintVals[1] = 0.0;//means no constraint
	contraintVals[2] = xVals[1]*(1.0+xVals[3]);// b(1+rho)>0 : total variance for infinite strikes is positive
	contraintVals[3] = xVals[1]*(1.0-xVals[3]);// b(1-rho)>0 : total variance for zero strikes is positive
	contraintVals[4] = 0.0;// means no constraint
	contraintVals[5] = xVals[0] + xVals[1]*xVals[2]*sqrt(1.0-(xVals[3]*xVals[3]));// minimal total variance ( constraint is very important )

	int	numSlices = getNumberOfSlices();

	for ( int i = 1 ; i < numSlices; i++ )
	{
		contraintVals[N*i+0]	= xVals[m_numSvlParams*i+2] - xVals[m_numSvlParams*(i-1)+2]; //dt(sig) > 0
		contraintVals[N*i+1]	= -(xVals[m_numSvlParams*i+4] - xVals[m_numSvlParams*(i-1)+4])* xVals[m_numSvlParams*i+3];//		// sign(dtm)=-sign(rho)
		contraintVals[N*i+2]	= xVals[m_numSvlParams*i+1] * (1 + xVals[m_numSvlParams*i+3]) - xVals[m_numSvlParams*(i-1)+1] * (1 + xVals[m_numSvlParams*(i-1)+3]);		//  d(t) [b*(1+rho)] > 0
		contraintVals[N*i+3]	= xVals[m_numSvlParams*i+1] * (1 - xVals[m_numSvlParams*i+3]) - xVals[m_numSvlParams*(i-1)+1] * (1 - xVals[m_numSvlParams*(i-1)+3]);		//  d(t) [b*(1-rho)] > 0
		contraintVals[N*i+4]	= ((xVals[m_numSvlParams*i+0] - xVals[m_numSvlParams*(i-1)+0]) - (xVals[m_numSvlParams*(i-1)+1] * (1 + xVals[m_numSvlParams*(i-1)+3])*(xVals[m_numSvlParams*i+4] - xVals[m_numSvlParams*(i-1)+4])));  // dt(a) - [b*(1+rho) * dt(m)] >0
		contraintVals[N*i+5]	= ((xVals[m_numSvlParams*i+0] - xVals[m_numSvlParams*(i-1)+0]) + (xVals[m_numSvlParams*(i-1)+1] * (1 - xVals[m_numSvlParams*(i-1)+3])*(xVals[m_numSvlParams*i+4] - xVals[m_numSvlParams*(i-1)+4])));  // dt(a) - [b*(1+rho) * dt(m)] >0
	}

}	


/****************************************************************
**	Class  : SvlFit 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/

void SvlFit::initialize(vector<CVector>& bidVols,vector<CVector>& midVols,vector<CVector>& askVols,vector < vector< MlEqStrikeHandle  > >&	strikes	,CVector& doubleTimes,double fitIdeals,double fitMidPoints,double fitSpread,CMatrix& initial_Guess,CMatrix& BoundsOnConstraints)
{

	if ( midVols.size() > 0 && askVols.size() > 0 && bidVols.size() > 0 )
	{
		if ( (midVols.size() != askVols.size()) || (midVols.size() != bidVols.size()) ){
			throw("inconststent bid/offer/mid vols entered");
		}
		m_numberOfSlices	=	midVols.size();
	}
	else if ( midVols.size() == 0 && askVols.size() > 0 && bidVols.size() > 0 )
	{
		if ( askVols.size() != bidVols.size() ){
			throw("inconststent bid/offer vols entered");
		}
		m_numberOfSlices	=	askVols.size();
	}
	else if ( midVols.size() > 0 && askVols.size() == 0 && bidVols.size() == 0 ){
			m_numberOfSlices	=	midVols.size();
	}

// other cases can easily be implemented

	m_bidVols = bidVols;
	m_midVols = midVols;
	m_askVols = askVols;
	m_fitIdeals = fitIdeals; 

	MlEqVolData::initialize(strikes,m_numberOfSlices);
	MlSvlVolData::initialize();

	int ii;
	if ( fitIdeals ){
		ii = 1;
	}else{
		ii = 0;
	}

	int numFitParams= m_numberOfSlices*m_numSvlParams+ii*3*m_numSvlParams;
	CVector initGuess,initGuess2;
	initGuess.resize(numFitParams);

	if ( initial_Guess.rows() > 0 )
	{
		int num= m_numberOfSlices*m_numSvlParams;

		if (!( (initial_Guess.rows() != num) || (initial_Guess.cols() != num)) ){
			throw("incorrect number of guesses entered");
		}

		int j, k = 0 ;
		for ( int i = 0 ; i < m_numberOfSlices; i++ )
		{
			for ( j = 0 ; j < m_numSvlParams; j++ ){
				initGuess[k++] = initial_Guess[i][j];
			}
		}

		if ( m_fitIdeals )
		{
			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = -1.0;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;
		}

	}
	else
	{
		initGuess2.resize(initGuess.getsize());

		int k = 0 ;
		int l = 0 ;
		for ( int i = 0 ; i < m_numberOfSlices; i++ )
		{
			initGuess[k++] = 0.1;
			initGuess[k++] = 0.1;
			initGuess[k++] = 1.0;
			initGuess[k++] = -0.5;
			initGuess[k++] = 0.1;
		}

		// input ideal guesses
		if ( m_fitIdeals )
		{
			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = -1.0;

			initGuess[k++] = 1.0;
			initGuess[k++] = 1.0;
			initGuess[k++] = 0.25;
		}

// second guess:

		for (int i = 0 ; i < m_numberOfSlices; i++ )
		{
			double	time	= doubleTimes[i];
			double	averageVar = 0.0;
			int count = 0;
			if ( midVols.size() > 0 )
			{
				if ( midVols.size() != m_numberOfSlices ){
					throw("not enough midVol slices have been supplied");
				}

				for ( int j = 0 ; j < midVols[i].getsize(); j++ ){
					averageVar += pow(midVols[i][j],2.0);
					count++;
				}
				averageVar /= count;
			}
			else
			{
				if ( (bidVols.size() != m_numberOfSlices) || (askVols.size() != m_numberOfSlices )){
					throw("not enough bidVol/askVol slices have been supplied");
				}

				for ( int j = 0 ; j < bidVols[i].getsize(); j++ ){
					averageVar += pow(0.5*(bidVols[i][j]+askVols[i][j]),2.0);
					count++;
				}
				averageVar /= count;
			}

			double	rho		= -0.5;
			double 	theta	= 5.0 / time * (1 - (1 - exp(-2*time))/(2*time));
			double 	a		= averageVar * (1 - rho * rho) / 2.0;
			double	b		= averageVar * theta / 2.0;
			double	m		= -rho/theta;
			double	sigma	= sqrt(1 - rho * rho) / theta;

			initGuess2[l++] = a;
			initGuess2[l++] = b;
			initGuess2[l++] = sigma;
			initGuess2[l++] = rho;
			initGuess2[l++] = m;
		}
		// input ideal guesses
		if ( m_fitIdeals )
		{
			initGuess[l++] = 1.0;
			initGuess[l++] = 1.0;
			initGuess[l++] = 0.25;

			initGuess[l++] = 1.0;
			initGuess[l++] = 1.0;
			initGuess[l++] = 0.25;

			initGuess[l++] = 1.0;
			initGuess[l++] = 1.0;
			initGuess[l++] = 0.25;

			initGuess[l++] = 1.0;
			initGuess[l++] = 1.0;
			initGuess[l++] = -1.0;

			initGuess[l++] = 1.0;
			initGuess[l++] = 1.0;
			initGuess[l++] = 0.25;
		}

	}

//	do some more work here!

	initialize(initGuess,strikes,doubleTimes,fitIdeals,fitMidPoints,fitSpread,BoundsOnConstraints,initGuess2);

}

/****************************************************************
**	Class  : SvlFit 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/

void SvlFit::initialize(CVector& initialGuess,vector < vector< MlEqStrikeHandle  > >&	strikes	,CVector& doubleTimes,double fitIdeals,double fitMidPoints,double fitSpread,CMatrix& BoundsOnConstraints,CVector& initialGuess2)
{	
	
// needs to be completed

	m_fitIdeals		=	fitIdeals;
	m_fitMidVols	=	fitMidPoints;
	m_fitSpread		=	fitSpread;
	m_doubleTime	=	doubleTimes;

	
// define bounds on constraints
	
	int i;
	int N = 6;
	double big = 1000;
	int	numSlices = getNumberOfSlices();
	int numConstraints = N*numSlices;	

	m_svl.resize(numSlices);
	m_totalNumVariables = numSlices*m_numSvlParams;
	if ( m_fitIdeals ){
		m_totalNumVariables += m_numSvlParams*3;
	}
	
	if ( BoundsOnConstraints.rows() == 0 )
	{
		BoundsOnConstraints.resize(numConstraints,2);// [iconstraint][2]	
		for ( i = 0 ; i < numConstraints; i++ )
		{
			BoundsOnConstraints[i][0] = -1000000*big;
			BoundsOnConstraints[i][1] =  1000000*big;
		}
//		for ( i = 0 ; i < numSlices; i++ ){
//			BoundsOnConstraints[i*N+(N-1)][0] = 0.0;// this is the contraint that the total variance must be greater than zero	
//		}
	}

// define bounds on variables
	
	CMatrix xValBounds(m_totalNumVariables,2);
	big=5.0;
	for ( i = 0 ; i < numSlices; i++ )
	{
		xValBounds[i*m_numSvlParams+0][0]	=	-big;
		xValBounds[i*m_numSvlParams+0][1]	=	 big;

		xValBounds[i*m_numSvlParams+1][0]	=	 0.00001;
		xValBounds[i*m_numSvlParams+1][1]	=	 big;

		xValBounds[i*m_numSvlParams+2][0]	=	 0.0001;
		xValBounds[i*m_numSvlParams+2][1]	=	 big;

		xValBounds[i*m_numSvlParams+3][0]	=	-1.0;
		xValBounds[i*m_numSvlParams+3][1]	=	 1.0;

		xValBounds[i*m_numSvlParams+4][0]	=	-big;
		xValBounds[i*m_numSvlParams+4][1]	=	 big;

	}

	int offset = numSlices*m_numSvlParams;
	if ( m_fitIdeals )
	{
		for ( i = offset; i < m_totalNumVariables; i++ )
		{
			xValBounds[i][0]	=	-1000.0*big;
			xValBounds[i][1]	=	 1000.0*big;
		}
	}

	iVector linearVariableIndex(2);
	CMatrix NonZeroPartialDerivativesSpecification;
	int outputFlag = 0;

	double InitialTolerance		=	1e-8;
	double FinalTolerance		=	1e-8;	
	double StoppingTolerance	=	1e-16;

	int maximizeFlag = 0;
	int returnCode;

	if ( initialGuess2.getsize() > 0 )
	{

		CMatrix boundsOnConstraints(numConstraints,2);// [iconstraint][2]	
		for ( i = 0 ; i < numConstraints; i++ )
		{
			boundsOnConstraints[i][0] = -1000000*big;
			boundsOnConstraints[i][1] =  1000000*big;
		}
//		for ( i = 0 ; i < numSlices; i++ ){
//			boundsOnConstraints[i*N+(N-1)][0] = 0.0;// this is the contraint that the total variance must be greater than zero	
//		}

		double fitIdeals	= m_fitIdeals;
		m_fitIdeals			= 0;

		double var1,var2;
//		do two unconstraint fits first as set up for main fitter

		LSGRGSolver::initialize(initialGuess,linearVariableIndex,
						   xValBounds,boundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);

		solve( initialGuess,maximizeFlag,returnCode);

		ObjectiveFcn(&var1,initialGuess.getPtr());

		LSGRGSolver::initialize(initialGuess2,linearVariableIndex,
						   xValBounds,boundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);

		solve( initialGuess2,maximizeFlag,returnCode);

		ObjectiveFcn(&var2,initialGuess2.getPtr());

// choose the better of the two

		if ( var2 < var1 ){
			initialGuess = initialGuess2;
		}

		m_fitIdeals			= fitIdeals;

	}

	LSGRGSolver::initialize(initialGuess,linearVariableIndex,
					   xValBounds,BoundsOnConstraints,
					   InitialTolerance,FinalTolerance,
					   StoppingTolerance,
					   NonZeroPartialDerivativesSpecification,outputFlag);



	solve( initialGuess,maximizeFlag,returnCode);

//	MlSvlVolData::initialize(strikes,m_svl)	;

	double xstrike = 0.0;//log(0.5);
	CVector params(6);
	for (  i = 0 ; i < 5;i++ ){
		params[i] = initialGuess[i];
	}
	params[i] = m_svl[0].t;

	SvlCoeffs xsvl(params);

	double test = getSvlVol(xstrike,xsvl);

}
		
/****************************************************************
**	Class  : SvlCoeffs 
**	Routine: SvlCoeffs
**	Returns: nothing
**	Action : constructor
**  
****************************************************************/

SvlCoeffs::SvlCoeffs() : bNewValid(false), bOldValid(false)
{
}
	
SvlCoeffs::SvlCoeffs(CVector& x)
{

	if (  x.getsize() != 6 && x.getsize() != 5  )
	{
		throw("incorrect number of parameters passed to initialize Svl parameters");
	}

	a	=	x[0];
	b	=	x[1];
	sig	=	x[2];
	rho	=	x[3];
	m	=	x[4];
	bNewValid = false;
	bOldValid = true;

	if ( x.getsize() == 6 ){
		t	=	x[5];
	}
	else{
		t	= 0.0;
	}

}

/****************************************************************
**	Class  : 
**	Routine: getHermitePrice
**	Returns: optionPrices
**	Action : 
**  
****************************************************************/

double HermiteCoeff::getHermiteCoeff(int n)
{
	if ( n >= hCoeff.getsize() ){
		return 0.0;
	}
	else{
		double x = hCoeff[n];
		return x;
	}
}

/****************************************************************
**	Class  : 
**	Routine: getHermitePrice
**	Returns: optionPrices
**	Action : 
**  
****************************************************************/


void getHermitePrice(CVector& callPrice,GVector< MlEqStrikeHandle >& strikes, const HermiteCoeff& sc)
{

	int nstrikes = strikes.getsize();
	CVector fixedStrikes(nstrikes);
	callPrice.resize(nstrikes);

	for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
	{	
		MlEqStrike xStrike;
		MlEqStrike::convertStrikes(xStrike,*strikes[istrike]);
		fixedStrikes[istrike] = xStrike.m_strike;
	}
	
	int returnVolFlag = 0;
	HermiteOptions(callPrice,fixedStrikes,sc.hCoeff,sc.mat,sc.fwd,returnVolFlag);

}

/****************************************************************
**	Class  : 
**	Routine: getHermitePrice
**	Returns: optionPrices
**	Action : 
**  
****************************************************************/

void getHermiteVol(CVector& vols, GVector< MlEqStrikeHandle >& strikes, const HermiteCoeff& sc)
{

	int nstrikes = strikes.getsize();
	CVector fixedStrikes(nstrikes);
	vols.resize(nstrikes);

	for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
	{	
		MlEqStrike xStrike;
		MlEqStrike::convertStrikes(xStrike,*strikes[istrike]);
		fixedStrikes[istrike] = xStrike.m_strike;
	}
	
	int returnVolFlag = 1;
	HermiteOptions(vols,fixedStrikes,sc.hCoeff,sc.mat,sc.fwd,returnVolFlag);

}


/****************************************************************
**	Class  : 
**	Routine: getHermitePrice
**	Returns: optionPrices
**	Action : 
**  
****************************************************************/

double  getHermiteVol(const MlEqStrike& strike, const HermiteCoeff& sc)
{
	MlEqStrike xStrike;
	MlEqStrike::convertStrikes(xStrike,strike);	
	int returnVolFlag = 1;
	CVector vols(1);

	CVector stk(1);
	stk[0] = xStrike.m_strike;
	HermiteOptions(vols,stk,sc.hCoeff,sc.mat,sc.fwd,returnVolFlag);

	double res = vols[0];
	return res;
}


/****************************************************************
**	Class  : MLHermiteVolData
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/

void MLHermiteVolData::initialize(  vector<HermiteCoeff> & hpol)
{


	m_numberOfSlices		=	hpol.size();
	m_HermiteCoeff.resize(m_numberOfSlices);

	for ( int i = 0 ; i < m_numberOfSlices; i++ )
	{
		m_HermiteCoeff[i] = hpol[i];
	}

//	MlEqVolData::initialize(strikes	,m_numberOfSlices);

}


/****************************************************************
**	Class  : MLHermiteVolData
**	Routine: MLHermiteVolData
**	Returns: nothing
**	Action : 
**  
****************************************************************/

MLHermiteVolData::MLHermiteVolData( vector<HermiteCoeff>& hpol)
{
	initialize(  hpol);
}



/****************************************************************
**	Class  : MLHermiteVolData
**	Routine: getHermiteData
**	Returns: MLHermiteVolData
**	Action : 
**  
****************************************************************/

const HermiteCoeff& MLHermiteVolData::getHermiteData(int whichSlice)const
{
	if ( whichSlice < 0  || whichSlice >= m_numberOfSlices ){
		throw("incorrect number of slices entered");
	}

	return m_HermiteCoeff[whichSlice];
}


/****************************************************************
**	Class  : MLHermiteVolData
**	Routine: getHermiteData
**	Returns: MLHermiteVolData
**	Action : 
**  
****************************************************************/


double MLHermiteVolData::getVol(const MlEqStrike& strike,int islice) 
{
	double res =  getHermiteVol(strike, m_HermiteCoeff[islice]);
	return res;
}



/****************************************************************
**	Class  : RamFit 
**	Routine: RamFit
**	Returns: nothing
**	Action : 
**  
****************************************************************/

RamFit::RamFit(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int fitBoundaries,int useVega,double fitSpread,double fitMidVols,vector < vector < MlEqStrikeHandle > > & upCutoff,vector < vector <  MlEqStrikeHandle > > & downCutoff)

{
	initialize(mktVols,FittingStrikes	,fittingDates, BoundsOnConstraints,fitBoundaries,useVega,fitSpread,fitMidVols,upCutoff,downCutoff);
}

/****************************************************************
**	Class  : RamFit 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/

void RamFit::initialize(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int fitBoundaries,int useVega,double fitSpread,double fitMidVols,vector < vector < MlEqStrikeHandle > > & upCutoff,vector < vector < MlEqStrikeHandle > > & downCutoff)
{	
	
	MlEqVolData::initialize(FittingStrikes	,FittingStrikes.size());
	
	m_mktVols		=	mktVols;
	m_fitBoundaries	=	fitBoundaries,
	m_useVega		=	useVega;
	m_fitSpread		=	fitSpread;
	m_fitMidVols	=	fitMidVols;

//  define bounds on constraints

	MlEqConstDateHandle dateToDouble =  m_mktVols->getDateToDouble();
	
	int i;
	double big = 1000;
	int	numSlices = getNumberOfSlices();
	int numConstraints = 0;
	
	m_RamData.resize(numSlices);

	int N = 3;
	int numVariables = N;
	if ( m_fitBoundaries ){
		numVariables += 2;
	}

	m_fittingDimension = numVariables;

	if (  downCutoff.size() == 0 && upCutoff.size() != 0 ){
		throw("ramfitting:either downCutoff AND upCutoffs are provided or none of them");
	}

	if ( ( downCutoff.size() == 0 || upCutoff.size() == 0 ) && numVariables == 3 ){
		throw("cutoffs are not provided in Ram vol parametrisation");
	}

	
	CMatrix xValBounds(numVariables,2);

	for ( i = 0 ; i < N; i++ )
	{
		xValBounds[0][0]	=	 0.0;
		xValBounds[0][1]	=	 1.5;
		xValBounds[1][0]	=	 -5;
		xValBounds[1][1]	=	 5;
		xValBounds[2][0]	=	 0;
		xValBounds[2][1]	=	 1;

	}

	for ( i = N+1; i < numVariables; i++ )
	{

// do not know how to do this sos

		xValBounds[3][0]	=	 0;
		xValBounds[3][1]	=	 1;

		xValBounds[4][0]	=	 1;
		xValBounds[4][1]	=	 3;
	}


	CVector initialGuess(numVariables);
	MlEqSpotBasedStrike strikeATMSpot(1.0,1.0,NULL);


	iVector linearVariableIndex(2);
	CMatrix NonZeroPartialDerivativesSpecification;
	int outputFlag = 0;

	double InitialTolerance		=	1e-6;
	double FinalTolerance		=	1e-4;	
	double StoppingTolerance	=	1e-16;


	int maximizeFlag = 0;
	int returnCode;

	MlEqStrikeHandle pStrike,pUpperStrike,pLowerStrike;


	for ( int idate = 0 ; idate < fittingDates.size(); idate++ )
	{

		initialGuess[0] = 	m_mktVols->getVol(strikeATMSpot,(double)fittingDates[idate]) ;

		if ( numVariables == 3){

			pUpperStrike	= getObjectFromVectorVector(upCutoff,idate,0,fittingDates.size());
			pLowerStrike	= getObjectFromVectorVector(downCutoff,idate,0,fittingDates.size());;

			if ( pUpperStrike->m_strike <= pLowerStrike->m_strike ){
				throw("upper strike limit is entered that is lower than lower strike limit");
			}

			initialGuess[1] = 	(m_mktVols->getVol(*pUpperStrike,(double)fittingDates[idate])-m_mktVols->getVol(*pLowerStrike,(double)fittingDates[idate]) )/(pUpperStrike->m_strike-pLowerStrike->m_strike);
			initialGuess[2] = 	(m_mktVols->getVol(*pUpperStrike,(double)fittingDates[idate])+m_mktVols->getVol(*pLowerStrike,(double)fittingDates[idate])-2.0*m_mktVols->getVol(strikeATMSpot,(double)fittingDates[idate]) )/pow(pUpperStrike->m_strike-pLowerStrike->m_strike,2.0);
		}
		else{

			MlEqSpotBasedStrike upperStrike(1.4,1.0,NULL);
			MlEqSpotBasedStrike lowerStrike(0.6,1.0,NULL);
			
			initialGuess[1] = 	(m_mktVols->getVol(upperStrike,(double)fittingDates[idate])-m_mktVols->getVol(lowerStrike,(double)fittingDates[idate]) )/(upperStrike.m_strike-lowerStrike.m_strike);
			initialGuess[2] = 	(m_mktVols->getVol(upperStrike,(double)fittingDates[idate])+m_mktVols->getVol(lowerStrike,(double)fittingDates[idate])-2.0*m_mktVols->getVol(strikeATMSpot,(double)fittingDates[idate]) )/pow(upperStrike.m_strike-lowerStrike.m_strike,2.0);

			initialGuess[3] =   lowerStrike.m_strike;
			initialGuess[4] =   upperStrike.m_strike;

		}

		pStrike = getpStrike(idate,0)->copy();

		m_RamData[idate].t				= dateToDouble->GetYearFraction(fittingDates[idate]);
		m_RamData[idate].currentDate	= fittingDates[idate];
		m_RamData[idate].m_strike		= pStrike;

		if ( numVariables ==  3 ){

			m_RamData[idate].DownCutoff		= pLowerStrike;
			m_RamData[idate].UpCutoff		= pUpperStrike;
			m_RamData[idate].VolLowStrike	= m_mktVols->getVol(*pLowerStrike,(double)fittingDates[idate]);
			m_RamData[idate].VolHighStrike	= m_mktVols->getVol(*pUpperStrike,(double)fittingDates[idate]);
		
		} else {
			m_RamData[idate].DownCutoff		= new MlEqSpotBasedStrike(1.4,1.0,NULL);
			m_RamData[idate].UpCutoff		= new MlEqSpotBasedStrike(0.6,1.0,NULL);
		}
		
		m_RamData[idate].Vol90Strike = m_mktVols->getVol(0.9, (double)fittingDates[idate]);
		m_RamData[idate].Vol110Strike = m_mktVols->getVol(1.1, (double)fittingDates[idate]);



		m_currentidate		= idate;

		LSGRGSolver::initialize(initialGuess,linearVariableIndex,
						   xValBounds,BoundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);

		solve( initialGuess,maximizeFlag,returnCode);

//		if ( returnCode ){
//			throw("solver did not coverge in ram vol fitting function");
//		}do this later sos
	}

}
		
/****************************************************************
**	Class  : RamFit 
**	Routine: Objective function
**	Returns: returns double
**	Action : 
**  
****************************************************************/
	
void RamFit::ObjectiveFcn(double* gVals,double* xVals)
{	

	int istrike;
	double time,d1,vega,sig,diff,vol,volMid,volAsk,volBid,fitvol;
	

	int idate = m_currentidate;
	time = 	m_RamData[idate].t;

	DATE currentDate = m_RamData[idate].currentDate;

	m_RamData[idate].BaseVol		= xVals[0];
	m_RamData[idate].Slope			= xVals[1];
	m_RamData[idate].Curvature		= xVals[2];

	if ( m_fittingDimension > 3 ){

		m_RamData[idate].DownCutoff->m_strike		= xVals[3];
		m_RamData[idate].UpCutoff->m_strike			= xVals[4];

		m_RamData[idate].VolLowStrike	=  	getVol(*m_RamData[idate].DownCutoff,idate) ;
		m_RamData[idate].VolHighStrike	=  	getVol(*m_RamData[idate].UpCutoff,idate) ;
	}

	int numberStrikesMid;  //,numberStrikesBid,numberStrikesAsk;
	m_mktVols->getMidVolData(numberStrikesMid,0,0);

	MlEqStrikeHandle strike;
	double objVal=0.0,vfit;	

	for ( istrike = 0 ;istrike < numberStrikesMid; istrike++ )
	{	

			vfit = 0.0;
			strike = getpStrike(idate,istrike);

			MlEqStrike::convertStrikes(*m_RamData[idate].m_strike,*strike);

			if ( !m_fitBoundaries )
			{
				if (  m_RamData[idate].m_strike->m_strike < m_RamData[idate].DownCutoff->m_strike || 
					  m_RamData[idate].m_strike->m_strike > m_RamData[idate].UpCutoff->m_strike ){
					  continue;
				}
			}

			if ( m_fitSpread )
			{
				vol			= m_mktVols->getVol(*strike,currentDate, Ask);
				volAsk		= vol;

				if ( volAsk < 0.0 && m_fitSpread ){
					throw("no ask vols have been provided");
				}

				vol			= m_mktVols->getVol(*strike,currentDate, Bid);
				volBid		= vol;

				if ( volBid < 0.0 && m_fitSpread ){
					throw("no bid vols have been provided");
				}
			}

			if ( m_fitMidVols )
			{
				vol			= m_mktVols->getVol(*strike,currentDate, Middle);
				volMid		=	vol;

			}
	
			if ( m_useVega )
			{
				sig		=	vol*sqrt(time);
				d1		=	(log(1.0/strike->m_strike)+0.5*sig*sig)/sig;
				vega	=	exp(-0.5*d1*d1);
			}

			fitvol  =	getVol(*strike,idate);

			if ( m_fitMidVols )
			{
				diff	=	fitvol-volMid;
				vfit	+=	m_fitMidVols*diff*diff;
			}

			if ( m_fitSpread )
			{
				diff	=	MlEqMaths::Max(volBid-fitvol,0.0);//use these to force inside bid/ask -- usually makes no difference
				vfit	+=	m_fitSpread*diff*diff;
			}

			if ( m_fitSpread )
			{
				diff	=	MlEqMaths::Max(fitvol-volAsk,0.0);
				vfit	+=	m_fitSpread*diff*diff;
			}

			if ( m_useVega ){
				vfit *= vega*vega;
			}

			objVal += vfit;
		}
	

	*gVals	=	objVal;

}

/****************************************************************
**	Class  : MlEqJumpWingFitVolData 
**	Routine: MlEqJumpWingFitVolData
**	Returns: nothing
**	Action : 
**  
****************************************************************/

MlEqJumpWingFitVolData::MlEqJumpWingFitVolData(vector<JumpWingCoeffs>& JWGuess,MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols)
{
	initialize(JWGuess,mktVols,FittingStrikes,fittingDates,BoundsOnConstraints,useVega,fitSpread,fitMidVols);
}

/****************************************************************
**	Class  : MlEqJumpWingFitVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/


void MlEqJumpWingFitVolData::initialize(MlEqJumpWingVolDataHandle guess,MlEqVolatilityStructureHandle& mktVols,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols)
{
	
	this->MlEqJumpWingVolData::MlEqJumpWingVolData(*guess);

//	initialize(guess.getJumpWingParameters(),mktVols,guess.getStrikes(),fittingDates,BoundsOnConstraints,useVega,fitSpread,fitMidVols,false);
	initialize(m_jw,mktVols,m_strikes,fittingDates,BoundsOnConstraints,useVega,fitSpread,fitMidVols,false);
}

/****************************************************************
**	Class  : MlEqJumpWingFitVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/


MlEqJumpWingFitVolData::MlEqJumpWingFitVolData(MlEqJumpWingVolDataHandle guess,MlEqVolatilityStructureHandle& mktVols,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols)
{
	initialize(guess,mktVols,fittingDates,BoundsOnConstraints,useVega,fitSpread,fitMidVols);
}

/****************************************************************
**	Class  : MlEqJumpWingFitVolData 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**  
****************************************************************/

void MlEqJumpWingFitVolData::initialize(vector<JumpWingCoeffs>& JWGuess,MlEqVolatilityStructureHandle& mktVols,const vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols,bool initializeBaseClass)
{	

	if ( initializeBaseClass ){
		MlEqVolData::initialize(FittingStrikes	,FittingStrikes.size());
	}

	
	m_mktVols		=	mktVols;
	m_useVega		=	useVega;
	m_fitSpread		=	fitSpread;
	m_fitMidVols	=	fitMidVols;

//  define bounds on constraints

	MlEqConstDateHandle dateToDouble =  m_mktVols->getDateToDouble();
	
	int i;
	double big = 1000;
	int	numSlices = getNumberOfSlices();
	int numConstraints = 0;

	if ( JWGuess.size() != numSlices ){
		throw("number of slices incorrect in initial guess for jumpwing parameters");
	}

	m_jw.resize(JWGuess.size());
	for ( i = 0 ; i < JWGuess.size(); i++ ){
		m_jw[i]	= JWGuess[i];
	}

//	JumpWingCoeffs.resize(numSlices);
	int numVariables = 5;
	m_fittingDimension = numVariables;

	
	CMatrix xValBounds(numVariables,2);

	xValBounds[0][0]	=	 0.0;
	xValBounds[0][1]	=	 1.5;
	xValBounds[1][0]	=	 -1.0;
	xValBounds[1][1]	=	 1.0;
	xValBounds[2][0]	=	 0;
	xValBounds[2][1]	=	 3;
	xValBounds[3][0]	=	 0;
	xValBounds[3][1]	=	 3;
	xValBounds[4][0]	=	 -2;
	xValBounds[4][1]	=	 2;



	CVector initialGuess(numVariables);
	MlEqSpotBasedStrike strikeATMSpot(1.0,1.0,NULL);

	iVector linearVariableIndex(2);
	CMatrix NonZeroPartialDerivativesSpecification;
	int outputFlag = 0;

	double InitialTolerance		=	1e-6;
	double FinalTolerance		=	1e-4;	
	double StoppingTolerance	=	1e-8;


	int maximizeFlag = 0;
	int returnCode;

	MlEqStrikeHandle pStrike,pUpperStrike,pLowerStrike;


	for ( int idate = 0 ; idate < fittingDates.size(); idate++ )
	{

		initialGuess[0] = JWGuess[idate].ATMVolF;     //m_mktVols->getVol(strikeATMSpot,(double)fittingDates[idate]) ;
		initialGuess[1] = JWGuess[idate].a;
		initialGuess[2] = JWGuess[idate].CallWing;
		initialGuess[3] = JWGuess[idate].PutWing;
		initialGuess[4] = JWGuess[idate].ASym;

		m_jw[idate].currentDate = fittingDates[idate];
		m_jw[idate].t			= dateToDouble->GetYearFraction(fittingDates[idate]);


		pStrike = getpStrike(idate,0)->copy();
		m_currentidate		= idate;

		LSGRGSolver::initialize(initialGuess,linearVariableIndex,
						   xValBounds,BoundsOnConstraints,
						   InitialTolerance,FinalTolerance,
						   StoppingTolerance,
						   NonZeroPartialDerivativesSpecification,outputFlag);


		solve( initialGuess,maximizeFlag,returnCode);

		setVolMin(m_jw[idate]);
		setSkewFactor(m_jw[idate]);


		JWGuess[idate] = m_jw[idate];


/*
//		test some results here
		for ( int i = 0 ; i < m_strikes[idate].size(); i++ ){
			
			double fitVol = getVol(*m_strikes[idate][i],idate);
			double mktVol = m_mktVols->getVol(*m_strikes[idate][i],idate) ;
			double diff = fitVol - mktVol;
		}	
*/


//		if ( returnCode ){
//			throw("solver did not coverge in ram vol fitting function");
//		}do this later sos
	}

}

/****************************************************************
**	Class  : MlEqJumpWingFitVolData 
**	Routine: Objective function
**	Returns: returns double
**	Action : 
**  
****************************************************************/

void MlEqJumpWingFitVolData::ObjectiveFcn(double* gVals,double* xVals)
{	

	int istrike;
	double time,d1,vega,sig,diff,vol,volMid,volAsk,volBid,fitvol;
	

	int idate	=	m_currentidate;
	time		= 	m_jw[idate].t;

	DATE currentDate = m_jw[idate].currentDate;

	m_jw[idate].ATMVolF			= xVals[0];
	m_jw[idate].a				= xVals[1];
	m_jw[idate].CallWing		= xVals[2];
	m_jw[idate].PutWing			= xVals[3];
	m_jw[idate].ASym			= xVals[4];


//	int numberStrikesMid;
//	m_mktVols->getMidVolData(numberStrikesMid,0,0);

	MlEqStrikeHandle strike;
	double objVal=0.0,vfit;	

	for ( istrike = 0 ;istrike < m_strikes[idate].size(); istrike++ )
	{	

			vfit = 0.0;
			strike = getpStrike(idate,istrike);

			MlEqStrike::convertStrikes((*m_strikes[idate][istrike]),*strike);

			if ( m_fitSpread )
			{
				vol			= m_mktVols->getVol(*strike,currentDate, Ask);
				volAsk		= vol;

				if ( volAsk < 0.0 && m_fitSpread ){
					throw("no ask vols have been provided");
				}

				vol			= m_mktVols->getVol(*strike,currentDate, Bid);
				volBid		= vol;

				if ( volBid < 0.0 && m_fitSpread ){
					throw("no bid vols have been provided");
				}
			}

			if ( m_fitMidVols )
			{
				vol			= m_mktVols->getVol(*strike,currentDate, Middle);
				volMid		=	vol;

			}
	
			if ( m_useVega )
			{
				sig		=	vol*sqrt(time);
				d1		=	(log(1.0/strike->m_strike)+0.5*sig*sig)/sig;
				vega	=	exp(-0.5*d1*d1);
			}

			fitvol  =	getVol(*strike,idate);

			if ( m_fitMidVols )
			{
				diff	=	fitvol-volMid;
				vfit	+=	m_fitMidVols*diff*diff;
			}

			if ( m_fitSpread )
			{
				diff	=	MlEqMaths::Max(volBid-fitvol,0.0);//use these to force inside bid/ask -- usually makes no difference
				vfit	+=	m_fitSpread*diff*diff;
			}

			if ( m_fitSpread )
			{
				diff	=	MlEqMaths::Max(fitvol-volAsk,0.0);
				vfit	+=	m_fitSpread*diff*diff;
			}

			if ( m_useVega ){
				vfit *= vega*vega;
			}

			objVal += vfit;
		}
	
	*gVals	=	objVal;

}


/****************************************************************
**	Class  : 
**	Routine: findATMFwdVol
**	Returns: int
**	Action : finds location of atm vols
**  
****************************************************************/

/*static*/ int MlEqVolatilityStructure::findATMFwdVol(int& result, vector< MlEqStrikeHandle >& Strikes, double fwd)
{
// find at the money vol

	int succeed = 0;
	MlEqStrike strike;

	for (int i = 0 ; i < Strikes.size(); i++ )
	{


		MlEqStrike::convertStrikes(strike,*Strikes[i]);
		if ( fabs(strike.m_strike/fwd-1.0 ) < 1e-6 )
		{
			succeed = 1;
			result  = i;
			return succeed;
		}
	}

	return succeed;
}

/****************************************************************
**	Class  : MLHullVolData 
**	Routine: MLHullVolData
**	Returns: nothing
**	Action : assignment operator
**           
****************************************************************/

MLHullVolData& MLHullVolData::operator=(const MLHullVolData& rhs)
{
	if ( this == &rhs)
		return *this;

	m_HullData			= rhs.m_HullData; //idate]
	m_numberOfSlices	= rhs.m_numberOfSlices;

	return *this;
}

/****************************************************************
**	Class  : MLRamVolData 
**	Routine: MLRamVolData
**	Returns: nothing
**	Action : assignment operator
**           
****************************************************************/

MLRamVolData& MLRamVolData::operator=(const MLRamVolData& rhs)
{
	if ( this == &rhs)
		return *this;

	m_RamData			= rhs.m_RamData; //idate]
	m_numberOfSlices	= rhs.m_numberOfSlices;

	return *this;
}


/****************************************************************
**	Class  : MLHullVolData 
**	Routine: getVol
**	Returns: double
**	Action : vol calc
**           
****************************************************************/


double MLHullVolData::getVol(const MlEqStrike& strike,int whichSlice) 
{

	const int nsize = getNumberOfSlices();
	double vol;
	
	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	
	if ( nsize == 1)
		whichSlice = 0;	
	
	MlEqStrikeHandle pStrike;
	
	double xVal;

	pStrike = getpStrike(whichSlice,0);
	
	xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,strike);

	double xstrike = pStrike->m_strike;
	vol =  getHullVol( xstrike, m_HullData[whichSlice]);
	pStrike->m_strike = xVal;//sos do something conceptual here!

	return vol;
}


/****************************************************************
**	Class  : 
**	Routine: getHullVol
**	Returns: double
**	Action : returns vols from Hull parametrization
**  
****************************************************************/

	
double getHullVol( double normStrike, const HullCoeffs& sc)
{

	double  z, temp;

	if ( normStrike > sc.NsAtmVoltMin )
	{
		z = ( 1 - exp( -(normStrike - sc.NsAtmVoltMin) / sc.RightInflectionPoint ) );
		temp = sc.VoltMin + sc.CallWing * z * z * z;
	}
	else
	{
		z = ( 1 - exp( -(sc.NsAtmVoltMin - normStrike) / sc.LeftInflectionPoint ) );
		temp = sc.VoltMin + sc.PutWing * z * z * z;
	}

	return( temp );
}

/****************************************************************
**	Class  : MLRamVolData 
**	Routine: getVol
**	Returns: double
**	Action : vol calc
**           
****************************************************************/


double MLRamVolData::getVol(const MlEqStrike& strike,int whichSlice) 
{

	const int nsize = getNumberOfSlices();
	double vol;
	
	if( ( nsize != 1 && whichSlice >= nsize )
		|| whichSlice < 0 ){
		throw("what's going on Franz?: wrong size for whichslice encountered!");
	}
	
	if ( nsize == 1)
		whichSlice = 0;	
	
	MlEqStrikeHandle pStrike;
	double xVal;

	pStrike = getpStrike(whichSlice,0);
	xVal = pStrike->m_strike;
	MlEqStrike::convertStrikes(*pStrike,strike);
	double xstrike = pStrike->m_strike;


	vol =  getRamVol( xstrike, m_RamData[whichSlice]);
	pStrike->m_strike = xVal;//sos do something conceptual here!

	return vol;
}




void MLRamVolData::getVolSpaceDerivs(CVector& res, MlEqStrike& strike,int whichSlice) 
{
	double xstrike = strike.m_strike;
	getRamVolDerivs( res, xstrike, m_RamData[whichSlice] );
}

/****************************************************************
**	Class  : 
**	Routine: getHullVol
**	Returns: double
**	Action : returns vols from Hull parametrization
**  
****************************************************************/

	
double getRamVol( double spotBasedStrike, const RamCoeffs& sc)
{

	double  z;
	if ( spotBasedStrike < sc.DownCutoff->m_strike ){
		return sc.VolLowStrike;}

	if ( spotBasedStrike > sc.UpCutoff->m_strike ){
		return sc.VolHighStrike;}


	if ( sc.m_strike->m_strikeType == SpotBased || sc.m_strike->m_strikeType == ForwardBased){
		spotBasedStrike -= 1.0;}

	z = sc.BaseVol + sc.Slope*spotBasedStrike + sc.Curvature*pow(spotBasedStrike,2.0);
	return z;
	
}


	
void getRamVolDerivs( CVector& res, double spotBasedStrike, const RamCoeffs& sc)
{

	if ( spotBasedStrike < sc.DownCutoff->m_strike )
	{
		res[0] = res[1] = 0.;
		res[2] = sc.VolLowStrike;
		return ;
	}

	if ( spotBasedStrike > sc.UpCutoff->m_strike )
	{
		res[0] = res[1] = 0.;
		res[2] = sc.VolHighStrike;
		return ;
	}


	if ( sc.m_strike->m_strikeType == SpotBased || sc.m_strike->m_strikeType == ForwardBased){
		spotBasedStrike -= 1.0;
	}

	res[0] = 2.* sc.Curvature * spotBasedStrike + sc.Slope;
	res[1] = 2.* sc.Curvature;
	res[2] = sc.BaseVol + sc.Slope*spotBasedStrike + sc.Curvature*pow(spotBasedStrike,2.0);

//	z = sc.BaseVol + sc.Slope*spotBasedStrike + sc.Curvature*pow(spotBasedStrike,2.0);
//	return z;
	
}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getMidVolData
**	Returns: MlEqVolDataHandle
**	Action : returns midvols
**           
****************************************************************/

void MlEqVolatilityStructure::reinitializeStrikes(const MlEqAsset& asset)
{

	if ( m_BidVolData.size() ){
		for ( int islice = 0 ; islice < m_BidVolData.size();islice++ ){
			m_BidVolData[islice]->reinitializeStrikes(asset);
		}
	}

	if ( m_MidVolData.size() ){
		for ( int islice = 0 ; islice < m_MidVolData.size();islice++ ){
			m_MidVolData[islice]->reinitializeStrikes(asset);
		}
	}

	if ( m_AskVolData.size() ){
		for ( int islice = 0 ; islice < m_AskVolData.size();islice++ ){
			m_AskVolData[islice]->reinitializeStrikes(asset);
		}
	}


}



/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getMidVolData
**	Returns: MlEqVolDataHandle
**	Action : returns midvols
**           
****************************************************************/

MlEqVolDataHandle   MlEqVolatilityStructure::getMidVolData(int& numberStrikes,int islice,int istrike)const
{
	if ( m_MidVolData.size() == 0 ){
		return NULL;
	}

	if ( m_MidVolData.size() >= islice && m_MidVolData.size()!= 1)
	{
		if ( m_MidVolData[islice]->getNumberOfSlices() != 1 ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_MidVolData[islice]->getNumberOfStikes(islice);
		return m_MidVolData[islice];
	}
	else if ( m_MidVolData.size() == 1 )
	{

		if ( m_MidVolData[islice]->getNumberOfSlices() <= islice ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_MidVolData[0]->getNumberOfStikes(islice);
		return m_MidVolData[0];
	}
	else
	{
			throw("incorrect voldata set up in volstructure encountered");
	}

}

/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getMidVolData
**	Returns: MlEqVolDataHandle
**	Action : returns midvols
**           
****************************************************************/

MlEqVolDataHandle   MlEqVolatilityStructure::getBidVolData(int& numberStrikes,int islice,int istrike)const 
{

	if ( m_BidVolData.size() == 0 ){
		return NULL;
	}

	if ( m_BidVolData.size() >= islice && m_BidVolData.size()!= 1)
	{
		if ( m_BidVolData[islice]->getNumberOfSlices() != 1 ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_BidVolData[islice]->getNumberOfStikes(islice);
		return m_BidVolData[islice];
	}
	else if ( m_BidVolData.size() == 1 )
	{

		if ( m_BidVolData[islice]->getNumberOfSlices() <= islice ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_BidVolData[0]->getNumberOfStikes(islice);
		return m_BidVolData[0];
	}
	else
	{
			throw("incorrect voldata set up in volstructure encountered");
	}

}

/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getAskVolData
**	Returns: MlEqVolDataHandle
**	Action : returns midvols
**           
****************************************************************/


MlEqVolDataHandle   MlEqVolatilityStructure::getAskVolData(int& numberStrikes,int islice,int istrike)const
{
	if ( m_AskVolData.size() == 0 ){
		return NULL;
	}

	if ( m_AskVolData.size() >= islice && m_AskVolData.size() != 1 )
	{
		if ( m_AskVolData[islice]->getNumberOfSlices() != 1 ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_AskVolData[islice]->getNumberOfStikes(islice);
		return m_AskVolData[islice];
	}
	else if ( m_AskVolData.size() == 1 )
	{

		if ( m_AskVolData[islice]->getNumberOfSlices() <= islice ){
			throw("incorrect voldata set up in volstructure encountered");
		}

		numberStrikes = m_AskVolData[0]->getNumberOfStikes(islice);
		return m_AskVolData[0];
	}
	else
	{
			throw("incorrect voldata set up in volstructure encountered");
	}

}


/****************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/


double MlEqVolatilityStructure::getNaiveFutureVol(const MlEqStrike& shortStrike,const MlEqStrike& longStrike,const DATE& startdate,const DATE& enddate,BidAskMidEnum bidOrAsk) 
{	
	double vol1,vol2;

	vol1 = getVol(shortStrike,startdate,bidOrAsk); 
	vol2 = getVol(longStrike,enddate,bidOrAsk); 
	double mat  = m_hDate->GetYearFraction(enddate)-m_hDate->GetYearFraction(startdate);
	double z = vol2*vol2*m_hDate->GetYearFraction(enddate)-vol1*vol1*m_hDate->GetYearFraction(startdate);

	if ( z < 0.0 ){
			throw("Negative forwad variance encountered in forward vol calculation");
	}

	if ( mat < 0 ){
			throw("Negative maturity encountered in getFutureVol function");
	}

	double res = sqrt(z/mat);
	return res;
}


/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

/*
double MlEqVolatilityStructure::getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,double mult,BidAskMidEnum bidOrAsk) 
{

//  this is the future vol for a constant notional forward starting option


	if ( m_hDate->GetDate() == startdate ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}


	DATE d =  enddate-startdate+m_hDate->GetDate();

	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);


	double vol = getVol(fixedStrike,d,bidOrAsk);
		   vol *= mult;

	if ( vol < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}

	if ( !m_asymptodicVolIsInitialized  ){

		double res = getNaiveFutureVol(fixedStrike,fixedStrike,startdate,enddate,bidOrAsk); 
		return res;
	}

//8888
	MlEqStrike fixedStrikeATM(spot);
	double volATMspot		= m_asymtodicVolStruct->getVol(fixedStrikeATM,d,bidOrAsk);
	double volATMspotNaive	= getNaiveFutureVol(fixedStrikeATM,fixedStrikeATM,startdate,enddate,bidOrAsk);
	double asmult;

	if ( fabs(volATMspotNaive) < 1e-3 ){
		asmult = 0.0;}
	else{
		asmult = volATMspotNaive/volATMspot;
	}


	double volAs = m_asymtodicVolStruct->getVol(fixedStrike,d,bidOrAsk);
	volAs *= asmult;

	if ( volAs < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}

	double t = m_hDate->GetYearFraction(startdate);
	vol = vol + (1.0-exp(-m_decayFactor*t))*(volAs-vol);
	return vol;

}
*/

/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,BidAskMidEnum bidOrAsk) 
{
//  returns  jim gatheral vol
	double res = getFutureVol(strike,startdate,enddate,spot,spot,spot,bidOrAsk) ;
	return res;

//  this is the future vol for a constant notional forward starting option


/*	if ( m_hDate->GetDate() == startdate ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	DATE d =  enddate-startdate+m_hDate->GetDate();

	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);


	MlEqStrike fixedStrikeATM(spot);
	double volATMspot = getVol(fixedStrikeATM,d,bidOrAsk);
	double volATMspotNaive = getNaiveFutureVol(fixedStrikeATM,fixedStrikeATM,startdate,enddate,bidOrAsk);
	double mult;

	if ( fabs(volATMspotNaive) < 1e-3 ){
		mult = 0.0;}
	else{
		mult = volATMspotNaive/volATMspot;
	}

	double res = MlEqVolatilityStructure::getFutureVol(strike,startdate,enddate,spot,mult,bidOrAsk) ;
	return res;
*/
}

/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,double fwdShort,double fwdLong,BidAskMidEnum bidOrAsk) 
{

//  this is the future vol for a constant notional forward starting option


	if ( m_hDate->GetDate() == startdate ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	DATE d =  enddate-startdate+m_hDate->GetDate();

	MlEqStrike fixedStrike;
	MlEqStrike::convertStrikes(fixedStrike,strike);

	MlEqStrike fixedStrikeATM(spot);
	double volATMspot = getVol(fixedStrikeATM,d,bidOrAsk);
	double volATMspotNaive = getNaiveFutureVol(MlEqStrike(fwdShort),MlEqStrike(fwdLong),startdate,enddate,bidOrAsk);
	double mult = volATMspotNaive/volATMspot;

	double vol = getVol(fixedStrike,d,bidOrAsk);
		   vol *= mult;

	if ( vol < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}

	if ( !m_asymptodicVolIsInitialized  ){

		double res = getNaiveFutureVol(fixedStrike,fixedStrike,startdate,enddate,bidOrAsk); 
		return res;
	}

	volATMspot		= m_asymtodicVolStruct->getVol(fixedStrikeATM,d,bidOrAsk);
	double asmult;

	if ( fabs(volATMspotNaive) < 1e-3 ){
		asmult = 0.0;}
	else{
		asmult = volATMspotNaive/volATMspot;
	}


	double volAs = m_asymtodicVolStruct->getVol(fixedStrike,d,bidOrAsk);
	volAs *= asmult;

	if ( volAs < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}

	double t = m_hDate->GetYearFraction(startdate);
	vol = vol + (1.0-exp(-m_decayFactor*t))*(volAs-vol);
	return vol;

}



/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startDate,long enddate,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) 
{

//	volref: target value of the renornmalization point

	long startdate  = startDate.GetDate();
	double fwdShort = asset.GetForward(startdate,false);
	double fwdLong  = asset.GetForward(enddate,false);

	double volRefNaive	= getNaiveFutureVol(MlEqStrike(fwdShort),MlEqStrike(fwdLong),startdate,enddate,bidOrAsk);

	MlEqDateHandle startDateHandle = new MlEqDate(startDate);
	MlEqStrikeHandle strikeInv;
	strikeInv = setUpDefaultStrikeInvariant(strikeInvariant,asset,startDateHandle,enddate);

	double vol = MlEqVolatilityStructure::getFutureVol(asset,strike,startDate,enddate,volRefNaive,strikeInv,bidOrAsk) ;
	return vol;
}


/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startDate,long enddate,double volRef,BidAskMidEnum bidOrAsk) 
{

//	volRef: target value of the renornmalization point

	MlEqStrikeHandle strikeInvariant;


	MlEqDateHandle startDateHandle = new MlEqDate(startDate);

//	double fwd = asset.GetForward(startDate.GetDate(),enddate,false);
//	strikeInvariant = new MlEqNormalizedStrike(0.0, enddate,volRef,fwd);

	strikeInvariant = setUpDefaultStrikeInvariant(strikeInvariant,asset, startDateHandle,enddate);

	double vol = getFutureVol(asset,strike,startDate,enddate,volRef,strikeInvariant, bidOrAsk) ;
	return vol;
}






/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startDate,long enddate,BidAskMidEnum bidOrAsk) 
{

//	volref: target value of the renornmalization point

	long startdate = startDate.GetDate();
	double fwdShort = asset.GetForward(startdate,false);
	double fwdLong  = asset.GetForward(enddate,false);

	double volRefNaive	= getNaiveFutureVol(MlEqStrike(fwdShort),MlEqStrike(fwdLong),startdate,enddate,bidOrAsk);

	double vol = getFutureVol(asset,strike,startDate,enddate,volRefNaive,bidOrAsk) ;

	return vol;
}


/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startDate,long enddate,double volRef,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) 
{

//  this is the future vol for a constant notional forward starting option

	if ( !strikeInvariant){
		throw("no strikeInvariant handle has been entered");
	}

	long startdate	= startDate.GetDate();
	
	if ( m_hDate->GetDate() == startdate ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;// that's the end of the story in
	}

	double fwdShort = asset.GetForward(startdate,false);
	double fwdLong  = asset.GetForward(enddate,false);

	if ( !m_asymptodicVolIsInitialized  )
	{
		double vol	 =  getNaiveFutureVol(fwdShort,fwdLong,startdate,enddate,bidOrAsk);
		double shift =  shiftFwdSkew(startDate.GetDate(), enddate, volRef, strikeInvariant);
		vol += shift;
		return vol; // that's the end of the story
	}

//  create strike variable that defines the parrallel transport

	MlEqStrikeHandle hStrike;
	hStrike = strikeInvariant->copy();
	MlEqStrike::convertStrikes(*hStrike,strike); // converts the input strike into the scale invariant strike type

	DATE d =  enddate-startdate+m_hDate->GetDate();
	hStrike->SetStartDateHandle(m_hDate);
	hStrike->SetEndDate(d);
	hStrike->reinitializeStrike(asset);// reset internals of hStrike to spot market data

	MlEqStrikeHandle hStrikeInv;
	hStrikeInv = strikeInvariant->copy();
	hStrikeInv->SetStartDateHandle(m_hDate);
	hStrikeInv->SetEndDate(d);
	hStrikeInv->reinitializeStrike(asset);// reset internals of hStrike to spot market data

//  start with spot vol contribution

	double volref	= getVol(*hStrikeInv,d,bidOrAsk);
	double mult		= volRef/volref;// this makes sure that fwd vols are centered around the naive forward

	double vol		= getVol(*hStrike,d,bidOrAsk);
		   vol		*= mult;

	if ( vol < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}


//  continue with asymptodic component first

//  the actual strike value in strikeInvariant defines the renormalization point of the asymtpdic skew
//  volAsRef the value at the renormalization point

	double volAsRef		= m_asymtodicVolStruct->getVol(*hStrikeInv,d,bidOrAsk);

	double asmult;
	if ( fabs(volAsRef) < 1e-3 ){
		asmult = 0.0;}
	else{
		asmult = volRef/volAsRef;
	}

	double volAs = m_asymtodicVolStruct->getVol(*hStrike,d,bidOrAsk);
	volAs *= asmult;

	if ( volAs < -0.1 ){
		throw("Negative vol encountered in getVol function");
		return -1;
	}


//  compose values to get asymprodic vols

	double t = m_hDate->GetYearFraction(startdate);
	vol = vol + (1.0-exp(-m_decayFactor*t))*(volAs-vol);

	double shift =  shiftFwdSkew( startDate.GetDate(),enddate,volRef,strikeInvariant); 
	vol += shift;

	return vol;


}


/*************************************************************
**	Class  : MlEqVolatilityStructure 
**	Routine: getVol
**	Returns: vol
**	Action : returns vol
**           
****************************************************************/

double MlEqVolatilityStructure::shiftFwdSkew(long startdate, long enddate, double volRef, const MlEqStrikeHandle& strikeInvariant)
{
	double f = m_scenarios.Apply(0.0, strikeInvariant, startdate, enddate, volRef, stForwardSkewShift);
	return f;	
}


/****************************************************************
**	Class  : none
**	Routine: setUpDefaultStrikeInvariant
**	Returns: MlEqStrikeHandle
**	Action : sets up default strike invariant handle for parralel transport
**           
****************************************************************/

/*static*/ MlEqStrikeHandle MlEqVolatilityStructure::setUpDefaultStrikeInvariant(MlEqStrikeHandle strikeInvariant,const MlEqAsset& asset,MlEqDateHandle startDateHandle,long enddate)
{	
	if ( !!strikeInvariant ){
		return strikeInvariant;
	}

	MlEqStrikeHandle strikeInv;

	double fwdLong = asset.GetForward(enddate,false);
	double fwdShort  = asset.GetForward(startDateHandle->GetDate(),false);

	double vol = asset.GetNaiveFutureVolatility(MlEqStrike(fwdShort),MlEqStrike(fwdLong),startDateHandle->GetDate(),enddate); 

	strikeInv = new MlEqNormalizedStrike(0.0, enddate,vol,fwdLong,
						 true,startDateHandle);

	return strikeInv;
}



/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: initialize
**	Returns: nothing
**	Action : returns vol
**           
****************************************************************/

MlEqVolatilityStructure::MlEqVolatilityStructure( 
						const MlEqVolatilityStructureHandle& spotVolStruct,
						const MlEqVolatilityStructureHandle& asymptodicVolStruct,
						const double					decayFactor)
{
	
	m_decayFactor		 = decayFactor;
//	m_asymtodicVolStruct = asymptodicVolStruct;sos
	MlEqVolatilityStructureHandle* pThis = (MlEqVolatilityStructureHandle*)this;
	pThis->MlEqVolatilityStructureHandle::operator=(spotVolStruct);
}	

	
	
void MlEqVolatilityStructure::initialize( 
						const MlEqVolatilityStructure& spotVolStruct,
						const MlEqVolatilityStructure& asymptodicVolStruct,
						double	decayFactor)
{	
	

	m_decayFactor		 = decayFactor;
	m_asymtodicVolStruct = new MlEqVolatilityStructure(asymptodicVolStruct);
	m_asymptodicVolIsInitialized	=	1;
	*this = spotVolStruct;

}

//	adjusts the volatility structure reference date
void MlEqVolatilityStructure::putDateToDouble(MlEqDateHandle hDate)
{		
	m_hDate = hDate;
	// ToDo - Alex - put any reinitialisation code here
};
void MlEqVolatilityStructure::PutDateNumber(long nDate)
{
	m_hDate->PutDate(nDate);
	// ToDo - Alex - put any reinitialisation code here
}




/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/

double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,BidAskMidEnum bidOrAsk)
{

	if ( m_hDate->GetDate() == startdate ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	double vol = MlEqVolatilityStructure::getFutureVol(strike,startdate,enddate,spot,bidOrAsk) ;
	vol *= getLevel(startdate,enddate);

	return vol;
}

/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/


double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdateHandle,long enddate,double volRef,BidAskMidEnum bidOrAsk) 
{


	if ( m_hDate->GetDate() == startdateHandle.GetDate() ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	double vol = MlEqVolatilityStructure::getFutureVol(asset,strike,startdateHandle,enddate,volRef,bidOrAsk) ;
	vol *= getLevel(startdateHandle.GetDate(),enddate);

	return vol;

}

/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/

double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdate,long enddate,BidAskMidEnum bidOrAsk) 
{

	if ( m_hDate->GetDate() == startdate.GetDate() ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	double vol = MlEqVolatilityStructure::getFutureVol( asset,strike,startdate,enddate,bidOrAsk) ;
	vol *= getLevel(startdate.GetDate(),enddate);
	return vol;
}


/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/

double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdate,long enddate,double volRef,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) 
{

	if ( m_hDate->GetDate() == startdate.GetDate() ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	double vol = MlEqVolatilityStructure::getFutureVol(asset,strike,startdate,enddate,volRef,strikeInvariant,bidOrAsk) ;
	vol *= getLevel(startdate.GetDate(),enddate);

	return vol;

}


/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/

double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdate,long enddate,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) 
{

	if ( m_hDate->GetDate() == startdate.GetDate() ){
		double vol = getVol(strike,(double)enddate,bidOrAsk) ;	
		return vol;
	}

	double vol = MlEqVolatilityStructure::getFutureVol(asset,strike,startdate,enddate,strikeInvariant,bidOrAsk) ;
	vol *= getLevel(startdate.GetDate(),enddate);

	return vol;
}


/****************************************************************
**	Class  : MlEqAsymptodicVolatilityStructure 
**	Routine: getFutureVol
**	Returns: double
**	Action : 
**           
****************************************************************/


double MlEqStochBetaVolatilityStructure::getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,double fwdShort,double fwdLong,BidAskMidEnum bidOrAsk) 
{
	double vol = MlEqVolatilityStructure::getFutureVol(strike,startdate,enddate,spot,fwdShort,fwdLong,bidOrAsk) ;
	vol *= getLevel(startdate,enddate);
	return vol;

}

/****************************************************************
**	Class  : MlEqStochBetaVolatilityStructure
**	Routine: getFrquencyRoundUp
**	Returns: double
**	Action : 
**           
****************************************************************/

// this should work if you don't use stochastic beta model to price 
// products with a frequency over 5 years...

/*static*/ double MlEqStochBetaVolatilityStructure::getFrequencyRoundUp(DATE dateStart, DATE dateEnd)
{
	long nbMonth = long((dateEnd - dateStart)/30.416667) ;	// very shitty
	MlEqDate start(dateStart), end(dateEnd);
	MlEqDate temp(	*(start.Add(0, nbMonth, 0, NoBusinessDayConvention, "", false)) );

	if((temp.Add(0,1,0, NoBusinessDayConvention, "", false))->GetDate() <= dateEnd)
	{
		temp = *(temp.Add(0,1,0, NoBusinessDayConvention, "", false));
		nbMonth++;
	}

	double nDays = (temp.Add(0,1, 0, NoBusinessDayConvention, "", false))->GetDate() - temp.GetDate() ;
	nDays = (dateEnd - temp.GetDate()) / nDays ;

	return nbMonth + nDays;
}

