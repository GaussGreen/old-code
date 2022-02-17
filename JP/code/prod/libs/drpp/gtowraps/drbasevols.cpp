// drbasevols.cpp: implementation of the drbasevols class.
//
//////////////////////////////////////////////////////////////////////

#include "drbasevols.h"
#include "drsymboltable.h"
extern "C" {
#include "lintrp.h"
#include "vtfmwrap.h"
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const DRString KAPITAL_BASE_VOL	= "basevol.dat";

DRBaseVols::DRBaseVols()
: DRCurve(), m_volMat(0), m_freq(0), m_fromCapVol(false), m_originalMat(-999), m_originalFreq(-999)
{}

DRBaseVols::DRBaseVols (DRDate& envDate, MbsIO& mbsIO)
: DRCurve(), m_freq(0), m_fromCapVol(false), m_originalMat(-999), m_originalFreq(-999)
{
	MbsIOMap& mbsIOMap = mbsIO.get_map();

	MbsIOVector& dateVec = mbsIOMap["BASEVOLDATES"].get_vector();
	MbsIOVector& rateVec = mbsIOMap["BASEVOLRATES"].get_vector();
	
	if (dateVec.size() != rateVec.size()) 
		throw DRException ("Size mismatch");

	DRDateList dates;
	DArray rates (dateVec.size());

	for (int i = 0; i < dateVec.size(); i++) {
		dates.AddDate (toDate(dateVec[i].get_string()));
		rates[i] = toNumber(rateVec[i].get_string()) / 100;
	}

	m_freq = toNumber(mbsIOMap.get("BaseVolFreq"));
	m_volMat = toNumber(mbsIOMap.get("BaseVolMat"));

	double basis = 12. / m_volMat;

	long dayCountConv = GTO_ACT_365F;

	m_tcurve = GtoMakeTCurve(envDate, dates, rates, dates.size(), basis, 
		dayCountConv);
}

DRBaseVols::~DRBaseVols()
{}

DRBaseVols::DRBaseVols (TCurve* volCurve, int volFreq, int volMat)
:DRCurve (volCurve), m_volMat(volMat), m_freq(volFreq), m_originalMat(-999), m_originalFreq(-999)
{
	m_fromCapVol = true;
	m_tcurve->fBasis = 12. / m_volMat;
}

DRBaseVols::DRBaseVols (FileType fileType, DRDate& envDate, 
						DRString filename, int mat)
						: DRCurve(), m_volMat(0), m_freq(0), m_originalMat(-999), m_originalFreq(-999)
{
	m_fromCapVol = true;
	if (fileType == SWAP) {
		LoadFromSwapFiles (envDate, filename, mat);
	}
	else {
		LoadFromKapitalFiles (envDate);
	}
}

DRBaseVols::DRBaseVols (DRDate& envDate, DRSwapVols& swapVols, double matForCalib)
					:DRCurve(), m_volMat(0), m_originalMat(-999), m_originalFreq(-999)
{
	// June 9, 1999 - MKS, Use todays date instead of envDate for vols
	// because KAPITAL is passing Vol dates based on today date

	DRString todayDateString("KAPITALTODAYDATE");
	DRDate todayDate = toDate(todayDateString);
	m_fromCapVol = false;
	m_freq = swapVols.freq();
	m_volMat= (int) matForCalib * 12;

	double basis = 12. / m_volMat;

	int numExp = swapVols.exps().size();
	int numMat = swapVols.mats().size();

	m_tcurve = GtoNewTCurve (envDate, numExp, basis, GTO_ACT_365F);

	DArray volMat (numExp);
	volMat = matForCalib;

	DArray matrix1D (numExp * numMat);
	MbsTDateArray volDates(numExp);
	double yearFrac;
	int iExp;

	for (iExp = 0; iExp < numExp; iExp++) {

		for (int iMat = 0; iMat < numMat; iMat++) {
			matrix1D[iExp + iMat * numExp] = swapVols.matrix()[iExp][iMat];
		}

		
		yearFrac = swapVols.exps()[iExp];
		if (yearFrac < 1./360)
			m_tcurve->fArray[iExp].fDate = (TDate) (todayDate + ONE_DAY);
		else 
			m_tcurve->fArray[iExp].fDate = (TDate) (todayDate + ((int) ceil(12*yearFrac)) * ONE_MONTH);
	}


	DArray vols(numExp);

	if (GtoLinInterpDoubleArray2 (
		swapVols.exps(),
		sizeof(double),
		numExp,
		swapVols.mats(),
		sizeof(double),
		numMat,
		matrix1D,
		sizeof(double),
		swapVols.exps(),
		sizeof(double),
		volMat,
		sizeof(double),
		numExp,
		NULL,
		sizeof(double),
		vols) IS FAILURE)
		throw DRException("Failed in GtoLinInterpDoubleArray2");

	for (iExp = 0; iExp < numExp; iExp++) {
		m_tcurve->fArray[iExp].fRate = vols[iExp];
	}
}

bool DRBaseVols::operator==(const DRBaseVols& a) const
{
	bool ans;
	if (m_volMat == a.m_volMat)
		ans = ((DRCurve&) *this == (DRCurve) a) && m_freq == a.m_freq;
	else if (m_originalMat == a.m_originalMat) {
		ans = (m_originalCurve == a.m_originalCurve) && m_originalFreq == a.m_originalFreq;
	}
	else if (m_volMat == a.m_originalMat) {
		ans = ((DRCurve&) *this == a.m_originalCurve) && m_freq == a.m_originalFreq;
	}
	else if (a.m_volMat == m_originalMat) {
		ans = (m_originalCurve == (DRCurve) a) && m_originalFreq == a.m_freq;
	}
	else {
		ans = false;
	}
	return ans;
}

void DRBaseVols::ChangeToSimpleBaseVolInLognormalSpace (long newMat,
														DRDate& envDate, 
														DRCurve& zeroCurve, 
														DArray& calibParams)
{
	if (newMat == maturity()) return;

	m_originalCurve = (DRCurve) *this;
	m_originalMat = m_volMat;
	m_originalFreq = m_freq;

	MbsTDateArray baseDatesL (3) ;
	baseDatesL[0] = 2;
	baseDatesL[1] = envDate;
	baseDatesL[2] = zeroCurve.baseDate();

	DRDateList zeroDates = zeroCurve.GetDates();
	DArray zeroRates = zeroCurve.GetRates();

	DRDateList volDates = GetDates();
	DArray vols = GetRates();

	MbsTDateArray tempZeroDates (zeroDates, zeroDates.size());
	MbsTDateArray zcDateL = MakeCountedArray(tempZeroDates);
	DArray zcRateL = MakeCountedArray(zeroRates);

	MbsTDateArray tempVolDates (volDates, volDates.size());
	MbsTDateArray volDatesL = MakeCountedArray(tempVolDates);
	DArray volsL = MakeCountedArray(vols);

	DArray tempCalibL(3);
	tempCalibL[0] = 2;
	tempCalibL[1] = calibParams[1];
	tempCalibL[2] = calibParams[2];

    LArray finalMatL(2);
	finalMatL[0] = 1;
	finalMatL[1] = 0;
 
    LArray volMatDateL(2);
	volMatDateL[0] = 1;
	volMatDateL[1] = 0;

    LArray volFreqL(2);
	volFreqL[0] = 1;
	volFreqL[1] = freq();

	DArray volMatL(2);
	volMatL[0] = 1;
	volMatL[1] = m_volMat / 12.;

	DArray bvMatL(2);
	bvMatL[0] = 1;
	bvMatL[1] = newMat / 12.;

	DArray bvRatesL (volDatesL[0] +1);
	bvRatesL[0] = volDatesL[0];

	if (VtfmVolConvert1VCurveL
            (baseDatesL,       // today and value date 
             zcDateL,          // Zero coupon dates
             zcRateL,          // Zero coupon rates (Ann,Act/365F) 
             tempCalibL,       // Calib parms (NDim,mrs,cors,wghts)
             volDatesL,        // Calibration dates 
             finalMatL,        // Input mat: 1=Use Date;0=Use Intval
             volMatDateL,      // Input vol maturity date 
             volMatL,          // Input vol maturity interval 
             volFreqL,         // Input vol frequency 
             volDatesL,        // Input vol dates 
             volsL,            // Input vols; interped at calDates 
             bvMatL,           // Output vol maturity 
             volDatesL,        // Output vol dates 
             bvRatesL)       // Output vols (freq=0) 
		IS FAILURE)
		throw DRException("Failed in VtfmVolConvert1VCurveL (covert vols)");

/*
	cout << GtoFormatDate(baseDatesL[1]) ;
	cout << zcDateL << zcRateL << tempCalibL << volDatesL ;
	cout << finalMatL << volMatDateL  << volMatL << volFreqL ;
	cout << volDatesL << volsL  << bvMatL  << volDatesL ;
	cout << bvRatesL ;
*/

	for (int i = 0; i < m_tcurve->fNumItems; i++) {
		m_tcurve->fArray[i].fRate = bvRatesL[i+1];
	}
	m_tcurve->fBasis = 12 / newMat;
	m_volMat = newMat;
	m_freq = 0;
}

void DRBaseVols::Print(ostream& stream) const
{
	DRCurve::Print(stream);
	stream << "Maturity:\t" << m_volMat << endl;
	stream << "Freq:\t" << m_freq << endl;
}

void DRBaseVols::LoadFromKapitalFiles (DRDate& envDate)
{
	ifstream inFile(KAPITAL_BASE_VOL.c_str()); 
	if (!inFile) 
		throw DRException ("Failure to Open Kapital File ") << KAPITAL_BASE_VOL;

	char trash[255];
	char freqLetter;

	inFile.getline(trash, 255);
	inFile >> freqLetter;

	if (freqLetter == 'Q') m_volMat = 3;
	else if (freqLetter == 'S') m_volMat = 6;
	int basis = 12 / m_volMat;

	inFile.getline(trash, 255);
	inFile.getline(trash, 255);

	int numRates;
	inFile >> numRates;

	inFile.getline(trash, 255);
	inFile.getline(trash, 255);

	MbsTDateArray dates (numRates);
	DArray rates (numRates);

	long tempDate;
	for (int i = 0 ; i < numRates ; i++) {
		inFile >> tempDate;
		dates[i] = toDate(tempDate);
		inFile >> rates[i];
	}

	rates /= 100.;

	m_tcurve = GtoMakeTCurve (envDate, dates, rates, numRates, basis, GTO_ACT_365F);	
}

void DRBaseVols::LoadFromSwapFiles (DRDate& envDate, DRString filename, int mat)
{
	int col;

	if (mat == 1) col = 0;
	else if (mat == 3) col = 1;
	else if (mat == 6) col = 2;
	else if (mat == 12) col = 3;
	else 
		throw DRException ("Illegal maturity");

	m_volMat = mat;
	int basis = 12/ mat;

	char temp[20];

	ifstream inFile(CheckPathDelimiters(theSymbolTable.get(filename)).c_str());
	if (!inFile) 
		throw DRException("Failed to open file: ") << filename;

	DArray rates(20);
	MbsTDateArray dates(20);

	int i;
	for (i = 0; i<6; i++) inFile >> temp;
	
	for (i = 0; i<20; i++) {
		DRString tempString;
		
		inFile >> temp;
		inFile >> temp;
		tempString = temp;
		
		tempString = left_substring(tempString, -1);
		tempString = right_substring(tempString, -1);
		
		TDate tempDate;
		char* t = const_cast<char*> (tempString.c_str());
		if (GtoStringToDate(t, &tempDate) IS FAILURE)
			throw DRException("Failure to convert string to date");
		
		dates[i] = tempDate;
		
		int j;
		for (j = 0; j < col; j++) inFile >> temp;
		
		inFile >> temp;
		tempString = temp;
		rates[i] = atof (left_substring(tempString, -1).c_str());
		
		int numRowsRemain = 3 - col;
		for (j = 0; j < numRowsRemain; j++) inFile >> temp;
	}

	rates /= 100;

	m_tcurve = GtoMakeTCurve (envDate, dates, rates, 20, basis, GTO_ACT_365F);
}

