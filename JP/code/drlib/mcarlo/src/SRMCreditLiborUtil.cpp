//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditLiborUtil.cpp
//
//   Description : Derived SRMCreditUtil class for CreditLibor model
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMCreditLiborUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/CRCalib.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"
#include "edginc/SwapTool.hpp"

#define NUMRESETDATES = 10

DRLIB_BEGIN_NAMESPACE

SRMCreditLiborUtil::SRMCreditLiborUtil(
    const DateTime&         baseDate,
    const string&           smileParamsKey,
    ICDSParSpreadsConstSP   stochCDSCurve): SRMCreditUtil(baseDate, smileParamsKey, stochCDSCurve)
{
    static const string method("SRMCreditLiborUtil::SRMCreditLiborUtil");
    try{
        // back out smile parameters
        CRCalib::SmileRequest smileRequest(smileParamsKey);
        CVolProcessed* vol = stochCDSCurve->getProcessedVol(&smileRequest);
        smartPtr<CRCalib::VolProcessed> volData(
            &dynamic_cast<CRCalib::VolProcessed&>(*vol));
        qLeft = volData->getQLeft();
        qRight = volData->getQRight();
		

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void SRMCreditLiborUtil::setTimeLine(DateTimeArrayConstSP simDates) // called when we know allDates
{
	int i;
    mLastSimTime       = 70.;   // in years
    mIntensityInterval = 12; 	// in months

	double temp;

    static const string method("SRMCreditLiborUtil::setTimeLine");
    
	if (initialized) 
        return; // nothing to do

    try{
        initialized = true;
    
        dates = simDates;

        // back out model parameters 
		// CHECK WE ARE GETTING THE RIGHT VOLS!!!!!!!!!!!
        MRSpotVolRequest betaRequest;
        CVolProcessed* volBeta = stochCDSCurve->getProcessedVol(&betaRequest);
        MRSpotVolProcessedSP volBetaData(
            &dynamic_cast<MRSpotVolProcessed&>(*volBeta));

        calcExtendedTimeLine(); // get 'extended' timeline

		// initialize dates: Accrued, Reset, etc...
	
		//****************************************************************
		//					WARNING: UNDER CONSTRUCTION
		//*****************************************************************
        DateTimeArray irregularSchedule(4);

	
		irregularSchedule[0] = baseDate.rollDateInMonths(3); 
		irregularSchedule[1] = baseDate.rollDateInMonths(6);
		irregularSchedule[2] = baseDate.rollDateInMonths(9);
		irregularSchedule[3] = baseDate.rollDateInMonths(12); 
		
		// regular schedule starts after one year
		DateTime startRegular = baseDate.rollDateInMonths(12); 

		DateTimeArray* regularSchedule(
		SwapTool::dateArray(
                            startRegular,							// start here
							mIntensityInterval,						// interval = count periods
							"M",									// e.g. Y, M, W, D
							1,										// 0=start @ basedate, 1=start @ baseDate + interval
							1,										// arrayIncrement, usually +1 or -1
							mLastSimTime * 12 / mIntensityInterval));	// how many dates, going out to 70 years

		vector<const DateTimeArray*> dtArrayVector(2);
		dtArrayVector[0] = &irregularSchedule;
		dtArrayVector[1] = regularSchedule;							   
		DateTimeArray TPDate(DateTime::merge(dtArrayVector));
		
		// sort dates and remove duplicates
		DateTime::doSortUniq(TPDate);

		mNumInt = TPDate.size() - 1;									// the last date is not a reset date
		
		// allocate memory for arrays of base class
		mResetDates.resize(mNumInt);
		mSpotSurvProb.resize(mNumInt+1);								// we need to consider the surv prob after the last reset date
		mAccruedFrac.resize(mNumInt);
		mInitialIntensities.resize(mNumInt);

	 
		temp = 0.0;
		mYearFracToFirstReset = baseDate.yearFrac(TPDate[0]);

		// initialize spot survival probabilities for each reset date	
		for (i = 0; i< mNumInt+1; i++)
		{	
			temp = logSurvProbRatio(baseDate,TPDate[i]);
			mSpotSurvProb[i] = exp(temp); 
		}
		

		for (i = 0; i< mNumInt; i++)
		{
			mResetDates[i]		   = baseDate.yearFrac(TPDate[i]);
			mAccruedFrac[i]		   = TPDate[i].yearFrac(TPDate[i+1]) ; 
			mInitialIntensities[i] = (1./mAccruedFrac[i]) *(mSpotSurvProb[i]/mSpotSurvProb[i+1]-1.);
		}
	
		DateTimeArray volDates;
		volDates = TPDate;
		volDates.pop_back();				// last TPDate is not a reset date


        // populate SpotVol (so far we only have FlatCDSSpotVol ...)
		// SpotVols are the discrete intensity vols
		// WARNING: mSpotVol[i] refers to the period T_i,T_{i+1}, where T_i is a reset date
		// CHECK THIS when you move to a non-flat vol term structure!!!!!!!!!!!!
        DoubleArray spotVolTemp(mNumInt);
        volBetaData->spotVol(baseDate, volDates, spotVolTemp);
        mSpotVols = vector<double>(spotVolTemp.begin(), spotVolTemp.end());


    } catch (exception& e){
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
