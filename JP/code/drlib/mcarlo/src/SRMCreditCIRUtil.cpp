//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditCIRUtil.cpp
//
//   Description : Derived SRMCreditUtil class for CIR model
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditCIRUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/CRCalib.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
SRMCreditCIRUtil::SRMCreditCIRUtil(
    const DateTime&         baseDate,
    const string&           smileParamsKey,
    ICDSParSpreadsConstSP   stochCDSCurve):
    SRMCreditUtil(baseDate, 
                  smileParamsKey, 
                  stochCDSCurve)
{
    static const string method("SRMCreditCIRUtil::SRMCreditCIRUtil");
    try {
        // back out smile parameters
        CRCalib::SmileRequest smileRequest(smileParamsKey);
        CVolProcessed* vol = stochCDSCurve->getProcessedVol(&smileRequest);
        smartPtr<CRCalib::VolProcessed> volData(
            &dynamic_cast<CRCalib::VolProcessed&>(*vol));
        qLeft = volData->getQLeft();
        qRight = volData->getQRight();

        // back out model parameters (in CR, there is only one model parameter, namely beta)
        MRSpotVolRequest betaRequest;
        CVolProcessed* volBeta = stochCDSCurve->getProcessedVol(&betaRequest);
        MRSpotVolProcessedSP volBetaData(&dynamic_cast<MRSpotVolProcessed&>(*volBeta));

        beta = volBetaData->meanReversion(); 

        DoubleArray spotVolTemp(1);
        volBetaData->spotVol(baseDate, DateTimeArray(1, baseDate), spotVolTemp);
        SpotVol = spotVolTemp[0];

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** destructor */
SRMCreditCIRUtil::~SRMCreditCIRUtil() {
}

/** returns spot vols */
double SRMCreditCIRUtil::getSpotVol() const { 
    return SpotVol; 
}

/** returns beta */
double SRMCreditCIRUtil::getBeta() const { 
    return beta; 
}

/** returns qRight */
double SRMCreditCIRUtil::getQRight() const { 
    return qRight; 
}

/** returns qLeft */
double SRMCreditCIRUtil::getQLeft() const { 
    return qLeft; 
}

/** finishes initialization */
// called when we know allDates
void SRMCreditCIRUtil::setTimeLine(DateTimeArrayConstSP simDates) {
    static const string method("SRMCreditCIRUtil::setTimeLine");
    if (initialized) 
        return; // nothing to do

    try {
        initialized = true;
    
        dates = simDates;

        calcExtendedTimeLine(); // get 'extended' timeline
    
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
