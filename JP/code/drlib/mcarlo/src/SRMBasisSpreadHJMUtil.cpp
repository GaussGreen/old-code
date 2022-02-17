//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : SRMBasisSpreadHJMUtil.cpp
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMBasisSpreadHJMUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SPCalib.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"

USING_DRLIB_NAMESPACE

///// constructs and populates SRMBasisSpreadUtil object
SRMBasisSpreadUtil::SRMBasisSpreadUtil(
    const DateTime&       _baseDate,
    IBasisIndexCurveConstSP _basisIndexCurve)
    :
        baseDate(_baseDate),
        basisCurve(_basisIndexCurve),
        beta(0.0),
        q(0.0),
        initialized(false)
{
}

SRMBasisSpreadUtil::~SRMBasisSpreadUtil() {}


void SRMBasisSpreadUtil::setTimeLine(DateTimeArrayConstSP allDates) // finishes initialization
{
    static const string method = "SRMBasisSpreadUtil::setTimeLine()";

    if (initialized) {
        if (allDates->size() != dates->size() || ! DateTime::isSubset(*allDates, *dates))
            throw ModelException(method, "Re-initialized with a different timeline");
        return;
    }
    initialized = true;
    try {

        dates = allDates;
        fwdSpreads.resize(dates->size()-1);
        for(size_t i=0; i<fwdSpreads.size(); ++i)
            fwdSpreads[i] = calcFwdSpread((*dates)[i]);

        // get q, vol, meanrev
        SPCalib::SPCalibRequest spVolRequest;
        CVolProcessed* vol = basisCurve->getProcessedVol(&spVolRequest);
        smartPtr<SPCalib::Processed> volData (
            &dynamic_cast<SPCalib::Processed&>(*vol));

        beta = volData->meanReversion();
        q = volData->getQSpread();

        // populate SpotVol (so far we only have FlatCDSSpotVol ...)
        DoubleArray spotVolTemp(dates->size());
        volData->spotVol(baseDate, *dates, spotVolTemp);
        spotVols.assign(spotVolTemp.begin(), spotVolTemp.end());

        //TODO: see if anything else is needed, e.g. for efficiency

    } catch (exception& e){
        throw ModelException(e, method, "For basis currency "+basisCurve->getName());
    }

}

/** Utility routine to calculate the cutoff rate and populate
    the cutoff rate arrays */
void SRMBasisSpreadUtil::calcEffRateLimit(
                        double            NbSigmasMax,    // (I) Number or sigmas to cut at
                        double            NbSigmasMin,    // (I) Number or sigmas to cut at
                        vector<double>&   MaxRate,        // (O) Cutoff forward spreads
                        vector<double>&   MinRate) const // (O) Cutoff forward spreads
{
    if (dates->empty())
        return ; // might throw ModelException instead...

    MaxRate.resize(dates->size()-1);
    MinRate.resize(dates->size()-1);
    double var = 0;

    // TODO: make sure all this (i+1) is correct -- maybe not, or QSRM might use different convention
    MaxRate[0] = fwdSpreads[0];
    MinRate[0] = Maths::min(fwdSpreads[0],0.0);
    for (size_t i = 0; i<MaxRate.size()-1; i++)
    {
        double logVol = spotVols[i] / fwdSpreads[i];

        var += SRMYearFrac((*dates)[i], (*dates)[i+1]) * logVol * logVol;

        if (var < 0.0) return;

        /* lognormal MAX-Cutoff */
        double dummy = NbSigmasMax * sqrt(var);
        dummy = Maths::min(dummy, 100.0); /* to avoid blow-up in exp() */
        MaxRate[i+1] = fwdSpreads[i+1] * exp(dummy);

        /* normal MIN-Cutoff */
        MinRate[i+1] = fwdSpreads[i+1] * (1 - NbSigmasMin * sqrt(var));
        MinRate[i+1] = Maths::min(MinRate[i+1], 0.0); /* to avoid cutting at positive rate */

    }
}


vector<double> SRMBasisSpreadUtil::computeSqrtYearFrac() const
{
    vector<double> sqrtYearFracs((*dates).size()-1);
    for (size_t i = 0; i < sqrtYearFracs.size(); ++i){
        sqrtYearFracs[i] = sqrt(SRMYearFrac((*dates)[i], (*dates)[i+1]));
    }
    return sqrtYearFracs;
}

double SRMBasisSpreadUtil::hFactor(const DateTime& d1, const DateTime& d2) const
{
    return exp(-beta*SRMYearFrac(d1, d2));
}

vector<double> SRMBasisSpreadUtil::computePartialH(const DateTime& today, const DateTimeArray& forwardDates) const
{
    vector<double> partialH(forwardDates.size());
    for(size_t i=0; i<partialH.size(); ++i)
        partialH[i] = hFactor(today,forwardDates[i]);
    return partialH;
}

vector<double> SRMBasisSpreadUtil::computeFwdParSpreads(const DateTimeArray& forwardDates) const
{
    //vector<double> spotSpreads(forwardDates.size()-1);
    vector<double> spotSpreads(forwardDates.size());
    for(size_t i=0; i<spotSpreads.size(); ++i)
        spotSpreads[i] = calcFwdSpread(forwardDates[i]);
    return spotSpreads;
}

double SRMBasisSpreadUtil::calcFwdSpread(const DateTime& resetDate) const
{
    //return basisCurve->parSpread(std::max(resetDate, basisCurve->getFirstCurveDate())); // never before the 1st curve date
    DateTime payDate = basisCurve->getRefPaymentDate(resetDate);
    return basisCurve->fwdSpread(std::max(resetDate, basisCurve->getFirstCurveDate()), payDate);
}


