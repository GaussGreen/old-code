//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditHJMUtil.cpp
//
//   Description : Derived SRMCreditUtil class for Markovian HJM / RS model
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMCreditHJMUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/CRCalib.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
SRMCreditHJMUtil::SRMCreditHJMUtil(
    const DateTime&         baseDate,
    const string&           smileParamsKey,
    ICDSParSpreadsConstSP   stochCDSCurve):
        SRMCreditUtil(baseDate, 
                  smileParamsKey, 
                  stochCDSCurve)
{
    static const string method("SRMCreditHJMUtil::SRMCreditHJMUtil");
    try {
        // back out smile parameters
        CRCalib::SmileRequest smileRequest(smileParamsKey);
        CVolProcessed* vol = stochCDSCurve->getProcessedVol(&smileRequest);
        smartPtr<CRCalib::VolProcessed> volData(
            &dynamic_cast<CRCalib::VolProcessed&>(*vol));
        qLeft = volData->getQLeft();
        qRight = volData->getQRight();
        fwdShift = volData->getFwdShift(); 
        if (Maths::equals(fwdShift, -1.0)) {
            throw ModelException(method,"Pivot ratio is 0");
        } 
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** destructor */
SRMCreditHJMUtil::~SRMCreditHJMUtil() {
};

/** finishes initialization */
// called when we know allDates
void SRMCreditHJMUtil::setTimeLine(DateTimeArrayConstSP simDates) {
    static const string method("SRMCreditHJMUtil::setTimeLine");
    if (initialized) 
        return; // nothing to do

    try {
        initialized = true;
    
        dates = simDates;

        // back out model parameters (in CR, there is only one model parameter, namely beta)
        MRSpotVolRequest betaRequest;
        CVolProcessed* volBeta = stochCDSCurve->getProcessedVol(&betaRequest);
        MRSpotVolProcessedSP volBetaData(
            &dynamic_cast<MRSpotVolProcessed&>(*volBeta));

        beta = volBetaData->meanReversion(); 

        calcExtendedTimeLine(); // get 'extended' timeline
    
        // populate SpotVol (so far we only have FlatCDSSpotVol ...)
        DoubleArray spotVolTemp(extendedTimeLine.size());
        volBetaData->spotVol(baseDate, extendedTimeLine, spotVolTemp);
        SpotVol = vector<double>(spotVolTemp.begin(), spotVolTemp.end());
        calcFwdIntensity();
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** computes 'KFactors' on srmUtil 'dates' 
    kFactor param is initialised to correct length. */
void SRMCreditHJMUtil::computeKFactor(vector<double>& kFactor) const {
    kFactor.resize(dates->size()-1);
    double lowerLBar = lBar(0);
    for (unsigned int i = 0; i < kFactor.size(); i++) {
        if (Maths::isZero(lowerLBar)) {
            throw ModelException("SRMCreditHJMUtil::computeKFactor",
                                 "Zero forward rate at "+
                                 (*dates)[i].toString());
        }
        double higherLBar = lBar(i+1);
        double ratio = higherLBar/lowerLBar;
        double del_t = (*dates)[i].yearFrac((*dates)[i+1]);
        kFactor[i] = ratio * exp(-beta * del_t);
        lowerLBar = higherLBar; // for next time
    }
}

/** computes 'GFactors' on srmUtil 'dates' 
    gFactor param  is initialised to correct length. */
void SRMCreditHJMUtil::computeGFactor(vector<double>& gFactor) const {
    gFactor.resize(dates->size()-1);
    for (unsigned int i = 0; i < gFactor.size(); i++) {
        double r_not = lBar(i);
        if (Maths::isZero(r_not)) {
            throw ModelException("SRMCreditHJMUtil::computeGFactor",
                                 "Zero forward rate at "+
                                 extendedTimeLine[i].toString());
        }
        double tFrom0 = dates->front().yearFrac((*dates)[i]);
        double denom = r_not * exp(-beta * tFrom0);
        double g = GfactorCR[i+1] - GfactorCR[i];
        gFactor[i]= g/denom;
    }
}

/** returns GfactorCR */
const vector<double>& SRMCreditHJMUtil::getGfactorCR() const { 
    assert(initialized);
    return GfactorCR;
}

/** returns the spot vols */
const vector<double>& SRMCreditHJMUtil::getSpotVols() const {
    assert(initialized);
    return SpotVol;
}

/** returns beta */
double SRMCreditHJMUtil::getBeta() const {
    assert(initialized);
    return beta;
}

/** returns qRight */
double SRMCreditHJMUtil::getQRight() const {
    return qRight;
}

/** returns qLeft */
double SRMCreditHJMUtil::getQLeft() const {
    return qLeft;
}

/** returns fwdShift */
double SRMCreditHJMUtil::getFwdShift() const {
    return fwdShift;
}

/** for calculating 'lbar' when moving FROM a date in dates whose index
    is given by futureDatesIndex */
double SRMCreditHJMUtil::lBar(int futureDatesIndex) const {
    assert(initialized);
    return FwdIntensity[futureDatesIndex];
}

/** for calculating 'Lbar' on supplied date */
double SRMCreditHJMUtil::lBar(const DateTime& date) const {
    int index = date.findLower(extendedTimeLine);
    if (index < 0) {
        throw ModelException("SRMCreditHJMUtil::lBar", "date is before today");
    }
    return FwdIntensity[index]; // FwdIntensity[index] = lambda(0,t), where t = date, computation in calcFwdIntensity()
}

/** Calculates and populates the GfactorCR[].
    This array is used in Gfactor() to streamline the calculation. From
    irdiffuse::PopulateGfactorCR */
void SRMCreditHJMUtil::populateGfactorCR() {
    GfactorCR.resize(extendedTimeLine.size());

    double day1 = extendedTimeLine.front().yearFrac(extendedTimeLine.back());
    double day_diff;
    for (int k = extendedTimeLine.size()-2; k >= 0; k--) {
        day_diff = day1;
        day1 = extendedTimeLine[0].yearFrac(extendedTimeLine[k]);
        day_diff -= day1;
        double lBarValue = lBar(extendedTimeLine[k]);
        
        GfactorCR[k] = GfactorCR[k+1] - SRMUtil::GFAC(day1, day_diff, beta) * lBarValue;
    }
}

/** this function is called on the edfForwardDates dates */
void SRMCreditHJMUtil::populatePartialIntCR(const DateTime& today,
                                     const DateTimeArray& myDatesIn, 
                                     vector<double>& partialIntegralOut)
{
    if (myDatesIn.empty())
        return;
    partialIntegralOut = vector<double>(myDatesIn.size());

    // FIXME: see if need to merge with RatesUtil's extended timeline
    DateTimeArray myDates = DateTime::merge(myDatesIn, extendedTimeLine); 
    vector<int> positions = DateTime::getIndexes(myDates, myDatesIn);

    double day1 = today.yearFrac(myDates.back());

    vector<double> pI(myDates.size());

    for (int k = myDates.size()-2; k >= 0; k--) {
        double day_diff = day1;
        day1 = today.yearFrac(myDates[k]);
        day_diff -= day1;
        double lBarValue = lBar(myDates[k]);
        pI[k] = pI[k+1]- SRMUtil::GFAC(day1, day_diff, beta) * lBarValue;
    }

    for(int i=0; i<myDatesIn.size(); ++i) {
        partialIntegralOut[i]=pI[positions[i]];
    }
}

/** this function is called on the requestedEDFdates */
void SRMCreditHJMUtil::populatePartialZeta(const DateTime& today,
                                    const DateTimeArray& myDates,
                                    vector<double>& zeta) 
{    
    zeta = vector<double>(myDates.size());
    
    for (int k = 0; k < myDates.size(); k++) {
        double day1 = today.yearFrac(myDates[k]);
        double lBarValue = lBar(myDates[k]);
        zeta[k] = exp(beta * day1) / lBarValue;
    }
}
   
/** populates FwdIntensity field on extendedTimeLine.
    requires 'dates' field */
void SRMCreditHJMUtil::calcFwdIntensity() {
    // use 1Y offset for last FwdIntensity
    DateTime endDate(MaturityPeriod::toDate(1,"A",extendedTimeLine.back())); // MaturityPeriod::toDate(count, interval (A=annual), aDate)
    int numDates = extendedTimeLine.size();
    FwdIntensity.resize(numDates);
    DoubleArray logOfRatio(numDates);
    DoubleArray dts(numDates);
    int i = 0;
    for (i = 0; i < numDates; i++) {
        const DateTime& start = extendedTimeLine[i];
        const DateTime& end = i==numDates-1? endDate: extendedTimeLine[i+1];        
        dts[i] = SRMYearFrac(start, end);
        logOfRatio[i] = logSurvProbRatio(start, end);
    }
    i = 0;
    while (i < numDates) {
        int lastDate = extendedTimeLine[i].getDate();
        double totalLogOfRatio = logOfRatio[i];
        double del_t = dts[i];
        int j = i + 1;
        while (j < numDates 
               && (lastDate == extendedTimeLine[j].getDate() || !Maths::isPositive(del_t))){
            totalLogOfRatio += logOfRatio[j];
            del_t += dts[j];
            ++j;
        }
        double lbar = - totalLogOfRatio / del_t;
        int start_i = i;
        i = j;
        while ((--j) >= start_i) {
            FwdIntensity[j] = lbar;
        }
    }
    populateGfactorCR(); 
}

/** computes KFactor between the 2 dates supplied */
double SRMCreditHJMUtil::kFactor(
    const DateTime& dateFrom,
    const DateTime& dateTo) const 
{
    double denominator = lBar(dateFrom);
    if (Maths::isZero(denominator)) {
        throw ModelException("SRMCreditHJMUtil::kFactor", "Zero forward intensity at "+dateFrom.toString());
    }
    double ratio = lBar(dateTo) / denominator;
    double del_t = dateFrom.yearFrac(dateTo);    
    return ratio*exp(-beta*del_t);
}

/** computes GFactor between the 2 dates supplied */
double SRMCreditHJMUtil::gFactor(const DateTime& dateFrom,
                                 const DateTime& dateTo) const 
{
    int idxM = dateFrom.findUpper(extendedTimeLine);
    /* if both date1 and date0 > FwdDate[NbFwdDates - 1] (last fwd
     * date), a flat extrapolation of the curve is done */
    if (idxM == extendedTimeLine.size()) {
        double del_t = dateFrom.yearFrac(dateTo);
        return SRMUtil::GFAC(0, del_t, beta); 
    } 

    int idxN = dateTo.findLower(extendedTimeLine);           
    if (idxM > idxN) {
        // both dateFrom and dateTo lie in the same interval
        double del_t = dateFrom.yearFrac(dateTo);
        return SRMUtil::GFAC(0, del_t, beta);
    }

    double r_not = lBar(dateFrom);
    if (Maths::isZero(r_not)) {
        throw ModelException("SRMCreditHJMUtil::gFactor", "Zero forward rate at "+
                             dateFrom.toString());
    }

    double tFrom0 = extendedTimeLine.front().yearFrac(dateFrom);
    double g;
    if (idxN > idxM) {
        g = GfactorCR[idxN] - GfactorCR[idxM];
        double denom = r_not * exp(-beta * tFrom0);
        g /= denom;
    } else {
        g = 0.0;
    }

    if (dateFrom != extendedTimeLine[idxM]) {
        double del_t = dateFrom.yearFrac(extendedTimeLine[idxM]);
        g += SRMUtil::GFAC(0, del_t, beta);
    }
        
    if (dateTo != extendedTimeLine[idxN]) {
        double tFrom0 = dateFrom.yearFrac(extendedTimeLine[idxN]);
        double del_t  = extendedTimeLine[idxN].yearFrac(dateTo);
        double r  = lBar(extendedTimeLine[idxN]);
        g += r * SRMUtil::GFAC(tFrom0, del_t, beta) / r_not ;
    }
    return g;                               
}

/** Utility routine to calculate the cutoff rate and populate
    the cutoff rate arrays */
void SRMCreditHJMUtil::calcEffRateLimit(
    double          NbSigmasMax,    // Number or sigmas to cut at
    double          NbSigmasMin,    // Number or sigmas to cut at
    vector<double>& MaxRate,        // Cutoff forward rates 
    vector<double>& MinRate) const  // Cutoff forward rates 
{
    // TO DO - Fix +10
    calcEffRateLimitCR(NbSigmasMax+10.0, NbSigmasMin+10.0, MaxRate, MinRate);
}

/** Utility routine to calculate the cutoff rate and populate the cutoff rate arrays 
    see crdiffuse::CalcEffRateLimitCR for details */
void SRMCreditHJMUtil::calcEffRateLimitCR(
    double           NbSigmasMax,    // (I) Number or sigmas to cut at
    double           NbSigmasMin,    // (I) Number or sigmas to cut at
    vector<double>&  MaxRate,        // (O) Cutoff forward rates
    vector<double>&  MinRate) const  // (O) Cutoff forward rates    
{
    MaxRate.resize(extendedTimeLine.size());
    MinRate.resize(extendedTimeLine.size());

    MaxRate[0] = FwdIntensity[0];
    MinRate[0] = Maths::min(FwdIntensity[0], 0.0);
    double var = 0.0;

    for (unsigned int i = 0; i < MinRate.size()-1; i++) {
        var += extendedTimeLine[i].yearFrac(extendedTimeLine[i+1]) * SpotVol[i] * SpotVol[i];
        
        /* lognormal MAX-Cutoff */
        double dummy = NbSigmasMax * sqrt(var);
        dummy = Maths::min(dummy, 100.0); /* to avoid blow-up in exp() */
        MaxRate[i+1] = FwdIntensity[i+1] * exp(dummy);

        /* normal MIN-Cutoff */
        MinRate[i+1] = FwdIntensity[i+1] * (1 - NbSigmasMin * sqrt(var));
        MinRate[i+1] = Maths::min(MinRate[i+1], 0.0);/* to avoid cutting at positive rate */
    }
}

DRLIB_END_NAMESPACE
