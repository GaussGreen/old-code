//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMRatesHJMUtil.cpp (from SRMRatesUtil)
//
//   Description : Helper for Ritchken-Sankarasubramanian model
//
//   Author      : Henrik Rasmussen (from earlier file from Mark Robson)
//
//   Date        : 2 May 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/SRMSwaption.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include <cassert>

DRLIB_BEGIN_NAMESPACE


///// constructs and populates SRMRatesUtil object
SRMRatesHJMUtil::SRMRatesHJMUtil(
    const DateTime&      baseDate,
    int                  numFactors,
    const string&        modelParamsKey,
    const string&        smileParamsKey,
    CVolProcessedSP      _processedVol,
    IYieldCurveConstSP   discYC,
    IYieldCurveConstSP   diffYC,
    bool                 _skipFlag,
    double               _flatVolIr,
    const string&        cutoffChoice,
    double               constCutoffValue,
    const string&        corrSwapStart, // eg 1Y  (offset to yc spot date)
    const string&        corrSwapMat,   // eg 10Y (offset to start)
    const string&        corrSwapDCC,   // eg Act/365F
    const string&        corrSwapFreq): // eg 6M
    SRMRatesUtil(
      baseDate,
      numFactors,
      modelParamsKey,
      smileParamsKey,
      _processedVol,
      discYC,
      diffYC,
      _skipFlag,
      _flatVolIr,
      cutoffChoice,
      constCutoffValue,
      corrSwapStart, // eg 1Y  (offset to yc spot date)
      corrSwapMat,   // eg 10Y (offset to start)
      corrSwapDCC,   // eg Act/365F
      corrSwapFreq)
{
    static const string method("SRMRatesHJMUtil::SRMRatesHJMUtil");
    try{
        // get alpha, beta, qLeft, qRight
        IRCalib::SmileRequest smileRequest(smileParamsKey);
        CVolProcessed* vol = diffYC->getProcessedVol(&smileRequest);
        smartPtr<IRCalib::VolProcessed> volData(
            &dynamic_cast<IRCalib::VolProcessed&>(*vol));
        const DoubleArray& smileParams = volData->getParams();
        if (smileParams.size() < 3){
            // internal error
            throw ModelException(method, "Number of ir vol smile params wrong");
        }
        qLeft = 1.0 - smileParams[0]; /* I guess somebody changed the meaning
                                         of q Left/Right at some point! */
        qRight = 1.0 - smileParams[1];
        fwdShift = smileParams[2];
        if (Maths::equals(fwdShift, -1.0)){
            throw ModelException(method, 
                                 "Pivot ratio is 0 for "+diffYC->getCcy());
        }
        IRCalib::ModelRequest modelRequest(modelParamsKey);
        vol = diffYC->getProcessedVol(&modelRequest);
        volData.reset(&dynamic_cast<IRCalib::VolProcessed&>(*vol));
        const DoubleArray& modelParams = volData->getParams();
        switch (numFactors){
        case 1:
            // one factor
            model = SRMRatesFactorModelSP(new SRMRatesFactorModel1F());
            model->factors[0].beta = modelParams[0];
            model->factors[0].alpha = modelParams[1];
            break;
        case 2:
            model = SRMRatesFactorModelSP(new SRMRatesFactorModel2F());
            model->factors[0].beta = modelParams[0];
            model->factors[1].beta = modelParams[1];
            model->factors[0].alpha = modelParams[2];
            model->factors[1].alpha = modelParams[3];
            model->rho[0] = modelParams[4];
            break;
        case 3:
            model = SRMRatesFactorModelSP(new SRMRatesFactorModel3F());
            model->factors[0].beta = modelParams[0];
            model->factors[1].beta = modelParams[1];
            model->factors[2].beta = modelParams[2];
            model->factors[0].alpha = modelParams[3];
            model->factors[1].alpha = modelParams[4];
            model->factors[2].alpha = modelParams[5];
            model->rho[0] = modelParams[6];
            model->rho[1] = modelParams[7];
            model->rho[2] = modelParams[8];
            break;
        default:
            throw ModelException(method, Format::toString(numFactors)+" factor "
                                 "IR model not supported");
        }
       
    } catch (exception& e){
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
}

/** Constructs and populates SRMRatesUtil object suitable using the 'VF' style
    calibration. qLeft and qRight are hard coded to 1, as is alpha. Note
    that only single factor IR is supported */
SRMRatesHJMUtil::SRMRatesHJMUtil(
    const DateTime&      baseDate,
    FXAssetConstSP       fx, /* optional - if specified will add FX smile
                                dates to timeline */
    const string&        modelParamsKey,
    VolProcessedBSIRSP   processedVol,
    IYieldCurveConstSP   discYC,
    IYieldCurveConstSP   diffYC,
    bool                 skipFlag,
    const string&        _corrSwapStart, // eg 1Y  (offset to today)
    const string&        _corrSwapMat,   // eg 10Y (offset to start)
    const string&        _corrSwapDCC,   // eg Act/365F
    const string&        _corrSwapFreq): // eg 6M
    SRMRatesUtil(
       baseDate,
       fx, /* optional - if specified will add FX smile
	   	      dates to timeline */
	   modelParamsKey,
	   processedVol,
	   discYC,
	   diffYC,
	   skipFlag,
	   _corrSwapStart, // eg 1Y  (offset to today)
	   _corrSwapMat,   // eg 10Y (offset to start)
	   _corrSwapDCC,   // eg Act/365F
	   _corrSwapFreq),
       qLeft(1.0), 
       qRight(1.0), 
       fwdShift(0.0) // hard coded values for q
{
    static const string method("SRMRatesHJMUtil::SRMRatesHJMUtil");
    try{
        if (!processedVol){
            throw ModelException(method, "IR Swaption vols must be available"
                                 " for VF style calibration");
        }
        // get hold of swaptionExpiries, swapStartDates, swapMatDates, vols
        processedVol->getBMDetails(swaptionExpiries, swapStartDates, 
                                   swapMatDates, swaptionVols);
        // nasty hack to avoid code falling over due to duplicate points per
        // date or EOD point followed by SOD point
        DateTime::setTimeOfDay(swaptionExpiries, baseDate.getTime());
        if (swaptionExpiries.front().equals(baseDate)){
            throw ModelException(method, "Swaption expiry on same day as "
                                 "today not supported");
        }
        // do the other two for safety
        DateTime::setTimeOfDay(swapStartDates, baseDate.getTime());
        DateTime::setTimeOfDay(swapMatDates, baseDate.getTime());
        swapFrequency = processedVol->getSwapFrequency();
        swapDCC = processedVol->getSwapDCC();

        IRCalib::ModelRequest modelRequest(modelParamsKey);
        CVolProcessed* vol = diffYC->getProcessedVol(&modelRequest);
        smartPtr<IRCalib::VolProcessed> volData(
            &dynamic_cast<IRCalib::VolProcessed&>(*vol));
        const DoubleArray& modelParams = volData->getParams();
        // one factor
        model = SRMRatesFactorModelSP(new SRMRatesFactorModel1F());
        model->factors[0].beta = modelParams[0];
        // shouldn't really need to hard code this to 1.0 I think as it should
        // cancel out - but it doesn't.
        model->factors[0].alpha = 1.0; // modelParams[1];
        // use today + swaption dates + input dates for our timeline
        dates = DateTimeArrayConstSP(new DateTimeArray(1, baseDate));
        DateTimeArray smileDates;
        if (fx.get()){
            // add fx smile dates to timeline
            smileDates = CriticalDateCollector::collectVolDates(
                fx, FXVolBase::TYPE, baseDate, swaptionExpiries.back());
            DateTime::setTimeOfDay(smileDates, baseDate.getTime());
        }

        vector<const DateTimeArray*> datesToMerge(3);
        datesToMerge[0] = &(*dates);
        datesToMerge[1] = &swaptionExpiries;
        datesToMerge[2] = &smileDates;
        dates = DateTimeArrayConstSP(new DateTimeArray(DateTime::merge(datesToMerge)));
        calcExtendedTimeLine(); // get 'extended' timeline
        // and then as for swaptionExpiries adjust time of day
        DateTime::setTimeOfDay(extendedTimeLine, baseDate.getTime());
        // and remove duplicates (could possibly be one)
        DateTime::removeDuplicates(extendedTimeLine, true);
        calcFwdRateAndRBar();

        // populate LastDate - first find benchmark on/after last sim date
        int bmIdx = (*dates).back().findUpper(swaptionExpiries);
        if (bmIdx == swaptionExpiries.size()){
            bmIdx--;
        }
        LastDate = swapMatDates[bmIdx].max((*dates).back());
        // calculate spot vol along swaption maturities using 'VF' methodology
        swaptionSpotVol = SRMSwaption::spotVolVF(*this, skipFlag);
        extendSpotVol(swaptionSpotVol);
    } catch (exception& e){
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
}

/** number of factors in model */
int SRMRatesHJMUtil::numFactors() const
{
    return model->factors.size();
}

void  SRMRatesHJMUtil::setTimeLine(DateTimeArrayConstSP simDates) // called when we know allDates
{
    static const string method("SRMRatesUtil::setTimeLine");
    if (initialized) {
        if (0 && (simDates->size() != dates->size() || ! DateTime::isSubset(*simDates, *dates)))
            throw ModelException(method, "Re-initialized with a different timeline");
        return;
    }
    initialized = true;
    try {
        dates = simDates;
        calcExtendedTimeLine(); // get 'extended' timeline
        calcFwdRateAndRBar();
        
//        if (!processedVol){
        if (! processedVol){
            // this wasn't in the original plan but it seems you can do stuff
            // even if you don't have the swaption vols
            spotVol(flatVolIr); // populate SpotVol with 1.0 at each point
            return;
        }
        else {
            if (VolProcessedBSIR::TYPE->isAssignableFrom(processedVol->getClass()))
            {
                VolProcessedBSIRSP processedVol = VolProcessedBSIRSP::dynamicCast(this->processedVol);
                swapFrequency = processedVol->getSwapFrequency();
                swapDCC = processedVol->getSwapDCC();
                // get hold of swaptionExpiries, swapStartDates, swapMatDates, vols
                processedVol->getBMDetails(swaptionExpiries, swapStartDates,
                                            swapMatDates, swaptionVols);
        
                // populate LastDate - first find benchmark on/after last sim date
                int bmIdx = (*dates).back().findUpper(swaptionExpiries);
                if (bmIdx == swaptionExpiries.size()){
                    bmIdx--;
                }
                LastDate = swapMatDates[bmIdx].max((*dates).back());
                spotVolBM(skipFlag); // calculate SpotVol and tau
            }
            else if (IRVolSpot::IRVolSpotProcessed::TYPE->isAssignableFrom(processedVol->getClass()))
            {
                // the vol is already calibrated:
                IRVolSpot::IRVolSpotProcessed* processedVol = dynamic_cast<IRVolSpot::IRVolSpotProcessed*>(this->processedVol.get());
                swaptionExpiries = processedVol->getSpotExpiryDates();
                swaptionSpotVol  = processedVol->getSpotVols();

                swaptionSpotVolAtSimDates.resize(swaptionSpotVol.getLength());
                SpotVol.resize(swaptionSpotVol.getLength());
                for (size_t i=0; i<SpotVol.size(); ++i)
                {
                    SpotVol[i] = swaptionSpotVol[i];
                    swaptionSpotVolAtSimDates[i] = swaptionSpotVol[i];
                }

                // hack here -- need to make it either automatically or user configurable
                // like in SimPathI (a.k->a. Driver Rate parameter)
                tau = getExtendedTimeLine();
                for(int i=0; i<tau.size(); ++i)
                    tau[i] = min(tau[i].rollDate(3650), tau.back());
                /* Interpolate spot vol throughout the timeline */
                extendSpotVol(SpotVol);

                // VERIFY THIS: also need to extend swaptionSpotVol and swaptionSpotVolAtSimDates
                vector<double> swaptionSpotVolNew(swaptionSpotVol.begin(), swaptionSpotVol.end());
                swaptionSpotVolAtSimDates = SRMUtil::extendVol(swaptionExpiries,
                    swaptionSpotVolNew, 
                    (*dates));

                /* Interpolate spot vol throughout the timeline */
                extendSpotVol(swaptionSpotVol);
            }
        }
        
    } catch (exception& e){
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
  
}

/** vol param */
double SRMRatesHJMUtil::getAlpha(int factorIdx) const
{ 
    return model->factors[factorIdx].alpha;
}

void SRMRatesHJMUtil::getAlpha(vector<double>& alphaByFactor) const
{
    alphaByFactor.resize(model->factors.size());
    for (unsigned int i = 0; i < alphaByFactor.size(); i++){
        alphaByFactor[i] = model->factors[i].alpha;
    }
}

void SRMRatesHJMUtil::getBeta(vector<double>& betaByFactor) const
{
    betaByFactor.resize(model->factors.size());
    for (unsigned int i = 0; i < betaByFactor.size(); i++){
        betaByFactor[i] = model->factors[i].beta;
    }
}

/** Returns beta for specified factor */
double SRMRatesHJMUtil::getBeta(int factorIdx) const
{
    return model->factors[factorIdx].beta;
}

/** Returns rho[rhoIndex] where rho are model parameters */
double SRMRatesHJMUtil::getRho(int rhoIndex) const
{
    return model->rho[rhoIndex];
}

/** Returns the model rho parameters - length corresponds to off diagonal
    elements of a symmetric matrix (size of which is number of factors) */
const vector<double>& SRMRatesHJMUtil::getRho() const
{
    return model->rho;
}

/** Returns 'variance' as determined by model parameters */
double SRMRatesHJMUtil::modelVariance(double delt_1, double delt_2) const
{
    return model->variance(delt_1, delt_2);
}

/* Exponential decay function for specified factor */
double  SRMRatesHJMUtil::expDecay(int factor, double t) const
{
    return model->factors[factor].expDecay(t);
}

/** Returns 'mean reversion integral' as determined by model parameters */
void SRMRatesHJMUtil::modelMeanReversionIntegral(
    double          delt_1, 
    double          delt_2,
    vector<double>& meanReversionIntegrals) const //O
{
    model->meanReversionIntegral(delt_1, delt_2, meanReversionIntegrals);
}
   
/** computes 'KFactors' on srmUtil 'dates' - kFactor field param  is
    initialised to correct length. From irdiffuse::Kfactor */
void SRMRatesHJMUtil::computeKFactor(vector<double>& kFactor, int factorIdx) const
{
    kFactor.resize((*dates).size()-1);
    double beta = model->factors[factorIdx].beta;
    double lowerRBar = rBar(0);
    for (unsigned int i = 0; i < kFactor.size(); i++){
        if (Maths::isZero(lowerRBar)){
            throw ModelException("SRMRatesHJMUtil::computeKFactor",
                                 "Zero forward rate at "+
                                 (*dates)[i].toString());
        }
        double higherRBar = rBar(i+1);
        double ratio = higherRBar/lowerRBar;
        double del_t = SRMYearFrac((*dates)[i], (*dates)[i+1]);
        kFactor[i] = ratio * exp(-beta * del_t);
        lowerRBar = higherRBar; // for next time
    }
}
/** computes 'GFactors' on srmUtil 'dates' - kFactor param  is
    initialised to correct length. From a subset of irdiffuse::Gfactor */
void SRMRatesHJMUtil::computeGFactor(vector<double>& gFactor, int factorIdx) const
{
    gFactor.resize((*dates).size()-1);
    const vector<double>& GfactorArrayIR = 
        model->factors[factorIdx].GfactorArrayIR;
    double beta = model->factors[factorIdx].beta;
    for (unsigned int i = 0; i < gFactor.size(); i++){
        double r_not = rBar(i);
        if (Maths::isZero(r_not)){
            throw ModelException("SRMRatesHJMUtil::computeGFactor",
                                 "Zero forward rate at "+
                                 extendedTimeLine[i].toString());
        }
        double tFrom0 = SRMYearFrac((*dates).front(), (*dates)[i]);
        double denom = r_not * exp(-beta * tFrom0);
        double g = GfactorArrayIR[i+1] - GfactorArrayIR[i];
        gFactor[i]= g/denom;
    }
}

/** Calculates and populates the GfactorArrayIR[].
    This array is used in Gfactor() to streamline  the calculation. From
    irdiffuse::PopulateGfactorArrayIR */
void SRMRatesHJMUtil::populateGfactorArrayIR()
{

    int numFactors = this->numFactors();
    for (int i = 0; i < numFactors; i++){
        model->factors[i].GfactorArrayIR.resize(extendedTimeLine.size());
    }
    double day1 = SRMYearFrac(extendedTimeLine.front(), extendedTimeLine.back());
    double day_diff;
    for (int k = extendedTimeLine.size()-2; k >= 0; k--) {
        day_diff = day1;
        day1 = SRMYearFrac(extendedTimeLine[0], extendedTimeLine[k]);
        day_diff -= day1;
        double rBarValue = rBar(extendedTimeLine[k]);
        for (int factorIdx = 0; factorIdx < numFactors; factorIdx++){
            vector<double>& GfactorArrayIR =
                model->factors[factorIdx].GfactorArrayIR;
            double beta = model->factors[factorIdx].beta;
            GfactorArrayIR[k] = GfactorArrayIR[k+1] 
                - SRMUtil::GFAC(day1, day_diff, beta) * rBarValue;
        }
    }
}

// this function is called on the edfForwardDates dates
void SRMRatesHJMUtil::populatePartialIntIR(const DateTime& today, 
                                     const DateTimeArray& myDatesIn, 
                                     vector<double> partialIntegralOut[])
{
    if (myDatesIn.empty())
    	return;
    int numFactors = this->numFactors();

    DateTimeArray myDates = DateTime::merge(myDatesIn, extendedTimeLine); 
    vector<int> positions = DateTime::getIndexes(myDates, myDatesIn);
    vector<vector<double> > partialIntegral(numFactors, vector<double>(myDates.size(), 0));

    double day1 = SRMYearFrac(today, myDates.back());
    
    for (int k = myDates.size()-2; k >= 0; k--) {
        double day_diff = day1;
        day1 = SRMYearFrac(today, myDates[k]);
        day_diff -= day1;
        double rBarValue = rBar(myDates[k]);
        for (int factorIdx = 0; factorIdx < numFactors; factorIdx++){
            vector<double>& pI = partialIntegral[factorIdx];
            double beta = model->factors[factorIdx].beta;
            pI[k] = pI[k+1]- SRMUtil::GFAC(day1, day_diff, beta) * rBarValue;
        }
    }

    for(int factorIdx=0; factorIdx < numFactors; ++factorIdx)
    { 
        partialIntegralOut[factorIdx].resize(myDatesIn.size());
        for(int i=0; i<myDatesIn.size(); ++i)
        {
            partialIntegralOut[factorIdx][i]=partialIntegral[factorIdx][positions[i]];
        }
    }

}

// this function is called on the requestedEDFdates
void SRMRatesHJMUtil::populatePartialZeta(const DateTime& today, 
									 const DateTimeArray& myDates, 
									 vector<double> zeta[])
{
    int numFactors = this->numFactors();
    
    for (int i = 0; i < numFactors; i++){
        zeta[i] = vector<double>(myDates.size());
    }
    
    for (int k = 0; k < myDates.size(); k++) {
        double day1 = SRMYearFrac(today, myDates[k]);
        double rBarValue = rBar(myDates[k]);
        for (int factorIdx = 0; factorIdx < numFactors; factorIdx++){
            double beta = model->factors[factorIdx].beta;
            zeta[factorIdx][k] = exp(beta * day1) / rBarValue;
        }
    }
}

/** Populates FwdRate and RBar fields. Requires calcExtendedTimeLine() to have
    been called.
    NB: This is a departure from 'native' srm3. This is so multiple time points
    in a day are supported. The definition of 'rBar' is Hitier's ie rBar is the 
    piecewise constant approximation of 'f' (fwdRate). The set of time points is 
    chosen to be that subset of extendedTimeLine such that every point in that
    subset belongs to a different day. No differences are expected for ported,
    native srm3 tests. */
void SRMRatesHJMUtil::calcFwdRateAndRBar()
{
    auto_ptr<IYieldCurve::IKey> key(diffYC->logOfDiscFactorKey());
    // use 1Y offset for last FwdRate
    DateTime endDate(MaturityPeriod::toDate(1,"A",extendedTimeLine.back()));
    int numDates = extendedTimeLine.size();
    FwdRate.resize(numDates);
    RBar.resize(numDates);
    DoubleArray logOfRatio(numDates);
    DoubleArray dts(numDates);
    int i = 0;
    for (; i < numDates; i++){
        const DateTime& start = extendedTimeLine[i];
        const DateTime& end = i==numDates-1? endDate: extendedTimeLine[i+1];
        dts[i] = SRMYearFrac(start, end);
        logOfRatio[i] = key->calc(start, end);
    }
    i = 0;
    while (i < numDates){
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
        double rbar = - totalLogOfRatio / del_t;
        int start_i = i;
        i = j;
        while ((--j) >= start_i){
            RBar[j] = rbar;
            FwdRate[j] = rbar;
        }
    }
    populateGfactorArrayIR();
}


double SRMRatesHJMUtil::rBar(int futureDatesIndex) const
{
	if (initialized == false) 
	{
        throw ModelException("Trying to call rBar(int) before initialized");
	}

	return RBar[futureDatesIndex]; 
}

// for calculating 'Rbar' on supplied date
double SRMRatesHJMUtil::rBar(const DateTime& date) const
{
    int index = date.findLower(extendedTimeLine);
    if (index < 0){
        throw ModelException("SRMRatesHJMUtil::rBar", 
                             "date is before today");
    }
    return RBar[index];
}

/** computes KFactor between the 2 dates supplied. From irdiffuse::Kfactor */
void SRMRatesHJMUtil::kFactor(const DateTime& dateFrom,
                        const DateTime& dateTo,
                        vector<double>& k) const // (0)
{
    double denominator = rBar(dateFrom);
    if (Maths::isZero(denominator)){
        throw ModelException("SRMRatesHJMUtil::kFactor", "Zero forward rate at "+
                             dateFrom.toString());
    }
    double ratio = rBar(dateTo) / denominator;
    double del_t = SRMYearFrac(dateFrom, dateTo);
    for (unsigned int i = 0; i < k.size(); i++){
        double beta = model->factors[i].beta;
        k[i] = ratio*exp(-beta*del_t);
    }
}
/** computes KFactor between the 2 dates supplied  for specified model
    factor index */
double SRMRatesHJMUtil::kFactor(const DateTime& dateFrom, 
                          const DateTime& dateTo,
                          int             factorIdx) const
{
    vector<double> kFac(factorIdx+1);
    kFactor(dateFrom, dateTo, kFac);
    return kFac[factorIdx];
}


/** computes GFactor between the 2 dates supplied. From irdiffuse::Gfactor */
void SRMRatesHJMUtil::gFactor(const DateTime& dateFrom,
                        const DateTime& dateTo,
                        vector<double>& g) const // (0)
{
    int idxM = dateFrom.findUpper(extendedTimeLine);
    /* if both date1 and date0 > FwdDate[NbFwdDates - 1] (last fwd
     * date), a flat extrapolation of the curve is done */
    if (idxM == extendedTimeLine.size()) {
        double del_t = SRMYearFrac(dateFrom, dateTo);
        for (unsigned int i = 0; i < g.size(); i++){
            double beta = model->factors[i].beta;
            g[i] = (1.0 / beta * (1.0 - exp(-beta * del_t)));   
            // think it should be (1.0 / (beta * (1.0 - exp(-beta * del_t))));
            // or even better: SRMUtil::GFAC(0, del_t, beta);
        }
        return;
    } 
    int idxN = dateTo.findLower(extendedTimeLine);   
            
    if (idxM > idxN) {
        // both dateFrom and dateTo lie in the same interval
        double del_t = SRMYearFrac(dateFrom, dateTo);
        for (unsigned int i = 0; i < g.size(); i++){
            double beta = model->factors[i].beta;
            g[i] = SRMUtil::GFAC(0, del_t, beta);
        }
        return;
    }
    double r_not = rBar(dateFrom);
    if (Maths::isZero(r_not)){
        throw ModelException("SRMRatesHJMUtil::gFactor", "Zero forward rate at "+
                             dateFrom.toString());
    }
    double tFrom0 = SRMYearFrac(extendedTimeLine.front(), dateFrom);
    for (unsigned int i = 0; i < g.size(); i++){
        if (idxN > idxM) {
            g[i] = model->factors[i].GfactorArrayIR[idxN] - 
                model->factors[i].GfactorArrayIR[idxM];
            double beta = model->factors[i].beta;
            double denom = r_not * exp(-beta * tFrom0);
            g[i] /= denom;
        } else {
            g[i] = 0.0;
        }
    }

    if (dateFrom != extendedTimeLine[idxM]) {
        double del_t = SRMYearFrac(dateFrom, extendedTimeLine[idxM]);
        for (unsigned int i = 0; i < g.size(); i++){
            double beta = model->factors[i].beta;
            g[i] += SRMUtil::GFAC(0, del_t, beta);
        }
    }
        
    if (dateTo != extendedTimeLine[idxN]) {
        double tFrom0 = SRMYearFrac(dateFrom, extendedTimeLine[idxN]);
        double del_t  = SRMYearFrac(extendedTimeLine[idxN], dateTo);
        double r  = rBar(extendedTimeLine[idxN]);
        for (unsigned int i = 0; i < g.size(); i++){
            double beta = model->factors[i].beta;
            g[i] += r * SRMUtil::GFAC(tFrom0, del_t, beta) / r_not ;
        }
    }
}

/** computes GFactor between the 2 dates supplied for specified model
    factor index */
double SRMRatesHJMUtil::gFactor(const DateTime& dateFrom,
                                const DateTime& dateTo,
                                int             factorIdx) const
{
    vector<double> gFac(factorIdx+1);
    gFactor(dateFrom, dateTo, gFac);
    return gFac[factorIdx];
}


/*  Calculate the 'stepped' tau given a date. From irdiffuse::STau */
const DateTime& SRMRatesHJMUtil::sTau(const DateTime& date) const
{
    int idx = date.findLower(extendedTimeLine);
    if (idx < 0){
        throw ModelException("SRMRatesHJMUtil::STau", "Date is before value date");
    }
    return (tau[idx]);
}

/** Utility routine to calculate the cutoff rate and populate
    the cutoff rate arrays */
void SRMRatesHJMUtil::calcEffRateLimit(
    double          NbSigmasMax, // (I) Number or sigmas to cut at
    double          NbSigmasMin, // (I) Number or sigmas to cut at
    vector<double>& MaxRate,     /* (O) Cutoff forward rates  */
    vector<double>& MinRate) const  /* (O) Cutoff forward rates      */
{
    if (swaptionExpiries.empty()){
        // no swaption vols so use old style. Why we add 10 here is undocumented
        calcEffRateLimitIR(NbSigmasMax+10.0, NbSigmasMin+10.0, MaxRate,MinRate);
    } else {
        // new style
        calcEffRateLimitIR_BM(NbSigmasMax, NbSigmasMin, MaxRate, MinRate);
    }
}

/***** From irdiffuse.c::CalcEffRateLimitIR  *******************************/
/*                                                                          */
/*       Utility routine to calculate the cutoff rate and populate          */
/*             the cutoff rate arrays                                       */
/*                                                                          */
/*       NOTE: arrays                                                       */
/*          Spot Volatilities (I)                                           */
/*          Forward Rates (I)                                               */
/*          Max Rates (O)                                                   */
/*       share the same set of Dates[0..NbDates-1]                          */
/*                                                                           */
void SRMRatesHJMUtil::calcEffRateLimitIR(
    double    NbSigmasMax,   // (I) Number or sigmas to cut at
    double    NbSigmasMin,   // (I) Number or sigmas to cut at
    vector<double>& MaxRate, /* (O) Cutoff forward rates  */
    vector<double>& MinRate) const  /* (O) Cutoff forward rates      */
{
    MaxRate.resize((*dates).size());
    MinRate.resize((*dates).size());

    MaxRate[0] = FwdRate[0];
    MinRate[0] = Maths::min(FwdRate[0], 0.0);

    double VolWeight = model->volWeight();
    double var       = 0.0;

    for (unsigned int i = 0; i < MinRate.size()-1; i++) {
        double SpotVolL = /* SpotVol[i] * */ VolWeight; // SpotVol[i] = 1.0
        var += SRMYearFrac((*dates)[i], (*dates)[i+1]) * SpotVolL * SpotVolL;
        /* lognormal MAX-Cutoff */
        double dummy = NbSigmasMax * sqrt(var);
        dummy = Maths::min(dummy, 100.0); /* to avoid blow-up in exp() */
        MaxRate[i+1] = FwdRate[i+1] * exp(dummy);

        /* normal MIN-Cutoff */
        MinRate[i+1] = FwdRate[i+1] * (1 - NbSigmasMin * sqrt(var));
        MinRate[i+1] = Maths::min(MinRate[i+1], 0.0);/* to avoid cutting at 
                                                        positive rate */
    }
}/* CalcEffRateLimitIR */


/** populates MaxEffRateIR and MinEffRateIR arrays. Requires tau. From
    irdiffuse::CalcEffRateLimitIR_BM */
void SRMRatesHJMUtil::calcEffRateLimitIR_BM(double          NbSigmasMax, 
                                            double          NbSigmasMin,
                                            vector<double>& MaxEffRateIR,
                                            vector<double>& MinEffRateIR) const
{
    static const string method("SRMRatesHJMUtil::calcEffRateLimitIR_BM");
    MaxEffRateIR.resize((*dates).size());
    MinEffRateIR.resize((*dates).size());
    /* At t= 0, the vol is deterministic  */
    MaxEffRateIR[0] = FwdRate[0];
    MinEffRateIR[0] = Maths::min(FwdRate[0], 0.0);
        
    /* Initialization */
    double sum_zeta = 0.0;
    double var = 0.0;
        
    /* single q approximation */
    double q = 0.5 * fabs(qLeft) + 0.5 * fabs(qRight);
    const DateTime& LastExpiry = getLastExpiry();
    vector<double> gfac(model->numFactors()); /* reserve some space for
                                                 gFactor calc */
    for (int i = 0; i < (*dates).size()-1; i++) {
        if (CString::equalsIgnoreCase(cutoffChoice, "constant")) {
            MaxEffRateIR[i+1] = constCutoffValue;
            MinEffRateIR[i+1] = 0.0; 
        } else {
            DateTime UsedTau;
            if ((*dates)[i] <= LastExpiry) {
                UsedTau = tau[i];
            } else {
                /* tau interpolation as implemented in calcSigma_r */
                const DateTime& firstT0 = sTau((*dates)[i]);
                const DateTime& lastT0 =  firstT0.max(LastDate);
                double denominator = SRMYearFrac(LastExpiry, LastDate);
                if (Maths::isZero(denominator)) {
                    throw ModelException(method, "Zero year fraction between "
                        + LastExpiry.toString()
                        + " and "
                        + LastDate.toString()); // not tested, more a theoretical case?
                }
                double slope = SRMYearFrac(firstT0, lastT0)
                    / denominator;
                int stepSize = (int) (slope * (*dates)[i].daysDiff(LastExpiry));
                UsedTau = firstT0.rollDate(stepSize);
            }
            
            double delt_1 = SRMYearFrac((*dates)[i], UsedTau);
            double delt_2 = SRMYearFrac((*dates)[i+1], UsedTau);
            // compute the variance
            double variance = model->variance(delt_1, delt_2);
            // then compute gFactor
            gFactor((*dates)[i], UsedTau, gfac);
            // and then the covariance between F and mu */
            // to do: see if we should move gFactor calc into covariance method
            double zeta = model->covariance(delt_1, delt_2,
                                            q * FwdRate[i] *
                                            Maths::square(SpotVol[i]), gfac);
            /* cumulative variance of F */
            var += Maths::square(q * SpotVol[i]) * variance;
        
            /* cumulative covariance between F and mu */
            sum_zeta += zeta;
        
            /* we build a confidence interval where exp(covariance) is
               the mean of the F under the risk-neutral measure (see
               report ). So that, MaxRate is applied to rbar * F under
               the risk-neutral measure */
            MaxEffRateIR[i+1] = FwdRate[i+1] * exp(-0.5 * var + NbSigmasMax * 
                                                   sqrt(var)) * exp(sum_zeta);
            MinEffRateIR[i+1] = FwdRate[i+1] * exp( - 0.5 * var - NbSigmasMin * 
                                                    sqrt(var)) * exp(sum_zeta);
            /* to avoid cutting at positive rate */
            MinEffRateIR[i+1] = Maths::min(MinEffRateIR[i+1], 0.0);
        }
    }   /* for i */
}

/* Calibration routine : populates SpotVol and tau array */
void SRMRatesHJMUtil::spotVolBM(bool skipFlag) // From swapvol::IR_SpotVol_BM
{
    try{
        swaptionSpotVol.resize(swaptionExpiries.size());
        tau = SRMSwaption::calibVol2Q(*this, skipFlag, swaptionSpotVol);

        vector<double> swaptionSpotVolNew(swaptionSpotVol.begin(), swaptionSpotVol.end());
        swaptionSpotVolAtSimDates = SRMUtil::extendVol(swaptionExpiries,
                                              swaptionSpotVolNew, 
                                              (*dates));

        /* Interpolate spot vol throughout the timeline */
        extendSpotVol(swaptionSpotVol);
    } catch (exception& e){
        throw ModelException(e, "SRMRatesHJMUtil::spotVolBM");
    }
}

/*
 *  Compute the xi factor as explained in the doc (section 3.1.2)
 *  Xi is the deterministic part of the measure change from the 
 *  tau-forward measure
 *  to the annuity measure. 
 *  From swapvol::XiFactor
 */
void SRMRatesHJMUtil::xiFactor(
    double&               xi,    /* (O) XiFactor */
    double&               der,   /* (O) derivative of xi with respect to tau */
    double&               maxXi, /* (O) maxXi is the vol of f(t,tau) */
    const vector<double>& weights,
    const DateTimeArray&  CpnPayDates,
    int                   NumCpns,
    const DateTime&       Tau,   /*(I) Tau value in the Newton-Raphson loop */
    const DateTime&       baseDate) const       /* (I) ZeroBaseDate */
{
    static const string method("SRMRatesUtil:xiFactor");
    /* Check that baseDate >= ir->FwdDate[0] */ 
    if (baseDate < getExtendedTimeLine().front()) {
        throw ModelException(method, " XiFactor : swap vol base date "
                             "must be >= today");
    }
    /* calculate rbar ratio */
    double denominator = rBar(baseDate);
    if (Maths::isZero(denominator)){
        throw ModelException(method, "Zero forward rate "
                             " on "+baseDate.toString());
    }

    double rbar_ratio = rBar(Tau) / denominator;
    /* Check that rbar_ratio is >= SRMConstants::SRM_TINY */
    if (fabs(rbar_ratio) < SRMConstants::SRM_TINY) {
        throw ModelException(method, "Ratio of forward rates < TINY between "+
                             Tau.toString()+" and "+baseDate.toString());
    }
    // reserve some space
    vector<double> xiByFactor(numFactors());
    vector<double> gfac(numFactors());
    /* calculate xi1 and its derivative with respect to tau for the
       Newton-Raphson algorithm */
    for (int j = 1; j <= NumCpns + 1;j++) {
        bool useTauDate = j == NumCpns + 1;
        gFactor(baseDate, useTauDate? Tau: CpnPayDates[j], gfac);
        for (unsigned int i = 0; i < xiByFactor.size(); i++){
            xiByFactor[i] -= (useTauDate? -1.0: weights[j]) * 
                gfac[i] * model->factors[i].alpha;
        }
    } /* for j*/
    vector<double> kfac(numFactors()); // reserve some space
    kFactor(baseDate, Tau, kfac);
    model->xiFactor(xi, der, maxXi, kfac, xiByFactor, rbar_ratio);
}

/** Calls bFactor for the swap used for correlation purposes */
vector<double> SRMRatesHJMUtil::bFactor() const
{
    return bFactor(corrSwapStart,
                   corrSwapMat,
                   corrSwapDCC,
                   MaturityPeriodSP(new MaturityPeriod(corrSwapFreq)));
}

/** Calls bFactor for the swap given by the specified index */
vector<double> SRMRatesHJMUtil::bFactor(int swapIndex) const
{
    return bFactor(swapStartDates[swapIndex], 
                   swapMatDates[swapIndex],
                   swapDCC, 
                   swapFrequency);
}

/*****  From swapvol.c: BFactor    *****************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*       this function is Flat Forward ready
*/
vector<double> SRMRatesHJMUtil::bFactor(const DateTime&      swapStart,
                                        const DateTime&      swapMat,
                                        DayCountConventionSP swapDayCC,
                                        MaturityPeriodSP     swapFreq) const
{
    /* MAR: Note that in the SRM3 code the zero curve passed in thinks that
       today is in fact value date (ie spot date) so all pv's etc are wrt
       spot date */
    static const string method("SRMRatesHJMUtil::bFactor");
    try{
        int freqCount;
        string freqInterval;
        swapFreq->decompose(freqCount, freqInterval);
        // This set of dates includes the swap start
        DateTimeArraySP dates(SwapTool::dateArray(swapStart, swapMat,
                                                  freqCount, freqInterval,
                                                  false/* stub at front */));
        const DateTime& spotDate = diffYC->getSpotDate();
        /* Calculate swap start in years from yield curve's spot date */
        double S = SRMYearFrac(spotDate, swapStart);
        double PrevT = S; // Same for previous coupon
        double Annuity  = 0.0; // Annuity price
        vector<double> B(numFactors()); // return value
        vector<double> A(B.size()); // reserve some space
        double ZerotoCpn = 0;     /* Zero to current cpn */

        /* Zero to swap start */
        auto_ptr<IYieldCurve::IKey> key(diffYC->logOfDiscFactorKey());
        double ZerotoS = exp(key->calc(spotDate, swapStart));

        double PrevZero = ZerotoS;

        for (int j = 1; j < dates->size(); j++) {
            const DateTime& CpnPmtDate = (*dates)[j];
            // compute the current coupon day count fraction
            double DCCFrac = swapDayCC->years((*dates)[j-1],CpnPmtDate);
            /* compute zero to current cpn */
            ZerotoCpn = exp(key->calc(spotDate, CpnPmtDate));
            /* compute the time to current coupon in years    */
            double TtoCpn   = SRMYearFrac(spotDate, CpnPmtDate);
            Annuity += DCCFrac * ZerotoCpn;
            for (unsigned int k = 0; k < B.size(); k++) {
                const SRMRatesFactor& factor = model->factors[k];
                /* cf Christian's memo for Vladimir's approximation */
                A[k] += exp (-factor.beta * (PrevT-S)) * 
                    log (PrevZero/ZerotoCpn) * factor.expDecay(TtoCpn-PrevT);
                B[k] += A[k] * DCCFrac * ZerotoCpn;
            }
            PrevT = TtoCpn;
            PrevZero = ZerotoCpn;
        }  /* for j */
        
        double FwdYield = (ZerotoS - ZerotoCpn) / Annuity; // Forward yield
        for (unsigned int k = 0; k < B.size(); k++) {
            B[k] *= FwdYield;
            B[k] += A[k] * ZerotoCpn;
            B[k] /= (ZerotoS - ZerotoCpn) ;
        }
        return B;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}  /* bFactor */

// returns the instantaneous vol by factor in the form (factor, time point)
void SRMRatesHJMUtil::instFactorVol(
			vector< vector<double> >& vol,				  
	        const vector<double>& DeltaTime,  // (I)  passed for convenience 
	        const DateTimeArray& TPDate,      // (I)
            const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
	        int expiryIndex,                  // (I)  index in TPDate
			int fwdMatIndex) const            // (I)  in case, the underlying has mat > expiry
{
	static const string method("SRMRatesHJMUtil::instFactorVol");

	try
	{
	    // checks 
		if ( ((int)DeltaTime.size() < TPDate.size() - 1) )
		{
			throw ModelException(method, "DeltaTime vector is too short");
		}

		int NbFac = this->numFactors();
		int NbTP = TPDate.size();

		if ((int)vol.size() < NbFac)
		{
			throw ModelException(method, "Factor dimension of vol is too low");
		}

		int i;
        for(i = 0; i < NbFac; i++)
		{
			if ( ((int)vol[i].size() < expiryIndex) )
			{
				throw ModelException(method, "Expiry dimension of vol is too low");
			}
		}

		// extend IR vols 
		// NB : one less spot vol than dates
		vector<double> irExtSpotVol(NbTP-1);
		const DateTimeArray& irDates = this->getExtendedTimeLine();
		int irNbPoint = irDates.size();
		const vector<double>& irSpotVol = this->getSpotVols();
        int k;

		for (i = 0, k = 1 ; i < NbTP - 1; i++) 
		{
			if (k < (irNbPoint - 1)) 
			{
				if (TPDate[i] >= irDates[k])
				{
					k++;
				}
			}
			irExtSpotVol[i] = irSpotVol[k-1];
		}

		vector< SRMRatesFactor > factor = (this->model)->factors;

		// reset every time
		vector<double> IrA(NbFac);

		for (i = fwdMatIndex - 1; i >= expiryIndex; i--) 
		{
			this->aFactor(IrFwdRate[i], DeltaTime[i], IrA);
		}

		for(i = expiryIndex - 1; i >= 0; i--)
		{
			this->aFactor(IrFwdRate[i], DeltaTime[i], IrA);	
			for (int l = 0; l < NbFac; l++) 
			{
				vol[l][i] = factor[l].alpha * IrA[l] * irExtSpotVol[i];
			}/* for l */
		}

		return;

	}
	catch (exception& e)
	{
		throw ModelException(e, method, "For currency "+discYC->getCcy());
	}
}

// Calculates a 'factor variance' as used by equity/FX vol calibration 
void SRMRatesHJMUtil::irAssetVariance(
    vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
    vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
	const vector<double>& rhoEqIr,     // (I)
	const vector<double>& DeltaTime,   // (I)  passed for convenience 
	const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
    const vector<double>& IrFwdRate,   // (I)  pre-computed for performance
	const vector<int>& periodIndex,    // (I)  indices in TPDate
 	int fwdMatIndex) const             // (I)  indices in TPDate
{
	static const string method("SRMRatesHJMUtil::irAssetVariance");

	try
	{
	    // checks 
		if (irVariance.size() != irAssetCovar.size())
		{
			throw ModelException(method, "irVariance and irAssetCovar have different lengths");
		}

		if (irVariance.size() < periodIndex.size())
		{
			throw ModelException(method, "irVariance vector is too short");
		}

		if ((int)DeltaTime.size() < TPDate.size() - 1)
		{
			throw ModelException(method, "DeltaTime vector is too short");
		}

		int IrNbFac = this->numFactors(); 
		int NbNoEqInt = IrNbFac + (IrNbFac - 1) * IrNbFac / 2;
		int NbWithEqInt = IrNbFac;

		if ((int)rhoEqIr.size() < IrNbFac)
		{
			throw ModelException(method, "Too few correlations between EQ and IR factors");
		}

		// extend IR vols 
		int NbTP = TPDate.size();
		vector<double> irExtSpotVol(NbTP-1);
		const DateTimeArray& irDates = this->getExtendedTimeLine();
		int irNbPoint = irDates.size();
		const vector<double>& irSpotVol = this->getSpotVols();
		int i, k, l;

		for (i = 0, k = 1 ; i <= NbTP - 2; i++) 
		{
			if (k < (irNbPoint - 1)) 
			{
				if (TPDate[i] >= irDates[k]) 
				{
					k++;
				}
			}
			irExtSpotVol[i] = irSpotVol[k-1];
		}

		vector<double> IrA(IrNbFac);

		int NbIndices = periodIndex.size();
		int beginIndex = periodIndex[NbIndices-1];
		int endIndex = beginIndex;

		if (fwdMatIndex < endIndex)
		{
			throw ModelException(method, "forward maturity < expiry");
		}

		for (i = fwdMatIndex - 1; i >= endIndex; i--) 
		{
			this->aFactor(IrFwdRate[i], DeltaTime[i], IrA);
		}

		// start calculation of variances
		for (k = NbIndices - 1; k >= 0; k--)  
		{   
			endIndex = beginIndex;
			beginIndex = (k>0) ? periodIndex[k-1] : 0; 

			// reset
			vector<double> IrVolOnly(NbNoEqInt);  // NbNoEqInt    = 1 (IR=1F), 3 (IR=2F), 6 (IR=3F) (s1,s2,s3,s12,s13,s23)
			vector<double> IrAsset(NbWithEqInt);  // NbWithEqInt  = 1 (IR=1F), 2 (IR=2F), 3 (IR=3F) (s1e,s2e,s3e)

			for (i = endIndex - 1; i >= beginIndex; i--) 
			{
				double deltaTime = DeltaTime[i]; // for ease
				double irVol = irExtSpotVol[i];

				// NB. SRM3	has	some adjustments due to	dividends, but also
				// requires	that there be no dollar	divs. This actually
				// means the adjustment	is 0. So, here we ignore the
				// adjustment code entirely.
				this->aFactor(IrFwdRate[i], deltaTime, IrA);

				// variance of IR factors
				vector<double>::iterator irVar = IrVolOnly.begin();
				irVar = 
					this->factorVariance(irVar, IrA,
					irVol, deltaTime);
				// covariance between IR factors
				irVar = 
					this->factorCovariance(irVar, IrA,
					irVol, deltaTime);
				// covariance between EQ/FX and IR
				vector<double>::iterator irEqCoVar = IrAsset.begin();
				irEqCoVar = 
					this->factorFXCovariance(true, // "+" and not "-" contribution
					irEqCoVar, IrA, rhoEqIr,
					irVol, deltaTime);
			}   /* for i*/

			irVariance[k] = 0.;
			irAssetCovar[k] = 0.;

			for(l = 0; l < NbNoEqInt; l++)
			{
				irVariance[k] += IrVolOnly[l];
			}

			for(l = 0; l < NbWithEqInt; l++)
			{
				irAssetCovar[k] += IrAsset[l];
			}
		}

		return;

	}
	catch (exception& e)
	{
		throw ModelException(e, method, "For currency "+discYC->getCcy());
	}
}

/*************** From util_s.c:Triangulation ******************************
* Produces the "usual" Lower triangular matrix for orthogonalising corrolated
* factors.  
* NOTE: If Nbfac < 3 then unused matrix elements are set to zero.
*
****************************************************************************/
DoubleMatrix SRMRatesHJMUtil::triangulation() const
{
    static const string method("SRMRatesHJMUtil::triangulation");
    // for ease
    int Nbfac = numFactors();   /* Number of factors */
    const vector<double>& Rho = model->rho;  /*(I) correlation coefficients */
    /* initialise matrix (return value) */
    DoubleMatrix TriangMtx(Nbfac, Nbfac);
    TriangMtx[0][0] = 1.;
    if (Nbfac > 1) {
        TriangMtx[1][0] = Rho[0];
        if (1.0- Rho[0] * Rho[0] < SRMConstants::SRM_TINY) {
            throw ModelException(method,
                                 "bad choice of correlation parameter Rho[0]");
        }
        TriangMtx[1][1] = sqrt(1.0 - Rho[0] * Rho[0]);
    }
    if (Nbfac > 2) {
        TriangMtx[2][0] = Rho[1];
        TriangMtx[2][1] = (Rho[2] - Rho[0] * Rho[1]) / (TriangMtx[1][1]);
        if ((1.0 - Rho[0] * Rho[0] - Rho[1] * Rho[1] - Rho[2] * Rho[2]
             + 2.0 * Rho[0] * Rho[1] * Rho[2]) < SRMConstants::SRM_TINY) {
            throw ModelException(method, "Bad choice of correlation "
                                 "parameters:\ncorrelation matrix is not "
                                 "positive definite!!\n");
        }
        TriangMtx[2][2] = sqrt(1.0 - Rho[0] * Rho[0] - 
                               Rho[1] * Rho[1] - Rho[2] * Rho[2]
                               + 2.0 * Rho[0] * Rho[1] * Rho[2]) 
            / TriangMtx[1][1];
    }
    return TriangMtx;
}/* Triangulation */


/*****  From util_s::Get3A **************************************************/
/*
 *       Calculate the decomposition of the Correlation in a 2,3-factor case
 *       given the swap start date, a tenor and a frequency
 */
DoubleMatrix SRMRatesHJMUtil::get3A() const
{
    static const string method("SRMRatesHJMUtil::get3A");
    if (!corrDecomposition.empty()){
        return corrDecomposition; // this is typically asked for many times
    }
    try{
        // ooh, this is classy
        double    t2Mat[] = {1,3,6,12,24,36,48,60,84,120,144,180,240,360};
        int       Nbt2Mat = sizeof(t2Mat)/sizeof(double);
        
        int NbFactor = numFactors(); //  for ease
        if (NbFactor == 1){
            DoubleMatrix OMtx(1, 1); // the results
            OMtx[0][0] = 1.0;
            return OMtx;
        }
        const vector<double>& rho = model->rho; // for ease
        vector<double>  beta(NbFactor);
        vector<double>  alpha(NbFactor);
        getBeta(beta); // for ease
        getAlpha(alpha); // for ease
        bool areBetasAllDiff;
        DoubleMatrix RHO(NbFactor, NbFactor);
        if (NbFactor == 2){
            RHO[0][0] = RHO[1][1] = 1.0;
            RHO[0][1] = RHO[1][0] = rho[0];
            areBetasAllDiff = (!Maths::equals(beta[0], beta[1]));
        } else {
            RHO[0][0] = RHO[1][1] = RHO[2][2] = 1.0;
            RHO[0][1] = RHO[1][0] = rho[0];
            RHO[0][2] = RHO[2][0] = rho[1];
            RHO[1][2] = RHO[2][1] = rho[2];
            areBetasAllDiff = ((!Maths::equals(beta[0], beta[1])) &&
                               (!Maths::equals(beta[0], beta[2])) &&
                               (!Maths::equals(beta[1], beta[2])));
        }
        if (!areBetasAllDiff) {
            /* Some betas are the same, therefore degeneration. */
            /* Use the lower triangular decomp instead          */
            throw ModelException(method, "Some betas levels are the same. "
                                 "Cannot complete decomposition.");
            // the next two lines are commented out in SRM3
            /* if (SqrtSymMatrix(NbFactor, RHO, OMtx) == FAILURE) goto RETURN;*/
            /* return (SUCCESS);   ???                                        */
        }

        /* create the rate covariance matrix */
        /* rCovMtx[1..Nbt2Mat][1..Nbt2Mat]   */
        DoubleMatrix rCovMtx(Nbt2Mat, Nbt2Mat);
        int i;
        for (i = 0; i < Nbt2Mat; i++) {
            for (int l = 0; l < Nbt2Mat; l++) {
                rCovMtx[i][l] = 0.0;
                for (int j = 0; j<NbFactor; j++) {
                    for (int k = 0; k < NbFactor; k++) {
                        rCovMtx[i][l] += 
                            (alpha[k] * exp(-beta[k] * t2Mat[l] / 12.0)) *
                            (alpha[j] * exp(-beta[j] * t2Mat[i] / 12.0)) * 
                            RHO[j][k];

                    } /* k */
                } /* j */
            } /* l */
        } /* i */

        /* eigenise the rate covariance matrix              */
        EigenVectorAnalysisSP eigenVecData(rCovMtx.computeEigenVectors());
        /* create the inverse of A */
        DoubleMatrix A(NbFactor, NbFactor);
        for (i = 0; i < NbFactor; i++) {
            // sqrt the eigenValues
            double& eigenValue = (*eigenVecData->eigenValues)[i];
            if (!Maths::isPositive(eigenValue)){
                // we do 1/sqrt so a small +ve is bad news too
                throw ModelException(method, "Zero eigenvalue!");
            }
            eigenValue = sqrt(eigenValue);
            for (int l = 0; l < NbFactor; l++) {
                A[i][l] = 0.0;
                for (int k = 0; k < Nbt2Mat; k++) {
                    A[i][l] += 1.0/(*eigenVecData->eigenValues)[i] *
                        (*eigenVecData->eigenVectors)[i][k] *
                        (alpha[l] * exp(-beta[l] * t2Mat[k] / 12.0));
                } /* k */
            } /* l */
        } /* i */

        /* find A by inverting inverse of A */
        corrDecomposition = A.computeInverse(); // cache
        return corrDecomposition;
    } catch (exception& e){
        throw ModelException(e,method, "Failed for currency "+discYC->getCcy());
    }
} /* Get3A */

// Calculates 'A factor' at next date using 'A factor' from previous date
// and supplied parameters 
// Note that fwdRateDt = fwdRate * dt
void SRMRatesHJMUtil::aFactor(double fwdRateDt,
                        double deltaTime,
                        vector<double>& a) const // (M)
{
    for (unsigned int l = 0; l < model->factors.size(); l++) {
        a[l] *= exp(-model->factors[l].beta * deltaTime);
        a[l] += fwdRateDt * model->factors[l].expDecay(deltaTime);
    } 
}

// Calculates a 'factor variance' as used by FX vol calibration 
vector<double>::iterator SRMRatesHJMUtil::factorVariance(
    vector<double>::iterator integral,  // (M)
    const vector<double>&    aFactor,
    double                   spotVol,
    double                   deltaTime) const
{
    // variance of domestic factors
    for (unsigned int l = 0;l < model->factors.size();l++) {   
        *integral += Maths::square(model->factors[l].alpha *
                                   spotVol * aFactor[l]) * 
            deltaTime;
        ++integral;
    }
    return integral;
}

// Calculates a 'factor covariance' as used by FX vol calibration 
vector<double>::iterator SRMRatesHJMUtil::factorCovariance(
    vector<double>::iterator integral,  // (M)
    const vector<double>&    aFactor,
    double                   spotVol,
    double                   deltaTime) const
{
    // covariance of domestic factors 
    // Methodology could be promoted to factors class
    int numFactors = model->factors.size();
    double spotVolSq = Maths::square(spotVol);
    if (numFactors > 1) {
        *integral += 2.0 * model->factors[0].alpha * model->factors[1].alpha * 
            aFactor[0] * aFactor[1]
            * spotVolSq *
            model->rho[0]*deltaTime;
        ++integral;
    }
    if (numFactors > 2) {   
        *integral += 2. * model->factors[0].alpha * model->factors[2].alpha *
            aFactor[0] * aFactor[2]
            * spotVolSq * 
            model->rho[1] * deltaTime;
        ++integral;
        *integral += 2.0 * model->factors[1].alpha * model->factors[2].alpha *
            aFactor[1] * aFactor[2]
            * spotVolSq * 
            model->rho[2] * deltaTime;
        ++integral;
    }
    return integral;
}


/** Calculates a 'factor covariance' as used by FX vol calibration */
vector<double>::iterator SRMRatesHJMUtil::factorFXCovariance(
    bool                     isDomestic,
    vector<double>::iterator integral,  // (M)
    const vector<double>&    aFactor,
    const vector<double>&    rhoFxIRFac, /* (I) correl IR/FX */
    double                   spotVol,
    double                   deltaTime) const
{
    /*covariance between Fx and IR */
    int numFac = model->factors.size();
    for (int l = 0; l < numFac; l++) {
        *integral += (isDomestic? 2.0: -2.0) * spotVol * 
            model->factors[l].alpha * aFactor[l] * rhoFxIRFac[l] * deltaTime;
        ++integral;
    }/* for l*/
    return integral;
}

/** this function returns product of spot vols and forward rates
    "interest rate basis point vol". The array returned is of the same
    length as the supplied array. 
    From CMLib:SpreadCurveAlgorithms.cpp:BasisPointVol */
DoubleArray SRMRatesHJMUtil::basisPointVol(
    const DateTimeArray& dates) const /* excludes today */
{
    const DateTimeArray& volDates = swaptionExpiries; // make port easier

    // currentVol is in effect until volEndDate
    // after that volEndDate is updated with volDates[nextVolIndex]
    // we check that the dates array contains all dates from volDates
    ASSERT(DateTime::isSubset(dates, volDates));

    DateTime volEndDate = baseDate;
    int nextVolIndex = 0;

    int numDates = dates.size(); // for ease
    DoubleArray rBpVol(numDates);
    double currentVol=0;
    Actual365F act365F;
    for (int n = 0; n < numDates; n++) {
        const DateTime& startDate = n == 0? baseDate: dates[n-1];
        if (startDate == volEndDate) {
            bool behindLastVol = nextVolIndex == volDates.size();

            // for dates past the last vol date we take overnight forward vol
            volEndDate = behindLastVol? 
                (volDates.back().rollDate(1)): volDates[nextVolIndex];
            double spotVol = behindLastVol? 
                swaptionSpotVol.back(): swaptionSpotVol[nextVolIndex++];
            currentVol = spotVol *
                diffYC->fwd(startDate, volEndDate, 
                            &act365F, CompoundBasis::ANNUAL); // why ANNUAL?

            if (behindLastVol){
                // if we are past the last vol 
                // use current vol for the rest of the timeline
                volEndDate = dates.back();
            }
        }
        ASSERT(startDate < volEndDate);
        
        rBpVol[n] = currentVol; // save the vol
    }
    return rBpVol;
}

DRLIB_END_NAMESPACE
