//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrCov.cpp
//
//   Description : Correlation/Covariance Swap
//
//   Date        : December 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/CorrSwapBasisAdj.hpp"
#include "edginc/CorrSwapSamplingAdj.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/MultiAsset.hpp"

DRLIB_BEGIN_NAMESPACE

/** CorrCov product */
class CorrCov: public GenericNFBase,
                           virtual public IMCIntoProduct{
private:
    // fields
    IDateBuilderSP          averageOutDates;
    DoubleArray             weights;
    double                  strike;
    bool                    isCorr;             // corr or cov
    bool                    isSpreadOfVarSwaps; // if isCorr = FALSE
    DoubleArray             spreadWeights;
    bool                    useMean;            // if isSpreadOfVarSwaps = FALSE
    bool                    useVolInWeights;    // if isCorr = FALSE and isSpreadOfVarSwaps = FALSE
    bool                    isCapped;
    double                  cap;
    bool                    isFloored;
    double                  floor;

    // transient fields
    // dates after ISDA adjustment
    // may have some dates omitted from original averageOutDates
    // may also contain duplicates
    DateTimeClusterConstSP  sampleDates;
    int nbDates;

public:
    static CClassConstSP const TYPE;
    friend class CorrCovMC;

    /** Get the asset and discount market data */
    virtual void GetMarket(const IModel*          model,
                           const CMarketDataSP    market) {
        static const string method("CorrCov::GetMarket");
        try {
            MarketObjectArraySP mo1 = market->GetAllDataWithType(CorrSwapBasisAdj::TYPE);
            MarketObjectArraySP mo2 = market->GetAllDataWithType(CorrSwapSamplingAdj::TYPE);
            if ( (mo1->size()>0) || (mo2->size()>0) ) {
                DateTimeSP maturity(new DateTime(averageOutDates->end()));
                MarketDataFetcherSP mdf(model->getMDF());
                // ISS - does this need to change with ISDA? Does it matter?
                mdf->setCorrSwapExpiry(maturity);
            }
            GenericNFBase::GetMarket(model,market);
            averageOutDates->getMarket(model,market.get());
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** returns the array of all sampling dates the instrument will ever need
        excluding (possibly) the ref level dates (handled in GenericNFBase)
        Note this may not be the final list as dates might be ISDA omitted/moved*/
    const DateTimeArray samplingDates() const {
        return *(averageOutDates->dates());
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // validation
    void validatePop2Object(){
        static const string routine("CorrCovMC::validatePop2Object");
        GenericNFBase::validatePop2Object();
        // check that we've got as many weights as there are assets
        // and if % weights that they sum to 100%
        AssetUtil::checkWeights(weights, assets->NbAssets());
        if (assets->NbAssets() != weights.size()){
            throw ModelException(routine,
                                 "Different number of assets ("+
                                 Format::toString(assets->NbAssets())+")"
                                 " to weights ("+
                                 Format::toString(weights.size())+")");
        }

        // validate dates are not empty and are in order
        DateTimeArraySP monDates = averageOutDates->dates();
        DateTime::ensureIncreasing(*monDates, "averageOutDates", true);

        // need at least two returns
        if( averageOutDates->size()<2 ) {
            throw ModelException(routine,
                                 "At least two dates required, "+
                                 Format::toString(averageOutDates->size())+
                                 " provided.");
        }

        // check flags
        if( isSpreadOfVarSwaps ) {
            throw ModelException(routine,
                                 "Spread Of Variance Swap is not supported.");
        }
        if( isCorr && isSpreadOfVarSwaps ) {
            throw ModelException(routine,
                                 "isCorr=TRUE and isSpreadOfVarSwaps=TRUE not supported.");
        }
        if( useVolInWeights && (isSpreadOfVarSwaps||!isCorr) ) {
            throw ModelException(routine,
                                 "useVolInWeights=TRUE and (isSpreadOfVarSwaps=TRUE or isCorr=FALSE) not supported.");
        }
        if( (isCapped || isFloored) && isSpreadOfVarSwaps ) {
            throw ModelException(routine,
                                 "isSpreadOfVarSwaps=TRUE and (isCapped=TRUE or isFloored=TRUE) not supported.");
        }

        // make sure that lurking values provided are not used if not required...
        if( !isCapped ) {
            cap = 1.e60;
        }
        if( !isFloored ) {
            floor = -1.e60;
        }

        // cap>=floor
        if( cap<floor ) {
            throw ModelException(routine,
                                 "Cap "+ Format::toString(cap) + " has to be below floor, " + Format::toString(floor) );
        }
    }

    /** Validate instrument having aquired market data */
    virtual void Validate() {
        static const string routine("CorrCovMC::Validate");
        // this call will handle the ISDA date generation and
        // validation of past values
        GenericNFBase::Validate();

        // some validation for ISDA itself
        // given that covariance swaps and spread of variance swaps involve
        // dt we don't know what this means in the presence of ISDA when dates
        // may roll independently
        if (!isdaDateAdjust->isUnadjusted() &&
                    (isSpreadOfVarSwaps || !isCorr)){
            throw ModelException(routine,
                        "Cannot have ISDA adjustment with Covariance swaps or"
                        " Spread of Variance Swaps. Please select None as the "
                        "ISDA adjustment rule");
        }
        sampleDates = DateTimeClusterConstSP(&(obsMap->getSampleDates()));
        // pointless checks hopefully
        nbDates = 0;
        for (int iAsset=0; iAsset<assets->NbAssets(); iAsset++ ) {
            // need at least two returns
            if((*sampleDates)[iAsset].size()<2 ) {
                throw ModelException(routine,
                                    "At least two dates required, "+
                                    Format::toString((*sampleDates)[iAsset].size())+
                                    " remains after ISDA adjustment for asset "+
                                    Format::toString(iAsset + 1) + ".");
            }
            if (iAsset == 0) {
                nbDates = (*sampleDates)[iAsset].size();
            } else {
                if (nbDates != (*sampleDates)[iAsset].size()) {
                    throw ModelException(routine,
                                        "There must be the same number of averaging"
                                        " dates for each asset after ISDA adjsutment. "
                                        "Asset 1 has " +
                                        Format::toString(nbDates)+
                                        " while asset "+
                                        Format::toString(iAsset+1) +
                                        "has " +
                                        Format::toString((*sampleDates)[iAsset].size()) +
                                        ".");
                }
            }
        }
    }

private:
    CorrCov(): GenericNFBase(TYPE), isCapped(false), cap(1.e60), isFloored(false),
               floor(-1.e60), nbDates(0) {}       // for reflection
    CorrCov(const CorrCov& rhs);            // not implemented
    CorrCov& operator=(const CorrCov& rhs); // not implemented

    static IObject* defaultCorrCov(){
        return new CorrCov();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CorrCov, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultCorrCov);
        FIELD(averageOutDates,   "averageOutDates");
        FIELD(weights,           "Weights");
        FIELD(strike,            "Strike");
        FIELD(isCorr,            "Corr (true) or Cov (false)");
        FIELD(isSpreadOfVarSwaps,"Cov as difference of bPath var and weighted sum of vars (true)");
        FIELD(spreadWeights,     "Weights of component variances (case isSpreadOfVarSwaps)");
        FIELD(useMean,           "Use sampled mean for covariance (true) or not (false)");
        FIELD(useVolInWeights,   "avg cov / (avg vol * avg vol) (true) or avg corr (false)");
        FIELD(isCapped,          "Flag: if true, expect cap");
        FIELD_MAKE_OPTIONAL(isCapped);
        FIELD(cap,               "cap on floating leg");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(isFloored,         "Flag: if true, expect floor");
        FIELD_MAKE_OPTIONAL(isFloored);
        FIELD(floor,             "floor on floating leg");
        FIELD_MAKE_OPTIONAL(floor);
        FIELD(sampleDates, "Average dates after ISDA adjustment");
        FIELD_MAKE_TRANSIENT(sampleDates);
        FIELD(nbDates, "Number of dates after ISDA adjustment");
        FIELD_MAKE_TRANSIENT(nbDates);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super CorrCov */
class CorrCovMC : public IMCProduct,
                        virtual public IMCProductLN {
private:
    const CorrCov*           inst;       // reference to original instrument
    int                      nbAssets;   // for convenience
    CDoubleArraySP           mean;
    CDoubleArraySP           meanSum;
    CDoubleArraySP           variance;
    CDoubleArraySP           varianceSum;
    CDoubleMatrixSP          covariance;
    CDoubleMatrixSP          covarianceSoFar;
    CDoubleMatrixSP          covarianceSum;
    int                      mcRuns;
    CDoubleArraySP           weightFactor; // working area - here to save alloc
    CDoubleArraySP           weightFactorNorm; // working area - here to save alloc
    int                      nbDates;
    double                   dt;
    double                   norm;
    double                   cleanPrice;
    double                   twoOverTime;
    CDoubleArraySP           logFwd;
    double                   logbFwd;
    DoubleArray              singleVarPast;
    DoubleArray              singleVarFuture;
    DoubleArray              singleVarTotal;
    DoubleArray              singleVarTotalSum;
    double                   basketVarPast;
    double                   basketVarFuture;
    double                   basketVarTotal;
    double                   basketVarTotalSum;
    double                   weightPast;
    int                      nbDatesPast;
    DoubleArray              logAssetReturn; // [nbAssets]

public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to CorrCov) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    CorrCovMC(  const CorrCov*      inst,
                const SimSeriesSP&  simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        mcRuns(0),
        logAssetReturn(nbAssets) {

            static const string routine = "CorrCovMC::CorrCovMC";
            mean            = CDoubleArraySP(new CDoubleArray(nbAssets));
            meanSum         = CDoubleArraySP(new CDoubleArray(nbAssets));
            variance        = CDoubleArraySP(new CDoubleArray(nbAssets));
            varianceSum     = CDoubleArraySP(new CDoubleArray(nbAssets));
            covariance      = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));
            covarianceSoFar = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));
            covarianceSum   = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));
            weightFactor    = CDoubleArraySP(new CDoubleArray(nbAssets));
            weightFactorNorm= CDoubleArraySP(new CDoubleArray(nbAssets));
            int iAsset,jAsset;

            for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                CAssetConstSP asset = CAssetConstSP::attachToRef(&getMultiFactors()->getAsset(iAsset));

                if (StruckEquity::TYPE->isInstance(asset) && !(AssetUtil::isBasket(asset))) {
                    throw ModelException(routine, "Struck assets are not supported");
                }

                if( inst->isCorr && !inst->useVolInWeights ) {  // isCorr && !useVolInWeights
                    (*covariance)[iAsset][iAsset] = 1;
                }
                else{                                           // !isCorr || useVolInWeights
                    (*weightFactor)[iAsset]       = 1;
                }
            }

            // ISS don't know what to do here in presence of ISDA
            // for now we validate against this case so it's safe to use sampleDates
            nbDates         = inst->nbDates;
            dt              = (double) (inst->averageOutDates->end().getDate()
                                        - inst->averageOutDates->start().getDate()        )
                              / (double)(nbDates*365); // use metric->yearFrac ?
            if( inst->useVolInWeights ) {
                norm = 1;
            }
            else {
                norm = 0;
                for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                    for( jAsset=iAsset+1; jAsset<nbAssets; jAsset++ ) {
                        norm    += inst->weights[iAsset] * inst->weights[jAsset];
                    }
                }
                if( norm!=0 ) {
                    norm = 1/norm;
                }
                if( inst->isSpreadOfVarSwaps ) {
                    norm = norm/2;
                }
                else if( !inst->isCorr ) {
                    norm = norm/dt;
                }
            }
            cleanPrice      = 0;
            // use metric->yearFrac ?

            //compute twoOverTime, logFwd, logbFwd only when inst has not expired, which are used only in futureVar...
            // ISS don't know what to do here in presence of ISDA
            // for now we validate against this case so it's safe to use averageDates
            if (inst->isSpreadOfVarSwaps) {
                twoOverTime = 0.0;
                logFwd      = CDoubleArraySP(new CDoubleArray(nbAssets,0.0));
                logbFwd     = 0;
                if (inst->averageOutDates->end().isGreater( getToday() )) {
                    twoOverTime
                        = 730.0 / (double)( inst->averageOutDates->end().getDate()
                        //- Maths::max(startDate.getDate(),getToday().getDate()) );
                            - getToday().getDate() );
                    double bFwd = 0;
                    for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                        double fwd  = getMultiFactors()->
                            assetFwdValue(iAsset,inst->averageOutDates->end());
                        (*logFwd)[iAsset] = log(fwd);
                        bFwd  += inst->weights[iAsset] * fwd;
                            bFwd  += inst->weights[iAsset] * fwd;
                    }
                    logbFwd = log(bFwd);
                }
            }        

            singleVarPast       = DoubleArray(nbAssets);
            singleVarFuture     = DoubleArray(nbAssets);
            singleVarTotal      = DoubleArray(nbAssets);
            singleVarTotalSum   = DoubleArray(nbAssets);
            basketVarPast       = 0;
            basketVarFuture     = 0;
            basketVarTotal      = 0;
            basketVarTotalSum   = 0;
            weightPast          = 0;
            nbDatesPast         = 0;
        }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*   pathGen,
                IMCPrices&              prices) {
        static const string method("CorrCovMC::payoff");
        try {
            int    iAsset, jAsset, iDate;

            if( inst->isSpreadOfVarSwaps ) {

                if (pathGen->doingPast()) {
                    historicalVar( pathGen, basketVarPast, singleVarPast, nbDatesPast);
                    weightPast = (double)nbDatesPast/(double)nbDates;
                }
                else {
                    futureVar( pathGen, basketVarFuture, singleVarFuture, nbDates);
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                        singleVarTotal[iAsset]
                            = weightPast * singleVarPast[iAsset]
                            + (1-weightPast) * singleVarFuture[iAsset];
                        singleVarTotalSum[iAsset]
                            += singleVarTotal[iAsset];
                    }
                    basketVarTotal
                        = weightPast * basketVarPast
                        + (1-weightPast) * basketVarFuture;
                    basketVarTotalSum += basketVarTotal;
                    mcRuns++;
                }

                double covariance = basketVarTotal;
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    covariance -= inst->spreadWeights[iAsset] * singleVarTotal[iAsset];
                }

                cleanPrice += covariance*norm;
                prices.add( inst->notional * (covariance*norm-inst->strike) );

            } // end if isSpreadOfVarSwaps
            else
            {
                CDoubleMatrixSP thisCovariance;
                int    beginIdx = 0;
                int    endIdx   = pathGen->end(0);
                // switch pointers past/future
                if (pathGen->doingPast()) {
                    thisCovariance = covarianceSoFar;

                    //if endIdx == 0, we are before the 'averaging' period so there is no need for any
                    // sampling adjustments
                    if (endIdx != 0) {
                        // work out the last completely fixed sample set
                        nbDatesPast = inst->obsMap->getLastSampleIndex(0/*assetidx*/,
                                                                    pathGen->end(0)-1);
                        for (iAsset = 1; iAsset <nbAssets; iAsset++) {
                            nbDatesPast = Maths::min(nbDatesPast,
                                                    inst->obsMap->getLastSampleIndex(iAsset,
                                                                                    pathGen->end(iAsset)-1));
                        }
                        nbDatesPast++; // to be consistent with isSpreadOfVarSwaps
                        endIdx = nbDatesPast;
                    }
                }
                else {
                    thisCovariance = covariance;
                    // initialize covariance (=sum of product of log returns with/w'out mean) by past sum
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                        for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                            (*covariance)[iAsset][jAsset] = (*covarianceSoFar)[iAsset][jAsset];
                        }
                    }
                    beginIdx = nbDatesPast;
                    endIdx = nbDates;
                }

                // cov w'out mean
                // We do this carefully minimising the number of log() calls since perf
                // can be seriously affected
                for( iDate=beginIdx; iDate<endIdx; iDate++ ) {
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                        logAssetReturn[iAsset] =
                            log( pathGen->Path(iAsset, 0)[inst->obsMap->getModellingIndex(iAsset, iDate)]
                                / ( (iDate>0) ?
                                    pathGen->Path(iAsset, 0)[inst->obsMap->getModellingIndex(iAsset, iDate-1)] :
                                    pathGen->refLevel(iAsset, 0) ) );
                    }
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                        // for ease 'n speed
                        double* iCov        = (*thisCovariance)[iAsset];

                        double logReturn1 = logAssetReturn[iAsset];
                        for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                            double oneCov = logReturn1 * logAssetReturn[jAsset];
                            iCov[jAsset] += oneCov;
                        }
                    }
                }
                double average  = 0;
                double denom    = inst->useVolInWeights ? 0 : 1;

                if (!pathGen->doingPast()||!hasFuture()) {

                    // diagonal
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {

                        // mean = log( lastSpot / refLevel ), get it wrong for past path...
                        (*mean)[iAsset]
                            = log( pathGen->Path(iAsset, 0)[inst->obsMap->getModellingIndex(iAsset, endIdx-1)]
                                    / pathGen->refLevel(iAsset, 0) ) / nbDates;
                        // for ease 'n speed
                        double* iCov        = (*thisCovariance)[iAsset];
                        double* iCovSum     = (*covarianceSum)[iAsset];

                        // variance
                        if( inst->useMean ) {
                            (*variance)[iAsset]
                                = ( iCov[iAsset] - nbDates * (*mean)[iAsset] * (*mean)[iAsset] ) / (nbDates-1);
                        }
                        else {
                            (*variance)[iAsset] = iCov[iAsset] / nbDates;
                        }

                        // weightFactor in denominator, weightFactorNorm in numerator...
                        if( inst->useVolInWeights ) {   // isCorr && useVolInWeights
                            iCov[iAsset]                = (*variance)[iAsset];
                            (*weightFactorNorm)[iAsset] = sqrt((*variance)[iAsset])*inst->weights[iAsset];
                        }
                        else if( inst->isCorr ) {       // isCorr && !useVolInWeights
                            (*weightFactor)[iAsset]     = 1 / sqrt((*variance)[iAsset]);
                            iCov[iAsset]                = 1;
                        }
                        else{                           // !isCorr
                            iCov[iAsset]                = (*variance)[iAsset];
                        }

                        // doing sums
                        (*meanSum)[iAsset]      += (*mean)[iAsset];
                        (*varianceSum)[iAsset]  += (*variance)[iAsset];
                        iCovSum[iAsset]         += iCov[iAsset];

                    } // end for iAsset diagonal

                    // off-diagonal
                    for (iAsset = 0; iAsset < nbAssets; iAsset++) {

                        // for ease 'n speed
                        double* iCov    = (*thisCovariance)[iAsset];
                        double* iCovSum = (*covarianceSum)[iAsset];

                        for (jAsset = iAsset+1; jAsset < nbAssets; jAsset++) {
                            if( inst->useMean ) {
                                iCov[jAsset] = (*weightFactor)[iAsset] * (*weightFactor)[jAsset]
                                            * ( iCov[jAsset] - nbDates * (*mean)[iAsset] * (*mean)[jAsset] )
                                            / (nbDates-1);
                            }
                            else {
                                iCov[jAsset] = (*weightFactor)[iAsset] * (*weightFactor)[jAsset]
                                            * iCov[jAsset]
                                            / nbDates;
                            }

                            // doing sums
                            iCovSum[jAsset] += iCov[jAsset];
                            average         += inst->weights[iAsset] * inst->weights[jAsset] * iCov[jAsset];
                            if( inst->useVolInWeights ) {
                                denom       += (*weightFactorNorm)[iAsset] * (*weightFactorNorm)[jAsset];
                            }

                        }  // end for jAsset off-diagonal

                    } // end for iAsset off-diagonal

                    cleanPrice += average*norm/denom;
                    mcRuns++;

                //prices.add should be inside the !doingPast loop
                prices.add( inst->notional * (Maths::max(inst->floor,Maths::min(inst->cap,average*norm/denom)) - inst->strike) );
                } // end if !doingPast

            } // end if !isSpreadOfVarSwaps
        } catch (exception& e) {            
            throw ModelException(e, method);
        }
    } // end payoff

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel;
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();

        if (fwdStarting){
            interpLevel = 1.0;
        } else {
            interpLevel = pathGen->refLevel(iAsset, 0);
        }
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                 startDate,
                                                                 lastSimDate,
                                                                 fwdStarting));
        return reqarr;
    }

    virtual void recordExtraOutput(CControl*        control,
                                   Results*         results,
                                   const IMCPrices& prices) const {
        static const string method("CorrCovMC::recordExtraOutput");
        try {
            int iAsset, jAsset;
            double* iCov;
            CDoubleArray&  meanReturn = *meanSum;
            CDoubleArray&  varReturn  = *varianceSum;
            CDoubleMatrix& covReturn  = *covarianceSum;

            // have to do averaging myself...
            if( inst->isSpreadOfVarSwaps ) {
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    varReturn[iAsset] = singleVarTotalSum[iAsset] / mcRuns;
                }
            }
            else {
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // annualise!
                    meanReturn[iAsset] /= mcRuns * dt;
                    varReturn[iAsset]  /= mcRuns * dt;

                    // for ease 'n speed
                    iCov = covReturn[iAsset];
                    for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                        iCov[jAsset] /= mcRuns * (inst->useVolInWeights||!inst->isCorr ? dt : 1);
                    }
                }
            }

            results->storeScalarGreek(cleanPrice/mcRuns,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName("cleanPrice")));
            results->storeGreek((IObjectSP)meanSum,
                                Results::DEBUG_PACKET,
                                OutputNameSP(new OutputName("mean")));
            results->storeGreek((IObjectSP)varianceSum,
                                Results::DEBUG_PACKET,
                                OutputNameSP(new OutputName("variance")));
            results->storeGreek((IObjectSP)covarianceSum,
                                Results::DEBUG_PACKET,
                                OutputNameSP(new OutputName("covariance")));
            results->storeScalarGreek(basketVarTotalSum/mcRuns,
                                Results::DEBUG_PACKET,
                                OutputNameSP(new OutputName("basketVar")));

            // calc realized correlation so far (not include between last sample to valueDate)
            // NB : CorrSwap is not able to decompose to past + future.  This is just for
            // checking the data of past sampling.
            // only called when past samples is more than 2 to avoid dividing 0.
            if (control) {
                OutputRequest* request =
                    control->requestsOutput(OutputRequest::CORR_REALIZED);
                if (request) {
                    try {
                        double realizedCorr = 0.0;
                        if (nbDatesPast>1 && inst->isCorr && !inst->useVolInWeights){
                            // diagonal
                            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                                double* iCov        = (*covarianceSoFar)[iAsset];
                                (*variance)[iAsset] = iCov[iAsset] / nbDatesPast;
                                if(Maths::isZero((*variance)[iAsset])) {
                                    (*weightFactor)[iAsset]     = 0.0;
                                } else {
                                    (*weightFactor)[iAsset]     = 1.0 / sqrt((*variance)[iAsset]);
                                }
                                iCov[iAsset]                = 1.;
                            }
                            // off-diagonal
                            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                                double* iCov    = (*covarianceSoFar)[iAsset];
                                for (jAsset = iAsset+1; jAsset < nbAssets; jAsset++) {
                                    iCov[jAsset] = (*weightFactor)[iAsset] * (*weightFactor)[jAsset]
                                                    * iCov[jAsset]
                                                    / nbDatesPast;

                                    realizedCorr += inst->weights[iAsset] * inst->weights[jAsset] * iCov[jAsset];
                                }
                            }
                            realizedCorr *= norm;                    
                        } 
                        results->storeRequestResult(request, realizedCorr); 
                    } catch (exception& e) {
                        IObjectSP error(new Untweakable(e));
                        results->storeRequestResult(request, error);
                    }
                }

                /** output request to compute "corr in future" = average implied corr at maturity */
                request = control->requestsOutput(OutputRequest::CORR_IMPLIED);
                if (request) {
                    try {
                        /** special treatment if fwdstarting -- only crucial in presence of CorrTS */
                        DateTime valueDate = getToday(); 
                        DateTime startDate = inst->averageOutDates->start(); 
                        bool fwdStarting = startDate.isGreater(valueDate);
                        
                        DateTimeArray toDates(2);
                        toDates[0] = valueDate.max(startDate); // if started, then toDates[0] = valueDate
                        toDates[1] = inst->averageOutDates->end();
                        
                        if (!toDates[1].isGreater(valueDate)) { // no contribution, since expired
                            results->storeRequestResult(request, 0.0); 
                        } else {
                            /** two cases -- yes CorrTS & no CorrTS */
                            MultiAsset* mAsset = dynamic_cast<MultiAsset*>(inst->assets.get());
                            if (!mAsset) { // just in case sthg goes wrong
                                throw ModelException("Internal error. This shouldn't happen.");
                            }
                            bool isCorrTS = (mAsset->asMultiFactors()->getCorrTermObjArray().size() > 0);
                            
                            if (isCorrTS) { /** YES CorrTS */
                                DoubleMatrix fwdVarAtDates(0,0); // dummy, since not necessary    
                                CorrTermDataSP corrTermData = 
                                    CorrelationTerm::getCorrelationTermSqueezesAndExpiries (
                                        valueDate, 
                                        nbAssets,
                                        mAsset->asMultiFactors()->getCorrObjArray(), 
                                        mAsset->asMultiFactors()->getCorrTermObjArray()); 
                                BoolArray dummy((nbAssets * nbAssets - nbAssets) / 2); // dummy
                                TimeMetricArray timeMetricArray = mAsset->asMultiFactors()->getTimeMetricArray();
                                DoubleMatrix deltaT(nbAssets,3);                            
                                for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                                    deltaT[iAsset][0] = // zero as soon as started
                                        timeMetricArray[iAsset]->yearFrac(valueDate, toDates[0]);
                                    deltaT[iAsset][1] = // thats the important one
                                        timeMetricArray[iAsset]->yearFrac(valueDate, toDates[1]);
                                    deltaT[iAsset][2] = // usually, deltaT[1] = deltaT[2]
                                        deltaT[iAsset][1]-deltaT[iAsset][0];
                                    if(Maths::isZero(deltaT[iAsset][2])) { // sthg goes wrong, shouldnt happen
                                        throw ModelException("Division by zero,");
                                    }
                                }                            
                                DoubleMatrixArraySP spotCorrelAtMat = 
                                    CorrelationTerm::CorrelationTermMatrix(valueDate,
                                                                        toDates, 
                                                                        DoubleMatrix(0,0),   // dummy
                                                                        timeMetricArray,
                                                                        corrTermData, 
                                                                        CorrelationTerm::toLast, 
                                                                        CorrelationTerm::EIGEN_VALUE_FLOOR, 
                                                                        CorrelationTerm::MAX_SQ_ERROR,
                                                                        DoubleArray(2),      // dummy
                                                                        dummy);              // dummy
                                double futureCorr = 0.0, sumWeights = 0.0;
                                for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                                    for (int jAsset= iAsset+1; jAsset<nbAssets; jAsset++) {
                                        double thisCorrel = 0.0;
                                        if (isCorrTS) {
                                            thisCorrel = 
                                                (*(*spotCorrelAtMat)[1])[iAsset][jAsset] 
                                                    * sqrt(deltaT[iAsset][1] * deltaT[jAsset][1])
                                                - (*(*spotCorrelAtMat)[0])[iAsset][jAsset]
                                                    * sqrt(deltaT[iAsset][0] * deltaT[jAsset][0]);
                                            thisCorrel /= sqrt(deltaT[iAsset][2] * deltaT[jAsset][2]);
                                        } else {
                                            thisCorrel = (*(*spotCorrelAtMat)[1])[iAsset][jAsset];
                                        }
                                        futureCorr += inst->weights[iAsset] 
                                                    * inst->weights[jAsset] 
                                                    * thisCorrel;
                                        sumWeights += inst->weights[iAsset] *inst->weights[jAsset];
                                    }
                                }
                                if (Maths::isPositive(sumWeights)) { // zero in one-asset case, thus, shouldnt happen
                                    futureCorr /= sumWeights;
                                }
                                results->storeRequestResult(request, futureCorr); 
                            
                            } else { /** NO CorrTS */
                                CDoubleMatrixConstSP correlMatrix = 
                                    mAsset->asMultiFactors()->factorsCorrelationMatrix(); 
                                double correl = 0.0, sumWeights = 0.0;
                                for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                                    for (int jAsset= iAsset+1; jAsset<nbAssets; jAsset++) {
                                        correl += inst->weights[iAsset]
                                               * inst->weights[jAsset] 
                                               * (*correlMatrix)[iAsset][jAsset];
                                       sumWeights += inst->weights[iAsset] *inst->weights[jAsset];
                                    }
                                }
                                if (Maths::isPositive(sumWeights)) { // zero in one-asset case, thus, shouldnt happen
                                    correl /= sumWeights;
                                }
                                results->storeRequestResult(request, correl); 
                            }
                        }
                    } catch (exception& e) {
                        IObjectSP error(new Untweakable(e));
                        results->storeRequestResult(request, error);
                    }
                }
            }
        } catch (exception& e) {            
            throw ModelException(e, method);
        }
    }

    /** calculate historical var */
    void historicalVar(
                    const IPathGenerator*   pathGen,
                    double&                 basketVar,
                    DoubleArray&            singleVar,
                    int&                    numReturns) const;
    /** calculate future var */
    void futureVar(
                    const IPathGenerator*   pathGen,
                    double&                 basketVar,
                    DoubleArray&            singleVar,
                    int&                    numReturns) const;
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CorrCov::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    if( isSpreadOfVarSwaps ) {
        // set up new date list with all past dates and only one future date
        int nbDates         = averageOutDates->size();
        int nbFutureDates   = valueDate.numFutureDates(*(averageOutDates->dates()));
        int iDate;
        DateTimeArray averageOutDatesUse(nbDates + 1 - nbFutureDates);
        for( iDate=0; iDate<nbDates-nbFutureDates; iDate++ ) {
            averageOutDatesUse[iDate] = (*averageOutDates)[iDate];
        }
        averageOutDatesUse[iDate] = averageOutDates->end();
        simSeries->addDates(averageOutDatesUse);
    }
    else {
        simSeries->addDates(obsMap->getModellingDates());
    }
    return new CorrCovMC(this, simSeries);
}

CClassConstSP const CorrCov::TYPE = CClass::registerClassLoadMethod(
    "CorrCov", typeid(CorrCov), CorrCov::load);

// * for class loading (avoid having header file) */
bool CorrCovSwapLoad() {
    return (CorrCov::TYPE != 0);
}

/** calculate future var */
void CorrCovMC::futureVar(
                const IPathGenerator*   pathGen,
                double&                 basketVar,
                DoubleArray&            singleVar,
                int&                    numReturns) const {
    static const string routine("CorrCov::historicalVol");
    try {
        int    beginIdx = pathGen->begin(0); // 0 <=same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset;

        numReturns = endIdx-beginIdx;
        if( numReturns<1 ) {
            throw ModelException(routine,
                                 "Number of returns = "+
                                 Format::toString(numReturns)+
                                 " greater 0 assumed.");
        }
        if( pathGen->doingPast() ) {
            throw ModelException(routine,
                                 "Doing future assumed.");
        }
        if( inst->weights.size()!=nbAssets || singleVar.size()!=nbAssets ) {
            throw ModelException(routine,
                                 "Lengths of weights = "+
                                 Format::toString(inst->weights.size())+
                                 " and singleVar = "+
                                 Format::toString(inst->weights.size())+
                                 " assumed to be equal to nbAssets = "+
                                 Format::toString(nbAssets) );
        }

        // do asset var
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            // var = 2/T * ( E[log(S_T)] - log(S_0) - F_T/S_0 )
            singleVar[iAsset]
                = twoOverTime * ( (*logFwd)[iAsset] - log( pathGen->Path(iAsset, 0)[endIdx-1] ) );
        }

        // do basket, use only last fixing
        double bPath = 0;
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            bPath += inst->weights[iAsset]*pathGen->Path(iAsset, 0)[endIdx-1];
        }

        // do basket var
        basketVar
            = twoOverTime * ( logbFwd - log( bPath ) );

    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}
/** calculate historical var */
void CorrCovMC::historicalVar(
                const IPathGenerator*   pathGen,
                double&                 basketVar,
                DoubleArray&            singleVar,
                int&                    numReturns) const {
    static const string routine("CorrCov::historicalVol");
    try {
        int    beginIdx = pathGen->begin(0); // 0 <=same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, iDate;

        numReturns = endIdx-beginIdx;
        if( numReturns<1 ) {
            throw ModelException(routine,
                                 "Number of returns = "+
                                 Format::toString(numReturns)+
                                 " greater 0 assumed.");
        }
        if( !pathGen->doingPast() ) {
            throw ModelException(routine,
                                 "Doing past assumed.");
        }
        if( inst->weights.size()!=nbAssets || singleVar.size()!=nbAssets ) {
            throw ModelException(routine,
                                 "Lengths of weights = "+
                                 Format::toString(inst->weights.size())+
                                 " and singleVar = "+
                                 Format::toString(inst->weights.size())+
                                 " assumed to be equal to nbAssets = "+
                                 Format::toString(nbAssets) );
        }

        // do asset var
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            // for ease 'n speed
            const double* iPath = &(pathGen->Path(iAsset, 0)[0]);
            const double spot = getMultiFactors()->assetGetSpot(iAsset);
            const double refLevel = pathGen->refLevel(iAsset, 0);
            const double  iMean = inst->useMean ? log( spot / refLevel ) / numReturns : 0;

            double logRelative;
            double sumSqr = 0.0;

            for (iDate = beginIdx; iDate < endIdx; iDate++) {
                logRelative = log(iPath[iDate]/(iDate>0 ? iPath[iDate-1] : refLevel)) - iMean;
                sumSqr      += logRelative * logRelative;
            }

            // ideally we want to calculate historic vol from the historic
            // samples and the future vol from the vol surface between value
            // date and the last sample date.  However, this leaves out the
            // vol between the last historic sample and value date.  To cover
            // this we set the sample after value date to the current spot
            // and calculate the historic vol up the sample after value date.
            // This essentially captures the vol between last historic sample
            // sample and value date. Set sample on or after value date to spot

            logRelative = log(spot/iPath[iDate-1]) - iMean;
            sumSqr      += logRelative * logRelative;

            // need at least 2 returns to get mean
            double var;
            if (!inst->useMean || numReturns < 2) {
                var = sumSqr/numReturns;
            }
            else {
                var = sumSqr/(numReturns-1);
            }
            singleVar[iAsset] = var;
        }

        // do basket
        DoubleArray bPath(endIdx);
        double basketRefLevel = 0.0;
        double basketSpot = 0.0;
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            const double* iPath = &(pathGen->Path(iAsset, 0)[0]);
            basketRefLevel  += inst->weights[iAsset]*pathGen->refLevel(iAsset, 0);
            basketSpot      += inst->weights[iAsset] *
                getMultiFactors()->assetGetSpot(iAsset);
            for (iDate = beginIdx; iDate < endIdx; iDate++) {
                bPath[iDate] += inst->weights[iAsset]*iPath[iDate];
            }
        }

        // do basket var
        const double  bMean = inst->useMean ?
            log( basketSpot / basketRefLevel ) / numReturns : 0;

        double logRelative;
        double sumSqr = 0.0;

        for (iDate = beginIdx; iDate < endIdx; iDate++) {
            logRelative = log( bPath[iDate]/(iDate==0 ? basketSpot : bPath[iDate-1]) ) - bMean;
            sumSqr      += logRelative * logRelative;
        }

        logRelative = log(basketSpot/bPath[iDate-1]) - bMean;
        sumSqr      += logRelative * logRelative;

        // need at least 2 returns to get mean
        if (!inst->useMean || numReturns < 2) {
            basketVar = sumSqr/numReturns;
        }
        else {
            basketVar = sumSqr/(numReturns-1);
        }

    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

DRLIB_END_NAMESPACE
