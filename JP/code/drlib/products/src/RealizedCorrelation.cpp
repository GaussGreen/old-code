//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RealizedCorrelation.cpp
//
//   Description : Computes the sample correlation matrix between N factors
//
//   Date        : December 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"


DRLIB_BEGIN_NAMESPACE

/** Computes the sample correlation matrix between N factors */
class RealizedCorrelation: public GenericNFBase, 
                           virtual public IMCIntoProduct{
private:
    /// fields ////////
    bool          isRealizedCorrelation;
    DateTimeArray averageOutDates;

public:
    static CClassConstSP const TYPE;
    friend class RealizedRealizedCorrelationMC;

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return averageOutDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    RealizedCorrelation(): GenericNFBase(TYPE) {} // for reflection
    RealizedCorrelation(const RealizedCorrelation& rhs); // not implemented
    RealizedCorrelation& operator=(const RealizedCorrelation& rhs); // not implemented

    static IObject* defaultRealizedCorrelation(){
        return new RealizedCorrelation();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RealizedCorrelation, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultRealizedCorrelation);
        FIELD(isRealizedCorrelation,     "is it a RealizedCorrelation");
        FIELD(averageOutDates,  "averageOutDates");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for RealizedCorrelation */
class RealizedRealizedCorrelationMC : public IMCProduct,
                                      virtual public IMCProductLN {
private:
    const RealizedCorrelation* inst;       // reference to original instrument
    int                        nbAssets;   // for convenience
    CDoubleMatrixSP            covariance;
    CDoubleArraySP             variance;
    CDoubleArraySP             mean;
    int                        mcRuns;
    
public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    RealizedRealizedCorrelationMC(const RealizedCorrelation* inst,
                                  const SimSeriesSP&         simSeries):
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
        mcRuns(0) { 
            covariance = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));
            variance   = CDoubleArraySP(new CDoubleArray(nbAssets));
            mean       = CDoubleArraySP(new CDoubleArray(nbAssets));
        }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int    iAsset, jAsset;
        int    endIdx   = pathGen->end(0);

        int date = endIdx - 1;
        CDoubleArray& meanReturn = *mean;

        // Any past samples
        if (pathGen->doingPast()) {
            throw ModelException("Payoff does not support past values");
        }

        // Continue using sumOut as a working area - now for perfs
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            double logReturn1 = log(pathGen->Path(iAsset, 0)[date] /
                                    pathGen->refLevel(iAsset, 0));
            meanReturn[iAsset] += logReturn1;
            // for ease 'n speed
            double* iCov = (*covariance)[iAsset];
            for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                iCov[jAsset] += log(pathGen->Path(jAsset, 0)[date] /
                                   pathGen->refLevel(jAsset, 0)) * logReturn1;
            }
        }
        
        prices.add(0.0); 
        mcRuns++;
    }

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

    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const {
        int iAsset, jAsset;
        double* iCov;
        CDoubleArray&  meanReturn = *mean;
        CDoubleArray&  varReturn  = *variance;
        CDoubleMatrix& covReturn  = *covariance;
        
        // Continue using sumOut as a working area - now for perfs
        for (iAsset = 0; iAsset < nbAssets; iAsset++) {
            meanReturn[iAsset] = meanReturn[iAsset] / mcRuns;
            varReturn[iAsset]  = (covReturn[iAsset][iAsset] - mcRuns * Maths::square(meanReturn[iAsset])) / 
                                 (mcRuns - 1);
        }

        for (iAsset = 0; iAsset < nbAssets; iAsset++) {    
            // for ease 'n speed
            iCov = covReturn[iAsset];
            for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                iCov[jAsset] = (iCov[jAsset] - mcRuns * meanReturn[iAsset] * meanReturn[jAsset]) / (mcRuns - 1);
            }
        }

        if(inst->isRealizedCorrelation) {
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {    
                // for ease 'n speed
                iCov = covReturn[iAsset];
                for (jAsset = iAsset; jAsset < nbAssets; jAsset++) {
                    iCov[jAsset] = iCov[jAsset] / sqrt(varReturn[iAsset] * varReturn[jAsset]);
                }
            }
        }

        results->storeGreek(IObjectSP(copy(mean.get())),
                            Results::DEBUG_PACKET,
                            OutputNameSP(new OutputName("mean")));
        results->storeGreek(IObjectSP(copy(variance.get())),
                            Results::DEBUG_PACKET,
                            OutputNameSP(new OutputName("variance")));
        results->storeGreek(IObjectSP(copy(covariance.get())),
                            Results::DEBUG_PACKET,
                            OutputNameSP(new OutputName("covariance")));
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RealizedCorrelation::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new RealizedRealizedCorrelationMC(this, simSeries);
}

CClassConstSP const RealizedCorrelation::TYPE = CClass::registerClassLoadMethod(
    "RealizedCorrelation", typeid(RealizedCorrelation), RealizedCorrelation::load);

// * for class loading (avoid having header file) */
bool RealizedCorrelationLoad() {
    return (RealizedCorrelation::TYPE != 0);
}

DRLIB_END_NAMESPACE
