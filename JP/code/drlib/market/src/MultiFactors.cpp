//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiFactors.cpp
//
//   Description : 
//
//   Date        : Oct 01
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_MULTIFACTORS_CPP
#include "edginc/MultiFactors.hpp"
#include "edginc/AsMultiFactors.hpp"
#include "edginc/GeneralAsset.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

/** Helper for centralised sampling/isda adjustment in the n-factor case
gets a list of dates for which the given sample will finally be known for 
each asset. Results are put into obsDates array and the BoolArray tracks
those assets for which date is omitted*/
void IMultiFactors::observationDates(const DateTime&           sampleDate,
                                const ObservationSourceArray&  sources,
                                const SamplingConvention*      sampleRule,
                                DateTimeArray&                 obsDates,
                                BoolArray&                     sampled) const {
    static const string routine("IMultiFactors::observationDates");
    try {
        // ISS - note this logic is fine if all assets roll in same
        // direction. Might need revisiting if we have asset level
        // overrides i.e. some roll fwd, some back
        bool rollAll = sampleRule->rollAssetsTogether();
        sampled[0] = (getAsset(0)).observationDate(sampleDate, 
                                                   sources[0].get(),
                                                   sampleRule,
                                                   &(obsDates[0]));
        bool omit = (rollAll && !sampled[0]);
        for (int i = 1; i < NbAssets() && !omit; i++) {
            sampled[i] = (getAsset(i)).observationDate(rollAll ? obsDates[i-1] : sampleDate, 
                                                       sources[i].get(),
                                                       sampleRule,
                                                       &(obsDates[i]));
            omit = (rollAll && !sampled[i]);
        }
        // if assets roll together we need to go back and adjust things
        if (rollAll) {
            if (omit) { //if one asset got omitted they all get omitted
                for (int i = 0; i < NbAssets(); i++) {
                    sampled[i] = false;
                }
            } else {
                // all assets move to the last date
                for (int i = 0; i < NbAssets()-1; i++) {
                    obsDates[i] = obsDates[NbAssets()-1];
                }
            }
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** Helper for parametrised (SPI) dates in the n-factor case
return a bool indicating whether any of the assets or their components
have a sampling holiday on the given date */
bool IMultiFactors::isHoliday(const DateTime&                sampleDate,
                              const ObservationSourceArray&  sources) const {
    static const string routine = "IMultiFactors::isHoliday";
    try {
        for (int i=0; i<NbAssets(); ++i) {
            if (getAsset(i).isHoliday(sampleDate, sources[i].get())) {
                return true;
            }
        }
        return false;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IAsMultiFactors::IAsMultiFactors(){}
IAsMultiFactors::~IAsMultiFactors(){}
IMultiMarketFactors::~IMultiMarketFactors(){}

static void multiFactorsLoad(CClassSP& clazz){
    REGISTER_INTERFACE(IMultiFactors, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IMultiFactors::TYPE = CClass::registerInterfaceLoadMethod(
    "IMultiFactors", typeid(IMultiFactors), multiFactorsLoad);

void IMultiMarketFactors::load(CClassSP& clazz){
    REGISTER_INTERFACE(IMultiMarketFactors, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IMultiMarketFactors::TYPE =
CClass::registerInterfaceLoadMethod(
    "IMultiMarketFactors", typeid(IMultiMarketFactors), load);

class IMultiMarketFactors::AsMulti: public CObject, 
                           public virtual IMultiMarketFactors{
    IMarketFactorConstSP marketFactor;
public:
    static CClassConstSP const TYPE; // in MultiFactors.cpp

    AsMulti(IMarketFactorConstSP marketFactor): 
        CObject(TYPE), marketFactor(marketFactor){}

    virtual ~AsMulti(){}

    /** Pull out the necessary (eg assets, correlations) from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market){
        // we'll assume market has already been obtained
    }
    /** Returns the number of factors in this collection (this does not include
        MarketFactors which are contained within top level MarketFactors) */
    virtual int numFactors() const{
        return 1;
    }

    /** For backwards compatibility: same as numFactors */
    virtual int NbAssets() const{
        return 1;
    }

    //// Returns the correlations between the immediate factors.
    //// This is numFactors() x numFactors() 
    virtual CDoubleMatrixConstSP factorsCorrelationMatrix() const{
        CDoubleMatrixSP matrix(new CDoubleMatrix(1,1));
        (*matrix)[0][0] = 1.0;
        return matrix;
    }

    /** Validate that the market data and instrument are consistent with
        each other etc */
    virtual void crossValidate(const DateTime&       startDate,
                               const DateTime&       valueDate,
                               const YieldCurve*     discCcy,
                               const CInstrument*    instrument) const{
    }

    /** Returns the name of the specified 'market factor' */
    virtual string getName(int index) const{
        return marketFactor->getName();
    }
    
    /** Returns the MarketFactor with the specifed index */
    virtual IMarketFactorConstSP getFactor(int index) const{
        return marketFactor;
    }

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    virtual IntArray getSensitiveFactors(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const{
        IntArray factors(1,0);
        return factors;
    }

    /** A utility to record 'forwards' at maturity (for each factor) as
        appropriate for each factor (if applicable for that factor) */
    virtual void recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const{
        if (IGeneralAsset::TYPE->isInstance(marketFactor.get())){
            // unclear exactly what methods we want on MarketFactor
            // at the moment
            IGeneralAsset* genAsset = STATIC_CAST(IGeneralAsset,
                                                  marketFactor.get());
            genAsset->recordFwdAtMat(request, results, maturityDate);
        }
    }

private:
    AsMulti(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new AsMulti();
    }

    static void load(CClassSP& clazz){
        REGISTER(AsMulti, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IMultiMarketFactors);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(marketFactor, "the market factor");
    }
};
CClassConstSP const IMultiMarketFactors::AsMulti::TYPE = 
CClass::registerClassLoadMethod(
    "IMultiMarketFactors::AsMulti", typeid(AsMulti), load);

/** Make a single market factor look like a IMarketFactor */
IMultiMarketFactorsSP IMultiMarketFactors::asMulti(
    IMarketFactorConstSP marketFactor){
    return IMultiMarketFactorsSP(new AsMulti(marketFactor));
}

DRLIB_END_NAMESPACE


