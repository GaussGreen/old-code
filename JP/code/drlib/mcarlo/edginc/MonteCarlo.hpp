//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MonteCarlo.hpp
//
//   Description : Monte Carlo model
//
//   Date        : May 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MONTECARLO_HPP
#define EDR_MONTECARLO_HPP
#include "edginc/Model.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/MCPricing.hpp"  // Could remain in .cpp directory?
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/LastSensDate.hpp"

DRLIB_BEGIN_NAMESPACE

class IMCQuickXGamma;
class IMCQuickGreeks;
class MonteCarloHelper;

class MCARLO_DLL MonteCarlo : public CModel,
                   virtual public IModelFamily,
                   virtual public LastProductSensDate {
public:
    static CClassConstSP const TYPE;
    friend class MonteCarloHelper;
    //friend class IMCProduct;

    /** Default maximum storage space used for caching paths */
    static const int MAX_STORAGE;
    static const string CACHE_PATH_DEFAULT;
    static const string CACHE_PROD_EXTENSION;

    static const string QUICK_GREEK_DEFAULT;
    static const string QUICK_GREEK_AND_X_GAMMA;
    static const string QUICK_GREEK_NO_X_GAMMA;
    static const string QUICK_X_GAMMA_ONLY;
    static const string LR_GREEK_DEFAULT;
    static const string SIM_MODE_DEFAULT;
    static const string DEBUG_GREEK_STD_ERR;
    virtual ~MonteCarlo();
 
    IMCPathConfigSP getPathConfig() const;
    
    /** override of main control - calculates price and sensitivities in
        blocks to limit memory usage. */
    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
                                     CControl* control);

     /** Implementation of method in IModel */
    virtual void Price(CInstrument*  instrument,
                       Control*      control, 
                       CResults*     results);

    /** Invoked after instrument has got its market data. Passes on the
        call to the MCPathConfig object */
    virtual void getMarket(const MarketData* market,
                           IInstrumentCollectionSP instrument);


    /** Override default clone to do shallow copy of 'origPrices' */
    IObject* clone() const;

    /** Override default createMDF in order to set the MDF returned
     * by pathConfig */
    virtual MarketDataFetcherSP createMDF() const;

    /** indicates whether the MonteCarlo infrastructure can support 
        Vega Matrix for supplied instrument */
    virtual bool vegaMatrixSupported(CInstrument*  instrument) const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns MCPathConfig::wantsRiskMapping() from the pathConfig.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** returns all strikes on the vol surface to which 
        the supplied instrument is sensitive. Use vegaMatrixSupported to see
        if this routine should work */
    DoubleArraySP getSensitiveStrikes(
        CInstrument*             instrument,
        const OutputNameConstSP& outputName) const;

    /** when to stop tweaking */
    DateTime endDate(const CInstrument*  instrument,
                     const Sensitivity*  sensitivity) const;

    virtual void validatePop2Object();

    /** Returns true if use of the state var framework is requested */
    bool stateVarUsed() const;

    /** Whether we are using stateless payoff */
    bool inStatelessMode() const;

    // For the IModelFamily interface
    const string* getFamilyName() const;
    IMCProduct* createProduct(const CInstrument*  instrument) const;
    void initializePathConfigPostPricing();

    // Accessors
    bool cachedQckGreeksOfType( int type ) const;
    IMCPricesSP getOrigPrices() const;
    int getNbIters(bool bypass = false) const;
    int getNbSubSamples() const;
    int getCachedQckGreeksModes() const;
    MCPricingSP getPricing() const;
    bool getCachedCachePaths() const; // FIX name of this!!

    void setOrigPrices( IMCPricesSP prices );
    void enableCachedQckGreeksMode(int mode);
    void disableCachedQckGreeksMode(int mode);

    // The following typdef is here to avoid breaking the public interface
    // It must be public since products extend this interface.
	typedef IMCIntoProduct IIntoProduct;

protected:
    MonteCarlo(CClassConstSP clazz);

    IMCPathConfigSP         pathConfig;
    int                     nbIter;
    int                     nbSubSamples;
    string                  quickGreeks; // yes/no/default/...
    string                  lrGreeks; // yes/no/default/...
    bool                    useStateVars;
    bool                    useCcyBasis;
    string                  simMode; // slice/path/pathOpt
    string                  cachePaths;  // yes/no/default/number
    int                     cachedQckGreeks; // bitwise values (transient)
    bool                    cachedCachePaths; // transient
    // following field not registered. clone() method shallow copies it
    IMCPricesSP     origPrices; /* from original pricing run $unregistered */
    /* Field is a clone of pathConfig taken when pricing (rather than greeks)
       and the sim is over - it is the pathConfig to use on the next block.
       This highlights the fact that we're mixing up user data with run time
       data */
    IMCPathConfigSP pathConfigPostPricing;

    //used to set pdfBoundaryProb in CVolProcessedBS::defaultStrikes() via a hard-coded scenario
    int                     pdfBoundaryFactor;    

private:
//    friend class Pricing;
//    class QGPricing;
//    friend class QGPricing;
//    class QXGamma;
//    friend class QXGamma;
//    class SubGamma;
//    friend class SubGamma;
//    class LRPricing;
//    friend class LRPricing;
//    class MCSlicePricing;
//    friend class MCSlicePricing;
//    class MCPathPricing;
//    friend class MCPathPricing;
//    class MCPastPricing;
//    friend class MCPastPricing;
//    class SimDates;
//    friend class SimDates;

    // internally configured specialised class for controlling pricing
    MCPricingSP               pricing;  // not registered: shallow copied $unregistered
    /* for reflection */
    MonteCarlo();

    void calculateSubGamma(
        IMCQuickXGamma*            qckXGamma,
        const MCPathGeneratorSP&  futurePathGen,
        Control*                            control, 
        CInstrument*                        instrument,
        CResults*                           results,
        const ScalarShiftArray&             theShifts);

     ResultsSP createXGammaResults(
        const string&                 ccyName,
        const OutputNameSP&           name1,
        const OutputNameSP&           name2,
        double                        subGamma12,
        double                        subGamma21,
        const ScalarShiftSP&          theShift);

    void priceXGamma(IMCQuickXGamma*            qckXGamma,
                     const MCPathGeneratorSP&  futurePathGen,
                     CInstrument*                        instrument,
                     Control*                            control, 
                     CResults*                           results,
                     const Sensitivity*                  xGammaSens);


    void quickGreeksChoice(Control*                 control,
                           IMCQuickGreeks* qckGreeksImpl,
                           IMCQuickXGamma* qckXGammaImp);

    /** Used to set pdfBoundaryProb for use in CProcessedVolBS */
    void updateMarketBoundaryProb(  CInstrumentSP instr, 
                                    TweakGroupSP tweakGroup, 
                                    CClassConstSP clazz);


protected:
    void createPricingObject(Control*   control,
                             IMCProduct* product);
    
};

DRLIB_END_NAMESPACE

#endif
