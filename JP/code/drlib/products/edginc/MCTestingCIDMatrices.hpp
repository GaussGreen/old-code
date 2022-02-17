#ifndef MC_TESTING_CID_MATRICES
#define MC_TESTING_CID_MATRICES


#include "edginc/Instrument.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/DECLARE.hpp"

#include "edginc/SimpathIInstrument.hpp"

DRLIB_BEGIN_NAMESPACE
class PRODUCTS_DLL MCTestingCIDMatrices:
            public GenericNFactor,
            public virtual IMCIntoProduct
{

public:
    static CClassConstSP const TYPE;

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct( const MonteCarlo* model ) const; // see below

    //  virtual void GetMarket(const IModel* model, const CMarketDataSP market);
    virtual void Validate();

    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load( CClassSP& clazz );

    static IObject* defaultMCTestingCIDMatrices();

private:
    MCTestingCIDMatrices();
    
    // will reuse SimpathIInstrument::classifyAssets()
    // the vector<..SP> is not stored in the instrument, instead it is passed to product
    void classifyAssets(
        vector<YieldCurveConstSP>& ycAsset,
        vector<FXAssetConstSP>& fxAsset,
        vector<ICDSParSpreadsConstSP>& cdsAsset,
        //              vector<SimpleEquityConstSP>& eqAsset,
        //              vector<EnergyFuturesCurveConstSP>& enrgAsset,
        //              vector<IBasisIndexCurveConstSP>& basisAsset,
        map<string, int>& fxAssetIDMap ) const;

    // Set via reflection
    DateTime today;
    DateTimeArray maturityDates; // maturity dates of the survival prob term structures
    DateTime measureDate; // measurement date
    bool doLog; // compute log
};

DRLIB_END_NAMESPACE

#endif
