//----------------------------------------------------------------------------
//
//   Description : Test product for the QSRM CID model
//
//   Author      :
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp" 
//#include "edginc/GenericNFactor.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/MCTestingCIDMatricesMC.hpp"
#include "edginc/SimpleEquity.hpp"

#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"

#include "edginc/MCTestingCIDMatrices.hpp"

DRLIB_BEGIN_NAMESPACE


// for class loading
bool MCTestingCIDMatricesLoad();

void MCTestingCIDMatrices::Validate()
{
    const string & method = "MCTestingCIDMatrices::Validate";
    // validate market data here:    

}
void MCTestingCIDMatrices::validatePop2Object()
{
    const string & method = "MCTestingCIDMatrices::validatePop2Object";
    // validate product details (not market data related) here:
}

IObject* MCTestingCIDMatrices::defaultMCTestingCIDMatrices()
{
    return new MCTestingCIDMatrices();
}

MCTestingCIDMatrices::MCTestingCIDMatrices() :
        GenericNFactor( TYPE )
{}

// void MCTestingCIDMatrices::GetMarket(const IModel* model, const CMarketDataSP market)
// {
//   // if we were to put tie an asset with this product, we would do something like
//   cdsCurveWrapper.getData(model, market);
//   yieldCurveWrapper.getData(model, market);

//   // for quanto adjustment, only if necessary...
//   if (!fxWrapper.isEmpty())
//  fxWrapper.getData(model, market);
// }

/** Invoked when Class is 'loaded' */
void MCTestingCIDMatrices::load( CClassSP& clazz )
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER( MCTestingCIDMatrices, clazz );
    SUPERCLASS( GenericNFactor );
    IMPLEMENTS( IMCIntoProduct );

    EMPTY_SHELL_METHOD( defaultMCTestingCIDMatrices );
    //        FIELD(cdsCurveWrapper, "cds curve wrapper");
    //         FIELD(yieldCurveWrapper, "yield curve wrapper");
    //         FIELD(fxWrapper, "fx wrapper");
    //         FIELD_MAKE_OPTIONAL(fxWrapper);
    FIELD( today, "today" );
    FIELD( maturityDates, "maturity dates of the survival prob term structures" );
    FIELD( measureDate, "measurement date" );
    FIELD( doLog, "compute log" );
}

CClassConstSP const MCTestingCIDMatrices::TYPE = CClass::registerClassLoadMethod(
            "MCTestingCIDMatrices", typeid( MCTestingCIDMatrices ), MCTestingCIDMatrices::load );

// for class loading
bool MCTestingCIDMatricesLoad()
{
    return ( MCTestingCIDMatrices::TYPE != 0 );
}

/// GenericNFactor contains a list of assets; we parse this list using SimpathIInstrument facilities;
void MCTestingCIDMatrices::classifyAssets( vector<YieldCurveConstSP>& ycAsset,
        vector<FXAssetConstSP>& fxAsset,
        vector<ICDSParSpreadsConstSP>& cdsAsset,
        //      vector<SimpleEquityConstSP>& eqAsset,
        //      vector<EnergyFuturesCurveConstSP>& enrgAsset,
        //        vector<IBasisIndexCurveConstSP>& basisAsset,
        map<string, int>& fxAssetIDMap ) const
{
    vector<SimpleEquityConstSP> eqAsset;
    vector<EnergyFuturesCurveConstSP> enrgAsset;
    vector<IBasisIndexCurveConstSP> basisAsset;
    
    // Call SimpathIInstrument to classify market objects
    SimpathIInstrument::classifyAssets( this->assets, ycAsset, fxAsset, cdsAsset, eqAsset, enrgAsset, basisAsset, fxAssetIDMap );
    
    QLIB_VERIFY( !cdsAsset.empty(), "Need at least one credit" );

    QLIB_VERIFY( eqAsset.empty(), "ERROR: SimpleEquity is not supported" );
    QLIB_VERIFY( enrgAsset.empty(), "ERROR: EnergyFuturesCurve is not supported" );
    QLIB_VERIFY( basisAsset.empty(), "ERROR: IBasisIndexCurve is not supported" );
}

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* MCTestingCIDMatrices::createProduct( const MonteCarlo* model ) const
{
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries( new SimSeries( 1 ) ); // create empty one
    //simSeries->addDates(DateTimeArray(1, maturityDate));
    simSeries->addDates( maturityDates ); // need fix ??

    // this is required
    InstrumentSettlementSP instSettle( new CashSettlePeriod( 0 ) );
    
    vector<ICDSParSpreadsConstSP> cdsAsset;

    vector<YieldCurveConstSP> ycAsset;
    vector<FXAssetConstSP> fxAsset;
    map<string, int> fxAssetIDMap;

    classifyAssets( ycAsset, fxAsset, cdsAsset, fxAssetIDMap );

    return new MCTestingCIDMatricesMC(
               assets,  // IMultiMarketFactors*
               ycAsset,
               cdsAsset,
               fxAsset,
               fxAssetIDMap,
               discount.get(),  // member of GenericNFactor
               today,
               maturityDates,
               measureDate,
               doLog );
}






DRLIB_END_NAMESPACE
