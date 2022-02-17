//----------------------------------------------------------------------------
//
//   Description : MC implementation of MCTestingCIDMatrices
//
//   Author      :
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"
#include "edginc/MCTestingCIDMatrices.hpp"

#include "edginc/MCTestingCIDMatricesMC.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
// #include "edginc/SVGenExpectedBasisFwdSpread.hpp"
// #include "edginc/SVExpEnergyFuture.hpp"

DRLIB_BEGIN_NAMESPACE
// sample class to accumulate complex payoff.
// we will have chance to postprocess after the simulation is finished and IMCProduct->recordExtraOutput is called

// This class captures everything we want to be passed for accumulation for this product as we generate path after path
// Exaples could be: vector of ESDFs between measurement date and all diffusion dates
// The internals of this class are improtant only to the relevant IMCPricesGeneric class that implements the actual accumulation (like MCMatrixPrices)

class MatrixPayoffEvent : public IPayoffEvent
{
    vector<vector<double> > esdf; // values of esdf[i][Tj] for all assets and maturity dates
public:
    MatrixPayoffEvent( const vector<vector<double> >& vals ) : esdf( vals )
    {}
    const vector<vector<double> >& getElems() const
    {
        return esdf;
    }
};


class MCMatrixPrices : public IMCPricesGeneric
{
    CDoubleMatrix m; // accumulator; convenient to have an IObject so we can just return it in getResult()

    void process( const MatrixPayoffEvent& mattr, double pathWeight )
    {
        const vector<vector<double> >& elems ( mattr.getElems() );

        QLIB_VERIFY( elems.size() == nbAssets, "Unexpected number of assets in input" );

        for ( size_t i = 0; i != nbAssets; ++i ) {
            QLIB_VERIFY( elems[ i ].size() == numMatDates, "MatrixPayoffEvent has size != numMatDates" );
            for ( size_t j = 0; j != numMatDates; ++j )
                m[ i ][ j ] += elems[ i ][ j ] * pathWeight;
        }
    }

    int nbIter;
    int nbSubSamples;
    int numMatDates;
    int nbAssets;
public:
    MCMatrixPrices(
        int NbIter,
        int NbSubSamples,
        int NbAssets,
        int NumMatDates // number of the maturity dates
    ) :
            nbIter( NbIter ),
            nbSubSamples( NbSubSamples ),
            numMatDates( NumMatDates ),
            nbAssets( NbAssets ),
            m( NbAssets, NumMatDates )
    {
        reset();
    }
    //////////////  IGenericPrices interface
    virtual void add
        ( const IPayoffEvent& ev, double weight )
    {
        const MatrixPayoffEvent & mattr =
            dynamic_cast<const MatrixPayoffEvent&> ( ev );
        process( mattr, weight );
    }

    // Will be used to return results
    virtual IObjectSP getResult() const
    {
        CDoubleMatrixSP res( new CDoubleMatrix( m ) );
        res->scale( 1.0 / nbIter );
        return res;
    }

    ////////////////// IGreeks interface ////////////
    virtual void reset()
    {
        m = CDoubleMatrix( nbAssets, numMatDates );
    }

    virtual int storagePerPath( IMCProduct* product ) const
    {
        return -1;
    }
    virtual void configureCache( const IntArray& changedAssets )
    {}

    ///////////////// IMCPrices interface ////////
public:
    /** Returns a deep copy of this object */
    virtual IMCPrices* clone() const
    {
        QLIB_VERIFY( 0, "Not implemented" );
        return NULL;
    }

protected:
    /** Ease cloning */
    virtual IMCPrices* emptyConstructor() const
    {
        return new MCMatrixPrices( 1, 1, 0, 0 );
    }
public:
    // to be called from MCProduct::recordExtraOutput
    void recordExtraOutput( CControl *control, Results * result ) const
    {
        // FIXME
    }
};
///////////////////////// end of MCMatrixPrices /////////////////

// Create custom Prices class
IMCPrices* MCTestingCIDMatricesMC::createOrigPrices( int nbIter,
        int nbSubSamples,
        int mode )
{
    size_t nbAssets = esdfGen.size();
    return new MCMatrixPrices( nbIter, nbSubSamples, nbAssets, nbMaturityDates );
}

/* At this point we know that we finished the simulation and are given a chance to output whatever we want.
   The prices here are what we created in createOrigPrices, so the actual type is MCMatrixPrices.
*/

void MCTestingCIDMatricesMC::recordExtraOutput( CControl* control,
        Results* results,
        const IMCPrices& prices ) const
{
    const MCMatrixPrices & matrixPrices = dynamic_cast<const MCMatrixPrices&>( prices );
    matrixPrices.recordExtraOutput( control, results );
}

/** Payoff now stores complex payoff via IHandlePayoffEvent interface of its prices object.
 * We could use static_cast to cast prices to IHandlePayoffEvent (as we created it in createOrigPrices() method)
 */

void MCTestingCIDMatricesMC::payoff(
    const IPathGenerator* pathGen,
    IMCPrices& prices )
{
    IGenericPrices & payoffHandler = dynamic_cast<IGenericPrices&>( prices );
    double weight = pathW->getWeight(0); 

    // only retreive the price at the measurement date along each path:
    if ( 0 /*dd->getDateOfDefault() <= today*/ )   // FIXME: revert when dd is impl.
    {
        // defaulted:
        //vector<double> esdfVals(nbMaturityDates, 0.0);
        //payoffHandler.add(MatrixPayoffEvent(esdfVals), weight); // TODO: modify this payoff!!
    } else {
        //             // not yet defaulted:
        //         double sum = 0.0;
        //         for (int i = 0; i < nbMaturityDates; ++i)
        //             sum += esdf->getExpSDF(i) * rr->getRecoveryRate(0); // TODO: Modify this! This payout is just for testing
        //
        //         sum = sum / (double)nbMaturityDates;
        //         prices.add(sum);

        // Pass matrix to event handler
        vector<vector<double> > esdfVals( nbAssets, vector<double>( nbMaturityDates, 0. ) );

        for ( int i = 0; i < nbAssets; ++i )
            for ( int j = 0; j < nbMaturityDates; ++j )
                esdfVals[ i ][ j ] = esdf[ i ] ->getExpSDF( j );

        payoffHandler.add( MatrixPayoffEvent( esdfVals ), weight );

    }
}

void MCTestingCIDMatricesMC::pathGenUpdated( IStateVariableGen::IStateGen* newPathGen )
{
    // create SV
    //        dd = ddGen->getSVDateOfDefault(newPathGen);
    esdf.resize( esdfGen.size() );
    for ( size_t i = 0; i != esdfGen.size(); ++i )
        esdf[ i ] = esdfGen[ i ] ->getSVExpectedSurvivalDiscFactor( newPathGen );
    domDf = domDfGen->getSVDiscFactor( newPathGen ); //  at least one DF has to be requested
    pathW = pathWGen->getSVPathWeight(newPathGen);
    //  domDf = domDfGen->getSVDiscFactor(newPathGen);

    //  if (needFx) {
    //      forDf = forDfGen->getSVDiscFactor(newPathGen);
    //      fx = fxGen->getSpotSV(newPathGen);
    //  }
}

void MCTestingCIDMatricesMC::collectStateVars( IStateVariableCollectorSP svCollector ) const
{
    // collect SV generator
    //        svCollector->append(ddGen.get());
    for ( size_t i = 0; i != esdfGen.size(); ++i )
        svCollector->append( esdfGen[ i ].get() );
    svCollector->append( domDfGen.get() );

    svCollector->append(pathWGen.get());
    //  svCollector->append(domDfGen.get());

    //  if (needFx) {
    //      svCollector->append(forDfGen.get());
    //      svCollector->append(fxGen.get());
    //  }
}

MCTestingCIDMatricesMC::MCTestingCIDMatricesMC( IMultiMarketFactorsSP assets,
        const vector<YieldCurveConstSP>& ycAsset,
        const vector<ICDSParSpreadsConstSP>& cdsAsset,
        const vector<FXAssetConstSP>& _fxAsset,
        const map<string, int>& fxAssetIDMap,
        const YieldCurve* discount,
        const DateTime& today,
        const DateTimeArray& _maturityDates,
        const DateTime& measureDate,
        bool doLog ) :
        MCProductClient( assets.get(), today, discount ),
        ycAsset( ycAsset ),
        cdsAsset( cdsAsset ),
        fxAsset( _fxAsset ),
        fxAssetIDMap( fxAssetIDMap ),
        today( today ),
        maturityDates( _maturityDates ),
        measureDate( measureDate ),
        doLog( doLog ),
        nbMaturityDates( _maturityDates.size() ),
        needFx( ! _fxAsset.empty() )

{
    cerr << "# of ycAsset: " << ycAsset.size() << endl;
    cerr << "# of cdsAsset: " << cdsAsset.size() << endl;
    needFx = ! fxAsset.empty();
    nbAssets = cdsAsset.size();
    // FIXME: For now only ESDFs
    // expected survival prob
    for ( size_t i = 0; i != cdsAsset.size(); ++i ) {
        SVGenExpectedSurvivalDiscFactorSP gen(
            new SVGenExpectedSurvivalDiscFactor(
                measureDate,
                measureDate,
                cdsAsset[ i ],
                maturityDates,
                doLog ) );
        QLIB_VERIFY( gen.get() != NULL, "Failed to create generator" );
        esdfGen.push_back( gen );
    }
    // domestic disc factor
    domDfGen = SVGenDiscFactorSP(
                   new SVGenDiscFactor(
                       today,
                       YieldCurveConstSP( discount ),   // FIXME
                       measureDate )
               );

}


#if 0
MCTestingCIDMatricesMC::MCTestingCIDMatricesMC( const MCTestingCIDMatrices* inst,
        const SimSeriesSP& simSeries,
        InstrumentSettlementSP instSettle )
        :
        MCProductClient(
            IMultiMarketFactors::asMulti( inst->cdsCurveWrapper.getSP() ).get(),
            inst->today,
            inst->yieldCurveWrapper.get(),
            IRefLevelSP( IRefLevel::Util::makeZero( inst->today ) ),   // TODO: fix!
            simSeries,
            IPastValuesSP( IPastValues::Util::makeTrivial( inst->today, 0.0 ) ),   // TODO: fix
            instSettle.get(),
            inst->maturityDates[ 0 ] ),   // TODO: FIX THIS!!
        today( inst->today )
{
    nbMaturityDates = inst->maturityDates.size();
    needFx = !inst->fxWrapper.isEmpty();

    // create SV generators
    // date of default:
    //         ddGen = SVGenDateOfDefaultSP(
    //             new SVGenDateOfDefault(
    //                 inst->today,
    //                 inst->cdsCurveWrapper.getSP())
    //             );

    // expected survival prob
    esdfGen = SVGenExpectedSurvivalDiscFactorSP(
                  new SVGenExpectedSurvivalDiscFactor(
                      inst->measureDate,
                      inst->measureDate,
                      inst->cdsCurveWrapper.getSP(),
                      inst->maturityDates,
                      inst->doLog )   // compute log
              );

    // path weight
    pathWGen.reset(
        new SVGenPathWeight(
            DateTimeArray(1,inst->measureDate))
        );

    // domestic disc factor
    domDfGen = SVGenDiscFactorSP(
                   new SVGenDiscFactor(
                       inst->today,
                       YieldCurveConstSP( inst->yieldCurveWrapper.get() ),
                       instSettle,
                       inst->measureDate )
               );

    if ( needFx ) {
        // create fx SV generator
        fxGen = SVGenSpotSP( new SVGenSpot( 1, inst->measureDate ) );
        forDfGen = SVGenDiscFactorSP(
                       new SVGenDiscFactor(
                           inst->today,
                           YieldCurveConstSP( inst->cdsCurveWrapper.get() ->getYieldCurve() ),
                           instSettle,
                           inst->measureDate )
                   );
    }
}
#endif


DRLIB_END_NAMESPACE
