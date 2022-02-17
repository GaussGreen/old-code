//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CDSCreditSupport.cpp
//
//   Description : Credit support for CDS instrument
//
//   Author      : Dmytro Zhuravytsky
//
//   Date        : 20 May 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CredDefSwap.hpp"
#include "edginc/CDSCreditSupport.hpp"
#include "edginc/ICreditEventOverrideName.hpp"

DRLIB_BEGIN_NAMESPACE


ParSpreadsUnd::ParSpreadsUnd( ICDSParSpreads * spreads, const CVolBase * vol ) :
    m_spreads( spreads ),
    m_vol( vol ),
    m_spot( 1. )
{
}
string ParSpreadsUnd::getName() const
{
    return m_spreads->getName();
}
string ParSpreadsUnd::getTrueName() const
{
    return m_spreads->getName();
}
CVolProcessed * ParSpreadsUnd::getProcessedVol(
    const CVolRequest * volRequest ) const
{
    return m_vol->getProcessedVol( volRequest, 0 );
}
void ParSpreadsUnd::fwdValue(
    const DateTimeArray & dateList,
    CDoubleArray        & result ) const
{
    for( int i = 0; i < dateList.size(); ++i )
        result[ i ] = 1.;
}
double ParSpreadsUnd::getSpot() const
{
    return m_spot;
}
void ParSpreadsUnd::setSpot(double spot)
{
    ParSpreadPropShift newSpot(spot/m_spot - 1.);
    m_spot = spot;
 
    OutputNameConstSP name(new OutputName(getTrueName()));
    if (!newSpot.findAndShift(IObjectSP::attachToRef(m_spreads), name)){
        throw ModelException( "ParSpreadsUnd::setSpot", "Could not change par spreads of '" + getTrueName() + "', par spreads may not be tweakable proportionally" );
    }
}

CDSCreditSupport::CDSCreditSupport(CInstrument* inst, CMarketDataSP market)
{
    // keep original
    instrCDSOrig = dynamic_cast<CredDefSwap*>(inst);
    // model created in the constructor
    model = IModelSP(new ClosedFormCDSPS());
    // call inst GetMarket to initiate market data selection.
    inst->GetMarket(model.get(), market);
    // call model to see if wants any extra data
    model->getMarket(market.get(), IInstrumentCollection::singleton(inst));
    // copy instrument
    copyInstrument(instrCDS, instrCDSOrig);

    //get vol data
    vol.setName( instrCDS->cdsParSpreads->getName() );
    vol.getData( model.get(), market.get() );
}

/** return model for this instrument */
IModelSP CDSCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string CDSCreditSupport::getInstCcyCode() const
{
    return instrCDS->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime CDSCreditSupport::getInstLastExposureDate() const
{
    return instrCDS->endDate(0);
}

/** return par spreads */
CreditUndSP CDSCreditSupport::getUnderlier() const
{
    return CreditUndSP( new ParSpreadsUnd( instrCDS->cdsParSpreads.get(), vol.get() ) );
}

/** preprocess instrument for a given set of path dates */
void CDSCreditSupport::preProcess(const DateTimeArray& dates, 
                                  const DoubleArray& atmFwd, 
                                  const DoubleArray& atmVar)
{
    const int NUM_VALUE_GRIDS = 21; // 10 on each side of strike
    // is this scaling good ?
    const double factor[NUM_VALUE_GRIDS] = 
    {-3.5, -2.5, -2.0, -1.5, -1.1, -0.8, -0.5, -0.3, -0.15, -0.05,
    0.0,
    0.05, 0.15, 0.3, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.5};

    // for spline
    const double y1 = 2e30; // default to natural spline
    const double yn = 2e30;

    int i, j;

    // keep value date
    DateTime valDate = instrCDS->getValueDate();

    CResults results;

    // wrap instrument
    IObjectSP inst = CredDefSwapSP::attachToRef(instrCDS.get());
    CreditUndSP und = getUnderlier();

    // set up spot grid
    setupSpotGrid(dates, atmFwd, atmVar, NUM_VALUE_GRIDS, factor, spotGrid);

    // init grid scaling
    valueGrid.resize(dates.size());
    y2.resize(dates.size());
    for (i=0; i<dates.size(); i++)
    {
        valueGrid[i].resize(NUM_VALUE_GRIDS);
        y2[i].resize(NUM_VALUE_GRIDS);
 
        // shift value date. need to restore atmFwd. useful for setting spot-at-forward-start
        und->setSpot( atmFwd[i] ); 
        CreditSupport::shiftValueDate(inst, instrCDS->getValueDate(), dates[i]);

        for (j=0; j<NUM_VALUE_GRIDS; j++)
        {
            und->setSpot(spotGrid[i][j]);
            model->Price(instrCDS.get(), smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &results);
            valueGrid[i][j] = results.retrievePrice();
        }
        // prepare cubic spline coefficients
        spline(&spotGrid[i][0]-1, &valueGrid[i][0]-1, NUM_VALUE_GRIDS, y1, yn, &y2[i][0]-1);
    }

    // restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(instrCDS, instrCDSOrig);
    // store atm fwd for fwd starting
    AtmFwd = atmFwd;
}

/** calculate values for a given path */
void CDSCreditSupport::calcPathValues(DoubleArray& results, 
                                      const DateTimeArray& dates, 
                                      const double* spots,
                                      double spotRef)
{
    // uses cubic spline 
    for (int i=0; i<dates.size(); i++)
    {
        // spline interp
        if (Maths::isZero(spotGrid[i][0] - spotGrid[i][spotGrid[i].size()-1]))
            results[i] = valueGrid[i][0]; // no vol or value date
        else
        {
            // adjusted spot floored at 0 and scaled if fwdstart. then price
            double spot = Maths::max( spots[i], MIN_EQUITY_PRICE ); 
            splint(&spotGrid[i][0]-1, &valueGrid[i][0]-1, &y2[i][0]-1, spotGrid[i].size(), spot, &results[i]);
        }
    }
}

DRLIB_END_NAMESPACE
