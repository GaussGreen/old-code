//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaCreditSupport.cpp
//
//   Description : Credit support for vanilla instrument
//
//   Author      : Ning Shen
//
//   Date        : 11 Sept 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/Vanilla.hpp"

DRLIB_BEGIN_NAMESPACE

VanillaCreditSupport::VanillaCreditSupport(CInstrument* inst, CMarketDataSP market){
    // keep original
    instrVanOrig = dynamic_cast<CVanilla*>(inst);
    // model created in the constructor
    if (instrVanOrig->canExerciseEarly)
        model = IModelSP(Tree1fLN::make("VolPreferred", false));
    else
        model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // copy instrument
    copyInstrument(instrVan, instrVanOrig);
}

/** return model for this instrument */
IModelSP VanillaCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string VanillaCreditSupport::getInstCcyCode() const
{
    return instrVan->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime VanillaCreditSupport::getInstLastExposureDate() const
{
    return instrVan->endDate(0);
}

/** return asset */
CreditUndSP VanillaCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( instrVan->asset.get() ) );
}

/** preprocess instrument for a given set of path dates */
void VanillaCreditSupport::preProcess(const DateTimeArray& dates, 
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
    DateTime valDate = instrVan->getValueDate();

    CResults results;

    // wrap instrument
    IObjectSP inst = CVanillaSP::attachToRef(instrVan.get());
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
        CreditSupport::shiftValueDate(inst, instrVan->getValueDate(), dates[i]);

        for (j=0; j<NUM_VALUE_GRIDS; j++)
        {
            und->setSpot(spotGrid[i][j]);
            model->Price(instrVan.get(), smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &results);
            valueGrid[i][j] = results.retrievePrice();
        }
        // prepare cubic spline coefficients
        spline(&spotGrid[i][0]-1, &valueGrid[i][0]-1, NUM_VALUE_GRIDS, y1, yn, &y2[i][0]-1);
	}

	// restore instrument to revert any changes (theta shift changes yield div)
	copyInstrument(instrVan, instrVanOrig);
	// store atm fwd for fwd starting
	AtmFwd = atmFwd;
}

/** calculate values for a given path */
void VanillaCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double* spots,
                                                  double spotRef)
{
	double scaleFactor = 1.0;
	bool alreadyComputed = !instrVan->fwdStarting;

    // uses cubic spline 
    for (int i=0; i<dates.size(); i++)
    {
		// need to compute scale factor for forward starting inclusive of fwd start date 
		if (!alreadyComputed && dates[i] >= instrVan->startDate) 
		{
			scaleFactor = CreditSupport::computeScaleFactor(instrVan->startDate, 
															dates, 
															spots,
															instrVan->getValueDate(),
															AtmFwd[i]);
			alreadyComputed = true;
		}

        // spline interp
        if (Maths::isZero(spotGrid[i][0] - spotGrid[i][spotGrid[i].size()-1]))
            results[i] = scaleFactor*valueGrid[i][0]; // no vol or value date
        else
        {
			// adjusted spot floored at 0 and scaled if fwdstart. then price
			double spot = Maths::max( scaleFactor *  spots[i], MIN_EQUITY_PRICE ); 
            splint(&spotGrid[i][0]-1, &valueGrid[i][0]-1, &y2[i][0]-1, spotGrid[i].size(), spot, &results[i]);
        }
        // floored to 0
        if(results[i] < 0.0)
            results[i] = 0.0;
    }
}

DRLIB_END_NAMESPACE

