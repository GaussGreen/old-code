//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AverageCreditSupport.cpp
//
//   Description : Credit support for Average instrument
//
//   Author      : Ning Shen
//
//   Date        : 12 Sept 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Average.hpp"

DRLIB_BEGIN_NAMESPACE

// grid size
static const int NUM_VALUE_GRIDS = 11;

AverageCreditSupport::AverageCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    instrAveOrig = dynamic_cast<Average*>(inst);
    // copy instrument
    copyInstrument(instrAve, instrAveOrig);
}

/** return model for this instrument */
IModelSP AverageCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string AverageCreditSupport::getInstCcyCode() const
{
    return instrAve->discount->getCcy();;
}

/** return instrument's last exposure date */
DateTime AverageCreditSupport::getInstLastExposureDate() const
{
    return instrAve->endDate(0);
}

/** return asset */
CreditUndSP AverageCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( instrAve->asset.get() ) );
}

// initialise sizes
void AverageCreditSupport::setupGrid(const DateTimeArray& dates, 
									 const DoubleArray& atmFwd, 
									 const DoubleArray& atmVar)
{
    // is this scaling good ?
    const double factor[NUM_VALUE_GRIDS] = 
                            {-3.5, -2.0, -1.1, -0.5, -0.15, 0.0,
                            0.15, 0.5, 1.1, 2.0, 3.5};
    int i, j;
    int numDate = dates.size();
    SampleListSP aveIn = instrAve->getAveIn();

	// set up the spot grid first
	setupSpotGrid(dates, atmFwd, atmVar, NUM_VALUE_GRIDS, factor, SpotGrid);

    // resize all arrays
    AveInSoFar.resize(numDate);
    AveOutSoFar.resize(numDate);
    ValueGrid.resize(numDate);
    Y2.resize(numDate);

    for (i=0; i<numDate; i++)
    {
        if (!!aveIn)
            AveInSoFar[i].resize(NUM_VALUE_GRIDS);
        else
            AveInSoFar[i].resize(1);

        AveOutSoFar[i].resize(NUM_VALUE_GRIDS);
        ValueGrid[i].resize(NUM_VALUE_GRIDS);
        Y2[i].resize(NUM_VALUE_GRIDS);
        for (int j=0; j<NUM_VALUE_GRIDS; j++)
        {
            // set value grid size aveOut and spot dimension
            ValueGrid[i][j].resize(NUM_VALUE_GRIDS);
            Y2[i][j].resize(NUM_VALUE_GRIDS);
            for (int k=0; k<NUM_VALUE_GRIDS; k++)
            {
                ValueGrid[i][j][k].resize(NUM_VALUE_GRIDS);
                Y2[i][j][k].resize(NUM_VALUE_GRIDS);
            }
        }
    }

    // set up grids
    for (i=0; i<numDate; i++)
    {
        if (!!aveIn)
        {
            // do ave in grids
            if (aveIn->getFirstDate() > dates[i])
            {
                AveInSoFar[i].resize(1); // only one point needed at 0
                AveInSoFar[i][0] = 0.0;
            }
            else
            {
                for (j=0; j<NUM_VALUE_GRIDS; j++)
                { //set ave in grid
                    AveInSoFar[i][j] = SpotGrid[i][j];
                }
            }
        }
        // do ave out grids
        if (instrAve->avgOut->getFirstDate() > dates[i])
        {
            AveOutSoFar[i].resize(1); // only one point needed at 0
            AveOutSoFar[i][0] = 0.0;
        }
        else
        {
            for (j=0; j<NUM_VALUE_GRIDS; j++)
            { //set ave in grid
                AveOutSoFar[i][j] = SpotGrid[i][j];
            }
        }
    }
}

// helper function
static inline void fillArray(DoubleArray& result, double v)
{
    for (int i=0; i<result.size(); i++)
        result[i] = v;
}

/** preprocess instrument for a given set of path dates */
void AverageCreditSupport::preProcess(const DateTimeArray& dates, 
									  const DoubleArray& atmFwd, 
									  const DoubleArray& atmVar)
{
    // for spline
  	const double y1 = 2e30; // default to natural spline
	const double yn = 2e30;

    int i, j, k, m;

    // keep value date
    DateTime valDate = instrAve->getValueDate();

    CResults results;

    // wrap instrument
    IObjectSP inst = AverageSP::attachToRef(instrAve.get());
    CreditUndSP und = getUnderlier();

    // set up grids
    setupGrid(dates, atmFwd, atmVar); 

    SampleListSP aveIn = instrAve->getAveIn();
    DoubleArray* aveInFixing = 0;

    for (i=0; i<dates.size(); i++)
    {
        if (!!aveIn){
            aveInFixing = &(const_cast<DoubleArray&>(aveIn->getValues()));
        }
        DoubleArray& aveOutFixing = const_cast<DoubleArray&>(instrAve->avgOut->getValues());
        // shift value date. restore atmFwd first
		und->setSpot( atmFwd[i] );
        CreditSupport::shiftValueDate(inst, instrAve->getValueDate(), dates[i]);

        // set sample to ave in so far grid
        for (j=0; j<AveInSoFar[i].size(); j++)
        {
            if (aveInFixing)
                fillArray(*aveInFixing, AveInSoFar[i][j]);
            // set sample to ave out so far grid
            for (k=0; k<AveOutSoFar[i].size(); k++)
            {
                fillArray(aveOutFixing, AveOutSoFar[i][k]);
                // set spot grid
                for (m=0; m<NUM_VALUE_GRIDS; m++)
                {
                    und->setSpot(SpotGrid[i][m]);
                    model->Price(instrAve.get(), smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &results);
                    ValueGrid[i][j][k][m] = results.retrievePrice();
                }
                // prepare cubic spline coefficients
                spline(&SpotGrid[i][0]-1, &ValueGrid[i][j][k][0]-1, 
                                NUM_VALUE_GRIDS, y1, yn, &Y2[i][j][k][0]-1);
            }
        }
    }

	// restore instrument to revert any changes (theta shift changes yield div)
	copyInstrument(instrAve, instrAveOrig);

	// store the atm fwd
	AtmFwd = atmFwd;
}

/** calculate values for a given path */
void AverageCreditSupport::calcPathValues(DoubleArray& results, 
                                          const DateTimeArray& dates, 
                                          const double* spots,
                                          double spotRef)
{
  	const double y1 = 2e30; // default to natural spline
	const double yn = 2e30;
	double scaleFactor = 1.0;
	DateTime startDate;
	bool alreadyComputed = !instrAve->getFwdStartDate( startDate );

    int i, j, k;
    double inSoFar, outSoFar;
//    double error;

    double ya[2][NUM_VALUE_GRIDS]; // values
    double yb[2][NUM_VALUE_GRIDS]; // 2nd deriv

    SampleListSP aveIn = instrAve->getAveIn();
    // uses cubic spline 
    for (i=0; i<dates.size(); i++)
    {
		// need to compute scale factor for forward starting inclusive of fwd start date 
		if (!alreadyComputed && dates[i] >= startDate) 
		{
			scaleFactor = CreditSupport::computeScaleFactor(startDate, 
															dates, 
															spots,
															instrAve->getValueDate(),
															AtmFwd[i]);
			alreadyComputed = true;
		}

        if (Maths::isZero(SpotGrid[i][0] - SpotGrid[i][SpotGrid[i].size()-1]))
            results[i] = ValueGrid[i][0][0][0]; // on value date
        else
        {
		    // roll over sample dates first, linearly interpolate sample fixing values
            if (i>0)
            {
                if (!!aveIn)
                {
                    aveIn->rollLinear(dates[i-1], spots[i-1],
                                      dates[i], spots[i]);
                }
                instrAve->avgOut->rollLinear(dates[i-1], spots[i-1],
                                              dates[i], spots[i]);
            }

			// adjusted spot floored at 0 and scaled if fwdstart. then price
			double spot = Maths::max( scaleFactor *  spots[i], MIN_EQUITY_PRICE ); 

			// calc average so far
            if (AveInSoFar[i].size() == 1)
                inSoFar = 0.0;
            else
                inSoFar = aveIn->averageToDate(dates[i]) * scaleFactor;

            if (AveOutSoFar[i].size() == 1)
                outSoFar = 0.0;
            else
                outSoFar = instrAve->avgOut->averageToDate(dates[i]) * scaleFactor;

/*
            // splint on spot but poly interpolate on average so far
            for (j=0; j<AveInSoFar[i].size(); j++)
            {
                for (k=0; k<AveOutSoFar[i].size(); k++)
                {
                    splint(SpotGrid[i].begin()-1, ValueGrid[i][j][k].begin()-1, Y2[i][j][k].begin()-1, 
                                        SpotGrid[i].size(), spots[i], &ya[0][k]);
                }
                if (AveOutSoFar[i].size() == 1)
                    ya[1][j] = ya[0][0];
                else
                {
                    // polynomial interp
                    polint(AveOutSoFar[i].begin()-1, ya[0]-1, AveOutSoFar[i].size(), outSoFar, &ya[1][j], &error);
                }
            }
            if (AveInSoFar[i].size() == 1)
                results[i] = ya[1][0];
            else
            {
                // polynomila interp
                polint(AveInSoFar[i].begin()-1, ya[1]-1, AveInSoFar[i].size(), inSoFar, &results[i], &error);
            }
*/

            // spline interpolate
            for (j=0; j<AveInSoFar[i].size(); j++)
            {
                for (k=0; k<AveOutSoFar[i].size(); k++)
                {
                    splint(&SpotGrid[i][0]-1, &ValueGrid[i][j][k][0]-1, &Y2[i][j][k][0]-1, 
                                        SpotGrid[i].size(), spot, &ya[0][k]);
                }
                if (AveOutSoFar[i].size() == 1)
                    ya[1][j] = ya[0][0];
                else
                {
                    spline(&*AveOutSoFar[i].begin()-1, ya[0]-1, AveOutSoFar[i].size(), y1, yn, yb[0]-1);
                    splint(&*AveOutSoFar[i].begin()-1, ya[0]-1, yb[0]-1, AveOutSoFar[i].size(), outSoFar, &ya[1][j]);
                }
            }
            if (AveInSoFar[i].size() == 1)
                results[i] = ya[1][0];
            else
            {
                spline(&*AveInSoFar[i].begin()-1, ya[1]-1, AveInSoFar[i].size(), y1, yn, yb[1]-1);
                splint(&*AveInSoFar[i].begin()-1, ya[1]-1, yb[1]-1, AveInSoFar[i].size(), inSoFar, &results[i]);
            }

            if (results[i] < 0.0)
                results[i] = 0.0;
        }
    }
}

DRLIB_END_NAMESPACE

