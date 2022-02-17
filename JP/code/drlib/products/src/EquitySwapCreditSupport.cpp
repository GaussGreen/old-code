//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquitySwapCreditSupport.cpp
//
//   Description : Credit support object for Equity Swap
//
//   Author      : Jay Blumenstein
//
//   Date        : 18 Sep 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EquitySwap.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


void EquitySwapCreditSupport::preProcess(const DateTimeArray& dates,
										 const DoubleArray& atmFwd,
										 const DoubleArray& atmVar)
{
    const double tweakSize = 1.01;
    if (!instrESW->isFixedNotional)
    {
        // for each date, need to first compute price array

        // wrap instrument
        IObjectSP inst = CEquitySwapSP::attachToRef(instrESW.get());

        // pre-calc prices with all historical levels populated with spot.
	    CreditSupport::computePriceCache(
		    inst,
		    dates,
		    this,
		    model,
		    priceCache);

        // at each date, now need to compute array of deltas, one for each sample point.

        /////////////// spot delta

        // restore instrument to revert any changes 
        copyInstrument(instrESW, instrESWOrig);
        // wrap instrument for shift
        inst = CEquitySwapSP::attachToRef(instrESW.get());

        CreditSupport::computeDeltaCache(
            inst,
            dates,
            this,
            model,
            tweakSize,
            priceCache,
            spotDelta);

        //  now compute delta caches for each sample list.  

        /////////////// avg in samples

        if (!!instrESW->eqLeg->avIn)
        {
            // restore instrument to revert any changes 
            copyInstrument(instrESW, instrESWOrig);
            // wrap instrument for shift
            inst = CEquitySwapSP::attachToRef(instrESW.get());
            
            computeSampleDeltaCache(
                inst,
                dates,
                 this,
                 model,
                 tweakSize,
                 priceCache,
                 instrESW->eqLeg->avIn->sampleDates,
                 instrESW->eqLeg->avIn->sampleLevels,                 
                 avInDelta);
        
        }

        ///////////// avg out samples

        if (!!instrESW->eqLeg->avOut)
        {
            // restore instrument to revert any changes 
            copyInstrument(instrESW, instrESWOrig);
            // wrap instrument for shift
            inst = CEquitySwapSP::attachToRef(instrESW.get());

            computeSampleDeltaCache(
                 inst,
                 dates,
                 this,
                 model,
                 tweakSize,
                 priceCache,
                 instrESW->eqLeg->avOut->sampleDates,
                 instrESW->eqLeg->avOut->sampleLevels,                 
                 avOutDelta);

        }
            
         //////////// start levels

        // restore instrument to revert any changes 
        copyInstrument(instrESW, instrESWOrig);
        // wrap instrument for shift
        inst = CEquitySwapSP::attachToRef(instrESW.get());

        DateTimeArray eqRefixDates = instrESW->eqLeg->eqRefixDates;
        
      
        DateTimeArray eqStartDates = eqRefixDates;
        eqStartDates.erase(eqStartDates.end() - 1);

        computeSampleDeltaCache(
            inst,
            dates,
            this,
            model,
            tweakSize,
            priceCache,
            eqStartDates,
            instrESW->eqLeg->eqStartLevel,
            eqStartDelta);

        ////////////////// end levels

        // restore instrument to revert any changes 
        copyInstrument(instrESW, instrESWOrig);
        // wrap instrument for shift
        inst = CEquitySwapSP::attachToRef(instrESW.get());

        
        DateTimeArray eqEndDates   = eqRefixDates;
        eqEndDates.erase(eqEndDates.begin());
        
        computeSampleDeltaCache(
            inst,
            dates,
            this,
            model,
            tweakSize,
            priceCache,
            eqEndDates,
            instrESW->eqLeg->eqEndLevel,
            eqEndDelta);
        
   }
    else 
    {
        // fixed Notional case, nothing to do.
    }
}

/** computes double array deltaCache[dates][samples] consisting of the delta 
    with respect to each historical sample */                                  
void EquitySwapCreditSupport::computeSampleDeltaCache(
     IObjectSP& inst,	
     const DateTimeArray& creditDates,
     CreditSupport *creditSupport,
     const IModelSP &model,
     double tweakSize,
     const DoubleArray &priceCache,
     const DateTimeArray& sampleDates,
     DoubleArray& sampleLevels,                 // levels to tweak.  remember to reset.
     vector<vector<double> > &deltaCache)
{
    const string method = "EquitySwapCreditSupport::computeSampleDeltaCache";

    // keep value date and spot
    CInstrument* instPtr = dynamic_cast<CInstrument*>(inst.get());
    if (!instPtr)
    {
        throw ModelException(method, "first parameter must be a smart pointer to an instrument.");
    }
    
    CreditUndSP und = creditSupport->getUnderlier();

    const DateTime valDate = instPtr->getValueDate();
    const double spot = und->getSpot();

    CResults results;

    // now loop through all credit dates, computing sensitivities to each sample, to be stored in the 2D double array
    // deltaCache[credit date][sample date].  Note: for each credit date, the Equity Swap price will be sensitive to
    // sample values on sample dates *only* from and excluding (excluding since there is no variation of spot price at 
    // value date) value date up to and including credit dates, and hence will have non zero sensitivity to to be 
    // recorded in the deltaCache.
    // Other sample dates (i.e. on or before value date, or strictly after credit date) will have deltaCache set to 0
    
    const int numSampleDates = sampleDates.size();
    deltaCache.resize(creditDates.size());

    for (int iCrdDt = 0; iCrdDt < creditDates.size(); iCrdDt++)
    {
        // shift to the credit date
        CreditSupport::shiftValueDate(inst, instPtr->getValueDate(), creditDates[iCrdDt]); 

        deltaCache[iCrdDt].resize(numSampleDates);
        for (int iSmpDt =0; iSmpDt < numSampleDates; iSmpDt++)
        {          
            // if the sample is after value date and on or before credit date, we need to compute delta
            if ( sampleDates[iSmpDt].isGreater(valDate) && 
                !sampleDates[iSmpDt].isGreater(creditDates[iCrdDt]))
            {
                // since we shifted to this date, the Theta method should have populated
                // all samples to spot.  Hence:
                if (sampleLevels.size() <= iSmpDt)
                {
                    throw ModelException(method, 
                        "sampleLevel size = " + Format::toString(sampleLevels.size()) + " but iSmpDt = " + Format::toString(iSmpDt));
                }
                if(!Maths::isZero(sampleLevels[iSmpDt] - spot))
                {
                    throw ModelException(method, "spot = " + Format::toString(spot) + "\n sampleLevels[iSmpDt] = " +
                        Format::toString(sampleLevels[iSmpDt]) + " and note iSmpDt = " + Format::toString(iSmpDt));
                }
                // tweak sample level, compute delta
                sampleLevels[iSmpDt] = spot * (1.0 + tweakSize);
                model->Price(instPtr, smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &results);
                deltaCache[iCrdDt][iSmpDt] = 
                    (results.retrievePrice() - priceCache[iCrdDt])/(tweakSize*spot);
                // now revert back to spot
                sampleLevels[iSmpDt] = spot;
            }
            else
            {
                deltaCache[iCrdDt][iSmpDt] = 0.0;
            }
        }
    }
}

double EquitySwapCreditSupport::changeFromSamples(
            const DateTime& rolledDate,
            const DateTimeArray& sampleDates,
            const DoubleArray& sampleLevels,
            const vector<double>& sampleDelta,
            const double spotRef)
{
    double netChange = 0.0;
    for (int i=0; i < sampleDates.size(); i++)
    {
        // only consider sensitivity to historical dates
        // compare dates to include today's EOD sensitivities.
        if (sampleDates[i].getDate() <= rolledDate.getDate())
        {
            netChange += sampleDelta[i] * (sampleLevels[i] - spotRef);
        }
    }

    return netChange;
}
/** calculate values for a given path */
void EquitySwapCreditSupport::calcPathValues(DoubleArray& results, const DateTimeArray& dates, 
											 const double* spots, double spotRef)
{
    if (!instrESW->isFixedNotional)
    {
        // simply uses a linear function
        for (int i=0; i<dates.size(); i++)
        {   
            // start with price when all samples are at spots
            results[i] = priceCache[i];

            // now add the effect of each sample.
            
            // first set all samples to new values
            if (i > 0)
            {
                rollLinearESWEquity(dates[i-1], spots[i-1],
                    dates[i], spots[i]);
            }

            // effect of spot
            results[i] += spotDelta[i]*(spots[i] - spotRef);


            // avg in samples
            if (!!instrESW->eqLeg->avIn)
            {
                results[i] += changeFromSamples(dates[i], instrESW->eqLeg->avIn->sampleDates, 
                    instrESW->eqLeg->avIn->sampleLevels, avInDelta[i], spotRef);
            }

            // avg out samples
            if (!!instrESW->eqLeg->avOut)
            {
                results[i] += changeFromSamples(dates[i], instrESW->eqLeg->avOut->sampleDates, 
                    instrESW->eqLeg->avOut->sampleLevels, avOutDelta[i], spotRef);
            }

            // need for start & end levels
            DateTimeArray eqRefixDates = instrESW->eqLeg->eqRefixDates;
        
            // start levels
            DateTimeArray eqStartDates = eqRefixDates;
            eqStartDates.erase(eqStartDates.end() - 1);

            results[i] += changeFromSamples(dates[i], eqStartDates, 
                instrESW->eqLeg->eqStartLevel, eqStartDelta[i], spotRef);
        
            // end levels
            DateTimeArray eqEndDates   = eqRefixDates;
            eqEndDates.erase(eqEndDates.begin());
        
            results[i] += changeFromSamples(dates[i], eqEndDates, 
                instrESW->eqLeg->eqEndLevel, eqEndDelta[i], spotRef);
            
        }
    }
    else
    {
        CResults resultObj;
        // keep value date
        DateTime valDate = instrESW->getValueDate();
    
        // wrap instrument
        IObjectSP inst = CEquitySwapSP::attachToRef(instrESW.get());
        CreditUndSP und = getUnderlier();
    
        for (int i=0; i<dates.size(); i++)
        {
        
            // shift value date
            CreditSupport::shiftValueDate(inst, instrESW->getValueDate(), dates[i]);
        
            // linearly interpolate the samples we rolled over
            if (i > 0)
            {
                rollLinearESWEquity(dates[i-1], spots[i-1],
                    dates[i], spots[i]);
            
                // set spot
                und->setSpot(spots[i]);
            }
        
            model->Price(instrESW.get(), smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &resultObj);
            results[i] = resultObj.retrievePrice();
        }
    
        // restore instrument to revert any changes (theta shift changes yield div)
        copyInstrument(instrESW, instrESWOrig);
    }
}


/** Populate any samples in ESWEquity between date1 and date2 using linear interpolation, assuming that the spot level was
value1 and value2 on these dates. */
void EquitySwapCreditSupport::rollLinearESWEquity(const DateTime& date1, double value1, 
                                                  const DateTime& date2, double value2)
{
    if (!!instrESW->eqLeg->avIn)
    {
        CreditSupport::populateSampleArray(instrESW->eqLeg->avIn->sampleDates,
            instrESW->eqLeg->avIn->sampleLevels,
            date1, value1, 
            date2, value2);
    }
    
    if (!!instrESW->eqLeg->avOut)
    {
        CreditSupport::populateSampleArray(instrESW->eqLeg->avOut->sampleDates,
            instrESW->eqLeg->avOut->sampleLevels,
            date1, value1, 
            date2, value2);
    }

    DateTimeArray eqRefixDates = instrESW->eqLeg->eqRefixDates;
    DateTimeArray eqStartDates = eqRefixDates;
    eqStartDates.erase(eqStartDates.end() - 1);

    DateTimeArray eqEndDates   = eqRefixDates;
    eqEndDates.erase(eqEndDates.begin());

    CreditSupport::populateSampleArray(eqStartDates,
            instrESW->eqLeg->eqStartLevel,
            date1, value1, 
            date2, value2);

    CreditSupport::populateSampleArray(eqEndDates,
        instrESW->eqLeg->eqEndLevel,
        date1, value1, 
        date2, value2);
}

EquitySwapCreditSupport::EquitySwapCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    instrESWOrig = dynamic_cast<CEquitySwap*>(inst); //todo right here try to curt
    // copy instrument
	copyInstrument(instrESW, instrESWOrig);
} 

/** return model for this instrument */
IModelSP EquitySwapCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string EquitySwapCreditSupport::getInstCcyCode() const
{
    return instrESW->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime EquitySwapCreditSupport::getInstLastExposureDate() const
{
    return instrESW->endDate(0);
}

/** return asset */
CreditUndSP EquitySwapCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( instrESW->asset.get() ) );
}

DRLIB_END_NAMESPACE
