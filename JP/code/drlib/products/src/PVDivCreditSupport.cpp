//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PVDivCreditSupport.cpp
//
//   Description : Credit support object for PVDiv
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PVDiv.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE

PVDivCreditSupport::PVDivCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    pvDivOriginal = dynamic_cast<CPVDiv*>(inst);
    // copy instrument
	copyInstrument(pvDiv, pvDivOriginal);
} 


/** return model for this instrument */
IModelSP PVDivCreditSupport::getModel()
{
    return model;
}

/** return instrument ccy ISO code */
string PVDivCreditSupport::getInstCcyCode() const
{
    return pvDiv->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime PVDivCreditSupport::getInstLastExposureDate() const
{
    return pvDiv->endDate(0);
}

/** return asset */
CreditUndSP PVDivCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( pvDiv->equityAsset.get() ) );
}

/** calculate values for a given path */
void PVDivCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double* spots,
                                                  double spotRef)
{
	/** if nominal case or classical forward starting, use linear func */
	if (!pvDiv->isFixedNotional || pvDiv->isLegacyFwdStarting())
	{
		double scaleFactor = 1.0;
		bool alreadyComputed = !pvDiv->isLegacyFwdStarting();  // so if it isn't fwd starting, we've indeed already computed (on the previous line).

		DateTime startDate = DateTime(0,0);
        if (pvDiv->isLegacyFwdStarting())
		{
			startDate = pvDiv->averageSched->getDates()[0];
		}

		// simply uses a linear function
        for (int i=0; i<dates.size(); i++)
        {
            // need to compute scale factor for forward starting
            if (!alreadyComputed
                && dates[i] >= startDate)
            {
                scaleFactor = CreditSupport::computeScaleFactor(startDate, 
                    dates, 
                    spots, 
                    pvDiv->getValueDate(),
                    spotRef);
                alreadyComputed = true;
            }

            if (pvDiv->isLegacyFwdStarting() && !alreadyComputed)
            {
                results[i] = priceCache[i] * spotRef / Maths::max(spots[i], MIN_EQUITY_PRICE);
            }
            else
            {
                // note that both delta and price are computed with 1/spotRef scaled for dates >= fwd starting date, 
                // hence we multiply by scaleFactor outside to adjust.
    			results[i] = priceCache[i]; // + deltaCache[i]*(spots [i] - spotRef);  add this when yield divs are properly handled.
            }
			results[i] *= scaleFactor;
		}
	}
	else
	{
        throw ModelException("PVDivCreditSupport::calcPathValues", "Fixed notional only supported for InvestOnWeightedAverage with one average in date");
    }

#if 0 
        // if this is necessary...
        vector<DateTimeArray> fixDatesArr;
		DateTimeArray allFixingDates;
		getFixingDates(allFixingDates);

		// wrap instrument for shift
		IObjectSP inst = PVDivSP::attachToRef(pvDiv.get());
		CreditUndSP und = getUnderlier();		
		
		for (int i=0; i<dates.size(); i++)
		{
			CreditSupport::shiftValueDate(inst, pvDiv->getValueDate(), dates[i]);
       
            // to do
			results[i] = 2;
		}
	}

#endif

}

/** preprocess instrument for a given set of path dates */
void PVDivCreditSupport::preProcess(const DateTimeArray& dates, 
									const DoubleArray& atmFwd, 
									const DoubleArray& atmVar)
{
    // wrap instrument for shift
    IObjectSP inst = PVDivSP::attachToRef(pvDiv.get());

	// pre-calc pv div price
	CreditSupport::computePriceCache(
		inst,
		dates,
		this,
		model,
		priceCache);
	
	// restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(pvDiv, pvDivOriginal);

    // until yield dividends are properly handled (rolling over yield dividends currently populates with spot, not the interpolated value 
    // based on simulated spot which bracket the yield div date.
#if 0
	// wrap instrument for shift
    inst = PVDivSP::attachToRef(pvDiv.get());

	double tweakSize = 0.01;

	CreditSupport::computeDeltaCache(
		inst,
		dates,
		this,
		model,
		tweakSize,
		priceCache,
		deltaCache);
    
    // restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(pvDiv, pvDivOriginal);
#endif
}

void PVDivCreditSupport::getFixingDates(DateTimeArray& allFixingDates)
{
    if (pvDiv->isFixedNotional)
    {
        // if type 1 avg need the dates in the contract schedule
        if (pvDiv->averagingType == "InvestOnContractDates")
        {
            allFixingDates = pvDiv->contractSched->getDates();
        }

        // if type 2 avg need all dates in average schedule
        if (pvDiv->averagingType == "InvestOnWeightedAverage")
        {
            allFixingDates = pvDiv->averageSched->getDates(); 
        }
    }


}
DRLIB_END_NAMESPACE

