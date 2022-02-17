//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContractCreditSupport.cpp
//
//   Description : Credit support object for ForwardContract
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ForwardContract.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

ForwardContractCreditSupport::ForwardContractCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    fwdOriginal = dynamic_cast<CForwardContract*>(inst);
    // copy instrument
    copyInstrument(fwdContract, fwdOriginal);

	if ( fwdContract->fwdStarting && fwdContract->oneContract)
    {
        throw ModelException("ForwardContractCreditSupport::ForwardContractCreditSupport", 
						"Can't be forward starting and one contract");
    }
} 

/** return model for this instrument */
IModelSP ForwardContractCreditSupport::getModel()
{
    return model;
}

/** return instrument ccy ISO code */
string ForwardContractCreditSupport::getInstCcyCode() const
{
    return fwdContract->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime ForwardContractCreditSupport::getInstLastExposureDate() const
{
    return fwdContract->endDate(0);
}

/** return asset */
CreditUndSP ForwardContractCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( fwdContract->asset.get() ) );
}

/** preprocess instrument for a given set of path dates */
void ForwardContractCreditSupport::preProcess(const DateTimeArray& dates,
											  const DoubleArray&,
											  const DoubleArray&)
{
    // wrap instrument for shift
    IObjectSP inst = ForwardContractSP::attachToRef(fwdContract.get());

	// pre-calc forwards for linear interpolation
	CreditSupport::computePriceCache(
		inst,
		dates,
		this,
		model,
		FwdCache);
	
	// restore instrument to revert any changes (theta shift changes yield div)
     copyInstrument(fwdContract, fwdOriginal);
    // wrap instrument for shift
    inst = ForwardContractSP::attachToRef(fwdContract.get());

	double tweakSize = 0.01;

	CreditSupport::computeDeltaCache(
		inst,
		dates,
		this,
		model,
		tweakSize,
		FwdCache,
		FwdDeltaCache);
    
    // restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(fwdContract, fwdOriginal);
}

/** calculate values for a given path */
void ForwardContractCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double* spots,
                                                  double spotRef)
{
	double scaleFactor=1.0;
	bool alreadyComputed = !fwdContract->fwdStarting;

    // simply uses a linear function
    for (int i=0; i<dates.size(); i++)
    {
		// need to compute scale factor for forward starting inclusive of fwd start date 
		if (!alreadyComputed && dates[i] >= fwdContract->startDate) 
		{
			scaleFactor = CreditSupport::computeScaleFactor(fwdContract->startDate, 
															dates, 
															spots, 
															fwdContract->getValueDate(),
															spotRef);
			alreadyComputed = true;
		}

		// adjusted spot floored at 0 and scaled if fwdstart. then price
		double spot = Maths::max( scaleFactor *  spots[i], MIN_EQUITY_PRICE ); 
        results[i] = FwdCache[i] + FwdDeltaCache[i]*(spot  - spotRef);
    }
}

DRLIB_END_NAMESPACE

