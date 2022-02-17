//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetValueCreditSupport.cpp
//
//   Description : Credit support object for AssetValue
//
//   Author      : Ning Shen
//
//   Date        : 19 Sept 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetValue.hpp"
#include "edginc/AssetUtil.hpp"

DRLIB_BEGIN_NAMESPACE

AssetValueCreditSupport::AssetValueCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new CClosedFormLN("VolPreferred"));

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    instrAVOriginal = dynamic_cast<AssetValue*>(inst);
    // copy instrument
    copyInstrument(instrAV, instrAVOriginal);
} 

/** return model for this instrument */
IModelSP AssetValueCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string AssetValueCreditSupport::getInstCcyCode() const
{
    return AssetUtil::assetCcy(instrAV->asset.get());
}

/** return instrument's last exposure date */
DateTime AssetValueCreditSupport::getInstLastExposureDate() const
{
    // asset does not have a maturity, so we should use the last requested exposure date.
    return instrAV->valueDate.rollDate(12000); // to do : need to be able to use last request date
}

/** return asset */
CreditUndSP AssetValueCreditSupport::getUnderlier() const
{
    return CreditUndSP( new AssetUnd( instrAV->asset.get() ) );
}

/** preprocess instrument for a given set of path dates */
void AssetValueCreditSupport::preProcess(const DateTimeArray& dates, const DoubleArray&, const DoubleArray&)
{
    // wrap instrument for shift
    IObjectSP inst = AssetValueSP::attachToRef(instrAV.get());

	// pre-calc forwards for linear interpolation
	CreditSupport::computePriceCache(
		inst,
		dates,
		this,
		model,
		FwdCache);
	
	// restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(instrAV, instrAVOriginal);
    // wrap instrument for shift
    inst = AssetValueSP::attachToRef(instrAV.get());

	double tweakSize = 1.01;

	CreditSupport::computeDeltaCache(
		inst,
		dates,
		this,
		model,
		tweakSize,
		FwdCache,
		FwdDeltaCache);
    
    // restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(instrAV, instrAVOriginal);
}

/** calculate values for a given path */
void AssetValueCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double* spots,
                                                  double spotRef)
{
    // simply uses a linear function
    for (int i=0; i<dates.size(); i++)
    {
        results[i] = FwdCache[i] + FwdDeltaCache[i]*(spots[i] - spotRef);
        if (results[i] < 0.0)
            results[i] = 0.0;
    }
}

DRLIB_END_NAMESPACE

