//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashFlowStreamCreditSupport.cpp
//
//   Description : Credit support object for CashFlowStream
//
//   Author      : Ning Shen
//
//   Date        : 19 Sept 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashFlowStream.hpp"

DRLIB_BEGIN_NAMESPACE

CashFlowStreamCreditSupport::CashFlowStreamCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new ClosedForm());

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    cfOriginal = dynamic_cast<CashFlowStream*>(inst);
    // copy instrument
    copyInstrument(cfStream, cfOriginal);
} 

/** return model for this instrument */
IModelSP CashFlowStreamCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string CashFlowStreamCreditSupport::getInstCcyCode() const
{
    return cfStream->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime CashFlowStreamCreditSupport::getInstLastExposureDate() const
{
    return cfStream->endDate(0);
}

/** return asset */
CreditUndSP CashFlowStreamCreditSupport::getUnderlier() const
{
    return CreditUndSP();
}

/** preprocess instrument for a given set of path dates */
void CashFlowStreamCreditSupport::preProcess(const DateTimeArray& dates, const DoubleArray&, const DoubleArray&)
{
    // wrap instrument for shift
    IObjectSP inst = CashFlowStreamSP::attachToRef(cfStream.get());

	// pre-calc forwards for linear interpolation
	CreditSupport::computePriceCache(
		inst,
		dates,
		this,
		model,
		cfValues);
	
	// restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(cfStream, cfOriginal);
}

/** calculate values for a given path */
void CashFlowStreamCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double*,
                                                  double )
{
    // loop through dates */
    for (int i = 0; i<dates.size(); i++)
        results[i] = cfValues[i];
}

DRLIB_END_NAMESPACE

