//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LiborStreamCreditSupport.cpp
//
//   Description : Credit support object for LiborStream
//
//   Author      : Ning Shen
//
//   Date        : 19 Sept 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LiborStream.hpp"

DRLIB_BEGIN_NAMESPACE

LiborStreamCreditSupport::LiborStreamCreditSupport(CInstrument* inst, CMarketDataSP market){
    // model created in the constructor
    model = IModelSP(new ClosedForm());

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    liborOriginal = dynamic_cast<LiborStream*>(inst);
    // copy instrument
    copyInstrument(libor, liborOriginal);
} 

/** return model for this instrument */
IModelSP LiborStreamCreditSupport::getModel()
{
    return model;
}

/** return instrument discount curve */
string LiborStreamCreditSupport::getInstCcyCode() const
{
    return libor->discount->getCcy();
}

/** return instrument's last exposure date */
DateTime LiborStreamCreditSupport::getInstLastExposureDate() const
{
    return libor->endDate(0);
}

/** return asset */
CreditUndSP LiborStreamCreditSupport::getUnderlier() const
{
    return CreditUndSP();
}

/** preprocess instrument for a given set of path dates */
void LiborStreamCreditSupport::preProcess(const DateTimeArray& dates, const DoubleArray&, const DoubleArray&)
{
    // wrap instrument for shift
    IObjectSP inst = LiborStreamSP::attachToRef(libor.get());

	// pre-calc forwards for linear interpolation
	CreditSupport::computePriceCache(
		inst,
		dates,
		this,
		model,
		cfValues);
	
	// restore instrument to revert any changes (theta shift changes yield div)
    copyInstrument(libor, liborOriginal);
}

/** calculate values for a given path */
void LiborStreamCreditSupport::calcPathValues(DoubleArray& results, 
                                                  const DateTimeArray& dates, 
                                                  const double*,
                                                  double )
{
    // loop through dates */
    for (int i = 0; i<dates.size(); i++)
        results[i] = cfValues[i];
}

DRLIB_END_NAMESPACE

