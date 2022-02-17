//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LiborStreamCreditSupport.hpp
//
//   Description : Credit support object for LiborStream
//
//   Author      : Ning Shen
//
//   Date        : 19 Sept 2002
//
//
//----------------------------------------------------------------------------

#ifndef LIBORSTREAM_CREDIT_HPP
#define LIBORSTREAM_CREDIT_HPP

#include "edginc/Asset.hpp"
#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class LiborStream;

/** LiborStream credit support object  */
class PRODUCTS_DLL LiborStreamCreditSupport : virtual public CreditSupport
{
public:

    /** preprocess instrument for a given set of path dates */
    virtual void preProcess(const DateTimeArray& dates, 
							const DoubleArray& atmFwd, 
							const DoubleArray& atmVar);

    /** return asset */
    virtual CreditUndSP getUnderlier() const;

    /** calculate values for a given path */
    virtual void calcPathValues(DoubleArray& results, const DateTimeArray& dates, 
                    const double* spots, double spotRef);

    /** return model for this instrument */
    virtual IModelSP getModel();

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const;

	/** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const;

    LiborStreamCreditSupport(CInstrument*, CMarketDataSP market);

private:
    IModelSP model;
    
    LiborStream* liborOriginal;

    smartPtr<LiborStream> libor;

    // pre-processed values for spline interpolation.
    DoubleArray cfValues;
};

DRLIB_END_NAMESPACE

#endif
