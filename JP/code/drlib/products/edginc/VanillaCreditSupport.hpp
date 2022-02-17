//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaCreditSupport.hpp
//
//   Description : Credit support for vanilla instrument
//
//   Author      : Ning Shen
//
//   Date        : 11 Sept 2002
//
//
//----------------------------------------------------------------------------

#ifndef VANILLA_CREDIT_SUPPORT_HPP
#define VANILLA_CREDIT_SUPPORT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class CVanilla;

/** Vanilla credit support object  */
class PRODUCTS_DLL VanillaCreditSupport : virtual public CreditSupport
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

    /** interpolate values from grid */
	void InterpValue(const CDoubleArray& spots, CDoubleArray& results);    

    /** store spot-value grid to be used for interpolation
        interp coeff (spline) pre-calculated */
    void SetGrid(const CDoubleArray& spots, const CDoubleArray& values);

    VanillaCreditSupport(CInstrument* inst, CMarketDataSP market);

private:

    IModelSP model;
    
    CVanilla* instrVanOrig;

    smartPtr<CVanilla> instrVan;

	DoubleArray		AtmFwd;

    // spline coeff
    DoubleArrayArray y2;
    // pre-processed values for spline interpolation.
    DoubleArrayArray spotGrid;
    DoubleArrayArray valueGrid;
};

DRLIB_END_NAMESPACE

#endif
