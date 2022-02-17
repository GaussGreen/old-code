//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AverageCreditSupport.hpp
//
//   Description : Credit support for Average instrument
//
//   Author      : Ning Shen
//
//   Date        : 11 Sept 2002
//
//
//----------------------------------------------------------------------------

#ifndef AVERAGE_CREDIT_SUPPORT_HPP
#define AVERAGE_CREDIT_SUPPORT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class Average;

/** Average credit support object  */
class PRODUCTS_DLL AverageCreditSupport : virtual public CreditSupport
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

    AverageCreditSupport(CInstrument* inst, CMarketDataSP market);

private:

    IModelSP model;
    
    Average* instrAveOrig;

    smartPtr<Average> instrAve;

    // pre-processed values for spline interpolation.
    DoubleArrayArray AveInSoFar;
    DoubleArrayArray AveOutSoFar;
    DoubleArrayArray SpotGrid;
	DoubleArray		 AtmFwd; // only used for fwd starting avg spot

    // [sample dates][AviIn grid][AveOut grid][spot grid]
    vector<vector<vector<vector<double> > > > ValueGrid;
    vector<vector<vector<vector<double> > > > Y2;

    // initialise sizes and set up grids
    void setupGrid(	const DateTimeArray& dates, 
					const DoubleArray& atmFwd, 
					const DoubleArray& atmVar);
};

DRLIB_END_NAMESPACE

#endif
