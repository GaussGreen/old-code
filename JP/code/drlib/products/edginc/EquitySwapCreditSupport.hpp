//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquitySwapCreditSupport.hpp
//
//   Description : Credit support object for Equity Swap
//
//   Author      : Jay Blumenstein
//
//   Date        : 19 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef EQUITY_SWAP_CREDIT_SUPPORT_HPP
#define EQUITY_SWAP_CREDIT_SUPPORT_HPP

#include "edginc/config.hpp"
#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class CEquitySwap;
/** simulated path values */
class PRODUCTS_DLL EquitySwapCreditSupport: virtual public CreditSupport
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

    /** return instrument discount curve */
    virtual string getInstCcyCode() const;

    /** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const;


    EquitySwapCreditSupport(CInstrument*, CMarketDataSP market);



private:
    IModelSP model;
    
    CEquitySwap* instrESWOrig;

    smartPtr<CEquitySwap> instrESW;

    // [credit dates][samples dates]
    // for each value date, arrays of deltas with respect to historical samples
    vector<vector<double> > avInDelta;
	vector<vector<double> > avOutDelta;
    vector<vector<double> > eqStartDelta;
    vector<vector<double> > eqEndDelta;

    // price array at each date
    DoubleArray priceCache;
    DoubleArray spotDelta;

    /** Populate any samples in ESWEquity between date1 and date2 using linear interpolation, assuming that the spot level was
    value1 and value2 on these dates. */
    void rollLinearESWEquity(const DateTime& date1, double value1, 
                             const DateTime& date2, double value2);

    /** computes double array deltaCache[dates][samples] consisting of the delta 
    with respect to each historical sample */                                  
    void computeSampleDeltaCache(
        IObjectSP& inst,	
        const DateTimeArray& creditDates,
        CreditSupport *creditSupport,
        const IModelSP &model,
        double tweakSize,
        const DoubleArray &priceCache,
        const DateTimeArray& sampleDates,
        DoubleArray& sampleLevels,                 // levels to tweak.  remember to reset
        vector<vector<double> > &deltaCache);
    
    /** computes change from tweaked historical levels */
    double changeFromSamples(
        const DateTime& rolledDate,
        const DateTimeArray& sampleDates,
        const DoubleArray& sampleLevels,
        const vector<double>& sampleDelta,
        const double spotRef);

};

typedef smartConstPtr<EquitySwapCreditSupport> EquitySwapCreditSupportConstSP;
typedef smartPtr<EquitySwapCreditSupport> EquitySwapCreditSupportSP;

DRLIB_END_NAMESPACE

#endif
