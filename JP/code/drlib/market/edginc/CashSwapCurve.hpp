//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashSwapCurve.hpp
//
//   Description : A classical yield curve of cash & swap rates
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CASHSWAPCURVE_HPP
#define EDR_CASHSWAPCURVE_HPP

#include "edginc/BootstrappedYieldCurve.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/**
 * A classical yield curve of cash, futures & swap rates.
 *
 * This class preserves the old (and deprecated) interface for bootstrapping
 * zero curves - all new users should use BootstrappedCurve instead.
 */

class MARKET_DLL CashSwapCurve : public BootstrappedYieldCurve
{
public:
    static CClassConstSP const TYPE;
    friend class CashSwapCurveHelper;

    // (deprecated)
    CashSwapCurve(
        const string&                ccy,
        const string&                name,
        const DateTime&              today, 
        int                          spotOffset,
        const HolidayWrapper&        hols,
        int                          moneyMarketDenom,
        const string&                swapDayCount,
        int                          swapFrequency,
        const string&                swapFloatDayCount,
        int                          swapFloatFrequency,
        const string&                swapBadDayConvention,
        const ExpiryArray&           expiries,
        const DoubleArray&           rates,
        const StringArray&           instruments,
        const IZeroCurveFactorySP    zcMethod,
        const CurrencyBasisWrapper&  ccyBasis,
        const string&                futuresMaturity,
        const IRVolBaseWrapper&      irVol,
        bool                         isIndexCurve = false);

    virtual ~CashSwapCurve();
    
    virtual IPublicObject* toPublicObject() const;

    /** Define the public interface to this class */
    class MARKET_DLL Interface : public CObject, public IPublicObject {
    public:
        static CClassConstSP const TYPE;
        static const int NULL_FREQUENCY;

        friend class CashSwapCurveInterfaceHelper;
        virtual IPrivateObject* toPrivateObject() const;

        Interface(
            const string&                ccy,
            const string&                name,
            const DateTime&              today, 
            int                          spotOffset,
            const HolidayWrapper&        hols,
            int                          moneyMarketDenom,
            const string&                swapDayCount,
            int                          swapFrequency,
            const string&                swapFloatDayCount,
            int                          swapFloatFrequency,
            const string&                swapBadDayConvention,
            const ExpiryArray&           expiries,
            const DoubleArray&           rates,
            const StringArray&           instruments,
            const IZeroCurveFactory*     interp,    //maintain old name for this public parameter
            const CurrencyBasisWrapper&  ccyBasis,
            const string&                futuresMaturity,
            const IRVolBaseWrapper&      irVol,
            bool                         isIndexCurve);

    private:
        Interface();

        string                  ccy;
        string                  name;
        DateTime                today; 
        int                     spotOffset;
        HolidayWrapper          hols;
        int                     moneyMarketDenom;
        string                  swapDayCount;
        int                     swapFrequency;
        string                  swapFloatDayCount;
        int                     swapFloatFrequency;
        string                  swapBadDayConvention;
        ExpiryArray             expiries;
        DoubleArray             rates;
        StringArray             instruments;
        IZeroCurveFactorySP     interp; //maintain old name for this public parameter
        CurrencyBasisWrapper    ccyBasis;
        string                  futuresMaturity;
        IRVolBaseWrapper        irVol;
        bool                    isIndexCurve;
    };

private:
    CashSwapCurve();
    CashSwapCurve(const CashSwapCurve &rhs);
    CashSwapCurve& operator=(const CashSwapCurve& rhs);

    // build zero curve for current yield curve
    virtual ZeroPairSP zeroCurve() const;    

    // build arrays for old public interface
    void toArrays(
        ExpiryArray&    expiries,
        DoubleArray&    rates,
        StringArray&    instruments,
        const DateTime& refDate,
        const Holiday&  holidays) const;
};

typedef smartPtr<CashSwapCurve> CashSwapCurveSP;
typedef smartConstPtr<CashSwapCurve> CashSwapCurveConstSP;


DRLIB_END_NAMESPACE
#endif
