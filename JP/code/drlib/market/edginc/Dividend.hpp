//   Filename    : Dividend.hpp
//
//   Description : Dividend representation
//
//   Author      : Stephen Hope
//
//   Date        : 08 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef DIVIDEND_HPP
#define DIVIDEND_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL Dividend : public CObject {
public:
    static CClassConstSP const TYPE;

    static const int DIVIDEND_EXDIV_TIME;

    enum TDividendType
    {
        AMOUNT = 0,
        PERCENT,
        CONTINUOUS
    };
    
    friend class DividendHelper;

  
    /** Constructor */
    Dividend(const DateTime&      exDivDate,
             const DateTime&      payDivDate,
             const TDividendType  divType,
             const double         divAmount);

    /** Public default constructor to allow creation of DividendArray's */ 
    Dividend();

    /** Overridden for performance */
    virtual IObject* clone() const;

    /** Destructor */
    ~Dividend(){} // inline for performance reasons
    
    /* Interface methods */

    /** return the ex-dividend date */
    const DateTime&  getExDate()const;

    /** return the payment date */
    const DateTime&  getPayDate()const;

    /** return the dividend amount */
    double  getDivAmount()const;

    /** return the dividend type */
    Dividend::TDividendType  getDivType()const;

    /** tweak the dividend amount by shift size */
    void tweakDividend(double shiftSize);

    /** set the dividend amount */
    void setDivAmount(double newAmount);

    /** Sets div type to 'dollar' (ie AMOUNT) and scales amount by
        scaling factor provided */
    void convertToDollar(double scalingFactor);

    /** Sets div type to PERCENT and scales amount by
        scaling factor provided. To convert a dollar dividend
        to yield one should pass a scaling factor of (1.0 / fwd) */
    void convertToYield(double scalingFactor);

    /** Scaled dividend by supplied factor 
        ie divAmount -> divAmount* scalingFactor */
    void scale(double scalingFactor);

    /** Sets the pay date to supplied value */
    void setPayDate(const DateTime& payDate);
private:
    virtual void validatePop2Object();

    DateTime       exDivDate;
    DateTime       payDivDate;
    int            divType;
    double         divAmount; 
};

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<Dividend>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const Dividend& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(Dividend& value);

    /** Turns the IObjectSP into a Dividend */
    static const Dividend& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class MARKET_DLL arrayClone<Dividend>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

typedef smartPtr<Dividend>              DividendSP;
typedef smartConstPtr<Dividend>         DividendConstSP;
typedef array<Dividend>                 DividendArray;
typedef smartPtr<DividendArray>         DividendArraySP;
typedef smartConstPtr<DividendArray>    DividendArrayConstSP;

typedef Dividend CDividend;
typedef DividendArray CDividendArray;
typedef DividendArraySP CDividendArraySP;
typedef DividendArrayConstSP CDividendArrayConstSP;


DRLIB_END_NAMESPACE

#endif


