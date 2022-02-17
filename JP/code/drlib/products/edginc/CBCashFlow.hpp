//----------------------------------------------------------------------------
//
//   Group       : Convertible/EDG Derivatives Research
//
//   Filename    : CBCashFlow.hpp
//
//   Description : CBCashFlow
//
//   Date        : 2 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CBCASHFLOW_HPP
#define EDR_CBCASHFLOW_HPP

#include "edginc/CashFlow.hpp"
#include "edginc/InstrumentSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CBCashFlow : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultCBCashFlow(){
        return new CBCashFlow();
    }

    /** returns cash flows after start date (inclusive) 
     optionally returns accrued to start date */
    virtual CashFlowArraySP getCashFlow(const DateTime& start, double* accrued=0, bool incRemdeption=true) const;
    
    /** calculate cash flows for the list of dates */
    virtual void getCashFlow(const DateTimeArray& dateList, DoubleArray& result) const;

    /** calculate accrued insterest for the list of dates */
    virtual void getAccrued(const DateTimeArray& dateList, DoubleArray& result) const;

    /** calc yield to maturity */
//    virtual double YieldToMat();

    // acceessors
    DateTime getMaturity() const {return MaturityDate;};
    double getFaceValue() const {return FaceValue;};

protected:

    CBCashFlow();

	double		        	FaceValue;
	double			        CouponPct;
	double		        	Redemption;
    string                  BondDCC;
	string                  CouponFreq;
    DateTime                IssueDate;
    DateTime                MaturityDate;
    string                  CashFlowType;

    bool                    useInputCashFlow;
    DateTimeArray           CashFlowDates;
    DoubleArray             CashFlowAmounts;

    friend class CConvBond;
    friend class CConvBondTree1fProd;
};

typedef smartPtr<CBCashFlow> CBCashFlowSP;

DRLIB_END_NAMESPACE
#endif
