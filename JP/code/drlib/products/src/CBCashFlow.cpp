//----------------------------------------------------------------------------
//
//   Group       : Convertible/EDG Derivatives Research
//
//   Filename    : CBCashFlow.cpp
//
//   Description : Convertible bond cash flow object
//
//   Date        : 2 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/CBCashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
void CBCashFlow::load(CClassSP& clazz){
	    clazz->setPublic();
        REGISTER(CBCashFlow, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCBCashFlow);
        FIELD(FaceValue, "face value of the bond");
        FIELD(CouponPct, "coupon rate as percentage of face value");
        FIELD(Redemption, "redemption value");
        FIELD(BondDCC, "bond day count convention");
        FIELD(CouponFreq, "coupon frequency: A, S, Q");
        FIELD(IssueDate, "bond issue date");
        FIELD(MaturityDate, "maturity date.");
        FIELD(CashFlowType, "cash flow type: BOND or PREFERRED");

        FIELD(useInputCashFlow, "true - use custom cash flow schedule, default=false");
        FIELD_MAKE_OPTIONAL(useInputCashFlow);
        FIELD(CashFlowDates, "user input cash flow dates");
        FIELD_MAKE_OPTIONAL(CashFlowDates);
        FIELD(CashFlowAmounts, "user input cash flow amounts");
        FIELD_MAKE_OPTIONAL(CashFlowAmounts);
}

CClassConstSP const CBCashFlow::TYPE = CClass::registerClassLoadMethod(
    "CBCashFlow", typeid(CBCashFlow), load);
bool  CBCashFlowLoad() {
    return (CBCashFlow::TYPE != 0);
   }


// constructor
CBCashFlow::CBCashFlow(): CObject(TYPE)
{
    useInputCashFlow = false;
};

/** returns cash flows after start date (inclusive) */
CashFlowArraySP CBCashFlow::getCashFlow(const DateTime& start, double* accrued, bool incRemdeption) const
{
    CashFlowArraySP cf (new CashFlowArray);

	int i;
    DateTime prevDate = IssueDate;

    DayCountConventionSP dcc(DayCountConventionFactory::make(BondDCC));

    if (start > MaturityDate || start < IssueDate)
        throw ModelException("CBCashFlow::getCashFlow", "cash flow start date cannot be after maturity or before issue date");

    if (useInputCashFlow)
	{// user supplied cash flow list
        if (CashFlowDates.size() != CashFlowAmounts.size())
            throw ModelException("CBCashFlow::getCashFlow", "cash flow dates and amounts array must be same size");

        int n = CashFlowDates.size();
		if (n == 0)
        {
            throw ModelException("CBCashFlow::getCashFlow", "no user cash flow schedule supplied");
        }
        else if (!(CashFlowDates[n-1] == MaturityDate))
        {
            throw ModelException("CBCashFlow::getCashFlow", "last cash flow date supplied must be maturity date");
        }
        else
        {
            for (i=n-1; i>=0; i--)
            {
                if (CashFlowDates[i] < start)
                {
                    prevDate = CashFlowDates[i];
                    break;
                }
                cf->push_back(CashFlow(CashFlowDates[i], CashFlowAmounts[i]));
            }
		}
	}
	else // calculate cashflows
	{
        DateTime prevDate = MaturityDate;
        do
        {// going from the back
            cf->insert(cf->begin(), CashFlow(prevDate, 0)); // leave amount 0 to be filled below
            prevDate = MaturityPeriod::toDate(-1, CouponFreq, prevDate);
        }while (prevDate >= start);
        
        if (prevDate <= IssueDate)
            prevDate = IssueDate;

        // calc cash flow amounts
        // front stub first
        (*cf)[0].amount = FaceValue*CouponPct*dcc->years(prevDate, (*cf)[0].date);
 		for (i=1; i<cf->size(); i++)
            (*cf)[i].amount = FaceValue*CouponPct*dcc->years((*cf)[i-1].date, (*cf)[i].date);

        if(!(CashFlowType == "PREFERRED") && incRemdeption)
		    (*cf)[cf->size()-1].amount += Redemption;
	}

    if(accrued!=0)
    {
        if(!(CashFlowType == "PREFERRED"))
            *accrued = FaceValue*CouponPct*dcc->years(prevDate, start);
	    else
		    *accrued = 0.0;
    }

    return cf;
}

/** calculate cash flows for the list of dates */
void CBCashFlow::getCashFlow(const DateTimeArray& dateList, DoubleArray& result) const
{
    if (dateList.size() == 0)
        throw ModelException("CBCashFlow::getCashFlow", "date list empty");
    if (dateList[dateList.size()-1] > MaturityDate)
        throw ModelException("CBCashFlow::getAccrued", "last date cannot be after maturity date");

    if (result.size() != dateList.size())
        result.resize(dateList.size());

    CashFlowArraySP cf = getCashFlow(dateList[0]);
    
    int j=0;
    for (int i=0; i<result.size(); i++)
    {
        result[i] = 0.0; // unmatched date will have 0 amount
        for (; j<cf->size(); j++)
        {
            if ((*cf)[j].date > dateList[i])
                break;
            else if ((*cf)[j].date == dateList[i])
            {
                result[i] = (*cf)[j].amount;
                j++;
                break;
            }
        }
    }
}

/** calculate accrued insterest for the list of dates */
void CBCashFlow::getAccrued(const DateTimeArray& dateList, DoubleArray& result) const
{
    if (dateList.size() == 0)
        throw ModelException("CBCashFlow::getAccrued", "date list empty");
    if (dateList[dateList.size()-1] > MaturityDate)
        throw ModelException("CBCashFlow::getAccrued", "last date cannot be after maturity date");

    if (result.size() != dateList.size())
        result.resize(dateList.size());

    CashFlowArraySP cf = getCashFlow(dateList[0], &*result.begin()); // first accrued done here
    bool isPreferred = (CashFlowType == "PREFERRED");

    if (isPreferred && !((*cf)[0].date == dateList[0]))
        result[0] = 0.0;

    DayCountConventionSP dcc(DayCountConventionFactory::make(BondDCC));

    DateTime prevDate = dateList[0];
    int j=0;
    for (int i=1; i<result.size(); i++)
    {
        for (; j<cf->size(); j++)
        {
            if ((*cf)[j].date > dateList[i])
                break;
        }
        if (isPreferred)
        {
            if (j>0 && (*cf)[j-1].date == dateList[i])
                result[i] = (*cf)[j-1].amount;
            else
                result[i] = 0.0;

        }
        else
        {
            if (j>0 && j<cf->size()-1)
                prevDate = (*cf)[j-1].date;

            result[i] = FaceValue*CouponPct*dcc->years(prevDate, dateList[i]);
        }
    }
}

DRLIB_END_NAMESPACE
