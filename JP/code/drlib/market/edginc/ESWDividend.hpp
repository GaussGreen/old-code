//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWDividend.hpp
//
//   Description   dividend class for equity swap.
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ESW_DIVIDEND_HPP
#define EDG_ESW_DIVIDEND_HPP


#include "edginc/DividendList.hpp"
#include "edginc/ESWEquity.hpp"
#include "edginc/ESWLibor.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ESWDividend : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWDividend(){
        return new ESWDividend();
    }

    ESWDividend();
        
    /** validate Div leg */
    virtual void Validate();

   /** initialisation */
    void init(ESWEquityConstSP eq,
              ESWLiborConstSP libor,
              const YieldCurveSP&     discount,
              const DateTime& callDate, 
              bool isCallable, 
              const DateTime& valueDate);

    double calcAccrueFactor(ESWLiborConstSP libor, const DateTime& valueDate,
                            const DateTime& divPayDate, const DateTime& instDivPayDate,
                            YieldCurveSP discCurve) const;

    void setFixing(const DateTime& valDate, const DateTime& rollDate, ESWLiborConstSP libor);

    /** price dividend leg */
    double priceDiv(YieldCurveWrapper discount, ESWLiborConstSP libor, const DateTime& valueDate, 
                    CAssetConstSP eq, const string& ccyTreatment,
                    const DateTime& endDate, const DateTime& eqCallSettleDate,
                    bool isCallable,const Control* control, CResults* results,
                    OutputRequestUtil::KnownCashFlows& knownCFs, 
                    DateTimeArray& payment_dates) const;

    string          LegName;
    bool            isReceiving;

    int             accrueDivType;
    bool            useSynthetic;
    int             payStyle;

    DateTimeArray   startDates;
    DateTimeArray   payDates;

    DividendListSP  synthDivs;
    DateTimeArray   divAccrueStartDates;
    DoubleArray     divAccrueDF;
    double          divWeight;

private:
    // unregistered
    DividendListSP  DivsInUse; // $unregistered
    DoubleArray     NumOfShares; // $unregistered
    bool            NoDivs; // $unregistered

    // rates for div accrual
    DoubleArray     DivAccrueRates; // $unregistered
    bool isCumulative();
    
    // This method returns the first historical ex date whose pay date is strictly in the future 
    // If none can be found, it returns valueDate.                                                 
    // Note: this method assumes the supplied divArray is sorted by ExDate.                      
    DateTime earliestUnpaidHistoricalExDate(const DateTime&      valueDate, 
                                            const DividendListSP& divList) const;

    void getDivInUse(ESWEquityConstSP         eq,
                     const YieldCurveSP&      discount,
                     const ESWLiborConstSP&   libor,
                     const DateTime&          start, 
                     const DateTime&          end, 
                     bool                     isCallable,
                     bool                     shouldValidate);
    DateTime actualDivPayDate(const Dividend& div) const;
    class DivAdjuster;
    friend class DivAdjuster;
    void adjustDividend(const YieldCurveSP&     discount,
                        const ESWLiborConstSP&  libor,
                        const ESWEquityConstSP& eq,
                        bool                    isCallable,
                        Dividend&               div);
};

typedef smartPtr<ESWDividend> ESWDividendSP;

DRLIB_END_NAMESPACE
#endif
