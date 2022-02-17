/*****************************************************************************
 *
 *    Group       : Equity Derivatives Research
 *
 *    Description : Defines a Payment stream
 *
 *    Author      : Stephen Hope
 *
 *    Date        : 23 May 2001
 *
 *
 *******************************************************************************/

#ifndef _EDG_PAYSTREAM_H
#define _EDG_PAYSTREAM_H

#include "edginc/Object.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/FloatRate.hpp"
#include "edginc/IXLInterfaceMap.hpp"
#include "edginc/AccrualCalendar.hpp"


DRLIB_BEGIN_NAMESPACE

class PayStreamXLInterface;

class MARKET_DLL PayStream: public CObject,
                            public IXLInterfaceMap{
public:
    static CClassConstSP const TYPE;
    friend class PayStreamHelper;
    friend class PaymentHelper;
    friend class RefixHelper;

    // * for class loading */
    static bool load();  

    /** overrides default */
    virtual void validatePop2Object(); 

    /** Map the object */
    virtual IObjectSP map()const;

    /** Get the floating rate curve from the market data cache */
    void getMarketData(const IModel* model, const CMarketDataSP market);

    /** Constructor with flat list type interface */
    PayStream(const DayCountConvention* payDayCountConv,
              bool  compoundFlat,
              const FloatRate* floatRate,
              const DoubleArray& notionals,
              const DateTimeArray& refixDates,
              const DateTimeArray& accruelDates,
              const DateTimeArray& payDates,
              const DoubleArray&   refixLevels,
              const DoubleArray&   spreads,
              const DoubleArray&   weights,
              const BoolArray&     compound);

    /** Returns the floatRate Ccy */
    string getCcy() const;

	/** Returns the floatRate YC Name */
	string getYCName() const;

    /** Curtail the PayStream for callable LiborStreams */
    void curtail(const DateTime& callDatePlusOffset);

    /** Return last accrue start date */
    const DateTime& getLastAccrueDate()const;
    
    /** Calculate the fwd rate between the callDatePlusOffset and
        the accrueEnd date of the relevant refix period. Note that
        the discount curve is passed since this may be a fixed PayStream
        i.e FloatRate is empty ! */
    const double callPlusOffsetFwdRate(const DateTime& callDatePlusOffset,
                                       const YieldCurveWrapper& discount)const;

    /** Calculate the cash flow list arising from a payment stream.
        The number of notionals must match the number of payments */
    CashFlowArrayConstSP cashFlowList(
              const DateTime& today,          /* (I) value date */
              const DoubleArray* notionals);  /* (I) Use NULL to use notionals inside payment
                                                 - otherwise these override it */

    /** As above but uses supplied fixings */
    void cashFlowList(const DateTime&    today,     // (I) 
                      const DoubleArray* notionals, // (I), optional
                      const DoubleArray& fixings,   // (I)
                      CashFlowArray&     cfl);      // (M) the result - done this way it avoids memory alloc

    /** Calculate the cash flow list arising from a payment stream.
        The number of notionals must match the number of payments.
        For use when the parent instrument has a call date */
    CashFlowArrayConstSP cashFlowList(
        const DateTime& today,                /* (I) value date */
        const DateTime& callDatePlusOffset,   /* (I) call date + offset */
        double calFwdRate,                    /* (I) fwd rate between callDatePlusOffset and accrueEnd */
        const DoubleArray* notionals);        /* (I) Use NULL to use notionals inside payment
                                                 - otherwise these override it */
    
    /** Calculate the cash flow list arising from a payment stream.,
        but only up to the cut off payment date passed in. i.e
        generate a cash flow list from a stream within a stream */
    CashFlowArrayConstSP cashFlowList(const DateTime& fromDate,
                                      const DateTime& lastPayDate)const;

    /** Calculate the cash flow list arising from a payment stream 
        containing payments after 'fromDate' */
    CashFlowArrayConstSP cashFlowListAfterDate(
        const DateTime& today,          /* (I) value date */
        const DateTime& fromDate);      /* (I) date from */
    
        /* Validation for libor flow. Checks that the payment dateTime is after
       the last refix dateTime for the payment period. */
    void checkRefixAndPayment()const;

    /** Invoked by THETA shift method. Roll the historical refixes 
        a day forward */
    void rollRefixes(const DateTime& valueDate,
                     const DateTime& newValueDate);

    /** return TRUE if floatingRate */
    bool isFloating()const;

    /** return the ACCRUED_INTEREST */
    double getAccruedInterest()const;

    // get the accrued as of 'when'
    double getAccruedInterest(const DateTime&    today,
                              const DateTime&    when,
                              const DoubleArray* notionals) const;

    /** when to stop tweaking */
    DateTime lastSensDate()const;

    /** when do payments occur ? */
    DateTimeArraySP paymentDates() const;

    /** and what are they ? Use NULL to use notionals inside payment */
    CashFlowArrayConstSP knownCashflows(const DateTime&    today,
                                        const DoubleArray* notionals,
                                        bool               isCallable,
                                        const DateTime&    cutoff,
                                        double             callFwdRate);

    /** returns cashflows corresponding to periods whose accrual start date is before cutoff */
    CashFlowArrayConstSP startedKnownCashflows(const DateTime&    today,
                                               const DoubleArray* notionals) const;

    AccrualCalendarArraySP couponDue(const DateTime&        fromDate,
                                     bool                   isCallable,
                                     const DateTime&        cutoff,
                                     double                 callFwdRate) const;

    AccrualCalendarArraySP accrualCalendar(const DateTime& date) const;

private:
    PayStream();
    PayStream(const PayStream &rhs);
    PayStream& operator=(const PayStream& rhs);
public:
    //////// Payment Class  //////////////
    class MARKET_DLL Payment: public CObject{
    public:
        static CClassConstSP const TYPE;
        friend class PayStream;
        friend class PaymentHelper;
        friend class RefixHelper;
    
        Payment();

        Payment(double notional,
                const DateTimeArray& refixDates,
                const DateTimeArray& accruelDates,
                const DateTime& payDate,
                const DoubleArray& refixLevels,
                const DoubleArray& spreads,
                const DoubleArray& weights,
                int   startIdx,
                int   endIdx);

        /** populate the passed cash flow array with the refix dates and
            levels within the payment */
        CashFlowArraySP getRefixes()const;

        /** Returns the floatRate Ccy */
        string getCcy()const;

        /** do we have all the data today to calculate this ? */
        bool isKnown(const DateTime& today, bool isFloating) const;

        /** is the payment known and started accruing? */
        bool isKnownAndStartedAccruing(const DateTime& today, bool isFloating) const;

        /** overrides default */
        virtual void validatePop2Object();

        /** Payment validation. Invoked by PayStream validatePop2Object */
        void validate(bool isFloating, const bool isAdjust)const;

        /** Calculate the cash flow arising from a payment.
            If we have a FloatRate then the future rates will be calculated, 
            otherwise the rates inside the payment will be used.
            For the notional, either pass NULL to use the notionals recorded inside
            the payments, otherwise pass a pointer to a double containing the single
            notional to use for that payment. */   
        CashFlow cashFlow(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const FloatRate* floatRate,                  /* (I) NULL unless payStream is a float type */
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional)const;                /* (I) Use NULL to use notional inside
                                                            payment otherwise this overrides it */

        /** Calculate the cash flow arising from a payment.
            If call date + offset is within the last refix period 
            then the following behaviour is invoked :
            (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
            (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
            Subtract (A) from (B) */
        CashFlow cashFlow(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const FloatRate* floatRate,                  /* (I) NULL unless payStream is a float type */
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional,                      /* (I) Use NULL to use notional inside
                                                            payment otherwise this overrides it */
            const DateTime&  callDatePlusOffset,         /* it is callable --> supply the date of call date + offset */ 
            const double callFwdRate)const;              /* fwd rate between callDatePlusOffset and final accrueEnd date */

        /** Invoked indirectly by THETA shift method. Rolls the 
            historical refixes for a payment */
        void rollRefixes(const DateTime& valueDate,
                         const DateTime& newValueDate,
                         const FloatRate& floatRate);

        // this is here just to allow ESW libor object access this class
        // alternatively add constructors, make fields public, or ?
        friend class ESWLibor;
    private:
        Payment(const Payment &rhs);
        Payment& operator=(const Payment& rhs);

        // like the other flavours except fixings are passed in explicitly

        CashFlow cashFlow(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const double* fixings,                       // (I) fix/flt fixings.
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional)const;                /* (I) Use NULL to use notional inside */

        // like the other flavours except fixings are passed in explicitly
        CashFlow cashFlowStandard(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const double* fixings,                       // (I) fix/flt fixings.
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional)const;                /* (I) Use NULL to use notional inside
                                                            payment otherwise this overrides it */

        CashFlow cashFlowBrazil(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const double* fixings,                       // (I) fix/flt fixings. "double*" allows more flexibility
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* theNotional)const;                 /* (I) Use NULL to use notional inside*/

        /** Calculate the cash flow arising from a payment.
            If call date + offset is within the last refix period 
            then the following behaviour is invoked :
            (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
            (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
            Subtract (A) from (B) */
        CashFlow cashFlowStandard(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const FloatRate* floatRate,                  /* (I) NULL unless payStream is a float type */
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional,                      /* (I) Use NULL to use notional inside
                                                            payment otherwise this overrides it */
            const DateTime&  callDatePlusOffset,         /* it is callable --> supply the date of call date + offset */ 
            const double callFwdRate)const;              /* fwd rate between callDatePlusOffset and final accrueEnd date */


        /** Calculate the cash flow arising from a payment.
            If call date + offset is within the last refix period 
            then the following behaviour is invoked :
            (A) Accrue over the Notional from accrueStart to accrueEnd as usual.
            (B) Accrue over the Notional from call date + offset to accrueEnd using the fwd rate at call date + offset.
            Subtract (A) from (B) */
        CashFlow cashFlowBrazil(
            const DayCountConvention* payDayCountConv,   /* (I) Pay date count convention */
            bool  compoundFlat,                          /* (I) TRUE = Index rate flat,
                                                            FALSE = Index rate + spread */  
            const FloatRate* floatRate,                  /* (I) NULL unless payStream is a float type */
            const DateTime& today,                       /* (I) value date */                       
            double* accInt,                              /* (M) accruedInterest so far */
            const double* notional,                      /* (I) Use NULL to use notional inside
                                                            payment otherwise this overrides it */
            const DateTime&  callDatePlusOffset,         /* it is callable --> supply the date of call date + offset */ 
            const double callFwdRate)const;              /* fwd rate between callDatePlusOffset and final accrueEnd date */

        
    public:
        /////// Refix Class //////////////
        class MARKET_DLL Refix: public CObject{
        public:
            static CClassConstSP const TYPE;
            friend class PayStream;
            friend class Payment;
            friend class RefixHelper;

            Refix();
            
            Refix(const DateTime& refixDate,
                  const DateTime& accrueStart,
                  const DateTime& accrueEnd,
                  double refixLevel,
                  double spread,
                  double weight);

            /** overrides default */
            virtual void validatePop2Object();

            /** Mainly validates against pricing in arrears for floating types */
            void validate(bool isFloating, const bool isAdjust)const;

            /** Invoked indirectly from THETA shift method. 
                Rolls historical refixes */
            void rollRefix(const DateTime& valueDate,
                           const DateTime& newValueDate,
                           const FloatRate& floatRate);

            // this is here just to allow ESW libor object access this class
            // alternatively add constructors, make fields public, or ?
            friend class ESWLibor;

        private:
            Refix(const Refix &rhs);
            Refix& operator=(const Refix& rhs);

            DateTime     refixDate;
            DateTime     accrueStart;
            DateTime     accrueEnd;
            double       refixLevel;
            double       spread;
            double       weight;
        };
        /////// End of Refix Class //////////////
        typedef smartPtr<Refix> RefixSP;
        typedef array<RefixSP, Refix> RefixArray;
        typedef smartPtr<RefixArray> RefixArraySP;
        typedef smartConstPtr<Refix> RefixConstSP;
  
      
    private:
        DateTime     payDate;
        RefixArraySP refixes;
        double       notional;
    };
    //////// End of Payment Class  //////////////
    typedef smartPtr<Payment> PaymentSP;
    typedef array<PaymentSP, Payment> PaymentArray;
    typedef smartPtr<PaymentArray> PaymentArraySP;
    typedef smartConstPtr<Payment> PaymentConstSP;
private:

private:
    PaymentArraySP            payments;
    FloatRateSP               floatRate;
    DayCountConventionConstSP dayCountConv;
    bool                      compoundFlat;

    // not registered
    double   accruedInterest;
};


typedef smartConstPtr<PayStream> PayStreamConstSP;
typedef smartPtr<PayStream> PayStreamSP;
typedef array<PayStreamSP, PayStream> PayStreamArray;

DRLIB_END_NAMESPACE

#endif
