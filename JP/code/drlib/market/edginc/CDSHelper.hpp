//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CDSHelper.cpp
//
//   Description : Credit default swap helper class
//
//   Author      : Andr�Segger
//
//   Date        : 14 March 2003
//
//
//----------------------------------------------------------------------------

#ifndef CDS_HELPER_HPP
#define CDS_HELPER_HPP

#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/DefaultRates.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FirmAsset);
FORWARD_DECLARE(CDSParSpreads);
FORWARD_DECLARE(ICDSParSpreads);

class MARKET_DLL CDSHelper {
public:

    FORWARD_DECLARE_REF_COUNT(CParSpreadDefaultRates);

    class CFirmAssetDefaultRates;
    /** Main purpose is to get default rates from par spreads. Can
        also be created by directly supplying default rates. The par
        spreads creator actually calls the direct default rates creator
        when iterating to find the default rates. This class is [almost]
        immutable ie should not be modified -- the only time this is allowed
        is by the bootstrapper which derives from this class. The
        immutability is required as references to this object are stored in
        caches. */
    class MARKET_DLL CParSpreadDefaultRates: public DefaultRates{
    public:

        //// Possible values for interpolation
        static const string FLAT_FORWARD;
        static const string LINEAR;

        CParSpreadDefaultRates(const DateTime&      valueDate);

        CParSpreadDefaultRates(
            const DateTime&      valueDate,
            const DateTimeArray& inDefaultDates,
            const CDoubleArray&  inDefaultRates);

        CParSpreadDefaultRates(
            const DateTime&      valueDate,
            const DateTimeArray& inDefaultDates,
            const CDoubleArray&  inDefaultRates,
            const string         interpStyle);

        CParSpreadDefaultRates(
            CDSParSpreadsConstSP     cdsParSpreads,
            YieldCurveConstSP        discount,
            const DateTime&          valueDate,
            const DateTime&          effDate,
            const BadDayConvention*  bdc,
            DayCountConventionSP     dcc);

        CParSpreadDefaultRates(
            const CreditSpreadCurve* creditSpreads,
            YieldCurveConstSP        discount,
            const DateTime&          valueDate,
            const double&            spreadRecovery,
            const DateTime*          maturityDate = 0);

        // Basically turns the CParSpreadDefaultRates into a cash flow
        // for the CLEAN_SPREAD_CURVE o/p request
        CashFlowArraySP getCleanSpreadCurve() const;

        double calcDefaultPV(
            const DateTime&      fromDate,
            const DateTime&      toDate) const;

        /** Returns a key used to optimise repeated calculations of
            default probabilities. The calc method for this key returns
            the natural logarithm of 1 - the default probability between
            the two dates (or equivalently the product of the forward
            default rate (continuous, Act/365F) and the negative year
            fraction (Act/365F) betweeen the two dates. Or another way of
            saying this is to say that it returns minus the integral of
            the forward default rate. The default implementation has no
            performance improvements. */
        virtual IKey* logOfDefaultPVKey() const;

        virtual ~CParSpreadDefaultRates();

        virtual DateTimeArraySP getDefaultDates() const;

        CDoubleArraySP getDefaultValues() const;

        CashFlowArraySP getDefaultRates() const;

        CParSpreadDefaultRatesSP createFwdContDefCurve(const DateTime& valueDate) const;

        CParSpreadDefaultRatesSP annualiseDefaultRates() const;

        /** Returns (forward) rates values */
        const CDoubleArray& getRates() const;

        /** Returns default dates values */
        const DateTimeArray& getDates() const;

        /** Returns the interpolation style */
        const string& getInterpolationStyle() const;


    protected:
        CParSpreadDefaultRates();

        /** Returns today */
        virtual const DateTime& getValueDate() const;
        /** Used for bootstrapping only: change value of default rate
            i to "rate". See comments at start of class about its
            use. Note that it must be used incrementally (ie index
            starts at 0 and then increases by 1)  */
        void setRate(double rate, int index);

        DateTime      valueDate;
        CDoubleArray  defaultRates;
        DateTimeArray defaultDates;
        string        interpolationStyle;

        void          cacheIntegratedRate();

    private:
        class         LogOfDefaultPVKey;
        friend class  LogOfDefaultPVKey;
        friend class  CFirmAssetDefaultRates;
        double        calcIntegratedRate(int  date, int& idx) const;

        CDoubleArray  integratedRate;
    };

    /** Used to get default rates from a firm asset market object.*/
    class MARKET_DLL CFirmAssetDefaultRates: public DefaultRates{
    public:
        // That constructor SHOULD BE REMOVED
        // Still used for compatibility reasons
        CFirmAssetDefaultRates(
            const FirmAsset*        inFirmAssetWrapper,
            CParSpreadDefaultRates  cleanLiquiditySpreads);

        /** Constructor */
        CFirmAssetDefaultRates(
            const FirmAsset*        inFirmAssetWrapper,
            DefaultRatesSP  cleanLiquiditySpreads);

        CFirmAssetDefaultRates(const CFirmAssetDefaultRates& rhs);
        CFirmAssetDefaultRates& operator=(const CFirmAssetDefaultRates& rhs);

        /* calculates the probability of a default due to the asset
           level hitting the barrier by diffusion */
        double calcDefaultPV(
            const DateTime&      fromDate,
            const DateTime&      toDate) const;

        // calculates the probability of a default due to a jump-to-zero event
        double calcJumpPV(
            const DateTime&      fromDate,
            const DateTime&      toDate) const;

        // Basically turns the CParSpreadDefaultRates into a cash flow
        // for the CLEAN_SPREAD_CURVE o/p request
        CashFlowArraySP getCleanSpreadCurve() const;

        DateTimeArraySP getDefaultDates() const;

        virtual ~CFirmAssetDefaultRates();

        /** Returns (forward) rates values */
        const CDoubleArray& getRates() const;

        /** Returns default dates values */
        const DateTimeArray& getDates() const;

    protected:
        /** Returns today */
        virtual const DateTime& getValueDate() const;

    private:
        // FIELDS
        FirmAssetSP firmAsset;
        DefaultRatesSP cleanLiquiditySpreads;
    };

    DECLARE_REF_COUNT(CFirmAssetDefaultRates);

    // ------------- static helper functions for CDS pricing -----------------

    /** converts a CDS par spread curve into a clean spread curve. */
    static CleanSpreadCurveSP getCleanSpreadCurve(
        const ICDSParSpreads& parSpreadCurve,
        const DateTime&       valueDate,
        bool                  returnFwdRates);

    /** create the CDS timeline for the numerical integration */
    static DateTimeArraySP createTimeLine(
        const DateTime&           valueDate,
        const DateTime&           effDate,
        const CashFlowArray&      feePayments,
        const DateTimeArray&      defaults,
        YieldCurveConstSP         discount);

    /** Create feePayment from frequency */
    static CashFlowArraySP calculateFeePayments(
        const DateTime&            startDate,
        const DateTime&            maturity,
        int                        frequency,
        const DayCountConvention*  dcc,
        double                     notional,
        double                     spread);

    /** Static function rather than a method because it is called by
        the par spreads default rate finder as well as the general
        pricing. */
    static double calculateRiskyCashFlowPV(
        const DateTime&                   valueDate,
        const DateTime&                   maturity,
        const CashFlowArraySP             feePayments,
        YieldCurveConstSP                 discount,
        const DefaultRates*               defRates,
        bool                              isE2C);

    /** Static function rather than a method because it is called by
        the par spreads default rate finder as well as the general
        pricing. */
    static void calculateDefaultPayments(
        const DateTime&                   valueDate,
        const DateTime&                   effDate,
        const CashFlowArraySP             feePayments,
        const DoubleArraySP               notionals, /* array of notional *
                                                        (1-recovery) */
        bool                              accrueFee,
        const DateTime&                   accruedEffectiveDate,
        const DateTime&                   maturity,
        const DefaultRates*               defRates,
        const DayCountConvention*         dcc,
        bool                              startFromEffDate,
        bool                              isE2C,
        const DateTimeArraySP             critDates,
        YieldCurveConstSP                 discount,
        double&                           accruedPayment,    // output
        double&                           defaultPayment);   // output

    /** Static function rather than a method because it is called by
        the par spreads default rate finder as well as the general
        pricing. */
    static void calculateDefaultPayments(
        const DateTime&                   valueDate,
        const DateTime&                   effDate,
        const CashFlowArraySP             feePayments,
        double                            notional,
        double                            recovery,
        bool                              accrueFee,
        const DateTime&                   accruedEffectiveDate,
        const DateTime&                   maturity,
        YieldCurveConstSP                 discount,
        const DefaultRates*               defRates,
        const DayCountConvention*         dcc,
        bool                              isE2C,
        bool                              startFromEffDate,
        bool                              fullTimeLine,
        double                            &accruedPayment,   // output
        double                            &defaultPayment,   // output
        double                            &cashFlowRiskyPV); // output

    /** Static function rather than a method because it is called by
        the par spreads default rate finder as well as the general
        pricing. */
    static void CredDefSwapFromDefRatesObject(
        const DateTime                 &valueDate,
        const DateTime                 &effDate,
        const CashFlowArraySP          feePayments,
        double                         notional,
        double                         recovery,
        bool                           accrueFee,
        const DateTime&                accruedEffectiveDate,
        const DateTime&                maturity,
        YieldCurveConstSP              discount,
        const DefaultRates*            defRates,
        const DayCountConvention*      dcc,
        bool                           isE2C,
        double*                        price);   // output

    /** Static function that calculates an implied par spread for a
     * spot or forward CDS */
    static double CDSParSpread(
        const DateTime&                valueDate,
        const DateTime&                effDate,
        int                            frequency,
        double                         notional,
        double                         recovery,
        bool                           accrueFee,
        const DateTime&                accruedEffectiveDate,
        const DateTime&                maturity,
        YieldCurveConstSP              discount,
        bool                           isE2C,
        const DefaultRates*            defRates,
        const DayCountConvention*      dcc);

private:
};

typedef CDSHelper::CParSpreadDefaultRates ParSpreadDefaultRates;
typedef CDSHelper::CParSpreadDefaultRatesSP ParSpreadDefaultRatesSP;
typedef CDSHelper::CParSpreadDefaultRatesConstSP ParSpreadDefaultRatesConstSP;

DRLIB_END_NAMESPACE
#endif
