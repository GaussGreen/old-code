//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : RiskyCDSCurve.hpp
//
//   Description : contains a yield curve and a spread curve.
//
//   Author      : André Segger
//
//
//----------------------------------------------------------------------------

#ifndef EDG_RISKY_CDS_CURVE_H
#define EDG_RISKY_CDS_CURVE_H

#include "edginc/Class.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/RiskyCurve.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL RiskyCDSCurve : public YieldCurve,
                      virtual public IRiskyCurve
{
public:
    static CClassConstSP const TYPE;
    friend class RiskyCurveHelper;

    RiskyCDSCurve(
		const string&         name, 
		const IYieldCurve&    ccyCurve, 
		const ICDSParSpreads& spreadCurve, 
		const DateTime&       today);

    RiskyCDSCurve(
        const string&           name, 
        const IYieldCurve&      ccyCurve, 
        const CleanSpreadCurve& spreadCurve, 
        double                  recovery);

    virtual string getCcy() const;
    virtual string getName() const;

    virtual DateTime getSpotDate() const;

    virtual DateTime getToday() const;

    /** Useful accessor methods */
    virtual ExpiryArrayConstSP getExpiries() const;
    virtual StringArrayConstSP getInstruments() const;

    /** Returns tradeDate + n business days where n = spot offset */
    virtual DateTime settles(const DateTime& tradeDate) const;

    virtual double pv(const DateTime& lodate, 
                      const DateTime& hidate) const;

    virtual double pv(const DateTime& date) const;

    virtual double zero(const DateTime& date) const;

    /** Returns rateMaturity->toDate(rateStart) bad day adjusted
        either using bad day convention supplied (if non null)
        together with the yield curve holidays or the yield curve's
        bad day convention and holidays */
    virtual DateTime rateMaturity(
        const DateTime&         rateStart,
        const Expiry*           rateMaturity,
        const BadDayConvention* bdc) const; // optional

    virtual double fwd(const DateTime&           lodate, 
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const;

    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

    virtual double fwd(const DateTime&           payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IKey* logOfDiscFactorKey() const;

    /** calculate the risky PV based on the assumption that the bond recovers 
        recovery * recoveryNotional in case of default */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const;

    /** this function calculates the discount factor based on the assumption that the
        on which the recovery is based is provided externally. This allows to use
        different methodologies (PV, face value + accrued etc.) to be included easily  -
        this function will use the externally given recovery rather than the underlying
        risky curve's recovery */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional,
                           bool            useAssetRecovery,
                           double          assetRecovery) const;

    /** make a risky curve from a credit spread curve */
    virtual IYieldCurveSP makeRiskyCurve(
        const CreditSpreadCurve& spreadCurve,
        const  DateTime*    maturityDate=NULL) const;

	virtual CreditSpreadCurveSP makeCreditSpreadCurve(
		const string&        name,
	    const CashFlowArray& defaultRates,
		double               recovery) const;

    virtual CashFlowArraySP getRatesAndDates() const;

    virtual IYieldCurveSP createForwardCurve(const DateTime& forwardDate) const;

    /** grab the dates used in the zero curve */
    virtual DateTimeArray zeroDates() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest. Currently this implementation 
        fails */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;
    
    /** drive which style of zero curve is used. Default is discounting */
    void setProjectionCurve(bool useEstimatingCurve = true) const {};

private:
    class LogOfDiscFactorKey;

    // data members
    IYieldCurveConstSP   riskFreeCurve;
    string               name;

    // private methods
    RiskyCDSCurve();

    // transient fields
    CleanSpreadCurveConstSP cleanSpreadCurve;
    double                  recovery;
};

typedef smartConstPtr<RiskyCDSCurve> RiskyCDSCurveConstSP;
typedef smartPtr<RiskyCDSCurve>      RiskyCDSCurveCurveSP;

DRLIB_END_NAMESPACE
#endif

