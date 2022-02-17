//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : QuasiVanillaModel.hpp
//
//   Description : First attempt at implementing quasi-vanilla in qlib
//                 subject to change in the future
//
//----------------------------------------------------------------------------

#ifndef QLIB_QuasiMQ_HPP
#define QLIB_QuasiMQ_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Model.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IRExoticParam.hpp"
#include "edginc/MarketTable.hpp"
#include "edginc/RateTree.hpp"   // ??? for linear/flat_fwd enum.  Move it when migrate to new enum style

DRLIB_BEGIN_NAMESPACE

// To avoid redundant file includes.
FORWARD_DECLARE(IndexSpec);
FORWARD_DECLARE(IndexSpecIR);
FORWARD_DECLARE(InstrumentSettlement);

/** QuasiMQ class. */
class TREE_DLL QuasiMQ : public CModel
{
public:
    static CClassConstSP const TYPE;

    /****************** methods **************/
    QuasiMQ(const CClassConstSP &type);
    static int errCallBack(const char *msg);
    static string popErrMsg();
    static ModelException makeException();

    struct VariateDebugData
    {
        VariateDebugData() : expiryDates(new DateTimeArray),
                             startDates(new DateTimeArray),
                             endDates(new DateTimeArray),
                             payDates(new DateTimeArray),
                             atmVols(new DoubleArray),
                             prices(new DoubleArray),
                             fwdRates(new DoubleArray)
        {}
        
        DateTimeArraySP expiryDates;
        DateTimeArraySP startDates;
        DateTimeArraySP endDates;
        DateTimeArraySP payDates;
        DoubleArraySP   atmVols;
        DoubleArraySP   prices;
        DoubleArraySP   fwdRates;
    };

    /** Univariate Pricer */
    void univarPricer(const IndexSpecIR&            rateSpec,        /*  1 (I) rate definition           */
                      const DateTime&               resetDate,       /*  2 (I) rate effective date       */
                      const DateTime&               payDate,         /*  3 (I) option payment date       */
                      double                        strike,          /*  4 (I) strike (0. for Q3_ADJ)    */
                      const DayCountConventionSP&   busDcc,          /*  5 (i) between today to reset date */
                      long                          optType,         /*  6 (I) option type               */
                      const DoubleArray&            payoffParams,    /*  7 (I) payoff coeffs             */
                      const InstrumentSettlementSP& setlType,        /*  8 (I) cash or phys settle       */
                      double&                       price,           /*  9 (O) fwd price and AA par rate */
                      double&                       fwdRate,         /* 10 (O) fwd price and AA par rate */
                      VariateDebugData&             debug) const;

    /** Univariate Pricer */
    void bivarPricer(const IndexSpecIR&            rate1Spec,       /*  1 (I) rate1 definition           */
                     const DateTime&               rate1Expiry,     /*  2 (I) rate1 effective date       */
                     const IndexSpecIR&            rate2Spec,       /*  3 (I) rate2 definition           */
                     const DateTime&               rate2Expiry,     /*  4 (I) rate2 effective date       */
                     const DateTime&               payDate,         /*  3 (I) option payment date       */
                     const DayCountConventionSP&   busDcc,          /*  5 (i) between today to reset date */
                     long                          optType,         /*  6 (I) option type               */
                     const DoubleArray&            payoffParams,    /*  7 (I) payoff coeffs             */
                     const InstrumentSettlementSP& setlType,        /*  8 (I) cash or phys settle       */
                     double&                       price,           /*  9 (O) fwd price and AA par rate */
                     double&                       fwdRate,         /* 10 (O) fwd price and AA par rate */
                     VariateDebugData&             debug1,
                     VariateDebugData&             debug2) const;

    // return discount factor
    double getZero(DateTime useDate, DateTime matDate, string curveName);

private:
    static void registerErrorCallback();

public:
    // CModel interface
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);
    virtual void Price(
        CInstrument* instrument,
        CControl*    control,
        CResults*    results);
    virtual WantsRiskMapping wantsRiskMapping() const { return IModel::riskMappingIrrelevant; }

    /****************** types **************/

    class Product
    {
    protected:
        QuasiMQ *model;
    public:
        virtual void price(Control* control, CResults* results) = 0;
        Product(QuasiMQ* model) : model(model) {}
        virtual ~Product() {};
    };
    typedef refCountPtr< Product > ProductSP;

    /****/

    // product creation interace
    class IIntoProduct : virtual public CModel::IModelIntoProduct
    {
    public:
        static CClassConstSP const TYPE;
    private:
        friend class QuasiMQ;
        virtual ProductSP createProduct(QuasiMQ* model) const = 0;
    };

    struct ModelMode {
        enum Enum {Q2Q, QMULTIQ};
    };
    typedef BoxedEnum<ModelMode::Enum> ModelModeBoxedEnum;

    /********************* variables *****************/
private:
    static string errMsg;

    /****************** exported fields **************/
public:
    IRCalibWrapper          irCalib;        // stuck with this for now
    string                  smileSet;       // same as tree IRCalib structure
    string                  modelSet;       // vnfm parameters
    ModelMode::Enum         modelMode;
    IRExoticParamTableWrapper smileTable;     // smile map collection
    IRExoticParamTableWrapper modelTable;     // model map collection
    ZeroInterpStyle::Enum   zeroInterpStyle;   // zero curve interpolation type

public:
    /****************** transient fields **************/
    // market data retrieved in getMarket, but has to be stored in the
    // model until the it is used in the eventual pricing call
    DateTime           today;
    DateTime           spotDate;
    YieldCurveWrapper  domesticYC;
    DoubleArray        vnfm;
    DateTimeArray      domDates;    // domesticYC split into
    DoubleArray        domRates;    // dates and rates
    IMarketFactorArray factors;     // market factors present in instrument
    IndexSpecArray     indexSpecs;  // array of index present in instrument
    bool               legacyPricingMode;  // discount to curveSpot date, don't adjust curves

private:
    static void load(CClassSP& );
    static IObject* defaultConstructor(void) { return new QuasiMQ(TYPE); }

    // Populates the smile data during each rib iteration.
    void populateSmileData(const DateTime& baseDate,
                           const DateTime& expiryDate,
                           const Expiry*   pMaturityTenor,
                           DoubleArray&    smileData) const; // output

    enum EnumSmileType
    {
        SMILE_MQ = 10,
        SMILE_SV = 20,
        SMILE_2Q = 30
    };
};

DRLIB_END_NAMESPACE

#endif
