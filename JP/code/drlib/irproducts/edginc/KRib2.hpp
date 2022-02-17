//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRib2.hpp
//
//   Description : Rib component
//
//----------------------------------------------------------------------------
#ifndef _KRib2_HPP
#define _KRib2_HPP

#include "edginc/KComponent.hpp"
#include "edginc/FlexDates.hpp"

DRLIB_BEGIN_NAMESPACE

/*------------------------------------------------------------------------*/

class IRPRODUCTS_DLL KRib2 : public KComponent,
                             virtual public FDModel::IIntoProduct 
{
public:
    static CClassConstSP const TYPE;    
    friend class KRib2Tree;

    /****************** exported fields ************/
    CModel::IProdCreatorArray obsUnds;  // can be set by shell instrument

    FlexDates obsDates;
    bool includeAccStartObs;

    /**************** transient fields ************/
    // These fields need to be set by the parent component
public:
    DateTimeArray accStart;
    DateTimeArray accEnd;
    DateTimeArray resetDates;
    DateTimeArray payDates;
    string discountCurveName;
private:
    IntArray obsToCpn;

    /****************** methods ************/
public:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);

    virtual double getRibFraction(int cpnIdx, CashflowInfo &cfi ) const;
protected:
    void assignObsDatesToCoupons(
            IntArray &convTable,
            DateTimeArray const &dates) const;

private:
    KRib2(void) : KComponent(TYPE), includeAccStartObs(true) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KRib2(); }
};
DECLARE(KRib2);

/*------------------------------------------------------------------------*/

class FDModel;
class KRib2Tree : public FDProduct {
  
    /****************** Compressor ************/
public:
    class Compressor;
    typedef smartPtr<Compressor> CompressorSP;

    class IRPRODUCTS_DLL Compressor : public CObject,
                                      virtual public FDModel::IProdModifier {
    public:
        static CClassConstSP const TYPE;    

        class Tester : public CObject {
        public:
            static CClassConstSP const TYPE;    

            DateTime today;
            DateTimeArray inDates;
            IntArray obsToCpn;
            CompressorSP compressor;
            
            IObjectSP run();
        private:
            Tester() : CObject(TYPE) {}
            static IObject* defaultConstructor(void) { return new Tester(); }
            static void load(CClassSP& clazz);
        };

        /****** exported fields ***/
    private:
        ExpirySP verbatimPeriod;
        MaturityPeriodSP targetFreq;
        MaturityPeriodSP decayHalfLife;

        /****** methods ***/
    public:
        void lop(DateTime today,
            DateTimeArray const &inDates, DateTimeArray &outDates, 
            IntArray const &obsToCpn);

    private:
        Compressor() : CObject(TYPE) {}
        static IObject* defaultConstructor(void) { return new Compressor(); }
        static void load(CClassSP& clazz);
    };

    /********************* methods ********************/
public:
    KRib2Tree(const KRib2ConstSP &inst, FDModel* model);

    virtual string getOutputName() const { return inst->outputName; }
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int &step, FDProduct::UpdateType updateType);

    void couponHasReset(TreeSliceSP couponFixing, int cpnIdx);
    // Returns the RIB coupon. 
    // Throws an error if all observations date have not been processed for this coupon
    // If couponHasReset() has not been called yet, the result must be multiplied by the reset to 
    // have the actual coupon value.
    void getRibCoupon(int step, int couponIdx, TreeSliceSP &ribFraction);

    // Stops DEV and deletes slices involved in a coupon (call when not needed anymore)
    void deleteCoupon(int couponIdx);

    // irrelevant functions
    virtual const TreeSlice& getValue(int step, DateTime eventDate) const {
        throw ModelException(__FUNCTION__, "This component does not return a price");
    }
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);

protected:
    void addObservation(int step, TreeSlice &accuSlice, int &nbObs, int couponIdx, int obsIdx);

    /******************** variables ********************/
public:
    // Store rib fraction (not divied by the number of observation) 
    // of the couponFixings[] (or 1 if not yet available)
    vector<TreeSliceSP> ribCoupons;  
    vector<TreeSliceSP> couponFixings;
    // number of observation made per coupon
    vector<int> nbObsPerCoupon;
private:
    KRib2ConstSP  inst;
    FDProductArray obsProd; // underlyings to observe
    DateTimeArray modelObsDates;
    IntArray modelObsToCpn;
    int nbCoupons;
};

typedef refCountPtr<KRib2Tree> KRib2TreeSP;


DRLIB_END_NAMESPACE

#endif


