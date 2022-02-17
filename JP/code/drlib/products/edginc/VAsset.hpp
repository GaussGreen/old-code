//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VAsset.hpp
//
//   Description : Implement volatility asset such as VIX 
//                 
//
//   Author      : xiaolan zhang
//
//   Date        : 24 Mar 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VAsset_HPP
#define EDG_VAsset_HPP

#include "edginc/VolBase.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/VolVarSwap.hpp" //to use var swap pricer

DRLIB_BEGIN_NAMESPACE

// VIXFModel class 
// same as VanVSModel with 1 additional input
class VIXFModel : public VanVSModel {
public:
    static CClassConstSP const TYPE;
            
    /** the class that the product must be able to create */
    class IProduct{
    public:
        virtual void price(VIXFModel* model,
                           Control*   control, 
                           Results*   results) const = 0;
        virtual ~IProduct(){};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(const VIXFModel* model) const = 0;
    };

    // inherited from CModel 
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
     

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(VIXFModel, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultVIXFModel);

        FIELD(pricingMethod, "BS_SV, SV, BS_SV_NO_BASIS, SV_NO_BASIS, VSW_BS, VSW_SV");
        FIELD_MAKE_OPTIONAL(pricingMethod);                                                   
    }

    // for VIXFModel::IIntoProduct  
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(VIXFModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultVIXFModel(){
        return new VIXFModel();
    }
    // constructor 
    VIXFModel(CClassConstSP type=TYPE): VanVSModel(type), pricingMethod("BS_SV"){}

    //returns pricingMethod used to compute prices VIX Future
    string getPricingMethod() const{return pricingMethod;};

public:
    static const string BS_SV ;        
    static const string SV;   
    static const string BS_SV_NO_BASIS;        
    static const string SV_NO_BASIS;   
    static const string VSW_BS;   
    static const string VSW_SV;   

private :
    /** registration */
    string                   pricingMethod; 

    /** BS_SV : BS + convexity + basis
        SV: E_t[sqrt(Fwd Var(T,T+ tau))] + basis
        BS_SV_NO_BASIS: BS + convexity
        SV_NO_BASIS: E_t[sqrt(Fwd Var(T,T+ tau))]
        VSW_BS: Var Swap using BS
        VSW_SV: Var Swap using SVJ
        */
    void validate();
};

typedef smartPtr<VIXFModel> VIXFModelSP;
typedef smartConstPtr<VIXFModel> VIXFModelConstSP;

/////////////////////////////////////////////////////////////////////////////////

class VAssetAlgorithm : public CObject {
public:
    static CClassConstSP const TYPE;
    VAssetAlgorithm() : CObject(TYPE) {}
    static IObject* defaultConstructor(){
        return new VAssetAlgorithm();
    }
    static void load(CClassSP& clazz){
        REGISTER(VAssetAlgorithm, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(model, "");
        FIELD(ctrl, "");
        FIELD(type, "");
    }

    VAssetAlgorithm(VIXFModelConstSP m, CControlSP c, const string& t)
        : CObject(TYPE), model(m), ctrl(c), type(t){}

    void  setType(const string calcType){type = calcType;} 
    void  setControl(CControlSP control){ctrl = control;}

    VIXFModelConstSP   model;
    CControlSP          ctrl;
    string              type;
};
typedef smartPtr<VAssetAlgorithm> VAssetAlgorithmSP;

////////////////////////////////////////////////////////////////////////////////////

/** VIX asset */
class VAsset: public IndexSpecEQ,
              //virtual public IRestorableWithRespectTo<Spot>,
              virtual public LastSensDate
{
//           virtual public DeltaSurface::IShift,

public:
    static const string PARVOL_SWAPVAR ;        
    static const string PARVOL_SWAPVAR_SV ;   
    static const string SQT_SWAPVAR_SV ;        
    static const string PARVOL_SWAPVAR_VIX_BASIS0 ;   

    static CClassConstSP const TYPE;
    friend class VAssetHelper;
    friend class VAssetFP;

    ~VAsset();

    /** Pull out the factor assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market) ;

    /** returns the spot price */
    /** if MTM = true, Market spot price; else, calculated spot price*/
    virtual double getSpot() const;

    /** returns the market spot price */
    double getMarketSpot() const;

    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** record forwards at maturity*/
    void recordFwdAtMat(OutputRequest*  request,
                        CResults*       results,
                        const DateTime& maturityDate) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Returns name identifying this object for Delta */
    virtual string sensName(Spot*) const;

    /** Shifts the object using given shift (see Delta::Shift)*/
    //virtual TweakOutcome sensShift(const PropertyTweak<Spot>& shift);

    /** Restores the object to its original form */
    //virtual void sensRestore(const PropertyTweak<Spot>& shift);
    
    /** Basis defined as the diff between Market value and model value for VAsset spot*/
    double getBasis(string basisType) const ;

    /** get fixing for a given date */
    double getFixing (DateTime& t, string source) const;

    VAssetAlgorithmSP getAlgorithm() const {
        return data;
    }

    ExpirySP getTenor() const {return tenor;}

    /** reprot sensitive strikes through var swap instrument*/
    DoubleArraySP sensitiveStrikes(OutputNameConstSP outputName, const VIXFModel* model) const;

    /** vector "out" contains BS, SV components, only has values after calling vAsset->fwdValue*/
    /** record "out" into results based on the control */
    const void recordExtraOutput(Control*   control,  CResults*  results) const;

private:
    VAsset();
    VAsset(const VAsset& rhs);
    VAsset& operator=(const VAsset& rhs);

    /** This is the method responsible for pricing the Fwd VarSwap 
        and the convexty adj using SVJ and aggregating the results
    
        compute fwd price according to model flag. all component values 
        are output to prices 
    */
    virtual double fwdValue(const DateTime& date, DoubleArray& out) const;

    /** price each components based on data type, need to set the data->type before calling*/
    double priceFwd(const DateTime& date) const;

    /** return BS model price */
    double priceBS(const DateTime& startDate) const;
    double calcPriceBS(const DateTime& startDate, const DateTime& matDate)const;

    /** return stochastic vol model price*/
    double priceSV(const DateTime& startDate) const;
    double calcPriceSV(const DateTime& startDate, const DateTime& matDate)const;

    double getBSSpot()const;
    double getSVSpot()const;

    /** create a var swap instrument*/
    VarianceSwapSP createVarSwap(const CAssetWrapper& asset, const YieldCurveWrapper& yc, 
                                const VAssetAlgorithm* data, const DateTime& startDate, 
                                const DateTime& matDate) const;

    //void recordExtraOutput( CResults*       results) const;

/////////////////////////////////////////////////////////
    //  exported fields
    double                  spot;               //VAsset's spot
    ExpirySP                tenor;			    // e.g. 30 days as default value
    bool                    useMTM;

    // transient fields
    CAssetWrapper           undSV;               //factor with SV
    VAssetAlgorithmSP       data;
    double                  basisBS;             //spot - varSwap price (BS)
    double                  basisSV;             //spot - varSwap price (SV)

    DoubleArray   out;                           // store components for output purpose
};

typedef smartPtr<VAsset> VAssetSP;
typedef smartConstPtr<VAsset> VAssetConstSP;

DRLIB_END_NAMESPACE
#endif
