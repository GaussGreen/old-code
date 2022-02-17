//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaGridCEVJ.hpp
//
//   Description : Grid of European Vanillas for CEV+Jump Tree
//
//   Author      : Keiji Kitazawa
//
//   Date        : 07 Janualy 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VANILLA_GRID_CEVJ_HPP
#define EDR_VANILLA_GRID_CEVJ_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/Asset.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Tree1fCEVJ.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interest rate future */
class VanillaGridCEVJ : public Generic1Factor, 
                        public CTree1fCEVJ::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    // override base implementation if required
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    virtual void Validate();
 
    /** create a tree1f payoff product */
    virtual CTree1fCEVJ::IProduct* createProduct(CTree1fCEVJ* model) const;

    // below are copied from Vanilla
    virtual DateTime getValueDate() const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
   
 
private:
    friend class VanillaGridCEVJHelper;
    friend class VanillaGridCEVJ1fProd;

protected:
    static void load(CClassSP& clazz);
    VanillaGridCEVJ();
    VanillaGridCEVJ(CClassConstSP clazz);

    // this block is the same as Vanilla
    //string                  PayoffMode; // registered
    bool                    isCall; // not registered $unregistered
    //bool                    canExerciseEarly;
    ScheduleSP              exerciseSchedule;

    DoubleArray     StrikeArray;    //Strike Grid
    DoubleArray     TargetPrices;   //Target Price Array (1D array to be matrix) $unregistered
    DoubleArray     TargetVegas;   //Target Vegas Array (1D array to be matrix) $unregistered
    
};

/////////////////////////////////////////////////////////
//           tree1f/fd product
/////////////////////////////////////////////////////////
/** VanillaGridCEVJ product payoff for a tree */
class VanillaGridCEVJ1fProd:  virtual public CTree1fCEVJ::IProduct
{
public:
//  typedef enum{CALL, PUT, BINARY, FORWARD} TPayoffMode;
//  typedef enum{NA=-1, KI=0, KO=1} TBarrier;
//  typedef enum{KI_KEEP_KO, KI_CANCEL_KO, ONCE_TOUCH, TWO_TOUCH} TBarDependence;

    VanillaGridCEVJ1fProd(const VanillaGridCEVJ* instr):instrCEVJ(instr){
        NumOfPriceArray= instrCEVJ->StrikeArray.size();
    }

    virtual ~VanillaGridCEVJ1fProd(){}

    virtual CAssetConstSP GetAssetRef();

    virtual bool GetFwdStartLV()
    {
           return instrCEVJ->fwdStarting;
    }

    virtual DateTime GetFwdStartDateLV()
    {
		   if (instrCEVJ->fwdStarting)
               return instrCEVJ->startDate;
		   else
               return instrCEVJ->valueDate;
    }

    virtual YieldCurveConstSP GetDiscCurveRef()
    {
            return YieldCurveConstSP::dynamicCast(
                (IObjectConstSP)instrCEVJ->discount.getSP());
    }

    virtual CVolRequestConstSP GetLNRequest();

    /** initialise tree1f - allow product customisation */
    virtual void InitTree(CControl* control);

    /** initialise product specific data */
    virtual void InitProd();

    /** calculate at barriers for tree */
    virtual void preCalcTree(int step, int idx);

    /** product payoff method at maturity */
    virtual void PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                             double * const * price);
    /** product payoff method at steps earlier than maturity */
    virtual void PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                 double * const * price);


    /** premium scaling */
//    virtual double scalePremium(const double& fairValue);

    /** extra output requests */
    virtual void recordOutputRequests(Control* control, Results* results, double fairValue);

    /** extra output requests */
    virtual bool Positive() { return true; }

    virtual string getCcyTreatment();

    /** make price refinement - control variate */
    virtual double RefinePrice(double basePrice, double discFactor, bool)
        {
            // scale premium for forward starting.
            //double price = discFactor*scalePremium(basePrice);
            double price = basePrice;
            return price;
        }

protected:
 
    int             NumOfPriceArray;

    // num of vol days per day, used for once a day barrier adjustment
    double          VolDaysPerYear;

    vector<double>  Strikes;

    int numKGrid;       // num of Strike Grid
    CDoubleArraySP     Prices;     //to store Results

    /** price dead instrument */
//    bool    CheckDeadInstr();
//    void    CalcDeadInstr(int step, int bot, int top, double * const * price);

private:
    const VanillaGridCEVJ* instrCEVJ;

};

typedef smartPtr<VanillaGridCEVJ> VanillaGridCEVJSP;

DRLIB_END_NAMESPACE
#endif

