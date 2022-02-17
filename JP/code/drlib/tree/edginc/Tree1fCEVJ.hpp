//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fCEVJ.hpp
//
//   Description : one factor trinomial tree for CEV+Jump process.
//
//----------------------------------------------------------------------------

#ifndef EDG_Tree1fCEVJ_HPP
#define EDG_Tree1fCEVJ_HPP

#include "edginc/Tree1f.hpp"
#include "edginc/CEVJProcessed.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/Model1F.hpp"
#include "edginc/Tree1fPayoff.hpp"

DRLIB_BEGIN_NAMESPACE

/** CEVJ tree */
class TREE_DLL CTree1fCEVJ : public CTree1f
{
public:

    /** old interface, just leave it to work */
    class IProduct : virtual public Model1F::Product1F,
                     public TreePayoff1F {
    public:

        CTree1f * tree1f;
        /** initialise tree1f, only be needed once for base price call */
        virtual void InitTree(CControl* control){
            Init(control);
        }

        /** initialising and setting product variables */
        // this is called per pricing call before tree sweep call (after InitTree)
        virtual void InitProdTree(){
            InitProd();
        }
        /** a chance for product to change inserted node before each roll() 
        default implementation returns false meaning that no node moving needed */
        virtual bool MoveInsNode(int idx, int iPrice){ return false;}

        IProduct() : tree1f(0){}
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        IIntoProduct(){};
        virtual ~IIntoProduct(){};
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(CTree1fCEVJ* model) const = 0;
    };

    static CClassConstSP const TYPE;
    static void load(CClassSP&);
    
    CTree1fCEVJ();
    virtual ~CTree1fCEVJ();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    int DEBUG_BSSmooth;  //BSSmoothing switch for debug use.
    int DEBUG_JumpFrom;  //JumpFrom switch for debug use.

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, even though at some level we would like
     * to do risk mapping on CEVJ if we were using it, because (a) we're not,
     * (b) CEVJProcessed is annoyingly not a proper reflection class.  See
     * IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

protected:
    CTree1fCEVJ(CClassConstSP const type);

    /** set up variance array */
    virtual    void PostSetup();

    /** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
        they must be either one per step or one per node (from -BotClip to TopClip)
        returns true if constant probabilities can be used for all nodes (fast roll) */
    virtual bool CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& dvar_dt,    vector<double>* drift);
    /** calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num */
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);
    /** calculate node spacing for each time segment, usually just one segment */
    virtual    void CalcNodeSpacing();
    /** set up tree nodes at the given step */
    virtual void NodeSetup(int step, int idx);
    /**   tree node space set up */
    virtual void NodeSetupSpace(double* s, int bot, int top, 
                  double* sLast, int botLast, int topLast, 
                  double f0, double spacing, int step, double driftNext);

    // total variance from step 0
    CDoubleArray Variance; // $unregistered
    CLatticeDouble* Spots;  //Stock Price on Tree Node $unregistered

    double    **JumpPrice[2]; // price arrays $unregistered
    double  **InsJumpPrice[2]; // inserted Jump Price Array $unregistered
    vector<double>  JumpUpNode;    // To store all Jump Up Node on Stock[][] $unregistered
    vector<double>  JumpDownNode;  // To store all Jump Down Node on Stock[][] $unregistered

    double  StepJumpRate;   // JumpRate at CurrStep $unregistered

//********* some helpers
    /** get processed vol for tree1f */
    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                        const DateTime& matDate, CTermStructure& v2);

    double    CalcNodeLevel(double x, double s_last, double time);

    void CalcJumpPrice(const vector<double>& d_node, 
                                const vector<double>& u_node,
                                const int bot,const int top, 
                                const bool isInsertNode = false, 
                                const int pPrice = 0);

    double CalcDrift(double x, double v_dt, double asset_drift, double t);

    virtual void RollTree(vector<double> *drift_out = 0, vector<double> *driftIns_out=0);

    virtual void Clear();

    virtual bool isTreeRebuilt(const CControl* control);

private:
    friend class CTree1fCEVJCalib;

    virtual void SetDefault();
    bool IsStoreJumpNode; // $unregistered

   // Interpolated data along to Terem Line.
    CTermStructure termDiffV2; // diff vol^2 term structure, scaled with s_ref^(1-CEVpower) $unregistered
    CTermStructure termCEVPower; // $unregistered
    CTermStructure termJumpRate; // $unregistered
    CTermStructure termJumpMean; // $unregistered
    CTermStructure termJumpWidth; // $unregistered

    //  For Store purpose at SetUpTree.
    CDateTime ValDate;           // $unregistered
    CDateTime StartDate; // $unregistered
    CDateTime MatDate; // $unregistered

    // CEV+Jump vol
    CEVJProcessedSP VolCEVJ; // $unregistered

    // for quanto and struck use
    FXAssetConstSP fxAsset; // $unregistered
    const Correlation* eqFXCorr; // $unregistered
    CVolProcessedBSSP volFXBS; // $unregistered
    FXVolBaseWrapper fxVol;  // $unregistered
};

DRLIB_END_NAMESPACE
#endif
