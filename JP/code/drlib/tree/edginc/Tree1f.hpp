//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1f.hpp
//
//   Description : one factor trinomial tree base class
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TREE1F_HPP
#define EDG_TREE1F_HPP

#include "edginc/Model.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/EventAssetMove.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/LatticeModelEDR.hpp"

DRLIB_BEGIN_NAMESPACE

//#define TREE_DEBUG

class ILatticeProdEDR;

class TREE_DLL CTree1f: public LatticeModelEDR,
               public virtual FDModel::IFDSolver
{
public:
    friend class CTree1fHelper;
    static CClassConstSP const TYPE;
    static void load(CClassSP& );

    //FDModel::IFDSolver interface
    virtual void roll();
    /** create the tree specific solver */
    virtual IFDSolverSP createSolver() { 
        return IFDSolverSP(this, NullDeleter());
        // return IFDSolverSP::attachToRef(this); 
    }

    // !!! ********** bridge to new interface
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;
    bool    equalTime; // $unregistered
    bool    isCall; // $unregistered
    int     noExerciseWindow; // $unregistered
    int        NumOfPrice; // number of prices in Price array $unregistered
    int        NumOfInsertNode; // number of insert nodes array $unregistered

    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    double getDrift(int step) const{
        double drift = 1.0; // used only for rolling boundary nodes at step < last step
        if (step < getLastStep())
            drift = StepForward[step+1]/StepForward[step]; // step drift

        if (RollDirection==-1)
            drift = 1.0/drift;
        return drift;
    }

    /**  products call this to set discount curve  */
    void setDiscountCurve( const YieldCurveConstSP & yc ) { discYC = yc; }

    /**  rebuild grid or not  */
    virtual bool isRebuilt(const CControl* control);

    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /** accept or reject factors */
    virtual bool acceptFactor( IMarketFactor * factor );

    /** retrieving market data from MDF, new fd/tree interface */    
    virtual void retrieveFactor();

    /**  model initialisation */
    virtual void initModel();

    /** add dividend critical dates */
    void addDivCritDates(const DateTimeArray& divCritDates);

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const;

    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const;

    typedef enum{NO_SMOOTHING=-1, DEFAULT, NODE_INSERTION, SMAX, INTEGRATION} TSmooth;
    typedef enum{U_SHIFT, M_SHIFT, D_SHIFT} TNodeDir;
    typedef enum{MAKE_CACHE, USE_CACHE, RECOMPUTE} TCacheMode;

    CTree1f(CClassConstSP clazz); // for derived classes

    virtual ~CTree1f();

    /** main model entry point */
    virtual void Price(CInstrument*  instrument, CControl* control, CResults* results);

#ifdef TREE_DEBUG
    CIntArraySP     stepCanExercise;
    int             DebugLevel;
#endif


    /** tree initialisation */
    virtual void Setup(const DateTime&      valDate, 
                       const DateTimeArray& segDates, 
                       const vector<int>&   density, 
                       const DateTimeArray* critDates,
                       double               minGap, 
                       bool                 equalTime, 
                       int                  numOfPrice, 
                       int                  numInsertNode,
                       bool                 isCall = true,
                       int                  noExerciseWindow = 0);

    /** price interpolation on a tree, public for use in tree stacks */
    double    TreeInterp(double s, bool useCurrentArray, int iPrice, int j, int *j_out=NULL);
    /** spline interpolation */
    void TreeInterpSpline(double* s, double* price, double* sInterp, 
                          int bot, int top, int iPrice, vector<double>& result);

    /** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
        they must be either one per step or one per node (from -BotClip to TopClip)
        returns true if constant probabilities can be used for all nodes (fast roll) */
    virtual bool CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& dvar_dt,    vector<double>* drift) = 0;
    
    /** determines if cached rebranching info should be used for same grid tweaks */
    virtual bool UseRebranchingCache(){
        return false;
    }
    
    ///*** access functions 
    // set initial tree input steps
    void SetStepsPerYear(int step_per_year);
    // get initial tree input steps
    int    GetStepsPerYear() const;
    // get sweep direction
    int    GetRollDirection() const;
    // get smoothing method chosen 
    TSmooth    GetSmoothMethod() const;
    // set smoothing method chosen 
    void SetSmoothMethod(TSmooth smoothMethod);
    // set smoothing method string
    void SetSmoothString(const string& smoothString);
    // get tree slice index
    int GetSliceIdx() const;
    /** set inserted node levels */
    void SetInsertNode(int idx, int insPt, double insStock, int insPriority, bool isActive = true);
    /** set inserted node level and price */
    void SetInsertNodeAndPrice(int idx, int insPt, double insStock, int insPriority, 
                               int pStart, int pEnd, double insPrice);
    bool getInsertNodePrice(int idx, int iPrice, int insPt, double* insLevel) const;

    /** reset tree node boudary within original limits, used for cases like KO */
    void ResetNodeBoundary(int step, int idx, int bot, int top);
    /** reset price array within original limit - eg. reduce it if some price arrays are no longer need */
    void ResetPriceArr(int step, int idx, int pStart, int pEnd);

    // penultimate smoothing also uses this
    double    TruncationStd;
    bool    TweakingSameGridDelta; // this is the flag to indicate if the price call is doing delta $unregistered
    // DEBUG_xxx are registered
    bool    DEBUG_SameGridDelta; // true => spot is tweaked, but same grid is used
    bool    DEBUG_SameGridVegaRho; // true => same grid for rho and vega 
    bool    DEBUG_UseCtrlVar; // true => same grid for rho and vega 
    //  DEBUG about Node Insert Method
    double  DEBUG_InsertWidth;      //Define the range of the node to be node inserted node at Roll
    int     DEBUG_SmoothWeigth;     //Define the weight of inserted node price at Roll.
    // for dollar div
    int     DEBUG_DollarDivMethod; // to be removed, not used
    int     RG_DEBUG_DollarDivMethod; // no longer used $unregistered
    // for preserving rebranching decisions for same grid tweaks
    bool    UseSameGridGeometry; // true => same rebranching / inserted node info for same grid tweaks
    bool    sameGridVega; // true => same grid (no rebuild tree) for vega tweak $unregistered
    bool    sameGridRho; // true => same grid (no rebuild tree) for vega tweak $unregistered
    bool    sameGridFwdStart; // true => same grid (no rebuild tree) even forward starting case.
    StringArray rebuildList;       // The list of Sensitivity when rebuild is required.  
                                   // 1st component Should be 'Default', 'DefaultPlus', 
                                   // 'DefaultMinus' or 'ListBelow'.
                                   // from 2nd, sensitivities packetName is expected.
    bool    useCcyBasis;       

//#ifdef NEW_TREE_TECH_TEST
    bool    useGammaNodes; // $unregistered
    int     gammaNodesInterval;
    double  gammaThreshold;         // unit is % of notioal is expected as input.
    double  gammaThresholdScaled;   //internal use, scaled by instruments.

    // set up gammaThresholdScaled.  
    void SetGammaThresholdScaled(double refStrike, double notional);

//#define TREE_THETA_CAP
#ifdef  TREE_THETA_CAP
    // activate theta smooth (via instrument).
    bool activateThetaCap(bool isPositiveThetaCap, double ThetaCapThrehold);

    // set up ThetaCapThreshold.  
    void ScaleThetaCapThreshold(double refStrike, double notional);
#endif

    virtual void NodeSetWithAddNode(double* s, int bot, int top, double f0, double spacing);

    // adding nodes.  
    // 1. call calcGammAndGridMap to know where to add new nodes.
    // 2. store the all stocks & prices.
    // 3. calulate new additional nodes and prices (interpolated by t+1).
    // 4. call products->update.
    // 5. re-arranging (combining) new nodes & already exising nodes.
    void addNodesByGamma(const FDProductArray & products, bool isMat = false);

    // calculate gamma on each tree nodes and return addGridMap.  
    // addGridMap will contain index array in original (big) measure, whose are going to be added.  
    // in this function, addMap is also updated.
    bool calcGammaAndGridMap( int step, vector<int> & addGridMap);

//#endif

    /** override a control shift (eg for delta on trees)
        returns new control if constructed else returns 0 */
    virtual SensControl* AlterControl(const SensControl* currSensControl) const;

    /** simple size adjustment for smoothing delta, return true if adjustment made, else return false */
    void DeltaSizeAdjust(double shiftSize, const DateTime& valueDate,
                         const DateTime& exerciseStartDate, const DateTime& matDate, 
                         bool isBermudan, bool hasDollarDivFromValueDateToMat,
                         bool isFwdStarting);

    bool DivsTreatedAsAbsolute() const;

    void SetDivAmountTreatment(bool divTreat);

    double GetStockFloor(int step) const;

    /** decide if a new tree should be built or use the same tree */
    virtual bool isTreeRebuilt(const CControl* control);

    /** overwrite control of sameGrid from instrument, mainly for fwdStarting */
    virtual void controlSameGridFwdStart(const string& ccyTreatment);

protected:
    // cached lattice prod (for performance reasons)
    ILatticeProdEDR * latticeProd; // $unregistered

    FDProductSP payoffIndex; // $unregistered

    // status info
    bool    BuildTree;  // if tree needs to be built in this sweep $unregistered

    // tree parameters, come in from ModelParam
    // registered fields
    int     StepsPerYear;  // num of tree steps per year, actual steps created in TimePts may be different
    double    TreeAlpha; // initial setting for sum of prob_up and prob_down
    bool    DivAmountTreatment; // if dividend treated as absolute amount
    TSmooth SmoothMethod; // not registered, translated from SmoothType string $unregistered

    // tree processing control
    int        RollDirection; // tree rolling direction, -1(default)=roll tree backward, 1=roll tree fwd to calc prob $unregistered
    bool    PerNodeProb; // true if recalculate probability at every node $unregistered

    double    Pu, Pm, Pd; // prob's for the time step (used only by model that gives constant prob's per step like log-normal case $unregistered $unregistered $unregistered

    int        CurrStep; // current tree step being processed $unregistered
    int        CurrIdx; // index for the current time slice, next time slice idx = 1- CurrIdx $unregistered
    int        CurrSeg; // current tree segment being processed, usually =0 $unregistered
    int        TopClip[2]; // top node clipping index at curr step and next step $unregistered
    int        BotClip[2]; // bottom node clipping index at curr step and next step $unregistered
    int        TopDim[2]; // top tree node with standard truncation $unregistered
    int        BotDim[2]; // bottom tree node with standard truncation $unregistered
    int        Ndim; // dimension of regular tree nodes on each side (2*Ndim+7 is total) $unregistered
    int        TopArrIdx[2]; // unpper bound index for price array $unregistered
    int        BotArrIdx[2]; // lower bound index for price array $unregistered

    // underlying node arrays
    double *             Stock[2]; // node level $unregistered
    TreeSliceGeneralSP   underlierPrice; // $unregistered
    double *             UnderlierPrice[2]; // same as underlierPrice $unregistered
    // price arrays, allow any NumOfPrice to be priced in one tree (can be used for state variable tree)
    vector< double * >   NodePrice[2]; // price arrays $unregistered

    TCacheMode CacheMode; // $unregistered
    int  **CacheBranching[4]; // Cache of branching decision made during pricing at each spot, t $unregistered
    bool **CacheInsNodeUsed[2]; // Cache of boolean (for both regular & ins nodes) representing whether an inserted node was used or not during pricing.   $unregistered
    int **CacheInsNodeIdx[2]; // Cache of index of inserted node  $unregistered
    double **CacheInsNodeValue[2]; // Cache of spot value for inserted node $unregistered
    int **CacheIdxUsed; // Cache of index j (index of middle node closest to spot) used for each j_curr  $unregistered
    int *CacheIdxLower; // Cache of lowest j_curr used for a CalcProb call for a given time slice $unregistered
    int *CacheIdxUpper; // Cache of highest j_curr used for a CalcProb call for a given time slice $unregistered
    StringArray usggList;       // The list of Sensitivity when useSameGridGeometry is active.
                                // 1st component Should be 'Default', 'DefaultPlus', 'DefaultMinus' or 'ListBelow'.
                                // from 2nd, sensitivities packetName is expected.

    // fwd stock at each step */
    CDoubleArray StepForward; // $unregistered
    CDoubleArray StockFloor; // $unregistered
    DateTimeArray divCritDates; // $unregistered
    // centre node */
    vector<double> CentreNode; // $unregistered
    // holds node space (in transformed unit) for all tree segments.
    vector<double>  NodeSpace; // $unregistered

    // for node insertion
    TreeSliceGeneral::RangeSP insRange; // inserted stock range $unregistered
    TreeSliceGeneralSP        insStock; // inserted stock levels $unregistered
    TreeSliceGeneralContSP    insNodePrice; // price arrays for inserted node $unregistered
    double *              InsStock[2]; // same as insStock $unregistered
    vector< double * >    InsNodePrice[2]; // same as insNodePrice $unregistered
    int *                 InsPriority[2]; //priority of inserted nodes if more than one within same range $unregistered
    bool                  isActiveInsNode;  //insert Node is active or not.   
    // This "isActiveInsNode" flag is added to treat the following case.
    // When UseCacheAtNode = true, and base price used the Inserted node, but it's not always true that
    // tweaked price calc (e.g. Spot Shifted) also use the Inserted Node.  
    // For example, the barrier could be breached after tweak.  In this case,
    // Inserted node is not used.  When there is no insert node, the insStock have very low or high number.
    // This shouldn't be used, but it was used in the UseCache framework and caused a bug.
    // The flag is used to avoid this.
    // In future, it may be better to use this flag widely, rather than giving the dummy insStock Level.
    // In this case, this should be Array which is corresponding to each inserted node levels!

    /** clean up */
    virtual void    Clear();
    
    /** called before exiting tree Setup */
    virtual void PostSetup();

    /** calculate probability given spot fwd=s1, var_dt wrt s1, current node index j,
        inserted node s_ins(0 if n/a)
        the calculated values are the probabilities and node branching  */
    inline void    CalcProb(double s1, double v_dt, int j, double s_ins, int &u_shift, int &m_shift, 
        int &d_shift, double &pu, double &pd, double &pm, bool tryRebranch);

    /** initialise for forward induction */
    virtual void InitFwdInduction(int priceStart, int priceEnd);

    /** calculate node spacing for each time segment, usually just one segment.*/
    virtual    void CalcNodeSpacing() = 0;

    /** set up tree nodes at the given step */
    virtual void NodeSetup(int step, int idx);
    /**   tree node space set up */
    virtual void NodeSetupSpace(double* s, int bot, int top, 
                                double* sLast, int botLast, int topLast, 
                                double f0, double spacing, int step, double driftNext) =0;

    /** forward or backward roll the tree at one time step */
    virtual void RollTree(vector<double> *drift_out = 0, vector<double> *driftIns_out=0);

    /** roll each node */
    virtual void RollNode(const vector< double * > & arr, double s1, int pStart, int pEnd,
                          double v_dt, double df, int j_curr, bool isInsertNode);
    // called by Roll() after looping through all nodes to interpolate two more nodes for roll
    void InterpTreeBoundary();

    // get segment index of a tree
    int GetTreeSeg (int step);
    /** changing time line segment of a different density */
    bool SwitchSegment();

    // probability stored
    int        NumOfProbMarking; // num of points that fwd prob's are to be stored $unregistered
    double    ***FwdProb; // forward probabilities of reaching a node $unregistered
    
    void convertSmoothTypeString();

    // true if we've computed cache info at this node, and should be using
    inline bool UseCacheAtNode(int j_curr)
    {
        return (CacheMode == USE_CACHE && 
            j_curr >= CacheIdxLower[CurrStep] &&
            j_curr <= CacheIdxUpper[CurrStep]);
    }
    // const refs to asset and discount curve
    CAssetConstSP               Underlier; // $unregistered

//********* some helpers
    /** get processed vol for tree1f */
    virtual void InitVol() = 0;
    /* calculate a term structure of vol^2 
       or simply one point at maturity */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2) = 0;

    // used to init pseudo equity in new interface
    virtual void finaliseModel(CControl*    control);

    //-- for the new tree technique, additional grids on Gamma --//
    // "original" means the ordinary grids geometory with this gamma grids method //
    int NLimit;                     // the original of Ndim.
    vector<int> addGridMap; // $unregistered
                                    // the index array in original Grids
                                    // of adding nodes (addMaps=2 nodes) 
    vector< vector<int> > addGridMap_cache; // for cache
    vector<double>  drift_GN;          // store the drift.

    //-- for the new tree technique, additional grids on Gamma --//
    class GammaNodeMap{
    public:
        GammaNodeMap();

        // initialize orgIdx and clean other int array
        void initIdx(int nSize);

        // store addGridMap, which include the OrgIndex to be added at this time step.
        void setAddGridMap(vector<int> addGridMap);

        void SetAtMat(int NLimit, int gammaInterval);
        
        void Set(const int botOrg, const int topOrg,const GammaNodeMap &gnmPre);

        // store the stocks and prices array.
        void storeStockAndPrice(int oldBot,int oldTop,int NumOfPrice,double* s,double** price);

        // integrate stored stock & prices and additional strikes & prices.
        void integrate(int numAdded, int NumOfPrice,
                        const vector<double>& addStock, 
                        const vector<vector<double> >& addPrice, 
                        double* s,                     //output
                        vector< double * > & price);    //output
        
        // update (set) idxArray.
        void updateIdxAfter();

        int calcNumOfActiveNodes() const;

        int calcLowerGrid(const vector<int> &idxArray) const;

        void turnOnByOrg(int idx);

        const vector<int>& getStatus() const;

        const vector<int>& getOrgIdx() const;

        int getBotBefore()  const;

        int getTopBefore() const;
        
        int getBotAfter() const;

        int getTopAfter() const;

        vector<int> getAddGridMap() const;

        void getOrgIdxOfAddLimit(int &addBotLimit, int &addTopLimit) const;

        // inline
        inline int getOrgIdxByBef(int currIdx) const{
            int iOffset = currIdx + baseBef;
            ASSERT(0<=iOffset && iOffset<(int)idxBefore.size());
            return idxBefore[iOffset];
        };

        inline int getStByBef(int currIdx) const{
            int iOrgOffset = getOrgIdxByBef(currIdx) + baseOrg;
            ASSERT(0<=iOrgOffset && iOrgOffset<(int)orgIdx.size());
            return status[iOrgOffset];
        };

        inline int getStByOrg(int idx) const{
            ASSERT(0<=idx+baseOrg && idx+baseOrg<(int)status.size());
            return status[idx+baseOrg];
        };

    private:
        vector<int> orgIdx;         //original node index.  It's used for the spot level in NodeSetUp
        vector<int> idxBefore;      //node index before additional gamma, but would be truncated.  
        vector<int> idxAfter;       //node index after additional gamma.  
        vector<int> status;         //node status 
                                    // 0 : always used, not additional nodes
                                    // 1 : the additional Nodes, and already used (activated).
                                    // 2 : the additional Nodes, and new additonal on the time slice.  (to be activated).
                                    // 3 : ceiling node (BotClip+1).
                                    // -1: the additional Nodes, but not used.
        vector<int> addGridMap;     // include the OrgIndex to be added at this time step.
        vector<double> S_store;
        vector<vector<double> > NP_store;
        int baseOrg;           // number of lower grid of OrigGrid.
        int baseBef;         // number of idxBefore
        int baseAft;        //number of lower grid of idxAfter
        typedef enum {NOT_USED_IDX = -999999} IdxStatus;
    };

    vector<GammaNodeMap> gammaNodeMaps;

    //----//

private:

    void SetDefault();
    string  SmoothMethodString; // registered string for smooth type

    double  TreeDeltaShift; // $unregistered

    friend class Tree1fSolver;

    void validateList();

    bool isActiveUSGG(const CControl* control);

    bool hasSensInList(const SensitivitySP sens, const StringArray *strgList);

    bool rebuildRequest(const SensitivitySP sens, bool  divAmountTreatment);

    // special class for smoothing theta 
    double  ThetaCapThrehold;
    double  threholdTC;         // internal scaled threshold
    bool    isPositiveThetaCap;
    bool    useThetaCap;
    vector <double> thetaArray;

    // overwrite price array
    void addThetaCap();

};
typedef smartPtr<CTree1f> CTree1fSP;

DRLIB_END_NAMESPACE
#endif
