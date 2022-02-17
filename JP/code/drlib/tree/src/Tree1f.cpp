//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1f.cpp
//
//   Description : one factor trinomial tree base class.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/FXCrossGamma.hpp"
#include "edginc/PseudoSimpleEquity.hpp"
#include "edginc/RollingTheta.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/VegaPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/IndexSpecEDR.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

static int numOfInterp = 0;

const double DefaultTruncation = 5.0; // 5 stdev truncation

/* XXX Not sure how to provide central definition of these */
#define TWK_DEFAULT "Default"
#define TWK_DEFAULT_PLUS "DefaultPlus"
#define TWK_DEFAULT_MINUS "DefaultMinus"
#define TWK_LIST_BELOW  "ListBelow"

//  tried following code to check how slow exp().
//  using VC++7.0, the results of both shows 22 sec!! (No diffs!!)
//#include   <time.h>
//#include   <stdio.h>
//
//int testtime()
//{
//    clock_t  start, end;
//
//
//    double x=0;
//    start = clock();
//    for (double i=0.5; i<1.5;i+=1.0/1e10){
//        x = exp(i);
//    }
//    end = clock();
//    cerr    <<  "calc time " << (end - start)/CLOCKS_PER_SEC <<" [sec] \n";
//
//    start = clock();
//    for (double i=0.5; i<1.5;i+=1.0/1e10){
//        x = i*i;
//    }
//    end = clock();
//    cerr    <<  "calc time " << (end - start)/CLOCKS_PER_SEC <<" [sec] \n";
//
//    return (0);
//}

/************ for debuging tree step outputs **************/
//#define TREE_DEBUG_FILE "c:\\temp\\mytree"
//#define TREE_GRID_FILE "c:\\temp\\mygrid"
//#define TREE_DEBUG_VOL_FILE "c:\\temp\\myvol"

#ifdef TREE_DEBUG_VOL_FILE

void Debug_OutPutVol(int step, double time, vector<double> var, double dt, double* s, int bot, int top, bool first)
{
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file

    string fileName = TREE_DEBUG_VOL_FILE;
    sprintf(fname, "%s.dat", fileName.c_str());
    if (first)
        DumpFile = fopen(fname, "w");
    else
        DumpFile = fopen(fname, "a+");
    if (DumpFile)
    {
        sprintf(buffer, "#step = %d, TradeTime, stock , Vol \n", step);
        fprintf(DumpFile, buffer);
        int vsize = var.size();
        for (int j_debug=0; j_debug < vsize; j_debug ++)
        {
            double vol = sqrt(var[j_debug]/dt);
            sprintf (buffer, "%f, %f, %f # %d \n", time, s[j_debug+bot], vol, j_debug);
            fprintf(DumpFile, buffer);
        }
        sprintf (buffer, "\n");
        fprintf (DumpFile, buffer);
        fclose(DumpFile);
        DumpFile = 0;
    }
}
#endif

#ifdef TREE_DEBUG_FILE

void Debug_OutPutTree(int step, double time, double *s, vector <double *> price, int bot, int top, int nPrice, const string& fileName, bool first,
                      int insSize, double *s_ins, vector <double *> insPrice)
{
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file

    for (int i=0; i<nPrice; i++)
    {
        sprintf(fname, "%s%d.dat", fileName.c_str(), i);
        if (first)
            DumpFile = fopen(fname, "w");
        else
            DumpFile = fopen(fname, "a+");
        if (DumpFile)
        {
            sprintf(buffer, "#step = %d, TradeTime, stock , Price \n", step);
            fprintf(DumpFile, buffer);
            int j_ins = 0;
            for (int j_debug=bot; j_debug<=top; j_debug ++)
            {
                sprintf (buffer, "%f, %f, %f # %d \n", time, s[j_debug], price[i][j_debug], j_debug);
                fprintf(DumpFile, buffer);
                if (j_ins <insSize && j_debug<top){
                    if (s[j_debug]<s_ins[j_ins] && s_ins[j_ins]<s[j_debug+1]){
                        sprintf (buffer, "%f, %f, %f # %d \n", time, s_ins[j_ins], insPrice[i][j_ins], 0);
                        fprintf(DumpFile, buffer);
                        j_ins ++;
                    }
                }            
            }
            sprintf (buffer, "\n");
            fprintf (DumpFile, buffer);
            fclose(DumpFile);
            DumpFile = 0;
        }
    }
}


void Debug_OutPutTreeProbs(double price, double s_next, double v_dt, double s_ins, double mean, double var, double skew, double pu, double pd, double pm, double df, int j_curr, int d_shift, int m_shift, int u_shift, bool priceInterplated, const string& fileName, bool top)
{
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file
    static bool first = true;
    int maxStepWantToSee = 70;
    static bool skip = true;
    if (!top)
    {
        if (!(j_curr > maxStepWantToSee))
        {
            skip = false;
        }
        else
        {
            skip = true;
        }
    }
    if (skip) return;

    sprintf(fname, "%s%s.dat", fileName.c_str(), "-probs");
    if (first)
    {
        DumpFile = fopen(fname, "w");
    }
    else
        DumpFile = fopen(fname, "a+");
    if (DumpFile)
    {
        if (!top)
        {
            sprintf (buffer, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d\n", price, s_next, s_ins, v_dt, mean, var, skew, pu, pm, pd, df, j_curr, d_shift, m_shift, u_shift, priceInterplated);
            fprintf(DumpFile, buffer);
            fclose(DumpFile);
            DumpFile = 0;

        }
        else
        {
            //hack: j_curr contains step
            if (!first)
            {
                sprintf (buffer, "\n");
                fprintf(DumpFile, buffer);
            }
            sprintf (buffer, "#step = %d, price  pu  pm  pd  df  j_curr  d_shift m_shift u_shift.\n", j_curr);
            fprintf(DumpFile, buffer);
                fclose(DumpFile);
                DumpFile = 0;

            }
            first = false;
        }
}
#endif

#ifdef TREE_GRID_FILE
void Debug_OutGrid(double *s, double *price, int bot, int top, const string& fileName, bool first, const string& sens)
{
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file

    sprintf(fname, "%s.dat", fileName.c_str());
    if (first)
        DumpFile = fopen(fname, "w");
    else
        DumpFile = fopen(fname, "a+");
    if (DumpFile)
    {
        sprintf(buffer, "#stock , Price @ %s \n",sens.c_str());
        fprintf(DumpFile, buffer);
        for (int j_debug=bot; j_debug<=top; j_debug ++)
        {
            sprintf (buffer, "%f, %f \n", s[j_debug], price[j_debug]);
            fprintf(DumpFile, buffer);
        }
        sprintf (buffer, "\n");
        fprintf (DumpFile, buffer);
        fclose(DumpFile);
        DumpFile = 0;
    }
}
#endif

/////////// special functions for smoothing theta /////////////
#ifdef  TREE_THETA_CAP
// overwrite price array
void CTree1f::addThetaCap(){
    if(useThetaCap){
        const TreeSlice & slice = payoffIndex->getValue( CurrStep );
        int bot, top;
        double nextPrice;
        slice.getCalcRange( bot, top );
        int sizeOfSlice =  top-bot+1;
        thetaArray.resize(sizeOfSlice, 0.0);
        ASSERT(CurrStep < timeLine->NumOfStep);
        double dt = timeLine->TradeTime[CurrStep+1]-timeLine->TradeTime[CurrStep];
        double tScale = dt*252.0;   // scaling the theta for 1 day (1/252 year).
        double thetaLevel = threholdTC * tScale;
        for (int i=bot+2; i<=top-2; i++){   // do not smooth for boundary area.
            nextPrice = TreeInterp(Stock[CurrIdx][i], false, 0, i);
            thetaArray[i-bot] = nextPrice - NodePrice[CurrIdx][0][i];
            if (thetaArray[i-bot] > thetaLevel){
                NodePrice[CurrIdx][0][i] += thetaArray[i-bot] - thetaLevel;
            }
        }
    }
}
#endif

/************ for debuging tree step outputs **************/

static bool checkTierTweak(const SensitivitySP sens)
{
    CClassConstSP clazz = CClass::forName("Tier");
    bool tierTweak = clazz->isInstance(sens.get());

    return tierTweak;
}

bool checkVegaTweak(const SensitivitySP sens)
{
    bool vegaTweak = VegaParallel::TYPE->isInstance(sens.get())
                || VegaPointwise::TYPE->isInstance(sens.get())
                || DDeltaDVol::TYPE->isInstance(sens.get())
                || VegaMatrix::TYPE->isInstance(sens.get());

    return vegaTweak;
}

bool checkDeltaTweak(const SensitivitySP sens)
{
    bool deltaTweak = Delta::TYPE->isInstance(sens.get())
                   || BasketDelta::TYPE->isInstance(sens.get())
                   || DeltaDDE::TYPE->isInstance(sens.get())
                   || DDeltaDVol::TYPE->isInstance(sens.get())
                   || CrossGamma::TYPE->isInstance(sens.get())
                   || FXCrossGamma::TYPE->isInstance(sens.get());

    return deltaTweak;
}

////////////////////////// CTree1f /////////////////////////

bool CTree1f::rebuildRequest(const SensitivitySP sens,
                             bool  divAmountTreatment)
{
    bool rebuild = Theta::TYPE->isInstance(sens.get())
                || MuParallel::TYPE->isInstance(sens.get())
                || MuPointwise::TYPE->isInstance(sens.get())
                || MuSpecial::TYPE->isInstance(sens.get())
                || FXDelta::TYPE->isInstance(sens.get())
                || DeltaNextDay::TYPE->isInstance(sens.get())
                || RollingTheta::TYPE->isInstance(sens.get());

    if (divAmountTreatment){
        rebuild = rebuild
                    || RhoParallel::TYPE->isInstance(sens.get())
                    || RhoPointwise::TYPE->isInstance(sens.get())
                    || RhoBorrowParallel::TYPE->isInstance(sens.get())
                    || RhoBorrowPointwise::TYPE->isInstance(sens.get());
    }

    return rebuild;
}

/** decide if a new tree should be built or use the same tree */
bool CTree1f::isTreeRebuilt(const CControl* control)
{
    bool buildTree = (!control || control->isPricing());
    if (!buildTree)
    {
        SensitivitySP sens = control->getCurrentSensitivity();
        if (checkDeltaTweak(sens))
        {
            buildTree = !DEBUG_SameGridDelta;
            TweakingSameGridDelta = DEBUG_SameGridDelta;
        }
        else if (checkVegaTweak(sens)){
            buildTree = !sameGridVega;
        }
        else if (checkTierTweak(sens)){
            buildTree = true;
        }
         else
        {// if sens is not vega or delta, it use Rho flag...  Maybe we want to change here.
            buildTree = !sameGridRho || rebuildRequest(sens, DivAmountTreatment);
        }
        // over control the rebuildTree by using rebuildList
        if (!rebuildList.empty()){
            if (hasSensInList(sens, &rebuildList)){                    
                if (rebuildList[0] == TWK_DEFAULT_PLUS || rebuildList[0] == TWK_LIST_BELOW){
                    buildTree  = true;
                }else if (rebuildList[0] == TWK_DEFAULT_MINUS){
                    buildTree  = false;
                }
            }
        }
    }
    return buildTree;
}

/** overwrite control of sameGrid from instrument, mainly for fwdStarting */
void CTree1f::controlSameGridFwdStart(const string& ccyTreatment)
{
    DEBUG_SameGridDelta = false;
    //DEBUG_SameGridVegaRho = false;
    sameGridVega = sameGridFwdStart;
    sameGridRho = false;                // we haven't tested, so false.
}


// validate the string array
void CTree1f::validateList()
{
    static const string method = "CTree1f::validateList";
    try{        
        for (int i=0; i<2; i++){
            StringArray stringList;
            string currentList;
            switch (i){
                case 0: 
                    if (!usggList.empty()){
                        stringList = usggList;
                        currentList = "usggList";
                        if (!UseSameGridGeometry && stringList[0] != TWK_DEFAULT)
                            throw ModelException(method, "UseSameGridGeometry = false, then stringList should be empty or 'Default'.");
                    }
                    break;
                case 1:
                    if (!rebuildList.empty()){
                        stringList = rebuildList;
                        currentList = "rebuildList";
                    }
                    break;
                default :
                    throw ModelException(method, "internal error.  Shouldn't be here");
            }
            if (!stringList.empty()){
                if (stringList[0] != TWK_DEFAULT
                    && stringList[0] != TWK_DEFAULT_PLUS
                    && stringList[0] != TWK_DEFAULT_MINUS
                    && stringList[0] != TWK_LIST_BELOW)
                {
                    throw ModelException(method, currentList + "'s first string is '"+stringList[0]+"'.\n"
                                        "Should be 'Default', 'DefaultPlus', 'DefaultMinus' or 'ListBelow'");
                }
                if (stringList.size() == 1 && stringList[0] != TWK_DEFAULT)
                    throw ModelException(method, currentList+"[0]='" + stringList[0]
                                                +"', so you need give sensitivities in list.");
                if (stringList.size() > 1 && stringList[0] == TWK_DEFAULT)
                    throw ModelException(method, currentList+"[0]='" + stringList[0]
                                                +"', so The rest should be empty.");            
            }
        }
    }
    catch(exception& e){
        throw ModelException(&e, method);
    }
}


// return the control found in stringArray or not.
bool CTree1f::hasSensInList(const SensitivitySP sens,const StringArray *strgList)
{
    static const string method = "CTree1f::hasSensInList";
    try{
        bool foundInList = false;
        string sensName = sens->getPacketName();
        for (int i=1;i<strgList->size();i++){
            if ((*strgList)[i] == sensName){
                foundInList = true; 
                break;
            }
        }
        return foundInList;
    }catch(exception& e){
        throw ModelException(&e, method);
    }
}

bool CTree1f::isActiveUSGG(const CControl* control)
{
    static const string method = "CTree1f::isActiveUSGG";
    try{        
        if (!UseSameGridGeometry)
            return false;

        bool isActive = false;
        if (usggList.size() == 0 || usggList[0] != TWK_LIST_BELOW){
            isActive = !isTreeRebuilt(control);
        }
        if (usggList.size() == 0 || usggList[0] == TWK_DEFAULT)
            return isActive;    //exit 

        if (usggList[0] == TWK_DEFAULT_PLUS && !isActive){
            throw ModelException(method, 
                                "You choose usggList[0] = 'DefaultPlus', but currently nothing you can add to default");
        }
        
        // search the string in List
        SensitivitySP sens = control->getCurrentSensitivity();
        bool foundInList = hasSensInList(sens, &usggList);
        if (usggList[0] == TWK_DEFAULT_MINUS && isActive && foundInList){
            isActive = false;
        }
        if (usggList[0] == TWK_LIST_BELOW && foundInList){
            isActive = !isTreeRebuilt(control);
        }
        return isActive;
    }
    catch(exception& e){
        throw ModelException(&e, method);
    }
}


/** main model entry point */
void CTree1f::Price(CInstrument* instrument,
                    CControl*    control,
                    CResults*    results)
{
    static const string method = "CTree1f::Price";

    // convert smooth type string to enum for state variable tree.  (only need at the first time)
    // SmoothType could be override by product, and the object could be used in later sens calc.
    // So only just need to at once. I don't think it's best to do at here.....
    if (control && control->isPricing())
        convertSmoothTypeString();

    // init old tree
    // rebuild tree nodes when required
    TweakingSameGridDelta = false;
    BuildTree = isTreeRebuilt(control);

    // if using same tree, use cached rebranching/ins node info to maintain same geometry
    CacheMode = RECOMPUTE;
    UseSameGridGeometry &= UseRebranchingCache();
    if (UseSameGridGeometry) // currently only Local Vol
    {
        if (control && control->isPricing())
        {
            CacheMode = MAKE_CACHE;
        }
        else if(control && isActiveUSGG(control))
        {
            CacheMode = USE_CACHE;
        }
    }

    CurrIdx = -999;  //just to check the inst is already dead or not.  When it's dead, CurrIdx is not set.  n
                        // of couse, need to review.
    FDModel::Price(instrument, control, results);
    if (control && control->isPricing() && CurrIdx != -999){
        if(DEBUG_SameGridDelta)
        {
            double deltaSize = control->getDeltaShiftSize();
//            if (Maths::isZero(deltaSize))
//                deltaSize = Delta::DEFAULT_SHIFT;
            TreeDeltaShift = 10.0*deltaSize*(Stock[CurrIdx][1] - Stock[CurrIdx][-1])/Stock[CurrIdx][0];

            // store some extra results
            results->storeScalarGreek(numOfInterp, Results::DEBUG_PACKET,
                                      OutputNameSP(new OutputName("TREE_NUM_INTERP")));
        }
    }
        //debug output (not RefinePrice is applied....)
#ifdef TREE_GRID_FILE
        string tmpSt;
        bool isPrice = false;
        if (control) {
            CDoubleMatrixSP stockAndPrice(new DoubleMatrix(TopClip[CurrIdx]+BotClip[CurrIdx]+7, 2));
            for (int iNode=-BotClip[CurrIdx]-3; iNode <=TopClip[CurrIdx]+3; iNode ++){
                (*stockAndPrice)[iNode+BotClip[CurrIdx]+3][0] = Stock[CurrIdx][iNode];
                (*stockAndPrice)[iNode+BotClip[CurrIdx]+3][1] = NodePrice[CurrIdx][0][iNode];
            }
            if (control->isPricing()){
                tmpSt = "price";
                isPrice = true;
            }else{
                SensitivitySP tmpSens = control->getCurrentSensitivity();
                tmpSt = tmpSens->getPacketName();
            }
            OutputNameSP tmpName = OutputNameSP(new OutputName(tmpSt));
            results->storeGreek(stockAndPrice, "DEBUG", tmpName);
        }
        Debug_OutGrid(Stock[CurrIdx], NodePrice[CurrIdx][0], -BotClip[CurrIdx]-3, TopClip[CurrIdx]+3, TREE_GRID_FILE, isPrice, tmpSt);
#endif
#ifdef TREE_DEBUG
        if (DebugLevel > 0)
        {
            OutputNameConstSP stepDatesName(new OutputName("StepDates"));
            OutputNameConstSP exerStepsName(new OutputName("ExerciseSteps"));

            results->storeGreek(IObjectSP(timeLine->StepDates.clone()), "DEBUG", stepDatesName);
            results->storeGreek(stepCanExercise, "DEBUG", exerStepsName);

            // copy results to original tree
            this->PriceEnd = treeToUse->PriceEnd;
            this->DeltaEnd = treeToUse->DeltaEnd;
            this->GammaEnd = treeToUse->GammaEnd;

        }
#endif
}

void CTree1f::load(CClassSP& clazz)
{
        REGISTER(CTree1f, clazz);
        SUPERCLASS(FDModel);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        FIELD(StepsPerYear, "num of tree steps per year");
        FIELD(TruncationStd, "num of stdev to truncate");
        FIELD(TreeAlpha, "initial setting for sum of prob up and prob down");
        FIELD(DivAmountTreatment, "true if use absolute amount dividend treatment in tree");
        FIELD_MAKE_OPTIONAL(DivAmountTreatment);
        FIELD(SmoothMethodString, "smoothing type: NO_SMOOTHING, DEFAULT, NODE_INSERTION, SMAX, INTEGRATION");
        FIELD_MAKE_OPTIONAL(SmoothMethodString);
        FIELD(DEBUG_DollarDivMethod, "For DR use only : 0=old tree interface, 1=new tree/fd interface");
        FIELD_MAKE_OPTIONAL(DEBUG_DollarDivMethod);

        FIELD(DEBUG_SameGridVegaRho, "For DR use only : true(default)=same grid tweak for vega and rho");
        FIELD_MAKE_OPTIONAL(DEBUG_SameGridVegaRho);
        FIELD(DEBUG_SameGridDelta, "For DR use only : true(default)=same grid tweak for delta/gamma");
        FIELD_MAKE_OPTIONAL(DEBUG_SameGridDelta);
        FIELD(DEBUG_UseCtrlVar, "For DR use only : true(default)=use control variate when possible");
        FIELD_MAKE_OPTIONAL(DEBUG_UseCtrlVar);

        FIELD(DEBUG_InsertWidth, "negative(-1:default) is same before.  positive dabule(e.g. 5%) is Insert Area Width")
        FIELD_MAKE_OPTIONAL(DEBUG_InsertWidth);
        FIELD(DEBUG_SmoothWeigth, "0 (default) no use, 1 is used.");
        FIELD_MAKE_OPTIONAL(DEBUG_SmoothWeigth);

        FIELD(UseSameGridGeometry, "true: use same rebranching and inserted node info as base price for same grid tweaks");
        FIELD_MAKE_OPTIONAL(UseSameGridGeometry);
        FIELD(sameGridFwdStart, "true: try to use same Grid Tree even Fwd Starting instruments.");
        FIELD_MAKE_OPTIONAL(sameGridFwdStart);
        FIELD(useCcyBasis, "use currency basis (default = false)");
        FIELD_MAKE_OPTIONAL(useCcyBasis);

        FIELD(gammaNodesInterval, "interval for gamma nodes.  0 or 1 is default (no use gamma nodes).");
        FIELD_MAKE_OPTIONAL(gammaNodesInterval);
        FIELD(gammaThreshold, "Gamma Threshhold in % of Notional, only used for gammaNodesInterval > 1.");
        FIELD_MAKE_OPTIONAL(gammaThreshold);

        FIELD(rebuildList, "The list of Sensitivity, in which rebuild Grid is required.");
        FIELD_MAKE_OPTIONAL(rebuildList);
        FIELD(usggList, "The list of Sensitivity, in which UseSameGridGeometry is active.");
        FIELD_MAKE_OPTIONAL(usggList);

#ifdef  TREE_THETA_CAP
        FIELD(useThetaCap, "flag for internal use.");
        FIELD_MAKE_TRANSIENT(useThetaCap);
        // others are not necessary to be TRANSIENT, because it's always set at initProd().
#endif

}

// helpers
CClassConstSP const CTree1f::TYPE = CClass::registerClassLoadMethod(
    "Tree1f", typeid(CTree1f), CTree1f::load);

CTree1f::CTree1f(CClassConstSP clazz) : LatticeModelEDR(clazz)
{
    useCcyBasis = false;
    StepsPerYear = 0;
    TruncationStd = DefaultTruncation;
    SmoothMethod = DEFAULT;
    SmoothMethodString = "DEFAULT";
    isActiveInsNode = true;
    CacheMode = RECOMPUTE;
    TreeAlpha = 1.0/3.0;
    DivAmountTreatment = false;
    Ndim = 0;
    equalTime = false;
    NumOfInsertNode = 0;
    NumOfPrice = 1;
    isCall = true;
    noExerciseWindow = 0;
    NumOfProbMarking = 0;
    RollDirection = -1;
    TreeDeltaShift = 0.0;
    DEBUG_DollarDivMethod = 1;
    Pu=Pm=Pd = 0.0; // stop UMR's for non LN tree read.
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = true;
    DEBUG_UseCtrlVar = true;

    DEBUG_InsertWidth = -1.0;
    DEBUG_SmoothWeigth = 0;

    TweakingSameGridDelta = false; // this must be false here !!!

    UseSameGridGeometry = false;  // currently smooths greeks but impacts performance
    sameGridVega  = true;   // same grid vega works generally, except LN quanto FwdStart case.
    sameGridRho  = true;   // same grid rho works generally, except Fwd Start Case.
    sameGridFwdStart = false; //same grid is usually not available when it's fwd starting.

    usggList.resize(1);
    usggList[0] = TWK_DEFAULT;

    NLimit = Ndim;
    addGridMap.clear();
    drift_GN.clear();
    useGammaNodes = false;
    gammaNodesInterval = 0;
    gammaThreshold = 0.;       
    gammaThresholdScaled = 0.;
    gammaNodeMaps.clear();

#ifdef  TREE_THETA_CAP
    ThetaCapThrehold = threholdTC = 0.0;
    isPositiveThetaCap = useThetaCap = false;
#endif

#ifdef TREE_DEBUG
    DebugLevel = 0;
#endif
    SetDefault();
}

CTree1f::~CTree1f()
{
    Clear();
}

/** clean up */
void CTree1f::Clear()
{
    int i, j;
    // cached rebranching/ins node info
    if (UseSameGridGeometry && CacheBranching[0] != 0)
    {
        int totalStep = timeLine->NumOfStep;        
        for (i=0; i<3; i++)
        {
            for (j=0; j<totalStep + 1; j++)
                delete [] (CacheBranching[i][j] -= (NLimit+3));
            delete [] CacheBranching[i];
        }

        for (i=0; i<2; i++)
        {
            for (j=0; j<totalStep + 1; j++)
            {
                delete [] (CacheInsNodeUsed[i][j] -= (NLimit+3));
                delete [] (CacheInsNodeIdx[i][j] -= (NLimit+3));
                delete [] (CacheInsNodeValue[i][j] -= (NLimit+3));
            }
                delete [] CacheInsNodeUsed[i];
                delete [] CacheInsNodeIdx[i];
                delete [] CacheInsNodeValue[i];
        }

        for (j=0; j<totalStep + 1; j++)
        {
            delete [] (CacheIdxUsed[j] -= (NLimit+3));

        }
        delete [] CacheIdxUsed;
        delete [] CacheIdxLower;
        delete [] CacheIdxUpper;
    }

    if (FwdProb != 0)
    {
        for (i=0; i<NumOfProbMarking; i++)
        {
            for (j=0; j<NumOfPrice; j++)
                delete [] (FwdProb[i][j] -= (NLimit+3));
            delete [] FwdProb[i];
        }
        delete [] FwdProb;
    }

    StepForward.clear();
    CentreNode.clear();
    NodeSpace.clear();

//#ifdef  NEW_TREE_TECH_TEST
    addGridMap.clear();
//#endif

    SetDefault();
}

void CTree1f::SetDefault()
{
    BuildTree = true;
    PerNodeProb = false;

    InsPriority[0] = InsPriority[1] = 0;

    int i = 0;
    for (i=0; i<4; i++)
    {
        CacheBranching[i] = 0;
    }

    for (i=0; i<2; i++)
    {
        CacheInsNodeUsed[i] = 0;
        CacheInsNodeIdx[i] = 0;
        CacheInsNodeValue[i] = 0;
    }

    FwdProb = 0;
}

bool CTree1f::DivsTreatedAsAbsolute() const{
    return DivAmountTreatment;
}

void CTree1f::SetDivAmountTreatment(bool divTreat){
    DivAmountTreatment = divTreat;
}

double CTree1f::GetStockFloor(int step) const{
    if (StockFloor.size() > 0){
        return StockFloor[step];    // no bound checking
    }
    return 0.0;
}

// convert input string to enum type
void CTree1f::convertSmoothTypeString()
{
    if (SmoothMethodString == "NO_SMOOTHING")
        SmoothMethod = NO_SMOOTHING;
    else if (SmoothMethodString == "DEFAULT")
        SmoothMethod = DEFAULT;
    else if (SmoothMethodString == "NODE_INSERTION")
        SmoothMethod = NODE_INSERTION;
    else if (SmoothMethodString == "SMAX")
        SmoothMethod = SMAX;
    else if (SmoothMethodString == "INTEGRATION")
        SmoothMethod = INTEGRATION;
    else
        throw ModelException("CTree1f::convertSmoothTypeString",
                             "unknown smooth type ("+SmoothMethodString+").\n"
                             "Type NO_SMOOTHING, DEFAULT, NODE_INSERTION, SMAX or INTEGRATION needed");
}

// ******* access functions ******
// set and return initial tree input steps
void CTree1f::SetStepsPerYear(int step_per_year)
{
    StepsPerYear = step_per_year;
}
// get and return initial tree input steps
int CTree1f::GetStepsPerYear() const
{
    return StepsPerYear;
}
// get sweep direction
int CTree1f::GetRollDirection() const
{
    return RollDirection;
}
// get current slice index
int CTree1f::GetSliceIdx() const
{
    return CurrIdx;
}
// get smoothing method chosen
CTree1f::TSmooth CTree1f::GetSmoothMethod() const
{
    return SmoothMethod;
}
// set smoothing method chosen
void CTree1f::SetSmoothMethod(TSmooth smoothMethod)
{
    SmoothMethod = smoothMethod;
}

// set smoothing method string
void CTree1f::SetSmoothString(const string& smoothString)
{
    SmoothMethodString = smoothString;
}

// called before exiting Setup
void CTree1f::PostSetup()
{
}

/** set inserted node level */
void CTree1f::SetInsertNode(int idx,int insPt,double insStock, int insPriority, bool isActive)
{
    if (insPt >= NumOfInsertNode) {
        throw ModelException("CTree1f::SetInsertNode",
                             "trying to set node beyond bound");
    }

    InsStock[idx][insPt] = insStock;
    InsPriority[idx][insPt] = insPriority;
    isActiveInsNode = isActive;
}

/** set inserted node level and price */
void CTree1f::SetInsertNodeAndPrice(int idx, int insPt, double insStock, int insPriority, int pStart, int pEnd, double insPrice)
{
    if (insPt >= NumOfInsertNode) {
        throw ModelException("CTree1f::SetInsertNodeAndPrice",
                             "trying to set node beyond bound");
    }

    InsStock[idx][insPt] = insStock;
    InsPriority[idx][insPt] = insPriority;
    for (int i=pStart; i<=pEnd; i++) {
        InsNodePrice[idx][i][insPt] = insPrice;
    }
}

bool CTree1f::getInsertNodePrice(int idx, int iPrice, int insPt, double* insLevel) const{
    if (NumOfInsertNode <1)
        return false;
    else{
        *insLevel = InsNodePrice[idx][iPrice][insPt];
        return true;
    }
}

/****** tree initialisation */
void CTree1f::Setup(const DateTime&      valDate,
                    const DateTimeArray& segDates,
                    const vector<int>&   density,
                    const DateTimeArray* critDatesIn,
                    double               minGap,
                    bool                 equalTime,
                    int                  numOfPrice,
                    int                  numInsertNode,
                    bool                 isCall,
                    int                  noExerciseWindow)
{
    static const string method = "CTree1f::Setup";
    try
    {
    int i, j;

        Clear();

        DateTimeArray critDates;

        // assign user supplied critical dates
        if (critDatesIn)
        {
            for (i=0; i<critDatesIn->size(); i++)
                critDates.push_back((*critDatesIn)[i]);
        }

        // add one step after start date seem to help exercise premium convergence for deep ITM case
        // and can be used to estimate theta, delta next day as well
        DateTime firstDate = AssetUtil::getHoliday(Underlier.get())->addBusinessDays(segDates[0], 1);
        critDates.push_back(DateTime(firstDate.getDate(),
                                            DateTime::END_OF_DAY_TIME));


        if (DivAmountTreatment) {
            this->isCall = isCall;
            this->noExerciseWindow = noExerciseWindow;
            Underlier = CAssetConstSP(
                PseudoSimpleEquity::create(Underlier.get(),
                                           divCritDates,
                                           isCall,
                                           noExerciseWindow));
        }

        InitVol();

        CTermStructure v_term;
        CalcV2Term(valDate, segDates[0], segDates[segDates.size()-1], v_term);
        // create time line
        timeLine = TimeLineSP(new CTimeLine());
        int effectiveStep = timeLine->CreateTimeLine(segDates, density, timeMetric, StepsPerYear, minGap, false, v_term,
                                                    critDates, equalTime);
        int totalStep = timeLine->NumOfStep;
        // number of nodes each side of central node
        Ndim = Maths::min((int)(TruncationStd*sqrt(effectiveStep*TreeAlpha)), effectiveStep);
//#ifdef NEW_TREE_TECH_TEST
        NLimit = Ndim;
        if (gammaNodesInterval>1)
            useGammaNodes=true;
        if (useGammaNodes){
            Ndim /= gammaNodesInterval;       
            // to avoid 0 Ndim.
            if (Ndim<1){
                int gNI_max = NLimit/2;
                throw ModelException(method, "gammaNodesInterval="
                                            + Format::toString(gammaNodesInterval)
                                            + " is too large.  Should be less than "
                                            + Format::toString(gNI_max)
                                            + ". Or increase StepsPerYear.");
            }
        }
//#endif
        NumOfInsertNode = numInsertNode;
        NumOfPrice = numOfPrice;

        // incl central node + 3 spare nodes on each side, make it [-n-3] to [n+3]
        // the pair of inner spare nodes assists branching along tree boundary
        // two outer spare nodes on each side close to floor (asset =0) and ceiling (some maximum asset level)
        // are for assisting extrapolations if needed

        if (useGammaNodes)
            range = TreeSliceGeneral::Range::create( 2, -NLimit-3, NLimit+3 );
        else
            range = TreeSliceGeneral::Range::create( 2, -Ndim-3, Ndim+3 );

        // other init
        if (DivAmountTreatment)
        {
            underlierPrice = DYNAMIC_POINTER_CAST< TreeSliceGeneral >( createSlice() );
            for( int i = 0; i < 2; ++i )
                UnderlierPrice[i] = (*underlierPrice)[i];
        }
        
        // allocate memory for original grid for GammaNodes
        if (useGammaNodes){
            gammaNodeMaps.clear();
            gammaNodeMaps.resize(totalStep+1);
        }

        // Allocate memory for cached shifts and nodes to preserve geometry across a tweak.
        if (UseSameGridGeometry && CacheMode == MAKE_CACHE)
        {
            int j_limit = 2*NLimit+7;
            if (NumOfInsertNode > j_limit){
                // this validation is to avoid accessing beyond of allocated array size.
                // It's rare to have such a big number of insert node, so leave it as this.
                // There is no critical reason to have such a limit.
                throw ModelException(method, "When UseSameGridGeometry is active, NumOfInsertNode["
                                            + Format::toString(NumOfInsertNode)
                                            + "] should be less than 2*NLimit+7 =["
                                            + Format::toString(j_limit)
                                            + "].");
            }

            CacheIdxUsed = new int*[totalStep + 1];
            CacheIdxLower = new int[totalStep + 1];
            CacheIdxUpper = new int[totalStep + 1];

            for (i=0; i<3; i++)
            {
                CacheBranching[i] = new int*[totalStep + 1];
                for (j=0; j<=totalStep; j++)
                {
                    CacheBranching[i][j] = new int[2*NLimit+7] +NLimit+3;
                }
            }

            for (i=0; i<2; i++)
            {
                CacheInsNodeUsed[i] = new bool*[totalStep + 1];
                CacheInsNodeIdx[i] = new int*[totalStep + 1];
                CacheInsNodeValue[i] = new double*[totalStep + 1];

                for (j=0; j<=totalStep; j++)
                {
                    CacheInsNodeUsed[i][j] = new bool[2*NLimit+7] +NLimit+3;
                    CacheInsNodeIdx[i][j] = new int[2*NLimit+7] +NLimit+3;
                    CacheInsNodeValue[i][j] = new double[2*NLimit+7] +NLimit+3;
                    // initialization of flag.  Formally, it's located in RollNode, 
                    // but it would be initialized many times when moveInsNode is true.
                    for (int k=-NLimit-3;k<=NLimit+3;k++)
                        CacheInsNodeUsed[i][j][k] = false;
                }
            }

            for (i=0; i<=totalStep; i++)
            {
                CacheIdxUsed[i] = new int[2*NLimit+7] +NLimit+3;
                CacheIdxLower[i] = 9999; // initial value
                CacheIdxUpper[i] = -9999; // initial value
            }
        }

        StepForward.resize(totalStep+1);
        CentreNode.resize(totalStep+1);

        // alloc storage for fwd induction prob's if required
        if (NumOfProbMarking>0)
        {
                int k;
                FwdProb = new double**[NumOfProbMarking];
                for (i=0; i<NumOfProbMarking; i++)
                {
                    FwdProb[i] = new double*[NumOfPrice];
                    for (j=0; j<NumOfPrice; j++)
                    {
                        FwdProb[i][j] = new double[2*NLimit+7] +NLimit+3;
                        for (k=-NLimit-3; k<=NLimit+3; k++)
                            FwdProb[i][j][k] = 0.0;
                    }
                }
        }

    }
    catch(exception& e){
        throw ModelException(&e, method, "tree initialisation failed");
    }
}

/**   tree nodes set up */
void CTree1f::NodeSetup(int step, int idx)
{
    BotDim[idx] = TopDim[idx] = Maths::min(step+1, Ndim);
    BotClip[idx] = Maths::min(BotClip[idx], BotDim[idx]);
    TopClip[idx] = Maths::min(TopClip[idx], TopDim[idx]);

    double drift = 1.0; // used only for rolling boundary nodes at step < last step
    if (step < timeLine->NumOfStep)
        drift = StepForward[step+1]/StepForward[step]; // step drift

    if (RollDirection==-1)
        drift = 1.0/drift;

    NodeSetupSpace(Stock[idx], -BotClip[idx]-1, TopClip[idx]+1,
                   Stock[1-idx], -BotClip[1-idx]-1, TopClip[1-idx]+1,
                   CentreNode[step], NodeSpace[GetTreeSeg(step)], step, drift);
}

//#ifdef NEW_TREE_TECH_TEST
// setting up Spot Level on Time Slice.
// It looks addedGridMap, whose node will be between the Normal Nodes.
// it's expected that bot = -BotClip-1 = range.bot-2) and top = TopClip+1 = range.top-2
void CTree1f::NodeSetWithAddNode(double* s, int bot, int top, double f0, double spacing){
    static const string method = "CTree1f::NodeSetWithAddNode";
    try
    {
        for (int j=bot;j<=top;j++){
            s[j] = f0*exp(spacing*(double)(gammaNodeMaps[CurrStep].getOrgIdxByBef(j)));
        }
    }
    catch(exception& e){
        throw ModelException(&e, method);
    }
}
//#endif

/**------------------------------------------------------------------------------
*   Description  :    calculate probabilities and node shifts required to match v_dt.
*                   An inserted node will be used if supplied and in such cases
*                   diffusion branches to 4 nodes.
*------------------------------------------------------------------------------*/
void CTree1f::CalcProb(double s1, double v_dt, int j, double s_ins,
                       int &u_shift, int &m_shift, int &d_shift,
                       double &pu, double &pd, double &pm, bool noRebranch)
{
    const int numOfNodeShift = 50; // allow 50 node shift operations
    double u, m, d;
    double a, pa=0.0, p1, p2, p3, temp[7];
    int i, iLast;

    double weight = 0.9;    // weight of inserted node.
    iLast = 1 - CurrIdx;

    // several weighting method can be choose
    if (DEBUG_InsertWidth<=0.0 || DEBUG_SmoothWeigth == 0)
        weight = 0.9;
    else if (s_ins > s1)
        weight = SmoothValue(0.9, s_ins, 0.0, s_ins*(1.0-DEBUG_InsertWidth), s1);
    else
        weight = SmoothValue(0.0, s_ins*(1.0+DEBUG_InsertWidth), 0.9, s_ins, s1);

    for (i=0; i<numOfNodeShift-1; i++)
    {
        m = Stock[iLast][j+m_shift]/s1 - 1.0;
        u = Stock[iLast][j+u_shift]/s1 - 1.0;
        d = Stock[iLast][j+d_shift]/s1 - 1.0;
        if ( Maths::isZero((u-m) * (u-d)) )
        {
            break;
        }
        pu = (v_dt + d*m)/(u-d)/(u-m);
        pd = (v_dt + u*m)/(d-u)/(d-m);
        if ((pu<1 && pu >=0 && pd<1 && pd>=0 && pu+pd<1) || (j+u_shift == TopClip[CurrIdx]+1 && pu > pd) ||
            (j+d_shift == -BotClip[CurrIdx]-1 && pd > pu) || i == numOfNodeShift)
        {// use inserted node if possible
            if (s_ins > 0 && (pu<1 && pu >=0 && pd<1 && pd>=0 && pu+pd<1))
            {// measure wrt m
                a = s_ins/s1 - 1.0 - m;
                u -= m;
                d -= m;
                temp[1] = u*(u-d);
                temp[2] = d*(d-u);
                temp[3] = a*(a-d);
                temp[4] = a*(a-u);
                temp[5] = v_dt+m*m+m*d;
                temp[6] = v_dt+m*m+m*u;
                p3 = (1-temp[5]/temp[1]-temp[6]/temp[2])/(1-temp[3]/temp[1]-temp[4]/temp[2]);
                if (temp[3]>0)
                    p1 = temp[5]/temp[3];
                else
                    p1 = (temp[5]-temp[1])/temp[3];
                if (temp[4]>0)
                    p2 = temp[6]/temp[4];
                else
                    p2 = (temp[6]-temp[2])/temp[4];
                p1 = (p1<0? 1:p1);
                p2 = (p2<0? 1:p2);
                p3 = (p3<0? 1:p3);
                pa = weight*Maths::min(Maths::min(p1, p2) ,p3); // 90% of maximum used to have less impact on 3rd moment, sensitive to barrier on dates only
                pu = (temp[5] - temp[3]*pa)/temp[1];
                pd = (temp[6] - temp[4]*pa)/temp[2];
                //ASSERT(pu>=0 && pa>=0 && pd>=0 && pu+pa+pd-FP_MIN <1.0);
            }
            break;
        }
        else if (noRebranch)
        {
            break; // rebranching details have already been computed
        }
        else if (Stock[iLast][j+d_shift] <= FP_MIN)
        {
           break; // stock = 0
        }
        else
        {
            bool reachedTop = (j+u_shift >= TopClip[CurrIdx]+1);
            bool reachedBot = (j+d_shift <= -BotClip[CurrIdx]-1);

            if ( pd < 0 || u < 0 )
            {// shift every node
                if (reachedTop)
                    break; // reached tree top edge
                else{
                    u_shift++; m_shift++; d_shift++;
                }
            }
            else if (pu < 0 || d >0)
            {
                if (reachedBot)
                    break; // reached tree bottom edge
                else{
                    u_shift--; m_shift--; d_shift--;
                }
            }
            // just use simple opening up. Single side open up is removed as it needs more refined decision test to avoid 3rd moment error.
            else if (reachedTop || reachedBot)
                break; // reached one edge and cannot open up
            else{
                u_shift ++; d_shift --;
            }
        }
    }
    pm = 1.0 - (pu + pd + pa);
}

/**------------------------------------------------------------------------------
*   Description  :      interpolate price for a given asset level in a tree.
*       Parameters       :      s                               stock level for interpolation
*                                       useCurrentArray false=use array at next step (t+1), true=use array at this step.
*                                       iPrice                  index array for Price
*                                       j                               estimated starting point in tree
*
*   Returns      :  price at given asset level if RollDirection==-1 or weight at j_out for prob calc
*------------------------------------------------------------------------------*/

// to add spline interpolation for coping tree step density change and possible tree stacks
double CTree1f::TreeInterp(double s, bool useCurrentArray,
                           int iPrice, int j, int *j_out)
{
    static const string method = "CTree1f::TreeInterp";
    try {
    int shift = 0;
    double result;
    int iArr, bot_dim, top_dim;

    if (useCurrentArray)
            iArr = CurrIdx;
    else
            iArr = 1-CurrIdx;

        bot_dim = BotClip[iArr];
    top_dim = TopClip[iArr];

    if (RollDirection == 1) // fwd rolling for prob
    {
            if(j_out == 0 || useCurrentArray)
                throw ModelException(method, "fwd rolling failed");
    }

    if (s <= FP_MIN)
    {
            if (RollDirection == -1) // pricing
                result = NodePrice[iArr][iPrice][-bot_dim-3]; // floor price for asset = 0
            else // fwd rolling for prob
            {
                *j_out = -bot_dim-3;
                result = 1.0;
            }
    }
    else
    {// decide branching
        while (j >= -bot_dim && j <= top_dim)
        {
            if (Stock[iArr][j] > s)
                j --;
            else if (Stock[iArr][j+1] < s)
                j ++;
            else
                break;
        }
        if (j < -bot_dim-1)
            shift = -j-bot_dim; // shift up node for extrapolation
        else if (j >= top_dim+1)
            shift = -j+top_dim-1; // shift down for extrapolation

        // linearly interpolate
        if (Stock[iArr][j+shift+1] <= Stock[iArr][j+shift]) {
            throw ModelException(method,
                                    "Asset level at [" + Format::toString(iArr) +
                                    "," + Format::toString(j+shift+1) + "] (" +
                                    Format::toString("%.12f", Stock[iArr][j+shift+1]) +
                                    ") is <= asset level at [" +
                                    Format::toString(iArr) + "," +
                                    Format::toString(j+shift) + "] (" +
                                    Format::toString("%.12f",Stock[iArr][j+shift]) +
                                    ") for asset " + Underlier->getName() +
                                    "  at time step = " + Format::toString(CurrStep));
        }

        if (RollDirection == -1) // pricing
            result = NodePrice[iArr][iPrice][j+shift]+(s-Stock[iArr][j+shift])*(NodePrice[iArr][iPrice][j+1+shift]-NodePrice[iArr][iPrice][j+shift])
                /(Stock[iArr][j+shift+1]-Stock[iArr][j+shift]);
        else // fwd rolling for prob
        {
            *j_out = j+shift;
            while (*j_out < top_dim+1 && s > Stock[iArr][*j_out+1])
                *j_out += 1;
            while (*j_out > -bot_dim-2 && s < Stock[iArr][*j_out])
                *j_out -= 1;
            result = (Stock[iArr][*j_out+1]-s)/(Stock[iArr][*j_out+1]-Stock[iArr][*j_out]); // approx
        }
        // refine tree boudary node extrapolation
        if (RollDirection == -1 // pricing
            && (s < Stock[iArr][-bot_dim] || s > Stock[iArr][top_dim])) // outside tree edge
        {
            double refined;
            int start = (s < Stock[iArr][-bot_dim]? -bot_dim-3 : top_dim);
            // try to improve with floor or ceiling price
            if (LinearInterpImprove(&Stock[iArr][start], &NodePrice[iArr][iPrice][start], 4, s, refined, false) > 0)
                result = refined;
        }
    }
    ASSERT((result >= 0 && result <=1) || RollDirection == -1);

    return result;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**------------------------------------------------------------------------------
*   Description  :      linear interpolate similar to TreeInterp, 
*                       but using explicit double array rather than Stock[] & NodePrice[].
*   Parameters   :      s           stock level for interpolation
*                       s_inp       input array of stock.
*                       p_inp       input array of price.
*                       j           estimated starting point in s_inp
*
*   Returns      :  price at given asset level.  
*------------------------------------------------------------------------------*/

// to add spline interpolation for coping tree step density change and possible tree stacks
double LinearTreeInterp(const double s, 
                        const vector<double> & s_inp, 
                        const vector<double> & p_inp,
                        int j)
{
    static const string method = "CTree1f::LinearTreeInterp";
    try {

    int botBuff = 3; //number of nodes between lowest to BotClip.
    int topBuff = 3; //number of nodes between top to TopClip.

    int arraySize = p_inp.size();
    //validation
    if (s_inp.size() != arraySize){
        ModelException(method, "s_inp [" + Format::toString(arraySize) 
                            + "] and p_inp["
                            + Format::toString(arraySize) + "] is not same!.");
    }
    if (botBuff < 0 || botBuff >= arraySize ||topBuff > arraySize){
        ModelException(method, "botBuff[" + Format::toString(botBuff) 
                             + "] or topBuff["
                             + Format::toString(topBuff) 
                             + "] would be out of array size["
                             + Format::toString(arraySize) + "].");
    }

    int shift = 0;
    double result;

    int bot_dim = botBuff;
    int top_dim = s_inp.size()-topBuff;

    if (s <= FP_MIN)
    {
        result = p_inp[0]; // floor price for asset = 0
    }
    else
    {// decide branching
        while (j >= bot_dim && j <= top_dim)
        {
            if (s_inp[j] > s)
                j --;
            else if (s_inp[j+1] < s)
                j ++;
            else
                break;
        }
        if (j < bot_dim-1)
            shift = -j+bot_dim; // shift up node for extrapolation
        else if (j >= top_dim+1)
            shift = -j+top_dim-1; // shift down for extrapolation

        // check j
        if (j+shift < 0 || arraySize <= j+shift+1){
            ModelException(method, "j(=" +Format::toString(j)
                                 + ") + shift(=" + Format::toString(shift) 
                                 + ") is out of array size = "
                                 + Format::toString(arraySize)+".");
        }

        // linearly interpolate
        if (s_inp[j+shift+1] <= s_inp[j+shift]) {
            throw ModelException(method,
                                    "s_inp is not increasing order at [" + Format::toString(j+shift+1) + "].");
        }

        result = p_inp[j+shift]+(s-s_inp[j+shift])*(p_inp[j+1+shift]-p_inp[j+shift])
            /(s_inp[j+shift+1]-s_inp[j+shift]);
        // refine tree boudary node extrapolation
        if ((s < s_inp[bot_dim] || s > s_inp[top_dim])) // outside tree edge
        {
            double refined;
            int start = (s < s_inp[bot_dim]? 0 : top_dim);
            // try to improve with floor or ceiling price
            if (start>=3 || s_inp.size()-top_dim >=3){
                if (LinearInterpImprove(&s_inp[start], &p_inp[start], 4, s, refined, false) > 0)
                    result = refined;
            }
        }
    }

    return result;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**------------------------------------------------------------------------------
*   Description  :    called by Roll() after looping through all nodes.
*                    interpolates 2 extra node prices for tree roll over use.
*------------------------------------------------------------------------------*/
void CTree1f::InterpTreeBoundary()
{
    static const string method = "CTree1f::InterpTreeBoundary";
    try {
        int i, j;
        // interpolate 2 extra node prices for roll
        // cannot use TreeInterp because it needs this point itself
        double result[2], x[4], y[4], *s, *p, node[2];
        node[0] = Stock[CurrIdx][-BotClip[CurrIdx]-1];
        node[1] = Stock[CurrIdx][TopClip[CurrIdx]+1];
        for (i=0; i<NumOfPrice; i++)
        {
            for (j=0; j<2; j++)
            {
                int start = (j==0? -BotClip[CurrIdx]-3 : TopClip[CurrIdx]-1);
                s = &Stock[CurrIdx][start];
                p = &NodePrice[CurrIdx][i][start];
                x[0] = s[0], y[0] = p[0];
                x[1] = s[1], y[1] = p[1];
                x[2] = s[3], y[2] = p[3];
                x[3] = s[4], y[3] = p[4];
                if (LinearInterpImprove(x, y, 4, node[j], result[j], false) < 0) {
                    throw ModelException(method, "interpolation failed.");
                }
            }
            NodePrice[CurrIdx][i][-BotClip[CurrIdx]-1] = result[0];
            NodePrice[CurrIdx][i][TopClip[CurrIdx]+1] = result[1];
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//------------------------------------------------------------------------------
/*** forward or backward roll the tree at one time step
*------------------------------------------------------------------------------*/
void CTree1f::RollTree(vector<double> *drift_out, vector<double> *driftIns_out)
{
    static const string method = "CTree1f::Roll";
    try {
     int j, iPrice;

        ASSERT(CurrStep < timeLine->NumOfStep);

        double df = (RollDirection == -1?
                     discYC->pv(timeLine->StepDates[CurrStep], timeLine->StepDates[CurrStep+1]) : 1.0);

        // start with usual array
        const vector< double * > & arr = NodePrice[CurrIdx];

        // indicies
        int pStart = BotArrIdx[CurrIdx];
        int pEnd = TopArrIdx[CurrIdx];
        int iLast = 1-CurrIdx;

        // fwd and variance array for this time slice
        vector<double> var_dt;
        // we may need drift outputs
        vector<double> mem_drift; // not to be returned
        vector<double>* ptr_drift = (drift_out? drift_out:&mem_drift);
        vector<double>* ptr_driftIns = (driftIns_out? driftIns_out:&mem_drift);

        // calculate the fwd and variance at this step
        bool fastRoll = CalcStepDriftAndVar(Stock[CurrIdx], -BotClip[CurrIdx], TopClip[CurrIdx], var_dt, ptr_drift);

#ifdef TREE_DEBUG_VOL_FILE
        if ((CurrStep == timeLine->NumOfStep-1) || (CurrStep < 10) )
            Debug_OutPutVol(CurrStep, timeLine->TradeTime[CurrStep+1], var_dt, timeLine->TradeYrFrac[CurrStep+1], Stock[CurrIdx], -BotClip[CurrIdx], TopClip[CurrIdx], CurrStep == timeLine->NumOfStep-1);
#endif


        int numDrift = ptr_drift->size();
        int numVar = var_dt.size();
        // must be 1 per step or one per node level for now
        if (numVar != 1 && numVar != BotClip[CurrIdx]+TopClip[CurrIdx]+1) {
            throw ModelException(method,
                                 "incorrect num of tree step variances.");
        }

        // loop 2 floor and 2 ceiling nodes first, vol and fwd not needed
        RollNode(arr, 0, pStart, pEnd, 0, df, -BotClip[CurrIdx]-3, false);
        RollNode(arr, 0, pStart, pEnd, 0, df, -BotClip[CurrIdx]-2, false);
        RollNode(arr, 0, pStart, pEnd, 0, df, TopClip[CurrIdx]+2, false);
        RollNode(arr, 0, pStart, pEnd, 0, df, TopClip[CurrIdx]+3, false);

        double v_dt = var_dt[0];
        double mu = (*ptr_drift)[0];

        if (fastRoll)
        {// fast rolling
            double* arrLast;
            for (int iPrice=pStart; iPrice<=pEnd; iPrice++)
            {
                arrLast = NodePrice[iLast][iPrice];
                for (j= -BotClip[CurrIdx]; j<=TopClip[CurrIdx]; j++)
                    arr[iPrice][j] = df*(Pd*arrLast[j-1] + Pm*arrLast[j] + Pu*arrLast[j+1]);
            }
        }
        else
        {// roll one node at a time
            // decide if each price array wants to move the inserted node
            bool moveInsNode = NumOfInsertNode > 0;
            moveInsNode = moveInsNode && latticeProd->moveInsertNode(CurrStep, pStart);

            for (j= -BotClip[CurrIdx]; j<=TopClip[CurrIdx]; j++)
            {
                if (numDrift > 1)
                    mu = (*ptr_drift)[j+BotClip[CurrIdx]];
                if (numVar >1)
                    v_dt = var_dt[j+BotClip[CurrIdx]]; // local vol

                // decide if needs to move inserted node for each price array
                if (moveInsNode)
                {
                    for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                    {
                        latticeProd->moveInsertNode(CurrStep, iPrice);
                        RollNode(arr, mu*Stock[CurrIdx][j], iPrice, iPrice, v_dt, df, j, false);
                    }
                }
                else
                    RollNode(arr, mu*Stock[CurrIdx][j], pStart, pEnd, v_dt, df, j, false);
            }
        }

        if (NumOfInsertNode > 0)
        {
            for (iPrice=pStart; iPrice<=pEnd; iPrice++)
            {
                latticeProd->moveInsertNode(CurrStep, iPrice); // a chance for product to change inserted node

                // reprocess nodes around inserted nodes if needed
                if (fastRoll)
                {// reprocess nodes around inserted nodes if needed
                    for (j= -BotClip[iLast]; j<=TopClip[iLast]; j++)
                    {
                        for (int k=0; k<NumOfInsertNode; k++)
                        {// use 3*node spacing to check for now
                            if (InsStock[iLast][k] > 0 &&
                                fabs(InsStock[iLast][k]-Stock[iLast][j])<3*(Stock[iLast][j+1]-Stock[iLast][j]))
                            {
                                RollNode(arr, mu*Stock[CurrIdx][j], iPrice, iPrice, v_dt, df, j, false);
                                break; // done for this j
                            }
                        }
                    }
                }
                // process inserted nodes itself
                if (numVar >1 || numDrift >1) // calculate drift and variance for inserted node if needed
                {
                    CalcStepDriftAndVar(InsStock[CurrIdx], 0, NumOfInsertNode-1, var_dt, ptr_driftIns);
                    v_dt = var_dt[0];
                    mu = (*ptr_driftIns)[0];
                }
                // loop inserted nodes
                for (j=0; j<NumOfInsertNode; j++)
                {   // do inserted node itself
                    if (ptr_driftIns->size() >1)
                        mu = (*ptr_driftIns)[j];
                    if (var_dt.size()> 1)
                        v_dt = var_dt[j]; // local vol

                    RollNode(InsNodePrice[CurrIdx], mu*InsStock[CurrIdx][j], iPrice, iPrice, v_dt, df, j, true);
                }
            }
        }

        // interpolate one more point each side
        InterpTreeBoundary();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//------------------------------------------------------------------------------
/*** rolling each node
*------------------------------------------------------------------------------*/
void CTree1f::RollNode(const vector< double * > & arr, double s1, int pStart, int pEnd,
                       double v_dt, double df, int j_org, bool isInsertNode)
{
    static const string method = "CTree1f::RollNode";
    try {
    // j , j_curr : Node Idx of usual Tree (-Bot-3 to Top+3).
    // j_org      : Node Idx of Tree or InsNode.  When insNode, could be above Top+3.
    int j=j_org, j_curr=j_org, iLast, d_shift =-1, m_shift=0, u_shift =1;
    int j_ins=0;
    double pu=Pu, pd=Pd, pm=Pm;
    double s_ins = 0.0;
    int ins_use = -1;    // index of inserted node found
    int iPrice;
    bool priceInterplated = false;
    bool probInterplated = false;
    bool nodeChangedForTweak = false;

    double s0;
    if(isInsertNode)
        s0 = InsStock[CurrIdx][j_org];
    else
        s0 = Stock[CurrIdx][j_org];

    iLast = 1- CurrIdx;

    // When isInsertNode = true, j_org is index of NumOfInsert Node.
    // But j_curr are assumed Tree Node's location (from -Bot to Top).
    // Two cases are not same at all.
    // Thus, when NumOfInsertNode is larger than TopClip[iLast],
    // it fails to find the "j" which correspoinding to InsStock level.
    // In such case, putting j and j_curr to 0 so as to start from centre node.
    // j_curr needs to be in Bot-Top,
    // othewise it's captured by "// floor and ceiling price simply carried over" part.
    if (isInsertNode && (j <= -BotClip[iLast] || j >= TopClip[iLast]) ){
                        // no need to check bottom, as j should be positive, but in case...
            j = j_curr = 0;
    }

    if ((s0 < FP_MIN || s0 > Stock[CurrIdx][TopClip[CurrIdx]+3]) && RollDirection == 1)
    {// prob sent to floor or ceiling
            for (iPrice=pStart; iPrice<=pEnd; iPrice++)
            {
                if (s0 < FP_MIN)
                    NodePrice[iLast][iPrice][-BotClip[iLast]-3] += arr[iPrice][j_org];
                else
                    NodePrice[iLast][iPrice][TopClip[iLast]+3] += arr[iPrice][j_org];
            }

            return;
    }

    // floor and ceiling price simply carried over
    if (j_curr <= -BotClip[CurrIdx]-2 || s0 < FP_MIN || j_curr >= TopClip[CurrIdx]+2)
    {
            if (j_curr == -BotClip[CurrIdx]-3 || s0 < FP_MIN)
                j = -BotClip[iLast]-3;
            else if (j_curr == -BotClip[CurrIdx]-2)
                j = -BotClip[iLast]-2;
            else if (j_curr == TopClip[CurrIdx]+3)
                j = TopClip[iLast]+3;
            else
                j = TopClip[iLast]+2;
            if (RollDirection == -1) // pricing
            {
                for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                    arr[iPrice][j_org] = df*NodePrice[1-CurrIdx][iPrice][j];
            }
            else // fwd prob calc
            {
                for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                {// df here is NOT disc factor but total diffusion prob
                    NodePrice[iLast][iPrice][j] += df*arr[iPrice][j_org];
                }
            }
            return;
    }

    if (PerNodeProb || (NumOfInsertNode>0 && SmoothMethod==NODE_INSERTION)
            || (pu<0.0 || pm<0.0 || pd <0.0)) // require prob calc per node
    {
            // decide branching
            while (j > -BotClip[iLast] && j < TopClip[iLast])
            {
                if ((Stock[iLast][j]+Stock[iLast][j-1])/2.0 > s1)
                    j --;
                else if ((Stock[iLast][j+1]+Stock[iLast][j])/2.0 < s1)
                    j ++;
                else // now between j and j+1, unless beyond tree boundary
                    break;
            }


            if (UseCacheAtNode(j_curr) && !isInsertNode)
            {
                int j_cache = CacheIdxUsed[CurrStep][j_curr];
                int topIndex = j < TopClip[iLast] ? j_cache + 1 : TopClip[iLast];
                int botIndex = j > -BotClip[iLast] ? j_cache - 1 : -BotClip[iLast];
                if (Stock[iLast][botIndex] <= s1 && s1 <= Stock[iLast][topIndex])
                {
                    // encompassed by outer nodes, can use cached j
                    j = j_cache;
                }
                else
                {
                    nodeChangedForTweak = true;
                }
            }

            if (Stock[iLast][j-1] <= FP_MIN || s1 < Stock[iLast][-BotClip[iLast]] || s1 > Stock[iLast][TopClip[iLast]])
            {// close to 0 or outside tree edge, simple interpolation used
                for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                {
                    if (RollDirection == -1)
                    {
                        arr[iPrice][j_org] = df*TreeInterp(s1, false, iPrice, j);
                    }
                    else
                    {// fwd prob calc
                        pd = TreeInterp(s1, false, iPrice, j, &j);
                        pm = 1.0 - pd;
                        // df here is NOT disc factor but total diffusion prob
                        NodePrice[iLast][iPrice][j] += df*pd*arr[iPrice][j_org];
                        NodePrice[iLast][iPrice][j+1] += df*pm*arr[iPrice][j_org];
                    }
                }
                priceInterplated = probInterplated = true;
            }
            else
            {
                // insert node method
                if (SmoothMethod == NODE_INSERTION) // use inserted node method
                {
                    // always try to compute, even if we are using cache.
                    for (j_ins=0; j_ins < NumOfInsertNode; j_ins ++)
                        {// find if an inserted level between relavent range exist
                            bool tmpbool;
                            if (DEBUG_InsertWidth <=0.0)
                                tmpbool = InsPriority[iLast][j_ins] >=0 && // if this inserted node is "on" at this step
                                    (InsStock[iLast][j_ins] > (3*Stock[iLast][j-1]-Stock[iLast][j])/2.0
                                     && InsStock[iLast][j_ins] < (3*Stock[iLast][j+1]-Stock[iLast][j])/2.0);
                            else
                                tmpbool =
                                    InsPriority[iLast][j_ins] >=0 && // if this inserted node is "on" at this step
                                    (Stock[iLast][j] > InsStock[iLast][j_ins]*(1.0-DEBUG_InsertWidth)
                                     && Stock[iLast][j] < InsStock[iLast][j_ins]*(1.0+DEBUG_InsertWidth) );
                            if (tmpbool)
                            {
    //                    if (InsPriority[iLast][j_ins] >=0 && // if this inserted node is "on" at this step
    //                        (InsStock[iLast])[j_ins] > (3*Stock[iLast][j-1]-Stock[iLast][j])/2.0
    //                         && InsStock[iLast])[j_ins] < (3*Stock[iLast][j+1]-Stock[iLast][j])/2.0))
    //                    {
                                if (ins_use < 0)
                                    ins_use = j_ins; // found first
                                else
                                {// found more than one
                                    if (InsPriority[iLast][j_ins] > InsPriority[iLast][ins_use])
                                        ins_use = j_ins; // replace with higher priority one
                                    else if (InsPriority[iLast][j_ins] == InsPriority[iLast][ins_use])
                                    {// same priority, choose closer node
                                        if (fabs(s1-InsStock[iLast][j_ins]) < fabs(s1-InsStock[iLast][ins_use]))
                                            ins_use = j_ins;
                                    }
                                }
                            }
                        }
                        if (ins_use>=0)
                        {// found inserted node to use
                            j_ins = ins_use;
                            s_ins = InsStock[iLast][j_ins];
                            // if pricing, save decision and node
                            if (CacheMode == MAKE_CACHE)
                            {
                                CacheInsNodeUsed[isInsertNode? 1:0][CurrStep][j_curr] = true;
                                CacheInsNodeIdx[isInsertNode? 1:0][CurrStep][j_curr] = j_ins;
                                CacheInsNodeValue[isInsertNode? 1:0][CurrStep][j_curr] = s_ins;
                            }
                        }
                }

                if (CacheMode == USE_CACHE && !nodeChangedForTweak) // need to index the insert node
                {  // For using cache, if we used an inserted node in the base, use the new inserted node level, if computed.  Otherwise,
                   // use original node
                    // only use the cache ins nodes, when isActiveInsNode = true.
                    if (CacheInsNodeUsed[isInsertNode? 1:0][CurrStep][j_curr] && isActiveInsNode)
                    {
                        if (Maths::isZero(s_ins))
                        {
                            j_ins = CacheInsNodeIdx[isInsertNode? 1:0][CurrStep][j_curr];
                            // check it's actually active or not, because moveInsNode could make it off.
                            if (InsPriority[iLast][j_ins] >=0)   
                                s_ins = CacheInsNodeValue[isInsertNode? 1:0][CurrStep][j_curr];
                            else
                                s_ins = 0.0;
                        }
                    }
                    else
                    {
                        // if we didn't use an inserted node for the base pricing, do not use for a tweak.
                        s_ins = 0.0;
                    }

                }

                // need to also check that j is still in the geometry
                bool noRebranch = UseCacheAtNode(j_curr) && !nodeChangedForTweak && !isInsertNode;

                // if doing a same grid tweak, precompute shifts
                if (noRebranch)
                {
                    u_shift = CacheBranching[U_SHIFT][CurrStep][j_curr];
                    m_shift    = CacheBranching[M_SHIFT][CurrStep][j_curr];
                    d_shift    = CacheBranching[D_SHIFT][CurrStep][j_curr];
                }


                if (s_ins != 0.0 || PerNodeProb || j != j_curr) // require prob calc
                    CalcProb(s1, v_dt, j, s_ins, u_shift, m_shift, d_shift, pu, pd, pm, noRebranch);
                bool badProb = (pd + pu > 1 || pu >1 || pd >1 || pd < 0 || pu < 0);

                // use record shifts to preserve geometry for same grid tweaks
                if (CacheMode == MAKE_CACHE && !isInsertNode)
                {
                    CacheBranching[U_SHIFT][CurrStep][j_curr] = u_shift;
                    CacheBranching[M_SHIFT][CurrStep][j_curr] = m_shift;
                    CacheBranching[D_SHIFT][CurrStep][j_curr] = d_shift;
                    CacheIdxUsed[CurrStep][j_curr] = j;
                    if (j_curr < CacheIdxLower[CurrStep]) CacheIdxLower[CurrStep] = j;
                    if (j_curr > CacheIdxUpper[CurrStep]) CacheIdxUpper[CurrStep] = j;
                }

                if (badProb)
                { // unable to obtain correct prob's
                    for (int iter = 0; iter < 2; iter++) // call this twice to guarantee pm = 0.000000
                    //(if pu, pd too large, will not get exactly 0 due to numerical precision issues)
                    {
                        if (pu >0.0 && pd > 0.0) // v_dt too large, too close to lower boundary
                        {    // maintain the mean, pm is -ve, so make it 0
                            double common = pm/(Stock[iLast][j+u_shift] - Stock[iLast][j+d_shift]);
                            pu += -common*(Stock[iLast][j+d_shift] - Stock[iLast][j+m_shift]) - FP_MIN;
                            pd += common*(Stock[iLast][j+u_shift] - Stock[iLast][j+m_shift]) - FP_MIN;
                            pm = 1.0 - (pu + pd);
                        }
                    }
                    if (pu <0.0 || pd <0.0){ // use interpolation at s1
                        pu = v_dt*s1*s1/(Stock[iLast][j+u_shift] - Stock[iLast][j+d_shift])/(Stock[iLast][j+u_shift]-s1);
                        pd = v_dt*s1*s1/(Stock[iLast][j+u_shift] - Stock[iLast][j+d_shift])/(s1-Stock[iLast][j+d_shift]);
                        if (pu + pd > 1.0){
                            pm = pu + pd + FP_MIN; // just to store sum to renormalize prob's
                            pu /= pm;
                            pd /= pm;
                        }
                        pm = 1.0 - (pu + pd);
                    }
                    badProb = (!Maths::finite(pd) || !Maths::finite(pu) 
                                   || pd + pu > 1 || pu >1 || pd >1 || pd < 0 || pu < 0);

                    // check the variance for same Geometry case.
                    // If it's not close to the v_dt, interploation could give some bias.
                    // Check the usual rebrunching availability (u_shift,m_shift,d_shift) = 1, 0, -1,  
                    // and use avoid interpolation, if it's possible.
                    // when it has insert node, use interpolation.
//                    bool doInterp = true;
//                    if (!badProb && noRebranch
//                        && (   !(j+u_shift ==  TopClip[CurrIdx]+1 && pu > pd)       // this two condition means 
//                            && !(j+d_shift == -BotClip[CurrIdx]-1 && pd > pu) )     // CalcProb return the same results.
//                        && Maths::finite(v_dt) )
//                        //&& Maths::isZero(s_ins) )
//                    {
//                        double u = Stock[iLast][j+u_shift]/s1-1.0;
//                        double d = Stock[iLast][j+d_shift]/s1-1.0;
//                        double var = pu*u*u+pd*d*d; // this is local variance if use interpolation.
//                        if (fabs(var/v_dt-1.0)> usggVarDiffAllowance)
//                        {// use relative diffs as criteria.  Some cases, it's necessary to give up 
//                            // otherwise it would generate bias, often be seen in vega_pointwise.
//                            double pu_save = pu, pm_save = pm, pd_save = pd;
//                            int u_shift_save = u_shift, m_shift_save = m_shift, d_shift_save = d_shift;
//                            u_shift = 1, m_shift = 0, d_shift = -1;
//                            CalcProb(s1, v_dt, j, s_ins, u_shift, m_shift, d_shift, pu, pd, pm, false);
//                            badProb = (pd + pu > 1 || pu >1 || pd >1 || pd < 0 || pu < 0);
//                            if (badProb){               // restore geometry and use Interp.
//                                u_shift = u_shift_save, m_shift = m_shift_save, d_shift = d_shift_save;
//                                pu = pu_save, pm = pm_save, pd = pd_save;
//                                badProb = false;        //before coming here, it must be false.
//                            }
//                            else
//                                doInterp = false;       // no more interpolate!
//                        }
//                    }

//                    if (doInterp){
                    if (!badProb){
                        for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                        {
                            if (RollDirection == -1)
                            {
                                arr[iPrice][j_org] = pm*TreeInterp(s1, false, iPrice, j); // price for node at s1
                                arr[iPrice][j_org] += pd*NodePrice[iLast][iPrice][j+d_shift] + pu*NodePrice[iLast][iPrice][j+u_shift];
                                arr[iPrice][j_org] *=df;
#ifdef TREE_DEBUG_FILE
                                if (iPrice==0 && DEBUG_InsertWidth < -2.0){
                                    string fileName = TREE_DEBUG_FILE;
                                    double u = Stock[iLast][j+u_shift]/s1-1.0;
                                    double m = 0.0;
                                    double d = Stock[iLast][j+d_shift]/s1-1.0;
                                    double a = s_ins/s1-1.0;
                                    double pa = 1.0 - (pu+pm+pd);
                                    double mean = pu*u+pm*m+pd*d+pa*a;
                                    double var = pu*u*u+pm*m*m+pd*d*d+pa*a*a;
                                    double skew = pu*u*u*u+pm*m*m*m+pd*d*d*d+pa*a*a*a;
                                    Debug_OutPutTreeProbs(Stock[iLast][j_org], s1 , v_dt, s_ins, mean, var, skew*1.0e5, pu, pd, pm,CurrStep,j_org, d_shift,m_shift,u_shift,true,fileName,false);
                                    //Debug_OutPutTreeProbs(Stock[iLast][j_org] , pu, pd, pm,BotClip[iLast],j_org, d_shift,m_shift,u_shift,fileName,false);
                                }
#endif
                            }
                            else
                            {
                                // more consideration on how implement here, approx ok for now
                            }
                        }
                    }
                    else{
                        // if still prob are not good, give up the same geometry.
                        // interpolate the value from 50% up and down stock level.  (Binomial Tree).
                        pu = 0.5; pd = 0.5; pm = 0.0;
                        double s_up = s1 * (1.0 + sqrt(v_dt));
                        double s_dn = s1 * (1.0 - sqrt(v_dt));
                        for (iPrice=pStart; iPrice<=pEnd; iPrice++)
                        {
                            if (RollDirection == -1)
                            {
                                arr[iPrice][j_org]  = TreeInterp(s_up, false, iPrice, j); // price of Up Node
                                arr[iPrice][j_org] += TreeInterp(s_dn, false, iPrice, j); // price of Down Node
                                arr[iPrice][j_org] *= df*0.5;
#ifdef TREE_DEBUG_FILE
                                if (iPrice==0 && DEBUG_InsertWidth < -2.0){
                                    string fileName = TREE_DEBUG_FILE;
                                    double u = s_up/s1-1.0;
                                    double m = 0.0;
                                    double d = s_dn/s1-1.0;
                                    double a = s_ins/s1-1.0;
                                    double pa = 1.0 - (pu+pm+pd);
                                    double mean = pu*u+pm*m+pd*d+pa*a;
                                    double var = pu*u*u+pm*m*m+pd*d*d+pa*a*a;
                                    double skew = pu*u*u*u+pm*m*m*m+pd*d*d*d+pa*a*a*a;
                                    Debug_OutPutTreeProbs(Stock[iLast][j_org], s1 , v_dt, s_ins, mean, var, skew*1.0e5, pu, pd, pm,CurrStep,j_org, d_shift,m_shift,u_shift,true,fileName,false);
                                    //Debug_OutPutTreeProbs(Stock[iLast][j_org] , pu, pd, pm,BotClip[iLast],j_org, d_shift,m_shift,u_shift,fileName,false);
                                }
#endif
                            }
                            else
                            {
                                // more consideration on how implement here, approx ok for now
                            }
                        }
                    }
                    priceInterplated = true;
//                    }
                }
            }
    }

    ASSERT((pd>-FP_MIN && pu>-FP_MIN && pm >-FP_MIN) || probInterplated);
    if (RollDirection == -1 && !priceInterplated)
    {
            for (iPrice=pStart; iPrice<=pEnd; iPrice++)
            {
                arr[iPrice][j_org] = df*(pd*NodePrice[iLast][iPrice][j+d_shift] +
                                          pm*NodePrice[iLast][iPrice][j+m_shift] +
                                          pu*NodePrice[iLast][iPrice][j+u_shift]);
                if (s_ins > 0)
                    arr[iPrice][j_org] += df*(1-pu-pm-pd)*InsNodePrice[iLast][iPrice][j_ins];

#ifdef TREE_DEBUG_FILE
                if (iPrice==0 && DEBUG_InsertWidth < -2.0)
                {
                    string fileName = TREE_DEBUG_FILE;
                    double u = Stock[iLast][j+u_shift]/s1-1.0;
                    double m = Stock[iLast][j+m_shift]/s1-1.0;
                    double d = Stock[iLast][j+d_shift]/s1-1.0;
                    double a = s_ins/s1-1.0;
                    double pa = 1.0 - (pu+pm+pd);
                    double mean = pu*u+pm*m+pd*d+pa*a;
                    double var = pu*u*u+pm*m*m+pd*d*d+pa*a*a;
                    double skew = pu*u*u*u+pm*m*m*m+pd*d*d*d+pa*a*a*a;
                    Debug_OutPutTreeProbs(Stock[iLast][j_org], s1 , v_dt, s_ins, mean, var, skew*1.0e5, pu, pd, pm,CurrStep,j_org, d_shift,m_shift,u_shift,priceInterplated,fileName,false);
                    //Debug_OutPutTreeProbs(Stock[iLast][j_org] , pu, pd, pm,BotClip[iLast],j_org, d_shift,m_shift,u_shift,fileName,false);

                }
#endif
            }
    }
    else if (RollDirection == 1 && !probInterplated)
    {
            for (iPrice=pStart; iPrice<=pEnd; iPrice++)
            {
                // df here is NOT disc factor but total diffusion prob
                NodePrice[iLast][iPrice][j+d_shift] += df*pd*arr[iPrice][j_org];
                NodePrice[iLast][iPrice][j+m_shift] += df*pm*arr[iPrice][j_org];
                NodePrice[iLast][iPrice][j+u_shift] += df*pu*arr[iPrice][j_org];
                if (s_ins > 0)
                    InsNodePrice[iLast][iPrice][j_ins] += df*(1-pu-pm-pd)*arr[iPrice][j_org];
            }
    }
    if (priceInterplated && j_org>-BotClip[CurrIdx]+3 && j_org<TopClip[CurrIdx]-3)
        numOfInterp ++;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** changing time line segment of a different density */
bool CTree1f::SwitchSegment()
{
    bool changed = false;
    int i;
    int idx = CurrIdx;

    if (CurrSeg>0 && timeLine->SegmentEnd[CurrSeg-1]>=CurrStep && RollDirection==-1)
    {
        // store nodes before segment change
        vector<double> s_arr(2*NLimit+7);
        vector<double>::iterator sBase =
            s_arr.begin()+ NLimit+3;  // for convenience
        for (i=-BotClip[idx]-3; i<=TopClip[idx]+3; i++)
        {// copy all nodes
            sBase[i] = Stock[idx][i];
        }

        // now set up nodes again using new segment
        CurrSeg --;
        NodeSetup(CurrStep, idx);

        // cubic spline
        vector<vector<double> > result(NumOfPrice);
        result.resize(BotClip[idx]+TopClip[idx]+7);

        for (int iPrice=0; iPrice<NumOfPrice; iPrice++)
        {
            TreeInterpSpline(&*sBase, NodePrice[idx][iPrice], Stock[idx],
                             -BotClip[idx]-3, TopClip[idx]+3, iPrice,
                             result[iPrice]);
            for (i=-BotClip[idx]-3; i<=TopClip[idx]+3; i++)
                NodePrice[idx][iPrice][i] = result[iPrice][i+BotClip[idx]+3];
        }


        changed = true;
    }
    else if (CurrSeg<(int)timeLine->SegmentEnd.size()-1 && timeLine->SegmentEnd[CurrSeg]<CurrStep && RollDirection==1)
    {
        CurrSeg ++;
        // to add spline here !!!
        changed = true;
    }

    return changed;
}

/** get current tree segment */
int CTree1f::GetTreeSeg(int step)
{
    int seg;

    if (step == 0)
        seg = CurrSeg = 0;
    else if (step == timeLine->NumOfStep)
        seg = CurrSeg = timeLine->SegmentEnd.size()-1;
    else if (step == CurrStep)
        seg = CurrSeg;
    else if (step == CurrStep+1)
    {
        seg = CurrSeg;
        if (seg<(int)timeLine->SegmentEnd.size()-1 && timeLine->SegmentEnd[seg]<CurrStep+1)
            seg ++;
    }
    else
        throw ModelException("CTree1f::GetTreeSeg", "unable to locate tree segment.");

    return seg;
}

/** initialise for forward induction */
void CTree1f::InitFwdInduction(int priceStart, int priceEnd)
{
    int j, iPrice;

    for (iPrice=priceStart; iPrice<=priceEnd; iPrice++)
    {
        {// init starting point
            for (j = 0; j<NumOfInsertNode; j++)
                InsNodePrice[CurrIdx][iPrice][j] = 0.0;
            for (j = -BotClip[CurrIdx]-3; j<=TopClip[CurrIdx]+3; j++)
                NodePrice[CurrIdx][iPrice][j] = 0.0;
            NodePrice[CurrIdx][iPrice][0] = 1.0; // prob starting point
        }
    }
}

/** reset tree node boudary within original limits, used for cases like KO */
void CTree1f::ResetNodeBoundary(int step, int idx, int bot, int top)
{
    int j, iPrice;

    if (bot > BotClip[idx])
    {
        bot =(bot>BotDim[idx] ? BotDim[idx] : bot);
        for (iPrice=0; iPrice<NumOfPrice; iPrice++)
        {// prices simply copied over
            for (j=BotClip[idx]+2; j<bot+3; j++)
            {
                NodePrice[idx][iPrice][-j-1] = NodePrice[idx][iPrice][-j];
            }
        }
    }
    if (top > TopClip[idx])
    {
        top =(top>TopDim[idx] ? TopDim[idx] : top);
        for (iPrice=0; iPrice<NumOfPrice; iPrice++)
        {// prices simply copied over
            for (j=TopClip[idx]+2; j<top+3; j++)
            {
                NodePrice[idx][iPrice][j+1] = NodePrice[idx][iPrice][j];
            }
        }
    }
    // reset clipping
    BotClip[idx] = BotDim[idx];
    TopClip[idx] = TopDim[idx];

    if (bot>BotClip[idx] || top>TopClip[idx])
        NodeSetup(step, idx);
}

/** reset price array index within original limit - eg. reduce it if some price arrays are no longer need */
void CTree1f::ResetPriceArr(int step, int idx, int pStart, int pEnd)
{
    if (pStart > pEnd)
        throw ModelException("CTree1f::ResetPriceArr", "invalid parameters: pStart>pEnd");

    if (pStart<0)
        pStart = 0;
    if (pEnd >= NumOfPrice)
        pEnd = NumOfPrice - 1;

    BotArrIdx[idx] = pStart;
    TopArrIdx[idx] = pEnd;
}


/** override a control shift (eg for delta on trees)
    returns new control if constructed else returns 0 */
SensControl* CTree1f::AlterControl(const SensControl* currSensControl) const
{
    if (!Delta::TYPE->isInstance(currSensControl)
        || Maths::isZero(TreeDeltaShift /* this is not a good test */)){
        return 0;
    }

    SensControlPerName* sensControl = new Delta(TreeDeltaShift);
    sensControl->setMarketDataName(currSensControl->getMarketDataName());

    return sensControl;
}

/** simple size adjustment for smoothing delta */
void CTree1f::DeltaSizeAdjust(double shiftSize, const DateTime& valueDate,
                            const DateTime& exerciseStartDate, const DateTime& matDate,
                            bool isBermudan, bool hasDollarDivFromValueDateToMat,
                            bool isFwdStarting)
{
    const double DEFAULT_SHIFT_SIZE = 0.005;
    double length;

    bool unstableInterp = (isBermudan &&
                           (exerciseStartDate.getDate()-valueDate.getDate())
                           < .3 * (matDate.getDate() - valueDate.getDate()));

    if (unstableInterp || isFwdStarting || hasDollarDivFromValueDateToMat)
    {
        length = valueDate.yearFrac(matDate);
        length = Maths::max(1.0/12.0, Maths::min(0.5,length));

        double tweakSize = DEFAULT_SHIFT_SIZE +
            (length - 1.0/12.0) * (0.02 - DEFAULT_SHIFT_SIZE) / ( 0.5 - 1.0/12.);

        TreeDeltaShift = Maths::max(tweakSize, shiftSize);
    }
}

/** spline interpolation */
void CTree1f::TreeInterpSpline(double* s_inp, double* price_inp, double* sInterp,
                               int bot, int top, int iPrice, vector<double>& result)
{
    static const string method = "CTree1f::TreeInterpSpline";
    try {
        int n_dim = top-bot+1;

        if ((int)result.size() != n_dim) {
            result.resize(n_dim);
        }

        vector<double>::iterator arr = result.begin()-bot; // for convenience
            // spline interpolation, NR routine, with guard
        double y1 = 2e30; // default to natural spline
        double yn = 2e30;

        // trim ndim for identical s values (at bottom)
        vector<double> s, price;
        s.push_back(s_inp[bot]);
        price.push_back(price_inp[bot]);
        for (int i=bot; i<top; i++) // remove duplicate values first
        {
            if (!Maths::isZero(s_inp[i] - s_inp[i+1]))
            {
                s.push_back(s_inp[i+1]);
                price.push_back(price_inp[i+1]);
            }
        }
        n_dim = s.size();

        vector<double> y2(n_dim);
    //    double y_left, y_right, d_left, d_right;
    //    int i_guard = 0;

        // call cubic spline routines
        spline(&*s.begin()-1, &*price.begin()-1, n_dim, y1, yn,
                   &*y2.begin()-1);
            // interpolate one result at a time
        for (int j=bot; j<=top; j++)
        {
                if (sInterp[j] <= s[2] || sInterp[j] >= s[n_dim-3])        // linear beyond inner boundary
                    arr[j] = LinearTreeInterp(sInterp[j], s, price, j-bot);
                else
                {
                    splint(&*s.begin()-1, &*price.begin()-1, &*y2.begin()-1,
                           n_dim, sInterp[j], &arr[j]);

                    // safeguard for spurious results
    /* to do : test this guard
       while (i_guard<ndim && s[i_guard+1]<s[j])
       i_guard ++;
       d_left = (price[i_guard]-price[i_guard-1])/(s[i_guard]-s[i_guard-1]);
       d_right = (price[i_guard+2]-price[i_guard+1])/(s[i_guard+2]-s[i_guard+1]);
       y_left = price[i_guard] + d_left*(sInterp[j]-s[i_guard]);
       y_right = price[i_guard+1] + d_right*(sInterp[j]-s[i_guard+1]);
       if ((d_right>=d_left && (arr[j]<=y_left || arr[j]<=y_right))
       || (d_right<=d_left && (arr[j]>=y_left || arr[j]>=y_right)))
       arr[j] = TreeInterp(sInterp[j], true, iPrice, j); // simple linear interp when failed the check
    */        }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**  rebuild grid or not  */
bool CTree1f::isRebuilt(const CControl* control)
{
    bool buildTree = isTreeRebuilt(control);
    // TweakingSameGridDelta was just set in isTreeRebuilt
    sameGridTweak = TweakingSameGridDelta;
    return buildTree;
}

/** Invoked after instrument has got its market data. */
void CTree1f::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
    FDModel::getMarket(market, instruments);

    if( factors.size() != 1 )
        throw ModelException( "CTree1f::getMarket", "only 1 factor suppported" );
}

/** accept or reject factors */
bool CTree1f::acceptFactor( IMarketFactor * factor )
{
    return dynamic_cast< CAsset * >( factor ) != 0;
}

/** retrieving market data from MDF , new fd/tree interface*/
void CTree1f::retrieveFactor()
{
    static const string method = "CTree1f::retrieveFactor";
    try
    {
         // param validation
        if (StepsPerYear < -1) {
            throw ModelException(method, "invalid tree steps per year (" +
                                 Format::toString(StepsPerYear) + ")");
        }

        if (TruncationStd < 2.0 || TruncationStd > 10.0) {
            throw ModelException(method,
                                 "tree truncation (" +
                                 Format::toString(TruncationStd) +
                                 ") must be between 2.0 and 10.0");
        }

        if (TreeAlpha < 0.1 || TreeAlpha > 0.9) {
            throw ModelException(method,
                                 "tree alpha (" +
                                 Format::toString(TreeAlpha) +
                                 ") must between 0.1 and 0.9");
        }

        if (gammaNodesInterval > 1 && Maths::isNegative(gammaThreshold)){
            throw ModelException(method,
                                 "gammaThreshold [=" 
                                 + Format::toString(gammaThreshold) 
                                 + "] should be positive.");
        }

        FDModel::retrieveFactor();

        Underlier = CAssetConstSP::dynamicCast( factors[0] );
        if( ! Underlier )
            throw ModelException( "CTree1f::retrieveFactor", "only asset underlier suppported" );

        latticeProd = dynamic_cast< ILatticeProdEDR * >( prod.get() );
        if( ! latticeProd )
            throw ModelException( "Product should support ILatticeProdEDR interface" );

        if (DivAmountTreatment && !AssetUtil::isSimpleEquity(Underlier)) {
            DivAmountTreatment = false;
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/**  model initialisation, base implementation is to set up time line
    overriden FDModel method */
void CTree1f::initModel()
{
    try
    {
        arrangeDates();

        TimeMetric::SortDate(true, divCritDates);

        vector<int> density(1,1);
        if (segDates.size()>1) // this is temp, review this !!!
        {
            density.resize(this->density.size());
            for (int i=0; i<(int)density.size(); i++){
                density[i] = this->density[i]; // this is a hack
            }
        }

        Setup(getValueDate(),
                segDates,
                density,
                &critDates,
                0.0,
                equalTime,
                NumOfPrice,
                NumOfInsertNode,
                isCall,
                noExerciseWindow );

        stepsPerYearFD = StepsPerYear;
        sameGridTweak = TweakingSameGridDelta;

        validateList();

        // this is just for initialization. Later in initProd(), 
        // inst, will change via SetGammaThresholdScaled.
        gammaThresholdScaled = gammaThreshold; 
    }
    catch(exception&e){
        throw ModelException(e, "CTree1f::initModel()");
    }
}

/** add critical dates to initData for timeline setup */
void CTree1f::addDivCritDates(const DateTimeArray& divCritDates)
{
    // add to critical dates
    addCritDates( divCritDates );

    this->divCritDates.insert(this->divCritDates.end(), divCritDates.begin(), divCritDates.end() );
}

/*------------------------------------------------------------------------------
*   Description  :    roll back or forward through the tree.
*                    probability is recalculated for each node if required.
*   Return      :   price
*------------------------------------------------------------------------------*/
void CTree1f::roll(){
    static string method = "CTree1f::roll";
    try {
        // starting step
        CurrStep = timeLine->NumOfStep;
        CurrIdx = getSliceIndex(CurrStep); // using the top prod, review this

        BotClip[CurrIdx] = TopClip[CurrIdx] = BotClip[1-CurrIdx] = TopClip[1-CurrIdx] = Ndim;
        BotDim[CurrIdx] = TopDim[CurrIdx] = BotDim[1-CurrIdx] = TopDim[1-CurrIdx] = Ndim;
        TopArrIdx[0] = TopArrIdx[1] = NumOfPrice-1;
        BotArrIdx[0] = BotArrIdx[1] = 0;

        if (DivAmountTreatment) {
            /* if we're here, Underlier really is a PseudoSimpleEquity. */
            PseudoSimpleEquityConstSP pseudoSimpleEquity = PseudoSimpleEquityConstSP::dynamicCast(Underlier);
            if (!pseudoSimpleEquity) {
                throw ModelException(method, "Should never happen.");
            }
            StockFloor.resize(timeLine->StepDates.size());
            pseudoSimpleEquity->calcStockFloor(timeLine->StepDates,
                                               StockFloor);
        }
        // calc fwd for each step in one go
        Underlier->fwdValue(timeLine->StepDates, StepForward);

        // calculate centre nodes and node spacing for all tree segments
        if (BuildTree)
            CalcNodeSpacing();

        // init segmenent index
        GetTreeSeg(CurrStep);
        // set up tree nodes
        NodeSetup(CurrStep, CurrIdx);

        //set current range index
        range->setIndex( CurrIdx );
        if( NumOfInsertNode > 0 )
            insRange->setIndex( CurrIdx );

        // set current compute range
        (*range)->limits.bot1 = -BotClip[CurrIdx]-3;
        (*range)->limits.top1 = TopClip[CurrIdx]+3;

        // calcuate payoff at maturity, taking care of 3 extra nodes
        bool floorStock = DivAmountTreatment && !Maths::isZero(StockFloor[CurrIdx]);
        if (floorStock) {
            if (NumOfInsertNode>0) {
                throw ModelException(method,
                                     "Node insertion is not supported when when DivAmountTreatment == TRUE.");
            }
            for (int j = - BotClip[CurrIdx] - 3; j <= TopClip[CurrIdx] + 3; j++) {
                Stock[CurrIdx][j] += StockFloor[CurrStep];
            }
        }

        int prodIndex;
        const FDProductArray & products = getProducts();

        // for each product
        for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
            products[ prodIndex ]->preCalc( CurrStep );

        // for each product
        for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
            products[ prodIndex ]->update( CurrStep, FDProduct::BWD_T );

        if (floorStock){ // restore
            for (int j = - BotClip[CurrIdx] - 3; j <= TopClip[CurrIdx] + 3; j++) {
                Stock[CurrIdx][j] -= StockFloor[CurrStep];
            }
        }

//#ifdef NEW_TREE_TECH_TEST
        addGridMap.clear();     // initialization.

//#endif
//#ifdef NEW_TREE_TECH_TEST
        int oldBot, oldTop;

        if (useGammaNodes){
            addNodesByGamma(products, true);
        }

        // move a step to start sweeping
        CurrStep += RollDirection;
#ifdef TREE_DEBUG_FILE
//        if (control->isPricing()){
            // open file and recored the price at Maturity
        if (DEBUG_InsertWidth < -2.0){
            bool openFile = true;
            string fileName = TREE_DEBUG_FILE;
            Debug_OutPutTree(CurrStep+1, timeLine->TradeTime[CurrStep+1], Stock[CurrIdx], NodePrice[CurrIdx], -BotClip[CurrIdx]-3, TopClip[CurrIdx]+3, NumOfPrice, fileName, openFile, NumOfInsertNode, InsStock[CurrIdx], InsNodePrice[CurrIdx]);
        }

        int stopTreeAt = 11;
#endif
        // sweep the tree
        for (; CurrStep >=0 && CurrStep <timeLine->NumOfStep; CurrStep += RollDirection)
        {            
            CurrIdx = 1- CurrIdx; // swap array index

#ifdef TREE_DEBUG_FILE
            if (stopTreeAt == CurrStep)
                double test = 0.0;
#endif
            // clear inserted node
            for (int j=0; j<NumOfInsertNode; j++)
                SetInsertNode(CurrIdx, j, -1.0, 0);

            NodeSetup(CurrStep, CurrIdx);
            if (useGammaNodes){
                oldBot = (*range)->limits.bot1;
                oldTop = (*range)->limits.top1;
            }

            //set current range index
            range->setIndex( CurrIdx );
            if( NumOfInsertNode > 0 )
                insRange->setIndex( CurrIdx );

            // set current compute range
            (*range)->limits.bot1 = -BotClip[CurrIdx]-3;
            (*range)->limits.top1 = TopClip[CurrIdx]+3;

            // for each product
            for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
                products[ prodIndex ]->preCalc( CurrStep );

            if (useGammaNodes)
                RollTree(&drift_GN);  //need to store the drift.
            else
                RollTree();

            
#ifdef  TREE_THETA_CAP
            // theta smoothing
            // putting here to avoid instrument's payoff...  correct?
            if (useThetaCap){
                addThetaCap();
            }
#endif

            if (RollDirection == -1)
            {// calcuate early exercise etc. taking care of 3 extra nodes
                bool floorStock = DivAmountTreatment && !Maths::isZero(StockFloor[CurrIdx]);
                if (floorStock) {
                    if (NumOfInsertNode>0) {
                        throw ModelException(method,
                                             "Node insertion is not supported when when DivAmountTreatment == TRUE.");
                    }
                    for (int j = - BotClip[CurrIdx] - 3; j <= TopClip[CurrIdx] + 3; j++) {
                        Stock[CurrIdx][j] += StockFloor[CurrStep];
                    }
                }

                // for each product
                for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
                    products[ prodIndex ]->update( CurrStep, FDProduct::BWD );

                if (floorStock) {// restore
                    for (int j = - BotClip[CurrIdx] - 3; j <= TopClip[CurrIdx] + 3; j++) {
                        Stock[CurrIdx][j] -= StockFloor[CurrStep];
                    }
                }
            }
            // now, I can move this after UpDate.
            if (useGammaNodes){
                addNodesByGamma(products);
            }

#ifdef TREE_DEBUG_FILE
//            if (CurrStep == 0){
            if (DEBUG_InsertWidth < -2.0){
                // only if fileName does not exists, we output new file, otherwise no new output
                string fileName = TREE_DEBUG_FILE;
                Debug_OutPutTree(CurrStep, timeLine->TradeTime[CurrStep], Stock[CurrIdx], NodePrice[CurrIdx], -BotClip[CurrIdx]-3, TopClip[CurrIdx]+3, NumOfPrice, fileName, false, NumOfInsertNode, InsStock[CurrIdx], InsNodePrice[CurrIdx]);
            }
#endif

            // switching segment if needed
            SwitchSegment();
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**  model parameter validation */
void CTree1f::finaliseModel(CControl*    control)
{
    static const string method = "CTree1f::finaliseModel";
    try
    {
        if( ! payoffIndex )
            throw ModelException( method, "no spot index was requested" );

        // map stock to spot
        TreeSliceGeneral & stock = static_cast< TreeSliceGeneral & >( *payoffIndex->getSlicesToDEV()[0] );
        // perform mapping for nodes. new interface products to old tree nodes
        const vector< TreeSliceSP > & nodePrice = prod->getSlicesToDEV();
        ASSERT( NumOfPrice <= (int)nodePrice.size() );
        for( int i = 0; i < 2; ++i )
        {
            Stock[i] = stock[i];
            NodePrice[i].resize( NumOfPrice );
            for( int j = 0; j < NumOfPrice; ++j )
                NodePrice[i][j] = static_cast< TreeSliceGeneral & >( *nodePrice[j] )[i];
        }

        if (NumOfInsertNode>0)
        {
            latticeProd->mapInsertNode(insRange, insStock, insNodePrice, InsPriority);
            for( int i = 0; i < 2; ++i ) // for Current and last (1-Current) array
            {
                InsStock[i] = (*insStock)[i];
                InsNodePrice[i].resize( NumOfPrice );
                for( int j = 0; j < NumOfPrice; ++j )
                    InsNodePrice[i][j] = (*insNodePrice)[i][j];
            }
        }

        // use it to init pseudo equity
        if (!BuildTree) {
            if( DivAmountTreatment) {
                Underlier = CAssetConstSP(
                    PseudoSimpleEquity::create(Underlier.get(),
                                               divCritDates,
                                               isCall,
                                               noExerciseWindow));
            }
            InitVol();
        }

        PostSetup();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double CTree1f::getPrice0( const TreeSlice & price ) const
{
    double * s = payoffIndex->getValue( 0 ).getValues();
    double * p = price.getValues();

    // if forward starting then adjust the spot
    double s0;
    if( getValueDate() < prod->getStartDate() )
        s0 = Underlier->fwdValue(prod->getStartDate());
    else
        s0 = Underlier->getSpot();

    return QuadraticInterp( s0,
        s[-1],
        s[0],
        s[1],
        p[-1],
        p[0],
        p[1]);
}

/** get the initial conditions and corresponding indexes in the slices
    for forward induction */

void CTree1f::getInitialConditions(IntArray& initialIndex,
                                   DoubleArray& initialValue) const{
}

FDProductSP CTree1f::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            return payoffIndex = FDModel::makeProduct( IProdCreatorSP( new IndexSpecEDR( spec->getFactor()->getName() ) ) );
        }
        else
        {
            throw ModelException( "CTree1f::makeProduct",
                "IndexSpec source is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' is expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

// adding nodes.  
// 1. call calcGammAndGridMap to know where to add new nodes.
// 2. store the all stocks & prices.
// 3. calulate new additional nodes and prices (interpolated by t+1).
// 4. call products->update.
// 5. re-arranging (combining) new nodes & already exising nodes.
void CTree1f::addNodesByGamma(const FDProductArray & products, bool isMat)
{
    static const string method = "CTree1f::addNodesByGamma";
    try
    {
        int i;
        int step = CurrStep;
        // 1. call calcGammAndGridMap to know where to add new nodes.
        if (CacheMode == USE_CACHE){            
            addGridMap = gammaNodeMaps[step].getAddGridMap();
            //addGridMap = addGridMap_cache[CurrStep];
        }
        else{
            calcGammaAndGridMap( CurrStep, addGridMap);
            if (CacheMode == MAKE_CACHE){   
                gammaNodeMaps[step].setAddGridMap(addGridMap);
            }
        }
        int numAdded = addGridMap.size();
        vector<double> addStock;

        if (numAdded>0){
            int j,iPrice;
            int oldBot = (*range)->limits.bot1;
            int oldTop = (*range)->limits.top1;

            // 2. store the all stocks & prices.
            // store the current slices to orginal grids.
            gammaNodeMaps[step].storeStockAndPrice(oldBot,oldTop,NumOfPrice,Stock[CurrIdx],&NodePrice[CurrIdx][0]);

            // 3. calulate new additional nodes levels and prices (interpolated by t+1).
            vector<vector<double> > addPrices;
            addPrices.resize(NumOfPrice);
            double df = (isMat ? 1.0 : discYC->pv(timeLine->StepDates[CurrStep], 
                                                  timeLine->StepDates[CurrStep+1]));
            double drft = 0.0;
            if (!isMat){
                if (drift_GN.size() == 1){
                    drft = drift_GN[0];
                }
                else if ((int)drift_GN.size() != TopClip[CurrIdx]+BotClip[CurrIdx]+1){
                    throw ModelException(method, "Internal Error.  drift_DN size["
                                                + Format::toString(drift_GN.size())
                                                +"]doesn't match to stock nodes size["
                                                + Format::toString(oldTop-oldBot+1)
                                                +"].");
                }
            }
            j=-BotClip[CurrIdx];   // search from bottom+3.  s must be increasing order!!
            for (i=0;i<(int)addGridMap.size();i++){
                double s = Stock[CurrIdx][0] * exp(NodeSpace[GetTreeSeg(CurrStep)]*(double)addGridMap[i]);
                // validation
                if (s<Stock[CurrIdx][-BotClip[CurrIdx]-1] ||
                    s>Stock[CurrIdx][TopClip[CurrIdx]+1]){
                    throw ModelException(method, "Internal Error.  Trying to add Nodes at "
                        + Format::toString(s)
                        + ", which is lower than BotClip+1 = "
                        + Format::toString(Stock[CurrIdx][-BotClip[CurrIdx]-1])
                        + " or higher than TopClip+1 = "
                        + Format::toString(Stock[CurrIdx][TopClip[CurrIdx]+1])
                        + ".");
                }
                addStock.push_back(s);
                double s_next = s;
                if (!isMat){
                    if (drift_GN.size() > 1){
                        if (s<Stock[CurrIdx][j]){
                            // between Stock[bot+2] and BotClip.  Use BotClip's drift.
                            drft = drift_GN[0]; 
                        }
                        else if (s>Stock[CurrIdx][TopClip[CurrIdx]]){
                            // between Stock[top-2] and TopClip.  Use TopClip's drift.
                            drft = drift_GN.back(); 
                        }
                        else{
                            while (j<TopClip[CurrIdx]){
                                if(Stock[CurrIdx][j] < s && s<Stock[CurrIdx][j+1])
                                    break;
                                else
                                    j++;
                            }
                            int j_drift = j+BotClip[CurrIdx];
                            ASSERT(0 <= j_drift && j_drift < (int)drift_GN.size());
                            drft = LinearInterp(s, Stock[CurrIdx][j],Stock[CurrIdx][j+1],
                                                   drift_GN[j_drift],drift_GN[j_drift+1]);
                        }
                    }
                    s_next *= drft;    
                }
                for (iPrice = 0; iPrice<NumOfPrice; iPrice++){
                    double p = df * TreeInterp(s_next, isMat, iPrice, addGridMap[i]/gammaNodesInterval);
                    addPrices[iPrice].push_back(p);
                }
            }
            
            // 4. call products->update.  
            // copy those new nodes to slice, and update.            
            (*range)->limits.bot1 = -NLimit-3; 
            (*range)->limits.top1 = -NLimit-3+numAdded-1; 
            for (i=0;i<numAdded;i++){
                Stock[CurrIdx][-NLimit-3+i] = addStock[i];
                for (iPrice=0; iPrice<NumOfPrice; iPrice++){
                    NodePrice[CurrIdx][iPrice][-NLimit-3+i] = addPrices[iPrice][i];
                }
            }

            // overwrite the prices by using instrument.
            for (int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
                products[ prodIndex ]->update( CurrStep, isMat ? FDProduct::BWD_T:FDProduct::BWD);
            
            // get out the update results
            for (i=0;i<numAdded;i++){
                for (iPrice=0; iPrice<NumOfPrice; iPrice++){
                    addPrices[iPrice][i] = NodePrice[CurrIdx][iPrice][-NLimit-3+i];
                }
            }

            // 5. re-arranging (combining) new nodes & already exising nodes.
            // put everything into tree slice
            // set orgIdx (re-new)
            if (CacheMode != USE_CACHE)
                gammaNodeMaps[step].updateIdxAfter();
            gammaNodeMaps[step].integrate(numAdded, NumOfPrice, addStock, addPrices, 
                                          Stock[CurrIdx],NodePrice[CurrIdx]);
            
            // check centre node
            if (fabs(Stock[CurrIdx][0] - CentreNode[CurrStep])>1.0e-5)
                throw ModelException(method, "Internal Error.  After added Nodes, centre nodes are away from CentreNode at ["
                                              +Format::toString(CurrStep)+ "].");
            (*range)->limits.bot1 = gammaNodeMaps[step].getBotAfter();
            (*range)->limits.top1 = gammaNodeMaps[step].getTopAfter();
            TopClip[CurrIdx] = (*range)->limits.top1-3;
            BotClip[CurrIdx] = -((*range)->limits.bot1+3);

            // InterpTreeBoundary was called before adding nodes, 
            // so the -BotClip-1 / TopClip+1 node could be poor.  
            // Here, recacluate them, if the nodes are added inside of -BotClip-1 or TopClip+1.  
            if (!isMat){
                const vector<int>& prevSt = gammaNodeMaps[step+1].getStatus();
                const vector<int>& currSt = gammaNodeMaps[step].getStatus();
                ASSERT(prevSt.size()>=7 && currSt.size()>=7);
                if (  (prevSt[3]<0 && currSt[3] >=0)
                    ||(prevSt[prevSt.size()-4]<0 && currSt[currSt.size()-4] >=0) ){                    
                    InterpTreeBoundary();
                }
            }

        }else{
            // need to update the idxAfter.
            gammaNodeMaps[step].updateIdxAfter();
        }

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// set up gammaThresholdScaled.  
void CTree1f::SetGammaThresholdScaled(double refStrike, double notional){
        gammaThresholdScaled = gammaThreshold * notional / refStrike / (0.01 * refStrike);
}

#ifdef  TREE_THETA_CAP
void CTree1f::ScaleThetaCapThreshold(double refStrike, double notional){
        threholdTC = ThetaCapThrehold * notional / refStrike ;
}

// activate theta smooth (via instrument).
bool CTree1f::activateThetaCap(bool isPositiveThetaCap, double ThetaCapThrehold)
{
    if (!Maths::isZero(ThetaCapThrehold))
        useThetaCap = true;
    else
        useThetaCap = false;
    this->ThetaCapThrehold = ThetaCapThrehold;
    this->isPositiveThetaCap = isPositiveThetaCap;
    return useThetaCap;
}
#endif

// calculate gamma on each tree nodes and return addGridMap.  
// addGridMap will contain index array in original (big) measure, whose are going to be added.  
// in this function, addMap is also updated.
bool CTree1f::calcGammaAndGridMap(int step,vector<int> & addGridMap)
{
    bool useFineGrid = false;
    if (useGammaNodes){
        const TreeSlice & slice = payoffIndex->getValue( step );
        int bot, top;
        slice.getCalcRange( bot, top );
        double * s = slice.getValues();

        int idx = CurrIdx; // using current Node's Price

        int addBotLimit;
        int addTopLimit;
        gammaNodeMaps[step].getOrgIdxOfAddLimit(addBotLimit, addTopLimit);

        double dSu, dSd, dSm;
        double gamma = 0.0;
        double tmpGamma = 0.0;
        bool wasBigGamma = false;

        // j_buff is to avoid the setting additional nodes around floor/ceiling boudary in tree.
        // tree usually have 3 nodes in boudary area.  Thus, 3 means avoid adding nodes around
        // those boundary nodes + lowest nodes in normal tree grids.
        // Also, to avoid trouble at shirinking time steps, adding node will start from 
        // gammaNodesInterval-1.
        // (NB) Additional node will be added below of the "addGridMap" nodes.
        //int j_buff = 3+tree1f->gammaNodesInterval-1;
        int j_buff = 3;

        addGridMap.clear();
        double peak = gammaThresholdScaled;
        //for (int i=0; i<tree1f->NumOfPrice; i++){
        for (int j=bot+j_buff;j<=top-j_buff;j++){ 
            if (gammaNodeMaps[step].getStByBef(j) == 0){
                //addMaps[CurrStep][orgIdx[j]] == 0){
                gamma = tmpGamma = 0.0;
                dSu = s[j+1]-s[j];
                dSd = s[j]-s[j-1];
                dSm = (s[j+1]-s[j-1])/2.0;
                for (int iPrice=0; iPrice<NumOfPrice; iPrice++){
                    // Use the biggest gamma among many time slices.
                    // gamma = {delta(+)-delta(-)} / {(S(+)+S(0))/2 - (S(0)-S(-))/2}
                    tmpGamma = fabs( ( ((NodePrice[idx][iPrice])[j+1]-(NodePrice[idx][iPrice])[j])/dSu
                                      -((NodePrice[idx][iPrice])[j]-(NodePrice[idx][iPrice])[j-1])/dSd) /dSm );
                    if (tmpGamma>gamma)
                        gamma = tmpGamma;
                }
                if (gamma > peak || wasBigGamma) {
                    int orgIdx = gammaNodeMaps[step].getOrgIdxByBef(j);
                    int l = orgIdx - gammaNodesInterval;
                    for (int k=1;k<gammaNodesInterval;k++){
                        l++; 
                        while (l<=addBotLimit){
                            k++,l++;
                        }
                        if (gammaNodeMaps[step].getStByOrg(l) < 0){
                            gammaNodeMaps[step].turnOnByOrg(l);
                            addGridMap.push_back(l);
                        }
                    }
                    if (j==top-j_buff){
                        // add node between top - j_buff to next.  
                        // l<=NLimit is to avoid adding node beyond original node.
                        for (int k=1;k<gammaNodesInterval && l<addTopLimit;k++){                                 
                            l = orgIdx+k;    // turn "ON" for up-side of top-j_buff.
                            if (gammaNodeMaps[step].getStByOrg(l) < 0){
                                gammaNodeMaps[step].turnOnByOrg(l);
                                addGridMap.push_back(l);
                            }
                        }
                    }
                    wasBigGamma = (fabs(gamma) > peak);
                    ASSERT(-NLimit-3<=l && l<=NLimit+3);    // check boundary
                }
            }
        }
        if (addGridMap.size()>0)
            useFineGrid = true;
    }

    return useFineGrid;
}


// class for Gamma Node
CTree1f::GammaNodeMap::GammaNodeMap():baseOrg(0),baseBef(0),baseAft(0){
    orgIdx.clear();
    idxBefore.clear();
    idxAfter.clear();
    status.clear();
    addGridMap.clear();
    S_store.clear();
    NP_store.clear();
}

void CTree1f::GammaNodeMap::initIdx(int nSize){
        int orgSize = nSize;
        orgIdx.clear();
        idxBefore.clear();
        status.clear();
        orgIdx.resize(orgSize, NOT_USED_IDX);
        status.resize(orgSize,-1);        
}

// bot&top is the lowest & highest node in all nodes.
void CTree1f::GammaNodeMap::SetAtMat(int NLimit, int gammaInterval){
    static const string method = "CTree1f::GammaNodeMap::SetAtMat";
    try
    {
        ASSERT(gammaInterval>0);
        initIdx(2*NLimit+7);       
        int orgBot = -(NLimit+3);
        int orgSize = orgIdx.size();
        for (int i=0; i<orgSize; i++){
            orgIdx[i] = i + orgBot;
            if (i<=2 || i>=orgSize-3 || orgIdx[i]%gammaInterval ==0){
                idxBefore.push_back(orgIdx[i]);
                status[i] = 0;
            }
        }
        baseBef = calcLowerGrid(idxBefore);
        baseOrg = calcLowerGrid(orgIdx); 
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// set up GammaNodeMap using topOrg.  
// return bot & top of idxBefore, which is set up.  
void CTree1f::GammaNodeMap::Set(const int botOrg, const int topOrg, const GammaNodeMap &gnmPre){
    static const string method = "CTree1f::GammaNodeMap::Set";
    try
    {
        int orgSize = topOrg-botOrg+1;
        initIdx(orgSize);
        for (int i=0; i<orgSize; i++){
            orgIdx[i] = i+botOrg;
        }

        const vector<int>& pre_orgIdx = gnmPre.getOrgIdx();
        const vector<int>& pre_status = gnmPre.getStatus();
        
        int pre_bot = pre_orgIdx[0];
        int pre_top = pre_orgIdx.back();
        
        // diff_bot is supposed to be positive, becuse botOrg & pre_bot are negative number and shirinking.
        int diff_bot = botOrg - pre_bot;    
        
        // update status.  
        for (int i=0; i<orgSize;i++){
            if (pre_status[i+diff_bot] == 2)
                status[i] = 1;
            else
                status[i] = pre_status[i+diff_bot];
        }
        // Force that top-2 & bot+2 node is always existing.
        status[2] = status[orgSize-3] = 0;

        // update idxBefore.  
        int bsize = calcNumOfActiveNodes();
        idxBefore.resize(bsize,NOT_USED_IDX);        
        for (int i=0, j=0; i<orgSize;i++){
            if (status[i]>=0){
                ASSERT(j<bsize);
                idxBefore[j] = orgIdx[i];
                j++;
            }
        }

        baseBef = calcLowerGrid(idxBefore);
        baseOrg = calcLowerGrid(orgIdx);         
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CTree1f::GammaNodeMap::storeStockAndPrice(int oldBot,int oldTop,int NumOfPrice,double* s, double** price){
    // initialization
    S_store.clear();
    NP_store.clear();
    int nSize = oldTop-oldBot+1;
    S_store.resize(nSize);
    NP_store.resize(NumOfPrice);
    for (int i=0;i<NumOfPrice;i++)
        NP_store[i].resize(nSize);
    // copy
    for (int i=oldBot; i<=oldTop; i++){
        S_store[i-oldBot] = s[i];
        for (int iPrice=0;iPrice<(int)NP_store.size();iPrice++)
            NP_store[iPrice][i-oldBot] = price[iPrice][i];
    }    
}

void CTree1f::GammaNodeMap::integrate(int numAdded, int NumOfPrice,
                                      const vector<double>& addStock, 
                                      const vector<vector<double> >& addPrice, 
                                      double* s,                     //output
                                      vector< double * > & price){    //output
    static const string method = "CTree1f::GammaNodeMap::integrate";
    try
    {
        int iOrg;             // index to access original idx
        int iNew;             // index to express new grid index.
        int j=0;    // index for new Array.  
        int k=0;    // index for additional node.
        int botOrg = orgIdx[0];
        int topOrg = orgIdx.back();
        int sizeAdd = addStock.size();
        int sizeStore = S_store.size();       
        if ((int)idxAfter.size() != sizeAdd + sizeStore){
            throw ModelException(method, "Internal Error.  IdxArray Size["
                                        + Format::toString(idxAfter.size())
                                        +"] is not equal sum of addStock size["
                                        + Format::toString(sizeAdd)
                                        +"] and S_store["
                                        + Format::toString(sizeStore)
                                        +"]. ");
        }
        for (int i=0;i<(int)idxAfter.size();i++){
            iOrg = idxAfter[i] + baseOrg;
            iNew = i-baseAft;
            ASSERT(0<=iOrg && iOrg<(int)orgIdx.size());
            ASSERT(botOrg<=iNew && iNew <=topOrg);
            if (status[iOrg] == 2){
                ASSERT(k<sizeAdd);
                s[iNew] = addStock[k];
                for (int iPrice=0;iPrice<NumOfPrice;iPrice++){
                    price[iPrice][iNew] = addPrice[iPrice][k];
                }
                k++;                
            }
            else if (status[iOrg] >= 0){
                ASSERT(0<=j && j<sizeStore);
                s[iNew] = S_store[j];
                for (int iPrice=0;iPrice<NumOfPrice;iPrice++){
                    price[iPrice][iNew] = NP_store[iPrice][j];
                }
                j++;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CTree1f::GammaNodeMap::updateIdxAfter(){
    static const string method = "CTree1f::GammaNodeMap::updateIdxAfter";
    try
    {
        int asize = calcNumOfActiveNodes();
        idxAfter.resize(asize,NOT_USED_IDX);        
        for (int i=0, j=0;i<(int)orgIdx.size();i++){
            if (status[i]>=0){
                ASSERT(j<asize);
                idxAfter[j] = orgIdx[i];
                j++;
            }
        }
        baseAft = calcLowerGrid(idxAfter);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

int CTree1f::GammaNodeMap::calcNumOfActiveNodes() const{
    int num=0;
    for (int i=0; i<(int)status.size(); i++){
        if(status[i]>=0)
            num++;
    }
    return num;
}

int CTree1f::GammaNodeMap::calcLowerGrid(const vector<int> &idxArray) const{
    int i=0, num=0;
    while (idxArray[i]<0){
        i++, num++;
    }
    return num;
}

void CTree1f::GammaNodeMap::turnOnByOrg(int idx){
    ASSERT(0<=idx+baseOrg && idx+baseOrg<(int)status.size());
    status[idx+baseOrg] = 2;
}

const vector<int>& CTree1f::GammaNodeMap::getStatus() const{
    return status;
}

const vector<int>& CTree1f::GammaNodeMap::getOrgIdx() const{
    return orgIdx;
}

int CTree1f::GammaNodeMap::getBotBefore()  const{
    return -baseBef;
}

int CTree1f::GammaNodeMap::getTopBefore() const{
    return idxBefore.size()-baseBef-1;
}

int CTree1f::GammaNodeMap::getBotAfter() const{
    return -baseAft;
}

int CTree1f::GammaNodeMap::getTopAfter() const{
    return idxAfter.size()-baseAft-1;
}

vector<int> CTree1f::GammaNodeMap::getAddGridMap() const{
    return addGridMap;
}
void CTree1f::GammaNodeMap::setAddGridMap(vector<int> new_addGridMap){
    addGridMap = new_addGridMap;
}

void CTree1f::GammaNodeMap::getOrgIdxOfAddLimit(int &addBotLimit, int &addTopLimit) const{
    addBotLimit = idxBefore[2];
    addTopLimit = idxBefore[idxBefore.size()-3];
}
DRLIB_END_NAMESPACE
