//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fCEVJ.cpp
//
//   Description : one factor trinomial tree for CEV+Jump process.
//
//----------------------------------------------------------------------------
//#define DEBUG_CEVJ_TREE

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Tree1fCEVJ.hpp"
#include "edginc/CEVJ.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/LatticeProdEDR.hpp"

#ifdef DEBUG_CEVJ_TREE
static FILE *DumpFile = NULL; // static file handle for dump file (keiji)
static int Dumpfilenum = 0;
static int iOut = 0;
#endif

DRLIB_BEGIN_NAMESPACE

//for reflection
CTree1fCEVJ::CTree1fCEVJ(CClassConstSP clazz):CTree1f(clazz)
{
    DEBUG_BSSmooth = 1;    //1 is "ON", else "Off"
    DEBUG_JumpFrom = -1;     //0 is jump at t, 1 is at t+1, else "Off"
    DEBUG_UseCtrlVar = false;
    IsStoreJumpNode = false;
    SetDefault();
    StepsPerYear = -1;  //avoid to leave 0 
}

CTree1fCEVJ::CTree1fCEVJ():CTree1f(TYPE)
{
    DEBUG_BSSmooth = 1;    //1 is "ON", else "Off"
    DEBUG_JumpFrom = -1;     //0 is jump at t, 1 is at t+1, else "Off"
    DEBUG_UseCtrlVar = false;
    IsStoreJumpNode = false;
    SetDefault();
    StepsPerYear = -1;  //avoid to leave 0 
}

CTree1fCEVJ::~CTree1fCEVJ()
{
    Clear();
}

/** clean up */
void CTree1fCEVJ::Clear()
{
    int i, j;
    if (JumpPrice[0] != 0)
    {
        for (i=0; i<2; i++)
        {
            for (j=0; j<NumOfPrice; j++)
                delete [] (JumpPrice[i][j] -= (Ndim+3));
            delete [] JumpPrice[i];
        }
    }
    if (InsJumpPrice[0] != 0)
    {
        for (i=0; i<2; i++)
        {
            for (j=0; j<NumOfPrice; j++)
                delete [] (InsJumpPrice[i][j]);
            delete [] (InsJumpPrice[i]);
        }        
    }

    SetDefault();
}
void CTree1fCEVJ::SetDefault()
{
    StepJumpRate = 0.0;
    JumpPrice[0] = JumpPrice[1] = 0;
    InsJumpPrice[0] = InsJumpPrice[1] = 0;
    PerNodeProb = true;
}

/** get processed vol*/
// need to examine if this is general
void CTree1fCEVJ::InitVol()
{
    // no vol request used - since an underlier can only hold one type of vol, request is redundent here
    VolCEVJ = CEVJProcessedSP::dynamicCast(
        CVolProcessedSP(Underlier->getProcessedVol((CVolRequest*)0)));
    // get time metric
    timeMetric = VolCEVJ->GetTimeMetric();

    if (DEBUG_UseCtrlVar)
        throw ModelException("CTree1fCEVJ", "CEVJ model does not support control variate technique.");

    // get vol and corr for quanto or struck
    CAssetConstSP plainAsset = Underlier;;

    if (StruckEquity::TYPE->isInstance(Underlier) && !(AssetUtil::isBasket(Underlier)))
    {
        plainAsset = StruckEquityConstSP::dynamicCast(Underlier)->getPlainAsset();
        fxAsset = StruckEquityConstSP::dynamicCast(Underlier)->getFX();
        eqFXCorr = StruckEquityConstSP::dynamicCast(Underlier)->getCorrelation();

        if (!fxAsset){
            throw ModelException("CTree1fLV::InitVol", "NULL fx asset");
        }
        if (!eqFXCorr){
            throw ModelException("CTree1fLV::InitVol", "NULL correlation");
        }
        // create an atm vol request for the fx side of things
        CVolRequestSP atmInterp(new ATMVolRequest());
        CVolProcessedSP interpFXVol(fxAsset->getProcessedVol(atmInterp.get()));
        // cast to the type of vol we need
        volFXBS = CVolProcessedBSSP::dynamicCast(interpFXVol);
    }
    else if (ProtEquity::TYPE->isInstance(Underlier) && !(AssetUtil::isBasket(Underlier)))
    {
        plainAsset = ProtEquityConstSP::dynamicCast(Underlier)->getPlainAsset();
        FXVolBaseWrapper fxVol = ProtEquityConstSP::dynamicCast(Underlier)->getFXVol();
        eqFXCorr = ProtEquityConstSP::dynamicCast(Underlier)->getCorrelation();
        /* interpolate the fx vols atm */
        ATMVolRequestSP fxVolRequest(new ATMVolRequest());
        // interpolate the vol - need asset to define ATM
        const Asset* fx = ProtEquityConstSP::dynamicCast(Underlier)->getFXAsset();
        CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), fx));
        // cast to the type of vol we're expecting
        volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
        Underlier = plainAsset; // assign it to plain asset for the tree
    }
}

MarketDataFetcherSP CTree1fCEVJ::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(CEVJ::TYPE->getName(),
                                               VolSurface::TYPE->getName()));   
}

IModel::WantsRiskMapping CTree1fCEVJ::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** calculate a term structure vol^2 or just one point at maturity 
    This vol^2 comes from Diffusion Vol only and not depend on StockLevel
    Thus, the v2 returned by here is assumed to only use for Build Tree, not RollBack.*/
void CTree1fCEVJ::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
                           const DateTime& matDate, CTermStructure& v2)
{
    static const string method = "CTree1fCEVJ::CalcV2Term";
    try {
        int i;
        DoubleArray vol_sq;
        
        // write down the bench marks
        const DateTimeArray& bmDates = VolCEVJ->GetBMDates();

        // retrieve variance for each date from the processed vol object
        vector<double> yrs;
        vol_sq.resize(bmDates.size());
        yrs.resize(bmDates.size());
        
        const DoubleArray& bmVol = VolCEVJ->GetParamArr("ATM_VOL");
        const DoubleArray& bmCEVPower = VolCEVJ->GetParamArr("CEV_POWER");
        for (i=0; i<vol_sq.size(); i++)
        {
            vol_sq[i] = bmVol[i]*bmVol[i]*::pow(VolCEVJ->GetSpotRef(), 2.0*(1.0-bmCEVPower[i]));
            yrs[i] = VolCEVJ->calcTradingTime(valDate, bmDates[i]);
            /* for Struck case adjustment : Need to be more thought
               Only Diffusion Vol are added adjust ment for building tree. */
            if (prod->getCcyTreatment() =="S")
            {
                double vol_fx = volFXBS->CalcVol(valDate, bmDates[i]);
                double corr = eqFXCorr->getCorrelation();
                double vol_eq = sqrt(vol_sq[i]);
                vol_sq[i] += vol_fx*vol_fx + 2.0*corr*vol_eq*vol_fx;
            }
        }

        // populate it to a vol^2 term structure
        v2.Populate(valDate, vol_sq.size(), bmDates.begin(), yrs.begin(), 
                    vol_sq.begin());

        // populate CEVJ params using year fractions.
        termDiffV2.Populate(valDate, vol_sq.size(), bmDates.begin(),
                            yrs.begin(), vol_sq.begin());
        termCEVPower.Populate(valDate, yrs.size(), bmDates.begin(), 
                              yrs.begin(),
                              VolCEVJ->GetParamArr("CEV_POWER").begin());
        termJumpRate.Populate(valDate, yrs.size(), bmDates.begin(), 
                              yrs.begin(),
                              VolCEVJ->GetParamArr("JUMP_RATE").begin());
        termJumpMean.Populate(valDate, yrs.size(), bmDates.begin(), 
                              yrs.begin(),
                              VolCEVJ->GetParamArr("JUMP_MEAN").begin());
        termJumpWidth.Populate(valDate, yrs.size(), bmDates.begin(), 
                               yrs.begin(),
                               VolCEVJ->GetParamArr("JUMP_WIDTH").begin());

        // For store purpose
        ValDate = valDate;
        MatDate = matDate;
        StartDate = startDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** set up variance array */
/*  I guess it would be better to build Variance depend on Stocks value. (keiji)*/
void CTree1fCEVJ::PostSetup()
{
    int i;
    Clear();
    // mem for jump price array
    for (i=0; i<2; i++)
    {// for Current and last (1-Current) array
        JumpPrice[i] = new double*[NumOfPrice];
        InsJumpPrice[i] = new double*[NumOfPrice];
        for (int j=0; j<NumOfPrice; j++)
        {
            JumpPrice[i][j] = new double[2*Ndim+7] +Ndim+3;
            InsJumpPrice[i][j] = new double[NumOfInsertNode+1];
        }
    }

    // Remake here for greeks.
    CTermStructure v_term; 
    CalcV2Term(ValDate, StartDate, MatDate, v_term);

    // set variance to each step
    Variance.resize(timeLine->NumOfStep+1);
    Variance[0] = 0.0; // start with 0
    for (i = 1; i <= timeLine->NumOfStep; i++){
        Variance[i] =  timeLine->TradeTime[i] * termDiffV2.InterpFlatFwd(0, timeLine->TradeTime[i]);
    }
}

/** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
    they must be either one per step or one per node (from -BotClip to TopClip)
    returns true if constant probabilities can be used for all nodes (fast roll) 
    For Jump Case, JumpNode can be store for pricing later. 
    NB!!!   Variance are only for diffusion part.  Lacking Jump Diffusion!! */
bool CTree1fCEVJ::CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& var_dt, vector<double>* drift)
{
    DateTime thisStepDate, nextStepDate;    
    static const string method = "CTree1fCEVJ::CalcStepDriftAndVar";
    try 
    {
        int step = CurrStep;
        if (CurrStep == timeLine->NumOfStep)
            step = CurrStep -1;     // return same vol between T-1 to T, e.g. for Barrier Adjust 
        //ASSERT(CurrStep < timeLine->NumOfStep);
        if (!drift)
            throw ModelException("CTree1fCEVJ::CalcStepDriftAndVar", "must supply valid array for drift.");

        thisStepDate = timeLine->StepDates[step];
        nextStepDate = timeLine->StepDates[step+1];
        var_dt.resize(end-start+1);
        drift->resize(end-start+1);

        int i;
        double drift_fwd = StepForward[step+1]/StepForward[step];
        StepJumpRate = VolCEVJ->CalcStepJumpRate(thisStepDate, nextStepDate);

        if (StepJumpRate == 0)        
        {// No Jump Case
            for (i=0; i<=end-start; i++)
                (*drift)[i] = drift_fwd;
        }
        else
        {// With Jump
            VolCEVJ->SetVarCalc(false);            
            vector<double> jump_drift(end-start+1);
            //avoid calc outside Tree.
            int bot=start, top=end;
            if (start < -BotClip[CurrIdx])
                bot = -BotClip[CurrIdx];
            if (end > TopClip[CurrIdx])
                top = TopClip[CurrIdx];
            if (IsStoreJumpNode){
                VolCEVJ->CalcJumpDrift(s, thisStepDate, jump_drift, JumpDownNode, JumpUpNode, bot, top);
                IsStoreJumpNode = false;
            }
            else{
                vector<double> jump_down;
                vector<double> jump_up;
                VolCEVJ->CalcJumpDrift(s, thisStepDate, jump_drift, jump_down, jump_up, bot, top);
            }

            for (i=0; i<=end-start; i++)
            {
                if (i+start < bot || i+start > top) //for boundary area of Tree.
                    (*drift)[i] = drift_fwd;
                else if (DEBUG_JumpFrom == 0 || DEBUG_JumpFrom == 1)
                    (*drift)[i] = drift_fwd / (1.0 - (-jump_drift[i])*StepJumpRate);
                else    //Re-normalized diffusion drift.
                    (*drift)[i] = 1.0 + (drift_fwd -1.0 - jump_drift[i] * StepJumpRate)/(1.0 - StepJumpRate);
            }
        }

        double cevpower = termCEVPower.InterpLinear(timeLine->TradeTime[step]);
        double dt = timeLine->TradeYrFrac[step+1];
        double var_diff = Variance[step+1] - Variance[step];

        if (DEBUG_JumpFrom == 0 || DEBUG_JumpFrom == 1)
            VolCEVJ->CalcLocVar(s, (*drift), cevpower, 0.0, dt, var_diff ,start, end, var_dt);  // No renormalize for variance
        else
            VolCEVJ->CalcLocVar(s, (*drift), cevpower, StepJumpRate, dt, var_diff ,start, end, var_dt);

        //quanto and struck adjustments 
        // no specific dollar div treatment
        
        string ccyTreatment = prod->getCcyTreatment();

        if ( ccyTreatment =="S" || ccyTreatment =="P")
        {
            double vol;
            double jump_vol, jump_volsq, jump_mean, jump_meansq;
            double vol_fx = volFXBS->CalcVol(timeLine->StepDates[step], timeLine->StepDates[step+1]);
            double corr = eqFXCorr->getCorrelation();   
            if (ccyTreatment =="S")
            {// struck
                for (i=0; i<=end-start; i++)
                {
                    vol = sqrt(var_dt[i]/dt);
                    if (StepJumpRate != 0)
                    {// Add Jump Vol for quanto adjustment vol
                        jump_vol = termJumpWidth.InterpLinear(timeLine->TradeTime[step]); //Get JumpWidth.
                        jump_volsq = jump_vol*jump_vol;
                        jump_mean = termJumpMean.InterpLinear(timeLine->TradeTime[step]); //Get JumpMean.
                        jump_meansq = jump_mean*jump_mean;
                        vol = sqrt( (var_dt[i]+ StepJumpRate*(jump_meansq+jump_volsq)/(*drift)[i]/(*drift)[i])/dt ); 
                    }
                    var_dt[i] += vol_fx*vol_fx*dt + 2.0*corr*vol*vol_fx*dt;
                }
            }
            if (ccyTreatment=="P") 
            {// quanto
                for (i=0; i<=end -start ; i++)
                {
                    vol = sqrt(var_dt[i]/dt);   //defualt no jump vol are included
                    if (StepJumpRate != 0)
                    {// Add Jump Vol for quanto adjustment vol
                        jump_vol = termJumpWidth.InterpLinear(timeLine->TradeTime[step]); //Get JumpWidth.
                        jump_volsq = jump_vol*jump_vol;
                        jump_mean = termJumpMean.InterpLinear(timeLine->TradeTime[step]); //Get JumpMean.
                        jump_meansq = jump_mean*jump_mean;
                        vol = sqrt( (var_dt[i]+ StepJumpRate*(jump_meansq+jump_volsq)/(*drift)[i]/(*drift)[i])/dt ); 
                    }
                    (*drift)[i] *= exp(-corr*vol_fx*vol*dt);
                }
            }
        }
 

/*  to do:  There may be a constant Pu, Pd and Pm by setting Tree cleverly

       // calculate prob for the step
        double mu = StepForward[step+1]/StepForward[step]; // step drift
        double m = Stock[1-CurrIdx][0]/Stock[CurrIdx][0]/mu - 1.0;
        double u = Stock[1-CurrIdx][1]/Stock[CurrIdx][0]/mu - 1.0;
        double d = Stock[1-CurrIdx][-1]/Stock[CurrIdx][0]/mu - 1.0;
        Pu = (var_dt[0] + d*m)/(u-d)/(u-m);
        Pd = (var_dt[0] + u*m)/(d-u)/(d-m);
        Pm = 1.0 - Pu - Pd;

        // protect rounding for ~0 vol interval
        if (var_dt[0]<1.0e-8)
        {
            if (Pu < 0.0 && Pu>-1.0e-10) Pu = 0.0;
            if (Pm < 0.0 && Pm>-1.0e-10) Pm = 0.0;
            if (Pd < 0.0 && Pd>-1.0e-10) Pd = 0.0;
        }
        if(Pu<0.0 || Pd<0.0 || Pm<0.0)
            fastRoll = false;     // -ve prob, needs per node recalc

        if (!fastRoll)
            PerNodeProb = true;*/
        
        PerNodeProb = true;
        bool fastRoll= false;
        return fastRoll;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** calculate vol for current step for a set of spot levels
    returns number of vol calculated - one (flat for all node) or num */
int CTree1fCEVJ::GetStepVol(int step, vector<double>& vol, const double* s, int bot, int top)
{//to do:  To use BSSmoothing, vol should be consider.
 
    vol.resize(top-bot+1);
    if (step >= timeLine->NumOfStep)
        step =timeLine->NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    
    vector<double> diff_drift;
    diff_drift.resize(top-bot+1);
    int i;
//    double jump_vol = 0.0;
    double jump_vol, jump_mean, jump_volsq, jump_meansq;

    double drift_fwd = StepForward[step+1]/StepForward[step];    
    for (i=bot; i<top+1; i++)
        diff_drift[i-bot] = drift_fwd;
    double step_jump_rate = VolCEVJ->CalcStepJumpRate(timeLine->StepDates[step], timeLine->StepDates[step+1]);

/*    if (DEBUG_BSSmooth ==0)
    {
        bool tmpStore = IsStoreJumpNode;
        IsStoreJumpNode = false;
        CalcStepDriftAndVar(s, bot, top, vol, &diff_drift);  //Get local variance.
        IsStoreJumpNode = tmpStore;
        jump_vol = termJumpWidth.InterpLinear(timeLine->TradeTime[step]); //Get JumpWidth.
        for (i=bot; i<top+1; i++)
            vol[i-bot] = sqrt( vol[i-bot]/timeLine->TradeYrFrac[step+1] + step_jump_rate*jump_vol*jump_vol/2.0 ); // not exactly right for jumpMean != 0
            //vol[i-bot] = sqrt( vol[i-bot] + step_jump_rate*jump_vol*jump_vol/2.0 ); // not exactly right for jumpMean != 0
    }
    else*/ if (DEBUG_BSSmooth == 1)
    {
        double cevpower = termCEVPower.InterpLinear(timeLine->TradeTime[step]);
        double dt = timeLine->TradeYrFrac[step+1];
        double var_diff = Variance[step+1] - Variance[step];
    
        if (dt == 0.0)
        {// avoid dividing zero.  It could be happen, T=date1(EOD) T+1 = date1+1day (SOD)
            for (i=bot; i<=top; i++)
                vol[i-bot] = 0.0;
        }
        else
        {
            VolCEVJ->CalcLocVar(s, diff_drift, cevpower, 0.0, dt, var_diff ,bot, top, vol);

            jump_vol = termJumpWidth.InterpLinear(timeLine->TradeTime[step]); //Get JumpWidth.
            jump_volsq = jump_vol*jump_vol;
            jump_mean = termJumpMean.InterpLinear(timeLine->TradeTime[step]); //Get JumpMean.
            jump_meansq = jump_mean*jump_mean;
            for (i=bot; i<top+1; i++)
                vol[i-bot] = sqrt( (vol[i-bot]+ step_jump_rate*(jump_meansq+jump_volsq)/diff_drift[i-bot]/diff_drift[i-bot])/timeLine->TradeYrFrac[step+1] ); // not exactly right for jumpMean != 0
        }
    }
    else
    {// just to turn off penultimate BS smoothing !!!
        for (i=bot; i<=top; i++)
            vol[i-bot] = 0.0;
    }
    return 1;
}

/** calculate node spacing for each time segment, usually just one segment */
void CTree1fCEVJ::CalcNodeSpacing()
{// to do:  It may be possible to set Node as Constant Pu, Pm and Pd.
    static const string method = "CTree1fCEVJ::CalcNodeSpacing";
    try {
        double v_dt = 0.0;
        double mean, spacing, lambda, assetSpot, m;
        int startStep = 0;
        int startStepNoCrit = 0;
        int i, j, initStep;

        int numOfSeg = timeLine->SegmentEnd.size();
        NodeSpace.resize(numOfSeg);
        for (i=0; i<numOfSeg; i++)
        {
            initStep = timeLine->SegmentEndNoCrit[i] - startStepNoCrit;
            v_dt = (Variance[timeLine->SegmentEnd[i]] - Variance[startStep])/initStep;
            //From Akasaksa's TreeCEVJ.cpp TreeSweep line 320.
            spacing = sqrt(3.0*v_dt); // node spacing

            NodeSpace[i] = spacing;
            startStep = timeLine->SegmentEnd[i];
            startStepNoCrit = timeLine->SegmentEndNoCrit[i];
        }

        //From Akasaksa's TreeCEVJ.cpp TreeSweep line 320.
        assetSpot = StepForward[0];

        // last v_dt is used if more than one segments ???
        spacing = sqrt(3.0*v_dt); // node spacing

        // tree centre node stored for same tree tweak
        for (j = timeLine->NumOfStep; j >=0; j--)
        {
            lambda = termCEVPower.InterpLinear(timeLine->TradeTime[j]);
            //if (CalcMethod == CEV_TREE || CalcMethod == CEV_MC)
            mean = ::pow(assetSpot, 1-lambda)/(1-lambda);
            /* to do: Jump Diffusion
            else
                mean = lambda*log(AssetSpot) + (1-lambda)*AssetSpot/VolTerm.SpotRef;*/
            m = StepForward[j];
            CentreNode[j] = CalcDrift(mean, Variance[j], log(m/assetSpot), timeLine->TradeTime[j]);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**   tree nodes set up */
void CTree1fCEVJ::NodeSetup(int step, int idx)
{
    int i, j, k;

    // NodeSetUp is assumed to call at once for one time slice.
    IsStoreJumpNode = true;

    //keep the same size
    BotDim[idx] = BotDim[1-idx];
    TopDim[idx] = TopDim[1-idx];
    BotClip[idx] = BotClip[1-idx];
    TopClip[idx] = TopClip[1-idx];
    double drift = 1.0; // used only for rolling boundary nodes at step < last step
    if (step < timeLine->NumOfStep)
        drift = StepForward[step+1]/StepForward[step]; // step drift

    if (RollDirection==-1)
        drift = 1.0/drift;

    NodeSetupSpace(Stock[idx], -BotClip[idx]-1, TopClip[idx]+1, 
                   Stock[1-idx], -BotClip[1-idx]-1, TopClip[1-idx]+1, 
                   CentreNode[step], NodeSpace[GetTreeSeg(step)], step, drift);
    
    // trim tree but keep nodes for jump interpolations
    //Because we don't know bot/top level before NodeSetupSpacing, we need those here.
    // node is not trimmed if it can be useful for jump interpolation.

    int target_bot=BotClip[idx], target_top=TopClip[idx];
    double ceiling3 = Stock[CurrIdx][TopClip[idx]+3];
    double ceiling2 = Stock[CurrIdx][TopClip[idx]+2];
    double floor3 = Stock[CurrIdx][-BotClip[idx]-3];
    double floor2 = Stock[CurrIdx][-BotClip[idx]-2];
    // make sure lowest tree boudary >0
    for (i=-BotClip[idx]-1; i<0; i++)       // BotClip shold be below of CentreNode.
    {
        if (Stock[idx][i] < FP_MIN)
            target_bot --; // reduce lower clip
        else
            break;
    }
    if (BotClip[idx] != target_bot)
    {//avoid 0 node at bot-1 and bot-2.
        BotClip[idx] = target_bot+1;  //Set smallest node as tree edg node -1.
        if (CurrStep == timeLine->NumOfStep)
           Stock[idx][-BotClip[idx]-2] = 0.01*Stock[idx][-BotClip[idx]-1]; // 1% of tree edg node
        else
        {
           Stock[idx][-BotClip[idx]-2] = floor2;
           Stock[idx][-BotClip[idx]-3] = floor3;
        }
    }

    // check jumps before trimming down tree dim
    j = Maths::min(step+1, Ndim); // standard tree boundary for LN
    double s[2];
    vector<double> jump_drift(2);
    vector<double> jumpUp_node(2);
    vector<double> jumpDown_node(2);

    // set it to new values, if changed
    target_bot = BotClip[idx];
    target_top = TopClip[idx];

    DateTime thisStepDate = timeLine->StepDates[CurrStep];
    k = (j<target_bot ? j : target_bot);
    s[0] = Stock[idx][-k];
    k = (j<target_top ? j : target_top);
    s[1] = Stock[idx][k];
    VolCEVJ->CalcJumpDrift(s, thisStepDate, jump_drift, jumpDown_node, jumpUp_node, 0, 1);

    if (jumpDown_node[0] > Stock[idx][-target_bot] && j < target_bot)
        target_bot --;
    if (jumpUp_node[1] < Stock[idx][target_top] && j < target_top)
        target_top --;
    if (RollDirection == -1)
    {  //Check all tree node is within floor to ceiling
        i = -target_bot-1;
        while (floor2 > Stock[CurrIdx][i])
        {
            target_bot --;
            i ++;
        }
        j = target_top+1;
        while (ceiling2 < Stock[CurrIdx][j])
        {
            target_top --;
            i ++;
        }
        //keep Floor and Ceiling Underlying Level.
        if (target_top != TopClip[idx])
        {
            Stock[CurrIdx][target_top+3] = ceiling3;
            Stock[CurrIdx][target_top+2] = ceiling2;
        }
        if (target_bot != BotClip[idx])
        {
            Stock[CurrIdx][-target_bot-3] = floor3;
            Stock[CurrIdx][-target_bot-2] = floor2;
        }
        BotDim[idx] = BotClip[idx] = target_bot;
        TopDim[idx] = TopClip[idx] = target_top;
    }
    else
    {// to do: need to check.
        BotDim[1-idx] = BotClip[1-idx] = target_bot;
        TopDim[1-idx] = TopClip[1-idx] = target_top;
    }
}

/**   tree node space set up */
void CTree1fCEVJ::NodeSetupSpace(double* s, int bot, int top, 
                          double* sLast, int botLast, int topLast, 
                          double f0, double spacing, int step, double driftNext)
{
    static const string method = "CTree1fCEVJ::NodeSetupSpace";
    try {
        // reset it to make sure default only calc one prob per step
        // PerNodeProb = false;

        // spot tree is created with proportional spacing
        double s_start = 0.0;
        for (int j = bot; j<=top; j++)
            s[j] = CalcNodeLevel(f0+spacing*j, s_start, timeLine->TradeTime[step]);

        // set floors and ceilings
        s[bot-2] = 0.0; // lowest floor set to zero
        if (step == timeLine->NumOfStep || (step == 0 && RollDirection==1))
        {
            int  n = Maths::min(step+1, Ndim);
            s[bot-1] = 0.01*s[-n-1]; // use 1% for bottom line, note that LN tree uses 10%
            s[top+2] = 2.0 * s[top];
             s[top+1] = s[n+1] +0.9*(s[n+3]-s[n+1]);
            if (step == 0)
            {
                double drift = StepForward[timeLine->NumOfStep]/StepForward[0];
                s[bot-1] /= drift;
                s[top+2] /= drift;
                s[top+1] /= drift;
            }
        }
        else
        {
            s[bot-1] = driftNext*sLast[botLast-1];
            s[top+2] = driftNext*sLast[topLast+2];
            s[top+1] = driftNext*sLast[topLast+1];
        }
    ASSERT(s[top+1] >= s[top]);
    ASSERT(s[top+2] >= s[top+1]);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/*------------------------------------------------------------------------------
*   Name         :  CTree1fCEV::CalcNodeLevel
*
*   Description  :    calculate node level
*
*   Returns      :  node level
*------------------------------------------------------------------------------*/
double CTree1fCEVJ::CalcNodeLevel(double x, double s_last, double time)
{
//    const double tol = 1e-5;
//    int i, n_iter;
    double result = s_last;
    double lambda;

    lambda = termCEVPower.InterpLinear(time);
    
/*    if (CalcMethod == CEV_TREE || CalcMethod == CEV_MC)
    {*/
        if (x < FP_MIN)
            result = 0.0;
        else
            result = ::pow((1-lambda)*x, 1/(1-lambda));
/*    }
    else if (CalcMethod == LAMBDA_TREE || CalcMethod == LAMBDA_MC)
    {// Lambda diffusion
        if (s_last == 0)
        {
            result = FP_MIN*AssetSpot; // initial estimate
            n_iter = 10; // 10 iterations if no initial estimate
        }
        else
            n_iter = 3; // 3 iter for the rest
        for (i=0; i<n_iter; i++)
        {
            x_below = lambda*log(result) + (1.0-lambda)/VolTerm.SpotRef*result;
            if (x_below >= x && result <= FP_MIN*AssetSpot + FP_MIN) // stop if too low
            {
                result = 0.0;
                break;
            }
            else
            {
                temp = 1.0/(lambda + (1.0-lambda)/VolTerm.SpotRef*result);
                temp = (x-x_below)*result*temp*(1.0 +0.5*(x-x_below)*lambda*temp*temp*(1.0
                        +(x-x_below)*temp/3.0*(lambda-2.0*(1.0-lambda)/VolTerm.SpotRef*result)));
                result += temp;
                if (result <= FP_MIN*AssetSpot + FP_MIN)
                {
                    result = 0.0;
                    break;
                }
                if (fabs(temp/result) < tol)
                    break;
            }
        }
    }
    else
        SaveMessage(ModelName, -1, "unknown calc method.");
*/
    return result;
}

/*------------------------------------------------------------------------------
*   Name         :  CTreeCEV::CalcDrift
*
*   Description  :    calculate drift from t=0 to t for setting up the tree
*
*   Returns      :  drift from 0 to t
*------------------------------------------------------------------------------*/
double CTree1fCEVJ::CalcDrift(double x, double v_dt, double asset_drift, double t)
{
    double dx, x0;
//    double temp1, temp2;
    double lambda;  //temp lambda for LambdaTS

    lambda = termCEVPower.InterpLinear(t);
//    if (CalcMethod == CEV_TREE || CalcMethod == CEV_MC)
//    {
        if (x > 0.0)
        {
            x0 = x;
            // 1st pass estimate
            dx = asset_drift*(1-lambda)*x - lambda*v_dt/2/x/(1-lambda);
            x = x0 + dx/2.0; // average value
            // 2nd pass estimate change
            if (x > 0.0)
            {
                dx = asset_drift*(1-lambda)*x - lambda*v_dt/2/x/(1-lambda);
                x = Maths::max(0.0, x0 + dx); // estimate drift
            }
            else
                x = 0.0;;
        }
        else
            x = 0.0;
/*    }
    else if (CalcMethod == LAMBDA_TREE || CalcMethod == LAMBDA_MC)
    {// lambda diffusion
        x0 = x;
        // 1st pass estimate
        temp2 = lambda+(1-lambda)/VolTerm.SpotRef*AssetSpot;
        //temp2 = VolTerm.Lambda+(1-VolTerm.Lambda)/VolTerm.SpotRef*AssetSpot;
        dx = asset_drift*(lambda+(1-lambda)/VolTerm.SpotRef*AssetSpot) - lambda*v_dt/2/temp2/temp2;
        //dx = asset_drift*(VolTerm.Lambda+(1-VolTerm.Lambda)/VolTerm.SpotRef*AssetSpot) - VolTerm.Lambda*v_dt/2/temp2/temp2;
        // 2nd pass estimate change
        x = x0 + dx/2.0; // average value
//        temp1 = CalcNodeLevel(x, AssetSpot);
        temp1 = CalcNodeLevel(x, AssetSpot, t);
        temp2 = lambda+(1-lambda)/VolTerm.SpotRef*temp1;
        //temp2 = VolTerm.Lambda+(1-VolTerm.Lambda)/VolTerm.SpotRef*temp1;
        dx = asset_drift*(lambda+(1-lambda)/VolTerm.SpotRef*temp1) - lambda*v_dt/2/temp2/temp2;
        //dx = asset_drift*(VolTerm.Lambda+(1-VolTerm.Lambda)/VolTerm.SpotRef*temp1) - VolTerm.Lambda*v_dt/2/temp2/temp2;
        x = x0 + dx; // estimate drift
    }
*/
    return x;
}

/*------------------------------------------------------------------------------
*   forward or backward roll the tree at one time step 
*------------------------------------------------------------------------------*/
void CTree1fCEVJ::RollTree(vector<double> *, vector<double> *)
{
    static const string method = "CTree1fCEVJ::Roll";
    try 
    {
        vector<double> drift; 
        vector<double> driftIns; 
        vector<double> insJumpUpNode, insJumpDownNode;
        vector<double> jump_drift(NumOfInsertNode);
        vector<vector<double> > insjump_drifts(NumOfPrice);

        if (DEBUG_JumpFrom == 1)
        {
            // Calculate Jump First at t+1
            StepJumpRate = VolCEVJ->CalcStepJumpRate(timeLine->StepDates[CurrStep], timeLine->StepDates[CurrStep+1]);
            IsStoreJumpNode = true;
            if (StepJumpRate != 0.0)
            { //If there is no jump, JumpUp/DownNode doesn't exist.    
                vector<double> tmp_jump_drift(TopClip[1-CurrIdx] + BotClip[1-CurrIdx] + 1 +6);
                VolCEVJ->CalcJumpDrift(Stock[1-CurrIdx], timeLine->StepDates[CurrStep+1], jump_drift, JumpDownNode, JumpUpNode, -BotClip[1-CurrIdx]-3, TopClip[1-CurrIdx]+3);
                CalcJumpPrice(JumpDownNode, JumpUpNode, -BotClip[1-CurrIdx]-3, TopClip[1-CurrIdx]+3);
                // Calculate Jump Inserted Price (re-calculate drift, jump_up/down node.)
                if (NumOfInsertNode > 0)
                {
                    for (int iPrice=0; iPrice<NumOfPrice; iPrice++)
                    {
                        insjump_drifts[iPrice] = vector<double>(NumOfPrice);
                        latticeProd->moveInsertNode(CurrStep, iPrice);

                        VolCEVJ->CalcJumpDrift(InsStock[1-CurrIdx], timeLine->StepDates[CurrStep+1], 
                                               insjump_drifts[iPrice], insJumpDownNode, insJumpUpNode, 0, NumOfInsertNode-1);
                        CalcJumpPrice(insJumpDownNode, insJumpUpNode,0, NumOfInsertNode-1, true, iPrice);
                    }
                }
            }

            // Combining Jump Price.
//            int numDrift = drift.size();
            int j;
//            double s1 =0.0, mu = 0.0; 
//            double step_rate = StepJumpRate;
            vector<vector<double> > storePrice(NumOfPrice);
            if (RollDirection == -1 && StepJumpRate!=0.0)
            {
                for (int iPrice=0; iPrice<NumOfPrice; iPrice++)
                {
                    storePrice[iPrice] = vector<double>(4);
                    int k=0;
                    for (j=-BotClip[1-CurrIdx]-3; j<=TopClip[1-CurrIdx]+3; j++)
                    {
                        if (j<-BotClip[1-CurrIdx]-1 || j>TopClip[1-CurrIdx]+1)
                            storePrice[iPrice][k++] = NodePrice[1-CurrIdx][iPrice][j];
                        NodePrice[1-CurrIdx][iPrice][j] = NodePrice[1-CurrIdx][iPrice][j]*(1 - StepJumpRate) 
                                                      + StepJumpRate * JumpPrice[1-CurrIdx][iPrice][j];
                      //  Stock[1-CurrIdx][j] *= 1.0+jump_drift[j + BotClip[1-CurrIdx]];   //Shift Node
                    }
                    for (j=0; j<NumOfInsertNode; j++)
                    {
                        InsNodePrice[1-CurrIdx][iPrice][j] = InsNodePrice[1-CurrIdx][iPrice][j]*(1 - StepJumpRate) 
                                                      + StepJumpRate * InsJumpPrice[1-CurrIdx][iPrice][j];
                       // InsStock[1-CurrIdx][j] *= 1.0+insjump_drifts[iPrice][j];   //Shift Ins Node
                    }
                }
            }
            else
            {
                //Nothing need to do.  Jump Prob are already added to NodePrice in CalcJumpPrice.
            }

            int store_CurrIdx = CurrIdx;
            CurrIdx = 1-CurrIdx;
            InterpTreeBoundary();
            CurrIdx = store_CurrIdx;

            CTree1f::RollTree(&drift, &driftIns);

            //Overwrite the boundary prices
            if (StepJumpRate != 0.0){
                double df = (RollDirection == -1? 
                             discYC->pv(timeLine->StepDates[CurrStep], timeLine->StepDates[CurrStep+1]) 
                             : 1.0);
                for (int iPrice=0; iPrice<NumOfPrice; iPrice++){
                    NodePrice[CurrIdx][iPrice][-BotClip[CurrIdx]-3] = df*storePrice[iPrice][0];
                    NodePrice[CurrIdx][iPrice][-BotClip[CurrIdx]-2] = df*storePrice[iPrice][1];
                    NodePrice[CurrIdx][iPrice][ TopClip[CurrIdx]+2] = df*storePrice[iPrice][2];
                    NodePrice[CurrIdx][iPrice][ TopClip[CurrIdx]+3] = df*storePrice[iPrice][3];
                }
                InterpTreeBoundary();
            }

        }
        else
        {
            // Calculate Diffusion Price and InsPrice
            IsStoreJumpNode = true;             // Store JumpUpNode and JumpDownNode -> Once stored, 
            CTree1f::RollTree(&drift, &driftIns);

            // Calculate Jump Price.  
            if (StepJumpRate != 0.0)
            { //If there is no jump, JumpUp/DownNode doesn't exist.    
                CalcJumpPrice(JumpDownNode, JumpUpNode, -BotClip[CurrIdx], TopClip[CurrIdx]);
                // Calculate Jump Inserted Price (re-calculate drift, jump_up/down node.)
                if (NumOfInsertNode > 0)
                {
                    for (int iPrice=0; iPrice<NumOfPrice; iPrice++)
                    {
                        latticeProd->moveInsertNode(CurrStep, iPrice);
                        
                        VolCEVJ->CalcJumpDrift(InsStock[CurrIdx], timeLine->StepDates[CurrStep], 
                                           jump_drift, insJumpDownNode, insJumpUpNode, 0, NumOfInsertNode-1);
                        CalcJumpPrice(insJumpDownNode, insJumpUpNode,0, NumOfInsertNode-1, true, iPrice);
                    }
                }
            }

            // Combining Jump Price.
            int numDrift = drift.size();
            int j;
            double s1 =0.0, mu = 0.0; 
            double step_rate = StepJumpRate;
            if (RollDirection == -1 && StepJumpRate!=0.0)
            {
                for (int iPrice=0; iPrice<NumOfPrice; iPrice++)
                {
                    for (j=-BotClip[CurrIdx]; j<=TopClip[CurrIdx]; j++)
                    {
                        if (numDrift > 1)
                            mu = drift[j+BotClip[CurrIdx]];
                        s1 = mu * Stock[CurrIdx][j];
                        if (s1 < FP_MIN)                  
                        {
                            s1 = 0.0;
                            double step_drift = StepForward[CurrStep+1]/StepForward[CurrStep];
                            double jump_drift = ( (step_drift - 1.0) - (mu -1.0)*(1.0 - StepJumpRate) )
                                                / StepJumpRate;         //calc back Jump_Drift
                            step_rate = Maths::min(1.0, Maths::max(0.0, step_drift/(1.0 + jump_drift)));
                        }
                        else
                            step_rate = StepJumpRate;

                        NodePrice[CurrIdx][iPrice][j] = NodePrice[CurrIdx][iPrice][j]*(1 - step_rate) 
                                                      + step_rate * JumpPrice[CurrIdx][iPrice][j];
                    }
                    for (j=0; j<NumOfInsertNode; j++)
                    {
                        mu = driftIns[j];
                        s1 = mu * InsStock[CurrIdx][j];
                        if (s1 < FP_MIN)
                        {
                            s1 = 0.0;
                            double step_drift = StepForward[CurrStep+1]/StepForward[CurrStep];
                            double jump_drift = ( (step_drift - 1.0) - (mu -1.0)*(1.0 - StepJumpRate) )
                                                / StepJumpRate;         //calc back Jump_Drift
                            step_rate = Maths::min(1.0, Maths::max(0.0, step_drift/(1.0 + jump_drift)));
                        }
                        else
                            step_rate = StepJumpRate;

                        InsNodePrice[CurrIdx][iPrice][j] = InsNodePrice[CurrIdx][iPrice][j]*(1 - step_rate) 
                                                      + step_rate * InsJumpPrice[CurrIdx][iPrice][j];
                    }
                }
            }
            else
            {
                //Nothing need to do.  Jump Prob are already added to NodePrice in CalcJumpPrice.
            }

            InterpTreeBoundary();
        }
#ifdef DEBUG_CEVJ_TREE
/****************************************/
//debug
    char buffer[2048];
    if (DumpFile == NULL && Dumpfilenum==0)
        DumpFile = fopen("C:\\Temp\\debugEDR.dat", "w");
    if (DumpFile != NULL && Dumpfilenum ==0)
    {
        if (CurrStep == timeLine->NumOfStep-1)
        {
            sprintf(buffer, "#step = %d, TradeTime, stock , Price \n", CurrStep+1);
            fprintf(DumpFile, buffer);
            for (int j_debug=-BotClip[CurrIdx]-3; j_debug<=TopClip[CurrIdx]+3; j_debug ++)
            {
                sprintf (buffer, "%f, %f, %f, # %d \n", timeLine->TradeTime[CurrStep+1], Stock[1-CurrIdx][j_debug], NodePrice[1-CurrIdx][iOut][j_debug], j_debug);
                fprintf(DumpFile, buffer);
            }
            sprintf (buffer, "\n");
            fprintf (DumpFile, buffer);
        }
        sprintf(buffer, "#step = %d, TradeTime, stock , Price \n", CurrStep);
        fprintf(DumpFile, buffer);
        for (int j_debug=-BotClip[CurrIdx]-3; j_debug<=TopClip[CurrIdx]+3; j_debug ++)
        {
            sprintf (buffer, "%f, %f, %f # %d \n", timeLine->TradeTime[CurrStep], Stock[CurrIdx][j_debug], NodePrice[CurrIdx][iOut][j_debug], j_debug);
            fprintf(DumpFile, buffer);
        }
        sprintf (buffer, "\n");
        fprintf (DumpFile, buffer);
    }
    if (CurrStep == 0)
    {
        fclose(DumpFile);
        Dumpfilenum=1;
    }
/***************************************/
#endif

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//------------------------------
//Calculate Jump Price
//Note:  Bot and Top should be from Tree ceil and fllor.  
//       Boundary area is forbidden due to fail TreeInterp.
//
//      isInsertNode and pPrice are optional for InsertNode case.
//      InsertNode need to use MoveInsNode outside this function, so
//      pPrice are used for that case.
//------------------------------
void CTree1fCEVJ::CalcJumpPrice(const vector<double>& d_node, 
                                const vector<double>& u_node,
                                const int bot,const int top, 
                                bool isInsertNode, int pPrice)
{   
    int bot_guess, top_guess;
    int iPrice;
    int i;

    int pStart, pEnd;
    if (isInsertNode)
    {
        pStart = pPrice;
        pEnd = pPrice+1;
    }
    else
    {
        pStart = 0;
        pEnd = NumOfPrice;
    }

    double df = (RollDirection == -1? 
                 discYC->pv(timeLine->StepDates[CurrStep], timeLine->StepDates[CurrStep+1]) : 1.0);
    bot_guess = int (bot/2);
    top_guess = int (top/2);

    if (RollDirection == -1) 
    {// price
        double d_price, u_price;
        if (DEBUG_JumpFrom == 0)
        {
            for (iPrice=pStart; iPrice<pEnd; iPrice++)
            {
                for (i=bot; i<=top; i++)
                {
//                    if (Stock[CurrIdx][i] < FP_MIN)
//                        JumpPrice[CurrIdx][iPrice][i] = 0.0;
//                    else
 //                   {
                        d_price = TreeInterp(d_node[i-bot], true, iPrice, i, &bot_guess);
                        u_price = TreeInterp(u_node[i-bot], true, iPrice, i, &top_guess);
                        if (isInsertNode)
                            InsJumpPrice[CurrIdx][iPrice][i] = (d_price+u_price)/2.0;
                        else
                            JumpPrice[CurrIdx][iPrice][i] = (d_price+u_price)/2.0;
//                    }
                }
            }
        }
        else if(DEBUG_JumpFrom == 1)
        {
            for (iPrice=pStart; iPrice<pEnd; iPrice++)
            {
                for (i=bot; i<=top; i++)
                {
//                    if (Stock[1-CurrIdx][i] < FP_MIN)
//                        JumpPrice[1-CurrIdx][iPrice][i] = 0.0;
//                    else
//                    {
                        d_price = TreeInterp(d_node[i-bot], false, iPrice, i, &bot_guess);
                        u_price = TreeInterp(u_node[i-bot], false, iPrice, i, &top_guess);
                        if (isInsertNode)
                            InsJumpPrice[1-CurrIdx][iPrice][i] = (d_price+u_price)/2.0;
                        else
                            JumpPrice[1-CurrIdx][iPrice][i] = (d_price+u_price)/2.0;
//                    }
                }
            }

        }
        else
        {
            for (iPrice=pStart; iPrice<pEnd; iPrice++)
            {
                for (i=bot; i<=top; i++)
                {
//                    if (Stock[CurrIdx][i] < FP_MIN)
//                        JumpPrice[CurrIdx][iPrice][i] = 0.0;
//                    else
//                    {
                        d_price = TreeInterp(d_node[i-bot], false, iPrice, i, &bot_guess);
                        u_price = TreeInterp(u_node[i-bot], false, iPrice, i, &top_guess);
                        if (isInsertNode)
                            InsJumpPrice[CurrIdx][iPrice][i] = df*(d_price+u_price)/2.0;
                        else
                            JumpPrice[CurrIdx][iPrice][i] = df*(d_price+u_price)/2.0;
//                    }
                }
            }
        }
    }
    else
    {// fwd prob calc       : Need to check
        double d_prob, u_prob;
        int j_d, j_u;
        for (i=bot; i<=top; i++)
        {
            for (iPrice=pStart; iPrice<pEnd; iPrice++)
            {
                d_prob = TreeInterp(d_node[i-bot], false, iPrice, bot_guess, &j_d);
                NodePrice[1-CurrIdx][iPrice][j_d] += 0.5*StepJumpRate*d_prob*NodePrice[CurrIdx][iPrice][i];
                NodePrice[1-CurrIdx][iPrice][j_d+1] += 0.5*StepJumpRate*(1.0-d_prob)*NodePrice[CurrIdx][iPrice][i];
                u_prob = TreeInterp(u_node[i-bot], false, iPrice, top_guess, &j_u);
                NodePrice[1-CurrIdx][iPrice][j_u] += 0.5*StepJumpRate*u_prob*NodePrice[CurrIdx][iPrice][i];
                NodePrice[1-CurrIdx][iPrice][j_u+1] += 0.5*StepJumpRate*(1.0-u_prob)*NodePrice[CurrIdx][iPrice][i];

                /* store Probality between two node for Asian  Keiji
                Prob_Path[iPrice][i][j_d]    += 0.5*StepJumpRate*d_prob*NodePrice[CurrIdx][iPrice][i];
                Prob_Path[iPrice][i][j_d+1] += 0.5*StepJumpRate*(1.0-d_prob)*NodePrice[CurrIdx][iPrice][i];
                Prob_Path[iPrice][i][j_u]    += 0.5*StepJumpRate*u_prob*NodePrice[CurrIdx][iPrice][i];
                Prob_Path[iPrice][i][j_u+1]    += 0.5*StepJumpRate*(1.0-u_prob)*NodePrice[CurrIdx][iPrice][i];*/
            }
        }
    }

}

/** decide if a new tree should be built or use the same tree */
bool CTree1fCEVJ::isTreeRebuilt(const CControl* control)
{
    SensitivitySP sens = control->getCurrentSensitivity();

// need to reconstruct center nodes for these sensitivities (specific to this model)
    return CTree1f::isTreeRebuilt( control )
        || CEVPowerParallel::TYPE->isInstance( sens.get() )
        || CEVPowerPointwise::TYPE->isInstance( sens.get() );
}

class Tree1fCEVJHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CTree1fCEVJ, clazz);
        SUPERCLASS(CTree1f);

        FIELD(DEBUG_BSSmooth, "0=On, else=Off");
        FIELD_MAKE_OPTIONAL(DEBUG_BSSmooth);
        FIELD(DEBUG_JumpFrom, "0=Jump@t+1, 1=Jump@t+1, else=Off(default)");
        FIELD_MAKE_OPTIONAL(DEBUG_JumpFrom); 

        EMPTY_SHELL_METHOD(defaultTree1fCEVJ);
    }

    static IObject* defaultTree1fCEVJ(){
        return new CTree1fCEVJ();
    }

    // for IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(CTree1fCEVJ::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }
};

CClassConstSP const CTree1fCEVJ::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("CTree1fCEVJ::IIntoProduct",
                                    typeid(CTree1fCEVJ::IIntoProduct),
                                    Tree1fCEVJHelper::loadIntoProduct);


CClassConstSP const CTree1fCEVJ::TYPE = CClass::registerClassLoadMethod(
    "Tree1fCEVJ", typeid(CTree1fCEVJ), Tree1fCEVJHelper::load);

DRLIB_END_NAMESPACE
