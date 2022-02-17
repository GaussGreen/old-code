//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fLV.cpp
//
//   Description : one factor trinomial tree for local vol process.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/StruckAsset.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ProtAsset.hpp"

DRLIB_BEGIN_NAMESPACE

const string CTree1fLV::SPEEDUP_METHOD_NONE     = "none";
const string CTree1fLV::SPEEDUP_METHOD_GRID     = "grid";

CTree1fLV::CTree1fLV(CClassConstSP type):CTree1f(type), volType(VolSurface::TYPE->getName()),
                       useTweakingForTimeDerivs(true),
                       tweakStrikeUnscaled(0.01),
                       tweakTimeUnscaled(0.001),
                       probDensRatioMin(0.01),
                       useMidPoint(false),
                       speedup(SPEEDUP_METHOD_NONE){}

/** get processed vol*/
// need to examine if this is general
void CTree1fLV::InitVol()
{
    static const string method("CTree1fLV::InitVol");
    try {
        DateTime startDate = prod->getStartDate();
        bool isFwdStart = startDate > getValueDate();

        // *** to do: just an example here, needs to review in actual use
        LocVolRequestSP volRequest(new LocVolRequest(startDate,
                                                    isFwdStart,
                                                    false,
                                                    useTweakingForTimeDerivs,
                                                    useMidPoint,
                                                    tweakStrikeUnscaled,
                                                    tweakTimeUnscaled,
                                                    probDensRatioMin,
                                                    speedup));
                               
        
        CAssetConstSP plainAsset = Underlier;

        if (StruckEquity::TYPE->isInstance(Underlier) && !(AssetUtil::isBasket(Underlier)))
        {
            plainAsset = StruckEquityConstSP::dynamicCast(Underlier)->getPlainAsset();
            fxAsset = StruckEquityConstSP::dynamicCast(Underlier)->getFX();
            eqFXCorr = StruckEquityConstSP::dynamicCast(Underlier)->getCorrelation();

            if (!fxAsset){
                throw ModelException(method, "NULL fx asset");
            }
            if (!eqFXCorr){
                throw ModelException(method, "NULL correlation");
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
        else if (StruckAsset::TYPE->isInstance(Underlier) || ProtAsset::TYPE->isInstance(Underlier)) {
            throw ModelException(method, "1-factor LV tree does not support struck or protected non-equity assets");
        }

        VolLV = CVolProcessedDVFSP::dynamicCast(
            CVolProcessedSP(plainAsset->getProcessedVol(volRequest.get())));

        // get time metric
        timeMetric = VolLV->GetTimeMetric();
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calculate a term structure vol^2 or just one point at maturity */
void CTree1fLV::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
                           const DateTime& matDate, CTermStructure& v2)
{
    static const string method = "CTree1fLV::CalcV2Term";
    try {
        int i;
        DateTimeArray benchMark;
        DoubleArray vol_sq;
        // a simple but very poor way of creating vol term structure
        // access to vol term structure denied

        // *** this is not exact for LV where param function is used to give smooth T,K values
        // but it is good enough - LV cannot have true equal variance step anyway.

        // write down the bench marks
        const int benchMarkSize = 19;
        string benchMarkString[benchMarkSize] = {"1D", "1W", "2W", "3W", "1M", "2M", "3M", "6M", "9M", 
                                                 "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "10Y", "20Y", "30Y"};
        // create bench mark dates
        DateTime valDatePM(valDate.getDate(), DateTime::END_OF_DAY_TIME);
        benchMark.resize(benchMarkSize*3-2);
        benchMark[0] = MaturityPeriod(benchMarkString[0]).toDate(valDatePM);
        for (i=1; i<benchMarkSize; i++)
        {
            benchMark[3*i-1] = MaturityPeriod(benchMarkString[i]).toDate(valDatePM);
            benchMark[3*i-2] = benchMark[3*i-1].rollDate(-1); // we won't be wrong by 1 day !
            benchMark[3*i] = benchMark[3*i-1].rollDate(1);
        }
        // use dates from start date onwards
        DateTimeArray benchMarkTrunc;
        for (i=0; i<benchMarkSize*3-2; i++)
        {
            if (benchMark[i] > startDate)
            {
                benchMarkTrunc.push_back(benchMark[i]);
            }
        }

        // retrieve variance for each date from the processed vol object
        vector<double> yrs;
        vol_sq.resize(benchMarkTrunc.size());
        yrs.resize(benchMarkTrunc.size());


        // *** this is only for setting up time line, using ATM vol should be ok ?
        double strike;
        strike = Underlier->fwdValue(startDate);

        for (i=0; i<vol_sq.size(); i++)
        {
            vol_sq[i] = VolLV->computeImpVol(benchMarkTrunc[i], strike);
            if (vol_sq[i] <= 0.0)
                throw ModelException(method, "received vol<=0 : cannot create tree with such vol");
            vol_sq[i] *= vol_sq[i];
            yrs[i] = VolLV->calcTradingTime(startDate, benchMarkTrunc[i]);
        }

        // populate it to a vol^2 term structure
        v2.Populate(startDate, vol_sq.size(), benchMarkTrunc.begin(),
                    yrs.begin(), vol_sq.begin());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
// set up variance array 
void CTree1fLV::PostSetup()
{
    Variance.resize(timeLine->NumOfStep+1);
    // *** requiring ATM vol should be ok, Variance is used for setting up nodes only
    double strike;

    DateTime startDate = timeLine->StepDates[0];

    string ccyTreatment = prod->getCcyTreatment();

    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
        strike = Underlier->fwdValue(startDate)/fxAsset->fwdValue(startDate);
    else
        strike = Underlier->fwdValue(startDate);

     Variance[0] = 0.0; // start with 0
    for (int i=1; i<=timeLine->NumOfStep; i++)
    {
        Variance[i] = VolLV->computeImpVol(timeLine->StepDates[i], strike);
        if (Variance[i] <= 0.0)
            throw ModelException("CTree1fLV::PostSetup", "received vol<=0 : cannot create tree with such vol");
        if (ccyTreatment==CAsset::CCY_TREATMENT_STRUCK)
        {
            double vol_fx = volFXBS->CalcVol(startDate, timeLine->StepDates[i]);
            double corr = eqFXCorr->getCorrelation();               
            Variance[i] = sqrt(Variance[i]*Variance[i] + vol_fx*vol_fx + 2.0*corr*Variance[i]*vol_fx);
        }
        Variance[i] *= Variance[i] * timeLine->TradeTime[i];
    }
}

/** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
    they must be either one per step or one per node (from -BotClip to TopClip)
    returns true if constant probabilities can be used for all nodes (fast roll) */
bool CTree1fLV::CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& var_dt,    vector<double>* drift)
{
    static const string method = "CTree1fLV::CalcStepDriftAndVar";
    try {
        ASSERT(CurrStep < timeLine->NumOfStep);

        double dt = timeLine->TradeYrFrac[CurrStep+1];
        double sqr_dt = sqrt(dt);

        // calc loc vol
        GetStepVol(CurrStep, var_dt, s, start, end);
        for (int j=0; j<(int)var_dt.size(); j++)
        {
            var_dt[j] *= var_dt[j]*dt;
        }
       
        // drift
        if (drift)
        {
            int i;
            // no specific dollar div treatment
            double base_drift = StepForward[CurrStep+1]/StepForward[CurrStep];
            // clear all data once, because resize will not clear the existing numbers!!
            // otherwise InserNode calculation will rememeber the previous number.
            drift->clear(); 
            drift->resize(end-start+1, base_drift);

            string ccyTreatment = prod->getCcyTreatment();

            if (ccyTreatment==CAsset::CCY_TREATMENT_PROTECTED) 
            {// quanto
                double vol_fx = volFXBS->CalcVol(timeLine->StepDates[CurrStep], 
                                                 timeLine->StepDates[CurrStep+1]);
                double vol_fx_sqdt = sqr_dt*vol_fx;
                double vol_eq_sqdt;
                double corr = eqFXCorr->getCorrelation();   
                for (i=0; i<=end -start ; i++)
                {
                    vol_eq_sqdt = sqrt(var_dt[i]);
                    (*drift)[i] *= exp(-corr * vol_fx_sqdt * vol_eq_sqdt);
                }
            }
        }
        // now make v_dt in unit of fwd
        for (int k=0; k<(int)var_dt.size(); k++)
            var_dt[k] = exp(var_dt[k])-1.0;

        PerNodeProb = true;

        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calculate vol for current step for a set of spot levels
    returns number of vol calculated - one (flat for all node) or num */
int CTree1fLV::GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end)
{
    vol.resize(end-start+1);

    if (step >= timeLine->NumOfStep)
        step =timeLine->NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    

    DateTimeArray t(2);
    t[0] = timeLine->StepDates[step];
    // go from now to next timepoint or to tomorrow whichever is further.
    // Local vol is not very accurate if you ask for it over too short an interval.
    if (timeLine->StepDates[step].rollDate(1) > timeLine->StepDates[step+1]) {
        t[1] = timeLine->StepDates[step].rollDate(1);
    } else {
        t[1] = timeLine->StepDates[step+1];
    }

    string ccyTreatment = prod->getCcyTreatment();

    if(ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
    {// struck
        vector<double> stmp;
        stmp.resize(end-start+1);
        vector<double>::iterator s = stmp.begin();
        for (int n_zero=0; n_zero<=end - start; n_zero++)
        {    
            s[n_zero] = s_inp[n_zero + start];
            s[n_zero] = s[n_zero]/fxAsset->fwdValue(t[0]);
            if (s[n_zero] < FP_MIN)
                s[n_zero] = FP_MIN;
        }
        CSliceDouble spots(&(*s), end - start + 1);
        CSliceDouble locVol(&vol[0], end - start + 1);
     
        // uses the version with cahing of forward values
        volCalculator->CalcLocVol(spots, 
                                  step,
                                  locVol);
        
        double vol_fx = volFXBS->CalcVol(t[0], t[1]);
        double corr = eqFXCorr->getCorrelation();   
        for (int i=0; i<=end - start; i++)
            vol[i] = sqrt(vol[i]*vol[i] + vol_fx*vol_fx + 2.0*corr*vol[i]*vol_fx);
    }
    else
    {
        double* s = const_cast<double*>(s_inp);
        for (int n_zero=start; n_zero<=end; n_zero++)
        {               
            if (s[n_zero] < FP_MIN)
                s[n_zero] = FP_MIN;
        }
        CSliceDouble spots(&*(s+start), end - start + 1);
        CSliceDouble locVol(&vol[0], end - start + 1);
        
        // use the version with caching of forward values
        volCalculator->CalcLocVol(spots, 
                                  step,
                                  locVol);
    }
    return (end-start+1);
}

/** calculate node spacing for each time segment, usually just one segment */
void CTree1fLV::CalcNodeSpacing()
{
    static const string method = "CTree1fLV::CalcNodeSpacing";
    try {
        double v_dt;
        int startStep = 0;
        int startStepNoCrit = 0;
        int i, j, initStep;
        double root1, root2, alpha;

        int numOfSeg = timeLine->SegmentEnd.size();
        NodeSpace.resize(numOfSeg);

        for (i=0; i<numOfSeg; i++)
        {
            initStep = timeLine->SegmentEndNoCrit[i] - startStepNoCrit;
            v_dt = (Variance[timeLine->SegmentEnd[i]] - Variance[startStep])/initStep;
            //NodeSpace[i] = sqrt(v_dt/(TreeAlpha/timeLine->PeakVarRatio)); // simple but does not work for high vol

            // search for max var in the segment, really for equal time step use
            // but keep it until equal variance step creation matches CalcVar() exactly
            for (j=startStep; j<timeLine->SegmentEnd[i]; j++)
            {
                if (v_dt < Variance[j+1] - Variance[j])
                    v_dt = Variance[j+1] - Variance[j];
            }

            double m = exp(v_dt/2.); // this is for variance adjusted centre node
            //double m =1.0;
            alpha = (-exp(v_dt)*m*m + 1-2*TreeAlpha)/(TreeAlpha-1+m);

            SolveQuadratic(1.0, alpha, 1.0, &root1, &root2);
            ASSERT(root2>1.0);
            NodeSpace[i] = log(root2);
            startStep = timeLine->SegmentEnd[i];
            startStepNoCrit = timeLine->SegmentEndNoCrit[i];
        } 
        // calculate centre nodes for all steps
        for (i=0; i<=timeLine->NumOfStep; i++)
        {
            // variance adjusted centre nodes give better sampling for prob space 
            // and thus slightly better convergence than fwds centre nodes.
            CentreNode[i] = StepForward[i]*exp(-Variance[i]/2.0);
 
            // but using fwds as centre nodes tend to range around usual boundaries (strikes etc.)
            // and may be more relvevant ?
            //CentreNode[i] = StepForward[i];
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** set up tree nodes at the given step */
void CTree1fLV::NodeSetup(int step, int idx)
{
    // leave 3 more nodes on each side of tree
    if (useGammaNodes)
        BotDim[idx] = TopDim[idx] = Maths::min(step+3, NLimit);
    else
        BotDim[idx] = TopDim[idx] = Maths::min(step+3, Ndim);
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

/**   tree node space set up */
void CTree1fLV::NodeSetupSpace(double* s, int bot, int top, 
                          double* sLast, int botLast, int topLast, 
                          double f0, double spacing, int step, double driftNext)
{
    static const string method = "CTree1fLV::NodeSetupSpace";
    try {
        // reset it to make sure default only calc one prob per step
        PerNodeProb = false;
        int j;

        if (useGammaNodes){
            // overwrite Dim & Clip 
            int idx = CurrIdx;
            // set up information
            if (CurrStep == timeLine->NumOfStep && CacheMode != USE_CACHE)
                gammaNodeMaps[CurrStep].SetAtMat(NLimit, gammaNodesInterval);
            else{
                if (CacheMode != USE_CACHE){
                    gammaNodeMaps[CurrStep].Set(-BotDim[idx]-3, TopDim[idx]+3, gammaNodeMaps[CurrStep+1]);
                }
                bot = gammaNodeMaps[CurrStep].getBotBefore() + 2;
                top = gammaNodeMaps[CurrStep].getTopBefore() - 2;
                ASSERT(bot<0 && top>0);     // at least, three nodes should be remained in both side.
                BotClip[idx] = -bot-1;
                TopClip[idx] = top-1;
            }
            BotDim[idx] = TopDim[idx] = Maths::min(step/gammaNodesInterval+1, Ndim);
            ASSERT(BotClip[idx]>=0 && TopClip[idx] >=0);
            NodeSetWithAddNode(s, bot, top, f0, spacing);
        }else{
            // spot tree is created with proportional spacing
            double spacingExp = exp(spacing);
            // starting node
            s[bot] = f0*exp(spacing*bot);
            for (j = bot+1; j<=top; j++)
                s[j] = s[j-1]*spacingExp;        
        }

        // set floors and ceilings
        s[bot-2] = 0.0; // lowest floor set to zero
        if (step == timeLine->NumOfStep || (step == 0 && RollDirection==1))
        {
            double scaleTop = 2.0;            
            // looking for the biggest scale among the segment
            double topSpacing = spacing;
            for (int k=0;k<(int)timeLine->SegmentEnd.size();k++){
                if (NodeSpace[k]>spacing)
                    topSpacing = NodeSpace[k];
            }
        
            s[bot-1] = 0.1*s[bot];
            if (useGammaNodes)
                s[top+2] = scaleTop*f0*exp(topSpacing*(Ndim*gammaNodesInterval+1)+Variance[timeLine->NumOfStep]/2.0);
            else
                s[top+2] = scaleTop*f0*exp(topSpacing*(Ndim+1)+Variance[timeLine->NumOfStep]/2.0);
            s[top+1] = s[top] +0.9*(s[top+2]-s[top]);
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
            if (useGammaNodes){
                botLast = gammaNodeMaps[CurrStep+1].getBotAfter() + 2;
                topLast = gammaNodeMaps[CurrStep+1].getTopAfter() - 2;
            }
            s[bot-1] = driftNext*sLast[botLast-1];
            s[top+2] = driftNext*sLast[topLast+2];
            s[top+1] = driftNext*sLast[topLast+1];
        }

        // check for the segment setting
        if (step < timeLine->NumOfStep){
            if (s[top+1] >= s[top] && GetTreeSeg(step) != GetTreeSeg(step+1))
                ModelException(method, "NodeSpace is changed to too big due to segment setting.");
        }

        ASSERT(s[top+1] >= s[top]);
        ASSERT(s[top+2] >= s[top+1]);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP CTree1fLV::createMDF() const{

    return MarketDataFetcherSP(new MarketDataFetcherLN(volType,useCcyBasis));
}

IModel::WantsRiskMapping CTree1fLV::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class Tree1fLVHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CTree1fLV, clazz);
        SUPERCLASS(CTree1f);
        EMPTY_SHELL_METHOD(defaultTree1fLV);
        FIELD(volType, "Type of log-normal vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(useTweakingForTimeDerivs, "useTweakingForTimeDerivs");
        FIELD_MAKE_OPTIONAL(useTweakingForTimeDerivs);
        FIELD(tweakStrikeUnscaled, "tweakStrikeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakStrikeUnscaled);
        FIELD(tweakTimeUnscaled, "tweakTimeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakTimeUnscaled);
        FIELD(probDensRatioMin, "probDensRatioMin");
        FIELD_MAKE_OPTIONAL(probDensRatioMin);
        FIELD(useMidPoint, "useMidPoint");
        FIELD_MAKE_OPTIONAL(useMidPoint);    
        FIELD(speedup, "speedup");
        FIELD_MAKE_OPTIONAL(speedup);        
    }

    static IObject* defaultTree1fLV(){
        return new CTree1fLV();
    }
};

// for new FDModel interface
void CTree1fLV::finaliseModel(CControl*    control) {
    CTree1f::finaliseModel(control);

    /* initialize the caching of forward values for the local volatility */
    volCalculator = refCountPtr<CVolProcessedDVF::IVolCalculator>(VolLV->CreateVolCalculator(timeLine->StepDates));
}

CClassConstSP const CTree1fLV::TYPE = CClass::registerClassLoadMethod(
    "Tree1fLV", typeid(CTree1fLV), Tree1fLVHelper::load);

DRLIB_END_NAMESPACE
