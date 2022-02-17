#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MDFAssetVol.hpp"

#include "edginc/ProtEquity.hpp"

#include "edginc/FD1DLV.hpp"
#include "edginc/UtilFuncs.hpp"

DRLIB_BEGIN_NAMESPACE

////////////////////////// FD1DLV /////////////////////////
void FD1DLV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD1DLV, clazz);
    SUPERCLASS(FD1D);
    EMPTY_SHELL_METHOD(defaultFD1DLV);

    FIELD(volType, "Type of log-normal vol to use");
    FIELD_MAKE_OPTIONAL(volType);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD1DLV::TYPE = CClass::registerClassLoadMethod(
    "FD1DLV", typeid(FD1DLV), load);

//----------------------------------------------------------------

FD1DLV::FD1DLV( const CClassConstSP & type ) :
    FD1D( type ),
    volType(VolSurface::TYPE->getName()),
    useTweakingForTimeDerivs(true),
    tweakStrikeUnscaled(0.01),
    tweakTimeUnscaled(0.001),
    probDensRatioMin(0.01),
    useMidPoint(false)
{}

//----------------------------------------------------------------

FD1DLV::~FD1DLV(){
}

//----------------------------------------------------------------

/** Override default createMDF in order to set the right MDF */

MarketDataFetcherSP FD1DLV::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(volType));
}

//----------------------------------------------------------------

IModel::WantsRiskMapping FD1DLV::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

//----------------------------------------------------------------

/** retrieving market data from MDF */
void FD1DLV::retrieveFactor()
{
    static const string method = "FD2DSVCJ::retrieveFactor";
    try{
        FD1D::retrieveFactor();
     
        // init vols    
        DateTime startDate = prod->getStartDate();
        DateTime valueDate = getValueDate();

        bool isFwdStart = startDate > valueDate;

        if(!!prod && !isFwdStart)
            startDate = getValueDate();

        // local vol request for each asset
        LocVolRequest volRequest(startDate,
                                 isFwdStart, // fwd start flag
                                 false,
                                 useTweakingForTimeDerivs,
                                 useMidPoint,
                                 tweakStrikeUnscaled,
                                 tweakTimeUnscaled,
                                 probDensRatioMin);


        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD1DLV::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }


/*      To add
        if (StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying)))
        {
            plainAsset = StruckEquityConstSP::dynamicCast(underlying)->getPlainAsset();
            fxAsset = StruckEquityConstSP::dynamicCast(underlying)->getFX();
            eqFXCorr = StruckEquityConstSP::dynamicCast(underlying)->getCorrelation();

            if (!fxAsset){
                throw ModelException("FD1DLV::InitVol", "NULL fx asset");
            }
            if (!eqFXCorr){
                throw ModelException("FD1DLV::InitVol", "NULL correlation");
            }
            // create an atm vol request for the fx side of things
            CVolRequestSP atmInterp(new ATMVolRequest());
            CVolProcessedSP interpFXVol(fxAsset->getProcessedVol(atmInterp.get()));
            // cast to the type of vol we need
            volFXBS = CVolProcessedBSSP::dynamicCast(interpFXVol);
        }
        else if (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying)))
        {
            plainAsset = ProtEquityConstSP::dynamicCast(underlying)->getPlainAsset();
            FXVolBaseWrapper fxVol = ProtEquityConstSP::dynamicCast(underlying)->getFXVol();
            eqFXCorr = ProtEquityConstSP::dynamicCast(underlying)->getCorrelation();
            /// interpolate the fx vols atm 
            ATMVolRequestSP fxVolRequest(new ATMVolRequest());
            // interpolate the vol
            CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), NULL));
            // cast to the type of vol we're expecting
            volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
            underlying = plainAsset; // assign it to plain asset for the tree
        }
*/
        VolLV = CVolProcessedDVFSP::dynamicCast(
            CVolProcessedSP(underlying->getProcessedVol(&volRequest)));

        // get time metric
        timeMetric = VolLV->GetTimeMetric();

//        int iAsset = 0;
//        CVolProcessedSP volProc(mAsset->factorGetProcessedVol(iAsset, &volRequest));
//
//        VolLV = CVolProcessedDVFSP::dynamicCast(volProc);

        //need to add quanto, then can remove preSetup

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------
// allow to overwrite new VolLV.
void FD1DLV::setVolLV(CVolProcessedDVFSP newVolProcessed)
{
    VolLV = newVolProcessed;
}

//----------------------------------------------------------------

/**-----------------------------------------------------
    model initialisation  (second part)
    most memory init should be here
    finalise model initialisation once after initModel() and product init() 
    called by validate 
-----------------------------------------------------*/

void FD1DLV::finaliseModel(CControl*    control)
{
    double volForBound;
    //calc vriance used for set the grid
    //need to call it everytime for tweek
    prepare(volForBound);

    // compute fwds
    int iAsset = 0;
    //stepForward[iAsset] = DoubleArray(timeLine->NumOfStep+1, 0.0);
    dynamic_cast<const CAsset *>(
        factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward); //stepForward[iAsset][step], iAsset = 1

    // call parent last
    FD1D::finaliseModel( control );
}

//----------------------------------------------------------------

//alpha is nb of std
void FD1DLV::setFdBounds(
    double alpha,
    double & outLowB,
    double & outUpB ) const
{
    //calc vriance used for set the grid
    double volForBound;
    prepare(volForBound);

    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];

    double temp1 = alpha * volForBound * sqrt(longest_T);

    /* double fwd1 = */ underlying->fwdValue(timeLine->StepDates[timeLine->NumOfStep]);
        
    double spot = underlying->fwdValue(timeLine->StepDates[0]);

    switch (changeOfVar) {
        case X:{
            //S
            //taking min only for the purpose to put S0 at the center
            double dis = Maths::min(1.0, 0.5 * (exp(temp1) - exp(-temp1)));

            outLowB = Maths::max(0.0, spot * (1 - dis));
            outUpB = spot * (1 + dis);
        }        
        break;
        case LOG_X: {
            //log(S)            
            outLowB =  log(spot) - temp1;
            outUpB = log(spot) + temp1;
//            outLowB =  log(fwd1) - temp1;
//            outUpB = log(fwd1) + temp1;

        }
        break;
        case LOG_FWDX: {
            //log(fwd)
            outLowB =  log(spot) - temp1;
            outUpB = log(spot) + temp1;
/*
            outLowB = log(fwd1) - temp1;
            outUpB = log(fwd1) + temp1;

            //make sure that S0 is always in the interval
            double tolerance = 0.00001; 
            outUpB = Maths::max(outLowB, log(spot + tolerance));
            outLowB = Maths::min(outLowB, log(spot - tolerance));
*/
        }
        break;
    }
}

//----------------------------------------------------------------

// set up variance array 
void FD1DLV::prepare(double& volForBound) const
{
    //Variance.resize(timeLine->NumOfStep+1);
    // *** requiring ATM vol should be ok, Variance is used for setting up nodes only
    
    double strike;
    DateTime startDate = timeLine->StepDates[0];
    //DateTime startDate = TimePts.StepDates[0];
    string ccyTreatment = prod->getCcyTreatment();
    
    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
        strike = underlying->fwdValue(startDate)/fxAsset->fwdValue(startDate);
    else
        strike = underlying->fwdValue(startDate);
    
    volForBound = VolLV->computeImpVol(timeLine->StepDates[timeLine->NumOfStep], strike);
    
    /* initialize the caching of forward values for the local volatility */
    volCalculator = refCountPtr<CVolProcessedDVF::IVolCalculator>(VolLV->CreateVolCalculator(timeLine->StepDates));
}

//----------------------------------------------------------------


/** get the initial conditions and corresponding indexes in the slices
    for forward induction */
void FD1DLV::getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{
    
    static const string method("FD1DLV::getInitialConditions");
    
    if (!isFwdInduction)
    {
        throw ModelException(method, "the initial conditions shouldn't be "
                                     "required for backward induction");
    }

    initialIndex.resize(1); // resize the array of initial indexes to the number of dimensions
    initialValue.resize(1); // resize the array of initial indexes to the number of dimensions

    double undSpot1 = underlying->fwdValue(timeLine->StepDates[0]); // initial spot
    initialValue[0] = undSpot1;

    // find the indexes corresponding to the initial values
    // mode = 1, position of the smallest item larger than or equal to target
    int bot, top;
    xValue->getCalcRange( bot, top );
    initialIndex[0] = Neighbour(initialValue[0], xValue, bot, top, 1); 
}

//----------------------------------------------------------------

/**----------------------------------------------------------------
    PDE 1 factor is
    U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
----------------------------------------------------------------*/
void FD1DLV::pdeCoeff(
    int step,
    TreeSliceEQ & a,
    TreeSliceEQ & c,
    TreeSliceEQ & f,
    TreeSliceEQ & g ) const
{
    //get current spot 
    //I need only 
    int bot1, top1;
    xValue->getCalcRange( bot1, top1 );
    double * s = xValue;

    int ix;
            
    double mu = 0;
    double r = 0;

    double dt_cal;
    double dt_trading ;
    double df ;
        
    if (isFwdInduction ){
        if (step == 0){
            throw ModelException("FD1DLV::pdeCoeff", 
                "pdeCoeff shouldn't be called at time step 0!.");
        }
    
        dt_trading =timeLine->TradeYrFrac[step];

//        if (!Maths::isZero(stepForward[step-1])){
//            mu = log(stepForward[step]/stepForward[step-1]);    //[iAsset][step]
//        }else{
//            throw ModelException("FD1DLV::pdeCoeff", 
//                                 "division by zero.");
//        }


        if(stepForward[step] ==0 ){ //stepForward[iAsset][step], iAsset = 1
            // hard coded for now !!!
            throw ModelException("FD1DLV::pdeCoeff", 
            "division by zero");
        }else{
            mu = log(stepForward[step]/stepForward[step-1]); //stepForward[iAsset][step], iAsset = 1
        }

        vector<double> volst;

        //GetStepVol(step, volst, s.begin(), bot1, top1);
        GetStepVol(step, volst, s, bot1, top1);

        vector<double> dvolst;
        vector<double> d2volst;
        int size = volst.size();

        //calc first, second derivatives
        dvolst.resize(size);
        d2volst.resize(size);

        for (ix = 1; ix <= size -2; ix++){
            double hi, hiplus1, h;
            
            hi = exp(s[ix]) - exp(s[ix-1]);
            hiplus1 = exp(s[ix+1]) - exp(s[ix]);
            h = hi + hiplus1;

            dvolst[ix] = (volst[ix+1] - volst[ix-1])/(h);
            d2volst[ix] = 2.0 * (volst[ix+1] * hi  + volst[ix-1] * hiplus1 - volst[ix] * h )/(hi * hiplus1* h );
        }

        dvolst[0] = dvolst[1];
        d2volst[0] = d2volst[1];
        dvolst[size -1] = dvolst[size -2];
        d2volst[size -1] = d2volst[size -2];

        //to think when dttrading = 0, matrix may not inversable, may switch to explicite???
        //to do
        //fill all coeff of PDE
        for (ix = bot1; ix <= top1; ix++){
            
            //due to "GetStepVol"'s size = top1-bot1+1
            int ixVol =  ix - bot1;
            double halfsig2 = 0.5 * volst[ixVol] * volst[ixVol] * dt_trading;
            double sigsig1 = volst[ixVol] * dvolst[ixVol] * dt_trading;
            double sigsig2 = volst[ixVol] * d2volst[ixVol] * dt_trading;
            double sig1_2 = dvolst[ixVol] * dvolst[ixVol] * dt_trading;

            c[ix] = halfsig2;  //c*dt, coeff of U_xx
            f[ix] = (2.0 * sigsig1 + (sig1_2 + sigsig2 )* exp(s[ix]) ) *exp(s[ix]);  //f*dt, coeff of U            
            f[ix] = 0.0;  //f*dt, coeff of U
            g[ix]= 0.0;  //g*dt, 

            switch (changeOfVar) {
                case X:{//S
                }        
                break;
                case LOG_X:{//log(S)
                    a[ix]  = - ((mu - halfsig2) - 2.0 * sigsig1 * exp(s[ix]));  //a *dt, coeff of U_x
                }
                break;
                case LOG_FWDX:{//log(fwd)
                }
                break;
                default:{
                    throw ModelException("FD1DLV::pdeCoeff", "This change of variable is not available");
                }
            }
        }          
/*
        //for test the fwd price
        for (ix = bot1; ix <= top1; ix++){
            
            //due to "GetStepVol"'s size = top1-bot1+1
            int ixVol =  ix - bot1;
            double halfsig2 = 0.5 * volst[ixVol] * volst[ixVol] * dt_trading;

            c[ix] = halfsig2;  //c*dt, coeff of U_xx

            f[ix] = -r + mu;  //f*dt, coeff of U
            
            g[ix]= 0.0;  //g*dt, 

            switch (changeOfVar) {
                case 1:{//S
//                    a[ix]= mu * s[ix];  //a *dt, coeff of U_x
//                    c[ix] = halfsig2 * s[ix] * s[ix];  //c*dt, coeff of U_xx
                }        
                break;
                case 2:{//log(S)
                    a[ix]= -(mu - halfsig2);  //a *dt, coeff of U_x
                }
                break;
                case 5:{//log(fwd)
                }
                break;
                default:{
                    throw ModelException("FD1DLV::pdeCoeff", "This change of variable is not available");
                }
            }
        }          
*/
        
    }else{//backward

        if (step == timeLine->NumOfStep){
            throw ModelException("FD1DLV::pdeCoeff", 
                "pdeCoeff shouldn't be called at last step!.");
        }
    
        dt_cal = timeLine->StepDates[step].yearFrac(timeLine->StepDates[step+1]);
        dt_trading =timeLine->TradeYrFrac[step+1];

/*
        if (!Maths::isZero(stepForward[step])){
            mu = log(stepForward[step+1]/stepForward[step]);    //[iAsset][step]
        }else{
            throw ModelException("FD1DLV::pdeCoeff", 
                                 "division by zero.");
        }
*/

        if(stepForward[step] ==0 ){  //stepForward[iAsset][step], iAsset = 1
            // hard coded for now !!!
            throw ModelException("FD1DLV::pdeCoeff", 
            "division by zero");
        }else{
            mu = log(stepForward[step+1]/stepForward[step]);
        }

        df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);

        r = -log(df);

        vector<double> volst;

        //GetStepVol(step, volst, s.begin(), bot1, top1);
        GetStepVol(step, volst, s, bot1, top1);

        //fill all coeff of PDE
        for (ix = bot1; ix <= top1; ix++){
            
            //due to "GetStepVol"'s size = top1-bot1+1
            int ixVol =  ix - bot1;
            double halfsig2 = 0.5 * volst[ixVol] * volst[ixVol] * dt_trading;

            c[ix] = halfsig2;  //c*dt, coeff of U_xx

            f[ix] = -r;  //f*dt, coeff of U
            
            g[ix]= 0.0;  //g*dt, 

            switch (changeOfVar) {
                case X:{//S
                    a[ix]= mu * s[ix];  //a *dt, coeff of U_x
                    c[ix] = halfsig2 * s[ix] * s[ix];  //c*dt, coeff of U_xx
                }        
                break;
                case LOG_X:{//log(S)
                    a[ix]= mu - halfsig2;  //a *dt, coeff of U_x
                }
                break;
                case LOG_FWDX:{//log(fwd)
                    a[ix] = - halfsig2;  //a *dt, coeff of U_x
                }
                break;
                default:{
                    throw ModelException("FD1DLV::pdeCoeff", "This change of variable is not available");
                }
            }
        }          
    }
}

//----------------------------------------------------------------

/** calculate vol for current step for a set of spot levels
    returns number of vol calculated - one (flat for all node) or num */
int FD1DLV::GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end) const
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
        CSliceDouble spots(s+start, end - start + 1);
        CSliceDouble locVol(&vol[0], end - start + 1);
        
        // use the version with caching of forward values
        volCalculator->CalcLocVol(spots, 
                                  step,
                                  locVol);
    }
    return (end-start+1);
}

//----------------------------------------------------------------

bool FD1DLVLoad(){
    return (FD1DLV::TYPE && true);
}


DRLIB_END_NAMESPACE
