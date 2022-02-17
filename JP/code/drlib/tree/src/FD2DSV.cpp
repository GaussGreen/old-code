//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSV.cpp
//
//   Description : two factor finite difference engine for SVCJ vol
//                 asset factor   dX = drift*dt + sqrt(V) dW + Merton jump in X
//                 vol factor  dV = k(V_0 - V)dt + vol*sqrt(V) dZ + Merton jump in V
//                 V_0 is allowed to be a gamma distribution of mean v0.
//
//   Author      : Ning Shen
//                   Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

//#include <math.h>
#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/FD2DSV.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/VolRequestRaw.hpp"

#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"



DRLIB_BEGIN_NAMESPACE

const double DEFAULT_LBOUND2 = 0.000025; 

const double DEFAULT_P_UND2 = 0.01;


////////////////////////// FD2DSV /////////////////////////
void FD2DSV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD2DSV, clazz);
    SUPERCLASS(FD2D);
    EMPTY_SHELL_METHOD(defaultFD2DSV);
    FIELD(allowRiskMapping, "Allow risk mapping");
    FIELD_MAKE_OPTIONAL(allowRiskMapping);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD2DSV::TYPE = CClass::registerClassLoadMethod(
    "FD2DSV", typeid(FD2DSV), load);

//----------------------------------------------------------------

FD2DSV::FD2DSV(const CClassConstSP & type):
    FD2D(type),
    allowRiskMapping(false)
{
    whichChgVarDim1 = LOG_FWDX;  //1: S; 2: log(S); 5: log(fwd)
    whichChgVarDim2 = X; //1 V, 2: log(V): 3: V^(1/2); 4: V^(1/4)
}

//----------------------------------------------------------------

FD2DSV::~FD2DSV(){}

//----------------------------------------------------------------

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP FD2DSV::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSV::TYPE->getName(),
                                               VolSurface::TYPE->getName()));
}

//----------------------------------------------------------------

IModel::WantsRiskMapping FD2DSV::wantsRiskMapping() const {
    return allowRiskMapping ? riskMappingAllowed : riskMappingDisallowed;
}

//----------------------------------------------------------------
/** retrieving market data from MDF */
void FD2DSV::retrieveFactor(){
    static const string method = "FD2DSV::retrieveFactor";
    try{
        FD2D::retrieveFactor();

        underlying = CAssetConstSP::dynamicCast( factors[0] );
        if( ! underlying )
            throw ModelException( "FD2DSV::retrieveFactor", "only asset underlier suppported" );

        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD2DSV::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

        // copy the vol that's in the multiasset
        vol = VolSVSP::dynamicCast(VolRequestRaw::copyVolBase(dynamic_cast<const CAsset &>(*factors[0])));
        // get time metric
        timeMetric = vol->GetTimeMetric();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double FD2DSV::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD2DSV::getPrice0";
    try {
        const TreeSlice & s = payoffIndex->getValue( 0 );
        int bot, top;
        s.getCalcRange( bot, top );

        vector<double> v(dim2);
        computeUndLevel2(0, dim2-1, 0, &*v.begin());

        double undSpot1 = underlying->getSpot();
        double undSpot2 = vol->initialVol; 

        undSpot2 *= undSpot2;

        return interpF2(undSpot1, 
                        undSpot2, 
                        s.getValues(),
                        &*v.begin(), 
                        price.getValues(),
                        top - bot + 1,
                        dim2);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------


void FD2DSV::finaliseModel(CControl*    control){

    // validation for SVCJ's params
    if (whichChgVarDim2 != X) {
        throw ModelException( "FD2DSV::finaliseModel", 
                             "whichChgVarDim2 should be 1 (working variable is V)!");
    }

    FD2D::finaliseModel(control);
    
    // compute fwds
    int iAsset = 0;
    stepForward[iAsset] = DoubleArray(timeLine->NumOfStep+1, 0.0);
    dynamic_cast<const CAsset *>(
        factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward[iAsset]);

        //set the flag hasJumps
    //double lambda_c = vol->commonCrashRate;
    //if (lambda_c == 0.0 ){
        hasJumps = false;
    //}else{
        //hasJumps = true;
    //}

    if (hasJumps) {//jump parts haven't been implemented yet
        numCoeffPDE = 6;
        throw ModelException( "FD2DSV::finaliseModel", 
                             "Jump part hasn't been implemented yet!");


    }else{
        numCoeffPDE = 6; 
        //for the schemes used now, it's always 6, since jumps are treated by FFT
    }
}

//----------------------------------------------------------------

/**---------------------------------------------------------------*/
/** calculate the stuff which specifies the grid */
/** calculate FD grid and its underlying
/// alpha : truncation1D = nb of std
/// volForBound1 : vol used to set the FD bound for first factor
/// volForBound2 : vol used to set the FD bound for second factor
///---------------------------------------------------------------*/

void FD2DSV::setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                double& volForBound2, double alpha2, double& outLowB2, double& outUpB2){

    //volForBound2 isn't used for SVCJ
    volForBound2 = vol->initialVol;    

    //for short term, use V0
    //for long term, use meanVol to set up the boundary of stock
    volForBound1 = vol->initialVol;    
//    if (longest_T > 1/k){
//        volForBound = Vol->meanVol;
//    }
    
    //set up default boundaries for FD, can be modified by product at init()    
    //for test
    double v_LR = vol->meanVol;
    v_LR = v_LR * v_LR ;
    double k = vol->meanReversRate;
    double sigma = vol->volVol; 
    
    //maybe need to combine model and products
    //hard code for now, to change!!!!!!!!!!!!!!!!
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];

    double temp1 = alpha1* volForBound1 * sqrt(longest_T);
    double temp2 = alpha2* volForBound1 * sqrt(longest_T);

    double undSpot1 = underlying->getSpot();
    double undSpot2 = vol->initialVol; 
    undSpot2 *= undSpot2;

    switch (whichChgVarDim1) {
        case X:{
            //S
            outLowB1 = undSpot1 * exp(-temp1);
            outUpB1 = undSpot1* exp(temp1);
        }        
        break;
        case LOG_X: {
            //log(S)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
        case LOG_FWDX: {
            //log(fwd)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
    }

    //---------
    //factor 2
    //---------

    //if it's < 1/k, it's a short T,
    //we centered at V0
    //otherwise, we use inverse gamma to set up the boundary for LT.
    //but, make sure to including V0 due to current structure,
    //since we output the price for (S0, V0)
    //Maybe, don't have to include V0 if we do distribution on V0 for final price!

//for test
    if (dim2 == 1){
        //lBound2 = uBound2 = undSpot2 ;
        outUpB2=outLowB2=undSpot2;
    }else{

        if (longest_T < 1/k){
            //temp2 = truncation2D * sigma * sqrt(undSpot2) * sqrt(longest_T);
            temp2 = alpha2* sigma * sqrt(undSpot2) * sqrt(longest_T);

            outLowB2 = max(undSpot2 - temp2, DEFAULT_LBOUND2);
            outUpB2 = undSpot2 + temp2;        
        }else{ //long term
            //hard code for now
            double guessVsup = 0.2; 
            //double guessVinf = DEFAULT_LBOUND2; 
            //double p = DEFAULT_P_UND2;
            double p_bar = 1 - DEFAULT_P_UND2;
            double deltaVup = 0.01;
            //double deltaVdown = 0.0001;

            double alpha = 2.0 * k * v_LR / (sigma * sigma);
            double beta = alpha / v_LR;

            //double lBound;
            double uBound;
            double cdf;

            //obtain the low bound
    //            cdf = 0.0;
    //            while (cdf < p ){
    //                guessVinf += deltaVdown;                
    //                lBound = guessVinf * beta;
    //                //cdf = gamma_cdf(lBound, alpha);
    //                cdf = imsl_d_gamma_cdf(lBound, alpha);
    //            }
    //            LBound2 = guessVinf;
            //if I decide the low boundary based on proba, 
            //since the idea is centered around v_LR, my V0 maybe out of the grid.
            //so, I'll have a pb to report the price at S0 and V0
            //So, hard code the low bound for now
            //can be removed if we want to add the distribution of V0
            outLowB2 = DEFAULT_LBOUND2; 
          
            //obtain the upper boundary
            cdf = 0.0;
            while (cdf < p_bar ){
                guessVsup += deltaVup;                
                uBound = guessVsup * beta;
                cdf = imsl_d_gamma_cdf(uBound, alpha);                
            }
            outUpB2 = guessVsup;        

            //due to same issue as mentioned above,
            //we include this part to make sure that V0 is in the Grid
            double middle = (outLowB2 + outUpB2) / 2.0;

            if ((undSpot2 < outLowB2) || (undSpot2 > outUpB2)){

                if(undSpot2 < middle){
                    double dis = middle - undSpot2;

                    //only make sure that we have at least some points around V0    
                    outLowB2 = Maths::max(outLowB2 - 2.0 * dis, 0.0);
                }else{
                    double dis = undSpot2 - middle ;
                    
                    outUpB2 = Maths::max(outUpB2 + 2.0 * dis, 0.0);
                }
            }
        }       
    }
}

/*
void FD2DSV::SetSpaceGridData(){

    double volForBound1;
    double volForBound2;

    double trunc5_1D = truncation1D;
    double trunc3_1D = 3.0 /5.0 * truncation1D;
    double trunc2_1D = 2.0/5.0 * truncation1D;

    double trunc5_2D = truncation2D;
    double trunc3_2D = 3.0 /5.0 * truncation2D;
    double trunc2_2D = 2.0/5.0 * truncation2D;

    double lBound5_1D, uBound5_1D, lBound3_1D, uBound3_1D, lBound2_1D, uBound2_1D;
    double lBound5_2D, uBound5_2D, lBound3_2D, uBound3_2D, lBound2_2D, uBound2_2D;


    setFdBounds(volForBound1, trunc2_1D, lBound2_1D, uBound2_1D,
                volForBound2, trunc2_2D, lBound2_2D, uBound2_2D);

    setFdBounds(volForBound1, trunc3_1D, lBound3_1D, uBound3_1D,
                volForBound2, trunc3_2D, lBound3_2D, uBound3_2D);
    
    setFdBounds(volForBound1, trunc5_1D, lBound5_1D, uBound5_1D,
                volForBound2, trunc5_2D, lBound5_2D, uBound5_2D);

    
    lBound1 = lBound5_1D;
    uBound1 = uBound5_1D;
    lBound2 = lBound5_2D;
    uBound2 = uBound5_2D;

    //GetGridDataFromProd();
    

    if (isVariableGrid){
        int n = dim1 -1; //to be review

        //weight for diff truncation 
//        double w5 = 0.1;
//        double w3 = 0.1;
//        double w2 = 0.6;


        //test
        double all = uBound1 - lBound1;
        double w5 = (lBound3_1D - lBound5_1D) / all;
        double w3 = (lBound2_1D - lBound3_1D) / all;
        double w2 = 1.0 - 2*w3 - 2 * w5;

        int n5 = int (w5 * n);
        int n3 = int (w3 * n);
        int n2 = n - 2 * n5 - 2 * n3;
        int i;

        //for X
        double dxTemp = (lBound3_1D - lBound5_1D)/(n5+1);
        for( i = 1; i <= n5; i++){
            v_dxM[i] = dxTemp;            
        }

        dxTemp = (lBound2_1D - lBound3_1D)/(n3);
        for( i = n5 + 1; i <= n3 + n5; i++){
            v_dxM[i] = dxTemp;            
        }

        dxTemp = (uBound2_1D - lBound2_1D)/(n2);
        for( i = n3 + n5 + 1; i <= n2 + n3 + n5; i++){
            v_dxM[i] = dxTemp;            
        }

        dxTemp = (uBound3_1D - uBound2_1D)/(n3);
        for( i = n2 + n3 + n5 + 1; i <= n3 + n2 + n3 + n5; i++){
            v_dxM[i] = dxTemp;            
        }

        dxTemp = (uBound5_1D - uBound3_1D)/(n5);
        for( i = n3 + n2 + n3 + n5 + 1; i <= n5 + n3 + n2 + n3 + n5; i++){
            v_dxM[i] = dxTemp;            
        }


        //for Y
        double dyM = (uBound2 - lBound2) / (dim2 - 1);

        for (i=0; i< (int)v_dyM.size(); i++){
            v_dyM[i] = dyM; //initial vector of dx
        }


//        n = dim2 -1; //to be review
//
//        //weight for diff truncation 
//        w5 = 0.05;
//        w3 = 0.1;
//        w2 = 0.7;
//
//        n5 = int (w5 * n);
//        n3 = int (w3 * n);
//        n2 = n - 2 * n5 - 2 * n3;
//
//        dxTemp = (lBound3_2D - lBound5_2D)/(n5+1);
//        for( i = 1; i <= n5; i++){
//            v_dyM[i] = dxTemp;            
//        }
//
//        dxTemp = (lBound2_2D - lBound3_2D)/(n3);
//        for( i = n5 + 1; i <= n3 + n5; i++){
//            v_dyM[i] = dxTemp;            
//        }
//
//        dxTemp = (uBound2_2D - lBound2_2D)/(n2);
//        for( i = n3 + n5 + 1; i <= n2 + n3 + n5; i++){
//            v_dyM[i] = dxTemp;            
//        }
//
//        dxTemp = (uBound3_2D - uBound2_2D)/(n3);
//        for( i = n2 + n3 + n5 + 1; i <= n3 + n2 + n3 + n5; i++){
//            v_dyM[i] = dxTemp;            
//        }
//
//        dxTemp = (uBound5_2D - uBound3_2D)/(n5);
//        for( i = n3 + n2 + n3 + n5 + 1; i <= n5 + n3 + n2 + n3 + n5; i++){
//            v_dyM[i] = dxTemp;            
//        }

    }else{
        //calculate dx
        // fill with constant for now, so only do once.
        double dxM = (uBound1 - lBound1) / (dim1 - 1);
        double dyM = (uBound2 - lBound2) / (dim2 - 1);
        int i;

        for (i=0; i< (int)v_dxM.size(); i++){
            v_dxM[i] = dxM; //initial vector of dx
        }

        for (i=0; i< (int)v_dyM.size(); i++){
            v_dyM[i] = dyM; //initial vector of dx
        }
    }


    //if the grid is variable, set up the v_dyM
//    if (isVariableGrid){
//        // need to added in a more general way        
//        throw ModelException("FD2DSV::FD2DSV", "put in the comments, need to think it and rewrite it!");
//    }

}
*/


/** ----------------------------------------------------------------
///  setup the coefficients of the PIDE needed to be solved
///  Assume the most general PDE is 
///  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
///  this ft is called at each time step, so they're time depd.
///----------------------------------------------------------------*/


void FD2DSV::pdeCoeff(int step, double*** coeff, 
                        int bot1, int top1, int bot2, int top2 ){

    double v0 = vol->initialVol; 
    v0 = v0 *v0; //transfer to inititial variance
    double v_LR = vol->meanVol;
    v_LR = v_LR * v_LR; // transfer to mean variance
    double sigma = vol->volVol;
    double sigma2 = sigma * sigma;
    double rho = vol->correlation;
    double k = vol->meanReversRate;
    
    int ix;
    int iy;

    if ((!isFwdInduction && step == timeLine->NumOfStep) ||
        (isFwdInduction && step == 0))
    {
        throw ModelException("FD2DSV::pdeCoeff", 
                             "pdeCoeff shouldn't be called at last step!");
    }
    
    double mu = 0;    
    double df;
    double r = 0;    
    double dt_trading;

    // backward induction
    if (!isFwdInduction)
    {
        dt_trading = timeLine->TradeYrFrac[step+1];

        if (!Maths::isZero(stepForward[0][step])){
            mu = log(stepForward[0][step+1]/stepForward[0][step]);    //[iAsset][step]
        }else{
            throw ModelException("FD2DSV::pdeCoeff", 
                                 "division by zero.");
        }
        

//        if(stepForward1[step] ==0 )
//        {
//            // hard coded for now !!!
//            throw ModelException("FD2DSV::pdeCoeff", 
//                                 "division by zero.");
//        }
//        else
//        {
//            mu = log(stepForward1[step+1]/stepForward1[step]);    
//        }

        df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
        r = -log(df);

        // fpr accessing und level, need review !!!
        //payoffIndex->update(step, FDProduct::BWD);
        double* und = payoffIndex->getValue( step ).getValues();

        //fill all coeff of PDE
        for (ix = bot1 ; ix <= top1 ; ix++)
        {
            for (iy = bot2 ; iy <= top2 ; iy++)
            {

                coeff[1][ix][iy] = k * (v_LR - gridLevel2[iy]) * dt_trading; //b * dt, coeff of U_y
                coeff[2][ix][iy] = 0.5 * gridLevel2[iy] * dt_trading; //c * dt, coeff of U_xx 
                coeff[3][ix][iy] = 0.5 * gridLevel2[iy] * sigma2 * dt_trading; //d * dt, coeff of U_yy
                coeff[4][ix][iy] = rho * sigma * gridLevel2[iy] * dt_trading; //e * dt, coeff of U_xy 
                coeff[5][ix][iy] = -r; //f * dt, coeff of U
                
                switch (whichChgVarDim1) {
                    case X:{ //S
                        coeff[0][ix][iy] = mu * und[ix]; //a * dt, coeff of U_x 
                        coeff[2][ix][iy] *= und[ix] * und[ix]; //c * dt, coeff of U_xx 
                        coeff[4][ix][iy] *= und[ix]; //e * dt, coeff of U_xy 
                    }        
                    break;
                    case LOG_X:{ //log(S)
                        coeff[0][ix][iy] = mu - 0.5 * gridLevel2[iy] * dt_trading; //a* dt, coeff of U_x 
                    }
                    break;
                    case LOG_FWDX:{ //log(fwd)
                        coeff[0][ix][iy] = - 0.5 * gridLevel2[iy] * dt_trading; //a * dt, coeff of U_x 
                    }
                    break;
                    default:{
                        throw ModelException("FD2DSV::pdeCoeff", "This change of variable is not available");
                    }
                }
            }
        }
    }
    // forward induction
    if (isFwdInduction)
    {
        dt_trading = timeLine->TradeYrFrac[step];


        if (!Maths::isZero(stepForward[0][step-1])){
            mu = log(stepForward[0][step]/stepForward[0][step-1]);    //[iAsset][step]
        }else{
            throw ModelException("FD2DSV::pdeCoeff", 
                                 "division by zero.");
        }


//        if (stepForward1[step-1] == 0)
//        {
//            // hard coded for now !!!
//            throw ModelException("FD2DSV::pdeCoeff", 
//                                 "division by zero.");
//        }
//        else
//        {
//            mu = log(stepForward1[step]/stepForward1[step-1]);    
//        }

        df = discYC->pv(timeLine->StepDates[step-1], timeLine->StepDates[step]);
        r = -log(df);

        // fpr accessing und level, need review !!!
        //int bot, top, dummy;
        //payoffIndex->update(step, FDProduct::FWD);
        //const vector< TreeSliceSP > & spots =
        //    payoffIndex->getValue(step, bot, top, dummy, dummy);

        //fill all coeff of PDE
        for (ix = bot1 ; ix <= top1 ; ix++)
        {
            for (iy = bot2 ; iy <= top2 ; iy++)
            {
                switch (whichChgVarDim1) {
                    case X:{ // S
                        coeff[0][ix][iy] = (gridLevel1[ix] * (2.0 * gridLevel2[iy] + rho * sigma - (mu / dt_trading))) * dt_trading; //a * dt, coeff of U_x
                        coeff[1][ix][iy] = (sigma2 + rho * sigma * gridLevel2[iy] - k * (v_LR - gridLevel2[iy])) * dt_trading; //b * dt, coeff of U_y
                        coeff[2][ix][iy] = (0.5 * gridLevel1[ix] * gridLevel1[ix] * gridLevel2[iy]) * dt_trading; //c * dt, coeff of U_xx
                        coeff[3][ix][iy] = (0.5 * sigma2 * gridLevel2[iy]) * dt_trading; //d * dt, coeff of U_yy
                        coeff[4][ix][iy] = (rho * sigma * gridLevel1[ix] * gridLevel2[iy]) * dt_trading; //e * dt, coeff of U_xy
                        coeff[5][ix][iy] = (gridLevel2[iy] + rho * sigma - (mu / dt_trading) + k) * dt_trading; //f * dt, coeff of U
                    }
                    case LOG_X: { // log(S)
                        coeff[0][ix][iy] = (rho * sigma - (mu / dt_trading) + 0.5 * gridLevel2[iy]) * dt_trading; //a * dt, coeff of U_x
                        coeff[1][ix][iy] = (sigma2 - k * (v_LR - gridLevel2[iy])) * dt_trading; //b * dt, coeff of U_y
                        coeff[2][ix][iy] = (0.5 * gridLevel2[iy]) * dt_trading; //c * dt, coeff of U_xx
                        coeff[3][ix][iy] = (0.5 * sigma2 * gridLevel2[iy]) * dt_trading; //d * dt, coeff of U_yy
                        coeff[4][ix][iy] = (rho * sigma * gridLevel2[iy]) * dt_trading; //e * dt, coeff of U_xy
                        coeff[5][ix][iy] = (k) * dt_trading; //f * dt, coeff of U
                    }
                    break;
                    case LOG_FWDX:{ // log(Fwd)
                        coeff[0][ix][iy] = (rho * sigma + 0.5 * gridLevel2[iy]) * dt_trading; //a * dt, coeff of U_x
                        coeff[1][ix][iy] = (sigma2 - k * (v_LR - gridLevel2[iy])) * dt_trading; //b * dt, coeff of U_y
                        coeff[2][ix][iy] = (0.5 * gridLevel2[iy]) * dt_trading; //c * dt, coeff of U_xx
                        coeff[3][ix][iy] = (0.5 * sigma2 * gridLevel2[iy]) * dt_trading; //d * dt, coeff of U_yy
                        coeff[4][ix][iy] = (rho * sigma * gridLevel2[iy]) * dt_trading; //e * dt, coeff of U_xy
                        coeff[5][ix][iy] = (k) * dt_trading; //f * dt, coeff of U
                    }
                    break;
                    default:{
                        throw ModelException("FD2DSV::pdeCoeff", "change of variable not available.");
                    }
                }
            }
        }
    }
}

/** get the initial conditions and corresponding indexes in the slices
    for forward induction */
void FD2DSV::getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{
    
    static const string method("FD2DSV::getInitialConditions");
    
    if (!isFwdInduction)
    {
        throw ModelException(method, "the initial conditions shouldn't be "
                                     "required for backward induction");
    }
    
    const TreeSlice & s = payoffIndex->getValue( 0 );
    int botS, topS;
    s.getCalcRange( botS, topS );

    vector<double> v(dim2);
    computeUndLevel2(0, dim2 - 1, 0, &*v.begin());

    initialIndex.resize(2); // resize the array of initial indexes to the number of dimensions
    initialValue.resize(2); // resize the array of initial indexes to the number of dimensions

    double undSpot1 = underlying->getSpot(); // initial spot
    initialValue[0] = undSpot1;
    double undSpot2 = vol->initialVol; 
    undSpot2 *= undSpot2; // initial variance
    initialValue[1] = undSpot2;

    // find the indexes corresponding to the initial values
    // mode = 1, position of the smallest item larger than or equal to target
    initialIndex[0] = Neighbour(initialValue[0], s.getValues(), botS, topS, 1);
    initialIndex[1] = Neighbour(initialValue[1], &*v.begin(), 0, dim2 - 1, 1);
}

//----------------------------------------------------------------

/*
void FD2DSV::pdeCoeff(int step, double*** coeff, int nbOfCoeffs, const FD2DSolver* solver){
 
//This ft can't be called at maturity due to timeLine
//Assuming that DimX is stock and DimY is Variance, or sqrt(Variance), or (Var)^(1/4)

    double v0 = vol->initialVol; 
    v0 = v0 *v0; //transfer to init variance
    double sigma = vol->volVol; 
    double rho   = vol->correlation;
    double v_LR = vol->meanVol;
    v_LR = v_LR * v_LR ;
    double k = vol->meanReversRate;
    
    int ix;
    int iy;

    int bot1 = 0;
    int top1= dim1-1;
    int bot2 = 0;
    int top2= dim2-1;

    //to update
    if (numCoeffPDE !=6 && nbOfCoeffs != 6){
        throw ModelException("FD2DSV::pdeCoeff", 
            "solver does not support 9 coefficients.");
    }

    if (step == timeLine->NumOfStep){
        throw ModelException("FD2DSV::pdeCoeff", 
            "PDE coeff ft can't be called at timeLine->NumOfStep!");
    }
        
    // make sure  it's not called at T
    double dt_cal = timeLine->StepDates[step].yearFrac(timeLine->StepDates[step+1]);
    double dt_trading =timeLine->TradeYrFrac[step+1];
    
    if (dt_trading ==0 || stepForward1[step] ==0 || dt_cal == 0){
        // hard coded for now !!!
        throw ModelException("FD2DSV::pdeCoeff", 
        "division by zero");
    }
    
    double mu = log(stepForward1[step+1]/stepForward1[step]);    
    mu = mu/dt_trading;

    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);

    double r = -log(df)/(dt_trading);




    double sigma2 = sigma * sigma;

    // fpr accessing und level, need review !!!
    int bot, top, dummy;
    payoffIndex->update(step, false);
    const vector< TreeSliceSP > & spots =
        payoffIndex->getValueArr(step, bot, top, dummy, dummy);
    double* und = spots[0];

    switch (whichChgVarDim2) {
        case 1:{
            //V
            //fill all coeff of PDE
            for (ix = bot1; ix <= top1; ix++){
                for (iy = bot2; iy <= top2; iy++){

                    solver->coeff[1][ix][iy] = k * (v_LR - gridLevel2[iy]);  //b

                    solver->coeff[2][ix][iy] = 0.5 * gridLevel2[iy];  //c

                    solver->coeff[3][ix][iy] = 0.5 * gridLevel2[iy] * sigma2;  //d

                    solver->coeff[4][ix][iy] = rho * sigma * gridLevel2[iy];  //e


                    if (hasDiscTerminPDE == true){
                        solver->coeff[5][ix][iy] = -r;  //f
                    }else{
                        solver->coeff[5][ix][iy] = 0.0;  //f
                    }

                    switch (whichChgVarDim1) {
                        case 1:{//S
                            solver->coeff[0][ix][iy] = mu * und[ix];  //a
                            solver->coeff[2][ix][iy] *= und[ix] * und[ix];  //c
                            solver->coeff[4][ix][iy] *= und[ix];  //e
                        }        
                        break;
                        case 2:{//log(S)
                            solver->coeff[0][ix][iy] = mu - 0.5 * gridLevel2[iy];  //a
                        }
                        break;
                        case 5:{//log(fwd)
                            solver->coeff[0][ix][iy] = - 0.5 * gridLevel2[iy];  //a
                        }
                        break;
                        default:{
                            throw ModelException("NewFD1FLV::pdeCoeff", "This change of variable is not available");
                        }
                    }

                }
            }
        }        
        break;
        case 3: {
           //V^(1/2)
            //fill all coeff of PDE
            for (ix = bot1; ix <= top1; ix++){
                for (iy = bot2; iy <= top2; iy++){
                    double y = gridLevel2[iy];
                    double y2 = y * y;

                    //solver->coeff[0][ix][iy] = mu - 0.5 * y2;  //a

                    solver->coeff[1][ix][iy] = 0.5 * (k * (v_LR - y2) - 0.25 * sigma2) / y ;  //b

                    solver->coeff[2][ix][iy] = 0.5 * y2;  //c

                    solver->coeff[3][ix][iy] = sigma2 / 8.0;  //d

                    solver->coeff[4][ix][iy] = 0.5 * rho * sigma * y;  //e

                    //solver->coeff[5][ix][iy] = -r;  //f

                    if (hasDiscTerminPDE == true){
                        solver->coeff[5][ix][iy] = -r;  //f
                    }else{
                        solver->coeff[5][ix][iy] = 0.0;  //f
                    }

                    //to optm, put it for test now
                    //log(fwd) on stock
                    if (whichChgVarDim1 == 5) {
                        solver->coeff[0][ix][iy] = - 0.5 * y2;  //a
                    }else if (whichChgVarDim1 == 2){
                        solver->coeff[0][ix][iy] = mu - 0.5 * y2;  //a
                    }else{
                        throw ModelException("NewFD1FLV::pdeCoeff", "This change of variable is not available");
                    }
                }
            }
        }
        break;
        case 4: {
           //V^(1/4)
            for (ix = bot1; ix <= top1; ix++){
                for (iy = bot2; iy <= top2; iy++){
                    double y = gridLevel2[iy];
                    double y2 = y*y;
                    
                    double y3 = y*y*y;

                    double y4 = y*y*y*y;

                    double sigma2 = sigma*sigma;


                    if (whichChgVarDim1 == 2) {
                    
                        solver->coeff[0][ix][iy] = mu - 0.5 * y4;  //a

                        solver->coeff[1][ix][iy] = 0.25 * (k * (v_LR - y4) - (3.0 / 8.0) * sigma2) / y3 ;  //b

                        solver->coeff[2][ix][iy] = 0.5 * y4;  //c

                        solver->coeff[3][ix][iy] = sigma2 / (32.0 * y2);  //d

                        solver->coeff[4][ix][iy] = 0.25 * rho * sigma * y2;  //e

                        solver->coeff[5][ix][iy] = -r;  //f
                    }else{
                        throw ModelException("NewFD1FLV::pdeCoeff", "This change of variable is not available");
                    }
                }
            }
        }
        break;

    }

}
*/
//----------------------------------------------------------------
//jump parts
//----------------------------------------------------------------

/*double FD2DSV::jumpProDen(double x, double y){
    double fy = 0;    // the pdf of y;
    double fx_y = 0;  // the pdf of x conditional on y;
    double fxy = 0;   // the joint pdf of x and y;

    double x_mean = log(1+ vol->commonStockCrashSizeMean) 
        - 0.5 * Maths::square(vol->commonStockCrashSizeUncertainty) 
                    + vol->stockVolCrashSizeCorrelation * y;
    double x_var  = Maths::square(vol->commonStockCrashSizeUncertainty);
    //double x_var  *= (Vol->commonStockCrashSizeUncertainty);

    double y_mean = vol->commonVolCrashSizeMean;

    fy = y_mean * exp(y_mean * (-y));
    fx_y = exp(-Maths::square(x - x_mean)/ (2*x_var))/sqrt(2*Maths::PI*x_var);

    fxy = fx_y * fy ;

    return fxy;
}
*/
//----------------------------------------------------------------

// the probability density function is provided in the write up.
/*
double FD2DSV::jumpProDen(double x, double y,FD2DSolver* solver){
    double fy = 0;    // the pdf of y;
    double fx_y = 0;  // the pdf of x conditional on y;
    double fxy = 0;   // the joint pdf of x and y;

    double x_mean = log(1+ solver->commonStockCrashSizeMean) 
                    - 0.5 * pow(solver->commonStockCrashSizeUncertainty, 2) 
                    - solver->stockVolCrashSizeCorrelation * y;
    double x_var  = pow(solver->commonStockCrashSizeUncertainty, 2);

    double y_mean = solver->commonVolCrashSizeMean;

    fy = y_mean * exp(y_mean * (-y));
    fx_y = exp(-pow(x - x_mean, 2)/ (2*x_var))/sqrt(2*Maths::PI*x_var);

    fxy = fx_y * fy ;

    return fxy;
}
*/


    
// cumulative prob ~ prob density * interval * interval.
/*
double FD2DSV::jumpPro(double x, double y, FD2DSolver* solver){
    return jumpProDen(x,y,solver) * solver->dx * solver->dy;
}
*/

// It is here the jump part being computed.
// Different stochastic process models have different jump

/*void FD2DSV::jumpComp(FD2DSolver *solver){
    double fxy = 0;
    double x = 0;
    double y = 0;
    int index = 0;
    int xNum = solver->xNum;
    int yNum = solver->yNum;

    int xMiddle = solver->xMiddle;
    int yMiddle = solver->yMiddle;

    double dx = solver->dx;
    double dy = solver->dy;

    double xCenter = solver->xCenter;
    double yCenter = solver->yCenter;
    
// Implement the jump part using the FFT.
    int idx, idy = 0; // the index for the jump probability array

    for (idx = 0; idx < xNum; idx++){
        for (idy = 0; idy < yNum; idy++){
            x = (double)(idx - xMiddle) * dx + xCenter;
            y = (double)(idy - yMiddle) * dy + yCenter;
    
            solver->jumpProb_slice[solver->Pos2Index(idx,idy, xNum, yNum)] 
                = jumpPro(x,y,solver);
        }
    }    
    solver->Conv2(xNum, yNum, solver->jumpProb_slice, 
                    solver->current_slice, solver->jumpPart);
}

*/

bool FD2DSVLoad(){
    return (FD2DSV::TYPE && true);
}

DRLIB_END_NAMESPACE
