//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSVCJ.cpp
//
//   Description : two factor finite difference engine for SVCJ vol
//                 asset factor   dX = drift*dt + sqrt(V) dW + Merton jump in X
//                 vol factor  dV = k(V_0 - V)dt + vol*sqrt(V) dZ + Merton jump in V
//                 V_0 is allowed to be a gamma distribution of mean v0.
//
//   Author      Xiaolan Zhang
//
//   Date        : November 29, 2006
//
//----------------------------------------------------------------------------

//#include <math.h>
#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/FD2DSVCJ.hpp"
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


////////////////////////// FD2DSVCJ /////////////////////////
void FD2DSVCJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD2DSVCJ, clazz);
    SUPERCLASS(FD2DSV);
    EMPTY_SHELL_METHOD(defaultFD2DSVCJ);

    FIELD(cCrashRatedt,"commonCrashRate*dt");
	FIELD_MAKE_TRANSIENT(cCrashRatedt);

    FIELD(whichM,"which method for jumps");
    FIELD_MAKE_OPTIONAL(whichM);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD2DSVCJ::TYPE = CClass::registerClassLoadMethod(
    "FD2DSVCJ", typeid(FD2DSVCJ), load);

//----------------------------------------------------------------

FD2DSVCJ::FD2DSVCJ(const CClassConstSP & type ):
    FD2DSV(type), whichM(1){}

//----------------------------------------------------------------

FD2DSVCJ::~FD2DSVCJ(){}

//----------------------------------------------------------------

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP FD2DSVCJ::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSVCJ::TYPE->getName(),
                                               VolSurface::TYPE->getName()));
}

//----------------------------------------------------------------

IModel::WantsRiskMapping FD2DSVCJ::wantsRiskMapping() const {
    return allowRiskMapping ? riskMappingAllowed : riskMappingDisallowed;
}

//----------------------------------------------------------------
/** retrieving market data from MDF */
void FD2DSVCJ::retrieveFactor(){
    static const string method = "FD2DSVCJ::retrieveFactor";
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
        volSVCJ = VolSVCJConstSP::dynamicCast(VolRequestRaw::copyVolBase(dynamic_cast<const CAsset &>(*factors[0])));
        VolSVCJSP volSV = VolSVCJSP( IObject::copy( volSVCJ.get() ) );
        
        //remove jumps in order to convert to volSV in FD2DSV
        volSV.get()->commonCrashRate = 0.0;

        //vol in FD2DSV.
        vol = volSV->convert( const_cast< VolSV *>( vol.get() ) );
        // get time metric
        timeMetric = volSVCJ->GetTimeMetric();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

void FD2DSVCJ::finaliseModel(CControl*    control){

    // validation for SVCJ's params
    if (whichChgVarDim2 != X) {
        throw ModelException( "FD2DSVCJ::finaliseModel", 
                             "whichChgVarDim2 should be 1 (working variable is V)!");
    }

    FD2DSV::finaliseModel(control);
    
    //set the flag hasJumps
    //double lambda_c = vol->commonCrashRate;
    //if (lambda_c == 0.0 ){
        hasJumps = true;
    //}else{
        //hasJumps = true;
    //}

//    if (hasJumps) {//jump parts haven't been implemented yet
//        numCoeffPDE = 6;
//        throw ModelException( "FD2DSVCJ::finaliseModel", 
//                             "Jump part hasn't been implemented yet!");
//    }else{
//        numCoeffPDE = 6; 
//        //for the schemes used now, it's always 6, since jumps are treated by FFT
//    }
}

/** ----------------------------------------------------------------
///  setup the coefficients of the PIDE needed to be solved
///  Assume the most general PDE is 
///  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
///  this ft is called at each time step, so they're time depd.
///----------------------------------------------------------------*/

void FD2DSVCJ::pdeCoeff(int step, double*** coeff, 
                        int bot1, int top1, int bot2, int top2 ){

    //test, treat jumps explicitly,
    //need to chg
    FD2DSV::pdeCoeff(step, coeff, bot1, top1, bot2, top2);
    
    if (hasJumps = true){
        double dt_trading = timeLine->TradeYrFrac[step+1];
        cCrashRatedt = dt_trading * volSVCJ->commonCrashRate;
    }
}

//----------------------------------------------------------------

//----------------------------------------------------------------
//jump parts
//----------------------------------------------------------------

double FD2DSVCJ::jumpProDen(double x, double y){
    double fy = 0;    // the pdf of y;
    double fx_y = 0;  // the pdf of x conditional on y;
    double fxy = 0;   // the joint pdf of x and y;

    double x_mean = log(1.0 + volSVCJ->commonStockCrashSizeMean) 
        - 0.5 * Maths::square(volSVCJ->commonStockCrashSizeUncertainty) 
                    + volSVCJ->stockVolCrashSizeCorrelation * y;
    double x_var  = Maths::square(volSVCJ->commonStockCrashSizeUncertainty);
    //double x_var  *= (Vol->commonStockCrashSizeUncertainty);

    double y_mean = volSVCJ->commonVolCrashSizeMean;

    fy = y_mean * exp(y_mean * (-y));
    fx_y = exp(-(Maths::square(x - x_mean))/ (2.0*x_var))/sqrt(2.0*Maths::PI*x_var);

    fxy = fx_y * fy ;

    return fxy;
}

//temporary, to remove
double FD2DSVCJ::getJumpSizeMean(double y){
   double x_mean = log(1.0 + volSVCJ->commonStockCrashSizeMean) 
        - 0.5 * Maths::square(volSVCJ->commonStockCrashSizeUncertainty) 
                    + volSVCJ->stockVolCrashSizeCorrelation * y;
    return x_mean; 
}
  
double FD2DSVCJ::getJumeSizeVar(){
    double x_var  = Maths::square(volSVCJ->commonStockCrashSizeUncertainty);
    return x_var; 
}

//----------------------------------------------------------------

// the probability density function is provided in the write up.
/*
double FD2DSVCJ::jumpProDen(double x, double y,FD2DSolver* solver){
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
double FD2DSVCJ::jumpPro(double x, double y, FD2DSolver* solver){
    return jumpProDen(x,y,solver) * solver->dx * solver->dy;
}
*/

// It is here the jump part being computed.
// Different stochastic process models have different jump

/*void FD2DSVCJ::jumpComp(FD2DSolver *solver){
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

bool FD2DSVCJLoad(){
    return (FD2DSVCJ::TYPE && true);
}

DRLIB_END_NAMESPACE
