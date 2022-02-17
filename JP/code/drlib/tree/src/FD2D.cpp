//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2D.cpp
//
//   Description : two factor finite difference engine
//
//   Author      : Ning Shen
//               : Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD2D.hpp"
#include "edginc/FD2DSolver.hpp"
#include "edginc/Maths.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

const double DefaultTruncation = 4.0; // 4 stdev truncation
const int DefaultDim1 = 51; // 51 point for dimension 1 (x direction)
const int DefaultDim2 = 51; // 51 point for dimension 2 (y direction)

//----------------------------------------------------------------

void FD2D::load(CClassSP& clazz){
    REGISTER(FD2D, clazz);
    SUPERCLASS(FDModel);
    clazz->setPublic(); // make visible to EAS/spreadsheet
    FIELD(truncation1D, "num of stdev to truncate for 1st dimention");
    FIELD(truncation2D, "num of stdev to truncate for 2nd dimention");
    FIELD(dim1, "num of point in FD grid for 1st dimension");
    FIELD(dim2, "num of point in FD grid for 2nd dimension");
    FIELD(stepsPerYear, "num of time points per year");

    //optional, for test purpose

    FIELD(whichChgVarDim1String, "For DR use only : 1: X, 2: log(X), 5: log(fwd)");
    FIELD_MAKE_OPTIONAL(whichChgVarDim1String);

    FIELD(whichChgVarDim2String, "For DR use only : 1: Y, 2: log(Y), 3: sqrt(Y), 4: Y^(1/4)");
    FIELD_MAKE_OPTIONAL(whichChgVarDim2String);

    FIELD(isVariableGrid, "For DR use only : true, variable grid in factor 2 ");
    FIELD_MAKE_OPTIONAL(isVariableGrid);

    FIELD(solveMethodString, "For DR use only : 1: LOD or ADI, 2: Advanced ADI");
    FIELD_MAKE_OPTIONAL(solveMethodString);
}

//----------------------------------------------------------------

CClassConstSP const FD2D::TYPE = CClass::registerClassLoadMethod(
    "FD2D", typeid(FD2D), load);

//----------------------------------------------------------------

FD2D::FD2D(CClassConstSP clazz) : LatticeModelEDR(clazz){
    // data initialised once and should not change
    truncation1D = DefaultTruncation;
    truncation2D = DefaultTruncation;
    dim1 = DefaultDim1;
    dim2 = DefaultDim2;
    sameGridTweak = false; // this must be false here !!!
    numCoeffPDE = 6;  //general PDE 2F has 6 coeffs, jumps model will be more.      
    hasJumps = false;
    isFwdInduction = false; //only backward for now
 
    //optional input
    whichChgVarDim1 = LOG_X;   
    whichChgVarDim2 = X;     
    solveMethod = ADVANCE_ADI; //1: ADI 2: Advanced ADI

    whichChgVarDim1String = "LOG_X";   
    whichChgVarDim2String = "X";      
    solveMethodString = "ADVANCE_ADI"; //1: ADI 2: Advanced ADI


    isVariableGrid = false;

    needRecalcUndValue1 = false;
    needRecalcUndValue2 = false;

    //default value for local variables
    botDim1= 0;
    topDim1= 0;
    botDim2= 0;
    topDim2= 0;

    nProd = -1;
    maxNumOfValue = -1;

    //barriers
    needSpecialFD = false;
}

//----------------------------------------------------------------

FD2D::~FD2D(){
     
    int j;

    if (botDim1 != 0){

        for (j=0; j < nProd; j++){
            // for each product
            delete [] botDim1[j];
            delete [] topDim1[j];
            delete [] botDim2[j];
            delete [] topDim2[j];
        }
        delete [] botDim1;
        delete [] topDim1;
        delete [] botDim2;
        delete [] topDim2;
    }

    botDim1 = 0;
    topDim1 = 0;
    botDim2 = 0;
    topDim2 = 0;
}

//----------------------------------------------------------------

// create a solver
FDModel::IFDSolverSP FD2D::createSolver(){
    return IFDSolverSP(new FD2DSolver(this));
}

//----------------------------------------------------------------

void FD2D::validatePop2Object()
{
    static const string method = "FD2D::validatePop2Object";
    try {

        FDModel::validatePop2Object();

        // param validation
        if (truncation1D < 1.0 || truncation1D > 10.0
            || truncation2D < 1.0 || truncation2D > 10.0) {
            throw ModelException(method, 
                                 "fd truncation must be between 1.0 and 10.0");
        }

        //to update
        if (numCoeffPDE !=6 ){
            throw ModelException("FD2D::validate", 
                "2F PDE only supports 6 coefficients.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked after instrument has got its market data. */
void FD2D::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
    FDModel::getMarket(market, instruments);

    if( factors.size() != 1 )
        throw ModelException( "FD2D::getMarket", "only 1 factor suppported" );
}

/** accept or reject factors */
bool FD2D::acceptFactor( IMarketFactor * factor )
{
    return dynamic_cast< CAsset * >( factor ) != 0;
}

/** call timeline builder usinf initData */
//can be used to prepare vars need to set the FD boundaries, such as total variance in LV
void FD2D::initModel()
{
    try{
        stepsPerYearFD = stepsPerYear; // set FDModel stepsPerYear

        // set up time line
        FDModel::initModel();

        stepForward.resize(factors.size());

        gridLevel1.resize(dim1);
        gridLevel2.resize(dim2);
        v_dxM.resize(dim1);
        v_dyM.resize(dim2);

        //convert the registered string to enum type variables
        convertRegisteredString();
        if (whichChgVarDim1 == LOG_FWDX){
            needRecalcUndValue1 = true;
        }

        if (whichChgVarDim2 == LOG_FWDX){
            needRecalcUndValue2 = true;
        }

        //create spot range
        range = TreeSliceEQ::Range::create( 0, dim1 - 1, 0, dim2 - 1 );
    }
    catch(exception& e){
        throw ModelException(e, "FD2D::initModel()");
    }
}

//----------------------------------------------------------------

/**-----------------------------------------------------*/
/** model initialisation  (second part)
/// most memory init should be here
/// finalise model initialisation once after initModel() and product init() 
/// called by validate 
///-----------------------------------------------------*/

/** set up the payoff */
void FD2D::finaliseModel(CControl*    control){
   
    static const string method = "FD2D::validate";
    try {
        bool rebuild = isRebuilt(control);
        if(rebuild){
            int i, j;

            const FDProductArray & products = getProducts();

            nProd = products.size();

            maxNumOfValue =0;
            for (j=0; j < nProd; j++){
                maxNumOfValue = Maths::max(maxNumOfValue, products[j]->getSlicesToDEV().size());
            }

            botDim1= new int*[nProd];
            topDim1= new int*[nProd];
            botDim2= new int*[nProd];
            topDim2= new int*[nProd];

            for (j=0; j<nProd; j++){
                // for each product
                //int numOfValue = (*products)[j]->getNumOfValue(0);

                botDim1[j]= new int [maxNumOfValue + 1];
                topDim1[j]= new int [maxNumOfValue + 1];
                botDim2[j]= new int [maxNumOfValue + 1];
                topDim2[j]= new int [maxNumOfValue + 1];

                for (i = 0; i <= maxNumOfValue; i++){
                    botDim1[j][i] = 0;
                    topDim1[j][i] = dim1-1;
                    botDim2[j][i] = 0;
                    topDim2[j][i] = dim2-1;
                }
            }

            // to be reviewed !!!
            SetSpaceGridData();

            //compute the FD GridLevel and its state variables
            computeFdGridLevel(0, dim1-1, v_dxM, lBound1, uBound1, gridLevel1);
            computeFdGridLevel(0, dim2-1, v_dyM, lBound2, uBound2, gridLevel2);
        }        
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

// convert input string to enum type
void FD2D::convertRegisteredString()
{
    if (whichChgVarDim1String == "X")
        whichChgVarDim1 = X;
    else if (whichChgVarDim1String == "LOG_X")
        whichChgVarDim1 = LOG_X;
    else if (whichChgVarDim1String == "LOG_FWDX")
        whichChgVarDim1 = LOG_FWDX;
    else
        throw ModelException("FD2D::convertRegisteredString", 
                             "unknown chg var type ("+whichChgVarDim1String+").\n"
                             "Type X, LOG_X, LOG_FWDX  needed");

    if (whichChgVarDim2String == "X")
        whichChgVarDim2 = X;
    else if (whichChgVarDim2String == "LOG_X")
        whichChgVarDim2 = LOG_X;
    else if (whichChgVarDim2String == "LOG_FWDX")
        whichChgVarDim2 = LOG_FWDX;
    else
        throw ModelException("FD2D::convertRegisteredString", 
                             "unknown chg var type ("+whichChgVarDim2String+").\n"
                             "Type X, LOG_X, LOG_FWDX  needed");

    if (solveMethodString == "ADI")
        solveMethod = ADI;
    else if (solveMethodString == "ADVANCE_ADI")
        solveMethod = ADVANCE_ADI;
    else
        throw ModelException("FD2D::convertRegisteredString", 
                             "unknown chg var type ("+solveMethodString+").\n"
                             "Type ADI, ADVANCE_ADI  needed");

}

//----------------------------------------------------------------

void FD2D::SetSpaceGridData(){

    double volForBound1;
    double volForBound2;

    setFdBounds(volForBound1, truncation1D, lBound1, uBound1,
                volForBound2, truncation2D, lBound2, uBound2);

    //GetGridDataFromProd();
    
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

    //if the grid is variable, set up the v_dyM
//    if (isVariableGrid){
//        // need to added in a more general way        
//        throw ModelException("FD2DSVCJ::FD2DSVCJ", "put in the comments, need to think it and rewrite it!");
//    }

}

//----------------------------------------------------------------
/**   fd stock level can call CDF mapping here */

void FD2D::computeFdGridLevel (int bot, int top, vector<double>& vdx, 
            double lBoundary, double uBoundary, vector<double>& s){

    static const string method = "FD2D::computeFdGridLevel";
    try {
        int j;

        //need to calc cum vdx
        vector<double > vcdx(vdx.size());

        vcdx[0] = 0;
        for (int i = 1; i < (int) vdx.size(); i++){
            vcdx[i] = vcdx[i-1] + vdx[i];
        }

        for (j = bot; j<=top; j++){
            s[j] = lBoundary + vcdx[j];
        }              
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------
/** update payoffIndex          */

void FD2D::payoffIndexUpdate (int& step, FDProduct::UpdateType type){
    payoffIndex->update(step, type);
}

//----------------------------------------------------------------

void FD2D::computeUndLevel1(int bot, int top,  int step, double*s) const{


    static const string method = "FD2D::computeUndLevel1";
    try {
        int j;

        //need to calc cum vdx
        vector<double > vcdx(v_dxM.size());

        vcdx[0] = 0;
        for (int i = 1; i < (int) v_dxM.size(); i++){
            vcdx[i] = vcdx[i-1] + v_dxM[i];
        }

        switch (whichChgVarDim1) {
            case X:{
                //S
                for (j = bot; j<=top; j++){
                    s[j] = lBound1 + vcdx[j];
                }
            }        
            break;
            case LOG_X: {
               //log(S)
                for (j = bot; j<=top; j++){
                    s[j] = exp(lBound1 + vcdx[j]);  
                }
            }
            break;
            case LOG_FWDX:{
                //log(fwd)


                double mu;
                if (!Maths::isZero(stepForward[0][0])){
                    mu = stepForward[0][step]/stepForward[0][0];    //[iAsset][step]
                }else{
                    throw ModelException(method, 
                                         "Forward of asset 1 is 0 at step = 0.");
                }

                for (j = bot; j<=top; j++){
                    s[j] = mu * exp(lBound1 + vcdx[j]);
                }
            }
            break;
        }     
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

//s: v in SVCJ
void FD2D::computeUndLevel2(int bot, int top, int step, double*s) const{


    static const string method = "FD2D::computeUndLevel2";
    try {
        int j;

        //need to calc cum vdx
        vector<double > vcdy(v_dyM.size());

        vcdy[0] = 0;
        for (int i = 1; i < (int) v_dyM.size(); i++){
            vcdy[i] = vcdy[i-1] + v_dyM[i];
        }

        switch (whichChgVarDim2) {
            case X:{
                //V
                for (j = bot; j<=top; j++){
                    s[j] = lBound2 + vcdy[j];
                }
            }        
            break;
            case LOG_X: {
               //log(V)
                for (j = bot; j<=top; j++){
                    s[j] = exp(lBound2 + vcdy[j]);  
                }
            }
            break;
            case LOG_FWDX:{
                //log(fwd)

                double mu;
                if (!Maths::isZero(stepForward[1][0])){
                    mu = stepForward[1][step]/stepForward[1][0];    //[iAsset][step]
                }else{
                    throw ModelException(method, 
                                         "Forward of asset 2 is 0 at step = 0.");
                }

                for (j = bot; j<=top; j++){
                    s[j] = mu * exp(lBound2 + vcdy[j]);
                }
            }
        }                              
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

void FD2D::Spot::Product::update( int & step, FDProduct::UpdateType type )
{
    // we assume just need the first und level for spot here for now
    const TreeSlice & spot = getValue( step );
    int bot, top;
    spot.getCalcRange( bot, top );
    static_cast< FD2D * >( model )->computeUndLevel1( bot, top, step, spot.getValues() );
}

//----------------------------------------------------------------

void FD2D::Spot2::Product::update( int & step, FDProduct::UpdateType type )
{
    // we assume just need the first und level for spot here for now
    const TreeSlice & spot = getValue( step );
    int bot, top;
    spot.getCalcRange( bot, top );
    static_cast< FD2D * >( model )->computeUndLevel2( bot, top, step, spot.getValues() );
}

//----------------------------------------------------------------

FDProductSP FD2D::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            return payoffIndex = FDModel::makeProduct( IProdCreatorSP( new Spot( spec->getFactor()->getName() ) ) );
        }
        else
        {
            throw ModelException( "FD2D::makeProduct",
                "IndexSpec source is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' is expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

bool FD2DLoad(){
    return (FD2D::TYPE !=0);
}

DRLIB_END_NAMESPACE
