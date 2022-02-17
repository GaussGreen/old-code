//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRet.cpp
//
//   Description : 1 factor finite difference engine
//
//   Author      : Xiaolan Zhang
//
//   Date        : Apr, 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FD1DRet.hpp"
#include "edginc/FD1DRetSolver.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

const double DefaultTruncation = 4.0; // 4 stdev truncation
const int DefaultDim1 = 50; // 50 point for dimension 1

const int minNe = 20;
const double oneD = 2.0/360.0;

//----------------------------------------------------------------

void FD1DRet::load(CClassSP& clazz){

    REGISTER(FD1DRet, clazz);
    SUPERCLASS(FDModel);

    FIELD(truncation1D, "num of stdev to truncate for 1st dimention");
    FIELD(dim1, "num of point in FD grid for 1st dimension");
    FIELD(stepsPerYear, "num of time points per year");

    //optional, for test purpose
    FIELD(whichChgVarDim1String, "X: underlying, LOG_X: log(underlying), LOG_FWDX: log(fwd underlying)");
    FIELD_MAKE_OPTIONAL(whichChgVarDim1String);

//    FIELD(hasDiscTerminPDE, "For DR use only : true term -ru in PDE, false chg var. no more -ru in PDE, only for equity");
//    FIELD_MAKE_OPTIONAL(hasDiscTerminPDE);

    FIELD(solveMethodString, "For DR use only : DEFAULT");
    FIELD_MAKE_OPTIONAL(solveMethodString);

    FIELD(whichBarrierMethodMString, "For DR use only, special barrier scheme on space: FIX_GRID: FD grid isn't changed for any t , VAR_GRID: FD grid change when there is barrier ");
    FIELD_MAKE_OPTIONAL(whichBarrierMethodMString);

    FIELD(minStepInAddedSegM, "For DR use only, special barrier scheme on space: min time step for the added segment around critical dates such as barrier");
    FIELD_MAKE_OPTIONAL(minStepInAddedSegM);

    FIELD(whichSetSpaceStepMString, "For DR use only, set var space dx: NONE: no special treatment on space; STD: space step set based on std; STD_BAR: STD + more pts around barriers");
    FIELD_MAKE_OPTIONAL(whichSetSpaceStepMString);

//    FIELD(isVariableGrid, "For DR use only : true, variable grid ");
//    FIELD_MAKE_OPTIONAL(isVariableGrid);

    FIELD(DEBUG_NODE_GEN, "For DR use only : true, false: no more than 2 barriers; false: more than 2 barriers ");
    FIELD_MAKE_OPTIONAL(DEBUG_NODE_GEN);

    FIELD(needSpecialFD, "For DR use only : true, special barrier ");
    FIELD_MAKE_OPTIONAL(needSpecialFD);

    FIELD(DEBUG_DumpToFile, "file name to dump slices to");
    FIELD_MAKE_OPTIONAL(DEBUG_DumpToFile);
}

//----------------------------------------------------------------

// convert input string to enum type
void FD1DRet::convertRegisteredString()
{
    if (whichChgVarDim1String == "X")
        whichChgVarDim1 = X;
    else if (whichChgVarDim1String == "LOG_X")
        whichChgVarDim1 = LOG_X;
    else if (whichChgVarDim1String == "LOG_FWDX")
        whichChgVarDim1 = LOG_FWDX;
    else
        throw ModelException("FD1DRet::convertRegisteredString", 
                             "unknown chg var type ("+whichChgVarDim1String+").\n"
                             "Type X, LOG_X, LOG_FWDX  needed");


    if (whichBarrierMethodMString == "FIX_GRID")
        whichBarrierMethodM = FIX_GRID;
    else if (whichBarrierMethodMString == "VAR_GRID")
        whichBarrierMethodM = VAR_GRID;
    else
        throw ModelException("FD1DRet::convertRegisteredString", 
                             "unknown chg var type ("+whichBarrierMethodMString+").\n"
                             "Type FIX_GRID, VAR_GRID  needed");

    if (whichSetSpaceStepMString == "NONE")
        whichSetSpaceStepM = NONE;
    else if (whichSetSpaceStepMString == "STD")
        whichSetSpaceStepM = STD;
    else if (whichSetSpaceStepMString == "STD_BAR")
        whichSetSpaceStepM = STD_BAR;
    else
        throw ModelException("FD1DRet::convertRegisteredString", 
                             "unknown chg var type ("+whichSetSpaceStepMString+").\n"
                             "Type DEFAULT, STD, STD_BAR  needed");


    if (solveMethodString == "DEFAULT")
        solveMethod = DEFAULT;
    else
        throw ModelException("FD1DRet::convertRegisteredString", 
                             "unknown chg var type ("+solveMethodString+").\n"
                             "DEFAULT  needed");

}

//----------------------------------------------------------------

// helpers
CClassConstSP const FD1DRet::TYPE = CClass::registerClassLoadMethod(
    "FD1DRet", typeid(FD1DRet), load);

//----------------------------------------------------------------

////////////////////////// FD1DRet /////////////////////////

FD1DRet::FD1DRet(CClassConstSP clazz) : LatticeModelEDR(clazz)
{
    // data initialised once and should not change
    truncation1D = DefaultTruncation;
    dim1 = DefaultDim1;
    numCoeffPDE = 4;  //general PDE 2F has 6 coeffs, jumps model will be more.  
    sameGridTweak = false; // this must be false here !!!
    hasJumps = false;
    isFwdInduction = false; //only backward for now

    //optional input
    whichChgVarDim1 = LOG_X;   //X: S; LOG_X: log(S); LOG_FWDX: log(fwd)
    whichChgVarDim1String = "LOG_X";   //X: S; LOG_X: log(S); LOG_FWDX: log(fwd)

    solveMethod = DEFAULT;
    solveMethodString = "DEFAULT";

    hasDiscTerminPDE = true;

    //default value for local variables
    needRecalcUndValue1 = false;
    
    //special barrier
    barrier = 0;

    whichBarrierMethodMString = "VAR_GRID";
    whichBarrierMethodM = VAR_GRID; //FIX_GRID: don't chg the grid between 0 to T
                                    // VAR_FIX: chg fd grid based on diff barrier levels
    minStepInAddedSegM = 15;

    whichSetSpaceStepMString = "NONE";
    whichSetSpaceStepM = NONE;   //NONE: don't use critical pts on space for setting theFD Grid
                                    //STD: set space steps based on diff std
                                    //STD_BAR: set space steps based on diff std and diff barrier levels

    //default value for local variables
    needSpecialFD = false;
    isVariableGrid = false;

    botDim1= 0;
    topDim1= 0;

    //for now
    nProd = -1;
    maxNumOfValue = -1;

    numOfInsertNodeM = 0; //only allow 2 inserted node at this moment

    DEBUG_NODE_GEN = false;
    DEBUG_forceNonNegative = false;
    DEBUG_UseCtrlVar = false;
}

//----------------------------------------------------------------

FD1DRet::~FD1DRet(){

    if (barrier){
        delete [] barrier;
    }
    barrier = 0;
    
    int j;
    if (botDim1 != 0){
        for (j=0; j < nProd; j++){
            // for each product
            delete [] botDim1[j];
            delete [] topDim1[j];
        }
        delete [] botDim1;
        delete [] topDim1;
    }

    botDim1 = 0;
    topDim1 = 0;
}

//----------------------------------------------------------------
/** add data for fd initialisation, timeline setup */
//temporary, need to chg
void FD1DRet::initSegments(const DateTimeArray& seg, 
                        const IntArray&   dens, 
                        DoubleArray& critSpacePts, //critical pts at space dimension
                        IntArray* isAddedSeg  //,  //                               
                        //bool needSpecialFD
                        )
{
    segDates = seg;
    density = dens;

    this->critSpacePts = critSpacePts;

   if (isAddedSeg)
        this->isAddedSeg = *isAddedSeg;

    //added->needSpecialFD = needSpecialFD;
}

void FD1DRet::setNbOfProd(int nProdIn){
    nProd = nProdIn + 2; // plus payoffIndex and payoffIndexOrig
}

void FD1DRet::setMaxNumOfValue(int maxNumOfValueIn){
    maxNumOfValue = maxNumOfValueIn;
}

void FD1DRet::setNumOfInsertNode(int numOfInsertNodeIn){
    numOfInsertNodeM = numOfInsertNodeIn;
}

int FD1DRet::getNumOfInsertNode(){
    return numOfInsertNodeM ;
}

/** Invoked after instrument has got its market data. */
void FD1DRet::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
    FDModel::getMarket(market, instruments);

    if( factors.size() != 1 )
        throw ModelException( "FD1DRet::getMarket", "only 1 factor suppported" );
}

/** accept or reject factors */
bool FD1DRet::acceptFactor( IMarketFactor * factor )
{
    return dynamic_cast< CAsset * >( factor ) != 0;
}

/** retrieving market data from MDF , new fd/tree interface*/    
void FD1DRet::retrieveFactor()
{
    static const string method = "FD1DRet::retrieveFactor";
    try
    {
        FDModel::retrieveFactor();

        underlying = CAssetConstSP::dynamicCast( factors[0] );
        if( ! underlying )
            throw ModelException( "FD1DRet::retrieveFactor", "only asset underlier suppported" );

        stepForward.resize(factors.size());  // factors.size() should be 1.
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** call timeline builder using initData */
/** can be used to prepare vars need to set the FD boundaries, such as total variance in LV*/
/** set the dates and FD set at time direction */

/**override the default FDModel::initModel(),*/
/**have more dada, in the case of DBlbarrier*/
/** to review*/
void FD1DRet::initModel(){
    static const string method = "FD1DRet::initModel";
    try{
        stepsPerYearFD = stepsPerYear; // set FDModel stepsPerYear

        int i;
        /**if more than 2 products we restrict:
         1: only has 1 barrier product and we assume it's the first product */

        if (critSpacePts.size()){

            arrangeDates();

            timeLine = TimeLineSP(new CTimeLine());
            //maybe need to call FDUtils::SetSegDates to set the criti seg, to review

            // create time line
            /* int effectiveStep = */ timeLine->CreateTimeLineSimple(segDates,
                                                                     timeMetric,
                                                                     stepsPerYear,
                                                                     critDates, 
                                                                     &isAddedSeg,
                                                                     &minStepInAddedSegM);

            //need a loop for nProd
            barrier = new FDUtilsBarrierSP[maxNumOfValue];

            for (i= 0; i < maxNumOfValue; i++ ){
                barrier[i] = FDUtilsBarrierSP(new FDUtilsBarrier());
            }
        }else{
            // set up time line
            FDModel::initModel();
        }

        stepForward[0].resize(timeLine->NumOfStep+1,0.0); //only has 1 underlying
        gridLevel1.resize(dim1);
        v_dxM.resize(dim1);
        v_dxMOrg.resize(dim1);

        // param validation
        if (truncation1D < 1.0 || truncation1D > 15.0) {
            throw ModelException(method, 
                                 "fd truncation must be between 1.0 and 15.0");
        }

        //to update
        if (numCoeffPDE !=4 ){
            throw ModelException("FD1DRet::validate", 
                "1 factor PDE only has 4 coefficients.");
        }

        //convert the registered string to enum type variables
        convertRegisteredString();

        if (whichChgVarDim1 == LOG_FWDX){
            needRecalcUndValue1 = true;
        }

        if((whichBarrierMethodM == VAR_GRID) && (needSpecialFD == true)){
            isVariableGrid = true;
        }

        //---------------
        //set space steps
        //---------------

        //create spot range
        range = TreeSliceEQ::Range::create( 0, dim1 - 1 );
    }
    catch(exception& e){
        throw ModelException(e, "FD1DRet::initModel()");
    }
}

//----------------------------------------------------------------

/**-----------------------------------------------------
    model initialisation  (second part)
    most memory init should be here
    finalise model initialisation once after initModel() and product init() 
    called by validate 
-----------------------------------------------------*/

/** set up the payoff */
void FD1DRet::finaliseModel(CControl*    control){
    static const string method = "FD1DRet::finaliseModel";
    try {
        int i, j;

        //to review
        //no barrier case, no extra data
        if (barrier ==0){
            needSpecialFD = false;

            const FDProductArray & products = getProducts();
            nProd = products.size();

            maxNumOfValue =0;
            for (j=0; j < nProd; j++){
                maxNumOfValue = Maths::max(maxNumOfValue, products[j]->getSlicesToDEV().size());
            }
        } 

        // these are sizes of solver
        // only reset when there is barriers

        bool rebuild = isRebuilt(control);
        if(rebuild){
            botDim1= new int*[nProd];
            topDim1= new int*[nProd];

            for (j=0; j<nProd; j++){
                // for each product
                botDim1[j]= new int [maxNumOfValue + 1];
                topDim1[j]= new int [maxNumOfValue + 1];

                for (i = 0; i <= maxNumOfValue; i++){
                    botDim1[j][i] = 0;
                    topDim1[j][i] = dim1-1;
                }
            }

            // set FD Grids. Use same grids to calc Greeks
            // to be reviewed !!!
            SetSpaceGridData();

            //calc the fd grids level
            computeFdGridLevel(0, dim1-1, v_dxM, lBound1, uBound1, gridLevel1);

        }
       
        // perform mapping for nodes. new interface products to old tree nodes
        ILatticeProdEDR * latticeProd = dynamic_cast< ILatticeProdEDR * >( prod.get() );
        if (latticeProd){
            int** dummyPriority = 0;
            latticeProd->mapInsertNode(insRangeM, insNodeM, insNodePriceM, dummyPriority);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

void FD1DRet::SetSpaceGridData(){

    double volForBound;

    setFdBounds(volForBound, truncation1D, lBound1, uBound1);

    //calculate dx
    // fill with constant for now, so only do once.
    double dxM = (uBound1 - lBound1) / (dim1-1);
    for (int i=0; i< (int)v_dxM.size(); i++){
        v_dxM[i] = dxM; //initial vector of dx
    }

    setVarSpaceSteps(volForBound);

    if (isVariableGrid == true){    
        //copy dxM to v_dxMOrg
        for (int i=0; i< (int)v_dxM.size(); i++){
            v_dxMOrg[i] = v_dxM[i]; //initial vector of dx
        }
    }
}

//----------------------------------------------------------------
/*
double FD1DRet::ConvertStoGridVar(int seg, double s){

    double x;
    
    switch (whichChgVarDim1) {
        case X:{
            //S
            x = s;
        }        
        break;
        case LOG_X: {
            //log(S)
            x = log(s);
        }
        break;
        case LOG_FWDX: {
            //log(fwd)                
            double mu = stepForward1[inSegEnd[seg]]/stepForward1[0];
            x = mu + log(s);
        }
        break;
    }

    return x;
}
*/
//----------------------------------------------------------------

void FD1DRet::setVarSpaceSteps(double volForBound){
        //only implemented when it's log(S), if not, need to modify addMoreSpacePts.

    switch (whichSetSpaceStepM) {
        case NONE:{

        }        
        break;
        case STD:{
            if (whichChgVarDim1 == LOG_X){
                isVariableGrid  = true;
                setVarSpaceStepsBasedonStd(volForBound, whichSetSpaceStepM, dim1-1);
            }else{
                throw ModelException("FD1DRet::setVarSpaceSteps", "only implemented for log(S)!.");
            }
        }
        break;
        case STD_BAR: { //case 1 + more pts arounds bars

            if ((whichChgVarDim1 == LOG_X) && (needSpecialFD == true)){
                isVariableGrid  = true;

                setVarSpaceStepsBasedonStd(volForBound, whichSetSpaceStepM, dim1-1);
            }else{
                throw ModelException("FD1DRet::setVarSpaceSteps", "only implemented for log(S)!.");
            }
        }
        break;
        default: {
            throw ModelException("FD1DRetLV::setVarSpaceSteps", "not implemented yet.");
        }
        break;
    }
}

//----------------------------------------------------------------

void FD1DRet::setVarSpaceStepsBasedonStd(double volForBound, TFdSetSpaceStepType whichSetSpaceStep, int n){

    int i;

    //hard code for now
    double trunc5 = truncation1D;
    double trunc3 = 3.0 /5.0 * truncation1D;
    double trunc2 = 2.0/5.0 * truncation1D;

    double lBound5, uBound5, lBound3, uBound3, lBound2, uBound2;

    setFdBounds(volForBound, trunc5, lBound5, uBound5);
    setFdBounds(volForBound, trunc3, lBound3, uBound3);
    setFdBounds(volForBound, trunc2, lBound2, uBound2);

    //weight for diff truncation 
    double w5 = 0.05;
    double w3 = 0.1;
    double w2 = 0.7;

    //2*w5 + 2 * w3 + w2 = 100%
    int n5 = int (w5 * n);
    int n3 = int (w3 * n);
    int n2 = n - 2 * n5 - 2 * n3;


    double dxTemp = (lBound3 - lBound5)/(n5+1);
    for( i = 1; i <= n5; i++){
        v_dxM[i] = dxTemp;
    }

    dxTemp = (lBound2 - lBound3)/(n3);
    for( i = n5 + 1; i <= n3 + n5; i++){
        v_dxM[i] = dxTemp;
    }

    dxTemp = (uBound2 - lBound2)/(n2);
    for( i = n3 + n5 + 1; i <= n2 + n3 + n5; i++){
        v_dxM[i] = dxTemp;
    }

    dxTemp = (uBound3 - uBound2)/(n3);
    for( i = n2 + n3 + n5 + 1; i <= n3 + n2 + n3 + n5; i++){
        v_dxM[i] = dxTemp;
    }

    dxTemp = (uBound5 - uBound3)/(n5);
    for( i = n3 + n2 + n3 + n5 + 1; i <= n5 + n3 + n2 + n3 + n5; i++){
        v_dxM[i] = dxTemp;
    }

    double oneDVolT = volForBound * sqrt(oneD);

    if (whichSetSpaceStep == STD_BAR){

        //test, to be generalized
        int minNeTemp = minNe;
        if (dim1 < 2.5 * minNe){//no more inserted pts around critical pts if dim is a small nb.
            minNeTemp = 0;
        }
        
        FDUtils::addMoreSpacePts(v_dxM, critSpacePts, n,  minNeTemp,
                    w5, w3, w2,
                    lBound5, uBound5, 
                    lBound3, uBound3, 
                    lBound2, uBound2, 
                    oneDVolT);
    }
}

//----------------------------------------------------------------
/**   fd stock level can call CDF mapping here */
//s : gridLevel

void FD1DRet::computeFdGridLevel (int bot, int top, vector<double>& vdx, 
            double lBoundary, double uBoundary, vector<double>& s){

    static const string method = "FD1DRet::computeFdGridLevel";
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

/**   fd stock level can call CDF mapping here */

void FD1DRet::computeUndLevel1(int bot, int top,  int step, double* s)const{

    //s = undValue
    static const string method = "FD1DRet::computeUndLevel1";
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
                double mu ;

                if (!Maths::isZero(stepForward[0][0])){  //stepForward[iAsset][step], iAsset = 1
                    mu = stepForward[0][step]/stepForward[0][0];
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

void FD1DRet::convertToGridBar(int step, int kValue){

    int k = kValue;

       switch (whichChgVarDim1) {
            case X:{
                //S
                if (barrier[k]->hasUpBarAtStep == true ){
                    barrier[k]->upBarrierGrid = barrier[k]->upBarrier;
                }

                if (barrier[k]->hasDownBarAtStep == true ){
                    barrier[k]->downBarrierGrid = barrier[k]->downBarrier;
                }
            }        
            break;
            case LOG_X: {
               //log(S)
                if (barrier[k]->hasUpBarAtStep == true ){
                    barrier[k]->upBarrierGrid = log(barrier[k]->upBarrier);
                }

                if (barrier[k]->hasDownBarAtStep == true ){
                    barrier[k]->downBarrierGrid = log(barrier[k]->downBarrier);
                }
            }
            break;
            case LOG_FWDX:{
                //log(fwd)

                double mu;

                if (!Maths::isZero(stepForward[0][0])){
                    mu = stepForward[0][step]/stepForward[0][0];
                }else{
                    throw ModelException("Forward of asset 1 is 0 at step = 0.");
                }

                if (barrier[k]->hasUpBarAtStep == true ){
                    barrier[k]->upBarrierGrid = log(barrier[k]->upBarrier) - log(mu);
                }

                if (barrier[k]->hasDownBarAtStep == true ){
                    barrier[k]->downBarrierGrid = log(barrier[k]->downBarrier) - log(mu);
                }
            }
            break;
        }
}

//----------------------------------------------------------------

bool FD1DRet::getDownBarrier(int iP, double s){
    
    barrier[iP]->downBarrier = s;
    
    return true;
} 

//----------------------------------------------------------------

bool FD1DRet::getUpBarrier(int iP, double s){

    barrier[iP]->upBarrier = s;
    
    return true;
}

//----------------------------------------------------------------

bool FD1DRet::hasBarrier(int step, int pStart, int pEnd){

    int k;

    bool hasBar = false;
    for (k = pStart; k<= pEnd; k++){
        if (barrier[k]->hasUpBarAtStep == true || barrier[k]->hasDownBarAtStep == true){
            hasBar = true;
            return hasBar ;
        }
    }
    return hasBar ;
}

//----------------------------------------------------------------
// get and return initial tree input steps
int FD1DRet::GetStepsPerYear() const
{
    return stepsPerYear;
}

//----------------------------------------------------------------

void FD1DRet::Spot::Product::update( int & step, FDProduct::UpdateType type )
{
    // we assume just need the first und level for spot here for now
    const TreeSlice & spot = getValue( step );
    int bot, top;
    spot.getCalcRange( bot, top );
    static_cast< FD1DRet * >( model )->computeUndLevel1( bot, top, step, spot.getValues() );
}

//----------------------------------------------------------------

void FD1DRet::Spot::Product::update( int & step, int & bot, int & top, FDProduct::UpdateType type )
{
    // we assume just need the first und level for spot here for now
    const TreeSlice & spot = getValue( step );
    spot.getCalcRange( bot, top );
    static_cast< FD1DRet * >( model )->computeUndLevel1( bot, top, step, spot.getValues() );
}

//----------------------------------------------------------------

FDProductSP FD1DRet::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            // avoid recursion (createProduct calls makeProduct back)
            if( ! dynamic_cast< const IndexSpecEDR * >( spec ) )
            {
                payoffIndex = FDModel::makeProduct( IProdCreatorSP( new Spot( spec->getFactor()->getName() ) ) );

                // need crateProduct (not makeProduct) to add product to the list
                payoffIndexOrig = FDModel::createProduct( IProdCreatorSP( new Spot( spec->getFactor()->getName() ) ) );

                return payoffIndex;
            }
        }
        else
        {
            throw ModelException( "FD1DRet::makeProduct",
                "IndexSpec source is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' is expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

FDModel::IFDSolverSP FD1DRet::createSolver()
{
    return IFDSolverSP(new FD1DRetSolver(this));
}

//----------------------------------------------------------------

/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double FD1DRet::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD1DRetLV::getPrice0";
    try {
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if( layer )
            return layer->getPrice0();

        double s0 = underlying->fwdValue(timeLine->StepDates[0]);
        const TreeSlice & s = payoffIndex->getValue( 0 );
        double * sVal = s.getValues();
        double * pVal = price.getValues();

        //int bot, top;
        //s->getCalcRange( bot, top );
        //int mid = Neighbour(s0, s, 0, top-bot,1);
        //price0[j]  = QuadraticInterp(s0, 
        //                        sVal[mid-1], 
        //                        sVal[mid], 
        //                        sVal[mid+1],
        //                        pVal[mid-1], 
        //                        pVal[mid], 
        //                        pVal[mid+1]);            

        double price0, deltaDummy, gammaDummy;
        int done = FDInterpolationD(dim1,sVal,pVal,1,&s0,&price0,&deltaDummy,&gammaDummy);

        //int done = FDInterpolationD(dim1,s[0],prices[j],1,&s0,&price0,&delta,&gamma);

        return price0;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

bool FD1DRetLoad(){
//    return (FD1DRet::TYPE && FD1DRet::IIntoProduct::TYPE);
    return (FD1DRet::TYPE != 0);
}

DRLIB_END_NAMESPACE
