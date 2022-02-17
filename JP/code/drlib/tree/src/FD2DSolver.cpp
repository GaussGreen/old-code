//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSolver.hpp
//
//   Description : two factor finite difference algorithm
//
//   Author      : Ning Shen
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD2DSolver.hpp"
#include "edginc/FDUtils.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE

/************ for debuging tree step outputs **************/
//#define FD2D_DEBUG_FILE "c:\\Temp\\Fd2f\\fd.txt"


/************ for debuging tree step outputs **************/
#define FD2D_DEBUG_FILE1 "c:\\temp\\fd.txt"

#ifdef FD2D_DEBUG_FILE1


void Debug_OutPutFD2DPrice(int step, double time, 
                      vector<double >& s1, vector<double >& s2, 
                      double*** price, 
                      int bot1, int top1, int bot2, int top2, 
                      int nPrice, 
                      const string& fileName, bool first){
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)

    for (int i=0; i<nPrice; i++){
        sprintf(fname, "%s%d.dat", fileName.c_str(), i);
        if (first)
            DumpFile = fopen(fname, "w");
        else
            DumpFile = fopen(fname, "a+");
        if (DumpFile)
        {
            sprintf(buffer, "#step = %d,  factor1,   price \n", step);
            fprintf(DumpFile, buffer);
            for (int ix =bot1; ix <=top1; ix ++)
            {
//                for (int iy =bot2; iy <=top2; iy ++)
//                {
//                    sprintf (buffer, "%f, %f, %f, %f # %d, %d \n", time, s1[ix], s2[iy], price[i][ix][iy], ix, iy);
//                    fprintf(DumpFile, buffer);
//                }

                sprintf (buffer, "%f, %f, %f, # %d \n", time, s1[ix],  price[i][ix][25], ix);
                fprintf(DumpFile, buffer);
            
            }
            sprintf (buffer, "\n");
            fprintf (DumpFile, buffer);
            fclose(DumpFile);
            DumpFile = 0;
        }
    }
}


void Debug_Output(int step, 
                      vector<double >& s1, 
                      double** price, 
                      int bot1, int top1,
                      int sDim,
                      int nPrice, 
                      const string& fileName, bool first){
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)

    for (int i=0; i<nPrice; i++){
        sprintf(fname, "%s%d.dat", fileName.c_str(), i);
        if (first)
            DumpFile = fopen(fname, "w");
        else
            DumpFile = fopen(fname, "a+");
        if (DumpFile)
        {
            sprintf(buffer, "#step = %d,  S1,   price \n", step);
            fprintf(DumpFile, buffer);
            for (int ix =bot1; ix <=top1; ix ++)
            {
//                for (int iy =bot2; iy <=top2; iy ++)
//                {
//                    sprintf (buffer, "%f, %f, %f, %f # %d, %d \n", time, s1[ix], s2[iy], price[i][ix][iy], ix, iy);
//                    fprintf(DumpFile, buffer);
//                }

                int s0Idx = ix % sDim; 

                sprintf (buffer, " %f, %f # %d, %d \n", s1[ix], price[i][s0Idx], s0Idx, ix);
                fprintf(DumpFile, buffer);
            
            }
            sprintf (buffer, "\n");
            fprintf (DumpFile, buffer);
            fclose(DumpFile);
            DumpFile = 0;
        }
    }
}

#endif

#ifdef FD2D_DEBUG_FILE




void Debug_OutPutFD2D(int step, double time, 
                      vector<double >& s1, vector<double >& s2, 
                      double*** price, 
                      int bot1, int top1, int bot2, int top2, 
                      int nPrice, 
                      const string& fileName, bool first){
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)

    for (int i=0; i<nPrice; i++){
        sprintf(fname, "%s%d.dat", fileName.c_str(), i);
        if (first)
            DumpFile = fopen(fname, "w");
        else
            DumpFile = fopen(fname, "a+");
        if (DumpFile)
        {
            sprintf(buffer, "#step = %d, TradeTime, factor1,  factor2, price \n", step);
            fprintf(DumpFile, buffer);
            for (int ix =bot1; ix <=top1; ix ++)
            {
//                for (int iy =bot2; iy <=top2; iy ++)
//                {
//                    sprintf (buffer, "%f, %f, %f, %f # %d, %d \n", time, s1[ix], s2[iy], price[i][ix][iy], ix, iy);
//                    fprintf(DumpFile, buffer);
//                }

                sprintf (buffer, "%f, %f, %f, %f # %d, %d \n", time, s1[ix], s2[ix], price[i][ix][10], ix, ix);
                fprintf(DumpFile, buffer);
            
            }
            sprintf (buffer, "\n");
            fprintf (DumpFile, buffer);
            fclose(DumpFile);
            DumpFile = 0;
        }
    }
}

void Debug_OutPutFD2DTime(int step, double time, double dt_trading, double dt_cal,
                      const string& fileName, bool first){
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)


    sprintf(fname, "%s%d.dat", fileName.c_str(), 4);
    if (first){
        DumpFile = fopen(fname, "w");

        if (DumpFile){
            sprintf(buffer, "#step = %d, TradeTime, , dt_trading, dt_cal,\n", step);
            fprintf(DumpFile, buffer);
        }
    }
    else{
        DumpFile = fopen(fname, "a+");
    }

    if (DumpFile){
//        sprintf(buffer, "#step = %d, TradeTime, , dt_trading, dt_cal,\n", step);
//        fprintf(DumpFile, buffer);

        sprintf (buffer, "%f, %f, %f,  # %d,  \n", time, dt_trading, dt_cal, step);
        fprintf(DumpFile, buffer);            
    }
    sprintf (buffer, "\n");
    fprintf (DumpFile, buffer);
    fclose(DumpFile);
    DumpFile = 0;        
}


#endif
/************ for debuging tree step outputs **************/
//-----------------------------------------------------------------------------

/** FD2DSolver algorithm class that supports FD2DSolver */
FD2DSolver::FD2DSolver(FD2D* model) : engine(model){

    int i;
    int k;

    //for now, theta is only used for ADI 1 without cross term
    theta = 0.5;

    //copy info. from model
    numCoeff = engine->numCoeffPDE;

    //1: ADI; 2: ADI advanced, ie: splitting
    solveMethod = engine->solveMethod;

    // Define x-axis as dim 1, and y-axis as dim 2;    
    xNum = engine->dim1;
    yNum = engine->dim2;
    //int tNum = xNum * yNum;

    v_dx.resize(xNum); 
    v_dy.resize(yNum); 
    for (k = 0; k < (int)(v_dx.size()); k++ ){
        v_dx[k] = engine->v_dxM[k];
    }

    for (k = 0; k < (int)(v_dy.size()); k++ ){
        v_dy[k] = engine->v_dyM[k];
    }

//    if ((solveMethod == ADI) && (engine->isVariableGrid)) {
//        throw ModelException("FD2DSolver::FD2DSolver", 
//            "variavle step isn't implemented for this scheme yet!.");
//    }

    //if (( !dynamic_cast< FD2DSVCJ*> (engine)->isPureHeston) && (engine->isVariableGrid)) {
    if (( (engine)->hasJumps) && (engine->isVariableGrid)) {
        throw ModelException("FD2DSolver::FD2DSolver", 
            "variavle step isn't implemented for this scheme yet!.");
    }

    coeff = 0;
    mSourceX = 0;
    mSourceY = 0;

    init_mem();

    a.resize(2);
    b.resize(2);
    c.resize(2);
    d.resize(2);
    e.resize(2);
    f.resize(2);
    for (k=0; k < 2; k++){
        a[k].resize(xNum);
        b[k].resize(xNum);
        c[k].resize(xNum);
        d[k].resize(xNum);
        e[k].resize(xNum);
        f[k].resize(xNum);

        for (i=0; i < xNum; i++){
            a[k][i].resize(yNum);
            b[k][i].resize(yNum);
            c[k][i].resize(yNum);
            d[k][i].resize(yNum);
            e[k][i].resize(yNum);
            f[k][i].resize(yNum);
        }
    }
    
//////////////////////////////

    alphaX.resize(xNum, 0.0);
    betaX.resize(xNum, 1.0);
    gammaX.resize(xNum, 0.0);
    aX.resize(xNum);
    bX.resize(xNum);
    cX.resize(xNum);

    alphaY.resize(yNum, 0.0);
    betaY.resize(yNum, 1.0);
    gammaY.resize(yNum, 0.0);
    aY.resize(yNum);
    bY.resize(yNum);
    cY.resize(yNum);
}

void FD2DSolver::init_mem(){

    int k, i, j, kk;

    coeff = new double**[numCoeff];

    // the size must be Dim1, Dim2 !
    for (i=0; i<numCoeff; i++){
        coeff[i] = new double* [xNum];
    
        for (j=0; j< xNum; j++){
            coeff[i][j] = new double [yNum];
        }               
    }


    // !!! temp solution, won't work with child product like this

    int tNum = xNum * yNum;

    int nbP = engine->maxNumOfValue;

    const FDProductArray & products = engine->getProducts();
    int nProd = products.size();
    
    mSourceX = new double** [nProd];
    mSourceY = new double** [nProd];

    // the size must be Dim1, Dim2 !
    for (kk=0; kk < nProd; kk++){
        mSourceX[kk] = new double* [nbP];
        mSourceY[kk] = new double* [nbP];

        for (k = 0; k < nbP; k++){
            mSourceX[kk][k] = new double [tNum];
            mSourceY[kk][k] = new double [tNum];
        }
    }

    //init value to 0
    for (kk=0; kk < nProd; kk++){
        for (k = 0; k < nbP; k++){
            for (j = 0; j < tNum; j++){
                mSourceX[kk][k][j] = 0.0;
                mSourceY[kk][k][j] = 0.0;
            }
        }
    }
}

//-----------------------------------------------------------------------------

FD2DSolver::~FD2DSolver(){
    int kk, k;

    if (coeff !=0){
        for (k = 0; k < numCoeff; k++){
            for (kk = 0; kk < xNum; kk++){
                delete [] coeff[k][kk];
            }
            delete [] coeff[k];
        }
        delete [] coeff;
    }
    coeff  =0;

    int maxNumOfValue = engine->maxNumOfValue;
    const FDProductArray & products = engine->getProducts();
    int nProd = products.size();
    
    if (mSourceX != 0){
        for (kk=0; kk < nProd; kk++){
            for (k = 0; k < maxNumOfValue; k++){
                delete [] mSourceX[kk][k] ;
                delete [] mSourceY[kk][k] ;
            }
            delete [] mSourceX[kk];
            delete [] mSourceY[kk];
        }
        delete [] mSourceX;
        delete [] mSourceY;
    }
    mSourceX = 0;
    mSourceY = 0 ;
}

//-----------------------------------------------------------------------------

// solve backward or forward induction through all steps
// return price
void FD2DSolver::roll(){
    static const string method = "FD2DSolver::Roll";

    // starting step, only backward for now, need update if adding forward
    int rollDirection = engine->isFwdInduction ? 1 : -1;
    int currStep = engine->isFwdInduction ? 0 : engine->timeLine->NumOfStep;
    
    // init segmenent index, not really used for now, to chg in the future if adds diff seg.
    //int currSeg = getfdSeg(currStep);

    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
        products[ prodIndex ]->preCalc( currStep );

    // update spot first is in this ft
    updateSolverInfo(currStep, 0/*engine->getSliceIndex(currStep)*/, rollDirection);

    tempSlices.resize( products.size() );
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( products[ prodIndex ]->isElementary() )
            continue;

        products[ prodIndex ]->update( currStep,
            engine->isFwdInduction ? FDProduct::FWD_0 : FDProduct::BWD_T );

        // initialize temporary slices
        const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();
        tempSlices[ prodIndex ].resize( slices.size() );
        for( int j = 0; j < (int)slices.size(); ++j )
            tempSlices[ prodIndex ][ j ] = slices[ j ]->clone( false );
    }

    // move a step to start sweeping
    currStep += rollDirection;

    // sweep the fd
    for (; currStep >= 0 && currStep <= engine->timeLine->NumOfStep; currStep += rollDirection)
    {
        //for barriers, need to update the info from products first
        //engine->prod->preCalcFD1D(currStep,  priceStart, priceEnd); 

        // for each product
        for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
            products[ prodIndex ]->preCalc( currStep );
            
        //pass new botDim, TopDim for x and y in case of barriers
        //call model to compute PDE coef
        updateSolverInfo(currStep, 0/*engine->getSliceIndex(currStep)*/, rollDirection);

        rollOneStep(currStep, rollDirection);

        // switching segment if needed
        switchSegment();
    }

#ifdef FD2D_DEBUG_FILE
    //fprintf(OUT_FILE, "%f6.3", nodePrice[idx][0][0][0]);
        //string fileName = FD2D_DEBUG_FILE;

        Debug_OutPutFD2D(currStep, engine->timeLine->TradeTime[0], 
                      engine->UndValue1, engine->UndValue2, 
                      nodePrice[currIdx], 
                      0, xNum-1, 0, yNum-1, 
                      engine->NumOfPrice, 
                      fileName, false);
#endif

}

//-----------------------------------------------------------------------------

//not sure it's useful for FD,
int FD2DSolver::getfdSeg(int step){

    int seg;

    if (step == 0){
        seg = 0;
    }
    else if (step == engine->timeLine->NumOfStep){
        seg = engine->timeLine->SegmentEnd.size()-1;
    }
    else{
        throw ModelException("FD1F::getfdSeg", "unable to locate fd segment.");
    }

    return seg;
}

//-----------------------------------------------------------------------------

//this ft is used to copy boundary and index info from model to solver at each segment
//for the momenet, we don't distinguish it and seg. to change in the future!!!!
void FD2DSolver::updateSolverInfo(int step, int idx, int rollDirection){

    if ((engine->isFwdInduction == false) ) {
        //calc the first slice
        if (step == engine->timeLine->NumOfStep ) {
            // update payoffIndex
            engine->payoffIndexUpdate(step, FDProduct::BWD_T);
        }

        //Underlying value will chg if we do log(fwd)
        //need to be before calc coeff, ex: LV    
        if (step < engine->timeLine->NumOfStep ) {//X=log(fwd), need to recal undValue                         
            if (step == engine->timeLine->NumOfStep -1 ) {
                //calc the second slice
                engine->payoffIndexUpdate(step, FDProduct::BWD);
            }else{
                //X=log(fwd), need to recal undValue 
                if ((engine->needRecalcUndValue1 == true)|| (engine->needRecalcUndValue2 == true)) {
                    // update payoffindex
                    engine->payoffIndexUpdate(step, FDProduct::BWD);
                }

//                if (engine->needRecalcUndValue2 == true){//X=log(fwd), need to recal undValue 
//                //engine->computeUndLevel2(0, yNum-1, dy, step);
//                //need to add if there is second state variable....        
//                }
            }            
        }
    }else{//fwd eq
        //calc the first slice
        if (step == 0 ) {
            // update spot first
            engine->payoffIndexUpdate(step, FDProduct::FWD_0);
        }

        //Underlying value will chg if we do log(fwd)
        //need to be before calc coeff, ex: LV    
        if (step > 0) {
                        
            if (step ==1){//calc the second slice
                engine->payoffIndexUpdate(step, FDProduct::FWD);
            }else{
                //X=log(fwd), need to recal undValue 
                if ((engine->needRecalcUndValue1 == true) || engine->needRecalcUndValue2 == true) {
                    // update spot 
                    engine->payoffIndexUpdate(step, FDProduct::FWD);               
                }
            }
        }
    }
            
    if ( (step < engine->timeLine->NumOfStep && engine->isFwdInduction == false) ||
        (step > 0 && engine->isFwdInduction == true)){

        engine->pdeCoeff(step, coeff, 0, xNum-1, 0, yNum-1);

        //keep 2 layers of coeffs in case we want to use it such as in LV case.
        passCoeff(step, idx, 0, xNum-1, 0, yNum-1);

        passCoeff(step, 1-idx, 0, xNum-1, 0,yNum-1);
    }
}

//-----------------------------------------------------------------------------
void FD2DSolver::passCoeff(int step, int idx, int low1, int top1, 
                           int low2, int top2){

// pass coeff to a,b,c,d,e,f)
//only to make the code easy to read
//for the final version, if we want,
//maybe, we can cut a,b,c,d,e,f
//but, need coefs at step and step + 1 anyway.
    int i = 0;
    int index_x;
    int index_y;

       for (index_x = low1; index_x <= top1; index_x++) {
           for (index_y = low2; index_y <= top2; index_y++) {
            a[idx][index_x][index_y] = coeff[i][index_x][index_y];
            b[idx][index_x][index_y] = coeff[i+1][index_x][index_y];
            c[idx][index_x][index_y] = coeff[i+2][index_x][index_y];
            d[idx][index_x][index_y] = coeff[i+3][index_x][index_y];
            e[idx][index_x][index_y] = coeff[i+4][index_x][index_y];
            f[idx][index_x][index_y] = coeff[i+5][index_x][index_y];
        }
    }
}

//-----------------------------------------------------------------------------

/** changing time line segment of a different density */
bool FD2DSolver::switchSegment(){
    
    bool changed = false;

    return changed;

}

//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*** forward or backward roll the fd at one time step 
*------------------------------------------------------------------------------*/

void FD2DSolver::rollOneStep(int step, int rollDirection)
{
    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // calc pv factor between steps
    double pv_x = 1.;
    if( ! engine->isFwdInduction ) // backward
    {
        pv_x = engine->discYC->pv(
            engine->timeLine->StepDates[step], engine->timeLine->StepDates[step+1] );
    }

    int idx = 0/*(engine->getSliceIndex(step)*/;

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( products[ prodIndex ]->isElementary() )
            continue;

        const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();
        const vector< TreeSliceSP > & slicesLast = tempSlices[ prodIndex ];

        // swap slices to preserve slicesLast
        for( int j = 0; j < (int)slices.size(); ++j )
        {
            static_cast< TreeSliceEQ & >( *slicesLast[ j ] ).swapValues(
                static_cast< TreeSliceEQ & >( *slices[ j ] ) );
        }

        int pStart = 0;
        int pEnd = slices.size() - 1;

#ifdef FD2D_DEBUG_FILE1
//    Debug_Output(step, 
//                      vector<double >& s1, 
//                      double** price, 
//                      int bot1, int top1,
//                      int sDim,
//                      int nPrice, 
//                      const string& fileName, bool first){
//
//        Debug_OutPut(currStep, engine->timeLine->TradeTime[0], 
//                      engine->UndValue1, engine->UndValue2, 
//                      nodePrice[currIdx], 
//                      0, xNum-1, 0, yNum-1, 
//                      engine->NumOfPrice, 
//                      fileName, false);
#endif

        nextStepFdBoundary(
            step,
            engine->botDim1[prodIndex],
            engine->topDim1[prodIndex],
            engine->botDim2[prodIndex],
            engine->topDim2[prodIndex],
            pStart, pEnd, pv_x, //engine->hasDiscTerminPDE ? pv_x : 1.,
            slices, slicesLast );

        euroOneStepwithSource2D(
            engine->maxNumOfValue,
            engine->botDim1[prodIndex],
            engine->topDim1[prodIndex],
            engine->botDim2[prodIndex],
            engine->topDim2[prodIndex],
            pStart, pEnd, idx,
            //(*products)[prodIndex]->engine->getSliceIndex(step-rollDirection),
            mSourceX[prodIndex], mSourceY[prodIndex],
            slices, slicesLast) ;

        // update the points at boundaries
        updateFdBoundary(
            step,
            engine->botDim1[prodIndex],
            engine->topDim1[prodIndex],
            engine->botDim2[prodIndex],
            engine->topDim2[prodIndex],
            pStart, pEnd,
            slices, slicesLast );

        products[ prodIndex ]->update( step,
            engine->isFwdInduction ? FDProduct::FWD : FDProduct::BWD );
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::nextStepFdBoundary(
    int step,
    int* low1, int* top1, int* low2, int* top2,
    int pStart, int pEnd, double pv_x,
    const vector< TreeSliceSP > & mcurrP,
    const vector< TreeSliceSP > & mlastP)
{
    int index_x, index_y;

    for( int k = pStart; k <= pEnd; ++k )
    {
        int l1 = low1[k];
        int t1 = top1[k];
        int l2 = low2[k];
        int t2 = top2[k];
        double * curr = (*mcurrP[k]).getValues();
        double * last = (*mlastP[k]).getValues();

        for (index_y = l2; index_y <= t2; index_y++){

            int i = index_y + l1 * yNum;            
            curr[i] = pv_x * last[i];
            
            i = index_y + t1 * yNum;            
            curr[i] = pv_x * last[i];
        }

        for (index_x = l1; index_x <= t1; index_x++){

            int i = l2 + index_x * yNum;
            curr[i] = pv_x * last[i];

            i = t2 + index_x * yNum;
            curr[i] = pv_x * last[i];
        }
    }
}

/**------------------------------------------------------------------------------
*   Description  :    called by Roll() after looping through all nodes.
*                    interpolates 2 extra node prices for fd roll over use.
*------------------------------------------------------------------------------*/

void FD2DSolver::updateFdBoundary(
    int step,
    int* low1, int* top1, int* low2, int* top2,
    int pStart, int pEnd,
    const vector< TreeSliceSP > & mcurrP,
    const vector< TreeSliceSP > & mlastP)
{
    int i;

    //the four corners
    //to review

    for( int k = pStart; k <= pEnd; ++k )
    {        
        int l1 = low1[k];
        int t1 = top1[k];
        int l2 = low2[k];
        int t2 = top2[k];
        double * curr = (*mcurrP[k]).getValues();
        double * last = (*mlastP[k]).getValues();

        i = l2 + l1 * yNum;
        curr[i] = last[i];

        i = t2 + l1 * yNum;
        curr[i] = last[i];

        i = l2 + t1 * yNum;
        curr[i] = last[i];

        i = t2 + t1 * yNum;
        curr[i] = last[i];

        if (engine->needSpecialFD == false){
            double ratioLow;
            double ratioTop;

            //direction X
            ratioLow = v_dx[l1+1] / v_dx[l1+2];
            ratioTop = v_dx[t1] / v_dx[t1-1];
            for (int index_y = l2 + 1; index_y < t2; index_y++){                
                i = index_y + l1 * yNum;
                curr[i] = (1.0 + ratioLow)*curr[i+yNum] - ratioLow* curr[i+2*yNum] ;
                curr[i] = Maths::max(curr[i]  ,0.0);

                i = index_y + t1 * yNum;
                curr[i] = (1.0 + ratioTop) *curr[i-yNum] - ratioTop * curr[i-2*yNum] ;
                curr[i] = Maths::max(curr[i]  ,0.0);
            }

            ratioLow = v_dy[l2+1] / v_dy[l2+2];
            ratioTop = v_dy[t2] / v_dy[t2-1];
            for (int index_x = l1 + 1; index_x < t1; index_x++){                
                i = l2 + index_x * yNum;
                curr[i] = (1.0 + ratioLow)*curr[i + 1] - ratioLow * curr[i + 2] ;

                curr[i] = Maths::max(curr[i]  ,0.0);

                i = t2 + index_x * yNum;
                curr[i] = (1.0 + ratioTop) *curr[i - 1] - ratioTop * curr[i - 2] ;

                curr[i] = Maths::max(curr[i]  ,0.0);
            }
        }else{//special barrier case
//                updateFdBoundaryBar(step, idx,
//                                    low1, top1, 
//                                    priceStart, priceEnd);
        }
    }

/*
    double pv; 
    //converting back to the true option value at current it.
    if ( engine->hasDiscTerminPDE == true){
        pv = 1;
    }else{
        pv = engine->pv(engine->timeLine->StepDates[step], 
                        engine->timeLine->StepDates[step+1]);
    }

    //think about index, for ex: with jumps
    //which is better: from 0 to xNum or from low to top???
    for( int k = pStart; k <= pEnd; ++k )
    {
        for (i = 0; i < xNum; i++){
            for (int j = 0; j < yNum; j++){
                int ij = j + i * yNum;
                curr[ij] = pv * curr[ij];
            }
        }        
    }
*/
}

//-----------------------------------------------------------------------------    
//classic ADI or LOD
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------    

void FD2DSolver::euroOneStepwithSource2D(
    int lastIndexTop,
    int* low1, int* top1, int* low2, int* top2,
    int pStart, int pEnd, int idx,
    double** mSrcX,
    double** mSrcY,
    const vector< TreeSliceSP > & mcurrP,
    const vector< TreeSliceSP > & mlastP)
{
    int index_x;
    int index_y;
    bool solveByLine = false;
    int l1, t1, l2, t2;

    //------------------------- 
    //compute the source terms
    //-------------------------
    for( int k = pStart; k <= pEnd; ++k )
    {
        l1 = low1[k];
        t1 = top1[k];
        l2 = low2[k];
        t2 = top2[k];

        calcSourceAll(solveMethod, l1, t1, l2, t2,idx, 
                        mSrcX[k], mSrcY[k], (*mlastP[k]).getValues());
    }

    //----------------
    //solvelution on X
    //----------------
    //loop on the biggest interval among the nb of value
    //should be the last element in the topDim, botDim, ie: max, or min
    int last = lastIndexTop;
    for (index_y = low2[last]; index_y <= top2[last]; index_y++ ) {

        //coeff of FD should be the same for loop k
        l1 = low1[last];
        t1 = top1[last];
        calcCoeffPdeAll_X(solveMethod, l1, t1, idx, index_y);      

        //loop on the nbOfPrice
        for( int k = pStart; k <= pEnd; ++k )
        {
            //loop only in the FD grid demanded by price[k]
            if ((index_y >= low2[k]) && (index_y <= top2[k])) {
                l1 = low1[k];
                t1 = top1[k];

                //calculate single price
                FDUtils::euroOneStepWithSource
                                    (l1, t1, 
                                    &alphaX[0], &betaX[0], &gammaX[0],
                                    &aX[0], &bX[0], &cX[0],
                                    (*mcurrP[k]).getValues(), (*mlastP[k]).getValues(), 
                                    mSrcX[k], solveByLine, index_y, yNum);  
            }
        }  
    }

    for( int k = pStart; k <= pEnd; ++k )
    {
        double * curr = (*mcurrP[k]).getValues();
        double * last = (*mlastP[k]).getValues();

        for (index_y = low2[k] ; index_y <= top2[k]; index_y++ ) {
            for (index_x = low1[k]; index_x <= top1[k]; index_x++ ) {
                int i = index_y + index_x * yNum;
                last[i] = curr[i];
            }
        }
    }    

    //----------------
    //solvelution on Y
    //----------------
    solveByLine = true;
    for (index_x = low1[last]; index_x <= top1[last]; index_x++ ) {

        //coeff of FD should be the same for loop k
        calcCoeffPdeAll_Y(solveMethod, low2[last], top2[last], idx,  index_x);

        //loop on the nbOfPrice
        for( int k = pStart; k <= pEnd; ++k )
        {

            //loop only in the FD grid demanded by price[k]
            if ((index_x >= low1[k]) && (index_x <= top1[k])) {

                l2 = low2[k];
                t2 = top2[k];
                //int yNumTemp = (t2 - l2) + 1;

                //calculate single price            
                FDUtils::euroOneStepWithSource
                                    (l2, t2, 
                                    &alphaY[0], &betaY[0], &gammaY[0],
                                    &aY[0], &bY[0], &cY[0],
                                    (*mcurrP[k]).getValues(), (*mlastP[k]).getValues(),
                                    mSrcY[k], solveByLine, index_x, yNum);
            }
        }  
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcSourceAll(const FD2D::TFdSolveType solveMethod,
                                    int low1, int top1,
                                    int low2, int top2,
                                    int idx,
                                    double* vSrcX, 
                                    double* vSrcY,
                                    double* M_price
                                    ){
    int i;
    int index_x, index_y;

    switch (solveMethod) {
        case FD2D::ADI:{
            calcSource(low1, top1, low2, top2,idx, vSrcX, M_price);

            //this is only for explicite method
            for (index_y = low2; index_y <= top2; index_y++ ) {
                for (index_x = low1; index_x <= top1; index_x++ ) {
                    i = index_y + index_x * yNum;
                    vSrcY[i] = 0;                        
                }
            }
        }        
        break;
        case FD2D::ADVANCE_ADI:{
            calcSourceAdvancedADI_X(low1, top1, low2, top2, 
                idx, vSrcX, M_price);
            calcSourceAdvancedADI_Y(low1, top1, low2, top2, 
                idx, vSrcY, M_price);      
        }
        break;
    }

    //adjusted source term if there is a jumps
    FD2DSVCJ* modelSVCJ = dynamic_cast< FD2DSVCJ*>(engine);
    //if (modelSVCJ->isPureHeston == false){

    if (engine->hasJumps == true){
        //put jumps explicitely when do ADI on X direction

        FD2DSolverJumpsSP fd2fsolverjumps(new FD2DSolverJumps(modelSVCJ, xNum, yNum));
        //FD2DSolverJumps fd2fsolverjumps(engine, xNum, yNum);

        if (!engine->isVariableGrid){
            //need constant dx and dy
            fd2fsolverjumps->compute_source_jumps(M_price, vSrcX, v_dx[1], v_dy[1]);
        }else{
            throw ModelException("FD2DSolver::calcSourceAll", 
                "variavle step isn't implemented for jumps solver yet!.");
        }
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeAll_X(const FD2D::TFdSolveType solveMethod, 
                         int low, int top, int idx,  int index_y){                                  

    switch (solveMethod) {
        case FD2D::ADI:{
            calcCoeffPdeADI_X(low, top, idx,  index_y);      
        }        
        break;
        case FD2D::ADVANCE_ADI:{
            calcCoeffPdeAdvancedADI_X(low, top, idx,  index_y);
        }
        break;
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeAll_Y(const FD2D::TFdSolveType solveMethod, 
                            int low, int top, int idx,  int index_x){                                  

    switch (solveMethod) {
        case FD2D::ADI:{
            calcCoeffPdeADI_Y(low, top, idx, index_x);      
        }        
        break;
        case FD2D::ADVANCE_ADI:{
            calcCoeffPdeAdvancedADI_Y(low, top, idx, index_x);
        }
        break;
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeADI_X(int low, int top, 
                                         int idx, int index_y ){
    int index_x;
    double Drift_Curr;
    double ItoTerm_Curr;
    double Drift_Last;
    double ItoTerm_Last;

    double opt1;
    double opt2;  

    int l1 = low +1;
    int t1 = top -1;
    int last_idx = 1- idx;

    if (engine->isVariableGrid){

        double dxi2;
        double dxiplus12;

        for (index_x = l1; index_x <= t1; index_x++) {
            dxi2 = v_dx[index_x] * v_dx[index_x];
            dxiplus12 = v_dx[index_x+1] * v_dx[index_x+1];

            opt1 = v_dx[index_x] * dxiplus12 + v_dx[index_x+1] * dxi2;
            //opt2 = dt / opt1;
            opt2 = 1 / opt1;

            Drift_Curr =  opt2 * a[idx][index_x][index_y];  
            ItoTerm_Curr = 2.0 * opt2 * c[idx][index_x][index_y]; 

            Drift_Last =  opt2 * a[last_idx][index_x][index_y];  
            ItoTerm_Last = 2.0 * opt2 * c[last_idx][index_x][index_y]; 

            alphaX[index_x] = theta * (Drift_Curr *  dxiplus12 - ItoTerm_Curr * v_dx[index_x+1]);
            aX[index_x]=(1.0 - theta) * ( - Drift_Last * dxiplus12 
                            + ItoTerm_Last * v_dx[index_x+1]);

            gammaX[index_x] = - theta * (Drift_Curr * dxi2 + ItoTerm_Curr * v_dx[index_x]);
            cX[index_x]= (1.0 - theta)* (Drift_Last * dxi2  + ItoTerm_Last * v_dx[index_x]);

    //        betaX[index_x] =1.0 + theta 
    //                        * ( Drift_Curr * (dxi2 - dxiplus12) +
    //                         ItoTerm_Curr * (v_dx[index_x] + v_dx[index_x+1] ) 
    //                        - dt * f[idx][index_x]);
    //        bX[index_x]= 1.0 - (1.0-theta) * ( Drift_Last * (dxi2 - dxiplus12) +
    //                        ItoTerm_Last *(v_dx[index_x] + v_dx[index_x+1])
    //                        - dt * f[last_idx][index_x]);             

            betaX[index_x] = 1 + theta 
                            * ( Drift_Curr * (dxi2 - dxiplus12) +
                             ItoTerm_Curr * (v_dx[index_x] + v_dx[index_x+1] ) 
                            -  f[idx][index_x][index_y]);
            bX[index_x]= 1 - (1.0-theta) * ( Drift_Last * (dxi2 - dxiplus12) +
                            ItoTerm_Last *(v_dx[index_x] + v_dx[index_x+1])
                            -  f[last_idx][index_x][index_y]);             
        }
    }else{
        //const grid

        //when const grid, use the first element
        double dx = v_dx[1];

    //    opt1 = 0.5 * dt / dx;
    //    opt2 = dt / (dx * dx);
        opt1 = 0.5 / dx;
        opt2 = 1.0 / (dx * dx);

        for (index_x = l1; index_x <= t1; index_x++) {
            Drift_Curr =  opt1 * a[idx][index_x][index_y];  
            ItoTerm_Curr =  opt2 * c[idx][index_x][index_y]; 

            Drift_Last =  opt1 * a[last_idx][index_x][index_y];  
            ItoTerm_Last =  opt2 * c[last_idx][index_x][index_y]; 

            alphaX[index_x] = theta * (Drift_Curr - ItoTerm_Curr);
            aX[index_x]=(1.0 - theta) * ( - Drift_Last 
                            + ItoTerm_Last);

            gammaX[index_x] = - theta * (Drift_Curr + ItoTerm_Curr);
            cX[index_x]= (1.0 - theta)* (Drift_Last + ItoTerm_Last);

            betaX[index_x] =1.0 + theta 
                            * ( 2.0 * ItoTerm_Curr - f[idx][index_x][index_y]);
            bX[index_x]= 1.0 - (1.0-theta) * (2.0 * ItoTerm_Last 
                            - f[last_idx][index_x][index_y]);             

    //        betaX[index_x] =1.0 + theta 
    //                        * ( 2.0 * ItoTerm_Curr - dt * f[idx][index_x][index_y]);
    //        bX[index_x]= 1.0 - (1.0-theta) * (2.0 * ItoTerm_Last 
    //                        - dt * f[last_idx][index_x][index_y]);             

        }
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeADI_Y(int low, int top, 
                                         int idx, int index_x){
    int index_y;
    double Drift_Curr;
    double ItoTerm_Curr;
    double Drift_Last;
    double ItoTerm_Last;

    double opt1;
    double opt2;  

    int l1 = low +1;
    int t1 = top -1;

    int last_idx = 1 - idx;  

    if (engine->isVariableGrid){

        double dyi2;
        double dyiplus12;

        for (index_y = l1; index_y <= t1; index_y++) {
            dyi2 = v_dy[index_y] * v_dy[index_y];
            dyiplus12 = v_dy[index_y+1] * v_dy[index_y+1];

            opt1 = v_dy[index_y] * dyiplus12 + v_dy[index_y+1] * dyi2;
            //opt2 = dt / opt1;
            opt2 = 1 / opt1;

            Drift_Curr =  opt2 * b[idx][index_x][index_y];  
            ItoTerm_Curr = 2.0 * opt2 * d[idx][index_x][index_y]; 

            Drift_Last =  opt2 * b[last_idx][index_x][index_y];  
            ItoTerm_Last = 2.0 * opt2 * d[last_idx][index_x][index_y]; 

            alphaY[index_y] = theta * (Drift_Curr *  dyiplus12 - ItoTerm_Curr * v_dy[index_y+1]);
            aY[index_y]=(1.0 - theta) * ( - Drift_Last * dyiplus12 
                            + ItoTerm_Last * v_dy[index_y+1]);

            gammaY[index_y] = - theta * (Drift_Curr * dyi2 + ItoTerm_Curr * v_dy[index_y]);
            cY[index_y]= (1.0 - theta)* (Drift_Last * dyi2  + ItoTerm_Last * v_dy[index_y]);

    //        betaX[index_x] =1.0 + theta 
    //                        * ( Drift_Curr * (dxi2 - dxiplus12) +
    //                         ItoTerm_Curr * (v_dx[index_x] + v_dx[index_x+1] ) 
    //                        - dt * f[idx][index_x]);
    //        bX[index_x]= 1.0 - (1.0-theta) * ( Drift_Last * (dxi2 - dxiplus12) +
    //                        ItoTerm_Last *(v_dx[index_x] + v_dx[index_x+1])
    //                        - dt * f[last_idx][index_x]);             

            betaY[index_y] = 1 + theta 
                            * ( Drift_Curr * (dyi2 - dyiplus12) +
                             ItoTerm_Curr * (v_dy[index_y] + v_dy[index_y+1] ) 
                            );
            bY[index_y]= 1 - (1.0-theta) * ( Drift_Last * (dyi2 - dyiplus12) +
                            ItoTerm_Last *(v_dy[index_y] + v_dy[index_y+1])
                            );             
        }
    }else{

        //when const grid, use the first element
        double dy = v_dy[1];

    //    opt1 = 0.5 * dt / dy;
    //    opt2 = dt / (dy * dy);

        opt1 = 0.5 / dy;
        opt2 = 1.0 / (dy * dy);

        for (index_y = l1; index_y <= t1; index_y++) {
            Drift_Curr =  opt1 * b[idx][index_x][index_y];  
            ItoTerm_Curr =  opt2 * d[idx][index_x][index_y]; 

            Drift_Last =  opt1 * b[last_idx][index_x][index_y];  
            ItoTerm_Last =  opt2 * d[last_idx][index_x][index_y]; 


            alphaY[index_y]= theta * (Drift_Curr - ItoTerm_Curr);
            aY[index_y] = 
                (1.0 - theta) * ( - Drift_Last + ItoTerm_Last);

            gammaY[index_y]= 
                    - theta * (Drift_Curr + ItoTerm_Curr);
            cY[index_y] =
                    (1.0 - theta) * (Drift_Last + ItoTerm_Last);


            betaY[index_y]= 
                    1.0 + theta * ( 2.0 * ItoTerm_Curr );
            bY [index_y] =
                    1.0 - (1.0-theta) * (2.0 * ItoTerm_Last);    
        }
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcSource(
                                int low1, int top1,
                                int low2, int top2,
                                int idx,
                                double*  vSrc, 
                                double* M_price
                                ){
    double opt;
    int index_x;
    int index_y;
    double source;
    int l1 = low1 +1;
    int t1 = top1 -1;

    int l2 = low2 +1;
    int t2 = top2 -1;


    if (engine->isVariableGrid){

        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {

                opt =  1.0 /( (v_dx[index_x] + v_dx[index_x+1] ) * (v_dy[index_y] + v_dy[index_y+1]) ) ; 

                int i1 = (index_y+1) + (index_x+1) * yNum;
                int i2 = (index_y+1) + (index_x-1) * yNum;
                int i3 = (index_y-1) + (index_x+1) * yNum;
                int i4 = (index_y-1) + (index_x-1) * yNum;

                source = M_price[i1] 
                        - M_price[i2]
                        - M_price[i3]
                        + M_price[i4];

                i1 = index_y + index_x * yNum;

                source *= opt * e[idx][index_x][index_y];  
                vSrc[i1] = source;
            }
        }
    }else{

        //when const grid, use the first element
        double dy = v_dy[1];
        double dx = v_dx[1];

    //    opt =  0.25 * dt/ (dx * dy) ; 
        opt =  0.25 / (dx * dy) ; 

        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {

                int i1 = (index_y+1) + (index_x+1) * yNum;
                int i2 = (index_y+1) + (index_x-1) * yNum;
                int i3 = (index_y-1) + (index_x+1) * yNum;
                int i4 = (index_y-1) + (index_x-1) * yNum;;

                source = M_price[i1] 
                        - M_price[i2]
                        - M_price[i3]
                        + M_price[i4];

                i1 = index_y + index_x * yNum;

                source *= opt * e[idx][index_x][index_y];  
                vSrc[i1] = source;
            }
        }
    }


/*
    for (index_y = l2; index_y <= t2; index_y++ ){
        for (index_x = l1; index_x <= t1; index_x++) {

            source = M_price[index_x + 1][index_y + 1] 
                    - M_price[index_x - 1][index_y + 1]
                    - M_price[index_x + 1][index_y - 1]
                    + M_price[index_x - 1][index_y - 1];

            source *= opt * e[idx][index_x][index_y];  
            vSrc[index_x][index_y] = source;
        }
    }

*/
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeAdvancedADI_X(int low, int top, int idx, int index_y){
    
    /** this scheme at direction X is the same as traditional LOD with theta =0.5 */
    theta = 0.5;
    calcCoeffPdeADI_X(low, top, idx, index_y );
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcCoeffPdeAdvancedADI_Y(int low, int top, int idx, 
                                                      int index_x){  
    int index_y;
    double Drift_Curr;
    double ItoTerm_Curr;

    double opt1;
    double opt2;  

    int l1 = low +1;
    int t1 = top -1;

    
    if (engine->isVariableGrid){
 
        for (index_y = l1; index_y <= t1; index_y++) {
            opt1 = v_dy[index_y] + v_dy[index_y + 1];                                  

//            Drift_Curr =  0.5 * dt / opt1 * b[idx][index_x][index_y];  
//            ItoTerm_Curr = dt / (opt1 * v_dy[index_y] * v_dy[index_y-1] ) 
//                            * d[idx][index_x][index_y]; 

            Drift_Curr =  0.5 / opt1 * b[idx][index_x][index_y];  
            ItoTerm_Curr = 1.0 / (opt1 * v_dy[index_y+1] * v_dy[index_y] ) 
                            * d[idx][index_x][index_y]; 


            alphaY[index_y]=  Drift_Curr - ItoTerm_Curr * v_dy[index_y+1] ;
            aY [index_y] = 0.0;

            gammaY[index_y] = - (Drift_Curr + ItoTerm_Curr * v_dy[index_y]);
            cY[index_y] = 0.0;
            
            betaY[index_y]=  1.0  +  opt1 * ItoTerm_Curr;
            bY[index_y] = 1.0;    
        }
    }else{
        //when const grid, use the first element
        double dy = v_dy[1];

//        opt1 = 0.5 * dt / dy;
//        opt1 *= 0.5; //1
//        opt2 = 0.5 * dt / (dy * dy);

        opt1 = 0.5  / dy;
        opt1 *= 0.5; //1
        opt2 = 0.5  / (dy * dy);

        
        for (index_y = l1; index_y <= t1; index_y++) {
            Drift_Curr =  opt1 * b[idx][index_x][index_y];  
            ItoTerm_Curr =  opt2 * d[idx][index_x][index_y]; 

            alphaY[index_y]=  Drift_Curr - ItoTerm_Curr;
            aY [index_y] = 0.0;

            gammaY[index_y] = - (Drift_Curr + ItoTerm_Curr);
            cY[index_y] = 0.0;
            
            betaY[index_y]=  1.0  +  2.0 * ItoTerm_Curr;
            bY[index_y] = 1.0;    
        }
    }
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcSourceAdvancedADI_X(int low1, int top1, int low2, int top2,  int idx, 
                                                   double* vSrc, double* M_price){

    double opt1;
    double opt2;
    double opt3;
    int index_x;
    int index_y;        
    double temp1;
    double temp2;
    double temp3;

    int l1 = low1 +1;
    int t1 = top1 -1;

    int l2 = low2 +1;
    int t2 = top2 -1;


    double dx = v_dx[1];

    if (engine->isVariableGrid){

        for (index_y = l2; index_y <= t2; index_y++ ){

            double opt0 = v_dy[index_y+1] + v_dy[index_y];

//            opt2 = 2.0 * dt/ (v_dy[index_y] * v_dy[index_y-1] * opt0);
//            opt1 = 0.5 * dt/ (dx * opt0) ; 
//            opt3 =  dt/ opt0;

            opt2 = 2.0 / (v_dy[index_y+1] * v_dy[index_y] * opt0);
//            opt1 = 0.5 / (dx * opt0) ; 
//            opt1 = 1.0 / ((v_dx[index_x] +v_dx[index_x+1]) * opt0) ; 
//
            opt3 =  1.0/ opt0;


            for (index_x = l1; index_x <= t1; index_x++) {

                opt1 = 1.0 / ((v_dx[index_x] +v_dx[index_x+1]) * opt0) ; 

                int i1 = (index_y+1) + index_x * yNum;
                int i2 = (index_y  ) + index_x * yNum;
                int i3 = (index_y-1) + index_x * yNum;

                temp1 = M_price[i1] * v_dy[index_y]
                            - M_price[i2] * opt0
                            + M_price[i3] * v_dy[index_y+1];             

                temp1 *= d[idx][index_x][index_y] * opt2;  

                int i4 = (index_y+1) + (index_x+1) * yNum;
                int i5 = (index_y-1) + (index_x+1) * yNum;
                int i6 = (index_y+1) + (index_x-1) * yNum;
                int i7 = (index_y-1) + (index_x-1) * yNum;

                temp2  = M_price[i4] 
                            - M_price[i5]
                            - M_price[i6]
                            + M_price[i7];

                temp2 *= e[idx][index_x][index_y] * opt1;
                temp3 =  M_price[i1] - M_price[i3];
                temp3 *=  b[idx][index_x][index_y] * opt3;
                vSrc[i2] =  temp1 + temp2 + temp3;
            }
        }
    }else{
        //when const grid, use the first element
        double dy = v_dy[1];

//        opt2 = dt/ (dy * dy) ; 
//        opt1 = 0.25 * dt/ (dx * dy) ; 
//        opt3 = 0.5 * dt/ dy;

        opt2 = 1.0/ (dy * dy) ; 
        opt1 = 0.25 / (dx * dy) ; 
        opt3 = 0.5 / dy;


        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {

                int i1 = (index_y+1) + index_x * yNum;
                int i2 = (index_y  ) + index_x * yNum;
                int i3 = (index_y-1) + index_x * yNum;


                temp1 = M_price[i1]
                        - 2.0 * M_price[i2] 
                        + M_price[i3];
                
                temp1 *= d[idx][index_x][index_y] * opt2;  

                int i4 = (index_y+1) + (index_x+1) * yNum;
                int i5 = (index_y-1) + (index_x+1) * yNum;
                int i6 = (index_y+1) + (index_x-1) * yNum;
                int i7 = (index_y-1) + (index_x-1) * yNum;

                temp2  = M_price[i4] 
                        - M_price[i5]
                        - M_price[i6]
                        + M_price[i7];

                temp2 *= e[idx][index_x][index_y] * opt1;

                temp3 =  M_price[i1] - M_price[i3];
                temp3 *=  b[idx][index_x][index_y] * opt3;
                vSrc[i2] =  temp1 + temp2 + temp3;
            }
        }
    }


/*
    if (engine->isVariableGrid){
        for (index_y = l2; index_y <= t2; index_y++ ){
            double opt0 = v_dy[index_y] + v_dy[index_y-1];
            opt2 = 2.0 * dt/ (v_dy[index_y] * v_dy[index_y-1] * opt0);

            opt1 = 0.5 * dt/ (dx * opt0) ; 
            opt3 =  dt/ opt0;

            for (index_x = l1; index_x <= t1; index_x++) {
                temp1 = M_price[index_x ][index_y+1] * v_dy[index_y-1]
                            - M_price[index_x][index_y] * opt0
                            + M_price[index_x ][index_y-1] * v_dy[index_y];             
                temp1 *= d[idx][index_x][index_y] * opt2;  
                temp2  = M_price[index_x + 1][index_y + 1] 
                            - M_price[index_x + 1][index_y - 1]
                            - M_price[index_x - 1][index_y + 1]
                            + M_price[index_x - 1][index_y - 1];

                temp2 *= e[idx][index_x][index_y] * opt1;
                temp3 =  M_price[index_x][index_y + 1] - M_price[index_x][index_y -1];
                temp3 *=  b[idx][index_x][index_y] * opt3;
                vSrc[index_x][index_y] =  temp1 + temp2 + temp3;
            }
        }
    }else{
        opt2 = dt/ (dy * dy) ; 
        opt1 = 0.25 * dt/ (dx * dy) ; 
        opt3 = 0.5 * dt/ dy;

        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {
                temp1 = M_price[index_x ][index_y+1]
                        - 2.0 * M_price[index_x][index_y] 
                        + M_price[index_x ][index_y-1];             
                temp1 *= d[idx][index_x][index_y] * opt2;  
                temp2  = M_price[index_x + 1][index_y + 1] 
                        - M_price[index_x + 1][index_y - 1]
                        - M_price[index_x - 1][index_y + 1]
                        + M_price[index_x - 1][index_y - 1];
                temp2 *= e[idx][index_x][index_y] * opt1;
                temp3 =  M_price[index_x][index_y + 1] - M_price[index_x][index_y -1];
                temp3 *=  b[idx][index_x][index_y] * opt3;
                vSrc[index_x][index_y] =  temp1 + temp2 + temp3;
            }
        }
    }
*/
}

//-----------------------------------------------------------------------------    

void FD2DSolver::calcSourceAdvancedADI_Y(
                            int low1, int top1, int low2, int top2,                                                    
                            int idx,
                            double* vSrc, double* M_price){

    int index_x;
    int index_y;

    double coeff1;
    double coeff2;
    double coeff3;
    double TermIto;
    double opt;
    double Drift_Curr;
    double opt1;

    int l1 = low1 +1;
    int t1 = top1 -1;
    int l2 = low2 +1;
    int t2 = top2 -1;

    if (engine->isVariableGrid){
        for (index_y = l2; index_y <= t2; index_y++ ){
            double opt0 = v_dy[index_y+1] + v_dy[index_y];

//            opt1 = 0.5 * dt / opt0;
//            opt = dt / (v_dy[index_y] * v_dy[index_y-1] *opt0);

            opt1 = 0.5 / opt0;
            opt = 1.0 / (v_dy[index_y+1] * v_dy[index_y] *opt0);


            for (index_x = l1; index_x <= t1; index_x++) {
                TermIto = opt * d[idx][index_x][index_y] ;
                Drift_Curr = opt1 * b[idx][index_x][index_y] ;

                coeff1 = - TermIto * v_dy[index_y] - Drift_Curr;
                coeff3 = - TermIto * v_dy[index_y+1] + Drift_Curr;
                coeff2 =  TermIto * opt0;


                int i1 = (index_y+1) + index_x * yNum;
                int i2 = (index_y  ) + index_x * yNum;
                int i3 = (index_y-1) + index_x * yNum;

                vSrc[i2] = coeff1 * M_price[i1] 
                                    + coeff2 * M_price[i2]
                                    + coeff3 * M_price[i3];
            }
        }
    }else{
        //when const grid, use the first element
        double dy = v_dy[1];

//        opt = 0.5 * dt / (dy * dy);
//        opt1 = 0.5 * dt /dy;
//        opt1 *= 0.5;

        opt = 0.5  / (dy * dy);
        opt1 = 0.5  /dy;
        opt1 *= 0.5;


        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {
     
                TermIto = opt * d[idx][index_x][index_y] ;
                Drift_Curr = opt1 * b[idx][index_x][index_y] ;

                coeff1 = - TermIto;
                coeff2 = - 2.0 * coeff1;
                coeff3 = coeff1 + Drift_Curr;
                coeff1 -= Drift_Curr;

                int i1 = (index_y+1) + index_x * yNum;
                int i2 = (index_y  ) + index_x * yNum;
                int i3 = (index_y-1) + index_x * yNum;

                vSrc[i2] = coeff1 * M_price[i1] 
                                    + coeff2 * M_price[i2]
                                    + coeff3 * M_price[i3];
            }                
        }
    }

/*
    if (engine->isVariableGrid){
        for (index_y = l2; index_y <= t2; index_y++ ){
            double opt0 = v_dy[index_y] + v_dy[index_y-1];
            opt1 = 0.5 * dt / opt0;
            opt = dt / (v_dy[index_y] * v_dy[index_y-1] *opt0);

            for (index_x = l1; index_x <= t1; index_x++) {
                TermIto = opt * d[idx][index_x][index_y] ;
                Drift_Curr = opt1 * b[idx][index_x][index_y] ;

                coeff1 = - TermIto * v_dy[index_y-1] - Drift_Curr;
                coeff3 = - TermIto * v_dy[index_y] + Drift_Curr;
                coeff2 =  TermIto * opt0;

                vSrc[index_x][index_y] = coeff1 * M_price[index_x][index_y + 1] 
                                    + coeff2 * M_price[index_x][index_y]
                                    + coeff3 * M_price[index_x][index_y - 1];
            }
        }
    }else{
        opt = 0.5 * dt / (dy * dy);
        opt1 = 0.5 * dt /dy;
        opt1 *= 0.5;

        for (index_y = l2; index_y <= t2; index_y++ ){
            for (index_x = l1; index_x <= t1; index_x++) {
     
                TermIto = opt * d[idx][index_x][index_y] ;
                Drift_Curr = opt1 * b[idx][index_x][index_y] ;

                coeff1 = - TermIto;
                coeff2 = - 2.0 * coeff1;
                coeff3 = coeff1 + Drift_Curr;
                coeff1 -= Drift_Curr;

                vSrc[index_x][index_y] = coeff1 * M_price[index_x][index_y + 1] 
                                    + coeff2 * M_price[index_x][index_y]
                                    + coeff3 * M_price[index_x][index_y - 1];
            }                
        }
    }
*/
}

DRLIB_END_NAMESPACE
