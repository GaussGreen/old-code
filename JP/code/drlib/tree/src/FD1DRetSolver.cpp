//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetSolver.cpp
//
//   Description : one factor finite difference algorithm
//
//   Author      : Xiaolan Zhang
//               : 
//   Date        : Apr, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1DRetSolver.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE
/************ for debuging tree step outputs **************/
//#define DEBUG_FILE1 "c:\\temp\\debugFD"

//#define FD1DRet_DEBUG_FILE1 "H:\\xlz\\temp\\fd1dOut.txt"

//#define DEBUG_FILE "c:\\temp\\debugFD"
#ifdef DEBUG_FILE

void Debug_OutPutFD(int step, double time, double **s, double **price, int bot, int top, int nPrice, const string& fileName, bool first)
{
    char fname[2048];
    char buffer[2048];

    static FILE *DumpFile = 0; // static file handle for dump file (keiji)

    for (int i=0; i<=nPrice; i++)
    {
        sprintf(fname, "%s%d.dat", fileName.c_str(), i);
        if (first)
            DumpFile = fopen(fname, "w");
        else
            DumpFile = fopen(fname, "a+");
        if (DumpFile)
        {
            sprintf(buffer, "#step = %d, %f, TradeTime, stock , Price \n", step, engine->timeLine->StepDates[step]);
            fprintf(DumpFile, buffer);
            for (int j_debug=-bot; j_debug<=top; j_debug ++)
            {
                sprintf (buffer, "%f, %f, %f # %d \n", time, s[i][j_debug], price[i][j_debug], j_debug);
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


#ifdef FD1DRet_DEBUG_FILE1
#include "stdio.h"  

void Debug_Output1D(DateTime&  TradeDate,             
                    int step, 
                      double* s1, 
                      double** lastPrice, 
                      double** price, 
                      int bot1, 
                      int top1,
                      int sDim,
                      int nPrice, 
                      const string& fileName, 
                      bool first){
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
            sprintf(buffer, "#StepDate , \n",  TradeDate);
            sprintf(buffer, "#step = %d, %d,  S1,   lastPrice, price, s0Idx, ix \n", step, TradeDate);
            fprintf(DumpFile, buffer);
            for (int ix =bot1; ix <=top1; ix ++)
            {
                //int s0Idx = ix % sDim; 
                int s0Idx = ix;
                sprintf (buffer, " %f, %f, %f # %d, %d \n", s1[s0Idx], lastPrice[i][s0Idx], price[i][s0Idx], s0Idx, ix);

//                sprintf (buffer, " %f, %f, %f # %d, %d \n", s1[ix], lastPrice[i][ix], price[i][ix], s0Idx, ix);
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
//-----------------------------------------------------------------------------

/** FD1DRetSolver algorithm class that supports FD1DRetSolver */
FD1DRetSolver::FD1DRetSolver(FD1DRet* engine) : engine(engine){

    int k;

    theta = 0.5;

    //copy info. from model
    numCoeff = engine->numCoeffPDE;
    solveMethod = engine->solveMethod;
    xNum = engine->dim1;  //vector form 0, to dim1
    whichBarrierMethod = engine->whichBarrierMethodM;
    
    v_dx.resize(xNum );
    for (k = 0; k < (int)(v_dx.size()); k++ ){
        v_dx[k] = engine->v_dxM[k];
    }

    init_mem();

    a.resize(2);
    c.resize(2);
    f.resize(2);
    g.resize(2);

    for (k=0; k < 2; k++){
        a[k].resize(xNum);
        c[k].resize(xNum);
        f[k].resize(xNum);
        g[k].resize(xNum);
    }
    
    alphaX.resize(xNum);
    betaX.resize(xNum);
    gammaX.resize(xNum);
    aX.resize(xNum);
    bX.resize(xNum);
    cX.resize(xNum);

    alphaX[0] = 0.0;
    betaX[0] = 1.0;
    gammaX[0] = 0.0;

    alphaX[xNum-1] = 0.0;
    betaX[xNum-1] = 1.0;
    gammaX[xNum-1] = 0.0;

    int maxBits = 1 << engine->range->nDim;
    sliceCache.resize( maxBits );
    for( int b = 0; b < maxBits; ++b )
        sliceCache[ b ] = DYNAMIC_POINTER_CAST< TreeSliceEQ >( engine->createSlice( b ) );

    //FD1DRetSolverJumpsSP FD1DRetSolverjumps(new FD1DRetSolverJumps(engine, xNum+1, yNum+1));   
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::init_mem(){

    int k, kk, j;

    coeff = 0;    
    mSourceX = 0;

    coeff = new double*[numCoeff];

    // the size must be Dim1, Dim2 !
    for (k =0; k<numCoeff; k++){
        coeff[k] = new double [xNum];    
    }

    // !!! temp solution, won't work with child product like this
    int nProd = engine->nProd;
    int maxNumOfValue =engine->maxNumOfValue;

    mSourceX = new double** [nProd];
    for (kk=0; kk < nProd; kk++){
        mSourceX[kk] = new double* [maxNumOfValue];

        for (k = 0; k < maxNumOfValue; k++){
            mSourceX[kk][k] = new double [xNum];
        }
    }

    //init value to 0
    for (kk=0; kk < nProd; kk++){
        for (k = 0; k < maxNumOfValue; k++){
            for (j = 0; j < xNum; j++){
                mSourceX[kk][k][j] = 0.0;
            }
        }
    }

    /**----------
        barriers
    ----------*/

    //0: low bar, 1: high bar
    numOfInsertNode = engine->numOfInsertNodeM; //2*nb of value

    insNode_FdIndex.resize(numOfInsertNode);
    insNodeActive.resize(numOfInsertNode, false);
    insNodePriceBak.resize(maxNumOfValue);

    for ( k= 0; k < maxNumOfValue; k++ ){
        insNodePriceBak[k].resize(numOfInsertNode);
        for (j = 0; j < numOfInsertNode; j++){
            insNodePriceBak[k][j] = 0;
        }
    }

    //hard code for now
    //need to add prod dimension
    insNode = engine->insNodeM;
    insNodePrice = engine->insNodePriceM;

    //to avoid UMR
    for (j = 0; j < numOfInsertNode; j++){
        (*insNode)[0][j] = 1.0;
    }
    
    for ( k= 0; k < maxNumOfValue; k++ ){
        for (j = 0; j < numOfInsertNode; j++){
            (*insNodePrice)[0][k][j] = 0;
        }
    }

    //FD1DRetSolverjumps = new FD1DRetSolverJumps(engine, xNum, yNum);
}

//-----------------------------------------------------------------------------

FD1DRetSolver::~FD1DRetSolver(){
    int k, kk;

    if (coeff !=0){
        for (k = 0; k < numCoeff; k++){
            delete [] coeff[k];
        }
        delete [] coeff;
    }
    coeff  =0;

    int nProd = engine->nProd;
    int nbP = engine->maxNumOfValue;

    if (mSourceX != 0){
        for (kk=0; kk < nProd; kk++){
            for (k = 0; k < nbP; k++){
                delete [] mSourceX[kk][k] ;
            }
            delete [] mSourceX[kk];
        }
        delete [] mSourceX;
    }
    mSourceX = 0;

    //delete [] FD1DRetSolverjumps;
}

//-----------------------------------------------------------------------------

// solve backward or forward induction through all steps
// return price
void FD1DRetSolver::roll(){

    static const string method = "FD1DRetSolver::Roll";

    // starting step
    int rollDirection = engine->isFwdInduction ? 1 : -1;
    int step = engine->isFwdInduction ? 0 : engine->timeLine->NumOfStep;
    
    // init segmenent index, 
    int  currSeg = getfdSeg(step);

    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
        products[ prodIndex ]->preCalc( step );

    // update spot first
    updateSolverInfo(step, 0/*engine->getSliceIndex(step)*/, rollDirection);

    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( ! products[ prodIndex ]->isElementary() )
        {
            products[ prodIndex ]->update( step,
                engine->isFwdInduction ? FDProduct::FWD_0 : FDProduct::BWD_T );
        }

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( engine->DEBUG_DumpToFile.length() )
        {
            string fileName =
                engine->DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "wt" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i ", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }
    }

#ifdef DEBUG_FILE
        string fileName = DEBUG_FILE;
        int bot, top, pStart, pEnd;

        //double** undValue1;
        double** undValue1Org;
        if (engine->needSpecialFD == false)
            updatePayoffIndex(step, FD1DRetSpotSP::dynamicCast(engine->payoffIndexOrig));
        engine->payoffIndexOrig->getValue(step);

        double** price;
        engine->prod->getValue(step);

        Debug_OutPutFD(step, engine->timeLine->TradeTime[step], undValue1Org, price, bot, top, pEnd-pStart, fileName, false);     
#endif
        
    // move a step to start sweeping
    step += rollDirection;

    // sweep the fd
    for (; step >=0 && step <= engine->timeLine->NumOfStep; step += rollDirection)
    {
        // for each product
        for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
            products[ prodIndex ]->preCalc( step );

        // update spot first
        updateSolverInfo(step, 0/*engine->getSliceIndex(step)*/, rollDirection);

        rollOneStep(step, rollDirection);

#ifdef DEBUG_FILE
        string fileName = DEBUG_FILE;
        int bot, top, pStart, pEnd;

        //double** undValue1;
        double** undValue1Org;
        if (engine->needSpecialFD == false)
            updatePayoffIndex(step, FD1DRetSpotSP::dynamicCast(engine->payoffIndexOrig));
        engine->payoffIndexOrig->getValue(step);

        double** price;
        engine->prod->getValue(step);

        Debug_OutPutFD(step, engine->timeLine->TradeTime[step], undValue1Org, price, bot, top, pEnd-pStart, fileName, false);     
#endif

        // switching segment if needed
        //has segment, but didn't chg space grid, so, shouldn't need interpol here
        switchSegment(step, 0/*engine->getSliceIndex(step)*/, currSeg, rollDirection);
    }
}

//-----------------------------------------------------------------------------

int FD1DRetSolver::getfdSeg(int step){
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
void FD1DRetSolver::updateSolverInfo(int step, int idx, int rollDirection)
{
    updatePayoffIndex(step, engine->payoffIndex);

    if (engine->needSpecialFD == true){
        updatePayoffIndex(step, engine->payoffIndexOrig);
        //set the new value for botDim, topDim
        //will calc the cloest index to barrier if isVariableGrid = true,
        //check if we have barrier at this step 
        setFDSolverDims(step);
    }

    /**----------
    // calc coeff
    //----------*/

    int pStart, pEnd;
    //at step, we use the same coeffs for coeff{idx] and coeff[1 - idx]
    if ( (step < engine->timeLine->NumOfStep && engine->isFwdInduction == false)  ||
        (step > 0 && engine->isFwdInduction == true) ){

        //set inserted node, reset v_dx, recalc underlying value for inserted node cases
        if (engine->needSpecialFD == true)
        {
            int prodIndex;
            const FDProductArray & products = engine->getProducts();

            // for each product
            for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
            {
                if( products[ prodIndex ]->isElementary() )
                    continue;

                const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();

                pStart = 0;
                pEnd = slices.size() - 1;

                setInsNodeInfo(pStart, pEnd);
                for( int k = pStart; k <= pEnd; ++k )
                {
                    TreeSliceLayer* layer = dynamic_cast< TreeSliceLayer* >( slices[ k ].get() );
                    int l = layer ? 0 : k;
                    int h = layer ? layer->getSlices().size() - 1 : k;
                    const vector< TreeSliceSP > & sCurr = layer ? layer->getSlices() : slices;

                    for( int j = l; j <= h; ++j )
                    {
                        TreeSliceEQ & curr = static_cast< TreeSliceEQ & >( *sCurr[ j ] );

                        calcBCDsBefore(
                            step,
                            engine->botDim1[prodIndex],
                            engine->topDim1[prodIndex],
                            pStart, pEnd, k,
                            curr );
                    }
                }

                //get approximated BCDs from the product at inserted nodes   
                engine->prod->update(step, FDProduct::BWD_NODE_INSERTION);

                for( int k = pStart; k <= pEnd; ++k )
                {
                    int pEndTemp = pEnd;
                    if (engine->DEBUG_NODE_GEN==false){
                        pEndTemp = pStart;
                    }
                    
                    for (int kk = pStart; kk <= pEndTemp; kk++ ){
                        int i = 2 * (kk-pStart);
                        if ((engine->barrier[kk]->hasDownBarAtStep == true) ){
                            engine->barrier[k]->valueAtDownBar = (*insNodePrice)[0][k][i];
                        }            

                        if ((engine->barrier[kk]->hasUpBarAtStep == true) ){
                            engine->barrier[k]->valueAtUpBar = (*insNodePrice)[0][k][i + 1];
                        }
                    }
                }

                //maybe need to put
                if (whichBarrierMethod == FD1DRet::VAR_GRID){
                    resetdx(v_dx, prodIndex, pStart, pEnd); //chg v_dx

                    //recalc undValue
                    for (int j=0; j < numOfInsertNode; j++ ){
                        if (insNodeActive[j] == true){ //has a bar here
                            int index = insNode_FdIndex[j];

                            if (engine->DEBUG_NODE_GEN == false){
                                if (j ==0) {//down bar
                                    index = insNode_FdIndex[j];
                                }else{
                                    index = insNode_FdIndex[j] + 1;
                                    index = Maths::min(index, xNum-1); //the size of payoffIndex is xNum
                                }
                            }
                            //recalc payoffIndex at inserted node via v_dx
                            static_cast< FD1DRet::Spot::Product & >(*engine->payoffIndex).update(step, index, index, FDProduct::BWD);  
                        }
                    }
                }

            }
            //the above code should be moved here if prod dim is added (use same fd grid for all prod)
        }

        engine->pdeCoeff(step, coeff, 0, xNum-1);
        passCoeff(idx, 0, xNum-1);
        //keep it in case we wants 2 layers coeff
        passCoeff(1-idx, 0, xNum-1);
    }

    if ((engine->isFwdInduction == false) && (step < engine->timeLine->NumOfStep )) {
        //set FD scheme
        if (step == engine->timeLine->NumOfStep -1 && engine->isFwdInduction == false){
            theta = 1;
        }else{
            theta = 0.5;
        }

        //if special barriers
        if (engine->needSpecialFD == true){
            if (whichBarrierMethod == FD1DRet::VAR_GRID){
                if (engine->hasBarrier(step, pStart, pEnd)){
                    theta = 1;                
                }
                else{
                    theta = 0.5;
                }
            }
        }    
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::updatePayoffIndex(int step, FDProductSP payoffIndex){

    if ((engine->isFwdInduction == false) ) {
        //update spot at first slice
        if (step == engine->timeLine->NumOfStep ){
            payoffIndex->update(step, FDProduct::BWD_T);
        }

        //Underlying value will chg if we do log(fwd) or var grid
        //need to be before calc coeff, ex: LV    
        if ((step == engine->timeLine->NumOfStep -1)  //calc the second slice
            || (step < engine->timeLine->NumOfStep - 1 && 
                engine->needRecalcUndValue1 == true) ){//X=log(fwd), need to recal undValue 
            payoffIndex->update(step, FDProduct::BWD);
        }
    }
    else{
        //update spot at first slices
        if (step == 0 ){
            payoffIndex->update(step, FDProduct::FWD_0);
        }

        //Underlying value will chg if we do log(fwd) or var grid
        //need to be before calc coeff, ex: LV    
        if ( (step == 1)   // calc the second slice
            || (step > 1 && 
                engine->needRecalcUndValue1 == true) ){//X=log(fwd), need to recal undValue 
            payoffIndex->update(step, FDProduct::FWD);
        }
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::postUpdateSolverInfo(int step){

    //come back to org grid
    if (numOfInsertNode && whichBarrierMethod == FD1DRet::VAR_GRID){
        const TreeSliceEQ & undValue1 =
            static_cast< const TreeSliceEQ & >( engine->payoffIndex->getValue( step ) );
        const TreeSliceEQ & undValue1Org =
            static_cast< const TreeSliceEQ & >( engine->payoffIndexOrig->getValue( step ) );

        //copy back recalc node in  undValue
        for (int j=0; j < numOfInsertNode; j++ ){
            if (insNodeActive[j] == true){
                int index = insNode_FdIndex[j];                
                if(engine->DEBUG_NODE_GEN == false){
                    if (j ==0) {
                        index = insNode_FdIndex[j];
                    }else{
                        index = insNode_FdIndex[j] + 1;
                        index = Maths::min(index, xNum-1); //the size of payoffIndex is xNum
                    }
                }

                undValue1[index] = undValue1Org[index];

                //we assume that we always use the same coeffs as at step
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::passCoeff(int idx, int bot1, int top1){

// pass coeff to a,b,c,d,e,f)
//only to make the code easy to read
//for the final version, if we want to speed,
//maybe, we can cut a,b,c,d,e,f
    int i = 0;
    int index_x;

       for (index_x = bot1; index_x <= top1; index_x++) {
        a[idx][index_x] = coeff[i][index_x];
        c[idx][index_x]= coeff[i+1][index_x];
        f[idx][index_x]= coeff[i+2][index_x];
        g[idx][index_x]= coeff[i+3][index_x];
    }
}

//-----------------------------------------------------------------------------

/** changing time line segment of a different density */
/** this function isn't really used*/
bool FD1DRetSolver::switchSegment(int step, int idx, int currSeg, int RollDirection){

    bool changed = false;
    return changed;
}

//------------------------------------------------------------------------------
/*** forward or backward roll the fd at one time step 
*------------------------------------------------------------------------------*/

//#define DEBUG
//#include <dprintf.h>
//#undef min
//#undef max

void FD1DRetSolver::rollOneStep(int step, int rollDirection){
    
    //if we want try diff method to solve, put it here    
    int pStart, pEnd;

    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // calc pv factor between steps
    double pv_x = 1.;
    if( ! engine->isFwdInduction ) // backward
    {
        pv_x = engine->discYC->pv(
            engine->timeLine->StepDates[step], engine->timeLine->StepDates[step+1] );
    }

    int idx = 0/*engine->getSliceIndex(step)*/;

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( products[ prodIndex ]->isElementary() )
            continue;

        const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();

        pStart = 0;
        pEnd = slices.size() - 1;

#ifdef FD1DRet_DEBUG_FILE1
        string fileName = FD1DRet_DEBUG_FILE1;

        double** undValue1;
        int nPrice = 1;
        int dummy, bot, top;

        engine->payoffIndex->getValue(step);

        Debug_Output1D(engine->timeLine->StepDates[step],
                        step, 
                      undValue1[0], 
                      lastPrice,
                      price, 
                      bot, 
                      top,
                      1,
                      nPrice, 
                      fileName, 
                      false);
#endif

        //----------------
        //solvelution on X
        //----------------
        //coeff of FD should be the same for loop k
        //always calc the whole coeff. due to 2 layer in FD
        //calcCoeffPdeAll_X(solveMethod, 0, xNum-1, idx);   
        int lastIndexTop = engine->maxNumOfValue;
        calcCoeffPdeAll_X(
            solveMethod,
            engine->botDim1[prodIndex][lastIndexTop],
            engine->topDim1[prodIndex][lastIndexTop],
            idx );

        for( int k = pStart; k <= pEnd; ++k )
        {
            TreeSliceLayer * layer = dynamic_cast< TreeSliceLayer * >( slices[ k ].get() );
            int l = layer ? 0 : k;
            int h = layer ? layer->getSlices().size() - 1 : k;
            const vector< TreeSliceSP > & sliceList = layer ? layer->getSlices() : slices;

            for( int j = l; j <= h; ++j )
            {
                TreeSliceEQ & currSlice = static_cast< TreeSliceEQ & >( *sliceList[ j ] );
                int dimBits = currSlice.dimBits();
                TreeSliceEQ & prevSlice = *sliceCache[ dimBits ];

                // preserve last step slice in prevSlice
                prevSlice.swapValues( currSlice );

                nextStepFdBoundary(
                    step,
                    engine->botDim1[prodIndex],
                    engine->topDim1[prodIndex],
                    k,
                    engine->hasDiscTerminPDE ? pv_x : 1.,
                    currSlice, prevSlice );

                euroOneStepwithSource1D(
                    engine->maxNumOfValue,
                    engine->botDim1[prodIndex],
                    engine->topDim1[prodIndex],
                    k, idx,
                    mSourceX[prodIndex],
                    currSlice, prevSlice );

                // update the points at boundaries
                updateFdBoundary(
                    step,
                    engine->botDim1[prodIndex],
                    engine->topDim1[prodIndex],
                    pStart, pEnd,
                    k, engine->hasDiscTerminPDE ? 1. : pv_x,
                    currSlice, prevSlice );

                //if true, set no negative value to prices
                forceNonNegative( k, currSlice );
            }
        }

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( engine->DEBUG_DumpToFile.length() && ! products[ prodIndex ]->isElementary() )
        {
            string fileName =
                engine->DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "at" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i+", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }

        if( products[ prodIndex ]->isElementary() )
            continue;

        products[ prodIndex ]->update( step, engine->isFwdInduction ? FDProduct::FWD : FDProduct::BWD );

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( engine->DEBUG_DumpToFile.length() && ! products[ prodIndex ]->isElementary() )
        {
            string fileName =
                engine->DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "at" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i-", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }
    }

    postUpdateSolverInfo(step);
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::forceNonNegative( int k, const TreeSliceEQ & price )
{
    if( engine->DEBUG_forceNonNegative )
    {
        for( int i = 0; i < xNum; ++ i)
            price[i] = Maths::max(0.0, price[i]);
    }
}
 
//-----------------------------------------------------------------------------

void FD1DRetSolver::nextStepFdBoundary(
    int step,
    int* bot1, int* top1,
    int k, double pv_x,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP )
{
    if( ! engine->needSpecialFD )
    {
        int l1 = bot1[k];
        int t1 = top1[k];
        mcurrP[l1] = pv_x * mlastP[l1];
        mcurrP[t1] =  pv_x * mlastP[t1];
    }
    else
        nextStepFdBoundaryBar(step,  bot1, top1, k, pv_x, mcurrP, mlastP);
}

/**------------------------------------------------------------------------------
*   Description  :    called by rollOneStep() after solving the equation
*------------------------------------------------------------------------------*/
void FD1DRetSolver::updateFdBoundary(
    int step,
    int* bot1, int* top1,
    int pStart, int pEnd,
    int k, double pv_x,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP )
{
    if( ! engine->needSpecialFD )
    {
        interpBoundaryLow( bot1[k], mcurrP );
        interpBoundaryUp( top1[k], mcurrP );
    }
    else // special barrier case
    {
        updateFdBoundaryBar(
            step,
            bot1, top1,
            pStart, pEnd,
            k, pv_x,
            mcurrP, mlastP );
    }

    if( ! engine->isFwdInduction )
    {
        for( int i = 0; i < xNum; ++i )
            mcurrP[i] = pv_x * mcurrP[i];
    }
}

//-----------------------------------------------------------------------------    
//---------------------------common fd scheme fts for all prods----------------
//-----------------------------------------------------------------------------

void FD1DRetSolver::interpBoundaryLow(int l, double* p){

    double ratio = 1;
    //using v_dx
    if ((engine->isVariableGrid == true ) ||
        ((engine->isVariableGrid == false) && (whichBarrierMethod == FD1DRet::VAR_GRID)) ) {
        ratio = v_dx[l+1] / v_dx[l+2];
    }

    p[l] = (1.0 + ratio)* p[l+1] - ratio * p[l+2];

    if (engine->DEBUG_forceNonNegative){
        p[l] = Maths::max(p[l], 0.0);
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::interpBoundaryUp(int t, double* p){

    double ratio = 1;
    //using v_dx
    if ((engine->isVariableGrid == true ) ||
        ((engine->isVariableGrid == false) && (whichBarrierMethod == FD1DRet::VAR_GRID)) ) {

        ratio = v_dx[t] / v_dx[t-1];
    }

    p[t] = (1.0 + ratio)* p[t-1] - ratio * p[t-2];

    if (engine->DEBUG_forceNonNegative){    
        p[t] = Maths::max(p[t], 0.0);
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::euroOneStepwithSource1D(
    int lastIndexTop,
    int* bot1, int* top1,
    int k, int idx,
    double** mSrcX,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP )
{
    bool solveByLine = true;
    
    //------------------------- 
    //compute the source terms
    //-------------------------
    int l1 = bot1[k];
    int t1 = top1[k];

    calcSourceAll( solveMethod, l1, t1, idx, mSrcX[k], mlastP );

    //calculate single price
    //xNum need to adjust if we deal with xNum+4       
    FDUtils::euroOneStepWithSource(
        l1, t1,
        &alphaX[0], &betaX[0], &gammaX[0],
        &aX[0], &bX[0], &cX[0],
        mcurrP, mlastP,
        mSrcX[k], solveByLine, 0, 0 );
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::calcSourceAll(const FD1DRet::TFdSolveType solveMethod,
                                int bot1, int top1,
                                int idx,
                                double* vSrc, 
                                double* vlastP
                                ){

    //init all source terms at boundaries to 0
    vSrc[bot1] = 0.0;
    vSrc[top1] = 0.0;

    switch (solveMethod) {
    case FD1DRet::DEFAULT:{
            //in fact, only calc low+1 to top1 -1
            calcSource(bot1, top1, idx, vSrc, vlastP);

            //this is only for explicite method
//            for (int index_x = bot1; index_x <= top1; index_x++ ) {
//                vSrc[index_x][0] = 0;                        
//            }

        }        
        break;
        default:{
            throw ModelException("FD1DRetSolver", "solveMethod != DEFAULT is not available");
        }
        break;
    }

    //adjusted source term if there is a jumps
    //FD2FSVCJ* modelSVCJ = dynamic_cast< FD2FSVCJ*>(engine);
    //if (modelSVCJ->isPureHeston == false){

    if (engine->hasJumps == true){
        //put jumps explicitely when do ADI on X direction

        //FD1DRetSolverJumpsSP FD1DRetSolverjumps(new FD1DRetSolverJumps(engine, xNum, yNum));
        //FD1DRetSolverJumps FD1DRetSolverjumps(engine, xNum, yNum);

        //FD1DRetSolverjumps->calcSource_jumps(vlastP, vSrc, dt, dx, dy);

    }
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::calcSource(int bot1, int top1,
                                int idx,
                                double*  vSrc, 
                                double* vlastP
                                ){
    int index_x;
    int l1 = bot1 +1;
    int t1 = top1 -1;

    //to double check, is my g = coupon term????aaaaaa
    for (index_x = l1; index_x <= t1; index_x++) {
        //vSrc[index_x][0] = dt * g[idx][index_x];
        //dt is moved to the model        
        vSrc[index_x] =  g[idx][index_x];
    }
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::calcCoeffPdeAll_X(const FD1DRet::TFdSolveType solveMethod, 
                         int low, int top, int idx){                                  

    switch (solveMethod) {
    case FD1DRet::DEFAULT:{
            if ((engine->isVariableGrid == true ) ||
                ((engine->isVariableGrid == false) && (whichBarrierMethod == FD1DRet::VAR_GRID)) ) {
                calcCoeffPdeVardx_X(low, top, idx);
            }else{
                calcCoeffPdeConstdx_X(low, top, idx);      
            }
        }        
        break;
        default:{
            throw ModelException("FD1DRetSolver", "solveMethod = Default is not available");
        }
        break;
    }
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::calcCoeffPdeConstdx_X(int low, int top, 
                                  int idx){
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
/*
    opt1 = 0.5 * dt / dx;
    opt2 = dt / (dx * dx);

    for (index_x = low; index_x <= top; index_x++) {
        Drift_Curr =  opt1 * a[idx][index_x];  
        ItoTerm_Curr =  opt2 * c[idx][index_x]; 

        Drift_Last =  opt1 * a[last_idx][index_x];  
        ItoTerm_Last =  opt2 * c[last_idx][index_x]; 

        alphaX[index_x] = theta * (Drift_Curr - ItoTerm_Curr);
        aX[index_x]=(1.0 - theta) * ( - Drift_Last 
                        + ItoTerm_Last);

        gammaX[index_x] = - theta * (Drift_Curr + ItoTerm_Curr);
        cX[index_x]= (1.0 - theta)* (Drift_Last + ItoTerm_Last);

        betaX[index_x] =1.0 + theta 
                        * ( 2.0 * ItoTerm_Curr - dt * f[idx][index_x]);
        bX[index_x]= 1.0 - (1.0-theta) * (2.0 * ItoTerm_Last 
                        - dt * f[last_idx][index_x]);             
    }

*/
    //mv dt back to model
//    opt1 = 0.5 * dt / dx;
//    opt2 = dt / (dx * dx);

    //to change

    if (engine->isVariableGrid == true ){// variable grid
        throw ModelException("FD1DRetSolver::calcCoeffPdeConstdx_X", "This function shouldn't be called when you need a variable FD Grid!");
    }

    //here, v_dx[i]= v_dx[j]=v_dx[1] =dx
    double dx = v_dx[1];

    opt1 = 0.5  / dx;
    opt2 = 1 / (dx * dx);

    for (index_x = l1; index_x <= t1; index_x++) {

        Drift_Curr =  opt1 * a[idx][index_x];  
        ItoTerm_Curr =  opt2 * c[idx][index_x]; 

        Drift_Last =  opt1 * a[last_idx][index_x];  
        ItoTerm_Last =  opt2 * c[last_idx][index_x]; 

        alphaX[index_x] = theta * (Drift_Curr - ItoTerm_Curr);
        aX[index_x]=(1.0 - theta) * ( - Drift_Last 
                        + ItoTerm_Last);

        gammaX[index_x] = - theta * (Drift_Curr + ItoTerm_Curr);
        cX[index_x]= (1.0 - theta)* (Drift_Last + ItoTerm_Last);

//        betaX[index_x] = 1 + theta 
//                        * ( 2.0 * ItoTerm_Curr - dt * f[idx][index_x]);
//        bX[index_x]= 1 - (1.0-theta) * (2.0 * ItoTerm_Last 
//                        - dt * f[last_idx][index_x]);             

        betaX[index_x] = 1 + theta 
                        * ( 2.0 * ItoTerm_Curr -  f[idx][index_x]);
        bX[index_x]= 1 - (1.0-theta) * (2.0 * ItoTerm_Last 
                        -  f[last_idx][index_x]);          
    }
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::calcCoeffPdeVardx_X(int low, int top, 
                                  int idx){
    int index_x;
    double Drift_Curr;
    double ItoTerm_Curr;
    double Drift_Last;
    double ItoTerm_Last;

    double opt1;
    double opt2;  
    double dxi2;
    double dxiplus12;

    int last_idx = 1- idx;
    int l1 = low +1;
    int t1 = top -1;


    for (index_x = l1; index_x <= t1; index_x++) {
        dxi2 = v_dx[index_x] * v_dx[index_x];
        dxiplus12 = v_dx[index_x+1] * v_dx[index_x+1];

        opt1 = v_dx[index_x] * dxiplus12 + v_dx[index_x+1] * dxi2;
        //opt2 = dt / opt1;
        opt2 = 1 / opt1;

        Drift_Curr =  opt2 * a[idx][index_x];  
        ItoTerm_Curr = 2.0 * opt2 * c[idx][index_x]; 


        Drift_Last =  opt2 * a[last_idx][index_x];  
        ItoTerm_Last = 2.0 * opt2 * c[last_idx][index_x]; 


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
                        -  f[idx][index_x]);
        bX[index_x]= 1 - (1.0-theta) * ( Drift_Last * (dxi2 - dxiplus12) +
                        ItoTerm_Last *(v_dx[index_x] + v_dx[index_x+1])
                        -  f[last_idx][index_x]);             
    }
}


//-----------------------------------------------------------------------------    

//-----------------------------------------------------------------------------
//----------------------special barriers---------------------------------------
//-----------------------------------------------------------------------------



/**------------------------------------------------------------
    set up the value for solver BotDim, TopDim
    for barriers: get barrier levels and looking for the cloest index
    it's called at updateSolverInfo 
------------------------------------------------------------*/

void FD1DRetSolver::setFDSolverDims(int step){
    //we need to recalc the botDim, topDim due to needSpecialFD == true
    //Our FD scheme size may be changing at each time step
    //reset the index to the initial dimension at each step
    //since the index may chg at previous step when there is a barrier
    int i;

    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        int nSlices = products[ prodIndex ]->getSlicesToDEV().size();
        for (i = 0; i < nSlices; i++){
            engine->botDim1[prodIndex][i] = 0;
            engine->topDim1[prodIndex][i] = xNum -1;
        }
        engine->botDim1[prodIndex][engine->maxNumOfValue ] = 0;  //should be the smallest of botDim1
        engine->topDim1[prodIndex][engine->maxNumOfValue ] = xNum -1;  //should be the higest of topDim1
    }
    setFDSolverDimsSpecial(step);
}

//----------------------------------------------------------------

void FD1DRetSolver::setFDSolverDimsSpecial(int step){
    int i;
    int minI = xNum-1;
    int maxI = 0;
    bool atNodeTop = false;
    bool atNodeDown = false;

    const TreeSliceEQ & undValue1Org =
        static_cast< const TreeSliceEQ & >( engine->payoffIndexOrig->getValue( step ) );
    int bot, top;
    undValue1Org.getCalcRange( bot, top );

    //barrier need be changed to nb of prod also, to review
    int prodIndex;
    const FDProductArray & products = engine->getProducts();

    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( products[ prodIndex ]->isElementary() )
            continue;

        int pStart = 0;
        int pEnd = products[ prodIndex ]->getSlicesToDEV().size() - 1;

        for (i = pStart; i <= pEnd; i++){
            int low = engine->botDim1[prodIndex][i];
            int high = engine->topDim1[prodIndex][i];

            //calc the closest index smaller or equal to bar, 
            //set the hasDownBarAtStep, hasUpBarAtStep based on the level of barrier 
            engine->barrier[i]->calcClosestIndex(low, high, undValue1Org, top - bot +1);

            //convert Bar to BarGrid based on chg of variable of model
            engine->convertToGridBar(step, i);

            //corrected numerical pb
            //low and high are all at the left side of barrier
            if (high +1 < xNum){
                // test UMR
                if(engine->barrier[i]->hasUpBarAtStep == true){
                    if (engine->gridLevel1[high +1] - engine->barrier[i]->upBarrierGrid < FP_MIN){
                        high = high +1; 
                        engine->barrier[i]->topDimBar = high;
                        atNodeTop = true; //barrier is at fd grid pts
                    }
                }
            }

            if (low +1 < xNum){

                // test UMR
                if(engine->barrier[i]->hasDownBarAtStep == true){
                    if (engine->gridLevel1[low +1] - engine->barrier[i]->downBarrierGrid < FP_MIN){
                        low = low +1; 
                        engine->barrier[i]->botDimBar = low;
                        atNodeDown= true; //barrier is at fd grid pts
                    }
                }
            }

            if (engine->barrier[i]->hasDownBarAtStep == true){
                engine->botDim1[prodIndex][i] = low;

                if  (  (engine->DEBUG_NODE_GEN == false)
                    && (whichBarrierMethod == FD1DRet::FIX_GRID) 
                    && (atNodeDown == false) ){
                    engine->botDim1[prodIndex][i] = low + 1;
                }
            }
            if (engine->barrier[i]->hasUpBarAtStep == true){
                engine->topDim1[prodIndex][i] = high;

                if  (  (engine->DEBUG_NODE_GEN == false)
                    && (whichBarrierMethod == FD1DRet::VAR_GRID) 
                    && (atNodeTop == false) ){
                    engine->topDim1[prodIndex][i] = high +1;
                }
            }

            //for test            
            engine->botDim1[prodIndex][i] = Maths::min(engine->botDim1[prodIndex][i], xNum-1);
            engine->topDim1[prodIndex][i] = Maths::min(engine->topDim1[prodIndex][i], xNum-1);
            minI = Maths::min(minI, engine->botDim1[prodIndex][i]);
            maxI = Maths::max(maxI, engine->topDim1[prodIndex][i]);
        }
        engine->botDim1[prodIndex][engine->maxNumOfValue] = minI;  //should be the smallest of botDim1
        engine->topDim1[prodIndex][engine->maxNumOfValue] = maxI;  //should be the higest of topDim1
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::setInsNodeInfo(int pStart, int pEnd){
    int k;
    int pEndTemp = pEnd;

    if (engine->DEBUG_NODE_GEN == false){
        //assume p[0] will be the price we want
        //FD grid and bar level are based on p[0]
        pEndTemp = pStart;
    }

    for (k = pStart; k <= pEndTemp; k++){
        int i = 2 * (k-pStart);

        //down
        if ((engine->barrier[k]->hasDownBarAtStep == true) ){
            insNode_FdIndex[i] = engine->barrier[k]->botDimBar;
            (*insNode)[0][i] = engine->barrier[k]->downBarrier;
            insNodeActive[i] = true;
        }else{
            //put value here only for the purpose to be able to call payoffBCD only one time
            insNode_FdIndex[i] = 0;
            (*insNode)[0][i] = engine->barrier[k]->downBarrier;
            insNodeActive[i] = false;
        }

        //up bar
        if ((engine->barrier[k]->hasUpBarAtStep == true) ){
            insNode_FdIndex[i+1] = engine->barrier[k]->topDimBar;
            (*insNode)[0][i+1] = engine->barrier[k]->upBarrier;
            insNodeActive[i+1] = true;
        }else{
            insNode_FdIndex[i+1] = xNum-1;
            (*insNode)[0][i+1] = engine->barrier[k]->upBarrier;
            insNodeActive[i+1]= false;
        }
    }
}

//-----------------------------------------------------------------------------

//#define DEBUG
//#include <dprintf.h>
//#undef min
//#undef max


void FD1DRetSolver::calcBCDsBefore(
    int step, 
    int* bot1, int* top1,
    int pStart, int pEnd, int k,
    const TreeSliceEQ & mlastP )
{
    int index;
    double x1, x2, y1, y2;

    const TreeSliceEQ & undValue1Org =
        static_cast< const TreeSliceEQ & >( engine->payoffIndexOrig->getValue( step ) );

    for (int kk=0; kk < int (insNode_FdIndex.size()); kk++){
        if (insNodeActive[kk] = true){
            index = insNode_FdIndex[kk];
            //for test URM
            if (index < xNum-1){
                x1 = undValue1Org[index];
                x2 = undValue1Org[index+1];
                y1= mlastP[index];
                y2= mlastP[index+1];
                (*insNodePrice)[0][k][kk] = LinearInterp((*insNode)[0][kk], x1, x2, y1, y2);
            }else{
                (*insNodePrice)[0][k][kk] = mlastP[xNum-1];
            }
        }
    }

    //seems only need these for var grid
    if (whichBarrierMethod == FD1DRet::VAR_GRID){        
        ajustLastPatInsNodes( bot1, top1, pStart, pEnd, k, insNodePriceBak, mlastP);
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::ajustLastPatInsNodes(
    int* bot1, int* top1,
    int pStart, int pEnd, int k,
    vector<vector<double> >& pBackup,
    const TreeSliceEQ & mlastP )
{
    int kk, i;
    int pEndTemp =pEnd;

    if (engine->DEBUG_NODE_GEN == false){
        pEndTemp = pStart;
    }

    for (kk = pStart; kk <= pEndTemp; kk++ ){
        i = 2 * (kk - pStart);
        if ((engine->barrier[kk]->hasDownBarAtStep == true) ){
            int downIndex = bot1[kk];

            pBackup[k][i] = mlastP[downIndex];
            //need to ajust the last vector price at the bar point if use VAR_GRID.
            mlastP[downIndex] = (*insNodePrice)[0][k][i];
        }

        if ((engine->barrier[kk]->hasUpBarAtStep == true) ){
            int upIndex = top1[kk];

            pBackup[k][i + 1] = mlastP[upIndex];
            mlastP[upIndex] = (*insNodePrice)[0][k][i + 1];
        }
    }
}

//-----------------------------------------------------------------------------
//this ft should be called every time when we want to reset dx 
//need to rearrange the index
//for now, only called when whichBarrierMethod ==VAR_GRID

//this ft need to be reviewed for the case if barriers are too close.
void FD1DRetSolver::resetdx(vector<double>& vdx, int kkProd, int pStart, int pEnd  ){

    int i, k, index;
    int pEndTemp = pEnd;

    vector<bool > vdxOverwrite;
    vector<double > vdxGridLevel;

    vdxGridLevel.resize(xNum, 0.0);
    vdxOverwrite.resize(xNum, false);

    vdx[0] = 0;
    for (i = 1; i < (int) vdx.size(); i++){
        vdx[i] = engine->v_dxMOrg[i];
    }

    //is set based on info obtained from p0
    if (engine->DEBUG_NODE_GEN == false){
        pEndTemp = pStart;
        for (k = pStart; k <= pEndTemp; k++){
            //is set based on info obtained from p0
            if (engine->barrier[k]->hasUpBarAtStep == true){
                //if more than one price, need to change v_dxM, one per price
                index = engine->topDim1[kkProd][k];            
                if ( index < xNum-1){
                    vdx[index] = engine->barrier[k]->upBarrierGrid - engine->gridLevel1[index-1];
                    vdx[index+1] = engine->gridLevel1[index+1] - engine->barrier[k]->upBarrierGrid; 
                }
            }

            if (engine->barrier[k]->hasDownBarAtStep == true){
                index = engine->botDim1[kkProd][k];
                if ( index > 0){
                    vdx[index] =  engine->barrier[k]->downBarrierGrid - engine->gridLevel1[index-1] ;
                    vdx[index+1] =  engine->gridLevel1[index+1] - engine->barrier[k]->downBarrierGrid ;
                }
            }
        }
    }else{

        for (k = pStart; k <= pEndTemp; k++){
            //is set based on info obtained from p0
            if (engine->barrier[k]->hasUpBarAtStep == true){
                //if more than one price, need to change v_dxM, one per price
                index = engine->topDim1[kkProd][k];            
                if ( index < xNum-1){

                    if ((vdxOverwrite[index] == true) && (vdxOverwrite[index+1] == true) ){ //two barriers within the same interval

                    }else{
                        if(vdxOverwrite[index] == false){ //haven't overwritten yet at this step
                            vdx[index] = engine->barrier[k]->upBarrierGrid - engine->gridLevel1[index-1];
                            vdxOverwrite[index] = true;
                            vdxGridLevel[index] = engine->barrier[k]->upBarrierGrid;
                        }else{
                            vdx[index] = engine->barrier[k]->upBarrierGrid - vdxGridLevel[index-1];
                            vdxGridLevel[index] = engine->barrier[k]->upBarrierGrid;
                        }

                        if(vdxOverwrite[index+1] == false){ //haven't overwritten yet at this step
                            vdx[index+1] = engine->gridLevel1[index+1] - engine->barrier[k]->upBarrierGrid; 
                            vdxOverwrite[index+1] = true;
                            vdxGridLevel[index] = engine->barrier[k]->upBarrierGrid;
                        }else{
                            vdx[index] = vdxGridLevel[index+1] - engine->barrier[k]->upBarrierGrid; 
                            vdxGridLevel[index] = engine->barrier[k]->upBarrierGrid;
                        }
                    }
                }
            }

            if (engine->barrier[k]->hasDownBarAtStep == true){
                index = engine->botDim1[kkProd][k];
                if ( index > 0){

                    if ((vdxOverwrite[index] == true) && (vdxOverwrite[index+1] == true) ){ //two barriers within the same interval

                    }else{
                        if(vdxOverwrite[index] == false){ //haven't overwritten yet at this step
                            vdx[index] =  engine->barrier[k]->downBarrierGrid - engine->gridLevel1[index-1] ;
                            vdxOverwrite[index] = true;
                            vdxGridLevel[index] = engine->barrier[k]->downBarrierGrid;
                        }else{
                            vdx[index] = engine->barrier[k]->downBarrierGrid - vdxGridLevel[index-1];
                            vdxGridLevel[index] = engine->barrier[k]->downBarrierGrid;
                        }

                        if(vdxOverwrite[index+1] == false){ //haven't overwritten yet at this step
                            vdx[index+1] = engine->gridLevel1[index+1] - engine->barrier[k]->downBarrierGrid ;
                            vdxOverwrite[index+1] = true;
                            vdxGridLevel[index] = engine->barrier[k]->downBarrierGrid;
                        }else{
                            vdx[index+1] = vdxGridLevel[index+1] - engine->barrier[k]->downBarrierGrid; 
                            vdxGridLevel[index] = engine->barrier[k]->downBarrierGrid;
                        }
                    }
                }
            }
        }
    }

    //also change model v_dxM
    for (i = 0; i < (int) vdx.size(); i++){
        engine->v_dxM[i] = vdx[i];
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::nextStepFdBoundaryBar(
    int step, 
    int* bot1, int* top1, 
    int k, double pv_x,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP)
{
    int upIndex ;
    int downIndex;

    //upBar
    if (engine->barrier[k]->hasUpBarAtStep == true){
        if (whichBarrierMethod == FD1DRet::VAR_GRID){
            upIndex = top1[k];
            mcurrP[upIndex] = engine->barrier[k]->valueAtUpBar;
        }else{
            upIndex = top1[k];
            downIndex = bot1[k];

            //don't chg grid, but need to calc nre rebate
            engine->barrier[k]->adjustCBDs(engine->gridLevel1, 
                        downIndex, upIndex,  //these two index are the cloest index to the critical points
                        mcurrP, mlastP );
        }
    }else{
        int t1 = top1[k];
        mcurrP[t1] =  pv_x * mlastP[t1];
    }

    //down Bar
    if (engine->barrier[k]->hasDownBarAtStep == true){
        if (whichBarrierMethod == FD1DRet::VAR_GRID){
            downIndex = bot1[k];
            mcurrP[downIndex] = engine->barrier[k]->valueAtDownBar;                                
        }else{
            //already done in "engine->getRebate", to see
            upIndex = top1[k];
            downIndex = bot1[k];
            engine->barrier[k]->adjustCBDs(engine->gridLevel1, 
                        downIndex, upIndex,  //these two index are the cloest index to the critical points
                        mcurrP, mlastP);
        }    
    }else{ //no bar 
        int l1 = bot1[k];
        mcurrP[l1] = pv_x * mlastP[l1];
    }                                        
}

//-----------------------------------------------------------------------------    

void FD1DRetSolver::updateFdBoundaryBar(
    int step,
    int* bot1, int* top1,
    int pStart, int pEnd,
    int k, double pv_x,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP )
{
    int upIndex ;
    int downIndex;

    //var grid at BCDs
    if (whichBarrierMethod == FD1DRet::VAR_GRID){
        calcBCDsAfterVarGrid( step, bot1, top1, pStart, pEnd, k, mcurrP, mlastP );
    }else{
        //recalc new rebate
        calcBCDsAfterFixGrid( step, bot1, top1, pStart, pEnd, k, mcurrP, mlastP );
    }
            
    //down barrier
    if (engine->barrier[k]->hasDownBarAtStep == true) {
        if (whichBarrierMethod == FD1DRet::FIX_GRID){//fixed grid
            upIndex = top1[k];
            downIndex = bot1[k];
            engine->barrier[k]->adjustCBDs(engine->gridLevel1, 
                                            downIndex, upIndex,  //these two index are the cloest index to the critical points
                                            mcurrP, mcurrP);
        }
    }else{
        int l1 = bot1[k];
        interpBoundaryLow(l1, mcurrP);
    }
    
    //up barrier
    if ((engine->barrier[k]->hasUpBarAtStep == true) ){                        
        if (whichBarrierMethod == FD1DRet::FIX_GRID){//fixed grid
            upIndex = top1[k];
            downIndex = bot1[k];
            engine->barrier[k]->adjustCBDs(engine->gridLevel1, 
                                            downIndex, upIndex,    //these two index are the cloest index to the critical points
                                            mcurrP, mcurrP);
        }            
    }else{
        int t1 = top1[k];
        interpBoundaryUp(t1, mcurrP);
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::calcBCDsAfterFixGrid(
    int step,
    int* bot1, int* top1,
    int pStart, int pEnd, int k,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP )
{
    //want use the current price vector to get the the current rebate
    //maybe need to rewrite this ft.
    //use calcBCDsBefore isn't opt. 

    calcBCDsBefore( step,  bot1, top1, pStart, pEnd, k, mcurrP );
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::calcBCDsAfterVarGrid(
    int step, 
    int* bot1, int* top1, 
    int pStart, int pEnd, int k,
    const TreeSliceEQ & mcurrP,
    const TreeSliceEQ & mlastP)
{

    //fisrt copy back the right old price to the borrowered index at last Step
    copyBackLastPatInsNodes( bot1, top1, pStart, pEnd, k, insNodePriceBak, mlastP );
    
    int h, l;
    int i, kk;
    int pEndTemp = pEnd;

    if (engine->DEBUG_NODE_GEN == false){
        pEndTemp = pStart;
    }

    for (kk = pStart; kk <= pEndTemp ; kk++){
        i = 2 * (kk-pStart);
        h= top1[kk];
        l = bot1[kk];

        (*insNodePrice)[0][k][i + 0] = mcurrP[l];  //low bar
        (*insNodePrice)[0][k][i + 1] = mcurrP[h];  //high bar
        engine->barrier[k]->valueAtDownBar = (*insNodePrice)[0][k][i];
        engine->barrier[k]->valueAtUpBar = (*insNodePrice)[0][k][i+1];
    }

    //get approximated BCDs from the product at inserted nodes   
/*    
    // for each product
    for( prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        if( products[ prodIndex ]->isElementary() )
            continue;

        products[ prodIndex ]->update( step, FDProduct::BWD_NODE_INSERTION );
    }

    for (kk = pStart; kk <= pEndTemp; kk++ ){
        i = 2 * (kk-pStart);

        if ((engine->barrier[kk]->hasDownBarAtStep == true) ){
            for (k = pStart; k <= pEnd; k++){
                engine->barrier[k]->valueAtDownBar =  (*insNodePrice)[0][k][i];
            }
        }

        if ((engine->barrier[kk]->hasUpBarAtStep == true) ){
            for (k = pStart; k <= pEnd; k++){
                engine->barrier[k]->valueAtUpBar =  (*insNodePrice)[0][k][i+1];
            }
        }        
    }
*/

    //get the ptr of spotOrg
/*
        if( bot1[k] >= 1 )
            mcurrP[ bot1[k] - 1 ] = 0.;
        if( bot1[k] >= 2 )
            mcurrP[ bot1[k] - 2 ] = 0.;
        if( top1[k] <= xNum - 1 )
            mcurrP[ top1[k] + 1 ] = 0.;
        if( top1[k] <= xNum - 2 )
            mcurrP[ top1[k] + 2 ] = 0.;
*/
    const TreeSliceEQ & undValue1Org =
        static_cast< const TreeSliceEQ & >( engine->payoffIndexOrig->getValue( step ) );
        
    double x, x1, x2, y1, y2;

    //need to set the correct value at the "borrowed index J, K"
    for (kk = pStart; kk <= pEndTemp; kk++){
        h = top1[kk];
        l = bot1[kk];
        i = 2 * (kk - pStart);

        if ((engine->barrier[kk]->hasDownBarAtStep == true) )
        {
            x = undValue1Org[l];
            x1 = (*insNode)[0][i];
            y1 = (*insNodePrice)[0][k][i];
            if (l ==0 ){
                x2 = undValue1Org[l+1];
                y2 = mcurrP[l+1];
            }else{
                x2 = undValue1Org[l-1];
                y2 = mcurrP[l-1];
            }
            mcurrP[l] = LinearInterp(x, x1, x2, y1, y2);
        }

        if ((engine->barrier[kk]->hasUpBarAtStep == true) )
        {
            x = undValue1Org[h];
            x1 = (*insNode)[0][i + 1];
            y1 = (*insNodePrice)[0][k][i + 1];

            if (engine->DEBUG_NODE_GEN == false){//DblBarrier case
                if(h+1 < xNum){
                    x2 = undValue1Org[h+1];
                    y2 = mcurrP[h+1];
                }else{
                    x2 = undValue1Org[h];
                    y2 = mcurrP[h];
                }
            }else{
                if (h ==0 ){
                    x2 = undValue1Org[h+1];
                    y2 = mcurrP[h+1];
                }else{
                    x2 = undValue1Org[h-1];
                    y2 = mcurrP[h-1];
                }                
            }
            mcurrP[h] = LinearInterp(x, x1, x2, y1, y2);
        }        
    }
}

//-----------------------------------------------------------------------------

void FD1DRetSolver::copyBackLastPatInsNodes(
    int* bot1, int* top1,
    int pStart, int pEnd, int k,
    vector<vector<double> >& pBackup,
    const TreeSliceEQ & mlastP )
{
    int i, kk;

    int pEndTemp = pEnd;

    if (engine->DEBUG_NODE_GEN == false){
        pEndTemp = pStart;
    }

    for (kk = pStart; kk <= pEndTemp; kk++ ){
        i = 2 * (kk-pStart);
        if ((engine->barrier[kk]->hasDownBarAtStep == true) ){
            int downIndex = bot1[kk];
            mlastP[downIndex] = pBackup[k][i];
        }

        if ((engine->barrier[kk]->hasUpBarAtStep == true) ){
            int upIndex = top1[kk];
            mlastP[upIndex] = pBackup[k][i + 1];
        }
    }   
}

DRLIB_END_NAMESPACE
