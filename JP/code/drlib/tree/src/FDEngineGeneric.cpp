//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FDEngineGeneric.cpp
//
//   Description : One factor generic finite difference engine
//
//   Author      : André Segger
//
//   Date        : 14 October 2003
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"            // needed for precompiled headers
#include "edginc/Maths.hpp"
#include "edginc/FDUtils.hpp"
#include "edginc/FD1FE2C.hpp"
#include "edginc/FDModel.hpp"

DRLIB_BEGIN_NAMESPACE


FDEngine1FGeneric::FDEngine1FGeneric(){
    stockArray              = 0;
    stockArrayNew           = 0;
    fwdGrid                 = 0;
    assetInterpLevel        = 0;

    dt                      = 0;
    forwards                = 0;
    pvs                     = 0;
    ir                      = 0;
    divy                    = 0;

    segEnd                  = 0;
    segMax                  = 0;
    segMin                  = 0;
    segStrike               = 0;

    upBarrierNew            = 0;
	upPayoutNew             = 0;
	upPayoutDeltaNew        = 0;
	upBarrierOld            = 0;
	upPayoutOld             = 0;
	upPayoutDeltaOld        = 0;
	downBarrierNew          = 0;
	downPayoutNew           = 0;
	downPayoutDeltaNew      = 0;
	downBarrierOld          = 0;
	downPayoutOld           = 0;
	downPayoutDeltaOld      = 0;

    useFwdGrid = false;

    numPriceArrays = -1;
}

FDEngine1FGeneric::~FDEngine1FGeneric()
{
    Clear();
}

/** clean up */
void FDEngine1FGeneric::Clear()
{
    if (stockArray != 0){
        delete [] stockArray;
        stockArray = 0;
    }    
    if (stockArrayNew != 0){
        delete [] stockArrayNew;
        stockArrayNew = 0;
    }
    if (assetInterpLevel != 0){
        delete [] assetInterpLevel;
        assetInterpLevel = 0;
    }
    if (fwdGrid != 0){
        delete [] fwdGrid;
        fwdGrid = 0;
    }
    if (dt != 0){
        delete [] dt;
        dt = 0;
    }
    if (forwards != 0){
        delete [] forwards;
        forwards = 0;
    }
    if (pvs != 0){
        delete [] pvs;
        pvs = 0;
    }
    if (ir != 0){
        delete [] ir;
        ir = 0;
    }
    if (divy != 0){
        delete [] divy;
        divy = 0;
    }    
    if (segEnd != 0){
        delete [] segEnd;
        segEnd = 0;
    }    
    if (segMax != 0){
        delete [] segMax;
        segMax = 0;
    }    
    if (segMin != 0){
        delete [] segMin;
        segMin = 0;
    }    
    if (segStrike != 0){
        delete [] segStrike;
        segStrike = 0;
    }    
    if (upBarrierNew != 0){
        delete [] upBarrierNew;
        upBarrierNew = 0;
    }  
    if (upPayoutNew != 0){
        delete [] upPayoutNew;
        upPayoutNew = 0;
    }   
    if (upPayoutDeltaNew != 0){
        delete [] upPayoutDeltaNew;
        upPayoutDeltaNew = 0;
    }  
    if (upBarrierOld != 0){
        delete [] upBarrierOld;
        upBarrierOld = 0;
    }   
    if (upPayoutOld != 0){
        delete [] upPayoutOld;
        upPayoutOld = 0;
    }   
    if (upPayoutDeltaOld != 0){
        delete [] upPayoutDeltaOld;
        upPayoutDeltaOld = 0;
    }   
    if (downBarrierNew != 0){
        delete [] downBarrierNew;
        downBarrierNew = 0;
    }   
    if (downPayoutNew != 0){
        delete [] downPayoutNew;
        downPayoutNew = 0;
    }   
    if (downPayoutDeltaNew != 0){
        delete [] downPayoutDeltaNew;
        downPayoutDeltaNew = 0;
    }
    if (downBarrierOld != 0){
        delete [] downBarrierOld;
        downBarrierOld = 0;
    }   
    if (downPayoutOld != 0){
        delete [] downPayoutOld;
        downPayoutOld = 0;
    }   
    if (downPayoutDeltaOld != 0){
        delete [] downPayoutDeltaOld;
        downPayoutDeltaOld = 0;
    }
}

void FDEngine1FGeneric::init(int timeSteps, int stockSteps, double *times, double *inForwards, double *inPVs,
                             int inNSegments, int *inSegEnd, double *inSegMax, double *inSegMin, double*inSegStrike,
                             bool inUseFwdGrid, double *inIr, double *inDivy, int inNumPriceArrays,int inGridType)
{
    int i;

    Clear();

    nTimeSteps  = timeSteps;
    nStockSteps = stockSteps;
    useFwdGrid  = inUseFwdGrid;
    dt = new double [nTimeSteps+1];
    for (i=0; i<=nTimeSteps; i++) {
        dt[i] = times[i];
    }

    if (useFwdGrid == true) {
        forwards = new double [nTimeSteps+1];
        pvs = new double [nTimeSteps+1];
        for (i=0; i<=nTimeSteps; i++) {
            forwards[i] = inForwards[i];
            pvs[i] = inPVs[i];
        }
    } else {
        forwards = new double [nTimeSteps+1];
        ir = new double [nTimeSteps+1];
        divy = new double [nTimeSteps+1];
        for (i=0; i<=nTimeSteps; i++) {
            forwards[i] = inForwards[i];
            ir[i] = inIr[i];
            divy[i] = inDivy[i];
        }
    }
   
    nSegments   = inNSegments;
    segEnd      = new int [nSegments];
    segMax      = new double [nSegments];
    segMin      = new double [nSegments];
    segStrike   = new double [nSegments];
    for (i=0; i<nSegments; i++) {
        segEnd[i] = inSegEnd[i];
        segMax[i] = inSegMax[i];
        segMin[i] = inSegMin[i];
        segStrike[i] = inSegStrike[i];
    }

    varMethod = 0;
    numPriceArrays = inNumPriceArrays;
    gridType = inGridType;

    FD1FGeneric* fdModel = dynamic_cast<FD1FGeneric*>(model);

    if ( fdModel->isE2C() && fdModel->hasEquityLayer() ) {
        // need one array for lambda unadjusted prices and one array for lambda adjusted prices
        // numPriceArrays = numPriceArrays+2;
    }


    AllocateArrays();
}


void FDEngine1FGeneric::AllocateArrays() {
    static const string method = "FDEngine1FGeneric::AllocateArrays";

    int i;

    stockArray       = new double[nStockSteps+1];
    stockArrayNew    = new double[nStockSteps+1];
    fwdGrid          = new double[nStockSteps+1];
    assetInterpLevel = new double[nStockSteps+1];
    
    range = TreeSliceGeneral::Range::create( 1, 0, nStockSteps );
    int dimBits = (1<<(*range)->nDim)-1;
    optionArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );
    optionOldArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );
    optionInterpArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );

    upBarrierNew            = new double [numPriceArrays];
	upPayoutNew             = new double [numPriceArrays];
	upPayoutDeltaNew        = new double [numPriceArrays];
	upBarrierOld            = new double [numPriceArrays];
	upPayoutOld             = new double [numPriceArrays];
	upPayoutDeltaOld        = new double [numPriceArrays];
	downBarrierNew          = new double [numPriceArrays];
	downPayoutNew           = new double [numPriceArrays];
	downPayoutDeltaNew      = new double [numPriceArrays];
	downBarrierOld          = new double [numPriceArrays];
	downPayoutOld           = new double [numPriceArrays];
	downPayoutDeltaOld      = new double [numPriceArrays];

    for (i=0; i<numPriceArrays; i++) {
        upBarrierNew[i]          = -1;
	    upPayoutNew[i]           =  0;
	    upPayoutDeltaNew[i]      =  0;
	    upBarrierOld[i]          = -1;
	    upPayoutOld[i]           =  0;
	    upPayoutDeltaOld[i]      =  0;
	    downBarrierNew[i]        = -1;
	    downPayoutNew[i]         =  0;
	    downPayoutDeltaNew[i]    =  0;
	    downBarrierOld[i]        = -1;
	    downPayoutOld[i]         =  0;
	    downPayoutDeltaOld[i]    =  0;
    }
}

void FDEngine1FGeneric::loop(FDPayoff1F *engineCBs, double stockNow, double *price, double *divPert, double *irPert) {
    static const string method = "FDEngine1F::loop";
    int i, j, k, l;
    FDSolver1FGeneric fdSolver;
    FDSolver1FGeneric fdAssetSolver;

    double *tmpPointer;
    double *assetArray              = NULL;

    double *bondFloors = NULL;
    double *bondDelta  = NULL;
    double *bondGamma  = NULL;

    double fwdNew;
    double fwdOld;
    double pvFact;
    int    currSeg;
    double delta, gamma;

    FDPayoff1FGeneric* genericPayoff = dynamic_cast<FDPayoff1FGeneric*>(engineCBs);

    if ( !genericPayoff) {
        throw ModelException(method, "Payoff must be of type FDPayoff1FGeneric.");
    }

    fdSolver.forceNonNegative(engineCBs->Positive());
    
    if (fdSolver.allocateMem(nStockSteps) == FAILURE) 
        throw ModelException(method, "fdSolver.allocateMem failure");

    if ( !FD1FGeneric::TYPE->isInstance(model)) {
        throw ModelException(method, "Model must be derived from FD1FGeneric");
    }

    FD1FGeneric* fdModel = dynamic_cast<FD1FGeneric*>(model);

    currSeg = nSegments-1;
    if (useFwdGrid == false) {
        if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
            throw ModelException(method, "fdSolver.createGrid failure");

        if ( fdModel->isE2C() && fdModel->hasEquityLayer() ) {
            if (fdAssetSolver.allocateMem(nStockSteps) == FAILURE) 
                throw ModelException(method, "fdAssetSolver.allocateMem failure");

            if (fdAssetSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
                throw ModelException(method, "fdSolver.createGrid failure");

            FD1FE2C* fdModelE2C = dynamic_cast<FD1FE2C*>(model);

            if (fdSolver.createGrid(fdModelE2C->equityInSegMin[currSeg], fdModelE2C->equityInSegMax[currSeg],
                                    nStockSteps,gridType,fdModelE2C->equityInSegStrike[currSeg]) == FAILURE) {
                throw ModelException(method, "fdAssetSolver.createGrid failure");
            }

            // get the stock price array
            tmpPointer = fdSolver.getSpotsPtr();
            for (j=0; j<=nStockSteps;j++){
                stockArray[j] = tmpPointer[j];
            }

            // calculate stock to asset mapping
            for (j=0 ; j<=nStockSteps ; ++j) {
                assetInterpLevel[j] = stockArray[j];
            }

            fdModelE2C->mapToEquitySpace(assetInterpLevel, nStockSteps+1);

            // allocate arrays
            bondFloors = new double[nStockSteps+1]; 
            bondDelta  = new double[nStockSteps+1];
            bondGamma  = new double[nStockSteps+1];

        } else {
            tmpPointer = fdSolver.getSpotsPtr();
            for (j=0; j<=nStockSteps;j++){
                stockArray[j] = tmpPointer[j];
            }
        }
    } else {
        if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
            throw ModelException(method, "fdSolver.createGrid failure");
        fwdNew = forwards[nTimeSteps];
        tmpPointer = fdSolver.getSpotsPtr();
        for (j=0; j<=nStockSteps;j++){
            fwdGrid[j] = tmpPointer[j];
            stockArray[j] = fwdNew*tmpPointer[j];
        }
    }

    if (!model) {
        throw ModelException(method, "The model must not be null");
    }

    int equityPriceArrays = numPriceArrays;
    int e2cPriceArrays    = 0;
    if ( fdModel->isE2C() && fdModel->hasEquityLayer() ) {
        double* grid = fdAssetSolver.getSpotsPtr();
        assetArray             = new double[nStockSteps+1];
        for (i=0 ; i <= nStockSteps ; ++i) {
            assetArray[i] = grid[i];
        }

        equityPriceArrays = equityPriceArrays-2;
        e2cPriceArrays = 2;
    }

    genericPayoff->preCalcFDGeneric(nTimeSteps, 0, 0, equityPriceArrays-1, (*optionArray)[0], (*optionOldArray)[0]);
    genericPayoff->PayoffAtMatFD(stockArray, nTimeSteps, 0, nStockSteps, 0, equityPriceArrays-1, (*optionArray)[0]);

    FDTermStructureSP driftTerm;
    FDTermStructureSP diffusionTerm;
    FDTermStructureSP discountTerm;
    FDTermStructureSP couponTerm;

    
    fdModel->preProcessGrid(stockArray, 0, nStockSteps);

    for (i = nTimeSteps; i>0; i--) {

        swapT( optionOldArray, optionArray );

        if (useFwdGrid == false) {
            genericPayoff->preCalcFDGeneric(i-1, 0, 0, equityPriceArrays-1, (*optionArray)[0], (*optionOldArray)[0]);

            if ( fdModel->isE2C() && fdModel->hasEquityLayer() ) {

                // pre-process E2C layer
                driftTerm       = fdModel->getDriftTerm    (i, assetArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], false);
                diffusionTerm   = fdModel->getDiffusionTerm(i, assetArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], false);
                couponTerm      = fdModel->getCouponTerm   (i, assetArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], false);
                discountTerm    = fdModel->getDiscountTerm (i, assetArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], false);

                fdSolver.forceNonNegative(true);

                if (fdAssetSolver.solveBarrier(nStockSteps, driftTerm, diffusionTerm, discountTerm, couponTerm, dt[i], 
                                               (*optionArray)[0][equityPriceArrays],(*optionOldArray)[0][equityPriceArrays], -1,0,0,
                                               -1,0,0,downBarrierNew[equityPriceArrays],downPayoutNew[equityPriceArrays],
                                               downPayoutDeltaNew[equityPriceArrays],downBarrierOld[equityPriceArrays],
                                               downPayoutOld[equityPriceArrays],downPayoutDeltaOld[equityPriceArrays],0,0) == FAILURE) {
                    throw ModelException(method, "fdAssetSolver.solveBarrier failure");
                } 

                // the lambda adjusted prices will be held in a separate array
                for (l=0 ; l<=nStockSteps ; ++l) {
                    (*optionArray)[0][equityPriceArrays+1][l] = (*optionArray)[0][equityPriceArrays][l];
                }

                if ( fdModel->doLambdaAdjust() ) {
                    double lambda = fdModel->getLambda();

                    if (Maths::isPositive(lambda)) {
	                    FDLambdaAdjustment (fdAssetSolver,nStockSteps,lambda,downBarrierNew[equityPriceArrays],downPayoutNew[equityPriceArrays],10,(*optionArray)[0][equityPriceArrays+1]);
                    }
                }

                // get the bond floors at the stock prices
                if (FDInterpolationD(nStockSteps+1,assetArray,(*optionArray)[0][equityPriceArrays+1],nStockSteps+1, assetInterpLevel,bondFloors,bondDelta,bondGamma) == FAILURE)
                    throw ModelException(method, "FDCubicSplineInterpD failure");

                for (l=0 ; l<=nStockSteps ; ++l) {
                    if ( assetInterpLevel[l] <= assetArray[nStockSteps] )
                        (*optionArray)[0][equityPriceArrays+1][l] = bondFloors[l];
                    else 
                        (*optionArray)[0][equityPriceArrays+1][l] = (*optionArray)[0][equityPriceArrays+1][nStockSteps];
                }


            }

            // run this again to get the correct bond floors at the boundaries - some scope for performance improvement here
            genericPayoff->preCalcFDGeneric(i-1, 0, 0, equityPriceArrays-1, (*optionArray)[0], (*optionOldArray)[0]);

            bool doEquityLayer = fdModel->hasEquityLayer();

            driftTerm       = fdModel->getDriftTerm    (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], doEquityLayer);
            diffusionTerm   = fdModel->getDiffusionTerm(i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], doEquityLayer);
            couponTerm      = fdModel->getCouponTerm   (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], doEquityLayer);
            discountTerm    = fdModel->getDiscountTerm (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], doEquityLayer);
            
            for (k=0; k<equityPriceArrays; k++) {
                if ( k == 0 ) {
                    fdSolver.forceNonNegative(engineCBs->Positive());
                } else {
                    fdSolver.forceNonNegative(true);
                }

                if (fdSolver.solveBarrier(nStockSteps, driftTerm, diffusionTerm, discountTerm, couponTerm, dt[i], 
                                          (*optionArray)[0][k],(*optionOldArray)[0][k], upBarrierNew[k],upPayoutNew[k],upPayoutDeltaNew[k],
                                          upBarrierOld[k],upPayoutOld[k],upPayoutDeltaOld[k],downBarrierNew[k],downPayoutNew[k],
                                          downPayoutDeltaNew[k],downBarrierOld[k],downPayoutOld[k],downPayoutDeltaOld[k],0,0) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }  
            }

            if (currSeg > 0 && i-1 == segEnd[currSeg-1]) {
                currSeg--;
                
                if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
                    throw ModelException(method, "fdSolver.createGrid failure");
                
                tmpPointer = fdSolver.getSpotsPtr();
                for (j=0; j<=nStockSteps;j++){
                   stockArrayNew[j] = tmpPointer[j];
                }

                
                for (k=0; k<equityPriceArrays; k++) {
                    // figure out the top stock level used
                    int nodeMax = nStockSteps;
                    for (j=nStockSteps; j>=0; j--) {
                        if(stockArray[j] < upBarrierNew[k]) {
                            nodeMax = j;
                            break;
                        }      
                    }

                    if (FDCubicSplineInterp(nodeMax+1,
                                            stockArray,
                                            (*optionArray)[0][k],
                                            nStockSteps+1,
                                            stockArrayNew,
                                            (*optionInterpArray)[0][k]) == FAILURE) 
                    {
                        throw ModelException(method, "FDCubicSplineInterp failure");
                    }
                }

                swapT( optionArray, optionInterpArray );

                tmpPointer = stockArray;
                stockArray = stockArrayNew;
                stockArrayNew = tmpPointer;

            }  

            engineCBs->PayoffBeforeMatFD(stockArray, i-1, 0, nStockSteps, 0, equityPriceArrays-1, (*optionArray)[0]);

        } else {
            fwdNew = forwards[i-1];
            fwdOld = forwards[i];
            for (j=0; j<=nStockSteps;j++) {
                stockArray[j] = fwdNew*fwdGrid[j];
            }
            pvFact = pvs[i];//DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);

            driftTerm       = fdModel->getDriftTerm    (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], true);
            diffusionTerm   = fdModel->getDiffusionTerm(i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], true);
            couponTerm      = fdModel->getCouponTerm   (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], true);
            discountTerm    = fdModel->getDiscountTerm (i, stockArray, 0, nStockSteps, useFwdGrid, irPert[i], divPert[i], true);

            diffusionTerm->scale(1./(fwdNew*fwdNew));
            driftTerm->scale(1./fwdNew);

            engineCBs->preCalcFD(i-1, 0, 0, equityPriceArrays-1);
            
            for (k=0; k<equityPriceArrays; k++) {
                if ( k == 0 ) {
                    fdSolver.forceNonNegative(engineCBs->Positive());
                } else {
                    fdSolver.forceNonNegative(true);
                }

                // to do: drift term should be ir and dividend perturbations only, as fwd otherwise are driftless
                if (fdSolver.solveBarrier(nStockSteps, driftTerm, diffusionTerm, discountTerm, couponTerm, dt[i], 
                                          (*optionArray)[0][k],(*optionOldArray)[0][k], upBarrierNew[k]/fwdNew,upPayoutNew[k],upPayoutDeltaNew[k]*fwdNew,
                                          upBarrierOld[k]/fwdOld,upPayoutOld[k],upPayoutDeltaOld[k]*fwdOld,downBarrierNew[k]/fwdOld,downPayoutNew[k],
                                          downPayoutDeltaNew[k]*fwdOld,downBarrierOld[k]/fwdOld,downPayoutOld[k],downPayoutDeltaOld[k]*fwdOld,0,0) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }  

                for (j=0; j<=nStockSteps;j++){
                    (*optionArray)[0][k][j] *= pvFact;
                }
            }
            
            if (currSeg > 0 && i-1 == segEnd[currSeg-1]) {
                currSeg--;                

               if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
                    throw ModelException(method, "fdSolver.createGrid failure");
                
                tmpPointer = fdSolver.getSpotsPtr();
                for (j=0; j<=nStockSteps;j++){
                    fwdGrid[j] = tmpPointer[j];
                    stockArrayNew[j] = fwdNew*tmpPointer[j];
                }  
                
                for (k=0; k<numPriceArrays; k++) {

                    // figure out the top stock level used
                    int nodeMax = nStockSteps;
                    for (j=nStockSteps; j>=0; j--) {
                        if(stockArray[j] < upBarrierNew[k]) {
                            nodeMax = j;
                            break;
                        }      
                    }
                

                    if (FDCubicSplineInterp(nodeMax+1,
                                            stockArray,
                                            (*optionArray)[0][k],
                   	                        nStockSteps+1,
                                            stockArrayNew,
                                            (*optionInterpArray)[0][k])== FAILURE) 
                    {
                        throw ModelException(method, "FDCubicSplineInterp failure");
                    }
                }

                swapT( optionArray, optionInterpArray );

                tmpPointer = stockArray;
                stockArray = stockArrayNew;
                stockArrayNew = tmpPointer;

            }  

            engineCBs->PayoffBeforeMatFD(stockArray, i-1, 0, nStockSteps, 0, equityPriceArrays-1, (*optionArray)[0]);
        }
    }

    // adjust the final bond price for lambda
    if ( fdModel->isE2C() ) {
        if ( !fdModel->hasEquityLayer()) {
            if ( fdModel->doLambdaAdjust() ) {
                double lambda = fdModel->getLambda();

                if (Maths::isPositive(lambda)) {
                    for (k=0; k<numPriceArrays; k++) {
	                    FDLambdaAdjustment (fdSolver,nStockSteps,lambda,downBarrierNew[k],downPayoutNew[k],10,(*optionArray)[0][k]);
                    }
                }
                stockNow = fdModel->getLambdaAdjustedSpot(stockNow);
            }
        } else {
            double lambda = fdModel->getLambda();
            // stockNow      = fdModel->getLambdaAdjustedSpot(stockNow);
            // stockNow      = stockNow - downBarrierNew[equityPriceArrays];
            FD1FE2C* fdModelE2C = dynamic_cast<FD1FE2C*>(model);

            stockNow      = (stockNow - downBarrierNew[equityPriceArrays]) * fdModelE2C->getFXRate();
        }
    }
    
    for (k=0; k<numPriceArrays; k++) {
        if ( stockNow > stockArray[nStockSteps] && ( upBarrierNew[k] >= 0 || 
             ( fdModel->isE2C() && fdModel->hasEquityLayer()))) {

		    // Correct values for spots above the barrier
			price[k] = (*optionArray)[0][k][nStockSteps] + 
                       upPayoutDeltaNew[k]*(stockNow - stockArray[nStockSteps]);
			delta = upPayoutDeltaNew[k];
			gamma = 0.;
        } else {
            if (FDInterpolationD(nStockSteps+1,stockArray,(*optionArray)[0][k],1,&stockNow,price+k,&delta,&gamma) == FAILURE)
                throw ModelException(method, "FDCubicSplineInterpD failure");
        }
    }

    fdModel->postProcessFD(price, equityPriceArrays);

    // de-allocate memory
    if (bondFloors)
        delete [] bondFloors;
    if (bondDelta)
        delete [] bondDelta;
    if (bondGamma)
        delete [] bondGamma;
    if (assetArray)
        delete [] assetArray;

    return;
}

void FDEngine1FGeneric::FDLambdaAdjustment (
    FDSolver1FGeneric&  fdSolver,
	int			        m,
	double		        lambda,
	double		        downBarrier,
	double		        downPayout,
	int			        n,
	double*		        V)
{
	double			delta_t = 1./n;

	fdSolver.solveQuick (m,0,0,lambda,0,delta_t,V,-1,0,downBarrier,downPayout,n);
}



DRLIB_END_NAMESPACE

