//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FDEngine.cpp
//
//   Description : One factor finite difference engine
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : February 7, 2002
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp" // needed for precompiled headers
#include "edginc/FDEngine.hpp"
#include "edginc/ZFDSolver1D.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FDModel.hpp"

DRLIB_BEGIN_NAMESPACE
//#define FD1F_DEBUG_FILE1 "c:\\Temp\\fd1f.txt"

#ifdef FD1F_DEBUG_FILE1


void Debug_OutPutFD1FPrice(int step, double time, 
                      double* s, double** price, 
                      int bot, int top, int nPrice, 
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
            sprintf(buffer, "#step = %d,  price \n", step);
            fprintf(DumpFile, buffer);
            for (int ix =bot; ix <=top; ix ++)
            {
                sprintf (buffer, "%f, %f, %f, # %d \n", time, s[ix],  price[i][ix], ix);
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


FDEngine1F::FDEngine1F(){
    stockArray = 0;
    stockArrayNew = 0;
 //   sigma = 0;
    fwdGrid = 0;

    dt = 0;
    forwards = 0;
    pvs = 0;
    ir = 0;
    divy = 0;

    segEnd = 0;
    segMax = 0;
    segMin = 0;
    segStrike = 0;

    upBarrierNew=0;
	upPayoutNew=0;
	upPayoutDeltaNew=0;
	upBarrierOld=0;
	upPayoutOld=0;
	upPayoutDeltaOld=0;
	downBarrierNew=0;
	downPayoutNew=0;
	downPayoutDeltaNew=0;
	downBarrierOld=0;
	downPayoutOld=0;
	downPayoutDeltaOld=0;
    
    useFwdGrid = false;

    numPriceArrays = -1;
}

FDEngine1F::~FDEngine1F()
{
    Clear();
}

/** clean up */
void FDEngine1F::Clear()
{
    if (stockArray != 0){
        delete [] stockArray;
        stockArray = 0;
    }    
    if (stockArrayNew != 0){
        delete [] stockArrayNew;
        stockArrayNew = 0;
    }
/*    if (sigma != 0){
        delete [] sigma;
        sigma = 0;
    }  */
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

void FDEngine1F::init(int timeSteps, int stockSteps, double *times, double *inForwards, double *inPVs,
                      int inNSegments, int *inSegEnd, double *inSegMax, double *inSegMin, double*inSegStrike,
                      bool inUseFwdGrid, double *inIr, double *inDivy, int inNumPriceArrays,int inGridType)
{
    int i;

    Clear();

    nTimeSteps = timeSteps;
    nStockSteps = stockSteps;
    useFwdGrid = inUseFwdGrid;
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
   
    nSegments = inNSegments;
    segEnd = new int [nSegments];
    segMax = new double [nSegments];
    segMin = new double [nSegments];
    segStrike = new double [nSegments];
    for (i=0; i<nSegments; i++) {
        segEnd[i] = inSegEnd[i];
        segMax[i] = inSegMax[i];
        segMin[i] = inSegMin[i];
        segStrike[i] = inSegStrike[i];
    }

    varMethod = 0;
    numPriceArrays = inNumPriceArrays;
    gridType = inGridType;

    AllocateArrays();
}


void FDEngine1F::AllocateArrays() {
    static const string method = "FDEngine1F::AllocateArrays";

    int i;

 //   sigma = new double [nStockSteps+1];

	sigma.resize(nStockSteps+1);
    
    stockArray = new double [nStockSteps+1];
    stockArrayNew = new double [nStockSteps+1];
    fwdGrid = new double [nStockSteps+1];

    range = TreeSliceGeneral::Range::create( 1, 0, nStockSteps );
    int dimBits = (1<<(*range)->nDim)-1;
    optionArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );
    optionOldArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );
    optionInterpArray = TreeSliceGeneralCont::create( *range, numPriceArrays, dimBits );

    upBarrierNew = new double [numPriceArrays];
	upPayoutNew = new double [numPriceArrays];
	upPayoutDeltaNew = new double [numPriceArrays];
	upBarrierOld = new double [numPriceArrays];
	upPayoutOld = new double [numPriceArrays];
	upPayoutDeltaOld = new double [numPriceArrays];
	downBarrierNew = new double [numPriceArrays];
	downPayoutNew = new double [numPriceArrays];
	downPayoutDeltaNew = new double [numPriceArrays];
	downBarrierOld = new double [numPriceArrays];
	downPayoutOld = new double [numPriceArrays];
	downPayoutDeltaOld = new double [numPriceArrays];

    for (i=0; i<numPriceArrays; i++) {
        upBarrierNew[i]=-1;
	    upPayoutNew[i]=0;
	    upPayoutDeltaNew[i]=0;
	    upBarrierOld[i]=-1;
	    upPayoutOld[i]=0;
	    upPayoutDeltaOld[i]=0;
	    downBarrierNew[i]=-1;
	    downPayoutNew[i]=0;
	    downPayoutDeltaNew[i]=0;
	    downBarrierOld[i]=-1;
	    downPayoutOld[i]=0;
	    downPayoutDeltaOld[i]=0;
    }
}

void FDEngine1F::loop(FDPayoff1F *engineCBs, double stockNow, double *price, double *divPert, double *irPert) {
    static const string method = "FDEngine1F::loop";
    int i, j, k;
    ZFDSolver1D fdSolver;
    double *tmpPointer;

//    double vol;
    double fwdNew;
    double fwdOld;
    double pvFact;
    int    currSeg;
    double delta, gamma;

    fdSolver.forceNonNegative(engineCBs->Positive());
    
    if (fdSolver.allocateMem(nStockSteps) == FAILURE) 
        throw ModelException(method, "fdSolver.allocateMem failure");

    currSeg = nSegments-1;
    if (useFwdGrid == false){
//        engineCBs->InitStockArrayFD(currSeg, stockArray);
        if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
            throw ModelException(method, "fdSolver.createGrid failure");
        tmpPointer = fdSolver.getSpotsPtr();
        for (j=0; j<=nStockSteps;j++){
            stockArray[j] = tmpPointer[j];
        }

    } else {
//        engineCBs->InitStockArrayFD(currSeg, fwdGrid);
        if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
            throw ModelException(method, "fdSolver.createGrid failure");
        fwdNew = forwards[nTimeSteps];
        tmpPointer = fdSolver.getSpotsPtr();
        for (j=0; j<=nStockSteps;j++){
            fwdGrid[j] = tmpPointer[j];
            stockArray[j] = fwdNew*tmpPointer[j];
        }
    }

    engineCBs->preCalcFD(nTimeSteps, 0, 0, numPriceArrays-1);
    engineCBs->PayoffAtMatFD(stockArray, nTimeSteps, 0, nStockSteps, 0, numPriceArrays-1, (*optionArray)[0]);
    
#ifdef FD1F_DEBUG_FILE1
        double time = 0.0;
        for (int kk=1; kk<=nTimeSteps;kk++)
            time += dt[kk];
#endif

    for (i = nTimeSteps; i>0; i--){

#ifdef FD1F_DEBUG_FILE1
        if (i==nTimeSteps){
            time = 0.0;        
            for (int kk=1; kk<=nTimeSteps;kk++)
                time += dt[kk];
        }
        else
            time -= dt[i+1];

        Debug_OutPutFD1FPrice(i, time, stockArray, optionArray, 
                              0, nStockSteps, numPriceArrays, FD1F_DEBUG_FILE1, i==nTimeSteps);
#endif

        swapT( optionOldArray, optionArray );

//        vol = sqrt((Variance[i] - Variance[i-1])/TimePts.TradeYrFrac[i]);
   //     vol = engineCBs->getVolFD(-1, i);
  //      for (j=0; j<=nStockSteps;j++){
   //         sigma[j] = vol;
   //     }


        
        if (useFwdGrid == false){
 //           throw ModelException(method, "useFwdGrid must be true for now");  
			

			engineCBs->getVolFD(i - 1, sigma, stockArray,0, nStockSteps);
            engineCBs->preCalcFD(i-1, 0, 0, numPriceArrays-1);
            
            for (k=0; k<numPriceArrays; k++) {
                if ( k == 0 ) {
                    fdSolver.forceNonNegative(engineCBs->Positive());
                } else {
                    fdSolver.forceNonNegative(true);
                }

                if (fdSolver.solveBarrier(nStockSteps,ir[i] + irPert[i],divy[i] + divPert[i],&*sigma.begin(),NULL,0.,dt[i],
                                          (*optionOldArray)[0][k],(*optionArray)[0][k],upBarrierNew[k],upPayoutNew[k],upPayoutDeltaNew[k],
                                          upBarrierOld[k],upPayoutOld[k],upPayoutDeltaOld[k],downBarrierNew[k],downPayoutNew[k],
                                          downPayoutDeltaNew[k],downBarrierOld[k],downPayoutOld[k],downPayoutDeltaOld[k],
                                          0,NULL,1.e-9,0,.5,0,varMethod) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }  
            }

            if (currSeg > 0 && i-1 == segEnd[currSeg-1]) {

                currSeg--;
             //   engineCBs->InitStockArrayFD(currSeg, stockArrayNew); 
                
             if (fdSolver.createGrid(segMin[currSeg],segMax[currSeg],nStockSteps,gridType,segStrike[currSeg]) == FAILURE) 
                    throw ModelException(method, "fdSolver.createGrid failure");
                
                tmpPointer = fdSolver.getSpotsPtr();
                for (j=0; j<=nStockSteps;j++){
                   stockArrayNew[j] = tmpPointer[j];
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
                                            (*optionInterpArray)[0][k]) == FAILURE) 
                    {
                        throw ModelException(method, "FDCubicSplineInterp failure");
                    } 
                }

                swapT( optionArray, optionInterpArray );
                
                tmpPointer = stockArray;
                stockArray = stockArrayNew;
                stockArrayNew = tmpPointer;

#ifdef FD1F_DEBUG_FILE1
        Debug_OutPutFD1FPrice(i, time, stockArray, optionArray, 
                              0, nStockSteps, numPriceArrays, FD1F_DEBUG_FILE1, i==nTimeSteps);
#endif

            }  

            engineCBs->PayoffBeforeMatFD(stockArray, i-1, 0, nStockSteps, 0, numPriceArrays-1, (*optionArray)[0]);

        } else {
            fwdNew = forwards[i-1];
            fwdOld = forwards[i];
            for (j=0; j<=nStockSteps;j++){
                stockArray[j] = fwdNew*fwdGrid[j];
            }
            pvFact = pvs[i];//DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);


			engineCBs->getVolFD(i - 1, sigma, stockArray,0, nStockSteps);

            engineCBs->preCalcFD(i-1, 0, 0, numPriceArrays-1);
            
            for (k=0; k<numPriceArrays; k++) {
                if ( k == 0 ) {
                    fdSolver.forceNonNegative(engineCBs->Positive());
                } else {
                    fdSolver.forceNonNegative(true);
                }

                if (fdSolver.solveBarrier(nStockSteps,irPert[i],divPert[i],&*sigma.begin(),NULL,0.,dt[i],
                                          (*optionOldArray)[0][k],(*optionArray)[0][k],upBarrierNew[k]/fwdNew,upPayoutNew[k],upPayoutDeltaNew[k]*fwdNew,
                                          upBarrierOld[k]/fwdOld,upPayoutOld[k],upPayoutDeltaOld[k]*fwdOld,downBarrierNew[k]/fwdNew,downPayoutNew[k],
                                          downPayoutDeltaNew[k]*fwdNew,downBarrierOld[k]/fwdOld,downPayoutOld[k],downPayoutDeltaOld[k]*fwdOld,
                                          0,NULL,1.e-9,0,.5,0,varMethod) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }


                for (j=0; j<=nStockSteps;j++){
                    (*optionArray)[0][k][j] *= pvFact;
                }
            }
            
            if (currSeg > 0 && i-1 == segEnd[currSeg-1]) {

                currSeg--;                
     //           engineCBs->InitStockArrayFD(currSeg, fwdGrid); 
     //           for (j=0; j<=nStockSteps;j++){
     //               stockArrayNew[j] = fwdNew*fwdGrid[j];
     //           }

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

                    // Include tolerance because forward prices (which stockArray[nStockSteps] can be calculated from)
                    // are computed with platform-dependent accuracy
                    if (Maths::areEqualWithinTol(stockArray[nStockSteps], upBarrierNew[k], DBL_THRESHHOLD)) {
                        nodeMax = nStockSteps - 1;
                    }
                    else
                    {
                        for (j=nStockSteps; j>=0; j--) {
                            if(stockArray[j] < upBarrierNew[k] ) {
                                nodeMax = j;
                                break;
                            }      
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

#ifdef FD1F_DEBUG_FILE1
        Debug_OutPutFD1FPrice(i, time, stockArray, optionArray, 
                              0, nStockSteps, numPriceArrays, FD1F_DEBUG_FILE1, i==nTimeSteps);
#endif

            }  

            engineCBs->PayoffBeforeMatFD(stockArray, i-1, 0, nStockSteps, 0, numPriceArrays-1, (*optionArray)[0]);
            
        }
    }   

    for (k=0; k<numPriceArrays; k++) {
        if ( stockNow > stockArray[nStockSteps] && upBarrierNew[k] >= 0) {
		    // Correct values for spots above the barrier
			price[k] = (*optionArray)[0][k][nStockSteps] + 
                       upPayoutDeltaNew[k]*(stockNow - stockArray[nStockSteps]);
			delta = upPayoutDeltaNew[k];
			gamma = 0.;
        } else {
            if (FDInterpolationD(nStockSteps+1,stockArray,(*optionArray)[0][k],1,&stockNow,price+k,&delta,&gamma) == FAILURE)
            // if (FDCubicSplineInterpD(nStockSteps+1,stockArray,(*optionArray)[0][k],1,&stockNow,price+k,&delta,&gamma) == FAILURE)
                throw ModelException(method, "FDCubicSplineInterpD failure");
        }
    }

#ifdef FD1F_DEBUG_FILE1
        Debug_OutPutFD1FPrice(i, 0, stockArray, optionArray, 
                              0, nStockSteps, numPriceArrays, FD1F_DEBUG_FILE1, i==nTimeSteps);
#endif

    return;
}


DRLIB_END_NAMESPACE

