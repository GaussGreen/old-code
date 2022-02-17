//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetStochGarfSolver.cpp
//
//   Description : mixing of one factor finite difference algorithm
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1DRetStochGarfSolver.hpp"
//#include "edginc/UtilFuncs.hpp"

DRLIB_BEGIN_NAMESPACE
class FD1DRetSolver;

//-----------------------------------------------------------------------------

/** FD1DRetStochGarfSolver algorithm class that supports FD1DRetStochGarfSolver 
FD1DRetStochGarfSolver::FD1DRetStochGarfSolver(vector<FD1DRet* > engines) : engines(engines){
    numOfEngines = engines.size();
}
*/

/** FD1DRetStochGarfSolver algorithm class that supports FD1DRetStochGarfSolver */
FD1DRetStochGarfSolver::FD1DRetStochGarfSolver(FD1DRetStochGarf*  engineStochGarf) : FD1DRetSolver(engineStochGarf->FD1DRetLVs[0].get()), engineStochGarf(engineStochGarf){
    numOfEngines = engineStochGarf->FD1DRetLVs.size();
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

FD1DRetStochGarfSolver::~FD1DRetStochGarfSolver(){
}

/*
// solve backward or forward induction through all steps
// return price
void FD1DRetStochGarfSolver::roll(){

    static const string method = "FD1DRetStochGarfSolver::Roll";

    // starting step, only backward for now, need update if adding forward
    int rollDirection ;  
    int i;

    // keep a copy to be restored
    //FD1DRetStochGarf* orgEngine = engine->get();
    
    if (engines[0]->isFwdInduction == false){ //backward
        rollDirection =-1;
    }else{
        rollDirection = 1;
    }

    int currStep = (rollDirection == 1 ? 0 : engines[0]->timeLine->NumOfStep);
    // init segmenent index, 
    int  currSeg = engines[0]->getfdSeg(currStep);


    //aaaaaaaaaaaaaaaa
    //how to add this for barrier
    //it's not right
    //engine->prod->preCalcFD1DRet(currStep,  priceStart, priceEnd); 

    for(i=0;i<numOfEngines;i++){
        engines[i]->prod->preCalc(currStep);
    }

    // update spot first
    for (i = 0;i<numOfEngines;i++){
        engines[i]->updateSolverInfo(currStep, getSliceIdx(currStep), rollDirection);
    }

    if (engine->isFwdInduction ){// fwd induction
        // calcuate payoff at maturity
        int firststep = 0;
        for (i=0;i<numOfEngines;i++){
            engines[i]->prod->update(firststep, FDProduct::FWD_0);
        }
    }
    else{
        // calcuate payoff at maturity
        for (i=0;i<numOfEngines;i++){
            engines[i]->prod->update(engine->timeLine->NumOfStep, FDProduct::BWD_T);
        }
    }

    // move a step to start sweeping
    currStep += rollDirection;

    // sweep the fd
    for (; currStep >=0 && currStep <= engines[0]->timeLine->NumOfStep; 
        currStep += rollDirection){
            
        for (i = 0;i<numOfEngines;i++){
            engines[i].updateSolverInfo(currStep, getSliceIdx(currStep), rollDirection);
            engines[i].rollOneStep(currStep, rollDirection);
        }

        // now, mixing the results before evalute product payoff
        // to do....

        for (i = 0;i<numOfEngines;i++){
            // update product info
            const FDProductArray & products = engines[i]->getProducts();
            int step = currStep;
            for(int kk=0; kk < (int)(*products).size(); kk++){
                if (engineToUse->isFwdInduction == false){ //backward
                    (*products)[kk]->update(step, FDProduct::BWD);           
                }
                else{
                    (*products)[kk]->update(step, FDProduct::FWD);
                }
            }
            engines[i].postUpdateSolverInfo(currStep, getSliceIdx(currStep));
        }


        for (i = 0;i<numOfEngines;i++){
            engines[i].postUpdateSolverInfo(currStep, getSliceIdx(currStep));
            engines[i].switchSegment(currStep, getSliceIdx(currStep), currSeg, rollDirection);
        }
    }
}
*/
//-----------------------------------------------------------------------------

// solve backward or forward induction through all steps
// return price
void FD1DRetStochGarfSolver::roll(){

    static const string method = "FD1DRetStochGarfSolver::Roll";

    // starting step, only backward for now, need update if adding forward
    int rollDirection ;  
    int i;

    // keep a copy to be restored
    FD1DRetLV* baseEngine = engineStochGarf->FD1DRetLVs[0].get();
    
    if (baseEngine->isFwdInduction == false){ //backward
        rollDirection =-1;
    }else{
        rollDirection = 1;
    }

    int currStep = (rollDirection == 1 ? 0 : baseEngine->timeLine->NumOfStep);
    // init segmenent index, 
    int  currSeg = getfdSeg(currStep);


    //aaaaaaaaaaaaaaaa
    //how to add this for barrier
    //it's not right
    //engine->prod->preCalcFD1DRet(currStep,  priceStart, priceEnd); 

    for(i=0;i<numOfEngines;i++){
//        engineStochGarf->FD1DRetLVs[i]->prod->preCalc(currStep);
        engine=engineStochGarf->FD1DRetLVs[i].get();
        engine->prod->preCalc(currStep); 
    }

    // update spot first
    for (i = 0;i<numOfEngines;i++){
        engine = engineStochGarf->FD1DRetLVs[i].get();
        updateSolverInfo(currStep, 0/*engine->getSliceIndex(currStep)*/, rollDirection);
    }

    if (baseEngine->isFwdInduction ){// fwd induction
        // calcuate payoff at maturity
        int firststep = 0;
        for (i=0;i<numOfEngines;i++){
            engine=engineStochGarf->FD1DRetLVs[i].get();
            engine->prod->update(firststep, FDProduct::FWD_0);
        }
    }
    else{
        // calcuate payoff at maturity
        for (i=0;i<numOfEngines;i++){
            engine=engineStochGarf->FD1DRetLVs[i].get();
            engine->prod->update(engine->timeLine->NumOfStep, FDProduct::BWD_T);
        }
    }

    // move a step to start sweeping
    currStep += rollDirection;

    // sweep the fd
    for (; currStep >=0 && currStep <= baseEngine->timeLine->NumOfStep; 
        currStep += rollDirection){
            
        for (i = 0;i<numOfEngines;i++){
            engine = engineStochGarf->FD1DRetLVs[i].get();
            updateSolverInfo(currStep, 0/*engine->getSliceIndex(currStep)*/, rollDirection);
            rollOneStep(currStep, rollDirection);
        }

        // now, mixing the results before evalute product payoff
        // to do....

//        for (i = 0;i<numOfEngines;i++){
//            engine = engineStochGarf->FD1DRetLVs[i].get();
//            // update product info
//            const FDProductArray & products = engine->getProducts();
//            int step = currStep;
//            for(int kk=0; kk < (int)(*products).size(); kk++){
//                if (engine->isFwdInduction == false){ //backward
//                    (*products)[kk]->update(step, FDProduct::BWD);           
//                }
//                else{
//                    (*products)[kk]->update(step, FDProduct::FWD);
//                }
//            }
//            postUpdateSolverInfo(currStep, getSliceIdx(currStep));
//        }


        for (i = 0;i<numOfEngines;i++){
            engine = engineStochGarf->FD1DRetLVs[i].get();
            //postUpdateSolverInfo(currStep, getSliceIdx(currStep));
            postUpdateSolverInfo(currStep);
            switchSegment(currStep, 0/*engine->getSliceIndex(currStep)*/, currSeg, rollDirection);
        }
    }

    engine = engineStochGarf->FD1DRetLVs[1].get();  // use middle as MV's results.
    //engineStochGarf = dynamic_cast<FD1DRetStochGarf*>(engine);                    // and copy it as Main Results
    engineStochGarf->prod = engine->prod;
}


//-----------------------------------------------------------------------------    

DRLIB_END_NAMESPACE
