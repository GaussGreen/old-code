/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treemethod.cpp
 *
 *  \brief object to implement the generic tree
 *	 first generation
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/treemethod.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/timer.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"		/// for matrix conversion

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"

//// gpnummethods
#include "gpnummethods/treetransition.h"

/// STL
#include <iomanip> /// for setprecision()
CC_USING_NS(std,dec);
CC_USING_NS(std,setw);
CC_USING_NS(std,fixed);
CC_USING_NS(std,setprecision);

#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

#define NB_PROBA    2
#define PROBA_UP    0
#define PROBA_DOWN  1

/// sqrt(3)
#define STD_MAPPING_FACTOR  1.73205080757

#define MAX_MAPPING_FACTOR  2.0

/// 2/sqrt(3)
#define MIN_MAPPING_FACTOR  1.15470053838

static double totalBackwardInductTime=0.0;

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_TreeMethod::CopyNoCleanUp(const ARM_TreeMethod& rhs)
{
    int i,nb=rhs.itsSpaceSteps.size();
    itsSpaceSteps.resize(nb);
    for(i=0;i<nb;++i)
        itsSpaceSteps[i] = rhs.itsSpaceSteps[i] ? (ARM_GP_TriangularMatrix*) rhs.itsSpaceSteps[i]->Clone() : NULL;

    nb=rhs.itsMaxIndex.size();
    itsMaxIndex.resize(nb);
    for(i=0;i<nb;++i)
        itsMaxIndex[i] = rhs.itsMaxIndex[i] ? (ARM_GridIndex*) rhs.itsMaxIndex[i]->Clone() : NULL;
    
    nb=rhs.itsNumMethodStateGlobalStdDevs.size();
    itsNumMethodStateGlobalStdDevs.resize(nb);
    for(i=0;i<nb;++i)
        itsNumMethodStateGlobalStdDevs[i] = rhs.itsNumMethodStateGlobalStdDevs[i] ? (ARM_GP_TriangularMatrix*) rhs.itsNumMethodStateGlobalStdDevs[i]->Clone() : NULL;
  
	nb=rhs.itsNumMethodStateGlobalVars.size();
    itsNumMethodStateGlobalVars.resize(nb);
    for(i=0;i<nb;++i)
        itsNumMethodStateGlobalVars[i] = rhs.itsNumMethodStateGlobalVars[i] ? (ARM_GP_Matrix*) rhs.itsNumMethodStateGlobalVars[i]->Clone() : NULL;
	
	itsStdDevRatio=rhs.itsStdDevRatio;
    itsMinStdDev=rhs.itsMinStdDev;
    itsNbMinSteps=rhs.itsNbMinSteps;

	nb=rhs.itsProbaSpotForward.size();
    itsProbaSpotForward.resize(nb);
    for(i=0;i<nb;++i)
        itsProbaSpotForward[i] = rhs.itsProbaSpotForward[i] ? (ARM_GP_Vector*) rhs.itsProbaSpotForward[i]->Clone(): NULL;

	nb=rhs.itsTransitionStates.size();
    itsTransitionStates.resize(nb);
    for(i=0;i<nb;++i)
	{
	 	ARM_TreeTransitions* transitionStates_i = rhs.itsTransitionStates[i];
		int nb_i = (*transitionStates_i).size();
		itsTransitionStates[i] = new vector<ARM_TreeTransition*>(nb_i);
		for(int j=0;j<nb_i;++j)
		    (*itsTransitionStates[i])[j] = (*transitionStates_i)[j] ? (ARM_TreeTransition*) ((*transitionStates_i)[j])->Clone(): NULL;
	}

    itsIsProbabilityComputation = rhs.itsIsProbabilityComputation;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_TreeMethod::CleanUp()
{
    size_t i;

    for(i=0;i<itsSpaceSteps.size();++i)
    {
        delete itsSpaceSteps[i];
        itsSpaceSteps[i]=NULL;
    }

    for(i=0;i<itsMaxIndex.size();++i)
    {
        delete itsMaxIndex[i];
        itsMaxIndex[i]=NULL;
    }   

    for(i=0;i<itsNumMethodStateGlobalStdDevs.size();++i)
    {
        delete itsNumMethodStateGlobalStdDevs[i];
        itsNumMethodStateGlobalStdDevs[i]=NULL;
    }

	for(i=0;i<itsNumMethodStateGlobalVars.size();++i)
    {
        delete itsNumMethodStateGlobalVars[i];
        itsNumMethodStateGlobalVars[i]=NULL;
    }
	
	for(i=0;i<itsProbaSpotForward.size();++i)
    {
        delete itsProbaSpotForward[i];
        itsProbaSpotForward[i]=NULL;
    }

	for(i=0;i<itsTransitionStates.size();++i)
    {
        if(itsTransitionStates[i] != NULL)
		{
			for(int j=0;j<(*itsTransitionStates[i]).size();++j)
			{
				delete (*itsTransitionStates[i])[j];
				(*itsTransitionStates[i])[j]=NULL;
			}
		}

		delete itsTransitionStates[i];
		itsTransitionStates[i] =NULL;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeMethod::ARM_TreeMethod(double stdDevRatio,double minStdDev,int nbMinStep,bool isProbabilityComputation)
:	itsStdDevRatio(stdDevRatio),itsMinStdDev(minStdDev),itsNbMinSteps(nbMinStep),itsIsProbabilityComputation(isProbabilityComputation)
{
	SetName(ARM_TREEMETHOD);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeMethod::ARM_TreeMethod(const ARM_TreeMethod& rhs)
: ARM_NumMethod(rhs)
{
    CopyNoCleanUp(rhs);
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_TreeMethod::~ARM_TreeMethod()
{
    CleanUp();
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: operator =
///	Returns: itself
///	Action : Affectation
////////////////////////////////////////////////////
ARM_TreeMethod& ARM_TreeMethod::operator = (const ARM_TreeMethod& rhs)
{
	if(this != &rhs)
	{
		ARM_NumMethod::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_TreeMethod
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_TreeMethod::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    
    os << "\n\n";
    os << indent << "ARM_TreeMethod\n";
    os << indent << "--------------\n";

    os << indent << "Pricing direction : " << GP_PricingDirectionTxt[(int)GetPricingDirection()] << "\n";
    os << indent << "Truncation Level : " << itsStdDevRatio << "\n";
    os << indent << "Min Annualised StdDev for space step resizing : " << itsMinStdDev << "\n";
    os << indent << "Required nb step before first event date : " << setw(4) << itsNbMinSteps << "\n";

    os << indent << "Nb Steps : " << GetNbSteps() << "\n";
    os << indent << "Diffusion Schedule :" << "\n";
    int i,j;

    if( GetTimeSteps() && GetTimeSteps()->size() ) // to get the pointor !!
    {
        os << indent << "Last Time Idx : " << GetLastTimeIdx() << "\n";
        os << indent << "  #     Time lag";
        if(!itsSpaceSteps.empty())
            os << "    Max Idx     Space Steps\n";
        for(i=0;i<GetTimeSteps()->size();++i)
        {
            os << indent << setw(4) << i << "  " << fixed << setprecision(2) << setw(7) << GetTimeStep(i);
            if(!itsSpaceSteps.empty())
            {
                os << "  I=( ";
                for(j=0;j<itsMaxIndex[i]->size();++j)
                    os << dec << setw(4) << (*(itsMaxIndex[i]))[j] << ",";
                os << ")  S=( ";
                for(j=0;j<itsSpaceSteps[i]->rows();++j)
                    os << fixed << setprecision(6) << setw(8) << (*(itsSpaceSteps[i]))(j,j) << ",";
                os << ")\n";
            }
        }
    }
    else
        os << indent << "not initialised !\n";

    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Clone
///	Returns: 
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_TreeMethod::Clone() const
{
	return new ARM_TreeMethod(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Init
///	Returns: 
///	Action : Initialiation of the tree
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeMethod::Init( ARM_PricingModel& model, double firstInductTime)
{
    // Restore the schedule initialised by the model
#ifdef __GP_STRICT_VALIDATION
    if((*GetTimeSteps())[0] != 0.0)
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "The model schedule is inconsistent");
#endif

    totalBackwardInductTime=0.0;

	ARM_Timer timer;
	timer.ClockStartTime();

    ComputeAndSetTimeSteps(model);
    int nbSteps = GetTimeSteps()->size();


    /// Last index setting for backward induction purpose
    SetLastTimeIdx(nbSteps-1);

	
	/// Free memory (due to object persistance)
    CleanUp();

	//ARM_GP_MatrixPtr relativeDrifts;
	//ARM_GP_MatrixPtr absoluteDrifts;
	//model.IntegratedLocalDrifts(*(GetTimeSteps()),relativeDrifts,absoluteDrifts);
	//SetNumMethodRelativDrift(relativeDrifts);
	//SetNumMethodAbsoluteDrift(absoluteDrifts);
	//model.ModelStateLocalVariancesAndStdDev(*(GetTimeSteps()));

	/// Call the model to compute NumMethodStates variances/covariances
	InitContinuousMapping(model);

	timer.ClockEndTime();
    double mappingInitTime = timer.GetDuration();

	timer.ClockStartTime();

	/// Construction of the tree (space steps, max grid Index, transion states, and probaforwards) 
	TreeConstruction(*GetTimeSteps(), itsProbaSpotForward, itsTransitionStates, model.NeedArrowDebreuPrices() || itsIsProbabilityComputation);

	timer.ClockEndTime();
    double forwardInductionTime = timer.GetDuration();

	model.PostInit();

    /// Compute initial NumMethodStates
    ARM_GridIndex maxIndex(GetMaxIndex(nbSteps-1));
    ARM_PricingStatesPtr initStates( model.FirstPricingStates( maxIndex.Range() ) );
	
    ComputeStateVariables(nbSteps-1,maxIndex,initStates);

	/// Compute initial ModelStates	
	model.TreeStatesToModelStates(initStates,nbSteps-1);

/***
FILE* f=fopen("c:\\temp\\dumpTree1G.txt","a");
for(size_t i=0;i<itsProbaSpotForward.size();++i)
{
    if(itsProbaSpotForward[i])
    {
        for(size_t j=0;j<itsProbaSpotForward[i]->size();++j)
            fprintf(f,"   (%3d,%3d) : %15.13lf\n",i,j,(*(itsProbaSpotForward[i]))[j]);
        fprintf(f,"\n");
    }
}
fclose(f);
***/

/***
FILE* f=fopen("c:\\temp\\dumpTree1G.txt","a");
fprintf(f,"mapping=%10.5lf\tforwardInduct=%10.5lf\n",mappingInitTime,forwardInductionTime);
fclose(f);
***/

    return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeCstTimeSchedule
///	Returns: 
///	Action : Initialiation of a schedule uniformly
///          spaced in time
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeCstTimeSchedule()
{
    const ARM_GP_Vector& modelSchedule = *GetTimeSteps();
    int nbModelSteps=modelSchedule.size()-1; // [0]=0

    int i,nbSteps = GetNbSteps();

    int nbToAdd=nbSteps-nbModelSteps;
    if(nbToAdd <= 0)
        return; /// keep the model schedule

    ARM_GP_Vector timeSteps(nbSteps+1);
    timeSteps[0]=0.0;
    int j,idx=1;

    int nbInSteps = nbToAdd/nbModelSteps;
    int nbExtraInSteps = nbToAdd - nbInSteps*nbModelSteps;
    int jmax;
    double lastModelTime=0.0,modelTime;
    double step,time;
    for(i=1;i<=nbModelSteps;++i)
    {
        modelTime = modelSchedule[i];
        jmax = (i-1<nbExtraInSteps ? nbInSteps+1 : nbInSteps);
        step = (modelTime-lastModelTime)/(jmax+1);
        for(j=0,time=lastModelTime+step;j<jmax;++j,time+=step,++idx)
            timeSteps[idx]=time;
        timeSteps[idx]=modelTime; // to avoid rounding pbs !
        ++idx;
        lastModelTime = modelTime;
    }

    SetTimeSteps(timeSteps);
 }

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeMixTimeVarSchedule
///	Returns: 
///	Action : Initialiation of a schedule uniformly
///          spaced in variance. It is based on the
///          variance sum of the state variables
///          of the model
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeMixTimeVarSchedule(const ARM_PricingModel& model)
{
    const ARM_GP_Vector& modelSchedule = *GetTimeSteps();
    int nbModelSteps=modelSchedule.size()-1; /// Actual model steps

    int nbSteps = GetNbSteps();
    if(nbSteps-nbModelSteps <= 0)
        return;

    ARM_MatrixVector localVCV;
    model.NumMethodStateLocalVariances(modelSchedule,localVCV);

    ARM_GP_Vector varIncr(localVCV.size());
    double varIncrSum=0.0;
    int i,j;
    for(i=0;i<localVCV.size();++i)
    {
        varIncr[i] = localVCV[i]->trace();
        varIncrSum += varIncr[i];
    }
    if(varIncrSum < itsMinStdDev*itsMinStdDev*modelSchedule[nbModelSteps]/K_YEAR_LEN)
    {
        /// Variance too low, switch to a time
        /// uniform spacing
        ComputeCstTimeSchedule();
        return;
    }

    
    /// 1st step : set steps before the 1st model time and
    /// fulfill the constraint on the minimum number
    int nbInSteps = (int)(floor(varIncr[0]/varIncrSum*nbSteps));
    if(nbInSteps<1)
        nbInSteps=1; // current event date
    double time=0.0;
    double lastModelTime=0.0,modelTime;
    double var = 0.0;
    int ifirst=1,nbAdded=0;
    double step;

    ARM_GP_Vector timeSteps(nbSteps+nbModelSteps);
    timeSteps[0]=0.0;
    int idx=1;

    ARM_IntVector modelTimeIdx(nbModelSteps+1);
    modelTimeIdx[0]=0;

    if(nbInSteps <= itsNbMinSteps)
    {
        /// Force itsNbMinSteps steps between 0 and
        /// the first model time

        var += varIncr[ifirst-1];
        modelTime = modelSchedule[ifirst];
        step = modelTime/itsNbMinSteps;
        for(;idx<itsNbMinSteps;++idx)
        {
            time += step;
            timeSteps[idx]=time;
        }
        time = modelTime;
        timeSteps[idx]=time;  // to avoid rounding pbs !
        modelTimeIdx[ifirst]=idx;
        ++idx;
        nbAdded = itsNbMinSteps;
        lastModelTime = modelTime;
        ifirst = 2;
    }


    /// 2nd step : insert steps between all model times
    for(i=ifirst;i<=nbModelSteps;++i)
    {
        var += varIncr[i-1] ;
        nbInSteps = (int)(floor(var/varIncrSum*nbSteps)) - nbAdded;
        if(nbInSteps<1)
            nbInSteps=1; // current event date
        if(nbAdded+nbInSteps > timeSteps.size())
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Inconsistency in the mix time/var tree schedule initialisation");
        }
        modelTime = modelSchedule[i];
        step = (modelTime-lastModelTime)/nbInSteps;
        for(j=0;j<nbInSteps-1;++j,++idx)
        {
            time += step;
            timeSteps[idx]=time;
        }
        time = modelTime;
        timeSteps[idx]=time;  // to avoid rounding pbs !
        modelTimeIdx[i]=idx;
        ++idx;
        nbAdded += nbInSteps;
        lastModelTime = modelTime;
    }

    /// 3rd step : review this time schedule using the
    /// local variance sum profile which is now
    /// strictly increasing (but not the variance profile)
    nbSteps=idx;
    ARM_GP_Vector X(nbSteps),Y(nbSteps);
    Y[0]=0.0;
    for(j=0;j<nbSteps;++j)
        X[j]=timeSteps[j];

	/// management of simplematrix to triangular matrix
    ARM_MatrixVector lastLocalVCV;
	model.NumMethodStateLocalVariances(X,lastLocalVCV);
    
    for(j=1;j<nbSteps;++j)
        Y[j] = Y[j-1] + lastLocalVCV[j-1]->trace();

    ARM_GP_Vector finalTimeSteps(nbSteps);
    finalTimeSteps[0]=0.0;
    idx=1;
    var=0.0;
    int nextIdx,lastIdx=0;
    double nextVar,varStep;
    for(int nextModelIdx=1;nextModelIdx<=nbModelSteps;++nextModelIdx)
    {
        nextIdx=modelTimeIdx[nextModelIdx];
        nextVar=Y[nextIdx];
        nbSteps=nextIdx-lastIdx;
        varStep=(nextVar-var)/nbSteps;
        var += varStep;
        for(int i=0;i<nbSteps-1;++i)
        {
            time=YtoX(X,Y,var,lastIdx,nextIdx);
            finalTimeSteps[idx]=time;
            var += varStep;
            ++idx;
        }
        if(idx != nextIdx)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Inconsistency in the final rescheduling in constant variance");
        }
        time=modelSchedule[nextModelIdx];
        finalTimeSteps[idx]=time;
        lastIdx=nextIdx;
        var=nextVar;
        ++idx;
    }

    SetTimeSteps(finalTimeSteps);

    /// Free memory
	DeletePointorVector<ARM_GP_Matrix>( localVCV );
	DeletePointorVector<ARM_GP_Matrix>( lastLocalVCV );

}


double ARM_TreeMethod::YtoX(ARM_GP_Vector& X,ARM_GP_Vector& Y,double y,int& lastIdx,int nextIdx)
{
    /// Search x such that Y(x)= y, Y() is linearly interpolated
    /// Y is strictly increasing and we have Y[lastIdx] < y <= Y[nextIdx]

    for(int idx=lastIdx+1;idx<=nextIdx;++idx)
    {
        if(y<=Y[idx])
            break;
    }

    lastIdx=idx-1; /// for next call by increasing y

    double xinf=X[lastIdx];
    double yinf=Y[lastIdx];
    double xsup=X[idx];
    double ysup=Y[idx];

    return xinf + (y-yinf)*(xsup-xinf)/(ysup-yinf);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeCstVarSchedule
///	Returns: 
///	Action : Initialiation of a schedule uniformly
///          spaced in variance. It is based on the
///          variance sum of the state variables
///          of the model
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeCstVarSchedule(const ARM_PricingModel& model)
{
    const ARM_GP_Vector& modelSchedule = *GetTimeSteps();
    int nbModelSteps=modelSchedule.size()-1; /// Actual model steps

    int nbSteps = GetNbSteps();
    int nbToAdd=nbSteps-(nbModelSteps-1); // -1 because last model steps will be the last step
    if(nbToAdd <= 0)
        return;

	/// management of simplematrix to triangular matrix
    ARM_MatrixVector localVCV;
	ARM_MatrixVector globalVCV;
   
    model.NumMethodStateLocalGlobalVariances(modelSchedule,localVCV,globalVCV);
	double firstModelTime=modelSchedule[1];
    double firstVar = globalVCV[1]->trace();
    double lastModelTime=modelSchedule[nbModelSteps];
    double lastVar = globalVCV[nbModelSteps]->trace();

    if(lastVar < itsMinStdDev*itsMinStdDev*lastModelTime/K_YEAR_LEN)
    {
        /// Variance too low, switch to a time
        /// uniform spacing
        ComputeCstTimeSchedule();
        return;
    }


    ARM_GP_Vector timeSteps(nbSteps+1); // +1 for [0]=0
    timeSteps[0]=0.0;
    int idx=1;

    double time,timeErr;
    double varErrRel=0.0001;

    double var,prevVar=0.0;

    double varStep=lastVar/nbToAdd;
    int nbStepBefore=(int)(floor(firstVar/varStep));
    time=model.VarianceToTime(nbStepBefore*varStep,0.0,firstModelTime);

    if(time < firstModelTime - varStep/10.0)
    {
        /// Don't merge additional step
        nbStepBefore++;
    }

    ARM_IntVector modelTimeIdx(nbModelSteps+1);
    modelTimeIdx[0]=0;

    /// 1st stage : set steps before the 1st model time with
    /// fulfilling the constraint on the minimum number
    if(nbStepBefore <= itsNbMinSteps)
    {
        /// Force itsNbMinSteps steps between 0 and
        /// the first model time

        double varStepBefore=firstVar/itsNbMinSteps;
        var=varStepBefore;
        time=model.VarianceToTime(var,0.0,firstModelTime);
        for(idx=1;idx<itsNbMinSteps;++idx)
        {
            timeSteps[idx]=time;
            var+=varStepBefore;
            time=model.VarianceToTime(var,time,firstModelTime);
        }
        timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-time;
        if(time > firstModelTime + timeErr)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Inconsistency before first model time in the tree schedule initialisation");
        }

        varStep=(lastVar-firstVar)/(nbToAdd-itsNbMinSteps+1);
    }
    else if(nbStepBefore > 1) 
    {
        /// Keep the same spacing in variance before and
        /// after the first model time
        var=varStep;
        time=model.VarianceToTime(var,0.0,firstModelTime);
        for(idx=1;idx<nbStepBefore;++idx)
        {
            timeSteps[idx]=time;
            prevVar=var;
            var+=varStep;
            time=model.VarianceToTime(var,time,firstModelTime);
        }
    }

    /// Next time step to add is the first model time step
    time=firstModelTime;


    /// 2nd stage : insert steps between all model times
    /// using the global constant variance step
    double nextVar=firstVar;
    double nextModelTime=firstModelTime;
    int nextModelIdx=1;
    timeErr=model.VarianceToTime(firstVar*(1+varErrRel),firstModelTime,firstModelTime+1.0)-
        firstModelTime;
    for(;idx<=nbSteps;++idx)
    {
        if(nextModelTime - timeErr <= time)
        {
            // After next model time then stop just on it
            timeSteps[idx]=nextModelTime;
            modelTimeIdx[nextModelIdx]=idx;

            if(nextModelIdx < nbModelSteps)
            {
                if(time <= nextModelTime + timeErr)
                {
                    // Merge with the current model time in fact
                    var = nextVar + varStep;
                    time=model.VarianceToTime(var,time,modelSchedule[nextModelIdx+1]);
                    timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-
                        time;
                }

                // Next model time is the new target to reach
                nextModelIdx++;
                prevVar=nextVar;
                nextModelTime=modelSchedule[nextModelIdx];
                nextVar=globalVCV[nextModelIdx]->trace();
            }
            else
            {
                // Last model time is reached
                break;
            }
        }
        else
        {
            // Strictly before next model time => a step is added
            timeSteps[idx]=time;
            prevVar=var;
            var+=varStep;
            time=model.VarianceToTime(var,time,nextModelTime);
            timeErr=model.VarianceToTime(var*(1+varErrRel),time,time+1.0)-time;
        }
    } // idx <= nbSteps


    /// 3rd stage : reschedule evenly in variance between each
    /// model time with the number of time steps computed in 2nd stage
    ARM_GP_Vector finalTimeSteps(idx+1);
    int lastIdx=0;
    var=0.0;
    time=0.0;
    finalTimeSteps[lastIdx]=time;
    idx=1;
    for(nextModelIdx=1;nextModelIdx<=nbModelSteps;++nextModelIdx)
    {
        nextVar=globalVCV[nextModelIdx]->trace();
        nbSteps=modelTimeIdx[nextModelIdx];
        nbSteps -= lastIdx;
        varStep=(nextVar-var)/nbSteps;
        nextModelTime=modelSchedule[nextModelIdx];
        var += varStep;
        for(int i=0;i<nbSteps-1;++i)
        {
            time=model.VarianceToTime(var,time,nextModelTime);
            finalTimeSteps[idx]=time;
            var += varStep;
            ++idx;
        }
        if(idx != modelTimeIdx[nextModelIdx])
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Inconsistency in the final rescheduling in constant variance");
        }
        time=nextModelTime;
        finalTimeSteps[idx]=time;
        lastIdx=idx;
        var=nextVar;
        ++idx;
    }

    SetTimeSteps(finalTimeSteps);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: InitContinousMapping
///	Returns: 
///	Action : Initialise the continuous mapping of
///          the state variables. Done from local
///          & global variances/covariances matrix
///          of state variables and using a Cholesky
///          transformation if necessary
///          Local variances are saved before
////////////////////////////////////////////////////
void ARM_TreeMethod::InitContinuousMapping( const ARM_PricingModel& model )
{
	ARM_MatrixVector tmpLocalVar;
	ARM_MatrixVector tmpStateLocalMappings;
	ARM_MatrixVector tmpStateMappings;

	/// computes in one go the global and local std dev and variances!
    /// step [i] at ti  => local variance ti -> ti+1
    ///                 => variance 0 -> ti (nothing to do at t0=0)

    /// Free memory
    DeletePointorVector<ARM_GP_Matrix>( itsNumMethodStateGlobalVars );
    itsNumMethodStateGlobalVars.resize(0);
    DeletePointorVector<ARM_GP_TriangularMatrix>( itsNumMethodStateGlobalStdDevs );
    itsNumMethodStateGlobalStdDevs.resize(0);

	model.NumMethodStateLocalGlobalVariancesAndStdDev(*GetTimeSteps(),tmpLocalVar,itsNumMethodStateGlobalVars,
		tmpStateLocalMappings,tmpStateMappings);

	//SetNumMethodStateLocalStdDevs( ConvertToTriangularMatrixVector( tmpStateLocalMappings ) );
	//itsNumMethodStateGlobalStdDevs	= ConvertToTriangularMatrixVector( tmpStateMappings );
	//SetNumMethodStateLocalVars(tmpLocalVar);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: InitDiscreteMapping
///	Returns: 
///	Action : Initialise the discrete mapping of
///          the state variables using the continuous
///          mapping of them
////////////////////////////////////////////////////
void ARM_TreeMethod::InitDiscreteMapping()
{
    int i,j=0;
    int nbSteps=GetTimeSteps()->size();
    itsSpaceSteps.resize(nbSteps);

    double mappingFactor,minTimeVar,localVar;
    double minVar = itsMinStdDev*itsMinStdDev/K_YEAR_LEN;

    ARM_GP_Matrix* stateLocalVar=NULL;

    size_t nbDir=0;// = GetNumMethodStateLocalStdDevs()[0]->GetRowsNb();

    itsSpaceSteps[0] = new ARM_GP_TriangularMatrix(nbDir,0.0);

    /// Space step [i] at ti, compute from local variance ti-1 -> ti
    for(int timeIdx=1;timeIdx<nbSteps;++timeIdx)
    {
        itsSpaceSteps[timeIdx] = new ARM_GP_TriangularMatrix(nbDir);
        //stateLocalVar=GetNumMethodStateLocalVars()[timeIdx-1];

        /// The mimimum value of space steps is computed to fufill
        /// the minimum local stdDev requirement
        for(i=0;i<nbDir;++i)
        {
            minTimeVar = minVar * (GetTimeStep(timeIdx)-GetTimeStep(timeIdx-1));
            localVar = (*stateLocalVar)(i,i);
            if(localVar < minTimeVar)
            {
                if(localVar < K_DOUBLE_TOL)
                {
                    throw Exception(__LINE__, __FILE__, ERR_TOL_PB,
                        "Can't resize space steps, local variance too low");
                }
                mappingFactor = STD_MAPPING_FACTOR * sqrt(minTimeVar / localVar);
            }
            else
                mappingFactor = STD_MAPPING_FACTOR;

            //for(j=0;j<=i;++j)
            //    (*(itsSpaceSteps[timeIdx]))(i,j) = mappingFactor * (*(GetNumMethodStateLocalStdDevs()[timeIdx-1]))(i,j);
        }
    }
}
////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeMaxIndex
///	Returns: 
///	Action : Compute the maximum value of the tree
///          index for each time step
///          At the moment, no truncation is  used
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeMaxIndex()
{
    int i,j,k,timeIdx,nbTimeStep=GetTimeSteps()->size();
    itsMaxIndex.resize(nbTimeStep);

    int nbDir=itsSpaceSteps[0]->GetRowsNb();

    double x,y,z,rz,sum;


    /// First step : compute truncatation indexes from
    /// the end to the begining of the tree
    ARM_IntVector truncIndex(nbDir);
    vector< ARM_IntVector* > truncIndexes(nbTimeStep);
    int index;
    ARM_GP_TriangularMatrix* stateMapping;
    ARM_GP_TriangularMatrix* spaceStep;
    for(timeIdx=nbTimeStep-1;timeIdx>=1;--timeIdx)
    {
        spaceStep=itsSpaceSteps[timeIdx];

        /// Inverse global continuous mapping to get
        /// Zi = Bi1.X1 + ... + Bii.Xi
        /// itsStateMappings is definitely inverted !
        stateMapping=itsNumMethodStateGlobalStdDevs[timeIdx];
        stateMapping->Inverse();


        /// Compute maximum index w.r.t. truncation strategy
        /// The tree limit is assumed to be a rectangle
        /// rather than the actual ellipsoid
        for(i=0;i<nbDir;++i)
        {
            x=0.0;
            for(k=i-1;k>=0;--k)
            {
                sum=0.0;
                for(j=k;j<=i;++j)
                    sum += (*stateMapping)(i,j) * (*spaceStep)(j,k);
                x += truncIndex[k] * (sum>0 ? sum : -sum);
            }
            y = (itsStdDevRatio + x) / (*stateMapping)(i,i);
            x = y / (*spaceStep)(i,i);
            truncIndex[i] = (int)(floor(x)+1);
        }

        truncIndexes[timeIdx]= new ARM_IntVector(truncIndex);
    }


    /// Second step : compute tree limit taking
    /// into account truncation indexes
    ARM_GP_TriangularMatrix* curDiscMapping = itsSpaceSteps[0];
    ARM_GP_TriangularMatrix* nextDiscMapping;

    ARM_GP_MatrixPtr stateLocalDrift;//=GetNumMethodRelativeDrifts();
    ARM_GP_Vector nextExpectedState(nbDir);

    int edgeIdx,edgePos,nbEdge = 1<<nbDir;
    ARM_GridIndexVector curStates(nbEdge);
    for(i=0;i<nbEdge;++i)
        curStates[i] = new ARM_GridIndex(nbDir);

    ARM_IntVector edgeNextIndex(nbDir);
    ARM_IntVector* truncNextIndex;

    ARM_GridIndex maxState(nbDir);

    double curDiscValue,nextDiscValue,contExpectedValue;

    bool isTrunc;

    for(timeIdx=0;timeIdx<nbTimeStep-1;++timeIdx)
    {
        itsMaxIndex[timeIdx] = new ARM_GridIndex(maxState);

        nextDiscMapping=itsSpaceSteps[timeIdx+1];

        truncNextIndex = truncIndexes[timeIdx+1];

        /// For each edge, compute next states
        for(edgeIdx=0;edgeIdx<nbEdge;++edgeIdx)
        {
            edgePos=edgeIdx;
            for(i=0;i<nbDir;++i)
            {
                curDiscValue=0.0;
                nextDiscValue=0.0;
                for(j=0;j<i;++j)
                {
                    curDiscValue += (*(curStates[edgeIdx]))[j] * (*curDiscMapping)(i,j);
                    nextDiscValue += nextExpectedState[j] * (*nextDiscMapping)(i,j);
                }
                curDiscValue += (*(curStates[edgeIdx]))[i] * (*curDiscMapping)(i,i);

                contExpectedValue = curDiscValue * (*stateLocalDrift)(timeIdx,i);

                z = (contExpectedValue - nextDiscValue)/(*nextDiscMapping)(i,i);

                rz = ROUND(z);

                index = (int)rz;

                isTrunc = index >= (*truncNextIndex)[i] || index <= - (*truncNextIndex)[i];

                if(isTrunc)
                {
                    /// Truncate to the nearest expected index
                    edgeNextIndex[i]=index;
                    if( (rz>=z && z>0.0) || (rz<=z && z<0.0) )
                        /// Continuous expectation will be fulfilled
                        nextExpectedState[i]=z;
                    else
                        /// Continuous expectation will NOT be fulfilled
                        nextExpectedState[i]=rz;
                }
                else
                {
                    /// No truncation
                    edgeNextIndex[i]=index + (edgePos&1 ? 1 : -1);
                    nextExpectedState[i]=z;
                }

                edgePos = edgePos>>1;
            }
            curStates[edgeIdx]->SetIndex(edgeNextIndex);
        }
        maxState.SetToMax(curStates);
        curDiscMapping=nextDiscMapping;

        /// Update current new tree limits
        for(edgeIdx=0;edgeIdx<nbEdge;++edgeIdx)
        {
            edgePos=edgeIdx;
            for(i=0;i<nbDir;++i)
            {
                (*(curStates[edgeIdx]))[i] = (edgePos&1 ? maxState[i] : -maxState[i]);
                edgePos = edgePos>>1;
            }
        }
    }

    itsMaxIndex[nbTimeStep-1] = new ARM_GridIndex(maxState);


    /// Free memory
    for(i=0;i<nbEdge;++i)
        delete curStates[i];

    for(i=1;i<nbTimeStep;++i)
        delete truncIndexes[i];

}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeBinomialTransition
///	Returns: 
///	Action : Compute a binomial transition by only
///          fitting a drift error.
///          The variance and the probabilities of
///          the transition are computed
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeBinomialTransition(int dirIdx,
                                               double driftErr,
                                               ARM_GP_Vector& proba,
                                               ARM_GP_Vector& nextStateVariance) const
{
    /// driftErr is the error between the continuous index (z)
    /// and the integer index (rz=round(z))
    if(driftErr < 0.0)
    {
        /// Possible transitions are (rz,pMid) & (rz-1,pDown)
        /// (Top edge if truncated or internal forced binomial transition)
        proba[PROBA_DOWN] = -driftErr;
        proba[PROBA_UP] = 0.0;
        /// then pMid = 1 - pDown

        nextStateVariance[dirIdx]=proba[PROBA_DOWN] * (1.0 - proba[PROBA_DOWN]);
    }
    else
    {
        /// Possible transitions are (rz+1,pUp) & (rz,pMid)
        /// (Bottom edge if truncated or internal forced binomial transition)
        proba[PROBA_UP] = driftErr;
        proba[PROBA_DOWN] = 0.0;
        /// then pMid = 1 - pUp

        nextStateVariance[dirIdx]=proba[PROBA_UP] * (1.0 - proba[PROBA_UP]);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeTransitionStates
///	Returns: A tranisition state list for the grid
///          located at the timeIdx
///	Action : Compute all connected states on the
///          next slice from a current one
////////////////////////////////////////////////////
ARM_TreeTransitions* ARM_TreeMethod::ComputeTransitionStates(
        int timeIdx) const
{
    /// Get datas of the current and next slices
    ARM_GP_TriangularMatrix* curDiscMapping=itsSpaceSteps[timeIdx];
    ARM_GP_TriangularMatrix* nextDiscMapping=itsSpaceSteps[timeIdx+1];

    ARM_GP_MatrixPtr stateLocalDrift;// = GetNumMethodRelativeDrifts();
    ARM_GP_TriangularMatrix* stateLocalMapping=NULL;// = const_cast<ARM_TreeMethod*>(this)->GetNumMethodStateLocalStdDevs()[timeIdx];
    ARM_GP_Matrix* stateLocalVar=NULL;//=const_cast<ARM_TreeMethod*>(this)->GetNumMethodStateLocalVars()[timeIdx];

    ARM_GridIndex curMaxState(GetMaxIndex(timeIdx));
    ARM_GridIndex curState(curMaxState);

    ARM_GridIndex nextMaxState(GetMaxIndex(timeIdx+1));
    ARM_GridIndex nextState(nextMaxState);

    int i,j,nbDir=nextState.size();
    ARM_GP_Vector nextExpectedState(nbDir);
    ARM_GP_Vector nextStateVariance(nbDir);

    /// Test if next space steps were resized beyong
    /// admissible bounds due to min stdDev requirement
    /// The flag below will force the variance correction
    /// in further proba computation
    double mappingFactor;
    ARM_GP_Vector invSpaceii(nbDir);
    bool forceVarCorrect=false;
    for(i=0;i<nbDir;++i)
    {
        invSpaceii[i] = 1.0 / (*nextDiscMapping)(i,i);
        if(!forceVarCorrect)
        {
            mappingFactor = (*nextDiscMapping)(i,i) / (*stateLocalMapping)(i,i);
            if(mappingFactor < MIN_MAPPING_FACTOR || mappingFactor > MAX_MAPPING_FACTOR)
                forceVarCorrect=true;
        }
    }

    ARM_VectorVector nextProba(nbDir);
    for(i=0;i<nbDir;++i)
        nextProba[i]=new ARM_GP_Vector(NB_PROBA,0.0);

    ARM_GP_Vector* proba;

    double curDiscValue,nextDiscValue,contExpectedValue;
    double nextDiscVar;
    double x,z,rz,driftErr,driftErr2,contVar;

    /// Create the transition state vector
    /// from the current slice to the next one
    ARM_TreeTransitions* newStates = new ARM_TreeTransitions(curState.Range());

    bool isVarCorrectAllowed,isVarCorrect,isTrunc,isOut;
    bool isTrinomial;

    /// Compute transition states for each state of the current slice
    int stateIdx=0;
    curState.Reset();
    for(;curState <= curMaxState;++curState)
    {
        isVarCorrectAllowed=true;       // to allow variance correction
        isVarCorrect=forceVarCorrect;   // to correct variance if truncation needed

        /// For each direction compute...
        for(i=0;i<nbDir;++i)
        {
            /// Approximated expected state in the next slice index counting
            curDiscValue=0.0;
            nextDiscValue=0.0;
            nextDiscVar=0.0;
            for(j=0;j<i;++j)
            {
                curDiscValue += curState[j] * (*curDiscMapping)(i,j);
                nextDiscValue += nextExpectedState[j] * (*nextDiscMapping)(i,j);
                if(isVarCorrect)
                    nextDiscVar += nextStateVariance[j] * (*nextDiscMapping)(i,j) * (*nextDiscMapping)(i,j);
            }
            curDiscValue += curState[i] * (*curDiscMapping)(i,i);

            contExpectedValue = curDiscValue * (*stateLocalDrift)(timeIdx,i);

            z = (contExpectedValue - nextDiscValue)*invSpaceii[i];

            rz = ROUND(z);

            /// Tree transition will be forced inside the range
            /// [-nextMaxIndex[i],nextMaxIndex[i]]
            nextState[i]=(int)(rz);

            proba = nextProba[i];

            isTrunc = (nextState[i] >= nextMaxState[i] || nextState[i] <= - nextMaxState[i]);

            if(isTrunc)
            {
                isOut = nextState[i] > nextMaxState[i] || nextState[i] < - nextMaxState[i];

                if(!isOut && ((rz>=z && z>0.0) || (rz<=z && z<0.0)) )
                {
                    /// Allow variance correction for further directions
                    if(isVarCorrectAllowed)
                        isVarCorrect=true;

                    /// Continuous expectation is still fulfilled but not variance...
                    nextExpectedState[i] = z;
                    driftErr = z - rz;

                    ///...because only both nearest nodes are connected through
                    /// a simple binomial transition
                    ComputeBinomialTransition(i,driftErr,*proba,nextStateVariance);

                } // if binomial transition

                else
                {
                    /// No more variance correction because the brute force truncation
                    /// may distrub probability computation
                    isVarCorrectAllowed=false;
                    isVarCorrect=false;

                    if(isOut)
                    {
                        // Limit the transition to max indexes
                        nextState[i] = (nextState[i] > nextMaxState[i] ? nextMaxState[i] : - nextMaxState[i]);
                        rz = (double) nextState[i];
                    }

                    /// Continuous expectation is no more fulfilled...
                    nextExpectedState[i] = rz;

                    ///...because the nearest central node is only connected (pMid=1) 
                    (*proba)[PROBA_UP] = 0.0;
                    (*proba)[PROBA_DOWN] = 0.0;
                    nextStateVariance[i] = 0.0;

                } // if straight transition

            } // if truncation

            else
            {
                /// Theoretical index expectation
                nextExpectedState[i] = z;

                /// Associated transition probabilities to fit
                /// residual drift error and variance
                driftErr = z - rz;
                driftErr2 = driftErr*driftErr;

                isTrinomial=true;
                if(isVarCorrect)
                {
                    /// Compute residual variance to be fitted
                    /// by the current direction...
                    contVar = ((*stateLocalVar)(i,i) - nextDiscVar)*invSpaceii[i]*invSpaceii[i];
                    if(contVar < 0.0)
                        /// Switch to a binomial transition
                        isTrinomial=false;
                    else if(contVar + driftErr2 > 1.0)
                        ///...but not too much else negative probas appear
                        contVar = 1.0 - driftErr2;
                }
                else if(isVarCorrectAllowed)
                {
                    /// No truncation on previous directions => only variance on
                    /// current direction is fitted (on others automatically fulfilled)
                    contVar = (*stateLocalMapping)(i,i) * invSpaceii[i];
                    contVar *= contVar;
                }
                else
                    isTrinomial=false;


                if(isTrinomial)
                {
                    x=contVar + driftErr2;

                    /// pUp
                    (*proba)[PROBA_UP]=0.5*(x + driftErr);

                    /// pDown and pMid=1-pUp-pDown...not saved !
                    (*proba)[PROBA_DOWN]=(*proba)[PROBA_UP] - driftErr;

                    if(x<0.0 || x>1.0 || (*proba)[PROBA_UP]<0.0 || (*proba)[PROBA_UP]>1.0 ||
                        (*proba)[PROBA_DOWN]<0.0 || (*proba)[PROBA_DOWN]>1.0)
                        /// Switch to a binomial transition
                        isTrinomial=false;
                    else
                        /// Actual index variance (may differ from the continuous value
                        /// if truncation occurs on previous directions)
                        nextStateVariance[i] = contVar;
                }

                if(!isTrinomial)
                {
                    /// Be less ambitious, index limitation is quite tricky
                    ComputeBinomialTransition(i,driftErr,*proba,nextStateVariance);
                }

            } // if no truncation

            x=(*proba)[PROBA_UP]+(*proba)[PROBA_DOWN];
            if(x<0.0 || x>1.0 || (*proba)[PROBA_UP]<0.0 || (*proba)[PROBA_UP]>1.0 ||
                (*proba)[PROBA_DOWN]<0.0 || (*proba)[PROBA_DOWN]>1.0)
            {
                throw Exception(__LINE__, __FILE__, ERR_PROPAGATION_PB,
                    "Negative probabilities in the tree");
            }

        } // loop nbDir


        /// Set a new transition state from current to next slice
        nextState.UpdatePosition();
        (*newStates)[stateIdx]=new ARM_TreeTransition(nextState,nextProba);
        ++stateIdx;
    }

    /// Free memory
    for(i=0;i<nbDir;++i)
        delete nextProba[i];

    return newStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeStateProba
///	Returns: Value of the transition probability
///	Action : Compute the probability of a state
///          referred by a relative index (to a
///          reference index)
///          Standard version (not optimised if
///          number of directions exceeds 2) 
////////////////////////////////////////////////////
double ARM_TreeMethod::ComputeStateProba(
        ARM_GridIndex& relIndex,
        const ARM_VectorVector& elemProbas) const
{
    double stateProba;
    ARM_GP_Vector* dirProba=elemProbas[0];
    if(relIndex[0]==-1)
        stateProba = (*dirProba)[PROBA_DOWN];
    else if(relIndex[0]==0)
        stateProba = 1.0-(*dirProba)[PROBA_UP]-(*dirProba)[PROBA_DOWN];
    else
        stateProba = (*dirProba)[PROBA_UP];

    if(stateProba > 0.0)
    {
        double p;
        for(int i=1;i<relIndex.size();++i)
        {
            dirProba=elemProbas[i];
            if(relIndex[i]==-1)
                p = (*dirProba)[PROBA_DOWN];
            else if(relIndex[i]==0)
                p = 1.0-(*dirProba)[PROBA_UP]-(*dirProba)[PROBA_DOWN];
            else
                p = (*dirProba)[PROBA_UP];

            if(p > 0.0)
                stateProba *= p;
            else
                return 0.0;
        }
    }
    else
        return 0.0;

    return stateProba;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeStateProba
///	Returns: Value of the transition probability
///	Action : Compute the probability of a state
///          referred by a relative index (to a
///          reference index)
///          Optimised version using datas from
///          the previous call : index and 1D
///          associated probabilities
////////////////////////////////////////////////////
double ARM_TreeMethod::ComputeStateProba(
        ARM_GridIndex& relIndex,
        const ARM_VectorVector& elemProbas, 
        ARM_GridIndex& lastIndex,
        ARM_GP_Vector& lastProbas) const
{
    int imax=relIndex.size()-1;

    double stateProba=1.0;
    int i=imax;
    if(relIndex[imax] == lastIndex[imax])
    {
        /// Proba already set, just use last proba
        // and loop back until indexes differs
        stateProba=lastProbas[imax];
        for(i=imax-1;i>=0;--i)
        {
            if(relIndex[i] == lastIndex[i])
                stateProba *= lastProbas[i];
            else
                break;
        }
    }

    /// Update probas & last index/proba
    ARM_GP_Vector* dirProba;
    for(;i>=0;--i)
    {
        dirProba=elemProbas[i];
        if(relIndex[i]==-1)
            lastProbas[i] = (*dirProba)[PROBA_DOWN];
        else if(relIndex[i]==0)
            lastProbas[i] = 1.0-(*dirProba)[PROBA_UP]-(*dirProba)[PROBA_DOWN];
        else
            lastProbas[i] = (*dirProba)[PROBA_UP];

        stateProba *= lastProbas[i];

        lastIndex[i] = relIndex[i];
    }

    return stateProba;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeStateVariables
///	Returns: List of state variables
///	Action : Compute the state variables for a
///          given position in the tree (defined
///          by its timeIdx)
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeStateVariables(
        int timeIdx,
        const ARM_GridIndex& maxIndex,
        ARM_PricingStatesPtr& states) const
{
    int i,j,nbDir=maxIndex.size();
    const ARM_GP_TriangularMatrix& discMapping = *(itsSpaceSteps[timeIdx]);

    ARM_GridIndex curIndex(maxIndex);
    double value;
    for(i=0;i<nbDir;++i)
    {
        for(curIndex.Reset();curIndex <= maxIndex;++curIndex)
        {
            value=0.0;
            for(j=0;j<=i;++j)
                value += curIndex[j] * discMapping(i,j);
            states->SetNumMethodState(curIndex.GetPosition(),i,value);
        }
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: Induct
///	Returns: List of pricing states
///	Action : backward induction from last time
///          to the next one
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeMethod::Induct(
        const ARM_PricingModel& model,
		ARM_PricingStatesPtr& states,
		double toTime)
{
	// Obsolete method 

    return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeExercise
///	Returns: An exercise boundary
///	Action : Compute an exercise boundary
///          of an exercise node
////////////////////////////////////////////////////
ARM_ExerciseBoundary* ARM_TreeMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: InductOneStep
///	Returns: List of payoffs
///	Action : backward induction from timeIdx+1
///          to timeIdx
////////////////////////////////////////////////////
void ARM_TreeMethod::InductOneStep(
        ARM_GP_Matrix& payoffs,
		ARM_TreeTransitions* transitionStates,
		int timeIdx)
{      
    
    int nbDir=itsMaxIndex[timeIdx]->size();
	int nbStates=transitionStates->size();
	int nbPayoffs=payoffs.GetColsNb();

	ARM_GP_Matrix nextPayoffs=payoffs;
	
	ARM_GP_Matrix currPayoffs(nbStates,nbPayoffs);
	
    
    ARM_GridIndex maxState(nbDir),curState(nbDir);
    ARM_GridIndex localRefState(nbDir);
    ARM_UnsymGridIndex relState(nbDir),relMinState(nbDir),relMaxState(nbDir),relLastState(nbDir);

    ARM_GP_Vector payoff(nbPayoffs);


    size_t curPosition,nextPosition;
    double stateProba;
	ARM_TreeTransition* transitionState;
	int i;
       
    maxState.SetIndex(GetMaxIndex(timeIdx));
    curState.SetIndex(maxState);  // better than...
    curState.Reset();        // ...curState = -maxState

    for(;curState<=maxState;++curState)
    {
        /// Get the transition of the current state
        curPosition=curState.GetPosition();
        transitionState=(*transitionStates)[curPosition];


        /// Get the central state of the transition...
        localRefState=transitionState->GetIndex();


        ///...and the upper right & lower left with
        /// truncation on the edge if necessary
        relMaxState.RelativeIndex(localRefState,1);
        relMinState.RelativeIndex(localRefState,-1);
        relState.SetRangeIndex(relMinState,relMaxState); // a reset is done


        /// Loop over connected states of the next slice
        for(i=0;i<nbPayoffs;++i)
            payoff[i]=0.0;
        relLastState.SetIndex(relMaxState);
        while(relState <= relMaxState)
        {
            nextPosition=relState.AbsolutePosition(localRefState);


            stateProba = ComputeStateProba(relState,transitionState->GetProbas());

            /// Compute conditionnal expectation for each payoff
            if(stateProba > 0.0)
            {
                for(i=0;i<nbPayoffs;++i)
                    payoff[i] += stateProba  * nextPayoffs(nextPosition,i);
            }
            ++relState;
        }

        /// Create a new pricing state (only payoffs are usefull at this stage)
        for(i=0;i<nbPayoffs;++i)
            currPayoffs(curPosition,i) = payoff[i];
    }

	nextPayoffs = currPayoffs;
	
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: InductForwardOneStep
///	Returns: List of payoffs
///	Action : backward induction from timeIdx+1
///          to timeIdx
////////////////////////////////////////////////////
void ARM_TreeMethod::InductForwardOneStep(
        ARM_GP_Vector& payoffs,
		ARM_TreeTransitions* transitionStates,
		int timeIdx)
{      
    
    int nbDir=itsMaxIndex[timeIdx]->size();
	int nbStates=transitionStates->size();
	
	int nbNextStates = itsMaxIndex[timeIdx+1]->Range();

	ARM_GP_Vector nextPayoffs(nbNextStates);
	
	ARM_GP_Vector currPayoffs = payoffs;
	
    
    ARM_GridIndex maxState(nbDir),curState(nbDir);
    ARM_GridIndex localRefState(nbDir);
    ARM_UnsymGridIndex relState(nbDir),relMinState(nbDir),relMaxState(nbDir),relLastState(nbDir);

    double payoff;


    size_t curPosition,nextPosition;
    double stateProba;
	
    maxState.SetIndex(GetMaxIndex(timeIdx));
    curState.SetIndex(maxState);  // better than...
    curState.Reset();        // ...curState = -maxState

    for(;curState<=maxState;++curState)
    {
        /// Get the transition of the current state
        curPosition=curState.GetPosition();

        /// Get the central state of the transition...
        localRefState=(*transitionStates)[curPosition]->GetIndex();


        ///...and the upper right & lower left with
        /// truncation on the edge if necessary
        relMaxState.RelativeIndex(localRefState,1);
        relMinState.RelativeIndex(localRefState,-1);
        relState.SetRangeIndex(relMinState,relMaxState); // a reset is done


        /// Loop over connected states of the next slice
        payoff=0.0;
        relLastState.SetIndex(relMaxState);
        while(relState <= relMaxState)
        {
            nextPosition=relState.AbsolutePosition(localRefState);


            stateProba = ComputeStateProba(relState,(*transitionStates)[curPosition]->GetProbas());

            /// Compute conditionnal expectation for each payoff
            if(stateProba > 0.0)
            {
               nextPayoffs.Elt(nextPosition) += stateProba  * currPayoffs.Elt(curPosition);
            }
            ++relState;
        }
      
    }

	payoffs = nextPayoffs;
	
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: InductForwardOneStep
///	Returns: List of payoffs
///	Action : backward induction from timeIdx+1
///          to timeIdx
////////////////////////////////////////////////////

void ARM_TreeMethod::InductForwardFromToNextStep(
        ARM_GP_Vector& payoffs,
		int timeIdx_1,
		int timeIdx_2)
{      
	int nbCurrStates = itsMaxIndex[timeIdx_1]->Range();
	int nbPayoffs=payoffs.size();

    if(nbCurrStates != nbPayoffs)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Inconsistency in Fwdpropagation ");
    }
	
	for (int timeIdx = timeIdx_1; timeIdx <timeIdx_2; timeIdx++)
	{
		ARM_TreeTransitions* transitionStates=ComputeTransitionStates(timeIdx);
		InductForwardOneStep(payoffs,transitionStates,timeIdx);	
		for(int i=0;i<transitionStates->size();++i)
            delete (*transitionStates)[i];
        delete transitionStates;
	}

}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ProbaSpotForward
///	Returns: List of payoffs
///	Action : 
////////////////////////////////////////////////////
void ARM_TreeMethod::TreeConstruction(const ARM_GP_Vector& timeSteps, vector< ARM_GP_Vector* >& probaSoptFwd, 
									  vector< ARM_TreeTransitions* >& transitionStates, bool needTreeProbaForward)
	 
{
	clock_t startTreeConstruction = clock();

	/// Compute the discrete mapping (space steps)
    InitDiscreteMapping();

	clock_t initDiscreteMapping = clock();

    /// Compute maximum value of tree indexes (tree shape)
    ComputeMaxIndex();

	clock_t computeMaxIndex = clock();

	//TransitionStates and Probaforward
	int lastTimeIdx=GetLastTimeIdx();

	if(needTreeProbaForward)
	{
		ARM_GP_Vector payoffs(1,1.0);	
		transitionStates.resize(lastTimeIdx);
		probaSoptFwd.resize(lastTimeIdx);
		for(int timeIdx=0;timeIdx<lastTimeIdx;++timeIdx)
		{
			transitionStates[timeIdx] = ComputeTransitionStates(timeIdx);
			InductForwardOneStep(payoffs,transitionStates[timeIdx],timeIdx);
			probaSoptFwd[timeIdx] = new ARM_GP_Vector(payoffs);

		}
	}
	else
	{
		transitionStates.resize(lastTimeIdx);
		for(int timeIdx=0;timeIdx<lastTimeIdx;++timeIdx)
		{
			transitionStates[timeIdx] = ComputeTransitionStates(timeIdx);
		}
	}
	
	clock_t endTreeConstruction = clock();

	CC_Ostringstream os;

	os.flush();
	os << "Init Discrete mapping = " << (double)(initDiscreteMapping-startTreeConstruction) / ((double)CLOCKS_PER_SEC) << std::endl;
	ARM_TheEventViewer.Instance()->AddToMessage(os.str());

	os.flush();
	os << "Compute max index = " << (double)(computeMaxIndex-initDiscreteMapping) / ((double)CLOCKS_PER_SEC) << std::endl;
	ARM_TheEventViewer.Instance()->AddToMessage(os.str());

	os.flush();
	os << "Compute transition states & probas = " << (double)(endTreeConstruction-computeMaxIndex) / ((double)CLOCKS_PER_SEC) << std::endl;
	ARM_TheEventViewer.Instance()->AddToMessage(os.str());
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ReInit
///	Returns: 
///	Action : reinitialisation for loop pricing! not supported
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeMethod::ReInit( const ARM_PricingModel& model)
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"tree method does not support multiple loops!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change (in case of pricing 
/// direction change)
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeTimeSteps
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_TreeMethod::ComputeAndSetTimeSteps( const ARM_PricingModel& model) const
{
	/// Switch for testing purpose !
    bool isCstTime=false;

    if(isCstTime)
    {
        /// For testing purpose, keep a schedule
        /// uniformly spaced in time
        const_cast< ARM_TreeMethod* >(this)->ComputeCstTimeSchedule();
    }
    else
    {
        /// Built a schedule based on time and variance considerations
        //ComputeCstVarSchedule(model);
        const_cast< ARM_TreeMethod* >(this)->ComputeMixTimeVarSchedule(model);

    }
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_TreeMethod::GetBuckets() const
{
	/// one bucket with size 1
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,1)); 
}

size_t ARM_TreeMethod::GetBucketIndex() const
{
	/// first bucket only
	return 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_TreeMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoop;
}



	
////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_TreeMethod::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    size_t i,j,nbStates=0,eventIdx=0,probaIdx;
    const ARM_GP_Vector& sched = *(GetTimeSteps());

    /// No probabilities available
    if(! &sched || !itsIsProbabilityComputation || itsProbaSpotForward.size()==0)
        return ARM_GP_MatrixPtr(NULL);

    if(eventTimes.size()==1 && eventTimes[0]==0.0)
        return ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,1,1.0));

    /// Time selection in schedule
    /// Take care : probabilities start at time step #1 (because #0 = 0 then proba=1)
    vector<int> timeIdx(eventTimes.size());
    for(i=0;i<timeIdx.size();++i) timeIdx[i]=-1;
    for(i=1;i<sched.size() && eventIdx < eventTimes.size();++i)
    {
        if(eventTimes[eventIdx] <= sched[i])
        {
            probaIdx = i-1;
            if(eventTimes[eventIdx] == sched[i] && itsProbaSpotForward[probaIdx])
            {
                timeIdx[eventIdx]=i;
                nbStates = (nbStates < itsProbaSpotForward[probaIdx]->size() ? itsProbaSpotForward[probaIdx]->size() : nbStates);
            }

            ++eventIdx;
        }
    }

    /// Save spot probabilities
    /// 1st line = event time
    ARM_GP_Matrix* probas = new ARM_GP_Matrix(nbStates+1,timeIdx.size());
    for(i=0;i<timeIdx.size();++i)
    {
        if(timeIdx[i]!=-1)
        {
            probaIdx = timeIdx[i]-1;
            (*probas)(0,i)=sched[timeIdx[i]];
            for(j=0;j<itsProbaSpotForward[probaIdx]->size();++j)
                (*probas)(j+1,i)=(*(itsProbaSpotForward[probaIdx]))[j];
        }
        else
        {
            /// Slice not found or probabilities not available
            (*probas)(0,i)=eventTimes[i];
            j=0;
        }

        for(;j<nbStates;++j)
            (*probas)(j+1,i)=0.0;
    }

    return ARM_GP_MatrixPtr(probas);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

