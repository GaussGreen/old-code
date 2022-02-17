/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file slice.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/slice.h"

/// gpbase
#include "gpbase/utilityport.h"  /// for CC_Round
#include "gpbase/timer.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

/// gpnumlib
#include "gpnumlib/newtonraphson.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/treeindex.h"
#include "gpnummethods/transitor.h"

///  define for code clarity
#define FIRST_STATE_VARIABLE 0


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_SliceBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_SliceBase
///	Routine: ARM_SliceBase
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////
ARM_SliceBase::ARM_SliceBase(const ARM_SliceBase& rhs)
:   ARM_RootObject(rhs),
    itsIndex(rhs.itsIndex)
{
    
    ARM_GP_VectorPtr null(NULL);
	itsSpotProbas = (rhs.itsSpotProbas != null
                    ? ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsSpotProbas->Clone()))
                    : null);
	itsArrowDebreuPrices =  (rhs.itsArrowDebreuPrices != null
                            ? ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsArrowDebreuPrices->Clone()))
                            : null);
}


////////////////////////////////////////////////////
///	Class  : ARM_SliceBase
///	Routine: ToSlice????
///	Returns: 
///	Action : Downcast to derivated object
////////////////////////////////////////////////////
ARM_Slice1DBase* ARM_SliceBase::ToSlice1DBase()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1DBase" );
}

const ARM_Slice1DBase* ARM_SliceBase::ToSlice1DBase() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1DBase" );
}

ARM_Slice1D* ARM_SliceBase::ToSlice1D()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1D" );
}

const ARM_Slice1D* ARM_SliceBase::ToSlice1D() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1D" );
}

ARM_Slice1DCstSymProba* ARM_SliceBase::ToSlice1DCstSymProba()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1DCstSymProba" );
}

const ARM_Slice1DCstSymProba* ARM_SliceBase::ToSlice1DCstSymProba() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Slice1DCstSymProba" );
}

ARM_SliceNDBase* ARM_SliceBase::ToSliceNDBase()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceNDBase" );
}

const ARM_SliceNDBase* ARM_SliceBase::ToSliceNDBase() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceNDBase" );
}

ARM_SliceND* ARM_SliceBase::ToSliceND()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceND" );
}

const ARM_SliceND* ARM_SliceBase::ToSliceND() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceND" );
}

ARM_SliceNDCstSymProba* ARM_SliceBase::ToSliceNDCstSymProba()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceNDCstSymProba" );
}

const ARM_SliceNDCstSymProba* ARM_SliceBase::ToSliceNDCstSymProba() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SliceNDCstSymProba" );
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Slice1DFunct
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : Slice1DFunct
///	Routine: Constructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
Slice1DFunct::Slice1DFunct( size_t timeIdx, double dt, const ARM_PricingModel* model, const ARM_GP_VectorPtr& xStates, const ARM_GP_VectorPtr& adPrices)
:   
    itsTimeIdx(timeIdx),
    itsDt(dt),
	itsModel(model),
    itsXStates(xStates), 
	itsADPrices(adPrices),
    itsStates(new ARM_PricingStates(xStates->size(),1,0,1))
{}

////////////////////////////////////////////////////
///	Class  : Slice1DFunct
///	Routine: operator()
///	Returns: double
///	Action : Compute the functor value
////////////////////////////////////////////////////
double Slice1DFunct::operator() (double x) const
{
    size_t i,nbStates = itsXStates->size();
    for(i=0;i<nbStates;++i)
        itsStates->SetNumMethodState(i,FIRST_STATE_VARIABLE,(*itsXStates)[i]+x);
	itsModel->TreeStatesToModelStates(itsStates,itsTimeIdx);

    ARM_GP_VectorPtr localDf = itsModel->LocalDiscounts(itsTimeIdx,itsDt,itsStates);
    double fx=0.0;
    for(i=0;i<nbStates;++i)
        fx += (*itsADPrices)[i] * (*localDf)[i];
    return fx;
}

////////////////////////////////////////////////////
///	Class  : Slice1DFunct
///	Routine: GetLocalDf
///	Returns: ARM_GP_VectorPtr
///	Action : Compute local discount associated to a
///        : deterministic shift x
////////////////////////////////////////////////////
ARM_GP_VectorPtr Slice1DFunct::GetLocalDf(double x) const
{
    size_t i,nbStates = itsXStates->size();
    for(i=0;i<nbStates;++i)
        itsStates->SetNumMethodState(i,FIRST_STATE_VARIABLE,(*itsXStates)[i]+x);
	itsModel->TreeStatesToModelStates(itsStates,itsTimeIdx);
    ARM_GP_VectorPtr localDf = itsModel->LocalDiscounts(itsTimeIdx,itsDt,itsStates);
    return localDf;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Slice1DBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Slice1DBase
///	Routine: ARM_Slice1DBase
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////
ARM_Slice1DBase::ARM_Slice1DBase(const ARM_Slice1DBase& rhs)
:   ARM_SliceBase(rhs),
    itsMin(rhs.itsMin), itsMax(rhs.itsMax), itsSpaceStep(rhs.itsSpaceStep),
    itsDriftCorrection(rhs.itsDriftCorrection)
{
    ARM_GP_VectorPtr null(NULL);
	itsXStates =    (rhs.itsXStates != null
                    ? ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsXStates->Clone()))
                    : null);
}

////////////////////////////////////////////////////
///	Class  : ARM_Slice1DBase
///	Routine: ComputeDriftCorrection
///	Returns: 
///	Action : Compute the correction to be added
///          to states to fit the discount factor
///          from 0 to next slice time
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_Slice1DBase::ComputeDriftCorrection(const ARM_PricingModel* model, ARM_SamplerBase* sampler, double nextTime, size_t driftOffset, double targetPayoff, bool isLocalDf)
{
    /// Get target discount price
    size_t timeIdx = GetIndex();
    double time = model->GetNumMethod()->GetTimeStep(timeIdx);

    double nrTol = 1.0e-8;;

    /// Generate X states if necessary
    if(itsXStates == ARM_GP_VectorPtr(NULL))
    {
        /// Not previously computed then compute Z states...
        ARM_GP_VectorPtr zStates = GetZStates();

        /// ... and restore X states through the sampler
        const ARM_Sampler1DBase* sampler1D = sampler->ToSampler1DBase();
        itsXStates = sampler1D->ComputeZtoXStates(timeIdx,zStates);
    }

    /// Compute the initial guess
    Slice1DFunct slice1DFct(timeIdx,nextTime-time,model,itsXStates,GetArrowDebreuPrices());
    UnaryFuncWithNumDerivative<double,double> fctToSolve( slice1DFct );
    double df0 = fctToSolve(0.0);
    double ddf0 = (*(fctToSolve.Derivative()))(0.0);
    double x0 = targetPayoff - df0;
    double driftCorrection=0.0;
    if(x0 < -nrTol || x0 > nrTol)
    {
        if(fabs(ddf0)<K_NEW_DOUBLE_TOL)
            ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME +
            " : can't find an initial solution to compute drift correction for current slice" );
        x0 /= ddf0;

        /// Create a NR solver then find the solution
		const double xtol = 1.0e-8;
		const size_t max_iter = 30;
        T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double,double> > solver(fctToSolve,targetPayoff,nrTol,xtol,max_iter);
		
		/// To initialize departure point
		solver.setInitialGuess(x0);

        driftCorrection=solver.Solve();
    }

    itsDriftCorrection = driftCorrection;

    if(isLocalDf)
        /// Compute local DF for further Arrow-Debreu updating
        return slice1DFct.GetLocalDf(driftCorrection);
    else
        return ARM_GP_VectorPtr(NULL);
}

////////////////////////////////////////////////////
///	Class  : ARM_Slice1DBase
///	Routine: UpdateProbasAndArrowDebreu
///	Returns: 
///	Action : Compute and set probabilities to reach
///          states of the slice. Idem for Arrow-Debreu
///          prices
////////////////////////////////////////////////////
void ARM_Slice1DBase::UpdateProbasAndArrowDebreu( bool needProbas, bool needADPrices, bool isLocalDf, ARM_GP_VectorPtr localDf, ARM_Slice1DBase* nextSlice1D)
{
    size_t stateIdx,nbStates = size();
    size_t nextNbStates = nextSlice1D->size();

    ARM_GP_Vector* spotProbas=NULL;
    if(needProbas) spotProbas = new ARM_GP_Vector(nextNbStates,0.0);
    ARM_GP_Vector* arrowDebreuPrices=NULL;
    if(needADPrices) arrowDebreuPrices = new ARM_GP_Vector(nextNbStates,0.0);
    double probaUp,probaDown,probaMid,proba,ADPrice;
    int nextStateIndex;
    for( stateIdx=0; stateIdx<nbStates; ++stateIdx )
    {
        nextStateIndex = GetNextNodeIndex(stateIdx) - nextSlice1D->GetMin();

        probaUp     = GetProbaUp(stateIdx);
        probaDown   = GetProbaDown(stateIdx);
        probaMid    = 1-probaDown-probaUp;

        if(needProbas)
        {
            proba = GetSpotProbas(stateIdx);
            if(probaUp>0)
				(*spotProbas)[nextStateIndex+1] += proba*probaUp;
            if(probaDown>0)
				(*spotProbas)[nextStateIndex-1] += proba*probaDown;
            (*spotProbas)[nextStateIndex] += proba*probaMid;
        }

        if(needADPrices)
        {
            ADPrice = (isLocalDf ? GetArrowDebreuPrices(stateIdx)*(*localDf)[stateIdx] :GetArrowDebreuPrices(stateIdx));
            if(probaUp>0)
				(*arrowDebreuPrices)[nextStateIndex+1] += ADPrice*probaUp;
            if(probaDown>0)
				(*arrowDebreuPrices)[nextStateIndex-1] += ADPrice*probaDown;
            (*arrowDebreuPrices)[nextStateIndex] += ADPrice*probaMid;
        }
    }

    nextSlice1D->SetSpotProbas(ARM_GP_VectorPtr(spotProbas));
    nextSlice1D->SetArrowDebreuPrices(ARM_GP_VectorPtr(arrowDebreuPrices));

#if defined( __GP_STRICT_VALIDATION)
    /// Check probabilities of next slice (AD prices without localDf <=> probabilities)
    if(needProbas || (needADPrices && !isLocalDf))
    {
        double probaSum=0.0;
        if(needProbas)
        {
            for( stateIdx=0; stateIdx<nextNbStates; ++stateIdx )
                probaSum += (*spotProbas)[stateIdx];
        }
        else
        {
            for( stateIdx=0; stateIdx<nextNbStates; ++stateIdx )
                probaSum += (*arrowDebreuPrices)[stateIdx];
        }

        if(probaSum - 1 < - K_NEW_DOUBLE_TOL || probaSum - 1 > K_NEW_DOUBLE_TOL)
            ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME + " : proba sum is not equal to 1 in the current slice" );
    }
#endif
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Slice1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Slice1D
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_Slice1D& ARM_Slice1D::operator = (const ARM_Slice1D& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_Slice1D();
		new (this) ARM_Slice1D (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_Slice1D
///	Routine: LinkedWithNextSlice
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Slice1D::LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff)
{
    const ARM_PricingModel* model = sampler->GetModel();

    ARM_Slice1DBase* nextSlice1D    = nextSlice->ToSlice1DBase();
    ARM_Truncator1D* truncator1D    = truncator->ToTruncator1D();
    ARM_Sampler1DBase* sampler1D    = sampler->ToSampler1DBase();
    size_t i,nbStates = size();

    size_t timeIdx = GetIndex();
    double sliceYf = model->GetNumMethod()->GetTimeStep(timeIdx)/K_YEAR_LEN;

    /// Get variable states on the current slice
    ARM_GP_VectorPtr zStates = GetZStates();

    /// Get conditional mean & variance and save X states
    ARM_GP_VectorPtr xStates;
    ARM_GP_VectorPtr mean   = sampler1D->ComputeConditionalMean( timeIdx,zStates,xStates );
    SetXStates(xStates);
    double var              = sampler1D->GetLocalVar( timeIdx );

    double meanError,meanRel,varRel;
    double spaceStepSquare  = nextSlice1D->GetSpaceStep() * nextSlice1D->GetSpaceStep();

    itsProbaUp.resize(nbStates);
    itsProbaDown.resize(nbStates);
    itsNextNodeIndex.resize(nbStates);

    /// Connection with possible truncation
    int minIndex=nextSlice1D->GetMax();
    int maxIndex=nextSlice1D->GetMin();
    for( i=0; i<nbStates; ++i )
    {
        meanRel             = (*mean)[i]/nextSlice1D->GetSpaceStep();
        itsNextNodeIndex[i] = CC_Round(meanRel);
        varRel              = var/spaceStepSquare;
        meanError           = meanRel-itsNextNodeIndex[i];
        itsProbaUp[i]       = 0.5 * (varRel+meanError*meanError+meanError);
        itsProbaDown[i]     = itsProbaUp[i]-meanError;

        /// Truncation if necessary
        truncator1D->UpdateNode(sliceYf,this,i,nextSlice,meanRel,varRel);

        /// Update index range
        if( itsProbaUp[i] > 0)
            maxIndex = CC_Max( maxIndex, itsNextNodeIndex[i]+1);
        else
            maxIndex = CC_Max( maxIndex, itsNextNodeIndex[i]);
        if( itsProbaDown[i] > 0)
            minIndex = CC_Min( minIndex, itsNextNodeIndex[i]-1);
        else
            minIndex = CC_Min( minIndex, itsNextNodeIndex[i]);
    }

    nextSlice1D->SetMin(minIndex);
    nextSlice1D->SetMax(maxIndex);

    /// Calibrate tree if necessary
    bool needADPrices = model->NeedArrowDebreuPrices();
    bool isLocalDf = model->NeedLocalDiscount();
    ARM_GP_VectorPtr localDf(NULL);
    if(isLocalDf && needADPrices)
    {
        double nextSliceTime = model->GetNumMethod()->GetTimeStep(nextSlice->GetIndex());
        localDf = ComputeDriftCorrection(model,sampler,nextSliceTime,driftOffset,targetPayoff);
    }


    /// Update probabilities and Arrow-Debreu prices
    if(needProbas || needADPrices)
        UpdateProbasAndArrowDebreu(needProbas,needADPrices,isLocalDf,localDf,nextSlice1D);


/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
if(needADPrices || needProbas)
{
    for( i=0; i<nbStates; ++i )
    {
        fprintf(f,"   (%3d,%3d) : ",timeIdx,i);
        if(needProbas)
            fprintf(f,"p=%15.13lf\t",GetSpotProbas(i));
        if(needADPrices)
            fprintf(f,"AD=%15.13lf\n",GetArrowDebreuPrices(i));
    }
    fprintf(f,"\n");
}
fclose(f);
****/

/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
if(needProbas)
{
    /// Check mean & variance of current slice (on mean reverting process X)
    ARM_GP_VectorPtr spotProbas = GetSpotProbas();
    ARM_GP_VectorPtr zStates = GetZStates();
    ARM_GP_VectorPtr xStates = sampler1D->ComputeStates(timeIdx,zStates);
    double treeMeanZ=0.0,treeVarZ=0.0,z;
    double treeMeanX=0.0,treeVarX=0.0,x;
    for( i=0; i<nbStates; ++i )
    {
        z = (*zStates)[i];
        treeMeanZ += z * (*spotProbas)[i];
        treeVarZ += z*z* (*spotProbas)[i];
        x = (*xStates)[i];
        treeMeanX += x * (*spotProbas)[i];
        treeVarX += x*x* (*spotProbas)[i];
    }
    treeVarZ -= treeMeanZ*treeMeanZ;
    treeVarX -= treeMeanX*treeMeanX;
    ARM_MatrixVector localVariances,variances;
    model->NumMethodStateLocalGlobalVariances(*(model->GetNumMethod()->GetTimeSteps()),localVariances,variances);
    double theoVar=(*(variances[timeIdx]))(0,0);
    DeletePointorVector<ARM_GP_Matrix>( localVariances );
    DeletePointorVector<ARM_GP_Matrix>( variances );

    double time = model->GetNumMethod()->GetTimeStep(timeIdx);

    fprintf(f,"1D : #%3d\tt=%6.2lf\tmin=%3d\tax=%3d\tmeanZ=%15.10lf\tmeanX=%15.10lf\tdriftCor=%15.10lf\tvarZ=%15.10lf\tvarX=%15.10lf\ttheo var=%15.10lf\n",
            timeIdx,time,GetMin(),GetMax(),treeMeanZ,treeMeanX,GetDriftCorrection(),treeVarZ,treeVarX,theoVar);
}
fclose(f);
****/

}


////////////////////////////////////////////////////
///	Class  : ARM_Slice1D
///	Routine: ComputeExpectation
///	Returns: double 
///	Action : 
////////////////////////////////////////////////////

double ARM_Slice1D::ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates ) const
{
    const ARM_Slice1D* nextSlice1D = nextSlice->ToSlice1D();

    int nextStateIdx = itsNextNodeIndex[stateIdx] - nextSlice1D->GetMin();

    double payoff = 0;
    if( itsProbaUp[stateIdx]>0 )
        payoff += itsProbaUp[stateIdx] * payoffFunc(nextStateIdx+1);

    if( itsProbaDown[stateIdx]>0 )
        payoff += itsProbaDown[stateIdx] * payoffFunc(nextStateIdx-1);

    payoff +=  payoffFunc(nextStateIdx) * 
        (1.0-itsProbaDown[stateIdx]-itsProbaUp[stateIdx]);

    return payoff;
}


////////////////////////////////////////////////////
///	Class  : ARM_Slice1D
///	Routine: toString
///	Returns: string 
///	Action : 
////////////////////////////////////////////////////

string ARM_Slice1D::toString( const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;
	os << indent << "1D Slice #" << dec << setw(3) << GetIndex() << "\t";
	os << "Min=" << dec << setw(5) << GetMin() << "\t";
	os << "Max=" << dec << setw(5) << GetMax() << "\t";
	os << "H=" << fixed << setw(10) << GetSpaceStep() << "\t";
	os << "Nb=" << dec << setw(6) << GetMax()-GetMin()+1;
	return os.str();
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Slice1DCstSymProba
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Slice1DCstSymProba
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_Slice1DCstSymProba& ARM_Slice1DCstSymProba::operator = (const ARM_Slice1DCstSymProba& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_Slice1DCstSymProba();
		new (this) ARM_Slice1DCstSymProba (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_Slice1DCstSymProba
///	Routine: LinkedWithNextSlice
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Slice1DCstSymProba::LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff)
{
    const ARM_PricingModel* model = sampler->GetModel();

    ARM_Slice1DBase* nextSlice1D    = nextSlice->ToSlice1DBase();
    ARM_Truncator1D* truncator1D    = truncator->ToTruncator1D();
    ARM_Sampler1DBase* sampler1D    = sampler->ToSampler1DBase();

    size_t nbStates = size();

    size_t timeIdx = GetIndex();
    double sliceYf = model->GetNumMethod()->GetTimeStep(timeIdx)/K_YEAR_LEN;

    /// Get variable states on the current slice
    ARM_GP_VectorPtr zStates = GetZStates();

    /// Save original X states
    SetXStates(sampler1D->ComputeZtoXStates(timeIdx,zStates));

    itsProbaUpDown.resize(nbStates);

    /// Reduce next slice range if truncation is needed
    /// Implicit straight connection
    int nextHedge;
    size_t stateIdx = 0;            // lower edge
    int nextNodeIndex = GetMin();   // lower edge
    if( truncator1D->NeedTruncation(sliceYf,this,stateIdx,nextSlice,nextNodeIndex) )
    {
        itsProbaUpDown[stateIdx]=0.0;
        nextSlice1D->SetMin(GetMin());
    }
    else
    {
        int nextMin = nextSlice1D->GetMin();
        nextHedge = GetMin()-1;
        nextSlice1D->SetMin(nextHedge > nextMin ? nextHedge : nextMin);
    }


    stateIdx = nbStates-1;   // upper edge
    nextNodeIndex = GetMax();       // upper edge
    if( truncator1D->NeedTruncation(sliceYf,this,stateIdx,nextSlice,nextNodeIndex) )
    {
        itsProbaUpDown[stateIdx]=0.0;
        nextSlice1D->SetMax(GetMax());
    }
    else
    {
        int nextMax = nextSlice1D->GetMax();
        nextHedge = GetMax()+1;
        nextSlice1D->SetMax(nextHedge < nextMax ? nextHedge : nextMax);
    }


    /// Calibrate tree if necessary
    bool needADPrices = model->NeedArrowDebreuPrices();
    bool isLocalDf = model->NeedLocalDiscount();
    ARM_GP_VectorPtr localDf(NULL);
    if(isLocalDf && needADPrices)
    {
        double nextSliceTime = model->GetNumMethod()->GetTimeStep(nextSlice->GetIndex());
        localDf = ComputeDriftCorrection(model,sampler,nextSliceTime,driftOffset,targetPayoff);
    }


    /// Update probabilities and Arrow-Debreu prices
    if(needProbas || needADPrices)
        UpdateProbasAndArrowDebreu(needProbas,needADPrices,isLocalDf,localDf,nextSlice1D);

}


////////////////////////////////////////////////////
///	Class  : ARM_Slice1DCstSymProba
///	Routine: ComputeExpectation
///	Returns: double 
///	Action : 
////////////////////////////////////////////////////

double ARM_Slice1DCstSymProba::ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates ) const
{
    const ARM_Slice1DBase* nextSlice1D = nextSlice->ToSlice1DBase();

    int nextNodeIdx = GetNextNodeIndex(stateIdx);
    int nextStateIdx = nextNodeIdx - nextSlice1D->GetMin();

    double payoff;
    if( itsProbaUpDown[stateIdx] == itsProba && itsProba>0)
    {
        /// No truncation
        payoff = itsProba * ( payoffFunc(nextStateIdx+1)
                              + payoffFunc(nextStateIdx-1)
                              - 2*payoffFunc(nextStateIdx) ) +
                 payoffFunc(nextStateIdx);
    }
    else
        /// On the upper or lower truncated edges => implicit straight connection
        payoff =  payoffFunc(nextStateIdx);


    return payoff;
}


////////////////////////////////////////////////////
///	Class  : ARM_Slice1DCstSymProba
///	Routine: toString
///	Returns: string 
///	Action : 
////////////////////////////////////////////////////

string ARM_Slice1DCstSymProba::toString( const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;
	os << indent << "1D CstSym Slice #" << dec << setw(3) << GetIndex() << "\t";
	os << "Min=" << dec << setw(5) << GetMin() << "\t";
	os << "Max=" << dec << setw(5) << GetMax() << "\t";
	os << "H=" << fixed << setw(10) << GetSpaceStep() << "\t";
	os << "Nb=" << dec << setw(6) << GetMax()-GetMin()+1;
	return os.str();
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	SliceNDFunct
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : SliceNDFunct
///	Routine: Constructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
SliceNDFunct::SliceNDFunct( size_t timeIdx, double dt, const ARM_PricingModel* model, const ARM_GP_MatrixPtr& xStates, const ARM_GP_VectorPtr& adPrices, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset)
:   
    itsTimeIdx(timeIdx),
    itsDt(dt),
	itsModel(model),
    itsXStates(xStates), 
	itsADPrices(adPrices),
	itsDriftCorrection(driftCorrection),
	itsDriftOffset(driftOffset),
    itsStates(new ARM_PricingStates(xStates->cols(),xStates->rows(),0,xStates->rows()))
{
    /// Copy states in ARM_PricingStates object for all directions except
    /// the driftOffset one that will be updated by the calibrated drift correction
    /// Non calibrated directions are adjusted by input drift corrections
    size_t i,nbStates = itsXStates->cols();
    size_t j,nbDims = itsXStates->rows();
    for(j=0;j<nbDims;++j)
    {
        if(j!=itsDriftOffset)
        {
            for(i=0;i<nbStates;++i)
                itsStates->SetNumMethodState(i,j,(*itsXStates)(j,i) + (*itsDriftCorrection)[j]);
        }
    }
}

////////////////////////////////////////////////////
///	Class  : SliceNDFunct
///	Routine: operator()
///	Returns: double
///	Action : Compute the functor value
///          The drift correction is arbitrarily set
///          to the first numethod variable
////////////////////////////////////////////////////
double SliceNDFunct::operator() (double x) const
{
    size_t i,nbStates = itsXStates->cols();

    /// Correct driftOffset direction only
    for(i=0;i<nbStates;++i)
        itsStates->SetNumMethodState(i,itsDriftOffset,(*itsXStates)(itsDriftOffset,i) + x);

	itsModel->TreeStatesToModelStates(itsStates,itsTimeIdx);

    ARM_GP_VectorPtr localPayoff = itsModel->LocalPayoffs(itsTimeIdx,itsDt,itsStates);

    double fx=0.0;
    for(i=0;i<localPayoff->size();++i)
        fx += (*itsADPrices)[i] * (*localPayoff)[i];
    return fx;
}

////////////////////////////////////////////////////
///	Class  : SliceNDFunct
///	Routine: GetLocalDf
///	Returns: ARM_GP_VectorPtr
///	Action : Compute local discount associated to a
///        : deterministic shift x
////////////////////////////////////////////////////
ARM_GP_VectorPtr SliceNDFunct::GetLocalDf(double x) const
{
    size_t i,nbStates = itsXStates->cols();

    /// Correct driftOffset direction only
    for(i=0;i<nbStates;++i)
        itsStates->SetNumMethodState(i,itsDriftOffset,(*itsXStates)(itsDriftOffset,i) + x);

	itsModel->TreeStatesToModelStates(itsStates,itsTimeIdx);

    ARM_GP_VectorPtr localPayoff = itsModel->LocalDiscounts(itsTimeIdx,itsDt,itsStates);
    return localPayoff;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_SliceNDBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_SliceNDBase
///	Routine: ARM_SliceNDBase
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////
ARM_SliceNDBase::ARM_SliceNDBase(const ARM_SliceNDBase& rhs)
:   ARM_SliceBase(rhs),
    itsMin(rhs.itsMin), itsMax(rhs.itsMax), itsSpaceStep(rhs.itsSpaceStep),
    itsNextNodeIndex(rhs.itsNextNodeIndex)
{
    itsTransitor = (rhs.itsTransitor ? static_cast< ARM_TransitorBase* >(rhs.itsTransitor->Clone()) : NULL);

    ARM_GP_VectorPtr nullV(NULL);
	itsDriftCorrection =    (rhs.itsDriftCorrection != nullV
                            ? ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsDriftCorrection->Clone()))
                            : nullV);

    ARM_GP_MatrixPtr nullM(NULL);
	itsXStates =    (rhs.itsXStates != nullM
                    ? ARM_GP_MatrixPtr(static_cast<ARM_GP_Matrix*>(rhs.itsXStates->Clone()))
                    : nullM);
}

////////////////////////////////////////////////////
///	Class  : ARM_SliceNDBase
///	Routine: ComputeExpectation
///	Returns: double 
///	Action : 
////////////////////////////////////////////////////

double ARM_SliceNDBase::ComputeExpectation( size_t stateIdx, const ARM_SliceBase* nextSlice, const PayoffFunc& payoffFunc, ARM_TransitionStates& transStates ) const
{
    /// Compute all actual transitions (proba > 0) from current node to next slice
    if(transStates.GetUsedSize()==0)
        GetTransitor().ComputeTransitions(this,stateIdx,nextSlice,transStates);

    /// Compute payoff expectation
    const ARM_IntVector& nextStateIndexes = transStates.GetStates();
    const ARM_GP_Vector& transProbas = transStates.GetProbas();
    double payoff = 0;
    for(size_t i=0;i<transStates.GetUsedSize();++i)
        payoff += transProbas[i] *(payoffFunc)(nextStateIndexes[i]);

    return payoff;
}


////////////////////////////////////////////////////
///	Class  : ARM_SliceNDBase
///	Routine: ComputeDriftCorrection
///	Returns: 
///	Action : Compute the correction to be added
///          to states to fit the discount factor
///          from 0 to next slice time
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_SliceNDBase::ComputeDriftCorrection(const ARM_PricingModel* model, ARM_SamplerBase* sampler, double nextTime, size_t driftOffset, double targetPayoff, bool isLocalDf)
{
    /// Get target discount price
    size_t timeIdx = GetIndex();
    double time = model->GetNumMethod()->GetTimeStep(timeIdx);

    double nrTol = 1.0e-8;

    /// Generate X states if necessary
    if(itsXStates == ARM_GP_MatrixPtr(NULL))
    {
        /// Not previously computed then compute Z states...
        ARM_GP_MatrixPtr zStates = GetZStates();

        /// ... and restore X states through the sampler
        const ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();
        itsXStates = samplerND->ComputeZtoXStates(timeIdx,zStates);
        
    }

    /// Compute the initial guess
    SliceNDFunct sliceNDFct(timeIdx,nextTime-time,model,itsXStates,GetArrowDebreuPrices(),itsDriftCorrection,driftOffset);
    UnaryFuncWithNumDerivative<double,double> fctToSolve( sliceNDFct );
    double df0 = fctToSolve(0.0);
    double ddf0 = (*(fctToSolve.Derivative()))(0.0);
    double x0 = targetPayoff - df0;
    double driftCorrection = 0.0;
    if(x0 < -nrTol || x0 > nrTol)
    {
        if(fabs(ddf0)<K_NEW_DOUBLE_TOL)
            ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME +
            " : can't find an initial solution to compute drift correction for current slice" );
        x0 /= ddf0;

        /// Create a NR solver then find the solution
		const double xtol = 1.0e-8;
		const size_t max_iter = 30;
        T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double,double> > solver(fctToSolve,targetPayoff,nrTol,xtol,max_iter);

		/// To initialize departure point
		solver.setInitialGuess(x0);

        driftCorrection=solver.Solve();

/****
df0=fctToSolve(driftCorrection);
ARM_GP_VectorPtr localPayoff = model->LocalPayoffs(timeIdx,nextTime-time,sliceNDFct.GetStates());
ARM_GP_VectorPtr localDf = model->LocalDiscounts(timeIdx,nextTime-time,sliceNDFct.GetStates());
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"\n\ntimeIdx=%2d\n",timeIdx);
for(size_t i=0;i<localDf->size();++i)
    fprintf(f,"spotFx=%15.10lf\t localDf=%15.13lf\n",(*localPayoff)[i],(*localDf)[i]);
fclose(f);
****/

    }

    (*itsDriftCorrection)[driftOffset] = driftCorrection;

    if(isLocalDf)
        /// Compute local DF for further Arrow-Debreu updating
        return sliceNDFct.GetLocalDf(driftCorrection);
    else
        return ARM_GP_VectorPtr(NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_SliceNDBase
///	Routine: UpdateProbasAndArrowDebreu
///	Returns: 
///	Action : Compute and set probabilities to reach
///          states of the slice. Idem for Arrow-Debreu
///          prices
////////////////////////////////////////////////////
void ARM_SliceNDBase::UpdateProbasAndArrowDebreu( bool needProbas, bool needADPrices, bool isLocalDf, ARM_GP_VectorPtr localDf, ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates)
{
    size_t i,stateIdx,nbStates = size();
    size_t nextNbStates = nextSlice->size();

    ARM_GP_Vector* spotProbas=NULL;
    if(needProbas) spotProbas = new ARM_GP_Vector(nextNbStates,0.0);

    ARM_GP_Vector* arrowDebreuPrices=NULL;
    if(needADPrices) arrowDebreuPrices = new ARM_GP_Vector(nextNbStates,0.0);

    double proba,ADPrice;
    ARM_TreeIndex index(GetMin(),GetMax());
    for( index.Reset(); (stateIdx=index.GetPosition()) < nbStates; ++index )
    {
        /// Compute all actual transitions (proba > 0) from current node to next slice
        GetTransitor().ComputeTransitions(this,stateIdx,nextSlice,transStates);

        const ARM_IntVector& nextNodeIndexes = transStates.GetStates();
        const ARM_GP_Vector& transProbas = transStates.GetProbas();

        if(needProbas)
        {
            /// Update states probas
            proba = GetSpotProbas(stateIdx);
            for(i=0;i<transStates.GetUsedSize();++i)
                (*spotProbas)[nextNodeIndexes[i]] += proba * transProbas[i];
        }

        if(needADPrices)
        {
            /// Update Arrow-Debreu prices
            ADPrice = (isLocalDf ? GetArrowDebreuPrices(stateIdx)*(*localDf)[stateIdx] :GetArrowDebreuPrices(stateIdx));
            for(i=0;i<transStates.GetUsedSize();++i)
                (*arrowDebreuPrices)[nextNodeIndexes[i]] += ADPrice * transProbas[i];
        }
    } // loop over states


/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
size_t nextIndex = GetIndex()+1;
fprintf(f,"TimeIdx=%3d\tNbStates=%7d\n",nextIndex,nextNbStates);
if(nextIndex == 41)
{
    ARM_IntVector minIndex(((ARM_SliceND*)nextSlice)->GetMin());
    ARM_IntVector maxIndex(((ARM_SliceND*)nextSlice)->GetMax());
    ARM_TreeIndex index(minIndex,maxIndex);
    double probaSum=0.0,proba,adPriceSum=0.0,adPrice;
    fprintf(f,"Min=(%3d,%3d,%3d)\tMax=(%3d,%3d,%3d)\n",minIndex[0],minIndex[1],minIndex[2],maxIndex[0],maxIndex[1],maxIndex[2]);
    for( index.Reset(); (stateIdx=index.GetPosition()) < nextNbStates; ++index )
    {
        proba = needProbas ? (*spotProbas)[stateIdx] : 0.0;
        probaSum += proba;
        adPrice = needADPrices ? (*arrowDebreuPrices)[stateIdx] : 0.0;
        adPriceSum += adPrice;
        fprintf(f,"  Idx=%7d\tI=(%3d,%3d,%3d)\tPr==%15.13lf\tAD=%15.13lf\n",stateIdx,index[0],index[1],index[2],proba,adPrice);
    }
    fprintf(f,"  PrSum==%15.13lf\tADSum=%15.13lf\n",probaSum,adPriceSum);
}
fclose(f);
****/

    nextSlice->SetSpotProbas(ARM_GP_VectorPtr(spotProbas));
    nextSlice->SetArrowDebreuPrices(ARM_GP_VectorPtr(arrowDebreuPrices));

#if defined( __GP_STRICT_VALIDATION)
    /// Check probabilities of next slice (AD prices without localDf <=> probabilities)
    if(needProbas || (needADPrices && !isLocalDf))
    {
        double probaSum=0.0;
        if(needProbas)
        {
            for( stateIdx=0; stateIdx<nextNbStates; ++stateIdx )
                probaSum += (*spotProbas)[stateIdx];
        }
        else
        {
            for( stateIdx=0; stateIdx<nextNbStates; ++stateIdx )
                probaSum += (*arrowDebreuPrices)[stateIdx];
        }

        if(probaSum - 1 < - K_NEW_DOUBLE_TOL || probaSum - 1 > K_NEW_DOUBLE_TOL)
            ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME + " : proba sum is not equal to 1 in the current slice" );
    }
#endif
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_SliceND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_SliceND
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_SliceND& ARM_SliceND::operator = (const ARM_SliceND& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_SliceND();
		new (this) ARM_SliceND (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_SliceND
///	Routine: LinkedWithNextSlice
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SliceND::LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff)
{
    const ARM_PricingModel* model = sampler->GetModel();

    ARM_SliceNDBase* nextSliceND    = nextSlice->ToSliceNDBase();
    ARM_TruncatorND* truncatorND    = truncator->ToTruncatorND();
    ARM_SamplerNDBase* samplerND    = sampler->ToSamplerNDBase();
    size_t i,nbDims = dim();
    size_t stateIdx,nbStates = size();

    /// Set previous drift corrections if any
    if(driftCorrection == ARM_GP_VectorPtr(NULL) || driftCorrection->size() != nbDims)
        SetDriftCorrectionVect(ARM_GP_VectorPtr( new ARM_GP_Vector(nbDims,0.0)));
    else
        SetDriftCorrectionVect(driftCorrection);

    size_t timeIdx = GetIndex();
    if(timeIdx + 2 == model->GetNumMethod()->GetTimeSteps()->size())
    {
        /// Initialise drift correction for next slice = last slice
        nextSliceND->SetDriftCorrectionVect(ARM_GP_VectorPtr( new ARM_GP_Vector(nbDims,0.0)));
    }
    double sliceYf = model->GetNumMethod()->GetTimeStep(timeIdx)/K_YEAR_LEN;


    /// Calibrate tree if necessary. Must be done before connecting to
    /// next slice because integrated version of the markovian drift
    /// may use the actual spot values at current slice
    bool needADPrices = model->NeedArrowDebreuPrices();
    bool isLocalDf = model->NeedLocalDiscount();
    ARM_GP_VectorPtr localDf(NULL);
    if(isLocalDf && needADPrices)
    {
        double nextSliceTime = model->GetNumMethod()->GetTimeStep(nextSlice->GetIndex());
        localDf = ComputeDriftCorrection(model,sampler,nextSliceTime,driftOffset,targetPayoff);
    }


    /// Get variable states on the current slice
    ARM_GP_MatrixPtr zStates = GetZStates();

    /// Get conditional mean & variance and save X states
    ARM_GP_MatrixPtr xStates;
    ARM_GP_MatrixPtr mean       = samplerND->ComputeConditionalMean( timeIdx,zStates,GetDriftCorrectionVect(),xStates );
    SetXStates(xStates);
    const ARM_GP_Vector& var    = samplerND->GetLocalVar( timeIdx );

    ARM_GP_Vector spaceStepSquare(nextSliceND->GetSpaceStep());
    for( i=0; i<nbDims; ++i )
        spaceStepSquare[i] *= spaceStepSquare[i];

    itsProbaUp.resize(nbStates,nbDims);
    itsProbaDown.resize(nbStates,nbDims);
    ResizeNextNodeIndex(nbStates,nbDims);

    ARM_IntVector minIndex(nbDims),maxIndex(nbDims),nextNodeIndex(nbDims);
    for( i=0; i<nbDims; ++i )
    {
        minIndex[i]=nextSliceND->GetMax(i);
        maxIndex[i]=nextSliceND->GetMin(i);
    }

    /// Connect to next slice
    double meanError;
    ARM_GP_Vector meanRel(nbDims),varRel(nbDims);
    double probaUp,probaDown;
    ARM_TreeIndex index(GetMin(),GetMax());
    for( index.Reset(); (stateIdx=index.GetPosition()) < nbStates; ++index )
    {
        /// Theorical transition center in the next slice and
        /// associated probabilities to fulfill mean & variance
        for( i=0; i<nbDims; ++i )
        {
            meanRel[i]                  = (*mean)(i,stateIdx)/nextSliceND->GetSpaceStep(i);
            nextNodeIndex[i]            = CC_Round(meanRel[i]);
            varRel[i]                   = var[i]/spaceStepSquare[i];
            meanError                   = meanRel[i]-nextNodeIndex[i];
            probaUp                     = 0.5 * (varRel[i]+meanError*meanError+meanError);
            probaDown                   = probaUp-meanError;
            itsProbaUp(stateIdx,i)      = probaUp;
            itsProbaDown(stateIdx,i)    = probaDown;

            SetNextNodeIndex(stateIdx,i,nextNodeIndex[i]);
        }

        /// Truncation if necessary
        truncatorND->UpdateNode(sliceYf,this,stateIdx,nextSlice,meanRel,varRel);

        /// Update multi-index range
        for( i=0; i<nbDims; ++i )
        {
            if( itsProbaUp(stateIdx,i) > 0)
                maxIndex[i] = CC_Max( maxIndex[i], GetNextNodeIndex(stateIdx,i)+1 );
            else
                maxIndex[i] = CC_Max( maxIndex[i], GetNextNodeIndex(stateIdx,i) );

            if( itsProbaDown(stateIdx,i) > 0)
                minIndex[i] = CC_Min( minIndex[i], GetNextNodeIndex(stateIdx,i)-1 );
            else
                minIndex[i] = CC_Min( minIndex[i], GetNextNodeIndex(stateIdx,i) );
        }
    }

    /// Update slice range
    for( i=0; i<nbDims; ++i )
    {
        nextSliceND->SetMin(i,minIndex[i]);
        nextSliceND->SetMax(i,maxIndex[i]);
    }

    /// Update probabilities and Arrow-Debreu prices
    size_t nextNbStates = nextSliceND->size();
    if(needProbas || needADPrices)
        UpdateProbasAndArrowDebreu(needProbas,needADPrices,isLocalDf,localDf,nextSlice,transStates);


/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
if(needProbas)
{
    /// Check mean & variance of current slice (on mean reverting process X)
    ARM_SliceNDBase* curSlice = ToSliceNDBase();
    bool isPrintable = true;
    while(isPrintable)
    {
        nbStates = curSlice->size();
        ARM_GP_VectorPtr spotProbas = curSlice->GetSpotProbas();
        ARM_GP_MatrixPtr zStates = curSlice->GetZStates();
        ARM_GP_MatrixPtr xStates = samplerND->ComputeZtoXStates(timeIdx,zStates);

        ARM_GP_Vector treeMeanZ(nbDims,0.0),treeVarZ(nbDims,0.0);
        ARM_GP_Vector treeMeanX(nbDims,0.0),treeVarX(nbDims,0.0);
        ARM_GP_Vector treeSkewX(nbDims,0.0),treeKurtX(nbDims,0.0);
        double z,x,y,probaSum=0.0;
        for( stateIdx=0; stateIdx<nbStates; ++stateIdx )
        {
            probaSum += (*spotProbas)[stateIdx];
            for( i=0; i<nbDims; ++i )
            {
                z = (*zStates)(i,stateIdx);
                treeMeanZ[i] += z * (*spotProbas)[stateIdx];
                treeVarZ[i] += z*z* (*spotProbas)[stateIdx];
                x = (*xStates)(i,stateIdx);
                y = x;
                treeMeanX[i] += y * (*spotProbas)[stateIdx];
                y *= x;
                treeVarX[i] += y * (*spotProbas)[stateIdx];
                y *= x;
                treeSkewX[i] += y * (*spotProbas)[stateIdx];
                y *= x;
                treeKurtX[i] += y * (*spotProbas)[stateIdx];
            }
        }

        ARM_MatrixVector localVariances,variances;
        model->NumMethodStateLocalGlobalVariances(*(model->GetNumMethod()->GetTimeSteps()),localVariances,variances);
        ARM_GP_Vector theoVar(nbDims),theoSkew(nbDims,0.0),theoKurt(nbDims);
        ARM_GP_Vector errRelVar(nbDims,0.0),errRelKurt(nbDims,0.0);

        double x2,x3;
        for( i=0; i<nbDims; ++i )
        {
            treeVarZ[i] -= treeMeanZ[i]*treeMeanZ[i];
            x = treeMeanX[i];
            x2 = treeVarX[i];
            x3 = treeSkewX[i];
            treeVarX[i] -= x*x;
            treeSkewX[i] += -3*x2*x + 2*x*x*x;
            treeKurtX[i] += -4*x3*x + 6*x2*x*x - 3*x*x*x*x;
            theoVar[i] = (*(variances[timeIdx]))(i,i);
            treeKurtX[i] /= (treeVarX[i]*treeVarX[i]);
            theoKurt[i] = 3.0;

            if(theoVar[i] != 0.0)
                errRelVar[i] = (treeVarX[i]-theoVar[i])/theoVar[i]*100;
            if(theoKurt[i] != 0.0)
                errRelKurt[i] = (treeKurtX[i]-theoKurt[i])/theoKurt[i]*100;
        }

        DeletePointorVector<ARM_GP_Matrix>( localVariances );
        DeletePointorVector<ARM_GP_Matrix>( variances );

        double time = model->GetNumMethod()->GetTimeStep(timeIdx);

        fprintf(f,"#%3d\tt=\t%6.2lf\tPrSum=\t%15.10lf\t",timeIdx,time,probaSum);
        for( i=0; i<nbDims; ++i )
            fprintf(f,"     %1dD : min=\t%3d\tmax=\t%3d\tstep=\t%8.6lf\tmeanZ=\t%8.6lf\tmeanX=\t%8.6lf\t\tvarZ=\t%8.6lf\tvarX=\t%8.6lf\ttheo var=\t%8.6lf\terrVarX=\t%8.6lf (%5.2lf%%)\tskewX=\t%8.6lf\ttheo skew=\t%8.6lf\terrSkewX=\t%8.6lf\tkurtX=\t%9.6lf\ttheo kurt=\t%8.6lf\terrKurtX=\t%9.6lf (%5.2lf%%)\t",
                i+1,GetMin(i),GetMax(i),GetSpaceStep(i),treeMeanZ[i],treeMeanX[i],treeVarZ[i],treeVarX[i],theoVar[i],treeVarX[i]-theoVar[i],errRelVar[i],treeSkewX[i],theoSkew[i],treeSkewX[i]-theoSkew[i],treeKurtX[i],theoKurt[i],treeKurtX[i]-theoKurt[i],errRelKurt[i]);
        fprintf(f,"\n");
        if(timeIdx+1 == model->GetNumMethod()->GetLastTimeIdx())
        {
            /// Idem for last slice
            curSlice = nextSliceND;
            ++timeIdx;
        }
        else
            isPrintable=false;
    }
}
fclose(f);
****/

}


////////////////////////////////////////////////////
///	Class  : ARM_SliceND
///	Routine: toString
///	Returns: string 
///	Action : 
////////////////////////////////////////////////////

string ARM_SliceND::toString( const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;
    size_t nbDims = dim();
	os << indent << dec << setw(1) << nbDims << "D Slice #" << setw(3) << GetIndex() << "\t";
    size_t j,nbNodes=1;
    os << "Min=(" << dec;
    for(j=0;j<nbDims;++j)
    {
	    os << setw(3) << GetMin(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }
    os << "Max=(" << dec;
    for(j=0;j<nbDims;++j)
    {
	    os << setw(3) << GetMax(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }
    os << "H=(" << fixed;
    for(j=0;j<nbDims;++j)
    {
        nbNodes *= GetMax(j) - GetMin(j) + 1;
	    os << setw(10) << GetSpaceStep(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }

	os << "Nb=" << dec << setw(6) << nbNodes;

    return os.str();
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_SliceNDCstSymProba
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_SliceNDCstSymProba
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_SliceNDCstSymProba& ARM_SliceNDCstSymProba::operator = (const ARM_SliceNDCstSymProba& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_SliceNDCstSymProba();
		new (this) ARM_SliceNDCstSymProba (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_SliceNDCstSymProba
///	Routine: LinkedWithNextSlice
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SliceNDCstSymProba::LinkedWithNextSlice( ARM_SliceBase* nextSlice, ARM_SamplerBase* sampler, ARM_TruncatorBase* truncator, bool needProbas, ARM_TransitionStates& transStates, const ARM_GP_VectorPtr& driftCorrection, size_t driftOffset, double targetPayoff)
{
    const ARM_PricingModel* model = sampler->GetModel();

    ARM_SliceNDBase* nextSliceND    = nextSlice->ToSliceNDBase();
    ARM_TruncatorND* truncatorND    = truncator->ToTruncatorND();
    ARM_SamplerNDBase* samplerND    = sampler->ToSamplerNDBase();

    size_t timeIdx = GetIndex();
    double sliceYf = model->GetNumMethod()->GetTimeStep(timeIdx)/K_YEAR_LEN;

    /// Get variable states on the current slice (and set implicit node connection)
    ARM_GP_MatrixPtr zStates = GetZStates();

    /// Save original X states
    SetXStates(samplerND->ComputeZtoXStates(timeIdx,zStates));

    /// Set drift corrections if any
    size_t nbDims = dim();
    size_t i,stateIdx,nbStates = size();
    if(driftCorrection == ARM_GP_VectorPtr(NULL) || driftCorrection->size() != nbDims)
        SetDriftCorrectionVect(ARM_GP_VectorPtr( new ARM_GP_Vector(nbDims,0.0)));
    else
        SetDriftCorrectionVect(driftCorrection);

    if(timeIdx + 2 == model->GetNumMethod()->GetTimeSteps()->size())
    {
        /// Initialise drift correction for next = last slice
        nextSliceND->SetDriftCorrectionVect(ARM_GP_VectorPtr( new ARM_GP_Vector(nbDims,0.0)));
    }

    /// Identification of nodes located on a truncated hedge
    itsTruncatedHedge.resize(nbStates,nbDims);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        for(i=0;i<nbDims;++i)
            itsTruncatedHedge(stateIdx,i)=false;


    ARM_IntVector maxIndex(GetMin());
    ARM_IntVector minIndex(GetMax());
    ARM_BoolVector isTruncated;
    ARM_TreeIndex index(GetMin(),GetMax());
    int nextHedgeNodeIdx;
    for( index.Reset(); (stateIdx=index.GetPosition()) < nbStates; ++index )
    {
        /// Straight connexion
        for(i=0;i<nbDims;++i)
            SetNextNodeIndex(stateIdx,i,index[i]);

        if(index.IsOnHedge())
        {
            /// Test if truncation is needed and then reduce next slice range
            /// nextNodeIndex = index because of straight connexion
            isTruncated = truncatorND->NeedTruncation(sliceYf,this,stateIdx,nextSlice,index);
            for(i=0;i<nbDims;++i)
            {
                if(isTruncated[i])
                {
                    itsTruncatedHedge(stateIdx,i) = true;
                    if(maxIndex[i] < index[i])
                        maxIndex[i] = index[i];
                    if(minIndex[i] > index[i])
                        minIndex[i] = index[i];
                }
                else
                {
                    if(maxIndex[i] < (nextHedgeNodeIdx=index[i]+1))
                        maxIndex[i] = nextHedgeNodeIdx;
                    if(minIndex[i] > (nextHedgeNodeIdx=index[i]-1))
                        minIndex[i] = nextHedgeNodeIdx;
                }
            }
        }
    } // for states of current slice


    /// Update next slice range
    nextSliceND->SetMin(minIndex);
    nextSliceND->SetMax(maxIndex);


    /// Calibrate tree if necessary
    bool needADPrices = model->NeedArrowDebreuPrices();
    bool isLocalDf = model->NeedLocalDiscount();
    ARM_GP_VectorPtr localDf(NULL);
    if(isLocalDf && needADPrices)
    {
        double nextSliceTime = model->GetNumMethod()->GetTimeStep(nextSlice->GetIndex());
        localDf = ComputeDriftCorrection(model,sampler,nextSliceTime,driftOffset,targetPayoff);
    }


    /// Update probabilities and Arrow-Debreu prices
    if(needProbas || needADPrices)
        UpdateProbasAndArrowDebreu(needProbas,needADPrices,isLocalDf,localDf,nextSlice,transStates);
}


////////////////////////////////////////////////////
///	Class  : ARM_SliceNDCstSymProba
///	Routine: toString
///	Returns: string 
///	Action : 
////////////////////////////////////////////////////

string ARM_SliceNDCstSymProba::toString( const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;
    size_t nbDims = dim();
	os << indent << dec << setw(1) << nbDims << "D CstSym Slice #" << setw(3) << GetIndex() << "\t";
    size_t j,nbNodes=1;
    os << "Min=(" << dec;
    for(j=0;j<nbDims;++j)
    {
	    os << setw(3) << GetMin(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }
    os << "Max=(" << dec;
    for(j=0;j<nbDims;++j)
    {
	    os << setw(3) << GetMax(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }
    os << "H=(" << fixed;
    for(j=0;j<nbDims;++j)
    {
	    os << setw(10) << GetSpaceStep(j);
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }

    os << "P=(" << fixed;
    for(j=0;j<nbDims;++j)
    {
        nbNodes *= GetMax(j) - GetMin(j) + 1;
	    os << setw(10) << itsProba[j];
        if(j<nbDims-1)
            os << ",";
        else
            os << ")\t";
    }

	os << "Nb=" << dec << setw(6) << nbNodes;

    return os.str();
}


CC_END_NAMESPACE()


#undef FIRST_STATE_VARIABLE

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

