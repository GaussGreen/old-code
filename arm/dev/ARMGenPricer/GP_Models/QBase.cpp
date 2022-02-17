/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QBase.cpp
 *
 *  \brief base class for the q model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpmodels/QBase.h"

/// gpbase
#include "gpbase/singleton.h" /// for numeraire factory

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"

///gpcalib
#include "gpcalib/calibmethod.h"
 
CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: CopyNoCleanUp
///	Returns: void 
///	Action : attribute object copy
////////////////////////////////////////////////////
void ARM_QModelBase::CopyNoCleanUp(const ARM_QModelBase& rhs)
{
    itsPrecomputedFwds = rhs.itsPrecomputedFwds? (ARM_GP_Vector*) rhs.itsPrecomputedFwds->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: CleanUp
///	Returns: void 
///	Action : attribute object destruction
////////////////////////////////////////////////////
void ARM_QModelBase::CleanUp()
{
    delete itsPrecomputedFwds;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBase
///	Routines: MappingFunctionDerivative
///	Returns :
///	Action  : mapping function derivative of the Q model
////////////////////////////////////////////////////

double ARM_QModelBase::MappingFunctionDerivative( double x, double f0, double q0 ) const
{
	if(q0<K_NEW_DOUBLE_TOL)
		return f0;
	else 
		return f0*exp(q0*x);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelBase::ARM_QModelBase(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params) 
:	ARM_PricingModel(zc,&params), itsPrecomputedFwds(NULL)
{
    ARM_NumerairePtr numerairePtr( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
	SetNumeraire( numerairePtr );
}



////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelBase::ARM_QModelBase(const ARM_QModelBase& rhs)
:	ARM_PricingModel(rhs), itsPrecomputedFwds(NULL)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QModelBase& ARM_QModelBase::operator=(const ARM_QModelBase& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModel::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_QModelBase::~ARM_QModelBase()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QModelBase::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL))
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numerical method not set in the Q mode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

        ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numeraire not set in the Qmode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
            
		}

        if( numeraire->GetType() != ARM_Numeraire::Cash )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            ": only Cash numeraire supported by Q model at the moment!");
		
		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_QModelBase::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb= FactorCount();
	size_t statesNb = states->size();
	double currentState;
	size_t modelNb	= GetModelNb();
	
	for( size_t i=0;i<statesNb; ++i )
	{
		for( size_t j=0;  j<factorsNb; ++j )
		{
			currentState   = 0.0;
			for( size_t k =0; k<=j; ++k )
			{
				double gaussian = states->GetNumMethodState(i,modelNb+k);
				currentState  += gaussian;
			}
			states->SetModelState(i,j+modelNb,currentState);
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_QModelBase::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}

////////////////////////////////////////////////////
///	Class   : ARM_QModelBase
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_QModelBase::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: PostInit
///	Returns: nothing
///	Action : Function to init after the numerical method
///			 has been initialised.. non const
////////////////////////////////////////////////////
void ARM_QModelBase::PostInit()
{
	delete itsPrecomputedFwds;
    ARM_GP_Vector* timeSteps = GetNumMethod()->GetTimeSteps();
	size_t timeStepSize = timeSteps->size();
	itsPrecomputedFwds = new ARM_GP_Vector(timeStepSize,0.0);
	
	for( size_t i=0; i<timeStepSize; ++i )
		(*itsPrecomputedFwds)[i] = ComputeFwdAtTime( (*timeSteps)[i] );
}

////////////////////////////////////////////////////
///	Class  : ARM_QModelBase
///	Routine: NeedMCIntegProcess
///	Returns: nothing
///	Action : This function tells to MC to return or not
/// an integrated process
////////////////////////////////////////////////////

ARM_BoolVector ARM_QModelBase::NeedMCIntegProcess() const
{
	size_t factorNb = FactorCount();

	ARM_BoolVector ret;

	for (size_t i = 0; i < factorNb; ++i)
		ret.push_back(false);

	return ret;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

