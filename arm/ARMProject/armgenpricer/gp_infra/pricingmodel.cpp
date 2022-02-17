/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodel.cpp
 *
 *  \brief pricing model is the base class for all models
 *
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */

/// the header of the file
#include "gpinfra/pricingmodel.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/curve.h"
#include "gpbase/surface.h"
#include "gpbase/datestrip.h"

///gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/vanillacap.h"

/// gpmodels
#include "gpmodels/BS_ModelParams.h"
#include "gpmodels/BS_Model.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/surfacemodelparam.h"

/// gpbase
#include <ccy/currency.h>

/// STL
#include <memory> /// for auto_ptr

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_PricingModel::CopyNoCleanUp(const ARM_PricingModel& rhs)
{
	itsModelName		= rhs.itsModelName;
	itsPayModelName		= rhs.itsPayModelName;
	itsZeroCurve		= rhs.itsZeroCurve;
    itsParams           = rhs.itsParams ? (ARM_ModelParams*) rhs.itsParams->Clone(): NULL;
    itsDensityFunctor   = rhs.itsDensityFunctor ? (ARM_DensityFunctor*) rhs.itsDensityFunctor->Clone(): NULL;
    itsNumeraire        = rhs.itsNumeraire;
    itsNumMethod        = rhs.itsNumMethod;

    itsDiscountFunctor  = rhs.itsDiscountFunctor->Duplicate(rhs,this);
    itsFixingFunctor    = rhs.itsFixingFunctor->Duplicate(rhs,this);
    itsMktDataManager   = rhs.itsMktDataManager ? (ARM_MarketData_ManagerRep*) rhs.itsMktDataManager->Clone() : NULL;
    itsMDMKeys          = rhs.itsMDMKeys;
	itsModelNb			= rhs.itsModelNb;
	itsModelRank		= rhs.itsModelRank;

	/// part on the std Dev and local vars
	itsModelStateLocalVarsIsShared		= rhs.itsModelStateLocalVarsIsShared;
	itsModelStateLocalStdDevsIsShared	= rhs.itsModelStateLocalStdDevsIsShared;
	itsNumMethodStateLocalVarsIsShared  = rhs.itsNumMethodStateLocalVarsIsShared;

	if ( itsModelStateLocalStdDevsIsShared )	
		itsModelStateLocalStdDevs = rhs.itsModelStateLocalStdDevs; 	
	else
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelStateLocalStdDevs, itsModelStateLocalStdDevs);	

	if ( itsModelStateLocalVarsIsShared )
		itsModelStateLocalVars = rhs.itsModelStateLocalVars;
	else
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelStateLocalVars, itsModelStateLocalVars);

	if (itsNumMethodStateLocalVarsIsShared) 
	{
		itsNumMethodStateLocalVars = rhs.itsNumMethodStateLocalVars;
		itsNumMethodStateLocalStdDevs = rhs.itsNumMethodStateLocalStdDevs;
	}
	else 
	{
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsNumMethodStateLocalVars, itsModelStateLocalVars);
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsNumMethodStateLocalStdDevs, itsNumMethodStateLocalStdDevs);
	}

	itsFromMultiFactor	= rhs.itsFromMultiFactor;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_PricingModel::CleanUp()
{
    delete itsParams;
    itsParams=NULL;

	delete itsDensityFunctor; 
	itsDensityFunctor = NULL;

	delete itsDiscountFunctor;
    itsDiscountFunctor=NULL;
    
    delete itsFixingFunctor;
    itsFixingFunctor=NULL;

    delete itsMktDataManager;
    itsMktDataManager=NULL;

	if( !itsModelStateLocalStdDevsIsShared ){ DeletePointorVector<ARM_GP_Matrix>( itsModelStateLocalStdDevs );}
	if( !itsModelStateLocalVarsIsShared) { DeletePointorVector<ARM_GP_Matrix>( itsModelStateLocalVars ); }

	if (!itsNumMethodStateLocalVarsIsShared) DeletePointorVector<ARM_GP_Matrix>( itsNumMethodStateLocalStdDevs );
	if (!itsNumMethodStateLocalVarsIsShared) DeletePointorVector<ARM_GP_Matrix>( itsNumMethodStateLocalVars );
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor by default a model is never calibrated
///			Note that it never check the sense of the model params!
////////////////////////////////////////////////////
ARM_PricingModel::ARM_PricingModel( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params, const ARM_DensityFunctor* densityFct)
:
	ARM_RootObject(),
	itsModelName( "NoName" ),
	itsPayModelName("NoName"),
	itsZeroCurve( zc ),
	itsParams( params ? static_cast<ARM_ModelParams*>(params->Clone()): NULL),
	itsDensityFunctor( densityFct ? static_cast<ARM_DensityFunctor*>(densityFct->Clone()): NULL),
	itsNumeraire(NULL), 
	itsNumMethod(NULL),
    itsDiscountFunctor(NULL),
    itsFixingFunctor(NULL),
    itsMktDataManager(NULL),
	itsMDMKeys(0),
	itsModelStateLocalVars(0),
	itsModelStateLocalVarsIsShared(false),
	itsModelStateLocalStdDevs(0),
	itsModelStateLocalStdDevsIsShared(false),
	itsNumMethodStateLocalVarsIsShared(false),
	itsNumMethodStateLocalVars(0),
	itsNumMethodStateLocalStdDevs(0),
	itsModelNb(0),
	itsModelRank(0),
	itsFromMultiFactor(false)
{
    itsDiscountFunctor = new ARM_ZeroCurveFunctor(this);
    itsFixingFunctor = new ARM_ZeroCurveFunctor(this);
	CC_ARM_SETNAME(ARM_PRICINGMODEL);

	/// by default model name is the name of the curve
	//tb
	/*if( itsZeroCurve != ARM_ZeroCurvePtr(NULL) )
		itsModelName = itsZeroCurve->GetCurrencyUnit()->GetCcyName()*/;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Constructor with market data objects
///	Returns: 
///	Action : Constructor that inialise the MDM with objects
////////////////////////////////////////////////////
ARM_PricingModel::ARM_PricingModel(const ARM_ObjectVector& marketDatas, const ARM_StringVector& mdmKeys, const ARM_ModelParams* params, const ARM_DensityFunctor* densityFct )
:	
	ARM_RootObject(),
	itsModelName( "NoName" ),
	itsPayModelName("NoName"),
	itsZeroCurve(NULL),
	itsParams( params ? static_cast<ARM_ModelParams*>(params->Clone()): NULL),
	itsDensityFunctor( densityFct ? static_cast<ARM_DensityFunctor*>(densityFct->Clone()): NULL),
	itsNumeraire(NULL), 
	itsNumMethod(NULL),
    itsDiscountFunctor(NULL),
    itsFixingFunctor(NULL),
    itsMktDataManager(NULL),
	itsMDMKeys(0),
	itsModelStateLocalVars(0),
	itsModelStateLocalVarsIsShared(false),
	itsModelStateLocalStdDevs(0),
	itsModelStateLocalStdDevsIsShared(false),
	itsNumMethodStateLocalVarsIsShared(false),
	itsModelNb(0),
	itsModelRank(0)
{
    /// set the discounting model
    itsDiscountFunctor  = new ARM_ZeroCurveFunctor(this);
    itsFixingFunctor    = new ARM_ZeroCurveFunctor(this);

    /// Save market datas in the MarketDataManager
    if(marketDatas.size() != mdmKeys.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Number of market datas input differs from keys" );

    ARM_Date asOfDate;
    size_t i;
   /* for(i=0;i<marketDatas.size();++i)
    {
        if(dynamic_cast< ARM_ZeroCurve* >(marketDatas[i]))
        {
            asOfDate = static_cast< ARM_ZeroCurve* >(marketDatas[i])->GetAsOfDate();
            break;
        }
    }*/

    itsMktDataManager = new ARM_MarketData_ManagerRep(asOfDate);

    for(i=0;i<mdmKeys.size();++i)
        itsMktDataManager->RegisterData(mdmKeys[i],marketDatas[i]);

	CC_ARM_SETNAME(ARM_PRICINGMODEL);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModel::ARM_PricingModel(const ARM_PricingModel& rhs)
:	
	ARM_RootObject(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PricingModel::~ARM_PricingModel()
{
    CleanUp();

}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PricingModel& ARM_PricingModel::operator=(const ARM_PricingModel& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PricingModel::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

	os << indent << "Model Name " << itsModelName << "\n\n";
	os << indent << "Pay Model Name " << itsPayModelName << "\n\n";

    os << indent << itsParams->toString(indent);
    
    if( itsNumeraire != ARM_NumerairePtr(NULL) )
        os << indent << itsNumeraire->toString(indent);
    else
        os << indent << "\nNumeraire : not initialised !\n";

    if(itsNumMethod != ARM_NumMethodPtr(NULL) )
        os << indent << itsNumMethod->toString(indent);
    else
        os << indent << "\nNumerical Method : not initialised !\n";

    if(itsMktDataManager)
        os << indent << itsMktDataManager->toString(indent);

	os << indent << "\nModel Nb " << itsModelNb << "\n";
    
    return os.str();
}



////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: SetModelParams
///	Returns :
///	Action  : Set model parameters
////////////////////////////////////////////////////
void ARM_PricingModel::SetModelParams(const ARM_ModelParams& params)
{
	if (itsParams != &params)
	{
		delete itsParams;
		if(!ValidateModelParams(params))
		{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		   "Bad model parameter types or wrong number of model parameters");
		}
		itsParams=(ARM_ModelParams*) params.Clone();
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: SetNumMethod
///	Returns :
///	Action  : Set numerical method of the model
////////////////////////////////////////////////////
void ARM_PricingModel::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
    if( numMethodPtr != ARM_NumMethodPtr(NULL) )
		if(		(numMethodPtr->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING && !SupportBackwardInduction() )
			||	(numMethodPtr->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING && !SupportForwardInduction() ) )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Model induction type and numerical method direction are not compatible");
    }
    itsNumMethod = numMethodPtr;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: SetMktDataManager
///	Returns :
///	Action  : Set market data manager
////////////////////////////////////////////////////
void ARM_PricingModel::SetMktDataManager(const ARM_MarketData_ManagerRep& mktDataManager,const ARM_StringVector& mdmKeys)
{
    delete itsMktDataManager;
    itsMktDataManager=(ARM_MarketData_ManagerRep*) mktDataManager.Clone();
    itsMDMKeys = mdmKeys;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Induct
///	Returns: States result of the induction
///	Action : Default implementation for
///          backward/forward induction
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PricingModel::Induct(
	ARM_PricingStatesPtr& states,
	double toTime )
{
    if(itsNumMethod == ARM_NumMethodPtr(NULL) )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Numerical method is missing in the model");
    }
	ResetDFMap();
    return itsNumMethod->Induct(*this,states,toTime);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: PricingTimeSteps
///	Returns :
///	Action  : returns PricingTimeSteps
////////////////////////////////////////////////////
std::vector<double>* ARM_PricingModel::PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos)
{
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
	std::vector<double>* ptimeSteps = discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this );
	return ptimeSteps;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: NeedLocalDiscount
///	Returns :
///	Action  : default implementation (subclasses
///           may redefine it)
////////////////////////////////////////////////////
bool ARM_PricingModel::NeedLocalDiscount() const
{
#if defined(__GP_STRICT_VALIDATION)
    if( itsNumeraire == ARM_NumerairePtr(NULL) )
    {
        ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME +
            " : numeraire is required for local discounting" );
    }
#endif

    return itsNumeraire->GetType()!= ARM_Numeraire::TerminalZc 
		&& itsNumeraire->GetType()!= ARM_Numeraire::TerminalEventZc;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: LocalDiscounts
///	Returns :
///	Action  : default implementation (subclasses
///           may redefine it)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
#if defined(__GP_STRICT_VALIDATION)
    if(    itsNumeraire == ARM_NumerairePtr(NULL)
        || itsNumeraire->GetType()!= ARM_Numeraire::TerminalZc 
		&& itsNumeraire->GetType()!= ARM_Numeraire::TerminalEventZc )
    {
		CC_Ostringstream os;
		os  << "Only " << ARM_Numeraire::NumeraireTypeTable[ARM_Numeraire::TerminalZc] << " or "
            << ARM_Numeraire::NumeraireTypeTable[ARM_Numeraire::TerminalEventZc]
			<< " supported by this model currently";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
#endif
    ARM_VectorPtr localDiscounts(NULL);  /// by default no discounting (terminal Zc numeraire !)
    return localDiscounts;
}



////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ProcessPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "numerisation"
////////////////////////////////////////////////////
void ARM_PricingModel::ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
#if defined( __GP_STRICT_VALIDATION )
	if( itsNumeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		ARM_USERNAME + ": numerisation impossible, numeraire is missing" );	
#endif
	itsNumeraire->ProcessPaidPayoffs(payModelName,payoffs,evalTime,states,*this);
	itsNumMethod->ProcessPaidPayoffs(payoffs,evalTime);
}




////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ProcessUnPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "un-numerisation"
////////////////////////////////////////////////////
void ARM_PricingModel::ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
#if defined( __GP_STRICT_VALIDATION )
	if( itsNumeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		ARM_USERNAME + ": numerisation impossible, numeraire is missing" );	
#endif
	itsNumeraire->ProcessUnPaidPayoffs(payModelName,payoffs,evalTime,states,*this);
	itsNumMethod->ProcessUnPaidPayoffs(payoffs,evalTime);
}

void ARM_PricingModel::ProcessUnPaidPayoffs(const string& payModelName, ARM_GP_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
#if defined( __GP_STRICT_VALIDATION )
	if( itsNumeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		ARM_USERNAME + ": numerisation impossible, numeraire is missing" );	
#endif
	itsNumeraire->ProcessUnPaidPayoffs(payModelName,payoffs,evalTime,states,*this);
	itsNumMethod->ProcessUnPaidPayoffs(payoffs,evalTime);
}
////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ComputestdDevMatrixVector
///	Returns : void
///	Action  : from a var covar matrix takes either the Cholesky
///				decompositon or simply the sqrt elementwise
////////////////////////////////////////////////////
void ARM_PricingModel::ComputestdDevMatrixVector( 
	ARM_MatrixVector& varCovarMatrixVector,
	ARM_MatrixVector& stdDevMatrixVector,
	bool needToCholeskyDecompose,
	bool skipFirst ) const
{
	/// second compute the std dev as the square root of the variance
	DeletePointorVector<ARM_GP_Matrix>(stdDevMatrixVector );

	stdDevMatrixVector.resize( varCovarMatrixVector.size() );
	if( needToCholeskyDecompose )
	{
		size_t i = 0;
		if( skipFirst )
		{
			stdDevMatrixVector[0] = new ARM_GP_TriangularMatrix;
			++i;
		}

		while( i<stdDevMatrixVector.size() )
		{
			stdDevMatrixVector[i] = static_cast<ARM_GP_Matrix*>( varCovarMatrixVector[i]->Clone() );

#if defined(__GP_STRICT_VALIDATION)
			if( dynamic_cast<ARM_GP_TriangularMatrix*>(stdDevMatrixVector[i] ) == NULL )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					ARM_USERNAME + ": trying to Cholesky decompose a matrix that is not !");
#endif
			( (ARM_GP_TriangularMatrix*) stdDevMatrixVector[i] )->CholeskyDecompose();
			++i;
		}
	}
	else
	{
		size_t i = 0;
		if( skipFirst )
		{
			stdDevMatrixVector[0] = static_cast<ARM_GP_Matrix*>( varCovarMatrixVector[0]->Clone() );
			++i;
		}

		while( i<stdDevMatrixVector.size() )
		{
#if defined(__GP_STRICT_VALIDATION)
			if( !varCovarMatrixVector[i] )
				ARM_THROW( ERR_INVALID_ARGUMENT, "varCovarMatrixVector[i] null!" );
#endif
			stdDevMatrixVector[i] = static_cast<ARM_GP_Matrix*>( varCovarMatrixVector[i]->Clone() );
			for(size_t j=0; j<stdDevMatrixVector[i]->GetRowsNb(); ++j )
			{
				for(size_t k=0; k<stdDevMatrixVector[i]->GetColsNb(); ++k )
					(*stdDevMatrixVector[i])(j,k) = sqrt((*stdDevMatrixVector[i])(j,k));
			}
			++i;
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_PricingModel::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	ARM_MatrixVector modelStateLocalVars;
	ARM_MatrixVector modelStateLocalStdDevs;
	
	/// computes the local variance and the std dev matrix
	ModelStateLocalVariances( timeSteps, modelStateLocalVars); 
	ComputestdDevMatrixVector( modelStateLocalVars, modelStateLocalStdDevs, NeedsToCholeskyDecomposeFactors() );

	/// set the result
	SetModelStateLocalVars(modelStateLocalVars);
	SetModelStateLocalStdDevs(modelStateLocalStdDevs);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: NumMethodsDrifts
///	Returns : void
///	Action  : computes absolute and relative drifts
////////////////////////////////////////////////////
void ARM_PricingModel::NumMethodsDrifts( const std::vector<double>& timeSteps)
{
	ARM_GP_MatrixPtr numMethodRelativeDrifts;
	ARM_GP_MatrixPtr numMethodAbsoluteDrifts;

	IntegratedLocalDrifts(
			timeSteps,
			numMethodRelativeDrifts, 
			numMethodAbsoluteDrifts);

	/// set the result
	SetNumMethodAbsoluteDrifts(numMethodAbsoluteDrifts);
	SetNumMethodRelativeDrifts(numMethodRelativeDrifts);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: NumMethodStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_PricingModel::NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps )
{
	ARM_MatrixVector numMethodStateLocalVars;
	ARM_MatrixVector numMethodStateLocalStdDevs;
	
	/// computes the local variance and the std dev matrix
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	ComputestdDevMatrixVector( numMethodStateLocalVars, numMethodStateLocalStdDevs, NeedsToCholeskyDecomposeFactors() );

	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	SetNumMethodStateLocalStdDevs(numMethodStateLocalStdDevs);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ReInit
///	Returns : ARM_PricingStatesPtr
///	Action  : ReInit is for multi-loop numerical method to init lightwise
///				the model and its numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PricingModel::ReInit()
{
#if defined(__GP_STRICT_VALIDATION)
	if( itsNumMethod == ARM_NumMethodPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method is NULL!" );
#endif
	ResetDFMap();
	return itsNumMethod->ReInit(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ReInitLoop
///	Returns : ARM_PricingStatesPtr
///	Action  : ReInit is useful for direction change pricing
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PricingModel::ReInitLoop()
{
#if defined(__GP_STRICT_VALIDATION)
	if( itsNumMethod == ARM_NumMethodPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method is NULL!" );
#endif
	ResetDFMap();
	return itsNumMethod->ReInitLoop(*this);
}



////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: NumMethodStateLocalGlobalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local and global variance and std dev 
////////////////////////////////////////////////////
void ARM_PricingModel::NumMethodStateLocalGlobalVariancesAndStdDev(
	const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& globalVariances,
	ARM_MatrixVector& localStDev,
    ARM_MatrixVector& globalStdDev ) const
{
	NumMethodStateLocalGlobalVariances( timeSteps, localVariances, globalVariances );
	/// computes the local std deviation
	ComputestdDevMatrixVector( localVariances, localStDev, NeedsToCholeskyDecomposeFactors() ); 
	/// skip the first argument for the global variances
	ComputestdDevMatrixVector( globalVariances, globalStdDev, NeedsToCholeskyDecomposeFactors(), true ); 
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routine : NumMethodStateLocalGlobalVariances
///	Returns : void
///	Action  : computes NumMethodStateLocalGlobalVariances (for the generic tree (1st and 2nd generation)
////////////////////////////////////////////////////

void ARM_PricingModel::NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& variances ) const
{
	NumMethodStateLocalVariances(timeSteps,localVariances);
	NumMethodStateGlobalVariances(timeSteps,variances);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: CreatePricerInfo
///	Returns : ARM_PricerInfo*
///	Action  : gives the pricer info corresponding to the numerical method
///             if asked to!
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_PricingModel::CreatePricerInfo() const
{
#if defined(__GP_STRICT_VALIDATION)
	if( itsNumMethod == ARM_NumMethodPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method is NULL!" );
#endif
    return itsNumMethod->CreatePricerInfo( *this );
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: CreatePricerInfo
///	Returns : ARM_PricerInfo*
///	Action  : gives the pricer info corresponding to the numerical method
///             if asked to!
////////////////////////////////////////////////////
void ARM_PricingModel::Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter)
{}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: SetZeroCurve
///	Returns : void
///	Action  : set the zero curve in a model
////////////////////////////////////////////////////
void ARM_PricingModel::SetZeroCurve( const ARM_ZeroCurvePtr& zc)
{ 
	/// not cloned because using smart pointor!
	itsZeroCurve = zc; 
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: EulerLocalDrifts
///	Returns : void
///	Action  : Computes the relative and absolute Euler local drifts
////////////////////////////////////////////////////
void ARM_PricingModel::EulerLocalDrifts(const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method EulerLocalDrifts" );
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: IntegratedLocalDrifts
///	Returns : void
///	Action  : Computes the relative and absolute Euler local drifts
////////////////////////////////////////////////////
void ARM_PricingModel::IntegratedLocalDrifts(const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method IntegratedLocalDrifts" );
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: IntegratedMarkovianDrift
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PricingModel::IntegratedMarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, const ARM_GP_VectorPtr& driftCorrection) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method IntegratedMarkovianDrift" );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PricingModel::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method MarkovianDrift" );
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: VolatilitiesAndCorrelations
///	Returns : void
///	Action  : Computes the vols, the vols derivatives and the correlation
////////////////////////////////////////////////////
void ARM_PricingModel::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method VolatilitiesAndCorrelations" );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: UpdatePDE3DCoeffs
///	Returns : void
///	Action  : Update the PDE 3D coefficients
////////////////////////////////////////////////////
void ARM_PricingModel::UpdatePDE3DCoeffs(
		size_t timeIdx,
		const ARM_PricingStatesPtr& states,
		ARM_GP_VectorPtr& qxx,
		ARM_GP_VectorPtr& qyy,
		ARM_GP_VectorPtr& qzz,
		ARM_GP_VectorPtr& qxy,
		ARM_GP_VectorPtr& qyz,
		ARM_GP_VectorPtr& qzx,
		ARM_GP_VectorPtr& px,
		ARM_GP_VectorPtr& py,
		ARM_GP_VectorPtr& pz,
		ARM_GP_VectorPtr& o,
		double lambda,
		bool IsInit
		) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented method UpdatePDE3DCoeffs" );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: NeedLocalDiscounts
///	Returns : bool
///	Action  : tells whether a model needs to compute local discount
////////////////////////////////////////////////////
bool ARM_PricingModel::NeedLocalDiscounts() const
{	return true; }


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: DefaultLibor
///	Returns : a vector of libor values
///	Action  : Default Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultLibor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{    
    /// Get Libor Zc through the fixing functor
    ARM_VectorPtr ZcStart   = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
    ARM_VectorPtr ZcEnd     = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		double libor = ((*ZcStart)[0]/(*ZcEnd)[0]-1.0)/period;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,libor) );
    }

    /// No payment convexity computed in default implementation
    int i,nbStates=ZcStart->size();
    ARM_VectorPtr values( new std::vector<double>(nbStates) );
    for(i=0;i<nbStates;++i)
        (*values)[i]=((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;

    return values;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: DefaultAnnuity
///	Returns : a vector of annuity
///	Action  : Default Annuity Computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultAnnuity(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_PricingStatesPtr& states) const
{
    size_t i, nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
    ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));

    /// Compute fixed leg values using the discounting zc functor
    size_t iFix,nbFixFlows=fixPayTimes.size();
    for(iFix=0;iFix<nbFixFlows;++iFix)
    {
        ARM_VectorPtr ZcFixFlow = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
        for(i=0;i<nbStates;++i)
            (*values)[i] += fixPayPeriods[iFix]*(*ZcFixFlow)[i];
    }

    return values;
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: DefaultAnnuityWithNominal
///	Returns : a vector of annuity
///	Action  : Default Annuity Computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultAnnuityWithNominal(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fixNominal,
		const ARM_PricingStatesPtr& states) const
{
	/// Check that the the nominal Vector size
	if( fixNominal.size() != fixPayTimes.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": The Nominal Vector Size is not equal the number of Payment times!" );
	
    size_t i, nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
    ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));

    /// Compute fixed leg values using the discounting zc functor
    size_t iFix,nbFixFlows=fixPayTimes.size();
    for(iFix=0;iFix<nbFixFlows;++iFix)
    {
        ARM_VectorPtr ZcFixFlow = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
        for(i=0;i<nbStates;++i)
            (*values)[i] += fixPayPeriods[iFix]*(*ZcFixFlow)[i]*fixNominal[iFix];
    }

    return values;
}



////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: DefaultSwapRateInPlaceWithComputedAnnuity
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
///				WARNING: need to clone the annuity as the computation
///				is done in place... I REPEAT
///				the annuity argument will be modified by the computation as it is in place
///				so to avoid side effect on the annuity (or if you want to keep the value of the annuity)
///				one needs to first clone the value of the annuity
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultSwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods,
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
		const std::vector<double>& margin,
        bool isDbleNotional,
		const ARM_VectorPtr& FixedComputedAnnuity,
		const ARM_VectorPtr& FloatComputedAnnuity,
		const ARM_PricingStatesPtr& states) const
{
    /// Compute the fixed leg price
    ARM_VectorPtr values = FixedComputedAnnuity;
    int i,nbStates=values->size();

    if( isDbleNotional && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
    {
        /// The double notional method is still an approximation but quite correct
        ARM_VectorPtr ZcFloatStart  =   GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
        ARM_VectorPtr ZcFloatEnd    =   GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatEndTime,states);

        ARM_VectorPtr floatAnnuity;
        if(margin[0] != 0.0)
        {
            for(i=0;i<nbStates;++i)
                (*values)[i] = ( (*ZcFloatStart)[i] - (*ZcFloatEnd)[i] + margin[0] * (*FloatComputedAnnuity)[i] ) / (*values)[i];
        }
        else
        {
            for(i=0;i<nbStates;++i)
                (*values)[i] = ((*ZcFloatStart)[i]-(*ZcFloatEnd)[i])/(*values)[i];
        }
    }
    else
    {
        /// We need to compute the floating leg by the forward method and no more by the double notional
        ARM_VectorPtr floatLegValues(new std::vector<double>(nbStates,0.0));
        size_t iFloat, nbFloatFlows = floatPayTimes.size();
        double marginCoef, fwdRatio;
        for(iFloat=0;iFloat<nbFloatFlows;++iFloat)
        {
            ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTimes[iFloat],states);
            ARM_VectorPtr ZcFwdEnd      = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTimes[iFloat],states);
            ARM_VectorPtr ZcPay         = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatPayTimes[iFloat],states);
            marginCoef = margin[0] * floatPayPeriods[iFloat];
            fwdRatio = floatPayPeriods[iFloat] / fwdPayPeriods[iFloat];
            for(i=0;i<nbStates;++i)
                (*floatLegValues)[i] += (fwdRatio*((*ZcFwdStart)[i] / (*ZcFwdEnd)[i] - 1.0) + marginCoef) *
                                        (*ZcPay)[i];
        }

        for(i=0;i<nbStates;++i)
            (*values)[i] = (*floatLegValues)[i] / (*values)[i];
    }
    return values;
}



////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: DefaultSwapRateInPlaceWithComputedAnnuityAndNotional
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
///				WARNING: need to clone the annuity as the computation
///				is done in place... I REPEAT
///				the annuity argument will be modified by the computation as it is in place
///				so to avoid side effect on the annuity (or if you want to keep the value of the annuity)
///				one needs to first clone the value of the annuity
///				
///				WARNING : The computed FixAnnuity must take into account the variable Nominal 
///				We have to compute the fixed Annuity
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultSwapRateInPlaceWithComputedAnnuityAndNominal(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods,
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
		const std::vector<double>& margin,
		const ARM_VectorPtr& FixedComputedAnnuity,
		const std::vector<double>& FloatNotional,
		const ARM_PricingStatesPtr& states) const
{
    /// Compute the fixed leg price
    ARM_VectorPtr values = FixedComputedAnnuity;
    int i,nbStates=values->size();


	/// We need to compute the floating leg by the forward method and no more by the double notional
	/// The function has to be optimised in case of float frequency greater than the fix freq and No variable margin
	ARM_VectorPtr floatLegValues(new std::vector<double>(nbStates,0.0));
	size_t iFloat, nbFloatFlows = floatPayTimes.size();
	double marginCoef, fwdRatio;
	for(iFloat=0;iFloat<nbFloatFlows;++iFloat)
	{
		ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTimes[iFloat],states);
		ARM_VectorPtr ZcFwdEnd      = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTimes[iFloat],states);
		ARM_VectorPtr ZcPay         = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatPayTimes[iFloat],states);
		marginCoef = margin[0] * floatPayPeriods[iFloat];
		fwdRatio = floatPayPeriods[iFloat] / fwdPayPeriods[iFloat];
		for(i=0;i<nbStates;++i)
			(*floatLegValues)[i] += (fwdRatio*((*ZcFwdStart)[i] / (*ZcFwdEnd)[i] - 1.0) + marginCoef) *	(*ZcPay)[i] * FloatNotional[iFloat];
	}
	
	for(i=0;i<nbStates;++i)
		(*values)[i] = (*floatLegValues)[i] / (*values)[i];
    
    return values;
}




////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: NPVSwap
///	Returns: a vector of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwap computation
///           using double notional method
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::DefaultNPVSwapWithComputedAnnuity(
		const string& curveName, 
		double evalTime,
		double floatStartTime,
		double floatEndTime, 
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& floatPayTimes, 
		const std::vector<double>& floatPayPeriods, 
		const std::vector<double>& margin, 
		bool isDbleNotional,
		const std::vector<double>& FixNotional, 
		const std::vector<double>& FloatNotional, 
		const ARM_GP_Matrix& strikesPerState,
		int payRec,
		const ARM_VectorPtr& floatAnnuity,
		const ARM_PricingStatesPtr& states) const
{
	/// Compute the fixed leg price
	int i, nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	std::vector<double>*  valuesNull = new std::vector<double>(nbStates);
	ARM_VectorPtr values(valuesNull);

	//// Calculation of the Fixed Leg
    ARM_VectorPtr FixedLegvalues(new std::vector<double>(nbStates,0.0));

    /// Compute fixed leg values using the discounting zc functor
    size_t iFix,nbFixFlows=fixPayTimes.size();
	if(nbFixFlows != strikesPerState.GetColsNb())
		ARM_THROW( ERR_INVALID_ARGUMENT, "Strikes Curve size should be equal to the number of Fix Flows!" );

	/// don't compute fixed leg if fixed notionals are = 0
	bool computeFixedLeg (false);
	for(iFix=0;iFix<nbFixFlows;++iFix)
	{
		if (FixNotional[iFix])
		{
			computeFixedLeg = true;
			break;
		}
	}

	if (computeFixedLeg)
	{
		for(iFix=0;iFix<nbFixFlows;++iFix)
		{
			ARM_VectorPtr ZcFixFlow = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
			for(i=0;i<nbStates;++i)
				(*FixedLegvalues)[i] += fixPayPeriods[iFix]*(*ZcFixFlow)[i]*strikesPerState(i,iFix) * (FixNotional)[iFix];
		}
	}


	///// Margin constant or Null
    if( isDbleNotional )
    {
        /// The double notional method is still an approximation but quite correct
        ARM_VectorPtr ZcFloatStart  =   GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
        ARM_VectorPtr ZcFloatEnd    =   GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatEndTime,states);

        
        if(margin[0] != 0.0)
		{
            for(i=0;i<nbStates;++i)
                (*values)[i] = payRec*( (*FixedLegvalues)[i]-( (*ZcFloatStart)[i] - (*ZcFloatEnd)[i] + margin[0] * (*floatAnnuity)[i] ) * (FloatNotional)[0]);
		}
        else
        {
            for(i=0;i<nbStates;++i)
                (*values)[i] = payRec*( (*FixedLegvalues)[i] -( (*ZcFloatStart)[i] - (*ZcFloatEnd)[i] )* (FloatNotional)[0] );
        }
    }

    else
    {
		/// Variable Margin
        /// We need to compute the floating leg by the forward method and no more by the double notional
        ARM_VectorPtr floatLegValues(new std::vector<double>(nbStates,0.0));
        size_t iFloat, nbFloatFlows = floatPayTimes.size();
        double marginCoef, fwdRatio;
        for(iFloat=0;iFloat<nbFloatFlows;++iFloat)
        {
            ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTimes[iFloat],states);
            ARM_VectorPtr ZcFwdEnd      = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTimes[iFloat],states);
            ARM_VectorPtr ZcPay         = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatPayTimes[iFloat],states);
            marginCoef = margin[iFloat] * floatPayPeriods[iFloat];
            fwdRatio = floatPayPeriods[iFloat] / fwdPayPeriods[iFloat];
            for(i=0;i<nbStates;++i)
                (*floatLegValues)[i] += (fwdRatio*((*ZcFwdStart)[i] / (*ZcFwdEnd)[i] - 1.0) + marginCoef) * (*ZcPay)[i] * (FloatNotional)[iFloat];
        }

        for(i=0;i<nbStates;++i)
            (*values)[i] = payRec*( (*FixedLegvalues)[i]-(*floatLegValues)[i]);
    }
    return values;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwapleg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PricingModel::DefaultNPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& payTimes, 
		const std::vector<double>& payPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const
{
	/// Compute the fixed leg price
	int i, nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// Variable Margin
    /// We need to compute the floating leg by the forward method and no more by the double notional
    size_t size = payTimes.size();
	ARM_GP_MatrixPtr values(new ARM_GP_Matrix(size,nbStates,0.0));
    double marginCoef, fwdRatio;
    for(size_t iFloat(0);iFloat<size;++iFloat)
    {
        ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTimes[iFloat],states);
        ARM_VectorPtr ZcFwdEnd      = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTimes[iFloat],states);
        ARM_VectorPtr ZcPay         = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTimes[iFloat],states);
        marginCoef = margin[iFloat] * payPeriods[iFloat];
        fwdRatio = payPeriods[iFloat] / fwdPayPeriods[iFloat];
        for(i=0;i<nbStates;++i)
            (*values)(iFloat,i) = (fwdRatio*((*ZcFwdStart)[i] / (*ZcFwdEnd)[i] - 1.0) + marginCoef) * (*ZcPay)[i] * (notional)[iFloat];
    }

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: DefaultNPVFixLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti)=0,K)
///	Action : 
/// Default: Default FixLeg computation 
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PricingModel::DefaultNPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const
{

	/// Compute the fixed leg price
	int i, nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	//// Calculation of the Fixed Leg

    /// Compute fixed leg values using the discounting zc functor
    size_t iFix,nbFixFlows=fixPayTimes.size();
	ARM_GP_MatrixPtr values(new ARM_GP_Matrix(nbFixFlows,nbStates,0.0));
	if(nbFixFlows != strikesPerState.GetColsNb())
		ARM_THROW( ERR_INVALID_ARGUMENT, "Strikes Curve size should be equal to the number of Fix Flows!" );

	for(iFix=0;iFix<nbFixFlows;++iFix)
	{
		ARM_VectorPtr ZcFixFlow = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
		for(i=0;i<nbStates;++i)
		{
			(*values)(iFix,i)= payRec*fixPayPeriods[iFix]*(*ZcFixFlow)[i]*strikesPerState(i,iFix) * (FixNotional)[iFix];
		}
	}
	return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: DefaultImpliedVol
///	Returns: double
///	Action : calculates implied volatility
/// by inverting B&S Formula
////////////////////////////////////////////////////
double ARM_PricingModel::DefaultImpliedVol(const ARM_VanillaArg& arg) const
{

   double target = arg.Price(const_cast<ARM_PricingModel*>(this));
   ARM_FlatSurface surface( 0.45 );
   ARM_SurfaceModelParam  surfaceModelParam(ARM_ModelParamType::Volatility, &surface, "ImpliedVol", 0.001, 2.00);
   ARM_BS_ModelParams  modelParams(ARM_ModelParamVector(1,&surfaceModelParam));
   ARM_BS_Model bsmodel(itsZeroCurve,modelParams);
   //ARM_CalibMethod calibMethod;

   /// FIX FIX should be finished after using only ARM_VanillaPortfolio
   /// to calibration
   
   switch( arg.GetType() )
	{
   case ARM_VanillaArg::VANILLA_CAP:
        {
        }
        break;
    case ARM_VanillaArg::VANILLA_SWAPTION:
        {
        }
        break;
    default:
        ARM_THROW( ERR_INVALID_ARGUMENT, "DefaultImpliedVol: Only either std caplet or std swaption is avalaible" );
		break;
   }
   

   return 0.0;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: FactorCount
///	Returns: size_t
///	Action : returns the nb of factor of the model
////////////////////////////////////////////////////

size_t ARM_PricingModel::FactorCount() const
{
#if defined(__GP_STRICT_VALIDATION)
	if( !itsParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "model params is not set!" );
#endif
	return itsParams->FactorCount();
}




////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: OffsetTimeIndexWithModelNb
///	Returns: size_t
///	Action : returns the offset for a given model nb
////////////////////////////////////////////////////

size_t ARM_PricingModel::OffsetTimeIndexWithModelNb( size_t timeIdx, size_t modelNb) const
{
	const ARM_NumMethodPtr numMethod = GetNumMethod();
#if defined(__GP_STRICT_VALIDATION)
	if( numMethod == ARM_NumMethodPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "numMethod == NULL!" );
#endif
	const ARM_GP_Vector* const timeSteps = numMethod->GetTimeSteps();
#if defined(__GP_STRICT_VALIDATION)
	if( !timeSteps )
		ARM_THROW( ERR_INVALID_ARGUMENT, "timeSteps == NULL!" );
#endif
	size_t timeStepsSize = timeSteps->size();
	
	/// warning, the individual size is timeStepsSize-1!
	return timeIdx+(modelNb*(timeStepsSize-1));
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: GetLocalMatrixElemWithModelNb
///	Returns: double
///	Action : return the element i,j of the matrixvector knowing the timeIdx and the model number
////////////////////////////////////////////////////

double ARM_PricingModel::GetLocalMatrixElemWithModelNb( const ARM_MatrixVector& matrix, size_t timeIdx, size_t stateIdx, size_t i, size_t j ) const
{
#if defined(__GP_STRICT_VALIDATION)
	size_t offsetIndex = OffsetTimeIndexWithModelNb(timeIdx, stateIdx );
	if( offsetIndex >= matrix.size() )
        ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME + " : out of range!" );
#endif
	return matrix[ OffsetTimeIndexWithModelNb(timeIdx, stateIdx ) ]->Elt(i,j);
}




////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: IsCurrentlyOnlyAnalyticalModel
///	Returns: bool
///	Action : return 1 if there is no numerical method
////////////////////////////////////////////////////

bool ARM_PricingModel::IsCurrentlyOnlyAnalyticalModel() const
{ 
	return itsNumMethod == ARM_NumMethodPtr(NULL); 
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: SetModelStateLocalVars
///	Returns: void
///	Action : set the model state local vars
////////////////////////////////////////////////////

void ARM_PricingModel::SetModelStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool shareStateLocalVars )
{
	if( !itsModelStateLocalVarsIsShared ) DeletePointorVector<ARM_GP_Matrix>(itsModelStateLocalVars); 
	itsModelStateLocalVars=stateLocalVars;
	itsModelStateLocalVarsIsShared = shareStateLocalVars;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: SetModelStateLocalStdDevs
///	Returns: void
///	Action : set the model state local vars
////////////////////////////////////////////////////

void ARM_PricingModel::SetModelStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool shareStateLocalStdDevs )
{
	if( !itsModelStateLocalStdDevsIsShared ) DeletePointorVector<ARM_GP_Matrix>(itsModelStateLocalStdDevs); 
	itsModelStateLocalStdDevs=stateLocalStdDevs;
	itsModelStateLocalStdDevsIsShared = shareStateLocalStdDevs;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: SetNumMethodLocalVars
///	Returns: void
///	Action : set the num method state local vars
////////////////////////////////////////////////////

void ARM_PricingModel::SetNumMethodStateLocalVars( const ARM_MatrixVector& stateLocalVars ,bool isShared )
{
	if(!isShared) DeletePointorVector<ARM_GP_Matrix>(itsNumMethodStateLocalVars);
	itsNumMethodStateLocalVars = stateLocalVars;
	itsNumMethodStateLocalVarsIsShared = isShared;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: SetNumMethodLocalStdDevs
///	Returns: void
///	Action : set the num method state local vars
////////////////////////////////////////////////////

void ARM_PricingModel::SetNumMethodStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool isShared )
{
	if(!isShared) DeletePointorVector<ARM_GP_Matrix>(itsNumMethodStateLocalStdDevs);
	itsNumMethodStateLocalStdDevs=stateLocalStdDevs;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: IntegratedRiskNeutralDrift
///	Returns: a vector
///	Action : computes the integrated risk neutral drift i.e.
///          Integ{s=t(timeIdx)->t(timeIdx+1), r(s)ds}
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::IntegratedRiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented method IntegratedRiskNeutralDrift()" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: RiskNeutralDrift
///	Returns: a vector
///	Action : computes the instanteneous risk neutral drift i.e.
///          r(t(timeIdx))
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModel::RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented method RiskNeutralDrift()" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: ModelFixTimeStep
///	Returns: give the number of days to use when computing a Euler Scheme
///	Action : the default is a discretisation of one year
////////////////////////////////////////////////////
int ARM_PricingModel::ModelFixTimeStep( int fixTimeStep ) const
{ 
	return K_YEAR_LEN/(fixTimeStep+1); 
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: Finalize
///	Returns: void
///	Action : finalize provided a numerical method
////////////////////////////////////////////////////
void ARM_PricingModel::Finalize()
{
	if(!itsNumMethod.IsNull())
		itsNumMethod->Finalize();

	if (!itsNumeraire.IsNull())
		GetNumeraire()->ResetNumDiscountMap();
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: GetOtherPayoffsFalg
///	Returns: bool
///	Action : does the nummethod advise to use otherpayoffs?
////////////////////////////////////////////////////
bool ARM_PricingModel::GetOtherPayoffsFlag() const 
{ 
	if( itsNumMethod == ARM_NumMethodPtr(NULL) )
		return true; 
	else 
		return itsNumMethod->GetOtherPayoffsFlag();
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: GetOtherPayoffsFalg
///	Returns: bool
///	Action : does the nummethod advise to use otherpayoffs?
////////////////////////////////////////////////////

double ARM_PricingModel::PartialDerivative( const ARM_ModelParam& modelParam, 
                         size_t number, 
                         size_t factorNb,
                         const ARM_VanillaArg& arg,
                         ARM_MktTargetType targetFuncType )
{
	if( HasClosedFormsDerivatives( modelParam.GetType(), factorNb ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, " model says that it has closed form derivatives but \
		did not implement it, pelase advise!" ); 

	ARM_THROW( ERR_INVALID_ARGUMENT, " does not make sense that you ask for closed form derivatives as \
		the model says that it does not have one!" ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: NeedMCIntegProcess
///	Returns: ARM_BoolVector
///	Action : With this function the model tells to the numerical method for each 
//  dimensions of pocess simulated if they are stored as the integrated
//  process or increments.
////////////////////////////////////////////////////

ARM_BoolVector ARM_PricingModel::NeedMCIntegProcess() const
{
	ARM_BoolVector ret(0);
	return ret;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: UnderlyingCorrelation
///	Returns: double
///	Action : Compute the correlation between instruments
///			  depending of input type
////////////////////////////////////////////////////
double ARM_PricingModel::UnderlyingCorrelation(  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const
{
	string className(typeid(*this).name());
	ARM_THROW( ERR_INVALID_ARGUMENT, "UnderlyingCorrelation not implemented for class " + className );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: UnderlyingCovariance
///	Returns: double
///	Action : Compute the correlation between instruments
///			  depending of input type
////////////////////////////////////////////////////
double ARM_PricingModel::UnderlyingCovariance (  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const
{
	string className(typeid(*this).name());
	ARM_THROW( ERR_INVALID_ARGUMENT, "UnderlyingCovariance not implemented for class " + className );
}


ARM_VectorPtr ARM_PricingModel::DiscountFactor( 
		const string& curveName,
        double evalTime, 
		const std::vector<double>&  maturityTime,
        const ARM_PricingStatesPtr& states) const
{
	return this->DiscountFactor(
								curveName, 
								evalTime, 
								maturityTime[0], 
								states);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
