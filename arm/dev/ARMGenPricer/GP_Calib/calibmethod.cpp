/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file calibmethod.cpp
 *
 *  \brief file for the object to describe the generic calibration
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/calibmethod.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/curve.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/modelnamemap.h"

/// gpcalib
#include "gpcalib/argconvdefault.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/modelfitter.h"


/// gpnumlib
#include "gpnumlib/solver.h"


/// gpmodel
#include "gpmodels/MultiAssets.h"

/// STL
#include <algorithm>
#include <functional>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	simple predictate to sort the Portfolio
////////////////////////////////////////////////////
typedef pair < double, ARM_CLASS_NAME> porfoliopair;
typedef vector <porfoliopair > porfoliopairVector;   

struct Functionporfoliopair : CC_NS( std, binary_function )< porfoliopair, porfoliopair, bool>
{  
	bool operator()( const porfoliopair& lhs, const porfoliopair& rhs) const
    {
        return lhs.first < rhs.first;
    }
}; 

////////////////////////////////////////////////////
///	simple predictate to find a pair in portfolio
////////////////////////////////////////////////////
struct FindSecurityWMktPriceAndClassNameUnaryVersion : public CC_NS( std, unary_function )< const porfoliopair&,bool>
{
	FindSecurityWMktPriceAndClassNameUnaryVersion( porfoliopair& unpair )
	:	itsPair( &unpair ) 
	{}
    
	bool operator()(const porfoliopair& unpair ) const
    {
        /// equality is only based on the MktPrice and Class_Name!
        return (fabs(unpair.first-itsPair->first) <ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE  && unpair.second == itsPair->second);
    }
private:
	porfoliopair* itsPair;
};

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Constuctor
///	Returns:    
///	Action : Constuctor
////////////////////////////////////////////////////

ARM_CalibMethod::ARM_CalibMethod( ARM_StdPortfolioPtr portfolio,
    const ARM_ModelParamVector&  calibParams, 
    ARM_MethodType methodType, 
    const size_t max_iter,
    ARM_MktTargetType targetFuncType,
    ARM_CalibMethod* linkedMethod,
    ARM_CalibMethod* previousMethod,
    bool isCalibMethodShared,
	size_t factorNb,
	size_t nbIter,
	bool validate,
	const ARM_DateStripPtr&					numSchedule,
	const ARM_VanillaSecDensityPtrVector&	numSecDensities)
:
	itsPortfolio(portfolio), 
    itsTargetFuncType( targetFuncType ),
    itsLinkedMethod(linkedMethod),
    itsPreviousMethod(previousMethod),
	itsNextMethod(NULL),
    itsCalibParams(),
	itsMethodType(methodType),
	itsModelFitterDes( NULL),
	itsModelFitter(NULL),
    itsIsCalibMethodShared(isCalibMethodShared),
	itsDoesValidation( validate ),
	itsFactorNb(factorNb ),
	itsNbIteration(nbIter),
	itsNumSchedule(numSchedule),
	itsNumSecDensities(numSecDensities)
{
    CC_ARM_SETNAME(ARM_CALIBMETHOD);

	switch(itsMethodType)
    {
	case ARM_CalibMethodType::Bootstrap1D:
		{
			itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterSolverType::SmoothNewthonRhaphson,max_iter,1.0e-6,SolverConstant::DefaultXTolerance);
		}
		break;
		
	case ARM_CalibMethodType::Optimize:
	case ARM_CalibMethodType::Optimize1D:
	case ARM_CalibMethodType::BootstrapND:
		{
			itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterOptimizerType::bounds_no_deriv,max_iter,1.0e-3);
		}
		break;
	case ARM_CalibMethodType::Numerical:
		{ /// FIXME ; enhancement needed
			/// itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterOptimizerType::bounds_no_deriv,max_iter,1.0e-3);
			itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterSolverType::SmoothNewthonRhaphson,max_iter,1.0e-6,SolverConstant::DefaultXTolerance);
		}
		break;
	case ARM_CalibMethodType::HW2FOnly:
		{ /// FIXME ; enhancement needed
			/// itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterOptimizerType::bounds_no_deriv,max_iter,1.0e-3);
			itsModelFitterDes = new ARM_ModelFitterDes(ARM_ModelFitterSolverType::SmoothNewthonRhaphson,max_iter,1.0e-6,SolverConstant::DefaultXTolerance);
		}
		break;
	
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"Calib Method with no calib method type! Cannot create model fitter descriptor" );
	}

    DuplicateCloneablePointorVectorInPlace<ARM_ModelParam>(calibParams,itsCalibParams);

    Validate();

    /// uses the same linked method to avoid extra copy cost!
    if(!itsIsCalibMethodShared)
    {
		itsPortfolio		= ARM_StdPortfolioPtr(portfolio!=ARM_StdPortfolioPtr(NULL) ?  (ARM_StdPortfolio*)portfolio->Clone() : NULL);
        itsLinkedMethod		= linkedMethod ? (ARM_CalibMethod*)linkedMethod->Clone() : NULL;
        itsPreviousMethod	= previousMethod ? (ARM_CalibMethod*)previousMethod->Clone() : NULL; 
    }
    itsCalibDirection = DefaultCalibDirection();
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: New Constuctor with descriptor
///	Returns:    
///	Action : Constuctor
////////////////////////////////////////////////////

ARM_CalibMethod::ARM_CalibMethod( ARM_StdPortfolioPtr portfolio,
    const ARM_ModelParamVector&  calibParams,
	ARM_MethodType methodType,
	ARM_ModelFitterDes* modelFitterDes,
    ARM_MktTargetType targetFuncType,
    ARM_CalibMethod* linkedMethod,
    ARM_CalibMethod* previousMethod,
    bool isCalibMethodShared,
	size_t FactorNb,
	size_t nbIter,
	bool validate,
	const ARM_DateStripPtr&					numSchedule,
	const ARM_VanillaSecDensityPtrVector&	numSecDensities)
:
	itsPortfolio(portfolio), 
	itsModelFitterDes( modelFitterDes? (ARM_ModelFitterDes*) modelFitterDes->Clone(): NULL ),
    itsTargetFuncType(targetFuncType ),
    itsLinkedMethod(linkedMethod),
    itsPreviousMethod(previousMethod),
	itsNextMethod(NULL),
    itsCalibParams(),
	itsMethodType(methodType),
    itsModelFitter(),
    itsIsCalibMethodShared(isCalibMethodShared),
	itsDoesValidation(validate),
	itsFactorNb(FactorNb),
	itsNbIteration(nbIter),
	itsNumSchedule(numSchedule),
	itsNumSecDensities(numSecDensities)
{
    CC_ARM_SETNAME(ARM_CALIBMETHOD);
    DuplicateCloneablePointorVectorInPlace<ARM_ModelParam>(calibParams,itsCalibParams);

    Validate();

    /// uses the same linked method to avoid extra copy cost!
    if(!itsIsCalibMethodShared)
    {
		itsPortfolio		= ARM_StdPortfolioPtr(portfolio!=ARM_StdPortfolioPtr(NULL) ?  (ARM_StdPortfolio*)portfolio->Clone() : NULL);
        itsLinkedMethod		= linkedMethod ? (ARM_CalibMethod*)linkedMethod->Clone() : NULL;
        itsPreviousMethod	= previousMethod ? (ARM_CalibMethod*)previousMethod->Clone() : NULL; 
    }
    itsCalibDirection = DefaultCalibDirection();
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Constuctor
///	Returns: 
///	Action : copy Constuctor
////////////////////////////////////////////////////

ARM_CalibMethod::ARM_CalibMethod(const ARM_CalibMethod& rhs)
: ARM_RootObject(rhs)
{
    CopyNoCleanUp(rhs);
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_CalibMethod:: ~ARM_CalibMethod()
{
    CleanUp();
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Assignment operator
///	Returns: 
///	Action : operator =
////////////////////////////////////////////////////
ARM_CalibMethod& ARM_CalibMethod::operator = (const ARM_CalibMethod& rhs)
{
    if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );

		/// Delete the old ones.
		CleanUp();
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : copies the additional member variables defined
///	in this object
////////////////////////////////////////////////////
void ARM_CalibMethod::CopyNoCleanUp( const ARM_CalibMethod& rhs )
{
    itsTargetFuncType	    = rhs.itsTargetFuncType;
	itsMethodType           = rhs.itsMethodType;
    itsCalibDirection	    = rhs.itsCalibDirection;    
    itsModelFitter		    = rhs.itsModelFitter;
	itsModelFitterDes		= rhs.itsModelFitterDes? (ARM_ModelFitterDes*) rhs.itsModelFitterDes->Clone(): NULL;
    itsIsCalibMethodShared  = rhs.itsIsCalibMethodShared;

	itsPortfolio			= itsIsCalibMethodShared ? rhs.itsPortfolio :
								ARM_StdPortfolioPtr(rhs.itsPortfolio!=ARM_StdPortfolioPtr(NULL) ?  (ARM_StdPortfolio*)rhs.itsPortfolio->Clone() : NULL);
    itsLinkedMethod         = itsIsCalibMethodShared ? rhs.itsLinkedMethod :
                                                       (rhs.itsLinkedMethod ? (ARM_CalibMethod*)rhs.itsLinkedMethod->Clone() :NULL); 
    itsPreviousMethod       = itsIsCalibMethodShared ? rhs.itsPreviousMethod :
                                                       (rhs.itsPreviousMethod ? (ARM_CalibMethod*)rhs.itsPreviousMethod->Clone() :NULL); 
    itsNextMethod       = itsIsCalibMethodShared ? rhs.itsNextMethod :
                                                       (rhs.itsNextMethod ? (ARM_CalibMethod*)rhs.itsNextMethod->Clone() :NULL); 
	itsCalibParams		    = DuplicateCloneablePointorVector<ARM_ModelParam>(rhs.itsCalibParams);
	itsDoesValidation		= rhs.itsDoesValidation;
	itsFactorNb				= rhs.itsFactorNb;
	itsNbIteration			= rhs.itsNbIteration;
    itsTargetFuncType       = rhs.itsTargetFuncType;
	
	itsNumSchedule			= CreateClonedPtr(&*rhs.itsNumSchedule);
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsNumSecDensities, itsNumSecDensities);
}


////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: CleanUp
///	Returns: 
///	Action : deletes any dynamic object defined at the object
///				level
////////////////////////////////////////////////////
void ARM_CalibMethod::CleanUp( )
{
    if( !itsIsCalibMethodShared )
    {
        delete itsLinkedMethod;
        itsLinkedMethod = NULL;
        delete itsPreviousMethod;
        itsPreviousMethod = NULL;
		delete itsNextMethod;
        itsNextMethod = NULL;
    }

    DeletePointorVector<ARM_ModelParam>(itsCalibParams);  
    delete itsModelFitterDes;
	itsModelFitterDes = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_CalibMethod::Clone() const
{
	return new ARM_CalibMethod( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Validate
///	Returns: 
///	Action : The goal of validate is to check that
///             1) we do not optimize on overlapping model params 
///             2) we have different nested calibration portfolios
///             3) there is no cycle in linked calibration
////////////////////////////////////////////////////

void ARM_CalibMethod::Validate() const
{
	/// checks that if we use limped Calibration
	/// we do not have overlapping model params
	if(itsDoesValidation&&itsLinkedMethod)
	{
		ARM_CalibMethod* nextLinkedMethod = itsLinkedMethod;
		if(itsPortfolio != ARM_StdPortfolioPtr(NULL))
		{
			/// 1) checks that the portfolios are different 
			/// ON MKT PRICES and CLASS_NAME pair by pair 
			size_t i;
			size_t sizepf1 = itsPortfolio->size();
			vector <porfoliopair> pf1(sizepf1);
			for(i =0; i< sizepf1; ++i)
			{
				double mktprice     = (*itsPortfolio->GetMktPrices())[i];
				ARM_CLASS_NAME name = itsPortfolio->GetAsset(i)->GetName();
				porfoliopair pairpf(mktprice, name);
				pf1[i] = pairpf;
			}
			/// sort the model calls for easy lookup of the corresponding dates!
			CC_NS( std, sort )( pf1.begin(), pf1.end(), Functionporfoliopair() );
			porfoliopairVector::iterator 
				searchBegin = pf1.begin(),
				searchEnd   = pf1.end();
			
			while( nextLinkedMethod && nextLinkedMethod->GetPortfolio() != ARM_StdPortfolioPtr(NULL))
			{
				size_t sizepf2 = nextLinkedMethod->GetPortfolio()->size();
				vector <pair<double, ARM_CLASS_NAME> > pf2(sizepf2);
				
				for(i =0; i< sizepf2; ++i)
				{
					double mktprice     = (*nextLinkedMethod->GetPortfolio()->GetMktPrices())[i];
					ARM_CLASS_NAME name = nextLinkedMethod->GetPortfolio()->GetAsset(i)->GetName();
					porfoliopair pairpf(mktprice, name);
					pf2[i] = pairpf;
				}
				/// sort the model calls for easy lookup of the corresponding dates!
				CC_NS( std, sort )( pf2.begin(), pf2.end(), Functionporfoliopair() );
				
				for(i=0; i<sizepf2; ++i)
				{
					vector <pair<double, ARM_CLASS_NAME> >::const_iterator found = CC_NS( std, find_if) ( searchBegin, searchEnd, 
						FindSecurityWMktPriceAndClassNameUnaryVersion( pf2[i]) );
					if( found != searchEnd )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						" There is at least a product in doubloon, please, check Pfs Construction " );
				}
				nextLinkedMethod = nextLinkedMethod->GetlinkedMethod();
			}
		}
		
		/// 2) checks that there is no model param overlap between linked method
		nextLinkedMethod = itsLinkedMethod;
		while( nextLinkedMethod)
		{
			ARM_ModelParamVector::iterator 
				searchBegin = nextLinkedMethod->GetCalibParams().begin(),
				searchEnd   = nextLinkedMethod->GetCalibParams().end();
			
			/// loop over the itsCalibParams and checks that they do not exist previously
			/// no need to clone and make superflous copies!
			/// warning the comparison is based on the model param type... and if the model param type
			/// are the same, we compare model name!
			/// we may need to test/revisit this in the case of multiple model params with same type
			
			size_t i;
			for( i=0; i<itsCalibParams.size(); ++i)
			{   
				if(itsCalibParams[i] && (*searchBegin))
				{
					ARM_ModelParamVector::iterator found = CC_NS( std, find_if) ( searchBegin, searchEnd, 
						CompareModelParamWEnumUnaryLeftVersion( itsCalibParams[i] ) );
					if( (found != searchEnd) && (nextLinkedMethod->GetFactorNb() == GetFactorNb()) )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"Overlap in calibration not allowed in the case of optimize calibration!" );
				}
			}
			/// go to the next linked method
			nextLinkedMethod = nextLinkedMethod->GetlinkedMethod();
		}
	}
}
////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Default Validate With Model
///	Returns: 
///	Action : The goal of validate is to check that
///             1) 1) tests that there is no nested cycle 
///				   2) checks that the portfolios are different ON MKT PRICES
////////////////////////////////////////////////////

void ARM_CalibMethod::DefaultValidateWithModel(const ARM_PricingModel& model) const
{
	if(itsDoesValidation && itsPreviousMethod)
	{		
		/// 2) checks that the portfolios are different ON MKT PRICES 
		ARM_CalibMethod* nextPreviousMethod = itsPreviousMethod;
		size_t i;
		size_t sizepf1 = itsPortfolio->size();
		vector <porfoliopair> pf1(sizepf1);
		for(i =0; i< sizepf1; ++i)
		{
			double mktprice     = (*itsPortfolio->GetMktPrices())[i];
			ARM_CLASS_NAME name = itsPortfolio->GetAsset(i)->GetName();
			porfoliopair pairpf(mktprice, name);
			pf1[i] = pairpf;
		}
		/// sort the model calls for easy lookup of the corresponding dates!
		CC_NS( std, sort )( pf1.begin(), pf1.end(), Functionporfoliopair() );
		porfoliopairVector::iterator 
			searchBegin = pf1.begin(),
			searchEnd   = pf1.end();
		while( nextPreviousMethod)
		{
			size_t sizepf2 = nextPreviousMethod->GetPortfolio()->size();
			vector <pair<double, ARM_CLASS_NAME> > pf2(sizepf2);
			
			for(i =0; i< sizepf2; ++i)
			{
				double mktprice     = (*nextPreviousMethod->GetPortfolio()->GetMktPrices())[i];
				ARM_CLASS_NAME name = nextPreviousMethod->GetPortfolio()->GetAsset(i)->GetName();
				porfoliopair pairpf(mktprice, name);
				pf2[i] = pairpf;
			}
			/// sort the model calls for easy lookup of the corresponding dates!
			CC_NS( std, sort )( pf2.begin(), pf2.end(), Functionporfoliopair() );
			
			for(i=0; i<sizepf2; ++i)
			{
				vector <pair<double, ARM_CLASS_NAME> >::const_iterator found = CC_NS( std, find_if) ( searchBegin, searchEnd, 
					FindSecurityWMktPriceAndClassNameUnaryVersion( pf2[i]) );
				if( found != searchEnd )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					" There is at least a product in doubloon, please, check Pfs Construction " );
			}
			nextPreviousMethod = nextPreviousMethod->GetPreviousMethod();
		}
	}
	
	if(GetlinkedMethod())
	{	           
		ARM_CalibMethod*  nextLinkedMethod = GetlinkedMethod();
		nextLinkedMethod->DefaultValidateWithModel(model);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Initialise
///	Returns: 
///	Action : Initialise the calibration schudle
///          and models fitter!
////////////////////////////////////////////////////

void ARM_CalibMethod::Initialise(ARM_PricingModel* model)
{
    if(itsLinkedMethod)
        itsLinkedMethod->Initialise(model);
    if(itsPreviousMethod)
        itsPreviousMethod->Initialise(model);
    
    InitialiseModelFitter( model );
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: ComputeCalibDirection
///	Returns: 
///	Action : computes the default calib direction ... model 
///     can overwritte this if necessary!
////////////////////////////////////////////////////
ARM_CalibDirection ARM_CalibMethod::DefaultCalibDirection() const
{
    switch( itsMethodType )
    {
    case ARM_CalibMethodType::Bootstrap1D: case ARM_CalibMethodType::BootstrapND:
        return CalibDirection_Forward;
    case  ARM_CalibMethodType::Optimize : case ARM_CalibMethodType::Unknown : case ARM_CalibMethodType::Optimize1D: case ARM_CalibMethodType::Numerical:case ARM_CalibMethodType::HW2FOnly:
        return CalibDirection_None;
    default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknonw calib method type" );
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: InitialiseModelFitter
///	Returns: 
///	Action : This initialises the corresponding model fitter
/// allowing for nested ARM_ModelFitter*
////////////////////////////////////////////////////
void ARM_CalibMethod::InitialiseModelFitter( ARM_PricingModel* model )
{
    ARM_ModelFitterPtr linkedmodelfitter = itsLinkedMethod ? itsLinkedMethod->GetModelFitter() : ARM_ModelFitterPtr(NULL);
	
	/// linked model fitter does not store current warning as this is only in the last loop that this is done
	if( linkedmodelfitter != ARM_ModelFitterPtr (NULL ) )
		linkedmodelfitter->SetKeepCurrentWarning( false );

    ARM_ModelFitterPtr previousmodelfitter = itsPreviousMethod ? itsPreviousMethod->GetModelFitter() : ARM_ModelFitterPtr(NULL);

	itsModelFitter = itsModelFitterDes->CreateModelFitter(model,
		itsPortfolio,
		itsCalibParams,
		itsMethodType,
		linkedmodelfitter,
		previousmodelfitter,
		itsTargetFuncType,
		itsFactorNb,
		itsNbIteration,
		&*itsNumSchedule,
		itsNumSecDensities);

	/// preprocessing of the model fitter
	model->PreProcessing( *itsModelFitter );
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the individual calibration method
////////////////////////////////////////////////////
void ARM_CalibMethod::Calibrate( ARM_PricingModel* model )
{
    ///Does some validations with model
	model->ValidateCalibMethod(*this);
	
	/// initialize calibmethod
    Initialise(model);

	/// does the calibration process
	itsModelFitter->DoCalibrationProcess();

	/// post processing of all the methods
	itsModelFitter->PostProcessing();

	/// set the refModel for multi-asset
	ARM_MultiAssetsModel* multiAssets = dynamic_cast<ARM_MultiAssetsModel*>(model);
	if( multiAssets )
	{
		ARM_ModelNameMap* modelMap = multiAssets->GetModelMap();
		multiAssets->SetRefModel( &*(*modelMap)[ multiAssets->GetModelNamePerFactor( itsFactorNb ) ]->Model() );
	}
	if(GetNextMethod())
		GetNextMethod()->Calibrate(model);
}



////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: toString
///	Returns: string
///	Action : browse the object and stringify it
////////////////////////////////////////////////////
string ARM_CalibMethod::toString( const string& indent, const string& nextIndent) const
{
	return toString( indent, true );
}


////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: toString
///	Returns: string
///	Action : browse the object and stringify it
////////////////////////////////////////////////////

string ARM_CalibMethod::toString( const string& indent, bool printPreviousMethods ) const
{
    /// standard string
    CC_Ostringstream os;
    os << indent << "\n";
    
	if( GetWarning() )
	{
		os << indent <<" #    #    ##    #####   #    #     #    #    #   #### \n";
		os << indent <<" #    #   #  #   #    #  ##   #     #    ##   #  #    #\n";
		os << indent <<" #    #  #    #  #    #  # #  #     #    # #  #  #	   \n";
		os << indent <<" # ## #  ######  #####   #  # #     #    #  # #  #  ###\n";
		os << indent <<" ##  ##  #    #  #   #   #   ##     #    #   ##  #    #\n";
		os << indent <<" #    #  #    #  #    #  #    #     #    #    #   #### \n";
		os << GetWarning()->toString(indent );
	}


    os << indent << "Method type \t: " << ARM_ArgConvReverse_CalibMethod.GetString( itsMethodType ) << "\n";
	os << indent << itsModelFitterDes->toString(indent)<<"\n";
	if(itsModelFitterDes->GetGetDetails())
	    os << indent << "Detailed mode   : On\n";

	size_t i;

	if (itsCalibParams.size() > 0)
	{
		/// calib params
		os << indent << "Calib params \t: ";
		os << 1 <<"- " << itsCalibParams[0]->GetTypeString() << " (" << itsCalibParams[0]->GetTypeString() << ") (Factor=" << itsFactorNb <<", Nb of Reiteration= "<< itsNbIteration <<")\n";

		/// handle other calib Params!
		string tmpIndent = indent + "  \t\t";
		for(i=1;i<itsCalibParams.size();++i)
		{
			os << tmpIndent << i+1<<"- " << itsCalibParams[i]->GetTypeString() << " (" << itsCalibParams[i]->GetTypeString() << ")\n";
		}
		os << indent << "\n";
	}

	if( printPreviousMethods )
	{
		/// browse and print the previous method(s)
		ARM_CalibMethod* previousCalibMethod = itsPreviousMethod;
		ARM_GP_T_Vector<ARM_CalibMethod*> previousCalibMethods(0);

		/// first determine if there are some previous methods and stores them backward
		if( previousCalibMethod )
		{
			previousCalibMethods.push_back(previousCalibMethod);
			while( previousCalibMethod->GetPreviousMethod() != NULL )
			{
				previousCalibMethod = previousCalibMethod->GetPreviousMethod();
				previousCalibMethods.push_back( previousCalibMethod);
			}
		}

		/// print them
		int j;
		for( j=previousCalibMethods.size()-1, i=0; j>=0; --j, ++i )
		{
			os << indent << "\n";
			os << indent << i+1 << ") Previous Method\n";
			os << previousCalibMethods[j]->toString( indent, false );
			os << indent << "\n";
			os << indent << "-------- End Of Previous Method (" << i+1 << ") <-----------\n";
		}
	}
    os << indent << "\n";

    /// browse and print the model fitter
    if( itsModelFitter != ARM_ModelFitterPtr(NULL))
    {
        os << itsModelFitter->toString( indent);
        os << indent << "\n";
    }

	/// browse and print the linked method
    if( itsLinkedMethod)
    {
        os << indent << "Linked Method(s):\n";
        os << itsLinkedMethod->toString( indent + "| \t" );
        os << indent << "|\n";
        os << indent << "======> End Of Linked Method(s) <==================\n";
        os << indent << "\n";
    }

    if( itsModelFitter != ARM_ModelFitterPtr(NULL) && itsModelFitterDes->GetGetDetails())
    {
        os << itsModelFitter->detailedString(indent);
        os << indent << "|\n"; 
        os << indent <<"======> End Of Model Fitter Details <==================\n";
        os << indent << "\n"; 
    }
    if(GetDuration() != ARM_Timer::NOT_CALC_DURATION)
        os << indent << "Time spent      : " << GetDuration() << " second\n";
    else
        os << indent << " Time spent : NOT_CALC_DURATION\n";
    
	if (itsNextMethod!=NULL)
	{
		os << indent << "\n";
		os << indent << "Next Method(s):\n";
        os << indent << "|\n";
		os << itsNextMethod->toString( indent);
	}
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: View
///	Returns: void
///	Action : standard view
////////////////////////////////////////////////////
void ARM_CalibMethod::View(char* id, FILE* ficOut) const
{
    FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;   

    fprintf(fOut, "\n GLOBAL CALIB METHOD OBJECT \n" );
    fprintf(fOut, "%s", toString( " | \t" ).c_str() ); 
    fprintf(fOut, " | \n" );

    fprintf(fOut, " ======> END OF INFOS CALIB METHOD <================== \n\n" );
    
    if ( ficOut == NULL )
       fclose(fOut);
}


////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: View
///	Returns: void
///	Action : standard view
////////////////////////////////////////////////////
double ARM_CalibMethod::GetDuration() const
{
	if( itsModelFitter != ARM_ModelFitterPtr(NULL) )
		return itsModelFitter->GetDuration();
	else
		return ARM_Timer::NOT_CALC_DURATION;
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibMethod
///	Routine: GetWarning
///	Action : returns the warning of the calibration
////////////////////////////////////////////////////
ARM_Warning* ARM_CalibMethod::GetWarning() const
{
	if( itsModelFitter != ARM_ModelFitterPtr(NULL) )
		return itsModelFitter->GetWarning();
	return 0;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
