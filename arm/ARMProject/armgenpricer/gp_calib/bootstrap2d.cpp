/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file bootstrap2d.cpp
 *
 *  \brief bootstrap 1 dimension model fitter
 *	\author  A. Schauly
 *	\version 1.0
 *	\date October 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/bootstrap2d.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/surface.h"

///gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/nagsolver.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"

/// gpcalib
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"

/// kernel
//#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_Bootstrap2D::ARM_Bootstrap2D(const ARM_Bootstrap2D& rhs)
:   
	ARM_ModelFitter(rhs),
	itsFunction(FunctionToSolveWithDerivativePtr(NULL)),
	itsSolver(NULL)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Bootstrap2D:: ~ARM_Bootstrap2D()
{
	CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Bootstrap2D& ARM_Bootstrap2D::operator = (const ARM_Bootstrap2D& rhs)
{
    if( this != &rhs )
	{
		CleanUp();
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: CleanUp
///	Returns: void 
///	Action : delete all pointors
////////////////////////////////////////////////////
void ARM_Bootstrap2D::CleanUp(  )
{}

////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : copies the additional member variables defined
///	in this object
////////////////////////////////////////////////////
void ARM_Bootstrap2D::CopyNoCleanUp( const ARM_Bootstrap2D& rhs )
{
    itsFunction			= FunctionToSolveWithDerivativePtr(rhs.itsFunction == FunctionToSolveWithDerivativePtr(NULL) ? NULL : new FunctionToSolveWithDerivative( *rhs.itsFunction )); 
    itsSolver			= ModifiedNRSolverPtr(rhs.itsSolver == ModifiedNRSolverPtr(NULL) ? NULL : rhs.itsSolver->Clone() );      
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_Bootstrap2D::Clone() const
{
	return new ARM_Bootstrap2D(*this); 
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_Bootstrap2D::ARM_Bootstrap2D( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    const size_t max_iter,
	bool printLevel,
	double fxTolerance,
	double xTolerance,
	const size_t dicho_max_iter,
	double dichoXTolerance,
	ARM_SolverType type,
	size_t factorNb,
	size_t NbCalibIter,
    ARM_MktTargetType  targetType)
:   
    ARM_ModelFitter( model, 
        portfolio, 
        modelParam,
        linkedModelFitter, 
        previousModelFitter,
        max_iter,
        printLevel, 
        factorNb, 
        NbCalibIter,
        targetType),
    itsFunction( NULL ),
	itsSolver( NULL )
{
	/// should we get some details!
	StoreDetailsInFileIfRequired();
    Validate();

	itsFunction = FunctionToSolveWithDerivativePtr(new FunctionToSolveWithDerivative( GetPricingModel(), this, GetCalibParam(0)->GetType(), new ARM_BinFuncMinus ));
    const double target   = 0.0;	

	switch( type )
	{
	case ARM_ModelFitterSolverType::NewtonRaphson:
		itsSolver   = ModifiedNRSolverPtr( new T_NewtonRaphsonSolver<FunctionToSolveWithDerivative>(
			*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
		break;
	case  ARM_ModelFitterSolverType::NewtonRaphsonWithRetrial:
		itsSolver   = ModifiedNRSolverPtr(new T_NewtonRaphsonSolverRetrial<FunctionToSolveWithDerivative>(
			*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
		break;
	case ARM_ModelFitterSolverType::SmoothNewthonRhaphson:
		itsSolver   = ModifiedNRSolverPtr( new T_SmoothNewtonRaphsonSolver<FunctionToSolveWithDerivative>(
			*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
		break;
	case ARM_ModelFitterSolverType::NewthonRhaphsonNoThrow:
		itsSolver   = ModifiedNRSolverPtr( new T_NewtonRaphsonSolverNoThrow<FunctionToSolveWithDerivative>(
			*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
        break;
	case ARM_ModelFitterSolverType::NewtonRaphsonWithDichotomy:
        {
		    ModifiedNRSolverPtr nrSolver= ModifiedNRSolverPtr( new T_SmoothNewtonRaphsonSolver<FunctionToSolveWithDerivative>(
			    *itsFunction,target, fxTolerance,dichoXTolerance, dicho_max_iter, printLevel ));
            itsSolver   = ModifiedNRSolverPtr( new T_DichotomySolver<FunctionToSolveWithDerivative>(
				    *itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel,NULL,nrSolver ));

		    /*itsSolver = ModifiedNRSolverPtr( new T_SolverWithDichotomy<FunctionToSolveWithDerivative,ModifiedNRSolverPtr>(
			    *itsFunction,nrSolver,target,fxTolerance,dichoXTolerance,dicho_max_iter,printLevel ));*/
        }
		break;
	case ARM_ModelFitterSolverType::Dichotomy:
		itsSolver   = ModifiedNRSolverPtr( new T_DichotomySolver<FunctionToSolveWithDerivative>(
				*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
		break;
	case ARM_ModelFitterSolverType::NagSolver:
		itsSolver   = ModifiedNRSolverPtr( new T_NagSolver<FunctionToSolveWithDerivative>(
				*itsFunction,target, fxTolerance,xTolerance, max_iter, printLevel ));
		break;
	default:
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,	"ARM_SolverType::unknown type!" );
	}
	Initialise();
};

////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D 
///	Routine: Validate
///	Returns: void
///	Action : initialise the object
////////////////////////////////////////////////////
void ARM_Bootstrap2D::Validate()
{
    /*if( GetPortfolio() == ARM_StdPortfolioPtr(NULL))
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"One dimensional bootstrap should have a non null portfolio!" );
    
	if(!GetPortfolio()->IsSameAssetsName())
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			" It is not authorised to calibarte a portofolio with differents  product Name\
            Classe. In case of need, please, use previous calibration principle or ask R&D team to help you!" );

	if(!GetPortfolio()->IsGrowingByExpiry())
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			" Only portfolio increasing by maturity is avalaible, please advise!" );
    
    if( GetCalibParams().size() != 1 )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"One dimensional bootstrap needs only one model params!" );*/
}



////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D 
///	Routine: Initialise
///	Returns: void
///	Action : initialise the object
////////////////////////////////////////////////////
void ARM_Bootstrap2D::Initialise()
{
    /// change the calibParam to the right type!
	GetPricingModel()->AdviseBreakPointTimes( GetPortfolio(), GetCalibParam(), GetFactorNb());
    GetPricingModel()->GetModelParams()->MergeModelParam(GetCalibParam(), GetFactorNb());
}

////////////////////////////////////////////////////
///	Class  : ARM_Calibrator
///	Routine: RootFinding
///	Returns: double
///	Action : RootFinding 1D
////////////////////////////////////////////////////
double ARM_Bootstrap2D::RootFinding(
	double mktPrice, 
	ARM_Security* sec,
	double guess_initial,
	double lowerbound,
	double upperbound,
	int index,
	double expiry,
	double tenor )
{
	/// then get the vanilla arg
    ARM_VanillaArg*  argVanilla = index < GetVanillaArgVector().size() ? (GetVanillaArgVector(index)!=NULL ? (ARM_VanillaArg*)GetVanillaArgVector(index)->Clone() : NULL) : NULL;
    if(argVanilla)
    {
        argVanilla->SetIndex(index);
        argVanilla->SetMktPrice(mktPrice);
        argVanilla->SetExpiry(expiry);
        argVanilla->SetTenor(tenor);
        itsFunction->SetVanillaProduct(ARM_VanillaArgPtr(argVanilla));
    }
    else 
        itsFunction->InitArgument(sec,mktPrice,index, expiry, tenor);
    itsSolver->setInitialGuess(guess_initial,lowerbound,upperbound);
    return itsSolver->Solve();    
}


////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D 
///	Routine: Initialise
///	Returns: void
///	Action : initialise the object
////////////////////////////////////////////////////
void ARM_Bootstrap2D::Calibrate()
{

	/*int size = GetPortfolio()->size();
    int start, end,increment;
	ARM_ModelFitterPtr linkedModelFitter( GetLinkedModelFitter() );

    switch( GetCalibDirection() )
    {
    case CalibDirection_Forward:
        {
            start       = 0;
            end         = size;
            increment   = 1;
        }
        break;
        
    case CalibDirection_None : case CalibDirection_Backward: default:
        {
            start       = size-1;
            end         = -1;
            increment   = -1;
        }
        break;
    }

	/// case of curve calib param
	if( ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam()) )
	{
		GetPricingModel()->Re_InitialiseCalibParams(*this);

		for(int secIndex=start; secIndex!=end; secIndex+=increment)
		{
			if(fabs( GetPortfolio()->GetWeights()->Elt(secIndex) ) > K_DOUBLE_TOL )
			{
                double root = (curveCalibParam->GetCurve()->GetOrdinates())[secIndex];
				double time = (curveCalibParam->GetCurve()->GetAbscisses())[secIndex];
		
				GetPricingModel()->AdviseCurrentCalibSecIndex(secIndex,*this);
				double mktprice		= GetPortfolio()->GetMktPrices()->Elt(secIndex); 
				double lowerbound	= (*(curveCalibParam->GetLowerBound()))[secIndex];
				double upperbound	= (*(curveCalibParam->GetUpperBound()))[secIndex];            
				
				linkedModelFitter->setStartAndEndCalibrationTimes( (secIndex == start)?0:(curveCalibParam->GetCurve()->GetAbscisses())[secIndex-increment],time);

				for( int j=0 ; j < 10 ; ++j)
				{
					if( ARM_ModelFitterPtr(NULL) == linkedModelFitter )
						linkedModelFitter->Calibrate();

					root = RootFinding(mktprice,GetPortfolio()->GetAsset(secIndex), root,lowerbound, upperbound, secIndex, time);
					if( root<lowerbound )
					{
						AddWarningMessage( "Bootstrap solution is lower than lowerbound, set the solution to the lower bound\n" );
						root = lowerbound;
					}
					else if( root > upperbound )
					{
						AddWarningMessage( "Bootstrap solution is upper than upperbound, set the solution to the upper bound\n" );
						root = upperbound;
					}
					if (	curveCalibParam->GetType() == ARM_ModelParamType::Volatility
						||  curveCalibParam->GetType() == ARM_ModelParamType::QVol )
						root = fabs(root);

					curveCalibParam->SetValueAtPoint(secIndex,root );
					GetPricingModel()->GetModelParams()->GetModelParam(curveCalibParam->GetType(),GetFactorNb()).SetValue(time,root);
				}
			}
		}
		
		GetPricingModel()->GetModelParams()->MergeModelParam(GetCalibParam(),GetFactorNb());
	}

	/// case of surface calib param
	else if( ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(GetCalibParam()) )
	{
		ARM_SurfaceWithInterpol* surface = dynamic_cast<ARM_SurfaceWithInterpol*>(surfaceCalibParam->GetSurface());
		if( !surface )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"Expected a discretised surface with an interpolator!" );
		
		GetPricingModel()->Re_InitialiseCalibParams(*this);
		const std::vector<double>& expiries = surface->GetX1();
		const std::vector<double>& tenors	= surface->GetX2();
		const ARM_GP_Matrix& values	= surface->GetX3();
		size_t secIndex = 0;

		for(size_t i=0; i<expiries.size(); ++i )
		{
			for( size_t j=0; j<tenors.size(); ++j )
			{
				
#ifdef __GP_STRICT_VALIDATION
				if(secIndex >= GetPortfolio()->GetSize())
					throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "out of bound on the portfolio!" );
#endif
				double root		= values(i,j);
				double expiry	= expiries[i];
				double tenor	= tenors[j];
				
				if(fabs( GetPortfolio()->GetWeights()->Elt(secIndex) ) > K_DOUBLE_TOL )
				{
					GetPricingModel()->AdviseCurrentCalibSecIndex(secIndex,*this);
					double mktprice		= GetPortfolio()->GetMktPrices()->Elt(secIndex); 
					double lowerbound	= surfaceCalibParam->GetLowerBound();
					double upperbound	= surfaceCalibParam->GetUpperBound();
					root = RootFinding(mktprice,GetPortfolio()->GetAsset(secIndex), root,lowerbound, upperbound, secIndex, expiry,tenor);
					if( root<lowerbound )
					{
						std::ostringstream os;
						os << "Warning: instrument " << i << " bootstrap lower than lowerbound, solution set to lowerbound\n";
						AddWarningMessage( os.str() );
						root = lowerbound;
					}
					else if( root > upperbound )
					{
						std::ostringstream os;
						os << "Warning: instrument " << i << " bootstrap upper than upperbound, solution set to upperbound\n";
						AddWarningMessage( os.str() );
						root = upperbound;
					}
					if(		ARM_ModelParamType::Volatility == surfaceCalibParam->GetType() 
						||	ARM_ModelParamType::QVol == surfaceCalibParam->GetType()  )
						root = fabs(root);

					surfaceCalibParam->SetValueAtPoint(i,j,root);
					GetPricingModel()->GetModelParams()->GetModelParam(surfaceCalibParam->GetType(),GetFactorNb()).SetValue(expiry,tenor,root);
				}
				++secIndex;
			}
		}
		GetPricingModel()->GetModelParams()->MergeModelParam(GetCalibParam(),GetFactorNb());
	}
	else
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "type not supported!" );
*/
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: detailedString
///	Returns: 
///	Action : when detail mode is on, helps to give more
///				details
////////////////////////////////////////////////////
string ARM_Bootstrap2D::detailedString( const string& indent ) const
{
		/*if(GetDetails()&&GetFileName())
		{
			CC_NS(std,ofstream) file(GetFileName()->c_str());
			file<< itsSolver->GetSolverDetails()->GetOs().str()<<CC_NS(std,endl);
			file.close();
		}*/
	return ARM_ModelFitter::detailedString(indent);
}

////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap2D
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Bootstrap2D::toString( const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
   
    size_t i,j;
    os << indent << "Bootstrap_1d Model Fitter \n";
    string tmpIndent = indent + "|\t ";
    os << tmpIndent << "Calib parameter     : "<< GetCalibParam(0)->GetTypeString() 
		<< " (" << GetCalibParam(0)->GetTypeString() << ")\n";		
    os << tmpIndent << "Corresponding errors:\n";
    os << tmpIndent << "\n";

    /// for unix compatibility
    /// we use printf
    /// more traditional C++ should use std::iostream manipulation
    char temp2[30];
    SetUpError();
   
    for(i=0; i<GetCalibParam()->rows(); ++i)
    {
        os << tmpIndent;
        for(j=0; j<GetCalibParam()->cols(); ++j)
        {
            CC_Ostringstream ostmp;
            if(GetCalibParam()->cols() == 0)
                ostmp <<CC_NS(std,right)<<"["<<i+1<<"]";
            else
                ostmp <<CC_NS(std,right)<<"["<<i+1<<","<<j+1<<"]";

            os<<CC_NS(std,setw)(7)<<ostmp.str();
            os<<CC_NS(std,setw)(2)<<CC_NS(std,left)<<"= ";
            sprintf( temp2, "%+6.3E (=%+4.4E %s)", (*GetError())(i+j,0), (*GetError())(i+j,1)*100.0, "%" );
            os <<temp2 << "  ";
        }
       os <<"\n";
    }
	
	os << tmpIndent << "\n";    
	os << indent << "======> End Of Infos <==================\n";
    return os.str();
}

CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

