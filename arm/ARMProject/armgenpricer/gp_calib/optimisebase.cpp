/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 *  \file optimisebase.cpp
 *  \brief optimisebase provides the skeleton for all optimiser function
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/optimisebase.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/surface.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/surfacetypedef.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/modelparams.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"

/// gpcalib
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"

/// gpnumlib
#include "gpnumlib/nagsolver.h"

/// kernel
//#include <inst/portfolio.h>

/// STL
#include <ctime>
#include <cstdio>
#include <fstream>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: static variables
////////////////////////////////////////////////////
const bool ARM_OptimiseBase::DefaultGetDetails	= false;
const double ARM_OptimiseBase::DefaultTolerance	= 1.0e-3;
const double ARM_OptimiseBase::DefaultStepMax	= 1.0e+002;
const bool ARM_OptimiseBase::LocalSearch		= false;
const double ARM_OptimiseBase::DefaultRelativeTolerance = 0.01;
const double ARM_OptimiseBase::DefaultAbsoluteTolerance = 0.0001;

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_OptimiseBase::ARM_OptimiseBase(const ARM_OptimiseBase& rhs)
:	
	ARM_ModelFitter(rhs), 
	itsTolerance(rhs.itsTolerance),
	itsStepMax(rhs.itsStepMax),
	itsLocalSearch(rhs.itsLocalSearch),
	itsInitialValuesVec( rhs.itsInitialValuesVec ? (ARM_GP_Vector*) rhs.itsInitialValuesVec->Clone(): NULL ),
	itsInitialValuesMat( rhs.itsInitialValuesMat ? (ARM_GP_Matrix*) rhs.itsInitialValuesMat->Clone(): NULL ),
	itsFunction( rhs.itsFunction? rhs.itsFunction->Clone() : NULL )
{}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_OptimiseBase:: ~ARM_OptimiseBase()
{
	delete itsInitialValuesVec;
	delete itsInitialValuesMat;
	delete itsFunction;
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_OptimiseBase& ARM_OptimiseBase::operator = (const ARM_OptimiseBase& rhs)
{
    if( this != &rhs )
	{
		ARM_Object::operator=( rhs );
		delete itsInitialValuesVec;
		delete itsInitialValuesMat;
		itsTolerance		= rhs.itsTolerance;
		itsStepMax			= rhs.itsStepMax;
		itsLocalSearch		= rhs.itsLocalSearch;
		itsInitialValuesVec = rhs.itsInitialValuesVec ? (ARM_GP_Vector*) rhs.itsInitialValuesVec->Clone() : NULL;
		itsInitialValuesMat = rhs.itsInitialValuesMat ? (ARM_GP_Matrix*) rhs.itsInitialValuesMat->Clone() : NULL;
		delete itsFunction;
		itsFunction			= rhs.itsFunction? rhs.itsFunction->Clone() : NULL;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_OptimiseBase::ARM_OptimiseBase( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    const size_t max_iter,
	bool getDetails,
	double tolerance,
	double stepMax,
	bool localSearch,
	size_t factorNb,
	size_t NbCalibIter,
    ARM_MktTargetType  targetType )
:   
    ARM_ModelFitter(
		model,
        portfolio,
        modelParam,
        linkedModelFitter,
        previousModelFitter,
        max_iter,
        getDetails,
        factorNb,
        NbCalibIter,
        targetType ), 
	
	itsTolerance(tolerance),
	itsStepMax(stepMax),
	itsLocalSearch(localSearch),
	itsInitialValuesVec(NULL),
	itsInitialValuesMat(NULL),
	itsFunction(NULL)
{
    Validate();
	Initialise();
};


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: Initialise
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_OptimiseBase::Initialise()
{
	size_t i;
	size_t calibparamsize	= GetCalibParams().size();
	for(i=0; i<calibparamsize ; ++i)
		if( GetCalibParam(i)->GetAdviseBreakPointTimes() )
			GetPricingModel()->AdviseBreakPointTimes( GetPortfolio(), GetCalibParam(i), GetFactorNb() );

    itsInitialValuesVec = new ARM_GP_Vector();

	size_t j,l,p;
	for(i=0; i<GetCalibParams().size(); ++i)
	{
		if( ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam(i)) )
		{
			for(j=0; j<curveCalibParam->size(); ++j)
				(*itsInitialValuesVec).push_back((curveCalibParam->GetCurve()->GetOrdinates())[j]);
		}
		else if( ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(GetCalibParam(i)) )
		{
			for(j=0; j<surfaceCalibParam->rows(); ++j)
				for(p=0; p<surfaceCalibParam->cols(); ++p)
					(*itsInitialValuesVec).push_back((surfaceCalibParam->GetSurface()->GetX3().Elt(j,p)));
		}
		else if( ARM_SurfaceListModelParam* surfaceListCalibParam = dynamic_cast<ARM_SurfaceListModelParam*>(GetCalibParam(i)))
		{
			for (l=0; l<surfaceListCalibParam->size(); ++l)
				for(j=0; j<surfaceListCalibParam->GetSurface(l)->GetX1().size(); ++j)
					for(p=0; p<surfaceListCalibParam->GetSurface(l)->GetX2().size(); ++p)
						(*itsInitialValuesVec).push_back(surfaceListCalibParam->GetSurface(l)->GetX3().Elt(j,p));
		}
		else
            throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "type not supported!" );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: Validate
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_OptimiseBase::Validate()
{
    if(GetPortfolio() == ARM_StdPortfolioPtr(NULL))
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"optimise calibration needs non null portfolio!" );

    if( GetCalibParams().empty() )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"optimise calibration needs at least one calib param!" );
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: MergeVariables
///	Returns: 
///	Action : Regroup all parameters to calibrate
////////////////////////////////////////////////////

ARM_VectorVector ARM_OptimiseBase::MergeVariables()
{
    size_t variablesize = 0;
    size_t calibparamsize = GetCalibParams().size();

    /// vector of vector to store variables to calibrate
    /// and its boundaries
    ARM_VectorVector variablesvector(3);

    size_t i,j,k,l;
    for(i=0; i<3 ; ++i)
		variablesvector[i] = new ARM_GP_Vector();

    for(i=0;i<calibparamsize; ++i)
    {
		if (ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam(i)))
        {
			for(j=0; j<curveCalibParam->size(); ++j)
			{
				variablesvector[0]->push_back( curveCalibParam->GetCurve()->GetOrdinate(j));
				variablesvector[1]->push_back((*(curveCalibParam->GetLowerBound()))[j]);
				variablesvector[2]->push_back((*(curveCalibParam->GetUpperBound()))[j]);
			}
		}
		else if (ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(GetCalibParam(i)))
		{
			for(j=0; j<surfaceCalibParam->rows(); ++j)
			{
				for(k=0; k<surfaceCalibParam->cols(); ++k )
				{
					variablesvector[0]->push_back( surfaceCalibParam->GetSurface()->GetX3().Elt(j,k));
					variablesvector[1]->push_back( surfaceCalibParam->GetLowerBound());
					variablesvector[2]->push_back( surfaceCalibParam->GetUpperBound());
				}
			}
		}
		else if( ARM_SurfaceListModelParam* surfaceListCalibParam = dynamic_cast<ARM_SurfaceListModelParam*>(GetCalibParam(i)))
		{
			ARM_SurfacePtrVector& surfacelist = surfaceListCalibParam->GetSurfaceList() ;
			size_t taille = surfaceListCalibParam->size();

			for (l=0; l<surfaceListCalibParam->size(); ++l)
			{
				for(j=0; j<surfacelist[l]->GetX1().size(); ++j)
				{
					for(k=0; k<surfacelist[l]->GetX2().size(); ++k)
					{
						variablesvector[0]->push_back( surfacelist[l]->GetX3().Elt(j,k));
						variablesvector[1]->push_back( surfaceListCalibParam->GetLowerBound());
						variablesvector[2]->push_back( surfaceListCalibParam->GetUpperBound());
					}
				}
			}
		}
		else throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "type not supported!" );
    }

    return variablesvector;
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: SplitVariables
///	Returns: 
///	Action : Split the parameters after optimisation
////////////////////////////////////////////////////
void ARM_OptimiseBase::SplitVariables( const ARM_GP_Vector& variables)
{
    for(size_t i=0, k=0; i<GetCalibParams().size() ; ++i)
    {
		if (ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam(i)))
		{
            ARM_GP_Vector value;
			for(size_t j=0; j < curveCalibParam->size(); ++j,++k) 
                value.push_back(variables[k]);
            
            curveCalibParam->SetAndUpDate(&value);
		}
		else if (ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(GetCalibParam(i)))
		{
            std::vector<double> value;
			for(size_t j=0; j < surfaceCalibParam->size(); ++j,++k) 
                value.push_back(variables[k]);

			surfaceCalibParam->SetAndUpDate(&value);
		}
		else if( ARM_SurfaceListModelParam* surfaceListCalibParam = dynamic_cast<ARM_SurfaceListModelParam*>(GetCalibParam(i)))
		{
			ARM_SurfacePtrVector& surfacelist = surfaceListCalibParam->GetSurfaceList() ;
            std::vector<double> value;
			for (size_t l=0; l<surfaceListCalibParam->size(); ++l)
                for(size_t j=0; j < surfacelist[l]->size(); ++j,++k) 
                    value.push_back(variables[k]);

			surfaceListCalibParam->SetAndUpDate(&value);
		}
		else
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "type not supported!" );

        GetPricingModel()->GetModelParams()->MergeModelParam(GetCalibParam(i),GetFactorNb());
    }
	GetPricingModel()->AdviseCurrentCalib(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_OptimiseBase::toString( const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    size_t i;
    os << indent << ModelFitterName() <<"\n";
    string tmpIndent = indent + "|\t ";
    os << tmpIndent << "Calib parameters  : ";
	for(i=0;i<GetCalibParams().size();++i)
	os << GetCalibParam(i)->GetTypeString() << " (" << GetCalibParam(i)->GetTypeString() << ")\n";
    os << tmpIndent << "Corresponding errors\n";
    os << tmpIndent << "\n";

    /// for unix compatibility
    /// we use printf
    /// more traditional C++ should use std::iostream manipulation
    char temp1[5];
    char temp2[30];
    SetUpError();
    for(i=0; i<GetError()->GetRowsNb(); ++i)
    {
        os << tmpIndent << "[";
        sprintf( temp1, "%2d", i );
        os << temp1 << "]=  ";
		sprintf( temp2, "%+6.6E (=%+4.4f %s)", (*GetError())(i,0), (*GetError())(i,1)*100.0, "%" );
        os << temp2 << "\n";
    }
	os << tmpIndent << "\n";
    os << indent << "======> End Of Infos <==================\n";
    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: StoreDetailsInFile
///	Returns: void
///	Action : store the details in the file
////////////////////////////////////////////////////
void ARM_OptimiseBase::StoreDetailsInFileIfRequired( Nag_E04_Opt& options, NagError& fail )
{
	/// should we get some details!
	if(GetDetails())
	{

		ARM_ModelFitter::StoreDetailsInFileIfRequired();
        /// copy the fileName for NAG
		const char* fileName = GetFileName()->c_str();
		strcpy(options.outfile, fileName );

        options.print_level = Nag_Soln_Iter_Full;
        options.list		= TRUE;
		fail.print			= TRUE;
		options.output_level= Nag_MPS_List;
		options.print_deriv	= Nag_D_Full;
		options.verify_grad = Nag_SimpleCheck;

	}
	else
	{
		fail.print					= FALSE;
        options.list				= FALSE;
        options.print_level			= Nag_NoPrint;
        options.output_level		= Nag_NoOutput;
		options.minor_print_level	= Nag_NoPrint;
		options.print_deriv			= Nag_D_NoPrint;

	}
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: SetFunctionNoClone
///	Returns: void
///	Action : set the function witout cloning it
////////////////////////////////////////////////////
void ARM_OptimiseBase::SetFunctionNoClone( MultiDimFunc* function )
{
	delete itsFunction;
	itsFunction = function;
}



////////////////////////////////////////////////////
///	Routine: ARM_OptimiseBase
///	Returns: ProcessNagWarning
///	Action : looks at nag warning and store it if flag set
////////////////////////////////////////////////////
void ARM_OptimiseBase::ProcessNagWarning( NagError fail )
{
	ProcessNagWarningOnObj( fail, *this );
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

