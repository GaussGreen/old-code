/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file modelfitter.cpp
 *
 *  \brief  base class for all model fitter
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/modelfitter.h"

/// gpbase
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/surface.h"

///gpmodel
#include "gpmodels/MultiAssets.h"



/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"

/// gpcalib
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"

/// kernel headers
//#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_WarningKeeper
///	Routine: AddModelFitterWarning
///	Returns: 
///	Action : Adds Model Fitter Warning
////////////////////////////////////////////////////
void ARM_ModelFitter::AddModelFitterWarning( const ARM_ModelFitterPtr& modelFitter )
{
	if( ARM_WarningKeeper::KeepWarning )
	{
		if( modelFitter->GetWarning() )
		{
			if( GetWarning() )
				GetWarning()->AddToMessage( modelFitter->GetWarning()->toString() );
			else 
				SetWarning( static_cast<ARM_Warning*>( modelFitter->GetWarning()->Clone()) );
		}
	}
}




////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_ModelFitter::ARM_ModelFitter()
:	
	ARM_WarningKeeper(),
	itsLinkedModelFitter(NULL ),
    itsPreviousModelFitter(NULL ),
    itsMax_iter (ARM_MAX_ITER),
	itsGetDetails (false),
	itsFactorNb (0), 
	itsNbCalibIter(1),
	itsTargetType(ARM_CalibrationTarget::PriceTarget)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_ModelFitter::ARM_ModelFitter(const ARM_ModelFitter& rhs)
: ARM_RootObject(rhs), ARM_Timer( rhs ), ARM_WarningKeeper(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_ModelFitter::ARM_ModelFitter( ARM_PricingModel* model,
     const ARM_StdPortfolioPtr portfolio,
     const ARM_ModelParamVector& calibParams,
     ARM_ModelFitterPtr& linkedModelFitter,
     ARM_ModelFitterPtr& previousModelFitter,
     size_t max_iter,
	 bool getDetailsFlag,
	 size_t factorNb,
	 size_t NbCalibIter,
     ARM_MktTargetType targetType)
:   
	ARM_RootObject(),
	ARM_Timer(),
	ARM_WarningKeeper(),
    itsPricingModel(model),
    itsPortfolio(portfolio),
    itsCalibDirection(CalibDirection_Forward),
    itsLinkedModelFitter(linkedModelFitter),
    itsPreviousModelFitter(previousModelFitter),
    itsArgsVector(),
    itsGetDetails(getDetailsFlag),
	itsFileName(NULL),
    itsMax_iter(max_iter),
	itsCalibParamsType(),
    itsError( ARM_GP_MatrixPtr(NULL)),
	itsFactorNb(factorNb),
	itsNbCalibIter(NbCalibIter),
    itsTargetType(targetType)
{
   	DuplicateCloneablePointorVectorInPlace<ARM_ModelParam>(calibParams,itsCalibParam);
    Validate();
	int size = itsCalibParam.size();
	for(int i=0; i<size; ++i)
		itsCalibParamsType.push_back(itsCalibParam[i]->GetType());
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: Validate
///	Returns: 
///	Action : to validate calibParms Boundaries
////////////////////////////////////////////////////
void ARM_ModelFitter::Validate()
{	
	/// conversion from model param to calib param type
	size_t i;
    for(i=0; i<itsCalibParam.size(); ++i)
    {
		/// nothing to convert
		if(		dynamic_cast<ARM_CurveModelParam*>(itsCalibParam[i]) 
			||	dynamic_cast<ARM_SurfaceModelParam*>(itsCalibParam[i]) 
			||  dynamic_cast<ARM_SurfaceListModelParam*>(itsCalibParam[i]) ) 
		{
			break;
		}
	
		else 
		{
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unsupported type! supported is curve, surface and surfacelist calib model param" );
        }
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: DoCalibrationProcess
///	Returns: 
///	Action : does the loop of calibration
////////////////////////////////////////////////////
void ARM_ModelFitter::DoCalibrationProcess()
{
	///// start timing
	//ClockStartTime();
	//int iter = itsNbCalibIter;

	///// calibration of previous methods
	//if(itsPreviousModelFitter != ARM_ModelFitterPtr(NULL ))
	//{
	//	itsPreviousModelFitter->DoCalibrationProcess();
	//	AddModelFitterWarning( itsPreviousModelFitter ); 
	//}
	//
	//while (iter>0)
	//{
	//	if(GetWarning())
	//		EraseWarningMessage();
	//	/// calibration of the current method
	//	Calibrate();
	//	
	//	/// once the calibration is done, relaunch the linked model fitter
	//	if(itsLinkedModelFitter != ARM_ModelFitterPtr(NULL))
	//	{
	//		itsLinkedModelFitter->DoCalibrationProcess();
	//		/// store in the last loop the warning
	//		itsLinkedModelFitter->SetKeepCurrentWarning( true );
	//		AddModelFitterWarning( itsLinkedModelFitter ); 
	//	}
	//	
	//	iter--;
	//	
	//	if(iter>0)
	//	{
	//		//Postprocessing()?
	//		std::auto_ptr<ARM_ModelParam> tmp(itsCalibParam[0]);
	//		SetCalibParam(static_cast<ARM_ModelParam*> (GetPricingModel()->GetModelParams()->GetModelParam(GetCalibParam()->GetType()).Clone()),0);
	//		this->Initialise();
	//		//Preprocessing()?
	//	}
	//}
	///// gets the end time
	//ClockEndTime();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: PostProcessing
///	Returns: 
///	Action : to post-processes the modelpricing
////////////////////////////////////////////////////

void ARM_ModelFitter::PostProcessing() 
{
    if(itsLinkedModelFitter != ARM_ModelFitterPtr(NULL ))
        itsLinkedModelFitter->PostProcessing();
    if(itsPreviousModelFitter != ARM_ModelFitterPtr(NULL ))
        itsPreviousModelFitter->PostProcessing();
  
    itsPricingModel->PostProcessing(*this);
};

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_ModelFitter:: ~ARM_ModelFitter()
{
	CleanUp();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: Assignment operator
///	Returns: 
///	Action : operator =
////////////////////////////////////////////////////

ARM_ModelFitter& ARM_ModelFitter::operator = (const ARM_ModelFitter& rhs)
{
    if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );
		ARM_Timer::operator=(rhs);
		ARM_WarningKeeper::operator=(rhs);

		/// Delete the old ones.
		CleanUp();
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : copies the additional member variables defined
///	in this object
////////////////////////////////////////////////////

void ARM_ModelFitter::CopyNoCleanUp( const ARM_ModelFitter& rhs )
{
    itsPricingModel			= rhs.itsPricingModel;
    itsMax_iter				= rhs.itsMax_iter;
    itsPortfolio			= rhs.itsPortfolio;
    itsError				= rhs.itsError ;
    itsGetDetails			= rhs.itsGetDetails;
    itsCalibDirection		= rhs.itsCalibDirection;
	itsFileName				= rhs.itsFileName? new string(*rhs.itsFileName): NULL;
	itsNbCalibIter			= rhs.itsNbCalibIter;
	itsFactorNb				= rhs.itsFactorNb;
	itsCalibParamsType		= rhs.itsCalibParamsType;
    itsTargetType			= rhs.itsTargetType;
    itsLinkedModelFitter	= ARM_ModelFitterPtr( rhs.itsLinkedModelFitter == ARM_ModelFitterPtr(NULL) ? NULL: (ARM_ModelFitter*)rhs.itsLinkedModelFitter->Clone());
    itsPreviousModelFitter	= ARM_ModelFitterPtr( rhs.itsLinkedModelFitter == ARM_ModelFitterPtr(NULL) ? NULL: (ARM_ModelFitter*)rhs.itsPreviousModelFitter->Clone());
	itsFactorNb				= rhs.itsFactorNb;
    itsTargetType			= rhs.itsTargetType;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: CleanUp
///	Returns: 
///	Action : deletes any dynamic object defined at the object
///				level
////////////////////////////////////////////////////
void ARM_ModelFitter::CleanUp( )
{
	//// delete the file
	if (itsFileName)
		remove( itsFileName->c_str() );
    delete itsFileName;
	itsFileName			= NULL;

    DeletePointorVector<ARM_ModelParam>(itsCalibParam);
    DeletePointorVector<ARM_VanillaArg>(itsArgsVector);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: SetVanillaArgVector method
///	Returns: args vector accessor
///	Action : Set the args vector
////////////////////////////////////////////////////
void ARM_ModelFitter::SetVanillaArgVector(ARM_VanillaArgVector argsVector)
{
	size_t i;
	for (i = 0; i < itsArgsVector.size(); ++i)
		delete itsArgsVector[i];
	itsArgsVector = argsVector;
}

///////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: SetUpError
///	Returns: void
///	Action : Update errors for all portfolio afer calibration
////////////////////////////////////////////////////
void ARM_ModelFitter::SetUpError() const
{
	/*if(!itsPricingModel)
		return;
	ARM_MultiAssetsModel* model = dynamic_cast<ARM_MultiAssetsModel*>(itsPricingModel);
	double asOfDate = model ? (*model->GetModelMap())[itsFactorNb]->Model()->GetAsOfDate().GetJulian() : itsPricingModel->GetAsOfDate().GetJulian(); 
    ARM_VanillaArgPtr VanillaArg;
    
    int size = GetPortfolio()->GetSize();
    int i;
    itsError = ARM_GP_MatrixPtr(new ARM_GP_Matrix(size,2));

    for(i=0; i<size; ++i)
    {      
        VanillaArg = ARM_VanillaArgPtr(ARM_ConverterFromKernel::ConvertSecuritytoArgObject( GetPortfolio()->GetAsset(i),asOfDate ));
		VanillaArg->SetCurveName( itsPricingModel->GetModelNamePerFactor( itsFactorNb ) );
		(*itsError)(i,0) = VanillaArg->Price( itsPricingModel ) - (*GetPortfolio()->GetMktPrices())[i];
		if( (*GetPortfolio()->GetMktPrices())[i]!= 0)
			(*itsError)(i,1) = (*itsError)(i,0) / (*GetPortfolio()->GetMktPrices())[i];
		else
			(*itsError)(i,1) = 0.0;
    }*/
}


///////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: SetUpError
///	Returns: void
///	Action : Update errors for all portfolio afer calibration
////////////////////////////////////////////////////
ARM_ModelParamVector::const_iterator ARM_ModelFitter::FindCalibParamWType( int type ) const
{
	return CC_NS( std, find_if) ( itsCalibParam.begin(), itsCalibParam.end(), FindModelParamWEnumUnaryVersion( ARM_ModelParamType::ParamNb ( type ) ) );
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: detailedString
///	Returns: 
///	Action : when detail mode is on, helps to give more
///				details
////////////////////////////////////////////////////
string ARM_ModelFitter::detailedString( const string& indent ) const
{
	CC_Ostringstream os;
	/*if(GetDetails()&&GetFileName())
    {
		CC_Ostringstream os;
		/// get file name and open the corresponding file
		CC_NS(std,ifstream) myFile(GetFileName()->c_str());
		if(!myFile)
		{
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,  
				"Could not open file " + *GetFileName() );
		}

		os << indent << "=============================================\n";
		os << indent << "========>   Model Fitter's Details    <======\n";
		os << indent << "=============================================\n";
		
		/// read 1024 characters or up to the newline
		const int lineSize = 1024;
		char buffer[lineSize];
		char firstChar;
		while(myFile.get(firstChar) )
		{
			myFile.getline(buffer,lineSize);
			os << indent<< indent<<firstChar << buffer << "\n";
		}
		os << indent <<" " << "\n";
	
		for(int i = 0; i < itsCalibParam.size(); ++i)
        os << itsCalibParam[i]->toString(indent).c_str() ;		
		return os.str();
	}
    else
        os << indent <<"  NO DETAILS!  \n"; */
	return os.str();
}

///	Routine: StoreDetailsInFile
///	Returns: void
///	Action : store the details in the file
////////////////////////////////////////////////////
void ARM_ModelFitter::StoreDetailsInFileIfRequired()
{
	/// should we get some details!
	if(GetDetails())
	{
		/// if there is no file name
		/// creates one
		if(!itsFileName)
		{
			/// filename include time to make it unique
			CC_Ostringstream os;

            char fOutName[200];

            ARM_GetTmpAbsFile("MethodDetails_", fOutName);

			os << fOutName;

			size_t i=0;
			for( i=0; i<GetCalibParams().size()-1; ++i)
				os << GetCalibParam(i)->GetTypeString() + string("_");
			os << GetCalibParam(i)->GetTypeString();
			
			/// to avoid overwritting the same file add the number of previous method!
			i=0;
			if( GetPreviousModelFitter() != ARM_ModelFitterPtr(NULL) )
			{
				ARM_ModelFitterPtr previousModelFitter = GetPreviousModelFitter();
				++i;
				while( (previousModelFitter = previousModelFitter->GetPreviousModelFitter()) != ARM_ModelFitterPtr(NULL) )
					++i;
			}
			os << "_" << i << ".txt";
			itsFileName = new string(os.str());

			/// take the fileName and erase the data
			/// open a file for output
			/// and clear it using the std ofstream
			//CC_NS(std,ofstream) myFile(itsFileName->c_str());
			//myFile.close();
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitter
///	Routine: SetCalibParamToInitialValue
///	Returns: void
///	Action : initialise the calib param to their initial value
////////////////////////////////////////////////////
void ARM_ModelFitter::SetCalibParamToInitialValue()
{
	for( size_t i=0; i<GetCalibParams().size(); ++i )
	{
		ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam(i));
		if( curveCalibParam )
			curveCalibParam->SetToInitialCurve();
		/// FIX: should be extended for surfaceCalibParam
	}
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

