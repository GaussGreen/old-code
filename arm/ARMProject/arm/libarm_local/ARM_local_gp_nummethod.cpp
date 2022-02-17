/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_nummethod.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_local_gp_nummethod.cpp
 *
 *  \brief file for the numerical method part of the generic pricer local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include "firstToBeIncluded.h"

///////////////////////////////////////////////
/// WARNING include these headers FIRST
/// because it uses the new style headers
/// CCxl uses cstyle headers and leads to conflict
/// if defined first
///////////////////////////////////////////////

#include <GP_Infra\gpinfra\typedef.h>

#include <GP_NumMethods\gpnummethods\typedef.h>
#include <GP_NumMethods\gpnummethods\binummethod.h>
#include <GP_NumMethods\gpnummethods\finummethod.h>
#include <GP_NumMethods\gpnummethods\mixtenummethod.h>
#include <GP_NumMethods\gpnummethods\treemethod.h>
#include <GP_NumMethods\gpnummethods\mcmethod.h>
#include <GP_NumMethods\gpnummethods\amc_andersen.h>
#include <GP_NumMethods\gpnummethods\amc_ls.h>
#include <GP_NumMethods\gpnummethods\treebase.h>
#include <GP_NumMethods\gpnummethods\treefactory.h>
#include <GP_NumMethods\gpnummethods\scheduler.h>
#include <GP_NumMethods\gpnummethods\schedulerfactory.h>
#include <GP_NumMethods\gpnummethods\sampler.h>
#include <GP_NumMethods\gpnummethods\samplerfactory.h>
#include <GP_NumMethods\gpnummethods\impsampler.h>
#include <GP_NumMethods\gpnummethods\impsamplerfactory.h>
#include <GP_NumMethods\gpnummethods\pathscheme.h>
#include <GP_NumMethods\gpnummethods\pathschemefactory.h>
#include <GP_NumMethods\gpnummethods\truncator.h>
#include <GP_NumMethods\gpnummethods\reconnector.h>
#include <GP_NumMethods\gpnummethods\smoother.h>
#include <GP_NumMethods\gpnummethods\amcmethod.h>
#include <GP_NumMethods\gpnummethods\amc_exercboundcalc.h>
#include <GP_NumMethods\gpnummethods\cfmethod.h>
#include <GP_NumMethods\gpnummethods\pdemethod.h>
#include <GP_NumMethods\gpnummethods\argconvdefault.h>
#include <GP_NumMethods\gpnummethods\impsampleropt.h>

#include <GP_NumLib\gpnumlib\randomgenfactory.h>
#include <GP_NumLib\gpnumlib\compositegen.h>
#include <GP_NumLib\gpnumlib\antitheticgen.h>
#include <GP_NumLib\gpnumlib\argconvdefault.h>

#include <GP_Base\gpbase\singleton.h>
#include <GP_Base\gpbase\gplinalgconvert.h>
#include <GP_Base\gpbase\gpvector.h>
#include <GP_Base\gpbase\gpmatrix.h>
#include <GP_Base\gpbase\typedef.h>
#include <GP_Base\gpbase\eventviewerfwd.h>
#include <GP_Base\gpbase\curve.h>

#include <gpbase\autocleaner.h>

#include <GP_Infra\gpinfra\numeraire.h>
#include <GP_Infra\gpinfra\numerairefactory.h>
#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\nummethod.h>
#include <GP_Infra\gpinfra\pricingstates.h>
#include <GP_Infra\gpinfra\gensecurity.h>
#include <GP_Infra\gpinfra\pricingadviser.h>
#include <GP_Infra\gpinfra\argconvdefault.h>

#include "refvalue.h"
#include "ARM_local_gp_nummethod.h"
#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"

#include "ARM_local_glob.h"
#include "ARM_local_wrapper.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "glob\expt.h"
#include "util\fromto.h"


//// using the namespace directive to access ARM object!
using ARM::ARM_IntVector;
using ARM::ARM_GP_T_Vector;
using ARM::std::vector<double>;
using ARM::ARM_BINumMethod;
using ARM::ARM_FINumMethod;
using ARM::ARM_MixteNumMethod;
using ARM::ARM_Numeraire;
using ARM::ARM_NumeraireFactory;
using ARM::ARM_TreeMethod;
using ARM::ARM_PricingModel;
using ARM::ARM_AutoCleaner;
using ARM::ARM_NumMethod;
using ARM::ARM_MCMethod;
using ARM::ARM_AMCMethod;
using ARM::ARM_AMCAndersen;
using ARM::ARM_AMCLongstaffSchwartz;
using ARM::ARM_ExerciseBoundaryCalc;
using ARM::ARM_RandomGenerator;
using ARM::ARM_RandomGeneratorPtr;
using ARM::ARM_RandomGeneratorPtrVector;
using ARM::ARM_RandGenFactory;
using ARM::ARM_RandGenFactoryImp;
using ARM::ARM_AntitheticOneGen;
using ARM::ARM_ArgConv_BaseGenAlgoType;
using ARM::ARM_ArgConv_TransformAlgoType;
using ARM::CreateARMGPVectorFromVECTOR;
using ARM::ARM_SchedulerBase;
using ARM::ARM_SchedulerFactory;
using ARM::ARM_SamplerBase;
using ARM::ARM_SamplerFactory;
using ARM::ARM_TreeBase;
using ARM::ARM_TreeFactory;
using ARM::ARM_TreeFactoryData;
using ARM::ARM_NumMethodPtr;
using ARM::ARM_NumerairePtr;
using ARM::ARM_CFMethod;
using ARM::ARM_ArgConv_CFmethodType;
using ARM::ARM_PDEMethod;
using ARM::ARM_PDENumericalScheme;
using ARM::ARM_ArgConv_PDENumSchemeType;
using ARM::ARM_ArgConv_PDEBoundConditionType;
using ARM::ARM_ArgConv_PDEGridType;
using ARM::ARM_TheEventViewer;
using ARM::ARM_MultiCurve;
using ARM::ARM_GP_MultiCurvePtr;
using ARM::ARM_LinInterpCstExtrapolVec;
using ARM::ARM_ImpSampler;
using ARM::ARM_ImpSamplerPtr;
using ARM::ARM_ImpSamplerFactory;
using ARM::ARM_PathScheme;
using ARM::ARM_PathSchemePtr;
using ARM::ARM_PathSchemeFactory;
using ARM::ARM_ArgConv_RandGenOrder;
using ARM::ARM_ImpSamplerOpt;
using ARM::ARM_GenSecurity;
using ARM::ARM_StepUpRightOpenCstExtrapolVec;
using ARM::ARM_Regression;
using ARM::ARM_ArgConv_RegMode;
using ARM::ARM_GP_Matrix;
using ARM::ARM_ArgConv_YesNo;


extern string ModelNumMethodToClass(long modelId)
{
    ARM_PricingModel* model = dynamic_cast<ARM_PricingModel *>(LOCAL_PERSISTENT_OBJECTS->GetObject(modelId));
	if (!model)
		return "";

    ARM_NumMethod* numMethod = &*model->GetNumMethod();

    if(dynamic_cast<ARM_BINumMethod *>(numMethod))
        return LOCAL_BINUMMETHOD_CLASS;

    else if(dynamic_cast<ARM_FINumMethod *>(numMethod))
        return LOCAL_FINUMMETHOD_CLASS;

	else if(dynamic_cast<ARM_MixteNumMethod *>(numMethod))
        return LOCAL_MIXTENUMMETHOD_CLASS;

    else if(dynamic_cast<ARM_MCMethod *>(numMethod))
        return LOCAL_MCMETHOD_CLASS;

    else if(dynamic_cast<ARM_TreeMethod *>(numMethod))
        return LOCAL_TREEMETHOD_CLASS;

    else if(dynamic_cast<ARM_TreeBase *>(numMethod) || dynamic_cast<ARM_TreeFactoryData *>(numMethod))
        return LOCAL_TREEND_CLASS;

    else
        return LOCAL_ANY_CLASS;
}


////////////////////////////////////////////
//// Function to create a backward induction num method
////////////////////////////////////////////
extern long ARMLOCAL_BINumMethod_Create(
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_BINumMethod* method= NULL;

	try
	{
		method = new ARM_BINumMethod;

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to create a backward induction num method
////////////////////////////////////////////
extern long ARMLOCAL_FINumMethod_Create(
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_FINumMethod* method= NULL;

	try
	{
		method = new ARM_FINumMethod;

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a mixte induction num method
////////////////////////////////////////////
extern long ARMLOCAL_MixteNumMethod_Create(
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_MixteNumMethod* method= NULL;

	try
	{
		method = new ARM_MixteNumMethod;

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}


}



////////////////////////////////////////////
//// Function to set a numerical method
////////////////////////////////////////////
extern long ARMLOCAL_SetNumMethod(
	long modelId,
	long numMethodId,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_PricingModel* mod= NULL;
    ARM_TreeBase* tree=NULL;

	try
	{

		ARM_PricingModel* oldMod= NULL;
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !GetObjectFromId( &oldMod, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};

		ARM_NumMethod* method=NULL;
		if( !GetObjectFromIdWithDynamicCastCheck( &method, numMethodId))
		{
            ARM_TreeFactoryData* treeFactoryData;
		    if( GetObjectFromIdWithDynamicCastCheck( &treeFactoryData, numMethodId) )
		    {
                /// Convert the TreeFactoryData to Tree1D, Tree2D or Tree3D depending of model number factors
		        mod = (ARM_PricingModel*) oldMod->Clone();
                tree = ARM_TreeFactory.Instance()->CreateTreeND(mod->FactorCount(),*treeFactoryData);
		        mod->SetNumMethod( ARM_NumMethodPtr( tree ) ); /// tree is cloned here
            }
            else
            {
			    result.setMsg ("ARM_ERR: numerical method is not of a good type");
			    return ARM_KO;
            }
		}
        else
        {
		    /// clone it to avoid side effect
		    mod = (ARM_PricingModel*) oldMod->Clone();
		    mod->SetNumMethod( ARM_NumMethodPtr( static_cast<ARM_NumMethod*>( method->Clone() ) ) );
        }


		/// beware that compared to standard assignObject macro, this one checks the
		/// root name and not the class name... so use it carefully
		/// and if you do not know what you are doing, please ask someone
		/// in fact this assign object checks that the object is really a derived object from pricing model!
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		delete tree;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete mod;
		delete tree;
		result.setMsg ("ARM_ERR: unrecognized failure in setting numerical method");
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to set numeraire to model
////////////////////////////////////////////
extern long ARMLOCAL_SetNumeraire(
	long modelId,
	long numeraireId,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_PricingModel* mod= NULL;

	try
	{

		ARM_PricingModel* oldMod= NULL;
		if( !GetObjectFromId( &oldMod, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};

		ARM_Numeraire* numeraire= NULL;
		if( !GetObjectFromId( &numeraire, numeraireId, ARM_NUMERAIRE ) )
		{
			result.setMsg ("ARM_ERR: numerical method is not of a good type");
			return ARM_KO;
		};

		mod = (ARM_PricingModel*) oldMod->Clone();
		mod->SetNumeraire( ARM_NumerairePtr( static_cast<ARM_Numeraire*>( numeraire->Clone() ) ) );

		/// beware that compared to standard assignObject macro, this one checks the
		/// root name and not the class name... so use it carefully
		/// and if you do not know what you are doing, please ask someone
		/// in fact this assign object checks that the object is really a derived object from pricing model!
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}





////////////////////////////////////////////
//// Function to create a Numeraire
////////////////////////////////////////////
extern long ARMLOCAL_Numeraire_Create(
    long                    numeraireType,
    const VECTOR<double>&   numeraireTimes,
	ARM_result&	            result, 
	long                    objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_Numeraire* numeraire= NULL;

	try
	{       
	    numeraire = ARM_NumeraireFactory.Instance()->CreateNumeraire( numeraireType );
		
		/// assign object
		if( !assignObject( numeraire, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete numeraire;
		x.DebugPrint();
		ARM_RESULT();
	}
}




////////////////////////////////////////////
//// Function to create a Tree method
////////////////////////////////////////////
extern long ARMLOCAL_TreeMethod_Create(
    int         nbSteps,
    double      stdDev,
    double      minStdDev,
    int         nbMinSteps,
    bool        isTree1GForced,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_TreeMethod* method= NULL;
    ARM_TreeFactoryData* treeFactoryData=NULL;

	try
	{
        if(isTree1GForced)
        {
		    method = new ARM_TreeMethod(stdDev,minStdDev,nbMinSteps);
            method->SetNbSteps(nbSteps);

		    /// assign object
		    if( !assignObject( method, result, objId ) ){
			    return ARM_KO; }
		    else{
			    return ARM_OK; }
        }
        else
        {
            std::vector<double> schedulerDatas(3);
            schedulerDatas[0]=nbSteps;
            schedulerDatas[1]=nbMinSteps;
            schedulerDatas[2]=minStdDev;
            treeFactoryData = new ARM_TreeFactoryData(
                ARM::ARM_SchedulerBase::ConstantVarianceMeanReverting,schedulerDatas,
                ARM::ARM_SamplerBase::MeanReverting, std::vector<double>(1,minStdDev),
                ARM::ARM_TruncatorBase::StandardDeviation,std::vector<double>(1,stdDev),false,
                ARM::ARM_ReconnectorBase::Mean,ARM::ARM_SmootherBase::DoNothing);

		    /// assign object
		    if( !assignObject( treeFactoryData, result, objId ) ){
			    return ARM_KO; }
		    else{
			    return ARM_OK; }
        }
	}
	
	catch(Exception& x)
	{
		delete method;
        delete treeFactoryData;
		x.DebugPrint();
		ARM_RESULT();
	}
}




////////////////////////////////////////////
//// Function to create a Monte Carlo method
////////////////////////////////////////////
extern long ARMLOCAL_MCMethod_Create(
    const size_t&			itersNb,
    const int&				fixStep,
	const VECTOR<long>&		randGensId,
	const long&				C_SamplerType,
	const VECTOR<double>&	C_SamplerDatas,
	const long&				C_SchedulerType,
    const VECTOR<double>&	C_SchedulerDatas,
	const long&				C_ImpSamplerType,
	const VECTOR<double>&	C_ImpSamplerData,
	const long&				C_PathSchemeType,
	const long&				ExercBoundCalcId,
	const size_t&			MaxBucketSize,
	ARM_result&				result, 
	long					objId )
{
	CC_Ostringstream os;

	os << "MaxBucketSize" << MaxBucketSize << endl;

	ARM_TheEventViewer.Instance()->AddToMessage(os.str());

	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_SchedulerBase* scheduler = NULL;
	ARM_SamplerBase* sampler = NULL;
	ARM_ImpSamplerPtr impSampler(NULL);
	ARM_PathSchemePtr pathScheme(NULL);
	ARM_MCMethod* mcMethod= NULL;
	ARM_MultiCurve* alpha;
	bool isAlpha = false;

	try
	{
		//ARM_RandomGeneratorPtr numRandGen;
		ARM_RandomGeneratorPtrVector randGenVector;
		if( randGensId.size()!=0)
		{
			size_t size = randGensId.size();
			VECTOR<ARM_RandomGenerator*> randGens;
			ARM_RandomGenerator* randGen = NULL;
			size_t i;
			for(i=0; i<size; ++i)
			{                
				if( !GetObjectFromId( &randGen, randGensId[i], ARM_RANDOMGEN) )
				{
					result.setMsg ("ARM_ERR: random generator is not of a good type");
					return ARM_KO;
				}
				randGens.push_back(randGen);
			}
			randGenVector.resize(size);
			for(i=0;i<size; ++i)
			{
					randGenVector[i] = ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) randGens[i]->Clone() );
			}
		}
		else
		{
			/// default is MRGK5 with InvNormCumFast + antithetic variates
			/// mrgk5
			ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::MRGK5,
				ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

			/// antithetic box muller
			ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::BoxMuller,
				pBaseRandomGen ) );

			/// antithetic variates!
			randGenVector.resize(1);
			randGenVector[0] = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::AntitheticOne,
				normRandGen ) );
		}

		// Scheduler
		int schedulerType = C_SchedulerType;
        std::vector<double> schedulerDatas;

		if (C_SchedulerDatas.size())
			schedulerDatas = std::vector<double>(C_SchedulerDatas);
		else
			schedulerDatas.push_back(fixStep);

		int nbSteps;

		scheduler = ARM_SchedulerFactory.Instance()->CreateScheduler(
			schedulerType,
			schedulerDatas,
			nbSteps);

		for (int i = 0; i < schedulerDatas.size(); ++i)
		{
			CC_Ostringstream os;
			os << "Data " << i << " : " << schedulerDatas[i] << "\n";
			ARM_TheEventViewer.Instance()->AddToMessage(os.str());
		}

		// Sampler
		// We need a multi dimension sampler
		int multiDim = 2;
		int samplerType = C_SamplerType;
		std::vector<double> samplerDatas(C_SamplerDatas);

		sampler = ARM_SamplerFactory.Instance()->CreateSampler(
			multiDim,
			samplerType,
			samplerDatas,
			scheduler);

		// Imp Sampler
		int impSamplerType;

		if (C_ImpSamplerData.size () == 1)
		{
			// By Default if there is no curve alpha equal to 0
			if(!(alpha = dynamic_cast<ARM_MultiCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_ImpSamplerData[0]))))
			{
				result.setMsg ("ARM_ERR: The Imp Sampler Data should contain the alpha multi curve!");
				return ARM_KO;
			}
		}
		else
		{
			std::vector<double> abscisses(1, 0.0);
			ARM_GP_T_Vector<std::vector<double> > ordinates(1, std::vector<double>(1,0.0));

			alpha = new ARM_MultiCurve(abscisses,ordinates,new ARM_LinInterpCstExtrapolVec);
			isAlpha = true;
		}

		impSamplerType = C_ImpSamplerType;

		impSampler = ARM_ImpSamplerPtr(ARM_ImpSamplerFactory.Instance()->CreateImpSampler(
			impSamplerType,
			alpha));

		// Path Scheme
		int pathSchemeType = C_PathSchemeType;
		pathScheme = ARM_PathSchemePtr(ARM_PathSchemeFactory.Instance()->CreatePathScheme(pathSchemeType));


		if( ExercBoundCalcId != ARM_NULL_OBJECT )
		{ // Our MCMethod is in fact an AMCMethod
			ARM_ExerciseBoundaryCalc * exerciseBoundaryCalc;
			if( (!GetObjectFromId( &exerciseBoundaryCalc, ExercBoundCalcId, ARM_AMCANDERSEN )) &&
				(!GetObjectFromId( &exerciseBoundaryCalc, ExercBoundCalcId, ARM_AMCLONGSTAFFSCHWARTZ )) )
			{
				result.setMsg ("ARM_ERR: Exercise Boundary Calculator is not of a good type");
				return ARM_KO;
			}

			mcMethod = new ARM_AMCMethod(
				itersNb,
				randGenVector,
				sampler,
				exerciseBoundaryCalc,
				MaxBucketSize,
				impSampler,
				pathScheme);
		}
		else
		{ // No ExerciseBoundaryCalculator -> Normal MCMethod
			mcMethod = new ARM_MCMethod(
				itersNb,
				randGenVector,
				sampler, 
				MaxBucketSize,
				impSampler,
				pathScheme);
		}

		delete scheduler;
		scheduler = NULL;
		delete sampler;
		sampler = NULL;
		if (isAlpha)
			delete alpha;

		/// assign object
		if( !assignObject( mcMethod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete scheduler;
		scheduler = NULL;
		delete sampler;
		sampler = NULL;
		delete mcMethod;
		if (isAlpha)
			delete alpha;
		mcMethod = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete scheduler;
		scheduler = NULL;
		delete sampler;
		sampler = NULL;
		delete mcMethod;
		mcMethod = NULL;
		if (isAlpha)
			delete alpha;
		result.setMsg ("ARM_ERR: unrecognized failure in creating Monte Carlo method");
		return ARM_KO;
	}
}

////////////////////////////////////////////
//// Function to create a PDE method
////////////////////////////////////////////

//Initial version
extern long ARMLOCAL_PDEMethod_Create(
	const string&	MethodName,				// Method Name
    const vector<double>& SchedulerData,	// Scheduler Data
	const int&		SchedulerDataNbRows,	// GridNbRows
	const int&		SchedulerDataNbCols,	// GridNbCols
	const int&		SpaceItersNb,			// Space Iter Nb 
	const string&	GridTypeName,			// Grid Type
	const vector<double>& GridData,			// Grid Data
	const int&		GridDataNbRows,			// GridNbRows
	const int&		GridDataNbCols,			// GridNbCols
	const int&		YGridItersNb,			// YGridIterNb
	const int&		ZGridItersNb,			// ZGridIterNb
	const double&	Theta1,					// Theta1
	const double&	Theta2,					// Theta2
	const double&	Theta3,					// Theta3
	const string&	BoundaryConditionName,	// Bound Cond
	const double&	lambda,					// Lambda
	ARM_result&		result,
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_PDEMethod* PDEMethod= NULL;

	try
	{
		/// Gets the right NumScheme from its name
		int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber(MethodName);
		int boundConditionType = ARM_ArgConv_PDEBoundConditionType.GetNumber(BoundaryConditionName);
		int GridType = ARM_ArgConv_PDEGridType.GetNumber(GridTypeName);
		ARM_GP_Matrix matGridData(GridDataNbRows,GridDataNbCols,GridData);
		ARM_GP_Matrix matSchedulerData(SchedulerDataNbRows,SchedulerDataNbCols,SchedulerData);
		ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType,SpaceItersNb,YGridItersNb,ZGridItersNb,Theta1,Theta2,Theta3,boundConditionType,lambda,GridType,matGridData,matSchedulerData);


		/// Builds PDE NumMethod With this numscheme
		double StdDevNb = matGridData.Elt(0,0);//Convention pour pouvoir reutiliser le CN1F
		double TimeItersNb = matSchedulerData.Elt(0,0);//Convention pour pouvoir reutiliser le CN1F
		PDEMethod = new ARM_PDEMethod(numScheme,TimeItersNb, SpaceItersNb, StdDevNb);

		/// assign object
		if( !assignObject( PDEMethod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete PDEMethod;
		PDEMethod = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete PDEMethod;
		PDEMethod = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating PDE method");
		return ARM_KO;
	}
}

//ND generalisation
extern long ARMLOCAL_PdeND_Create(
	string			C_MethodName,
    const long&                C_SchedulerType,
	long			C_SchulerData,
	long			C_SpaceData,
	long			C_SchemeData,
	string			C_BoundCondName,
	ARM_result&		result, 
	long			objId )
{
	return ARM_OK;//for the moment
}

////////////////////////////////////////////
//// Function to create a random nb generator
////////////////////////////////////////////
extern long ARMLOCAL_RandGen_Create(
	const string& genType,
	const string& algo,
	const long& baseGen1Id,
	const long& baseGen2Id,
	const int& seed,
	const int& dim,
	const int& factorDim,
	const int& nbOfPoints,
	const double& nbStdDevs,
	const int& firstNbTimes,
	const int& firstNbDims,
	const string& order,
	const int& firstSimulations,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_RandomGenerator* randGen= NULL;

	try
	{
		// The 2 bas generators for the algorithm

		ARM_RandomGeneratorPtr pBaseGen1, pBaseGen2;

		if( baseGen1Id != ARM_NULL_OBJECT )
		{
			ARM_RandomGenerator* baseGen1 = NULL;
			if( !GetObjectFromId( &baseGen1, baseGen1Id, ARM_RANDOMGEN ) )
			{
				result.setMsg ("ARM_ERR: base random generator 1 is not of a good type");
				return ARM_KO;
			};
			pBaseGen1 = ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) baseGen1->Clone() );
		}

		if( baseGen2Id != ARM_NULL_OBJECT )
		{
			ARM_RandomGenerator* baseGen2 = NULL;
			if( !GetObjectFromId( &baseGen2, baseGen2Id, ARM_RANDOMGEN ) )
			{
				result.setMsg ("ARM_ERR: base random generator 2 is not of a good type");
				return ARM_KO;
			};
			pBaseGen2 = ARM_RandomGeneratorPtr( (ARM_RandomGenerator*) baseGen2->Clone() );
		}

		ARM_TheEventViewer.Instance()->AddToMessage("Object creation start");

		/// should be an argument of the Monte Carlo
		randGen = ARM_RandGenFactory.Instance()->CreateRandGen(
			(ARM_RandGenFactoryImp::BaseGenType)	 ARM_ArgConv_BaseGenAlgoType.GetNumber(genType),
			(ARM_RandGenFactoryImp::TransformAlgo)   ARM_ArgConv_TransformAlgoType.GetNumber(algo),
			pBaseGen1,
			pBaseGen2,
			seed, dim, factorDim, ARM_GP_T_Vector<size_t>(1,nbOfPoints), nbStdDevs, firstNbDims, firstNbTimes,
			(ARM_RandGenFactoryImp::RandGenOrder) ARM_ArgConv_RandGenOrder.GetNumber(order), firstSimulations );

		ARM_TheEventViewer.Instance()->AddToMessage("Object creation end");

		/// assign object
		if( !assignObject( randGen, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete randGen;
		randGen = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete randGen;
		randGen = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating a random generator" );
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to create a random nb generator
////////////////////////////////////////////

extern long ARMLOCAL_SimpleRandGen_Create(
	const string& genType1,
	const string& genType2,
	const string& algo1,
	const string& algo2,
	const int& firstNbTimes,
	const int& firstNbDims,
	const string& isAntithetic,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_RandomGenerator* randGen= NULL;

	try
	{
		bool isAntitheticFlag;


		if (isAntithetic == "Y")
			isAntitheticFlag = true;
		else if (isAntithetic == "N")
			isAntitheticFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: IsAntithetic is Y or N.");
			return ARM_KO;
		}
			

		randGen = ARM_RandGenFactory.Instance()->CreateSimpleRandGen( 
			genType1,
			genType2,
			algo1,
			algo2,
			"PathOrder",
			firstNbDims,
			firstNbTimes,
			isAntitheticFlag);


		/// assign object
		if( !assignObject( randGen, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete randGen;
		randGen = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete randGen;
		randGen = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating a random generator" );
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to draw nbs from a random nb generator
////////////////////////////////////////////
extern long ARMLOCAL_RandGen_DrawVector(
	long randomGenId,
	int size,
	VECTOR<double>& data,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_RandomGenerator* randGen= NULL;

	try
	{
		ARM_RandomGenerator* randGen = NULL;
		if( !GetObjectFromId( &randGen, randomGenId, ARM_RANDOMGEN ) )
		{
			result.setMsg ("ARM_ERR: base random generator is not of a good type");
			return ARM_KO;
		};

		std::vector<double> result(size);
		randGen->draw(result);

		data = VECTOR<double>(size);
		std::copy(result.begin(), result.end(), data.begin() );

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete randGen;
		randGen = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete randGen;
		randGen = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating a random generator" );
		return ARM_KO;
	}
}


///////////////////////////////////////////////
//// Get the numerical method of a GP model
///////////////////////////////////////////////
extern long ARMLOCAL_GetNumMethodFromModel(
        long modelId,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_PricingModel* model;
    ARM_Object* object=NULL;

	try
	{
		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

        object =  static_cast< ARM_Object* >(model->GetNumMethod()->Clone());

        /// Assign the object in the ARM cache
		if( !assignObject( object, result, objId ) )
        {
            delete object;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Set the spot proba computation flag of a tree method
////////////////////////////////////////////
extern long ARMLOCAL_Local_Tree_SetProbaFlag(
	long treeId,
	bool isSpotProba,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_TreeMethod* treeMethod = NULL;
		ARM_TreeBase* treeBase = NULL;
		ARM_TreeFactoryData* treeFactory = NULL;
        ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(treeId);

		if( (treeMethod=dynamic_cast<ARM_TreeMethod*>( armObj )) )
		    treeMethod->SetProbilityComputation( isSpotProba );

		else if( (treeBase=dynamic_cast<ARM_TreeBase*>( armObj )) )
            treeBase->SetComputeSpotProbas(isSpotProba);

		else if( (treeFactory=dynamic_cast<ARM_TreeFactoryData*>( armObj )) )
            treeFactory->SetProbasFlag(isSpotProba);

        else
		{
			result.setMsg ("ARM_ERR: Tree is not of a good type");
			return ARM_KO;
		};

		string txt( "Spot Probas : " );
		txt += isSpotProba? "Forced" : "Free";
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an Andersen Method
////////////////////////////////////////////
extern long ARMLOCAL_AMCAndersen_Create(
	const double& ItersNb,
	const bool& sortedMaximization,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_AMCAndersen * AMCAndersen = NULL;

	try
	{

		AMCAndersen = new ARM_AMCAndersen( ( size_t )ItersNb, sortedMaximization );

		/// assign object
		if( !assignObject( AMCAndersen, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete AMCAndersen;
		AMCAndersen = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete AMCAndersen;
		AMCAndersen = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating Monte Carlo method");
		return ARM_KO;
	}
}

////////////////////////////////////////////
//// Function to create an LongstaffSchwartz Method
////////////////////////////////////////////
extern long ARMLOCAL_AMCLongstaffSchwartz_Create(
	const double& ItersNb,
	const string& regMode,
	const double& span,
	const string& isAutomatic,
	const int& degree,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_AMCLongstaffSchwartz * AMCLongstaffSchwartz = NULL;

	try
	{
		int isAutomaticInt = ARM_ArgConv_YesNo.GetNumber(isAutomatic);

		AMCLongstaffSchwartz = new ARM_AMCLongstaffSchwartz( ( size_t )ItersNb, (ARM_Regression::RegressionMode) ARM_ArgConv_RegMode.GetNumber(regMode), span, isAutomaticInt, degree );

		/// assign object
		if( !assignObject( AMCLongstaffSchwartz, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete AMCLongstaffSchwartz;
		AMCLongstaffSchwartz = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete AMCLongstaffSchwartz;
		AMCLongstaffSchwartz = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating Monte Carlo method");
		return ARM_KO;
	}
}

////////////////////////////////////////////////////////
//// Function to create a Tree ND method
///////////////////////////////////////////////////////
extern long ARMLOCAL_TreeND_Create(
         const long&                C_NbDims,
         const long&                C_SchedulerType,
         const VECTOR<double>&      C_SchedulerDatas,
         const long&                C_SamplerType,
         const VECTOR<double>&      C_SamplerDatas,
         const long&                C_TruncatorType,
         const VECTOR<double>&      C_TruncatorDatas,
         const bool&                C_ProbasFlag,
         const long&                C_ReconnectorType,
         const long&                C_SmootherType,
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_TreeBase* tree                  = NULL;
	try
	{
        int nbDims = C_NbDims;
        int schedulerType = C_SchedulerType;
        std::vector<double> schedulerDatas(C_SchedulerDatas);
        int samplerType = C_SamplerType;
        std::vector<double> samplerDatas(C_SamplerDatas);
        int truncatorType = C_TruncatorType;
        std::vector<double> truncatorDatas(C_TruncatorDatas);
        bool probasFlag = C_ProbasFlag;
        int reconnectorType = C_ReconnectorType;
        int smootherType = C_SmootherType;

        tree = ARM_TreeFactory.Instance()->CreateTreeND(nbDims,
            schedulerType,schedulerDatas,samplerType,samplerDatas,truncatorType,truncatorDatas,probasFlag,reconnectorType,smootherType);

		/// assign object
		if( !assignObject( tree, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete tree;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete tree;
		result.setMsg ("ARM_ERR: unrecognized failure while creating a tree ND method");
		return ARM_KO;
	}
}

////////////////////////////////////////////
//// Function to create a CF method
////////////////////////////////////////////
extern long ARMLOCAL_CFMethod_Create(
		 const string&	MethodeName,
		 const long&	matrixId,	
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CFMethod* cfMethod= NULL;
	ARM_GP_Matrix * matrixParams= new ARM_GP_Matrix(0);
    ARM_AutoCleaner< ARM_GP_Matrix > HoldMP(matrixParams);

	try
	{
		//fundLevrage
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &matrixParams, matrixId,"Matrix of integral params",result ) ) return ARM_KO;

		ARM_CFMethod::ARM_CFMethodType methodType = (ARM_CFMethod::ARM_CFMethodType) ARM_ArgConv_CFmethodType.GetNumber(MethodeName);

		cfMethod = new ARM_CFMethod( methodType , *matrixParams );

		/// assign object
		if( !assignObject( cfMethod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete cfMethod;
		cfMethod = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete cfMethod;
		cfMethod = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating Closed form method");
		return ARM_KO;
	}
}

////////////////////////////////////////////////////////
//// Function to create a Tree ND method
///////////////////////////////////////////////////////
extern long ARMLOCAL_ImpSampler_Optimize(
         const long&                C_GenSecId,
         const long&                C_ModelId,
         const long&                C_InitGuessId,
		 const double&				C_InitGuess,
         const long&                C_LowerBoundId,
		 const double&              C_LowerBound,
		 const long&                C_UpperBoundId,
		 const double&              C_UpperBound,
		 const string&              C_WithMC,
		 const long&				C_NbSteps,
		 const string&              C_Bootstrap,
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_PricingModel* newMod                  = NULL;
	ARM_GenSecurity* sec                  = NULL;
	ARM_PricingModel* mod                  = NULL;

	ARM_MultiCurve* initGuess;
	ARM_MultiCurve* lowerBound;
	ARM_MultiCurve* upperBound;

	bool isInitGuess=false;
	bool isLowerBound=false;
	bool isUpperBound=false;

	try
	{
		if(!(sec = dynamic_cast<ARM_GenSecurity*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_GenSecId))))
		{
			result.setMsg ("ARM_ERR: The gen sec id should contains a generic security!");
			return ARM_KO;
		}

		if(!(mod = dynamic_cast<ARM_PricingModel*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_ModelId))))
		{
			result.setMsg ("ARM_ERR: The model id should contains a pricing model!");
			return ARM_KO;
		}

		if (C_InitGuessId != ARM_NULL_OBJECT)
		{
			if(!(initGuess = dynamic_cast<ARM_MultiCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_InitGuessId))))
			{
				result.setMsg ("ARM_ERR: The init guess should contains a multi curve!");
				return ARM_KO;
			}
		}
		else
		{
			std::vector<double> abscisses(1, 0.0);
			ARM_GP_T_Vector<std::vector<double> > ordinates(1, std::vector<double>(1,C_InitGuess));

			initGuess = new ARM_MultiCurve(abscisses,ordinates,new ARM_StepUpRightOpenCstExtrapolVec);
			isInitGuess = true;
		}
		
		if (C_LowerBoundId != ARM_NULL_OBJECT)
		{
			if(!(lowerBound = dynamic_cast<ARM_MultiCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_LowerBoundId))))
			{
				result.setMsg ("ARM_ERR: The lower bound should contains a multi curve!");
				return ARM_KO;
			}
		}
		else
		{
			std::vector<double> abscisses(1, 0.0);
			ARM_GP_T_Vector<std::vector<double> > ordinates(1, std::vector<double>(1,C_LowerBound));

			lowerBound = new ARM_MultiCurve(abscisses,ordinates,new ARM_StepUpRightOpenCstExtrapolVec);
			isLowerBound = true;
		}

		if (C_UpperBoundId != ARM_NULL_OBJECT)
		{
			if(!(upperBound = dynamic_cast<ARM_MultiCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_UpperBoundId))))
			{
				result.setMsg ("ARM_ERR: The upper bound should contains a multi curve!");
				return ARM_KO;
			}
		}
		else
		{
			std::vector<double> abscisses(1, 0.0);
			ARM_GP_T_Vector<std::vector<double> > ordinates(1, std::vector<double>(1,C_UpperBound));

			upperBound = new ARM_MultiCurve(abscisses,ordinates,new ARM_StepUpRightOpenCstExtrapolVec);
			isUpperBound = true;
		}

		bool withMCFlag;
		if (C_WithMC == "Y")
			withMCFlag = true;
		else if (C_WithMC == "N")
			withMCFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: With MC should contain 'Y' or 'N'!");
			return ARM_KO;
		}

		bool bootstrapFlag;
		if (C_Bootstrap == "Y")
			bootstrapFlag = true;
		else if (C_Bootstrap == "N")
			bootstrapFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: With MC should contain 'Y' or 'N'!");
			return ARM_KO;
		}

		newMod = ARM_ImpSamplerOpt::ComputeModel(
					sec,
					mod,
					*initGuess,
					*lowerBound,
					*upperBound,
					withMCFlag,
					C_NbSteps,
					bootstrapFlag);

		/// assign object
		if( !assignObject( newMod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		if (isInitGuess)
			delete initGuess;
		if (isLowerBound)
			delete lowerBound;
		if (isUpperBound)
			delete upperBound;
	}
	
	catch(Exception& x)
	{
		delete newMod;

		if (isInitGuess)
			delete initGuess;
		if (isLowerBound)
			delete lowerBound;
		if (isUpperBound)
			delete upperBound;

		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete newMod;

		if (isInitGuess)
			delete initGuess;
		if (isLowerBound)
			delete lowerBound;
		if (isUpperBound)
			delete upperBound;

		result.setMsg ("ARM_ERR: unrecognized failure in the optimization of the importance sampling.");
		return ARM_KO;
	}
}