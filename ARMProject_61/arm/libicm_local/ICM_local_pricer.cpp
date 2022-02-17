
#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARMKernel\inst\security.h"
#include "ICMKernel\inst\icm_gen.h"

/// remove warnings on va_start and va_end
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
/// and do not change the order as this one has to be precisely here
/// if you do not know what you are doing, please ask someone who does!
#include <ARM\libarm_local\undef_va_vars.h>
#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ICMKernel/pricer/icm_pricer_tree_binomial_spreadoption.h"
#include "ICMKernel/pricer/icm_pricer_cds.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous.h"
#include "ICMKernel\mod\modelmulticurves.h"

long ICMLOCAL_Pricer (const ARM_Date*  asof_,
					  long idSecurity, 
					  long idModel,
					  int PricerType,
					  int nbpaths,
					  long idParameters,
					  // double AsOfdate,
					  // long idMktDataManager,
					  // ARM_result& result, 
					  long prevId)
{

	ARM_Security* security; 
	LocalPersistent::get().convert(idSecurity,security ); 

	ARM_Object * model; 
	LocalPersistent::get().convert(idModel,model); 
 
	ICM_Parameters* Parameters; 
	LocalPersistent::get().convert(idParameters,Parameters); 

	// ICM_MktDataMng* MktDataManager = NULL;
	// LocalPersistent::get().convert(idMktDataManager,MktDataManager); 

	ICM_Pricer_Advisor advisor = ICM_Pricer_Advisor();
	
	ARM_Date asof ;
	if (asof_) asof=*asof_ ;
	else 
	{
		ARM_Model* mod_ = dynamic_cast<ARM_Model*>(model); 
		if (!mod_)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Pricer: AsOf required (can't cast Model to ARM_Model)"); 
		asof = mod_->GetStartDate(); 
	}
	ICM_Pricer* Pricer = advisor.GeneratePricer(security, model, (ARM_CLASS_NAME)PricerType,nbpaths,Parameters,asof);

	if (!Pricer) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Pricer: Can't create pricer."); 

	return LocalPersistent::get().adopt(Pricer,prevId); 
}


/** 
double ICMLOCAL_Debug_Function (long pricerId, double Double, VECTOR<double>& Data, ARM_result& result)
{
	ICM_Pricer_Tree_Binomial_SpreadOption* pricer=NULL;
	double Aux = 0.;

	double* pDatas = NULL;

	pDatas = new double[Data.size()];

	for	(int j=0; j<Data.size(); j++)
		pDatas[j] = Data[j];
	

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{		
		pricer = (ICM_Pricer_Tree_Binomial_SpreadOption*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);


		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		Aux = pricer->GenericLegPrice((double) Double, (double*) pDatas);
		
		result.setDouble(Aux);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Debug_Function : unrecognized failure");
		return ARM_KO;
	}

}

**/ 
long ICMLOCAL_PricerDefaultCdsNew (const ARM_Date* asof_,
								   long idSecurity, 
								long idModel,
								// ARM_result& result, 
								long prevId)
{

 	ARM_Security* security; 
	LocalPersistent::get().convert(idSecurity,security); 
	ARM_Model* model; 
	LocalPersistent::get().convert(idModel,model); 

	ARM_Date asof ; 
	if (asof_) asof=*asof_; 
	else asof = model->GetStartDate(); 
	ICM_Pricer_Cds *Pricer = new ICM_Pricer_Cds ; 
	Pricer->Set(security,model,ICM_Parameters(),asof);

	if (Pricer == NULL)
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_PricerDefaultCdsNew: Can't create pricer."); 

	return LocalPersistent::get().adopt(Pricer,prevId); 
}


long ICMLOCAL_SetVolatility(long pricerId,
							long idvolcurve,
							ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	ARM_VolCurve* pVolCurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pVolCurve = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(idvolcurve);

		ARM_Date NewAsOfDate;

		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		ICM_ModelMultiCurves* Model = (ICM_ModelMultiCurves*) pricer->GetModel();
		Model->SetVolCurve(pVolCurve) ;
		// pricer->SetPriceFlg(false);
		pricer->unsetFlg(qCMPPRICE); 
		pricer->SetInitialPriceFlg(false);


		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetVolatility : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_GetPricer_DataFromLabel(
							long				PricerId,
							CCString			DataLabel,
							ARM_result&			result)
{
	double	dResult=0.;
	string	sResult;

	ICM_Pricer* ThePricer = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// ---------------------------------------------------------------------------------
		ThePricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(PricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ThePricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// discrimate between STRING or NUMERICS DATA from LABEL
		string	MyString(DataLabel);

		if (MyString.find("STR_") == string::npos)
		{
			// NUMERICS
			ThePricer->GetDataFromLabel(MyString, dResult);
			result.setDouble(dResult);
		}
		else
		{
			// STRING
			ThePricer->GetDataFromLabel(MyString, sResult);
			result.setMsg(sResult.c_str());
		}

		// ---------------------------------------------------------------------------------
		

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}

long ICMLOCAL_GenerateImpliedCurve (long pricerId, ARM_result& result,long objId)
{
	ICM_DefaultCurve* curve=NULL;
	ICM_DefaultCurve* prevcurve=NULL;
	ICM_Pricer_Distrib* pricer=NULL;
	double price=0.0;
	long curveId;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	CCString msg ("");

	try
	{
		pricer = (ICM_Pricer_Distrib *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		curve = pricer->GenerateImpliedCurve();

		if (curve == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)curve);

			if (curveId == RET_KO)
			{
				if (curve)
					delete curve;
				curve = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);
			return ARM_OK;
		}
		else
		{
			prevcurve = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(prevcurve, ICM_DEFAULTCURVE) == 1)
			{
				if (prevcurve)
				{
					delete prevcurve;
					prevcurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)curve, objId);

				return ARM_OK;
			}
			else
			{
				if (curve)
					delete curve;
				curve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Price : unrecognized failure");
		return ARM_KO;
	}

}
