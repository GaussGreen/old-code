

#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libicm_local\icm_local_pf.h>
#include <ICMKernel\util\icm_matrix.h>

#include <ARMKernel\inst\security.h>
#include <ICMKernel\inst\icm_leg_basket.h>
#include <ICMKernel\inst\icm_collateral.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>




long ICMLOCAL_Portfolio (const VECTOR<long>& SecuritiesID,
						 const long& CashFlowId,
						 ARM_result& result,
						 long objId )
{
	long pfId = 0;
	int i = 0, j=0;

	ICM_Portfolio* prevpf = NULL;
	ICM_Portfolio* newpf = NULL;

	ARM_Security** pSecurities = NULL;
	ARM_Security* pSecurity = NULL;

	ICM_Parameters* params = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg(""); 

	int size = SecuritiesID.size();
	int nbassets = 0;

	try
	{
		//On filtre les bads inputs
		for (j = 0; j <size; j++)
		{
			if (SecuritiesID[j]>=0)
				nbassets++;
		}

		pSecurities = new ARM_Security*[nbassets];

		for (j = 0; j < nbassets; j++)
		{

			pSecurity = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(SecuritiesID[j]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pSecurity, ARM_SECURITY) == 0)
			{
				result.setMsg ("ARM_ERR: one Security is not of a good type");
			
				return ARM_KO;
			}

			pSecurities[j] = pSecurity;
		}	

		params = (ICM_Parameters*) LOCAL_PERSISTENT_OBJECTS->GetObject(CashFlowId);

		if (params)
			newpf = new ICM_Portfolio(pSecurities,nbassets,params);
		else
			newpf = new ICM_Portfolio(pSecurities,nbassets);

		if (newpf == NULL)
		{
			result.setMsg ("ARM_ERR:Portfolio is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			pfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Portfolio*)newpf);

			if (pfId == RET_KO)
			{
				if (newpf)
					delete newpf;
				newpf = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(pfId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* PfPf = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prevpf = (ICM_Portfolio*) (PfPf);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PfPf, ARM_SECURITY) == 1)
			{
				if (prevpf)
				{
					delete prevpf;
					prevpf = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Portfolio*)(newpf), objId);

				return ARM_OK;
			}
			else
			{
				if (newpf)
					delete newpf;
				newpf = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newpf)
			delete newpf;
		newpf = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_Leg_Basket (const VECTOR<long>& SecuritiesID,
						  const long& CashFlowId,
						  ARM_result& result,
						  long objId )
{
	long pfId = 0;
	int i = 0, j=0;

	CCString msg(""); 
	ICM_Leg_Basket* prevpf = NULL;
	ICM_Leg_Basket* newpf = NULL;

	ARM_Security** pSecurities = NULL;
	ARM_Security* pSecurity = NULL;

	ICM_Parameters* params = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	

	int size = SecuritiesID.size();
	int nbassets = 0;

	try
	{
		//On filtre les bads inputs
		for (j = 0; j <size; j++)
		{
			if (SecuritiesID[j]>=0)
				nbassets++;
		}

		pSecurities = new ARM_Security*[nbassets];

		for (j = 0; j < nbassets; j++)
		{

			pSecurity = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(SecuritiesID[j]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pSecurity, ARM_SECURITY) == 0)
			{
				result.setMsg ("ARM_ERR: one Security is not of a good type");
			
				return ARM_KO;
			}

			pSecurities[j] = pSecurity;
		}	

		params = (ICM_Parameters*) LOCAL_PERSISTENT_OBJECTS->GetObject(CashFlowId);

		newpf = new ICM_Leg_Basket(pSecurities,nbassets,params);

		if (newpf == NULL)
		{
			result.setMsg ("ARM_ERR:Portfolio is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			pfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Leg_Basket*)newpf);

			if (pfId == RET_KO)
			{
				if (newpf)
					delete newpf;
				newpf = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(pfId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* PfPf = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prevpf = (ICM_Leg_Basket*) (PfPf);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PfPf, ARM_SECURITY) == 1)
			{
				if (prevpf)
				{
					delete prevpf;
					prevpf = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Leg_Basket*)(newpf), objId);

				return ARM_OK;
			}
			else
			{
				if (newpf)
					delete newpf;
				newpf = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newpf)
			delete newpf;
		newpf = NULL;

		ARM_RESULT();
	}
}

//	-----------------------------------------------------------------------------------------------------------------------
long 
ICMLOCAL_Collateral (const vector<string>& Labels,
					 const ARM_Vector & Notionals,
					 // const vector<bool>& IsInDefault,
					 long prevId)
{
	long pfId = 0;
	int i = 0, j=0;

	ICM_Collateral* prevpf = NULL;
	ICM_Collateral* newpf = NULL;

	int size = Labels.size();


	newpf = new ICM_Collateral(Labels,
							   Notionals
							   // IsInDefault
							   );
	if (!newpf) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Collateral: can't create ICM_Collateral") ;

	return LocalPersistent::get().adopt(newpf,prevId) ; 
}
//	-----------------------------------------------------------------------------------------------------------------------
long 
ICMLOCAL_VariableCollateral (const vector<string>& Labels,
							 const ICM_Vector<long>& notionalIds,
							 long prevId)
{
	

	int size = Labels.size();
	ICM_Vector<ARM_ReferenceValue> notionals(notionalIds.size()); 
	for(int i=0;i<notionalIds.size();i++) 
	{
		ARM_ReferenceValue* item=0; 
		LocalPersistent::get().convert(notionalIds[i],item); 
		if (!item) ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_VariableCollateral: can't retrieve RefValue "<<i); 
		notionals[i]=*item; 
	}

	ICM_Collateral* newpf = new ICM_Collateral(); 
	if (!newpf) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Collateral: can't create ICM_Collateral") ;
	newpf->Set(Labels,notionals);

	return LocalPersistent::get().adopt(newpf,prevId) ; 
}