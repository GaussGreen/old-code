
#pragma warning(disable :4786 4541 4800 4250)

#include "firstToBeIncluded.h"
#include "CCdate.h"
#include "CCstring.h"

#include <math.h>

#include <ICMKernel\inst\icm_convertible.h>
#include <ICMKernel\inst\icm_cbasw.h>

#include "ARMKernel\inst\irindex.h"

#include <ARM\libarm_local\ARM_local_ccy.h>
#include <ARM\libarm_local\ARM_local_swap.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

extern long ICMLOCAL_CBOPTION(double		FirstDateIn,
							  double		LastDateIn,
							  int			OptionTypeIn,
							  double		StrikeIn,
							  int			StrikeTypeIn,
							  int			AccruedOnExIn,
							  int			BarrierTypeIn,
							  double		BarrierStrikeIn,
							  int			BarrierStrikeType,
							  ARM_result&	result,
							  long			objId)
{
	long cboptionId;
	ICM_CBOption* cboption	=NULL;
	ICM_CBOption* cbO_old	=NULL;

	char*pFirstDate=new char[11];
	char*pLastDate=new char[11];

	Local_XLDATE2ARMDATE(FirstDateIn,pFirstDate);
	Local_XLDATE2ARMDATE(LastDateIn,pLastDate);

	//Creation de l'objet ICM_CBOption************************
	
	CCString msg("");

	try
	{
		cboption = new ICM_CBOption((ARM_Date)				pFirstDate,
									(ARM_Date)				pLastDate,
									(qCBOPTIONTYPE)			OptionTypeIn,
															StrikeIn,
									(qCBOPTIONSTRIKETYPE)	StrikeTypeIn,
									(qCBOPTIONACCRUED)		AccruedOnExIn,
									(qCBOPTIONBARRIERTYPE)	BarrierTypeIn,
															BarrierStrikeIn,
									(qCBOPTIONSTRIKETYPE)   BarrierStrikeType);
		if(pFirstDate)
			delete[]pFirstDate;
		if(pLastDate)
			delete[]pLastDate;

		if(cboption==NULL)
		{
			result.setMsg("ARM_ERR: CB Option is null");
			return ARM_KO;
		};

		if(objId==-1)
		{
			CREATE_GLOBAL_OBJECT();

			cboptionId=LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cboption);

			if(cboptionId==RET_KO)
			{
				if(cboption)
					delete cboption;
				cboption=NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			};

			result.setLong(cboptionId);

			return ARM_OK;
		}
		else
		{
			cbO_old=(ICM_CBOption*)LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cbO_old,ICM_CBOPTION)==1)
			{
				if(cbO_old)
				{
					delete cbO_old;
					cbO_old=NULL;
				};

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cboption,objId);

				return ARM_OK;
			}
			else
			{
				if(cboption)
					delete cboption;
				cboption=NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			};
		};
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if(cboption)
			delete cboption;
		cboption=NULL;

		if(pFirstDate)
			delete[]pFirstDate;
		pFirstDate=NULL;
		
		if(pLastDate)
			delete[]pLastDate;
		pLastDate=NULL;

		ARM_RESULT();
	};

};


extern long ICMLOCAL_CONVERTIBLE(long		DefBondId,
								 long		StockId,
								 double		YieldIn,
								 const VECTOR<long>& CBOptionIds,
								 ARM_result&result,
								 long		objId)
{
	long convertId;
	int i=0;
	ICM_CBOption*pCBOption=NULL;
	vector<ICM_CBOption*> CBOptionsIn;

	ICM_Convertible* Convertible=NULL;
	ICM_Convertible* newConvertible=NULL;

	ICM_Bond* pBondId=NULL;
	ICM_Stock*pStockId=NULL;

	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS)==0)
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	};

	CCString msg("");

	try
	{
		//On check le bond sous-jacent***************************

		pBondId=(ICM_Bond*)LOCAL_PERSISTENT_OBJECTS->GetObject(DefBondId);
		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pBondId,ICM_BOND)==0)
		{
			result.setMsg("ARM_ERR: Defbond is not of a good type");
			return ARM_KO;
		};

		pStockId=(ICM_Stock*)LOCAL_PERSISTENT_OBJECTS->GetObject(StockId);
		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pStockId,ICM_STOCK)==0)
		{
			result.setMsg("ARM_ERR: Stock is not of a good type");
			return ARM_KO;
		};

		int n=CBOptionIds.size();
		//On va récupérer une à une les options*********
		for(i=0;i<n;i++)
		{
			pCBOption=(ICM_CBOption*)LOCAL_PERSISTENT_OBJECTS->GetObject(CBOptionIds[i]);
			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pCBOption,ICM_CBOPTION)==0)
			{
				result.setMsg("ARM_ERR: One CB Option is not of a good type");
				return ARM_KO;
			};
			CBOptionsIn.push_back(pCBOption);
		};

		//On crée l'OC à partir du DefBondId et des options*********************
		newConvertible=new ICM_Convertible(pBondId,pStockId,YieldIn,CBOptionsIn);

		if(newConvertible==NULL)
		{
			result.setMsg("Convertible is null");
			return ARM_KO;
		};

		if(objId==-1)
		{
			CREATE_GLOBAL_OBJECT();

			convertId=LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newConvertible);

			if(convertId==RET_KO)
			{
				if(newConvertible)
					delete newConvertible;
				newConvertible=NULL;

				CBOptionsIn.clear();

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			};

			result.setLong(convertId);

			return ARM_OK;
		}
		else
		{
			Convertible=(ICM_Convertible*)LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Convertible,ICM_CONVERTIBLE)==0)
			{
				if(newConvertible)
					delete newConvertible;
				newConvertible=NULL;

				CBOptionsIn.clear();

				result.setMsg("ARM_ERR: Previous object is not of a good type");
				return ARM_KO;
			}
			else
			{
				if(Convertible)
					delete Convertible;
				Convertible=NULL;

				CBOptionsIn.clear();

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newConvertible,objId);

				return ARM_OK;
			};

		};
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if(newConvertible)
			delete newConvertible;
		newConvertible=NULL;

		CBOptionsIn.clear();

		ARM_RESULT();
	};
	
};

extern long ICMLOCAL_CBCALLABLEASW(long		ConvertId,
								 const VECTOR<long>& CBOptionIds,
								 double		SwapStart,
								 long		irIndexId,
								 double		Spread,
								 double		EarlyCall,
								 double		RecallSpread,
								 bool		hasEuribidClause,
								 bool		CallIfITM,
								 ARM_result&result,
								 long		objId)
{
	long cbaswId;
	int i=0;

	ARM_result C_result;
	long retCode=ARM_KO;

	ICM_CbAsw* CbAsw=NULL;
	ICM_CbAsw* newCbAsw=NULL;

	ICM_Convertible* Convertible=NULL;

	ARM_IRIndex* IrIndex=NULL;

	ICM_CBOption* pCBOption=NULL;
	vector<ICM_CBOption*> CallFeatures;

	char* myStartDate=new char[11];
	char* myEarlyCallDate=new char[11];
	ARM_Date SwapStartDate;
	ARM_Date EarlyCallDate;

	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS)==0)
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	};

	CCString msg("");

	try
	{
		Convertible=(ICM_Convertible*)LOCAL_PERSISTENT_OBJECTS->GetObject(ConvertId);

		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Convertible,ICM_CONVERTIBLE)==0)
		{
			result.setMsg("ARM_ERR: Convertible is not of a good type");
			return ARM_KO;
		};

		IrIndex=(ARM_IRIndex*)LOCAL_PERSISTENT_OBJECTS->GetObject(irIndexId);
		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IrIndex,ARM_IRINDEX)==0)
		{
			result.setMsg("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		};

		int n=CBOptionIds.size();
		for(i=0;i<n;i++)
		{
			pCBOption=(ICM_CBOption*)LOCAL_PERSISTENT_OBJECTS->GetObject(CBOptionIds[i]);
			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pCBOption,ICM_CBOPTION)==0)
			{
				result.setMsg("ARM_ERR: Call Feature is not of a good type");
				return ARM_KO;
			};
			CallFeatures.push_back(pCBOption);
		};

		Local_XLDATE2ARMDATE(SwapStart,myStartDate);
		SwapStartDate=(ARM_Date)myStartDate;

		if(EarlyCall==0)
		{
			EarlyCallDate=((ARM_Date)myStartDate).AddMonths(6);
		}
		else
		{
			Local_XLDATE2ARMDATE(EarlyCall,myEarlyCallDate);
			EarlyCallDate=(ARM_Date)myEarlyCallDate;
		};

		newCbAsw=new ICM_CbAsw(Convertible,SwapStartDate,IrIndex,Spread,CallFeatures,RecallSpread,EarlyCallDate,hasEuribidClause,CallIfITM);

		if(objId==-1)
		{
			CREATE_GLOBAL_OBJECT();
			cbaswId=LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCbAsw);
			if(cbaswId==RET_KO)
			{
				
				if(newCbAsw)
					delete newCbAsw;
				newCbAsw=NULL;

				CallFeatures.clear();

				if(myStartDate)
					delete[] myStartDate;
				myStartDate=NULL;

				if(myEarlyCallDate)
					delete[] myEarlyCallDate;
				myEarlyCallDate=NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			};

			if(CbAsw)
				delete CbAsw;
			CbAsw=NULL;

			CallFeatures.clear(); 

			if(myStartDate)
				delete[] myStartDate;
			myStartDate=NULL;

			if(myEarlyCallDate)
				delete[] myEarlyCallDate;
			myEarlyCallDate=NULL;

			result.setLong(cbaswId);
			return ARM_OK;
		}
		else
		{
			CbAsw=(ICM_CbAsw*)LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CbAsw,ICM_CBASW)==0)
			{
				if(newCbAsw)
					delete newCbAsw;
				newCbAsw=NULL;

				CallFeatures.clear();

				if(myStartDate)
					delete[] myStartDate;
				myStartDate=NULL;

				if(myEarlyCallDate)
					delete[] myEarlyCallDate;
				myEarlyCallDate=NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			if(CbAsw)
				delete CbAsw;
			CbAsw=NULL;

			CallFeatures.clear();

			if(myStartDate)
				delete[] myStartDate;
			myStartDate=NULL;

			if(myEarlyCallDate)
				delete[] myEarlyCallDate;
			myEarlyCallDate=NULL;

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCbAsw,objId);
			return ARM_OK;
		};
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if(newCbAsw)
			delete newCbAsw;
		newCbAsw=NULL;

		CallFeatures.clear(); //Call Feature is cleared because construction failed

		if(myStartDate)
			delete[] myStartDate;
		myStartDate=NULL;

		if(myEarlyCallDate)
			delete[] myEarlyCallDate;
		myEarlyCallDate=NULL;

		ARM_RESULT();
	};

};
