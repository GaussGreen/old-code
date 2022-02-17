

#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <crv\zerocurv.h>
#include <ICMKernel\crv\icm_defaultcurve.h>
#include <ICMKernel\mod\modelmulticurves.h>
#include <ICMKernel\glob\icm_mktdatamng.h>
#include <ICMKernel\mod\icm_defcurvemodel.h>
#include <ICMKernel\mod\icm_meta_model.h>
#include <ICMKernel\inst\icm_gen.h>
#include "ICMKernel\mod\icm_customized_credit_multicurves.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ICMKernel\glob\icm_correlation.h>
#include "gpinflation\infcurv.h"
#include <memory>

using namespace ARM;
using namespace std;

long ICMLOCAL_ModelMultiCurves (int NbDefCurves, 
								VECTOR<long> DefCurvesID,
								long DiscCurveId,
								VECTOR<double> RecoveryRates ,
								long CorrelationId,
								long idvolcurve,
								bool CloneOrNot, 
								long idCpnInflationCurve,
								long idCpnIRCurve,
								ARM_result& result, 
								long objId)
{
	long modId = 0;
	int i = 0, j=0;
	int nbdefcurves = DefCurvesID.size();
	int sizereal = nbdefcurves;
	int modif = 0;

	ARM_ZeroCurve* pZCPY = NULL;
	ICM_ModelMultiCurves* prevModel = NULL;
	ICM_ModelMultiCurves* pmccm = NULL;
	// ICM_DefaultCurve** pZCCDS = NULL;
	ICM_DefaultCurve* pCdsZC = NULL;
	ICM_Correlation* pCorrelation = NULL;
	ARM_VolCurve* pVolCurve = NULL;
	ARM_InfCurv* pInfCurv = NULL;
	ARM_ZeroCurve* pCpnIRCurve = NULL;

	ARM_Vector pRecoveryRates ; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idvolcurve);
		pInfCurv = (ARM_InfCurv *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCpnInflationCurve);
		pCpnIRCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCpnIRCurve);

		pZCPY = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscCurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pZCPY, ARM_ZERO_CURVE) == 0) 
		{
			result.setMsg ("ARM_ERR: Zero Curve is not of a good type");
			return ARM_KO;
		}

		LocalPersistent::get().convert(CorrelationId,pCorrelation); 

		for (j = 0; j < nbdefcurves; j++) {if (DefCurvesID[j] == -1) sizereal-=1;}

		std::vector<const ICM_DefaultCurve*> pZCCDS (sizereal) ; 
		// ICM_DefaultCurve** pZCCDS = new ICM_DefaultCurve*[sizereal];
		
		modif = 0;
		for (j = 0; j < nbdefcurves; j++)
		{
			if (DefCurvesID[j] == -1) continue;

			// on recupere les courbes ****************************************************************
			pCdsZC = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DefCurvesID[j]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCdsZC, ICM_DEFAULTCURVE) == 0)
			{
				result.setMsg ("ARM_ERR: one CDS Zc Curve is not of a good type");
			
				return ARM_KO;
			}

			pZCCDS[modif] = pCdsZC;
			modif++;
		}	

		if (RecoveryRates.size()==0)
			pRecoveryRates.clear(); 
		else if (RecoveryRates[0] == -999.)
			pRecoveryRates.clear(); 
		else
		{
			pRecoveryRates.Resize(RecoveryRates.size()); 
			modif = 0;
			for	(j=0; j<RecoveryRates.size(); j++)
			{ if (DefCurvesID[j] != -1) pRecoveryRates[modif] = RecoveryRates[j];modif++; }
		}

		pmccm = new ICM_ModelMultiCurves(// sizereal,
								   pZCCDS,
								   pZCPY,
								   pRecoveryRates,
								   pCorrelation,
								   pVolCurve,
								   CloneOrNot,
								   pInfCurv,
								   pCpnIRCurve);



		if (pmccm == NULL)
		{
			result.setMsg ("ARM_ERR:ACVMC Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_ModelMultiCurves *)pmccm);

			if (modId == RET_KO)
			{
				if (pmccm)
					delete pmccm;
				pmccm = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* ModMod = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prevModel = (ICM_ModelMultiCurves*) (ModMod);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ModMod, ICM_MODELMULTICURVES ) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_ModelMultiCurves*)(pmccm), objId);

				return ARM_OK;
			}
			else
			{
				if (pmccm)
					delete pmccm;
				pmccm = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (pmccm)
			delete pmccm;
		pmccm = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_ModelMultiCurves ( int NbDefCurves, 
									VECTOR<long> DefCurvesID,
									long DiscCurveId,
									VECTOR<double> LossRates ,
									long CorrelationId,
									long marketDataMngId,
									long idvolcurve,
									bool  CloneOrNot,
									long objId = -1)
{
	long modId = 0;
	int i = 0, j=0;
	int nbdefcurves = DefCurvesID.size();
	int sizereal = nbdefcurves;
	int modif = 0;

	ICM_ModelMultiCurves* pmccm = NULL;
	
	ARM_ZeroCurve* pZCPY = NULL;	
	
	ICM_DefaultCurve* pZCCDS = NULL;
	ARM_VolCurve* pVolCurve = NULL;
	ICM_MktDataMng* aICM_MktDataMng = NULL;
	// ICM_DefaultCurve** PpZCCDS = NULL;
	std::vector<const ICM_DefaultCurve*> PpZCCDS ; 
	ICM_Correlation* pCorrelation = NULL;

	// double* pRecoveryRates = NULL;
	ARM_Vector pRecoveryRates ; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves : Pb with accessing objects");
	}
	try {

		pCorrelation = dynamic_cast<ICM_Correlation *>(LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelationId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCorrelation, ICM_CORRELATION) == 0) 
		{
			ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves : CorrelationId is not of the good type");
		}

		pVolCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(idvolcurve));
		aICM_MktDataMng = dynamic_cast<ICM_MktDataMng *>( LOCAL_PERSISTENT_OBJECTS->GetObject(marketDataMngId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(aICM_MktDataMng, ICM_MKTDATAMNG) ==0)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves : MktDataMng is not of a good type");
		}

		pZCPY = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pZCPY, ARM_ZERO_CURVE) == 0) 
		{
			ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves : Zero Curve is not of a good type");
		}

		
			for (j = 0; j < nbdefcurves; j++) {if (DefCurvesID[j] == -1) sizereal-=1;}
			// PpZCCDS = new ICM_DefaultCurve*[sizereal];
			PpZCCDS.resize(sizereal); 
			
			
			modif = 0;
			for (j = 0; j < nbdefcurves; j++)
			{
				if (DefCurvesID[j] == -1) continue;
				// on recupere les courbes ****************************************************************
				pZCCDS = dynamic_cast<ICM_DefaultCurve *>((dynamic_cast<ICM_DefaultCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(DefCurvesID[j])))->Clone());
				if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pZCCDS, ICM_DEFAULTCURVE) == 0)
				{
					ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves : one CDS Def Curve is not of a good type");
				}

				PpZCCDS[modif] = pZCCDS;
				modif++;
			}	

			if (LossRates.size()!=0 && LossRates[0] != -999.)
			{	
				if (LossRates.size() != 0)	pRecoveryRates.Resize(sizereal);
				modif = 0;
				for	(j=0; j<LossRates.size(); j++)
				{ 
					if (DefCurvesID[j] != -1) pRecoveryRates[modif] = LossRates[j];
					modif++; 
				}
			}

			pmccm = new ICM_ModelMultiCurves(// sizereal,
									   PpZCCDS,
									   pZCPY,
									   pRecoveryRates,
									   pCorrelation,
									   aICM_MktDataMng,
									   pVolCurve,
									   CloneOrNot);


		if (! pmccm)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves can't construct a ModelMulticurve with MktDataMng");
		}
		long QId = LocalPersistent::get().adopt(pmccm,objId );
		//if (pRecoveryRates) { delete []  pRecoveryRates; pRecoveryRates = NULL; }
		if (CloneOrNot) { // ="Y"
			// if (PpZCCDS)
			// {
				for (i=0;i<NbDefCurves;i++)
				{
					if (PpZCCDS[i])
						delete  (PpZCCDS[i]);
				}
				// delete[] PpZCCDS;
				// PpZCCDS = NULL;
			// }
		}
		return QId;
	}catch(Exception& )
	{
		// if (pRecoveryRates) { delete [] pRecoveryRates; pRecoveryRates = NULL; }
		if (CloneOrNot) { // ="Y"
			// if (PpZCCDS)
			// {
				for (i=0;i<NbDefCurves;i++)
				{
					if (PpZCCDS[i]){
						delete  PpZCCDS[i]; PpZCCDS[i] = NULL;}
				}
			// 	delete[] PpZCCDS;
			// 	PpZCCDS = NULL;
			// }
		}
	//	x.DebugPrint();

	//	ARM_RESULT();
	}
	 catch(...) // only for pRecoveryRates that can't be construct with a auto_ptr
	{
		// if (pRecoveryRates) { delete [] pRecoveryRates; pRecoveryRates = NULL; }
		if (CloneOrNot) { // ="Y"
			// if (PpZCCDS)
			// {
				for (i=0;i<NbDefCurves;i++)
				{
					if (PpZCCDS[i]){
						delete  PpZCCDS[i]; PpZCCDS[i] = NULL;}
				}
				// delete[] PpZCCDS;
				// PpZCCDS = NULL;
			// }
		}
		ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_ModelMultiCurves fatal error");
	}
	 return ARM_KO;
}


long ICMLOCAL_SetPropMixCopule (long ModelId, double PropIndep,double PropFullCorrel, ARM_result& result)
{
	long modId = 0;
	ICM_ModelMultiCurves* pmccm = NULL;

	ICMMSG(WARN,"Using ICMLOCAL_SetPropMixCopule"); 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pmccm = (ICM_ModelMultiCurves *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pmccm, ICM_MODELMULTICURVES) == 0) 
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		pmccm->SetIndependantPart(PropIndep);
		pmccm->SetFullCorrelPart(PropFullCorrel);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetPropMixCopule : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_DefProbModel ( long idDefProb, 
							 long idIRcurve,
							 long idvolcurve,
							 ARM_result& result, 
							 long objId)
{
	long modId =0;

	ICM_DefaultCurveModel* Model = NULL;
	ICM_DefaultCurveModel* prevModel = NULL;

	ARM_ZeroCurve* IRCurve = NULL;
	ICM_DefaultCurve* pDefProbCurve=NULL;
	ARM_VolCurve* pVolCurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		// on recupere les courbes ****************************************************************
		pDefProbCurve = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idDefProb);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pDefProbCurve, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Curve is not of a good type");
			
			return ARM_KO;
		}

		IRCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idIRcurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(IRCurve, ARM_ZERO_CURVE) == 0)
		{

			result.setMsg ("ARM_ERR: Zero Curve is not of a good type");
			return ARM_KO;
		}

		pVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idvolcurve);

		// creation du model de default ****************************************************************
		Model = new ICM_DefaultCurveModel(pDefProbCurve,IRCurve,pVolCurve);

		if (Model == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}


		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Model);

			if (modId == RET_KO)
			{
				if (Model)
					delete Model;
				Model = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ICM_DefaultCurveModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ICM_DEFAULT_CURVE_MODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Model, objId);

				return ARM_OK;
			}
			else
			{
				if (Model)
					delete Model;
				Model = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (Model)
			delete Model;
		Model = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_MetaModel (int Nbmodel, 
						 const VECTOR<long>& ModelID,
						 VECTOR<int>& PricerType,
						 ARM_result& result, 
						 long objId)
{
	long modId = 0;
	int i = 0, j=0;
	int nbdefcurves = ModelID.size();
	int sizereal = nbdefcurves;
	int modif = 0;

	ICM_Meta_Model* prevModel = NULL;
	ICM_Meta_Model* pmccm = NULL;
	ARM_Model** pZCCDS = NULL;
	ARM_Model* pCdsZC = NULL;

	VECTOR<ARM_CLASS_NAME> PricerTypeReal;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		for (j = 0; j < nbdefcurves; j++) {if (ModelID[j] == -1) sizereal-=1;}

		ARM_Model** pZCCDS = new ARM_Model*[sizereal];
		PricerTypeReal.resize(sizereal);
		
		modif = 0;
		for (j = 0; j < nbdefcurves; j++)
		{
			if (ModelID[j] == -1) continue;

			// on recupere les courbes ****************************************************************
			pCdsZC = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelID[j]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCdsZC, ARM_MODEL) == 0)
			{
				result.setMsg ("ARM_ERR: one Model is not of a good type");
			
				return ARM_KO;
			}

			pZCCDS[modif] = pCdsZC;
			PricerTypeReal[modif] = (ARM_CLASS_NAME) PricerType[j];

			modif++;
		}	

		pmccm = new ICM_Meta_Model(pZCCDS,
								   sizereal,
								   PricerTypeReal);


		if (pmccm == NULL)
		{
			result.setMsg ("ARM_ERR:ACVMC Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Meta_Model *)pmccm);

			if (modId == RET_KO)
			{
				if (pmccm)
					delete pmccm;
				pmccm = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* ModMod = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prevModel = (ICM_Meta_Model*) (ModMod);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ModMod, ICM_META_MODEL ) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Meta_Model*)(pmccm), objId);

				return ARM_OK;
			}
			else
			{
				if (pmccm)
					delete pmccm;
				pmccm = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (pmccm)
			delete pmccm;
		pmccm = NULL;

		ARM_RESULT();
	}
}
/*
long ICMLOCAL_MarketDataMng (const	VECTOR<long>& Data,
									long objId = -1)
{
	long modId = 0;
	int i = 0, j=0;
	int nbdatas = Data.size();
	int sizereal = nbdatas;
	int modif = 0;

	ICM_MktDataMng* pMktDataMng = new ICM_MktDataMng;

	
	if (! pMktDataMng)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_MarketDataMng can't construct a pMktDataMng");
	}

	//for (j = 0; j < nbdatas; j++) {if (Data[j] == -1) sizereal-=1;}
	for (j = 0; j < nbdatas; j++)
	{
		if (Data[j] == -1) continue;
		// on recupere les courbes ****************************************************************
		ARM_Object* object = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(Data[j]);
		// creation d'une clef particulière pour que cela fonctionne 
		//	avec la recuperation figée du modelmulticurves.
		ARM_CLASS_NAME name = object->GetName();
		char TMP[1000];
		switch (name)
		{
		case ICM_DEFAULTCURVE:
		case ICM_CST_PIECEWISE:
		case ICM_LINEAR_PIECEWISE:
		case ICM_DEFAULT_CURVE_MODEL:
		{//DEFPROB_MO_20050302
		ICM_DefaultCurve* curve = (ICM_DefaultCurve*) object;
		string plabel = curve->GetLabel();
		char* date = curve->GetAsOfDate().GetStrDate();
		sprintf(TMP,"DEF:%s_%s",plabel.c_str(),date);
		delete[] date;
		}
		break;
		case ARM_ZERO_CURVE:
		case ARM_ZERO_INTERPOLATION:
		case ARM_ZERO_LIN_INTERPOL:
		case ARM_ZERO_CUBDIFF:
		case ARM_ZERO_SPLINES:
		case ARM_ZERO_SPLICUB:
		case ARM_ZERO_VASICEK:
		case ARM_ZERO_FLAT:
		{//EUR_EURIB_MO_20050302
			ARM_ZeroCurve* curve = (ARM_ZeroCurve*) object;
			string index(IRCURVE); 
			sprintf(TMP,"%s",index.c_str());
		}
		break;
		case ICM_CORRMATRIX:
		case ICM_SMILE_CORRMATRIX:
		case ICM_BETA_CORRMATRIX:
		{//LABEL_MO_20050302
			ICM_Smile_Correlation* curve = (ICM_Smile_Correlation*) object;
			string structname(CORREL);
			sprintf(TMP,"%s",structname.c_str());
		}
		break;
		
		case ARM_INFCURV :
		{	ARM_InfCurv* infla = dynamic_cast<ARM_InfCurv*>(object);
			string structname(INFCURV);
			sprintf(TMP,"%s",structname.c_str());
		}
		break;
		case ICM_IMPLIED_LOSS_TREE :
		{//LABEL_MO_20050302
			ICM_ImpLossTree* tree = (ICM_ImpLossTree*) object;
			ICM_Smile_Correlation* correl = tree->GetCorrelation();
			string structname = correl->GetStructName();
			char* date = correl->GetAsOfDate().GetStrDate();
			sprintf(TMP,"TREE:%s_%s",structname.c_str(),date);
			delete[] date;
		}
		break;
		case ARM_SMILE_CURVE:
		case ARM_VOL_CURVE:
		case ARM_VOL_FLAT:
		case ARM_VOL_LIN_INTERPOL:
		case ARM_VOL_CUBE:
		{//EUR_MO_20050302
			ARM_VolCurve* curve = (ARM_VolCurve*) object;
			char* name = curve->GetIndexName(); 
			char* ccy = curve->GetCurrency()->GetCcyName();
			char* date = curve->GetAsOfDate().GetStrDate();
			sprintf(TMP,"VOL:%s_%s_%s",name,ccy,date);
			delete[] date;
		}
		break;
			
		}
		string label(TMP);
		pMktDataMng->associate(object, label);
	}	

	long QId = LocalPersistent::get().adopt(pMktDataMng,objId );

	return QId; 
}

*/
long ICMLOCAL_MarketDataMng (const	VECTOR<long>& Data,
									long objId = -1)
{
	long modId = 0;
	int i = 0, j=0;
	int nbdatas = Data.size();
	int sizereal = nbdatas;
	int modif = 0;

	ICM_MktDataMng* pMktDataMng = new ICM_MktDataMng;

	
	if (! pMktDataMng)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_MarketDataMng can't pMktDataMng");
	}

	//for (j = 0; j < nbdatas; j++) {if (Data[j] == -1) sizereal-=1;}
	for (j = 0; j < nbdatas; j++)
	{
		if (Data[j] == -1) continue;
		// on recupere les courbes ****************************************************************
		ARM_Object* pObj = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(Data[j]);
		pMktDataMng->associate(pObj);
	}	

	long QId = LocalPersistent::get().adopt(pMktDataMng,objId );

	return QId; 
}

long ICMLOCAL_Customized_Credit_MultiCurves(
									long				DiscountCurveId,
									VECTOR<CCString>	labels,
									vector<double>		spreads,
									VECTOR<CCString>	maturities,								
									long				Data_DescriptionId,
									long				Market_ParametersId,
									long				CDO_Square_ParametersId,
									long				CDO_Square_DataId,
									long				CorrelationId,
									long				ModelMultiCurvesId,
									ARM_result& result, 
									long objId)
{
	long	modId	=	0;
	double	dResult	=	0.;

	ARM_ZeroCurve*			pZCPY = NULL;

	ARM_Model*				pModel = NULL;

	ICM_Correlation*		pCorrelation	=	NULL;
	ICM_ModelMultiCurves*	pMultiCurves	=	NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Description = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_CDOSquare_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_CDOSquare_Data = NULL;

	ICM_Customized_Credit_MultiCurves*	pccmc		=	NULL;
	ICM_ModelMultiCurves*				prevModel	=	NULL;

	ICM_GenCF* TheCashFlows = NULL;

	// ---------------------------------------------------------------------------------
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	// ---------------------------------------------------------------------------------

	CCString msg ("");

	int i,j,il;
	int nbissuers = 0;
	int nbmat= 0;
	char** pmaturities = NULL;

	vector<string> AllCreditLabels;
	vector<string> AllMaturities;

	try
	{
		// ---------------------------------------------------------------------------------
		// ZERO COUPON CURVE
		// ---------------------------------------------------------------------------------
		pZCPY = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscountCurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pZCPY, ARM_ZERO_CURVE) == 0) 
		{
			result.setMsg ("ARM_ERR: Zero Curve is not of a good type");
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// CORRELATION
		// ---------------------------------------------------------------------------------
		pCorrelation = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelationId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCorrelation, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Object is not of a good type");
			return ARM_KO;
		}
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		pModel = (ARM_Model*) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelMultiCurvesId);
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// LABELS
		// ---------------------------------------------------------------------------------
		nbissuers = labels.size();
		AllCreditLabels.resize(nbissuers);

		if (nbissuers)
		{
			for (j=0; j<nbissuers; j++)
				AllCreditLabels[j] = (string)(labels[j]);
		}
		/*
		char** psLabels = new char*[nbissuers];

		for (j=0;j<nbissuers;j++)
		{
			psLabels[j] = new char[60];
			sprintf(psLabels[j],(const char*)labels[j]);
		}*/
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// MATURITIES
		// ---------------------------------------------------------------------------------
		nbmat = maturities.size();
		AllMaturities.resize(nbmat);

		if (nbmat)
		{
			for (il=0; il<nbmat; il++)
			{
				AllMaturities[il]	=	(string)(maturities[il]);
			}
		}
		/*
		char** psMatu = new char*[nbmat];

		for (j=0;j<nbmat;j++)
		{
			psMatu[j] = new char[6];
			sprintf(psMatu[j],(const char*)maturities[j]);
		}*/		
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// SPREADS
		// ---------------------------------------------------------------------------------

		ICM_QMatrix<double>* CreditSpreads = new ICM_QMatrix<double>(nbissuers,nbmat,CREDIT_DEFAULT_VALUE);
			
		for (i=0;i<nbissuers;i++)
			for (j=0;j<nbmat;j++)
				CreditSpreads->SetValue(i,j, spreads[i*nbmat+j] / 10000.0);		// in bps.


		// ---------------------------------------------------------------------------------
		// Credit Data Description
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(Data_DescriptionId); 

		if (TheCashFlows)	Matrix_Description = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Credit Data Parameters
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(Market_ParametersId); 

		if (TheCashFlows)	Matrix_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// CDO^2 Credit Data Parameters
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CDO_Square_ParametersId); 

		if (TheCashFlows)	Matrix_CDOSquare_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// CDO^2 Credit Data Data
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CDO_Square_DataId); 

		if (TheCashFlows)	Matrix_CDOSquare_Data = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		if (pModel == NULL)
			pccmc	=	new ICM_Customized_Credit_MultiCurves(
													pZCPY,
													Matrix_Description,
													Matrix_CDOSquare_Parameters,
													Matrix_CDOSquare_Data,
													Matrix_Parameters,
													CreditSpreads,
													AllCreditLabels,
													AllMaturities,
													pCorrelation
								   );
		else
			pccmc	=	new ICM_Customized_Credit_MultiCurves((ICM_ModelMultiCurves*) pModel);

		// ---------------------------------------------------------------------------------

		if (pccmc == NULL)
		{
			result.setMsg ("ARM_ERR: Customized Multi Curves Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_ModelMultiCurves *)pccmc);

			if (modId == RET_KO)
			{
				if (pccmc)
					delete pccmc;
				pccmc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object: Customized Multi Curves");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* ModMod = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prevModel = (ICM_Customized_Credit_MultiCurves*) (ModMod);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ModMod, ICM_CUSTOMIZED_CREDIT_MULTICURVES ) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Customized_Credit_MultiCurves*)(pccmc), objId);

				return ARM_OK;
			}
			else
			{
				if (pccmc)
					delete pccmc;
				pccmc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type: Customized Multi Curves");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}

long ICMLOCAL_SetVolCurve(long ModelId,
							long idvolcurve,
							ARM_result& result)
{
	ARM_Model* Model=NULL;
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

		Model = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model, ARM_MODEL) == 0) 
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		ICM_DefaultCurveModel* DefModel = dynamic_cast<ICM_DefaultCurveModel*>(Model);
		if (!DefModel) {
			ICM_MktDataMng* MktModel = dynamic_cast<ICM_MktDataMng*>(Model);
			if(!MktModel) {
				result.setMsg ("ARM_ERR: Model is nor ICM_DefaultCurveModel nor ICM_MktDataMng ");
				return ARM_KO;
			}
			// Clone VolCurve
			MktModel->adopt(pVolCurve);
		} else {
			DefModel->SetVolCurve(pVolCurve) ;
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
		result.setMsg ("ARM_ERR: ICMLOCAL_SetVolatility : unrecognized failure");
		return ARM_KO;
	}

}