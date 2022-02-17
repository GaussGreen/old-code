#pragma warning(disable : 4541)
#pragma warning(disable : 4250)

#include "firstToBeIncluded.h"
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 


#include "ARM_local_glob.h"
#include <glob\armglob.h>
#include <inst\security.h>

#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_nthtd.h>

#include <libCCdate\CCdate.h>

#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>

#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_parsexml_for_icm.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <util\fromto.h> 


/** 
long ICMLOCAL_GetModelFromSummit (const long DiscountCurveId,
								  const CCString& idSummit,
								  const CCString& type,
								  const CCString& CurveId,
								  const CCString& CorrCurveId,
								  ARM_result& result,
								  long objId)
{

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	ARM_Model* mmc = NULL;
	ARM_Object* newObject = NULL;
	ARM_Object* oldObject = NULL;
	ARM_ZeroCurve * ircurve = NULL; 

	long Id;
	ARM_CLASS_NAME thisClass;

	CCString AssetType;

	if (type == "MMC")
	{
		AssetType = "CDO";
		thisClass = ICM_MODELMULTICURVES;
	}
	if (type == "CDSOPT")
	{
		AssetType = "CDSOPT";
		thisClass = ICM_DEFAULT_CURVE_MODEL;
	}
	else
		thisClass = ICM_MODELMULTICURVES;

	CCString msg (" ");

	// we catch the exceptions
	try
	{

		ircurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscountCurveId);

		 if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ircurve, ARM_ZERO_CURVE) == 0)
			{
			result.setMsg ("ARM_ERR: PWC Curve is not of a good type");
			return ARM_KO;
			}

		CCString xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,AssetType);

		switch (thisClass)
		{
		case ICM_MODELMULTICURVES :
			{
			newObject = ICMLOCAL_ParseModelCredit(ircurve, 
											  xmlResponse,
											  idSummit,
											  type,
											  CurveId,
											  CorrCurveId);
			break;
			}
		case ICM_DEFAULT_CURVE_MODEL :
			{
			newObject = ICMLOCAL_ParseModelForSpreadOptions(ircurve, 
															xmlResponse);
			break;
			}
		default:
			;
		}

		if (newObject == NULL)
		{
			result.setMsg ("ARM_ERR: MMC is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObject);

			if (Id == RET_KO)
			{
				if (newObject)
					delete newObject;
				newObject = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(Id);

			return ARM_OK;
		}
		else
		{
			if (type == "MMC")
			{
				mmc = (ICM_ModelMultiCurves*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
				oldObject = mmc;
			}

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldObject, thisClass) == 1)
			{
				if (oldObject)
				{
					delete oldObject;
					oldObject = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObject, objId);

				return ARM_OK;
			}

			else
			{
				if (newObject)
					delete newObject;
				newObject = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	// process the exceptions 
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newObject)
			delete newObject;
		newObject = NULL;

		ARM_RESULT();
	}

}
**/ 

//	-----------------------------------------------------------------------------------------
void 
ICMLOCAL_GetBasketCorrelMkDataFromCalypso(const std::string& pricingEnv,const ARM_Date& date,
										  const std::string& forceCurveName,const std::string& xmlFileName,
										  std::vector<std::string>& matus,
										  std::vector<double>&tenors,
										  ICM_QMatrix<double>& correls)
{

	std::string xmlOutput; 
	ARM_CalypsoToolkit::GetBasketCorrel(forceCurveName,pricingEnv,date,xmlFileName,xmlOutput); 
	ICMLOCAL_ParseBasketCorrelMkDataFromCalypso(xmlOutput,matus,tenors,correls); 
}
long ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(const std::string& pricingEnv,const ARM_Date& date,
												  const std::string& forceCurveName,const std::string Ccy,
												  const std::string& xmlFilename,long IndexId ,long objId)
{

	std::string xmlOutput;
	ARM_CalypsoToolkit::GetBasketCorrel(forceCurveName,pricingEnv,date, xmlFilename, xmlOutput); 

	std::vector<std::string> matus;
	std::vector<double>tenors;
	ICM_QMatrix<double> correls;
	ICMLOCAL_ParseBasketCorrelMkDataFromCalypso(xmlOutput,matus,tenors,correls); 
	// getting year fraction from maturities

	// here we have the correl datas 
	ARM_result C_result;
	long volId;

	ARM_VolLInterpol* newVolCrv = NULL;

	ARM_Vector* vMatu = NULL;
	ARM_Vector* vStrikes = NULL;
	ARM_Matrix* mVol = NULL;
	long strikeType;
	long volType;
	ARM_Date* AsOf = new ARM_Date(date);
	ARM_Date* tmpdate ;
	vMatu = new ARM_Vector(matus.size());

	for (int k=0;k<matus.size();k++)
	{

		tmpdate =new ARM_Date(matus[k].c_str(),"YYYYMMDD");

		vMatu->Elt(k) = ( tmpdate->GetJulian() - AsOf->GetJulian()) /365.;
		if (tmpdate)
			delete tmpdate;
	}
	
	vStrikes = new ARM_Vector(tenors.size());

	for ( k =0;k<tenors.size();k++)
	{
		vStrikes ->Elt(k) = tenors[k]/100;
	}
	//vStrikes = CreateARMVectorFromVECTOR(tenors);

	ARM_Currency* currency = new ARM_Currency(Ccy.c_str());
	double * pdVols = NULL;

	pdVols = new double[correls.Getnbrows()*correls.Getnbcols()];
	for (int i=0 ; i< correls.Getnbrows();i++)
	{
		for (int j=0;j<correls.Getnbcols();j++)
		{
			pdVols[j+i*correls.Getnbcols()] = correls(i,j);
		}
	}

	mVol = new ARM_Matrix(correls.Getnbrows(),correls.Getnbcols(),pdVols);
	
	if (pdVols)
		delete [] pdVols;
	pdVols = NULL;

	volType = ARM_ConvVolType("SMILE", C_result) ;
	strikeType = StrikeCode("S", C_result);

	newVolCrv = new ARM_VolLInterpol((ARM_Date) date, vMatu,
									     vStrikes, mVol, strikeType, 
                                         volType,
                                         currency);

	newVolCrv ->SetInterpType(K_LINEAR);

	ARM_IRIndex* theIndex; 
	LocalPersistent::get().convert(IndexId,theIndex); 
	if (theIndex) newVolCrv->SetIndex(*theIndex); 

	newVolCrv ->SetIndexName("");
	// setting matus labels
	
	for ( i=0; i<ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
	}

	for (i=0; i<matus.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsX[i], (const char *) matus[i].c_str());
	}

	if (objId=-1)
	{
		CREATE_GLOBAL_OBJECT();

		volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv);
	}
	else
	{
		ARM_VolLInterpol* VolCrv = (ARM_VolLInterpol*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (VolCrv)
			delete VolCrv;
		VolCrv=NULL;

		LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv, objId);
	}

	return volId;
	
}