//disable STL warnings in Windows NT

#include "firstToBeIncluded.h"
#include "ARM_local_glob.h"

#include <glob\armglob.h>
#include <inst\security.h>
#include <ICMKernel\inst\icm_cds.h>
#include <ICMKernel\pricer\icm_pricer_basket.h>
#include <ICMKernel\inst\icm_option.h>
#include <ICMKernel\glob\icm_corrmatrix.h>
#include <ICMKernel\glob\icm_smile_correlation.h>
#include <ICMKernel\glob\icm_index_correlation.h>
#include <ICMKernel\util\icm_schedule_info.h>
#include <ICMKernel\glob\icm_flatcorrel.h>

#include <ICMKernel\util\icm_superseniorleverage.h>
#include "ICMKernel/pricer/icm_pricer_analytic_cdo2_smile.h"
#include "ICMKernel\inst\icm_collateral.h"

/// remove warnings on va_start and va_end
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
/// and do not change the order as this one has to be precisely here
/// if you do not know what you are doing, please ask someone who does!
#include <ARM\libarm_local\undef_va_vars.h>
#include "ICMKernel\pricer\icm_pricer_adviser.h"
// #include "ICMKernel\util\icm_ftdbacktest.h"

#include <libCCdate\CCdate.h>

#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>
#include <ARM\libarm_frometk\arm_local_etoolkit_for_ICM.h>
#include <ARM\libarm_frometk\arm_local_parsexml_for_ICM.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <ARM\libicm_local\ICM_local_pwccurve.h>
#include <ARM\libicm_local\ICM_local_leg.h>
#include "ICMKernel\inst\icm_gen.h"
#include "ICMKernel\glob\icm_correlation_sector.h"
#include <ICMKernel/pricer/icm_pricer_mc_cdo2.h>
#include <ICMKernel/pricer/icm_pricer_option.h>
#include <ICMKernel/pricer/icm_pricer_tree_binomial_spreadoption.h>
#include <ICMKernel\glob\icm_maths.h>
#include "ICMKernel\crv\icm_distriblosscalculator.h"
#include <ICMKernel\glob\icm_cubicspline.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_credit_index.h>
#include <ICMKernel\pricer\icm_pricer_mc_cdo.h>
#include <ICMKernel\glob\icm_calibrator.h>
#include <ICMKernel\random\icm_RandomGenerator.h>
#include <ICMKernel\random\icm_RandomNag.h>
#include <ICMKernel\random\icm_RandomRan1.h>
#include <ICMKernel\random\icm_RandomRan2.h>
#include <ICMKernel\random\icm_RandomRanDef.h>
#include <ICMKernel\random\icm_RandomRanmar.h>
#include <ICMKernel\random\icm_randomKISS.h>
#include <ICMKernel\random\icm_randomLaws.h>
#include <ICMKernel\random\icm_randomInvNormAcklam.h>
#include <ICMKernel\random\icm_randomInvNormMoro.h>
#include <ICMKernel\\random\icm_RandomRNG-Str.h>

#include <stdio.h>
#define itoa _itoa

std::string ICMLOCAL_Version()
{
	return VersionInfoSingleton().longversion(); 
}


long ARM_Credit_DisplayScheduleDates (long instId,
										long datesType,
										long recId,
										ARM_result& result)
{
	ARM_Security* sec = NULL;

	// ARM_SwapLeg* leg = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);

		ARM_SwapLeg* swapLeg = dynamic_cast<ARM_SwapLeg*>(sec); 
		if (swapLeg) { swapLeg->DisplayScheduleDates(datesType,K_YES,"123");  return ARM_OK ; }
		ICM_Cds* cds = dynamic_cast<ICM_Cds*>(sec); 
		if (cds) { cds->DisplayScheduleDates(datesType,"123"); return ARM_OK ; } 
		result.setMsg ("ARM_ERR: function not implemented for this object");
		return ARM_KO;

	
    }
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	return ARM_OK;
}

long ARM_Credit_DisplayScheduleValues (long pricerId,
										 long valuesType,
										 long recId,
										 ARM_result& result)
{
	ICM_Pricer* pricer = NULL;
	ARM_Security* sec = NULL;

	// ARM_SwapLeg* leg = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		pricer = (ICM_Pricer*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}
		
		pricer->Price(qCMPPRICE);

		sec = pricer->GetSecurity();

		ARM_SwapLeg* swapLeg = dynamic_cast<ARM_SwapLeg*>(sec); 
		if (sec) { swapLeg->DisplayScheduleValues(valuesType,"123"); return ARM_OK ;}
		ICM_Cds * cds = dynamic_cast<ICM_Cds*>(sec);  
		if (cds) { cds->DisplayScheduleValues(valuesType,"123"); return ARM_OK ;}
		
		result.setMsg ("ARM_ERR: function not implemented for this object");
		return ARM_KO;

	
    }
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	catch (...)
	{
		result.setMsg ("ARM_ERR: DisplaySV : unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


long ICMLOCAL_Spread(long pricerId, double MtM, ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	double price=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ARM_Date NewAsOfDate;

		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		price = pricer->ComputeSpread(MtM);

		result.setDouble(price);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Spread : unrecognized failure");
		return ARM_KO;
	}

}



long ICMLOCAL_Price (long pricerId, double AsOfDate, ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	double price=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	char pAsOfDate[11];
		
	try
	{
		ARM_Date NewAsOfDate;

		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		if (AsOfDate == -1.0)
			strcpy(pAsOfDate,"NULL");
		else
		{
			Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate);
		}

		price = pricer->Price(qCMPPRICE);

		result.setDouble(price);

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


long ICMLOCAL_GetLabel (long curveId, ARM_result& result)
{
	ICM_DefaultCurve* crv = NULL;
	// char* label = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		crv = dynamic_cast<ICM_DefaultCurve*> ( LOCAL_PERSISTENT_OBJECTS->GetObject(curveId) );

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(crv, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Curve is not of a good type");
			return ARM_KO;
		}

		std::string label = crv->GetLabel();


		result.setMsgString(label);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: GetLabel : unrecognized failure");
		return ARM_KO;
	}

}



long ICMLOCAL_SetLabel (long curveId, CCString label, ARM_result& result)
{
	ICM_DefaultCurve* crv = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		char* tmp = (char*)label;

		crv = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(crv, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Curve is not of a good type");
			return ARM_KO;
		}

		crv->SetLabel(tmp);

		result.setMsg(label);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ICMLOCAL_DTR (long pricerId , int NbDefaults,CCString S_or_L , VECTOR<CCString>& Labels, VECTOR<double> RecoveryRates , ARM_result& result )
{
	ICM_Pricer_MC_Cdo2* pricer=NULL;
	double sensitivity=0.0;
	char** labels =  NULL;
	double* pRecoveryRates = NULL;


	int size = Labels.size();
	int i=0, j=0;

	CCString TmpStr;
	CCString msg ("");
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{

		labels = new char*[size];
		
		for (i=0; i<size; i++)
		{
			TmpStr = Labels[i];
			TmpStr.trim_right();
			labels[i] = TmpStr.GetStr();
		}

		pRecoveryRates = new double[size];
		
		for	(j=0; j<size; j++)
			pRecoveryRates[j] = RecoveryRates[j];
		
		
		pricer = (ICM_Pricer_MC_Cdo2*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		sensitivity = pricer->ComputeDTR(NbDefaults, S_or_L, labels,pRecoveryRates);		

		result.setDouble(sensitivity);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Sensitivity : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_Sensitivity (long pricerId ,
						   qSENSITIVITY_TYPE curvtype, 
						   const std::string&  plot, 
						   const std::string& label, 
						   double epsilon,
						   double epsilonGamma,
						   ARM_result& result)
{

	
// 17783 	if (curvtype==ICM_FAST_SPREAD_TYPE) ICMMSG(WARN,"Using ICMLOCAL_Sensitivity with  ICM_FAST_SPREAD_TYPE"); 
// 17783 	if (curvtype==ICM_FAST_RECOVERY_TYPE) ICMMSG(WARN,"Using ICMLOCAL_Sensitivity with  ICM_FAST_RECOVERY_TYPE"); 

	ICM_Pricer* pricer=NULL;
	double sensitivity=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		sensitivity = pricer->Hedge(curvtype,plot,label,epsilon, epsilonGamma);
		
		result.setDouble(sensitivity);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Sensitivity : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_DisplayMatrix (long pricerId, 
							   VECTOR<CCString>& VectorOut,
							   int& NbCol,  
							   VECTOR<CCString>& VectorName,
							   ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	ARM_Security* sec = NULL;

	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		pricer = (ICM_Pricer*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER)) 
		{
			retCode = ICMLOCAL_Price(pricerId,-1.,result);

			ICM_Matrix<ARM_Vector>* matrix = NULL; //(ICM_Matrix*) pricer->GetPricingMatrix();
			
			if (matrix == NULL)
			{
				result.setMsg ("ARM_ERR: No model assigned to this Security ! ");
				return ARM_KO;
			}

			int nbcols = matrix->GetNbCol();
			int nbrows = matrix->GetCol(0)->GetSize();
			NbCol = nbcols;

			for (int j = 0; j<nbcols; j++)
			{
				VectorName.push_back(matrix->GetColName(j));
				ARM_Vector* vect = matrix->GetCol(j);

				for (int i = 0; i<nbrows; i++)
				{
					char Tmp[30];
					if (j == 0)
					{
					ARM_Date TmpDate = (ARM_Date) vect->Elt(i);
					TmpDate.GetCompleteStrDate(Tmp);
					VectorOut.push_back(Tmp);
					}
					else
					{
					sprintf(Tmp,"%f",vect->Elt(i));
					VectorOut.push_back(Tmp);
					}

				}
			}	
		}
		else
		{
			result.setMsg ("ARM_ERR: function not implemented for this object");
			return ARM_KO;
		}
    }
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	catch (...)
	{
		result.setMsg ("ARM_ERR: DisplayMatrix : unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


//	----------------------------------------------------------------------------------------------
double 
ICMLOCAL_RiskyPV01 (long DefCurveId, const ARM_Date* date1, const ARM_Date& date2) 
{
	ICM_DefaultCurve* defprobcurve=NULL;
	LocalPersistent::get().convert(DefCurveId,defprobcurve); 
	if (!defprobcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RiskyPV01: no defcurve specified"); 
	if (date1) return defprobcurve->RiskyPV01(*date1,date2); 
	else return defprobcurve->RiskyPV01(defprobcurve->GetAsOfDate(),date2); 
}
//	----------------------------------------------------------------------------------------------
double 
ICMLOCAL_RiskyPV01 (long DefCurveId, const ARM_Date* date1, const std::string& tenor) 
{
	ICM_DefaultCurve* defprobcurve=NULL;
	LocalPersistent::get().convert(DefCurveId,defprobcurve); 
	if (!defprobcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RiskyPV01: no defcurve specified"); 
	if (date1) return defprobcurve->RiskyPV01(*date1,tenor); 
	else return defprobcurve->RiskyPV01(defprobcurve->GetAsOfDate(),tenor); 
}
//	----------------------------------------------------------------------------------------------
double 
ICMLOCAL_RiskyPV01AsSensitivity (long DefCurveId, const std::string& tenor) 
{
	ICM_DefaultCurve* defprobcurve=NULL;
	LocalPersistent::get().convert(DefCurveId,defprobcurve); 
	if (!defprobcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RiskyPV01: no defcurve specified"); 
	return defprobcurve->RiskyPV01AsSensitivity(tenor); 
} 
 
//	----------------------------------------------------------------------------------------------
double 
ICMLOCAL_RiskyDuration (long DefCurveId, const std::string& tenor) 
{
	ICM_DefaultCurve* defprobcurve=NULL;
	LocalPersistent::get().convert(DefCurveId,defprobcurve); 
	if (!defprobcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RiskyDuration: no defcurve specified"); 
	return defprobcurve->RiskyDuration(tenor); 
} 
//	----------------------------------------------------------------------------------------------
double 
ICMLOCAL_RiskyDuration (long DefCurveId, const ARM_Date& date) 
{
	ICM_DefaultCurve* defprobcurve=NULL;
	LocalPersistent::get().convert(DefCurveId,defprobcurve); 
	if (!defprobcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RiskyDuration: no defcurve specified"); 
	return defprobcurve->RiskyDuration(date); 
} 


long ICMLOCAL_RiskyDuration (long ModelId,CCString issuer, double date, CCString Tenor, ARM_result& result)
{
	const ICM_DefaultCurve* defprobcurve=NULL;
	char* pdate =new char[11];
	double duration = 0.;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ARM_Date NewAsOfDate;

		defprobcurve = ((ICM_ModelMultiCurves*)
		LOCAL_PERSISTENT_OBJECTS->GetObject(ModelId))->GetDefaultCurve((char*)(const char*)issuer);

		if (Tenor=="NONE")
		{
	 		Local_XLDATE2ARMDATE(date,pdate);
			ARM_Date edate = (ARM_Date) pdate;

			duration = defprobcurve->RiskyDuration(edate);
		}
		else
			duration = defprobcurve->RiskyDuration(CCSTringToSTLString(Tenor));

		result.setDouble(duration);

		if (pdate)
			delete[] pdate;
		pdate = NULL;

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		if (pdate)
			delete[] pdate;
		pdate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: RiskyDuration : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_Delivery (double AsOfDate,CCString maturity,ARM_result& result)
{

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	char* pAsOfDate=new char[11];
	char* nDate = new char[30];	
	ARM_Currency* Currency = NULL;
	int m, y;

	try
	{
		Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate);

		ARM_Date NewAsOfDate = (ARM_Date)(pAsOfDate);

		GetMonthYearFromExpiryDate((char*)maturity, &m ,&y);
		NewAsOfDate.Calc_MATIF_IRFutureDelivery(m, y);

		NewAsOfDate.GetCompleteStrDate(nDate);
		result.setString(nDate);

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


double ICMLOCAL_GetBeta (long pricerId, CCString label, ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	double beta=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ICM_Pricer_Basket* pricer = dynamic_cast<ICM_Pricer_Basket*>( LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId) );

		if (!pricer) ICMTHROW(ERR_INVALID_ARGUMENT,"Not an ICM_Pricer_Basket"); 

		beta = pricer->GetBeta((char*)label);

		result.setDouble(beta);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: GetBeta : unrecognized failure");
		return ARM_KO;
	}

}


double ICMLOCAL_GetDefProbTranche (long pricerId, double yearterm, ARM_result& result)
{
	ICM_Pricer_Security * pricer=NULL;
	double DPT=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pricer = dynamic_cast<ICM_Pricer_Security*> ( LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId) );

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		vector<double> losses;

		DPT = pricer->ExpectedLossTranche(yearterm,losses);

		result.setDouble(DPT);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: GetDefProbTranche : unrecognized failure");
		return ARM_KO;
	}

}


double ICMLOCAL_GetDuration (long pricerId, ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	double Durat=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		Durat = pricer->Price(qCMPDURATION);

		result.setDouble(Durat);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetDuration : unrecognized failure");
		return ARM_KO;
	}

}
long ICMLOCAL_CORRMATRIX (const ARM_Date& AsOf,
						  const string& Name,		
						  const std::vector<std::string>& labels,
						  const VECTOR<double>& coefs, 
						  ARM_result& result, 
						  long objId)
{

	ICM_CorrMatrix* PrevMatrix = NULL;
	ICM_CorrMatrix* Matrix = NULL;
	
	// char** labels2 = NULL;
	std::vector<std::string> labels2 ;
	// double** dmatrix = NULL;
	ICM_QMatrix<double> dmatrix; 
	int size ;
	int il=0,il2=0;

	// char* pAsOfDate=NULL;

	long CorrId = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	CCString TmpStr;

	try
	{


		size = labels.size();

		// labels2 = new char*[size];
		labels2.resize(size); 
		// dmatrix = new double*[size];
		dmatrix.Resize(size,size); 

		for (il=0; il<size; il++)
		{
		// dmatrix[il] = new double[size];
		TmpStr = labels[il].c_str();
		TmpStr.trim_right();
		labels2[il] = TmpStr.GetStr();	
		}

		for (il2=0;il2<size;il2++)
			for (il=0;il<size;il++)
			{
				if ((fabs(coefs[il+il2*size])>=1) && (il!=il2))
				{
				result.setMsg ("ARM_ERR: one correlation is >= 1");
				return ARM_KO;
				}

				dmatrix(il,il2) = 	coefs[il+il2*size];
			}


        for (il=0;il<size;il++)
		{
			for (il2=il;il2<size;il2++)
				if (dmatrix(il,il2) != dmatrix(il2,il) )
				{
					char buffer[20];
					itoa(il,buffer,9);
					msg = "Correlation Matrix is not symetric ! - pb (row,col)=(";
					msg += (CCString)buffer;
					itoa(il2,buffer,9);
					msg = msg + "," + (CCString)buffer+")";
					result.setMsg(msg);
					return ARM_KO;
				}
		}



		Matrix = new ICM_CorrMatrix(AsOf,Name, labels2, dmatrix );

	
		
		if (Matrix == NULL)
		{
			result.setMsg ("ARM_ERR: Correaltion Matrix is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix);

			if (CorrId == RET_KO)
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			PrevMatrix = (ICM_CorrMatrix*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevMatrix, ICM_CORRELATION) == 1)
			{
				if (PrevMatrix)
				{
					delete PrevMatrix;
					PrevMatrix = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix, objId);

				return ARM_OK;
			}
			else
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: CorrMatrix : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_EXTRACTCORRMATRIX (long CorrMatrixId,
								 const std::vector<std::string>& labels,
								 VECTOR<double>& coefsout, 
								 ARM_result& result)
{

	ICM_CorrMatrix* ExtMatrix = NULL;
	ICM_CorrMatrix* Matrix = NULL;

	coefsout.clear();
	// char** labels2 = NULL;
	std::vector<std::string> labels2; 
	double** dmatrix = NULL;
	int size = labels.size();
	int il=0,il2=0;

	long CorrId = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	CCString TmpStr;

	try
	{

		// labels2 = new char*[size];
		labels2.resize(size); 
		dmatrix = new double*[size];

		for (il=0; il<size; il++)
		{
		dmatrix[il] = new double[size];
		TmpStr = labels[il].c_str();
		TmpStr.trim_right();
		labels2[il] = TmpStr.GetStr();	
		}

		Matrix = (ICM_CorrMatrix *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrMatrixId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Matrix, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Matrix is not of a good type");
			return ARM_KO;
		}

		ExtMatrix = Matrix->ExtractCorrs(labels2);

		if (ExtMatrix == NULL)
			return ARM_KO;

		const ICM_QMatrix<double>& matrix = ExtMatrix->GetMatrix();

		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
				coefsout.push_back(matrix.Getvalue(i,j));

		if (ExtMatrix)
			delete ExtMatrix;
		ExtMatrix = NULL;

		if (dmatrix)
		{
		for (il=0;il<size;il++)
		{		
			delete[] dmatrix[il];
		}

		delete[] dmatrix;
		}

		return ARM_OK;;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (dmatrix)
		{
		for (il=0;il<size;il++)
		{		
			delete[] dmatrix[il];
		}

		delete[] dmatrix;
		}

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: CorrMatrix : unrecognized failure");
		return ARM_KO;
	}

}

//	----------------------------------------------------------------------------------------------------
double  ICMLOCAL_NPV(long pricerId, qCMPMETH measure)
{
	ICM_Pricer* pricer=NULL;
	LocalPersistent::get().convert(pricerId,pricer); 
	return pricer->Price(measure); 	
}


long ICMLOCAL_GenSchedule (double	EffectiveDate,
							double EndDate,
							double ReferenceDate,
							int	Frequency,
							int	DayCount,
							CCString Currency,
							long datesType,
							int modfoll,
							int gapcredit,
							ARM_result& result)
{

	ICM_Leg* leg = NULL;
	ICM_Cds* cds = NULL;
	// ARM_Currency* ccy = NULL;
	
	int stubrule = 3;// stub rule is shortend

	int type = (qTYPEGENODATES) datesType;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];

	CCString msg ("");

	try
	{
	
		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate; 
		if (ReferenceDate == -1.0)
		{	
			stubrule = 1; // stub rule is shortstart
		}
		else
		Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		
		cds = new ICM_Cds((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
							  0,
							  (ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							   0.01,
							   100,100,// 0,0, // notionals 
							   Frequency,
							   DayCount,
							   qCONTINUE_TO_MATURITY,
							   CCSTringToSTLString(Currency), 
							   stubrule,
								DEFAULT_CREDIT_LAG, // const double& CreditLag /*=  DEFAULT_CREDIT_LAG */,
								DEFAULT_FRQ_DEFLEG,// const int& FrequencyDefLeg /* DEFAULT_FRQ_DEFLEG */,
								K_ADJUSTED ,// const int& intRule /* K_ADJUSTED */,
								INCLUDE_MATURITY,// const bool& includematurity /* INCLUDE_MATURITY*/ ,
								K_ADJUSTED,// const int& adjStartDate /* K_ADJUSTED*/ ,
								std::string(), // const std::string& /* char* payCalName = NULL */ ,
								qRunning_Leg,// const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/ ,
								qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/ ,
								ISSUER_UNDEFINE , // const string& name /* = ISSUER_UNDEFINE */ ,
								CREDIT_DEFAULT_VALUE // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ 
							   );

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		// if (pRefDate)
		//	delete [] pRefDate;
		//pRefDate = NULL;

		if (cds == NULL)
		{
			result.setMsg ("ARM_ERR: CDS is null");
			return ARM_KO;
		}

		leg = (ICM_Leg*)cds->GetFeeLeg();

		switch(type)
		{
		case qACC_START_DATE :
			leg->GenScheduleDates(qACC_START_DATE,K_YES,modfoll,gapcredit,"123");
			break;
		case qACC_END_DATE :
			leg->GenScheduleDates(qACC_END_DATE,K_YES,modfoll,gapcredit,"123");
			break;
		case qPAY_DATE :
			leg->GenScheduleDates(qPAY_DATE,K_YES,modfoll,gapcredit,"123");
			break;
		case qOBS_START_DATE :
			leg->GenScheduleDates(qOBS_START_DATE,K_YES,modfoll,gapcredit,"123");
			break;
		case qOBS_END_DATE :
			leg->GenScheduleDates(qOBS_END_DATE,K_YES,modfoll,gapcredit,"123");
			break;
		}

		if (cds)
			delete cds;
		cds = NULL;

	}

	catch(Exception& x)
    {
		x.DebugPrint();

		//if (ccy)
		//	delete ccy;
		//ccy = NULL;

		if (cds)
			delete cds;
		cds = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		// if (pRefDate)
		//	delete [] pRefDate;
		// pRefDate = NULL;

		ARM_RESULT();
	}

	return ARM_OK;
}


/**
extern long ICMLOCAL_ZCPrice(double ValoSpread,
							  double Maturity,
							  double Recovery,
							  double NotionalAmount,
							  ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(ICM_ZCPrice(ValoSpread,Maturity,Recovery,NotionalAmount));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCPrice : unrecognized failure");
		return ARM_KO;
	}
};

extern long ICMLOCAL_ZCFDPrice(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(ICM_ZCFDPrice(ValoSpreads,Maturity,Beta,Recovery,NotionalAmount));
	//	result.setDouble(ICM_ZCFDPrice(ValoSpreads,Shifts,Maturity,Beta,Recovery,NotionalAmount));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCFDPrice : unrecognized failure");
		return ARM_KO;
	}
};

extern long ICMLOCAL_ZCNDPrice(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(ICM_ZCNDPrice(ValoSpreads,Maturity,Beta,Counter,Recovery,NotionalAmount));
	//	result.setDouble(ICM_ZCNDPrice(ValoSpreads,Shifts,Maturity,Beta,Counter,Recovery,NotionalAmount));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCNDPrice : unrecognized failure");
		return ARM_KO;
	}
};

extern long ICMLOCAL_ZCNDDelta(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   int Counter,
							   int Underl,
							   double Recovery,
							   ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(ICM_ZCNDDelta(ValoSpreads,Maturity,Beta,Counter,Underl,Recovery));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCNDPrice : unrecognized failure");
		return ARM_KO;
	}
};

extern long ICMLOCAL_ZCHdgNDGamma(VECTOR<double>& ValoSpreads,
							   VECTOR<double>& SpreadsAnnVol, 
							   double SpreadBeta,
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(NotionalAmount*ICM_ZCNDwithHedgeGammaExp(ValoSpreads,SpreadsAnnVol,SpreadBeta,Maturity,Beta,Counter,Recovery));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCHdgNDGamma : unrecognized failure");
		return ARM_KO;
	}
};

extern long ICMLOCAL_ZCNDTimeShift(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result)
{
	CCString msg ("");

	try
	{
		result.setDouble(NotionalAmount*ICM_ZCNDTimeshift(ValoSpreads,Maturity,Beta,Counter,Recovery));
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch(...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_ZCNDTimeShift : unrecognized failure");
		return ARM_KO;
	}
};
**/ 

long ICMLOCAL_SetCorrelationMatrix (long ModelMultiCurvesId,
									long CorrMatrixId,
									ARM_result& result)
{

	ICM_CorrMatrix* Matrix = NULL;
	ICM_ModelMultiCurves* MMC = NULL;


	CCString msg ("");
	CCString TmpStr;

	try
	{

		MMC = (ICM_ModelMultiCurves *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelMultiCurvesId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(MMC, ICM_MODELMULTICURVES) == 0) 
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		Matrix = (ICM_CorrMatrix *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrMatrixId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Matrix, ICM_CORRMATRIX) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Matrix is not of a good type");
			return ARM_KO;
		}

		MMC->SetCorrelation(Matrix);

		return ARM_OK;;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: SetCorrelationMatrix : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_CloneCorrMatrixBary(int CorrMatrixId,
								  double Beta,
								  int UpOrDown,	
								  ARM_result& result,
								  long objId)
{
	long mtxId;

	ICM_CorrMatrix* pmtx = NULL;
	ICM_CorrMatrix* pnewmtx = NULL;
	ICM_CorrMatrix* pmtxold = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pmtx = (ICM_CorrMatrix *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrMatrixId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pmtx, ICM_CORRMATRIX) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Matrix is not of a good type");
			return ARM_KO;
		}

		pnewmtx = (ICM_CorrMatrix*) pmtx->Clone();
		pnewmtx->ModifyCorrMatrixForBary(Beta,UpOrDown);

		if (pnewmtx == NULL)
		{
			result.setMsg ("ARM_ERR: Correlation Matrix is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			mtxId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pnewmtx);

			if (mtxId == RET_KO)
			{
				if (pnewmtx)
					delete pnewmtx;
				pnewmtx = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mtxId);

			return ARM_OK;
		}
		else
		{
			pmtxold = (ICM_CorrMatrix *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pmtxold, ICM_CORRMATRIX) == 1)
			{
				if (pmtxold)
				{
					delete pmtxold;
					pmtxold = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pnewmtx, objId);

				return ARM_OK;
			}
			else
			{
				if (pnewmtx)
					delete pnewmtx;
				pnewmtx = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (pnewmtx)
			delete pnewmtx;
		pnewmtx = NULL;

		ARM_RESULT();
	}
}

extern long ICMLOCAL_BSGreeks (long pricerId,
							   long greektype,
							   ARM_result& result)
{
	ICM_Pricer_Option* pricer=NULL;
	double greek=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer_Option *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		int greektype2 = (int)greektype ;

		greek = pricer->ComputeBSGreeks(greektype2);

		result.setDouble(greek);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Sensitivity : unrecognized failure");
		return ARM_KO;
	}

}

extern long ICMLOCAL_Greeks (long pricerId,
							   long greektype,
							   ARM_result& result)
{
	ICM_Pricer_Tree_Binomial_SpreadOption* pricer=NULL;
	double greek=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer_Tree_Binomial_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		int greektype2 = (int)greektype ;

		greek = pricer->ComputeGreeks(greektype2);

		result.setDouble(greek);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Sensitivity : unrecognized failure");
		return ARM_KO;
	}

}

extern long ICMLOCAL_ImpliedVol(long pricerId,
							    double Price,
								ARM_result& result)
{
	
	ICM_Pricer* pricer=NULL;
	double Vol=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		Vol = pricer->ComputeImpliedVol(Price);

		result.setDouble(Vol);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Arm_Credit_ImpliedVol : unrecognized failure");
		return ARM_KO;
	}


}

extern long ICMLOCAL_FwdSpreadPricer(long pricerId,
									 const ARM_Date& date1, // double matu1,
									 const ARM_Date& date2, // double matu2,
									 ARM_result& result)
{
	ICM_Pricer_Cds* pricer=NULL;
	double fwdspread=0.0;



	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// ARM_Date NewAsOfDate;
		
		pricer = dynamic_cast<ICM_Pricer_Cds *> (LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId));

		if(!pricer){
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}
		double dur = 0.;
		fwdspread = pricer->Compute_Fwd_Spread(date1, date2, dur);
		result.setDouble(fwdspread);	
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
	

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Compute Forward Spread : unrecognized failure");
		return ARM_KO;
	}

}

double ICMLOCAL_VirtualCdsSpread (long pricerId, double MaturityDate, ARM_result& result)
{
	ICM_Pricer_CDSIndex* pricer=NULL;
	double Spread=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* pMaturityDate=new char[11];

	CCString msg ("");



	try
	{
		pricer = (ICM_Pricer_CDSIndex*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);
		
		Local_XLDATE2ARMDATE(MaturityDate,pMaturityDate);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		Spread = pricer->GetVirtualCDSSpread((ARM_Date) pMaturityDate);
		result.setDouble(Spread);

		if (pMaturityDate)
			delete [] pMaturityDate;
		pMaturityDate = NULL;

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: VirtualCdsSpread : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_CORRELATION_STRIKE(const ARM_Date& AsOf,
								 const std::string& Name,
								 const std::vector<std::string>& labels,
								const std::vector<long>& volcurves, 
								const std::vector<double>& proportions,
								const std::vector<double>& smilestrikelow,
								const std::vector<double>& smilestrikehight,
								const ICM_QMatrix<double>&  fullStrikeLow,
								const ICM_QMatrix<double>&  fullStrikeUp,
								const std::vector<long>& IndexIds,
								ARM_result& result, 
								long objId)
{

	// ARM_Date AsOfModified;AsOfModified.Today();
	// char* pAsOfDate=NULL;

	ICM_Smile_Correlation* PrevMatrix = NULL;
	ICM_Smile_Correlation* Matrix = NULL;
	
	// char** labels2 = NULL;
	std::vector<std::string> labels2; 
	int size = 0;
	int il=0,il2=0;

	vector<const ARM_VolCurve*> VVolcurves;
	ARM_Vector Vproportions;
	ARM_Vector* Vsmilestrikelow=NULL;
	ARM_Vector* Vsmilestrikehight=NULL;
	vector<const ICM_Credit_Index*> VVIndex;

	long CorrId = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	CCString TmpStr;

	try
	{


		size = labels.size();

		// labels2 = new char*[size];
		labels2.resize(size); 
		
		VVolcurves.resize(size);
		Vproportions.Resize(size);
		if (smilestrikelow.size()>0) Vsmilestrikelow=new ARM_Vector(size,0.);
		if (smilestrikehight.size()>0) Vsmilestrikehight=new ARM_Vector(size,0.);

	  if (IndexIds.size()>0) 
		{	if (IndexIds[0] != -1) 
				VVIndex.resize(size);  
		}

		for (il=0; il<size; il++)
		{
			TmpStr = labels[il].c_str();
			TmpStr.trim_right();
			labels2[il] = TmpStr.GetStr();	

			ARM_VolCurve* Volcurves = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(volcurves[il]);

			VVolcurves[il] = Volcurves;
			Vproportions.InitElt(il,proportions[il]);
			if (smilestrikelow.size()>0) Vsmilestrikelow->InitElt(il,smilestrikelow[il]);
			if (smilestrikehight.size()>0) Vsmilestrikehight->InitElt(il,smilestrikehight[il]);

			if (IndexIds.size()>0)
			{
				if (IndexIds[il] != -1)
				{
				ICM_Credit_Index* pIndexIds = (ICM_Credit_Index*) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexIds[il]);
				VVIndex[il] = pIndexIds;
				}
			}	
		}


		if (IndexIds.size()>0 && IndexIds[0] != -1)  
		{
			Matrix = new ICM_Smile_Correlation(AsOf,Name, &(*VVolcurves.begin()), labels2,&Vproportions,&(*VVIndex.begin()));
		}
		else
		{
			if (fullStrikeLow.IsEmpty() || fullStrikeUp.IsEmpty() ) 

// FIXMEFRED: mig.vc8 (30/05/2007 18:04:14):cast
Matrix = new ICM_Smile_Correlation(AsOf,Name, &(*VVolcurves.begin()), labels2,&Vproportions,Vsmilestrikelow,Vsmilestrikehight);
			else 
			{
				Matrix = new ICM_Smile_Correlation ; 
				Matrix->Set(AsOf,Name,&(* VVolcurves.begin()), labels2,Vproportions,fullStrikeLow,fullStrikeUp);
			}
		}

		/** if (labels2)
		{
			for (il=0;il<size;il++) {
				delete labels2[il];
				//delete VVolcurves[il];
				//delete VVIndex[il];
			}
		}
		**/ 
		//if (VVIndex) delete[] VVIndex;
		//if (VVolcurves) delete[] VVolcurves;
		//if (labels2) delete[] labels2;

		// if (Vproportions) delete Vproportions;
		if (Vsmilestrikelow) delete Vsmilestrikelow;
		if (Vsmilestrikehight) delete Vsmilestrikehight;

		if (Matrix == NULL)
		{
			result.setMsg ("ARM_ERR: Correlation Object is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix);

			if (CorrId == RET_KO)
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			PrevMatrix = (ICM_Smile_Correlation*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevMatrix, ICM_CORRELATION) == 1)
			{
				if (PrevMatrix)
				{
					delete PrevMatrix;
					PrevMatrix = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix, objId);

				return ARM_OK;
			}
			else
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		
		// if (Vproportions) delete Vproportions;
		if (Vsmilestrikelow) delete Vsmilestrikelow;
		if (Vsmilestrikehight) delete Vsmilestrikehight;

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: CorrMatrix : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_CORRELATION_SMILE_STRIKE(const double& AsOfDate,
									   const string& Name,
									   const vector<string>& labels,
									   const vector<long>& volcurves, 
									   const ARM_Vector& proportions,
									   const ARM_Vector& smilestrikelow,
									   const ARM_Vector& smilestrikehight,
									   const ICM_QMatrix<double>&  fullStrikeLow,
									   const ICM_QMatrix<double>&  fullStrikeUp,
									   const vector<long>& IndexIds,
									   long objId)
{

	ARM_Date AsOfModified;AsOfModified.Today();
	ICM_Smile_Correlation* Matrix = NULL;	
	int size = 0;
	int il=0;

	// char** labels2 = NULL;
	std::vector<std::string> labels2 ;
	vector<const ARM_VolCurve*> VVolcurves;
	vector<const ICM_Credit_Index*> VVIndex;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
	//	result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	try
	{
		if (AsOfDate != -1.0)
		{
			auto_ptr<char> pAsOfDate(new char[11]);
			Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate.get());
			AsOfModified = (ARM_Date) pAsOfDate.get();
		}
		size = labels.size();
		// labels2 = new char*[size];
		labels2.resize(size); 
		VVolcurves.resize(size);
		if (IndexIds.size()>0) 
		{	if (IndexIds[0] != -1) VVIndex.resize(size);  }

		for (il=0; il<size; il++)
		{
			labels2[il] = labels[il] ;
			ARM_VolCurve* Volcurves = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(volcurves[il]);
			VVolcurves[il] = Volcurves;
			if (IndexIds.size()>0)
			{
				if (IndexIds[il] != -1)
				{
				ICM_Credit_Index* pIndexIds = (ICM_Credit_Index*) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexIds[il]);
				VVIndex[il] = pIndexIds;
				}
			}	
		}

		if (IndexIds.size()>0 && IndexIds[0] != -1)  
		{
			Matrix = new ICM_Smile_Correlation(AsOfModified,Name, &(*VVolcurves.begin()), labels2,&proportions,&(*VVIndex.begin()));
		}
		else
		{
			if (fullStrikeLow.IsEmpty() || fullStrikeUp.IsEmpty() ) 
				Matrix = new ICM_Smile_Correlation(AsOfModified,Name, &(*VVolcurves.begin()), labels2,&proportions,&smilestrikelow,&smilestrikehight);
			else 
			{
				Matrix = new ICM_Smile_Correlation ; 				
				Matrix->Set(AsOfModified,Name, &(*VVolcurves.begin()), labels2,proportions,fullStrikeLow,fullStrikeUp);
			}
		}

		// if (labels2) delete[] labels2;

		if (Matrix == NULL)
		{	//result.setMsg ("ARM_ERR: Correlation Object is null");
			return ARM_KO;
		}

		long QId = LocalPersistent::get().adopt(Matrix,objId);
		return QId; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		return ARM_KO;
	}
	catch (...)
	{
		return ARM_KO;
	}

}

long ICMLOCAL_BETA_CORRELATION(const ARM_Date &AsOf,
						  const string& Name,
						  vector<string>& labels,
						  VECTOR<double>& betas,
						  long idIndex1,
						  long idIndex2,
						  ARM_result& result, 
						  long objId)
{

	ICM_Beta_Correlation* PrevMatrix = NULL;
	ICM_Beta_Correlation* Matrix = NULL;
	
	// ARM_Date AsOfModified;AsOfModified.Today();
	// char* pAsOfDate=NULL;

	std::vector<std::string> labels2 ; 
	int size = 0;
	int il=0,il2=0;

	ARM_Vector* Vbetas=NULL;

	long CorrId = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	CCString TmpStr;

	try
	{

		ARM_IRIndex *index1 ; LocalPersistent::get().convert(idIndex1,index1); 
		ARM_IRIndex *index2 ; LocalPersistent::get().convert(idIndex2,index2); 
		size = labels.size();

		// labels2 = new char*[size];
		labels2.resize(size); 
		
		ARM_Vector Vbetas(size,0.);

		for (il=0; il<size; il++)
		{
			TmpStr = labels[il].c_str();
			TmpStr.trim_right();
			labels2[il] = TmpStr.GetStr();	

			Vbetas.InitElt(il,betas[il]);
		}

		if (labels.size()==1)
			Matrix = new ICM_Beta_Correlation(AsOf,betas[0],Name,index1,index2);
		else
			Matrix = new ICM_Beta_Correlation(AsOf,Name,Vbetas,labels2,index1,index2);


		if (Matrix == NULL)
		{
			result.setMsg ("ARM_ERR: Correlation Object is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix);

			if (CorrId == RET_KO)
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			PrevMatrix = (ICM_Beta_Correlation*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevMatrix, ICM_CORRELATION) == 1)
			{
				if (PrevMatrix)
				{
					delete PrevMatrix;
					PrevMatrix = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Matrix, objId);

				return ARM_OK;
			}
			else
			{
				if (Matrix)
					delete Matrix;
				Matrix = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		// for (il=0;il<size;il++)
		// {			
		// 	delete labels2[il];
		// }
		
		if (Vbetas) delete Vbetas;

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: CorrMatrix : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_SetCorrelation (long ModelMultiCurvesId,
							  long CorrelationId,
							  ARM_result& result)
{
	ICM_Correlation* Correlation = NULL;
	ICM_ModelMultiCurves* MMC = NULL;

	CCString msg ("");
	CCString TmpStr;

	try
	{

		MMC = (ICM_ModelMultiCurves *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelMultiCurvesId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(MMC, ICM_MODELMULTICURVES) == 0) 
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		Correlation = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelationId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Correlation, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation is not of a good type");
			return ARM_KO;
		}

		MMC->SetCorrelation(Correlation);

		return ARM_OK;;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: SetCorrelationMatrix : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_GetEqStrikeDown (long objId,
							   CCString indexname,
							  ARM_result& result)
{
	
	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		// if we are working with a model : retrieve the correlation from the pricer. 
		ICM_Smile_Correlation* correl=0;
		ICM_Pricer* pricer; 
		ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId)) ;
		if (model) 
			correl = dynamic_cast<ICM_Smile_Correlation *>( model->GetCorrelation() ) ;
		else if 
			( correl = dynamic_cast<ICM_Smile_Correlation *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  )  ;
		else if ( pricer = dynamic_cast<ICM_Pricer *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  ) 
		{
			pricer->Price(qCMPPRICE); 
			correl = dynamic_cast<ICM_Smile_Correlation*> ( dynamic_cast<ICM_ModelMultiCurves&>(*pricer->GetModel()).GetCorrelation()); 
		}
			
		if (!correl) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_GetEqStrikeDown: No ICM_Smile_Correlation found"); 
		// if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correl, ICM_CORRELATION) == 0) 
		// {
		// 	result.setMsg ("ARM_ERR: Correlation is not of a good type");
		// 	return ARM_KO;
		// }

		Result = correl->GetEqStrikeDown((char*)indexname);

		result.setDouble(Result);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetEqStrikeDown : unrecognized failure");
		return ARM_KO;
	}
}

long ICMLOCAL_GetEqStrikeUp (long objId,
							 CCString indexname,
							  ARM_result& result)
{
	ICM_Correlation* correl=NULL;
	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ICM_Smile_Correlation* correl=0;
		ICM_Pricer* pricer; 
		ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId)) ;
		if (model) 
			correl = dynamic_cast<ICM_Smile_Correlation *>( model->GetCorrelation() ) ;
		else if 
			( correl = dynamic_cast<ICM_Smile_Correlation *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  )  ;
		else if ( pricer = dynamic_cast<ICM_Pricer *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  ) 
		{
			pricer->Price(qCMPPRICE); 
			correl = dynamic_cast<ICM_Smile_Correlation*> ( dynamic_cast<ICM_ModelMultiCurves&>(*pricer->GetModel()).GetCorrelation()); 
		}
		if (!correl) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_GetEqStrikeDown: No ICM_Smile_Correlation found"); 


		Result = correl->GetEqStrikeUp((char*)indexname);

		result.setDouble(Result);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetEqStrikeUp : unrecognized failure");
		return ARM_KO;
	}
}

long ICMLOCAL_GetCorrelStrikeDown (long CorrelId,
						      double yf_maturity,
							  ARM_result& result)
{
	ICM_Correlation* correl=NULL;
	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		correl = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correl, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation is not of a good type");
			return ARM_KO;
		}

		Result = correl->GetCorrelStrikeDown(yf_maturity);

		result.setDouble(Result);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetCorrelStrikeDown : unrecognized failure");
		return ARM_KO;
	}
}

long ICMLOCAL_GetCorrelStrikeUp (long CorrelId,
								 double yf_maturity,
								ARM_result& result)
{

	ICM_Correlation* correl=NULL;
	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		correl = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correl, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation is not of a good type");
			return ARM_KO;
		}

		Result = correl->GetCorrelStrikeUp(yf_maturity);

		result.setDouble(Result);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetCorrelStrikeUp : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_GetCorrelation (long ModelId,
							  ARM_result& result,
							  long objId)
{

	ICM_ModelMultiCurves* model=NULL;
	ICM_Correlation* correlation = NULL;
	ICM_Correlation* PREVcorrelation = NULL;

	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	long CorrId = -1;

	try
	{

		model = (ICM_ModelMultiCurves *) LOCAL_PERSISTENT_OBJECTS->GetObject(ModelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(model, ARM_MODEL) == 0) 
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		correlation = (ICM_Correlation*) model->GetCorrelation()->Clone();

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)correlation);

			if (CorrId == RET_KO)
			{
				if (correlation)
					delete correlation;
				correlation = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			PREVcorrelation = (ICM_Correlation*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PREVcorrelation, ICM_CORRELATION) == 1)
			{
				if (PREVcorrelation)
				{
					delete PREVcorrelation;
					PREVcorrelation = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)correlation, objId);

				return ARM_OK;
			}
			else
			{
				if (correlation)
					delete correlation;
				correlation = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

		result.setLong(CorrId);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetCorrelation : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_GetExpectedLoss (long PricerId,
							   double yearterm,
							   ARM_result& result)
{

	ICM_Pricer_Security* pricer = NULL;
	
	double ExpectedLoss =0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(PricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}
		vector<double> losses;

		ExpectedLoss = pricer->ExpectedLossTranche(yearterm,losses);

		result.setDouble(ExpectedLoss);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetExpectedLoss : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_FwdSpreadAsIndex (const long& defCurve,
								const double& StartDate,
								const double& EndDate,	
								ARM_result& result)
{
	ICM_DefaultCurve* pDefCurve = NULL;
	double FwdSpread = 0.;
	double FlatDur =0.;
	double DurDiff=0.;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	CCString msg ("");
	char* startDate=new char[11];
	char* endDate=new char[11];
	pDefCurve = dynamic_cast<ICM_DefaultCurve*>( LOCAL_PERSISTENT_OBJECTS->GetObject(defCurve) );
	if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pDefCurve, ICM_DEFAULTCURVE) == 0)
	{
		result.setMsg ("ARM_ERR: Model is not of a good type");
		return ARM_KO;
	}
	Local_XLDATE2ARMDATE(StartDate,startDate);
	Local_XLDATE2ARMDATE(EndDate,endDate);
	ARM_Date Start = (ARM_Date)startDate;
	ARM_Date End = (ARM_Date)endDate;
	if (startDate) delete [] startDate;
	if (endDate)   delete [] endDate;

	FwdSpread = pDefCurve->FwdSpread_AsIndex(Start, End, FlatDur, DurDiff);
	result.setDouble(FwdSpread);
	return ARM_OK;
}


long ICMLOCAL_SetProportionsInfos (const long& CorrelId,
								   const CCString& IndexName,
								   const double& proportion,
								   const double& forcedstrikelow,
								   const double& forcedstrikehigh,
								   ARM_result& result)
{

	ICM_Smile_Correlation* correlation = NULL;

	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		correlation = (ICM_Smile_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlation, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation is not of a good type");
			return ARM_KO;
		}

		correlation->SetProportionsInfos(CCSTringToSTLString(IndexName),
										 proportion,
										 forcedstrikelow,
										 forcedstrikehigh);

		result.setLong(0);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetProportionsInfos : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_ComputeImplicitCurveForCDO2(const long& PricerId, 
								   const CCString& Name,
								   const CCString& Tenor,
								   ARM_result& result)
{

	ICM_Pricer_Analytic_Cdo2_Smile* pricer = NULL;

	CCString msg ("");

	try
	{
		pricer = (ICM_Pricer_Analytic_Cdo2_Smile*) LOCAL_PERSISTENT_OBJECTS->GetObject(PricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pricer, ICM_PRICER_ANALYTIC_CDO2_STRIKE) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		double spread = pricer->ComputeImplicitSpread((string)(const char*)Name,(string)(const char*)Tenor);

		result.setDouble(spread);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetProportionsInfos : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_Credit_AddPeriod(double AsOfDate , 
									  CCString Maturity, 
									  CCString Currency , 
									  bool Adj, 
									  qCDS_ADJ AdjCDS, 
									  ARM_result& result)


{

	char* pAsOfDate = new char[11] ;
	Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate);
	
	char buf[30];
	double dDate;

	// ARM_Currency* ccy =NULL;

	ARM_Date dMatu; 




	try
	{	dMatu = AddPeriod((ARM_Date) pAsOfDate, 
							CCSTringToSTLString(Maturity),
							CCSTringToSTLString(Currency),
							Adj,
							AdjCDS);	
		
		dMatu.JulianToStrDate(buf);
		dDate = Local_ARMDATE2XLDATE(buf);

		result.setDouble(dDate);

		if (pAsOfDate) 
			delete[] pAsOfDate;
		pAsOfDate=NULL; 



		return ARM_OK; 

	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Credit_AddPeriod : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_GetBaseCorrelFromSummit(const ARM_Date&  AsOfDate , 
									  CCString Index, 
									  CCString CurveType, 
									  VECTOR<CCString> Currency, 
									  VECTOR<CCString> CvIssuerName, 
									  VECTOR<double> Proportions,
									  VECTOR<double> smilestrikelow,
									  VECTOR<double> smilestrikehight,
									  CCString CorrelName,
									  // long indexId,
									  ARM_result& result,
									  int ObjectId)


{


	// VECTOR<CCString> Names;Names.resize(Currency.size());
	std::vector<std::string> Names ; Names.resize(Currency.size());
	VECTOR<long> volid;volid.resize(Currency.size());
	VECTOR<long> crvid;crvid.resize(Currency.size());
	VECTOR<long> indexid;

	if (smilestrikelow.size()==0) {indexid.resize(Currency.size());}

	ARM_result Aresult_vol;
	long result_vol;
	
	try
	{

	for (int i=0;i<Currency.size();i++)
	{
	
		// double xlDate ;
		result_vol = ARMLOCAL_GetVolFromSummit(Index,
										 Currency[i],
										 CurveType,
										 JulianToXLDate(AsOfDate.GetJulian()),
									     "IRG",
										 "S:5Y",
									     "IRFWDVOL",
										 -1, // no index
									     Aresult_vol);

		if (result_vol==ARM_KO) 
		{	result.setMsg ("ARM_ERR: ICMLOCAL_GetBaseCorrelFromSummit : pb with vol");
			return ARM_KO;
		}
		volid[i] = Aresult_vol.getLong();

		result_vol = ICMLOCAL_GetDPFromSummit (AsOfDate,
										   CvIssuerName[i],
										   CurveType,
										   -1,	
										   	CvIssuerName[i],
											Aresult_vol);

		if (result_vol==ARM_KO) 
		{	result.setMsg ("ARM_ERR: ICMLOCAL_GetBaseCorrelFromSummit : pb with defprob");
			return ARM_KO;
		}
		crvid[i] = Aresult_vol.getLong();

		int size_index=0;
		Names[i] = (Index + (CCString)"_" + Currency[i]);
			
		if (Currency[i]=="EUR") 
		{size_index=TRX_EUR_NBNAMES;}
		else
		{size_index=TRX_USD_NBNAMES;}

		vector<string> labels;labels.resize(size_index);
		for (int k=0;k<size_index;k++) {labels[k]=CvIssuerName[i];}
		ARM_Vector YT(1); YT.Elt(0) = 5.0;
		ARM_Vector spread(1); spread.Elt(0) = 0.;
		if (smilestrikelow.size()==0)
		{
		result_vol = ICMLOCAL_Index( labels[0], labels,1,4,4,YT,spread,"EUR",0,390000.,-1,K_FOLLOWING,0,0,0,0,(qCDS_ADJ)0,qCredit_Adjust20,
			0,1		// roll day and occurence for CMS: not significant
			); 
		if (result_vol==ARM_KO) 
		{	result.setMsg ("ARM_ERR: ICMLOCAL_GetBaseCorrelFromSummit : pb with index");
			return ARM_KO;
		}
		indexid[i] = result_vol;
		}
	}

	result_vol = ICMLOCAL_CORRELATION_STRIKE(AsOfDate,
											(const char*) CorrelName,
											Names,
											crvid, 
											Proportions,
											smilestrikelow,
											smilestrikehight,
											ICM_QMatrix<double>(),
											ICM_QMatrix<double>(),
											indexid,
											Aresult_vol, 
											ObjectId);

	for (int j=0;j<Currency.size();j++)
	{
	LOCAL_PERSISTENT_OBJECTS->FreeObject(volid[j]);
	LOCAL_PERSISTENT_OBJECTS->FreeObject(crvid[j]);
	LOCAL_PERSISTENT_OBJECTS->FreeObject(indexid[j]);
	}

	if (result_vol==ARM_KO) 
	{	result.setMsg ("ARM_ERR: ICMLOCAL_GetBaseCorrelFromSummit : pb with object correlation");
		return ARM_KO;
	}

	result.setLong(Aresult_vol.getLong());

	return ARM_OK; 

	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetBaseCorrelFromSummit : unrecognized failure");
		return ARM_KO;
	}

}


long ICMLOCAL_INDEX_CORRELATION(const double& AsOfDate,
								const CCString& Name,
								const long& CalMethod,
								const long& CreditIndexId,
								VECTOR<double>& StrikeLow,
								VECTOR<double>& StrikeHigh,
								VECTOR<double>& MktBid,
								VECTOR<double>& MktAsk,
								VECTOR<double>& UpfBid,
								VECTOR<double>& UpfAsk,
								VECTOR<double>& InitialCorrelation,
								double RFLBeta0,
								const long ParamId,
								const long MMCId,
								ARM_result& result, 
								long objId)
{

	ARM_Date AsOfModified;AsOfModified.Today();
	char* pAsOfDate=NULL;

	ICM_Parameters* PrevObject = NULL;
	ICM_Parameters* Object = NULL;
	ICM_Parameters* InParams = NULL;
	
	ARM_Vector* VStrikeLow=NULL;
	ARM_Vector* VStrikeHigh=NULL;
	ARM_Vector* VMktBid=NULL;
	ARM_Vector* VMktAsk=NULL;
	ARM_Vector* VInitialCorrelation=NULL;
	ARM_Vector* VUpfBid=NULL;
	ARM_Vector* VUpfAsk=NULL;
	ARM_Vector VLeverages;
	//ICM_Parameters* Object = NULL;

	ICM_ModelMultiCurves* mmc = NULL;
	ICM_Credit_Index* Index = NULL;

	long CorrId = 0;
	int size = 0,il=0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	CCString TmpStr;

	try
	{

		if (AsOfDate != -1.0)
		{
			pAsOfDate=new char[11];
			Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate);
			AsOfModified = (ARM_Date) pAsOfDate;
			delete[] pAsOfDate;
		}

		size = StrikeLow.size();
		if (size)
		{VStrikeLow = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VStrikeLow->Elt(il)=StrikeLow[il];}}

		size = StrikeHigh.size();
		if (size)
		{VStrikeHigh = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VStrikeHigh->Elt(il)=StrikeHigh[il];}}

		size = MktBid.size();
		if (size)
		{VMktBid = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VMktBid->Elt(il)=MktBid[il];}}

		size = MktAsk.size();
		if (size)
		{VMktAsk = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VMktAsk->Elt(il)=MktAsk[il];}}

		size = UpfBid.size();
		if (size)
		{VUpfBid = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VUpfBid->Elt(il)=UpfBid[il];}}

		size = UpfAsk.size();
		if (size)
		{VUpfAsk = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VUpfAsk->Elt(il)=UpfAsk[il];}}

		size = InitialCorrelation.size();
		if (size)
		{VInitialCorrelation = new ARM_Vector(size,0.);
		for (il=0; il<size; il++) {VInitialCorrelation->Elt(il)=InitialCorrelation[il];}}

		mmc = (ICM_ModelMultiCurves*) LOCAL_PERSISTENT_OBJECTS->GetObject(MMCId);
		Index = (ICM_Credit_Index*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditIndexId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Index, ICM_CREDIT_INDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Credit Index is not of a good type");
			return ARM_KO;
		}
		
		InParams = (ICM_Parameters*) LOCAL_PERSISTENT_OBJECTS->GetObject(ParamId);

		ICM_Index_Correlation ObjectIdx(AsOfModified,(string)(const char*)Name,(const qCAL_INDEX_CORR_TYPE)CalMethod,Index,
										VStrikeLow,VStrikeHigh,VMktBid,VMktAsk,VUpfBid,VUpfAsk,&VLeverages,VInitialCorrelation,mmc); 

		ObjectIdx.SetParams(InParams);
		ObjectIdx.SetRFLBeta0(RFLBeta0);
		ObjectIdx.GlobalCalibrate();
		Object = new ICM_Parameters();
		ARM_Vector* values = new ARM_Vector(ObjectIdx.itsImpliedParameters[0].size(),&(*ObjectIdx.itsImpliedParameters[0].begin()));
		Object->Push(values,"RF_PARAMS");

		if (VStrikeLow) delete VStrikeLow;
		if (VStrikeHigh) delete VStrikeHigh;
		if (VMktBid) delete VMktBid;
		if (VMktAsk) delete VMktAsk;
		if (VUpfBid) delete VUpfBid;
		if (VUpfAsk) delete VUpfAsk;
		if (VInitialCorrelation) delete VInitialCorrelation;

		if (Object == NULL)
		{
			result.setMsg ("ARM_ERR: Correlation Object is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Object);

			if (CorrId == RET_KO)
			{
				if (Object)
					delete Object;
				Object = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			PrevObject = (ICM_Parameters*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PrevObject, ICM_PARAMETERS) == 1)
			{
				if (PrevObject)
				{
					delete PrevObject;
					PrevObject = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Object, objId);

				return ARM_OK;
			}
			else
			{
				if (Object)
					delete Object;
				Object = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (VStrikeLow) delete VStrikeLow;
		if (VStrikeHigh) delete VStrikeHigh;
		if (VMktBid) delete VMktBid;
		if (VMktAsk) delete VMktAsk;
		if (VInitialCorrelation) delete VInitialCorrelation;
		if (VUpfBid) delete VUpfBid;
		if (VUpfAsk) delete VUpfAsk;


		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: IndexCorrelation : unrecognized failure");
		return ARM_KO;
	}
}



long ICMLOCAL_CPT_BASE_CORRELATION(const double& AsOfDate,
							//	const CCString& Name,
								const long& CalMethod,
								const long& CreditIndexId,
								VECTOR<double>& StrikeLow,
								VECTOR<double>& StrikeHigh,
								VECTOR<double>& MktBid,
								VECTOR<double>& MktAsk,
								VECTOR<double>& UpfBid,
								VECTOR<double>& UpfAsk,
								VECTOR<double>& InitialCorrelation,
								VECTOR<double>& Leverages,
								VECTOR<double>& basesCorrelation,
								const long& mmcId,
								const int& intstep,
								const int& startlag,
								const int& creditlag,
								VECTOR<long>& PrevCreditIndexId,
								VECTOR<double>& MatrixPrevBC,
								const double step,
								const qOPTIMIZE_TYPE method,
								const long ParamId,
								ARM_result& result)
{

	ARM_Date AsOfModified;AsOfModified.Today();
	char* pAsOfDate=NULL;

	ICM_Index_Correlation* PrevObject = NULL;
	ICM_Index_Correlation* Object = NULL;
	
	ARM_Vector VStrikeLow;
	ARM_Vector VStrikeHigh;
	ARM_Vector VMktBid;
	ARM_Vector VMktAsk;
	ARM_Vector VInitialCorrelation;
	ARM_Vector VUpfBid;
	ARM_Vector VUpfAsk;
	ARM_Vector VLeverages;

	int row = PrevCreditIndexId.size();
	int col = 0;
	if (row) col = MatrixPrevBC.size()/row;

	ICM_QMatrix<double> MPrevBC(row,col);
	vector<ICM_Credit_Index*> VPrevIndex;

	ICM_Credit_Index* Index = NULL;
	//ICM_Credit_Index* PrevIndex = NULL;

	ICM_ModelMultiCurves* mmc = NULL;

	VECTOR<double> BC;

	long CorrId = 0;
	int size = 0,size2 = 0,il=0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	
	try
	{
		for (int il1=0;il1<row;il1++)
		for (int il2=0;il2<col;il2++)
		{MPrevBC(il1,il2) = MatrixPrevBC[il1+il2*row];}
					
		for (int il3=0;il3<PrevCreditIndexId.size();il3++)
		{VPrevIndex.push_back((ICM_Credit_Index*) LOCAL_PERSISTENT_OBJECTS->GetObject(PrevCreditIndexId[il3]));}

		if (AsOfDate != -1.0)
		{
			pAsOfDate=new char[11];
			Local_XLDATE2ARMDATE(AsOfDate,pAsOfDate);
			AsOfModified = (ARM_Date) pAsOfDate;
			delete[] pAsOfDate;
		}

		size = StrikeLow.size();
		if (size)
		{VStrikeLow.Resize(size);
		for (il=0; il<size; il++) {VStrikeLow.Elt(il)=StrikeLow[il];}}

		size = StrikeHigh.size();
		if (size)
		{VStrikeHigh.Resize(size);
		for (il=0; il<size; il++) {VStrikeHigh.Elt(il)=StrikeHigh[il];}}

		size = MktBid.size();
		if (size)
		{VMktBid.Resize(size);
		for (il=0; il<size; il++) {VMktBid.Elt(il)=MktBid[il];}}

		size = MktAsk.size();
		if (size)
		{VMktAsk.Resize(size);
		for (il=0; il<size; il++) {VMktAsk.Elt(il)=MktAsk[il];}}

		size = UpfBid.size();
		if (size)
		{VUpfBid.Resize(size);
		for (il=0; il<size; il++) {VUpfBid.Elt(il)=UpfBid[il];}}

		size = Leverages.size();
		if (size)
		{VLeverages.Resize(size);
		for (il=0; il<size; il++) {VLeverages.Elt(il)=Leverages[il];}}

		size = UpfAsk.size();
		if (size)
		{VUpfAsk.Resize(size);
		for (il=0; il<size; il++) {VUpfAsk.Elt(il)=UpfAsk[il];}}

		size2 = InitialCorrelation.size();
		if (size2)
		{VInitialCorrelation.Resize(size2);
		for (il=0; il<size2; il++) {VInitialCorrelation.Elt(il)=InitialCorrelation[il];}}

		Index = (ICM_Credit_Index*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditIndexId);
		ICM_Parameters* InParams = (ICM_Parameters*) LOCAL_PERSISTENT_OBJECTS->GetObject(ParamId);

		mmc = (ICM_ModelMultiCurves*) LOCAL_PERSISTENT_OBJECTS->GetObject(mmcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Index, ICM_CREDIT_INDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Credit Index is not of a good type");
			return ARM_KO;
		}


		ICM_Index_Correlation ObjectIdx(AsOfModified,"CORR"/*(string)(const char*)Name*/,(const qCAL_INDEX_CORR_TYPE)CalMethod,Index,
									&VStrikeLow,&VStrikeHigh,&VMktBid,&VMktAsk,&VUpfBid,&VUpfAsk,&VLeverages,&VInitialCorrelation,mmc,
									intstep,startlag,creditlag,VPrevIndex,MPrevBC,step,method); 

		if (InParams)
		{ObjectIdx.SetParams(InParams);}
		ObjectIdx.GlobalCalibrate();

		//basesCorrelation.resize(size);
		basesCorrelation.resize(ObjectIdx.GetImpliedCorrelation()->GetSize());

		for (il=0;il<size;il++) 
		{basesCorrelation[il]=ObjectIdx.GetImpliedCorrelation()->Elt(il);}


/*		FILE* pFile = NULL;
		pFile = fopen("c:\\temp\\ICM_Index_Correlation.txt", "w");
		ICM_Index_Correlation Object(AsOfModified,(string)(const char*)Name,(const qCAL_INDEX_CORR_TYPE)CalMethod,Index,
										&VStrikeLow,&VStrikeHigh,&VMktBid,&VMktAsk,&VUpfBid,&VUpfAsk,&VLeverages,&VInitialCorrelation,mmc,
										intstep,startlag,creditlag,VPrevIndex,MPrevBC,step,method); 
		
		fprintf(pFile, " ICM_Index_Correlation before calibration \n");
		Object.View("", pFile);
		fprintf(pFile, " \nIndex \n");
		Index->View("",pFile);
		fprintf(pFile, " \n mmc \n");
		//mmc->View("",pFile);

		Object.GlobalCalibrate();
		fprintf(pFile, " ICM_Index_Correlation after calibration \n");
		Object.View("", pFile);
		if ( pFile ) fclose (pFile);

		basesCorrelation.resize(size);
		for (il=0;il<size;il++) 
		{basesCorrelation[il]=Object.GetImpliedCorrelation()->Elt(il);}
		Object.GenImplyCorrelation("123");
*/
		return ARM_OK;
	}
	catch (Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	catch (...)
	{
		//result.setMsg ("ARM_ERR: IndexCorrelation : unrecognized failure");
		return ARM_KO;
	}

}


double ICMLOCAL_GetDataFromLabel (long pricerId, 
								  CCString Label,
								  ARM_result& result)
{
	ICM_Pricer* pricer=NULL;
	double value=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{result.setMsg ("ARM_ERR: Pb with accessing objects");
	 return ARM_KO;}

	CCString msg ("");

	try
	{
		pricer = (ICM_Pricer *) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		pricer->GetDataFromLabel((string)(const char*)Label,value);
		result.setDouble(value);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GetDuration : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_GetEqStrike (long objId,
							CCString indexname,
							int isUp,
							std::vector<double>& matu,
							std::vector<double>& strikes,
							ARM_result& result)
{
	double Result=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// if we are working with a model : retrieve the correlation from the pricer. 
		ICM_Smile_Correlation* correl=0;
		ICM_Pricer* pricer; 
		ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId)) ;
		if (model) 
			correl = dynamic_cast<ICM_Smile_Correlation *>( model->GetCorrelation() ) ;
		else if 
			( correl = dynamic_cast<ICM_Smile_Correlation *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  )  ;
		else if ( pricer = dynamic_cast<ICM_Pricer *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId))  ) 
		{
			pricer->Price(qCMPPRICE); 
			correl = dynamic_cast<ICM_Smile_Correlation*> ( dynamic_cast<ICM_ModelMultiCurves&>(*pricer->GetModel()).GetCorrelation()); 
		}
		if (!correl) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_GetEqStrike: No ICM_Smile_Correlation found"); 

		
		if (isUp == 1) 
		{
		correl->GetSmileStrikeUp((char*)indexname,matu,strikes);
		}

		if (isUp == 0) 
		{
		correl->GetSmileStrikeDown((char*)indexname,matu,strikes);
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
		result.setMsg ("ARM_ERR: ICMLOCAL_GetEqStrike : unrecognized failure");
		return ARM_KO;
	}
}


extern long ICMLOCAL_CptLeverageLevels(long pricerId,
									 long ParametersId,
									 long	TriggerCorrelationId,
									 VECTOR<double>		MatrixMultiplesInput,
									 VECTOR<double>		MatrixFlagsInput,
									 VECTOR<double>& VectOfLosses,
									 VECTOR<double>& VectOfMaturitiesInYF,
									 ICM_QMatrix<double>*&	MatrixMultiples,
									 ARM_result& result)
{
	ICM_Pricer_Security* pricer = NULL;
	
	ICM_SuperSeniorLeverage* PrevObject = NULL;
	ICM_SuperSeniorLeverage* Object = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Description = NULL;
	ICM_GenCF* TheCashFlows = NULL;
	ICM_QMatrix<double>*	MatrixOutputs = NULL;

	ICM_Correlation* Correlation = NULL;

	double	Value;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pricer = (ICM_Pricer_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Credit Parameters
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(ParametersId); 

		if (TheCashFlows)	Matrix_Description = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		Object = new ICM_SuperSeniorLeverage();
		
		if (Object == NULL)
		{
			result.setMsg ("ARM_ERR: Super Senior Leverage Object is null");
			return ARM_KO;
		}

		// set attributes
		Object->SetPricer(pricer);
		Object->SetDataParameters(Matrix_Description);
		
		// TriggerCorrelationId,
		Correlation = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(TriggerCorrelationId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Correlation, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation is not of a good type");
			return ARM_KO;
		}

		Object->SetTriggeredCorrelation(Correlation);

		int	NbRows	=	Object->Get_Nb_Flags_Rows();
		int	NbCols	=	Object->Get_Nb_Flags_Cols();
		
		// ---------------------------------------------------------------------------------
		// Matrice de MatrixMultiplesInput	
		// should check the size
		ICM_QMatrix<double>* MyMultiplesFlags = new ICM_QMatrix<double>(NbRows, NbCols, 0.0);
		ICM_QMatrix<double>* MyMultiplesInput = new ICM_QMatrix<double>(NbRows, NbCols, 0.0);
		
		int	i, j;

		for (i=0; i<NbRows; i++)
			for (j=0; j<NbCols; j++)
			{
				MyMultiplesFlags->SetValue(i, j, MatrixFlagsInput[i*NbCols+j]);
				MyMultiplesInput->SetValue(i, j, MatrixMultiplesInput[i*NbCols+j]);
			}

		Object->Set_Multiples_Input(MyMultiplesInput);
		Object->Set_Matrix_Flags(MyMultiplesFlags),

		// computes
		Value	=	Object->ComputeLeverage();

		// get values
		MatrixOutputs		=	Object->GetMatrix_Outputs();
		MatrixMultiples		=	(ICM_QMatrix<double>*) MatrixOutputs->Clone();

		Object->Get_Vect_of_Used_Losses(VectOfLosses);
		Object->Get_Vect_of_Used_MaturitiesInYF(VectOfMaturitiesInYF);


		result.setDouble(Value);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_CptLeverageLevels : unrecognized failure");
		return ARM_KO;
	}
}


long ICMLOCAL_SetMatuLabel (long CurveId, 
								  VECTOR<CCString>& MatuLabels,
								  ARM_result& result)
{
	ARM_VolLInterpol* VolCrv = NULL;
	double value=0.0;
	int j = 0; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{result.setMsg ("ARM_ERR: Pb with accessing objects");
	 return ARM_KO;}

	CCString msg ("");

	try
	{
		VolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(CurveId);

		/*if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Vol Curve is not of a good type");
			return ARM_KO;
		}*/

		for (j=0; j<ARM_NB_TERMS; j++) 
		{
			strcpy(VolCrv->itsYearTermsX[j], "X");
		}

		for (j=0; j<MatuLabels.size(); j++)
		{	
			strcpy(VolCrv->itsYearTermsX[j], (const char *) MatuLabels[j]);
		}

		result.setDouble(value);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetMatuLabel : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_SetRecovCoef (long SecId, 
							const double& RecovCoef,
							ARM_result& result)
{
	ICM_Mez* mezz = NULL;
	double value=0.0;
	int j = 0; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{result.setMsg ("ARM_ERR: Pb with accessing objects");
	 return ARM_KO;}

	CCString msg ("");

	try
	{
		mezz = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(SecId);
		
		//Mise a jour du Recov Coef associé au Collateral
		mezz->GetCollateral()->SetRecovCoef(RecovCoef);

		result.setDouble(value);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetRecovCoef : unrecognized failure");
		return ARM_KO;
	}
}

extern long ICMLOCAL_SetFees(long securityId,
							 long RefvalueId,
							 ARM_result& result)
{
	ARM_Security* sec=NULL;
	ARM_ReferenceValue* ref=NULL;


	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(securityId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0) 
		{	result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;}

		ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(RefvalueId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_REFERENCE_VALUE) == 0) 
		{	result.setMsg ("ARM_ERR: Reference value is not of a good type");
			return ARM_KO;}

		sec->SetFee(ref);
		result.setMsg("ok");

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_SetFees : unrecognized failure");
		return ARM_KO;
	}

}

extern long ICMLOCAL_SetInterpolationType (long VolCurveId,
										   long InterpolType,
										   ARM_result& result)
{
	ARM_VolCurve* volcurve=NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		volcurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolCurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volcurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object VolCurve is not of a good type");
				return ARM_KO;
			}

		volcurve->SetInterpType(InterpolType);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: SetInterpolationType : unrecognized failure");
		return ARM_KO;
	}

}

extern long ICMLOCAL_SECTORIAL_CORRELATION( const ARM_Date& AsOf,
										   const std::string& structName,
										   qTWO_FACTORS_CORRELATION_TYPE		Sectorial_Correlation_Id,
											const vector<string>&	Labels,
											const VECTOR<int>&		Sector_Membership,
											double				Intra_Sector_Correlation,
											double				Inter_Sector_Correlation,
											const VECTOR<double>&		Sector_Betas,
											const VECTOR<double>&		Sector_Lambdas,
											const VECTOR<double>&		Sector_Betas_Down,
											const VECTOR<double>&		Sector_Lambdas_Down,
											ARM_result&			result, 
											long				objId)
{
	ICM_Correlation_Sector*		Prev_Sectorial_Correlation	= NULL;
	ICM_Correlation_Sector*		Sectorial_Correlation		= NULL;
	
	int i;
	int	Nb_Sectors;
	int	Nb_Names;

	long CorrId = 0;

	// char**	The_Labels	=  NULL;
	std::vector<std::string> The_Labels ;


	CCString TmpStr;
	CCString msg ("");
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects - ICMLOCAL_SECTORIAL_CORRELATION");
		return ARM_KO;
	}

	try
	{
		Nb_Names	=	Labels.size();
		// _Labels	=	new char*[Nb_Names];
		The_Labels.resize(Nb_Names); 
		
		for (i=0; i<Nb_Names; i++)
		{
			TmpStr = Labels[i].c_str();
			TmpStr.trim_right();
			The_Labels[i] = TmpStr.GetStr();
		}
			
		Nb_Sectors	=	Sector_Betas.size();

		Sectorial_Correlation	= new ICM_Correlation_Sector(
				AsOf,
				"whatstructname?",
				Sectorial_Correlation_Id,
				The_Labels,
				Nb_Sectors,
				Sector_Membership,
				Intra_Sector_Correlation,
				Inter_Sector_Correlation,
				Sector_Betas,
				Sector_Lambdas,
				Sector_Betas_Down,
				Sector_Lambdas_Down);

		if (Sectorial_Correlation == NULL)
		{
			result.setMsg ("ARM_ERR: Sectorial Correlation Object is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			CorrId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Sectorial_Correlation);

			if (CorrId == RET_KO)
			{
				if (Sectorial_Correlation)
					delete Sectorial_Correlation;
				Sectorial_Correlation = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object - Sectorial Correlation");				
				return ARM_KO;
			}

			result.setLong(CorrId);

			return ARM_OK;
		}
		else
		{
			Prev_Sectorial_Correlation = (ICM_Correlation_Sector*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Prev_Sectorial_Correlation, ICM_CORRELATION) == 1)
			{
				if (Prev_Sectorial_Correlation)
				{
					delete Prev_Sectorial_Correlation;
					Prev_Sectorial_Correlation = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Sectorial_Correlation, objId);

				return ARM_OK;
			}
			else
			{
				if (Sectorial_Correlation)
					delete Sectorial_Correlation;
				Sectorial_Correlation = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type - ICM_Correlation_Sector");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICM_Correlation_Sector : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_GenTsCorrfromBaseCorr (long irCurveId, 
											long defCurveId, 
											long volcurveId, 
											int creditlag,
											ICM_QMatrix<double>*& Correls,
											ARM_result& result)
{
	ARM_ZeroCurve* ircurve=NULL;
	ICM_DefaultCurve* defcurve=NULL;
	ARM_VolCurve* volcurve=NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		ircurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(irCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ircurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object ircurve is not of a good type");
				return ARM_KO;
			}

		defcurve = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(defCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(defcurve, ICM_DEFAULTCURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object default is not of a good type");
				return ARM_KO;
			}

		volcurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(volcurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volcurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object VolCurve is not of a good type");
				return ARM_KO;
			}

		TSBaseCorrelCallib call;

		Correls = call.GenerateTSBaseCorrelation(ircurve,defcurve,*(ICM_VolInterpol*)volcurve,creditlag);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_GenTsCorrfromBaseCorr : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_Get_IR_Curve_Moved_In_Time (long IRCurveId,
												 double MoveDate,
												ARM_result&		result, 
												long			objId)
{
	ARM_ZeroCurve* ircurve=NULL;

	ARM_ZeroCurve* moved_ircurve = NULL;
	ARM_ZeroCurve* Prev_moved_ircurve = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	long Moved_IR_Id= 0;
	
	char* pMoveDate=new char[11];
	
	try
	{

		ircurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(IRCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ircurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object ircurve is not of a good type");
				return ARM_KO;
			}


		Local_XLDATE2ARMDATE(MoveDate, pMoveDate);
		
		moved_ircurve	=	GenerateIRCurveMovedInTime(ircurve, (ARM_Date) pMoveDate);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			Moved_IR_Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)moved_ircurve);

			if (Moved_IR_Id == RET_KO)
			{
				if (moved_ircurve)
					delete moved_ircurve;
				moved_ircurve = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(Moved_IR_Id);

			return ARM_OK;
		}
		else
		{
			Prev_moved_ircurve = (ARM_ZeroCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Prev_moved_ircurve, ARM_ZERO_CURVE) == 1)
			{
				if (Prev_moved_ircurve)
				{
					delete Prev_moved_ircurve;
					Prev_moved_ircurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)moved_ircurve, objId);

				return ARM_OK;
			}
			else
			{
				if (moved_ircurve)
					delete moved_ircurve;
				moved_ircurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Get_IR_Curve_Moved_In_Time : unrecognized failure");
		return ARM_KO;
	}

}
// EOF %M%					   


long ICMLOCAL_Math_Bivariate_normale(double x,double y, double rho, ARM_result& result)
{

	CCString msg ("");

	try
	{

		double price = NAG_bivariate_normal_dist(x,y,rho);

		result.setDouble(price);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Math_Bivariate_normale : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_Math_random_uniform(int seed, ARM_result& result)
{

	CCString msg ("");

	try
	{
		if (seed>=0) {NAG_random_init_repeatable(seed);}

		double price = NAG_random_continuous_uniform();

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Math_Bivariate_normale : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_Math_random_normal(double a,double b,int seed, ARM_result& result)
{

	CCString msg ("");

	try
	{
		if (seed>=0) {NAG_random_init_repeatable(seed);}

		double price = NAG_random_normal(a,b);

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Math_Bivariate_normale : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_QMatrix(const ICM_QMatrix<double>& aQMatrix, long objId)
{
	ICM_QMatrix<double>* perQMatrix = new ICM_QMatrix<double>(aQMatrix);

	if ( !perQMatrix)
		ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_QMatrix can't create QMatrix");
	// insertion in the cache :
	long QId = LocalPersistent::get().adopt(perQMatrix,objId );

	return QId; 
}


long ICMLOCAL_QuickELoss(	   double pdefault,
							   double recovery,
							   int nbnames,
							   double strikedw,
							   double strikeup,
							   double correldw,
							   double correlup,
							   int lossesno,
							   int intstep,
							   bool lhp,
							   ARM_result& result)
{

	CCString msg ("");

	try
	{

		vector<double> losses;

	 double	price = cpt_quick_el_full_homog(pdefault,
							   recovery,
							   nbnames,
							   strikedw,
							   strikeup,
							   correldw,
							   correlup,
							   losses,
							   intstep,
							   lhp);

		if ((lossesno>=0)&&(lossesno<losses.size()))
		{price = losses[lossesno];}
		else if (lossesno>=0)
		{price = -999.;}

		result.setDouble(price);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_QuickELoss : unrecognized failure");
		return ARM_KO;
	}

}

long ICMLOCAL_QuickCDOPV(long discid,
						 long pdefid,
						 double startdate,
						 double enddate,
						 int frequency,
						 double rate,
						 double strikedw,
					     double strikeup,
					     double correldw,
					     double correlup,
					     double notfeeleg,
					     double notdefleg,
					     double recovery,
					     int nbnames,
					     int intstep,
					     bool LHP,
						 int pvtype,
						 long correlid,
						 ARM_result& result)
{
	ARM_VolCurve* pcorrel=NULL;
	ICM_DefaultCurve* pdef=NULL;
	ARM_ZeroCurve* disc=NULL;
	ARM_Date asof;

	char pstartdate[11];
	char penddate[11];
	Local_XLDATE2ARMDATE(startdate,pstartdate);
	Local_XLDATE2ARMDATE(enddate,penddate);

	ARM_Date StartDate(pstartdate);
	ARM_Date EndDate(penddate);
	double feepv;
	double defpv;

	double price=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		double NPV = 0.;

		disc = (ARM_ZeroCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(discid);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(disc, ARM_ZERO_CURVE) == 0) 
		{
			result.setMsg ("ARM_ERR: zerocurve is not of a good type");
			return ARM_KO;
		}

		asof = disc->GetAsOfDate();

		pdef = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(pdefid);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pdef, ICM_DEFAULTCURVE) == 0) 
		{
			result.setMsg ("ARM_ERR: default curve is not of a good type");
			return ARM_KO;
		}

		if (correlid == -1){
		NPV = FastCDOPricing(asof,StartDate,EndDate,frequency,rate,strikedw,strikeup,
					  correldw,correlup,notfeeleg,notdefleg,pdef,disc,recovery,
					  nbnames,intstep,feepv,defpv,LHP);}
		else {

			pcorrel = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(correlid);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pcorrel, ARM_VOL_CURVE) == 0) 
			{
			result.setMsg ("ARM_ERR: correlation is not of a good type");
			return ARM_KO;
			}

			NPV = FastCDOPricing_TSR(asof,StartDate,EndDate,frequency,rate,strikedw,strikeup,
					  notfeeleg,notdefleg,pdef,disc,recovery,nbnames,pcorrel,intstep,feepv,defpv,LHP);
		}

		switch (pvtype)
		{
		case qCMPFEELEGPV:
			result.setDouble(feepv);
			break;
		case qCMPDEFLEGPV:
			result.setDouble(defpv);
			break;
		default :
			result.setDouble(NPV);
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
		result.setMsg ("ARM_ERR: QuickCDOPV : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_Math_Interpol(VECTOR<double>& a,
								   VECTOR<double>& b,
								   double value,
								   int type, 
								   double smooth,
								   VECTOR<double>& weights,
								   int modespline,
								   int withC1condition,
								   double leftSlope,
								   double rightSlope,
								   ARM_result& result)
{
	CCString msg ("");

	try
	{

	 double	price = -999.;

	 if (a.size()==b.size())
	 {	
		if (type==0)
			price=LinearVectorInterpol(a,b,value);
		 else if (type==1)
		 {
			if ((value>=a[0]) && (value<=a[a.size()-1]))
				{price=SplineInterpol(a,b,value,smooth,weights);}
		 }
		 else if (type==2)
		 {
			if ((value>=a[0]) && (value<=a[a.size()-1]))
			{
				vector<double> alpha;alpha.resize(a.size());
				for (int i=0;i<a.size();i++) {alpha[i]=0.5;}

// FIXMEFRED: mig.vc8 (30/05/2007 18:05:05):cast
				CCubicspline c(&(*a.begin()),&(*b.begin()),a.size(),(eEMMODE)modespline,(eSPLINE_EMMODE)withC1condition,leftSlope,rightSlope,alpha);
				price = c(value);
			}
		 }	
		 else if (type==3)
		 { price = HermiteInterpol(a,b,value); }
	 }

	result.setDouble(price);

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_QuickELoss : unrecognized failure");
		return ARM_KO;
	}

}

extern long ICMLOCAL_SCHEDULE_INFO(const double	EffectiveDate,
									const double MaturityDate,			
									const int	payFrequency,
									const int ResetFreq ,
									const int	DayCount,
									const int	Stubrule,
									const int	intRule,
									const std::string& payCalName,
									const int PayTiming,
									const int ResetTiming,
									const int fwdRule,
									const bool	IncludeMaturity,
									const int adj,
									const int	intStartAdj,
									const int AccDayCount,
									const double ReferenceDate, // roll date
									const double FirstCpnEffDate,
									const int CDSAdj,
									long objId = -1)
{
		CCString msg ("");
		char pstartdate[11];
		char penddate[11];
		Local_XLDATE2ARMDATE(EffectiveDate,pstartdate);
		Local_XLDATE2ARMDATE(MaturityDate,penddate);
		ARM_Date StartDate(pstartdate);
		ARM_Date EndDate(penddate);

		ARM_Date * pRollDate = NULL;
		if( ReferenceDate != -1) {
				char cRolldate[11];
				Local_XLDATE2ARMDATE(ReferenceDate,cRolldate);
				pRollDate = new ARM_Date(cRolldate);
		}

		ARM_Date * pFCEDDate = NULL;
		if( ReferenceDate != -1) {
				char cFCEDdate[11];
				Local_XLDATE2ARMDATE(FirstCpnEffDate,cFCEDdate);
				pFCEDDate = new ARM_Date(cFCEDdate);
		}

		ICM_Schedule_Info* pScheduleInfo = NULL;
		pScheduleInfo= new ICM_Schedule_Info(StartDate,
							EndDate,			
							payFrequency,
							ResetFreq ,
							DayCount,
							Stubrule,
							intRule,
							payCalName,
							PayTiming,
							ResetTiming,
							fwdRule,
							IncludeMaturity,
							adj,
							intStartAdj,
							AccDayCount,
							pRollDate,
							pFCEDDate,
							(qCDS_ADJ) CDSAdj);
		if (pFCEDDate) delete pFCEDDate; 
		pFCEDDate = NULL;

		if (pRollDate) delete pRollDate; 
		pRollDate = NULL;

		long QId = LocalPersistent::get().adopt(pScheduleInfo,objId );

		return QId; 
	
}
long ICMLOCAL_FLAT_CORRELATION(const ARM_Date&AsOf,
							   const std::string& structName,
							   double correlValue,
							   long idIndex1,
							   long idIndex2,
							   long prevId)
{
	ARM_IRIndex* idx1 ; LocalPersistent::get().convert(idIndex1,idx1); 
	ARM_IRIndex* idx2 ; LocalPersistent::get().convert(idIndex2,idx2); 
	ICM_FlatCorrel* item = new ICM_FlatCorrel(AsOf,structName,idx1,idx2,correlValue); 
	long newId = LocalPersistent::get().adopt(item,prevId); 
	return newId; 
}


void ICMLOCAL_PriceVector(long pricerId, const std::string& measure, ARM_Vector&output)
{
	
	ICM_Pricer* pricer=NULL;
	

	LocalPersistent::get().convert(pricerId,pricer); 
	if (!pricer) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_PriceVector: Pricer not defined") ;

	const ARM_Vector& pricevector = pricer->PriceVector(measure) ;
	output=pricevector ; 
}


void ICMLOCAL_GenPrice(long pricerId, const std::string& measure, double& price)
{
	
	ICM_Pricer* pricer=NULL;

	LocalPersistent::get().convert(pricerId,pricer); 
	if (!pricer) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_GenPrice: Pricer not defined") ;

	price = pricer->Price(measure) ;
		
}


extern long ICMLOCAL_Math_CF_SpreadOption(double& yearterm,
										double& strike,
										double& correlation,
										VECTOR<double>& coefs,
										VECTOR<double>& spots,
										VECTOR<double>& vol,
										int intstep,
										ARM_result& result)
{
	CCString msg ("");

	try
	{

	double	price = -999.;

	LinearLognProcess LP(yearterm,strike,correlation,coefs,spots,vol,intstep);
	price = LP.IntBND2(-10.,10.);

	result.setDouble(price);

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Math_CF_SpreadOption : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_Register(long address,ARM_result& result)
{
	CCString msg ("");

	try
	{

	Register(address);	

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Register : unrecognized failure");
		return ARM_KO;
	}

}


extern long ICMLOCAL_Calibrator(VECTOR<long>& security_vector,
								VECTOR<double>& price_vectorBid,
								VECTOR<double>& price_vectorAsk,
								VECTOR<CCString>& iparameters_vector,
								VECTOR<double>& tsparams_vector,
								long l_model,
								long ipricerType,
								long l_parameters,
								long l_parameters_inf,
								long l_parameters_sup,
								int iPricingType,
								long l_Optimparameters,
								long prevId)
{
	CCString msg ("");

	int i=0;
	vector<string> parameters_vector;
	vector<ICM_Security*>	SecurityVector;SecurityVector.resize(security_vector.size());
	ARM_Model*				Model;
	ICM_Parameters*			Parameters;
	ICM_Parameters*			Parameters_inf;
	ICM_Parameters*			Parameters_sup;
	ICM_Parameters*			Parameters_out;
	ICM_Parameters*			OptimParameters;

	qCMPMETH				PricingType = (qCMPMETH) iPricingType;
	ARM_CLASS_NAME			pricerType = (ARM_CLASS_NAME) ipricerType;

	for (i=0;i<iparameters_vector.size();i++)
	{
		parameters_vector.push_back((const char*)iparameters_vector[i]);
	}

	for (i=0;i<security_vector.size();i++)
	{
	LocalPersistent::get().convert(security_vector[i],SecurityVector[i]); 
	if (!SecurityVector[i]) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: security not defined") ;
	}

	LocalPersistent::get().convert(l_model,Model); 
	if (!Model) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: Model not defined") ;

	LocalPersistent::get().convert(l_parameters,Parameters); 
	if (!Parameters) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: Parameters not defined") ;

	LocalPersistent::get().convert(l_parameters_inf,Parameters_inf); 
	if (!Parameters_inf) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: Parameters_inf not defined") ;

	LocalPersistent::get().convert(l_parameters_sup,Parameters_sup); 
	if (!Parameters_sup) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: Parameters_sup not defined") ;

	LocalPersistent::get().convert(l_Optimparameters,OptimParameters); 
	if (!OptimParameters) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_Calibrator: Optim Parameters not defined") ;


	ICM_Calibrator cal(SecurityVector,price_vectorBid,price_vectorAsk,parameters_vector,tsparams_vector,
			Model,pricerType,Parameters,Parameters_inf,Parameters_sup,PricingType,OptimParameters);


	Parameters_out = cal.Optimize();

	long newId = LocalPersistent::get().adopt(Parameters_out,prevId); 
	
	return newId; 

}

extern long ICMLOCAL_RandomGenerator(const qRAN_GEN& RandomType,
									 long parameterId, long prevId)
{

	ICM_Parameters* parameter = NULL;
	ICM_RandomGenerator* RandomGenerator = NULL;
	int iIsRandom; 
	bool IsRandom = false; 
	int InitialSeed; 
	ICM_QMatrix<string>* randGenerator;
	ICM_RandomGenerator* RandomGeneratorForLaws = NULL;

	LocalPersistent::get().convert(parameterId,parameter); 
	if (!parameter) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Parameter not defined") ;
	if(RandomType < q_INV_CUM_NORM_ACKLAM ) {	
		iIsRandom = parameter->getLong("IsRandom");
		if (iIsRandom == 0) 
			IsRandom = true; 
		InitialSeed = parameter->getLong("InitialSeed");
	} else { // case laws
		randGenerator = parameter->GetColVectStr("STR_UNIFORM_GENERATOR");
		
		if(!randGenerator)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Parameter RandomGeneratorForLaws not defined") ;
		}
		string strIDRandomGenerator = randGenerator->Elt(0);
		long lUniGenerator = LocalGetNumObjectId (CCString(strIDRandomGenerator.c_str()));
		LocalPersistent::get().convert(lUniGenerator,RandomGeneratorForLaws); 
		if (!RandomGeneratorForLaws)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Parameter RandomGeneratorForLaws not defined") ;
		}
	}
		
	switch(RandomType)
	{	
	case q_NAG  :
	{
		RandomGenerator = new ICM_RandomNag(InitialSeed, IsRandom); 
		break;
	}
	case q_RAN1 :
	{
		RandomGenerator = new ICM_RandomRan1(InitialSeed); 
		break;
	}
	case q_RAN2 :
	{
		RandomGenerator = new ICM_RandomRan2(InitialSeed); 
		break;
	}
	case q_DEF :
	{
		RandomGenerator = new ICM_RandomRanDef(InitialSeed, IsRandom); 
		break;
	}
	case q_RANMAR:
	{
		int InitialSeed1; 
		int InitialSeed2; 
		InitialSeed1 = parameter->getLong("InitialSeed1");
		InitialSeed2 = parameter->getLong("InitialSeed2");
		RandomGenerator = new ICM_RandomRanmar(InitialSeed1, InitialSeed2); 
		break;
	}
	case q_RNG_STR :
	{
		RandomGenerator = new ICM_RandomRNG_Str(InitialSeed); 
		break;
	}
	case q_KISS :
	{
		unsigned long InitialSeed_X = 0; 
		unsigned long InitialSeed_Y = 0; 
		unsigned long InitialSeed_Z = 0;
		unsigned long InitialSeed_C = 0;

		InitialSeed_X = parameter->getLong("X");
		InitialSeed_Y = parameter->getLong("Y");
		InitialSeed_Z = parameter->getLong("Z");
		InitialSeed_C = parameter->getLong("C");
		RandomGenerator = new ICM_RandomKISS(InitialSeed_X, InitialSeed_Y, InitialSeed_Z, InitialSeed_C); 
		break;
	}
	case q_INV_CUM_NORM_ACKLAM :
		RandomGenerator = new ICM_RandomInvNormAcklam(*RandomGeneratorForLaws);
		break;
	case q_INV_CUM_NORM_MORO : 
	{
		RandomGenerator = new ICM_RandomInvNormMoro(*RandomGeneratorForLaws);
		break;
	}
	default :
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: type of random not defined") ;
	}
	// for all generators
	long RandId = LocalPersistent::get().adopt(RandomGenerator, prevId);
	return RandId ; 
}

void ICMLOCAL_GenerateOneRandom(long RandomGenId, double& Random)
{

	ICM_RandomGenerator* RandomGenerator = NULL;

	LocalPersistent::get().convert(RandomGenId,RandomGenerator); 
	if (!RandomGenerator) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Random Generator not defined") ;

	Random = RandomGenerator->GenerateOneRandom();

}
void ICMLOCAL_GenerateRandoms(long RandomGenId, ARM_Vector& RandomVector)
{

	ICM_RandomGenerator* RandomGenerator = NULL;

	LocalPersistent::get().convert(RandomGenId,RandomGenerator); 
	if (!RandomGenerator) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Random Generator not defined") ;

	RandomGenerator->GenerateRandoms(RandomVector);

}
void ICMLOCAL_ResetRandom(long RandomGenId)
{

	ICM_RandomGenerator* RandomGenerator = NULL;

	LocalPersistent::get().convert(RandomGenId,RandomGenerator); 
	if (!RandomGenerator) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RandomGenerator: Random Generator not defined") ;

	RandomGenerator->reset();

}



