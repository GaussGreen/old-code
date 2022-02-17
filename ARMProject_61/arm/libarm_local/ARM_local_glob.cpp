#include "firstToBeIncluded.h"

#include <gpmodels\typedef.h>

#include <glob\armglob.h>
#include <glob\securityflows.h>
#include <inst\security.h>
#include <inst\bond.h>
#include <inst\bondtec.h>
#include <inst\xccyconvert.h>
#include <inst\option.h>
#include <inst\swaption.h>
#include <inst\swap.h>
#include <pricer\ipricer.h>
#include <mod\frmana.h>
#include <inst\powrev.h>
#include <inst\fixleg.h>
#include <inst\swapleg.h>
#include <inst\spreadoption.h>
#include <inst\armdigital.h>
#include <mod\bsflexible.h>
#include <crv\volint.h>
#include <util\interpol.h>
#include <util\fromto.h>
#include <crv\volcube.h>
#include <crv\hypercube.h>
#include <crv\correlmanager.h>
#include <mod\convadjustmanager.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>
#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
#include <ARM\libarm_local\ARM_local_wrapper.h>
#include <ARM\libarm_local\ARM_local_glob.h>


/*
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_nthtd.h>
#include <ICMKernel\inst\icm_option.h>
#include <ICMKernel\inst\icm_cdo2.h>
*/
#include <ARM\libarm_local\undef_va_vars.h>

#include <windows.h>

#ifdef GetObject
#undef GetObject
#endif


/// part of the generic pricer
#include <GP_Base\gpbase\autocleaner.h>
#include <GP_Base\gpbase\gpmatrix.h>
#include <GP_Base\gpbase\singleton.h>
#include <GP_Base\gpbase\stringmanip.h>

#include <GP_Infra\gpinfra\gensecurity.h>
#include <GP_Infra\gpinfra\genpricer.h>
#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\typedef.h>
#include <GP_Infra\gpinfra\modelnrefcall.h>
#include <GP_Infra\gpinfra\gramnode.h>
#include <GP_Infra\gpinfra\pricingadviser.h>
#include <GP_Infra\gpinfra\dealdescription.h>
#include <GP_Infra\gpinfra\mktdatamanagerrep.h>
#include <GP_Infra\gpinfra\mktdatamanager.h>
#include <GP_Infra\gpinfra\curvemodelparam.h>
#include <GP_Infra\gpinfra\argconvdefault.h>

#include <GP_Calib\gpcalib\kerneltogp.h>
#include <GP_Calib\gpcalib\vanillaarg.h>
#include <GP_Calib\gpcalib\calibmethod.h>

#include <GP_Calculators\gpcalculators\pricerfactory.h>
#include <GP_Calculators\gpcalculators\gencalculator.h>
#include <GP_Calculators\gpcalculators\crfcalculator.h>

#include <GP_Calculators\gpcalculators\craspreadcalculator.h>
#include <GP_Calculators\gpcalculators\localcsocalculator.h>
#include <GP_Calculators\gpcalculators\tarncalculator.h>
#include <GP_Calculators\gpcalculators\tarnfxcalculator.h>
#include <GP_Calculators\gpcalculators\maturitycapcalculator.h>
#include <GP_Calculators\gpcalculators\prdccalculator.h>
#include <GP_Calculators\gpcalculators\bermudaswaptioncalculator.h>
#include <GP_Calculators\gpcalculators\captioncalculator.h>
#include <GP_Calculators\gpcalculators\argconvdefault.h>
#include <GP_Calculators\gpcalculators\globalcapcalculator.h>
#include <GP_Calculators\gpcalculators\callablesnowballcalculator.h>

#include <GP_Models\gpmodels\MarketIRModel.h>
#include <GP_Models\gpmodels\Mixture_FX.h>

#include <GP_Inflation\gpinflation\resetmanager.h>

#include <GP_NumLib\gpnumlib\argconvdefault.h>

/// gphelp
#include <GP_Help\gphelp\crmcookies.h>

/// Objects are in namespace, hence the using directive!
using ARM::ARM_GenSecurity;
using ARM::ARM_GenPricer;
using ARM::ARM_GP_Matrix;
using ARM::ARM_GP_Vector;
using ARM::ARM_VanillaArg;
using ARM::ARM_PricingModel;
using ARM::ARM_PricingModelPtr;
using ARM::ARM_GenCalculator;
using ARM::ARM_CRFCalculator;
using ARM::ARM_MaturityCapCalculator;
using ARM::ARM_PRDCCalculator;
using ARM::ARM_PricerFactory;
using ARM::stringGetUpper;
using ARM::ARM_CRASpreadCalculator;
using ARM::ARM_LocalCSOCalculator;
using ARM::ARM_TARNCalculator;
using ARM::ARM_TARNFXCalculator;
using ARM::ARM_BermudaSwaptionCalculator;
using ARM::ARM_CaptionCalculator;
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_MarketIRModel;
using ARM::ARM_ArgConv_CaptionCalibMode;
using ARM::ARM_ModelType;
using ARM::ARM_PricingModelType;
using ARM::ARM_ResetManager;
using ARM::ARM_ArgConv_PricingModelType;
using ARM::ARM_ArgConv_VnsPricingMethod;
using ARM::ARM_CurveModelParam;
using ARM::ARM_GlobalCapCalculator;
using ARM::ARM_CallableSnowBallCalculator;
using ARM::ARM_CRMCookies;
using ARM::ARM_ParamsMixture_Fx;
using ARM::ARM_ArgConv_BaseGenAlgoType;
using ARM::ARM_ArgConv_TransformAlgoType;
using ARM::ARM_ArgConv_YesNo;
using ARM::ARM_ArgConv_TARNFXModelType;
using ARM::ARM_ArgConv_MMCorrelType;

//extern ARMLOCAL_Init* armlocal_init;

/*********************************************
 Déclaration de la liste des objets persistents
**********************************************/
LocalPersistent* LOCAL_PERSISTENT_OBJECTS;

/******************************************/


int DBL_INT(double x)
{
    char   buf[50];
    double v;
    int    res;
  
    v = floor(x);
 
    sprintf(buf, "%lf", v);
 
    sscanf(buf, "%d", &res);
 
    return(res);
}



ARM_Matrix* CreateARMMatrixFromVECTOR(const VECTOR<double>& param, int nbrows, int nbcolumns)
{
    long size = param.size();
    if(size == 0)
        return NULL;

	if( size != nbrows * nbcolumns )
	{
		char msg[255];
		sprintf( msg, "Matrix set with incorrect nb of rows %d and nb of columns %d, total size is %d",
			nbrows, nbcolumns, param.size() );
    	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}

	ARM_Matrix* res = new ARM_Matrix( nbrows, nbcolumns );

	int i,j;
	for ( i=0;i<nbrows; ++i)
		for( j=0; j<nbcolumns; ++j )
			res->Elt(i,j) = param[i*nbcolumns+j];
	return res;
}



ARM_GP_Matrix* CreateARM_GP_MatrixFromVECTOR(const VECTOR<VECTOR<double> >& param)
{
    long nbrows = param.size();
    long nbcolumns;
    if(nbrows == 0)
        return NULL;
    else
        nbcolumns = param[0].size();	

    ARM_GP_Matrix* res = new ARM_GP_Matrix( nbrows, nbcolumns );

	int i,j;
	for ( i=0;i<nbrows; ++i)
		for( j=0; j<nbcolumns; ++j )
			(*res)(i,j) = param[i][j];
	return res;

}


long ARMLOCAL_GetCurrency(long ObjId, ARM_result& result)
{
	ARM_Object* sec = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return(ARM_KO);
	}

	CCString msg("");

	try
	{
		sec = LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY))
		{
			result.setString(((ARM_Security*)sec)->GetCurrencyUnit()->GetCcyName());
			return ARM_OK;
		}
		else
		{
			result.setString((((ARM_AbstractMarketClass*)sec)->GetStrCurrency()).c_str());
			return ARM_OK;
		}
	
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_ARM_Price");
        return(ARM_KO);
	}
}

long ARMLOCAL_ARM_GetDefaultCurrency (ARM_result& result)
{
	result.setString(ARM_DEFAULT_COUNTRY);
	result.setRetCode(ARM_OK);
	return ARM_OK;
}

long ARMLOCAL_ARM_SetDefaultCurrency (const CCString& isoccy, ARM_result& result)
{
	char buf[20];

	CCString msg (" ");
	char* sISO = (char*) isoccy;
	try
	{
		char* curDefCountry = ARM_SetDefaultCountry(sISO);

		if (sISO)
			delete sISO;
		sISO = NULL;

		strcpy(buf,curDefCountry);
	}
	catch(Exception& x)
	{
		if (sISO)
			delete [] sISO;
		sISO = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}

	result.setDouble(0.);
	result.setRetCode(ARM_OK);
	return ARM_OK;
}

long ARMLOCAL_ARM_GetFRFCurrency (ARM_result& result)
{
	result.setString("FRF");
	result.setRetCode(ARM_OK);
	return ARM_OK;
}


double Local_ARMDATE2XLDATE (const CCString& armdate)
{
	int y, m, d;

	if (armdate.GetLen() == 8)
	{
		sscanf (armdate, "%2d/%2d/%2d", &d, &m, &y);
		y += 2000;
	}
	else
	{
		sscanf (armdate, "%2d.%2d.%4d", &d, &m, &y);
	}
	
	return (DAT_struct_to_ssdate (y, m, d));
}


ARM_Date ConvertToARMDATE (double xldate)
{
	char  charDate[11];
	Local_XLDATE2ARMDATE( xldate, charDate );
	return ARM_Date( charDate );
}

long ARMLOCAL_bsflexible (double F, double V,double B, double K,double CallPutFlag,
						  ARM_result& result)
{
	double price;
if (F<0)
		{
			result.setMsg ("ARM_ERR: negative Forward");
			return ARM_KO;
		}
if (V<0)
		{
			result.setMsg ("ARM_ERR: negative Volatility");
			return ARM_KO;
		}
if (K<0)
		{
			result.setMsg ("ARM_ERR: negative Strike");
			return ARM_KO;
		}
		price=bsflexible(F,V,B,K,CallPutFlag);

		result.setDouble(price);

		return ARM_OK;
}


long ARMLOCAL_ARM_Price (long secId, long modId, ARM_result& result)
{
	ARM_Object* sec		= NULL;
   	ARM_Portfolio* pf	= NULL;
	ARM_Object* mod		= NULL;

	double price		= 0.0;

    bool isPf			= false; // A security is a Portfolio?



	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");

		return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		// get the security
		sec = LOCAL_PERSISTENT_OBJECTS->GetObject(secId);
		
        // get the model
		if ( modId != ARM_NULL_OBJECT )
		   mod = LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		// validation: checks whether we are pricing a security that is a portfolio or not
		if ( !dynamic_cast<ARM_GenCalculator*>(sec) 
			&& 
             !dynamic_cast<ARM_GenSecurity*>(sec) 
           )
		{
            if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0 )
		    {
			    result.setMsg("ARM_ERR: ARM Security or Portfolio expected");

			    return(ARM_KO);
		    }

            pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);
			
            if ( pf->GetName() == ARM_PORTFOLIO )
               isPf = true;
		}

        
		// validation on a model
        if ( mod != NULL )
        {
            // Price with a model
            if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0
			     && 
                 !dynamic_cast<ARM_PricingModel*>(mod) 
               )
		    {
			    result.setMsg("ARM_ERR: Expected Model is not of a good type");
			    
                return(ARM_KO);
		    }
		}
		else
		{
			if( !dynamic_cast<ARM_GenCalculator*>(sec) )
		    {
			    result.setMsg("ARM_ERR: only calculators could be priced without a model, input either a calculator or, a security and a model!");
			    
                return(ARM_KO);
		    }
		}

		// call factory
		CC_NS(std,pair)<bool,double> pricingResult = ARM_PricerFactory::Price(sec, mod);

		if (pricingResult.first)
        {
		   price = pricingResult.second;
        }
        else
		{
            if ( ((ARM_Model *) mod)->PricerType() == PT_NONE )
            {
			   /// this is a simple pricer!

			   ARM_IFPricer pricer((ARM_Security*) sec, (ARM_Model *) mod);

			   price = pricer.Price();
            }
            else
            {
                // other cases are the standard ARM cases
			    // Is it a portfolio?
			    
                if (isPf)
			    {
				    long pfSize				= pf->GetSize();
				    ARM_Vector*  pfWeights	= pf->GetWeights();
				    ARM_Security* pfSecurity;

				    for (long i = 0; i < pfSize; i++)
				    {
					    pfSecurity = pf->GetAsset(i);

					    ARM_IFPricer pricer(pfSecurity, (ARM_Model *) mod);

					    price += pricer.Price() * pfWeights->Elt(i);
				    }
			    }
			    else
			    {
				    /// this is also a simple pricer!
				    ARM_IFPricer pricer((ARM_Security*) sec, (ARM_Model*) mod);

				    price = pricer.Price();
			    }
            }
        }

		// save the result!
		result.setDouble(price);

		return ARM_OK;
	}

	// first catch ARM type Exception
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_ARM_Price");
		
        return(ARM_KO);
	}
}


long ARMLOCAL_ARM_Cover (long secId, long modId, ARM_result& result)
{
	ARM_Security* sec=NULL;
	ARM_Model* mod=NULL;
	ARM_Vector*   Cover = NULL;
	long retour = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		if (sec->GetName()!=ARM_BONDTEC)
		{
			result.setMsg ("ARM_ERR: Only available for bonds");
			return ARM_KO;
		}

		if (!((ARM_BondTEC*)sec)->GetPfTEC())
		{
			result.setMsg ("ARM_ERR: Cover unavailable without portfolio specified for TEC Bond");
			return ARM_KO;
		}
		
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		if (mod->GetZeroCurve()->GetName()!=ARM_ZERO_VASICEK)
		{
			result.setMsg ("ARM_ERR: Only available for Vasicek Model");
			return ARM_KO;
		}

		((ARM_BondTEC *)sec)->CptTecHdg(mod);
		Cover = ((ARM_BondTEC *)sec)->GetTecHdg();

		result.setLong(((ARM_BondTEC *)sec)->GetTecHdg()->GetSize());

		for (long j=0; j<Cover->GetSize(); j++)
			result.setArray(Cover->Elt(j),j);

		//if (Cover) //DAM
		//	delete Cover;
		//Cover = NULL;

		return ARM_OK;
	}
 
	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_FreeObject (long secId, ARM_result& result)
{
    if (LOCAL_PERSISTENT_OBJECTS)
    {
		LOCAL_PERSISTENT_OBJECTS->FreeObject(secId);

		result.setLong(secId);

		return ARM_OK;
    }
	else
		return ARM_KO;
}


long ARMLOCAL_FreeAllObjects (ARM_result& result)
{
    if (LOCAL_PERSISTENT_OBJECTS)
    {
		LOCAL_PERSISTENT_OBJECTS->FreeAllObjects();

		result.setLong(0);

		return ARM_OK;
    }
	else
		return ARM_KO;
}

long ARMLOCAL_NextBusinessDay (double date, const CCString& cal, long days, ARM_result& result)
{
	char nDate[30];	
	char sDate[11];
	char* ccy = cal.GetStr();

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(date,sDate);

        char buf[30];
        ARM_Date myDate(sDate);

		if ( days > 0 )
		{
			myDate.NextBusinessDay(days, ccy).JulianToStrDate(buf);
		}
		else if ( days == 0)
		{
			if (myDate.IsBusinessDay(ccy))
			{
				myDate.JulianToStrDate(buf);
			}
			else
			{
				myDate.NextBusinessDay(1, ccy).JulianToStrDate(buf);
			}
		}
		else
		{
			myDate.PreviousBusinessDay(-days, ccy).JulianToStrDate(buf);
		}

		strcpy(nDate, (char *) buf);

		result.setString(nDate);

		if (ccy)
			free(ccy);
		ccy = NULL;

		return ARM_OK;
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (ccy)
			free(ccy);
		ccy = NULL;

        ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (ccy)
			free(ccy);
		ccy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_IsBusinessDay (double date, const CCString& isoccyname, ARM_result& result)
{
	long flag;

	char sDate[11];
	char* ccy = isoccyname.GetStr();

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

		flag = ((ARM_Date) sDate).IsBusinessDay(ccy);

		result.setLong(flag);

		if (ccy)
			free(ccy);
		ccy = NULL;

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (ccy)
			free(ccy);
		ccy = NULL;

        ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (ccy)
			free(ccy);
		ccy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_FutDelivery (const CCString& fut,
						   const CCString& currency,
						   ARM_result& result)
{
	char sDate[11];
	char buf[30];

	CCString msg ("");

	try
	{
		ARM_Date matDate;
		int month, year;

		GetMonthYearFromExpiryDate(fut, &month, &year);
		matDate.ChangeDate(1, month, year);
		matDate.PiborDelivery((const char*) currency);
		matDate.JulianToStrDate(buf),
		strcpy(sDate, (char *) buf);

		result.setString(sDate);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

        ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_ARM_Accrued (long secId, double fwdDate, long modId, ARM_result& result)
{
	double accrued;
	ARM_Security* sec = NULL;
	ARM_Model* mod = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(fwdDate,sDate);

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}
 
		if ( modId == ARM_NULL_OBJECT )
		{
		   mod = NULL;
		}
		else
		{
			mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
			{
				result.setMsg ("ARM_ERR: Model is not of a good type");
				return ARM_KO;
			}
		}
 
		if ( modId != ARM_NULL_OBJECT )
		{
			sec->SetModel(mod);
		}

		sec->SetSettlement((ARM_Date) sDate);

		accrued = sec->ComputeAccrued();

		result.setDouble(accrued);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_IMPLIEDVOL (long instId, long modelId, double price, bool LnVol, ARM_result& result)
{
	ARM_Security *security = NULL;

	ARM_Model* mod = NULL;
	double vol = 0.; 
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
    
	try
    {
		security = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(security, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}
  
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}
 
        security->SetModel(mod);
 
		if( LnVol )
			vol = security->ComputeImpliedVol(price);

		else
			vol = security->ComputeImpliedVolNor(price);

		result.setDouble(vol);

		return ARM_OK;
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_ARM_View (long instId, ARM_result& result, bool aXMLResult)
{
	ARM_Object* sec = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		sec = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);

		if (sec == NULL)
		{
			result.setMsg ("ARM_ERR: Unknown or Null object");
			return ARM_KO;
		}

		char username[100];
		DWORD nbChar = sizeof(username);

		GetUserName(username,&nbChar);

		CCString	vFileId((CCString)"123" + (CCString)username);
		if(aXMLResult)
			vFileId += ".xml";

		char* sFile = (char*)vFileId;
		sec->View(sFile);
		if (sFile)
			delete sFile;
    }

	catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


long ARMLOCAL_SetNotional (long secId,
						   long rId,
						   double percentRemainder,
						   ARM_result& result,
						   int interpolDates) //default : K_PAY_DATES (=3)
{
	ARM_Security* sec = NULL;
	ARM_ReferenceValue* ref = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(rId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ref, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		sec->SetAmount(ref,percentRemainder);
		//dates type on which notional will be interpolated
		sec->SetInterpolDates(interpolDates);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_GetExpiry (long secId, ARM_result& result)
{
	ARM_Security* sec = NULL;
	ARM_Date expiry;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		expiry = sec->GetExpiryDate();

		char strDate[30];

		expiry.JulianToStrDate(strDate);

		result.setString(strDate);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_Sensitivity (long secId, long modId, long paramId, ARM_result& result)
{
	ARM_Security* sec = NULL;
	ARM_Model* mod = NULL;
	double sensitivity;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

	   sec->SetModel(mod);

	   sensitivity = sec->ComputeSensitivity(paramId);

	   result.setDouble(sensitivity);

	   return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_GetFRMShortRateVols(long modId, VECTOR<double>* matu, VECTOR<double>* rate, 
                                    ARM_result& result)
{
	ARM_Model* mod = NULL;
	ARM_FrmAna* frmmod = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		if ( !IsFrmFamily(mod->GetName()) || IsFrmNumeric(mod->GetName()) )
    		throw Exception (__LINE__, __FILE__,
				ERR_INVALID_ARGUMENT, "Only available for FRM analytical models");

		frmmod = dynamic_cast<ARM_FrmAna *> (mod);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(frmmod, ARM_FRM_ANA) == 0)
		{
			result.setMsg ("ARM_ERR: FRMANA Model is not of a good type");
			return ARM_KO;
		}

		ARM_Vector*   dates = NULL;
		ARM_Vector*   vols = NULL;

		vols = frmmod->GetSpotVol();
		dates = frmmod->GetResetYearTerms();

		if (dates)
		{
			long size = dates->GetSize();

			for ( int i = 0; i < size; i++)
			{
				matu->push_back((*dates)[i]);
				rate->push_back((*vols)[i]);
			}
		}
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




long ARMLOCAL_XCccyAdjustment( long startDate,
                          long endDate, 
				          long payFreq,
                          const CCString& domCcy,
				          long forIndexType,
				          const CCString& forCcy,
                          long spreadsId,
                          long zcDomId,
                          long discDomId,
                          long zcForId,
                          long discForId,
                          double FX, long couponId, 
                          long domDc, long forDc,
                          ARM_result& result, long objId)
{
	long refvalId;

    ARM_ReferenceValue *createdAdj = NULL, *refv = NULL;
	ARM_ReferenceValue* spreads = NULL;
    ARM_ZeroCurve* zcDom = NULL;
	ARM_ZeroCurve* discDom = NULL;
	ARM_ZeroCurve* zcFor = NULL;
	ARM_ZeroCurve* discFor = NULL;
    ARM_ReferenceValue* coupon = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];

	char* sDomCcy = NULL;
	char* sForCcy = NULL;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		sDomCcy = (char*) domCcy;
		sForCcy = (char*) forCcy;

		zcDom = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcDomId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcDom, ARM_ZERO_CURVE) == 0)
		{
			if (sDomCcy)
				delete sDomCcy;
			sDomCcy = NULL;

			if (sForCcy)
				delete sForCcy;
			sForCcy = NULL;

			result.setMsg ("ARM_ERR: zcDom is not of a good type");
			return ARM_KO;
		}

		discDom = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(discDomId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(discDom, ARM_ZERO_CURVE) == 0)
		{
			if (sDomCcy)
				delete sDomCcy;
			sDomCcy = NULL;

			if (sForCcy)
				delete sForCcy;
			sForCcy = NULL;

			result.setMsg ("ARM_ERR: discDom is not of a good type");
			return ARM_KO;
		}

		zcFor = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcForId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcFor, ARM_ZERO_CURVE) == 0)
		{
			if (sDomCcy)
				delete sDomCcy;
			sDomCcy = NULL;

			if (sForCcy)
				delete sForCcy;
			sForCcy = NULL;

			result.setMsg ("ARM_ERR: zcFor is not of a good type");
			return ARM_KO;
		}

		discFor = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(discForId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(discFor, ARM_ZERO_CURVE) == 0)
		{
			if (sDomCcy)
				delete sDomCcy;
			sDomCcy = NULL;

			if (sForCcy)
				delete sForCcy;
			sForCcy = NULL;

			result.setMsg ("ARM_ERR: discFor is not of a good type");
			return ARM_KO;
		}

		spreads = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spreads, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDomCcy)
				delete sDomCcy;
			sDomCcy = NULL;

			if (sForCcy)
				delete sForCcy;
			sForCcy = NULL;

			result.setMsg ("ARM_ERR: spreads is not of a good type");
			return ARM_KO;
		}

		if (couponId >=0)
		{
			coupon = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(couponId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(coupon, ARM_REFERENCE_VALUE) == 0)
			{
				if (sDomCcy)
					delete sDomCcy;
				sDomCcy = NULL;

				if (sForCcy)
					delete sForCcy;
				sForCcy = NULL;

				result.setMsg ("ARM_ERR: coupon is not of a good type");
				return ARM_KO;
			}
		}

		createdAdj = XCcyAdjustment((ARM_Date)sStartDate, (ARM_Date)sEndDate,
									payFreq, sDomCcy,
									(ARM_INDEX_TYPE)forIndexType,
									sForCcy,
									spreads,
									zcDom,
									discDom,
									zcFor,
									discFor,
									FX,
									coupon,
									domDc,
									forDc);

		if (sDomCcy)
			delete sDomCcy;
		sDomCcy = NULL;

		if (sForCcy)
			delete sForCcy;
		sForCcy = NULL;

		if (createdAdj == NULL)
		{
			result.setMsg ("ARM_ERR: refValue is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			refvalId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdAdj);

			if (refvalId == RET_KO)
			{
				if (createdAdj)
					delete createdAdj;
				createdAdj = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(refvalId);

			return ARM_OK;
		}
		else
		{
			refv = (ARM_ReferenceValue  *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refv, ARM_REFERENCE_VALUE) == 1)
			{
				if (refv)
				{
					delete refv;
					refv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdAdj, objId);

				return ARM_OK;
			}
			else
			{
				if (createdAdj)
					delete createdAdj;
				createdAdj = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdAdj)
			delete createdAdj;
		createdAdj = NULL;

		if (sDomCcy)
			delete sDomCcy;
		sDomCcy = NULL;

		if (sForCcy)
			delete sForCcy;
		sForCcy = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdAdj)
			delete createdAdj;
		createdAdj = NULL;

		if (sDomCcy)
			delete sDomCcy;
		sDomCcy = NULL;

		if (sForCcy)
			delete sForCcy;
		sForCcy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ADJUSTTOBUSDATE (double date, const CCString& currency, long ruleId, ARM_result& result)
{
	char nDate[30];
	char sDate[11];
	char* ccy = currency.GetStr();

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

        char buf[30];

		ARM_Date adjDate;
		ARM_Date lDate(sDate);
		lDate.AdjustToBusDate(ccy, ruleId);

		adjDate = lDate;

		adjDate.JulianToStrDate(buf);

		strcpy(nDate, (char *) buf);

		result.setString(nDate);

		if (ccy)
			free(ccy);
		ccy = NULL;

		return ARM_OK;
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (ccy)
			free(ccy);
		ccy = NULL;

        ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (ccy)
			free(ccy);
		ccy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_FwdPrice (long secId, long modId, double fwdDate, ARM_result& result)
{
	ARM_Security* sec=NULL;
	ARM_Model* mod=NULL;
	double fprice;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sFwdDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(fwdDate,sFwdDate);

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		sec->SetModel(mod);

		sec->SetForwardDate((ARM_Date) sFwdDate);

		fprice = sec->ComputeForwardPrice();

		result.setDouble(fprice);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_CvSensitivity (long secId,
							 long modId,
							 long paramId,
							 ARM_result& result)
{
	ARM_Security* sec=NULL;
	ARM_Model* mod=NULL;
	double sensitivity;
	double bpShift=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}
		
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		sec->SetModel(mod);

		sensitivity = sec->ComputeSensitivityCurve(paramId, 0, NULL, bpShift);
		
		result.setDouble(sensitivity);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_ARM_BetweenDates (long date1,
								long date2,
								long daycountId,
								long isYearFrac,
								ARM_result& result)
{
	double diffDates(0.);

	char sDate1[11];
	char sDate2[11];

	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(date1,sDate1);
		Local_XLDATE2ARMDATE(date2,sDate2);

		ARM_Date d1(sDate1);
		ARM_Date d2(sDate2);

		if (isYearFrac)
		{
		   diffDates = CountYears(daycountId, d1, d2);
		}
		else
		{
		   diffDates = DaysBetweenDates(daycountId, d1, d2);
		}

		result.setDouble(diffDates);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_CountBusinessDays (long date1,
									 long date2,
									 CCString calendar,
									 ARM_result& result)
{
	double diffDates(0.);

	char sDate1[11];
	char sDate2[11];
	
	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(date1,sDate1);
		Local_XLDATE2ARMDATE(date2,sDate2);

		const ARM_Date d1(sDate1);
		const ARM_Date d2(sDate2);

	    diffDates = CountBusinessDays(d1, d2, (char*) calendar);
		   
		result.setDouble(diffDates);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_ADDYEARS  (long date,
							 long nb,
							 long ruleId,
							 const CCString& Ccy,
							 ARM_result& result)
{
	char sDate[11];
	CCString msg (" ");

    char buf[30];
    char nDate[30];

	char* sCCY = NULL;

    try
    {
 		Local_XLDATE2ARMDATE(date,sDate);
        ARM_Date date(sDate);

		sCCY = (char*)Ccy;

		if ( nb != 0 )
		{
			date.AddYears(nb);
			if (ruleId != 0)
				date.AdjustToBusDate(sCCY,ruleId);
		}

		delete sCCY;

		date.JulianToStrDate(buf);

		strcpy(nDate, (char *) buf);
		result.setString(nDate);

		return ARM_OK;
    }

    catch(Exception& x)
    {
        x.DebugPrint();

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_ADDMONTHS  (long date,
							  long nb,
							  long ruleId,
							  const CCString& Ccy,
							  ARM_result& result)
{
	char sDate[11];
	CCString msg (" ");

    char buf[30];
    char nDate[30];

	char* sCCY = NULL;

    try
    {
 		Local_XLDATE2ARMDATE(date,sDate);
        ARM_Date date(sDate);

		sCCY = (char*)Ccy;

		if ( nb != 0 )
		{
			date.AddMonths(nb);
			if (ruleId != 0)
				date.AdjustToBusDate(sCCY,ruleId);
		}

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		date.JulianToStrDate(buf);

        strcpy(nDate, (char *) buf);
		result.setString(nDate);

		return ARM_OK;
    }

    catch(Exception& x)
    {
        x.DebugPrint();

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_FxConvert (const CCString& ccy1,
						 const CCString& ccy2,
						 double asOfDate,
						 double amount,
						 const CCString& cvname,
						 ARM_result& result)
{
	CCString msg (" ");

	if (ccy1 == ccy2)
	{
		result.setDouble(1.0);
		return ARM_OK;
	}

	if (GetDataRetrieverVersion() >= ETKRETRIEVER)
	{
		char sDate[11];
		Local_XLDATE2ARMDATE(asOfDate,sDate);
		ARM_Date asOf(sDate);

		try
		{
			double val = ARMLOCAL_GetFxCurve(ccy1,
											 ccy2,
											 asOf,
											 amount,
											 cvname);

			result.setDouble(val);
			return ARM_OK;
		}
		catch (Exception& x)
		{
			x.DebugPrint();

			ARM_RESULT();
		}
		/// catch the rest
		catch (...)
		{
			result.setMsg ("ARM_ERR: unrecognized failure");
			return ARM_KO;
		}
	}
	else
	{
		FILE *Fp = NULL;

		double val1=1.; // valeur devise1/USD
		double val2=1.; // valeur USD/devise2

		CCString myRepertory ((CCString)(armlocal_init->data_folder.c_str()) + FX_SUMMIT_FILE_LOCATION);
		CCString myRepertory2 ((CCString)(armlocal_init->data_folder.c_str()) + FX_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

		CCString FileName1 = myRepertory + (CCString)"FX_";
		CCString FileNameHisto1 = myRepertory2 + (CCString)"FX_";
		CCString FileName2 = myRepertory + (CCString)"FX_";
		CCString FileNameHisto2 = myRepertory2 + (CCString)"FX_";

		CCString Suffixe;
		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			Suffixe = CCString ("_.");
		}
		else
		{
			Suffixe = CCString ("_") + cvname + CCString (".");
		}

		char sDate[11];

		Local_XLDATE2ARMDATE(asOfDate,sDate);

		ARM_Date myDate(sDate);

		ARM_Date dateToday; // Constructeur par Défaut donne la date du jour
	//	dateToday.SysToday();// donne la date systeme du jour

		char buffer[50];

		if (myDate > dateToday)
		{
			result.setMsg ("ARM_ERR: Invalid Date");
			return ARM_KO;
		}

		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			_ltoa(myDate.GetYear(),buffer,10);
			Suffixe += (CCString) (buffer);

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				Suffixe += "0";
			}
			Suffixe += (CCString) buffer;

			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				Suffixe += "0";
			}
			Suffixe += (CCString) buffer;
		}
		else
		{
			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				Suffixe += "0";
			}
			Suffixe += (CCString) buffer;

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				Suffixe += "0";
			}
			Suffixe += (CCString) buffer;

			_ltoa(myDate.GetYear(),buffer,10);
			Suffixe += (CCString) (buffer+2);
		}

		Suffixe += (CCString) ".000";
		
		vector<string> currencyList = armlocal_init->currencyList;
		vector<int> indirectList = armlocal_init->indirectList;
		int indirect1, indirect2;

		if (strcmp((const char*) ccy1,"USD") != 0)
		{
			for (int i=0; i< currencyList.size(); i++)
			{
				if (strcmp((const char*) currencyList[i].c_str(),(const char*) ccy1) == 0)
				{
					indirect1 = indirectList[i];
					i = currencyList.size();
				}
			}

			if (indirect1 == 1)
			{
				FileName1 = FileName1 + ccy1 + "_USD" + Suffixe;
				FileNameHisto1 = FileNameHisto1 + ccy1 + "_USD" + Suffixe;
			}
			else
			{
				FileName1 = FileName1 + "USD_" + ccy1 + Suffixe;
				FileNameHisto1 = FileNameHisto1 + "USD_" + ccy1 + Suffixe;
			}

			if ((Fp = fopen(FileName1,"r")) == NULL)
			{
				if ((Fp = fopen(FileNameHisto1,"r")) == NULL)
				{
					if (GetFallBackDataRetrieverVersion() == 0)
					{
						result.setMsg ("ARM_ERR: File not found");
						return ARM_KO;
					}
					else
					{
						char sDate[11];
						Local_XLDATE2ARMDATE(asOfDate,sDate);
						ARM_Date asOf(sDate);

						try
						{
							double val = ARMLOCAL_GetFxCurve(ccy1,
															 ccy2,
															 asOf,
															 amount,
															 cvname);

							result.setDouble(val);
							return ARM_OK;
						}
						catch (Exception& x)
						{
							x.DebugPrint();

							ARM_RESULT();
						}
						/// catch the rest
						catch (...)
						{
							result.setMsg ("ARM_ERR: unrecognized failure");
							return ARM_KO;
						}
					}
				}
			}

			int rc = 0;

			char sEch[50];

			rc = fscanf(Fp, "%s",sEch);
			rc = fscanf(Fp, "%s",buffer);

			fclose (Fp);

			if (indirect1 == 1)
				val1 = atof((const char*) buffer);
			else
				val1 = 1. / atof((const char*) buffer);
		}

		if (strcmp((const char*) ccy2,"USD") != 0)
		{
			for (int i=0; i< currencyList.size(); i++)
			{
				if (strcmp((const char*) currencyList[i].c_str(),(const char*) ccy2) == 0)
				{
					indirect2 = indirectList[i];
					i = currencyList.size();
				}
			}

			if (indirect2 == 1)
			{
				FileName2 = FileName2 + ccy2 + "_USD" + Suffixe;
				FileNameHisto2 = FileNameHisto2 + ccy2 + "_USD" + Suffixe;
			}
			else
			{
				FileName2 = FileName2 + "USD_" + ccy2 + Suffixe;
				FileNameHisto2 = FileNameHisto2 + "USD_" + ccy2 + Suffixe;
			}

			if ((Fp = fopen(FileName2,"r")) == NULL)
			{
				if ((Fp = fopen(FileNameHisto2,"r")) == NULL)
				{
					if (GetFallBackDataRetrieverVersion() == 0)
					{
						result.setMsg ("ARM_ERR: File not found");
						return ARM_KO;
					}
					else
					{
						char sDate[11];
						Local_XLDATE2ARMDATE(asOfDate,sDate);
						ARM_Date asOf(sDate);

						try
						{
							double val = ARMLOCAL_GetFxCurve(ccy1,
															 ccy2,
															 asOf,
															 amount,
															 cvname);

							result.setDouble(val);
							return ARM_OK;
						}
						catch (Exception& x)
						{
							x.DebugPrint();

							ARM_RESULT();
						}
						/// catch the rest
						catch (...)
						{
							result.setMsg ("ARM_ERR: unrecognized failure");
							return ARM_KO;
						}
					}
				}
			}

			int rc = 0;

			char sEch[50];

			rc = fscanf(Fp, "%s",sEch);
			rc = fscanf(Fp, "%s",buffer);

			fclose (Fp);

			if (indirect2 == 0)
				val2 = atof((const char*) buffer);
			else
				val2 = 1. / atof((const char*) buffer);
		}

		result.setDouble(val1 * val2 *amount);
		return ARM_OK;
	}
}



long ARMLOCAL_ARM_ADDPERIOD  (double date,
							  long freq,
							  const CCString& ccy,
							  long nbPeriods,
							  long adjRuleId,
							  long goToEndOfMonth,
							  ARM_result& result)
{
	char sDate[11];
	CCString msg (" ");

    char buf[30];
    char nDate[30];

	char* sCcy = NULL;

    try
    {
		if (freq == K_DAILY)
		{
			long retCode = ARMLOCAL_NextBusinessDay (date, ccy, nbPeriods, result);
			return retCode;
		}
		else
		{
			Local_XLDATE2ARMDATE(date,sDate);
			ARM_Date myDate(sDate);

			sCcy = (char*) ccy;

			myDate.AddPeriodMult(freq,nbPeriods,sCcy,goToEndOfMonth);

			if (adjRuleId != (int)0)
				myDate.AdjustToBusDate(sCcy,adjRuleId);

			myDate.JulianToStrDate(buf);

			if (sCcy)
				delete sCcy;
			sCcy = NULL;

			strcpy(nDate, (char *) buf);
			result.setString(nDate);

			return ARM_OK;
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();

		if (sCcy)
			delete sCcy;
		sCcy = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (sCcy)
			delete sCcy;
		sCcy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_Price_OptUnder (long secId,
								  long modId,
								  ARM_result& result)
{
	ARM_Object* security = NULL;
	ARM_Security* sec = NULL;

	ARM_Model* mod = NULL;
	ARM_PricingModel* pricingModel = NULL;
	ARM_IFPricer *pricer = NULL;
	double price;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
			
		/// gets the security
		security = LOCAL_PERSISTENT_OBJECTS->GetObject(secId);
		/// get model
		pricingModel = (ARM_PricingModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);
		if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(security, ARM_GENCALCULATOR))
        {
			ARM_GenCalculator* genCalc=dynamic_cast<ARM_GenCalculator*>(security);
			if( !genCalc)
				ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into generic calculator failed!" );

            /// force the model of the calculator then price
			if(pricingModel && LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pricingModel, ARM_PRICINGMODEL))
			{
				ARM_PricingModel* pmod=dynamic_cast<ARM_PricingModel*>(pricingModel);
				if( !pmod )
					ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into pricing model failed!" );
				genCalc->SetPricingModel( ARM_PricingModelPtr((ARM_PricingModel*)pmod->Clone()) );
			}

            ARM_Vector* prices	= genCalc->ComputeAll();
			//CC_NS(ARM, ARM_AutoCleaner)<ARM_Vector> Hold(prices);

		   for (int i = 0; i < prices->GetSize(); i++)
			   result.setArray(prices->Elt(i), i);

		   result.setDouble(prices->GetSize());

		   delete prices;

		   return(ARM_OK);
        }
		
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);
		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0 )
		{
		   result.setMsg ("ARM_ERR: Security is not of a good type");			
		   return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);	
		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0 )
		{
		   result.setMsg ("ARM_ERR: Model is not of a good type");		
		   return ARM_KO;
		}

		

        if ( mod->GetName() == ARM_FRMMARKOVTREE )
        {
		    pricer = new ARM_IFPricer(sec, mod);

            ARM_Vector* resPrices = pricer->ComputePrices();

            for (int i = 0; i < resPrices->GetSize(); i++)
            {
		        result.setArray(resPrices->Elt(i),i);
            }

            result.setDouble(resPrices->GetSize());

		    if (pricer)
		       delete pricer;
		    pricer = NULL;

		    if (resPrices)
		       delete resPrices;
		    resPrices = NULL;

		    return(ARM_OK);
        }
        else if ( mod->GetName() == ARM_TREE3F )
        {
           sec->SetModelVariable(mod);

           ARM_Vector* res = ((ARM_PowerReverse *) sec)->ComputeAll();

           for (int i = 0; i < res->GetSize(); i++)
           {
		       result.setArray(res->Elt(i), i);
           }

           result.setDouble(res->GetSize());

		   if (res)
		      delete res;
		   res = NULL;

		   return(ARM_OK);
        }
        else
        {
		    pricer = new ARM_IFPricer(sec, mod);
            int i=0;

		    price = pricer->Price();

		    result.setArray(price,i);
            ++i;

            if ( mod->GetName() == ARM_FRMMCMODEL)
            {
                ARM_MCModel *mcmod=dynamic_cast<ARM_MCModel *>(mod);

                double stdDev = mcmod->GetStdDev();

                result.setArray(stdDev,i);
                ++i;
                result.setDouble(double (i));
            
            }

		    if ( sec->GetName() == ARM_OPTION )
		    {           
		    
			   ARM_Option* TmpSec = (ARM_Option*) sec;

			   result.setArray(TmpSec->GetUnderlying()->GetPrice(),i);
               ++i;
               
		       result.setDouble(double (i));
		    }
            else
            
		        result.setDouble(double (i));

		    if (pricer)
		       delete pricer;
		    pricer = NULL;

		    result.setDouble(2.);

		    return(ARM_OK);
        }
	}
 
	catch(Exception& x)
	{
		if (pricer)
		   delete pricer;
		pricer = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		if (pricer)
		   delete pricer;
		pricer = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

}


long ARMLOCAL_ParallelShift (long crvId,
							 double value,
							 ARM_result& result,
							 long objId)
{
	long curveId;

	ARM_Object* inCurve = NULL;
	ARM_ZeroCurve* createdZcCurve = NULL;
	ARM_VolCurve* createdVolCurve = NULL;
	ARM_Object* oldCurve = NULL;

	ARM_Object* object = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		inCurve = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(crvId);

		if (LocalPersistent::LOCAL_IS_CURVE_CLASS_OK(inCurve) == 0)
		{
			result.setMsg ("ARM_ERR: In Curve is not a curve");
			return ARM_KO;
		}

		if (inCurve->GetRootName() == ARM_ZERO_CURVE)
		{
			createdZcCurve = (ARM_ZeroCurve*) inCurve->Clone();
			createdZcCurve->ParallelShift(value);
			object = createdZcCurve;
		}
		else if (inCurve->GetRootName() == ARM_VOL_CURVE)
		{
			createdVolCurve = (ARM_VolCurve*) inCurve->Clone();
			createdVolCurve->ParallelShift(value);
			object = createdVolCurve;
		}
		else
		{
			result.setMsg ("ARM_ERR: Pb");
			return ARM_KO;
		}

		if (object == NULL)
		{
			result.setMsg ("ARM_ERR: created Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)object);

			if (curveId == RET_KO)
			{
				if (object)
					delete object;
				object = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			oldCurve = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_CURVE_CLASS_OK(oldCurve) == 1)
			{
				if (oldCurve)
				{
					delete oldCurve;
					oldCurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)object, objId);

				return ARM_OK;
			}

			else
			{
				if (object)
					delete object;
				object = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (object)
			delete object;
		object = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (object)
			delete object;
		object = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ClonedAndSetNotional (long secId,
									long rId,
									double percentRemainder,
									ARM_result& result,
									int interpolDates, // default value = K_PAY_DATES (=3)
									long objId) 
{
	long clonedId;

	ARM_Security* sec = NULL;
	ARM_Security* oldSec = NULL;
	ARM_Security* clonedSec = NULL;

	ARM_ReferenceValue* ref= NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec,ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Object to clone is not a security");
			return ARM_KO;
		}

		ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(rId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ref,ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Notional is not a Reference Value");
			return ARM_KO;
		}

		int calcMethod = ref->GetCalcMethod();
		int extrpMeth = ref->GetExtrapolMeth();

		ref->SetCalcMethod(K_STEPUP_RIGHT);
		ref->SetExtrapolMeth(1);

		clonedSec = (ARM_Security *) sec->Clone();
		clonedSec->SetModel(NULL);
		clonedSec->SetAmount(ref,percentRemainder);
		clonedSec->SetInterpolDates(interpolDates);

		ref->SetCalcMethod(calcMethod);
		ref->SetExtrapolMeth(extrpMeth);

		if (clonedSec == NULL)
		{
			result.setMsg ("ARM_ERR: created security is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			clonedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec);

			if (clonedId == RET_KO)
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clonedId);

			return ARM_OK;
		}
		else
		{
			oldSec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(oldSec,ARM_SECURITY) == 1)
			{
				if (oldSec)
				{
					delete oldSec;
					oldSec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec, objId);

				return ARM_OK;
			}

			else
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




long ARMLOCAL_ClonedAndSet (long secId,
							long valtosetType,
							double valtoset,
							const CCString typetoset,
							ARM_result& result,
							long objId)
{
	long clonedId;

	ARM_SpreadOption* sec = NULL;
	ARM_SpreadOption* oldSec = NULL;
	ARM_SpreadOption* clonedSec = NULL;

	ARM_Portfolio* port = NULL;
	ARM_Portfolio* oldPort = NULL;
	ARM_Portfolio* clonedPort = NULL;

	ARM_Object* newObj;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		sec = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SPREADOPTION) == 0)
		{
			sec = NULL;

			port = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(port,ARM_PORTFOLIO) == 0)
			{
				result.setMsg ("ARM_ERR: Object to clone has to be a spreadoption or a portfolio");
				return ARM_KO;
			}
		}

		if (sec)
		{
			clonedSec = (ARM_SpreadOption *) sec->Clone();
			newObj = clonedSec;
		}
		else
		{
			clonedPort = (ARM_Portfolio *) port->DeepClone();
			newObj = clonedPort;
		}

		if (typetoset == "STRIKE")
		{
			if (valtosetType == 0)
			{
				if (clonedSec)
					clonedSec->SetStrike(valtoset);
				else
					clonedPort->SetStrike(valtoset);
			}
			else
			{
				ARM_ReferenceValue* strike;

				strike = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(valtoset));

				if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(strike, ARM_REFERENCE_VALUE) == 0)
				{
					result.setMsg ("ARM_ERR: strike is not of a good type");
					return ARM_KO;
				}

				if (clonedSec)
					clonedSec->SetStrikes(strike);
				else
					clonedPort->SetStrike(valtoset);
			}
		}
		if (typetoset == "STRIKE_SHIFT")
		{
			if (valtosetType == 0)
			{
				if (clonedSec)
					clonedSec->SetShiftStrike(valtoset);
				else
					clonedPort->SetShiftStrike(valtoset);
			}
			else
			{
				result.setMsg ("ARM_ERR: Shift must be a numeric");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: invalid Type; Valid Is STRIKE or STRIKE_SHIFT");
			return ARM_KO;
		}

		if (newObj == NULL)
		{
			result.setMsg ("ARM_ERR: created security is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			clonedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObj);

			if (clonedId == RET_KO)
			{
				if (newObj)
					delete newObj;
				newObj = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clonedId);

			return ARM_OK;
		}
		else
		{
			oldSec = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSec,ARM_SPREADOPTION) == 1)
			{
				if (oldSec)
				{
					delete oldSec;
					oldSec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObj, objId);

				return ARM_OK;
			}
			else
			{
				oldPort = (ARM_Portfolio*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

				if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldPort,ARM_PORTFOLIO) == 1)
				{
					if (oldPort)
					{
						delete oldPort;
						oldPort = NULL;
					}

					LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObj, objId);

					return ARM_OK;
				}
				else
				{
					if (newObj)
						delete newObj;
					newObj = NULL;

					result.setMsg ("ARM_ERR: previous object is not of a good type");
					return ARM_KO;
				}
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newObj)
			delete newObj;
		newObj = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newObj)
			delete newObj;
		newObj = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_Clone (long objectId,
					 ARM_result& result,
					 long objId)
{
	long clonedId;

	ARM_Object* sec = NULL;
	ARM_Object* oldSec = NULL;
	ARM_Object* clonedSec = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		sec = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(objectId);

		if (sec)
			clonedSec = (ARM_Object *) sec->Clone();

		if (clonedSec == NULL)
		{
			result.setMsg ("ARM_ERR: created security is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			clonedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec);

			if (clonedId == RET_KO)
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clonedId);

			return ARM_OK;
		}
		else
		{
			oldSec = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec, objId);

			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




long ARMLOCAL_DoPastReset (long secId,
						   VECTOR<CCString>& resetMgrIds,
						   long AsOf,
						   ARM_result& result,
						   long objId)
{
	long clonedId;

	ARM_Security* sec = NULL;
	ARM_Security* oldSec = NULL;
	ARM_Security* clonedSec = NULL;

	vector<ARM_ResetManager*> resetMgrs;
	ARM_ResetManager* resetMgr = NULL;

	char sDate[11];

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(AsOf,sDate);
        ARM_Date myDate(sDate);

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Object to clone is not a security");
			return ARM_KO;
		}

		if (resetMgrIds.size() < 1)
		{
			result.setMsg ("ARM_ERR: size of resetMgrs must be at least 1");
			return ARM_KO;
		}

		for (int i = 0; i < resetMgrIds.size(); i++)
		{
			resetMgr = (ARM_ResetManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(resetMgrIds[i]));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(resetMgr, ARM_RESETMANAGER) == 0)
			{
				result.setMsg ("ARM_ERR: Object is not a resetMgr");
				return ARM_KO;
			}
			resetMgrs.push_back(resetMgr);
		}

		clonedSec = (ARM_Security*) sec->Clone();
		clonedSec->DoPastReset(resetMgrs,myDate);

		if (clonedSec == NULL)
		{
			result.setMsg ("ARM_ERR: created security is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			clonedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec);

			if (clonedId == RET_KO)
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clonedId);

			return ARM_OK;
		}
		else
		{
			oldSec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(oldSec,ARM_SECURITY) == 1)
			{
				if (oldSec)
				{
					delete oldSec;
					oldSec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec, objId);

				return ARM_OK;
			}

			else
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_FIXRATES (double secId,
						VECTOR<double>& rate,
						ARM_result& result,
						long objId)
{
	long clonedId;

	ARM_Security* sec = NULL;
	ARM_Security* oldSec = NULL;
	ARM_Security* clonedSec = NULL;

	ARM_SwapLeg* leg = NULL;

	ARM_Vector* vFixing = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Object to clone is not a security");
			return ARM_KO;
		}

		vFixing = CreateARMVectorFromVECTOR(rate);

		switch ( sec->GetName() )
		{
		    case ARM_SWAPLEG:
		    case ARM_CMSLEG:
		    case ARM_CMTLEG:
		    case ARM_TMLEG:
			{
				ARM_SwapLeg* clonedLeg = (ARM_SwapLeg*) sec->Clone();
				clonedSec = clonedLeg;
				leg = clonedLeg;
			};
			break;

		    case ARM_SWAPTION:
			{
				ARM_Swaption* clonedSwaption = (ARM_Swaption*) sec->Clone();
				clonedSec = clonedSwaption;
				leg = ( (ARM_Swap*)clonedSwaption)->GetFloatLeg();
			};
			break;

		    case ARM_SWAP:
			{
				ARM_Swap* clonedSwap = (ARM_Swap*) sec->Clone();
				clonedSec = clonedSwap;
				leg = clonedSwap->GetFloatLeg();
			};
			break;

		    case ARM_CAPFLOOR:
			{
				ARM_CapFloor* clonedCF = (ARM_CapFloor*) sec->Clone();
				clonedSec = clonedCF;
				leg = clonedCF->GetSwapLeg();
			};
			break;

		    default:
			{
				if (vFixing)
					delete vFixing;
				vFixing = NULL;

				result.setMsg ("ARM_ERR: Object to clone is not a instance of expected class");
				return ARM_KO;
			}
			break;
		}

		leg->SetFixRates(vFixing);

		if (vFixing)
			delete vFixing;
		vFixing = NULL;

		if (clonedSec == NULL)
		{
			result.setMsg ("ARM_ERR: created security is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			clonedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec);

			if (clonedId == RET_KO)
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clonedId);

			return ARM_OK;
		}
		else
		{
			oldSec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(oldSec,ARM_SECURITY) == 1)
			{
				if (oldSec)
				{
					delete oldSec;
					oldSec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)clonedSec, objId);

				return ARM_OK;
			}

			else
			{
				if (clonedSec)
					delete clonedSec;
				clonedSec = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (vFixing)
			delete vFixing;
		vFixing = NULL;

		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (vFixing)
			delete vFixing;
		vFixing = NULL;

		if (clonedSec)
			delete clonedSec;
		clonedSec = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_DisplayScheduleDates (long instId,
										long datesType,
										long recId,
										long viewInitExchId,
										ARM_result& result)
{
	ARM_Security* sec = NULL;
	ARM_SwapLeg* leg = NULL;
	ARM_ReferenceValue* refval = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(sec) == 1)
		{
			leg = (ARM_SwapLeg*)sec;

			leg->DisplayScheduleDates(datesType,viewInitExchId,"123");
		}

		else if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAP) == 1) ||
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAPTION) == 1)
			)
		{
			if (((ARM_Swap*) sec)->Get1stLeg()->GetRcvOrPay() == recId)
				leg = ((ARM_Swap*) sec)->Get1stLeg();
			else
				leg = ((ARM_Swap*) sec)->Get2ndLeg();
				
			leg->DisplayScheduleDates(datesType,viewInitExchId,"123");
		}

		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_CAPFLOOR) == 1)
		{
			leg = ((ARM_CapFloor*)sec)->GetSwapLeg();
				
			leg->DisplayScheduleDates(datesType,viewInitExchId,"123");
		}

		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_DIGITAL) == 1)
		{				
			leg = ((ARM_Digital*)sec)->GetSwapLeg();
				
			leg->DisplayScheduleDates(datesType,viewInitExchId,"123");
		}
		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_REFERENCE_VALUE) == 1)
		{				
			refval = (ARM_ReferenceValue*)sec;
				
			refval->DisplayScheduleDates(datesType,viewInitExchId,"123");
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

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}

long ARMLOCAL_ARM_DisplayScheduleValues (long instId,
										 long valuesType,
										 long recId,
										 long modelId,
										 ARM_result& result)
{
	ARM_Security* sec = NULL;
	ARM_SwapLeg* leg = NULL;
	ARM_ReferenceValue* refval = NULL;


	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(sec) == 1)
		{
			if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

			leg = (ARM_SwapLeg*)sec;

			leg->DisplayScheduleValues(valuesType,"123");
		}

		else if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAP) == 1) ||
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAPTION) == 1)
			)
		{
			if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

			if (((ARM_Swap*) sec)->Get1stLeg()->GetRcvOrPay() == recId)
				leg = ((ARM_Swap*) sec)->Get1stLeg();
			else
				leg = ((ARM_Swap*) sec)->Get2ndLeg();
				
			leg->DisplayScheduleValues(valuesType,"123");
		}
		
		else if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_CAPFLOOR) == 1) ||
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_MATCAPFLOOR) == 1))
		{
			if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

            if (valuesType == K_VOL_CAP)
            {
                ((ARM_CapFloor*)sec)->DisplayScheduleValues(valuesType,"123");
            }
            else
            {
			    leg = ((ARM_CapFloor*)sec)->GetSwapLeg();
				    
			    leg->DisplayScheduleValues(valuesType,"123");
            }
		}
	
		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_DIGITAL) == 1)
		{				
			if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

			leg = ((ARM_Digital*)sec)->GetSwapLeg();
				
			leg->DisplayScheduleValues(valuesType,"123");
		}
		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_REFERENCE_VALUE) == 1)
		{				
			refval = (ARM_ReferenceValue*)sec;
				
			refval->DisplayScheduleValues(valuesType,"123");
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

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


long ARMLOCAL_ARM_DisplayReplicPortfolio(long instId,
										 long WeightOrStrike,
                                         long PayoffOrSensi,
                                         long recId,
                                         long modelId,
                                         VECTOR<double>& DataResult,
										 ARM_result& result)
{
	ARM_Security* sec = NULL;
	ARM_SwapLeg* leg = NULL;
    ARM_CapFloor* capFloor = NULL;
    ARM_ReplicPortfolio* replicPort = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    long retCode;

    try
    {
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(instId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(sec) == 1)
		{
            if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

			leg = (ARM_SwapLeg*)sec;
    
            replicPort = ARM_ReplicPortfolio::GetReplicPortfolio(*(leg->GetReplicPortKey()));
		}

		else if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAP) == 1) ||
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_SWAPTION) == 1)
			)
		{
            if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

			if (((ARM_Swap*) sec)->Get1stLeg()->GetRcvOrPay() == recId)
				leg = ((ARM_Swap*) sec)->Get1stLeg();
			else
				leg = ((ARM_Swap*) sec)->Get2ndLeg();
				
			replicPort = ARM_ReplicPortfolio::GetReplicPortfolio(*(leg->GetReplicPortKey()));
		}
		
		else if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec,ARM_CAPFLOOR) == 1)
		{
            if (modelId != ARM_NULL_OBJECT)
				retCode = ARMLOCAL_ARM_Price(instId,modelId,result);

            capFloor = (ARM_CapFloor*)sec;

		    replicPort = ARM_ReplicPortfolio::GetReplicPortfolio(*(leg->GetReplicPortKey()));
		}

		else
		{
			result.setMsg ("ARM_ERR: function not implemented for this object");
			return ARM_KO;
		}

        if (!PayoffOrSensi)
        {
            replicPort = replicPort->GetSensiPortfolio();
        }

        if (replicPort)
        {
            const std::vector<double>& weightVector = replicPort->GetWeightVector();
            const std::vector<double>& strikeVector = replicPort->GetStrikeVector();

            result.setLong(weightVector.size());

            for (int i = 0; i < weightVector.size(); ++i)
            {
                if (WeightOrStrike)
                {
                    DataResult.push_back(weightVector[i]);
                }
                else
                {
                    DataResult.push_back(strikeVector[i]);
                }
            }
        }
    }
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}



long ARMLOCAL_INTERPOL (const VECTOR<double>& vecX,
						const VECTOR<double>& vecY,
						double X,
						long interpId,
						ARM_result& result)
{
	double val;

	if ( vecX.size () != vecY.size () )
	{
		result.setMsg("ARM_ERR: vector X and vector Y must the have same size");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		ARM_Vector vX(vecX.size());
		ARM_Vector vY(vecY.size());

		for (int i = 0; i < vecX.size(); i++)
		{
			vX.Elt(i) = vecX[i];
			vY.Elt(i) = vecY[i];
		}

        if ( interpId == K_SPLINE )
        {
           val = SplineInterpolateFunc(&vX, &vY,
                                       X,
                                       NULL, // SecondDerivCalc = NULL,
                                       1);   // keep2Der = 0
        }
        else
        {
		   val = linInterpol2(&vX,X,&vY);
        }

		result.setDouble(val);
	}

    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


long ARMLOCAL_TRIANGULARINTERPOL (const VECTOR<double>& vecX,
								  const VECTOR<double>& vecY,
								  const VECTOR<double>& matZ,
								  double X,
								  double Y,
								  ARM_result& result)
{
	double val;

	ARM_Vector vX(vecX.size());
	ARM_Vector vY(vecY.size());

	if (matZ.size() != (vecX.size() * vecY.size()))
	{
		result.setMsg ("ARM_ERR: check your matrix dimension");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		ARM_Matrix mZ(vecX.size(),vecY.size());

		for (int i=0;i<vecX.size();i++)
		{
			vX.Elt(i) = vecX[i];

			for (int j=0;j<vecY.size();j++)
			{
				mZ.Elt(i,j) = matZ[i*vecY.size()+j];
			}
		}

		for (int j=0;j<vecY.size();j++)
		{
			vY.Elt(j) = vecY[j];
		}

		val = triangularInterpol(&vX,&vY,&mZ,X,Y);

		result.setDouble(val);
	}

    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}


long ARMLOCAL_KImp (long secId,
					long modId,
					double price,
					long param,
					ARM_result &result)
{
	ARM_Security* opt=NULL;
	ARM_Model* mod=NULL;
	double strike=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		opt = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(opt, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		opt->SetModel(mod);

		strike = opt->ComputeImpliedStrike(price,param);

		result.setDouble(strike);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_ARM_DisplayZC(long zcId,
							ARM_result &result)
{
	ARM_ZeroCurve* zc=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		result.setLong(zc->GetYearTerms()->GetSize());

		for (int i=0;i<zc->GetYearTerms()->GetSize();i++)
		{
			result.setArray(zc->GetYearTerms()->Elt(i),i);
			result.setArray(zc->GetDiscountFactors()->Elt(i),zc->GetYearTerms()->GetSize()+i);
			result.setArray(zc->GetZeroRates()->Elt(i),2*zc->GetYearTerms()->GetSize()+i);
		}

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_BSSpot (long secId,
					  long modId,
					  double date,
					  ARM_result& result)
{
	ARM_Security* sec=NULL;
	ARM_Model* mod=NULL;
	double spot=0.0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		sec->SetModel(mod);

		spot = sec->ComputeBSSpot((ARM_Date)sDate);

		result.setDouble(spot);

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




ARM_Vector* CreateARMMatuVectorFromStrVECTOR(const VECTOR<CCString>& param)
{
	long size = param.size();

	if (size == 0)
		return NULL;

	ARM_Vector* res = new ARM_Vector(size);

	for (int i=0;i<size;i++)
	{
		size_t l = strspn( param[i], " +-.0123456789");
		if( l == strlen( param[i] ) )
		{
			double tmp =  atof( param[i] );
			res->Elt(i) = atof( param[i] );
		}
		else
		{
			char* tmpChar = param[i];
			double tmp =  StringMatuToYearTerm( tmpChar );
			res->Elt(i) = StringMatuToYearTerm( tmpChar );
			delete tmpChar;
		}
	}

	return res;
}




long ARMLOCAL_GETINSTRUMENTFROMSUMMIT (const CCString& idSummit,
									   CCString& typeId,
									   double asOf,
									   const CCString& filter,
									   ARM_result& result,
									   long objId)
{
	if ( (GetDataRetrieverVersion() == FFRETRIEVER) && (GetFallBackDataRetrieverVersion() == FFRETRIEVER) )
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	ARM_Object* newObject = NULL;
	ARM_Object* oldObject = NULL;

	long Id;
	ARM_CLASS_NAME thisClass;

	char sDate[11];

	CCString msg (" ");

	try
	{
		ARM_Date date(MINDATE);

		if (asOf != -1)
		{
			Local_XLDATE2ARMDATE(asOf,sDate);
			ARM_Date tmpdate(sDate);
			date = tmpdate;
		}

		//Temporary modification to test bermuda calculator efficiency.
		CCString xmlResponse;
		if (typeId == "BERM") 
		{	
			typeId = CCString("SWOPT");
			xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,typeId);
			typeId = CCString("BERM");
		}
		else if(typeId == "SWAP.RA")
			xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,"SWAP");
		else if(typeId == "EXOTIC.RA")
			xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,"EXOTIC");
		else
		{
			xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,typeId);		
		}
		if (xmlResponse.GetLen() == 0)
		{
			msg.Set("No trade available with Id=");
			msg += idSummit;
			msg += " and type=";
			msg += typeId;
			result.setMsg(msg);
			return ARM_KO;
		}

#ifdef _DEBUG
		char filename[50];
		if (filter == "ALL")
			sprintf(filename, "C:\\%s - %s.xml", (const char*)typeId, idSummit);
		else
			sprintf(filename, "C:\\%s.%s - %s.xml", (const char*)typeId, (const char*)filter, idSummit);

		FILE* file = fopen(filename, "w");

		if (file)
		{
			fprintf(file, "%s", (const char*)xmlResponse);
			fclose(file);
		}
#endif

		CCString bookName;
		CCString structureId;
		CCString custId;
		CCString dealId;

		newObject = ARMLOCAL_ParseObject(xmlResponse,typeId,date,filter,bookName,structureId,custId,dealId);

		if (newObject == NULL)
		{
			result.setMsg ("ARM_ERR: Object created from Summit is null");
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
			oldObject = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			thisClass = newObject->GetName();
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
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newObject)
			delete newObject;
		newObject = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newObject)
			delete newObject;
		newObject = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



//	--------------------------------------------------------------------------------------------
ARM_Object* 
obj_getTradeFromCalypso (const std::string & IdCalypso,
					  const ARM_Date& date,
					  const std::string & type,
					  const std::string& modelType)
{
	ARM_Object* trade = NULL;
	std::string xmlContent ;
	ARM_CalypsoToolkit::GetTrade(IdCalypso,date,xmlContent); 
//	trade = ARMLOCAL_ParseXMLForCalypsoObj(xmlContent,type,AsOf);
	trade =ARMLOCAL_ParseCalypsoObject(xmlContent,type,modelType,date);
	return trade;
}



long ARMLOCAL_GETINSTRUMENTFROMCALYPSO (CCString& idCalypso,
										CCString& typeId,
										CCString& modelType,
										double asOf,
										ARM_result& result,
										long objId)
{
	ARM_Object* newObject = NULL;
	ARM_Object* oldObject = NULL;

	long Id;
	ARM_CLASS_NAME thisClass;

	char sDate[11];

	CCString msg (" ");

	try
	{
		//ARM_Date date(MINDATE);
        ARM_Date date;

		if (asOf != -1)
		{
			Local_XLDATE2ARMDATE(asOf,sDate);
			ARM_Date tmpdate(sDate);
			date = tmpdate;
		}

		newObject = obj_getTradeFromCalypso((const char*)idCalypso,date,(const char*)typeId,(const char*) modelType);

		if (newObject == NULL)
		{
			result.setMsg ("ARM_ERR: Object created from Calypso is null");
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
			oldObject = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			thisClass = newObject->GetName();
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
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newObject)
			delete newObject;
		newObject = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newObject)
			delete newObject;
		newObject = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_GetMeanRevFromSummit (const CCString& C_ccy,
										const CCString& C_index,
										const CCString& C_cvname,
										double C_date,
										const CCString& C_NumFactor,
										ARM_result& result)
{
	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		char tmpDate[11];

		Local_XLDATE2ARMDATE(C_date,tmpDate);

		double res = etoolkit_getXMLMeanRevFromSummit(C_ccy,
													  C_index,
													  C_cvname,
													  (ARM_Date) tmpDate,
													  C_NumFactor);

		if (res == -10000.)
		{
			result.setMsg("ARM_ERR: error in getting Mean Reversion");
			return ARM_KO;
		}

		result.setDouble(res);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_GetCutOffFromSummit (const CCString& C_ccy,
									   const CCString& C_index,
									   const CCString& C_cvname,
									   const CCString& C_NumFactor,
									   double C_date,
									   ARM_result& result)
{
	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		char tmpDate[11];

		Local_XLDATE2ARMDATE(C_date,tmpDate);

		double res = etoolkit_getXMLCutOffFromSummit(C_ccy,
													 C_index,
													 C_cvname,
													 C_NumFactor,
													 (ARM_Date) tmpDate);

		if (res == -10000.)
		{
			result.setMsg("ARM_ERR: error in getting CutOff");
			return ARM_KO;
		}

		result.setDouble(res);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_HEDGE (long secId,
					 long hedgeId,
					 ARM_result& result)
{
	ARM_Matrix* res = NULL;

	ARM_PowerReverse* Prcs = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		Prcs = (ARM_PowerReverse*) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (PRCS for the moment");
			return ARM_KO;
		}

		res = Prcs->Hedge(hedgeId);

		if (res != NULL)
		{
			result.setLong(res->GetNumLines());
			result.setDouble((double) (res->GetNumCols()));

			for (int i=0;i<res->GetNumLines();i++)
				for (int j=0;j<res->GetNumCols();j++)
					result.setArray(res->Elt(i,j),i+j*res->GetNumLines());

			delete res;
		}
		else
		{
			result.setMsg ("ARM_ERR: Hedge Matrix is NULL");
			return ARM_KO;
		}

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_GetFxSmileFromSummit (const CCString& C_ccy1,
										const CCString& C_ccy2,
										const CCString& C_cvname,
										double C_date,
										ARM_result& result,
										long objId)
{
	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	long volId;

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(C_date,sDate);

	try
	{
		newVolCrv = etoolkit_GetXMLFxSmileFromSummit(C_ccy1,
													 C_ccy2,
													 C_cvname,
													 (ARM_Date) sDate);

		if (newVolCrv == NULL)
		{
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv);

			if (volId == RET_KO)
			{
				if (newVolCrv)
					delete newVolCrv;
				newVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			vc = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK((ARM_Object*)vc, ARM_VOL_LIN_INTERPOL) == 1)
			{
				if (vc)
				{
					delete vc;
					vc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (newVolCrv)
					delete newVolCrv;
				newVolCrv = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newVolCrv)
			delete newVolCrv;
		newVolCrv = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newVolCrv)
			delete newVolCrv;
		newVolCrv = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_GETINFOFROMPRCS(long prcsId,
							  const CCString& datatype,
							  ARM_result& result)
{
	ARM_PowerReverse* Prcs = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		Prcs = (ARM_PowerReverse*) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Prcs, ARM_POWERREVERSE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Security is not of a good type (PRCS for the moment");
			
           return(ARM_KO);
		}

		if ( strcmp((const char*) datatype,"FXUNDCCY") == 0 )
		   result.setString(Prcs->GetItsFxUnderLeg()->GetCurrencyUnit()->GetCcyName());
		else if ( strcmp((const char*) datatype,"FXNUMCCY") == 0 )
		   result.setString(Prcs->GetItsFxNumLeg()->GetCurrencyUnit()->GetCcyName());
		else if	( strcmp((const char*) datatype,"FUNDINGCCY") == 0 )
		   result.setString(Prcs->GetItsRealFundLeg()->GetCurrencyUnit()->GetCcyName());
		else if ( strcmp((const char*) datatype,"FXNUMNOT") == 0 )
		{
		   ARM_ReferenceValue* amount = Prcs->GetItsFxNumLeg()->GetAmount();
			
           result.setDouble(amount->GetDiscreteValues()->Elt(0));
		}
		else if ( strcmp((const char*) datatype,"FUNDINGNOT") == 0 )
		{
			ARM_ReferenceValue* amount = Prcs->GetItsRealFundLeg()->GetAmount();

			result.setDouble(amount->GetDiscreteValues()->Elt(0));
		}
        else if (( strcmp((const char *) datatype, "FUNDINGPV") == 0 )
                 ||
                 ( strcmp((const char*) datatype,"FUNDPV") == 0 )
                )
        {
            double fundPV = Prcs->GetRealFundingLegPV();

            result.setDouble(fundPV);
        }
		else
		{
			result.setMsg("Wrong data type: Valid are FXUNDCCY,FXNUMCCY,FUNDINGCCY,FXNUMNOT,FUNDINGNOT,FUNDPV or FUNDINGPV");
			return ARM_KO;
		}

		return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_GETOBJINFOFROMPRCS (long prcsId,
								  const CCString& datatype,
								  ARM_result& result,
								  long objId)
{
	ARM_ReferenceValue* newRefVal = NULL;
	ARM_ReferenceValue* oldRefVal = NULL;

	ARM_SwapLeg* newSwapleg = NULL;
	ARM_SwapLeg* oldSwapleg = NULL;

	ARM_Object* newObject = NULL;
	ARM_Object* oldObject = NULL;

	ARM_PowerReverse* Prcs = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	long newobjId;

	try
	{
		Prcs = (ARM_PowerReverse*) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (PRCS for the moment");
			return ARM_KO;
		}

		if (strcmp((const char*) datatype,"FXNUMCPN") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(((ARM_FixLeg*)Prcs->GetItsFxNumLeg())->GetVarCoupons())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"FXUNDCPN") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(((ARM_FixLeg*)Prcs->GetItsFxUnderLeg())->GetVarCoupons())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"FX0") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(Prcs->GetItsFxVariable())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"CAP") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(Prcs->GetItsCap())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"FLOOR") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(Prcs->GetItsFloor())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"FUNDSPREAD") == 0)
		{
			newRefVal = (ARM_ReferenceValue*)(Prcs->GetItsRealFundLeg()->GetVariableSpread())->Clone();
			newObject = newRefVal;
		}
		else if (strcmp((const char*) datatype,"FXNUM") == 0)
		{
			newSwapleg = (ARM_SwapLeg*)(Prcs->GetItsFxNumLeg())->Clone();
			newObject = newSwapleg;
		}
		else if (strcmp((const char*) datatype,"FIXINIT") == 0)
		{
			newSwapleg = (ARM_SwapLeg*)(Prcs->GetInitFixedLeg())->Clone();
			newObject = newSwapleg;
		}
		else if (strcmp((const char*) datatype,"FXUND") == 0)
		{
			newSwapleg = (ARM_FixLeg*)(Prcs->GetItsFxUnderLeg())->Clone();
			newObject = newSwapleg;
		}
		else if (strcmp((const char*) datatype,"FUNDING") == 0)
		{
			newSwapleg = (ARM_FixLeg*)(Prcs->GetItsRealFundLeg())->Clone();
			newObject = newSwapleg;
		}
		else
		{
			result.setMsg("Wrong data type: valid are FXNUMCPN,FXUNDCPN,FX0,CAP,FLOOR,FUNDSPREAD,FXNUM,FIXINIT,FXUND or FUNDING");
			return ARM_KO;
		}

		if (newObject == NULL)
		{
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newobjId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObject);

			if (newobjId == RET_KO)
			{
				if (newObject)
					delete newObject;
				newObject = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newobjId);

			return ARM_OK;
		}
		else
		{
			if ( (strcmp(datatype,"FXNUMCPN") == 0)
				|| (strcmp(datatype,"FXUNDCPN") == 0)
				|| (strcmp(datatype,"FX0") == 0)
				|| (strcmp(datatype,"CAP") == 0)
				|| (strcmp(datatype,"FLOOR") == 0)
				|| (strcmp(datatype,"FUNDSPREAD") == 0)
				)
			{
				oldObject = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			}
			else
			{
				oldObject = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			}

/*			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK((ARM_Object*)oldObject, ARM_REFERENCE_VALUE) == 1)
			{
*/				if (oldObject)
				{
					delete oldObject;
					oldObject = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newObject, objId);

				return ARM_OK;
/*			}
			else
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
*/		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newObject)
			delete newObject;
		newObject = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newObject)
			delete newObject;
		newObject = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




long ARMLOCAL_INITCRF(long crfId,
					  long zcCpnId,
					  long swVolId,
					  long capVolId,
					  long rhoCapId,
					  long nuCapId,
					  long betaCapId,
					  long rhoSwaptId,
					  long nuSwaptId,
					  long betaSwaptId,
					  long zcFundId,
					  long zcCpnBasisId,
					  long zcFundBasisId,
					  double fxSpot,
					  long isUpdate,
					  CCString& C_modelType,
					  long meanReversionId,
                      long skewRecalibFlag,
                      long inSABRSigmaOrAlpha,
					  ARM_result& result,
					  long objId)
{
	long cloneId;

	ARM_CRFCalculator* security  = NULL;
	ARM_CRFCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zcCpn       = NULL;
	ARM_ZeroCurve* zcFund      = NULL;
	ARM_ZeroCurve* zcCpnBasis  = NULL;
	ARM_ZeroCurve* zcFundBasis = NULL;

	ARM_VolCurve* swVol  = NULL;
	ARM_VolCurve* capVol = NULL;

	ARM_VolLInterpol* rhoCap = NULL;
	ARM_VolLInterpol* nuCap  = NULL;
	ARM_VolLInterpol* betaCap = NULL;

	ARM_VolLInterpol* rhoSwapt = NULL;
	ARM_VolLInterpol* nuSwapt = NULL;
	ARM_VolLInterpol* betaSwapt = NULL;

	ARM_ReferenceValue*	meanReversion = NULL;

    int SABRSigmaOrAlpha = int(inSABRSigmaOrAlpha);


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	ARM_result C_result;

	try
	{
		security = (ARM_CRFCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(crfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CALLREVFLOATER) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (CRF)");
			return ARM_KO;
		}

		zcCpn = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcCpn, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		swVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: swVol is not of a good type");
			return ARM_KO;
		}

		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: capVol is not of a good type");
			return ARM_KO;
		}

		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		if (rhoSwaptId != ARM_NULL_OBJECT)
		{
			rhoSwapt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwaptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoSwapt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoSwapt is not of a good type");
				return ARM_KO;
			}
		}

		if (nuSwaptId != ARM_NULL_OBJECT)
		{
			nuSwapt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwaptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuSwapt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuSwapt is not of a good type");
				return ARM_KO;
			}
		}

		if (betaSwaptId != ARM_NULL_OBJECT)
		{
			betaSwapt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwaptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaSwapt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaSwapt is not of a good type");
				return ARM_KO;
			}
		}


		if (zcFundId != ARM_NULL_OBJECT)
		{
			zcFund = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcFundId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcFund, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcFund is not of a good type");
				return ARM_KO;
			}
		}

		if (zcCpnBasisId != ARM_NULL_OBJECT)
		{
			zcCpnBasis = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnBasisId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcCpnBasis, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcCpnBasis is not of a good type");
				return ARM_KO;
			}
		}

		if (zcFundBasisId != ARM_NULL_OBJECT)
		{
			zcFundBasis = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcFundBasisId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcFundBasis, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcFundBasis is not of a good type");
				return ARM_KO;
			}
		}

		if (meanReversionId != ARM_NULL_OBJECT)
		{
			meanReversion = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(meanReversionId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(meanReversion, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: meanReversion is not of a good type");
				return ARM_KO;
			}
		}

		sec_clone = (ARM_CRFCalculator *) security->Clone();
		sec_clone->SetModelType(CCSTringToSTLString(C_modelType));
		sec_clone->SetMRS(meanReversion);
		sec_clone->InitCRFForSummit(zcCpn, swVol, capVol, 
                                    rhoCap, nuCap, betaCap, 
                                    rhoSwapt, nuSwapt, betaSwapt, 
									zcFund, zcCpnBasis, zcFundBasis, 
                                    fxSpot, 
                                    isUpdate,
                                    skewRecalibFlag,
                                    SABRSigmaOrAlpha);	

		if ( sec_clone == NULL )
		{
		   return(ARM_KO);
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if ( cloneId == RET_KO )
			{
				if (sec_clone)
				   delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_CRFCalculator* oldSec = (ARM_CRFCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSec, ARM_CALLREVFLOATER) == 0)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg ("ARM_ERR: Previous object is not of a good type (CRF)");
				return ARM_KO;
			}

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



// Initialise maturity cap
long ARMLOCAL_INITMATCAP(long matcapId,
					   long zcId,
					   long capVolId,
					   long rhoCapId,
					   long nuCapId,
					   long betaCapId,
					   long nbpas,
                       long inSABRSigmaOrAlpha,
					   ARM_result& result,
					   long objId)
{
    long cloneId;

	ARM_MaturityCapCalculator* security  = NULL;
    ARM_MaturityCapCalculator* sec_clone = NULL;

	ARM_ZeroCurve*         zc = NULL;
	ARM_VolCurve*      capVol = NULL;
	ARM_VolLInterpol*  rhoCap = NULL;
	ARM_VolLInterpol*   nuCap = NULL;
	ARM_VolLInterpol* betaCap = NULL;

    int SABRSigmaOrAlpha = int(inSABRSigmaOrAlpha);


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_MaturityCapCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(matcapId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_MATURITYCAP) == 0 )
		{
		   result.setMsg("ARM_ERR: Security is not of a good type (MATURITY CAP)");
			
           return(ARM_KO);
		}

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zc is not of a good type");

            return(ARM_KO);
		}

		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg("ARM_ERR: capVol is not of a good type");

			return(ARM_KO);
		}

		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

        sec_clone = (ARM_MaturityCapCalculator *) security->Clone();

		sec_clone->SetNbIterations(nbpas);

		sec_clone->InitMaturityCapForSummit(zc, capVol, rhoCap, nuCap, betaCap, 
                                            SABRSigmaOrAlpha);

		if ( sec_clone == NULL )
		{
		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) sec_clone);

			if ( cloneId == RET_KO )
			{
				if (sec_clone)
				   delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");

				return(ARM_KO);
			}
			
			result.setLong(cloneId);

			return(ARM_OK);
		}
		else
		{
			ARM_MaturityCapCalculator* oldSec = (ARM_MaturityCapCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSec, ARM_MATURITYCAP) == 0 )
			{
				if (sec_clone)
				   delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Previous object is not of a good type (MaturityCap)");

				return(ARM_KO);
			}

			if (oldSec)
			{
			   delete oldSec;
			   oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) sec_clone, objId);

			result.setLong(objId);

			return(ARM_OK);
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_INITTARN(long tarnId,
					   long zcCpnId,
					   long swVolId,
					   long capVolId,
					   long rhoCapId,
					   long nuCapId,
					   long betaCapId,
					   long rhoSwoptId,
					   long nuSwoptId,
					   long betaSwoptId,
					   long zcFundId,
					   long zcCpnBasisId,
					   long zcFundBasisId,
					   double fxSpot,
                       CCString& inModelType,
                       double betaCorrel,
                       double hump,
                       long   SABRSigmaOrAlpha,
					   ARM_result& result,
					   long objId)
{
	long cloneId;

	ARM_TARNCalculator* security  = NULL;
	ARM_TARNCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zcCpn         = NULL;
	ARM_ZeroCurve* zcFund        = NULL;
	ARM_ZeroCurve* zcCpnBasis    = NULL;
	ARM_ZeroCurve* zcFundBasis   = NULL;

	ARM_VolCurve* swVol          = NULL;
	ARM_VolCurve* capVol         = NULL;

	ARM_VolLInterpol* rhoCap     = NULL;
	ARM_VolLInterpol* nuCap      = NULL;
	ARM_VolLInterpol* betaCap    = NULL;

	ARM_VolLInterpol* rhoSwopt   = NULL;
	ARM_VolLInterpol* nuSwopt    = NULL;
	ARM_VolLInterpol* betaSwopt  = NULL;


  
       
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_TARNCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(tarnId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_TARN) == 0 )
		{
		   result.setMsg("ARM_ERR: Security is not of a good type (RFTARN)");
		
           return(ARM_KO);
		}

		zcCpn = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcCpn, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		//Volatlity (ATM or Vol Cube).
		swVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: swVol is not of a good type");
			return ARM_KO;
		}

		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: capVol is not of a good type");
			return ARM_KO;
		}
		
		// SABR contributions
		//Short term
		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		//Long term
		if (rhoSwoptId != ARM_NULL_OBJECT)
		{
			rhoSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (nuSwoptId != ARM_NULL_OBJECT)
		{
			nuSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (betaSwoptId != ARM_NULL_OBJECT)
		{
			betaSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaSwopt is not of a good type");
				return ARM_KO;
			}
		}

		//Additional infos
		if (zcFundId != ARM_NULL_OBJECT)
		{
			zcFund = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcFundId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcFund, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcFund is not of a good type");
				return ARM_KO;
			}
		}

		if (zcCpnBasisId != ARM_NULL_OBJECT)
		{
			zcCpnBasis = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnBasisId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcCpnBasis, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcCpnBasis is not of a good type");
				return ARM_KO;
			}
		}

		if (zcFundBasisId != ARM_NULL_OBJECT)
		{
			zcFundBasis = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcFundBasisId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcFundBasis, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: zcFundBasis is not of a good type");
				return ARM_KO;
			}
		}

        ARM_ModelType modelType;	

        if ( inModelType == "SBGM" )
        {
           modelType = ARM_PricingModelType::SBGM;
        }
        else
        {
           modelType = ARM_PricingModelType::SFRM2F;
        }

		sec_clone = (ARM_TARNCalculator *) security->Clone();

		sec_clone->InitTARNForSummit(zcCpn, swVol, capVol, rhoCap, nuCap, betaCap, 
                                     rhoSwopt, 
                                     nuSwopt, 
                                     betaSwopt, zcFund, zcCpnBasis, zcFundBasis,
                                     fxSpot,
                                     0, // hedgeUpdate
                                     modelType,
                                     betaCorrel,
                                     hump,
                                     int(SABRSigmaOrAlpha));

		if ( sec_clone == NULL )
		{
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_TARNCalculator* oldSec = (ARM_TARNCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			/*if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSec, ARM_TARN) == 0)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg ("ARM_ERR: Previous object is not of a good type (RFTARN)");
				return ARM_KO;
			}*/

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_INITCSB(long csbId,
					  long zcCpnId,
					  long swVolId,
					  long capVolId,
					  long rhoCapId,
					  long nuCapId,
					  long betaCapId,
					  long rhoSwoptId,
					  long nuSwoptId,
					  long betaSwoptId,
					  double hump,
					  double betaCorrel,
					  double reCorrel,
					  long inSABRSigmaOrAlpha,
					  ARM_result& result,
					  long objId)
{
	long cloneId;

	ARM_CallableSnowBallCalculator* security  = NULL;
	ARM_CallableSnowBallCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zcCpn         = NULL;

	ARM_VolCurve* swVol          = NULL;
	ARM_VolCurve* capVol         = NULL;

	ARM_VolLInterpol* rhoCap     = NULL;
	ARM_VolLInterpol* nuCap      = NULL;
	ARM_VolLInterpol* betaCap    = NULL;

	ARM_VolLInterpol* rhoSwopt   = NULL;
	ARM_VolLInterpol* nuSwopt    = NULL;
	ARM_VolLInterpol* betaSwopt  = NULL;


  
       
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_CallableSnowBallCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(csbId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CALLABLE_SNOWBALL) == 0 )
		{
		   result.setMsg("ARM_ERR: Security is not of a good type (CSB)");
		
           return(ARM_KO);
		}

		zcCpn = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zcCpn, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		//Volatlity (ATM or Vol Cube).
		swVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: swVol is not of a good type");
			return ARM_KO;
		}

		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: capVol is not of a good type");
			return ARM_KO;
		}
		
		// SABR contributions
		//Short term
		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		//Long term
		if (rhoSwoptId != ARM_NULL_OBJECT)
		{
			rhoSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (nuSwoptId != ARM_NULL_OBJECT)
		{
			nuSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (betaSwoptId != ARM_NULL_OBJECT)
		{
			betaSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaSwopt is not of a good type");
				return ARM_KO;
			}
		}

		sec_clone = (ARM_CallableSnowBallCalculator*) security->Clone();

		sec_clone->InitCSBFromSummit(zcCpn, swVol, capVol, rhoCap, nuCap, betaCap, 
                                     rhoSwopt, 
                                     nuSwopt, 
                                     betaSwopt,
                                     hump,
									 betaCorrel,
                                     reCorrel,
									 inSABRSigmaOrAlpha);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_CallableSnowBallCalculator* oldSec = (ARM_CallableSnowBallCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			/*if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSec, ARM_TARN) == 0)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg ("ARM_ERR: Previous object is not of a good type (RFTARN)");
				return ARM_KO;
			}*/

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_INITCAPTION(long captionId,
						  long mktDataManagerId,
						  vector<string> mdmKeys,
						  int nbFactors,
						  int SFRMVolType,
						  string SwoptCalibrationMode,
						  string BetaCalibrationMode,
						  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
						  ARM_result& result,
						  long objId)
{
	long cloneId;

	ARM_CaptionCalculator* security  = NULL;
	ARM_CaptionCalculator* sec_clone = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	try
	{
		security = (ARM_CaptionCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(captionId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CAPTION) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (CAPTION)");
			return ARM_KO;
		}

		//Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}
		
		ARM_CaptionCalculator::CalibrationMode SWOPTCalibMode = (ARM_CaptionCalculator::CalibrationMode) ARM::ARM_ArgConv_CaptionCalibMode.GetNumber(SwoptCalibrationMode); 
		ARM_CaptionCalculator::CalibrationMode BETACalibMode = (ARM_CaptionCalculator::CalibrationMode) ARM::ARM_ArgConv_CaptionCalibMode.GetNumber(BetaCalibrationMode);

		sec_clone = (ARM_CaptionCalculator *) security->Clone();

		sec_clone->InitCaptionForSummit(mdmKeys,
										mktDataManager, 
										nbFactors,
										SFRMVolType,
										SWOPTCalibMode,
										BETACalibMode,
										productsToPrice);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_CaptionCalculator* oldSec = (ARM_CaptionCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_INITCSO(long csoId,
					  long zcId,
					  long capVolId,
					  long swoptVolId,
					  long inSABRSigmaOrAlpha,
					  long rhoCapId,
					  long nuCapId,
					  long betaCapId,
					  long rhoSwoptId,
					  long nuSwoptId,
					  long betaSwoptId,
					  long flatVolId,
					  long convAdjustVolId,
					  long convAdjustManagerId,
					  long correlDiagCapId,
					  long correlDiagSwoptId,
					  long correlCorrId,
					  long mrsId,
					  long correlId,
					  long volRatioId,
					  long mrsSpreadId,
					  vector<double> modelParams,
					  vector<string> calibParams,
					  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
					  long forexId,
					  long fundZcId,
					  long domBasisZcId,
					  long fundBasisZcId,
					  ARM_result& result,
					  long objId)
{
	long cloneId;

	ARM_LocalCSOCalculator* security  = NULL;
	ARM_LocalCSOCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zc				= NULL;
	ARM_VolCurve* capVol			= NULL;
	ARM_VolCurve* swoptVol			= NULL;
	ARM_VolCurve* flatVol			= NULL;

	ARM_VolLInterpol* rhoCap		= NULL;
	ARM_VolLInterpol* nuCap			= NULL;
	ARM_VolLInterpol* betaCap		= NULL;
	ARM_VolLInterpol* rhoSwopt		= NULL;
	ARM_VolLInterpol* nuSwopt		= NULL;
	ARM_VolLInterpol* betaSwopt		= NULL;

	ARM_VolCurve* convAdjustVol		= NULL;
	ARM_ConvAdjustManager* convAdjustManager = NULL;

	ARM_VolCube* correlDiagCap		= NULL;
	ARM_VolCube* correlDiagSwopt	= NULL;
	ARM_CorrelManager* correlCorr	= NULL;

	ARM_CurveModelParam* mrs		= NULL;
	ARM_CurveModelParam* correl		= NULL;
	ARM_CurveModelParam* volRatio	= NULL;
	ARM_CurveModelParam* mrsSpread	= NULL;
       
	ARM_Forex* forex				= NULL;
	ARM_ZeroCurve* fundZc			= NULL;
	ARM_ZeroCurve* domBasisZc		= NULL;
	ARM_ZeroCurve* fundBasisZc		= NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_LocalCSOCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(csoId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CALLABLE_SPREADOPTION) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (callable spread option)");
			return ARM_KO;
		}

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		//Volatlity (ATM or Vol Cube).
		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: cap volatility is not of a good type");
			return ARM_KO;
		}
		
		swoptVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swoptVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: swaption volatility is not of a good type");
			return ARM_KO;
		}

		if (flatVolId != ARM_NULL_OBJECT)
		{
			flatVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(flatVolId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(flatVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: flat volatility is not of a good type");
				return ARM_KO;
			}
		}

		// SABR contributions
		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		if (rhoSwoptId != ARM_NULL_OBJECT)
		{
			rhoSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (nuSwoptId != ARM_NULL_OBJECT)
		{
			nuSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (betaSwoptId != ARM_NULL_OBJECT)
		{
			betaSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (convAdjustVolId != ARM_NULL_OBJECT)
		{
			convAdjustVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(convAdjustVolId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(convAdjustVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: convexity adjustment volatility is not of a good type");
				return ARM_KO;
			}
		}

		if (convAdjustManagerId != ARM_NULL_OBJECT)
		{
			convAdjustManager = (ARM_ConvAdjustManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(convAdjustManagerId);
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjustManager, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjustManager, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjustManager, ARM_MAPCONVADJUST) == 0))
			{
				result.setMsg ("ARM_ERR: convexity adjustment manager is not of a good type");
				return ARM_KO;
			}
		}

		//Correlation Cubes
		correlDiagCap = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlDiagCapId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlDiagCap, ARM_VOL_CUBE) == 0)
		{
			result.setMsg ("ARM_ERR: correlation cap is not of a good type");
			return ARM_KO;
		}
		
		if (correlDiagSwoptId != ARM_NULL_OBJECT)
		{
			correlDiagSwopt = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlDiagSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlDiagSwopt, ARM_VOL_CUBE) == 0)
			{
				result.setMsg ("ARM_ERR: correlation swaption is not of a good type");
				return ARM_KO;
			}
		}

		if (correlCorrId != ARM_NULL_OBJECT)
		{
			correlCorr = (ARM_CorrelManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlCorrId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlCorr, ARM_CORRELMANAGER) == 0)
			{
				result.setMsg ("ARM_ERR: correl corr is not of a good type");
				return ARM_KO;
			}
		}

		//Model parameters
		mrs = (ARM_CurveModelParam *) LOCAL_PERSISTENT_OBJECTS->GetObject(mrsId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mrs, ARM_MODELPARAM) == 0)
		{
			result.setMsg ("ARM_ERR: mrs is not of a good type");
			return ARM_KO;
		}

		if (correlId != ARM_NULL_OBJECT)
		{
			correl = (ARM_CurveModelParam *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correl, ARM_MODELPARAM) == 0)
			{
				result.setMsg ("ARM_ERR: correl is not of a good type");
				return ARM_KO;
			}
		}

		if (volRatioId != ARM_NULL_OBJECT)
		{
			volRatio = (ARM_CurveModelParam *) LOCAL_PERSISTENT_OBJECTS->GetObject(volRatioId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volRatio, ARM_MODELPARAM) == 0)
			{
				result.setMsg ("ARM_ERR: volRatio is not of a good type");
				return ARM_KO;
			}
		}

		if (mrsSpreadId != ARM_NULL_OBJECT)
		{
			mrsSpread = (ARM_CurveModelParam *) LOCAL_PERSISTENT_OBJECTS->GetObject(mrsSpreadId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mrsSpread, ARM_MODELPARAM) == 0)
			{
				result.setMsg ("ARM_ERR: mrs spread is not of a good type");
				return ARM_KO;
			}
		}

		// calibration and model parameters
		ARM_LocalCSOCalculator::CalibrationType calibType;
		ARM_ModelType modelType;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;

		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		if (calibParams.size()>4) //switch to new calib hw2f
		{
			if (calibParams.size()==6)
			{
				//MODEL
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[0]);

				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (calibParams[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_LONG;
					else if (calibParams[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_SHORT;
					else if (calibParams[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::SHORT_LONG_SPREAD;
					else if (calibParams[1] == "DIAG")
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (calibParams[1] == "BASKET")
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (calibParams[1] == "DIAG_BASKET")
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (calibParams[1] == "DIAG") 
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (calibParams[1] == "BASKET") 
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (calibParams[1] == "DIAG_BASKET") 
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}

				if (calibParams[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (calibParams[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);
				else if (modelType==ARM_PricingModelType::HWM1F &&
						 calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 calibParams[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (calibParams[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (calibParams[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ZERO);
				else if (calibParams[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::CAP);
				else if (calibParams[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (calibParams[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}
			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[5]);
		}
		else
		{
			if (calibParams[0] == "DIAG") 
				calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;

			else if (calibParams[0] == "BASKET") 
				calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;

			else if (calibParams[0] == "DIAG_BASKET") 
				calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;

			else
			{
				result.setMsg ("ARM_ERR: invalid calib type");
				return ARM_KO;
			}

			if (calibParams.size()>1)
			{
				if (calibParams[1] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);

				else if (calibParams[1] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);

				else if (calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 calibParams[1] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (calibParams.size()>2)
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[2]);

			if (calibParams.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[3]);
		}

		if (forexId != ARM_NULL_OBJECT)
		{
			forex = (ARM_Forex *) LOCAL_PERSISTENT_OBJECTS->GetObject(forexId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(forex, ARM_FOREX) == 0 )
			{
				result.setMsg ("ARM_ERR: forex is not of a good type");
				return ARM_KO;
			}
		}

		if (fundZcId != ARM_NULL_OBJECT)
		{
			fundZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fundZcId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fundZc, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: funding ZC is not of a good type");
				return ARM_KO;
			}
		}

		if (domBasisZcId != ARM_NULL_OBJECT)
		{
			domBasisZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(domBasisZcId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(domBasisZc, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: domestic basis ZC is not of a good type");
				return ARM_KO;
			}
		}

		if (fundBasisZcId != ARM_NULL_OBJECT)
		{
			fundBasisZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fundBasisZcId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fundBasisZc, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: funding basis ZC is not of a good type");
				return ARM_KO;
			}
		}

		sec_clone = (ARM_LocalCSOCalculator *) security->Clone();

		sec_clone->InitCSOFromSummit(zc, capVol, swoptVol, inSABRSigmaOrAlpha,
									 rhoCap, nuCap, betaCap, rhoSwopt, nuSwopt, betaSwopt, 
									 flatVol, convAdjustVol, convAdjustManager,
									 correlDiagCap, correlDiagSwopt, correlCorr,
                                     mrs, correl, volRatio, mrsSpread,
									 productsToPrice,
									 modelParams,
									 calibType,
									 calibStrikeTypeVec,
									 modelType,
									 vnsMethod,
									 forex,
									 fundZc, 
									 domBasisZc,
									 fundBasisZc);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("Init CSO: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_LocalCSOCalculator* oldSec = (ARM_LocalCSOCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("InitCSO: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_INITCSO(long csoId,
					  long mktDataManagerId,
					  // vector<string> mdmKeys, : Now irrelevant
					  vector<double> modelParams,
					  vector<string> calibParams,
					  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
					  ARM_result& result,
					  long objId)
{
	long cloneId;

	ARM_LocalCSOCalculator* security  = NULL;
	ARM_LocalCSOCalculator* sec_clone = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	try
	{
		security = (ARM_LocalCSOCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(csoId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CALLABLE_SPREADOPTION) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (callable spread option)");
			return ARM_KO;
		}

		//Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		ARM_LocalCSOCalculator::CalibrationType calibType;
		ARM_ModelType modelType;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;

		/*ARM_LocalCSOCalculator::CalibStrikeType calibStrikeType;
		
		if (calibParams[0] == "DIAG") 
			calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
		else if (calibParams[0] == "BASKET") 
			calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
		else if (calibParams[0] == "DIAG_BASKET") 
			calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
		else
		{
			result.setMsg ("ARM_ERR: invalid calib type");
			return ARM_KO;
		}
		
		if (calibParams.size()>1)
		{
			if (calibParams[1] == "ATM") 
				calibStrikeType = ARM_LocalCSOCalculator::ATM;

			else if (calibParams[1] == "EQUIVALENT") 
				calibStrikeType = ARM_LocalCSOCalculator::EQUIVALENT;
			else
			{
				result.setMsg ("ARM_ERR: invalid calib strike type");
				return ARM_KO;
			}
		}

		if (calibParams.size()>2)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[2]);

		if (calibParams.size()>3)
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[3]);

		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(1,calibStrikeType);*/

				/*
		ARM_LocalCSOCalculator::CalibStrikeType calibStrikeType;
		
		 if (calibParams[0] == "DIAG") 
			calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
		else if (calibParams[0] == "BASKET") 
			calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
		else if (calibParams[0] == "DIAG_BASKET") 
			calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
		else
		{
			result.setMsg ("ARM_ERR: invalid calib type");
			return ARM_KO;
		}
		
		if (calibParams.size()>1)
		{
			if (calibParams[1] == "ATM") 
				calibStrikeType = ARM_LocalCSOCalculator::ATM;
			else if (calibParams[1] == "EQUIVALENT") 
				calibStrikeType = ARM_LocalCSOCalculator::EQUIVALENT;
			else
			{
				result.setMsg ("ARM_ERR: invalid calib strike type");
				return ARM_KO;
			}
		}

		if (calibParams.size()>2)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[2]);

		if (calibParams.size()>3)
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[3]);
		
		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(1,calibStrikeType);
		*/

		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		if (calibParams.size()>4) //switch to new calib hw2f
		{
			if (calibParams.size()==6)
			{
				//MODEL
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[0]);

				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (calibParams[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_LONG;
					else if (calibParams[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_SHORT;
					else if (calibParams[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::SHORT_LONG_SPREAD;
					else if (calibParams[1] == "DIAG")
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (calibParams[1] == "BASKET")
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (calibParams[1] == "DIAG_BASKET")
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (calibParams[1] == "DIAG") 
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (calibParams[1] == "BASKET") 
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (calibParams[1] == "DIAG_BASKET") 
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}

				if (calibParams[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (calibParams[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);
				else if (modelType==ARM_PricingModelType::HWM1F &&
						 calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 calibParams[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (calibParams[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (calibParams[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ZERO);
				else if (calibParams[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::CAP);
				else if (calibParams[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (calibParams[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}
					
			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[5]);
		}
		else
		{
			if (calibParams[0] == "DIAG") 
				calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;

			else if (calibParams[0] == "BASKET") 
				calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;

			else if (calibParams[0] == "DIAG_BASKET") 
				calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;

			else
			{
				result.setMsg ("ARM_ERR: invalid calib type");
				return ARM_KO;
			}

			if (calibParams.size()>1)
			{
				if (calibParams[1] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);

				else if (calibParams[1] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);

				else if (calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 calibParams[1] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (calibParams.size()>2)
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibParams[2]);

			if (calibParams.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(calibParams[3]);
		}


		sec_clone = (ARM_LocalCSOCalculator *) security->Clone();

		sec_clone->InitCSOFromSummit(mktDataManager, 
									 productsToPrice,
									 modelParams,
									 calibType,
									 calibStrikeTypeVec,
									 modelType,
									 vnsMethod);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_LocalCSOCalculator* oldSec = (ARM_LocalCSOCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_INITBERMUDASWAPTION (long bsId,
								   long mktDataManagerId,
								   vector<string> mdmKeys,
								   vector<int>* controlVariates,
								   vector<double>* controlPrices,
								   vector<string> modelParams,
								   bool mrsCalibFlag,
								   bool atmCalibFlag,
								   int numMethodType,
								   int amcIter,
								   int mcIter,
								   int maxBucketSize,
								   string genType,
								   int treeSteps,
								   vector<int> portfolioMode,
								   bool boundaryFlag,
								   bool approxMarginFlag,
								   bool freezeBetasFlag,
								   int modelType,
								   bool calculateProbaFlag,
								   ARM_result& result,
								   long objId)
{
	long cloneId;

	ARM_BermudaSwaptionCalculator* security  = NULL;
	ARM_BermudaSwaptionCalculator* sec_clone = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	try
	{
		security = (ARM_BermudaSwaptionCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(bsId);

		//Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		// model params : convert into Reference Values
		vector<ARM_ReferenceValue*>* modelParamsRef = new vector<ARM_ReferenceValue*>(7);
		for (int i=0; i<modelParams.size(); i++)
		{
			CCString objId(modelParams[i].c_str());
			if (LocalGetNumObjectId(objId) == ARM_KO)
			{
				double value = 0.0;
				if (modelParams[i] == "DEFAULT")
				{
					if (i == 5)
						value = -0.1;
					if (i == 6)
						value = 0.1;
				}
				else
					value = atof(modelParams[i].c_str());
		
				(*modelParamsRef)[i] = new ARM_ReferenceValue(value);
			}
			else
			{
				(*modelParamsRef)[i] = dynamic_cast<ARM_ReferenceValue*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(objId))->Clone());
			}
		}

		sec_clone = (ARM_BermudaSwaptionCalculator *) security->Clone();
		sec_clone->InitBermudaSwaptionForSummit( mdmKeys,
												 mktDataManager,
												 controlVariates,
												 controlPrices,
												 modelParamsRef,
												 mrsCalibFlag,
												 atmCalibFlag,
												 numMethodType,
												 amcIter,
												 mcIter,
												 maxBucketSize,
												 genType,
												 treeSteps,
												 portfolioMode,
												 boundaryFlag,
												 approxMarginFlag,
												 freezeBetasFlag,
												 modelType,
												 calculateProbaFlag);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_BermudaSwaptionCalculator* oldSec = (ARM_BermudaSwaptionCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_INITSWAPTIONBERMUDA ( long bsId,
								    vector<int>* controlVariates,
								    vector<double>* controlPrices,
								    vector<string> modelParams,
								    bool mrsCalibFlag,
								    bool atmCalibFlag,
								    int numMethodType,
								    int amcIter,
								    int mcIter,
								    int maxBucketSize,
								    string genType,
								    int treeSteps,
								    vector<int> portfolioMode,
								    bool boundaryFlag,
								    bool approxMarginFlag,
								    bool freezeBetasFlag,
								    int modelType,
								    bool calculateProbaFlag,
								    long zcCpnId,
									long swoptVcId,
									long capVcId,
									long capRoId,
									long capNuId,
									long capBetaId,
									long swoptRoId,
									long swoptNuId,
									long swoptBetaId,
									long normalModelId,
									long inSABRSigmaOrAlpha,
									ARM_result& result,
								    long objId)
{
	long cloneId;

	ARM_BermudaSwaptionCalculator* security  = NULL;
	ARM_BermudaSwaptionCalculator* sec_clone = NULL;

	ARM_ZeroCurve*		zcCpn		= NULL;
	ARM_VolCurve*		swoptVc		= NULL;
	ARM_VolCurve*		capVc		= NULL;
	ARM_VolLInterpol*	capRo		= NULL;
	ARM_VolLInterpol*	capNu		= NULL;
	ARM_VolLInterpol*	capBeta		= NULL;
	ARM_VolLInterpol*	swoptRo		= NULL;
	ARM_VolLInterpol*	swoptNu		= NULL;
	ARM_VolLInterpol*	swoptBeta	= NULL;
	ARM_MarketIRModel*	normalModel	= NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	try
	{
		security = dynamic_cast<ARM_BermudaSwaptionCalculator*> (LOCAL_PERSISTENT_OBJECTS->GetObject(bsId));

		//Zc Cpn Id
		if (zcCpnId != ARM_NULL_OBJECT)
		{
			zcCpn = dynamic_cast<ARM_ZeroCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(zcCpnId));

			if (!zcCpn)
			{
				result.setMsg ("ARM_ERR: ZcCpn is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: ZcCpn is not of a good type");
			return ARM_KO;
		}

		//Swopt Vol Curve Id
		if (swoptVcId != ARM_NULL_OBJECT)
		{
			swoptVc = dynamic_cast<ARM_VolCurve*> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptVcId));

			if (!swoptVc)
			{
				result.setMsg ("ARM_ERR: Swopt Vc is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: Swopt Vc is not of a good type");
			return ARM_KO;
		}

		//Cap Vol Curve Id
		if (capVcId != ARM_NULL_OBJECT)
		{
			capVc = dynamic_cast<ARM_VolCurve*> (LOCAL_PERSISTENT_OBJECTS->GetObject(capVcId));

			if (!capVc)
			{
				result.setMsg ("ARM_ERR: Cap Vc is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: Cap Vc is not of a good type");
			return ARM_KO;
		}

		//Cap Ro Curve Id
		if (capRoId != ARM_NULL_OBJECT)
		{
			capRo = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(capRoId));

			if (!capRo)
			{
				result.setMsg ("ARM_ERR: Cap Ro is not of a good type");
				return ARM_KO;
			}
		}
	
		//Cap Nu Curve Id
		if (capNuId != ARM_NULL_OBJECT)
		{
			capNu = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(capNuId));

			if (!capNu)
			{
				result.setMsg ("ARM_ERR: Cap Nu is not of a good type");
				return ARM_KO;
			}
		}
	
		//Cap Beta Curve Id
		if (capBetaId != ARM_NULL_OBJECT)
		{
			capBeta = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(capBetaId));

			if (!capBeta)
			{
				result.setMsg ("ARM_ERR: Cap Beta is not of a good type");
				return ARM_KO;
			}
		}

		//Swopt Ro Curve Id
		if (swoptRoId != ARM_NULL_OBJECT)
		{
			swoptRo = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptRoId));

			if (!swoptRo)
			{
				result.setMsg ("ARM_ERR: Swopt Ro is not of a good type");
				return ARM_KO;
			}
		}

		//Swopt Nu Curve Id
		if (swoptNuId != ARM_NULL_OBJECT)
		{
			swoptNu = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptNuId));

			if (!swoptNu)
			{
				result.setMsg ("ARM_ERR: Swopt Nu is not of a good type");
				return ARM_KO;
			}
		}

		//Swopt Beta Curve Id
		if (swoptBetaId != ARM_NULL_OBJECT)
		{
			swoptBeta = dynamic_cast<ARM_VolLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptBetaId));

			if (!swoptBeta)
			{
				result.setMsg ("ARM_ERR: Swopt Beta is not of a good type");
				return ARM_KO;
			}
		}
	
		//Normal model Id
		if (normalModelId != ARM_NULL_OBJECT)
		{
			normalModel = dynamic_cast<ARM_MarketIRModel*> (LOCAL_PERSISTENT_OBJECTS->GetObject(normalModelId));

			if (!normalModel)
			{
				result.setMsg ("ARM_ERR: Normal Model is not of a good type");
				return ARM_KO;
			}
		}

		// model params : convert into Reference Values
		vector<ARM_ReferenceValue*>* modelParamsRef = new vector<ARM_ReferenceValue*>(7);
		for (int i=0; i<modelParams.size(); i++)
		{
			CCString objId(modelParams[i].c_str());
			if (LocalGetNumObjectId(objId) == ARM_KO)
			{
				double value = 0.0;
				if (modelParams[i] == "DEFAULT")
				{
					if (i == 5)
						value = -0.1;
					if (i == 6)
						value = 0.1;
				}
				else
					value = atof(modelParams[i].c_str());
				
				(*modelParamsRef)[i] = new ARM_ReferenceValue(value);
			}
			else
			{
				(*modelParamsRef)[i] = dynamic_cast<ARM_ReferenceValue*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(objId))->Clone());
			}
		}

		sec_clone = (ARM_BermudaSwaptionCalculator *) security->Clone();

		sec_clone->InitBermudaSwaptionForSummit(controlVariates,
												controlPrices,
												modelParamsRef,
												mrsCalibFlag,
												atmCalibFlag,
												numMethodType,
												amcIter,
												mcIter,
												maxBucketSize,
												genType,
												treeSteps,
												portfolioMode,
												boundaryFlag,
												approxMarginFlag,
												freezeBetasFlag,
												modelType,
												calculateProbaFlag,
												zcCpn, 
												swoptVc, 
												capVc, 
												capRo, 
												capNu,
												capBeta,
												swoptRo,
												swoptNu,
												swoptBeta,
												normalModel,
												0, //hedge update
												inSABRSigmaOrAlpha);

		delete modelParamsRef;

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_BermudaSwaptionCalculator* oldSec = (ARM_BermudaSwaptionCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_GetFixing (const CCString& source,
							 const CCString& index,
							 const CCString& tenor,
							 const CCString& ccy,
							 double asof,
							 ARM_result& result)
{
	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		char tmpDate[11];

		Local_XLDATE2ARMDATE(asof,tmpDate);

		double res = etoolkit_getXMLFixingFromSummit(source,
													 index,
													 tenor,
													 ccy,
													 (ARM_Date) tmpDate);

		if (res == -10000.)
		{
			result.setMsg("ARM_ERR: error in getting Fixing");
			return ARM_KO;
		}

		result.setDouble(res);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_SetSecurityData (long secId,
							   const VECTOR<double>& data,
					           ARM_result& result,
					           long objId)
{
    CCString msg(""); // used in macro ARM_RESULT()

    ARM_Swaption *newOSW=NULL;
    ARM_Option *newOpt=NULL;
	ARM_CapFloor* newCF=NULL;
    ARM_Security *newSec=NULL;
	try
	{
		if (data.size()<1)
		{
			result.setMsg("ARM_ERR: Data to set are missing");
			return ARM_KO;
		}

		if((newOSW=dynamic_cast<ARM_Swaption *>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId))))
		{
			newOSW = static_cast<ARM_Swaption *>(newOSW->Clone());

			double strike=data[0];
			newOSW->UpdateStrike(strike);

			newSec = newOSW;
		}
		else if((newOpt=dynamic_cast<ARM_Option *>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId))))
		{
			newOpt = static_cast<ARM_Option *>(newOpt->Clone());

			double strike=data[0];
			newOpt->SetStrike(strike);

		   newSec = newOpt;
		}
		else if((newCF=dynamic_cast<ARM_CapFloor *>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId))))
		{
			/// Cap/floor or spread option
			newCF = static_cast<ARM_CapFloor *>(newCF->Clone());

			double strike=data[0];
			newCF->SetStrike(strike);

			if(data.size()>1)
				newCF->SetCapFloorType(data[1]>0 ? K_CAP : K_FLOOR);

		   newSec = newCF;
		}
		else
		{
			result.setMsg("ARM_ERR: SetData available only for ARM_SWAPTION, ARM_OPTION, ARM_CAPFLOOR or ARM_SPREADOPTION (strike, option type)");
			return ARM_KO;
		}


        /// Assign the new swaption in ARM cache
		if(!assignObject( newSec, result, objId ) )
        {
            delete newSec;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}
	catch(Exception& x)
	{
        delete newSec;

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_GetSecurityData (long secId,
							   const CCString& data,
					           ARM_result& result)
{
    CCString msg(""); // used in macro ARM_RESULT()

    ARM_Swaption *OSW;
    ARM_Option *Opt;
	try
	{
        string typestr = CCSTringToSTLString(data);
        string typeToGet(stringGetUpper(typestr));
		if (typeToGet != "STRIKE")
		{
			result.setMsg("ARM_ERR: Data is not valid");
			return ARM_KO;
		}

        double strike;
		if((OSW = dynamic_cast<ARM_Swaption *>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId))))
            strike = OSW->GetStrike();
        else if((Opt = dynamic_cast<ARM_Option *>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId))))
            strike = Opt->GetStrike();
		else
        {
			result.setMsg("ARM_ERR: GetData available only for ARM_SWAPTION or ARM_OPTION");
			return ARM_KO;
        }


        result.setDouble(strike);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_GetLastDateWarm(ARM_result& result)
{
	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		CCString res = ARMLOCAL_ParseXMLForGetDateLastWarm();

		result.setString(res);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_SecurityFlows(vector<CCString>& C_labels,
							vector<double>& C_values,
							ARM_result& result,
							long objId)
{
	int	nbCol = C_labels.size();
	int	nbLin = C_values.size() / nbCol;
	long securityFlowsId;

	vector<string>		labelList(nbCol);
	vector<ARM_Vector*>	dataList(nbCol);
	ARM_Vector*			dataVector = NULL;
	double*	data = new double[nbLin];
	
	ARM_SecurityFlows*	createdSecurityFlows = NULL;
	ARM_SecurityFlows*	prevSecurityFlows = NULL;

	if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	try
	{
		for(int i=0; i<nbCol; i++)
		{
			labelList.at(i) = C_labels.at(i);

			for(int j=0; j<nbLin; j++)
				data[j] = C_values.at(i + j*nbCol);

			dataVector = new ARM_Vector(nbLin, data);
			
			dataList.at(i) = dataVector;
		}

		delete [] data;

		createdSecurityFlows = new ARM_SecurityFlows(labelList, dataList);

		if( createdSecurityFlows == NULL )
		{
		   result.setMsg("ARM_ERR: createdSecurityFlows is null");

		   return	ARM_KO;
		}

		if( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			securityFlowsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurityFlows);

			if( securityFlowsId == RET_KO )
			{
				delete	createdSecurityFlows;

				result.setMsg("ARM_ERR: Pb with inserting object");				

				return	ARM_KO;
			}

			result.setLong(securityFlowsId);

			return	ARM_OK;
		}
		else
		{
			prevSecurityFlows = (ARM_SecurityFlows*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevSecurityFlows, ARM_SECURITY_FLOWS) == 1)
			{
				delete	prevSecurityFlows;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurityFlows, objId);

				return ARM_OK;
			}
			else
			{
				delete	createdSecurityFlows;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				
				return	ARM_KO;
			}
		}		
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}

long ARMLOCAL_INITCRASPREAD(  long ccsoId,
							  long zcId,
							  long capVolId,
							  long rhoCapId,
							  long nuCapId,
							  long betaCapId,
							  long swoptVolId,
							  long rhoSwoptId,
							  long nuSwoptId,
							  long betaSwoptId,
							  long sigmaOrAlpha,
							  long convAdjustVolCapId,
							  long convAdjustVolSwoptId,
							  long convAdjustType,
							  long correlCorrId,
							  long correlDiagCapId,
							  long correlDiagSwoptId,
							  double mrs,
							  double volRatio,
							  double mrsSpread,
							  double correl,
							  vector<string> modelParams,
							  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
							  const vector<string>& localCalibFlags,
							  ARM_result& result,
							  long objId)
{
	long cloneId;

	ARM_CRASpreadCalculator* security  = NULL;
	ARM_CRASpreadCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zc				= NULL;
	ARM_VolCurve* capVol			= NULL;
	ARM_VolLInterpol* rhoCap		= NULL;
	ARM_VolLInterpol* nuCap			= NULL;
	ARM_VolLInterpol* betaCap		= NULL;

	ARM_VolCurve* swoptVol			= NULL;
	ARM_VolLInterpol* rhoSwopt		= NULL;
	ARM_VolLInterpol* nuSwopt		= NULL;
	ARM_VolLInterpol* betaSwopt		= NULL;

	ARM_VolCurve* convAdjustVolCap		= NULL;
	ARM_VolCurve* convAdjustVolSwopt	= NULL;

	ARM_CorrelManager* correlCorr	= NULL;
	ARM_VolCurve* correlDiagCap		= NULL;
	ARM_VolCurve* correlDiagSwopt	= NULL;
       
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_CRASpreadCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccsoId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CALLABLE_CORRIDOR_SPREADOPTION) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (callable spread option)");
			return ARM_KO;
		}

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		//Cap Volatility (ATM or Vol Cube).
		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: cap volatility is not of a good type");
			return ARM_KO;
		}
		
		if (rhoCapId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuCapId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaCapId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		//Swaption Volatility (ATM or Vol Cube).
		swoptVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swoptVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: swaption volatility is not of a good type");
			return ARM_KO;
		}
		
		if (rhoSwoptId != ARM_NULL_OBJECT)
		{
			rhoSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (nuSwoptId != ARM_NULL_OBJECT)
		{
			nuSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (betaSwoptId != ARM_NULL_OBJECT)
		{
			betaSwopt = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaSwopt, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}

		if (convAdjustVolCapId != ARM_NULL_OBJECT)
		{
			convAdjustVolCap = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(convAdjustVolCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(convAdjustVolCap, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: convexity adjustment volatility for cap is not of a good type");
				return ARM_KO;
			}
		}

		if (convAdjustVolSwoptId != ARM_NULL_OBJECT)
		{
			convAdjustVolSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(convAdjustVolSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(convAdjustVolSwopt, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: convexity adjustment volatility for swaption is not of a good type");
				return ARM_KO;
			}
		}

		//Correlation
		if (correlDiagCapId != ARM_NULL_OBJECT)
		{
			correlDiagCap = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlDiagCapId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlDiagCap, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: correlation cap is not of a good type");
				return ARM_KO;
			}
		}
		
		if (correlDiagSwoptId != ARM_NULL_OBJECT)
		{
			correlDiagSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlDiagSwoptId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlDiagSwopt, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: correlation swaption is not of a good type");
				return ARM_KO;
			}
		}

		if (correlCorrId != ARM_NULL_OBJECT)
		{
			correlCorr = (ARM_CorrelManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlCorrId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlCorr, ARM_CORRELMANAGER) == 0)
			{
				result.setMsg ("ARM_ERR: correl corr is not of a good type");
				return ARM_KO;
			}
		}

		// calibration and model parameters
		ARM_CRASpreadCalculator::CalibrationType calibType;
		ARM_ModelType modelType;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;
	
		vector<ARM_CRASpreadCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		if (modelParams.size()>4) //switch to new calib hw2f
		{
			if (modelParams.size()==6)
			{
				//MODEL
		if (modelParams[0] == "HW1F")
			modelType = ARM_PricingModelType::HWM1F;
		else if (modelParams[0] == "HW2F")
			modelType = ARM_PricingModelType::HWM2F;
		else
		{
			result.setMsg ("ARM_ERR: invalid model type");
			return ARM_KO;
		}

				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (modelParams[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_LONG;
					else if (modelParams[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_SHORT;
					else if (modelParams[1] == "DIAG_INDEX_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_INDEX;
					else if (modelParams[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::SHORT_LONG_SPREAD;
					else if (modelParams[1] == "DIAG")
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (modelParams[1] == "BASKET")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (modelParams[1] == "BASKET_SIMPLE")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (modelParams[1] == "DIAG") 
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (modelParams[1] == "BASKET") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (modelParams[1] == "BASKET_SIMPLE") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
		
				if (modelParams[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (modelParams[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);
				else if (modelParams[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (modelParams[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (modelParams[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ZERO);
				else if (modelParams[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::CAP);
				else if (modelParams[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (modelParams[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}

			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(modelParams[5]);
		}
		else
		{
			if (modelParams[0] == "HW1F")
				modelType = ARM_PricingModelType::HWM1F;
			else if (modelParams[0] == "HW2F")
				modelType = ARM_PricingModelType::HWM2F;
			else
			{
				result.setMsg ("ARM_ERR: invalid model type");
				return ARM_KO;
			}

		if (modelParams.size()>1)
		{
			if (modelParams[1] == "DIAG") 
				calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;

			else if (modelParams[1] == "BASKET") 
				calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;

			else if (modelParams[1] == "BASKET_SIMPLE") 
				calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;

				else if (modelParams[1] == "DIAG_LONG_SPREAD")
					calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_LONG;

				else if (modelParams[1] == "DIAG_SHORT_SPREAD")
					calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_SHORT;

				else if (modelParams[1] == "SHORT_LONG_SPREAD")
					calibType = ARM_CRASpreadCalculator::SHORT_LONG_SPREAD;

			else
			{
				result.setMsg ("ARM_ERR: invalid calib type");
				return ARM_KO;
			}
		}		

		if (modelParams.size()>2)
		{
			if (modelParams[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);

				else if (modelParams[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);

				else if (modelParams[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (modelParams.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(modelParams[3]);
		}

		sec_clone = (ARM_CRASpreadCalculator *) security->Clone();

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ localCalibSelectors(localCalibFlags.size(), false);
		for (size_t i = 0; i < localCalibFlags.size(); ++i)
		{
			localCalibSelectors[i] = ((localCalibFlags[i] == "Y") || (localCalibFlags[i] == "YES"));
		}

		sec_clone->InitCRASpreadFromSummit( zc,
											capVol,
											rhoCap,
											nuCap,
											betaCap,
											swoptVol,
											rhoSwopt,
											nuSwopt,
											betaSwopt,
											sigmaOrAlpha,
											convAdjustVolCap,
											convAdjustVolSwopt,
											convAdjustType,
											correlCorr,
											correlDiagCap,
											correlDiagSwopt,
											mrs,
											volRatio,
											mrsSpread,
											correl,
											modelType,
											calibType,
											calibStrikeTypeVec,
											vnsMethod,
											productsToPrice,
											localCalibSelectors);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("Init CSO: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_CRASpreadCalculator* oldSec = (ARM_CRASpreadCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("InitCSO: unrecognized failure");
		return ARM_KO;
	}
}

extern long ARMLOCAL_INITGLOBALCAP( long calculatorId,
									long zcId,
									long capVolId,
									long rhoId,
									long nuId,
									long betaId,
									double mrs,
									ARM_Vector* calibParams,
									int nbSteps,
									vector<string> randGenerator,
									int samplerType,
									std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
									ARM_result& result,
									long objId)
{
	
	long cloneId;

	ARM_GlobalCapCalculator* security  = NULL;
	ARM_GlobalCapCalculator* sec_clone = NULL;

	ARM_ZeroCurve* zc				= NULL;
	ARM_VolCurve* capVol			= NULL;

	ARM_VolLInterpol* rhoCap		= NULL;
	ARM_VolLInterpol* nuCap			= NULL;
	ARM_VolLInterpol* betaCap		= NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = (ARM_GlobalCapCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(calculatorId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_GLOBALCAP_CALCULATOR) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (global cap calculator)");
			return ARM_KO;
		}

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: zcCpn is not of a good type");
			return ARM_KO;
		}

		//Volatlity (ATM or Vol Cube).
		capVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: cap volatility is not of a good type");
			return ARM_KO;
		}
		
		// SABR contributions
		if (rhoId != ARM_NULL_OBJECT)
		{
			rhoCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(rhoCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: rhoCap is not of a good type");
				return ARM_KO;
			}
		}

		if (nuId != ARM_NULL_OBJECT)
		{
			nuCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nuCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: nuCap is not of a good type");
				return ARM_KO;
			}
		}

		if (betaId != ARM_NULL_OBJECT)
		{
			betaCap = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(betaCap, ARM_VOL_LIN_INTERPOL) == 0)
			{
				result.setMsg ("ARM_ERR: betaCap is not of a good type");
				return ARM_KO;
			}
		}
		
		// Global Cap Value
		//globalCapParams->Elt(2) /= 100;

		sec_clone = (ARM_GlobalCapCalculator *) security->Clone();
		
		sec_clone->InitGlobalCapFromSummit(zc, capVol, rhoCap, nuCap, betaCap, 
                                     mrs, calibParams, nbSteps, randGenerator,
									 samplerType, productsToPrice);

		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("Init Global cap: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_GlobalCapCalculator* oldSec = (ARM_GlobalCapCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("InitGlobalCap: unrecognized failure");
		return ARM_KO;
	}
	
}

extern long ARMLOCAL_INITTARNFX(long calculatorId,
								vector<CCString> MCParams,
								int factorNb,
								int timeStepNb,
								int spaceStepNb,
								double stdDevNb,
								int skipPDE,
								int rescalling,
								string modelType,
								int smileFlag,
								int mixCalib,
								int oneFactorFlag,
								string correlType,
								vector<CCString> zeroCurvesId,
								vector<CCString> basisCurvesId,
								vector<CCString> forexId,
								vector<CCString> ATMVolId, //for swopt BSGen
								vector<CCString> fxVolId, //for BS fx models
								vector<CCString> mixtureParamsId, //for mixture fx models
								vector<CCString> mrsParamsId,
								vector<CCString> QParamsId,
								long correlMatrixId,
								ARM_result& result,
								long objId)
{
	long cloneId;

	ARM_TARNFXCalculator* security  = NULL;
	ARM_TARNFXCalculator* sec_clone = NULL;

	vector<ARM_ZeroCurve*> zeroCurves(0);
	vector<ARM_ZeroCurve*> basisCurves(0);
	vector<ARM_Forex*> forex(0);
	vector<ARM_VolCurve*> ATMVol(0);
	vector<ARM_VolCurve*> fxVol(0);
	vector<ARM_ParamsMixture_Fx*> mixtureParams(0);
	vector<ARM_CurveModelParam*> mrsParams(0);
	vector<ARM_CurveModelParam*> QParams(0);
	ARM_GP_Matrix* correlMatrix = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
       return(ARM_KO);
	}

	CCString msg("");

	ARM_result C_result;

	try
	{
		security = dynamic_cast<ARM_TARNFXCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(calculatorId));
		if (!security)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type (global cap calculator)");
			return ARM_KO;
		}

		int i;
		
		for (i=0; i<zeroCurvesId.size(); i++)
		{
			ARM_ZeroCurve* zc = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(zeroCurvesId[i])));
			if ( !zc )
			{
				result.setMsg ("ARM_ERR: zc is not of a good type");
				return ARM_KO;
			}
			zeroCurves.push_back(zc);
		}

		for (i=0; i<basisCurvesId.size(); i++)
		{
			ARM_ZeroCurve* basis = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(basisCurvesId[i])));
			if ( !basis )
			{
				result.setMsg ("ARM_ERR: basis curve is not of a good type");
				return ARM_KO;
			}
			basisCurves.push_back(basis);
		}

		for (i=0; i<forexId.size(); i++)
		{
			ARM_Forex* fx = dynamic_cast<ARM_Forex*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(forexId[i])));
			if ( !fx )
			{
				result.setMsg ("ARM_ERR: forex is not of a good type");
				return ARM_KO;
			}
			forex.push_back(fx);
		}
		
		for (i=0; i<ATMVolId.size(); i++)
		{
			ARM_VolCurve* vol = dynamic_cast<ARM_VolCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(ATMVolId[i])));
			if ( !vol )
			{
				result.setMsg ("ARM_ERR: ATM vol is not of a good type");
				return ARM_KO;
			}
			ATMVol.push_back(vol);
		}
		
		for (i=0; i<fxVolId.size(); i++)
		{
			ARM_VolCurve* vol = dynamic_cast<ARM_VolCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(fxVolId[i])));
			if ( !vol )
			{
				result.setMsg ("ARM_ERR: FX vol is not of a good type");
				return ARM_KO;
			}
			fxVol.push_back(vol);
		}
		
		for (i=0; i<mixtureParamsId.size(); i++)
		{
			ARM_ParamsMixture_Fx* param = dynamic_cast<ARM_ParamsMixture_Fx*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(mixtureParamsId[i])));
			if ( !param )
			{
				result.setMsg ("ARM_ERR: Mxiture param is not of a good type");
				return ARM_KO;
			}
			mixtureParams.push_back(param);
		}
		
		for (i=0; i<mrsParamsId.size(); i++)
		{
			ARM_CurveModelParam* param = dynamic_cast<ARM_CurveModelParam*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(mrsParamsId[i])));
			if ( !param )
			{
				result.setMsg ("ARM_ERR: MRS param is not of a good type");
				return ARM_KO;
			}
			mrsParams.push_back(param);
		}

		for (i=0; i<QParamsId.size(); i++)
		{
			ARM_CurveModelParam* param = dynamic_cast<ARM_CurveModelParam*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(QParamsId[i])));
			if ( !param )
			{
				result.setMsg ("ARM_ERR: Q param is not of a good type");
				return ARM_KO;
			}
			QParams.push_back(param);
		}

		correlMatrix = dynamic_cast<ARM_GP_Matrix*>(LOCAL_PERSISTENT_OBJECTS->GetObject(correlMatrixId));
		if ( !correlMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of a good type");
			return ARM_KO;
		}

		sec_clone = (ARM_TARNFXCalculator *) security->Clone();
		
		sec_clone->Init(atoi(MCParams[0].c_str()),
						atoi(MCParams[1].c_str()), 
						ARM_ArgConv_BaseGenAlgoType.GetNumber(MCParams[2].c_str()),
						ARM_ArgConv_TransformAlgoType.GetNumber(MCParams[3].c_str()),
						ARM_ArgConv_BaseGenAlgoType.GetNumber(MCParams[4].c_str()),
						ARM_ArgConv_TransformAlgoType.GetNumber(MCParams[5].c_str()),
						atof(MCParams[6].c_str()), 
						atof(MCParams[7].c_str()), 
						factorNb, timeStepNb, spaceStepNb, stdDevNb, skipPDE, rescalling, 
						ARM_ArgConv_TARNFXModelType.GetNumber(modelType), 
						smileFlag, mixCalib, oneFactorFlag, 
						ARM_ArgConv_MMCorrelType.GetNumber(correlType),
						zeroCurves, basisCurves, forex, ATMVol, fxVol, mixtureParams, mrsParams, QParams, correlMatrix);


		if (sec_clone == NULL)
		{
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			cloneId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone);

			if (cloneId == RET_KO)
			{
				if (sec_clone)
					delete sec_clone;
				sec_clone = NULL;

				result.setMsg("Init TARN FX: Pb with inserting object");
				return ARM_KO;
			}
			
			result.setLong(cloneId);

			return ARM_OK;
		}
		else
		{
			ARM_TARNFXCalculator* oldSec = (ARM_TARNFXCalculator *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldSec)
			{
				delete oldSec;
				oldSec = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)sec_clone, objId);

			result.setLong(objId);

			return ARM_OK;
		}

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		ARM_RESULT();
	}
	catch (...)
	{
		if (sec_clone)
			delete sec_clone;
		sec_clone = NULL;

		result.setMsg ("InitTARNFX: unrecognized failure");
		return ARM_KO;
	}
}

extern long ARMLOCAL_ViewCell(
	const long&				ObjectId,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	   
	try
	{ 		
		ARM_RootObject* object =(ARM_RootObject *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjectId);

		if (object == NULL)
		{
			result.setMsg ("ARM_ERR: Unknown or Null object");
			return ARM_KO;
		}
		
		string txt = object->toString("","");
		result.setString(txt.c_str());
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}




extern long ARMLOCAL_MatrixVectorViewer(
	 const long&		ObjId,
	 long&				nbRows,
     long&		nbCols,
	 vector<double>&    values,
     ARM_result&		result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "MatrixVector Info" );

		ARM_Object* armObj		= LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
		ARM_GP_Matrix* mat	    = dynamic_cast<ARM_GP_Matrix*>(armObj);
		ARM_GP_Vector* vec	    = dynamic_cast<ARM_GP_Vector*>(armObj);

		if ( mat)
		{
			nbRows = mat->GetRowsNb();
			nbCols = mat->GetColsNb();
			values = mat->GetValues();
		}
		else if ( vec)
		{
			nbRows = vec->size();
			nbCols = 0;
			values = vec->GetValues();
		}
		else
		{
			result.setMsg ("ARM_ERR: local object must be a matrix or a vector ID");
			return ARM_KO;
		}
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_GetFixingsFromInstrument (const CCString& idSummit,
										ARM_result& result,
										long objId)
{
	ARM_ReferenceValue* newRefVal = NULL;
	ARM_ReferenceValue* oldRefVal = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_result C_result;

	long newobjId;

	try
	{
		switchToETK();

		CCString xmlResponse = etoolkit_getXMLObjectFromSummit(idSummit,"EXOTIC");
		
		newRefVal = ARMLOCAL_ParseFixingFromInstrument(xmlResponse);

		if (newRefVal == NULL)
		{
			result.setMsg ("ARM_ERR: no past fixing");				
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newobjId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal);

			if (newobjId == RET_KO)
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newobjId);

			return ARM_OK;
		}
		else
		{
			oldRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldRefVal)
			{
				delete oldRefVal;
				oldRefVal = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal, objId);

			return ARM_OK;
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newRefVal)
			delete newRefVal;
		newRefVal = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newRefVal)
			delete newRefVal;
		newRefVal = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}
