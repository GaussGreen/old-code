#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "fixingsched.h"

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_zccurve.h>

#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_eToolkit.h>

#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\ARM_local_init.h>
#include <ARM\libarm_local\ARM_local_wrapper.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
#include <ccy\currency.h>
#include <crv\zerocurv.h>
#include <crv\oldzccurve.h>
#include <crv\zerointimp.h>
#include <crv\zerointspreaded.h>
#include <crv\zeroflat.h>
#include <crv\zerocbdf.h>
#include <crv\zerovsk.h>
#include <crv\zerospl3.h>
#include <crv\zerospli.h>
#include <crv\zerosplsum.h>
#include <inst\bond.h>
#include <util\fromto.h>
#include <stdlib.h>
#include <string>

#include <gpinflation\infcurv.h>
#include <gpbase\gplinalgconvert.h>
#include <gpinflation\resetmanager.h>

#include <libCCTools++\CCString.h>



using ARM::ARM_InfCurv;
using ARM::To_ARM_Vector;
using ARM::CreateARMGPVectorFromVECTOR;

long ARMLOCAL_CreateZCSwapInt (double Date,
							   VECTOR<CCString>& matu,
							   VECTOR<double>& rate,
							   long MMVsFut,
							   long SwapVsFut,
							   long Raw,
							   long interp,
							   const CCString& Ccy,
							   long swapFrqId,
							   long fixDayCount,
							   ARM_result& result,
							   long objId)
{
	long curveId;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	ARM_ZeroLInterpol* prevZcLin = NULL;
	ARM_Currency* aCcy = NULL;
    ARM_Vector* mktData = NULL;
	int i=0;
	int real_size = matu.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* tmp = NULL;
	char myCurveDate[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(Date,myCurveDate);

		tmp = (char*) Ccy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			delete tmp;;
		}
		tmp = NULL;

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j],"%s", "X");

		for (int j = 0; j < real_size; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZcLin = new ARM_ZeroLInterpol((ARM_Date) myCurveDate, psMatu,
							                 mktData, MMVsFut, SwapVsFut, 
							                 Raw, interp, 0, aCcy, swapFrqId, fixDayCount);

		if (mktData)
		   delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZcLin == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin);

			if (curveId == RET_KO)
			{
				if (createdZcLin)
					delete createdZcLin;
				createdZcLin = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmp)
			delete tmp;
		tmp = NULL;

		if (createdZcLin)
			delete createdZcLin;
		createdZcLin = NULL;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_zcflat(const double zeroFlat,
					 double date,
                     const CCString& ccy,
					 ARM_result& result,
					 long objId)
{
	long curveId;

    ARM_Currency currency((const char *) ccy);

	ARM_ZeroFlat* createdZeroCurveFlat = NULL;
	ARM_ZeroFlat* prevZeroCurveFlat = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char* myCurveDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,myCurveDate);
	
		createdZeroCurveFlat = new ARM_ZeroFlat((ARM_Date) myCurveDate, 
                                                zeroFlat,
                                                &currency);

		if (myCurveDate)
		   delete [] myCurveDate;
		myCurveDate = NULL;

		if (createdZeroCurveFlat == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurveFlat);

			if (curveId == RET_KO)
			{
				if (createdZeroCurveFlat)
				   delete createdZeroCurveFlat;
				createdZeroCurveFlat = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZeroCurveFlat = (ARM_ZeroFlat *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZeroCurveFlat, ARM_ZERO_FLAT) == 1)
			{
				if (prevZeroCurveFlat)
				{
					delete prevZeroCurveFlat;
					prevZeroCurveFlat = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurveFlat, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg("ARM_ERR: previous object is not of a good type");
				
                return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdZeroCurveFlat)
			delete createdZeroCurveFlat;
		createdZeroCurveFlat = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_DiscountPrice(long idCurve,
							double matu,
							ARM_result& result)
{
	double dResult;
	ARM_ZeroCurve* zc = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = zc->DiscountPrice(matu);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_DiscountYield(long idCurve,
							double matu,
							long meth,
							ARM_result& result)
{
	double dResult;
	ARM_ZeroCurve* zc=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = zc->DiscountYield(matu, meth);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_ForwardPrice (long idCurve,
							double matu1,
							double matu2,
							ARM_result& result)
{
	double dResult;
	ARM_ZeroCurve* zc=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = zc->ForwardPrice(matu1, matu2);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_ForwardYield (long idCurve,
							double matu1,
							double matu2,
							long meth,
							long adjDaycountId,
							long decompFreqId,
							long daycountId,
							ARM_result& result)
{
	double dResult;
	ARM_ZeroCurve* zc=NULL;
	double pMatu1(matu1), pMatu2(matu2);

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		if ( daycountId != -1 )
		{
			if ( (matu1 < 30000) || (matu2 < 30000) )
			{
				ARM_Date d1 = zc->GetAsOfDate();
				ARM_Date d2 = zc->GetAsOfDate();

				d1.AddDays(matu1*365);
				d2.AddDays(matu2*365);
				
				dResult = zc->ForwardYieldWithDayCount(d1,
													   d2, 
													   -1,
													   daycountId);
			}
			else
			{
				char myDate1[11];
				char myDate2[11];

				Local_XLDATE2ARMDATE(matu1, myDate1);
				Local_XLDATE2ARMDATE(matu2, myDate2);

				// pMatu1 = CountYears(daycountId,zc->GetAsOfDate(),(ARM_Date) myDate1);

				// pMatu2 = CountYears(daycountId, zc->GetAsOfDate(), (ARM_Date) myDate2);

				dResult = zc->ForwardYieldWithDayCount(ARM_Date(myDate1),
													   ARM_Date(myDate2), 
													   -1,
													   daycountId);
			}
		}
        else
        {           
		   dResult = zc->ForwardYield(pMatu1, pMatu2, meth);
        }

		if (( adjDaycountId == K_YES )
            &&
            ( daycountId == -1 )
           )
		{
			const char* isoccy = zc->GetCurrencyUnit()->GetCcyName();

			if (   (strcmp(isoccy,"ZAR") != 0)
				&& (strcmp(isoccy,"PLN") != 0)
				&& (strcmp(isoccy,"BEF") != 0)
				&& (strcmp(isoccy,"GBP") != 0)
				&& (strcmp(isoccy,"PTE") != 0)
				&& (strcmp(isoccy,"AUD") != 0)
				)
			{
				dResult *= 360./365.;
			}
		}

		if (decompFreqId != -9999)
		{
			dResult = FromRateToRate(dResult, 1.0, K_COMP_ANNUAL, decompFreqId);
		}
		
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_zclint(const  VECTOR<double>& matu,
					 const  VECTOR<double>& rate,
					 long   meth,
					 double aDate,
					 const  CCString& sCcy,
					 long   interpMeth,
                     int    matuAreDoubles, // 1: YES, O: NO
                     ARM_CRV_TERMS& sTerms,
					 ARM_result& result,
					 long objId)
{
	long curveId;

	ARM_ZeroLInterpol* zc		= NULL;
	ARM_ZeroLInterpol* newZc	= NULL;
	ARM_Vector* vMatu			= NULL;
	ARM_Vector* vRates			= NULL;

	int real_size = matu.size ();

	if ( rate.size() != real_size )
	{
	   result.setMsg("ARM_ERR: rate and matu must have the same size");
		
       return(ARM_KO);
	}

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}
 
	char* myCurveDate = new char[11];

	double* pdRate = NULL;
	double* pdMatu = NULL;

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(aDate,myCurveDate);

		pdRate = new double[real_size];
		pdMatu = new double[real_size];

		for (int i = 0; i < real_size; i++)
		{
			pdRate[i] = rate[i];
			pdMatu[i] = matu[i];
		}

        vMatu = new ARM_Vector(real_size, pdMatu);
 
        vRates = new ARM_Vector(real_size, pdRate);

            char* strCcy = (char *) sCcy;

        ARM_Currency ccy(strCcy);

        delete strCcy;
		

        if ( matuAreDoubles == 1 ) // maturities are doubles
        {
           newZc = new ARM_ZeroLInterpol((ARM_Date) myCurveDate, vMatu, 
                                         vRates, meth, 0, interpMeth, &ccy);
        }
        else
        {
           newZc = new ARM_ZeroLInterpol((ARM_Date) myCurveDate, vRates, (ARM_CRV_TERMS) sTerms, 
                                         meth, 0, interpMeth, &ccy);
        }  

        if (vMatu)
           delete vMatu;
		vMatu = NULL;

        if (vRates)
           delete vRates;
		vRates = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;

		if ( newZc == NULL )
		{
		   result.setMsg("ARM_ERR: Curve is null");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc);

			if ( curveId == RET_KO )
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc, objId);

				return ARM_OK;
			}
			else
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();

		if (newZc)
			delete newZc;
		newZc = NULL;

        if (vMatu)
           delete vMatu;
		vMatu = NULL;

        if (vRates)
           delete vRates;
		vRates = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;
 
		ARM_RESULT();
    }
}

long ARMLOCAL_shiftzclint (double value,
                          long nbplot,
                          const VECTOR<double>& matu,
					      const VECTOR<double>& rate,
					      long meth,
					      double aDate,
					      const CCString& sCcy,
					      long interpMeth,
					      ARM_result& result,
					      long objId)
{
	long curveId;

	ARM_ZeroLInterpol* zc = NULL;

	ARM_ZeroLInterpol* newZc = NULL;
	ARM_Vector* vMatu        = NULL;
	ARM_Vector* vRates       = NULL;

	int real_size = matu.size ();

	if ( rate.size() != real_size )
	{
	   result.setMsg("ARM_ERR: rate and matu must have the same size");
		
       return(ARM_KO);
	}

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}
 
	char* myCurveDate = new char[11];

	double* pdRate = NULL;
	double* pdMatu = NULL;

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(aDate,myCurveDate);

		pdRate = new double[real_size];
		pdMatu = new double[real_size];

		for (int i = 0; i < real_size; i++)
		{
            if(i <nbplot)
                pdRate[i] = rate[i]+value;
            else
                pdRate[i] = rate[i];

			pdMatu[i] = matu[i];
		}

        vMatu = new ARM_Vector(real_size, pdMatu);
 
        vRates = new ARM_Vector(real_size, pdRate);

        newZc = new ARM_ZeroLInterpol((ARM_Date) myCurveDate, vMatu, vRates, meth, 0, interpMeth);
		
        char* strCcy = (char *) sCcy;

        ARM_Currency ccy(strCcy);

        delete strCcy;
		
        newZc->SetCurrencyUnit(&ccy);


        if (vMatu)
           delete vMatu;
		vMatu = NULL;

        if (vRates)
           delete vRates;
		vRates = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;

		if ( newZc == NULL )
		{
		   result.setMsg("ARM_ERR: Curve is null");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc);

			if ( curveId == RET_KO )
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc, objId);

				return ARM_OK;
			}
			else
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();

		if (newZc)
			delete newZc;
		newZc = NULL;

        if (vMatu)
           delete vMatu;
		vMatu = NULL;

        if (vRates)
           delete vRates;
		vRates = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;
 
		ARM_RESULT();
    }
}


long ARMLOCAL_zcspreaded (long zcSprId,
						  long zcInitId,
						  double date,
						  long MMFreq,
						  long SwapFreq,
						  bool ccyIsObject,
						  const CCString& ccyName,
						  ARM_result& result,
						  long objId)
{
	long curveId;
	ARM_ZeroLInterpol *zc = NULL;
	ARM_BasisCurve *newZc = NULL;
	ARM_ZeroCurve *ZCSpread = NULL, *ZCInit = NULL;
	ARM_Currency* ccy= NULL;
 

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return(ARM_KO);
	}

	char myCurveDate[20];

	CCString msg("");

    try
    {
		Local_XLDATE2ARMDATE(date, myCurveDate);

		ZCSpread = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcSprId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ZCSpread, ARM_ZERO_LIN_INTERPOL) == 0)
		{
		   result.setMsg ("ARM_ERR: ZCSpread is not of a good type");
			
           return(ARM_KO);
		}

		ZCInit = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcInitId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZCInit, ARM_ZERO_CURVE) == 0)
		{
		   result.setMsg("ARM_ERR: ZCInit is not of a good type");
			
           return(ARM_KO);
		}
 
		if ( ccyName == "DEFAULT" )
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}

        newZc = new ARM_BasisCurve(	(ARM_Date) myCurveDate, ZCSpread, ZCInit,
									MMFreq, SwapFreq, ccy);

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newZc == NULL)
		{
			result.setMsg ("ARM_ERR: ZcSpread is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc);

			if (curveId == RET_KO)
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = dynamic_cast<ARM_BasisCurve*> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_BASIS_CURVE) == 1)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc, objId);

				return ARM_OK;
			}
			else
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
		if (newZc)
			delete newZc;
		newZc = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_CreateZCSwapIntSmooth (double Date,
									 VECTOR<CCString>& matu,
									 VECTOR<double>& rate,
									 long MMVsFut,
									 long SwapVsFut,
									 long Raw,
									 long interp,
									 const CCString& Ccy,
									 double lambda,
									 long prec,
									 ARM_result& result,
									 long objId)
{
	long curveId;

    ARM_ZeroInterpolation* prevZcLin = NULL;
    ARM_ZeroInterpolation* createdZeroCurve = NULL;

    ARM_Vector* mktData=NULL;
 
	int real_size = matu.size ();
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char asOfDate[11];

	ARM_Currency* aCcy = NULL;

	char* tmp = NULL;
	ARM_CRV_TERMS psMatu;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(Date,asOfDate);

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "X");

		for (int j = 0; j < real_size; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}
 
		tmp = Ccy.GetStr();

		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			if (tmp)
				delete tmp;
		}

		tmp = NULL;

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZeroCurve = new ARM_ZeroInterpolation((ARM_Date) asOfDate, psMatu,
									mktData, MMVsFut, SwapVsFut, 
									Raw, interp, aCcy, lambda, prec);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroInterpolation *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_INTERPOLATION) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}

			else
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (mktData)
		   delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GetMaturitiesFromZC(ARM_result& result,
								  long objId)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	ARM_ZeroCurve* ZC = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

	if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZC, ARM_ZERO_CURVE) == 1)
	{
		int sizeMkt = ZC->GetMktData()->itsMktValue->GetSize();
		for(int i = 0; i < sizeMkt; i++)
			result.setStringInVect(ZC->GetMktData()->itsMktTerms[i]);

		return ARM_OK;
	}

	else
	{
		result.setMsg ("ARM_ERR: object is not of a ZeroCurve type");
		return ARM_KO;
	}


}


long ARMLOCAL_CreateZCSwapFutInt(double Date,
								 VECTOR<CCString>& matu,
								 VECTOR<double>& rate,
								 long MMVsFut,
								 long SwapVsFut,
								 long Raw,
								 long interp,
								 const CCString& Ccy,
								 ARM_result& result,
								 long objId)
{
	long curveId;

    ARM_ZeroLInterpol* prevZcLin = NULL;
    ARM_ZeroLInterpol* createdZeroCurve = NULL;
 
    ARM_Vector* mktData=NULL;

	int real_size = matu.size ();
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char asOfDate[11];
	ARM_Currency* aCcy = NULL;
	char* tmp = NULL;
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(Date,asOfDate);

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "XX");

		for(int j = 0; j < real_size; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		tmp = Ccy.GetStr();

		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			if (tmp)
				delete tmp;
		}

		tmp = NULL;

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZeroCurve = new ARM_ZeroLInterpol(psMatu,(ARM_Date) asOfDate,
									             mktData, MMVsFut, SwapVsFut, 
									             Raw, interp, 0, aCcy);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}

			else
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (mktData)
		   delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_CreateZCSwapFutIntSmooth (double Date,
										VECTOR<CCString>& matu,
										VECTOR<double>& rate,
										long MMVsFut,
										long SwapVsFut,
										long Raw,
										long interp,
										const CCString& Ccy,
										double lambda,
										long prec,
										ARM_result& result,
										long objId)
{
	long curveId;

    ARM_ZeroInterpolation* prevZcLin = NULL;
    ARM_ZeroInterpolation* createdZeroCurve = NULL;

    ARM_Vector* mktData=NULL;
 
	int real_size = matu.size ();
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char asOfDate[11];
	ARM_Currency* aCcy = NULL;
	char* tmp = NULL;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(Date,asOfDate);

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "XX");

		for( int j = 0; j < real_size; j++)
		{
			sprintf(psMatu[j], "%s", (const char*)matu[j]);
		}
 
		tmp = Ccy.GetStr();

		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			if (tmp)
				delete tmp;
		}

		tmp = NULL;

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZeroCurve = new ARM_ZeroInterpolation(psMatu, (ARM_Date) asOfDate, 
									mktData, MMVsFut, SwapVsFut, 
									Raw, interp, aCcy, lambda, prec);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroInterpolation *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_INTERPOLATION) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}

			else
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (mktData)
		   delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_ZCINTSMOOTH (long inCvId,
						   double lambda,
						   long prec,
						   ARM_result& result,
						   long objId)
{
	long curveId;

	ARM_ZeroInterpolation* zc = NULL;
	ARM_ZeroInterpolation* createdZeroCurve = NULL;
	ARM_ZeroLInterpol* inCurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		inCurve = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(inCvId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inCurve, ARM_ZERO_LIN_INTERPOL) == 0)
		{
			result.setMsg ("ARM_ERR: In Curve is not of a good type");
			return ARM_KO;
		}

		createdZeroCurve = new ARM_ZeroInterpolation(inCurve, lambda, prec);

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_ZeroInterpolation *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_INTERPOLATION) == 1)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}

			else
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

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
}




long ARMLOCAL_zcswapcubdiff (double date,
							 VECTOR<CCString>& matu,
							 VECTOR<double>& rate,
							 long mmVsFut,
							 long swapVsFut,
							 long raw,
							 long interp,
							 const CCString& ccy,
							 ARM_result& result,
							 long objId)
{
	long curveId;

	if(matu.size () != rate.size ())
	{
		result.setMsg ("ARM_ERR: maturities and rates must have same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	long nblines = matu.size ();

	ARM_ZeroCubDiff* createdZeroCurve = NULL;
	ARM_ZeroCubDiff* prevZcCubDiff = NULL;
	ARM_Currency* aCcy = NULL;
    ARM_Vector* mktData = NULL;

	char* tmp = NULL;
	char myCurveDate[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(date,myCurveDate);

		tmp = (char*) ccy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			free(tmp);
		}
		tmp = NULL;

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "X");

		for(int j = 0; j < nblines; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZeroCurve = new ARM_ZeroCubDiff((ARM_Date) myCurveDate, psMatu,
												mktData, mmVsFut, swapVsFut,
												raw, interp, 0, aCcy);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcCubDiff = (ARM_ZeroCubDiff *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcCubDiff, ARM_ZERO_CUBDIFF) == 1)
			{
				if (prevZcCubDiff)
				{
					delete prevZcCubDiff;
					prevZcCubDiff = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmp)
			delete tmp;
		tmp = NULL;

		if (createdZeroCurve)
			delete createdZeroCurve;
		createdZeroCurve = NULL;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
    }
}




long ARMLOCAL_ZCINTSMOOTH (const VECTOR<double>& matu,
						   const VECTOR<double>& rate,
						   double aDate,
						   long meth,
						   double lambda,
						   long prec,
						   ARM_result& result,
						   long objId)
{
	long zcId;

	/*--- parameters checking ---*/
	if(matu.size () != rate.size ())
	{
		result.setMsg ("ARM_ERR: maturities and rates must have same size");
		return ARM_KO;
	}

	int real_size = matu.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_ZeroInterpolation* zc = NULL;
	ARM_ZeroInterpolation* newZc = NULL;

	ARM_Vector* vMatu = NULL;
	ARM_Vector* vRates = NULL;

	double * pdRate = NULL;
	double * pdMatu = NULL;
	char myCurveDate[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(aDate,myCurveDate);

		pdRate = new double[real_size];
		pdMatu = new double[real_size];
		for(int i = 0; i < real_size; i++)
		{
			pdRate[i] = rate [i];
			pdMatu[i] = matu [i];
		}

		vRates = new ARM_Vector(real_size,pdRate);
		vMatu = new ARM_Vector(real_size,pdMatu);

		newZc = new ARM_ZeroInterpolation((ARM_Date) myCurveDate, vMatu, 
													vRates, meth, lambda, prec);

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;

		if (vRates)
			delete [] vRates;
		vRates = NULL;

		if (vMatu)
			delete [] vMatu;
		vMatu = NULL;

		if (newZc == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			zcId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc);

			if (zcId == RET_KO)
			{
				if (newZc)
					delete newZc;
				newZc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(zcId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_ZeroInterpolation *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_INTERPOLATION) == 1)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newZc)
			delete newZc;
		newZc = NULL;

		if (pdRate)
			delete [] pdRate;
		pdRate = NULL;

		if (pdMatu)
			delete [] pdMatu;
		pdMatu = NULL;

		if (vRates)
			delete [] vRates;
		vRates = NULL;

		if (vMatu)
			delete [] vMatu;
		vMatu = NULL;

		ARM_RESULT();
    }
}



ARM_ZeroLInterpol* GetZCFromSummitNoETK(const CCString& index,
							           const CCString& currency,
							           const CCString& cvName,
							           ARM_Date asof,
									   long interpId,
                                       ARM_result& result)
{
    ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	ARM_Currency* aCcy = NULL;

	FILE *Fp = NULL;

	ARM_Date myDateFromFile=asof;
	ARM_Date dateToday;

	if (myDateFromFile > dateToday)
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName, FileName2;

	double *pdRate = NULL;
	double *pdMatu = NULL;
	
	CCString msg (" ");

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + ZC_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + ZC_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);


	FileName = myRepertory + (CCString)"ZC_" + currency + "_" + index + "_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"ZC_" + currency + "_" + index + "_" + cvName + ".";

	char sEch[50];
	char buffer[50];

	if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
	{
		_ltoa(myDateFromFile.GetYear(),buffer,10);
		FileName += (CCString) (buffer);
		FileName2 += (CCString) (buffer);

		_ltoa(myDateFromFile.GetMonth(),buffer,10);
		if (myDateFromFile.GetMonth() < 10)
		{
			FileName += "0";
			FileName2 += "0";
		}
		FileName += (CCString) buffer;
		FileName2 += (CCString) buffer;

		_ltoa(myDateFromFile.GetDay(),buffer,10);
		if (myDateFromFile.GetDay() < 10)
		{
			FileName += "0";
			FileName2 += "0";
		}
		FileName += (CCString) buffer;
		FileName2 += (CCString) buffer;
	}
	else
	{
		_ltoa(myDateFromFile.GetDay(),buffer,10);
		if (myDateFromFile.GetDay() < 10)
		{
			FileName += "0";
			FileName2 += "0";
		}
		FileName += (CCString) buffer;
		FileName2 += (CCString) buffer;

		_ltoa(myDateFromFile.GetMonth(),buffer,10);
		if (myDateFromFile.GetMonth() < 10)
		{
			FileName += "0";
			FileName2 += "0";
		}
		FileName += (CCString) buffer;
		FileName2 += (CCString) buffer;

		_ltoa(myDateFromFile.GetYear(),buffer,10);
		FileName += (CCString) (buffer + 2);
		FileName2 += (CCString) (buffer + 2);
	}

	FileName += (CCString) ".000";
	FileName2 += (CCString) ".000";

	double val;
	int rc = 0;
	int j = 0;

	// latest curves
	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		// historical curves
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg( CCString("ARM_ERR: Could not open the following files:\n" )
						+ FileName + "\n" + FileName2 + "\nCheck parameters ..." );
			return NULL;
		}
	}

	char* tmp = (char*) currency;

	pdRate = new double[ARM_NB_TERMS];
	pdMatu = new double[ARM_NB_TERMS];

	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",buffer);
		rc = fscanf(Fp, "%s",sEch);
		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			ARM_Date tmpDate(sEch,"MM/DD/YYYY");
			pdMatu[j] = ((double)(tmpDate-asof))/365.;
		}
		else
		{
			ARM_Date tmpDate(sEch,"DD/MM/YYYY");
			pdMatu[j] = ((double)(tmpDate-asof))/365.;
		}
		rc = fscanf(Fp, "%lf",&val);
		pdRate[j] = val;

		j++;
	}

	fclose(Fp);

	if (tmp)
	{
		aCcy = new ARM_Currency(tmp);
		delete tmp;
	}
	tmp = NULL;

	mktData = new ARM_Vector(j-1,pdRate);
	mktMatu = new ARM_Vector(j-1,pdMatu);

	ARM_ZeroLInterpol* createdZcLin = new ARM_ZeroLInterpol(asof, mktMatu,
						mktData, 0, 0, interpId);

	createdZcLin->SetCurrencyUnit(aCcy);
	
	if (pdRate)
		delete [] pdRate;
	pdRate = NULL;

	if (mktData)
		delete mktData;
	mktData = NULL;

	if (pdMatu)
		delete [] pdMatu;
	pdMatu = NULL;

	if (mktMatu)
		delete mktMatu;
	mktMatu = NULL;

	if (aCcy)
		delete aCcy;
	aCcy = NULL;

    return createdZcLin;
}


ARM_ZeroLInterpol* obj_getZcFromSummit (const CCString& index,
										const CCString& currency,
										const CCString& cvName,
										double aSdate,
										long interpId,
										ARM_result& result)
{
	ARM_ZeroLInterpol* zc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(aSdate,sDate);

	if (GetDataRetrieverVersion () >= ETKRETRIEVER)
	{
		CCString xmlResponse;

		xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,currency,cvName,(ARM_Date)sDate);

		int interp = ARMLOCAL_ParseXMLForInterpolator(xmlResponse);

		ARM_Date myDate(sDate);

		xmlResponse = etoolkit_getXMLZCFromSummit(index,currency,cvName,myDate);

		zc = ARMLOCAL_ParseXMLForZC(xmlResponse, myDate, (const char*) currency, interp);
	}
	else
	{
        ARM_Date asOfDate(sDate);
	    zc = GetZCFromSummitNoETK(index,
								  currency,
								  cvName,
								  asOfDate,
								  interpId,
								  result);

		if ( (!zc) && (GetFallBackDataRetrieverVersion () >= ETKRETRIEVER) )
		{			
			CCString xmlResponse;

			xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,currency,cvName,(ARM_Date)sDate);

			int interp = ARMLOCAL_ParseXMLForInterpolator(xmlResponse);

			xmlResponse = etoolkit_getXMLZCFromSummit(index,currency,cvName,asOfDate);

			zc = ARMLOCAL_ParseXMLForZC(xmlResponse, asOfDate,(const char*) currency, interp);
		}
	}

	return zc;
}

//	--------------------------------------------------------------------------------------------
ARM_ZeroLInterpol* 
obj_getZcFromCalypso (const ARM_Date& AsOf,
					const std::string & index,
					const std::string & ccy,
					const std::string & term,
					const std::string & PricingEnv,
					const std::string & ForceCurveName,
					long interpId,
					const std::string & xmlFileName)
{
	ARM_ZeroLInterpol* zc = NULL;
	std::string xmlContent ;
	ARM_CalypsoToolkit::GetCurveZero(index,ccy,term,ForceCurveName,PricingEnv,AsOf,xmlFileName,xmlContent); 
	return ARMLOCAL_ParseXMLForCalypsoZC(xmlContent,AsOf,ccy,interpId);	
}



long ARMLOCAL_GetInitialCurveFromSummit (const CCString& index,
										 const CCString& currency,
										 const CCString& cvName,
										 double aSdate,
										 long adjOrNotId,
										 VECTOR<CCString>* matu,
										 VECTOR<double>* yield,
										 ARM_result& result)
{
	CCString msg(" ");

	char* sDate = new char[11];

	Local_XLDATE2ARMDATE(aSdate,sDate);

	ARM_Date myDate(sDate);

	if (sDate)
		delete [] sDate;
	sDate = NULL;

	if (GetDataRetrieverVersion() >= ETKRETRIEVER)
	{
		CCString xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,currency,cvName,myDate);

		long retCode;

		try
		{
			retCode = ARMLOCAL_ParseXMLForMY(xmlResponse,adjOrNotId,matu,yield);
			result.setLong(GetDataRetrieverVersion());
		}
	
        catch(Exception& x)
		{
			x.DebugPrint();

			ARM_RESULT();
		}

		return retCode;
	}
	else
	{
		FILE *Fp = NULL;

		CCString myRepertory ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION);
		CCString myRepertory2 ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

		CCString FileName = myRepertory + (CCString)"YLD_" + currency + "_" + index + "_" + cvName + ".";
		CCString FileName2 = myRepertory2 + (CCString)"YLD_" + currency + "_" + index + "_" + cvName + ".";

		ARM_Date dateToday;

	//	dateToday.SysToday();// donne la date systeme du jour

		char buffer[50];

		int trouve = 0;
		int compteur = 0;

		if (myDate > dateToday)
		{
			result.setMsg ("ARM_ERR: Invalid Date");
			return ARM_KO;
		}

		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			_ltoa(myDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer);
			FileName2 += (CCString) (buffer);

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;
		}
		else
		{
			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer + 2);
			FileName2 += (CCString) (buffer + 2);
		}

		FileName += (CCString) ".000";
		FileName2 += (CCString) ".000";

		if ((Fp = fopen(FileName,"r")) == NULL)
		{
			if ((Fp = fopen(FileName2,"r")) == NULL)
			{
				if (GetFallBackDataRetrieverVersion() != 0)
				{
					CCString msg (" ");

					CCString xmlResponse;
					CCString msgList;

					char sDate[11];
					Local_XLDATE2ARMDATE(aSdate,sDate);
					ARM_Date myDate(sDate);

					sprintf(sDate, "%04d%02d%02d", myDate.GetYear(), myDate.GetMonth(), myDate.GetDay());

					CCString myMarket("s_market:MktYieldRead");
					CCString myRequest;
					myRequest = (CCString)"<Request><Ccy>" + currency + (CCString)"</Ccy><Index>" + index + (CCString)"</Index><AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate><CurveId>" + cvName + (CCString)"</CurveId></Request>";

					if (etoolkit_execute(myMarket,myRequest,xmlResponse,msgList)==ARM_KO)
						return ARM_KO;

					const int buffSize = 1024;

					long retCode;

					try
					{
						retCode = ARMLOCAL_ParseXMLForMY(xmlResponse,adjOrNotId,matu,yield);
						result.setLong(GetFallBackDataRetrieverVersion());
					}
				
					catch(Exception& x)
					{
						x.DebugPrint();

						ARM_RESULT();
					}

					return retCode;
				}
				else
				{
					result.setMsg( CCString("ARM_ERR: Could not open the following files:\n" )
					+ FileName + "\n" + FileName2 + "\nCheck parameters ..." );
					return ARM_KO;
				}
			}
			else
			{
				fclose(Fp);

				result.setString(FileName2);
				result.setLong(GetDataRetrieverVersion());
				return ARM_OK;
			}
		}
		else
		{
			fclose(Fp);

			result.setString(FileName);
			result.setLong(GetDataRetrieverVersion());
			return ARM_OK;
		}
	}
}

long ARMLOCAL_GetZCFromSummit (const CCString& index,
							   const CCString& currency,
							   const CCString& cvName,
							   double aSdate,
							   long interpId,
							   ARM_result& result,
							   long objId)
{
	long curveId;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	ARM_ZeroLInterpol* prevZcLin = NULL;

	CCString msg (" ");
	try
	{
		createdZcLin = obj_getZcFromSummit(index, currency, cvName, aSdate, interpId, result);
		if ( createdZcLin == NULL )
		{
		   result.setMsg("Object is Null");
			
           return(ARM_KO);
		}

        // Update Mkt data characteristics
		string	vType("ZERO");
		string	vIndex((const char*) index);
		string	vCurrency((const char*) currency);
		string	vCrvId((const char*) cvName);

        if( index == "BS" )
			vType = "BASIS USD";

        createdZcLin->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);
	}

    catch (Exception& x)
	{
		x.DebugPrint();

		if (createdZcLin)
			delete createdZcLin;
		createdZcLin = NULL;

		ARM_RESULT();
	}

	if(objId == -1)
	{
		CREATE_GLOBAL_OBJECT();

		curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin);

		if (curveId == RET_KO)
		{
			if (createdZcLin)
				delete createdZcLin;
			createdZcLin = NULL;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return ARM_KO;
		}

		result.setLong(curveId);

		return ARM_OK;
	}
	else
	{
		prevZcLin = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_LIN_INTERPOL) == 1)
		{
			if (prevZcLin)
			{
				delete prevZcLin;
				prevZcLin = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin, objId);

			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}
	}

}


//	--------------------------------------------------------------------------
long
ARMLOCAL_GetZCFromCalypso (const ARM_Date& AsOf,
							const std::string& index,
							const std::string& ccy,
							const std::string& term,
							const std::string& PricingEnv,
							const std::string& ForceCurveName,
							long interpId,
							const std::string& xmlFileName,
							long objId)
{
	// long curveId;
	// ARM_ZeroLInterpol* createdZcLin = NULL;
	// ARM_ZeroLInterpol* prevZcLin = NULL;



	//	creates the instance & acquires ownership	
	std::auto_ptr<ARM_ZeroLInterpol> createdZcLin (  obj_getZcFromCalypso(AsOf,
													index,
													ccy,
													term,
													PricingEnv,
													ForceCurveName,
													interpId,
													xmlFileName) ) ;
	
	// creation might return ull pointer
	if (!createdZcLin.get())
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_GetZCFromCalypso: Can't create the curve"); 

	// Update Mkt data characteristics : really located here ??? 
	string	vType("ZERO");
	if( index == "BS" ) vType = "BASIS USD";
	createdZcLin->SetMktExternalCharacteristics(unconst(index), unconst(ccy), unconst(PricingEnv), unconst(vType));

	//	insertion in the cache
	if (objId==-1) 
	{
		CREATE_GLOBAL_OBJECT();
		long curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdZcLin.get());
		if (curveId==RET_KO) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_GetZCFromCalypso:: Can't insert new object"); 
		createdZcLin.release(); 
		return curveId; 
	}
	
	ARM_ZeroLInterpol * prevZcLin = dynamic_cast<ARM_ZeroLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId) );
	if (prevZcLin) delete prevZcLin ; 

	LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdZcLin.get(), objId);
	createdZcLin.release(); 
	return objId; 
}

long ARMLOCAL_CreateZCFromSummit (const CCString& index,
								  const CCString& currency,
								  const CCString& cvName,
								  double aSdate,
								  long adjOrNotId,
								  const CCString& raw,
								  long swapFrqId,
								  ARM_result& result,
								  long objId)
{
/*	if (GetETKVersion () == 0)
	{
		result.setMsg ("ARM_ERR: This function is not implemented without ETK");
		return ARM_KO;
	}
*/
	long curveId;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	ARM_ZeroLInterpol* prevZcLin = NULL;

	CCString msg (" ");
	
	CCString xmlResponse;

	char sDate[11];
	Local_XLDATE2ARMDATE(aSdate,sDate);
	ARM_Date myDate(sDate);

	long retCode;

	try
	{
		createdZcLin = ARMLOCAL_CreateZC(index, currency, cvName, 
                                         myDate, raw, adjOrNotId, swapFrqId);

		if ( createdZcLin == NULL )
		{
		   result.setMsg("ARM_ERR: Curve is null");
			
           return(ARM_KO);
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin);

			if (curveId == RET_KO)
			{
				if (createdZcLin)
					delete createdZcLin;
				createdZcLin = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
    catch (Exception& x)
	{
		if (createdZcLin)
			delete createdZcLin;
		createdZcLin = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}

	return retCode;
}


long ARMLOCAL_CreateZCSpreadedFromSummit(const CCString& index,
										 const CCString& currency,
										 const CCString& cvName,
										 double aSdate,
										 long adjOrNotId,
										 const CCString& raw,
										 long swapFreq,
										 long mmFreq,
										 long interpId,
										 long zcInitId,
										 ARM_result& result,
										 long objId)
{
	long curveId;
	ARM_ZeroLInterpol*	zc = NULL;
	ARM_BasisCurve*	newZc = NULL;
	ARM_ZeroCurve* ZCInit = NULL;

	if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return(ARM_KO);
	}

	char myCurveDate[20];

	CCString msg("");

	long retCode;

    try
    {
		Local_XLDATE2ARMDATE(aSdate, myCurveDate);

		if (zcInitId != -1)
		{
			ZCInit = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcInitId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZCInit, ARM_ZERO_CURVE) == 0)
			{
			   result.setMsg("ARM_ERR: ZCInit is not of a good type");
				
			   return(ARM_KO);
			}
		}
		if( currency == "USD" )
		{
			result.setMsg ("ARM_ERR: Can not create a ZcSpreaded curve for USD");
			return ARM_KO;
		}

		newZc = ARMLOCAL_CreateZCSpreaded(index, currency, cvName, myCurveDate, adjOrNotId, raw, swapFreq, mmFreq, interpId, ZCInit);

		if (newZc == NULL)
		{
			result.setMsg ("ARM_ERR: ZcSpread is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
				
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc);

			if (curveId == RET_KO)
			{
				delete newZc;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = dynamic_cast<ARM_ZeroLInterpol*> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_BASIS_CURVE) == 1)
			{
				delete zc;
				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newZc, objId);
				return ARM_OK;
			}
			else
			{
				delete newZc;
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
		delete newZc;
		ARM_RESULT();
    }

	return retCode;
}

/**********/
void ARMLOCAL_GetInitialCurveFromCalypso(const ARM_Date& AsOf,
										 const std::string & pricingEnv,
										 const std::string & index,
										 const std::string & currency,
										 const std::string & term,
										 const std::string & forceCurveName,
										 const std::string & xmlFileName,
										 bool doAdjust ,
										 std::vector<std::string>& matus,
										 std::vector<double>&yields)
{
	std::string xml ;
	ARM_CalypsoToolkit::GetCurveZero(index,currency,term,forceCurveName,pricingEnv,AsOf,xmlFileName,xml); 
	ARMLOCAL_ParseXMLForMYCalypso(xml,doAdjust,matus,yields); 
/**
	if (GetDataRetrieverVersion() >= ETKRETRIEVER)
	{
		CCString msg (" ");

		CCString xmlResponse;
		CCString msgList;

		char sDate[11];
		Local_XLDATE2ARMDATE(aSdate,sDate);
		ARM_Date myDate(sDate);

		sprintf(sDate, "%04d%02d%02d", myDate.GetYear(), myDate.GetMonth(), myDate.GetDay());

		CCString myMarket("s_market:MktYieldRead");
		CCString myRequest;
		myRequest = (CCString)"<Request><Ccy>" + currency + (CCString)"</Ccy><Index>" + index + (CCString)"</Index><AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate><CurveId>" + cvName + (CCString)"</CurveId></Request>";

		if (etoolkit_execute(myMarket,myRequest,xmlResponse,msgList)==ARM_KO)
			return ARM_KO;

		const int buffSize = 1024;

		long retCode;

		try
		{
			retCode = ARMLOCAL_ParseXMLForMY(xmlResponse,adjOrNotId,matu,yield);
			result.setLong(GetDataRetrieverVersion());
		}
	
        catch(Exception& x)
		{
			x.DebugPrint();

			ARM_RESULT();
		}

		return retCode;
	}
	else
	{
		FILE *Fp = NULL;

		CCString myRepertory ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION);
		CCString myRepertory2 ((CCString)(armlocal_init->data_folder.c_str()) + YLD_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

		CCString FileName = myRepertory + (CCString)"YLD_" + currency + "_" + index + "_" + cvName + ".";
		CCString FileName2 = myRepertory2 + (CCString)"YLD_" + currency + "_" + index + "_" + cvName + ".";
		char* sDate = new char[11];

		Local_XLDATE2ARMDATE(aSdate,sDate);

		ARM_Date myDate(sDate);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		ARM_Date dateToday;

	//	dateToday.SysToday();// donne la date systeme du jour

		char buffer[50];

		int trouve = 0;
		int compteur = 0;

		if (myDate > dateToday)
		{
			result.setMsg ("ARM_ERR: Invalid Date");
			return ARM_KO;
		}

		if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		{
			_ltoa(myDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer);
			FileName2 += (CCString) (buffer);

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;
		}
		else
		{
			_ltoa(myDate.GetDay(),buffer,10);
			if (myDate.GetDay() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetMonth(),buffer,10);
			if (myDate.GetMonth() < 10)
			{
				FileName += "0";
				FileName2 += "0";
			}
			FileName += (CCString) buffer;
			FileName2 += (CCString) buffer;

			_ltoa(myDate.GetYear(),buffer,10);
			FileName += (CCString) (buffer + 2);
			FileName2 += (CCString) (buffer + 2);
		}

		FileName += (CCString) ".000";
		FileName2 += (CCString) ".000";

		if ((Fp = fopen(FileName,"r")) == NULL)
		{
			if ((Fp = fopen(FileName2,"r")) == NULL)
			{
				if (GetFallBackDataRetrieverVersion() != 0)
				{
					CCString msg (" ");

					CCString xmlResponse;
					CCString msgList;

					char sDate[11];
					Local_XLDATE2ARMDATE(aSdate,sDate);
					ARM_Date myDate(sDate);

					sprintf(sDate, "%04d%02d%02d", myDate.GetYear(), myDate.GetMonth(), myDate.GetDay());

					CCString myMarket("s_market:MktYieldRead");
					CCString myRequest;
					myRequest = (CCString)"<Request><Ccy>" + currency + (CCString)"</Ccy><Index>" + index + (CCString)"</Index><AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate><CurveId>" + cvName + (CCString)"</CurveId></Request>";

					if (etoolkit_execute(myMarket,myRequest,xmlResponse,msgList)==ARM_KO)
						return ARM_KO;

					const int buffSize = 1024;

					long retCode;

					try
					{
						retCode = ARMLOCAL_ParseXMLForMY(xmlResponse,adjOrNotId,matu,yield);
						result.setLong(GetFallBackDataRetrieverVersion());
					}
				
					catch(Exception& x)
					{
						x.DebugPrint();

						ARM_RESULT();
					}

					return retCode;
				}
				else
				{
					result.setMsg( CCString("ARM_ERR: Could not open the following files:\n" )
					+ FileName + "\n" + FileName2 + "\nCheck parameters ..." );
					return ARM_KO;
				}
			}
			else
			{
				fclose(Fp);

				result.setString(FileName2);
				result.setLong(GetDataRetrieverVersion());
				return ARM_OK;
			}
		}
		else
		{
			fclose(Fp);

			result.setString(FileName);
			result.setLong(GetDataRetrieverVersion());
			return ARM_OK;
		}
	} 
	**/ 
}




long ARMLOCAL_zcvsk (const VECTOR<double>& param,
					 double date,
					 ARM_result& result,
					 long objId)
{
	long curveId=0;

    ARM_ZeroVasicek* zc = NULL;
    ARM_ZeroVasicek* newzc = NULL;
    ARM_ZeroVasicek* createdZeroCurveVasicek = NULL;

	int param_size = param.size ();

	double* pparam = new double[param_size];

	char* sDate = new char[11];
	CCString msg (" ");

	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	if (sDate)
		delete [] sDate;
	sDate = NULL;
	

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	try
	{
		for(int i = 0; i < param_size; i++)
		{
			pparam[i] = param[i];
		}

		newzc = new ARM_ZeroVasicek(myDate, pparam);

		if (pparam)
			delete [] pparam;
		pparam = NULL;

		if (objId != -1)
		{
			zc = (ARM_ZeroVasicek *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_VASICEK) == 0)
			{
				result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
				return ARM_KO;
			}

			(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc, objId);
			
			result.setLong(objId);

			if (zc)
				delete zc;

			return ARM_OK;

		}
		else
		{
			curveId=LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc);

			result.setLong(curveId);
		}

		return ARM_OK;
	}

    catch(Exception& x)

    {
		x.DebugPrint();
		
		if (pparam)
			delete [] pparam;
		pparam = NULL;

		if (newzc)
			delete newzc;
		newzc = NULL;

		ARM_RESULT();
    }

}



long ARMLOCAL_zcsplicub (VECTOR<double>& matu,
						 VECTOR<double>& rate,
						 long meth,
						 double date,
						 long lastBucket,
						 ARM_result& result,
						 long objId)
{
	long curveId=0;

    ARM_ZeroSpliCub* zc = NULL;
    ARM_ZeroSpliCub* newzc = NULL;

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	int size = matu.size ();

	ARM_Vector* vMatu = new ARM_Vector(size);
	ARM_Vector* vRate = new ARM_Vector(size);

	char* sDate = new char[11];
	CCString msg (" ");

	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	if (sDate)
		delete [] sDate;
	sDate = NULL;
	

	try
	{
		for(int i = 0; i < size; i++)
		{
			vMatu->Elt(i) = matu[i];
			vRate->Elt(i) = rate[i];
		}

		newzc = new ARM_ZeroSpliCub(myDate, vMatu, vRate, meth, lastBucket);

		if (lastBucket == 0)
		{
			if (vMatu)
				delete vMatu;
			vMatu = NULL;

			if (vRate)
				delete vRate;
			vRate = NULL;
		}

		if (objId != -1)
		{
			zc = (ARM_ZeroSpliCub *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_SPLICUB) == 0)
			{
				result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
				return ARM_KO;
			}

			(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc, objId);
			
			result.setLong(objId);

			if (zc)
				delete zc;
			zc = NULL;

			return ARM_OK;
		}
		else
		{
			curveId=LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc);

			result.setLong(curveId);
		}

		return ARM_OK;
	}

    catch(Exception& x)

    {
		x.DebugPrint();

/*		if (vMatu)
			delete vMatu;
		vMatu = NULL;

		if (vRate)
			delete vRate;
		vRate = NULL;
*/
		if (newzc)
			delete newzc;
		newzc = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_zcspl (const VECTOR<double>& param,
					 double date,
					 ARM_result& result,
					 long objId)
{
	long curveId=0;

    ARM_ZeroSplines* zc = NULL;
    ARM_ZeroSplines* newzc = NULL;

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	int size = param.size ();

	double* dParam = new double[size];

	char* sDate = new char[11];
	CCString msg (" ");

	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	if (sDate)
		delete [] sDate;
	sDate = NULL;
	

	try
	{
		for(int i = 0; i < size; i++)
			dParam[i] = param[i];

		newzc = new ARM_ZeroSplines(myDate,size,dParam);

		if (objId != -1)
		{
			zc = (ARM_ZeroSplines *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_ZERO_SPLINES) == 0)
			{
				result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
				return ARM_KO;
			}

			(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc, objId);
			
			result.setLong(objId);

			if (zc)
				delete zc;
			zc = NULL;

			return ARM_OK;
		}
		else
		{
			curveId=LOCAL_PERSISTENT_OBJECTS->SetPersistent(newzc);

			result.setLong(curveId);
		}

		return ARM_OK;
	}

    catch(Exception& x)

    {
		x.DebugPrint();

		if (newzc)
			delete newzc;
		newzc = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_CreateTOYNYZCSwapInt (double Date,
									VECTOR<CCString>& matu,
									VECTOR<double>& rate,
									long MMVsFut,
									long SwapVsFut,
									long Raw,
									long interp,
									const CCString& Ccy,
									long frq,
									ARM_result& result,
									long objId)
{
	long curveId;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	ARM_ZeroLInterpol* prevZcLin = NULL;
	ARM_Currency* aCcy = NULL;
    ARM_Vector* mktData = NULL;
	int i=0;
	int real_size = matu.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* tmp = NULL;
	char myCurveDate[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(Date,myCurveDate);

		tmp = (char*) Ccy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			free(tmp);
		}
		tmp = NULL;

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "X");

		for(int j = 0; j < real_size; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZcLin = new ARM_ZeroLInterpol((ARM_Date)myCurveDate, psMatu,
							mktData, MMVsFut, SwapVsFut, 
							Raw, interp, 0, frq, aCcy);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZcLin == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin);

			if (curveId == RET_KO)
			{
				if (createdZcLin)
					delete createdZcLin;
				createdZcLin = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcLin = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcLin, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (prevZcLin)
				{
					delete prevZcLin;
					prevZcLin = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmp)
			delete tmp;
		tmp = NULL;

		if (createdZcLin)
			delete createdZcLin;
		createdZcLin = NULL;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_BumpCurve (long cvId,
						 VECTOR<CCString>& matu,
						 VECTOR<double>& epsilon,
						 long meth,
						 ARM_result& result,
						 long objId)
{
	long curveId;

	ARM_Object* zc = NULL;
	ARM_InfCurv* infZc = NULL;
	ARM_ZeroCurve* createdZeroCurve = NULL;
	ARM_InfCurv* createdInfZeroCurve = NULL;
	ARM_ZeroCurve* inCurve = NULL;
	ARM_InfCurv* inInfCurve = NULL;
	ARM_Vector* vEpsilon = NULL;
	ARM_Object* newCurv = NULL;

	int i=0;
	int real_size = matu.size ();
	int value_size = epsilon.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	vector<string> vMatu;
	vector<double> vInfEpsilon;

	CCString msg (" ");

	try
	{
		inCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(cvId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inCurve, ARM_ZERO_CURVE) == 0)
		{
			inInfCurve = (ARM_InfCurv *) LOCAL_PERSISTENT_OBJECTS->GetObject(cvId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inCurve, ARM_INFCURV) == 0)
			{
				result.setMsg ("ARM_ERR: In Curve is not a Zc Curve");
				return ARM_KO;
			}
		}

		ARM_CRV_TERMS psMatu;

		if ( (real_size == 0) && (value_size == 1) )
		{
			if (inInfCurve)
			{
				for(int j = 0; j < inInfCurve->GetMktTerms().size(); j++)
				{
					vMatu.push_back(inInfCurve->GetMktTerms()[j]);
					vInfEpsilon.push_back(epsilon[0]);
				}
			}
			else
			{
				for(int j = 0; j < inCurve->GetMktData()->itsMktValue->GetSize(); j++)
				{
					sprintf(psMatu[j], "%s", (const char*)(inCurve->GetMktData()->itsMktTerms[j]));
				}

				vEpsilon = new ARM_Vector(inCurve->GetMktData()->itsMktValue->GetSize(),epsilon[0]);
			}
		}
		else
		{
			for (int j = 0; j < ARM_NB_TERMS; j++)
				sprintf(psMatu[j], "%s", "X");

			for(int j = 0; j < real_size; j++)
			{
				if (inInfCurve)
				{
					/// convert dates to Julian type
					string maturity( "yYMmWwDd" );
					string s;
					char buffer[20];

					for( i = 0; i < matu.size(); i++)
					{
						s = matu[i];
						if( s.find_first_of( maturity ) == string::npos )
						{
							double julianDate = atof( s.c_str() );
							sprintf( buffer, "%f", XLDateToJulian( julianDate ) );
							s = buffer;
						}
						vMatu.push_back(s);
					}
					vInfEpsilon.push_back(epsilon[j]);
				}
				else
				{
					sprintf(psMatu[j], (const char *) matu[j]);
				}
			}

			vEpsilon = CreateARMVectorFromVECTOR(epsilon);
		}

		// BUMP d'une courbe Inflation
		if (inInfCurve != NULL)
		{
			createdInfZeroCurve = inInfCurve->GenerateShiftCurve(vMatu,vInfEpsilon);
			newCurv = createdInfZeroCurve;
		}
		else
		{
			if (meth == 0)
				createdZeroCurve = inCurve->GenerateShiftCurve(psMatu,vEpsilon);
			else
				createdZeroCurve = inCurve->GenerateShiftCurveFwd(psMatu,vEpsilon);

			newCurv = createdZeroCurve;

			if (vEpsilon)
				delete vEpsilon;
			vEpsilon = NULL;
		}

		if (newCurv == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurv);

			if (curveId == RET_KO)
			{
				if (newCurv)
					delete newCurv;
				newCurv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 1)
				|| (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(zc, ARM_INFCURV) == 1)
				)
			{
				if (zc)
				{
					delete zc;
					zc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurv, objId);

				return ARM_OK;
			}

			else
			{
				if (newCurv)
					delete newCurv;
				newCurv = NULL;

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
}


long ARMLOCAL_BumpSpreadedCurve (long cvId,
										VECTOR<CCString>& matu,
										VECTOR<double>& epsilon,
										CCString& curveToBump,
										long meth,
										ARM_result& result,
										long objId )
{
	long curveId;

	ARM_Object* zc = NULL;
	ARM_BasisCurve* createdZeroCurve = NULL;
	ARM_BasisCurve* inCurve = NULL;
	ARM_Vector* vEpsilon = NULL;
	ARM_Object* newCurv = NULL;

	int i=0;
	int real_size = matu.size ();
	int value_size = epsilon.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	vector<string> vMatu;

	CCString msg (" ");

	try
	{
		inCurve = (ARM_BasisCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(cvId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inCurve, ARM_BASIS_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In Curve is not an ARM_BASIS_CURVE");
			return ARM_KO;
		}

		ARM_CRV_TERMS psMatu;

		if ( (real_size == 0) && (value_size == 1) )
		{
			for(int j = 0; j < inCurve->GetMktData()->itsMktValue->GetSize(); j++)
			{
				sprintf(psMatu[j], "%s", (const char*)(inCurve->GetMktData()->itsMktTerms[j]));
			}

			vEpsilon = new ARM_Vector(inCurve->GetMktData()->itsMktValue->GetSize(), epsilon[0]);
		}
		else
		{
			for (int j = 0; j < ARM_NB_TERMS; j++)
				sprintf(psMatu[j], "%s", "X");

			for(int j = 0; j < real_size; j++)
			{
				sprintf(psMatu[j], (const char *) matu[j]);
			}

			vEpsilon = CreateARMVectorFromVECTOR(epsilon);
		}

		if (meth == 0)
			createdZeroCurve = inCurve->GenerateShiftCurve(psMatu,vEpsilon,curveToBump);
		else
			createdZeroCurve = inCurve->GenerateShiftCurveFwd(psMatu,vEpsilon,curveToBump);

		newCurv = createdZeroCurve;

		delete vEpsilon;


		if (newCurv == NULL)
		{
			result.setMsg ("ARM_ERR: createdZeroCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurv);

			if (curveId == RET_KO)
			{
				delete newCurv;

				result.setMsg ("ARM_ERR: Pb with inserting object");				

				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			zc = (ARM_Object *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 1)
			{
				delete zc;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurv, objId);

				return ARM_OK;
			}

			else
			{
				delete newCurv;

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
}


long ARMLOCAL_CreateZCCashInt (double date,
							   VECTOR<CCString>& matu,
							   VECTOR<double>& rate,
							   VECTOR<long>& bondsId,
							   VECTOR<double>& yields,
							   long MMVsFut,
							   const CCString& Ccy,
							   ARM_result& result,
							   long objId)
{
	long curveId;

	ARM_ZeroLInterpol* createdZeroCurve = NULL;
	ARM_ZeroLInterpol* oldZeroCurve = NULL;

	ARM_Vector* mktData = NULL;
	ARM_Container* Bonds = NULL;
	ARM_Vector* Yields = NULL;
	ARM_Bond* Bond = NULL;

	ARM_Currency* aCcy = NULL;

	int rateSize = rate.size ();

	int szBd = bondsId.size();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* tmp = NULL;
	char myCurveDate[11];

	CCString msg (" ");

	try
	{
		Local_XLDATE2ARMDATE(date,myCurveDate);

		tmp = (char*) Ccy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			free(tmp);
		}
		tmp = NULL;

		mktData = CreateARMVectorFromVECTOR(rate);

		Bonds = new ARM_Container();

		for (int i = 0; i < szBd; i++)
		{
			Bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondsId[i]);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Bond, ARM_BOND) == 0)
			{
				result.setMsg ("ARM_ERR: Bond is not of a good type");
				return ARM_KO;
			}

			Bonds->Append(Bond);
		}

		Yields = CreateARMVectorFromVECTOR(yields);

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "X");

		for(int j = 0; j < rateSize; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		
		createdZeroCurve = new ARM_ZeroLInterpol((ARM_Date) myCurveDate,
												 psMatu,
												 mktData,
												 Bonds,
												 Yields,
												 MMVsFut,
												 0,
												 aCcy);

		if (aCcy)
			delete aCcy;
		aCcy = NULL;


		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			oldZeroCurve = (ARM_ZeroLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldZeroCurve, ARM_ZERO_LIN_INTERPOL) == 1)
			{
				if (oldZeroCurve)
				{
					delete oldZeroCurve;
					oldZeroCurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmp)
			delete tmp;
		tmp = NULL;

		if (createdZeroCurve)
			delete createdZeroCurve;
		createdZeroCurve = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_GenerateBasisAdjCurve(long Crv1Id,
									long BSCrv1Id,
									long Crv2Id,
									long BSCrv2Id,
									long flagInputAsSprds,
									long flagRetSprds,
									VECTOR<CCString>& matu,
									ARM_result& result,
									long objId)
{	
	ARM_ZeroCurve* Crv1=NULL;
	ARM_ZeroCurve* BSCrv1=NULL;
	ARM_ZeroCurve* Crv2=NULL;
	ARM_ZeroCurve* BSCrv2=NULL;
	
	ARM_ZeroCurve* createBasisCrv = NULL;

	ARM_ZeroCurve* oldBasisCrv = NULL;
	long basisCrvId;

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		Crv1 = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Crv1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Crv1, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR:  curve1 is not of a good type");
			return ARM_KO;
		}

		BSCrv1 = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(BSCrv1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(BSCrv1, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR:  BS curve1 is not of a good type");
			return ARM_KO;
		}

		Crv2 = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Crv2Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Crv2, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: curve2 is not of a good type");
			return ARM_KO;
		}

		BSCrv2 = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(BSCrv2Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(BSCrv2, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: BS curve2 is not of a good type");
			return ARM_KO;
		}

		int real_size = matu.size ();

		ARM_CRV_TERMS psMatu;

		if ( real_size > 0 )
		{
			for (int j = 0; j < ARM_NB_TERMS; j++)
				sprintf(psMatu[j],"%s", "X");

			for (int j = 0; j < real_size; j++)
			{
				sprintf(psMatu[j], "%s", (const char *) matu[j]);
			}
		}

        createBasisCrv = GenerateTwoCcyBSAdjusted(Crv1, // Exp: EUR: The reference CCY
                                                  BSCrv1,
                                                  Crv2, // Exp: JPY
                                                  BSCrv2,
                                                  real_size,
                                                  psMatu,
												  flagInputAsSprds,
                                                  flagRetSprds);

		if (createBasisCrv == NULL)
		{
			result.setMsg ("ARM_ERR: Returned Spreaded Curve is null!!");
			return ARM_KO;			
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			basisCrvId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createBasisCrv);

			if (basisCrvId == RET_KO)
			{
				if (createBasisCrv)
					delete createBasisCrv;
				createBasisCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(basisCrvId);

			return ARM_OK;
		}
		else
		{
			oldBasisCrv = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(oldBasisCrv, ARM_ZERO_CURVE) == 1)
			{
				if (oldBasisCrv)
				{
					delete oldBasisCrv;
					oldBasisCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createBasisCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (createBasisCrv)
					delete createBasisCrv;
				createBasisCrv = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

    catch(Exception& x)
    {
		if (createBasisCrv)
			delete createBasisCrv;
		createBasisCrv = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GenForwardYield (	long idCurve,
								double matu1,
								double matu2,
								long isSwapRate,
								long decompFreqId,
								long daycountId,
								ARM_result& result )
{
	double dResult;
	ARM_ZeroCurve* zc=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = ApproximatedForwardRate(matu1,matu2,zc,isSwapRate,daycountId);

		if ( decompFreqId != -9999 )
		{
			dResult = FromRateToRate(dResult, 1.0, K_COMP_ANNUAL, decompFreqId);
		}

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_OldZcCurve(
	const long& zcId,
	const double& asOf,
	ARM_result& result,
	long objId)
{
	ARM_ZeroCurve* zc = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date asOfDate;
		Local_XLDATE2ARMDATE(asOf,asOfDate);

		ARM_OldZcCurve* oldZcCurve = new ARM_OldZcCurve(asOfDate,zc);

		// assign object
		if ( !assignObject( oldZcCurve, result, objId ) )
		{
			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_zcswapsplsum (double date,
							 VECTOR<CCString>& matu,
							 VECTOR<double>& rate,
							 long mmVsFut,
							 long swapVsFut,
							 long raw,
							 long interp,
							 const CCString& ccy,
							 ARM_result& result,
							 long objId)
{
	long curveId;

	if(matu.size () != rate.size ())
	{
		result.setMsg ("ARM_ERR: maturities and rates must have same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	long nblines = matu.size ();

	ARM_ZeroSpliSum* createdZeroCurve = NULL;
	ARM_ZeroSpliSum* prevZcSplSum = NULL;
	ARM_Currency* aCcy = NULL;
    ARM_Vector* mktData = NULL;

	char* tmp = NULL;
	char myCurveDate[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(date,myCurveDate);

		tmp = (char*) ccy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			free(tmp);
		}
		tmp = NULL;

		ARM_CRV_TERMS psMatu;

		for (int j = 0; j < ARM_NB_TERMS; j++)
			sprintf(psMatu[j], "%s", "X");

		for(int j = 0; j < nblines; j++)
		{
			sprintf(psMatu[j], "%s", (const char *) matu[j]);
		}

		mktData = CreateARMVectorFromVECTOR(rate);

		createdZeroCurve = new ARM_ZeroSpliSum((ARM_Date) myCurveDate, psMatu,
												mktData, mmVsFut, swapVsFut,
												raw, interp, 0, aCcy);

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (createdZeroCurve == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve);

			if (curveId == RET_KO)
			{
				if (createdZeroCurve)
					delete createdZeroCurve;
				createdZeroCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevZcSplSum = (ARM_ZeroSpliSum *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevZcSplSum, ARM_ZERO_SPLSUM) == 1)
			{
				if (prevZcSplSum)
				{
					delete prevZcSplSum;
					prevZcSplSum = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZeroCurve, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmp)
			delete tmp;
		tmp = NULL;

		if (createdZeroCurve)
			delete createdZeroCurve;
		createdZeroCurve = NULL;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		ARM_RESULT();
    }
}

long ARMLOCAL_FixingSched(const double&		asOf,
						  vector<string>	LiborKeys,
						  vector<string>	FXKeys,
						  vector<long>		LiborCurveId,
						  vector<long>		FXCurveId,
						  ARM_result&		result,
						  long				objId)
{
	CCString msg ("");

	try
	{
		char* sDate = new char[11];

		Local_XLDATE2ARMDATE(asOf,sDate);
		ARM_Date myDate(sDate);

		if (sDate)
			delete [] sDate;
		sDate = NULL;
	

		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		map <string,ARM_Curve *> LiborFixingMap;
		map <string,ARM_Curve *> FXFixingMap;

		ARM_FixingSched * FixingSched = new ARM_FixingSched(myDate);

		for (int i = 0; i < LiborCurveId.size(); i++ )
		{
			ARM_Curve * curve = (ARM_Curve *) LOCAL_PERSISTENT_OBJECTS->GetObject(LiborCurveId[i]);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(curve, ARM_GENERIC_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			FixingSched->AddLiborFixing(LiborKeys[i],curve);
		}

		for (int i = 0; i < FXCurveId.size(); i++ )
		{
			ARM_Curve * curve = (ARM_Curve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FXCurveId[i]);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(curve, ARM_GENERIC_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			FixingSched->AddFxFixing(FXKeys[i],curve);
		}

		long FixingSchedId;

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			FixingSchedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FixingSched);

			if (FixingSchedId == RET_KO)
			{
				if (FixingSched)
					delete FixingSched;
				FixingSched = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(FixingSchedId);

			return ARM_OK;
		}
		else
		{
			ARM_FixingSched * previous_FixingSched = (ARM_FixingSched *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(previous_FixingSched, ARM_FIXING_SCHED) == 1)
			{
				if (previous_FixingSched)
				{
					delete previous_FixingSched;
					previous_FixingSched = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FixingSched, objId);

				return ARM_OK;
			}
			else
			{
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

	return ARM_OK;
}


long ARMLOCAL_GetFixingSchedFromSummit(	VECTOR<string>&		ListOfKeys,
										const double&		AsOf,
										const CCString&		Source,
										long&				DateStripId,
										ARM_result&			result,
										long				objId)
{
	CCString			msg ("");
	ARM_result			C_result;
	ARM_FixingSched*	FixingSched = NULL;

	try
	{
		char* sDate = new char[11];

		Local_XLDATE2ARMDATE(AsOf,sDate);
		ARM_Date myDate(sDate);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		FixingSched = new ARM_FixingSched(myDate);

		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		ARM_ResetManager* resetMgr;

		for (int i = 0; i < ListOfKeys.size(); i++)
		{
			if ( ListOfKeys[i].size() == 7 )
			{
				resetMgr = ARMLOCAL_ParseXMLForResetMgr(myDate,
														"FX",
														"MO",
														CCString(ListOfKeys[i].c_str()),
														K_NO,
														""
													 );
			}
			else
			{
				resetMgr = ARMLOCAL_ParseXMLForResetMgr(myDate,
														CCString(ListOfKeys[i].substr(0,5).c_str()),
														Source,
														CCString(ListOfKeys[i].substr(6,3).c_str()),
														K_NO,
														CCString(ListOfKeys[i].substr(10,(ListOfKeys[i].size() - 10)).c_str())
													 );
			}

			ARM_DateStrip * DS	= (ARM_DateStrip *) LOCAL_PERSISTENT_OBJECTS->GetObject(DateStripId);
			
			std::vector<double> * ResetDates		= DS->GetResetDates();
			std::vector<double> * PaymentDates	= DS->GetPaymentDates();

			if ( myDate.GetJulian() < ARM_Date(ResetDates->Elt(0)).GetJulian() )
			{
				result.setMsg ("ARM_ERR: Date Strip must start forward. No past fixing");
				return ARM_KO;
			}
			
			int			j = 0;
			ARM_Date	ResetDateResearched;
	
			while ( ARM_Date(PaymentDates->Elt(j)).GetJulian() < myDate.GetJulian() ) { j++; }

			if ( myDate.GetJulian() > ARM_Date(	ResetDates->Elt(j)).GetJulian() )
				ResetDateResearched = ARM_Date(	ResetDates->Elt(j)	);
			else
			{
				result.setMsg ("ARM_ERR: AsOf inferior to the Reset date. This case should never happen");
				return ARM_KO;
			}

			double FixingResearched = resetMgr->GetReset(ResetDateResearched.GetJulian());

			ARM_Curve * curve;

			curve = new ARM_Curve(	std::vector<double>( 1, (ResetDateResearched - myDate) ),	/* abscisses */
									std::vector<double>( 1, FixingResearched ),				/* ordinates */
									new ARM_StepUpRightOpenCstExtrapolDble);			/* interpolator */

			if ( ListOfKeys[i].size() == 7 )
			{
				FixingSched->AddFxFixing(ListOfKeys[i],curve);
			}			

			else
			{
				FixingSched->AddLiborFixing(ListOfKeys[i],curve);
			}
		}

		if ( FixingSched->GetNbFixings() == 0 )
		{
			delete FixingSched;
			FixingSched = NULL;
		}

		long FixingSchedId;

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			FixingSchedId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FixingSched);

			if (FixingSchedId == RET_KO)
			{
				if (FixingSched)
					delete FixingSched;
				FixingSched = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(FixingSchedId);

			return ARM_OK;
		}
		else
		{
			ARM_FixingSched * previous_FixingSched = (ARM_FixingSched *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(previous_FixingSched, ARM_FIXING_SCHED) == 1)
			{
				if (previous_FixingSched)
				{
					delete previous_FixingSched;
					previous_FixingSched = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FixingSched, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();

		if (FixingSched)
			delete FixingSched;
		FixingSched = NULL;
	}

	return ARM_OK;


}


/*---- End Of File ----*/


