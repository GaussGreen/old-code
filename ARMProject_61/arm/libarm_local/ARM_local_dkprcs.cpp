#pragma warning(disable : 4786)
#pragma warning(disable : 4541)
#pragma warning(disable : 4250)


#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <crv\zerocurv.h>
#include <crv\zeroflat.h>
#include <crv\volflat.h>

#include <mod\ycmodel.h>
#include <mod\bsmodel.h>
#include <mod\xg2ycmod.h>
#include <mod\mcfnhw1fsv.h>
#include <mod\frmtree.h>
#include <mod\dftreehwsigvar.h>
#include <mod\bootstrapcalibration.h>
#include <mod\smc_frm.h>
#include <mod\irtbk.h>
#include <mod\smiledmcrnldc.h>
#include <mod\irthw2.h>
#include <mod\irthwsc.h>
#include <mod\irthwsigvar.h>
#include <mod\xfhwfx.h>
#include <mod\dfirthw.h>
#include <mod\bmc_frm.h>
#include <mod\bssmiled.h>
#include <mod\crrtree.h>
#include <mod\bscorrmodel.h>
#include <mod3F\DK_prcs.h>
#include <mod3F\tree3f.h>
#include <mod\xbsfx.h>
#include <glob\armdef.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\libarm_local\ARM_local_dkprcs.h>



long ARMLOCAL_PRCS3F_Lattice_HWVFDK_Pricing(VECTOR<double>& dLatticeGeometryDataIn,
		 								    double dNumTimeLinesBeforeFirstNotice,
										    double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
										    double dNumTimeLines, 
										    double evalDate, 
										    long dBaseYieldCurveId,
										    long dForeignYieldCurveId,
											long dBaseRatesNoBasisCurveId,
											long dForeignRatesNoBasisCurveId,
										    long volSwopBaseId,
                                            long volSwopForeignId,
                                            long volFxId,
										    VECTOR<double>& dNoticeDatesIn,
										    double dStrike,
										    double dType,
										    double dOptionExpiryIn, 
										    double dMeanReversionBase,
										    double dMeanReversionForeign,  
										    double dSpotFX,
										    long dBaseForeignCorrelationId,
										    long dBaseSpotFXCorrelationId,
										    long dForeignSpotFXCorrelationId,
										    double dProductModelCodeId,
										    double dX1Limit,
										    double dX2Limit,
										    double dX3Limit,
										    double dI1Limit,
										    double dI2Limit,
										    double dI3Limit,
										    double dADPLimit,
										    double dOptimal,
										    double dTimeBoost,
										    double dDeltaFlag,
										    double dStringModel,
										    long nbCol,
										    VECTOR<double>& dBoosterDataIn,
										    ARM_result& result)
{
    CCString msg ("");

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 ) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");		
		return(NULL);
	}

	char tmpDate[11];

    Local_XLDATE2ARMDATE(evalDate, tmpDate);

    ARM_Date AsOfDate(tmpDate);


    Local_XLDATE2ARMDATE(dOptionExpiryIn, tmpDate);
    
	ARM_Date dOptionExpiry(tmpDate);

    // Vectors and Matrixes

	ARM_Vector* dLatticeGeometryData = CreateARMVectorFromVECTOR(dLatticeGeometryDataIn);

    ARM_Vector* dNoticeDates         = CreateARMVectorFromVECTOR(dNoticeDatesIn);


    int i, j;

    int lsz = dLatticeGeometryData->GetSize();

    for (j = 1; j < lsz; j++)
	{
		Local_XLDATE2ARMDATE(dLatticeGeometryData->Elt(j), tmpDate);

        ARM_Date buffDate;
        
        try
        {
            ARM_Date curDate(tmpDate);

            buffDate = curDate;
        }

        catch(Exception& x)
        {
		    ARM_RESULT();
        }

        dLatticeGeometryData->Elt(j) = buffDate.GetJulian();
	}

	int dsz = dNoticeDates->GetSize();

	for (j = 0; j < dsz; j++)
	{
		Local_XLDATE2ARMDATE(dNoticeDates->Elt(j), tmpDate);

        ARM_Date buffDate;

        try
        {
           ARM_Date curDate(tmpDate);

           buffDate = curDate;
        }

        catch(Exception& x)
        {
		    ARM_RESULT();
        }

        dNoticeDates->Elt(j) = buffDate.GetJulian();
	}

	ARM_Vector* dBoosterData = CreateARMVectorFromVECTOR(dBoosterDataIn);

    
	long nbLines = dBoosterDataIn.size();

	long matrixNbLines = nbLines/nbCol;

    ARM_Matrix dataMatrix(matrixNbLines, nbCol);

	for (i = 0; i < matrixNbLines; i++)
	{
		for (j = 0; j < nbCol; j++)
		{
            if (( 0 <= j ) && ( j <= 5 ))
            {
               Local_XLDATE2ARMDATE(dBoosterDataIn[i*nbCol+j], tmpDate);

               ARM_Date		buffDate;

               try
               {
                   ARM_Date curDate(tmpDate);
               
                   buffDate = curDate;
               }

               catch(Exception& x)
               {
		            ARM_RESULT();
               }

               dataMatrix.Elt(i, j) = buffDate.GetJulian();
            }
            else
            {
               dataMatrix.Elt(i, j) = dBoosterDataIn[i*nbCol+j];
            }
		}
	}


	if (dBoosterData)
	   delete dBoosterData;


    // Yield curves
	ARM_ZeroCurve* dBaseYieldCurve = NULL;
	ARM_ZeroCurve* dForeignYieldCurve = NULL;

	ARM_ZeroCurve* dBaseRatesNoBasisCurve;
    ARM_ZeroCurve* dForeignRatesNoBasisCurve;
	
	// Volatity objects	
	ARM_VolCurve* volSwopBase = NULL;
    ARM_VolCurve* volSwopForeign = NULL;	
	ARM_VolCurve* volFx = NULL;
	
	// Correlations
    ARM_VolCurve* dBaseForeignCorrelation = NULL;	
    ARM_VolCurve* dBaseSpotFXCorrelation = NULL;
    ARM_VolCurve* dForeignSpotFXCorrelation = NULL;

	ARM_Vector* res = NULL;

	try
	{
		dBaseYieldCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBaseYieldCurveId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBaseYieldCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Base Yield Curve is expected");			
		   return(NULL);
		}

		dForeignYieldCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dForeignYieldCurveId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dForeignYieldCurve, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign Yield Curve is expected");			
			return(NULL);
		}

		dBaseRatesNoBasisCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBaseRatesNoBasisCurveId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBaseRatesNoBasisCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg("ARM_ERR: Base no basis Yield Curve is expected");
		   
		   return(NULL);
		}

		dForeignRatesNoBasisCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dForeignRatesNoBasisCurveId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dForeignRatesNoBasisCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg("ARM_ERR: Base no basis Yield Curve is expected");
		   
		   return(NULL);
		}

		volSwopBase = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volSwopBaseId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volSwopBase, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. Swopt. Base is not a Vol Curve");
			return(NULL);
		}

		volSwopForeign = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volSwopForeignId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volSwopForeign, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. Swopt. Foreign is not a Vol Curve");
			return(NULL);
		}

		volFx = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volFxId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volFx, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. FX is not a Vol Curve Object");	
			return(NULL);
		}

		dBaseForeignCorrelation = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBaseForeignCorrelationId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBaseForeignCorrelation, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Base Foreign correlation is not a Vol Curve");
			return(NULL);
		}

		dBaseSpotFXCorrelation = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBaseSpotFXCorrelationId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBaseSpotFXCorrelation, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Base Spot FX correlation is not a Vol Curve");			
			return(NULL);
		}
		
		dForeignSpotFXCorrelation = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dForeignSpotFXCorrelationId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dForeignSpotFXCorrelation, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign Spot FX correlation is not a Vol Curve");			
			return(NULL);
		}

		res = PRCS3F_Lattice_HWVFDK_Pricing(dLatticeGeometryData,
		  							         dNumTimeLinesBeforeFirstNotice,
									         dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
									         dNumTimeLines, 
									         AsOfDate, 
									         dBaseYieldCurve,
									         dForeignYieldCurve,
                                             dBaseRatesNoBasisCurve,
                                             dForeignRatesNoBasisCurve,
									         volSwopBase,
									         volSwopForeign,
									         volFx,
											 dNoticeDates,
										 	 dStrike,
											 dType,
											 dOptionExpiry, 
											 dMeanReversionBase, 
											 dMeanReversionForeign,  
											 dSpotFX,
											 dBaseForeignCorrelation,
											 dBaseSpotFXCorrelation,
											 dForeignSpotFXCorrelation,
											 dProductModelCodeId,
											 dX1Limit,
											 dX2Limit,
											 dX3Limit,
											 dI1Limit,
											 dI2Limit,
											 dI3Limit,
											 dADPLimit,
											 dOptimal,
											 dTimeBoost,
											 dDeltaFlag,
											 0.0,
                                             0,
											 0.0,
											 0.0,
											 0,
											 0.0,
											 dStringModel,
											 &dataMatrix);
	}
	
	catch(Exception& x)
    {
        if (dLatticeGeometryData)
	       delete dLatticeGeometryData;

	    if (dNoticeDates)
           delete dNoticeDates;

		ARM_RESULT();
    }

	catch(...)
    {
        if (dLatticeGeometryData)
	       delete dLatticeGeometryData;

	    if (dNoticeDates)
           delete dNoticeDates;	
		
		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal PRCS pricing failure!?");

        ARM_RESULT();
	}

	if ( res == NULL )
	{
       if (dLatticeGeometryData)
	      delete dLatticeGeometryData;

	    if (dNoticeDates)
           delete dNoticeDates;	
		
		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
						"Internal PRCS pricing failure!?");

		ARM_RESULT();	
	}

    if (dLatticeGeometryData)
	   delete dLatticeGeometryData;

	if (dNoticeDates)
       delete dNoticeDates;

    // Build the Result Structure

	long sz = res->GetSize();
	
	result.setLong(sz);

	for (int k = 0; k < sz; k++)
	{
	    result.setArray(res->Elt(k), k);
	}

	if (res)
	   delete res;

	return(ARM_OK);
}


long ARMLOCAL_TREE3FACT (double asof,
						 long xbsfxId,
						 long volSwopBaseId,
						 long volSwopForeignId,
						 long dBaseForeignCorrelationId,
						 VECTOR<double>& dLatticeGeometryDataIn,
						 double dNumTimeLinesBeforeFirstNotice,
						 double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
						 double dNumTimeLines,
						 double dMeanReversionBase,
						 double dMeanReversionForeign,
						 double dX1Limit,
						 double dX2Limit,
						 double dX3Limit,
						 double dI1Limit,
						 double dI2Limit,
						 double dI3Limit,
						 double dADPLimit,
						 double dOptimal,
						 double dTimeBoost,
						 double dDeltaFlag,
						 double dSmoothingValue,
						 long calcProbaSurvOrNotId,
						 double QBaseSmile,
				         double QForeignSmile,
						 long spotVolOrFwdVol,
						 double histoVolLongTerm,
						 long convInVolFwd,
						 long calibBasisIncluded,
						 ARM_result& result,
						 long objId)
{
	long modId;

	ARM_Tree3F* newModel = NULL;
	ARM_Tree3F* oldModel = NULL;

	ARM_DFBSModel* xbsfx = NULL;
	ARM_VolCurve* volSwopBase = NULL;
    ARM_VolCurve* volSwopForeign = NULL;	
    ARM_VolCurve* dBaseForeignCorrelation = NULL;	

    // Vectors and Matrixes
	ARM_Vector* dLatticeGeometryData = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char tmpDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(asof, tmpDate);
		ARM_Date AsOfDate(tmpDate);

		xbsfx = (ARM_DFBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(xbsfxId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xbsfx, ARM_DFBSMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: xbsfx model is not of a good type");
			return ARM_KO;
		}

		volSwopBase = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volSwopBaseId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volSwopBase, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. Swopt. Base is not a Vol Curve");
			return(ARM_KO);
		}

		volSwopForeign = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volSwopForeignId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volSwopForeign, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. Swopt. Foreign is not a Vol Curve");
			return(ARM_KO);
		}

		dBaseForeignCorrelation = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBaseForeignCorrelationId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBaseForeignCorrelation, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Base Foreign correlation is not a Vol Curve");
			return(ARM_KO);
		}

		dLatticeGeometryData = CreateARMVectorFromVECTOR(dLatticeGeometryDataIn);

		int dsz;
		int j;

		dsz = dLatticeGeometryData->GetSize();

		for (j = 1; j < dsz; j++)
		{
			Local_XLDATE2ARMDATE(dLatticeGeometryData->Elt(j), tmpDate);
			ARM_Date buffDate;

			try
			{
			   ARM_Date curDate(tmpDate);
			   buffDate = curDate;
			}

			catch(Exception& x)
			{
				ARM_RESULT();
			}

			dLatticeGeometryData->Elt(j) = buffDate.GetJulian();
		}

		newModel = new ARM_Tree3F(AsOfDate,
								  xbsfx,
								  volSwopBase,
								  volSwopForeign,
								  dBaseForeignCorrelation,
								  dLatticeGeometryData,
								  dNumTimeLinesBeforeFirstNotice,
								  dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
								  dNumTimeLines,
								  dMeanReversionBase,
								  dMeanReversionForeign,
								  dX1Limit,
								  dX2Limit,
								  dX3Limit,
								  dI1Limit,
								  dI2Limit,
								  dI3Limit,
								  dADPLimit,
								  dOptimal,
								  dTimeBoost,
								  dDeltaFlag,
								  dSmoothingValue,
								  calcProbaSurvOrNotId,
								  QBaseSmile,
				                  QForeignSmile,
								  spotVolOrFwdVol,
								  histoVolLongTerm,
								  convInVolFwd,
								  0,
								  calibBasisIncluded);

		if (dLatticeGeometryData)
	       delete dLatticeGeometryData;

		if (newModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newModel);

			if (modId == RET_KO)
			{
				if (newModel)
					delete newModel;
				newModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldModel = (ARM_Tree3F *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldModel, ARM_TREE3F) == 1)
			{
				if (oldModel)
				{
					delete oldModel;
					oldModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newModel, objId);

				return ARM_OK;
			}
			else
			{
				if (newModel)
					delete newModel;
				newModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
        if (dLatticeGeometryData)
	       delete dLatticeGeometryData;

		ARM_RESULT();
    }

	catch(...)
    {
        if (dLatticeGeometryData)
	       delete dLatticeGeometryData;

		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal PRCS pricing failure!?");

        ARM_RESULT();
	}

}


long  ARMLOCAL_Bootstrapping_VFDK_HW1To3F ( long volId,
                                            long zcId,
                                            VECTOR<double>& noticeDatesIn,
                                            VECTOR<double>& swapStartDates,
                                            VECTOR<double>& swapEndDates,
                                            VECTOR<double>& HW3FParameters,
                                            double observationDate,
                                            ARM_result& result)
{
    ARM_VolCurve* volCurve = NULL;
    ARM_ZeroCurve* zeroCurve = NULL;	
    ARM_Vector* noticeDates = NULL;
    ARM_Vector* startDates = NULL;
    ARM_Vector* endDates = NULL;
	ARM_Vector* HW3FParam = NULL;

    ARM_Vector* sigma = NULL;

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char obsDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(observationDate, obsDate);
		ARM_Date ARMObsDate(obsDate);

		volCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurve, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. is not a Vol Curve");
			return(ARM_KO);
		}

        zeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) zcId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zeroCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Zero Curve is expected");			
		   return(NULL);
		}
        
        char tmpDate[11];
        ARM_Date buffDate;
        
        noticeDates = CreateARMVectorFromVECTOR(noticeDatesIn);
        int dsz = noticeDates->GetSize();
        int j;
	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(noticeDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
		        ARM_RESULT();
            }

            noticeDates->Elt(j) = buffDate.GetJulian();
	    }
        startDates = CreateARMVectorFromVECTOR(swapStartDates);
        dsz = startDates->GetSize();

	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(startDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
		        ARM_RESULT();
            }

            startDates->Elt(j) = buffDate.GetJulian();
	    }
        endDates = CreateARMVectorFromVECTOR(swapEndDates);
        dsz = endDates->GetSize();

	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(endDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
		        ARM_RESULT();
            }

            endDates->Elt(j) = buffDate.GetJulian();
	    }
        HW3FParam = CreateARMVectorFromVECTOR(HW3FParameters);

         	
        

       sigma = PRCS3F_Bootstrapping_VFDK_HW1To3F(volCurve,
		 										 zeroCurve,
												 noticeDates,
                                                 startDates,
                                                 endDates,
                                                 HW3FParam,
                                                 ARMObsDate);

	    delete noticeDates;
        delete startDates;
        delete endDates;
        delete HW3FParam;

        // Build the Result Structure

	    long sz = sigma->GetSize();
	    
	    result.setLong(sz);

	    for (int k = 0; k < sz; k++)
	    {
	        result.setArray(sigma->Elt(k), k);
	    }

	    if (sigma)
	       delete sigma;


		return ARM_OK;
	}


	catch(...)
    {
        
	    delete noticeDates;
        delete startDates;
        delete endDates;
        delete HW3FParam;
        delete sigma;

		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal Bootstrapping_VFDK_HW1To3F failure!?");

        ARM_RESULT();
	}


}

long  ARMLOCAL_SwaptionPrice_VFDK_HW1To3F ( double dSwaptionExpiryInYears,
                                            double dSwaptionTenorInYears,
                                            double dNoticePeriodInDays,
                                            double dStrike,
                                            double dCallPut,
                                            long zcId,
                                            VECTOR<double>& dNoticeDates,
                                            VECTOR<double>& dSigma,
                                            VECTOR<double>& HW3FParameters,
                                            double observationDate,
                                            ARM_result& result)
{
    ARM_ZeroCurve* zeroCurve = NULL;	
    ARM_Vector* noticeDates = NULL;
    ARM_Vector* sigma = NULL;
    ARM_Vector* HW3FParams = NULL;

    double Price;

 

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char obsDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(observationDate, obsDate);
		ARM_Date ARMObsDate(obsDate);

        zeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) zcId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zeroCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Zero Curve is expected");			
		   return(NULL);
		}
        
        char tmpDate[11];
        ARM_Date buffDate;
        
        noticeDates = CreateARMVectorFromVECTOR(dNoticeDates);
        int dsz = noticeDates->GetSize();
        int j;
	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(noticeDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
		        ARM_RESULT();
            }

            noticeDates->Elt(j) = buffDate.GetJulian();
	    }
        

       HW3FParams = CreateARMVectorFromVECTOR(HW3FParameters);
       sigma = CreateARMVectorFromVECTOR(dSigma);

         	
        

       Price = PRCS3F_SwaptionPrice_VFDK_HW1To3F(   dSwaptionExpiryInYears,
                                                    dSwaptionTenorInYears,
                                                    dNoticePeriodInDays,
                                                    dStrike,
                                                    dCallPut,
		 										    zeroCurve,
												    noticeDates,
                                                    sigma,
                                                    HW3FParams,
                                                    ARMObsDate);

	    delete noticeDates;
        delete HW3FParams;
        delete sigma;


        // Build the Result Structure

	    
	    result.setDouble(Price);

		return ARM_OK;
	}


	catch(...)
    {
        
	    delete noticeDates;
        delete sigma;
        delete HW3FParams;

		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal SwaptionPrice_VFDK_HW1To3F failure!?");

        ARM_RESULT();
	}


}


long  ARMLOCAL_ImpliedFwdCorrelation_VFDK_HW1To3F ( double dSwaptionExpiryInYears,
                                                    double dSwaptionTenorInYears,
                                                    double dSwaptionTenor2InYears,
                                                    double dNoticePeriodInDays,
                                                    long zcId,
                                                    VECTOR<double>& dNoticeDates,
                                                    VECTOR<double>& dSigma,
                                                    VECTOR<double>& HW3FParameters,
                                                    double observationDate,
                                                    ARM_result& result)
{
    ARM_ZeroCurve* zeroCurve = NULL;	
    ARM_Vector* noticeDates = NULL;
    ARM_Vector* sigma = NULL;
    ARM_Vector* HW3FParams = NULL;

    double Price;

 

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char obsDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(observationDate, obsDate);
		ARM_Date ARMObsDate(obsDate);

        zeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) zcId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zeroCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Zero Curve is expected");			
		   return(NULL);
		}
        
        char tmpDate[11];
        ARM_Date buffDate;
        
        noticeDates = CreateARMVectorFromVECTOR(dNoticeDates);
        int dsz = noticeDates->GetSize();
        int j;
	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(noticeDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
		        ARM_RESULT();
            }

            noticeDates->Elt(j) = buffDate.GetJulian();
	    }
        

       HW3FParams = CreateARMVectorFromVECTOR(HW3FParameters);
       sigma = CreateARMVectorFromVECTOR(dSigma);

         	
        

       Price = PRCS3F_ImpliedFwdCorrelation_VFDK_HW1To3F(   dSwaptionExpiryInYears,
                                                            dSwaptionTenorInYears,
                                                            dSwaptionTenor2InYears,
                                                            dNoticePeriodInDays,
                                                            zeroCurve,
												            noticeDates,
                                                            sigma,
                                                            HW3FParams,
                                                            ARMObsDate);

	    delete noticeDates;
        delete HW3FParams;
        delete sigma;


        // Build the Result Structure

	    
	    result.setDouble(Price);

		return ARM_OK;
	}


	catch(...)
    {
        
	    delete noticeDates;
        delete sigma;
        delete HW3FParams;

		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal ImpliedFwdCorrelation_VFDK_HW1To3F failure!?");

        ARM_RESULT();
	}


}



long  ARMLOCAL_HW3F_CALIBRATION(double dSpotDate,
                                long zcId,
                                VECTOR<double>& HW3FParametersIn,
                                long volId,
                                long correlationId,
                                long volWeightsId,
                                long correlationWeightsId,
                                VECTOR<double>& dNoticeDates,
                                VECTOR<double>& dSwapStartDates,
                                VECTOR<double>& dSwapEndDates,
                                ARM_result& result)
{
    ARM_ZeroCurve* zeroCurve = NULL;
    ARM_VolCurve* volCurve = NULL;
    ARM_VolCurve* volWeightCurve = NULL;
    ARM_VolCurve* correlationCurve = NULL;
    ARM_VolCurve* correlationWeightCurve = NULL;

    ARM_Vector* noticeDates    = NULL;
    ARM_Vector* swapStartDates = NULL;
    ARM_Vector* swapEndDates   = NULL;
    
	ARM_Vector* HW3FParamsIn  = NULL;
    ARM_Vector* HW3FParamsOut = NULL;

    

    if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
	   return(ARM_KO);
	}

	char spotDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(dSpotDate, spotDate);
		ARM_Date ARMSpotDate(spotDate);

        zeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) zcId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zeroCurve, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: Zero Curve is expected");			
		   return(NULL);
		}
        
        volCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurve, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Vol. is not a Vol Curve");
			return(ARM_KO);
		}

        if ( volWeightsId != ARM_NULL_OBJECT )
		{
           volWeightCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) volWeightsId);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volWeightCurve, ARM_VOL_CURVE) == 0 )
		   {
			  result.setMsg("ARM_ERR: VolWeights is not a Vol Curve");
			
			  return(ARM_KO);
		   }
		}
		else
		{
			// Create a default object vol initialized to 1

			int nbLines = volCurve->GetVolatilities()->GetNumLines();
			int nbCol   = volCurve->GetVolatilities()->GetNumCols();

			ARM_Vector* yearTerms  = new ARM_Vector(volCurve->GetExpiryTerms());
			ARM_Vector* maturities = new ARM_Vector(((ARM_VolLInterpol *) volCurve)->GetStrikes());

			ARM_Matrix* elements = new ARM_Matrix(nbLines, nbCol, 1.0);

            volWeightCurve = new ARM_VolLInterpol(ARMSpotDate, yearTerms, maturities, elements);
		}

		if ( correlationId != ARM_NULL_OBJECT )
		{
           correlationCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) correlationId);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlationCurve, ARM_VOL_CURVE) == 0 )
		   {
              if ( volWeightsId == ARM_NULL_OBJECT )
			  {				   
			     delete volWeightCurve;
			  }
			
		      result.setMsg("ARM_ERR: CorrelationCurve is not a Vol Curve");
			
		      return(ARM_KO);
		   }
		}
		else
		{
			// Create a default object vol initialized to 0.0

			int nbLines = 1;
			int nbCol   = 1;

			ARM_Vector* yearTerms  = new ARM_Vector(1, 0.0);
			ARM_Vector* maturities = new ARM_Vector(1, 0.0);

			ARM_Matrix* elements = new ARM_Matrix(nbLines, nbCol, 0.0);

            correlationCurve = new ARM_VolLInterpol(ARMSpotDate, yearTerms, maturities, elements);
		}

		if ( correlationWeightsId != ARM_NULL_OBJECT )
		{
           correlationWeightCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) correlationWeightsId);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlationWeightCurve, ARM_VOL_CURVE) == 0 )
		   {
		      if ( volWeightsId == ARM_NULL_OBJECT )
			  {				   
			     delete volWeightCurve;
			  }

			  if ( correlationId == ARM_NULL_OBJECT )
			  {
				 delete correlationCurve;
			  }

		      result.setMsg("ARM_ERR: CorrelationWeights is not a Vol Curve");
			
		      return(ARM_KO);
		   }
		}
		else
		{
			// Create a default object vol initialized to 0.0

			int nbLines = 1;
			int nbCol   = 1;

			ARM_Vector* yearTerms  = new ARM_Vector(1, 0.0);
			ARM_Vector* maturities = new ARM_Vector(1, 0.0);

			ARM_Matrix* elements = new ARM_Matrix(nbLines, nbCol, 0.0);

            correlationWeightCurve = new ARM_VolLInterpol(ARMSpotDate, yearTerms, maturities, elements);
		}
        
		char tmpDate[11];
        ARM_Date buffDate;
        
        noticeDates = CreateARMVectorFromVECTOR(dNoticeDates);
        int dsz = noticeDates->GetSize();
        int j;
	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(noticeDates->Elt(j), tmpDate);

            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
                if ( volWeightsId == ARM_NULL_OBJECT )
				{				   
				   delete volWeightCurve;
				}

				if ( correlationId == ARM_NULL_OBJECT )
				{
				   delete correlationCurve;
				}

				if ( correlationWeightsId == ARM_NULL_OBJECT )
				{
				   delete correlationWeightCurve;
				}

				if (noticeDates)
				   delete noticeDates;

				ARM_RESULT();
            }

            noticeDates->Elt(j) = buffDate.GetJulian();
	    }

        swapStartDates = CreateARMVectorFromVECTOR(dSwapStartDates);
        dsz = swapStartDates->GetSize();
	    for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(swapStartDates->Elt(j), tmpDate);
            try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
				if ( volWeightsId == ARM_NULL_OBJECT )
				{
				   delete volWeightCurve;
				}

				if ( correlationId == ARM_NULL_OBJECT )
				{
				   delete correlationCurve;
				}

				if ( correlationWeightsId == ARM_NULL_OBJECT )
				{
				   delete correlationWeightCurve;
				}

				if (noticeDates)
				   delete noticeDates;

				if (swapStartDates)
				   delete swapStartDates;

		        ARM_RESULT();
            }

            swapStartDates->Elt(j) = buffDate.GetJulian();
	    }


        swapEndDates = CreateARMVectorFromVECTOR(dSwapEndDates);
        dsz = swapEndDates->GetSize();
        for (j = 0; j < dsz; j++)
	    {
		    Local_XLDATE2ARMDATE(swapEndDates->Elt(j), tmpDate);
            
			try
            {
               ARM_Date curDate(tmpDate);

               buffDate = curDate;
            }

            catch(Exception& x)
            {
				if ( volWeightsId == ARM_NULL_OBJECT )
				{
				   delete volWeightCurve;
				}

				if ( correlationId == ARM_NULL_OBJECT )
				{
				   delete correlationCurve;
				}

				if ( correlationWeightsId == ARM_NULL_OBJECT )
				{
				   delete correlationWeightCurve;
				}

				if (noticeDates)
				   delete noticeDates;

				if (swapStartDates)
				   delete swapStartDates;
				
		        ARM_RESULT();
            }

            swapEndDates->Elt(j) = buffDate.GetJulian();
	    }
        

        HW3FParamsIn = CreateARMVectorFromVECTOR(HW3FParametersIn);


        HW3FParamsOut = HW3F_CALIBRATION(ARMSpotDate,
                                         zeroCurve,
                                         HW3FParamsIn,
                                         volCurve,
                                         correlationCurve,
                                         volWeightCurve,
                                         correlationWeightCurve,
                                         noticeDates,
                                         swapStartDates,
                                         swapEndDates);

  
        delete noticeDates;
        delete swapStartDates;
        delete swapEndDates;

        delete HW3FParamsIn;       

        // Build the Result Structure

	    long sz = HW3FParamsOut->GetSize();
	    
	    result.setLong(sz);

	    for (int k = 0; k < sz; k++)
	    {
	        result.setArray(HW3FParamsOut->Elt(k), k);
	    }

		// Clean Up

	    if (HW3FParamsOut)
	       delete HW3FParamsOut;

		if ( volWeightsId == ARM_NULL_OBJECT )
		{
		   delete volWeightCurve;
		}

		return (ARM_OK);
	}


	catch(...)
    {
	    if ( volWeightsId == ARM_NULL_OBJECT )
		{
		   delete volWeightCurve;
		}

		if ( correlationId == ARM_NULL_OBJECT )
		{
		   delete correlationCurve;
		}

		if ( correlationWeightsId == ARM_NULL_OBJECT )
		{
		   delete correlationWeightCurve;
		}

	    delete noticeDates;
        delete swapStartDates;
        delete swapEndDates;

        delete HW3FParamsIn;

		Exception x(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Internal HW3F_Calibration failure!?");

        ARM_RESULT();
	}
}



/*-------------------------------------------------------------------------------------------------------*/
/*----- End Of File ----*/

