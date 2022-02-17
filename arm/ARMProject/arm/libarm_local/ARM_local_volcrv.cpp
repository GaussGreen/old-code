#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <libCCtools++\CCstring.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include "ARM_local_wrapper.h"

#include <crv\SABRVol.h>
#include <crv\volint.h>
#include <crv\volspline.h>
#include <crv\volflat.h>
#include <crv\volcube.h>
#include <crv\oldvolcurve.h>
#include <crv\hypercube.h>
#include <crv\indexindexcorrelcube.h>
#include <crv\correlmanager.h>
#include <crv\volfxspinterp.h>
#include <util\interpol.h>
#include <mod3f\tree3f.h>
#include <mod3f\DK_prcs.h>
#include <inst\powrev.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 

#include "gpmodels\Mixture_FX.h"
#include "gpbase\gplinalgconvert.h"
using ARM::ARM_ParamsMixture_Fx;
using ARM::To_ARM_GP_Vector;




//reading from file
#include <iostream>
#include <fstream>
using namespace std;
// end reading from file

/*----------------------------------------------------------------------*/





long ARMLOCAL_volcurv(const VECTOR<double>& matu,
					  const VECTOR<double>& strikes,
					  const VECTOR<double>& vols,
					  double date,
					  long strikeType,
					  long volType,
                      long interpolType,
                      const CCString& ccy,
					  const CCString& indexName,
					  long indexId,
					  ARM_result& result,
					  long objId)
{
	long volId;

	ARM_VolLInterpol* vc = NULL;
	ARM_VolLInterpol* newVolCrv = NULL;

	ARM_Vector* vMatu = NULL;
	ARM_Vector* vStrikes = NULL;
	ARM_Matrix* mVol = NULL;


    ARM_Currency currency((const char *) ccy);

	int matu_size = matu.size ();
	int strikes_size = strikes.size ();
	int vols_size = vols.size();

	if ( vols_size != (matu_size * strikes_size))
	{
	   result.setMsg ("ARM_ERR: check your volatility matrix dimension");
		
       return ARM_KO;
	}

	
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	double * pdVols = NULL;

	char startDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,startDate);
	
		pdVols = new double[vols_size];
		for(int i = 0; i < vols_size; i++)
		{
			pdVols[i] = vols [i];
		}
		
		vMatu = CreateARMVectorFromVECTOR(matu);

		vStrikes = CreateARMVectorFromVECTOR(strikes);

		mVol = new ARM_Matrix(matu_size, strikes_size, pdVols);

		if (pdVols)
		   delete [] pdVols;
		pdVols = NULL;

		if(interpolType == K_SPLINE)
			newVolCrv = new ARM_VolSplineInterpol((ARM_Date) startDate, vMatu,
									     vStrikes, mVol, strikeType, 
                                         volType,
                                         &currency);
		else
			newVolCrv = new ARM_VolLInterpol((ARM_Date) startDate, vMatu,
									     vStrikes, mVol, strikeType, 
                                         volType,
                                         &currency);


		if (indexId!=-1) 
		{
			ARM_IRIndex* index = NULL;
			LocalPersistent::get().convert(indexId,index); 
			if (index) newVolCrv->SetIndex(*index); 
		}
		// Pas de destruction : les matrices et les vecteurs ne sont pas clonés

		if (newVolCrv == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		newVolCrv->SetInterpType(int(interpolType));

		newVolCrv->SetIndexName(CCSTringToSTLString(indexName));

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv);

			if ( volId == RET_KO )
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

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
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

		if (vMatu)
		   delete vMatu;
		vMatu = NULL;

		if (vStrikes)
		   delete vStrikes;
		vStrikes = NULL;

		if (mVol)
		  delete mVol;
		mVol = NULL;

		if (newVolCrv)
			delete newVolCrv;
		newVolCrv = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_FxVolcurv(double date,
                        const VECTOR<double>& matu,
					    const VECTOR<double>& fxFwds,
                        const VECTOR<double>& pivotVols,
                        const VECTOR<double>& pivotTypes,
                        const VECTOR<double>& deltaPuts,
                        const VECTOR<double>& deltaCalls,
					    const VECTOR<double>& volsPuts,
                        const VECTOR<double>& volsCalls,
                        const VECTOR<double>& interpolTypes,
						long whatIsInterpolated,
						long correctSplineWithLinear,
                        double fxSpot,
                        long domZcCurveId,
                        long forZcCurveId,
                        int inRRSTR,
                        int isATM,
					    ARM_result& result,
					    long objId)
{
	long volId;

	ARM_FXVolCurve*          vc = NULL;
	ARM_FXVolCurve* newFxVolCrv = NULL;

	ARM_Vector* vMatu           = NULL;
    ARM_Vector* vFxFwds         = NULL;

    ARM_Vector* vPivotVols      = NULL;
    ARM_Vector* vPivotTypes     = NULL;

    ARM_Vector* vDeltaCalls     = NULL;
    ARM_Vector* vDeltaPuts      = NULL;
 
    ARM_Vector* vInterpolTypes  = NULL;
    ARM_ZeroCurve* domZcCurve   = NULL;
    ARM_ZeroCurve* forZcCurve   = NULL;

	int matu_size   = matu.size ();
	int fxFwds_size = fxFwds.size ();

    int pivotVols_size  = pivotVols.size();
    int pivotTypes_size = pivotTypes.size();

	int testSizeOK = ( pivotVols_size == pivotTypes_size )
                     && ( matu_size == pivotTypes_size );

	if (!(testSizeOK))
	{
	   result.setMsg ("ARM_ERR: Check your Vectors sizes: must be the same!");
		
       return(ARM_KO);
	}

	
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb in accessing objects");
		
       return(ARM_KO);
	}

	
	char AsOfDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date, AsOfDate);

        if ( domZcCurveId != ARM_NULL_OBJECT )
        {
           domZcCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(domZcCurveId);

		    if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(domZcCurve, ARM_ZERO_CURVE) == 0 )
		    {
		       result.setMsg ("ARM_ERR: In ZC Curve is not of a good type");
			    
               return(ARM_KO);
		    }
        }

        if ( forZcCurveId != ARM_NULL_OBJECT )
        {
            forZcCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(forZcCurveId);

		    if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(forZcCurve, ARM_ZERO_CURVE) == 0 )
		    {
		       result.setMsg ("ARM_ERR: In ZC Curve is not of a good type");
			    
               return(ARM_KO);
		    }
        }

        vMatu = CreateARMVectorFromVECTOR(matu);
	
        // Process vol Matrixes

            // PUT VOLS
        int volsPut_size = volsPuts.size();

        double dVolsPut[500];
        
        int i;

		for (i = 0; i < volsPut_size; i++)
		{
			dVolsPut[i] = volsPuts[i];
		}
		
        int deltaPuts_size  = deltaPuts.size();

	    ARM_Matrix mVolPuts(matu_size, deltaPuts_size, dVolsPut);


                 // CALL VOLS
        int volsCall_size = volsCalls.size();

        double dVolsCall[500];

		for (i = 0; i < volsCall_size; i++)
		{
			dVolsCall[i] = volsCalls[i];
		}
		
        int deltaCalls_size  = deltaCalls.size();

	    ARM_Matrix mVolCalls(matu_size, deltaCalls_size, dVolsCall);

        // End Process Vols
       

		vFxFwds     = CreateARMVectorFromVECTOR(fxFwds);
        vPivotVols  = CreateARMVectorFromVECTOR(pivotVols);
        vPivotTypes = CreateARMVectorFromVECTOR(pivotTypes);
        vDeltaCalls = CreateARMVectorFromVECTOR(deltaCalls);
        vDeltaPuts  = CreateARMVectorFromVECTOR(deltaPuts);

        vInterpolTypes = CreateARMVectorFromVECTOR(interpolTypes);

		newFxVolCrv = new ARM_FXVolCurve((ARM_Date) AsOfDate,
                                         *vMatu, 
					                     *vPivotVols, 
					                     *vPivotTypes,
					                     *vDeltaCalls,
					                     mVolCalls,
					                     *vDeltaPuts,
					                     mVolPuts,
					                     ARM_Vector(vFxFwds),
					                     *vInterpolTypes,
										 whatIsInterpolated,
										 correctSplineWithLinear,
                                         fxSpot,
                                         domZcCurve,
                                         forZcCurve,
                                         inRRSTR,
                                         isATM);

		if (vMatu)
			delete vMatu;
		vMatu = NULL;

		if (vFxFwds)
			delete vFxFwds;
		vFxFwds = NULL;

		if (vPivotVols)
			delete vPivotVols;
		vPivotVols = NULL;

		if (vPivotTypes)
			delete vPivotTypes;
		vPivotTypes = NULL;

		if (vDeltaCalls)
			delete vDeltaCalls;
		vDeltaCalls = NULL;
		
		if (vDeltaPuts)
			delete vDeltaPuts;
		vDeltaPuts = NULL;

		if (vInterpolTypes)
			delete vInterpolTypes;
		vInterpolTypes = NULL;

		if ( newFxVolCrv == NULL )       
		{
		   result.setMsg ("ARM_ERR: FX VolCurve construction failed!");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
		   CREATE_GLOBAL_OBJECT();

		   volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newFxVolCrv);

		   if ( volId == RET_KO )
           {
			  if (newFxVolCrv)
				 delete newFxVolCrv;
			  newFxVolCrv = NULL;

			  result.setMsg("ARM_ERR: Pb with  inserting FX Vol object");				
				
              return(ARM_KO);
           }

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			vc = (ARM_FXVolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_FX_VOLAT) == 1 )
			{
				if (vc)
				{
				   delete vc;
					
                   vc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) newFxVolCrv, objId);

				return(ARM_OK);
			}
			else
			{
				if (newFxVolCrv)
				   delete newFxVolCrv;
				newFxVolCrv = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
				
                return(ARM_KO);
			}
		}
	}

	catch(Exception& x)
	{
		if (vMatu)
		   delete vMatu;
		vMatu = NULL;

		if (vFxFwds)
           delete vFxFwds;
        vFxFwds = NULL;

        if (vPivotVols)
           delete vPivotVols;
        vPivotVols = NULL;

        if (vPivotTypes)
           delete vPivotTypes;
        vPivotTypes = NULL;

        if (vDeltaCalls)
           delete vDeltaCalls;
        vDeltaCalls = NULL;
        
        if (vDeltaPuts)
           delete vDeltaPuts;
        vDeltaPuts = NULL;

        if (vInterpolTypes)
           delete vInterpolTypes;
        vInterpolTypes = NULL;

		if (newFxVolCrv)
		   delete newFxVolCrv;
		newFxVolCrv = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_NewComputeFxVolatility(long idCurve,
								     double moneyness, 
                                     double maturity,
								     ARM_result& result)
{
	double dVol;
	ARM_FXVolCurve* vc;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		vc = (ARM_FXVolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_FX_VOLAT) == 0 )
		{
		   result.setMsg("ARM_ERR: Fx VolCurve is not of a good type");
			
           return(ARM_KO);
		}

		dVol = vc->computeVol(moneyness, maturity);

		result.setDouble(dVol);

		return(ARM_OK);
	}

	catch(Exception& x)
	{
		ARM_RESULT();
    }
}



long ARMLOCAL_GetInfoFromFxVolatility(long idCurve, double*& res, int& nbColumns, 
                                      int& nbRows,
								      ARM_result& result)
{
	ARM_Matrix mRes;
	ARM_FXVolCurve* vc;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		vc = (ARM_FXVolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_FX_VOLAT) == 0 )
		{
		   result.setMsg("ARM_ERR: Fx VolCurve is not of a good type");
			
           return(ARM_KO);
		}

		mRes      = vc->get_info();
		nbColumns = mRes.GetNumCols();
		nbRows    = mRes.GetNumLines();
		
        res = new double[nbColumns*nbRows];

		for (int i = 0; i < nbRows; i++)
			for (int j = 0; j < nbColumns; j++)
				res[i*nbColumns+j] = mRes.Elt(i,j);
		

		result.setDouble(0);

		return(ARM_OK);
	}

	catch(Exception& x)
	{
		ARM_RESULT();
    }
}



long ARMLOCAL_volflat(double volFlat,
					  double date,
                      const CCString& ccy,
					  ARM_result& result,
					  long objId)
{
	long volId;

    ARM_Currency currency((const char *) ccy);


	ARM_VolFlat* createdVolCurveFlat = NULL;
	ARM_VolFlat* prevVolCurveFlat = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* myCurveDate=new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,myCurveDate);

		createdVolCurveFlat = new ARM_VolFlat((ARM_Date) myCurveDate, 
                                              volFlat,
                                              &currency);

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		if (createdVolCurveFlat == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCurveFlat);

			if (volId == RET_KO)
			{
				if (createdVolCurveFlat)
					delete createdVolCurveFlat;
				createdVolCurveFlat = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			prevVolCurveFlat = (ARM_VolFlat *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevVolCurveFlat, ARM_VOL_FLAT) == 1)
			{
				if (prevVolCurveFlat)
				{
					delete prevVolCurveFlat;
					prevVolCurveFlat = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCurveFlat, objId);

				return ARM_OK;
			}
			else
			{
				if (createdVolCurveFlat)
				   delete createdVolCurveFlat;
				createdVolCurveFlat = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdVolCurveFlat)
			delete createdVolCurveFlat;
		createdVolCurveFlat = NULL;
		
		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_ComputeVolatility(long idCurve,
								double matu,
								double strike,
								double tenor,
								ARM_result& result)
{
	double dVol;
	ARM_VolCurve* vc;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		vc = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}
	
		dVol = vc->ComputeVolatility(matu, strike, tenor);

		result.setDouble(dVol);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}



long ARMLOCAL_ComputeModelVolatility (long idModel,
									  double matu,
									  double tenor,
									  double fwd,
									  double strike,
									  ARM_result& result,
									  int useSabr)
{
	double dVol;
	ARM_Model* model;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		model = (ARM_Model*) LOCAL_PERSISTENT_OBJECTS->GetObject(idModel);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(model, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}
	
		dVol = model->ComputeVol(matu, tenor, fwd, strike, useSabr);

		result.setDouble(dVol);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}



long ARMLOCAL_VolCube(long ATMVolId,
					  const VECTOR<long>& volCurveIds,
					  const VECTOR<double>& tenors,
                      long volType,
					  long checkCcy,
					  ARM_result& result,
					  long objId)
{
	long volId;

	ARM_VolCube*  createdVolCurve = NULL;
	ARM_VolCube*  prevVolCube     = NULL;

	ARM_VolCurve* ATMVolCrv = NULL;
	ARM_VolCurve* volCv     = NULL;
	ARM_VolCurve* inVols[200];
	double dTenors[200];



	int nbVolCurveIds = volCurveIds.size();

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		ATMVolCrv = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ATMVolId);

		if( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ATMVolCrv, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: ATM VolCurve is not of a good type");
			return ARM_KO;
		}

		for (int i = 0; i < nbVolCurveIds; i++)
		{
			volCv = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(volCurveIds[i]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCv, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: VolCurve is not of a good type");
				return ARM_KO;
			}

			inVols[i] = volCv;

			dTenors[i] = tenors[i];
		}

		ARM_Vector tenorsVect(nbVolCurveIds, dTenors);

		createdVolCurve = new ARM_VolCube(ATMVolCrv,
										  (ARM_VolCurve** ) inVols,
										  nbVolCurveIds,
										  &tenorsVect,
                                          (int) volType,
										  (int) checkCcy);

		if ( createdVolCurve == NULL )
		{
		   result.setMsg("ARM_ERR: VolCurve is null");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCurve);

			if ( volId == RET_KO )
			{
				if (createdVolCurve)
					delete createdVolCurve;
				createdVolCurve = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				return(ARM_KO);
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			prevVolCube = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevVolCube, ARM_VOL_CUBE) == 1)
			{
				if (prevVolCube)
				{
					delete prevVolCube;
					prevVolCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCurve, objId);

				return ARM_OK;
			}
			else
			{
				if (createdVolCurve)
					delete createdVolCurve;
				createdVolCurve = NULL;

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



ARM_VolLInterpol* ARMLOCAL_GetVolATMFromSummit(const CCString& index,
											   const CCString& currency,
											   const CCString& cvName,
											   ARM_Date date,
											   const CCString& vtype,
				   							   ARM_result& result)
{
	ARM_VolLInterpol *newVolCrv = NULL;
	char* sCurrency = currency;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	/// A conversion to char* creates a new... should be released with a delete!
	delete sCurrency;

	FILE *Fp = NULL;

	ARM_result C_result;

	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;


	char sDate[11];
	date.JulianToStrDate(sDate);
	double dDate = Local_ARMDATE2XLDATE(sDate);

	if (myDateFromFile > dateToday)
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + VOL_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + VOL_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"IRFWDVOL_" + currency + "_" + index + "_" + vtype + "_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"IRFWDVOL_" + currency + "_" + index + "_" + vtype + "_" + cvName + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg ("ARM_ERR: File not found");
			return NULL;
		}
	}

	long spotDays = sCCY->GetSpotDays();
	ARM_Date settleDate = asOfDate;
    char* currencyChar = currency.c_str();
	settleDate.NextBusinessDay(spotDays,currencyChar);
    delete currencyChar;
	settleDate.JulianToStrDate(sDate);
	double dSettleDate = Local_ARMDATE2XLDATE(sDate);

	vector<CCString> sStrikes;
	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* vVols = NULL;

	int indiceI;
	int indiceJ;
	int i, j;
	int compteur(0);

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;

		for (i = 0; i < sStrikes.size(); i++)
		{
			if ( strcmp((const char*) sStrikes[i], sEch) == 0 )
			{
				indiceI = i;
				i = sStrikes.size();
			}
		}

		if ( indiceI == -1 )
		{
			indiceI = sStrikes.size();
			sStrikes.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		indiceJ = -1;

		for (j = 0; j < sYearTerms.size(); j++)
		{
			if ( strcmp((const char*)sYearTerms[j], sEch) == 0 )
			{
				indiceJ = j;
				j = sYearTerms.size();
			}
		}

		if ( indiceJ == -1 )
		{
			indiceJ = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	vStrikes = new ARM_Vector(sStrikes.size());
	vYearTerms = new ARM_Vector(sYearTerms.size());
	vVols = new ARM_Matrix(sYearTerms.size(),sStrikes.size());

	long Nb;
	char cMatu;
	long freqId;

	for (i = 0; i < sStrikes.size(); i++)
	{
		long isDate = 0;

		sscanf(sStrikes[i], "%ld%c", &Nb, &cMatu);

		cMatu = toupper(cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
			freqId = K_DAILY;
		else if ( cMatu == 'W' )  
			freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
			freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			freqId = K_ANNUAL;
		else
		{
			isDate = 1;
			char sTmpDate[11];
			char sTmpDate1[7];

			strncpy(sTmpDate1, sStrikes[i], 6);
			sTmpDate1[6] = '\0';
			
            sprintf(sTmpDate,"%s%s%s",sTmpDate1,"20",((const char*) sStrikes[i])+6);

			if ( strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0 )
			{
			   ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
				
               vStrikes->Elt(i) = (tmpDate-asOfDate)/365.0;
            }
			else
			{
			   ARM_Date tmpDate(sTmpDate);

			   vStrikes->Elt(i) = (tmpDate-asOfDate)/365.0;
			}
		}
	
		if ( isDate == 0 )
		{
            /* OLD CODE: MA
			if (freqId == K_DAILY)
				long retCode = ARMLOCAL_NextBusinessDay(dDate,(const char*)currency,
                                                        Nb,
                                                        C_result);
			else
				long retCode = ARMLOCAL_ARM_ADDPERIOD(dDate, freqId, (const char*)currency, (long) Nb, 0, C_result);

			double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			vStrikes->Elt(i) = (myDate - dDate) /365.;
            */

            vStrikes->Elt(i) = FromStrMatuToDouble((const char *) sStrikes[i], &date);
		}
	}

	for (i = 0; i < sYearTerms.size(); i++)
	{
		long isDate = 0;

		// Contrat
		if ( strlen((const char*) sYearTerms[i]) == 5 )
		{
			int month, year;
			ARM_Date matDate;

			GetMonthYearFromExpiryDate(sYearTerms[i], &month, &year);
			matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery();
			vYearTerms->Elt(i) = (matDate.GetJulian()-date.GetJulian())/365.0;
		}
		else
		{
			sscanf(sYearTerms[i], "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
			{
				isDate = 1;
				char sTmpDate[11];
				char sTmpDate1[11];
                // TMP
				strncpy(sTmpDate1,sYearTerms[i],6);
				sTmpDate1[6] = '\0';
				
                sprintf(sTmpDate,"%s%s%s",sTmpDate1,"20",((const char*) sYearTerms[i])+6);

				if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
				{
					ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
				else
				{
					ARM_Date tmpDate(sTmpDate);
					vYearTerms->Elt(i) = (tmpDate-asOfDate)/365.0;
				}
			}

			if ( isDate == 0 )
			{
				double myDate;

				if ( freqId == K_DAILY )
				{
					ARM_Date tmpDate(settleDate);
					tmpDate.AddDays(Nb);
					char* ccy = currency.c_str();
					tmpDate.AdjustToBusDate(ccy,K_MOD_FOLLOWING);
					delete ccy;
					tmpDate.JulianToStrDate(sDate);
					myDate = Local_ARMDATE2XLDATE(sDate);
					//long retCode = ARMLOCAL_NextBusinessDay(dSettleDate,(const char*)currency,Nb,C_result);
				}
				else
				{
					long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)currency, (long) Nb, K_MOD_FOLLOWING, 0, C_result);
					myDate = Local_ARMDATE2XLDATE(C_result.getString());
				}
				vYearTerms->Elt(i) = (myDate - dDate) /365.;
/*				long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)currency, (long) Nb, K_MOD_FOLLOWING, C_result);

				double myDate = Local_ARMDATE2XLDATE(C_result.getString());

				vYearTerms->Elt(i) = (myDate - dDate) /365.;
*/			}
		}	
	}

	for (i=0;i<vVols->GetNumLines();i++)
		for (j=0;j<vVols->GetNumCols();j++)
			vVols->Elt(i,j) = Vols[j*vVols->GetNumLines()+i];


	newVolCrv = new ARM_VolLInterpol( asOfDate, vYearTerms,
								vStrikes, vVols, 1, K_ATMF_VOL);

	newVolCrv->SetCurrencyUnit(sCCY);
	// newVolCrv->SetIndexName((char*)(const char*)index);
	newVolCrv->SetIndexName(CCSTringToSTLString(index));

	if (Vols)
		delete [] Vols;
	Vols = NULL;

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	return newVolCrv;
}


ARM_VolLInterpol* ARMLOCAL_FF_GetCorrelFromSummit(const CCString& ccy1,
												  const CCString& index1,
												  const CCString& ccy2,
												  const CCString& index2,
												  ARM_Date date,
												  const CCString& cvname)
{
	ARM_VolLInterpol *newVolCrv = NULL;
	
	FILE *Fp = NULL;

	int isTranspose = 0;

	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;

	char sDate[11];
	date.JulianToStrDate(sDate);
	double dDate = Local_ARMDATE2XLDATE(sDate);

	if (myDateFromFile > dateToday)
	{
		return NULL;
	}

	CCString tmpCcy1 = ccy1;

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + CORREL_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + CORREL_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"INDEXCOR_ZERO_" + ccy1 + "_" + index1 + "_ZERO_" + ccy2 + "_" + index2 + "_" + cvname + ".";
	FileName2 = myRepertory2 + (CCString)"INDEXCOR_ZERO_" + ccy1 + "_" + index1 + "_ZERO_" + ccy2 + "_" + index2 + "_" + cvname + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			FileName.Replace(ccy1 + "_" + index1 + "_ZERO_" + ccy2 + "_" + index2,ccy2 + "_" + index2 + "_ZERO_" + ccy1 + "_" + index1);
			FileName2.Replace(ccy1 + "_" + index1 + "_ZERO_" + ccy2 + "_" + index2,ccy2 + "_" + index2 + "_ZERO_" + ccy1 + "_" + index1);
			if ((Fp = fopen(FileName,"r")) == NULL)
			{
				if ((Fp = fopen(FileName2,"r")) == NULL)
				{
					return NULL;
				}
				tmpCcy1 = ccy2;
				isTranspose = 1;
			}
			tmpCcy1 = ccy2;
			isTranspose = 1;
		}
	}

	vector<CCString> sStrikes;
	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* vVols = NULL;

	int indiceI;
	int indiceJ;
	int i, j;
	int compteur(0);

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;

		for (i = 0; i < sStrikes.size(); i++)
		{
			if ( strcmp((const char*) sStrikes[i], sEch) == 0 )
			{
				indiceI = i;
				i = sStrikes.size();
			}
		}

		if ( indiceI == -1 )
		{
			indiceI = sStrikes.size();
			sStrikes.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		indiceJ = -1;

		for (j = 0; j < sYearTerms.size(); j++)
		{
			if ( strcmp((const char*)sYearTerms[j], sEch) == 0 )
			{
				indiceJ = j;
				j = sYearTerms.size();
			}
		}

		if ( indiceJ == -1 )
		{
			indiceJ = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	vStrikes = new ARM_Vector(sStrikes.size());
	vYearTerms = new ARM_Vector(sYearTerms.size());
	vVols = new ARM_Matrix(sYearTerms.size(),sStrikes.size());

	ARM_Currency cCcy1((const char*)tmpCcy1);
	long spotDays = cCcy1.GetSpotDays();
    char* payCalTmp = cCcy1.GetPayCalName(cCcy1.GetVanillaIndexType());
    ARM_Date settleDate(date);
	settleDate.NextBusinessDay(spotDays, payCalTmp);

    char payCal[30];
    strcpy(payCal, payCalTmp);
    delete [] payCalTmp;

	for (i = 0; i < sStrikes.size(); i++)
	{
		vStrikes->Elt(i) = convPlotInYearTerm((const char *) sStrikes[i], date, 
                                                settleDate,
                                                payCal);
	}

	for (i = 0; i < sYearTerms.size(); i++)
	{
		vYearTerms->Elt(i) = convPlotInYearTerm((const char *) sYearTerms[i], date, 
                                                settleDate,
                                                payCal);
	}
	
	for (i=0;i<vVols->GetNumLines();i++)
		for (j=0;j<vVols->GetNumCols();j++)
			vVols->Elt(i,j) = Vols[j*vVols->GetNumLines()+i];

	if (isTranspose == 1)
	{
		vVols->Transpose();

		newVolCrv = new ARM_VolLInterpol(date,vStrikes,vYearTerms,
										 vVols, 1, K_ATMF_VOL);
	}
	else
	{
		newVolCrv = new ARM_VolLInterpol(date,vYearTerms,
										 vStrikes, vVols, 1, K_ATMF_VOL);

	}

	for (i = 0; i < ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}
	
	if(ccy1==ccy2)
	{
		char* tmpccy = (char*)ccy1;
		ARM_Currency* myCcy = new ARM_Currency(tmpccy);
		delete [] tmpccy;

		if(myCcy != ARM_DEFAULT_CURRENCY)
			newVolCrv->SetCurrencyUnit(myCcy);

		if (myCcy)
			delete myCcy;
		myCcy = NULL;
	}

	// Correl case
	// we need to fill itsYearTerms
	for (i = 0; i < newVolCrv->GetExpiryTerms()->GetSize(); i++)
	{
		sprintf(newVolCrv->itsYearTermsX[i], "%s", (const char*) sYearTerms[i]);
	}

	for (i = 0; i < newVolCrv->GetStrikes()->GetSize(); i++)
	{
		sprintf(newVolCrv->itsYearTermsY[i], "%s", (const char*) sStrikes[i]);
	}

	string	vIndex((const char *) ccy2);
	string	vCurrency((const char *) ccy1);
	string	vCrvId((const char *) cvname);
	string	vType("IRIR CORR");
	string	vNbFact((const char *) index1);
	if(vNbFact.substr(0, 3) == "COR")
	{
		vNbFact = vNbFact.substr(3, 2);
		vType += " " + vNbFact;
	}

	newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	if (Vols)
		delete [] Vols;
	Vols = NULL;

	return newVolCrv;

}



ARM_VolLInterpol* ARMLOCAL_GetSmileFromSummit(const CCString& index,
											  const CCString& currency,
											  const CCString& cvName,
											  ARM_Date date,
											  const CCString& vtype,
											  const CCString& matuIndex,
											  int smileType,
											  ARM_result& result)
{
	ARM_VolLInterpol *newVolCrv = NULL;
	char* sCurrency = currency;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	
	FILE *Fp = NULL;

	char sDate[11];
	ARM_result C_result;

	date.JulianToStrDate(sDate);
	double dDate = Local_ARMDATE2XLDATE(sDate);

	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;


	long spotDays = sCCY->GetSpotDays();
	ARM_Date settleDate = asOfDate;
	settleDate.NextBusinessDay(spotDays,sCurrency);
	settleDate.JulianToStrDate(sDate);
	double dSettleDate = Local_ARMDATE2XLDATE(sDate);

	delete sCurrency;

	if ( myDateFromFile > dateToday )
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + SMILE_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + SMILE_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"SMILE_" + currency + "_" + index + "_" + vtype + "_" + matuIndex + "_C_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"SMILE_" + currency + "_" + index + "_" + vtype + "_" + matuIndex + "_C_" + cvName + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg ("ARM_ERR: File not found");
			return NULL;
		}
	}

	vector<CCString> sStrikes;
	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* vVols = NULL;

	int indiceI;
	int indiceJ;
	int i, j;
	int compteur(0);

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<sYearTerms.size();i++)
		{
			if (strcmp((const char*)sYearTerms[i],sEch) == 0)
			{
				indiceI = i;
				i = sYearTerms.size();
			}
		}
		if (indiceI == -1)
		{
			indiceI = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		indiceJ = -1;
		for (j=0;j<sStrikes.size();j++)
		{
			if (strcmp((const char*)sStrikes[j],sEch) == 0)
			{
				indiceJ = j;
				j = sStrikes.size();
			}
		}
		if (indiceJ == -1)
		{
			indiceJ = sStrikes.size();
			sStrikes.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	if ( (sStrikes.size() == 0) || (sYearTerms.size() == 0) )
	{
		if (Vols)
			delete [] Vols;
		Vols = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		return NULL;
	}

	vStrikes = new ARM_Vector(sStrikes.size());
	vYearTerms = new ARM_Vector(sYearTerms.size());
	vVols = new ARM_Matrix(sYearTerms.size(),sStrikes.size());

	long Nb;
	char cMatu;
	long freqId;

	for (i=0;i<sStrikes.size();i++)
	{
		long isDouble = 0;

		sscanf(sStrikes[i], "%ld%c", &Nb, &cMatu);

		cMatu = toupper(cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
			freqId = K_DAILY;
		else if ( cMatu == 'W' )  
			freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
			freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			freqId = K_ANNUAL;
		else
		{
			isDouble = 1;

			vStrikes->Elt(i) = atof(sStrikes[i]);
		}
	
		if (isDouble == 0)
		{
			if (freqId == K_DAILY)
				long retCode = ARMLOCAL_NextBusinessDay(dDate,(const char*)currency,Nb,C_result);
			else
				long retCode = ARMLOCAL_ARM_ADDPERIOD(dDate, freqId, (const char*)currency, (long) Nb, 0, 0, C_result);

			double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			vStrikes->Elt(i) = (myDate - dDate) /365.;
		}
	}

	for (i=0;i<sYearTerms.size();i++)
	{
		long isDate = 0;

		// Contrat
		if (strlen((const char*) sYearTerms[i]) == 5)
		{
			int month, year;
			ARM_Date matDate;

			GetMonthYearFromExpiryDate(sYearTerms[i], &month, &year);
			matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery();
			vYearTerms->Elt(i) = (matDate.GetJulian() - date.GetJulian()) /365.;
		}
		else
		{
			sscanf(sYearTerms[i], "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
			{
				isDate = 1;
				char sTmpDate[11];
				char sTmpDate1[7];
				strncpy(sTmpDate1,sYearTerms[i],6);
				sTmpDate1[6] = '\0';
				sprintf(sTmpDate,"%s%s%s",sTmpDate1,"20",((const char*) sYearTerms[i])+6);
				if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
				{
					ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
				else
				{
					ARM_Date tmpDate(sTmpDate);
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
			}
		
			if (isDate == 0)
			{
				double myDate;

				if (freqId == K_DAILY)
				{
					ARM_Date tmpDate(settleDate);
					tmpDate.AddDays(Nb);
					char* ccy = currency.c_str();
					tmpDate.AdjustToBusDate(ccy,K_MOD_FOLLOWING);
					delete ccy;
					tmpDate.JulianToStrDate(sDate);
					myDate = Local_ARMDATE2XLDATE(sDate);
					//long retCode = ARMLOCAL_NextBusinessDay(dSettleDate,(const char*)currency,Nb,C_result);
				}
				else
				{
					long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)currency, (long) Nb, K_MOD_FOLLOWING, 0, C_result);
					myDate = Local_ARMDATE2XLDATE(C_result.getString());
				}
				vYearTerms->Elt(i) = (myDate - dDate) /365.;
			}
		}
	}

	for (i=0;i<vVols->GetNumLines();i++)
		for (j=0;j<vVols->GetNumCols();j++)
			vVols->Elt(i,j) = Vols[i*vVols->GetNumCols()+j];


	if(smileType == 0)
	{
		newVolCrv = new ARM_VolLInterpol( asOfDate, vYearTerms,
									vStrikes, vVols, 1, K_SMILE_VOL);
	}
	else
	{
		newVolCrv = new ARM_VolSplineInterpol(asOfDate, vYearTerms, vStrikes, vVols, 1, K_SMILE_VOL);
	}

	newVolCrv->SetCurrencyUnit(sCCY);
	newVolCrv->SetIndexName(CCSTringToSTLString(index));

	if (Vols)
		delete [] Vols;
	Vols = NULL;

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	return newVolCrv;
}




ARM_VolLInterpol* ARMLOCAL_GetVolSmileFromSummit(const CCString& index,
												 const CCString& currency,
												 const CCString& cvName,
												 ARM_Date date,
												 const CCString& vtype,
												 const CCString& matuIndex,
												 ARM_result& result)
{

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* VolATM = NULL;
	ARM_VolLInterpol* Smile = NULL;



	if ( GetDataRetrieverVersion() >= ETKRETRIEVER )
	{
		VolATM = etoolkit_GetVolATMFromSummit(index,
											  currency,
											  cvName,
											  date,
											  vtype);

		if ( VolATM == NULL )
		   return NULL;

		try
		{
			Smile = etoolkit_GetSmileFromSummit(index,
												currency,
												cvName,
												date,
												vtype,
												matuIndex,
												0 // interpol = LINEAR
												);
		}
		catch(...)
		{
			delete VolATM;
			return NULL;
		}

	}
	else
	{
		VolATM = ARMLOCAL_GetVolATMFromSummit(index,
											  currency,
											  cvName,
											  date,
											  vtype,
											  result);

		if( (VolATM == NULL) && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER))
		{
			VolATM = etoolkit_GetVolATMFromSummit(index,
												  currency,
												  cvName,
												  date,
												  vtype);

			if ( VolATM == NULL )
			   return NULL;
		}

		Smile = ARMLOCAL_GetSmileFromSummit(index,
											currency,
											cvName,
											date,
											vtype,
											matuIndex,
											0, // interpol = LINEAR
											result);

		if( (Smile == NULL)  && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER) )
		{
			try
			{
				Smile = etoolkit_GetSmileFromSummit(index,
													currency,
													cvName,
													date,
													vtype,
													matuIndex,
													0 // interpol = LINEAR
													);
			}
			catch(...)
			{
				delete VolATM;
				return NULL;
			}
		}
	}

	char* sCurrency = currency;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	delete sCurrency;

	long Nb;
	char cMatu;
	long freqId;

	sscanf(matuIndex, "%ld%c", &Nb, &cMatu);

	cMatu = toupper(cMatu);

	if ( cMatu == 'D' ) // Ex : "1D"
		freqId = K_DAILY;
	else if ( cMatu == 'W' )  
		freqId = K_WEEKLY;
	else if ( cMatu == 'M' ) 
		freqId = K_MONTHLY;
	else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
		freqId = K_ANNUAL;

	long retCode;

	char sDate[11];
	date.JulianToStrDate(sDate);
	double myAsOfDate = Local_ARMDATE2XLDATE(sDate);

	double myDate;

   


	if (freqId == K_DAILY)
		retCode = ARMLOCAL_NextBusinessDay(myAsOfDate, (const char*) currency,
                     Nb,result);
	else
		retCode = ARMLOCAL_ARM_ADDPERIOD(myAsOfDate, freqId, (const char*)currency, (long) Nb, 0, 0, result);

	if ( retCode == ARM_OK )
		myDate = Local_ARMDATE2XLDATE(result.getString());
	else
		return NULL;

	//double dStrike = (myDate-myAsOfDate)/365.0;


    double dStrike = FromStrMatuToDouble((const char *) matuIndex, &date);

    int i = 0;
	
	while (( dStrike > VolATM->GetStrikes()->Elt(i) ) 
           && 
           ( i < VolATM->GetStrikes()->GetSize() )
          )
	{
		i++;
	}
	
	if  ( i == VolATM->GetStrikes()->GetSize() )
		return NULL;

	ARM_pCol colStrike = VolATM->GetVolatilities()->GetCol(i);

	double x1, x2, y1, y2, x;
	double valueSmile;

	ARM_Matrix* mVolATM = VolATM->GetVolatilities();
	ARM_Matrix* mSmile = Smile->GetVolatilities();

	ARM_Vector* vYearTermsATM = (ARM_Vector*) VolATM->GetExpiryTerms()->Clone();
	ARM_Vector* vStrikesSmile = (ARM_Vector*) Smile->GetStrikes()->Clone();
	ARM_Vector* vYearTermsSmile = Smile->GetExpiryTerms();
	
	ARM_Matrix* mVolSmile = new ARM_Matrix(vYearTermsATM->GetSize(),vStrikesSmile->GetSize());

	int nblgSm = Smile->GetExpiryTerms()->GetSize();
	int nbcolSm = Smile->GetStrikes()->GetSize();
	int nblgATM = VolATM->GetExpiryTerms()->GetSize();

	ARM_result C_result;

	for (int j = 0; j < nblgATM; j++)
	{
		x = vYearTermsATM->Elt(j);
		long k = 0;

		while(k < nblgSm && vYearTermsSmile->Elt(k) < x)
			k++;

		if (x >= vYearTermsSmile->Elt(nblgSm-1))
		{
			for (i=0;i<nbcolSm;i++)
			{
				valueSmile = mSmile->Elt(nblgSm-1,i);
				mVolSmile->Elt(j,i) = valueSmile + colStrike[j];
			}
		}
		else if (x <= vYearTermsSmile->Elt(0))
		{
			for (i=0;i<nbcolSm;i++)
			{
				valueSmile = mSmile->Elt(0,i);
				mVolSmile->Elt(j,i) = valueSmile + colStrike[j];
			}
		}
		else
		{
			for (i = 0; i < nbcolSm; i++)
			{
				x1 = vYearTermsSmile->Elt(k-1);
				x2 = vYearTermsSmile->Elt(k);

				y1 = mSmile->Elt(k-1,i);
				y2 = mSmile->Elt(k,i);

				valueSmile = linInterpol(x,x1,y1,x2,y2);
				mVolSmile->Elt(j,i) = valueSmile + colStrike[j];
			}
		}
	}

	newVolCrv = new ARM_VolLInterpol(date, vYearTermsATM,
								     vStrikesSmile, mVolSmile, 
                                     1, K_SMILE_VOL);

	newVolCrv->SetCurrencyUnit(sCCY);
	newVolCrv->SetIndexName((char*)(const char*)index);


	if (VolATM)
		delete VolATM;
	VolATM = NULL;

	if (Smile)
		delete Smile;
	Smile = NULL;

	if (sCCY)
		delete (sCCY);
	sCCY = NULL;

	return newVolCrv;
}


long ARMLOCAL_GetVolFromSummit(const CCString& index,
							   const CCString& currency,
							   const CCString& cvName,
							   double date,
							   const CCString& vtype,
							   const CCString& matuIndex,
							   const CCString& impOrHist,
							   long indexId,
							   ARM_result& result,
							   long objId)
{
	long volId;

    string volImpOrHist("");

    if (!( impOrHist == "IRFWDVOL" ))
    {
       volImpOrHist = "HISTO";  
    }

    string smileOrATM("");

    string volType((const char *) vtype);

    int optionType = K_IRG; // OR K_SWOPT

    int volEnumType = K_ATMF_VOL; // K_ATMF_VOL(ATM), K_SMILE_VOL(SMILE), K_FX_VOL_SP_INTERP

    if ( vtype == "IRG" )
    {
       volType = "CAP";
    }
    else
    {
		if( !( volImpOrHist.empty() ) )
			volType += string(" ");

       volType = volType + volImpOrHist;

       optionType = K_SWOPT;
    }

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	const char* sMatuIndex = (const char*) matuIndex;

	char sDate[11];
	Local_XLDATE2ARMDATE(date, sDate);
	ARM_Date myDate(sDate);

	CCString msg(" ");
	
	try
	{
		ARM_IRIndex* theIndex; 
		LocalPersistent::get().convert(indexId,theIndex); 
		
		if ( strcmp(sMatuIndex,"ATM") == 0 )
		{
			if ( GetDataRetrieverVersion () >= ETKRETRIEVER )
			   newVolCrv = etoolkit_GetVolATMFromSummit(index,currency,cvName,myDate,vtype,impOrHist);
			else
			{
				newVolCrv = ARMLOCAL_GetVolATMFromSummit(index,currency,cvName,myDate,vtype,result);

				if ( (!newVolCrv) && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER))
					newVolCrv = etoolkit_GetVolATMFromSummit(index,currency,cvName,myDate,vtype,impOrHist);
			}

		
            smileOrATM = ""; // ATM VOL

            volEnumType = K_ATMF_VOL;
        }
		else if (( sMatuIndex[0] == 'S' ) || ( sMatuIndex[0] == 's' ))
		{
			char smileMat[20];
			char onechar;
			char dummy[20];

			sscanf(sMatuIndex,"%[^:]%c%s",dummy,&onechar,smileMat);

			if (GetDataRetrieverVersion () >= ETKRETRIEVER)
				newVolCrv = etoolkit_GetSmileFromSummit(index,currency,cvName,myDate,vtype,smileMat,
														0, // interpol = LINEAR
														impOrHist);
			else
			{
				newVolCrv = ARMLOCAL_GetSmileFromSummit(index,currency,cvName,myDate,vtype,smileMat,
														0, // interpol = LINEAR
														result);
				
				if (!newVolCrv)
					newVolCrv = etoolkit_GetSmileFromSummit(index,currency,cvName,myDate,vtype,smileMat,
															0, // interpol = LINEAR
															impOrHist);
			}
		
            smileOrATM = "SMILE";
            volEnumType = K_SMILE_VOL;
        }
		else
		{
			newVolCrv = ARMLOCAL_GetVolSmileFromSummit(index,currency,cvName,myDate,vtype,sMatuIndex,result);
		
            smileOrATM = "SMILED"; // Notice the "D" for a smiled VOL

            volEnumType = K_SMILE_VOL;
        }

		if ( newVolCrv == NULL )
		{
		   return(ARM_KO);
		}

		if (theIndex) newVolCrv->SetIndex(*theIndex); 
        newVolCrv->SetOptionType(optionType);
        newVolCrv->SetVolType(volEnumType);

        // Update Mkt data characteristics
        string	vType = string("VOL ") + volType;
		if( !(smileOrATM.empty()) )
			vType += string(" ") + smileOrATM;

		string	vIndex((const char*) index);
		string	vCurrency((const char*) currency);
		string	vCrvId((const char*) cvName);
		
		if (( index == "ROLIB" ) || ( index == "ROEUR" ))
        {
			vIndex = "RO";
        }
		else if (( index == "NULIB" ) || ( index == "NUEUR" ))
        {
			vIndex = "NU";
        }
		else if( vIndex.substr(0, 2) == "CO" )
		{
			vType = "IRIR CORR";
			vIndex = "";
		}
 
        newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		if ( objId == -1 )
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

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
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
}


long ARMLOCAL_GetVolCubeFromSummit(const CCString& index,
								   const CCString& currency,
								   const CCString& cvName,
								   double date,
								   const CCString& vtype,
								   VECTOR<CCString>& tenors,
								   const CCString& smileOrNot,
								   long indexId,
								   int smileType,
								   ARM_result& result,
								   long objId)
{
	long volId;

    VECTOR<CCString> smileFallBacks;

	ARM_VolCube*  createdVolCube = NULL;
	ARM_VolCube*  prevVolCube = NULL;

	ARM_VolCurve* ATMVolCrv = NULL;
	ARM_VolCurve* volCv     = NULL;
	ARM_VolCurve* inVols[200];
	double dTenors[200];

	char* ccy = (char*) currency;
	ARM_Currency* sCCY = new ARM_Currency(ccy);
	delete ccy;

	int nbCrv;
	char matIndex[20];

	char sDate[11];
	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

    string volType((const char *) vtype);

    if ( vtype == "IRG" )
    {
       volType = "CAP";
    }

    // Temporary code for treating MO40 FALL BACKS

    if ( cvName == "MO40" )
    {
       smileFallBacks.push_back("MO");
    }

    CCString msg (" ");

	

	try
	{
			//	safely convert ID to ARM_IRIndex. throw if bad type,
			//	convert to 0 if -1 
			ARM_IRIndex* theIndex ; 
			LocalPersistent::get().convert(indexId,theIndex); 

		if ( (smileOrNot[0] == 'A') ||(smileOrNot[0] == 'a') )
		{
			if ( GetDataRetrieverVersion () >= ETKRETRIEVER )
			{
			   ATMVolCrv = etoolkit_GetVolATMFromSummit(index, currency, cvName, myDate, vtype);
			}
			else
			{
				ATMVolCrv = ARMLOCAL_GetVolATMFromSummit(index, currency, cvName, myDate, vtype, result);

				if( (ATMVolCrv == NULL) && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER))
				{
					ATMVolCrv = etoolkit_GetVolATMFromSummit(index, currency, cvName, myDate, vtype);
				}
			}

			if ( ATMVolCrv == NULL )
			{
				result.setMsg("ARM_ERR: Pb with ATM VolCurve");				
				
                return(ARM_KO);
			}
		}

		if (theIndex) ATMVolCrv->SetIndex(*theIndex); 

		nbCrv = 0;

		VECTOR<CCString> listTenors = tenors;
		if ( (listTenors.size() == 0) && ( (GetDataRetrieverVersion () > FFRETRIEVER) || GetFallBackDataRetrieverVersion() > FFRETRIEVER) )
		{
			listTenors = ARMLOCAL_ParseListTenors(index, currency, cvName, myDate, vtype);
		}

		for (int i = 0; i < listTenors.size();i++)
		{
			strcpy(matIndex, listTenors[i]);

			if (( smileOrNot[0] == 'S' ) || ( smileOrNot[0] == 's' )
				|| 
                ( smileOrNot[0] == 'A') ||( smileOrNot[0] == 'a'))
			{
				if ( GetDataRetrieverVersion () >= ETKRETRIEVER )
				{
					try
					{
						volCv = etoolkit_GetSmileFromSummit(index,currency,cvName,myDate,vtype,matIndex, smileType);
					}

					catch(...)
					{
						volCv = NULL;

                        // Try fall backs

                        int size = smileFallBacks.size();
                            
                        int found = 0;

                        for (int h = 0; (( h < size ) && (!(found))); h++)
                        {
                            try
                            {
                           
                                volCv = etoolkit_GetSmileFromSummit(index, currency,
                                                                    smileFallBacks[h],
                                                                    myDate,
                                                                    vtype,
                                                                    matIndex, smileType);

                                if ( volCv != NULL )
                                {
                                   found = 1;
                                }
                            }

                            catch(...)
                            {
                                volCv = NULL;
                            }
                        }
					}
				}
				else
				{
					volCv = ARMLOCAL_GetSmileFromSummit(index,currency,cvName,myDate,vtype,matIndex,smileType,result);

					if( (volCv == NULL) && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER))
					{
						volCv = etoolkit_GetSmileFromSummit(index,currency,cvName,myDate,vtype,matIndex, smileType);
					}
				}
			}
			else
			{
			   volCv = ARMLOCAL_GetVolSmileFromSummit(index,currency,cvName,myDate,vtype,matIndex,result);
			}

			if ( volCv != NULL )
			{
				if (theIndex) volCv->SetIndex(*theIndex); 
				dTenors[nbCrv] = FromStrMatuToDouble(matIndex);			
				inVols[nbCrv] = volCv;
				nbCrv++;
			}
		}

		if ( nbCrv == 0 )
		{
			result.setMsg ("ARM_ERR: No Smile");
			return ARM_KO;
		}

		ARM_Vector existingTenors(nbCrv,dTenors);

		if ( (smileOrNot[0] == 'A') ||(smileOrNot[0] == 'a') )
		{
			createdVolCube = new ARM_VolCube(ATMVolCrv,inVols,nbCrv,&existingTenors);
		}
		else
		{
			createdVolCube = new ARM_VolCube(inVols,nbCrv,&existingTenors);
		}

		if (createdVolCube == NULL)
		{
			return ARM_KO;
		}

		if (theIndex) createdVolCube->SetIndex(*theIndex); 
		createdVolCube->SetCurrencyUnit(sCCY);

       // Update Mkt data characteristics
        string	vType = string("VOL ") + volType + " SMILE";
		string	vIndex((const char*) index);
		string	vCurrency((const char*) currency);
		string	vCrvId((const char*) cvName);

        createdVolCube->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		for (int i = 0; i < nbCrv; i++)
		{
			if (inVols[i])
				delete inVols[i];
			inVols[i] = NULL;
		}

		if (ATMVolCrv)
			delete ATMVolCrv;
		ATMVolCrv = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCube);

			if ( volId == RET_KO )
			{
				if (createdVolCube)
					delete createdVolCube;
				createdVolCube = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			prevVolCube = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevVolCube, ARM_VOL_CUBE) == 1)
			{
				if (prevVolCube)
				{
					delete prevVolCube;
					prevVolCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCube, objId);

				return ARM_OK;
			}
			else
			{
				if (createdVolCube)
					delete createdVolCube;
				createdVolCube = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		for (int i = 0; i < nbCrv; i++)
		{
			if (inVols[i])
				delete inVols[i];
			inVols[i] = NULL;
		}

		if (ATMVolCrv)
			delete ATMVolCrv;
		ATMVolCrv = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (createdVolCube)
			delete createdVolCube;
		createdVolCube = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_ARM_BumpVolatility(long VolId,
								 double valueToBump,
								 long nthLine,
								 long nthCol,
								 long cumulId,
								 long absoluteId,
								 ARM_result& result,
								 long objId)
{
	long newVolId;

	ARM_VolCurve* inVolCrv = NULL;
	ARM_VolCurve* outVolCrv = NULL;
	ARM_VolCurve* newOutVolCrv = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inVolCrv = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inVolCrv,ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In VolCurve is not of a good type");
			return ARM_KO;
		}

		newOutVolCrv = (ARM_VolCurve *)inVolCrv->Clone();
		
		if (newOutVolCrv == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		newOutVolCrv->BumpVolatility(valueToBump,nthLine,nthCol,cumulId,absoluteId);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if (newVolId == RET_KO)
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outVolCrv = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(outVolCrv,ARM_VOL_CURVE) == 1)
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}		
	}

	catch(Exception& x)
	{
		if (newOutVolCrv)
			delete newOutVolCrv;
		newOutVolCrv = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
}


long ARMLOCAL_ARM_BumpSmile(long VolId,
							double valueToBump,
							double tenor,
							long nthLine,
							long nthCol,
							long cumulId,
							long absoluteId,
							ARM_result& result,
							long objId)
{
	long newVolId;

	ARM_VolCube* inVolCube = NULL;
	ARM_VolCube* outVolCube = NULL;
	ARM_VolCube* newOutVolCube = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inVolCube = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inVolCube, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In VolCube is not of a good type");
			return ARM_KO;
		}

		newOutVolCube = (ARM_VolCube *)inVolCube->Clone();
		
		if (newOutVolCube == NULL)
		{
			result.setMsg ("ARM_ERR: VolCube is null");
			return ARM_KO;
		}

		newOutVolCube->BumpSmile(valueToBump,tenor,nthLine,nthCol,cumulId,absoluteId);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCube);

			if (newVolId == RET_KO)
			{
				delete newOutVolCube;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outVolCube = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(outVolCube, ARM_VOL_CURVE) == 1)
			{
				delete outVolCube;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCube, objId);

				return ARM_OK;
			}
			else
			{
				delete newOutVolCube;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}		
	}

	catch(Exception& x)
	{
		delete newOutVolCube;

		x.DebugPrint();

		ARM_RESULT();
    }
}


long ARMLOCAL_ARM_BumpHyperCubeSmile(long VolId,
									double valueToBump,
									double cubeTenor,
									double smileTenor,
									long nthLine,
									long nthCol,
									long cumulId,
									long absoluteId,
									ARM_result& result,
									long objId)
{
	long newVolId;

	ARM_HyperCube* inHyperCube = NULL;
	ARM_HyperCube* outHyperCube = NULL;
	ARM_HyperCube* newOutHyperCube = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inHyperCube = (ARM_HyperCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inHyperCube, ARM_HYPER_CUBE) == 0)
		{
			result.setMsg ("ARM_ERR: In HyperCube is not of a good type");
			return ARM_KO;
		}

		newOutHyperCube = (ARM_HyperCube *)inHyperCube->Clone();
		
		if (newOutHyperCube == NULL)
		{
			result.setMsg ("ARM_ERR: HyperCube is null");
			return ARM_KO;
		}

		newOutHyperCube->BumpSmile(valueToBump,cubeTenor,smileTenor,nthLine,nthCol,cumulId,absoluteId);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutHyperCube);

			if (newVolId == RET_KO)
			{
				delete newOutHyperCube;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outHyperCube = (ARM_HyperCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inHyperCube, ARM_HYPER_CUBE) == 1)
			{
				delete outHyperCube;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutHyperCube, objId);

				return ARM_OK;
			}
			else
			{
				delete newOutHyperCube;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}		
	}

	catch(Exception& x)
	{
		delete newOutHyperCube;

		x.DebugPrint();

		ARM_RESULT();
    }
}


long ARMLOCAL_ARM_FXBumpRRorSTR(long VolId,
								double valueToBump,
								long nthLine,
								long nthCol,
								double spotFX,
								long isCumul,
								long isAbsolute,
								long isRR,
								ARM_result& result,
								long objId)
{
	long newVolId;

	ARM_FXVolCurve* inVolCrv = NULL;
	ARM_FXVolCurve* outVolCrv = NULL;
	ARM_FXVolCurve* newOutVolCrv = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inVolCrv = (ARM_FXVolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inVolCrv,ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In VolCurve is not of a good type");
			return ARM_KO;
		}

		newOutVolCrv = (ARM_FXVolCurve *)inVolCrv->Clone();
		
		if (newOutVolCrv == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		newOutVolCrv->FXBumpRRorSTR(valueToBump,nthLine,nthCol,spotFX,isRR,isCumul,isAbsolute);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if (newVolId == RET_KO)
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outVolCrv = (ARM_FXVolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(outVolCrv,ARM_VOL_CURVE) == 1)
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}		
	}

	catch(Exception& x)
	{
		if (newOutVolCrv)
			delete newOutVolCrv;
		newOutVolCrv = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
}


ARM_VolLInterpol* ARMLOCAL_GetFileFxSmileFromSummit(const CCString& ccy1,
													const CCString& ccy2,
													const CCString& cvName,
													ARM_Date date,
													ARM_result& result)
{

	ARM_VolLInterpol* newVolCrv = NULL;
	
	FILE *Fp = NULL;

	char* sCurrency = ccy1;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	delete sCurrency;

	char sDate[11];
	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;
//	dateToday.SysToday();// donne la date systeme du jour

	if (myDateFromFile > dateToday)
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + FXSMILE_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + FXSMILE_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"FXSMILE_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"FXSMILE_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg ("ARM_ERR: File not found");
			return NULL;
		}
	}

	ARM_Date settleDate = asOfDate;
	settleDate.JulianToStrDate(sDate);
	double dSettleDate = Local_ARMDATE2XLDATE(sDate);

	vector<double> sStrikes;
	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;

	int indiceI;
	int indiceJ;
	int i;
	int compteur(0);
	double tmpStk;

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<sYearTerms.size();i++)
		{
			if (strcmp((const char*)sYearTerms[i],sEch) == 0)
			{
				indiceI = i;
				i = sYearTerms.size();
			}
		}

		if (indiceI == -1)
		{
			indiceI = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		tmpStk = atof(sEch);

		indiceJ = -1;
		for (i=0;i<sStrikes.size();i++)
		{
			if (sStrikes[i] == tmpStk)
			{
				indiceJ = i;
				i = sStrikes.size();
			}
		}

		if (indiceJ == -1)
		{
			indiceJ = sStrikes.size();
			sStrikes.push_back(tmpStk);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	vYearTerms = new ARM_Vector(sYearTerms.size());
	vStrikes = new ARM_Vector(sStrikes.size());

	long Nb;
	char cMatu;
	long freqId;


	for (i=0;i<sYearTerms.size();i++)
	{
		long isDate = 0;

		// Contrat
		if (strlen((const char*) sYearTerms[i]) == 5)
		{
			int month, year;
			ARM_Date matDate;

			GetMonthYearFromExpiryDate(sYearTerms[i], &month, &year);
			matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery();
			vYearTerms->Elt(i) = (matDate.GetJulian() - settleDate.GetJulian()) /365.;
		}
		else
		{
			sscanf(sYearTerms[i], "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
			{
				isDate = 1;
				char sTmpDate[11];
				char sTmpDate1[11];

				strncpy(sTmpDate1, sYearTerms[i], 6);
				sTmpDate1[6] = '\0';
				sprintf(sTmpDate, "%s%s%s",sTmpDate1,"20",
                        ((const char*) sYearTerms[i])+6);

				if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
				{
					ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
				else
				{
					ARM_Date tmpDate(sTmpDate);
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
			}
			if (isDate == 0)
			{
				long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)ccy1, (long) Nb, 0, 0, result);

				ARM_Date myDate(result.getString());

				vYearTerms->Elt(i) = (myDate - date) /365.;
			}

		}	
	}

	ARM_Matrix* mVols = new ARM_Matrix(vYearTerms->GetSize(),vStrikes->GetSize());

	for (i=0;i<vYearTerms->GetSize();i++)
	{
		for (int j = 0; j < vStrikes->GetSize(); j++)
		{
			if (i == 0)
				vStrikes->Elt(j) = sStrikes[j];

			mVols->Elt(i,j) = Vols[i*vStrikes->GetSize()+j];
		}
	}

	newVolCrv = new ARM_VolLInterpol( asOfDate, vYearTerms,
								vStrikes, mVols, K_STK_TYPE_PRICE, K_FX_VOL_SP_INTERP);



	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	return newVolCrv;
}

/*
long ARMLOCALFF_GetInitialFXVolFromSummit (const CCString& ccy1,
										   const CCString& ccy2,
										   const CCString& cvName,
										   const ARM_Date& date,
										   VECTOR<CCString>* maturities,
										   ARM_Vector* vVolATM);
{

	ARM_VolLInterpol* newVolCrv = NULL;
	
	FILE *Fp = NULL;

	char* sCurrency = ccy1;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	delete sCurrency;

	char sDate[11];
	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;
//	dateToday.SysToday();// donne la date systeme du jour

	if (myDateFromFile > dateToday)
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + FXVOL_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + FXVOL_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"FXVOL_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"FXVOL_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg ("ARM_ERR: File not found");
			return NULL;
		}
	}

	ARM_Date settleDate = asOfDate;
	settleDate.JulianToStrDate(sDate);
	double dSettleDate = Local_ARMDATE2XLDATE(sDate);

	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVols = NULL;

	int indiceI;
	int i;
	int compteur(0);

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<sYearTerms.size();i++)
		{
			if (strcmp((const char*)sYearTerms[i],sEch) == 0)
			{
				indiceI = i;
				i = sYearTerms.size();
			}
		}

		if (indiceI == -1)
		{
			indiceI = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	vYearTerms = new ARM_Vector(sYearTerms.size());

	long Nb;
	char cMatu;
	long freqId;


	for (i=0;i<sYearTerms.size();i++)
	{
		long isDate = 0;

		// Contrat
		if (strlen((const char*) sYearTerms[i]) == 5)
		{
			int month, year;
			ARM_Date matDate;

			GetMonthYearFromExpiryDate(sYearTerms[i], &month, &year);
			matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery();
			vYearTerms->Elt(i) = (matDate.GetJulian() - settleDate.GetJulian()) /365.;
		}
		else
		{
			sscanf(sYearTerms[i], "%d%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
			{
				isDate = 1;
				char sTmpDate[11];
				char sTmpDate1[6];
				strncpy(sTmpDate1,sYearTerms[i],6);
				sTmpDate1[6] = '\0';
				sprintf(sTmpDate,"%s%s%s",sTmpDate1,"20",((const char*) sYearTerms[i])+6);

				if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
				{
					ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
				else
				{
					ARM_Date tmpDate(sTmpDate);
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
			}
			if (isDate == 0)
			{
				long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)ccy1, (long) Nb, 0, result);

				ARM_Date myDate(result.getString());

				vYearTerms->Elt(i) = (myDate - date) /365.;
			}

		}	
	}

	ARM_Vector vVols(vYearTerms->GetSize());

	vStrikes = new ARM_Vector(1,100.);

	for (i=0;i<vVols.GetSize();i++)
		vVols.Elt(i) = Vols[i];

	mVols = new ARM_Matrix(vVols); 

	newVolCrv = new ARM_VolLInterpol( asOfDate, vYearTerms,
								vStrikes, mVols, K_STK_TYPE_PRICE, K_ATMF_VOL);

	if (Vols)
		delete [] Vols;
	Vols = NULL;

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	return newVolCrv;
}

*/

ARM_VolLInterpol* ARMLOCAL_GetFXVolATMFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist,
												 ARM_result& result)
{

	ARM_VolLInterpol* newVolCrv = NULL;
	
	FILE *Fp = NULL;

	char* sCurrency = ccy1;
	ARM_Currency* sCCY = new ARM_Currency(sCurrency);
	delete sCurrency;

	char sDate[11];
	ARM_Date asOfDate(date);
	ARM_Date myDateFromFile(date);

	ARM_Date dateToday;


	if (myDateFromFile > dateToday)
	{
		result.setMsg ("ARM_ERR: Invalid Date");
		return NULL;
	}

	CCString myRepertory;
	CCString myRepertory2;

	CCString FileName;
	CCString FileName2;

	myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + FXVOL_SUMMIT_FILE_LOCATION);
	myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + FXVOL_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

	FileName = myRepertory + (CCString)"FXVOL_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";
	FileName2 = myRepertory2 + (CCString)"FXVOL_" + ccy1 + "_" + ccy2 +  "_" + cvName + ".";

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

	if ((Fp = fopen(FileName,"r")) == NULL)
	{
		if ((Fp = fopen(FileName2,"r")) == NULL)
		{
			result.setMsg ("ARM_ERR: File not found");
			return NULL;
		}
	}

	ARM_Date settleDate = asOfDate;
	settleDate.JulianToStrDate(sDate);
	double dSettleDate = Local_ARMDATE2XLDATE(sDate);

	vector<CCString> sYearTerms;
	double* Vols = new double[1000];

	ARM_Vector* vStrikes = NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVols = NULL;

	int indiceI;
	int i;
	int compteur(0);

	rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<sYearTerms.size();i++)
		{
			if (strcmp((const char*)sYearTerms[i],sEch) == 0)
			{
				indiceI = i;
				i = sYearTerms.size();
			}
		}

		if (indiceI == -1)
		{
			indiceI = sYearTerms.size();
			sYearTerms.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		Vols[compteur] = val;

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	vYearTerms = new ARM_Vector(sYearTerms.size());

	long Nb;
	char cMatu;
	long freqId;


	for (i = 0; i < sYearTerms.size(); i++)
	{
		long isDate = 0;

		// Contrat
		if ( strlen((const char*) sYearTerms[i]) == 5 )
		{
			int month, year;
			ARM_Date matDate;

			GetMonthYearFromExpiryDate(sYearTerms[i], &month, &year);
			matDate.ChangeDate(1, month, year);

			matDate.PiborDelivery();
			vYearTerms->Elt(i) = (matDate.GetJulian() - settleDate.GetJulian()) /365.;
		}
		else
		{
			sscanf(sYearTerms[i], "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
			{
				isDate = 1;
				char sTmpDate[11];
				char sTmpDate1[7];
				strncpy(sTmpDate1,sYearTerms[i],6);
				sTmpDate1[6] = '\0';
				sprintf(sTmpDate,"%s%s%s",sTmpDate1,"20",((const char*) sYearTerms[i])+6);

				if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
				{
					ARM_Date tmpDate(sTmpDate,"MM/DD/YYYY");
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
				else
				{
					ARM_Date tmpDate(sTmpDate);
					vYearTerms->Elt(i) = (tmpDate - asOfDate) /365.;
				}
			}
			if (isDate == 0)
			{
				long retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, (const char*)ccy1, (long) Nb, 0, 0, result);

				ARM_Date myDate(result.getString());

				vYearTerms->Elt(i) = (myDate - date) /365.;
			}

		}	
	}

	ARM_Vector vVols(vYearTerms->GetSize());

	vStrikes = new ARM_Vector(1,100.);

	for (i=0;i<vVols.GetSize();i++)
		vVols.Elt(i) = Vols[i];

	mVols = new ARM_Matrix(vVols); 

	newVolCrv = new ARM_VolLInterpol( asOfDate, vYearTerms,
								vStrikes, mVols, K_STK_TYPE_PRICE, K_ATMF_VOL);

	if (Vols)
		delete [] Vols;
	Vols = NULL;

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	return newVolCrv;
}



ARM_VolLInterpol* ARMLOCAL_CreateFXVolFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist,
												 const CCString& VolType,
												 ARM_result& result)
{
	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* FXVolATM  = NULL;
	ARM_VolCube* newVolCub      = NULL;

	/// voltype is only used when ETK is used
	/// hence included here
	int volType;
	ARM_Vector* vStrikes    =  NULL;
	ARM_Vector* vYearTerms  = NULL;
	ARM_Matrix* mVolATM     = NULL;
	ARM_Matrix* mVols       = NULL;
	ARM_VolLInterpol* smile = NULL;

    string is3F("");

	if ( GetDataRetrieverVersion() >= ETKRETRIEVER )
	{
	   FXVolATM = etoolkit_GetFXVolATMFromSummit(ccy1, ccy2, date, cvName, impOrHist);
	}
	else
	{
		FXVolATM = ARMLOCAL_GetFXVolATMFromSummit(ccy1, ccy2, date, cvName, impOrHist, result);

		if ( (FXVolATM == NULL) && (GetFallBackDataRetrieverVersion() != 0))
			FXVolATM = etoolkit_GetFXVolATMFromSummit(ccy1, ccy2, date, cvName, impOrHist);

		if ( FXVolATM == NULL )
			return NULL;
	}


	if ( strcmp(VolType,"ATM") == 0 )	
	{
	   smile   = NULL;
		
       volType = K_ATMF_VOL;
	}
	else if (( strcmp(VolType,"SMILE") == 0 )
			 || ( strcmp(VolType,"FXSPI") == 0 )
			 || ( strcmp(VolType,"CSMILE") == 0 )
			 || ( strcmp(VolType,"CFXSPI") == 0 ) 
            )
	{
		if ( GetDataRetrieverVersion() >= ETKRETRIEVER )
		{
		   smile = etoolkit_GetXMLFxSmileFromSummit(ccy1, ccy2, cvName, date);
		}
		else
		{
			smile = ARMLOCAL_GetFileFxSmileFromSummit(ccy1, ccy2, cvName, date, result);

			if ( (smile == NULL) && (GetFallBackDataRetrieverVersion() != 0) )
			{
			   smile = etoolkit_GetXMLFxSmileFromSummit(ccy1, ccy2, cvName, date);
			}
		}

		volType = K_FX_VOL_SP_INTERP;

        is3F = "3F";
	}
	else 
	   throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			            "Invalid VolType Argument for getting FX Vol Smile");

	// on recupere les strikes de l'objet smile;
	if ( (smile) && ( volType == K_FX_VOL_SP_INTERP ) 
         && (( strcmp(VolType, "SMILE") == 0 ) || ( strcmp(VolType, "FXSPI") == 0 )) 
       )
	{
		vStrikes= new ARM_Vector(smile->GetStrikes()->GetSize()+1);
		vStrikes->Elt(0) = 0.0;
		for (int i = 1; i < vStrikes->GetSize(); i++)
		{
			vStrikes->Elt(i)=((ARM_Vector*) smile->GetStrikes())->Elt(i-1);
		}

		mVols = MergeMatrix(FXVolATM->GetVolatilities(), smile->GetVolatilities());
		
	}
	else if ( !(smile) && ( volType == K_FX_VOL_SP_INTERP ))
	{
		if (FXVolATM)
		   delete FXVolATM;
		FXVolATM = NULL;

		return(NULL);
	}
	else
	{
		vStrikes = new ARM_Vector(1, 0.0);

		mVols = (ARM_Matrix *) FXVolATM->GetVolatilities()->Clone();
	}

	if (( strcmp(VolType,"CSMILE") == 0 ) || ( strcmp(VolType,"CFXSPI") == 0 ))
	{
		ARM_Vector* noStrikes = new ARM_Vector(1, 0.0);
		
		ARM_VolCurve* inVols[1];
		inVols[0] = smile;
		
		ARM_Vector* underlying = new ARM_Vector(1, 0.0);

		newVolCub = new ARM_VolCube((ARM_VolCurve *) FXVolATM, 
                                    inVols, 1, underlying, volType);

		delete FXVolATM;
		delete underlying;

		newVolCrv = newVolCub;
	}
	else
       newVolCrv = new ARM_VolLInterpol(date, FXVolATM->GetExpiryTerms(),
									    vStrikes, mVols, K_STK_TYPE_PRICE, volType);

	if (smile)
	   delete smile;
	smile = NULL;

   // Update Mkt data characteristics
    string	vType;
	string	vIndex((const char*) ccy2);
	string	vCurrency((const char*) ccy1);
	string	vCrvId((const char*) cvName);

    if ( impOrHist == "HISTVOL" )
    {
		vType = string("FXVOL HISTO");
    }
    else if (( impOrHist == "FXVOL" ) && ( VolType == "ATM" ))
    {
		vType = string("FXVOL");
    }
    else if (( impOrHist == "FXVOL" ) && ( VolType == "SMILE" ))
    {
		vType = string("FXVOL SMILE");
    }
    else if (( impOrHist == "FXVOL" ) && ( VolType == "FXSPI" ))
    {
		vType = string("FXVOL SMILE 3F");
    }

    newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	return(newVolCrv);

}


long ARMLOCAL_GetFXVolFromSummit(const CCString& ccy1,
								 const CCString& ccy2,
								 double date,
								 const CCString& cvName,
								 const CCString& impOrHist,
								 const CCString& volType,
								 ARM_result& result,
								 long objId)
{
	ARM_result C_result;

	long volId;

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	CCString msg("");

	try
	{
		newVolCrv = ARMLOCAL_CreateFXVolFromSummit(ccy1,ccy2,myDate,cvName,impOrHist,volType,result);

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
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 1)
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
}


long ARMLOCAL_GetNewFXVolFromSummit(const CCString& ccy1,
									const CCString& ccy2,
									double date,
									const CCString& cvName,
									long domZcId,
									long forZcId,
									double fxSpot,
									const VECTOR<double>& forwards,
									long WhatIsInterpolated,
									long correctSplineWithLinear,
									long isATM,
									CCString& curClass,
									ARM_result& result,
									long objId)
{
	ARM_result C_result;

	long volId;

	ARM_VolCurve* newVolCrv = NULL;
	ARM_VolCurve* vc = NULL;

	ARM_ZeroCurve* domZc = NULL;
	ARM_ZeroCurve* forZc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	CCString msg("");

	try
	{
		if (GetDataRetrieverVersion() == FFRETRIEVER)
		{
			result.setMsg ("ARM_ERR: Function not Implemented without ETK");
			return ARM_KO;
		}

		if (domZcId != ARM_NULL_OBJECT)
		{
			domZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(domZcId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(domZc, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Domestic Zc Curve is not of a good type");
				return ARM_KO;
			}

			forZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(forZcId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(forZc, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Foreign Zc Curve is not of a good type");
				return ARM_KO;
			}

			newVolCrv = ARMLOCAL_CreateNewFXVolFromSummitWithCurves(ccy1,ccy2,myDate,cvName,domZc,forZc,fxSpot,WhatIsInterpolated,correctSplineWithLinear,isATM);

			if ( newVolCrv->GetName() == ARM_FX_VOLAT )
			   curClass = LOCAL_VOLAT_FX_CLASS;
			else
			   curClass = LOCAL_VOL_CURVE_LIN_CLASS;
		}
		else
		{
			newVolCrv = ARMLOCAL_CreateNewFXVolFromSummitWithForwards(ccy1,ccy2,myDate,cvName,forwards,WhatIsInterpolated,correctSplineWithLinear,isATM);
		}

		if ( newVolCrv == NULL )
		{
		   return(ARM_KO);
		}

        // Update Mkt data characteristics
		string	vType("FXVOL SMILE IMPROVED");
		string	vIndex((const char*) ccy2);
		string	vCurrency((const char*) ccy1);
		string	vCrvId((const char*) cvName);

        if(isATM)
			vType = string("FXVOL SMILE ATM");

        newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		if ( objId == -1 )
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
			vc = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_FX_VOLAT) == 1)
				|| (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
				)
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
}



long ARMLOCAL_GetInitialFXVolFromSummit(const CCString& ccy1,
										const CCString& ccy2,
										double date,
										const CCString& cvName,
										const CCString& impOrHist,
										const CCString& volType,
										VECTOR<CCString> *maturities,
										VECTOR<double> *tenors,
										VECTOR<double> *vol,
										ARM_result& result)
{
	ARM_Vector vVolATM(0);
	ARM_Matrix mSmile(0);
	ARM_Vector vTenors(0);

	char sDate[11];

	long retcode;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);
		ARM_Date myDate(sDate);

		if ( (GetDataRetrieverVersion() == FFRETRIEVER) && (GetFallBackDataRetrieverVersion() == FFRETRIEVER) )
		{
			result.setMsg ("ARM_ERR: Function not implemented without ETK");
			return ARM_KO;
		}

		retcode = etoolkit_GetInitialFXVolFromSummit (ccy1,
													  ccy2,
													  cvName,
													  myDate,
													  impOrHist,
													  maturities,
													  &vVolATM);

		if (retcode == ARM_KO)
		{
			result.setMsg("ARM_ERR : error in getting FX VOL ATM");
			return ARM_KO;
		}


		if ( (strcmp(volType,"SMILE")==0)
			|| (strcmp(volType,"FXSPI")==0)
			|| (strcmp(volType,"CSMILE")==0)
			|| (strcmp(volType,"CFXSPI")==0) )
		{
			retcode = etoolkit_GetInitialFXSmileFromSummit (ccy1,
															ccy2,
															cvName,
															myDate,
															impOrHist,
															&mSmile,
															&vTenors);
			if (retcode == ARM_KO)
			{
				result.setMsg("ARM_ERR : error in getting FX SMILE");
				return ARM_KO;
			}
		}

		// Renseignement des vecteurs de sortie
		tenors->push_back(0.0);

		for (int i = 0; i < vTenors.GetSize(); i++)
			tenors->push_back(vTenors.Elt(i));

		for ( int i = 0; i < maturities->size(); i++)
			vol->push_back(vVolATM.Elt(i));

		for (int j = 0; j < vTenors.GetSize(); j++)
		{
			for (int i = 0; i < maturities->size(); i++)
			{
				vol->push_back(mSmile.Elt(i,j));
			}
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_GetInitialVolFromSummit (const CCString& index,
									   const CCString& currency,
									   const CCString& cvName,
									   double date,
									   const CCString& vtype,
									   const CCString& matuIndex,
									   VECTOR<CCString> *maturities,
									   VECTOR<CCString> *tenors,
									   VECTOR<double> *vol,
									   ARM_result& result)
{
	if (GetDataRetrieverVersion () >= ETKRETRIEVER)
	{
		long retcode = etoolkit_GetInitialVolFromSummit (index,
														 currency,
														 cvName,
														 date,
														 vtype,
														 matuIndex,
														 maturities,
														 tenors,
														 vol);

		if (retcode == ARM_KO)
		{
			result.setMsg("ARM_ERR : error in getting Volatilities");
		}

		result.setLong(GetDataRetrieverVersion ());


		return retcode;
	}
	else
	{
		FILE *Fp = NULL;

		const char* sMatuIndex = (const char*) matuIndex;
		CCString myRepertory;
		CCString myRepertory2;
		CCString FileName;
		CCString FileName2;

		char sDate[11];
			
		if (strcmp(sMatuIndex,"ATM") == 0)
		{
			myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + VOL_SUMMIT_FILE_LOCATION);
			myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + VOL_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

			FileName = myRepertory + (CCString)"IRFWDVOL_" + currency + "_" + index + "_" + vtype + "_" + cvName + ".";
			FileName2 = myRepertory2 + (CCString)"IRFWDVOL_" + currency + "_" + index + "_" + vtype + "_" + cvName + ".";
		}
		else if ( (sMatuIndex[0] == 'S') || (sMatuIndex[0] == 's') )
		{
			char smileMat[20];
			char onechar;
			char dummy[20];

			sscanf(sMatuIndex,"%[^:]%c%s",dummy,&onechar,smileMat);

			myRepertory = ((CCString)(armlocal_init->data_folder.c_str()) + SMILE_SUMMIT_FILE_LOCATION);
			myRepertory2 = ((CCString)(armlocal_init->data_folder.c_str()) + SMILE_SUMMIT_FILE_LOCATION + SUMMIT_FILE_LOCATION_HISTO);

			FileName = myRepertory + (CCString)"SMILE_" + currency + "_" + index + "_" + vtype + "_" + smileMat + "_C_" + cvName + ".";
			FileName2 = myRepertory2 + (CCString)"SMILE_" + currency + "_" + index + "_" + vtype + "_" + smileMat + "_C_" + cvName + ".";
		}
		else
		{
			result.setMsg( "ARM_ERR: Check parameter matu Index" );
			return ARM_KO;
		}

		Local_XLDATE2ARMDATE(date,sDate);

		ARM_Date asOfDate(sDate);
		ARM_Date myDateFromFile(sDate);

		ARM_Date dateToday;
	//	dateToday.SysToday();// donne la date systeme du jour

		if (myDateFromFile > dateToday)
		{
			result.setMsg ("ARM_ERR: Invalid Date");
			return ARM_KO;
		}

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


		if ((Fp = fopen(FileName,"r")) == NULL)
		{
			if ((Fp = fopen(FileName2,"r")) == NULL)
			{
				if (GetFallBackDataRetrieverVersion() != 0)
				{
					long retcode = etoolkit_GetInitialVolFromSummit (index,
																	 currency,
																	 cvName,
																	 date,
																	 vtype,
																	 matuIndex,
																	 maturities,
																	 tenors,
																	 vol);

					result.setLong(GetFallBackDataRetrieverVersion());

					if (retcode == ARM_KO)
					{
						result.setMsg("ARM_ERR : error in getting Volatilities");
					}

					return retcode;
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
				result.setLong(GetDataRetrieverVersion ());
				return ARM_OK;
			}
		}
		else
		{
			fclose(Fp);

			result.setString(FileName);
			result.setLong(GetDataRetrieverVersion ());
			return ARM_OK;
		}
	}
}


long ARMLOCAL_GetFXCorrelFromSummit(const CCString& ccy1,
									const CCString& index,
									const CCString& ccy2,
									double date,
									const CCString& cvName,
									const VECTOR<CCString>& tenors,
									ARM_result& result,
									long objId)
{
/*	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}
*/
	long volId;

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	CCString msg("");

	try
	{
		newVolCrv = etoolkit_GetFXCorrelFromSummit(ccy1,index,ccy2,myDate,cvName,tenors);

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

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
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
}



long ARMLOCAL_GetCorrelFromSummit(const CCString& ccy1,
								  const CCString& index1,
								  const CCString& ccy2,
								  const CCString& index2,
								  double date,
								  const CCString& cvName,
								  ARM_result& result,
								  long objId)
{
	long volId;

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(date,sDate);
	ARM_Date myDate(sDate);

	CCString msg("");

	try
	{
		if ( GetDataRetrieverVersion() >= ETKRETRIEVER )
		{
			newVolCrv = etoolkit_GetCorrelFromSummit(ccy1,index1,ccy2,index2,myDate,cvName);
		}
		else
		{
			newVolCrv = ARMLOCAL_FF_GetCorrelFromSummit(ccy1,index1,ccy2,index2,myDate,cvName);

			if( (newVolCrv == NULL) && (GetFallBackDataRetrieverVersion() != 0))
			{
				newVolCrv = etoolkit_GetCorrelFromSummit(ccy1,index1,ccy2,index2,myDate,cvName);
			}
		}

		if (newVolCrv == NULL)
		{
			return ARM_KO;
		}

        newVolCrv->SetInterpType(K_DIAG_INTERPOL);

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

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
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
}



long ARMLOCAL_GetNthMaturity (long idCurve,
							  long nLine,
							  ARM_result& result)
{
	double dMat;
	ARM_VolCurve* vc;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		vc = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}

		dMat = vc->GetNthMaturity(nLine);

		result.setDouble(dMat);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}


long ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL(long volId,
										   long zcId,
										   long isSwoptVol,
										   long inPct,
										   ARM_result& result,
										   long objId)
{
	long newVolId;

	ARM_ZeroCurve* inZCCrv = NULL;
	ARM_VolLInterpol* inVolCrv = NULL;
	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inZCCrv = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inZCCrv,ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In ZC Curve is not of a good type");
			return ARM_KO;
		}

		inVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inVolCrv,ARM_VOL_LIN_INTERPOL) == 0)
		{
			result.setMsg ("ARM_ERR: In VolCurve is not of a good type");
			return ARM_KO;
		}

		newOutVolCrv = inVolCrv->ConvertToNormalVol(inZCCrv,
													isSwoptVol,
													inPct);

		if (newOutVolCrv == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if (newVolId == RET_KO)
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1)
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

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



long ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL(long volId,
										   long zcId,
										   long isSwoptVol,
										   long inPct,
										   long outPct,
										   ARM_result& result,
										   long objId)
{
	long newVolId;

	ARM_ZeroCurve* inZCCrv = NULL;
	ARM_VolLInterpol* inVolCrv = NULL;
	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		inZCCrv = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inZCCrv,ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: In ZC Curve is not of a good type");
			return ARM_KO;
		}

		inVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inVolCrv,ARM_VOL_LIN_INTERPOL) == 0)
		{
			result.setMsg ("ARM_ERR: In VolCurve is not of a good type");
			return ARM_KO;
		}

		newOutVolCrv = inVolCrv->ConvertToBSVol(inZCCrv,
												isSwoptVol,
												inPct,
												outPct);

		if (newOutVolCrv == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			newVolId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if (newVolId == RET_KO)
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(newVolId);

			return ARM_OK;
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1)
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return ARM_OK;
			}
			else
			{
				if (newOutVolCrv)
					delete newOutVolCrv;
				newOutVolCrv = NULL;

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



long ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(long tree3fId,
										 long PrcsId,
                                         const VECTOR<double>& dForwardVolDates,
										 ARM_result& result,
										 long objId)
{
	long volId;

	ARM_Tree3F* tree3fMod = NULL;
    ARM_PowerReverse* PRCSObject = NULL;

	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;


    char strDate[20];

    ARM_Vector* vForwardVolDates  = CreateARMVectorFromVECTOR(dForwardVolDates);

    if (vForwardVolDates!=NULL)
    {
        for (int i = 0; i < vForwardVolDates->GetSize(); i++)
        {
            Local_XLDATE2ARMDATE((*vForwardVolDates)[i], strDate);

            vForwardVolDates->Elt(i) = ((ARM_Date) strDate).GetJulian();
        }
    }
	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		tree3fMod = (ARM_Tree3F *) LOCAL_PERSISTENT_OBJECTS->GetObject(tree3fId);
        PRCSObject= (ARM_PowerReverse * ) LOCAL_PERSISTENT_OBJECTS->GetObject(PrcsId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(tree3fMod,ARM_TREE3F) == 0 )
		{
		   result.setMsg("ARM_ERR: tree 3f model is not of a good type");

		   return(ARM_KO);
		}

        if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PRCSObject,ARM_POWERREVERSE) == 0 )
		{
		   result.setMsg("ARM_ERR: PRCS Object is not of a good type");

		   return(ARM_KO);
		}


    	// MSG_printf_message (MSG_TRACE, "avant PRCS3F_ConvertObjSpotVolToFwdVol \n");

        try
        {
		     newOutVolCrv = PRCS3F_ConvertObjSpotVolToFwdVol(   tree3fMod->itsAsOfDate,
														        tree3fMod->itsDFBSModel->GetDomYieldCurve(),
														        tree3fMod->itsDFBSModel->GetDBsCrv(),
														        tree3fMod->itsDFBSModel->GetForeignYieldCurve(),
														        tree3fMod->itsDFBSModel->GetFBsCrv(),
														        tree3fMod->itsVolSwopBase,
														        tree3fMod->itsVolSwopForeign,
						                                        (ARM_VolLInterpol *) tree3fMod->itsDFBSModel->GetFxVol(),
														        tree3fMod->itsMeanReversionBase, 
														        tree3fMod->itsMeanReversionForeign,
														        tree3fMod->itsDFBSModel->GetdFxCorr(),// BaseSpotFXCorrel
														        tree3fMod->itsDFBSModel->GetfFxCorr(),// ForgnSpotFXCorrel
														        tree3fMod->itsBaseForeignCorrelation,
                                                                PRCSObject->GetItsNoticeDates(),
                                                                PRCSObject->GetItsFxUnderLeg()->GetResetDates(),
                                                                PRCSObject->GetItsFxNumLeg()->GetPaymentDates(),
                                                                tree3fMod->itsCutOff,
                                                                tree3fMod->itsLongDatedSpotFxVol,
														        tree3fMod->itsCalibSwoptBasis,
                                                                vForwardVolDates);

        }

        catch(Exception a3FExpt)
        {
            result.setMsg ("ARM_ERR: Pb in: PRCS3F_ConvertObjSpotVolToFwdVol");				
			
            return ARM_KO;
        }

       catch(...)
       {
           result.setMsg ("ARM_ERR: Failure in: PRCS3F_ConvertObjSpotVolToFwdVol");				
			
           return ARM_KO;
       }

             // MSG_printf_message (MSG_TRACE, "apres PRCS3F_ConvertObjSpotVolToFwdVol \n");

		if ( newOutVolCrv == NULL )
		{
		   result.setMsg("ARM_ERR: PRCS3F_ConvertObjSpotVolToFwdVol Failed!");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if ( volId == RET_KO )
			{
				if (newOutVolCrv)
					delete newOutVolCrv;

				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return(ARM_OK);
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1 )
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return(ARM_OK);
			}
			else
			{
				if (newOutVolCrv)
				   delete newOutVolCrv;
				newOutVolCrv = NULL;

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

 

long ARMLOCAL_InterpolInStrikeFwdTime (long idCurve,
									  double forward,
									  double strike,		
									  double matu,
									  double precision,
									  double sigmaATMF,
									  long y2NULL,
									  ARM_result& result)
{
	double dVol;
	ARM_VolLInterpol* vc;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		vc = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}

		ARM_FXVolSmileInterpol* fxvolsmile=vc->GetFxVolSmileInterpolation();
		if (fxvolsmile)
		{
			dVol = fxvolsmile->InterpolInStrikeFwdTime(forward, strike, matu, precision, sigmaATMF, y2NULL);
			
		}
		else
		{
			result.setMsg ("ARM_ERR: FxVolCurve is not of a good type");
			return ARM_KO;
		}

		result.setDouble(dVol);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}



long ARMLOCAL_ComputeFxVol(long idCurve,
						   double asof,
						   double matu,
						   double calcmatu,
						   double fxspot,
						   double strike,
						   long idDiscCrv,
						   long idDivCrv,
						   ARM_result& result)
{
	double dVol;
	ARM_VolLInterpol* vc = NULL;
	ARM_ZeroCurve* discount = NULL;
	ARM_ZeroCurve* dividend = NULL;



	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
        char myCurveDate[11];

		Local_XLDATE2ARMDATE(asof, myCurveDate);
        ARM_Date asOfDate(myCurveDate);

        Local_XLDATE2ARMDATE(matu, myCurveDate);
        ARM_Date matuDate(myCurveDate);


		vc = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 0 )
		{
		   result.setMsg("ARM_ERR: VolCurve is not of a good type");
			
           return(ARM_KO);
		}

		discount = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idDiscCrv);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(discount, ARM_ZERO_CURVE) == 0 )
		{
		   result.setMsg("ARM_ERR: Discount Zc is not of a good type");
			
           return(ARM_KO);
		}
	
        dividend = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idDivCrv);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dividend, ARM_ZERO_CURVE) == 0)
		{
		   result.setMsg("ARM_ERR: Dividend Zc is not of a good type");
			
           return(ARM_KO);
		}

		ARM_FXVolSmileInterpol* fxvolsmile = vc->GetFxVolSmileInterpolation();
		
        if (fxvolsmile)
		{
		   dVol = fxvolsmile->ComputeFxVol(asOfDate, matuDate, calcmatu, 
                                           fxspot, strike, 
                                           discount, dividend);			
		}
		else
		{
		   result.setMsg("ARM_ERR: FxVolCurve is not of a good type");
			
           return(ARM_KO);
		}

		result.setDouble(dVol);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}




long ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(   double asOf,
										    long dZcId,
                                            long fZcId,
                                            long dBSZcId,
				                            long fBSZcId,
                                            long volSwopBaseId,
						                    long volSwopForeignId,
                                            long fxVolId,
                                            double dMeanReversionBase,
						                    double dMeanReversionForeign,
                                            long dFxRdCorrId,
			                                long dFxRfCorrId,
                                            long dRdRfCorrId,
                                            double dCutOff,
                                            double dVolLongTerm,
                                            long calibBasisIncluded,
                                            const VECTOR<double>& dForwardVolDates,
										    ARM_result& result,
										    long objId)
{
	long volId;


	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;

    ARM_ZeroCurve* dZc = NULL;
    ARM_ZeroCurve* fZc = NULL;
    ARM_ZeroCurve* dBSZc = NULL;
    ARM_ZeroCurve* fBSZc = NULL;
    ARM_VolCurve* volSwopBase = NULL;
    ARM_VolCurve* volSwopForeign = NULL;	
    ARM_VolCurve* fxVol = NULL;
    ARM_VolCurve* dFxRdCorr = NULL;	
    ARM_VolCurve* dFxRfCorr = NULL;	
    ARM_VolCurve* dRdRfCorr = NULL;	
    


    char strDate[20];

    ARM_Vector* vForwardVolDates  = CreateARMVectorFromVECTOR(dForwardVolDates);

    for (int i = 0; i < vForwardVolDates->GetSize(); i++)
    {
        Local_XLDATE2ARMDATE((*vForwardVolDates)[i], strDate);

        vForwardVolDates->Elt(i) = ((ARM_Date) strDate).GetJulian();
    }


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char tmpDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(asOf, tmpDate);
		ARM_Date AsOfDate(tmpDate);

        ARM_Vector* NoticesDatesFictives = new ARM_Vector(30);
        for( int i = 1; i<=30;i++)
        {
            NoticesDatesFictives->Elt(i-1)=AsOfDate.GetJulian()+i*365;
        }

        dZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Domestic Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}
        
        fZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}

        dBSZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBSZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBSZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Domestic BS Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}

        fBSZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fBSZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fBSZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign BS Zero Curve is not a Zero Curve");
			return(ARM_KO);
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

        fxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fxVolId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVol, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Vol is not a Vol Curve");
			return(ARM_KO);
		}

        dFxRdCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dFxRdCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxRdCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Rd correlation is not a Vol Curve");
			return(ARM_KO);
        }

        dFxRfCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dFxRfCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxRfCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Rf correlation is not a Vol Curve");
			return(ARM_KO);
        }

		dRdRfCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dRdRfCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dRdRfCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Base Foreign correlation is not a Vol Curve");
			return(ARM_KO);
        }


        try
        {
		     newOutVolCrv = PRCS3F_ConvertObjSpotVolToFwdVol(   AsOfDate,
                                                                dZc,
                                                                dBSZc,
                                                                fZc,
                                                                fBSZc,
                                                                volSwopBase,
                                                                volSwopForeign,
                                                                (ARM_VolLInterpol *) fxVol,
                                                                dMeanReversionBase,
                                                                dMeanReversionForeign,
                                                                dFxRdCorr,
                                                                dFxRfCorr,
                                                                dRdRfCorr,
                                                                NoticesDatesFictives,
                                                                NoticesDatesFictives,
                                                                NoticesDatesFictives,
                                                                dCutOff,
                                                                dVolLongTerm,
                                                                calibBasisIncluded,
                                                                vForwardVolDates);

        }

        catch(Exception a3FExpt)
        {
            result.setMsg ("ARM_ERR: Pb in: PRCS3F_ConvertObjSpotVolToFwdVol");
            
            if(NoticesDatesFictives)
                delete NoticesDatesFictives;
			
            return ARM_KO;
        }

       catch(...)
       {
           result.setMsg ("ARM_ERR: Failure in: PRCS3F_ConvertObjSpotVolToFwdVol");
           
           if(NoticesDatesFictives)
                delete NoticesDatesFictives;
			
           return ARM_KO;
       }

             // MSG_printf_message (MSG_TRACE, "apres PRCS3F_ConvertObjSpotVolToFwdVol \n");

		if ( newOutVolCrv == NULL )
		{
		   result.setMsg("ARM_ERR: PRCS3F_ConvertObjSpotVolToFwdVol Failed!");

           if(NoticesDatesFictives)
                delete NoticesDatesFictives;

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

            if(NoticesDatesFictives)
                delete NoticesDatesFictives;

			if ( volId == RET_KO )
			{
				if (newOutVolCrv)
					delete newOutVolCrv;

				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return(ARM_OK);
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

            if(NoticesDatesFictives)
                delete NoticesDatesFictives;

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1 )
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return(ARM_OK);
			}
			else
			{
				if (newOutVolCrv)
				   delete newOutVolCrv;
				newOutVolCrv = NULL;

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



long ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(long tree3fId,
										 long PrcsId,
										 ARM_result& result,
										 long objId)
{
	long volId;

	ARM_Tree3F* tree3fMod = NULL;
    ARM_PowerReverse* PRCSObject = NULL;

	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;



    if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		tree3fMod = (ARM_Tree3F *) LOCAL_PERSISTENT_OBJECTS->GetObject(tree3fId);
        PRCSObject= (ARM_PowerReverse * ) LOCAL_PERSISTENT_OBJECTS->GetObject(PrcsId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(tree3fMod,ARM_TREE3F) == 0 )
		{
		   result.setMsg("ARM_ERR: tree 3f model is not of a good type");

		   return(ARM_KO);
		}

        if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PRCSObject,ARM_POWERREVERSE) == 0 )
		{
		   result.setMsg("ARM_ERR: PRCS Object is not of a good type");

		   return(ARM_KO);
		}


        try
        {
		     newOutVolCrv = PRCS3F_ConvertObjFwdVolToSpotVol(   tree3fMod->itsAsOfDate,
														        tree3fMod->itsDFBSModel->GetDomYieldCurve(),
														        tree3fMod->itsDFBSModel->GetDBsCrv(),
														        tree3fMod->itsDFBSModel->GetForeignYieldCurve(),
														        tree3fMod->itsDFBSModel->GetFBsCrv(),
														        tree3fMod->itsVolSwopBase,
														        tree3fMod->itsVolSwopForeign,
						                                        (ARM_VolLInterpol *) tree3fMod->itsDFBSModel->GetFxVol(),
														        tree3fMod->itsMeanReversionBase, 
														        tree3fMod->itsMeanReversionForeign,
														        tree3fMod->itsDFBSModel->GetdFxCorr(),// BaseSpotFXCorrel
														        tree3fMod->itsDFBSModel->GetfFxCorr(),// ForgnSpotFXCorrel
														        tree3fMod->itsBaseForeignCorrelation,
                                                                PRCSObject->GetItsNoticeDates(),
                                                                PRCSObject->GetItsFxUnderLeg()->GetResetDates(),
                                                                PRCSObject->GetItsFxNumLeg()->GetPaymentDates(),
                                                                tree3fMod->itsCutOff,
                                                                tree3fMod->itsLongDatedSpotFxVol,
														        tree3fMod->itsCalibSwoptBasis);

        }

        catch(Exception a3FExpt)
        {
            result.setMsg ("ARM_ERR: Pb in: PRCS3F_ConvertObjFwdVolToSpotVol");				
			
            return ARM_KO;
        }

        catch(...)
        {
            result.setMsg ("ARM_ERR: Failure in: PRCS3F_ConvertObjFwdVolToSpotVol");				
			
            return ARM_KO;
        }

	    if ( newOutVolCrv == NULL )
		{
		   result.setMsg("ARM_ERR: PRCS3F_ConvertObjFwdVolToSpotVol Failed!");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

			if ( volId == RET_KO )
			{
				if (newOutVolCrv)
					delete newOutVolCrv;

				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return(ARM_OK);
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1 )
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return(ARM_OK);
			}
			else
			{
				if (newOutVolCrv)
				   delete newOutVolCrv;
				newOutVolCrv = NULL;

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


long ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(   double asOf,
										    long dZcId,
                                            long fZcId,
                                            long dBSZcId,
				                            long fBSZcId,
                                            long volSwopBaseId,
						                    long volSwopForeignId,
                                            long fxVolId,
                                            double dMeanReversionBase,
						                    double dMeanReversionForeign,
                                            long dFxRdCorrId,
			                                long dFxRfCorrId,
                                            long dRdRfCorrId,
                                            double dCutOff,
                                            double dVolLongTerm,
                                            long calibBasisIncluded,
                                            ARM_result& result,
										    long objId)
{
	long volId;


	ARM_VolLInterpol* outVolCrv = NULL;
	ARM_VolLInterpol* newOutVolCrv = NULL;

    ARM_ZeroCurve* dZc = NULL;
    ARM_ZeroCurve* fZc = NULL;
    ARM_ZeroCurve* dBSZc = NULL;
    ARM_ZeroCurve* fBSZc = NULL;
    ARM_VolCurve* volSwopBase = NULL;
    ARM_VolCurve* volSwopForeign = NULL;	
    ARM_VolCurve* fxVol = NULL;
    ARM_VolCurve* dFxRdCorr = NULL;	
    ARM_VolCurve* dFxRfCorr = NULL;	
    ARM_VolCurve* dRdRfCorr = NULL;	
    
    

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char tmpDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(asOf, tmpDate);
		ARM_Date AsOfDate(tmpDate);

        ARM_Vector* NoticesDatesFictives = new ARM_Vector(30);
        for( int i = 1; i<=30;i++)
        {
            NoticesDatesFictives->Elt(i-1)=AsOfDate.GetJulian()+i*365;
        }

        dZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Domestic Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}
        
        fZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}

        dBSZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dBSZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBSZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Domestic BS Zero Curve is not a Zero Curve");
			return(ARM_KO);
		}

        fBSZc= (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fBSZcId);
        
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fBSZc, ARM_ZERO_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Foreign BS Zero Curve is not a Zero Curve");
			return(ARM_KO);
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

        fxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) fxVolId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVol, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Vol is not a Vol Curve");
			return(ARM_KO);
		}

        dFxRdCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dFxRdCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxRdCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Rd correlation is not a Vol Curve");
			return(ARM_KO);
        }

        dFxRfCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dFxRfCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxRfCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: FX Rf correlation is not a Vol Curve");
			return(ARM_KO);
        }

		dRdRfCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) dRdRfCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dRdRfCorr, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Base Foreign correlation is not a Vol Curve");
			return(ARM_KO);
        }


        try
        {
		     newOutVolCrv = PRCS3F_ConvertObjFwdVolToSpotVol(   AsOfDate,
                                                                dZc,
                                                                dBSZc,
                                                                fZc,
                                                                fBSZc,
                                                                volSwopBase,
                                                                volSwopForeign,
                                                                (ARM_VolLInterpol *) fxVol,
                                                                dMeanReversionBase,
                                                                dMeanReversionForeign,
                                                                dFxRdCorr,
                                                                dFxRfCorr,
                                                                dRdRfCorr,
                                                                NoticesDatesFictives,
                                                                NoticesDatesFictives,
                                                                NoticesDatesFictives,
                                                                dCutOff,
                                                                dVolLongTerm,
                                                                calibBasisIncluded);

        }

        catch(Exception a3FExpt)
        {
            result.setMsg ("ARM_ERR: Pb in: PRCS3F_ConvertObjFwdVolToSpotVol");
            
            if(NoticesDatesFictives)
                delete NoticesDatesFictives;
			
            return ARM_KO;
        }

       catch(...)
       {
           result.setMsg ("ARM_ERR: Failure in: PRCS3F_ConvertObjFwdVolToSpotVol");
           
           if(NoticesDatesFictives)
                delete NoticesDatesFictives;
			
           return ARM_KO;
       }

		if ( newOutVolCrv == NULL )
		{
		   result.setMsg("ARM_ERR: PRCS3F_ConvertObjFwdVolToSpotVol Failed!");

           if(NoticesDatesFictives)
                delete NoticesDatesFictives;

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv);

            if(NoticesDatesFictives)
                delete NoticesDatesFictives;

			if ( volId == RET_KO )
			{
				if (newOutVolCrv)
					delete newOutVolCrv;

				newOutVolCrv = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return(ARM_OK);
		}
		else
		{
			outVolCrv = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

            if(NoticesDatesFictives)
                delete NoticesDatesFictives;

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(outVolCrv,ARM_VOL_LIN_INTERPOL) == 1 )
			{
				if (outVolCrv)
				{
					delete outVolCrv;
					outVolCrv = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOutVolCrv, objId);

				return(ARM_OK);
			}
			else
			{
				if (newOutVolCrv)
				   delete newOutVolCrv;
				newOutVolCrv = NULL;

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


extern double ARMLOCAL_SetVolCurveName(long idCurve,
									   CCString c_name,
									   ARM_result& result)
{
	ARM_VolCurve* vc;
	c_name.toUpper();
	char* name = (char*)(const char*)(c_name);


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		vc = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vc, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}
	
		vc->SetIndexName(name);

		result.setDouble(0.);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}


// HyperCube
extern long ARMLOCAL_HyperCube(vector<long>& aVolCurveIds, 
							   vector<CCString>& aKeys, 
							   ARM_result& result, 
							   long objId = -1)
{
	long	hyperCubeId;
	ARM_HyperCube*  createdHyperCube = NULL;
	ARM_HyperCube*  prevHyperCube     = NULL;

	ARM_VolCurve*	volCurve = NULL;

	int	nbVolCurveIds = aVolCurveIds.size();

	vector<string>			keyList(nbVolCurveIds);
	vector<ARM_VolCurve*>	volCurveList(nbVolCurveIds);

	if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	try
	{
		for(int i=0; i < nbVolCurveIds; i++)
		{
			volCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(aVolCurveIds[i]);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: VolCurve is not of a good type");

				return ARM_KO;
			}

			volCurveList[i] = volCurve;
			keyList[i] = aKeys[i];
		}

		createdHyperCube = new ARM_HyperCube(volCurveList, keyList);

		if( createdHyperCube == NULL )
		{
		   result.setMsg("ARM_ERR: HyperCube is null");

		   return	ARM_KO;
		}

		if( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			hyperCubeId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdHyperCube);

			if( hyperCubeId == RET_KO )
			{
				if(createdHyperCube)
					delete	createdHyperCube;

				createdHyperCube = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				

				return	ARM_KO;
			}

			result.setLong(hyperCubeId);

			return	ARM_OK;
		}
		else
		{
			prevHyperCube = (ARM_HyperCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevHyperCube, ARM_HYPER_CUBE) == 1)
			{
				if(prevHyperCube)
				{
					delete	prevHyperCube;

					prevHyperCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdHyperCube, objId);

				return ARM_OK;
			}
			else
			{
				if(createdHyperCube)
					delete	createdHyperCube;

				createdHyperCube = NULL;

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


extern long ARMLOCAL_CreateCorrelCubeByExpiry(long C_hyperCubeId,
											  vector<CCString>& aTenorList,
											  vector<CCString>& aExpiryList,
											  CCString& aIntersurfaceInterpol,
											  ARM_result& result,
											  long objId)
{
	long	vCorrelCubeId;
	int		vNbTenors = aTenorList.size();
	int		vNbExpiries = aExpiryList.size();
	char*	vTenorStr;
	double	vTenor;
	bool	vIntersurfaceInterpol = false;
	if( aIntersurfaceInterpol == "YES" )
	{
		vIntersurfaceInterpol = true;
	}

	ARM_VolCube*	vPreviousCorrelCube = NULL;
	ARM_VolCube*	vCreatedCorrelCube = NULL;

	vector<double>	vTenor1List(vNbTenors);
	vector<double>	vTenor2List(vNbTenors);
	vector<double>	vExpiryList(vNbExpiries);

	for( int i=0; i<vNbTenors; i++)
	{
		vTenorStr = aTenorList[i].c_str();
		vTenor = StringMatuToYearTerm(vTenorStr);
		vTenor1List[i] = vTenor;
		vTenor2List[i] = vTenor;

		delete vTenorStr;
	}

	for(int j=0; j<vNbExpiries; j++)
	{
		vTenorStr = aExpiryList[j].c_str();
		vTenor = StringMatuToYearTerm(vTenorStr);
		vExpiryList[j] = vTenor;

		delete vTenorStr;
	}

	if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	try
	{			
		ARM_HyperCube*	vHyperCube = (ARM_HyperCube*) LOCAL_PERSISTENT_OBJECTS->GetObject(C_hyperCubeId);
		vHyperCube->SetIntersurfaceInterpol(vIntersurfaceInterpol);

		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vHyperCube, ARM_HYPER_CUBE) == 0)
		{
			result.setMsg ("ARM_ERR: HyperCube is not of a good type");

			return ARM_KO;
		}

		vCreatedCorrelCube = vHyperCube->CreateCorrelCubeByExpiry(vTenor1List, vTenor2List, vExpiryList);

		if( vCreatedCorrelCube == NULL )
		{
		   result.setMsg("ARM_ERR: CreatedCorrelCube is null");

		   return	ARM_KO;
		}

		if( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			vCorrelCubeId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vCreatedCorrelCube);

			if( vCorrelCubeId == RET_KO )
			{
				if(vCreatedCorrelCube)
					delete	vCreatedCorrelCube;

				vCreatedCorrelCube = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				

				return	ARM_KO;
			}

			result.setLong(vCorrelCubeId);

			return	ARM_OK;
		}
		else
		{
			vPreviousCorrelCube = (ARM_VolCube*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vPreviousCorrelCube, ARM_VOL_CUBE) == 1)
			{
				if(vPreviousCorrelCube)
				{
					delete	vPreviousCorrelCube;

					vPreviousCorrelCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vCreatedCorrelCube, objId);

				return ARM_OK;
			}
			else
			{
				if(vCreatedCorrelCube)
					delete	vCreatedCorrelCube;

				vCreatedCorrelCube = NULL;

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


extern long ARMLOCAL_ComputeCorrelByExpiry(long C_correlCubeId,
										   double C_Expiry,
										   double C_Tenor1,
										   double C_Tenor2,
										   ARM_result& result)
{
	ARM_VolCube*	vCorrelCube = NULL;

	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return	ARM_KO;
	}

	CCString	msg("");

	try
	{
		vCorrelCube = (ARM_VolCube*)LOCAL_PERSISTENT_OBJECTS->GetObject(C_correlCubeId);

		if( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vCorrelCube, ARM_VOL_CUBE) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vCorrelCube, ARM_HYPER_CUBE) == 0) )
		{
			result.setMsg ("ARM_ERR: CorrelCube is not of a good type");
			
			return	ARM_KO;
		}

		double	vCorrel = vCorrelCube->ComputeCorrelByExpiry(C_Expiry, C_Tenor1, C_Tenor2);

		result.setDouble(vCorrel);

		return	ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}


extern long ARMLOCAL_ComputeHyperCorrel(long C_correlCubeId,
										   double C_Tenor1,
										   double C_Tenor2,
										   double C_Expiry,
										   double C_Moneyness,
										   ARM_result& result)
{
	ARM_VolCube*	vCorrelCube = NULL;

	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return	ARM_KO;
	}

	CCString	msg("");

	try
	{
		vCorrelCube = (ARM_VolCube*)LOCAL_PERSISTENT_OBJECTS->GetObject(C_correlCubeId);

		if( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vCorrelCube, ARM_VOL_CUBE) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vCorrelCube, ARM_HYPER_CUBE) == 0) )
		{
			result.setMsg ("ARM_ERR: CorrelCube is not of a good type");
			
			return	ARM_KO;
		}

		double	vCorrel = vCorrelCube->ComputeHyperCorrel( C_Expiry, C_Tenor1, C_Tenor2, C_Moneyness );

		result.setDouble(vCorrel);

		return	ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}

// IndexIndexVolCube
extern long ARMLOCAL_IndexIndexCorrelCube(vector<long>& aVolCurveIds,
										  vector<CCString>& aTenors1List,
										  vector<CCString>& aTenors2List,
										  CCString& aIntersurfaceInterpol,
										  ARM_result& result,
										  long objId = -1)
{
	long	vIndexIndexCorrelCubeId;
	ARM_IndexIndexCorrelCube*  createdIndexIndexCorrelCube = NULL;
	ARM_IndexIndexCorrelCube*  prevIndexIndexCorrelCube    = NULL;

	ARM_VolCurve*	volCurve = NULL;

	int	nbVolCurveIds = aVolCurveIds.size();

	vector<string>			vTenors1List(nbVolCurveIds);
	vector<string>			vTenors2List(nbVolCurveIds);
	vector<ARM_VolCurve*>	volCurveList(nbVolCurveIds);

	if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg("");

	try
	{
		for(int i=0; i < nbVolCurveIds; i++)
		{
			volCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(aVolCurveIds[i]);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: VolCurve is not of a good type");

				return ARM_KO;
			}

			volCurveList[i] = volCurve;
			vTenors1List[i] = aTenors1List[i];
			vTenors2List[i] = aTenors2List[i];
		}

		createdIndexIndexCorrelCube = new ARM_IndexIndexCorrelCube(volCurveList, vTenors1List, vTenors2List);
		createdIndexIndexCorrelCube->SetIntersurfaceInterpol(aIntersurfaceInterpol == "YES" ? true: false);

		if( createdIndexIndexCorrelCube == NULL )
		{
		   result.setMsg("ARM_ERR: IndexIndexCorrelCube is null");

		   return	ARM_KO;
		}

		if( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			vIndexIndexCorrelCubeId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIndexIndexCorrelCube);

			if( vIndexIndexCorrelCubeId == RET_KO )
			{
				if(createdIndexIndexCorrelCube)
					delete	createdIndexIndexCorrelCube;

				createdIndexIndexCorrelCube = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				

				return	ARM_KO;
			}

			result.setLong(vIndexIndexCorrelCubeId);

			return	ARM_OK;
		}
		else
		{
			prevIndexIndexCorrelCube = (ARM_IndexIndexCorrelCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevIndexIndexCorrelCube, ARM_INDEX_INDEX_CORREL_CUBE) == 1)
			{
				if(prevIndexIndexCorrelCube)
				{
					delete	prevIndexIndexCorrelCube;

					prevIndexIndexCorrelCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIndexIndexCorrelCube, objId);

				return ARM_OK;
			}
			else
			{
				if(createdIndexIndexCorrelCube)
					delete	createdIndexIndexCorrelCube;

				createdIndexIndexCorrelCube = NULL;

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


extern long ARMLOCAL_ComputeIndexIndexCorrel(long C_correlCubeId,
											 double C_Tenor1,
											 double C_Tenor2,
											 double C_Expiry1,
											 double C_Expiry2,
											 ARM_result& result)
{
	ARM_IndexIndexCorrelCube*	vCorrelCube = NULL;

	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return	ARM_KO;
	}

	CCString	msg("");

	try
	{
		vCorrelCube = (ARM_IndexIndexCorrelCube*)LOCAL_PERSISTENT_OBJECTS->GetObject(C_correlCubeId);

		if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vCorrelCube, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: CorrelCube is not of a good type");
			
			return	ARM_KO;
		}
	
		double	vCorrel = vCorrelCube->ComputeIndexIndexCorrel(C_Tenor1, C_Expiry1, C_Tenor2, C_Expiry2);

		result.setDouble(vCorrel);

		return	ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }
}

extern long ARMLOCAL_BumpVolatilityCorrelManager(long correlManagerId,
										  string& C_ccy,
										  long TypeCorrel,
										  double value,
										  long nthLine,
										  long nthCol,
										  long isCumul,
										  long isAbsolute,
										  long isToClone,
										  ARM_result& result,
										  long ObjId)
{
	long newCorrelManagerId;

	ARM_CorrelatorManager* theCorrelManager = NULL;
	ARM_CorrelatorManager* newCorrelManager = NULL;
	ARM_CorrelatorManager* oldCorrelManager = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		theCorrelManager = (ARM_CorrelatorManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlManagerId);
		if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(theCorrelManager,ARM_CORRELMANAGER) == 0 )
		{
			result.setMsg ("ARM_ERR: In VolCurve is not of a good type");
			return ARM_KO;
		}

		if(isToClone)
			newCorrelManager = (ARM_CorrelatorManager* )theCorrelManager->Clone();
		else
			newCorrelManager = theCorrelManager;

		if( newCorrelManager == NULL)
		{
			result.setMsg ("ARM_ERR: CorrelManager is null");
			return ARM_KO;
		}

		newCorrelManager->BumpVolatility(C_ccy,ARM_CorrelatorManager::MapType(TypeCorrel),value,nthLine,
										 nthCol,isCumul,isAbsolute);

		if( ObjId == -1 )
		{
			if( isToClone )
			{
				CREATE_GLOBAL_OBJECT();

				newCorrelManagerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object* )newCorrelManager);

				if(newCorrelManagerId == RET_KO)
				{
					delete newCorrelManager;
					result.setMsg("ARM_ERR: Pb with inserting object");
					return ARM_KO;
				}
				result.setLong(newCorrelManagerId);
				return ARM_OK;
			}
			else 
			{
				result.setLong(correlManagerId);
				return ARM_OK;
			}
		}
		else
		{
			if( isToClone )
			{
			//	if( ObjId != correlManagerId )
			//	{
					oldCorrelManager = (ARM_CorrelatorManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
					delete oldCorrelManager;

					LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCorrelManager, ObjId);
					
					result.setLong(ObjId);
					return ARM_OK;
		/*		}
				else
				{
					CREATE_GLOBAL_OBJECT();

					newCorrelManagerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object* )newCorrelManager);

					if(newCorrelManagerId == RET_KO)
					{
						delete newCorrelManager;
						result.setMsg("ARM_ERR: Pb with inserting object");
						return ARM_KO;
					}
					result.setLong(newCorrelManagerId);
					return ARM_OK;
				}*/
			}
			else
			{
				// trop risqué, si on fait un premier Bump d'un objet sans Cloner et après sur la même cellule on 
				// fait un bump sans clone d'un autre objet, celà va détruire le premier objet, donc 
				// on va mettre en commentaire la distruction de l'ancien objet.
				
/*				if( ObjId != correlManagerId )
				{
					oldCorrelManager = (ARM_CorrelatorManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
					delete oldCorrelManager;
				}
*/				
				result.setLong(correlManagerId);
				return ARM_OK;
			}
		}
	}
	
	catch(Exception& x)
	{
		delete newCorrelManager;

		x.DebugPrint();

		ARM_RESULT();
    }

}


extern long ARMLOCAL_BumpVolatilityCorrelManager( long correlManagerId,
												  vector<string> mktTag,
												  vector<string> intraMktTag,
												  double value,
												  long nthLine,
												  long nthCol,
												  long isCumul,
												  long isAbsolute,
												  long isToClone,
												  ARM_result& result,
												  long ObjId)
{
	long newCorrelManagerId;

	ARM_CorrelManager* theCorrelManager = NULL;
	ARM_CorrelManager* newCorrelManager = NULL;
	ARM_CorrelManager* oldCorrelManager = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		theCorrelManager = (ARM_CorrelManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlManagerId);
		if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(theCorrelManager,ARM_CORRELMANAGER) == 0 )
		{
			result.setMsg ("ARM_ERR: Correl Manager is not of a good type");
			return ARM_KO;
		}

		if(isToClone)
			newCorrelManager = (ARM_CorrelManager* )theCorrelManager->Clone();
		else
			newCorrelManager = theCorrelManager;

		if( newCorrelManager == NULL)
		{
			result.setMsg ("ARM_ERR: CorrelManager is null");
			return ARM_KO;
		}

		newCorrelManager->BumpVolatility(mktTag,
										 intraMktTag,
										 value,
										 nthLine,
										 nthCol,
										 isCumul,
										 isAbsolute);

		if( ObjId == -1 )
		{
			if( isToClone )
			{
				CREATE_GLOBAL_OBJECT();

				newCorrelManagerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object* )newCorrelManager);

				if(newCorrelManagerId == RET_KO)
				{
					delete newCorrelManager;
					result.setMsg("ARM_ERR: Pb with inserting object");
					return ARM_KO;
				}
				result.setLong(newCorrelManagerId);
				return ARM_OK;
			}
			else 
			{
				result.setLong(correlManagerId);
				return ARM_OK;
			}
		}
		else
		{
			if( isToClone )
			{
				oldCorrelManager = (ARM_CorrelManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
				delete oldCorrelManager;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCorrelManager, ObjId);
				
				result.setLong(ObjId);
				return ARM_OK;
			}
			else
			{
				result.setLong(correlManagerId);
				return ARM_OK;
			}
		}
	}
	
	catch(Exception& x)
	{
		delete newCorrelManager;

		x.DebugPrint();

		ARM_RESULT();
    }
}


extern long ARMLOCAL_GetCorrelDiag(long IdIdCorrelId,
							string& C_Tenor1,
							vector<CCString>& C_Tenor,
							ARM_result& result,
							long ObjId)
{

	ARM_IndexIndexCorrelCube* vIndexCube;
	ARM_VolCurve* outCurve = NULL;
	ARM_VolCurve* oldCurve = NULL;
	
	int tenors_size = C_Tenor.size();
	vector<string> tenors(tenors_size);
	long CurveId;


	if(CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return	ARM_KO;
	}

	CCString	msg("");

	try
	{
		vIndexCube = (ARM_IndexIndexCorrelCube*)LOCAL_PERSISTENT_OBJECTS->GetObject(IdIdCorrelId);

		if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vIndexCube, ARM_INDEX_INDEX_CORREL_CUBE) == 0)
		{
			result.setMsg ("ARM_ERR: CorrelCube is not of a good type");
			
			return	ARM_KO;
		}
		
		for( int i = 0; i < tenors_size; i++)
			tenors[i] = CCSTringToSTLString(C_Tenor[i]);

		outCurve = vIndexCube->GetCorrelDiag(C_Tenor1,tenors);
		if(outCurve == NULL)
		{
			result.setMsg("ARM_ERR: CurveDiagCorrel is null");

		   return	ARM_KO;
		}

		if( ObjId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			CurveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object* )outCurve);

			if(CurveId == RET_KO)
			{
				delete outCurve;
				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}
			result.setLong(CurveId);
			return ARM_OK;
		}
		else
		{
			oldCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
			delete oldCurve;

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)outCurve, ObjId);
			
			result.setLong(ObjId);
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete outCurve;

		x.DebugPrint();

		ARM_RESULT();
    }

}

extern long ARMLOCAL_CreateGenCorrelatorManager(VECTOR<CCString>& mktTags,
										 vector<long>& HyperDiagVolIds,
										 vector<long>& IndexIndexVolIds, 
										 vector<long>& CorrelVolIds,
										 vector<long>& IndexVolIds,
										 vector<long>& IRVolHyperCubeIds,
										 vector<long>& VolVolHyperCubeIds,
										 vector<long>& FXVolHyperCubeIds,
										 ARM_result& result, 
										 long objId)
{
	long CorrelatorManagerId;
	int mkt_size = mktTags.size();
	int IndexIndexId_size = IndexIndexVolIds.size();
	int HyperDiagId_size = HyperDiagVolIds.size();
	int CorrelId_size = CorrelVolIds.size();
	int IndexId_size = IndexVolIds.size();
	int IRVolId_size = IRVolHyperCubeIds.size();
	int VolVolId_size = VolVolHyperCubeIds.size();
	int FXVolId_size = FXVolHyperCubeIds.size();
	ARM_CorrelatorManager* correlatorManager = NULL;
	ARM_CorrelatorManager* correlatorManagerOld = NULL;
	ARM_Object* CorrelCurve = NULL;

	if( (( mkt_size == 0 ) )
		|| ( ( HyperDiagId_size != 0 ) && ( HyperDiagId_size != mkt_size ) )
		|| ( ( IndexIndexId_size != 0 ) && ( IndexIndexId_size != mkt_size ) )
		|| ( ( CorrelId_size != 0 ) && ( CorrelId_size != mkt_size ) )
		|| ( ( IndexId_size != 0 ) && ( IndexId_size != mkt_size ) )
		|| ( ( IRVolId_size != 0 ) && ( IRVolId_size != mkt_size ) )
		|| ( ( VolVolId_size != 0 ) && ( VolVolId_size != mkt_size ) )
		|| ( ( FXVolId_size != 0 ) && ( FXVolId_size != mkt_size ) )
		)
	{
		result.setMsg ("ARM_ERR: check your matrix dimension");
		return ARM_KO;		
	}

	if( (HyperDiagId_size == 0) && (IndexIndexId_size == 0) && (CorrelId_size == 0) && (IndexId_size == 0)  &&
		(IRVolId_size == 0) && (VolVolId_size == 0) )
	{
		result.setMsg ("ARM_ERR: Empty correl manager");
		return ARM_KO;
	}

	CCString msg("");
	char MsgError[150];

	try
	{

		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		////////// Create CorrelatorManager ////////////
		correlatorManager = new ARM_CorrelatorManager();
		int i;
		for( i = 0; i < HyperDiagId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( HyperDiagVolIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( CorrelCurve, ARM_HYPER_CUBE ) == 0 )
			{

				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not Hyper Cube VolCurve",HyperDiagVolIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve, ARM_CorrelatorManager::MapType(0) );

		}

		for( i = 0; i < IndexIndexId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( IndexIndexVolIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( CorrelCurve, ARM_INDEX_INDEX_CORREL_CUBE ) == 0 )
			{
				
				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not IndexIndex Cube VolCurve",IndexIndexVolIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve,ARM_CorrelatorManager::MapType(1) );

		}

		for( i = 0; i < CorrelId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( CorrelVolIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK( CorrelCurve, ARM_VOL_CURVE ) == 0 )
			{
				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not VolCurve",CorrelVolIds[i]);
				result.setMsg( MsgError );
	

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve,ARM_CorrelatorManager::MapType(2) );

		}

		for( i = 0; i < IndexId_size; i++ )
		{				
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( IndexVolIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK( CorrelCurve, ARM_VOL_CURVE ) == 0 )
			{
				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not VolCurve",IndexVolIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve,ARM_CorrelatorManager::MapType(3));

		}

		for( i = 0; i < IRVolId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( IRVolHyperCubeIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( CorrelCurve, ARM_HYPER_CUBE ) == 0 )
			{

				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not an Hyper Cube",IRVolHyperCubeIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve, ARM_CorrelatorManager::MapType(4) );

		}

		for( i = 0; i < VolVolId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( VolVolHyperCubeIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( CorrelCurve, ARM_HYPER_CUBE ) == 0 )
			{

				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not an Hyper Cube", VolVolHyperCubeIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve, ARM_CorrelatorManager::MapType(5) );

		}

		for( i = 0; i < FXVolId_size; i++ )
		{
			CorrelCurve = LOCAL_PERSISTENT_OBJECTS->GetObject( FXVolHyperCubeIds[i] );

			if( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( CorrelCurve, ARM_VOL_CUBE ) == 0 )
			{

				sprintf(MsgError,"ARM_ERR: CorrelCurve: %ld is not an Hyper Cube", FXVolHyperCubeIds[i]);
				result.setMsg( MsgError );

				delete correlatorManager;

				return(ARM_KO);
			}
			correlatorManager->Add( (const char*)mktTags[i], CorrelCurve, ARM_CorrelatorManager::MapType(6) );

		}


		if( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			CorrelatorManagerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) correlatorManager);
			
			if ( CorrelatorManagerId == RET_KO )
			{
				if (correlatorManager)
					delete correlatorManager;
				correlatorManager = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(CorrelatorManagerId);

			return ARM_OK;	
		}
		else
		{
			correlatorManagerOld = (ARM_CorrelatorManager *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlatorManagerOld, ARM_CORRELMANAGER) == 1)
			{
				if (correlatorManagerOld)
				{
					delete correlatorManagerOld;
					correlatorManagerOld = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) correlatorManager, objId);

				return(ARM_OK);
			}
			else
			{
				if (correlatorManager)
					delete correlatorManager;
				correlatorManager = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}
	catch(Exception& x)
	{


		if (correlatorManager)
			delete correlatorManager;
		correlatorManager = NULL;


		ARM_RESULT();
	}

}


extern long ARMLOCAL_ComputeIdIdCorrelFromCorrelatorManager(const string& ccy,
															const string& tenor1,
															const string& tenor2,
															double expiry1,
															double expiry2,
															long correlManagerId,
															ARM_result& result )
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_CorrelatorManager* correlManager = NULL;
	
	CCString msg ("");

	try
	{
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		if( !GetObjectFromId( &correlManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		}

		double dResult = correlManager->ComputeIndexIndexCorrel(ccy,tenor1,tenor2,expiry1,expiry2);
		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
	//	delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_ComputeHyperCorrelFromCorrelatorManager(const string& ccy,
															 const string& tenor1,
															 const string& tenor2,
															 double expiry,
															 long correlManagerId,
															 ARM_result& result,
															 bool isByExpiry)
{

	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_CorrelatorManager* correlManager = NULL;
	
	CCString msg ("");

	try
	{
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		if( !GetObjectFromId( &correlManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		}

		double dResult;

		if( isByExpiry )
			dResult = correlManager->ComputeCorrelByExpiry(ccy,tenor1,tenor2,expiry);
		else
			dResult = correlManager->ComputeHyperCorrel(ccy,tenor1,tenor2,expiry);

		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
	//	delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_ComputeCorrSimplFromCorrelatorManager(const string& ccy,
													 const string& tenor,
													 double expiry,
													 long correlManagerId,
													 ARM_result& result,
													 bool isCorr)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_CorrelatorManager* correlManager = NULL;
	
	CCString msg ("");

	try
	{
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		if( !GetObjectFromId( &correlManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		}

		double dResult;

		if( isCorr )
			dResult = correlManager->ComputeCorrModeCorrel(ccy,tenor,expiry);
		else
			dResult = correlManager->ComputeSimpleModeCorrel(ccy,tenor,expiry);

		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
	//	delete correlManager;

		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_ComputeFXCorrelFromCorrelatorManager(const string& ccy,
														  const string& ccy2,
														  const string& tenor,
														  double expiry,
														  long correlManagerId,
														  ARM_result& result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_CorrelatorManager* correlManager = NULL;
	
	CCString msg ("");

	try
	{
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}

		if( !GetObjectFromId( &correlManager, correlManagerId, ARM_CORRELMANAGER ) )
		{
			result.setMsg ("ARM_ERR: correl manager is not of a good type");
			return ARM_KO;
		}

		double dResult = correlManager->ComputeFXCorrel(ccy, ccy2, tenor, expiry);

		result.setDouble(dResult);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_GetMixtureParamsFromSummit ( const CCString& index,
										   const CCString& currency,
										   const CCString& cvName,
										   double date,
										   const CCString& interpolMethod,
										   ARM_result& result,
										   long objId )
{
	long mixtureParamsId;

    string smileOrATM("ATM");
	CCString vType("IRG");
    string volType("CAP");
    int optionType = K_IRG; // OR K_SWOPT
    int volEnumType = K_ATMF_VOL; // K_ATMF_VOL(ATM), K_SMILE_VOL(SMILE), K_FX_VOL_SP_INTERP

	ARM_VolLInterpol* volCrv = NULL;

	ARM_ParamsMixture_Fx* paramsMixture_Fx = NULL;
	ARM_ParamsMixture_Fx* oldParamsMixture_Fx = NULL;
	
	char sDate[11];
	Local_XLDATE2ARMDATE(date, sDate);
	ARM_Date myDate(sDate);

	CCString msg(" ");
	
	try
	{
		if ( GetDataRetrieverVersion () >= ETKRETRIEVER )
		    volCrv = etoolkit_GetVolATMFromSummit(index,currency,cvName,myDate,vType);
		else
		{
			volCrv = ARMLOCAL_GetVolATMFromSummit(index,currency,cvName,myDate,vType,result);

			if ( (!volCrv) && (GetFallBackDataRetrieverVersion() >= ETKRETRIEVER))
				volCrv = etoolkit_GetVolATMFromSummit(index,currency,cvName,myDate,vType);
		}

		ARM_Vector* lags = volCrv->GetExpiryTerms();
		*lags *= 365.0;

		ARM_Vector* volATM = volCrv->GetVolatilities()->GetColumn(0);
		ARM_Vector* decVol = volCrv->GetVolatilities()->GetColumn(1);
		ARM_Vector* shift  = volCrv->GetVolatilities()->GetColumn(2);
		ARM_Vector* lambda = volCrv->GetVolatilities()->GetColumn(3);

		paramsMixture_Fx = new ARM_ParamsMixture_Fx(To_ARM_GP_Vector(*lags),
													To_ARM_GP_Vector(*volATM),
													To_ARM_GP_Vector(*decVol),
													To_ARM_GP_Vector(*shift),
													To_ARM_GP_Vector(*lambda),
													string(interpolMethod) );

		delete lags;
		delete volATM;
		delete decVol;
		delete shift;
		delete lambda;

		if ( paramsMixture_Fx == NULL )
		{
		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			mixtureParamsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)paramsMixture_Fx);

			if (mixtureParamsId == RET_KO)
			{
				if (paramsMixture_Fx)
					delete paramsMixture_Fx;
				paramsMixture_Fx = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mixtureParamsId);

			return ARM_OK;
		}
		else
		{
			oldParamsMixture_Fx = (ARM_ParamsMixture_Fx*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldParamsMixture_Fx, paramsMixture_Fx->GetName()) == 1)
			{
				if (oldParamsMixture_Fx)
				{
					delete oldParamsMixture_Fx;
					oldParamsMixture_Fx = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)paramsMixture_Fx, objId);

				return ARM_OK;
			}
			else
			{
				if (paramsMixture_Fx)
					delete paramsMixture_Fx;
				paramsMixture_Fx = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (paramsMixture_Fx)
			delete paramsMixture_Fx;
		paramsMixture_Fx = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_GetVolFromCalypso(const string& index,
									  const string& currency,
									  const string& cvName,
									  ARM_Date date,
									  const string& vtype,
									  const string& matuIndex,
									  const string& impOrHist,
									  long indexId,
									  ARM_result& result,
									  long objId ){

	

	long volId;

    string volImpOrHist("");

    if (!( impOrHist == "IRFWDVOL" ))
    {
       volImpOrHist = "HISTO";  
    }

    string smileOrATM("");

    string volType( vtype.c_str());

    int optionType = K_IRG; // OR K_SWOPT

    int volEnumType = K_ATMF_VOL; // K_ATMF_VOL(ATM), K_SMILE_VOL(SMILE), K_FX_VOL_SP_INTERP

    if ( vtype == "IRG" )
    {
       volType = "CAP";
    }
    else
    {
		if( !( volImpOrHist.empty() ) )
			volType += string(" ");

       volType = volType + volImpOrHist;

       optionType = K_SWOPT;
    }

	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* vc = NULL;

	const char* sMatuIndex = (const char*) matuIndex.c_str();


	ARM_Date myDate(date);

	CCString msg(" ");
	
	try
	{
		ARM_IRIndex* theIndex; 
		LocalPersistent::get().convert(indexId,theIndex); 
		
		if ( strcmp(sMatuIndex,"ATM") == 0 )
		{
            smileOrATM = ""; // ATM VOL
            volEnumType = K_ATMF_VOL;
        }
		else
		{
            smileOrATM = "SMILED"; 
            volEnumType = K_SMILE_VOL;
        }



		string xmlInput;
		string volSurfName=cvName+"_"+vtype+"_"+currency+"_"+index;
		ARM_CalypsoToolkit::GetVolatilitySurface(volSurfName, cvName, date,	xmlInput);
	/*	string line;
		ifstream myfile ("C:\\hbelefquih\\vol\\test_vol.xml");
		if (myfile.is_open()){
			while (! myfile.eof()){
				getline (myfile,line);
				xmlInput+=line;
			}
			myfile.close();
		}*/
		newVolCrv = GetVolSurfaceFromCalypso(date,xmlInput);

		if ( newVolCrv == NULL )
		{
		   return(ARM_KO);
		}

		if (theIndex) newVolCrv->SetIndex(*theIndex); 
        newVolCrv->SetOptionType(optionType);
        newVolCrv->SetVolType(volEnumType);

        // Update Mkt data characteristics
        string	vType = string("VOL ") + volType;
		if( !(smileOrATM.empty()) )
			vType += string(" ") + smileOrATM;

		string	vIndex(index);
		string	vCurrency( currency);
		string	vCrvId( cvName);
		
		if (( index == "ROLIB" ) || ( index == "ROEUR" ))
        {
			vIndex = "RO";
        }
		else if (( index == "NULIB" ) || ( index == "NUEUR" ))
        {
			vIndex = "NU";
        }
		else if( vIndex.substr(0, 2) == "CO" )
		{
			vType = "IRIR CORR";
			vIndex = "";
		}
 
        newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		if ( objId == -1 )
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

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vc, ARM_VOL_LIN_INTERPOL) == 1)
			{
				if (vc)
				{
					delete vc;
					vc = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVolCrv, objId);
				result.setLong(objId);
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

}


long ARMLOCAL_GetVolCubeFromCalypso(const string& index,            //EURIB
										  const string& currency,   //EUR
										  const string& cvName,     //MO
										  ARM_Date date,            //date
										  const string& vtype,      //IRG
                                          const string& type,
                                          const string suffix,
										  VECTOR<string>& tenors,
										  const string& smileOrNot,
										  long indexId,
										  ARM_result& result,
										  long objId = -1){

	long volId;

    VECTOR<CCString> smileFallBacks;

	ARM_VolCube*  createdVolCube = NULL;
	ARM_VolCube*  prevVolCube = NULL;

	
	ARM_Currency* sCCY = new ARM_Currency(currency.c_str());

	ARM_Date myDate(date);

    string volType= vtype;

    if ( vtype == "IRG" )
    {
       volType = "CAP";
    }

    // Temporary code for treating MO40 FALL BACKS

    if ( cvName == "MO40" )
    {
       smileFallBacks.push_back("MO");
    }

    CCString msg (" ");

	

	try
	{
			//	safely convert ID to ARM_IRIndex. throw if bad type,
			//	convert to 0 if -1 
			ARM_IRIndex* theIndex ; 
			LocalPersistent::get().convert(indexId,theIndex); 

		

		
		string volSurfName;
		string xmlInput;
        if (type == "SABR")
    		volSurfName=cvName+"_"+vtype+"_"+currency+"_"+index+"_SABR";
		else if(suffix !="")
            volSurfName=cvName+"_"+vtype+"_"+currency+"_"+index+"_"+suffix;
        else
            volSurfName=cvName+"_"+vtype+"_"+currency+"_"+index;

        ARM_CalypsoToolkit::GetVolatilitySurface(volSurfName, cvName, date,	xmlInput);
		
        /*string line;
												
		ifstream myfile ("C:\\hbelefquih\\vol\\test_vol.xml");
		if (myfile.is_open()){
			while (! myfile.eof()){
				getline (myfile,line);
				xmlInput+=line;
			}
			myfile.close();
		}
*/
		
		createdVolCube =GetVolCubeFromCalypso(date,cvName,xmlInput,suffix);

		if (createdVolCube == NULL)
		{
			return ARM_KO;
		}

		if (theIndex) createdVolCube->SetIndex(*theIndex); 
		createdVolCube->SetCurrencyUnit(sCCY);

       // Update Mkt data characteristics
        string	vType = string("VOL ") + volType + " SMILE";
		string	vIndex= index;
		string	vCurrency=currency;
		string	vCrvId= cvName;

        createdVolCube->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		

	

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCube);

			if ( volId == RET_KO )
			{
				if (createdVolCube)
					delete createdVolCube;
				createdVolCube = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			prevVolCube = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevVolCube, ARM_VOL_CUBE) == 1)
			{
				if (prevVolCube)
				{
					delete prevVolCube;
					prevVolCube = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdVolCube, objId);
                result.setLong(objId);
				return ARM_OK;
			}
			else
			{
				if (createdVolCube)
					delete createdVolCube;
				createdVolCube = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}catch(Exception& x)
	{
		x.DebugPrint();

	
		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (createdVolCube)
			delete createdVolCube;
		createdVolCube = NULL;

		ARM_RESULT();
	}
}



//--------------------------------------------------------------------------//
// Construction of a SABR volatility curve from 4 volatility curves :		//
// SigmaOrAlpha, Rho, Beta, Nu												//
//--------------------------------------------------------------------------//
long ARMLOCAL_SABRVol(long SigmaOrAlphaId,
					  long RhoId,
					  long BetaId,
					  long NuId,
					  long SigmaOrAlphaFlag,
					  long ModelType,
					  double Weight,
					  ARM_result& result,
					  long objId)
{
	long previousId;

	ARM_SABRVol* SABRVol	= NULL;
	ARM_SABRVol* oldSABRVol	= NULL;

	ARM_VolCurve* SigmaOrAlpha;
	ARM_VolCurve* Rho;
	ARM_VolCurve* Nu;

	ARM_VolCurve* Beta		= NULL;

	int SOrAFlag = SigmaOrAlphaFlag;
	int modType  = ModelType;

	int IRGorSWOPT, tmpIRGorSWOPT;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// Recovery of the ATM volatility curve
		//----------------------------------------------------

		SigmaOrAlpha = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(SigmaOrAlphaId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(SigmaOrAlpha, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg("ARM_ERR: SigmaOrAlpha curve is not of good type");

			return ARM_KO;
		}

		IRGorSWOPT	  = SigmaOrAlpha->GetOptionType();


		// Recovery of the Rho curve
		//----------------------------------------------------

		Rho = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(RhoId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Rho, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg("ARM_ERR: Rho curve is not of good type");

			return ARM_KO;
		}

		tmpIRGorSWOPT = Rho->GetOptionType();

		if ( tmpIRGorSWOPT != IRGorSWOPT )
		{
			result.setMsg("ARM_ERR: Sigma option type and Rho option type are not compliant");

			return ARM_KO;
		}


		// Recovery of the Beta curve
		//----------------------------------------------------

		if ( BetaId != ARM_NULL_OBJECT )
		{
			Beta = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(BetaId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Beta, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg("ARM_ERR: Beta curve is not of good type");

				return ARM_KO;
			}

			tmpIRGorSWOPT = Beta->GetOptionType();

			if ( tmpIRGorSWOPT != IRGorSWOPT )
			{
				result.setMsg("ARM_ERR: Sigma option type and Beta option type are not compliant");

				return ARM_KO;
			}
		}


		// Recovery of the Nu curve
		//----------------------------------------------------

		Nu = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(NuId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Nu, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg("ARM_ERR: Nu curve is not of good type");

			return ARM_KO;
		}

		tmpIRGorSWOPT = Nu->GetOptionType();

		if ( tmpIRGorSWOPT != IRGorSWOPT )
		{
			result.setMsg("ARM_ERR: Sigma option type and Nu option type are not compliant");

			return ARM_KO;
		}


		// Construction of the ARM_SABRVol structure
		//----------------------------------------------------

		SABRVol = new ARM_SABRVol(SigmaOrAlpha,
								  Beta,
								  Rho,
								  Nu,
								  SOrAFlag,
								  modType,
								  Weight);

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			previousId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) SABRVol);

			if ( previousId == RET_KO )
			{
				if ( SABRVol )
				{
					delete SABRVol;
					SABRVol = NULL;
				}

				result.setMsg("ARM_ERR: Pb with inserting object");
				return ARM_KO;
			}

			result.setLong(previousId);

			return ARM_OK;
		}
		else
		{
			oldSABRVol = (ARM_SABRVol*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSABRVol, ARM_VOL_SABR) == 1 )
			{
				if ( oldSABRVol )
				{
					delete oldSABRVol;
					oldSABRVol = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)SABRVol, objId);

				return ARM_OK;
			}
			else
			{
				if ( SABRVol )
				{
					delete SABRVol;
					SABRVol = NULL;
				}

				result.setMsg("ARM_ERR: previous object is not of good type");

				return ARM_KO;
			}
		}
	}

	// Catch Exceptions
	//----------------------------------------------------
	
	catch(Exception& x)
	{
		x.DebugPrint();

		if ( SABRVol )
		{
			delete SABRVol;
			SABRVol = NULL;
		}

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");

		if ( SABRVol )
		{
			delete SABRVol;
			SABRVol = NULL;
		}

		return ARM_KO;
	}
}


//--------------------------------------------------------------------------//
// Construction of a SABR volatility structure thanks to Summit				//
//--------------------------------------------------------------------------//
long ARMLOCAL_GetSABRVolFromSummit(const CCString& index,
								   const CCString& currency,
								   const CCString& cvName,
								   double date,
								   const CCString& vtype,
								   const CCString& matuIndex,
								   const CCString& impOrHist,
								   long indexId,
								   long SigmaOrAlphaFlag,
								   long ModelType,
								   double Weight,
								   ARM_result& result,
								   long objId)
{
	CCString SABRCurve[3];
	ARM_result volCurve, roCurve, nuCurve, btCurve;
	long volCurveId, roCurveId, nuCurveId, btCurveId;
	
	CCString msg ("");
	
	try
	{	
		// Two possible cases:	if currency is Euro, then choose the euribor like curves
		//						if currency is different from Euro, choose the libor like curves
		//---------------------------------------------------------------------------------------

		if ( currency == "EUR" )
		{
			SABRCurve[0] = "ROEUR";
			SABRCurve[1] = "NUEUR";
			SABRCurve[2] = "BTEUR";
		}
		else
		{
			SABRCurve[0] = "ROLIB";
			SABRCurve[1] = "NULIB";
			SABRCurve[2] = "BTLIB";
		}

		// Recovery of the volatility curves from Summit
		//---------------------------------------------------------------------------------------

		volCurveId = ARMLOCAL_GetVolFromSummit(index, currency, cvName, date, vtype, matuIndex, impOrHist, indexId, volCurve,-1);
		roCurveId  = ARMLOCAL_GetVolFromSummit(SABRCurve[0], currency, cvName, date, vtype, matuIndex, impOrHist, indexId, roCurve,-1);
		nuCurveId  = ARMLOCAL_GetVolFromSummit(SABRCurve[1], currency, cvName, date, vtype, matuIndex, impOrHist, indexId, nuCurve,-1);
		btCurveId  = ARMLOCAL_GetVolFromSummit(SABRCurve[2], currency, cvName, date, vtype, matuIndex, impOrHist, indexId, btCurve,-1);

		volCurveId = volCurve.getLong();
		roCurveId  = roCurve.getLong();
		nuCurveId  = nuCurve.getLong();
		btCurveId  = btCurve.getLong();

		// Construction of the SABR volatility 
		//---------------------------------------------------------------------------------------

		return ARMLOCAL_SABRVol(volCurveId, roCurveId, btCurveId, nuCurveId, SigmaOrAlphaFlag, ModelType, Weight, result, objId);
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_OldVolCurve(
	const long& volId,
	const double& asOf,
	ARM_result& result,
	long objId)
{
	ARM_VolCurve* vol = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Vol Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date asOfDate;
		Local_XLDATE2ARMDATE(asOf,asOfDate);

		ARM_OldVolCurve* oldVolCurve = new ARM_OldVolCurve(asOfDate,vol);

		// assign object
		if ( !assignObject( oldVolCurve, result, objId ) )
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



long ARMLOCAL_ConvIndexInYearTerm(const CCString& index,
								 double AsOf, 
								 const CCString& Ccy,
								 ARM_result& result)
{
	CCString msg ("");

	try
	{


		char sDate[11];
		Local_XLDATE2ARMDATE(AsOf,sDate);
		ARM_Date asof(sDate);

		char* sCcy = Ccy.GetStr();
		ARM_Currency ccy (sCcy);
		delete sCcy;

		char* payCalTmp = ccy.GetPayCalName(ccy.GetVanillaIndexType());
		char payCal[30];
		strcpy(payCal, payCalTmp);
		delete payCalTmp;

		long spotDays = ccy.GetSpotDays();

		ARM_Date settleDate = asof;
		settleDate.NextBusinessDay(spotDays, payCal);

		//Convert index into plot
		double d_indexterm; 
		long indexId; 
		CCString inPlot;
		indexId = ARM_ConvIrType(index); 

		if (IsCMSIndex((ARM_INDEX_TYPE)indexId) == 1)
		{
			d_indexterm = (indexId - K_CMS1 +1 );
		}	
		else
		{
			d_indexterm = 1./FromLiborTypeToFrequency(indexId); 
		}
		
		inPlot = CCString(ConvertYearTermToStringMatu(d_indexterm).c_str()); //+ CCString("Y");
		double yearTerm = convPlotInYearTerm(inPlot, asof, settleDate, payCal);

		result.setDouble(yearTerm);

		return ARM_OK; 
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    }

}
/*----- End Of File ----*/
