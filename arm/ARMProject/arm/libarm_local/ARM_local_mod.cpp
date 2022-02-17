#pragma warning(disable : 4541)
#pragma warning(disable : 4250)

#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <crv\zerocurv.h>
#include <crv\zeroflat.h>
#include <crv\volflat.h>
#include <crv\SABRVol.h>
#include <mod\ycmodel.h>
#include <mod\y2cmodel.h>
#include <mod\bsmodel.h>
#include <mod\xg2ycmod.h>
#include <mod\mcfnhw1fsv.h>
#include <mod\mcfnhw1f.h>
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
#include <mod\xbsfx.h>
#include <glob\armdef.h>
#include <mod\armfrmmodel.h>
#include <mod\armfrmmodelmixture.h>
#include <mod\armfrmmcmodelmixture.h>
#include <mod\armfrmmarkovvol.h>
#include <mod\armfrmhwvol.h>
#include <mod\armfrmmcmodel.h>
#include <mod\markovtree.h>
#include <mod\qmodel.h>
#include <crv\correlmanager.h>
#include <mod\crossmodel.h>
#include <mod\bsconvadjust.h>
#include <mod\replicconvadjust.h>
#include <mod\mapconvadjust.h>
#include <mod\replicmod.h>
#include <inst\irindex.h>
#include <mod\tribsmodel.h>
#include <mod\calibratorfrm.h>
#include <mod\tribsdual.h>
#include <mod\bsxtic.h>
#include <mod\trixbsmodel.h>
#include <mod\bsconvadjustrep.h>


#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"

#include <GP_Calculators\gpcalculators\cracalculator.h>
#include <GP_Calculators\gpcalculators\cralocalcalculator.h>
#include <gpinfra/gensecurity.h>
#include <gpcalib/calibmethod.h>
#include <gpinfra/pricingmodel.h>
#include <gpinfra/pricingadviser.h>
using namespace ARM;

long ARMLOCAL_ycmod (long idCurve,
					 long discCurveId,
					ARM_result& result,
					long objId)
{
	long modId;

	ARM_YCModel* myCreatedModel = NULL;
	ARM_YCModel* myCreatedY2CModel = NULL;
	ARM_YCModel* prevModel = NULL;
	ARM_ZeroCurve* myCurve = NULL;
	ARM_ZeroCurve* discountCurve = NULL;
	ARM_YCModel* tmpModel = NULL;
	

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		myCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(myCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		if (discCurveId == ARM_NULL_OBJECT)
		{
			myCreatedModel = new ARM_YCModel(myCurve);

			tmpModel = myCreatedModel;

			if (myCreatedModel == NULL)
			{
				result.setMsg ("ARM_ERR: Model is null");
				return ARM_KO;
			}
		}
		else
		{
			discountCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(discCurveId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(myCurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: discount Curve is not of a good type");
				return ARM_KO;
			}

			myCreatedY2CModel = new ARM_Y2CModel(myCurve,discountCurve);

			tmpModel = myCreatedY2CModel;

			if (myCreatedY2CModel == NULL)
			{
				result.setMsg ("ARM_ERR: Model is null");
				return ARM_KO;
			}
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)tmpModel);

			if (modId == RET_KO)
			{
				if (tmpModel)
					delete tmpModel;
				tmpModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_YCModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(prevModel, ARM_MODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)tmpModel, objId);

				return ARM_OK;
			}
			else
			{
				if (tmpModel)
					delete tmpModel;
				tmpModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (tmpModel)
			delete tmpModel;
		tmpModel = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_bsmodel(double date,
					  double spot,
					  long dividend_type,
					  double dividend,
					  long discrate_type,
					  double discrate,
					  long volat_type,
					  double volat,
					  long typstk,
					  ARM_result& result,
					  long objId)
{
	long modId;

	ARM_BSModel* BSmod = NULL;
	ARM_BSModel* newBSmod = NULL;
	ARM_ZeroCurve* divid = NULL;
	ARM_ZeroCurve* drate = NULL;
	ARM_VolCurve*  vol = NULL;

	ARM_ZeroCurve* tmpdivid = NULL;
	ARM_ZeroCurve* tmpdrate = NULL;
	ARM_VolCurve*  tmpvol = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * startDate = new char[11];
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,startDate);

		if ( dividend_type == 1 )
		{
			divid = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)dividend);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(divid, ARM_ZERO_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				result.setMsg ("ARM_ERR: Dividend is not a Curve");
				return ARM_KO;
			}

		}
		else
		{
			divid = new ARM_ZeroFlat((ARM_Date) startDate, dividend);

			tmpdivid = divid;
		}


		if ( discrate_type == 1 )
		{
			drate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)discrate);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(drate, ARM_ZERO_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (tmpdivid)
					delete tmpdivid;
				tmpdivid = NULL;

				result.setMsg ("ARM_ERR: discount rate is not a Curve");
				return ARM_KO;
			}
		}
		else
		{
			drate = new ARM_ZeroFlat((ARM_Date) startDate, discrate);

			tmpdrate = drate;
		}

		if (volat_type == 1)
		{
			vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)volat);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (tmpdivid)
					delete tmpdivid;
				tmpdivid = NULL;

				if (tmpdrate)
					delete tmpdrate;
				tmpdrate = NULL;

				result.setMsg ("ARM_ERR: volatility is not a Vol Curve");
				return ARM_KO;
			}

		}
		else
		{
		   vol = new ARM_VolFlat((ARM_Date) startDate, volat);

		   tmpvol = vol;
		}

		newBSmod = (ARM_BSModel *) new ARM_BSModel((ARM_Date) startDate, 
												spot, divid,
												drate, vol, typstk);

		if (tmpdivid)
		   delete tmpdivid;

		if (tmpdrate)
		   delete tmpdrate;

		if (tmpvol)
		   delete tmpvol;

		if (startDate)
			delete [] startDate;

		if (newBSmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod);

			if (modId == RET_KO)
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			BSmod = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

/*
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BSmod, ARM_BSMODEL) == 1)
			{
*/
				if (BSmod)
				{
					delete BSmod;
					BSmod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod, objId);

				result.setLong(objId);

				return ARM_OK;
/*			}
			else
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
*/
		}
	}

	catch(Exception& x)
	{
		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmpvol)
			delete tmpvol;
		tmpvol = NULL;

		if (newBSmod)
			delete newBSmod;
		newBSmod = NULL;

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_NewBSModel(long yieldCurveId,
					          long volatilityId,
							  long convadjustVolatId,
                              long correlmanagerId,
							  long convadjmanagerId,
							  long spreadLockId,
                              long discountCurveId,
					          ARM_result& result,
					          long objId)
{
	long modId;

	ARM_BSModel* BSmod = NULL;
	ARM_BSModel* newBSmod = NULL;
	ARM_ZeroCurve* YielCurve = NULL;
    ARM_ZeroCurve* DiscountCurve = NULL;
	ARM_VolCurve*  Volatility = NULL;
	ARM_VolCurve*  ConvAdjustVolatility = NULL;
    ARM_CorrelManager* CorrelManager = NULL;
    ARM_VolLInterpol* CorrelCurve = NULL;
    ARM_VolFlat* CorrelFlat = NULL;
	ARM_ConvAdjustManager* ConvAdjManager = NULL;
	ARM_VolCurve* SpreadLockCurve = NULL;	

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		YielCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)yieldCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(YielCurve, ARM_ZERO_CURVE) == 0)
		{			
			result.setMsg ("ARM_ERR: YielCurve is not a Curve");
			return ARM_KO;
		}
        ARM_Date startDate = YielCurve->GetAsOfDate();
		
		Volatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)volatilityId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Volatility, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Swaption Volatility is not a Volatility Curve");
			return ARM_KO;
		}

		if (convadjustVolatId != ARM_NULL_OBJECT)
		{
			ConvAdjustVolatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)convadjustVolatId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ConvAdjustVolatility, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Caplet Volatility is not a Volatility Curve");
				return ARM_KO;
			}
		}

        if ( correlmanagerId != ARM_NULL_OBJECT )
        {
            CorrelManager = (ARM_CorrelManager *)(LOCAL_PERSISTENT_OBJECTS->GetObject(correlmanagerId));
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CorrelManager, ARM_CORRELMANAGER) == 0)
		    { 
                CorrelCurve = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)correlmanagerId);
                
                if (CorrelCurve->GetName() != ARM_VOL_LIN_INTERPOL)
                {
                    CorrelFlat = (ARM_VolFlat*) LOCAL_PERSISTENT_OBJECTS->GetObject((long)correlmanagerId);
                    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CorrelCurve, ARM_VOL_CURVE) == 0)
		            {
                        result.setMsg ("ARM_ERR: CorrelManager  is not of a good type");
			            return ARM_KO;
                    }
                    CorrelManager= new ARM_CorrelManager("X/X" , "XXX_XXX", new ARM_CorrelMatrix(CorrelFlat));
                }
                else
                {
                    /*string indexName	= ARM_IRIndex(YielCurve->GetCurrencyUnit()->GetVanillaIndexType()).GetIndexName();
				    string ccy			= YielCurve->GetCurrencyUnit()->GetCcyName();
				    string intraMktTag	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;*/
                    CorrelManager= new ARM_CorrelManager("" , "", new ARM_CorrelMatrix(CorrelCurve));
                }
                
		    }
        }

		if (convadjmanagerId != ARM_NULL_OBJECT)
		{
			ConvAdjManager = (ARM_ConvAdjustManager*) (LOCAL_PERSISTENT_OBJECTS->GetObject(convadjmanagerId));
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_MAPCONVADJUST) == 0))
		    {

			    result.setMsg ("ARM_ERR: ConAdjManager  is not of a good type");
			    return ARM_KO;
		    }
		}

		if(spreadLockId !=  ARM_NULL_OBJECT)
		{
			SpreadLockCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadLockId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(SpreadLockCurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Spread Lock Curve  is not of a good type");
			    return ARM_KO;
			}
		}

        if (discountCurveId != ARM_NULL_OBJECT)
        {
                DiscountCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)discountCurveId);
 		        if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DiscountCurve, ARM_ZERO_CURVE) == 0)
 		        {			
 			        result.setMsg ("ARM_ERR: DiscountCurve is not a Curve");
 			        return ARM_KO;
 		        }
        }
		
		newBSmod = new ARM_BSModel((ARM_Date) startDate, YielCurve, SpreadLockCurve,
                                    ConvAdjustVolatility,Volatility,CorrelManager,ConvAdjManager, DiscountCurve);

				
		if (newBSmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod);

			if (modId == RET_KO)
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			BSmod = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
/*
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BSmod, ARM_BSMODEL) == 1)
			{
*/
				if (BSmod)
				{
					delete BSmod;
					BSmod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod, objId);

				result.setLong(objId);

				return ARM_OK;
/*
			}
			else
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
*/
		}
	}

	catch(Exception& x)
	{
		if (newBSmod)
			delete newBSmod;
		newBSmod = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_BSMODELGEN(long yieldCurveId,
						 long spreadLockId,
						 long convadjustVolatId, 
						 long volatilityId,
						 long correlmanagerId,
						 long convadjmanagerId,//convxMngrId
						 long discountCurveId,
						 long correlId,
						 long cashVolId,
						 long spreadVolId,
						 long modelTypeId,
						 long spreadVolTypeId,
						 long sabrModId,
                         bool isLnVol,
						 long numSteps,
						 ARM_result& result,
						 long objId)
{
	long modId;

	ARM_BSModel* BSmod = NULL;
	ARM_BSModel* newBSmod = NULL;
	ARM_ZeroCurve* YielCurve = NULL;
    ARM_ZeroCurve* DiscountCurve = NULL;
	ARM_VolCurve*  Volatility = NULL;
	ARM_VolCurve*  ConvAdjustVolatility = NULL;
    ARM_CorrelManager* CorrelManager = NULL;
    ARM_VolLInterpol* CorrelCurve = NULL;
    ARM_VolFlat* CorrelFlat = NULL;
	ARM_ConvAdjustManager* ConvAdjManager = NULL;
	ARM_VolCurve* SpreadLockCurve = NULL;
	ARM_VolCurve* SpreadVolCurve = NULL;

	ARM_VolCurve* Correlation = NULL;
	ARM_VolCurve* CashVol = NULL;

	ARM_BSSmiledModel* sabrMod = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		YielCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)yieldCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(YielCurve, ARM_ZERO_CURVE) == 0)
		{			
			result.setMsg ("ARM_ERR: YielCurve is not a Curve");
			return ARM_KO;
		}
        ARM_Date startDate = YielCurve->GetAsOfDate();
		
		Volatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)volatilityId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Volatility, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Swaption Volatility is not a Volatility Curve");
			return ARM_KO;
		}

		if (convadjustVolatId != ARM_NULL_OBJECT)
		{
			ConvAdjustVolatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)convadjustVolatId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ConvAdjustVolatility, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Caplet Volatility is not a Volatility Curve");
				return ARM_KO;
			}
		}

        if (correlmanagerId != ARM_NULL_OBJECT)
        {
            CorrelManager = (ARM_CorrelManager *)(LOCAL_PERSISTENT_OBJECTS->GetObject(correlmanagerId));
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CorrelManager, ARM_CORRELMANAGER) == 0)
		    {
			    
                CorrelCurve = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)correlmanagerId);
                
                if (CorrelCurve->GetName() != ARM_VOL_LIN_INTERPOL)
                {
                    CorrelFlat = (ARM_VolFlat*) LOCAL_PERSISTENT_OBJECTS->GetObject((long)correlmanagerId);
                    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CorrelCurve, ARM_VOL_CURVE) == 0)
		            {
                        result.setMsg ("ARM_ERR: CorrelManager  is not of a good type");
			            return ARM_KO;
                    }
                    CorrelManager= new ARM_CorrelManager("X/X" , "XXX_XXX", new ARM_CorrelMatrix(CorrelFlat));
                }
                else
                {
                    /*string indexName	= ARM_IRIndex(YielCurve->GetCurrencyUnit()->GetVanillaIndexType()).GetIndexName();
				    string ccy			= YielCurve->GetCurrencyUnit()->GetCcyName();
				    string intraMktTag	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;*/
                    CorrelManager= new ARM_CorrelManager("" , "", new ARM_CorrelMatrix(CorrelCurve));
                }
                
		    }
        }

		if (convadjmanagerId != ARM_NULL_OBJECT)
		{
			ConvAdjManager = (ARM_ConvAdjustManager*) (LOCAL_PERSISTENT_OBJECTS->GetObject(convadjmanagerId));
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ConvAdjManager, ARM_MAPCONVADJUST) == 0))
		    {

			    result.setMsg ("ARM_ERR: ConvAdjManager  is not of a good type");
			    return ARM_KO;
		    }
		}

		if(spreadLockId !=  ARM_NULL_OBJECT)
		{
			SpreadLockCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadLockId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(SpreadLockCurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Spread Lock Curve  is not of a good type");
			    return ARM_KO;
			}
		}

        if (discountCurveId != ARM_NULL_OBJECT)
        {
            DiscountCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)discountCurveId);
 		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DiscountCurve, ARM_ZERO_CURVE) == 0)
 		    {			
 				result.setMsg ("ARM_ERR: DiscountCurve is not a Curve");
 				return ARM_KO;
 		    }
        }
		
		
		if (correlId != ARM_NULL_OBJECT)
		{
			Correlation = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Correlation, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: correlations is not of a good type");
				return ARM_KO;
			}
		}

		if(cashVolId != ARM_NULL_OBJECT)
		{
			CashVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)cashVolId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CashVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Cash volatility is not a Volatility Curve");
				return ARM_KO;
			}
		}

		if(spreadVolId != ARM_NULL_OBJECT)
		{
			SpreadVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)spreadVolId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(SpreadVolCurve, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Spread volatility is not a Volatility Curve");
				return ARM_KO;
			}
		}

		if(spreadVolId == ARM_NULL_OBJECT && spreadVolTypeId == K_INPUTED)
		{
			result.setMsg ("ARM_ERR: SpreadVol is necessaire because the SpreadVolType is: Input");
			return ARM_KO;
		}

		if(sabrModId != ARM_NULL_OBJECT)
		{
			sabrMod = (ARM_BSSmiledModel*) LOCAL_PERSISTENT_OBJECTS->GetObject((long)sabrModId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sabrMod, ARM_BSSMILEDMODEL) == 0)
			{
				result.setMsg ("ARM_ERR: sabr model is not a bssmile model");
				return ARM_KO;
			}
		}

        if(isLnVol)
		    newBSmod = new ARM_BSModel((ARM_Date) startDate, YielCurve, SpreadLockCurve,
                                        ConvAdjustVolatility,Volatility,CorrelManager,ConvAdjManager, DiscountCurve,
									    Correlation, CashVol,SpreadVolCurve, modelTypeId,spreadVolTypeId,sabrMod,NULL,numSteps);
        else
		    newBSmod = new ARM_BSNorModel((ARM_Date) startDate, YielCurve, SpreadLockCurve,
                                        ConvAdjustVolatility,Volatility,CorrelManager,ConvAdjManager, DiscountCurve,
									    Correlation, CashVol,SpreadVolCurve, modelTypeId,spreadVolTypeId,sabrMod);

				
		if (newBSmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod);

			if (modId == RET_KO)
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			BSmod = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
/*
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BSmod, ARM_BSMODEL) == 1)
			{
*/
				if (BSmod)
				{
					delete BSmod;
					BSmod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newBSmod, objId);

				result.setLong(objId);

				return ARM_OK;
/*
			}
			else
			{
				if (newBSmod)
					delete newBSmod;
				newBSmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
*/
		}
	}

	catch(Exception& x)
	{
		if (newBSmod)
			delete newBSmod;
		newBSmod = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_BSPricingModel(long yieldCurveId,
							 long convexityManagerId,
							 long capModelId, 
							 long swoptModelId,
							 long correlManagerId,
							 long modelTypeId,
							 long spreadVolCurveId,
							 long discountCurveId,
                             long adjConvVolCurveId,
							 const VECTOR<double>& calibInfos,
							 const VECTOR<double>& numInfos,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_BSPricingModel* BSmod    = NULL;
	ARM_BSPricingModel* newBSmod = NULL;

	ARM_ZeroCurve* yieldCurve    = NULL;
    ARM_ZeroCurve* discountCurve = NULL;
	ARM_BSModel* capModel        = NULL;
	ARM_BSModel* swoptModel      = NULL;
    ARM_CorrelManager* correlManager      = NULL;
	ARM_ConvAdjustManager* convAdjManager = NULL;

	ARM_VolCurve* spreadVolCurve          = NULL;
    ARM_VolCurve* adjConvVol              = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}
	
	CCString msg("");

	try
	{
		yieldCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)yieldCurveId);
		
        if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(yieldCurve, ARM_ZERO_CURVE) == 0 )
		{			
		   result.setMsg ("ARM_ERR: YielCurve is not a Curve");
			
           return(ARM_KO);
		}

        ARM_Date startDate = yieldCurve->GetAsOfDate();
		
		if ( convexityManagerId != ARM_NULL_OBJECT )
		{
			convAdjManager = (ARM_ConvAdjustManager *) (LOCAL_PERSISTENT_OBJECTS->GetObject(convexityManagerId));
			
            if (( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjManager, ARM_BSCONVADJUST) == 0 )
				&& 
                ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjManager, ARM_REPLICCONVADJUST) == 0 )
				&&
                ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(convAdjManager, ARM_MAPCONVADJUST) == 0 ) 
               )
		    {
			   result.setMsg("ARM_ERR: ConvAdjManager is not of a good type");
			    
               return(ARM_KO);
		    }
		}

		if ( capModelId != ARM_NULL_OBJECT )
		{
			capModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) capModelId);

			if (( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(capModel, ARM_BSMODEL) == 0 )
                &&
                ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(capModel, ARM_BSSMILEDMODEL) == 0 )
               )
			{
			   result.setMsg("ARM_ERR: cap model is not a BS or BSSMILED model");
				
               return(ARM_KO);
			}
		}

		if ( swoptModelId != ARM_NULL_OBJECT )
		{
			swoptModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)swoptModelId);
			
            if (( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swoptModel, ARM_BSMODEL) == 0 )
                &&
                ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swoptModel, ARM_BSSMILEDMODEL) == 0 )
               )
			{
			   result.setMsg("ARM_ERR: swaption model is not a BS or BSSMILED model");
				
               return(ARM_KO);
			}
		}

        if ( correlManagerId != ARM_NULL_OBJECT )
        {
            correlManager = (ARM_CorrelManager *)(LOCAL_PERSISTENT_OBJECTS->GetObject(correlManagerId));

		    if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(correlManager, ARM_CORRELMANAGER) == 0 )
		    {
                ARM_VolLInterpol* CorrelCurve = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)correlManagerId);
                
                if ( CorrelCurve->GetName() != ARM_VOL_LIN_INTERPOL )
                {
                   ARM_VolFlat* CorrelFlat = (ARM_VolFlat *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) correlManagerId);
                    
                   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CorrelCurve, ARM_VOL_CURVE) == 0 )
		           {
                      result.setMsg ("ARM_ERR: CorrelManager  is not of a good type");
			            
                      return(ARM_KO);
                   }
                    
                   correlManager= new ARM_CorrelManager("X/X" , "XXX_XXX", new ARM_CorrelMatrix(CorrelFlat));
                }
                else
                {
                    correlManager= new ARM_CorrelManager("" , "", new ARM_CorrelMatrix(CorrelCurve));
                }
		    }
        }

		if ( spreadVolCurveId != ARM_NULL_OBJECT )
		{
		   spreadVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) spreadVolCurveId);
			
           if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(spreadVolCurve, ARM_VOL_CURVE) == 0 )
           {
			  result.setMsg("ARM_ERR: Spread volatility is not a Volatility Curve");
			
              return(ARM_KO);
           }
		}

        if ( adjConvVolCurveId != ARM_NULL_OBJECT )
		{
		   adjConvVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long) adjConvVolCurveId);
			
           if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(adjConvVol, ARM_VOL_CURVE) == 0 )
           {
			  result.setMsg("ARM_ERR: Spread volatility is not a Volatility Curve");
			
              return(ARM_KO);
           }
		}

        if ( discountCurveId != ARM_NULL_OBJECT )
        {
           discountCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)discountCurveId);
 		   
           if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(discountCurve, ARM_ZERO_CURVE) == 0 )
 		   {			
 			  result.setMsg("ARM_ERR: DiscountCurve is not a Curve");
 				
              return(ARM_KO);
 		   }
        }
		
		ARM_Vector* calibInfosVect = CreateARMVectorFromVECTOR(calibInfos);
		ARM_Vector* numInfosVect   = CreateARMVectorFromVECTOR(numInfos);

		newBSmod = new ARM_BSPricingModel(startDate, 
										  yieldCurve, 
										  capModel,
										  swoptModel,
										  discountCurve,
										  spreadVolCurve,
                                          adjConvVol,
										  modelTypeId,
									      calibInfosVect,
										  numInfosVect,
										  correlManager,
										  convAdjManager);

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newBSmod);

			if ( modId == RET_KO )
			{
				if (newBSmod)
				   delete newBSmod;
				newBSmod = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(modId);

			return(ARM_OK);
		}
		else
		{
			BSmod = (ARM_BSPricingModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			
			if (BSmod)
			{
			   delete BSmod;
				
               BSmod = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newBSmod, objId);

			result.setLong(objId);

			return(ARM_OK);
		}
	}

	catch(Exception& x)
	{
		if (newBSmod)
			delete newBSmod;
		newBSmod = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_GTWOYC(double dMeanRevSpeed,
					 double dSigma,
					 long dZcId,
					 long fZcId,
					 double ratesCorr,
					 double fMeanRevSpeed,
					 double fSigma,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_ZeroCurve* dZeroCurve = NULL;
	ARM_ZeroCurve* fZeroCurve = NULL;

	ARM_G2YCModel* createdG2YCModel = NULL;;
	ARM_G2YCModel* oldG2YCModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign Zc Curve is not of a good type");
			return ARM_KO;
		}

		createdG2YCModel = new ARM_G2YCModel(dMeanRevSpeed, dSigma, 
											 dZeroCurve, fZeroCurve, 
											 fMeanRevSpeed,fSigma,ratesCorr);

		if (createdG2YCModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdG2YCModel);

			if (modId == RET_KO)
			{
				if (createdG2YCModel)
					delete createdG2YCModel;
				createdG2YCModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldG2YCModel = (ARM_G2YCModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldG2YCModel, ARM_G2YCMODEL) == 1)
			{
				if (oldG2YCModel)
				{
					delete oldG2YCModel;
					oldG2YCModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdG2YCModel, objId);

				return ARM_OK;
			}

			else
			{
				if (createdG2YCModel)
					delete createdG2YCModel;
				createdG2YCModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdG2YCModel)
			delete createdG2YCModel;
		createdG2YCModel = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GYCMODEL(long zcId,
					   double a,
					   double sigma,
					   ARM_result& result,
					   long objId)
{
	long modId;

	ARM_GYCModel*  createdMod = NULL;
	ARM_GYCModel*  mod = NULL;
	ARM_ZeroCurve* zc = NULL;
	double a_;

	if (fabs(a)<=10e-5)
		a_ = 10e-5;
	else
		a_ = a;

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
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

        createdMod = new ARM_GYCModel(a_, sigma, zc);

		if (createdMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMod);

			if (modId == RET_KO)
			{
				if (createdMod)
					delete createdMod;
				createdMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_GYCModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_GYCMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdMod)
					delete createdMod;
				createdMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}		
	}
    catch(Exception& x)
    {
		x.DebugPrint();

		if (createdMod)
			delete createdMod;
		createdMod = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_HWFNMONTECARLOSV (long zcId,
								double horizon,
								double a,
								const VECTOR<double>& sigmaDate,
								const VECTOR<double>& sigmaVal,
								long nbTraj,
								ARM_result& result,
								long objId)
{
	long modId;

    ARM_MCFNHullWhite1FSigVar* HWMCmod = NULL;
    ARM_Object *omod = NULL;
    ARM_MCFNHullWhite1FSigVar*  mod = NULL;
    ARM_ZeroCurve* zc = NULL;

    double x[ARM_NB_TERMS];
    double y[ARM_NB_TERMS];

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sHorizon = new char [11];
	char* sDate = new char [11];

	CCString msg ("");

    try
    {
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			if (sDate)
				delete sDate;
			sDate = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		if(sigmaDate.size () != sigmaVal.size ())
		{
			result.setMsg ("ARM_ERR: dates and values array must have same size");
			return ARM_KO;
		}

		int size = sigmaDate.size ();

		Local_XLDATE2ARMDATE(horizon,sHorizon);

        ARM_Date beginDate = zc->GetAsOfDate();
 
        for (int i = 0; i < size; i++)
        {
			Local_XLDATE2ARMDATE(sigmaDate[i],sDate);

            ARM_Date tmpDate(sDate);
 
            double dDate=(tmpDate.GetJulian()-beginDate.GetJulian())/K_YEAR_LEN;
 
            x[i] = dDate;

			y[i] = sigmaVal[i];
        }


		HWMCmod = new ARM_MCFNHullWhite1FSigVar(zc, (ARM_Date) sHorizon,
												a, size, x, y, nbTraj);

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (HWMCmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(HWMCmod);

			if (modId == RET_KO)
			{
				if (HWMCmod)
					delete HWMCmod;
				HWMCmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			omod = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			mod = dynamic_cast<ARM_MCFNHullWhite1FSigVar *> (omod);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_MCFNHW1FSIGVAR) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWMCmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWMCmod)
					delete HWMCmod;
				HWMCmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWMCmod)
			delete HWMCmod;
		HWMCmod = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_HWFNMONTECARLO (long zcId,
							  double horizon,
							  double a,
							  double sigma,
							  long nbTraj,
							  ARM_result& result,
							  long objId)
{
	long modId;

    ARM_MCFNHullWhite1F* HWMCmod = NULL;
    ARM_Object *omod = NULL;
    ARM_MCFNHullWhite1F* mod = NULL;
    ARM_ZeroCurve* zc = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sHorizon = new char [11];

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		HWMCmod = new ARM_MCFNHullWhite1F(zc, (ARM_Date) sHorizon,
												a, sigma, nbTraj);

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (HWMCmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(HWMCmod);

			if (modId == RET_KO)
			{
				if (HWMCmod)
					delete HWMCmod;
				HWMCmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			omod = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			mod = dynamic_cast<ARM_MCFNHullWhite1F *> (omod);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_MCFNHW1FMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWMCmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWMCmod)
					delete HWMCmod;
				HWMCmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWMCmod)
			delete HWMCmod;
		HWMCmod = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_FRMTREE_AUTO(long zcId,
						   long volId,
						   long smileId,
						   long autoModeId,
						   double horizon,
						   long fineMonthId,
						   double shapeDecay,
						   double shapeSlope,
						   double shapeAsymptote,
						   long nbFactor,
						   const VECTOR<double>& corrMatu,
						   const VECTOR<double>& corrMatrix,
						   const VECTOR<double>& corr2Matu,
						   ARM_result& result,
						   long objId)
{
	long modId;

	ARM_FrmTree1	*createdFrmTree, *mod = NULL;
	ARM_Security	*Security	= NULL;
	ARM_VolCurve	*Vol		= NULL;
	ARM_VolCurve	*Smile		= NULL;
	ARM_ZeroCurve	*zc			= NULL;    

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;

	ARM_Matrix* MCorrMatrix = NULL;

	char* sHorizon = new char[11];
	Local_XLDATE2ARMDATE(horizon,sHorizon);

	CCString msg ("");

	try
	{
		zc  = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		Vol = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(volId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Vol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Volatility is not of a good type");
			return ARM_KO;
		}

		Smile = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(smileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Smile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Smile is not of a good type");
			return ARM_KO;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];
/*
		switch(autoModeId)
		{
			case K_SWOPT     : autoMode = PT_SWOPT; break;
			case K_IRG       : autoMode = PT_IRG  ; break;
			case K_SWOPT_IRG : autoMode = PT_SWOPT_IRG ; break;
			case K_IRG_SWOPT : autoMode = PT_IRG_SWOPT ; break;
			default          : autoMode = PT_SWOPT;
		}
*/
		createdFrmTree = new ARM_FrmTree1(zc, Vol, Smile, fineMonthId, 
											autoModeId, (ARM_Date)sHorizon, shapeDecay,
											shapeSlope, shapeAsymptote, nbFactor,
											VCorr2Matu, VCorrMatu, MCorrMatrix);

		if (sHorizon)
			delete sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (createdFrmTree == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdFrmTree);

			if (modId == RET_KO)
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = dynamic_cast<ARM_FrmTree1 *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_TREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmTree, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmTree)
			delete createdFrmTree;
		createdFrmTree = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FRMTREE_AUTO_G(long zcId,
							 long zdId,
							 long volId,
							 long smileId,
							 long irgvolId,
							 long irgsmileId,
							 long autoModeId,
							 double horizon,
							 long fineMonthId,
							 double shapeDecay,
							 double shapeSlope,
							 double shapeAsymptote,
							 long nbFactor,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_FrmTree1	*createdFrmTree = NULL;
	ARM_FrmTree1	*mod		= NULL;
	ARM_Security	*Security	= NULL;
	ARM_VolCurve	*Vol		= NULL;
	ARM_VolCurve	*Smile		= NULL;
	ARM_VolCurve	*irgVol		= NULL;
	ARM_VolCurve	*irgSmile	= NULL;
	ARM_ZeroCurve	*zc			= NULL; 
	ARM_ZeroCurve	*zd			= NULL; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Vector* VCorrMatu;
	ARM_Vector* VCorr2Matu;

	ARM_Matrix* MCorrMatrix;

	char* sHorizon = new char[11];
	Local_XLDATE2ARMDATE(horizon,sHorizon);

	CCString msg ("");

	try
	{
		zc  = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Curve is not of a good type");
			return ARM_KO;
		}

		zd  = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zdId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zd, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Basis Curve is not of a good type");
			return ARM_KO;
		}

		Vol = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(volId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Vol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Volatility is not of a good type");
			return ARM_KO;
		}

		Smile = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(smileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Smile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Smile is not of a good type");
			return ARM_KO;
		}

		irgVol = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(irgvolId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgVol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Irg Volatility is not of a good type");
			return ARM_KO;
		}

		irgSmile = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(irgsmileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgSmile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Irg Smile is not of a good type");
			return ARM_KO;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

/*
		switch(autoModeId)
		{
			case K_SWOPT     : autoMode = PT_SWOPT; break;
			case K_IRG       : autoMode = PT_IRG  ; break;
			case K_SWOPT_IRG : autoMode = PT_SWOPT_IRG ; break;
			case K_IRG_SWOPT : autoMode = PT_IRG_SWOPT ; break;
			default          : autoMode = PT_SWOPT;
		}
*/
		createdFrmTree = ARM_FrmTree1::CreateFrmTree1(zc, zd, Vol, Smile, irgVol, irgSmile, fineMonthId, 
											autoModeId, (ARM_Date)sHorizon, shapeDecay,
											shapeSlope, shapeAsymptote, nbFactor,
											VCorr2Matu, VCorrMatu, MCorrMatrix);

		if (sHorizon)
			delete sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (createdFrmTree == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdFrmTree);

			if (modId == RET_KO)
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = dynamic_cast<ARM_FrmTree1 *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_TREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmTree, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmTree)
			delete createdFrmTree;
		createdFrmTree = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_DFHWSIGVARTREE (double startDate,
							  double horizon,
							  long numSteps,
							  long dZcId,
							  long fZcId,
							  double dMeanRevSpeed,
							  const VECTOR<double>& dDate,
							  const VECTOR<double>& dSigma,
							  double fMeanRevSpeed,
							  const VECTOR<double>& fDate,
							  const VECTOR<double>& fSigma,
							  long dFxCorrId,
							  long fFxCorrId,
							  long fxVolId,
							  double ratesCorr,
							  double fxSpotRate,
							  ARM_result& result,
							  long objId)
{
	long modId;

	int dsize = dDate.size();
	int fsize = fDate.size();

	if(dsize != dSigma.size ())
	{
		result.setMsg ("ARM_ERR: domestic date and sigma array must have same size");
		return ARM_KO;
	}

	if(fsize != fSigma.size ())
	{
		result.setMsg ("ARM_ERR: foreign date and sigma array must have same size");
		return ARM_KO;
	}

	double dx[ARM_NB_TERMS];
	double fx[ARM_NB_TERMS];
	double ddSigma[ARM_NB_TERMS];
	double dfSigma[ARM_NB_TERMS];

	ARM_DFHWSigVarTree* mod = NULL;
	ARM_DFHWSigVarTree* createdDFTreeMod = NULL;
	ARM_VolCurve* dfxCorrCurve = NULL;
	ARM_VolCurve* ffxCorrCurve = NULL;
	ARM_VolCurve* fxVolCurve = NULL;
	ARM_ZeroCurve* dZeroCurve = NULL;
	ARM_ZeroCurve* fZeroCurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	char* sStartDate = new char[11];
	char* sHorizon = new char[11];
	char* sTmpDate;

	VECTOR<CCString> dDate_str;
	VECTOR<CCString> fDate_str;

	CCString msg ("");

	try
	{
		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: fxVolCurve is not of a good type");
			return ARM_KO;
		}

		dfxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dfxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: dfxVolCurve is not of a good type");
			return ARM_KO;
		}

		ffxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ffxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: ffxVolCurve is not of a good type");
			return ARM_KO;
		}

		Local_XLDATE2ARMDATE(startDate,sStartDate);		
		Local_XLDATE2ARMDATE(horizon,sHorizon);		

		for(int i = 0; i < dDate.size (); i++)
		{
			Local_XLDATE2ARMDATE(dDate[i],sDate);		
			dDate_str.push_back (sDate);
		}

		for(i = 0; i < fDate.size (); i++)
		{
			Local_XLDATE2ARMDATE(fDate[i],sDate);
			fDate_str.push_back (sDate);
		}

		ARM_Date dateDeb = dZeroCurve->GetAsOfDate();

		for ( i = 0; i < dsize; i++)
		{
			sTmpDate = dDate_str[i];
			ARM_Date tmpDate(sTmpDate);

			if (sTmpDate)
				free(sTmpDate);

			double date = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;

			dx[i] = date;

			ddSigma[i] = dSigma[i];
		}

		for ( i = 0; i < fsize; i++)
		{
			sTmpDate = fDate_str[i];
			ARM_Date tmpDate(sTmpDate);

			if (sTmpDate)
				free(sTmpDate);

			double date = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;

			fx[i] = date;

			dfSigma[i] = fSigma[i];
		}

		sTmpDate = NULL;
		RealStepFct dStepSigma(dsize,dx,ddSigma);

		RealStepFct fStepSigma(fsize,fx,dfSigma);

		createdDFTreeMod = new ARM_DFHWSigVarTree((ARM_Date)sStartDate,
												  (ARM_Date)sHorizon, 
												  numSteps, dMeanRevSpeed,
												  fMeanRevSpeed,
												  dStepSigma, fStepSigma,
												  dfxCorrCurve,ffxCorrCurve,
												  fxVolCurve, ratesCorr,
												  fxSpotRate,
												  dZeroCurve, fZeroCurve);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdDFTreeMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFTreeMod);

			if (modId == RET_KO)
			{
				if (createdDFTreeMod)
					delete createdDFTreeMod;
				createdDFTreeMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFTreeMod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdDFTreeMod)
					delete createdDFTreeMod;
				createdDFTreeMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdDFTreeMod)
			delete createdDFTreeMod;
		createdDFTreeMod = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT()
	}	
}


long ARMLOCAL_DFHWBASISTREE (double startDate,
							 double horizon,
							 long numSteps,
							 long dZcId,
							 long fZcId,
							 double dMeanRevSpeed,
							 const VECTOR<double>& dDate,
							 const VECTOR<double>& dSigma,
							 double fMeanRevSpeed,
							 const VECTOR<double>& fDate,
							 const VECTOR<double>& fSigma,
							 long dFxCorrId,
							 long fFxCorrId,
							 long fxVolId,
							 double ratesCorr,
							 double fxSpotRate,
							 long dNonBasisZcId,
							 long fNonBasisZcId,
							 long SwaptionBSVolId,
							 long SwaptionBSSmileId,
							 long calageType,
							 long pfType,
							 double amin,
							 double amax,
							 double volmin,
							 double volmax,
							 ARM_result& result,
							 long objId)
{
	long modId;

	int dsize = dDate.size();
	int fsize = fDate.size();

	if(dsize != dSigma.size ())
	{
		result.setMsg ("ARM_ERR: domestic date and sigma array must have same size");
		return ARM_KO;
	}

	if(fsize != fSigma.size ())
	{
		result.setMsg ("ARM_ERR: foreign date and sigma array must have same size");
		return ARM_KO;
	}

	double dx[ARM_NB_TERMS];
	double fx[ARM_NB_TERMS];
	double ddSigma[ARM_NB_TERMS];
	double dfSigma[ARM_NB_TERMS];

	ARM_DFHWSigVarTree* mod = NULL;
	ARM_DFHWSigVarTree* createdDFTreeMod = NULL;
	ARM_VolCurve* dfxCorrCurve = NULL;
	ARM_VolCurve* ffxCorrCurve = NULL;
	ARM_VolCurve* fxVolCurve = NULL;
	ARM_ZeroCurve* dZeroCurve = NULL;
	ARM_ZeroCurve* fZeroCurve = NULL;
	ARM_ZeroCurve* dNonBasisZc = NULL;
	ARM_ZeroCurve* fNonBasisZc = NULL;
	ARM_VolCurve* dSwaptionBsVol = NULL;
	ARM_VolCube* dSwaptionBsSmile = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	char* sStartDate = new char[11];
	char* sHorizon = new char[11];
	char* sTmpDate;

	VECTOR<CCString> dDate_str;
	VECTOR<CCString> fDate_str;

	CCString msg ("");

	try
	{
		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: fxVolCurve is not of a good type");
			return ARM_KO;
		}

		dfxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dfxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: dfxVolCurve is not of a good type");
			return ARM_KO;
		}

		ffxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ffxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: ffxVolCurve is not of a good type");
			return ARM_KO;
		}

		dNonBasisZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dNonBasisZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dNonBasisZc, ARM_ZERO_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: domestic non basis Zc Curve is not of a good type");
			return ARM_KO;
		}

		fNonBasisZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fNonBasisZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fNonBasisZc, ARM_ZERO_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: foreign non basis Zc Curve is not of a good type");
			return ARM_KO;
		}

		dSwaptionBsVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(SwaptionBSVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dSwaptionBsVol, ARM_VOL_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swaption BS Vol is not of a good type");
			return ARM_KO;
		}

		dSwaptionBsSmile = (ARM_VolCube *) LOCAL_PERSISTENT_OBJECTS->GetObject(SwaptionBSSmileId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(dSwaptionBsSmile, ARM_VOL_CUBE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swaption BS Smile is not of a good type");
			return ARM_KO;
		}

		Local_XLDATE2ARMDATE(startDate,sStartDate);		
		Local_XLDATE2ARMDATE(horizon,sHorizon);		

		for(int i = 0; i < dDate.size (); i++)
		{
			Local_XLDATE2ARMDATE(dDate[i],sDate);		
			dDate_str.push_back (sDate);
		}

		for(i = 0; i < fDate.size (); i++)
		{
			Local_XLDATE2ARMDATE(fDate[i],sDate);
			fDate_str.push_back (sDate);
		}

		ARM_Date dateDeb = dZeroCurve->GetAsOfDate();

		for ( i = 0; i < dsize; i++)
		{
			sTmpDate = dDate_str[i];
			ARM_Date tmpDate(sTmpDate);

			if (sTmpDate)
				delete sTmpDate;

			double date = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;

			dx[i] = date;

			ddSigma[i] = dSigma[i];
		}

		for ( i = 0; i < fsize; i++)
		{
			sTmpDate = fDate_str[i];
			ARM_Date tmpDate(sTmpDate);

			if (sTmpDate)
				delete sTmpDate;

			double date = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;

			fx[i] = date;

			dfSigma[i] = fSigma[i];
		}

		sTmpDate = NULL;
		RealStepFct dStepSigma(dsize,dx,ddSigma);

		RealStepFct fStepSigma(fsize,fx,dfSigma);

		createdDFTreeMod = new ARM_DFHWSigVarTree((ARM_Date)sStartDate,
												  (ARM_Date)sHorizon, 
												  numSteps, dMeanRevSpeed,
												  fMeanRevSpeed,
												  dStepSigma, fStepSigma,
												  dfxCorrCurve,ffxCorrCurve,
												  fxVolCurve, ratesCorr,
												  fxSpotRate,
												  dZeroCurve, fZeroCurve,
												  dNonBasisZc, fNonBasisZc,
												  dSwaptionBsVol,dSwaptionBsSmile,
												  (ARM_DFHWSigVarTree_CalibrationAction) calageType,
												  (ARM_Calibration_PORTTYPE) pfType,
												  amin, amax, volmin, volmax);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdDFTreeMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFTreeMod);

			if (modId == RET_KO)
			{
				if (createdDFTreeMod)
					delete createdDFTreeMod;
				createdDFTreeMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFTreeMod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdDFTreeMod)
					delete createdDFTreeMod;
				createdDFTreeMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdDFTreeMod)
			delete createdDFTreeMod;
		createdDFTreeMod = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT()
	}
}



long ARMLOCAL_GetParameter(long modId,
						   long paraId,
						   ARM_result& result)
{
	ARM_Model*  mod=NULL;
	double factor;

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
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		factor = mod->GetParameter(paraId);

		result.setDouble(factor);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_FRMTREE(long FrmAnaId, 
                 double horizon, 
                 long fineMonth,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
{
	long modId;

    ARM_FrmTree1  *createdFrmTree=NULL;
	ARM_FrmTree1   *mod = NULL;
    ARM_FrmAna* FrmAnaMod = NULL;

    ARM_Vector* VCorrMatu = NULL;
    ARM_Vector* VCorr2Matu = NULL;
    ARM_Matrix* MCorrMatrix = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sHorizon = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		FrmAnaMod = dynamic_cast<ARM_FrmAna *> (LOCAL_PERSISTENT_OBJECTS->GetObject(FrmAnaId));
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FrmAnaMod, ARM_FRM_ANA) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

		createdFrmTree = new ARM_FrmTree1(FrmAnaMod, (ARM_Date) sHorizon, fineMonth,
                                          VCorr2Matu, VCorrMatu, MCorrMatrix);

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdFrmTree == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmTree);

			if (modId == RET_KO)
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_FrmTree1 *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_TREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmTree, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmTree)
			delete createdFrmTree;
		createdFrmTree = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_FRMANA(long zcId,
					 const VECTOR<CCString>& resetDates,
					 const VECTOR<double>&   spotVols,
					 long shapeType,
					 double shapeDecay,
					 double shapeSlope,
					 double shapeAsymptote,
					 long nbFactor,
					 const VECTOR<double>& corrMatu,
					 const VECTOR<double>& corrMatrix,
					 const VECTOR<double>& corr2Matu,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_FrmAna* mod = NULL;
	ARM_FrmAna* Fanamod = NULL;
	ARM_ZeroCurve* zc = NULL;
	double dx[ARM_NB_TERMS];

    ARM_Vector* dates = NULL;
	ARM_Vector* VSpotVols = NULL;
    ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;

    ARM_Matrix* MCorrMatrix = NULL;

	ARM_Date       date;

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
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date dateDeb = zc->GetAsOfDate();

		for (int i=0;i<resetDates.size();i++)
		{
			ARM_Date tmpDate(resetDates[i]);

		  //double dDate = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;
			double dDate = tmpDate.GetJulian();

			dx[i] = dDate;
		}

		dates = new ARM_Vector(resetDates.size(), dx);

		VSpotVols = CreateARMVectorFromVECTOR(spotVols);

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

        Fanamod = new ARM_FrmAna(zc, dates, VSpotVols, 0,
								shapeType, shapeDecay,
								shapeSlope, shapeAsymptote,
								nbFactor, VCorr2Matu, VCorrMatu,
								MCorrMatrix);

		if (dates)
			delete dates;
		dates = NULL;

		if (VSpotVols)
			delete VSpotVols;
		VSpotVols = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (Fanamod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Fanamod);

			if (modId == RET_KO)
			{
				if (Fanamod)
					delete Fanamod;
				Fanamod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_FrmAna *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_ANA) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Fanamod, objId);

				return ARM_OK;
			}
			else
			{
				if (Fanamod)
					delete Fanamod;
				Fanamod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (Fanamod)
			delete Fanamod;
		Fanamod = NULL;

		if (dates)
			delete dates;
		dates = NULL;

		if (VSpotVols)
			delete VSpotVols;
		VSpotVols = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_FRMANA_PORT(long zcId,
						  long portId,
						  const VECTOR<CCString>& resetDates,
						  double precision,
						  double min_paras,
						  double max_paras,
						  long   max_iters,
						  long shapeType,
						  double shapeDecay,
						  double shapeSlope,
						  double shapeAsymptote,
						  long nbFactor,
						  const VECTOR<double>& corrMatu,
						  const VECTOR<double>& corrMatrix,
						  const VECTOR<double>& corr2Matu,
						  ARM_result& result,
						  long objId)
{
	long modId;

	ARM_FrmAna* mod = NULL;
	ARM_FrmAna* Fanamod = NULL;
	ARM_ZeroCurve* zc = NULL;
	ARM_Portfolio* portfolio = NULL;  
	ARM_Security *curAsset = NULL;

    ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;

    ARM_Matrix* MCorrMatrix = NULL;

    double MatsCurve[200];
    double SigmaCurve[200];

    MEMSET(MatsCurve,  0.0, sizeof(double)*200);
    MEMSET(SigmaCurve, 0.0, sizeof(double)*200);

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
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		portfolio = (ARM_Portfolio *)LOCAL_PERSISTENT_OBJECTS->GetObject(portId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(portfolio, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: portfolio is not of a good type");
			return ARM_KO;
		}

		ARM_Date dateDeb = zc->GetAsOfDate();

		for (int i = 0; i < resetDates.size(); i++)
		{
			ARM_Date tmpDate(resetDates[i]);

			double dDate = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;

			MatsCurve[i] = dDate;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

		Fanamod = new ARM_FrmAna(zc, (ARM_Object*) portfolio, 0,
							   shapeType, shapeDecay,
							   shapeSlope, shapeAsymptote,
							   nbFactor, VCorr2Matu,
							   VCorrMatu, MCorrMatrix);


		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;
		
		if (Fanamod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}
		
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Fanamod);

			if (modId == RET_KO)
			{
				if (Fanamod)
					delete Fanamod;
				Fanamod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			int PFSize = portfolio->GetSize();
 
			for (i = 0; i < PFSize; i++)
			{
				curAsset = portfolio->GetAsset(i);
 
				//curAsset->SetModel(createdFrmAna);
				//curAsset->PrepareToPrice(dateDeb);
				curAsset->SetSettlement(dateDeb);
			}


			if (shapeType != 0)
				ComputeParasCurve(Fanamod, portfolio, resetDates.size(), MatsCurve,
								  precision, min_paras, max_paras, max_iters,
								  SigmaCurve, K_DIAG_TYPE);
			else
				ComputeParasCurve(Fanamod, portfolio, resetDates.size(), MatsCurve,
								  precision, min_paras, max_paras, max_iters,
								  SigmaCurve);

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_FrmAna *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_ANA) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)Fanamod, objId);

				if (shapeType != 0)
					ComputeParasCurve(Fanamod, portfolio, resetDates.size(), MatsCurve, precision,
										min_paras, max_paras, max_iters, SigmaCurve, K_DIAG_TYPE);
				else
					ComputeParasCurve(Fanamod, portfolio, resetDates.size(), MatsCurve, precision,
										min_paras, max_paras, max_iters, SigmaCurve);

				return ARM_OK;
			}
			else
			{
				if (Fanamod)
					delete Fanamod;
				Fanamod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (Fanamod)
			delete Fanamod;
		Fanamod = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_FRMLSMC(long FrmAnaId,
					  double horizon,
					  long nbTraj,
					  long fineMonth,
					  long mcMethod,
					  const VECTOR<double>& corrMatu,
					  const VECTOR<double>& corrMatrix,
					  const VECTOR<double>& corr2Matu,
					  ARM_result& result,
					  long objId)
{
	long modId;

	ARM_SMCFrm  *createdSMCFrm = NULL, *mod = NULL;
	ARM_FrmAna* FrmAnaMod = NULL;

	ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;

	ARM_Matrix* MCorrMatrix = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* sHorizon = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);
		
		FrmAnaMod = dynamic_cast<ARM_FrmAna *> (LOCAL_PERSISTENT_OBJECTS->GetObject(FrmAnaId));

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FrmAnaMod, ARM_FRM_ANA) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: FRMANA model is not of a good type");
			return ARM_KO;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

        createdSMCFrm = new ARM_SMCFrm(FrmAnaMod, (ARM_Date)sHorizon, nbTraj, mcMethod,
										VCorr2Matu, VCorrMatu, MCorrMatrix);

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;
		
		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdSMCFrm == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSMCFrm);

			if (modId == RET_KO)
			{
				if (createdSMCFrm)
					delete createdSMCFrm;
				createdSMCFrm = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_SMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_SMCFRM) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSMCFrm, objId);

				return ARM_OK;
			}
			else
			{
				if (createdSMCFrm)
					delete createdSMCFrm;
				createdSMCFrm = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdSMCFrm)
			delete createdSMCFrm;
		createdSMCFrm = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FRMLSMC_AUTO(long zcId,
						   long volId,
						   long smileId,
						   long intAutoMode,
						   double horizon,
						   long fineMonth,
						   long   nbTraj,
						   long mcMethod,
						   long noControl,
						   double shapeDecay,
						   double shapeSlope,
						   double shapeAsymptote,
						   long nbFactor,
						   const VECTOR<double>& corrMatu,
						   const VECTOR<double>& corrMatrix,
						   const VECTOR<double>& corr2Matu,
						   ARM_result& result,
						   long objId)
{
	long modId;

	ARM_SMCFrm  *createdFrmMc = NULL, *mod = NULL;
	ARM_Security*	Security	= NULL;
	ARM_VolCurve*	Vol			= NULL;
	ARM_VolCurve*	Smile		= NULL;
	ARM_ZeroCurve*	zc			= NULL;    

//	int autoMode;

	ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;

	ARM_Matrix* MCorrMatrix = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* sHorizon = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		Vol = dynamic_cast<ARM_VolCurve   *> (LOCAL_PERSISTENT_OBJECTS->GetObject(volId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Vol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Volatility is not of a good type");
			return ARM_KO;
		}

		Smile = dynamic_cast<ARM_VolCurve   *> (LOCAL_PERSISTENT_OBJECTS->GetObject(smileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Smile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Smile is not of a good type");
			return ARM_KO;
		}

		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];
/*
        switch(intAutoMode)
        {
            case K_SWOPT     : autoMode = AUTO_MODE::IDX; break;
            case K_IRG       : autoMode = AUTO_MODE::IRG  ; break;
            case K_SWOPT_IRG : autoMode = AUTO_MODE::SWOPT_IRG ; break;
            case K_IRG_SWOPT : autoMode = AUTO_MODE::IRG_SWOPT ; break;
            default          : autoMode = AUTO_MODE::IDX;
        }
*/
		createdFrmMc = new ARM_SMCFrm(zc, Vol, Smile, intAutoMode, (ARM_Date)sHorizon, nbTraj, mcMethod,
									shapeDecay, shapeSlope, shapeAsymptote, nbFactor, 
									VCorr2Matu, VCorrMatu, MCorrMatrix, noControl);

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;
		
		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdFrmMc == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc);

			if (modId == RET_KO)
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_SMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_SMCFRM) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmMc)
			delete createdFrmMc;
		createdFrmMc = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FRMLSMC_AUTO_2(long zcId,
							 long swoptVolId,
							 long swoptSmileId,
							 long irgVolId,
							 long irgSmileId,
							 long intAutoMode,
							 double horizon,
							 long fineMonth,
							 long nbTraj,
							 long mcMethod,
							 long noControl,
							 double shapeDecay,
							 double shapeSlope,
							 double shapeAsymptote,
							 long nbFactor,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_SMCFrm  *createdFrmMc = NULL;
	ARM_SMCFrm  *mod = NULL;
	ARM_Security* Security   = NULL;
	ARM_VolCurve*   swoptVol = NULL;
	ARM_VolCurve*   swoptSmile = NULL;
	ARM_VolCurve*   irgVol   = NULL;
	ARM_VolCurve*   irgSmile = NULL;
	ARM_ZeroCurve      *zc   = NULL;

	ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;
	ARM_Matrix* MCorrMatrix = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sHorizon = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		swoptVol      = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptVolId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptVol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: swopt Vol is not of a good type");
			return ARM_KO;
		}

		swoptSmile    = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptSmileId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptSmile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swopt Smile is not of a good type");
			return ARM_KO;
		}

		irgVol        = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgVolId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgVol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Vol is not of a good type");
			return ARM_KO;
		}

		irgSmile      = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgSmileId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgSmile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Smile is not of a good type");
			return ARM_KO;
		}

/*        switch(intAutoMode)
        {
			case K_SWOPT     : autoMode = AUTO_MODE::IDX; break;
			case K_IRG       : autoMode = AUTO_MODE::IRG  ; break;
			case K_SWOPT_IRG : autoMode = AUTO_MODE::SWOPT_IRG ; break;
			case K_IRG_SWOPT : autoMode = AUTO_MODE::IRG_SWOPT ; break;
			default          : autoMode = AUTO_MODE::IDX;
        }
*/
		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

		createdFrmMc = new ARM_SMCFrm(zc, swoptVol, swoptSmile,
									  irgVol, irgSmile, intAutoMode,
									  (ARM_Date)sHorizon, nbTraj, mcMethod,
									  shapeDecay, shapeSlope,
									  shapeAsymptote, nbFactor, 
									  VCorr2Matu, VCorrMatu,
									  MCorrMatrix, noControl);

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;
		
		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdFrmMc == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc);

			if (modId == RET_KO)
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_SMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_SMCFRM) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmMc)
			delete createdFrmMc;
		createdFrmMc = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FRMLSMC_AUTO_G(long zcId,
							 long zdId,
							 long swoptVolId,
							 long swoptSmileId,
							 long irgVolId,
							 long irgSmileId,
							 long intAutoMode,
							 double horizon,
							 long fineMonth,
							 long nbTraj,
							 long mcMethod,
							 long noControl,
							 double shapeDecay,
							 double shapeSlope,
							 double shapeAsymptote,
							 long nbFactor,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_SMCFrm  *createdFrmMc = NULL;
	ARM_SMCFrm  *mod = NULL;
	ARM_Security* Security   = NULL;
	ARM_VolCurve*   swoptVol = NULL;
	ARM_VolCurve*   swoptSmile = NULL;
	ARM_VolCurve*   irgVol   = NULL;
	ARM_VolCurve*   irgSmile = NULL;
	ARM_ZeroCurve      *zc   = NULL;
	ARM_ZeroCurve      *zd   = NULL;

	ARM_Vector* VCorrMatu = NULL;
	ARM_Vector* VCorr2Matu = NULL;
	ARM_Matrix* MCorrMatrix = NULL;

//	int autoMode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sHorizon = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Curve is not of a good type");
			return ARM_KO;
		}

		zd = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zdId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zd, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic BS Curve is not of a good type");
			return ARM_KO;
		}

		swoptVol      = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptVolId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptVol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: swopt Vol is not of a good type");
			return ARM_KO;
		}

		swoptSmile    = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptSmileId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptSmile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swopt Smile is not of a good type");
			return ARM_KO;
		}

		irgVol        = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgVolId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgVol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Vol is not of a good type");
			return ARM_KO;
		}

		irgSmile      = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgSmileId));
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgSmile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Smile is not of a good type");
			return ARM_KO;
		}
/*
        switch(intAutoMode)
        {
			case K_SWOPT     : autoMode = AUTO_MODE::IDX; break;
			case K_IRG       : autoMode = AUTO_MODE::IRG  ; break;
			case K_SWOPT_IRG : autoMode = AUTO_MODE::SWOPT_IRG ; break;
			case K_IRG_SWOPT : autoMode = AUTO_MODE::IRG_SWOPT ; break;
			default          : autoMode = AUTO_MODE::IDX;
        }
*/
		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corrMatu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

		createdFrmMc = ARM_SMCFrm::CreateSmcFrm(zc, zd, swoptVol, swoptSmile,
									  irgVol, irgSmile, intAutoMode,
									  (ARM_Date)sHorizon, nbTraj, mcMethod,
									  shapeDecay, shapeSlope,
									  shapeAsymptote, nbFactor, 
									  VCorr2Matu, VCorrMatu,
									  MCorrMatrix, noControl);

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;
		
		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (createdFrmMc == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc);

			if (modId == RET_KO)
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = dynamic_cast<ARM_SMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_SMCFRM) == 1)
				|| (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BMCFRM) == 1) )
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmMc, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmMc)
					delete createdFrmMc;
				createdFrmMc = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmMc)
			delete createdFrmMc;
		createdFrmMc = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_bsslmodel(double startDate,
						long zc_type,
						double zc,
						long volSpreadLock_type,
						double volSpreadLock,
						long capVol_type,
						double capVol,
						long indexVol_type,
						double indexVol,
						ARM_result& result,
						long objId)
{
	long modId;

	ARM_BSModel* createdBSSLModel = NULL;
	ARM_BSModel* mod = NULL;

	ARM_ZeroCurve* ZcCurve = NULL;
	ARM_ZeroCurve* tmpZc = NULL;

	ARM_VolCurve* volSL = NULL;
	ARM_VolCurve* tmpVolSL = NULL;

	ARM_VolCurve*  cvCapVol = NULL;
	ARM_VolCurve*  tmpCvCapVol = NULL;

	ARM_VolCurve*  cvIndexVol = NULL;
	ARM_VolCurve*  tmpCvIndexVol = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sStartDate = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);

		// si type = 1 (obj) on recherche objet
		if ( zc_type == 1 )
		{
			ZcCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(zc));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZcCurve, ARM_ZERO_CURVE) == 0)
			{
				if (sStartDate)
					delete [] sStartDate;
				sStartDate = NULL;

				result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ZcCurve = new ARM_ZeroFlat((ARM_Date) sStartDate, zc);

			tmpZc = ZcCurve;
		}

		if ( volSpreadLock_type == 1 )
		{
			volSL = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(volSpreadLock));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volSL, ARM_VOL_CURVE) == 0)
			{
				if (sStartDate)
					delete [] sStartDate;
				sStartDate = NULL;
				
				if (tmpZc)
					delete tmpZc;
				tmpZc = NULL;

				result.setMsg ("ARM_ERR: volSpreadLock is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			volSL = new ARM_VolFlat((ARM_Date) sStartDate, volSpreadLock);

			tmpVolSL = volSL;
		}

		if ( capVol_type == 1 )
		{
			cvCapVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(capVol));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(cvCapVol, ARM_VOL_CURVE) == 0)
			{
				if (sStartDate)
					delete [] sStartDate;
				sStartDate = NULL;
				
				if (tmpZc)
					delete tmpZc;
				tmpZc = NULL;

				if (tmpVolSL)
					delete tmpVolSL;
				tmpVolSL = NULL;

				result.setMsg ("ARM_ERR: cvCapVol is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			cvCapVol = new ARM_VolFlat((ARM_Date) sStartDate, capVol);

			tmpCvCapVol = cvCapVol;
		}

		if ( indexVol_type == 1 )
		{
			cvIndexVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(indexVol));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(cvIndexVol, ARM_VOL_CURVE) == 0)
			{
				if (sStartDate)
					delete [] sStartDate;
				sStartDate = NULL;
				
				if (tmpZc)
					delete tmpZc;
				tmpZc = NULL;

				if (tmpVolSL)
					delete tmpVolSL;
				tmpVolSL = NULL;

				if (tmpCvCapVol)
					delete tmpCvCapVol;
				tmpCvCapVol = NULL;

				result.setMsg ("ARM_ERR: cvIndexVol is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			cvIndexVol = new ARM_VolFlat((ARM_Date) sStartDate, indexVol);

			tmpCvIndexVol = cvIndexVol;
		}

		createdBSSLModel = new ARM_BSModel((ARM_Date) sStartDate,ZcCurve, volSL,cvCapVol,cvIndexVol);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;
		
		if (tmpZc)
			delete tmpZc;
		tmpZc = NULL;

		if (tmpVolSL)
			delete tmpVolSL;
		tmpVolSL = NULL;

		if (tmpCvCapVol)
			delete tmpCvCapVol;
		tmpCvCapVol = NULL;

		if (tmpCvCapVol)
			delete tmpCvCapVol;
		tmpCvCapVol = NULL;

		if (createdBSSLModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBSSLModel);

			if (modId == RET_KO)
			{
				if (createdBSSLModel)
					delete createdBSSLModel;
				createdBSSLModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
		    mod = (ARM_BSModel *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BSMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBSSLModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdBSSLModel)
					delete createdBSSLModel;
				createdBSSLModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;
		
		if (tmpZc)
			delete tmpZc;
		tmpZc = NULL;

		if (tmpVolSL)
			delete tmpVolSL;
		tmpVolSL = NULL;

		if (tmpCvCapVol)
			delete tmpCvCapVol;
		tmpCvCapVol = NULL;

		if (tmpCvCapVol)
			delete tmpCvCapVol;
		tmpCvCapVol = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_BKIRTREE (long zcId,
						double startDate,
						double endDate,
						long pas,
						double a,
						double sigma,
						ARM_result& result,
						long objId)
{
	long modId;

	ARM_BKIRTree* BKTmod = NULL;
	ARM_BKIRTree*  mod = NULL;
	ARM_ZeroCurve* zc = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sStartDate = new char[11];
	char * sEndDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		BKTmod = (ARM_BKIRTree *) new ARM_BKIRTree((ARM_Date) sStartDate,
													(ARM_Date) sEndDate,
													pas,
													a,
													sigma,
													zc);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (BKTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)BKTmod);

			if (modId == RET_KO)
			{
				if (BKTmod)
					delete BKTmod;
				BKTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_BKIRTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BKIRTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)BKTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (BKTmod)
					delete BKTmod;
				BKTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (BKTmod)
			delete BKTmod;
		BKTmod = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;
		
		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_LDCMC_FROM_ANA(long anaModId,
							 double horizon,
							 long nbTraj,
							 long mcmethod,
							 long pricerTypeId,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_MCRNLogDecal* createdMCRNLogDecal=NULL;
	ARM_MCRNLogDecal* precMCRNLogDecal=NULL;
	ARM_LogDecalANA* LDC_ANA = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sHorizon = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		LDC_ANA = dynamic_cast<ARM_LogDecalANA *> (LOCAL_PERSISTENT_OBJECTS->GetObject(anaModId));

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(LDC_ANA, ARM_LOGDECALANA) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

        createdMCRNLogDecal = new ARM_MCRNLogDecal(LDC_ANA, (ARM_Date)sHorizon,
                                                   nbTraj, mcmethod, (ARM_PRICER_TYPE) pricerTypeId);

		if (createdMCRNLogDecal == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMCRNLogDecal);

			if (modId == RET_KO)
			{
				if (createdMCRNLogDecal)
					delete createdMCRNLogDecal;
				createdMCRNLogDecal = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			precMCRNLogDecal = dynamic_cast<ARM_MCRNLogDecal *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(precMCRNLogDecal, ARM_MCRNLOGDECAL) == 1)
			{
				if (precMCRNLogDecal)
				{
					delete precMCRNLogDecal;
					precMCRNLogDecal = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMCRNLogDecal, objId);

				result.setLong(objId);

				return ARM_OK;
			}
			else
			{
				if (createdMCRNLogDecal)
					delete createdMCRNLogDecal;
				createdMCRNLogDecal = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdMCRNLogDecal)
			delete createdMCRNLogDecal;
		createdMCRNLogDecal = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;
		
		ARM_RESULT();
	}
}



long ARMLOCAL_HWTREE(long zcId,
					 double begDate,
					 double endDate,
					 long nbSteps,
					 double a,
					 double sigma,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_ZeroCurve* zc=NULL;
	ARM_IRTreeHW* HWTmod=NULL;
	ARM_IRTreeHW*  mod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sStartDate = new char[11];
	char * sEndDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(begDate,sStartDate);		
		Local_XLDATE2ARMDATE(endDate,sEndDate);

	    zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		HWTmod = new ARM_IRTreeHW((ARM_Date) sStartDate,
										(ARM_Date) sEndDate,
										nbSteps,
										a,
										sigma,
										zc);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (HWTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod);

			if (modId == RET_KO)
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_IRTreeHW *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_IRTREEHW) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWTmod)
			delete HWTmod;
		HWTmod = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_HWSIGCONST(long zcId,
						 double begDate,
						 double endDate,
						 long nbSteps,
						 double a,
						 double sigma,
						 ARM_result& result,
						 long objId)
{
	long modId;

	double a_;

	ARM_ZeroCurve* zc=NULL;
	ARM_HWSigCstTree* HWTmod=NULL;
	ARM_HWSigCstTree*  mod=NULL;

	if (fabs(a)<=10e-5)
		a_ = 10e-5;
	else
		a_ = a;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sStartDate = new char[11];
	char * sEndDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(begDate,sStartDate);		
		Local_XLDATE2ARMDATE(endDate,sEndDate);

	    zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		HWTmod = new ARM_HWSigCstTree((ARM_Date) sStartDate,
									 (ARM_Date) sEndDate,
									 nbSteps, a_, sigma,
									 zc);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (HWTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod);

			if (modId == RET_KO)
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_HWSigCstTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_HWSIGCSTTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWTmod)
			delete HWTmod;
		HWTmod = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		ARM_RESULT();
	}
}




long ARMLOCAL_HWSIGVAR(long zcId,
					   long begDate,
					   long endDate,
					   long nbSteps,
					   double a,
					   long real_size,
					   const VECTOR<CCString>& dates,
					   const VECTOR<double>& sigmas,
					   ARM_result& result,
					   long objId)
{
	long modId;

	ARM_HWSigVarTree* HWTmod=NULL;
	ARM_HWSigVarTree* mod=NULL;
	ARM_ZeroCurve* zc=NULL;
	double a_;

	if (fabs(a)<=10e-5)
	   a_ = 10e-5;
	else
	   a_ = a;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	double* frac = new double [real_size];
	double* y = new double [real_size];

	char * sStartDate = new char[11];
	char * sEndDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(begDate,sStartDate);		
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			if (frac)
				delete [] frac;
			frac = NULL;

			if (y)
				delete [] y;
			y = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date beginDate(sStartDate);

		char* sDate;

		for (int i = 0; i < real_size; i++)
		{
			sDate = dates[i];
			ARM_Date tmpDate(sDate);
			if (sDate)
				free(sDate);
			double dDate = (tmpDate.GetJulian()-beginDate.GetJulian())/K_YEAR_LEN;
			frac[i] = dDate;
			y[i] = sigmas[i];
		}

		sDate = NULL;

		HWTmod = new ARM_HWSigVarTree((ARM_Date) sStartDate,
                                     (ARM_Date) sEndDate,
                                     nbSteps,
                                     a_,
                                     real_size, frac, y,
                                     zc);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (frac)
			delete [] frac;
		frac = NULL;

		if (y)
			delete [] y;
		y = NULL;

		if (HWTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod);

			if (modId == RET_KO)
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_HWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_HWSIGVARTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWTmod)
			delete HWTmod;
		HWTmod = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (frac)
			delete [] frac;
		frac = NULL;

		if (y)
			delete [] y;
		y = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_DFGYC (double dMeanRevSpeed,
					 double fMeanRevSpeed,
					 double dSigma,
					 double fSigma,
					 double fxCorr,
					 double fxVol,
					 double ratesCorr,
					 long dZcId,
					 long fZcId,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_DFGYCModel* mod=NULL;
	ARM_DFGYCModel* DFGYCmod=NULL;
	ARM_ZeroCurve* dZeroCurve=NULL;
	ARM_ZeroCurve* fZeroCurve=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Foreign Zc Curve is not of a good type");
			return ARM_KO;
		}

		DFGYCmod = new ARM_DFGYCModel(dMeanRevSpeed, fMeanRevSpeed,
									  dSigma, fSigma, fxCorr, fxVol, ratesCorr,
									  dZeroCurve, fZeroCurve);

		if (DFGYCmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)DFGYCmod);

			if (modId == RET_KO)
			{
				if (DFGYCmod)
					delete DFGYCmod;
				DFGYCmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_DFGYCModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFGYCMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)DFGYCmod, objId);

				return ARM_OK;
			}
			else
			{
				if (DFGYCmod)
					delete DFGYCmod;
				DFGYCmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (DFGYCmod)
			delete DFGYCmod;
		DFGYCmod = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_IDTREEHW(double startDate,
					   double horizon,
					   long nbSteps,
					   double dMeanRevSpeed,
					   double fMeanRevSpeed,
					   double dSigma,
					   double fSigma,
					   double prtyCorr,
					   double prtyVol,
					   double ratesCorr,
					   long dZcId,
					   long fZcId,
					   ARM_result& result,
					   long objId)
{
	long modId;

    ARM_DFIRTreeHW* mod=NULL;
    ARM_DFIRTreeHW* IR3DTmod=NULL;
    ARM_ZeroCurve* dZeroCurve=NULL;
    ARM_ZeroCurve* fZeroCurve=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sStartDate = new char[11];
	char * sHorizon = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Foreign Zc Curve is not of a good type");
			return ARM_KO;
		}

		IR3DTmod = new ARM_DFIRTreeHW((ARM_Date) sStartDate, 
									  (ARM_Date) sHorizon, nbSteps,
									  dMeanRevSpeed, fMeanRevSpeed,
									  dSigma, fSigma, prtyCorr, prtyVol,
									  ratesCorr, dZeroCurve, fZeroCurve);

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (IR3DTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)IR3DTmod);

			if (modId == RET_KO)
			{
				if (IR3DTmod)
					delete IR3DTmod;
				IR3DTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_DFIRTreeHW *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFIRTREEHW) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)IR3DTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (IR3DTmod)
					delete IR3DTmod;
				IR3DTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (IR3DTmod)
			delete IR3DTmod;
		IR3DTmod = NULL;

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_HWANALYTICSIGVAR(long zcId,
							   double a,
							   long real_size,
							   const VECTOR<CCString>& dates,
							   const VECTOR<double>& sigmas,
							   ARM_result& result,
							   long objId)
{
	long modId;

    ARM_GYCSigVarModel* HWTmod=NULL;
    ARM_GYCSigVarModel* mod=NULL;
    ARM_ZeroCurve* zc=NULL;
    double a_;

    if (fabs(a)<=10e-5)
       a_ = 10e-5;
    else
       a_ = a;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	double* frac = new double [real_size];
	double* y = new double [real_size];

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date dateDeb = zc->GetAsOfDate();

		char* sDate;

		for (int i = 0; i < real_size; i++)
		{
			sDate = dates[i];
			ARM_Date tmpDate(sDate);
			if (sDate)
				delete sDate;
			double dDate = (tmpDate.GetJulian()-dateDeb.GetJulian())/K_YEAR_LEN;
			frac[i] = dDate;
			y[i] = sigmas[i];
		}

		sDate = NULL;

		HWTmod = new ARM_GYCSigVarModel(a_, real_size, frac, y, zc);

		if (frac)
			delete [] frac;
		frac = NULL;

		if (y)
			delete [] y;
		y = NULL;

		if (HWTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod);

			if (modId == RET_KO)
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_GYCSigVarModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_GYCSIGVARMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWTmod)
			delete HWTmod;
		HWTmod = NULL;

		if (frac)
			delete [] frac;
		frac = NULL;

		if (y)
			delete [] y;
		y = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_BASISMCFRM(long zcId,
						 long baZcId,
						 long volId,
						 long smileId,
						 long productType,
						 double horizon,
						 long nbTraj,
						 long MCGeneratorType,
						 double shapeDecay,
						 double shapeSlope,
						 double shapeAsymptote,
						 long nbFactor,
						 const VECTOR<double>& indexes,
						 const VECTOR<double>& correlatedIndex,
						 const VECTOR<double>& corrMatrix,
						 long control,
						 long seed,
						 ARM_result& result,
						 long objId)
{
	long modId;

	ARM_ZeroCurve* zCurve		= NULL;
	ARM_ZeroCurve* baZCurve		= NULL;
	ARM_VolCurve*  volCurve		= NULL;
	ARM_VolCurve*  smileCurve	= NULL;

	ARM_Vector* Vindexes		= NULL;
	ARM_Vector* VCorrelatedIndexes = NULL;

	ARM_Matrix* Mcorrelations	= NULL;

	ARM_BMCFrm*    createdModel = NULL;
	ARM_BMCFrm*    mod			= NULL;

	char * sHorizon = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zCurve   = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		baZCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(baZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(baZCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Basis Zc Curve is not of a good type");
			return ARM_KO;
		}

		volCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(volId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Vol Curve is not of a good type");
			return ARM_KO;
		}

		smileCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(smileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(smileCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Smile Curve is not of a good type");
			return ARM_KO;
		}

		Vindexes = CreateARMVectorFromVECTOR(indexes);

		VCorrelatedIndexes = CreateARMVectorFromVECTOR(correlatedIndex);

		Mcorrelations = new ARM_Matrix(indexes.size(), correlatedIndex.size());
		for (int i=0;i<indexes.size();i++)
			for (int j=0;j<correlatedIndex.size();j++)
				Mcorrelations->Elt(i,j) = corrMatrix[i*correlatedIndex.size()+j];

		createdModel = new ARM_BMCFrm(zCurve, baZCurve,
									  volCurve, smileCurve,
									  productType, (ARM_Date)sHorizon,
									  nbTraj, MCGeneratorType,
									  shapeDecay, shapeSlope,
									  shapeAsymptote, nbFactor,
									  VCorrelatedIndexes, Vindexes,
									  Mcorrelations,
									  control, seed);

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorrelatedIndexes)
			delete VCorrelatedIndexes;
		VCorrelatedIndexes = NULL;

		if (Vindexes)
			delete Vindexes;
		Vindexes = NULL;

		if (Mcorrelations)
			delete Mcorrelations;
		Mcorrelations = NULL;

		if (createdModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel);

			if (modId == RET_KO)
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = dynamic_cast<ARM_BMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BMCFRM) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdModel)
			delete createdModel;
		createdModel = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorrelatedIndexes)
			delete VCorrelatedIndexes;
		VCorrelatedIndexes = NULL;

		if (Vindexes)
			delete Vindexes;
		Vindexes = NULL;

		if (Mcorrelations)
			delete Mcorrelations;
		Mcorrelations = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_BASISMCFRM2CR(long zcId,
							long baZcId,
							long swoptvolId,
							long swoptsmileId,
							long irgvolId,
							long irgsmileId,
							long productType,
							double horizon,
							long nbTraj,
							long MCGeneratorType,
							double shapeDecay,
							double shapeSlope,
							double shapeAsymptote,
							long nbFactor,
							const VECTOR<double>& indexes,
							const VECTOR<double>& correlatedIndex,
							const VECTOR<double>& corrMatrix,
							long control,
							long seed,
							ARM_result& result,
							long objId)
{
	long modId;

	ARM_ZeroCurve* zCurve		= NULL;
	ARM_ZeroCurve* baZCurve		= NULL;
	ARM_VolCurve*  irgvolCurve		= NULL;
	ARM_VolCurve*  irgsmileCurve	= NULL;
	ARM_VolCurve*  swoptvolCurve		= NULL;
	ARM_VolCurve*  swoptsmileCurve	= NULL;

	ARM_Vector* Vindexes		= NULL;
	ARM_Vector* VCorrelatedIndexes = NULL;

	ARM_Matrix* Mcorrelations	= NULL;

	ARM_BMCFrm*    createdModel = NULL;
	ARM_BMCFrm*    mod			= NULL;

	char * sHorizon = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zCurve   = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		baZCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(baZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(baZCurve, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Basis Zc Curve is not of a good type");
			return ARM_KO;
		}

		irgvolCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgvolId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgvolCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Vol Curve is not of a good type");
			return ARM_KO;
		}

		irgsmileCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(irgsmileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(irgsmileCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: IRG Smile Curve is not of a good type");
			return ARM_KO;
		}

		swoptvolCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptvolId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptvolCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swopt Vol Curve is not of a good type");
			return ARM_KO;
		}

		swoptsmileCurve = dynamic_cast<ARM_VolCurve *> (LOCAL_PERSISTENT_OBJECTS->GetObject(swoptsmileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(swoptsmileCurve, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete [] sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Swopt Smile Curve is not of a good type");
			return ARM_KO;
		}

		Vindexes = CreateARMVectorFromVECTOR(indexes);

		VCorrelatedIndexes = CreateARMVectorFromVECTOR(correlatedIndex);

		Mcorrelations = new ARM_Matrix(indexes.size(), correlatedIndex.size());
		for (int i=0;i<indexes.size();i++)
			for (int j=0;j<correlatedIndex.size();j++)
				Mcorrelations->Elt(i,j) = corrMatrix[i*correlatedIndex.size()+j];

		createdModel = new ARM_BMCFrm(zCurve, baZCurve,
									  swoptvolCurve, swoptsmileCurve,
									  irgvolCurve, irgsmileCurve,
									  productType, (ARM_Date)sHorizon,
									  nbTraj, MCGeneratorType,
									  shapeDecay, shapeSlope,
									  shapeAsymptote, nbFactor,
									  VCorrelatedIndexes, Vindexes,
									  Mcorrelations,
									  control, seed);

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorrelatedIndexes)
			delete VCorrelatedIndexes;
		VCorrelatedIndexes = NULL;

		if (Vindexes)
			delete Vindexes;
		Vindexes = NULL;

		if (Mcorrelations)
			delete Mcorrelations;
		Mcorrelations = NULL;

		if (createdModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel);

			if (modId == RET_KO)
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = dynamic_cast<ARM_BMCFrm *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BMCFRM) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdModel)
			delete createdModel;
		createdModel = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorrelatedIndexes)
			delete VCorrelatedIndexes;
		VCorrelatedIndexes = NULL;

		if (Vindexes)
			delete Vindexes;
		Vindexes = NULL;

		if (Mcorrelations)
			delete Mcorrelations;
		Mcorrelations = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_GetCalibrationOutputPFSize(long modId,
										 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		if (calHWOutput->PFSize < 0)
			return ARM_KO;

		result.setDouble(calHWOutput->PFSize);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_GetCalibrationOutputHWVolSize(long modId,
											ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		if (calHWOutput->HWVolSize < 0)
			return ARM_KO;

		result.setDouble(calHWOutput->HWVolSize);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_GetCalibrationOutputMeanRev(long modId,
										  ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->meanrev);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_GetCalibrationOutputDate(long modId,
									   long i,
									   ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();
		ARM_Date tmp = calHWOutput->datesched[i];

		char sTmp[11];
		tmp.JulianToStrDate(sTmp);

		result.setString(sTmp);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationVolSched(long modId,
									 long i,
									 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->volsched[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputSecExeDate(long modId,
											 long i,
											 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		ARM_Date tmp = calHWOutput->SecurityExeDate[i];

		char sTmp[11];
		tmp.JulianToStrDate(sTmp);

		result.setString(sTmp);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputSecStartDate(long modId,
											   long i,
											   ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		ARM_Date tmp = calHWOutput->SecurityStartDate[i];

		char sTmp[11];
		tmp.JulianToStrDate(sTmp);

		result.setString(sTmp);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_GetCalibrationOutputSecEndDate(long modId,
											 long i,
											 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		ARM_Date tmp = calHWOutput->SecurityEndDate[i];

		char sTmp[11];
		tmp.JulianToStrDate(sTmp);

		result.setString(sTmp);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_GetCalibrationOutputSecStrike(long modId,
											long i,
											ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->SecurityStrike[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputInputVol(long modId,
										   long i,
										   ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->InputVol[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputOutputVol(long modId,
											long i,
											ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->OutputVol[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputErrorVol(long modId,
										   long i,
										   ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->ErrorVol[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputInputPrice(long modId,
											 long i,
											 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->InputPrice[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputOutputPrice(long modId,
											  long i,
											  ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->OutputPrice[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}

long ARMLOCAL_GetCalibrationOutputErrorPrice(long modId,
											 long i,
											 ARM_result& result)
{
	ARM_CalibrationHWSV_OUTPUT* calHWOutput=NULL;

	ARM_DFHWSigVarTree* mod	= (ARM_DFHWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_DFHWSIGVARTREE) == 1)
	{
		calHWOutput = mod->GetCalibrationOutput();

		result.setDouble(calHWOutput->ErrorPrice[i]);

		return ARM_OK;
	}
	else
	{
		result.setMsg ("ARM_ERR: Model object is not of a good type");
		return ARM_KO;
	}
}


long ARMLOCAL_SMILEDLDCANA(long anaModId,
						   double dVolSpotUp,
						   double Asymp_RowUp,
						   double Asymp_ColUp,
						   double pUp,
						   double dVolSpotDo,
						   double Asymp_RowDo,
						   double Asymp_ColDo,
						   double pDo,
						   ARM_result& result,
						   long objId)
{
    ARM_SmiledLDCANA* createdSmiledLDCANA=NULL;
    ARM_LogDecalANA* LDC_ANA = NULL;
	int id;

    LDC_ANA = dynamic_cast<ARM_LogDecalANA *>(LOCAL_PERSISTENT_OBJECTS->GetObject(anaModId));
   

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(LDC_ANA, ARM_LOGDECALANA) == 0)
		{
			return ARM_KO;
		}

	//IS_OBJECT_CLASS_OK(LDC_ANA, ARM_LOGDECALANA, "LogDecalANA MODEL"); 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
		}

	CCString msg ("");
	
	try
	{
			//stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1) // aprés 1er appel ******************************************
			{
   
		       createdSmiledLDCANA = new ARM_SmiledLDCANA(LDC_ANA,
                                        dVolSpotUp, Asymp_RowUp,
                                        Asymp_ColUp, pUp,
                                        dVolSpotDo, Asymp_RowDo,
                                        Asymp_ColDo, pDo);

				CREATE_GLOBAL_OBJECT();

				if ( anaModId < 0 )
				{
				id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledLDCANA);
				}
				else
				{
				id = anaModId;
	           (void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledLDCANA, 
                                                    anaModId);
				}


				result.setLong(id);
				return ARM_OK;
				

			} // Fin aprés 1er appel *****************************************************
			else
			{ 
		       createdSmiledLDCANA = new ARM_SmiledLDCANA(LDC_ANA,
                                                   dVolSpotUp, Asymp_RowUp,
                                                   Asymp_ColUp, pUp,
                                                   dVolSpotDo, Asymp_RowDo,
                                                   Asymp_ColDo, pDo);

				CREATE_GLOBAL_OBJECT();

				if ( objId < 0 )
				{
				id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledLDCANA);
				}
				else
				{
				id = objId;
	           (void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledLDCANA, 
                                                    objId);
			    }

				result.setLong(id);
				return ARM_OK;
			} 
	}

    catch(Exception &x)
	{
		x.DebugPrint();

		if (createdSmiledLDCANA)
			delete createdSmiledLDCANA;
		createdSmiledLDCANA = NULL;

		ARM_RESULT();
	}

	return ARM_KO;
}


long ARMLOCAL_SMILEDLDCFROMANA(long anaModId,
							   double horizon,
							   long nbTraj,
							   long mcmethod,
							   ARM_result& result,
							   long objId)
{
    int id;

    ARM_SmiledMCRNLDC* createdSmiledMCRNLDC=NULL;
    ARM_SmiledLDCANA* SMLDC_ANA = NULL;

    //ARM_BEGIN(modTrace);

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

    SMLDC_ANA = dynamic_cast<ARM_SmiledLDCANA *>(LOCAL_PERSISTENT_OBJECTS->GetObject(anaModId));

	//IS_OBJECT_CLASS_OK(SMLDC_ANA, ARM_SMILEDLDCANA, "SmiledLogDecalANA MODEL"); 

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(SMLDC_ANA, ARM_SMILEDLDCANA) == 0)
	{
		return ARM_KO;
	}

	try
	{

			if (objId != -1)
			{  
				// ********************************************************
 
				ARM_Date dateFin(horizon);
			    createdSmiledMCRNLDC = new ARM_SmiledMCRNLDC(SMLDC_ANA, dateFin, nbTraj, mcmethod);
	
		        CREATE_GLOBAL_OBJECT();
	
				if ( objId < 0 )
				{
				id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledMCRNLDC);
				}
				else
				{
				id = objId;
 
				(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledMCRNLDC, 
                                                    objId);
				}

				// ********************************************************

				result.setLong(id);
				return ARM_OK;
    
			}
			else
			{
			ARM_Date dateFin(horizon);
		    createdSmiledMCRNLDC = new ARM_SmiledMCRNLDC(SMLDC_ANA, dateFin, nbTraj, mcmethod);

	        CREATE_GLOBAL_OBJECT();
 
			if ( anaModId < 0 )
			{
			id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledMCRNLDC);
			}
			else
		    {
			   id = anaModId;
 
			(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdSmiledMCRNLDC, 
                                                    anaModId);
			}
				result.setLong(id);
				return ARM_OK;
		}

	}

    catch(Exception &x)
	{
		x.DebugPrint();

		if (createdSmiledMCRNLDC)
			delete createdSmiledMCRNLDC;
		createdSmiledMCRNLDC = NULL;

		//ARM_RESULT();
	}

	return ARM_KO;
}



long ARMLOCAL_bssmiledmodel(double date,
							double spot,
							long dividend_type,
							double dividend,
							long discrate_type,
							double discrate,
							long volat_type,
							double volat,
							long typstk,
							const VECTOR<double>& matu,
							long rhoType,
                            long rhoObj,
                            const VECTOR<double>& rho,
                            long nuType,
                            long nuObj,
							const VECTOR<double>& nu,
							long isSABR,
							long beta_type,
							double beta,
							const VECTOR<double>& beta_vect,
							double weight,
                            long sigmaOrAlphaFlag,
                            long correlManagerId,
							long realdiscrate_type,
							double rdiscrate,
							long adjCvxVolId,
							long convToAdjVolWhenAlpha,
							ARM_result& result,
							long objId)
{
	long modId;

	int size = matu.size();

	ARM_BSSmiledModel* createdBSModel = NULL;
	ARM_BSSmiledModel* oldModel = NULL;

	ARM_ZeroCurve* divid = NULL;
	ARM_VolCurve*  vol = NULL;
	ARM_ZeroCurve* drate = NULL;
	ARM_ZeroCurve* realDiscountRate = NULL;
	ARM_ZeroCurve* tmpdivid = NULL;
	ARM_ZeroCurve* tmpdrate = NULL;
	ARM_ZeroCurve* tmprdrate = NULL;
	ARM_VolCurve*  tmpvol = NULL;
	ARM_VolCurve*  betaCurve = NULL;
	ARM_VolCurve*  tmpbeta = NULL;
	ARM_VolCurve*  adjCvxVol = NULL;

    ARM_CorrelManager* CorrelManager = NULL;

	ARM_Vector matuVect(size);
	ARM_Vector rhoVect(size);
	ARM_Vector nuVect(size);
	ARM_Vector betaVect(beta_vect.size());

    if (( rhoType == 0 ) && ( nuType == 0 ))
    {
        for (int i = 0; i < size; i++)
        {
	        matuVect.Elt(i) = matu[i];
	        rhoVect.Elt(i) = rho[i];
			nuVect.Elt(i) = nu[i];
        }
    }

	if ( beta_type == 0 )
	{
		if ( beta_vect.size() > 1 )
		{
			for (int i = 0; i < beta_vect.size(); i++)
			{
                matuVect.Elt(i) = matu[i];
				betaVect.Elt(i) = beta_vect[i];
			}
		}
		else // if ( beta_vect.size() == 1 )
		{
			betaVect.Elt(0) = beta_vect[0];
		}
	}

	char sDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

		if ( adjCvxVolId != ARM_NULL_OBJECT )
		{
			adjCvxVol = (ARM_VolCurve *) (LOCAL_PERSISTENT_OBJECTS->GetObject(adjCvxVolId));

           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(adjCvxVol, ARM_VOL_CURVE) == 0 )
		   {
              result.setMsg("ARM_ERR: A volatility curve expected");
				
              return(ARM_KO);
           }
		}

        if ( correlManagerId != ARM_NULL_OBJECT )
        {
           CorrelManager = (ARM_CorrelManager *) (LOCAL_PERSISTENT_OBJECTS->GetObject(correlManagerId));
        
           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CorrelManager, ARM_CORRELMANAGER) == 0 )
		   {
              result.setMsg("ARM_ERR: A correl Manager expected");
				
              return(ARM_KO);
           }
        }

		if ( dividend_type == 1 )
		{
			divid = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(dividend));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(divid, ARM_ZERO_CURVE) == 0)
			{
			   result.setMsg ("ARM_ERR: dividend Zc Curve is not of a good type");
				
               return(ARM_KO);
			}
		}
		else
		{
			divid = new ARM_ZeroFlat((ARM_Date) sDate, dividend);

			tmpdivid = divid;
		}

		if ( discrate_type == 1 )
		{
			drate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(discrate));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(drate, ARM_ZERO_CURVE) == 0)
			{
			   result.setMsg ("ARM_ERR: Forecast rates Zc Curve is not of a good type");
				
               return(ARM_KO);
			}
		}
		else
		{
		   drate = new ARM_ZeroFlat((ARM_Date) sDate, discrate);

		   tmpdrate = drate;
		}

		if ( volat_type == 1 )
		{
		   vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(volat));

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0 )
           {
			  result.setMsg ("ARM_ERR: Vol Curve is not of a good type");
				
              return(ARM_KO);
			}
		}
		else
		{
		   vol = new ARM_VolFlat((ARM_Date) sDate, volat);

		   tmpvol = vol;
		}

        if ( beta_type == 1 )
		{
		   if (beta != ARM_NULL_OBJECT)
		   {
			   
			   betaCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(beta));

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(betaCurve, ARM_VOL_CURVE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: Beta Curve is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else
		{
			if (betaVect.GetSize() == 1)
			{
				betaCurve = new ARM_VolFlat((ARM_Date) sDate, betaVect.Elt(0));
				tmpbeta = betaCurve;
			}
			else
			{
				if (betaVect.GetSize() == matuVect.GetSize())
				{
					ARM_Vector* yearTermsBeta = (ARM_Vector *) matuVect.Clone();
					ARM_Vector* undMatuBeta   = new ARM_Vector(1);
					ARM_Matrix* theBeta       = new ARM_Matrix(betaVect);

					betaCurve = new ARM_VolLInterpol((ARM_Date) sDate,
													 yearTermsBeta,
													 undMatuBeta,
													 theBeta);
				}
				else
				{
				   result.setMsg ("ARM_ERR: Beta and Matu are not of the same size");

				   return(ARM_KO);
				}
			}

			tmpbeta = betaCurve;
		}
		
		if (realdiscrate_type >= 0 && rdiscrate >= 0)
		{
			if ( realdiscrate_type == 1 )
			{
				realDiscountRate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(rdiscrate));

				if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(realDiscountRate, ARM_ZERO_CURVE) == 0)
				{
				   result.setMsg ("ARM_ERR: Discount rates Zc Curve is not of a good type");
					
				   return(ARM_KO);
				}
			}
			else
			{
				   realDiscountRate = new ARM_ZeroFlat((ARM_Date) sDate, rdiscrate);

				   tmprdrate = realDiscountRate;
			}
		}
		
		if (( rhoType == 0 ) && ( nuType == 0 ))
        {
		   createdBSModel = new ARM_BSSmiledModel((ARM_Date) sDate,
			  									  spot, divid, drate,
												  vol, typstk,
												  &matuVect,
												  &rhoVect,
												  &nuVect,
												  isSABR,
                                                  CorrelManager,
												  realDiscountRate,
												  adjCvxVol,
												  convToAdjVolWhenAlpha);
        }
        else if (( rhoType == 1 ) && ( nuType == 1))
        {
           ARM_VolCurve* rho = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoObj);

           if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(rho, ARM_VOL_CURVE) == 0 )
           {
			  result.setMsg ("ARM_ERR: Rho Curve is not of a good type");
				
              return(ARM_KO);
           }

           ARM_VolCurve* nu  = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuObj);
        
           if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(nu, ARM_VOL_CURVE) == 0 )
           {
			  result.setMsg ("ARM_ERR: Nu Curve is not of a good type");
				
              return(ARM_KO);
           }
           
           createdBSModel = new ARM_BSSmiledModel((ARM_Date) sDate,
			  									  spot, divid, drate,
												  vol, typstk,
												  rho,
												  nu,
												  isSABR, NULL, weight, 
                                                  sigmaOrAlphaFlag,
                                                  CorrelManager,
												  realDiscountRate,
												  adjCvxVol,
												  convToAdjVolWhenAlpha);

        }
        else
        {
           result.setMsg ("Rho, Nu input parameters problem?");

		   return(ARM_KO);
        }

		createdBSModel->SetSigmaOrAlphaInput(sigmaOrAlphaFlag);

		createdBSModel->SetBeta(betaCurve);

		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmpvol)
			delete tmpvol;
		tmpvol = NULL;

		if (tmpbeta)
			delete tmpbeta;
		tmpbeta = NULL;
		
		if (tmprdrate)
			delete tmprdrate;
		tmprdrate = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		if ( createdBSModel == NULL )
		{
		   result.setMsg ("ARM_ERR: Model is null");

		   return(ARM_KO);
		}

		CREATE_GLOBAL_OBJECT();

		if ( objId == -1 )
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) createdBSModel);

			if ( modId == RET_KO )
			{
				if (createdBSModel)
					delete createdBSModel;
				createdBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return(ARM_OK);
		}
		else
		{
			oldModel = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
/*
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldModel, ARM_BSSMILEDMODEL) == 1)
			{
*/
				if (oldModel)
				{
					delete oldModel;
					oldModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBSModel, objId);

				return ARM_OK;
/*
            }
			else
			{
				if (createdBSModel)
					delete createdBSModel;
				createdBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
*/
		}
	}

	catch(Exception& x)
	{
		// x.DebugPrint();

		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmpvol)
			delete tmpvol;
		tmpvol = NULL;

		if (tmpbeta)
			delete tmpbeta;
		tmpbeta = NULL;

		if (createdBSModel)
		   delete createdBSModel;
		createdBSModel = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		ARM_RESULT();
	}

	catch(...)
	{
		Exception x(__LINE__, __FILE__, ERR_INVALID_DATA, 
			        "Unrecognized failure in : ARMLOCAL_bssmiledmodel");

 		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmpvol)
			delete tmpvol;
		tmpvol = NULL;

		if (tmpbeta)
			delete tmpbeta;
		tmpbeta = NULL;

		if (createdBSModel)
		   delete createdBSModel;
		createdBSModel = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		ARM_RESULT();
	}
}

//------------------------------------------------------------------------------//
// BSSmiled interface with SABR volatility structure (ARM_SABRVol)				//
//------------------------------------------------------------------------------//
long ARMLOCAL_bssmiledmodel_SABR(double date,
								 double spot,
								 long dividend_type,
								 double dividend,
								 long discrate_type,
								 double discrate,
								 long SABRVolId,
								 long typstk,
								 long correlManagerId,
								 long realdiscrate_type,
								 double rdiscrate,
								 long adjCvxVolId,
								 long convToAdjVolWhenAlpha,
								 ARM_result& result,
								 long objId)
{
	long modId;

	ARM_BSSmiledModel* createdBSModel = NULL;
	ARM_BSSmiledModel* oldModel = NULL;

	ARM_ZeroCurve* divid = NULL;
	ARM_ZeroCurve* drate = NULL;
	ARM_ZeroCurve* realDiscountRate = NULL;
	ARM_ZeroCurve* tmpdivid = NULL;
	ARM_ZeroCurve* tmpdrate = NULL;
	ARM_ZeroCurve* tmprdrate = NULL;
	ARM_VolCurve*  adjCvxVol = NULL;

    ARM_CorrelManager* CorrelManager = NULL;
	ARM_SABRVol* SABRVol = NULL;

	ARM_Date StartDate;
	char sDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

		// Recovery of the SABR volatility curve
		//----------------------------------------------------------

		if ( SABRVolId != ARM_NULL_OBJECT )
		{
			SABRVol = (ARM_SABRVol *) (LOCAL_PERSISTENT_OBJECTS->GetObject(SABRVolId));

           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(SABRVol, ARM_VOL_SABR) == 0 )
		   {
              result.setMsg("ARM_ERR: A SABR volatility curve expected");
				
              return(ARM_KO);
           }
		}

		// Recovery of the convexity adjustments volatility curve
		//----------------------------------------------------------

		if ( adjCvxVolId != ARM_NULL_OBJECT )
		{
			adjCvxVol = (ARM_VolCurve *) (LOCAL_PERSISTENT_OBJECTS->GetObject(adjCvxVolId));

           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(adjCvxVol, ARM_VOL_CURVE) == 0 )
		   {
              result.setMsg("ARM_ERR: A volatility curve expected");
				
              return(ARM_KO);
           }
		}

		// Recovery of the correlation manager
		//----------------------------------------------------------

        if ( correlManagerId != ARM_NULL_OBJECT )
        {
           CorrelManager = (ARM_CorrelManager *) (LOCAL_PERSISTENT_OBJECTS->GetObject(correlManagerId));
        
           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CorrelManager, ARM_CORRELMANAGER) == 0 )
		   {
              result.setMsg("ARM_ERR: A correl Manager expected");
				
              return(ARM_KO);
           }
        }

		// Recovery of the dividend curve
		//----------------------------------------------------------

		if ( dividend_type == 1 )
		{
			divid = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(dividend));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(divid, ARM_ZERO_CURVE) == 0)
			{
			   result.setMsg ("ARM_ERR: dividend Zc Curve is not of a good type");
				
               return(ARM_KO);
			}
		}
		else
		{
			divid = new ARM_ZeroFlat((ARM_Date) sDate, dividend);

			tmpdivid = divid;
		}

		// Recovery of the forecast zero rates curve
		//----------------------------------------------------------

		if ( discrate_type == 1 )
		{
			drate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(discrate));

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(drate, ARM_ZERO_CURVE) == 0)
			{
			   result.setMsg ("ARM_ERR: Forecast rates Zc Curve is not of a good type");
				
               return(ARM_KO);
			}
		}
		else
		{
		   drate = new ARM_ZeroFlat((ARM_Date) sDate, discrate);

		   tmpdrate = drate;
		}

		// Recovery of the real zero rates curve
		//----------------------------------------------------------

		if (realdiscrate_type >= 0 && rdiscrate >= 0)
		{
			if ( realdiscrate_type == 1 )
			{
				realDiscountRate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(rdiscrate));

				if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(realDiscountRate, ARM_ZERO_CURVE) == 0)
				{
				   result.setMsg ("ARM_ERR: Discount rates Zc Curve is not of a good type");
					
				   return(ARM_KO);
				}
			}
			else
			{
				   realDiscountRate = new ARM_ZeroFlat((ARM_Date) sDate, rdiscrate);

				   tmprdrate = realDiscountRate;
			}
		}


		// Start Date
		//----------------------------------------------------------
		
		StartDate      = (ARM_Date) sDate;
	
		
		// Construction of the model 
		//----------------------------------------------------------

		createdBSModel = new ARM_BSSmiledModel(StartDate,
		   									   spot, 
											   divid, 
											   drate,
											   SABRVol,
											   typstk,
                                               CorrelManager,
											   realDiscountRate,
											   adjCvxVol,
											   convToAdjVolWhenAlpha);



		// Objects destruction
		//----------------------------------------------------------

		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmprdrate)
			delete tmprdrate;
		tmprdrate = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		if ( createdBSModel == NULL )
		{
		   result.setMsg ("ARM_ERR: Model is null");

		   return(ARM_KO);
		}

		// Global object recovery
		//----------------------------------------------------------

		CREATE_GLOBAL_OBJECT();

		if ( objId == -1 )
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) createdBSModel);

			if ( modId == RET_KO )
			{
				if (createdBSModel)
					delete createdBSModel;
				createdBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return(ARM_OK);
		}
		else
		{
			oldModel = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (oldModel)
			{
				delete oldModel;
				oldModel = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBSModel, objId);

			return ARM_OK;
		}
	}

	// Catch Exceptions
	//----------------------------------------------------------

	catch(Exception& x)
	{
		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (createdBSModel)
		   delete createdBSModel;
		createdBSModel = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		ARM_RESULT();
	}

	catch(...)
	{
		Exception x(__LINE__, __FILE__, ERR_INVALID_DATA, "Unrecognized failure in : ARMLOCAL_bssmiledmodel");

 		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (createdBSModel)
		   delete createdBSModel;
		createdBSModel = NULL;

		if (adjCvxVol)
			delete adjCvxVol;
		adjCvxVol = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_CRRTREE(double begDate,
					  double endDate,
					  long nbSteps,
					  double spot,
					  double dividend,
					  double discrate,
					  double volat,
					  ARM_result& result,
					  long objId)
{
    int modId;

    ARM_CRRTree* mod=NULL;
    ARM_CRRTree* createdCRRTmod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sDate1 = new char[11];
	char * sDate2 = new char[11];

	CCString msg ("");
	
	try
	{
		Local_XLDATE2ARMDATE(begDate,sDate1);
		Local_XLDATE2ARMDATE(endDate,sDate2);

		createdCRRTmod = new ARM_CRRTree((ARM_Date)sDate1,
										 (ARM_Date)sDate2,
										 nbSteps,
										 spot,
										 dividend,
										 discrate,
										 volat);

		if (sDate1)
			delete [] sDate1;
		sDate1 = NULL;

		if (sDate2)
			delete [] sDate2;
		sDate2 = NULL;

		if (createdCRRTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCRRTmod);

			if (modId == RET_KO)
			{
				if (createdCRRTmod)
					delete createdCRRTmod;
				createdCRRTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_CRRTree *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_CRRTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCRRTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdCRRTmod)
					delete createdCRRTmod;
				createdCRRTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdCRRTmod)
			delete createdCRRTmod;
		createdCRRTmod = NULL;

		if (sDate1)
			delete [] sDate1;
		sDate1 = NULL;

		if (sDate2)
			delete [] sDate2;
		sDate2 = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_BSCORRMODEL(double startDate,
						  long zcId,
						  long spreadlockId,
						  long capirgvolId,
						  long capcashvolId,
						  long indexadjvolId,
						  long spreadvolId,
						  long correlationsId,
						  long modelTypeId,
						  long volTypeId,
						  ARM_result& result,
						  long objId)
{
    int modId;

    ARM_BSCorrModel* mod = NULL;
    ARM_BSCorrModel* createdbscorrmod = NULL;

	ARM_ZeroCurve* zc = NULL;
	ARM_VolCurve* spreadLock = NULL;
	ARM_VolCurve* capIRGVol = NULL;
	ARM_VolCurve* capCashVol = NULL;
	ARM_VolCurve* tmpCapCashVol = NULL;
	ARM_VolCurve* indexVol = NULL;
	ARM_VolCurve* spreadVol = NULL;
	ARM_VolCurve* correlations = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sDate);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		spreadLock = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadlockId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(spreadLock, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: spread Lock is not of a good type");
			return ARM_KO;
		}

		capIRGVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capirgvolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capIRGVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Cap IRG vol is not of a good type");
			return ARM_KO;
		}

		if (capcashvolId != ARM_NULL_OBJECT)
		{
			capCashVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capcashvolId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(capCashVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Cap Cash vol is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			capCashVol = new ARM_VolFlat(zc->GetAsOfDate(),0.0);

			tmpCapCashVol = capCashVol;
		}


		if (indexadjvolId != ARM_NULL_OBJECT)
		{
			indexVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(indexadjvolId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(indexVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Index vol is not of a good type");
				return ARM_KO;
			}
		}

		if (spreadvolId != ARM_NULL_OBJECT)
		{
			spreadVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadvolId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(spreadVol, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: spread vol is not of a good type");
				return ARM_KO;
			}
		}

		if (correlationsId != ARM_NULL_OBJECT)
		{
			correlations = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(correlationsId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(correlations, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: correlations is not of a good type");
				return ARM_KO;
			}
		}

		createdbscorrmod = new ARM_BSCorrModel((ARM_Date) sDate,
												zc,
												spreadLock,
												capIRGVol,
												capCashVol,
												indexVol,
												spreadVol,
												correlations,
												modelTypeId,
												volTypeId);

		if (tmpCapCashVol)
			delete tmpCapCashVol;
		tmpCapCashVol = NULL;

		if (createdbscorrmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdbscorrmod);

			if (modId == RET_KO)
			{
				if (createdbscorrmod)
					delete createdbscorrmod;
				createdbscorrmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_BSCorrModel *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BSMODEL) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdbscorrmod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdbscorrmod)
					delete createdbscorrmod;
				createdbscorrmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdbscorrmod)
			delete createdbscorrmod;
		createdbscorrmod = NULL;

		if (tmpCapCashVol)
			delete tmpCapCashVol;
		tmpCapCashVol = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_FRMTREE_AUTO_B(long zcId,
							 long zdId,
							 long volId,
							 long smileId,
							 long autoModeId,
							 double horizon,
							 long fineMonthId,
							 double shapeDecay,
							 double shapeSlope,
							 double shapeAsymptote,
							 long nbFactor,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId)
{
	long modId;

	ARM_FrmTree1	*createdFrmTree, *mod = NULL;
	ARM_Security	*Security	= NULL;
	ARM_VolCurve	*Vol		= NULL;
	ARM_VolCurve	*Smile		= NULL;
	ARM_ZeroCurve	*zc			= NULL; 
	ARM_ZeroCurve	*zd			= NULL; 

//	ARM_PRODUCT_TYPE autoMode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Vector* VCorrMatu;
	ARM_Vector* VCorr2Matu;

	ARM_Matrix* MCorrMatrix;

	char* sHorizon = new char[11];
	Local_XLDATE2ARMDATE(horizon,sHorizon);

	CCString msg ("");

	try
	{
		zc  = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Curve is not of a good type");
			return ARM_KO;
		}

		zd  = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zdId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zd, ARM_ZERO_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Domestic Basis Curve is not of a good type");
			return ARM_KO;
		}

		Vol = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(volId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Vol, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Volatility is not of a good type");
			return ARM_KO;
		}

		Smile = dynamic_cast<ARM_VolCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(smileId));

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Smile, ARM_VOL_CURVE) == 0)
		{
			if (sHorizon)
				delete sHorizon;
			sHorizon = NULL;

			result.setMsg ("ARM_ERR: Smile is not of a good type");
			return ARM_KO;
		}


		VCorrMatu = CreateARMVectorFromVECTOR(corrMatu);

		VCorr2Matu = CreateARMVectorFromVECTOR(corr2Matu);

		MCorrMatrix = new ARM_Matrix(corrMatu.size(), corr2Matu.size());
		for (int i=0;i<corrMatu.size();i++)
			for (int j=0;j<corr2Matu.size();j++)
				MCorrMatrix->Elt(i,j) = corrMatrix[i*corr2Matu.size()+j];

/*		switch(autoModeId)
		{
			case K_SWOPT     : autoMode = PT_SWOPT; break;
			case K_IRG       : autoMode = PT_IRG  ; break;
			case K_SWOPT_IRG : autoMode = PT_SWOPT_IRG ; break;
			case K_IRG_SWOPT : autoMode = PT_IRG_SWOPT ; break;
			default          : autoMode = PT_SWOPT;
		}
*/
		createdFrmTree = ARM_FrmTree1::CreateFrmTree1(zc, zd, Vol, Smile, fineMonthId, 
											autoModeId, (ARM_Date)sHorizon, shapeDecay,
											shapeSlope, shapeAsymptote, nbFactor,
											VCorr2Matu, VCorrMatu, MCorrMatrix);

		if (sHorizon)
			delete sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		if (createdFrmTree == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdFrmTree);

			if (modId == RET_KO)
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = dynamic_cast<ARM_FrmTree1 *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_FRM_TREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFrmTree, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFrmTree)
					delete createdFrmTree;
				createdFrmTree = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdFrmTree)
			delete createdFrmTree;
		createdFrmTree = NULL;

		if (sHorizon)
			delete [] sHorizon;
		sHorizon = NULL;

		if (VCorr2Matu)
			delete VCorr2Matu;
		VCorr2Matu = NULL;

		if (VCorrMatu)
			delete VCorrMatu;
		VCorrMatu = NULL;

		if (MCorrMatrix)
			delete MCorrMatrix;
		MCorrMatrix = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_DFFXBS (long dVolId,
					  long fVolId,
					  long dZcId,
					  long fZcId,
					  long dFxCorrId,
					  long fFxCorrId,
					  long fxVolId,
					  double ratesCorr,
					  long dBSZcId,
					  long fBSZcId,
					  double spot,
                      long FundZcId,
                      long FundBSZcId,
                      double fundFxSpot,
					  long fxSmileId,
                      bool isLnVols,
					  long RhoDomId,
					  long NuDomId,
					  long BetaDomId,
					  long RhoForId,
					  long NuForId,
					  long BetaForId,
					  long fxVolModelId,
					  ARM_result& result,
					  long objId)
{
	long modId;

	ARM_VolCurve* dVolatility = NULL;
	ARM_VolCurve* fVolatility = NULL;

	ARM_VolCurve* RhoDom = NULL;
	ARM_VolCurve* NuDom = NULL;
	ARM_VolCurve* BetaDom = NULL;

	ARM_VolCurve* RhoFor = NULL;
	ARM_VolCurve* NuFor = NULL;
	ARM_VolCurve* BetaFor = NULL;

	ARM_VolCurve* dFxCorrCurve = NULL;
	ARM_VolCurve* fFxCorrCurve = NULL;

	ARM_VolCurve* fxVolCurve = NULL;
	ARM_VolCurve* fxSmileCurve = NULL;

	ARM_ZeroCurve* dZeroCurve = NULL;
	ARM_ZeroCurve* fZeroCurve = NULL;

	ARM_ZeroCurve* dBSZeroCurve = NULL;
	ARM_ZeroCurve* fBSZeroCurve = NULL;

	ARM_ZeroCurve* FundZeroCurve = NULL;
	ARM_ZeroCurve* FundBSZeroCurve = NULL;

	ARM_DFBSModel* createdDFFXBSModel = NULL;
	ARM_DFBSModel* prevModel = NULL;

	ARM_PricingModel* fxVolModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		dZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic Zc Curve is not of a good type");
			return ARM_KO;
		}

		fZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fZcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fZeroCurve, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign Zc Curve is not of a good type");
			return ARM_KO;
		}

		dVolatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dVolatility, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic Volatility Curve is not of a good type");
			return ARM_KO;
		}

		fVolatility = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fVolatility, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign Volatility Curve is not of a good type");
			return ARM_KO;
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: fx Volatility Curve is not of a good type");
			return ARM_KO;
		}
		
		if ( fxSmileId != ARM_NULL_OBJECT)
		{
			fxSmileCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxSmileId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxSmileCurve, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: fx Smile Curve is not of a good type");

				return(ARM_KO);
			}
		}

		
		dFxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic fx Correlation Curve is not of a good type");
			return ARM_KO;

		}

		fFxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fFxCorrId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fFxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign fx Correlation Curve is not of a good type");
			return ARM_KO;
		}

		if (dBSZcId != ARM_NULL_OBJECT)
		{
			dBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dBSZcId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBSZeroCurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");
				return ARM_KO;
			}
		}

		if (fBSZcId != ARM_NULL_OBJECT)
		{
			fBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fBSZcId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fBSZeroCurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: foreign BS Zc Curve is not of a good type");
				return ARM_KO;
			}
		}

		if ( FundZcId != ARM_NULL_OBJECT)
		{
			FundZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FundZcId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FundZeroCurve, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");

				return(ARM_KO);
			}
		}

		if ( FundBSZcId != ARM_NULL_OBJECT)
		{
			FundBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FundBSZcId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FundBSZeroCurve, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");

				return(ARM_KO);
			}
		}

		if ( RhoDomId != ARM_NULL_OBJECT)
		{
			RhoDom = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(RhoDomId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(RhoDom, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho Dom is not of a good type");

				return(ARM_KO);
			}
		}

		if ( NuDomId != ARM_NULL_OBJECT)
		{
			NuDom = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(NuDomId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(NuDom, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu Dom is not of a good type");

				return(ARM_KO);
			}
		}

		if ( BetaDomId != ARM_NULL_OBJECT)
		{
			BetaDom = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(BetaDomId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(BetaDom, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Beta Dom is not of a good type");

				return(ARM_KO);
			}
		}

		if ( RhoForId != ARM_NULL_OBJECT)
		{
			RhoFor = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(RhoForId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(RhoFor, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho For is not of a good type");

				return(ARM_KO);
			}
		}

		if ( NuForId != ARM_NULL_OBJECT)
		{
			NuFor = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(NuForId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(NuFor, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu For is not of a good type");

				return(ARM_KO);
			}
		}

		if ( BetaForId != ARM_NULL_OBJECT)
		{
			BetaFor = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(BetaForId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(BetaFor, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Beta For is not of a good type");

				return(ARM_KO);
			}
		}

		if ( fxVolModelId != ARM_NULL_OBJECT)
		{
			fxVolModel = (ARM_PricingModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolModelId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolModel, ARM_PRICINGMODEL) == 0 )
			{
				result.setMsg ("ARM_ERR: fx vol model is not of a good type");

				return(ARM_KO);
			}
		}

		createdDFFXBSModel = new ARM_DFBSModel(dVolatility,
											   fVolatility,
											   dZeroCurve,
											   fZeroCurve,
											   dFxCorrCurve,
											   fFxCorrCurve,
											   fxVolCurve,
											   ratesCorr,
											   K_YIELD,
											   K_YIELD,
											   dBSZeroCurve,
											   fBSZeroCurve,
											   spot,
                                               FundZeroCurve,
                                               FundBSZeroCurve,
                                               fundFxSpot,
											   fxSmileCurve,
                                               isLnVols,
											   RhoDom,
											   NuDom,
											   BetaDom,
											   RhoFor,
											   NuFor,
											   BetaFor,
											   fxVolModel);

		if (createdDFFXBSModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFFXBSModel);

			if (modId == RET_KO)
			{
				if (createdDFFXBSModel)
					delete createdDFFXBSModel;
				createdDFFXBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_DFBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ARM_DFBSMODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFFXBSModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdDFFXBSModel)
					delete createdDFFXBSModel;
				createdDFFXBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdDFFXBSModel)
			delete createdDFFXBSModel;
		createdDFFXBSModel = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_DFFXBS(long dBSId,
					 long fBSId,
					 long dFxCorrId,
					 long fFxCorrId,
					 long fxVolId,
					 double ratesCorr,
					 long dBSZcId,
					 long fBSZcId,
					 double spot,
					 long FundZcId,
					 long FundBSZcId,
					 double fundFxSpot,
					 long fxSmileId,
					 long fxVolModelId,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_BSSmiledModel* dBS = NULL;
	ARM_BSSmiledModel* fBS = NULL;

	ARM_VolCurve* dFxCorrCurve = NULL;
	ARM_VolCurve* fFxCorrCurve = NULL;

	ARM_VolCurve* fxVolCurve = NULL;
	ARM_VolCurve* fxSmileCurve = NULL;

	ARM_ZeroCurve* dZeroCurve = NULL;
	ARM_ZeroCurve* fZeroCurve = NULL;

	ARM_ZeroCurve* dBSZeroCurve = NULL;
	ARM_ZeroCurve* fBSZeroCurve = NULL;

	ARM_ZeroCurve* FundZeroCurve = NULL;
	ARM_ZeroCurve* FundBSZeroCurve = NULL;

	ARM_DFBSModel* createdDFFXBSModel = NULL;
	ARM_DFBSModel* prevModel = NULL;

	ARM_PricingModel* fxVolModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		dBS = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(dBSId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBS, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: domestic BS model is not of a good type");
			return ARM_KO;
		}

		fBS = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(fBSId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fBS, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: foreign BS model is not of a good type");
			return ARM_KO;
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: fx Volatility Curve is not of a good type");
			return ARM_KO;
		}
		
		if ( fxSmileId != ARM_NULL_OBJECT)
		{
			fxSmileCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxSmileId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxSmileCurve, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: fx Smile Curve is not of a good type");
				return(ARM_KO);
			}
		}

		dFxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dFxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dFxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic fx Correlation Curve is not of a good type");
			return ARM_KO;

		}

		fFxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fFxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fFxCorrCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign fx Correlation Curve is not of a good type");
			return ARM_KO;
		}

		if (dBSZcId != ARM_NULL_OBJECT)
		{
			dBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dBSZcId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(dBSZeroCurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");
				return ARM_KO;
			}
		}

		if (fBSZcId != ARM_NULL_OBJECT)
		{
			fBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fBSZcId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fBSZeroCurve, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: foreign BS Zc Curve is not of a good type");
				return ARM_KO;
			}
		}

		if ( FundZcId != ARM_NULL_OBJECT)
		{
			FundZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FundZcId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FundZeroCurve, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");
				return(ARM_KO);
			}
		}

		if ( FundBSZcId != ARM_NULL_OBJECT)
		{
			FundBSZeroCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FundBSZcId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FundBSZeroCurve, ARM_ZERO_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: domestic BS Zc Curve is not of a good type");
				return(ARM_KO);
			}
		}

		if ( fxVolModelId != ARM_NULL_OBJECT)
		{
			fxVolModel = (ARM_PricingModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolModelId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolModel, ARM_PRICINGMODEL) == 0 )
			{
				result.setMsg ("ARM_ERR: fx vol model is not of a good type");

				return(ARM_KO);
			}
		}

		createdDFFXBSModel = new ARM_DFBSModel(dBS,
											   fBS,
											   dFxCorrCurve,
											   fFxCorrCurve,
											   fxVolCurve,
											   ratesCorr,
											   dBSZeroCurve,
											   fBSZeroCurve,
											   spot,
                                               FundZeroCurve,
                                               FundBSZeroCurve,
                                               fundFxSpot,
											   fxSmileCurve,
											   fxVolModel);

		if (createdDFFXBSModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFFXBSModel);

			if (modId == RET_KO)
			{
				if (createdDFFXBSModel)
					delete createdDFFXBSModel;
				createdDFFXBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_DFBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ARM_DFBSMODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDFFXBSModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdDFFXBSModel)
					delete createdDFFXBSModel;
				createdDFFXBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdDFFXBSModel)
			delete createdDFFXBSModel;
		createdDFFXBSModel = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_GetComputedSigmaSABR(long modId,
								   ARM_result& result)
{
	ARM_BSSmiledModel*  mod=NULL;
	double sigmaSABR;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		mod = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_BSSMILEDMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		sigmaSABR = mod->GetSABRSigma();

		result.setDouble(sigmaSABR);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_CALIBRATIONHWSV(double asof,
							  long zcId,
							  long volId,
							  long secId,
							  double amin,
							  double amax,
							  double volmin,
							  double volmax,
							  const VECTOR<double>& dates,
							  long pfId,
							  long RhoSwoptId,
							  long NuSwoptId,
							  long BetaSwoptId,
							  ARM_result& result,
							  long objId)
{
	long calibId;

	ARM_CalibrationHWSV* newCalib = NULL;
	ARM_CalibrationHWSV* oldCalib = NULL;

	ARM_ZeroCurve* zc = NULL;
	ARM_VolCurve* vol = NULL;
	ARM_Security* sec = NULL;
	ARM_Portfolio* pf = NULL;

	ARM_VolCurve* RhoSwopt = NULL;
	ARM_VolCurve* NuSwopt  = NULL;
	ARM_VolCurve* BetaSwopt= NULL;

	ARM_Vector* vDates = NULL;

	char sDate[11];

	CCString msg (" ");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR:Zc Curve is not of a good type");
			return ARM_KO;
		}

		vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: security is not of a good type");
			return ARM_KO;
		}

		if (pfId != ARM_NULL_OBJECT)
		{
			pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
			{
				result.setMsg ("ARM_ERR: portfolio is not of a good type");
				return ARM_KO;
			}
		}

		if (RhoSwoptId != ARM_NULL_OBJECT)
		{
			RhoSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(RhoSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(RhoSwopt, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: RhoSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (NuSwoptId != ARM_NULL_OBJECT)
		{
			NuSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(NuSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(NuSwopt, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: NuSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (BetaSwoptId != ARM_NULL_OBJECT)
		{
			BetaSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(BetaSwoptId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(BetaSwopt, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: BetaSwopt is not of a good type");
				return ARM_KO;
			}
		}

		if (dates.size () > 0)
			vDates = new ARM_Vector(dates.size ());

		for(int i = 0; i < dates.size (); i++)
		{
			Local_XLDATE2ARMDATE(dates[i],sDate);		
			ARM_Date tmpDate = (ARM_Date) sDate;
			vDates->Elt(i) = tmpDate.GetJulian();
		}

		Local_XLDATE2ARMDATE(asof,sDate);

		newCalib = new ARM_CalibrationHWSV((ARM_Date) sDate,
											zc,
											vol,
											sec,
											amin,
											amax,
											volmin,
											volmax,
											vDates,
											pf,
											RhoSwopt,
											NuSwopt,
											BetaSwopt);

		if (newCalib == NULL)
		{
			result.setMsg ("ARM_ERR: Calibration is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			calibId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCalib);

			if (calibId == RET_KO)
			{
				if (newCalib)
					delete newCalib;
				newCalib = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(calibId);

			return ARM_OK;
		}
		else
		{
			oldCalib = (ARM_CalibrationHWSV *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldCalib, ARM_HWSVCALIBRATION) == 1)
			{
				if (oldCalib)
				{
					delete oldCalib;
					oldCalib = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCalib, objId);

				return ARM_OK;
			}
			else
			{
				if (newCalib)
					delete newCalib;
				newCalib = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newCalib)
			delete newCalib;
		newCalib = NULL;

		if (vDates)
			delete vDates;
		vDates = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_CALIBRATE (long calibId,
						 long calibType,
						 ARM_result& result,
						 long objId)
{
	long modId;

	ARM_GYCSigVarModel* newModel = NULL;
	ARM_GYCSigVarModel* oldModel = NULL;

	ARM_CalibrationHWSV* calib = NULL;

	ARM_CalibrationHWSV_OUTPUT Output;

	char msg2[20];

	CCString msg (" ");

	try
	{
		calib = (ARM_CalibrationHWSV *) LOCAL_PERSISTENT_OBJECTS->GetObject(calibId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(calib, ARM_HWSVCALIBRATION) == 0)
		{
			result.setMsg ("ARM_ERR: calibrator is not of a good type");
			return ARM_KO;
		}

		newModel = calib->Calibrate(msg2,(ARM_Calibration_TYPE) calibType,Output);

		if (newModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
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
			oldModel = (ARM_GYCSigVarModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldModel, ARM_GYCSIGVARMODEL) == 1)
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
		x.DebugPrint();

		if (newModel)
			delete newModel;
		newModel = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_HWSIGVARFROMANA( long AnaModId,
							   long endDate,
							   long nbSteps,
							   ARM_result& result,
							   long objId)
{
	long modId;

	ARM_HWSigVarTree* HWTmod=NULL;
	ARM_HWSigVarTree* mod=NULL;
	ARM_GYCSigVarModel* anaMod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sEndDate = new char[11];
	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(endDate,sEndDate);

		anaMod = (ARM_GYCSigVarModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(AnaModId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(anaMod, ARM_GYCSIGVARMODEL) == 0)
		{
			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			result.setMsg ("ARM_ERR: HW Ana Mod is not of a good type");
			return ARM_KO;
		}

		HWTmod = new ARM_HWSigVarTree(anaMod,
                                      (ARM_Date) sEndDate,
                                      nbSteps);

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		if (HWTmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod);

			if (modId == RET_KO)
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			mod = (ARM_HWSigVarTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mod, ARM_HWSIGVARTREE) == 1)
			{
				if (mod)
				{
					delete mod;
					mod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)HWTmod, objId);

				return ARM_OK;
			}
			else
			{
				if (HWTmod)
					delete HWTmod;
				HWTmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (HWTmod)
			delete HWTmod;
		HWTmod = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_BSSMILEDCALIBRATE(long modelId,
								long pfId,
								long calVolOrNotId,
								long calRhoOrNotId,
								long calNuOrNotId,
                                long calBetaOrNotId,
								double volTenor,
								double timeStep,
								double minSig,
								double maxSig,
								double minRho,
								double maxRho,
								double minNu,
								double maxNu,
                                double minBeta,
                                double maxBeta,
								long interpMethodId,
								double tol,
								long maxIter,
								long gradCalcId,
								double lambda,
								long globOrBootstrapId,
								double beginSmoothMatu,
								ARM_result& result,
								long objId)
{
	long modId;

	ARM_BSSmiledModel* createdModel = NULL;
	ARM_BSSmiledModel* oldModel     = NULL;
	ARM_BSSmiledModel* model        = NULL;
	ARM_Portfolio* pf               = NULL;

	CCString msg ("");

	try
	{
		model = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(model, ARM_BSSMILEDMODEL) == 0 )
            &&
            ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(model, ARM_BSMODEL) == 0 )
           )
		{
		   result.setMsg("ARM_ERR: model is not of a good type");

		   return(ARM_KO);
		}

        if ( model->GetName() == ARM_BSMODEL )
        {
           ARM_BSSmiledModel* sabrModel = ((ARM_BSModel *) model)->GetSabrModel();

           model = sabrModel;
        }

		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0 )
		{
		   result.setMsg("ARM_ERR: portfolio is not of a good type");
			
           return(ARM_KO);
		}

		createdModel = Calibrate(model,
								 pf,
								 calVolOrNotId,
								 calRhoOrNotId,
								 calNuOrNotId,
								 volTenor,
								 timeStep,
								 minSig,
								 maxSig,
								 minRho,
								 maxRho,
								 minNu,
								 maxNu,
								 interpMethodId,
								 tol,
								 maxIter,
								 gradCalcId,
								 lambda,
								 globOrBootstrapId,
								 beginSmoothMatu,
                                 calBetaOrNotId,
                                 minBeta,
                                 maxBeta);

		if (createdModel == NULL)
		{
			result.setMsg ("ARM_ERR: calibration failed");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel);

			if (modId == RET_KO)
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldModel = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldModel, ARM_BSSMILEDMODEL) == 1)
			{
				if (oldModel)
				{
					delete oldModel;
					oldModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdModel, objId);

				return ARM_OK;
			}
			else
			{
				if (createdModel)
					delete createdModel;
				createdModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createdModel)
			delete createdModel;
		createdModel = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GETCALIBRATED_SIGRHONU(long modelId,
									 long paramId,
									 ARM_result& result,
									 long objId)
{
	long volId;

	ARM_VolCurve* newVol = NULL;
	ARM_VolCurve* oldVol = NULL;
	ARM_BSSmiledModel* model = NULL;

	CCString msg ("");

	try
	{
		model = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(model, ARM_BSSMILEDMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		switch (paramId)
		{
		case ARM_SIGMA :
			newVol = (ARM_VolCurve *)(model->GetVolatility())->Clone();
			break;

		case ARM_RHO :
			newVol = (ARM_VolCurve *)(model->GetRho())->Clone();
			break;

		case ARM_NU :
			newVol = (ARM_VolCurve *)(model->GetNu())->Clone();
			break;

		case ARM_BETA :
			newVol = (ARM_VolCurve *)(model->GetBeta())->Clone();
			break;
		}

		if (newVol == NULL)
		{
			result.setMsg ("ARM_ERR: VolCurve is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			volId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVol);

			if (volId == RET_KO)
			{
				if (newVol)
					delete newVol;
				newVol = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(volId);

			return ARM_OK;
		}
		else
		{
			oldVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(oldVol, ARM_VOL_CURVE) == 1)
			{
				if (oldVol)
				{
					delete oldVol;
					oldVol = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newVol, objId);

				return ARM_OK;
			}
			else
			{
				if (newVol)
					delete newVol;
				newVol = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newVol)
			delete newVol;
		newVol = NULL;

		ARM_RESULT();
	}
}



//______________________Ezzine's FRM Calibration________________________________________________________


long ARMLOCAL_CalibrationFRMModel(long zcId,
								  long pf1Id,
                                  long pf2Id,
								  long nbFactors,
								  long liborTypeId,
								  const VECTOR<double>& powers,
								  const VECTOR<double>& smoothParams,
								  const VECTOR<double>& ACalibrationSchedule,
								  const VECTOR<double>& KCalibrationSchedule,
								  long initA,
								  long initK,
								  const VECTOR<double>& meanRevA,
								  const VECTOR<double>& meanRevK,
                                  const VECTOR<double>& SchedPrice,
                                  const VECTOR<double>& mdec,
								  const VECTOR<double>& Bounds,
								  const VECTOR<double>& optimizerParams,
								  long calibrate,
                                  ARM_result& result,
								  long objId)
{
	long modId;

	ARM_FRMModel* FRMModel = NULL;
	ARM_FRMModel* oldFRMModel = NULL;

	ARM_ZeroCurve* zc = NULL;
	ARM_Vector* MCurve = NULL;// Attention, il faut refaire
	ARM_Portfolio* pf1 = NULL;
    ARM_Portfolio* pf2 = NULL;
 
	ARM_ReferenceValue* vrefA = NULL;
	ARM_ReferenceValue* vrefK = NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Vector* vVectA = NULL;    

	ARM_Vector* vVectK = NULL;

    ARM_Vector* vVectSchedPrice = NULL;

	ARM_Vector* vMeanRA = NULL;

	ARM_Vector* vMeanRK = NULL;

	ARM_Vector* vmdec = NULL;

	double * pdBounds = NULL;
	ARM_Matrix* vBounds = NULL;

	ARM_Vector* vPowers = NULL;

	ARM_Vector* vSmoothParams = NULL;

	ARM_Vector* vOptimizerParams = NULL;


	CCString msg ("");

    try
    {
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		pf1 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf1Id);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf1, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Pf Curve is not of a good type");
			return ARM_KO;
		}

        if (pf2Id != ARM_NULL_OBJECT)
        {
            pf2 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf2Id);

		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf2, ARM_PORTFOLIO) == 0)
		    {
			    result.setMsg ("ARM_ERR: Pf Curve is not of a good type");
			    return ARM_KO;
		    }
        }
        
		vrefA = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(initA);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vrefA, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: ACurve refValue is not of a good type");
			return ARM_KO;
		}

		vrefK = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(initK);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vrefK, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: KCurve refValue is not of a good type");
			return ARM_KO;
		}

		int i;       

            //_____________________________AScheduleCalibration______________________________________
		// Detruit plus bas
		vVectSchedPrice = CreateARMVectorFromVECTOR(SchedPrice);

		//_____________________________AScheduleCalibration______________________________________
		// Detruit plus bas
		vVectA = CreateARMVectorFromVECTOR(ACalibrationSchedule);
		
		//_____________________________KScheduleCalibration______________________________________
		// Detruit plus bas
		vVectK = CreateARMVectorFromVECTOR(KCalibrationSchedule);
		
		//_____________________________MeanReversionA____________________________________________
		// Affecter sans clone dans le constructeur , detruit dans le destructeur
		vMeanRA = CreateARMVectorFromVECTOR(meanRevA);

		//_____________________________MeanReversionK____________________________________________
		// Affecter sans clone dans le constructeur , detruit dans le destructeur
		vMeanRK = CreateARMVectorFromVECTOR(meanRevK);

        //_____________________________Mdec______________________________________________________
		// Detruit dans le constructeur
		vmdec = CreateARMVectorFromVECTOR(mdec);

		//_____________________________Bounds____________________________________________________
		// Detruit plus bas
        if(Bounds.size()>0)
        {
		    pdBounds = new double[Bounds.size()];
		    for(i = 0; i < Bounds.size (); i++)
			    pdBounds[i] = Bounds [i];

		    vBounds  = new ARM_Matrix(Bounds.size () / 2, 2, pdBounds);
        }
		
		//_____________________________Powers___________________________________________________
		// Detruit plus bas
		vPowers = CreateARMVectorFromVECTOR(powers);

		//_____________________________smoothParams_____________________________________________
		// Detruit plus bas
		vSmoothParams = CreateARMVectorFromVECTOR(smoothParams);

		//_____________________________optimizerParams___________________________________________
		// Detruit plus bas
		vOptimizerParams = CreateARMVectorFromVECTOR(optimizerParams);
		
		//_______________________Objects Initialization_________________________________________
		ARM_Date AsOfDate = zc->GetAsOfDate();
   
        char sDate[11];

        if(vVectSchedPrice)
        {
            int size = vVectSchedPrice->GetSize();
        
            for(i = 0; i < size; i++)
            {
                Local_XLDATE2ARMDATE(vVectSchedPrice->Elt(i),sDate);
                ARM_Date tmpDate = (ARM_Date) sDate;
                vVectSchedPrice->Elt(i)   =  tmpDate.GetJulian();   
            }
        }

        if(vVectA)
        {
            int sizeA = vVectA->GetSize();
        
            for(i = 0; i < sizeA; i++)
            {
                Local_XLDATE2ARMDATE(vVectA->Elt(i),sDate);
                ARM_Date tmpDate = (ARM_Date) sDate;
                vVectA->Elt(i)   =  tmpDate.GetJulian();   
            }
        }

        if(vVectK)
        {
            int sizeK = vVectK->GetSize();
        
            for(i = 0; i < sizeK; i++)
            {
                Local_XLDATE2ARMDATE(vVectK->Elt(i),sDate);
                ARM_Date tmpDate = (ARM_Date) sDate;
                vVectK->Elt(i)   =  tmpDate.GetJulian();   
            }
        }
        
        ARM_IRIndex* IRIndex = ((ARM_CapFloor*)pf1->GetAsset(0))->GetSwapLeg()->GetIRIndex();

		ARM_FRMMarkovVol* FRMMarkovVol = new ARM_FRMMarkovVol(
				nbFactors,
				AsOfDate,
				vVectA,
				vVectK,
				(ARM_ReferenceValue*) vrefA->Clone(),
				(ARM_ReferenceValue*) vrefK->Clone(),
				vMeanRA,
				vMeanRK,
				vVectSchedPrice,
				vmdec,
				vBounds,
				vPowers,
				vSmoothParams,
				IRIndex);

	    FRMModel = new ARM_FRMModel(AsOfDate,zc,FRMMarkovVol,pf1,pf2);
        
        if(vOptimizerParams)
        {
            double tol    = vOptimizerParams->Elt(1);
		    long maxIter  = long (vOptimizerParams->Elt(0));
		    long gradCalc = long (vOptimizerParams->Elt(2));

            if(calibrate)
                FRMModel->GlobalCalibrate(tol,maxIter, gradCalc);
        }

        

		// Delete

		if (vVectA)
			delete vVectA;
		vVectA = NULL;

		if (vVectK)
			delete vVectK;
		vVectK = NULL;

		if (vBounds)
			delete vBounds;
		vBounds = NULL;

		if (vPowers)
			delete vPowers;
		vPowers = NULL;

		if (vSmoothParams)
			delete vSmoothParams;
		vSmoothParams = NULL;

		if (vOptimizerParams)
			delete vOptimizerParams;
		vOptimizerParams = NULL;

		//______________________________________________________________________________________

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(FRMModel);

			if (modId == RET_KO)
			{
				if (FRMModel)
					delete FRMModel;
				FRMModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldFRMModel = (ARM_FRMModel*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFRMModel, ARM_FRMMODEL) == 1)
			{
				if (oldFRMModel)
				{
					delete oldFRMModel;
					oldFRMModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FRMModel, objId);

				return ARM_OK;
			}
			else
			{
				if (FRMModel)
					delete FRMModel;
				FRMModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();

		if (vVectA)
			delete vVectA;
		vVectA = NULL;

		if (vVectK)
			delete vVectK;
		vVectK = NULL;

		if (vBounds)
			delete vBounds;
		vBounds = NULL;

		if (vPowers)
			delete vPowers;
		vPowers = NULL;

		if (vSmoothParams)
			delete vSmoothParams;
		vSmoothParams = NULL;

		if (vOptimizerParams)
			delete vOptimizerParams;
		vOptimizerParams = NULL;

		if (FRMModel)
			delete FRMModel;
		FRMModel = NULL;

		ARM_RESULT();
	}
}


// Ezzine's Version
/*________________________________Boostrap Calibration Method ________________________________________________*/
long ARMLOCAL_BootstCalibFRMModel(long zcId,
								  long pf1Id,
								  long pf2Id,
                                  long pf3Id,
                                  long marketmodelId,
								  long liborTypeId,
								  const VECTOR<double>& CalibParams,
								  const VECTOR<double>& initCurve,
								  const VECTOR<double>& mdec,
								  double meanRev,
                                  long   nbfactor,
								  long   nbrows,
								  long nbcolumns,
                                  const VECTOR<double>& correlmatrix,
                                  long volType,
                                  long calibrate, 
                                  long PreInitialise,  
                                  double Presicion,
                                  double VegaLevel,
								  ARM_result& result,
								  long objId)
{
	long modId;
    int i;
    char sStartDate[11];

	ARM_FRMModel* FRMModel    = NULL;
	ARM_FRMModel* oldFRMModel = NULL;

	ARM_ZeroCurve* zc  = NULL;
	ARM_Portfolio* pf1  = NULL;
	ARM_Portfolio* pf2  = NULL;
    ARM_Portfolio* pf3  = NULL;
    ARM_Model* model  = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}  

    double * pdcorrelmatrix = NULL;
	ARM_Matrix* vcorrelmatrix = NULL;

	CCString msg ("");

    try
    {
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

        if (marketmodelId != ARM_NULL_OBJECT)
        {
            model = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(marketmodelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(model, ARM_MODEL) == 0)
		    {
			    result.setMsg ("ARM_ERR: Market Model is not of a good type");
			    return ARM_KO;
		    }
        }

		pf1 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf1, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Pf1 Curve is not of a good type");
			return ARM_KO;
		}
        if (pf2Id != ARM_NULL_OBJECT)
        {
		    pf2 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf2Id);
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf2, ARM_PORTFOLIO) == 0)
		    {
			    result.setMsg ("ARM_ERR: Pf2 Curve is not of a good type");
			    return ARM_KO;
		    }
        }
        
        if (pf3Id != ARM_NULL_OBJECT)
        {
		    pf3 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf3Id);
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf3, ARM_PORTFOLIO) == 0)
		    {
			    result.setMsg ("ARM_ERR: Pf2 Curve is not of a good type");
			    return ARM_KO;
		    }
        }

        //_____________________________Calib Params_____________________________________________		
		// Affecte sans clone dans le constructeur de ARM_FRMModel, Detruit dans le destructeur
		int SizeCalibParams = CalibParams.size();
        double MFlag = CalibParams[0];
        ARM_Vector* vcalibParams = NULL;              
        if(SizeCalibParams > 1)
        {
		    vcalibParams = new ARM_Vector(SizeCalibParams);
		    for(i = 1; i < SizeCalibParams; i++)
			    vcalibParams->Elt(i-1) = CalibParams[i];
        }        

		//_____________________________Initial Curve_____________________________________________		
		// Detruit dans le constructeur
		ARM_Vector* vInitCurve = CreateARMVectorFromVECTOR(initCurve);

		//_____________________________Initial Mdec_______________________________________________		
		// Detruit dans le constructeur
		ARM_Vector* vmdec = CreateARMVectorFromVECTOR(mdec);

        //_____________________________Correlation Matrix________________________________________

		// Detruit plus bas
        if(nbfactor > 1)
        {
		    pdcorrelmatrix = new double[correlmatrix.size()];
		    for(i = 0; i < correlmatrix.size (); i++)
			    pdcorrelmatrix[i] = correlmatrix [i];
		    vcorrelmatrix  = new ARM_Matrix(nbrows, nbcolumns, pdcorrelmatrix);

            delete [] pdcorrelmatrix;

            for(i = 0; i < nbrows; i++)
            {
                double startDate = vcorrelmatrix->Elt(i,0);
			    Local_XLDATE2ARMDATE(startDate,sStartDate);
                vcorrelmatrix->Elt(i,0) = ((ARM_Date)sStartDate).GetJulian();
            }
        }
				
		//_______________________Objects Initialization_________________________________________		
		ARM_FRMHWVol* FRMWHVol = new ARM_FRMHWVol(zc->GetAsOfDate(),
                                                  (ARM_INDEX_TYPE)liborTypeId,
                                                  pf1,
                                                  pf3,
												  vInitCurve,
												  vmdec,
                                                  MFlag,
												  meanRev,
                                                  nbfactor,
                                                  vcorrelmatrix,
                                                  volType);
			
		FRMModel = new ARM_FRMModel(zc->GetAsOfDate(),
                                    zc,FRMWHVol,
                                    pf1,
                                    pf2,
                                    pf3,
                                    vcalibParams,
                                    MFlag,
                                    model,
                                    PreInitialise,
                                    Presicion,VegaLevel);

        if(vcorrelmatrix)                                  
            delete vcorrelmatrix;        
		
        if(calibrate)
            FRMModel->Calibrate();        
		//______________________________________________________________________________________

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT(); 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(FRMModel);
			if (modId == RET_KO)
			{
				if (FRMModel)
					delete FRMModel;
				FRMModel = NULL;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);
			return ARM_OK;
		}
		else
		{
			oldFRMModel = (ARM_FRMModel*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFRMModel, ARM_FRMMODEL) == 1)
			{
				if (oldFRMModel)
				{
					delete oldFRMModel;
					oldFRMModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FRMModel, objId);
				return ARM_OK;
			}
			else
			{
				if (FRMModel)
					delete FRMModel;
				FRMModel = NULL;
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		if (FRMModel)
			delete FRMModel;
		FRMModel = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_GETFROMFRMMODEL(long modelId,
							  const CCString& param,
							  ARM_result& result)
{
	ARM_FRMModel* FRMModel = NULL;
	ARM_ReferenceValue* myRefval = NULL;
	ARM_ReferenceValue* ACurveToPrice = NULL;
	ARM_ReferenceValue* KCurveToPrice = NULL;
	ARM_Vector* myVector = NULL;
	ARM_Vector* AMean = NULL;
	ARM_Vector* KMean = NULL;

    double myValue;

	CCString msg (" ");

	try
	{
		FRMModel = (ARM_FRMModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FRMModel, ARM_FRMMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: FRM model is not of a good type");
			return ARM_KO;
		}

		ARM_FRMVol* FRMVol = FRMModel->GetFRMVol();
		long NbF = FRMVol->GetNbFactors();

		if (!strcmp((const char*)param, "ACURVE"))
		{
			if(FRMVol->GetName() == ARM_FRMHWVOL )
			{
				myRefval = FRMVol->GetCurrentCurve();
			}
			else if(FRMVol->GetName() == ARM_FRMMARKOVVOL)
			{
				FRMVol->GetCurvesPrice(ACurveToPrice, KCurveToPrice);
				myRefval = ACurveToPrice;
			}
		}

        else if (!strcmp((const char*)param, "SHIFTCURVE"))
		{
            myRefval = FRMVol->GetMCurve();
        }

        else if (!strcmp((const char*)param, "MEANREV"))
		{
			if(FRMVol->GetName() == ARM_FRMHWVOL )
			{
				myValue = FRMVol->GetMeanRevParams();
			}
			else if(FRMVol->GetName() == ARM_FRMMARKOVVOL)
			{
			    result.setMsg ("ARM_ERR: Unimplemented method");	
			}
		}
        else if (!strcmp((const char*)param, "PROBA_LOG"))
		{
			if(FRMVol->GetName() == ARM_FRMHWVOL )
			{
				myRefval = FRMVol->GetLogProbaCurve();
			}
			else if(FRMVol->GetName() == ARM_FRMMARKOVVOL)
			{
			    result.setMsg ("ARM_ERR: Unimplemented method");	
			}
		}
        

		else if (!strcmp((const char*)param,"KCURVE"))
		{
			FRMVol->GetCurvesPrice(ACurveToPrice, KCurveToPrice);
			myRefval = KCurveToPrice;
		}
		else if (!strcmp((const char*)param,"AMEAN"))
		{
			((ARM_FRMMarkovVol*)FRMVol)->GetMeanRevParams(AMean, KMean);
			myVector = AMean;
		}

		else if (!strcmp((const char*)param,"KMEAN"))
		{
			((ARM_FRMMarkovVol*)FRMVol)->GetMeanRevParams(AMean, KMean);
			myVector = KMean;
		}

		else if (!strcmp((const char*)param,"ERROR"))
		{
			myVector = FRMModel->GetResult();
		}

		else
		{
			return ARM_KO;
		}

		if ( (myRefval == NULL) && (myVector == NULL)
             && strcmp((const char*)param, "MEANREV"))
		{
			result.setMsg ("returned Object is NULL");
			return ARM_KO;
		}

		if (myRefval != NULL)
		{
			for (int i=0;i<myRefval->GetSize();i++)
			{
				result.setArray(myRefval->GetDiscreteDates()->Elt(i),i);
				result.setArray(myRefval->GetDiscreteValues(0)->Elt(i),i+myRefval->GetSize());
				
				if (myRefval->NbCurves() == 2)
					result.setArray(myRefval->GetDiscreteValues(1)->Elt(i),i+2*myRefval->GetSize());
				else
					result.setArray(myRefval->GetDiscreteValues(0)->Elt(i),i+2*myRefval->GetSize());
			}

			result.setLong(myRefval->GetSize());
		}
		else if (myVector != NULL)
		{
			for (int i=0;i<myVector->GetSize();i++)
			{
				result.setArray(myVector->Elt(i),i);
			}

			result.setLong(myVector->GetSize());
		}
        else
        {
            result.setLong(1);
            result.setArray(myValue,0);
        }

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_BootstCalibFRMModelMixture(int N,
										long zcId,
										long pf1Id,
										long pf2Id,
										int nbProducts,
										long liborTypeId,
										const VECTOR<double>& lambdas,
										double meanRev,
										long   nbfactor,
										long   nbrows,
										long nbcolumns,
										const VECTOR<double>& correlmatrix,
										long VolType,
										const VECTOR<double>& lowerBound,
										const VECTOR<double>& upperBound,
										long calibrate,
										ARM_result& result,
										long objId)
{
	long modId;
    int i;
    char sStartDate[11];

	ARM_FRMModelMixture* FRMModelMixture    = NULL;
	ARM_FRMModelMixture* oldFRMModelMixture = NULL;

	ARM_ZeroCurve* zc  = NULL;
	ARM_Portfolio* pf1  = NULL;
	ARM_Portfolio* pf2  = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}  

    double * pdcorrelmatrix = NULL;
	ARM_Matrix* vcorrelmatrix = NULL;

	CCString msg ("");

    try
    {
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		pf1 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf1, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Pf is not of a good type");
			return ARM_KO;
		}

		if (pf2Id != ARM_NULL_OBJECT)
        {
		    pf2 = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pf2Id);
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf2, ARM_PORTFOLIO) == 0)
		    {
			    result.setMsg ("ARM_ERR: Pf2 Curve is not of a good type");
			    return ARM_KO;
		    }
        }

		//_____________________________Initial Curve _____________________________________________		
		ARM_Matrix* mInitCurve = new ARM_Matrix(N, nbProducts);

		//_____________________________Initial Mdec_______________________________________________		
		ARM_Vector* vmdec = new ARM_Vector(nbProducts);
		
		//_____________________________Initial Spread_______________________________________________
		ARM_Matrix* spreadInitCurve = new ARM_Matrix(N, nbProducts);
		
        //_____________________________Correlation Matrix________________________________________

		pdcorrelmatrix = new double[correlmatrix.size()];
		for(i = 0; i < correlmatrix.size (); i++)
			pdcorrelmatrix[i] = correlmatrix [i];
		vcorrelmatrix  = new ARM_Matrix(nbrows, nbcolumns, pdcorrelmatrix);

        for(i = 0; i < nbrows; i++)
        {
            double startDate = vcorrelmatrix->Elt(i,0);
			Local_XLDATE2ARMDATE(startDate,sStartDate);
            vcorrelmatrix->Elt(i,0) = ((ARM_Date)sStartDate).GetJulian();
        }

		//_____________________________Initial Lambdas_________________________________________
		int SizeLambdas = lambdas.size();
		ARM_Vector vlambdas(SizeLambdas);
		for(i = 0; i < SizeLambdas; i++)
		{
			vlambdas.Elt(i) =lambdas[i]; 
		}


		//_____________________________Initial Bound Low_________________________________________
		int SizeBoundLow = lowerBound.size();
		ARM_Vector vboundlow(SizeBoundLow);
		for(i = 0; i < SizeBoundLow; i++)
		{
			vboundlow.Elt(i) =lowerBound[i]; 
		}

		//_____________________________Initial Bound Up__________________________________________
		int SizeBoundUp = upperBound.size();
		ARM_Vector vboundup(SizeBoundUp);
		for(i = 0; i < SizeBoundUp; i++)
		{
			vboundup.Elt(i) = upperBound[i];
		}
		
		std::vector<ARM_FRMVol*> FRMWHVols(N);

		int nbProductsPerMaturity = pf1->GetSize()/nbProducts;

		//_______________________Objects Initialization_________________________________________		
		for (int i = 0; i < N; ++i)
		{
			ARM_Vector* vmdecCopy = (ARM_Vector*) vmdec->Clone();
			FRMWHVols[i] = new ARM_FRMHWVol(zc->GetAsOfDate(),
                                            (ARM_INDEX_TYPE)liborTypeId,
                                                  pf1,
												  nbProductsPerMaturity,
												  mInitCurve->GetRow(i),
												  vmdecCopy,
												  spreadInitCurve->GetRow(i),
												  meanRev,
                                                  nbfactor,
                                                  vcorrelmatrix,
                                                  VolType);
		}

		if (vmdec)
		{
			delete vmdec;
		}
		
		FRMModelMixture = new ARM_FRMModelMixture(
			N,
			zc->GetAsOfDate(),
			zc,
			FRMWHVols,
			vlambdas,
			pf1,
			pf2,
			nbProducts);

        if (vcorrelmatrix)
		{
			delete vcorrelmatrix;
			vcorrelmatrix = NULL;
		}

        if(calibrate)
            FRMModelMixture->BootstCalibrate(
			vboundlow,
			vboundup
			);

		//______________________________________________________________________________________

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT(); 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(FRMModelMixture);
			if (modId == RET_KO)
			{
				if (FRMModelMixture)
					delete FRMModelMixture;
				FRMModelMixture = NULL;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);
			return ARM_OK;
		}
		else
		{
			oldFRMModelMixture = (ARM_FRMModelMixture*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFRMModelMixture, ARM_FRMMODELMIXTURE) == 1)
			{
				if (oldFRMModelMixture)
				{
					delete oldFRMModelMixture;
					oldFRMModelMixture = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)FRMModelMixture, objId);
				return ARM_OK;
			}
			else
			{
				if (FRMModelMixture)
					delete FRMModelMixture;
				FRMModelMixture = NULL;
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}
	
	catch(Exception& x)
	{
		if (vcorrelmatrix)
			delete vcorrelmatrix;

		x.DebugPrint();
		if (FRMModelMixture)
			delete FRMModelMixture;
		FRMModelMixture = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_GETFROMFRMMODELMIXTURE(long modelId,
							  int modelNumber,
							  const CCString& param,
							  ARM_result& result)
{
	ARM_FRMModelMixture* FRMModelMixture = NULL;
	ARM_FRMModel* FRMModel = NULL;
	ARM_ReferenceValue* myRefval = NULL;
	ARM_ReferenceValue* ACurveToPrice = NULL;
	ARM_ReferenceValue* KCurveToPrice = NULL;
	ARM_Vector* myVector = NULL;
	ARM_Vector* AMean = NULL;
	ARM_Vector* KMean = NULL;

    double myValue;

	CCString msg (" ");

	try
	{
		FRMModelMixture = (ARM_FRMModelMixture *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FRMModelMixture, ARM_FRMMODELMIXTURE) == 0)
		{
			result.setMsg ("ARM_ERR: FRM Mixture model is not of a good type");
			return ARM_KO;
		}

		if (modelNumber >= FRMModelMixture->GetN())
		{
			result.setMsg ("ARM_ERR: Model number is too high");
			return ARM_KO;
		}

		FRMModel = FRMModelMixture->GetModel(modelNumber);

		ARM_FRMVol* FRMVol = FRMModel->GetFRMVol();
		long NbF = FRMVol->GetNbFactors();

		if (!strcmp((const char*)param, "ACURVE"))
		{
			if(FRMVol->GetName() == ARM_FRMHWVOL )
			{
				myRefval = FRMVol->GetCurrentCurve();
			}
			else if(FRMVol->GetName() == ARM_FRMMARKOVVOL)
			{
				FRMVol->GetCurvesPrice(ACurveToPrice, KCurveToPrice);
				myRefval = ACurveToPrice;
			}
		}
		else if (!strcmp((const char*)param, "SPREADCURVE"))
		{
            myRefval = FRMVol->GetSpreadCurve();
        }
        else if (!strcmp((const char*)param, "MEANREV"))
		{
			if(FRMVol->GetName() == ARM_FRMHWVOL )
			{
				myValue = FRMVol->GetMeanRevParams();
			}
			else if(FRMVol->GetName() == ARM_FRMMARKOVVOL)
			{
			    result.setMsg ("ARM_ERR: Unimplemented method");	
			}
		}
		else
		{
			return ARM_KO;
		}

		if ( (myRefval == NULL) && (myVector == NULL) && strcmp((const char*)param, "MEANREV") )
		{
			result.setMsg ("returned Object is NULL");
			return ARM_KO;
		}

		if (myRefval != NULL)
		{
			for (int i=0;i<myRefval->GetSize();i++)
			{
				result.setArray(myRefval->GetDiscreteDates()->Elt(i),i);
				result.setArray(myRefval->GetDiscreteValues(0)->Elt(i),i+myRefval->GetSize());
				
				if (myRefval->NbCurves() == 2)
					result.setArray(myRefval->GetDiscreteValues(1)->Elt(i),i+2*myRefval->GetSize());
				else
					result.setArray(myRefval->GetDiscreteValues(0)->Elt(i),i+2*myRefval->GetSize());
			}

			result.setLong(myRefval->GetSize());
		}
		else if (myVector != NULL)
		{
			for (int i=0;i<myVector->GetSize();i++)
			{
				result.setArray(myVector->Elt(i),i);
			}

			result.setLong(myVector->GetSize());
		}
        else
        {
            result.setLong(1);
            result.setArray(myValue,0);
        }

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_MCFRMMODEL(long FRMModId,
						 long nbTraj,
						 long nbStepIn,
                         long pricerTypeId,
                         long mcMethod,
						 ARM_result& result,
						 long objId)
{
	long modId;

	ARM_FRMMCModel* newFRMMCModel=NULL;
	ARM_FRMMCModel* oldFRMMCModel=NULL;
	ARM_FRMModel* FRMModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		FRMModel = (ARM_FRMModel *)(LOCAL_PERSISTENT_OBJECTS->GetObject(FRMModId));
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FRMModel, ARM_FRMMODEL) == 0)
		{			
			result.setMsg ("ARM_ERR: FRMModel  is not of a good type");
			return ARM_KO;
		}

        newFRMMCModel = new ARM_FRMMCModel(FRMModel,nbTraj, nbStepIn,(ARM_PRICER_TYPE)pricerTypeId,mcMethod);
		if (newFRMMCModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFRMMCModel);

			if (modId == RET_KO)
			{
				if (newFRMMCModel)
					delete newFRMMCModel;
				newFRMMCModel = NULL;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldFRMMCModel = dynamic_cast <ARM_FRMMCModel *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFRMMCModel, ARM_FRMMCMODEL) == 1)
			{
				if (oldFRMMCModel)
				{
					delete oldFRMMCModel;
					oldFRMMCModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFRMMCModel, objId);
				result.setLong(objId);

				return ARM_OK;
			}
			else
			{
				if (newFRMMCModel)
					delete newFRMMCModel;
				newFRMMCModel = NULL;
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newFRMMCModel)
			delete newFRMMCModel;
		newFRMMCModel = NULL;
		
		ARM_RESULT();
	}
}

long ARMLOCAL_MCFRMMODELMIXTURE(long FRMModMixId,
						 long nbTraj,
						 long nbStepIn,
                         long pricerTypeId,
						 ARM_result& result,
						 long objId)
{
	long modMixId;

	ARM_FRMMCModelMixture* newFRMMCModelMixture=NULL;
	ARM_FRMMCModelMixture* oldFRMMCModelMixture=NULL;
	ARM_FRMModelMixture* FRMModelMixture = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		FRMModelMixture = (ARM_FRMModelMixture *)(LOCAL_PERSISTENT_OBJECTS->GetObject(FRMModMixId));
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FRMModelMixture, ARM_FRMMODELMIXTURE) == 0)
		{			
			result.setMsg ("ARM_ERR: FRMModel  is not of a good type");
			return ARM_KO;
		}

        newFRMMCModelMixture = new ARM_FRMMCModelMixture(FRMModelMixture,nbTraj, nbStepIn,(ARM_PRICER_TYPE)pricerTypeId);
		if (newFRMMCModelMixture == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modMixId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFRMMCModelMixture);

			if (modMixId == RET_KO)
			{
				if (newFRMMCModelMixture)
					delete newFRMMCModelMixture;
				newFRMMCModelMixture = NULL;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modMixId);

			return ARM_OK;
		}
		else
		{
			oldFRMMCModelMixture = dynamic_cast <ARM_FRMMCModelMixture *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFRMMCModelMixture, ARM_FRMMCMODELMIXTURE) == 1)
			{
				if (oldFRMMCModelMixture)
				{
					delete oldFRMMCModelMixture;
					oldFRMMCModelMixture = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFRMMCModelMixture, objId);
				result.setLong(objId);

				return ARM_OK;
			}
			else
			{
				if (newFRMMCModelMixture)
					delete newFRMMCModelMixture;
				newFRMMCModelMixture = NULL;
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newFRMMCModelMixture)
			delete newFRMMCModelMixture;
		newFRMMCModelMixture = NULL;
		
		ARM_RESULT();
	}
}
 
long ARMLOCAL_FRMMARKOVTREE(double startDate,
                            double horizon,
                            long zcId,
                            long PathNumber,
                            const VECTOR<double>& Params,
                            long FRMModelId,
                            ARM_result& result,
                            long objId)
{
	long modId;

	ARM_FRMMarkovTree* newMod = NULL;
	ARM_FRMMarkovTree* oldMod = NULL;
	ARM_ZeroCurve* zc = NULL;
    ARM_FRMModel* FRMModel = NULL;    
    ARM_Vector* Parameters = NULL;

	char sStartDate[11];
	char sHorizon[11];

    CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: zc is not of a good type");
			return ARM_KO;
		}
        
        FRMModel = (ARM_FRMModel*) LOCAL_PERSISTENT_OBJECTS->GetObject(FRMModelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(FRMModel, ARM_FRMMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: FRMModel is not of a good type");
			return ARM_KO;
		}

		Parameters = CreateARMVectorFromVECTOR(Params);

        newMod = new ARM_FRMMarkovTree( (ARM_Date)sStartDate,
                                        (ARM_Date)sHorizon,
                                        zc,
                                        PathNumber,
                                        Parameters,
                                        FRMModel);

        if (Parameters)
			delete Parameters;
		Parameters = NULL;

        if (newMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod);

			if (modId == RET_KO)
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldMod = (ARM_FRMMarkovTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldMod, ARM_FRMMARKOVTREE) == 1)
			{
				if (oldMod)
				{
					delete oldMod;
					oldMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod, objId);

				return ARM_OK;
			}
			else
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newMod)
			delete newMod;
		newMod = NULL;

        if (Parameters)
			delete Parameters;
		Parameters = NULL;             

		ARM_RESULT();
	}
}



long ARMLOCAL_LOGDECANA(long zcId,
						long frequency,
						const VECTOR<double>& resetMatu,
						const VECTOR<double>& shifts,
						const VECTOR<double>& fwdVols,
						long isWeighted,
						ARM_result& result,
						long objId)
{
	long modId;

	ARM_ZeroCurve* zc = NULL;

    ARM_Vector* vResetMatu = NULL;
    ARM_Vector* vShifts = NULL;
    ARM_Matrix* mFwdVols = NULL;

	ARM_LogDecalANA* newMod=NULL;
	ARM_LogDecalANA* oldMod=NULL;

	if (resetMatu.size() != shifts.size())
	{
		result.setMsg ("ARM_ERR: reset Dates and shifts must have the same size");
		return ARM_KO;
	}

    CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: zc is not of a good type");
			return ARM_KO;
		}

   		vResetMatu = new ARM_Vector(resetMatu.size());
		for (int i=0;i<resetMatu.size();i++)
			vResetMatu->Elt(i) = resetMatu[i];

		vShifts = new ARM_Vector(shifts.size());
		for (i=0;i<shifts.size();i++)
			vShifts->Elt(i) = shifts[i];

		mFwdVols = new ARM_Matrix(shifts.size(),shifts.size());
		for (i=0;i<shifts.size();i++)
			for (int j=0;j<shifts.size();j++)
				mFwdVols->Elt(i,j) = fwdVols[i*shifts.size()+j];

		newMod = new ARM_LogDecalANA(zc,
									 frequency,
									 (ARM_Vector*)vResetMatu->Clone(),
									 (ARM_Vector*)vShifts->Clone(),
									 (ARM_Matrix*)mFwdVols->Clone(),
									 isWeighted);
;

		if (vResetMatu)
			delete vResetMatu;
		vResetMatu = NULL;

        if (vShifts)
			delete vShifts;
		vShifts = NULL;   

        if (mFwdVols)
			delete mFwdVols;
		mFwdVols = NULL;   

		if (newMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod);

			if (modId == RET_KO)
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldMod = dynamic_cast<ARM_LogDecalANA *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldMod, ARM_LOGDECALANA) == 1)
			{
				if (oldMod)
				{
					delete oldMod;
					oldMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod, objId);

				return ARM_OK;
			}
			else
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newMod)
			delete newMod;
		newMod = NULL;

		if (vResetMatu)
			delete vResetMatu;
		vResetMatu = NULL;

        if (vShifts)
			delete vShifts;
		vShifts = NULL;   

        if (mFwdVols)
			delete mFwdVols;
		mFwdVols = NULL;   

		ARM_RESULT();
	}
}



long ARMLOCAL_QMODEL(double date,
					 double spot,
					 long dividend_type,
					 double dividend,
					 long discrate_type,
					 double discrate,
					 long volat_type,
					 double volat,
					 double q0,
					 double q1,
					 long typstk,
					 ARM_result& result,
					 long objId)
{
	long modId;

	ARM_Q_Model* Qmod = NULL;
	ARM_Q_Model* newQmod = NULL;
	ARM_ZeroCurve* divid = NULL;
	ARM_ZeroCurve* drate = NULL;
	ARM_VolCurve*  vol = NULL;

	ARM_ZeroCurve* tmpdivid = NULL;
	ARM_ZeroCurve* tmpdrate = NULL;
	ARM_VolCurve*  tmpvol = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * startDate = new char[11];
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,startDate);

		if ( dividend_type == 1 )
		{
			divid = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)dividend);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(divid, ARM_ZERO_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				result.setMsg ("ARM_ERR: Dividend is not a Curve");
				return ARM_KO;
			}

		}
		else
		{
			divid = new ARM_ZeroFlat((ARM_Date) startDate, dividend);

			tmpdivid = divid;
		}


		if ( discrate_type == 1 )
		{
			drate = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)discrate);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(drate, ARM_ZERO_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (tmpdivid)
					delete tmpdivid;
				tmpdivid = NULL;

				result.setMsg ("ARM_ERR: discount rate is not a Curve");
				return ARM_KO;
			}
		}
		else
		{
			drate = new ARM_ZeroFlat((ARM_Date) startDate, discrate);

			tmpdrate = drate;
		}

		if (volat_type == 1)
		{
			vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)volat);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (tmpdivid)
					delete tmpdivid;
				tmpdivid = NULL;

				if (tmpdrate)
					delete tmpdrate;
				tmpdrate = NULL;

				result.setMsg ("ARM_ERR: volatility is not a Vol Curve");
				return ARM_KO;
			}

		}
		else
		{
		   vol = new ARM_VolFlat((ARM_Date) startDate, volat);

		   tmpvol = vol;
		}

		newQmod = new ARM_Q_Model((ARM_Date) startDate,
								   spot,
								   divid,
								   drate,
								   vol,
								   q0,
								   q1,
								   typstk);

		if (tmpdivid)
		   delete tmpdivid;

		if (tmpdrate)
		   delete tmpdrate;

		if (tmpvol)
		   delete tmpvol;

		if (startDate)
			delete [] startDate;

		if (newQmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newQmod);

			if (modId == RET_KO)
			{
				if (newQmod)
					delete newQmod;
				newQmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			Qmod = (ARM_Q_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Qmod, ARM_Q_MODEL) == 1)
			{
				if (Qmod)
				{
					delete Qmod;
					Qmod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newQmod, objId);

				return ARM_OK;
			}
			else
			{
				if (newQmod)
					delete newQmod;
				newQmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (tmpdivid)
			delete tmpdivid;
		tmpdivid = NULL;

		if (tmpdrate)
			delete tmpdrate;
		tmpdrate = NULL;

		if (tmpvol)
			delete tmpvol;
		tmpvol = NULL;

		if (newQmod)
			delete newQmod;
		newQmod = NULL;

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_CROSSMODEL(double date,
						long DbsmodId,
						long FbsmodId,
						long fxvolId,
						long dfxcorrelId,
						long ffxcorrelId,
						long dfcorrelId,
						long discountZcId,
						double correlforAdj,
						double adjustFlag,
						double slopeFlag,
						ARM_result& result,
						long objId)
{
	long modId;

	ARM_VolCurve* DfxCorr = NULL;
	ARM_VolCurve* FfxCorr = NULL;
	ARM_VolCurve* DomIdxFgnIdxCorrel = NULL;

	ARM_VolCurve* fxVolCurve = NULL;

	ARM_BSModel* DBSModel = NULL;
	ARM_BSModel* FBSModel = NULL;

	ARM_ZeroCurve* discountZc = NULL;

	ARM_CrossModel* newcrossmod = NULL;
	ARM_CrossModel* oldcrossmod = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	char * startDate = new char[11];

	try
	{
		Local_XLDATE2ARMDATE(date,startDate);

		DBSModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(DbsmodId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DBSModel, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: domestic BS Model is not of a good type");
			return ARM_KO;
		}

		FBSModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(FbsmodId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FBSModel, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: foreign BS Model is not of a good type");
			return ARM_KO;
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxvolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Forex Volatility Curve is not of a good type");
			return ARM_KO;
		}

		DfxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dfxcorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DfxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: domestic index/ forex Correlation Curve is not of a good type");
			return ARM_KO;
		}

		FfxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ffxcorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FfxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: foreign index/ forex Correlation Curve is not of a good type");
			return ARM_KO;
		}

		DomIdxFgnIdxCorrel = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(dfcorrelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DomIdxFgnIdxCorrel, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Domestic/foreign index Correlation Curve is not of a good type");
			return ARM_KO;
		}
	
		
		if (discountZcId == ARM_NULL_OBJECT)
		{
			newcrossmod = new ARM_CrossModel((ARM_Date) startDate,
											   DBSModel,
											   FBSModel,
											   fxVolCurve,
											   DfxCorr,
											   FfxCorr,
											   DomIdxFgnIdxCorrel,
											   NULL,
											   correlforAdj,
											   (int) adjustFlag,
											   (int) slopeFlag);

		}
		else
		{
			discountZc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(discountZcId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(discountZc, ARM_ZERO_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Discount Zero Coupon Curve is not of a good type");
				return ARM_KO;
			}
		
			newcrossmod = new ARM_CrossModel((ARM_Date) startDate,
											   DBSModel,
											   FBSModel,
											   fxVolCurve,
											   DfxCorr,
											   FfxCorr,
											   DomIdxFgnIdxCorrel,
											   discountZc,
											   correlforAdj,
											   (int) adjustFlag,
											   (int) slopeFlag);
		}
		if (startDate)
			delete [] startDate;

		if (newcrossmod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}
		
		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcrossmod);

			if (modId == RET_KO)
			{
				if (newcrossmod)
					delete newcrossmod;
				newcrossmod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldcrossmod = (ARM_CrossModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldcrossmod, ARM_CROSSMODEL) == 1)
			{
				if (oldcrossmod)
				{
					delete oldcrossmod;
					oldcrossmod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcrossmod, objId);

				return ARM_OK;
			}
			else
			{
				if (newcrossmod)
					delete newcrossmod;
				newcrossmod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newcrossmod)
			delete newcrossmod;
		newcrossmod = NULL;

		if (oldcrossmod)
			delete oldcrossmod;
		oldcrossmod = NULL;

		ARM_RESULT();
	}




}
long ARMLOCAL_GLOBDFBS(long DomBSId,
					   long DomCurrId,
					   long FrgBSId,
					   long FrgCurrId,
					   long fxVolCrvId,
					   long FFxCorrId,
					   long RatesCorrId,
					   long FxVolModelId,
					   ARM_result& result,
					   long objId)
{
	long modId;

	ARM_BSModel* dbs = NULL;
	ARM_BSModel* fbs = NULL;
	ARM_Currency* dCurr = NULL;
	ARM_Currency* fCurr = NULL;
	ARM_VolCurve* fxVolCurve = NULL;
	ARM_VolCurve* ffxCorrCurve = NULL;
	ARM_VolCurve* RatesCorrCurve = NULL;

	ARM_DFBSModel* createdMod = NULL;
	ARM_DFBSModel* oldMod = NULL;

	ARM_PricingModel* fxVolModel = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	try
	{
		dbs = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(DomBSId);

		if (dbs == NULL)
		{
		   result.setMsg ("ARM_ERR: Domestic model is NULL");
		
           return(ARM_KO);
		}

		if (!dbs->IsBSLikeModel())
		{
		   result.setMsg ("ARM_ERR: Domestic model is not a BSModel like");
		
           return(ARM_KO);
		}

		fbs = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(FrgBSId);

		if (fbs == NULL)
		{
		   result.setMsg ("ARM_ERR: Foreign model is NULL");
		
           return(ARM_KO);
		}

		if (!fbs->IsBSLikeModel())
		{
		   result.setMsg ("ARM_ERR: Foreign model is not a BSModel like");
		
           return(ARM_KO);
		}

		dCurr = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(DomCurrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(dCurr, ARM_CURRENCY) == 0 )
		{
		   result.setMsg ("ARM_ERR: Domestic Ccy is not a Currency");
			
           return(ARM_KO);
		}

		fCurr = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(FrgCurrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(fCurr, ARM_CURRENCY) == 0 )
		{
		   result.setMsg ("ARM_ERR: Foreign Ccy is not a Currency");
		
           return(ARM_KO);
		}

		fxVolCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxVolCrvId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolCurve, ARM_VOL_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: FX volcurve is not a VolCurve");

		   return(ARM_KO);
		}

		ffxCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(FFxCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ffxCorrCurve, ARM_VOL_CURVE) == 0 )
		{
		   result.setMsg("ARM_ERR: ffxCorrCurve is not a VolCurve");
			
           return(ARM_KO);
		}

		RatesCorrCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(RatesCorrId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(RatesCorrCurve, ARM_VOL_CURVE) == 0 )
		{
		   result.setMsg ("ARM_ERR: RatesCorrCurve is not a VolCurve");

		   return(ARM_KO);
		}

		if (FxVolModelId != ARM_NULL_OBJECT)
		{
			fxVolModel = (ARM_PricingModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(FxVolModelId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fxVolModel, ARM_PRICINGMODEL) == 0 )
			{
			   result.setMsg ("ARM_ERR: Fx vol model is not a pricing model");

			   return(ARM_KO);
			}
		}

		createdMod = new ARM_DFBSModel(dbs,
									   dCurr,
									   fbs,
									   fCurr,
									   fxVolCurve,
									   ffxCorrCurve,
									   RatesCorrCurve,
									   fxVolModel);

		if (createdMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMod);

			if (modId == RET_KO)
			{
				if (createdMod)
					delete createdMod;
				createdMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldMod = (ARM_DFBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldMod, ARM_DFBSMODEL) == 1)
			{
				if (oldMod)
				{
					delete oldMod;
					oldMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMod, objId);

				return ARM_OK;
			}
			else
			{
				if (createdMod)
					delete createdMod;
				createdMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createdMod)
			delete createdMod;
		createdMod = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_BSConvAdjust_Create(int SUMMITFormulaeUsed, int UseSabrCMS, ARM_result& result, long objId)
{
	long convAdjustId;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	ARM_BSConvAdjust* createConvAdjutst = NULL;

	try
	{
		createConvAdjutst = new ARM_BSConvAdjust(SUMMITFormulaeUsed, UseSabrCMS);

		if (createConvAdjutst == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			convAdjustId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst);

			if (convAdjustId == RET_KO)
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(convAdjustId);

			return ARM_OK;
		}
		else
		{
			ARM_BSConvAdjust* oldConvAdjutst = (ARM_BSConvAdjust *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldConvAdjutst, ARM_BSCONVADJUST) == 1)
			{
				if (oldConvAdjutst)
				{
					delete oldConvAdjutst;
					oldConvAdjutst = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst, objId);

				return ARM_OK;
			}
			else
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createConvAdjutst)
			delete createConvAdjutst;
		createConvAdjutst = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_BSConvAdjustRep_Create(long UsedModelId, long swoptVolCurveId, 
											const VECTOR<double>& Stddev, 
											int NbPtsForRepliq,
											long MRId,
											bool FullRepliq,
											double upperProba, 
											double lowerProba,
											ARM_result& result,
											long objId)
{
	long convAdjustId;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	ARM_BSConvAdjustRep* createConvAdjutst = NULL;

	try
	{
		ARM_Model * UsedModel;

		if (UsedModelId != ARM_NULL_OBJECT)
        {
            UsedModel = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(UsedModelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(UsedModel, ARM_MODEL) == 0)
		    {			
			    result.setMsg ("ARM_ERR: Used  is not of a good type");
			    return ARM_KO;
		    }
        }
		else
		{
			result.setMsg("ARM_ERR : Used Model is NULL");
			return ARM_KO;
		}
		
		ARM_VolLInterpol * Volatility;

		if(swoptVolCurveId != ARM_NULL_OBJECT)
		{
			Volatility = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)swoptVolCurveId);
			if(!Volatility)
			{
				result.setMsg ("ARM_ERR: Swaption Volatility is not a Volatility Curve");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg("ARM_ERR : Swopt Volatility is NULL");
			return ARM_KO;
		}

		ARM_VolLInterpol * MR;

		if(MRId != ARM_NULL_OBJECT)
		{
			MR = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)MRId);
			if(!MR)
			{
				result.setMsg ("ARM_ERR: Mean Reversion Surface is not a Volatility Curve");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg("ARM_ERR : Mean Reversion Surface is NULL");
			return ARM_KO;
		}


		ARM_Vector stddev(Stddev.size());
		for(int k = 0; k < Stddev.size(); k++) stddev[k] = Stddev[k];

		createConvAdjutst = new ARM_BSConvAdjustRep(UsedModel, *Volatility, stddev, *MR, NbPtsForRepliq, FullRepliq, upperProba, lowerProba);

		if (createConvAdjutst == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			convAdjustId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst);

			if (convAdjustId == RET_KO)
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(convAdjustId);

			return ARM_OK;
		}
		else
		{
			ARM_BSConvAdjustRep* oldConvAdjutst = (ARM_BSConvAdjustRep *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldConvAdjutst, ARM_BSCONVADJUSTREP) == 1)
			{
				if (oldConvAdjutst)
				{
					delete oldConvAdjutst;
					oldConvAdjutst = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst, objId);

				return ARM_OK;
			}
			else
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createConvAdjutst)
			delete createConvAdjutst;
		createConvAdjutst = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_ReplicConvAdjust_Create(
int Payoff_ReplicMode,
double Payoff_StepOrReplicPrecision,
int Payoff_StopMode,
double Payoff_StopThreshold,
int Sensi_ReplicMode,
double Sensi_StepOrReplicPrecision,
int Sensi_StopMode,
double Sensi_StopThreshold,
long UsedModelId,
double StrikeMinReplic,
double StrikeMaxReplic,
ARM_result& result,
long objId)
{
	long convAdjustId;

	ARM_Model* UsedModel = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	ARM_ReplicConvAdjust* createConvAdjutst = NULL;

	try
	{
		if (UsedModelId != ARM_NULL_OBJECT)
        {
            UsedModel = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(UsedModelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(UsedModel, ARM_MODEL) == 0)
		    {			
			    result.setMsg ("ARM_ERR: Used  is not of a good type");
			    return ARM_KO;
		    }
        }	

		createConvAdjutst = new ARM_ReplicConvAdjust(
			UsedModel,
			Payoff_ReplicMode,
			Payoff_StepOrReplicPrecision,
			Payoff_StopMode,
			Payoff_StopThreshold,
			Sensi_ReplicMode,
			Sensi_StepOrReplicPrecision,
			Sensi_StopMode,
			Sensi_StopThreshold,
			StrikeMinReplic,
			StrikeMaxReplic);

		if (createConvAdjutst == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			convAdjustId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst);

			if (convAdjustId == RET_KO)
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(convAdjustId);

			return ARM_OK;
		}
		else
		{
			ARM_ReplicConvAdjust* oldConvAdjutst = (ARM_ReplicConvAdjust *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldConvAdjutst, ARM_REPLICCONVADJUST) == 1)
			{
				if (oldConvAdjutst)
				{
					delete oldConvAdjutst;
					oldConvAdjutst = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createConvAdjutst, objId);

				return ARM_OK;
			}
			else
			{
				if (createConvAdjutst)
					delete createConvAdjutst;
				createConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createConvAdjutst)
			delete createConvAdjutst;
		createConvAdjutst = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_MapConvAdjust_Create(
                            long LiborArrearAdjId,
                            long NaturalCMSAdjId,
                            long PaymentLagAdjId,
                            ARM_result& result,
                            long objId)
{
	long mapConvAdjustId;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	ARM_MapConvAdjust* createMapConvAdjutst = NULL;

	ARM_ConvAdjustManager* liborArrearConvAdj = NULL;
	ARM_ConvAdjustManager* naturalCMSConvAdj = NULL;
	ARM_ConvAdjustManager* paymentLagConvAdj = NULL;

	try
	{
		if (LiborArrearAdjId != ARM_NULL_OBJECT)
		{
			liborArrearConvAdj = (ARM_ConvAdjustManager*) (LOCAL_PERSISTENT_OBJECTS->GetObject(LiborArrearAdjId));
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(liborArrearConvAdj, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(liborArrearConvAdj, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(liborArrearConvAdj, ARM_BSCONVADJUSTREP) == 0))
		    {

			    result.setMsg ("ARM_ERR: Libor Arrear Adjustment Manager  is not of a good type");
			    return ARM_KO;
		    }
		}
		
		if (NaturalCMSAdjId != ARM_NULL_OBJECT)
		{
			naturalCMSConvAdj = (ARM_ConvAdjustManager*) (LOCAL_PERSISTENT_OBJECTS->GetObject(NaturalCMSAdjId));
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(naturalCMSConvAdj, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(naturalCMSConvAdj, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(naturalCMSConvAdj, ARM_BSCONVADJUSTREP) == 0))
		    {

			    result.setMsg ("ARM_ERR: Natural CMS Adjustment Manager  is not of a good type");
			    return ARM_KO;
		    }
		}

		if (PaymentLagAdjId != ARM_NULL_OBJECT)
		{
			paymentLagConvAdj = (ARM_ConvAdjustManager*) (LOCAL_PERSISTENT_OBJECTS->GetObject(PaymentLagAdjId));
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(paymentLagConvAdj, ARM_BSCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(paymentLagConvAdj, ARM_REPLICCONVADJUST) == 0)
				&& (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(paymentLagConvAdj, ARM_BSCONVADJUSTREP) == 0))
		    {

			    result.setMsg ("ARM_ERR: Payment Lag Adjustment Manager  is not of a good type");
			    return ARM_KO;
		    }
		}

		createMapConvAdjutst = new ARM_MapConvAdjust(liborArrearConvAdj,
													naturalCMSConvAdj,
													paymentLagConvAdj);

		if (createMapConvAdjutst == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			mapConvAdjustId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createMapConvAdjutst);

			if (mapConvAdjustId == RET_KO)
			{
				if (createMapConvAdjutst)
					delete createMapConvAdjutst;
				createMapConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mapConvAdjustId);

			return ARM_OK;
		}
		else
		{
			ARM_MapConvAdjust* oldMapConvAdjutst = (ARM_MapConvAdjust *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);


			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldMapConvAdjutst, ARM_MAPCONVADJUST) == 1)
			{
				if (oldMapConvAdjutst)
				{
					delete oldMapConvAdjutst;
					oldMapConvAdjutst = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createMapConvAdjutst, objId);

				return ARM_OK;
			}
			else
			{
				if (createMapConvAdjutst)
					delete createMapConvAdjutst;
				createMapConvAdjutst = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createMapConvAdjutst)
			delete createMapConvAdjutst;
		createMapConvAdjutst = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_ReplicModel_Create(int ReplicMode,
								 double StepOrReplicPrecision,
								 int StopMode,
								 double StopThreshold,
								 int SensiReplicMode,
								 double SensiStepOrReplicPrecision,
								 int SensiStopMode,
								 double SensiStopThreshold,
								 long UsedModelId,
								 int SUMMITFormulaeUsed,
								 double StrikeMinReplic,
								 double StrikeMaxReplic,
								 ARM_result& result,
								 long objId)
{
	long replicModelId;

	ARM_Model* UsedModel = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");

	   return(ARM_KO);
	}

	CCString msg ("");

	ARM_ReplicModel* createReplicModel = NULL;

	try
	{
        UsedModel = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(UsedModelId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(UsedModel, ARM_MODEL) == 0)
		{			
			result.setMsg ("ARM_ERR: Used  is not of a good type");
			return ARM_KO;
		}	

		createReplicModel = new ARM_ReplicModel(
			UsedModel,
			ReplicMode,
			StepOrReplicPrecision,
			StopMode,
			StopThreshold,
            SensiReplicMode,
			SensiStepOrReplicPrecision,
			SensiStopMode,
			SensiStopThreshold,
            SUMMITFormulaeUsed,
			StrikeMinReplic,
			StrikeMaxReplic);

		if (createReplicModel == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			replicModelId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createReplicModel);

			if (replicModelId == RET_KO)
			{
				if (createReplicModel)
					delete createReplicModel;
				createReplicModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(replicModelId);

			return ARM_OK;
		}
		else
		{
			ARM_ReplicModel* oldReplicModel = (ARM_ReplicModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldReplicModel, ARM_REPLICMOD) == 1)
			{
				if (oldReplicModel)
				{
					delete oldReplicModel;
					oldReplicModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createReplicModel, objId);

				return ARM_OK;
			}
			else
			{
				if (oldReplicModel)
					delete oldReplicModel;
				oldReplicModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		if (createReplicModel)
			delete createReplicModel;
		createReplicModel = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_SetReplicDebugMode(int DebugMode,
								 ARM_result& result)
{
    CCString msg ("");

	try
	{
        ARM_ReplicPortfolio::SetDebugMode(DebugMode);

        string txt( "Replic Debug Mode:" );
		txt += DebugMode? "On" : "Off";
		result.setString(txt.c_str());

        return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_TriBSModel (long Model1Id,
								 long Model2Id,
								 long DiscModelId,
								 long Fx1DiscVolId,
								 long Fx2DiscVolId,
								 long Idx1Idx2CorrId,
								 long Idx1DiscIdxCorrId,
								 long Idx2DiscIdxCorrId,
								 long Idx1FxCorrId,
                                 long Idx2FxCorrId,
								 int quantoadjflag,
								 ARM_result& result,
								 long objId)
{
	long modId;

    ARM_VolCurve* Idx1Idx2Corr = NULL;
    ARM_VolCurve* Idx1DiscIdxCorr = NULL;
    ARM_VolCurve* Idx2DiscIdxCorr = NULL;
    
    ARM_VolCurve* Idx1FxCorr = NULL;    // Correlation between Index1 and forex Ccy1/DiscountCcy
    ARM_VolCurve* Idx2FxCorr = NULL;    // Correlation between Index2 and forex Ccy2/DiscountCcy
    // Vol Forex    
    ARM_VolCurve* Ccy1DiscFxVol = NULL;
    ARM_VolCurve* Ccy2DiscFxVol = NULL;

	ARM_BSModel* Model1 = NULL;
    ARM_BSModel* Model2 = NULL;
    ARM_BSModel* DiscountModel = NULL;
	

	ARM_TRIBSModel* TriBSModel = NULL;
	ARM_TRIBSModel* prevModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
                                 
        Ccy1DiscFxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fx1DiscVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Ccy1DiscFxVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Forex (BSModel1's ccy/Discount ccy) Volatility Curve is not of a good type");
			return ARM_KO;
		}
		
		                         
        Ccy2DiscFxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fx2DiscVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Ccy2DiscFxVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Forex (BSModel2's ccy/Discount ccy) Volatility Curve is not of a good type");
			return ARM_KO;
		}

		
        Idx1Idx2Corr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1Idx2CorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1Idx2Corr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/Index2 Correlation Curve is not of a good type");
			return ARM_KO;
		}

		Idx1DiscIdxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1DiscIdxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1DiscIdxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/DiscountBSModel'sIndex Correlation Curve is not of a good type");
			return ARM_KO;
		}

        Idx2DiscIdxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2DiscIdxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx2DiscIdxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index2/DiscountBSModel'sIndex Correlation Curve is not of a good type");
			return ARM_KO;
        }

		Idx1FxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1FxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1FxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/forex (ccy1/Disc) Correlation Curve is not of a good type");
			return ARM_KO;
		}

        Idx2FxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2FxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx2FxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index2/forex (ccy2/Disc) Correlation Curve is not of a good type");
			return ARM_KO;
		}

		Model1 = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Model1, ARM_BSMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: First BS Model is not of a good type");
			return ARM_KO;
		}
		
		Model2 = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model2Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Model2, ARM_BSMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Second BS Model is not of a good type");
			return ARM_KO;
		}
        
        if ( DiscModelId != ARM_NULL_OBJECT)
		{
		    DiscountModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscModelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(DiscountModel, ARM_BSMODEL) == 0)
		    {
			    result.setMsg ("ARM_ERR: Discount BS Model is not of a good type");
			    return ARM_KO;
		    }
        }
	

		TriBSModel = new ARM_TRIBSModel(Model1, Model2, DiscountModel, 
                                        Ccy1DiscFxVol, Ccy2DiscFxVol, Idx1Idx2Corr,
                                        Idx1DiscIdxCorr, Idx2DiscIdxCorr, Idx1FxCorr, Idx2FxCorr,
                                        quantoadjflag);
		if (TriBSModel == NULL)
		{
			result.setMsg ("ARM_ERR: TriBSModel is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriBSModel);

			if (modId == RET_KO)
			{
				if (TriBSModel)
					delete TriBSModel;
				TriBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_TRIBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ARM_TRIBSMODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriBSModel, objId);

				return ARM_OK;
			}
			else
			{
				if (TriBSModel)
					delete TriBSModel;
				TriBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (TriBSModel)
			delete TriBSModel;
		TriBSModel = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_SetDiscountPricingMode(long modelId,
									 long discPriceMode,
									 ARM_result& result)
{
    ARM_Model* model = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{                                 
		model = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(model, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}
		
		model->SetDiscPricingMode(discPriceMode);

		result.setLong(modelId);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}




extern long ARMLOCAL_TriBSDualModel (long Model1Id,
									 long Model2Id,
									 long DiscModelId,
									 long Fx1DiscVolId,
									 long Fx2DiscVolId,
									 long Idx1Idx2CorrId,
									 long Idx1DiscIdxCorrId,
									 long Idx2DiscIdxCorrId,
									 long Idx1FxCorrId,
									 long Idx2FxCorrId,
									 int quantoadjflag,
									 double correlforadj,
									 int withslopeflag,
									 ARM_result& result,
									 long objId)
{
    long modId;

    ARM_VolCurve* Idx1Idx2Corr = NULL;
    ARM_VolCurve* Idx1DiscIdxCorr = NULL;
    ARM_VolCurve* Idx2DiscIdxCorr = NULL;
    
    ARM_VolCurve* Idx1FxCorr = NULL;    // Correlation between Index1 and forex Ccy1/DiscountCcy
    ARM_VolCurve* Idx2FxCorr = NULL;    // Correlation between Index2 and forex Ccy2/DiscountCcy
    // Vol Forex    
    ARM_VolCurve* Ccy1DiscFxVol = NULL;
    ARM_VolCurve* Ccy2DiscFxVol = NULL;

    ARM_BSModel* Model1 = NULL;
    ARM_BSModel* Model2 = NULL;
    ARM_BSModel* DiscountModel = NULL;
	

    ARM_TRIBSDualModel* TriBSModel = NULL;
	ARM_TRIBSDualModel* prevModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
                                 
        Ccy1DiscFxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fx1DiscVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Ccy1DiscFxVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Forex (BSModel1's ccy/Discount ccy) Volatility Curve is not of a good type");
			return ARM_KO;
		}
		
		                         
        Ccy2DiscFxVol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fx2DiscVolId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Ccy2DiscFxVol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Forex (BSModel2's ccy/Discount ccy) Volatility Curve is not of a good type");
			return ARM_KO;
		}

		
        Idx1Idx2Corr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1Idx2CorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1Idx2Corr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/Index2 Correlation Curve is not of a good type");
			return ARM_KO;
		}

		Idx1DiscIdxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1DiscIdxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1DiscIdxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/DiscountBSModel'sIndex Correlation Curve is not of a good type");
			return ARM_KO;
		}

        Idx2DiscIdxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2DiscIdxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx2DiscIdxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index2/DiscountBSModel'sIndex Correlation Curve is not of a good type");
			return ARM_KO;
        }

		Idx1FxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1FxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1FxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index1/forex (ccy1/Disc) Correlation Curve is not of a good type");
			return ARM_KO;
		}

        Idx2FxCorr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2FxCorrId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx2FxCorr, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Index2/forex (ccy2/Disc) Correlation Curve is not of a good type");
			return ARM_KO;
		}

		Model1 = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model1, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: First BS Model is not of a good type");
			return ARM_KO;
		}
		
		Model2 = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model2Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model2, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Second BS Model is not of a good type");
			return ARM_KO;
		}
        
        if ( DiscModelId != ARM_NULL_OBJECT)
		{
		    DiscountModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(DiscModelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DiscountModel, ARM_MODEL) == 0)
		    {
			    result.setMsg ("ARM_ERR: Discount BS Model is not of a good type");
			    return ARM_KO;
		    }
        }
	

		TriBSModel = new ARM_TRIBSDualModel(Model1, Model2, DiscountModel, 
                                        Ccy1DiscFxVol, Ccy2DiscFxVol, Idx1Idx2Corr,
                                        Idx1DiscIdxCorr, Idx2DiscIdxCorr, Idx1FxCorr, Idx2FxCorr,
                                        quantoadjflag,correlforadj,withslopeflag);
		if (TriBSModel == NULL)
		{
			result.setMsg ("ARM_ERR: TriBSModel is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriBSModel);

			if (modId == RET_KO)
			{
				if (TriBSModel)
					delete TriBSModel;
				TriBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_TRIBSDualModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ARM_TRIBSDUALMODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriBSModel, objId);

				return ARM_OK;
			}
			else
			{
				if (TriBSModel)
					delete TriBSModel;
				TriBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (TriBSModel)
			delete TriBSModel;
		TriBSModel = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_CALIBRATORSFRM(long secId,
							 long capModelId,
							 long swoptModelId,
							 double meanRev,
							 const VECTOR<double>& calibParams,
							 long preInitFlagId,
							 const VECTOR<double>& initSigmaCrv,
							 const VECTOR<double>& initBetaOrShift,
							 long sigmaPfId,
							 long betaPfId,
							 long meanRevPfId,
							 long voltypeId,
							 long correlId,
                             const VECTOR<double>& SecurityParams,
                             const CCString c_tocalswaptATM,
							 ARM_result& result,
							 long objId)
{
    /// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

    /// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibratorSFRM* sfrmcalibrator = NULL;
	try
	{
        ARM_Model* capModel = NULL;
	    if( !GetObjectFromId( &capModel, capModelId, ARM_MODEL) )
	    {
            result.setMsg ("ARM_ERR: Cap Model is not of a good type");
			return ARM_KO;
	    }

        ARM_Model* swoptModel = NULL;
	    if( !GetObjectFromId( &swoptModel, swoptModelId, ARM_MODEL) )
	    {
            result.setMsg ("ARM_ERR: Swaption Model is not of a good type");
			return ARM_KO;
	    }

        ARM_Security* security = NULL;
        if( !GetObjectFromId( &security, secId, ARM_SECURITY) )
	    {
            result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
	    }

        ARM_StdPortfolio* sigmaPF   = NULL;
        if ( sigmaPfId != ARM_NULL_OBJECT)
		{
            if( !GetObjectFromId( &sigmaPF, sigmaPfId, ARM_PORTFOLIO) )
	        {
                result.setMsg ("ARM_ERR: sigma PF is not of a good type");
			    return ARM_KO;
	        }
        }

        ARM_StdPortfolio* betaPF   = NULL;
        if ( betaPfId != ARM_NULL_OBJECT)
		{
            if( !GetObjectFromId( &betaPF, betaPfId, ARM_PORTFOLIO) )
	        {
                result.setMsg ("ARM_ERR: beta PF is not of a good type");
			    return ARM_KO;
	        }
        }

        ARM_StdPortfolio* meanRevPF   = NULL;
        if ( meanRevPfId != ARM_NULL_OBJECT)
		{
            if( !GetObjectFromId( &meanRevPF, meanRevPfId, ARM_PORTFOLIO) )
	        {
                result.setMsg ("ARM_ERR: MRS PF is not of a good type");
			    return ARM_KO;
	        }
        }

        ARM_VolLInterpol* mCorrel = NULL;
        if ( correlId != ARM_NULL_OBJECT)
		{
            if( !GetObjectFromId( &mCorrel, correlId, ARM_VOL_LIN_INTERPOL) )
	        {
                result.setMsg ("ARM_ERR: correlation matrix is not of a good type");
			    return ARM_KO;
	        }
        }
        ARM_Matrix* matrixCorrel = NULL;
		if (mCorrel)
			matrixCorrel = mCorrel->GetVolatilities();
	
		if (calibParams.size() == 0)
		{
			result.setMsg ("ARM_ERR: calib param is a null size");
			return ARM_KO;
		}

        ARM_Vector* vCalibParams = CreateARMVectorFromVECTOR(calibParams);
        /// use auto_ptr for exception safety!
		CC_NS(std,auto_ptr)<ARM_Vector> autoptrCalibParams(vCalibParams);

        ARM_Vector* vInitSigmaCurv = NULL;	    
		if (initSigmaCrv.size() != 0)
			vInitSigmaCurv = CreateARMVectorFromVECTOR(initSigmaCrv);
        /// use auto_ptr for exception safety!
		CC_NS(std,auto_ptr)<ARM_Vector> autoptrInitSigmaCurv(vInitSigmaCurv);

        ARM_Vector* vInitBeta = NULL;
		if (initBetaOrShift.size() != 0)
			vInitBeta = CreateARMVectorFromVECTOR(initBetaOrShift);
        /// use auto_ptr for exception safety!
		CC_NS(std,auto_ptr)<ARM_Vector> autoptrInitBeta(vInitBeta);

        ARM_Vector* vSecurityParams= NULL;
        if (SecurityParams.size() != 0)
			vSecurityParams = CreateARMVectorFromVECTOR(SecurityParams);
        /// use auto_ptr for exception safety!
		CC_NS(std,auto_ptr)<ARM_Vector> autoptrSecurityParams(vSecurityParams);

        string calswaptATM = CCSTringToSTLString(c_tocalswaptATM);

		sfrmcalibrator = new ARM_CalibratorSFRM(security,
										  capModel,
										  swoptModel,
										  meanRev,
										  vCalibParams,
										  preInitFlagId,
										  vInitSigmaCurv,
										  vInitBeta,
										  sigmaPF,
										  betaPF,
										  meanRevPF,
										  voltypeId,
										  matrixCorrel,
                                          vSecurityParams,
                                          calswaptATM);

        /// assign object
		if( !assignObject( sfrmcalibrator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
        delete sfrmcalibrator;
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_SFRMCALIBRATE (long calibratorSFRMId,
                             long mktcapmodelId,
                             long mktswaptmodelId,
                             CCString tocalibrateBeta,
                             CCString tocalibrateMR,
                             long kerneltoGP,
							 ARM_result& result,
							 long objId)
{
	long frmModelId;

	ARM_CalibratorSFRM* Calibrator = NULL;
    ARM_Model* MktCapModel = NULL;
    ARM_Model* MktSwaptModel = NULL;

	ARM_FRMModel* newModel = NULL;
	ARM_FRMModel* oldModel = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		Calibrator = (ARM_CalibratorSFRM *) LOCAL_PERSISTENT_OBJECTS->GetObject(calibratorSFRMId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Calibrator, ARM_CALIBRATORSFRM) == 0)
		{
			result.setMsg ("ARM_ERR: Calibrator SFRM is not of a good type");
			return ARM_KO;
		}

        if ( mktcapmodelId != ARM_NULL_OBJECT)
        {
            MktCapModel = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(mktcapmodelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(MktCapModel, ARM_MODEL) == 0)
		    {
			    result.setMsg ("ARM_ERR: Cap Model is not of a good type");
			    return ARM_KO;
		    }
        }

        if ( mktswaptmodelId != ARM_NULL_OBJECT)
        {
		    MktSwaptModel = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(mktswaptmodelId);
		    if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(MktSwaptModel, ARM_MODEL) == 0)
		    {
			    result.setMsg ("ARM_ERR: Swaption Model is not of a good type");
			    return ARM_KO;
		    }
        }

        string calibrateBeta = CCSTringToSTLString(tocalibrateBeta);
        string calibrateMR = CCSTringToSTLString(tocalibrateMR);

        Calibrator->SetParamsToHedge(MktCapModel,MktSwaptModel,calibrateBeta,calibrateMR);
        newModel = Calibrator->Calibrate(NULL,MktCapModel==NULL || MktSwaptModel==NULL);
        
		if (newModel == NULL)
		{
			result.setMsg ("ARM_ERR: FRM Model is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			frmModelId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newModel);

			if (frmModelId == RET_KO)
			{
				if (newModel)
					delete newModel;
				newModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(frmModelId);

			return ARM_OK;
		}
		else
		{
			oldModel = (ARM_FRMModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldModel, ARM_FRMMODEL) == 1)
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
		x.DebugPrint();

		if (newModel)
			delete newModel;
		newModel = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_CRACALCULATOR(long secId,
							long zcId,
							long capVolATMId,
							long rhoId,
							long nuId,
							long swptVolId,
							long betaSABRId,
							double meanRev,
							const VECTOR<double>& calibParams,
							long preInitFlagId,
							const VECTOR<double>& initSigmaCrv,
							const VECTOR<double>& initBetaOrShift,
							long voltypeId,
                            const VECTOR<double>& SecurityParams,
							double horizon,
							long PathNumber,
							const VECTOR<double>& TreeParams,
							const VECTOR<double>& ReCalibFlags,
							CCString CalswaptATM,
							long rhoSwoptId,
							long nuSwoptId,
							long betaSwoptId,
                            long SABRSigmaOrAlphaFlag,
							ARM_result& result,
							long objId)
{
	long modId;

	ARM_CRAFRMMarkovTree* newMod = NULL;
	ARM_CRAFRMMarkovTree* oldMod = NULL;

    // What type of SABR
    int sigmaOrAlpha = int(SABRSigmaOrAlphaFlag);

	char sHorizon[11];
	
	ARM_ZeroCurve* zc = NULL;
	ARM_VolCurve* volcapATM = NULL;
	ARM_VolCurve* volswpt = NULL;
	ARM_VolCurve* rho = NULL;
	ARM_VolCurve* nu = NULL;
	ARM_VolCurve* betaSABR = NULL;
	ARM_VolCurve* rhoSwopt = NULL;
	ARM_VolCurve* nuSwopt = NULL;
	ARM_VolCurve* betaSwopt = NULL;

	ARM_Security* sec = NULL;	

	ARM_Vector* vInitSigmaCurv = NULL;
	ARM_Vector* vInitBetaOrShift = NULL;

	ARM_Vector* vCalibParams = NULL;
	ARM_Vector* vSecurityParams= NULL;
    ARM_Vector* vTreeParameters = NULL;
	ARM_Vector* vReCalibFlags = NULL;

	ARM_StdPortfolio* sigmaPF   = NULL;
	ARM_StdPortfolio* betaPF    = NULL;
	ARM_StdPortfolio* meanRevPF = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(horizon,sHorizon);

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
		   result.setMsg ("ARM_ERR: rate Zc Curve is not of a good type");
			
           return(ARM_KO);
		}
	
		volcapATM = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolATMId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volcapATM, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Cap Vol Curve is not of a good type");

			return(ARM_KO);
		}

		if ( rhoId != ARM_NULL_OBJECT)
        {
			rho = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(rho, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( nuId != ARM_NULL_OBJECT)
        {
			nu = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(nu, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( betaSABRId != ARM_NULL_OBJECT)
        {
			betaSABR = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSABRId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(betaSABR, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: beta Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		if ( rhoSwoptId != ARM_NULL_OBJECT)
        {
			rhoSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(rhoSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( nuSwoptId != ARM_NULL_OBJECT)
        {
			nuSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(nuSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( betaSwoptId != ARM_NULL_OBJECT)
        {
			betaSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(betaSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Beta Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		volswpt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swptVolId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volswpt, ARM_VOL_CURVE) == 0 )
		{
			result.setMsg ("ARM_ERR: Swaption Vol Curve is not of a good type");

			return(ARM_KO);
		}

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec, ARM_OPTIONPORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not an OptionPortfolio type ");
			return ARM_KO;
		}

		if (calibParams.size() == 0)
		{
			result.setMsg ("ARM_ERR: calib param is a null size");
			return ARM_KO;
		}

		vCalibParams = CreateARMVectorFromVECTOR(calibParams);

		if (initSigmaCrv.size() != 0)
			vInitSigmaCurv = CreateARMVectorFromVECTOR(initSigmaCrv);

		if (initBetaOrShift.size() != 0)
			vInitBetaOrShift = CreateARMVectorFromVECTOR(initBetaOrShift);

        
        if (SecurityParams.size() != 0)
			vSecurityParams = CreateARMVectorFromVECTOR(SecurityParams);

		
		vTreeParameters = CreateARMVectorFromVECTOR(TreeParams);

		if (ReCalibFlags.size() == 0)
		{
			result.setMsg ("ARM_ERR: Re calib flag is a null size");
			return ARM_KO;
		}

		vReCalibFlags = CreateARMVectorFromVECTOR(ReCalibFlags);

		string ArgcalswaptATM = CCSTringToSTLString(CalswaptATM);

        newMod = new ARM_CRAFRMMarkovTree(sec,
										zc->GetAsOfDate(),
										zc,
										volcapATM,
										rho,
										nu,
										volswpt,
										meanRev,
										vCalibParams,
										preInitFlagId,
										vInitSigmaCurv,
										vInitBetaOrShift,
										voltypeId,
										vSecurityParams,
										(ARM_Date)sHorizon,
										PathNumber,
										vTreeParameters,
										vReCalibFlags,
										ArgcalswaptATM,
										betaSABR,
										rhoSwopt,
										nuSwopt,
										betaSwopt,
                                        sigmaOrAlpha);

        delete vCalibParams;
		vCalibParams=NULL;
		delete vSecurityParams;
		vSecurityParams=NULL;
        delete vTreeParameters;
		vTreeParameters = NULL;
		delete vReCalibFlags;
		vReCalibFlags = NULL;


        if (newMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if(objId == -1)
		{
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod);

			if (modId == RET_KO)
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			oldMod = (ARM_CRAFRMMarkovTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldMod, ARM_FRMMARKOVTREE) == 1)
			{
				if (oldMod)
				{
					delete oldMod;
					oldMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod, objId);

				return ARM_OK;
			}
			else
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		delete newMod;
		newMod = NULL;

		delete vCalibParams;
		vCalibParams=NULL;
		
		delete vSecurityParams;
		vSecurityParams=NULL;
		
		delete vTreeParameters;
		vTreeParameters = NULL;

		delete vReCalibFlags;
		vReCalibFlags = NULL;


		ARM_RESULT();
	}
}



long ARMLOCAL_BUMPCRACALCULATOR(long cracalcId,
								long zcId,
								long capVolATMId,
								long rhoId,
								long nuId,
								long swptVolId,
								long betaSABRId,
								long rhoSwoptId,
								long nuSwoptId,
								long betaSwoptId,
                                long SABRSigmaOrAlphaFlag,
								ARM_result& result,
								long objId)
{

	ARM_ZeroCurve* zc = NULL;
	ARM_VolCurve* volcapATM = NULL;
	ARM_VolCurve* volswpt = NULL;
	ARM_VolCurve* rho = NULL;
	ARM_VolCurve* nu = NULL;
	ARM_VolCurve* betaSABR = NULL;
	ARM_VolCurve* rhoSwopt = NULL;
	ARM_VolCurve* nuSwopt = NULL;
	ARM_VolCurve* betaSwopt = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

        // What type of SABR
    int sigmaOrAlpha = int(SABRSigmaOrAlphaFlag);

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
		   result.setMsg ("ARM_ERR: rate Zc Curve is not of a good type");
			
           return(ARM_KO);
		}
	
		volcapATM = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(capVolATMId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volcapATM, ARM_VOL_CURVE) == 0 )
		{
		  result.setMsg ("ARM_ERR: Cap Vol Curve is not of a good type");
			
          return(ARM_KO);
		}

		if ( rhoId != ARM_NULL_OBJECT)
        {
			rho = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(rho, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( nuId != ARM_NULL_OBJECT)
        {
			nu = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(nu, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		if ( betaSABRId != ARM_NULL_OBJECT)
        {
			betaSABR = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSABRId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(betaSABR, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: beta Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		if ( rhoSwoptId != ARM_NULL_OBJECT)
        {
			rhoSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(rhoSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(rhoSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Rho Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( nuSwoptId != ARM_NULL_OBJECT)
        {
			nuSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(nuSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(nuSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Nu Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}
		if ( betaSwoptId != ARM_NULL_OBJECT)
        {
			betaSwopt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(betaSwoptId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(betaSwopt, ARM_VOL_CURVE) == 0 )
			{
				result.setMsg ("ARM_ERR: Beta Swaption Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		volswpt = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(swptVolId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volswpt, ARM_VOL_CURVE) == 0 )
		{
		  result.setMsg ("ARM_ERR: Swaption Vol Curve is not of a good type");
			
          return(ARM_KO);
		}

		ARM_CRAFRMMarkovTree* Cracalculator = (ARM_CRAFRMMarkovTree *) LOCAL_PERSISTENT_OBJECTS->GetObject(cracalcId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Cracalculator, ARM_FRMMARKOVTREE) != 0)
		{
			Cracalculator->Bump(zc,
								volcapATM,
								rho,
								nu,
								volswpt,
								betaSABR,
								rhoSwopt,
								nuSwopt,
								betaSwopt,
                                sigmaOrAlpha);
		}
		else
		{
            // We have to manage SABR/Alpha for the GP calculator later!!!!

            ARM_CRACalculator* Cracalculator = dynamic_cast<ARM_CRACalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cracalcId));
			
            if (Cracalculator == NULL)
			{
				result.setMsg ("ARM_ERR: CRA calculator is not of a good type");

				return(ARM_KO);
			}

			Cracalculator->Bump(zc,
								volcapATM,
								rho,
								nu,
								volswpt,
								betaSABR,
								rhoSwopt,
								nuSwopt,
								betaSwopt);
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



long ARMLOCAL_COMPUTEFWDBSDELTA(double fwd, double strike, double vol,
                                double T, 
                                int CallPut, 
								ARM_result& result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		double delta = 100000; 
        
        delta = ComputeFwdBSDelta(fwd, strike, vol/100.0, T, CallPut)*100.0;

		result.setDouble(delta);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



long ARMLOCAL_ComputeSplinedSigmaATMF(VECTOR<double>& deltas, VECTOR<double>& sigmas, 
									  double matu, double SigmaZDS, double Precision,
                                      double FX_SPOT,
									  ARM_result& result)
{
    ARM_Vector* mktDeltas = NULL;
	ARM_Vector* mktSigmas = NULL;
	
	/// input checks
	if (!GlobalPersistanceOk( result ) )
	   return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
	   return ARM_KO;
	}
	
	try
	{
		mktDeltas = CreateARMVectorFromVECTOR(deltas);
		mktSigmas = CreateARMVectorFromVECTOR(sigmas);

		double sigmaATMF = 100000.0; 
		
        sigmaATMF = ComputeSplinedSigmaATMF(mktDeltas, mktSigmas,
                                            matu, 
                                            SigmaZDS, 
                                            Precision,
                                            FX_SPOT);

		result.setDouble(sigmaATMF);

        if (mktDeltas)
           delete mktDeltas;

        if (mktSigmas)
           delete mktSigmas;

		return(ARM_OK);

	}
		
	catch(Exception& x)
	{	
        if (mktDeltas)
           delete mktDeltas;

        if (mktSigmas)
           delete mktSigmas;
		
        x.DebugPrint();
		
		ARM_RESULT();
	}

}



long ARMLOCAL_ComputeDeltaFwdFromDeltaWP(double AsOf, 									
                                         double matu,
                                         double sigma,
                                         double fxSpot,
                                         double deltaWithPremium,
									     long domCrvId,
                                         long foreignCrvId,
									     ARM_result& result)
{			
	/// input checks
	if ( !GlobalPersistanceOk( result ) )
	   return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( AsOf, charDate );
		ARM_Date dateIn( charDate );
		
		ARM_ZeroCurve* domCrv=NULL;
		ARM_ZeroCurve* foreignCrv=NULL;

		if ( !GetObjectFromId( &domCrv, domCrvId, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Zero curve is not of a good type");
			return ARM_KO;
		}

		if ( !GetObjectFromId( &foreignCrv, foreignCrvId, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Zero curve is not of a good type");
			return ARM_KO;
		}

		double value = 100000; 

		value = ComputeDeltaFwdFromDeltaWP(dateIn, matu, sigma,
                                           fxSpot, 
                                           deltaWithPremium,
										   domCrv, 
                                           foreignCrv);
		
		result.setDouble(value);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

long ARMLOCAL_CalculateImpliedStrikeFromDeltaWithPremium(double AsOf, 									
                                         double matu,
                                         double sigma,
                                         double fxSpot,
                                         double deltaWithPremium,
									     long domCrvId,
                                         long foreignCrvId,
									     ARM_result& result)
{			
	/// input checks
	if ( !GlobalPersistanceOk( result ) )
	   return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		char  charDate[11];
		Local_XLDATE2ARMDATE( AsOf, charDate );
		ARM_Date dateIn( charDate );
		
		ARM_ZeroCurve* domCrv=NULL;
		ARM_ZeroCurve* foreignCrv=NULL;

		if ( !GetObjectFromId( &domCrv, domCrvId, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Zero curve is not of a good type");
			return ARM_KO;
		}

		if ( !GetObjectFromId( &foreignCrv, foreignCrvId, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Zero curve is not of a good type");
			return ARM_KO;
		}

		double value = 100000; 

		value = CalculateImpliedStrikeFromDeltaWithPremium
										(dateIn, matu, sigma,
                                           fxSpot, 
                                           deltaWithPremium,
										   domCrv, 
                                           foreignCrv);
		
		result.setDouble(value);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

long ARMLOCAL_CalcFwdFXSpot(double AsofDate,
							double Spot,
                            double aFwdDate,
                            long NumDiscountCurveID, //domestic curve
                            long UndDiscountCurveID, //foreign curve
							ARM_result& result)
{			
	/// input checks
	if ( !GlobalPersistanceOk( result ) )
	   return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		double value = 0;
		char  charDateAsOf[11];
		char  charDateFwd[11];

		Local_XLDATE2ARMDATE(AsofDate, charDateAsOf);
		Local_XLDATE2ARMDATE(aFwdDate, charDateFwd);

		ARM_Date AsOfDate(charDateAsOf);
		ARM_Date FwdDate(charDateFwd);
		
		ARM_ZeroCurve* NumDiscountCurve=NULL;
		ARM_ZeroCurve* UndDiscountCurve=NULL;

		if ( !GetObjectFromId(&NumDiscountCurve, NumDiscountCurveID, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Domestic curve is not of a good type");
			return ARM_KO;
		}

		if ( !GetObjectFromId( &UndDiscountCurve, UndDiscountCurveID, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Foreign curve is not of a good type");
			return ARM_KO;
		} 

		value = CalcFwdFXSpot(AsOfDate,
							  Spot,
                              FwdDate,
                              NumDiscountCurve,
                              UndDiscountCurve);
		
		result.setDouble(value);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}

}




long ARMLOCAL_GetQModelQVolatility(long QModelId,
	ARM_result& result)
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	CCString msg ("");
	ARM_Q_Model* Qmod = NULL;

	try
	{
		Qmod = (ARM_Q_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(QModelId);

		/// An inflation curve should be of the base type ZERO_CURVE
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Qmod, ARM_Q_MODEL ) == 0)
		{
			result.setMsg ("ARM_ERR: Q Model is not of a good type");
			return ARM_KO;
		}

		result.setDouble(Qmod->GetQVolatility());
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_FromSigmaToAlpha(long bssmilemodId,
							   long sigmaCurveId,
							   double strike,
							   ARM_result& result,
							   long objId)
{
	long curveId;

	ARM_VolCurve* newCurve = NULL;
	ARM_VolCurve* oldCurve = NULL;


	ARM_BSSmiledModel* bssmilemod = NULL;
	ARM_VolCurve* sigmaCurve = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg ("ARM_ERR: Pb with accessing objects");
		
       return ARM_KO;
	}

	CCString msg ("");

	try
	{
		bssmilemod = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(bssmilemodId);
		
        if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bssmilemod, ARM_BSSMILEDMODEL) == 0 )
		{
			result.setMsg ("ARM_ERR: BS Smiled Model Curve is not of a good type");
				
			return(ARM_KO);
		}

		if ( sigmaCurveId != ARM_NULL_OBJECT)
        {
			sigmaCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(sigmaCurveId);
			
            if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sigmaCurve, ARM_VOL_LIN_INTERPOL) == 0 )
			{
				result.setMsg ("ARM_ERR: Sigma Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		newCurve = bssmilemod->FromSigmaToAlpha(sigmaCurve, strike);

		if ( newCurve == NULL )
		{
			result.setMsg ("ARM_ERR: FromSigmaToAlpha returned a NULL Curve");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) newCurve);

			if ( curveId == RET_KO )
			{
				if (newCurve)
					delete newCurve;
				newCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			oldCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldCurve, ARM_VOL_LIN_INTERPOL) == 1 )
			{
				if (oldCurve)
				{
				   delete oldCurve;
					
                   oldCurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newCurve, objId);

				return ARM_OK;
			}
			else
			{
				if (newCurve)
				   delete newCurve;
				newCurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newCurve)
		   delete newCurve;
		newCurve = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FromAlphaToSigma(long bssmilemodId,
							   long alphaCurveId,
							   double strike,
							   ARM_result& result,
							   long objId)
{
	long curveId;

	ARM_BSSmiledModel* bssmilemod = NULL;

	ARM_VolCurve* newCurve = NULL;
	ARM_VolCurve* oldCurve = NULL;
	
	ARM_VolCurve* alphaCurve = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		bssmilemod = (ARM_BSSmiledModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(bssmilemodId);

        if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bssmilemod, ARM_BSSMILEDMODEL) == 0 )
		{
			result.setMsg ("ARM_ERR: BS Smiled Model Curve is not of a good type");
				
			return(ARM_KO);
		}

		if ( alphaCurveId != ARM_NULL_OBJECT)
        {
			alphaCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(alphaCurveId);
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(alphaCurve, ARM_VOL_LIN_INTERPOL) == 0 )
			{
				result.setMsg ("ARM_ERR: alpha Curve is not of a good type");
					
				return(ARM_KO);
			}
		}

		newCurve = bssmilemod->FromAlphaToSigma(alphaCurve, strike);

		if ( newCurve == NULL )
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		CREATE_GLOBAL_OBJECT();

		if ( objId == -1 )
		{
			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurve);

			if (curveId == RET_KO)
			{
				if (newCurve)
					delete newCurve;
				newCurve = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			oldCurve = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldCurve, ARM_VOL_LIN_INTERPOL) == 1)
			{
				if (oldCurve)
				{
					delete oldCurve;
					oldCurve = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCurve, objId);

				return ARM_OK;
			}
			else
			{
				if (newCurve)
					delete newCurve;
				newCurve = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if(newCurve)
			delete newCurve;
		newCurve = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_TriXBSModel (	long Model1Id,
							long Model2Id,
							long Model3Id,
							long Idx1Idx2CorrCurveId,
							long Idx1Idx3CorrCurveId,
							long Idx2Idx3CorrCurveId,
							ARM_result& result,
							long objId)
{
	long modId;

	ARM_Model* Model1 = NULL;
    ARM_Model* Model2 = NULL;
    ARM_Model* Model3 = NULL;

    ARM_VolCurve* Idx1Idx2Corr = NULL;
    ARM_VolCurve* Idx1Idx3Corr = NULL;
    ARM_VolCurve* Idx2Idx3Corr = NULL;

	ARM_TriXBSModel* prevModel = NULL;
	ARM_TriXBSModel* TriXBSModel = NULL;
    
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
        Model1 = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model1Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model1, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: First Model is not of a good type");
			return ARM_KO;
		}
		
		Model2 = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model2Id);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model2, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Second Model is not of a good type");
			return ARM_KO;
		}
		
		if ( Model3Id != ARM_NULL_OBJECT)
		{
			Model3 = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(Model3Id);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Model3, ARM_MODEL) == 0)
			{
				result.setMsg ("ARM_ERR: Third Model is not of a good type");
				return ARM_KO;
			}
		}
		
		if ( Idx1Idx2CorrCurveId != ARM_NULL_OBJECT)
		{
			Idx1Idx2Corr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1Idx2CorrCurveId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1Idx2Corr, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1/Index2 Correlation Curve is not of a good type");
				return ARM_KO;
			}
		}

		if ( Idx1Idx3CorrCurveId != ARM_NULL_OBJECT)
		{
			Idx1Idx3Corr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1Idx3CorrCurveId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx1Idx3Corr, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1/Index3 Correlation Curve is not of a good type");
				return ARM_KO;
			}
		}

		if ( Idx2Idx3CorrCurveId != ARM_NULL_OBJECT)
		{
			Idx2Idx3Corr = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2Idx3CorrCurveId);
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Idx2Idx3Corr, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Index2/Index3 Correlation Curve is not of a good type");
				return ARM_KO;
			}
		}

		TriXBSModel = new ARM_TriXBSModel(	Model1, Model2, Model3, 
											Idx1Idx2Corr, Idx1Idx3Corr, Idx2Idx3Corr);

		if (TriXBSModel == NULL)
		{
			result.setMsg ("ARM_ERR: TriXBSModel is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriXBSModel);

			if (modId == RET_KO)
			{
				if (TriXBSModel)
					delete TriXBSModel;
				TriXBSModel = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevModel = (ARM_TriXBSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevModel, ARM_TRIXBSMODEL) == 1)
			{
				if (prevModel)
				{
					delete prevModel;
					prevModel = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)TriXBSModel, objId);

				return ARM_OK;
			}
			else
			{
				if (TriXBSModel)
					delete TriXBSModel;
				TriXBSModel = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (TriXBSModel)
			delete TriXBSModel;
		TriXBSModel = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GetAsOf(long modelId,
					  ARM_result& result)
{			
	/// input checks
	if ( !GlobalPersistanceOk( result ) )
	   return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	ARM_Model* model;

	try
	{
		if ( !GetObjectFromId(&model, modelId, ARM_MODEL) )
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}
		
		result.setDouble(model->GetStartDateJul());

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}

}


/*--------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
