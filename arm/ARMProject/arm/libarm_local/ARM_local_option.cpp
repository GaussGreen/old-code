#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_option.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"

#include <inst\security.h>
#include <inst\option.h>
#include <inst\portfolio.h>
#include <inst\optionportfolio.h>
#include <inst\sumopt.h>
#include <inst\swaption_smile.h>

#include <util\exercise.h>

#include <mod\bsxtic.h>
#include <mod\bscorrmodel.h>
#include <mod\calibratorfrm.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"

#include "gpbase/datestrip.h"
#include "gpbase/curve.h"
#include "gpbase/curveconvert.h"
#include "stripoption.h"
#include "stripdigitaloption.h"
#include "gpbase/interpolator.h"

using ARM::ARM_LinearInterpolatorCstExtrapol;
using ARM::RefValueToCurve;

long ARMLOCAL_EXOPTION (long underId,
						long optionType,
						long styleId,
						long kRefValId,
						ARM_result& result,
						long objId)
{
	long optId;

    ARM_Security* sec=NULL;
    ARM_Option* option=NULL;
    ARM_Option* newOption=NULL;
    ARM_ExerciseStyle *xStyle=NULL;
    ARM_ReferenceValue *kRefVal=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(styleId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: xStyle is not of a good type");
			return ARM_KO;
		}

		kRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(kRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		newOption = new ARM_Option(sec, optionType, kRefVal, xStyle);

		if (newOption == NULL)
		{
			result.setMsg ("ARM_ERR: Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOption);

			if (optId == RET_KO)
			{
				if (newOption)
					delete newOption;
				newOption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(optId);

			return ARM_OK;
		}
		else
		{
		    option = (ARM_Option *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(option, ARM_OPTION) == 1)
			{
				if (option)
					delete option;
				option = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOption, objId);

				return ARM_OK;
			}
			else
			{
				if (newOption)
					delete newOption;
				newOption = NULL;

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



long ARMLOCAL_OPTION(long underId,
					 double maturityDate,
					 double strike,
					 long optionType,
					 long exerciseType,
					 long strikeType,
					 double FstXDate,
					 double PayDate,
					 ARM_result& result,
					 long objId)
{
	long optId;

	ARM_Security* sec      = NULL;
	ARM_Option* option     = NULL;
    ARM_Option* prevOption = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char sMaturDate[11];
	char sFstXDate[11];
	char sPayDate[11];

	CCString msg("");

	try
	{
		Local_XLDATE2ARMDATE(maturityDate, sMaturDate);

		ARM_Date myFstXDate(ARM_DEFAULT_DATE);

		if ( FstXDate != -1 )
		{
		   Local_XLDATE2ARMDATE(FstXDate,sFstXDate);
			
           myFstXDate = (ARM_Date) sFstXDate;
		}

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

		if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0 )
		{
		   result.setMsg("ARM_ERR: Underlying Security is not of a good type");
			
           return(ARM_KO);
		}

		ARM_Date myPayDate(ARM_DEFAULT_DATE);

		if ( PayDate != -1 )
		{
		   Local_XLDATE2ARMDATE(PayDate, sPayDate);
			
           myPayDate = (ARM_Date) sPayDate;
		}

        option = new ARM_Option(sec, (ARM_Date) sMaturDate,
								strike, optionType, exerciseType, strikeType,
								myFstXDate);

        if ( option == NULL )
		{
		   result.setMsg("ARM_ERR: Option is null");
				
           return(ARM_KO);
		}

        if ( PayDate != -1 )
           option->SetPayDate(myPayDate);

        if ( objId == -1 )
		{
		   CREATE_GLOBAL_OBJECT();

		   optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) option);

		   if ( optId == RET_KO )
           {
			  if (option)
				 delete option;
			  option = NULL;

			  result.setMsg("ARM_ERR: Pb with inserting object");				
				
              return(ARM_KO);
           }

		   result.setLong(optId);

		   return(ARM_OK);
		}
		else
		{
			prevOption = (ARM_Option *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
 
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevOption, ARM_OPTION) == 1 )
			{
			   if (prevOption)
               {
                  delete prevOption;

                  prevOption = NULL;
               }

               	LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) option, objId);

				return(ARM_OK);
			}
            else
            {
               if (option)
				  delete option;
			   option = NULL;

			   result.setMsg ("ARM_ERR: previous object is not of a good type");
				
               return(ARM_KO);
            }
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (option)
			delete option;
		option = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_VolImp (long secId,
					  long modId,
					  double price,
					  ARM_result& result)
{
	double vol;

	ARM_Model* mod=NULL;
	ARM_Security* option=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		option = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(option, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: security is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		option->SetModel(mod);

		vol = option->ComputeImpliedVol(price);

		result.setDouble(vol);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_bsOption (double spot,
						double strike,
						double volatility,
						double dividend,
						double discountRate,
						double maturity,
						double CallPut,
						ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;

	try
	{
		bsOpt = bsOption(spot/100.,
						 strike/100.,
						 volatility/100.,
						 dividend/100.,
						 discountRate/100.,
						 maturity,
						 CallPut);

			result.setDouble(bsOpt);

			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_bsDelta (double spot,
					   double strike,
					   double volatility,
					   double dividend,
					   double discountRate,
					   double maturity,
					   double CallPut,
					   ARM_result& result)
{

	CCString msg ("");
	double bsDel=0;

	try
	{
		bsDel = bsDelta(spot/100.,
						strike/100.,
						volatility/100.,
						dividend/100.,
						discountRate/100.,
						maturity,
						CallPut);

			result.setDouble(bsDel);

			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_bsVega (double spot,
					  double strike,
					  double volatility,
					  double dividend,
					  double discountRate,
					  double maturity,
					  double CallPut,
					  ARM_result& result)
{

	CCString msg ("");
	double bsVeg=0;

	try
	{
		bsVeg = bsVega(spot/100.,
					   strike/100.,
					   volatility/100.,
					   dividend/100.,
					   discountRate/100.,
					   maturity,
					   CallPut);

			result.setDouble(bsVeg);

			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_bsGamma (double spot,
					   double strike,
					   double volatility,
					   double dividend,
					   double discountRate,
					   double maturity,
					   double CallPut,
					   ARM_result& result)
{

	CCString msg ("");
	double bsGam=0;

	try
	{
		bsGam = bsGamma(spot/100.,
						strike/100.,
						volatility/100.,
						dividend/100.,
						discountRate/100.,
						maturity,
						CallPut);

			result.setDouble(bsGam);

			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_bsTheta (double spot,
					   double strike,
					   double volatility,
					   double dividend,
					   double discountRate,
					   double maturity,
					   double CallPut,
					   ARM_result& result)
{

	CCString msg ("");
	double bsThe=0;

	try
	{
		bsThe = bsTheta(spot/100.,
						strike/100.,
						volatility/100.,
						dividend/100.,
						discountRate/100.,
						maturity,
						CallPut);

			result.setDouble(bsThe);

			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_ARM_OPTIONPORTFOLIO (long portfolioId,
								   long styleId,
								   long kRefValId,
								   long optionType,
								   ARM_result& result,
								   long objId)
{
	long optId;

    ARM_Portfolio* pf=NULL;
    ARM_OptionPortfolio* option=NULL;
    ARM_OptionPortfolio* newOption=NULL;
    ARM_ExerciseStyle *xStyle=NULL;
    ARM_ReferenceValue *kRefVal=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(portfolioId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}

		xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(styleId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: xStyle is not of a good type");
			return ARM_KO;
		}

		if (kRefValId != ARM_NULL_OBJECT)
		{
			kRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kRefValId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(kRefVal, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: refValue is not of a good type");
				return ARM_KO;
			}
		}

		newOption = new ARM_OptionPortfolio(pf, xStyle, kRefVal, optionType);

		if (newOption == NULL)
		{
			result.setMsg ("ARM_ERR: Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOption);

			if (optId == RET_KO)
			{
				if (newOption)
					delete newOption;
				newOption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(optId);

			return ARM_OK;
		}
		else
		{
		    option = (ARM_OptionPortfolio *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(option, ARM_OPTIONPORTFOLIO) == 1)
			{
				if (option)
					delete option;
				option = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newOption, objId);

				return ARM_OK;
			}
			else
			{
				if (newOption)
					delete newOption;
				newOption = NULL;

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

long ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (long calibratorId,
									            CCString portfolioType,
									            ARM_result& result,
									            long objId )

{
	long pfId;

    ARM_CalibratorSFRM * calibrator=NULL;
    ARM_StdPortfolio* pf=NULL;
    ARM_StdPortfolio* newpf=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		calibrator = (ARM_CalibratorSFRM *) LOCAL_PERSISTENT_OBJECTS->GetObject(calibratorId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(calibrator, ARM_CALIBRATORSFRM) == 0)
		{
			result.setMsg ("ARM_ERR: SFRM calibrator of CRA is not of a good type");
			return ARM_KO;
		}


        if(portfolioType == "Volatility")
        {
            newpf = calibrator->GetSigmaPortfolio() ? (ARM_StdPortfolio *)calibrator->GetSigmaPortfolio()->Clone() : NULL;
            if(!newpf)
            {
                	result.setMsg ("ARM_ERR: SFRM Calibrator contains a null portfolio of Volatility");
			    return ARM_KO;
            }

        }
        else if(portfolioType == "MeanReversion")
        {
            newpf = calibrator->GetMRSPortfolio() ? (ARM_StdPortfolio *)calibrator->GetMRSPortfolio()->Clone() : NULL; 
            if(!newpf)
            {
                	result.setMsg ("ARM_ERR: SFRM Calibrator contains a null portfolio of MeanReversion");
			    return ARM_KO;
            }
        }
        else if(portfolioType == "Beta")
        {
            newpf = calibrator->GetBetaPortfolio() ? (ARM_StdPortfolio *)calibrator->GetBetaPortfolio()->Clone() : NULL;  
            if(!newpf)
            {
                	result.setMsg ("ARM_ERR: SFRM Calibrator contains a null portfolio of Beta");
			    return ARM_KO;
            }
        }
        else
        {
			result.setMsg("ARM_ERR: CalibType is not of a good type(Volatility, MeanReversion, Beta)");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
            pfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newpf);

			if ( pfId == RET_KO )
			{
                delete newpf;
				newpf = NULL;

				result.setMsg ("ARM_ERR: Pb when inserting an object");				
				return ARM_KO;
			}

			result.setLong(pfId);

			return ARM_OK;
		}
		else
		{
		    pf = (ARM_StdPortfolio *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 1)
			{
				delete pf;
				pf = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newpf, objId);

				return ARM_OK;
			}
			else
			{
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

		ARM_RESULT();
	}
}




/// straight call to the formula ARM_BSCorrModel::SpreadOptionFormula
long ARMLOCAL_SpreadOptionFormula(
	double fwd1, 
	double fwd2, 
	double vol1, 
	double vol2, 
	double Correl, 
	double strike, 
	double optMat, 
	int optType, 
	int modelType, 
	double spreadVol,
	ARM_result& result)
{
	CCString msg ("");

	try
	{
		double price=ARM_BSCorrModel::SpreadOptionFormula(fwd1,fwd2,vol1,vol2,Correl,strike,optMat,optType,modelType,spreadVol);
		result.setDouble(price);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


// Build a sum option object
long ARMLOCAL_SumOption(
		const long& DateStripId,
		const double& strike,
		const long& capFloor,
		const double& coeff,
		const long& coeffCurveId,
		const double& payDate,
		const long& indexDayCount,
		ARM_result& result,
		long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_SumOpt* sumOpt	= NULL;

	ARM_Curve* coeffCurve = NULL;
	bool isCoeff = false;
    
	try
	{ 
		ARM_DateStrip* DateStrip= NULL;

		if( !GetObjectFromId( &DateStrip, DateStripId, ARM_DATESTRIP ) )
		{
			result.setMsg ("ARM_ERR: datestrip is not of a good type");
			return ARM_KO;
		}

		if (payDate != -1)
		{
			double julianPayDate = XLDateToJulian(payDate);
		}

		if(coeffCurveId != ARM_NULL_OBJECT)
        {
		    coeffCurve = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(coeffCurveId));

		    if (!coeffCurve)
		    {
			    result.setMsg ("ARM_ERR: coeffs curve  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			std::vector<double> times(1,0.0);
			std::vector<double> values(1,coeff);
            coeffCurve = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCoeff=true;
        }

        sumOpt = new ARM_SumOpt(
			*DateStrip, 
			payDate,
			strike,
			capFloor,
			indexDayCount,
			*coeffCurve);

		if (isCoeff)
			delete coeffCurve;

		/// assign object
		if( !assignObject( sumOpt, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (isCoeff)
			delete coeffCurve;

		delete sumOpt;

		x.DebugPrint();
		ARM_RESULT();
	}
}
// Build a smiled swaption object
long ARMLOCAL_SmiledSwaption(
		const long& BaseSwaptionId,
		const VECTOR<double>& Data,
		ARM_result& result,
		long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_SmiledSwaption* smiledSwaption	= NULL;
	ARM_Swaption* baseSwaption = NULL;
	int sizeVector = Data.size();
	std::vector<double> smileData(sizeVector);
	try
	{ 
		for(int i=0;i<sizeVector;++i)
		{
			smileData[i] = Data[i];
		}
		if( !GetObjectFromId( &baseSwaption, BaseSwaptionId, ARM_SWAPTION ) )
		{
			result.setMsg ("ARM_ERR: swaption is not of a good type");
			return ARM_KO;
		}

        smiledSwaption = new ARM_SmiledSwaption(
			*baseSwaption, 
			smileData);

		/// assign object
		if( !assignObject( smiledSwaption, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete baseSwaption;
		delete smiledSwaption;

		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_STRIPOPTION (	double asOfDate,
						    // Option data
							long underId,
							long optionType,
							long strikesId,
							long schedId,
							CCString PorS,
							long fxFixingsId,
							long leverageId,
							double leverageValue,
							// Results
							ARM_result& result,
							long objId)
{
	long optId;

    ARM_Security* sec					= NULL;
	ARM_DateStrip* schedule				= NULL;
	ARM_Object* tmpCurve				= NULL;
	ARM_Curve* strikes		            = NULL;
	ARM_Curve* fxFixings	            = NULL;
	ARM_Curve* leverage		            = NULL;
    ARM_StripOption* newStripOption		= NULL;
	ARM_StripOption* StripOption		= NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// Getting underlying
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		// Getting schedule
		schedule = (ARM_DateStrip *) LOCAL_PERSISTENT_OBJECTS->GetObject(schedId);

		// Getting strikes
		tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(strikesId);
		if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
		{
			strikes = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
		}
		else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
		{
			ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
			strikes = RefValueToCurve(*refValue, asOfDate);
		}
		else
		{
			result.setMsg ("ARM_ERR: Strikes curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
			return ARM_KO;
		}

		double vPorS = K_RCV;
		if (PorS == "PAY")
			vPorS = K_PAY;

		if (fxFixingsId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(fxFixingsId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				fxFixings = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				fxFixings = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: FX fixings curve is not of a good type");
				return ARM_KO;
			}
		}

		if (leverageId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				leverage = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				leverage = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: Leverage curve is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			// create a leverage constant to leverageValue
			leverage = new ARM_Curve ( std::vector<double>(1, -1), //lag
									   std::vector<double>(1, leverageValue), //value
									   new ARM_LinearInterpolatorCstExtrapol<double,double> );
		}

		newStripOption = new ARM_StripOption(ARM_Date(asOfDate),
											 *schedule, 
											 *strikes, 
											 sec, 
											 optionType, 
											 fxFixings, 
											 0.7, //correl
											 ARM_DEFAULT_COUNTRY, //payCcy
											 vPorS,
											 leverage);

		if (strikes)
		{
			delete strikes;
			strikes = NULL;
		}
		if (fxFixings)
		{
			delete fxFixings;
			fxFixings = NULL;
		}
		if (leverage)
		{
			delete leverage;
			leverage = NULL;
		}

		if (newStripOption == NULL)
		{
			result.setMsg ("ARM_ERR: Strip Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newStripOption);

			if (optId == RET_KO)
			{
				if (newStripOption)
					delete newStripOption;
				newStripOption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(optId);

			return ARM_OK;
		}
		else
		{
		    StripOption = (ARM_StripOption *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(StripOption, ARM_STRIPOPTION) == 1)
			{
				if (StripOption)
					delete StripOption;
				StripOption = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newStripOption, objId);

				return ARM_OK;
			}
			else
			{
				if (newStripOption)
					delete newStripOption;
				newStripOption = NULL;

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



long ARMLOCAL_STRIPDIGITALOPTION (	double asOfDate,
								    // Option data
									long underId,
									long optionType,
									long strikesId,
									long schedId,
									double correl,
									CCString PorS,
									long fxFixingsId,
									int callSpreadFlag,
									double epsilon,
									long payoffCurveId,
									long leverageId,
									double leverageValue,
									// Results
									ARM_result& result,
									long objId)
{
	long optId;

	ARM_Security* sec								= NULL;
	ARM_DateStrip* schedule							= NULL;
	ARM_Object* tmpCurve							= NULL;
	ARM_Curve* strikes								= NULL;
	ARM_Curve* fxFixings							= NULL;
	ARM_Curve* leverage								= NULL;
	ARM_Curve* payoff								= NULL;
    ARM_StripDigitalOption* newStripDigitalOption	= NULL;
	ARM_StripDigitalOption* StripDigitalOption		= NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// Getting underlying
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		// Getting schedule
		schedule = (ARM_DateStrip *) LOCAL_PERSISTENT_OBJECTS->GetObject(schedId);

		// Getting strikes
		tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(strikesId);
		if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
		{
			strikes = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
		}
		else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
		{
			ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
			strikes = RefValueToCurve(*refValue, asOfDate);
		}
		else
		{
			result.setMsg ("ARM_ERR: Strikes curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
			return ARM_KO;
		}

		double vPorS = K_RCV;
		if (PorS == "PAY")
			vPorS = K_PAY;

		if (fxFixingsId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(fxFixingsId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				fxFixings = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				fxFixings = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: FX fixings curve is not of a good type");
				return ARM_KO;
			}
		}

		if (payoffCurveId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(payoffCurveId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				payoff = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				payoff = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: Payoff curve is not of a good type");
				return ARM_KO;
			}
			payoff->GetOrdinates() *= 100.0;
		}

		if (leverageId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				leverage = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				leverage = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: Leverage curve is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			// create a leverage constant to leverageValue
			leverage = new ARM_Curve ( std::vector<double>(1, -1), //lag
									   std::vector<double>(1, leverageValue), //value
									   new ARM_LinearInterpolatorCstExtrapol<double,double> );
		}

		newStripDigitalOption = new ARM_StripDigitalOption( ARM_Date(asOfDate),
															*schedule,
															*strikes, 
															sec, 
															optionType, 
															correl,
															ARM_DEFAULT_COUNTRY,
															vPorS,
															callSpreadFlag,
															epsilon,
															fxFixings,
															payoff,
															leverage);

		if (strikes)
		{
			delete strikes;
			strikes = NULL;
		}
		if (fxFixings)
		{
			delete fxFixings;
			fxFixings = NULL;
		}
		if (leverage)
		{
			delete leverage;
			leverage = NULL;
		}
		if (payoff)
		{
			delete payoff;
			payoff = NULL;
		}

		if (newStripDigitalOption == NULL)
		{
			result.setMsg ("ARM_ERR: Strip Digital Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newStripDigitalOption);

			if (optId == RET_KO)
			{
				if (newStripDigitalOption)
					delete newStripDigitalOption;
				newStripDigitalOption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(optId);
			return ARM_OK;
		}
		else
		{
		    StripDigitalOption = (ARM_StripDigitalOption *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(StripDigitalOption, ARM_STRIPDIGITALOPTION) == 1)
			{
				if (StripDigitalOption)
					delete StripDigitalOption;
				StripDigitalOption = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newStripDigitalOption, objId);
				return ARM_OK;
			}
			else
			{
				if (newStripDigitalOption)
					delete newStripDigitalOption;
				newStripDigitalOption = NULL;

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



long ARMLOCAL_FxOptionStrip( double asOfDate,
							 long C_underlyingId,
							 long C_strikesCurveId,
							 long optionType,
							 double C_startDate,
							 double C_endDate,
							 long C_notionalId,
							 CCString C_paymentCcy,
							 long C_resetFreq,
							 long C_dayCount,
							 CCString C_resetCalendar,
							 long C_fwdRule,
							 long C_intRule,
							 long C_stubRule,
							 long C_resetGap,
							 long C_payFreq, 
							 long C_payGap, 
							 CCString C_payCalendar, 
							 long C_resetTiming, 
							 long C_payTiming,
							 CCString C_PorS,
							 long C_fxFixingsId,
							 bool isDigital,
							 int C_callSpreadFlag,
							 double C_epsilon,
							 long C_payoffCurveId,
							 long C_leverageId,
							 double C_leverageValue,
							 ARM_result& result,
							 long objId)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	long vStripId;
	ARM_Security*	vUnderlying = NULL;
	ARM_Object*		tmpCurve = NULL;
	ARM_Curve*		vStrikes = NULL;
	ARM_Curve*		vNotional = NULL;
	ARM_Curve*		vFxFixings = NULL;
	ARM_Curve*		vPayoff = NULL;
	ARM_Curve*		vLeverage = NULL;
	int i=0;
	ARM_StripOption*	vStripOption = NULL;
	ARM_StripOption*	vNewStripOption = NULL;

	try
	{
		vUnderlying = (ARM_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(C_underlyingId)->Clone();
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vUnderlying, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(C_strikesCurveId);
		if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
		{
			vStrikes = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
		}
		else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
		{
			ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
			vStrikes = RefValueToCurve(*refValue, asOfDate);
		}
		else
		{
			result.setMsg ("ARM_ERR: Strikes curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
			return ARM_KO;
		}

		tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(C_notionalId);
		if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
		{
			vNotional = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
		}
		else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
		{
			ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
			vNotional = RefValueToCurve(*refValue, asOfDate);
		}
		else
		{
			result.setMsg ("ARM_ERR: Notional curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
			return ARM_KO;
		}

		double PorS = K_RCV;
		if (C_PorS == "PAY")
			PorS = K_PAY;

		if (C_fxFixingsId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(C_fxFixingsId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				vFxFixings = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				vFxFixings = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: FX fixings curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
				return ARM_KO;
			}
		}

		if (C_leverageId != ARM_NULL_OBJECT)
		{
			tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(C_leverageId);
			if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
			{
				vLeverage = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
			}
			else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
			{
				ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
				vLeverage = RefValueToCurve(*refValue, asOfDate);
			}
			else
			{
				result.setMsg ("ARM_ERR: Leverage curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
				return ARM_KO;
			}
		}
		else
		{
			// create a leverage constant to leverageValue
			vLeverage = new ARM_Curve ( std::vector<double>(1, -1), //lag
									   std::vector<double>(1, C_leverageValue), //value
									   new ARM_LinearInterpolatorCstExtrapol<double,double> );
		}

		if (!isDigital)
		{
			vNewStripOption = new ARM_StripOption(ARM_Date(asOfDate),
												  (ARM_Forex*)vUnderlying, 
												  *vStrikes, 
												  optionType, 
												  ConvertToARMDATE(C_startDate), 
												  ConvertToARMDATE(C_endDate), 
												  C_resetFreq, C_dayCount, 
												  *vNotional, 
												  C_resetCalendar, 
												  C_fwdRule, 
												  C_intRule, 
												  C_stubRule, 
												  C_resetGap,
												  C_payFreq, 
												  C_payGap, 
												  C_payCalendar, 
												  C_resetTiming, 
												  C_payTiming,
												  C_paymentCcy,
												  PorS,
												  vFxFixings,
												  vLeverage);
		}
		else
		{
			if (C_payoffCurveId != ARM_NULL_OBJECT)
			{
				tmpCurve = LOCAL_PERSISTENT_OBJECTS->GetObject(C_payoffCurveId);
				if (tmpCurve->GetName() == ARM_GENERIC_CURVE)
				{
					vPayoff = dynamic_cast<ARM_Curve*>(tmpCurve->Clone());
				}
				else if (tmpCurve->GetName() == ARM_REFERENCE_VALUE)
				{
					ARM_ReferenceValue* refValue = dynamic_cast<ARM_ReferenceValue*>(tmpCurve);
					vPayoff = RefValueToCurve(*refValue, asOfDate);
				}
				else
				{
					result.setMsg ("ARM_ERR: Payoff curve is not of a good type (ARM_Curve or ARM_ReferenceValue)");
					return ARM_KO;
				}
				vPayoff->GetOrdinates() *= 100.0;
			}

			vNewStripOption = new ARM_StripDigitalOption( ARM_Date(asOfDate),
														  (ARM_Forex*)vUnderlying, 
														  *vStrikes, 
														  optionType, 
														  ConvertToARMDATE(C_startDate), 
														  ConvertToARMDATE(C_endDate), 
														  C_resetFreq, C_dayCount, 
														  *vNotional, 
														  C_resetCalendar, 
														  C_fwdRule, 
														  C_intRule, 
														  C_stubRule, 
														  C_resetGap,
														  C_payFreq, 
														  C_payGap, 
														  C_payCalendar, 
														  C_resetTiming, 
														  C_payTiming,
														  C_paymentCcy,
														  PorS,
														  C_callSpreadFlag,
														  C_epsilon,
														  vFxFixings,
														  vPayoff,
														  vLeverage);
		}

		if (vNotional)
		{
			delete vNotional;
			vNotional = NULL;
		}
		if (vStrikes)
		{
			delete vStrikes;
			vStrikes = NULL;
		}
		if (vFxFixings)
		{
			delete vFxFixings;
			vFxFixings = NULL;
		}
		if (vLeverage)
		{
			delete vLeverage;
			vLeverage = NULL;
		}
		if (vPayoff)
		{
			delete vPayoff;
			vPayoff = NULL;
		}

		if(vNewStripOption == NULL)
		{
			result.setMsg ("ARM_ERR: FX Option Strip Object is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			vStripId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vNewStripOption);

			if (vStripId == RET_KO)
			{
				delete vNewStripOption;
				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(vStripId);
			return ARM_OK;
		}
		else
		{
			vStripOption = dynamic_cast<ARM_StripOption*>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
			
			if ( (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vStripOption, ARM_STRIPOPTION) == 1)
				|| (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vStripOption, ARM_STRIPDIGITALOPTION) == 1) )
			{
				delete vStripOption;				
				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vNewStripOption, objId);
				return ARM_OK;
			}
			else
			{
				delete vNewStripOption;
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