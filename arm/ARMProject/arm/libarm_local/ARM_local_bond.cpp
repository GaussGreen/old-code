#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <inst\bond.h>
#include <inst\bondtec.h>
#include <inst\portfolio.h>

#include <util\fromto.h>


long ARMLOCAL_bond (double issueDate, double maturityDate,
					  double firstCouponDate, double couponRate,
					  double redemptionPrice, long periodicity,
					  long dayCount, long settleGap, long couponDateFlag, long ccyId, 
					  ARM_result& result, long objId)
{
	long bondId;

    ARM_Bond* createdBond=NULL;
    ARM_Bond* bond=NULL;
	ARM_Currency* pccy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char dateEmiss[11];
	char dateEch[11];
	char datePremCoupon[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(issueDate,dateEmiss);
		Local_XLDATE2ARMDATE(maturityDate,dateEch);
		Local_XLDATE2ARMDATE(firstCouponDate,datePremCoupon);

		if (ccyId == -1) 
		  pccy = ARM_DEFAULT_CURRENCY;
		else
		{	
		  pccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);	
		}

        createdBond = new ARM_Bond((ARM_Date) dateEmiss, (ARM_Date) dateEch,
                                   (ARM_Date) datePremCoupon, couponRate,
                                   redemptionPrice, periodicity, dayCount, settleGap,
                                   couponDateFlag, pccy);

		if (createdBond == NULL)
		{
			result.setMsg ("ARM_ERR: Bond is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			bondId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBond);

			if (bondId == RET_KO)
			{
				if (createdBond)
					delete createdBond;
				createdBond = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(bondId);

			return ARM_OK;
		}
		else
		{

			bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);


			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 1)
			{
				if (bond)
				{
					delete bond;
					bond = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBond, objId);
				
				return ARM_OK;
			}
			else
			{
				if (createdBond)
					delete createdBond;
				createdBond = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdBond)
			delete createdBond;
		createdBond = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_RiskyBond (double issueDate, double maturityDate,
						double firstCouponDate, double couponRate,
						double redemptionPrice, long periodicity,
						long dayCount, long settleGap, long couponDateFlag, long ccyId, 
						double sRepo, double ssl, double recoveryRate,
						ARM_result& result, long objId)
{
	long bondId;

    ARM_Bond* createdRiskyBond=NULL;
    ARM_Bond* riskyBond=NULL;
	ARM_Currency* pccy = NULL;
	bool IsRiskyBond = true;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char dateEmiss[11];
	char dateEch[11];
	char datePremCoupon[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(issueDate,dateEmiss);
		Local_XLDATE2ARMDATE(maturityDate,dateEch);
		Local_XLDATE2ARMDATE(firstCouponDate,datePremCoupon);

		if (ccyId == -1) 
		  pccy = ARM_DEFAULT_CURRENCY;
		else
		{	
		  pccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);	
		}
		

		if (sRepo == -1 && ssl== -1 && recoveryRate == -1)
		{
			IsRiskyBond = false;	
			createdRiskyBond = new ARM_Bond((ARM_Date) dateEmiss, (ARM_Date) dateEch,
											(ARM_Date) datePremCoupon, couponRate,
											redemptionPrice, periodicity, dayCount, settleGap,
											couponDateFlag, pccy);

		}
		else
		{
			
			if (sRepo == -1)
			{
				result.setMsg ("ARM_ERR: srepo must have a value");
				return ARM_KO;
			}
			
			if (ssl == -1)
			{
				result.setMsg ("ARM_ERR: ssl must have a value");
				return ARM_KO;
			}

			if (recoveryRate == -1)
			{
				result.setMsg ("ARM_ERR: Recovery rate must have a value");
				return ARM_KO;
			}
			
			createdRiskyBond = new ARM_RiskyBond((ARM_Date) dateEmiss, (ARM_Date) dateEch,
												(ARM_Date) datePremCoupon, couponRate,
												redemptionPrice, periodicity, dayCount, settleGap,
												couponDateFlag, pccy, sRepo, ssl, recoveryRate);

		}

		if (createdRiskyBond == NULL)
		{
			result.setMsg ("ARM_ERR: Bond is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			bondId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRiskyBond);

			if (bondId == RET_KO)
			{
				if (createdRiskyBond)
					delete createdRiskyBond;
				createdRiskyBond = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(bondId);

			return ARM_OK;
		}
		else
		{
			if(IsRiskyBond)
			{
				riskyBond = (ARM_RiskyBond *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			}
			else
			{
				riskyBond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			}

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(riskyBond, ARM_RISKYBOND) == 1 || 
				LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(riskyBond, ARM_BOND) == 1)
			{
				if (riskyBond)
				{
					delete riskyBond;
					riskyBond = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRiskyBond, objId);
				
				return ARM_OK;
			}
			else
			{
				if (createdRiskyBond)
					delete createdRiskyBond;
				createdRiskyBond = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdRiskyBond)
			delete createdRiskyBond;
		createdRiskyBond = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_RiskyBondWithCF (double asOfDate, double redemptionPrice, long periodicity,
							   long dayCount, long settleGap, long couponDateFlag, long ccyId, 
							   VECTOR<double>& yearTerms, VECTOR<double>& cashFlows, 
							   double sRepo, double ssl, double recoveryRate,
							   ARM_result& result, long objId)
{
	long bondId;

    ARM_RiskyBondWithCF* createdRiskyBondWithCF = NULL;
    ARM_RiskyBondWithCF* riskyBond = NULL;
	ARM_Currency* pccy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char asOf[11];
	char tmpDate[11];
	vector<ARM_Date> yearTermsVect;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(asOfDate, asOf);

		
		for (int j = 0; j <yearTerms.size(); j++) 
		{
			Local_XLDATE2ARMDATE(yearTerms[j],tmpDate);
			yearTermsVect.push_back(tmpDate);
		}

		if (ccyId == -1) 
		  pccy = ARM_DEFAULT_CURRENCY;
		else
		{	
		  pccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);	
		}
		
			
		if (sRepo == -1)
		{
			result.setMsg ("ARM_ERR: srepo must have a value");
			return ARM_KO;
		}
		
		if (ssl == -1)
		{
			result.setMsg ("ARM_ERR: ssl must have a value");
			return ARM_KO;
		}

		if (recoveryRate == -1)
		{
			result.setMsg ("ARM_ERR: Recovery rate must have a value");
			return ARM_KO;
		}

		if ( yearTerms.size() > 0 && cashFlows.size() > 0)
		{
		   ARM_Vector* cashFlowsVect = CreateARMVectorFromVECTOR(cashFlows);

		
		   createdRiskyBondWithCF = new ARM_RiskyBondWithCF((ARM_Date)asOf, redemptionPrice, periodicity, dayCount, 
															&yearTermsVect, ARM_Vector(cashFlowsVect), settleGap,
															couponDateFlag, pccy,
															sRepo, ssl, recoveryRate);
		   delete cashFlowsVect;
		}
		else 
		{
			result.setMsg ("ARM_ERR: Year Terms and Cash Flows : array expected");
			return ARM_KO;
		}
		

		if (createdRiskyBondWithCF == NULL)
		{
			result.setMsg ("ARM_ERR: Bond is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			bondId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) createdRiskyBondWithCF);

			if (bondId == RET_KO)
			{
				if (createdRiskyBondWithCF)
					delete createdRiskyBondWithCF;
				createdRiskyBondWithCF = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(bondId);

			return ARM_OK;
		}
		else
		{

			riskyBond = (ARM_RiskyBondWithCF *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);


			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(riskyBond, ARM_RISKYBONDWITHCF) == 1 || 
				LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(riskyBond, ARM_BOND) == 1)
			{
				if (riskyBond)
				{
					delete riskyBond;
					riskyBond = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRiskyBondWithCF, objId);
				
				return ARM_OK;
			}
			else
			{
				if (createdRiskyBondWithCF)
					delete createdRiskyBondWithCF;
				createdRiskyBondWithCF = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdRiskyBondWithCF)
			delete createdRiskyBondWithCF;
		createdRiskyBondWithCF = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_bondTEC (double issueDate, double maturityDate,
					  double firstCouponDate, double couponRate,
					  double redemptionPrice, long periodicity,
					  long dayCount, long settleGap, long couponDateFlag, long ccyId, double tec, long pfTECId, long modTECid ,
					  ARM_result& result, long objId)
{
	long bondId;

    ARM_Bond* createdBond=NULL;
    ARM_Bond* bond=NULL;
	ARM_BondTEC* bondtec=NULL;
	ARM_Currency* pccy = NULL;
	ARM_Portfolio* ppftecId = NULL;
	ARM_Model* pmodTECid = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char dateEmiss[11];
	char dateEch[11];
	char datePremCoupon[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(issueDate,dateEmiss);
		Local_XLDATE2ARMDATE(maturityDate,dateEch);
		Local_XLDATE2ARMDATE(firstCouponDate,datePremCoupon);

		if (ccyId == -1) 
		  pccy = ARM_DEFAULT_CURRENCY;
		else
		{	
		  pccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);	
		}


	    ppftecId = (ARM_Portfolio*) LOCAL_PERSISTENT_OBJECTS->GetObject(pfTECId);
		pmodTECid = (ARM_Model*) LOCAL_PERSISTENT_OBJECTS->GetObject(modTECid);

		if ((tec) && (ppftecId))  // Pour la TEC le portfolio doit être un portfolio de bonds
		{

		if (ppftecId->GetAsset(0)->GetName()!=ARM_BOND)
			{
			result.setMsg ("ARM_ERR: TEC cover portfolio is not available for others securities !");
			return ARM_KO;
			}	
		}

        createdBond = new ARM_BondTEC((ARM_Date) dateEmiss, (ARM_Date) dateEch,
                                   (ARM_Date) datePremCoupon, couponRate,
                                   redemptionPrice, periodicity, dayCount, settleGap,
                                   couponDateFlag, pccy ,tec, ppftecId, pmodTECid);


		if (createdBond == NULL)
		{
			result.setMsg ("ARM_ERR: Bond is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			bondId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBond);

			if (bondId == RET_KO)
			{
				if (createdBond)
					delete createdBond;
				createdBond = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(bondId);

			return ARM_OK;
		}
		else
		{
			bondtec = (ARM_BondTEC *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bondtec, ARM_BONDTEC) == 1) && (tec))
			{
				if (bondtec)
				{
					delete bondtec;
					bondtec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdBond, objId);
				
				return ARM_OK;
			}
			else
			{
				if (createdBond)
					delete createdBond;
				createdBond = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdBond)
			delete createdBond;
		createdBond = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_YTOPRICE (long bondId, double settlement, double yield, ARM_result& result)
{
	double price;
	ARM_Bond* bond=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetSettlement((ARM_Date) sSettlement);
		bond->SetYield(yield); 

		price = bond->YieldToPrice((ARM_Date) sSettlement, yield);

		result.setDouble(price);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_PTOYIELD (long bondId, double settlement, double price, ARM_result& result)
{
	double yield;
	ARM_Bond* bond=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BONDTEC) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetSettlement((ARM_Date) sSettlement);
		bond->SetPrice(price); 

		yield = bond->PriceToYield((ARM_Date) sSettlement, price);

		result.setDouble(yield);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_BDFAPRICE (long bondId, double settlement, double actuPrice, double forwardDate, 
					       double repoRate, ARM_result& result)
{
    ARM_Bond* bond=NULL;
	double price;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];
	char sFwdDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);
		Local_XLDATE2ARMDATE(forwardDate,sFwdDate);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetModel(NULL);

		bond->SetSettlement((ARM_Date) sSettlement);

		bond->SetPrice(actuPrice);

		bond->SetForwardDate((ARM_Date) sFwdDate);

		bond->SetRepoRate(repoRate);

		price = bond->ComputeForwardPrice();

		result.setDouble(price);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



long ARMLOCAL_BDREPORATE (long bondId, double settlement, 
					        double actuPrice, double forwardDate,
					        double forwardPrice,
					        ARM_result& result)
{
    ARM_Bond* bond=NULL;
	double repoRate;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];
	char sFwdDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);
		Local_XLDATE2ARMDATE(forwardDate,sFwdDate);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetModel(NULL);

		bond->SetSettlement((ARM_Date) sSettlement);

		bond->SetPrice(actuPrice);

		bond->SetForwardDate((ARM_Date) sFwdDate);

		bond->SetForwardPrice(forwardPrice);

		repoRate = bond->ComputeRepoRate();

		result.setDouble(repoRate);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_YTODURATION (long bondId,
						   double settlement,
						   double actuRate,
						   long flagCpn,
						   ARM_result& result)
{
	double duration;
	
	ARM_Bond* bond=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetModel(NULL);

		duration = bond->YieldToDuration((ARM_Date) sSettlement, actuRate, flagCpn);

		result.setDouble(duration);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_YTOCONVEXITY (long bondId,
							double settlement,
							double actuRate,
							ARM_result& result)
{
	double convexity;
	
	ARM_Bond* bond=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sSettlement[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		bond = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) &&
			(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_RISKYBOND) == 0))
		{
			result.setMsg ("ARM_ERR: bond is not of a good type");
			return ARM_KO;
		}

		bond->SetModel(NULL);

		convexity = bond->YieldToConvexity((ARM_Date) sSettlement, actuRate);

		result.setDouble(convexity);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
