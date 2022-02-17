#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <ARM\libarm\ARM_result.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_ccy.h>
#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_assetswap.h>
#include <ARM\libarm_local\ARM_local_swap.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_zccurve.h>
#include <ARM\libarm_local\ARM_local_irindex.h>

#include <ccy\currency.h>
#include <util\refvalue.h>
#include <inst\bond.h>
#include <inst\irindex.h>

#include <ARM\armdate\ARM_date.h>

#include <util\ournag.h>
#include <nag_stdlib.h>

#include <nage04.h>

#include <math.h>



/************* Fin des Include *************/

#define CURRENCY_USD	"USD"

#define c_a				0.00001
#define c_s				0.00001
#define c_cor			0.0001

static long g_index1Id = ARM_NULL_OBJECT_ID;
static long g_index2Id = ARM_NULL_OBJECT_ID;
static long g_model1Id = ARM_NULL_OBJECT_ID;
static long g_model2Id = ARM_NULL_OBJECT_ID;

static long leg1Id = ARM_NULL_OBJECT_ID;
static long leg2Id = ARM_NULL_OBJECT_ID;

static double g_price1, g_price2;
static double g_minimum;
static int g_nbiter;
static double sensi_std, sensi_GBP;
static double startXNL1, endXNL1, startXNL2, endXNL2;

static long LASMC_ret = ARM_KO;
static long CDP_ret = ARM_KO;

static ARM_result* g_result = new ARM_result ();


/**************** Fin des define et static ************/

struct ARMLOCAL_CDP_arg
{
	double margin1;
	double delivery;
	double maturity;
	CCString* ccy1;
	long ccy1Id;
	CCString* ccy2;
	long ccy2Id;
};
typedef struct ARMLOCAL_CDP_arg ARMLOCAL_CDP_arg;



double ASP_getStartXNL1 ()
{
	return startXNL1;
}
double ASP_getEndXNL1 ()
{
	return endXNL1;
}
double ASP_getStartXNL2 ()
{
	return startXNL2;
}
double ASP_getEndXNL2 ()
{
	return endXNL2;
}
double ASP_getPrice1 ()
{
	return g_price1;
}
double ASP_getPrice2 ()
{
	return g_price2;
}
double ASP_getMinimum ()
{
	return g_minimum;
}
double ASP_getNbIter ()
{
	return g_nbiter;
}


static long Local_computePrice1 (ARMLOCAL_CDP_arg* arg, double margin1, double* price1)
{
	ARM_result C_result;
	double liborLegPrice = 0.0;

	double delivery = arg->delivery;
	double maturity = arg->maturity;
	
	*price1 = 0.0;
	
	// if(ARM_SWAPLEG (g_index1Id, delivery, maturity, K_RCV, margin1 / 100.0, g_dccy1Id, -1, C_result, leg1Id) == ARM_OK)
	if(ARMLOCAL_SWAPLEG (g_index1Id,
						 delivery,
						 maturity,
						 K_RCV,
						 0L,
						 margin1 / 100.0,
						 //arg->ccy1Id,
						 false,
						 *(arg->ccy1),
						 -1,
						 10000,
						 "NULL",
						 "NULL",
						 1,
						 K_NX_NONE,
						 K_SHORTSTART,
						 -1.,
						 K_YES,
						 -1,
						 C_result,
						 leg1Id) == ARM_OK)
	{
		if(leg1Id == ARM_NULL_OBJECT_ID)
		{
			leg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "Local_computePrice1 (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (leg1Id,
						   g_model1Id,
						   C_result) == ARM_OK)
	{
		liborLegPrice = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "Local_computePrice1 (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	*price1 = liborLegPrice + startXNL1 + endXNL1;

// Rajout correction base sensibilit� si devise = GBP

	if(*(arg->ccy1) == "GBP")
    {
		*price1 = *price1 + margin1 * (sensi_GBP - sensi_std);
	}

	return ARM_OK;
}

static long Local_computePrice2 (ARMLOCAL_CDP_arg* arg, double margin2, double* price2)
{
	ARM_result C_result;
	double liborLegPrice = 0.0;

	double delivery = arg->delivery;
	double maturity = arg->maturity;

	*price2 = 0.0;

	// if(ARM_SWAPLEG (g_index2Id, delivery, maturity, K_RCV, margin2 / 100.0, g_dccy2Id, -1, C_result, leg2Id) == ARM_OK)
	if(ARMLOCAL_SWAPLEG (g_index2Id,
						 delivery,
						 maturity,
						 K_RCV,
						 0L,
						 margin2 / 100.0,
						 //arg->ccy2Id,
						 false,
						 *(arg->ccy2),
						 -1,
						 10000,
						 "NULL",
						 "NULL",
						 1,
						 K_NX_NONE,
						 K_SHORTSTART,
						 -1.,
						 K_YES,
						 -1,
						 C_result,
						 leg2Id) == ARM_OK)
	{
		if(leg2Id == ARM_NULL_OBJECT_ID)
		{
			leg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "Local_computePrice2 (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (leg2Id,
						   g_model2Id,
						   C_result) == ARM_OK)
	{
		liborLegPrice = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "Local_computePrice2 (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	*price2 = liborLegPrice + startXNL2 + endXNL2;

// Rajout correction base sensibilit� si devise = GBP

	if(*(arg->ccy2) == "GBP")
    {	
		*price2 = *price2 + margin2 * (sensi_GBP - sensi_std);
	}

	return ARM_OK;
}


static void __stdcall Local_computeDiffPrice (double margin2, double* liborASW, Nag_Comm* comm)
{
	double l_price1, l_price2;

	ARMLOCAL_CDP_arg* arg = (ARMLOCAL_CDP_arg*)comm->p;
	
	if(Local_computePrice1 (arg, arg->margin1, &l_price1) == ARM_KO)
	//if(computePrice1 (arg->delivery, arg->maturity, arg->margin1, &l_price1) == ARM_KO)
	{
		comm->p = (void*)g_result;
		*liborASW = 0.0;
		CDP_ret = ARM_KO;
		return;
	}
	if(Local_computePrice2 (arg, margin2, &l_price2) == ARM_KO)
	//if(computePrice2 (arg->delivery, arg->maturity, margin2, &l_price2) == ARM_KO)
	{
		comm->p = (void*)g_result;
		*liborASW = 0.0;
		CDP_ret = ARM_KO;
		return;
	}

	*liborASW = fabs (l_price1 - l_price2) * 1000000.0;
	CDP_ret = ARM_OK;
	return;
}

static void Local_destroySwapLegs ()
{
	ARM_result C_result;

	ARMLOCAL_FreeObject (leg1Id, C_result);
	leg1Id = ARM_NULL_OBJECT_ID;
	ARMLOCAL_FreeObject (leg2Id, C_result);
	leg2Id = ARM_NULL_OBJECT_ID;
}





/**********************************************************

***********************************************************/
long ARMLOCAL_fixing (double asOfDate,
					  double startDate,
					  double endDate,
					  long frequencyId,
					  long forwardCurveId,
					  const CCString& ccy,
					  double fixing,
					  double spread,
					  ARM_result& result)
{
	ARM_result C_result;

	double dfs = 0.0;
	double dfn = 0.0;
	double dfnAdj = 0.0;
	double value = 0.0;
	double daycount = 360.0;
	double forward = 0.0;
	double temp = 0.0;

	if(ccy == "GBP")
	{
		daycount = 365.0;
	}

	double nextDate, prevDate, nextDateAdj, resetDate;
		
	if(ARMLOCAL_ARM_NextCpnDate (startDate,
								 endDate,
								 frequencyId,
								 K_MOD_FOLLOWING,
								 ccy,
								 C_result) == ARM_OK)
	{
		nextDate = C_result.getDouble ();
	}
	
	if(ARMLOCAL_ARM_NextCpnDate (startDate,
								 endDate,
								 frequencyId,
								 K_MOD_FOLLOWING,
								 ccy,
								 C_result,
								 K_ADJUSTED) == ARM_OK)
	{
		nextDateAdj = C_result.getDouble ();
	}

	if(ARMLOCAL_ARM_PrevCpnDate (startDate,
								 endDate,
								 frequencyId,
								 K_MOD_FOLLOWING,
								 ccy,
								 C_result) == ARM_OK)
	{
		prevDate = C_result.getDouble ();
	}

	if(ARMLOCAL_NextBusinessDay (prevDate,
								 ccy,
								 -2,
								 C_result) == ARM_OK)
	{
		resetDate = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	if(ARMLOCAL_DiscountPrice (forwardCurveId,
							   (startDate - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfs = result.getDouble ();

	if(ARMLOCAL_DiscountPrice (forwardCurveId,
							   (nextDate - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfn = result.getDouble ();

	if(ARMLOCAL_DiscountPrice (forwardCurveId,
							   (nextDateAdj - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfnAdj = result.getDouble ();

	forward = ((dfs / dfn) - 1.0) * (daycount * 100.0) / (nextDate - startDate);

	// value = (fixing - spread - ((dfs / dfn) - 1.0) * daycount / (nextDate - startDate)) * (nextDate - startDate) / daycount * (dfn/dfs);

	if (resetDate <= asOfDate)
		temp = (fixing - forward) * (nextDate - startDate) / daycount * (dfnAdj / dfs);

	value = temp + ((fixing + spread/100.) * ((startDate - prevDate) / daycount) * ((dfnAdj / dfs) - 1.0));

	MSG_printf_message (MSG_TRACE, "FRN Fixing adjustment : (nextDate = %lf) ; (prevDate = %lf) ; (resetDate = %lf)", nextDate,prevDate,resetDate);
	MSG_printf_message (MSG_TRACE, "FRN Fixing adjustment : (Fwd = %lf) ; (dfs = %lf) ; (dfn = %lf) - (first adj : Fixing fwd = %lf) ; (2nd adj : Accrued = %lf) ; (Total adj = %lf)", forward, dfs, dfn, temp, value - temp, value);

	result.setDouble (value);
	result.setRetCode (ARM_OK);

	return ARM_OK;
}



/**************************************************************

***************************************************************/
long ARMLOCAL_FRNPriceNew (double asOfDate,
						   double delivery,
						   double maturity,
						   const CCString& ccy1,
						   long liborType1Id,
						   long forwardCurve1Id,
						   long discountCurve1Id,
						   double facialMargin,
						   double valoMargin,
						   long frequencyId,
						   const CCString& ccy2,
						   long liborType2Id,
						   long forwardCurve2Id,
						   long discountCurve2Id,
						   long amortizationId,
						   long frequencyId2,
						   double fixing,
						   double spread,
						   long solve,
						   ARM_result& result)
{
	ARM_result C_result;
	
	double sensibility = 0.0;
	double sensibility2 = 0.0;
	double price = 0.0;
	double ccy1margin = 0.0;

	long retCode = ARM_KO;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}

	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
		return ARM_KO;
	}

    int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	if(ARMLOCAL_CptBPVNew (asOfDate,
						   delivery,
						   maturity,
						   forwardCurve1Id,
						   frequencyId,
						   inputDayCount,
						   ccy1Id,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Sensibility : (%lf)", sensibility);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(sensibility == 0.0)
	{
		result.setRetCode (ARM_KO);
		result.setMsg ("ARM_ERR: sensibility is nul");
		return ARM_KO;
	}

	ccy1margin = valoMargin;

/*	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		retCode = ARM_OK;
		ccy1margin = valoMargin;
		
		if (frequencyId != frequencyId2)
		{
			if(ARMLOCAL_CptBPVNew (asOfDate, delivery, maturity, forwardCurve1Id, frequencyId2, inputDayCount, ccy1Id, amortizationId, C_result) == ARM_OK)
			{
				sensibility2 = C_result.getDouble ();
				MSG_printf_message (MSG_TRACE, "FRN Sensibility2 : (%lf)", sensibility2);
			}
			else
			{
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if(sensibility2 == 0.0)
			{
				result.setRetCode (ARM_KO);
				result.setMsg ("ARM_ERR: sensibility2 is nul");
				return ARM_KO;
			}

			ccy1margin *= sensibility2 / sensibility;
		}
//		price = facialMargin - ccy1margin;
	}
	else
	{
*/	if ((ccy1 != ccy2) && (!(ccy2 == "NONE")))
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
			return ARM_KO;
		}

		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 maturity,
										 valoMargin,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 amortizationId,
										 solve,
										 C_result);
		
		if(retCode == ARM_OK)
		{
			ccy1margin = C_result.getDouble ();
			//price = facialMargin - ccy1margin;
			MSG_printf_message (MSG_TRACE, "FRN Basis Conversion Margin : valoMargin(2) = (%lf) ; implyMargin(1) = (%lf) ; method = (%ld)", valoMargin, C_result.getDouble (), solve);
		}
		else
		{
			result.setMsg (C_result.getMsg ());
		}
	}

	double dfs = 1.0;

	if(ARMLOCAL_DiscountPrice (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfs = result.getDouble ();
    
/*
	if(ARMLOCAL_CptBPVNew (asOfDate, delivery, maturity, forwardCurve1Id, frequencyId2, inputDayCount, ccy1Id, C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
*/
//	modif 10/10/2001
//	price = 100.0 + price * sensibility / dfs;

	long liborfltId;
	long liborswapId;

	// Modif du 28/11/01
//	if (ARMLOCAL_LIBORLEG(delivery,maturity,liborType1Id,K_RCV,facialMargin/100.,frequencyId,frequencyId,K_ADVANCE,K_ARREARS,ccy1Id,K_ADJUSTED,10,"NULL","NULL",1,C_result) == ARM_OK)
	if (ARMLOCAL_LIBORLEG(delivery,
						  maturity,
						  liborType1Id,
						  K_RCV,
						  0,
						  facialMargin/100.,
						  frequencyId,
						  frequencyId,
						  K_ADVANCE,
						  K_ARREARS,
						  //ccy1Id,
						  false,
						  ccy1,
						  K_UNADJUSTED,
						  10,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.,
						  K_YES,-1,
						  C_result) == ARM_OK)
	{
		liborfltId = C_result.getLong ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (liborfltId,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_FRNPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if (ARMLOCAL_LIBORLEG(delivery,
						  maturity,
						  liborType2Id,
						  K_RCV,
						  0,
						  ccy1margin/100.,
						  frequencyId2,
						  frequencyId2,
						  K_ADVANCE,
						  K_ARREARS,
						  //ccy1Id,
						  false,
						  ccy1,
						  K_ADJUSTED,
						  10,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.,
						  K_YES,-1,
						  C_result) == ARM_OK)
	{
		liborswapId = C_result.getLong ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (liborswapId,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_FRNPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	long ycModId;

	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result) == ARM_OK)
	{
		ycModId = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	double price1, price2;

	if(ARMLOCAL_ARM_Price (liborfltId,
						   ycModId,
						   C_result) == ARM_OK)
	{
		price1 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Price Details : price float leg = %lf\n", price1);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (liborswapId,
						   ycModId,
						   C_result) == ARM_OK)
	{
		price2 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Price Details : price swap leg = %lf\n", price2);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	ARMLOCAL_FreeObject (liborfltId, C_result);
	ARMLOCAL_FreeObject (liborswapId, C_result);
	ARMLOCAL_FreeObject (ycModId, C_result);

	price = 100. + (price1 - price2) / dfs;

// fin de modif
	if(fixing != 0.0)
	{
		if(ARMLOCAL_fixing (asOfDate,
							delivery,
							maturity,
							frequencyId,
							forwardCurve1Id,
							ccy1,
							fixing,
							facialMargin,
							C_result) == ARM_OK)
		{
			price += C_result.getDouble ();
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}
	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	MSG_printf_message (MSG_TRACE, "FRN Price Details : (DPStart = %lf) ; (Sensi1 = %lf) (valoMargin2 = %lf) ; (ccy1margin = %lf) ; (Price = %lf)", dfs, sensibility, valoMargin, ccy1margin, price);

	result.setRetCode (retCode);
	result.setDouble (price);

	return retCode;
}



long ARMLOCAL_FRNPrice (double asOfDate,
						double delivery,
						double maturity,
						const CCString& ccy1,
						long liborType1Id,
						long forwardCurve1Id,
						long discountCurve1Id,
						double facialMargin,
						double valoMargin,
						long frequencyId,
						const CCString& ccy2,
						long liborType2Id,
						long forwardCurve2Id,
						long discountCurve2Id,
						long amortizationId,
						long frequencyId2,
						double fixing,
						double spread,
						long solve,
						ARM_result& result)
{
	ARM_result C_result;
	
	double sensibility = 0.0;
	double sensibility2 = 0.0;
	double price = 0.0;
	double ccy1margin = 0.0;

	long retCode = ARM_KO;

    int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	if(ARMLOCAL_CptBPV (asOfDate,
						delivery,
						maturity,
						forwardCurve1Id,
						frequencyId,
						inputDayCount,
						ccy1,
						amortizationId,
						C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Sensibility : (%lf)", sensibility2);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(sensibility == 0.0)
	{
		result.setRetCode (ARM_KO);
		result.setMsg ("ARM_ERR: sensibility is nul");
		return ARM_KO;
	}

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		retCode = ARM_OK;
		ccy1margin = valoMargin;

		if (frequencyId != frequencyId2)
		{
			if(ARMLOCAL_CptBPV (asOfDate,
								delivery,
								maturity,
								forwardCurve1Id,
								frequencyId2,
								inputDayCount,
								ccy1,
								amortizationId,
								C_result) == ARM_OK)
			{
				sensibility2 = C_result.getDouble ();
				MSG_printf_message (MSG_TRACE, "FRN Sensibility : (%lf)", sensibility2);
			}
			else
			{
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if(sensibility2 == 0.0)
			{
				result.setRetCode (ARM_KO);
				result.setMsg ("ARM_ERR: sensibility2 is nul");
				return ARM_KO;
			}

			ccy1margin *= sensibility2 / sensibility;
		}
		price = facialMargin - ccy1margin;
	}
	else
	{
		retCode = ARMLOCAL_BasisSwap (asOfDate,
									  delivery,
									  maturity,
									  valoMargin,
									  ccy2,
									  liborType2Id,
									  forwardCurve2Id,
									  discountCurve2Id,
									  ccy1,
									  liborType1Id,
									  forwardCurve1Id,
									  discountCurve1Id,
									  amortizationId,
									  solve,
									  C_result);
		
		if(retCode == ARM_OK)
		{
			ccy1margin = C_result.getDouble ();
			price = facialMargin - ccy1margin;
			MSG_printf_message (MSG_TRACE, "FRN Basis Conversion Margin : valoMargin(2) = (%lf) ; implyMargin(1) = (%lf) ; method = (%ld)", valoMargin, C_result.getDouble (), solve);
		}
		else
		{
			result.setMsg (C_result.getMsg ());
		}
	}

	double dfs = 1.0;

	if(ARMLOCAL_DiscountPrice (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfs = result.getDouble ();
    

/*	if(ARMLOCAL_CptBPV (asOfDate, delivery, maturity, forwardCurve1Id, frequencyId2, inputDayCount, ccy1, C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
*/
	price = 100.0 + price * sensibility / dfs;

	if(fixing != 0.0)
	{
		if(ARMLOCAL_fixing (asOfDate,
							delivery,
							maturity,
							frequencyId,
							forwardCurve1Id,
							ccy1,
							fixing,
							facialMargin,
							C_result) == ARM_OK)
		{
			price += C_result.getDouble ();
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	MSG_printf_message (MSG_TRACE, "FRN Price Details : (DPStart = %lf) ; (Sensi1 = %lf) (valoMargin2 = %lf) ; (ccy1argin = %lf) ; (Price = %lf)", dfs, sensibility, valoMargin, ccy1margin, price);
		
	result.setRetCode (retCode);
	result.setDouble (price);

	return retCode;
}






/*****************************************************************

******************************************************************/
long ARMLOCAL_ASWPriceNew (double bondMaturity,
						   long cpnType,
						   double bondCoupon,
						   long bondFrequency,
						   long bondBase,
						   double bondMargin,
						   double bondRedemptionPrice,
						   double asOfDate,
						   double delivery,
						   long fixDecompFrequency,
						   long floatResetFreq,
						   long floatPayFreq,
						   const CCString& ccy1,
						   long liborType1Id,
						   long forwardCurve1Id,
						   long discountCurve1Id,
						   const CCString& ccy2,
						   long liborType2Id,
						   long forwardCurve2Id,
						   long discountCurve2Id,
						   long amortizationId,
						   long solve,
						   double minValue,
						   double maxValue,
						   ARM_result& result)
{
	ARM_result C_result;
	long ycModId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
		return ARM_KO;
	}

	if(!(ccy2 == "NONE") && (ccy1 != ccy2))
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
			return ARM_KO;
		}
	}

	double startDate = 27659.0;	// 22/09/1975!


	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		oneCurrency = true;
	}

	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double vFlag = 0.0;
	CCString id ("");
	double prevCpnDate = 0.0;

	double l_bondMargin = bondMargin;

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 bondMaturity,
										 bondMargin,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 amortizationId,
										 solve,
										 C_result);

		if(retCode == ARM_KO)
		{	
			MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}

		l_bondMargin = C_result.getDouble ();

		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (bondMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", bondMargin, l_bondMargin, solve);

	}

	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_PrevCpnDate (delivery,
								 bondMaturity,
								 bondFrequency,
								 fixIntRule,
								 ccy1,
								 C_result) == ARM_OK)
	{
		prevCpnDate = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): prevCpnDate = %lf", prevCpnDate);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  prevCpnDate,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  K_ANNUAL,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg3Id) == ARM_OK)
	{
		if(fixedLeg3Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg3Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg3Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if( (solve == BasisSwap_Method_Alg) || (solve == BasisSwap_Method_Alg_Spr))
	{
		retCode = ARMLOCAL_LiborAssetSwapPriceAlgNew (ycModId,
													  forwardCurve1Id,
													  fixedLeg1Id,
													  fixedLeg2Id,
													  fixedLeg3Id,
													  asOfDate,
													  delivery,
													  bondMaturity,
													  floatPayFreq,
													  l_bondMargin,
													  ccy1Id,
													  amortizationId,
													  bondRedemptionPrice,
													  result);
	}
	else
	{
		retCode = ARMLOCAL_LiborAssetSwapPrice (ycModId,
												forwardCurve1Id,
												fixedLeg1Id,
												fixedLeg2Id,
												fixedLeg3Id,
												asOfDate,
												delivery,
												bondMaturity,
												cpnType,
												bondCoupon,
												fixReceiveOrPay,
												floatPayFreq,
												K30_360,
												bondFrequency,
												fixIntRule,
												l_bondMargin,
												ccy1Id,
												amortizationId,
												bondRedemptionPrice,
												minValue,
												maxValue,
												result);
	}

	ARMLOCAL_FreeObject (ycModId, C_result);

	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return retCode;
}




long ARMLOCAL_BondASWPrice (long bondId,
							double bondMargin,
							double asOfDate,
							double delivery,
							long floatResetFreq,
							long floatPayFreq,
							const CCString& ccy1,
							long forwardCurve1Id,
							long discountCurve1Id,
							const CCString& ccy2,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							long solve,
							double minValue,
							double maxValue,
							ARM_result& result)
{
	ARM_result C_result;
	long ycModId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;

	ARM_Bond* bond = NULL;
	ARM_Currency* Ccy1Obj = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	bond = (ARM_Bond*) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0)
	{
		result.setMsg ("ARM_ERR: bond is not of a good type");
		return ARM_KO;
	}

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
		return ARM_KO;
	}

	Ccy1Obj = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Ccy1Obj, ARM_CURRENCY) == 0)
	{
		result.setMsg ("ARM_ERR: ccy1 is not of a good type");
		return ARM_KO;
	}

	char* ccy1name = ccy1.GetStr();

	long liborType1Id = GetDefaultIndexFromCurrency(ccy1name);

	if (ccy1name)
		delete ccy1name;
	ccy1name = NULL;

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		oneCurrency = true;
	}
	else
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_BondASWPrice (): %s", C_result.getMsg ());
			return ARM_KO;
		}
	}

	double startDate = 27659.0;	// 22/09/1975!

	char strDate[11];
	bond->GetMaturity().JulianToStrDate(strDate);
	double bondMaturity = Local_ARMDATE2XLDATE(strDate);
	double bondCoupon = bond->GetCoupon();
	double bondBase = bond->GetDayCount();
	double bondFrequency = bond->GetFrequency();
	double bondRedemptionPrice = bond->GetRedemptionValue();

	// Finesse du Swap (dixit AR)
	long fixDecompFrequency = 1;

	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double vFlag = 0.0;
	CCString id ("");
	double prevCpnDate = 0.0;

	double l_bondMargin = bondMargin;

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 bondMaturity,
										 bondMargin,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 amortizationId,
										 solve,
										 C_result);

		if(retCode == ARM_KO)
		{	
			MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}

		l_bondMargin = C_result.getDouble ();

		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (bondMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", bondMargin, l_bondMargin, solve);

	}

	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_PrevCpnDate (delivery,
								 bondMaturity,
								 bondFrequency,
								 fixIntRule,
								 ccy1,
								 C_result) == ARM_OK)
	{
		prevCpnDate = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): prevCpnDate = %lf", prevCpnDate);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  prevCpnDate,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  K_ANNUAL,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg3Id) == ARM_OK)
	{
		if(fixedLeg3Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg3Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg3Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if( (solve == BasisSwap_Method_Alg) || (solve == BasisSwap_Method_Alg_Spr))
	{
		retCode = ARMLOCAL_LiborAssetSwapPriceAlgNew (ycModId,
													  forwardCurve1Id,
													  fixedLeg1Id,
													  fixedLeg2Id,
													  fixedLeg3Id,
													  asOfDate,
													  delivery,
													  bondMaturity,
													  floatPayFreq,
													  l_bondMargin,
													  ccy1Id,
													  amortizationId,
													  bondRedemptionPrice,
													  result);
	}
	else
	{
		retCode = ARMLOCAL_LiborAssetSwapPrice (ycModId,
												forwardCurve1Id,
												fixedLeg1Id,
												fixedLeg2Id,
												fixedLeg3Id,
												asOfDate,
												delivery,
												bondMaturity,
												0,
												bondCoupon,
												fixReceiveOrPay,
												floatPayFreq,
												K30_360,
												bondFrequency,
												fixIntRule,
												l_bondMargin,
												ccy1Id,
												amortizationId,
												bondRedemptionPrice,
												minValue,
												maxValue,
												result);
	}

	ARMLOCAL_FreeObject (ycModId, C_result);

	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return retCode;
}


long ARMLOCAL_ASWPrice (double bondMaturity,
						double bondCoupon,
						long bondFrequency,
						long bondBase,
						double bondMargin,
						double bondRedemptionPrice,
						double asOfDate,
						double delivery,
						long fixDecompFrequency,
						long floatResetFreq,
						long floatPayFreq,
						const CCString& ccy1,
						long liborType1Id,
						long forwardCurve1Id,
						long discountCurve1Id,
						const CCString& ccy2,
						long liborType2Id,
						long forwardCurve2Id,
						long discountCurve2Id,
						long amortizationId,
						long solve,
						double minValue,
						double maxValue,
						ARM_result& result)
{
	ARM_result C_result;
	long ycModId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	
	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWMargin (): %s", C_result.getMsg ());
		return ARM_KO;
	}

	double startDate = 27659.0;	// 22/09/1975!
	
	if(ARMLOCAL_FIXEDLEG (startDate,
						  bondMaturity, 
						  K_RCV,
						  0,					
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  fixDecompFrequency,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								amortizationId,
								100.,
								C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_FIXEDLEG (startDate,
						  bondMaturity,
						  K_RCV,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  fixDecompFrequency,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		oneCurrency = true;
	}
	
	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

/*	double accrued1, accrued2;
	if(ARMLOCAL_ARM_Accrued (fixLeg1Id, delivery, ycModId, C_result) == ARM_OK)
	{
		accrued1 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}
	if(ARMLOCAL_ARM_Accrued (fixLeg2Id, delivery, ycModId, C_result) == ARM_OK)
	{
		accrued2 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}
	double supplFee = accrued1 - accrued2;

	if (supplFee != 0.0)
	{
		MSG_printf_message (MSG_TRACE, "ARM_ASWMargin - Correction de coupon : (Accr1 = %lf) ; (Accr2 = %lf)", accrued1, accrued2);
	}
	
	ARMLOCAL_FreeObject (fixLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixLeg2Id, C_result);
*/	
	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double vFlag = 0.0;
	CCString id ("");

	double l_bondMargin = bondMargin;

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwap (asOfDate,
									  delivery,
									  bondMaturity,
									  bondMargin,
									  ccy2,
									  liborType2Id,
									  forwardCurve2Id,
									  discountCurve2Id,
									  ccy1,
									  liborType1Id,
									  forwardCurve1Id,
									  discountCurve1Id,
									  amortizationId,
									  solve,
									  C_result);

		if(retCode == ARM_KO)
		{	
			MSG_printf_message (MSG_ERROR, "ARM_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}

		l_bondMargin = C_result.getDouble ();

		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (bondMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", bondMargin, l_bondMargin, solve);

	}

	if( (solve == BasisSwap_Method_Alg) || (solve == BasisSwap_Method_Alg_Spr))
	{
		retCode = ARMLOCAL_LiborAssetSwapPriceAlg (forwardCurve1Id,
												   asOfDate,
												   delivery,
												   bondMaturity,
												   bondCoupon,
												   bondBase,						  
												   K_RCV,
												   bondFrequency,
												   floatPayFreq,
												   K_UNADJUSTED,
												   l_bondMargin,
												   ccy1,
												   amortizationId,
												   bondRedemptionPrice,
												   result);
	}
	else
	{
		retCode = ARMLOCAL_LiborAssetSwapPrice (ycModId,
												forwardCurve1Id,
												fixedLeg1Id,
												fixedLeg2Id,
												fixedLeg3Id,
												asOfDate,
												delivery,
												bondMaturity,
												0,
												bondCoupon,
												fixReceiveOrPay,
												floatPayFreq,
												K30_360,
												bondFrequency,
												fixIntRule,
												l_bondMargin,
												ccy1Id,
												amortizationId,
												bondRedemptionPrice,
												minValue,
												maxValue,
												result);
	}

	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);

	ARMLOCAL_FreeObject (ycModId, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);

	return retCode;
}















/************************************************************************

*************************************************************************/
struct ARMLOCAL_LASM_arg
{
	long modelId;
	long zcId;
	long fixedLeg1Id;	
	long fixedLeg2Id;
	long fixedLeg3Id;
	double startDate;
	double delivery;
	double endDate;
	long cpnType;
	double fixedRate;
	long fixDayCountId;
	long fixReceiveOrPayId;
	long fixFrequencyId;
	long indexFrequencyId;
	long fixIntRuleId;
	double price;
	long discountCcyId;
	long amortizationId;
	double redemptionPrice;
	double margin;
}; 

typedef struct ARMLOCAL_LASM_arg ARMLOCAL_LASM_arg;


static void __stdcall LocalLiborAssetSwapMarginCaller (double price, double* diff, Nag_Comm* comm)
{
	ARM_result C_result;

	ARMLOCAL_LASM_arg* arg = (ARMLOCAL_LASM_arg*)comm->p;

	if (ARMLOCAL_LiborAssetSwapMarginAlgNew (arg->modelId,
											 arg->zcId,
											 arg->fixedLeg1Id,
											 arg->fixedLeg2Id,
											 arg->fixedLeg3Id,
											 arg->startDate,
											 arg->delivery,
											 arg->endDate,
											 arg->cpnType,
											 arg->fixedRate,
											 arg->fixDayCountId,
											 arg->fixReceiveOrPayId,
											 arg->fixFrequencyId,
											 arg->indexFrequencyId,
											 arg->fixIntRuleId,
											 price,
											 arg->discountCcyId,
											 arg->amortizationId,
											 arg->redemptionPrice,
											 C_result) == ARM_OK)
	{
		*diff = fabs (C_result.getDouble () - arg->margin) * 1000000.0;
		LASMC_ret = ARM_OK;
	}
	else
	{
		g_result->setMsg (C_result.getMsg ());
		g_result->setRetCode (C_result.getRetCode ());
		comm->p = (void*)g_result;
		*diff = 0.0;
		LASMC_ret = ARM_KO;
	}

	return;
}



static void __stdcall LocalLiborAssetSwapPriceCaller (double margin, double* diff, Nag_Comm* comm)
{
	ARM_result C_result;

	ARMLOCAL_LASM_arg* arg = (ARMLOCAL_LASM_arg*)comm->p;

	if (ARMLOCAL_LiborAssetSwapPriceAlgNew (arg->modelId,
											arg->zcId,
											arg->fixedLeg1Id,
											arg->fixedLeg2Id,
											arg->fixedLeg3Id,
											arg->startDate,
											arg->delivery,
											arg->endDate,
											arg->indexFrequencyId,
											margin,
											arg->discountCcyId,
											arg->amortizationId,
											arg->redemptionPrice,
											C_result) == ARM_OK)
	{
		*diff = fabs (C_result.getDouble () - arg->price) * 1000000.0;
		LASMC_ret = ARM_OK;
	}
	else
	{
		g_result->setMsg (C_result.getMsg ());
		g_result->setRetCode (C_result.getRetCode ());
		comm->p = (void*)g_result;
		*diff = 0.0;
		LASMC_ret = ARM_KO;
	}

	return;
}




/**************************************************************************************

***************************************************************************************/
long ARMLOCAL_LiborAssetSwapPrice (long ycModId,
								   long forwardCurve1Id,
								   long fixedLeg1Id,
								   long fixedLeg2Id,
								   long fixedLeg3Id,
								   double asOfDate,
								   double delivery,
								   double bondMaturity,
								   long cpnType,
								   double bondCoupon,
								   long fixReceiveOrPayId,
								   long floatPayFreq,
								   long fixDayCountId,
								   long fixFrequencyId,
								   long fixIntRuleId,
								   double bondMargin,
								   long ccy1Id,
								   long amortizationId,
								   double bondRedemptionPrice,
								   double minValue,
								   double maxValue,
								   ARM_result& result)
{
	// solver
	double a, b, e1, e2, x, f;
	Integer max_fun;
	Nag_Comm comm;
	static NagError fail;

	e1 = e2 = 0.0;
	a = minValue;
	b = maxValue;
	max_fun = 10000;
	fail.print = TRUE;

	ARMLOCAL_LASM_arg* l_arg = (ARMLOCAL_LASM_arg*)malloc (sizeof (ARMLOCAL_LASM_arg));
	l_arg->zcId = forwardCurve1Id;
	l_arg->modelId = ycModId;
	l_arg->fixedLeg1Id = fixedLeg1Id;
	l_arg->fixedLeg2Id = fixedLeg2Id;
	l_arg->fixedLeg3Id = fixedLeg3Id;
	l_arg->startDate = asOfDate;
	l_arg->delivery = delivery;
	l_arg->endDate = bondMaturity;
	l_arg->cpnType = cpnType;
	l_arg->fixedRate = bondCoupon;
	l_arg->fixReceiveOrPayId = fixReceiveOrPayId;
	l_arg->fixDayCountId = fixDayCountId;
	l_arg->fixFrequencyId = fixFrequencyId;
	l_arg->indexFrequencyId = floatPayFreq;
	l_arg->fixIntRuleId = fixIntRuleId;
	l_arg->discountCcyId = ccy1Id;
	l_arg->redemptionPrice = bondRedemptionPrice;
	l_arg->amortizationId = amortizationId;
	l_arg->margin = bondMargin;

	comm.p = (void*)l_arg;
   
	e04abc (LocalLiborAssetSwapMarginCaller, e1, e2, &a, &b, max_fun, &x, &f, &comm, &fail);

	if(LASMC_ret == ARM_OK)
	{
		result.setDouble (x);
		result.setRetCode (ARM_OK);
	}
	else
	{
		ARM_result* l_result = (ARM_result*)(comm.p);
		result.setRetCode (ARM_KO);
		result.setMsg (l_result->getMsg ());
	}

	if(l_arg)
	{
		free (l_arg);
	}

	return result.getRetCode ();
}




long ARMLOCAL_LiborAssetSwapMarginNumNew (long ycModId,
										  long forwardCurve1Id,
										  long fixedLeg1Id,
										  long fixedLeg2Id,
										  long fixedLeg3Id,
										  double asOfDate,
										  double delivery,
										  long bondMaturity,
										  long floatPayFreq,
										  double bondPrice,
										  long ccy1Id,
										  long amortizationId,
										  double bondRedemptionPrice,
										  double minValue,
										  double maxValue,
										  ARM_result& result)
{
	// solver
	double a, b, e1, e2, x, f;
	Integer max_fun;
	Nag_Comm comm;
	static NagError fail;

	e1 = e2 = 0.0;
	a = minValue;
	b = maxValue;
	max_fun = 100;
	fail.print = TRUE;

	ARMLOCAL_LASM_arg* l_arg = (ARMLOCAL_LASM_arg*)malloc (sizeof (ARMLOCAL_LASM_arg));
	l_arg->zcId = forwardCurve1Id;
	l_arg->modelId = ycModId;
	l_arg->fixedLeg1Id = fixedLeg1Id;
	l_arg->fixedLeg2Id = fixedLeg2Id;
	l_arg->fixedLeg3Id = fixedLeg3Id;
	l_arg->startDate = asOfDate;
	l_arg->delivery = delivery;
	l_arg->endDate = bondMaturity;
	l_arg->indexFrequencyId = floatPayFreq;
	l_arg->discountCcyId = ccy1Id;
	l_arg->redemptionPrice = bondRedemptionPrice;
	l_arg->amortizationId = amortizationId;
	l_arg->price = bondPrice;

	comm.p = (void*)l_arg;
   
	e04abc (LocalLiborAssetSwapPriceCaller, e1, e2, &a, &b, max_fun, &x, &f, &comm, &fail);

	if ( (fail.code == NE_NOERROR) || (fail.code == NW_COND_MIN) )
	{
		result.setDouble (x);
		result.setRetCode (ARM_OK);
	}
	else
	{
		ARM_result* l_result = (ARM_result*)(comm.p);
		result.setRetCode (ARM_KO);
		result.setMsg (l_result->getMsg ());
	}

	if(l_arg)
	{
		free (l_arg);
	}

	return result.getRetCode ();
}



/*******************************************************************************

********************************************************************************/
long ARMLOCAL_BasisSwapNew (double asOfDate,
							double delivery,
							double maturity,
							double margin1,
							long ccy1Id,
							long liborType1Id,
							long forwardCurve1Id,
							long discountCurve1Id,
							long ccy2Id,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							long solve,
							ARM_result& result)
{
	long retCode = ARM_KO;

	ARM_result C_result;

	long zcspreaded1Id = ARM_NULL_OBJECT_ID;
	long zcspreaded2Id = ARM_NULL_OBJECT_ID;
//	long ccy1Id = ARM_NULL_OBJECT_ID;
//	long ccy2Id = ARM_NULL_OBJECT_ID;
	char* tmp1 = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id))->GetCcyName();
	CCString ccy1 (tmp1);
	char* tmp2 = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy2Id))->GetCcyName();
	CCString ccy2 (tmp2);


	// Pour la m�thode analytique 'BasisSwap_Method_Alg' on construit la courbe discount
	// � partir de la courbe forward et la courbe de cotation des spreads
	// Dans les autres cas (it�rative, � la Summit ou variantes analytique), on utilise directement la courbe discount
	// qui est pass�e en param�tre !
	// les m�thodes 2 et 3 ne sont pas appelables dans les fonctions ASW ou FRN

	if(solve == BasisSwap_Method_Alg)
	{
		if(ARMLOCAL_zcspreaded (discountCurve1Id,
								forwardCurve1Id,
								asOfDate,
								K_MONTHLY,
								K_SEMIANNUAL,
								//ccy1Id,
								false,
								ccy1,
								C_result,
								zcspreaded1Id) == ARM_OK)
		{
			zcspreaded1Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	
		if(ARMLOCAL_zcspreaded (discountCurve2Id,
								forwardCurve2Id,
								asOfDate,
								K_MONTHLY,
								K_SEMIANNUAL,
								//ccy2Id,
								false,
								ccy2,
								C_result,
								zcspreaded2Id) == ARM_OK)
		{
			zcspreaded2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	switch(solve)
	{
		case BasisSwap_Method_Num:// m�thode it�rative (comme Summit) avec discount pass�s en param
			retCode = ARMLOCAL_BasisSwapNumNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id, // Libor pour les Fwd du swapleg
												discountCurve1Id, // Discount BS pass� en param
												ccy2Id,
												liborType2Id,
												forwardCurve2Id,
												discountCurve2Id,
												amortizationId,
												result);
			break;
		case BasisSwap_Method_Alg: // m�thode analytique avec discounts reconstruits avec Fwd et discount
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												discountCurve1Id, // Cotation BS
												zcspreaded1Id, // Discount BS recalcul�
												ccy2Id,
												liborType2Id,
												discountCurve2Id,
												zcspreaded2Id,
												amortizationId,
												BasisSwap_Method_Alg,
												result);
			break;
		case BasisSwap_Method_Num_Spr:// m�thode analytique avec discount (libor spread�) pass�s en param
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id,
												discountCurve1Id,
												ccy2Id,
												liborType2Id,
												forwardCurve2Id, // Cotation BS
												discountCurve2Id, // Discount BS pass� en param
												amortizationId,
												BasisSwap_Method_Num_Spr,
												result);
			break;
		case BasisSwap_Method_Alg_Spr:// m�thode analytique avec discount (libor) pass�s en param 
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id, // Cotation BS
												discountCurve1Id, // Discount Libor pass� en param
												ccy2Id,
												liborType2Id,
												forwardCurve2Id,
												discountCurve2Id,
												amortizationId,
												BasisSwap_Method_Alg_Spr,
												result);
			break;
		default:
			retCode = ARM_KO;
			break;
	}

	if(solve == BasisSwap_Method_Alg)
	{
		ARMLOCAL_FreeObject (zcspreaded1Id, C_result);
		ARMLOCAL_FreeObject (zcspreaded2Id, C_result);
/*		ARMLOCAL_FreeObject (ccy1Id, C_result);
		ARMLOCAL_FreeObject (ccy2Id, C_result);*/
	}

	return retCode;
}


long ARMLOCAL_BasisSwap (double asOfDate,
						 double delivery,
						 double maturity,
						 double margin1,
						 const CCString& ccy1,
						 long liborType1Id,
						 long forwardCurve1Id,
						 long discountCurve1Id,
						 const CCString& ccy2,
						 long liborType2Id,
						 long forwardCurve2Id,
						 long discountCurve2Id,
						 long amortizationId,
						 long solve,
						 ARM_result& result)
{
	long retCode = ARM_KO;

	ARM_result C_result;

	long zcspreaded1Id = ARM_NULL_OBJECT_ID;
	long zcspreaded2Id = ARM_NULL_OBJECT_ID;
	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;


	// Pour la m�thode analytique 'BasisSwap_Method_Alg' on construit la courbe discount
	// � partir de la courbe forward et la courbe de cotation des spreads
	// Dans les autres cas (it�rative, � la Summit ou variantes analytique), on utilise directement la courbe discount
	// qui est pass�e en param�tre !
	// les m�thodes 2 et 3 ne sont pas appelables dans les fonctions ASW ou FRN

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy2 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy2,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy2Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy2,
								   C_result,
								   ccy2Id);
	}
	if(retCode == ARM_OK)
	{
		ccy2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(solve == BasisSwap_Method_Alg)
	{
		if(ARMLOCAL_zcspreaded (discountCurve1Id,
								forwardCurve1Id,
								asOfDate,
								K_MONTHLY,
								K_SEMIANNUAL,
								//ccy1Id,
								false,
								ccy1,
								C_result,
								zcspreaded1Id) == ARM_OK)
		{
			zcspreaded1Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	
		if(ARMLOCAL_zcspreaded (discountCurve2Id,
								forwardCurve2Id,
								asOfDate,
								K_MONTHLY,
								K_SEMIANNUAL,
								//ccy2Id,
								false,
								ccy2,
								C_result,
								zcspreaded2Id) == ARM_OK)
		{
			zcspreaded2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARM_BasisSwap (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	switch(solve)
	{
		case BasisSwap_Method_Num:// m�thode it�rative (comme Summit) avec discount pass�s en param
			retCode = ARMLOCAL_BasisSwapNumNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id, // Libor pour les Fwd du swapleg
												discountCurve1Id, // Discount BS pass� en param
												ccy2Id,
												liborType2Id,
												forwardCurve2Id,
												discountCurve2Id,
												amortizationId,
												result);
			break;
		case BasisSwap_Method_Alg: // m�thode analytique avec discounts reconstruits avec Fwd et discount
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												discountCurve1Id, // Cotation BS
												zcspreaded1Id, // Discount BS recalcul�
												ccy2Id,
												liborType2Id,
												discountCurve2Id,
												zcspreaded2Id,
												amortizationId,
												BasisSwap_Method_Alg,
												result);
			break;
		case BasisSwap_Method_Num_Spr:// m�thode analytique avec discount (libor spread�) pass�s en param
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id,
												discountCurve1Id,
												ccy2Id,
												liborType2Id,
												forwardCurve2Id, // Cotation BS
												discountCurve2Id, // Discount BS pass� en param
												amortizationId,
												BasisSwap_Method_Num_Spr,
												result);
			break;
		case BasisSwap_Method_Alg_Spr:// m�thode analytique avec discount (libor) pass�s en param 
			retCode = ARMLOCAL_BasisSwapAlgNew (asOfDate,
												delivery,
												maturity,
												margin1,
												ccy1Id,
												liborType1Id,
												forwardCurve1Id, // Cotation BS
												discountCurve1Id, // Discount Libor pass� en param
												ccy2Id,
												liborType2Id,
												forwardCurve2Id,
												discountCurve2Id,
												amortizationId,
												BasisSwap_Method_Alg_Spr,
												result);
			break;
		default:
			retCode = ARM_KO;
			break;
	}

/*	if(solve == BasisSwap_Method_Alg)
	{
*/		ARMLOCAL_FreeObject (zcspreaded1Id, C_result);
		ARMLOCAL_FreeObject (zcspreaded2Id, C_result);
//	}

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return retCode;
}






long ARMLOCAL_LiborAssetSwapPriceAlg (long zcId,
									  double startDate,
									  double delivery,
									  double endDate,
									  double fixedRate,
									  long fixDayCountId,
									  long fixReceiveOrPayId, 
									  long fixFrequencyId,
									  long indexFrequencyId,
									  long fixIntRuleId,
									  double margin,
									  const CCString& discountCcy,
									  long amortizationId,
									  double redemptionPrice,
									  ARM_result& result)
{
	ARM_result C_result;
	long retCode = ARM_OK;

	double prevCpnDate = 0.0;
	double endDateAdjusted = 0.0;

	double sensibility = 0.0;
	double price = 0.0;
	
	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;
	long discountCcyId = ARM_NULL_OBJECT_ID;
	long modelId = ARM_NULL_OBJECT_ID;

	if(discountCcy == "GBP")
	{
		retCode = ARMLOCAL_CCY (discountCcy,
								zcId,
								1.0,
								KACTUAL_365,
								C_result,
								discountCcyId);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (discountCcy,
								   C_result,
								   discountCcyId);
	}
	if(retCode == ARM_OK)
	{
		if(discountCcyId == ARM_NULL_OBJECT_ID)
		{
			discountCcyId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
	
	if(ARMLOCAL_FIXEDLEG (27659.0,
						  endDate,
						  fixReceiveOrPayId,
						  0,
						  fixedRate,
						  K30_360,
						  fixFrequencyId,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //discountCcyId,
						  false,
						  discountCcy,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_PrevCpnDate (delivery,
								 endDate,
								 fixFrequencyId,
								 fixIntRuleId,
								 discountCcy,
								 C_result) == ARM_OK)
	{
		prevCpnDate = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): prevCpnDate = %lf", prevCpnDate);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  prevCpnDate,
						  fixReceiveOrPayId,
						  0,
						  fixedRate,
						  K30_360,
						  fixFrequencyId,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //discountCcyId,
						  false,
						  discountCcy,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ycmod (zcId,
					   ARM_NULL_OBJECT,
					   C_result,
					   modelId) == ARM_OK)
	{
		if(modelId == ARM_NULL_OBJECT_ID)
		{
			modelId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLeg1Id,
						   modelId,
						   C_result) == ARM_OK)
	{
		price += C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 1 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLeg2Id,
						   modelId,
						   C_result) == ARM_OK)
	{
		price -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 2 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (endDate,
								 discountCcy,
								 1,
								 C_result) == ARM_OK)
	{
		endDateAdjusted = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId == -1)
	{
		if(ARMLOCAL_DiscountPrice (zcId,
								   (endDateAdjusted - startDate) / 365.0,
								   C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			price += C_result.getDouble () * redemptionPrice;
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble () * redemptionPrice);
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}
	else
	{
		if(ARMLOCAL_DiscountPriceRefvalue(zcId,
										  amortizationId,
										  discountCcyId,
										  startDate,
										  endDateAdjusted,
										  C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			price += C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble ());
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

    int inputDayCount = KACTUAL_360;

    if(discountCcy == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	if(ARMLOCAL_CptBPV (startDate,
						delivery,
						endDate,
						zcId,
						indexFrequencyId,
						inputDayCount,
						discountCcy,
						amortizationId,
						C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): sensibility = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	price -= margin * sensibility;

	if((ARMLOCAL_DiscountPrice (zcId,
								(delivery - startDate) / 365.0,
								C_result) == ARM_OK) && (C_result.getDouble () != 0.0))
	{
		price /= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): correction J-3 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  endDate,
						  fixReceiveOrPayId,
						  0,
						  fixedRate,
						  fixDayCountId,
						  fixFrequencyId,
						  K_ANNUAL,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //discountCcyId,
						  false,
						  discountCcy,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg3Id) == ARM_OK)
	{
		if(fixedLeg3Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg3Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg3Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_Accrued (fixedLeg3Id,
							 delivery,
							 modelId,
							 C_result) == ARM_OK)
	{
		price -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): accrued = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);
	ARMLOCAL_FreeObject (discountCcyId, C_result);
	ARMLOCAL_FreeObject (modelId, C_result);
	
	result.setRetCode (retCode);
	result.setDouble (price);

	return retCode;
}




long ARMLOCAL_LiborAssetSwapPriceAlgNew (long modelId,
										 long zcId,
										 long fixedLeg1Id,
										 long fixedLeg2Id,
										 long fixedLeg3Id,
										 double startDate,
										 double delivery,
										 double endDate,
										 long indexFrequencyId,
										 double margin,
										 long discountCcyId,
										 long amortizationId,
										 double redemptionPrice,
										 ARM_result& result)
{
	ARM_result C_result;
	long retCode = ARM_OK;

	double prevCpnDate = 0.0;
	double endDateAdjusted = 0.0;

	double sensibility = 0.0;
	double price = 0.0;

//	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
//	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
//	long fixedLeg3Id = ARM_NULL_OBJECT_ID;
//	long discountCcyId = ARM_NULL_OBJECT_ID;
//	long modelId = ARM_NULL_OBJECT_ID;

	// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
	// Modif du 28/11/01
//	char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId))->GetCcyName();
	ARM_Currency* ccy = (ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId);

	char* tmp = ccy->GetCcyName();
	CCString discountCcy(tmp);

    ARM_INDEX_TYPE K_INDEX_TYPE = GetDefaultIndexFromCurrency(ccy->GetCcyName());
    char* payCalName   = ccy->GetPayCalName(K_INDEX_TYPE);
	CCString cPayCalName(payCalName);

	if (payCalName)
		free(payCalName);
	payCalName = NULL;

/*	if(discountCcy == "GBP")
	{
		retCode = ARMLOCAL_CCY (discountCcy, zcId, 1, 2, C_result, discountCcyId);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (discountCcy, C_result, discountCcyId);
	}
	if(retCode == ARM_OK)
	{
		if(discountCcyId == ARM_NULL_OBJECT_ID)
		{
			discountCcyId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
*/


	if(ARMLOCAL_ARM_Price (fixedLeg1Id,
						   modelId,
						   C_result) == ARM_OK)
	{
		price += C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 1 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLeg2Id,
						   modelId,
						   C_result) == ARM_OK)
	{
		price -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 2 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (endDate,
								 cPayCalName,
								 K_MOD_FOLLOWING,
								 C_result) == ARM_OK)
	{
		endDateAdjusted = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId == -1)
	{
		if(ARMLOCAL_DiscountPrice (zcId,
								   (endDateAdjusted - startDate) / 365.0,
								   C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			price += C_result.getDouble () * redemptionPrice;
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble () * redemptionPrice);
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}
	else
	{
		if(ARMLOCAL_DiscountPriceRefvalue(zcId,
										  amortizationId,
										  discountCcyId,
										  startDate,
										  endDateAdjusted,
										  C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			price += C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble ());
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}
/*
    int inputDayCount = KACTUAL_360;

    if(discountCcy == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }
*/
	if(ARMLOCAL_CptBPVNew (startDate,
						   delivery,
						   endDate,
						   zcId,
						   indexFrequencyId,
						   ccy->GetLiborIndexDayCount(),
						   discountCcyId,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): sensibility = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	price -= margin * sensibility;

	if((ARMLOCAL_DiscountPrice (zcId,
								(delivery - startDate) / 365.0,
								C_result) == ARM_OK)
								&&
								(C_result.getDouble () != 0.0))
	{
		price /= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): correction J-3 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Accrued (fixedLeg3Id,
							 delivery,
							 modelId,
							 C_result) == ARM_OK)
	{
		price -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): accrued = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}


//	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
//	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
//	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);
//	ARMLOCAL_FreeObject (discountCcyId, C_result);
//	ARMLOCAL_FreeObject (modelId, C_result);
	
	result.setRetCode (retCode);
	result.setDouble (price);

	return retCode;
}



long ARMLOCAL_ARM_NextCpnDate (double asOfDate,
							   double maturity,
							   long frequencyId,
							   long ruleId,
							   const CCString& ccy,
							   ARM_result& result,
							   long intruleId)
{
	ARM_result C_result;

	if (frequencyId == K_ZEROCOUPON)
	{
		result.setRetCode (ARM_OK);
		result.setDouble (asOfDate);

		return ARM_OK;
	}

	int d, m, y;

	DAT_ssdate_to_struct (asOfDate, &y, &m, &d);
	ARM_Local_Date C_asOfDate (d, m, y);

	DAT_ssdate_to_struct (maturity, &y, &m, &d);
	ARM_Local_Date C_nextDate (d, m, y);
	ARM_Local_Date C_tempDate (d, m, y);

	double l_nextDate;

	int i = 1;

	// while((double)C_nextDate >= (double)C_asOfDate) Modif du 22/06/01 pour corriger bug sur date anniversaire
	if ((double)C_asOfDate<0)
		return ARM_KO;

	while((double)C_nextDate > (double)C_asOfDate)
	{
		d = C_nextDate.GetDay ();
		m = C_nextDate.GetMonth ();
		y = C_nextDate.GetYear ();

		l_nextDate = DAT_struct_to_ssdate (y, m, d);
	
		C_nextDate = C_tempDate;
		C_nextDate.AddPeriodMult (-1 * frequencyId,i);
		i++;
	}

	if (intruleId == K_ADJUSTED)
	{
		if(ARMLOCAL_ADJUSTTOBUSDATE (l_nextDate,
									 ccy,
									 ruleId,
									 C_result) == ARM_OK)
		{
			l_nextDate = Local_ARMDATE2XLDATE(C_result.getString ());
		}
	}

	result.setRetCode (ARM_OK);
	result.setDouble (l_nextDate);

	return ARM_OK;
}



long ARMLOCAL_ARM_PrevCpnDate (double asOfDate,
							   double maturity,
							   long frequencyId,
							   long ruleId,
							   const CCString& ccy,
							   ARM_result& result,
							   long intruleId)
{
	ARM_result C_result;

	if (frequencyId == K_ZEROCOUPON)
	{
		result.setRetCode (ARM_OK);
		result.setDouble (asOfDate);

		return ARM_OK;
	}

	int d, m, y;

	DAT_ssdate_to_struct (asOfDate, &y, &m, &d);
	ARM_Local_Date C_asOfDate (d, m, y);

	DAT_ssdate_to_struct (maturity, &y, &m, &d);
	ARM_Local_Date C_nextDate (d, m, y);
	ARM_Local_Date C_tempDate (d, m, y);
	
	int i(1);

	if ((double)C_asOfDate<0)
		return ARM_KO;

	// while((double)C_nextDate >= (double)C_asOfDate) Modif du 12/06/01 pour corriger bug sur date anniversaire
	while((double)C_nextDate > (double)C_asOfDate)
	{
		C_nextDate = C_tempDate;
		C_nextDate.AddPeriodMult (-1 * frequencyId,i);
	
		i++;
	}

	d = C_nextDate.GetDay ();
	m = C_nextDate.GetMonth ();
	y = C_nextDate.GetYear ();

	double l_nextDate = DAT_struct_to_ssdate (y, m, d);
	
	if (intruleId == K_ADJUSTED)
	{
		if(ARMLOCAL_ADJUSTTOBUSDATE (l_nextDate,
									 ccy,
									 ruleId,
									 C_result) == ARM_OK)
		{
			l_nextDate = Local_ARMDATE2XLDATE(C_result.getString ());
		}
	}

	result.setRetCode (ARM_OK);
	result.setDouble (l_nextDate);

	return ARM_OK;
}



long ARMLOCAL_CptBPV (double asOfDate,
					  double delivery,
					  double maturity,
					  long zcId,
					  long frequencyId,
					  long dayCountId,
					  const CCString& ccy,
					  long amortizationId,
					  ARM_result& result)
{
	ARM_result C_result;
	long retCode = ARM_OK;
	double sensibility = 0.0;
	
	long fixedLegId = ARM_NULL_OBJECT_ID;
	long ccyId = ARM_NULL_OBJECT_ID;
	long ycModId = ARM_NULL_OBJECT_ID;
	
	if(ccy == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy,
								zcId,
								1.0,
								KACTUAL_365,
								C_result,
								ccyId);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy,
								   C_result,
								   ccyId);
	}
	if(retCode == ARM_OK)
	{
		if(ccyId == ARM_NULL_OBJECT_ID)
		{
			ccyId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (delivery,
						  maturity,
						  K_RCV,
						  0,
						  0.01,
						  dayCountId,
						  frequencyId,
						  K_ANNUAL,
						  K_ARREARS,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  //ccyId,
						  false,
						  ccy,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLegId) == ARM_OK)
	{
		if(fixedLegId == ARM_NULL_OBJECT_ID)
		{
			fixedLegId = C_result.getLong ();
		}
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLegId,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ycmod (zcId,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLegId,
						   ycModId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
	
	ARMLOCAL_FreeObject (fixedLegId, C_result);
	ARMLOCAL_FreeObject (ccyId, C_result);
	ARMLOCAL_FreeObject (ycModId, C_result);

	result.setRetCode (retCode);
	result.setDouble (sensibility);

	return retCode;
}


long ARMLOCAL_CptBPVNew (double asOfDate,
						 double delivery,
						 double maturity,
						 long zcId,
						 long frequencyId,
						 long dayCountId,
						 long ccyId,
						 long amortizationId,
						 ARM_result& result)
{
	if (delivery == maturity)
	{
		result.setDouble(0.0);
		return ARM_OK;
	}

	ARM_result C_result;
	long retCode = ARM_OK;
	double sensibility = 0.0;
	
	long fixedLegId = ARM_NULL_OBJECT_ID;
//	long ccyId = ARM_NULL_OBJECT_ID;
	long ycModId = ARM_NULL_OBJECT_ID;
	
/*	if(ccy == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy, zcId, 1, 2, C_result, ccyId);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy, C_result, ccyId);
	}
	if(retCode == ARM_OK)
	{
		if(ccyId == ARM_NULL_OBJECT_ID)
		{
			ccyId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
*/
	char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId))->GetCcyName();
	CCString ccy (tmp);

	if(ARMLOCAL_FIXEDLEG (delivery,
						  maturity,
						  K_RCV,
						  0,
						  0.01,
						  dayCountId,
						  frequencyId,
						  K_ANNUAL,
						  K_ARREARS,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  //ccyId,
						  false,
						  ccy,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLegId) == ARM_OK)
	{
		if(fixedLegId == ARM_NULL_OBJECT_ID)
		{
			fixedLegId = C_result.getLong ();
		}
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLegId,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ycmod (zcId,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLegId,
						   ycModId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}
	
	ARMLOCAL_FreeObject (fixedLegId, C_result);
//	ARMLOCAL_FreeObject (ccyId, C_result);
	ARMLOCAL_FreeObject (ycModId, C_result);

	result.setRetCode (retCode);
	result.setDouble (sensibility);

	return retCode;
}


long ARMLOCAL_BasisSwapNum (double asOfDate,
							double delivery,
							double maturity,
							double margin1,
							const CCString& ccy1,
							long liborType1Id,
							long forwardCurve1Id,
							long discountCurve1Id,
							const CCString& ccy2,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							ARM_result& result)
{
	long retCode = ARM_KO;

	static int first_call_flag = 0;
	static long dccyDefaultId = ARM_NULL_OBJECT_ID;
	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;
	double adjustedEnd1;
	double adjustedEnd2;
	long inputDayCount;
	
	ARM_result C_result;

	if(first_call_flag == 0)
	{
		if(ARMLOCAL_ISOCCY (CURRENCY_USD,
							C_result,
							dccyDefaultId) == ARM_OK)
		{
			dccyDefaultId = C_result.getLong ();
		}
		else
		{
			g_price1 = -1.0;
			g_price2 = -1.0;
			g_nbiter = -1.0;
			g_minimum = -1.0;

			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		first_call_flag = 1;
		sensi_std = 0.0;
		sensi_GBP = 0.0;
	}

	/*g_dccy1Id = dccyDefaultId;
	g_dccy2Id = dccyDefaultId;*/
	
	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy2 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy2,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy2Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy2,
								   C_result,
								   ccy2Id);
	}
	if(retCode == ARM_OK)
	{
		ccy2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	
	if(ARMLOCAL_DiscountPrice (discountCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		startXNL1 = -100.0 * C_result.getDouble ();
	}
	else
	{
		g_price1 = -1.0;
		g_price2 = -1.0;
		g_nbiter = -1.0;
		g_minimum = -1.0;

		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy1,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd1 = Local_ARMDATE2XLDATE (C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}


	if(ARMLOCAL_DiscountPrice (discountCurve1Id,
							   (adjustedEnd1 - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		endXNL1 = 100.0 * C_result.getDouble ();
	}
	else
	{
		g_price1 = -1.0;
		g_price2 = -1.0;
		g_nbiter = -1.0;
		g_minimum = -1.0;

		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy1 == "GBP")
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   maturity,
								   discountCurve1Id,
								   2, //floatPayFreq,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);

		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			// MSG_printf_message (MSG_TRACE, "ARM_LiborAssetSwapPriceAlg (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		
		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   maturity,
								   discountCurve1Id,
								   2, //floatPayFreq,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			// MSG_printf_message (MSG_TRACE, "ARM_LiborAssetSwapPriceAlg (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_computeprice1 (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	if(ARMLOCAL_DiscountPrice (discountCurve2Id,
							   (delivery - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		startXNL2 = -100.0 * C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}


	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy2,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd2 = Local_ARMDATE2XLDATE (C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountPrice (discountCurve2Id,
							   (adjustedEnd2 - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		endXNL2 = 100.0 * C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy2 == "GBP")
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   maturity,
								   discountCurve2Id,
								   2, //floatPayFreq,
								   inputDayCount,
								   ccy2,
								   amortizationId,
								   C_result);

		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		
		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   maturity,
								   discountCurve2Id,
								   2, //floatPayFreq,
								   inputDayCount,
								   ccy2,
								   amortizationId,
								   C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_computeprice1 (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	if(g_index1Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_index1Id, C_result);
		g_index1Id = ARM_NULL_OBJECT_ID;
	}

	bool dccyDefaultIsObject = false;
	CCString dccyDefault = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(dccyDefaultId))->GetCcyName();

	if(ARMLOCAL_LIBOR (liborType1Id, dccyDefaultIsObject, dccyDefault, -1, -1, -1 ,K_ADJUSTED, C_result, g_index1Id) == ARM_OK)
	{
		g_index1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_index2Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_index2Id, C_result);
		g_index2Id = ARM_NULL_OBJECT_ID;
	}

	if(ARMLOCAL_LIBOR (liborType2Id, dccyDefaultIsObject, dccyDefault, -1, -1,-1, K_ADJUSTED,C_result, g_index2Id) == ARM_OK)
	{
		g_index2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_model1Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_model1Id, C_result);
		g_model1Id = ARM_NULL_OBJECT_ID;
	}
	if(ARMLOCAL_GTWOYC (c_a, c_s, discountCurve1Id, forwardCurve1Id, c_cor, c_a, c_s, C_result, g_model1Id) == ARM_OK)
	{
		g_model1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_model2Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_model2Id, C_result);
		g_model2Id = ARM_NULL_OBJECT_ID;
	}
	if(ARMLOCAL_GTWOYC (c_a, c_s, discountCurve2Id, forwardCurve2Id, c_cor, c_a, c_s, C_result, g_model2Id) == ARM_OK)
	{
		g_model2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	// solver
	double a, b, e1, e2, x, f;
	Integer max_fun;
	Nag_Comm comm;
	static NagError fail;

	e1 = e2 = 0.0;
	a = -1000.0;
	b = 2000.0;
	max_fun = 100;
	fail.print = TRUE;

	ARMLOCAL_CDP_arg* l_arg = (ARMLOCAL_CDP_arg*)malloc (sizeof (ARMLOCAL_CDP_arg));
	l_arg->margin1 = margin1;
	l_arg->delivery = delivery;
	l_arg->maturity = maturity;
	l_arg->ccy1 = new CCString (ccy1);
	l_arg->ccy1Id = ccy1Id;
	l_arg->ccy2 = new CCString (ccy2);
	l_arg->ccy2Id = ccy2Id;
	
	comm.p = (void*)l_arg;
      
	e04abc (Local_computeDiffPrice, e1, e2, &a, &b, max_fun, &x, &f, &comm, &fail);
     
	if(CDP_ret == ARM_OK)
	{
		Local_computePrice1 (l_arg, margin1, &g_price1);
		Local_computePrice2 (l_arg, x, &g_price2);
		g_nbiter = comm.nf;
		g_minimum = f / 1000000.0;

		result.setDouble (x);
		result.setRetCode (ARM_OK);

		MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum - LEG1 Infos : (ccy = %s) ; (margin = %lf) ; (StartXNL = %lf) ; (EndXNL = %lf) ; (LegValue = %lf)", (const char*)ccy1, margin1, startXNL1, endXNL1, g_price1);
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum - LEG2 Infos : (ccy = %s) ; (margin = %lf) ; (StartXNL = %lf) ; (EndXNL = %lf) ; (LegValue = %lf)", (const char*)ccy2, x, startXNL2, endXNL2, g_price2);

	}
	else
	{
		ARM_result* l_result = (ARM_result*)(comm.p);
		result.setRetCode (ARM_KO);
		result.setMsg (l_result->getMsg ());
	}

	if(l_arg)
	{
		free (l_arg->ccy1);
		free (l_arg->ccy2);
		free (l_arg);
	}

	Local_destroySwapLegs ();

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return result.getRetCode ();
}


long ARMLOCAL_BasisSwapNumNew (double asOfDate,
							   double delivery,
							   double maturity,
							   double margin1,
							   long ccy1Id,
							   long liborType1Id,
							   long forwardCurve1Id,
							   long discountCurve1Id,
							   long ccy2Id,
							   long liborType2Id,
							   long forwardCurve2Id,
							   long discountCurve2Id,
							   long amortizationId,
							   ARM_result& result)
{
	long retCode = ARM_KO;

	static int first_call_flag = 0;
	static long dccyDefaultId = ARM_NULL_OBJECT_ID;
//	long ccy1Id = ARM_NULL_OBJECT_ID;
//	long ccy2Id = ARM_NULL_OBJECT_ID;
	double adjustedEnd1;
	double adjustedEnd2;
	long inputDayCount;
	
	ARM_result C_result;

	if(first_call_flag == 0)
	{
		if(ARMLOCAL_ISOCCY (CURRENCY_USD,
							C_result,
							dccyDefaultId) == ARM_OK)
		{
			dccyDefaultId = C_result.getLong ();
		}
		else
		{
			g_price1 = -1.0;
			g_price2 = -1.0;
			g_nbiter = -1.0;
			g_minimum = -1.0;

			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		first_call_flag = 1;
		sensi_std = 0.0;
		sensi_GBP = 0.0;
	}

	/*g_dccy1Id = dccyDefaultId;
	g_dccy2Id = dccyDefaultId;*/
	
/*	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1, forwardCurve1Id, 1, 2, C_result, ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1, C_result, ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy2 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy2, forwardCurve1Id, 1, 2, C_result, ccy2Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy2, C_result, ccy2Id);
	}
	if(retCode == ARM_OK)
	{
		ccy2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

*/	
	if(ARMLOCAL_DiscountPrice (discountCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		startXNL1 = -100.0 * C_result.getDouble ();
	}
	else
	{
		g_price1 = -1.0;
		g_price2 = -1.0;
		g_nbiter = -1.0;
		g_minimum = -1.0;

		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
	char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id))->GetCcyName();
	CCString ccy1 (tmp);

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy1,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd1 = Local_ARMDATE2XLDATE (C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}


	if(ARMLOCAL_DiscountPrice (discountCurve1Id,
							   (adjustedEnd1 - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		endXNL1 = 100.0 * C_result.getDouble ();
	}
	else
	{
		g_price1 = -1.0;
		g_price2 = -1.0;
		g_nbiter = -1.0;
		g_minimum = -1.0;

		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy1 == "GBP")
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
									  delivery,
									  maturity,
									  discountCurve1Id,
									  2, //floatPayFreq,
									  inputDayCount,
									  ccy1Id,
									  amortizationId,
									  C_result);

		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			// MSG_printf_message (MSG_TRACE, "ARM_LiborAssetSwapPriceAlg (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		
		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
									  delivery,
									  maturity,
									  discountCurve1Id,
									  2, //floatPayFreq,
									  inputDayCount,
									  ccy1Id,
									  amortizationId,
									  C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			// MSG_printf_message (MSG_TRACE, "ARM_LiborAssetSwapPriceAlg (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_computeprice1 (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	if(ARMLOCAL_DiscountPrice (discountCurve2Id,
							   (delivery - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		startXNL2 = -100.0 * C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
	tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy2Id))->GetCcyName();
	CCString ccy2(tmp);

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy2,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd2 = Local_ARMDATE2XLDATE (C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountPrice (discountCurve2Id,
							   (adjustedEnd2 - asOfDate) / 365.0,
							   C_result) == ARM_OK)
	{
		endXNL2 = 100.0 * C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ccy2 == "GBP")
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
									  delivery,
									  maturity,
									  discountCurve2Id,
									  2, //floatPayFreq,
									  inputDayCount,
									  ccy2Id,
									  amortizationId,
									  C_result);

		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
		
		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
									  delivery,
									  maturity,
									  discountCurve2Id,
									  2, //floatPayFreq,
									  inputDayCount,
									  ccy2Id,
									  amortizationId,
									  C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_computeprice1 (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	if(g_index1Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_index1Id, C_result);
		g_index1Id = ARM_NULL_OBJECT_ID;
	}
	
	bool ccyIsObject = false;

	if(ARMLOCAL_LIBOR (liborType1Id,
					   ccyIsObject,
					   ccy1,
					   -1,
					   -1,
                       -1,
					   K_ADJUSTED,
					   C_result,
					   g_index1Id) == ARM_OK)
	{
		g_index1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_index2Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_index2Id, C_result);
		g_index2Id = ARM_NULL_OBJECT_ID;
	}

	if(ARMLOCAL_LIBOR (liborType2Id,
					   ccyIsObject,
					   ccy2,
					   -1,
					   -1,
                       -1,
					   K_ADJUSTED,
					   C_result,
					   g_index2Id) == ARM_OK)
	{
		g_index2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_model1Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_model1Id, C_result);
		g_model1Id = ARM_NULL_OBJECT_ID;
	}
	/*
	if(ARMLOCAL_GTWOYC (c_a,
						c_s,
						discountCurve1Id,
						forwardCurve1Id,
						c_cor,
						c_a,
						c_s,
						C_result,
						g_model1Id) == ARM_OK)
	*/
	if(ARMLOCAL_ycmod(forwardCurve1Id, discountCurve1Id, C_result) == ARM_OK)
	{
		g_model1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARM_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(g_model2Id != ARM_NULL_OBJECT_ID)
	{
		ARMLOCAL_FreeObject (g_model2Id, C_result);
		g_model2Id = ARM_NULL_OBJECT_ID;
	}
	/*
	if(ARMLOCAL_GTWOYC (c_a,
						c_s,
						discountCurve2Id,
						forwardCurve2Id,
						c_cor,
						c_a,
						c_s,
						C_result,
						g_model2Id) == ARM_OK)
	*/
	if(ARMLOCAL_ycmod(forwardCurve2Id, discountCurve2Id, C_result) == ARM_OK)
	{
		g_model2Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapNum (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	// solver
	double a, b, e1, e2, x, f;
	Integer max_fun;
	Nag_Comm comm;
	static NagError fail;

	e1 = e2 = 0.0;
	a = -1000.0;
	b = 2000.0;
	max_fun = 100;
	fail.print = TRUE;

	ARMLOCAL_CDP_arg* l_arg = (ARMLOCAL_CDP_arg*)malloc (sizeof (ARMLOCAL_CDP_arg));
	l_arg->margin1 = margin1;
	l_arg->delivery = delivery;
	l_arg->maturity = maturity;
	l_arg->ccy1 = new CCString (ccy1);
	l_arg->ccy1Id = ccy1Id;
	l_arg->ccy2 = new CCString (ccy2);
	l_arg->ccy2Id = ccy2Id;
	
	comm.p = (void*)l_arg;
      
	e04abc (Local_computeDiffPrice, e1, e2, &a, &b, max_fun, &x, &f, &comm, &fail);
     
	if(CDP_ret == ARM_OK)
	{
		Local_computePrice1 (l_arg, margin1, &g_price1);
		Local_computePrice2 (l_arg, x, &g_price2);
		g_nbiter = comm.nf;
		g_minimum = f / 1000000.0;

		result.setDouble (x);
		result.setRetCode (ARM_OK);

		MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum - LEG1 Infos : (ccy = %s) ; (margin = %lf) ; (StartXNL = %lf) ; (EndXNL = %lf) ; (LegValue = %lf)", (const char*)ccy1, margin1, startXNL1, endXNL1, g_price1);
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapNum - LEG2 Infos : (ccy = %s) ; (margin = %lf) ; (StartXNL = %lf) ; (EndXNL = %lf) ; (LegValue = %lf)", (const char*)ccy2, x, startXNL2, endXNL2, g_price2);

	}
	else
	{
		ARM_result* l_result = (ARM_result*)(comm.p);
		result.setRetCode (ARM_KO);
		result.setMsg (l_result->getMsg ());
	}

	if(l_arg)
	{
		free (l_arg->ccy1);
		free (l_arg->ccy2);
		free (l_arg);
	}

	Local_destroySwapLegs ();

/*	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);
*/
	return result.getRetCode ();
}


long ARMLOCAL_BasisSwapAlg (double asOfDate,
							double delivery,
							double maturity,
							double margin1,
							const CCString& ccy1,
							long liborType1Id,
							long forwardCurve1Id,
							long discountCurve1Id,
							const CCString& ccy2,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							long mode,
							ARM_result& result)
{
	double adjusted = 0.0;
	double basisFlat1 = 0.0;
	double basisATMStart1 = 0.0;
	double basisATMEnd1 = 0.0;
	double basisATMFwd1 = 0.0;
	double basisFlat2 = 0.0;
	double basisATMStart2 = 0.0;
	double basisATMEnd2 = 0.0;
	double basisATMFwd2 = 0.0;
	double sensibility1 = 0.0;
	double sensibilityStart1 = 0.0;
	double sensibilityEnd1 = 0.0;
	double sensibilityFwd1 = 0.0;
	double sensibility2 = 0.0;
	double sensibilityStart2 = 0.0;
	double sensibilityEnd2 = 0.0;
	double sensibilityFwd2 = 0.0;
	double margin2 = 0.0;

	ARM_result C_result;

	long sensiCurve1Id;
	long sensiCurve2Id;
	long spotdays1 = 2;
	long spotdays2 = 2;
	double spotDate1 = 0.0;
	double spotDate2 = 0.0;
	double adjustedEnd1 = 0.0;
	double adjustedEnd2 = 0.0;
	long retCode = ARM_OK;


	// sensis sont toujours calcules avec discountcurve quelle que soit la methode
	// c'est la nature du discount qui change avec methode par le moyen des parametres d'appel
	// c-a-d si mode=alg alors sensi=BS (pass� en param ou construit avant l'appel) sinon sensi=Libor

	/*
	if(mode == BasisSwap_Method_Alg_Spr)
	{
		sensiCurve1Id = discountCurve1Id;
		sensiCurve2Id = discountCurve2Id;
	}
	else
	{
		sensiCurve1Id = forwardCurve1Id;
		sensiCurve2Id = forwardCurve2Id;
	}
	*/

	sensiCurve1Id = discountCurve1Id;
	sensiCurve2Id = discountCurve2Id;

	long indexFrequency1Id = ARM_ConvIrIndIdToFreqId (liborType1Id, C_result);
	long indexFrequency2Id = ARM_ConvIrIndIdToFreqId (liborType2Id, C_result);
	

	if(ccy1 == "GBP")
	{
		spotdays1 = 0;
	}
	if(ccy2 == "GBP")
	{
		spotdays2 = 0;
	}

	// Computing Sensi Fwd 2 & Basis ATM Fwd 2

	if(ARMLOCAL_NextBusinessDay (asOfDate,
								 ccy1,
								 spotdays1,
								 C_result) == ARM_OK)
	{
		spotDate1 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy1,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd1 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMStart1 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve1Id,
							   (adjustedEnd1 - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMEnd1 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	
	if (delivery >= spotDate1)
	{
		retCode = ARMLOCAL_CptBPV (asOfDate, 
								   spotDate1,
								   delivery,
								   sensiCurve1Id,
								   indexFrequency1Id,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart1 = C_result.getDouble ();
		}
	}
	else
	{
		retCode = ARMLOCAL_CptBPV (asOfDate, 
								   delivery,
								   spotDate1,
								   sensiCurve1Id,
								   indexFrequency1Id,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart1 = -1.0 * C_result.getDouble ();
		}
	}

	if (retCode == ARM_KO)
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_CptBPV (asOfDate,
						spotDate1,
						// adjustedEnd1,
						maturity, // l'ajustement est fait directement dans le swapleg
						sensiCurve1Id,
						indexFrequency1Id,
						inputDayCount,
						ccy1,
						amortizationId,
						C_result) == ARM_OK)
	{
		sensibilityEnd1 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	sensibilityFwd1 = sensibilityEnd1 - sensibilityStart1;

	basisATMFwd1 = (basisATMEnd1*sensibilityEnd1 - basisATMStart1*sensibilityStart1) / sensibilityFwd1;


	// Computing Sensi Fwd 2 & Basis ATM Fwd 2

	if(ARMLOCAL_NextBusinessDay (asOfDate,
								 ccy2,
								 spotdays2,
								 C_result) == ARM_OK)
	{
		spotDate2 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy2,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd2 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve2Id,
							   (delivery - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMStart2 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve2Id,
							   (adjustedEnd2 - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMEnd2 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	inputDayCount = KACTUAL_360;

    if(ccy2 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }
	
	if (delivery >= spotDate2)
	{
		retCode = ARMLOCAL_CptBPV (asOfDate, 
								   spotDate2,
								   delivery,
								   sensiCurve2Id,
								   indexFrequency2Id,
								   inputDayCount,
								   ccy2,
								   amortizationId,
								   C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart2 = C_result.getDouble ();
		}
	}
	else
	{
		retCode = ARMLOCAL_CptBPV (asOfDate, 
								   delivery,
								   spotDate2,
								   sensiCurve2Id,
								   indexFrequency2Id,
								   inputDayCount,
								   ccy2,
								   amortizationId,
								   C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart2 = -1.0 * C_result.getDouble ();
		}
	}

	if (retCode == ARM_KO)
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_CptBPV (asOfDate,
						spotDate2,
						// adjustedEnd2,
						maturity, // l'ajustement est fait directement dans le swapleg
						sensiCurve2Id,
						indexFrequency2Id,
						inputDayCount,
						ccy2,
						amortizationId,
						C_result) == ARM_OK)
	{
		sensibilityEnd2 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	sensibilityFwd2 = sensibilityEnd2 - sensibilityStart2;

	basisATMFwd2 = (basisATMEnd2*sensibilityEnd2 - basisATMStart2*sensibilityStart2) / sensibilityFwd2;

	// Compute Basis Fwd2 From Basis Fwd1 & ATM Basis Fwds & Sensi Fwds

	margin2 = (margin1 - basisATMFwd1) * sensibilityFwd1 / sensibilityFwd2 + basisATMFwd2;

	result.setDouble (margin2);
	result.setRetCode (ARM_OK);

	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (spotdate1 = %lf) ; (spotdate2 = %lf) ; (adjend1 = %lf) ; (adjend2 = %lf)", mode, spotDate1, spotDate2, adjustedEnd1, adjustedEnd2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSEnd1 = %lf) ; (BSStart1 = %lf) ; (BPVEnd1 = %lf) ; (BPVStart1 = %lf)", mode, basisATMEnd1, basisATMStart1, sensibilityEnd1, sensibilityStart1);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSEnd2 = %lf) ; (BSStart2 = %lf) ; (BPVEnd2 = %lf) ; (BPVStart2 = %lf)", mode, basisATMEnd2, basisATMStart2, sensibilityEnd2, sensibilityStart2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSATMF1 = %lf) ; (BSATMF2 = %lf) ; (BPVF1 = %lf) ; (BPVF2 = %lf)", mode, basisATMFwd1, basisATMFwd2, sensibilityFwd1, sensibilityFwd2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (margin1 = %lf) ; (margin2 = %lf)", mode, margin1, margin2);

	return result.getRetCode ();
}


long ARMLOCAL_BasisSwapAlgNew (double asOfDate,
							   double delivery,
							   double maturity,
							   double margin1,
							   long ccy1Id,
							   long liborType1Id,
							   long forwardCurve1Id,
							   long discountCurve1Id,
							   long ccy2Id,
							   long liborType2Id,
							   long forwardCurve2Id,
							   long discountCurve2Id,
							   long amortizationId,
							   long mode,
							   ARM_result& result)
{
	double adjusted = 0.0;
	double basisFlat1 = 0.0;
	double basisATMStart1 = 0.0;
	double basisATMEnd1 = 0.0;
	double basisATMFwd1 = 0.0;
	double basisFlat2 = 0.0;
	double basisATMStart2 = 0.0;
	double basisATMEnd2 = 0.0;
	double basisATMFwd2 = 0.0;
	double sensibility1 = 0.0;
	double sensibilityStart1 = 0.0;
	double sensibilityEnd1 = 0.0;
	double sensibilityFwd1 = 0.0;
	double sensibility2 = 0.0;
	double sensibilityStart2 = 0.0;
	double sensibilityEnd2 = 0.0;
	double sensibilityFwd2 = 0.0;
	double margin2 = 0.0;

	ARM_result C_result;

	long sensiCurve1Id;
	long sensiCurve2Id;
	long spotdays1 = 2;
	long spotdays2 = 2;
	double spotDate1 = 0.0;
	double spotDate2 = 0.0;
	double adjustedEnd1 = 0.0;
	double adjustedEnd2 = 0.0;
	long retCode = ARM_OK;

	// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
	char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id))->GetCcyName();
	CCString ccy1(tmp);

	tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(ccy2Id))->GetCcyName();
	CCString ccy2(tmp);

	// sensis sont toujours calcules avec discountcurve quelle que soit la methode
	// c'est la nature du discount qui change avec methode par le moyen des parametres d'appel
	// c-a-d si mode=alg alors sensi=BS (pass� en param ou construit avant l'appel) sinon sensi=Libor

	/*
	if(mode == BasisSwap_Method_Alg_Spr)
	{
		sensiCurve1Id = discountCurve1Id;
		sensiCurve2Id = discountCurve2Id;
	}
	else
	{
		sensiCurve1Id = forwardCurve1Id;
		sensiCurve2Id = forwardCurve2Id;
	}
	*/

	sensiCurve1Id = discountCurve1Id;
	sensiCurve2Id = discountCurve2Id;

	long indexFrequency1Id = ARM_ConvIrIndIdToFreqId (liborType1Id, C_result);
	long indexFrequency2Id = ARM_ConvIrIndIdToFreqId (liborType2Id, C_result);
	

	if(ccy1 == "GBP")
	{
		spotdays1 = 0;
	}
	if(ccy2 == "GBP")
	{
		spotdays2 = 0;
	}

	// Computing Sensi Fwd 2 & Basis ATM Fwd 2

	if(ARMLOCAL_NextBusinessDay (asOfDate,
								 ccy1,
								 spotdays1,
								 C_result) == ARM_OK)
	{
		spotDate1 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy1,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd1 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMStart1 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve1Id,
							   (adjustedEnd1 - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMEnd1 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	
	if (delivery >= spotDate1)
	{
		retCode = ARMLOCAL_CptBPVNew (asOfDate, 
									  spotDate1,
									  delivery,
									  sensiCurve1Id,
									  indexFrequency1Id,
									  inputDayCount,
									  ccy1Id,
									  amortizationId,
									  C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart1 = C_result.getDouble ();
		}
	}
	else
	{
		retCode = ARMLOCAL_CptBPVNew (asOfDate, 
									  delivery,
									  spotDate1,
									  sensiCurve1Id,
									  indexFrequency1Id,
									  inputDayCount,
									  ccy1Id,
									  amortizationId,
									  C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart1 = -1.0 * C_result.getDouble ();
		}
	}

	if (retCode == ARM_KO)
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_CptBPVNew (asOfDate,
						   spotDate1,
						   // adjustedEnd1,
						   maturity, // l'ajustement est fait directement dans le swapleg
						   sensiCurve1Id,
						   indexFrequency1Id,
						   inputDayCount,
						   ccy1Id,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibilityEnd1 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	sensibilityFwd1 = sensibilityEnd1 - sensibilityStart1;

	basisATMFwd1 = (basisATMEnd1*sensibilityEnd1 - basisATMStart1*sensibilityStart1) / sensibilityFwd1;


	// Computing Sensi Fwd 2 & Basis ATM Fwd 2

	if(ARMLOCAL_NextBusinessDay (asOfDate,
								 ccy2,
								 spotdays2,
								 C_result) == ARM_OK)
	{
		spotDate2 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy2,
								 1,
								 C_result) == ARM_OK)
	{
		adjustedEnd2 = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve2Id,
							   (delivery - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMStart2 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_DiscountYield (forwardCurve2Id,
							   (adjustedEnd2 - asOfDate) / 365.0,
							   -1,
							   C_result) == ARM_OK)
	{
		basisATMEnd2 = -1.0 * C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	inputDayCount = KACTUAL_360;

    if(ccy2 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }
	
	if (delivery >= spotDate2)
	{
		retCode = ARMLOCAL_CptBPVNew (asOfDate, 
									  spotDate2,
									  delivery,
									  sensiCurve2Id,
									  indexFrequency2Id,
									  inputDayCount,
									  ccy2Id,
									  amortizationId,
									  C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart2 = C_result.getDouble ();
		}
	}
	else
	{
		retCode = ARMLOCAL_CptBPVNew (asOfDate, 
									  delivery,
									  spotDate2,
									  sensiCurve2Id,
									  indexFrequency2Id,
									  inputDayCount,
									  ccy2Id,
									  amortizationId,
									  C_result);
		
		if (retCode == ARM_OK)
		{
			sensibilityStart2 = -1.0 * C_result.getDouble ();
		}
	}

	if (retCode == ARM_KO)
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_CptBPVNew (asOfDate,
						   spotDate2,
						   // adjustedEnd2,
						   maturity, // l'ajustement est fait directement dans le swapleg
						   sensiCurve2Id,
						   indexFrequency2Id,
						   inputDayCount,
						   ccy2Id,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibilityEnd2 = C_result.getDouble ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_BasisSwapAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	sensibilityFwd2 = sensibilityEnd2 - sensibilityStart2;

	basisATMFwd2 = (basisATMEnd2*sensibilityEnd2 - basisATMStart2*sensibilityStart2) / sensibilityFwd2;

	// Compute Basis Fwd2 From Basis Fwd1 & ATM Basis Fwds & Sensi Fwds

	margin2 = (margin1 - basisATMFwd1) * sensibilityFwd1 / sensibilityFwd2 + basisATMFwd2;

	result.setDouble (margin2);
	result.setRetCode (ARM_OK);

	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (spotdate1 = %lf) ; (spotdate2 = %lf) ; (adjend1 = %lf) ; (adjend2 = %lf)", mode, spotDate1, spotDate2, adjustedEnd1, adjustedEnd2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSEnd1 = %lf) ; (BSStart1 = %lf) ; (BPVEnd1 = %lf) ; (BPVStart1 = %lf)", mode, basisATMEnd1, basisATMStart1, sensibilityEnd1, sensibilityStart1);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSEnd2 = %lf) ; (BSStart2 = %lf) ; (BPVEnd2 = %lf) ; (BPVStart2 = %lf)", mode, basisATMEnd2, basisATMStart2, sensibilityEnd2, sensibilityStart2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (BSATMF1 = %lf) ; (BSATMF2 = %lf) ; (BPVF1 = %lf) ; (BPVF2 = %lf)", mode, basisATMFwd1, basisATMFwd2, sensibilityFwd1, sensibilityFwd2);
	MSG_printf_message (MSG_TRACE, "ARMLOCAL_BasisSwapAlg - (meth = %ld) ; (margin1 = %lf) ; (margin2 = %lf)", mode, margin1, margin2);

	return result.getRetCode ();
}


long ARMLOCAL_FRNMargin (double asOfDate,
						 double delivery,
						 double maturity,
						 const CCString& ccy1,
						 long liborType1Id,
						 long forwardCurve1Id,
						 long discountCurve1Id,
						 double facialMargin,
						 double price,
						 long frequencyId,
						 const CCString& ccy2,
						 long liborType2Id,
						 long forwardCurve2Id,
						 long discountCurve2Id,
						 long amortizationId,
						 long frequencyId2,
						 double fixing,
						 double spread,
						 long solve,
						 ARM_result& result)
{
	ARM_result C_result;
	
	double sensibility = 0.0;
	double sensibility2 = 0.0;
	double ccy1margin = 0.0;
	double margin = 0.0;


	if(fixing != 0.0)
	{
		if(ARMLOCAL_fixing (asOfDate,
							delivery,
							maturity,
							frequencyId,
							forwardCurve1Id,
							ccy1,
							fixing,
							facialMargin,
							C_result) == ARM_OK)
		{
			price -= C_result.getDouble ();
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	double dfs = 1.0;

	if(ARMLOCAL_DiscountPrice (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}

	dfs = result.getDouble ();

    int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	// Modif du 14/06/01 if(ARM_CptBPV (asOfDate, delivery, maturity, discountCurve1Id, frequencyId, inputDayCount, ccy1, C_result) == ARM_OK)

	if(ARMLOCAL_CptBPV (asOfDate,
						delivery,
						maturity,
						forwardCurve1Id,
						frequencyId,
						inputDayCount,
						ccy1,
						amortizationId,
						C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(sensibility == 0.0)
	{
		result.setRetCode (ARM_KO);
		result.setMsg ("ARM_ERR: sensibility is nul");
		return ARM_KO;
	}

	long retCode = ARM_KO;

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		retCode = ARM_OK;
		ccy1margin = facialMargin - (price - 100.0) * dfs / sensibility;
		margin = ccy1margin;

		if (frequencyId2 != frequencyId)
		{
			if(ARMLOCAL_CptBPV (asOfDate,
								delivery,
								maturity,
								forwardCurve1Id,
								frequencyId2,
								inputDayCount,
								ccy1,
								amortizationId,
								C_result) == ARM_OK)
			{
				sensibility2 = C_result.getDouble ();
				MSG_printf_message (MSG_TRACE, "FRN Sensibility 2 : (%lf)", sensibility2);
			}
			else
			{
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if(sensibility2 == 0.0)
			{
				result.setRetCode (ARM_KO);
				result.setMsg ("ARM_ERR: sensibility is nul");
				return ARM_KO;
			}
		
			margin *= sensibility / sensibility2;
		}
	}
	else
	{
		ccy1margin = facialMargin - (price - 100.0) * dfs / sensibility;
		retCode = ARMLOCAL_BasisSwap (asOfDate,
									  delivery,
									  maturity,
									  ccy1margin,
									  ccy1,
									  liborType1Id,
									  forwardCurve1Id,
									  discountCurve1Id,
									  ccy2,
									  liborType2Id,
									  forwardCurve2Id,
									  discountCurve2Id,
									  amortizationId,
									  solve,
									  C_result);
		
		if(retCode == ARM_OK)
		{
			margin = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "FRN Basis Conversion Margin : valoMargin(1) = (%lf) ; implyMargin(2) = (%lf) ; method = (%ld)", ccy1margin, margin, solve);
		}
		else
		{
			margin = -1.0;
			result.setMsg (C_result.getMsg ());
		}
	}

	MSG_printf_message (MSG_TRACE, "FRN Margin Details : (DPStart = %lf) ; (Sensi1 = %lf) ccy1Margin = (%lf) ; valoMargin(2) = (%lf)", dfs, sensibility, ccy1margin, margin);

	result.setRetCode (retCode);
	result.setDouble (margin);

	return retCode;
}

long ARMLOCAL_FRNMarginNew (double asOfDate,
							double delivery,
							double maturity,
							const CCString& ccy1,
							long liborType1Id,
							long forwardCurve1Id,
							long discountCurve1Id,
							double facialMargin,
							double price,
							long frequencyId,
							const CCString& ccy2,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							long frequencyId2,
							double fixing,
							double spread,
							long solve,
							ARM_result& result)
{
	ARM_result C_result;
	
	double sensibility = 0.0;
	double sensibility2 = 0.0;
	double ccy1margin = 0.0;
	double margin = 0.0;

	long retCode = ARM_KO;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;
	
	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_FRNMargin (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(fixing != 0.0)
	{
		if(ARMLOCAL_fixing (asOfDate,
							delivery,
							maturity,
							frequencyId,
							forwardCurve1Id,
							ccy1,
							fixing,
							facialMargin,
							C_result) == ARM_OK)
		{
			price -= C_result.getDouble ();
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	margin = 100. - price;

	double dfs = 1.0;
	double dfe = 1.0;
	double endDateAdjusted = 0.0;

	if(ARMLOCAL_DiscountPrice (forwardCurve1Id,
							   (delivery - asOfDate) / 365.0,
							   result) != ARM_OK)
	{
		return ARM_KO;
	}
	
	dfs = result.getDouble ();
	margin *= dfs;

	if(ARMLOCAL_ADJUSTTOBUSDATE (maturity,
								 ccy1,
								 K_MOD_FOLLOWING,
								 C_result) == ARM_OK)
	{
		endDateAdjusted = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId == -1)
	{
		if(ARMLOCAL_DiscountPrice (forwardCurve1Id,
								   (endDateAdjusted - asOfDate) / 365.0,
								   result) != ARM_OK)
		{
			return ARM_KO;
		}

		dfe = result.getDouble ();

		margin -= 100. * (dfs - dfe);
	}
	else
	{
		if(ARMLOCAL_DiscountPriceRefvalue(forwardCurve1Id,
										  amortizationId,
										  ccy1Id,
										  delivery,
										  endDateAdjusted,
										  C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			margin -= 100 * dfs - C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_FRNMarginNew (): amortissement = %lf", C_result.getDouble ());
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

	long liborfltId;
//	Modif du 28/11/01
//	if (ARMLOCAL_LIBORLEG(delivery,maturity,liborType1Id,K_RCV,facialMargin/100.,frequencyId,frequencyId,K_ADVANCE,K_ARREARS,ccy1Id,K_ADJUSTED,10,"NULL","NULL",1,C_result) == ARM_OK)
	if (ARMLOCAL_LIBORLEG(delivery,
						  maturity,
						  liborType1Id,
						  K_RCV,
						  0,
						  facialMargin/100.,
						  frequencyId,
						  frequencyId,
						  K_ADVANCE,
						  K_ARREARS,
						  //ccy1Id,
						  false,
						  ccy1,
						  K_UNADJUSTED,
						  10,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.,
						  K_YES,-1,
						  C_result) == ARM_OK)
	{
		liborfltId = C_result.getLong ();
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (liborfltId,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_FRNPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	long ycModId;

	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result) == ARM_OK)
	{
		ycModId = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "CptBPV (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	double price1 (0.0);

	if(ARMLOCAL_ARM_Price (liborfltId,
						   ycModId,
						   C_result) == ARM_OK)
	{
		price1 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Price Details : price float leg = %lf\n", price1);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	margin += price1;

	ARMLOCAL_FreeObject (ycModId, C_result);
	ARMLOCAL_FreeObject (liborfltId, C_result);

    int inputDayCount = KACTUAL_360;

    if(ccy1 == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }

	// Modif du 14/06/01 if(ARM_CptBPV (asOfDate, delivery, maturity, discountCurve1Id, frequencyId, inputDayCount, ccy1, C_result) == ARM_OK)
	if(ARMLOCAL_CptBPVNew (asOfDate,
						   delivery,
						   maturity,
						   forwardCurve1Id,
						   frequencyId2,
						   inputDayCount,
						   ccy1Id,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "FRN Sensibility : (%lf)", sensibility);
	}
	else
	{
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(sensibility == 0.0)
	{
		result.setRetCode (ARM_KO);
		result.setMsg ("ARM_ERR: sensibility is nul");
		return ARM_KO;
	}

	margin /= sensibility;

	ccy1margin = margin;

/*	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		retCode = ARM_OK;
		ccy1margin = facialMargin - (price - 100.0) * dfs / sensibility;
		margin = ccy1margin;

		if (frequencyId2 != frequencyId)
		{
			if(ARMLOCAL_CptBPVNew (asOfDate, delivery, maturity, forwardCurve1Id, frequencyId2, inputDayCount, ccy1Id, amortizationId, C_result) == ARM_OK)
			{
				sensibility2 = C_result.getDouble ();
				MSG_printf_message (MSG_TRACE, "FRN Sensibility 2 : (%lf)", sensibility2);
			}
			else
			{
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if(sensibility2 == 0.0)
			{
				result.setRetCode (ARM_KO);
				result.setMsg ("ARM_ERR: sensibility is nul");
				return ARM_KO;
			}

			margin *= sensibility / sensibility2;
		}

	}
	else
	{
*/
	if ((ccy1 != ccy2) && (!(ccy2 == "NONE")))
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		ccy1margin = margin;
		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 maturity,
										 margin,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 amortizationId,
										 solve,
										 C_result);

		if(retCode == ARM_OK)
		{
			margin = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "FRN Basis Conversion Margin : valoMargin(1) = (%lf) ; implyMargin(2) = (%lf) ; method = (%ld)", ccy1margin, margin, solve);
		}
		else
		{
			margin = -1.0;
			result.setMsg (C_result.getMsg ());
		}
	}

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	MSG_printf_message (MSG_TRACE, "FRN Margin Details : (DPStart = %lf) ; (Sensi1 = %lf) ; endDateAdj = %lf ; dfe = %lf ; ccy1Margin = (%lf) ; valoMargin(2) = (%lf)", dfs, sensibility, endDateAdjusted, dfe, ccy1margin, margin);

	result.setRetCode (retCode);
	result.setDouble (margin);

	return retCode;
}




long ARMLOCAL_LiborAssetSwapMarginAlgNew (long ycmodId,
										  long zcId,
										  long fixedLeg1Id,
										  long fixedLeg2Id,
										  long fixedLeg3Id,
										  double startDate,
										  double delivery,
										  double endDate,
										  long cpnType,
										  double fixedRate,
										  long fixDayCountId,
										  long fixReceiveOrPayId,
										  long fixFrequencyId,
										  long indexFrequencyId,
										  long fixIntRuleId,
										  double price,
										  long discountCcyId,
										  long amortizationId,
										  double redemptionPrice,
										  ARM_result& result)
{
	ARM_result C_result;
	long retCode = ARM_OK;

	double prevCpnDate = 0.0;
	double endDateAdjusted = 0.0;

	double sensibility = 0.0;
	double margin = -price;
	
//	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
//	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
//	long fixedLeg3Id = ARM_NULL_OBJECT_ID;
//	long discountCcyId = ARM_NULL_OBJECT_ID;
//	long modelId = ARM_NULL_OBJECT_ID;

	// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
//	char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId))->GetCcyName();
//	CCString discountCcy(tmp);
	
	ARM_Currency* ccy = (ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId);

	char* tmp = ccy->GetCcyName();
	CCString discountCcy(tmp);

    ARM_INDEX_TYPE K_INDEX_TYPE = GetDefaultIndexFromCurrency(ccy->GetCcyName());
    char* payCalName   = ccy->GetPayCalName(K_INDEX_TYPE);
	CCString cPayCalName(payCalName);

	if (payCalName)
		free(payCalName);
	payCalName = NULL;


	if(ARMLOCAL_ARM_Accrued (fixedLeg3Id,
							 delivery,
							 ycmodId,
							 C_result) == ARM_OK)
	{
		margin -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): accrued = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}
	
	
	if((ARMLOCAL_DiscountPrice (zcId,
								(delivery - startDate) / 365.0,
								C_result) == ARM_OK)
								&&
								(C_result.getDouble () != 0.0))
	{
		margin *= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): correction J-3 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLeg1Id,
						   ycmodId,
						   C_result) == ARM_OK)
	{
		margin += C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 1 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ARM_Price (fixedLeg2Id,
						   ycmodId,
						   C_result) == ARM_OK)
	{
		margin -= C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 2 = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_ADJUSTTOBUSDATE (endDate,
								 cPayCalName,
								 K_MOD_FOLLOWING,
								 C_result) == ARM_OK)
	{
		endDateAdjusted = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId == -1)
	{
		if(ARMLOCAL_DiscountPrice (zcId,
								   (endDateAdjusted - startDate) / 365.0,
								   C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			margin += C_result.getDouble () * redemptionPrice;
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble () * redemptionPrice);
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}
	else
	{
		if(ARMLOCAL_DiscountPriceRefvalue(zcId,
										  amortizationId,
										  discountCcyId,
										  startDate,
										  endDateAdjusted,
										  C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			margin += C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): price 3 = %lf", C_result.getDouble ());
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}

/*    int inputDayCount = KACTUAL_360;

    if(discountCcy == "GBP")
    {
       inputDayCount = KACTUAL_365;
    }
*/
	if(ARMLOCAL_CptBPVNew (startDate,
						   delivery,
						   endDate,
						   zcId,
						   indexFrequencyId,
						   ccy->GetLiborIndexDayCount(),
						   discountCcyId,
						   amortizationId,
						   C_result) == ARM_OK)
	{
		sensibility = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): sensibility = %lf", C_result.getDouble ());
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	margin /= sensibility;


//	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
//	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
//	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);
//	ARMLOCAL_FreeObject (discountCcyId, C_result);
//	ARMLOCAL_FreeObject (modelId, C_result);
	
	result.setRetCode (retCode);
	result.setDouble (margin);

	return retCode;
}



long ARMLOCAL_ASWMargin (double bondMaturity,
						 double bondCoupon,
						 long bondFrequency,
						 long bondBase,
						 double bondPrice,
						 double bondRedemptionPrice,
						 double asOfDate,
						 double delivery,
						 long fixDecompFrequency,
						 long floatResetFreq,
						 long floatPayFreq,
						 const CCString& ccy1,
						 long liborType1Id,
						 long forwardCurve1Id,
						 long discountCurve1Id,
						 const CCString& ccy2,
						 long liborType2Id,
						 long forwardCurve2Id,
						 long discountCurve2Id,
						 long amortizationId,
						 long solve,
						 long viewFlag,
						 const CCString& id,
						 ARM_result& result)
{
	ARM_result C_result;
	long ycModId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	
	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	long fixLeg1Id = ARM_NULL_OBJECT_ID;
	long fixLeg2Id = ARM_NULL_OBJECT_ID;

	if (amortizationId != -1)
	{
/*		refval = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(amortizationId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval, ARM_REFERENCE_VALUE) == 0)
		{
			ARMLOCAL_FreeObject (ccy1Id, C_result);

			result.setMsg ("ARM_ERR: amortization is not of a good type");
			return ARM_KO;
		}*/
	}
	
	double startDate = 27659.0;	// 22/09/1975!
	
	if(ARMLOCAL_FIXEDLEG (startDate,
						  bondMaturity,
						  K_RCV,
						  0,
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  fixDecompFrequency,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixLeg1Id) == ARM_OK)
	{
		if(fixLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixLeg1Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): 3");

	if(ARMLOCAL_FIXEDLEG (startDate,
						  bondMaturity,
						  K_RCV,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  fixDecompFrequency,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixLeg2Id) == ARM_OK)
	{
		if(fixLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixLeg2Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		oneCurrency = true;
	}
	
	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycModId) == ARM_OK)
	{
		if(ycModId == ARM_NULL_OBJECT_ID)
		{
			ycModId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): 4");

	double accrued1, accrued2;
	if(ARMLOCAL_ARM_Accrued (fixLeg1Id,
							 delivery,
							 ycModId,
							 C_result) == ARM_OK)
	{
		accrued1 = C_result.getDouble ();
		// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): accrued1 = %lf", accrued1);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}
	if(ARMLOCAL_ARM_Accrued (fixLeg2Id,
							 delivery,
							 ycModId,
							 C_result) == ARM_OK)
	{
		accrued2 = C_result.getDouble ();
		// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): accrued2 = %lf", accrued2);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}
	double supplFee = accrued1 - accrued2;

	if (supplFee != 0.0)
	{
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMargin - Correction de coupon : (Accr1 = %lf) ; (Accr2 = %lf)", accrued1, accrued2);
	}

	ARMLOCAL_FreeObject (fixLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixLeg2Id, C_result);
	
	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double margin1 = 0.0;
	long asw_solve = 0; // m�thode analytique (en sensis)

//	en attendant la correction du bug de ARM_liborAssetSwapMargin sur la m�thode it�rative quand 
//	previous coupon se situe entre AsOf et Delivery
//
	if((solve == BasisSwap_Method_Num) || (solve == BasisSwap_Method_Num_Spr))
	{
		asw_solve = 1; // m�thode it�rative
	}
//

	retCode = ARMLOCAL_LiborAssetSwapMargin (ycModId,
											 asOfDate,
											 bondMaturity,
											 bondCoupon,
											 fixReceiveOrPay,
											 K30_360,
											 bondFrequency,
											 fixDecompFrequency,
											 fixPayTiming,
											 fixIntRule,
											 liborType1Id,
											 spread,
											 floatResetFreq,
											 floatPayFreq,
											 assetGap,
											 viewFlag,
											 bondPrice,
											 ccy1,
											 bondRedemptionPrice,
											 supplFee,
											 asw_solve,
											 id,
											 C_result);

	if(retCode == ARM_OK)
	{
		margin1 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapMargin (): (margin = %lf)", margin1);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	// Si asw_solve = 0 et que la devise1 est du GBP alors il faut faire une correction 
	// de base de sensibilit� (en dur � KACTUAL_360 dans la fonction liborassetswapmargin !)

    int inputDayCount = KACTUAL_360;
	double sensi_std = 0.0;
	double sensi_GBP = 0.0;

    if((ccy1 == "GBP") && (asw_solve == 0))
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   bondMaturity,
								   forwardCurve1Id,
								   floatPayFreq,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);
		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPV (asOfDate,
								   delivery,
								   bondMaturity,
								   forwardCurve1Id,
								   floatPayFreq,
								   inputDayCount,
								   ccy1,
								   amortizationId,
								   C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		margin1 = margin1 * (sensi_std / sensi_GBP);
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg () - Correction sensi : (sensiSTD = %lf) ; (sensiGBP = %lf) ; (margin corrected = %lf)", sensi_std, sensi_GBP, margin1);
	}


	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): 6");

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): (margin = %lf)", margin1);

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwap (asOfDate,
									  delivery,
									  bondMaturity,
									  margin1,
									  ccy1,
									  liborType1Id,
									  forwardCurve1Id,
									  discountCurve1Id,
									  ccy2,
									  liborType2Id,
									  forwardCurve2Id,
									  discountCurve2Id,
									  amortizationId,
									  solve,
									  result);
		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (valoMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", margin1, C_result.getDouble (), solve);
	}
	else
	{
		// result.setDouble (C_result.getDouble ());
		result.setDouble (margin1);
		result.setRetCode (ARM_OK);
	}
	
	ARMLOCAL_FreeObject (ycModId, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);

	return retCode;
}



long ARMLOCAL_ASWMarginNew (double bondMaturity,
							long cpnType,
							double bondCoupon,
							long bondFrequency,
							long bondBase,
							double bondPrice,
							double bondRedemptionPrice,
							double asOfDate,
							double delivery,
							long fixDecompFrequency,
							long floatResetFreq,
							long floatPayFreq,
							const CCString& ccy1,
							long liborType1Id,
							long forwardCurve1Id,
							long discountCurve1Id,
							const CCString& ccy2,
							long liborType2Id,
							long forwardCurve2Id,
							long discountCurve2Id,
							long amortizationId,
							long solve,
							double minValue,
							double maxValue,
							ARM_result& result)
{
	ARM_result C_result;
	long ycmodId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;

	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1,
								forwardCurve1Id,
								1.0,
								KACTUAL_365,
								C_result,
								ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1,
								   C_result,
								   ccy1Id);
	}

	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if((ccy2 == "NONE") || (ccy1 == ccy2))
	{
		oneCurrency = true;
	}
	else
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}	
	
	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double margin1 = 0.0;
	long asw_solve = 0; // m�thode analytique (en sensis)

	double prevCpnDate = 0.0;
	
	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycmodId) == ARM_OK)
	{
		if(ycmodId == ARM_NULL_OBJECT_ID)
		{
			ycmodId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  K_ANNUAL,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg3Id) == ARM_OK)
	{
		if(fixedLeg3Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg3Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg3Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_PrevCpnDate (delivery,
								 bondMaturity,
								 bondFrequency,
								 fixIntRule,
								 ccy1,
								 C_result) == ARM_OK)
	{
		prevCpnDate = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): prevCpnDate = %lf", prevCpnDate);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  prevCpnDate,
						  fixReceiveOrPay,
						  cpnType,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

//	en attendant la correction du bug de ARM_liborAssetSwapMargin sur la m�thode it�rative quand 
//	previous coupon se situe entre AsOf et Delivery
//
	if((solve == BasisSwap_Method_Num) || (solve == BasisSwap_Method_Num_Spr))
	{
		retCode = ARMLOCAL_LiborAssetSwapMarginNumNew (ycmodId,
													   forwardCurve1Id,
													   fixedLeg1Id,
													   fixedLeg2Id,
													   fixedLeg3Id,
													   asOfDate,
													   delivery,
													   bondMaturity,
													   floatPayFreq,
													   bondPrice,
													   ccy1Id,
													   amortizationId,
													   bondRedemptionPrice,
													   minValue,
													   maxValue,
													   C_result);

	}
	else
	{
		// Nouvelle m�thode pour enlever le lien vers liborassetswapmargin
		retCode = ARMLOCAL_LiborAssetSwapMarginAlgNew (ycmodId,
													   forwardCurve1Id,
													   fixedLeg1Id,
													   fixedLeg2Id,
													   fixedLeg3Id,
													   asOfDate,
													   delivery,
													   bondMaturity,
													   cpnType,
													   bondCoupon,
													   bondBase,
													   K_RCV,
													   bondFrequency,
													   floatPayFreq,
													   K_UNADJUSTED,
													   bondPrice,
													   ccy1Id,
													   amortizationId,
													   bondRedemptionPrice,
													   C_result);
	}

	if(retCode == ARM_OK)
	{
		margin1 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): (margin = %lf)", margin1);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	// Si asw_solve = 0 et que la devise1 est du GBP alors il faut faire une correction 
	// de base de sensibilit� (en dur � KACTUAL_360 dans la fonction liborassetswapmargin !)

/*    int inputDayCount = KACTUAL_360;
	double sensi_std = 0.0;
	double sensi_GBP = 0.0;

    if((ccy1 == "GBP") && (asw_solve == 0))
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
				   delivery,
				   bondMaturity,
				   forwardCurve1Id,
				   floatPayFreq,
				   inputDayCount,
				   ccy1Id,
				   amortizationId,
				   C_result);
		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
				   delivery,
				   bondMaturity,
				   forwardCurve1Id,
				   floatPayFreq,
				   inputDayCount,
				   ccy1Id,
				   amortizationId,
				   C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		margin1 = margin1 * (sensi_std / sensi_GBP);
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew () - Correction sensi : (sensiSTD = %lf) ; (sensiGBP = %lf) ; (margin corrected = %lf)", sensi_std, sensi_GBP, margin1);
	}

*/
	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): 6");

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): (margin = %lf)", margin1);

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 bondMaturity,
										 margin1,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 amortizationId,
										 solve,
										 result);

		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (valoMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", margin1, C_result.getDouble (), solve);
	}
	else
	{
		// result.setDouble (C_result.getDouble ());
		result.setDouble (margin1);
		result.setRetCode (ARM_OK);
	}
	
	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);

	ARMLOCAL_FreeObject (ycmodId, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return retCode;
}



long ARMLOCAL_BondASWMargin (long bondId,
							 double bondPrice,
							 double asOfDate,
							 double delivery,
							 long floatResetFreq,
							 long floatPayFreq,
							 const CCString& ccy1,
							 long forwardCurve1Id,
							 long discountCurve1Id,
							 const CCString& ccy2,
							 long liborType2Id,
							 long forwardCurve2Id,
							 long discountCurve2Id,
							 long amortizationId,
							 long solve,
							 double minValue,
							 double maxValue,
							 ARM_result& result)
{
	ARM_result C_result;
	long ycmodId = ARM_NULL_OBJECT_ID;
	bool oneCurrency = false;
	long retCode = ARM_KO;

	long ccy1Id = ARM_NULL_OBJECT_ID;
	long ccy2Id = ARM_NULL_OBJECT_ID;

	long fixedLeg1Id = ARM_NULL_OBJECT_ID;
	long fixedLeg2Id = ARM_NULL_OBJECT_ID;
	long fixedLeg3Id = ARM_NULL_OBJECT_ID;

	ARM_Bond* bond = NULL;
	ARM_Currency* Ccy1Obj = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	bond = (ARM_Bond*) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ARM_BOND) == 0) 
	{
		result.setMsg ("ARM_ERR: bond is not of a good type");
		return ARM_KO;
	}

	if(ccy1 == "GBP")
	{
		retCode = ARMLOCAL_CCY (ccy1, forwardCurve1Id, 1, 2, C_result, ccy1Id);
	}
	else
	{
		retCode = ARMLOCAL_ISOCCY (ccy1, C_result, ccy1Id);
	}
	if(retCode == ARM_OK)
	{
		ccy1Id = C_result.getLong ();
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMargin (): %s", C_result.getMsg ());
		g_result->setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	Ccy1Obj = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id);

	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Ccy1Obj, ARM_CURRENCY) == 0)
	{
		result.setMsg ("ARM_ERR: ccy1 is not of a good type");
		return ARM_KO;
	}

	char* ccy1name = ccy1.GetStr();

	long liborType1Id = GetDefaultIndexFromCurrency(ccy1name);

	if (ccy1name)
		delete ccy1name;
	ccy1name = NULL;

	if(ccy1 == ccy2)
	{
		oneCurrency = true;
	}
	else
	{
		if(ccy2 == "GBP")
		{
			retCode = ARMLOCAL_CCY (ccy2,
									forwardCurve2Id,
									1.0,
									KACTUAL_365,
									C_result,
									ccy2Id);
		}
		else
		{
			retCode = ARMLOCAL_ISOCCY (ccy2,
									   C_result,
									   ccy2Id);
		}
		if(retCode == ARM_OK)
		{
			ccy2Id = C_result.getLong ();
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			g_result->setMsg (C_result.getMsg ());
			return ARM_KO;
		}
	}	
	
	char strDate[11];
	bond->GetMaturity().JulianToStrDate(strDate);
	double bondMaturity = Local_ARMDATE2XLDATE(strDate);
	double bondCoupon = bond->GetCoupon();
	double bondBase = bond->GetDayCount();
	double bondFrequency = bond->GetFrequency();
	double bondRedemptionPrice = bond->GetRedemptionValue();

	// Finesse du Swap (dixit AR)
	long fixDecompFrequency = 1;
	
	long fixReceiveOrPay = K_RCV;
	long fixPayTiming = K_ARREARS;
	long fixIntRule = K_UNADJUSTED;
	double spread = 0.0;
	long assetGap = 3;
	double margin1 = 0.0;
	long asw_solve = 0; // m�thode analytique (en sensis)

	double prevCpnDate = 0.0;
	
	if(ARMLOCAL_ycmod (forwardCurve1Id,
					   ARM_NULL_OBJECT,
					   C_result,
					   ycmodId) == ARM_OK)
	{
		if(ycmodId == ARM_NULL_OBJECT_ID)
		{
			ycmodId = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  bondBase,
						  bondFrequency,
						  K_ANNUAL,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg3Id) == ARM_OK)
	{
		if(fixedLeg3Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg3Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg3Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  bondMaturity,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg1Id) == ARM_OK)
	{
		if(fixedLeg1Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg1Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg1Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

	if(ARMLOCAL_ARM_PrevCpnDate (delivery,
								 bondMaturity,
								 bondFrequency,
								 fixIntRule,
								 ccy1,
								 C_result) == ARM_OK)
	{
		prevCpnDate = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_LiborAssetSwapPriceAlg (): prevCpnDate = %lf", prevCpnDate);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if(ARMLOCAL_FIXEDLEG (27659.0,
						  prevCpnDate,
						  fixReceiveOrPay,
						  0,
						  bondCoupon,
						  K30_360,
						  bondFrequency,
						  1,
						  K_ARREARS,
						  K_UNADJUSTED,
						  K_SHORTSTART,
						  //ccy1Id,
						  false,
						  ccy1,
						  "NULL",
						  K_NX_NONE,
						  -1.,
						  K_YES,
						  -1.,
						  C_result,
						  fixedLeg2Id) == ARM_OK)
	{
		if(fixedLeg2Id == ARM_NULL_OBJECT_ID)
		{
			fixedLeg2Id = C_result.getLong ();
		}
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
		result.setRetCode (ARM_KO);
		result.setMsg (C_result.getMsg ());
		return ARM_KO;
	}

	if (amortizationId != -1)
	{
		if(ARMLOCAL_SetNotional (fixedLeg2Id,
								 amortizationId,
								 100.,
								 C_result) != ARM_OK)
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWPrice (): %s", C_result.getMsg ());
			result.setMsg (C_result.getMsg ());
			result.setRetCode (ARM_KO);
			return ARM_KO;
		}
	}

//	en attendant la correction du bug de ARM_liborAssetSwapMargin sur la m�thode it�rative quand 
//	previous coupon se situe entre AsOf et Delivery
//
	if((solve == BasisSwap_Method_Num) || (solve == BasisSwap_Method_Num_Spr))
	{
		retCode = ARMLOCAL_LiborAssetSwapMarginNumNew (ycmodId,
													   forwardCurve1Id,
													   fixedLeg1Id,
													   fixedLeg2Id,
													   fixedLeg3Id,
													   asOfDate,
													   delivery,
													   bondMaturity,
													   floatPayFreq,
													   bondPrice,
													   ccy1Id,
													   amortizationId,
													   bondRedemptionPrice,
													   minValue,
													   maxValue,
													   C_result);

	}
	else
	{
		// Nouvelle m�thode pour enlever le lien vers liborassetswapmargin
		retCode = ARMLOCAL_LiborAssetSwapMarginAlgNew (ycmodId,
													   forwardCurve1Id,
													   fixedLeg1Id,
													   fixedLeg2Id,
													   fixedLeg3Id,
													   asOfDate,
													   delivery,
													   bondMaturity,
													   0,
													   bondCoupon,
													   bondBase,
													   K_RCV,
													   bondFrequency,
													   floatPayFreq,
													   K_UNADJUSTED,
													   bondPrice,
													   ccy1Id,
													   amortizationId,
													   bondRedemptionPrice,
													   C_result);
	}

	if(retCode == ARM_OK)
	{
		margin1 = C_result.getDouble ();
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): (margin = %lf)", margin1);
	}
	else
	{
		MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
		result.setMsg (C_result.getMsg ());
		result.setRetCode (ARM_KO);
		return ARM_KO;
	}

	// Si asw_solve = 0 et que la devise1 est du GBP alors il faut faire une correction 
	// de base de sensibilit� (en dur � KACTUAL_360 dans la fonction liborassetswapmargin !)

/*    int inputDayCount = KACTUAL_360;
	double sensi_std = 0.0;
	double sensi_GBP = 0.0;

    if((ccy1 == "GBP") && (asw_solve == 0))
    {
		inputDayCount = KACTUAL_360;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
				   delivery,
				   bondMaturity,
				   forwardCurve1Id,
				   floatPayFreq,
				   inputDayCount,
				   ccy1Id,
				   amortizationId,
				   C_result);
		if (retCode == ARM_OK)
		{
			sensi_std = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): sensi_std = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		inputDayCount = KACTUAL_365;
    
		retCode = ARMLOCAL_CptBPVNew (asOfDate,
				   delivery,
				   bondMaturity,
				   forwardCurve1Id,
				   floatPayFreq,
				   inputDayCount,
				   ccy1Id,
				   amortizationId,
				   C_result);
		if (retCode == ARM_OK)
		{
			sensi_GBP = C_result.getDouble ();
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew (): sensi_GBP = %lf", C_result.getDouble ());
		}
		else
		{
			MSG_printf_message (MSG_ERROR, "ARMLOCAL_ASWMarginNew (): %s", C_result.getMsg ());
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		margin1 = margin1 * (sensi_std / sensi_GBP);
		MSG_printf_message (MSG_TRACE, "ARMLOCAL_ASWMarginNew () - Correction sensi : (sensiSTD = %lf) ; (sensiGBP = %lf) ; (margin corrected = %lf)", sensi_std, sensi_GBP, margin1);
	}
*/

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): 6");

	// MSG_printf_message (MSG_TRACE, "ARM_ASWMargin (): (margin = %lf)", margin1);

	if(oneCurrency == false)
	{
		retCode = ARMLOCAL_BasisSwapNew (asOfDate,
										 delivery,
										 bondMaturity,
										 margin1,
										 ccy1Id,
										 liborType1Id,
										 forwardCurve1Id,
										 discountCurve1Id,
										 ccy2Id,
										 liborType2Id,
										 forwardCurve2Id,
										 discountCurve2Id,
										 amortizationId,
										 solve,
										 result);
		MSG_printf_message (MSG_TRACE, "ASW Basis Conversion : (valoMargin(1) = %lf) ; (implyMargin(2) = %lf) ; (method = %ld)", margin1, C_result.getDouble (), solve);
	}
	else
	{
		// result.setDouble (C_result.getDouble ());
		result.setDouble (margin1);
		result.setRetCode (ARM_OK);
	}
	
	ARMLOCAL_FreeObject (fixedLeg1Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg2Id, C_result);
	ARMLOCAL_FreeObject (fixedLeg3Id, C_result);

	ARMLOCAL_FreeObject (ycmodId, C_result);

	ARMLOCAL_FreeObject (ccy1Id, C_result);
	ARMLOCAL_FreeObject (ccy2Id, C_result);

	return retCode;
}



long ARMLOCAL_DiscountPriceRefvalue(long zcId,
									long refvalId,
									long discountCcyId,
									double startDate,
									double endDate,
									ARM_result& result)
{
	ARM_result C_result;

	double price = 0;
	ARM_ReferenceValue* refval = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{		
		// Ne pas d�truire le tmp, il est d�truit � la destruction de ARM_Currency
		char* tmp = ((ARM_Currency*)LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId))->GetCcyName();
		CCString discountCcy(tmp);

		refval = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: RefValue is not of a good type");
			return ARM_KO;
		}

		ARM_Vector * refDates = refval->GetDiscreteDates();
		ARM_Vector * refValues = refval->GetDiscreteValues();
		int i = refDates->GetSize() - 1;
		double dAmort, dFrac;
		double datei (0.);
		long lTest = 0;

		while ( (i >= 0 ) && (lTest == 0 ) )
		{
			char buf[30];
			ARM_Date aDate (refDates->Elt(i));
			aDate.JulianToStrDate(buf);

			datei = Local_ARMDATE2XLDATE(buf);

			if(ARMLOCAL_ADJUSTTOBUSDATE (datei,
										 discountCcy,
										 K_MOD_FOLLOWING,
										 C_result) == ARM_OK)
			{
				datei = Local_ARMDATE2XLDATE(C_result.getString ());
			}
			else
			{
				MSG_printf_message (MSG_ERROR, "ARMLOCAL_LiborAssetSwapPriceAlg (): %s", C_result.getMsg ());
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if (datei <= endDate)
				lTest = 1;
			else
				i--;
		}

		dAmort = refValues->Elt(i);
		dFrac = (endDate - startDate) / 365.0;

		if(ARMLOCAL_DiscountPrice (zcId,
								   dFrac,
								   C_result) == ARM_OK)
		{
			// price += C_result.getDouble () * redemptionPrice * -100.0;
			price += C_result.getDouble () * dAmort;
			MSG_printf_message (MSG_TRACE, "ARMLOCAL_DiscountPriceRefvalue (): amortissement %lf : %lf, price = %lf", endDate, dAmort, C_result.getDouble ()*dAmort );
		}
		else
		{
			result.setRetCode (ARM_KO);
			result.setMsg (C_result.getMsg ());
			return ARM_KO;
		}

		while ( (i >  0) && (datei > startDate) )
		{
			char buf[30];
			ARM_Date aDate (refDates->Elt(i-1));
			aDate.JulianToStrDate(buf);

			datei = Local_ARMDATE2XLDATE(buf);
			if(ARMLOCAL_ADJUSTTOBUSDATE (datei,
										 discountCcy,
										 1,
										 C_result) == ARM_OK)
			{
				datei = Local_ARMDATE2XLDATE(C_result.getString ());
			}
			else
			{
				MSG_printf_message (MSG_ERROR, "ARMLOCAL_DiscountPriceRefvalue (): %s", C_result.getMsg ());
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}

			if (datei <= startDate)
				break;

			dAmort = refValues->Elt(i-1) - refValues->Elt(i);
			dFrac = (datei - startDate) / 365.0;

			if(ARMLOCAL_DiscountPrice (zcId,
									   dFrac,
									   C_result) == ARM_OK)
			{
				// price += C_result.getDouble () * redemptionPrice * -100.0;
				price += C_result.getDouble () * dAmort;
				MSG_printf_message (MSG_TRACE, "ARMLOCAL_DiscountPriceRefvalue (): amortissement %lf: %lf, price = %lf", datei, dAmort, C_result.getDouble()*dAmort);
			}
			else
			{
				result.setRetCode (ARM_KO);
				result.setMsg (C_result.getMsg ());
				return ARM_KO;
			}
			i--;
		}

		result.setDouble(price);
		return ARM_OK;
	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
        ARM_RESULT();
    }
}