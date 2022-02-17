
/*! \file ARM_local_gp_genericaddin.cpp
 *
 *  \brief file for the generic addin
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2006
 */


#include "ARM_local_gp_genericaddin.h"

#include "ARM_local_gp_calculators.h"
#include "ARM_local_gp_inflation.h"

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"
#include "ARM_local_class.h"
#include "ARM_local_gp_model.h"
#include "ARM_local_xxxproject.h"
#include "ARM_local_gp_calib.h"

#include "ARM_local_gp_closedforms.h"
#include <algorithm>

/// gphelp
#include <GP_Help\gphelp\crmcookies.h>

#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
using ARM::stringTrim;
using ARM::stringToUpper;
using ARM::ARM_CRMCookies;



ARM_GenericParamValue DefZero(0.0);
ARM_GenericParamValue DefOne(1.0);
ARM_GenericParamValue DefNullObject((long)ARM_NULL_OBJECT); //is defined as #define ARM_NULL_OBJECT -11111

ARM_GenericParamValue DefNbSimul(50000.0);
ARM_GenericParamValue DefNbFactors(6.0);
ARM_GenericParamValue DefTimeSteps(1000.0);
ARM_GenericParamValue DefSpaceSteps(901.0);
ARM_GenericParamValue DefNbStdDevs(6.0);
ARM_GenericParamValue DefEpsilon(0.00001);
ARM_GenericParamValue DefNbGaussLeg(60.0);
ARM_GenericParamValue DefPtGaussLeg(5.0);
ARM_GenericParamValue DefResetGap(-2.0);
ARM_GenericParamValue DefParam(-999.0);
ARM_GenericParamValue DefObjParam((long)-1);

ARM_GenericParamValue DefUnknown("Unknown");
ARM_GenericParamValue DefNone("NONE");
ARM_GenericParamValue DefFreq("S");
ARM_GenericParamValue DefIsNothing("IsNothing");
ARM_GenericParamValue DefYes("Y");
ARM_GenericParamValue DefNo("N");
ARM_GenericParamValue DefNO("NO");

ARM_GenericParamValue DefPayRec("PAY");
ARM_GenericParamValue DefPayTiming("ARR");
ARM_GenericParamValue DefResetTiming("ADV");
ARM_GenericParamValue DefDayCount("A360");
ARM_GenericParamValue DefStub("SS");
ARM_GenericParamValue DefFwdRule("MF");
ARM_GenericParamValue DefIntRule("UNADJ");
ARM_GenericParamValue DefRedeemType("STANDARD");
ARM_GenericParamValue DefTerminalZc("TerminalZc");
ARM_GenericParamValue DefRandGenType1("Sobol");
ARM_GenericParamValue DefRandGenAlgo1("InvNormCum");
ARM_GenericParamValue DefRandGenType2("NR_Ran2");
ARM_GenericParamValue DefRandGenAlgo2("InvNormCum");
ARM_GenericParamValue DefCorrelType("Beta");
ARM_GenericParamValue DefTARNFXModel("1IRFX");
ARM_GenericParamValue DefInfInterType("CPILINEAR");
ARM_GenericParamValue DefModel("GAUSS");
ARM_GenericParamValue DefCurrency(ARM_DEFAULT_COUNTRY);
ARM_GenericParamValue DefSabrFlag("SABR_IMPLNVOL");
ARM_GenericParamValue DefCorrType("CorrelMatrix");
ARM_GenericParamValue DefCallPut("C");
ARM_GenericParamValue DefMinMax("Max");
ARM_GenericParamValue DefDigitType("CENTRED");
ARM_GenericParamValue DefTARNFXChoise("Worst");
ARM_GenericParamValue DefOptionType("Vanilla");
ARM_GenericParamValue DefDom("Dom");
ARM_GenericParamValue DefPayoffType("TARNFX");

ARM_PRDCCalculator_CreateFunctor			ThePRDCCalculator_CreateFunctor;
ARM_PRDCCalculator_InitFunctor				ThePRDCCalculator_InitFunctor;
ARM_TARNFXCalculator_CreateFunctor			TheTARNFXCalculator_CreateFunctor;
ARM_TARNFXCalculator_InitFunctor			TheTARNFXCalculator_InitFunctor;
ARM_TARNCalculatorIndian_CreateFunctor		TheTARNCalculatorIndian_CreateFunctor;
ARM_CCSCalculator_CreateFunctor				TheCCSCalculator_CreateFunctor;
ARM_CCSCalculator_InitFunctor				TheCCSCalculator_InitFunctor;

ARM_FXVanillaCalculator_CreateFunctor	    TheFXVanillaCalculator_CreateFunctor;
ARM_FXVanillaCalculator_InitFunctor		    TheFXVanillaCalculator_InitFunctor;
ARM_FXRACalculator_CreateFunctor			TheFXRACalculator_CreateFunctor;
ARM_FXRACalculator_InitFunctor				TheFXRACalculator_InitFunctor;

ARM_ParamsMixtureFx_CreateFunctor			TheParamsMixtureFx_CreateFunctor;
ARM_MixtureModelFx_CreateWithParamsFunctor	TheMixtureModelFx_CreateWithParamsFunctor;
ARM_ModelBumpParams_CreateFunctor			TheModelBumpParams_CreateFunctor;
ARM_GP_XXXMktDataFromMktDataMgerFunctor		TheXXXMktDataFromMktDataMgerFunctor;
ARM_FXModelDensity_CreateFunctor			TheFXModelDensity_CreateFunctor;
ARM_QDensityFunctor_Create					TheQDensityFunctor_CreateFunctor;


ARM_SABR_SmileCalibration_ParamFunctor		TheSABRSmileCalib_CreateFunctor;
ARM_Heston_SmileCalibration_ParamFunctor	TheHestonSmileCalib_CreateFunctor;
ARM_Heston2b_SmileCalibration_ParamFunctor	TheHeston2bSmileCalib_CreateFunctor;
ARM_BiSABR_SmileCalibration_ParamFunctor	TheBiSABRSmileCalib_CreateFunctor;
XXX_MktData_ApplyScenarioFunctor			TheXXX_MktData_ApplyScenarioFunctor;
ARM_VanillaDensity_CreateFunctor			TheVanillaDensityFunctor;
ARM_SABR2B_SmileCalibration_ParamFunctor	TheSABR2BSmileCalib_CreateFunctor;
ARM_Merton_SmileCalibration_ParamFunctor	TheMertonSmileCalib_CreateFunctor;
ARM_2IRFX_ComputeTimeLagFunctor				The2IRFX_ComputeTimeLagFunctor;
ARM_2IRFX_ComputeFwdFxVolFunctor			The2IRFX_ComputeFwdFxVolFunctor;
ARM_2IRFX_ComputeFwdFxModelParamFunctor		The2IRFX_ComputeFwdFxModelParamFunctor;


ARM_Corridor_InfIr_CreateFunctor			TheCorridor_InfIr_CreateFunctor;
ARM_HybridInfIrLeg_CreateFunctor			TheHybridInfIrLeg_CreateFunctor;
ARM_HybridInfIrPayOff_CreateFunctor			TheHybridInfIrPayOff_CreateFunctor;
ARM_HybridInfIrModel_CreateFunctor			TheHybridInfIrModel_CreateFunctor;

/**************************************************************************************************/

const GenericParamStruct HybridInfIrModel_CreateParams[] = 
{
	///////////////			////////////////////	//////////////////////////////	////////////////
	//	Param Type			Param Name				ParamDescription				DefaultValue

	{	GA_STRING,			"ModelName",			"Model Name",						NULL			},
	{	GA_DOUBLE,			"Discretisation",		"Discretisation",					NULL			},
	{	GA_DOUBLE,			"Domain",				"Domain",							NULL			},
	{	GA_DOUBLE,			"Center",				"Center",							&DefZero		},
	{	GA_DOUBLE,			"Epsilon",				"Epsilon",							NULL			},

	{	GA_TERMINAL }
};

const GenericParamStruct HybridInfIrPayOff_CreateParams[] = 
{
	///////////////			////////////////////	//////////////////////////////	////////////////
	//	Param Type			Param Name				ParamDescription				DefaultValue

	{	GA_STRING,			"MainCpnName",			"Main Coupon  Index Name",			&DefNO			},
	{	GA_DOUBLE,			"MainCpnCoef",			"Main Coupon  Index Coef",			NULL			},
	{	GA_DOUBLE,			"SubCpnName",			"Sub  Coupon  Index Name",			&DefNO			},
	{	GA_DOUBLE,			"SubCpnCoef",			"Sub  Coupon  Index Coef",			NULL			},
	{	GA_STRING,			"SupCpnName",			"Sup  Coupon  Index Name",			&DefNO			},
	{	GA_OBJECT,			"SupCpnCoef",			"Sup  Coupon  Index Coef",			NULL			},
	{	GA_OBJECT,			"CstCpnCoef",			"Cst  Coupon  Index Coef",			NULL			},
	{	GA_STRING,			"MainOptName",			"Main Option  Index Name",			&DefNO			},
	{	GA_OBJECT,			"MainOptCoef",			"Main Option  Index Coef",			NULL			},
	{	GA_STRING,			"SubOptName",			"Sub  Option  Index Name",			&DefNO			},
	{	GA_OBJECT,			"SubOptCoef",			"Sub  Option  Index Coef",			NULL			},
	{	GA_STRING,			"SupOptName",			"Sup  Option  Index Name",			&DefNO			},
	{	GA_OBJECT,			"SupOptCoef",			"Sup  Option  Index Coef",			NULL			},
	{	GA_OBJECT,			"CstOptCoef",			"Cst  Option  Index Coef",			NULL			},

	{	GA_TERMINAL }
};


const GenericParamStruct HybridInfIrLeg_CreateParams[] = 
{
	///////////////			////////////////////	//////////////////////////////	////////////////
	//	Param Type			Param Name				ParamDescription				DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"As of Date",									NULL			},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",						NULL			},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",							NULL			},
	{	GA_STRING,			"MainIndex",			"MainIndex",									NULL			},
	{	GA_STRING,			"SubIndex",				"SubIndex",										&DefNo			},
	{	GA_STRING,			"SupIndex",				"SupIndex",										&DefNo			},
	{	GA_STRING,			"InfInterType",			"Inflation Interpolation Type",					&DefInfInterType},
	{	GA_STRING,			"Currency",				"Currency",										&DefCurrency	},
	{	GA_STRING,			"ResetFreq",			"Reset Frequency (A/S/Q/M)",					NULL			},
	{	GA_STRING,			"ResetTiming",			"Reset Timing (ARR/ADV)",						&DefResetTiming	},
	{	GA_DOUBLE,			"ResetGap",				"Reset Gap",									&DefResetGap	},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",								&DefCurrency	},
	{	GA_DOUBLE,			"ResetNumGap",			"Reset Numerator Gap",							NULL		},
	{	GA_DOUBLE,			"ResetDemGap",			"Reset Denominator Gap",						NULL		},
	{	GA_STRING,			"PayFreq",				"Pay Frequency (A/S/Q/M)",						NULL			},
	{	GA_STRING,			"PayTiming",			"Payment Timing (ARR/ADV)",						&DefPayTiming	},
	{	GA_DOUBLE,			"PayGap",				"Payment Gap",									&DefZero		},	
	{	GA_STRING,			"PayCal",				"Payment Calendar",								&DefCurrency	},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",					&DefIntRule		},
	{	GA_STRING,			"StubRule",				"Stub Type (SS/LS)",							&DefStub		},
	{	GA_STRING,			"FwdRule",				"Reset Rule (F/MF/P/MP)",						&DefFwdRule		},
	{	GA_STRING,			"AdjFirstRule",			"Adjust First Rule (ADJ/UNADJ)",				&DefIntRule		},
	{	GA_STRING,			"DayCount",				"Day Count (A360/A365/30/360)",					&DefDayCount	},

	{	GA_TERMINAL }
};


const GenericParamStruct Corridor_InfIr_CreateParams[] = 
{
	///////////////			////////////////////	//////////////////////////////	////////////////
	//	Param Type			Param Name				ParamDescription				DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"As of Date",									NULL			},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",						NULL			},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",							NULL			},
	{	GA_STRING,			"IrIndex",				"Rate Index",									NULL			},
	{	GA_STRING,			"InfIndex",				"Inflation Index",								NULL			},
	{	GA_STRING,			"InfInterType",			"Inflation Interpolation Type",					&DefInfInterType},
	{	GA_STRING,			"Currency",				"Currency",										&DefCurrency	},
	{	GA_STRING,			"ResetFreq",			"Reset Frequency (A/S/Q/M)",					NULL			},
	{	GA_STRING,			"ResetTiming",			"Reset Timing (ARR/ADV)",						&DefResetTiming	},
	{	GA_DOUBLE,			"ResetGap",				"Reset Gap",									&DefResetGap	},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",								&DefCurrency	},
	{	GA_DOUBLE,			"ResetNumGap",			"Reset Numerator Calendar",						&DefCurrency	},
	{	GA_DOUBLE,			"ResetDemGap",			"Reset Denominator Calendar",					&DefCurrency	},
	{	GA_STRING,			"PayFreq",				"Pay Frequency (A/S/Q/M)",						NULL			},
	{	GA_STRING,			"PayTiming",			"Payment Timing (ARR/ADV)",						&DefPayTiming	},
	{	GA_DOUBLE,			"PayGap",				"Payment Gap",									&DefZero		},	
	{	GA_STRING,			"PayCal",				"Payment Calendar",								&DefCurrency	},
	{	GA_STRING,			"IrIntRule",			"Ir Index Interest Rule (ADJ/UNADJ)",			&DefIntRule		},
	{	GA_STRING,			"InfIntRule",			"Inf IndexInterest Rule (ADJ/UNADJ)",			&DefIntRule		},
	{	GA_STRING,			"StubRule",				"Stub Type (SS/LS)",							&DefStub		},
	{	GA_STRING,			"FwdRule",				"Reset Rule (F/MF/P/MP)",						&DefFwdRule		},
	{	GA_STRING,			"AdjFirstRule",			"Adjust First Rule (ADJ/UNADJ)",				&DefIntRule		},
	{	GA_STRING,			"RecOrPay",				"Receive or Pay",								&DefPayRec		},
	{	GA_STRING,			"IrDayCount",			"Rate Index Day Count (A360/A365/30/360)",		&DefDayCount	},
	{	GA_STRING,			"InfDayCount",			"Inflation Index Day Count (A360/A365/30/360)",	&DefDayCount	},
	{	GA_STRING,			"CpnDayCount",			"Coupon Index Day Count (A360/A365/30/360)",	&DefDayCount	},
	{	GA_STRING,			"IsIrCritera",			"Corridor Critera on Ir Index (YES/NO)",		&DefYes			},	
	{	GA_STRING,			"IsModulable",			"Ratio Modulation on Index (YES/NO)",			&DefNo			},	
	{	GA_STRING,			"Model",				"Kind of Copula (NAIVE/GAUSS)",					&DefModel		},	
	{	GA_DOUBLE,			"Epsilon",				"Epsilon",										&DefEpsilon		},
	{	GA_DOUBLE,			"NbGaussLeg",			"NbGaussLeg",									&DefNbGaussLeg	},
	{	GA_DOUBLE,			"PtGaussLeg",			"PtGaussLeg",									&DefPtGaussLeg	},
	{	GA_OBJECT,			"Notional",				"Notional",										NULL			},
	{	GA_OBJECT,			"IrLeverage",			"IrLeverage",									NULL			},
	{	GA_OBJECT,			"InfLeverage",			"InfLeverage",									NULL			},
	{	GA_OBJECT,			"Constant",				"Constant",										NULL			},
	{	GA_OBJECT,			"MultipleUp",			"Multiple Up",									NULL			},
	{	GA_OBJECT,			"MultipleDown",			"Multiple Down",								NULL			},
	{	GA_OBJECT,			"RangeUp",				"Range Up",										NULL			},
	{	GA_OBJECT,			"RangeDown",			"Range Down",									NULL			},

	{	GA_TERMINAL }
};

/**************************************************************************************************/
const GenericParamStruct PRDCCalculator_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"FixEndDate",			"End Date of the fixed period",									NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",											NULL				},
	{	GA_STRING,			"CpnCcy",				"Coupon currency",												NULL				},
	{	GA_STRING,			"FgnCcy",				"Foreign currency,",											NULL				},
	{	GA_STRING,			"FundCcy",				"Funding currency",												NULL				},
	{	GA_STRING,			"CpnFreq",				"Coupon Frequency (A/S/Q/M)",									&DefFreq			},
	{	GA_STRING,			"CpnDayCount",			"Coupon DayCount (A360;A365;30/360)",							&DefDayCount		},
	{	GA_DOUBLE,			"CpnResetGap",			"Coupon Reset Gap",												&DefResetGap		},
	{	GA_STRING,			"PayCal",				"Payment Calendar",												NULL				},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",												NULL				},
	{	GA_STRING,			"StubType",				"Stub Type (SS/LS)",											&DefStub			},
	{	GA_STRING,			"ResetTiming",			"Fixing Timing (ARR/ADV)",										&DefPayTiming		},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",									&DefIntRule			},
	{	GA_OBJECT,			"CpnNominal",			"Coupon Nominal",												NULL				},
	{	GA_OBJECT,			"DomCpn",				"Domestic Coupon",												NULL				},
	{	GA_OBJECT,			"FgnCpn",				"Foreign Coupon",												NULL				},
	{	GA_OBJECT,			"InitFX",				"Initial FX",													NULL				},
	{	GA_OBJECT,			"MinCpn",				"Minimum Coupon",												NULL				},
	{	GA_OBJECT,			"MaxCpn",				"Maximum Coupon",												NULL				},
	{	GA_OBJECT,			"Barrier",				"Barrier of the product in case of PRDKO",						&DefNullObject		},
	{	GA_OBJECT,			"NoticeType",			"Falgs to specify Either CALL or KO",							&DefNullObject		},
	{	GA_STRING,			"FundFreq",				"Funding Freq (A/S/Q/M)",										&DefFreq			},
	{	GA_STRING,			"FundDayCount",			"Funding DayCount (A360/A365/30/360)",							&DefDayCount		},
	{	GA_STRING,			"CompFreq",				"Compounding frequency (A,S,Q,M,W,D,ZC,None)",					&DefNone			},
	{	GA_STRING,			"CompType",				"Compounding type (SpreadInc, SpreadExc, Flat, None)",			&DefNone			},
	{	GA_OBJECT,			"FundNominal",			"Funding Nominal",												NULL				},
	{	GA_OBJECT,			"FundSpread",			"Funding Spread",												NULL				},
	{	GA_STRING,			"NoticeFreq",			"Notification  Frequency (A/S/Q/M)",							&DefFreq			},
	{	GA_DOUBLE,			"NoticeGap",			"Notification  Gap",											&DefResetGap		},
	{	GA_STRING,			"PayRec",				"Pay or Receive the coupon leg",								&DefPayRec			},
	{	GA_DOUBLE,			"NbNoCall",				"Nb no Call period",											&DefZero			},
	{	GA_STRING,			"RedemType",			"Redemption Type (STANDARD/MANDATORY/DUALOPTION)",				&DefRedeemType		},
	{	GA_DOUBLE,			"RedemGap",				"Redemption Gap",												&DefZero			},
	{	GA_DOUBLE,			"RedemStrike",			"Redemption Strike",											NULL				},
	{	GA_OBJECT,			"Fees",					"Call Fees",													NULL				},
	{	GA_OBJECT,			"ProductsToPrice",		"Products to price ",											NULL				},
	{	GA_STRING,			"FXLocalFlag",			"Flag to active/disactive the local vol adjustment(Y/N)",		&DefYes				},
	{	GA_STRING,			"BasisIRFlag",			"Flag to take in account the basis effect in IR calib(Y/N)",	&DefYes				},



	{	GA_TERMINAL }
};


const GenericParamStruct PRDCCalculator_InitParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"PRDCId",				"PRDC Calculator Object",										NULL			},
	{	GA_OBJECT,			"MktDataMger",			"market data manager",											NULL			},
	{	GA_OBJECT,			"SchedulerTree",		"Schedule of 3 factors tree",									NULL			},
	{	GA_OBJECT,			"TruncatorTree",		"Truncator of 3 factors tree",						            NULL			},
	{	GA_STRING,			"MarkovianDrift",		"Markovian drift sampler",										&DefYes			},
	{	GA_STRING,			"CalibType",			"Type of calibration",											&DefUnknown		},
	{	GA_OBJECT,			"CalibDatas",			"Datas of calibration",											NULL			},

	{	GA_TERMINAL }
};


/**************************************************************************************************/
const GenericParamStruct TARNFXCalculator_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",											NULL				},
	{	GA_STRING,			"CpnCcy",				"Coupon currency",												NULL				},
	{	GA_STRING_VECTOR,	"FgnCcy",				"Foreign currency (vector in case of Tarn Chooser)",			NULL				},
	{	GA_STRING,			"FundCcy",				"Funding currency",												NULL				},
	{	GA_STRING,			"PayRec",				"Pay or Receive the coupon leg",								&DefPayRec			},
	{	GA_STRING,			"CpnFreq",				"Coupon Frequency (A/S/Q/M)",									NULL				},
	{	GA_STRING,			"CpnDayCount",			"Coupon DayCount (A360/A365/30/360)",							NULL				},
	{	GA_DOUBLE,			"CpnResetGap",			"Coupon Reset Gap",												NULL				},
	{	GA_STRING,			"PayCal",				"Payment Calendar",												NULL				},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",												NULL				},
	{	GA_STRING,			"StubType",				"Stub Type (SS/LS)",											&DefStub			},
	{	GA_STRING,			"Timing",				"Fixing Timing (ARR/ADV)",										NULL				},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",									&DefIntRule			},
	{	GA_OBJECT,			"CpnNominal",			"Coupon Nominal",												NULL				},
	{	GA_OBJECT,			"DomCpn",				"Domestic Coupon",												NULL				},
	{	GA_OBJECT,			"FgnCpn",				"Foreign Coupon",												NULL				},
	{	GA_OBJECT,			"InitFX",				"Initial FX",													NULL				},
	{	GA_OBJECT,			"MinCpn",				"Minimum Coupon",												NULL				},
	{	GA_OBJECT,			"MaxCpn",				"Maximum Coupon",												NULL				},
	{	GA_STRING,			"FundFreq",				"Funding Freq (A/S/Q/M)",										NULL				},
	{	GA_STRING,			"FundDayCount",			"Funding DayCount (A360/A365/30/360)",							NULL				},
	{	GA_OBJECT,			"FundNominal",			"Funding Nominal",												NULL				},
	{	GA_OBJECT,			"FundSpread",			"Funding Spread",												NULL				},
	{	GA_STRING,			"CompFreq",				"Compounding frequency (A,S,Q,M,W,D,ZC)",						&DefUnknown			},
	{	GA_STRING,			"CompType",				"Compounding type (SpreadInc, SpreadExc, Flat)",				&DefUnknown			},
	{	GA_OBJECT,			"Target",				"Target of the product",										NULL				},
	{	GA_STRING,			"FXChoice",				"Condition on the currencies (For TARN Chooser)",				&DefTARNFXChoise	},
	{	GA_STRING,			"RedeemType",			"Redemption Type (STANDARD/MANDATORY/DUALOPTION)",				&DefRedeemType		},
	{	GA_DOUBLE,			"RedeemGap",			"Redemption Gap",												NULL				},
	{	GA_DOUBLE_VECTOR,	"RedeemStrike",			"Redemption Strike (vector in case of Tarn Chooser)",			NULL				},
	{	GA_OBJECT,			"Fees",					"Trigger Fees",													NULL				},
	{	GA_STRING,			"PayoffType",			"Payoff Type (TARNFX/PRDKO)",									&DefPayoffType		},
	{	GA_STRING,			"IntermPrices",			"Intermediate Prices Flag (Y/N)",								&DefYes				},
	{	GA_STRING_VECTOR,	"ProductsToPrice",		"Vector with names of columns to price",						&DefObjParam		},
	{	GA_OBJECT,			"FixingSched",			"Fixing Sched",													&DefObjParam		},

	{	GA_TERMINAL }
};

const GenericParamStruct TARNFXCalculator_InitParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"TARNFXId",				"TARN FX Calculator Object",									NULL				},
	{	GA_DOUBLE,			"NbSimul",				"Number of Monte Carlo Simulation",								&DefNbSimul			},
	{	GA_DOUBLE,			"BucketSize",			"Size of Monte Carlo Simulation Bucket",						&DefNbSimul			},
	{	GA_STRING,			"RandGenType1",			"Random Generator 1 (MRGK5/NR_RAN2 ...)",						&DefRandGenType1	},
	{	GA_STRING,			"RandGenAlgo1",			"Random Generator Algorithm 1 (BoxMuller/InvNormCum ...)",		&DefRandGenAlgo1	},
	{	GA_STRING,			"RandGenType2",			"Random Generator 2 (MRGK5/NR_RAN2 ...)",						&DefRandGenType2	},
	{	GA_STRING,			"RandGenAlgo2",			"Random Generator Algorithm 2 (BoxMuller/InvNormCum ...)",		&DefRandGenAlgo2	},
	{	GA_DOUBLE,			"FirstNbDims",			"First Nb Dimensions used in the Mixed Generator",				&DefZero			},
	{	GA_DOUBLE,			"FirstNbTimes",			"First Nb Times used in the Mixed Generator",					&DefZero			},
	{	GA_DOUBLE,			"FactorNb",				"Model Number of factors",										&DefNbFactors		},
	{	GA_DOUBLE,			"TimeStepNb",			"PDE Number of Time Steps",										&DefTimeSteps		},
	{	GA_DOUBLE,			"SpaceStepNb",			"PDE Number of Space Steps",									&DefSpaceSteps		},
	{	GA_DOUBLE,			"StdDevNb",				"PDE Number of Standard Deviation",								&DefNbStdDevs		},
	{	GA_STRING,			"SkipPDE",				"Skip PDE flag (Y/N)",											&DefNo				},
	{	GA_STRING,			"Rescalling",			"Rescalling flag (Y/N)",										&DefNo				},
	{	GA_STRING,			"ModelType",			"Choice of the model (2IRFX/1IRFX)",							&DefTARNFXModel		},
	{	GA_STRING,			"Smile",				"Smile Flag (Y/N)",												&DefYes				},
	{	GA_STRING,			"MixCalib",				"Mixture Model Calibraton Flag (Y/N)",							&DefNo				},
	{	GA_STRING,			"OneFactor",			"One Factor Flag (Deterministic IR) (Y/N)",						&DefNo				},
	{	GA_STRING,			"CorrelType",			"Type of correlation (Beta/CorrelMatrix/Fwd)",					&DefCorrType		},
	{	GA_OBJECT,			"MDM",					"Market Data Manager",											NULL				},

	{	GA_TERMINAL }
};


/**************************************************************************************************/
const GenericParamStruct TARNCalculatorIndian_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",											NULL				},
	{	GA_STRING,			"CpnCcy",				"Coupon currency",												NULL				},
	{	GA_STRING,			"FgnCcy",				"Foreign currency",												NULL				},
	{	GA_STRING,			"PayRec",				"Pay or Receive the coupon leg",								&DefPayRec			},
	{	GA_STRING,			"CpnFreq",				"Coupon Frequency (A/S/Q/M)",									NULL				},
	{	GA_STRING,			"CpnDayCount",			"Coupon DayCount (A360/A365/30/360)",							NULL				},
	{	GA_DOUBLE,			"CpnResetGap",			"Coupon Reset Gap",												NULL				},
	{	GA_STRING,			"PayCal",				"Payment Calendar",												NULL				},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",												NULL				},
	{	GA_STRING,			"StubType",				"Stub Type (SS/LS)",											&DefStub			},
	{	GA_STRING,			"Timing",				"Fixing Timing (ARR/ADV)",										NULL				},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",									&DefIntRule			},
	{	GA_OBJECT,			"CpnNominal",			"Coupon Nominal",												NULL				},
	{	GA_OBJECT,			"Strike",				"Strike",														NULL				},
	{	GA_OBJECT,			"BarrierUp",			"Barrier Up",													NULL				},
	{	GA_OBJECT,			"BarrierDown",			"Barrier Down",													&DefObjParam		},
	{	GA_OBJECT,			"Target",				"Target of the product",										NULL				},
	{	GA_DOUBLE,			"Epsilon",				"Coefficient on strike",										NULL				},
	{	GA_OBJECT,			"Fees",					"Trigger Fees",													NULL				},
	{	GA_STRING,			"OptionType",			"Call or Put",													&DefCallPut			},
	{	GA_STRING,			"IndianType",			"DownUp, Trigger, DigitalTrigger",								NULL				},
	{	GA_STRING,			"IntermPrices",			"Intermediate Prices Flag (Y/N)",								&DefYes				},
	{	GA_STRING_VECTOR,	"ProductsToPrice",		"Vector with names of columns to price",						&DefObjParam		},

	{	GA_TERMINAL }
};


const GenericParamStruct ParamsMixtureFx_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"Lags",					"Lags in number of days (vector of double)",					NULL				},
	{	GA_OBJECT,			"VolATM",				"ATM volatility (vector of double)",							NULL				},
	{	GA_OBJECT,			"DecVol",				"Gap between vol ATM and the vol1 (vector of double)",			NULL				},
	{	GA_OBJECT,			"Shift",				"Shift of the fwd for the regime 1 (vector of double)",			NULL				},
	{	GA_OBJECT,			"Lambda",				"Probability of the regime 1 (vector of double)",				NULL				},
	{	GA_STRING,			"Interpol",				"Interpolation mode (STEPUPRIGHT/STEPUPRIGHT/LINEAR)",			NULL				},

	{	GA_TERMINAL }
};


const GenericParamStruct MixtureModelFx_CreateWithParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"DomCrv",				"Domestic Yield Curve",											NULL				},
	{	GA_OBJECT,			"ForCrv",				"Foreign Yield Curve",											NULL				},
	{	GA_DOUBLE,			"Spot",					"FX spot rate",													NULL				},
	{	GA_OBJECT,			"MixParams",			"Mixture parameters",											NULL				},

	{	GA_TERMINAL }
};

const GenericParamStruct FXModelDensity_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"Pricing model object",											NULL				},
	{	GA_STRING,			"DensityType",			"Normale Density or InvNormale Density",							&DefUnknown			},
	{	GA_DOUBLE,			"Expiry",				"Expiry of the forward in DAYS",								&DefZero			},
	{	GA_DOUBLE,			"Xmin",					"xmin",															&DefZero			},
	{	GA_DOUBLE,			"Xmax",					"xmax",															&DefZero			},
	{	GA_DOUBLE,			"NbPoints",				"Nb de Points between xmin and xmax",							&DefZero			},

	{	GA_TERMINAL }
};

const GenericParamStruct ModelBumpParams_CreateWithParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"Pricing model object",											NULL				},
	{	GA_STRING,			"ParamType",			"Model arameter type (volatility,Correlation,...)",				&DefUnknown			},
	{	GA_DOUBLE,			"RowNumber",			"Number of rows",												&DefZero			},
	{	GA_DOUBLE,			"ColumnNumber",			"Number of Columns ",											&DefZero			},
	{	GA_DOUBLE,			"Shift",				"Shift ( 0.01 for 1%)",											&DefZero			},
	{	GA_STRING,			"IsCumulative",			"Type of bump ( IsCumulative/ IsPerturbative)",					&DefIsNothing		},

	{	GA_TERMINAL }
};

const GenericParamStruct XXXMktDataFromMktDataMgerWithParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"MktDataMgerId",		"Mkt Data Manger object",											NULL				},

	{	GA_TERMINAL }

};

/**************************************************************************************************/
const GenericParamStruct CCSCalculator_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"Asof Date of the deal",										NULL				},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",											NULL				},
	{	GA_STRING,			"DomCcy",				"Domestic currency",											NULL				},
	{	GA_STRING,			"DomFreq",				"Domestic Frequency (A/S/Q/M)",									&DefFreq			},
	{	GA_STRING,			"DomDayCount",			"Domestic DayCount (A360/A365/30/360)",							&DefDayCount		},
	{	GA_STRING,			"PayCal",				"Payment Calendar",												&DefCurrency		},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",												&DefCurrency		},
	{	GA_OBJECT,			"DomNominal",			"Domestic Nominal",												NULL				},
	{	GA_OBJECT,			"DomSpread",			"Domestic Spread",												NULL				},
	{	GA_STRING,			"ForCcy",				"Foreign currency (A360/A365/30/360)",							NULL				},
	{	GA_STRING,			"ForFreq",				"Foreign Frequency (A/S/Q/M)",									&DefFreq			},
	{	GA_STRING,			"ForDayCount",			"Foreign DayCount (A360/A365/30/360)",							&DefDayCount		},
	{	GA_STRING,			"ForPayCal",			"Foreign Payment Calendar",										&DefCurrency		},
	{	GA_STRING,			"ForResetCal",			"Foreign Reset Calendar",										&DefCurrency		},
	{	GA_OBJECT,			"ForNominal",			"Foreign Nominal",												NULL				},
	{	GA_OBJECT,			"ForSpread",			"Foreign Spread",												NULL				},
	{	GA_DOUBLE,			"FxResetGap",			"Forex Reset Gap",												&DefResetGap		},
	{	GA_STRING,			"PayRec",				"Pay or Receive the Domestic leg",								&DefPayRec			},
	{	GA_DOUBLE,			"NbNoCall",				"Number of no call periods",									&DefZero			},
	{	GA_STRING,			"NoticeFreq",			"Notification Frequency (A/S/Q/M)",								&DefFreq			},
	{	GA_DOUBLE,			"NoticeGap",			"Notification Gap",												&DefResetGap		},	
	{	GA_STRING,			"StubType",				"Stub Type (SS/LS)",											&DefStub			},
	{	GA_STRING,			"ResetTiming",			"Fixing Timing (ARR/ADV)",										NULL				},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",									&DefIntRule			},
	{	GA_OBJECT,			"Fees",					"Fees",															NULL				},
	{	GA_OBJECT,			"ProductsToPrice",		"payoffs to price",												NULL				},

	{	GA_TERMINAL }
};

const GenericParamStruct CCSCalculator_InitParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"CCSCalculator",		"CCS Calculator Object",										NULL			},
	{	GA_OBJECT,			"MktDataMger",			"market data manager",											NULL			},
	{	GA_OBJECT,			"SchedulerTree",		"Schedule of 3 factors tree",									NULL			},
	{	GA_OBJECT,			"TruncatorTree",		"Truncator of 3 factors tree",						            NULL			},
	{	GA_STRING,			"MarkovianDrift",		"Markovian drift sampler",										&DefYes			},
	{	GA_STRING,			"CalibType",			"Type of calibration",											&DefUnknown		},
	{	GA_OBJECT,			"CalibDatas",			"Datas of calibration",											NULL			},

	{	GA_TERMINAL }
};


const GenericParamStruct QDensityFunctor_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_DOUBLE,			"vol",					"double",														NULL			},
	{	GA_DOUBLE,			"q",					"double",														NULL			},

	{	GA_TERMINAL }
};

const GenericParamStruct SABR_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"alpha",				"double",														&DefParam				},
	{	GA_DOUBLE,			"beta",					"double",														&DefParam			},
	{	GA_DOUBLE,			"rho",					"double",														&DefParam			},
	{	GA_DOUBLE,			"nu",					"double",														&DefParam			},
	{	GA_STRING,			"sabrflag",				"string",														&DefSabrFlag		},
	{	GA_DOUBLE,			"shift",				"double (only if beta = 1)",									&DefOne				},
	{	GA_STRING,			"localRhoCalib",		"string",														&DefNo				},
	{	GA_TERMINAL }
};

const GenericParamStruct SABR2B_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"alpha",				"double",														&DefParam			},
	{	GA_DOUBLE,			"beta1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"beta2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rho",					"double",														&DefParam			},
	{	GA_DOUBLE,			"nu",					"double",														&DefParam			},
	{	GA_DOUBLE,			"zero",					"double",														&DefParam			},
	{	GA_DOUBLE,			"lambda",				"double",														&DefParam			},

	{	GA_TERMINAL }
};


const GenericParamStruct Heston_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"v0",					"double",														NULL				},
	{	GA_DOUBLE,			"kappa",				"double",														&DefParam			},
	{	GA_DOUBLE,			"theta",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rho",					"double",														&DefParam			},
	{	GA_DOUBLE,			"vvol",					"double",														&DefParam			},
	{	GA_DOUBLE,			"shift",				"double",														&DefOne				},
	{	GA_DOUBLE,			"level",				"double",														&DefOne				},
	{	GA_STRING,			"localRhoCalib",		"string",														&DefNo				},
	{	GA_STRING,			"bootstraplevel",		"string",														&DefNo				},
	{	GA_STRING,			"normalheston",			"string",														&DefNo				},
	{	GA_OBJECT,			"sigma",				"double",														&DefObjParam		},
	{	GA_DOUBLE,			"weight",				"double",														&DefZero			},
	{	GA_TERMINAL }
};

const GenericParamStruct Heston2b_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"v01",					"double",														NULL				},
	{	GA_DOUBLE,			"kappa1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"theta1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rho1",					"double",														&DefParam			},
	{	GA_DOUBLE,			"vvol1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"v02",					"double",														NULL				},
	{	GA_DOUBLE,			"kappa2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"theta2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rho2",					"double",														&DefParam			},
	{	GA_DOUBLE,			"vvol2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"shift",				"double",														&DefOne				},
	{	GA_DOUBLE,			"level",				"double",														&DefOne				},
	{	GA_STRING,			"localRhoCalib1",		"string",														&DefNo				},
	{	GA_STRING,			"localRhoCalib2",		"string",														&DefNo				},
	{	GA_STRING,			"bootstraplevel",		"string",														&DefNo				},
	{	GA_TERMINAL }
};

const GenericParamStruct BiSABR_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"forward1",				"double",														NULL				},
	{	GA_DOUBLE,			"alpha1",				"double",														NULL				},
	{	GA_DOUBLE,			"beta1",				"double",														NULL				},
	{	GA_DOUBLE,			"rho1",					"double",														NULL				},
	{	GA_DOUBLE,			"nu1",					"double",														NULL				},
	{	GA_DOUBLE,			"forward2",				"double",														NULL				},
	{	GA_DOUBLE,			"alpha2",				"double",														NULL				},
	{	GA_DOUBLE,			"beta2",				"double",														NULL				},
	{	GA_DOUBLE,			"rho2",					"double",														NULL				},
	{	GA_DOUBLE,			"nu2",					"double",														NULL				},
	{	GA_DOUBLE,			"rhoS1S2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rhoS1V2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rhoS2V1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"rhoV1V2",				"double",														&DefParam			},

	{	GA_TERMINAL }
};

const GenericParamStruct Merton_SmileCalibration_Params[] =
{
	{	GA_DOUBLE,			"sigma",				"double",														&DefParam			},
	{	GA_DOUBLE,			"lambda1",				"double",														&DefParam			},
	{	GA_DOUBLE,			"U1",					"double",														&DefParam			},
	{	GA_DOUBLE,			"lambda2",				"double",														&DefParam			},
	{	GA_DOUBLE,			"U2",					"double",														&DefParam			},

	{	GA_TERMINAL }
};

/**************************************************************************************************/

const GenericParamStruct FXVanillaCalculator_CreateParams[] = 
{
	///////////////			////////////////////	//////////////////////////////////////////////////////////////////////////////			///////////////
	//	Param Type			Param Name				ParamDescription																		DefaultValue

	{	GA_DOUBLE,			"AsOfDate",				"Asof Date of the deal",														NULL				},
	{	GA_OBJECT,			"DateStrip",			"DateStrip of the deal ",														&DefNullObject		},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",														NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",															NULL				},
	{	GA_DOUBLE,			"ExpiryGap",			"Forex Expiry Gap toward the EndDate",											&DefResetGap		},
	{	GA_DOUBLE,			"SetlmentGap",		    "Forex Settlement Gap toward the Expiry",										&DefResetGap		},
	{	GA_STRING,			"Frequency",			"Exercise Frequency (A/S/Q/M)",													&DefFreq			},
	{	GA_STRING,			"DayCount",			    "DayCount (A360/A365/30/360)",													&DefDayCount		},
	{	GA_STRING,			"PayCal",				"Payment Calendar",																NULL				},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",																NULL				},
	{	GA_STRING,			"FX1name",				"Name of the first FX (FOR/DOM)",												NULL				},
	{	GA_STRING,			"FX2name",				"Name of the second FX (FOR/DOM)",												NULL				},
	{	GA_STRING,			"PayCcy",			    "Payment currency",																&DefDayCount		},
	{	GA_OBJECT,			"Nominal",				"Nominal",																		NULL				},
	{	GA_STRING,			"CallPut",				"Call or Put for the option on FX1  ",											&DefCallPut			},
	{	GA_OBJECT,			"Strike",				"Common strike",																NULL				},
	{	GA_OBJECT,			"Alpha",				"Alpha (Spread)",																NULL				},
	{	GA_OBJECT,			"Beta",					"Beta (Spread)",																&DefNullObject		}, // &DefObjParam	},
	{	GA_OBJECT,			"Strike2",				"Strike of the second performance (Basket)",									&DefNullObject		},
	{	GA_OBJECT,			"Leverage",				"Leverage of the option",														NULL                },
	{	GA_STRING,			"CallPut2",				"Call or Put for the option on FX2 (Basket)  ",									&DefCallPut			},
	{	GA_STRING,			"MinMax",				"BestOf or WorstOf (Basket)",													&DefMinMax			},
	{	GA_STRING,			"DigitType",			"ANALYTIC,CENTRED,BACKWARD,FORWARD",											&DefDigitType		},
	{	GA_DOUBLE,			"Epsilon",				"Epsilon if the call spread for the digit",										&DefEpsilon			},
	{	GA_STRING,			"OptionType",			"Vanilla,Spread,Basket,Digit,Perf,DigitSpread,Quotient,FXBall,FXBallPerf",      &DefOptionType		},
	{	GA_OBJECT,			"CouponMin",			"Lower bound for MinMax coupon",												&DefNullObject		},
	{	GA_OBJECT,			"CouponMax",			"Upper bound for MinMax coupon",												&DefNullObject		},

	{	GA_TERMINAL }
};
const GenericParamStruct FXVanillaCalculator_InitParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"FXCalculatorId",		"FX Vanilla Calculator Object",									NULL				},
	{	GA_OBJECT,			"Kmin1",				"Matrix of Kmin for the first FX",								NULL				},
	{	GA_OBJECT,			"Kmax1",				"Matrix of Kmax for the first FX",								NULL				},
	{	GA_OBJECT,			"Kmin2",				"Matrix of Kmin for the second FX",								NULL				},
	{	GA_OBJECT,			"Kmax2",				"Matrix of Kmin for the second FX",								NULL				},
	{	GA_DOUBLE,			"Nleft1",				"Number of points for the left part of the FX1 smile",			NULL				},
	{	GA_DOUBLE,			"Ncenter1",				"Number of points for the center part of the FX1 smile",		NULL				},
	{	GA_DOUBLE,			"Nright1",				"Number of points for the right part of the FX1 smile",			NULL				},
	{	GA_DOUBLE,			"Nleft2",				"Number of points for the left part of the FX2 smile",			NULL				},
	{	GA_DOUBLE,			"Ncenter2",				"Number of points for the center part of the FX2 smile",		NULL				},
	{	GA_DOUBLE,			"Nright2",				"Number of points for the right part of the FX2 smile",			NULL				},
	{	GA_OBJECT,			"MktDataManager",		"Market Data Manager",											NULL				},
	{	GA_TERMINAL }
};

/**************************************************************************************************/

const GenericParamStruct FXRACalculator_CreateParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue


	{	GA_DOUBLE,			"AsOfDate",				"Asof Date of the deal",										NULL				},
	{	GA_OBJECT,			"DateStrip",			"DateStrip of the deal ",									    NULL				},
	{	GA_DOUBLE,			"StartDate",			"Start Date of the deal",										NULL				},
	{	GA_DOUBLE,			"EndDate",				"End Date of the deal",											NULL				},
	{	GA_DOUBLE,			"ExpiryGap",			"Forex Expiry Gap toward the EndDate",							&DefResetGap		},
	{	GA_DOUBLE,			"SetlmentGap",		    "Forex Settlement Gap toward the Expiry",						&DefResetGap		},
	{	GA_DOUBLE,			"PaymentGap",			"Forex Payment Gap toward the Settlement",						&DefResetGap		},
	{	GA_STRING,			"Frequency",			"Exercise Frequency (A/S/Q/M)",									&DefFreq			},
	{	GA_STRING,			"DayCount",			    "DayCount (A360/A365/30/360)",							        &DefDayCount		},
	{	GA_STRING,			"PayCal",				"Payment Calendar",												NULL				},
	{	GA_STRING,			"ResetCal",				"Reset Calendar",												NULL				},
	{	GA_STRING,			"FX1name",				"Name of the first FX (FOR/DOM)",								NULL				},
	{	GA_STRING,			"FX2name",				"Name of the second FX (FOR/DOM)",								NULL				},
	{	GA_STRING,			"PayCcy",			    "Payment currency",												&DefDayCount		},
	{	GA_OBJECT,			"Nominal",				"Nominal",												        NULL				},
	{	GA_STRING,			"CallPut",				"Call or Put for the option on FX1  ",							&DefCallPut			},
	{	GA_OBJECT,			"Strike",				"Common strike",										        NULL				},
	{	GA_OBJECT,			"Alpha",				"Alpha (Spread)",												NULL				},
	{	GA_OBJECT,			"Beta",					"Beta (Spread)",										        NULL				},
	{	GA_OBJECT,			"Strike2",				"Strike of the second performance (Basket)",			        NULL				},
	{	GA_OBJECT,			"Leverage",				"Leverage of the option",		                                NULL                },
	{	GA_STRING,			"CallPut2",				"Call or Put for the option on FX2 (Basket)  ",					&DefCallPut			},
	{	GA_STRING,			"MinMax",				"BestOf or WorstOf (Basket)",									&DefMinMax			},
	{	GA_STRING,			"DigitType",			"ANALYTIC,CENTRED,BACKWARD,FORWARD",							&DefDigitType		},
	{	GA_DOUBLE,			"Epsilon",				"Epsilon if the call spread for the digit",						&DefEpsilon			},
	{	GA_STRING,			"OptionType",			"Vanilla,Spread,Basket,Digit,Perf,DigitSpread,Quotient",        NULL				},
	{	GA_STRING,			"IntRule",				"Interest Rule (ADJ/UNADJ)",									&DefIntRule			},
	{	GA_STRING,			"StubType",				"Stub Type (SS/LS)",											&DefStub			},
	{	GA_STRING,			"ResetTiming",			"Fixing Timing (ARR/ADV)",										NULL				},
	{	GA_STRING,			"FixingFrequency",		"Fixing Frequency (A/S/Q/M/W/D)",								&DefFreq			},
	{	GA_STRING,			"PAYidx",				"Index of the paid rate (FIXED/LIBOR/CMS)",						NULL				},
	{	GA_DOUBLE,			"PAYidxSpread",			"Spread on the Payment Index",									NULL				},
	{	GA_STRING,			"PAYidxIT",				"Paid Index Term = nb followed by d/w/m/y like 6m or 1y",		NULL				},
	{	GA_STRING,			"IRidx",				"Type of the Ir fixing index (LIBOR/CMS)",						NULL				},
	{	GA_STRING,			"IRidxIT",				"Paid Index Term  = nb followed by d/w/m/y like 6m or 1y",		NULL				},
	{	GA_OBJECT,			"FxDownBarrier",		"Down Barrier for FX",										    NULL				},
	{	GA_OBJECT,			"FxUpBarrier",			"Up Barrier for FX",											NULL				},
	{	GA_OBJECT,			"IrDownBarrier",		"Down Barrier for IR",											NULL				},
	{	GA_OBJECT,			"IrUpBarrier",			"Up Barrier for IR",											NULL				},
	{	GA_TERMINAL }
};

const GenericParamStruct FXRACalculator_InitParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"FXCalculatorId",		"FX Vanilla Calculator Object",									NULL				},
	{	GA_DOUBLE,			"NbPoints1",			"Number of points for GL integration of the first ind",			NULL				},
	{	GA_DOUBLE,			"NbPoints2",			"Number of points for GL integration of the second ind",		NULL				},
	{	GA_OBJECT,			"MktDataManager",		"Market Data Manager",											NULL				},
	{	GA_TERMINAL }
};
const GenericParamStruct XXX_MktData_ApplyScenarioParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"MktDataId",			"MktData Object",												NULL				},
	{	GA_OBJECT,			"ScenarioId",			"Scenario Object",												NULL				},
	{	GA_DOUBLE,			"NbShift",				"NbShift",														&DefOne				},
	{	GA_TERMINAL }
};

const GenericParamStruct VanillaDensityFunctorParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"pricing model Id",												NULL				},
	{	GA_DOUBLE,			"ExpiryDate",			"Expiry date",													NULL				},
	{	GA_STRING,			"IsDirect",				"bool to choose th sens of density",							&DefYes 			},
	{	GA_TERMINAL }
};


const GenericParamStruct ARM_2IRFX_ComputeTimeLagFunctorParams[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"2IRFX model Id",												NULL				},
	{	GA_DOUBLE,			"ExpiryTime",			"Expiry time",													NULL				},
	{	GA_DOUBLE,			"PayTime",				"Pay time",														NULL 				},
	{	GA_STRING,			"DomForFlag",			"Dom or Foreign Flag",											&DefDom 				},
	{	GA_TERMINAL }
};

const GenericParamStruct ARM_2IRFX_ComputeFwdFxVolFunctorParams[] =
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"2IRFX model Id",												NULL				},
	{	GA_DOUBLE_VECTOR,	"ResetTimes",			"Reset times",													NULL				},
	{	GA_DOUBLE_VECTOR,	"SettlTimes",			"Settlment times",												NULL 				},
	{	GA_TERMINAL }
};

const GenericParamStruct ARM_2IRFX_ComputeFwdFxModelParamFunctor[] = 
{
	///////////////			////////////////////	////////////////////////////////////////////////////////////	///////////////
	//	Param Type			Param Name				ParamDescription												DefaultValue
	
	{	GA_OBJECT,			"ModelId",				"2IRFX model Id",												NULL				},
	{	GA_DOUBLE,			"EvalTime",				"Eval time",													NULL				},
	{	GA_DOUBLE,			"SettlementTime",		"Settlement time",														NULL 				},
	{	GA_TERMINAL }
};

const GenericAddinStruct GenericAddinsTable[] = 
{
//	Function Name								Function Description									Function Params							Functor
	{"ARM_PRDCCalculator_Create",				"Build a global PRDC calculator object.",				PRDCCalculator_CreateParams,			&ThePRDCCalculator_CreateFunctor,				true },
	{"ARM_PRDCCalculator_Init",					"Init a PRDC calculator object.",						PRDCCalculator_InitParams,				&ThePRDCCalculator_InitFunctor,					true },
	{"TARNFXCalculator_Create",					"Build a TARN FX calculator object.",					TARNFXCalculator_CreateParams,			&TheTARNFXCalculator_CreateFunctor,				true },
	{"TARNFXCalculator_Init",					"Init a TARN FX calculator object.",					TARNFXCalculator_InitParams,			&TheTARNFXCalculator_InitFunctor,				true },
	{"TARNCalculatorIndian_Create",				"Build an indian TARN calculator object.",				TARNCalculatorIndian_CreateParams,		&TheTARNCalculatorIndian_CreateFunctor,			true },
	{"ARM_GP_CCSCalculator_Create",				"Build a CCS calculator object.",						CCSCalculator_CreateParams,				&TheCCSCalculator_CreateFunctor,				true },
	{"ARM_GP_CCSCalculator_Init",				"Init a CCS  calculator object.",						CCSCalculator_InitParams,				&TheCCSCalculator_InitFunctor,					true },
	{"ARM_GP_FXVanillaCalculator_Create",	    "Build a FXVanilla calculator object.",					FXVanillaCalculator_CreateParams,		&TheFXVanillaCalculator_CreateFunctor,			true },
	{"ARM_GP_FXVanillaCalculator_Init",			"Init a FXVanilla  calculator object.",					FXVanillaCalculator_InitParams,			&TheFXVanillaCalculator_InitFunctor,			true },
	{"ARM_GP_FXRACalculator_Create",			"Build a FXRA calculator object.",						FXRACalculator_CreateParams,			&TheFXRACalculator_CreateFunctor,				true },
	{"ARM_GP_FXRACalculator_Init",				"Init a FXRA  calculator object.",						FXRACalculator_InitParams,				&TheFXRACalculator_InitFunctor,					true },	
	{"ParamsMixtureFx_Create",					"Create the  Mixture Fx  parameter.",					ParamsMixtureFx_CreateParams,			&TheParamsMixtureFx_CreateFunctor,				true },
	{"MixtureModelFx_CreateWithParams",			"Create the  Mixture Model Fx with parameters.",		MixtureModelFx_CreateWithParams,		&TheMixtureModelFx_CreateWithParamsFunctor,		true },
	{"ARM_GP_ModelBumpParams",					"Bump the parameters of a given  GP Model.",			ModelBumpParams_CreateWithParams,   	&TheModelBumpParams_CreateFunctor,				true },
	{"ARM_GP_XXXMktDataFromMktDataMger",		"Create a XXX Mkt Data from Mkt Data Manager.",			XXXMktDataFromMktDataMgerWithParams,	&TheXXXMktDataFromMktDataMgerFunctor,			true },
	{"ARM_GP_FXModelDensity_Create",			"Create the Normale or Invnormale forward density.",	FXModelDensity_CreateParams,  			&TheFXModelDensity_CreateFunctor,				true },

	{"Corridor_InfIr_Create",					"Build a Inflation Corridor Leg object.",				Corridor_InfIr_CreateParams,			&TheCorridor_InfIr_CreateFunctor,				true },
	{"Hybrid_Inf/IrLeg_Create",					"Build a Inf/Ir leg instrument.",						HybridInfIrLeg_CreateParams,			&TheHybridInfIrLeg_CreateFunctor,				true },
	{"Hybrid_Inf/IrPayOff_Create",				"Build a Inf/Ir PayOff.",								HybridInfIrPayOff_CreateParams,			&TheHybridInfIrPayOff_CreateFunctor,			true },
	{"Hybrid_Inf/IrModel_Create",				"Build a Inf/Ir Model.",								HybridInfIrModel_CreateParams,			&TheHybridInfIrModel_CreateFunctor,				true },

	{"ParamsMixtureFx_Create",					"Create the  Mixture Fx  parameter.",					ParamsMixtureFx_CreateParams,			&TheParamsMixtureFx_CreateFunctor,				true },
	{"MixtureModelFx_CreateWithParams",			"Create the  Mixture Model Fx with parameters.",		MixtureModelFx_CreateWithParams,		&TheMixtureModelFx_CreateWithParamsFunctor,		true },
	{"ARM_GP_ModelBumpParams",					"Bump the parameters of a given  GP Model.",			ModelBumpParams_CreateWithParams,   	&TheModelBumpParams_CreateFunctor,				true },
	{"ARM_GP_XXXMktDataFromMktDataMger",		"Create a XXX Mkt Data from Mkt Data Manager.",			XXXMktDataFromMktDataMgerWithParams,	&TheXXXMktDataFromMktDataMgerFunctor,			true },
	{"ARM_GP_FXModelDensity_Create",			"Create the Normale or Invnormale forward density.",	FXModelDensity_CreateParams,  			&TheFXModelDensity_CreateFunctor,				true },

	{"SABR_SmileCalib_Params_Create",			"Create a structure of sabr params",					SABR_SmileCalibration_Params,			&TheSABRSmileCalib_CreateFunctor,				true },
	{"SABR2B_SmileCalib_Params_Create",			"Create a structure of sabr 2 beta params",				SABR2B_SmileCalibration_Params,			&TheSABR2BSmileCalib_CreateFunctor,				true },
	{"Heston_SmileCalib_Params_Create",			"Create a structure of heston params",					Heston_SmileCalibration_Params,			&TheHestonSmileCalib_CreateFunctor,				true },
	{"Heston2b_SmileCalib_Params_Create",		"Create a structure of heston2b  params",				Heston2b_SmileCalibration_Params,		&TheHeston2bSmileCalib_CreateFunctor,			true },
	{"BiSABR_SmileCalib_Params_Create",			"Create a structure of bi-sabr params",					BiSABR_SmileCalibration_Params,			&TheBiSABRSmileCalib_CreateFunctor,				true },
	{"Merton_SmileCalib_Params_Create",			"Create a structrue of merton params",					Merton_SmileCalibration_Params,			&TheMertonSmileCalib_CreateFunctor,				true },
	{"QDensityFunctor_Create",					"Create a Q density functor.",							QDensityFunctor_CreateParams,			&TheQDensityFunctor_CreateFunctor,				true },
	{"XXX_MktData_ApplyScenario",				"Apply a scenario to the market data.",					XXX_MktData_ApplyScenarioParams,		&TheXXX_MktData_ApplyScenarioFunctor,			true },
	{"ARM_GP_VanillaDensityFunctor",			"Create a vanilla Security Density.",					VanillaDensityFunctorParams,			&TheVanillaDensityFunctor,						true },
	{"2IRFX_ComputeTimeLag",					"Compute a time lag for a forward FX.",					ARM_2IRFX_ComputeTimeLagFunctorParams,	&The2IRFX_ComputeTimeLagFunctor,				false },
	{"2IRFX_ComputeFwdFxModelParam",			"Compute forward FX Q parameter and volatility.",		ARM_2IRFX_ComputeFwdFxModelParamFunctor,&The2IRFX_ComputeFwdFxModelParamFunctor,		false },
	{"2IRFX_ComputeFwdFxVolFunctor",			"Compute the theorical forward FX vol.",				ARM_2IRFX_ComputeFwdFxVolFunctorParams,	&The2IRFX_ComputeFwdFxVolFunctor,				false },
	

};

void ARM_GenericParamValue::Init()
{
	itsType = GA_DOUBLE;
	itsDoubleValue = 0.0;
	itsStringValue = "";
	itsObjectValue = 0;
}

ARM_GenericParamValue::ARM_GenericParamValue(double value)
{
	Init();
	itsType = GA_DOUBLE;
	itsDoubleValue = value;
}

ARM_GenericParamValue::ARM_GenericParamValue(const string& value)
{
	Init();
	itsType = GA_STRING;
	itsStringValue = value;
}

ARM_GenericParamValue::ARM_GenericParamValue(long value)
{
	Init();
	itsType = GA_OBJECT;
	itsObjectValue = value;
};

double ARM_GenericParamValue::GetDouble() const
{
	if (itsType != GA_DOUBLE)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "This parameter value is not a double.");;

	return itsDoubleValue;
}

const string &ARM_GenericParamValue::GetString() const
{
	if (itsType != GA_STRING)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "This parameter value is not a string.");;

	return itsStringValue;
}

long ARM_GenericParamValue::GetObjectId() const
{
	if ((itsType != GA_OBJECT) && (itsType != GA_DOUBLE_VECTOR) && (itsType != GA_STRING_VECTOR))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "This parameter value is not an object.");;

	return itsObjectValue;
}

ARM_GenericParamsDesc::ARM_GenericParamsDesc(const string& functionName, const GenericParamStruct* paramsStruct)
: itsFunctionName(functionName)
{
	while (paramsStruct->Type != GA_TERMINAL)
	{
		string nameUC(paramsStruct->ParamName);
		stringToUpper(nameUC);
		stringTrim(nameUC);

		ARM_GenericParamDesc param(*paramsStruct);

		GenericParamMap::value_type p(nameUC,param);

		itsParams.insert(p);
		itsParamNames.push_back(paramsStruct->ParamName);
		++paramsStruct;
	}
}

const ARM_GenericParamDesc& ARM_GenericParamsDesc::GetParam(const string& name) const
{
	string nameUC(name);
	stringToUpper(nameUC);

	GenericParamMap::const_iterator it = itsParams.find(nameUC);
	
	if (it == itsParams.end())
	{
		CC_Ostringstream os;

		os << "The parameter " << name << " is missing in the description in " << itsFunctionName << ".";

		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());;
	}

	return it->second;
}

void ARM_GenericParamsDesc::fillDefaultParams(ARM_GenericParams* genericParams)
{
	for (GenericParamMap::iterator it = itsParams.begin(); it != itsParams.end(); ++it)
	{
		string paramName(it->second.GetParamName());
		stringToUpper(paramName);
		stringTrim(paramName);

		if (!genericParams->IsExists(paramName))
		{
			if (it->second.GetDefaultParamValue() != NULL)
				genericParams->SetParamValue(paramName, (*it->second.GetDefaultParamValue()));
			else
			{
				CC_Ostringstream os;
				os << "Default value is missing for the parameter " << paramName << " in " << itsFunctionName << ".";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());;
			}
		}
	}
}

vector<string> ARM_GenericParamsDesc::GetParamNames(bool withDefault) const
{
	vector<string> paramNames;

	vector<string>::const_iterator it;
	for (it = itsParamNames.begin(); it != itsParamNames.end(); ++it)
	{
		if (!GetParam(*it).GetDefaultParamValue() || withDefault)
			paramNames.push_back(GetParam(*it).GetParamName());
	}

	return paramNames;
}

string ARM_GenericParamsDesc::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << std::setiosflags(std::ios::left);
	os << indent;
	os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15) << "Type";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(20) << "Name";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(60) << "Description";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15) << "Default" << endl;

	os << indent << "**************************************************************************************************************" << endl;

	vector<string>::const_iterator it;
	for (it = itsParamNames.begin(); it != itsParamNames.end(); ++it)
	{
		os << indent;
		if (GetParam(*it).GetType() == GA_STRING)
		{
			os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15)<< "STRING";
		}
		else if (GetParam(*it).GetType() == GA_DOUBLE)
		{
			os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15)<< "DOUBLE";
		}
		else
		{
			os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15)<< "OBJECT";
		}

		os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(20)<< GetParam(*it).GetParamName();
		os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(60)<< GetParam(*it).GetParamDescription();
		if ( GetParam(*it).GetDefaultParamValue())
		{
			if (GetParam(*it).GetDefaultParamValue()->GetType() == GA_STRING)
			{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15)<< GetParam(*it).GetDefaultParamValue()->GetString();
			}
			else if (GetParam(*it).GetDefaultParamValue()->GetType() == GA_DOUBLE)
			{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15)<< GetParam(*it).GetDefaultParamValue()->GetDouble();
			}
			else
			{
				os << CC_NS(std,setw)(15) << "DEF";
			}
		}
		else
		{
			os << CC_NS(std,setfill(' ')) << CC_NS(std,setw)(15) << "";
		}

		os << endl;
	}

	return os.str();
}

void ARM_GenericParams::SetParamValue(const string& name, const ARM_GenericParamValue& value)
{
	string nameUC(name);
	stringToUpper(nameUC);
	stringTrim(nameUC);

	GenericParamValueMap::iterator it = itsMap.find(nameUC);

	if (it == itsMap.end())
	{
		pair<string,ARM_GenericParamValue> p(nameUC, value);
		itsMap.insert(p);		
	}
	else
	{
		it->second = value;
	}
}

const ARM_GenericParamValue& ARM_GenericParams::GetParamValue(const string& name) const
{
	string nameUC(name);
	stringToUpper(nameUC);
	stringTrim(nameUC);

	GenericParamValueMap::const_iterator it = itsMap.find(nameUC);

	if (it == itsMap.end())
	{
		CC_Ostringstream os;

		os << "Missing Parameter " << name << "in "<< itsFunctionName << ".";

		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());;
	}

	return it->second;
}

bool ARM_GenericParams::IsExists(const string& name) const
{
	string nameUC(name);
	stringToUpper(nameUC);

	GenericParamValueMap::const_iterator it = itsMap.find(nameUC);

	if (it == itsMap.end())
		return false;
	else
		return true;
}


ARM_GenericAddinDesc::ARM_GenericAddinDesc(const GenericAddinStruct& addinStruct)
: itsFunctionName(addinStruct.FunctionName),
itsFunctionDescription(addinStruct.FunctionDescription),
itsParams(addinStruct.FunctionName, addinStruct.params),
itsFunctor(addinStruct.Functor),
itsRetObj(addinStruct.RetObj)
{
}

string ARM_GenericAddinDesc::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "**************************************************************************************************************" << endl;
	os << indent << "** Function Name : " << GetFunctionName() << endl;
	os << indent << "** Function Description : " << GetFunctionDescription() << endl;
	os << indent << "**************************************************************************************************************" << endl;
	os << itsParams.toString(indent);

	return os.str();
}


ARM_GenericAddinsDesc* ARM_GenericAddinsDesc::itsAddinsTable = NULL;

ARM_GenericAddinsDesc::ARM_GenericAddinsDesc(const GenericAddinStruct* addinStruct, size_t nbAddins)
{
	for (int i = 0; i < nbAddins; ++i)
	{
		string nameUC(addinStruct[i].FunctionName);
		stringToUpper(nameUC);

		ARM_GenericAddinDesc addinDesc(addinStruct[i]);

		pair<string,ARM_GenericAddinDesc> p(nameUC,addinDesc);

		itsMap.insert(p);
	}
}

ARM_GenericAddinDesc ARM_GenericAddinsDesc::GetAddin(const string& name)
{
	string nameUC(name);
	stringToUpper(nameUC);

	GenericAddinMap::const_iterator it = itsMap.find(nameUC);

	if (it == itsMap.end())
	{
		CC_Ostringstream os;

		os << "The addin " << name << " is missing";

		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}

	return it->second;
}

ARM_GenericAddinsDesc ARM_GenericAddinsDesc::ExtractAddinsDesc(const string& name)
{
	ARM_GenericAddinsDesc addinsDesc;

	string nameUC(name);
	stringToUpper(nameUC);

	pair<string,ARM_GenericAddinDesc> p(nameUC,GetAddin(name));

	addinsDesc.itsMap.insert(p);

	return addinsDesc;
}

void ARM_GenericAddinsDesc::CreateTheAddingTable()
{
	itsAddinsTable = new ARM_GenericAddinsDesc(GenericAddinsTable,sizeof(GenericAddinsTable)/sizeof(GenericAddinStruct));
}

void ARM_GenericAddinsDesc::ReleaseTheAddingTable()
{
	delete itsAddinsTable;
}


const string AdditionalHelpText[] = 
{
	"",
	"",
	"In this file we describe how to use the generic addin.",
	"To call a function you should use the addin ARM_GP_GenericAddin.",
	"It has 3 parameters:",
	"FunctionName: name of the function",
	"ParamNames: name of the parameters",
	"ParamValues: values of the parameters",
	"",
	"You can use the addin ARM_GP_GenericAddin_Helper to get information about",
	"any function.",
	"With the addin ARM_GP_GenericAddin_ParamNames you can get the list of all ",
	"the parameter of the function."
	"",
	"",
};

string ARM_GenericAddinsDesc::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	int i;
	int AdditionalHelpTextSize = sizeof(AdditionalHelpText)/sizeof(AdditionalHelpText[0]);

	string space = " ";

	for(i=0; i<AdditionalHelpTextSize; ++i)
		os  << space << AdditionalHelpText[i] << "\n";

	GenericAddinMap::const_iterator it;

	for (it = itsMap.begin(); it != itsMap.end(); ++it)
	{
		os << it->second.toString(space);
		os << endl;
	}

	return os.str();
}

void ARM_GenericAddinFunctor::ResizeRetValues(int nbRows, int nbCols)
{
	itsNbRows = nbRows;
	itsNbCols = nbCols;

	itsRetValues.resize(nbRows*nbCols);
}

void ARM_GenericAddinFunctor::SetValue(int row, int col, double val)
{
	RetStruct ret;

	ret.type = DOUBLE;
	ret.dblVal = val;
	ret.strVal = "";

	itsRetValues[itsNbRows*row+col] = ret;
}

void ARM_GenericAddinFunctor::SetValue(int row, int col, const string& str)
{
	RetStruct ret;

	ret.type = STRING;
	ret.dblVal = 0.0;
	ret.strVal = str;

	itsRetValues[itsNbRows*row+col] = ret;
}

long ARMLOCAL_GenericAddin_Helper(
		const string& functionName,
		ARM_result&	result, 
        long objId)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GenericAddinsDesc* genericAddinsDesc = NULL;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Generic Addin Helper" );

		if (functionName != "")
		{
			genericAddinsDesc = static_cast<ARM_GenericAddinsDesc*>(ARM_GenericAddinsDesc::TheAddinsTable()->ExtractAddinsDesc(functionName).Clone());
		}
		else
		{
			genericAddinsDesc = static_cast<ARM_GenericAddinsDesc*>(ARM_GenericAddinsDesc::TheAddinsTable()->Clone());
		}
	
		// assign object
		if ( !assignObject( genericAddinsDesc, result, objId ) )
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

long ARMLOCAL_GenericAddin_ParamNames(
		const string& functionName,
		bool withDefault,
		vector<string>& paramNames,
		ARM_result&	result)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GenericAddinsDesc* genericAddinsDesc = NULL;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Generic Addin Param Names" );

		ARM_GenericAddinDesc addin = ARM_GenericAddinsDesc::TheAddinsTable()->GetAddin(functionName);
	
		paramNames = addin.GetParams().GetParamNames(withDefault);
		
		return ARM_OK; 
		
	}
	
	catch(Exception& x)
	{	
		x.DebugPrint();
		ARM_RESULT();
	}
}
