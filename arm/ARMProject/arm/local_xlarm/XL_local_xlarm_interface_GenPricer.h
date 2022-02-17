/////////////////////////////////////////
///
/// Copyright (c) CDC IXIS CM July 2003 Paris
///
/// \brief contains all the addins for the generic pricer
/// \author E. Benhamou, JM Prie, E. Mostafa Ezzine
/// this is part of a general table that contains all
/// the excel addins
///
///	\version 1.0
///	\date December 2003
//////////////////////////////////////////


	/// Gen Pricer Section
{
        	" Local_GenericAddin",	/// name of the C++ function
            " RRRR",				/// 4 parametres = 3 d'entree + 1 de retour 
            " ARM_GP_GenericAddin_Create",
 			" FunctionName,ParamNames,ParamValues",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Call the generic addin.",
            " Function Name",
			" Param Names",
			" Param Values"
    },
	{
        	" PXL_Local_GenericAddin",	/// name of the C++ function
            " RRRR",				/// 4 parametres = 3 d'entree + 1 de retour 
            " PXL_GP_GenericAddin_Create",
 			" FunctionName,ParamNames,ParamValues",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Call the generic addin.",
            " Function Name",
			" Param Names",
			" Param Values"
    },
	{
        	" Local_GenericAddin_Helper",	/// name of the C++ function
            " RR",						/// 2 parametres = 1 d'entree + 1 de retour 
            " ARM_GP_GenericAddin_Helper_Create",
 			" FunctionName",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create the generic addin helper.",
            " Function Name",
    },
	{
        	" PXL_Local_GenericAddin_Helper",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 de retour 
            " PXL_GP_GenericAddin_Helper_Create",
 			" FunctionName",
            " 0",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create the generic addin helper.",
            " Function Name",
    },
	{
        	" Local_GenericAddin_ParamNames",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GP_GenericAddin_ParamNames_Create",
 			" FunctionName, WithDefault",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create the generic addin param names.",
            " Function Name",
			" With Default Parameters (N)",
    },

	{
        	" Local_GetWarning_OnObject",		/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_Warning_GetOnObj",
			" Message",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " gets a warning object from object",
			" obj",
    },
	{
        	" Local_PXL_GetWarning_OnObject",	/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_Warning_GetOnObj",
			" MEssaage",
            " 0",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " gets a warning object from object",
			" obj",
    },
	{
        	" Local_Warning_SetPopUpFlag",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_Warning_SetPopUpFlag",
			" Flag",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " set the popup warning flag",
			" boolean",
    },
	{
        	" Local_Warning_Activate",		/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_Warning_Activate",
			" Flag",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " set warning flag",
			" boolean",
    },
	{
        	" Local_DealDes_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_DealDes_Create",
			" DealDes1,[DealDes2],[PricedColumns]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a deal description object",
			" first input as a matrix of cash flows",
			" second input as a matrix of cash flows",
			" priced columns",
    },





	{
        	" Local_DealDes_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_DealDes_Create",
			" DealDes1,[DealDes2],[PricedColumns]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a deal description object",
			" first input as a matrix of cash flows",
			" second input as a matrix of cash flows",
			" priced columns",
    },
    {
        	" Local_PXL_DealDes_Create",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DealDes_Create",
			" DealDes1,[DealDes2],[PricedColumns]",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a deal description object",
			" first input as a matrix of cash flows",
			" second input as a matrix of cash flows",
			" priced columns",
    },
    {
        	" Local_GenSec_Create",			/// name of the C++ function
            " RRRRRRRRR",					/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_Create",
			" DealDes1,[DealDes2],[Discounting],[CstManagerId],[ExercBoundaryResetFlag],[OtherPayoffsFlag],[PricedColumns],[IVFlag]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
			" Creates a deal description object",
            " first input as a matrix of cash flows",
			" second input as a matrix of cash flows (Optional)",
			" Discounting curve: model name (Optional)",
			" Cst manager object (Optional)",
			" Exercise Boundary Reset Flag (Optional ; default = Y)",
			" Intermediate Payoffs & Snapshots Flag (Optional ; default = Y)"
			" priced columns",
			" Intermediate Values Flag (Optional ; default = Y)"
    },
    {
        	" Local_PXL_GenSec_Create",		/// name of the C++ function
            " RRRRRRRRR",					/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenSec_Create",
			" DealDes1,[DealDes2],[Discounting],[CstManagerId],[ExercBoundaryResetFlag],[OtherPayoffsFlag],[PricedColumns],[IVFlag]",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a generic security object",
			" first input is another generic security object",
			" second input as a matrix of cash flows (Optional)",
			" Discounting curve: model name (Optional)",
			" Cst manager object (Optional)",
			" Exercise Boundary Reset Flag (Optional ; default = Y)",
			" Intermediate Payoffs & Snapshots Flag (Optional ; default = Y)"
			" priced columns",
			" Intermediate Values Flag (Optional ; default = Y)"
    },
    {
        	" Local_GenSec_ChangeAmericanIntoTrigger",			/// name of the C++ function
            " RRR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_ChangeAmericanIntoTrigger",
			" GenSec,PricingModel",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a trigger generic security object from an american one",
			" first input is a generic security object", 
			" second input is a pricing model "
    },
    {
        	" Local_PXL_GenSec_ChangeAmericanIntoTrigger",		/// name of the C++ function
            " RRR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenSec_ChangeAmericanIntoTrigger",
			" GenSec,PricingModel",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a trigger generic security object from an american one",
			" first input is another generic security object", 
			" second input is a pricing model "
    },
	{
        	" Local_ChangeSecurityIntoGenSec",			/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_ChangeSecurityIntoGenSec",
			" SecurityId,AsOfDate,ModelName",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a generic security object from a std ARM Security",
			" security object", 
			" AsOfDate",
			" ModelName "
    },
    {
        	" Local_PXL_ChangeSecurityIntoGenSec",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_ChangeSecurityIntoGenSec",
			" SecurityId,AsOfDate,ModelName",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a generic security object from a std ARM Security",
			" security object", 
			" AsOf",
			" ModelName "
    },
    {
        	" Local_GenSec_ChangeMAXPVIntoExercise",			/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_ChangeMAXPVIntoExercise",
			" GenSec",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a new generic security from another one.",
			" first input is another generic security object. ", 
			" Its MAX(PV), keyworkds will be changed to Exercise keywords."
    },
    {
        	" Local_PXL_GenSec_ChangeMAXPVIntoExercise",		/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenSec_ChangeMAXPVIntoExercise",
			" GenSec",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a new generic security from another one.",
			" first input is another generic security object. ", 
			" Its MAX(PV), keyworkds will be changed to Exercise keywords."
    },
	{
        	" Local_DealDes_ExtractSubDealDes",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_ExtractSubDealDes",
			" GenSec,ColName,[Cols],[ccy]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Create a sub generic security based on the name of the last column",
			" generic security Id",
			" column Name",
			" multi column names",
			" ccy",
    },
    {
        	" Local_PXL_DealDes_ExtractSubDealDes",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenSec_ExtractSubDealDes",
			" GenSec,ColName,[Cols],[ccy]",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Create a sub generic security based on the name of the last column",
			" generic security Id",
			" column Name",
			" multi column names",
			" ccy",
    },
    {
        	" Local_GenSec_SetPTFlag",		/// name of the C++ function
            " RRR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_SetPTFlag",
			" generic security, parse tree flag",
            " 1",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " set on and off the parse tree flag on the generic security",
			" generic security object",
			" boolean flag",
    },
    {
        	" Local_GenSec_GetDealDesTable",    /// name of the C++ function
            " RR",							    /// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_GetDealDesTable",
			" Generic Security Id",
            " 1",							    /// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
			" Get the deal description table of a generic security",
            " Generic security object",
    },
	{
        	" Local_GenSec_GetCstManager",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GenSec_GetCstManager",
			" GenSec Id",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Extract the cst manager from a generic security",
			" Generic Security Id",
    },
    {
        	" Local_PXL_GenSec_GetCstManager",	/// name of the C++ function
            " RR",							    /// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenSec_GetCstManager",
			" GenSec Id",
            " 0",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Extract the cst manager from a generic security",
			" Generic Security Id",
    },
    {
        	" Local_IRFwdMod_Create",		/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_IRFwdMod_Create",
			" ZeroCurve",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an interest rate fwd model",
			" Object",
    },
    {
        	" Local_PXL_IRFwdMod_Create",		/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_IRFwdMod_Create",
			" ZeroCurve",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an interest rate fwd model",
			" Object",
    },
    {
        	" Local_InfFwdMod_Create",		/// name of the C++ function
            " RR",							/// 3 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_InfFwdMod_Create",
			" InfCurve",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an inflation fwd model",
			" Object",
    },
    {
        	" Local_PXL_InfFwdMod_Create",		/// name of the C++ function
            " RR",							/// 3 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_InfFwdMod_Create",
			" InfCurve",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an inflation fwd model",
			" Object",
    },
    {
        	" Local_InflationEquityModel_Create",		/// name of the C++ function
            " RRRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_InflationEquityModel_Create",
			" InfCurve, PublicationLag, ModelParam1, ModelParam2",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an inflation equity model",
			" Object (Zc Curve)",
			" Publication Lag (in days)",
			" Object (Volatility or Multiplier type)",
			" Object (Volatility or Multiplier type)",
    },
    {
        	" Local_InflationEquityModel_Create",		/// name of the C++ function
            " RRRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_InflationEquityModel_Create",
			" InfCurve, PublicationLag, ModelParam1, ModelParam2",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an inflation equity model",
			" Object (Zc Curve)",
			" Publication Lag (in days)",
			" Object (Volatility or Multiplier type)",
			" Object (Volatility or Multiplier type)",
    },
	{
        	" Local_SABREquityModel_Create",		/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_SABREquityModel_Create",
			" ZcCurve,Spot,ModelParams",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SABR equity model",
			" Object (Zc Curve)",
			" double (spot value)",
			" Vector of model params (Alpha,Beta,Rho,VolOfVol,Div)",
    },
	{
        	" Local_PXL_SABREquityModel_Create",		/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_SABREquityModel_Create",
			" ZcCurve,Spot,ModelParams",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SABR equity model",
			" Object (Zc Curve)",
			" double (spot value)",
			" Vector of model params (Alpha,Beta,Rho,VolOfVol,Div)",
    },
    {
        	" Local_BINumMethod_Create",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_BiNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a backward induction method object",
			" number of steps",
			" truncation policy"
    },
    {
        	" Local_PXL_BINumMethod_Create",/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BiNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a backward induction method object",
			" number of steps",
			" truncation policy"
    },
    {
        	" Local_FINumMethod_Create",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_FiNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a forward induction method object",
			" number of steps",
			" truncation policy"
    },
    {
        	" Local_PXL_FINumMethod_Create",/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_FiNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a forward induction method object",
			" number of steps",
			" truncation policy"
    },
	{
        	" Local_MixteNumMethod_Create",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MixteNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a mixte induction method object",
			" number of steps",
			" truncation policy"
    },
    {
        	" Local_PXL_MixteNumMethod_Create",/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MixteNum_Create",
			" [Steps Nb],[Truncation Policy]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a mixte induction method object",
			" number of steps",
			" truncation policy"
    },
    {
        	" Local_SetNumMethodtoModel",	/// name of the C++ function
            " RRR",							/// 3 parametres = 1 d'entree + 2 parametre de retour 
            " ARM_GP_SetNumMethodtoModel",
			" Pricing Mod, Numerical Method",
            " 1",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Sets the numerical method of the pricing model",
			" Object (Generic pricing model)",
			" Object (Numerical method)",
    },
    {
        	" Local_PXL_SetNumMethodtoModel",/// name of the C++ function
            " RRR",							/// 3 parametres = 1 d'entree + 2 parametre de retour 
            " PXL_ARM_GP_SetNumMethodtoModel",
			" Pricing Mod, Numerical Method",
            " 0",							/// not visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Sets the numerical method of the pricing model",
			" Object (Generic pricing model)",
			" Object (Numerical method)",
    },
    {
        	" Local_SetNumerairetoModel",	/// name of the C++ function
            " RRR",							/// 3 parametres = 1 d'entree + 2 parametre de retour 
            " ARM_GP_SetNumerairetoModel",
			" Pricing Mod, Numeraire",
            " 1",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Sets the numeraire of the pricing model",
			" Object (Generic pricing model)",
			" Object (Numeraire)",
    },
    {
        	" Local_PXL_SetNumerairetoModel",/// name of the C++ function
            " RRR",							/// 3 parametres = 1 d'entree + 2 parametre de retour 
            " PXL_ARM_GP_SetNumerairetoModel",
			" Pricing Mod, Numeraire",
            " 0",							/// not visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Sets the numeraire of the pricing model",
			" Object (Generic pricing model)",
			" Object (Numeraire)",
    },
    {
        	" Local_GetNumMethodFromModel",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GetNumMethodFromModel",
			" Pricing Mod",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Get the numerical method of a pricing model",
			" Object (Generic pricing model)",
    },
    {
        	" Local_PXL_GetNumMethodFromModel",/// name of the C++ function
            " RR",							   /// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GetNumMethodFromModel",
			" Pricing Mod",
            " 0",							/// not visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Get the numerical method of a pricing model",
			" Object (Generic pricing model)",
    },
    {
        	" Local_GramHelper_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_Helper",
			" [FuncName]",
            " 1",							/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Help information on generic pricer",
			" Function Name (case insensitive) as a string",
    },
    {
        	" Local_PXL_GramHelper_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Helper",
			" [FuncName]",
            " 0",							/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Help information on generic pricer",
			" Function Name (case insensitive) as a string",
    },
    {
        	" Local_GramHelper_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_Helper_Create",
			" [FuncName]",
            " 1",							/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Help information on generic pricer",
			" Function Name (case insensitive) as a string",
    },
    {
        	" Local_PXL_GramHelper_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Helper_Create",
			" [FuncName]",
            " 0",							/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Help information on generic pricer",
			" Function Name (case insensitive) as a string",
    },
    {
        	" Local_Numeraire_Create",		/// name of the C++ function
            " RRR",						    /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Numeraire_Create",
			" NumeraireType, [Time lags Vector]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Numeraire",
			" String for numeraire type : Cash,TerminalZc,TerminalEventZC,RollingEvent,RollingPayment",
			" Vector of time lags (optional)",
    },
    {
        	" Local_PXL_Numeraire_Create",  /// name of the C++ function
            " RRR",					    	/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Numeraire_Create",
			" Numeraire Type, [Time Lags Vector]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Numeraire",
			" String for numeraire type : Cash,TerminalZc,,TerminalEventZC,RollingEvent,RollingPayment",
			" Vector of time lags (optional)",
    },
    {
        	" Local_ModelParam_Create",		/// name of the C++ function
            " RRRRRRRRRR",					/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_ModelParam_Create",
			" Param Type, Timelags,Values, [Param Name],[LoweBoundary],[UpperBoundary],[InterpolMethod],[AdviseBreakPoint],[Currency]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Model Parameter",
			" string: Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread,Correlation, Beta, Shift",
			" vector of time lags",
			" vector of values",
			" string for param name (optional)",
			" vector of values (optional vector = NULL)",
			" vector of values (optional vector = NULL)",
            " string (optional method = STEPUPRIGHT),LINEAR, STEPUPLEFT,STEPUPRIGHT  ",
			" boolean true or false (default=false)",
			" string for currency",
    },
    {
        	" Local_PXL_ModelParam_Create", /// name of the C++ function
            " RRRRRRRRRR",					/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ModelParam_Create",
			" Param Type, Timelags, Values, [Param Name],[LoweBoundary],[UpperBoundary],[InterpolMethod],[AdviseBreakPoint],[Currency]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Model Parameter",
			" string: Volatility,MeanReversion,Beta,...please check the correct param type before building princing model",
			" vector of time lags",
			" vector of values",
			" vtring for param name (optional)",
			" vector of values (optional vector = NULL)",
			" vector of values (optional vector = NULL)",
            " string (optional method = STEPUPRIGHT),LINEAR, STEPUPLEFT,STEPUPRIGHT  ",
			" boolean true or false (default=false)",
			" string for currency",
    },
    {
        	" Local_HW1FModel_Create",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_HW1FModel_Create",
			" Zero Curve,Model Param1,Model Param2,[flags]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull White 1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility or MeanReversion type)",
			" Object (Volatility or MeanReversion type)",
			" Flags : [0]=ApproxSO(1)/ExactSO(0), [1]=NoITMOptim(1)/ITMOptim(0)",
    },
    {
        	" Local_PXL_HW1FModel_Create",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HW1FModel_Create",
			" Zero Curve,Model Param1,Model Param2,[flags]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull White 1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility or MeanReversion type)",
			" Object (Volatility or MeanReversion type)",
			" Flags : [0]=ApproxSO(1)/ExactSO(0), [1]=NoITMOptim(1)/ITMOptim(0)",
    },
    {
        	" Local_HW2FModel_Create",		/// name of the C++ function
            " RRRRRRRR",					/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_HW2FModel_Create",
			" Zero Curve, Model Param1,Model Param2,Model Param3,Model Param4,Model Param5,[flags]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull White 2F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Flags : [0]=ApproxSO(1)/ExactSO(0), [1]=NoITMOptim(1)/ITMOptim(0)",
    },
    {
        	" Local_PXL_HW2FModel_Create",	/// name of the C++ function
            " RRRRRRRR",					/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HW2FModel_Create",
			" Zero Curve, Model Param1,Model Param2,Model Param3,Model Param4,Model Param5,[flags]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull White 2F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Object (Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread or Correlation type)",
			" Flags : [0]=ApproxSO(1)/ExactSO(0), [1]=NoITMOptim(1)/ITMOptim(0)",
    },
    {
        	" Local_MF1FModel_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_MarkovFunctionalModel_Create",
			" ZeroCurve,ModelParam1,[ModelParam2]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a 1F Markov Functional Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility or MR type)",
			" Object (Volatility or MR type)",
    },
    {
        	" Local_PXL_MF1FModel_Create",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MarkovFunctionalModel_Create",
			" ZeroCurve,ModelParam1,[ModelParam2]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create a 1F Markov Functional Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility or MR type)",
			" Object (Volatility or MR type)",
    },
    {
        	" Local_OldTreeMethod_Create",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_OldTreeMethod_Create",
			" Time Step Nb, [Trunc Ratio], [Min StdDev], [Min Step Nb]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree method object (1st generation)",
			" Number of time steps",
			" Number of global StdDev for truncation (optional, default=5)",
			" Min annual local StdDev for node reduction (optional, default=1e-3)",
			" Min number of time steps before 1st event date (optional, default=0)"
    },
    {
        	" Local_PXL_OldTreeMethod_Create",	 /// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_OldTreeMethod_Create",
			" Time Step Nb, [Trunc Ratio], [Min StdDev], [Min Step Nb]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree method object (1st generation)",
			" Number of steps",
			" Number of global StdDev for truncation (optional, default=5)",
			" Min annual local StdDev for node reduction (optional, default=1e-3)",
			" Min number of time steps before 1st event date (optional, default=0)"
    },
    {
        	" Local_TreeMethod_Create",	    /// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_TreeMethod_Create",
			" Time Step Nb, [Trunc Ratio], [Min StdDev], [Min Step Nb]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree method or Tree ND object depending of internal flag",
			" Number of time steps",
			" Number of global StdDev for truncation (optional, default=5)",
			" Min annual local StdDev for node reduction (optional, default=1e-3)",
			" Min number of time steps before 1st event date (optional, default=0)"
    },
    {
        	" Local_PXL_TreeMethod_Create",	    /// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_TreeMethod_Create",
			" Time Step Nb, [Trunc Ratio], [Min StdDev], [Min Step Nb]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree method or Tree ND object depending of internal flag",
			" Number of steps",
			" Number of global StdDev for truncation (optional, default=5)",
			" Min annual local StdDev for node reduction (optional, default=1e-3)",
			" Min number of time steps before 1st event date (optional, default=0)"
    },
    {
        	" Local_Tree_SetProbaFlag",	        	/// name of the C++ function
            " RRR",							        /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Tree_SetProbaFlag",
			" Tree method, proba computation flag",
            " 1",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " set on and off the proba computation flag on the tree method",
			" Tree object",
			" boolean flag",
    },
	{
        	" Local_GenSec_SetCRFlag",				/// name of the C++ function
            " RRR",							        /// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_SetCRFlag",
			" Circular Ref flag",
            " 1",									/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " set flag to look for circular reference",
			" circular reference flag (true or false)"
    },
    {
        	" Local_SFRMModel_Create",		/// name of the C++ function
            " RRRRRRRR",					/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_SFRMModel_Create",
			" Zero Curve, Model Params,volType, factorsNb,IRIndex,[ShiftConvPort],[NonParamDrift]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SFRM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Model Params Objects (Volatility,MeanReversion,Correlation,Shift)",
			" volType either ROW or DIAG",
			" factorsNb = 1..3",
			" IRIndex object Id"
			" Portfolio used to convert beta to shift using the fwds",
			" Non Parametric Diffusion (Y,N))"
    },
	{
        	" Local_SFRMModelVolSwapVolFRADump",		/// name of the C++ function
            " RRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_SFRMModel_VolSwapVolFRADump",
			" Swaption, SFRM Model",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Dump a Vol Swap Vol FRA information from a swaption",
			" Swaption Id",
			" SFRM Id"
    },
    {
        	" Local_PXL_SFRMModel_Create",		/// name of the C++ function
            " RRRRRRRR",						/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SFRMModel_Create",
			" Zero Curve, Model Params,volType, factorsNb,IRIndex,[ShiftConvPort],[NonParamDrift]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SFRM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Model Params Objects (Volatility,MeanReversion,Correlation,Shift)",
			" volType either ROW or DIAG",
			" factorsNb = 1..3",
			" IRIndex object Id",
			" Portfolio used to convert beta to shift using the fwds",
			" Non Parametric Diffusion (Y,N))"
    },
	{
        	" Local_SetSFRMFixScheduler",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_SetSFRMFixScheduler",
			" SFRMModId, StartDate, EndDate",
            " 1",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Set a fix scheduler on a SFRM Model.",
			" SFRM Model Id",
			" Start Date",
			" End Date",
    },
    {
        	" Local_PXL_SetSFRMFixScheduler",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SetSFRMFixScheduler",
			" SFRMModId, StartDate, EndDate",
            " 0",							/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Set a fix scheduler on a SFRM Model.",
			" SFRM Model Id",
			" Start Date",
			" End Date",
    },
	{
        	" Local_TrigoMatrix",		/// name of the C++ function
            " RRRR",					/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_TrigoMatrix",
			" N, Alpha",
            " 1",						/// visible in excel
            XLLOCALARM_GLOBAL_GROUP,
            " ",
            " ",
            " Computes the trigonometric matrix",
			" size",
			" constant",
			
    },
    {
        	" Local_LocalCovariance",	/// name of the C++ function
            " RRRRRRRRRRRRR",			/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_LocalCovariance",
			" Model, underlyingType, FromDate, ToDate, StartDate1, EndDate1, StartDate2, EndDate2, [StartDate3], [EndDate3], [StartDate4], [EndDate4]",
            " 1",						/// visible in excel
            XLLOCALARM_GLOBAL_GROUP,
            " ",
            " ",
            " Computes the local covariance (From -> To) of 2 underlyings",
			" Model (HW, SFRM, 2IR+FX)",
			" underlying type (ZC, FWD, CMS, FX, SO/CMS)",
			" Integration From Date",
			" Integration To Date",
			" Underlying #1 Start Date",
			" Underlying #1 End Date",
			" Underlying #2 Start Date",
			" Underlying #2 End Date",
			" If SO = S1-S2 vs CMS : Underlying #3 Start Date",
			" If SO = S1-S2 vs S3 : Underlying #3 End Date",
			" If SO = S1-S2 vs S3-S4 : Underlying #4 Start Date",
			" If SO = S1-S2 vs S3-S4 : Underlying #4 End Date"
    },
    {
        	" Local_LocalCorrelation",	/// name of the C++ function
            " RRRRRRRRRRRRR",			/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_LocalCorrelation",
			" Model, underlyingType, FromDate, ToDate, StartDate1, EndDate1, StartDate2, EndDate2, [StartDate3], [EndDate3], [StartDate4], [EndDate4]",
            " 1",						/// visible in excel
            XLLOCALARM_GLOBAL_GROUP,
            " ",
            " ",
            " Computes the local correlation (From -> To) of 2 underlyings",
			" Model (HW or SFRM)",
			" underlying type (ZC, FWD, CMS, SO/CMS, SO/SO)",
			" Integration From Date",
			" Integration To Date",
			" Underlying #1 Start Date",
			" Underlying #1 End Date",
			" Underlying #2 Start Date",
			" Underlying #2 End Date",
			" If SO = S1-S2 vs S3 : Underlying #3 Start Date",
			" If SO = S1-S2 vs S3 : Underlying #3 End Date",
			" If SO = S1-S2 vs S3-S4 : Underlying #4 Start Date",
			" If SO = S1-S2 vs S3-S4 : Underlying #4 End Date"
	},
    {
        	" Local_MCMethod_Create",		/// name of the C++ function
            " RRRRRRRRRRRRR",				/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_MCMethod_Create",
			" [Iteration Nb], [Fix Time Step],[RandGenId],[SamplerType],[SamplerDatas],[SchedulerType],[ScheduerDatas],[ExerciseBoundaryCalculator],[Max Paths per Bucket],[ImpSampler],[ImpSamplerDatas],[PathScheme]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the MC Method",
			" Nb of iteration (Default=2)",
			" Nb of fix time step in days for Euler scheme (Default=1)",
			" object default MRGK5 BoxMuller Antithetic",
			" Sampler type, default=NormalCentred",
			" Sampler Datas [0]=MinStdDev",
			" Scheduler type, default=TimeStepPerYear",
			" Scheduler Datas, default time step=1",
			" exercise boundary calculator",
			" Bucket Size limit ",
			" Importance Sampler type, default=Dummy",
			" Importance Sampler Datas (Empty)",
			" Path Scheme type, default=Incremental",
    },
    {
        	" Local_PXL_MCMethod_Create",	/// name of the C++ function
            " RRRRRRRRRRRRR",				/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MCMethod_Create",
			" [Iteration Nb],[Fix Time Step],[RandGenId],[SamplerType],[SamplerDatas],[SchedulerType],[ScheduerDatas],[ExerciseBoundaryCalculator],[Max Paths per Bucket],[ImpSampler],[ImpSamplerDatas],[PathScheme]",
            " 0",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the MC Method",
			" Nb of iteration (Default=2)",
			" Nb of fix time step in days for Euler scheme (Default= 1)",
			" object default MRGK5 BoxMuller Antithetic",
			" Sampler type, default=NormalCentred",
			" Sampler Datas (Empty)",
			" Scheduler type, default=TimeStepPerYear",
			" Scheduler Datas, default time step=1",
			" exercise boundary calculator",
			" Bucket Size limit ",
			" Importance Sampler type, default=Dummy",
			" Importance Sampler Datas (Empty)",
			" Path Scheme, default=Incremental",
    },
    {
        	" Local_PDEMethod_Create",		/// name of the C++ function
            " RRRRRRRRRRRRR",							/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_PDEMethod_Create",
			" [MethodName], [TimeStepsNb], [SpaceStepsNb], [StdDevNb], [YGridNb], [ZGridNb], [Theta1], [Theta2], [Theta3], [BCName], [Lambda], [GridType]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the PDE Method",
			" Method Name : Default = Explicit",
			" Nb of Additional Discr Time Steps which can be extend to Matrix of Scheduler (Default = 0)",
			" Nb of Discr Space Steps (Default = 101)",
			" Nb Of StdDev for Space Extension whcih can be extend to Matrix of GridData (Default = 6)",
			" Nb of YGrid Steps (Default = 3)",
			" Nb of ZGrid Steps (Default = 3)",
			" Theta1: Implicitness of the Diffusion term (Theta1=0 is the explicit case) (Default = 0.5)",
			" Theta2: Implicitness of the Convection term (Theta2=0 is the explicit case) (Default = 0.5)",
			" Theta3: Implicitness of the Actualisation term (Theta3=0 is the explicit case) (Default = 0.5)",
			" Boundary Condition Name (Default = VonNeumann)",
			" Lambda: PDE switch, ie number which enables to switch of PDE according to the model",
			" GridType: Fixed for choosing the extreme value, StdDev for letting the model do it, BiReg for an StdDev with 2 regims (Default=StdDev)"
    },
    {
        	" Local_PXL_PDEMethod_Create",	/// name of the C++ function
            " RRRRRRRRRRRRR",							/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_PDEMethod_Create",
			" [MethodName], [TimeStepsNb], [SpaceStepsNb], [StdDevNb], [YGridNb], [ZGridNb], [Theta1], [Theta2], [Theta3], [BoundaryConditionName], [Lambda], [GridType]",
            " 0",						/// visible in excel
     XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the PDE Method",
			" Method Name : Default = Explicit",
			" Nb of Additional Discr Time Steps (Default = 0)",
			" Nb of Discr Space Steps (Default = 100)",
			" Nb Of StdDev for Space Extension (Default = 6)",
			" Nb of YGrid Steps (Default = 2)",
			" Nb of ZGrid Steps (Default = 2)",
			" Theta1: Implicitness of the Diffusion term (Theta1=0 is the explicit case) (Default = 0.5)",
			" Theta2: Implicitness of the Convection term (Theta2=0 is the explicit case) (Default = 0.5)",
			" Theta3: Implicitness of the Actualisation term (Theta3=0 is the explicit case) (Default = 0.5)",
			" Boundary Condition Name : Default = VonNeumann",
			" Lambda: Weight of the correction term in the Predictor-Corrector (Lambda=0 is the non corrected scheme) (Default = 0.0)",
			" GridType (Default=StdDev)"
    },

	{
        	" Local_PdeND_Create",		/// name of the C++ function
            " RRRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_PdeND_Create",
			" MethodName, [SchedulerType], [SchedulerData], [SpacesData], [SchemeData], [BoundCond]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Pde ND method object",
			" Method Name",
			" Scheduler type, default=ConstantVarianceMeanRevertingScheduler",
			" Scheduler Datas",
			" Nb Of StdDev for Space Extension, a matrix with as many columns as the nb of factors",
			" Scheme Datas: matrix 4*1 with M11=Implicitness of the Diffusion term (default=0.5), M21=Implicitness of the Convection term(default=0.5), M31: Implicitness of the Actualisation term(default=0.5), M41: Weight of the correction term in the Predictor-Corrector, default=0 is the non corrected scheme)",
			" Bound Cond Name: default = VonNeumann",
			
    },
    {
        	" Local_PXL_PdeND_Create",		/// name of the C++ function
            " RRRRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_PdeND_Create",
			" MethodName, [SchedulerType], [SchedulerData], [SpacesData], [SchemeData], [BoundCond]",
            " 0",								/// visible in excel
     XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the PDE Method",
			" Method Name : Default = Explicit",
			" Nb of Additional Discr Time Steps (Default = 0)",
			" Nb of Discr Space Steps (Default = 100)",
			" Sampler data"
			" Nb of YGrid Steps (Default = 2)",
			" Nb of ZGrid Steps (Default = 2)",
			" Theta1: Implicitness of the Diffusion term (Theta1=0 is the explicit case) (Default = 0.5)",
			" Theta2: Implicitness of the Convection term (Theta2=0 is the explicit case) (Default = 0.5)",
			" Theta3: Implicitness of the Actualisation term (Theta3=0 is the explicit case) (Default = 0.5)",
			" Boundary Condition Name : Default = VonNeumann",
			" Lambda: Weight of the correction term in the Predictor-Corrector (Lambda=0 is the non corrected scheme) (Default = 0.0)"
    },




    {
            " RRRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_PDEMethod_Create",
			" [MethodName], [Discr TimeSteps Nb], [Discr SpaceSteps Nb]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the PDE Method",
			" Method Name : Default = Explicit",
			" Nb of Additional Discr Time Steps (Default = 0)",
			" Nb of Discr Space Steps (Default = 100)",
    },
    {
        	" PXL_Local_PXL_PDEMethod_Create",	/// name of the C++ function
            " RRRR",						/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_PDEMethod_Create",
			" [MethodName], [Discr TimeSteps Nb], [Discr SpaceSteps Nb]",
            " 0",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the PDE Method",
			" Method Name : Default = Explicit",
			" Nb of Additional Discr Time Steps (Default = 0)",
			" Nb of Discr Space Steps (Default = 100)",
    },
    {
        	" Local_GenPricer_Create",	/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_GenPricer_Create",
			" GenSecurity,Pricing Model,[CVColumnNames],[CVPrices],[RefColumn],[Beta]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a generic pricer",
			" Generic Security",
			" Pricing Model",
			" Control Variate Column Names",
			" Control Variate Reference Prices",
			" Ref Column Name",
			" Vector of regression values",
    },
    {
        	" Local_PXL_GenPricer_Create",	/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenPricer_Create",
			" GenSecurity,Pricing Model,[CVColumnNames],[CVPrices],[RefColumn],[Beta]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a generic pricer",
			" Generic Security",
			" Pricing Model",
			" Control Variate Column Names",
			" Control Variate Reference Prices",
			" Ref Column Name",
			" Vector of regression values",
    },
    {
        	" Local_GetData_FromGenPricer",	/// name of the C++ function
            " RRRRR",						/// 4 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_GenPricer_GetData",
			" Gen Pricer, key, column name, int",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get data from pricer info of a generic pricer",
			" Generic Pricer Object",
			" string",
			" string",
			" int",
    },
    {
        	" Local_PXL_GetData_FromGenPricer",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GenPricer_GetData",
			" Gen Pricer, key, column name, int",
            " 0",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get data from pricer info of a generic pricer",
			" Generic Pricer Object",
			" string",
			" string",
			" int",

    },
    {
        	" Local_GenPricer_SetDetailMode",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GenPricer_SetDetailMode",
			" GenPricer,flag",
            " 1",								/// visible in excel
            XLLOCALARM_SETTOID_GROUP,
            " ",
            " ",
            " Sets detail mode on a genpricer view",
			" Obj",
			" flag true or false"
    },
    {
        	" Local_PricingModel_GetModelParam",	/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_Model_GetModelParam",
			" ModelId,paramType,dataType,[index],[factorNb]",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get model param data type from a model",
			" model object",
			" string: Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread,Correlation, Beta, Alpha,.....",
			" string :values, breakPointTimes, tenors (if surface), strikes( if surfaceListModelParm)",
            " index to identify within a surface model param",
			" factorNb to identify within a MultiAsset model",
    },
    {
        	" Local_RandGen_Create",			/// name of the C++ function
            " RRRRRRRRRRRRRR",						/// 14 parametres = 13 d'entree + 1 parametre de retour 
            " ARM_GP_RandGen_Create",
			" [genType],[algoType],[baseGen1Id],[seed],[dim],[factorDim],[nbOfPaths],[NbStdDevs],[baseGen2Id],[FirstNbTimes],[Order], [firstSimulations],[FirstNbDims]",
            " 1",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " creates a random nb generator",
			" default=UnknownBaseGenAlgorithm, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=UnknownTransformAlgo, BoxMuller, InvNormCum CondInvNormCum AntitheticOne InvnormCumFast Transpose Skipper MixteGen...",
			" object",
			" number, default=-1",
			" number, default= 1",
			" number, default= 1",
			" number, default=-1",
			" nbStdDevs, default=4.0",
			" object",
			" firstNbTimes, default=0",
			" order, default=BucketOrder",
			" firstSimulations default=0",
			" firstNbDims, default=0",
    },
    {
        	" Local_PXL_RandGen_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",					/// 14 parametres = 12 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_RandGen_Create",
            " [genType],[algoType],[baseGen1Id],[seed],[dim],[factorDim],[nbOfPaths],[NbStdDevs],[baseGen2Id],[FirstNbTimes],[Order],[firstSimulations],[firstNbDims]",
            " 0",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " creates a random nb generator",
			" default=UnknownBaseGenAlgorithm, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=UnknownTransformAlgo, BoxMuller, InvNormCum CondInvNormCum AntitheticOne InvnormCumFast Transpose Skipper MixteGen...",
			" object",
			" number, default=-1",
			" number, default= 1",
			" number, default= 1",
			" number, default=-1",
			" nbStdDevs, default=4.0",
			" object",
			" firstNbDims, default=0",
			" order, default=BucketOrder",
			" firstSimulations default=0",
			" firstNbTimes, default=0",
    },
	{
        	" Local_SimpleRandomGen_Create",		/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_SimpleRandGen_Create",
			" [genType1],[genType2],[algoType1],[algoType2],[FirstNbTimes],[FirstNbDims],[IsAntithetic]",
            " 1",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " creates a random nb generator",
			" default=NR_Ran2, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=Sobol, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=BoxMuller, BoxMuller, InvNormCum CondInvNormCum",
			" default=InvNormCum, BoxMuller, InvNormCum CondInvNormCum",
			" firstNbTimes, default=0",
			" firstNbDims, default=0",
			" isAntithetic, default=Y"
    },
    {
        	" Local_PXL_SimpleRandomGen_Create",	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SimpleRandGen_Create",
			" [genType1],[genType2],[algoType1],[algoType2],[FirstNbTimes],[FirstNbDims],[IsAntithetic]",
            " 0",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " creates a random nb generator",
			" default=NR_Ran2, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=Sobol, MRGK5/3, NR_Ran1/2/3/4, Faure Halton Niederreiter Sobol Knuth Lecuyer Mersene MerseneStd ParkMiller Tausworthe...",
			" default=BoxMuller, BoxMuller, InvNormCum CondInvNormCum",
			" default=InvNormCum, BoxMuller, InvNormCum CondInvNormCum",
			" firstNbTimes, default=0",
			" firstNbDims, default=0",
			" isAntithetic, default=Y"
    },
    {
        	" Local_RandomGen_DrawVector",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_RandGen_Draw",
			" randGenId,[size]",
            " 1",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " draw a vector of given size",
            " random nb object",
			" number, default=10"
    },
    {
        	" Local_HybridBasisFwdIRModel_Create",		/// name of the C++ function
            " RRRRRR",						            /// 6 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_BasisFwdIRModel_Create",
			" IRRefModel, IRMarginZc, BasisMarginZc, Forex, ModelNames",
            " 1",							        /// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create an Hybrid Basis Forward & IR Model",
			" Object (IR Reference Model)",
			" Object (Zc Curve for IR Forward Margin Model)",
			" Object (Zc Curve for Basis Forward Margin Model)",
            " Object (Forex)",
			" String vector for model names ([0]=RefModel,[1]=IRMarginModel,[2]=BasisMarginModel)",
    },
    {
        	" Local_PXL_HybridBasisFwdIRModel_Create",	/// name of the C++ function
            " RRRRRR",						            /// 6 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BasisFwdIRModel_Create",
			" IRRefModel, IRMarginZc, BasisMarginZc, Forex, ModelNames",
            " 0",							            /// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create an Hybrid Basis Forward & IR Model",
			" Object (IR Reference Model)",
			" Object (Zc Curve for IR Forward Margin Model)",
			" Object (Zc Curve for Basis Forward Margin Model)",
            " Object (Forex)",
			" String vector for model names ([0]=RefModel,[1]=IRMarginModel,[2]=BasisMarginModel)",
    },
    {
        	" Local_GenericCurve_Create",		/// name of the C++ function
            " RRRRRR",							/// 6 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_GCurve_Create",
			" abscisses,ordinates,[interpolatorType],[sortAbscisses],[alwaysMulti]",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " creates a generic curve",
			" vector of abscisses",
			" vector of ordinates (matrix for multi curve)",
			" LINEAR, STEPUPRIGHT, STEPUPLEFT, default=LINEAR",
			" boolean to sort the abscisses, default = false",
			" boolean to always create a multi curve, default = false",
    },
    {
        	" Local_PXL_GenericCurve_Create",	/// name of the C++ function
            " RRRRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GCurve_Create",
			" abscisses,ordinates,[interpolatorType],[sortAbscisses],[alwaysMulti]",
            " 0",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " creates a generic curve",
			" vector of abscisses",
			" vector of ordinates (matrix for multi curve)",
			" LINEAR, STEPUPRIGHT, STEPUPLEFT, default=LINEAR",
			" boolean to sort the abscisses, default = false",
			" boolean to always create a multi curve, default = false",
    },
	{
        	" Local_GenericCurve_Interpolate",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GCurve_Interpolate",
			" GCurveId,abscisse",
            " 1",								/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " interpolate a generic curve",
            " generic curve object",
			" abscisse"
    },
    {
        	" Local_GenericCurve_CptCurve",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GCurve_CptCurve",
			" GCurveId,new abscisses",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " from a generic curve, get new curve based on new absicces",
            " generic curve object",
			" new abscisses"
    },
    {
        	" Local_PXL_GenericCurve_CptCurve",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GCurve_CptCurve",
			" GCurveId,new abscisses",
            " 0",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " from a generic curve, get new curve based on new absicces",
            " generic curve object",
			" new abscisses"
    },
    {
        	" Local_GenericCurve_Insert",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_GCurve_Insert",
			" GCurveId,new abscisse, new ordinate",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " insert a new value to a generic curve",
            " generic curve object",
			" new abscisse"
			" new ordonate"
    },
    {
        	" Local_PXL_GenericCurve_Insert",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GCurve_Insert",
			" GCurveId,new abscisse, new ordinate",
            " 0",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " insert a new value to a generic curve",
            " generic curve object",
			" new abscisse"
			" new ordonate"
    },
	{
        	" Local_CurveMatrix_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_CurveMatrix_Create",
			" CorrelsId, NbRows, NbCols",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " create a matrix curve",
            " vector of curves",
			" number of rows",
			" number of columns"
    },
    {
        	" Local_PXL_CurveMatrix_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CurveMatrix_Create",
			" CorrelsId, NbRows, NbCols",
            " 0",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " create a matrix curve",
            " vector of curves",
			" number of rows",
			" number of columns"
    },
	{
        	" Local_Model_ZCCurveSet",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Model_ZCCurveSet",
			" ModelId,zcCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Sets a zc curve object to an IR Model",
			" IR Model Object",
			" ZC Curve Object",
    },
	{
        	" Local_PXL_Model_ZCCurveSet",		/// name of the C++ function
            " RRRRRRRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Model_ZCCurveSet",
			" ModelId,zcCurveId",
            " 0",									/// not visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Sets a zc curve object to an IR Model",
			" IR Model Object",
			" ZC Curve Object",
    },
    {
			" Local_DateStripCreate",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",			// 19 = 18 parametre d'entree + 1 parametre de retour
            " ARM_DateStrip_Create",
            " startDate, endDate,[resetFreq],[dayCount],[resetCalendar],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCalendar],[resetTiming],[payTiming],[adjFirstDate],[refDate],[AccruedOrFull],[FirstDateFwdRule]",
            " 1",
            XLLOCALARM_UTIL_GROUP,
            " ",
            IDH_ARM_DateStripCreate,
            " Creates a date strip", 
            " date",
			" date or single maturity",
			" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default=(A)nnual",
			" ACTUAL, A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is ACTUAL",
			" ccy mainly",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
			" A S, Q, B, M, W, D, Z", 
			" in days",
			" ccy mainly",
			" ADV(ance), ARR(ears)",
			" ADV(ance), ARR(ears)",
			" true or false",
			" reference date",
            " K_ACCURUED, K_FULL, default(K_ACCURUED) : Accrud mode or full mode to generat FwdDates",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous). Default: fwdRule value",
    },
	{
        	" Local_PXL_DateStripCreate",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",			// 19 = 18 parametre d'entree + 1 parametre de retour
            " PXL_ARM_DateStrip_Create",
            " startDate, endDate,[resetFreq],[dayCount],[resetCalendar],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCalendar],[resetTiming],[payTiming],[adjFirstDate],[refDate],[AccruedOrFull],[FirstDateFwdRule]",
            " 0",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a date strip", 
            " date",
			" date or single maturity",
			" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default=(A)nnual",
			" ACTUAL, A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is ACTUAL",
			" ccy mainly"
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
			" A S, Q, B, M, W, D, Z", 
			" in days",
			" ccy mainly",
			" ADV(ance), ARR(ears)",
			" ADV(ance), ARR(ears)",
			" true or false",
			" reference date",
            " K_ACCURUED, K_FULL, default(K_ACCURUED) : Accrud mode or full mode to generat FwdDates.",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous). Default: fwdRule value",
    },
	{
			" Local_DataFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_DateStrip_GetData",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataSizeFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_DateStrip_GetDataSize",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets data size From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataSizeFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_DateStrip_GetDataSize",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_DataFromDateStrip",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DateStripVecCreate",
            " RRRRRRRRR",						// 9 parametres = 8 d'entree + 1 parametre de retour
            " ARM_DateStripFromVec",
            " FlowStartDates, FlowEndDates, FwdStartDates, FwdEndDates, ResetDates, PaymentDates, InterestDays, InterestTerms",
            " 1",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a datestrip from vector inputs",
            " FlowStartDates as a vector",
			" FlowEndDates as a vector",
            " FwdStartDates as a vector",
			" FwdEndDates as a vector",
			" ResetDates as a vector",
			" PaymentDates as a vector",
			" InterestDays as a vector",
			" InterestTerms as a vector",
	},
	{
			" Local_PXL_DateStripVecCreate",
            " RRRRRRRRR",						// 9 parametres = 8 d'entree + 1 parametre de retour
            " PXL_ARM_DateStripFromVec",
            " FlowStartDates, FlowEndDates, FwdStartDates, FwdEndDates, ResetDates, PaymentDates, InterestDays, InterestTerms",
            " 0",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a datestrip from vector inputs",
            " FlowStartDates as a vector",
			" FlowEndDates as a vector",
            " FwdStartDates as a vector",
			" FwdEndDates as a vector",
			" ResetDates as a vector",
			" PaymentDates as a vector",
			" InterestDays as a vector",
			" InterestTerms as a vector",
	},
	/// GP Version
	{
			" Local_DateStripCreate",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",			// 19 = 18 parametre d'entree + 1 parametre de retour
            " ARM_GP_DateStrip_Create",
            " startDate,endDate,[resetFreq],[dayCount],[resetCalendar],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCalendar],[resetTiming],[payTiming],[adjFirstDate],[refDate],[AccruedOrFull],[FirstDateFwdRule]",
            " 1",
            XLLOCALARM_UTIL_GROUP,
            " ",
            IDH_ARM_DateStripCreate,
            " Creates a date strip", 
            " date",
			" date or single maturity",
			" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default=(A)nnual",
			" ACTUAL, A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is ACTUAL",
			" ccy mainly",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
			" A S, Q, B, M, W, D, Z", 
			" in days",
			" ccy mainly",
			" ADV(ance), ARR(ears)",
			" ADV(ance), ARR(ears)",
			" true or false",
			" reference date",
            " K_ACCURUED, K_FULL, default(K_ACCURUED) : Accrud mode or full mode to generat FwdDates.",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous). Default: fwdRule value",
    },
	{
        	" Local_PXL_DateStripCreate",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",			// 19 = 18 parametre d'entree + 1 parametre de retour
            " PXL_ARM_GP_DateStrip_Create",
            " startDate,endDate,[resetFreq],[dayCount],[resetCalendar],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCalendar],[resetTiming],[payTiming],[adjFirstDate],[refDate],[AccruedOrFull],[FirstDateFwdRule]",
            " 0",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a date strip", 
            " date",
			" date or single maturity",
			" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default=(A)nnual",
			" ACTUAL, A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is ACTUAL",
			" ccy mainly"
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
			" A S, Q, B, M, W, D, Z", 
			" in days",
			" ccy mainly",
			" ADV(ance), ARR(ears)",
			" ADV(ance), ARR(ears)",
			" true or false",
			" reference date",
            " K_ACCURUED, K_FULL, default(K_ACCURUED) : Accrud mode or full mode to generat FwdDates.",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous). Default: fwdRule value",
    },
	{
			" Local_DataFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_GP_DateStrip_GetData",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataSizeFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_GP_DateStrip_GetDataSize",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets data size From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataSizeFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_GP_DateStrip_GetDataSize",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DataFromDateStrip",
            " RRR",							// 3 parametres = 2 d'entree + 1 parametre de retour
            " ARM_GP_DataFromDateStrip",
            " Datestrip, DataType",
            " 1",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets Data From DateStrip", 
            " Date Strip object", 
            " STARTDATE(SD), ENDDATE(ED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)"
	},
	{
			" Local_DateStripVecCreate",
            " RRRRRRRRR",						// 9 parametres = 8 d'entree + 1 parametre de retour
            " ARM_GP_DateStripFromVec",
            " FlowStartDates, FlowEndDates, FwdStartDates, FwdEndDates, ResetDates, PaymentDates, InterestDays, InterestTerms",
            " 1",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a datestrip from vector inputs",
            " FlowStartDates as a vector",
			" FlowEndDates as a vector",
            " FwdStartDates as a vector",
			" FwdEndDates as a vector",
			" ResetDates as a vector",
			" PaymentDates as a vector",
			" InterestDays as a vector",
			" InterestTerms as a vector",
	},
	{
			" Local_PXL_DateStripVecCreate",
            " RRRRRRRRR",						// 9 parametres = 8 d'entree + 1 parametre de retour
            " PXL_ARM_GP_DateStripFromVec",
            " FlowStartDates, FlowEndDates, FwdStartDates, FwdEndDates, ResetDates, PaymentDates, InterestDays, InterestTerms",
            " 0",						
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Creates a datestrip from vector inputs",
            " FlowStartDates as a vector",
			" FlowEndDates as a vector",
            " FwdStartDates as a vector",
			" FwdEndDates as a vector",
			" ResetDates as a vector",
			" PaymentDates as a vector",
			" InterestDays as a vector",
			" InterestTerms as a vector",
	},
    {
        	" Local_QGM1FModel_Create",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_QGM1FModel_Create",
			" Zero Curve, Model Param1, Model Param2, Model Param3",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an QGM 1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility, MeanReversion or Skew type)",
			" Object (Volatility, MeanReversion or Skew type)",
			" Object (Volatility, MeanReversion or Skew type)",
    },
    {
        	" Local_PXL_QGM1FModel_Create",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_QGM1FModel_Create",
			" Zero Curve, Model Param1, Model Param2, Model Param3",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an QGM 1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Object (Volatility, MeanReversion or Skew type)",
			" Object (Volatility, MeanReversion or Skew type)",
			" Object (Volatility, MeanReversion or Skew type)",
    },
    {
        	" Local_QGM1F_Test",		/// name of the C++ function
            " RRRRRRRRRRRR",			/// 12 parametres = 11 d'entree + 1 parametre de retour 
            " ARM_GP_QGM1F_Test",
			" QGM1F",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " QGM1F Test",
			" Object (QGM1F)",
			" Time",
			" Maturity or Expiry Time",
			" Numeraire or Pay time",
			" X(Time)",
			" Start time",
			" End time",
			" Strike",
			" CapOrPay(1)/FloorOrRec(-1)",
			" Fixed pay times",
			" Fixed periods",
    },
    {
        	" Local_QModel1F_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_Q1FModel_Create",
			" Zero Curve, Model Params,HW Flag",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Q1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Model Params Vector",
			" Flag to degenerate in HW",
    },
    {
        	" Local_PXL_QModel1F_Create",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Q1FModel_Create",
			" Zero Curve, Model Params,HW Flag",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Q1F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Model Params Vector",
			" Flag to degenerate in HW",
    },
    {
        	" Local_QModel1FAna_Create",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_Q1FModelAna_Create",
			" Zero Curve, Model Param1, Model Param2",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Q1F Analytical Model",
			" Object (Zc Curve)",
			" Object (Volatility)",
			" Object (QParameter)",
    },
    {
        	" Local_PXL_QModel1FAna_Create",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Q1FModelAna_Create",
			" Zero Curve, Model Param1, Model Param2",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Q1F Analytical Model",
			" Object (Zc Curve)",
			" Object (Volatility)",
			" Object (QParameter)",
    },
    {
        	" Local_CstModelParam_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_CstModelParam_Create",
			" Param Type, Value,[keep Times]",
            " 1",							/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Creates a Constant Model Parameter",
			" String for param type : Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread,Correlation, Beta, Shift,Q...",
			" Constant value",
			" boolean (default is false)",
    },
    {
        	" Local_PXL_CstModelParam_Create", /// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CstModelParam_Create",
			" Param Type, Value,[keep Times]",
            " 0",							/// not visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
    },
	{
        	" Local_TrigoCorrelParam_Create",	/// name of the C++ function
            " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_TrigoCorrelParam_Create",
			" AsOfDate, DateStrip, Theta,[lowerBound],[upperBound],[interpolatorName]",
            " 1",							/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Create a Trigo Correl Param cos(theta*PI*(i-j)/(n-1))",
			" Date",
			" DateStrip Id/Vector Id",
			" Theta value",
			" lowerBound (size 1)",
			" upperBound (size 1)",
			" interpolatorName default = STEPUPRIGHT",
    },
    {
        	" Local_PXL_TrigoCorrelParam_Create",	/// name of the C++ function
            " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_TrigoCorrelParam_Create",
			" AsOfDate, DateStrip, Theta,[lowerBound],[upperBound],[interpolatorName]",
            " 0",							/// visible in excel
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Create a Trigo Correl Param cos(theta*PI*(i-j)/(n-1))",
			" Date",
			" DateStrip Id",
			" Theta value",
			" lowerBound (size 1)",
			" upperBound (size 1)",
			" interpolatorName default = STEPUPRIGHT",
    },
	{
        	" Local_CstManager_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_CstManager_Create",
			" names,values",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a cst manager object",
			" vector of string",
			" vector of values",
    },
    {
        	" Local_PXL_CstManager_Create",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CstManager_Create",
			" names,values",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a cst manager object",
			" vector of string",
			" vector of values",
    },
	{
        	" Local_ObjManager_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_ObjManager_Create",
			" names,Objects",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates an objects manager",
			" vector of string",
			" vector of Object Ids",
    },
    {
        	" Local_PXL_ObjManager_Create",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ObjManager_Create",
			" names,Objects",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates an objects manager",
			" vector of string",
			" vector of Object Ids",
    },
	{
        	" Local_FlatSurface_Create",	/// name of the C++ function
            " RR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_FlatSurface_Create",
			" value",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a flat surface",
			" value",
    },
    {
        	" Local_PXL_FlatSurface_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_FlatSurface_Create",
			" value",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a flat surface",
			" value",
    },
	{
        	" Local_LinSurface_Create",	/// name of the C++ function
            " RRRRR",					/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_LinSurface_Create",
			" X1,X2,X3,interpMethod",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a linear Surface",
			" array of X1",
			" array of X2",
			" matrix of X3",
            " string: default( LINEAR_COLUMN),LINEAR_COLUMN_ROW,LINEAR_ROW,LINEAR_ROW_COLUMN,STEPUP_RIGHT_COLUMN,STEPUP_RIGHT_ROW,STEPUP_LEFT_COLUMN,STEPUP_LEFT_ROW,CONSTANT",
    },
    {
        	" Local_PXL_LinSurface_Create",	/// name of the C++ function
            " RRRRR",					/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_LinSurface_Create",
			" X1,X2,X3,interpMethod",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a linear Surface",
    },
	{
        	" Local_Surface_Interpolate",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_Surface_Interpolate",
			" SurfaceId,X1,X2",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " interpolate a surface",
			" surface id",
			" x1 value",
			" x2 value",
    },
    {
        	" Local_Surface_Insert",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_Surface_Insert",
			" SurfaceId,X1,X2,X3",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Inserts a point in a surface",
			" surface id",
			" x1 value",
			" x2 value",
			" x3 value",
    },
	{
        	" Local_PXL_Surface_Insert",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Surface_Insert",
			" SurfaceId,X1,X2,X3",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Inserts a point in a surface",
			" surface id",
			" x1 value",
			" x2 value",
			" x3 value",
    },

	{
        	" Local_FromVolSummitToSurfce_Create",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_FromVolSummitToSurface_Create",
			" vol id, interpMethod",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a surface from Summit volatility",
			" string: volatility Id",
			" string: default( LINEAR_COLUMN),LINEAR_COLUMN_ROW,LINEAR_ROW,LINEAR_ROW_COLUMN,STEPUP_RIGHT_COLUMN,STEPUP_RIGHT_ROW,STEPUP_LEFT_COLUMN,STEPUP_LEFT_ROW,CONSTANT",
    },
	{
        	" Local_PXL_FromVolSummitToSurfce_Create",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_FromVolSummitToSurface_Create",
			" vol id, interpMethod",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
    },
	{
        	" Local_FromVolSummitToCurve_Create",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_FromVolSummitToCurve_Create",
			" vol id, interpMethod",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a surface from Summit volatility",
			" string: volatility Id",
			" string: default (LIN), STEPUP_LEFT, STEPUP_RIGHT",
    },
	{
        	" Local_PXL_FromVolSummitToCurve_Create",	/// name of the C++ function
            " RRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_FromVolSummitToCurve_Create",
			" vol id, interpMethod",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
    },


	{
        	" Local_SurfaceModelParam_Create",		/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_SurfaceModelParam_Create",
			" ParamType,SurfaceId,[ParamName],[LowerBound],[UpperBound],[adviseBreakPoint]",
            " 1",									/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a surface model param",
			" param type",
			" surface id",
			" param name",
			" Lower Bound",
			" Upper Bound",
			" boolean (default=false)"
    },
	{
        	" Local_PXL_SurfaceModelParam_Create",	/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SurfaceModelParam_Create",
			" ParamType,SurfaceId,[ParamName],[LowerBound],[UpperBound],[adviseBreakPoint]",
            " 0",									/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a surface model param",
			" param type",
			" surface id",
			" param name",
			" Lower Bound",
			" Upper Bound",
			" boolean (default=false)"
    },
	{
        	" Local_SurfaceListModelParam_Create",		/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_SurfaceListModelParam_Create",
			" ParamType,IndexRange, SurfaceListIds,[ParamName]",
            " 1",									/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a List of surface model param",
			" param type",
			" Range of Index",
			" Range of surface ids",
			" param name",
    },
	{
        	" Local_PXL_SurfaceListModelParam_Create",	/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SurfaceListModelParam_Create",
			" ParamType,IndexRange, SurfaceListIds,[ParamName]",
            " 0",									/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a surface model param",
			" param type",
			" Range of Index",
			" Range of surface ids",
			" param name",
    },
    {
        	" Local_Heston_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_HestonModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Heston model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_Heston_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HestonModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Heston model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_ShiftedHeston_Model_Create",		/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_ShiftedHestonModel_Create",
			" CurveID,ModelParamsId,[MC],[nbSteps],[nbSimulations],[nbIntegrationSteps]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Heston model",
			" curve id",
			" ParamIds: InitialVol,LongTermVol,VolMeanReversion,VolOfVol,Correlation,Shift",
			" MC Method",
			" nbSteps",
			" nbSimulations",
			" nbIntegrationSteps",
    },
	{
        	" Local_PXL_ShiftedHeston_Model_Create",	/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ShiftedHestonModel_Create",
			" CurveID,ModelParamsId,[MC],[nbSteps],[nbSimulations],[nbIntegrationSteps]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Heston model",
			" curve id",
			" ParamIds: InitialVol,LongTermVol,VolMeanReversion,VolOfVol,Correlation,Shift",
			" MC Method",
			" nbSteps",
			" nbSimulations",
			" nbIntegrationSteps",
    },
	{
        	" Local_MSV1FModel_Create",		/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_MSV1FModel_Create",
			" CurveID,ModelParamsId,ForwardTerm,[IRIndex],[IsSwapRate]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Markov Stochastic 1F Model",
			" curve id",
			" ParamIds: MeanReversion,InitialVol,VolMeanReversion,VolOfVol,Shift,LongTermVol",
			" ForwardTerm (in Years)",
			" IRIndex",
			" IsSwapRate",
    },
	{
        	" Local_PXL_MSV1FModel_Create",	/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MSV1FModel_Create",
			" CurveID,ModelParamsId,ForwardTerm,[IRIndex],[IsSwapRate]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Markov Stochastic 1F model",
			" curve id",
			" ParamIds: MeanReversion,InitialVol,VolMeanReversion,VolOfVol,Shift,LongTermVol",
			" ForwardTerm (in Years)",
			" IRIndex",
			" IsSwapRate",
    },
	{
        	" Local_FRMSVModel_Create",		/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_FRMSVModel_Create",
			" CurveID,ModelParamsId,ModelParamsId2,IRIndex",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a FRM with Stochastic Volatility model",
			" curve id",
			" ParamIds: MeanReversion,Volatility,Shift",
			" ParamIds: InitialVol,VolMeanReversion,VolOfVol,LongTermVol",
			" IRIndex",
    },
	{
        	" Local_PXL_FRMSVModel_Create",	/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_FRMSVModel_Create",
			" CurveID,ModelParamsId,ModelParamsId2,IRIndex",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a FRM with Stochastic Volatility model",
			" curve id",
			" ParamIds: MeanReversion,Volatility,Shift",
			" ParamIds: InitialVol,VolMeanReversion,VolOfVol,LongTermVol",
			" IRIndex",
    },
	{
        	" Local_HWSV1FModel_Create",		/// name of the C++ function
            " RRRRRRRRRRR",						/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_HWSV1FModel_Create",
			" CurveId,ModelParamsId,[SolverType],[SolverDatas],[FormulaType],[FormulaDatas],[MaxDecay],[FormulaTypeSO],[FormulaDatasSO],[MaxDecaySO]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull & White 1F model with stochastic volatility",
			" Zc curve",
			" Vol,MRS,VoV,Rho,VolMRS",
			" RK4Constant or RK5Adaptative",
			" NbPY,Min,Max,RefK,CK,CT,CN or Precision,TinyLevel",
			" Heston or Lewis for CF/OSW",
			" 1stLimStp,2ndLimStp,NextLimStp,1stNb,2ndNb,NextNb + (Heston case) 1stLimNb,2ndLimNb,NextLimNb,",
			" MaxDecay (default=0)",
			" Heston or Lewis for SO",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb  or (Lewis + SO) Inv FT pt",
			" MaxDecay for SO (default=0)",
    },
	{
        	" Local_PXL_HWSV1FModel_Create",	/// name of the C++ function
            " RRRRRRRRRRR",						/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HWSV1FModel_Create",
			" CurveId,ModelParamsId,[SolverType],[SolverDatas],[FormulaType],[FormulaDatas],[MaxDecay],[FormulaTypeSO],[FormulaDatasSO],[MaxDecaySO]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull & White 1F model with stochastic volatility",
			" Zc curve",
			" Vol,MRS,VoV,Rho,VolMRS",
			" RK4Constant or RK5Adaptative",
			" NbPY,Min,Max,RefK,CK,CT,CN or Precision,TinyLevel",
			" Heston or Lewis for CF/OSW",
			" 1stLimStp,2ndLimStp,NextLimStp,1stNb,2ndNb,NextNb + (Heston case) 1stLimNb,2ndLimNb,NextLimNb,",
			" MaxDecay (default=0)",
			" Heston or Lewis for SO",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb  or (Lewis + SO) Inv FT pt",
			" MaxDecay for SO (default=0)",
    },
	{
        	" Local_HWSV2FModel_Create",		/// name of the C++ function
            " RRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GP_HWSV2FModel_Create",
			" CurveId,ModelParamsId,[SolverDatas],[FormulaType],[FormulaDatas],[MaxDecay],[FormulaTypeSO],[FormulaDatasSO],[MaxDecaySO]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull & White 2F model with stochastic volatility",
			" Zc curve",
			" Vol,MRS,VolRatio,Rhos,VoV,VolMRS",
			" Precision,TinyLevel",
			" Heston or Lewis for CF/OSW",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb",
			" MaxDecay (default=0)",
			" Heston or Lewis for SO",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb or (Lewis + SO) Inv FT pt",
			" MaxDecay for SO (default=0)",
    },
	{
        	" Local_PXL_HWSV2FModel_Create",	/// name of the C++ function
            " RRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HWSV2FModel_Create",
			" CurveId,ModelParamsId,[SolverDatas],[FormulaType],[FormulaDatas],[MaxDecay],[FormulaTypeSO],[FormulaDatasSO],[MaxDecaySO]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Hull & White 2F model with stochastic volatility",
			" Zc curve",
			" Vol,MRS,VolRatio,Rhos,VoV,VolMRS",
			" Precision,TinyLevel",
			" Heston or Lewis for CF/OSW",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb",
 			" MaxDecay (default=0)",
			" Heston or Lewis for SO",
			" 1stLStp,2ndLStp,NxtLStp,1stNb,2ndNb,NxtNb,IntegPrec + (Heston case) 1stLNb,2ndLNb,NxtLNb or (Lewis + SO) Inv FT pt",
			" MaxDecay for SO (default=0)",
   },

	{
        	" Local_EQHWSV_ModelParamsCreate",		/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_EQHWSV_ModelParamsCreate",
			" ModelParamsId",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Generate the parameter set of a HW 1F model with stochastic volatility adapted to an equity model",
			" CompoundVol, MeanReversion, ScalingVol, VolOfVol, VolMeanReversion",
    },
	{
        	" Local_PXL_EQHWSV_ModelParamsCreate",	/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_EQHWSV_ModelParamsCreate",
			" ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Generate the parameter set of a HW 1F model with stochastic volatility adapted to an equity model",
			" CompoundVol, MeanReversion, ScalingVol, VolOfVol, VolMeanReversion",
    },
	{
        	" Local_EQHWSV_NumMethodsCreate",		/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_EQHWSV_NumMethodsCreate",
			" [IntStep],[ImAxis],[MaxDecay]",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Generate the parameter set of Numerical valorization",
			" IntStep	= GaussLegendre space discterization ( default value = 20) ",
			" ImAxis	= axis of complex integration u= z+I.ImAxis ( default value = 0.5)",
			" MaxDecay  = critera allows the converstion of piecewize function of the time dependent parameters ( default value = 0)",
    },
	{
        	" Local_PXL_EQHWSV_NumMethodsCreate",	/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_EQHWSV_NumMethodsCreate",
			" [IntStep],[ImAxis],[MaxDecay]",
            " 0",									/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Generate the parameter set of Numerical valorization",
			" IntStep	= GaussLegendre space discterization ( default value = 20) ",
			" ImAxis	= axis of complex integration u= z+I.ImAxis ( default value = 0.5)",
			" MaxDecay  = critera allows the converstion of piecewize function of the time dependent parameters ( default value = 0)",
    },
	{
        	" Local_EQHWSV_Create",			
            " RRRRR",							
            " ARM_GP_EQHWSV_Create",
			" ZeroCurve,ModelParamsId,NumMethodsId,[dilatation]",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create Equity Stochastic Volatility Model",
			" Zero Curve",
			" Model Params",
			" Numeric Methods",
			" Dilatation Factor",
    },
	{
        	" Local_PXL_EQHWSV_Create",			
            " RRRRR",								 
            " PXL_ARM_GP_EQHWSV_Create",
			" ZeroCurve,ModelParamsId,NumMethodsId,[dilatation]"
            " 0",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create Equity Stochastic Volatility Model",
			" Zero Curve",
			" Model Params",
			" Numeric Methods",
			" Dilatation Factor",
    },
	{
        	" Local_AMCAndersen_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_AMCAndersen_Create",
			" Frontier Steps, [sorted flag]",
            " 1",								/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates an Andersen exercice frontier calculator",
			" Integer",
			" Boolean",
    },
	{
        	" Local_PXL_AMCAndersen_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_AMCAndersen_Create",
			" Frontier Steps, [sorted flag]",
            " 0",								/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates an Andersen exercice frontier calculator",
			" Integer",
			" Boolean",
    },
	{
        	" Local_AMCLongstaffSchwartz_Create",		/// name of the C++ function
            " RRRRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_AMCLongstaffSchwartz_Create",
			" Frontier Steps, [RegMode], [Span], [IsAutomatic], [Degree]",
            " 1",								/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a LongstaffSchwartz exercice frontier calculator"
			" Integer",
			" Mode (LS/LOESS)",
			" Span (0.8)",
			" IsAutomatic (N)",
			" Degree (N)"
    },
	{
        	" Local_PXL_AMCLongstaffSchwartz_Create",		/// name of the C++ function
            " RRRRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_AMCLongstaffSchwartz_Create",
			" Frontier Steps, [RegMode], [Span], [IsAutomatic], [Degree]",
            " 0",								/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates an LongstaffSchwartz exercice frontier calculator",
			" Integer",
			" Mode (LS/LOESS)",
			" Span (0.8)",
			" IsAutomatic (N)",
			" Degree (N)"
    },
	{
        	" Local_BS_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_BSModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an BS model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_BS_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BSModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an BS model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_CEV_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_CEVModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an CEV model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_CEV_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CEVModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an CEV model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_Merton_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MertonModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Merton model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_Merton_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MertonModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Merton model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_Normal_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_NormalModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Normal model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_Normal_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_NormalModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an Normal model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_SABR_Model_Create",		/// name of the C++ function
            " RRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_SABRModel_Create",
			" Zc CurveID,ModelParamsId,ImpliedVolType",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an SABR model",
			" Zc curve id",
			" model param ids (InitialVol,Beta,Correlation and VolOfVol)",
			" SABR implied vol (default(DIRECTEXACT),ANALYTIC,NINTEGRATION,AUTOMATIC,DIRECTEXACTSTRIKE,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2,SABR_IMPLNVOL,SABR_IMPLNVOL2,SABR_A,SABR_G)",
    },
	{
        	" Local_PXL_SABR_Model_Create",	/// name of the C++ function
            " RRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SABRModel_Create",
			" Zc CurveID,ModelParamsId,ImpliedVolType",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an SABR model",
			" Zc curve id",
			" model param ids (InitialVol,Beta,Correlation and VolOfVol)",
			" SABR implied vol (default(DIRECTEXACT),ANALYTIC,NINTEGRATION,AUTOMATIC,DIRECTEXACTSTRIKE,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2,SABR_IMPLNVOL,SABR_IMPLNVOL2,SABR_A,SABR_G)",

    },
	{
        	" Local_SLN_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_SLNModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an SLN model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_SLN_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SLNModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an SLN model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_GetVarianceSqueeze",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GetVarianceSqueezeFromLocalModel",
			" model id, detail",
            " 1",								/// not visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " Get variance squeeze info from a given local model",
			" model id",
			" detail (default(NO))",
    },
    {
        	" Local_TreeND_Create",	        /// name of the C++ function
            " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_TreeND_Create",
			" Dim Nb, [SchedulerType], [SchedulerData], [SamplerType], [SamplerData], [TruncatorType], [TruncatorData], [ProbasFlag], [ReconnectorType], [SmootherType]",
            " 1",							/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree ND method object",
			" Number of dimensions",
			" Scheduler type, default=ConstantVarianceMeanRevertingScheduler",
			" Scheduler Datas",
			" Sampler type, default=MeanRevertingSampler",
			" Sampler Datas [0]=MinStdDev",
			" Truncator type, default=StdDevTruncator",
			" Truncator datas",
			" State probas : 0 (no computation=default), 1=computation",
			" Reconnector type, default=MeanReconnector",
			" Smoother type, default=DoNothingSmoother",
    },
    {
        	" Local_PXL_TreeND_Create",	    /// name of the C++ function
            " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_TreeND_Create",
			" Dim Nb, [SchedulerType], [SchedulerData], [SamplerType], [SamplerData], [TruncatorType], [TruncatorData], [ProbasFlag], [ReconnectorType], [SmootherType]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates a Tree ND method object",
			" Number of dimensions",
			" Scheduler type, default=ConstantVarianceMeanRevertingScheduler",
			" Scheduler Datas",
			" Sampler type, default=MeanRevertingSampler",
			" Sampler Datas [0]=MinStdDev",
			" Truncator type, default=StdDevTruncator",
			" Truncator datas",
			" State probas : 0 (no computation=default), 1=computation",
			" Reconnector type, default=MeanReconnector",
			" Smoother type, default=DoNothingSmoother",
    },
    {
        	" Local_HW1FModelParam_Create",	        /// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_HW1FModelParam_Create",
			" ModelParams Vector",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a HW1F Model param object",
			" Vector of model params object (vol+mean reversion)",
    },
    {
        	" Local_PXL_HW1FModelParam_Create",	        /// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HW1FModelParam_Create",
			" ModelParams Vector",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a HW1F Model param object",
			" Vector of model params object (vol+mean reversion)",
    },
	{
        	" Local_QNFModelParam_Create",	        /// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_QNFModelParam_Create",
			" Q Parameter,HW1F Model Params Vector, Correlation matrix",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a QNF Model param object",
			" Q Model Param object",
			" Vector of HW1F ModelParam",
			" Correlation matrix",
    },
    {
        	" Local_PXL_QNFModelParam_Create",	        /// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_QNFModelParam_Create",
			" Q Parameter,HW1F Model Params Vector, Correlation matrix",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a QNF Model  object",
			" Q Model Param object",
			" Vector of HW1F ModelParam",
			" Correlation matrix",
    },
	{
        	" Local_QNFModel_Create",				/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_QNFModel_Create",
			" Curve Id,QNF Model Param ID,HW Flag",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a QNF Model object",
			" Curve object id",
			" QNF Model Param Id",
			" Flag to degenerate in HW",
    },
    {
        	" Local_PXL_QNFModel_Create",	        /// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_QNFModel_Create",
			" Curve Id,QNF Model Param ID,HW Flag",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a QNF Model object",
			" Curve object id",
			" QNF Model Param Id",
			" Flag to degenerate in HW",
    },
	{
        	" Local_GPVector_Create",				/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Vector_Create",
			" Values, type",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a GP Vector object",
			" range: vector values",
			" string: type (default(Real), String)",
    },
	{
        	" Local_PXL_GPVector_Create",				/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Vector_Create",
			" Values, type",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
    },
	{
        	" Local_GPMatrix_Create",				/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GPMatrix_Create",
			" Values",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a GP Matrix object",
			" matrix values",
    },
	{
        	" Local_PXL_GPMatrix_Create",				/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GPMatrix_Create",
			" Values",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a GP Matrix object",
			" matrix values",
    },
	{
        	" Local_GPMatrix_Create",				/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_Matrix_Create",
			" Values",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a GP Matrix object",
			" matrix values",
    },
	{
        	" Local_PXL_GPMatrix_Create",				/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Matrix_Create",
			" Values",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a GP Matrix object",
			" matrix values",
    },
	{
        	" Local_GP_LeastSquareRegression",		/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_LeastSquareRegression",
			" X,Y",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calculate a regression",
			" matrix values",
			" vector values"
    },
	{
        	" Local_GP_ACP",						/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_ACP",
			" Matrix",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calculate an ACP",
			" matrix values"
    },
	{
        	" Local_GP_Regression",					/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_Regression",
			" Y,X,XInter,RegMode,[Span]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calculate a Regression",
			" vector values",
			" matrix values",
			" matrix values",
			" string value",
			" double value"
    },
	{
        	" Local_Q1FModel_FX_Create",				/// name of the C++ function
            " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_Q1FFXModel_Create",
			" Curve Id,Q1F Model Param IDs, Spot, For Curve",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Q1F FX Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
			" Foreign Curve Id",
    },
    {
        	" Local_PXL_Q1FModel_FX_Create",	        /// name of the C++ function
            " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Q1FFXModel_Create",
			" Curve Id,Q1F Model Param IDs, Spot, For Curve",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Q1F FX Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
			" Foreign Curve Id",
    },
	{
        	" Local_CEVModel_FX_Create",				/// name of the C++ function
            " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_CEVFXModel_Create",
			" Curve Id,CEV Model Param IDs, Spot, For Curve",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a CEV FX Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
			" Foreign Curve Id",
    },
    {
        	" Local_PXL_CEVModel_FX_Create",	        /// name of the C++ function
            " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CEVFXModel_Create",
			" Curve Id,CEV Model Param IDs, Spot, For Curve",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Q1F FX Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
			" Foreign Curve Id",
    },
	{
        	" Local_BSModel_FX_Create",				/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_BSFXModel_Create",
			" Curve Id,BS Model Param IDs, Spot, For Curve",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a BS FX Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
			" Foreign Curve Id",
    },
    {
        	" Local_PXL_BSModel_FX_Create",	        /// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BSFXModel_Create",
			" Curve Id,BS Model Param IDs, Spot, For Curve",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
    },
	{
        	" Local_HestonModel_FX_Create",				/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_HestonFXModel_Create",
			" Domestic curve Id, ModelParam IDs, SpotValue, Foreign curve Id, [MCScheme]",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Heston FX Model object",
			" string: curve object id",
			" vector of string: Vector of modelParam Ids",
			" value: spot value (1 forccy = X dom Ccy)",
			" string: foreign curev Id",
			" MCScheme: (EULER/ANDREASEN)"
    },
    {
        	" Local_PXL_HestonModel_FX_Create",	        /// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HestonFXModel_Create",
			" Domestic curve Id,ModelParam IDs, Spot, Foreign curve Id, [MCScheme]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
			" Creates a Heston FX Model object",
			" string: curve object id",
			" vector of string: Vector of modelParam Ids",
			" value: spot value (1 forccy = X dom Ccy)",
			" string: foreign curev Id",
			" MCScheme: (EULER/ANDREASEN)"
    },
	{
        	" Local_SABRModel_FX_Create",				/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_SABRFXModel_Create",
			" Domestic curve Id, ModelParam IDs, SpotValue, Foreign curve Id",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SABR FX Model object",
			" string: curve object id",
			" vector of string: Vector of modelParam Ids",
			" value: spot value (1 forccy = X dom Ccy)",
			" string: foreign curev Id",
    },
    {
        	" Local_PXL_SABRModel_FX_Create",	        /// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SABRFXModel_Create",
			" Domestic curve Id,ModelParam IDs, Spot, Foreign curve Id",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
    },
	{
        	" Local_MixtureModel_FX_Create",				/// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_MixtureFXModel_Create",
			" Domestic curve Id, ModelParam IDs, SpotValue, Foreign curve Id",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Mixture FX Model object",
			" string: curve object id",
			" vector of string: Vector of modelParam Ids",
			" value: spot value (1 forccy = X dom Ccy)",
			" string: foreign curve Id",
    },
    {
        	" Local_PXL_MixtureModel_FX_Create",	        /// name of the C++ function
            " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MixtureFXModel_Create",
			" Domestic curve Id,ModelParam IDs, Spot, Foreign curve Id",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
    },
	{
        	" Local_MixtureFXModel_Calibrate",				/// name of the C++ function
            " RRRRRRRRR",									/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_MixtureFXModel_Calibrate",
			" Fwd, Expiry, CallPut, Strikes, Vols, DecVol, Alpha, Lambda",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calibrate a Mixture FX Model",
			" value: Fwd FX value",
			" value: Option Expiry",
			" string: (C/F)",
			" vector of value: vector of strikes",
			" vector of value: vector of volatilities",
			" vector of value: DecVol Init/LB/UB",
			" vector of value: Alpha Init/LB/UB",
			" vector of value: Lambda Init/LB/UB",

    },
	{
        	" Local_HestonModel_Eq_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_HestonEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Heston Eq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
    {
        	" Local_PXL_HestonModel_Eq_Create",	        /// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HestonEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Heston Eq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
	{
        	" Local_Q1FModel_Eq_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_Q1FEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Q1F Eq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
    {
        	" Local_PXL_Q1FModel_Eq_Create",	        /// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Q1FEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Q1F Eq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
	{
        	" Local_BSModel_Eq_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_BSEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 1",						/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a BS Eq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
    {
        	" Local_PXL_BSModel_Eq_Create",	        /// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BSEqModel_Create",
			" Curve Id,BS Model Param IDs, Spot",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a BSEq Model object",
			" Curve object id",
			" Vector of Model Param Ids",
			" Spot value",
    },
	{
        	" Local_ModelNameMap_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_ModelNameMap_Create",
			" Names, Models Objs, Other Model Names",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a model name map object",
			" names",
			" model objects",
			" other model names",
    },
	{
        	" Local_PXL_ModelNameMap_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ModelNameMap_Create",
			" Names, Models Objs, Other Model Names",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a model name map object",
			" names",
			" model objects",
			" other model names",
    },
	{
        	" Local_ModelNameMap_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_ModelNameMap_Create",
			" Names, Models Objs, Other Model Names",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a model name map object",
			" names",
			" model objects",
			" other model names",
    },
	{
        	" Local_PXL_ModelNameMap_Create",				/// name of the C++ function
            " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ModelNameMap_Create",
			" Names, Models Objs, Other Model Names",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a model name map object",
			" names",
			" model objects",
			" other model names",
    },
	{
        	" Local_MultiAssetsModel_Create",			/// name of the C++ function
            " RRRR",										/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_MultiAssetsModel_Create",
			" Model Name Map, correlation Matrix, [MultiAssetsName]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates the multi asset model whose name is MultiAssets Name",
			" string: model map id",
			" string: correlation matrix id",
			" string: MultiAsset model name (Default:unknown) "
    },
	{
        	" Local_PXL_MultiAssetsModel_Create",			/// name of the C++ function
            " RRRR",										/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MultiAssetsModel_Create",
			" Model Name Map, correlation Matrix, MultiAssets Name",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a multi asset model",
			" model name map",
			" correlation matrix",
			" MultiAsset model name "
    },
	{
        	" Local_PricingModel_GetModelParamId",			/// name of the C++ function
            " RRRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Model_GetModelParamId",
			" Model, ParamType, [FactorNb]",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get model param data type from a model",
			" model object",
			" string for param type : Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread,Correlation",
			" Factor Nb"
    },
	{
        	" Local_PXL_PricingModel_GetModelParamId",			/// name of the C++ function
            " RRRR",											/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Model_GetModelParamId",
			" Model, ParamType, [FactorNb]",
            " 0",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get model param data type from a model",
			" model object",
			" string for param type : Volatility,MeanReversion,VolatilityRatio,MeanReversionSpread,Correlation"
			" Factor Nb"
    },
	{
        	" Local_PricingModel_SetModelMap",			/// name of the C++ function
            " RRRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Model_SetModelMap",
			" Model,ModelMap",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " set modelMap to model",
			" model object",
			" model Name object",
			" Factor Nb"
    },
	{
        	" Local_PXL_PricingModel_SetModelMap",			/// name of the C++ function
            " RRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Model_SetModelMap",
			" Model,ModelMap ",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " set modelMap to model",
			" model object",
			" model Name object"
    },
	{
        	" Local_PricingModel_GetModelMap",			/// name of the C++ function
            " RR",										/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_Model_GetModelMap",
			" Model",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get modelMap from a model",
			" model object"
    },
	{
        	" Local_PXL_PricingModel_GetModelMap",			/// name of the C++ function
            " RR",										/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Model_GetModelMap",
			" Model",
            " 0",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get modelMap from a model",
			" model object"
    },
	{
        	" Local_GetModelFromModelMap",			/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GetModelFromModelMap",
			" ModelMap,Name",
            " 1",							/// visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get model from model map",
			" model map object",
			" name"
    },
	{
        	" Local_PXL_GetModelFromModelMap",			/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GetModelFromModelMap",
			" ModelMap,Name",
            " 0",							/// not visible in excel
            XLLOCALARM_GETFROMID_GROUP,
            " ",
            " ",
            " get model from model map",
			" model map object",
			" name"
    },

	{
        	" Local_SetModelToModelMap",			/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_SetModelToModelMap",
			" ModelMap,Name,ModelId",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " get model from model map",
			" model map object",
			" name",
			" model object",
    },
	{
        	" Local_PXL_SetModelToModelMap",			/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SetModelToModelMap",
			" ModelMap,Name,ModelId",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " get model from model map",
			" model map object",
			" name",
			" model object",
    },
	{
        	" Local_Create2IRFXModel",			/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_Create2IRFXModel_Create",
			" names,models,correlation",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 2IR+FX Model",
			" names",
			" models",
			" correlation",
    },
	{
        	" Local_PXL_Create2IRFXModel",			/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Create2IRFXModel_Create",
			" names,models,correlation",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 2IR+FX Model",
			" names",
			" models",
			" correlation",
    },
	{
        	" Local_Create1IRFXModel",			/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_Create1IRFXModel_Create",
			" names,models,correlation,[Model2IRFX]",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 1IR+FX Model",
			" names",
			" models",
			" correlation",
			" Model2IRFX inputed for the correlation (NULL)"
    },
	{
        	" Local_PXL_Create1IRFXModel",			/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Create1IRFXModel_Create",
			" names,models,correlation,[Model2IRFX]",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 1IR+FX Model",
			" names",
			" models",
			" correlation",
			" Model2IRFX inputed for the correlation (NULL)"
    },
	{
        	" Local_NP1IRNFXModel_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_NP1IRNFXModel_Create",
			" names,models,correlation",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 1IR+FX Model",
			" names",
			" models",
			" correlation"
    },
	{
        	" Local_PXL_NP1IRNFXModel_Create",			/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_NP1IRNFXModel_Create",
			" names,models,correlation",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid model 1IR+FX Model",
			" names",
			" models",
			" correlation"
    },
	{
        	" Local_2IRFXSV_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_2IRFXSV_Create",
			" names,models,correlation",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid 2IR+FX SV Model",
			" names",
			" models",
			" correlation"
    },
	{
        	" Local_PXL_2IRFXSV_Create",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_2IRFXSV_Create",
			" names,models,correlation",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates an hybrid 2IR+FX SV Model",
			" names",
			" models",
			" correlation"
    },
	{
        	" Local_NP1IRNFX_CalibrateFunctional",		/// name of the C++ function
            " RRRRRRR",									/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_NP1IRNFX_CalibrateFunctional",
			" NP1IRNFXId,ResetDates,Densities,[GridSize],[StdDevNb],[Rescaling]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " NP1IRNFX model",
			" Reset Dates",
			" Density Functors",
			" GridSize (501)",
			" StdDevNb (6)",
			" Rescaling (Y)"
    },
	{
        	" Local_PXL_NP1IRNFX_CalibrateFunctional",		/// name of the C++ function
            " RRRRRRR",										/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_NP1IRNFX_CalibrateFunctional",
			" NP1IRNFXId,ResetDates,Densities,[GridSize],[StdDevNb],[Rescaling]",
            " 0",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " NP1IRNFX model",
			" Reset Dates",
			" Density Functors",
			" GridSize (501)",
			" StdDevNb (6)",
			" Rescaling (Y)"
    },
	{
        	" Local_2IRFXSV_CalibrateFunctional",		/// name of the C++ function
            " RRRRRRR",									/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_2IRFXSV_CalibrateFunctional",
			" 2IRFXSVId,ResetDates,Densities,[GridSize],[StdDevNb]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " 2IRFXSV model",
			" Reset Dates",
			" Density Functors",
			" GridSize (501)",
			" StdDevNb (6)"
    },
	{
        	" Local_PXL_2IRFXSV_CalibrateFunctional",		/// name of the C++ function
            " RRRRRRR",										/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_2IRFXSV_CalibrateFunctional",
			" 2IRFXSVId,ResetDates,Densities,[GridSize],[StdDevNb]",
            " 0",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " 2IRFXSV model",
			" Reset Dates",
			" Density Functors",
			" GridSize (501)",
			" StdDevNb (6)"
    },
	{
        	" Local_HWHWQtoModel_Create",			/// name of the C++ function
            " RRRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_HWHWQtoModel_Create",
			" names,models,correlation,FxFlag",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates the qto model",
			" names",
			" models",
			" correlation",
			" FxFlag",
    },
	{
        	" Local_PXL_HWHWQtoModel_Create",		/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HWHWQtoModel_Create",
			" names,models,correlation",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates the qto model",
			" names",
			" models",
			" correlation",
    },
	{
        	" Local_HWHW2FQtoModel_Create",			/// name of the C++ function
            " RRRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_HWHW2FQtoModel_Create",
			" names,models,correlation,FxFlag",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates the qto model",
			" names",
			" models",
			" correlation",
			" FxFlag",
    },
	{
        	" Local_PXL_HWHW2FQtoModel_Create",		/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HWHW2FQtoModel_Create",
			" names,models,correlation",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates the qto model",
			" names",
			" models",
			" correlation",
    },
	{
        	" Local_CreateFwdMarginModel",		/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_CreateFwdMarginModell_Create",
			" basis curve",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates a fwd margin Model",
			" basis curve",
    },
	{
        	" Local_PXL_CreateFwdMarginModel",		/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CreateFwdMarginModell_Create",
			" basis curve",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " creates a fwd margin Model",
			" basis curve",
    },
	{
        	" Local_SetRefModelNameToMultiAsset",		/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_SetRefModelNameToMultiAsset",
			" name,multi aset model",
            " 1",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " set the reference model in a multi asset model",
			" name",
			" multi asset object",
    },
	{
        	" Local_PXL_SetRefModelNameToMultiAsset",		/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SetRefModelNameToMultiAsset",
			" name,multi aset model",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " set the reference model in a multi asset model",
			" name",
			" multi asset object",
    },
    {
        	" Local_IntegratedCorrelation",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_IntegratedCorrelation",
			" Model,Tenor1,Tenors,Expiries",
            " 1",							/// visible in excel
            XLLOCALARM_GLOBAL_GROUP,
            " ",
            " ",
            " Computes the integrated correlation",
			" Model",
			" Tenor1",
			" Tenors",
			" Expiries",
    },
    {
        	" Local_PXL_IntegratedCorrelation",	/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_IntegratedCorrelation",
			" Model,Tenor1,Tenors,Expiries",
            " 0",							/// visible in excel
            XLLOCALARM_GLOBAL_GROUP,
            " ",
            " ",
            " Computes the integrated correlation",
			" Model",
			" Tenor1",
			" Tenors",
			" Expiries",
    },
	{
        	" Local_LocalNormal_Model_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_LocalNormalModel_Create",
			" CurveID,ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Local Normal model",
			" curve id",
			" model param ids",
    },
	{
        	" Local_PXL_LocalNormal_Model_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_LocalNormalModel_Create",
			" CurveID,ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Local Normal model",
			" model param ids",
    },
	{
        	" Local_LocalSLN_Model_Create",		/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_LocalSLNModel_Create",
			" ModelParamsId",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Local Shifted LogNormal model",
			" Model param ids (Fwd adj & Vol)",
    },
	{
        	" Local_PXL_LocalSLN_Model_Create",	/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_LocalSLNModel_Create",
			" ModelParamsId",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Local Shifted LogNormal model",
			" Model param ids (Fwd adj & Vol)",
    },
	{
        	" Local_Local_Model_Calibrate",		/// name of the C++ function
            " RRRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_LocalModel_Calibrate",
			" MultiAssetsId,LocalModelName,PortfolioId,EvalDates",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calibrates a Local model",
			" MultiAssets model id",
			" Local model name",
			" Portfolio id",
			" Eval dates",
    },
	{
        	" Local_PXL_Local_Model_Calibrate",	/// name of the C++ function
            " RRRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_LocalModel_Calibrate",
			" MultiAssetsId,LocalModelName,PortfolioId,EvalDates",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calibrates a Local model",
			" MultiAssets model id",
			" Local model name",
			" Portfolio id",
			" Eval dates",
    },
	{
        	" Local_LocalModel_CalibrateFunctional",		/// name of the C++ function
            " RRRRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_LocalModel_CalibrateFunctional",
			" MultiAssetsId,LocalModelName,Securities,Densities,[Rescaling]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calibrates a Local model to densities",
			" MultiAssets model id",
			" Local model name",
			" Security Ids",
			" Density Ids",
			" Rescaling (Y/N, default N)",
    },
	{
        	" Local_PXL_LocalModel_CalibrateFunctional",	/// name of the C++ function
            " RRRRRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_LocalModel_CalibrateFunctional",
			" MultiAssetsId,LocalModelName,Securities,Densities,[Rescaling]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Calibrates a Local model to densities",
			" MultiAssets model id",
			" Local model name",
			" Security Ids",
			" Density Ids",
			" Rescaling (Y/N, default N)",
	},
	{
        	" Local_Local_Model_VarSqueeze",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_LocalModel_VarSqueeze",
			" MultiAssetsId,LocalModelName",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Var Squeeze in Local model?",
			" MultiAssets model id",
			" Local model name",
    },

	{	
        	" Local_CFMethod_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_CFMethod_Create",
			" [MethodeType],[IntegralParameters]",
            " 1",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the CF Method",
			" string: Type of Closed Form, analytic or integral (Default:unknown)",
			" string: ID of the MatrixParameters (Default:Matrix null)",

    },
    {
        	" Local_PXL_CFMethod_Create",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CFMethod_Create",
			" ",
            " 0",						/// visible in excel
            XLLOCALARM_GENNUMMETHODS_GROUP,
            " ",
            " ",
            " Creates the CF Method",
			" string: Type of Closed Form, analytic or integral (Default:unknown)",
			" vector: When MethodType=INTEGRAL,[xmin,xmax,NbPoints] for GL approximation (Default:[-6.0,+6.0,121])",
    },
	{
        	" Local_MarketIRModel_Create",		/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MarketIRModel_Create",
			" MarketDatas,[VNS_PricingMethod]",
            " 1",								/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a IR Market Model",
			" Mkt Data Manager + Keys",
			" Var. notio. swaption: MONEYNESS or ATM (default: MONEYNESS)",
    },
	{
        	" Local_PXL_MarketIRModel_Create",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MarketIRModel_Create",
			" MarketDatas,[VNS_PricingMethod]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a IR Market Model",
			" Mkt Data Manager + Keys",
			" Var. notio. swaption: MONEYNESS or ATM (default: MONEYNESS)",

    },
	{
        	" Local_SmiledFRMModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",					/// 12 parametres = 11 d'entree + 1 parametre de retour 
            " ARM_GP_SmiledFRMModel_Create",
			" Zero Curve, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [allowInterpol], [swaptionApprox], [recorrel], [rescalling]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled FRM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" allowInterpol",
			" swaptionApprox (Atm, Moment, Local)",
			" recorrel",
			" rescalling (N)",
    },
	{
        	" Local_PXL_SmiledFRMModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",						/// 12 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SmiledFRMModel_Create",
			" Zero Curve, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [allowInterpol], [swaptionApprox], [recorrel], [rescalling]",
             " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled FRM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" allowInterpol",
			" swaptionApprox (Atm, Moment, Local)",
			" recorrel",
			" rescalling (N)",
    },
	{
        	" Local_ImpSampler_Optimize",			/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_ImpSampler_Optimize",
			" GenSec, Model, InitGuess, LowerBound, UpperBound, [WithMC], [NbSteps], [Bootstrap]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Optimize a model for importance",
			" GenSec",
			" Model",
			" InitGuess",
			" LowerBound",
			" UpperModel",
			" WithMC: Y,N (Y)",
			" NbSteps",
			" Bootstrap: Y,N (N)",
    },
	{
        	" Local_PXL_ImpSampler_Optimize",		/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ImpSampler_Optimize",
			" GenSec, Model, InitGuess, LowerBound, UpperBound, [WithMC], [NbSteps], [Bootstrap]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Optimize a model for importance",
			" GenSec",
			" Model",
			" InitGuess",
			" LowerBound",
			" UpperModel",
			" WithMC: Y,N (Y)",
			" NbSteps",
			" Bootstrap: Y,N (N)",
    },
	{
        	" Local_SmiledMarketModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",						/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_SmiledMarketModel_Create",
			" Calib Pattern, Zero Curve, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [allowInterpol], [swaptionApprox], [recorrel]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled Market Model with default numeraire = Terminal ZC",
			" Calib Pattern (Libor,SwapCol,SwapDiag)",
			" Object (Zc Curve)",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" allowInterpol",
			" swaptionApprox (Atm, Moment, Local)",
			" recorrel",
    },
	{
        	" Local_SmiledMarketModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",							/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SmiledMarketModel_Create",
			" Calib Pattern, Zero Curve, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [allowInterpol], [swaptionApprox], [recorrel]",
             " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled Market Model with default numeraire = Terminal ZC",
			" Calib Pattern (Libor,SwapCol,SwapDiag)",
			" Object (Zc Curve)",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" allowInterpol",
			" swaptionApprox (Atm, Moment, Local)",
			" recorrel",
    },
	{
        	" Local_SmiledMarketModelDS_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",						/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " ARM_GP_SmiledMarketModelDS_Create",
			" CalibPattern, StartDate, EndDate, ResetFreq, IndexFreq, IndexType, ResetTiming, [dayCount], [resetCalendar], [fwdRule], [intRule], [stubRule], [resetGap]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled Market Model DateStrip",
			" Calib Pattern (Libor,CMS,VMS)",
			" StartDate",
			" EndDate",
			" ResetFreq",
			" IndexFreq",
			" IndexType (tenor in months, for CMS only)",
			" ADV(ance), ARR(ears)",
			" dayCount",
			" resetCalendar",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
    },
	{
        	" Local_SmiledMarketModelDS_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRR",							/// 13 parametres = 12 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SmiledMarketModelDS_Create",
			" CalibPattern, StartDate, EndDate, ResetFreq, IndexFreq, IndexType, ResetTiming, [dayCount], [resetCalendar], [fwdRule], [intRule], [stubRule], [resetGap]",
             " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled Market Model DateStrip",
			" Calib Pattern (Libor,CMS,VMS)",
			" StartDate",
			" EndDate",
			" ResetFreq",
			" IndexFreq",
			" IndexType (tenor in months, for CMS only)",
			" ADV(ance), ARR(ears)",
			" dayCount",
			" resetCalendar",
			" F(ollowing), MF (modified Following), P(revious), MP(modified Previous)",
			" ADJ, UNADJ, default ADJ",
			" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd) ",
			" in days",
    },
	{
        	" Local_SVBGMModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_SVBGMModel_Create",
			" Zero Curve, Shift, Alpha, Nu, Rho, BetaCorrel, RateVolCorrel, VolVolCorrel, [recorrel], [factorsNb], [min ratio], [Proxy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Shift",
			" Alpha",
			" Nu",
			" Rho",
			" BetaCorrel",
			" RateVolCorrel",
			" VolVolCorrel",
			" recorrel",
			" factorsNb",
			" minratio",
			" Proxy",
    },
	{
        	" Local_PXL_SVBGMModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SVBGMModel_Create",
			" Zero Curve, Shift, Alpha, Nu, Rho, BetaCorrel, RateVolCorrel, VolVolCorrel, [recorrel], [factorsNb], [min ratio], [Proxy]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Shift",
			" Alpha",
			" Nu",
			" Rho",
			" BetaCorrel",
			" RateVolCorrel",
			" VolVolCorrel",
			" recorrel",
			" factorsNb",
			" minratio",
			" Proxy",
    },
	{
        	" Local_BGMSV1FModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_BGMSV1FModel_Create",
			" Zero Curve, Shift, Level, InitialVar, LongTermVar, VarVol, VarMeanRev, RateVolCorrel, BetaCorrel, LocalCalibration, [recorrel], [factorsNb], [min ratio], [LocalRhoCalib], [StdDevForCalib], [Proxy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Shift",
			" Level",
			" InitialVar",
			" LongTermVar",
			" VarVol",
			" VarMeanRev",
			" RateVolCorrel",
			" BetaCorrel",
			" Local Or Global Calib"
			" recorrel",
			" factorsNb",
			" minratio",
			" global calib",
			" Proxy",
    },
	{
        	" Local_PXL_BGMSV1FModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BGMSV1FModel_Create",
			" Zero Curve, Shift, Level, InitialVar, LongTermVar, VarVol, VarMeanRev, RateVolCorrel, BetaCorrel, LocalCalibration, [recorrel], [factorsNb], [min ratio], [LocalRhoCalib], [StdDevForCalib], [Proxy]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Shift",
			" Level",
			" InitialVar",
			" LongTermVar",
			" VarVol",
			" VarMeanRev",
			" RateVolCorrel",
			" BetaCorrel",
			" recorrel",
			" factorsNb",
			" minratio",
			" global calib",
			" Proxy",
    },
	{
        	" Local_BGMSV2FModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_BGMSV2FModel_Create",
			" Zero Curve, BetaCorrel, v01, kappa1, v02, kappa2, [rho1], [rho2], [LocalRho1Calib], [LocalRho2Calib], [Shift], [recorrel], [factorsNb], [min ratio], [StdDevForCalib], [Proxy]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" "
    },
	{
        	" Local_PXL_BGMSV2FModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_BGMSV2FModel_Create",
			" Zero Curve, BetaCorrel, v01, kappa1, v02, kappa2, [rho1], [rho2], [LocalRho1Calib], [LocalRho2Calib], [Shift], [recorrel], [factorsNb], [min ratio], [StdDevForCalib], [Proxy]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" "
    },
	{
        	" Local_SVMMSpreadModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_SVMMSpreadModel_Create",
			" Zero Curve, Level, InitialVar, LongTermVar, VarVol, VarMeanRev, RateVolCorrel, BetaCorrel, [recorrel], [factorsNb], [min ratio]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Level",
			" InitialVar",
			" LongTermVar",
			" VarVol",
			" VarMeanRev",
			" RateVolCorrel",
			" BetaCorrel",
			" recorrel",
			" factorsNb",
			" minratio",
    },
	{
        	" Local_PXL_SVMMSpreadModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRR", /// 11 parametres = 10 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SVMMSpreadModel_Create",
			" Zero Curve, Level, InitialVar, LongTermVar, VarVol, VarMeanRev, RateVolCorrel, BetaCorrel, [recorrel], [factorsNb], [min ratio]",
            " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a SVBGM Model with default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Level",
			" InitialVar",
			" LongTermVar",
			" VarVol",
			" VarMeanRev",
			" RateVolCorrel",
			" BetaCorrel",
			" recorrel",
			" factorsNb",
			" minratio",
    },
	{
        	" Local_QGM2FModel_Create",		/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_QGM2FModel_Create",
			" Zero Curve, ModelParams 1F, ModelParams 2F",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an QGM 2F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Vector ModelParams (Volatility, MeanReversion, Skew and CrossFactors)",
			" Vector ModelParams (Volatility, MeanReversion, Skew and CrossFactors)",
    },
    {
        	" Local_PXL_QGM2FModel_Create",	/// name of the C++ function
            " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_QGM2FModel_Create",
			" Zero Curve, ModelParams 1F, ModelParams 2F",
            " 0",							/// not visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates an QGM 2F Model with a default numeraire = Terminal ZC",
			" Object (Zc Curve)",
			" Vector ModelParams (Volatility, MeanReversion, Skew and CrossFactors)",
			" Vector ModelParams (Volatility, MeanReversion, Skew and CrossFactors)",
    },
	{
        	" Local_SmiledFXModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRR",					/// 15 parametres = 14 d'entree + 1 parametre de retour 
            " ARM_GP_SmiledFXModel_Create",
			" Dom Zero Curve, Foreign Zero Curve, FXSpot, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [recorrel], [rescalling], [Model2IRFX]",
            " 1",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled FRM Model with default numeraire = Terminal ZC",
			" Object (Dom Zc Curve)",
			" Object (For Zc Curve)",
			" FXSpot",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" recorrel",
			" rescalling (N)",
			" Model2IRFX inputed for the correlation (NULL)"
    },
	{
        	" Local_PXL_SmiledFXModel_Create",		/// name of the C++ function
            " RRRRRRRRRRRRRRR",					/// 15 parametres = 14 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_SmiledFXModel_Create",
			" Dom Zero Curve, Foreign Zero Curve, FXSpot, BetaCorrel, Hump, factorsNb, [timeStepsNb], [gridSize], [stdDevNb], [skipPDE], [switchTheta], [recorrel], [rescalling], [Model2IRFX]",
             " 0",							/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Creates a Smiled FRM Model with default numeraire = Terminal ZC",
			" Object (Dom Zc Curve)",
			" Object (For Zc Curve)",
			" FXSpot",
			" BetaCorrel",
			" Hump",
			" factorsNb",
			" timeStepsNb",
			" gridSize",
			" stdDevNb",
			" skipPDE",
			" switchTheta",
			" recorrel",
			" rescalling (N)",
			" Model2IRFX inputed for the correlation (NULL)"
    },
	{
			" Local_BiSVMM_Create",
			" RRRR",
			" ARM_GP_BiSVMM_Create",
			" Names, Models, Correlation",
			" 1",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
			" Creates a cross svmm model"
	},
	{
			" Local_PXL_BiSVMM_Create",
			" RRRR",
			" PXL_ARM_GP_BiSVMM_Create",
			" Names, Models, Correlation",
			" 0",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
			" Creates a cross svmm model"
	},
	{
			" Local_HWxSVMMSpread_Create",
			" RRRRRR",
			" ARM_GP_HWxSVMMSpread_Create",
			" modelNames, modelIds, [HW2Fmodel], [CorrIdxEndTimes], [cstXCorrel]",
			" 1",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
			" Creates an hybrid hw x svmm spread model",
			" ",
			" ",
			" HW2F needed to compute cross correlation",
			" Index End Times to compute cross correlation",
			" Constant Cross Correlation (if not computed with HW2F) (default = 0)"
	},
	{
			" Local_PXL_HWxSVMMSpread_Create",
			" RRRRRR",
			" PXL_ARM_GP_HWxSVMMSpread_Create",
			" modelNames, modelIds, [HW2Fmodel], [CorrIndexEndTimes], [constantCrossCorrel]",
			" 0",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
			" Creates an hybrid hw x svmm spread model"
	},
	{
			" Local_HWSBGMQtoModel_Create",
			" RRRRR",
			" ARM_GP_HWSBGMQtoModel_Create",
			" modelNames, modelIds, correlationMatrix",
			" 1",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
	},
	{
			" Local_PXL_HWSBGMQtoModel_Create",
			" RRRRR",
			" PXL_ARM_GP_HWSBGMQtoModel_Create",
			" modelNames, modelIds, correlationMatrix",
			" 0",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
	},
	{
			" Local_HWSVBGMQtoModel_Create",
			" RRRRR",
			" ARM_GP_HWSVBGMQtoModel_Create",
			" modelNames, modelIds, correlationMatrix",
			" 1",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
	},
	{
			" Local_PXL_HWSVBGMQtoModel_Create",
			" RRRRR",
			" PXL_ARM_GP_HWSVBGMQtoModel_Create",
			" modelNames, modelIds, correlationMatrix",
			" 0",
			XLLOCALARM_GENMODELS_GROUP,
			" ",
			" ",
	},
	/// END GENERIC PRICER
