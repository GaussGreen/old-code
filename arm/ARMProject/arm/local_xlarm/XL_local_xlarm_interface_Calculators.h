/////////////////////////////////////////
///
/// Copyright (c) CDC IXIS CM July 2003 Paris
///
/// \brief contains all the addins for the calculators
/// \author E. Benhamou, JM Prie, E. Mostafa Ezzine
/// this is part of a general table that contains all
/// the excel addins
///
///	\version 1.0
///	\date December 2003
//////////////////////////////////////////


	/// Calculators Section
    {
        	" Local_DateStripCombiner_Create",	/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DateStripCombiner_Create",
			" dateStrips,[merge func]",
            " 1",						/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Combines datestrip",
			" vector of dateStrip objects",
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD)"
    },

    {
        	" Local_PXL_DateStripCombiner_Create",	/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DateStripCombiner_Create",
			" dateStrips,[merge func]",
            " 0",						/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Combines datestrip",
			" objects",
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD)"
    },
    {
        	" Local_DataStripCombiner_GetData",		/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_DateStripCombiner_GetData",
			" DateStripCombiner,dataType,dateStripNb",
            " 1",						/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Get Data from a Datestrip Combiner",
			" objects",
            " STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE (FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD)"
			" dateStrip Nb",
    },
    {
        	" Local_DataStripCombiner_GetMergeData",	/// name of the C++ function
            " RR",										/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_DateStripCombiner_GetMergeData",
			" DateStripCombiner",
            " 1",						/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Get Data from a Datestrip Combiner",
			" object"
    },
    {
        	" Local_CRFCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_CRFCalculator_Create",
			" Start,End,Strike,P/R,CpnDatas,MktDatas,[FixEnd],[FixDC],[ResetGap],[Lvge],[CpnMin],[CpnMax],[FundSpd],[FundDatas],[Nominal],[ExerGap],[NbNCall],[ExerFee],[Flags],[FundingNominal]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CRF Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq(M,Q,S,A),DC,Timing(ADV), Idxterm(1M,3M,6M,1Y),IndexDC(A360),ResetCal(EUR),PayCal (EUR),StubRule (SS)",
            " Market Datas = Mkt Data Mger Id, Cpn Yc Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key, Fund Yc Key, Basis Yc Key, Forex Key",
            " Fix leg end date (start date)",
            " Fix leg day count (coupon DC)",
            " Coupon reset lag (2 days)",
            " Leverage (1.0)",
            " Coupon min (0.0)",
            " Coupon max (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
            " Notification lag (reset Lag)",
            " Number of non call period (0)",
            " Exercise fees (0.0)",
            " Autocal Flags = Sigma(Y),MRSPfType(Y),cap(Y),floor(Y),product(CRF),modelType(HWM1F),Skew(Y),MRSStrikeType(Kequi)",
            " FundingNominal (Nominal*1/FxSpot)"
    },
    {
        	" Local_PXL_CRFCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_CRFCalculator_Create",
			" Start,End,Strike,P/R,CpnDatas,MktDatas,[FixEnd],[FixDC],[ResetGap],[Lvge],[CpnMin],[CpnMax],[FundSpd],[FundDatas],[Nominal],[ExerGap],[NbNCall],[ExerFee],[Flags],[FundingNominal]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CRF Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq (M,Q,S,A), DC, Timing (ADV), Idx term (1M,3M,6M,1Y), Index DC (A360), Reset Cal (EUR), Pay Cal (EUR)",
            " Market Datas = Mkt Data Mger Id, Cpn Yc Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key, Fund Yc Key, Basis Yc Key, Forex Key",
            " Fix leg end date (start date)",
            " Fix leg day count (coupon DC)",
            " Coupon reset lag (2 days)",
            " Leverage (1.0)",
            " Coupon min (0.0)",
            " Coupon max (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
            " Notification lag (reset Lag)",
            " Number of non call period (0)",
            " Exercise fees (0.0)",
            " Autocal Flags = diag OSW (Y), short term vanilla (Y), cap (Y), floor (N), prod to price (CRF),QGMFlag,SkewAutoCalFlag, SkewShift",
            " FundingNominal (Nominal*1/FxSpot)"
    },
    {
        	" Local_Crude_CRFCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",	    /// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_GC_Crude_CRFCalculator_Create",
			" Start,End,Strike,P/R,CpnDatas,MktDatas,[FixEnd],[FixDC],[ResetGap],[Lvge],[CpnMin],[CpnMax],[FundSpd],[FundDatas],[Nominal],[ExerGap],[NbNCall],[ExerFee],[Flags],[FundingNominal]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CRF Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq(M,Q,S,A),DC,Timing(ADV), Idxterm(1M,3M,6M,1Y),IndexDC(A360),ResetCal(EUR),PayCal (EUR),StubRule (SS)",
            " Fix leg end date (start date)",
            " Fix leg day count (coupon DC)",
            " Coupon reset lag (2 days)",
            " Leverage (1.0)",
            " Coupon min (0.0)",
            " Coupon max (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
            " Notification lag (reset Lag)",
            " Number of non call period (0)",
            " Exercise fees (0.0)",
            " FundingNominal (Nominal*1/FxSpot)"
    },
    {
        	" Local_CRFCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CRFCalculator_GetData",
			" CRFObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRF Calculator",
			" CRF Id",
            " Type = SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_PXL_CRFCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CRFCalculator_GetData",
			" CRFObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRF Calculator",
			" CRF Id",
            " Type = SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO"
    },
    {
        	" Local_CRFCalculator_GetMRS",	/// name of the C++ function
            " RR",	                        /// 2 parametres = 1 d'entree + 1 de retour 
            " ARM_GC_CRFCalculator_GetMRS",
			" CRFObject",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get mean reversion from a CRF Calculator",
			" CRF Id"
    },
    {
        	" Local_PXL_CRFCalculator_GetMRS",	/// name of the C++ function
            " RR",	                            /// 2 parametres = 1 d'entree + 1 de retour 
            " PXL_ARM_GC_CRFCalculator_GetMRS",
			" CRFObject",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get mean reversion from a CRF Calculator",
			" CRF Id"
    },
    {
        	" Local_CRFCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CRFCalculator_SetData",
			" CRFObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRF Calculator",
			" CRF Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
    {
        	" Local_PXL_CRFCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_CRFCalculator_SetData",
			" CRFObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRF Calculator",
			" CRF Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
    {
        	" Local_CRFCalculator_SetAutoCalFlags", /// name of the C++ function
            " RRR",							        /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_CRFCalculator_SetAutoCalFlags",
			" CRFObject, Autocal flags",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " CRF Autocalibration flags management",
			" CRF object Id",
			" Autocal flags (Y or N) = diag OSW (Y), short term vanilla (Y), cap (Y), floor (N), ProductName (CRF), QGMFlag (False), SkewFlag (False), Shift (-0.5)",
    },
	{
        	" Local_CRFCalculator_SetAutoCalFlagsAndClone", /// name of the C++ function
            " RRR",							        /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_CRFCalculator_SetAutoCalFlagsAndClone",
			" CRFObject, Autocal flags",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " CRF Autocalibration flags management",
			" CRF object Id",
			" Autocal flags (Y or N) = diag OSW (Y), short term vanilla (Y), cap (Y), floor (N), ProductName (CRF), QGMFlag (False), SkewFlag (False), Shift (-0.5)",
    },
	{
        	" Local_PXL_CRFCalculator_SetAutoCalFlagsAndClone", /// name of the C++ function
            " RRR",							        /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_GC_CRFCalculator_SetAutoCalFlagsAndClone",
			" CRFObject, Autocal flags",
            " 0",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " CRF Autocalibration flags management",
			" CRF object Id",
			" Autocal flags (Y or N) = diag OSW (Y), short term vanilla (Y), cap (Y), floor (N), ProductName (CRF), QGMFlag (False), SkewFlag (False), Shift (-0.5)",
    },
	{
        	" Local_CRFCalculator_SetOneCalFlag", /// name of the C++ function
            " RRRR",							  /// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GC_CRFCalculator_SetOneCalFlag",
			" CRFObject, Flag, Type",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " CRF calibration flag management",
			" CRF object Id",
			" string: flag (default = NO)",
			" string: flag type (default = SKEW)",
    },
    {
        	" Local_CRFCalculator_SetProductToPrice",   /// name of the C++ function
            " RRR",							            /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_CRFCalculator_SetProductToPrice",
			" CRFObject, ProductName",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " CRF implicit product selection for pricing",
			" CRF object Id",
			" Product Name = CRF, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA",
    },
	{
        	" Local_MktDataManager_Create",	/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_Create",
			" asOfDate",
            " 1",							/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Creates a mkt data manager object",
			" date",
    },
    {
        	" Local_PXL_MktDataManager_Create",		/// name of the C++ function
            " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_Create",
			" asOfDate",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Creates a mkt data manager object",
			" date",
    },
	{
        	" Local_MktDataManager_ZCCurveGet",		/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ZCCurveGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,source",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a zc curve from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" string (anything for instance SUMMIT or USER)",
    },
    {
        	" Local_PXL_MktDataManager_ZCCurveGet",	/// name of the C++ function
            " RRRRRRR",								/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_ZCCurveGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,source",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a zc curve from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_ZCCurveSet",		/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ZCCurveSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,source,zcCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a zc curve object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
	{
        	" Local_PXL_MktDataManager_ZCCurveSet",		/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_ZCCurveSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,source,zcCurveId",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a zc curve object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
	{
        	" Local_MktDataManager_ZCCurveGetKey",	/// name of the C++ function
            " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ZCCurveGetKey",
			" indexName,ccy,cvName,asOf,source,zcCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets the market data manager key for a zc curve",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_VolCurveGet",	/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolCurveGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,volType,source",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a vol curve object from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, ATMANDSMILE",
			" string (anything for instance SUMMIT or USER)",
    },
    {
        	" Local_PXL_MktDataManager_VolCurveGet",/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_VolCurveGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,volType,source",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a vol curve from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, ATMANDSMILE",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_VolCurveSet",	/// name of the C++ function
            " RRRRRRRRRR",							/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolCurveSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,volType,source,volCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol curve object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, SMILE",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
	{
        	" Local_MktDataManager_VolCurveSet",	/// name of the C++ function
            " RRRRRRRRRR",							/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolCurveSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,volType,source,volCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol curve object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, SMILE",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
	{
        	" Local_PXL_MktDataManager_VolCurveSet",	/// name of the C++ function
            " RRRRRRRRRR",								/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_VolCurveSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,volType,source,volCurveId",
            " 0",										/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol curve object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, SMILE",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
    {
        	" Local_MktDataManager_VolCurveGetKey",/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolCurveGetKey",
			" indexName,ccy,cvName,asOf,volMktType,volType,source",
            " 1",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets the market data manager key for a vol curve ",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" volType like ATM, ATMANDSMILE",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_VolMktModelGet",	/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolMktModelGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,source",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a vol mkt model object from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" string (anything for instance SUMMIT or USER)",
    },
    {
        	" Local_PXL_MktDataManager_VolMktModelGet",/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_VolMktModelGet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,source",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets a vol mkt model object from a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_VolMktModelSet",	/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolMktModelSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,source,volCurveId",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol mkt model object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
	{
        	" Local_PXL_MktDataManager_VolMktModelSet",	/// name of the C++ function
            " RRRRRRRRR",								/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_VolMktModelSet",
			" mktDataManagerId,indexName,ccy,cvName,asOf,volMktType,source,volCurveId",
            " 0",										/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol mkt model object to a mkt data manager",
			" Object",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" string (anything for instance SUMMIT or USER)",
			" Object",
    },
    {
        	" Local_MktDataManager_VolMktModelGetKey",/// name of the C++ function
            " RRRRRRRR",							/// 8 parametres = 7 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_VolMktModelGetKey",
			" indexName,ccy,cvName,asOf,volMktType,source",
            " 1",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets the market data manager key for a vol mkt model ",
			" string like EURIB",
			" string like EUR",
			" string like MO",
			" date like today()",
			" vol Mkt Type like IRG,SWOPT",
			" string (anything for instance SUMMIT or USER)",
    },
	{
        	" Local_MktDataManager_MktDataGet",		/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_MktDataGet",
			" mktDataManagerId,XL_objectKey",
            " 1",									/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets an object registered with its key from a mkt data manager",
			" Object",
			" string representing the referencing key",
    },
    {
        	" Local_PXL_MktDataManager_MktDataGet",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_MktDataGet",
			" mktDataManagerId,objectKey",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Gets an object registered with its key from a mkt data manager",
			" Object",
			" string representing the referencing key",
    },
	{
        	" Local_MktDataManager_MktDataSet",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_MktDataSet",
			" mktDataManagerId,objectKey,object",
            " 1",								/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol mkt model object to a mkt data manager",
			" Object",
			" string representing the referencing key",
			" object",
    },
	{
        	" Local_PXL_MktDataManager_MktDataSet",	/// name of the C++ function
            " RRRR",								/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_MktDataMger_MktDataSet",
			" mktDataManagerId,objectKey,object",
            " 0",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets a vol mkt model object to a mkt data manager",
			" Object",
			" string representing the referencing key",
			" object",
    },
    {
        	" Local_MktDataManager_MktDataFill",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_MktDataFill",
			" mktDataManagerId,objectKey,object",
            " 1",								/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Fills a vol mkt model object to a mkt data manager",
			" Object",
			" string representing the referencing key",
			" object",
    },
	{
        	" Local_MktDataManager_SetDetailMode",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_SetDetailMode",
			" mktDataManagerId,detailMode",
            " 1",								/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Sets the detail mode on a mkt data manager (for the view method)",
			" Object",
			" bool value",
    },
	{
        	" Local_EventViewer_Create",		/// name of the C++ function
            " R",								/// 1 parametres = 0 d'entree + 1 parametre de retour 
            " ARM_GP_EventViewer_Get",
			" ",
            " 1",								/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets the event viewer",
    },
	{
        	" Local_PXL_EventViewer_Create",	/// name of the C++ function
            " R",								/// 1 parametres = 0 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_EventViewer_Get",
			" ",
            " 0",								/// not visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets the event viewer",
    },
	{
        	" Local_EventViewer_ResetMssg",		/// name of the C++ function
            " RR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_EventViewer_ResetMssg",
			" EvtViewer",
            " 1",								/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Resets message in the event viewer",
			" Object",
    },
	{
        	" Local_EventViewer_SetVerboseMode",		/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_EventViewer_SetVerboseMode",
			" verboseMode",
            " 1",								/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Sets verbose mode in the event viewer",
			" boolean",
    },
	{
        	" Local_ErrViewer_Create",			/// name of the C++ function
            " RR",								/// 1 parametres = 0 d'entree + 1 parametre de retour 
            " ARM_GP_ErrViewer_Get",
			" [Reset]",
            " 1",								/// visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets the last error viewer",
			" reset flag"
    },
	{
        	" Local_PXL_ErrViewer_Create",		/// name of the C++ function
            " RR",								/// 1 parametres = 0 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_ErrViewer_Get",
			" [Reset]",
            " 0",								/// not visible in excel
            XLLOCALARM_UTIL_GROUP,
            " ",
            " ",
            " Gets the last error viewer",
			" reset flag"
    },
	{
        	" Local_MktDataManager_ChangeAsOf",	/// name of the C++ function
            " RRR",								/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ChangeAsOf",
			" mktDataManagerId, as of date",
            " 1",								/// visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Changes the as of date in the mkt data manager",
			" object",
			" as Of Date"
    },
	{
        	" Local_MktDataManager_ResetMyData",	/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ResetMyData",
			" mktDataManagerId",
            " 1",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Resets the mkt data manager",
			" object"
    },
	{
        	" Local_MktDataManager_ResetAllData",	/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_MktDataMger_ResetAllData",
			" mktDataManagerId",
            " 1",									/// not visible in excel
            XLLOCALARM_MKTMANAGER_GROUP,
            " ",
            " ",
            " Resets all data in the mkt data manager",
			" object"
    },
	{
        	" Local_GenCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_GenCalculator_GetPricingData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a generic calculator",
			" Generic Calculator Object",
			" string"
    },
	{
        	" Local_GenCalculator_GetData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_GenCalculator_GetData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from any calculator",
			" Calculator Object",
			" string : SECURITY, MODEL, CALIBMETHOD, MDM"
    },
	{
        	" Local_PXL_GenCalculator_GetData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GC_GenCalculator_GetData",
			" Gen Calc, key",
            " 0",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from any calculator",
			" Calculator Object",
			" string : SECURITY, MODEL, CALIBMETHOD, MDM"
	},
	{
        	" Local_GenCalculator_SetData",	/// name of the C++ function
            " RRRR",	                    /// 4 parametres = 3 d'entree + 1 de retour 
            " ARM_GC_GenCalculator_SetData",
			" Calculator Object,Object to set,[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Set object to a calculator",
			" Calculator Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM",
			" Market data keys"
	},
	{
        	" Local_PXL_GenCalculator_SetData",	/// name of the C++ function
            " RRRR",	                    /// 4 parametres = 3 d'entree + 1 de retour 
            " PXL_ARM_GC_GenCalculator_SetData",
			" Calculator Object,Object to set,[MktDataKeys]",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Set object to a calculator",
			" Calculator Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM",
			" Market data keys"
	},
    {
        	" Local_ARM_INITCRF",				/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",				/// 20 parametres = 19 d'entree + 1 de retour 
            " ARM_INITCRF",
			" CRF,ZcCpn,swaptionVol,capVol,[rhoCap],[nuCap],[ZcFund],[ZcCpnBasis],[ZcFundBasis],[Fx],[isUpdate],[betaCap],[rhoSwopt],[nuSwopt],[betaSwopt],[modelType],[meanRev],[skewRecalFlag],[SigmaOrAlpha]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CRF Calculator",
			" CRF",
            " Zc Coupon Leg",
            " swoptVol",
            " ATM cap Vol",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" Zc Funding Leg (default NULL)",
			" Zc Coupon Basis (default NULL)",
			" Zc Funding Basis (default NULL)",
			" Forex (default NULL)",
			" is Update or Init (default Init)",
			" beta Cap (default NULL)",
			" rho Swopt (default NULL)",
			" nu Swopt (default NULL)",
			" beta Swopt (default NULL)",
			" Model type (default HWM1F), QGM1F",
			" Mean reversion : ReferenceValue (default : value in Summit DB)",
            " Skew Recalib Flag:Yes,Y,NO,N, By default irrelevant",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
        	" Local_PXL_ARM_INITCRF",				/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",				/// 20 parametres = 19 d'entree + 1 de retour 
            " PXL_ARM_INITCRF",
			" CRF,ZcCpn,swaptionVol,capVol,[rhoCap],[nuCap],[ZcFund],[ZcCpnBasis],[ZcFundBasis],[Fx],[isUpdate],[betaCap],[rhoSwopt],[nuSwopt],[betaSwopt],[modelType],[meanRev],[skewRecalFlag],[SigmaOrAlpha]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CRF Calculator",
			" CRF",
            " Zc Coupon Leg",
            " swoptVol",
            " ATM cap Vol",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" Zc Funding Leg (default NULL)",
			" Zc Coupon Basis (default NULL)",
			" Zc Funding Basis (default NULL)",
			" Forex (default NULL)",
			" is Update or Init (default Init)",
			" beta Cap (default NULL)",
			" rho Swopt (default NULL)",
			" nu Swopt (default NULL)",
			" beta Swopt (default NULL)",
			" Model type (default HWM1F), QGM1F",
            " Mean reversion : ReferenceValue (default: value in Summit DB)",
            " Skew Recalib Flag:Yes,Y,NO,N, By default irrelevant",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
        	" PXL_Local_CRFCalculator_Initialize",	/// name of the C++ function
            " RRRRRRRRRRRR",						    /// 10 parametres = 9 d'entree + 1 parametre de retour 
            " PXL_ARM_GC_CRFCalculator_Initialize",
            " calculatorId, mktManagerId,[toCalSigma],[toCalMrs],[MRSStrikeType],[toAdjKcap],[toAdjKfloor],[modelType],[toCalSkew],[Kshift],[frontier]",
			
            " 0",									/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " "
    },
	{
        	" Local_CRFCalculator_Initialize",	/// name of the C++ function
            " RRRRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GC_CRFCalculator_Initialize",
			" calculatorId, mktManagerId,[toCalSigma],[toCalMrs],[MRSStrikeType],[toAdjKcap],[toAdjKfloor],[pricingModel],[toCalSkew],[Kshift],[frontier]",
            " 1",									
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialise a GC_Calculator with a given mkt Data Manager",
			" string: Generic Calculator",
			" string: Market Data Manager",
            " string: to calibrate instantanous volatility, default(YES),NO",
            " string: to calibrate Mean Reversion, default(NO),YES",
			" string: strike type to calibrate MRS, default(EQUIVALENT),ATM,MOYENESS",
            " string: to adjust caplet strikes , default(YES),NO",
            " string: to adjust floorlet strikes , default(NO),YES",
            " string: to choose pricing model, default(HWM1F), QGM1F",
            " string: to calibrate QGM Skew, default(NO),YES",
            " value:  default(25%), Kskew=(1-value)Keq+value*Katm",
			 "value:  default(0), Nb of tree running to calculate frontier"
    },
	{
        	" Local_CRFCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CRFCalculator_Update",
			" CRFObject, Object, [PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a CRF Calculator by the given object",
			" CRF object Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
			" Market Data Keys" 
    },
	{
        	" Local_GenCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_Calculator_Update",
			" CalObject, Object, [PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a Gen Calculator by the given object",
			" string: Calculator Id",
            " string: Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " string: for portfolio, type = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO,....",
			" vector of string: Market Data Keys" 
    },
	{
        	" Local_ARM_INITMATCAP",	/// name of the C++ function
            " RRRRRRRRR",				/// 9 parametres = 8 d'entree + 1 de retour 
            " ARM_INITMATCAP",
			" MatCap,Zc,capVol,[rhoCap],[nuCap],[betaCap],[nbpas],[SigmaOrAlpha]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Maturity Cap Calculator",
			" MATCAP",
            " Zc",
            " ATM cap Vol",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" beta Cap (default NULL)",
			" nb of iterations (default: 10000)",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
        	" Local_PXL_ARM_INITMATCAP",/// name of the C++ function
            " RRRRRRRRR",				/// 9 parametres = 8 d'entree + 1 de retour 
            " PXL_ARM_INITMATCAP",
			" MatCap,Zc,capVol,[rhoCap],[nuCap],[betaCap],[nbpas],[SigmaOrAlpha]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Maturity Cap Calculator",
			" MATCAP",
            " Zc",
            " ATM cap Vol",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" beta Cap (default NULL)",
			" nb of iterations (default: 10000)",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
    {
        	" Local_TARNCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_TARNCalculator_Create",
 			" Start,End,Strike,P/R,CpnDatas,[MktDatas],[ResetGap],[IntRule],[Lvge],[LifeTimeCapDatas],[LifeTimeFloor],[FundSpd],[FundDatas],[Nominal],[NbIterations],[calibFlags],[outputFlags],[FundNominal],[Fees],[AsOf]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a TARN Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq (M,Q,S,A), DC, Timing (ADV), Idx term (1M,3M,6M,1Y), Index DC (A360), Reset Cal (EUR), Pay Cal (EUR)",
            " Market Datas = Mkt Data Mger Id, Crv Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key",
            " Coupon reset lag (2 days)",
			" Interest Term Rules (ADJ)",
            " Coupon Curves : Leverage, Cpn Min (0.0), Cpn Max (0.0)",
			" life time cap datas (0.0,Y)",
            " life time floor (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
			" Nb Iterations (10000)",
            " AutoCal Flags = cap floor  (Y), digital(N), diag OSW (N), exer strikes (MC)",
            " Pricing Flags = TARN (Y), Swap (N), LifeTimeCap (N), LifeTimeFloor (N), DigitalFunding (N), Funding (N), Exer Strikes (N), Exer Probas (N), ExerciseTimes (N)",
 			" Funding Nominal (Nominal)",
			" Fees (1)",
            " AsOf"
    },
    {
        	" Local_PXL_TARNCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_TARNCalculator_Create",
 			" Start,End,Strike,P/R,CpnDatas,[MktDatas],[ResetGap],[IntRule],[Lvge],[LifeTimeCapDatas],[LifeTimeFloor],[FundSpd],[FundDatas],[Nominal],[NbIterations],[calibFlags],[outputFlags],[FundNominal],[Fees],[AsOf]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a TARN Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq (M,Q,S,A), DC, Timing (ADV), Idx term (1M,3M,6M,1Y), Index DC (A360), Reset Cal (EUR), Pay Cal (EUR)",
            " Market Datas = Mkt Data Mger Id, Crv Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key",
            " Coupon reset lag (2 days)",
			" Interest Term Rules (ADJ)",
            " Coupon Curves : Leverage, Cpn Min (0.0), Cpn Max (0.0)",
			" life time cap datas (0.0,Y)",
            " life time floor (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
			" Nb Iterations (10000)",
            " AutoCal Flags = cap floor  (Y), digital(N), diag OSW (N), exer strikes (MC)",
            " Pricing Flags = TARN (Y), Swap (N), LifeTimeCap (N), LifeTimeFloor (N), DigitalFunding (N), Funding (N), Exer Strikes (N), Exer Probas (N), ExerciseTimes (N)",
 			" Funding Nominal (Nominal)",
			" Fees (1)",
            " AsOf"
    },
	{
        	" Local_ARM_INITTARN",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",				/// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_INITTARN",
			" TARN,ZcCpn,swaptionVol,capVol,[rhoCap],[nuCap],[betaCap],[rhoSwopt],[nuSwopt],[betaSwopt],[ZcFund],[ZcCpnBasis],[ZcFundBasis],[Fx],[modelType],[betaCorrel],[hump],[SigmaOrAlpha]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialize a TARN Calculator with market curves",
			" TARN",
            " Zc Coupon Leg",
            " swoptVol",
            " cap Vol (ATM if SABR)",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" beta Cap (default NULL)",
			" rho Swopt (default NULL)",
			" nu Swopt (default NULL)",
			" beta Swopt (default NULL)",
			" Zc Funding Leg (default NULL)",
			" Zc Coupon Basis (default NULL)",
			" Zc Funding Basis (default NULL)",
			" Forex (default NULL)",
            " modelType: SFRM2F(Default),SBGM"
            " betaCorrel : in SBGM case",
            " hump: in SBGM case",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
        	" Local_PXL_ARM_INITTARN",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",				/// 19 parametres = 18 d'entree + 1 de retour 
            " PXL_ARM_INITTARN",
			" TARN,ZcCpn,swaptionVol,capVol,[rhoCap],[nuCap],[betaCap],[rhoSwopt],[nuSwopt],[betaSwopt],[ZcFund],[ZcCpnBasis],[ZcFundBasis],[Fx],[modelType],[betaCorrel],[hump],[SigmaOrAlpha]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialize a TARN Calculator with market curves",
			" TARN",
            " Zc Coupon Leg",
            " swoptVol",
            " cap Vol (ATM if SABR)",
			" rho Cap (default NULL)",
			" nu Cap (default NULL)",
			" beta Cap (default NULL)",
			" rho Swopt (default NULL)",
			" nu Swopt (default NULL)",
			" beta Swopt (default NULL)",
			" Zc Funding Leg (default NULL)",
			" Zc Coupon Basis (default NULL)",
			" Zc Funding Basis (default NULL)",
			" Forex (default NULL)",
            " modelType: SFRM2F(Default),SBGM",
            " betaCorrel : in SBGM case",
            " hump: in SBGM case",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
        	" Local_ARM_INITTARNFX",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",				/// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_INITTARNFX",
			" TARN,Zc,Basis,Forex,ATM,[BSFxVol],[MixtureParams],MRS,Q,CorrelMatrix,MCParams,ModelType,[nbFactor],[PDEParams],[rescaling],[smileFlag],[mixCalib],[oneFactor],[correlType]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialize a TARN Calculator with market curves",
			" TARN",
            " vector of ZC curves (1 dom + N fgn)",
            " vector of basis curves (1 dom + N fgn)",
            " vector of forex (1 dom + N fgn)",
			" vector of ATM swopt (1 dom + N fgn)",
			" vector of BS FX vol (empty if Mixture, N fx)",
			" vector of mixture params (empty if BS, N fx)",
			" vector of MRS params (1 dom + N fgn)",
			" vector of Q params (N fgn)",
			" ",
			" nb simul, bucket size, random generator, 1st nb dims, 1st nb times",
			" ",
			" default : 6",
			" timeStepNb (1000), spaceStepNb (901), stdDevNb (6), skipPDE (N)",
            " default : N",
            " default : Y",
            " default : N",
            " default : N",
            " default : CorrelMatrix"
    },
	{
        	" Local_PXL_ARM_INITTARNFX",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",				/// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_INITTARNFX",
			" TARN,Zc,Basis,Forex,ATM,[BSFxVol],[MixtureParams],MRS,Q,CorrelMatrix,MCParams,ModelType,[nbFactor],[PDEParams],[rescaling],[smileFlag],[mixCalib],[oneFactor],[correlType]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " "
    },
    {
        	" Local_TARNSnowBallCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_TARNSnowBallCalculator_Create",
			" Start,End,Strike,C0,P/R,CpnDatas,MktDatas,[ResetGap],[IntRule],[Leverage],[LevPrev],[LifeTimeCapDatas],[LifeTimeFloor],[FundSpd],[FundDatas],[Nominals],[NbIterations],[calibFlags],[outputFlags],[FundCcy]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a TARN SnowBall Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " first coupon",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq (M,Q,S,A), DC, Timing (ADV), Idx term (1M,3M,6M,1Y), Index DC (A360), Reset Cal (EUR), Pay Cal (EUR)",
            " Market Datas = Mkt Data Mger Id, Crv Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key",
            " Coupon reset lag (2 days)",
			" Interest Term Rules (ADJ)",
            " Coupon Curves : Leverage, Cpn Min (0.0), Cpn Max (0.0)",
			" LevPrev (1.0)",
			" life time cap datas (0.0,Y)",
            " life time floor (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " vector of string :Nominal (100), funding notional (default = Nominal)",
			" Nb Iterations (10000)",
            " AutoCalFlags = cap floor  (Y), digital(N), diag OSW (N), exer strikes (MC)",
            " Pricing Flags = TARN (Y), Swap (N), LifeTimeCap (N), LifeTimeFloor (N), DigitalFunding (N), Funding (N), Exer Strikes (N), Exer Probas (N), ExerciseTimes (N)"
			" Funding Ccy (Cpn Ccy)"
    },
    {
        	" Local_PXL_TARNSnowBallCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_TARNSnowBallCalculator_Create",
			" Start,End,Strike,C0,P/R,CpnDatas,MktDatas,[ResetGap],[IntRule],[Lvge],[LevPrev],[LifeTimeCapDatas],[LifeTimeFloor],[FundSpd],[FundDatas],[Nominal],[NbIterations],[calibFlags],[outputFlags],[FundCcy]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a TARN SnowBall Calculator",
			" Start date",
            " End date",
            " Fix rate or strike",
            " first Coupon",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq (M,Q,S,A), DC, Timing (ADV), Idx term (1M,3M,6M,1Y), Index DC (A360), Reset Cal (EUR), Pay Cal (EUR)",
            " Market Datas = Mkt Data Mger Id, Crv Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key",
            " Coupon reset lag (2 days)",
			" Interest Term Rules (ADJ)",
            " Coupon Curves : Leverage, Cpn Min (0.0), Cpn Max (0.0)",
			" LevPrev (1.0)",
			" life time cap datas (0.0,Y)",
            " life time floor (10000.0)",
            " Funding leg spread (0.0)",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " Nominal (100)",
			" Nb Iterations (10000)",
            " AutoCalFlags = cap floor  (Y), digital(N), diag OSW (N), exer strikes (MC)",
            " Pricing Flags = TARN (Y), Swap (N), LifeTimeCap (N), LifeTimeFloor (N), DigitalFunding (N), Funding (N), Exer Strikes (N), Exer Probas (N), ExerciseTimes (N)",
			" Funding Ccy (Cpn Ccy)"
    },
	{
        	" Local_TARNCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_TARNCalculator_GetData",
			" TARNObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a TARN Calculator",
			" TARN Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO"
    },
    {
        	" Local_PXL_TARNCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_TARNCalculator_GetData",
			" TARNObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a TARN Calculator",
			" TARN Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO"
    },
    {
        	" Local_TARNCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_TARNCalculator_SetData",
			" TARNObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a TARN Calculator",
			" TARN Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
    {
        	" Local_PXL_TARNCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 5 d'entree + 1 de retour 
            " PXL_ARM_GC_TARNCalculator_SetData",
			" TARNObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a TARN Calculator",
			" TARN Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_TARNCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_TARNCalculator_Update",
			" TARNObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a TARN Calculator by the given object",
			" TARN object Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_TARNCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_TARNCalculator_GetPricingData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a tarn calculator",
			" TARN Calculator Object",
			" string"
    },
	{
        	" Local_MaturityCapCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_MaturityCapCalculator_Create",
			" Start, End, UnderlyingEnd, [L/S], [C/F], LoanDatas, [Spread], [InitNominal], InitTRI, Annuity, [MatCapMode], [Coeff], [Amortizing], [ResetGap], [ResetCal], [PayCal], [CalibMode], [NbIterations], [Flags], MktDatas",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Maturity Cap Calculator",
			" Start date",
            " End date",
			" Underlying End date",
			" Long Short (L)",
			" Cap Floor (C)",
			" Loan Datas = Reset Freq, PayFreq, IndexTerm, DayCount, IntRule",
            " Spread (0.0)",
            " Init Nominal (100.0)",
			" Init TRI",
            " Annuity",
            " Maturity Cap Mode (TRI)",
            " Coeff (1.0)",
			" Amortizing (1.0)",
            " Reset Gap (-2)",
            " Reset Calendar (EUR)",
            " Payment Calendar (EUR)",
            " Calibration mode (ATM)",
			" Nb Iterations (10000)",
            " Pricing Flags = maturity cap  (Y), ref std cap (N), esim TRI (N), estim nom (N)",
			" Market Datas = Mkt Data Mger Id, CF Mkt Model Key, MRS Key, Beta Key"
    },
    {
        	" Local_PXL_MaturityCapCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_MaturityCapCalculator_Create",
			" Start, End, UnderlyingEnd, [L/S], [C/F], LoanDatas, [Spread], [InitNominal], InitTRI, Annuity, [MatCapMode], [Coeff], [Amortizing], [ResetGap], [ResetCal], [PayCal], [CalibMode], [NbIterations], [Flags], MktDatas",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Maturity Cap Calculator",
			" Start date",
            " End date",
			" Underlying End date",
			" Long Short (L)",
			" Cap Floor (C)",
			" Loan Datas = Reset Freq, PayFreq, IndexTerm, DayCount (A360), IntRule (ADJ)",
            " Spread (0.0)",
            " Init Nominal (100.0)",
			" Init TRI",
            " Annuity",
            " Maturity Cap Mode (TRI)",
            " Coeff (0.0)",
			" Amortizing (1.0)",
            " Reset Gap (-2)",
            " Reset Calendar (EUR)",
            " Payment Calendar (EUR)",
            " Calibration mode (ATM)",
			" Nb Iterations (10000)",
            " Pricing Flags = maturity cap  (Y), ref std cap (N), esim TRI (N), estim nom (N)",
			" Market Datas = Mkt Data Mger Id, CF Mkt Model Key, MRS Key, Beta Key"
    },
	{
        	" Local_MaturityCapCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_MaturityCapCalculator_GetData",
			" MaturityCapObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Maturity Cap Calculator",
			" Maturity Cap Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO"
    },
    {
        	" Local_PXL_MaturityCapCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_MaturityCapCalculator_GetData",
			" MaturityCapObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Maturity Cap Calculator",
			" Maturity Cap Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO"
    },
    {
        	" Local_MaturityCapCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_MaturityCapCalculator_SetData",
			" MaturityCapObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Maturity Cao Calculator",
			" Maturity Cap Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
    {
        	" Local_PXL_MaturityCapCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 5 d'entree + 1 de retour 
            " PXL_ARM_GC_MaturityCapCalculator_SetData",
			" MaturityCapObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Maturity Cap Calculator",
			" Maturity Cap Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_MaturityCapCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_MaturityCapCalculator_Update",
			" MaturityCapObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a Maturity Cap Calculator by the given object",
			" Maturity Cap object Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_MaturityCapCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_MaturityCapCalculator_GetPricingData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a maturity cap calculator",
			" Maturity Cap calculator Object",
			" string"
    },
	{
        	" Local_CaptionCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",	    /// 17parametres = 16 d'entree + 1 de retour 
            " ARM_GC_CaptionCalculator_Create",
			" Start,End,MktDatas,RequiredCpnData,Coupon,FundIdxTerm,[ExerciseData],[ExerStyle],[Notional],[CpnData],[CpnSpread],[FundData],[FundSpread],[FactorNb],[CalibModData],[Flags]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Caption Calculator",
			" Start date",
            " End date",
            " MktDataMgerId, Keys(ZCCrv,CF MktModel,OSW MktModel,CF MktModel,MRS,beta,Correl,Fund ZCCrv,Basis ZCCrv)",
            " RequiredCpnData : Idxterm (6M,1Y,...), P/R, C/F",
			" Coupon of the CpnLeg(double or curve)",
			" Funding Idxterm(6M,1Y,..)",
			" ExerciseData = NotifDays (CpnCCy resetGap), ProbCalculIdx (0)",
			" Exercise fees (0.0)",
			" Notional (100)",
			" Coupon Datas = DayCount,ResetTiming, ResetCal, PayCal (cpnCcy)",
			" Spread of the CpnLeg (double or curve, 0)",
			" Funding Datas = DayCount,ResetCal, PayCal (FundCcy)",
			" Spread of the FundLeg (0)",
			" Number of Factors ofthe SFRM model",
			" Vol Type (DIAG), Calibration Mode Swaption(EXSWOPT,AllSWOPT,N,Y)",
			" Pricing Flags(Y/N) = CapPrice,CouponLegPrice,FundingLegPrice,CaptionStrikes,ExerciseProbas(N)"
    },
    {
        	" Local_PXL_CaptionCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",	    /// 17parametres = 16 d'entree + 1 de retour
            " PXL_ARM_GC_CaptionCalculator_Create",
			" Start,End,MktDatas,RequiredCpnData,Coupon,FundIdxTerm,[ExerciseData],[ExerStyle],[Notional],[CpnData],[CpnSpread],[FundData],[FundSpread],[FactorNb],[CalibModData],[Flags]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Caption Calculator",
			" Start date",
            " End date",
            " MktDataMgerId, Keys(ZCCrv,CF MktModel,OSW MktModel,CF MktModel,MRS,beta,Correl,Fund ZCCrv,Basis ZCCrv)",
            " RequiredCpnData : Idxterm (6M,1Y,...), P/R, C/F",
			" Coupon of the CpnLeg(double or curve)",
			" Funding Idxterm(6M,1Y,..)",
			" ExerciseData = NotifDays (CpnCCy resetGap), ProbCalculIdx (0)",
			" Exercise fees (0.0)",
			" Notional (100)",
			" Coupon Datas = DayCount,ResetTiming, ResetCal, PayCal (cpnCcy)",
			" Spread of the CpnLeg (double or curve, 0)",
			" Funding Datas = DayCount,ResetCal, PayCal (FundCcy)",
			" Spread of the FundLeg (0)",
			" Number of Factors ofthe SFRM model",
			" Vol Type (DIAG), Calibration Mode Swaption(EXSWOPT,AllSWOPT,N,Y)",
			" Pricing Flags(Y/N) = CapPrice,CouponLegPrice,FundingLegPrice,CaptionStrikes,ExerciseProbas(N)"

    },
	{
			" Local_ARM_INITCAPTION",	/// name of the C++ function
			" RRRRR",				/// 5 parametres = 4 d'entree + 1 de retour 
			" ARM_INITCAPTION",
			" Caption, MktDataManager, CalibParams, ProductsToPrice",
			" 1",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialise an 'empty' Caption Calculator",
			" Caption Id : the 'empty' caption leg",
            " MktDataManager : vector containing MarketDataManager and market data keys",
			" CalibParams : Factor number, Vol type, Calib swopt mode, Calib beta mode",
			" Products to price : vector of string (Y/N)"
    },
	{
			" Local_PXL_ARM_INITCAPTION",	/// name of the C++ function
			" RRRRR",				/// 5 parametres = 4 d'entree + 1 de retour 
			" PXL_ARM_INITCAPTION",
			" Caption, MktDataManager, CalibParams, ProductsToPrice",
			" 0",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialise an 'empty' Caption Calculator",
			" Caption Id : the 'empty' caption leg",
            " MktDataManager : vector containing MarketDataManager and market data keys",
			" CalibParams : Factor number, Vol type, Calib swopt mode, Calib beta mode",
			" Products to price : vector of string (Y/N)"
    },
	{
        	" Local_CaptionCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CaptionCalculator_GetData",
			" CaptionObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Caption Calculator",
			" Caption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_PXL_CaptionCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CaptionCalculator_GetData",
			" CaptionObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Caption Calculator",
			" Caption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_CaptionCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CaptionCalculator_SetData",
			" CaptionObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Caption Calculator",
			" Caption Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys",
    },
    {
        	" Local_PXL_CaptionCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 5 d'entree + 1 de retour 
            " PXL_ARM_GC_CaptionCalculator_SetData",
			" CaptionObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Caption Calculator",
			" Caption Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_CaptionCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CaptionCalculator_Update",
			" CaptionObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a Caption Calculator by the given object",
			" Caption object Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_CaptionCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_CaptionCalculator_GetPricingData",
			" Gen Pricer, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a caption calculator",
			" Caption Calculator Object",
			" string"
    },
	{
        	" Local_PXL_PRDCCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRR",	            /// 13 parametres = 12 d'entree + 1 de retour 
            " PXL_ARM_GC_PRDCCalculator_Create",
			" PRDC Obj, Model Obj, OtherMktObjs, [SchedDatas], [TruncDatas], [ColNames], [MDSFlag], [LocalFxFlag], [CalType], [CalDatas], [BasisCalibFlag], [MarginFlag]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a PRDC Calculator",
    },
    {
        	" Local_PRDCCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRR",	            /// 13 parametres = 12 d'entree + 1 de retour 
            " ARM_GC_PRDCCalculator_Create",
			" PRDC Obj, Model Obj, OtherMktObjs, [SchedDatas], [TruncDatas], [ColNames], [MDSFlag], [LocalFxFlag], [CalType], [CalDatas], [BasisCalibFlag], [MarginFlag]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a PRDC Calculator",
			" PRDC Id",
            " Model Id",
            " Other Mkt Datas Ids : [0]=DomMRS, [1]=ForMRS, [ [2]=QFx, [3]=QDom, [4]=QFor ]",
            " [Scheduler Datas] : [0]=MinBef1st, [1]=MaxBef1st, [2]=NbPYBef1st, [3]=OptDate, [4]=NbPYBefOpt, [5]=NbPYAftOpt",
            " [Truncator Datas] : [0]=MaxStdDev, [1]=MaxIdx, [2]=MinADP, [3]=MinADPTime",
            " [ColumnsToPrice] : (default = PRDCOption)",
            " [MarkovianDriftFlag] : Y/N (default = Y)",
            " [LocalFxModelFlag] : Y/N (default = N)",
            " [CalibType] : ATM,ATSFX,ATSFXMIXED,ATSFXPROFILE,ATSFXMONEYNESS,ATSFXSHIFTED (default = ATM)",
            " [CalibDatas] : (default = [])",
            " [BasisIRCalibFlag] : Y/N (default = Y)",
			" [MarginFlag] :to convert funLeg to dom Ccy (FlowByFlow, Average)(default = FlowByFlow)",
    },
	{
        	" Local_PXL_PRCSCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",		/// 20 parametres = 19 d'entree + 1 de retour 
            " PXL_ARM_GC_PRCSCalculator_Create",
			" Start, FixEnd, End, Ccys, MktDataManager ,StructDatas, FundDatas, ExerDatas, RedemDatas, Notionals, InitialFx,Cpns, [FundCurves], [Fees], [calibTypes],[calibDatas], [ProductFlags],[SchedulerDatas],[TruncatorDatas] ",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Power reverse dual currencies swap calculator",
	},
	{
        	" Local_PRCSCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",		/// 20 parametres = 19 d'entree + 1 de retour 
            " ARM_GC_PRCSCalculator_Create",
			" Start, FixEnd, End, Ccys, MktDataManager ,StructDatas, FundDatas, ExerDatas, RedemDatas, Notionals, InitialFx,Cpns, [FundCurves], [Fees], [calibTypes],[calibDatas], [ProductFlags],[SchedulerDatas],[TruncatorDatas] ",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Power reverse dual currencies swap calculator",
			" date: Start date",
            " date: End date of fix period",
			" date: End date",
			" vector of strings: Coupon Ccy, Foriegn Ccy,Funding Ccy = default(Coupon Currency)",
			" string: Mkt dataManager ",
            " vector of strings: Struc Datas =freq,FxResretLag, ResetTiming(K_ADVANCE),RestCalendar(""),Pay Calendar(""), Stub Rule (SS)",
            " vector of strings: Funding datas = Freq (cpn freq), DC (A360)",
            " vector of strings: ExerDatas = freq, notification gap, Pay/Rec(REC),nbNoCall (0)",
			" vector of strings: RedemptionDatas type,gap, strike",
            " vector of strings: Cpn notional curve id, Fund notional = default(Cpn notional)",
			" string: initial forex Id",
            " vector of strings: domestic Cpn curve id, foreign Cpn curve id , CpnMin curve id,CpnMax curve id",
            " string: Funding margin curve id, Funding Levrage id ( default(1.0))", 
            " string: Fees curve id ",
            " vector of strings: CalibTypes = default (localFxModel(N), calib type(ATM)) ",
			" vector of strings: CalibDatas ",
            " vector of strings: ProductNames = default (PRDCOption) ",
            " vector of double : SchedulerDatas = default( min steps 1st date=4,max steps 1st date=20, Nb/year 1st date=4, optimal date=7, nb/year before opt date=2,nb/year before opt date=1)",
			" vector of double : TruncatorDatas = default( MasStdDev=4, MaxMaxIndex=20, Min Arrow-Debreu=1.0e-5, Timetotrancate=4 yf)",
    }, 
	{
        	" Local_PXL_PRDKOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",		/// 17 parametres = 16 d'entree + 1 de retour 
            " PXL_ARM_GC_PRDKOCalculator_Create",
			" Dates, Ccys, MktDataManager ,StructDatas, FundDatas, ExerDatas, RedemDatas, Notionals, CpnCurves, [FundCurves], [Fees], [calibTypes],[calibDatas], [ProductFlags],[SchedulerDatas],[TruncatorDatas] ",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Power reverse dual currencies swap calculator",
	},
	{
        	" Local_PRDKOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",		/// 17 parametres = 16 d'entree + 1 de retour 
            " ARM_GC_PRDKOCalculator_Create",
			" Dates, Ccys, MktDataManager ,StructDatas, FundDatas, ExerDatas, RedemDatas, Notionals, CpnCurves, [FundCurves], [Fees], [calibTypes],[calibDatas], [ProductFlags],[SchedulerDatas],[TruncatorDatas] ",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Power reverse dual currencies swap calculator",
			" vector of dates: Start date, End date, End date of fix period, switch date ( PRDC to PRDKO)",
			" vector of strings: Coupon Ccy, Foriegn Ccy,Funding Ccy = default(Coupon Currency)",
			" string: Mkt dataManager ",
            " vector of strings: Struc Datas =freq,FxResretLag, ResetTiming(K_ADVANCE),RestCalendar(""),Pay Calendar(""), Stub Rule (SS)",
            " vector of strings: Funding datas = Freq (cpn freq), DC (A360)",
            " vector of strings: ExerDatas = freq, notification gap, Pay/Rec(REC),nbNoCall (0)",
			" vector of strings: RedemptionDatas type,gap, strike",
            " vector of strings: Cpn notional curve id, Fund notional = default(Cpn notional)",
            " vector of strings: initial fx Curve Id, domestic Cpn curve id, foreign Cpn curve id, CpnMin curve id,CpnMax curve id, barrier curve Id",
            " string: Funding margin curve id, Funding Levrage id ( default(1.0))", 
            " string: Fees curve id ",
            " vector of strings: CalibTypes = default (localFxModel(N), calib type(ATM)) ",
			" vector of strings: CalibDatas ",
            " vector of strings: ProductNames = default (PRDCOption) ",
            " vector of double : SchedulerDatas = default( min steps 1st date=4,max steps 1st date=20, Nb/year 1st date=4, optimal date=7, nb/year before opt date=2,nb/year before opt date=1)",
			" vector of double : TruncatorDatas = default( MasStdDev=4, MaxMaxIndex=20, Min Arrow-Debreu=1.0e-5, Timetotrancate=4 yf)",
    },    
	{
        	" Local_PRDCCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_PRDCCalculator_GetData",
			" PRDC Object, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a PRDC Calculator",
			" PRDC Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, DOM(FOR)OSWPORTFOLIO, FXPORTFOLIO, FLOOREDFXPORTFOLIO, CAPPEDFXPORTFOLIO",
    },
	{
        	" Local_PXL_PRDCCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_PRDCCalculator_GetData",
			" PRDC Object, TypeToGet",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a PRDC Calculator",
			" PRDC Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, DOM(FOR)OSWPORTFOLIO, FXPORTFOLIO, FLOOREDFXPORTFOLIO, CAPPEDFXPORTFOLIO",
    },
    {
        	" Local_CallableSBCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRR",			/// 15 parametres = 14 d'entree + 1 de retour 
            " ARM_GC_CallableSBCalculator_Create",
			" Start, End, MktDatas, P/R, CpnDatas, Notional, CpnCurves, C/F, FundDatas, ExerDatas, fees, calibFlags, ProductFlags, modelDatas",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SnowBall Calculator",
			" Start date",
            " End date",
            " Market Datas = Mkt Data Mger Id, Cpn Yc Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key, Fund Yc Key, Basis Yc Key, Forex Key",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq(M,Q,S,A),DC,Timing(ADV), Idxterm(1M,3M,6M,1Y),IndexDC(A360),ResetCal(EUR),PayCal (EUR),IntRule(ADJ/UNADJ),ResetLag",
            " Notional",
            " CpnCurves (Const,LevPrev,LevNew,StrikeOpt,CpnMin,CpnMax)",
            " Cap or Floor",
            " Funding datas = Freq (cpn freq), DC (A360), FundMargin, FundCoeff",
            " ExerDatas",
            " fees",
			" calibFlags",
            " ProductFlags",
			" modelDatas"
    },
    {
        	" Local_PXL_CallableSBCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRR",		/// 15 parametres = 14 d'entree + 1 de retour 
            " PXL_ARM_GC_CallableSBCalculator_Create",
			" Start, End, MktDatas, P/R, CpnDatas, Notional, CpnCurves, C/F, FundDatas, FundMargin, ExerDatas, fees, calibFlags, ProductFlags, modelDatas",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SnowBall Calculator",
			" Start date",
            " End date",
            " Market Datas = Mkt Data Mger Id, Cpn Yc Key, OSW Mkt Model Key, CF Mkt Model Key, MRS Key, Fund Yc Key, Basis Yc Key, Forex Key",
            " Coupon Payer(P) or Receiver (R)",
            " Coupon Datas = Freq(M,Q,S,A),DC,Timing(ADV), Idxterm(1M,3M,6M,1Y),IndexDC(A360),ResetCal(EUR),PayCal (EUR),IntRule(ADJ/UNADJ),ResetLag",
            " vector of string :Nominal (100), funding notional (default = Nominal)",
            " CpnCurves (Const,LevPrev,LevNew,StrikeOpt,CpnMin,CpnMax)",
            " Cap or Floor",
            " Funding datas = Freq (cpn freq), DC (A360)",
            " funding margin",
            " ExerDatas",
            " fees",
			" calibFlags",
            " ProductFlags",
			" modelDatas"
    },
	{
        	" Local_CallableSBCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CallableSBCalculator_GetData",
			" callableSBObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Callable SnowBall Calculator",
			" Callable SnowBall Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO",
    },
	{
        	" Local_PXL_CallableSBCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CallableSBCalculator_GetData",
			" callableSBObject, TypeToGet",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Callable SnowBall Calculator",
			" Callable SnowBall Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO",
    },
	{
        	" Local_CallableSBCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CallableSBCalculator_SetData",
			" CSBObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CSB Calculator",
			" CSB Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys",
    },
    {
        	" Local_PXL_CallableSBCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 5 d'entree + 1 de retour 
            " PXL_ARM_GC_CallableSBCalculator_SetData",
			" CSBObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CSB Calculator",
			" CSB Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_CallableSBCalculator_Update",	/// name of the C++ function
            " RRRRR",				        /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CallableSBCalculator_Update",
			" CSBObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update a CSB Calculator by the given object",
			" CSB object Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type = OSWPORTFOLIO, CFPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_CallableSBCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_CallableSBCalculator_GetPricingData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a callable SnowBall calculator",
			" callable SnowBall Calculator Object",
			" string"
    },
	{
			" Local_ARM_INITCSB",	/// name of the C++ function
			" RRRRRRRRRRRRRRR",				/// 15 parametres = 14 d'entree + 1 de retour 
			" ARM_INITCSB",
			" CSB,ZcCpn,SwoptVc,CapVc,[CapRo],[CapNu],[CapBeta],[SwoptRo],[SwoptNu],[SwoptBeta],[hump],[betaCorrel],[reCorrel],[SigmaOrAlpha]",
			" 1",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSB with Market Data",
			" CSB Id",
			" Yield Curve",
			" Swopt Vol (vol cube or ATM vol)",
			" Cap Vol (vol cube or ATM vol)",
			" Cap Ro (optional)",
			" Cap Nu (optional)",
			" Cap Beta (optional)",
			" Swopt Ro (optional)",
			" Swopt Nu (optional)",
			" Swopt Beta (optional)",
			" hump",
			" betaCorrel",
			" reCorrel",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
	{
			" Local_PXL_ARM_INITCSB",	/// name of the C++ function
			" RRRRRRRRRRRRRRR",				/// 15 parametres = 14 d'entree + 1 de retour 
			" PXL_ARM_INITCSB",
			" CSB,ZcCpn,SwoptVc,CapVc,[CapRo],[CapNu],[CapBeta],[SwoptRo],[SwoptNu],[SwoptBeta],[hump],[betaCorrel],[reCorrel],[SigmaOrAlpha]",
			" 0",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSB with Market Data",
			" CSB Id",
			" Yield Curve",
			" Swopt Vol (vol cube or ATM vol)",
			" Cap Vol (vol cube or ATM vol)",
			" Cap Ro (optional)",
			" Cap Nu (optional)",
			" Cap Beta (optional)",
			" Swopt Ro (optional)",
			" Swopt Nu (optional)",
			" Swopt Beta (optional)",
			" hump",
			" betaCorrel",
			" reCorrel",
            " SABR model mode (SIGMA or ALPHA) def: SIGMA)"
    },
    {
        	" Local_BermudaSwaptionCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRR",	    /// 12 parametres = 11 d'entree + 1 de retour 
            " ARM_GC_BermudaSwaptionCalculator_Create",
			" Ccy,Start,End,NotionalCurve,StrikeCurve,FeesCurve,SpreadCurve,Pay/Receive,CallDatas,FixDatas,VarDatas,GenSecType,MktDataManager,CalibParams,[ControlVariates]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator",
			" Currency", 
			" Swap Start date",
            " Swap End date",
            " Notional Curve",
            " Strike Curve",
			" Fees Curve",
			" Spread Curve",
			" Pay / Receive",
			" Call Datas",
			" Fix Datas",
			" Var Datas",
			" GenSecType",
			" MktDataManager",
			" Calib Parameters",
			" ControlVariates"
   },
   {
        	" Local_PXL_BermudaSwaptionCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRR",	    /// 12 parametres = 11 d'entree + 1 de retour 
            " PXL_ARM_GC_BermudaSwaptionCalculator_Create",
			" Ccy,Start,End,NotionalCurve,StrikeCurve,FeesCurve,SpreadCurve,Pay/Receive,CallDatas,FixDatas,VarDatas,GenSecType,MktDataManager,CalibParams,[ControlVariates]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a BermudaSwaption Calculator",
			" Currency", 
			" Swap Start date",
            " Swap End date",
            " Notional Curve",
            " Strike Curve",
			" Fees Curve",
			" Spread Curve",
			" Pay / Receive",
			" Call Datas",
			" Fix Datas",
			" Var Datas",
			" GenSecType",
			" MktDataManager",
			" Calib Parameters",
			" ControlVariates"
	},
	{
        	" Local_BermudaSwaptionCalc_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRR",	    /// 12 parametres = 11 d'entree + 1 de retour 
            " ARM_GC_BermudaSwaptionCalc_Create",
			" Ccy,Start,End,NotionalCurve,StrikeCurve,FeesCurve,SpreadCurve,Pay/Receive,CallDatas,FixDatas,VarDatas,GenSecType,MktDataManager,CalibParams,[ControlVariates]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator with default parameters",
			" Currency", 
			" Swap Start date",
            " Swap End date",
            " Notional Curve",
            " Strike Curve",
			" Fees Curve",
			" Spread Curve",
			" Pay / Receive",
			" Call Datas",
			" Fix Datas",
			" Var Datas",
			" GenSecType",
			" MktDataManager",
			" Calib Parameters",
			" ControlVariates"
	},
	{
        	" Local_PXL_BermudaSwaptionCalc_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRR",	    /// 12 parametres = 11 d'entree + 1 de retour 
            " PXL_ARM_GC_BermudaSwaptionCal_Create",
			" Ccy,Start,End,NotionalCurve,StrikeCurve,FeesCurve,SpreadCurve,Pay/Receive,CallDatas,FixDatas,VarDatas,GenSecType,MktDataManager,CalibParams,[ControlVariates]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a BermudaSwaption Calculator with default parameters",
			" Currency", 
			" Swap Start date",
            " Swap End date",
            " Notional Curve",
            " Strike Curve",
			" Fees Curve",
			" Spread Curve",
			" Pay / Receive",
			" Call Datas",
			" Fix Datas",
			" Var Datas",
			" GenSecType",
			" MktDataManager",
			" Calib Parameters",
			" ControlVariates"
	},
	{
        	" Local_BermudaSwaptionCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_BermudaSwaptionCalculator_GetData",
			" BermudaSwaptionObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, OSWPORTFOLIO, STMPORTFOLIO, EXOSWAPTION, SWAP"
    },
    {
        	" Local_PXL_BermudaSwaptionCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_BermudaSwaptionCalculator_GetData",
			" BermudaSwaptionObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, OSWPORTFOLIO, STMPORTFOLIO, EXOSWAPTION, SWAP"
    },
    {
        	" Local_BermudaSwaptionCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_BermudaSwaptionCalculator_SetData",
			" BermudaSwaptionObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
			" Object to force",
            " Type = SECURITY, MODEL, CALIBMETHOD, OSWPORTFOLIO, STMPORTFOLIO",
			" Market datas keys"
    },
    {
        	" Local_PXL_BermudaSwaptionCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_BermudaSwaptionCalculator_SetData",
			" BermudaSwaptionObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
			" Object to force",
            " Type = SECURITY, MODEL, CALIBMETHOD, OSWPORTFOLIO, STMPORTFOLIO",
			" Market datas keys"
    },
	{
        	" Local_BermudaSwaptionCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GC_BermudaSwaptionCalculator_GetPricingData",
			" Gen Calc, key",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a bermuda swaption calculator",
			" Bermuda Swaption Calculator Object",
			" string"
    },
	{
        	" Local_BermudaSwaptionCalculator_GetPricingData",	/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GC_BermudaSwaptionCalculator_GetPricingData",
			" Gen Calc, key",
            " 0",							/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " get data from a bermuda swaption calculator",
			" Bermuda Swaption Calculator Object",
			" string"
    },
	{
			" Local_ARM_INITBERMUDASWAPTION",	/// name of the C++ function
			" RRRRR",				/// 5 parametres = 4 d'entree + 1 de retour 
			" ARM_INITBERMUDASWAPTION",
			" BermudaSwaption,MktDataManager,CalibParams,[ControlVariates]",
			" 1",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
            " MktDataManager",
			" CalibParams",
			" ControlVariates"
    },
	{
			" Local_PXL_ARM_INITBERMUDASWAPTION",	/// name of the C++ function
			" RRRRR",				/// 5 parametres = 4 d'entree + 1 de retour 
			" PXL_ARM_INITBERMUDASWAPTION",
			" BermudaSwaption,MktDataManager,CalibParams,[ControlVariates]",
			" 0",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
            " MktDataManager",
			" CalibParams",
			" ControlVariates"
    },
	{
			" Local_ARM_INITSWAPTIONBERMUDA",	/// name of the C++ function
			" RRRRRRRRRRRRRRR",			/// 15 parametres = 14 d'entree + 1 de retour 
			" ARM_INITSUMMITSWAPTIONBERMUDA",
			" Bermuda,CalibParams,ZcCpn,SwoptVc,CapVc,[CapRo],[CapNu],[CapBeta],[SwoptRo],[SwoptNu],[SwoptBeta],[NormalModel],[ControlVariates],[SigmaOrAlpha]",
			" 1",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
			" CalibParams",
			" Yield Curve",
			" Swopt Vol (vol cube or ATM vol)",
			" Cap Vol (vol cube or ATM vol)",
			" Cap Ro (optional)",
			" Cap Nu (optional)",
			" Cap Beta (optional)",
			" Swopt Ro (optional)",
			" Swopt Nu (optional)",
			" Swopt Beta (optional)",
			" Normal Model (optional)",
			" ControlVariates (optional)",
			" SABR Vol type : S (default), SIGMA, A, ALPHA"
    },
	{
			" Local_PXL_ARM_INITSWAPTIONBERMUDA",	/// name of the C++ function
			" RRRRRRRRRRRRRRR",			/// 15 parametres = 14 d'entree + 1 de retour 
			" PXL_ARM_INITSUMMITSWAPTIONBERMUDA",
			" Bermuda,CalibParams,ZcCpn,SwoptVc,CapVc,[CapRo],[CapNu],[CapBeta],[SwoptRo],[SwoptNu],[SwoptBeta],[NormalModel],[ControlVariates],[SigmaOrAlpha]",
			" 0",						/// visible in excel
        	XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Bermuda Swaption Calculator",
			" Bermuda Swaption Id",
			" CalibParams",
			" Yield Curve",
			" Swopt Vol (vol cube or ATM vol)",
			" Cap Vol (vol cube or ATM vol)",
			" Cap Ro (optional)",
			" Cap Nu (optional)",
			" Cap Beta (optional)",
			" Swopt Ro (optional)",
			" Swopt Nu (optional)",
			" Swopt Beta (optional)",
			" Normal Model (optional)",
			" ControlVariates (optional)",
			" SABR Vol type : S (default), SIGMA, A, ALPHA"
    },
	{
        	" Local_BermudaSwaptionCalculator_RootMrs",	/// name of the C++ function
            " RRRRR",									/// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_BermudaSwaptionCalculator_RootMrs",
			" BermudaSwaptionObject, TargetPrice, Tolerance, MaxIter",
            " 1 ",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Calculate MRS level corresponding to a target price",
			" Bermuda Swaption Id",
			" Bermuda Target Price",
			" F Tolerance",
			" Maximum Iterations"
    },
    {
        	" Local_PXL_BermudaSwaptionCalculator_RootMrs",	/// name of the C++ function
            " RRRRR",										/// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_BermudaSwaptionCalculator_RootMrs",
			" BermudaSwaptionObject, TargetPrice, Tolerance, MaxIter",
            " 0 ",											/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Calculate MRS level corresponding to a target price",
			" Bermuda Swaption Id",
			" Bermuda Target Price",
			" F Tolerance",
			" Maximum Iterations"
    },
	{		" Local_CRALocalCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 20 parametres = 19 d'entree + 1 de retour 
            " ARM_GC_CRALocalCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,CpnSpread,BoostedFix,BarrierDown,BarrierUp,LocalModelParams,MRSBeta,calibSecPFParams, NbSteps,flagToGenerateOSWATM,MktDataManager,ProductsToPrice,[IsSfrmStdCalib]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" General Datas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDates / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Local Model Parameters",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, [Recalibrate MRS flag, Recalibrate Beta flag]",
			" Calibration parameters",
			" Number of steps",
			" ATM swaption flag",
			" MktDataManager",
			" Products to price",
			" Is a standard SFRM calibration"
   },
   {
        	" Local_PXL_CRALocalCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 20 parametres = 19 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,CpnSpread,BoostedFix,BarrierDown,BarrierUp, LocalModelParameters, MRSBeta,calibSecPFParams,NbSteps,flagToGenerateOSWATM,MktDataManager,ProductsToPrice,[IsSfrmStdCalib]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" General Datas", 
		   	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Local Model Parameters",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, Recalibrate MRS flag,Recalibrate Beta flag",
			" Calibration parameters",
			" Number of steps",
			" ATM swaptions flag",
			" MktDataManager",
			" Products to price"
			" Is a standard SFRM calibration"
	},
	{		" Local_CRACalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_CRACalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,CpnSpread,BoostedFix,BarrierDown,BarrierUp,MRSBeta,calibSecPFParams,nbSteps,flagToGenerateOSWATM,MktDataManager,ProductsToPrice,[IsSfrmStdCalib]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, [Recalibrate MRS flag, Recalibrate Beta flag]",
			" Calibration parameters",
			" Tree steps number",
			" ATM swaptions flag",
			" MktDataManager",
			" Products to price",
			" Is a standard SFRM calibration"
   },
   {
        	" Local_PXL_CRACalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	    /// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,CpnSpread,BoostedFix,BarrierDown,BarrierUp,MRSBeta,calibSecPFParams,nbSteps,flagToGenerateOSWATM,MktDataManager,ProductsToPrice,[IsSfrmStdCalib]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" General Datas", 
		   	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, [Recalibrate MRS flag, Recalibrate Beta flag]",
			" Calibration parameters",
			" Tree steps number",
			" ATM swaptions flag",
			" MktDataManager",
			" Products to price",
			" Is a standard SFRM calibration"
	},
	{		
			" Local_CRACalculator_CreateFromPf",	/// name of the C++ function
            " RRRRRRRRRR",	    /// 10 parametres = 9 d'entree + 1 de retour 
            " ARM_GC_CRACalculator_CreateFromPf",
			" OptionPf,[MRSBeta],[calibSecPFParams],[nbSteps],[flagToGenerateOSWATM],MktDataManager,ProductsToPrice, LocalModelParams,[IsSfrmStdCalib]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator from a given option portfolio",
			" Option portfolio object",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, [Recalibrate MRS flag, Recalibrate Beta flag] (def: taken from portfolio)",
			" Calibration parameters (def: taken from portfolio)",
			" Tree steps number (def: taken from portfolio)",
			" ATM swaption flag (def: taken from portfolio)",
			" MktDataManager",
			" Products to price",
			" Local Model Parameters",
			" Is a standard SFRM calibration"
   },
   {
			" Local_PXL_CRACalculator_CreateFromPf",	/// name of the C++ function
            " RRRRRRRRRR",	    /// 10 parametres = 9 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_CreateFromPf",
			" OptionPf,[MRSBeta],[calibSecPFParams],[nbSteps],[flagToGenerateOSWATM],MktDataManager,ProductsToPrice, LocalModelParams,[IsSfrmStdCalib]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator from a given option portfolio",
			" Option portfolio object",
			" Vector : MRS Min, MRS Max, Beta Min, Beta Max, [Recalibrate MRS flag, Recalibrate Beta flag] (def: taken from portfolio)",
			" Calibration parameters (def: taken from portfolio)",
			" Tree steps number (def: taken from portfolio)",
			" ATM swaption flag (def: taken from portfolio)",
			" MktDataManager",
			" Products to price",
			" Local Model Parameters",
			" Is a standard SFRM calibration"
	},
    {
        	" Local_CRACalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CRACalculator_GetData",
			" CRAObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRA Calculator",
			" CRA Id",
            " Type = SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_PXL_CRACalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_GetData",
			" CRAObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRA Calculator",
			" CRA Id",
            " Type = SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",

    },
    {
        	" Local_CRACalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_CRACalculator_SetData",
			" CRAObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRA Calculator",
			" CRA Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
    {
        	" Local_PXL_CRACalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_SetData",
			" CRAObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRA Calculator",
			" CRA Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
	{
        	" Local_LocalCRACalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_LocalCRACalculator_GetData",
			" LocalCRAObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Local CRA Calculator",
			" Local CRA Id",
            " Type = LOCALPORTFOLIO, SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_PXL_LocalCRACalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_LocalCRACalculator_GetData",
			" LocalCRAObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a local CRA Calculator",
			" Local CRA Id",
            " Type = LOCALPORTFOLIO, SECURITY, CAP, FLOOR, FUNDING, STDLEG, RFLEG, STDSWAP, RFSWAP, BERMUDA, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
    },
	{
        	" Local_LocalCRACalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_LocalCRACalculator_SetData",
			" LocalCRAObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a local CRA Calculator",
			" Local CRA Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = LOCALPORTFOLIO, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
    {
        	" Local_PXL_LocalCRACalculator_SetData",	/// name of the C++ function
            " RRRRR",	                        /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_LocalCRACalculator_SetData",
			" LocalCRAObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a local CRA Calculator",
			" Local CRA Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = LOCALPORTFOLIO, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
	{
        	" Local_CRACalculator_SetRecalibFlags",	/// name of the C++ function
            " RRRR",	                    /// 4 parametres = 3 d'entree + 1 de retour 
            " ARM_GC_CRACalculator_SetRecalibFlags",
			" CRAObject, MRS flag, Beta flag",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Set new calibration flags (MRS ans Beta) to a CRA Calculator",
			" Local CRA Id",
            " recalibrate = 1 / do not recalibrate = 0",
            " recalibrate = 1 / do not recalibrate = 0"
    },
    {
        	" Local_PXL_CRACalculator_SetRecalibFlags",	/// name of the C++ function
            " RRRR",	                        /// 4 parametres = 3 d'entree + 1 de retour 
            " PXL_ARM_GC_CRACalculator_SetRecalibFlags",
			" CRAObject, MRS flag, Beta flag",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Set new calibration flags (MRS ans Beta) to a CRA Calculator",
			" Local CRA Id",
            " recalibrate = 1 / do not recalibrate = 0",
            " recalibrate = 1 / do not recalibrate = 0"
    },
	{		" Local_CRASpreadCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_GC_CRASpreadCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Ref Coeff 1",
			" Ref Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
   },
   {
        	" Local_PXL_CRASpreadCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " PXL_ARM_GC_CRASpreadCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" General Datas", 
		   	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Ref Coeff 1",
			" Ref Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
	},
	{		" Local_CRAVMSCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_GC_CRAVMSCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Ref Coeff 1",
			" Ref Tenor 1",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
   },
   {
        	" Local_PXL_CRAVMSCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " PXL_ARM_GC_CRAVMSCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefTenor1,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" General Datas", 
		   	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Ref Coeff 1",
			" Ref Tenor 1",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
	},
	{		
			" Local_CRASpreadCalculator_CreateFromPf",	/// name of the C++ function
            " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 de retour 
            " ARM_GC_CRASpreadCalculator_CreateFromPf",
			" OptionPf,ModelData,payIndexMult,MktDataManager,ProductsToPrice,[refIndexResetFreq]",
            " 1",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Spread Calculator (not a CCSO!) from a given option portfolio",
			" Option portfolio object",
            " vector of strings {model, calib type, calib strike type, vol basket}",
			" payIndexMult",
			" MktDataManager",
			" Products to price",
			" modify reset frequency of corridor's reference index (if empty, keep original frequency)"
	},
	{		
            " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 de retour 
            " ARM_GC_CRASpreadCalculator_CreateFromPf",
			" OptionPf,ModelData,payIndexMult,MktDataManager,ProductsToPrice,[refIndexResetFreq]",
            " 0",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Spread Calculator (not a CCSO!) from a given option portfolio",
			" Option portfolio object",
            " vector of strings {model, calib type, calib strike type, vol basket}",
			" payIndexMult",
			" MktDataManager",
			" Products to price",
			" modify reset frequency of corridor's reference index (if empty, keep original frequency)"
	},
	{		
			" Local_BasicCRASpreadCalculator_CreateFromPf",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_BasicCRASpreadCalculator_CreateFromPf",
			" AsOfDate,OptionPf,payIndexMult,[refIndexResetFreq]",
            " 1",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CCSO degenerated into CRA (without market data)",
			" as of date",
			" Option portfolio object",
			" payIndexMult",
			" modify reset frequency of corridor's reference index (if empty, keep original frequency)"
	},
	{		
			" Local_PXL_BasicCRASpreadCalculator_CreateFromPf",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_BasicCRASpreadCalculator_CreateFromPf",
			" AsOfDate,OptionPf,payIndexMult,[refIndexResetFreq]",
            " 0",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CCSO degenerated into CRA (without market data)",
			" as of date",
			" Option portfolio object",
			" payIndexMult",
			" modify reset frequency of corridor's reference index (if empty, keep original frequency)"
	},
	{		
			" Local_BasicCRASpreadCalculator_CreateFromSwaption",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_BasicCRASpreadCalculator_CreateFromSwaption",
			" AsOfDate,Swaption,payIndexMult,[refIndexResetFreq]",
            " 1",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a degenerated CCSO from Swaption (without market data)",
			" as of date",
			" Swaption object",
			" payIndexMult",
			" modify reset frequency of reference index (if empty, keep original frequency)"
	},
	{		
			" Local_PXL_BasicCRASpreadCalculator_CreateFromSwaption",	/// name of the C++ function
            " RRRRR",							/// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_BasicCRASpreadCalculator_CreateFromSwaption",
			" AsOfDate,Swaption,payIndexMult,[refIndexResetFreq]",
            " 0",										/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a degenerated CCSO from Swaption (without market data)",
			" as of date",
			" Swaption object",
			" payIndexMult",
			" modify reset frequency of reference index (if empty, keep original frequency)"
	},
	{		" Local_CRASpreadCalculator_CreateWithoutMarketData",	/// name of the C++ function
            " RRRRRRRRRRRRRR",										/// 14 parametres = 13 d'entree + 1 de retour 
            " ARM_GC_CRASpreadCalculator_CreateWithoutMarketData",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Ref Coeff 1",
			" Ref Coeff 2"
	},
	{
        	" Local_ARM_INITCRASPREAD",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRRR",	/// 19 inputs + 1 output
            " ARM_INITCRASPREAD",
			" Calculator,Zc,CapVol,[CapSabr],SwoptVol,[SwoptSabr],[SigmaOrAlpha],[ConvAdjVolCap],[ConvAdjVolSwopt],[ConvAdjType],[CorrelCorr],[CorrelCap],[CorrelSwopt],Mrs,VolRatio,MrsSpread,Correl,ModelParams,[ProductsToPrice],[LocalCalibFlags]",
            " 1",					/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialise a CRA Calculator with Market Data and model parameters",
			" CRA Spread object Id",
            " Zero curve object",
			" Cap volatility object",
			" Swaption volatility object",
			" Vector of objects for Rho, Nu and Beta Cap volatilities (default: NULL)",
			" Vector of objects for Rho, Nu and Beta Swaption volatilities (default: NULL)",
			" SABR Vol type : S (default), SIGMA, A, ALPHA",
			" Convexity adjustment volatility (cap) object",
			" Convexity adjustment volatility (swaption) object",
			" string (default: SUMEXP)",
			" Correlation manager (CORR correlation) object",
			" Correlation cube (index-index IRG) object",
			" Correlation cube (index-index SWOPT) object",
			" MRS value",
			" VolRatio value",
			" MRS spread value",
			" Correl value",
            " vector of strings {model, calib type, calib strike type, vol basket}",
            " vector of strings: default {CSO(N),Funding(N),Fix(N),Floor(N),Cap(N)}",
			" vector of strings: default {N,N,N,N,N}"
    },
	{
        	" Local_PXL_ARM_INITCRASPREAD",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	/// 19 inputs + 1 output
            " PXL_ARM_INITCRASPREAD",
			" Calculator,Zc,CapVol,[CapSabr],SwoptVol,[SwoptSabr],[SigmaOrAlpha],[ConvAdjVolCap],[ConvAdjVolSwopt],[ConvAdjType],[CorrelCorr],[CorrelCap],[CorrelSwopt],Mrs,VolRatio,MrsSpread,Correl,ModelParams,[ProductsToPrice]",
            " 0",					/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Initialise a CRA Calculator with Market Data and model parameters"
    },
	{
        	" Local_CRASpreadCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CRASpreadCalculator_GetData",
			" CRASpreadObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRA Calculator",
			" CRA Spread Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",
    },
    {
        	" Local_PXL_CRASpreadCalculator_GetData",	/// name of the C++ function
            " RRR",	                            /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CRASpreadCalculator_GetData",
			" CRAObject, TypeToGet",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a CRA Calculator",
			" CRA Spread Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO",

    },
    {
        	" Local_CRASpreadCalculator_SetData",	/// name of the C++ function
            " RRRRRR",	                    /// 6 parametres = 5 d'entree + 1 de retour 
            " ARM_GC_CRASpeadCalculator_SetData",
			" CRASpreadObject,Object,[PortfolioType],[MktDataKeys],[Update]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRA Calculator",
			" CRA Spread Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys",
			" Update Calculator (Y/N)"
    },
    {
        	" Local_PXL_CRASpreadCalculator_SetData",	/// name of the C++ function
            " RRRRRR",	                    /// 6 parametres = 5 d'entree + 1 de retour 
            " PXL_ARM_GC_CRASpeadCalculator_SetData",
			" CRAObject,Object,[PortfolioType],[MktDataKeys],[Update]",
            " 0",						        /// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a CRA Calculator",
			" CRA Spread Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys",
			" Update Calculator (Y/N)"
    },
	{
        	" Local_ARM_GetCcyFromGenCalculator",	/// name of the C++ function
            " RRR",	                                /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_GetCcyNameFromGenCalculator",
			" Gen Calculator Object,Ccy type",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Return currency name of a generic calculator",
			" string:GenCalculator Id",
            " string: Currency type (DOMESTIC,FOREIGN,FUNDING,COUPON)"
    },
	{		" Local_GlobalCapCalculator_CreateFromPf",	/// name of the C++ function
            " RRRRRRRR",	    /// 8 parametres = 7 d'entree + 1 de retour 
            " ARM_GC_GlobalCapCalculator_CreateFromPf",
			" GlobalCap,FundLev,CapLev,ModelParams,CalibParams,MktDataManager,ProductsToPrice",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Global Cap Calculator",
			" Global Cap : a cap/floor security", 
			" Fund Lev : Double or Curve",
			" Cap Lev  : Double or Curve",
			" Model Params : NbSteps, GeneratorType, InversionMethod, Antithetic, SampleType",
			" Calib Params : TimeStepNb, SpaceStepNb, StdDevNb, HkVolStart, HkVolAdd",
			" MktDataManager",
			" Products to price : Product, Swap, Coupon, Fund, GlobalCap"
	},
	{		" Local_PXL_GlobalCapCalculator_CreateFromPf",	/// name of the C++ function
            " RRRRRRRR",	    /// 8 parametres = 7 d'entree + 1 de retour 
            " PXL_ARM_GC_GlobalCapCalculator_CreateFromPf",
			" GlobalCap,FundLev,CapLev,ModelParams,CalibParams,MktDataManager,ProductsToPrice",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Global Cap Calculator",
			" Global Cap : a cap/floor security", 
			" Fund Lev : Double or Curve",
			" Cap Lev  : Double or Curve",
			" Model Params : NbSteps, GeneratorType, InversionMethod, Antithetic, SampleType",
			" Calib Params : TimeStepNb, SpaceStepNb, StdDevNb, HkVolStart, HkVolAdd",
			" MktDataManager",
			" Products to price : Product, Swap, Coupon, Fund, GlobalCap"
	},
	{		" Local_GlobalCapCalculator_CreateFromPfWithoutMktData",	/// name of the C++ function
            " RRRRR",	    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_GlobalCapCalculator_CreateFromPfWithoutMktData",
			" AsOfDate, GlobalCap,FundLev,CapLev",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Global Cap Calculator from Summit Securities",
			" As Of Date : Pricing date",
			" GlobalCap :  a cap/floor security", 
			" Fund Lev : Double or Curve",
			" Cap Lev  : Double or Curve"
	},
	{		" Local_GlobalCapCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRR",	    /// 15 parametres = 14 d'entree + 1 de retour 
            " ARM_GC_GlobalCapCalculator_Create",
			" General,Notional,Fund,FundLev,GlobalCapParams,[PastFixings],CapLev,CapFixed,CapStrike,CapSpread,ModelParams,CalibParams,MktDataManager,ProductsToPrice",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Global Cap Calculator",
			" GeneralData : Ccy, StartDate, EndDate, Receive/Pay", 
           	" Notional : Double or Curve",
			" Fund Data : IndexType, DayCount, Frequency, ResetGap, PayGap, PayTiming, PayCalendar, AdjRule, IntRule",
			" Fund Lev : Double or Curve",
			" Global Cap Params : Nb Period, Final Ratio, Global Cap Value",
			" Past Fixings : Vector of past fixings (in %)",
			" Cap Lev  : Double or Curve",
			" Cap Fixed  : Double or Curve",
			" Cap Strike  : Double or Curve",
			" Cap Spread  : Double or Curve",
			" Model Params : NbSteps, GeneratorType, InversionMethod, TransformAlgo[, GeneratorType, InversionMethod, TransformAlgo], SampleType",
			" Calib Params : TimeStepNb, SpaceStepNb, StdDevNb, HkVolStart, HkVolAdd",
			" MktDataManager",
			" Products to price : Product, Swap, Coupon, Fund, GlobalCap"
	},
	{		" Local_PXL_GlobalCapCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRR",	    /// 15 parametres = 14 d'entree + 1 de retour 
            " PXL_ARM_GC_GlobalCapCalculator_Create",
			" General,Notional,Fund,FundLev,GlobalCapParams,[PastFixings],CapLev,CapFixed,CapStrike,CapSpread,ModelParams,CalibParams,MktDataManager,ProductsToPrice",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Global Cap Calculator",
			" GeneralData : Ccy, StartDate, EndDate, Receive/Pay", 
           	" Notional : Double or Curve",
			" Fund Data : IndexType, DayCount, Frequency, ResetGap, PayGap, PayTiming, PayCalendar, AdjRule, IntRule",
			" Fund Lev : Double or Curve",
			" Global Cap Params : Nb Period, Final Ratio, Global Cap Value",
			" Past Fixings : Vector of past fixings (in %)",
			" Cap Lev  : Double or Curve",
			" Cap Fixed  : Double or Curve",
			" Cap Strike  : Double or Curve",
			" Cap Spread  : Double or Curve",
			" Model Params : NbSteps, GeneratorType, InversionMethod, TransformAlgo[, GeneratorType, InversionMethod, TransformAlgo], SampleType",
			" Calib Params : TimeStepNb, SpaceStepNb, StdDevNb, HkVolStart, HkVolAdd",
			" MktDataManager",
			" Products to price : Product, Swap, Coupon, Fund, GlobalCap"
	},
	{
        	" Local_GlobalCapCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_GlobalCapCalculator_SetData",
			" GlobalCapObject,Object,[PortfolioType],[MktDataKeys]",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Global Cap Calculator",
			" Global Cap Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
	{
        	" Local_PXL_GlobalCapCalculator_SetData",	/// name of the C++ function
            " RRRRR",	                    /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_GlobalCapCalculator_SetData",
			" GlobalCapObject,Object,[PortfolioType],[MktDataKeys]",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Force internal objects of a Global Cap Calculator",
			" Global Cap Id",
            " Object Id = a GenSecurity, a PricingModel, a CalibMethod, a MDM or a Portfolio",
            " For portfolio, type to set = OSWPORTFOLIO, STMPORTFOLIO, CFPORTFOLIO, OSWSECONDPORTFOLIO",
			" Market Data Keys"
    },
	{
			" Local_InitGlobalCap",		// name of the C++ function
			" RRRRRRRRRRR",				// 11 parametres = 10 d'entree + 1 de retour
			" ARM_INIT_GLOBALCAP",
			" GCId, ZcId, CapVolId, RhoCapId, NuCapId, BetaCapId, Mrs, CalibParams, ModelParams, ProductsToPrice",
			" 1",
			XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
			" Initialise a Global Cap calculator with Market Data",
			" Global Cap Calculator Id",
			" Zero Curve Id",
			" Cap Volatility Curve Id",
			" Rho Cap Id",
			" Nu Cap Id",
			" Beta Cap Id",
			" Mean Rerversion",
			" CalibParams : TimeStepNb, SpaceStepNb, StdDevNb, HkVolStart, HkVolAdd"
			" Model Params : NbSteps, GeneratorType, InversionMethod, TransformAlgo[, GeneratorType, InversionMethod, TransformAlgo], SampleType",
			" Products to price : Product, Swap, Coupon, Fund, GlobalCap"
	},
	{
        	" Local_TarnSetOutputFlags",	/// name of the C++ function
            " RRR",	                            /// 2 parametres = 1 d'entree + 1 de retour 
            " ARM_TarnSetOutputFlags",
			" Tarn Calculator string Id, Flag Vector",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Update the calculated outputs for a TARN security",
			" Calculator Id"
    },
	{
        	" Local_CallOnMepiVanillaArg_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",                      /// 16 parametres = 15 d'entree + 1 de retour 
            " ARM_GP_VanillaCallOnMepiCreate",
			" CurveName,EquityName,StartDate,EndDate,ResetFreq,RiskFactor,Strike,MaxBorrow,ProtectionCurveStart,ProtectionCurveEnd,StartingPf,StartingCash,MinInvested,LvgCost,CashSpread,Fees, AlreadyAsianed, AsianDatesNb",
            " 1 ",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Call On Mepi Vanilla Arg ",
			" CurveName",
			" EquityName", 
			" StartDate",
			" EndDate",
			" ResetFreq",
			" RiskFactor",
			" Strike",
			" MaxBorrow",
			" ProtectionCurveStart", 
			" ProtectionCurveEnd", 
			" StartingPf",
			" StartingCash", 
			" MinInvested",
			" LvgCost", 
			" CashSpread", 
			" Fees",
			" AlreadyAsianed",
			" AsianDatesNb"
    },
    {
        	" Local_CSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",			/// 17 parametres = 16 d'entree + 1 de retour 
            " ARM_GC_CSOCalculator_Create",
			" Start, End, MktDatas, CpnDatas, FundDatas, ExerDatas, Notional, CpnMin, CpnMax, Leverage, FixCpnCurve, FundSpread, Fees, calibFlags, ProductFlags, ModelDatas",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
			" Start date",
            " End date",
            " Market Datas",
            " Coupon Datas",
            " Funding datas = Freq (cpn freq), DC (A360), FundMargin, FundCoeff",
            " ExerDatas",
            " Notional",
            " CpnMin",
            " CpnMax",
			" Leverage",
			" Fix Coupon Curve",
            " Funding Spread",
            " fees",
            " calibFlags",
            " ProductFlags",
            " ModelDatas"
    },
    {
        	" Local_PXL_CSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",			/// 17 parametres = 16 d'entree + 1 de retour 
            " PXL_ARM_GC_CSOCalculator_Create",
			" Start, End, MktDatas, CpnDatas, FundDatas, ExerDatas, Notional, CpnMin, CpnMax, Leverage, FixCpnCurve, FundSpread, Fees, calibFlags, ProductFlags, ModelDatas",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
			" Start date",
            " End date",
            " Market Datas",
            " Coupon Datas",
            " Funding datas = Freq (cpn freq), DC (A360), FundMargin, FundCoeff",
            " ExerDatas",
            " Notional",
            " CpnMin",
            " CpnMax",
			" Leverage",
			" Fix Coupon Curve",
            " Funding Spread",
            " fees",
            " calibFlags",
            " ProductFlags",
            " ModelDatas"
    },
	{
        	" Local_PXL_ExtendedCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 20 parametres = 19 d'entree + 1 de retour 
            " PXL_ARM_GC_ExtendedCSOCalculator_Create",
            " Start,End,MktDatas,CpnDatas,FundDatas,ExerDatas,[CpnNotional],[CpnMin],[CpnMax],[LeverageLong],[LeverageShort],[Strike],[FixCpnCurve],[FundSpread],[Fees],[calibFlags],[ProductFlags],[ModelDatas],[FundLeverage]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
    },
	{
        	" Local_ExtendedCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 20 parametres = 19 d'entree + 1 de retour 
            " ARM_GC_ExtendedCSOCalculator_Create",
            " Start,End,MktDatas,CpnDatas,FundDatas,ExerDatas,[Notional],[CpnMin],[CpnMax],[LvgeLong],[LvgeShort],[Strike],FixCpnCurve,[FundSpread],[Fees],[CalibFlags],[ProductFlags],[ModelDatas],[FundLeverage]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
			" date: Start date",
            " date: End date",
            " vector of strings: Mkt dataManager and Mkt DataKeys",
            " vector of strings: Coupon Datas =freq, Long index, short index, day count, ResetTiming(K_ADVANCE)",
            " vector of strings: Funding datas = Freq (cpn freq), DC (A360), ResetTiming (K_ADVANCE)",
            " vector of strings: ExerDatas = freq, notification gap",
            " string: Cpn notional curve id  = default(100)",
            " string: CpnMin curve id = default(0.0)",
            " string: CpnMax curve id = default(1.0e+10)",
			" string: LeverageLong curve id = default(1.0)",
			" string: LeverageShort curve id = default(1.0)",
			" string: Strike curve id = default(0.0)",
			" string: Fix Coupon curve id",
            " string: Funding Spread curve id = default(0.0)", 
            " string: Fees curve id ",
            " vector of strings: CalibFlags = default (calib type(DIAG), strike type(ATM), model type(HWM1F)) ",
            " vector of strings: ProductFlags = default (CSO (N),Funding(N),Fix(N),Floor(N),Cap(N)) ",
            " vector of double : ModelDatas = default(vector(1,0))",
			" string: Funding Leverage id (default = flat curve = 1.0)",
    },
		/// it is ugly FIX FIX FIX, we have to merge all interface.
	{
        	" Local_PXL_BasisCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRR",			/// 18 parametres = 17 d'entree + 1 de retour 
            " PXL_ARM_GC_BasisCSOCalculator_Create",
            " Start, FixEnd, End, Ccys, MktDataManager,CpnDatas, FundDatas, ExerDatas, Notionals, Cpns, Leverages, Strike, [FundSpread], [Fees], [CalibDatas], [ProductFlags],[ModelDatas] ",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator bi currencies",
    },
	{
        	" Local_BasisCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRR",			/// 18 parametres = 17 d'entree + 1 de retour 
            " ARM_GC_BasisCSOCalculator_Create",
			" Start, FixEnd, End, Ccys, MktDataManager ,CpnDatas, FundDatas, ExerDatas, Notionals, Cpns, Leverages, Strike, [FundSpread], [Fees], [CalibDatas], [ProductFlags],[ModelDatas] ",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator bi currencies",
			" date: Start date",
            " date: End date of fix period",
			" date: End date",
			" vector of strings: Coupon Ccy, Funding Ccy= default(Coupon Currency)",
			" string: Mkt dataManager ",
            " vector of strings: Coupon Datas =freq, Long index, short index, DC, resetTiming(ADV),restCalendar(), payCalendar(), stubRule(SS)",
            " vector of strings: Funding datas = freq(cpn freq), DC (A360), resetTiming (ADV)",
            " vector of strings: ExerDatas = freq, notification gap, pay/rec(REC),nbNoCall (0)",
            " vector of strings: Cpn notional curve id, Fund notional = default(Cpn notional)",
            " vector of strings: CpnMin curve id,CpnMax curve id",
			" vector of strings: LeverageLong curve id, LeverageShort curve id",
			" string: Strike curve id ",
            " string: Funding Spread curve id = default(0.0)", 
            " string: Fees curve id = default(0.0)",
            " vector of strings: CalibDatas = default (calib type(DIAG), strike type(ATM), model type(HWM1F), vns method(MOYENESS), moyeness level(1.0)) ",
            " vector of strings: ProductFlags = default (CSO (N),Funding(N),Fix(N),Floor(N),Cap(N)) ",
            " vector of double : ModelDatas = default(0)",
    },
	{
        	" Local_BasicCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 input + 1 output 
            " ARM_GC_BasicCSOCalculator_Create",
            " AsOf,Start,End,CpnDatas,FundDatas,ExerDatas,[Notional],[CpnMin],[CpnMax],[LvgeLong],[LvgeShort],[Strike],FixCpnCurve,[FundSpread],[Fees],[FundLeverage],[CpnCcy],[FundCcy],[fixEndDate]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
			" date : as of date",
			" date: Start date",
            " date: End date",
            " vector of strings: Coupon Datas =freq, Long index, short index, day count, ResetTiming(K_ADVANCE)",
            " vector of strings: Funding datas = Freq (cpn freq), DC (A360), ResetTiming (K_ADVANCE)",
            " vector of strings: ExerDatas = freq, notification gap",
            " string: Cpn notional curve id  = default(100)",
            " string: CpnMin curve id = default(0.0)",
            " string: CpnMax curve id = default(1.0e+10)",
			" string: LeverageLong curve id = default(1.0)",
			" string: LeverageShort curve id = default(1.0)",
			" string: Strike curve id = default(0.0)",
			" string: Fix Coupon curve id",
            " string: Funding Spread curve id = default(0.0)", 
            " string: Fees curve id ",
			" string: Funding Leverage id (default = flat curve = 1.0)",
            " string: Coupon Currency",
            " string: Funding Currency",
            " date:   Fix End Date"
    },
	{
        	" Local_PXL_BasicCSOCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 input + 1 output 
            " PXL_ARM_GC_BasicCSOCalculator_Create",
            " AsOf,Start,End,CpnDatas,FundDatas,ExerDatas,[Notional],[CpnMin],[CpnMax],[LvgeLong],[LvgeShort],[Strike],FixCpnCurve,[FundSpread],[Fees],[FundLeverage],[CpnCcy],[FundCcy],[fixEndDate]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable SpreaOption Calculator",
			" date : as of date",
			" date: Start date",
            " date: End date",
            " vector of strings: Coupon Datas =freq, Long index, short index, day count, ResetTiming(K_ADVANCE)",
            " vector of strings: Funding datas = Freq (cpn freq), DC (A360), ResetTiming (K_ADVANCE)",
            " vector of strings: ExerDatas = freq, notification gap",
            " string: Cpn notional curve id  = default(100)",
            " string: CpnMin curve id = default(0.0)",
            " string: CpnMax curve id = default(1.0e+10)",
			" string: LeverageLong curve id = default(1.0)",
			" string: LeverageShort curve id = default(1.0)",
			" string: Strike curve id = default(0.0)",
			" string: Fix Coupon curve id",
            " string: Funding Spread curve id = default(0.0)", 
            " string: Fees curve id ",
			" string: Funding Leverage id (default = flat curve = 1.0)",
            " string: Coupon Currency",
            " string: Funding Currency",
            " date:   Fix End Date"
    },
	{
        	" Local_ARM_INITCSO",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	/// 20 inputs + 1 output
            " ARM_INITCSO",
            " CsoId,Zc,CapVol,SwoptVol,[SabrCap],[SabrSwopt],[SigmaOrAlpha],[FlatVol],[ConvAdjVol],[ConvAdjMgr],[CorrelCap],[CorrelSwopt],[CorrelCorr],CurveModelParams,[ModelData],CalibParams,ProductsToPrice,[fundZc],[domBasisZc],[fundBasisZc]",
            " 1",					/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with Market Data",
			" CSO object",
            " Zero curve object",
			" Cap volatility object",
			" Swaption volatility object",
			" Vector of objects : Rho, Nu and Beta volatilities (default: NULL)",
			" Vector of objects : Rho, Nu and Beta volatilities (default: NULL)",
			" SABR Vol type : S (default), SIGMA, A, ALPHA",
 			" Flat volatility object (for spread lock)",
			" Convexity adjustment volatility object",
			" Convexity adjustment manager object",
			" Correlation cube (index-index IRG) object",
			" Correlation cube (index-index SWOPT) object",
			" Correlation manager (CORR correlation) object",
			" MRS, Correl, VolRatio, MRS spread objects",
            " vector of double: default {0}",
            " vector of string: calib(DIAG),calib strike(ATM),model(HW1F),vns(MONEYNESS)",
            " vector of string: default {CSO(N),Funding(N),Fix(N),Floor(N),Cap(N)}",
			" Zero curve object (default: NULL)",
 			" Zero curve object (default: NULL)",
			" Zero curve object (default: NULL)"
   },
	{
        	" Local_PXL_ARM_INITCSO",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",	/// 20 inputs + 1 output
            " PXL_ARM_INITCSO",
            " CsoId,Zc,CapVol,SwoptVol,[SabrCap],[SabrSwopt],[SigmaOrAlpha],[FlatVol],[ConvAdjVol],[ConvAdjMgr],[CorrelCap],[CorrelSwopt],[CorrelCorr],CurveModelParams,[ModelData],CalibParams,ProductsToPrice,[fundZc],[domBasisZc],[fundBasisZc]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with Market Data",
    },
	{
        	" Local_ARM_INITCSO",		/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",		/// 19 inputs + 1 output
            " ARM_INITCSO",
            " CsoId,Zc,CapVol,SwoptVol,[SabrCap],[SabrSwopt],[SigmaOrAlpha],[FlatVol],[ConvAdjVol],[ConvAdjMgr],[CorrelCap],[CorrelSwopt],[CorrelCorr],CurveModelParams,[ModelData],CalibParams,ProductsToPrice,[Forex],[Fund&BasisCurves]",
            " 1",					/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with Market Data",
			" CSO object",
            " Zero curve object",
			" Cap volatility object",
			" Swaption volatility object",
			" Vector of objects : Rho, Nu and Beta volatilities (default: NULL)",
			" Vector of objects : Rho, Nu and Beta volatilities (default: NULL)",
			" SABR Vol type : S (default), SIGMA, A, ALPHA",
 			" Flat volatility object (for spread lock)",
			" Convexity adjustment volatility object",
			" Convexity adjustment manager object",
			" Correlation cube (index-index IRG) object",
			" Correlation cube (index-index SWOPT) object",
			" Correlation manager (CORR correlation) object",
			" MRS, Correl, VolRatio, MRS spread objects",
            " vector of double: default {0}",
            " vector of string: calib(DIAG),calib strike(ATM),model(HW1F),vns(MONEYNESS)",
            " vector of string: default {CSO(N),Funding(N),Fix(N),Floor(N),Cap(N)}",
			" Forex object (default: NULL)",
			" Vector of objects : FundZc, DomBasisZc, FundBasisZc (default: NULL)"
	},
	{
        	" Local_PXL_ARM_INITCSO",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",		/// 19 inputs + 1 output
            " PXL_ARM_INITCSO",
            " CsoId,Zc,CapVol,SwoptVol,[SabrCap],[SabrSwopt],[SigmaOrAlpha],[FlatVol],[ConvAdjVol],[ConvAdjMgr],[CorrelCap],[CorrelSwopt],[CorrelCorr],CurveModelParams,[ModelData],CalibParams,ProductsToPrice,[Forex],[Fund&BasisCurves]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with Market Data",
    },
	{
        	" Local_ARM_INITCSO_withMDM",	/// name of the C++ function
            " RRRRRR",				/// 5 inputs + 1 output
            " ARM_INITCSO_withMDM",
            " CsoId, MktData, ModelParams, CalibParams, ProductsToPrice",
            " 1",					/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with a Market Data Manager",
			" CSO object Id",
            " vector of strings: MktDataManager and MktDataKeys",
            " vector of double: default {0}",
            " vector of strings: calib type(def DIAG), calib strike type(def ATM), model type(def HW1F), vns method(def MONEYNESS)",
            " vector of strings: default {CSO(N),Funding(N),Fix(N),Floor(N),Cap(N)} ",
    },
	{
        	" Local_PXL_ARM_INITCSO_withMDM",	/// name of the C++ function
            " RRRRRR",					/// 5 inputs + 1 output
            " PXL_ARM_INITCSO_withMDM",
            " CsoId, MktData, ModelParams, CalibParams, ProductsToPrice",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a CSO calculator and initialise it with a Market Data Manager",
    },
	{
        	" Local_CSOCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GC_CSOCalculator_GetData",
			" CSOObject, TypeToGet",
            " 1",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Callable SpreadOption Calculator",
			" Callable SpreadOption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO,CSOCSTMANAGER",
    },
	{
        	" Local_PXL_CSOCalculator_GetData",	/// name of the C++ function
            " RRR",	                        /// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GC_CSOCalculator_GetData",
			" CSOObject, TypeToGet",
            " 0",						    /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Get internal objects from a Callable SpreadOption Calculator",
			" Callable SpreadOption Id",
            " Type = SECURITY, MODEL, CALIBMETHOD, MDM, OSWPORTFOLIO, CFPORTFOLIO",
    },
	{		" Local_SnowRangeCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",	    /// 15 inputs + 1 output
            " ARM_GC_SnowRangeCalculator_Create",
			" General,Notional,FundingData,CouponData,Spread,Strike,Ratchet,CashFlow,FixedRate,Leverage,SnowRangeParams,CalibParams,ModelParams,MCParams,MktDataManager,ProductsToPrice",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Snow Range Calculator",
			" GeneralData : Ccy, StartDate, EndDate, Receive/Pay", 
           	" Notional : Double or Curve",
			" Funding Data : Index Term, DayCount",
			" Coupon Data : Index Term, DayCount, ResetFreq, PayFreq, ResetTiming, PayTiming, ResetCal, PayCal, AdjRule, IntRule",
			" Spread : Double or Curve",
			" Strike  : Double or Curve",
			" Ratchet  : Double or Curve",
			" Cash Flow  : Double or Curve",
			" Fixed Rate  : Double or Curve",
			" Leverage  : Double or Curve",
			" Snow Range Params : BarrierShift (%), IsCallable (YES/NO), UpFront (%)",
			" Calibration Params : Time Step Nb, Space Step Nb, Std Dev Nb, Calib Swaption (YES/NO)",
			" Model Params : Model Name + corresponding model parameters",
			" Monte Carlo Params : NbSteps, GeneratorType, InversionMethod, Antithetic, SampleType",
			" MktDataManager",
			" Products to price : Funding, Coupon, StdCoupon, PlainCoupon, Swap, CallOption, CallSwap"
	},
	{		" Local_PXL_SnowRangeCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRR",	    /// 15 inputs + 1 output
            " PXL_ARM_GC_SnowRangeCalculator_Create",
			" General,Notional,FundingData,CouponData,Spread,Strike,Ratchet,CashFlow,FixedRate,Leverage,SnowRangeParams,CalibParams,ModelParams,MCParams,MktDataManager,ProductsToPrice",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Snow Range Calculator",
			" GeneralData : Ccy, StartDate, EndDate, Receive/Pay", 
           	" Notional : Double or Curve",
			" Funding Data : Index Term, DayCount",
			" Coupon Data : Index Term, DayCount, ResetFreq, PayFreq, ResetTiming, PayTiming, ResetCal, PayCal, AdjRule, IntRule",
			" Spread : Double or Curve",
			" Strike  : Double or Curve",
			" Ratchet  : Double or Curve",
			" Cash Flow  : Double or Curve",
			" Fixed Rate  : Double or Curve",
			" Leverage  : Double or Curve",
			" Snow Range Params : BarrierShift (%), IsCallable (YES/NO), UpFront (%)",
			" Calibration Params : Time Step Nb, Space Step Nb, Std Dev Nb, Calib Swaption (YES/NO)",
			" Model Params : Model Name + corresponding model parameters",
			" Monte Carlo Params : NbSteps, GeneratorType, InversionMethod, Antithetic, SampleType",
			" MktDataManager",
			" Products to price : Funding, Coupon, StdCoupon, PlainCoupon, Swap, CallOption, CallSwap"
	},
	{		" Local_CRADoubleCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",			/// 21 parametres = 20 d'entree + 1 de retour 
            " ARM_GC_CRADoubleCalculator_Create",
			" General,Call,Fund,Cpn,Notio,CallFees,FundData,BoostFix,PayIdxMult,RBarDown,RBarUp,SpdBarDown,SpdBarUp,SpdCoef1,SpdCoef2,ModDatas,MktDatMngr,ProdToPrice,[LocCalFlag],[MiscData]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Double Condition Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Datas : [0] = Spread Curve, [1] = Fund Lev",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Rate Condition Barrier Down Curve",
			" Rate Condition Barrier Up Curve",
			" Spread Barrier Down Curve",
			" Spread Barrier Up Curve",
			" Spread Coeff 1",
			" Spread Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
   },
	{		" Local_PXL_CRADoubleCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRRR",			/// 21 parametres = 20 d'entree + 1 de retour 
            " PXL_ARM_GC_CRADoubleCalculator_Create",
			" General,Call,Fund,Cpn,Notio,CallFees,FundData,BoostFix,PayIdxMult,RBarDown,RBarUp,SpdBarDown,SpdBarUp,SpdCoef1,SpdCoef2,ModDatas,MktDatMngr,ProdToPrice,[LocCalFlag],[MiscData]",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Double Condition Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Datas : [0] = Spread Curve, [1] = Fund Lev",
			" Boosted Fix Curve",
			" PayIndexMult Curve",
			" Rate Condition Barrier Down Curve",
			" Rate Condition Barrier Up Curve",
			" Spread Barrier Down Curve",
			" Spread Barrier Up Curve",
			" Spread Coeff 1",
			" Spread Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
   },
	{		" Local_ConvertToVarNotionalSwaption",	/// name of the C++ function
            " RRR",									/// 3 parametres = 2 d'entree + 1 de retour 
            " ARM_GP_ConvertToVarNotionalSwaption",
			" YieldCurveId, SwaptionId",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Convert a general swaption to a variable notional swaption",
			" Yield curve", 
           	" Swaption to convert"
   },
	{		" Local_PXL_ConvertToVarNotionalSwaption",	/// name of the C++ function
            " RRR",										/// 3 parametres = 2 d'entree + 1 de retour 
            " PXL_ARM_GP_ConvertToVarNotionalSwaption",
			" YieldCurveId, SwaptionId",
            " 0",										/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Convert a general swaption to a variable notional swaption",
			" Yield curve", 
           	" Swaption to convert"
   },
	{		" Local_CMRASpreadCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_GC_CMRASpreadCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Triple Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curves",
			" PayIndexMult Curve",
			" Barrier Down Curves",
			" Barrier Up Curves",
			" Ref Coeff 1",
			" Ref Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
   },
   {
        	" Local_PXL_CMRASpreadCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " PXL_ARM_GC_CMRASpreadCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,PayIndexMult,BarrierDown,BarrierUp,RefCoeff1,RefCoeff2,ModelDatas,MktDataManager,ProdToPrice,[LocCalFlag],[MiscData]",
            " 0",						/// not visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Triple Range Accrual Calculator",
			" General Datas", 
		   	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Cpn Spread Curve",
			" Boosted Fix Curves",
			" PayIndexMult Curve",
			" Barrier Down Curves",
			" Barrier Up Curves",
			" Ref Coeff 1",
			" Ref Coeff 2",
			" ModelDatas",
			" MktDataManager",
			" Products to price",
			" Calib Flags",
			" ExerProba:[0]=y(1)/n(0),[1]=prIdx, OptimData:[2]=y(1)/n(0),[3]=DLim,[4]=WLim,[5]=MLim"
	},
	{
        	" Local_BasisConverter",	/// name of the C++ function
            " RRRRRRRRRRRRRRRRRRR",				/// 19 parametres = 18 d'entree + 1 de retour 
            " ARM_GP_BasisConverter",
 			" AsOfDate,DomCcy,ForCcy,DomDateStrip,ForDateStrip,FundDateStrip,DomDayCount,DomFreq,ForDayCount,ForFreq,DomCurve,ForCurve,DomDiscCurve,ForDiscCurve,Forex,DomNotional,ForNotional,ForSpread",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Convert the funding spread by taking into account the basis.",
			" AsOfDate",
            " Domestic Currency",
			" Foreign Currency",
			" Domestic DateStrip",
			" Foreign DateStrip",
			" Funding DateStrip",
            " Domestic DayCount",
			" Domestic Freq",
			" Foreign DayCount",
			" Foreign Freq",
			" Domestic Zero Curve",
            " Foreign Zero Curve",
			" Domestic Zero Discount Curve",
            " Foreign Zero Discount Curve",
			" Forex",
			" Domestic Notional",
			" Foreign Notional",
			" Foreign Spread"
    },
	{
        	" Local_GenericAddin",	/// name of the C++ function
            " RRRR",				/// 4 parametres = 3 d'entree + 1 de retour 
            " ARM_GA",
 			" FunctionName,ParamNames,ParamValues",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Call the generic addin.",
			" AsOfDate",
            " Function Name",
			" Param Names",
			" Param Values"
    },
	{
        	" PXL_Local_GenericAddin",	/// name of the C++ function
            " RRRR",				/// 4 parametres = 3 d'entree + 1 de retour 
            " PXL_ARM_GA",
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
            " ARM_GA_H",
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
            " PXL_ARM_GA_H",
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
            " ARM_GA_P",
 			" FunctionName, WithDefault",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create the generic addin helper.",
            " Function Name",
			" With Default Parameters (N)",
    },
	{		" Local_CRAQuantoCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRR",			/// 12 parametres = 11 d'entree + 1 de retour 
            " ARM_GC_CRAQuantoCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,BarrierDown,BarrierUp,ProdToPrice",
            " 1",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Products to price"
   },
   {
        	" Local_PXL_CRAQuantoCalculator_Create",	/// name of the C++ function
            " RRRRRRRRRRRR",			/// 19 parametres = 18 d'entree + 1 de retour 
            " PXL_ARM_GC_CRAQuantoCalculator_Create",
			" General,Call,Fund,Cpn,Notional,CallFees,FundSpread,BoostedFix,BarrierDown,BarrierUp,ProdToPrice",
            " 0",						/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a Callable Range Accrual Calculator",
			" GeneralDatas", 
           	" Call Datas",
			" Fund Datas",
			" Cpn Datas",
			" Notional Curve",
			" CallDate / Fees Curve",
			" Fund Spread Curve",
			" Boosted Fix Curve",
			" Barrier Down Curve",
			" Barrier Up Curve",
			" Products to price"
	},
    {
        	" Local_FXVanillaCalculator_CreateFromSecurity",	/// name of the C++ function
            " RRRRR",											/// 5 parametres = 4 d'entree + 1 de retour 
            " ARM_GC_FXVanillaCalculator_CreateFromSecurity",
			" Security,[BasketType],[DigitType],[VanillaType]",
            " 1",												/// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a FXVanilla Calculator from a security",
			" Security",
            " Basket Type (default: max)",
            " Digital Type (default: analytic)",
            " Vanilla Type (default: spread)",
    },
    {
        	" Local_PXL_FXVanillaCalculator_CreateFromSecurity", /// name of the C++ function
            " RRRRR",											 /// 5 parametres = 4 d'entree + 1 de retour 
            " PXL_ARM_GC_FXVanillaCalculator_CreateFromSecurity",
			" Security,[BasketType],[DigitType],[VanillaType]",
            " 0",												 /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a FXVanilla Calculator from a security",
			" Security",
            " Basket Type (default: max)",
            " Digital Type (default: analytic)",
            " Vanilla Type (default: spread)",
    },
	{
        	" Local_VolBondCalculator_Create",					 /// name of the C++ function
            " RRRRRRRRRRRRRRRRRRRR",						 /// 18 parametres = 17 d'entree + 1 de retour 
            " ARM_GC_VolBondCalculator_Create",
			" General,Leverage,PayFreq,ResetFreq,DayCount,Tenor,IntRule,StubRule,ResetGap,PayCalendar,ResetCalendar,OdeSolvers,RKParameters,MCParameters,RandomGenerator,PayOffType,MktDataMgr,MktDataMgrKeys,ProdToPrice",
            " 1",												 /// visible in excel
            XLLOCALARM_GENCALCULATOR_GROUP,
            " ",
            " ",
            " Create a VolBond Calculator",
			" a",
			" a",
			" a",
			" a",
			" a",
			" a",
			" a",
			" a",			
			" a",
			" a",			
    },	
/// END CALCULATOR

     
