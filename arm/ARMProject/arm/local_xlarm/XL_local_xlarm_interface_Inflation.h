	{
        		" Local_CreateInfCurv",		/// name of the C++ function
                " RRRRRRRRRRRRRR",			/// 14 parametres = 13 d'entree + 1 parametre de retour beware to include also Returns parameter
                " ARM_INF_CreateCurv",
				" asOfDate,indexName,CPIRefVal,CPIRefDate,maturities,values,[MonthlyInterp],[DailyInterp],[DCFMonthly],[DCFDaily],[ExtrapolType],[ResetManager],[SeasonManager]",
											/// first line in excel function wizard
                " 1",						/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a forward inflation curve",
											/// description in excel in the wizard
                " Date",				//// description of the various parameters
				" index Name like CPALEMU",	
				" CPI Index Value",
				" CPI Reference Date",
				" Array of Maturities",
                " Array of Values",
				" Default per index: CPILINEAR, ZCLINEAR, ZCCTFWD",
				" Default per index: STEPWISESTART, STEPWISEMID, STEPWISEEND, STEPWISE, CPILINEAR, ZCLINEAR, ZCCTFWD",
				" Default per index: ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE",
				" Default per index: ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE",
				" Default per index: LASTTWO, MIDDLE, FIRSTTWO",
				" ResetManager Object",
				" SeasonManager Object"
        },
		{
        		" Local_GetInfZcFromSummit_Create",		/// name of the C++ function
                " RRRRRRR",			/// 7 parametres = 6 d'entree + 1 parametre de retour beware to include also Returns parameter
                " ARM_INF_GetZcFromSummit",
				" indexName,ccy,cvName,asOfDate,[seasonAdj],[seasonAdjMode]",
											/// first line in excel function wizard
                " 1",						/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a forward inflation curve from Summit Datas",
											/// description in excel in the wizard
                " index Name like EMU",//// description of the various parameters
				" currency",
				" Summit curve Name like MO",
				" As Of date",
				" Season Adjustment? (default : NO)",
				" Seasonality Correction Mode (default : PLUS)"
        },
		{
        		" Local_PXL_GetInfZcFromSummit_Create",		/// name of the C++ function
                " RRRRRRR",			/// 7 parametres = 6 d'entree + 1 parametre de retour beware to include also Returns parameter
                " PXL_ARM_INF_GetZcFromSummit",
				" indexName,ccy,cvName,asOfDate,[seasonAdj],[seasonAdjMode]",
											/// first line in excel function wizard
                " 0",						/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a forward inflation curve from Summit Datas",
											/// description in excel in the wizard
                " index Name like EMU",//// description of the various parameters
				" currency",
				" Summit curve Name like MO",
				" As Of date",
				" Season Adjustment? (default : NO)",
				" Seasonality Correction Mode (default : PLUS)"
        },
        {
        		" Local_PXL_CreateInfCurv",
                " RRRRRRRRRRRRRR",			/// 14 parametres = 13 d'entree + 1 parametre de retour beware to include also Returns parameter
                " PXL_ARM_INF_CreateCurv",
				" asOfDate,indexName,CPIRefVal,CPIRefDate,maturities,values,[MonthlyInterp],[DailyInterp],[DCFMonthly],[DCFDaily],[ExtrapolType],[ResetManager],[SeasonManager]", 
                " 0",						//// visible
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a forward inflation curve",
                " Date",				//// description of the various parameters
				" index Name like CPALEMU",	
				" CPI Index Value",
				" CPI Reference Date",
				" Array of Maturities",
                " Array of Values",
				" Default per index: CPILINEAR, ZCLINEAR, ZCCTFWD",
				" Default per index: STEPWISESTART, STEPWISEMID, STEPWISEEND, STEPWISE, CPILINEAR, ZCLINEAR, ZCCTFWD",
				" Default per index: ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE",
				" Default per index: ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE",
				" Default per index: LASTTWO, MIDDLE, FIRSTTWO",
				" ResetManager Object",
				" SeasonManager Object"
        },
        {
        		" Local_InfCurv_CPIInterp",
                " RRRRRRR",					///  7 parametres 6 d'entree + 1 parametre de retour
                " ARM_INF_InterpCPI",
				" Crvid,CPIDate,[DCFLag],[DailyInterp],[ResetLag],[weight]",
                " 1",						//// non visible
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_INF_InterpCPI,
                " Interpolates an inflation curve", 
				" Curve id", 
                " CPI Date", 
				" DCF Lag, default per index, like 3m",
				" Default from the curve, STEPWISESTART, STEPWISEMID, STEPWISEEND, STEPWISE, CPILINEAR, ZCLINEAR, ZCCTFWD",
                " CPI Lag, default per index, like 0m",
                " weight if using CPISTEPWISE, default 0, between 0 1"
        },
        {
        		" Local_InfCurv_ZCRateInterp",
                " RRRRRRR",					///  7 parametres 6 d'entree + 1 parametre de retour
                " ARM_INF_InterpZCRate",
				" Crvid,CPIDate,[DCFLag],[DailyInterp],[ResetLag],[weight]",
                " 1",						//// non visible
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_INF_InterpZCRate,
                " Interpolates an inflation curve", 
				" Curve id", 
                " CPI Date", 
				" DCF Lag, default per index, like 3m",
				" Default from the curve, STEPWISESTART, STEPWISEMID, STEPWISEEND, STEPWISE, CPILINEAR, ZCLINEAR, ZCCTFWD",
                " CPI Lag, default per index, like 0m",
                " weight if using CPISTEPWISE, default 0 between 0 1"
        },
		{
				" Local_ARM_LIVRETACURVE",
                " RRRRRRRRRR",					// 9 parametres d'entree + 1 parametre de retour
                " ARM_LivretACurve",
                " asOfDate,infCurvId, euribCurvId, [flagRounding], [resetmanager],[fixingLivretAId],[fixingEuribId],[monthForAugust],[monthForFebruary]",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a livretA curve",
                " As of date",
				" Inflation curve id ",
                " Euribor curve id",
                " rounded flag",
				" reset Manager for inflation curve",
				" Fixing taux Livret A id",
				" Fixing taux euribor id",
				" reset month for August LivretA Reset (default JUN)",
				" reset month for February LivretA Reset (default NOV)"
        },
		{
				" Local_ARM_LIVRETACURVE",
                " RRRRRRRRRR",					// 9 parametres d'entree + 1 parametre de retour
                " ARM_INF_LivretACurve",
                " asOfDate,infCurvId, euribCurvId, [flagRounding], [resetmanager],[fixingLivretAId],[fixingEuribId],[monthForAugust],[monthForFebruary]",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a livretA curve",
                " As of date",
				" Inflation curve id ",
                " Euribor curve id",
                " rounded flag",
				" reset Manager for inflation curve",
				" Fixing taux Livret A id",
				" Fixing taux euribor id",
				" reset month for August LivretA Reset (default JUN)",
				" reset month for February LivretA Reset (default NOV)"
        },
		{
				" Local_PXL_ARM_LIVRETACURVE",
                " RRRRRRRRRR",					// 9 parametres d'entree + 1 parametre de retour
                " PXL_ARM_LivretACurve",
                " asOfDate,infCurvId, euribCurvId, [flagRounding], [resetmanager],[fixingLivretAId],[fixingEuribId],[monthForAugust],[monthForFebruary]",
                " 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a livretA curve",
                " As of date",
				" Inflation curve id ",
                " Euribor curve id",
                " rounded flag",
				" reset Manager for inflation curve",
				" Fixing taux Livret A id",
				" Fixing taux euribor id",
				" reset month for August LivretA Reset (default JUN)",
				" reset month for February LivretA Reset (default NOV)"
        },
		{
				" Local_PXL_ARM_LIVRETACURVE",
                " RRRRRRRRRR",					// 9 parametres d'entree + 1 parametre de retour
                " PXL_ARM_INF_LivretACurve",
                " asOfDate,infCurvId, euribCurvId, [flagRounding], [resetmanager],[fixingLivretAId],[fixingEuribId],[monthForAugust],[monthForFebruary]",
                " 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a livretA curve",
                " As of date",
				" Inflation curve id ",
                " Euribor curve id",
                " rounded flag",
				" reset Manager for inflation curve",
				" Fixing taux Livret A id",
				" Fixing taux euribor id",
				" reset month for August LivretA Reset (default JUN)",
				" reset month for February LivretA Reset (default NOV)"
        },
		{
				" Local_ARM_LivreACurveGetRateDate",
                " RRR",					// 2 parametres d'entree + 1 parametre de retour
                " ARM_LivretACurveGetRateDate",
                " livretACurvId, date",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Get the begin date of a period in which we take the rate for the date given",
                " Livret A id",
				" Date"
		},
		{
        		" Local_InfIdxCreate",		/// name of the C++ function
                " RRRRRR",					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_INF_CreateIdx",
				" indexName,[resetLag],[DCFLag],[ccyObj]", 
                " 1",						/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_CreateIdx,
                " Creates an inflation index",
				" index Name like CPALEMU",	
				" resetLag like 0m, Default per index",
				" DCFLag like 3m, Default per index",
				" currency Object, Default per index"
        },
        {
        		" Local_PXL_InfIdxCreate",	/// name of the C++ function
                " RRRRRR",					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_CreateIdx",
				" indexName,[resetLag],[DCFLag],[ccyObj]", 
                " 0",						/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates an inflation index",
				" index Name like CPALEMU",	
				" resetLag like 0m, Default per index",
				" DCFLag like 3m, Default per index",
				" currency Object,Default per index"
		},
        {
        		" Local_InfLegwDateStripCreate",/// name of the C++ function
                " RRRRRRRRRR",					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_INF_CreateGenericLeg",
				" infidx,rcvOrPay,[interpType],[multiple],[constant],[notionalType],[NumDateStrip],[DenomDateStrip]",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_CreateGenericLeg,
                " Creates an inflation leg generic",
				" index Name like CPALEMU",	
				" REC or PAY, default PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default 1",
				" numeric, default 0",
				" NXNONE, NXEND NXINFEND, default NXNONE",
				" date strip object for the numerator, default No Object",
				" date strip object for the denominator, default No Object"
        },
        {
        		" Local_PXL_InfLegwDateStripCreate",/// name of the C++ function
                " RRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_CreateGenericLeg",
				" infidx,rcvOrPay,[interpType],[multiple],[constant],[NotionalType],[NumDateStrip],[DenomDateStrip]",
                " 0",								/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates an inflation leg generic",
				" index Name like CPALEMU",	
				" REC or PAY, default PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default 1",
				" numeric, default 0",
				" NXNONE, NXEND NXINFEND, default NXNONE",
				" date strip object for the numerator, default No Object",
				" date strip object for the denominator, default No Object"
        },
        {
        		" Local_InfLegYtYCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",		/// 20 parametres =  19 d'entree + 1 parametre de retour 
                " ARM_INF_CreateYtYLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[leverage],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetNumGap],[resetDenomGap],[payFreq],[payGap],[payCal],[adjFirstDate],[CoLeverage]",
                " 1",						
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_CreateYtYLeg,
                " Creates an inflation leg year to year",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",
				" numeric, default is 0",
				" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default is A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE default is 30/360",
				" string: ccy mainly, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default is UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default is SS",
				" in days (or roll date), default is 0",
				" in days (or roll date), default is numResetGap",
				" A, S, Q, B, M, W, D, Z, default is resetFreq", 
				" in days, default is 0",
				" string: ccy mainly, default is ccy of index",
				" true or false, default false"
				" numeric, default is 1",
       },
       {
        		" Local_PXL_InfLegYtYCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRRR",		/// 20 parametres =  19 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_CreateYtYLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[leverage],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetNumGap],[resetDenomGap],[payFreq],[payGap],[payCal],[adjFirstDate],[Coleverage]",
                " 0",							/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates an inflation leg year to year",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",

				" numeric, default is 0",
				" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default is A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE default is 30/360",
				" string: ccy mainly, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default is UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default is SS",
				" in days (or roll date), default is 0",
				" in days (or roll date), default is numResetGap",
				" A, S, Q, B, M, W, D, Z, default is resetFreq", 
				" in days, default is 0",
				" string: ccy mainly, default is ccy of index",
				" true or false, default false"
				" numeric, default is 1",
        },
        {
        		" Local_InfLegZCCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRR",				/// 15 parametres =  14 d'entree + 1 parametre de retour 
                " ARM_INF_CreateZCLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[leverage],[spread],[resetCal],[resetNumGap],[resetDenomGap],[payGap],[payCal],[adjFirstDate],[firstReset]",
                " 1",						
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_CreateZCLeg,
                " Creates an inflation leg year to year",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",
				" numeric, default is 0",
				" string: ccy mainly, default is INF",
				" in days (or roll date), default is 0",
				" in days (or roll date), default is numResetGap",
				" A, S, Q, B, M, W, D, Z, default is resetFreq", 
				" in days, default is 0",
				" string: ccy mainly, default is ccy of index",
				" true or false, default false",
				" value of the first Reset on 100 basis, default is taken from the curve"
       },
       {
        		" Local_PXL_InfLegZCCreate",	/// name of the C++ function
                " ",				/// 15 parametres =  14 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_CreateZCLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[leverage],[spread],[resetCal],[resetNumGap],[resetDenomGap],[payGap],[payCal],[adjFirstDate],[firstReset]",
                " 0",							/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates an inflation leg year to year",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",
				" numeric, default is 0",
				" string: ccy mainly, default is INF",
				" in days (or roll date), default is 0",
				" in days (or roll date), default is numResetGap",
				" A, S, Q, B, M, W, D, Z, default is resetFreq",
				" in days (or roll date), default is 0",
				" in days (or roll date), default is numResetGap",
				" string: ccy mainly, default is ccy of index",
				" true or false, default false",
				" value of the first Reset on 100 basis, default is taken from the curve"
        },
		{
        		" Local_InfLegOATCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",	/// 21 parametres =  20 d'entree + 1 parametre de retour 
                " ARM_INF_CreateOATLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[coupon],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetNumGap],[resetDenomGap],[payFreq],[payGap],[payCal],[notionalType],[firstReset]",
                " 1",						
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_CreateOATLeg,
                " Creates an OAT Type inflation leg",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" coupon of the oat in 100 basis, default 1 percent",
				" in 100 basis, default 0",
				" A, S, Q, B, M, W, D, Z, default A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE default 30/360",
				" string, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default  SS",
				" in days, default 0",
				" A, S, Q, B, M, W, D, Z, default resetFreq", 
				" in days (or roll date), default 0",
				" in days (or roll date), default is numResetGap",
				" string: ccy mainly, default is ccy of index",
				" NXNONE, NXEND NXINFEND, NXASINFEND, default NXINFEND",
				" value of first Reset, default taken from the curve"
       },
       {
        		" Local_PXL_InfLegOATCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",		/// 21 parametres =  20 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_CreateOATLeg",
				" start,end,infIdxId,rcvOrPay,[interpTyp],[coupon],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetNumGap],[resetDenomGap],[payFreq],[payGap],[payCal],[notionalType],[firstReset]",
                " 0",							/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates an OAT Type inflation leg",
                " date",
                " date",
				" index Name like CPALEMU",	
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" coupon of the oat in 100 basis, default 1 percent",
				" in 100 basis, default 0",
				" A, S, Q, B, M, W, D, Z, default A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE default 30/360",
				" string, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default  SS",
				" in days, default 0",
				" A, S, Q, B, M, W, D, Z, default resetFreq", 
				" in days (or roll date), default 0",
				" in days (or roll date), default is numResetGap",
				" string: ccy mainly, default is ccy of index",
				" NXNONE, NXEND NXINFEND, NXASINFEND, default NXINFEND",
				" value of first Reset, default taken from the curve"
        },
        {
        		" Local_FixZC",			/// name of the C++ function
                " RRRRRRRRRRR",			/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_FixZCCreate",
				" startDate, endDate, fixRate, rcvOrPay, [dayCount],[intRule],[stubRule],[payGap],[payCalendar],[discountCcy]",
                " 1",					/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_FixZCCreate,
                " Creates a fix zero coupon leg",
                " date",
                " date",
				" rate in 100 basis",
				" REC or PAY",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE default is 30/360",
				" ADJ, UNADJ, default UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default  SS",
				" in days (or roll date), default is 0",
				" string: ccy mainly, default is ccy of index",
				" currency object"
        },
        {
        		" Local_PXL_FixZC",		/// name of the C++ function
                " RRRRRRRRRRR",			/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " PXL_ARM_FixZCCreate",
				" startDate, endDate, fixRate, rcvOrPay, [dayCount],[intRule],[stubRule],[payGap],[payCalendar],[discountCcy]",
                " 0",					/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a fix zero coupon leg",
                " date",
                " date",
				" rate in 100 basis",
				" REC or PAY",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is 30/360",
				" ADJ, UNADJ, default UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default  SS",
				" in days (or roll date), default is 0",
				" string: ccy mainly, default is ccy of index",
				" currency object"
        },
        {
        		" Local_ResetManager",			/// name of the C++ function
                " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " ARM_INF_ResetManagerCreate",
				" resetData",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                IDH_ARM_INF_ResetManagerCreate,
                " Creates a reset manager object",
				" should be an array"
        },
        {
        		" Local_PXL_ResetManager",		/// name of the C++ function
                " RR",							/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " PXL_ARM_INF_ResetManagerCreate",
				" resetData",
                " 0",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a reset manager object",
				" should be an array"
        },

		{
        		" Local_CorrelMat_Create",		/// name of the C++ function
                " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CorrelMat_Create",
				" asOfDate,X,Y,Z",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a simple correlation object",
				" date",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics"
        },
        {
        		" Local_PXL_CorrelMat_Create",	/// name of the C++ function
                " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " PXL_ARM_CorrelMat_Create",
				" asOfDate,X,Y,Z",
                " 0",							/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a simple correlation object",
				" date",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics"
        },
		{
        		" Local_ComputeCorrel",		/// name of the C++ function
                " RRRR",				/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_ComputeCorrel",
				" correlMat,X,Y",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Get the correl from a correl matrix",
				" Correl Object",
				" double",
				" double"
        },
		{
        		" Local_CorrelManager_Create",		/// name of the C++ function
                " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CorrelManager_Create",
				" mktTag,intraMktTag,asOfDate,X,Y,Z",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" date",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics"
        },
        {
        		" Local_PXL_CorrelManager_Create",	/// name of the C++ function
                " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " PXL_ARM_CorrelManager_Create",
				" mktTag,intraMktTag,asOfDate,X,Y,Z",
                " 0",							/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" date",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics"
        },
		{
        		" Local_CorrelManagerFromMat_Create",/// name of the C++ function
                " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CorrelManagerFromMat_Create",
				" mktTag,intraMktTag,correlMatrix",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager from a correlation matrix",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" correlation matrix object"
        },
        {
        		" Local_PXL_CorrelManagerFromMat_Create",	/// name of the C++ function
                " RRRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " PXL_ARM_CorrelManagerFromMat_Create",
				" mktTag,intraMktTag,correlMatrix",
                " 0",										/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager from a correlation matrix",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" correlation matrix object"
        },
		{
        		" Local_CreateGenCorrelManager",/// name of the C++ function
                " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CreateGenCorrelManager",
				" mktTags,intraMktTags,correlCurves",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager",
				" tags of the type mkt1/mkt2 sorted in alphabetical order! (array)",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER (array)",
				" correl curves (array)"
		},
		{
        		" Local_PXL_CreateGenCorrelManager",/// name of the C++ function
                " RRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " PXL_ARM_CreateGenCorrelManager",
				" mktTags,intraMktTags,correlCurves",
                " 0",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Creates a correlation manager",
				" tags of the type mkt1/mkt2 sorted in alphabetical order! (array)",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER (array)",
				" correl curves (array)"
		},
		{
        		" Local_CorrelManager_Fill",		/// name of the C++ function
                " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CorrelManager_Fill",
				" mktTag,intraMktTag,X,Y,Z,correlManagerId",
                " 1",								/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Fills a correlation manager",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics",
				" correl manager object"
        },
        {
        		" Local_PXL_CorrelManager_Fill",	/// name of the C++ function
                " RRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " PXL_ARM_CorrelManager_Fill",
				" mktTag,intraMktTag,X,Y,Z,correlManagerId",
                " 0",								/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Fills a correlation manager",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" array of numerics or maturities",
				" array of numerics or maturities",
				" matrix of numerics",
				" correl manager object"
        },
		{
        		" Local_CorrelManagerFromMat_Fill",			/// name of the C++ function
                " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CorrelManagerFromMat_Fill",
				" mktTag,intraMktTag,correlMatId,correlManagerId",
                " 1",										/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Fills a correlation manager with a correlation matrix",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" correl manager object"
        },
        {
        		" Local_PXL_CorrelManagerFromMat_Fill",		/// name of the C++ function
                " RRRRR",									/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " PXL_ARM_CorrelManagerFromMat_Fill",
				" mktTag,intraMktTag,correlMatId,correlManagerId",
                " 0",										/// not visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " Fills a correlation manager with a correlation matrix",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" correl manager object"
        },
		{
        		" Local_ComputeCorrelFromManager",		/// name of the C++ function
                " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_ComputeCorrelFromManager",
				" mktTag,intraMktTag,X,Y,correlManagerId",
                " 1",									/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
                " ",
                " computes correlation from the correlation manager",
				" tag of the type mkt1/mkt2 sorted in alphabetical order!",
				" can be ccy1_Ten1_Ccy2_DIAG, ccy1_Ten1_Ccy2_Ten2_COMP | NOCHECK, or ccy|ind1_ccy|ind2 in ALPHABETICAL ORDER",
				" double",
				" double",
				" correl manager object"
        },
		{
        		" Local_Vol_to_Cor",	/// name of the C++ function
                " RRRR",				/// 4 parametres = " d'entree + 1 parametre de retour 
                " ARM_INF_VolToCor",	  /// nom s/ Excel
				" ZCVol, YtYVol, Maturity",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes the correlation from volatility data",
				" Zero Coupon volatilities vector",
				" Year to Year volatilities vector",
				" maturities vector"
        },
		{
        		" Local_YtYCor_to_ZC",	/// name of the C++ function
                " RRRR",				/// 4 parametres = " d'entree + 1 parametre de retour 
                " ARM_INF_YtYCor_to_ZC",	  /// nom s/ Excel
				" YtYVol, Cor, Maturity",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes zero zoupon volatilities from these following data",
				" the correlation vector",
				" year to year volatilities vector",
				" maturities vector"
        },
		{
        		" Local_ZCCor_to_YtY",	/// name of the C++ function
                " RRRRRRR",				/// 8 parametres = " d'entree + 3 parametre de retour 
				" ARM_INF_ZCCor_to_YtY",	  /// nom s/ Excel
				" ZCVol, Cor, Maturity",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes year to year volatilities from these following data",
				" the correlation vector",
				" zero zoupon volatilities vector",
				" maturities vector"
        },
		{
        		" Local_Bounds",	/// name of the C++ function
                " RRRRRR",				/// 5 parametres = " d'entree + 1 parametre de retour 
                " ARM_INF_Bounds",	  /// nom s/ Excel
				" ZCVol, YtYVol, Cor, Maturity, Choice",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " provides confidence interval for zero coupon and year to year volatilities and correlations",
				" zero coupon volatlities", 
                " year to year volatlities",
				" correlations",
				" Maturities",
				" Choose between ZC, YtY and Cor"
        },
		{
        		" Local_HmgVol_to_Cor",	/// name of the C++ function
                " RRRRR",				/// 4 parametres = " 3 d'entree + 1 parametre de retour 
                " ARM_INF_HmgVolToCor",	  /// nom s/ Excel
				" ZCVol, YtYVol, Maturity, length",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes the correlation from volatility data",
				" Zero Coupon volatilities vector",
				" Year to Year volatilities vector",
				" maturities vector",
				" length for interpoling zero coupon volatilities"
        },
		{
        		" Local_HmgZCCor_to_YtY",	/// name of the C++ function
                " RRRRR",					/// 5 parametres = " 4 d'entree + 1 parametre de retour 
				" ARM_INF_HmgZCCor_to_YtY",	  /// nom s/ Excel
				" ZCVol, Cor, Maturity, length",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes year to year volatilities from these following data",
				" the correlation vector",
				" zero zoupon volatilities vector",
				" maturities vector",
				" length for interpoling zero coupon volatilities"
        },
		{
        		" Local_HmgYtYCor_to_ZC",	/// name of the C++ function
                " RRRRR",				/// 5 parametres = " 4 d'entree + 1 parametre de retour 
                " ARM_INF_HmgYtYCor_to_ZC",	  /// nom s/ Excel
				" YtYVol, Cor, Maturity, length",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes zero zoupon volatilities from these following data",
				" the correlation vector",
				" year to year volatilities vector",
				" maturities vector",
				" length for interpoling zero coupon volatilities"
        },
		{
        		" Local_VolYoY_to_VolSwp",	/// name of the C++ function
                " RRRRRRRRR",				/// 9 parametres = " 8 d'entree + 1 parametre de retour 
                " ARM_INF_VolYoY_to_VolSwp",	  /// nom s/ Excel
				" DFactor, FwdCPI, Vol_DF, Vol_YoY, AvgCor, Dates, Tenors, SwpRate",
                " 1",							/// visible in excel
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Computes the implied vol of the inflation swap from the YoY vols and the DF's ones",
				" Discount Factor vector",
				" Fwd CPI vector",
				" The Discount Factor's volatilities vector",
				" The year to year volatilities vector",
				" The averages correlations betweens the DF and the fwd  CPI vector",
				" Maturities vector",
				" Tenors vector",
				" The inf swap rate"
        },
		{
        		" Local_SeasonalityManager_Create",		/// name of the C++ function
                " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour
                " ARM_INF_SeasonalityManager_Create",
				" monthList,seasonSpreadList,[horizonList], [mode]",
                " 1",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
				" ",
                " Creates a seasonality manager object",
				" should be an array of months (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)",
				" should be an array of mumeric (1, 2, ..., Last InfCurve Expiry)",
				" should be an arrayof numeric",
				" correction mode : PLUS or MULT"
        },
        {
        		" Local_PXL_SeasonalityManager_Create",		/// name of the C++ function
                " RRRRR",								/// 5 parametres = 4 d'entree + 1 parametre de retour
                " PXL_ARM_INF_SeasonalityManager_Create",
				" monthList,seasonSpreadList,[horizonList], [mode]",
                " 0",							/// visible in excel
                XLLOCALARM_SEC_GROUP,
                " ",
				" ",
                " Creates a seasonality manager object",
				" should be an array of months (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)",
				" should be an array of mumeric (1, 2, ..., Last InfCurve Expiry)",
				" should be an arrayof numeric",
				" correction mode : PLUS or MULT"
        },
		{
        		" Local_InfCapFloorCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",	/// 21 parametres =  20 d'entree + 1 parametre de retour
                " ARM_INF_CapFloorCreate",
				" start,end,infIdxId,capOrFloor,strike,swaptype,rcvOrPay,[interpTyp],[leverage],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCal]",
                " 1",						
                XLLOCALARM_OPTION_GROUP,
                " ",
                IDH_ARM_INF_CapFloorCreate,
                " Creates an inflation cap floor option",
                " date",
                " date",
				" index Name like CPALEMU",	
				" CAP or FLOOR",
				" double",
				" ZC YTY OAT",
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",
				" numeric, default is 0",
				" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default is A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is 30/360",
				" string: ccy mainly, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default is UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default is SS",
				" in days, default is 0",
				" A, S, Q, B, M, W, D, Z, default is resetFreq", 
				" in days, default is 0",
				" string: ccy mainly, default is ccy of index"
       },
       {
        		" Local_PXL_InfCapFloorCreate",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",	/// 21 parametres =  20 d'entree + 1 parametre de retour
                " PXL_ARM_INF_CapFloorCreate",
				" start,end,infIdxId,capOrFloor,strike,swaptype,rcvOrPay,[interpTyp],[leverage],[spread],[resetFreq],[dayCount],[resetCal],[fwdRule],[intRule],[stubRule],[resetGap],[payFreq],[payGap],[payCal]",
                " 0",							/// not visible in excel
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation cap floor option",
                " date",
                " date",
				" index Name like CPALEMU",	
				" CAP or FLOOR",
				" double",
				" ZC YTY OAT",
				" REC or PAY",
				" CPILINEAR, STEPWISESTART, STEPWISEMID, STEPWISEEND, Default CPILINEAR",
				" numeric, default is 1",
				" numeric, default is 0",
				" A(nnual), S(emiAnnual), Q(uarterly), B(imonthly), M(onthly), W(eekly), D(aily), Z(eroCoupon), default is A",
				" ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, default is 30/360",
				" string: ccy mainly, default is INF",
				" F(ollowing), MF (modified Following), P(revious), MP(modified Previous), default is MF",
				" ADJ, UNADJ, default is UNADJ",
				" SS (shortStart), LS (LongStart), SE(ShortEnd), LE(LongEnd), default is SS",
				" in days, default is 0",
				" A, S, Q, B, M, W, D, Z, default is resetFreq", 
				" in days, default is 0",
				" string: ccy mainly, default is ccy of index"
        },
        {
				" Local_InfYCMOD",
                " RRR",						// 3 parametres = 2 d'entree + 1 parametre de retour
                " ARM_INF_YCMOD",
                " discountCurve,CPIForwardCurve",
                " 1",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_YCMOD,
                " Model for the inflation curve: takes a discounting curve and a forward CPI Curve",
                " ARM_ZEROCURVE object,used to compute the discount factors",
                " ARM_INFCURVE object,used to compute the CPI forward"
        },
		{
				" Local_PXL_InfYCMOD",
                " RRR",						// 3 parametres = 2 d'entree + 1 parametre de retour
                " PXL_ARM_INF_YCMOD",
                " discountCurve,CPIForwardCurve",
                " 0",
                XLLOCALARM_MODEL_GROUP,
                " ",
                " ",
                " Model for the inflation curve: takes a discounting curve and a forward CPI Curve",
                " ARM_ZEROCURVE object,used to compute the discount factors",
                " ARM_INFCURVE object,used to compute the CPI forward"
        },
   		{
				" Local_InfBSMOD",
                " RRRRRRRRR",						/// 8 parametres = 7 d'entree + 1 parametre de retour
                " ARM_INF_BSMOD",
                " asOfDate,discountCurve,CPIForwardCurve,volData,[correlManager],[IRBSmodel],[InflationSwaptionCurve],[IRSwaptionCurve]",
                " 1",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Model for the inflation curve: takes a discounting curve and a forward CPI Curve",
				" Date",
                " ARM_ZEROCURVE object,used to compute the discount factors",
                " ARM_INFCURVE object,used to compute the CPI forward",
				" ARM_VOLCURVE object,either volcurve or volcube object",
				" ARM_CORRELMANAGER object",
				" ARM_IRBSModel type object",
				" Inflation swaption vol curve",
				" Interest Rate swaption vol curve"
        },
		{
				" Local_PXL_InfBSMOD",
                " RRRRRRRRR",						/// 8 parametres = 7 d'entree + 1 parametre de retour
                " PXL_ARM_InfBSMOD",
                " asOfDate,discountCurve,CPIForwardCurve,volData,[correlManager],[IRBSmodel],[InflationSwaptionCurve],[IRSwaptionCurve]",
                " 0",
                XLLOCALARM_MODEL_GROUP,
                " ",
                " ",
                " Model for the inflation curve: takes a discounting curve and a forward CPI Curve",
				" Date",
                " ARM_ZEROCURVE object,used to compute the discount factors",
                " ARM_INFCURVE object,used to compute the CPI forward",
				" ARM_VOLCURVE object,either volcurve or volcube object",
				" ARM_CORRELMANAGER object",
				" ARM_IRBSModel type object",
				" Inflation swaption vol curve",
				" Interest Rate swaption vol curve"
        },
		{
				" Local_SparseVolCube_CreateNFill",
                " RRRRRRRRRRRR",						// 12= 11 parametres d'entree + 1 parametre de retour
                " ARM_SparseVolCube_CreateNFill",
                " asOfDate,lastKnownDate,indexName,dim1Type,dim1Value,dim2Type,dim2Value,strikes,vols,[strikeType],[volType]",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a sparse volatility cube",
                " date",
                " date (last known published index)",
				" string",
				" either Tenor or TimeToStart or TimeToExpiry",
				" array of doubles or maturity",
				" either Tenor or TimeToStart or TimeToExpiry",
                " corresponding double value",
				" array of doubles",
				" matrix of data",
				" Y(ield) or P(rice) (default:Price)",
				" S(mile) or A(TM) (default:ATM)"
	    },
		{
				" Local_PXL_SparseVolCube_CreateNFill",
                " RRRRRRRRRRRR",						// 12= 11 parametres d'entree + 1 parametre de retour
                " PXL_ARM_SparseVolCube_CreateNFill",
                " asOfDate,lastKnownDate,indexName,dim1Type,dim1Value,dim2Type,dim2Value,strikes,vols,[strikeType],[volType]",
				" 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
				" ",
                " Creates a sparse volatility cube",
                " date",
                " date (last known published index)",
				" string",
				" either Tenor or TimeToStart or TimeToExpiry",
				" array of doubles or maturity",
				" either Tenor or TimeToStart or TimeToExpiry",
                " corresponding double value",
				" array of doubles",
				" matrix of data",
				" Y(ield) or P(rice) (default:Price)",
				" S(mile) or A(TM) (default:ATM)"
	    },
		{
				" Local_SparseVolCube_Fill",
                " RRRRRRRR",							// 8= 7 parametres d'entree + 1 parametre de retour
                " ARM_SparseVolCube_Fill",
                " sparseVolCubeId,dim1Type,dim1Value,otherDimType,otherDimValue,strikes,vols",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_SparseVolCube_Fill,
                " Creates a sparse volatility cube", 
				" Sparse Vol Cube Object",
				" either Tenor or TimeToStart or TimeToExpiry",
				" array of doubles or maturity",
				" either Tenor or TimeToStart or TimeToExpiry",
                " corresponding double value",
				" array of doubles",
				" matrix of data"
	    },
		{
				" Local_PXL_SparseVolCube_Fill",
                " RRRRRRRR",							// 8= 7 parametres d'entree + 1 parametre de retour
                " PXL_ARM_SparseVolCube_Fill",
                " sparseVolCubeId,dim1Type,dim1Value,otherDimType,otherDimValue,strikes,vols",
				" 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a sparse volatility cube", 
				" Sparse Vol Cube Object",
				" either Tenor or TimeToStart or TimeToExpiry",
				" array of doubles or maturity",
				" either Tenor or TimeToStart or TimeToExpiry",
                " corresponding double value",
				" array of doubles",
				" matrix of data"
	    },
		{
				" Local_VolCubeFromSparseVolCube",
                " RR",							// 2= 1 parametres d'entree + 1 parametre de retour
                " ARM_VolCubeFromSparseVolCube",
                " sparseVolCubeId",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Creates a volatility cube from a sparse volatility cube", 
				" Sparse Vol Cube Object"
	    },
		{
				" Local_PXL_VolCubeFromSparseVolCube",
                " RR",							// 2= 1 parametres d'entree + 1 parametre de retour
                " PXL_ARM_VolCubeFromSparseVolCube",
                " sparseVolCubeId",
				" 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a volatility cube from a sparse volatility cube", 
				" Sparse Vol Cube Object"
	    },
		{
				" Local_InfSwoVolCurveFromModel_Create",
                " RRRRRR",							// 6= 5 parametres d'entree + 1 parametre de retour
                " ARM_INF_InfSwoVolCurv_Create",
                " asOfDate,inflation/IR BS ModelId,[tenors],[expiries],[method]",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Computes the inflation swopt volatility surface from a model", 
				" date",
				" object",
				" vector",
				" vector",
				" Std, Equal Weight, DF Weight, DF Weight Square"
	    },
		{
				" Local_PXL_InfSwoVolCurveFromModel_Create",
                " RRRRRR",							// 6= 5 parametres d'entree + 1 parametre de retour
                " PXL_ARM_INF_InfSwoVolCurv_Create",
                " asOfDate,inflation/IR BS ModelId,[tenors],[expiries],[method]",
				" 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Computes the inflation swopt volatility surface from a model", 
				" date",
				" object",
				" vector",
				" vector",
				" Std, Equal Weight, DF Weight, DF Weight Square"
	    },
		{
				" Local_InfSwoVolCubeFromModel_Create",
                " RRRRRRRR",							// 8= 7 parametres d'entree + 1 parametre de retour
                " ARM_INF_InfSwoVolCube_Create",
                " asOfDate,inflation/IR BS ModelId,[tenors],[expiries],[smiledTenors],[strikes],[method]",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Computes the inflation swopt volatility Cube from a model", 
				" date",
				" object",
				" vector",
				" vector",
				" vector",
				" vector",
				" Std, Equal Weight, DF Weight, DF Weight Square"
	    },
		{
				" Local_PXL_InfSwoVolCubeFromModel_Create",
                " RRRRRRRR",							// 8= 7 parametres d'entree + 1 parametre de retour
                " PXL_ARM_INF_InfSwoVolCube_Create",
                " asOfDate,inflation/IR BS ModelId,[tenors],[expiries],[smiledTenors],[strikes],[method]",
				" 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Computes the inflation swopt volatility Cube from a model", 
				" date",
				" object",
				" vector",
				" vector",
				" vector",
				" vector",
				" Std, Equal Weight, DF Weight, DF Weight Square"
	    },
		{
				" Local_InfOATSwoVolCurveFromModel_Create",
                " RRRRRRRR",							// 7= 6 parametres d'entree + 1 parametre de retour
                " ARM_INF_InfOATSwoVolCurv_Create",
                " asOfDate,inflation/IR BS ModelId,[tenors],[expiries],[method],[coupon],[choice]",
				" 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                IDH_ARM_VolCubeFromSparseVolCube,
                " Computes the inflation OAT swopt volatility surface from a model", 
				" date",
				" object",
				" vector",
				" vector",
				" Std, Equal Weight, DF Weight, DF Weight Square",
				" double",
				" 0 or 1"
	    },
		{
				" Local_InfMultiBSMOD",
                " RR",						/// 2 parametres = 1 d'entree + 1 parametre de retour
                " ARM_INF_MultiBSMOD",
                " InfBSModels",
                " 1",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates a multi inflation bs model",
				" Vector of inflation bs models"
        },
   		{
				" Local_PXL_InfMultiBSMOD",
                " RR",						/// 2 parametres = 1 d'entree + 1 parametre de retour
                " PXL_ARM_INF_MultiBSMOD",
                " InfBSModels",
                " 0",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates a multi inflation bs model",
				" Vector of inflation bs models"
        },
   		{
				" Local_InfCurv_SetResetManager",
                " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour
                " ARM_INF_InfCurv_SetResetManager",
                " InfCurvId,ResetManagerId",
                " 1",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates an inflation curve with a reset manager",
				" inflation curve object",
				" reset manager object"
        },
   		{
				" Local_PXL_InfCurv_SetResetManager",
                " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour
                " PXL_ARM_INF_InfCurv_SetResetManager",
                " InfCurvId,ResetManagerId",
                " 0",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates an inflation curve with a reset manager",
				" inflation curve object",
				" reset manager object"
        },
   		{
				" Local_InfCurv_SetSeasonalityManager",
                " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour
                " ARM_INF_InfCurv_SetSeasonalityManager",
                " InfCurvId,SeasonalityManagerId",
                " 1",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates an inflation curve with a seasonality manager",
				" inflation curve object",
				" seasonality manager object"
        },
   		{
				" Local_PXL_InfCurv_SetSeasonalityManager",
                " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour
                " PXL_ARM_INF_InfCurv_SetSeasonalityManager",
                " InfCurvId,SeasonalityManagerId",
                " 0",
                XLLOCALARM_MODEL_GROUP,
                " ",
                IDH_ARM_INF_BSMOD,
                " Creates an inflation curve with a seasonality manager",
				" inflation curve object",
				" seasonality manager object"
        },
   		{
				" Local_ARM_GetSeasonMgrFromSummit",
                " RRRRRR",						/// 6 parametres = 5 d'entree + 1 parametre de retour
                " ARM_INF_GetSeasonMgrFromSummit",
                " Index,Ccy,CvName,AsOf,[Mode]",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a seasonality manager from Summit",
				" Index",
				" Currency",
				" Curve Name",
				" AsOfDate",
				" Seasonality Correction Mode (default : PLUS)"
        },
   		{
				" Local_PXL_ARM_GetSeasonMgrFromSummit",
                " RRRRRR",						/// 6 parametres = 5 d'entree + 1 parametre de retour
                " PXL_ARM_INF_GetSeasonMgrFromSummit",
                " Index,Ccy,CvName,AsOf,[Mode]",
                " 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a seasonality manager from Summit",
				" Index",
				" Currency",
				" Curve Name",
				" AsOfDate",
				" Seasonality Correction Mode (default : PLUS)"
        },
   		{
				" Local_ARM_GetResetMgrFromSummit",
                " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour
                " ARM_INF_GetResetMgrFromSummit",
                " AsOf, Index, [Source], [Ccy], [IsInflatIndex], [Term]",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a reset manager from Summit",
				" AsOfDate",
				" Index(EMU, IFRF, EURIB...) or 'FX' for spot FX",
				" Summit Source (Default MO)",
				" Currency (Default EUR) or couple of currencies if FX (like 'EUR_USD')",
				" Is it for an Inflation Index? (default Y, empty for FX)",
				" Term of Index (default 1D for inflation Index, mandatory for IR Index, empty for FX)"
        },
		{
				" Local_PXL_ARM_GetResetMgrFromSummit",
                " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour
                " PXL_ARM_INF_GetResetMgrFromSummit",
                " AsOf, Index, [Source], [Ccy], [IsInflatIndex], [Term]",
                " 0",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Creates a reset manager from Summit",
				" AsOfDate",
				" Index(EMU, IFRF, EURIB...) or 'FX' for spot FX",
				" Summit Source (Default MO)",
				" Currency (Default EUR) or couple of currencies if FX (like 'EUR_USD')",
				" Is it for an Inflation Index? (default Y, empty for FX)",
				" Term of Index (default 1D for inflation Index, mandatory for IR Index, empty for FX)"
        },
		{
        		" Local_ARM_GP_INFCAPFLOOR",	/// name of the C++ function
                " RRRR",	/// 4 parametres =  3d'entree + 1 parametre de retour
                " ARM_GP_INF_CapFloor",
				" swapId,capOrFloor,strike",
                " 1",						
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation cap floor option From a swap the two legs must have the same Freq",
                " swapId",
                " Cap or Floor",
				" strike double or RefValue"
       },
       {
        		" Local_PXL_ARM_GP_INFCAPFLOOR",	/// name of the C++ function
                " RRRR",	/// 4 parametres =  3d'entree + 1 parametre de retour
                " PXL_ARM_GP_INF_CapFloor",
				" swapId,capOrFloor,strike",
                " 0",							/// not visible in excel
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation digitak option from a swap the two legs must have the same Freq",
                " swapId",
                " Cap or Floor",
				" strike double or RefValue"
		},
		{
        		" Local_ARM_GP_INFCALLSPREAD",	/// name of the C++ function
                " RRRR",	/// 4 parametres =  3d'entree + 1 parametre de retour
                " ARM_GP_INF_CallSpread",
				" swapId,capOrFloor,strike",
                " 1",						
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation digital option from a swap the two legs must have the same Freq",
                " swapId",
                " Cap or Floor",
				" strike double or RefValue"
       },
       {
        		" Local_PXL_ARM_GP_INFCALLSPREAD",	/// name of the C++ function
                " RRRR",	/// 4 parametres =  3d'entree + 1 parametre de retour
                " PXL_ARM_GP_INF_CallSpread",
				" swapId,capOrFloor,strike",
                " 0",							/// not visible in excel
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation cap floor option From a swap the two legs must have the same Freq",
                " swapId",
                " Cap or Floor",
				" strike double or RefValue"
		},
		{
				" Local_ARM_GetReset",
                " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour
                " ARM_GetReset",
                " ResetMgr, date",
                " 1",
                XLLOCALARM_CURVE_GROUP,
                " ",
                " ",
                " Gets a reset from a reset manager",
				" Reset Manager Id",
				" date"
        },
		{
        		" Local_ARM_GP_INFDIGITAL",	/// name of the C++ function
                " RRRRRRR",	/// 7 parametres =  6d'entree + 1 parametre de retour
                " ARM_GP_INF_DIGITAL",
				" PaySwapId,DigitSwapId,PayOffType,Barrier,CapOrFloor,PayOrRec",
                " 1",						
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation Digital for Two Payoff Types: (PayLeg)*I_{DigitLeg>Barrier} (PayLeg*DigitLeg)*I_{DigitLeg>Barrier}",
                " PaySwapLegId : FIX or Libor/CMS or INF leg",
                " DigitSwapLegId : Libor/CMS or INF leg",
				" PayOffType : SIMPLE or PRODUCT",
				" Barrier double or RefValue",
				" CapOrFloor 'C' for Down and 'F' for Up",
				" PayOrRec for the PaySwapleg"
       },
       {
        		" Local_PXL_ARM_GP_INFDIGITAL",	/// name of the C++ function
                " RRRRRRR",	/// 7 parametres =  6d'entree + 1 parametre de retour
                " PXL_ARM_GP_INF_DIGITAL",
				" PaySwapId,DigitSwapId,PayOffType,Barrier,CapOrFloor,PayOrRec",
                " 0",							/// not visible in excel
                XLLOCALARM_OPTION_GROUP,
                " ",
                " ",
                " Creates an inflation Digital for Two Payoff Types: (PayLeg)*I_{DigitLeg>Barrier} (PayLeg*DigitLeg)*I_{DigitLeg>Barrier}",
                " PaySwapLegId : FIX or Libor/CMS or INF leg",
                " DigitSwapLegId : Libor/CMS or INF leg",
				" PayOffType : SIMPLE or PRODUCT",
				" Barrier double or RefValue",
				" CapOrFloor 'C' for Down and 'F' for Up",
				" PayOrRec for the PaySwapleg"
		},
		{
        	" Local_ARM_INF_HybridInfIrMkt_Create",		/// name of the C++ function
            " RRRR",										/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_INF_HybridInfIrMkt_Create",
			" asOfDate, mktDatasKeys, mktDatasId",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Hybrid Mkt Datas object",
			" as of date",
			" list of keys",
			" list of objects id",
		},
		{
        	" Local_PXL_ARM_INF_HybridInfIrMkt_Create",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_INF_HybridInfIrMkt_Create",
			" asOfDate, mktDatasKeys, mktDatasId",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Hybrid Mkt Datas object",
			" as of date",
			" list of keys",
			" list of objects id",
		},
		{
        	" Local_ARM_INF_HybridInfIrPayOff_Create",				/// name of the C++ function
            " RRRRRRRRRRRRRRR",										/// 15 parametres = 14 d'entree + 1 parametre de retour 
            " ARM_INF_HybridInfIrPayOff_Create",
			" CstCpnCoef, CstOptCoef, MainCpnName, MainCpnCoef, MainOptName, MainOptCoef, [SubCpnName], [SubCpnCoef], [SubOptName], [SubOptCoef], [SupCpnName], [SupCpnCoef], [SupOptName], [SupOptCoef]",
            " 1",												/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Hybrid PayOff",
			" Coef of Cst  Index Coupon",
			" Coef of Cst  Index Option",
			" Name of Main Index Coupon",
			" Coef of Main Index Coupon",
			" Name of Main Index Option",
			" Coef of Main Index Option",
			" Name of Sub  Index Coupon",
			" Coef of Sub  Index Coupon",
			" Name of Sub  Index Option",
			" Coef of Sub  Index Option",
			" Name of Sup  Index Coupon",
			" Coef of Sup  Index Coupon",
			" Name of Sup  Index Option",
			" Coef of Sup  Index Option",
		},
		{
        	" Local_PXL_ARM_INF_HybridInfIrPayOff_Create",	/// name of the C++ function
            " RRRRRRRRRRRRRRR"								/// 15 parametres = 14 d'entree + 1 parametre de retour 
            " PXL_ARM_INF_HybridInfIrPayOff_Create",
			" CstCpnCoef, CstOptCoef, MainCpnName, MainCpnCoef, MainOptName, MainOptCoef, [SubCpnName], [SubCpnCoef], [SubOptName], [SubOptCoef], [SupCpnName], [SupCpnCoef], [SupOptName], [SupOptCoef]",
			" 0",											/// visible in vba
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Hybrid PayOff",
			" Coef of Cst  Index Coupon",
			" Coef of Cst  Index Option",
			" Name of Main Index Coupon",
			" Coef of Main Index Coupon",
			" Name of Main Index Option",
			" Coef of Main Index Option",
			" Name of Sub  Index Coupon",
			" Coef of Sub  Index Coupon",
			" Name of Sub  Index Option",
			" Coef of Sub  Index Option",
			" Name of Sup  Index Coupon",
			" Coef of Sup  Index Coupon",
			" Name of Sup  Index Option",
			" Coef of Sup  Index Option",
		},
		{
        	" Local_ARM_INF_HybridInfIr_Load",					/// name of the C++ function
            " RRRRR",											/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_INF_HybridInfIr_Load",
			" Instrument, [MktDatas], [Model], [Payoff]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " load all objects fro pricing",
			" Intrument ( Leg_Create )",
			" Market Data ( Mkt_Create )",
			" Model  ( Mod_Create )",
			" Payoff ( Payoff_Create )",
		},
		{
        	" Local_PXL_ARM_INF_HybridInfIr_Load",			/// name of the C++ function
            " RRRRR",										/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_INF_HybridInfIr_Load",
			" Instrument, [MktDatas], [Model], [Payoff]",
            " 0",											/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " load all objects fro pricing",
			" Intrument ( Leg_Create )",
			" Market Data ( Mkt_Create )",
			" Model  ( Mod_Create )",
			" Payoff ( Payoff_Create )",
		},
		{
			" Local_InfBSSmiledModel",
            " RRRRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour
            " ARM_INF_BSSMILEDMODEL",
            " asOfDate,discountCurve,CPIForwardCurve,VolSigma,VolNu,VolRho,VolBeta,VolAtmCap,Correl,[CorrelAdj]",
            " 1",
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Model for Smiled Inflation volatility",
			" ARM_DATE object which contains the asofdate",
            " ARM_ZEROCURVE object,used to compute the discount factors",
            " ARM_INFCURVE object,used to compute the CPI forward",
			" ARM_VOLCURVE the sigma Sabr parameter",
			" ARM_VOLCURVE the nu Sabr parameter",
			" ARM_VOLCURVE the rho Sabr parameter",
			" ARM_VOLCURVE the beta Sabr parameter",
			" ARM_VOLCURVE the bs atm ir vol",
			" ARM_VOLCURVE the correlation",
			" ARM_VOLCURVE the adj conv correlation",
        },
		{
			" Local_PXL_InfBSSmiledModel",
            " RRRRRRRRRRR",						/// 10 parametres = 9 d'entree + 1 parametre de retour
            " PXL_ARM_INF_BSSMILEDMODEL",
            " asOfDate,discountCurve,CPIForwardCurve,VolSigma,VolNu,VolRho,VolBeta,VolAtmCap,Correl,[CorrelAdj]",
            " 0",
            XLLOCALARM_MODEL_GROUP,
            " ",
            " ",
            " Model for Smiled Inflation volatility",
			" ARM_DATE object which contains the asofdate",
            " ARM_ZEROCURVE object,used to compute the discount factors",
            " ARM_INFCURVE object,used to compute the CPI forward",
			" ARM_VOLCURVE the sigma Sabr parameter",
			" ARM_VOLCURVE the nu Sabr parameter",
			" ARM_VOLCURVE the rho Sabr parameter",
			" ARM_VOLCURVE the beta Sabr parameter",
			" ARM_VOLCURVE the bs atm ir vol",
			" ARM_VOLCURVE the correlation",
			" ARM_VOLCURVE the adj conv correlation",
        },
		{
        	" Local_ARM_INF_GetPrice",	/// name of the C++ function
            " RRR",						/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_INF_GetPrice",
			" Inf Pricer,[Key]",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Pricing Function",
			" Inf Pricer",
			" Key ( Price ( by def)  or Flows)",
		},
		{
        	" Local_PXL_ARM_INF_GetPrice",	/// name of the C++ function
            " RRR",						/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_INF_GetPrice",
			" Inf Pricer,[Key]",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Pricing Function",
			" Inf Pricer",
			" Key ( Price ( by def)  or Flows)",
		},
		{
        	" Local_ARM_INF_GetSchedule",	/// name of the C++ function
            " RRR",						/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_INF_GetSchedule",
			" Inf Leg,Key",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Edit Schedule",
			" Inf Pricer",
			" Key ( StartDates, EndDates, ResDates, PayDates)",
		},
		{
        	" Local_PXL_ARM_INF_GetSchedule",	/// name of the C++ function
            " RRR",								/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_INF_GetSchedule",
			" Inf Leg,Key",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Pricing Function",
			" Inf Pricer",
			" Key (StartDates, EndDates, ResDates, PayDates)",
		},
//=======> PAYOFF
		{
        	" Local_ARM_INF_SpreadCap_Create",				
            " RRRRRRRRR",										
            " ARM_INF_SpreadCap_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, Strike, Notional ",
            " 1",	
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Cap PayOff",
			" Main Index Name",
			" Main Index Type (INF, IR)",
			" Main Leverage Value",
			" Sub  Index Name",
			" Sub  Index Type (INF, IR)",
			" Sub  Leverage Value",
			" Strike Value",
			" Notional Value",
		},
		{
        	" Local_PXL_ARM_INF_SpreadCap_Create",
            " RRRRRRRRR"								 
            " PXL_ARM_INF_SpreadCap_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, Strike, Notional ",
            " 0",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Cap PayOff",
			" Main Index Name",
			" Main Index Type (INF, def: IR)",
			" Main Leverage Value",
			" Sub Index Name",
			" Sub Index Type (INF, def: IR)",
			" Sub  Leverage Value",
			" Strike Value",
			" Notional Value",
		},
		{
        	" Local_ARM_INF_SpreadDigital_Create",				
            " RRRRRRRRR",										
            " ARM_INF_SpreadDigital_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, Strike, Notional ",
            " 1",	
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Digital PayOff",
			" Main Index Name",
			" Main Index Type (INF, IR)",
			" Main Leverage Value",
			" Sub  Index Name",
			" Sub  Index Type (INF, IR)",
			" Sub  Leverage Value",
			" Strike Value",
			" Notional Value",
		},
		{
        	" Local_PXL_ARM_INF_SpreadDigital_Create",
            " RRRRRRRRR"								 
            " PXL_ARM_INF_SpreadDigital_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, Strike, Notional ",
            " 0",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Digital PayOff",
			" Main Index Name",
			" Main Index Type (INF, def: IR)",
			" Main Leverage Value",
			" Sub Index Name",
			" Sub Index Type (INF, def: IR)",
			" Sub  Leverage Value",
			" Strike Value",
			" Notional Value",
		},
		{
        	" Local_ARM_INF_DoubleDigital_Create",				
            " RRRRRRRRRRR",										
            " ARM_INF_DoubleDigital_Create",
			" MainIndex, [MainType], MainLeverage, MainStrike, SubIndex, [SubType], SubLeverage, SubStrike, Spread,  Notional ",
            " 1",	
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Double Digital PayOff",
			" Main Index Name",
			" Main Index Type (INF, IR)",
			" Main Leverage Value",
			" Main Strike Value",
			" Sub  Index Name",
			" Sub  Index Type (INF, IR)",
			" Sub  Leverage Value",
			" Sub  Strike Value",
			" Spread Value",
			" Notional Value",
		},
		{
        	" Local_PXL_ARM_INF_DoubleDigital_Create",
            " RRRRRRRRRRR",									 
            " PXL_ARM_INF_DoubleDigital_Create",
			" MainIndex, [MainType], MainLeverage, MainStrike, SubIndex, [SubType], SubLeverage, SubStrike, Spread, Notional ",
            " 0",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Double Digital PayOff",
			" Main Index Name",
			" Main Index Type (INF, def: IR)",
			" Main Leverage Value",
			" Main Strike Value",
			" Sub Index Name",
			" Sub Index Type (INF, def: IR)",
			" Sub  Leverage Value",
			" Sub  Strike Value",
			" Spread Value",
			" Notional Value",
		},
		{
        	" Local_ARM_INF_Corridor_Create",				
            " RRRRRRRRRR",										
            " ARM_INF_Corridor_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, StrikeInf, StrikeSup, Notional ",
            " 1",	
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Cap PayOff",
			" Main Index Name",
			" Main Index Type (INF, IR)",
			" Main Leverage Value",
			" Sub  Index Name",
			" Sub  Index Type (INF, IR)",
			" Sub  Leverage Value",
			" Strike Inf Value",
			" Strike Sup Value",
			" Notional Value",
		},
		{
        	" Local_PXL_ARM_INF_Corridor_Create",
            " RRRRRRRRR"								 
            " PXL_ARM_INF_Corridor_Create",
			" MainIndex, [MainType], MainLeverage, SubIndex, [SubType], SubLeverage, StrikeInf, StrikeSup, Notional ",
            " 0",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a Spread Cap PayOff",
			" Main Index Name",
			" Main Index Type (INF, def: IR)",
			" Main Leverage Value",
			" Sub Index Name",
			" Sub Index Type (INF, def: IR)",
			" Sub  Leverage Value",
			" Strike Inf Value",
			" Strike Sup Value",
			" Notional Value",
		},
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		{
        	" Local_ARM_INF_GetAdjCorrel",				
            " RR",										
            " ARM_INF_GetAdjCorrel",
			" InfBsSmiledModel ",
            " 1",	
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " get the adj conv correlation",
			" Inf BS Smiled Model",
		},
		{
        	" Local_PXL_ARM_INF_GetAdjCorrel",
            " RR",									 
            " RR",										
            " PXL_ARM_INF_GetAdjCorrel",
			" InfBsSmiledModel ",
            " 0",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " get the adj conv correlation",
			" Inf BS Smiled Model",
		},
		{
        	" Local_ARM_INF_AdjConvexity",
            " RRRR",									 
            " ARM_INF_AdjConvexity",
			" InfIndex, YoY_t, YoY_T",
            " 1",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Compute the adjustment of convexity",
			" Main Index Name",
			" Main Index Type (INF, def: IR)",
			" Main Leverage Value",
			" Main Strike Value",
			" Sub Index Name",
			" Sub Index Type (INF, def: IR)",
			" Sub  Leverage Value",
			" Sub  Strike Value",
			" Spread Value",
			" Notional Value",
		},
		{
        	" Local_ARM_INF_EQHWSV_Laplace",				
            " RRRRRRRRRR",							
            " ARM_INF_EQHWSV_Laplace",
			" EQHWSVId,evalTime, startTime, endTime,[xt],[vt],[k_Re],[k_Im],[isReal]",
            " 1",										/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute E( exp( -k.x_{S}^{E} ) | x_{t}^{E}=xt, v_{t}=vt ) return its real or imaginary part",
			" Equity HW vol sto model",
			" t: evalTime",
			" S: startTime",
			" E: endTime",
			" initial value of xt(p), [ default xt=0.0 ]",
			" initial value of vt(),  [ default vt=1.0 ]",
			" k_Re: parameter	[ default k_Re = 0.0]",
			" k_Im: parameter	[ default k_Im = 1.0]",
			" isReal: Y or N [ default Y ]",
		},
		{
        	" Local_ARM_INF_EQHWSV_Density",				
            " RRRRRRRRRR",							
            " ARM_INF_EQHWSV_Density",
			" EQHWSVId,evalTime, startTime, endTime,[xt],[vt],[x],[Period],[frequency]",
            " 1",										/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute the density using the inverse Laplace",
			" Equity HW vol sto model",
			" t: evalTime",
			" S: startTime",
			" E: endTime",
			" initial value of xt(p)",
			" initial value of vt()",
			" x: parameter",
			" Period: Laplace pade frequency",
			" Freq: Laplace pade frequency",
		},