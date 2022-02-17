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

	/// Generic Calibration Section
	{
        	" Local_CalibMethod_Create",			/// name of the C++ function
            " RRRRRRRRRR",							/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GP_CalibMethod_Create",
			" Portfolio Id, CalibParms (ids),CalibType,[Max_iter],[TargetFunType],[LinkCalMethod],[PrevCalMethod],[factorNb],[validate]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Calib method object",
			" string:Portfolio Id",
			" vector of string: CalibParms (Ids) Vector of CalibParms",
			" string: type of calibration (BOOTSTRAP1D, OPTIMIZE, BOOTSTRAPND, OPTIMIZE1D)",
			" integer: max_iterations (optional 100)",
			" string: whitch target used  (optional(PRICE_TAR),IMPVOL_TAR)",
			" string: Linked CalibMethod Id ( Optional(NULL))",
            " string: Previuos CalibMethod Id ( Optional(NULL))",
			" Long: for multiAsset, witch model do you want to calibrate",
			" boolean: to validate or not calib Method ( default true)",
    },
	{
        	" Local_PXL_CalibMethod_Create",	    /// name of the C++ function
            " RRRRRRRRRR",							/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CalibMethod_Create",
			" Portfolio Id, CalibParms (ids), CalibType, [Max_iter], [TargetFuncType],[LinkCalibMethod],[PrevCalibMethod],[factorNb],[validate]",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
    }, 
	{
        	" Local_PXL_CalibMethodWithDesc_Create",	    /// name of the C++ function
            " RRRRRRRRRR",							/// 10 parametres = 9 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CalibMethodWithDesc_Create",
			" MethodType, Portfolio Id, CalibParms (ids), ModelFitterDes Id,[TargetFuncType],[LinkCalibMethod],[PrevCalibMethod],[FactorNb],[validate]",
            " 0",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
           
    },  
	{
        	" Local_CalibMethodWithDesc_Create",	    /// name of the C++ function
            " RRRRRRRRRRR",							/// 11 parametres = 10 d'entree + 1 parametre de retour 
            " ARM_GP_CalibMethodWithDesc_Create",
			" MethodType, Portfolio Id,CalibParms (ids),ModelFitterDes Id,[TargetFunType],[LinkCalMethod],[PrevCalMethod],[FactorNb],[NbIteration],[validate]",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Calib Method object with a given Model Fitter Id",
			" string:type of calibration (BOOTSTRAP1D, OPTIMIZE, BOOTSTRAPND, OPTIMIZE1D)",
			" string: Portfolio Id",
			" vector of string: model Param ids",
			" string: Model Fitter Id",
			" string: whitch target used  (optional(PRICE_TAR),IMPVOL_TAR)",
			" string: Linked CalibMethod Id default(NULL)",
            " string: Previuos CalibMethod Id default(NULL)",
			" Long: for multiAsset, witch model do you want to calibrate",
			" long: to repeat same routine n time ( default n=1)",
			" boolean: to validate or not calib Method ( default true)",
    },
	{
        	" Local_CalibMethod2D_Create",	    /// name of the C++ function
            " RRRRRRRRRRR",							/// 11 parametres = 9 d'entree + 1 parametre de retour 
            " ARM_GP_CalibMethod2D_Create",
			" Portfolio1 Id,Portfolio2 Id, CalibParms1 (ids), CalibParms2 (ids), ModelFitterDes1 Id,ModelFitterDes2 Id, [TargetFuncType], [LinkCalibMethod], [PrevCalibMethod], [CalibDirection]",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Calib method2D object with description",
			" Portfolio1 Id",
			" Portfolio2 Id",
			" CalibParms1 (Ids) Vector of CalibParms",
			" CalibParms2 (Ids) Vector of CalibParms",
			" Model Fitter 1 Description (Id)",
			" Model Fitter 2 Description (Id)",
			" TargetFuncType (optional(PriceCalibration))",
			" LinkedCalibMethod Id   ( Optional(NULL))",
            " PreviuosCalibMethod Id ( Optional(NULL))",
			" Calib Direction ( Optional(Forward))"
    },
    {
        	" Local_GetDataFromCalibMethod",	    /// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_GetDataFromCalibMethod",
			" CalibMethod Id, DataType",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Get data from Calibration method Id",
			" CalibMethod Id",
			" Datatype (PORTFOLIO, LINKEDCALIBMETHOD, PREVIOUSCALIBMATHOD, CALIBMETHODTYPE, CALIBPARAMS)"
    },
    {
        	" Local_PXL_GetDataFromCalibMethod",	    /// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_GetDataFromCalibMethod",
			" CalibMethod Id, DataType",
            " 0",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Get data from Calibration method Id",
			" CalibMethod Id",
			" Datatype (PORTFOLIO, LINKEDCALIBMETHOD, PREVIOUSCALIBMATHOD, CALIBMETHODTYPE, CALIBPARAMS)"
    },
	
	{
        	" Local_SetDetailFlagToCalibMethod",	    /// name of the C++ function
            " RRR",							            /// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_SetDetailToCalibMethod",
			" CalibMethod Id,Detail flag",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Set flag for detail view on the calibration",
			" CalibMethod Id",
			" detail flag (true or false)"
    },
	{
        	" Local_GetDurationFromCalibMethod",	    /// name of the C++ function
            " RR",							            /// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_GetDurationFromCalibMethod",
			" CalibMethod Id",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " get duration of a calibration",
			" CalibMethod Id",
    },
    {
        	" Local_Calibrate",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_Calibrate",
			" Model Id,Calib Method Id",
            " 1",							/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " To calibrate a Model with a portfolio",
			" string: Pricing Model Id",
			" string: Calib Method Id",

    },
    {
        	" Local_PXL_Calibrate",	/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Calibrate",
			" Model Id,Calib Method Id",
            " 0",							/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_ARM_Optimizer_Create",		/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_Optimizer_Create",
			" AlgorithmType,[Max_Iter],[Tolerance],[Step_Max],[LocalSearch],[PrintLevel]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a NAG Optimizer",
			" string: Nag routine (BOUNDS_NO_DERIV, NLIN_LSQ, LSQ_DERIV, LSQ_CHECK_DERIV)",
			" long: number of iteration ( default 100)",
			" double: precision (default power(machine_precsion, 0.72))",
			" double: step Max  (default 2.0)",
			" boolean: local search (default false)",
			" boolean: to store Detail ( default false)",

    },
	{
        	" Local_PXL_Optimizer_Create",		/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Optimizer_Create",
			" algorithmType,[max_Iter],[xTolerance],[step_Max],[localSearch],[printLevel]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_PXL_ARM_Solver_Create",		/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_Solver_Create",
			" algorithmType,[max_Iter],[xTolerance],[fxTolerance],[gradTolerance],[printLevel]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_ARM_Solver_Create",		/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_Solver_Create",
			" algorithmType,[max_Iter],[xTolerance],[fxTolerance],[gradTolerance],[printLevel]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a mono-dim Solver. Attention!!!, the last selection is a NAG optimizer mono_dim",
			" string:(STD_NR,NR_SMOOTHING,NR_RETRIAL,NR_NOTHROW,NAG_SOLVER,DICHOTOMY)",
			" integer:number of iteration ( default 100)",
			" float: (default power(machine_precsion, 0.72))",
			" float: (default 1.0e-14)",
			" float: (default 1.0e-10)",
			" boolean:to get Detail ( default false)",
    },
	{
        	" Local_DensityFunctor_CallOption",		/// name of the C++ function
            " RRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctor_CallOption",
			" functor, forward, strike, maturity",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Price of std call option", 
			" "
			" "
    },
	{
        	" Local_DensityFunctor_Quantile",		/// name of the C++ function
            " RRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctor_Quantile",
			" functor, forward, proba, maturity",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Quantile", 
			" "
			" "
    },
	{/// Density functors for Markov Functional Calibration
        	" Local_PXL_SLNDensityFunctor",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorSLN_Create",
			" Volatility, Shift",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_SLNDensityFunctor",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorSLN_Create",
			" Volatility, Shift",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a shifted Log Normal density functor", 
			" Volatility (20% = 0.2)"
			" Shift "
    },
	{/// Density functors for Markov Functional Calibration
        	" Local_PXL_MixtureDensityFunctor",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorMixture_Create",
			" Volatility1, Volatility2, Alpha, Lambda",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
			" ",
            " ",
            " Creates a Mixture density functor", 
			" Volatility1 (20% = 0.2)",
			" Volatility2 (20% = 0.2)",
			" Alpha ",
			" Lambda "
    },
	{
        	" Local_MixtureDensityFunctor",		/// name of the C++ function
            " RRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorMixture_Create",
			" Volatility1, Volatility2, Alpha, Lambda",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Mixture density functor", 
			" Volatility1 (20% = 0.2)",
			" Volatility2 (20% = 0.2)",
			" Alpha ",
			" Lambda "
    },
	{/// Density functors for Markov Functional Calibration
        	" Local_PXL_MixtureDensityFunctorWithATMVol",		/// name of the C++ function
            " RRRRRRR",						/// 7 parametres = 4 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorMixtureWithATMVol_Create",
			" Fwd, Maturity, VolATM, DecVol, Alpha, Lambda",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
			" ",
            " ",
            " Creates a Mixture density functor with ATM volatility", 
			" Fwd",
			" Maturity",
			" VolATM (20% = 0.2)",
			" DecVol",
			" Alpha ",
			" Lambda "
    },
	{
        	" Local_MixtureDensityFunctorWithATMVol",		/// name of the C++ function
            " RRRRRRR",						/// 5 parametres = 4 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorMixtureWithATMVol_Create",
			" Fwd, Maturity, VolATM, DecVol, Alpha, Lambda",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Mixture density functor with ATM volatility", 
			" Fwd",
			" Maturity",
			" VolATM (20% = 0.2)",
			" DecVol",
			" Alpha ",
			" Lambda "
    },
	{/// Density functors for Markov Functional Calibration
        	" Local_PXL_HestonDensityFunctor",		/// name of the C++ function
            " RRRRRRRRR",							/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_HestonDensityFunctor_Create",
			" V0, Kappa, Theta, VVol, Rho, Shift, Level, [Sigma]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
			" ",
            " ",
            " Creates a Heston density functor.", 
			" V0",
			" Kappa",
			" Theta",
			" VVol",
			" Rho",
			" Shift",
			" Level",
			" Sigma",
    },
	{
        	" Local_HestonDensityFunctor",		/// name of the C++ function
            " RRRRRRRRR",						/// 9 parametres = 8 d'entree + 1 parametre de retour 
            " ARM_GP_HestonDensityFunctor_Create",
			" V0, Kappa, Theta, VVol, Rho, Shift, Level, [Sigma]",
            " 1",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
			" ",
            " ",
            " Creates a Heston density functor.", 
			" V0",
			" Kappa",
			" Theta",
			" VVol",
			" Rho",
			" Shift",
			" Level",
			" Sigma"
    },
	{/// Density functors for Markov Functional Calibration
        	" Local_PXL_SplineDensityFunctor",		/// name of the C++ function
            " RRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorSpline_Create",
			" moneyness, vol, voltype, smile",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_SplineDensityFunctor",		/// name of the C++ function
            " RRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorSpline_Create",
			" moneyness, vol, voltype, [smile]",
           " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Takes either the triplet (moneyness,vol,voltype) or smile", 
			" moneyness",
			" vol",
			" vol type (GAUSS/BLACK), default=GAUSS",
			" smile viewer id",
    },
	{
        	" Local_PXL_SABRDensityFunctor",		/// name of the C++ function
            " RRRRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorSABR_Create",
			" Alpha,Beta,Rho,Nu,[SabrType],[GridSize]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_SABRDensityFunctor",		/// name of the C++ function
            " RRRRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorSABR_Create",
			" Alpha,Beta,Rho,Nu,[SabrType],[GridSize]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a SABR density functor", 
			" SABR Initial Vol",
			" SABR Backbone",
			" SABR Correlation",
			" SABR Vol of Vol",
			" SABR_IMPLNVOL, SABR_IMPLNVOL2, SABR_A, SABR_G, DIRECTEXACT, DIRECTGEOMETRIC, (Default = SABR_IMPLNVOL)",
			" Size of Precomputed Grid of Quantiles",
    },
	{
        	" Local_PXL_BiSABRDensityFunctor",		/// name of the C++ function
            " RRRRRRRRRRRRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorBiSABR_Create",
			" Alpha1,Beta1,Rho1,Nu1,Alpha2,Beta2,Rho2,Nu2,RhoS1S2,RhoS1V2,RhoS2V1,RhoV1V2,[SabrType],[GridSize]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_BiSABRDensityFunctor",		/// name of the C++ function
            " RRRRRRRRRRRRRRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_GP_DensityFunctorBiSABR_Create",
			" Alpha1,Beta1,Rho1,Nu1,Alpha2,Beta2,Rho2,Nu2,RhoS1S2,RhoS1V2,RhoS2V1,RhoV1V2,[SabrType],[GridSize]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a SABR density functor", 
			" SABR Initial Vol (1st Rate)",
			" SABR Backbone (1st Rate)",
			" SABR Correlation (1st Rate)",
			" SABR Vol of Vol (1st Rate)",
			" SABR Initial Vol (2nd Rate)",
			" SABR Backbone (2nd Rate)",
			" SABR Correlation (2nd Rate)",
			" SABR Vol of Vol (2nd Rate)",
			" BiSABR Rho S1/S2",
			" BiSABR Rho S1/V2",
			" BiSABR Rho S2/V1",
			" BiSABR Rho V1/V2",
			" SABR_IMPLNVOL, SABR_IMPLNVOL2, SABR_A, SABR_G, DIRECTEXACT, DIRECTGEOMETRIC, (Default = SABR_IMPLNVOL)",
			" Size of Precomputed Grid of Quantiles",
    },
	{
			" Local_NormalHestonDensityFunctor",
			" RRRRRRRR",
			" ARM_GP_NormalHestonDensityFunctor_Create",
			" Forward, V0, Kappa, Theta, VVol, Rho, [Level]",
			" 1",
			XLLOCALARM_GENCALIB_GROUP,
			" ",
			" ",
			" Creates a Normal Heston density functor (for spread)",
			" Fwd Spread",
			" Normal Heston Initial Var",
			" Normal Heston Var Mean Reversion",
			" Normal Heston Long Term Var", 
			" Normal Heston Vol of Var",
			" Normal Heston Correlation",
			" Normal Heston Level of Var",
	},
	{
			" Local_PXL_NormalHestonDensityFunctor",
			" RRRRRRRR",
			" PXL_ARM_GP_NormalHestonDensityFunctor_Create",
			" Forward, V0, Kappa, Theta, VVol, Rho, [Level]",
			" 0",
			XLLOCALARM_GENCALIB_GROUP,
			" ",
			" ",
			" Creates a Normal Heston density functor (for spread)",
			" Fwd Spread",
			" Normal Heston Initial Var",
			" Normal Heston Var Mean Reversion",
			" Normal Heston Long Term Var", 
			" Normal Heston Vol of Var",
			" Normal Heston Correlation",
			" Normal Heston Level of Var",
	},
	{
        	" Local_PXL_IrFwdDensityFunctor",		/// name of the C++ function
            " R",						/// 1 parametres = 1 parametre de retour 
            " PXL_ARM_GP_DensityFunctorNoVol_Create",
			" ",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_IrFwdDensityFunctor",		/// name of the C++ function
            " R",						/// 1 parametres = 1 parametre de retour 
            " ARM_GP_DensityFunctorNoVol_Create",
			" ",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates an IrFwd density functor", 
    },
	{
        	" Local_PXL_VanillaSecurityDensity",		/// name of the C++ function
            " RRRRRRRRRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_VanillaSecurityDensity_Create",
			" ResetDate,StartDate,EndDate,DensityFunctor,[Frequency],[DayCount],[StubRule], [Weight],[FwdAdd],[FwdMult]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_VanillaSecurityDensity",		/// name of the C++ function
            " RRRRRRRRRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_VanillaSecurityDensity_Create",
			" ResetDate,StartDate,EndDate,DensityFunctor,[Frequency],[DayCount],[StubRule], [Weight],[FwdAdd],[FwdMult]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Libor / SwapRate description (Hunt Kennedy calibration)", 
			" Rate Reset Date",
			" Rate Start Date",
			" Rate End Date",
			" Density Functor Id",
			" Frequency (default = calib schedule frequency)",
			" Day count convention (default = calib schedule day count convention)",
			" Stub rule (default = short start)",
			" Weight",
			" add adjustment to forward (default = 0)",
			" mult adjustment to forward (default = 1)"

    },
	{
        	" Local_PXL_VanillaSecurityDensitySpread",		/// name of the C++ function
            " RRRRRRRRRRRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_VanillaSecurityDensitySpread_Create",
			" ResetDate,StartDate1,EndDate1,StartDate2,EndDate2,DensityFunctor,[Frequency1],[DayCount1],[Frequency2],[DayCount2],[StubRule],[Weight]",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
    },
	{
        	" Local_VanillaSecurityDensitySpread",		/// name of the C++ function
            " RRRRRRRRRRRRR",						/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_GP_VanillaSecurityDensitySpread_Create",
			" ResetDate,StartDate1,EndDate1,StartDate2,EndDate2,DensityFunctor,[Frequency1],[DayCount1],[Frequency2],[DayCount2],[StubRule], [Weight]",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Libor / SwapRate description (Hunt Kennedy calibration)", 
			" Rate Reset Date",
			" Rate Start Date (1st rate)",
			" Rate End Date (1st rate)",
			" Rate Start Date (2nd rate)",
			" Rate End Date (2nd rate)",
			" Density Functor Id",
			" Frequency (1st rate) (default = calib schedule frequency)",
			" Day count convention (1st rate) (default = calib schedule day count convention)",
			" Frequency (2nd rate) (default = calib schedule frequency)",
			" Day count convention (2nd rate) (default = calib schedule day count convention)",
			" Stub rule (default = short start)",
			" Weight",

    },
	{
        	" Local_PXL_VanillaSecurityDensityFX",		/// name of the C++ function
            " RRRRRR",									/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_VanillaSecurityDensityFX_Create",
			" ResetDate,DensityFunctor,DomCurve,ForCurve,FXSpot",
            " 0",								/// not visible in excel
            XLLOCALARM_GENCALIB_GROUP,
			" ",
            " ",
            " Smiled FX calibtation data", 
			" Rate Reset Date",
			" Density Functor Id",
			" DomCurve",
			" ForCurve",
			" FXSpot"
    },
	{
        	" Local_VanillaSecurityDensityFX",		/// name of the C++ function
			" RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_VanillaSecurityDensityFX_Create",
			" ResetDate,DensityFunctor,DomCurve,ForCurve,FXSpot",
            " 1",								/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Smiled FX calibtation data", 
			" Rate Reset Date",
			" Density Functor Id",
			" DomCurve",
			" ForCurve",
			" FXSpot"
    },
	{ /// Numerical CalibMethod
        	" Local_NumericalCalibMethod_Create",			/// name of the C++ function
            " RRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_GP_NumericalCalibMethod_Create",
			" CalibDateStrip,Securities,[Portfolio]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Numerical Calib method object",
			" DateStrip Id",
			" Vanilla Security + Densities description",
			" Portfolio Id (bootstrap calibration)",
    },
	{
        	" Local_PXL_NumericalCalibMethod_Create",	    /// name of the C++ function
            " RRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_NumericalCalibMethod_Create",
			" CalibDateStrip,Securities,[Portfolio]",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
    },
	{ /// HW2F CalibMethod
        	" Local_CalibMethodHW2F_Create",			/// name of the C++ function
            " RRRRRRRR",							/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_CalibMethodHW2F_Create",
			" Portfolio1,Param1,[Portfolio2],[Param2],[Portfolio3],[Param3],[WithOptim]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Creates a Calib Method object dedicated to HW2F",
			" ",
			" ",
			" ",
			" ",
			" ",
			" ",
			" WithOptim? (0=false, 1=yes)",
    },
	{
        	" Local_PXL_CalibMethodHW2F_Create",	    /// name of the C++ function
            " RRRRRRRR",						/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_NumericalCalibMethod_Create",
			" Portfolio1,Param1,[Portfolio2],[Param2],[Portfolio3],[Param3],[WithOptim]",
           " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
			" ",
			" ",
			" ",
			" ",
			" ",
			" ",
			" ",
			" WithOptim? (0=false, 1=yes)",
    },
	{ /// Basket
        	" Local_BasketDecomp_Create",			/// name of the C++ function
            " RRRRRRRRRRR",							
            " ARM_GP_BasketCalibration_Create",
			" SecuritiesId,ModelsId,Weights,Notional,ExerDatestrip,Fees,[Side],PricingModel,[Method],[Strike]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Computes Basket Decomposition of a Security",
			" Vector of securities",
			" Vector of bs models (to price each security)",
			" Vector of weights",
			" Notional (numeric of reference value)",
			" Exercise Datestrip",
			" Call fees (numeric of reference value)",
			" Side (P/R)",
			" Pricing model (for VN swaptions)",
			" Method (BASKET / BASKET_SIMPLE) (default BASKET)",
			" Strike (EQUIVALENT / ATM) (default EQUIVALENT)",
    },
	{
        	" Local_PXL_BasketDecomp_Create",	    /// name of the C++ function
            " RRRRRRRRR",						
            " PXL_ARM_GP_BasketCalibration_Create",
			" SecuritiesId,ModelsId,Weights,Notional,ExerDatestrip,Fees,[Side],PricingModel,[Method],[Strike]",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
    },
	{ 
        	" Local_GetData_FromBasket",			/// name of the C++ function
            " RRR",							
            " ARM_GP_BasketCalibration_GetData",
			" BasketId,Key",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Get datas from basket",
			" Basket Id",
			" Key (PORTFOLIO/BASKET)",
    },
	{
        	" Local_PXL_GetData_FromBasket",	    /// name of the C++ function
            " RRR",						
            " PXL_ARM_GP_BasketCalibration_GetData",
			" BasketId,Key",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
    },
	{ /// Smile Viewer
        	" Local_SmileViewer_Create",			/// name of the C++ function
            " RRRRRR",							
            " ARM_GP_SmileViewer_Create",
			" SecurityId,ModelId,Moneyness,[MoneyType],[ExtraStrikes]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Computes Smile Viewer of a Security",
			" Security",
			" Model",
			" Vector of moneyness",
			" Moneyness Type (GAUSS/BLACK), default=GAUSS",
			" Extra Strikes, default=empty",
			" Tol",
    },
	{
        	" Local_PXL_SmileViewer_Create",	    /// name of the C++ function
            " RRRRRR",						
            " PXL_ARM_GP_SmileViewer_Create",
			" SecurityId,ModelId,Moneyness,[MoneyType],[ExtraStrikes]",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
    },
	{ /// Smile Viewer
        	" Local_DensityFunctorGen_Create",			/// name of the C++ function
            " RRRRRRR",							
            " ARM_GP_DensityFunctorGen_Create",
			" SecurityId,ModelId,DecStrike,[IsDirect],[MinProba],[MaxProba]",
            " 1",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Computes Implied Density Functor of a Security priced with a model",
			" Security",
			" Model",
			" Dec Strike (for digital)",
			" Direct (Y/N, default N)",
			" Min Proba (default 0.00001)",
			" Max Proba (default 0.99999)",
    },
	{ /// Smile Viewer
        	" Local_DensityFunctorGen_Create",			/// name of the C++ function
            " RRRRRRR",							
            " PXL_ARM_GP_DensityFunctorGen_Create",
			" SecurityId,ModelId,DecStrike,[IsDirect],[MinProba],[MaxProba]",
            " 0",									/// visible in excel
            XLLOCALARM_GENCALIB_GROUP,
            " ",
            " ",
            " Computes Implied Density Functor of a Security priced with a model",
			" Security",
			" Model",
			" Dec Strike (for digital)",
			" Direct (Y/N, default N)",
			" Min Proba (default 0.00001)",
			" Max Proba (default 0.99999)",
    },
