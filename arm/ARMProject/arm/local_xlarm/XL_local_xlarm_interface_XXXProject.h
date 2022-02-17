///	 ==> XXXProject
	{
        	" Local_MktDatas_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_XXX_MktDatas_Create",
			" AsOfDate,mktDatasKeys, mktDatasId",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a mkt datas object",
			" As of Date",
			" list of keys",
			" list of objects id",
    },
    {
        	" Local_PXL_MktDatas_Create",	/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " PXL_ARM_XXX_MktDatas_Create",
			" AsOfDate,mktDatasKeys, mktDatasId",
            " 0",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a mkt datas object",
			" As of Date",
			" list of keys",
			" list of objects id",
    },
	{
        	" Local_ARM_XXX_Price",		/// name of the C++ function
            " RRR",							/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_XXX_Price",
 			" security, MktDatas",
            " 1",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Prices a Vanilla Security with BS model", 
            " Security id", 
            " Market Datas",
    },	
	{
        	" Local_ARM_XXX_Hedge_Create",		/// name of the C++ function
            " RRRR",							/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_XXX_Hedge_Create",
			" Security, Scenario, MktDatas",
            " 1",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Compute the global hedge",
			" Security Id",
            " Scenario Id",
			" Market Datas Id",
	},
	{
        	" Local_ARM_XXX_Scenari_Compose",	
            " RRR",								
            " ARM_XXX_Scenari_Compose",
			" Scenario1, Scenario2",
            " 1",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Compose 2 distinct scenari",
			" Scenario Id1",
            " Scenario Id2",
	},
	{
			" Local_ARM_XXX_Hedge_GetData",
			" RRR",									/// 4 parametres = 3 d'entree + 1 parametre de retour 
            " ARM_XXX_Hedge_GetData",
			" HedgeId, Key",
            " 1",
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Get Result from hedging computation",
			" Hedge Id",
            " Key Scenario (DELTA_BS_EUR, VEGA_CAP_RHO_JPY )"
	},
	{
        	" Local_ARM_XXX_Scenario",			/// name of the C++ function
            " RRRRRRRRR",							/// 6 parametres = 5 d'entree + 1 parametre de retour 
            " ARM_XXX_Scenario",				/// Excel Name
			" Shift, Currency, Type Scenario, [SubType Scenario], [Stress Order], [Is Relative], [Is Backward], [Is Perturbative]",
            " 1",								/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
            " Creates a scenario of hedge",
			" Value of the shift",
			" Currency (EUR, USD, ...)",
			" Type of scenario (DELTA_ZC, DELTA_BS, VEGA_CAP, VEGA_OSW,...)",
			" Sub type of scenario (ATM, RHO, NU, BETA,...)",
			" Specification of the stress hedge ( T;E, T=3M;E=6M,2Y..)",
			" Specification Relative (Y) or Absolute (N) ( default N )",
			" Specification Backward (Y) or Forward behaviour (N)( default N )",
			" Specification Perturbative (Y) or Cumulative (N)( default N )"
    },
	
