/// ==> MERCURE (5 functions)
{
	" Mercure_Hedge",					// Nom de la fonction qu'on va appeler
    " RRRRRRRRRRRR",					// 11 parametres d'entree + 1 parametre de retour
    " Mercure_Hedge",					// Nom qui apparaît dans le wizard
    " TradeId,TradeType,AsOfDate,FallBack,[SensiType],[MarketDataManager],[HedgeRatioFilesDirectory],[ConfigFile],[MarketDataFile],[Trace],[HedgeRatioFile]",// Paramètres
    " 1",								// 1 = visible sur XL, 0 = non visible
    XLLOCALARM_MERCURE_GROUP,			// Nom du groupe dans lequel apparaîtra la fonction
    " ",
    " ",								// Chemin du fichier d'aide
    " Computes sensis for a deal",		// Description générale
    " Summit id (or ARM object)",		// Description du 1er param
    " Must be empty if trade Id is an ARM object",	// Description du 2ème param
	" As Of Date",
	" FallBack",
	" Corresponds to the eventual subdirectory where the ratio file is stored",
	" To add other market data into dictionnary",
	" HedgeRatioFiles Directory",
	" Config File",
	" MarketData File",
	" flag to get a view of the security and the model : YES or NO (default value)",
	" Take this file instead of the default one"
},
{
	" Local_View_XML",
    " RR",								// 1 parametre d'entree + 1 parametre de retour
    " ARM_View_XML",
    " Mercure result id",
    " 1",
    XLLOCALARM_MERCURE_GROUP,
    " ",
    " ",
    " Views Mercure result XML file", 
    " Mercure result (LMRCR) id"
},
{
	" Mercure_Hedge_Array",				// Nom de la fonction qu'on va appeler
    " RRRRRRRRRRRR",					// 11 parametres d'entree + 1 parametre de retour
    " Mercure_Hedge_Array",				// Nom qui apparaît dans le wizard
    " TradeId,TradeType,AsOfDate,FallBack,[SensiType],[MarketDataManager],[HedgeRatioFilesDirectory],[ConfigFile],[MarketDataFile],[Trace],[HedgeRatioFile]",// Paramètres
    " 1",								// 1 = visible sur XL, 0 = non visible
    XLLOCALARM_MERCURE_GROUP,			// Nom du groupe dans lequel apparaîtra la fonction
    " ",
    " ",								// Chemin du fichier d'aide
    " Computes sensis for a deal",		// Description générale
    " Summit id (or ARM object)",		// Description du 1er param
    " Must be empty if trade Id is an ARM object",	// Description du 2ème param
	" As Of Date",
	" FallBack",
	" Corresponds to the eventual subdirectory where the ratio file is stored",
	" To add other market data into dictionnary",
	" HedgeRatioFiles Directory",
	" Config File",
	" MarketData File",
	" flag to get a view of the security and the model : YES or NO (default value)",
	" Take this file instead of the default one"
},
{
	" Mercure_ARM_Hedge",				// Nom de la fonction qu'on va appeler
    " RRRRRRRRRRR",						// 10 parametres d'entree + 1 parametre de retour
    " Mercure_ARM_Hedge",				// Nom qui apparaît dans le wizard
    " SecurityId,AsOfDate,FallBack,[SensiType],[MarketDataManager],[HedgeRatioFilesDirectory],[ConfigFile],[MarketDataFile],[Trace],[HedgeRatioFile]",// Paramètres
    " 1",								// 1 = visible sur XL, 0 = non visible
    XLLOCALARM_MERCURE_GROUP,			// Nom du groupe dans lequel apparaîtra la fonction
    " ",
    " ",								// Chemin du fichier d'aide
    " Computes sensis for an ARM security",	// Description générale
    " Security",						// Description du 1er param
	" As Of Date",						// Description du 2ème param
	" FallBack",
	" Corresponds to the eventual subdirectory where the ratio file is stored",
	" To add other market data into dictionnary",
	" HedgeRatioFiles Directory",
	" Config File",
	" MarketData File",
	" flag to get a view of the security and the model : YES or NO (default value)",
	" Take this file instead of the default one"
},
{
	" CreateMarketDataManager",			// Nom de la fonction qu'on va appeler
    " RRRRR",							// 4 parametres d'entree + 1 parametre de retour
    " Mercure_CreateMarketDataManager",	// Nom qui apparaît dans le wizard
    " MarketDataIds,AsOfDate,[FallBack],[SwitchToETK]",// Paramètres
    " 1",								// 1 = visible sur XL, 0 = non visible
    XLLOCALARM_MERCURE_GROUP,			// Nom du groupe dans lequel apparaîtra la fonction
    " ",
    " ",								// Chemin du fichier d'aide
    " Create MarketDataManager",		// Description générale
    " MarketData Ids",					// Description du 1er param
	" As Of Date",						// Description du 2ème param
	" FallBack",
	" SwitchToETK : YES or NO (default value)"
},
{
	" CreateARMScalarData",				// Nom de la fonction qu'on va appeler
    " RRRRRRRR",						// 7 parametres d'entree + 1 parametre de retour
    " Mercure_CreateARMScalarData",		// Nom qui apparaît dans le wizard
    " Value,Type,Currency,Index,AsOfDate,[CurveId],[SwitchToETK]",	// Paramètres
    " 1",								// 1 = visible sur XL, 0 = non visible
    XLLOCALARM_MERCURE_GROUP,			// Nom du groupe dans lequel apparaîtra la fonction
    " ",
    " ",								// Chemin du fichier d'aide
    " Create ARMScalarData",			// Description générale
    " Value",							// Description du 1er param
	" Type : FX, CUTOFF 2F, CUTOFF 3F, MEANREV 2F, MEANREV 3F, QMOD0, QMOD1, FXIR CORR 2F, FXIR CORR 3F, IRIR CORR 2F, IRIR CORR 3F" // Description du 2ème param
	" Currency",
	" Index : LIBOR, EURIB, ...",
	" As Of Data"
	" Curve Id : MO (default value), MOSMILE, MO3F, ...",
	" SwitchToETK : YES or NO (default value)"
},
{
	" Mercure_GetPostProcessedData",
    " RRRR",	// 3 parametres d'entree + 1 parametre de retour
    " Mercure_GetPostProcessedData",
    " MercureResultObject, DataName, [PlotName]",
    " 1",
    XLLOCALARM_MERCURE_GROUP,
    " ",
    " ",
    " Extract post processing value from a hedge result object",
    " Mercure result object (LMRCR)",
	" Data name (PV, or market data name like 'Zero_Curve/EUR/EURIB')",
	" Plot ('2D', '1Y'...) or underlying when PV (like 'Funding', 'Option'...). Depends on processor name in config file."
},
{
	" Mercure_ViewModelParams",
    " RR",	// 1 parametre d'entree + 1 parametre de retour
    " Mercure_ViewModelParams",
    " MetaModelName",
    " 1",
    XLLOCALARM_MERCURE_GROUP,
    " ",
    " ",
    " Create an object that displays the extra model parameters for a given MetaModel",
	" 'BS', 'BSSmiled', 'BSGen', 'CRACalculator', '2IRFXModel', 'Multi3F', 'QModel', 'Mixture', 'FXOption3F' (empty to display all models)"
},
/// END MERCURE