		/// ==> UTILITIES (91 functions)
		{
				" Local_ISOCCY",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ISOCCY",
                " Name",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an ISO currency object",
                " Currency ISO name"
        },
		{
				" Local_ISOCCY",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_ISOCCY",
                " Name",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ISOCCY,
                " Creates an ISO currency object",
                " Currency ISO name"
        },
		{
				" Local_PXL_ISOCCY",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " PXL_ISOCCY",
                " Name",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an ISO currency object",
                " Currency ISO name"
        },
		{
				" Local_PXL_ISOCCY",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " PXL_ARM_ISOCCY",
                " Name",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an ISO currency object",
                " Currency ISO name"
        },
		{
				" Local_CCY",
                " RRRRR",					// 4 parametres d'entree + 1 parametre de retour
                " CCY",
                " Name,CurveId,CrossValue,[DayCount]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a currency object",
                " Currency name",
                " Curve id",
				" Cross value",
				" Day count (default : default daycount of ccy name)"
        },
		{
				" Local_CCY",
                " RRRRR",					// 4 parametres d'entree + 1 parametre de retour
                " ARM_CCY",
                " Name,CurveId,CrossValue,[DayCount]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_CCY,
                " Creates a currency object",
                " Currency name",
                " Curve id",
				" Cross value",
				" Day count (default : default daycount of ccy name)"
        },
		{
				" Local_PXL_CCY",
                " RRRRR",					// 4 parametres d'entree + 1 parametre de retour
                " PXL_CCY",
                " Name,CurveId,CrossValue,[DayCount]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a currency object",
                " Currency name",
                " Curve id",
				" Cross value",
				" Day count (default : default daycount of ccy name)"
        },
		{
				" Local_PXL_CCY",
                " RRRRR",					// 4 parametres d'entree + 1 parametre de retour
                " PXL_ARM_CCY",
                " Name,CurveId,CrossValue,[DayCount]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a currency object",
                " Currency name",
                " Curve id",
				" Cross value",
				" Day count (default : default daycount of ccy name)"
        },
		{
				" Local_ARM_GetInfoFromCcy",
                " RRR",					// 2 parametres d'entree + 1 parametre de retour
                " ARM_GetInfoFromCcy",
                " CcyId,type",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets Information from currency object",
                " Currency object Id",
                " Type (RESETCAL, PAYCAL)"
        },
		{
				" Local_ARM_Price",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_Price",
                " Security id,Model id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_Price,
                " Prices a security with the given model", 
                " Security id", 
                " Model id"
        },
		{
				" Local_bsflexible",
                " RRRRRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_BSFlexible",
                " Forward,Total Volatility,Bond Price,Strike,Call=1 or Put=-1",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ", //IDH_ARM_BSFlexible,
                " Return a Black and Scholes Price", 
                " Forward Price of the Underlying at Reset Time", 
                " Volatility of the underlying at reset time",
				" Price of Bond maturing at payment time",
				" Strike Price",
				" Call = 1 or Put = -1"
        },
		{
				" Local_FreeObject",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " FreeObject",
                " Security id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Frees an object on ARM server", 
                " Security id"
        },
		{
				" Local_FreeObject",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_FreeObject",
                " Security id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FreeObject,
                " Free an object on ARM server", 
                " Security id"
        },
		{
				" Local_FreeAllObjects",
                " R",						// 0 parametre d'entree + 1 parametre de retour
                " FreeAllObjects",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Frees all objects on ARM server"
        },
		{
				" Local_FreeAllObjects",
                " R",						// 0 parametre d'entree + 1 parametre de retour
                " ARM_FreeAllObjects",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FreeAllObjects,
                " Frees all objects on ARM server"
        },
		{
				" Local_NextBusinessDay",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " NextBusinessDay",
                " Date,[Currency],[Days]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gives the next business day", 
                " Date", 
                " Currency",
				" Number of days (default: 1)"
        },
		{
				" Local_NextBusinessDay",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_NextBusinessDay",
                " Date,[Currency],[Days]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_NextBusinessDay,
                " Gives the next business day", 
                " Date", 
                " Currency",
				" Number of days (default: 1)"
        },
		{
				" Local_IsBusinessDay",
                " RRR",							// 2 parametres d'entree + 1 parametre de retour
                " IsBusinessDay",
                " Date,[Currency]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Tests if a day is a business one", 
                " Date", 
                " Currency"
        },
		{
				" Local_IsBusinessDay",
                " RRR",							// 2 parametres d'entree + 1 parametre de retour
                " ARM_IsBusinessDay",
                " Date,[Currency]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_IsBusinessDay,
                " Tests if a day is a business one", 
                " Date", 
                " Currency"
        },
		{
				" Local_ARM_FutDelivery",
                " RRR",							// 2 parametres d'entree + 1 parametre de retour
                " ARM_FutDelivery",
                " Future,Currency",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gives the delivery date of the Future Contract", 
                " Future (MAR06, DEC07...)",
				" Currency (EUR, USD...)"
        },
		{
				" Local_ARM_Accrued",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_Accrued",
                " Security id, a Date,[Model id]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_Accrued,
                " Returns accrued of a security with the given model", 
                " Security id", 
                " a Date",
				" Model id"
        },
		{
				" Local_IMPLIEDVOL",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " IMPLIEDVOL",
                " Security id,Model id,Price,[TypeVol]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Returns implied volatility of a security with the given model", 
                " Security id", 
                " Model id",
				" Price",
				" Volatility type LnNor or Nor (Y/N) (default LnNor:Y)"
        },
		{
				" Local_IMPLIEDVOL",
                " RRRRR",						// 4  parametres d'entree + 1 parametre de retour
                " ARM_IMPLIEDVOL",
                " Security id,Model id,Price,[TypeVol]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_IMPLIEDVOL,
                " Returns implied volatility of a security with the given model", 
                " Security id", 
                " Model id",
				" Price",
				" Volatility type LnNor or Nor (Y/N) (default LnNor:Y)"
        },
		{
				" Local_ARM_View",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_View",
                " Security id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_View,
                " Views an instrument", 
                " Instrument id"
        },
		{
				" Local_SetNotional",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " SetNotional",
                " Security id,RefVal Id,[percentRemainder],[interpolationDates]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Modifies the security notional", 
                " Security id", 
                " Reference value id",
				" percentage of Remainder (default 100.)",
				" Dates on which notional will be interpolated (start dates:0, end dates:1, reset dates:2, payment dates:3)"
        },
		{
				" Local_SetNotional",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_SetNotional",
                " Security id,RefVal Id,[percentRemainder],[interpolationDates]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_SetNotional,
                " Modifies the security notional", 
                " Security id", 
                " Reference value id",
				" percentage of Remainder (default 100.)",
				" Dates on which notional will be interpolated (start dates:0, end dates:1, reset dates:2, payment dates:3)"
        },
		{
				" Local_ARM_today",
                " R",							// 0 parametre d'entree + 1 parametre de retour
                " ARM_today",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_today,
                " Returns the date of the current day"
        },
		{
				" Local_GetExpiry",
                " RR",							// 1 parametre d'entree + 1 parametre de retour
                " GetExpiry",
                " Security id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Returns the security expiry", 
                " Security id"
        },
		{
				" Local_GetExpiry",
                " RR",							// 1 parametre d'entree + 1 parametre de retour
                " ARM_GetExpiry",
                " Security id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_GetExpiry,
                " Returns the security expiry", 
                " Security id"
        },
		{
				" Local_Sensitivity",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " Sensitivity",
                " Security id,Model id,Param",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Returns the sensitivity of a security with the given model", 
                " Security id", 
                " Model id",
				" Sensitivity parameter"
        },
		{
				" Local_Sensitivity",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_Sensitivity",
                " Security id,Model id,Param",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_Sensitivity,
                " Returns the sensitivity of a security with the given model", 
                " Security id", 
                " Model id",
				" Sensitivity parameter"
        },
		{
				" Local_ADJUSTTOBUSDATE",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ADJUSTTOBUSDATE",
                " Date,Currency,Rule",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Adjusts to business date", 
                " Date", 
                " Currency",
				" Rule (previous or forward)"
        },
		{
				" Local_ADJUSTTOBUSDATE",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_ADJUSTTOBUSDATE",
                " Date,Currency,Rule",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ADJUSTTOBUSDATE,
                " Adjusts to business date", 
                " Date", 
                " Currency",
				" Rule (previous or forward)"
        },
		{
				" Local_FwdPrice",
                " RRRR",							// 3 parametres d'entree + 1 parametre de retour
                " FwdPrice",
                " Security id,Model id,Forward date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Computes a forward price of a security with the given model", 
                " Security id", 
                " Model id",
				" Forward date"
        },
		{
				" Local_FwdPrice",
                " RRRR",							// 3 parametres d'entree + 1 parametre de retour
                " ARM_FwdPrice",
                " Security id,Model id,Forward date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FwdPrice,
                " Computes a forward price of a security with the given model", 
                " Security id", 
                " Model id",
				" Forward date"
        },
		{
				" Local_CvSensitivity",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " CvSensitivity",
                " Security id,Model id,Param",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Returns the security sensitivity with the given model", 
                " Security id", 
                " Model id",
				" Sensitivity parameter"
		},
		{
				" Local_CvSensitivity",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_CvSensitivity",
                " Security id,Model id,Param",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_CvSensitivity,
                " Returns the security sensitivity with the given model", 
                " Security id", 
                " Model id",
				" Sensitivity parameter"
		},
		{
				" Local_ARM_BetweenDates",
                " RRRRR",						// 4 parametre d'entree + 1 parametre de retour
                " ARM_BetweenDates",
                " date1,date2,daycount,[isYearFrac]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_BetweenDates,
                " Calculates nb of days or year fraction between 2 dates", 
                " date 1", 
                " date 2",
				" daycount basis: ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, NOBASE",
				" isYearFrac: if set to 1 the result is a year fraction, else a number of days"
        },
		{
				" Local_ARM_CountBusinessDays",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_CountBusinessDays",
                " date1,date2,calendar",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Calculates nb of business days between 2 dates", 
                " date 1", 
                " date 2",
				" calendar",
        },
		{
				" Local_ARM_ADDMONTHS",
                " RRRRR",						// 4 parametre d'entree + 1 parametre de retour
                " ARM_ADDMONTHS",
                " date,nbmonths,[adjRule],[ccy]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ADDMONTHS,
                " Adds nbmonths months to the date ", 
                " date", 
                " nbmonths",
				" adjustment rule (default NONE)",
				" currency"
        },
		{
				" Local_ARM_ADDYEARS",
                " RRRRR",						// 4 parametre d'entree + 1 parametre de retour
                " ARM_ADDYEARS",
                " date,nbyears,[adjRule],[ccy]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ADDYEARS,
                " Adds nbyears years to the date ", 
                " date", 
                " nbyears",
				" adjustment rule (default NONE)",
				" currency"
        },
		{
				" Local_FxConvert",
                " RRRRRR",						// 5 parametres d'entree + 1 parametre de retour
                " FxConvert",
                " Currency 1,Currency 2,As of date,[Amount],[cvname]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Forex conversion", 
                " Currency 1", 
                " Currency 2",
				" As of date",
				" Amount (default: 1)",
				" Curve name (default MO)"
        },
		{
				" Local_FxConvertFromCalypso",
                " RRRRR",						// 5 parametres d'entree + 1 parametre de retour
                " FxConvertFromCalypso",
                " Currency 1,Currency 2,As of date,[cvname]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Forex conversion from Calypso", 
                " Currency 1", 
                " Currency 2",
				" As of date",
				" Curve name (default MO)"
        },
		{
				" Local_ARM_GetCurrency",
                " RR",						// 1 parametres d'entree + 1 parametre de retour
                " ARM_GetCurrency",
                " Security",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get currency name from Security", 
                " Security"
        },
		{
				" Local_FxConvert",
                " RRRRRR",						// 5 parametres d'entree + 1 parametre de retour
                " ARM_FxConvert",
                " Currency 1,Currency 2,As of date,[Amount],[cvname]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FxConvert,
                " Forex conversion", 
                " Currency 1", 
                " Currency 2",
				" As of date",
				" Amount (default: 1)",
				" Curve name (default MO)"
        },
		{
				" Local_ARM_GetDefaultCurrency",
                " R",						// 0 parametre d'entree + 1 parametre de retour
                " ARM_GetDefaultCurrency",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_GetDefaultCurrency,
                " Gets the default currency"
        },
		{
				" Local_ARM_SetDefaultCurrency",
                " RR",							// 1 parametre d'entree + 1 parametre de retour
                " ARM_SetDefaultCurrency",
                " Currency",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_SetDefaultCurrency,
                " Sets the default currency", 
                " Currency"
        },
		{
				" Local_ARM_ADDPERIOD",
                " RRRRRRR",						// 6 parametre d'entree + 1 parametre de retour
                " ARM_ADDPERIOD",
                " date,frequency,[ccy],[nbperiods],[adjRule],[goToEndOfMonth]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ADDPERIOD,
                " Adds a period (or more) to the date", 
                " date",
                " frequency",
				" ccy (dafault : default currency)",
				" number of Periods (default 1)",
				" adjusting Rule (default NONE)",
				" Go to End of Month (default: 0 for false)"
        },
		{
				" Local_ARM_Cover",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_Cover",
                " Security id,Model id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_Cover,
                " Covers a security with the given model", 
                " Security id ", 
                " Model id "
        },
		{
				" Local_ARM_Price_OptUnder",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_Price_OptUnder",
                " Security id,Model id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_Price_OptUnder,
                " Prices a security with the given model and returns an array containing option and underlying prices", 
                " Security id",
                " Model id"
        },
		{
				" Local_ARM_GetPID",
                " R",						// 0 parametre d'entree + 1 parametre de retour
                " ARM_GetPID",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_GetPID,
                " Gets the last date of compilation of this version of ARM"
        },
		{
				" Local_ARM_ClonedAndSetNotional",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_ClonedAndSetNotional",
                " Security id,RefVal Id,[percentRemainder],[interpolationDates]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_ClonedAndSetNotional,
                " Creates a cloned security and sets the notional", 
                " Security id", 
                " Reference value id",
				" percentage of Remainder (default 100)",
				" Dates on which notional will be interpolated (start dates:0, end dates:1, reset dates:2, payment dates:3)"
        },
		{
				" Local_PXL_ARM_ClonedAndSetNotional",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " PXL_ARM_ClonedAndSetNotional",
                " Security id,RefVal Id,[percentRemainder],[interpolationDates]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned security and sets the notional", 
                " Security id", 
                " Reference value id",
				" percentage of Remainder (default 100)",
				" Dates on which notional will be interpolated (start dates:0, end dates:1, reset dates:2, payment dates:3)"
        },
		{		" Local_ARM_ClonedAndSet",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " ARM_ClonedAndSet",
                " Security id, value to set, [type to set]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned spreadoption and changes an attribute", 
                " SpreadOption id",
                " value to set set (refvalue or numeric)",
				" type to set (default STRIKE)"
        },
		{		" Local_PXL_ARM_ClonedAndSet",
                " RRRR",						// 3 parametres d'entree + 1 parametre de retour
                " PXL_ARM_ClonedAndSet",
                " Security id, value to set, [type to set]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned spreadoption and changes an attribute", 
                " SpreadOption id",
                " value to set set (refvalue or numeric)",
				" type to set (default STRIKE)"
        },
		{
				" Local_PXL_ARM_Clone",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " PXL_ARM_Clone",
                " object id",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned object", 
                " source object", 
        },
		{
				" Local_PXL_ARM_Clone",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_Clone",
                " object id",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned object", 
                " source object", 
        },
		{
				" Local_ARM_DoPastReset",
                " RRRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_DoPastReset",
                " Security id,Reset Manager Id, AsOf",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Reset Manager Id (or list for spreadoption)",
				" AsOfDate"
        },
		{
				" Local_PXL_ARM_DoPastReset",
                " RRRR",						// 2 parametres d'entree + 1 parametre de retour
                " PXL_ARM_DoPastReset",
                " Security id,Reset Manager Id, AsOf",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Reset Manager Id (or list for spreadoption)",
				" AsOfDate"
        },
		{
				" Local_FIXRATES",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " FIXRATES",
                " Security id,Fixing Rates",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FIXRATES,
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Fixing Rates"
        },
		{
				" Local_FIXRATES",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_FIXRATES",
                " Security id,Fixing Rates",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_FIXRATES,
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Fixing Rates"
        },
		{
				" Local_PXL_FIXRATES",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " PXL_FIXRATES",
                " Security id,Fixing Rates",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Fixing Rates"
        },
		{
				" Local_PXL_FIXRATES",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " PXL_ARM_FIXRATES",
                " Security id,Fixing Rates",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a cloned security and sets the fixing rates", 
                " Security id", 
                " Fixing Rates"
        },
		{
				" Local_ARM_DisplayScheduleValues",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_DisplayScheduleValues",
                " Instrument id,typeValue,[RcvOrPay],[modelId]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_DisplayScheduleValues,
                " Returns a vector containing a type of values of a swapleg schedule", 
                " object id", 
                " type of values (FWDRATE(FR), RAWFWDRATE(RF), INTERESTDAYS(ID), INTERESTTERMS(IT), FLOWVALUE(FV), NOTIONALVALUE(NV), AMORTVALUE(AV), VOLFWD(VF), FLOWVALUEPV(FPV), VOLCAP(VLC))",
				" Rcv or Pay Leg in case of swap (default R)",
				" ModelId to use for pricing the leg before getting the value (default NULL)"
        },
		{
				" Local_ARM_DisplayScheduleDates",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_DisplayScheduleDates",
                " Instrument id,typeDate,[RcvOrPay],[viewInitExch]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_DisplayScheduleDates,
                " Returns a vector containing a type of Dates of a swapleg schedule", 
                " object id", 
                " type of Dates (STARTDATE(SD), ENDDATE(ED), RESETDATE(RD), PAYDATE(PD), FWDSTARTDATE(FSD), FWDENDDATE(FED))",
				" Rcv or Pay Leg in case of swap (default R)",
				" view the inial Exchange Date? (default Yes)"
        },
        {
				" Local_ARM_DisplayReplicPort",
                " RRRRRR",						// 6 parametres d'entree + 1 parametre de retour
                " ARM_DisplayReplicPort",
                " Instrument id,WeightOrStrike,PayoffOrSensi,[RcvOrPay],[ModelId]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Returns a vector containing replic portfolio informations of a swapleg schedule", 
                " object id", 
                " Weight Or Strike (WEIGHT or STRIKE)",
                " Payoff Or Sensi (PAYOFF or SENSI)",
				" Rcv or Pay Leg in case of swap (default R)"
        },
		{
				" Local_ARM_INTERPOL",
                " RRRRR",						// 4 parametre d'entree + 1 parametre de retour
                " ARM_INTERPOL",
                " vecX,vecY,X,[typeInterpol]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_INTERPOL,
                " Interpolation function in 2 dimensions", 
                " X vector",
                " Y vector",
				" X data",
				" Interpolation type: Linear, Continuous, Spline (default LINEAR)"
        },
		{
				" Local_ARM_TRIANGULARINTERPOL",
                " RRRRRR",						// 5 parametre d'entree + 1 parametre de retour
                " ARM_TRIANGULARINTERPOL",
                " vecX,vecY,matZ,X,Y",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_TRIANGULARINTERPOL,
                " Triangular interpolation function in 3 dimensions", 
                " X vector",
                " Y vector",
				" Z matrix",
				" X data",
				" Y data"
        },
		{
				" Local_KImp",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " KImp",
                " SecId,ModId,Price,[Param]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Implicit K calculation", 
                " Security id", 
                " Model id",
				" Price",
				" Parameter (default: PRICE)"
        },
		{
				" Local_KImp",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_KImp",
                " SecId,ModId,Price,[Param]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_KImp,
                " Implicit K calculation", 
                " Security id", 
                " Model id",
				" Price",
				" Parameter (default: PRICE)"
        },
		{
				" Local_ARM_DisplayZC",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_DisplayZC",
                " ZcId",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_DisplayZC,
                " Display yearterms, discount factors and zero rates from a zc curve", 
                " Zc id"
        },
		{
				" Local_BSSpot",
                " RRRR",							// 3 parametres d'entree + 1 parametre de retour
                " BSSpot",
                " Security id,Model id,As of Date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Computes a spot of a security with the given model", 
                " Security id", 
                " Model id",
				" As of date"
        },
		{
				" Local_BSSpot",
                " RRRR",							// 3 parametres d'entree + 1 parametre de retour
                " ARM_BSSpot",
                " Security id,Model id,As of Date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_BSSpot,
                " Computes a spot of a security with the given model", 
                " Security id", 
                " Model id",
				" As of date"
        },
		{
				" Local_CONSTREFVALUE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " CONSTREFVALUE",
                " Value",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a constant reference value", 
                " Value"
        },
		{
				" Local_CONSTREFVALUE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " ARM_CONSTREFVALUE",
                " Value",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_CONSTREFVALUE,
                " Creates a constant reference value", 
                " Value"
        },
		{
				" Local_PXL_CONSTREFVALUE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " PXL_CONSTREFVALUE",
                " Value",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a constant reference value", 
                " Value"
        },
		{
				" Local_PXL_CONSTREFVALUE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " PXL_ARM_CONSTREFVALUE",
                " Value",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a constant reference value", 
                " Value"
        },
		{
				" Local_REFVALUE",
                " RRRRRRR",								// 6 parametres d'entree + 1 parametre de retour
                " REFVALUE",
                " Dates,Values,[ValueType],[Conversion],[Calculation_Method],[Values2]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value", 
                " Dates array", 
                " Values array",
				" Type of value (default: Y)",
				" Conversion (default: 1)",
                " Calculation Method (default : LIN)",
				" values of the 2nd dimension (default NULL)"
        },
		{
				" Local_REFVALUE",
                " RRRRRRR",								// 6 parametres d'entree + 1 parametre de retour
                " ARM_REFVALUE",
                " Dates,Values,[ValueType],[Conversion],[Calculation_Method],[Values2]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_REFVALUE,
                " Creates a reference value", 
                " Dates array", 
                " Values array",
				" Type of value (default: Y)",
				" Conversion (default: 1)",
                " Calculation Method (default : LIN)",
				" values of the 2nd dimension (default NULL)"
        },
		{
				" Local_PXL_REFVALUE",
                " RRRRRRR",								// 6 parametres d'entree + 1 parametre de retour
                " PXL_REFVALUE",
                " Dates,Values,[ValueType],[Conversion],[Calculation_Method],[Values2]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value", 
                " Dates array", 
                " Values array",
				" Type of value (default: Y)",
				" Conversion (default: 1)",
                " Calculation Method (default : LIN)",
				" values of the 2nd dimension (default NULL)"
        },
		{
				" Local_PXL_REFVALUE",
                " RRRRRRR",								// 6 parametres d'entree + 1 parametre de retour
                " PXL_ARM_REFVALUE",
                " Dates,Values,[ValueType],[Conversion],[Calculation_Method],[Values2]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value", 
                " Dates array", 
                " Values array",
				" Type of value (default: Y)",
				" Conversion (default: 1)",
                " Calculation Method (default : LIN)",
				" values of the 2nd dimension (default NULL)"
        },
		{
				" Local_ARM_DisplayRefValue",
                " RRR",								// 2 parametres d'entree + 1 parametre de retour
                " ARM_DisplayRefValue",
                " refvalue, [isDate]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_REFVALUE,
                " Display the content of a reference value", 
                " reference value object",
				" is the abscisses date or not (default 1)"
        },
		{
				" Local_IATHREELEVREFVAL",
                " RRRRRRRR",								// 7 parametres d'entree + 1 parametre de retour
                " IATHREELEVREFVAL",
                " Level 1,Amort 1,Level 2,Amort 2,Level 3,Amort 3,[Notional]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value with three levels", 
                " Level 1", 
                " Amortisement 1",
				" Level 2",
				" Amortisement 2",
                " Level 3",
                " Amortisement 3",
                " Notional value (default: 100)"
		},
		{
				" Local_IATHREELEVREFVAL",
                " RRRRRRRR",								// 7 parametres d'entree + 1 parametre de retour
                " ARM_IATHREELEVREFVAL",
                " Level 1,Amort 1,Level 2,Amort 2,Level 3,Amort 3,[Notional]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_IATHREELEVREFVAL,
                " Creates a reference value with three levels", 
                " Level 1", 
                " Amortisement 1",
				" Level 2",
				" Amortisement 2",
                " Level 3",
                " Amortisement 3",
                " Notional value (default: 100)"
		},
		{
				" Local_PXL_IATHREELEVREFVAL",
                " RRRRRRRR",								// 7 parametres d'entree + 1 parametre de retour
                " PXL_IATHREELEVREFVAL",
                " Level 1,Amort 1,Level 2,Amort 2,Level 3,Amort 3,[Notional]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value with three levels", 
                " Level 1", 
                " Amortisement 1",
				" Level 2",
				" Amortisement 2",
                " Level 3",
                " Amortisement 3",
                " Notional value (default: 100)"
		},
		{
				" Local_PXL_IATHREELEVREFVAL",
                " RRRRRRRR",								// 7 parametres d'entree + 1 parametre de retour
                " PXL_ARM_IATHREELEVREFVAL",
                " Level 1,Amort 1,Level 2,Amort 2,Level 3,Amort 3,[Notional]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a reference value with three levels", 
                " Level 1", 
                " Amortisement 1",
				" Level 2",
				" Amortisement 2",
                " Level 3",
                " Amortisement 3",
                " Notional value (default: 100)"
		},
		{
				" Local_BERMUDANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " BERMUDANXSTYLE",
                " XDates,[ExpiryDates]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a bermudan exercise style type", 
                " Exercise dates",
				" Expiry dates"
        },
		{
				" Local_BERMUDANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " ARM_BERMUDANXSTYLE",
                " XDates,[ExpiryDates]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_BERMUDANXSTYLE,
                " Creates a bermudan exercise style type", 
                " Exercise dates",
				" Expiry dates"
        },
		{
				" Local_PXL_BERMUDANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " PXL_BERMUDANXSTYLE",
                " XDates,[ExpiryDates]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a bermudan exercise style type", 
                " Exercise dates",
				" Expiry dates"
        },
		{
				" Local_PXL_BERMUDANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " PXL_ARM_BERMUDANXSTYLE",
                " XDates,[ExpiryDates]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a bermudan exercise style type", 
                " Exercise dates",
				" Expiry dates"
        },
		{
				" Local_EUROPEANXSTYLE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " EUROPEANXSTYLE",
                " XDate",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a european exercise style type", 
                " Exercise date"
        },
		{
				" Local_EUROPEANXSTYLE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " ARM_EUROPEANXSTYLE",
                " XDate",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_EUROPEANXSTYLE,
                " Creates a european exercise style type", 
                " Exercise date"
        },
		{
				" Local_PXL_EUROPEANXSTYLE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " PXL_EUROPEANXSTYLE",
                " XDate",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a european exercise style type", 
                " Exercise date"
        },
		{
				" Local_PXL_EUROPEANXSTYLE",
                " RR",									// 1 parametre d'entree + 1 parametre de retour
                " PXL_ARM_EUROPEANXSTYLE",
                " XDate",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates a european exercise style type", 
                " Exercise date"
        },
		{
				" Local_AMERICANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " AMERICANXSTYLE",
                " XStartDate,XEndDate",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an american exercise style type", 
                " Start exercise date", 
                " End exercise date"
        },
		{
				" Local_AMERICANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " ARM_AMERICANXSTYLE",
                " XStartDate,XEndDate",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                IDH_ARM_AMERICANXSTYLE,
                " Creates an american exercise style type", 
                " Start exercise date", 
                " End exercise date"
        },
		{
				" Local_PXL_AMERICANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " PXL_AMERICANXSTYLE",
                " XStartDate,XEndDate",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an american exercise style type", 
                " Start exercise date", 
                " End exercise date"
        },
		{
				" Local_PXL_AMERICANXSTYLE",
                " RRR",									// 2 parametres d'entree + 1 parametre de retour
                " PXL_ARM_AMERICANXSTYLE",
                " XStartDate,XEndDate",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Creates an american exercise style type", 
                " Start exercise date", 
                " End exercise date"
        },
		{
				" Local_GetCorrelInst",
                " RRRRRRRRRRRRRRRRRRRRR",							// 20 parametres d'entree + 1 parametre de retour
                " GetCorrelInst",
                " date1,date2,ccy1,index1,fixing1,tenor1,curve1_ccy1,[curve2_ccy1],[nbmonths_curve1_ccy1],ccy2,index2,fixing2,tenor2,curve1_ccy2,[curve2_ccy2],[nbmonths_curve1_ccy2],type,[lambda],[prec],[ccy]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Calculates instantaneous correlation", 
                " startDate", 
                " endDate",
				" currency1",
				" index1",
				" expiry1",
				" tenor1",
				" curve1 used with ccy1",
				" curve2 used with ccy1, default is curve1",
				" nb of months for ccy1 (after nbmonths, curve2 is used instead of curve1), default is 0",
				" currency2",
				" index2",
				" expiry2",
				" tenor2",
				" curve1 used with ccy2",
				" curve2 used with ccy2, default is curve1",
				" nb of months for ccy2 (after nbmonths, curve2 is used instead of curve1), default is 0",
				" type: NOR or LOGNOR",
				" lambda",
				" precision",
				" calendar currency"
		},
		{
				" Local_GetMoyCorrel",
                " RRRRRRRRRRRRRRRRRRRRR",							// 20 parametres d'entree + 1 parametre de retour
                " GetMoyCorrel",
                " date1,date2,ccy1,index1,fixing1,tenor1,curve1_ccy1,[curve2_ccy1],[nbmonths_curve1_ccy1],ccy2,index2,fixing2,tenor2,curve1_ccy2,[curve2_ccy2],[nbmonths_curve1_ccy2],type,[lambda],[prec],[ccy]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Calculates integrated correlation", 
                " startDate", 
                " endDate",
				" currency1",
				" index1",
				" expiry1",
				" tenor1",
				" curve1_ccy1",
				" curve2 used with ccy1, default is curve1",
				" nb of months for ccy1 (after nbmonths, curve2 is used instead of curve1), default is 0",
				" currency2",
				" index2",
				" expiry2",
				" tenor2",
				" curve1 used with ccy2",
				" curve2 used with ccy2, default is curve1",
				" nb of months for ccy2 (after nbmonths, curve2 is used instead of curve1), default is 0",
				" type: NOR or LOGNOR",
				" lambda",
				" precision",
				" calendar currency"
		},
		{
				" Local_GetCorrelQuanto",
                " RRRRRRRRRRRRRRRRRRR",							// 18 parametres d'entree + 1 parametre de retour
                " GetCorrelQuanto",
                " date1,date2,ccy,index,fixing,tenor,cvname1,[cvname2],[switchinmonth],domccy,domindex,forccy,forindex,type,[lambda],[prec],[calccy],[fwdornot]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Calculates instantaneous quanto correlation", 
                " startDate", 
                " endDate",
				" currency",
				" index",
				" fixing",
				" tenor",
				" curve1 used with ccy",
				" curve2 used with ccy, default is curve1",
				" nb of months for ccy (after nbmonths, curve2 is used instead of curve1), default is 0",
				" domestic currency",
				" domestic index",
				" foreign currency",
				" foreign index",
				" type: NOR or LOGNOR",
				" lambda",
				" precision",
				" calendar currency",
				" Fx Spot or Forward"
		},
		{
				" Local_GetMoyCorrelQuanto",
                " RRRRRRRRRRRRRRRRRRR",							// 18 parametres d'entree + 1 parametre de retour
                " GetMoyCorrelQuanto",
                " date1,date2,ccy,index,fixing,tenor,cvname1,[cvname2],[switchinmonth],domccy,domindex,forccy,forindex,type,[lambda],[prec],[calccy],[fwdornot]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Calculates integrated quanto correlation", 
                " startDate", 
                " endDate",
				" currency",
				" index",
				" fixing",
				" tenor",
				" curve1 used with ccy",
				" curve2 used with ccy, default is curve1",
				" nb of months for ccy (after nbmonths, curve2 is used instead of curve1), default is 0",
				" domestic currency",
				" domestic index",
				" foreign currency",
				" foreign index",
				" type: NOR or LOGNOR",
				" lambda",
				" precision",
				" calendar currency",
				" Fx Spot or Forward"
		},
		{
				" Local_ARM_Hedge",
                " RRR",							// 2 parametres d'entree + 1 parametre de retour
                " ARM_Hedge",
                " security,type",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Computes hedge ratio on PRCS", 
                " PRCS Id", 
                " Hedge's type (IRDELTA,IRVEGA...)"
		},
		{
				" Local_ARM_GetMeanRevFromSummit",
                " RRRRRR",							// 5 parametres d'entree + 1 parametre de retour
                " ARM_GetMeanRevFromSummit",
                " ccy,index,cvname,date,[2or3Factor]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Mean Reversion from Summit", 
                " Currency", 
                " Index",
				" Summit Cv Name",
				" date",
				" 2F or 3F (default 3F)"
		},
		{
				" Local_ARM_GetCutOffFromSummit",
                " RRRRRR",							// 5 parametres d'entree + 1 parametre de retour
                " ARM_GetCutOffFromSummit",
                " ccy,index,cvname,numfactor,date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get CutOff from Summit", 
                " Currency", 
                " Index",
				" Summit Cv Name",
				" Number of Factor 2F or 3F",
				" date"
		},
		{
				" Local_ARM_GETINFOFROMPRCS",
                " RRR",						// 2 parametre d'entree + 1 parametre de retour
                " ARM_GETINFOFROMPRCS",
                " PRCS,DataType",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets information from a PRCS",
                " PRCS Id",
				" Wanted data:(FUNDINGCCY,FXNUMCCY,FXUNDCCY,FXNUMNOT, FUNDPV or FUNDINGPV)"
        },
		{
				" Local_ARM_GETOBJINFOFROMPRCS",
                " RRR",						// 2 parametre d'entree + 1 parametre de retour
                " ARM_GETOBJINFOFROMPRCS",
                " PRCS,DataType",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets object information from a PRCS",
                " PRCS Id",
				" string corresponding to the data you want to get (FXNUMCPN,FXUNDCPN,FX0,CAP or FLOOR)"
        },
		{
				" Local_ARM_GETOPTIONDATES",
                " RRR",						// 2 parametre d'entree + 1 parametre de retour
                " ARM_GETOPTIONDATES",
                " PRCS,DateType",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets option dates (Notice or Cancel) from a PRCS",
                " PRCS Id",
				" string corresponding to the data you want to get (N for notice, C for cancel)"
        },
		{
				" Local_ARM_GETDUALOPTIONSTRIKE",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_GETDUALOPTIONSTRIKE",
                " PRCS",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets dual option strike from a PRCS",
                " PRCS Id"
        },
		{
				" Local_ARM_CptRefValue",
                " RRR",						// 2 parametre d'entree + 1 parametre de retour
                " ARM_CptRefValue",
                " RefValue,date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " computes in a refvalue",
                " Reference Value Id",
				" date"
        },
		{
				" Local_ARM_GetFixing",
                " RRRRRR",							// 5 parametres d'entree + 1 parametre de retour
                " ARM_GetFixing",
                " source,index,term,ccy,date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Fixing from Summit", 
                " Summit source of fixing (T248,...)", 
                " Index",
				" Term (3M,...)",
				" Currency",
				" Date"
		},
		{
				" Local_ARM_GetFixingFromCalypso",
                " RRRRRRR",							// 6 parametres d'entree + 1 parametre de retour
                " ARM_GetFixingFromCalypso",
                " index,term,ccy,source,cvname,date",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Fixing from Calypso", 
                " Index",
				" Term (3M,...)",
				" Currency",
				" Source of fixing (T248,...)", 
				" Pricing Env (default MO)", 
				" Date"
		},
		{

				" Local_ARM_GetDealsFromSummitFilter",
                " RR",						// 1 parametre d'entree + 1 parametre de retour
                " ARM_GetDealsFromSummitFilter",
                " filter",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " retrieves list of deals from a summit filter",
                " Summit filter name"
        },
        {
				" Local_ARM_GetAsOfVolOrRate",
                " RRRRRRRRRR",						// 9 parametres d'entree + 1 parametre de retour
                " GetAsOfVolOrRate",
                " AsOfDate,Currency,Index,CvName,Expiry,Maturity,[YieldOrVol],[CalcMod],[VolType]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Retrieve a historic rate or Vol. from Summit",
                " AsOfDate",
                " Currency (EUR, USD, ...)",
				" Index (EURIB, PIBOR, ...)", 
				" CvName (MO, ...)",
				" Expiry (2M, 3Y, ...)",
                " Maturity (3M, 6Y, ...)",
				" Yield or volatility (Default : Y)",
				" Calculation mode (Default : LOGNOR)",
				" Volatility type (Default : IRG)"
        },
        {
				" Local_ARM_GetAsOfVolOrRate",
                " RRRRRRRRRR",						// 9 parametres d'entree + 1 parametre de retour
                " ARM_GetAsOfVolOrRate",
                " AsOfDate,Currency,Index,CvName,Expiry,Maturity,[YieldOrVol],[CalcMod],[VolType]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Retrieve a historic rate or Vol. from Summit",
                " AsOfDate",
                " Currency (EUR, USD, ...)",
				" Index (EURIB, PIBOR, ...)", 
				" CvName (MO, ...)",
				" Expiry (2M, 3Y, ...)",
                " Maturity (3M, 6Y, ...)",
				" Yield or volatility (Default : Y)",
				" Calculation mode (Default : LOGNOR)",
				" Volatility type (Default : IRG)"
        },
		{
				" Local_ARM_SumRefValue",
                " RRRR",						// 3 parametre d'entree + 1 parametre de retour
                " ARM_SumRefValue",
                " RefValue1,[RefValue2],[coef]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " create a refvalue = R1 + coef * R2 or coef * R1 if R2 is NULL",
                " Reference Value Id 1",
				" Reference Value Id 2 (default NULL)",
				" coef (default 1)"
        },
		{
				" Local_PXL_ARM_SumRefValue",
                " RRRR",						// 3 parametre d'entree + 1 parametre de retour
                " PXL_ARM_SumRefValue",
                " RefValue1,[RefValue2],[coef]",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " create a refvalue = R1 + coef * R2 or coef * R1 if R2 is NULL",
                " Reference Value Id 1",
				" Reference Value Id 2 (default NULL)",
				" coef (default 1)"
        },
		{
				" Local_ARM_GetLastDateWarm",
                " R",						// 0 parametre d'entree + 1 parametre de retour
                " ARM_GetLastDateWarm",
                " ",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Gets the last hour of the summit database",
        },
        {
				" Local_ARM_GetFwdRatesMatrix",
                " RRRRR",						// 4 parametres d'entree + 1 parametre de retour
                " ARM_GetFwdRatesMatrix",
                " AsOfDate,Currency,Index,CvName",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Retrieve Swap Fwd Rate Matrix",
                " AsOfDate",
                " Currency (EUR, USD, ...)",
				" Index (EURIB, PIBOR, ...)", 
				" CvName (MO, ...)",
        },
        {
				" Local_ARM_GetInfo",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_GetInfo",
                " secId,type",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Info from a security",
                " Instrument Id",
                " type"
        },
        {
				" Local_ARM_GetInfo",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_GetInfo",
                " secId,type",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Info from a security",
                " Instrument Id",
                " type"
        },
        {
				" Local_PXL_ARM_GetInfo",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " PXL_ARM_GetInfo",
                " secId,type",
                " 0",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Info from a security",
                " Instrument Id",
                " type"
        },
        {
				" Local_ARM_SecurityFlows",
                " RRR",						// 2 parametres d'entree + 1 parametre de retour
                " ARM_SecurityFlows",
                " labels,values",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Security flows",
                " labels : vector of string",
                " values : vector of double"
        },
        {
				" Local_ARM_GetModelFactorFromSummit",
                " RRRRRRRRR",				// 8 parametres d'entree + 1 parametre de retour
                " ARM_GetModelFactorFromSummit",
                " AsOfDate,Model,InstType,FactorName,Currency,Index,CvName,[Calculation_Method]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " Get Model Factor",
                " as of date",
                " model : string",
				" inst type : IRG, SWOPT, ...",
				" factor name : REVERSION, FACTOR1, ...",
                " currency : EUR, USD, ...",
				" index : EURIB, LIBOR, ...",
				" cvName : MO, ...",
				" calculation method (default LINEAR)"
        },
		{
				" Local_ARM_MatrixVectorViewer",
                " RR",				// 1 parametre d'entree + 1 parametre de retour
                " ARM_MatrixVectorViewer",
                " Object",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " View the Matrix or Vector ID in Excel cells",
                " 1"
        },
		{
				" Local_ARM_EtkConnect",
                " RRRRRR",						// 5 parametre d'entree + 1 parametre de retour
                " ARM_AnyEtkConnect",
                " username,passwd,context,itconfigdomainsdir,itdomainname",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " for connecting to any Summit Database",
                " username",
                " passwd",
				" context",
				" itconfigdomainsdir",
                " itdomainname"
        },
		{
        		" Local_ARM_GP_CQSO_Create",			/// name of the C++ function
				" RRRRRRRRRRRRRRRRR",					/// 17 parametres = 16 d'entree + 1 parametre de retour 
				" ARM_GP_CQSO_Create",
				" startDate,endDate,FundingCurrency,UnderlyingCurrency,1stIndex,2ndIndex,ResetDates,Notionals,Margins,Strikes,Leverages1,Leverages2,CpnMin,CpnMax,Fees,ScheduleArguments",
				" 1",									/// visible in excel
				XLLOCALARM_UTIL_GROUP,
				" ",
				" ",
				" Create a generic security corresponding to a callable quanto spread option",
				" startDate",
				" endDate",
				" FundingCurrency",
				" UnderlyingCurrency",
				" 1stIndex",
				" 2ndIndex",
				" ResetDates",
				" Notionals",
				" Margins",
				" Strikes",
				" Leverages1",
				" Leverages2",
				" CpnMin",
				" CpnMax",
				" Fees",
				" ScheduleArguments"
		},
		{
	
				" Local_ARM_NumericalRegression",
                " RRRRR",						// 4 parametre d'entree + 1 parametre de retour
                " NumericalRegression",
                " refValue,newValue,[tolerance],[epsilon]",
                " 1",
                XLLOCALARM_UTIL_GROUP,
                " ",
                " ",
                " check if newValue is equal to refValue according to tolerance ",
                " reference value",
                " new value",
				" tolerance defaulted to 10E-8",
				" technical value for testing around zero, defaulted to 10E-14 "
        },

        /// END UTILITIES
