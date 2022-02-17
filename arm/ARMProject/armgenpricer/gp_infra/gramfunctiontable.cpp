/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctiontable.cpp
 *
 *  \brief this is more a configuration file
 *  we made lots of efforts to make this as simple as possible..
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/gramfunctiontable.h"
#include "gpinfra/gramfunctordef.h"
#include "gpinfra/gramnodebuilder.h"
#include "gpbase/gpmatrix.h"

CC_BEGIN_NAMESPACE( ARM )

/// some constants
const char* DEFAULT_CCY_CHAR	= "DEFAULT_PER_CURRENCY";
double		DEFAULT_CCY_DOUBLE	= GETDEFAULTVALUE;

///	list of default argument ... in order to make it really built
/// we need to declare it as static variables!
///															itsChar					itsDouble	
const ARM_GramFunctionArgDef DefaultPerCcyCharStruct	= { "DEFAULT_PER_CURRENCY"	}; 
const ARM_GramFunctionArgDef DefaultPerCcyDbleStruct	= {	"",						GETDEFAULTVALUE };
const ARM_GramFunctionArgDef DefaultDbleZeroStruct		= {	"",						0 };
const ARM_GramFunctionArgDef DefaultDbleOnePerCentStruct= { "",						0.01 };
const ARM_GramFunctionArgDef DefaultDbleOneStruct		= { "",						1 };
const ARM_GramFunctionArgDef DefaultDbleTwoStruct		= { "",						2.0 };
const ARM_GramFunctionArgDef DefaultDbleInfiniteStruct	= { "",						1e20 };
const ARM_GramFunctionArgDef DefaultDbleTwentyStruct	= { "",					20 };
const ARM_GramFunctionArgDef DefaultDbleFiftyStruct		= { "",					50 };
const ARM_GramFunctionArgDef DefaultDbleMagic54Struct	= { "",					0.54 };
const ARM_GramFunctionArgDef DefaultDbleMinusOneStruct	= { "",						-1 };
const ARM_GramFunctionArgDef DefaultResetTimingStruct	= { "ADV"		}; 
const ARM_GramFunctionArgDef DefaultPaymentTimingStruct	= { "ARR"		}; 
const ARM_GramFunctionArgDef Default30Per360			= { "30/360"	}; 
const ARM_GramFunctionArgDef DefaultActualActual		= { "ACTUAL"	}; 
const ARM_GramFunctionArgDef DefaultOneYear				= { "1Y"		}; 
const ARM_GramFunctionArgDef DefaultThreeMonths			= { "3M"		}; 
const ARM_GramFunctionArgDef DefaultZeroMonths			= { "0M"		}; 
const ARM_GramFunctionArgDef DefaultCountryChar			= { ARM_DEFAULT_COUNTRY }; 
const ARM_GramFunctionArgDef DefaultNothing			    = { ""			}; 
const ARM_GramFunctionArgDef DefaultLN					= { "LN"		}; 
const ARM_GramFunctionArgDef DefaultIntRule				= { "ADJ"		};
const ARM_GramFunctionArgDef DefaultIndexType			= { "LIBOR"		};
const ARM_GramFunctionArgDef DefaultIndexTypeCms		= { "CMS"		};
const ARM_GramFunctionArgDef DefaultFixedIndex			= { "FIXED"		};
const ARM_GramFunctionArgDef DefaultStubType			= { "SS"		};
const ARM_GramFunctionArgDef DefaultEXNotionalType		= { "BOTH"		};
const ARM_GramFunctionArgDef DefaultCFType				= { "CAP"		};
const ARM_GramFunctionArgDef DefaultWeeklyFreqType		= { "W"			};
const ARM_GramFunctionArgDef DefaultDigitTypeStruct		= { "ANALYTIC"	}; 
const ARM_GramFunctionArgDef DefaultEpsilonStruct		= { "",					1e-4 }; 

/// Argument lists for the various functions.
static ARM_GramFunctionArgStruct IfArg[] = 
{ 
	/// Arg type			Required?	VaArg?		Name				Description                             Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"Condition    ",	"Condition to evaluate for each path",  true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"Value If True",	"Value returned if Condition is true",  true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"Value If False",	"Value returned if Condition is false", true                    },
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct StatAverageArg[] = 
{ 
	/// Arg type			Required?	VaArg?		Name				Description                             Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"New Value     ",	"New Value Per States",                 true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"Initial Values",	"Initial Values to compute average",    true                    },
	{ GFAT_DOUBLE_TYPE,		true,		false,		"Average size  ",	"Number of Elements to use in Average", true                    },
	{ GFAT_DOUBLE_TYPE,		false,		false,		"Reset         ",	"Reset initial values",					true ,	&DefaultDbleZeroStruct},

	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct StatStdDevArg[] = 
{ 
	/// Arg type			Required?	VaArg?		Name				Description                             Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"New Value     ",	"New Value Per States",                 true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"Average Value ",	"Average Value Per States",             true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"Initial Values",	"Initial Values to compute average",    true                    },
	{ GFAT_DOUBLE_TYPE,		true,		false,		"Average size  ",	"Number of Elements to use in Average", true                    },
	{ GFAT_DOUBLE_TYPE,		false,		false,		"Reset         ",	"Reset initial values",					true ,	&DefaultDbleZeroStruct},
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct MinMaxArg[] = 
{ 
	/// Arg type			Required?	VaArg?		Name				Description	        Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"A",				"First argument",   true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"B",				"Second argument",  true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};



static ARM_GramFunctionArgStruct DoubleArg[] =
{
	/// Arg type			Required?	VaArg?		Name				Description	    Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"x",				"Operand",      true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct PowArg[] =
{
	/// Arg type			Required?	VaArg?		Name				Description     Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"x",				"Base",         true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"y",				"Exponent",     true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SumSerieArg[] = 
{
	/// Arg type			Required?	VaArg?		Name				Description     Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"x",				"Base",         true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"y",				"Exponent",     true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct BinArg[] =
{
	/// Arg type			Required?	VaArg?		Name				Description         Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"lhs",				"left hand side",   true                    },
	{ GFAT_VECTOR_TYPE,		true,		false,		"rhs",				"right hand side",  true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct DCFArg[] =
{
	/// Arg type			        Required?   VaArg?		Name				Description												                Visible?    Default Value
	{ GFAT_DATE_OR_VECTOR_TYPE,		true,		false,		"From     ",		"From this date",                                                       true                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		true,		false,		"To       ",		"To this date",                                                         true                    },
	{ GFAT_MULTITOKENSTRING_TYPE,   true,	    false,		"Day Count",		"Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",  true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct DFArg[] = 
{ 
	/// Arg type			Required?	VaArg?		Name				 Description                                                                                    Visible?    Default Value
	{ GFAT_MODEL_TYPE,			true,		false,		"Model            ",  "The model to compute the discount factor with",         true                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,	true,	 	false,		"Vector or End Date", "Vector of Julian or End Date of the dicount Factor",    true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct LiborArg[] = 
{ 
	/// Arg type			        Required?	VaArg?		Name					Description												                                            Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model               ",	"The model to compute the Libor on",                                                                true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date          ",	"Start Date of the rate",                                                                           true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date",	"Maturity or End Date of the cap/floor : a maturity is nb followed by d/w/m/y like 6m or 1y",       true                                    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Days count          ",	"Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",					true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Reset Gap or Date   ",	"Reset date or number of days between fixing and start date",										true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Pay Gap  or Date    ",	"Payment date or number of days between start date(ADV)/end date(ARR) and payment date",			true,       &DefaultDbleZeroStruct      },
	{ GFAT_STRING_TYPE,		        false,		false,		"Payment Timing      ",	"Payment in-advance (ADV) or in-arrears (ARR), default is ARR",										true,       &DefaultPaymentTimingStruct },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct AnnuityArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name				 Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model            ", "The model to compute with",                                                                       true                                    },
	{ GFAT_DATE_TYPE,		        true,       false,		"Start Date       ", "Start date of the annuity",                                                                       true                                    }, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date", "Tenor or End Date of Annuity, a tenor is nb followed by d/w/m/y like 6m or 1y",                   true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period ", "Fix Leg Period = nb followed by d/w/m/y like 6m or 1y",									        true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Days count       ", "Fix Leg Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Notional         ", "Fix Leg Notional",							                                                    true,       &DefaultDbleOneStruct		},
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Date Strip       ", "Date strip describing the schedule of reset, settlement & payment dates",				            true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Date Strip Offset", "Schedule offset in the date strip, default is 0",										            true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Date Strip Width ",  "Schedule offset in the date strip, default is -1",												true,		&DefaultDbleMinusOneStruct  },
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct YTPArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name				 Description												                                        Visible?    Default Value
	{ GFAT_DATE_TYPE,		        true,		false,		"Settlement Date  ", "SettlementDate",                                                                                  true                                },
	{ GFAT_DATE_TYPE,		        true,		false,		"Coupon Date      ", "First coupon date",                                                                               true                                }, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date", "Tenor or End Date, a tenor is nb followed by d/w/m/y like 6m or 1y",                              true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Coupon           ", "Coupon in real terms",                                                                            true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Notional         ", "Redemption Value in real terms",                                                                  true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Yield            ", "Yield Value in real terms",                                                                       true                                },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period ", "Coupon Frequency Period = nb followed by d/w/m/y like 6m or 1y",									true,       &DefaultOneYear         },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,		false,		"Days count       ", "Coupon Days Count Convention: ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultActualActual    },
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct PTYArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name				 Description												                                        Visible?    Default Value
	{ GFAT_DATE_TYPE,		        true,		false,		"Settlement Date  ", "SettlementDate",                                                                                  true                                }, 
	{ GFAT_DATE_TYPE,		        true,		false,		"Coupon Date      ", "First coupon date",                                                                               true                                }, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date", "Tenor or End Date, a tenor is nb followed by d/w/m/y like 6m or 1y",                              true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Coupon           ", "Coupon in real terms",                                                                            true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Notional         ", "Redemption Value in real terms",                                                                  true                                },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Price            ", "Price in real terms",                                                                             true                                },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period ", "Coupon Frequency Period = nb followed by d/w/m/y like 6m or 1y",									true,       &DefaultOneYear         },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,		false,		"Days count       ", "Coupon Days Count Convention: ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultActualActual    },
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct SwapRateArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name						 Description												                                Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model						", "The model to compute with",															true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date					", "Start Date of the rate",															true									}, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date			", "Tenor or End Date of Swap, a tenor is nb followed by d/w/m/y like 6m or 1y",		true									},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixed Frequency Period		", "Fixed Leg Period = nb followed by d/w/m/y like 6m or 1y",							true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Fixed Days count			", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",				true,       &DefaultPerCcyCharStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Floating Frequency Period	", "Floating Leg Period = nb followed by d/w/m/y like 6m or 1y",						true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Floating Days count		", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",				true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Margin						", "Floating leg margin in real terms",													true,       &DefaultDbleZeroStruct		},
	{ GFAT_STRING_TYPE,		        false,		false,		"StubType					",	"Stub Type SS/LS/LS/LE",															true,       &DefaultStubType			},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Dble notional flag			", "****** Internal argument : default flag is NO double notional use ******",			false,      &DefaultDbleZeroStruct		},
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"floating Date Strip		",  "Date strip describing the schedule of reset, settlement & payment dates",			false,      &DefaultDbleZeroStruct		},
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"fixed Date Strip			",  "Date strip describing the schedule of reset, settlement & payment dates",			false,      &DefaultDbleZeroStruct		},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Float Date Strip Offset	",  "Schedule offset in the date strip, default is 0",									false,      &DefaultDbleZeroStruct		},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Fix Date Strip Offset		",	"Schedule offset in the date strip, default is 0",									false,      &DefaultDbleZeroStruct		},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Float Date Strip Width		",  "Schedule offset in the date strip, default is -1",									false,      &DefaultDbleMinusOneStruct  },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Fix Date Strip Width		",  "Schedule offset in the date strip, default is -1",									false,      &DefaultDbleMinusOneStruct  },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SpreadArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name						 Description												                                Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model						", "The model to compute with",                                                             true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date					", "Start Date of the rate",                                                                true									}, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date Rate 1	", "Tenor or End Date of Swap 1, (tenor is nb followed by d/w/m/y like 6m or 1y)",          true									},
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date Rate 2	", "Tenor or End Date of Swap 1, (tenor is nb followed by d/w/m/y like 6m or 1y)",          true									},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Coeff1						", "Coeff for CMS1 (default is 1)",															true,	&DefaultDbleOneStruct			},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Coeff2						", "Coeff for CMS2 (default is 1)",															true,	&DefaultDbleOneStruct			},

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct MaxRateArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name						 Description												                                Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                      ", "The model to compute with",                                                             true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date                 ", "Start Date of the rate",                                                                true									},
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date Rate     ", "Tenor or End Date of Swap, (tenor is nb followed by d/w/m/y like 6m or 1y)",			true									},
	{ GFAT_DATE_TYPE,		        true,		false,		"first Reset Date           ", "MinMax computed between first reset date and present reset date",						true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"first Start Date           ", "Start Date of the rate that begins minmax calculation",		                            true									},
	{ GFAT_VECTOR_TYPE,				true,		false,		"first Rate                 ", "PreComputed first rate that begins minmax calculation",									true									},
	{ GFAT_STRING_TYPE,				true,		false,		"Min/Max/MaxMin             ", "Min Or Max Or Max-Min Option Calculation",												true									},
	{ GFAT_STRING_TYPE,				false,		false,		"Reset Frequency            ", "Reset Frequency in observation period (D/W/M/B/Q/S/A)",									true,	&DefaultWeeklyFreqType			},
	{ GFAT_VECTOR_TYPE,				false,		false,		"Strike                     ", "Strike in option case",																	true,	&DefaultDbleZeroStruct			},
	{ GFAT_STRING_TYPE,				false,		false,		"Cap/Floor                  ", "Cap or Floor in option case",															true,	&DefaultCFType					},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Min vs Max correl          ", "Min vs Max correlation in option case",													true,	&DefaultDbleMagic54Struct		},
	{ GFAT_VECTOR_TYPE,		        false,		false,		"Min/Max Accrued            ", "[0]=Min and [1]=Max accrued rate",														true,	&DefaultDbleZeroStruct			},
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SwapArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name						 Description												                                Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model						", "The model to compute with",                                                             true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date					", "Start Date of the rate",                                                                true									}, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date			", "Tenor or End Date of Swap, a tenor is nb followed by d/w/m/y like 6m or 1y",            true									},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	true,		false,		"Fixed Rate					", "Swap's fixed rate in real terms",                                                       true									},
	{ GFAT_STRING_TYPE,		        true,		false,		"Pay/rec					", "Swap's type PAY or RCV",                                                                true									},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixed Frequency Period		", "Fixed Leg Period = nb followed by d/w/m/y like 6m or 1y",								true,		&DefaultPerCcyCharStruct	},
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Fixed Days count			", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",					true,		&DefaultPerCcyCharStruct	},
	{ GFAT_STRING_TYPE,		        false,		false,		"Floating Frequency Period	", "Floating Leg Period = nb followed by d/w/m/y like 6m or 1y",							true,		&DefaultPerCcyCharStruct	},
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Floating Days count		", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",					true,		&DefaultPerCcyCharStruct	},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Margin						", "Floating leg margin in real terms",							                            true,       &DefaultDbleZeroStruct      },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Swap's notional			", "Swap's notional in real terms, default is 1",											true,		&DefaultDbleOneStruct		},
	{ GFAT_STRING_TYPE,		        false,		false,		"StubType					", "Stub Type SS/LS/LS/LE",																    true,       &DefaultStubType			},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Dble notional flag			", "****** Internal argument : default flag is NO double notional use ******",              false,		&DefaultDbleZeroStruct		},
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"floating Date Strip		",  "Date strip describing the schedule of reset, settlement & payment dates",				false,      &DefaultDbleZeroStruct      },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"fixed Date Strip			",  "Date strip describing the schedule of reset, settlement & payment dates",				false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Float Date Strip Offset	",  "Schedule offset in the date strip, default is 0",										false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Fix Date Strip Offset		",  "Schedule offset in the date strip, default is 0",										false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Float Date Strip Width		",  "Schedule offset in the date strip, default is -1",										false,      &DefaultDbleMinusOneStruct  },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Fix Date Strip With 		",  "Schedule offset in the date strip, default is -1",										false,      &DefaultDbleMinusOneStruct  },


	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct BasisSwapArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name											 Description																								Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"DomModel										", "The domestic model to compute floating leg with",														true									},
	{ GFAT_MODEL_TYPE,		        true,		false,		"ForModel										", "The foreign model to compute fixed leg with",															true									},
	{ GFAT_MODEL_TYPE,		        true,		false,		"FxModel										", "The Forex model to compute Notional Exchange with",														true									},
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date										", "Start Date of the rate",																				true									}, 
	{ GFAT_DATEORMATU_TYPE,         true,		false,		"Tenor or End Date								", "Tenor or End Date of Swap, a tenor is nb followed by d/w/m/y like 6m or 1y",							true									},
	{ GFAT_STRING_TYPE,		        true,		false,		"Pay/rec										", "Swap's type PAY or RCV domestic flows",																	true									},
	{ GFAT_STRING_TYPE,		        false,		false,		"Domestic frequency Period						", "Domestic Leg Period = nb followed by d/w/m/y like 6m or 1y",											true,       &DefaultPerCcyCharStruct	},
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Domestic days count							", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",									true,       &DefaultPerCcyCharStruct	},																																																	
	{ GFAT_STRING_TYPE,		        false,		false,		"Foreign frequency Period						", "Foreign Leg Period = nb followed by d/w/m/y like 6m or 1y",												true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Foreign days count								", "Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",									true,       &DefaultPerCcyCharStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Domestic Type									", "Type of domestic rate to flottant or fix rate, FIXED,FLOTTANT",											true,       &DefaultFixedIndex			},
	{ GFAT_STRING_TYPE,		        false,		false,		"Foreign Type									", "Type of foreign rate to flottant or fix rate, FIXED,FLOTTANT",											true,       &DefaultFixedIndex			},
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Fx Reset Gap  or Date							", "Forex reset date or number of days between start date(ADV)/end date(ARR) and forex fixing date",		true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Fx Settlment Gap  or Date						", "Forex settlement date or number of days between forex expiry date and forex settlment date",			true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Domestic margin								", "Domestic leg margin in real terms",																		true,       &DefaultDbleZeroStruct		},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Domestic Fixed Rate							", "Swap's domestic fixed rate in real terms",																true,		&DefaultDbleZeroStruct		},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Foreign margin									", "Foreign leg margin in real terms",																		true,       &DefaultDbleZeroStruct      },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Foreign Fixed Rate								", "Swap's foreign fixed rate in real terms",																true,		&DefaultDbleZeroStruct	    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Swap's domestic notional						", "Swap's notional for domestic leg in real terms, default is 1 if it is pay Ccy",							true,       &DefaultDbleOneStruct       },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Swap's foreign notional						", "Swap's notional for foreign leg in real terms, default is 1 if it is pay Ccy",							true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,				false,		false,		"ExchangeNotionalType							", "Swap's notional Type, default is BOTH ( START, END)",													true,       &DefaultEXNotionalType      },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Domestic Date Strip							", "Date strip describing the schedule of reset, settlement & payment dates",								false,      &DefaultDbleZeroStruct      },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Foreign Date Strip								", "Date strip describing the schedule of reset, settlement & payment dates",								false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Domestic Date Strip Offset						", "Schedule offset in the date strip, default is 0",														false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Domestic Date Strip width						", "Schedule offset in the date strip, default is 0",														false,      &DefaultDbleMinusOneStruct  },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Foreign Date Strip Offset						", "Schedule offset in the date strip, default is 0",														false,      &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Foreign Date Strip width						", "Schedule offset in the date strip, default is 0",														false,      &DefaultDbleMinusOneStruct  },
																											
																											
	/// says that this is the end
	{ GFAT_TERMINATOR }
};




static ARM_GramFunctionArgStruct PVArg[] =
{
	/// Arg type			Required?		VaArg?		Name				Description									Visible?    Default Value
	{ GFAT_FUTUREREF_TYPE,	true,			false,		"Future Payoff",	"Reference to a Deal Description Cell",     true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct UnPayArg[] =
{
	/// Arg type			Required?	VaArg?		Name				Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,		true,		false,		"Payoff to unpay",	"Cancels the discounting on the payoff",    true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct ExerciseArg[] =
{
	/// Arg type			            Required?	VaArg?		Name								Description														Visible?	Default Value		AddRefNode
	{ GFAT_VECTOR_TYPE,			        true,		false,		"Payoff1",							"Always Paid Payoff",											true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_TYPE,			        true,		false,		"Payoff2",							"Exercise Payoff",												true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_OR_FUTUREREF_TYPE,	true,		false,		"Continuation",					    "Continuation Payoff",							                true,		&DefaultNothing,			},
	{ GFAT_VECTOR_TYPE,			        false,		true,		"Regression Variables",				"Reference to Regression Variables",							true,		&DefaultNothing,	true	},

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct FrontierArg[] =
{
	/// Arg type			            Required?	VaArg?		Name								Description														Visible?	Default Value		AddRefNode
	{ GFAT_VECTOR_TYPE,			        true,		false,		"Intrinsic",						"Intrinsic Value",												true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_OR_FUTUREREF_TYPE,	true,		false,		"Continuation",					    "Continuation Payoff",							                true,		&DefaultNothing,			},
	{ GFAT_VECTOR_TYPE,			        true,		false,		"Index",							"Index to be evaluated on frontier",							true,		&DefaultNothing,	true	},
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct TriggerArg[] =
{
	/// Arg type			Required?		VaArg?		Name								Description																	Visible?    Default Value		AddRefNode
	{ GFAT_VECTOR_TYPE,			true,		false,		"VarTest1",							"First test variable",														true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_TYPE,			true,		false,		"VarTest2",							"Second test variable",														true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_TYPE,			true,		false,		"Payoff1",							"Reference to what is always paid",											true,		&DefaultNothing,	true	},
	{ GFAT_VECTOR_TYPE,			true,		false,		"Payoff2",							"Reference to what is paid if trigger occurs (VarTest1>=VarTest2)",			true,		&DefaultNothing,	true	},
	{ GFAT_FUTUREREF_TYPE,		true,		false,		"Future Payoff",					"Reference to what happens if trigger does not occur",						true									},

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CapletArg[] =
{
	/// Arg type			        Required?	VaArg?		Name					Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model               ", "The model to compute the capfloor with",                                                       true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date          ", "Start Date of the cap/floor",                                                                  true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date", "Maturity or End Date of the cap/floor : a maturity is nb followed by d/w/m/y like 6m or 1y",   true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike              ", "Option's strike in real terms",                                                                true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Cap/Floor           ", "Option's type CAP or FlOOR",                                                                   true                                    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Cpn Days count      ", "Coupon Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Reset Gap or Date   ",	"Reset date or number of days between fixing and start date",									true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Pay Gap  or Date    ",	"Payment date or number of days between start date(ADV)/end date(ARR) and payment date",		true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Option's notional   ", "Option's notional in real terms, default is 1",												true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing        ", "Reset in-advance (ADV) or in-arrears (ARR), default is ADV",							        true,       &DefaultResetTimingStruct   },
	{ GFAT_MATURITY_TYPE,           false,	    false,		"Index Term          ", "Index term : a nb followed by d/w/m/y like 6m or 1y",                            		        true,       &DefaultNothing             },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Index Days count    ", "Index Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		    true,       &DefaultPerCcyCharStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct SwaptionArg[] =
{
	/// Arg type			        Required?	VaArg?		Name			 	 Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model            ", "The model to compute with",                                                                       true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date       ", "Date which the option is paid (= start of the swap rate)",                                        true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Tenor or End Date", "Tenor or End Date of the Swap, a tenor is nb followed by d/w/m/y like 6m or 1y",                  true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike           ", "Swaption's strike in real terms",                                                                 true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Pay/rec          ", "Swaption's type PAY or RCV",                                                                      true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period ", "Fix Leg Period = nb followed by d/w/m/y like 6m or 1y",											true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Days count       ", "Fix Leg Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Reset Gap or Date", "Reset date or number of days between fixing and start date",									    true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Option's notional", "Swaption's notional in real terms, default is 1",													true,       &DefaultDbleOneStruct       },


	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct CapArg[] =
{
	/// Arg type			        Required?	VaArg?		Name			    	   Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                  ", "The model to compute with",                                                                     true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date             ", "Date which the option is paid (= start of the swap rate)",                                      true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date   ", "Tenor or End Date of the Cap, a tenor is nb followed by d/w/m/y like 6m or 1y",                 true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                 ", "Cap's strike in real terms",                                                                    true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Cap/Floor              ", "Option's type CAP or FlOOR",                                                                    true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period       ", "Cpn Frequency Period = nb followed by d/w/m/y like 6m or 1y",									true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Cpn Days count         ", "Coupon Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Reset Gap              ", "Number of days between fixing and start date",											        true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Pay Gap                ", "Number of days between end and payment date",										            true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Option's notional      ", "Option's notional in real terms, default is 1",													true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing           ", "Reset in-advance (ADV) or in-arrears (ARR), default is ADV",							        true,       &DefaultResetTimingStruct   },
	{ GFAT_MATURITY_TYPE,           false,	    false,		"Index Term          ", "Index term : a nb followed by d/w/m/y like 6m or 1y",                            		            true,       &DefaultNothing             },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Index Days count       ", "Index Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		    true,       &DefaultPerCcyCharStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SumOptionArg[] =
{
	/// Arg type			        Required?	VaArg?		Name			    	   Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                  ", "The model to compute with",                                                                     true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date             ", "Date which the option is paid (= start of the swap rate)",                                      true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date   ", "Tenor or End Date of the Cap, a tenor is nb followed by d/w/m/y like 6m or 1y",                 true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Pay Date				", "Date which the option is paid",																	true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                 ", "Sum Option's strike in real terms",                                                             true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Cap/Floor              ", "Option's type CAP or FlOOR",                                                                    true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency Period       ", "Cpn Frequency Period = nb followed by d/w/m/y like 6m or 1y",									true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Reset Gap              ", "Number of days between fixing and start date",											        true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing           ", "Reset in-advance (ADV) or in-arrears (ARR), default is ADV",							        true,       &DefaultResetTimingStruct   },
	{ GFAT_MATURITY_TYPE,           false,	    false,		"Index Term             ", "Index term : a nb followed by d/w/m/y like 6m or 1y",                            		        true,       &DefaultNothing             },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Index Days count       ", "Index Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		    true,       &DefaultPerCcyCharStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,	    false,		"Coeff                  ", "Sum Coefficient",																				true,       &DefaultDbleOneStruct		},
	{ GFAT_DOUBLE_TYPE,				false,	    false,		"First Idx              ", "First Idx",																				true,       &DefaultDbleMinusOneStruct	},
	{ GFAT_DOUBLE_TYPE,				false,	    false,		"Last Idx               ", "Last Idx",																				true,       &DefaultDbleMinusOneStruct	},

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SpreadOptionArg[] =
{
	/// Arg type			        Required?	VaArg?		Name			    	   Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                  ", "The model to compute with",                                                                     true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date             ", "Start date of the spread options schedule",														true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date   ", "Tenor or End Date of the spread options schedule, a tenor is nb followed by d/w/m/y like 6m or 1y",true                                 },
	{ GFAT_MATURITY_TYPE,	        true,		false,		"Index #1 Tenor         ", "Tenor of 1st index",																		    true                                    },
	{ GFAT_MATURITY_TYPE,	        true,		false,		"Index #2 Tenor         ", "Tenor of 2nd index",																		    true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    true,		false,		"Index #1 Gearing       ", "Gearing(s) of 1st index",																	    true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    true,		false,		"Index #2 Gearing       ", "Gearing(s) of 2nd index",																	    true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    true,		false,		"Strike                 ", "Spread option strike(s) in real terms",                                                         true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Cap/Floor              ", "Option's type CAP or FlOOR",                                                                    true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Cpn Frequency          ", "Cpn Frequency Period = nb followed by d/w/m/y like 6m or 1y",									true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Cpn Days count         ", "Coupon Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Reset Gap              ", "Number of days between fixing and start date",											        true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Pay Gap                ", "Number of days between end and payment date",										            true,       &DefaultDbleZeroStruct      },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    false,		false,		"Option's notional      ", "Option's notional(s) in real terms, default is 1",												true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing           ", "Reset in-advance (ADV) or in-arrears (ARR)",                                                    true,       &DefaultResetTimingStruct   },
	{ GFAT_STRING_TYPE,		        false,		false,		"Index #1 Type          ", "Type of the 1st index (LIBOR/CMS)",                                                             true,		&DefaultIndexTypeCms		},
	{ GFAT_STRING_TYPE,		        false,		false,		"Index #2 Type          ", "Type of the 2nd index (LIBOR/CMS)",                                                             true,		&DefaultIndexTypeCms		},
	{ GFAT_STRING_TYPE,		        false,		false,		"Payment Frequency      ", "Payment Frequency = nb followed by d/w/m/y like 6m or 1y or of type A, S, Q etc...",            true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Previous Leverage      ", "Previous Leverage = leverage on the previous coupon",											true,       &DefaultDbleZeroStruct		},
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CorridorArg[] =
{
	/// Arg type			        Required?	VaArg?		Name			    	   Description												                                        Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                      ", "The model to compute with",                                                                     true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date                 ", "Start date of the corridor schedule",															true									}, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date       ", "Tenor or End Date of the corrdior schedule, a tenor is nb followed by d/w/m/y like 6m or 1y",   true									},
	{ GFAT_STRING_TYPE,		        true,		false,		"Pay/rec                    ", "Corridor's type PAY or RCV",																	true									},
	{ GFAT_STRING_TYPE,		        true,		false,		"Payment Index Type         ", "Index of the paid rate (FIXED/LIBOR/CMS)",                                                      true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	true,		false,		"Fix Value                  ", "Value of the fixed paid index",																	true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	true,		false,		"Payment Index Mult         ", "Value of the PayIndexMult when the the paid rate is floating.",									true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Payment Frequency          ", "Payment Frequency (A/S/Q/M/W/D/Z)",																true,									},
	{ GFAT_STRING_TYPE,		        false,		false,		"Payment Index Term	        ", "Paid Index Term = nb followed by d/w/m/y like 6m or 1y",										true,       &DefaultPerCcyCharStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Payment Index Timing       ", "Paid Index Timing in-advance (ADV) or in-arrears (ARR)",										true,       &DefaultResetTimingStruct   },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Payment Day count          ", "Payment Day Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct	},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Reset Gap                  ", "Reset Gap",																						true,		&DefaultDbleZeroStruct      },
	{ GFAT_STRING_TYPE,		        false,		false,		"Int Rule                   ", "Interest Rule (ADJ/UNADJ)",																		true,		&DefaultIntRule             },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Spread Value               ", "Spread on the Payment Index",																	true,		&DefaultDbleZeroStruct      },
	{ GFAT_STRING_TYPE,		        true,		false,		"Fixing Frequency           ", "Fixing Frequency (A/S/Q/M/W/D/ZC)",																true,									},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Timing              ", "Fixing in-advance (ADV) or in-arrears (ARR), default is ARR",							        true,       &DefaultResetTimingStruct },
	{ GFAT_STRING_TYPE,		        true,		false,		"Fixing Index Type 1        ", "Type of the second fixing index (LIBOR/CMS)",													true                                    },
	{ GFAT_STRING_OR_CURVE_TYPE,    false,		false,		"Fixing Index Term 1        ", "Paid Index Term 1 = nb followed by d/w/m/y like 6m or 1y",										true,       &DefaultPerCcyCharStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Coeff 1 Value              ", "Value of the coeff for the first index",														true,       &DefaultDbleOneStruct		},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Barrier Down Value         ", "Value of the barrier down",																		true,       &DefaultDbleZeroStruct      },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Barrier Up Value           ", "Value of the barrier up",																		true,       &DefaultDbleInfiniteStruct  },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Notional Value             ", "Value of the notional",																			true,       &DefaultDbleOneStruct		},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Index Type 2        ", "Type of the first fixing index (LIBOR/CMS)",													true,		&DefaultIndexType			},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Index Term 2        ", "Paid Index Term 2 = nb followed by d/w/m/y like 6m or 1y",										true,       &DefaultPerCcyCharStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Coeff 2 Value              ", "Value of the coeff for the second fixing index",												true,       &DefaultDbleOneStruct		},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Index Type 3        ", "Type of the single rate condition index (LIBOR/CMS)",											true,		&DefaultFixedIndex			},
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Index Term 3        ", "Term of the single rate condition index = nb followed by d/w/m/y like 6m or 1y",				true,       &DefaultPerCcyCharStruct    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Barrier Down 3             ", "Value of the barrier down for the single rate condition",										true,       &DefaultDbleZeroStruct		},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Barrier Up 3               ", "Value of the barrier up for the single rate condition",											true,       &DefaultDbleInfiniteStruct  },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Date Strip Struct			", "Date strip for underlying structure",															true,       &DefaultDbleZeroStruct      },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Date Strip Pay				", "Date strip for Payment Index",																	true,       &DefaultDbleZeroStruct      },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Date Strip Fix				", "Date strip for Ref Index",																		true,       &DefaultDbleZeroStruct		},
	{ GFAT_VECTOR_OR_CURVE_TYPE,	false,		false,		"Index						", "Vector of indexes",																				true,       &DefaultDbleZeroStruct		},
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Date Strip Offset			", "Schedule offset in the date strip, default is 0",												true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Nb of flows				", "Nb of flows in the date strip, default is 0",													true,       &DefaultDbleZeroStruct      },	

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct DoubleDigitalArg[] = 
{
	/// Arg type			        Required?	VaArg?		Name						 Description																					Visible?    Default Value
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model                      ", "The model to compute with",																		true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"First Rate                 ", "PreComputed rate used for the first digital condition",											true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"First Strike Down          ", "Strike Down for the first digital range condition",												true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"First Strike Up            ", "Strike Up for the first digital range condition",												true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"Second Rate                ", "PreComputed rate used for the second digital condition",										true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"Second Strike Down         ", "Strike Down for the second digital range condition",											true},
	{ GFAT_VECTOR_TYPE,				true,		false,		"Second Strike Up           ", "Strike Up for the second digital range condition",												true},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"First Strike Spread        ", "Strike Spread (>0) for the first digital condition",											true,      &DefaultEpsilonStruct      },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Second Strike Spread       ", "Strike Spread (>0) for the second digital condition",											true,      &DefaultEpsilonStruct      },
	
	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct SpotArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description														        Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",	"The model to compute the spot value with",						        true                                    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Settlement Gap or Date ",	"Settlement date or number of days between expiry and settlement date", true,       &DefaultPerCcyDbleStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct FwdArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description														        Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",  "The model to compute the spot value with",						        true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Expiry Date            ",	"Expiry date of the forward",									        true                                    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Settlement Gap or Date ",	"Settlement date or number of days between expiry and settlement date", true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Payment Gap or Date    ",	"Payment date or number of days between expiry and payment date",	    true,       &DefaultPerCcyDbleStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description														        Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",  "The model to compute the spot value with",						        true                                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		true,		false,		"Expiry Date            ",	"Expiry Date of the option",									        true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                 ",	"Call strike in real terms",									        true					                },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put               ",  "Option's type Call or Put",									        true                                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		false,		false,		"Settlement Gap or Date ",	"Settlement date or number of days between expiry and settlement date", true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		false,		false,		"Payment Gap or Date    ",	"Payment date or number of days between expiry and payment date",	    true,       &DefaultPerCcyDbleStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct GreekArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description														        Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",  "The model to compute the spot value with",						        true                                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		true,		false,		"Expiry Date            ",	"Expiry Date of the option",									        true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                 ",	"Call strike in real terms",									        true					                },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put               ",  "Option's type Call or Put",									        true                                    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		false,		false,		"Settlement Gap or Date ",	"Settlement date or number of days between expiry and settlement date", true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_VECTOR_TYPE,		false,		false,		"Payment Gap or Date    ",	"Payment date or number of days between expiry and payment date",	    true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Greek Type             ",	"Delta, Vega, etc...",													true									},

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct CallSpreadArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description																				Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model1                  ", "The model1 to compute the spot1 value with",											  true                                    },
	{ GFAT_MODEL_TYPE,				true,		false,		"Model2                  ", "The model2 to compute the spot2 value with",											  true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Expiry Date             ",	"Expiry Date of the option",															  true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                  ",	"CallSpread strike in real terms",														  true					                  },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Alpha					 ",	"Real Coef of the Spot1",																  true					                  },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Beta					 ",	"Real Coef of the Spot2",																  true					                  },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Settlement1 Gap or Date ",	"Settlement date or number of days between expiry and settlement date for the first FX",  true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Settlement2 Gap or Date ",	"Settlement date or number of days between expiry and settlement date for the second FX", true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Payment Gap or Date     ",	"Payment date or number of days between expiry and payment date",						  true,       &DefaultPerCcyDbleStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct EqDigArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description														            Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",  "The model to compute the spot value with",						        true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Expiry Date            ",	"Expiry Date of the forward",									        true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike                 ",	"Call strike in real terms",									        true					                },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put               ",  "Option's type CALL or PUT (Digital Case)",						        true                                    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Notional               ",	"Notional of the digital",										        true,		&DefaultDbleOneStruct       },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Settlement Gap or Date ",	"Settlement date or number of days between expiry and settlement date", true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Payment Gap or Date    ",	"Payment date or number of days between expiry and payment date",	    true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"DigitalType            ",	"ANALYTIC, CENTRED, BACKWARD, FORWARD",									false,      &DefaultDigitTypeStruct     },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Epsilon                ",	"Epsilon for the callspread",									        false,		&DefaultEpsilonStruct	    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct CallStripArg[] =
{
	/// Arg type					Required?	VaArg?		Name					    Description														                                Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Model                  ",  "The model to compute the option strip value with",						                        true                                    },
	{ GFAT_DATE_TYPE,	  	        true,		false,		"Start Date             ",  "Start date of the option strip",                                                               true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date   ",  "Maturity or End Date of the option strip, a maturity is nb followed by d/w/m/y like 6m or 1y", true                                    },
	{ GFAT_VECTOR_OR_CURVE_TYPE,	true,		false,		"Strike                 ",	"Call strike in real terms",									                                true					                },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put               ",  "Option's type Call or Put",									                                true                                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Frequency              ",  "Reset Frequency = nb followed by d/w/m/y like 6m or 1y",									    true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Reset Gap              ",  "Number of days between reset and settlement date",											    true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Pay Gap                ",  "Number of days between settlement and payment date",										    true,       &DefaultDbleZeroStruct      },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    false,		false,		"Notional               ",  "Notional in real terms, default is 1",											                true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing           ",  "Reset in-advance (ADV) or in-arrears (ARR), default is ADV",							        true,       &DefaultResetTimingStruct   },
	{ GFAT_VECTOR_OR_CURVE_TYPE,    false,		false,		"Leverage               ",  "Leverage, default is 1",											                            true,       &DefaultDbleOneStruct       },
	{ GFAT_DATESTRIP_TYPE,		    false,		false,		"Date Strip             ",  "Date strip describing the schedule of reset, settlement & payment dates",				        true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Date Strip Offset      ",  "Schedule offset in the date strip, default is 0",										        true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Date Strip Width		",  "Schedule offset in the date strip, default is -1",												true,      &DefaultDbleMinusOneStruct  },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct RangeAccrualArg[] =
{
	/// Arg type					Required?	VaArg?		Name							Description																						Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"IRModel					",  "The model to compute IRindex and PAYindex with",												true,                                   },//0
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date					",	"Start Date of the Range Accrual schedule",														true,                                   },//1
	{ GFAT_DATE_TYPE,		        true,		false,		"End Date					",	"End Date of the Range Accrual schedule",														true,                                   },//2
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Payment Gap or Date		",	"Payment date or number of days between end date and payment date",								true,       &DefaultPerCcyDbleStruct    },//3
	{ GFAT_MULTITOKENSTRING_TYPE,   false,      false,		"Payment Day count          ",  "Payment Day Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		true,       &DefaultPerCcyCharStruct	},//4
	{ GFAT_STRING_TYPE,		        true,		false,		"PAYindexType				",  "Index of the paid rate (FIXED/LIBOR/CMS)",                                                     true,								    },//5
	{ GFAT_STRING_TYPE,		        false,		false,		"PAYindexTerm		        ",  "Paid Index Term = nb followed by d/w/m/y like 6m or 1y",										true,       &DefaultPerCcyCharStruct    },//6
	{ GFAT_STRING_TYPE,		        false,		false,		"PAYindexTiming				",  "Paid Index Timing in-advance (ADV) or in-arrears (ARR)",										true,       &DefaultResetTimingStruct   },//7
	{ GFAT_DOUBLE_TYPE,				false,		false,		"PAYindexSpread             ",  "Spread on the Payment Index",																	true,		&DefaultDbleZeroStruct      },//8
	{ GFAT_STRING_TYPE,		        true,		false,		"Fixing Frequency			",  "Fixing Frequency (A/S/Q/M/W/D/ZC)",															true,                                   },//9
	{ GFAT_MODEL_TYPE,		        true,		false,		"FXModel                    ",  "The model to compute the FX with",																true,                                   },//10
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Barrier Down Value for FX  ",	"Value of the barrier down for the FX",															true,		&DefaultDbleZeroStruct      },//11
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Barrier Up Value for FX	",	"Value of the barrier up for the FX",															true,		&DefaultDbleInfiniteStruct  },//12
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Notional Value				",	"Notional in real terms",																		true,		&DefaultDbleOneStruct	    },//13
	{ GFAT_STRING_TYPE,		        true,		false,		"IRFixing Index Type        ", "Type of the ir fixing index (LIBOR/CMS)",														true,                                   },//14
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixing Index Term 1        ",  "Paid Index Term 1 = nb followed by d/w/m/y like 6m or 1y",										true,       &DefaultPerCcyCharStruct    },//15
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Barrier Down Value for IR  ",	"Value of the barrier down for the IR index",													true,		&DefaultDbleZeroStruct		},//16
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Barrier Up Value for IR	",	"Value of the barrier up for the IR index",														true,		&DefaultDbleInfiniteStruct	},//17

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct ModelFactorArg[] =
{
	/// Arg type					Required?	VaArg?		Name				Description										Visible?    Default Value
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Model        ",	"Returns a model factor",					     true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CFCallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Fwd               ",	"forward of the call",					    true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike            ",	"strike",									true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Volatility        ",	"volatility",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Time to Expiry    ",	"time to expiry",							true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Discount Factor   ",	"discount factor",							true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put          ",	"Option's type Call or Put",				true                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Distribution Type ",	"Normal, Lognormal",						true,	&DefaultLN      },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CFGreekArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Fwd               ",	"forward of the call",					    true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike            ",	"strike",									true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Volatility        ",	"volatility",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Time to Expiry    ",	"time to expiry",							true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Discount Factor   ",	"discount factor",							true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Call/Put          ",	"Option's type Call or Put",				true                    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Distribution Type ",	"Normal, Lognormal",						true,		&DefaultLN  },
	{ GFAT_STRING_TYPE,		        true,		false,		"Greek Type        ",	"Delta, Vega, etc...",						true                    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct ImpliedVolArg[] =
{
	{ GFAT_MODEL_TYPE,		        true,		false,		"Model               ", "The model to compute the capfloor with",                                                       true                                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date          ", "Start Date of the cap/floor",                                                                  true                                    }, 
	{ GFAT_DATEORMATU_TYPE,	        true,		false,		"Maturity or End Date", "Maturity or End Date of the cap/floor : a maturity is nb followed by d/w/m/y like 6m or 1y",   true                                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Strike              ", "Option's strike in real terms",                                                                true                                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Cap/Floor           ", "Option's type CAP or FlOOR",                                                                   true                                    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Cpn Days count      ", "Coupon Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Reset Gap or Date   ",	"Reset date or number of days between fixing and start date",									true,       &DefaultPerCcyDbleStruct    },
	{ GFAT_DATE_OR_DOUBLE_TYPE,		false,		false,		"Pay Gap  or Date    ",	"Payment date or number of days between start date(ADV)/end date(ARR) and payment date",		true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Option's notional   ", "Option's notional in real terms, default is 1",												true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Reset Timing        ", "Reset in-advance (ADV) or in-arrears (ARR), default is ADV",							        true,       &DefaultResetTimingStruct   },
	{ GFAT_MATURITY_TYPE,           false,	    false,		"Index Term          ", "Index term : a nb followed by d/w/m/y like 6m or 1y",                            		        true,       &DefaultNothing             },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Index Days count    ", "Index Days Count Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",		    true,       &DefaultPerCcyCharStruct    },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


///////////////////////////////////////////////////////////////////////////


////////////////// Option on Funds keywords ///////////////////////////////


///////////////////////////////////////////////////////////////////////////


static ARM_GramFunctionArgStruct CF_MepiCallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Spot               ",	"Initial value of the portfolio",			true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike             ",	"strike of the option",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Relative Maturity  ",	"Relative Maturity of the option",			true                    },		
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Risk Factor        ",	"Risk Factor",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Min Exposure       ",	"Exposure minimum",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Max Exposure       ",	"exposure Maximum",							true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Initial Protection ",	"Initial protection",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Final Protection   ",	"Final Protection (1.0 for capital garanty)",true                   },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Borrowing Spread   ",	"Borrowing spread",							true                    },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Yearly Fees        ",	"Yearly Fees",								true,	&DefaultLN      },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Cash Spread        ",	"Cash Spread							  ",true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Period Nb          ",	"Number of rehedging periods",				true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Avg Period Nb      ",	"Number of periods for averaging",			true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Asian Reset		",	"Already Computed Average				  ",true                    },
	{ GFAT_STRING_TYPE,				true,		false,		"Call/Put           ",	"Call or Put",								true                    },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Portef Min          ",	"Lower boundary for the portfolio",			true				    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Portef Max         ",	"High boundary for th portfolio",			true                    },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Port Discret Size  ",	"Size of the portfolio vector",				true ,&DefaultDbleFiftyStruct },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Hermite Size       ",	"Number of Integration points",				true ,&DefaultDbleTwentyStruct},
	{ GFAT_STRING_TYPE,		        true,		false,		"Curve Name         ",	"Currency Name                        ",    true					},
	{ GFAT_STRING_TYPE,		        true,		false,		"Equity Name        ",	"Equity Name                          ",	true					},
	{ GFAT_STRING_TYPE,		        true,		false,		"Algorithm Type     ",	"Type of Algorithm or model",				true                    },
	

	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CF_MepiGreekCallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Spot               ",	"Initial value of the portfolio",				true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike             ",	"strike of the option",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Relative Maturity  ",	"Relative Maturity of the option",			true                    },		
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Risk Factor        ",	"Risk Factor",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Min Exposure       ",	"Exposure minimum",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Max Exposure       ",	"exposure Maximum",							true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Initial Protection ",	"Initial protection",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Final Protection   ",	"Final Protection (1.0 for capital garanty)",true                   },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Borrowing Spread   ",	"Borrowing spread",							true                    },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Yearly Fees        ",	"Yearly Fees",								true,	&DefaultLN      },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Cash Spread        ",	"Cash Spread							  ",true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Period Nb          ",	"Number of rehedging periods",				true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Avg Period Nb      ",	"Number of periods for averaging",			true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Asian Reset		",	"Already Computed Average				  ",true                    },
	{ GFAT_STRING_TYPE,				true,		false,		"Call/Put           ",	"Call or Put",								true                    },
	{ GFAT_DOUBLE_TYPE,		        true,		false,		"Portef Min          ",	"Lower boundary for the portfolio",			true ,&DefaultDbleZeroStruct},
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Portef Max         ",	"High boundary for the portfolio",			true  ,&DefaultDbleTwoStruct},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Port Discret Size  ",	"Size of the portfolio vector",				true ,&DefaultDbleFiftyStruct },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Hermite Size       ",	"Number of Integration points",				true ,&DefaultDbleTwentyStruct},
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Shift Size       ",	"Relative Shift for the numerical derivation",true ,&DefaultDbleOnePerCentStruct},
	{ GFAT_STRING_TYPE,		        true,		false,		"Curve Name         ",	"Currency Name                        ",    true					},
	{ GFAT_STRING_TYPE,		        true,		false,		"Equity Name        ",	"Equity Name                          ",	true					},
	{ GFAT_STRING_TYPE,		        true,		false,		"Algorithm Type     ",	"Type of Algorithm or model",				true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Greek Type        ",	"Delta, Vega, etc...",						true                    },
	

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct NumericalCallOnMepiDeltaArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_STRING_TYPE,				true,		false,		"Curve Name         ",  "Discount Curve Name",						true                    },
	{ GFAT_STRING_TYPE,				true,		false,		"Equity Model Name  ",	"Name of the Equity Model",					true                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date         ",	"Start Date of the Call On Mepi",			true					},
	{ GFAT_DATE_TYPE,		        true,		false,		"End Date           ",	"End Date of the Call On Mepi",				true					},
	{ GFAT_STRING_TYPE,		        true,		false,		"Reset Frequency    ",  "Fund reset Period = nb followed by d/w/m/y like 6m or 1y",	true	},
	{ GFAT_DOUBLE_TYPE,				true,		false,		"RiskFactor         ",	"RiskFactor",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Strike		        ",	"Strike of the Call",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"MaxBorrow          ",	"MaxBorrow",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"ProtectionCurveStart",	"ProtectionCurveStart",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"ProtectionCurveEnd ",	"ProtectionCurveEnd",						true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"StartingPortFolio  ",	"StartingPortfolio",						true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"StartingCash       ",	"StartingCash",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Min Invested       ",	"Min Invested",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Leverage Cost      ",	"Leverage Cost",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"CashSpread	        ",	"CashSpread",								true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Fees    	        ",	"Fees",										true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"AlreadyAsianed	    ",	"AlreadyAsianed",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"AsianingPeriodNb   ",	"AsianingPeriodNb",							true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"NbSimuls           ",	"NbSimuls For Mc",							true		            },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Pf Bump           ",	"Pf Bump",									true		            },


	/// says that this is the end
	{ GFAT_TERMINATOR }
};


static ARM_GramFunctionArgStruct CF_StochVolCallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Fwd               ",	"forward price of the underlying",			true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike             ",	"strike of the option",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Maturity           ",	"Maturity of the option",					true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"ZC Rate            ",	"Zero coupon rate at maturity of the option",true                   },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Init Volatility    ",	"Initial Volatility Stae of the underlying",true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Vol Drift          ",	"Drift of the Volatility",						true,	&DefaultLN      },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Vol of Vol         ",	"Volatility of the Volatility",					    true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Avg Period	        ",	"averaging period",							true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Reset				",	"reset value",								true                    },
	{ GFAT_STRING_TYPE,				true,		false,		"Call/Put           ",	"Call or Put",								true                    },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Hermite Size       ",	"Number of Integration points",				true ,&DefaultDbleTwentyStruct},
	{ GFAT_STRING_TYPE,		        true,		false,		"Algorithm Type     ",	"Type of Algorithm or model",				true                    },
	

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct CF_StochVolGreekCallArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description									Visible?    Default Value
	{ GFAT_VECTOR_TYPE,				true,		false,		"Fwd               ",	"forward price of the underlying",			true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Strike             ",	"strike of the option",						true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Maturity           ",	"Maturity of the option",					true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"ZC Rate            ",	"Zero coupon rate at maturity of the option",true                   },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Init Volatility    ",	"Initial Volatility Stae of the underlying",true                    },
	{ GFAT_VECTOR_TYPE,		        true,		false,		"Vol Drift          ",	"Drift of the Volatility",						true,	&DefaultLN      },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Vol of Vol         ",	"Volatility of the Volatility",					    true                    },
	{ GFAT_DOUBLE_TYPE,				true,		false,		"Avg Period	        ",	"averaging period",							true                    },
	{ GFAT_VECTOR_TYPE,				true,		false,		"Reset				",	"reset value",								true                    },
	{ GFAT_STRING_TYPE,				true,		false,		"Call/Put           ",	"Call or Put",								true                    },
	{ GFAT_DOUBLE_TYPE,				false,		false,		"Hermite Size       ",	"Number of Integration points",				true ,&DefaultDbleTwentyStruct},
	{ GFAT_STRING_TYPE,		        true,		false,		"Algorithm Type     ",	"Type of Algorithm or model",				true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"Greek Type        ",	"Delta, Vega, etc...",						true                    },
	

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

/////////////////////////////////////////////////////////
////////////////// Inflation keywords //////////////////
/////////////////////////////////////////////////////////

static ARM_GramFunctionArgStruct CPIArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description														Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"InfModel          ",	"The model to compute the CPI spot value with",					true									},
	{ GFAT_DATE_TYPE,		        true,		false,		"CPIDate           ",	"CPIDate",														true									},
	{ GFAT_MATURITY_TYPE,	        false,		false,		"DCFLag            ",	"Difference between CPIDate and CPIFixing Date",				true,		&DefaultThreeMonths			},
	{ GFAT_STRING_TYPE,		        false,		false,		"Daily Interp Type ",	"The way to do daily interpolation",							true,		&DefaultPerCcyCharStruct    },
	{ GFAT_MATURITY_TYPE,			false,		false,		"Reset Lag         ",	"Reset Lag of the CPI",											true,       &DefaultZeroMonths	        },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct InfSwapRateArg[] =
{
	/// Arg type					Required?	VaArg?		Name							Description																Visible?    Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Currency                   ",	"The IR model to compute the the SwapRate with",						true                    },
	{ GFAT_MODEL_TYPE,				true,		false,		"InfModel                   ",	"The Inf model to compute the SwapRate with",							true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"YoY/OAT                    ",  "InfSwap's type YoY or OAT",											true                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date                 ",	"Start Date of the swap",												true                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"End Date                   ",	"End Date of the swap ",												true                    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Spread                     ",	"Spread",																true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Floating's Reset Gap       ",	"Reset Gap (in days) of the floating",									true,       &DefaultDbleOneStruct       },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Floating's Payment Gap     ",	"Payment Gap (in days) of the floating",								true,		&DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Floating Frequency Period  ",	"Floating Leg Period = nb followed by d/w/m/y like 6m or 1y",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Floating Days count        ",	"Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",	true,       &DefaultPerCcyCharStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixed Frequency Period     ",	"Fixed Leg Period = nb followed by d/w/m/y like 6m or 1y",				true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Fixed Days count           ",	"Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",	true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Coupon                     ",	"Coupon of the OAT type",												true,       &DefaultDbleOneStruct       },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};

static ARM_GramFunctionArgStruct InfSwapArg[] =
{
	/// Arg type					Required?	VaArg?		Name					Description																				Visible?		Default Value
	{ GFAT_MODEL_TYPE,				true,		false,		"Currency                   ",	"The IR model to compute the the SwapRate with",								true                    },
	{ GFAT_MODEL_TYPE,				true,		false,		"InfModel                   ",	"The Inf model to compute the SwapRate with",									true                    },
	{ GFAT_STRING_TYPE,		        true,		false,		"YoY/OAT                    ",  "InfSwap's type YoY or OAT",													true                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"Start Date                 ",	"Start Date of the swap",														true                    },
	{ GFAT_DATE_TYPE,		        true,		false,		"End Date                   ",	"End Date of the swap ",														true                    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Strike                     ",	"Swap's Strike in real terms",													true,       &DefaultDbleOneStruct       },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Spread                     ",	"Spread",																		true,       &DefaultDbleZeroStruct      },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Floating's Reset Gap       ",	"Reset Gap (in days) of the floating",											true,       &DefaultDbleOneStruct       },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Floating's Payment Gap     ",	"Payment Gap (in days) of the floating",										true,       &DefaultDbleOneStruct       },
	{ GFAT_STRING_TYPE,		        false,		false,		"Floating Frequency Period  ",	"Floating Leg Period = nb followed by d/w/m/y like 6m or 1y",					true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Floating Days count        ",	"Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_STRING_TYPE,		        false,		false,		"Fixed Frequency Period     ",	"Fixed Leg Period = nb followed by d/w/m/y like 6m or 1y",						true,       &DefaultPerCcyCharStruct    },
	{ GFAT_MULTITOKENSTRING_TYPE,   false,	    false,		"Fixed Days count           ",	"Convention ACTUAL, A365, A360, 30/360, ACTREAL, ACT29, ISMA, NOBASE",			true,       &DefaultPerCcyCharStruct    },
	{ GFAT_DOUBLE_TYPE,		        false,		false,		"Coupon                     ",	"Coupon of the OAT type",														true,       &DefaultDbleOneStruct       },

	/// says that this is the end
	{ GFAT_TERMINATOR }
};



/// The table of all functions.
/// should be revisited when completing the framework


/// currently missing function definition!!!
const ARM_GramFunctionStruct FunctionsTable[] =
{
	/// Financial functions
	/// interest rates functions
	/// Name			Return type			Description				Arguments		Functor			Category		Special GramNode Builder
	{ "Annuity",		GFAT_VECTOR_TYPE,	"Annuity (PVBP)",		AnnuityArg,		&AnnuityFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "Cap",			GFAT_VECTOR_TYPE,	"Cap",				    CapArg,			&CapFctor,	    "Financial",	&GramNodeFuncBuilder			},
	{ "Caplet",			GFAT_VECTOR_TYPE,	"Cap/Floor-let Option",	CapletArg,		&CapletFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "Digital",		GFAT_VECTOR_TYPE,	"Digital Cap/Floor-let",CapletArg,		&DigitalFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "DF",				GFAT_VECTOR_TYPE,	"Discount factor",		DFArg,			&DFFctor,		"Financial",	&GramNodeFuncBuilder			},
	{ "Libor",			GFAT_VECTOR_TYPE,	"Libor Rate",			LiborArg,		&LiborFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "SwapRate",		GFAT_VECTOR_TYPE,	"Swap Rate ",			SwapRateArg,	&SwapRateFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "SpreadCMS",		GFAT_VECTOR_TYPE,	"Spread ",				SpreadArg,		&SpreadFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "Swap",			GFAT_VECTOR_TYPE,	"Swap",				    SwapArg,		&SwapFctor,	    "Financial",	&GramNodeFuncBuilder			},
	{ "BasisSwap",		GFAT_VECTOR_TYPE,	"BasisSwap",			BasisSwapArg,	&BasisSwapFctor,"Financial",	&GramNodeFuncBuilder			},
	{ "Swaption",		GFAT_VECTOR_TYPE,	"Swaption",				SwaptionArg,	&SwaptionFctor,	"Financial",  	&GramNodeFuncBuilder			},
	{ "SumOption",		GFAT_VECTOR_TYPE,	"SumOption",			SumOptionArg,	&SumOptionFctor,"Financial",  	&GramNodeFuncBuilder			},
	{ "SpreadOption",	GFAT_VECTOR_TYPE,	"SpreadOption",			SpreadOptionArg,&SpreadOptionFctor,"Financial",	&GramNodeFuncBuilder			},
	{ "Corridor",		GFAT_VECTOR_TYPE,	"Corridor",				CorridorArg,	&CorridorFctor,"Financial",		&GramNodeFuncBuilder			},
	{ "DoubleDigital",	GFAT_VECTOR_TYPE,	"Double Digital",		DoubleDigitalArg,&DoubleDigitalFctor,"Financial",&GramNodeFuncBuilder			},
	{ "MinMax",			GFAT_VECTOR_TYPE,	"MinMax ",				MaxRateArg,		&MaxRateFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "ImpliedVol",		GFAT_VECTOR_TYPE,	"ImpliedVol",			ImpliedVolArg,	&ImpliedVolFctor, "Financial",	&GramNodeFuncBuilder			},

	/// inflation functions
	{ "CPI",			GFAT_VECTOR_TYPE,	"CPI",					CPIArg,			&CPIFctor,		"Financial",  	&GramNodeFuncBuilder			},
	{ "InfSwapRate",	GFAT_VECTOR_TYPE,	"Inflation Swap Rate",	InfSwapRateArg,	&InfSwapRateFctor,"Financial",  	&GramNodeFuncBuilder		},
	{ "InfSwap",		GFAT_VECTOR_TYPE,	"Inflation Swap",		InfSwapArg,		&InfSwapFctor,	"Financial",  	&GramNodeFuncBuilder			},

	/// equity part
	{ "Call",			GFAT_VECTOR_TYPE,	"Call or put",			CallArg,		&CallFctor,	        "Financial",	&GramNodeFuncBuilder			},
	{ "CallSpread",		GFAT_VECTOR_TYPE,	"CallSpread",			CallSpreadArg,	&CallSpreadFctor,   "Financial",	&GramNodeFuncBuilder			},
	{ "EqDigital",		GFAT_VECTOR_TYPE,	"Equity digital",		EqDigArg,		&EqDigitalFctor,    "Financial",	&GramNodeFuncBuilder			},
	{ "Fwd",			GFAT_VECTOR_TYPE,	"Fwd",				    FwdArg,			&FwdFctor,	        "Financial",	&GramNodeFuncBuilder			},
	{ "Spot",			GFAT_VECTOR_TYPE,	"Spot",				    SpotArg,		&SpotFctor,	        "Financial",	&GramNodeFuncBuilder			},
	{ "CallStrip",		GFAT_VECTOR_TYPE,	"Option strip",		    CallStripArg,   &CallStripFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "RangeAccrual",	GFAT_VECTOR_TYPE,	"RangeAccrual",		    RangeAccrualArg,&RangeAccrualFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "Greek",			GFAT_VECTOR_TYPE,	"Greek",				GreekArg,		&GreekFctor,	    "Financial",	&GramNodeFuncBuilder			},


	/// Payment function
	{ "LinearPV",		GFAT_VECTOR_TYPE,	"Linear Present value",	PVArg,		&PVFctor,		"Financial",	&GramNodeLinearPVFuncBuilder	},
	/// Name			Return type			Description				Arguments		Functor			Category
	{ "Exercise",		GFAT_VECTOR_TYPE,	"Exercise",				ExerciseArg,	&ExerciseFctor,	"Financial",	&GramNodeExerciseFuncBuilder	},
	{ "Frontier",		GFAT_VECTOR_TYPE,	"Frontier",				FrontierArg,	&FrontierFctor,	"Financial",	&GramNodeFuncBuilder	},
	{ "PV",				GFAT_VECTOR_TYPE,	"Present value",		PVArg,			&PVFctor,		"Financial",	&GramNodeFuncBuilder			},
	{ "Trigger",		GFAT_VECTOR_TYPE,	"Trigger",				TriggerArg,		&TriggerFctor,	"Financial",	&GramNodeTriggerFuncBuilder		},
	{ "UnPay",			GFAT_VECTOR_TYPE,	"Remove discounting",	UnPayArg,		&UnPayFctor,	"Financial",	&GramNodeFuncBuilder			},

	/// Closed formed function
	{ "CFCall",			GFAT_VECTOR_TYPE,	"Closed Form Call",		CFCallArg,		&CFCallFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "CFGreek",		GFAT_VECTOR_TYPE,	"Closed Form Greek",	CFGreekArg,		&CFGreekFctor,	"Financial",	&GramNodeFuncBuilder			},


		/// Option on Mepi Funds function
	{ "CFMepiCall",		GFAT_VECTOR_TYPE,	"Closed Form Mepi Option",		CF_MepiCallArg,		&CF_MepiCallFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "CFMepiCallGreek",	GFAT_VECTOR_TYPE,	"Closed Form Mepi Greek Option",CF_MepiGreekCallArg,&CF_MepiGreekFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "NumericalCallOnMepiDelta", GFAT_VECTOR_TYPE,	"Pricing of Call On Mepi with MC", NumericalCallOnMepiDeltaArg, &NumericalCallOnMepiDeltaFctor,	"Financial",	&GramNodeFuncBuilder			},

		/// Option on Stochastic Volatility Funds function
	{ "CFStochVolCall",		GFAT_VECTOR_TYPE,	"Closed Form Stochastic Vol Option",		CF_StochVolCallArg,		&CF_StochVolCallFctor,	"Financial",	&GramNodeFuncBuilder			},
	{ "CFStochVolCallGreek",	GFAT_VECTOR_TYPE,	"Closed Form Stochastic Vol Greek Option",CF_StochVolGreekCallArg,&CF_StochVolGreekFctor,	"Financial",	&GramNodeFuncBuilder			},

	/// bond :yield function
	{ "PriceToYield",	GFAT_VECTOR_TYPE,	"Price to yield.",		PTYArg,		&PTYFctor,		"Financial",	&GramNodeFuncBuilder			},
	{ "YieldToPrice",	GFAT_VECTOR_TYPE,	"Yield To Price.",		YTPArg,		&YTPFctor,		"Financial",	&GramNodeFuncBuilder			},

	/// General model functions
	{ "ModelFactor",	GFAT_VECTOR_TYPE,	"Model Factor",			ModelFactorArg,	&MFactorFctor,	"Financial",	&GramNodeFuncBuilder			},

	/// Conditionals															
	/// Name			Return type			Description				Arguments	Functor			Category
	{ "If",				GFAT_VECTOR_TYPE,	"If/Then/Else",			IfArg,		&IfFctor,		"Conditional",	&GramNodeFuncBuilder			},
	{ "Min",			GFAT_VECTOR_TYPE,	"Minimum of two values",MinMaxArg,	&MinFctor,		"Conditional",	&GramNodeFuncBuilder			},
	{ "Max",			GFAT_VECTOR_TYPE,	"Maximum of two values",MinMaxArg,	&MaxFctor,		"Conditional",  &GramNodeFuncBuilder			},

	/// Date Algebra			
	/// Name			Return type			Description				Arguments	Functor			Category
	{ "DCF",			GFAT_DOUBLE_TYPE,	"Day Count Fraction",	DCFArg,		&DCFFctor,		"Dates", 		&GramNodeFuncBuilder			},

	/// Name			Return type			Description				Arguments	Functor			Category
	{ "Exp",			GFAT_VECTOR_TYPE,	"Exponential",			DoubleArg,	&ExpVecFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "Log",			GFAT_VECTOR_TYPE,	"Natural Logarithm",	DoubleArg,	&LogVecFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "Sqrt",			GFAT_VECTOR_TYPE,	"Square Root",			DoubleArg,	&SqrtVecFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "Pow",			GFAT_VECTOR_TYPE,	"Power (x^y) with the second argument being a deterministic double!",
																	PowArg,		&PowFctor,		"Mathematical", &GramNodeFuncBuilder			},
	{ "CumNorm",		GFAT_VECTOR_TYPE,	"Cumulative Normal!",	DoubleArg,	&CumNormFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "Abs",			GFAT_VECTOR_TYPE,	"Absolute value",		DoubleArg,	&AbsVecFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "SumSerie",		GFAT_VECTOR_TYPE,	"Absolute value",		SumSerieArg, &SumSerieFctor,	"Mathematical", &GramNodeFuncBuilder			},
	{ "StatAverage",	GFAT_VECTOR_TYPE,	"Statistics Average",	StatAverageArg, &StatAverageFctor,	"Mathematical", &GramNodeFuncBuilder	},
	{ "StatStdDev",		GFAT_VECTOR_TYPE,	"Statistics Standard Deviation",	StatStdDevArg, &StatStdDevFctor,	"Mathematical", &GramNodeFuncBuilder	},
};

const size_t FuncTableSize = sizeof(FunctionsTable)/sizeof(FunctionsTable[0]);


/// this is only for the help of the grammar function since these functions have no functor
/// and ARE BUILT-IN FUNCTION OF THE LANGUAGE!

const ARM_GramFunctionStruct BuiltInFunctionsTable[] =
{
	/// builtin token do not have functor since the functor called is done directly in the grammar parsing!
	/// Name			Return type			Description				Arguments	Functor			Category			
	{ "+",				GFAT_VECTOR_TYPE,	"Plus",					BinArg,		NULL,			"Arithmetic" },
	{ "-",    			GFAT_VECTOR_TYPE,	"Minus",				BinArg,		NULL,			"Arithmetic" },
	{ "*",				GFAT_VECTOR_TYPE,	"Multiplies",			BinArg,		NULL,			"Arithmetic" },
	{ "/",				GFAT_VECTOR_TYPE,	"Divides",				BinArg,		NULL,			"Arithmetic" },
	{ "%",				GFAT_VECTOR_TYPE,	"Modulo (casts the 2 arguments in integer!)",				
																	BinArg,		NULL,			"Arithmetic" },   
	{ "==",    			GFAT_VECTOR_TYPE,	"Equal to",				BinArg,		NULL,			"Comparison" },
	{ "!=",				GFAT_VECTOR_TYPE,	"Not Equal to",			BinArg,		NULL,			"Comparison" },
	{ "<",				GFAT_VECTOR_TYPE,	"Less",					BinArg,		NULL,			"Comparions" },
	{ ">",				GFAT_VECTOR_TYPE,	"Greater",				BinArg,		NULL,			"Comparions" },
	{ "<=",				GFAT_VECTOR_TYPE,	"Less or Equal",		BinArg,		NULL,			"Comparions" },  
	{ ">=",    			GFAT_VECTOR_TYPE,	"Greater or Equal",		BinArg,		NULL,			"Comparions" }, 
	{ "!", 				GFAT_VECTOR_TYPE,	"Not",					BinArg,		NULL,			"Logical"	 },  
	{ "&&",				GFAT_VECTOR_TYPE,	"And (logical)",		BinArg,		NULL,			"Logical"	 },  
	{ "||",				GFAT_VECTOR_TYPE,	"Or (logical)",			BinArg,		NULL,			"Logical"	 }
};

const size_t BuiltInFuncTableSize = sizeof(BuiltInFunctionsTable)/sizeof(BuiltInFunctionsTable[0]);


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

