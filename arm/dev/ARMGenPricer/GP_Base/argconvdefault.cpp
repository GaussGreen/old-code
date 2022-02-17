/*!
 *
 * Copyright (c) IXIS CIB July 2003 Paris
 *
 *	\file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#include "gpbase/argconvdefault.h"
#include "gpbase/argconv.h"
#include "gpbase/enumbase.h"

/// ARM Kernel
#include <glob/armdef.h>

CC_BEGIN_NAMESPACE( ARM ) 	

ARGConvRevTable LgNameDayCountTable[]=
	{
		/// methodFlag		methodName
		{	KACTUAL_ACTUAL,		"ACTUAL ACTUAL"	}, 
		{	KACTUAL_365,		"A365"	},  
		{	KACTUAL_360,		"A360"	},  
		{	K30_360,			"30/360"		},  
		{	KACTUAL_REAL,		"ACTUAL REAL"	},  
		{	KACTUAL_FEB29,		"ACTUAL FEB29"	},  
		{	KACTUAL_ISMA,		"ACTUAL ISMA"	},
        {   K30_360E,           "30E"           },
		{	KNOBASE,			"NO BASE"		},

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_LgNameDayCount( LgNameDayCountTable, "LgNameDayCount Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_LgNameDayCount( LgNameDayCountTable, "LgNameDayCount Table" );

ARGConvRevTable LongNameFrequencyTable[]=
{
		/// methodFlag		methodName
		{	K_ANNUAL,		"ANNUAL"		}, 
		{	K_SEMIANNUAL,	"SEMIANNUAL"	}, 
		{	K_QUARTERLY,	"QUARTERLY"		}, 
		{	K_BIMONTHLY,	"BIMONTHLY"		}, 
		{	K_MONTHLY,		"MONTHLY"		}, 
		{	K_WEEKLY,		"WEEKLY"		}, 
		{	K_DAILY,		"DAILY"			}, 
		{	K_ZEROCOUPON,	"ZEROCOUPON"	}, 
		{	K_ANNUAL,		"A"				}, 
		{	K_SEMIANNUAL,	"S"				}, 
		{	K_QUARTERLY,	"Q"				}, 
		{	K_MONTHLY,		"M"				}, 
		{	K_WEEKLY,		"W"				}, 
		{	K_DAILY,		"D"				}, 
		{	K_DEF_FREQ,		"DEF_FREQ, default value used differently in various contexts "	},
		
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
};

ARGConvRevTable LongNameFrequencyTableReverse[]=
{
		/// methodFlag		methodName
		{	K_ANNUAL,		"ANNUAL"	}, 
		{	K_SEMIANNUAL,	"SEMIANNUAL"	}, 
		{	K_QUARTERLY,	"QUARTERLY"	}, 
		{	K_BIMONTHLY,	"BIMONTHLY"	}, 
		{	K_MONTHLY,		"MONTHLY"	}, 
		{	K_WEEKLY,		"WEEKLY"	}, 
		{	K_DAILY,		"DAILY"	}, 
		{	K_ZEROCOUPON,	"ZEROCOUPON"	}, 
		{	K_DEF_FREQ,		"DEF_FREQ, default value used differently in various contexts "	},
		
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
};

const ARM_ArgConv ARM_ArgConv_LgNameFrequency( LongNameFrequencyTable, "Lg Name Frequency" );
const ARM_ArgConvReverse ARM_ArgConvReverse_LgNameFrequency( LongNameFrequencyTableReverse, "Lg Name Frequency" );


ARGConvRevTable FwdRulesTable[]=
	{
		/// methodFlag				methodName
		{	K_PREVIOUS,				"PREVIOUS"	}, 
		{	K_MOD_PREVIOUS,			"MOD_PREVIOUS"	}, 
		{	K_FOLLOWING,			"FOLLOWING"	}, 
		{	K_MOD_FOLLOWING,		"MOD_FOLLOWING"	}, 
		
		/// do not forget to put this following end of line
        /// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_FwdRules( FwdRulesTable, "Fwd Rules" );
const ARM_ArgConvReverse ARM_ArgConvReverse_FwdRules( FwdRulesTable, "Fwd Rules" );

ARGConvRevTable CompoundingTypeTable[]=
	{
		/// methodFlag				methodName
		{	K_COMP_NONE,			"NONE"	}, 
		{	K_SPREAD_INC,			"SPREADINC"	}, 
		{	K_SPREAD_EXC,			"SPREADEXC"	}, 
		{	K_FLAT,					"FLAT"		}, 
		
		/// do not forget to put this following end of line
        /// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_CompoundingType( CompoundingTypeTable, "Compounding Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CompoundingType( CompoundingTypeTable, "Compounding Type" );

ARGConvRevTable CompoundingFrequencyTable[]=
{
		/// methodFlag		methodName
		{	K_COMP_CONT,		"CONT"				}, 
		{	K_COMP_PROP,		"NONE"				}, 
		{	K_COMP_ANNUAL,		"ANNUAL"			}, 
		{	K_COMP_SEMIANNUAL,	"SEMIANNUAL"		}, 
		{	K_COMP_QUARTERLY,	"QUARTERLY"			}, 
		{	K_COMP_MONTHLY,		"MONTHLY"			}, 
		{	K_COMP_BIMONTHLY,	"BIMONTHLY"			}, 
		{	K_COMP_DAILY_360,	"DAILY360"			}, 
		{	K_COMP_DAILY_365,	"DAILY365"			}, 
		{	K_COMP_ANNUAL,		"A"					}, 
		{	K_COMP_SEMIANNUAL,	"S"					}, 
		{	K_COMP_QUARTERLY,	"Q"					}, 
		{	K_COMP_MONTHLY,		"M"					}, 
		{	K_COMP_DAILY_365,	"D"					}, 
		
		{	K_DEF_FREQ,		"DEF_FREQ, default value used differently in various contexts "	},
		
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
};

const ARM_ArgConv ARM_ArgConv_CompoundingFrequency( CompoundingFrequencyTable, "Compounding Frequency" );
//const ARM_ArgConvReverse ARM_ArgConvReverse_CompoundingFrequency( CompoundingFrequencyTable, "Compounding Frequency" );

ARGConvRevTable LgTimingModTable[]=
	{
		/// methodFlag		methodName
		{	K_ADVANCE,		"ADVANCE"	}, 
		{	K_ARREARS,		"ARREARS"	}, 
		
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_LgTimingMod( LgTimingModTable, "LgTimingModTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_LgTimingMod( LgTimingModTable, "LgTimingModTable" );


ARGConvRevTable InterestRulesTable[]=
	{
		/// methodFlag				methodName
		{	K_ADJUSTED,				"ADJUSTED"	}, 
		{	K_UNADJUSTED,			"UNADJUSTED"	}, 
		
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_InterestRules( InterestRulesTable, "InterestRules Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InterestRules( InterestRulesTable, "InterestRules Table" );


ARGConvRevTable StubRulesTable[]=
	{
		/// methodFlag				methodName
		{	K_SHORTSTART,			"SHORTSTART"}, 
		{	K_LONGSTART,			"LONGSTART"	}, 
		{	K_SHORTEND,				"SHORTEND"	}, 
		{	K_LONGEND,				"LONGEND"	}, 
		{	K_SHORTSTART,			"SS"		}, 
		{	K_LONGSTART,			"LS"		}, 
		{	K_SHORTEND,				"SE"		}, 
		{	K_LONGEND,				"LE"		}, 

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};


ARGConvRevTable StubRulesTableReverse[]=
	{
		/// methodFlag				methodName
		{	K_SHORTSTART,			"SHORTSTART"}, 
		{	K_LONGSTART,			"LONGSTART"	}, 
		{	K_SHORTEND,				"SHORTEND"	}, 
		{	K_LONGEND,				"LONGEND"	}, 

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_StubRules( StubRulesTable, "StubRules Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_StubRules( StubRulesTableReverse, "StubRules Table" );

ARGConvRevTable ReceiveOrPayTable[]=
	{
		/// methodFlag		methodName
		{	K_RCV,		"RCV"	}, 
		{	K_PAY,		"PAY"	}, 

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_RcvOrPay( ReceiveOrPayTable, "Receive or Pay Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_RcvOrPay( ReceiveOrPayTable, "Receive or Pay Table" );

ARGConvTable InterpolationTypeTable[]=
	{
		/// methodFlag		methodName
		{"LINEAR_COLUMN",				ARM_InterpolationType::linear_column_extrapoleCst					}, 
		{"LINEAR_ROW",					ARM_InterpolationType::linear_row_extrapoleCst						},  
		{"LINEAR_COLUMN_ROW",			ARM_InterpolationType::linear_column_row_extrapoleCst_column_row	},  
		{"LINEAR_ROW_COLUMN",			ARM_InterpolationType::linear_row_column_extrapoleCst_row_column	},  
		{"STEPUP_RIGHT_COLUMN",			ARM_InterpolationType::stepup_right_column							},  
		{"STEPUP_RIGHT_ROW",			ARM_InterpolationType::stepup_right_row								},  
		{"STEPUP_LEFT_COLUMN",			ARM_InterpolationType::stepup_left_column							},
		{"STEPUP_LEFT_ROW",				ARM_InterpolationType::stepup_left_row								},
		{"CONSTANT",					ARM_InterpolationType::unknown										},

		/// very important as it tells that this is the end
		ENDOFLINE_CHAR	

	};

const ARM_ArgConv ARM_ArgConv_InterpolType( InterpolationTypeTable, "InterpolationType Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InterpolType( InterpolationTypeTable, "InterpolationType Table" );

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

