/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file argconvdefault.cpp
 *  \brief file for string conversion
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#include "xxxproject/argconvdefault.h"
#include "gpbase/argconv.h"
#include "xxxproject/mktdatas.h"

CC_BEGIN_NAMESPACE( ARM )

ARGConvRevTable MktDataTypeTable[]=
{
	/// methodFlag		methodName
	{	YC,			"YC"		}, 
	{	BSMOD,		"YC_BASIS"	}, 
	{	CAPMOD,		"CAPMOD"	},  
	{	OSWMOD,		"OSWMOD"	},  
	{	SPOTFX,		"FOREX"		},  
	{	FXMOD,		"FXMOD"		}, 
	{	SOMOD,		"SOMOD"		},

	/// do not forget to put this following end of line
	/// used as an equivalent of '\0' to stop a loop
	{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
};

ARGConvRevTable MktVolTypeTable[]=
{
/* 0  */	{	ATM,		"ATM"		}, 
/* 1  */	{	RHO,		"RHO"		}, 
/* 2  */	{	NU,			"NU"		},  
/* 3  */	{	BETA,		"BETA"		}, 
/* 4  */	{	PIV,		"PIV"		}, 
/* 5  */	{	RR,			"RR"		}, 
/* 6  */	{	STR,		"STR"		}, 
/* 7  */	{	VOL,		"VOL"		},
/* 8  */	{	SMILE,		"SMILE"		},
/* 9  */	{	SHIFT,		"SHIFT"		},
/* 10 */	{	Q,			"Q"			},
/* 11 */	{	ADJ,		"ADJ"		},	
/* 12 */	{	CPI,		"CPI"		},
/* 13 */	{	YOY,		"YOY"		},
			{	ENDOFLINE_INT,		ENDOFLINE_CHAR	}
};




const ARM_ArgConv ARM_ArgConv_MktDataType	( MktDataTypeTable, "MktDataType	Table" );
const ARM_ArgConv ARM_ArgConv_MktVolType	( MktVolTypeTable,	"MktVolType		Table" );

const ARM_ArgConvReverse ARM_ArgConvReverse_MktDataType	( MktDataTypeTable, "MktDataType	Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MktVolType	( MktVolTypeTable,	"MktVolType		Table" );

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

