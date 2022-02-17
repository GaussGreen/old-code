/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvddefault.cpp,v $
 * Revision 1.1  2004/09/09 16:39:43  ebenhamou
 * Initial revision
 *
 */


/*! \file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#include "gpinflation/argconvdefault.h"
#include "gpinflation/infswopfactory.h"

CC_BEGIN_NAMESPACE( ARM )

ARGConvTable InfSwoptComputationMethodTable[] =
{	
    /// Type Name       /// number
    "Std",				ARM_InfVolComputation_Factory::INF_SWOPT_STD,
    "Equal Weight",     ARM_InfVolComputation_Factory::INF_SWOPT_EQUAL_WEIGHT,
    "DF Weight",		ARM_InfVolComputation_Factory::INF_SWOPT_DF_WEIGHT,
    "DF Weight Square",	ARM_InfVolComputation_Factory::INF_SWOPT_DF_WEIGHT_SQUARE,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_InfSwoptComputationMethod( InfSwoptComputationMethodTable, "InfSwoptComputationMethod" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InfSwoptComputationMethod( InfSwoptComputationMethodTable, "InfSwoptComputationMethod" );


ARGConvRevTable CopulaTable[]=
	{
		/// methodFlag		methodName
		{	K_NAIVE,		"NAIVE"	}, 
		{	K_GAUSS,		"GAUSS"	}, 

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_Copula( CopulaTable, "Copula Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_Copula( CopulaTable, "Copula Table" );

ARGConvRevTable InfIndexTable[]=
	{
		/// methodFlag		methodName
		{	EMU,		"EMU"	}, 
		{	IFRF,		"IFRF"	}, 

		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_InfIndex( InfIndexTable, "Inf Index Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InfIndex( InfIndexTable, "Inf Index Table" );

ARGConvRevTable SubordIndexTable[]=
	{
		/// methodFlag		methodName
		{	MAI,		"Main"	}, 
		{	SUB,		"Sub"	}, 
		{	SUP,		"Sup"	},
		/// do not forget to put this following end of line
		/// used as an equivalent of '\0' to stop a loop
		{	ENDOFLINE_INT	, ENDOFLINE_CHAR	}
	};

const ARM_ArgConv ARM_ArgConv_SubordIndex( SubordIndexTable, "Subord Index Table" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SubordIndex( SubordIndexTable, "Subord Index Table" );



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

