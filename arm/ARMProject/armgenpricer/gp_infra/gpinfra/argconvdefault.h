/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.h
 *
 *  \brief default table for interface reading
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFRA_ARGCONVDEFAUIT_H
#define _INGPINFRA_ARGCONVDEFAUIT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include <gpbase/argconv.h>

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_DayCount;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_DayCount;

extern const ARM_ArgConv ARM_ArgConv_CalibMethod;

extern const ARM_ArgConv ARM_ArgConv_TargetFuncMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TargetFuncMethod;

extern const ARM_ArgConv ARM_ArgConv_CalibDirectionMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CalibDirectionMethod;

extern const ARM_ArgConv ARM_ArgConv_CalibDirection2DMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CalibDirection2DMethod;

extern const ARM_ArgConv ARM_ArgConv_ModelParam;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ModelParam;

extern const ARM_ArgConv ARM_ArgConv_ModelParamDataType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ModelParamDataType;

extern const ARM_ArgConv ARM_ArgConv_Numeraire;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_Numeraire;

extern const ARM_ArgConv ARM_ArgConv_CapFloor;

extern const ARM_ArgConv ARM_ArgConv_CallPut;

extern const ARM_ArgConv ARM_ArgConv_InOut;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InOut;

extern const ARM_ArgConv ARM_ArgConv_PayRec;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PayRec;

extern const ARM_ArgConv ARM_ArgConv_StdFrequency;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_StdFrequency;

extern const ARM_ArgConv ARM_ArgConv_MatFrequency;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MatFrequency;

extern const ARM_ArgConv ARM_ArgConv_IndexType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_IndexType;

extern const ARM_ArgConv ARM_ArgConv_Timing;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_Timing;

extern const ARM_ArgConv ARM_ArgConv_VolType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_VolType;

extern const ARM_ArgConv ARM_ArgConv_VolMktType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_VolMktType;

extern const ARM_ArgConv ARM_ArgConv_interpolCurveType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_interpolCurveType;

extern const ARM_ArgConv ARM_ArgConv_cfDistType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_cfDistType;

extern const ARM_ArgConv ARM_ArgConv_cfMepiDistType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_cfMepiDistType;

extern const ARM_ArgConv ARM_ArgConv_cfGreekType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_cfGreekType;

extern const ARM_ArgConv ARM_ArgConv_IndexClass;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_IndexClass;

extern const ARM_ArgConv ARM_ArgConv_IntRule;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_IntRule;

extern const ARM_ArgConv ARM_ArgConv_YesNo;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_YesNo;

extern const ARM_ArgConv ARM_ArgConv_BumpParam;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_BumpParam;

extern const ARM_ArgConv ARM_ArgConv_InterpolInfType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InterpolInfType;

extern const ARM_ArgConv ARM_ArgConv_FwdRule;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_FwdRule;

extern const ARM_ArgConv ARM_ArgConv_DigitType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_DigitType;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
