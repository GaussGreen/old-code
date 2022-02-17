/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file argconvdefault.h
 *  \brief file for string conversion
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */


#ifndef _INXXXPROJECT_ARGCONVDEFAULT_H
#define _INXXXPROJECT_ARGCONVDEFAULT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include "gpbase/argconv.h"

typedef enum { 	ATM,
				RHO,
				NU,  
				BETA,
				PIV,
				RR,
				STR,
				VOL,
				SMILE,
				SHIFT,
				Q,
				ADJ,
				CPI,
				YOY		}	MktVolType;


CC_BEGIN_NAMESPACE( ARM )


typedef enum { 	YC,		
				BSMOD,		
				CAPMOD,		
				OSWMOD,		
				SPOTFX,		
				FXMOD,	
				FXMIX,		
				SOMOD	}	MktModType;


extern const ARM_ArgConv				ARM_ArgConv_MktDataType;
extern const ARM_ArgConvReverse			ARM_ArgConvReverse_MktDataType;

extern const ARM_ArgConv				ARM_ArgConv_MktVolType;
extern const ARM_ArgConvReverse			ARM_ArgConvReverse_MktVolType;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
