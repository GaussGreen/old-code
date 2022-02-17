/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM for the gpinflation project
 *	\author  EBenhamou
 *	\version 1.0
 *	\date September 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_TYPEDEF_H
#define _INGPINFLATION_TYPEDEF_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpbase/gpvector.h"
#include "glob/linalg.h"

#include "gpbase/typedef.h"


/// Start ARM namespace
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration in global namespace (ARM kernel)
class ARM_InfBSModel;

/// reference counted pointor
typedef ARM_CountedPtr< ARM_InfBSModel >  ARM_InfBSModelPtr;

#define K_NAIVE     0
#define K_GAUSS     1
#define K_TEST		2 

typedef enum { MAI, SUB, SUP }	IndexType;

typedef enum { EMU, IFRF }	InfIndexName;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
