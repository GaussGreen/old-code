/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM
 *	\author  J-M Prié, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPNUMMETHODS_TYPEDEF_H
#define _INGPNUMMETHODS_TYPEDEF_H

#include "gpinfra/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_SliceBase;
class ARM_GridIndex;
class ARM_TreeTransition;
struct ARM_ReconnectorBase;
struct ARM_TruncatorBase;
class ARM_ExerciseBoundaryCalc;
class ARM_ImpSampler;
class ARM_PathScheme;

/// vector
typedef vector< ARM_SliceBase* >                        ARM_SliceVector;
typedef vector< ARM_TreeTransition* >                   ARM_TreeTransitions;
typedef vector< ARM_TreeTransitions* >                  ARM_TreeTransitionsVector;
typedef vector< ARM_GridIndex* >						ARM_GridIndexVector;

													
/// reference counted pointor
typedef ARM_CountedPtr< ARM_SliceVector >			    ARM_SliceVectorPtr;
typedef ARM_CountedPtr< ARM_ExerciseBoundaryCalc >      ARM_ExerciseBoundaryCalcPtr;
typedef ARM_CountedPtr< ARM_ImpSampler >				ARM_ImpSamplerPtr;
typedef ARM_CountedPtr< ARM_PathScheme >				ARM_PathSchemePtr;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
