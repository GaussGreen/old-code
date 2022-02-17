/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file surfacetypedef.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date Ocotber 2004
 */

#ifndef _INGPBASE_SURFACETYPEDEF_H
#define _INGPBASE_SURFACETYPEDEF_H

/// this header comes firts as it includes some preprocessor constants!
#include "removeidentifiedwarning.h"

/// header for portable macro
#include "port.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >		class ARM_T_SurfaceBase;
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >		class ARM_T_FlatSurface;
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >		class ARM_T_SurfaceWithInterpol;
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >		class ARM_T_Interporlator2DLin2CstExtrapol;
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >		class ARM_T_Interporlator2DLin1CstExtrapol;
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >    	class ARM_T_Interporlator2D;


/// surface definition
typedef ARM_T_SurfaceBase<double,double,double>								ARM_Surface;
typedef ARM_T_FlatSurface<double,double,double>								ARM_FlatSurface;
typedef ARM_T_SurfaceWithInterpol<double,double,double>						ARM_SurfaceWithInterpol;
typedef ARM_T_Interporlator2D<double,double,double>		                    ARM_Interporlator2D;
typedef ARM_T_Interporlator2DLin1CstExtrapol<double,double,double>		    ARM_2DLin1Interpol;
typedef ARM_T_Interporlator2DLin2CstExtrapol<double,double,double>		    ARM_2DLin2Interpol;



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
