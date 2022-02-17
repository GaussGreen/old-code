/*
***************************************************************************
** HEADER FILE: surface.h
**
** Generic surface functions.
**
** $Header$
***************************************************************************
*/

#ifndef CX_SURFACE_H
#define CX_SURFACE_H

#include "cxutils.h"  /* basic data types */

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Constructor for CXSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMakeFromInputs(
    int             numItems,            /* (I) */
    double*         x,                   /* (I) [numItems] */
    double*         y,                   /* (I) [numItems] */
    double*         z,                   /* (I) [numItems] */
    TBoolean*       include,             /* (I) [numItems] optional */
    TBoolean        xRegular,
    TBoolean        yRegular);

/*
***************************************************************************
** Interpolate a single value from a surface
***************************************************************************
*/
int CxSurfaceInterp(
    CxTSurface*           surface,
    CxTSurfaceInterp*     interp,
    double               x,
    double               y,
    double              *result);

/*
***************************************************************************
** Interpolate a matrix from a surface
***************************************************************************
*/
TMatrix2D* CxSurfaceInterpMatrix(
    CxTSurface*           surface,
    CxTSurfaceInterp*     interp,
    int                  xSize,
    double              *x,
    int                  ySize,
    double              *y);

/*
***************************************************************************
** Interpolate a data table from a surface
***************************************************************************
*/
TDataTable* CxSurfaceInterpTable(
    CxTSurface*           surface,
    CxTSurfaceInterp*     interp,
    int                  xSize,
    double              *x,
    char               **xHeading,
    int                  ySize,
    double              *y,
    char               **yHeading);

#ifdef __cplusplus
}
#endif

#endif
