/*
***************************************************************************
** HEADER FILE: cxutils.h
**
** Defines data structures used in the cxutils library.
***************************************************************************
*/

#ifndef _CX_CXUTILS_H
#define _CX_CXUTILS_H

#include "alib.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Bad day conventions used for converting bad days (weekends and holidays)
    into good days (working days where the market is open for business).
   
    Defined as an enumerated data type. Only the first letter is significant
    in the text format.
   
    Designed to be cast as long into ALIB type. */
typedef enum
{
/** No adjustment is made for bad days. */
    CX_BAD_DAY_NONE = 78,                    /* None */
/** Bad days are adjusted to the first good day after the bad day. */
    CX_BAD_DAY_FOLLOW = 70,                  /* Following */
/** Bad days are adjusted to the first good day before the bad day. */
    CX_BAD_DAY_PREVIOUS = 80,                /* Previous */
/** Bad days are adjusted to the first good day after the bad day, unless
    this takes you into a new month. In the latter case, the adjustment is
    then made to the first good day before the bad day. */
    CX_BAD_DAY_MODIFIED = 77                 /* Modified Following */
} CxTBadDayConv;

/** Day count conventions are used to calculate the time elapsed between
    two dates. These are often used for the calculation of accrued interest
    payments and coupon payments.
   
    Designed to be cast as long into equivalent ALIB type. */
typedef enum
{
/** The ISDA ACT/ACT convention. Not to be confused
    with the Bond ACT/ACT convention which is different! */
    CX_ACT_ACT = 1,                          /* ACT/ACT */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{date2-date1}{365}
    \end{equation*} */
    CX_ACT_365F = 2,                         /* ACT/365F */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{date2-date1}{360}
    \end{equation*}
    This is a very common convention for CDS coupon calculation. */
    CX_ACT_360 = 3,                          /* ACT/360 */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{DaysDiff(date1,date2)}{360}
    \end{equation*}
    where the days difference is calculated as
    $360.(year2-year1) + 30.(month2-month1) + (day2 - day1)$
    where we convert day1 from 31 to 30 and we convert day2 from 31 to 30
    but only if day1 is 30 (or 31).
    This is a very common convention for swap coupon calculation. */
    CX_B30_360 = 4,                          /* 30/360 */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{DaysDiff(date1,date2)}{360}
    \end{equation*}
    where the days difference is calculated as
    $360.(year2-year1) + 30.(month2-month1) + (day2 - day1)$
    where we convert day1 from 31 to 30 and we convert day2 from 31 to 30
    in all cases. */
    CX_B30E_360 = 5                          /* 30E/360 */
} CxTDayCountConv;

/** This describes stubs for cash flow lists where the start date and end
    date of the cash flow list are not on cycle. A stub is a period which is
    either shorter or longer than the usual interval. For example, cash flows
    at a 3M interval for a total period of 8M will have a short stub of 2M
    or a long stub of 5M.
   
    Designed to be cast as long into ALIB type. */
typedef enum
{
/** Short front stub.
    Dates are on cycle with the end date and a stub (if required) is
    shorter than the usual period. */
    CX_SHORT_FRONT_STUB = 0,                 /* SF */
/** Short back stub.
    Dates are on cycle with the start date and a stub (if required) is
    shorter than the usual period. */
    CX_SHORT_BACK_STUB = 1,                  /* SB */
/** Long front stub.
    Dates are on cycle with the end date and a stub (if required) is longer
    than the usual period. */
    CX_LONG_FRONT_STUB = 2,                  /* LF */
/** Long back stub.
    Dates are on cycle with the start date and a stub (if required)
    is longer than the usual period. */
    CX_LONG_BACK_STUB = 3                    /* LB */
} CxTStubType;

/** Defines the regularity of the points in a surface.
   
    See the Surface structure for more details. */
typedef enum
{
/** A regular surface is a full and complete grid. One dimensional
    interpolation can take place in both directions. */
    CX_SURFACE_REGULAR_GRID,                 /* Grid */
/** A semi-regular surface which can be represented as an array (in x)
    of arrays (in y). You can do one-dimensional interpolation in the
    y-direction for particular values of x. */
    CX_SURFACE_REGULAR_X,                    /* X */
/** A semi-regular surface which can be represented as an array (in y)
    of arrays (in x). You can do one-dimensional interpolation in the
    x-direction for particular values of y. */
    CX_SURFACE_REGULAR_Y,                    /* Y */
/** An irregular surface could not be represented as an array of arrays.
    Only a two-dimensional algorithm would be possible to interpolate from
    such a surface. */
    CX_SURFACE_IRREGULAR                     /* Irregular */
} CxTSurfaceRegularity;

/** Surface interpolation methods - describes how to seek an interior point
    within the surface. Interpolation methods apply in both directions on
    the surface and can be different. */
typedef enum
{
/** No interpolation in this direction - for example for an irregular grid. */
    CX_SURFACE_INTERP_NONE,                  /* None */
/** Linear interpolation. Extrapolation is supported. */
    CX_SURFACE_INTERP_LINEAR,                /* Linear */
/** Traditional spline interpolation. Extrapolation is not supported. */
    CX_SURFACE_INTERP_SPLINE,                /* Spline */
/** Splice interpolation. This method does linear interpolation on parabolas.
    You also need some assumption on how the first and last segment of the
    curve behaves - we interpolate between the parabola and the linearly
    extrapolated line on the edge segments. Extrapolation is not supported. */
    CX_SURFACE_INTERP_SPLICE,                /* Splice */
/** Base correlation interpolation. This is spline with linear extrapolation.
    In addition if the first x-value is negative (x-value is strike for base
    correlation), then the interpolation switches to linear between $x_1$ and
    $x_2$ using the values $z_0$ and $z_2$ instead of $z_1$ and $z_2$. */
    CX_SURFACE_INTERP_BASE_CORR              /* BaseCorr */
} CxTSurfaceInterpMethod;

/** Surface extrapolation methods - describes how to seek a point which is
    not within the bounds of points on the surface. Extrapolation methods
    apply in both directions on the surface and can be different. */
typedef enum
{
/** No extrapolation. */
    CX_SURFACE_EXTRAP_NONE,                  /* N */
/** Extrapolates to the last point on the surface. */
    CX_SURFACE_EXTRAP_FLAT,                  /* F */
/** Extrapolates by continuing the interpolation to points outside the
    surface. */
    CX_SURFACE_EXTRAP_EXTRAP                 /* X */
} CxTSurfaceExtrapMethod;

/** Define the calendar structure. The constructor is clever - it removes
    weekends from the list of dates.
   
    Designed for compatibility with ALIB holiday file functions - the data
    type is the same as one of the low-level data representation in the ALIB.
   
    Note that an undefined calendar = weekends only */
typedef struct _CxTCalendar
{
    TDateList*      dateList;
    long            weekends;
} CxTCalendar;

/** Defines a generic surface object. The surface can have degrees of
    regularity. The values are defined by linked arrays of x-values,
    y-values and z-values where z = f(x,y) where f is a function on the
    surface.
   
    The surface might be regular, semi-regular or irregular.
   
    A regular surface is a full and complete grid. One dimensional
    interpolation can take place in both directions.
   
    A semi-regular surface could realistically be represented as an array
    of arrays. One-dimensional interpolation can take place in one
    direction.
   
    An irregular surface could not be represented as an array of arrays.
    Only a two-dimensional algorithm would be possible to interpolate from
    such a surface.
   
    Interpolation rules will be defined separately. */
typedef struct _CxTSurface
{
    /** Defines array size for x, y, z. */
    int             numItems;
    /** Array of size numItems.
        x-values for the surface function z=f(x,y) */
    double*         x;
    /** Array of size numItems.
        y-values for the surface function z=f(x,y) */
    double*         y;
    /** Array of size numItems.
        z-values for the surface function z=f(x,y) */
    double*         z;
    /** Regularity of the x,y points. Used to determine what types of interpolation
        are possible. */
    CxTSurfaceRegularity regularity;
} CxTSurface;

/** Overall surface interpolation structure - defines interpolation and
    extrapolation in both directions, or whether to use surface interpolation
    instead of two one-dimensional interpolations. */
typedef struct _CxTSurfaceInterp
{
    /** Use a 2-dimensional surface interpolation. Not actually implemented. */
    TBoolean        interpSurface;
    /** Interpolation method in the x-direction. */
    CxTSurfaceInterpMethod xInterp;
    /** Extrapolation method in the x-direction. */
    CxTSurfaceExtrapMethod xExtrap;
    /** Interpolation method in the y-direction. */
    CxTSurfaceInterpMethod yInterp;
    /** Extrapolation method in the y-direction. */
    CxTSurfaceExtrapMethod yExtrap;
} CxTSurfaceInterp;

/** Defines the current state of a multi-dimensional minimization. Designed so
    that when this structure is returned from the multi-dimensional solver
    that it can be used as the initial state for a subsequent call. */
typedef struct _CxTMultiDimMinState
{
    /** Defines array size for x.
        Number of dimensions */
    int             n;
    /** Array of size n.
        Initial starting point (on initialization).
        Solution (on exit). */
    double*         x;
    /** Initial directions for search (on initialization). Final directions for
        search (on exit). */
    TMatrix2D*      direction;
    /** Number of iterations. */
    int             iter;
    /** Minimal value. */
    double          value;
    /** Last step-size in the algorithm. */
    double          vdiff;
} CxTMultiDimMinState;

/**
***************************************************************************
** Converts BadDayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxBadDayConvToString (CxTBadDayConv value);

/**
***************************************************************************
** Converts BadDayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxBadDayConvFromString (const char* str, CxTBadDayConv *val);

/**
***************************************************************************
** Converts DayCountConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxDayCountConvToString (CxTDayCountConv value);

/**
***************************************************************************
** Converts DayCountConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxDayCountConvFromString (const char* str, CxTDayCountConv *val);

/**
***************************************************************************
** Converts StubType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxStubTypeToString (CxTStubType value);

/**
***************************************************************************
** Converts StubType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxStubTypeFromString (const char* str, CxTStubType *val);

/**
***************************************************************************
** Converts SurfaceRegularity to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxSurfaceRegularityToString (CxTSurfaceRegularity value);

/**
***************************************************************************
** Converts SurfaceRegularity from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceRegularityFromString (const char* str, CxTSurfaceRegularity *val);

/**
***************************************************************************
** Converts SurfaceInterpMethod to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxSurfaceInterpMethodToString (CxTSurfaceInterpMethod value);

/**
***************************************************************************
** Converts SurfaceInterpMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceInterpMethodFromString (const char* str, CxTSurfaceInterpMethod *val);

/**
***************************************************************************
** Converts SurfaceExtrapMethod to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxSurfaceExtrapMethodToString (CxTSurfaceExtrapMethod value);

/**
***************************************************************************
** Converts SurfaceExtrapMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceExtrapMethodFromString (const char* str, CxTSurfaceExtrapMethod *val);

/**
***************************************************************************
** Constructor for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarMake(
TDateList*      dateList,
long            weekends
);

/**
***************************************************************************
** Memory allocator for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarCopy(CxTCalendar* src);

/**
***************************************************************************
** Destructor for CxTCalendar
***************************************************************************
*/
void CxCalendarFree(CxTCalendar *p);

/**
***************************************************************************
** Constructor for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMake(
/** Defines array size for x, y, z. */
int             numItems,
/** Array of size numItems.
    x-values for the surface function z=f(x,y) */
double*         x,
/** Array of size numItems.
    y-values for the surface function z=f(x,y) */
double*         y,
/** Array of size numItems.
    z-values for the surface function z=f(x,y) */
double*         z,
/** Regularity of the x,y points. Used to determine what types of interpolation
    are possible. */
CxTSurfaceRegularity regularity
);

/**
***************************************************************************
** Memory allocator for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMakeEmpty(
/** Defines array size for x, y, z. */
int             numItems
);

/**
***************************************************************************
** Copy constructor for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceCopy(CxTSurface* src);

/**
***************************************************************************
** Destructor for CxTSurface
***************************************************************************
*/
void CxSurfaceFree(CxTSurface *p);

/**
***************************************************************************
** Constructor for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpMake(
/** Use a 2-dimensional surface interpolation. Not actually implemented. */
TBoolean        interpSurface,
/** Interpolation method in the x-direction. */
CxTSurfaceInterpMethod xInterp,
/** Extrapolation method in the x-direction. */
CxTSurfaceExtrapMethod xExtrap,
/** Interpolation method in the y-direction. */
CxTSurfaceInterpMethod yInterp,
/** Extrapolation method in the y-direction. */
CxTSurfaceExtrapMethod yExtrap
);

/**
***************************************************************************
** Memory allocator for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpCopy(CxTSurfaceInterp* src);

/**
***************************************************************************
** Destructor for CxTSurfaceInterp
***************************************************************************
*/
void CxSurfaceInterpFree(CxTSurfaceInterp *p);

/**
***************************************************************************
** Constructor for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateMake(
/** Defines array size for x.
    Number of dimensions */
int             n,
/** Array of size n.
    Initial starting point (on initialization).
    Solution (on exit). */
double*         x,
/** Initial directions for search (on initialization). Final directions for
    search (on exit). */
TMatrix2D*      direction,
/** Number of iterations. */
int             iter,
/** Minimal value. */
double          value,
/** Last step-size in the algorithm. */
double          vdiff
);

/**
***************************************************************************
** Memory allocator for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateMakeEmpty(
/** Defines array size for x.
    Number of dimensions */
int             n
);

/**
***************************************************************************
** Copy constructor for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateCopy(CxTMultiDimMinState* src);

/**
***************************************************************************
** Destructor for CxTMultiDimMinState
***************************************************************************
*/
void CxMultiDimMinStateFree(CxTMultiDimMinState *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _CX_CXUTILS_H */
