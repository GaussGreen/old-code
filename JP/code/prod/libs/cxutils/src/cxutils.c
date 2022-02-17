/*
***************************************************************************
** SOURCE FILE: cxutils.c
**
** Defines data structures used in the cxutils library.
***************************************************************************
*/

#include "cxutils.h"
#include "cxmacros.h"
#include "calendar.h"
#include "zerocurve.h"
#include "dateutils.h"
#include "surface.h"
#include "mdmin.h"

#include <ctype.h>

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
char* CxBadDayConvToString (CxTBadDayConv value)
{
    static char routine[] = "CxBadDayConvToString";

    switch (value)
    {
    case CX_BAD_DAY_NONE:               return "None";
    case CX_BAD_DAY_FOLLOW:             return "Following";
    case CX_BAD_DAY_PREVIOUS:           return "Previous";
    case CX_BAD_DAY_MODIFIED:           return "Modified Following";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts BadDayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxBadDayConvFromString (const char* str, CxTBadDayConv *val)
{
    static char routine[] = "CxBadDayConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "F") == 0) *val = CX_BAD_DAY_FOLLOW;
        else goto done;
        break;
    case 'M':
        if (strcmp (buf, "M") == 0) *val = CX_BAD_DAY_MODIFIED;
        else goto done;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = CX_BAD_DAY_NONE;
        else goto done;
        break;
    case 'P':
        if (strcmp (buf, "P") == 0) *val = CX_BAD_DAY_PREVIOUS;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* CxDayCountConvToString (CxTDayCountConv value)
{
    static char routine[] = "CxDayCountConvToString";

    switch (value)
    {
    case CX_ACT_ACT:                    return "ACT/ACT";
    case CX_ACT_365F:                   return "ACT/365F";
    case CX_ACT_360:                    return "ACT/360";
    case CX_B30_360:                    return "30/360";
    case CX_B30E_360:                   return "30E/360";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts DayCountConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxDayCountConvFromString (const char* str, CxTDayCountConv *val)
{
    static char routine[] = "CxDayCountConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[9];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case '3':
        if (strcmp (buf, "30/360") == 0) *val = CX_B30_360;
        else if (strcmp (buf, "30E/360") == 0) *val = CX_B30E_360;
        else goto done;
        break;
    case 'A':
        if (strcmp (buf, "ACT/360") == 0) *val = CX_ACT_360;
        else if (strcmp (buf, "ACT/365F") == 0) *val = CX_ACT_365F;
        else if (strcmp (buf, "ACT/ACT") == 0) *val = CX_ACT_ACT;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* CxStubTypeToString (CxTStubType value)
{
    static char routine[] = "CxStubTypeToString";

    switch (value)
    {
    case CX_SHORT_FRONT_STUB:           return "SF";
    case CX_SHORT_BACK_STUB:            return "SB";
    case CX_LONG_FRONT_STUB:            return "LF";
    case CX_LONG_BACK_STUB:             return "LB";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts StubType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxStubTypeFromString (const char* str, CxTStubType *val)
{
    static char routine[] = "CxStubTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[3];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'L':
        if (strcmp (buf, "LB") == 0) *val = CX_LONG_BACK_STUB;
        else if (strcmp (buf, "LF") == 0) *val = CX_LONG_FRONT_STUB;
        else goto done;
        break;
    case 'S':
        if (strcmp (buf, "SB") == 0) *val = CX_SHORT_BACK_STUB;
        else if (strcmp (buf, "SF") == 0) *val = CX_SHORT_FRONT_STUB;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* CxSurfaceRegularityToString (CxTSurfaceRegularity value)
{
    static char routine[] = "CxSurfaceRegularityToString";

    switch (value)
    {
    case CX_SURFACE_REGULAR_GRID:       return "Grid";
    case CX_SURFACE_REGULAR_X:          return "X";
    case CX_SURFACE_REGULAR_Y:          return "Y";
    case CX_SURFACE_IRREGULAR:          return "Irregular";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts SurfaceRegularity from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceRegularityFromString (const char* str, CxTSurfaceRegularity *val)
{
    static char routine[] = "CxSurfaceRegularityFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'G':
        if (strcmp (buf, "G") == 0) *val = CX_SURFACE_REGULAR_GRID;
        else goto done;
        break;
    case 'I':
        if (strcmp (buf, "I") == 0) *val = CX_SURFACE_IRREGULAR;
        else goto done;
        break;
    case 'X':
        if (strcmp (buf, "X") == 0) *val = CX_SURFACE_REGULAR_X;
        else goto done;
        break;
    case 'Y':
        if (strcmp (buf, "Y") == 0) *val = CX_SURFACE_REGULAR_Y;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* CxSurfaceInterpMethodToString (CxTSurfaceInterpMethod value)
{
    static char routine[] = "CxSurfaceInterpMethodToString";

    switch (value)
    {
    case CX_SURFACE_INTERP_NONE:        return "None";
    case CX_SURFACE_INTERP_LINEAR:      return "Linear";
    case CX_SURFACE_INTERP_SPLINE:      return "Spline";
    case CX_SURFACE_INTERP_SPLICE:      return "Splice";
    case CX_SURFACE_INTERP_BASE_CORR:   return "BaseCorr";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts SurfaceInterpMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceInterpMethodFromString (const char* str, CxTSurfaceInterpMethod *val)
{
    static char routine[] = "CxSurfaceInterpMethodFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[9];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTSurfaceInterpMethod)(CX_SURFACE_INTERP_NONE);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'B':
        if (strcmp (buf, "BASECORR") == 0) *val = CX_SURFACE_INTERP_BASE_CORR;
        else goto done;
        break;
    case 'L':
        if (strcmp (buf, "LINEAR") == 0) *val = CX_SURFACE_INTERP_LINEAR;
        else goto done;
        break;
    case 'N':
        if (strcmp (buf, "NONE") == 0) *val = CX_SURFACE_INTERP_NONE;
        else goto done;
        break;
    case 'S':
        if (strcmp (buf, "SPLICE") == 0) *val = CX_SURFACE_INTERP_SPLICE;
        else if (strcmp (buf, "SPLINE") == 0) *val = CX_SURFACE_INTERP_SPLINE;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* CxSurfaceExtrapMethodToString (CxTSurfaceExtrapMethod value)
{
    static char routine[] = "CxSurfaceExtrapMethodToString";

    switch (value)
    {
    case CX_SURFACE_EXTRAP_NONE:        return "N";
    case CX_SURFACE_EXTRAP_FLAT:        return "F";
    case CX_SURFACE_EXTRAP_EXTRAP:      return "X";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts SurfaceExtrapMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxSurfaceExtrapMethodFromString (const char* str, CxTSurfaceExtrapMethod *val)
{
    static char routine[] = "CxSurfaceExtrapMethodFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTSurfaceExtrapMethod)(CX_SURFACE_EXTRAP_NONE);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "F") == 0) *val = CX_SURFACE_EXTRAP_FLAT;
        else goto done;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = CX_SURFACE_EXTRAP_NONE;
        else goto done;
        break;
    case 'X':
        if (strcmp (buf, "X") == 0) *val = CX_SURFACE_EXTRAP_EXTRAP;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Constructor for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarMake(
TDateList*      dateList,            /* (I) */
long            weekends             /* (I) */
)
{
    static char routine[] = "CxCalendarMake";
    int status = FAILURE;

    CxTCalendar* p = NULL;


    p = NEW(CxTCalendar);
    if (p==NULL) goto done;

    if (dateList != NULL)
    {
        p->dateList = GtoCopyDateList(dateList);
        if (p->dateList == NULL) goto done;
    }
    else
    {
        p->dateList = NULL;
    }

    p->weekends        = weekends;

    if (CxCalendarValidate(p) != SUCCESS) goto done; /* failure */

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCalendarFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarMakeEmpty(void)
{
    static char routine[] = "CxCalendarMakeEmpty";
    int status = FAILURE;

    CxTCalendar* p = NULL;


    p = NEW(CxTCalendar);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCalendarFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCalendar
***************************************************************************
*/
CxTCalendar* CxCalendarCopy(CxTCalendar* src)
{
    CxTCalendar* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCalendarMake(src->dateList,
                         src->weekends);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCalendar
***************************************************************************
*/
void CxCalendarFree(CxTCalendar *p)
{
    if (p != NULL)
    {
        GtoFreeDateList(p->dateList);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMake(
int             numItems,            /* (I) */
double*         x,                   /* (I) [numItems] */
double*         y,                   /* (I) [numItems] */
double*         z,                   /* (I) [numItems] */
CxTSurfaceRegularity regularity           /* (I) */
)
{
    static char routine[] = "CxSurfaceMake";
    int status = FAILURE;

    CxTSurface* p = NULL;

    REQUIRE(numItems > 0);
    REQUIRE(x != NULL);
    REQUIRE(y != NULL);
    REQUIRE(z != NULL);

    p = NEW(CxTSurface);
    if (p==NULL) goto done;

    p->numItems        = numItems;
    p->x = NEW_ARRAY(double, p->numItems);
    if (p->x == NULL) goto done;
    COPY_ARRAY (p->x, x, double, p->numItems);

    p->y = NEW_ARRAY(double, p->numItems);
    if (p->y == NULL) goto done;
    COPY_ARRAY (p->y, y, double, p->numItems);

    p->z = NEW_ARRAY(double, p->numItems);
    if (p->z == NULL) goto done;
    COPY_ARRAY (p->z, z, double, p->numItems);

    p->regularity      = regularity;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxSurfaceFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMakeEmpty(
int             numItems             /* (I) */
)
{
    static char routine[] = "CxSurfaceMakeEmpty";
    int status = FAILURE;

    CxTSurface* p = NULL;

    REQUIRE(numItems > 0);

    p = NEW(CxTSurface);
    if (p==NULL) goto done;

    p->numItems        = numItems;

    p->x = NEW_ARRAY(double, p->numItems);
    if (p->x == NULL) goto done;

    p->y = NEW_ARRAY(double, p->numItems);
    if (p->y == NULL) goto done;

    p->z = NEW_ARRAY(double, p->numItems);
    if (p->z == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxSurfaceFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTSurface
***************************************************************************
*/
CxTSurface* CxSurfaceCopy(CxTSurface* src)
{
    CxTSurface* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxSurfaceMake(src->numItems,
                        src->x,
                        src->y,
                        src->z,
                        src->regularity);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTSurface
***************************************************************************
*/
void CxSurfaceFree(CxTSurface *p)
{
    if (p != NULL)
    {
        FREE(p->x);
        FREE(p->y);
        FREE(p->z);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpMake(
TBoolean        interpSurface,       /* (I) */
CxTSurfaceInterpMethod xInterp,             /* (I) */
CxTSurfaceExtrapMethod xExtrap,             /* (I) */
CxTSurfaceInterpMethod yInterp,             /* (I) */
CxTSurfaceExtrapMethod yExtrap              /* (I) */
)
{
    static char routine[] = "CxSurfaceInterpMake";
    int status = FAILURE;

    CxTSurfaceInterp* p = NULL;


    p = NEW(CxTSurfaceInterp);
    if (p==NULL) goto done;

    p->interpSurface   = interpSurface;
    p->xInterp         = xInterp;
    p->xExtrap         = xExtrap;
    p->yInterp         = yInterp;
    p->yExtrap         = yExtrap;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxSurfaceInterpFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpMakeEmpty(void)
{
    static char routine[] = "CxSurfaceInterpMakeEmpty";
    int status = FAILURE;

    CxTSurfaceInterp* p = NULL;


    p = NEW(CxTSurfaceInterp);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxSurfaceInterpFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTSurfaceInterp
***************************************************************************
*/
CxTSurfaceInterp* CxSurfaceInterpCopy(CxTSurfaceInterp* src)
{
    CxTSurfaceInterp* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxSurfaceInterpMake(src->interpSurface,
                              src->xInterp,
                              src->xExtrap,
                              src->yInterp,
                              src->yExtrap);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTSurfaceInterp
***************************************************************************
*/
void CxSurfaceInterpFree(CxTSurfaceInterp *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateMake(
int             n,                   /* (I) */
double*         x,                   /* (I) [n] */
TMatrix2D*      direction,           /* (I) */
int             iter,                /* (I) */
double          value,               /* (I) */
double          vdiff                /* (I) */
)
{
    static char routine[] = "CxMultiDimMinStateMake";
    int status = FAILURE;

    CxTMultiDimMinState* p = NULL;

    REQUIRE(n > 0);
    REQUIRE(x != NULL);

    p = NEW(CxTMultiDimMinState);
    if (p==NULL) goto done;

    p->n               = n;
    p->x = NEW_ARRAY(double, p->n);
    if (p->x == NULL) goto done;
    COPY_ARRAY (p->x, x, double, p->n);

    if (direction != NULL)
    {
        p->direction = GtoMatrixCopy(direction);
        if (p->direction == NULL) goto done;
    }
    else
    {
        p->direction = NULL;
    }

    p->iter            = iter;
    p->value           = value;
    p->vdiff           = vdiff;

    if (CxMultiDimMinStateValidate(p) != SUCCESS) goto done; /* failure */

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxMultiDimMinStateFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateMakeEmpty(
int             n                    /* (I) */
)
{
    static char routine[] = "CxMultiDimMinStateMakeEmpty";
    int status = FAILURE;

    CxTMultiDimMinState* p = NULL;

    REQUIRE(n > 0);

    p = NEW(CxTMultiDimMinState);
    if (p==NULL) goto done;

    p->n               = n;

    p->x = NEW_ARRAY(double, p->n);
    if (p->x == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxMultiDimMinStateFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTMultiDimMinState
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinStateCopy(CxTMultiDimMinState* src)
{
    CxTMultiDimMinState* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxMultiDimMinStateMake(src->n,
                                 src->x,
                                 src->direction,
                                 src->iter,
                                 src->value,
                                 src->vdiff);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTMultiDimMinState
***************************************************************************
*/
void CxMultiDimMinStateFree(CxTMultiDimMinState *p)
{
    if (p != NULL)
    {
        FREE(p->x);
        GtoMatrixFree(p->direction);
        FREE(p);
    }
}

