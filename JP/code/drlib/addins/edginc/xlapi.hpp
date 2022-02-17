///////////////////////////////////////////////////////////////////////
//    Group       : EDG DR
//
//    Filename    : xlapi.hpp
//
//    Description : Interface for excel. Includes Microsoft's published file
//                  together with definition of OPER structure. For non
//                  windows, provides alternative definitions of the 
//                  necessary structures
//
//    Author      : Mark A Robson
//
//    Date        : 17 May 2000
//
//
///////////////////////////////////////////////////////////////////////

#ifndef EDG_XLAPI_H 
#define EDG_XLAPI_H
             
#if defined(_WIN32) || defined (__CYGWIN32__)
#if 0
/* I think the latest version of gcc (Win32) doesn't need this */
#ifdef __CYGWIN32__
#undef PASCAL
#undef far
#undef near
#define far
#define near
#endif
#endif
/* MS-Windows API - avoid including as it's incredibly slow */
//#include "windows.h"
typedef char CHAR;
typedef unsigned short WORD;
typedef CHAR *LPSTR;
typedef unsigned char BYTE;
typedef unsigned long DWORD;
typedef void* HANDLE;
#define FAR
#define far
#define pascal __stdcall

/* MS-Windows Excel API (unfortunately missing the OPER structure) */   
#include "edginc/xlcall.hpp"    /* Must come after windows.h */
/*
** EDG_XLFUNC macro defines the type of function which can be called from
** Excel. For example, XLFUNC (LPXLOPER) expands to
**     __declspec(dllexport) LPXLOPER
*/
#define EDG_XLFUNC(retType) __declspec(dllexport) retType

/* try and undo some of the dubious defines that MS have set up */
#ifdef interface
#undef interface
#endif

/*
** OPER structure - same shape as XLOPER but less types defined in the union
**
** This structure was originally defined in an appendix to one of the Excel manuals
** before the release of the Excel development kit (which was for Excel4)
*/

/*
** OPER structure - same shape as XLOPER but less types defined in the union
**
** This structure was originally defined in an appendix to one of the Excel
** manuals before the release of the Excel development kit (which was for
** Excel4)
*/

typedef struct _oper
{
    union
    {
        double num;                     /* xltypeNum */
        LPSTR str;                      /* xltypeStr */
        WORD boolVal;                   /* xltypeBool */
        WORD err;                       /* xltypeErr */
        struct 
        {
            struct _oper *lparray;
            WORD rows;
            WORD columns;
        } xlarray;                        /* xltypeMulti */
    } val;
    WORD type; /* one of xltypeNum, xltypeStr, xltypeBool, xltypeErr,
                  xltypeMulti, xltypeMissing, xltypeNil according to MS docs */

} OPER, *LPOPER ;

#else
/* allow compilation so tests can run on unix */
#define EDG_XLFUNC(retType) retType
typedef struct _oper
{
    union
    {
        double num;
        char*  str;
        unsigned short int boolVal;
        unsigned short int err;
        struct
        {
            struct _oper *lparray;
            unsigned short int rows;
            unsigned short int columns;
        } xlarray;
    } val;
    unsigned short int type;
} OPER, *LPOPER;

/* XLREF structure 
**
** Describes a single rectangular reference
*/

typedef struct xlref 
{
    unsigned short int rwFirst;
    unsigned short int rwLast;
    unsigned char      colFirst;
    unsigned char      colLast;
} XLREF, *LPXLREF;


/*
** XLMREF structure
**
** Describes multiple rectangular references.
** This is a variable size structure, default 
** size is 1 reference.
*/

typedef struct xlmref 
{
    unsigned short int count;
    XLREF reftbl[1];                        /* actually reftbl[count] */
} XLMREF, *LPXLMREF;


/*
** XLOPER structure 
**
** Excel's fundamental data type: can hold data
** of any type. Use "R" as the argument type in the 
** REGISTER function.
**/

typedef struct xloper 
{
    union 
    {
        double num;                     /* xltypeNum */
        char*  str;                     /* xltypeStr */
        unsigned short int boolVal;     /* xltypeBool */
        unsigned short int err;         /* xltypeErr */
        short int w;                    /* xltypeInt */
        struct 
        {
            unsigned short int count;   /* always = 1 */
            XLREF ref;
        } sref;                         /* xltypeSRef */
        struct 
        {
            XLMREF*         lpmref;
            unsigned long   idSheet;
        } mref;                         /* xltypeRef */
        struct 
        {
            struct xloper*      lparray;
            unsigned short int  rows;
            unsigned short int  columns;
        } xlarray;                        /* xltypeMulti */
        struct 
        {
            union
            {
                short int     level;        /* xlflowRestart */
                short int     tbctrl;       /* xlflowPause */
                unsigned long idSheet;      /* xlflowGoto */
            } valflow;
            unsigned short int rw;          /* xlflowGoto */
            unsigned char      col;         /* xlflowGoto */
            unsigned char      xlflow;
        } flow;                         /* xltypeFlow */
        struct
        {
            union
            {
                unsigned char*  lpbData;      /* data passed to XL */
                void *          hdata;        /* data returned from XL */
            } h;
            long cbData;
        } bigdata;                      /* xltypeBigData */
    } val;
    unsigned short int type;
} XLOPER, *LPXLOPER;

/*
** XLOPER data types
**
** Used for xltype field of XLOPER structure
*/
#define xltypeNum        0x0001
#define xltypeStr        0x0002
#define xltypeBool       0x0004
#define xltypeRef        0x0008
#define xltypeErr        0x0010
#define xltypeFlow       0x0020
#define xltypeMulti      0x0040
#define xltypeMissing    0x0080
#define xltypeNil        0x0100
#define xltypeSRef       0x0400
#define xltypeInt        0x0800
#define xlbitXLFree      0x1000
#define xlbitDLLFree     0x4000

#define xltypeBigData    (xltypeStr | xltypeInt)

/*
** Error codes
**
** Used for val.err field of XLOPER structure
** when constructing error XLOPERs
*/

#define xlerrNull    0
#define xlerrDiv0    7
#define xlerrValue   15
#define xlerrRef     23
#define xlerrName    29
#define xlerrNum     36
#define xlerrNA      42

/** Return codes */
#define xlretSuccess    0     /* success */ 
#define xlretAbort      1     /* macro halted */
#define xlretInvXlfn    2     /* invalid function number */ 
#define xlretInvCount   4     /* invalid number of arguments */ 
#define xlretInvXloper  8     /* invalid OPER structure */  
#define xlretStackOvfl  16    /* stack overflow */  
#define xlretFailed     32    /* command failed */  
#define xlretUncalced   64    /* uncalced cell */

/* Function number bits */ 
#define xlCommand    0x8000
#define xlSpecial    0x4000
#define xlIntl       0x2000
#define xlPrompt     0x1000

/* Auxiliary function numbers */
#define xlFree             (0  | xlSpecial)
#define xlStack            (1  | xlSpecial)
#define xlCoerce           (2  | xlSpecial)
#define xlSet              (3  | xlSpecial)
#define xlSheetId          (4  | xlSpecial)
#define xlSheetNm          (5  | xlSpecial)
#define xlAbort            (6  | xlSpecial)
#define xlGetInst          (7  | xlSpecial)
#define xlGetHwnd          (8  | xlSpecial)
#define xlGetName          (9  | xlSpecial)
#define xlEnableXLMsgs     (10 | xlSpecial)
#define xlDisableXLMsgs    (11 | xlSpecial)
#define xlDefineBinaryName (12 | xlSpecial)
#define xlGetBinaryName    (13 | xlSpecial)

#endif

// All files are written using XL_OPER  - allows easier switch between 
// OPER and XLOPER
typedef XLOPER XL_OPER;

#endif /* XLAPI_H */
