/*
***************************************************************************
** HEADER FILE: memutils.h
**
** Some basic and some complex memory utilities.
**
** $Header$
***************************************************************************
*/

#ifndef CX_MEMUTILS_H
#define CX_MEMUTILS_H

#include <alib/cgeneral.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef void* (CxTMemAllocFunc) (size_t);
typedef void (CxTMemFreeFunc) (void*);

#define CX_NEW(T)                (T*)(CxMemAlloc(sizeof(T)))
#define CX_NEW_ARRAY(T,n)        (T*)(CxMemAlloc(sizeof(T)*(n)))
#define CX_NEW_ARRAY_2D(T,nr,nc) (T**)(CxArray2DNew((nr),(nc),sizeof(T)))
#define CX_FREE(p)               CxMemFree(p)

/* frees an array p which starts at s instead of 0 */
#define CX_FREE_ARRAY(p,s) if ((p)!=NULL) CX_FREE((p)+(s))

/* frees an array of pointers with a free function for each pointer */
#define CX_FREE_PTR_ARRAY(array,size,freeFunc) \
do {\
    if (array != NULL && size > 0)\
    {\
        long ARRAY_I;\
        long ARRAY_N = (long)(size);\
        for (ARRAY_I = 0; ARRAY_I < ARRAY_N; ++ARRAY_I)\
        {\
            freeFunc((array)[ARRAY_I]);\
        }\
    }\
    CX_FREE(array);\
} while(0)


/*f
***************************************************************************
** Allocates memory from the heap and initialises the allocated memory to
** zero.
**
** If there is no memory then returns NULL and tries to log an error.
***************************************************************************
*/
void* CxMemAlloc(size_t theSize);
    
/*f
***************************************************************************
** Frees memory allocated via CxMemAlloc. Does nothing if the provided
** pointer is NULL.
***************************************************************************
*/
void CxMemFree(void *ptr);

/*
***************************************************************************
** Allocates a matrix which can be accessed by x[i][j].
**
** The matrix is 0-based for both rows and columns, and cannot be done
** otherwise.
**
** The memory is allocated in a single block and can be free'd by a single
** FREE statement.
**
** Example use where you want a matrix of doubles (the usual case):
**
** double **matrix = (double**)CxArray2DNew(numRows, numCols, sizeof(double))
**
** Alternatively you can use the CX_NEW_ARRAY_2D macro which does this
***************************************************************************
*/
void** CxArray2DNew
(size_t numRows,
 size_t numCols,
 size_t elementSize);

/*f
***************************************************************************
** Allows the application to override the default memory allocation and
** free routines.
**
** The idea is that when plugged into an interface layer we may wish to
** use common memory allocation and free routines.
**
** Returns 0 (CX_SUCCESS) or -1 (CX_FAILURE).
***************************************************************************
*/
int CxMemFuncsSet
(CxTMemAllocFunc* memAllocFunc,
 CxTMemFreeFunc*  memFreeFunc);

/*t
***************************************************************************
** Defines a memory cache. We use a memory cache when we want to allocate
** lots of small blocks for a particular data type. The data type needs to
** be a fixed size.
***************************************************************************
*/
typedef struct _CxTMemCache CxTMemCache;

/*f
***************************************************************************
** Initialises the memory cache.
***************************************************************************
*/
CxTMemCache* CxMemCacheInit
(size_t elemSize,      /* (I) Size of a single element */
 int    elemsPerBlock  /* (I) Number of elements that we allocate when
                          needing a new memory block */
);

/*f
***************************************************************************
** Returns a block of memory.
***************************************************************************
*/
void* CxMemCacheAlloc(CxTMemCache* memCache);

/*f
***************************************************************************
** Release a block of memory so that it can be re-used.
***************************************************************************
*/
void CxMemCacheFree(CxTMemCache* memCache, void *dataPtr);

/*f
***************************************************************************
** Frees all memory allocated within a memory cache.
***************************************************************************
*/
void CxMemCacheFreeAll(CxTMemCache* memCache);

/*f
***************************************************************************
** Frees unused memory. Typically a free list consumes memory as required
** to satisfy allocation requests, but never releases it.
***************************************************************************
*/
void CxMemCachePurge(CxTMemCache* memCache);

/*f
***************************************************************************
** Frees unused memory from all caches.
***************************************************************************
*/
void CxMemCachePurgeAll(void);


/*
***************************************************************************
** If you want to use the memory cache for a particular data type, then
** use this macro within your source file for that particular data type.
**
** For example:
**
** CX_DECLARE_AUTO_MEMCACHE(someType, CXSomeType)
**
** You will get the following static functions in your source file:
**
** static CXSomeType* someTypeNew(void);
** static void someTypeFree(CXSomeType* ptr);
**
** Memory allocations are done in blocks of 64 thus reducing the number
** of mallocs and making memory allocation more efficient.
***************************************************************************
*/
#define CX_DECLARE_AUTO_MEMCACHE(objIdentifier, objType)           \
                                                                   \
static CxTMemCache* g_##objIdentifier##FreeList = NULL;             \
static int g_##objIdentifier##FreeListNumAlloced = 0;              \
                                                                   \
static objType *objIdentifier##New(void)                           \
{                                                                  \
    objType *obj = NULL;                                           \
                                                                   \
    if ( g_##objIdentifier##FreeList == NULL )                     \
    {                                                              \
        g_##objIdentifier##FreeList =                              \
            CxMemCacheInit(sizeof(objType) /* element size */,     \
                           64              /* granularity */);     \
        if ( g_##objIdentifier##FreeList == NULL )                 \
        {                                                          \
            goto done;                                             \
        }                                                          \
    }                                                              \
                                                                   \
    obj = (objType *)CxMemCacheAlloc(g_##objIdentifier##FreeList); \
    if ( obj != NULL )                                             \
    {                                                              \
        g_##objIdentifier##FreeListNumAlloced++;                   \
    }                                                              \
                                                                   \
done:                                                              \
    /* normally when there are no allocations, free list is     */ \
    /* NULL; if we are recovering from an error, we may violate */ \
    /* this convention, so remedy below:                        */ \
                                                                   \
    if ( g_##objIdentifier##FreeListNumAlloced == 0 &&             \
         g_##objIdentifier##FreeList != NULL )                     \
    {                                                              \
        CxMemCacheFreeAll(g_##objIdentifier##FreeList);            \
        g_##objIdentifier##FreeList = NULL;                        \
    }                                                              \
                                                                   \
    return obj;                                                    \
}                                                                  \
                                                                   \
static void objIdentifier##Free(objType *obj)                      \
{                                                                  \
    if ( obj != NULL )                                             \
    {                                                              \
        CxMemCacheFree(g_##objIdentifier##FreeList, obj);          \
        --g_##objIdentifier##FreeListNumAlloced;                   \
        if ( g_##objIdentifier##FreeListNumAlloced == 0 )          \
        {                                                          \
            CxMemCacheFreeAll(g_##objIdentifier##FreeList);        \
            g_##objIdentifier##FreeList = NULL;                    \
        }                                                          \
    }                                                              \
}

#ifdef __cplusplus
}
#endif

#endif
