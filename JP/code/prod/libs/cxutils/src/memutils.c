/*
***************************************************************************
** SOURCE FILE: memutils.c
**
** Some basic and some complex memory utilities.
**
** Note that the memory cache code was written by such experts in software
** engineering as Doug Gallager and Simon Meldrum and hence is beyond all
** criticism!
**
** $Header$
***************************************************************************
*/

#include "memutils.h"

#include <string.h>  /* memset */

#include "cxmacros.h"

#include <alib/cmemory.h>

static CxTMemAllocFunc* g_memAllocFunc = GtoMallocSafe;
static CxTMemFreeFunc*  g_memFreeFunc  = GtoFreeSafe;

/*
***************************************************************************
** Allocates memory from the heap and initialises the allocated memory to
** zero.
**
** If there is no memory then returns NULL and tries to log an error.
***************************************************************************
*/
void* CxMemAlloc(size_t theSize)
{
    static char routine[]="CxMemAlloc";

    void *ptr;

    if (theSize <= 0)
    {
        GtoErrMsg ("%s: Number of bytes (%lu) must be at least 1.\n",
                 routine, (unsigned long) theSize);
        return NULL;
    }
    
    /* by default the next line calls GtoMallocSafe - but this can be
       overridden by calling CxMemFuncsSet */
    ptr = g_memAllocFunc(theSize);

    if (ptr == NULL)
    {
        GtoErrMsg("%s: Insufficient memory to allocate %lu bytes.\n", 
                routine, (unsigned long) theSize);
    }
    else
    {
        memset((char *)ptr, (size_t)0, theSize);
    }
    return(ptr);
}
    
/*
***************************************************************************
** Frees memory allocated via CxMemAlloc. Does nothing if the provided
** pointer is NULL.
***************************************************************************
*/
void CxMemFree(void *ptr)
{
    if (ptr != NULL)
    {
        g_memFreeFunc(ptr);
    }
}

/*
***************************************************************************
** Allocates a matrix which can be accessed by x[i][j].
** i can be in the range [0,numRows]
** j can be in the range [0,numCols]
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
 size_t elemSize)
{
    static char routine[] = "CxArray2DNew";

    size_t idx;
    size_t ptrBytes  = numRows * sizeof(void *);
    size_t dataBytes = numCols * numRows * elemSize;
    size_t totBytes;                       /* total # bytes to allocate */
    size_t rowSize;                        /* # bytes in a row */
    char **rowPtr;
    char *elemPtr;

/*
** We are going to allocate the memory in one big block.
**
** Note: we allocate enough pointers so that the pointers take up a multiple 
** of 8 bytes.
**
** This ensures that doubles can be accessed on 8 byte boundaries. On some
** platforms this is potentially very important.
**
** Hence we advance to the next 8-byte boundary when finding out how much
** space we need for the initial block of pointers.
*/

#undef DOUBLE_MASK
#define DOUBLE_MASK (sizeof(double)-1)
    ptrBytes = (ptrBytes+DOUBLE_MASK) & ~DOUBLE_MASK;
    totBytes = ptrBytes + dataBytes;
    
    /* Allocate all memory together-this reduces paging. */
    rowPtr = (char **)CxMemAlloc(totBytes);
    if (rowPtr == NULL)
    {
        GtoErrMsgFailure(routine);
        return NULL;
    }

    elemPtr = (char *)rowPtr + ptrBytes;
    rowSize = numCols * elemSize;
    for (idx=0; idx < numRows; idx++)
    {
        rowPtr[idx] = elemPtr + idx*rowSize;
    }

    return ((void **)rowPtr);
} 

/*
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
 CxTMemFreeFunc*  memFreeFunc)
{
    static char routine[] = "CxMemFuncsSet";

    if (memAllocFunc == NULL && memFreeFunc == NULL)
    {
        g_memAllocFunc = malloc;
        g_memFreeFunc  = free;
    }
    else if (memAllocFunc == NULL || memFreeFunc == NULL)
    {
        GtoErrMsg ("%s: Must set either both memAllocFunc and memFreeFunc or neither", routine);
        return FAILURE;
    }
    else
    {
        g_memAllocFunc = memAllocFunc;
        g_memFreeFunc  = memFreeFunc;
    }
    return SUCCESS;
}


/* Note that filler has been inserted here so that data ends up
 * on an 8 byte boundary (which is how long doubles are.) Otherwise
 * doubles cannot be accessed on the HP without using the +u4
 * compiler switch. This sacrifices a little memory for speed.
 */
#undef RESERVED_PTR_SIZE
#define RESERVED_PTR_SIZE sizeof(double) /* For HP alignment */

typedef struct _CACHE_ELEM              /* Used to put a memory elem */
{                                       /* on the free list. */
    struct _CACHE_ELEM *nextp;
    char filler[RESERVED_PTR_SIZE-sizeof(struct _CACHE_ELEM *)];
    char data[4];                       /* Real data goes here */
} CACHE_ELEM;  

typedef struct _CACHE_BLOCK             /* Needed only for freeing memory */
{                                       /* in the end. */
    struct _CACHE_BLOCK *nextp;
    char filler[RESERVED_PTR_SIZE-sizeof(struct _CACHE_BLOCK *)];
    char data[4];                       /* Real data goes here */
}  CACHE_BLOCK;  

struct _CxTMemCache
{
    size_t elemSize;                    /* # bytes in cache element */
    int elemsPerBlock;                  /* # blocks allocated at once */    
    CACHE_ELEM  *freeList;              /* List of available elements */
    CACHE_BLOCK *blockList;             /* List of allocated blocks */
    struct _CxTMemCache *masterListPrev;
    struct _CxTMemCache *masterListNext;
};

static void insertInFreeList(CxTMemCache *memCache, CACHE_ELEM *elemp);

static CACHE_ELEM inUse;
static CxTMemCache *masterListHead = NULL;
static CxTMemCache *masterListTail = NULL;

/* as a matter of academic curiosity, it's possible to construct a
   doubly-linked list for which the insert/delete operations involve
   no conditionals (use a circular list of link structures with the
   head/tail pointers being the next & prev of the always-present base
   link) */

/*
***************************************************************************
** Initializes memory cache.
***************************************************************************
*/
CxTMemCache* CxMemCacheInit(size_t elemSize, int elemsPerBlock)
{
    CxTMemCache *memCache = CX_NEW(CxTMemCache);
    if (memCache == NULL) return NULL; /* failure - no memory */

    memCache->freeList  = (CACHE_ELEM *)NULL;
    memCache->blockList = (CACHE_BLOCK *)NULL;
    memCache->elemSize  = elemSize + RESERVED_PTR_SIZE;
    memCache->elemsPerBlock = elemsPerBlock;

    if (masterListTail != NULL)
    {
        masterListTail->masterListNext = memCache; /* attach self as successor
                                                      of old tail */
    }
    else
    {
        masterListHead = memCache;  /* we are both head and tail */
    }
    memCache->masterListPrev = masterListTail;  /* predecessor is old tail 
                                                   (could be NULL) */
    memCache->masterListNext = NULL;            /* no successor */
    masterListTail = memCache;                  /* we become tail */

    return memCache;
}


/*
***************************************************************************
** Returns a free memory block
***************************************************************************
*/
void* CxMemCacheAlloc(CxTMemCache* memCache)
{
    CACHE_ELEM *elemp;
    CACHE_BLOCK *blockp;
    int count;                          /* Block counter */
    char *ptr;                          /* Generic pointer */
    
    /* If no memory on free list, allocate some. */
    if (memCache->freeList == (CACHE_ELEM *)NULL)
    {
        blockp = (CACHE_BLOCK*) CxMemAlloc(
            memCache->elemsPerBlock * memCache->elemSize + RESERVED_PTR_SIZE);
        if (blockp == NULL) return NULL; /* failure */

        /* Insert block in blockList */
        blockp->nextp = memCache->blockList;
        memCache->blockList = blockp;

        /* Now insert cache elements in freeList. */
        for (count = memCache->elemsPerBlock, ptr = (char *)blockp->data;
             count > 0;
             count--, ptr += memCache->elemSize)
        {
             insertInFreeList(memCache, (CACHE_ELEM *)ptr);
        }
    }
    
    elemp = memCache->freeList;            /* Get first in list. */
    memCache->freeList = elemp->nextp;     /* Remove from list */
    
    elemp->nextp = &inUse;        /* enables identification of used elements */
    return ((void *)elemp->data); /* Return allocated memory */
}


/*
***************************************************************************
** Puts a memory block back in free list.
***************************************************************************
*/
void CxMemCacheFree(CxTMemCache* memCache, void *dataPtr)
{
    insertInFreeList(memCache,
                     (CACHE_ELEM *)((char *)dataPtr - RESERVED_PTR_SIZE));
}

/*
***************************************************************************
** Frees all memory allocated.
***************************************************************************
*/
void CxMemCacheFreeAll(CxTMemCache* memCache)
{
    CACHE_BLOCK *blockp, *nextBlockp;

    if (memCache == NULL) return;

    /* Deallocate all blocks. */
    blockp = memCache->blockList;
    nextBlockp = blockp;
    for (;
         blockp != (CACHE_BLOCK *)NULL;
         blockp = nextBlockp)
    {
        nextBlockp = blockp->nextp;
        CX_FREE(blockp);
    }
         
    /* drop from master list */
    if ( memCache->masterListPrev == NULL )
    { /* we were head */
        masterListHead                          = memCache->masterListNext;
    }
    else
    {
        (memCache->masterListPrev)->masterListNext = memCache->masterListNext;
    }
    
    if ( memCache->masterListNext == NULL )
    { /* we were tail */
        masterListTail                          = memCache->masterListPrev;
    }
    else
    {
        (memCache->masterListNext)->masterListPrev = memCache->masterListPrev;
    }

    /* Now free the Head block. */
    CX_FREE(memCache);     
}


/*
***************************************************************************
** Frees unused memory. Typically a free list consumes memory as required
** to satisfy allocation requests, but never releases it.
***************************************************************************
*/
void CxMemCachePurge(CxTMemCache* memCache)
{
    int blocksDeleted = 0;

    /* (1) remove unused blocks */
    CACHE_BLOCK **blkRef;
    
    /* iterate over blocks, but using an indirect iterator which
       permits easy deletion */
    blkRef = &memCache->blockList;
    while ( *blkRef != NULL )
    {
        int n;
        char *el;
        
        /* iterate over elements */
        for (n = 0, el = (*blkRef)->data;
             n < memCache->elemsPerBlock;
             n++, el += memCache->elemSize)
        {
            if ( ((CACHE_ELEM *) el)->nextp == &inUse )
            {
                break;
            }
        }
        
        /* no elements were in use? */
        if ( n == memCache->elemsPerBlock )
        {
            CACHE_BLOCK *nextp = (*blkRef)->nextp;
            
            CX_FREE(*blkRef);
            
            /* drop current block from list, by patching
               its "next" pointer into "prev->next" */
            *blkRef = nextp;
            
            blocksDeleted++;
        }
        else
        {
            blkRef = &((*blkRef)->nextp); /* step on */
        }
    }

    /* (2) reconstruct free element linked list */
    if ( blocksDeleted )
    {
        CACHE_BLOCK *blk;

        memCache->freeList = NULL;

        /* iterate over blocks */
        for (blk = memCache->blockList; blk != NULL; blk = blk->nextp)
        {
            int n;
            char *el;
            
            /* iterate over elements */
            for (n = 0, el = blk->data;
                 n < memCache->elemsPerBlock;
                 n++, el += memCache->elemSize)
            {
                if ( ((CACHE_ELEM *) el)->nextp != &inUse )
                {
                    insertInFreeList(memCache, (CACHE_ELEM *) el);
                }
            }
        }
    }
}

/*
***************************************************************************
** Frees unused memory from all caches.
***************************************************************************
*/
void CxMemCachePurgeAll(void)
{
    CxTMemCache *memCache;

    for (memCache = masterListHead;
         memCache != NULL;
         memCache = memCache->masterListNext)
    {
        CxMemCachePurge(memCache);
    }
}

/*
***************************************************************************
** insertInFreeList(): Inserts an element in free list.
***************************************************************************
*/
static void insertInFreeList(CxTMemCache *memCache, CACHE_ELEM *elemp)
{
    elemp->nextp = memCache->freeList;
    memCache->freeList = elemp;
}



