/*
***************************************************************************
** SOURCE FILE: hash.c
**
** Defines a generic hash table, plus some utilities for manipulating
** simple hash tables where the hash key is a string or an integer (long).
**
** Written by such luminaries from J.P.Morgan's past such as Doug Gallager
** and Rob Norman, so this code is beyond all criticism.
***************************************************************************
*/

#include "irx/hash.h"

#include <string.h>

#include "irx/error.h"
#include "irx/macros.h"
#include "irx/memutils.h"    

/*
** Not sure why this stuff is not in a header file yet
*/
typedef struct _IrxTHashElement
{
    void* key;
    void* data;
    struct _IrxTHashElement *next;
} IrxTHashElement;

IRX_DECLARE_AUTO_MEMCACHE(hashElement, IrxTHashElement)

struct _IrxTHashTable
{
    long                 bucketIdx;  /* Used for irxHashFirst, irxHashNext*/
    IrxTHashElement*      nextElem;   /* Used for irxHashFirst, irxHashNext*/
    IrxTHashElement**     buckets;    /* Array of pointers to IrxTHashElements */
    IrxTBool             uniqueKey;  /* Keys must be unique */
    long                 numBuckets;   /* Number of buckets */
    IrxTHashFunc*         hash;         /* Hash function */
    IrxTHashKeyCmpFunc*   hashKeyCmp;   /* Hash key comparison function */
    IrxTHashKeyCopyFunc*  hashKeyCopy;  /* Hash key copy function */
    IrxTHashKeyFreeFunc*  hashKeyFree;  /* Hash key free function */
    IrxTHashDataFreeFunc* hashDataFree; /* Hash data free function */
};


/* 
 * Used if caller specifies a NULL print function on print.
 */ 

static long  hashString (char*, long);
static char* hashCopyString (char*);
static void  printPointer(void* pointer);
static long  hashLong (void*, long);
static int   hashCompareLong (void*, void*);
static void* hashCopyLong (void*);
static void  hashFreeLong (void*);

/*f
***************************************************************************
** Creates a hash table where the key used is a long.
**
** Note that this is possible since long and void* are the same size,
** and although normally the hash table expects pointer keys we can
** therefore safely use long keys instead.
***************************************************************************
*/
IrxTHashTable* irxHashTableMakeLongKey
(long                numBuckets,   /* (I) Length of the table */
 IrxTBool              uniqueKey,    /* (I) Keys need to be unique */
 IrxTHashDataFreeFunc* hashDataFree  /* (I) Hash data free function 
                                      (NULL=don't free)*/
)
{
    static char routine[] = "irxHashTableMakeLongKey";
    
    IrxTHashTable* hashTable = NULL;

    if (sizeof(long) != sizeof(void*))
    {
        irxError ("%s: Cannot create a hash table with long keys because "
                   "sizeof(long) %ld is not equal to sizeof(void*) %ld\n",
                   routine, (long)(sizeof(long)), (long)(sizeof(void*)));
        return NULL;
    }

    hashTable = irxHashTableMake (numBuckets,
                                 uniqueKey,
                                 hashLong,
                                 hashCompareLong,
                                 hashCopyLong,
                                 hashFreeLong,
                                 hashDataFree);
    return hashTable;
}

/*f
***************************************************************************
** Creates a hash table where the key used is a string.
***************************************************************************
*/
IrxTHashTable* irxHashTableMakeStringKey
(long                numBuckets,   /* (I) Length of the table */
 IrxTBool              uniqueKey,    /* (I) Keys need to be unique */
 IrxTHashDataFreeFunc* hashDataFree  /* (I) Hash data free function 
                                      (NULL=don't free)*/
)
{
    IrxTHashTable* hashTable = NULL;

    hashTable = irxHashTableMake (numBuckets,
                                 uniqueKey,
                                 (IrxTHashFunc*)hashString,
                                 (IrxTHashKeyCmpFunc*)strcmp,
                                 (IrxTHashKeyCopyFunc*)hashCopyString,
                                 irxMemFree,
                                 hashDataFree);

    return hashTable;
}

/*f
***************************************************************************
** Creates a generic hash table.
**
** The internal data structure is exceedingly complex, which is why the
** internal structure of struct _IrxTHashTable is not revealed in this header
** file.
**
** You can think of a hash table as a large array where the index of the
** the array is an arbitrary data type (not necessarily a number).
**
** In practice a hash table is implemented (as its name suggests) by
** using a hashing algorithm to convert the key into a particular bucket,
** and then doing a linear search through the bucket to find the key.
**
** Therefore for the hash table to be efficient you need to define the
** number of buckets to be a reasonable number compared to the expected
** number of entries in the hash table, and also for the hashing algorithm
** to assign the buckets in a relatively uniform manner.
**
** The key is an arbitrary data type. The hash table needs to be able to
** compare keys, hash keys, duplicate keys and free keys.
**
** The data is also an arbitrary data type. When the hash table is deleted
** you may need to free the data within the hash table (this depends on
** your implied ownership policy for the data within the hash table).
**
** If you need to free the data, then supply a free function for the data.
** Otherwise the free function for the data can be NULL.
***************************************************************************
*/
IrxTHashTable* irxHashTableMake
(long                numBuckets,   /* (I) Length of the table */
 IrxTBool              uniqueKey,    /* (I) Keys need to be unique */
 IrxTHashFunc*         hashFunc,     /* (I) Hash function */
 IrxTHashKeyCmpFunc*   hashKeyCmp,   /* (I) Hash key compare func */
 IrxTHashKeyCopyFunc*  hashKeyCopy,  /* (I) Hash key copy func */
 IrxTHashKeyFreeFunc*  hashKeyFree,  /* (I) Hash key free function */
 IrxTHashDataFreeFunc* hashDataFree  /* (I) Hash data free function 
                                      (NULL=don't free)*/
)
{
    static char routine[] = "irxHashInit";
    IrxTHashTable *hashTable = NULL;

    int status = FAILURE;

    REQUIRE(numBuckets > 0);
    REQUIRE(hashFunc != NULL);
    REQUIRE(hashKeyCmp != NULL);
    REQUIRE(hashKeyCopy != NULL);
    REQUIRE(hashKeyFree != NULL);
    
    hashTable = IRX_NEW(IrxTHashTable);
    if (hashTable == NULL) goto RETURN;

    hashTable->hash         = hashFunc;
    hashTable->uniqueKey    = uniqueKey;
    hashTable->hashKeyCmp   = hashKeyCmp;
    hashTable->hashKeyCopy  = hashKeyCopy;
    hashTable->hashKeyFree  = hashKeyFree;
    hashTable->hashDataFree = hashDataFree;

    /* Allocate array of buckets. */
    hashTable->numBuckets = numBuckets;
    hashTable->buckets  = IRX_NEW_ARRAY(IrxTHashElement*, numBuckets);
    if (hashTable->buckets == NULL) goto RETURN;

    status = SUCCESS;

  RETURN:

    if (status != SUCCESS)
    {
        irxError ("%s: Failed.\n", routine);
        irxHashTableFree (hashTable);
        hashTable = NULL;
    }

    return hashTable;
}


/*f
***************************************************************************
** Inserts a new key,value pair into the hash table.
**
** The data is shallow copied. The key is copied using the register
** key copy function.
**
** If keys need to be unique, then it is an error to insert the same key
** twice. If keys do not need to be unique, then inserting the same key
** twice will put duplicate keys into the hash table.
***************************************************************************
*/
int irxHashTableInsert
(IrxTHashTable* hashTable,   /* (I) Hash table */
 void*  key,              /* (I) Identifier */
 void*  data)             /* (I) Data to be stored */
{
    static char routine[] = "irxHashTableInsert";
    int status = FAILURE;        /* Status */

    IrxTHashElement  *listElement;      /* Pointer to elements in list */
    IrxTHashElement  *newElement;       /* Pointer to new element to add */
    long bkt;                         /* Bucket index */

    REQUIRE(key != NULL);
    
    /* Get the hash bucket. */
    bkt = (hashTable->hash)(key, hashTable->numBuckets);

    /* Check for uniqueness */
    if (hashTable->uniqueKey)
    {
        for (listElement = hashTable->buckets[bkt];
             listElement != NULL;
             listElement = listElement->next)
        {
            if ((hashTable->hashKeyCmp)(listElement->key, key) == 0)
            {
                irxError ("%s: Key is not unique\n", routine);
                goto RETURN; /* failure */
            }
        }
    }

    /* Point to 1st element in list. */
    listElement = hashTable->buckets[bkt];

    /* Allocate a new element to insert in hash table. */
    newElement = hashElementNew();
    if (newElement == NULL) goto RETURN; /* failure */

    /* Link new element into the head of the list. */
    hashTable->buckets[bkt] = newElement;
    newElement->next = listElement;
    
    /* Set new element */
    newElement->key = (hashTable->hashKeyCopy)(key);
    newElement->data = data;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure (routine);
    return status;
}


/*f
***************************************************************************
** Returns data given the key.
**
** If there is more than one entry with the same key, then the last one
** entered into the hash table is returned.
**
** Returns -1 (IRX_FAILURE) if the element cannot be found - but does not
** call the error handler in this case, since it may be completely normal
** to look-up in a hash table without knowing if the item is in the hash
** table. Otherwise returns 0 (IRX_SUCCESS).
**
** The output value is returned as is - i.e. as it was passed into the
** hash table using irxHashTableInsert.
***************************************************************************
*/
int irxHashTableSearch
(IrxTHashTable*  hashTable,        /* (I) Hash table identifier */
 void*         key,              /* (I) Key */
 void**        data)             /* (O) Address of pointer to data */
{
    IrxTHashElement *listElement;  /* Pointer to elements in list */
    long bkt;                    /* Bucket index */

    /* Get the hash bucket. */
    bkt = (hashTable->hash)(key, hashTable->numBuckets);

    /* Look for the element in the list with a matching key. */
    for (listElement = hashTable->buckets[bkt];
         listElement != NULL;
         listElement = listElement->next)
    {
        if ((hashTable->hashKeyCmp)(listElement->key, key) == 0)
        {
            /* found it */
            *data = listElement->data;
            return SUCCESS;
        }
    }

    return FAILURE; /* not found */
}


/*f
***************************************************************************
** Deletes an entry in the hash table for the given key.
**
** Note that the corresponding data item will be deleted if a free function
** was registered for it.
**
** The key is always deleted.
**
** Returns -1 (IRX_FAILURE) if the entry is not in the hash table, and in
** this case the error handler is called. Otherwise returns 0 (IRX_SUCCESS).
***************************************************************************
*/
int irxHashTableDeleteEntry
(IrxTHashTable* hashTable,               /* (I) Hash table */
 void*        key)                     /* (I) Key */
{
    static char routine[] = "irxHashTableDeleteEntry";

    void *data;

    if (irxHashTableExtract(hashTable, key, &data) != SUCCESS)
    {
        irxErrorFailure(routine);
        return FAILURE;
    }

    if (hashTable->hashDataFree != NULL && data != NULL)
        hashTable->hashDataFree(data);

    return SUCCESS;
}


/*f
***************************************************************************
** Returns data and removes it from the hash table.
**
** Returns -1 (IRX_FAILURE) if the entry is not in the hash table, and in
** this case the error handler is called. Otherwise returns 0 (IRX_SUCCESS).
**
** The key is deleted within the hash table, but the data is not deleted
** (otherwise we could not return it).
***************************************************************************
*/
int irxHashTableExtract
(IrxTHashTable* hashTable,            /* (I) Hash table */
 void*        key,                  /* (I) Key */
 void**       data)                 /* (O) Pointer to data */
{
    char    routine[] = "irxHashTableExtract";

    IrxTHashElement *listElement;        /* Pointer to elements in list */
    IrxTHashElement *previousElement;    /* Pointer to previous one in list */
    long bkt;                          /* Bucket index */

    /* Get the hash bucket.*/
    bkt = (hashTable->hash)(key, hashTable->numBuckets);

    previousElement = hashTable->buckets[bkt];

    /* The first element in the list is a special case. */
    if ((hashTable->hashKeyCmp)(previousElement->key, key) == 0 )
    {
        hashTable->buckets[bkt] = previousElement->next;
        hashTable->hashKeyFree(previousElement->key);
        *data = previousElement->data;
        hashElementFree(previousElement);
        return SUCCESS;
    }
 
    /* Look for the element in the list with a matching key. */
    for (listElement = previousElement->next;
         listElement != NULL;
         previousElement = listElement, listElement = listElement->next)
    {
        if ((hashTable->hashKeyCmp)(listElement->key, key) == 0)
        {
            /* Found it - so remove it from the list */
            previousElement->next = listElement->next;
            hashTable->hashKeyFree(listElement->key);
            *data = listElement->data;
            hashElementFree(listElement);
            return SUCCESS;
        }
    }

    irxError("%s: Unable to find key.\n", routine);
    return FAILURE;
}

/*f
***************************************************************************
** Prints out the contents of the hash table using the provided print
** functions for key and data.
***************************************************************************
*/
void irxHashTablePrint
(IrxTHashTable*         hashTable,      /* (I) Hash table identifier */ 
 IrxTHashKeyPrintFunc*  hashKeyPrint,   /* (I) Hash key print function */ 
 IrxTHashDataPrintFunc* hashDataPrint   /* (I) Hash data print function */ 
)
{
    IrxTHashElement *element;            /* Element in the list */
    long bkt;                          /* Bucket index */

    if (hashKeyPrint == NULL) hashKeyPrint = printPointer;
    if (hashDataPrint == NULL) hashDataPrint = printPointer;

    if (hashTable == NULL)
    {
        irxError ("\nHash Table: (NULL)\n");
        return;
    }
    irxError("\nHash Table:\n");

    /* For each bucket... */
    for (bkt=0; bkt < hashTable->numBuckets; bkt++)
    {
        if (hashTable->buckets[bkt] != NULL)
            irxError("Bucket: %ld\n", bkt);

        for (element = hashTable->buckets[bkt];
             element != NULL;
             element = element->next)            
        {
            irxError("\tKEY: ");
            hashKeyPrint(element->key);

            irxError(" DATA: ");
            hashDataPrint(element->data);

            irxError("\n");
        }
    }
}

/*f
***************************************************************************
** Returns the 1st element in the hash table. Returns TRUE if there are
** any elements in the hash table, and FALSE if there are no elements in
** the hash table.
**
** Designed to be part of an iteration to be used in conjunction with
** irxHashTableNext.
**
** For example:
**
** IrxTBool more = irxHashTableFirst (ht, &key, &data);
** while (more)
** {
**     process key and data
**     more = irxHashTableNext (ht, &key, &data);
** }
**
** This function does not affect the contents of the hash table.
**
** The key copy function is not called when setting the output of
** the key. In fact the user should treat the key and data returned with
** great reverence - they are effectively const although not declared
** as such.
**
** Note that irxHashTableFirst and irxHashTableNext are not a thread-safe
** combination - you should not iterate through the same hash table in
** different threads.
***************************************************************************
*/
IrxTBool irxHashTableFirst
(IrxTHashTable* hashTable,             /* (I) Hash table Id */
 void**       key,                   /* (O) Key for next elem */
 void**       data)                  /* (O) Data for 1st elem */
{
    if (hashTable == NULL) return FALSE;

    hashTable->bucketIdx = 0;
    hashTable->nextElem  = hashTable->buckets[0];

    return irxHashTableNext(hashTable, key, data);
}


/*f
***************************************************************************
** Returns the next element in the hash table. Returns TRUE if there is
** an another element, and FALSE if there are no more elements.
**
** See irxHashTableFirst for further details and caveats - the function
** has the same prototype and the outputs should be treated in the same
** manner.
***************************************************************************
*/
IrxTBool irxHashTableNext
(IrxTHashTable* hashTable,             /* (I) Hash table Id */
 void**       key,                   /* (O) Key for next elem */
 void**       data)                  /* (O) Data for next elem */
{
    IrxTHashElement *nextElem;
    long bkt;                     /* Bucket index */

    if (hashTable == NULL) return FALSE;

    if (hashTable->nextElem == NULL)
    {
        /* Look through the buckets. */
        for (bkt = hashTable->bucketIdx+1;
             bkt < hashTable->numBuckets && hashTable->buckets[bkt] == NULL;
             bkt++)
        {
            /* Just searching for non-empty bucket */;
        }
        
        /* Check if no more left. */
        if (bkt == hashTable->numBuckets)
        {
            return FALSE;
        }

        hashTable->nextElem  = hashTable->buckets[bkt];
        hashTable->bucketIdx = bkt;
    }

    /* For convenience. */
    nextElem = hashTable->nextElem;

    /* Set outputs. */
    *key  = nextElem->key;
    *data = nextElem->data;

    /* Update for next time. */
    hashTable->nextElem = nextElem->next;

    return TRUE;
}


/*f
***************************************************************************
** Deallocates the hash table, including all keys and data.
**
** Note that the data is only free'd if the user provided a free function
** when constructing the hash table.
***************************************************************************
*/
void irxHashTableFree (IrxTHashTable* hashTable)
{
    if (hashTable != NULL)
    {
        IrxTHashElement *element;  /* Element in the list */
        long bkt;                /* Bucket index */

        for (bkt=0; bkt < hashTable->numBuckets; ++bkt)
        {
            while (hashTable->buckets[bkt] != NULL)
            {
                element = hashTable->buckets[bkt];
                hashTable->hashKeyFree(element->key);
                if (hashTable->hashDataFree)
                    hashTable->hashDataFree(element->data);
                /* tidy-up as we go in case we are interrupted */
                /* admittedly that seems somewhat unlikely :-) */
                hashTable->buckets[bkt] = element->next;
                hashElementFree(element);
            }
        }
        IRX_FREE(hashTable->buckets);
        IRX_FREE(hashTable);
    }
}

/*
***************************************************************************
** Used if caller specifies a NULL print function on print.
***************************************************************************
*/ 
static void printPointer(void* pointer)
{
    irxError("0x%lX\n", (long)pointer);
}

/*
***************************************************************************
** hash function for strings
***************************************************************************
*/
static long hashString(char* key, long numBuckets)
{
    long hash = 0;

    if (key == NULL) return 0;

    while (*key) 
    {
        /* Extract the significant ascii bits from the current character
         * and shift the into the hash key */
        hash = (hash^(((*key)&0x001f) | (((*key)&0x0040)>>1)))<<5;
        key++;
    }
    hash = hash>>5;

    return (unsigned)hash % numBuckets;
}

/*
***************************************************************************
** copy function for strings - we use irxMemAlloc for the usual reasons
***************************************************************************
*/
static char* hashCopyString (char* key)
{
    if (key == NULL) return NULL;
    return strcpy(IRX_NEW_ARRAY(char,strlen(key)+1), key);
}


/*
***************************************************************************
** Hash function for long keys
***************************************************************************
*/
static long hashLong (void* key, long numBuckets)
{
    return (long)key % numBuckets;
}

/*
***************************************************************************
** Compare function for long keys
***************************************************************************
*/
static int hashCompareLong (void* key1, void* key2)
{
    return (long)key1 - (long)key2;
}

/*
***************************************************************************
** Copy function for long keys
***************************************************************************
*/
static void* hashCopyLong (void* key)
{
    return key;
}

/*
***************************************************************************
** Free function for long keys
***************************************************************************
*/
static void hashFreeLong(void* key) { key=key;/* EMPTY */ }

/*f
***************************************************************************
** Convenience function for printing string keys and data.
***************************************************************************
*/
void irxHashPrintString(void* str)
{
    irxError("%-20s", (char*)str);
}

/*f
***************************************************************************
** Convenience function for printing double data values.
***************************************************************************
*/
void irxHashPrintDouble(void* data)
{
    irxError("%f", *(double*)data);
}




