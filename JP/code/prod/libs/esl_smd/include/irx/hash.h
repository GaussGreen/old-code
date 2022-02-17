/*
***************************************************************************
** HEADER FILE: hash.h
**
** Defines a generic hash table with an arbitrary key type and data type,
** plus some concrete implementations where the hash key is a string or
** an integer (long) - although the data type is still arbitrary.
**
** $Header$
***************************************************************************
*/
#ifndef IRX_HASH_H
#define IRX_HASH_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Data typedefs */

typedef struct _IrxTHashTable IrxTHashTable;

/* Function typedefs */

/* Hash key comparator function pointer for comparing hash keys */
typedef int (IrxTHashKeyCmpFunc) (void* pvKey1, void* pvKey2);

/* Hash key duplicator function pointer for duplicating hash keys */
typedef void* (IrxTHashKeyCopyFunc) (void* pvKey);

/* Hash key free function pointer for freeing keys */
typedef void (IrxTHashKeyFreeFunc) (void* pvKey);

/* Hashing function pointer */
typedef long (IrxTHashFunc) (void* pKey, long numBuckets);

/* Hash key print function pointer */
typedef void (IrxTHashKeyPrintFunc) (void* pKey);

/* Hash table data print function pointer */
typedef void (IrxTHashDataPrintFunc) (void* pData);

/* Hash table data free function pointer */
typedef void (IrxTHashDataFreeFunc) (void* pData);

/*
***************************************************************************
** To save an excessive amount of casting, we provide a number of macros
** for manipulating a hash table.
***************************************************************************
*/
#define IRX_HASH_INSERT(ht, key, data) \
irxHashTableInsert((ht), (void*)(key), (void*)(data))

#define IRX_HASH_SEARCH(ht, key, address_of_data) \
irxHashTableSearch((ht), (void*)(key), (void**)(address_of_data))

#define IRX_HASH_DELETE_ENTRY(ht, key) \
irxHashTableDeleteEntry((ht), (void*)(key))

#define IRX_HASH_EXTRACT(ht, key, address_of_data) \
irxHashTableExtract((ht), (void*)(key), (void**)(address_of_data))

#define IRX_HASH_FIRST(ht, address_of_key, address_of_data) \
irxHashTableFirst((ht), (void**)(address_of_key), (void**)(address_of_data))

#define IRX_HASH_NEXT(ht, address_of_key, address_of_data) \
irxHashTableNext((ht), (void**)(address_of_key), (void**)(address_of_data))

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
);

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
);

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
** compare keys, hash keys, copy keys and free keys.
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
(long                  numBuckets,   /* (I) Length of the table */
 IrxTBool              uniqueKey,    /* (I) Keys need to be unique */
 IrxTHashFunc*         hashFunc,     /* (I) Hash function */
 IrxTHashKeyCmpFunc*   hashKeyCmp,   /* (I) Hash key compare func */
 IrxTHashKeyCopyFunc*  hashKeyCopy,  /* (I) Hash key copy func */
 IrxTHashKeyFreeFunc*  hashKeyFree,  /* (I) Hash key free function */
 IrxTHashDataFreeFunc* hashDataFree  /* (I) Hash data free function 
                                      (NULL=don't free)*/
);

/*f
***************************************************************************
** Inserts a new key,value pair into the hash table.
**
** The data is shallow copied. The key is copied using the hashKeyCopy
** function.
**
** If keys need to be unique, then it is an error to insert the same key
** twice. If keys do not need to be unique, then inserting the same key
** twice will put duplicate keys into the hash table.
***************************************************************************
*/
int  irxHashTableInsert
(IrxTHashTable* hashTable,    /* (I) Hash table                */
 void*          key,          /* (I) Key (to be copied)        */
 void*          data          /* (I) Data (to be stored as is) */
);

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
(IrxTHashTable* hashTable,             /* (I) Hash table identifier */
 void*        key,                   /* (I) Key */
 void**       data);                 /* (O) Pointer to data */

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
 void*        key);                    /* (I) Key */

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
 void**       data);                /* (O) Pointer to data */


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
);

/*f
***************************************************************************
** Deallocates the hash table, including all keys and data.
**
** Note that the data is only free'd if the user provided a free function
** when constructing the hash table.
***************************************************************************
*/
void irxHashTableFree (IrxTHashTable* hashTable);

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
** The key copy function is not called when setting the output of the key.
** In fact the user should treat the key and data returned with great
** reverence - they are effectively const although not declared as such.
**
** Note that irxHashTableFirst and irxHashTableNext are not a thread-safe
** combination - you should not iterate through the same hash table in
** different threads.
***************************************************************************
*/
IrxTBool irxHashTableFirst
(IrxTHashTable* hashTable,             /* (I) Hash table Id */
 void**       key,                   /* (O) Key for next elem */
 void**       data);                 /* (O) Data for 1st elem */


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
 void**       data                   /* (O) Data for next elem */
);

/*f
***************************************************************************
** Convenience function for printing string keys and data.
***************************************************************************
*/
void irxHashPrintString(void* str);

/*f
***************************************************************************
** Convenience function for printing double data values.
***************************************************************************
*/
void irxHashPrintDouble(void* data);

#ifdef __cplusplus
}
#endif

#endif    /* IRX_HASH_H */

