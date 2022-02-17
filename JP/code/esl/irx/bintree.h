/*m -----------------------------------------------------------------
   SOURCE FILE: bintree.h

   PURPOSE:     Public functions create and search a binary tree.

   CREATED BY:  Lawrence Miller (25 August 1998)
   MODIFIED BY: Charles Irwin (22 March 1999) - added red-black tree
                algorithm to maintain balance in the tree. Also added
                deleting capabilities.

 --------------------------------------------------------------------*/
#include "cgeneral.h"

#ifndef IRX_BINTREE_H
#define IRX_BINTREE_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct _IrxTBinaryTree IrxTBinaryTree;

typedef void (IrxTBinaryTreePrintFunc) (
    void * key,      /* (I) Key */
    void * data);    /* (I) Data */

typedef int (IrxTBinaryTreeCompareFunc) (
    void * k1,   /* (I) Key 1 */
    void * k2);  /* (I) Key 2 */

typedef void (IrxTBinaryTreeFreeFunc) (
    void * p);  /* (I) Pointer to data structure */

/* returns SUCCESS or FAILURE */
typedef int (IrxTBinaryTreeCopyFunc) (
    void  * iKey,   /* (I) Pointer to key structure */
    void  * iData,  /* (I) Pointer to data structure */
    void ** oKey,   /* (O) Pointer to key structure */
    void ** oData); /* (O) Pointer to data structure */

typedef int (IrxTBinaryTreeForEachFunc) (
    IrxTBinaryTree *T, /* (I) Red-Black tree being recursed */
    void        *client_data,   /* (I/O) data the application registered */
    void        *call_data);    /* (I/O) callback specific data */

void irxBinaryTreePrint(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to print from */
    IrxTBinaryTreePrintFunc * PRINT);      /* (I) Printing function */

void irxBinaryTreePrintNoIndent(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to print from */
    IrxTBinaryTreePrintFunc * PRINT);      /* (I) Printing function */

/*
 * Sequentially calls function fn with each node's data in tree T.
 */
int irxBinaryTreeForEach
(IrxTBinaryTree            *T,     /* (I/O) Red-Black tree to print from */
 IrxTBinaryTreeForEachFunc *fn,    /* (I/O) function to apply to each node
                                 * call data will be a pointer to node's key */
 void             *client_data  /* (I/O) passed to each call of fn */
);

/*
 * Sequentially calls function fn with each node's key in tree T.
 */
int irxBinaryTreeForEachKey
(IrxTBinaryTree            *T,     /* (I/O) Red-Black tree to print from */
 IrxTBinaryTreeForEachFunc *fn,    /* (I/O) function to apply to each node
                                 * call data will be a pointer to node's key */
 void             *client_data  /* (I/O) passed to each call of fn */
);

int irxBinaryTreeDelete(
    IrxTBinaryTree            *T,             /* (I) Tree to delete from */
    void                   *key,           /* (I) Key to match */
    IrxTBinaryTreeCompareFunc *COMPARE,       /* (I) Comparison function */
    void                  **returnedKey,   /* (O) Key found */
    void                  **returnedData); /* (O) Data found. */

int irxBinaryTreeInsert(
    IrxTBinaryTree *T,                       /* (I) Tree to insert into */
    void    *key,                         /* (I) Key to insert */
    void    *data,                        /* (I) Data to insert */
    IrxTBinaryTreeCompareFunc * COMPARE);    /* (I) Comparison function */

void * irxBinaryTreeSearch(
    IrxTBinaryTree *T,                         /* (I) Tree to search */
    void    *key,                           /* (I) Key to match */
    IrxTBinaryTreeCompareFunc * COMPARE);      /* (I) Comparison function */

void * irxBinaryTreeMinimum(
    IrxTBinaryTree *T);                  /* (I) Red-Black tree to find minimum */

int irxBinaryTreeMaximum(
    IrxTBinaryTree  *T,        /* (I) Red-Black tree to find maximum */
    void        **data);    /* (O) maximum data */

void irxBinaryTreeFree(
    IrxTBinaryTree *T,                      /* (I) Red-Black tree to deallocate */
    IrxTBinaryTreeFreeFunc *FREE_DATA,      /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY);      /* (I) Key free function */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeCopy
   Created by: Lawrence Miller 5/18/99
   Purpose: Copy a binary tree. (deep)
            The key and data function pointers are used to clean up if the
            copy fails midway through creation of the new tree. The free
            functions are treated in the same way as in irxBinaryTreeFree.
   -------------------------------------------------------------------------- */
IrxTBinaryTree * irxBinaryTreeCopy(
    IrxTBinaryTree           *T,             /* (I) Red-Black tree */
    IrxTBinaryTreeCopyFunc   *copyDataKey,   /* (I) Copy data & key */
    IrxTBinaryTreeFreeFunc   *freeData,      /* (I) Data free function */
    IrxTBinaryTreeFreeFunc   *freeKey);      /* (I) Key free function */

IrxTBinaryTree * irxBinaryTreeNew(void);

/*
 * List strings that are used for keys in a binary tree.
 * The array and its contents that is returned is owned by the caller.
 */
int irxBinaryTreeListKeys
(IrxTBinaryTree    *T,         /* (I) Red-Black tree to list keys from */
 char         ***listOfKeys,/* (O) The list of strings */
 int            *numKeys    /* (O) Number of keys in list */
);

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeInsertStr
   Created by: Lawrence Miller 8/25/98
   Purpose: Insert data into the tree using a null terminated string as the
            key, and interpreting all keys in the tree as null terminated
            strings.

            The tree is ordered based on strcmp of the new key and the key
            associated with each node.  If two identical keys are passed
            into the tree, both will be stored.

            A given data block should not be inserted into the tree more than
            once, unless you use the irxBinaryTreeFreeNodes routine to destroy
            the tree, and handle freeing the data blocks yourself. (the
            irxBinaryTreeFreeData routine frees both data and nodes, eliminating
            the need to keep track of the data pointers seperately)

            Returns FAILURE if unable to insert the key
   -------------------------------------------------------------------------- */
int irxBinaryTreeInsertStr(
                  char * key,             /* (I) pointer to key string */
                  void * data,            /* (I) pointer to data */
                  IrxTBinaryTree ** root);   /* (I/O) root of tree (pointer to
                                                   NULL for empty tree) */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeSearchStr
   Created by: Lawrence Miller 8/25/98
   Purpose: Insert data into the tree using a null terminated string as the
            key, and interpreting all keys in the tree as null terminated
            strings.

            The tree is ordered based on strcmp of the new key and the key
            associated with each node.  If two identical keys are passed
            into the tree, both will be stored.

            Returns pointer to data, or NULL if key was not found
   -------------------------------------------------------------------------- */
void * irxBinaryTreeSearchStr(
                  char * key,             /* (I) pointer to key string */
                  IrxTBinaryTree * root);   /* (I) root of tree */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFreeData
   Created by: Lawrence Miller 8/25/98
   Purpose: Free a binary tree.  Frees tree nodes as well as calling FREE on
            all data pointers inserted into the tree.  Does not call FREE
            on keys.

            Returns
   -------------------------------------------------------------------------- */
int irxBinaryTreeFreeData(
                  IrxTBinaryTree * root);   /* (I) root of tree */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFreeNodes
   Created by: Lawrence Miller 8/25/98
   Purpose: Free a binary tree.  Frees tree nodes only, not data nodes

            Returns
   -------------------------------------------------------------------------- */
int irxBinaryTreeFreeNodes(
                  IrxTBinaryTree * root);   /* (I) root of tree */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFreeNodes
   Created by: Lawrence Miller 8/25/98
   Purpose: Free a binary tree.  Frees tree nodes only, not data nodes

            Returns
   -------------------------------------------------------------------------- */
int irxBinaryTreeFreeNodes(
                  IrxTBinaryTree * root);   /* (I) root of tree */

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeGetLength
   Created by: Noel Yap
   Purpose: Returns the number of nodes inside a binary tree
   -------------------------------------------------------------------------- */
int irxBinaryTreeGetLength(
    IrxTBinaryTree *root,  /* (I) node of tree */
    size_t *length);    /* (O) number of nodes in the tree */

#ifdef __cplusplus
}
#endif


#endif
