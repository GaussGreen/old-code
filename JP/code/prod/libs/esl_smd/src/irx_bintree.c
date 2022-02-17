/*m -----------------------------------------------------------------
   SOURCE FILE: bintree.c

   PURPOSE:     Functions to create and search a binary tree.

   CREATED BY:  Lawrence Miller (25 August 1998)
   MODIFIED BY: Charles Irwin (22 March 1999) - added red-black tree
                algorithm to maintain balance in the tree. Also added
                deleting capabilities.

   NOTE: The algorithm is derived from "Introduction to Algorithms", by
   Thomas Cormen, Charles Leiserson, and Ronald Rivest.

   $Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/src/bintree.c,v 1.36 2004/03/02 14:57:19 cbauer Exp $
 --------------------------------------------------------------------*/

#include "irx/bintree.h"

#include <string.h>

#include "irx/macros.h"         /* NEW */
#include "irx/error.h"          /* IrxErrMsg */
#include "irx/memutils.h"
#include "irx/strutils.h"

typedef enum _TColorType {BLACK, RED} TColorType;

typedef struct _TTreeNode
{
    void              *key;
    void              *data;
    TColorType         color;
    struct _TTreeNode *parent;
    struct _TTreeNode *left;
    struct _TTreeNode *right;
} TTreeNode;

struct _IrxTBinaryTree
{
    TTreeNode *root;
    TTreeNode *nil;
    size_t     length;
};

typedef struct _TBinTreeKeysBuffer
{
    char  **list;
    int     maxKeys;
    int     index;
} TBinTreeKeysBuffer;

/*
 * These macros are to be used in the source code, to retrieve or set the
 * properties of tree nodes.
 */
#define ROOT(T)   (T->root)
#define NIL(T)    (T->nil)
#define LENGTH(T) (T->length)

#define COLOR(x)  (x->color)
#define LEFT(x)   (x->left)
#define RIGHT(x)  (x->right)
#define PARENT(x) (x->parent)
#define DATA(x)   (x->data)
#define KEY(x)    (x->key)

static TTreeNode * newEmptyRBNode(void);

static TTreeNode * newRBNode(
    IrxTBinaryTree *T,                     /* (I) Tree to create new node in */
    void    *key,                       /* (I) Key to fill into node */
    void    *data);                     /* (I) Data to fill into node */

static void freeRBNode(
    TTreeNode           *x,           /* (I) RB Node to free */
    IrxTBinaryTreeFreeFunc *FREE_DATA,   /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY);   /* (I) Key free function */

static void freeRBSubTree(
    IrxTBinaryTree   *T,                 /* (I) x is the subtree of T to free. */
    TTreeNode *x,                     /* (I) Subtree of T to deallocate.*/
    IrxTBinaryTreeFreeFunc *FREE_DATA,   /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY);   /* (I) Key free function */

static TTreeNode * rbTreeSearch(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to search */
    TTreeNode *x,                       /* (I) Subtree on which to search */
    void *key,                          /* (I) Key to match */
    IrxTBinaryTreeCompareFunc *COMPARE);   /* (I) Comparison function for keys */

static TTreeNode * rbTreeMinimum(
    IrxTBinaryTree *T,                           /* (I) Red-Black tree */
    TTreeNode   *x);                          /* (I) Subtree to find minimum */

static TTreeNode * rbTreeMaximum(
    IrxTBinaryTree *T,                           /* (I) Red-Black tree */
    TTreeNode   *x);                          /* (I) Subtree to find maximum */

static TTreeNode * rbTreeSuccessor(
    IrxTBinaryTree   *T,                             /* Red-Black tree */
    TTreeNode *x);                            /* Node to find successor of */

static int rbTreeInsert(
    IrxTBinaryTree            *T,          /* (I) Tree to insert into */
    TTreeNode              *z,          /* (I) Node to insert into tree */
    IrxTBinaryTreeCompareFunc *COMPARE);   /* (I) Comparison function for keys */

static int rbLeftRotate(
    IrxTBinaryTree *T,                           /* (I) Tree to rotate */
    TTreeNode   *x);                          /* (I) Node to rotate at */

static int rbRightRotate(
    IrxTBinaryTree *T,                           /* (I) Tree to rotate */
    TTreeNode   *x);                          /* (I) Node to rotate at */

static int rbInsert(
    IrxTBinaryTree            *T,          /* (I) Tree to insert into */
    TTreeNode              *x,          /* (I) Node to insert */
    IrxTBinaryTreeCompareFunc *COMPARE);   /* (I) Comparison function for keys */

static int rbDeleteFixup(
    IrxTBinaryTree *T,                    /* (I) Red-Black tree to fix up */
    TTreeNode   *x);                   /* (I) Node to start fixing up at */

static TTreeNode * rbDelete(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to delete from */
    TTreeNode   *z);                    /* (I) Node to delete */

static void rbTreePrint(
    IrxTBinaryTree          *T,          /* (I) Red-Black tree to print from */
    TTreeNode            *x,          /* (I) Subtree to print */
    IrxTBinaryTreePrintFunc *PRINT,      /* (I) Printing function */
    long                  level,      /* (I) Level of subtree */
    IrxTBool              indent);    /* (I) Should indenting be done? */

static int rbTreeForEach(
    IrxTBinaryTree             *T,            /* (I/O) Red-Black tree to print from */
    TTreeNode               *x,            /* (I/O) Subtree to print */
    IrxTBinaryTreeForEachFunc  *fn,           /* (I/O) generate function */
    void                    *client_data); /* (I/O) client data passed into fn */
static int rbTreeForEachKey(
    IrxTBinaryTree             *T,            /* (I/O) Red-Black tree to print from */
    TTreeNode               *x,            /* (I/O) Subtree to print */
    IrxTBinaryTreeForEachFunc  *fn,           /* (I/O) generate function */
    void                    *client_data); /* (I/O) client data passed into fn */

static TTreeNode * copyRBNode(
    TTreeNode           * in,            /* (I) Source node */
    IrxTBinaryTreeCopyFunc * copyDataKey);  /* (I) Function to copy cata & key */

static int rbCopySubTree(
    IrxTBinaryTree         * inT,          /* (I) Source tree */
    IrxTBinaryTree         * outT,         /* (I) Target tree */
    IrxTBinaryTreeCopyFunc * copyDataKey,  /* (I) Key and data copy function */
    IrxTBinaryTreeFreeFunc * freeData,     /* (I) Function to free data */
    IrxTBinaryTreeFreeFunc * freeKey,      /* (I) Function to free key */
    TTreeNode           * in,           /* (I) Source node */
    TTreeNode           * outParent,    /* (I) Parent of target node */
    TTreeNode          ** newOut);      /* (O) Target node */

IRX_DECLARE_AUTO_MEMCACHE(treeNode, TTreeNode)

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeDelete
   Created by: Charles Irwin 3/15/1999
   Purpose: Deletes a node which matches key.
   The key and data found are outputs. The user is responsible for
   freeing the returned key and data.
   --------------------------------------------------------------------------*/

int irxBinaryTreeDelete(
    IrxTBinaryTree            *T,             /* (I) Tree to delete from */
    void                   *key,           /* (I) Key to match */
    IrxTBinaryTreeCompareFunc *COMPARE,       /* (I) Comparison function */
    void                  **returnedKey,   /* (O) Key found */
    void                  **returnedData)  /* (O) Data found. */
{
    static char routine[] = "irxBinaryTreeDelete";
    int            status = FAILURE;

    TTreeNode *foundNode = NULL;         /* Node found which matches key */
    TTreeNode *deleteNode = NULL;        /* Node to deallocate from tree */

    if ((T == NULL) ||
        (key == NULL) ||
        (COMPARE == NULL) ||
        (returnedKey == NULL) ||
        (returnedData == NULL))
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    foundNode = rbTreeSearch(T, ROOT(T), key, COMPARE);
    if (foundNode == NULL)
        goto RETURN;

    if (foundNode != NIL(T))
    {
        *returnedKey  = KEY(foundNode);
        *returnedData = DATA(foundNode);

        deleteNode = rbDelete(T, foundNode);
        if (deleteNode == NULL)
            goto RETURN;
    }
    else
    {
        *returnedKey  = NULL;
        *returnedData = NULL;
    }

    status = SUCCESS;

RETURN:
    treeNodeFree(deleteNode);
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeInsert
   Created by: Charles Irwin 3/15/1999
   Purpose: Inserts a node containing a key and data into a binary tree.
   --------------------------------------------------------------------------*/

int irxBinaryTreeInsert(
    IrxTBinaryTree            *T,                 /* (I) Tree to insert into */
    void                   *key,               /* (I) Key to insert */
    void                   *data,              /* (I) Data to insert */
    IrxTBinaryTreeCompareFunc *COMPARE)           /* (I) Comparison function */
{
    static char routine[] = "irxBinaryTreeInsert";
    int            status = FAILURE;

    TTreeNode *newNode = NULL;            /* New node to insert */

    if ((T == NULL) ||
        (key == NULL) ||
        (data == NULL) ||
        (COMPARE == NULL))
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    newNode = newRBNode(T, key, data);
    if (newNode == NULL)
        goto RETURN;

    if (rbInsert(T, newNode, COMPARE) != SUCCESS)
        goto RETURN;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        treeNodeFree(newNode);
        irxErrorFailure(routine);
    }
    return (status);
}


/* --------------------------------------------------------------------------
   Function: irxBinaryTreeSearch
   Created by: Charles Irwin 3/15/1999
   Purpose: Searches a tree for a node which matches key. Upon finding a node
            returns the data associated with it. Returns NULL, if key is not
            found.
   --------------------------------------------------------------------------*/

void * irxBinaryTreeSearch(
    IrxTBinaryTree            *T,             /* (I) Tree to search */
    void                   *key,           /* (I) Key to match */
    IrxTBinaryTreeCompareFunc *COMPARE)       /* (I) Comparison function */
{
    static char routine[] = "irxBinaryTreeSearch";
    int            status = FAILURE;

    TTreeNode *foundNode = NULL;            /* Found node */
    void      *foundData = NULL;            /* Found data */

    if ((T == NULL) ||
        (key == NULL) ||
        (COMPARE == NULL))
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    foundNode = rbTreeSearch(T, ROOT(T), key, COMPARE);
    if (foundNode == NULL)
        goto RETURN;

    /*
     * Note: Even if foundNode = NIL(T), the data in NIL(T) is NULL, so
     * this works out OK.
     */

    foundData = DATA(foundNode);
    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        foundData = NULL;
        irxErrorFailure(routine);
    }
    return (foundData);
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeMinimum
   Created by: Charles Irwin 3/15/1999
   Purpose: Finds the minimum of a red-black tree.
   --------------------------------------------------------------------------*/
void * irxBinaryTreeMinimum(
    IrxTBinaryTree *T)                    /* (I) Red-Black tree to find minimum */
{
    static char routine[] = "irxBinaryTreeMinimum";
    int            status = FAILURE;

    TTreeNode *foundNode = NULL;
    void      *foundData = NULL;

    if (T == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    foundNode = rbTreeMinimum(T, ROOT(T));
    if (foundNode == NULL)
        goto RETURN;

    foundData = DATA(foundNode);

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        foundData = NULL;
        irxErrorFailure(routine);
    }
    return (foundData);
}
/* --------------------------------------------------------------------------
   Function: irxBinaryTreeMaximum
   Created by: Charles Irwin 3/15/1999
   Purpose: Finds the maximum of a red-black tree.
   --------------------------------------------------------------------------*/
int irxBinaryTreeMaximum(
    IrxTBinaryTree  *T,        /* (I) Red-Black tree to find maximum */
    void        **data)     /* (O) maximum data */
{
    TTreeNode       *foundNode = NULL;
    static  char     routine[] = "irxBinaryTreeMaximum";
    int              status = FAILURE;

    *data = NULL;

    if (T == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    foundNode = rbTreeMaximum(T, ROOT(T));
    if (foundNode == NULL)
        goto RETURN;

    *data = DATA(foundNode);

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }

    return status;
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFree
   Created by: Charles Irwin 3/15/1999
   Purpose: Deallocates all the nodes in a red-black tree including the NIL(T)
            node.
   --------------------------------------------------------------------------*/

void irxBinaryTreeFree(
    IrxTBinaryTree *T,                      /* (I) Red-Black tree to deallocate */
    IrxTBinaryTreeFreeFunc *FREE_DATA,      /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY)       /* (I) Key free function */
{
    if (T != NULL)
    {
        freeRBSubTree(T, ROOT(T),FREE_DATA,FREE_KEY);
        treeNodeFree(NIL(T));
        FREE(T);
    }
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeNew
   Created by: Charles Irwin 3/15/1999
   Purpose: Creates a new empty Red-Black tree. The only node in this tree
            is the NIL(T) node, which is black. Returns NULL upon failure.
   --------------------------------------------------------------------------*/

IrxTBinaryTree * irxBinaryTreeNew(void)
{
    static char routine[] = "irxBinaryTreeNew";
    int            status = FAILURE;

    IrxTBinaryTree *newTree = NULL;         /* New empty tree */

    newTree = NEW(IrxTBinaryTree);
    if (newTree == NULL)
    {
        irxError("%s: Cannot construct new empty RB tree.\n", routine);
        goto RETURN;
    }

    NIL(newTree) = NULL;
    ROOT(newTree) = NULL;

    NIL(newTree) = newEmptyRBNode();
    if (NIL(newTree) == NULL)
    {
        irxError("%s: Cannot construct new empty NIL node.\n", routine);
        goto RETURN;
    }

    ROOT(newTree) = NIL(newTree);
    PARENT(NIL(newTree)) = ROOT(newTree);
    LENGTH(newTree) = 0;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxBinaryTreeFree(newTree,NULL,NULL);
        newTree = NULL;
        irxErrorFailure(routine);
    }
    return (newTree);
}

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
    IrxTBinaryTreeFreeFunc   *freeKey)       /* (I) Key free function */
{
    static char routine[] = "irxBinaryTreeCopy";
    int            status = FAILURE;
    IrxTBinaryTree  *newTree = NULL;

    if ((T == NULL) ||
        (copyDataKey == NULL))
    {
        irxError("%s: Cannot process NULL pointers.\n", routine);
        goto RETURN;
    }

    newTree = irxBinaryTreeNew();

    if (newTree == NULL)
        goto RETURN;

    /*
     * Note: If there is failure here, ROOT(newTree) does NOT get set.
     */
    if (rbCopySubTree(T,
                      newTree,
                      copyDataKey,
                      freeData,
                      freeKey,
                      ROOT(T),
                      NIL(newTree),
                      &(ROOT(newTree))) != SUCCESS)
    {
        goto RETURN;
    }

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxBinaryTreeFree(newTree, freeData, freeKey);
        newTree = NULL;
        irxErrorFailure(routine);
    }
    return (newTree);
}

/* --------------------------------------------------------------------------
   Function: irxRBTreePrint
   Created by: Charles Irwin 3/15/1999
   Purpose: Prints the contents of a binary tree.
   --------------------------------------------------------------------------*/

void irxBinaryTreePrint(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to print from */
    IrxTBinaryTreePrintFunc * PRINT)       /* (I) Printing function */
{
    if (T != NULL)
    {
        rbTreePrint(T, ROOT(T), PRINT, 0, TRUE);
    }
}

/* --------------------------------------------------------------------------
   Function: irxRBTreePrint
   Created by: Lawrence Miller 5/12/1999
   Purpose: Prints the contents of a binary tree, but without auto-indenting.
   --------------------------------------------------------------------------*/

void irxBinaryTreePrintNoIndent(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to print from */
    IrxTBinaryTreePrintFunc * PRINT)       /* (I) Printing function */
{
    if (T != NULL)
    {
        rbTreePrint(T, ROOT(T), PRINT, 0, FALSE);
    }
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeForEach
   Created by: Noel Yap
   Purpose: Sequentially calls function fn with each node's data in tree T.
   --------------------------------------------------------------------------*/
int irxBinaryTreeForEach(
    IrxTBinaryTree             *T,             /* (I/O) Red-Black tree to print from */
    IrxTBinaryTreeForEachFunc  *fn,            /* (I/O) function
                                             * call data will be a pointer to node's data */
    void                    *client_data)   /* (I/O) passed to each call of fn */
{
    int status = FAILURE;

    if (T != NULL)
    {
        if (rbTreeForEach(T, ROOT(T), fn, client_data) != SUCCESS)
        {
            goto RETURN;
        }
    }

    status = SUCCESS;

RETURN:
    return status;
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeForEachKey
   Created by: Noel Yap
   Purpose: Sequentially calls function fn with each node's key in tree T.
   --------------------------------------------------------------------------*/
int irxBinaryTreeForEachKey
(IrxTBinaryTree            *T,   /* (I/O) Red-Black tree to print from */
 IrxTBinaryTreeForEachFunc *fn,  /* (I/O) function
                               * call data will be a pointer to node's key */
 void         *client_data)   /* (I/O) passed to each call of fn */
{
    int status = FAILURE;

    if (T != NULL)
    {
        if (rbTreeForEachKey(T, ROOT(T), fn, client_data) != SUCCESS)
        {
            goto RETURN;
        }
    }

    status = SUCCESS;

RETURN:
    return status;
}

static int AppendBinTreeKeyToBuffer
(IrxTBinaryTree *T,
 void        *client_data,
 void        *call_data)
{
    static char routine[] = "AppendBinTreeKeyToBuffer";
    int         status = FAILURE;

    TBinTreeKeysBuffer *buf = (TBinTreeKeysBuffer*) client_data; 
    char               *key = (char*) call_data;

    if (buf->index >= buf->maxKeys)
    {
        irxError("%s: Binary tree keys buffer overflow.\n", routine);
        goto RETURN; /* failure */
    }

    buf->list[buf->index] = STRDUP(key);
    if (buf->list[buf->index] == NULL)
    {
        goto RETURN; /* failure */
    }
    buf->index++;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return(status);
}


/* --------------------------------------------------------------------------
   Function: irxBinaryTreeListKeys
   Created by: Chris Bauer
   Purpose: List strings that are used for keys in a binary tree.
            The array and its contents that is returned is owned by the caller.
   --------------------------------------------------------------------------*/
int irxBinaryTreeListKeys
(IrxTBinaryTree    *T,         /* (I) Red-Black tree to list keys from */
 char         ***listOfKeys,/* (O) The list of strings */
 int            *numKeys    /* (O) Number of keys in list */
)
{
    static char routine[] = "irxBinaryTreeListKeys";
    int         status = FAILURE;
    char      **list = NULL;
    size_t      length;

    TBinTreeKeysBuffer buf;

    REQUIRE(listOfKeys != NULL);
    REQUIRE(numKeys != NULL);

    if (irxBinaryTreeGetLength(T, &length) != SUCCESS)
    {
        goto RETURN; /* failure */
    }
    list = NEW_ARRAY(char*, length);
    if (list == NULL)
    {
        goto RETURN; /* failure */
    }

    buf.list = list;
    buf.maxKeys = (int)length;
    buf.index = 0;

    if (irxBinaryTreeForEachKey(T, 
                         (IrxTBinaryTreeForEachFunc*) AppendBinTreeKeyToBuffer,
                         (void*)&buf) != SUCCESS)
    {
        goto RETURN; /* failure */
    }

    *listOfKeys = list;
    *numKeys = (int) length;
            
    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        FREE_PTR_ARRAY(list, (int)length, irxMemFree);
        irxErrorFailure(routine);
    }
    return(status);
}

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
            irxBinaryTreeFreeData routine frees both data and nodes,
            eliminating the need to keep track of the data pointers
            separately)

            Returns FAILURE if unable to insert the key
   ------------------------------------------------------------------------- */
int irxBinaryTreeInsertStr(
                  char * key,             /* (I) pointer to key string */
                  void * data,            /* (I) pointer to data */
                  IrxTBinaryTree ** root)    /* (I/O) root of tree (pointer to
                                              NULL for empty tree) */
{
    static char routine[] = "irxBinaryTreeInsertStr";
    int         status    = FAILURE;

    int neededNewTree = FALSE;    /* Whether we allocated a new tree */

    if (root == NULL)
    {
        irxError("%s: Cannot process NULL inputs.\n", routine);
        goto RETURN;
    }

    if (*root == NULL)
    {
        /*
         * Allocate a new empty binary tree.
         */
        *root = irxBinaryTreeNew();
        if (*root == NULL)
        {
            irxError("%s: Unable to create new tree\n",routine);
            goto RETURN;
        }
        neededNewTree = TRUE;   /* We needed to allocate a new tree */
    }

    if (irxBinaryTreeInsert(*root,
                            key,
                            data,
                            (IrxTBinaryTreeCompareFunc*) strcmp) != SUCCESS)
    {
        irxError("%s: Unable to insert key \"%s\" into tree.\n",
                  routine,
                  key);
        goto RETURN;
    }

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        if (neededNewTree)
        {
            /*
             * If we needed to allocate a new tree, and for some reason
             * couldn't insert afterwards, then we should clean up
             * the memory used in allocated a new empty tree.
             */
            irxBinaryTreeFree(*root, NULL, NULL);
            *root = NULL;
        }
        irxError("%s: Failed.\n", routine);
    }

    return (status);

}

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
                  IrxTBinaryTree * root)     /* (I) root of tree */
{
    void *retVal = NULL;

    retVal = irxBinaryTreeSearch(root,
                                 key,
                                 (IrxTBinaryTreeCompareFunc*)strcmp);
    return (retVal);
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFreeData
   Created by: Lawrence Miller 8/25/98
   Purpose: Free a binary tree.  Frees tree nodes as well as calling FREE on
            all data pointers inserted into the tree.  Does not call FREE
            on keys.

            Returns SUCCESS.
   -------------------------------------------------------------------------- */
int irxBinaryTreeFreeData(
                  IrxTBinaryTree * root)    /* (I) root of tree */
{
    irxBinaryTreeFree(root, irxMemFree, NULL);
    return (SUCCESS);
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeFreeNodes
   Created by: Lawrence Miller 8/25/98
   Purpose: Free a binary tree.  Frees tree nodes only, not data nodes

            Returns SUCCESS.
   -------------------------------------------------------------------------- */
int irxBinaryTreeFreeNodes(
                  IrxTBinaryTree * root)    /* (I) root of tree */
{
    irxBinaryTreeFree(root, NULL, NULL);
    return (SUCCESS);
}

/* --------------------------------------------------------------------------
   Function: newEmptyRBNode
   Created by: Charles Irwin 3/15/1999
   Purpose: Creates a new empty black treenode. Returns NULL on failure.
   --------------------------------------------------------------------------*/

static TTreeNode * newEmptyRBNode(void)
{
    static char routine[] = "newEmptyRBNode";
    int            status = FAILURE;

    TTreeNode *newNode = NULL;

    newNode = treeNodeNew();
    if (newNode == NULL)
    {
        irxError("%s: Cannot construct new empty RB node.\n", routine);
        goto RETURN;
    }

    KEY(newNode)    = NULL;
    DATA(newNode)   = NULL;
    COLOR(newNode)  = BLACK;
    PARENT(newNode) = NULL;
    LEFT(newNode)   = NULL;
    RIGHT(newNode)  = NULL;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        treeNodeFree(newNode);
        newNode = NULL;
        irxErrorFailure(routine);
    }
    return (newNode);
}

/* --------------------------------------------------------------------------
   Function: newRBNode
   Created by: Charles Irwin 3/15/1999
   Purpose: Creates a new black treenode with key and data fields filled in.
            The left and right pointers are set to NIL(T). Returns NULL
            upon failure.
   --------------------------------------------------------------------------*/

static TTreeNode * newRBNode(
    IrxTBinaryTree *T,                         /* (I) Tree to create new node in */
    void    *key,                       /* (I) Key to fill into node */
    void    *data)                      /* (I) Data to fill into node */
{
    static char routine[] = "newRBNode";
    int            status = FAILURE;

    TTreeNode *newNode = NULL;

    newNode = newEmptyRBNode();
    if (newNode == NULL)
    {
        irxError("%s: Cannot create new empty RB node.\n", routine);
        goto RETURN;
    }

    LEFT(newNode)  = NIL(T);
    RIGHT(newNode) = NIL(T);
    DATA(newNode)  = data;
    KEY(newNode)   = key;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        treeNodeFree(newNode);
        newNode = NULL;
        irxErrorFailure(routine);
    }
    return(newNode);
}

/* --------------------------------------------------------------------------
   Function: freeRBNode
   Created by: Charles Irwin 7/8/1999
   Purpose: Deallocates a node in a tree.
   --------------------------------------------------------------------------*/

static void freeRBNode(
    TTreeNode           *x,           /* (I) RB Node to free */
    IrxTBinaryTreeFreeFunc *FREE_DATA,   /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY)    /* (I) Key free function */
{
    if (x != NULL)
    {
        if (FREE_DATA != NULL)
        {
            FREE_DATA(x->data);
        }
        if (FREE_KEY != NULL)
        {
            FREE_KEY(x->key);
        }
        treeNodeFree(x);
    }
}

/* --------------------------------------------------------------------------
   Function: freeRBSubTree
   Created by: Charles Irwin 3/15/1999
   Purpose: Deallocates the nodes in a subtree.
   --------------------------------------------------------------------------*/

static void freeRBSubTree(
    IrxTBinaryTree         *T,           /* (I) x is the subtree of T to free. */
    TTreeNode           *x,           /* (I) Subtree of T to deallocate.*/
    IrxTBinaryTreeFreeFunc *FREE_DATA,   /* (I) Data free function */
    IrxTBinaryTreeFreeFunc *FREE_KEY)    /* (I) Key free function */
{
    if ((T != NULL) &&
        (x != NULL))
    {
        if (x != NIL(T))
        {
            /*
             * Take care not to free the NIL(T) node.
             */
            freeRBSubTree(T, LEFT(x),FREE_DATA,FREE_KEY);
            freeRBSubTree(T, RIGHT(x),FREE_DATA,FREE_KEY);
            freeRBNode(x, FREE_DATA, FREE_KEY);
        }
    }
}

/* --------------------------------------------------------------------------
   Function: rbTreeSearch
   Created by: Charles Irwin 3/15/1999
   Purpose: Searches a subtree of a red black tree for a node whose key
            field matches the key specified. COMPARE is a function to compare
            two keys. Returns NULL upon failure, the node if key is found, and
            NIL(T) if key is not found.
   --------------------------------------------------------------------------*/

static TTreeNode * rbTreeSearch(
    IrxTBinaryTree            *T,         /* (I) Red-Black tree to search */
    TTreeNode              *x,         /* (I) Subtree on which to search */
    void                   *key,       /* (I) Key to match */
    IrxTBinaryTreeCompareFunc *COMPARE)   /* (I) Comparison functions for keys */
{
    while ((x != NIL(T)) &&
           (COMPARE(key, KEY(x)) != 0))
    {
        if (COMPARE(key, KEY(x)) < 0)
        {
            x = LEFT(x);
        }
        else
        {
            x = RIGHT(x);
        }
    }

    return x;
}

/* --------------------------------------------------------------------------
   Function: rbTreeMinimum
   Created by: Charles Irwin 3/15/1999
   Purpose: Finds the minimum node in a binary search tree.
            Returns NULL upon failure, or NIL(T) if the subtree is empty.
   --------------------------------------------------------------------------*/

static TTreeNode * rbTreeMinimum(
    IrxTBinaryTree   *T,                       /* (I) Red-Black tree */
    TTreeNode     *x)                       /* (I) Subtree to find minimum */
{
    if (x != NIL(T))
    {
        while (LEFT(x) != NIL(T))
        {
            x = LEFT(x);
        }
    }

    return x;
}

/* --------------------------------------------------------------------------
   Function: rbTreeMaximum
   Created by: Charles Irwin 3/15/1999
   Purpose: Finds the maximum node in a binary search tree.
            Returns NULL upon failure, or NIL(T) if the subtree is empty.
   --------------------------------------------------------------------------*/

static TTreeNode * rbTreeMaximum(
    IrxTBinaryTree   *T,                        /* (I) Red-Black tree */
    TTreeNode     *x)                        /* (I) Subtree to find maximum */
{
    if (x != NIL(T))
    {
        while (RIGHT(x) != NIL(T))
        {
            x = RIGHT(x);
        }
    }

    return x;
}

/* --------------------------------------------------------------------------
   Function: rbTreeSuccessor
   Created by: Charles Irwin 3/15/1999
   Purpose: Finds the successor node in a binary search tree.
            Returns NULL upon failure, or if the subtree is empty. Returns
            NIL(T) if there is no successor to x.
   --------------------------------------------------------------------------*/

static TTreeNode * rbTreeSuccessor(
    IrxTBinaryTree   *T,                         /* Red-Black tree */
    TTreeNode     *x)                         /* Node to find successor of */
{
    static char routine[] = "rbTreeSuccessor";
    int            status = FAILURE;

    TTreeNode *y = NULL;

    if (x == NIL(T))
    {
        irxError("%s: Cannot find successor of a NIL node.\n", routine);
        goto RETURN;
    }

    if (RIGHT(x) != NIL(T))
    {
        y = rbTreeMinimum(T, RIGHT(x));
        if (y == NULL)
        {
            irxError("%s: Cannot find tree minimum.\n", routine);
            goto RETURN;
        }
    }
    else
    {
        y = PARENT(x);
        while ((y != NIL(T)) &&
               (x == RIGHT(y)))
        {
            x = y;
            y = PARENT(y);
        }
    }

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        y = NULL;
        irxErrorFailure(routine);
    }
    return y;
}

/* --------------------------------------------------------------------------
   Function: rbTreeInsert
   Created by: Charles Irwin 3/15/1999
   Purpose: Inserts a node into a binary search tree. Returns SUCCESS or
            FAILURE.
   --------------------------------------------------------------------------*/

static int rbTreeInsert(
    IrxTBinaryTree            *T,       /* (I) Tree to insert into */
    TTreeNode              *z,       /* (I) Node to insert into tree */
    IrxTBinaryTreeCompareFunc *COMPARE) /* (I) Comparison function for keys */
{
    static char routine[] = "rbTreeInsert";
    int            status = FAILURE;

    TTreeNode *x = NULL;           /* Node we are processing */
    TTreeNode *y = NULL;           /* Parent node of x */

    if (z == NIL(T))
    {
        irxError("%s: Cannot insert NIL node into tree.\n", routine);
        goto RETURN;
    }

    y = NIL(T);
    x = ROOT(T);

    while (x != NIL(T))
    {
        y = x;
        if (COMPARE(KEY(z), KEY(x)) < 0)
        {
            x = LEFT(x);
        }
        else
        {
            x = RIGHT(x);
        }
    }
    PARENT(z) = y;
    if (y == NIL(T))
    {
        ROOT(T) = z;
    }
    else
    {
        if (COMPARE(KEY(z), KEY(y)) < 0)
        {
            LEFT(y) = z;
        }
        else
        {
            RIGHT(y) = z;
        }
    }

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: rbLeftRotate
   Created by: Charles Irwin 3/15/1999
   Purpose: Performs a left rotation in the binary tree at a node.
   --------------------------------------------------------------------------*/

static int rbLeftRotate(
    IrxTBinaryTree *T,                           /* (I) Tree to rotate */
    TTreeNode   *x)                           /* (I) Node to rotate at */
{
    static char routine[] = "rbLeftRotate";
    int            status = FAILURE;

    TTreeNode *y = NULL;

    if (x == NIL(T))
    {
        irxError("%s: Cannot left rotate at a NIL node.\n", routine);
        goto RETURN;
    }

    y = RIGHT(x);
    RIGHT(x) = LEFT(y);

    if (LEFT(y) != NIL(T))
    {
        PARENT(LEFT(y)) = x;
    }

    PARENT(y) = PARENT(x);

    if (PARENT(x) == NIL(T))
    {
        ROOT(T) = y;
    }
    else
    {
        if (x == LEFT(PARENT(x)))
        {
            LEFT(PARENT(x)) = y;
        }
        else
        {
            RIGHT(PARENT(x)) = y;
        }
    }

    LEFT(y) = x;
    PARENT(x) = y;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: rbRightRotate
   Created by: Charles Irwin 3/15/1999
   Purpose: Performs a right rotation in the binary tree at a node.
   --------------------------------------------------------------------------*/

static int rbRightRotate(
    IrxTBinaryTree *T,                           /* (I) Tree to rotate */
    TTreeNode   *x)                           /* (I) Node to rotate at */
{
    static char routine[] = "rbRightRotate";
    int            status = FAILURE;

    TTreeNode *y = NULL;

    if (x == NIL(T))
    {
        irxError("%s: Cannot right rotate at a NIL node.\n", routine);
        goto RETURN;
    }

    y = LEFT(x);
    LEFT(x) = RIGHT(y);

    if (RIGHT(y) != NIL(T))
    {
        PARENT(RIGHT(y)) = x;
    }

    PARENT(y) = PARENT(x);

    if (PARENT(x) == NIL(T))
    {
        ROOT(T) = y;
    }
    else
    {
        if (x == RIGHT(PARENT(x)))
        {
            RIGHT(PARENT(x)) = y;
        }
        else
        {
            LEFT(PARENT(x)) = y;
        }
    }

    RIGHT(y) = x;
    PARENT(x) = y;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: rbInsert
   Created by: Charles Irwin 3/15/1999
   Purpose: Performs an insertion into a red-black tree. After inserting,
            rotations are performed to preserve balance in the tree.
   --------------------------------------------------------------------------*/

static int rbInsert(
    IrxTBinaryTree            *T,                /* (I) Tree to insert into */
    TTreeNode              *x,                /* (I) Node to insert */
    IrxTBinaryTreeCompareFunc *COMPARE)          /* (I) Comparision function */
{
    static char routine[] = "rbInsert";
    int            status = FAILURE;

    TTreeNode *y = NULL;

    if (x == NIL(T))
    {
        irxError("%s: Cannot insert NIL node into tree.\n", routine);
        goto RETURN;
    }

    if (rbTreeInsert(T, x, COMPARE) != SUCCESS)
    {
        irxError("%s: Cannot insert into binary tree.\n", routine);
        goto RETURN;
    }

    COLOR(x) = RED;
    while ((x != ROOT(T)) &&
           (COLOR(PARENT(x)) == RED))
    {
        if (PARENT(x) == LEFT(PARENT(PARENT(x))))
        {
            y = RIGHT(PARENT(PARENT(x)));
            if (COLOR(y) == RED)
            {
                COLOR(PARENT(x)) = BLACK;
                COLOR(y) = BLACK;
                COLOR(PARENT(PARENT(x))) = RED;
                x = PARENT(PARENT(x));
            }
            else
            {
                if (x == RIGHT(PARENT(x)))
                {
                    x = PARENT(x);
                    if (rbLeftRotate(T,x) != SUCCESS)
                    {
                        irxError("%s: Cannot left rotate tree.\n",
                                  routine);
                        goto RETURN;
                    }
                }
                COLOR(PARENT(x)) = BLACK;
                COLOR(PARENT(PARENT(x))) = RED;
                if (rbRightRotate(T,PARENT(PARENT(x))) != SUCCESS)
                {
                    irxError("%s: Cannot right rotate tree.\n",
                              routine);
                    goto RETURN;
                }
            }
        }
        else
        {
            y = LEFT(PARENT(PARENT(x)));
            if (COLOR(y) == RED)
            {
                COLOR(PARENT(x)) = BLACK;
                COLOR(y) = BLACK;
                COLOR(PARENT(PARENT(x))) = RED;
                x = PARENT(PARENT(x));
            }
            else
            {
                if (x == LEFT(PARENT(x)))
                {
                    x = PARENT(x);
                    if (rbRightRotate(T,x) != SUCCESS)
                    {
                        irxError("%s: Cannot left rotate tree.\n",
                                  routine);
                        goto RETURN;
                    }
                }
                COLOR(PARENT(x)) = BLACK;
                COLOR(PARENT(PARENT(x))) = RED;
                if (rbLeftRotate(T,PARENT(PARENT(x))) != SUCCESS)
                {
                    irxError("%s: Cannot right rotate tree.\n",
                              routine);
                    goto RETURN;
                }
            }
        }
    }
    COLOR(ROOT(T)) = BLACK;
    ++LENGTH(T);

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: rbDeleteFixup
   Created by: Charles Irwin 3/15/1999
   Purpose: Restores the balance in a binary tree after a delete operation.
   --------------------------------------------------------------------------*/

static int rbDeleteFixup(
    IrxTBinaryTree *T,                    /* (I) Red-Black tree to fix up */
    TTreeNode   *x)                    /* (I) Node to start fixing up at */
{
    static char routine[] = "rbDeleteFixup";
    int            status = FAILURE;

    TTreeNode *w = NULL;

    while ((x != ROOT(T)) &&
           (COLOR(x) == BLACK))
    {
        if (x == LEFT(PARENT(x)))
        {
            w = RIGHT(PARENT(x));
            if (COLOR(w) == RED)
            {
                COLOR(w) = BLACK;
                COLOR(PARENT(x)) = RED;
                if (rbLeftRotate(T, PARENT(x)) != SUCCESS)
                {
                    irxError("%s: Cannot left rotate tree.\n", routine);
                    goto RETURN;
                }
                w = RIGHT(PARENT(x));
            }
            if ((COLOR(LEFT(w)) == BLACK) &&
                (COLOR(RIGHT(w)) == BLACK))
            {
                COLOR(w) = RED;
                x = PARENT(x);
            }
            else
            {
                if (COLOR(RIGHT(w)) == BLACK)
                {
                    COLOR(LEFT(w)) = BLACK;
                    COLOR(w) = RED;
                    if (rbRightRotate(T, w) != SUCCESS)
                    {
                        irxError("%s: Cannot right rotate tree.\n", routine);
                        goto RETURN;
                    }
                    w = RIGHT(PARENT(x));
                }
                COLOR(w) = COLOR(PARENT(x));
                COLOR(PARENT(x)) = BLACK;
                COLOR(RIGHT(w)) = BLACK;
                if (rbLeftRotate(T, PARENT(x)) != SUCCESS)
                {
                    irxError("%s: Cannot left rotate tree.\n", routine);
                    goto RETURN;
                }
                x = ROOT(T);
            }

        }
        else
        {
            w = LEFT(PARENT(x));
            if (COLOR(w) == RED)
            {
                COLOR(w) = BLACK;
                COLOR(PARENT(x)) = RED;
                if (rbRightRotate(T, PARENT(x)) != SUCCESS)
                {
                    irxError("%s: Cannot left rotate tree.\n", routine);
                    goto RETURN;
                }
                w = LEFT(PARENT(x));
            }
            if ((COLOR(RIGHT(w)) == BLACK) &&
                (COLOR(LEFT(w)) == BLACK))
            {
                COLOR(w) = RED;
                x = PARENT(x);
            }
            else
            {
                if (COLOR(LEFT(w)) == BLACK)
                {
                    COLOR(RIGHT(w)) = BLACK;
                    COLOR(w) = RED;
                    if (rbLeftRotate(T, w) != SUCCESS)
                    {
                        irxError("%s: Cannot right rotate tree.\n", routine);
                        goto RETURN;
                    }
                    w = LEFT(PARENT(x));
                }
                COLOR(w) = COLOR(PARENT(x));
                COLOR(PARENT(x)) = BLACK;
                COLOR(LEFT(w)) = BLACK;
                if (rbRightRotate(T, PARENT(x)) != SUCCESS)
                {
                    irxError("%s: Cannot left rotate tree.\n", routine);
                    goto RETURN;
                }
                x = ROOT(T);
            }
        }
    }

    COLOR(x) = BLACK;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: rbDelete
   Created by: Charles Irwin 3/15/1999
   Purpose: Deletes a node out of a binary tree. Returns a node the calling
            procedure may free.
   --------------------------------------------------------------------------*/

static TTreeNode * rbDelete(
    IrxTBinaryTree *T,                     /* (I) Red-Black tree to delete from */
    TTreeNode   *z)                     /* (I) Node to delete */
{
    static char routine[] = "rbDelete";
    int            status = SUCCESS;

    TTreeNode *x = NULL;
    TTreeNode *y = NULL;

    if (z == NIL(T))
    {
        irxError("%s: Cannot delete NIL node out of tree.\n", routine);
        goto RETURN;
    }

    if ((LEFT(z) == NIL(T)) ||
        (RIGHT(z) == NIL(T)))
    {
        y = z;
    }
    else
    {
        y = rbTreeSuccessor(T, z);
        if (y == NULL)
        {
            irxError("%s: Cannot find tree successor.\n", routine);
            goto RETURN;
        }
    }
    if (LEFT(y) != NIL(T))
    {
        x = LEFT(y);
    }
    else
    {
        x = RIGHT(y);
    }
    PARENT(x) = PARENT(y);
    if (PARENT(y) == NIL(T))
    {
        ROOT(T) = x;
    }
    else
    {
        if (y == LEFT(PARENT(y)))
        {
            LEFT(PARENT(y)) = x;
        }
        else
        {
            RIGHT(PARENT(y)) = x;
        }
    }
    if (y != z)
    {
        KEY(z)  = KEY(y);
        DATA(z) = DATA(y);
    }
    if (COLOR(y) == BLACK)
    {
        if (rbDeleteFixup(T,x) != SUCCESS)
        {
            irxError("%s: Cannot fixup tree.\n", routine);
            goto RETURN;
        }
    }

    --LENGTH(T);

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        y = NULL;
        irxErrorFailure(routine);
    }
    return y;
}

/* --------------------------------------------------------------------------
   Function: rbTreePrint
   Created by: Charles Irwin 3/15/1999
   Purpose: Prints the contents of a binary tree.
   --------------------------------------------------------------------------*/

static void rbTreePrint(
    IrxTBinaryTree          *T,            /* (I) Red-Black tree to print from */
    TTreeNode            *x,            /* (I) Subtree to print */
    IrxTBinaryTreePrintFunc *PRINT,        /* (I) Printing function */
    long                  level,        /* (I) Level of subtree */
    IrxTBool              indent)       /* (I) Auto indent? */
{
    long i;
    if (x != NULL)
    {
        if (x != NIL(T))
        {
            rbTreePrint(T, LEFT(x), PRINT, level + 1, indent);
            if (indent)
                for (i=0; i< level; i++)
                    irxError(" ");
            PRINT(KEY(x), DATA(x));
            rbTreePrint(T, RIGHT(x), PRINT, level + 1, indent);
        }
    }
}

/* --------------------------------------------------------------------------
   Function: rbTreeForEach
   Created by: Noel Yap
   Purpose: Sequentially calls function fn with each node's data in tree T.
   --------------------------------------------------------------------------*/

static int rbTreeForEach(
    IrxTBinaryTree             *T,            /* (I/O) Red-Black tree to print from */
    TTreeNode               *x,            /* (I/O) Subtree to print */
    IrxTBinaryTreeForEachFunc  *fn,           /* (I/O) generate function */
    void                    *client_data)  /* (I/O) client data passed into fn */
{
    int status = FAILURE;

    if (x != NULL)
    {
        if (x != NIL(T))
        {
            if (rbTreeForEach(T, LEFT(x), fn, client_data) != SUCCESS)
            {
                goto RETURN;
            }
            if (fn(T, client_data, DATA(x)) != SUCCESS)
            {
                goto RETURN;
            }
            if (rbTreeForEach(T, RIGHT(x), fn, client_data) != SUCCESS)
            {
                goto RETURN;
            }
        }
    }

    status = SUCCESS;

RETURN:
    return status;
}

/* --------------------------------------------------------------------------
   Function: rbTreeForEachKey
   Created by: Noel Yap
   Purpose: Sequentially calls function fn with each node's key in tree T.
   --------------------------------------------------------------------------*/

static int rbTreeForEachKey(
    IrxTBinaryTree             *T,            /* (I/O) Red-Black tree to print from */
    TTreeNode               *x,            /* (I/O) Subtree to print */
    IrxTBinaryTreeForEachFunc  *fn,           /* (I/O) generate function */
    void                    *client_data)  /* (I/O) client data passed into fn */
{
    int status = FAILURE;

    if (x != NULL)
    {
        if (x != NIL(T))
        {
            if (rbTreeForEachKey(T, LEFT(x), fn, client_data) != SUCCESS)
            {
                goto RETURN;
            }
            if (fn(T, client_data, (void*)KEY(x)) != SUCCESS)
            {
                goto RETURN;
            }
            if (rbTreeForEachKey(T, RIGHT(x), fn, client_data) != SUCCESS)
            {
                goto RETURN;
            }
        }
    }

    status = SUCCESS;

RETURN:
    return status;
}

/* --------------------------------------------------------------------------
   Function: copyRBNode
   Created by: Charles Irwin 7/8/1999
   Purpose: Copies a node in a binary tree.
   --------------------------------------------------------------------------*/
static TTreeNode * copyRBNode(
    TTreeNode           * in,            /* (I) Source node */
    IrxTBinaryTreeCopyFunc * copyDataKey)   /* (I) Function to copy cata & key */
{
    static char routine[] = "copyRBNode";
    int         status    = FAILURE;

    void      *newKey = NULL;
    void      *newData = NULL;
    TTreeNode *out = NULL;

    out = newEmptyRBNode();
    if (out == NULL)
        goto RETURN;

    if (copyDataKey(DATA(in),
                    KEY(in),
                    &newData,
                    &newKey) != SUCCESS)
    {
        goto RETURN;
    }

    COLOR(out) = COLOR(in);
    DATA(out) = newData;
    KEY(out) = newKey;

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
    {
        treeNodeFree(out);
        out = NULL;
        irxErrorFailure(routine);
    }
    return (out);
}

/* --------------------------------------------------------------------------
   Function: rbCopySubTree
   Created by: Lawrence Miller 5/18/1999
   Purpose: Copies a binary tree.
   --------------------------------------------------------------------------*/
static int rbCopySubTree(
    IrxTBinaryTree         * inT,          /* (I) Source tree */
    IrxTBinaryTree         * outT,         /* (I) Target tree */
    IrxTBinaryTreeCopyFunc * copyDataKey,  /* (I) Key and data copy function */
    IrxTBinaryTreeFreeFunc * freeData,     /* (I) Function to free data */
    IrxTBinaryTreeFreeFunc * freeKey,      /* (I) Function to free key */
    TTreeNode           * in,           /* (I) Source node */
    TTreeNode           * outParent,    /* (I) Parent of target node */
    TTreeNode          ** newOut)       /* (O) Target node */
{
    static char routine[] = "rbCopySubTree";
    int         status = FAILURE;
    TTreeNode * out = NULL;

    if (in == NIL(inT))
    {
        out = NIL(outT);
    }
    else
    {
        out = copyRBNode(in, copyDataKey);
        if (out == NULL)
            goto RETURN;

        /*
         * Set the links to NIL(T) for the new node, just in case
         * we have to free the whole subtree below.
         */
        PARENT(out) = outParent;
        LEFT(out) = NIL(outT);
        RIGHT(out) = NIL(outT);

        /*
         * NOTE: If there was a failure here, LEFT(out) would NOT be set,
         * and would still be equal to NIL(outT). Thus, we can free the
         * subtree underneath the out node in a consistent manner.
         */
        if (rbCopySubTree(inT,
                          outT,
                          copyDataKey,
                          freeData,
                          freeKey,
                          LEFT(in),
                          out,
                          &(LEFT(out))) != SUCCESS)
        {
            goto RETURN;
        }

        /*
         * NOTE: If there was a failure here, RIGHT(out) would NOT be set,
         * and would still be equal to NIL(outT). Thus, we can free the
         * subtree underneath the out node in a consistent manner.
         */
        if (rbCopySubTree(inT,
                          outT,
                          copyDataKey,
                          freeData,
                          freeKey,
                          RIGHT(in),
                          out,
                          &(RIGHT(out))) != SUCCESS)
        {
            goto RETURN;
        }
    }

    status = SUCCESS;

    /*
     * Only set the output variable if completely successful.
     */
    *newOut = out;
    out = NULL;

 RETURN:
    /*
     * Will free the entire subtree if there was ANY failure above.
     */
    freeRBSubTree(outT, out, freeData, freeKey);
    if (status != SUCCESS)
    {
        irxErrorFailure(routine);
    }
    return (status);
}

/* --------------------------------------------------------------------------
   Function: irxBinaryTreeGetLength
   Created by: Noel Yap
   Purpose: Returns the number of nodes inside a binary tree
   -------------------------------------------------------------------------- */
int irxBinaryTreeGetLength(
    IrxTBinaryTree *root,  /* (I) node of tree */
    size_t *length) /* (O) number of nodes in the tree */
{
    static  char    routine[] = "irxBinaryTreeGetLength";
    auto    int     status = FAILURE;

    if(root == NULL)
    {
        irxError("%s: Cannot process NULL input\n", routine);

        goto RETURN;
    }

    *length = LENGTH(root);

    status = SUCCESS;

RETURN:
    return status;
}
