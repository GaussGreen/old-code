/*
***************************************************************************
** SOURCE FILE: iterator.c
** CREATED BY:  Peter Taylor (April 2003)
**
** Iterator functions which work well with container classes.
**
** $Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/src/iterator.c,v 1.1 2003/04/02 16:23:54 ptaylor Exp $
***************************************************************************
*/

#include "irx/iterator.h"

/*f
***************************************************************************
** Iterates through a hash table.
**
** Returns SUCCESS if each element of the hash table was successfully
** processed by the iterator action function. Otherwise returnes the
** reply code for the first use of the action function which did not return
** success.
***************************************************************************
*/
int irxIterateHashTable
(IrxTHashTable       *hashTable,  /* (I) Hash table */
 IrxTIteratorActionFunc *actionFunc, /* (I) Action function */
 void                *data)       /* (I) Passed to each call of actionFunc */
{
    IrxTBool  more;
    void*  key;
    void* item;

    more = irxHashTableFirst (hashTable, &key, &item);
    while (more)
    {
        int status = actionFunc (item, data);
        if (status != SUCCESS)
            return status;
        more = irxHashTableNext (hashTable, &key, &item);
    }
    return SUCCESS;
}

/*
***************************************************************************
** The existing iterator functionality of the TBinaryTree is not quite
** correct. Therefore we use our own structure and binary tree for each
** function to fix these problems.
***************************************************************************
*/
struct BIN_TREE_ITERATOR_DATA
{
    IrxTIteratorActionFunc *actionFunc;
    void                *data;
    int                  status;
};

/*
***************************************************************************
** This a TBinaryTreeForEach function which for some unknown reason (poor
** design from a previous developer) has a binaryTree as its first 
** parameter.
**
** We use this function as the TBinaryTreeForEach function whenever we
** call irxIterateBinaryTree.
***************************************************************************
*/
static int BinTreeForEachFunc
(IrxTBinaryTree *binaryTree,
 void        *item,
 void        *data)
{
    struct BIN_TREE_ITERATOR_DATA *x = (struct BIN_TREE_ITERATOR_DATA*)data;
    int status = x->actionFunc (item, x->data);

    if (status != SUCCESS)
    {
        x->status = status;
        return FAILURE;
    }

    return SUCCESS;
}

/*f
***************************************************************************
** Iterates through a binary tree.
**
** Returns SUCCESS if each element of the binary tree was successfully
** processed by the iterator action function. Otherwise returnes the
** reply code for the first use of the action function which did not return
** success.
***************************************************************************
*/
int irxIterateBinaryTree
(IrxTBinaryTree         *binaryTree, /* (I) Binary tree */
 IrxTIteratorActionFunc *actionFunc, /* (I) Action function */
 void                *data)       /* (I) Passed to each call of actionFunc */
{
    struct BIN_TREE_ITERATOR_DATA x;

    x.actionFunc = actionFunc;
    x.data       = data;
    x.status     = SUCCESS;

    irxBinaryTreeForEach (binaryTree, BinTreeForEachFunc, (void*)&x);

    return x.status;
}
