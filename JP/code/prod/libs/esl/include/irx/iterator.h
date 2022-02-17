/*
***************************************************************************
** HEADER FILE: iterator.h
** CREATED BY:  Peter Taylor (April 2003)
**
** Iterator functions which work well with container classes.
**
** $Header: /home/alibcvs/alibcvs/.cvsroot-alib/utils/include/iterator.h,v 1.1 2003/04/02 16:23:54 ptaylor Exp $
***************************************************************************
*/

#ifndef IRX_ITERATOR_H
#define IRX_ITERATOR_H

#include "bintree.h"
#include "hash.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*t
***************************************************************************
** An iterator action function acts upon all items within a container.
**
** Container classes provided are hash table and binary tree.
**
** The iterator function will have three parameters - the container, the
** iteraction action function, and data (as void*) which is passed as
** a parameter in each iterator function.
**
** If the iterator action function does not return SUCCESS, then the loop
** will stop, and the return code from the action function will be returned
** by the iterator function. This enables you to distinguish between cases
** where the iteractor action failed unexpectedly and when you want to
** end the loop in a controlled fashion.
***************************************************************************
*/
typedef int (IrxTIteratorActionFunc) (
    void *item,  /* (I/O) The current item from the container */
    void *data); /* (I/O) Data provided by high level function */

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
(IrxTHashTable          *hashTable,  /* (I) Hash table */
 IrxTIteratorActionFunc *actionFunc, /* (I) Action function */
 void                   *data);      /* (I) Passed to each call of actionFunc */

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
 void                   *data);      /* (I) Passed to each call of actionFunc */

#ifdef __cplusplus
}
#endif


#endif /* IRX_ITERATOR_H */

