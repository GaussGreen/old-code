/*
 * $Log: tree.h,v $
 * Revision 1.3  2003/06/30 09:29:48  ebenhamou
 * added the correct header for exception
 *
 * Revision 1.2  2003/06/30 08:54:51  ebenhamou
 * initial version log + robust code
 *
 *
 */

/*-----------------------------------------------------------------------------*/
/*                                                                             */
/*   File		    : tree.h                                                   */
/*                                                                             */
/*   Description	: Binary Tree class include file                           */
/*                                                                             */
/*                                                                             */
/*-----------------------------------------------------------------------------*/
#ifndef _TREE_H
#define _TREE_H





#include "armglob.h"
#include "armval.h"
#include "expt.h"






/*----------------------------------------------------------------------------*/


inline int ARM_TREE_MAX(int a, int b)
{
    if ( a < b )
    {
       return(b);
    }

    return(a);
}


/*---- Only necessary methods are implemented ----*/


class ARM_BinaryTree : public ARM_Object
{
                   ARM_BinaryTree* leftChild;
                   ARM_BinaryTree* rightChild; 

    public:

                   ARM_BinaryTree(void) 
                   { 
                       leftChild  = (ARM_BinaryTree *) NULL;

                       rightChild = (ARM_BinaryTree *) NULL;
                   }
 
                   ARM_BinaryTree(ARM_BinaryTree* l, ARM_BinaryTree* r)
                   {
                       leftChild = l;

                       rightChild = r;      
                   }

      virtual     ~ARM_BinaryTree(void);

                   short Height(void)
                   {
                       if ( ( leftChild == (ARM_BinaryTree *) NULL ) &&
                            ( rightChild == (ARM_BinaryTree *) NULL )
                          )
                       {
                          return(0);
                       }                                    
 
                       if ( leftChild == (ARM_BinaryTree *) NULL )
                       {
                          return(1 + rightChild->Height());
                       } 

                       if ( rightChild == (ARM_BinaryTree *) NULL )
                       {
                          return(1 + leftChild->Height());
                       } 

                       return(1 + ARM_TREE_MAX(leftChild->Height(), 
                                          rightChild->Height()));
                   }

                   void InsertLeft(ARM_BinaryTree* left)    
                   { 
                       leftChild = left;
                   }

                   void InsertRight(ARM_BinaryTree* right)
                   { 
                       rightChild = right;
                   }

                   ARM_BinaryTree* GetLeftChild(void)  
                   { 
                       return(leftChild);
                   }

                   ARM_BinaryTree* GetRightChild(void)
                   { 
                       return(rightChild); 
                   }
};



/*-----------------------------------------------------------------------------*/
/* The ARM_BinaryTree class destructor                                         */
/*-----------------------------------------------------------------------------*/

inline ARM_BinaryTree::~ARM_BinaryTree(void)
{
    if (leftChild)
       delete leftChild;

    if (rightChild)
       delete rightChild;
}



/*-----------------------------------------------------------------------------*/



/*--- The Balanced tree (AVL tree) ---*/


class ARM_ValueTable : public ARM_BinaryTree
{
      private:

                      ARM_Val* var;

      public:

                      ARM_ValueTable(void)         
                      { 
                          var = (ARM_Val *) NULL; 
                      }

                      ARM_ValueTable(ARM_Val* v)   
                      { 
                          var = v; 
                      }

                      ARM_ValueTable(ARM_Val* v, ARM_ValueTable* l, 
                                     ARM_ValueTable* r) : ARM_BinaryTree(l,r)
                      {
                          var = v;
                      }

                     ~ARM_ValueTable(void) 
                      {
                          var = (ARM_Val *) NULL;
                      }

                      ARM_Val* GetValue(void)  
                      { 
                          return(var); 
                      }

                      void SetValue(ARM_Val* v) 
                      { 
                          var = v; 
                      }

                      ARM_ValueTable* RotateRight(void);

                      ARM_ValueTable* RotateLeft(void);

                      ARM_ValueTable* RotateLeftRight(void);

                      ARM_ValueTable* RotateRightLeft(void);

                      ARM_ValueTable* BalanceTree(void);

                      ARM_ValueTable* InsertValue(ARM_Val* v);

                      ARM_Val* SearchValue(ARM_Val* v);

                      // TMP ARM_Val* SearchValueWithOrder(short order);

                      void CountVars(short* nbVars);
};



/*-----------------------------------------------------------------------------*/
/* calculate Height difference (leftChild, rightChild)                         */
/*-----------------------------------------------------------------------------*/

inline short deseq(ARM_ValueTable* t)
{
    if ( t == (ARM_ValueTable *) NULL )
    {
       return(0);
    }

    if (( t->GetLeftChild() == (ARM_BinaryTree *) NULL ) 
        &&
        ( t->GetRightChild() == (ARM_BinaryTree *) NULL )
       )
    {
       return(0);
    }

    if ( t->GetLeftChild() == (ARM_BinaryTree *) NULL )
    {
       return(0 - t->GetRightChild()->Height());
    }

    if ( t->GetRightChild() == (ARM_BinaryTree *) NULL )
    {
       return(0 + t->GetLeftChild()->Height());
    }

    return(t->GetLeftChild()->Height() - t->GetRightChild()->Height());
}



/*-----------------------------------------------------------------------------*/
/* rotate right a binary tree                                                  */
/*-----------------------------------------------------------------------------*/

inline ARM_ValueTable* ARM_ValueTable::RotateRight(void)
{
    ARM_ValueTable* aux;


    aux = (ARM_ValueTable *) this->GetLeftChild();

    if (this->GetLeftChild())
    {
       this->InsertLeft(aux->GetRightChild());

       aux->InsertRight(this);
   
       return(aux);
    }
    else
    {
       return(this);
    }
}
 


/*-----------------------------------------------------------------------------*/
/* rotate left a binary tree                                                   */
/*-----------------------------------------------------------------------------*/

inline ARM_ValueTable* ARM_ValueTable::RotateLeft(void)
{
    ARM_ValueTable* aux;
 
    aux = (ARM_ValueTable *) this->GetRightChild();
 
    if (this->GetRightChild())
    {
       this->InsertRight(aux->GetLeftChild());
 
       aux->InsertLeft(this); 
 
       return(aux);
    }
    else
    {
       return(this);
    }
}



/*-----------------------------------------------------------------------------*/
/* rotate left and right                                                       */
/*-----------------------------------------------------------------------------*/
 
inline ARM_ValueTable* ARM_ValueTable::RotateLeftRight(void)
{
    if (this->GetLeftChild())
    {
       this->InsertLeft(((ARM_ValueTable *) this->GetLeftChild())->RotateLeft());

       return(this->RotateRight());
    }

    return(this);
} 



/*-----------------------------------------------------------------------------*/
/* rotate right and left                                                       */
/*-----------------------------------------------------------------------------*/
 
inline ARM_ValueTable* ARM_ValueTable::RotateRightLeft(void)
{
    if (this->GetRightChild())
    {
       this->InsertRight(((ARM_ValueTable *) 
                          this->GetRightChild())->RotateRight());

       return(this->RotateLeft());
    }

    return(this);
}



/*-----------------------------------------------------------------------------*/
/* Balance an AVL tree                                                         */
/*-----------------------------------------------------------------------------*/

inline ARM_ValueTable* ARM_ValueTable::BalanceTree(void)
{
    if (( deseq(this) == 0 ) || ( deseq(this) == 1 ) || ( deseq(this) == -1 ))
    {
       return(this);
    }
    
    if (( deseq(this) == +2 ) 
        && 
        ( deseq((ARM_ValueTable *) this->GetLeftChild()) == +1 ) 
       )
    {
       return(this->RotateRight());
    }

    if (( deseq(this) == -2 ) 
        && 
        ( deseq((ARM_ValueTable *) this->GetRightChild()) == -1 )
       )
    {
       return(this->RotateLeft());
    }

  
    if (( deseq(this) == +2 ) 
        && 
        ( deseq((ARM_ValueTable *) this->GetLeftChild()) == -1 ) 
       )
    {
       return(this->RotateLeftRight());
    }

    if (( deseq(this) == -2 ) 
        && 
        ( deseq((ARM_ValueTable *) this->GetRightChild()) == +1 )
       )    
    { 
       return(this->RotateRightLeft()); 
    } 

	// other case
    throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
         "Problem with reorganising the tree");


}



/*-----------------------------------------------------------------------------*/
/* walk int the tree and search an object                                      */
/*-----------------------------------------------------------------------------*/

inline ARM_Val* ARM_ValueTable::SearchValue(ARM_Val* v)
{
    if ( *v == *var )
    {
       return(var);
    }

    if ( *v < *var )
    {
       if ( this->GetLeftChild() == (ARM_ValueTable *) NULL )   
       {
          return((ARM_Val *) NULL);
       }

       return(((ARM_ValueTable *)this->GetLeftChild())->SearchValue(v));
    }
    else
    {
       if ( this->GetRightChild() == (ARM_ValueTable *) NULL )   
       {
          return((ARM_Val *) NULL);
       } 

       return(((ARM_ValueTable *) this->GetRightChild())->SearchValue(v));
    }
}




/*-----------------------------------------------------------------------------*/
/* count the total number of variables                                         */
/*-----------------------------------------------------------------------------*/

inline void ARM_ValueTable::CountVars(short* nbVars)
{
    (*nbVars)++;

    if ( this->GetLeftChild() != (ARM_ValueTable *) NULL )
    {
       ((ARM_ValueTable *)this->GetLeftChild())->CountVars(nbVars);
    }
 
    if ( this->GetRightChild() != (ARM_ValueTable *) NULL ) 
    {    
       ((ARM_ValueTable *)this->GetRightChild())->CountVars(nbVars);
    } 
}





/*-----------------------------------------------------------------------------*/
/* insert in AVL tree                                                          */
/*-----------------------------------------------------------------------------*/

inline ARM_ValueTable* InsertInAVLtree(ARM_ValueTable* t, ARM_Val* v)
{
    if ( t == (ARM_ValueTable *) NULL )
    {
       return(new ARM_ValueTable(v));
    }

    if ( *v < *(t->GetValue()) )
    {
       t->InsertLeft(InsertInAVLtree((ARM_ValueTable *) t->GetLeftChild(), v));

       return(t->BalanceTree());
    }

    if ( *(t->GetValue()) < *v )
    {
       t->InsertRight(InsertInAVLtree((ARM_ValueTable *) t->GetRightChild(), v));

       return(t->BalanceTree());
    }

    return(t);
}



/*-----------------------------------------------------------------------------*/
/* insert a value in the values table(a tree)                                  */
/*-----------------------------------------------------------------------------*/

inline ARM_ValueTable* ARM_ValueTable::InsertValue(ARM_Val* v)
{
    return(InsertInAVLtree(this, v));
}







#endif /* end of tree.h */
/*-----------------------------------------------------------------------------*/
