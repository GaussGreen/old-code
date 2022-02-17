/*
 * $Log: containr.h,v $
 * Revision 1.2  2003/07/31 16:36:18  mab
 * RCS comment added
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : containr.h                                                   */
/*                                                                            */
/* DESCRIPTION : A Container class Header File                                */
/*                                                                            */
/* DATE        : Wed Dec  4 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _CONTAINR_H
#define _CONTAINR_H


#include "armglob.h"

#include "expt.h"


 
struct ARM_Link
{
    ARM_Object* item;

    ARM_Link*   _next;    



    
    ARM_Link(void)
    {
        item  = (ARM_Object *) NULL;

        _next = (ARM_Link *) NULL;
    }

    ARM_Link(ARM_Object* obj)
    {
        item = obj;

        _next = (ARM_Link *) NULL;
    }

   ~ARM_Link(void)
    {
        item  = (ARM_Object *) NULL;

        _next = (ARM_Link *) NULL;
    }
   
    void Destroy(void)
    {
        if (item)
        {
           delete item;

           item = (ARM_Object *) NULL;
        }

        _next = (ARM_Link *) NULL;

        delete this;
    } 
};





class ARM_Container :  public ARM_Object
{
     private:

                 ARM_Link* last;

                 ARM_Link* curItem; 

                 int       size;

         void Init(void)
         {
             last = NULL;

             curItem = NULL; 

             size = 0;
         }

         void FreeLinks(void)
         {
             if ( last == (ARM_Link *) NULL )
             {
                return;
             }

             ARM_Link* next;
             ARM_Link* cur  = last->_next;

             while ( cur != last )
             {
                 next = cur->_next;
               
                 delete cur;                   

                 cur = next;
             }

             delete last;

             last = (ARM_Link *) NULL;

             size = 0;
         }

         // k = 0 --> return the first element

         ARM_Link* GetNthLink(int k)
         {
              if (( size == 0 ) || ( last == (ARM_Link *) NULL ))
              {
                 throw Exception(__LINE__, __FILE__, ERR_CONTAINER_EMPTY,
                         "ARM_Container Is Empty");

                 return((ARM_Link *) NULL); 
              }
 
              if (( k > (size-1) ) || ( k < 0 ))
              {
                 throw Exception(__LINE__, __FILE__, ERR_CONTAINER_INDEX,
                         "ARM_Container : NOT VALID INDEX");
 
                 return((ARM_Link *) NULL);
              }

              ARM_Link* cur = last;
              ARM_Link* nthItem = (ARM_Link *) NULL;
 
              for (int i = 0; i <= k; i++)
              {
                  cur = cur->_next;
              }
 
              nthItem = cur;

              return(nthItem);
         }

     public:

         ARM_Container(void)
         {
             SetName(ARM_CONTAINER);

             last = (ARM_Link *) NULL;

             curItem = (ARM_Link *) NULL;

             size = 0;
         }

         ARM_Container(int sz)
         {
             ARM_Link* cur = (ARM_Link *) NULL;

 
             Init();

             SetName(ARM_CONTAINER);

             if ( sz <=  0 )
             {
                throw Exception(__LINE__, __FILE__, 
                      ERR_CONTAINER_SIZE_NOT_VALID,
                         "SIZE NOT VALID FOR ARM_Container");
             }

             size = sz;

             for (int i = 0; i < size; i++)
             {
                 if ( i != 0 )
                 {
                    cur = new ARM_Link();

                    cur->_next = last->_next;

                    last->_next = cur;

                    last = cur;
                 }
                 else
                 {    /*--- first item ---*/

                    last = new ARM_Link();

                    last->_next = last; 
                 }
             } 
         }

         ARM_Container(const ARM_Container& container);


         ARM_Container& operator = (const ARM_Container& container);
         

         void BitwiseCopy(const ARM_Object* srcCont);


         void Copy(const ARM_Object* srcCont)
         {
            ARM_Object::Copy(srcCont);

            BitwiseCopy(srcCont);
         }


         ARM_Object* Clone(void)
         {
            ARM_Container* theClone = new ARM_Container();

            theClone->Copy(this);
 
            return(theClone);
         }

         ARM_Container* DuplicateWithContent(void);

        ~ARM_Container(void)
         {
             FreeLinks();
         }

         void Destroy(void)
         {
             if ( last == (ARM_Link *) NULL )
             {
                return;
             }
 
             ARM_Link* next;
             ARM_Link* cur  = last->_next;
 
             while ( cur != last )
             {
                 next = cur->_next;
               
                 cur->Destroy();                   
 
//                 delete cur;

                 cur = next;
             }
 
             last->Destroy();
 
             last = (ARM_Link *) NULL;

             size = 0;

             delete this;
         }

         friend void Destroy(ARM_Container* cont)
         {
             if (cont)
             {
                cont->Destroy();
             }
         }

         ARM_CLASS_NAME GetRootName(void)
         {
              return(ARM_CONTAINER);
         }

         int GetSize(void)
         {
             return(size);
         }

         int IsEmpty(void)
         {
             return(( size == 0 ) || ( last == (ARM_Link *) NULL));
         }

         ARM_Object* GetNthItem(int k)
         {
             ARM_Link* nthItem = (ARM_Link *) NULL;


             nthItem = GetNthLink(k);

             if (nthItem)
             {
                return(nthItem->item);
             }
             else
             {
                return((ARM_Object *) NULL);
             }
         }

         void Append(ARM_Object* obj); 

         void SetItem(ARM_Object* obj, int i);

         ARM_Object* GetFirstItem(void)
         {
             if ( ( size == 0 ) || ( last == (ARM_Link *) NULL ) )
             {
                throw Exception(__LINE__, __FILE__, ERR_CONTAINER_EMPTY,
                         "ARM_Container Is Empty");
             }

             return(last->_next->item);
         }

         ARM_Object* GetLastItem(void)
         {
             if ( ( size == 0 ) || ( last == (ARM_Link *) NULL ) )
             {
                throw Exception(__LINE__, __FILE__, ERR_CONTAINER_EMPTY,
                         "ARM_Container Is Empty");
             }
             
             return(last->item);
         }

         /*--- Some Iteration Methods ----*/

         ARM_Object* Start(void)
         {
            if ( ( size == 0 ) || ( last == (ARM_Link *) NULL ) )
            {
               return((ARM_Object *) NULL); /*--- Empty Container ---*/
            }

            curItem = last->_next;

            return(curItem->item);
         }

         ARM_Object* Next(void)
         {
             if (( size == 0 ) || ( last == (ARM_Link *) NULL )
                 || ( curItem == NULL )
                )
             {
                return((ARM_Object *) NULL); /*--- Empty Container ---*/
             }

             if ( curItem == last )
             {
                curItem = (ARM_Link *) NULL;

                return((ARM_Object *) NULL); /*--- No more objects ---*/
             }
          
             curItem = curItem->_next; 

             return(curItem->item);
         }
};












#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
