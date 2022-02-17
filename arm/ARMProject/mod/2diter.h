
/** \file 2diter.h */

#ifndef _2DITER_H
#define _2DITER_H


#include "geniter.h"


/**
	\class ARM_2DIter
	
	desc
*/
class ARM_2DIter : public ARM_GenIter
{
    protected:

		/*!
			initialisation method of the object
		*/
		void Init(void)
        {
             SetName(ARM_2DITER);
        }

    public:

           
        /*!
			default constructor
		*/
        ARM_2DIter(void)
        {
            ;
        }
    
		/*!
			destructor
		*/
        ~ARM_2DIter(void)
         {
            ;
         }

		/*!
			copy constructor
		*/
        ARM_2DIter(const ARM_2DIter& obj /*! object to be copied */): ARM_GenIter(obj)
        {
 
            SetName(ARM_2DITER);

            this->BitwiseCopy(&obj);
        }

		/*!
			Affectation operator
		*/
        ARM_2DIter& operator = (const  ARM_2DIter &iter /*! object which is affected to \e this */)
        {
            (*this).ARM_GenIter::operator = (iter);
 
            this->BitwiseCopy(&iter);
 
            return(*this);
        }

		/*!
			Method which copies all the attributes of oiter which are specific to the class ARM_2DIter to the pointer \a this
		*/ 
        void BitwiseCopy(const ARM_Object* oiter /*! object to be copied */)
        {
           ;
        }


		/*!
			Method which copies all the attributes of oiter to the pointer \a this
		*/ 
        void Copy(const ARM_Object* iter /*! object to be copied */)
        {
            ARM_GenIter::Copy(iter);

            BitwiseCopy(iter);
        }

		/*!
			cloning method
			
			\return the cloned object
		*/
        ARM_Object* Clone(void)
        {
            ARM_2DIter* theClone = new ARM_2DIter();

            theClone->Copy(this);

            return(theClone);
        }


        
       
};

#endif
/*---------------------------------------------------------------------*/
/*---- End of File ----*/
