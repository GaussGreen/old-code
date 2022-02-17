
/** \file 1diter.h */

#ifndef _1DITER_H
#define _1DITER_H


#include "geniter.h"


/**
	\class ARM_1DIter

  desc
*/

class ARM_1DIter : public ARM_GenIter
{
    protected:


		/*!
			initialisation method of the object
		*/
         void Init(void)
         {
             SetName(ARM_1DITER);
             
         }

    public:

   
        /*!
			default constructor
		*/
        ARM_1DIter(void)
        {
            ;
        }
    
		/*!
			destructor
		*/
        ~ARM_1DIter(void)
         {
           ;
         }
         
		/*!
			copy constructor
		*/
        ARM_1DIter(const ARM_1DIter& obj /*! object to be copied */):ARM_GenIter(obj)
        {
            SetName(ARM_1DITER);

			BitwiseCopy(&obj);
        }


		/*!
			Affectation operator
		*/
        ARM_1DIter& operator = (const ARM_1DIter& iter /*! object which is affected to \a this */)
        {
            (*this).ARM_GenIter::operator = (iter);
 
            this->BitwiseCopy(&iter);
 
            return(*this);
        }
		
		/*!
			Method which copies all the attributes of oiter which are specific to the class ARM_1DIter to the pointer \a this
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
            ARM_1DIter* theClone = new ARM_1DIter();

            theClone->Copy(this);

            return(theClone);
        }


        
        /*!
        	accessor method of the proba Up
        	\return 0
        */
        virtual double GetProbaUp(int j /*! index */) {return 0.0;};

        /*!
        	accessor method of the proba Down
        	\return 0
        */
        virtual double GetProbaDown(int j /*! index */) {return 0.0;};

        /*!
        	accessor method of the proba Mid
        	\return 0
        */
        virtual double GetProbaMiddle(int j /*! index */){return 0.0;};
};

#endif
/*---------------------------------------------------------------------*/
/*---- End of File ----*/
