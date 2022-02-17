/*
 * $Log: 1dlattice.h,v $
 * Revision 1.1  2003/03/28 15:19:02  mab
 * Initial revision
 *
 * Revision 1.2  2002/04/18 13:36:34  jpriaude
 * Modif pour doxygen
 *
 * Revision 1.1  1999/11/24 17:16:57  vberger
 * Initial revision
 *
 *
 */

/** \file 1dlattice.h */

/*----------------------------------------------------------------------*/
#ifndef _1DLATTICE_H
#define _1DLATTICE_H




#include "gentree.h"

/**
	\class ARM_1DLattice
	
  	desc
*/

class ARM_1DLattice : public ARM_GenTree
{

    public:

		/*!
			constructor
		*/
        ARM_1DLattice( 	ARM_Date startree /*! start Date of the tree */
        				,ARM_Date endtree /*! end Date of the tree */
        				,int nbStep /*! number of steps */);
        
		/*!
			default constructor
		*/
        ARM_1DLattice(void)
        {;}

		/*!
			copy constructor
		*/
        ARM_1DLattice(const ARM_1DLattice& obj /*! object to be copied */): ARM_GenTree(obj)
        {
            SetName(ARM_1DLATTICE);

            BitwiseCopy(&obj);
        }

		/*!
			destructor
		*/
       ~ARM_1DLattice(void)
        {;}

		/*!
			Method which copies all the attributes of ogentree which are specific to the class ARM_1DLattice to the pointer \a this
		*/ 
        void BitwiseCopy(const ARM_Object* ogentree/*! object to be copied */)
        {;}


		/*!
			Method which copies all the attributes of gentree to the pointer \a this
		*/ 
        void Copy(const ARM_Object* gentree/*! object to be copied */)
        {
            ARM_GenTree::Copy(gentree);

            BitwiseCopy(gentree);
        }

		/*!
			cloning method
			
			\return the cloned object
		*/
        ARM_Object* Clone(void)
        {
            ARM_1DLattice* theClone = new ARM_1DLattice();

            theClone->Copy(this);

            return(theClone);
        }

		/*!
			caracteristic of the tree
			\return 0
		*/
        virtual int Width(int tindex /*! index */) {return 0;};

		/*!
			caracteristic of the tree
			\return 0
		*/
        virtual int Jmax(int tindex /*! index */) {return 0;};

		/*!
			caracteristic of the tree
			\return 0
		*/
        virtual int MaxJMax(void) {return 0;};

		/*!
			caracteristic of the tree
			\return 0
		*/
        virtual int WidthMax(void) {return 0;};


		/*!
			Conversion function
			\return 0
		*/
        int DateToIndex(ARM_Date tdate /*! date */) {return 0;};
        
		/*!
			Conversion function
			\return 0
		*/
        double IndexToTimeIndex(int tindex /*! index */) {return 0.0;};

		/*!
			Conversion function
			\return 0
		*/
        double DistFromPrevNode(ARM_Date date /*! date */) {return 0.0;};


       
};


#endif
/*-----------------------------------------------------------------------*/
