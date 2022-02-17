/*
 * $Log: 2dlattice.h,v $
 * Revision 1.1  2003/03/28 15:21:06  mab
 * Initial revision
 *
 * Revision 1.3  2002/04/18 13:42:35  jpriaude
 * Modif pour Doxygen
 *
 * Revision 1.2  1999/12/10 08:25:33  vberger
 *  modifs liees a l'ajout d'un lattice 2 dim a pas non constant
 *
 * Revision 1.1  1999/11/24 17:18:35  vberger
 * Initial revision
 *
 *
 */

/** \file 2dlattice.h */

/*----------------------------------------------------------------------*/
#ifndef _2DLATTICE_H
#define _2DLATTICE_H



#include "gentree.h"


/**
	\class ARM_2DLattice
	
	desc
*/
class ARM_2DLattice : public ARM_GenTree
{

    public:

		/*!
			constructor
		*/
        ARM_2DLattice( 	ARM_Date startree /*! start Date of the tree */
        				,ARM_Date endtree /*! end Date of the tree */
        				,int nbStep /*! number of steps */);
        /*!
        	default constructor
        */
        ARM_2DLattice(void)
        {;}

		/*!
			copy constructor
		*/
        ARM_2DLattice(const ARM_2DLattice& obj /*! object to be copied */): ARM_GenTree (obj)
        {
            SetName(ARM_2DLATTICE);

            BitwiseCopy(&obj);
        }

		/*!
			destructor
		*/
       ~ARM_2DLattice(void)
        {;}

		/*!
			Method which copies all the attributes of olattice which are specific to the class ARM_2DLattice to the pointer \a this
		*/ 
         void BitwiseCopy(const ARM_Object* olattice /*! object to be copied */)
        {;}

		/*!
			Method which copies all the attributes of lattice to the pointer \a this
		*/ 
        void Copy(const ARM_Object* lattice /*! object to be copied */)
        {
            ARM_GenTree::Copy(lattice);

            BitwiseCopy(lattice);
        }

		/*!
			cloning method
			
			\return the cloned object
		*/
        ARM_Object* Clone(void)
        {
            ARM_2DLattice* theClone = new ARM_2DLattice();

            theClone->Copy(this);

            return(theClone);
        }


		/*!
			caracteristic of the tree
			\return 0
		*/
        double StepLen(int index) {return 0.0;};

		/*!
			caracteristic of the tree
			\return 0
		*/
        double GetDt(int index = 0) {return 0.0;};

		/*!
			caracteristic of the tree
			\return \b itsNbStep attribute
		*/
        int NbSteps(void) {return itsNbStep;};

        /*!
        	\throw Exception Unimplemented method
        */
        virtual double 	GetDeltaX1 (int index) 
        {
			throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <GetDeltaX1> method");
		}

        /*!
        	\throw Exception Unimplemented method
        */
		virtual double 	GetDeltaX2 (int index) 
		{
			throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <GetDeltaX2> method");
		}
        
        /*!
        	\throw Exception Unimplemented method
        */
		virtual int GetNbNode1(int index)
		{
			throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <GetNbNode1> method");
		}
 
        /*!
        	\throw Exception Unimplemented method
        */
		virtual int GetNbNode2(int index) 
		{
		
			throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <GetNbNode2> method");
		}
 
		/*!
			conversion function
			\return 0
		*/
        int DateToIndex(ARM_Date tdate /*! date */) {return 0;};
        
		/*!
			conversion function
			\return 0
		*/
        double IndexToTimeIndex(int tindex /*! index */) {return 0.0;};

		/*!
			conversion function
			\return 0
		*/
        double DistFromPrevNode(ARM_Date date /*! date */) {return 0.0;};
       
};


#endif
/*-----------------------------------------------------------------------*/
