/*
 * $Log $
 *
 */
/** \file 2dconstlattice.h */

#ifndef _2DCONSTLATTICE_H
#define _2DCONSTLATTICE_H
#include "armglob.h"
#include "dates.h"
#include "2dlattice.h"

/**
	\class ARM_2DConstLattice
	
	desc
*/
class ARM_2DConstLattice : public ARM_2DLattice
{
	private :
		double itsStepLen; // longueur du pas exprimé en nbre d'années
		double itsDeltaX1 ;
		double itsDeltaX2 ;

	public:
		/*!
			constructor
		*/
		ARM_2DConstLattice( ARM_Date startree /*! start Date of the tree */
							, ARM_Date endtree /*! end Date of the tree */
							, int nbStep /*! number of steps */
							, double var1 /*! var1 */
							, double var2 /*! var2 */);
		
		/*!
			constructor
		*/
		ARM_2DConstLattice( ARM_Date startree /*! start Date of the tree */
							, ARM_Date endtree /*! end Date of the tree */
							, int nbStep /*! number of steps */);

		/*!
			default constructor
		*/
		ARM_2DConstLattice(void)
		{
			itsStepLen = 0.0 ;
			itsDeltaX1 = 0.0 ;
			itsDeltaX2 = 0.0 ;
			
		}

		/*!
			copy constructor
		*/
		ARM_2DConstLattice(const ARM_2DConstLattice& obj /*! object to be copied */):ARM_2DLattice (obj)
		{
			SetName(ARM_2DCONSTLATTICE);
			BitwiseCopy(&obj);
		}


		/*!
			Method which copies all the attributes of o2dconstlattice which are specific to the class ARM_2DConstLattice to the pointer \a this
		*/ 
		void BitwiseCopy(const ARM_Object* o2dconstlattice /*! object to be copied */)
		{
			ARM_2DConstLattice* dconstlattice = (ARM_2DConstLattice*) o2dconstlattice;
			itsDeltaX1 = dconstlattice->itsDeltaX1;
			itsDeltaX2 = dconstlattice->itsDeltaX2;
			itsStepLen = dconstlattice->itsStepLen;
		}

		/*!
			Method which copies all the attributes of dconstlattice to the pointer \a this
		*/ 
		void Copy(const ARM_Object* dconstlattice /*! object to be copied */)
		{
			ARM_Object::Copy(dconstlattice);
			BitwiseCopy(dconstlattice);
		}

		/*!
			cloning method
			
			\return the cloned object
		*/
		ARM_Object* Clone(void)
		{
			ARM_2DConstLattice* theClone = new ARM_2DConstLattice();
			theClone->Copy(this);
			return(theClone);
		}

		/*!
			Accessor method
			\return \b itsStepLen attribute corresponding to the length of the step in years
		*/
		double GetDt(int index = 0 /*! index */)
		{
			return itsStepLen;
		}

		/*!
			Same as GetDt
		*/
		double StepLen(int index /*! index */) {return itsStepLen;};

		/*!
			Accessor method
			\return \b itsDeltaX1 attribute
		*/
		double GetDeltaX1 (int index =0 /*! index */) {return itsDeltaX1 ;};

		/*!
			Accessor method
			\return \b itsDeltaX2 attribute
		*/
		double GetDeltaX2 (int index = 0 /*! index */) {return itsDeltaX2;};


	
		/*! 
			Conversion method
		 	\return index which preceds the date
		 */
		int DateToIndex(ARM_Date tdate /*! date */)
		{
			double tol = 10e-8;
			return(int(((tdate.GetJulian()-itsStartTree.GetJulian())/K_YEAR_LEN)/itsStepLen + tol));
		}

		/*!
			\return the distance (in years) between the \b tindex index in the tree and the start date of the tree
		*/
		double IndexToTimeIndex(int tindex /*! index */)
		{
			return tindex* itsStepLen;
		}

		/*!
			\return the distance (in years) between \b datecour and the node which preceds \b datecour
		*/
		double DistFromPrevNode(ARM_Date datecour /*! date */)
		{
			double distfromstart = (datecour - itsStartTree)/K_YEAR_LEN;
			int idxprev = DateToIndex(datecour);
			return distfromstart - (idxprev*itsStepLen);
		}

		/*!
			Method which updates \b itsDeltaX1 (\f$=\sqrt(MaxVar1\f$) and \b itsDeltaX2 (\f$=\sqrt(MaxVar1\f$) attributes
		*/
		void CompleteDeltaX( double MaxVarl, double MaxVar2);
	    
	    /*!
	    	see the code
	    */
	    int GetNbNode1(int index)
		{
			return ( 2 * index +1 )  ;
		}
 
	    /*!
	    	see the code
	    */
	    int GetNbNode2(int index) 
		{
			return ( 2 * index +1 )  ;
		}
	};


#endif
/*-----------------------------------------------------------------------*/
