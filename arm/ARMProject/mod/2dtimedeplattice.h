/*
 * $Log: 2dtimedeplattice.h,v $
 * Revision 1.1  2003/03/28 15:22:34  mab
 * Initial revision
 *
 * Revision 1.6  2002/04/18 13:44:21  jpriaude
 * Modif pour Doxygen
 *
 * Revision 1.5  2001/06/21 18:03:12  vberger
 * correction du destructeur : reinitialisation des vecteur à NULL
 *
 * Revision 1.4  2000/11/20 15:34:13  vberger
 * Init passe dans le .cpp
 * ajout d'un arg dans CompleteLattice
 *
 * Revision 1.3  2000/09/07 13:16:50  mab
 * *** empty log message ***
 *
 * Revision 1.2  2000/02/10 13:34:00  vberger
 * correction pour le log message
 *
 *
 */

/** \file 2dtimedeplattice.h */

#ifndef _2DTDEPTLATTICE_H
#define _2DTDEPTLATTICE_H
#include "armglob.h"
#include "dates.h"
#include "2dlattice.h"

#define SIMPLECHG 1
#define UPCHG 2
#define DOWNCHG 3
#define NOCHG 4


/**
	\class ARM_2DtDepLattice
	
	desc
*/
class ARM_2DtDepLattice : public ARM_2DLattice
{
	private :
		double itsDeltaX1 ;
		
		int itsFreqNbStep ;
		double itsDeltaX2ChangeFreq ;
		
		int 		* itsNbNode 	;
		int 		* itsDeltaX2Mode 	;
		ARM_Vector  * itsDt			; 
		ARM_Vector  * itsCumulDt	;	
		ARM_Vector 	* itsDeltaX2 	;
		ARM_Vector 	* itsMinX2 		;
		

	public:
	
		/*!
			constructor
		*/
		ARM_2DtDepLattice( 	ARM_Date startree /*! start Date of the tree */
							,ARM_Date endtree /*! end Date of the tree */
							,int nbStep /*! number of steps */);

		/*!
			default constructor
		*/
		ARM_2DtDepLattice(void)
		{
			Init();	
		}
		

		/*!
			copy constructor
		*/
		ARM_2DtDepLattice(const ARM_2DtDepLattice& obj /*! object to be copied */):ARM_2DLattice (obj)
		{
			Init();
			SetName(ARM_2DTDEPLATTICE);
			BitwiseCopy(&obj);
		}
		
		/*!
			initialisation method
		*/
		void Init();
		
		/*!
			destructor
		*/
		 ~ARM_2DtDepLattice(void)
        {
            
            if (itsDeltaX2Mode)
            {
               delete itsDeltaX2Mode;
               itsDeltaX2Mode = NULL;
			}
               
            if (itsNbNode)
            {
               delete itsNbNode;
               itsNbNode = NULL;
			}
               
            if (itsDt)
            {
               delete itsDt;
               itsDt = NULL;
			}

			if (itsCumulDt)
			{
               delete itsCumulDt;
               itsCumulDt = NULL;
			}

			if (itsDeltaX2)
			{
               delete itsDeltaX2;
               itsDeltaX2 = NULL;
			}

			if (itsMinX2)
			{
               delete itsMinX2;
               itsMinX2 = NULL;
			}
        }
        
		/*!
			Method which copies all the attributes of o2dtdeplattice which are specific
			to the class ARM_2DtDepLattice to the pointer \a this
		*/ 
		void BitwiseCopy(const ARM_Object* o2dtdeplattice /*! object to be copied */)
		{
			ARM_2DtDepLattice* dtdeplattice = (ARM_2DtDepLattice*) o2dtdeplattice;
			
			itsDeltaX1 = dtdeplattice->itsDeltaX1;
			
			itsFreqNbStep =dtdeplattice->itsFreqNbStep;
			
			itsDeltaX2ChangeFreq =dtdeplattice->itsDeltaX2ChangeFreq;
			
			
			
			if (itsDeltaX2Mode)
               delete itsDeltaX2Mode;

            if (dtdeplattice->itsDeltaX2Mode)
			{
				itsDeltaX2Mode = new int[dtdeplattice->NbSteps()];
				for (int i = 0 ; i< dtdeplattice->NbSteps() ; i++)
					itsDeltaX2Mode[i] = dtdeplattice->GetDeltaX2Mode(i);
			}
			
			if (itsNbNode)
               delete itsNbNode;

            if (dtdeplattice->itsNbNode)
			{
				itsNbNode = new int[dtdeplattice->NbSteps()];
				for (int i = 0 ; i< dtdeplattice->NbSteps() ; i++)
					itsNbNode[i] = dtdeplattice->GetNbNode2(i);
			}
			
                                    
			if (itsDt)
               delete itsDt;

            if (dtdeplattice->itsDt)
               itsDt = (ARM_Vector *) 
                                   dtdeplattice->itsDt->Clone();
			if (itsCumulDt)
               delete itsCumulDt;

            if (dtdeplattice->itsCumulDt)
               itsCumulDt = (ARM_Vector *) 
                                   dtdeplattice->itsCumulDt->Clone();
                                   
			if (itsDeltaX2)
               delete itsDeltaX2;

            if (dtdeplattice->itsDeltaX2)
               itsDeltaX2 = (ARM_Vector *) 
                                   dtdeplattice->itsDeltaX2->Clone();
			if (itsMinX2)
               delete itsMinX2;

            if (dtdeplattice->itsMinX2)
               itsMinX2 = (ARM_Vector *) 
                                   dtdeplattice->itsMinX2->Clone();
                                   
		}

		/*!
			Method which copies all the attributes of dtdeplattice to the pointer \a this
		*/ 
		void Copy(const ARM_Object* dtdeplattice /*! object to be copied */)
		{
			ARM_Object::Copy(dtdeplattice);
			BitwiseCopy(dtdeplattice);
		}

		/*!
			cloning method
			
			\return the cloned object
		*/
 		ARM_Object* Clone(void)
		{
			ARM_2DtDepLattice* theClone = new ARM_2DtDepLattice();
			theClone->Copy(this);
			return(theClone);
		}

		/*!
			accessor method
			\return \a itsFreqNbStep attribute
		*/
	    int 	FreqNbStep(void) {return itsFreqNbStep;};
	    
		/*!
			accessor method
			\return \a itsDeltaX2ChangeFreq attribute
		*/
	    double 	GetDeltaX2ChangeFreq (void) {return itsDeltaX2ChangeFreq ;};

		/*!
			accessor method
			\return \a itsDeltaX2Mode[index] attribute
		*/
	    int	 	GetDeltaX2Mode(int index) {return itsDeltaX2Mode[index];};

		/*!
			\return 2*index+1
		*/
	    int GetNbNode1(int index)
		{
			return ( 2 * index +1 )  ;
		}
	    
		/*!
			\return the \f$index\f$ème element of \b itsNbNode attribute
		*/
	    int GetNbNode2(int index) 
		{
			return itsNbNode[index]  ;
		}
		
		/*!
			\return the \f$index\f$ème element of \b itsDt attribute
		*/
		double 	GetDt(int index) {return itsDt->Elt(index);};

		/*!
			\return the \f$index\f$ème element of \b itsCumulDt attribute
		*/
		double 	GetCumulDt(int index) {return itsCumulDt->Elt(index);};

		/*!
			\return \b itsDeltaX1 attribute
		*/
		double 	GetDeltaX1 (int index) {return itsDeltaX1 ;};

		/*!
			\return the \f$index\f$ème element of \b itsDeltaX2 attribute
		*/
		double 	GetDeltaX2 (int index) {return itsDeltaX2->Elt(index);};

		/*!
			\return the \f$index\f$ème element of \b itsMinX2 attribute
		*/
		double 	GetMinX2 (int index) {return itsMinX2->Elt(index);};



		/*!
			conversion function
			\return the index preceding the date (in years) in the tree
		*/
		int DateInYearToIndex(double date) ;

		/*!
			conversion function
			\return the index preceding the date in the tree
		*/
		int DateToIndex(ARM_Date tdate)
		{
			double date =  (tdate.GetJulian()-itsStartTree.GetJulian())/K_YEAR_LEN;
			return ClosestIndex( date );	
		}
		
		/*!
			conversion function
			\return the closest index in the tree corresponding to the date 
		*/
		int ClosestIndex(ARM_Date tdate);

		/*!
			conversion function
			\return the closest index in the tree corresponding to the date (in years)
		*/
		int ClosestIndex(double dateinyear);
		// retourne la distance en fraction d'annee
		// entre le noeud tindex et date
		// de debut de I'arbre

		/*!
			conversion function
			\return the distance (in years) between the date corresponding to \b tindex node and the start date tree
		*/
		double IndexToTimeIndex(int tindex)
		{
			return itsCumulDt->Elt(tindex);
		}
	
		/*!
			conversion function
			\return the distance (in years) between \b datecour and the node preceding \b f=datecour
		*/
		double DistFromPrevNode(ARM_Date datecour)
		{
			double distfromstart = (datecour - itsStartTree)/K_YEAR_LEN;
			int idxprev = DateToIndex(datecour);
			return distfromstart - itsCumulDt->Elt(idxprev);
		}

		/*!
			method which sets \b itsDeltaX1
		*/
		void SetDeltaX1( double deltax1)
		{
			itsDeltaX1 = deltax1;
		}
	
		/*!
			method which updates \b itsDeltaX2, \b itsMinX2, \b itsNbNode and \b itsDeltaX2Mode attributes 
		*/
		void CompleteDeltaX2( double deltax2 , ARM_Vector * locvar2) ;
		
		void CompleteLattice(double dtprime, int nbdtprime, double  dt, int nbstep,double meanrev);

	};


#endif
/*-----------------------------------------------------------------------*/
