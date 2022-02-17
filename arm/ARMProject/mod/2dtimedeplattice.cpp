/*
 * $Log: 2dtimedeplattice.cpp,v $
 * Revision 1.2  2003/06/23 08:05:27  ebenhamou
 * remove unused var
 *
 * Revision 1.1  2003/03/28 15:21:31  mab
 * Initial revision
 *
 * Revision 1.7  2001/10/24 08:52:36  vberger
 * gestions des leaks : ajout de "delete locvar2"
 *
 * Revision 1.6  2001/01/29 17:15:10  vberger
 * modif iostream et utilisation de la stl
 *
 * Revision 1.4  2000/12/27 12:09:41  vberger
 * suppression de #include "iostream.h" remplace par #include <iostream>
 *
 */
#include "2dtimedeplattice.h"
#include <iostream>
#ifndef __ARM_NO_NAMESPACE
using namespace std;
#endif

ARM_2DtDepLattice::ARM_2DtDepLattice(ARM_Date startree, ARM_Date endtree,
				int nbStep)
:ARM_2DLattice(startree,  endtree,  nbStep)

{
Init();
itsFreqNbStep = nbStep ;
}
void ARM_2DtDepLattice::Init()
{
		itsDeltaX1 = 0.0;
		itsDeltaX2ChangeFreq =3.0 ;
		itsFreqNbStep = 1 ;
			
		itsDeltaX2Mode = NULL;
		itsNbNode = NULL;
		itsDt = NULL ;
		itsCumulDt = NULL ;	
		itsDeltaX2 = NULL ;
		itsMinX2 = NULL;
			
			
}

void ARM_2DtDepLattice::CompleteLattice(double dtprime, int nbdtprime, double  dt, int nbstep, double meanrev)
{
	int i ;
	//VBMODIF1
	itsDeltaX2ChangeFreq =-log(0.4)/2.0/meanrev ;
//	cout<<"***********itsDeltaX2ChangeFreq in CompleteLattice "<<itsDeltaX2ChangeFreq<<endl;
	SebtNbStep (nbstep);
	if (itsDt)
	   delete itsDt;
	itsDt = new ARM_Vector( nbstep ) ;
	
	if (itsCumulDt)
	   delete itsCumulDt;
	itsCumulDt = new ARM_Vector( nbstep ) ;  
	
  	itsDt->Elt(0)= dtprime ;
	itsCumulDt->Elt(0) = 0;
	   
	for (i = 1 ; i < nbdtprime; i++ )
	{
		itsDt->Elt(i)= dtprime ;
		itsCumulDt->Elt(i)= itsCumulDt->Elt(i-1) + itsDt->Elt(i-1) ;
	}

	
	for (i = nbdtprime ; i < nbstep; i++ )
	{
		itsDt->Elt(i)= dt;
		itsCumulDt->Elt(i)= itsCumulDt->Elt(i-1) + itsDt->Elt(i-1)   ;
	}
		

}


void ARM_2DtDepLattice::CompleteDeltaX2( double deltax2 , ARM_Vector * locvar2)
{
	int i = 0;
	
	if (itsDeltaX2)
	   delete itsDeltaX2;
	itsDeltaX2 = new ARM_Vector( NbSteps() ) ;	
		
	if (itsMinX2)
	   delete itsMinX2;
	itsMinX2 = new ARM_Vector( NbSteps() ) ; 
	
	if (itsNbNode)
	   delete itsNbNode;
	   
	itsNbNode 		= new int[NbSteps()] ;   
	itsDeltaX2Mode 	= new int[NbSteps()] ;  

//premiere slice 

	itsNbNode[0] = 1;
	itsMinX2->Elt(i) = 0.0;
    
	while ( i < NbSteps()-1 && itsCumulDt->Elt(i)< itsDeltaX2ChangeFreq ) 
	{
		itsDeltaX2->Elt(i) = deltax2 ;
		itsMinX2->Elt(i+1) = -((double)(i+1)) * deltax2;
		itsNbNode[i+1] =  2*(i+1) + 1 ;	
		itsDeltaX2Mode[i] = NOCHG ;
		i++;
	}
	int nextchg = DOWNCHG;
	double currentdeltaX2 = deltax2 ;
	
	while ( i <  NbSteps() - 1) 
	{
		currentdeltaX2 += currentdeltaX2 ;
		//cas de la slice de transition
		if ( locvar2->Elt(i) <  currentdeltaX2*currentdeltaX2  )
		{	
			if (  itsNbNode[i]%2 == 1  )//cas ou le nb de noeud de la slice est impair
			{
			
				itsNbNode[i+1] = (itsNbNode[i] - 1) / 2  + 3 ;
				itsMinX2->Elt(i+1) = itsMinX2->Elt(i) - currentdeltaX2 ;
				itsDeltaX2Mode[i] = SIMPLECHG ;
			}
			else 
			{
				itsNbNode[i+1] = itsNbNode[i] / 2  + 2 ;
				if ( nextchg == DOWNCHG )
					{
					itsMinX2->Elt(i+1) = itsMinX2->Elt(i) - currentdeltaX2 ;
					itsDeltaX2Mode[i] = nextchg ;
					nextchg = UPCHG ;
					}
				else
					{
					itsMinX2->Elt(i+1)= itsMinX2->Elt(i) - currentdeltaX2/2.0 ;
					itsDeltaX2Mode[i] = nextchg ;
					nextchg = DOWNCHG ;
					}
			}
			itsDeltaX2->Elt(i) = currentdeltaX2 ;
			i++;
		}
		//cas general
		
		while ( i <  NbSteps()-1 && ( locvar2->Elt(i) <  currentdeltaX2*currentdeltaX2  ) )
		{
			itsDeltaX2->Elt(i)= currentdeltaX2 ;
			itsMinX2->Elt(i+1)= itsMinX2->Elt(i) - currentdeltaX2 ;
			itsNbNode[i+1] = itsNbNode[i] + 2 ;
			itsDeltaX2Mode[i] = NOCHG ;
			i++;
		}
		
					
	}

	delete locvar2;

/*	for ( i = 0 ; i <  NbSteps()-1 ; i++)
		cout<<i<<"  nbr de noeud  "<<itsNbNode[i]<<" "<<
				itsDeltaX2Mode[i]<<"   "<<itsDeltaX2->Elt(i)<<"  "<<itsMinX2->Elt(i)
				<<" "<<itsDt->Elt(i)<<" "<<itsCumulDt->Elt(i)<<endl  ;*/
	
}

int ARM_2DtDepLattice::DateInYearToIndex(double date)
{
	int maxindex = itsDt->GetSize()- 1;
	int minindex = 0 ;
	if (date < 0.0 || date > itsCumulDt->Elt(maxindex) + 0.1/365.0)
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
			 "la date excede les limite du lattice");
	int i = maxindex / 2 ;
	int test = 0 ;
	if (date >= itsCumulDt->Elt(i) &&  date < itsCumulDt->Elt(i+1))
		test = 1 ;
	while ( test == 0 )
	{
	if (date > itsCumulDt->Elt(i))
		{
		minindex = i ;
		i = minindex + MAX ( (  maxindex - minindex ) / 2  ,1);
		}
	else 
		{
		maxindex  = i;
		i = ( maxindex + minindex ) / 2;
		}

	if (i == maxindex  ||  ( date >= itsCumulDt->Elt(i)  &&  date < itsCumulDt->Elt(i+1) ))
		test = 1 ;	
	}
	return i;

}
int ARM_2DtDepLattice::ClosestIndex(ARM_Date tdate)
{
	int i = DateToIndex(tdate);
	double date =  (tdate.GetJulian()-itsStartTree.GetJulian())/K_YEAR_LEN ;
	if (i < itsDt->GetSize()-1)
		if (fabs(itsCumulDt->Elt(i) - date) > fabs(itsCumulDt->Elt(i+1) - date))
			i ++ ;
	return i ;


}

int ARM_2DtDepLattice::ClosestIndex(double dateinyear)
{
	int i = DateInYearToIndex(dateinyear);
	if (i < itsDt->GetSize()-1)
		if (fabs(itsCumulDt->Elt(i) - dateinyear) > fabs(itsCumulDt->Elt(i+1) - dateinyear))
			i ++ ;
	return i ;


}
