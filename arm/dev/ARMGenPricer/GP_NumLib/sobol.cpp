/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file quasirandom.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnumlib/sobol.h"
#include "gpnumlib/sobolinit.h"
#include "gpbase/ostringstream.h"
#include <cstdlib>
#include <ctime>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////////////////
///             ARM_Sobol multidimentionnal sequence		 ///
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: constructor
///	Returns: 
///	Action : build the object
////////////////////////////////////////////////////
ARM_Sobol::ARM_Sobol(int firstSimulations)
:	ARM_QuasiRandom(firstSimulations),
	itsIndex(false),
	ix(), 
	iv(),
	itsFactor(0),
	itsInitRand(true)
{
	ixsobol = NULL;
	iusobol = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: DrawOne
///	Returns: 
///	Action : draw one set of values
////////////////////////////////////////////////////
void ARM_Sobol::DrawAll()
{
	newDrawAll();
	return;

	/// this routine calculates the next random number 
	/// and writes it into x
	int j, k;
	unsigned long im = itsIndex++;

	for (j=0;j< CC_NS(ARM_SOBOL,MAXBIT);j++) {
		if (!(im & 1)) break;
		im >>= 1;
	}

	/// this sequence is to do the
	/// XOR for each component
	im=j *itsDim;
	for (k=0;k< itsDim;k++) 
	{
		ix[k] ^= iv[im+k];
	
		/// we force it to be a double
		/// if for memory management we
		/// wanted to avoid this, we could use float!
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][k]=double( ix[k] ) *itsFactor;
	}
}

void ARM_Sobol::newDrawAll()
{
	int min;
	int j,k,im;

	im = insobol;

	for(j = 1;j <= MAXBITsobol; j++)
	{
		if (!(im & 1)) break;
		im >>=1;
	}

	im = (j-1)*MAXDIMsobol;

	min = itsDim < MAXDIMsobol ? itsDim : MAXDIMsobol;

	for(k = 1; k <= min; k++)
	{
		ixsobol[k]	^=	ivsobol[im+k];
		itsCurrentValues[itsCurrentBucketPos][itsCurrentPos][k-1]	= ixsobol[k]*facsobol;
	}

	insobol ++;	
}

////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: Init
///	Returns:  
///	Action : initilisation routine
////////////////////////////////////////////////////
void ARM_Sobol::Init()
{
	newInit();
	return;

	int j,k,l;
	unsigned long i,ipp, sum = 0;
	vector<unsigned long> mdeg(itsDim);
	vector<unsigned long> ip(itsDim );

	// for the initialisation 
	// we use for high dimension
	// random seed (combination of monte carlo 
	// and quasi monte carlo methods
	for(k=0;k< itsDim;k++)
	{
		mdeg[k]=CC_NS(ARM_SOBOL,MDEG)[k];
		ip[k]=CC_NS(ARM_SOBOL,IP)[k];
		for (j=0;j< (signed) mdeg[k];j++) 
			if( k+j+sum < sizeof(CC_NS(ARM_SOBOL,IV)) /sizeof( long ) )
				iv[k+j*itsDim]=CC_NS(ARM_SOBOL,IV)[k+j+sum];
			else
				iv[k+j*itsDim]= int( RandomNb() * ( 1L << j ) ) | 1;
		sum+=mdeg[k]-1;
	}

	// we compute the other direction by the recurrence
	// based on the primitive polynomials
	for (k=0;k<itsDim;k++) {
		for (j=0;j< (signed) mdeg[k];j++) 
			iv[j*itsDim+k] <<= (CC_NS(ARM_SOBOL,MAXBIT)-j-1);
		for (j=mdeg[k];j<CC_NS(ARM_SOBOL,MAXBIT);j++) {
			ipp=ip[k];
			i=iv[(j-mdeg[k])*itsDim+k];
			i ^= (i >> mdeg[k]);
			for (l=mdeg[k]-1;l>=1;l--) {
				if (ipp & 1) i ^= iv[(j-l)*itsDim+k];
					ipp >>= 1;
			}
			iv[j*itsDim+k]=i;
		}
	}

	//since the code obtained a number between 0 and 2^CC_NS(ARM_SOBOL,MAXBIT)
	// we will need to divide by 2^CC_NS(ARM_SOBOL,MAXBIT)!
	itsFactor=1.0/(1L << CC_NS(ARM_SOBOL,MAXBIT));

	// the index is number of random
	// numbers already generated
	// the draw is still not done
	itsIndex=0;
}

void ARM_Sobol::newInit()
{
	int j,k,n;
	int i,ipp;

	// allocation
	if(ixsobol != NULL) delete [] ixsobol;
	if(iusobol != NULL) delete [] iusobol;

	ixsobol = new int [MAXDIMsobol + 1];
	iusobol = new int*[MAXBITsobol + 1];

	// initialisation du générateur	
	for(i = 1; i <= MAXDIMsobol; ixsobol[i++] = 0); 
	insobol = 0;
	
	// génération aléatoire  des impairs pour initialiser le processus
	srand(9685);
	for(i = 1; i <= MAXBITsobol; i++)
	{
		for(j = 1; j <= MAXDIMsobol; j++)
		{
			ivsobol[j+(i-1)*MAXDIMsobol] = 2*(rand()%(1<<(i-1)))+1;
		}
	}
	
	// pour les plus petites dimensions (<16) on utilise les nombres de sobol issues des donnees statiques
	i=1;
	for(j = 1;j <= 16; j++)
	{
		for(k = 1; k <= mdegsobol[j]; k++)
		{
			ivsobol[j+(k-1)*MAXDIMsobol] = iwsobol[i];
		}
	}

	// retour à l'initialisation de numerical
	if(ivsobol[1] != 1) return;

	facsobol = 1.0/(1L << MAXBITsobol);

	for(j = 1,k = 0; j <= MAXBITsobol; j++, k += MAXDIMsobol)
	{
		iusobol[j] = &ivsobol[k];
	}

	for(k = 1; k <= MAXDIMsobol; k++)
	{
		for(j = 1; j <= mdegsobol[k]; j++)
		{
			iusobol[j][k] <<= (MAXBITsobol - j);
		}
		for(j = mdegsobol[k]+1; j <= MAXBITsobol; j++)
		{
			ipp	= ipsobol[k];
			i	= iusobol[j - mdegsobol[k]][k];
			i	^= (i >> mdegsobol[k]);

			for(n = mdegsobol[k]-1; n >= 1; n--)
			{
				if(ipp & 1) 
				{
					i^= iusobol[j-n][k];
				}
				ipp >>= 1;
			}
			iusobol[j][k] = i;
		}
	}

	itsIndex = 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: RandomNb
///	Returns: 
///	Action : standard random number generator
////////////////////////////////////////////////////

double ARM_Sobol::RandomNb()
{
	if( itsInitRand == false) 
	{
		itsInitRand = true;
		srand( 1 );
	}
	return rand()/(RAND_MAX+1.0);
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: SetDim
///	Returns: 
///	Action : reset the dimension!
////////////////////////////////////////////////////
void ARM_Sobol::SetDim( size_t dim )
{
	if(dim > MAXDIMsobol)
	{
		CC_Ostringstream os;
		os  << "sobol sequence with dim " << dim
			<< " but max dim for sobol is 366";
	    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());		
	}

	itsDim = dim;
	ix = vector<unsigned long>(dim,0);
	iv = vector<unsigned long>(dim*CC_NS(ARM_SOBOL,MAXBIT));
	itsInitRand = false;
	Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: reset
///	Returns: 
///	Action : reset the generator
////////////////////////////////////////////////////
void ARM_Sobol::reset( size_t dim, size_t nbOfPoints )
{
	SetDim(dim);
	itsInitRand = false;
	Init();
}





////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: Copy constructor and operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Sobol::ARM_Sobol(const ARM_Sobol& rhs)
:	ARM_QuasiRandom( rhs), 
	itsIndex(	rhs.itsIndex	),
	ix(			rhs.ix			),
	iv(			rhs.iv			),
	itsFactor(	rhs.itsFactor	),
	itsInitRand(rhs.itsInitRand )
{
	iusobol = NULL;
	ixsobol = NULL;
}

ARM_Sobol& ARM_Sobol::operator=(const ARM_Sobol& rhs )
{
	if(this!=&rhs)
	{
		ARM_QuasiRandom::operator=( rhs );
		itsIndex	= rhs.itsIndex;
		ix			= rhs.ix;		
		iv			= rhs.iv;
		itsFactor	= rhs.itsFactor;
		itsInitRand	= rhs.itsInitRand;

		if(iusobol != NULL) delete [] iusobol;
		if(ixsobol != NULL) delete [] ixsobol;

		iusobol		= NULL;
		ixsobol		= NULL;
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Sobol::~ARM_Sobol()
{
	delete [] ixsobol;
	delete [] iusobol;

	ixsobol = NULL;
	iusobol = NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_Sobol::Clone() const
{
	return new ARM_Sobol(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_Sobol
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Sobol::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "Sobol Sequence with dim " << itsDim;
	return os.str();
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


