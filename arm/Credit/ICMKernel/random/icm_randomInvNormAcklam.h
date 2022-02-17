#ifndef _ICM_RANDOM_INVNORMACK_H_
#define _ICM_RANDOM_INVNORMACK_H_

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\random\icm_RandomLaws.h"

/* from 
 ep::MathSrv::cumNorm*/


class ICM_RandomInvNormAcklam : public ICM_RandomLaws
{
private :
	// unsetted : CONST
	double its_a1;
	double its_a2;
	double its_a3;
	double its_a4;
	double its_a5;
	double its_a6;

	double its_b1;
	double its_b2;
	double its_b3;
	double its_b4;
	double its_b5;

	double its_c1;
	double its_c2;
	double its_c3;
	double its_c4;
	double its_c5;
	double its_c6;

	double its_d1;
	double its_d2;
	double its_d3;
	double its_d4;

	double its_u_low;
	double its_u_high;

protected :
	void Init();
	
public : 

	ICM_RandomInvNormAcklam(); 
	ICM_RandomInvNormAcklam(const ICM_RandomGenerator& pRandomUnif);
	ICM_RandomInvNormAcklam(const ICM_RandomInvNormAcklam& ref); 
	~ICM_RandomInvNormAcklam();
	ICM_RandomInvNormAcklam& operator=(const ICM_RandomInvNormAcklam& ref);
	ARM_Object * Clone(void);

	// Pure virtual functions to be defined in derived classes : 

	// function to generate  a vector
	void GenerateRandomInvNormAcklam(ARM_Vector& RandomVector) const ;
	/// function to generate a vector
	virtual void GenerateRandoms(ARM_Vector& RandomVector) ;
};

inline void ICM_RandomInvNormAcklam::GenerateRandoms(ARM_Vector& RandomVector) 
{	
	ICM_RandomLaws::GenerateRandoms(RandomVector);
	GenerateRandomInvNormAcklam(RandomVector);
}



inline void ICM_RandomInvNormAcklam::GenerateRandomInvNormAcklam(ARM_Vector& RandomVector) const 
{
	double * pRandVect = RandomVector.GetElt();
	//---------------------------------------------------------------------------------------------
	// The Inverse cumulative normal distribution function
	// Peter J. Acklam Method.
	// URL: http://www.math.uio.no/~jacklam/notes/invnorm


	  //msgTools::assert ( fabs(u-0.5)>=0.5, "u should belong to (0,1)");
	  double z=0, r=0;

	  for ( int i=0; i<RandomVector.size(); i++)
	  {
		  // Rational approximation for the lower region. ( 0 < u < its_u_low )
		  if( pRandVect[i] < its_u_low ){
			z = sqrt(-2.0*log(pRandVect[i]));
			z = (((((its_c1*z+its_c2)*z+its_c3)*z+its_c4)*z+its_c5)*z+its_c6) / ((((its_d1*z+its_d2)*z+its_d3)*z+its_d4)*z+1.0);
		  }
  
		  // Rational approximation for the central region. ( its_u_low <= u <= its_u_high )
		  else if( pRandVect[i] <= its_u_high ){
			z = pRandVect[i] - 0.5;
			r = z*z;
			z = (((((its_a1*r+its_a2)*r+its_a3)*r+its_a4)*r+its_a5)*r+its_a6)*z / (((((its_b1*r+its_b2)*r+its_b3)*r+its_b4)*r+its_b5)*r+1.0);
		  }
  
		  // Rational approximation for the upper region. ( its_u_high < u < 1 )
		  else {
			z = sqrt(-2.0*log(1.0-pRandVect[i]));
			z = -(((((its_c1*z+its_c2)*z+its_c3)*z+its_c4)*z+its_c5)*z+its_c6) /  ((((its_d1*z+its_d2)*z+its_d3)*z+its_d4)*z+1.0);
		  }

		  // The relative error of the approximation has absolute value less
		  // than 1.15e-9.  One iteration of Halley's rational method (third
		  // order) gives full machine precision.

		  //r = (cumNorm(z) - u) * sqrt(2*PI) * exp( 0.5 * z * z );	//	f(z)/df(z)
		  //z -= r/(1+0.5*z*r);							//	Halley's method
		   pRandVect[i] = z;
	  }
}

#endif 