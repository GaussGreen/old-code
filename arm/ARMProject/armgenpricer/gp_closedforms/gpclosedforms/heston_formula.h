/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=theta*(lgtvol^2-V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///			V(0)=sig^2
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
//////
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_HESTON_FORMULA_H
#define _GP_CF_HESTON_FORMULA_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result


CC_BEGIN_NAMESPACE(ARM)

/// forward declaration
class ArgumentList;

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Heston_JumpDiffusion_Formula
///			Purpose : Evaluation of avanilla option in the Heston with jumps model
///			Assumptions: lognormal hypothesis with stochastic volatility and normal jumps of the log
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Heston_JumpDiffusion_Formula
{
	enum ArgumentType
	{
		INDEX = 0,					// i=0
		STRIKE,						// i=1
		INITIALVOL,			// i=2
		TIMETOMATURITY,				// i=3
		LONGTERMVOL,				// i=4  omega
		VOLSQUAREREVERTINGSPEED,	// i=5  theta
		VOLSQUAREVOLATILITY,		// i=6  ksi
		CORRELATION,				// i=7  rho
		JUMPPROBABILITY,			// i=8	lambda	
		JUMPMEAN,					// i=9  muJ
		JUMPVOLATILITY,				// i=10  sigmaJ
		CALLORPUT,					// i=11
		NBTERMS						// i=12
	};
	enum Vector_Interpolation_Method
	{
		LINEARINTERPOLATION=0,
		PARABOLICINTERPOLATION
	};
	enum 
	{ 
		Nb_Parameters =13
	};
	enum
	{ 
		Nb_Derivable_Parameters =11
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};


///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_Heston_Formula
///			Purpose : Evaluation of avanilla option in the Heston with jumps model
///			Assumptions: lognormal hypothesis with stochastic volatility 
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_Heston_Formula
{
	enum ArgumentType
	{
		INDEX = 0,					// i=0
		STRIKE,						// i=1
		INITIALVOLSQUARE,					// i=2
		TIMETOMATURITY,				// i=3
		LONGTERMVOLSQUARE,				// i=4
		VOLSQUAREREVERTINGSPEED,	// i=5
		VOLSQUAREVOLATILITY,		// i=6
		CORRELATION,				// i=7
		CALLORPUT,					// i=8
		NBTERMS						// i=9
	};
	enum Vector_Interpolation_Method
	{
		LINEARINTERPOLATION=0,
		PARABOLICINTERPOLATION
	};
	enum 
	{ 
		Nb_Parameters =10
	};
	enum
	{ 
		Nb_Derivable_Parameters =8
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};



///////////////////////////////////////////////////////////////////////
///  
///			class : FFT
///			Purpose : Evaluation of a vanilla option with FFT
///			Assumptions: lognormal hypothesis with stochastic volatility 
///
///////////////////////////////////////////////////////////////////////

class ARM_FFT
{

public:	
	enum FFTMethodType
	{
		forward =   +1,
		inverse =   -1,
	};
	typedef ARM_FFT::FFTMethodType ARM_FFTType;
	
	ARM_FFT(const int n, const ARM_FFTType& type) 
		:
		itsNSize(n),
		itsFFTType(type),
		itsRealData()
	{
	}

	~ARM_FFT() {}

	//! \brief computeFFT interfaces the fft algorithm
	/*! 
	 *
	 *  \author vanhoovecl
	 *
	 *  \param a is a vector of real data
	 *  \param dir if is equal to 1 computes fft, if is equal to -1 computes inverse fft.
	 *
	 *  \return the vector of complex fft
	 */
	std::vector<double> computeFFT(const std::vector<double>& realData, int dir);
private :

	int itsNSize; /// should be power of 2, otherwise, it will be resized

	ARM_FFTType  itsFFTType;

	std::vector<double> itsRealData;
	//! \brief completeToNextPowerOfTwo completes a vector with 0 until his length reachs a power of 2.
	 /*  \author vanhoovecl
	 *  \param data is a vector of data on which FFT will be applied
	 *  \return a vector containing values of <code>data</code> and 0. Its length is a power of 2. 
	 */
	std::vector<double> completeToNextPowerOfTwo(const std::vector<double>& realData);

	//! \brief fromRealDataToComplexData transforms a real vector into a complex vector
	 /*  \author vanhoovecl
	 *  \param realData is a vector of real data
	 *  \return a vector of complex data such as real part and imaginary part alternates, 
	 *		full filled with <code>realData</code> values.
	 */
	std::vector<double> fromRealDataToComplexData(const std::vector<double>& realData);

	//! \brief complexTransform computes FFT (Danielson-Lanczos)
	/*  One-dimensional discrete complex fourier transform
	*  Replaces a[ 0...2*len ] by its discrete Fourier transform.
	*  In the INVERSE operation, a gain normalisierung by 1/<code>len</code> is
	*  applied automatically.
	*  <p>
	*  The routine (Danielson-Lanczos) was adapted from 'Numerical Recipes in C',
	*  optimized, and given readable variable names.
	*
	*  \author vanhoovecl
	*  \param  a       complex vector with real part in
	*                  <code>a</code>[ 0, 2, 4, ... 2*len - 2 ],
	*                  imaginary part in <code>a</code>[ 1, 3, ... 2 * len -1 ]
	*  \param  len     MUST be an integer power of 2
	*  \param  dir     use <code>INVERSE</code> or <code>FORWARD</code>
	*/
	static void complexTransform( std::vector<double> realData, int len, int dir );

};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
