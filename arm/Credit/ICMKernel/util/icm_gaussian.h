#ifndef _ICM_GAUSSIAN_H_
#define _ICM_GAUSSIAN_H_

// ----------------------------------------------------------------------------------
//	class		ICM_gaussian
//  author		L. Jacquel
//	version		1.0
//	date		September 2004
//	file		ICM_gaussian.h

//	brief		Useful functions for Standard Gaussian handling
// ----------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------
// PI
# ifndef PI
# define PI					3.14159265358979
# endif
# define TWOPI				6.28318530717959
# define SQRT2PI			2.50662827463100
# define ONEOVERSQRT2PI		0.39894228040143267793994605993438

// ----------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------
// for NCum

#define NORCUMB1	0.319381530
#define NORCUMB2	-0.356563782
#define NORCUMB3	1.781477937
#define NORCUMB4	-1.821255978
#define NORCUMB5	1.330274429
#define NORCUMA		0.2316419

// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// for InvNCum

#define INVNORCUMA0		2.50662823884												
#define INVNORCUMA1		-18.61500062529
#define INVNORCUMA2		41.39119773534
#define INVNORCUMA3		-25.44106049637
#define INVNORCUMB1		 -8.47351093090
#define INVNORCUMB2		23.08336743743
#define INVNORCUMB3		-21.06224101826 
#define INVNORCUMB4		3.13082909833
#define INVNORCUMC0		0.3374754822726147
#define INVNORCUMC1		0.9761690190917186
#define INVNORCUMC2	 	0.1607979714918209
#define INVNORCUMC3		0.0276438810333863
#define INVNORCUMC4	 	0.0038405729373609
#define INVNORCUMC5	 	0.0003951896511919
#define INVNORCUMC6		0.0000321767881768
#define INVNORCUMC7		0.0000002888167364
#define INVNORCUMC8		0.0000003960315187
#define INVNORCUMEPS	1.0e-9
#define INVNORCUMDEF	0.0
// ----------------------------------------------------------------------------------


	// Standard Gaussian Cumulative Function
	double NCum(double	x);

	// Standard Inverse Gaussian Cumulative Function
	double InvNCum(double InValue);

	// Standard Gaussian density
	double StdNorDensity(double x);				// with /Rac(2 Pi)
	double ExpMinusHalfXSquare(double x);		// exp(-x2/2)

//---------------------- binormal cumulée----------------------------------------------------------------------------------

double IntermFunction(double x,
				      double y);

double BinormalCumuleNegativeorzerocase(double a,
										  double b,
										  double rho);


double BinormalGenericCase(double a,
						   double b,
						   double rho);

double BinormalCumule(double a,
					  double b,
					  double rho);
#endif