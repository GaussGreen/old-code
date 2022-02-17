/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_PROBABILITY_DENSITY_GAUSSIAN.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for a Gaussian Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef _ICM_PROBABILITY_DENSITY_GAUSSIAN_H_
#define _ICM_PROBABILITY_DENSITY_GAUSSIAN_H_

# include "ICMKernel\cair\icm_probability_density.h"

class ICM_Probability_Density_Gaussian : public ICM_Probability_Density
{
	// -----------------------------------------------------------
	// CONSTRUCTORS AND DESTRUCTORS
	// -----------------------------------------------------------

	public:

		// constructors
		ICM_Probability_Density_Gaussian() {Init();}
		ICM_Probability_Density_Gaussian(double mean, double varaiance) {Init(); SetParameters(mean, variance);}
		ICM_Probability_Density_Gaussian(const ICM_Probability_Density_Gaussian& Data) {Init(); Copy(Data);}
		
		// destructor
		~ICM_Probability_Density_Gaussian()	{Reset(); }

	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		virtual ICM_Probability_Density*	Clone() const;
		virtual void	Destroy() {delete this;};
		virtual void	Copy(const ICM_Probability_Density& Data);
		virtual void	Reset() {}; 

	private:
	
		void Init();
		void SetParameters(double data1, double data2);

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// CLASS USE
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public: 

		double	GetMean() {return mean;}
		void	SetMean(double data) {mean = data;}

		double	GetVariance() {return variance;}
		void	SetVariance(double data) {variance = data;}

	public : 
		
		// -----------------------------------------------------------
		// Function Returning Values
		// ----------------------------------------------------------- 
		
		// Density Function
		double	Density_Function(double x);

		// Cumulative Density Function
		double	Cumulative_Density_Function(double x);

		// Inverse Cumulative Density Function
		double	Inverse_Cumulative_Density_Function(double x);
	
		// -----------------------------------------------------------
		// DATA
		// ----------------------------------------------------------- 

	protected:
		double	mean;
		double	variance;

};

#endif
	
