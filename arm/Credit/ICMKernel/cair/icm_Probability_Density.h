/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_PROBABILITY_DENSITY.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _ICM_PROBABILITY_DENSITY_H_
# define _ICM_PROBABILITY_DENSITY_H_

#include "ICMKernel\cair\types.h"

class ICM_Probability_Density
{
	// -----------------------------------------------------------
	// CONSTRUCTORS AND DESTRUCTORS
	// -----------------------------------------------------------

	public:

		// constructors
		ICM_Probability_Density() {Init();}
		ICM_Probability_Density(const ICM_Probability_Density& Data) {Init(); Copy(Data);}
		
		// destructor
		~ICM_Probability_Density()	{Reset();}

	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		virtual ICM_Probability_Density*	Clone() const;
		virtual void	Destroy() {delete this;};
		virtual void	Copy(const ICM_Probability_Density& Data);
		virtual void	Reset() {}; 

	private:

		void Init() {};

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// CLASS USE
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		// -----------------------------------------------------------
		// Function Returning Values
		// ----------------------------------------------------------- 
		
		// Density Function
		virtual double	Density_Function(double x) = 0;

		// Cumulative Density Function
		virtual double	Cumulative_Density_Function(double x) = 0;

		// Inverse Cumulative Density Function
		virtual double	Inverse_Cumulative_Density_Function(double x) = 0;

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// DATA
	// -----------------------------------------------------------
	// -----------------------------------------------------------

		
};

# endif