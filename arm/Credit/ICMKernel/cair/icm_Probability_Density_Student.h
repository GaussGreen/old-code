/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_PROBABILITY_DENSITY_STUDENT.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for a Student Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef _ICM_PROBABILITY_DENSITY_STUDENT_H_
#define _ICM_PROBABILITY_DENSITY_STUDENT_H_

# include "icm_probability_density.h"


class ICM_Probability_Density_Student : public ICM_Probability_Density
{
	// -----------------------------------------------------------
	// CONSTRUCTORS AND DESTRUCTORS
	// -----------------------------------------------------------

	public:

		// constructors
		ICM_Probability_Density_Student() {Init();}
		ICM_Probability_Density_Student(double degree) {Init(); SetParameters(degree);}
		ICM_Probability_Density_Student(const ICM_Probability_Density_Student& Data) {Init(); Copy(Data);}
		
		// destructor
		~ICM_Probability_Density_Student()	{Reset(); }

	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		virtual ICM_Probability_Density*	Clone() const;
		virtual void	Destroy() {delete this;};
		virtual void	Copy(const ICM_Probability_Density& Data);
		virtual void	Reset() {}; 

	private:
	
		void Init();
		void SetParameters(double data_degree);

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// CLASS USE
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public: 

		double	GetDegree() {return degree;}
		void	SetDegree(double data) {degree = data;}

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
		double	degree;

};

#endif
	
