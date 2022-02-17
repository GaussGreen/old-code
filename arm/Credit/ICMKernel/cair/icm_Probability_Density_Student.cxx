#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Probability_Density_Student.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for Student Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICM_Probability_Density_Student.h"


#include "ICMKernel\glob\icm_maths.h"



void ICM_Probability_Density_Student :: Init()
{
}

ICM_Probability_Density* ICM_Probability_Density_Student :: Clone() const
{
	ICM_Probability_Density_Student *New;
	New = new ICM_Probability_Density_Student(*this);
	return New;
}

void ICM_Probability_Density_Student :: Copy(const ICM_Probability_Density& Data)
{
	degree		=	((ICM_Probability_Density_Student&)Data).degree;
}


void ICM_Probability_Density_Student :: SetParameters(double data_degree)
{
	degree		=	data_degree;
}


		
// ----------------------------------------------------------- 
// Density Function
// ----------------------------------------------------------- 

double	ICM_Probability_Density_Student :: Density_Function(double x)
{
	return StudentDensity(degree, x);
}


// -----------------------------------------------------------
// Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Student :: Cumulative_Density_Function(double x)
{
	// nag_prob_students_t
	return	NAG_prob_students_t( x, degree );
}


// -----------------------------------------------------------
// Inverse Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Student :: Inverse_Cumulative_Density_Function(double x)
{
	return	NAG_deviates_students_t( x, degree );
}
	
