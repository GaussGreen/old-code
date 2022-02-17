/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CONSTANTS.H
	PROJECT:	GLOB
	
	DESCRIPTION:	list of personnalized constants

  -----------------------------------------------------------------

 	CREATION:	October 8, 2004

	LAST MODIF:	October 8, 2004
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _ICM_CONSTANTS_H_
# define _ICM_CONSTANTS_H_

	# define	INFINITY_MINUS	-1e25
	# define	INFINITY_PLUS	1e+25

	// Code Error

	# define	RetNotOk	0
	# define	RetOk		1

	// Maths Constants
	# define ONEOVER365			2.73972602739726027397260273972e-3
	# define ONEOVER365SQUARE	7.50609870519797335334959654719e-6
//	# define PI					3.14159265358979
	# define TWOPI				6.28318530717959
	# define SQRT2PI			2.50662827463100
	# define ONEOVERSQRTPI		0.56418958354775628694807945156077
	# define ONEOVERPI			0.31830988618379

	# define SQRT2				1.4142135623730950488016887242097
	
	// Str Constant
	// # define AsOf_Str			2451545		// Correspond au julian du 1 janvier 2000.

	// String constants
	#define DayKeyWord			"D"
	#define BusinessDayKeyWord	"BD"
	#define WeekKeyWord			"W"
	#define MonthKeyWord		"M"
	#define YearKeyWord			"Y"


# endif
