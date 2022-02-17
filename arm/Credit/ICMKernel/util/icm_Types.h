/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		TYPES.H
	PROJECT:	KERNEL
	
	DESCRIPTION:	list of personnalized types

   -----------------------------------------------------------------
   
	CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */


# ifndef _TYPES
# define _TYPES

# include "enums.h"


	typedef		unsigned short		ReturnCode;

	typedef		unsigned long 		Index;
	typedef		double				RelativeDate;


	typedef		double	Real;
	typedef		int		Integer;
	typedef		bool	Boolean;


	
	typedef int Day;
	typedef int Month;
	typedef int Year;
/*
	typedef struct 
	{
		Day		the_day;
		Month	the_month;
		Year	the_year;
	}		CLDMY;

/*
	typedef struct 
	{
		Day		the_day;
		Month	the_month;
		Year	the_year;
	}		CAIR_DMY;

	typedef		CAIR_DMY	CLDMY;
*/
	// Frequency = Duration
	struct Frequency
	{
		long		Nb;
		TimeUnit	U;
	};
	
	typedef Frequency Duration;


// TO DEAL WITH
	
	extern bool IsLargerThan(Duration D1, Duration D2);

# endif