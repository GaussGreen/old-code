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



#include <vector>
#include <string>
#include <map>

// 17783 using namespace std;

#include "enums.h"
#include "ICMKernel\glob\icm_constants.h"


	typedef		unsigned short		ReturnCode;

	typedef		unsigned long 		IndexId;
	typedef		double				RelativeDate;

	typedef		double	Double;
	// defined in NAG
//	typedef		int		Integer;
//	typedef		bool	Boolean;
	
	typedef int Day;
	typedef int Month;
	typedef int Year;

	typedef struct 
	{
		Day		the_day;
		Month	the_month;
		Year	the_year;
	}		DateDMY;

	typedef std::vector<int>				IntVector;
	typedef std::vector<double>			DoubleVector;
	typedef std::vector<RelativeDate>	DateVector;
	typedef std::vector<std::string>			StringVector;
	typedef std::vector<double*>			DoubleMatrix;


	// Kind of map

class Loss_Key
{
	public:

		double	Time_T;		// time T
		double	K_Min;		// Strike Min
		double	K_Max;		// Strike Max
		
	public:

		Loss_Key()
		{
			Init();
		}

		Loss_Key(const Loss_Key& data)
		{
			Time_T	=	data.Time_T;
			K_Min	=	data.K_Min;
			K_Max	=	data.K_Max;
		}

		Loss_Key(double d_time, double d_k_min = 0.0, double d_k_max = 0.0)
		{
			Time_T	=	d_time;
			K_Min	=	d_k_min;
			K_Max	=	d_k_max;
		}

		~Loss_Key(void)
		{
		}

	public:

		bool operator < (const Loss_Key & rhs) const {
			if (Time_T < rhs.Time_T)
				return true;
			else
				if (Time_T > rhs.Time_T)
					return false;
				else
					if (K_Min < rhs.K_Min)
						return true;
					else
						if (K_Min > rhs.K_Min)
							return false;
						else
							return (K_Max < rhs.K_Max);

		}
		
		bool operator == (const Loss_Key & rhs) const {
			return ((Time_T == rhs.Time_T) && (K_Min == rhs.K_Min) && (K_Max == rhs.K_Max));
		}
		
	protected:

		void Init()
		{
			Time_T	=	0.0;
			K_Min	=	0.0;
			K_Max	=	0.0;

		}
};


typedef	std::map<Loss_Key, double> CreditManager_Map_Loss;


# endif