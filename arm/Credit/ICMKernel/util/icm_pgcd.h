
#if !defined(ICM_PGCD_CLASS)
#define ICM_PGCD_CLASS

#include <string>
#include <vector>
#include "ARMKernel\glob\linalg.h"


/*********************************************************************************/
/*! \class  pgcd icm_pgcd.h "icm_pgcd.h"
 *  \author F. Dezormes
 *	\version 1.0
 *	\date   March 2004
 *	\file   icm_pgcd.h
 *	\brief for pgcd computing 
**********************************************************************************/

class pgcd  
{
public:
	pgcd();
	virtual ~pgcd();
	static long		solve(long a, long b);
	static long		solve (const std::vector<long>& a);
	static double	solve (const std::vector<double>& a);
	static double	solve (const ARM_Vector& a);
	static double	round(double value);

	static bool test(std::string& errStr);
};

#endif // !defined(ICM_PGCD_CLASS)
