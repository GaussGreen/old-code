#ifndef _FIRSTOBEINCLUDED_H
#define _FIRSTOBEINCLUDED_H

/*
 * File to be included EVERYWHERE
 * disable STL warnings in Windows NT
 * to get back STL warning use the flag _SHOW_STL_WARNING
 */

#ifdef WIN32
	#ifndef _SHOW_STL_WARNING
		#pragma warning(disable : 4786)
		#pragma warning(disable : 4503)
		#pragma warning(disable : 4275)
	    #pragma warning(disable : 4006)
		#pragma warning(disable : 4005)	// macro redefinition
	#endif
#endif


struct ARM_CF_StochVolDispatcher
{
	enum DistributionType
	{
		K_STOCHASTIC_BLACKSCHOLES_DIST,
		K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST,
		K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST,
		K_SABR_DIST,
		K_SABR_MC_DIST,
	};
	
	enum GreekType
	{
		K_Delta,
		K_Vega,
		K_Gamma
	};
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
