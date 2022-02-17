#ifndef _ICM_PRICERMC_H_
#define _ICM_PRICERMC_H_

#include "ICMKernel/pricer/icm_pricer_security.h"
#include "ICMKernel/util/icm_root_generator.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_MC icm_pricer_mc.h "icm_pricer_mc.h"
 *  \author Ruben Marciano 
 *	\version 1.0
 *	\date   April 2004
 *	\brief  generate random variables or default times  */
/***********************************************************************************/

class ICM_Pricer_MC : public ICM_Pricer_Security
{
	
	private :

	ICM_Root_Generator * itsGenerator ;
	int	itsNbPaths ;

public :

	ICM_Pricer_MC() {Init();}

	void Init(void) 
	{
		SetName(ICM_MCPRICER);
		itsGenerator=NULL;
		itsNbPaths=-1;
	};

	void Set(ARM_Security * sec, 
			 ARM_Model * mod, 
			 const ICM_Parameters& params,
			 const ARM_Date& asof,
			 int nbpaths
			 ) 
	{
		ICM_Pricer_Security::Set(sec,mod,params,asof);

		SetNbPaths(nbpaths);
	}

	ICM_Pricer_MC(ARM_Security * sec, 
				  ARM_Model * mod, 
				  const ICM_Parameters&params,
				  const ARM_Date&asof,
				  int nbpaths) 
	{
		Init();
		Set(sec,mod,params,asof,nbpaths);
	}

	~ICM_Pricer_MC()
	{
		if (itsGenerator)
			delete itsGenerator;
		itsGenerator = NULL;
	}

	virtual double Random(unsigned int i) { return itsGenerator->generateRandom(i);}
	
	void SetGenerator(ICM_Root_Generator * Generator) 
	{ 
		if (itsGenerator)
			delete itsGenerator;
		itsGenerator = Generator; 
	}
	
	ICM_Root_Generator * GetGenerator(void) { return itsGenerator; }

	void SetNbPaths(int NbPaths) { itsNbPaths = NbPaths;}
	
	int GetNbPaths(void) { return itsNbPaths;}

	void ResetPricer(void)  
	{
		ICM_Pricer_Security::ResetPricer();

		if (itsGenerator) itsGenerator->Reset();
	}

} ;

#endif
