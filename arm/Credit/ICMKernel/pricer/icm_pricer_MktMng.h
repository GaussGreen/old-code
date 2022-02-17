
#ifndef _ICM_Pricer_MKTMNG_H
#define _ICM_Pricer_MKTMNG_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/crv/icm_distribloss.h"
#include "ICMKernel/inst/icm_cds.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_MktMng ICM_Pricer_MktMng.h "ICM_Pricer_MktMng.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004 
 *	\file   ICM_Pricer_MktMng.h
 *	\brief  Pricer for all security using ICM_Security */
/***********************************************************************************/


class ICM_Pricer_MktMng : public ICM_Pricer
{
private:

public :
	ICM_Pricer_MktMng(void) 
	{ Init(); }

	~ICM_Pricer_MktMng()
	{
	}

	void Init()
	{
		//SetName(ICM_PRICERSECURITY);
	}

	void Set(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters&,const ARM_Date&);


public:
	// Pricer MktMng Specifics

	virtual void Reset(void)
	{}

	void View(char* id = NULL, FILE* ficOut = NULL);
protected :
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
											 const std::string& plot, 
											 const std::string& labelKey, 
											 double epsvalue, double epsvalueGamma = 0);
	

private:
	//	Unimplemented
	ICM_Pricer_MktMng(const ICM_Pricer_MktMng&ref); 
	ICM_Pricer_MktMng& operator=(const ICM_Pricer_MktMng&); 
	void BeforePrice(const std::string& label,qSENSITIVITY_TYPE /** type = ICMSPREAD_TYPE **/ );
	//	Throw
};


#endif

