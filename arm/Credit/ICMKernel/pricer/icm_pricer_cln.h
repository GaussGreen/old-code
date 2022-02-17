#ifndef _ICM_Pricer_CLN_H
#define _ICM_Pricer_CLN_H

#include "ICMKernel/pricer/icm_pricer_cds.h"
#include "ICMKernel/inst/icm_cln.h"
#include "ICMKernel/util/icm_polynoms.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Cln icm_pricer_cln.h "icm_pricer_cln.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2004 
 *	\file   icm_pricer_cds.h
 *	\brief  Credit Linked Note pricer */
/***********************************************************************************/


class ICM_Pricer_Cln : public ICM_Pricer_Cds
{

private :

public :
	ICM_Pricer_Cln(void) ; 

	/** ICM_Pricer_Cln(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&param,const ARM_Date&) ;
	**/ 

	~ICM_Pricer_Cln() ; 

	void Init()	; 

	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&param,const ARM_Date&) ;

private:
	virtual double FeeLegPV () ; 
	virtual double Accrued() ;
	virtual double DefLegPV () ; 
public:
	virtual void Reset(void); 
	

};

#endif

