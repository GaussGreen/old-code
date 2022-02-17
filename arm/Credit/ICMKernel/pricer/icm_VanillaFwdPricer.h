
#pragma warning( disable : 4150 ) 
#ifndef _ICM_VANILLA_FWD_PRICER_H
#define _ICM_VANILLA_FWD_PRICER_H

#include "ICMKernel/pricer/icm_pricer_basket.h"
#include "ICMKernel/inst/icm_ftd.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Distrib icm_pricer_homogeneous.h "icm_pricer_homogeneous.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  ce pricer est basé sur l'utilisation d'un cache de distribution */
/***********************************************************************************/

class ICM_VanillaFwdPricer : public ICM_Pricer_Basket
{

private :

	ICM_Pricer*		itsFwdStartPricer;	//Pricer for (itsFwdStartSec,Model)
	ICM_Pricer*		itsFwdEndPricer;	//Pricer for (itsFwdEndSec,Model)

	ICM_Security*	itsFwdStartSec;		//Security(AsOf,Start)
	ICM_Security*	itsFwdEndSec;		//Security(AsOf,Maturity)

	ICM_Correlation*	itsFwdStartCorr;//Correlation linked to itsFwdStartSec
	ICM_Correlation*	itsFwdEndCorr;	//Correlation linked to itsFwdEndSec
	bool				itsIsSensiLinkedtoCorr;
	bool				itsForcedRescaling; //ForcedRescaling

public :
	ICM_VanillaFwdPricer(void) { Init();}

	ICM_VanillaFwdPricer(ICM_Pricer *Initpricer)
	{
		Init();
		Set(Initpricer);
	}

	void Init() 
	{
		SetName(ICM_PRICER_VANILLA_FWD);

		itsFwdStartPricer=NULL;
		itsFwdEndPricer=NULL;

		itsFwdStartSec=NULL;
		itsFwdEndSec=NULL;

		itsFwdStartCorr=NULL;
		itsFwdEndCorr=NULL;
		itsIsSensiLinkedtoCorr = false;
		itsForcedRescaling = false;
	}

	~ICM_VanillaFwdPricer(void)	
	{
	if (itsFwdStartPricer) delete itsFwdStartPricer;itsFwdStartPricer=NULL;
	if (itsFwdEndPricer) delete itsFwdEndPricer;itsFwdEndPricer=NULL;
	if (itsFwdStartSec) delete itsFwdStartSec;itsFwdStartSec=NULL;
	if (itsFwdEndSec) delete itsFwdEndSec;itsFwdEndSec=NULL;
	if (itsFwdStartCorr) delete itsFwdStartCorr;itsFwdStartCorr=NULL;
	if (itsFwdEndCorr) delete itsFwdEndCorr;itsFwdEndCorr=NULL;
	}

	void Set(ICM_Pricer *Initpricer);

	virtual double ExpectedLossTranche(const double& yearterm,vector<double>& losses){return 0.;}

	virtual double ComputePrice(qCMPMETH measure);



	inline void SetFwdStartPricer(ICM_Pricer* pricer) {itsFwdStartPricer = pricer;}
	inline void SetFwdEndPricer(ICM_Pricer* pricer) {itsFwdEndPricer = pricer;}
	inline void SetFwdStartSec(ICM_Security* sec) {itsFwdStartSec = sec;}
	inline void SetFwdEndSec(ICM_Security* sec) {itsFwdEndSec = sec;}
	inline void SetFwdStartCorr(ICM_Correlation* corr) {itsFwdStartCorr = corr;}
	inline void SetFwdEndCorr(ICM_Correlation* corr) {itsFwdEndCorr = corr;}

	inline ICM_Pricer* GetFwdStartPricer() {return itsFwdStartPricer;}
	inline ICM_Pricer* GetFwdEndPricer() {return itsFwdEndPricer;}
	inline ICM_Security* GetFwdStartSec() {return itsFwdStartSec;}
	inline ICM_Security* GetFwdEndSec() {return itsFwdEndSec;}
	inline ICM_Correlation* GetFwdStartCorr() {return itsFwdStartCorr;}
	inline ICM_Correlation* GetFwdEndCorr() {return itsFwdEndCorr;}

	virtual void BeforePrice(const std::string& label,qSENSITIVITY_TYPE type ) 
	{
		itsIsSensiLinkedtoCorr = false;
		itsForcedRescaling  = false;

		if ((type==ICMCORRELATION_TYPE)||
			(type==ICMCORREL_STRIKE_DOWN_TYPE)||
			(type==ICMCORREL_STRIKE_UP_TYPE))
		{itsIsSensiLinkedtoCorr = true;}	

		if (type==ICM_INDX_SPREAD_RESCALING)
		{itsForcedRescaling = true;}
	} 

	virtual void PropagateModel(ARM_Model *model) 
    {
        itsFwdStartPricer->SetModel(model);
		itsFwdEndPricer->SetModel(model);
    }

	virtual void Reset(void)
	{
		ICM_Pricer_Basket::Reset();
		itsFwdStartPricer->ResetPricer();
		itsFwdEndPricer->ResetPricer();
	}

};

#endif

