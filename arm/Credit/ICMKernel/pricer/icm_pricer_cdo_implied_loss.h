
#ifndef _ICM_PRICER_CDO_IMPLIED_LOSS_H
#define _ICM_PRICER_CDO_IMPLIED_LOSS_H

#include "ICMKernel\pricer\icm_pricer.h"
// #include "ICMKernel\crv\icm_implied_loss_tree.h"
//#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\crv\icm_distribloss.h"
//#include "ICMKernel/glob/icm_mktdatamng.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_ImpliedLoss icm_pricer_cdo_implied_loss.h "icm_pricer_cdo_implied_loss.h"
 *  \author Damien Pouponneau
 *	\version 1.0
 *	\date   June 2006
 *	\brief  Pricer for Credit Tree  */
/***********************************************************************************/

class ICM_Pricer_Security; 
class ICM_ImpLossTree; 

class ICM_Pricer_ImpliedLoss : public ICM_Pricer
{
private:
	ICM_ImpLossTree*		itsTree;
	vector<double>			itsYearTerms;
	ICM_DistribLoss			itsDistribLoss;
	ICM_Pricer_Security*	itsPricer;

public :
	ICM_Pricer_ImpliedLoss(void) { Init();}

	ICM_Pricer_ImpliedLoss(ARM_Security* sec, 
						   ARM_Object* datamkt,
						   const ICM_Parameters& parameters,
						   const ARM_Date&asof)
	{
		Init();
		Set(sec, datamkt,parameters,asof);
	}

	~ICM_Pricer_ImpliedLoss() ;

	void Init()
	{
		SetName(ICM_IMPLIED_LOSS_TREE_PRICER);
		itsTree=NULL;
		itsDistribLoss.clear();
		itsYearTerms.clear();
		itsPricer = NULL;
	}

	void Set(ARM_Security *option, 
			 ARM_Object *mod,
			 const ICM_Parameters& parameters ,const ARM_Date&asof) ;
/**	{
		ICM_Pricer::Set(option,mod,parameters,&asof);
		itsTree = (ICM_ImpLossTree*)((ICM_MktDataMng*)mod)->GetLossTree("USERDEF",asof);
	}**/ 

	virtual void Reset(void)
	{
		ICM_Pricer::Reset();
	}



	virtual double ComputePrice(qCMPMETH measure) ;

	void CptFwdLoss(int idxstart,
					int idxend,
					int dw_loss ,
					int up_loss,
					double notional_dw,
					double notional_up);

	void CptFwdLoss(double dw_loss ,double up_loss);

	void View(char* id, FILE* ficOut);
protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon,double epsilonGamma =0 ) ; 

private:
	double Accrued(void) ; 
	double FeeLegPV(void) ;
	double DefLegPV(void) ;
	double Compute_Fwd_Spread(const  ARM_Date &,const  ARM_Date &, double& dur); 
	double ComputeDuration(void) ;
	double ComputeSpread(const double &) ;
	double ComputeImpliedVol(const double &) ;

};

#endif

