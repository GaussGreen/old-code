
#ifndef _ICM_Pricer_Corridor_H
#define _ICM_Pricer_Corridor_H

#pragma warning (disable:4786)

#include "ARMKernel/crv/zerocurv.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/pricer/icm_pricer_MktMng.h"


class ICM_Fixing_Curve;
class ICM_Correlation;
class ICM_CorridorLeg;
class ICM_MktDataMng;

/*********************************************************************************/
/*! \class  ICM_Pricer_Corridor ICM_Pricer_Corridor.h "ICM_Pricer_Corridor.h"
 *  \author Laurent Jacquel
 *	\version 1.0
 *	\date   April 2005
 *	\file   ICM_Pricer_Corridor.h
 *	\brief  Pricer for Corridor */
/***********************************************************************************/


class ICM_Pricer_Corridor : public ICM_Pricer_MktMng
{
private:

	ARM_ZeroCurve*		itsZeroCurve;
	ICM_DefaultCurve*	itsDefCurveIndex;
	ARM_VolCurve*		itsIRVolatility;
	ARM_VolCurve*		itsIdxVolatility;
	ICM_Correlation*	itsCorrelation;
	ICM_Fixing_Curve*	itsFixingCurveCredit;
	ICM_Fixing_Curve*	itsFixingCurveIR;

public :
	ICM_Pricer_Corridor(void) { Init();}

	/** ICM_Pricer_Corridor(ARM_Security *option, ARM_Model *mod, const ICM_Parameters& parameters , const ARM_Date&asof )
	{
		Init();
		Set(option, mod, parameters,asof);
	}
	**/ 
	~ICM_Pricer_Corridor()
	{}

	void Init()
	{
		SetName(ICM_PRICER_CORRIDOR);
		
		itsZeroCurve=NULL;
		itsDefCurveIndex=NULL;
		itsIRVolatility=NULL;
		itsIdxVolatility=NULL;
		itsCorrelation=NULL;
	}

	void Set(ARM_Security *sec, ARM_Object *mod, const ICM_Parameters& parameters , const ARM_Date&asof) ;
 
	virtual double ComputePrice(qCMPMETH measure);

	double	CorridorSubPeriodPV(double TExer,
								  double BarrierDown,
								  double BarrierUp,
								  double FwdStartDate,
								  double FwdEndDate,
								  double MaturityIdxYF,
								  double FixingCredit,
								  double& voldown,
								  double& volup,
								  double& FwdSpreadcvx);

	double	CorridorSubPeriodPVFloat(double TExer,
								  double BarrierDown,
								  double BarrierUp,
								  double VolIrCurve,
								  double FwdStartDate,
								  double FwdEndDate,
								  double MaturityIdxYF,
								  double FixingCredit,
								  double& voldown,
								  double& volup,
								  double& FwdSpreadcvx);

	void View(char* id = NULL, FILE* ficOut = NULL);

private:

	double Accrued(void) ; 
	double FeeLegPV(void) ;
	double DefLegPV(void) ;
	double Compute_Fwd_Spread(const class ARM_Date &,const class ARM_Date &, double& dur ); 
	double ComputeDuration(void) ;
	double ComputeSpread(const double &) ;
	double ComputeImpliedVol(const double &) ;
	double CptPrice_Corridor(qCMPMETH measure);
};



#endif

