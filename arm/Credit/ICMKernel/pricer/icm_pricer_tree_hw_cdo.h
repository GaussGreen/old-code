
#ifndef _ICM_Pricer_tree_hw_cdo_H
#define _ICM_Pricer_tree_hw_cdo_H

#include "ICMKernel/pricer/icm_pricer_homogeneous.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel/mod/icm_hw_credit_tree_dl.h"
#include "ICMKernel\random\icm_randomnag.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_tree_hw_cdo ICM_Pricer_tree_hw_cdo.h "ICM_Pricer_tree_hw_cdo.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   February 2007 
 *	\file   ICM_Pricer_tree_hw_cdo.h
 *	\brief  Hull & White model for term structure loss */
/***********************************************************************************/

class ICM_Pricer_tree_hw_cdo : public ICM_Pricer_Distrib
{

static  ICM_QMatrix<double> itsCombinations; //Pascal Combinations for full homogeneous case

private:

	ICM_HWTree_intensity* itsTree;
	double		itsTimeStep;
	ARM_Vector	itsSchedule;
	ARM_Vector	itsResetSchedule;
	ARM_Vector	itsMergedSchedule;

	ARM_Vector	itsParamsSchedule;
	ARM_Vector	itsIntensity;
	ARM_Vector	itsFlatIntensity;
	ARM_Vector	itsJumpSize;
	ARM_Vector	itsAlpha;
	ARM_Vector	itsBeta;
	ARM_Vector	itsFuncType;
	
	double		itsResetSubPercent;	
	double		itsResetSub;	
	double		itsMaturitySub;	

	//calibration only
	int			its_call_noperiod;
	double		its_call_proba;
	double		its_call_lastTheta;

	//Monte-Carlo
	double			its_mc_pricing_activate;
	int				its_mc_nbsimuls;
	ICM_RandomNag	its_mc_RandGen;

	//Pricing Only
	double			itsPtfLosse;
	double			itsPercentTrigger;
	double			itsPercentLoss;
	double			itsPercentLosstrigger;
	double			itsStrike;
	double			itsExerciceDate;

	//Matrix
	ICM_QMatrix<double>*  itsFwdMatrix;

public :
	ICM_Pricer_tree_hw_cdo(void) 
	{ Init(); }

	ICM_Pricer_tree_hw_cdo(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters, const ARM_Date&AsOf,long nbsimuls=-1) 
	{
		Init();
		Set(sec,mod,parameters,AsOf,nbsimuls);
	}

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters, const ARM_Date&AsOf,long nbsimuls=-1)
	{	
		ICM_Pricer::Set(sec,mod,parameters,&AsOf);

		ARM_Vector* INTENSITY = parameters.GetColVect("INTENSITY");
		if (INTENSITY) itsIntensity = *INTENSITY;

		ARM_Vector* SCHEDULE = parameters.GetColVect("SCHEDULE");
		if (SCHEDULE) itsParamsSchedule = *SCHEDULE;

		ARM_Vector* JUMP = parameters.GetColVect("JUMP");
		if (JUMP) itsJumpSize = *JUMP;

		ARM_Vector* ALPHA = parameters.GetColVect("ALPHA");
		if (ALPHA) itsAlpha = *ALPHA;

		ARM_Vector* BETA = parameters.GetColVect("BETA");
		if (BETA) itsBeta = *BETA;

		ARM_Vector* TimeStep = parameters.GetColVect("TIME_STEP");
		if (TimeStep) itsTimeStep = TimeStep->Elt(0);
		
		ARM_Vector* mc_pricing = parameters.GetColVect("MC_PRICING");
		if (mc_pricing) its_mc_pricing_activate = mc_pricing->Elt(0);

		ARM_Vector* mc_nbsimuls = parameters.GetColVect("MC_NBSIMULS");
		if (mc_nbsimuls) its_mc_nbsimuls = (int) mc_nbsimuls->Elt(0);

		ARM_Vector* resetSub = parameters.GetColVect("RESETSUB");
		if (resetSub) itsResetSubPercent = (double) resetSub->Elt(0);

		ARM_Vector* MaturitySub = parameters.GetColVect("MATURITYSUB");
		if (MaturitySub) itsMaturitySub = (double) MaturitySub->Elt(0);

		ARM_Vector* Strike = parameters.GetColVect("STRIKE");
		if (Strike) itsStrike = (double) Strike->Elt(0);

		ARM_Vector* ExerciceDate = parameters.GetColVect("EXDATE");
		if (ExerciceDate) itsExerciceDate = (double) ExerciceDate->Elt(0);

		ARM_Vector* FuncType = parameters.GetColVect("FUNCTYPE");
		if (FuncType) itsFuncType = *FuncType;

		ARM_Vector* FlatIntensity = parameters.GetColVect("FLAT_INTENSITY");
		if (FlatIntensity) itsFlatIntensity = *FlatIntensity;

	}

	virtual void RefreshParameters(void) 
	{
		ICM_Parameters parameters = GetParameters();

		ARM_Vector* INTENSITY = parameters.GetColVect("INTENSITY");
		if (INTENSITY) itsIntensity = *INTENSITY;

		ARM_Vector* SCHEDULE = parameters.GetColVect("SCHEDULE");
		if (SCHEDULE) itsParamsSchedule = *SCHEDULE;

		ARM_Vector* JUMP = parameters.GetColVect("JUMP");
		if (JUMP) itsJumpSize = *JUMP;

		ARM_Vector* ALPHA = parameters.GetColVect("ALPHA");
		if (ALPHA) itsAlpha = *ALPHA;

		ARM_Vector* BETA = parameters.GetColVect("BETA");
		if (BETA) itsBeta = *BETA;

		ARM_Vector* TimeStep = parameters.GetColVect("TIME_STEP");
		if (TimeStep) itsTimeStep = TimeStep->Elt(0);

		ARM_Vector* FlatIntensity = parameters.GetColVect("FLAT_INTENSITY");
		if (FlatIntensity) itsFlatIntensity = *FlatIntensity;

	}

	~ICM_Pricer_tree_hw_cdo()
	{
		if (itsTree)
			delete itsTree;
		itsTree = NULL;

		if (itsFwdMatrix)
			delete itsFwdMatrix;
		itsFwdMatrix = NULL;
	}

	void Init()
	{
		SetName(ICM_PRICER_TREE_HW_CDO);
		itsTimeStep=0.;
		itsTree = NULL;
		its_call_lastTheta=0.;
		its_mc_nbsimuls=10;
		its_mc_pricing_activate=0.;
		itsResetSub=0.;
		itsResetSubPercent=0.;	
		itsPtfLosse=0.;
		itsMaturitySub=1000.;
		itsPercentTrigger=0.;
		itsPercentLoss=0.;
		itsPercentLosstrigger=0.;
		itsFwdMatrix=NULL;
		itsStrike=0.;
		itsExerciceDate=0.;
	}

	virtual void CptExpectedLossTranche();
	ARM_Vector& GenCDOSchedule(bool cdo=true);

	void BuildTree();
	//virtual double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
	void CalibrateTheta(int& period);
	void HomogEL(int& period,int nbnames,const double& LossUnit,const double& tranche_down,const double& tranche_up);
	double Evaluate(const double& x);

	void GeneratePaths(vector<double>&	pdefaults,
					   vector<double>&	spreads,
					   vector<int>&		namesindefault,
					   vector<double>&  cumulativeloss);

	void CDOLoss(vector<double>&  cumulativeloss,
				vector<double>&  cdoloss,
				vector<double>&  spreads,
				vector<double>&  pdefaults,
				bool& trigger,
				int& triggerperiod);

	void OCLoss(vector<double>&  cumulativeloss,
				vector<double>&  cdoloss,
				vector<double>&  spreads,
				vector<double>&  pdefaults,
				bool& trigger,
				int& triggerBegin,
				int& triggerend);

	void PayOffCDO(const vector<double>&  cumulativeloss,qCMPMETH mode,double& feepev,double& defpev,int endperiod);
	void PayOffCDOFwd(const vector<double>&  cumulativeloss,qCMPMETH mode,double& feepev,double& defpev,int start,int end);
	double MCPricing_OC(qCMPMETH mode);
	double MCPricing_CDO(qCMPMETH mode);
	double CallOptionPricing(qCMPMETH mode);

	virtual double ComputePrice(qCMPMETH mode) ;
	void View(char* id, FILE* ficOut);

	double Pdef2Spread(double& date1,double& date2,vector<double>&  pdef,bool spreadok = true);
	void Pdef2Spreads(vector<double>&  pdef,vector<double>&  spreads);

	ICM_QMatrix<double>* GetFwdMatrix() {return itsFwdMatrix;}
	void ResetFwdMatrix()
	{
		if (itsFwdMatrix)
			delete itsFwdMatrix;
		itsFwdMatrix=NULL;
	}
		
	void ComputeFwdMatrix();
};


#endif

