
#ifndef _ICM_PRICER_LSS_GAP_OPTION_H
#define _ICM_PRICER_LSS_GAP_OPTION_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/crv/icm_distribloss.h"
#include "ICMKernel/inst/icm_mez.h"

/*********************************************************************************/
/*! \class  ICM_PRICER_LSS_GAP_OPTION icm_pricer_lss_gap_option.h "icm_pricer_lss_gap_option.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   July 2006 
 *	\file   icm_pricer_lss_gap_option.h
 *	\brief  Pricer for LSS gap option */
/***********************************************************************************/

class ICM_Pricer_LssGapOption : public ICM_Pricer
{
private:

	ICM_Pricer*	its_SS_Pricer;
	long		itsnbsimuls;

	double		its_drift_spreads;
	double		its_vol_spreads;

	double		its_drift_correl;
	double		its_vol_correl;

	double		its_proba_trigger;
	double		its_proba_loss;

	//vector<double> its_MTM;

public :
	ICM_Pricer_LssGapOption(void) 
	{ Init(); }

	/** ICM_Pricer_LssGapOption(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters , const ARM_Date&asof,long nbsimuls=-1)
	{
		Init();
		Set(sec, mod,parameters,nbsimuls);
	} **/ 

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters, const ARM_Date&AsOf,long nbsimuls=-1)
	{	
		ICM_Pricer::Set(sec,mod,parameters,&AsOf);

		if (nbsimuls>0) 
		{itsnbsimuls = nbsimuls;}
	}

	~ICM_Pricer_LssGapOption()
	{
		if (its_SS_Pricer)
			delete its_SS_Pricer;
		its_SS_Pricer=NULL;
	}

	void Init()
	{
		SetName(ICM_PRICER_LSS_GAP_OPTION);
		its_SS_Pricer = NULL;
		itsnbsimuls = 100;

		its_drift_spreads=0.;
		its_vol_spreads=0.4;

		its_drift_correl=0.;
		its_vol_correl=0.1;
	}

	void LoadingParameters(long& seed,long& nbdef,
						double& initial_price,double& Matu_Spread,double& loss,
						double& Simply,double& Time2Matu,	double& Leverage,double& Fast,
						double& volspreads, double& volcorrel,double& VarReduction);

	virtual double ComputePrice(qCMPMETH mode);

	//build pricer for underlying & price it
	void PriceUnderlying();

	virtual double ComputeDuration() {return CREDIT_DEFAULT_VALUE;}	
	virtual double ComputeSpread(const double& MtM = 0.) {return CREDIT_DEFAULT_VALUE;}	
 	virtual	double ComputeImpliedVol(const double& Price) {return CREDIT_DEFAULT_VALUE;}	
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) {return CREDIT_DEFAULT_VALUE;}	
private:
	virtual double Accrued() {return CREDIT_DEFAULT_VALUE;}	
	virtual double FeeLegPV () {return CREDIT_DEFAULT_VALUE;}	
	virtual double DefLegPV () {return CREDIT_DEFAULT_VALUE;}	
public:
	void FillInitalValues(ICM_Mez* cdo,
						  double matu_spread,
						  vector<double>& InitialSpreads,
						  double& Correldw0,
						  double& CorrelUp0,
						  double simplify);

	void FillInitalCorrels(ICM_Mez* mez,
							double& Correldw0,
							double& CorrelUp0);

	//generate spreads & correlations following brownian motion
	void GenerateValues(vector<double> InitialSpreads,
											double Maturity,
											double drift_spreads,
											double vol_spreads,
											vector<double>& spreads,
											double Correldw0,
											double CorrelUp0,
											double drift_correls,
											double vol_correls,
											double& correls_dw,
											double& correls_up,
											double& avgspreads);

	//generate spreads & correlations following brownian motion
	void GenerateValuesLight(vector<double>& InitialSpreads,
											vector<double>& Maturity,
											vector<double>& Triggers_spreads,
											double drift_spreads,
											double vol_spreads,
											vector<double>& spreads,
											double Correldw0,
											double CorrelUp0,
											double drift_correls,
											double vol_correls,
											double& correls_dw,
											double& correls_up,
											double& avgspreads,
											double& MatuTriggered,
											bool& istriggered,
											double VarReduction);


	ICM_Mez* GenerateCDOWithNNamesInDefault(ICM_Mez* initialCDO,
											int nbdefaults,
											double& strike_dw,
											double& strike_up);

	ICM_Mez* GenerateCDOImpactInLoss(ICM_Mez* initialCDO,
											double loss,
											double& strike_dw,
											double& strike_up);

	ICM_ModelMultiCurves* BuildFastModel(ARM_ZeroCurve* ircurve,
										 vector<double> spreads,
										 vector<double> recovery,
										 double Strike_Dw,
										 double Strike_Up,
										 double Correl_Dw,
										 double Correl_Up,
										 double T2Matu,
										 double quick = 0.);

	double GetSpreadLevel(double time2matu,double loslevel,long& row, long& col);
	inline double St(double drift, double vol,double init,double t1,double t2,double rand)
	{ 
		double value =0.;
		value = init*exp((drift-0.5*vol*vol)*(t2-t1)+vol*sqrt(t2-t1)*rand);
		return value;
	}

	double CptCorrForCDO(ICM_Mez* initialCDO,double strike);
										 

public:

	void View(char* id, FILE* ficOut)
	{
		FILE* fOut;
		char fOutName[200];

		if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		   (void) unlink(fOutName);
	       fOut = fopen(fOutName, "w"); 
		}
		else
		{
			fOut = ficOut;
		} 

		int size =0;
		fprintf(fOut, "\t\t\t ----------------- Gap Option LSS Pricer ----------------- \n\n");
		fprintf(fOut, "proba trigger activé =%2.10lf\n",its_proba_trigger);
		fprintf(fOut, "proba Loss =%2.10lf\n",its_proba_loss);
		fprintf(fOut, "\n");

		if ( ficOut == NULL )
		{
		fclose(fOut);
		}

	}
protected:
	double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double  epsilonGamma) 
	{
		return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon, epsilonGamma); 
	}


};


#endif

