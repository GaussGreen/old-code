#ifndef _ICM_PRICER_CPDO_H
#define _ICM_PRICER_CPDO_H

//#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel\pricer\icm_pricer.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_CPDO icm_pricer_cpdo.h "icm_pricer_cpdo.h"
 *  \author VR
 *	\version 1.0
 *	\date   Nov 2006
 *	\file   icm_pricer_cpdo.h
 *	\brief  Simulations for CPDO  */
/***********************************************************************************/

class ICM_Security ;

class ICM_Pricer_CPDO : public ICM_Pricer
{
private:
	//Inputs 

	// Random parameters 
	// IR
	double itsa; 
	double itsIRVol; 
	//Spreads
	double itsSpreadVol; 
	
	//Simulation
	long itsNbSimulations;		// Nb of Monte Carlo simulations
	long itsMaturityDays;		// Nb of days between AsOf and Maturity
	int itsNbPath;				// Number of paths of one simulation 
	int itsSimulationPath;		// Path of the initial Monte Carlo schedule 
								// The initial MC schedule can be modified to include roll dates and coupon dates 

	double itsRunning;			// Running spread
	double itsBidMid;			// Bid Mid costs
	double itsBidMidFloor;		// Bid Mid Floor 
	double itsBidMidCap;			// Bid Mid Cap


	double itsInitialLoss;		// Initial Loss
	double itsSpreadSlope;		// Spread Slope 
	double itsStressRoll;		// stress Roll

	bool itsIsDetailedSimulation;	// To have the details of a simulation 
	long itsNumDetailedSimulation; // Number of the detailed simulation 

	bool itsIsRandomIR;				// To generate random paths for Interest rates 
	bool itsIsRandomSpreads;		// To generate random paths for spreads 
	bool itsIsRandomDefaults;		// To generate random paths for Defaults
	bool itsIsSRSched;				

	bool itsIsCashOut;				// To know if there was a cash out (desactivation)
	bool itsIsCashIn;				// To know if there was a cash in (target reached) 

	double itsDuration;			// value of duration 
	double itsRecovery;			// value of recovery
	// End of inputs
	


	// Outputs 
	ARM_Vector itsDF0;				// Vector of Discount Factors at AsOf used for all simulations - size : NbNoRiskyCpns
	ARM_Vector itsDF;				// Vector of Discount Factors at AsOf used for all simulations - size : itsNbPath
	ARM_Vector itsFwdRate;			// Vector of forward rates used for all simulation - size : itsNbPath
	
	//ARM_Vector* itsInitSimulationSchedule;  // Vector with the days used to simulate (Nb of days from AsOf) - size : : itsNbPath
	ARM_Vector* itsSimulationSchedule;  // Vector with the days used to simulate (Nb of days from AsOf) - size : : itsNbPath
	ARM_Vector itsIsRollDate;			// Vector of Roll Dates - size : itsNbPath
	ARM_Vector itsIsRiskyCouponDate;	// Vector of Risky coupon Dates - size : itsNbPath
	ARM_Vector itsIsNoRiskyCouponDate;	// Vector of Non risky Dates - size : itsNbPath
	ARM_Vector** itsDFs;					// Vector of DF curve - size : itsNbPath
	ARM_Vector** itsDates;					// Vector of dates curve - size : itsNbPath

	ARM_Vector itsRollDate;				// Vector of Roll Dates - size : Nb of Roll Dates
	ARM_Vector itsRiskyCouponDate;		// Vector of Risky coupon Dates - size : Nb of Risky Cpn Dates
	ARM_Vector itsNoRiskyCouponDate;	// Vector of Non risky Dates - size : Nb of No Risky Cpn Dates

	ARM_Vector itsDetImpliedSpread;	// Vector of implied spreads for one simulation - size : itsNbPath
	ARM_Vector itsDetInstRate;			// Vector of instantaneous rates for one simulation - size : itsNbPath
	ARM_Vector itsDetDefaults;			// Vector of number of defaults for one simulation - size : itsNbPath 
	ARM_Vector itsDetXt;			// Vector of number of defaults for one simulation - size : itsNbPath 

	ARM_Vector itsImpliedSpread;	// Vector of implied spreads for one simulation - size : itsNbPath
	ARM_Vector itsXt;				// Vector of Xt for one simulation - size : itsNbPath
	ARM_Vector itsInstRate;			// Vector of instantaneous rates for one simulation - size : itsNbPath
	ARM_Vector itsDefaults;			// Vector of number of defaults for one simulation - size : itsNbPath 

	ARM_Vector itsRiskyValue;		// Vector of Risky value == itsMTM + itsCash - size : itsNbPath
	ARM_Vector itsTargetValue;		// Vector of target value - size : itsNbPath
	ARM_Vector itsMTM;				// Vector of MTM - size : itsNbPath
	ARM_Vector itsCash;				// Vector of cah - size : itsNbPath
	ARM_Vector itsPNL;				// Vector of PNL - size : itsNbPath
	ARM_Vector itsNotional;			// Vector of Notionals	- size : itsNbPath

	ARM_Vector itsSimulDayTarget;	// Vector of the day of the target - size : itsNbSimul
	ARM_Vector itsSimulPNL;			// Vector of the PNL at maturity - size : itsNbSimul
	ARM_Vector itsSimulIsCashIn;	// Vector of the nb of CashIn - size : itsNbSimul
	ARM_Vector itsSimulIsCashOut;	// Vector of the nb of CahsOut - size : itsNbSimul
	ARM_Vector itsValo;				// Vector of 

	long itsNbCashIn;				// Nb of Cash In
	long itsNbCashOut;				// Nb of Cash Out
	long itsNbNoInNoOut;			// Nb of No cash in no cash out

	double itsCountRebalance;		// Nb of Rebalancement
	double itsCashInPerc;			// Percentage of Cash In
	double itsCashOutPerc;			// Percentage of Cash Out
	double itsAverageCashInDay;		// Average of Cash In day
	double itsAverageCashOutDay;	// Average of Cash In day
	double itsAverageFinalVal;		// Average of Valo if No cash in, no cash out
	// End of Outputs

public :
	void Init() ;
	
	ICM_Pricer_CPDO(void) {Init();}

	ICM_Pricer_CPDO(ARM_Security* sec, ARM_Object* mod, const ICM_Parameters &params, const ARM_Date &asof)
	{
		Init();
		Set(sec, mod, params, asof);
	}

	void Set(ARM_Security *sec, ARM_Object *mod, const ICM_Parameters &params, const ARM_Date &asof) ;

	virtual ~ICM_Pricer_CPDO()
	{
		if (itsSimulationSchedule)
			delete itsSimulationSchedule;
		itsSimulationSchedule = NULL;

		if (itsDFs)
		{
			for (int i = 0; i < itsNbPath ; i++) 
			{
				if (itsDFs[i])
					delete itsDFs[i];
				itsDFs[i] = NULL; 
			}
			delete [] itsDFs;
		}

		if (itsDates)
		{
			for (int i = 0; i < itsNbPath ; i++) 
			{
				if (itsDates[i])
					delete itsDates[i];
				itsDates[i] = NULL; 
			}
			delete [] itsDates;
		}


		//if (itsInitSimulationSchedule)
		//	delete itsInitSimulationSchedule;
		//itsInitSimulationSchedule = NULL;
	}

	void View(char* id, FILE* ficOut) ;

	void AddSchedule(const ICM_Security& sec, const ARM_Date& dAsOf);

	void GenerateCpnRollDates(const ICM_Security& sec, 
							const ARM_Date& dAsOf, 
						    const ARM_Date& dPastFixingDate, 
							ARM_Vector& IsDate,
							ARM_Vector& CpnDate); 

	void Simulate();

	void SimulateOnePath(const int& NumSimu, 
						const double& Target,
						const double& UFFees,
						const double& RunningFees,
						const double& InitialValo,
						const double& ExpoOnV,
						const double& ExpoOnV0,
						const double& Alpha,
						const double& Beta,
						const double& Desactivation,
						const int& NbAssets);


	void ComputeDiscountFactorCurve(const double& NumPath, 
									 ARM_Vector& DateVect, 
									 ARM_Vector& DFVect);

	double ComputeTarget(const int& FirstPeriodLen,
						  const bool& IsCpnDate,
						  const double& V0, 
						  const double& Target, 
						  const double& RunningFees,
						  const ARM_Vector& VecDate,
						  const ARM_Vector& VecDF);

	double SetEURIB3M(const ARM_Vector& VecDate,
						const ARM_Vector& VecDF);

	int GetSizeDFCurve(const double& NumPath);

	void ComputeInstRate();
	void ComputeXt();

	void InitSimulation();
	
	void ComputeGlobalOutputs();

	void CreateSimulationSchedule();

	void PopulateVectors();

	void TestSizeNonRandomVectors();

protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ) ; 

	virtual void DoPriceVector(const std::string& measure) ;
	virtual double DoPrice(const std::string& measure); 


private:
	virtual double Accrued() ; 
	virtual double FeeLegPV () ;
	virtual double DefLegPV () ;
	virtual double ComputeDuration(void) ;
	virtual double ComputeSpread(const double& MtM = 0.) ;
 	virtual	double ComputeImpliedVol(const double& Price) ;
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate,double& dur) ;

};

#endif

