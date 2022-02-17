#pragma warning (disable : 4786 )

#if !defined(_ICM_CPDO_H_)
#define _ICM_CPDO_H_

#include "ICMKernel\inst\icm_leg.h"

/*********************************************************************************/
/*! \class  ICM_cpdo icm_cpdo.h "icm_cpdo.h"
 *  \author VR  
 *	\version 1.0
 *	\date   Nov 2006
 *	\file   icm_cpdo.h
 *	\brief  Simulation CPDO 
/***********************************************************************************/

//class ICM_Credit_Index;

class ICM_cpdo : public ICM_Security  
{
private: 
	ICM_Leg*			itsRiskyLeg;			// Risky Leg  	
	ICM_Leg*			itsRollLeg;				// Roll Leg (to know the roll dates) 
	ICM_Leg*			itsNoRiskyLeg;			// Non risky Leg	

	ARM_Date			itsCPDOMaturity;		// Maturity of the CPDO
	double				itsInitialValo;			// Initial Valo to invest in CPDO
	double				itsTarget;				// Target coupon (ex; 200bps)
	std::string			itsCpnType;				// Type of coupon for non risky asset (EUR3M - INFL ...)
	
	double				itsUFFees;				// Up Front fees
	double				itsRunningFees;			// Running fees
	
	double				itsVExpo;				// Exposition on V
	double				itsV0Expo;				// Exposition on V0
	
	double				itsAlpha;				// Rebalancing parameters
	double				itsBeta;				// Rebalancing parameters
	
	double				itsDesactivation;		// Value of desactivation 
	
	int					itsNbAssets;			// Nb of Assets 

public:

	void Init();

	ICM_cpdo() {Init();}	
	ICM_cpdo(const ICM_cpdo& cpdo) : ICM_Security(cpdo)
	{ 
		Init();
	    BitwiseCopy(&cpdo);  
	}	

	ICM_cpdo(ICM_Leg* RiskyLeg, 
			ICM_Leg* RollLeg,
			ICM_Leg* NoRiskyLeg,
			const double& InitialValo,
			const double& Target,
			const ARM_Date& Maturity,
			const std::string& CpnType,
			const double& UFFees,
			const double& RunningFees,
			const double& VExpo,
			const double& V0Expo,
			const double& Alpha,
			const double& Beta,
			const double& Desactivation, 
			const int& NbAssets);


	void Set(ICM_Leg* RiskyLeg, 
			ICM_Leg* RollLeg,
			ICM_Leg* NoRiskyLeg,
			const double& InitialValo,
			const double& Target,
			const ARM_Date& Maturity,
			const std::string& CpnType,
			const double& UFFees,
			const double& RunningFees,
			const double& VExpo,
			const double& V0Expo,
			const double& Alpha,
			const double& Beta,
			const double& Desactivation,
			const int& NbAssets);

	virtual ~ICM_cpdo()
	{
		if (itsRiskyLeg)
			delete itsRiskyLeg;
		itsRiskyLeg = NULL;

		if (itsRollLeg)
			delete itsRollLeg;
		itsRollLeg = NULL;

		if (itsNoRiskyLeg)
			delete itsNoRiskyLeg;
		itsNoRiskyLeg = NULL;

	}	

	// -----------------------------------
	// Copy and Clone Methode
	// -----------------------------------
	void BitwiseCopy(const ARM_Object* srcCPDO);

	void Copy(const ARM_Object* srcCPDO);

	ARM_Object* Clone(void);

	// -----------------------------------
	// View Methode
	// -----------------------------------
	void View(char* id, FILE* ficOut);

	// -----------------------------------
	// Get Methodes
	// -----------------------------------
	double GetInitialValo()	const 			{ return itsInitialValo;}
	double GetTarget() const				{ return itsTarget;}
	ARM_Date GetCPDOMaturity() const		{ return itsCPDOMaturity;}
	std::string GetCpnType() const			{ return itsCpnType;}
	double GetUFFees() const				{ return itsUFFees;}
	double GetRunningFees()	const			{ return itsRunningFees;}
	double GetVExpo() const					{ return itsVExpo;}
	double GetV0Expo() const				{ return itsV0Expo;}
	double GetAlpha() const					{ return itsAlpha;}
	double GetBeta() const					{ return itsBeta;}
	double GetDesactivation() const			{ return itsDesactivation;}
	int GetNbAssets() const					{ return itsNbAssets;}

	const ICM_Leg& GetRiskyLeg() const 
	{
		if (!itsRiskyLeg) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"itsRiskyLeg: risky leg does not exists"); 
		return *itsRiskyLeg; 
	}

	const ICM_Leg& GetRollLeg() const 
	{
		if (!itsRollLeg) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"itsRollLeg: roll leg does not exists"); 
		return *itsRollLeg; 
	}

	const ICM_Leg& GetNoRiskyLeg() const 
	{
		if (!itsNoRiskyLeg) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"itsNoRiskyLeg: non risky leg does not exists"); 
		return *itsNoRiskyLeg; 
	}


};

#endif 