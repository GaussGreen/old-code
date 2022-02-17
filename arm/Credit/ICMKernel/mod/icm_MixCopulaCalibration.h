#error no longer part of the project

/*----------------------------------------------------------------------------*/

/* \file icm_MixCopulaCalibration.h
 *
 *  \brief Calibration du Modèle de Mix Copula
 *  \author Fakher Ben Atig
 *  \version 1.0
 *  \date April 2004
 */

/*----------------------------------------------------------------------------*/

#ifndef _MIXCOPULACALIBRATION_H
#define _MIXCOPULACALIBRATION_H


#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\pricer\icm_pricer_basket.h"
#include "ICMKernel\inst\icm_pf.h"

/*! \class MixCopulaCalibration
 *  \brief Cette classe est utilisée pour la calibration
 *  du Modèle de Mix Copula
 *  Elle utilise la fonction  de minimisation de NAG c05tbc.
 */

class MixCopulaCalibration : public ARM_Object
{

private:

	int itsBCType ; // type of Base Correlation (Beta Smile or Mean Correaltion)
	double itsBaseCorrelation1 ; // Base Correlation or Equivalent Beta Smile of the two tranches
	double itsBaseCorrelation2 ;
	double itsSeed1 ; // Initial guess of solutions (seeds)
	double itsSeed2 ;
	double itsImplied_Proportion_independant ;
	double itsImplied_Proportion_FullCorrel ;
	double itsSmiledPrice1 ;
	double itsSmiledPrice2 ;
	ICM_Portfolio* itsPtf;
	ICM_Pricer* itsPricer;
	
	
	double itsFirstPrice ;
	double itsImplied_UniqueBeta ;
	double itsSeedBeta;

public:
	void Init();

	MixCopulaCalibration()
	{
		Init();
	}	

	MixCopulaCalibration(ICM_Pricer* pricer, 
						ICM_Portfolio* ptf, 
						int BCType, 
						double BC1, 
						double BC2, 
						double Seed1, 
						double Seed2);

	/*void BitwiseCopy(const ARM_Object* src)*/;
	
	void Copy(const ARM_Object* src);
	
	/*ARM_Object* Clone(void);*/

	virtual ~MixCopulaCalibration(void);

public:
	void SetitsImplied_Proportion_independant(double value){itsImplied_Proportion_independant = value;}
	void SetitsImplied_Proportion_FullCorrel (double value){itsImplied_Proportion_FullCorrel = value;}

	void SetitsImplied_UniqueBeta (double value){itsImplied_UniqueBeta = value;}
	

	void SetitsSmiledPrice1 (double value){itsSmiledPrice1 = value;}
	void SetitsSmiledPrice2 (double value){itsSmiledPrice2 = value;}

	void SetitsFirstPrice (double value){itsFirstPrice = value;}
	
	virtual double Get_MixCopula_Factor(int FactorType)
	{
		if((!itsImplied_Proportion_independant) || (!itsImplied_Proportion_FullCorrel))
			Calibrate();
		if (FactorType == 0)
			return itsImplied_Proportion_independant;
		else
			return itsImplied_Proportion_FullCorrel;
	}

	virtual double CptPrice_BC (ARM_Security* security, double BC, int BCType);
	virtual double CptPrice_MixCopula (ARM_Security* security, double proportion_independant,double proportion_full_correl);
	virtual double CptPrice_SameBeta_Indep(ARM_Security* Sec, double proportion_independant,double Beta);

	virtual double Get_ReducedMixCopula_Factor(int FactorType)
	{
		if((!itsImplied_Proportion_independant) || (!itsImplied_UniqueBeta))
			Calibrate3();
		if (FactorType == 0)
			return itsImplied_Proportion_independant;
		else
			return itsImplied_UniqueBeta;
	}

public:

	ARM_Vector* FunctionToMinimize  (double* x);
	ARM_Vector* FunctionToMinimize3 (double* x);

	// Fonction qui utilise l'algorithme de NAG
	virtual void Calibrate (void);
	virtual void Calibrate3(void);
	
};

#endif
