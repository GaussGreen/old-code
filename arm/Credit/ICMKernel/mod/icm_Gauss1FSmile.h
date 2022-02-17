
#error no longer part of the project

/*----------------------------------------------------------------------------*/

/* \file icm_Gauss1FSmile.h
 *
 *  \brief Calul du smile de correlation dans  un modèle de Gauss à 1 facteur
 *  \author Fakher Ben Atig
 *  \version 1.0
 *  \date April 2004
 */

/*----------------------------------------------------------------------------*/

#ifndef _GAUSS1FSMILE_H
#define _GAUSS1FSMILE_H


#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\pricer\icm_pricer_basket.h"

/*! \class Gauss1FSmile
 *  \brief Cette classe est utilisée pour le calcul du smile
 *  du Modèle de Gauss à 1 facteur
 *  Elle utilise la fonction  de minimisation de NAG c05tbc.
 */

class Gauss1FSmile : public ARM_Object
{
	
private:

	ICM_Pricer* itsPricer;
	int itsDataType ;
	double itsMktData ;
	double itsUpfrontPay;
	double itsSeed ;
	int itsSmileType ;
	double itsSmile ;
	double itsFirstPrice ;

public:
	void Init();

	Gauss1FSmile()
	{
		Init();
	}	

	Gauss1FSmile(ICM_Pricer* pricer,
				 int DataType, 
				 double MktData,
				 double UpfrontPay,
				 double Seed,
				 int SmileType);

	/*void BitwiseCopy(const ARM_Object* src)*/;
	
	void Copy(const ARM_Object* src);
	
	/*ARM_Object* Clone(void);*/

	virtual ~Gauss1FSmile(void);

public:
	void SetitsSmile (double value){itsSmile = value;}
	void SetitsFirstPrice (double value){itsFirstPrice = value;}
	
	virtual double GetitsSmile(void)
	{
		if(!itsSmile)
			Calibrate2();
		return itsSmile;
	}
	
public:

	ARM_Vector* FunctionToMinimize2(double* x);

	// Fonction qui utilise l'algorithme de NAG
	virtual void Calibrate2(void);
};

#endif
