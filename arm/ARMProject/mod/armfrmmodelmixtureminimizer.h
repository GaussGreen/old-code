/*
 *
 * Copyright (c) CDC IXIS CM October 2003 Paris
 *
 * $Log: armfrmmodelmixtureminimizer.h,v $
 * Revision 1.4  2004/01/26 13:38:35  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.3  2003/10/21 07:26:11  rguillemot
 * Ajout de commentaires
 *
 *
 */

/*----------------------------------------------------------------------------*/

/*! \file armfrmmodelmixtureminimizer.h
 *
 *  \brief Modèle mixture Calibration
 *  \author Richard GUILLEMOT
 *  \version 1.0
 *  \date October 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _ARMMINIMIZE_H
#define _ARMMINIMIZE_H

#include "calibration.h"

class  ARM_FRMModelMixture;
struct ModelAndIndex;

/*! \class ARM_FRMModelMixture
 *  \brief Cette classe est utilisée par la classe ARM_FRMModelMixture
 * pour le callage des volatilités et du spread, ainsi que de la 
 * mean reversion. Elle s'appuie sur un algorithme de minimisation
 * NAg e04jbc. 
 * Les deux fonctions principales (les plus commentées) sont:
 * _ CalibrateMixture : qui calibre le modèle mixture pour les 
 * options d'une même maturité (cf. Index)
 * _ CalibrateMeanReversion : qui calibre la mean reversion 
 * à partir d'un portefeuille d'options
 */
class ARM_FRMModelMixtureMinimizer : public ARM_Calibration
{
public:
	// Constructors and destructor
	ARM_FRMModelMixtureMinimizer();

	ARM_FRMModelMixtureMinimizer(ARM_FRMModelMixture* model,
		const ARM_Vector& lowerBound,
		const ARM_Vector& upperBound);

	ARM_FRMModelMixtureMinimizer(const ARM_FRMModelMixtureMinimizer& rhs);

	virtual ~ARM_FRMModelMixtureMinimizer(void);

	ARM_FRMModelMixtureMinimizer& operator=(const ARM_FRMModelMixtureMinimizer& rhs);

	// Services
	virtual ARM_Object* Clone(void);
		
public:
	// Calibration methods
	void CalibrateMixture(ARM_Vector& x, int index);
	void CalibrateMeanReversion(double& x);

private:
	// Private function needed for calibration
	void StochasticMinimisation(ARM_Vector& x, int index);
	void sortVolParam(ARM_Vector& x);
	bool testValidSolution(const ARM_Vector& x);
	// Fonction qui encapsule l'appelle à NAg
	void Minimize(ModelAndIndex* modelAndIndex, int n, ARM_Vector& x);

private :
	ARM_Vector itsLowerBound;
	ARM_Vector itsUpperBound;
	ARM_FRMModelMixture* itsFRMModelMixture;
};

#endif
