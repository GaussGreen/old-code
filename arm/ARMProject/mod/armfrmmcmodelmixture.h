/*
 *
 * Copyright (c) CDC IXIS CM October 2003 Paris
 *
 * $Log: armfrmmcmodelmixture.h,v $
 * Revision 1.3  2004/01/26 13:36:56  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.2  2003/10/21 07:26:41  rguillemot
 * Ajout de commentaires
 *
 *
 */

/*----------------------------------------------------------------------------*/

/*! \file armfrmmcmodelmixture.h
 *
 *  \brief Monte Carlo Method for Mixture Model
 *
 *	\author  Richard GUILLEMOT
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef ARMFRMMCMODELMIXTURE_H
#define ARMFRMMCMODELMIXTURE_H

#include "mcmodel.h"

#include <vector>

#include <mcmodel.h>

class ARM_FRMModelMixture;
class ARM_FRMMCModel;

/*! \class   ARM_FRMMCModelMixture
 *  \brief This class aggregates frew ARM_FRMMCModel. Its main interest
 *  is to resuse random variables to simulate forward in each  model.
 *
 *  \author Richard GUILLEMOT
 *  \version 1.0
 */
class ARM_FRMMCModelMixture : public ARM_MCModel
{
	private:
		std::vector<ARM_FRMMCModel*> itsFRMMCModels;
		ARM_FRMModelMixture* itsFRMModelMixture;
	public : 
		// Constructors
		ARM_FRMMCModelMixture(void);

	    ARM_FRMMCModelMixture(
			ARM_FRMModelMixture* FRMModelMixture,
			long nbTraj,
			long nbStepIn,
			ARM_PRICER_TYPE PricerType);

		ARM_FRMMCModelMixture(const ARM_FRMMCModelMixture& ARM_FRMMCModelMixture);

	   ~ARM_FRMMCModelMixture(void);

		ARM_FRMMCModelMixture& operator=(const ARM_FRMMCModelMixture& rhs);

        void FreeMemory(void);

		// Services
        virtual ARM_Object* Clone(void);

		// MC Model functions
        void BeFittedTo(ARM_Security* sec);
        
        void PropagateFwds(double * bmInnov);

        ARM_PathGenerator* Generator(void);

        void CptOnePath(double * bmInnov);

        ARM_PRICER_TYPE PricerType(void);

		double PayoffFN(double * bmInnov, ARM_Security* sec);
};

#endif
