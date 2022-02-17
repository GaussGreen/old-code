/*
 *
 * Copyright (c) CDC IXIS CM October 2003 Paris
 *
 * $Log: armfrmmcmodelmixture.cpp,v $
 * Revision 1.3  2004/01/26 13:37:50  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.2  2003/10/21 07:27:37  rguillemot
 * Ajout de commentaires
 *
 *
 */


#include "armfrmmcmodelmixture.h"

#include "armfrmmodel.h"
#include "armfrmmcmodel.h"
#include "armfrmmodelmixture.h"


/*!
 * By default constructor
 */
ARM_FRMMCModelMixture:: ARM_FRMMCModelMixture(void)
: ARM_MCModel(), itsFRMMCModels()
{
	SetName(ARM_FRMMCMODELMIXTURE);
}

/*!
 * Main Constructor
 */
ARM_FRMMCModelMixture::ARM_FRMMCModelMixture(
ARM_FRMModelMixture* FRMModelMixture,
long nbTraj,
long nbStepIn,
ARM_PRICER_TYPE PricerType)
: ARM_MCModel()
{
	SetZeroCurve(FRMModelMixture->GetZeroCurve());
    SetStartDate(GetZeroCurve()->GetAsOfDate());
    SetNbIters(nbTraj);
	Init();

	int nbModels = FRMModelMixture->GetN();

	int i;

	for (i = 0; i < nbModels; ++i)
	{
		ARM_FRMMCModel* FRMMCModel = new ARM_FRMMCModel(
			FRMModelMixture->GetModel(i),
			nbTraj,
			nbStepIn,
			PricerType);

		itsFRMMCModels.push_back(FRMMCModel);
	}

	itsFRMModelMixture = FRMModelMixture;
}

/*!
 * Copy Constructor
 */
ARM_FRMMCModelMixture::ARM_FRMMCModelMixture(const ARM_FRMMCModelMixture& rhs)
: ARM_MCModel(rhs)
{
	SetName(ARM_FRMMCMODELMIXTURE);
	
	if (this != &rhs)
	{
		ARM_MCModel::operator =(*this);
		itsFRMModelMixture = rhs.itsFRMModelMixture;
	}
}

/*!
 * Destructor
 */
ARM_FRMMCModelMixture::~ARM_FRMMCModelMixture()
{
	FreeMemory();
}

/*!
 * Fonction to free the memory allocated by the class. It is used by
 * the destructor.
 */
void ARM_FRMMCModelMixture::FreeMemory(void)
{
	int i;

	for (i = 0; i < itsFRMMCModels.size(); ++i)
	{
		if (itsFRMMCModels[i])
		{
			delete itsFRMMCModels[i];
		}
	}
}

/*!
 * Affectation operator
 */
ARM_FRMMCModelMixture& ARM_FRMMCModelMixture::operator=(const ARM_FRMMCModelMixture& rhs)
{
	if (this != &rhs)
	{
		ARM_Model::operator =(rhs);

		int i;

		for (i = 0; i < itsFRMMCModels.size(); ++i)
		{
			ARM_Model* Model_i = static_cast<ARM_Model*>(rhs.itsFRMMCModels[i]->Clone());
			ARM_FRMMCModel* FRMMCModel_i = dynamic_cast<ARM_FRMMCModel *>(Model_i);
			itsFRMMCModels.push_back(FRMMCModel_i);
		}
	}
	return *this;
}

/*!
 * Function, which return an instance copy of the current object.
 */
ARM_Object* ARM_FRMMCModelMixture::Clone(void)
{
	return new ARM_FRMMCModelMixture(*this);
}

/*!
 * Function to prepare model before the pricing of specific
 * security.
 */
void ARM_FRMMCModelMixture::BeFittedTo(ARM_Security* sec)
{
	int i;

	for (i = 0; i < itsFRMMCModels.size(); ++i)
	{
		itsFRMMCModels[i]->BeFittedTo(sec);
	}
}

/*!
 * Compute random variables needed to simulate forward.
 */
ARM_PathGenerator* ARM_FRMMCModelMixture::Generator(void)
{
	if (itsFRMMCModels.size() >= 1)
	{
		return itsFRMMCModels[0]->Generator();
	}
	else
	{
		return NULL;
	}
}

/*!
 * Compute forward vector values with the vector of random variables.
 */
void ARM_FRMMCModelMixture::CptOnePath(double*bmInnov)
{
	PropagateFwds(bmInnov);
}

/*!
 * Compute forward vector values with the vector of random variables.
 */
void ARM_FRMMCModelMixture::PropagateFwds(double * bmInnov)
{
	int i;

	for (i = 0; i < itsFRMMCModels.size(); ++i)
	{
		itsFRMMCModels[i]->PropagateFwds(bmInnov);
	}
}

/*!
 * Compute the payoff of the security for a specific trajectory for each model
 * of the Mixture and aggregate them with the probability of each model.
 */
double ARM_FRMMCModelMixture::PayoffFN(double * bmInnov, ARM_Security* sec)
{
	double payoff = 0.0;

	int i;

	for (i = 0; i < itsFRMMCModels.size(); ++i)
	{
		double lambda = itsFRMModelMixture->GetLambdas().Elt(i);
		payoff += lambda*itsFRMMCModels[i]->PayoffFN(bmInnov,sec);
	}

	return payoff;
}

/*!
 * Function which returns the pricer type of the model.
 */
ARM_PRICER_TYPE ARM_FRMMCModelMixture::PricerType(void)
{
	if (itsFRMMCModels.size() >= 1)
	{
		return itsFRMMCModels[0]->PricerType();
	}
	else
	{
		return PT_NONE;
	}
}
