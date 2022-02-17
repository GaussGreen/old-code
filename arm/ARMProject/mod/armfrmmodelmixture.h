/*
 *
 * Copyright (c) CDC IXIS CM October 2003 Paris
 *
 * $Log: armfrmmodelmixture.h,v $
 * Revision 1.8  2004/02/16 13:58:40  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.7  2004/02/09 08:54:19  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.6  2004/02/04 15:17:46  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.5  2004/01/26 13:36:05  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.4  2004/01/12 07:14:24  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.3  2003/10/21 07:24:38  rguillemot
 * Ajout de commentaires
 *
 *
 */

#ifndef _ARMFRMMODELMIXTURE_H
#define _ARMFRMMODELMIXTURE_H

#include "model.h"

#include <vector>

class ARM_FRMVol;
class ARM_FRMModel;

/*----------------------------------------------------------------------------*/

/*! \file armfrmmodelmixture.h
 *
 *  \brief Modèle mixture
 *  \author Richard GUILLEMOT
 *  \version 1.0
 *  \date October 2003
 */

/*----------------------------------------------------------------------------*/

/*! \class ARM_FRMModelMixture
 *  \brief La classe modèle ARM_FRMModelMixture est un agrégat de plusieurs
 * ARM_FRMModel. Dans ce modèle les prix sont une somme pondérées de 
 * prix dans un modèle FRM.
 * La principale fonctionalité de cette classe est la calibration réalisée par
 * la méthode BootsCalibrate(cf. commentaires dans .cpp).
 *
 *  \author Richard GUILLEMOT
 *  \version 1.0
 */
class ARM_FRMModelMixture : public ARM_Model  
{
public:
	ARM_FRMModelMixture::ARM_FRMModelMixture();

	ARM_FRMModelMixture(
						int n,
						ARM_Date& asOfDate,
						ARM_ZeroCurve* zc,
						const std::vector<ARM_FRMVol*>& FRMVols,
						const ARM_Vector& lambdas,
						ARM_Portfolio* Pf1,
						ARM_Portfolio* Pf2,
						int nbProducts);
	
	~ARM_FRMModelMixture(void);

	// Affectation operator
	ARM_FRMModelMixture& operator=(const ARM_FRMModelMixture& rhs);

	// Services
	virtual ARM_Object* Clone(void);
	virtual void View(char* id, FILE* ficOut);
	
public:
	//Accesseurs
	int GetN() const;
	ARM_FRMModel* GetModel(int i) const;
	const ARM_Vector& GetLambdas() const;

	// Méthodes de calibration
	void BootstCalibrate(
						const ARM_Vector& lowerBound,
						const ARM_Vector& upperBound);
	
	void ARM_FRMModelMixture::CalibrateMeanReversion(
		const ARM_Vector& lowerBound,
		const ARM_Vector& upperBound);

	void VolCalibrate();
	void SpreadCalibrate();


	// Fonctions appélées par la classe ARM_FRMModelMixtureMinimizer
	double FuncToMinimizeForMixture(const ARM_Vector& x,size_t Idx);
	double FuncToMinimizeForMeanReversion(double x);

	// Fonctions de pricing
	double Price(ARM_Security* ProductToPrice);

	double ExpectedFwdYield(
							double startDate,
							double endDate,
							double payDate,
							int compMeth,
							int dayCount,
							int domOrFrg,
							int discYC,
							int YieldDecomp = K_COMP_PROP,
							double Margin = 0.0);

	double EuroCaplet(
						double settlement,
						double resetMaturity,
						double startMaturity,
						double endMaturity,
						double payMaturity,
						int optionType,
						double fwd,
						double strike,
						int domOrFrg,
						int compMeth,
						int dayCount,
						bool IsCMS = false,
						int FixFrequency = 1,
						double UnderlyingTenor = 0.0,
						int YieldDecomp = K_COMP_PROP,
						double Margin = 0.0,
						ARM_Currency* ccy = NULL,
						StoreFwdRateAndCapletInfo* StoreInfo = NULL);


	double EuroSwaption(
						ARM_Security* sec,
						double startswap,
						double matswaption, 
						int optionType,
						double swapfwd,
						double Strike,
						ARM_Vector* Forwds,
						ARM_Vector* poids,
						int DomOrFrg=1);

private:
	void deleteFRMModels();
	void SetMu(ARM_Portfolio* FRMPortfolio, int NbProductsPerMaturity);
	void SetOneMuWithIndex(int index);

private:
	int itsN;
	int itsNbProducts;
	int itsNbProductsPerMaturity;
	std::vector<ARM_FRMModel*> itsFRMModels;
	ARM_Vector itsLambdas;
	bool itsIndependentMode;
	ARM_Matrix itsIndependentVol;
	ARM_Matrix itsIndependentSpread;
	int itsCurrentPos;
	ARM_Portfolio* itsFRMPortfolio1;
	ARM_Matrix itsFRMPortfolio1Error;
	ARM_Portfolio* itsFRMPortfolio2;
	ARM_Vector itsFRMPortfolio2Error;
};

#endif
