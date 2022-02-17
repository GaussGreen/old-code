/*
 *
 * Copyright (c) CDC IXIS CM October 2003 Paris
 *
 * $Log: armfrmmodelmixture.cpp,v $
 * Revision 1.9  2004/02/16 13:58:32  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.8  2004/02/09 08:54:15  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.7  2004/02/04 15:18:10  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.6  2004/01/26 13:36:26  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.5  2004/01/12 07:14:58  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.4  2003/10/24 16:28:44  rguillemot
 * One factor for Mixture
 *
 * Revision 1.3  2003/10/21 07:25:14  rguillemot
 * Ajout de commetaires
 *
 *
 */

#include "armfrmmodelmixture.h"

#include "armfrmmodel.h"
#include "armfrmhwvol.h"
#include "portfolio.h"
#include "armfrmmodelmixtureminimizer.h"
#include "swaption.h"
#include "ipricer.h"


/*!
 * By default constructor
 */
ARM_FRMModelMixture::ARM_FRMModelMixture()
: ARM_Model(), itsN(0), itsFRMModels(), itsLambdas(), itsFRMPortfolio1(NULL), itsFRMPortfolio2(NULL)
{
	SetName(ARM_FRMMODELMIXTURE);
}

/*!
 * Main constructor
 */
ARM_FRMModelMixture::ARM_FRMModelMixture(
										int n,
										ARM_Date& asOfDate,
										ARM_ZeroCurve* YieldCurve,
										const std::vector<ARM_FRMVol*>& FRMVols,
										const ARM_Vector& lambdas,
										ARM_Portfolio* Pf1,
										ARM_Portfolio* Pf2,
										int nbProducts)
: ARM_Model(YieldCurve), itsN(n), itsFRMPortfolio1(Pf1), itsFRMPortfolio2(Pf2), itsNbProducts(nbProducts)										
{	
	SetName(ARM_FRMMODELMIXTURE);
	if (n > 2)
	{
        throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
                        "FRM Mixture Model allows two states maximum.");
	}

	if (FRMVols.size() != n)
	{
		throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
                        "Mismatch between the number of element in the volatility vector and the number of mixture states.");
	}

	if (lambdas.GetSize() != n-1)
	{
		throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
                        "Mismatch between the number of element in the lambda vector and the number of mixture states.");
	}

	if (((int)Pf1->GetSize()%nbProducts != 0) || (Pf1->GetSize() == 0))
	{
		throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
                        "The first portfolio don't have a valid number of produtcs.");
	}

	int i = 0;

	for (i = 0; i < itsN; ++i)
	{
		if (FRMVols[i]->GetNbFactors() >= 2)
		{
			throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
                        "Just one factor model is allowed.");
		}
	}

	for (i = 0; i < itsN; ++i)
	{
		itsFRMModels.push_back(new ARM_FRMModel(
												asOfDate,
												YieldCurve,
												FRMVols[i],
												NULL,
												NULL));
	}

	itsIndependentMode = false;
	itsIndependentVol.Resize(itsNbProducts,itsN);
	itsIndependentSpread.Resize(itsNbProducts,itsN);
	itsNbProductsPerMaturity = Pf1->GetSize()/itsNbProducts;
	itsFRMPortfolio1Error.Resize(itsNbProducts,itsNbProductsPerMaturity);
	if (itsFRMPortfolio2)
	{
		itsFRMPortfolio2Error.Resize(itsNbProducts);
	}

	itsLambdas.Resize(itsN);

	double sum = 0.0;

	for (i = 0; i < itsN-1; ++i)
	{
		itsLambdas[i] = lambdas[i];

		sum += itsLambdas[i];
	}
	
	itsLambdas[itsN-1] = 1.0-sum;
}

/*!
 * Destructor
 */
ARM_FRMModelMixture::~ARM_FRMModelMixture()
{
	deleteFRMModels();
}

/*!
 * operator= function
 */
ARM_FRMModelMixture& ARM_FRMModelMixture::operator=(const ARM_FRMModelMixture& rhs)
{
	// Test pour verifier l'auto affectation
	if (this != &rhs)
	{
		ARM_Object::operator =(rhs);

		deleteFRMModels();

		itsN = rhs.itsN;

		itsFRMModels.resize(itsN);

		int i;

		for (i = 0; i < itsN; ++i)
		{
			if (itsFRMModels[i])
			{
				itsFRMModels[i] = static_cast<ARM_FRMModel*>(rhs.itsFRMModels[i]->Clone());
			}
		}
    
		itsLambdas = rhs.itsLambdas;

		itsFRMPortfolio1 = NULL;
		if (rhs.itsFRMPortfolio1)
			itsFRMPortfolio1 = rhs.itsFRMPortfolio1;

		itsFRMPortfolio2 = NULL;
		if (rhs.itsFRMPortfolio2)
			itsFRMPortfolio2 = rhs.itsFRMPortfolio2;
	}

	return *this;
}

/*!
 * Function, which return an instance copy of the current object.
 */
ARM_Object* ARM_FRMModelMixture:: Clone(void)
{
	return new ARM_FRMModelMixture(*this);
}

/*!
 * Function which returns the number of models in the mixture
 */
int ARM_FRMModelMixture::GetN() const
{
	return itsN;
}

/*!
 * Function which returns the ith FRM Model in the mixture
 */
ARM_FRMModel* ARM_FRMModelMixture::GetModel(int i) const
{
	if ((i >= 0) || (i < itsN))
	{
		return itsFRMModels[i];
	}
	else
	{
		return NULL;
	}
}

/*!
 * Function which returns the probability vector of the model.
 */
const ARM_Vector& ARM_FRMModelMixture::GetLambdas() const
{
	return itsLambdas;
}

/*!
 * Function used by ARM_FRMModelMixtureMinimizer to calibrate the model.
 */
double ARM_FRMModelMixture::FuncToMinimizeForMixture(const ARM_Vector& x,size_t Idx)
{
	int i;

	for (i = 0; i < itsN; ++i)
	{
		if (itsIndependentMode)
		{
			itsIndependentVol.Elt(itsCurrentPos,i) = x[i];
		}
		else
		{
			itsFRMModels[i]->GetFRMVol()->SetEltToPrice(x[i]*100, Idx);
		}
	}

	double fwd = 0;

	ARM_Security* sec = itsFRMPortfolio1->GetAsset(itsNbProductsPerMaturity*Idx);

	if (sec->GetName() == ARM_CAPFLOOR)
	{
		ARM_CapFloor* capfloor = (ARM_CapFloor*) sec;
		
		fwd = capfloor->GetSwapLeg()->GetFwdRates()->Elt(0);
	}
	else if (sec->GetName() == ARM_SWAPTION)
	{
		ARM_Swaption* swaption = (ARM_Swaption*) sec;

		fwd = swaption->ComputeBSSpot(swaption->GetExpiryDate());
	}

	double sumSpread = 0.0;

	for (i = 0; i < itsN-1; ++i)
	{
		if (itsIndependentMode)
		{
			itsIndependentSpread.Elt(itsCurrentPos,i) = x[itsN+i]*100;
		}
		else
		{
			itsFRMModels[i]->GetFRMVol()->GetSpreadCurve()->GetDiscreteValues()->Elt(Idx)=x[itsN+i]*100;
		}
		sumSpread = itsLambdas[i]*(fwd+x[itsN+i]*100);
	}

	if (itsIndependentMode)
	{
		itsIndependentSpread.Elt(itsCurrentPos,itsN-1) = (fwd-sumSpread)/itsLambdas[itsN-1]-fwd;;
	}
	else
	{
		itsFRMModels[itsN-1]->GetFRMVol()->GetSpreadCurve()->GetDiscreteValues()->Elt(Idx)=(fwd-sumSpread)/itsLambdas[itsN-1]-fwd;
	}

	double MarketError = 0.0;

	for (i = 0; i < itsNbProductsPerMaturity; ++i)
	{
		double MarketValue = 0.0, ModelValue = 0.0;

		MarketValue = itsFRMPortfolio1->GetMktPrices()->Elt(itsNbProductsPerMaturity*Idx+i);
		ModelValue  = itsFRMPortfolio1->GetAsset(itsNbProductsPerMaturity*Idx+i)->ComputePrice();

		itsFRMPortfolio1Error.Elt(Idx,i) = MarketValue-ModelValue;

		MarketError += (MarketValue-ModelValue)*(MarketValue-ModelValue);
	}
	
    return MarketError;
}

/*!
 * Function used by ARM_FRMModelMixtureMinimizer to calibrate mean reversion in the model.
 */
double ARM_FRMModelMixture::FuncToMinimizeForMeanReversion(double x)
{
	int i;

	double MarketError = 0.0;

	for (i = 0; i < GetN(); ++i)
	{
		itsFRMModels[i]->GetFRMVol()->SetMeanRevParams(x);
	}

	VolCalibrate();	

	for (i = 0; i < itsFRMPortfolio2->GetSize(); ++i)
	{
		ARM_Security* sec = itsFRMPortfolio2->GetAsset(i);

			double MarketValue = 0.0, ModelValue = 0.0;

			MarketValue = itsFRMPortfolio2->GetMktPrices()->Elt(i);

			SetOneMuWithIndex(i);
			ModelValue  = itsFRMPortfolio2->GetAsset(i)->ComputePrice();

			itsFRMPortfolio2Error.Elt(i) = MarketValue-ModelValue;

			MarketError += (MarketValue-ModelValue)*(MarketValue-ModelValue);
	}
		
	return MarketError;
}

/*!
 * Fonction principale de calibration du modèle Mixture
 * Elle se compose de 2 étapes:
 * Etape 1: Calibration des options de même maturité du Portefuille 1, c'est à dire
 * de l'intégrale de la fonction sigma et des spreads.
 * Etape 2: Calcul analytique de la fonction sigma à partir des intégrales précédemment 
 * calculées. On peut aussi calller la mean reversion à partir des options du portefeuille 2.
 */
void ARM_FRMModelMixture::BootstCalibrate(
const ARM_Vector& lowerBound,
const ARM_Vector& upperBound)
{	
	// Tests préliminaires sur la tailles des bornes

	if (itsFRMPortfolio2 != NULL)
	{
		if (lowerBound.GetSize() != itsN+1)
		{
			throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
							"The lower bound size must be N (Number of mixture states) + 1 (for mean reversion).");
		}

		if (upperBound.GetSize() != itsN+1)
		{
			throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
							"The upper bound size must be N (Number of mixture states) + 1 (for mean reversion).");
		}
	}
	else
	{
		if (lowerBound.GetSize() != itsN)
		{
			throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
							"The lower bound size must be N (Number of mixture states).");
		}

		if (upperBound.GetSize() != itsN)
		{
			throw Exception(__LINE__, __FILE__,ARM_FRMMODELMIXTURE,
							"The upper bound size must be N (Number of mixture states).");
		}
	}

	size_t i, j;

	// Précalculs préliminaires

	// Le flag itsCalibrate permet de ne pas recalculer les Mu des swaptions
	// pour chaque calcul de prix.
	for (i = 0; i < itsN; ++i)
	{
		itsFRMModels[i]->SetFlagToCalibrate(K_TRUE);
	}

	// Précalcule les Mu du Portefeuille 1
	SetMu(itsFRMPortfolio1,itsNbProductsPerMaturity);

	// Etape 1

	// Pour les volatilités les bornes sont fournies par l'utilisateur.

	ARM_Vector lowerBoundWithSpread(2*itsN-1);
	ARM_Vector upperBoundWithSpread(2*itsN-1);

	for (i = 0; i < itsN; ++i)
	{
		lowerBoundWithSpread[i] = lowerBound[i];
		upperBoundWithSpread[i] = upperBound[i];
	}

	// Ce flag permet de caller les intégrales des volatilité
	// La fonction sigma sera calculer gràce à la fonction VolCalibrate
	itsIndependentMode = true;

	for (i = 0; i < itsNbProducts; i++)
	{	
		itsCurrentPos = i;

		// Pour le spread le spread les bornes sont calculées automatiquement afin que:
		// _ le spread 1 soit positif.
		// _ F + spread2 soit positif
		// NB: le spread 1 positif est associé à la volatilité 1, la plus basse, et le 
		// spread 2 est associé à la volatilité 2, la plus haute.

		for (j = 0; j < itsN-1; ++j)
		{
			ARM_Security* sec = itsFRMPortfolio1->GetAsset(itsNbProductsPerMaturity*i);

			double fwd = 0.0;

			if (sec->GetName() == ARM_CAPFLOOR)
			{
				ARM_CapFloor* capfloor = (ARM_CapFloor*) sec;
				
				fwd = capfloor->GetSwapLeg()->GetFwdRates()->Elt(0);
			}
			else if (sec->GetName() == ARM_SWAPTION)
			{
				ARM_Swaption* swaption = (ARM_Swaption*) sec;

				fwd = swaption->ComputeBSSpot(swaption->GetExpiryDate());
			}

			lowerBoundWithSpread[itsN+j] = 0.0;
			upperBoundWithSpread[itsN+j] = (fwd/itsLambdas[j]-fwd)/100;
		}

		// Minimisation de la somme des carrés de prix des produits pour chaque maturité
		// Cette minimisation est assurée par la classe ARM_FRMModelMixtureMinimizer
		
		// On appelle SetModel avant la minimisation et on appelle ensuite ComputePrice
		// afin de ne pas reaculer les flux de l'échéancier à chaque calcul de prix.
		for (j = 0; j < itsNbProductsPerMaturity; ++j)
		{
			itsFRMPortfolio1->GetAsset(itsNbProductsPerMaturity*i+j)->SetModel(this);
		}

		ARM_FRMModelMixtureMinimizer* Minimizer = new ARM_FRMModelMixtureMinimizer
			(this,
			lowerBoundWithSpread,
			upperBoundWithSpread);

		ARM_Vector x0(2*itsN-1);

		Minimizer->CalibrateMixture(x0,i);
		FuncToMinimizeForMixture(x0,i);

		delete Minimizer;
	}	

	itsIndependentMode = false;

	// Cette fonction transforme (dans le cas d'un callage de swaptions) les spreads de taux de
	// swap en spreads de taux forward.
	SpreadCalibrate();

	// Si le Portefeuille2 n'est pas vide on calibre la Mean Reversion sinon on déduit seulement
	// la fonction sigma à partir de la mean reversion entré par l'utilisateur.
	if (itsFRMPortfolio2)
	{
		CalibrateMeanReversion(lowerBound, upperBound);
	}
	else
	{
		VolCalibrate();
	}

	for (i = 0; i < itsN; ++i)
	{
		itsFRMModels[i]->SetOneMu(NULL);
	}

	for (i = 0; i < itsN; ++i)
	{
		itsFRMModels[i]->SetFlagToCalibrate(K_FALSE);
	}

	return;
}


/*!
 * Calibration de la mean reversion 
 * Cette fonction s'appuie sur les intégrales de la fonction sigma
 * calculée précédemment. Elle minimise la somme des carrés de tous 
 * les prix des instruments de Portefeuille 2.
 */
void ARM_FRMModelMixture::CalibrateMeanReversion(
const ARM_Vector& lowerBound,
const ARM_Vector& upperBound)
{
	// Précalcule les Mu du portefeuille 2
	SetMu(itsFRMPortfolio2,1);

	int i;

	for (i = 0; i < itsFRMPortfolio2->GetSize(); ++i)
	{
		itsFRMPortfolio2->GetAsset(i)->SetModel(this);
	}
	
	// Les bornes sont entrées par l'utilisateur
	ARM_Vector lowerBoundMeanRev(1);
	ARM_Vector upperBoundMeanRev(1);

	lowerBoundMeanRev.Elt(0) = lowerBound.Elt(itsN);
	upperBoundMeanRev.Elt(0) = upperBound.Elt(itsN);

	// La minimisation est assurée par la classe ARM_Minimize
	ARM_FRMModelMixtureMinimizer* Minimizer = new ARM_FRMModelMixtureMinimizer
			(this,
			lowerBoundMeanRev,
			upperBoundMeanRev);

	double MeanRev = 0.0;

	Minimizer->CalibrateMeanReversion(MeanRev);
	FuncToMinimizeForMeanReversion(MeanRev);

	delete Minimizer;
}

/*
 * Calcul analytique de la fonction sigma à partir de ses intégrales
 */
void ARM_FRMModelMixture::VolCalibrate()
{
	int i, j, k;

	for (i = 0; i < itsN; ++i)
	{
		ARM_FRMHWVol* FRMVol = (ARM_FRMHWVol*)itsFRMModels[i]->GetFRMVol();

		ARM_ReferenceValue* Curve = FRMVol->GetCurrentCurve();

		ARM_Vector* dates = Curve->GetDiscreteDates();

		ARM_Vector* sigma = new ARM_Vector(itsNbProducts);

		double meanRev = FRMVol->GetMeanRevParams();

		if (itsFRMPortfolio1->GetAsset(0)->GetName() == ARM_CAPFLOOR)
		{
			if (FRMVol->GetVolType() == K_ROW)
			{
				double sum = 0.0;

				if (meanRev != 0)
				{
					double exp2mrTj_1 = 1.0;
					for (j = 0; j < itsNbProducts; ++j)
					{
						double exp2mrTj = exp(2*meanRev*(*dates)[j]/K_YEAR_LEN);
						double volImp = itsIndependentVol.Elt(j,i);
						(*sigma)[j] = (2*meanRev*(*dates)[j]/K_YEAR_LEN*volImp*volImp-sum/exp2mrTj)/(1-exp2mrTj_1/exp2mrTj);
						if ((*sigma)[j] > 0)
						{
							(*sigma)[j] = sqrt((*sigma)[j]);
						}
						else
						{
							(*sigma)[j] = 0;
						}
						sum += (*sigma)[j]*(*sigma)[j]*(exp2mrTj-exp2mrTj_1);
						exp2mrTj_1 = exp2mrTj;
					}
				}
				else
				{
					for (j = 0; j < itsNbProducts; ++j)
					{
						double Tj_Tj_1;
						if (j == 0)
						{
							Tj_Tj_1 = (*dates)[j]/K_YEAR_LEN;
						}
						else
						{
							Tj_Tj_1 = (*dates)[j]/K_YEAR_LEN-(*dates)[j-1]/K_YEAR_LEN;
						}
						double volImp = itsIndependentVol.Elt(j,i);
						(*sigma)[j] = ((*dates)[j]/K_YEAR_LEN*volImp*volImp-sum)/(Tj_Tj_1);
						if ((*sigma)[j] > 0)
						{
							(*sigma)[j] = sqrt((*sigma)[j]);
						}
						else
						{
							(*sigma)[j] = 0;
						}
						sum += Tj_Tj_1*(*sigma)[j]*(*sigma)[j];
					}
				}
			}
			else if (FRMVol->GetVolType() == K_DIAG)
			{
				for (j = 0; j < itsNbProducts; ++j)
				{
					double volImp = itsIndependentVol.Elt(j,i);
					(*sigma)[j] = sqrt(2*meanRev*(*dates)[j]/K_YEAR_LEN/(1-exp(-2*meanRev*(*dates)[j]/K_YEAR_LEN)))*volImp;
				}
			}
		}
		else if (itsFRMPortfolio1->GetAsset(0)->GetName() == ARM_SWAPTION)
		{
			double sum = 0.0;

			ARM_Vector** mu = itsFRMModels[i]->GetMu();

			if (FRMVol->GetVolType() == K_ROW)
			{
				if (meanRev != 0)
				{
					double exp2mrTj_1 = 1.0;

					for (j = 0; j < itsNbProducts; ++j)
					{
						double coeffMu = 0.0;

						for (k = 0; k < (*mu[j]).GetSize(); ++k)
						{
							coeffMu += exp(-meanRev*(*dates)[j+k]/K_YEAR_LEN)*(*mu[j]).Elt(k);
						}

						coeffMu *= coeffMu;

						double exp2mrTj = exp(2*meanRev*(*dates)[j]/K_YEAR_LEN);
						double volImp = itsIndependentVol.Elt(j,i);
						(*sigma)[j] = (2*meanRev*(*dates)[j]/K_YEAR_LEN*volImp*volImp-sum*coeffMu)/coeffMu/(exp2mrTj-exp2mrTj_1);
						if ((*sigma)[j] > 0)
						{
							(*sigma)[j] = sqrt((*sigma)[j]);
						}
						else
						{
							(*sigma)[j] = 0;
						}
						sum += (*sigma)[j]*(*sigma)[j]*(exp2mrTj-exp2mrTj_1);
						exp2mrTj_1 = exp2mrTj;
					}
				}
				else
				{
					for (j = 0; j < itsNbProducts; ++j)
					{
						double Tj_Tj_1;
						if (j == 0)
						{
							Tj_Tj_1 = (*dates)[j]/K_YEAR_LEN;
						}
						else
						{
							Tj_Tj_1 = (*dates)[j]/K_YEAR_LEN-(*dates)[j-1]/K_YEAR_LEN;
						}
						double volImp = itsIndependentVol.Elt(j,i);
						(*sigma)[j] = ((*dates)[j]/K_YEAR_LEN*volImp*volImp-sum)/(Tj_Tj_1);
						if ((*sigma)[j] > 0)
						{
							(*sigma)[j] = sqrt((*sigma)[j]);
						}
						else
						{
							(*sigma)[j] = 0;
						}
						sum += Tj_Tj_1*(*sigma)[j]*(*sigma)[j];
					}
				}
			}
			else if (FRMVol->GetVolType() == K_DIAG)
			{
				if (meanRev != 0)
				{
					for (j = itsNbProducts-1; j >=0; --j)
					{
						double coeffMu = 0.0;

						for (k = 1; k < (*mu[j]).GetSize(); ++k)
						{
							coeffMu += exp(-meanRev*(*dates)[j+k]/K_YEAR_LEN)*(*mu[j]).Elt(k)*(*sigma)[j+k];
						}

						double a = 0.0, b = 0.0, c = 0.0;

						double expmrTj = exp(meanRev*(*dates)[j]/K_YEAR_LEN);
						double exp2mrTj = expmrTj*expmrTj;

						double volImp = itsIndependentVol.Elt(j,i);

						a = ((*mu[j]).Elt(0))*((*mu[j]).Elt(0))/exp2mrTj;

						b = 2*((*mu[j]).Elt(0))*coeffMu/expmrTj;

						c = coeffMu*coeffMu-2*meanRev*(*dates)[j]/K_YEAR_LEN*volImp*volImp/(exp2mrTj-1.0);
						
						double sqrDelta = sqrt(b*b-4*a*c);

						(*sigma)[j] = (-b+sqrDelta)/2/a;
					}
				}
				else
				{
					for (j = itsNbProducts-1; j >=0; --j)
					{
						double coeffMu = 0.0;

						for (k = 1; k < (*mu[j]).GetSize(); ++k)
						{
							coeffMu += (*mu[j]).Elt(k)*(*sigma)[j+k];
						}

						double a = 0.0, b = 0.0, c = 0.0;

						double volImp = itsIndependentVol.Elt(j,i);

						a = ((*mu[j]).Elt(0))*((*mu[j]).Elt(0));

						b = 2*((*mu[j]).Elt(0))*coeffMu;

						c = coeffMu*coeffMu-volImp*volImp;
						
						double sqrDelta = sqrt(b*b-4*a*c);

						(*sigma)[j] = (-b+sqrDelta)/2/a;
					}
				}
			}
				 
		}

		for (j = 0; j < itsNbProducts; ++j)
		{
			(*sigma)[j] *= 100;
		}

		Curve->SetDiscreteValues(sigma);

		for (j = 0; j < itsNbProducts; j++)
		{
			FRMVol->UpdateCurves(j);
		}
	}
}

/*
 * Transformation des spreads de taux de swap en spread de taux forward
 */
void ARM_FRMModelMixture::SpreadCalibrate()
{
	if (itsFRMPortfolio1->GetAsset(0)->GetName() == ARM_CAPFLOOR)
	{
		for (int i = 0; i < itsN; ++i)
		{
			ARM_FRMHWVol* FRMVol = (ARM_FRMHWVol*)itsFRMModels[i]->GetFRMVol();
			for (int j = 0; j < itsNbProducts; ++j)
			{
				FRMVol->GetSpreadCurve()->GetDiscreteValues()->Elt(j) = itsIndependentSpread.Elt(j,i);
			}
		}
	}
	else if (itsFRMPortfolio1->GetAsset(0)->GetName() == ARM_SWAPTION)
	{
		ARM_Date AsOfDate = GetStartDate();

		ARM_Swaption* Swaption0 = (ARM_Swaption *) itsFRMPortfolio1->GetAsset(0);

		ARM_Date StartDate = Swaption0->GetStartDate();
		ARM_Date EndDate   = Swaption0->GetEndDate();
        
		ARM_SwapLeg* SwaplegVar = Swaption0->GetFloatLeg();

		ARM_Vector* Fwds = SwaplegVar->GetFwdRates();
		ARM_Vector* Terms = SwaplegVar->GetInterestTerms();

		for (int i = 0; i < itsN; ++i)
		{
			double Sij = 0.0, Siplus1j = 0.0;
			double mij = 0.0, miplus1j = 0.0;
			double Aiplus1j = 0.0;
			double DFiplus1j = 1.0;

			ARM_FRMHWVol* FRMVol = (ARM_FRMHWVol*)itsFRMModels[i]->GetFRMVol();

			for (int j = itsNbProducts-1; j >=0; --j)
			{
				double mi = 0.0;

				ARM_Swaption* Swaption = (ARM_Swaption *) itsFRMPortfolio1->GetAsset(j*itsNbProductsPerMaturity);
				Sij = Swaption->PriceToRate(AsOfDate,0.0);
				mij = itsIndependentSpread.Elt(j,i);

				mi = Sij+mij-Fwds->Elt(j)+1/Terms->Elt(j)*(Sij-Siplus1j+mij-miplus1j)*Aiplus1j;
				DFiplus1j *= 1/(1+(Terms->Elt(j)*(Fwds->Elt(j)+mi))/100);
				Aiplus1j += Terms->Elt(j)*DFiplus1j;

				FRMVol->GetSpreadCurve()->GetDiscreteValues()->Elt(j) = mi;

				Siplus1j = Sij;
				miplus1j = mij;
			}
		}
	}
}

/*!
 * Function which return the price of the security using the right pricer
 */
double ARM_FRMModelMixture::Price(ARM_Security* ProductToPrice)
{
	ARM_IFPricer pricer(ProductToPrice,this);
    
    return (pricer.Price());
}

/*!
 * Function which return the swap rate of the model 
 */
double ARM_FRMModelMixture::ExpectedFwdYield(
											double startDate,
											double endDate,
											double payDate,
											int compMeth,
											int dayCount,
											int domOrFrg,
											int discYC,
											int YieldDecomp,
											double Margin)
{
	if (itsN > 0)
	{
		return itsFRMModels[0]->ExpectedFwdYield(
		startDate,
		endDate,
		payDate,
		compMeth,
		dayCount,
		domOrFrg,
		discYC,
		YieldDecomp,
		Margin);
	}

	return 0.0;
}

/*!
 * Function which return the caplet price 
 */
double ARM_FRMModelMixture::EuroCaplet(
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
										bool IsCMS,
										int FixFrequency, double UnderlyingTenor,
										int YieldDecomp, double Margin,
										ARM_Currency* ccy,
										StoreFwdRateAndCapletInfo* StoreInfo )
{
	double optRate = 0.0;

	int i = 0;

	if (itsIndependentMode)
	{
		for (i = 0; i < itsN; ++i)
		{
			double m = itsFRMModels[i]->GetFRMVol()->Mdec(resetMaturity*K_YEAR_LEN);;
			double spread = itsIndependentSpread.Elt(itsCurrentPos,i);
			double vol = itsIndependentVol.Elt(itsCurrentPos,i);

			optRate+=itsLambdas[i]*bsOption(fwd+spread+m,strike+m,vol,0.0,0.0,resetMaturity,optionType);
		}
	}
	else
	{
		for (i = 0; i < itsN; ++i)
		{
			optRate+=itsLambdas[i]*itsFRMModels[i]->EuroCaplet(
										settlement,
										resetMaturity,
										startMaturity,
										endMaturity,
										optionType,
										fwd,
										strike,
										domOrFrg,
										compMeth,
										dayCount);
		}
	}

	return optRate;
}

/*!
 * Function which return the swaption price 
 */
double ARM_FRMModelMixture::EuroSwaption(
										ARM_Security* sec,
										double startswap,
										double matswaption, 
										int optionType,
										double swapfwd,
										double Strike,
										ARM_Vector* Forwds,
										ARM_Vector* poids,
										int DomOrFrg)
{
	double optRate = 0.0;

	int i;

	if (itsIndependentMode)
	{
		for (i = 0; i < itsN; ++i)
		{
			double mij = itsFRMModels[i]->Mdec(sec);
			double spread = itsIndependentSpread.Elt(itsCurrentPos,i);
			double vol = itsIndependentVol.Elt(itsCurrentPos,i);

			optRate+=itsLambdas[i]*bsOption(swapfwd+spread+mij,Strike+mij,vol,0.0,0.0,matswaption,optionType);
		}
	}
	else
	{
		for (i = 0; i < itsN; ++i)
		{
			optRate+=itsLambdas[i]*itsFRMModels[i]->EuroSwaption(
										sec,
										startswap,
										matswaption, 
										optionType,
										swapfwd,
										Strike,
										Forwds,
										poids,
										DomOrFrg);
		}
	}
    
    return optRate;
}

/*!
 * View function, which displays all the contents of the object.
 */
void ARM_FRMModelMixture::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
    char strDate[30];
     
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

    GetStartDate().JulianToStrDate(strDate);

	fprintf(fOut, "Number of models %d:\n", itsN);

	int i,j;

	for (i = 0; i < itsN; ++i)
	{
		fprintf(fOut,"\n===> Model %d <===\n", i);
		fprintf(fOut,"Lambda %d = %f\n", i, itsLambdas[i]);

		fprintf(fOut, "Independent Vol:\n");	
		for (j = 0; j < itsIndependentVol.GetNumLines(); ++j)
		{
			fprintf(fOut, "%f\n", itsIndependentVol.Elt(j,i));	
		}

		fprintf(fOut, "Independent Spread:\n");	
		for (j = 0; j < itsIndependentSpread.GetNumLines(); ++j)
		{
			fprintf(fOut, "%f\n", itsIndependentSpread.Elt(j,i));	
		}

		ARM_FRMVol* FRMVol = itsFRMModels[i]->GetFRMVol();
		ARM_ReferenceValue* vol = FRMVol->GetCurrentCurve();
		int NumVols = vol->GetSize();
		ARM_ReferenceValue* spread = FRMVol->GetSpreadCurve();
		int NumSpreads = spread->GetSize();

		fprintf(fOut, "Mean Reversion %d: %f\n",i,FRMVol->GetMeanRevParams());	

		fprintf(fOut, "Sigma function:\n");	
		for (j = 0; j < NumVols; ++j)
        {
            double datelag = (int)vol->GetDiscreteDates()->Elt(j);
            double value   = vol->GetDiscreteValues()->Elt(j);

            fprintf(fOut,"\t%6.0lf\t\t%lf \n",datelag, value);
        }

		fprintf(fOut, "Spread function:\n");	
		for (j = 0; j < NumSpreads; ++j)
        {
            double datelag = (int)spread->GetDiscreteDates()->Elt(j);
            double value   = spread->GetDiscreteValues()->Elt(j);

            fprintf(fOut,"\t%6.0lf\t\t%lf \n",datelag, value);
        }
	}

	fprintf(fOut, "\nPortfolio 1 Error:\n");

	for (i = 0; i < itsFRMPortfolio1Error.GetNumLines(); ++i)
	{
		for (j = 0; j < itsFRMPortfolio1Error.GetNumCols(); ++j)
		{
			fprintf(fOut,"%f",itsFRMPortfolio1Error.Elt(i,j));

			if (j < itsNbProductsPerMaturity-1)
			{
				fprintf(fOut,"\t");
			}
			else
			{
				fprintf(fOut,"\n");
			}
		}
	}

	if (itsFRMPortfolio2)
	{
		fprintf(fOut, "\nPortfolio 2 Error:\n");

		for (i = 0; i < itsFRMPortfolio2Error.GetSize(); ++i)
		{
			fprintf(fOut,"%f\n",itsFRMPortfolio2Error.Elt(i));
		}
	}
         
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}

void ARM_FRMModelMixture::deleteFRMModels()
{
	int i;

	for (i = 0; i < itsN; ++i)
	{
		if (itsFRMModels[i])
		delete itsFRMModels[i];
		itsFRMModels[i] = NULL;
	}
}

/*!
 * Fonction qui calcule le vecteur des Mu pour chacun des modèles de la mixture
 * Cette fonction n'est pas un doublon de la fonction setMu du ARM_FRMModel car elle
 * prend en compte le fait d'avoir plusieurs produits de même maturité dans le 
 * portfolio.
 */
void ARM_FRMModelMixture::SetMu(ARM_Portfolio* FRMPortfolio, int NbProductsPerMaturity) 
{
	int i,j;

	for (i = 0; i < itsN; ++i)
	{
		if (itsFRMModels[i]->GetMu())
			itsFRMModels[i]->FreeVectors();

		ARM_Security* Product;
		ARM_Vector* Mu = NULL;

		ARM_Vector** MuArray = new ARM_Vector* [itsNbProducts];

		for (j = 0; j < itsNbProducts; j++)
		{
			MuArray[i] = NULL;
		}

		for (j = 0; j < itsNbProducts; j++)
		{
			Product= FRMPortfolio->GetAsset(j*NbProductsPerMaturity);
			
			// Dans le cas des caplets le vecteur des Mu vaut 1.0
			if ( Product->GetName() == ARM_CAPFLOOR )            
			{
				Mu = new ARM_Vector(1,1.0);
				MuArray[j] = Mu;
			}
        
			// Dans le cas des swaptions on calcule le vecteur des Mu à partir de
			// la relation vol SWAP / vol FRA
			if ( Product->GetName() == ARM_SWAPTION )            
			{
				Mu = itsFRMModels[i]->ComputeMu(Product);
				MuArray[j] = Mu;
			}
		}

		itsFRMModels[i]->SetGlobalMu(MuArray,itsNbProducts);
	}
}


/*
 * Cette fonction est utilisée dans pour la calibration de la Mean Reversion afin de ne pas 
 * recalculer les Mu à pour chaque calcul de prix.
 */
void ARM_FRMModelMixture::SetOneMuWithIndex(int index)
{
	int i;

	for (i = 0; i < itsN; ++i)
	{
		itsFRMModels[i]->SetOneMu(itsFRMModels[i]->GetMu()[index]);
	}
}