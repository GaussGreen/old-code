
#include "ICMKernel\util\icm_CorrelFunction.h"
#include "ICMKernel\util\icm_gaussian.h"
#include "ICMKernel/util/icm_Rootfinder1D.h"
#include "ICMKernel\glob\icm_maths.h"
#include <math.h>
#include <algorithm>


// Constructor

ICM_CorrelFunction::ICM_CorrelFunction(int FactorSize, qCopula_TYPE copule, qCAL_INDEX_CORR_TYPE type, vector<double> &param,double Volatility)
{
	Init();
	Set(FactorSize, copule, type, param, Volatility);
	InitVAInter();
}

void ICM_CorrelFunction::Set(int FactorSize, qCopula_TYPE copule, qCAL_INDEX_CORR_TYPE type, vector<double> &param,double Volatility)
{
	SetFactorSize(FactorSize);
	SetCorrelFunctionType(type);
	SetCopulaType(copule);
	
	//Size
	int i=0;
	int size_param=param.size();
	itsVolatility = Volatility;


	//Init Param : 
		// PWC : n levels, n-1 thresholds (sorted)
		// PWL : n (lin ceof, level), n-1 threshold (sorted)

	//itsNbLevels = 2;

	switch (itsCorrelFunctionType)
	{
	case qCAL_PWC_CORREL:
		
		//Size

		itsNbLevels = (size_param+1)/2;
/*		if ((size_param+1)%2 != 0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"ERROR : PWC Parameters Size not OK (n levels + n-1 threshold)");
		}
			
		else
*/
		{
			itsLinCoef.resize(itsNbLevels);
			itsLevels.resize(itsNbLevels);
			itsThreshold.resize(itsNbLevels-1);
			
			//Init level
			for (i=0; i<itsNbLevels-1; i++)
			{
				itsLinCoef[i]=0.;
				itsLevels[i]=param[i];
				itsThreshold[i]=param[itsNbLevels+i];
			}
			itsLinCoef[itsNbLevels-1]=0.;
			itsLevels[itsNbLevels-1]=param[itsNbLevels-1];
			//Tri du vecteur des seuils
			std::sort(itsThreshold.begin(), itsThreshold.end());
		}
		break;
	case qCAL_PWL_CORREL:
		//Size
	
		itsNbLevels = (size_param+1)/2;
/*		if ((size_param+1)%3 != 0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"ERROR : PWL Parameters Size not OK (n (coef,levels) + n-1 threshold)");
		}
		else
*/
		{
			itsLinCoef.resize(itsNbLevels);
			itsLevels.resize(itsNbLevels);
			itsThreshold.resize(itsNbLevels-1);
			
			//Init level
			for (i=0; i<itsNbLevels-1; i++)
			{
				itsLinCoef[i]=param[i];
				itsLevels[i]=param[itsNbLevels+i];
				itsThreshold[i]=param[2*itsNbLevels+i];
			}
			itsLinCoef[itsNbLevels-1]=param[itsNbLevels-1];
			itsLevels[itsNbLevels-1]=param[2*itsNbLevels-1];
			//Tri du vecteur des seuils
			std::sort(itsThreshold.begin(), itsThreshold.end());
		}
		break;
	default :
		//Size
/*
		itsNbLevels = (size_param+1)/2;
		if ((size_param+1)%2 != 0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"ERROR : PWC Model Parameters Size not OK (n levels + n-1 threshold)");//Error
		}
		else
*/
		{
			itsLinCoef.resize(itsNbLevels);
			itsLevels.resize(itsNbLevels);
			itsThreshold.resize(itsNbLevels-1);
			
			//Init level
			for (i=0; i<itsNbLevels-1; i++)
			{
				itsLinCoef[i]=0.;
				itsLevels[i]=param[i];
				itsThreshold[i]=param[itsNbLevels+i];
			}
			itsLinCoef[itsNbLevels-1]=0.;
			itsLevels[itsNbLevels-1]=param[itsNbLevels-1];
			//Tri du vecteur des seuils
			std::sort(itsThreshold.begin(), itsThreshold.end());
		}
		break;
	}	
}

void ICM_CorrelFunction::InitVAInter(void)
{
	//Estimation du M et NU : pour l'instant uniquement modèle à 2 et 3 états PWC
	itsM=0.;
	itsNU=0.;

	switch (itsCopulaType)
	{
	case qGAUSSIAN:
		switch(itsCorrelFunctionType)
		{
		case qCAL_PWC_CORREL:
			if (itsNbLevels == 2)
			{
				//Va Interm
				double alpha=itsLevels[0];
				double beta=itsLevels[1];
				double theta=itsThreshold[0];
				itsM = alpha*StdNorDensity(theta) - beta*StdNorDensity(theta);
				double VI = (alpha*alpha) * (NAG_cumul_normal(theta) - theta * StdNorDensity(theta)) 
					+ (beta*beta) * (theta * StdNorDensity(theta) + 1 - NAG_cumul_normal(theta)) - (itsM) * (itsM);
				if (VI>1)
					itsNU = 0;
				else
					itsNU = sqrt(1-VI);
			}
			else if (itsNbLevels == 3)
			{
				//Va Interm
				double alpha=itsLevels[0];
				double beta=itsLevels[1];
				double gamma=itsLevels[2];
				double theta1=itsThreshold[0];
				double theta2=itsThreshold[1];
				itsM = StdNorDensity(theta1) * (alpha - beta) + StdNorDensity(theta2) * (beta - gamma);

				double VI = (alpha*alpha) * (NAG_cumul_normal(theta1) - theta1 * StdNorDensity(theta1));
				VI += (beta*beta) * ( (NAG_cumul_normal(theta2) - NAG_cumul_normal(theta1)) + (theta1 * StdNorDensity(theta1) - theta2 * StdNorDensity(theta2) ));
				VI += (gamma*gamma) * ( theta2 * StdNorDensity(theta2) + ( 1 - NAG_cumul_normal(theta2) ) );
				VI -= (itsM) * (itsM);
				if (VI>1)
					itsNU = 0;
				else
					itsNU = sqrt(1-VI);
			}
			break;
		case qCAL_PWL_CORREL:
			break;
		}
		
	case qSTUDENT:
		break;

	default:
		break;
	}
}

// Return Beta coef given the factor
double ICM_CorrelFunction::BetaCond(const double &factor)
{
	int i=0;
	double value=0.,carre =0.;
	
	// 1 correl level
	if (itsVolatility==-999.)
	{
	if (itsNbLevels==1)
		return itsLinCoef[i]*factor + itsLevels[i];
	
	while ( (factor > itsThreshold[i]) && (i<itsNbLevels-1) )
		i++;

	value = itsLinCoef[i]*factor + itsLevels[i];
	}
	else
	{
		if (itsNbLevels==1)
			return itsLinCoef[i]*factor + itsLevels[i];
	
		while ( (factor > itsThreshold[i]) && (i<itsNbLevels-1) )
			i++;

		if (i)
		{	
		carre = (factor-itsThreshold[i-1]);
		carre *= carre;

		value = itsLevels[i] + (itsLevels[i-1]-itsLevels[i]) 
			* exp(-carre/(2.*itsVolatility*itsVolatility));
		}
		else
		{
		value = itsLevels[0];
		}
	}

	return (value);

}


// Return the Marginal Default Probabilitie
double ICM_CorrelFunction::DiffMarginalDefprob(const double &seuil)
{
	double a=0;
	double b=0;
	double rho=0;
	double res=0.;
	
	switch (itsCopulaType)
	{
	case qGAUSSIAN:
		switch(itsCorrelFunctionType)
		{
		case qCAL_PWC_CORREL:
			if (itsNbLevels == 2)
			{
				double alpha=itsLevels[0];
				double beta=itsLevels[1];
				double theta=itsThreshold[0];

				a = (seuil - itsM) / (sqrt(itsNU*itsNU+alpha*alpha));
				b = theta;
				rho = alpha / (sqrt(itsNU*itsNU+alpha*alpha));

				res += NAG_bivariate_normal_dist(a,b,rho);

				a = (seuil - itsM) / (sqrt(itsNU*itsNU+beta*beta));
				res += NAG_cumul_normal(a);
				
				rho = beta / (sqrt(itsNU*itsNU+beta*beta));
				res += -NAG_bivariate_normal_dist(a,b,rho);
			}
			else if (itsNbLevels == 3)
			{
				//Va Interm
				double alpha=itsLevels[0];
				double beta=itsLevels[1];
				double gamma=itsLevels[2];
				double theta1=itsThreshold[0];
				double theta2=itsThreshold[1];
			
				//First part : [-inf, theta1]
				a = (seuil - itsM) / (sqrt(itsNU*itsNU+alpha*alpha));
				b = theta1;
				rho = alpha / (sqrt(itsNU*itsNU+alpha*alpha));

				res += NAG_bivariate_normal_dist(a,b,rho);

				//Seconv part : [theta1,theta2]
				a = (seuil - itsM) / (sqrt(itsNU*itsNU+beta*beta));
				rho = beta / (sqrt(itsNU*itsNU+beta*beta));
				res += -NAG_bivariate_normal_dist(a,b,rho);
				
				b = theta2;
				res += NAG_bivariate_normal_dist(a,b,rho);

				//Third part : [theta2, +inf]
				a = (seuil - itsM) / (sqrt(itsNU*itsNU+gamma*gamma));
				res += NAG_cumul_normal(a);

				rho = gamma / (sqrt(itsNU*itsNU+gamma*gamma));
				res += -NAG_bivariate_normal_dist(a,b,rho);
			}
			break;
		case qCAL_PWL_CORREL:
			break;
		}
		
	case qSTUDENT:
		break;

	default:
		break;
	}
	return itsMarginalProb - res;
}

// Find the Normal threshold given marginal prob
double ICM_CorrelFunction::EstimeSeuil(const double &MarginalDefprob)
{
	//Méthode Root Finder : on estime le seuil qui donne la proba marginal recherchée
	double seuil=0.;
	SetMarginalProb(MarginalDefprob);
	seuil = RootFinder1D(ff1::mem_call(&ICM_CorrelFunction::DiffMarginalDefprob,(*this))).Dichotomy(-6,6,100,1.E-5,100.);
	return seuil;
}