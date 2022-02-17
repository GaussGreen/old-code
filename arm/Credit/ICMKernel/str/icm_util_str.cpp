/// include this to avoid annoying warnings
#include "ARMKernel/glob/firsttoinc.h"

#include "icm_util_str.h"
#include "icm_MathSrv.h"
#include "ICMKernel\util\icm_fft.h"
#include "ICMKernel\util\icm_ComplexHermiteIntegration.h"
#include "ARMKernel\glob\linalg.h"

// #include "ICMKernel/str/icm_fft_distrib.h"
#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel/crv/icm_default_str.h"
#include "ICMKernel/util/icm_gaussian.h"
#include "ICMKernel/glob/icm_corrmatrix.h"
#include "ICMKernel\\random\icm_RandomRNG-Str.h"
#include "ICMKernel\\random\icm_RandomKISS.h"
#include "ICMKernel\\random\icm_randomRanDef.h"
#include "ICMKernel\\random\ICM_RandomInvNormAcklam.h"

#include <time.h>
#include <sys/timeb.h>

// --------------------------------------------------------------------
// Outils
// --------------------------------------------------------------------
void CompleteVector(vector<double>& Vinit, int size)
{
	int taille = Vinit.size();
	if (taille < size)
	{
		double temp = Vinit[0];
		Vinit.resize(size);
		for (int i=0;i<size;i++)
		Vinit[i] = temp;
	}
	if (taille > size)
		Vinit.resize(size);
}


double RapportNotional(double price, double InitNotional, double NewNotional)
{
	double result = 0.;
	double tol = 1e-10;
	if (fabs(InitNotional)<tol)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Init Notional = 0!");
	}

	result = (NewNotional/InitNotional) * price;
	return result;
}

//  Fonction utilisée dans le cas ou les 2 vecteurs VectToConcat et RefVect sont liés (plot d'une courbe et Spread au niveau du plot par exemple)
//  Alors s'il faut concatené le vecteur VectToConcat, il faut également concatené le vecteur RefVect
//  Comme le vecteur RefVect etant souvent définit pour plusieurs courbes, on ne le modifie pas mais on renvoie son equivalent concatené
vector<double> ConcatVects(vector<double>& VecttoConcat, vector<double> RefVect)
{
	vector<double> RefVectConcat;
	vector<double> VectTemp;
	int i =0;
	int sizeVecttoConcat= VecttoConcat.size();
	for (i=0;i<sizeVecttoConcat;i++)
	{
		if (fabs(VecttoConcat[i] - CREDIT_DEFAULT_VALUE)>1E-8)
		{
			VectTemp.push_back(VecttoConcat[i]);
			RefVectConcat.push_back(RefVect[i]);
		}
	}

	if (VectTemp.size() < sizeVecttoConcat)
	{
		VecttoConcat.clear();
		VecttoConcat = VectTemp;
	}

	return RefVectConcat;
}


vector<double> ConcatVects(ARM_Vector*& VecttoConcat, ARM_Vector* RefVect)
{
	vector<double> RefVectConcat;
	vector<double> VectTemp;
	int i =0;
	int sizeVecttoConcat= VecttoConcat->GetSize();
	for (i=0;i<sizeVecttoConcat;i++)
	{
		if (fabs(VecttoConcat->Elt(i) - CREDIT_DEFAULT_VALUE)>1E-8)
		{
			VectTemp.push_back(VecttoConcat->Elt(i));
			RefVectConcat.push_back(RefVect->Elt(i));
		}
	}

	int sizeVectTemp = VectTemp.size();
	if ( sizeVectTemp< sizeVecttoConcat)
	{
		VecttoConcat->Resize(sizeVectTemp);
		for (i=0;i<sizeVectTemp;i++)
			VecttoConcat->Elt(i) = VectTemp[i];
	}

	return RefVectConcat;
}

// Gestion des indices d'une List
bool IsIn(int indice, vector<int> List)
{
	int i =0;
	bool isIn = false;
	int sizeList = List.size();
	while ( (indice != List[i]) && (i<sizeList))
		i++;

	if (i < sizeList)
		isIn = true;

	return isIn;
}
// ------------------------------------------------------------------------------------------------------
// Fast Stripping ? 
bool FastStripping(string stripping)
{
	bool result = true;
    transform(stripping.begin(), stripping.end(), stripping.begin(), toupper);
	
	if (stripping == "SLOW")
		result = false;

	return result;		
}
// ------------------------------------------------------------------------------------------------------
// Swap
void Swap(double& a, double& b)
{
	double c = a;
	a = b;
	b = c;
}
// --------------------------------------------------------------------
// IF (TargetValue = ValDefaut) then (TargetValue = Refvalue)
// --------------------------------------------------------------------
void CompleteValue(double Refvalue, double& TargetValue, double ValDefaut)
{
	double tol = 1E-5;
	if (fabs(TargetValue - ValDefaut)<tol)
	{
		TargetValue = Refvalue;
	}
}

// --------------------------------------------------------------------
// Calcule le Nb de Flux
// --------------------------------------------------------------------
int Nb_Flux(double Maturity, double Frequency)
{
	int NbFlows = 0;
	if (CHECK_EQUAL(Frequency,0.))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Frequence NULLE");
	}

	double TmpMaturity = Maturity;
	double YF_Step = 1./Frequency;
	while (TmpMaturity > DB_TOL)
	{
		NbFlows ++;
		TmpMaturity -= YF_Step;
	}
	return NbFlows;
}


int Nb_Flux(double Maturity, double Frequency, double& FirstPeriod)
{
	int NbFlows = 0;
	if (CHECK_EQUAL(Frequency,0.))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Frequence NULLE");
	}

	double TmpMaturity = Maturity;
	double YF_Step = 1./Frequency;
	while (TmpMaturity > DB_TOL)
	{
		NbFlows ++;
		TmpMaturity -= YF_Step;
	}

	// On retourne la valeur de la premiere periode de paiement
	FirstPeriod = TmpMaturity + YF_Step;

	return NbFlows;
}

// ------------------------------------------------------------------------------------------------------
// Cree un vecteur de date correspondant aux dates de paiement
// determinées a l'aide de la maturité de de la Frequence
void GenerateCalendar(double Maturity, double Frequency, vector<double>& calendar)
{
	if (CHECK_EQUAL(Frequency,0.))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Frequence NULLE");
	}

	double FirstPeriod = 0.;
	int NbFlow = Nb_Flux(Maturity,Frequency,FirstPeriod);
	calendar.clear();
	calendar.resize(NbFlow);
	calendar[0] = FirstPeriod;
	for (int i = 1;i<NbFlow;i++)
		calendar[i] = FirstPeriod + i*1./Frequency;
}
// ------------------------------------------------------------------------------------------------------
// Cette fonction Retourne un entier correpsondant au type énuméré définissant la fréquence de paiement d'un CDO
int FindFrequency(double Freq, double Maturity)
{
	int Resultat = 0;
	double tol = 1e-6;
	
	if (fabs(Freq - 1)<tol)
	{	
		Resultat = K_ANNUAL;
		return Resultat;
   	}
	else if (fabs(Freq - 2)<tol)
	{	
		Resultat = K_SEMIANNUAL;
		return Resultat;
	} 
	else if (fabs(Freq - 4)<tol)
	{	
		Resultat = K_QUARTERLY;
		return Resultat;
	} 
	else if (fabs(Freq - 6)<tol)
	{	
		Resultat = K_BIMONTHLY;
		return Resultat;
	} 
	else if (fabs(Freq - 12)<tol)
	{	
		Resultat = K_MONTHLY;
		return Resultat;
	} 
	else if (fabs(Freq - 52)<tol)
	{	
		Resultat = K_WEEKLY;
		return Resultat;
	} 
	else if (fabs(Freq - 365)<tol)
	{	
		Resultat = K_DAILY;
		return Resultat;
	} 
	else if ( (fabs(Freq)<tol) || (fabs(Freq - Maturity)<tol) )
	{	
		Resultat = K_ZEROCOUPON;
		return Resultat;
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Frequence de paiement non autorisée");
	}
}

// --------------------------------------------------------------------
// Retourne element du vecteur inférieur ou égale à value (le vecteur est classé par odre croissant !!)
// --------------------------------------------------------------------
int Find_InfEqual(double value, const ARM_Vector* arm_vect)
{
	int i = 0;
	int size = arm_vect->GetSize();
	double tol = 1E-8;
	while ((i<size) && (fabs(value - arm_vect->Elt(i))>tol) && (value > arm_vect->Elt(i)))
		i++;

	if ( (i==0) && (fabs(value - arm_vect->Elt(i))>tol) )		{return CREDIT_DEFAULT_VALUE;}
	else if ((i<size) && (fabs(value - arm_vect->Elt(i))<tol))	{return i;}
	else														{return i-1;}
}

int Find_InfEqual(double value, vector<double> vect)
{
	int i = 0;
	int size = vect.size();
	double tol = 1E-8;
	while ((i<size) && (fabs(value - vect[i])>tol) && (value > vect[i]))
		i++;

	if ( (i==0) && (fabs(value - vect[i])>tol) )		{return CREDIT_DEFAULT_VALUE;}
	else if ((i<size) && (fabs(value - vect[i])<tol))	{return i;}
	else												{return i-1;}
}

int Find_SupEqual(double value, ARM_Vector* arm_vect)
{
	int i = 0;
	int size = arm_vect->GetSize();
	double tol = 1E-8;
	while ((i<size) && (fabs(value - arm_vect->Elt(i))>tol) && (value > arm_vect->Elt(i)))
		i++;

	if ( (i==size) && (fabs(value - arm_vect->Elt(i))>tol) )	{return CREDIT_DEFAULT_VALUE;}
	else if ((i<size) && (fabs(value - arm_vect->Elt(i))<tol))	{return i;}
	else														{return i;}
}

int Find_SupEqual(double value, vector<double> vect)
{
	int i = 0;
	int size = vect.size();
	double tol = 1E-8;
	while ((i<size) && (fabs(value - vect[i])>tol) && (value > vect[i]))
		i++;

	if ( (i==size) && (fabs(value - vect[i])>tol) )		{return CREDIT_DEFAULT_VALUE;}
	else if ((i<size) && (fabs(value - vect[i])<tol))	{return i;}
	else												{return i;}
}

// --------------------------------------------------------------------
// Interpole les spreads au niveau des plots
// --------------------------------------------------------------------
ICM_QMatrix<double>* InterpoleSpread(vector<double> YF_plot,	// Dates (YF) auxquelles on souhaite calculer les Prob de défaut
									 vector<double> Spread,	// Spread des noms au temps T
									 double DateSpread,		// Plot correspondant au spread
									 double pente)			// Pente annuelle en % permettant d'interpoler les spreads au niveau des autres plots
{
	int i,j = 0;
	int Nb_issuer = Spread.size();
	int Nb_plot   = YF_plot.size();
	double deltaYF = 0.;

	ICM_QMatrix<double>* Spread_Mat = new ICM_QMatrix<double>(Nb_issuer,Nb_plot,0.);
	
	// On remplit la matrice
	for (j=0;j<Nb_plot;j++)
	{
		deltaYF = YF_plot[j] - DateSpread;

		for (i=0;i<Nb_issuer;i++)
			Spread_Mat->SetValue(i,j,Spread[i] * pow((1 + pente),deltaYF));
	}

	return Spread_Mat;
}

ARM_Vector* InterpoleOneSpread(ARM_Vector* YF_plot,		
							   double SpreadInit,
							   double DateSpread,			
							   double pente)
{
	int i = 0;
	int size = YF_plot->GetSize();
	ARM_Vector* Spread = new ARM_Vector(size,0.);
	double deltaYF = 0.;
	for (i=0;i<size;i++)
	{
		deltaYF = YF_plot->Elt(i) - DateSpread;
		Spread->Elt(i) = SpreadInit*pow(1+pente,deltaYF);
	}

	return Spread;
}



/** 
void InterpolLineaire(vector<double> Xref,
					  vector<double> Yref,
					  vector<double> Xtarget,
					  vector<double>& Ytarget)
{
	int sizeXref = Xref.size();
	int sizeYref = Yref.size();
	if (sizeYref != sizeXref )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Size of Y and X reference are not compatible");
	}
	
	int sizeXtarget = Xtarget.size();
	Ytarget.clear();
	Ytarget.resize(sizeXtarget);
	
	// On recupere l'indice du Xref <= Xtarget
	vector<int> PosX;
	PosX.resize(sizeXtarget);

	for (int i = 0;i<sizeXtarget;i++)
		PosX[i] = Find_InfEqual(Xtarget[i],Xref);
	
	// On calcule le Ytarget	
	for (i=0;i<sizeXtarget;i++)
	{
		if (PosX[i] == -999)		// Xtarget[i]<Xref[j] pour tout j 
		{
			Ytarget[i] = Yref[0];	// On considere que le Y est constant avant le premier Xref
		}
		else if ( (PosX[i] != -999) && (PosX[i] <sizeYref-1))
		{
			double deltaY = Yref[PosX[i]+1] - Yref[PosX[i]];
			double deltaX = Xref[PosX[i]+1] - Xref[PosX[i]];

			double pente = 0.;
			if (deltaX != 0)
				pente = deltaY/deltaX;

			double cte = Yref[PosX[i]] - pente*Xref[PosX[i]];

			Ytarget[i] = pente * Xtarget[i] + cte;
		}
		else if (PosX[i] == sizeYref-1)
		{
			Ytarget[i] = Yref[sizeYref-1];
		}
	}
}



ARM_Vector* InterpolLineaire (ARM_Vector* Xref,
							  ARM_Vector* Yref,
							  ARM_Vector* Xtarget)
{
	int sizeXref = Xref->GetSize();
	int sizeYref = Yref->GetSize();
	if (sizeYref != sizeXref )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Size of Y and X reference are not compatible");
	}
	
	int sizeXtarget = Xtarget->GetSize();
	if (sizeYref < sizeXtarget )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Size of target > size of reference");
	}
	
	ARM_Vector* Ytarget = new ARM_Vector(sizeXtarget,1.);
	
	// On recupere l'indice du Xref <= Xtarget
	vector<int> PosX;
	PosX.resize(sizeXtarget);

	for (int i = 0;i<sizeXtarget;i++)
		PosX[i] = Find_InfEqual(Xtarget->Elt(i),Xref);
	
	// On calcule le Ytarget	
	for (i=0;i<sizeXtarget;i++)
	{
		if (PosX[i] == -999)		// Xtarget[i]<Xref[j] pour tout j 
		{
			Ytarget->Elt(i) = Yref->Elt(0);	// On considere que le Y est constant avant le premier Xref
		}
		else if ( (PosX[i] != -999) && (PosX[i] <sizeYref-1))
		{
			double deltaY = Yref->Elt(PosX[i]+1) - Yref->Elt(PosX[i]);
			double deltaX = Xref->Elt(PosX[i]+1) - Xref->Elt(PosX[i]);

			double pente = 0.;
			if (deltaX != 0)
				pente = deltaY/deltaX;

			double cte = Yref->Elt(PosX[i]) - pente*Xref->Elt(PosX[i]);

			Ytarget->Elt(i) = pente * Xtarget->Elt(i) + cte;
		}
		else if (PosX[i] == sizeYref-1)
		{
			Ytarget->Elt(i) = Yref->Elt(sizeYref-1);
		}
	}

	return Ytarget;
}
**/ 
// --------------------------------------------------------------------
// Generateur de nombre Aleatoire : Generateur de Lehmer
// --------------------------------------------------------------------
double GenerateurAlea(int& seed)
{
	double resultat = 0.;
	double x,u,a, test = 0.;
	a = seed;

	test = a * 16807./2147483647.;
	if (test > 0)   // Fix = Floor
	{	x = (a * 16807.) - floor(test) * 2147483647;}
	else			// Fix = Floor + 1 
	{	x = (a * 16807.) - (floor(test) + 1) * 2147483647;}

	u = x/2147483647.;
	seed = x;
	resultat = u;

	return resultat;
}



// --------------------------------------------------------------------
// Correl Inter-Tranche (correl evenement de défaut) et Correl Tps de défaut équivalent
// Cette fonction renvoie aussi un vecteur contenant les probas de défaut de chaque tranche et 
// Un vecteur contenant la variance propre a chaque tranche
// Utilisée lorsque l'on étudie un CDO2
// --------------------------------------------------------------------
void MatriceCorrelInterTranche(const long& NbSimul,					
							   const int& NbTranche,
							   const int& NbName,
							   const ARM_Vector& Notionals,
							   const ARM_Vector&  LossRecovery,					
							   const double& BetaIssuers, 
							   ICM_QMatrix<double>* TupTdown,			
							   const ICM_QMatrix<int>& MatAppartenance,	
							   const ARM_Vector&  PdefIssuers,
							   const std::vector<std::string>& labels, 
							   ICM_CorrMatrix*& CorrelEvtDefaut,	
							   ICM_CorrMatrix*& CorrelTpsDefaut,
							   ARM_Vector&  ProbaDefautTranche,
							   ARM_Vector&  VarianceTranche,
							   bool populateCorrelEvtDefaut
							   )
{
	int i = 0;
	int j = 0;

	int TupTdownCol = TupTdown->Getnbcols();
	int TupTdownRow = TupTdown->Getnbrows();

	int MatAppartenanceCol = MatAppartenance.Getnbcols();
	int MatAppartenanceRow = MatAppartenance.Getnbrows();

	double Sqrt_un_moins_beta2 =  sqrt(1-BetaIssuers*BetaIssuers);

	ProbaDefautTranche.Resize(0);
	ProbaDefautTranche.Resize(NbTranche);

	VarianceTranche.Resize(0);
	VarianceTranche.Resize(NbTranche);
	// ---------------------------------------------------------------------------
	// Messages d'erreurs
	if (NbSimul <= 1)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR not enought simulation");
	}
	if (((TupTdownRow%2) != 0) || (TupTdownCol<NbTranche))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR on size of TupTdown array");
	}
	if ((MatAppartenanceRow < NbName) || (MatAppartenanceCol<NbTranche))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR on size of MatAppartenance array");
	}
	if (PdefIssuers.size() < NbName)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR on size of SpreadIssuers vector");
	}
	// ---------------------------------------------------------------------------
	// On stocke toutes les probas de défaut à Maturité
	vector<double> BarriereUnivers;			// Barriere propre à chaque nom. 
	BarriereUnivers.resize(NbName);

	for (i =0;i<NbName;i++)
		BarriereUnivers[i] = ep::MathSrv::invCumNorm(PdefIssuers[i]); 

	
	ICM_QMatrix<double> CorTepsDef_(NbTranche,NbTranche) ; 

// Boucle sur le nombre de simulation
	// On cree les elements dans lesquels on va stocker les données temporaires
	vector<double> SimulAleaUnivers;		// Aléa propre à chaque nom. Modifié pour chaque simulation
	SimulAleaUnivers.resize(NbName);

	vector<double> LossName;				// Loss de chaque nom (1-R si défaut, 0 sinon)
	LossName.resize(NbName);

	vector<double> LossTranche;				// Loss de chaque Tranche
	LossTranche.resize(NbTranche);

	vector<double> PdefCumuleT;				// PDefautTranche cumulées (sur l'ensemble des simulations)
	PdefCumuleT.resize(NbTranche);

	vector<double> PdefSquareCumuleT;		// PDefautTranche² cumulées 
	PdefSquareCumuleT.resize(NbTranche);
 
	ICM_QMatrix<double> CoVariance_cumule_(NbTranche,NbTranche); 


	int seed = 1000;	
	double Uc = 0.;
	double Nc = 0.;
	double Ui = 0.;
	double Ni = 0.;
	double BetaXNc = 0.;
	double sizePtf = 0.;

	// On recupere la taille des tranches pour pouvoir rescaler par la taille de la tranche
	vector<double> SizeTranche;
	SizeTranche.resize(NbTranche);
	for (int l=0;l<NbTranche;l++)
		SizeTranche[l] = TupTdown->Getvalue(1,l) - TupTdown->Getvalue(0,l) ;		// = T_UP - T_Down


	//std::vector<double> allAlea(NbName) ; 
	//ICM_RandomRNG_Str aNewRNG(seed);
	//ICM_RandomRanDef aNewRNG;
	std::vector<double> allAlea(NbName) ;
	//ARM_Vector armAlea(NbName);
	//ICM_RandomInvNormAcklam aNormalGenerator(aNewRNG);
	for (i=0;i<NbSimul;i++)
	{
		// On construit un Alea commun à tous les noms
		Uc = ep::MathSrv::RNG(seed);
		Nc = ep::MathSrv::invCumNorm(Uc);
		BetaXNc = BetaIssuers * Nc; //aNormalGenerator.GenerateOneRandom();
		
		// En sortie de cette boucle on a dans le vecteur LossName la loss de chaque nom 
		// 0 si le nom n'a pas fait défaut
		// 1- R si le nom a fait défaut.
		ep::MathSrv::RNG(seed,&(*allAlea.begin()),allAlea.size()) ; 
		//aNormalGenerator.GenerateRandoms(armAlea);
		long iAlea(0); 
		for (j=0;j<NbName;j++)
		{
			// On reinitialise le vecteur des loss par nom;
			LossName[j] = 0.;

			// JLA Ui = ep::MathSrv::RNG(seed);
			Ni = ep::MathSrv::invCumNorm(allAlea[iAlea++]);
			double SimulAleaUnivers = BetaXNc + Sqrt_un_moins_beta2 * Ni; //armAlea.Elt(j);
			
			// On recupere dans le vecteur LossName la loss correspondant au nom j 
			if (SimulAleaUnivers <BarriereUnivers[j])
				LossName[j] = (1.-LossRecovery[j])*Notionals[j];
			else
				LossName[j] = 0.;
		}
		
		for (int l = 0; l<NbTranche; l++)
		{
			// On reinitialise les données
			sizePtf = 0.;
			double nbnom = 0.;
			LossTranche[l] = 0.;
			// On calcule la loss de la tranche
			for (int k = 0;k<NbName;k++)
			{
				if (MatAppartenance(k,l) == 1)
				{
					LossTranche[l] += LossName[k];
					sizePtf += Notionals[k];
				}
			}

			LossTranche[l] /= sizePtf;

			// On en déduit la proba de défaut de la tranche
			double T_Up = TupTdown->Getvalue(1,l); 
			double Tdown = TupTdown->Getvalue(0,l); 

			if (LossTranche[l]>T_Up)
			{
				ProbaDefautTranche[l] = T_Up - Tdown;
				ProbaDefautTranche[l] /= SizeTranche[l];
			}
			else if (LossTranche[l]>Tdown)
			{
				ProbaDefautTranche[l] = LossTranche[l] - Tdown;
				ProbaDefautTranche[l] /= SizeTranche[l];
			}
			else
				ProbaDefautTranche[l] = 0.;

			PdefCumuleT[l] += ProbaDefautTranche[l];
			PdefSquareCumuleT[l] += ProbaDefautTranche[l]*ProbaDefautTranche[l];
		}

		for (l=0;l<NbTranche;l++)
			for (int k=l+1;k<NbTranche;k++)
				CoVariance_cumule_(l,k)  += ProbaDefautTranche[l]*ProbaDefautTranche[k];			
	}



	// Variance
	for (i=0;i<NbTranche;i++)
	{
		VarianceTranche[i] = (PdefSquareCumuleT[i] - PdefCumuleT[i]*PdefCumuleT[i]/NbSimul)/(NbSimul-1);
		
		if (fabs(VarianceTranche[i])<1E-8)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Variance == 0");
		}

		if (VarianceTranche[i] < 0)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Variance < 0");
		}
	}

	// Proba de défaut tranche
	for (i=0;i<NbTranche;i++)
		ProbaDefautTranche[i] = PdefCumuleT[i] / NbSimul;
	
	
	if (populateCorrelEvtDefaut) 
	{
		ICM_QMatrix<double> CorEvtDef_(NbTranche,NbTranche) ; 
		// On remplit le tableau avec lec correlations, l'esperance et la variace
		for (i=0;i<NbTranche;i++)
		{
			for (j=i+1;j<NbTranche;j++)
			{
 
				CorEvtDef_(i,j) = (CoVariance_cumule_(i,j)/NbSimul - (PdefCumuleT[i]/NbSimul)*
					(PdefCumuleT[j]/NbSimul))/sqrt(VarianceTranche[i]*VarianceTranche[j]);
				CorEvtDef_(j,i) = CorEvtDef_(i,j);

				if (fabs(CorEvtDef_(i,j))>1.)
				{	
					throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"CorEvtDef > 1.");
				}
			}
			CorEvtDef_(i,i)  = 1.;	// Afin d'éviter les erreurs liées a l'utilisation de double
		}

		// On remplit l'objet CorrMatrix qui contiendra la correl Evenement de défaut.
		ARM_Date AsOf;
		string name = "STRCORR";
		CorrelEvtDefaut = new ICM_CorrMatrix(AsOf,name,labels,CorEvtDef_);
	}
	

	double Rho = 0.;
	double test = 0.;
	double Barriere_i = 0.;
	double Barriere_j = 0.;
	double PjointeEquivalente = 0.;
	double ProbJointeTarget = 0.;

	double min = -1.;
	double max = 1.;
	int nbpas = 0;

	for (i=0;i<NbTranche;i++)
	{	
		Rho = 0;
		Barriere_i = ep::MathSrv::invCumNorm(PdefCumuleT[i]/NbSimul);

		for (j=i+1;j<NbTranche;j++)
		{
			max = 1.;
			min = -1.;

			ProbJointeTarget = CoVariance_cumule_(i,j)/NbSimul;
			Barriere_j = ep::MathSrv::invCumNorm(PdefCumuleT[j]/NbSimul);
			
			PjointeEquivalente = BinormalCumule(Barriere_i,Barriere_j,Rho);
			test = ProbJointeTarget - PjointeEquivalente;
			nbpas = 0;
			while ( (fabs(test)>1E-14) && (nbpas <= 100))
			{
				if (PjointeEquivalente>ProbJointeTarget)
					max = Rho;
				else if (PjointeEquivalente<ProbJointeTarget)
					min = Rho;

				Rho = (max + min)/2.;

				PjointeEquivalente = BinormalCumule(Barriere_i,Barriere_j,Rho);
				test = ProbJointeTarget - PjointeEquivalente;

				nbpas ++;
			}
			
			if (nbpas ==  100)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"Dichotomie n'a pas convergé");
			}
			if (fabs(Rho)> 1)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"Correlation temps de défaut > 1");
			}

 			CorTepsDef_(i,j) = Rho;
			CorTepsDef_(j,i) = CorTepsDef_(i,j) ; // [i][j];
		}

 		CorTepsDef_(i,i) = 1.;	// Afin d'éviter les erreurs liées a l'utilisation de double
	}

	// Matrice de correl temps de défaut
	// CorrelTpsDefaut = new ICM_CorrMatrix(AsOf,name,NbTranche,CorTepsDef,Label);
	ARM_Date AsOf;
	string name = "STRCORR";

	CorrelTpsDefaut = new ICM_CorrMatrix(AsOf,name,labels,CorTepsDef_);

 
	// View the Default Times Correlation Matrix 
	#ifdef _DEBUG
		FILE *stream = fopen("c:\\temp\\CorrelDefaultTimes.txt", "w+");
		CorrelTpsDefaut->View("",stream);
		fclose(stream);
	#endif

}


// --------------------------------------------------------------------
// Correl Inter-Tranche (correl evenement de défaut) et Correl Tps de défaut équivalent
// Cette fonction renvoie aussi un vecteur contenant les probas de défaut de chaque tranche et 
// Un vecteur contenant la variance propre a chaque tranche
// Utilisée lorsque l'on étudie un CDO2
// --------------------------------------------------------------------
void MatriceCorrelInterTranche(const long& NbSimul,					
							   const int& NbTranche,
							   const int& NbName,
							   const ARM_Vector & Notionals,
							   const ARM_Vector & IssuersLossRecovery,					
							   const double& BetaIssuers, 
							   const vector<int>& NbDef,			
							   const ICM_QMatrix<int>& MatAppartenance,	
							   const ARM_Vector& PdefIssuers,
							   const std::vector<std::string>& labels,
							   ICM_CorrMatrix*& CorrelEvtDefaut,	
							   ICM_CorrMatrix*& CorrelTpsDefaut,
							   ARM_Vector& ProbaDefautTranche,
							   ARM_Vector& VarianceTranche,
							   bool populateCorrelEvtDefaut)
{
	int i = 0;
	int j = 0;

	int MatAppartenanceCol = MatAppartenance.Getnbcols();
	int MatAppartenanceRow = MatAppartenance.Getnbrows();

	double Sqrt_un_moins_beta2 =  sqrt(1-BetaIssuers*BetaIssuers);

	ProbaDefautTranche.Resize(0);
	ProbaDefautTranche.Resize(NbTranche);

	VarianceTranche.Resize(0);
	VarianceTranche.Resize(NbTranche);
	// ---------------------------------------------------------------------------
	// Messages d'erreurs
	if (NbSimul <= 1)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR not enought simulation");
	}
	if ((MatAppartenanceRow < NbName) || (MatAppartenanceCol<NbTranche))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR on size of MatAppartenance array");
	}
	if (PdefIssuers.size() < NbName)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR on size of SpreadIssuers vector");
	}
	// ---------------------------------------------------------------------------
	// On stocke toutes les probas de défaut à Maturité
	vector<double> BarriereUnivers;			// Barriere propre à chaque nom. 
	BarriereUnivers.resize(NbName);

	for (i =0;i<NbName;i++)
		BarriereUnivers[i] = ep::MathSrv::invCumNorm(PdefIssuers[i]); 

	// On définit les double** qui contiendront les probas de défaut des tranches
 	
	ICM_QMatrix<double> CorTepsDef_(NbTranche,NbTranche); 

// Boucle sur le nombre de simulation
	// On cree les elements dans lesquels on va stocker les données temporaires
	vector<double> SimulAleaUnivers;		// Aléa propre à chaque nom. Modifié pour chaque simulation
	SimulAleaUnivers.resize(NbName);

	vector<double> LossName;				// Loss de chaque nom (1-R si défaut, 0 sinon)
	LossName.resize(NbName);

	// std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ DefaultName;				// Defaut de chaque nom 
	// DefaultName.resize(NbName);

	vector<double> LossTranche;				// Loss de chaque Tranche
	LossTranche.resize(NbTranche);

	// vector<int> NbDefTranche;				// Loss de chaque Tranche
	// NbDefTranche.resize(NbTranche);

	vector<double> PdefCumuleT;				// PDefautTranche cumulées (sur l'ensemble des simulations)
	PdefCumuleT.resize(NbTranche);

	vector<double> PdefSquareCumuleT;		// PDefautTranche² cumulées 
	PdefSquareCumuleT.resize(NbTranche);
	
 
	ICM_QMatrix<double> CoVariance_cumule_(NbTranche,NbTranche); 


 
	int seed = 1000;	
	double Uc = 0.;
	double Nc = 0.;
	double Ui = 0.;
	double Ni = 0.;
	double BetaXNc = 0.;
	double sizePtf = 0.;
	double AvgNot = 0.;
	double AvgRec = 0.;
	double T_Up = 0.;
	double Tdown = 0.;
	
	ICM_RandomRNG_Str aNewRNG(seed);
	//ICM_RandomRanDef aNewRNG;
	//std::vector<double> allAlea(NbName) ;
	ARM_Vector armAlea(NbName);
	ICM_RandomInvNormAcklam aNormalGenerator(aNewRNG);
	for (i=0;i<NbSimul;i++)
	{
		// On construit un Alea commun à tous les noms
		//Uc = ep::MathSrv::RNG(seed);
		//Uc = aNewRNG.GenerateOneRandom();
		//Nc = ep::MathSrv::invCumNorm(Uc);
		
		Nc = aNormalGenerator.GenerateOneRandom();
		BetaXNc = BetaIssuers * Nc;
		
		// En sortie de cette boucle on a dans le vecteur LossName la loss de chaque nom 
		// 0 si le nom n'a pas fait défaut
		// 1- R si le nom a fait défaut.


// FIXMEFRED: mig.vc8 (28/05/2007 15:06:25):cast
		//ep::MathSrv::RNG(seed,&(*allAlea.begin()),allAlea.size()); 
		//ep::MathSrv::RNG(seed,allAlea.begin(),allAlea.size()); 
		//aNewRNG.GenerateRandoms(armAlea);
		aNormalGenerator.GenerateRandoms(armAlea);
		long iAlea(0); 
		for (j=0;j<NbName;j++)
		{
			// On reinitialise le vecteur des loss par nom;
			LossName[j] = 0.;
			//DefaultName[j] = false;

			//JLA Ui = ep::MathSrv::RNG(seed);
			//Ni = ep::MathSrv::invCumNorm(allAlea[iAlea++]);
			//
			
			//Ni = ep::MathSrv::invCumNorm(armAlea[iAlea++]);
			//SimulAleaUnivers[j] = BetaXNc + Sqrt_un_moins_beta2 * Ni;
			SimulAleaUnivers[j] = BetaXNc + Sqrt_un_moins_beta2 * armAlea.Elt(j);
			// On recupere dans le vecteur LossName la loss correspondant au nom j 
			if (SimulAleaUnivers[j]<BarriereUnivers[j])
			{
				LossName[j] = (1.-IssuersLossRecovery[j])*Notionals[j];
				//DefaultName[j] = true;
			}
			else
				LossName[j] = 0.;
		}
		
		for (int l = 0; l<NbTranche; l++)
		{
			// On reinitialise les données
			double nbnom = 0.;
			LossTranche[l] = 0.;
			//NbDefTranche[l] = 0;
			ProbaDefautTranche[l]=0.;
			sizePtf = 0.;
			AvgNot = 0.;
			AvgRec = 0.;

			// On calcule la loss de la tranche
			for (int k = 0;k<NbName;k++)
			{ if (MatAppartenance(k,l)==1) 
				{LossTranche[l] += LossName[k];
				//if (DefaultName[k])	NbDefTranche[l]+=1;		
				sizePtf += Notionals[k];
				nbnom++;
				AvgRec += IssuersLossRecovery[k];
				}
			}

			LossTranche[l] /= sizePtf;
			AvgNot = sizePtf/(double)nbnom;
			AvgRec /= (double)nbnom;

			// On en déduit la proba de défaut de la tranche Version 1
			//if (NbDefTranche[l]>=NbDef[l]) {ProbaDefautTranche[l]=1.;}

			T_Up = (1.-AvgRec)*AvgNot*(double)(NbDef[l])/sizePtf;
			Tdown = (1.-AvgRec)*AvgNot*(double)(NbDef[l]-1)/sizePtf;

			if (LossTranche[l]>T_Up)
			{
				ProbaDefautTranche[l] = 1.;
			}
			else if (LossTranche[l]>Tdown)
			{
				ProbaDefautTranche[l] = LossTranche[l] - Tdown;
				ProbaDefautTranche[l] /= (T_Up-Tdown);
			}


			PdefCumuleT[l] += ProbaDefautTranche[l];
			PdefSquareCumuleT[l] += ProbaDefautTranche[l]*ProbaDefautTranche[l];
		}

		for (l=0;l<NbTranche;l++)
			for (int k=l+1;k<NbTranche;k++)
				CoVariance_cumule_(l,k) += ProbaDefautTranche[l]*ProbaDefautTranche[k];			
	}


	// Variance
	for (i=0;i<NbTranche;i++)
	{
		VarianceTranche[i] = (PdefSquareCumuleT[i] - PdefCumuleT[i]*PdefCumuleT[i]/NbSimul)/(NbSimul-1);
		
		if (fabs(VarianceTranche[i])<1E-8)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Variance == 0");
		}

		if (VarianceTranche[i] < 0)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Variance < 0");
		}
	}

	// Proba de défaut tranche
	for (i=0;i<NbTranche;i++)
		ProbaDefautTranche[i] = PdefCumuleT[i] / NbSimul;
	
	
	if (populateCorrelEvtDefaut) 
	{
		ICM_QMatrix<double> CorEvtDef_(NbTranche,NbTranche); 
		// On remplit le tableau avec lec correlations, l'esperance et la variace
		for (i=0;i<NbTranche;i++)
		{
			for (j=i+1;j<NbTranche;j++)
			{
				CorEvtDef_(i,j) = (CoVariance_cumule_(i,j)/NbSimul - (PdefCumuleT[i]/NbSimul)*
					(PdefCumuleT[j]/NbSimul))/sqrt(VarianceTranche[i]*VarianceTranche[j]);
				CorEvtDef_(j,i) = CorEvtDef_(i,j);

				if (fabs(CorEvtDef_(i,j))>1.)
				{	
					throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"CorEvtDef > 1.");
				}
			}
			CorEvtDef_(i,i) = 1.;	// Afin d'éviter les erreurs liées a l'utilisation de double
		}

		// On remplit l'objet CorrMatrix qui contiendra la correl Evenement de défaut.
		ARM_Date AsOf;
		string name = "STRCORR";
		CorrelEvtDefaut = new ICM_CorrMatrix(AsOf,name,labels,CorEvtDef_);
	} 

	double Rho = 0.;
	double test = 0.;
	double Barriere_i = 0.;
	double Barriere_j = 0.;
	double PjointeEquivalente = 0.;
	double ProbJointeTarget = 0.;

	double min = -1.;
	double max = 1.;
	int nbpas = 0;

	for (i=0;i<NbTranche;i++)
	{	
		Rho = 0;
		Barriere_i = ep::MathSrv::invCumNorm(PdefCumuleT[i]/NbSimul);

		for (j=i+1;j<NbTranche;j++)
		{
			max = 1.;
			min = -1.;

			ProbJointeTarget = CoVariance_cumule_(i,j)/NbSimul;
			Barriere_j = ep::MathSrv::invCumNorm(PdefCumuleT[j]/NbSimul);
			
			PjointeEquivalente = BinormalCumule(Barriere_i,Barriere_j,Rho);
			test = ProbJointeTarget - PjointeEquivalente;
			nbpas = 0;
			while ( (fabs(test)>1E-14) && (nbpas <= 100))
			{
				if (PjointeEquivalente>ProbJointeTarget)
					max = Rho;
				else if (PjointeEquivalente<ProbJointeTarget)
					min = Rho;

				Rho = (max + min)/2.;

				PjointeEquivalente = BinormalCumule(Barriere_i,Barriere_j,Rho);
				test = ProbJointeTarget - PjointeEquivalente;

				nbpas ++;
			}
			
			if (nbpas ==  100)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"Dichotomie n'a pas convergé");
			}
			if (fabs(Rho)> 1)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
					"Correlation temps de défaut > 1");
			}

			CorTepsDef_(i,j) = Rho;
			CorTepsDef_(j,i) = CorTepsDef_(i,j);
		}

		CorTepsDef_(i,i) = 1.;	// Afin d'éviter les erreurs liées a l'utilisation de double
	}

	// Matrice de correl temps de défaut
	ARM_Date AsOf;
	string name = "STRCORR";
	CorrelTpsDefaut = new ICM_CorrMatrix(AsOf,name,labels,CorTepsDef_);


}


// ------------------------------------------------------------------------------------------------------
// Fonction calculant la year fraction séparant deux dates en prenant en compte la convention souhaitée
double YFfromDate(const ARM_Date& AsOfDate,
				  const ARM_Date& TargetDate,
				  const double& YearConvention)
{
	if (YearConvention <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : YearConvention <= 0");
	}
	
	double result = (TargetDate.GetJulian() - AsOfDate.GetJulian())/YearConvention;

	return result;
}

double YFfromDate(const double& AsOfDate,
				  const double& TargetDate,
				  const double& YearConvention)
{
	if (YearConvention <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : YearConvention <= 0");
	}
	if ( (AsOfDate < 1000) || (TargetDate < 1000))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
		    "ERROR : Not true ExcelDate");
	}

	double result = (TargetDate - AsOfDate)/YearConvention;

	return result;
}

double YFfromDate(const ARM_Date& AsOfDate,
				  const double& TargetDate,
				  const double& YearConvention)
{
	if (YearConvention <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : YearConvention <= 0");
	}
	if (TargetDate < 1000)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"ERROR : Not true ExcelDate");
	}

	// La constante 2415019 est la constante entre les julian et les excel date
	double AsOfExcelDate = AsOfDate.GetJulian() - 2415019;
	double result = (TargetDate - AsOfExcelDate)/YearConvention;

	return result;
}
// ------------------------------------------------------------------------------------------------------
// Convertir une maturité ExcelDate en maturité YF
void MaturityYF(double& maturity,
				const double& JulianAsOfDate,
				const double& Convention)
{
	// Cte séparant les julian des Excel Date
	double Cte = 2415019;
	double ExcelAsOf = 0.;

	if (maturity > 100)
	{
		if (fabs(JulianAsOfDate)<1e-6)
		{
			ARM_Date AsOf;
			AsOf.Today();
			ExcelAsOf = AsOf.GetJulian() - 2415019;
		}
		else
			ExcelAsOf = JulianAsOfDate - 2415019;
		
		maturity -= ExcelAsOf;
		maturity *= 1./Convention;
	}
}
// --------------------------------------------------------------------
// ********************************************************************
// --------------------------------------------------------------------
// Duration d'un CDS
// --------------------------------------------------------------------
double DurationCds(const ICM_DefaultCurve* Dcurve,
				   double Maturity,
				   double PasCds)
{
	if (PasCds <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : PasCds <= 0");
	}
	int j = 0;
	vector<double> calendar;
	GenerateCalendar(Maturity,1./PasCds,calendar);

	vector<double> Psurvie;
	ProbaSurvieCdS(Dcurve,calendar,Psurvie);

	ARM_ZeroCurve* CrbTaux = Dcurve->GetZeroCurve();
	vector<double> DFactor;
	DiscountFactor(CrbTaux,calendar,DFactor);
	double Dur = DurationCds(Psurvie,Maturity,DFactor,PasCds);

	return Dur;
}
// --------------------------------------------------------------------
// FeeLeg d'un CDS
// --------------------------------------------------------------------
double FeeLegCds(ICM_DefaultCurve* Dcurve,
				 double Maturity, 
				 double Spread, 
				 double Nominal, 
				 double PasCds)
{
	int j = 0;
	double FL = DurationCds(Dcurve,Maturity,PasCds);
	
	FL *= Spread * Nominal;

	return FL;
}
// --------------------------------------------------------------------
// DefLeg d'un CDS
// --------------------------------------------------------------------
double DefLegCds(const ICM_DefaultCurve* Dcurve,
				 double Maturity,
				 double Nominal,
				 double PasCds,
				 double RecovLoss)
{
	int j = 0;
	double DL = 0.;
	if (RecovLoss == CREDIT_DEFAULT_VALUE)
		RecovLoss = Dcurve->GetRecovery();
	ARM_ZeroCurve* CrbTaux = Dcurve->GetZeroCurve();

	vector<double> calendar;
	GenerateCalendar(Maturity,1./PasCds,calendar);

	vector<double> PSurvie;
	ProbaSurvieCdS(Dcurve,calendar,PSurvie);
	vector<double> DFactor;
	DiscountFactor(CrbTaux,calendar,DFactor);

	DL = DefLegCds(PSurvie,RecovLoss,DFactor,Nominal);
	return DL;
}
// --------------------------------------------------------------------
// NPV d'un CDS
// --------------------------------------------------------------------
double NPVCds(const ICM_DefaultCurve* Dcurve,
			  double Maturity,
			  double Spread,
			  double Nominal,
			  double PasCds,
			  double RecovLoss)
{	
	if (PasCds <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : PasCds <= 0");
	}
	vector<double> calendar;
	GenerateCalendar(Maturity,1./PasCds,calendar);

	vector<double> PSurvie;
	ProbaSurvieCdS(Dcurve,calendar,PSurvie);

	ARM_ZeroCurve* CrbTaux = Dcurve->GetZeroCurve();
	vector<double> DFactor;
	DiscountFactor(CrbTaux,calendar,DFactor);
	
	if (RecovLoss == CREDIT_DEFAULT_VALUE)
		RecovLoss = Dcurve->GetRecovery();

	double NPV = NPVCds(PSurvie,Maturity,Spread,RecovLoss,DFactor,Nominal,PasCds);
	return NPV;
}
// --------------------------------------------------------------------
// Break Event d'un CDS
// --------------------------------------------------------------------
double SpreadCds(const ICM_DefaultCurve* Dcurve,
				 double Maturity,
				 double PasCds,
				 double RecovLoss)
{
	double DL  = DefLegCds(Dcurve, Maturity, 1., PasCds,RecovLoss);
	double Dur = DurationCds(Dcurve, Maturity, PasCds);
	double SpreadCds = DL/Dur;

	return SpreadCds;
}
// --------------------------------------------------------------------
// SpreadForward
// --------------------------------------------------------------------
double SpreadForward(const ICM_DefaultCurve* Dcurve,
					 double YFDateAsOf,
					 double Duree,
					 double PasCds,
					 double RecovLoss)
{
	double SForward = 0.;
	double S1,S2 = 0.;
	double Dur1, Dur2 = 0.;

	if (YFDateAsOf == 0)
	{
		SForward = SpreadCds(Dcurve, Duree, PasCds);
		return SForward;
	}

	S1 = SpreadCds(Dcurve, YFDateAsOf, PasCds ,RecovLoss);
	S2 = SpreadCds(Dcurve, YFDateAsOf+Duree, PasCds ,RecovLoss);

	Dur1 = DurationCds(Dcurve, YFDateAsOf, PasCds);
	Dur2 = DurationCds(Dcurve, Duree + YFDateAsOf, PasCds);

	SForward = (Dur2 * S2 - Dur1 * S1)/(Dur2 - Dur1);

	return SForward;
}

// --------------------------------------------------------------------
// Duration d'un CDS
// --------------------------------------------------------------------
double DurationCds(const vector<double>& Psurvie,
				   double Maturity,
				   const vector<double>& DFactor,
				   double PasCds)
{
	if (PasCds <= 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : PasCds <= 0");
	}
	int size = Psurvie.size();
	if (size != DFactor.size())
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Size Error <= 0");
	}

	double Dur = 0.;
	double DurInit = 0.;
	double FirstPeriod = Maturity - (size-1)*PasCds;
	
	// Premiere periode
	DurInit = Psurvie[0]*DFactor[0]*FirstPeriod;
	
	for (int j=1;j<size;j++)
		Dur += Psurvie[j]*DFactor[j];
	
	Dur *= PasCds;
	Dur += DurInit;

	return Dur;
}
// --------------------------------------------------------------------
// FeeLeg d'un CDS
// --------------------------------------------------------------------
double FeeLegCds(const vector<double>& Psurvie,
				 double Maturity,
				 double Spread,
				 const vector<double>& DFactor,
				 double Nominal,
				 double PasCds)
{
	int size = Psurvie.size();

	double FL = DurationCds(Psurvie,Maturity,DFactor,PasCds);
	FL *= Spread*Nominal;

	return FL;
}
// --------------------------------------------------------------------
// DefLeg d'un CDS
// --------------------------------------------------------------------
double DefLegCds(const vector<double>& Psurvie,
				 double Recov,
				 const vector<double>& DFactor,
				 double Nominal)
{
	int size = Psurvie.size();
	if (size != DFactor.size())
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Size Error <= 0");
	}
	
	double DL = (1-Psurvie[0])*DFactor[0];

	for (int j = 1;j<size;j++)
		DL += (Psurvie[j-1] - Psurvie[j])*DFactor[j];
	
	DL *= (1-Recov)*Nominal;

	return DL;
}
// --------------------------------------------------------------------
// NPV d'un CDS
// --------------------------------------------------------------------
double NPVCds(const vector<double>& Psurvie,
			  double Maturity,
			  double Spread,
			  double Recov,
			  const vector<double>& DFactor,
			  double Nominal,
			  double PasCds)
{	
	double FL = FeeLegCds(Psurvie, Maturity,Spread, DFactor, Nominal, PasCds);
	double DL = DefLegCds(Psurvie,Recov,DFactor,Nominal);
	double NPV = DL - FL;

	return NPV;
}
// --------------------------------------------------------------------
// Break Event d'un CDS
// --------------------------------------------------------------------
double SpreadCds(const vector<double>& P_survie, double Maturity, double Recov,const vector<double>& DFactor, double PasCds)
{
	double DL  = DefLegCds(P_survie,Recov,DFactor,1.);
	double Dur = DurationCds(P_survie, Maturity,DFactor,PasCds);
	double SpreadCds = DL/Dur;

	return SpreadCds;
}
// ---------------------------------------------------------------------------
// Calcule la proba de survie d'un CDS aux dates situées dans le vecteur Dates
// ---------------------------------------------------------------------------
void ProbaSurvieCdS(const ICM_DefaultCurve* DCurve, const vector<double>& Dates, vector<double>& PSurvie)
{
	int size = Dates.size();
	PSurvie.clear();
	PSurvie.resize(size);
	
	for (int i =0;i<size;i++)
		PSurvie[i] = DCurve->SurvivalProba(Dates[i]);
}
// ---------------------------------------------------------------------------
// Calcule la proba de survie d'un CDS aux dates situées dans le vecteur Dates
// ---------------------------------------------------------------------------
void DiscountFactor(ARM_ZeroCurve* CrbTaux, const vector<double>& Dates, vector<double>& DFactor)
{
	int size = Dates.size();
	DFactor.clear();
	DFactor.resize(size);
	
	for (int i =0;i<size;i++)
		DFactor[i] = CrbTaux->DiscountPrice(Dates[i]);
}
