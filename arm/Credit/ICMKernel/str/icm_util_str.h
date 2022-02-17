#ifndef _UTIL_STR_H_
#define _UTIL_STR_H_

#include <vector>
// #include "ICMKernel/str/icm_fft_distrib.h"
#include "ICMKernel/util/icm_qmatrix.h"
// #include "ICMKernel/glob/icm_corrmatrix.h"
#include <math.h>
// using namespace std;
class ARM_Vector;
class ICM_DefCurvStr;
class ICM_DefaultCurve;
class ICM_CorrMatrix ;
// Fonction lib�rant la m�moire et deletant les objets
template <class T>
inline void MyDelete(T*& object)
{
	if (object)
		delete object;
	object = NULL;
}

template <class T>
inline void MyDelete(T**& object, int size)
{
	if (object)
	{
		for (int i=0;i<size;i++)
		{
			if (object[i])
				MyDelete(object[i]);
		}
		delete[] object;
		object = NULL;
	}
}

template <class T>
inline void MyDeleteTab(T*& object)
{
	if (object)
		delete[] object;
	object = NULL;
}

template <class T>
inline void MyDeleteTab(T**& object, int size)
{
	if (object)
	{
		for (int i=0;i<size;i++)
		{
			if (object[i])
				MyDeleteTab(object[i]);
		}
		delete[] object;
		object = NULL;
	}
}
// ------------------------------------------------------------------------------------------------------
// Fast Stripping ? 
bool FastStripping(string stripping);
// ------------------------------------------------------------------------------------------------------
// Algorithme de Tri rapide
void Swap(double& a, double& b);
// ------------------------------------------------------------------------------------------------------
// Transforme le vecteur Vinit si celui ci
// est de dimension diff�rente de celle pass�e en parametre.
// Alors Vinit devient un vecteur de taille :"size" dont toutes les composantes valent Vinit[0] si Vinit.size < size
// juste amputer des valeurs superflues si Vinit.size > size
void CompleteVector(vector<double>& Vinit, int size);
// ------------------------------------------------------------------------------------------------------
// Renvoie : (NewNotional/InitNotional) * price
double RapportNotional(double price, double InitNotional, double NewNotional);

// ------------------------------------------------------------------------------------------------------
// Concatenation de vecteur. 
//  Fonction utilis�e dans le cas ou les 2 vecteurs VectToConcat et RefVect
//  sont li�s (plot d'une courbe et Spread au niveau du plot par exemple)
//  Alors s'il faut concaten� le vecteur VectToConcat, il faut �galement concaten� le vecteur RefVect
//  Comme le vecteur RefVect etant souvent d�finit pour plusieurs courbes,
//  on ne le modifie pas mais on renvoie son equivalent concaten�
vector<double> ConcatVects(vector<double>& VecttoConcat, vector<double> RefVect);
vector<double> ConcatVects(ARM_Vector*& VecttoConcat, ARM_Vector* RefVect);

// ------------------------------------------------------------------------------------------------------
// Renvoie le nombre de flux lors d'un Deal de Maturity et de Frequency (4. pour quaterly)
int Nb_Flux(double Maturity, double Frequency);
int Nb_Flux(double Maturity, double Frequency, double& FirstPeriod);

// ------------------------------------------------------------------------------------------------------
// Cree un vecteur de date correspondant aux dates de paiement
// determin�es a l'aide de la maturit� de de la Frequence
void GenerateCalendar(double Maturity, double Frequency, vector<double>& calendar);


// ------------------------------------------------------------------------------------------------------
// Retourne un entier correpsondant au type �num�r� d�finissant la fr�quence de paiement d'un CDO
int FindFrequency(double Freq, double Maturity = CREDIT_DEFAULT_VALUE);


// ------------------------------------------------------------------------------------------------------
// IF (TargetValue = ValDefaut) then (TargetValue = Refvalue); 
void CompleteValue(double Refvalue, double& TargetValue, double ValDefaut);

// ------------------------------------------------------------------------------------------------------
// Renvoie l'element du vecteur imm�diatement inf�rieur ou �gale � la value.
// Retourne - 999 si tous les elements du tableau sont sup�rieurs � value
int Find_InfEqual(double value, const ARM_Vector* arm_vect);
int Find_InfEqual(double value, vector<double> vect);
int Find_SupEqual(double value, const ARM_Vector* arm_vect);
int Find_SupEqual(double value, vector<double> vect);

// ------------------------------------------------------------------------------------------------------
// Interpole les spreads aux niveaux des diff�rents plots
ICM_QMatrix<double>* InterpoleSpread(vector<double> YF_plot,	// Dates (YF) auxquelles on souhaite calculer les Prob de d�faut
									 vector<double> Spread,		// Spread des noms au temps T
									 double DateSpread,			// Plot correspondant au spread
									 double pente);				// Pente annuelle en % permettant d'interpoler les spreads au niveau des autres plots


ARM_Vector* InterpoleOneSpread(ARM_Vector* YF_plot,		
							   double SpreadInit,	
							   double DateSpread,			
							   double pente);				





// ------------------------------------------------------------------------------------------------------
// Permet une interpolation lin�aire lorsque l'on connait les valeurs d'une fonction sur un nuage de points.
/**
void InterpolLineaire(vector<double> Xref,
					  vector<double> Yref,
					  vector<double> Xtarget,
					  vector<double>& Ytarget);


ARM_Vector* InterpolLineaire (ARM_Vector* Xref,
							  ARM_Vector* Yref,
							  ARM_Vector* Xtarget);
							  **/ 

// ------------------------------------------------------------------------------------------------------
// Generateur d'Alea
double GenerateurAlea(int& seed);	// Generateur de Lehmer
// ------------------------------------------------------------------------------------------------------
// Correl Inter-Tranche (correl evenement de d�faut) et Correl Tps de d�faut �quivalent
// Utilis�e lorsque l'on �tudie un CDO2

void MatriceCorrelInterTranche(const long& NbSimul,							// nb de simulation
							   const int& NbTranche,						// nb de CDO fille
							   const int& NbName,							// Nbr de nom total 
							   const ARM_Vector& Notionals,			// vecteur de notionels par issuers
							   const ARM_Vector&  LossRecovery,			// taux de recouvrement determinant la loss
							   const double& BetaIssuers,
							   ICM_QMatrix<double>* TupTdown,				// Matrice des strikes des CDO filles
							   const ICM_QMatrix<int>& MatAppartenance,			// Matrice de 0 et de 1 definissant les tranches filles
							   const ARM_Vector&  ProbDefIssuers,		// vecteur des probas de d�faut des noms � maturit�
							   const std::vector<std::string>& labels,
							   ICM_CorrMatrix*& CorrelEvtDefaut,	
							   ICM_CorrMatrix*& CorrelTpsDefaut,
							   ARM_Vector&  ProbaDefautTranche,
							   ARM_Vector&  VarianceTranche,
							   bool populateCorrelEvtDefaut);


void MatriceCorrelInterTranche(const long& NbSimul,					
							   const int& NbTranche,
							   const int& NbName,
							   const ARM_Vector&  Notionals,
							   const ARM_Vector&  IssuersLossRecovery,					
							   const double& BetaIssuers, 
							   const vector<int>& NbDef,			
							   const ICM_QMatrix<int>& MatAppartenance,	
							   const ARM_Vector&  PdefIssuers,
							   const std::vector<std::string>& labels,
							   ICM_CorrMatrix*& CorrelEvtDefaut,	
							   ICM_CorrMatrix*& CorrelTpsDefaut,
							   ARM_Vector&  ProbaDefautTranche,
							   ARM_Vector&  VarianceTranche,
							   bool populateCorrelEvtDefaut);

// ------------------------------------------------------------------------------------------------------
// Fonction calculant la year fraction s�parant deux dates en prenant en compte la convention souhait�e
double YFfromDate(const ARM_Date& AsOfDate,
				  const ARM_Date& TargetDate,
				  const double& YearConvention = 365.);

// On consid�re qu'une excel date est un double n�cessairement sup�rieur � 1000
double YFfromDate(const double& AsOfDate,
				  const double& TargetDate,
				  const double& YearConvention = 365.);

double YFfromDate(const ARM_Date& AsOfDate,
				  const double& TargetDate,
				  const double& YearConvention = 365.);
// ------------------------------------------------------------------------------------------------------
// Convertir une maturit� ExcelDate en maturit� YF
void MaturityYF(double& maturity,
				const double& JulianAsOfDate = 0.,
				const double& Convetion = 365.25); 
// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
// Fonction CDS
double NPVCds(const ICM_DefaultCurve* Dcurve, double Maturity, double Spread, double Nominal = 10000000, double PasCds = 0.25, double RecovLoss = CREDIT_DEFAULT_VALUE);
double SpreadCds(const ICM_DefaultCurve* Dcurve, double Maturity, double PasCds = 0.25, double RecovLoss = CREDIT_DEFAULT_VALUE);
double SpreadForward(const ICM_DefaultCurve* Dcurve, double YFDateAsOf,double Duree, double PasCds = 0.25, double RecovLoss = CREDIT_DEFAULT_VALUE);

double DefLegCds(const vector<double>& Psurvie,  double Recovery, const vector<double>& DFactor,double Nominal = 10000000);
double NPVCds   (const vector<double>& Psurvie,  double Maturity, double Spread, double Recovery,const vector<double>& DFactor, double Nominal = 10000000, double PasCds = 0.25);
double DurationCds(const vector<double>& Psurvie,double Maturity, const vector<double>& DFactor,double PasCds = 0.25);

void ProbaSurvieCdS(const ICM_DefaultCurve* DCurve, const vector<double>& Dates, vector<double>& PSurvie);


void DiscountFactor(ARM_ZeroCurve* CrbTaux, const vector<double>& Dates, vector<double>& DFactor);


#endif