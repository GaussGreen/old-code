#ifndef _ICM_PRICER_MC_CDO2_H
#define _ICM_PRICER_MC_CDO2_H


#include "ICMKernel/pricer/icm_pricer_mc.h"
#include <set>

typedef std::vector<int> IntVector ;
typedef std::vector<IntVector> IntMatrix ;

class ICM_DefaultCurve; 
class ICM_ModelMultiCurves; 
/*********************************************************************************/
/*! \class  ICM_Pricer_MC icm_pricer_mc.h "icm_pricer_mc.h"
 *  \author Ruben Marciano 
 *	\version 1.0
 *	\date   April 2004
 *	\brief  generate random variables or default times  */
/***********************************************************************************/


class Y
{
	public :
		int id ;
		double tau ;
	public :
		Y(int i, double d) : 
          id(i), tau(d)
		{}
		bool operator < (const Y & rhs) const {
			return tau < rhs.tau ;
		}
} ;



class ICM_Pricer_MC_Cdo2 : public ICM_Pricer_MC
{
	
	private :

	std::vector<std::string> itsUnionNames;
	ICM_QMatrix<int> * itsField ;				// matrice d'appartenance (contient des 1 et des 0)  
	ICM_QMatrix<double>* itsCollatLoss ;		// matrice des losses

	std::set<Y>	itsDefaultTimes_Standart;			// temps de défauts ordonnés

	IntMatrix itsIndCDO ;						// correspondance nom <-> numéro de tranche 
	std::vector<double> itsZcPay ;					// zéro-coupons aux dates de paiement
	std::vector<double> itsSpreads ;                 // Spreads du CDO2 
	std::vector<double> itsAccrualDates ;            // YF AccrualDates du CDO2
	std::vector<double> itsCouponsDates ;            // YF Coupons du CDO2
	std::vector<double> itsYFPayDates ;              // YFPayDates du CDO2
	std::vector<double> itsStrikesDown ;             // vecteur des sub amount
	std::vector<double> itsStrikesUp ;               // vecteur des sub amount + nominal
	std::vector<double> itsBarriers_Standart;		// vecteur des N(-1) (Fi(T))

	std::vector<double> itsPortfolioLosses ;			// vecteur des losses pour chaque portefeuille
	std::vector<double> itsTrancheLosses ;			// vecteur des losses pour haque tranche
	double itsStrikeDown ;
	double itsStrikeUp ;
	double itsCdo2PortfolioLoss ;               // losse du portefeuille global
	double itsCdo2TrancheLoss ;					// losse de la tranche de tranche
	int itsNbNames ;
	int itsNbTranches ;
	int itsNbPeriods ;

	std::vector<ARM_CLASS_NAME>  itsInstrumentType;
	std::vector<int> its_NTD_NumDef;		            // vecteur des nb de défaults par NTD
	std::vector<int> its_NTD_ActualDef;              // vecteur des nb de défaults activés par NTD

	ICM_QMatrix<double>  itsCholeskyMatrix ;	// matrice de cholesky

	int its_sz_dt_i;
	int its_sz_dt_j;

	//Fast Sensitivity Computing
	// 17783 ICM_QMatrix<set<Y>*>* itsDefaultTimes ;		// temps de défauts ordonnés
	// 17783 ICM_QMatrix<vector<double>*>* itsBarriers ;	// vecteur des N(-1) (Fi(T))
	// 17783 ICM_QMatrix<ICM_DefaultCurve*>* itsPerturbDefaultCurves; //Courbes de défaults perturbées
	// 17783 std::vector<int>	itsCorrespondanceUnionNames_Sensis; //Correspondance entre UnionNames & Sensis


	// 17783 ARM_Model*	itsModel_Temp;

	void Init() ;

public :

	ICM_Pricer_MC_Cdo2() {Init();}
 
	void Set(ARM_Security * sec, ARM_Model * mod, const ICM_Parameters& params,const ARM_Date&asof,int nbpath) 
	{

		ICM_Pricer_MC::Set(sec, mod, params,asof,nbpath) ; 
	}

	void BeforePrice(ICM_ModelMultiCurves* Model,qSENSITIVITY_TYPE type = ICMSPREAD_TYPE) ; 
	

	void CptIntermediate(qSENSITIVITY_TYPE typesensi); 
 
	void SetTrancheLosses(double l)
	{
		itsTrancheLosses.resize(itsNbTranches) ;
		for(int i = 0; i < itsNbTranches ; i++)
			itsTrancheLosses[i] = 0. ;
	}

	void SetPortfolioLosses(double l)
	{
		itsPortfolioLosses.resize(itsNbTranches) ;
		for(int i = 0; i < itsNbTranches ; i++)
			itsPortfolioLosses[i] = 0. ;
	}

	void SetCdo2PortfolioLoss(double l)
	{
		itsCdo2PortfolioLoss = l ;
	}

	void SetCdo2TrancheLoss(double l)
	{
		itsCdo2TrancheLoss = l ;
	}

	void ComputeZc() ;
	void ComputeSpreads() ;
	void ComputeAccrualDates() ;
	void ComputeYFPayDates() ;

	void ComputeUnionNames() ;

// 17783 	 	void CptCorrespondanceUnionNames_Sensis() ;
	void ComputeBarriers_Complete() ;

	void GenerateDefaultTimes_Betas_Complete();
	void GenerateDefaultTimes_Betas_Standart();
	void GenerateDefaultTimes_troughCorrelMatrix_Standart();

// 17783 	void PerturbDefaultCurves(void);

	void ComputeFieldAndCollatLoss() ;

	void ComputeIndCDO() ;

	int Eta(const double& tau) ;

	virtual double ComputePrice(qCMPMETH measure );

	virtual double Accrued() { return (0.);}

	double ComputePrice_Basket(qCMPMETH measure );
	//JLA double ComputePrice_Basket_Complete(const double& initialprice);

	virtual double ComputeSpread(const double& MtM =0.);

	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double  epsilonGamma = 0); 
	public:

		

	double ComputeDTR(int NbDefaults, char* S_or_L, char** labels, double* Recoveries) ;
	
	~ICM_Pricer_MC_Cdo2() ;

	void View(char* id, FILE* ficOut);

	void ResetPricer(void)  ;


// 17783 		void ResetSensiManager() ;


	void AffectLossForCrossSub(const int& NoName,
								const int& NoCdo,
								ICM_QMatrix<double>* CollatLoss,
								const vector<double>& Strike_Down,
								const vector<double>& Strike_Up,
								vector<double>& Loss,
								const int& FirstNoCdo  = CREDIT_DEFAULT_VALUE,
								const double& residu = CREDIT_DEFAULT_VALUE);

	double ComputePrice_Basket_CrossSub(qCMPMETH measure);
};


// ----------------------------------------------------------------------------
// Generation of default times with Betas
// ----------------------------------------------------------------------------

#endif
