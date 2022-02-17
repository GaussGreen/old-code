
#ifndef _ICM_SuperSeniorLeverage_H_
#define _ICM_SuperSeniorLeverage_H_

/*********************************************************************************/
/*! \class  ICM_SuperSeniorLeverage ICM_SuperSeniorLeverage.h "ICM_SuperSeniorLeverage.h"
 *  \author L Jacquel
 *	\version 1.0
 *	\date   September 2005
 *	\file   ICM_SuperSeniorLeverage.h
 *		\brief Creates a Super Senior Leverage 
/***********************************************************************************/

#include "ICMKernel/util/icm_qmatrix.h"

// #include "ICMKernel/pricer/icm_pricer.h"
//#include "ICMKernel/inst/icm_mez.h"

class ICM_Correlation; 
class ICM_ModelMultiCurves;
class ICM_Pricer;  
class ICM_Mez ;

template  <class ARM_Vector> class ICM_Matrix ; 

enum	SSL_PricingType{SSL_PV, SSL_DEF_LEG, SSL_DEF_LEG_ZC, SSL_SPREAD};
enum	SSL_PricingMethod{SSL_PM_DICHOTOMY, SSL_PM_BRENT};
enum	SSL_ComputationChoice{SSL_CC_MULTIPLES, SSL_CC_TRIGGERS};
enum	SSL_GetData{SSL_GD_NONE, SSL_GD_MULTIPLES, SSL_GD_TRIGGERS};

class ICM_SuperSeniorLeverage : public ARM_Object
{

	private:
	
		ICM_Pricer* its_Pricer;
	
		// PARAMETERS from EXCEL
		double	its_leverage;
		
		double	its_trigger_pct;		// x% for MtM Delta
		double	its_trigger_tol;	// tolerance in the dichotomy

		double	its_target_spread;
		double	its_initial_mtm;		// should be 0.0

		double	its_loss_level;
		double	its_loss_step;
		
		double	its_maturity_step;

		double	its_min_down_multiple;
		double	its_max_down_multiple;

		double	its_min_high_multiple;
		double	its_max_high_multiple;

		int		its_cutoff_loss;
		int		its_cutoff_maturity;

		double	its_loss_trigger;
		double	its_correlation_triggered;

		SSL_PricingType			its_Pricing_Type;
		SSL_PricingMethod		its_Pricing_Method;
		SSL_ComputationChoice	its_Computation_Type;

		SSL_GetData		its_GetData_Type;

		// INTERNAL PARAMETERS
		int	its_NbCredits;
		ARM_Date	its_ValDate;
		ENUM_CMPMETH	its_PricingMode;
		double		its_Tranche_Size;
		ICM_ModelMultiCurves*	its_MultiCurvesModel;
		
		ICM_Correlation* its_Triggered_Correlation;

		double	its_NbMat;
		double	its_NbLosses;

		double	its_Notional;
		double	its_Init_Detach;
		double	its_Init_Attach;

		int		its_Nb_Flags_Rows;
		int		its_Nb_Flags_Cols;

		FILE	*its_fOut;

	void Init()
	{
		SetName(ICM_LEVERAGE);

		its_Pricer	=	NULL;
		its_Matrix_Multiples	=	NULL;
		
		Vect_of_MaturitiesInYF.clear();

		its_leverage	=	0.0;
		its_trigger_pct	=	0.0;
		its_trigger_tol	=	0.0;
		its_target_spread	=	0.0;
		its_loss_level	=	0.0;
		its_loss_step	=	0.0;
		its_maturity_step	=	0.0;;

		its_initial_mtm	=	0.0;

		its_min_down_multiple	=	0.0;
		its_max_down_multiple	=	0.0;
		its_min_high_multiple	=	0.0;
		its_max_high_multiple	=	0.0;

		its_cutoff_loss		=	0;
		its_cutoff_maturity	=	0;

		its_loss_trigger			=	0.0;
		its_correlation_triggered	=	0.0;;

		its_Pricing_Type	=	SSL_PV;
		its_Pricing_Method	=	SSL_PM_DICHOTOMY;
		its_Computation_Type	=	SSL_CC_MULTIPLES;

		its_MultiCurvesModel	=	NULL;
		its_NbCredits	=	0;
		its_ValDate	=	"01/01/2001";
		its_PricingMode	=	qCMPDEFLEGPV;
		its_Tranche_Size	=	0.0;
		its_fOut	= NULL;
		its_Initial_Spreads	=	NULL;

		its_Min_Matrix_Multiples	=	NULL;
		its_Max_Matrix_Multiples	=	NULL;

		its_NbMat		=	0;
		its_NbLosses	=	0;

		its_Notional	=	0.0;
		its_Init_Detach	=	0.0;
		its_Init_Attach	=	0.0;

		its_Nb_Flags_Rows	=	0;
		its_Nb_Flags_Cols	=	0;

		its_VectofProducts.clear();
		its_Vect_of_Maturities.clear();
		its_Vect_of_MaturitiesInYF.clear();
		
		its_Vect_of_Used_Losses.clear();
		its_Vect_of_Used_MaturitiesInYF.clear();

		its_Matrix_Multiples_Flags			=	NULL;
		its_Matrix_Multiples_Products		=	NULL;
		its_Min_Matrix_Control_Multiples	=	NULL;
		its_Max_Matrix_Control_Multiples	=	NULL;
		its_Min_Matrix_Control_Values		=	NULL;
		its_Max_Matrix_Control_Values		=	NULL;

		its_Matrix_Triggers		=	NULL;
		its_Matrix_Flags		=	NULL;
		its_Matrix_MultiplesInput	=	NULL;

		its_Triggered_Correlation	=	NULL;
	}


public:

	ICM_SuperSeniorLeverage() {Init();}
	
	ICM_SuperSeniorLeverage(ICM_SuperSeniorLeverage& data);
	
	~ICM_SuperSeniorLeverage() ;

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object * src) ;



	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
 
		BitwiseCopy(src) ;
	}

	// --------------
	//	Clone Method
	// --------------

	ARM_Object * Clone(void)
	{
		ICM_SuperSeniorLeverage * theClone = new ICM_SuperSeniorLeverage() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

	// --------------------------------------------------------------------
	// DATA
	// --------------------------------------------------------------------

	// Set Cash Flows Matrix
	void	SetDataParameters(ICM_Matrix<ARM_Vector>* parameters);

	inline	ICM_Pricer* GetPricer() {return its_Pricer;}
	void	SetPricer(ICM_Pricer* data) ;

	double	Get_Trigger_Pct()	{return its_trigger_pct;}
	void	Set_Trigger_Pct(double data)	{its_trigger_pct	=	data;}

	double	Get_Trigger_Tolerance()	{return its_trigger_tol;}
	void	Set_Trigger_Tolerance(double data)	{its_trigger_tol	=	data;}

	double	Get_Leverage()	{return its_leverage;}
	void	Set_Leverage(double data)	{its_leverage	=	data;}

	double	Get_Target_Spread()	{return its_target_spread;}
	void	Set_Target_Spread(double data)	{its_target_spread	=	data;}

	double	Get_Loss_Level()	{return its_loss_level;}
	void	Set_Loss_Level(double data)	{its_loss_level	=	data;}

	double	Get_Loss_Step()	{return its_loss_level;}
	void	Set_Loss_Step(double data)	{its_loss_level	=	data;}

	double	Get_Maturity_Step()	{return its_maturity_step;}
	void	Set_Maturity_Step(double data)	{its_maturity_step	=	data;}

	double	Get_Min_Down_Multiple()	{return its_min_down_multiple;}
	void	Set_Min_Down_Multiple(double data)	{its_min_down_multiple	=	data;}

	double	Get_Max_Down_Multiple()	{return its_max_down_multiple;}
	void	Set_Max_Down_Multiple(double data)	{its_max_down_multiple	=	data;}

	double	Get_Min_High_Multiple()	{return its_min_high_multiple;}
	void	Set_Min_High_Multiple(double data)	{its_min_high_multiple	=	data;}

	double	Get_Max_High_Multiple()	{return its_max_high_multiple;}
	void	Set_Max_High_Multiple(double data)	{its_max_high_multiple	=	data;}

	int		Get_CutOff_Loss()	{return its_cutoff_loss;}
	void	Set_CutOff_Loss(double data)	{its_cutoff_loss	=	data;}

	double	Get_CutOff_Maturity()	{return its_cutoff_maturity;}
	void	Set_CutOff_Maturity(double data)	{its_cutoff_maturity	=	data;}

	double	Get_Initial_Mtm()	{return its_initial_mtm;}
	void	Set_Initial_Mtm(double data)	{its_initial_mtm	=	data;}

	SSL_PricingType	Get_PricingType()	{return its_Pricing_Type;}
	void	Set_PricingType(SSL_PricingType data)	{its_Pricing_Type	=	data;}

	SSL_PricingMethod	Get_PricingMethod()	{return its_Pricing_Method;}
	void	Set_PricingMethod(SSL_PricingMethod data)	{its_Pricing_Method	=	data;}

	SSL_ComputationChoice	Get_ComputationType()	{return its_Computation_Type;}
	void	Set_ComputationType(SSL_ComputationChoice data)	{its_Computation_Type	=	data;}


	int		Get_Nb_Flags_Rows()	{return its_Nb_Flags_Rows;}
	void	Set_Nb_Flags_Rows(double data)	{its_Nb_Flags_Rows	=	data;}

	int		Get_Nb_Flags_Cols()	{return its_Nb_Flags_Cols;}
	void	Set_Nb_Flags_Cols(double data)	{its_Nb_Flags_Cols	=	data;}

	
	void Set_Matrix_Flags(ICM_QMatrix<double>* value) ;
	 

	void Set_Multiples_Input(ICM_QMatrix<double>* value) ;
 

	virtual void SetTriggeredCorrelation(ICM_Correlation* correl) ;

	// --------------------------------------------------------------------
	// SUPER SENIOR LEVERAGE COMPUTATION
	// --------------------------------------------------------------------

	public:

		void	Get_Vect_of_Used_Losses(vector<double>& data) {data = its_Vect_of_Used_Losses;}
		void	Get_Vect_of_Used_MaturitiesInYF(vector<double>& data) {data = its_Vect_of_Used_MaturitiesInYF;}

		ICM_QMatrix<double>* GetMatrix_Outputs() ;
		 

		ICM_QMatrix<double>* GetMatrix_Multiples() {return its_Matrix_Multiples;}
		ICM_QMatrix<double>* GetMatrix_Triggers() {return its_Matrix_Triggers;}

	private:
		// INTERNAL DATA

		// storage
		ICM_QMatrix<double>*	its_Initial_Spreads;
		ICM_QMatrix<double>*	its_Matrix_Multiples;
		vector<double>			Vect_of_MaturitiesInYF;

		ICM_QMatrix<double>*	its_Min_Matrix_Multiples;
		ICM_QMatrix<double>*	its_Max_Matrix_Multiples;

		vector<ICM_Mez*>		its_VectofProducts;
		vector<ARM_Date>		its_Vect_of_Maturities;
		vector<double>			its_Vect_of_MaturitiesInYF;
		
		// used
		vector<double>			its_Vect_of_Used_Losses;
		vector<double>			its_Vect_of_Used_MaturitiesInYF;

		ICM_QMatrix<bool>*		its_Matrix_Multiples_Flags;
		ICM_QMatrix<ICM_Mez*>*	its_Matrix_Multiples_Products;
		ICM_QMatrix<double>*	its_Min_Matrix_Control_Multiples;
		ICM_QMatrix<double>*	its_Max_Matrix_Control_Multiples;
		ICM_QMatrix<double>*	its_Min_Matrix_Control_Values;
		ICM_QMatrix<double>*	its_Max_Matrix_Control_Values;

		ICM_QMatrix<double>*	its_Matrix_Triggers;

		// --------------------------------------------------
		ICM_QMatrix<double>*	its_Matrix_Flags;
		ICM_QMatrix<double>*	its_Matrix_MultiplesInput;

	public:

		double	ComputeLeverage();

		// void	TestCalibrationWithGuess();

	private:
	
		void	ComputeMultipleRanges(bool IsDecreasingMaturity, bool IsIncresingLoss, int iloss, int iMat, int NbLosses, int NbMat, double& min_multiple, double& max_multiple);
		void	GetOneMultiple(double min_multiple, double max_multiple, int& nbiter, double& multiple);

		void	UpdateAllPrices(int curr_iloss, int curr_iMat, double multiple);

		void	PrepareProductForPricing(int iloss, int iMat);
		void	UpdateMaxMultiples(int from_Mat, int to_Mat, int from_Loss, int to_Loss, double max_multiple);
		void	UpdateMinMultiples(int from_Mat, int to_Mat, int from_Loss, int to_Loss, double min_multiple);


		void	GetOneMultiple_Brent(double min_multiple, double max_multiple, int& nbiter, double& multiple, double& f_multiple);
		double	ZBrent(double a, double fa, double b, double fb, double c, double fc) const;
		void	GetOnePrice(double multiple, double& result, bool Imposed_NPV_Flag = false);

		void	ComputeTriggersFromMultiples();
};


#endif 
