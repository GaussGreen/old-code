#ifndef _ICM_IMPLIED_LOSS_TREE_H
#define _ICM_IMPLIED_LOSS_TREE_H

/*********************************************************************************/
/*! \class  ICM_ImpLossTree icm_implied_loss_tree.h "icm_implied_loss_tree.h"
 *  \author : Damien pouponneau
 *	\version 0.0
 *	\date   30 may 2006
 *	\file   icm_implied_loss_tree.h
 *	\brief  build Index Loss tree 
/**********************************************************************************/


#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\mod\icm_binomialtree.h"


class ImpliedNode
{
public:
	ImpliedNode() : fwd_state(0.),p_def(0.),state(0.),p_nodef(0.) {}
	double p_def;
	double p_nodef;
	double state;
	double fwd_state;
};

class NagSplineHolder; 

class ICM_ImpLossTree : public ARM_Object 
{
	
	ICM_DefaultCurve*			itsIdxDefaultCurve;
	ICM_Smile_Correlation*		itsCorrelation;
	qCAL_INDEX_CORR_TYPE		itsCorreltype;

	double	itsLossUnit;

	ICM_QMatrix<double>			itsTimeStep;
	vector< BinomialTree<ImpliedNode> >	itsTreeVector;
	vector< ICM_QMatrix<double> > itsELoss;

	vector<ICM_QMatrix<double> > itsPricesBefore;
	vector<ICM_QMatrix<double> > itsPricesAfter;
	vector< ICM_QMatrix<double> > itsELossAfter;
	
	// vector<NagSplineHolder*>			itsSpleens;
	NagSplineHolder * itsSpleens; 
	vector<double>				itsStrikes;
	vector<double>				itsSchedCDO;

	int		its_cal_noindex;
	int		its_cal_slice;
	int		its_cal_noloss;

	bool	its_IsCalibrate;
	bool    its_IsAdjusted;

	int itsStep;

private:    

	void Init() ;


public:
	ICM_ImpLossTree() {Init();};

	ICM_ImpLossTree(ICM_DefaultCurve* defcurve,
					ICM_Smile_Correlation* correl,
					qCAL_INDEX_CORR_TYPE correltype,
					bool Adjusted,
					int step = 30) 
	{Init();
	 Set(defcurve,correl,correltype,Adjusted,step);}

	double InterpolEL(double strike, double maturity);

	void Set(ICM_DefaultCurve* defcurve,
					ICM_Smile_Correlation* correl,
					qCAL_INDEX_CORR_TYPE correltype,
					bool Adjusted,
					int step = 30)
	{
		itsTreeVector.clear();
		if (itsIdxDefaultCurve) delete itsIdxDefaultCurve;
		itsIdxDefaultCurve = (ICM_DefaultCurve*)defcurve->Clone();

		if (itsCorrelation) delete itsCorrelation;
		itsCorrelation = (ICM_Smile_Correlation*)correl->Clone();

		itsCorreltype=correltype;
		its_IsAdjusted = Adjusted;
		itsStep = step;

		//calibrate();
	}

	~ICM_ImpLossTree(); 

	void CptImplicitEL(double Maturity = 0.,bool tree = true);			//compute implicit EL
	void CptImplicitELforsplines(double Maturity = 0.);			//compute implicit EL
	
	void CreateTrees();				//create empty trees
	void CptStateLosses(bool adjust = false);			//compute state losses

	void CptTransitionProbaForNode(const int noindex,
								   const int slice,
								   const int noloss,
								   bool adjust);

	void CptTransitionProbaForSlice(const int noindex,
									const int slice,
									bool adjust);

	void CptTransitionProba(bool adjust = false);	//compute transition probabilities

	void CptStateLossesFromProbaForNode(const int noindex,
										const int slice,
										const int noloss,
										const double proba);

	void CptNextStateLossesFromProbaForNode(const int noindex,
											const int slice,
											const int noloss,
											const double proba);

	double CptNextStateLossesFromProbaForNodeRoot(double proba);
	double CptNextStateLossesFromProbaForNodeRootDF(double proba);

	void CalibrateTree();			//calibrate tree with adjusted probabilities
	void CptPrices();				//compute new EL

	void calibrate(double Maturity = 0.,bool tree = true)
	{

		if (tree) 
		{	CptImplicitEL(Maturity,tree);
			BuildTree();}
		else 
		{	CptImplicitELforsplines(Maturity);
			BuildSplines();}

	}

	void BuildTree()
	{
		CreateTrees();
		CptStateLosses(its_IsAdjusted);
		CptTransitionProba(its_IsAdjusted);
		//CalibrateTree();
		its_IsCalibrate=true;
	}

	void BuildSplines();

	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	double CptCumEL(int IndMatu,
						 int dw_IndLoss,
						 int up_IndLoss,
						 int noindex,
						 bool fwdstatus=false,
						 double notional_dw = 0.,
						 double notional_up = 0.);
	
	double CptCumELSpline(double strikedw,
		 			      double strikeup,
						  double maturity);


	bool IsCalibrated() {return its_IsCalibrate;}
	void SetCalibrated(bool status) {its_IsCalibrate=true;}

	void ResetTree() 
	{
		itsLossUnit=-1;
		itsTreeVector.clear();
		itsELoss.clear();

		itsPricesBefore.clear();
		itsPricesAfter.clear();
		itsELossAfter.clear();

		its_cal_noindex=-1;
		its_cal_slice=-1;
		its_cal_noloss=-1;

		its_IsCalibrate=false;
	}

	void SetTimeStep(ICM_QMatrix<double>& times) {itsTimeStep=times;}
	ICM_QMatrix<double> GetTimeStep() {return itsTimeStep;}

	BinomialTree<ImpliedNode>& GetTreeVector(int i) {return itsTreeVector[i];}
	vector<BinomialTree<ImpliedNode> >& GetTreeVector() {return itsTreeVector;}

	void SetValuesCal(int noindex,int slice,int noloss)
	{
	its_cal_noindex=noindex;
	its_cal_slice=slice;
	its_cal_noloss=noloss;
	}

	void ResetValuesForFwdStates(int noindex,
								 int slice,
								 int noloss,
								 double value=1.)
	{
		for (int idxslice=0;idxslice<itsTreeVector[noindex].depth();idxslice++)
			for (int idxnoloss=0;idxnoloss<=idxslice;idxnoloss++)
			{
				itsTreeVector[noindex].data(idxslice,idxnoloss).fwd_state=0.;
				if ((idxslice==slice) && (idxnoloss==noloss)) 
				{itsTreeVector[noindex].data(idxslice,idxnoloss).fwd_state=value;}

			}

	}	
	
	void CptNextStateLossesFromProbaForNode(const int noindex,const int slice,
											const int noloss,bool fwd=false);

	void CptNextStateLossesForSlice(const int noindex,const int slice,bool fwd=false);
	void DiffuseStateLosses(int begin,int end,bool fwd=false);

	ICM_Smile_Correlation*	GetCorrelation() {return itsCorrelation;}

	inline double GetLossUnit() {return itsLossUnit;}
private:
	ICM_ImpLossTree(const ICM_ImpLossTree&); //NA 
	ICM_ImpLossTree& operator=(const ICM_ImpLossTree&); //NA 
};


#endif
