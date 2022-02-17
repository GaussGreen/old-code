#ifndef _ICM_IMPLIED_LOSS_TREE_H
#define _ICM_IMPLIED_LOSS_TREE_H

/*********************************************************************************/
/*! \class  ICM_HWTree_intensity icm_implied_loss_tree.h "icm_implied_loss_tree.h"
 *  \author : Damien pouponneau
 *	\version 0.0
 *	\date   February 2007
 *	\file   icm_hw_credit_tree_dl.h
 *	\brief  build Hull White tree following paper Dynamic Models of portfolio credit risk: Simplified approach 
/**********************************************************************************/


#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\mod\icm_binomialtree.h"

class HW_node
{
public:
	HW_node() : theta(0.),jumpsize(0.),lambda(0.),proba(0.),EP_i_j(0.),FwdSpread(CREDIT_DEFAULT_VALUE) {}
	double theta;
	double jumpsize;
	double lambda;
	double proba;
	double EP_i_j;
	double FwdSpread;
	double duration;
};

class ICM_HWTree_intensity : public ARM_Object 
{
	ARM_Vector				itsTimeStep;
	BinomialTree<HW_node> 	itsHWTree;

private:    
	bool view_all;
	bool view_theta_only;
	bool view_jumpsize_only;
	bool view_lambda_only;
	bool view_proba_only;

	void Init()	
	{
		view_all=true;
		view_theta_only=true;
		view_jumpsize_only=true;
		view_lambda_only=true;
		view_proba_only=true;
	}

public:
	ICM_HWTree_intensity() {Init();};
	ICM_HWTree_intensity(ARM_Vector& time) 
	{
		Init();
		SetTimeStep(time);
		BuildTree();	};


	~ICM_HWTree_intensity(){} 

	void BuildTree()
	{
		ICM_QMatrix<double> disctimes(itsTimeStep.GetSize(),1);

		for (int i=0;i<itsTimeStep.GetSize();i++)
		{disctimes(i,0)=itsTimeStep.Elt(i);}

		itsHWTree.setTimes(disctimes);
	}

	void SetTimeStep(ARM_Vector& times) {itsTimeStep=times;}
	ARM_Vector& GetTimeStep() {return itsTimeStep;}

	BinomialTree<HW_node>& GetTree() {return itsHWTree;}

	void DiffuseTree(bool corp = true);
	void DiffuseTree(int i0,int j0,bool corp);
	double CorpEsp(int& i,bool corp = true);
	void SetTheta(int& period,double value);
	double SumProba(int& i);
	double SumJumps(int& i);

	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	double EC_Node_until_maturity(int i,int j,vector<double>& Epdef);

};


#endif
