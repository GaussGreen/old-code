
#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"

void ICM_Gauss1FLossDistrib_LJ::Set_VN(CtxtDistrib* c)
{
	SetCtxt(c);

	int size = c->m_IntegrationStep;

	int lup=floor(c->m_TrancheUp/c->m_LossUnit);
	int ldown=floor(c->m_TrancheDown/c->m_LossUnit);

	if (c->m_ProbaCond_Down) delete c->m_ProbaCond_Down;
	c->m_ProbaCond_Down= new ICM_QCubix<double>(c->m_nbnames+1,size,ldown+1,-999.);

	if (c->m_ProbaCond_Up) delete c->m_ProbaCond_Up;
	c->m_ProbaCond_Up= new ICM_QCubix<double>(c->m_nbnames+1,size,lup+1,-999.);

	vector<CtxtDistrib>* vdistrib = &(c->m_vn_contexts);

	for (int t=0;t<vdistrib->size();t++)
	{
		CtxtDistrib* c_ = NULL;

		c_ = &((*vdistrib)[t]);
		if (c_->m_ProbaCond_Down) delete c_->m_ProbaCond_Down;
		c_->m_ProbaCond_Down= new ICM_QCubix<double>(c->m_nbnames+1,size,ldown+1,-999.);

		if (c_->m_ProbaCond_Up) delete c_->m_ProbaCond_Up;
		c_->m_ProbaCond_Up= new ICM_QCubix<double>(c->m_nbnames+1,size,lup+1,-999.);
	}

	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ((c->m_IntegrationStep % 2 == 0) && (c->m_IntegrationMethod == qGAUSS_LEGENDRE))
		c->m_IntegrationMethod = qGAUSS_HERMITE;
	else if ((c->m_IntegrationStep % 2 == 1) && (c->m_IntegrationMethod == qGAUSS_HERMITE))
		c->m_IntegrationMethod = qGAUSS_LEGENDRE;
	
	c->ComputeAll();
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_VN(const int& lup)
{
	switch (GetCtxt()->m_IntegrationMethod)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	case qTRAPEZE:
	default:
		compute_distrib_Gauss_Legendre_VN(lup);
	break;
	}
}

// -------------------------------------------------------------------
// GAUSSIAN COPULA
// -------------------------------------------------------------------

//	Loss Level is 0.0

extern "C" void 
compute_cond_distribZeroLoss_Gaussian_Integrator_VN(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distribZeroLoss_Gaussian_VN(x,Param);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_Gaussian_VN(double x, void* params)
{
	int t;
	double	pk=0.,pk_suiv = 0.,pk_suiv_cond = 0.;
	double	tmp_barrier=0.;

	its_ind_x+=1;

	unsigned long sizemax=-1;
	(*(AddressVector*)params).GetCount(sizemax);

	ICM_QCubix<double>* PROBASCUM = (ICM_QCubix<double>*)((*(AddressVector*)params)[sizemax-4]);
	ICM_QCubix<double>* PROBAS = (ICM_QCubix<double>*)((*(AddressVector*)params)[sizemax-3]);
	ICM_QMatrix<double>* matrix = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)params)[sizemax-1]);

//	double CondProbaZeroLoss=1.;

	CtxtDistrib* c0 = GetCtxt();
	CtxtDistrib* c = NULL;
	CtxtDistrib* c_suiv = NULL;
	vector<CtxtDistrib>* vdistrib = &(c0->m_vn_contexts);

//	ICM_QCubix<double>* pProbCond_Up = c0->m_ProbaCond_Up;
//	pProbCond_Up->SetElt(0,its_ind_x,0,CondProbaZeroLoss);

	ICM_Gauss1FLossDistrib_LJ*	TheModel =	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)params)[0]);
	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod ==	qGAUSS_HERMITE)
		x	*=	SQRT2;


	for (t=0;t<vdistrib->size()-1;t++)
	{
			c_suiv = &((*vdistrib)[t+1]);
			c = &((*vdistrib)[t]);		

			ICM_QCubix<double>* pProbCond_Up_ = c->m_ProbaCond_Up;
			pProbCond_Up_->SetElt(0,its_ind_x,0,1. /*CondProbaZeroLoss*/ );

			double CondProbaZeroLoss_ = 1.;
//			pk_suiv_cond = 1.;

			for (int k=1;k<=c0->m_nbnames;k++) 
			{

				tmp_barrier	=(c_suiv->m_barriers_maturity[k-1]-c_suiv->m_beta_up_maturity[k-1]*x)/sqrt(1.-c_suiv->m_beta_up_maturity[k-1]*c_suiv->m_beta_up_maturity[k-1]);
				pk_suiv=NAG_cumul_normal(tmp_barrier);

				tmp_barrier=(c->m_barriers_maturity[k-1]-c->m_beta_up_maturity[k-1]*x)/sqrt(1.-c->m_beta_up_maturity[k-1]*c->m_beta_up_maturity[k-1]);
				pk=NAG_cumul_normal(tmp_barrier);

				pk_suiv_cond = (pk_suiv - pk);

				double LossRates_inf = floor(c_suiv->m_LossRates[k-1]);
				if (!CHECK_NULL(LossRates_inf))
				{CondProbaZeroLoss_ *= (1.-pk_suiv_cond);}

				pProbCond_Up_->SetElt(k,its_ind_x,0,CondProbaZeroLoss_);
			}

			(*PROBAS).SetElt(its_ind_x,t+1,0,CondProbaZeroLoss_);

	}

	double result = 1.; 

	for (t=1;t<vdistrib->size();t++)
	{
		result *= (*PROBAS).Elt(its_ind_x,t,0);
		(*PROBASCUM).SetElt(its_ind_x,t,0,result);
	}

	return (result);
}


//	Loss Level is 0.0

extern "C" void 
compute_cond_distribZeroLoss_Gaussian_Smile_VN(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	int k;
	double	pk=0.,pk_down=0.;
	double	tmp_barrier=0.;

	double CondProbaZeroLoss	=	1.;
	double CondProbaZeroLoss_Down=	1.;

	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	unsigned long sizemax=-1;
	(*(AddressVector*)Param).GetCount(sizemax);

	ICM_QCubix<double>* PROBASCUM_UP = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-7]);
	ICM_QCubix<double>* PROBASCUM_DW = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-6]);
	ICM_QCubix<double>* PROBAS_UP = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-5]);
	ICM_QCubix<double>* PROBAS_DW = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-4]);

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	CtxtDistrib* c0 = TheModel->GetCtxt();
	CtxtDistrib* c = NULL;
	CtxtDistrib* c_suiv = NULL;
	vector<CtxtDistrib>* vdistrib = &(c0->m_vn_contexts);

	ICM_QCubix<double>* pProbCond_Up = c0->m_ProbaCond_Up;
	ICM_QCubix<double>* pProbCond_Down = c0->m_ProbaCond_Down;

	pProbCond_Up->SetElt(0, *its_ind_x,0,CondProbaZeroLoss);
	pProbCond_Down->SetElt(0, *its_ind_x,0,CondProbaZeroLoss_Down);

	double	pk_suiv=0.,pk_down_suiv=0.;

	for (int t=0;t<vdistrib->size()-1;t++)
	{
			c_suiv = &((*vdistrib)[t+1]);
			c = &((*vdistrib)[t]);

			ICM_QCubix<double>* pProbCond_Up_ = c->m_ProbaCond_Up;
			pProbCond_Up_->SetElt(0,*its_ind_x,0,CondProbaZeroLoss);

			ICM_QCubix<double>* pProbCond_Down_ = c->m_ProbaCond_Down;
			pProbCond_Down_->SetElt(0,*its_ind_x,0,CondProbaZeroLoss);

			double CondProbaZeroLoss_Up = 1.,CondProbaZeroLoss_Down = 1.;
			double pk_suiv_cond = 1.,pk_down_suiv_cond = 1.;

			for (k=0;k<c0->m_nbnames;k++) 
			{

				tmp_barrier	=	c_suiv->m_a_up_maturity[k] - c_suiv->m_b_up_maturity[k] * x;
				pk_suiv=NAG_cumul_normal(tmp_barrier);

				tmp_barrier	=	c_suiv->m_a_dw_maturity[k] - c_suiv->m_b_dw_maturity[k] * x;
				pk_down_suiv=NAG_cumul_normal(tmp_barrier);

				tmp_barrier	=	c->m_a_up_maturity[k] - c->m_b_up_maturity[k] * x;
				pk=NAG_cumul_normal(tmp_barrier);
				pk_suiv_cond = (pk_suiv - pk);

				tmp_barrier	=	c->m_a_dw_maturity[k] - c->m_b_dw_maturity[k] * x;
				pk_down=NAG_cumul_normal(tmp_barrier);
				pk_down_suiv_cond = (pk_down_suiv - pk_down);

				double LossRates_inf = floor(c_suiv->m_LossRates[k]);
				if (!CHECK_NULL(LossRates_inf))
				{
				CondProbaZeroLoss_Up *= (1.-pk_suiv_cond);
				CondProbaZeroLoss_Down *= (1.-pk_down_suiv_cond);
				}
			
				pProbCond_Up_->SetElt(k+1, *its_ind_x,0,CondProbaZeroLoss_Up);
				pProbCond_Down_->SetElt(k+1, *its_ind_x,0,CondProbaZeroLoss_Down);
			}

			(*PROBAS_UP).SetElt(*its_ind_x,t+1,0,CondProbaZeroLoss_Up);
			(*PROBAS_DW).SetElt(*its_ind_x,t+1,0,CondProbaZeroLoss_Down);

	}

	double result_up = 1.,result_dw = 1.; 

	for (t=1;t<vdistrib->size();t++)
	{
		result_up *= (*PROBAS_UP).Elt(*its_ind_x,t,0);
		(*PROBASCUM_UP).SetElt(*its_ind_x,t,0,result_up);

		result_dw *= (*PROBAS_DW).Elt(*its_ind_x,t,0);
		(*PROBASCUM_DW).SetElt(*its_ind_x,t,0,result_dw);
	}

	ProbCond_Up		=	result_up; /* pProbCond_Up->Elt(c0->m_nbnames,*its_ind_x,0);	*/
	ProbCond_Down	=	result_dw; /* pProbCond_Down->Elt(c0->m_nbnames, *its_ind_x,0);	*/
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------

extern "C" void 
compute_cond_distrib_Gaussian_Integrator_VN(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distrib_Gaussian_VN(x,Param);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_Gaussian_VN(double x,void* params)
{

	its_ind_x+=1;

	int t;
	int	k;
	int j;
	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp1_=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0. /*,tmp3=0.*/; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	double	pk_suiv=0.;
	double	pk_suiv_cond=0.;
	int ind_lastloss=0;

	double tmp_barrier=0./*,tmp4=0.*/;

	unsigned long sizemax=-1;
	(*(AddressVector*)params).GetCount(sizemax);

	ICM_QCubix<double>* PROBASCUM = (ICM_QCubix<double>*)((*(AddressVector*)params)[sizemax-4]);
	ICM_QCubix<double>* PROBAS = (ICM_QCubix<double>*)((*(AddressVector*)params)[sizemax-3]);
	ICM_QMatrix<double>* matrix = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)params)[sizemax-1]);

	CtxtDistrib* c0 = GetCtxt();
	CtxtDistrib* c = NULL;
	CtxtDistrib* c_suiv = NULL;
	vector<CtxtDistrib>* vdistrib = &(c0->m_vn_contexts);

	// HERMITE
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	its_sqrt2;

		for (t=0;t<vdistrib->size()-1;t++)
		{
				c_suiv = &((*vdistrib)[t+1]);
				c = &((*vdistrib)[t]);

				ICM_QCubix<double>* pProbCond_Up_ = c->m_ProbaCond_Up;
				pProbCond_Up_ ->SetElt(0,its_ind_x,its_ind_Loss,0.);

				double prob = (*PROBAS).Elt(its_ind_x,t+1,its_ind_Loss);
				
				if (CHECK_EQUAL(prob, (CREDIT_DEFAULT_VALUE)))
				{
				for (k=1;k<=c0->m_nbnames;k++) 
				{
	
					tmp_barrier=(c_suiv->m_barriers_maturity[k-1]-c_suiv->m_beta_up_maturity[k-1]*x)/sqrt(1.-c_suiv->m_beta_up_maturity[k-1]*c_suiv->m_beta_up_maturity[k-1]);
					pk_suiv=NAG_cumul_normal(tmp_barrier);

					// current level of loss - loss rate of names k
					double LossRates_inf = floor(c_suiv->m_LossRates[k-1]);
					double LossRates_delta = c_suiv->m_LossRates[k-1] - LossRates_inf;

					ind_lastloss = its_ind_Loss - (int)LossRates_inf;

					if CHECK_NULL(LossRates_inf)
					{
						pProbCond_Up_->SetElt(k,its_ind_x, its_ind_Loss,pProbCond_Up_->Elt(k-1,its_ind_x, its_ind_Loss));
						continue;
					}

					// barrier
					tmp_barrier=(c->m_barriers_maturity[k-1]-c->m_beta_up_maturity[k-1]*x)/sqrt(1.-c->m_beta_up_maturity[k-1]*c->m_beta_up_maturity[k-1]);
					// conditional default probability
					pk=NAG_cumul_normal(tmp_barrier);

					if (!(fabs(1.-pk)<1.e-15))
					{pk_suiv_cond = (pk_suiv - pk)/(1.-pk);}
					else
					{double test=0.;}

//					if (c->m_IsBegin)
					{pk_suiv_cond = (pk_suiv - pk);}


					tmp1 = pProbCond_Up_->Elt(k-1, its_ind_x, its_ind_Loss);

					// Loss Distribution for a basket of size k knowing factor value its_ind_x
					tmp2	=	tmp1 * (1.-pk_suiv_cond);
			
					if (ind_lastloss >= 0)
					{
						tmp1_ = pProbCond_Up_->Elt(k-1,its_ind_x,ind_lastloss); 
						tmp2	+=	pk_suiv_cond * tmp1_ * (1.-LossRates_delta);

						if ((ind_lastloss>0)&&(LossRates_delta))
						{	
						tmp1_ = pProbCond_Up_->Elt(k-1,its_ind_x,ind_lastloss-1); 
						tmp2	+=	pk_suiv_cond * tmp1_ * LossRates_delta;
						}
					}

					// else: the loss of name k is too big
					pProbCond_Up_->SetElt(k,its_ind_x,its_ind_Loss, tmp2);

				}

				(*PROBAS).SetElt(its_ind_x,t+1,its_ind_Loss,pProbCond_Up_->Elt(c0->m_nbnames,its_ind_x,its_ind_Loss));
				}
		}

	if CHECK_EQUAL((*PROBASCUM).Elt(its_ind_x,1,its_ind_Loss),(CREDIT_DEFAULT_VALUE))
	{
			for (k=0;k<=its_ind_Loss;k++)
			{(*PROBASCUM).Elt(its_ind_x,1,k)=(*PROBAS).Elt(its_ind_x,1,k);}
	}

	for (t=1;t<vdistrib->size()-1;t++)
	{

		if (CHECK_EQUAL((*PROBASCUM).Elt(its_ind_x,t+1,its_ind_Loss),(CREDIT_DEFAULT_VALUE)))
			{
				double intermed=0.;
				for (j=0;j<=its_ind_Loss;j++)
				{intermed += (*PROBASCUM).Elt(its_ind_x,t,j)*(*PROBAS).Elt(its_ind_x,t+1,its_ind_Loss-j);}

				(*PROBASCUM).Elt(its_ind_x,t+1,its_ind_Loss)=intermed;
			}
	}

	double result = (*PROBASCUM).Elt(its_ind_x,vdistrib->size()-1,its_ind_Loss);

	return (result);
}


extern "C" void 
compute_cond_distrib_Gaussian_Smile_VN(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	unsigned long sizemax=-1;
	(*(AddressVector*)Param).GetCount(sizemax);

	ICM_QCubix<double>* PROBASCUM_UP = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-7]);
	ICM_QCubix<double>* PROBASCUM_DW = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-6]);
	ICM_QCubix<double>* PROBAS_UP = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-5]);
	ICM_QCubix<double>* PROBAS_DW = (ICM_QCubix<double>*)((*(AddressVector*)Param)[sizemax-4]);
	ICM_QMatrix<double>* matrix_up = (ICM_QMatrix<double>*)((*(AddressVector*)Param)[sizemax-3]);
	ICM_QMatrix<double>* matrix_down = (ICM_QMatrix<double>*)((*(AddressVector*)Param)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)Param)[sizemax-1]);

	int	k,j;
	double tmp1,tmp1_; // proba de loss=L pour un panier de taille k
	double tmp2,tmp2_down; // proba de loss=L-1 pour un panier de taille k
	double pk=0., pk_down=0.;     // proba cond de défaut du nom k
	double	pk_suiv=0.,pk_down_suiv=0.;
	int ind_lastloss;

	double tmp_barrier, tmp_barrier_down;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	CtxtDistrib* c0 = TheModel->GetCtxt();
	CtxtDistrib* c = NULL;
	CtxtDistrib* c_suiv = NULL;
	vector<CtxtDistrib>* vdistrib = &(c0->m_vn_contexts);

//	ICM_QCubix<double>* pProbCond_Up = c0->m_ProbaCond_Up;
//	ICM_QCubix<double>* pProbCond_Down = c0->m_ProbaCond_Down;

//	pProbCond_Up->SetElt(0,*its_ind_x,*its_ind_Loss,0.);
//	pProbCond_Down->SetElt(0,*its_ind_x,*its_ind_Loss,0.);


	for (int t=0;t<vdistrib->size()-1;t++)
	{
			// next probabilities
			c_suiv = &((*vdistrib)[t+1]);
			// previous probabilities
			c = &((*vdistrib)[t]);	

			ICM_QCubix<double>* pProbCond_Up_ = c->m_ProbaCond_Up;
			pProbCond_Up_->SetElt(0,*its_ind_x,*its_ind_Loss,0.);

			ICM_QCubix<double>* pProbCond_Down_ = c->m_ProbaCond_Down;
			pProbCond_Down_->SetElt(0,*its_ind_x,*its_ind_Loss,0.);

			double pk_suiv_cond = 1.,pk_down_suiv_cond = 1;
				
			for (k=0;k<c0->m_nbnames;k++) 
			{

				tmp_barrier	=	c_suiv->m_a_up_maturity[k] - c_suiv->m_b_up_maturity[k] * x;
				pk_suiv=NAG_cumul_normal(tmp_barrier);

				tmp_barrier	=	c_suiv->m_a_dw_maturity[k] - c_suiv->m_b_dw_maturity[k] * x;
				pk_down_suiv=NAG_cumul_normal(tmp_barrier);

				double LossRates_inf = floor(c_suiv->m_LossRates[k]);
				double LossRates_delta = c_suiv->m_LossRates[k] - LossRates_inf;

				ind_lastloss = *its_ind_Loss - (int)LossRates_inf;

				if CHECK_NULL(LossRates_inf)
				{
					pProbCond_Up_->SetElt(k+1,*its_ind_x, *its_ind_Loss,pProbCond_Up_->Elt(k,*its_ind_x, *its_ind_Loss));
					pProbCond_Down_->SetElt(k+1,*its_ind_x, *its_ind_Loss,pProbCond_Down_->Elt(k,*its_ind_x, *its_ind_Loss));
					continue;
				}

				// barrier
				tmp_barrier = c->m_a_up_maturity[k] - c->m_b_up_maturity[k] * x;
				tmp_barrier_down = c->m_a_dw_maturity[k] - c->m_b_dw_maturity[k] * x;
		
				// conditional default probability
				pk=NAG_cumul_normal(tmp_barrier);
				pk_down=NAG_cumul_normal(tmp_barrier_down);

				pk_suiv_cond = 1.;
				pk_down_suiv_cond = 1.;

//				if (!(fabs(1.-pk)<1.e-15))
//				{pk_suiv_cond = (pk_suiv - pk)/(1.-pk);}

				pk_suiv_cond = (pk_suiv - pk);

//				if (!(fabs(1.-pk_down)<1.e-15))
//				{pk_down_suiv_cond = (pk_down_suiv - pk_down)/(1.-pk_down);}

				pk_down_suiv_cond = (pk_down_suiv - pk_down);

				{
					//strike up
					tmp2 = pProbCond_Up_->Elt(k, *its_ind_x, *its_ind_Loss);
					tmp1 = tmp2 * (1. - pk_suiv_cond);

					if (ind_lastloss >= 0)
					{	
						tmp1_ = pProbCond_Up_->Elt(k,*its_ind_x,ind_lastloss); 
						tmp1	+=	pk_suiv_cond * tmp1_ * (1.-LossRates_delta);
						
						if ((ind_lastloss > 0)&&(LossRates_delta))
						{	
						tmp1_ = pProbCond_Up_->Elt(k,*its_ind_x,ind_lastloss-1); 
						tmp1	+=	pk_suiv_cond * tmp1_ * LossRates_delta;
						}
					}

					pProbCond_Up_->SetElt(k+1,*its_ind_x, *its_ind_Loss,tmp1);
					//tmp2 = tmp1;

					//strike down
					tmp2_down = pProbCond_Down_->Elt(k, *its_ind_x, *its_ind_Loss);
					tmp1 = tmp2_down * (1. - pk_down_suiv_cond);

					if (ind_lastloss >= 0)
					{	
						tmp1_ = pProbCond_Down_->Elt(k,*its_ind_x,ind_lastloss); 
						tmp1	+=	pk_down_suiv_cond * tmp1_ * (1.-LossRates_delta);	

						if ((ind_lastloss > 0)&&(LossRates_delta))
						{	
						tmp1_ = pProbCond_Down_->Elt(k,*its_ind_x,ind_lastloss-1); 
						tmp1	+=	pk_down_suiv_cond * tmp1_ * LossRates_delta;
						}
					}

					pProbCond_Down_->SetElt(k+1,*its_ind_x, *its_ind_Loss,tmp1);
					//tmp2_down = tmp1;
				}

			}

			(*PROBAS_UP).SetElt(*its_ind_x,t+1,*its_ind_Loss,pProbCond_Up_->Elt(c0->m_nbnames,*its_ind_x,*its_ind_Loss));
			(*PROBAS_DW).SetElt(*its_ind_x,t+1,*its_ind_Loss,pProbCond_Down_->Elt(c0->m_nbnames,*its_ind_x,*its_ind_Loss));

	}

	if CHECK_EQUAL((*PROBASCUM_UP).Elt(*its_ind_x,1,*its_ind_Loss),(CREDIT_DEFAULT_VALUE))
	{
			for (k=0;k<=*its_ind_Loss;k++)
			{(*PROBASCUM_UP).Elt(*its_ind_x,1,k)=(*PROBAS_UP).Elt(*its_ind_x,1,k);
			(*PROBASCUM_DW).Elt(*its_ind_x,1,k)=(*PROBAS_DW).Elt(*its_ind_x,1,k);
			}
	}

	for (t=1;t<vdistrib->size()-1;t++)
	{
	//	for (k=0;k<=*its_ind_Loss;k++)
	//	{
	//		if (!CHECK_EQUAL((*PROBASCUM_UP).Elt(*its_ind_x,t+1,k),(CREDIT_DEFAULT_VALUE)))
	//		{continue;}

			if (CHECK_EQUAL((*PROBASCUM_UP).Elt(*its_ind_x,t+1,*its_ind_Loss),(CREDIT_DEFAULT_VALUE)))
			{
			double intermed_up=0.,intermed_dw=0.;
			for (j=0;j<=*its_ind_Loss;j++)
			{
				intermed_up += (*PROBASCUM_UP).Elt(*its_ind_x,t,j)*(*PROBAS_UP).Elt(*its_ind_x,t+1,*its_ind_Loss-j);
				intermed_dw += (*PROBASCUM_DW).Elt(*its_ind_x,t,j)*(*PROBAS_DW).Elt(*its_ind_x,t+1,*its_ind_Loss-j);
			}

			(*PROBASCUM_UP).Elt(*its_ind_x,t+1,*its_ind_Loss)=intermed_up;
			(*PROBASCUM_DW).Elt(*its_ind_x,t+1,*its_ind_Loss)=intermed_dw;
			}
	//	}
	}

	//double result = (*PROBASCUM).Elt(*its_ind_x,vdistrib->size()-1,*its_ind_Loss);


	ProbCond_Up		=	 (*PROBASCUM_UP).Elt(*its_ind_x,vdistrib->size()-1,*its_ind_Loss);
	ProbCond_Down	=	 (*PROBASCUM_DW).Elt(*its_ind_x,vdistrib->size()-1,*its_ind_Loss);
}

// -------------------------------------------------------------------
// STUDENT COPULA
// -------------------------------------------------------------------
//	Loss Level is 0.0

double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_Student_VN(const double& x)
{
	double	pk=0.,pk_suiv = 0.,pk_suiv_cond = 0.;
	double	tmp_barrier=0.;

	its_ind_x+=1;

	double CondProbaZeroLoss=1.;

	CtxtDistrib* c_suiv = NULL;
	CtxtDistrib* c = GetCtxt();
	ICM_QCubix<double>* pProbCond_Up = c->m_ProbaCond_Up;

	vector<CtxtDistrib>* vdistrib = &(c->m_vn_contexts);

	pProbCond_Up->SetElt(0,its_ind_x,0,CondProbaZeroLoss);

	for (int k=1;k<=c->m_nbnames;k++) 
	{
		for (int t=0;t<vdistrib->size()-1;t++)
		{
			c_suiv = &((*vdistrib)[t+1]);
			c = &((*vdistrib)[t]);

			tmp_barrier=(c->m_barriers_maturity[k-1]-c->m_beta_up_maturity[k-1]*x)/sqrt(1.-c->m_beta_up_maturity[k-1]*c->m_beta_up_maturity[k-1]);
			pk = NAG_prob_students_t( tmp_barrier, itsFreedomDegree );

			tmp_barrier=(c_suiv->m_barriers_maturity[k-1]-c_suiv->m_beta_up_maturity[k-1]*x)/sqrt(1.-c_suiv->m_beta_up_maturity[k-1]*c_suiv->m_beta_up_maturity[k-1]);
			pk_suiv = NAG_prob_students_t( tmp_barrier, itsFreedomDegree );

			pk_suiv_cond = 1.;
			if (!(fabs(1.-pk)<1.e-15))
			{pk_suiv_cond = (pk_suiv - pk)/(1.-pk);}

			if (c->m_IsBegin)
			{
				CondProbaZeroLoss = pProbCond_Up->Elt(k-1,its_ind_x,0);
				CondProbaZeroLoss *= (1.-pk_suiv);
			}
			else
			{
				CondProbaZeroLoss *= (1.-pk_suiv_cond);
			}

			pProbCond_Up->SetElt(k,its_ind_x,0,CondProbaZeroLoss);
		}
	}

	return pProbCond_Up->Elt(c->m_nbnames,its_ind_x,0);
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------

double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_Student_VN(const double& x)
{

	//TO DO : modification non effectuée 
	its_ind_x+=1;

	int	k;
	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0.,tmp3=0.; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	int wk=0,ind_lastloss=0;

	double tmp_barrier=0.,tmp4=0.;

	CtxtDistrib* c = GetCtxt();
	ICM_QCubix<double>* pProbCond_Up = c->m_ProbaCond_Up;
	ICM_QCubix<double>* pProbCond_Down = c->m_ProbaCond_Down;

	pProbCond_Up->SetElt(0,its_ind_x,its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1	=	its_ProbCond->Elt(0, its_ind_x, its_ind_Loss);


	for (k=1;k<=c->m_nbnames;k++) 
	{
			// current level of loss - loss rate of names k
			ind_lastloss = its_ind_Loss - c->m_LossRates[k-1];

			// barrier
			tmp_barrier=(c->m_barriers_maturity[k-1]-c->m_beta_up_maturity[k-1]*x)/sqrt(1.-c->m_beta_up_maturity[k-1]*c->m_beta_up_maturity[k-1]);

			// conditional default probability
			// Cumulative Student distribution
			pk = NAG_prob_students_t( tmp_barrier, itsFreedomDegree);

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2	=	tmp1 * (1.-pk);
			
			if (ind_lastloss >= 0)
				tmp2	+=	pk * pProbCond_Up->Elt(k-1,its_ind_x,ind_lastloss);
			// else: the loss of name k is too big

			pProbCond_Up->SetElt(k,its_ind_x,its_ind_Loss, tmp2);

			tmp1 = tmp2;
	}


	return tmp1;
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_VN(const int& lup,const int& ldown)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	// ---------------------
	its_lup	=	0;
//	max_lup = MAX(lup,its_lup);
	
	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;
	AddressVector	TheParameterVector2;

	TheParameterVector.Append(this);
	TheParameterVector2.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector2.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);
	TheParameterVector2.Append(&its_ind_Loss);

	CtxtDistrib* c = GetCtxt();

	vector<CtxtDistrib>* vdistrib = &(c->m_vn_contexts);

	ICM_QCubix<double> PROBASCUM_UP(c->m_IntegrationStep,vdistrib->size(),lup+1,-999.);
	TheParameterVector.Append(&PROBASCUM_UP);
	TheParameterVector2.Append(&PROBASCUM_UP);

	ICM_QCubix<double> PROBASCUM_DW(c->m_IntegrationStep,vdistrib->size(),lup+1,-999.);
	TheParameterVector.Append(&PROBASCUM_DW);

	ICM_QCubix<double> PROBAS_UP(c->m_IntegrationStep,vdistrib->size(),lup+1,-999.);
	TheParameterVector.Append(&PROBAS_UP);
	TheParameterVector2.Append(&PROBAS_UP);

	ICM_QCubix<double> PROBAS_DW(c->m_IntegrationStep,vdistrib->size(),lup+1,-999.);
	TheParameterVector.Append(&PROBAS_DW);

	ICM_QMatrix<double> COEFS_up(c->m_nbnames,c->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_up);
	TheParameterVector2.Append(&COEFS_up);

	ICM_QMatrix<double> COEFS_down(c->m_nbnames,c->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_down);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);
	TheParameterVector2.Append(&nointegrationcoef);

	c->m_ProbasLossesUp.resize(lup+1);
	c->m_ProbasLossesDown.resize(lup+1);

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (l=0/*its_lup*/;l<=ldown;l++) 
			{
				if (l ==0)
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_Gaussian_Smile_VN, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					c->m_ProbasLossesUp[0] = tmpPrbLoss;
					c->m_ProbasLossesDown[0] = tmpPrbLoss_down;
					its_ind_Loss=1;
				}
				else
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Smile_VN, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					c->m_ProbasLossesUp[l] = tmpPrbLoss;
					c->m_ProbasLossesDown[l] = tmpPrbLoss_down;

					its_ind_Loss+=1;
				}
			}

			COEFS_down.RAZ(-999.);

			for (l=ldown+1;l<=lup;l++) 
			{
				its_ind_x=-1;

				TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Integrator_VN, &TheParameterVector2, tmpPrbLoss); 

				c->m_ProbasLossesUp[l] = tmpPrbLoss;
				its_ind_Loss+=1;
			}
			break;

		case qSTUDENT:

			break;
	}

	its_ind_Loss--;

	for (l=0;l<=ldown;l++)
	{
		cumul_distrib += c->m_ProbasLossesUp[l];
		cumul_distrib_down += c->m_ProbasLossesDown[l];
	}
	for (l=ldown+1;l<=lup;l++)
	{
		cumul_distrib += c->m_ProbasLossesUp[l];
	}

	c->m_TailProbaLosseUp	=	1.0 - cumul_distrib;
	c->m_TailProbaLosseDown	=	1.0 - cumul_distrib_down;

}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_Gauss_Legendre_VN(const int& lup)
{ 
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	// ------------------------
	its_lup		=	0;
	max_lup		=	lup;
//	max_lup = MAX(lup,its_lup);
	// ------------------------

	// ------------------------
	// Parameters
	AddressVector	TheParameterVector;

	CtxtDistrib* c = GetCtxt();

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);
	// ------------------------

	vector<CtxtDistrib>* vdistrib = &(c->m_vn_contexts);

	ICM_QCubix<double> PROBASCUM(c->m_IntegrationStep,vdistrib->size(),max_lup+1,-999.);
	TheParameterVector.Append(&PROBASCUM);

	ICM_QCubix<double> PROBAS(c->m_IntegrationStep,vdistrib->size(),max_lup+1,-999.);
	TheParameterVector.Append(&PROBAS);

	ICM_QMatrix<double> COEFS(c->m_nbnames,c->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	c->m_ProbasLossesUp.resize(max_lup+1);
	c->m_ProbasLossesDown.resize(max_lup+1);

	switch (c->m_CopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (l=0/*its_lup*/;l<=max_lup;l++) 
			{
				if (l ==0)
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_Gaussian_Integrator_VN, &TheParameterVector, tmpPrbLoss); 

					c->m_ProbasLossesUp[0] = tmpPrbLoss;
					its_ind_Loss=1;

				}
				else
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Integrator_VN, &TheParameterVector, tmpPrbLoss); 
					c->m_ProbasLossesUp[l] = tmpPrbLoss;

					its_ind_Loss+=1;
				}
			}

			break;

		case qSTUDENT:
			break;

	}

	its_ind_Loss--;

	for (l=0;l<=lup;l++)
		cumul_distrib += c->m_ProbasLossesUp[l];

	c->m_TailProbaLosseUp	=	1.0 - cumul_distrib;

}


double
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_VN(CtxtDistrib* c)
{
	double exp_loss_tranche_down = 0.;
	double exp_loss_tranche_up = 0.;
	
	int lup=floor(c->m_TrancheUp/c->m_LossUnit);
	int ldown=floor(c->m_TrancheDown/c->m_LossUnit);

	int i=0;	
	
	c->m_ProbasLossesUp.resize(lup+1);
	c->m_ProbasLossesDown.resize(ldown+1);
	// Compute the Distribution (by default, ldown = 0)
	if (c->m_TrancheDown)
		compute_distrib_VN(lup, ldown);
	else
		compute_distrib_VN(lup);

	int l	=	0;
	double	loss_level;

	for (l=1; l<=ldown; l++) 
	{
		loss_level	=	l*c->m_LossUnit;
		exp_loss_tranche_down	+=	c->m_ProbasLossesDown[l] * loss_level;		
		exp_loss_tranche_up		+=	c->m_ProbasLossesUp[l] * loss_level;
	}
		
	// second loop, deal only with beat_up
	for (l=ldown+1;l<=lup;l++)
	{
		loss_level	=	l*c->m_LossUnit;
		exp_loss_tranche_up		+=	c->m_ProbasLossesUp[l] * loss_level;
	}

	c->m_OutLosses = c->m_ProbasLossesUp;

	// add the tail
	if (c->m_TrancheDown)
		exp_loss_tranche_down	+=	c->m_TrancheDown * c->m_TailProbaLosseDown;
	exp_loss_tranche_up		+=	c->m_TrancheUp * c->m_TailProbaLosseUp;

	return (exp_loss_tranche_up - exp_loss_tranche_down);
}


void ICM_Gauss1FLossDistrib_LJ :: SetIntegratorType_VN(const qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep)
{
	qIntegratorChoice itsIntegrationMethod	=	TheIntegratorType;
	int itsIntegrationStep		=	TheStep;

	TheIntegrator.SetIntegrationType(TheIntegratorType);
	
	// Rajout du || sur la méthode qTrapeze afin que le step de l'integrator soit setter a la bonne valeur
	if ((TheIntegratorType == qGAUSS_LEGENDRE) || (TheIntegratorType == qGAUSS_HERMITE) || (TheIntegratorType == qTRAPEZE))
		TheIntegrator.SetIntegrationStep(TheStep);
}



