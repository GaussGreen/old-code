// AbstractCopulaCalculator.cpp: implementation of the AbstractCopulaCalculator class.
//
//////////////////////////////////////////////////////////////////////

#include "ARMKernel\glob\firsttoinc.h" 

#include "ICMKernel\mod\icm_Homogene_FastLossDistrib.h"
#include "ICMKernel\glob\icm_maths.h"

ICM_QMatrix<double> ICM_Gauss1FLossDistrib::itsCombinations=ICM_QMatrix<double>();
ICM_QMatrix<double> ICM_LightDistrib::itsCombinations=ICM_QMatrix<double>();


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void
ICM_Gauss1FLossDistrib::Init()
{
its_ProbCond_indep = NULL;
its_proportion_full_correl = 0.;
its_proportion_independant = 0.;
its_FH_Barrier=0.;
its_FH_i =0;
its_FH_Beta=0.;
its_ProbCond_indep = NULL;
}



//----------------------------------------------------------------------------------------------------------------------
ICM_Gauss1FLossDistrib::ICM_Gauss1FLossDistrib(const int& nbnames,
											   const ARM_Vector&  pdef,
											   const ARM_Vector&  beta,
											   const ARM_Vector&  LossRates)
{
	Init();
	Set(nbnames,pdef,beta,LossRates);
}

void ICM_Gauss1FLossDistrib::Set(const int& nbnames,
											   const ARM_Vector&  pdef,
											   const ARM_Vector&  beta,
											   const ARM_Vector&  LossRates)
{
	Set(nbnames);

	SetUniqueBeta(beta);
	SetPdefAtMaturity(pdef);
	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);
	compute_barrier();		// Copula dependant

	if (its_ProbCond_indep) delete its_ProbCond_indep;
	its_ProbCond_indep= new ICM_QMatrix<double>(nbnames+1,nbnames+1);
}

ICM_Gauss1FLossDistrib::~ICM_Gauss1FLossDistrib()
{

	if(its_ProbCond_indep) 
		delete its_ProbCond_indep; 
	its_ProbCond_indep = NULL;

}

double ICM_Gauss1FLossDistrib::compute_cond_distribZeroLoss(const double& x)
{
	double pk=0.,tmp_barrier=0.;
	its_ind_x+=1;

	double CondProbaZeroLoss=1.,CondProbaZeroLoss_indep=1.;
	const std::vector<double>& beta = GetUniqueBeta();
	const std::vector<double>& barrier = GetBarrier();
	const std::vector<double>& pdef = GetPdefAtMaturity();

	const std::vector<double>& beta_t2 = GetUniqueBeta_t2();
	const std::vector<double>& barrier_t2 = GetBarrier_t2();
	const std::vector<double>& pdef_t2 = GetPdefAtMaturity_t2();

	its_ProbCond->SetElt(0,its_ind_x,0,CondProbaZeroLoss);
	its_ProbCond_indep->SetValue(0,0,CondProbaZeroLoss);

	for(int k=1;k<=its_nbnames;k++) 
	{
		tmp_barrier=(barrier[k-1]-beta[k-1]*x*its_sqrt2)/sqrt(1.-beta[k-1]*beta[k-1]);
		pk=NAG_cumul_normal(tmp_barrier);

		if (IsTSR() && GetT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	(barrier_t2[k-1]-beta_t2[k-1]*x*its_sqrt2)/sqrt(1.-beta_t2[k-1]*beta_t2[k-1]);
			pk += NAG_cumul_normal(tmp_barrier_t2);

			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
		}

		CondProbaZeroLoss*=(1-pk);
		CondProbaZeroLoss_indep*=(1-pdef[k-1]);
		its_ProbCond->SetElt(k,its_ind_x,0,CondProbaZeroLoss);
		its_ProbCond_indep->SetValue(k,0,CondProbaZeroLoss_indep);
	}

	return (its_ProbCond->Elt(its_nbnames,its_ind_x,0));
}


double ICM_Gauss1FLossDistrib::compute_cond_distrib(const double& x)
{

	its_ind_x+=1;

	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0.,tmp3=0.; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	int wk=0,ind_lastloss=0;
	double tmp_barrier=0.,tmp4=0.;

	const std::vector<double>& beta = GetUniqueBeta();
	const std::vector<double>& barrier = GetBarrier();
	const std::vector<double>& pdef = GetPdefAtMaturity();

	const std::vector<double>& beta_t2 = GetUniqueBeta_t2();
	const std::vector<double>& barrier_t2 = GetBarrier_t2();
	const std::vector<double>& pdef_t2 = GetPdefAtMaturity_t2();

	its_ProbCond->SetElt(0,its_ind_x,its_ind_Loss,0.);
	its_ProbCond_indep->SetValue(0,its_ind_Loss,0.);

	if (its_ishomogeneous)
	{

		for(int k=1;k<=its_nbnames;k++) 
			{
			tmp_barrier=(barrier[k-1]-beta[k-1]*x*its_sqrt2)/sqrt(1.-beta[k-1]*beta[k-1]);
			pk=NAG_cumul_normal(tmp_barrier);

			if (IsTSR() && GetT2())
			{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	(barrier_t2[k-1]-beta_t2[k-1]*x*its_sqrt2)/sqrt(1.-beta_t2[k-1]*beta_t2[k-1]);
			pk += NAG_cumul_normal(tmp_barrier_t2);

			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
			}

			tmp1=its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss);

			its_ProbCond->SetElt(k,its_ind_x,its_ind_Loss, its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss)*(1-pk)+pk*its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss-1));
			its_ProbCond_indep->SetValue(k,its_ind_Loss, its_ProbCond_indep->Getvalue(k-1,its_ind_Loss)*(1-pdef[k-1])+pdef[k-1]*its_ProbCond_indep->Getvalue(k-1,its_ind_Loss-1));

			tmp3=its_ProbCond->Elt(k,its_ind_x,its_ind_Loss);

			}
	}
	else
	{
		for(int k=1;k<=its_nbnames;k++) 
			{

			ind_lastloss=its_ind_Loss-its_int_lossrates[k-1];
			tmp_barrier=(barrier[k-1]-beta[k-1]*x*its_sqrt2)/sqrt(1.-beta[k-1]*beta[k-1]);
			pk=NAG_cumul_normal(tmp_barrier);

			if (IsTSR() && GetT2())
			{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	= (barrier_t2[k-1]-beta_t2[k-1]*x*its_sqrt2)/sqrt(1.-beta_t2[k-1]*beta_t2[k-1]);
			pk += NAG_cumul_normal(tmp_barrier_t2);

			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
			}

			tmp1=its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss);

			if(ind_lastloss<0)
			its_ProbCond->SetElt(k,its_ind_x,its_ind_Loss, its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss)*(1-pk));
			else
			its_ProbCond->SetElt(k,its_ind_x,its_ind_Loss, its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss)*(1-pk)+pk*its_ProbCond->Elt(k-1,its_ind_x,ind_lastloss));

			tmp3=its_ProbCond->Elt(k,its_ind_x,its_ind_Loss);
			}
	}

	return its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss);
}


void ICM_Gauss1FLossDistrib::compute_distrib(const int& lup)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	its_lup	=	0;
	max_lup	=	lup;
	
	for (l=its_lup;l<=max_lup;l++) 
	{
		if (l ==0)
		{
			its_ind_x=-1;
			tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::compute_cond_distribZeroLoss, (*this))).Integrate_20();	
			tmpPrbLoss*=its_normPi;

			its_lossdistrib[0]=(1-its_proportion_independant-its_proportion_full_correl)*tmpPrbLoss+its_proportion_full_correl*(1-its_min_pdef)+its_proportion_independant*its_ProbCond_indep->Getvalue(its_nbnames,l);
			its_ind_Loss=1;
		}
		else
		{
			its_ind_x=-1;

			tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::compute_cond_distrib, (*this))).Integrate_20();
			tmpPrbLoss*=its_normPi;

			its_lossdistrib[l]=(1-its_proportion_independant-its_proportion_full_correl)*tmpPrbLoss+its_proportion_independant*its_ProbCond_indep->Getvalue(its_nbnames,l);

			its_ind_Loss+=1;
		}
	}

	its_ind_Loss--;

	for (l=0;l<=lup;l++) cumul_distrib+=its_lossdistrib[l];

	its_taildistrib=1-cumul_distrib;
}


double ICM_Gauss1FLossDistrib::compute_expectedlosstranche(const double& tranche_up, 
														   const double& tranche_down, 
														   const double& lossunit,
														   //ICM_QMatrix<double>* ShiftMatrix,
														   // const int& Tenor,
														   // const int& Issuer,
														   vector<double>& losses)
{

	// if ((GetShifts())&&(Tenor!=-1)&&(Issuer!=-1)) 
	// 	return GetShifts()->Getvalue(Issuer,Tenor);


	double exp_loss_tranche=0.;
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	int i=0;	

	SetLup(0);

	if (lup>0)
	{
		its_ProbCond->ResizeWithCopy(lup+1,0.);
		its_ProbCond_indep->Resize(its_nbnames+1,lup+1);

		std::vector<double>  Intermediate;
		Intermediate.resize(its_lossdistrib.size());

		for (i=0; i<its_lup; i++)
			Intermediate[i] = its_lossdistrib[i];
		
		its_lossdistrib.resize(lup+1);
		for (i=0; i<its_lup; i++)
			its_lossdistrib[i] = Intermediate[i];
	}

	compute_distrib(lup);

	for(int l=ldown+1;l<=lup;l++) 
		exp_loss_tranche += its_lossdistrib[l]*(l*lossunit-tranche_down);

	exp_loss_tranche += (tranche_up-tranche_down)*its_taildistrib;

	return (exp_loss_tranche);
}


/**
17783 
void ICM_Gauss1FLossDistrib::compute_expectedlosstranche_fast_spread_hedge(const double& tranche_up, 
																	  const double& tranche_down, 
																	  const double& lossunit,
																	  ICM_QMatrix<double>* ShiftMatrix)
{
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);

	its_lup	=	0;

	its_ProbCond->ResizeWithCopy(lup+1,0.);
	its_ProbCond_indep->Resize(its_nbnames+1,lup+1);
	its_lossdistrib.resize(lup+1);

	compute_distrib(lup);
	compute_expectedlosstranche_perturb(tranche_up,tranche_down,lossunit,ShiftMatrix);
}
**/ 
// for hedge computation
/**
17783 
double ICM_Gauss1FLossDistrib::compute_cond_distribZeroLoss_shift_k(const double& x)
{

	double pk=0.,pshiftk=0.,tmp_barrier_perturb=0.,tmp_barrier=0.,CondProbaZeroLossshiftk=1.;

	its_ind_x+=1;

	double	beta_value;
	double	SQRT_one_minus_beta_square;

	const std::vector<double>& beta = GetUniqueBeta();

	beta_value	=	beta[its_ind_name];
	SQRT_one_minus_beta_square	=	1.0	- beta_value * beta_value;
	SQRT_one_minus_beta_square	=	sqrt(SQRT_one_minus_beta_square);

	tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-beta_value*x*its_sqrt2)/SQRT_one_minus_beta_square;
	pshiftk=NAG_cumul_normal(tmp_barrier_perturb);

	tmp_barrier=(its_barrier[its_ind_name]-beta_value*x*its_sqrt2)/SQRT_one_minus_beta_square;
	pk=NAG_cumul_normal(tmp_barrier);

	CondProbaZeroLossshiftk=its_ProbCond->Elt(its_nbnames,its_ind_x,0) / (1-pk);

	its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,0,CondProbaZeroLossshiftk);

	return (CondProbaZeroLossshiftk*(1-pshiftk));
}

double ICM_Gauss1FLossDistrib::compute_cond_distrib_shift_k(const double& x)
{

	its_ind_x+=1;

	double pk=0.,pk_perturb=0.,tmp=0.;     // proba cond de défaut du nom k
	double tmp_barrier=0.,tmp_barrier_perturb=0.;
	double Prb_Cond_l_k=0.;  // proba cond de loss=l pour un panier sans le name k
	double Prb_Cond_l_shiftk=0.;  // proba cond de loss=l pour un panier où le name k a été shifté
	double Prb_Cond_last_l_k=0.;  // proba cond de loss=l-1 pour un panier sans le name k
	int ind_lastloss=0;

	const std::vector<double>& beta = GetUniqueBeta();
	const std::vector<double>& barrier = GetBarrier();
	const std::vector<double>& pdef = GetPdefAtMaturity();

	if (its_ishomogeneous)
	{
		tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-beta[its_ind_name]*x*its_sqrt2)/sqrt(1.-beta[its_ind_name]*beta[its_ind_name]);
		pk_perturb=NAG_cumul_normal(tmp_barrier_perturb);
		tmp_barrier=(barrier[its_ind_name]-beta[its_ind_name]*x*its_sqrt2)/sqrt(1.-beta[its_ind_name]*beta[its_ind_name]);
		pk=NAG_cumul_normal(tmp_barrier);
		tmp=(its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss)-its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,its_ind_Loss-1)*pk)/(1-pk);
		its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,its_ind_Loss,tmp);
		Prb_Cond_l_shiftk=tmp*(1-pk_perturb)+its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,its_ind_Loss-1)*pk_perturb;
	}
	else
	{
		ind_lastloss=its_ind_Loss-its_int_lossrates[its_ind_name];
		tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-beta[its_ind_name]*x*its_sqrt2)/sqrt(1.-beta[its_ind_name]*beta[its_ind_name]);
		pk_perturb=NAG_cumul_normal(tmp_barrier_perturb);
		tmp_barrier=(barrier[its_ind_name]-beta[its_ind_name]*x*its_sqrt2)/sqrt(1.-beta[its_ind_name]*beta[its_ind_name]);
		pk=NAG_cumul_normal(tmp_barrier);

		if (ind_lastloss<0)
		{
			tmp=(its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss))/(1-pk);
			its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,its_ind_Loss,tmp);
			Prb_Cond_l_shiftk=tmp*(1-pk_perturb);

		}
		else
		{
			tmp=(its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss)-its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,ind_lastloss)*pk)/(1-pk);
			its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,its_ind_Loss,tmp);
			Prb_Cond_l_shiftk=tmp*(1-pk_perturb)+its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,ind_lastloss)*pk_perturb;
		}

	}

	return (Prb_Cond_l_shiftk);
}


void ICM_Gauss1FLossDistrib::compute_distrib_perturb(const int& lup)
{

	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;

	its_lup	=	0;

	for (int k=0;k<its_nbnames;k++) 
	{
		its_ind_x=-1;
		its_ind_name=k;
		max_lup	=	lup;
		tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::compute_cond_distribZeroLoss_shift_k, (*this))).Integrate_20();	
		tmpPrbLoss*=its_normPi;
		its_lossdistrib_perturb->SetValue(0,k,tmpPrbLoss);
		cumul_distrib=tmpPrbLoss;
		its_ind_Loss=1;

		for (int l=its_lup+1;l<=max_lup;l++) 
		{
			 its_ind_x=-1;
			 tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::compute_cond_distrib_shift_k, (*this))).Integrate_20();
			 tmpPrbLoss*=its_normPi;
			 cumul_distrib+=tmpPrbLoss;
	 		 its_lossdistrib_perturb->SetValue(l,k,tmpPrbLoss);
			 its_ind_Loss+=1;
		}

		its_taildistrib_perturb[k]=1-cumul_distrib;
	}
}
**/ 
/**
17783 
void ICM_Gauss1FLossDistrib::compute_expectedlosstranche_perturb(const double& tranche_up, 
															  const double& tranche_down, 
															  const double& lossunit,
															  ICM_QMatrix<double>* ShiftMatrix)
{

	int nbscenario=ShiftMatrix->Getnbrows();		// inversion
	ICM_QMatrix<double>* result= new ICM_QMatrix<double>(its_nbnames,nbscenario);

	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	int s,k,l;	
	double tmp=0.;
	its_lup	=	0;

	for (s=0;s<nbscenario;s++)  //boucle sur les scenarios de hedges par buckets
	{
		ARM_Vector * V = ShiftMatrix->RowAsVector(s) ; // ->DupRowAsTab(s);
		SetPdefPerturb(*V);delete V;		// LJ, isn't it its_nbnames?
		compute_barrier_perturb();
		its_ProbCond_Perturb->ResizeWithCopy(lup+1,0.);
		its_lossdistrib_perturb->Resize(lup+1,its_nbnames);

		compute_distrib_perturb(lup);

		for (k=0;k<its_nbnames;k++)  //boucle sur les scenarios de hedges par labels
		{
			result->SetValue(k,s,0.);
			tmp=0.;

			for (l=ldown+1;l<=lup;l++) 
				tmp+=(its_lossdistrib_perturb->Getvalue(l,k))*(l*lossunit-tranche_down);

			tmp+=(tranche_up-tranche_down)*its_taildistrib_perturb[k];
			result->SetValue(k,s,tmp);
		}
	}

	// 17783 SetShifts(result);
}
**/ 

void ICM_Gauss1FLossDistrib::BitwiseCopy(const ARM_Object* src)
{
    ICM_Gauss1FLossDistrib* srcdistrib = (ICM_Gauss1FLossDistrib *) src;

	if (its_ProbCond_indep) delete its_ProbCond_indep;
	its_ProbCond_indep = NULL;

	if (srcdistrib->its_ProbCond_indep)
		its_ProbCond_indep = (ICM_QMatrix<double>*) srcdistrib->its_ProbCond_indep->Clone();
	
	its_proportion_full_correl = srcdistrib->its_proportion_full_correl;
	its_proportion_independant = srcdistrib->its_proportion_independant;

}

// -------------
//	Copy Method 
// -------------
void ICM_Gauss1FLossDistrib::Copy(const ARM_Object* src)
{
	ICM_Distribution::Copy(src);
    BitwiseCopy(src);
}

ARM_Object* ICM_Gauss1FLossDistrib::Clone(void)
{
     ICM_Gauss1FLossDistrib* theClone = new ICM_Gauss1FLossDistrib();

     theClone->Copy(this);
 
     return(theClone);
}


double ICM_Gauss1FLossDistrib::ComputeEL_FullHomog(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& beta_down,
												const double& beta_up,
												const double& Pdef,
												const int& IntStep,
												vector<double>& losses)
{

	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;
	losses.clear();

	int lup=MIN(floor(tranche_up/LossUnit),its_nbnames);
	int ldown=MIN(floor(tranche_down/LossUnit),its_nbnames);

	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
	{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
	}

	double its_FH_Pdef = Pdef;

	if (fabs(Pdef) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
	else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
	else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

	vector<double> tmpPrbLoss_down;tmpPrbLoss_down.resize(lup+1);
	vector<double> tmpPrbLoss_up;tmpPrbLoss_up.resize(lup+1);
	double CumtmpPrbLoss=0.;
	double Loss = 0.;
	double Cnp = 0;
	int i=0;

	for (i=0;i<=lup;i++)
	{
		its_FH_i = i;
		Cnp = itsCombinations(its_nbnames,i);
		tmpPrbLoss_down[i]=0.;

		its_FH_Beta = beta_up;
		tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
		tmpPrbLoss_up[i]*=Cnp*its_normPi;
		CumtmpPrbLoss+=tmpPrbLoss_up[i];
	}

	for (i=0;i<=ldown;i++)
	{
		its_FH_i = i;
		Cnp = itsCombinations(its_nbnames,i);

		its_FH_Beta = beta_down;
		tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
		tmpPrbLoss_down[i]*=Cnp*its_normPi;
	}

	for (int i2=0;i2<=lup;i2++) 
	{
		Loss+=(tmpPrbLoss_up[i2]-tmpPrbLoss_down[i2])*(i2*LossUnit-tranche_down);
		losses.push_back(tmpPrbLoss_up[i2]);
	}
	Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);

	return (Loss);
}


double ICM_Gauss1FLossDistrib::ZeroLossFullHomog(const double& x)
{
	double pk=0.,tmp_barrier=0.;
	double retour=0;

	tmp_barrier=(its_FH_Barrier-its_FH_Beta*x*its_sqrt2)/sqrt(1.-its_FH_Beta*its_FH_Beta);
	pk=NAG_cumul_normal(tmp_barrier);

	retour = pow(pk,its_FH_i)*pow(1.-pk,static_cast<int>(its_nbnames-its_FH_i));

	return (retour);
}


// ----------------------------------------------------------------------------------------------
// Light Distrib
// ----------------------------------------------------------------------------------------------

double ICM_LightDistrib::ZeroLossFullHomog(const double& x)
{
	double pk=0.,tmp_barrier=0.;
	double retour=0;

	if (its_IsUp)
	{

		tmp_barrier=(its_FH_Barrier- its_FH_Beta*x*its_sqrt2)/sqrt(1.-its_FH_Beta*its_FH_Beta);
		pk=NAG_cumul_normal(tmp_barrier);

	//forward collateral
	if (its_FH_collat_Barrier)
	{
		tmp_barrier=(its_FH_collat_Barrier-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier);
	}

	}
	else
	{

	tmp_barrier=(its_FH_Barrier_down-its_FH_Beta*x*its_sqrt2)/sqrt(1.-its_FH_Beta*its_FH_Beta);
	pk=NAG_cumul_normal(tmp_barrier);

	//forward collateral
	if (its_FH_collat_Barrier)
	{
		tmp_barrier=(its_FH_collat_Barrier_down-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier);
	}

	}

	if (its_FullProba)
		retour = pow(pk,its_FH_i)*pow(1-pk,its_nbnames-its_FH_i);
	else
		retour=pk;

	return (retour);
}

double ICM_LightDistrib::ZeroLossFullHomog_TSR(const double& x)
{
	double pk=0.,tmp_barrier=0.,tmp_barrier1=0.,tmp_barrier2=0.;
	double retour=0;

	if (CHECK_EQUAL(its_T2,0.)){
		return ZeroLossFullHomog(x);}

	if (its_IsUp)
	{

	tmp_barrier1=(its_FH_Barrier-its_FH_Beta*x*its_sqrt2)/sqrt(1.-its_FH_Beta*its_FH_Beta);
	pk=NAG_cumul_normal(tmp_barrier1);

	tmp_barrier2=(its_FH_Barrier_T-its_FH_Beta_T*x*its_sqrt2)/sqrt(1.-its_FH_Beta_T*its_FH_Beta_T);
	pk+=NAG_cumul_normal(tmp_barrier2);

	pk-=NAG_bivariate_normal_dist(tmp_barrier1,tmp_barrier2,RHO_DIST(its_T1,its_T2));

	//forward collateral Term Structure
	if (its_FH_collat_Barrier && its_collat_T2)
	{
		tmp_barrier1=(its_FH_collat_Barrier-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier1);

		tmp_barrier2=(its_FH_collat_Barrier_T-its_FH_collat_Beta_T*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta_T*its_FH_collat_Beta_T);
		pk-=NAG_cumul_normal(tmp_barrier2);

		pk+=NAG_bivariate_normal_dist(tmp_barrier1,tmp_barrier2,RHO_DIST(its_collat_T1,its_collat_T2));
	} else if (its_FH_collat_Barrier) {
		tmp_barrier=(its_FH_collat_Barrier-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier);
	}

	}else
	{

	tmp_barrier1=(its_FH_Barrier_down-its_FH_Beta*x*its_sqrt2)/sqrt(1.-its_FH_Beta*its_FH_Beta);
	pk=NAG_cumul_normal(tmp_barrier1);

	tmp_barrier2=(its_FH_Barrier_down_T-its_FH_Beta_T*x*its_sqrt2)/sqrt(1.-its_FH_Beta_T*its_FH_Beta_T);
	pk+=NAG_cumul_normal(tmp_barrier2);

	pk-=NAG_bivariate_normal_dist(tmp_barrier1,tmp_barrier2,RHO_DIST(its_T1,its_T2));

	//forward collateral Term Structure
	if (its_FH_collat_Barrier && its_collat_T2)
	{
		tmp_barrier1=(its_FH_collat_Barrier_down-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier1);

		tmp_barrier2=(its_FH_collat_Barrier_down_T-its_FH_collat_Beta_T*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta_T*its_FH_collat_Beta_T);
		pk-=NAG_cumul_normal(tmp_barrier2);

		pk+=NAG_bivariate_normal_dist(tmp_barrier1,tmp_barrier2,RHO_DIST(its_collat_T1,its_collat_T2));
	} else if (its_FH_collat_Barrier) {
		tmp_barrier=(its_FH_collat_Barrier_down-its_FH_collat_Beta*x*its_sqrt2)/sqrt(1.-its_FH_collat_Beta*its_FH_collat_Beta);
		pk-=NAG_cumul_normal(tmp_barrier);
	}

	}

	if (its_FullProba)
		retour = pow(pk,its_FH_i)*pow(1-pk,its_nbnames-its_FH_i);
	else
		retour = pk;

	return (retour);
}

//-------------------------------------------------------------------
// Full Homogeneous case with integration
//-------------------------------------------------------------------
double ICM_LightDistrib::ComputeEL_FullHomog(const double& LossUnit,
							   const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   const int& IntStep,
							   vector<double>& losses)
{

		double Limit_case_Minus =	-10.;
		double Limit_case_Plus	=	10.;
		losses.clear();

		int lup=MIN(floor(tranche_up/LossUnit),its_nbnames);
		int ldown=MIN(floor(tranche_down/LossUnit),its_nbnames);

		if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
		{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
		}

		if (fabs(Pdef) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
		else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
		else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

		vector<double> tmpPrbLoss_down;tmpPrbLoss_down.resize(lup+1);
		vector<double> tmpPrbLoss_up;tmpPrbLoss_up.resize(lup+1);
		double CumtmpPrbLoss=0.;
		double Loss = 0.;
		double Cnp = 0;
		int i=0;

		for (i=0;i<=lup;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);
			tmpPrbLoss_down[i]=0.;

			its_FH_Beta = beta_up;
			tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
			tmpPrbLoss_up[i]*=Cnp*its_normPi;
			CumtmpPrbLoss+=tmpPrbLoss_up[i];
		}

		for (i=0;i<=ldown;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);

			its_FH_Beta = beta_down;
			tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
			tmpPrbLoss_down[i]*=Cnp*its_normPi;
		}

		for (int i2=0;i2<=lup;i2++) 
		{
			Loss+=(tmpPrbLoss_up[i2]-tmpPrbLoss_down[i2])*(i2*LossUnit-tranche_down);
			losses.push_back(tmpPrbLoss_up[i2]);
		}
		Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);

		return (Loss);
}

//-------------------------------------------------------------------
// Full Homogeneous case with integration TSR
//-------------------------------------------------------------------
double ICM_LightDistrib::ComputeEL_FullHomog_TSR(const double& LossUnit,
										const double& tranche_down,
										const double& tranche_up,
										const double& beta_down_T1,
										const double& beta_down_T2,
										const double& beta_up_T1,
										const double& beta_up_T2,
										const double& Pdef_T1,
										const double& Pdef_T2,
										const double& Pdef_down_T1,
										const double& Pdef_down_T2,
										const double& T1,
										const double& T2,
										const int& IntStep,
										vector<double>& losses)
{

	its_T1=T1;
	its_T2=T2;

	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;
	losses.clear();

	int lup=MIN(floor(tranche_up/LossUnit),its_nbnames);
	int ldown=MIN(floor(tranche_down/LossUnit),its_nbnames);

	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
	{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
	}

	if (fabs(Pdef_T1) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
	else if (fabs(Pdef_T1-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
	else its_FH_Barrier = NAG_deviates_normal_dist(Pdef_T1);

	if (fabs(Pdef_T2) < DB_TOL) its_FH_Barrier_T = Limit_case_Minus;
	else if (fabs(Pdef_T2-1.0) < DB_TOL) its_FH_Barrier_T = Limit_case_Plus;
	else its_FH_Barrier_T = NAG_deviates_normal_dist(Pdef_T2);

	if (fabs(Pdef_down_T1) < DB_TOL) its_FH_Barrier_down = Limit_case_Minus;
	else if (fabs(Pdef_down_T1-1.0) < DB_TOL) its_FH_Barrier_down = Limit_case_Plus;
	else its_FH_Barrier_down = NAG_deviates_normal_dist(Pdef_down_T1);

	if (fabs(Pdef_down_T2) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Minus;
	else if (fabs(Pdef_down_T2-1.0) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Plus;
	else its_FH_Barrier_down_T = NAG_deviates_normal_dist(Pdef_down_T2);

	vector<double> tmpPrbLoss_down;tmpPrbLoss_down.resize(lup+1);
	vector<double> tmpPrbLoss_up;tmpPrbLoss_up.resize(lup+1);
	double CumtmpPrbLoss=0.;
	double Loss = 0.;
	double Cnp = 0;
	int i=0;

	its_IsUp=true;

		for (i=0;i<=lup;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);
			tmpPrbLoss_down[i]=0.;

			its_FH_Beta = beta_up_T1;
			its_FH_Beta_T = beta_up_T2;
			tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(IntStep);
			tmpPrbLoss_up[i]*=Cnp*its_normPi;
			CumtmpPrbLoss+=tmpPrbLoss_up[i];
		}

	its_IsUp=false;

		for (i=0;i<=ldown;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);

			its_FH_Beta = beta_down_T1;
			its_FH_Beta_T = beta_down_T2;
			tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(IntStep);
			tmpPrbLoss_down[i]*=Cnp*its_normPi;
		}

		for (int i2=0;i2<=lup;i2++) 
		{
			Loss+=(tmpPrbLoss_up[i2]-tmpPrbLoss_down[i2])*(i2*LossUnit-tranche_down);
			losses.push_back(tmpPrbLoss_up[i2]);
		}
		Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);
		losses.push_back(1.-CumtmpPrbLoss);

		return (Loss);
}

//-------------------------------------------------------------------
// LHP
//-------------------------------------------------------------------

double ICM_LightDistrib::ComputeEL_LHP(const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   double recovery,
							   vector<double>& losses)
{
		losses.clear();

		if (fabs(Pdef) < DB_TOL) its_FH_Barrier = -10.;
		else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = 10.;
		else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

		double ExpectedLoss = 0.,PrbLoss_down=0.,PrbLoss_up=0.;
		double K1=0.,K2=0.;
		double DevNorm1=0.,DevNorm2=0.;

		//calcul des l'expected loss de la tranche
		if (tranche_up>tranche_down)
		{
		K1 = tranche_down/(1.-recovery);
		if (CHECK_EQUAL(K1, 0.0)) {DevNorm1 = _MINUS_INFINITY_;}
		else if (CHECK_EQUAL(K1, 1.0)) {DevNorm1 = _PLUS_INFINITY_;}
		else DevNorm1 = NAG_deviates_normal(K1);

		K2 = tranche_up/(1.-recovery);
		if (CHECK_EQUAL(K2, 0.0)) {DevNorm2 = _MINUS_INFINITY_;}
		else if (CHECK_EQUAL(K2, 1.0)) {DevNorm2 = _PLUS_INFINITY_;}
		else DevNorm2 = NAG_deviates_normal(K2);

		PrbLoss_up = NAG_bivariate_normal_dist(-DevNorm2,
								its_FH_Barrier,-sqrt(1.-beta_up*beta_up));


		PrbLoss_down = NAG_bivariate_normal_dist(-DevNorm1,
								its_FH_Barrier,-sqrt(1.-beta_down*beta_down));


		ExpectedLoss=(PrbLoss_down-PrbLoss_up)/(K2-K1);
		}
		else //calcul de l'expected loss du portfeuille
		{
		ExpectedLoss = NAG_cumul_normal(
								(1./beta_down)*
								(sqrt(1.-beta_down*beta_down)*NAG_deviates_normal_dist(tranche_down)
								-its_FH_Barrier));

		}

		return (ExpectedLoss);
}

//-------------------------------------------------------------------
// Full Homogeneous case with integration
//-------------------------------------------------------------------
double ICM_LightDistrib::ComputeEL_FullHomog_collat(const double& LossUnit,
							   const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   const double& collat_beta_down,
							   const double& collat_beta_up,
							   const double& collat_Pdef,
							   const int& IntStep,
							   vector<double>& losses)
{

		double Limit_case_Minus =	-10.;
		double Limit_case_Plus	=	10.;
		losses.clear();

		int lup=MIN(floor(tranche_up/LossUnit),its_nbnames);
		int ldown=MIN(floor(tranche_down/LossUnit),its_nbnames);

		if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
		{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
		}

		if (fabs(Pdef) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
		else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
		else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

		if (fabs(collat_Pdef) < DB_TOL) its_FH_collat_Barrier = Limit_case_Minus;
		else if (fabs(collat_Pdef-1.0) < DB_TOL) its_FH_collat_Barrier = Limit_case_Plus;
		else its_FH_collat_Barrier = NAG_deviates_normal_dist(collat_Pdef);

		vector<double> tmpPrbLoss_down;tmpPrbLoss_down.resize(lup+1);
		vector<double> tmpPrbLoss_up;tmpPrbLoss_up.resize(lup+1);
		double CumtmpPrbLoss=0.;
		double Loss = 0.;
		double Cnp = 0;
		int i=0;

		for (i=0;i<=lup;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);
			tmpPrbLoss_down[i]=0.;

			its_FH_Beta = beta_up;
			its_FH_collat_Beta = collat_beta_up;
			tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
			tmpPrbLoss_up[i]*=Cnp*its_normPi;
			CumtmpPrbLoss+=tmpPrbLoss_up[i];
		}

		for (i=0;i<=ldown;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);

			its_FH_Beta = beta_down;
			its_FH_collat_Beta = collat_beta_down;
			tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog, (*this))).Integrate(IntStep);
			tmpPrbLoss_down[i]*=Cnp*its_normPi;
		}

		for (int i2=0;i2<=lup;i2++) 
		{
			Loss+=(tmpPrbLoss_up[i2]-tmpPrbLoss_down[i2])*(i2*LossUnit-tranche_down);
			losses.push_back(tmpPrbLoss_up[i2]);
		}
		Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);

		return (Loss);
}

//-------------------------------------------------------------------
// Full Homogeneous case with integration TSR
//-------------------------------------------------------------------
double ICM_LightDistrib::ComputeEL_FullHomog_collat_TSR(const double& LossUnit,
										const double& tranche_down,
										const double& tranche_up,
										const double& beta_down_T1,
										const double& beta_down_T2,
										const double& beta_up_T1,
										const double& beta_up_T2,
										const double& Pdef_T1,
										const double& Pdef_T2,
										const double& Pdef_down_T1,
										const double& Pdef_down_T2,
										const double& T1,
										const double& T2,
										const double& collat_beta_down_T1,
										const double& collat_beta_down_T2,
										const double& collat_beta_up_T1,
										const double& collat_beta_up_T2,
										const double& collat_Pdef_T1,
										const double& collat_Pdef_T2,
										const double& collat_Pdef_down_T1,
										const double& collat_Pdef_down_T2,
										const double& collat_T1,
										const double& collat_T2,
										const int& IntStep,
										vector<double>& losses)
{

		its_T1=T1;
		its_T2=T2;

		its_collat_T1=collat_T1;
		its_collat_T2=collat_T2;

		double Limit_case_Minus =	-10.;
		double Limit_case_Plus	=	10.;
		losses.clear();

		int lup=MIN(floor(tranche_up/LossUnit),its_nbnames);
		int ldown=MIN(floor(tranche_down/LossUnit),its_nbnames);

		if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
		{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
		}

		if (fabs(Pdef_T1) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
		else if (fabs(Pdef_T1-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
		else its_FH_Barrier = NAG_deviates_normal_dist(Pdef_T1);

		if (fabs(Pdef_T2) < DB_TOL) its_FH_Barrier_T = Limit_case_Minus;
		else if (fabs(Pdef_T2-1.0) < DB_TOL) its_FH_Barrier_T = Limit_case_Plus;
		else its_FH_Barrier_T = NAG_deviates_normal_dist(Pdef_T2);

		//bariers down
		if (fabs(Pdef_down_T1) < DB_TOL) its_FH_Barrier_down = Limit_case_Minus;
		else if (fabs(Pdef_down_T1-1.0) < DB_TOL) its_FH_Barrier_down = Limit_case_Plus;
		else its_FH_Barrier_down = NAG_deviates_normal_dist(Pdef_down_T1);

		if (fabs(Pdef_down_T2) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Minus;
		else if (fabs(Pdef_down_T2-1.0) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Plus;
		else its_FH_Barrier_down_T = NAG_deviates_normal_dist(Pdef_down_T2);

		if (fabs(collat_Pdef_T1) < DB_TOL) its_FH_collat_Barrier = Limit_case_Minus;
		else if (fabs(collat_Pdef_T1-1.0) < DB_TOL) its_FH_collat_Barrier = Limit_case_Plus;
		else its_FH_collat_Barrier = NAG_deviates_normal_dist(collat_Pdef_T1);

		if (fabs(collat_Pdef_T2) < DB_TOL) its_FH_collat_Barrier_T = Limit_case_Minus;
		else if (fabs(collat_Pdef_T2-1.0) < DB_TOL) its_FH_collat_Barrier_T = Limit_case_Plus;
		else its_FH_collat_Barrier_T = NAG_deviates_normal_dist(collat_Pdef_T2);

		//bariers down
		if (fabs(collat_Pdef_down_T1) < DB_TOL) its_FH_collat_Barrier_down= Limit_case_Minus;
		else if (fabs(collat_Pdef_down_T1-1.0) < DB_TOL) its_FH_collat_Barrier_down = Limit_case_Plus;
		else its_FH_collat_Barrier_down = NAG_deviates_normal_dist(collat_Pdef_down_T1);

		if (fabs(collat_Pdef_down_T2) < DB_TOL) its_FH_collat_Barrier_down_T = Limit_case_Minus;
		else if (fabs(collat_Pdef_down_T2-1.0) < DB_TOL) its_FH_collat_Barrier_down_T = Limit_case_Plus;
		else its_FH_collat_Barrier_down_T = NAG_deviates_normal_dist(collat_Pdef_down_T2);

		vector<double> tmpPrbLoss_down;tmpPrbLoss_down.resize(lup+1);
		vector<double> tmpPrbLoss_up;tmpPrbLoss_up.resize(lup+1);
		double CumtmpPrbLoss=0.;
		double Loss = 0.;
		double Cnp = 0;
		int i=0;

		its_IsUp=true;

		for (i=0;i<=lup;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);
			tmpPrbLoss_down[i]=0.;

			its_FH_Beta = beta_up_T1;
			its_FH_Beta_T = beta_up_T2;
			its_FH_collat_Beta = collat_beta_up_T1;
			its_FH_collat_Beta_T = collat_beta_up_T2;
			tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(IntStep);
			tmpPrbLoss_up[i]*=Cnp*its_normPi;
			CumtmpPrbLoss+=tmpPrbLoss_up[i];
		}

		its_IsUp=false;

		for (i=0;i<=ldown;i++)
		{
			its_FH_i = i;
			Cnp = itsCombinations(its_nbnames,i);

			its_FH_Beta = beta_down_T1;
			its_FH_Beta_T = beta_down_T2;
			its_FH_collat_Beta = collat_beta_down_T1;
			its_FH_collat_Beta_T = collat_beta_down_T2;
			tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(IntStep);
			tmpPrbLoss_down[i]*=Cnp*its_normPi;
		}

		for (int i2=0;i2<=lup;i2++) 
		{
			Loss+=(tmpPrbLoss_up[i2]-tmpPrbLoss_down[i2])*(i2*LossUnit-tranche_down);
			losses.push_back(tmpPrbLoss_up[i2]);
		}
		Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);
		losses.push_back(1.-CumtmpPrbLoss);

		return (Loss);
}


double ICM_LightDistrib::ComputeEL_FullHomog(CtxtDistrib* ctxt)
{
	double output=0.;

	SetCtxt(ctxt);
	GetCtxt()->ComputeALLBarriers();

	switch (ctxt->m_DistribType)
	{
		case qDISTRIB_COLLFWD:
			{
				output = ComputeEL_FullHomog_collat(ctxt->m_LossUnit,
							   ctxt->m_TrancheDown,
							   ctxt->m_TrancheUp,
							   ctxt->m_beta_dw_maturity[0],
							   ctxt->m_beta_up_maturity[0],
							   ctxt->m_pdef_maturity[0],
							   ctxt->m_collat_beta_dw_start[0],
							   ctxt->m_collat_beta_up_start[0],
							   ctxt->m_collat_pdef_start[0],
							   ctxt->m_IntegrationStep,
							   ctxt->m_OutLosses);
				break;
			}
		case qDISTRIB_STD_TSR:
			{
				output = ComputeEL_FullHomog_TSR(ctxt->m_LossUnit,
							   ctxt->m_TrancheDown,
							   ctxt->m_TrancheUp,
							   ctxt->m_ts_beta_dw_maturity_t1[0],
							   ctxt->m_ts_beta_dw_maturity_t2[0],
							   ctxt->m_ts_beta_up_maturity_t1[0],
							   ctxt->m_ts_beta_up_maturity_t2[0],
							   ctxt->m_ts_pdef_up_maturity_t1[0],
							   ctxt->m_ts_pdef_up_maturity_t2[0],
							   ctxt->m_ts_pdef_dw_maturity_t1[0],
							   ctxt->m_ts_pdef_dw_maturity_t2[0],
							   ctxt->m_ts_t1_corr,
							   ctxt->m_ts_t2_corr,
							   ctxt->m_IntegrationStep,
							   ctxt->m_OutLosses);
				break;
			}
		case qDISTRIB_COLLFWD_TSR:
			{
				output = ComputeEL_FullHomog_collat_TSR(ctxt->m_LossUnit,
							   ctxt->m_TrancheDown,
							   ctxt->m_TrancheUp,
							   ctxt->m_ts_beta_dw_maturity_t1[0],
							   ctxt->m_ts_beta_dw_maturity_t2[0],
							   ctxt->m_ts_beta_up_maturity_t1[0],
							   ctxt->m_ts_beta_up_maturity_t2[0],
							   ctxt->m_ts_pdef_up_maturity_t1[0],
							   ctxt->m_ts_pdef_up_maturity_t2[0],
							   ctxt->m_ts_pdef_dw_maturity_t1[0],
							   ctxt->m_ts_pdef_dw_maturity_t2[0],
							   ctxt->m_ts_t1_corr,
							   ctxt->m_ts_t2_corr,
							   ctxt->m_ts_collat_beta_dw_start_t1[0],
							   ctxt->m_ts_collat_beta_dw_start_t2[0],
							   ctxt->m_ts_collat_beta_up_start_t1[0],
							   ctxt->m_ts_collat_beta_up_start_t2[0],
							   ctxt->m_ts_collat_pdef_up_start_t1[0],
							   ctxt->m_ts_collat_pdef_up_start_t2[0],
							   ctxt->m_ts_collat_pdef_dw_start_t1[0],
							   ctxt->m_ts_collat_pdef_dw_start_t2[0],
							   ctxt->m_ts_collat_start_t1_corr,
							   ctxt->m_ts_collat_start_t2_corr,	
							   ctxt->m_IntegrationStep,
							   ctxt->m_OutLosses);
				break;
			}
		case qDISTRIB_STD:
		default :
			{
				output = ComputeEL_FullHomog(ctxt->m_LossUnit,
							   ctxt->m_TrancheDown,
							   ctxt->m_TrancheUp,
							   ctxt->m_beta_dw_maturity[0],
							   ctxt->m_beta_up_maturity[0],
							   ctxt->m_pdef_maturity[0],	
							   ctxt->m_IntegrationStep,
							   ctxt->m_OutLosses);
				break;
			}
	}

	return output;
}



double ICM_LightDistrib::ZeroLossFullHomog_TSR_stepup(const double& x)
{
	its_FullProba=false;
	
	int lup=MIN(floor(GetCtxt()->m_TrancheUp/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	int ldown=MIN(floor(GetCtxt()->m_TrancheDown/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	double Z_v1_v2_j=0.;
	double Pi_j_u=0.,Pi_j_s=0.;
	double Pu=0.,Ps=0.;
	double result = 0.;
	int j=0;

	if (its_IsUp) {
	its_FH_Beta = GetCtxt()->m_ts_beta_up_maturity_t1[0];
	its_FH_Beta_T = GetCtxt()->m_ts_beta_up_maturity_t2[0];
	its_FH_Barrier = GetCtxt()->m_ts_barriers_up_maturity_t1[0];
	its_FH_Barrier_T= GetCtxt()->m_ts_barriers_up_maturity_t2[0];
	Pi_j_u=ZeroLossFullHomog_TSR(x);
	
	its_FH_Beta = GetCtxt()->m_ts_stepup_beta_up_maturity_t1[0];
	its_FH_Beta_T = GetCtxt()->m_ts_stepup_beta_up_maturity_t2[0];
	its_FH_Barrier = GetCtxt()->m_ts_stepup_barriers_up_maturity_t1[0];
	its_FH_Barrier_T= GetCtxt()->m_ts_stepup_barriers_up_maturity_t2[0];
	Pi_j_s=ZeroLossFullHomog_TSR(x);
	}
	else 
	{
	its_FH_Beta = GetCtxt()->m_ts_beta_dw_maturity_t1[0];
	its_FH_Beta_T = GetCtxt()->m_ts_beta_dw_maturity_t2[0];
	its_FH_Barrier_down= GetCtxt()->m_ts_barriers_dw_maturity_t1[0];
	its_FH_Barrier_down_T= GetCtxt()->m_ts_barriers_dw_maturity_t2[0];
	Pi_j_u=ZeroLossFullHomog_TSR(x);

	its_FH_Beta = GetCtxt()->m_ts_stepup_beta_dw_maturity_t1[0];
	its_FH_Beta_T = GetCtxt()->m_ts_stepup_beta_dw_maturity_t2[0];
	its_FH_Barrier_down= GetCtxt()->m_ts_stepup_barriers_dw_maturity_t1[0];
	its_FH_Barrier_down_T= GetCtxt()->m_ts_stepup_barriers_dw_maturity_t2[0];
	Pi_j_s=ZeroLossFullHomog_TSR(x);
	}

	Pu=Pi_j_u;
	Ps=Pi_j_s;

	double value = 0.,cum=0.;
	int nbnames = GetCtxt()->m_nbnames;

	for (j=1;j<nbnames+1;j++)
	{
		for (int v1=0;v1<= MIN(j,GetCtxt()->m_losses_reset_beg);v1++)
			for (int v2=v1;v2<= MIN(j,GetCtxt()->m_losses_reset_end) ;v2++)
		{
			Z_v1_v2_j = 0.;

			if (v1>=1) 
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1) + 
						Ps*its_ProbCond_stepup->Elt(v1-1,v2-1,j-1) +
						(Pu - Ps)*(value=its_ProbCond_stepup->Elt(v1,v2-1,j-1));
			}
			else if (v2<1)
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1);
			}	
			else if ((v1<1) && (1<=v2))
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1) + 
						(Pu - Ps)*its_ProbCond_stepup->Elt(v1,v2-1,j-1);
			}

		its_ProbCond_stepup->SetElt(v1,v2,j,Z_v1_v2_j);
		}
	} 

	result = its_ProbCond_stepup->Elt(GetCtxt()->m_losses_reset_beg,GetCtxt()->m_losses_reset_end,nbnames);

	its_FullProba=true;

	return (result);
}


double ICM_LightDistrib::ComputeStepUp(CtxtDistrib* ctxt)
{
	SetCtxt(ctxt);
	GetCtxt()->ComputeALLBarriers();

	its_T1=ctxt->m_ts_t1_corr;
	its_T2=ctxt->m_ts_t2_corr;
	its_nbnames=GetCtxt()->m_nbnames;

	vector<double> tmpPrbLoss_up;
	vector<double> tmpPrbLoss_down;
	double CumtmpPrbLoss=0.;
	double CumtmpPrbLoss_=0.;
	double Loss = 0.,Loss1 = 0.,Loss2 = 0.;
	int i=0;
	double Cnp=0.;

	//-------------------------------------------------------------------
	// partie non step up
	//-------------------------------------------------------------------
	int lup=MIN(floor(GetCtxt()->m_TrancheUp/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	int ldown=MIN(floor(GetCtxt()->m_TrancheDown/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);

	double stepup = GetCtxt()->m_TrancheUp_stepup;
	double stepdw = GetCtxt()->m_TrancheDown_stepup;

	int lup_stepup=MIN(floor(stepup/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	int ldown_stepup=MIN(floor(stepdw/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);

	tmpPrbLoss_up.resize(lup_stepup+1,0.);
	tmpPrbLoss_down.resize(lup_stepup+1,0.);

	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup_stepup+1)))
	{	itsCombinations.Resize(MAX(its_nbnames+1,lup_stepup+1),MAX(GetCtxt()->m_nbnames+1,lup_stepup+1));
		itsCombinations.GenCombinations();	}

	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;

	if (fabs(GetCtxt()->m_ts_stepup_pdef_up_maturity_t1[0]) < DB_TOL) its_FH_Barrier = Limit_case_Minus;
	else if (fabs(GetCtxt()->m_ts_stepup_pdef_up_maturity_t1[0]-1.0) < DB_TOL) its_FH_Barrier = Limit_case_Plus;
	else its_FH_Barrier = NAG_deviates_normal_dist(GetCtxt()->m_ts_stepup_pdef_up_maturity_t1[0]);

	if (fabs(GetCtxt()->m_ts_stepup_pdef_up_maturity_t2[0]) < DB_TOL) its_FH_Barrier_T = Limit_case_Minus;
	else if (fabs(GetCtxt()->m_ts_stepup_pdef_up_maturity_t2[0]-1.0) < DB_TOL) its_FH_Barrier_T = Limit_case_Plus;
	else its_FH_Barrier_T = NAG_deviates_normal_dist(GetCtxt()->m_ts_stepup_pdef_up_maturity_t2[0]);

	if (fabs(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t1[0]) < DB_TOL) its_FH_Barrier_down = Limit_case_Minus;
	else if (fabs(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t1[0]-1.0) < DB_TOL) its_FH_Barrier_down = Limit_case_Plus;
	else its_FH_Barrier_down = NAG_deviates_normal_dist(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t1[0]);

	if (fabs(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t2[0]) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Minus;
	else if (fabs(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t2[0]-1.0) < DB_TOL) its_FH_Barrier_down_T = Limit_case_Plus;
	else its_FH_Barrier_down_T = NAG_deviates_normal_dist(GetCtxt()->m_ts_stepup_pdef_dw_maturity_t2[0]);

	its_IsUp=true;
	for (i=0;i<=lup_stepup;i++)
	{
		its_FH_i = i;
		Cnp = itsCombinations(its_nbnames,i);
		tmpPrbLoss_down[i]=0.;
		its_FH_Beta = GetCtxt()->m_ts_stepup_beta_up_maturity_t1[0];
		its_FH_Beta_T = GetCtxt()->m_ts_stepup_beta_up_maturity_t2[0];
		tmpPrbLoss_up[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		tmpPrbLoss_up[i]*=Cnp*its_normPi;
		CumtmpPrbLoss+=tmpPrbLoss_up[i];
		if (i<=lup) {CumtmpPrbLoss_=CumtmpPrbLoss;}
	}

	its_IsUp=false;
	for (i=0;i<=ldown_stepup;i++)
	{
		its_FH_i = i;
		Cnp = itsCombinations(its_nbnames,i);
		its_FH_Beta = GetCtxt()->m_ts_stepup_beta_dw_maturity_t1[0];
		its_FH_Beta_T = GetCtxt()->m_ts_stepup_beta_dw_maturity_t2[0];
		tmpPrbLoss_down[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		tmpPrbLoss_down[i]*=Cnp*its_normPi;
	}

	//-------------------------------------------------------------------
	// partie step up
	//-------------------------------------------------------------------
	double LossUp1 = 0.,LossDown1 = 0.;
	vector<double> tmpPrbLoss_up_stepup;
	vector<double> tmpPrbLoss_down_stepup;
	double CumtmpPrbLossUp_stepup=0.;
	double CumtmpPrbLossDown_stepup=0.;

	if (its_ProbCond_stepup)
		delete its_ProbCond_stepup;
	its_ProbCond_stepup= new ICM_QCubix<double>(GetCtxt()->m_nbnames+1,GetCtxt()->m_nbnames+1,GetCtxt()->m_nbnames+1,0.);
	its_ProbCond_stepup->SetElt(0,0,0,1.);

	vector<double> CheckLoss1;
	vector<double> CheckLoss2;

	for (i=0;i<=lup;i++)
	{
		double proba_dw = tmpPrbLoss_down[i];
		if (i>ldown) 
			{proba_dw=0.;}
		Loss2 +=(tmpPrbLoss_up[i]-proba_dw)*(i*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown);
		CheckLoss2.push_back(Loss2);
	}
	Loss2 += (1.-CumtmpPrbLoss_)*(GetCtxt()->m_TrancheUp-GetCtxt()->m_TrancheDown);
	Loss2 /= (GetCtxt()->m_TrancheUp-GetCtxt()->m_TrancheDown);

	for (GetCtxt()->m_losses_reset_beg=0;GetCtxt()->m_losses_reset_beg<=lup_stepup;GetCtxt()->m_losses_reset_beg++)
	{
		CumtmpPrbLossUp_stepup=0.;
		CumtmpPrbLossDown_stepup=0.;
		LossDown1 = LossUp1 =0.;

		double CDOloss = MAX((GetCtxt()->m_TrancheUp - GetCtxt()->m_TrancheDown)-
						MAX(GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown,0.),0.)/
						(GetCtxt()->m_TrancheUp - GetCtxt()->m_TrancheDown);

		double stepdw_ = MAX(stepdw-GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit,0.);

		double stepup_ = stepdw_ 
				+ MAX(
						(GetCtxt()->m_TrancheUp_stepup - GetCtxt()->m_TrancheDown_stepup)-
						MAX(GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown_stepup,0.)
				,0.);

		if ((!stepup_)||(!CDOloss))  
			{continue;}

		double LossProbaUp = tmpPrbLoss_up[GetCtxt()->m_losses_reset_beg];
		double LossProbaDown = tmpPrbLoss_down[GetCtxt()->m_losses_reset_beg];

		tmpPrbLoss_up_stepup.resize(lup_stepup+1,0.);
		tmpPrbLoss_down_stepup.resize(lup_stepup+1,0.);
	
		its_IsUp=true;
		for (i=GetCtxt()->m_losses_reset_beg;i<=lup_stepup;i++)
		{
		GetCtxt()->m_losses_reset_end=i;
		tmpPrbLoss_up_stepup[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR_stepup, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		tmpPrbLoss_up_stepup[i]*=its_normPi;
		CumtmpPrbLossUp_stepup+=tmpPrbLoss_up_stepup[i];
		}

		its_IsUp=false;
		for (i=GetCtxt()->m_losses_reset_beg;i<=ldown_stepup;i++)
		{
		GetCtxt()->m_losses_reset_end=i;
		tmpPrbLoss_down_stepup[i]=HermiteIntegration(ff1::mem_call(&ICM_LightDistrib::ZeroLossFullHomog_TSR_stepup, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		tmpPrbLoss_down_stepup[i]*=its_normPi;
		CumtmpPrbLossDown_stepup+=tmpPrbLoss_down_stepup[i];
		}

		for (i=GetCtxt()->m_losses_reset_beg;i<=ldown_stepup;i++) 
		{
			LossDown1+=tmpPrbLoss_down_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit; 
			LossUp1+=tmpPrbLoss_up_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit;
		}

		for (i=MAX(GetCtxt()->m_losses_reset_beg,ldown_stepup+1);i<=lup_stepup;i++) 
		{LossUp1+=tmpPrbLoss_up_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit;}

		LossDown1+=(LossProbaDown-CumtmpPrbLossDown_stepup)*stepdw_;
		LossUp1+=(LossProbaUp-CumtmpPrbLossUp_stepup)*stepup_;

		double StepUpLoss = LossUp1-LossDown1;

		Loss1 += CDOloss*StepUpLoss/(stepup_-stepdw_);

		//forward case
		//Loss1 += (LossUp1-LossDown1)/(GetCtxt()->m_TrancheUp_stepup - GetCtxt()->m_TrancheDown_stepup);
		CheckLoss1.push_back(CDOloss*StepUpLoss/(stepup_-stepdw_));

	}

	if (Loss1<0.)
		{Loss1 = 0.;}

	Loss = Loss2 + Loss1 ;

	return (Loss);
}



