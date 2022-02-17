
#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"


void ICM_Gauss1FLossDistrib_LJ::Set_APPROX(CtxtDistrib* c)
{
	SetCtxt(c);

	int size = c->m_IntegrationStep;

	int lup=floor(c->m_TrancheUp/c->m_LossUnit);
	int ldown=floor(c->m_TrancheDown/c->m_LossUnit);

	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ((c->m_IntegrationStep % 2 == 0) && (c->m_IntegrationMethod == qGAUSS_LEGENDRE))
		c->m_IntegrationMethod = qGAUSS_HERMITE;
	else if ((c->m_IntegrationStep % 2 == 1) && (c->m_IntegrationMethod == qGAUSS_HERMITE))
		c->m_IntegrationMethod = qGAUSS_LEGENDRE;
	
	c->ComputeAll();
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------

extern "C" void 
compute_cond_distrib_Gaussian_Integrator_APPROX(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distrib_Gaussian_APPROX(x,Param);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_Gaussian_APPROX(double x,void* params)
{

	CtxtDistrib* c0 = GetCtxt();

	double seuil = 15.;
	double result = 0.;
	int k;
	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;
	double Mu_up = 0.,Mu_down = 0.;
	double Sigma_up = 0.,Sigma_down = 0.;
	int nbnames = c0->m_nbnames;
	double tmp_barrier=0.;
	double sum_pdef_up=0.,sum_pdef_down=0.;
	double sum_pdef_up_square=0.,sum_pdef_down_square=0.;

	vector<double> Pdef_up;Pdef_up.resize(c0->m_nbnames);
	vector<double> Pdef_down;Pdef_down.resize(c0->m_nbnames);

	double TrancheUp = c0->m_TrancheUp/c0->m_TotalNotionals;
	double TrancheDown = c0->m_TrancheDown/c0->m_TotalNotionals;

	ICM_Gauss1FLossDistrib_LJ*	TheModel =	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)params)[0]);
	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod ==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	double EL_tranche_up=0.,EL_tranche_down=0.;
		
	for (k=0;k<c0->m_nbnames;k++) 
	{
		tmp_barrier=c0->m_a_up_maturity[k] - c0->m_b_up_maturity[k] * x;
		Pdef_up[k] = NAG_cumul_normal(tmp_barrier);
		sum_pdef_up+=Pdef_up[k];
		sum_pdef_up_square += Pdef_up[k]*Pdef_up[k];

		tmp_barrier=c0->m_a_dw_maturity[k] - c0->m_b_dw_maturity[k] * x;
		Pdef_down[k] = NAG_cumul_normal(tmp_barrier);
		sum_pdef_down+=Pdef_down[k];
		sum_pdef_down_square += Pdef_down[k]*Pdef_down[k];

		double Mui_up= (1.-c0->m_Recoveries[k])*Pdef_up[k]/(double)nbnames;
		Mu_up += Mui_up;

		double Mui_down= (1.-c0->m_Recoveries[k])*Pdef_down[k]/(double)nbnames;
		Mu_down += Mui_down;
	
		double Sigmai_up= Mui_up*Mui_up/Pdef_up[k]*(1.-Pdef_up[k]);
		Sigma_up += Sigmai_up;

		double Sigmai_down= Mui_down*Mui_down/Pdef_down[k]*(1.-Pdef_down[k]);
		Sigma_down += Sigmai_down;
	}

	// --------------------------------------------------------------
	//strike UP
	// --------------------------------------------------------------
	if (true /*sum_pdef_up<=seuil*/)
	{
		int m_up_m = (int) floor(TrancheUp*nbnames/(1.-c0->m_Recoveries[0])) + 1;
		int m_max = (int) floor(nbnames/(1.-c0->m_Recoveries[0])) ;

		double probaequalk,probaequalk_bef;

		probaequalk = exp(-sum_pdef_up);
		for (int i3=1;i3<m_up_m;i3++)
		{
			probaequalk *=sum_pdef_up;
		    probaequalk /= (double)i3;	
		}
		probaequalk_bef=probaequalk;

		EL_tranche_up =0.;

		for (int i1=m_up_m;i1<=m_max;i1++)
		{
			probaequalk *= sum_pdef_up;
			probaequalk /= (double)i1;
			//NAG_poisson_dist(sum_pdef_up, i1, probainfequalk, probasupk, probaequalk);
			EL_tranche_up += ((1./(double)nbnames)*(1.-c0->m_Recoveries[0])*((double)i1)-TrancheUp) * probaequalk;
		}

		//NAG_poisson_dist(sum_pdef_up, m_up_m-1, probainfequalk, probasupk, probaequalk);
		EL_tranche_up -= 0.5 * sum_pdef_up_square * (1./(double)nbnames)*(1.-c0->m_Recoveries[0])*probaequalk_bef;
	}
	else
	{
		Sigma_up = sqrt(Sigma_up);
		double tranche_up_tild = TrancheUp - Mu_up;
		EL_tranche_up = -tranche_up_tild*NAG_cumul_normal(-tranche_up_tild/Sigma_up);

		double sum_down = 0.,sum_up = 0.;

		for (k=0;k<c0->m_nbnames;k++)  
		{
			double Minus_Rec = (1.-c0->m_Recoveries[k])*(1.-c0->m_Recoveries[k])*(1.-c0->m_Recoveries[k]);
			sum_up += Minus_Rec/(nbnames*nbnames*nbnames)*Pdef_up[k]*(1.-Pdef_up[k])*(1.-2*Pdef_up[k]);
		}

		double ND_up = (1./(sqrt(2.*PI)*Sigma_up)) * exp(-0.5*(tranche_up_tild/Sigma_up)*(tranche_up_tild/Sigma_up));
		EL_tranche_up += (1./6.)*(1./(Sigma_up*Sigma_up))*tranche_up_tild*ND_up*sum_up;
	}

	// --------------------------------------------------------------
	//strike DOWN
	// --------------------------------------------------------------
	if (true /*sum_pdef_down<=seuil*/)
	{
		int m_down_m = (int) floor(TrancheDown*nbnames/(1.-c0->m_Recoveries[0])) + 1;
		int m_max = (int) floor(nbnames/(1.-c0->m_Recoveries[0])) ;

		double probaequalk,probaequalk_bef;

		probaequalk = exp(-sum_pdef_down);
		for (int i3=1;i3<m_down_m;i3++)
		{
			probaequalk *=sum_pdef_down;
		    probaequalk /= (double)i3;	
		}
		probaequalk_bef=probaequalk;

		EL_tranche_down =0.;
		for (int i2=m_down_m;i2<=m_max;i2++)
		{
			//NAG_poisson_dist(sum_pdef_down, i2, probainfequalk, probasupk, probaequalk);
			probaequalk *= sum_pdef_down;
			probaequalk /= (double)i2;
			EL_tranche_down += ((1./(double)nbnames)*(1.-c0->m_Recoveries[0])*((double)i2)-TrancheDown) * probaequalk;
		}

		//NAG_poisson_dist(sum_pdef_down, m_down_m-1, probainfequalk, probasupk, probaequalk);
		EL_tranche_down -= 0.5 * sum_pdef_down_square * (1./(double)nbnames)*(1.-c0->m_Recoveries[0])*probaequalk_bef;
	}
	else
	{
		Sigma_down = sqrt(Sigma_down);

		double tranche_down_tild = TrancheDown - Mu_down;
		EL_tranche_down = -tranche_down_tild*NAG_cumul_normal(-tranche_down_tild/Sigma_down);
		double sum_down = 0.,sum_up = 0.;

		for (k=0;k<c0->m_nbnames;k++)  
		{
			double Minus_Rec = (1.-c0->m_Recoveries[k])*(1.-c0->m_Recoveries[k])*(1.-c0->m_Recoveries[k]);
			sum_down += Minus_Rec/(nbnames*nbnames*nbnames)*Pdef_down[k]*(1.-Pdef_down[k])*(1.-2*Pdef_down[k]);
		}

		double ND_down = (1./(sqrt(2.*PI)*Sigma_down))*exp(-0.5*(tranche_down_tild/Sigma_down)*(tranche_down_tild/Sigma_down));
		EL_tranche_down += (1./6.)*(1./(Sigma_down*Sigma_down))*tranche_down_tild*ND_down*sum_down;
	}

	result = (EL_tranche_down-EL_tranche_up);
	result *= c0->m_TotalNotionals;

	return (result);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_distrib_APPROX(const int& lup,const int& ldown)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	double result =0.;

	// ---------------------
	its_lup	=	0;
//	max_lup = MAX(lup,its_lup);
	
	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;
	TheParameterVector.Append(this);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	CtxtDistrib* c = GetCtxt();

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:
		default :
		{
			TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Integrator_APPROX, &TheParameterVector, tmpPrbLoss); 
			result = tmpPrbLoss;
			//result = TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_APPROX, &TheParameterVector, tmpPrbLoss);
		} 
		break;
	}

	return (result);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_APPROX(CtxtDistrib* c)
{
	return	compute_distrib_APPROX(0., 0.);
}


void ICM_Gauss1FLossDistrib_LJ :: SetIntegratorType_APPROX(const qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep)
{
	qIntegratorChoice itsIntegrationMethod	=	TheIntegratorType;
	int itsIntegrationStep		=	TheStep;

	TheIntegrator.SetIntegrationType(TheIntegratorType);
	
	// Rajout du || sur la méthode qTrapeze afin que le step de l'integrator soit setter a la bonne valeur
	if ((TheIntegratorType == qGAUSS_LEGENDRE) || (TheIntegratorType == qGAUSS_HERMITE) || (TheIntegratorType == qTRAPEZE))
		TheIntegrator.SetIntegrationStep(TheStep);
}



