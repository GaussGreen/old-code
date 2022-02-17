#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_RF.h"
#include "ICMKernel\glob\icm_maths.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void
ICM_Gauss1FLossDistrib_RF::Init()
{
	//Integration
	itsIntegrationMethod	=	qGAUSS_HERMITE;		// Gauss-Legendre, Hermite, etc.
	itsIntegrationStep		=	20;
	SetIntegratorType(itsIntegrationMethod, itsIntegrationStep);

	//Copula
	itsCopulaType		=	qNO_COPULA;
	itsFreedomDegree	=	4;

	//CorrelFunction
	itsCorrelationType = qCAL_PWC_CORREL;
	itsCorrelFunctionParam.clear();
	itsCorrelFunction=NULL;

	//Precompute coef
	its_coeff_a = NULL;
	its_coeff_b = NULL;
}

//----------------------------------------------------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------------------------------------------------

ICM_Gauss1FLossDistrib_RF::ICM_Gauss1FLossDistrib_RF(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  beta,
							const ARM_Vector&  LossRates,
							const std::vector<int>& SortedIndice,
							qCAL_INDEX_CORR_TYPE CorrelType,
							const std::vector<double>& CorrelParam,
							const int& discretizationstep,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep)
{
	Init();
	Set(nbnames,pdefault,beta,LossRates, SortedIndice, CorrelType, CorrelParam, discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);
}

ICM_Gauss1FLossDistrib_RF::ICM_Gauss1FLossDistrib_RF(const int& nbnames,
							qCAL_INDEX_CORR_TYPE CorrelType,
							const std::vector<double>& CorrelParam,
							const int& discretizationstep,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep)
{
	Init();
	Set(nbnames, CorrelType, CorrelParam, discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);
}

void ICM_Gauss1FLossDistrib_RF::Set(const int& nbnames,
			 const ARM_Vector&  pdef,
			 const ARM_Vector&  beta,
			 const ARM_Vector&  LossRates,
			 const std::vector<int>& SortedIndice,
			 qCAL_INDEX_CORR_TYPE CorrelType,
			 const std::vector<double>& CorrelParam,
			 const int& discretizationstep,
			 const int& CopulaType,
			 const qIntegratorChoice&	IntegrationMethod,
			 const int& IntegrationStep)
{
	
	//Precompute coef size
	int size;
	if (discretizationstep	==	0) size	=	20;
	else size	=	discretizationstep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		if (its_coeff_a==NULL) its_coeff_a = new ICM_QMatrix<double>();
		if (its_coeff_b==NULL) its_coeff_b = new ICM_QMatrix<double>();

		its_coeff_a->Resize(nbnames, IntegrationStep);
		its_coeff_b->Resize(nbnames, IntegrationStep);
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) ||  (!its_ProbCond))
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		// 17783 if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
		// 17783 its_ProbCond_Perturb= new ICM_QCubix<double>(nbnames+1,size,1,0.);

	}

	SetNbNames(nbnames);

	//MaJ des paramètres Beta : dans ce modèle elle ne sert à rien (mais indispensable de les mettre 
	//à jour)
	SetUniqueBeta(beta);

	//Copule
	itsCopulaType = CopulaType;
	
	//Integration
	itsIntegrationStep = IntegrationStep;
	itsIntegrationMethod = IntegrationMethod;
	
	//Correl Function (always 1 factor model )
	if (itsCorrelFunction) delete itsCorrelFunction;
	int FactorSize = 1;
	itsCorrelationType = CorrelType;
	itsCorrelFunctionParam = CorrelParam;
	std::vector<double> TmpCorrel = itsCorrelFunctionParam;
	TmpCorrel.resize(TmpCorrel.size()-1);

	double vol = 0.;
	if  ((itsCorrelFunctionParam.size()%2)==0)
	{
		vol = itsCorrelFunctionParam[itsCorrelFunctionParam.size()-1];
		TmpCorrel.resize(TmpCorrel.size()-1);
	}
	else
	{	vol=-999.;	}


	itsCorrelFunction = new ICM_CorrelFunction(FactorSize, (qCopula_TYPE)itsCopulaType,  itsCorrelationType, TmpCorrel,vol);

	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegrationMethod == qGAUSS_LEGENDRE))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegrationMethod == qGAUSS_HERMITE))
		itsIntegrationMethod = qGAUSS_LEGENDRE;

	//MaJ de l'Integrator
	SetIntegratorType(itsIntegrationMethod, itsIntegrationStep);
	
	//Marginal Default Prob
	SetPdefAtMaturity(pdef); 
	compute_min_pdef(GetPdefAtMaturity());

	//Loss Amounts
	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);

	//Copula & Model Dependant
	compute_barrier();		// Copula dependant
	precompute_coeffs();

	//VA intermediaire si Long Short ds le portefeuille
	ComputeNbNameShort(LossRates);
	SetSortedIndices(SortedIndice);
}


void ICM_Gauss1FLossDistrib_RF::BitwiseCopy(const ARM_Object* src)
{

    ICM_Gauss1FLossDistrib_RF* srcdistrib = (ICM_Gauss1FLossDistrib_RF *) src;

	//Copula
	itsCopulaType = srcdistrib->itsCopulaType;
	itsFreedomDegree = srcdistrib->itsFreedomDegree;

	//Integration
	itsIntegrationStep = srcdistrib->itsIntegrationStep;
	itsIntegrationMethod = srcdistrib->itsIntegrationMethod;

	//MaJ de l'Integrator
	SetIntegratorType(itsIntegrationMethod, itsIntegrationStep);
	
	//Correlation Function
	itsCorrelationType=srcdistrib->itsCorrelationType;
	itsCorrelFunctionParam=srcdistrib->itsCorrelFunctionParam;
	if (srcdistrib->itsCorrelFunction)
		itsCorrelFunction = (ICM_CorrelFunction*) srcdistrib->itsCorrelFunction->Clone();

	//Precompute coef
	if (srcdistrib->its_coeff_a)
		its_coeff_a = (ICM_QMatrix<double>*) srcdistrib->its_coeff_a->Clone();
	if (srcdistrib->its_coeff_b)
		its_coeff_b = (ICM_QMatrix<double>*) srcdistrib->its_coeff_b->Clone();
}


void ICM_Gauss1FLossDistrib_RF::Set(const int& nbnames,
			 qCAL_INDEX_CORR_TYPE CorrelType,
			 const std::vector<double>& CorrelParam,
			 const int& discretizationstep,
			 const int& CopulaType,
			 const qIntegratorChoice&	IntegrationMethod,
			 const int& IntegrationStep)
{
	//Precompute coef size
	int size;
	if (discretizationstep	==	0) size	=	20;
	else size	=	discretizationstep;

	//Precompute coef are not used in this case
	SetNbNames(nbnames);
	itsCopulaType = CopulaType;
	
	//Integration
	itsIntegrationStep = IntegrationStep;
	itsIntegrationMethod = IntegrationMethod;
	
	//Correl Function (always 1 factor model )
	if (itsCorrelFunction) delete itsCorrelFunction;
	int FactorSize = 1;
	itsCorrelationType = CorrelType;
	itsCorrelFunctionParam = CorrelParam;
	std::vector<double> TmpCorrel = itsCorrelFunctionParam;

	double vol = 0.;
	if  ((itsCorrelFunctionParam.size()%2)==0)
	{
		vol = itsCorrelFunctionParam[itsCorrelFunctionParam.size()-1];
		TmpCorrel.resize(TmpCorrel.size()-1);
	}
	else
	{	vol=-999.;	}

	itsCorrelFunction = new ICM_CorrelFunction(FactorSize, (qCopula_TYPE)itsCopulaType,  itsCorrelationType, TmpCorrel,vol);

	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegrationMethod == qGAUSS_LEGENDRE))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegrationMethod == qGAUSS_HERMITE))
		itsIntegrationMethod = qGAUSS_LEGENDRE;

	//MaJ de l'Integrator
	SetIntegratorType(itsIntegrationMethod, itsIntegrationStep);
}

// -------------
//	Copy Method 
// -------------
void ICM_Gauss1FLossDistrib_RF::Copy(const ARM_Object* src)
{
	ICM_Distribution::Copy(src);
    BitwiseCopy(src);
}


ARM_Object* ICM_Gauss1FLossDistrib_RF::Clone(void)
{
     ICM_Gauss1FLossDistrib_RF* theClone = new ICM_Gauss1FLossDistrib_RF();

     theClone->Copy(this);
 
     return(theClone);
}



// -------------------------------------------------------------------
// GAUSSIAN COPULA
// -------------------------------------------------------------------

void
ICM_Gauss1FLossDistrib_RF::compute_barrier()
{
	double Limit_case_Minus = -10.;
	double Limit_case_Plus	= 10.;

	int k;

	its_barrier.resize(its_nbnames);

	switch (itsCopulaType)
	{
	case qNO_COPULA:	
	case qGAUSSIAN:
		if (check_homogeneous_Pdef(its_pdef_at_maturity))
		{
			double seuil = 0.;
			if (fabs(its_pdef_at_maturity[0]) < DB_TOL)
				seuil = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity[0]-1.0) < DB_TOL)
				seuil = Limit_case_Plus;
			else
				seuil = itsCorrelFunction->EstimeSeuil(its_pdef_at_maturity[0]);
			for (k=0;k<its_nbnames;k++)
				its_barrier[k] = seuil;
		}
		else
		{
			for (k=0;k<its_nbnames;k++)
			{
				if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
					its_barrier[k] = Limit_case_Minus;
				else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
					its_barrier[k] = Limit_case_Plus;
				else
					its_barrier[k] = itsCorrelFunction->EstimeSeuil(its_pdef_at_maturity[k]);
			}
		}
		
		break;
	case qSTUDENT:
		for (k=0;k<its_nbnames;k++)
			// Inverse Cumulative Student: nag_t_deviate
			its_barrier[k]	=	itsCorrelFunction->EstimeSeuil(its_pdef_at_maturity[k]);
		break;
	}
}

void ICM_Gauss1FLossDistrib_RF :: precompute_coeffs()
{
	int i=0, j=0;
	double	beta_value, a, b, factor, m, nu;

	for (i=0; i<its_nbnames; i++)
	{
		//Init Correl Function Value
		m = itsCorrelFunction->GetM();
		nu = itsCorrelFunction->GetNU();

		for (j=0; j<itsIntegrationStep; j++)
		{
			//Init Beta value  (always Gauss Hermite 40)
			//factor = TheIntegrator.GetAbscissa(j)*SQRT2;
			factor = TheIntegrator.GetRealAbscissa(j);
			beta_value = itsCorrelFunction->BetaCond(factor);

			if (fabs(beta_value) >= 1.0 || nu == 0)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
					"Beta values must be strictly between -1.0 and 1.0");
			}
			
			a	=	(its_barrier[i]-m) / nu;
			b	=	beta_value / nu;
			
			//Update coef
			its_coeff_a->SetValue(i,j,a);
			its_coeff_b->SetValue(i,j,b);
		}
	}
}




void ICM_Gauss1FLossDistrib_RF :: SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, 												
													const int&	TheStep, 
													const double& lbound,
													const double& ubound)
{
	itsIntegrationMethod	=	TheIntegratorType;
	itsIntegrationStep		=	TheStep;

	TheIntegrator.SetIntegrationType(TheIntegratorType);
	
	// Rajout du || sur la méthode qTrapeze afin que le step de l'integrator soit setter a la bonne valeur
	if ((TheIntegratorType == qGAUSS_LEGENDRE) || (TheIntegratorType == qGAUSS_HERMITE) || (TheIntegratorType == qTRAPEZE))
		TheIntegrator.SetIntegrationStep(TheStep);
	
	//MaJ des bornes d'intégration
	TheIntegrator.SetLowBound(lbound);
	TheIntegrator.SetUpBound(ubound);
}



// -------------------------------------------------------------------
//	Long Short CDO
// -------------------------------------------------------------------

double
ICM_Gauss1FLossDistrib_RF::compute_expectedlosstranche_LongShort(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  // const double& minloss
//												  ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  )
{
	double exp_loss_tranche = 0.;
	
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	int lmin= floor(fabs(its_LSminloss));

	int i=0;	
	
	its_ProbCond->ResizeWithCopy(lmin+lup+1,0.);

	its_lossdistrib.resize(lmin+lup+1);
	
	// Compute the Distribution
	compute_distrib_LongShort(lup, lmin);
	
	//View the loss ditribution
	/*#ifdef _DEBUG
	FILE *stream = fopen("c:\\temp\\LossDistribLongShort.txt", "w+");
	this->View("",stream);
	fclose(stream);
	#endif */

	int l	=	0;
	double	loss_level;

	//Long part of the portfolio (short part : no loss)
	for (l=ldown+lmin+1;l<=lup+lmin;l++)
	{
		loss_level	=	((l-lmin)*lossunit) - tranche_down;
		exp_loss_tranche		+=	its_lossdistrib[l] * loss_level;
	}
	
	exp_loss_tranche		+=	(tranche_up - tranche_down) * its_taildistrib;

	return exp_loss_tranche;
}

//Distribution

void
ICM_Gauss1FLossDistrib_RF::compute_distrib_LongShort(const int& lup, const int& lmin)
{
	switch (itsIntegrationMethod)
	{
		case qGAUSS_LEGENDRE:	
			compute_distrib_LongShortLegendre(lup, lmin);
			break;
		case qGAUSS_HERMITE:
			compute_distrib_LongShortHermite(lup, lmin);
			break;
		default:
			compute_distrib_LongShortLegendre(lup, lmin);
			break;
	}
}


void
ICM_Gauss1FLossDistrib_RF::compute_distrib_LongShortLegendre(const int& lup, const int& lmin)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0, i=0;

	int ind_x = 0;
	int ind_name = 0;
	int ind_loss = 0;
	int ind_sortedname = 0;
	max_lup		=	lup+lmin;
	its_ind_Loss = 0;
	std::vector<double> conddist;

	int nblevels = itsCorrelFunction->GetNbLevels();
	vector<double> lbound(nblevels);
	vector<double> ubound(nblevels);

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

		//Autant d'intégrations que de seuils différents si LEGENDRE
		for (i=0; i<nblevels; i++)
		{
			//Init Loss Distrib
			for (l=0;l<=max_lup;l++)
				its_lossdistrib[l]=0.;

			//First Loop : underlying factor (integration)
			for (i=0; i<nblevels; i++)
			{
				//Seuils
				if (i==0)
					lbound[i]=-6;
				else
					lbound[i] = itsCorrelFunction->GetThreshold(i-1);
				
				if (i==nblevels-1)
					ubound[i] = 6;
				else
					ubound[i] = itsCorrelFunction->GetThreshold(i);
				
				//Modif Integrateur
				SetIntegratorType(itsIntegrationMethod, itsIntegrationStep, lbound[i],ubound[i]);

				//Precompute Coef
				precompute_coeffs();

				for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
				{
					//Second Loop : portfolio size
					InitProbCondPortSize0(lmin, ind_x);
					
					//Short Issuer Part
					for (ind_name=1; ind_name<=this->its_nbnames_short; ind_name++)
					{
						//Sorted Index
						ind_sortedname = this->its_sorted_indices[ind_name-1];

						//Conditionnal loss distribution on short part
						for (ind_loss=0; ind_loss<=lmin; ind_loss++)
						{
							ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname);
						}
					}
				
					//Long Issuer Part
					for (ind_name=its_nbnames_short+1; ind_name<=this->its_nbnames; ind_name++)
					{
						//Sorted Index
						ind_sortedname = this->its_sorted_indices[ind_name-1];

						for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
						{
							ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname);
						}
					}					
				}

				//Loos distribution computation
				for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
				{
					//MaJ du vecteur comprenant les valeurs de la fonction à intégrer
					conddist = its_ProbCond->GetColVectorV(its_nbnames, ind_loss);
					
					//Intégration numérique
					its_lossdistrib[ind_loss] += TheIntegrator.IntegrateVector(conddist);
				}
			}
		}
	}

	for (l=0;l<=max_lup;l++)
		cumul_distrib += its_lossdistrib[l];

	its_taildistrib			=	1.0 - cumul_distrib;
}


void ICM_Gauss1FLossDistrib_RF::compute_distrib_LongShortHermite(const int& lup, const int& lmin)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	int ind_x = 0;
	int ind_name = 0;
	int ind_loss = 0;
	int ind_sortedname = 0;
	max_lup		=	lup+lmin;
	its_ind_Loss = 0;
	std::vector<double> conddist;

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Prob Cond computation

			//First Loop : underlying factor (integration)
			for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
			{
				//Second Loop : portfolio size
				InitProbCondPortSize0(lmin, ind_x);
				
				//Short Issuer Part
				for (ind_name=1; ind_name<=this->its_nbnames_short; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];

					//Conditionnal loss distribution on short part
					for (ind_loss=0; ind_loss<=lmin; ind_loss++)
					{
						ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname);
					}
				}
			
				//Long Issuer Part
				for (ind_name=its_nbnames_short+1; ind_name<=this->its_nbnames; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];

					for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
					{
						ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname);
					}
				}					
			}

			//Loos distribution computation
			for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
			{
				//MaJ du vecteur comprenant les valeurs de la fonction à intégrer
				conddist = its_ProbCond->GetColVectorV(its_nbnames, ind_loss);
				
				//Intégration numérique
				its_lossdistrib[ind_loss] = TheIntegrator.IntegrateVector(conddist);
			}
	}

	for (l=0;l<=max_lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}

	its_taildistrib			=	1.0 - cumul_distrib;
}


//Util Long Short
void ICM_Gauss1FLossDistrib_RF::InitProbCondPortSize0(const int& lmin, const int& ind_x)
{
	its_ProbCond->Elt(0,ind_x,lmin) = 1;
}


void ICM_Gauss1FLossDistrib_RF::ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname)
{
	double tmp1=0.;
	double tmp2=0.;
	double tmp_barrier=0.;
	double pk=0.;
	
	//double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	double x=TheIntegrator.GetRealAbscissa(ind_x);
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier
	tmp_barrier = (*its_coeff_a)(ind_sortedname,ind_x) - (*its_coeff_b)(ind_sortedname,ind_x) * x;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
	pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss >= 0)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}



void ICM_Gauss1FLossDistrib_RF::ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname)
{
	double tmp1=0.;
	double tmp2=0.;
	double tmp_barrier=0.;
	double pk=0.;
	
	//double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	double x=TheIntegrator.GetRealAbscissa(ind_x);
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier
	tmp_barrier = (*its_coeff_a)(ind_sortedname,ind_x) - (*its_coeff_b)(ind_sortedname,ind_x) * x;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
	pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss <= lmin)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}

//Compute ELoss if Full Homogene (index calibration) 
//Dépend de la méthode d'intégration

double ICM_Gauss1FLossDistrib_RF::compute_expectedlosstrancheFullHomog(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef)
{
	double res=0;
	switch (itsIntegrationMethod)
	{
		case qGAUSS_LEGENDRE:	
			res = compute_expectedlosstrancheFullHomogLegendre(LossUnit,tranche_down,tranche_up,Pdef);
			break;
		case qGAUSS_HERMITE:
			res =compute_expectedlosstrancheFullHomogHermite(LossUnit,tranche_down,tranche_up,Pdef);
			break;
		default:
			res = compute_expectedlosstrancheFullHomogLegendre(LossUnit,tranche_down,tranche_up,Pdef);
			break;
	}
	return res;
}

double ICM_Gauss1FLossDistrib_RF::compute_expectedlosstrancheFullHomogHermite(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef)
{
	//Compute Barrier
	double Limit_case_Minus = -10.;
	double Limit_case_Plus	= 10.;

	double barrier;
	if (fabs(Pdef) < DB_TOL)
		barrier = Limit_case_Minus;
	else if (fabs(Pdef-1.0) < DB_TOL)
		barrier = Limit_case_Plus;
	else
		barrier = itsCorrelFunction->EstimeSeuil(Pdef);

	//Compute coef a b
	vector<double> a(itsIntegrationStep);
	vector<double> b(itsIntegrationStep);
	vector<double> x(itsIntegrationStep);

	int i=0, j=0;
	double	beta_value, m, nu;

	//Init Correl Function Value
	m = itsCorrelFunction->GetM();
	nu = itsCorrelFunction->GetNU();

	for (j=0; j<itsIntegrationStep; j++)
	{
		x[j] = TheIntegrator.GetRealAbscissa(j);
		beta_value = itsCorrelFunction->BetaCond(x[j]);

		if (fabs(beta_value) >= 1.0 || nu == 0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		a[j]	=	(barrier-m) / nu;
		b[j]	=	beta_value / nu;
	}

	//Tranche def
	int lup=floor(tranche_up/LossUnit);
	int ldown=floor(tranche_down/LossUnit);

	//Combin computation
	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
	{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
	}

	//Compute Distrib
	vector<double> tmpPrbLoss;
	tmpPrbLoss.resize(itsIntegrationStep);
	its_lossdistrib.resize(lup+1);
	
	double cumul_distrib=0.;
	double Cnp = 0;
	double tmp_barrier=0., pk=0.;

	for (i=0;i<=lup;i++)
	{
		Cnp = itsCombinations(its_nbnames,i);
		for (j=0;j<itsIntegrationStep;j++)
		{
			tmp_barrier = a[j] - b[j] * x[j];
			pk=NAG_cumul_normal(tmp_barrier);
			tmpPrbLoss[j] = Cnp * pow(pk,i)*pow(1.-pk,static_cast<int>(its_nbnames-i));
		}
		its_lossdistrib[i] = TheIntegrator.IntegrateVector(tmpPrbLoss);
		cumul_distrib += its_lossdistrib[i];
	}

	its_taildistrib			=	1.0 - cumul_distrib;

	//View the loss ditribution
	/*#ifdef _DEBUG
	FILE *stream = fopen("c:\\temp\\LossDistribLongShort.txt", "w+");
	this->View("",stream);
	fclose(stream);
	#endif */

	//Compute Expected Loss
	int l	=	0;
	double	loss_level, exp_loss_tranche=0.;

	//Long part of the portfolio (short part : no loss)
	for (l=ldown+1;l<=lup;l++)
	{
		loss_level	=	(l*LossUnit) - tranche_down;
		exp_loss_tranche		+=	its_lossdistrib[l] * loss_level;
	}
	
	exp_loss_tranche		+=	(tranche_up - tranche_down) * its_taildistrib;

	return exp_loss_tranche;
}


double ICM_Gauss1FLossDistrib_RF::compute_expectedlosstrancheFullHomogLegendre(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef)
{
	//Compute Barrier
	double Limit_case_Minus = -10.;
	double Limit_case_Plus	= 10.;

	double barrier;
	if (fabs(Pdef) < DB_TOL)
		barrier = Limit_case_Minus;
	else if (fabs(Pdef-1.0) < DB_TOL)
		barrier = Limit_case_Plus;
	else
		barrier = itsCorrelFunction->EstimeSeuil(Pdef);

	int i=0, j=0, k=0;
	double	beta_value, m, nu;
	int nblevels = itsCorrelFunction->GetNbLevels();

	//Compute coef a b
	vector<vector<double> > a(nblevels);
	vector<vector<double> > b(nblevels);
	vector<vector<double> > x(nblevels);
	vector<double> lbound(nblevels);
	vector<double> ubound(nblevels);

	//Init Correl Function Value
	m = itsCorrelFunction->GetM();
	nu = itsCorrelFunction->GetNU();

	//Autant d'intégrations que de seuils différents si LEGENDRE
	for (i=0; i<nblevels; i++)
	{
		//Size
		a[i].resize(itsIntegrationStep);
		b[i].resize(itsIntegrationStep);
		x[i].resize(itsIntegrationStep);

		//Seuils
		if (i==0)
			lbound[i]=-6;
		else
			lbound[i] = itsCorrelFunction->GetThreshold(i-1);
		
		if (i==nblevels-1)
			ubound[i] = 6;
		else
			ubound[i] = itsCorrelFunction->GetThreshold(i);
		
		//Modif Integrateur
		SetIntegratorType(itsIntegrationMethod, itsIntegrationStep, lbound[i],ubound[i]);

		//Compute
		for (j=0; j<itsIntegrationStep; j++)
		{
			x[i][j] = TheIntegrator.GetRealAbscissa(j);
			beta_value = itsCorrelFunction->BetaCond(x[i][j]);

			if (fabs(beta_value) >= 1.0 || nu == 0)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
					"Beta values must be strictly between -1.0 and 1.0");
			}
			
			a[i][j]	=	(barrier-m) / nu;
			b[i][j]	=	beta_value / nu;
		}
	}

	//Tranche def
	int lup=floor(tranche_up/LossUnit);
	int ldown=floor(tranche_down/LossUnit);

	//Combin computation
	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(its_nbnames+1,lup+1)))
	{
		itsCombinations.Resize(MAX(its_nbnames+1,lup+1),MAX(its_nbnames+1,lup+1));
		itsCombinations.GenCombinations();
	}

	//Compute Distrib
	vector<double> tmpPrbLoss;
	tmpPrbLoss.resize(itsIntegrationStep);
	its_lossdistrib.resize(lup+1);
	
	double cumul_distrib=0.;
	double Cnp = 0;
	double tmp_barrier=0., pk=0.;

	for (i=0;i<=lup;i++)
	{
		Cnp = itsCombinations(its_nbnames,i);
		its_lossdistrib[i] = 0;
		for (k=0; k<nblevels; k++)
		{
			SetIntegratorType(itsIntegrationMethod, itsIntegrationStep, lbound[k],ubound[k]);
			for (j=0;j<itsIntegrationStep;j++)
			{
				tmp_barrier = a[k][j] - b[k][j] * x[k][j];
				pk=NAG_cumul_normal(tmp_barrier);
				tmpPrbLoss[j] = Cnp * pow(pk,i)*pow(1.-pk,static_cast<int>(its_nbnames-i));
			}
			its_lossdistrib[i] += TheIntegrator.IntegrateVector(tmpPrbLoss);
		}
		cumul_distrib += its_lossdistrib[i];
	}

	its_taildistrib			=	1.0 - cumul_distrib;

	//View the loss ditribution
	/*#ifdef _DEBUG
	FILE *stream = fopen("c:\\temp\\LossDistribLongShort.txt", "w+");
	this->View("",stream);
	fclose(stream);
	#endif */

	//Compute Expected Loss
	int l	=	0;
	double	loss_level, exp_loss_tranche=0.;

	//Long part of the portfolio (short part : no loss)
	for (l=ldown+1;l<=lup;l++)
	{
		loss_level	=	(l*LossUnit) - tranche_down;
		exp_loss_tranche		+=	its_lossdistrib[l] * loss_level;
	}
	
	exp_loss_tranche		+=	(tranche_up - tranche_down) * its_taildistrib;

	return exp_loss_tranche;
}


// View 

void ICM_Gauss1FLossDistrib_RF::View(char* id, FILE* ficOut)
{	
	FILE* fOut;	
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 

	int size =0; 

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Fast Loss Dsitrib Gaussian  ----------------- \n");
	int ind_loss=0;
	for (ind_loss=0; ind_loss<its_lossdistrib.size(); ind_loss++)
		fprintf(fOut,"ind_x : %i\t%f\n",ind_loss, its_lossdistrib[ind_loss]);
	
	fprintf(fOut,"tail : \t%f\n",its_taildistrib);

	/* fprintf(fOut, "\t\t\t ----------------- Prob Cond  ----------------- \n");
	for (int ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
	{
		for (ind_loss=0; ind_loss<its_lossdistrib.size(); ind_loss++)
		{
			for (int ind_name=0; ind_name<=this->its_nbnames; ind_name++)
				fprintf(fOut,"%f\t",its_ProbCond->Elt(ind_name,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}
		
		fprintf(fOut,"\n\n");
	} */
	
	
	
	
	ICM_Distribution::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
