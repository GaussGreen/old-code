#include "ARMKernel\glob\firsttoinc.h" 


#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_2F.h"
#include "ICMKernel\glob\icm_maths.h"



//#include "ARMKernel\util\interpol.h"

//#include "ICMKernel\glob\icm_enums.h"
//#include "ICMKernel\glob\icm_constants.h"
//#include "ICMKernel\glob\icm_addressvector.h"


//#include "ICMKernel\util\icm_HermiteIntegration.h"
//#include "ICMKernel\util\icm_GaussLegendreIntegration.h"


//#include "ICMKernel\util\icm_qmatrix.h"
//#include "ICMKernel\mod\icm_homogene_fastlossdistrib.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void
ICM_Gaussian_LossDistrib_2F::Init()
{
	its_CopulaType	=	qGAUSSIAN;

	its_IntegrationStep_Common_Factor	=	20;
	its_IntegrationStep_Sector_Factor	=	20;

	itsIntegrationMethod_Common_Factor	=	qGAUSS_HERMITE;
	itsIntegrationMethod_Sector_Factor	=	qGAUSS_HERMITE;

	its_CorrelationType		=	TFCT_FULL;

	its_coeff_a.clear();
	its_coeff_b.clear();
	its_coeff_c.clear();
	its_coeff_den.clear();

	TheIntegrator.Reset();

	// --------------------------------------
	// SECTORIAL APPROACH
	its_sector_membership.clear();

	its_Single_intra_sector_correlation	=	0.0;
	its_Single_inter_sector_correlation	=	0.0;

	its_Betas.clear();
	its_Lambdas.clear();

	its_Nb_Sectors	=	0;

	its_SQRT_OneMinus_intra_sector_correlation	=	0.0;
	its_SQRT_inter_sector_correlation	=	0.0;
	its_SQRT_inter_minus_intra_sector_correlation	=	0.0;

	// INTEGRATION	
	its_factor_state	=	0.0;
	its_sector_state	=	0.0;

	its_index_loss		=	0;

	its_Current_DefProb.clear();

	its_lossdistrib.clear();
	its_taildistrib		=	0.0;

	its_first_passing_flag	=	true;

	
}

//----------------------------------------------------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------------------------------------------------
ICM_Gaussian_LossDistrib_2F::ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  betas_sector,
							const ARM_Vector&  lambdas_sector,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{
	Init();

	Set(nbnames, pdefault, betas_sector, lambdas_sector, sector_membership, nb_sectors, LossRates, 
				CopulaType, IntegrationMethod, IntegrationStep, IntegrationMethod_2 ,IntegrationStep_2);

	Set_First_Passing();
}

ICM_Gaussian_LossDistrib_2F::ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  betas_sector,
							double	intra_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{
	Init();

	Set(nbnames, pdefault, betas_sector, intra_sector_correl, sector_membership, nb_sectors, LossRates, 
				CopulaType, IntegrationMethod, IntegrationStep, IntegrationMethod_2 ,IntegrationStep_2);

	Set_First_Passing();
}

ICM_Gaussian_LossDistrib_2F::ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector&  pdefault,
							double	intra_sector_correl,
							double	inter_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{
	Init();

	Set(nbnames, pdefault, intra_sector_correl, inter_sector_correl, sector_membership, nb_sectors, LossRates, 
				CopulaType, IntegrationMethod, IntegrationStep, IntegrationMethod_2 ,IntegrationStep_2);

	Set_First_Passing();
}

/*********************************************************************/

// SET FOR:		TFCT_FULL

/*********************************************************************/

void ICM_Gaussian_LossDistrib_2F::Set(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  betas_sector,
							const ARM_Vector&  lambdas_sector,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{

	int i, size;
	
	if (IntegrationStep	==	0)
		size	=	20;
	else
		size	=	IntegrationStep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		its_sector_membership.resize(nbnames);
		its_Betas.resize(nbnames);
		its_Lambdas.resize(nbnames);
		
		for (i=0; i<nbnames; i++)
		{
			its_sector_membership[i]	=	sector_membership[i];
			its_Betas[i]	=	betas_sector[i];
			its_Lambdas[i]	=	lambdas_sector[i];
		}

		// ALLOCATIONS
		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_c.clear();
		its_coeff_den.clear();

		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_c.resize(nbnames);
		its_coeff_den.resize(nbnames);

		its_Barriers.clear();
		its_Barriers.resize(nbnames);

		its_Current_DefProb.clear();
		its_Current_DefProb.resize(nbnames);

		its_lossdistrib.clear();
	}

	// betas...
	// double*	betas;
	// betas	=	new double [nbnames];
	ARM_Vector betas(nbnames); 
	ICM_Gauss1FLossDistrib::Set(nbnames, pdefault, betas, LossRates);

	// modified by preceeding line
//	if ((nbnames != its_nbnames) || (its_IntegrationStep_Common_Factor != size))
//	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond	=	new ICM_QCubix<double>(nbnames + 1, IntegrationStep * IntegrationStep_2, 1,0.);
//	}

	// NB NAMES
	SetNbNames(nbnames);
	SetNbSectors(nb_sectors);

	its_CopulaType	=	qGAUSSIAN;	//	CopulaType
	
	its_IntegrationStep_Common_Factor	=	IntegrationStep;
	its_IntegrationStep_Sector_Factor	=	IntegrationStep_2;

	itsIntegrationMethod_Common_Factor	=	IntegrationMethod;
	itsIntegrationMethod_Sector_Factor	=	IntegrationMethod_2;

	its_CorrelationType	=	TFCT_FULL;
	
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	// COMPUTE BARRIERS and alike
	SetPdefAtMaturity(pdefault); 
//	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);

//	itsCopulaType	=	its_CopulaType;
	Compute_Barriers();			// Copula dependant

	Compute_Conditional_Probabilities_Coeffs();

	//delete	[] betas;
}


/*********************************************************************/

// SET FOR:		TFCT_SAME_INTER_DIFF_INTRA

/*********************************************************************/

void ICM_Gaussian_LossDistrib_2F::Set(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  betas_sector,
							double	intra_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{

	int i, size;
	
	if (IntegrationStep	==	0)
		size	=	20;
	else
		size	=	IntegrationStep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		its_sector_membership.resize(nbnames);
		its_Betas.resize(nbnames);
		its_Lambdas.resize(nbnames);
		
		for (i=0; i<nb_sectors; i++)
		{
			its_sector_membership[i]	=	sector_membership[i];
			its_Betas[i]	=	betas_sector[i];
		}

		// ALLOCATIONS
		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_c.clear();
		its_coeff_den.clear();

		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_c.resize(nbnames);
		its_coeff_den.resize(nbnames);

		its_Barriers.clear();
		its_Barriers.resize(nbnames);

		its_Current_DefProb.clear();
		its_Current_DefProb.resize(nbnames);

		its_lossdistrib.clear();
	}

	// betas...
	//double*	betas;

	//betas	=	new double [nbnames];

	ARM_Vector betas(nbnames);
	ICM_Gauss1FLossDistrib::Set(nbnames, pdefault, betas, LossRates);

	// modified by preceeding line
//	if ((nbnames != its_nbnames) || (its_IntegrationStep_Common_Factor != size))
//	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond	=	new ICM_QCubix<double>(nbnames + 1, IntegrationStep * IntegrationStep_2, 1,0.);
//	}

	// NB NAMES
	SetNbNames(nbnames);
	SetNbSectors(nb_sectors);

	its_Single_intra_sector_correlation	=	intra_sector_correl;

	its_CopulaType	=	qGAUSSIAN;	//	CopulaType
	
	its_IntegrationStep_Common_Factor	=	IntegrationStep;
	its_IntegrationStep_Sector_Factor	=	IntegrationStep_2;

	itsIntegrationMethod_Common_Factor	=	IntegrationMethod;
	itsIntegrationMethod_Sector_Factor	=	IntegrationMethod_2;

	its_CorrelationType	=	TFCT_SAME_INTER_DIFF_INTRA;
	
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	// COMPUTE BARRIERS and alike
	SetPdefAtMaturity(pdefault); 
//	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);

//	itsCopulaType	=	its_CopulaType;
	Compute_Barriers();			// Copula dependant

	Compute_Conditional_Probabilities_Coeffs();

	//delete	[] betas;
}


/*********************************************************************/

// SET FOR:		TFCT_SAME_INTER_DIFF_INTRA

/*********************************************************************/

void ICM_Gaussian_LossDistrib_2F::Set(const int& nbnames,
							const ARM_Vector&  pdefault,
							double	intra_sector_correl,
							double	inter_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep,
							const qIntegratorChoice&	IntegrationMethod_2,
							const int& IntegrationStep_2
							)
{

	int i, size;
	
	if (IntegrationStep	==	0)
		size	=	20;
	else
		size	=	IntegrationStep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		its_sector_membership.resize(nbnames);
		its_Betas.resize(nbnames);
		its_Lambdas.resize(nbnames);
		
		for (i=0; i<nb_sectors; i++)
		{
			its_sector_membership[i]	=	sector_membership[i];
		}

		// ALLOCATIONS
		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_c.clear();
		its_coeff_den.clear();

		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_c.resize(nbnames);
		its_coeff_den.resize(nbnames);

		its_Barriers.clear();
		its_Barriers.resize(nbnames);

		its_Current_DefProb.clear();
		its_Current_DefProb.resize(nbnames);

		its_lossdistrib.clear();
	}

	// betas...
	// double*	betas;
	// betas	=	new double [nbnames];

	ARM_Vector betas(nbnames); 
	ICM_Gauss1FLossDistrib::Set(nbnames, pdefault, betas, LossRates);

	// modified by preceeding line
//	if ((nbnames != its_nbnames) || (its_IntegrationStep_Common_Factor != size))
//	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond	=	new ICM_QCubix<double>(nbnames + 1, IntegrationStep * IntegrationStep_2, 1,0.);
//	}

	// NB NAMES
	SetNbNames(nbnames);
	SetNbSectors(nb_sectors);

	its_Single_intra_sector_correlation	=	intra_sector_correl;
	its_Single_inter_sector_correlation	=	inter_sector_correl;
	
	its_SQRT_OneMinus_intra_sector_correlation	=	sqrt(1.0 - intra_sector_correl);
	its_SQRT_inter_sector_correlation	=	sqrt(1.0 - inter_sector_correl);
	if (inter_sector_correl > intra_sector_correl)
		ICMTHROW(ERR_INVALID_DATA,"Invalid values Intra and Inter Correlation must satisfy Inter <= Intra! - Inter: " << inter_sector_correl << " - Intra: " << intra_sector_correl);

	its_SQRT_inter_minus_intra_sector_correlation	=	sqrt(inter_sector_correl - intra_sector_correl);

	its_CopulaType	=	qGAUSSIAN;	//	CopulaType
	
	its_IntegrationStep_Common_Factor	=	IntegrationStep;
	its_IntegrationStep_Sector_Factor	=	IntegrationStep_2;

	itsIntegrationMethod_Common_Factor	=	IntegrationMethod;
	itsIntegrationMethod_Sector_Factor	=	IntegrationMethod_2;

	its_CorrelationType	=	TFCT_SAME_INTER_SAME_INTRA;
	
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	// COMPUTE BARRIERS and alike
	SetPdefAtMaturity(pdefault); 
//	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);

//	itsCopulaType	=	its_CopulaType;
	Compute_Barriers();			// Copula dependant

	Compute_Conditional_Probabilities_Coeffs();

	//delete	[] betas;
}

/*********************************************************************/


void ICM_Gaussian_LossDistrib_2F::BitwiseCopy(const ARM_Object* src)
{

    ICM_Gaussian_LossDistrib_2F* srcdistrib = (ICM_Gaussian_LossDistrib_2F *) src;

	its_CopulaType			=	srcdistrib->its_CopulaType;

	its_IntegrationStep_Common_Factor	=	srcdistrib->its_IntegrationStep_Common_Factor;
	its_IntegrationStep_Sector_Factor	=	srcdistrib->its_IntegrationStep_Sector_Factor;
	
	its_CorrelationType		=	srcdistrib->its_CorrelationType;
	
	TheIntegrator			=	srcdistrib->TheIntegrator;

	its_coeff_a				=	srcdistrib->its_coeff_a;
	its_coeff_b				=	srcdistrib->its_coeff_b;
	its_coeff_c				=	srcdistrib->its_coeff_c;
	its_coeff_den			=	srcdistrib->its_coeff_den;

	its_sector_membership	=	srcdistrib->its_sector_membership;
	its_Single_intra_sector_correlation	=	srcdistrib->its_Single_intra_sector_correlation;
	its_Single_inter_sector_correlation	=	srcdistrib->its_Single_inter_sector_correlation;

	its_Betas	=	srcdistrib->its_Betas;
	its_Lambdas	=	srcdistrib->its_Lambdas;

	its_Nb_Sectors		=	srcdistrib->its_Nb_Sectors;

	its_SQRT_inter_sector_correlation	=	srcdistrib->its_SQRT_inter_sector_correlation;
	its_SQRT_OneMinus_intra_sector_correlation	=	srcdistrib->its_SQRT_OneMinus_intra_sector_correlation;

	its_factor_state	=	srcdistrib->its_factor_state;
	its_sector_state	=	srcdistrib->its_sector_state;

	its_index_loss		=	srcdistrib->its_index_loss;

	its_Barriers		=	srcdistrib->its_Barriers;
	its_Current_DefProb	=	srcdistrib->its_Current_DefProb;

	its_lossdistrib		=	srcdistrib->its_lossdistrib;
	its_taildistrib		=	srcdistrib->its_taildistrib;

}

// -------------
//	Copy Method 
// -------------
void ICM_Gaussian_LossDistrib_2F::Copy(const ARM_Object* src)
{
	ICM_Gauss1FLossDistrib::Copy(src);
    BitwiseCopy(src);
}


ARM_Object* ICM_Gaussian_LossDistrib_2F::Clone(void)
{
     ICM_Gaussian_LossDistrib_2F* theClone = new ICM_Gaussian_LossDistrib_2F();

     theClone->Copy(this);
 
     return(theClone);
}



void ICM_Gaussian_LossDistrib_2F::Compute_Barriers()
{
	double Limit_case_Minus = -10.;
	double Limit_case_Plus	= 10.;

	int k;

	its_Barriers.resize(its_nbnames);

	switch (its_CopulaType)
	{
	case qNO_COPULA:	
	case qGAUSSIAN:

		for (k=0; k<its_nbnames; k++)
		{
			if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
				its_Barriers[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
				its_Barriers[k] = Limit_case_Plus;
			else
				its_Barriers[k] =	NAG_deviates_normal_dist(its_pdef_at_maturity[k]);
		}
		break;

	case qSTUDENT:
		break;
	}
}

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

extern "C" void  The_Function_To_Integrate(void* Param, double x, double y, double& Result)
{
	int		id_param;
	double	*LossMin;
	double	*LossMax;
	int		*factorId;
	int		*sectorId;
	qIntegratorChoice	*TheIntegratorChoice;

	ICM_Gaussian_LossDistrib_2F	*TheModel = NULL;

//	qIntegratorChoice	*TheIntegrationMethod;
	// -------------------------------------------
	// Get parameters from Stack
	// -------------------------------------------
	
	id_param	=	0;

	TheModel			=	(ICM_Gaussian_LossDistrib_2F*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	TheIntegratorChoice	=	(qIntegratorChoice*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	LossMin		=	(double*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	LossMax		=	(double*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	factorId		=	(int*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	sectorId		=	(int*)((*(AddressVector*)Param)[id_param]);
	id_param++;

	// -------------------------------------------

	// -------------------------------------------
	// STEP 0: SET THE DATA into THE MODEL
	// -------------------------------------------
	TheModel->SetFactorState(*factorId);
	TheModel->SetSectorState(*sectorId);

	// -------------------------------------------
	// STEP 1: COMPUTE COEFFICIENTS
	// -------------------------------------------

	// -------------------------------------------
	// STEP 2: COMPUTE CONDITIONAL PROBABILITIES
	// -------------------------------------------
	TheModel->ComputeConditionalDefaultProbabilities(x, y);

	// -------------------------------------------
	// STEP 3: COMPUTE EXPECTED LOSS TRANCHE
	// -------------------------------------------
	Result	=	TheModel->ComputeExpectedLossTranche(*LossMin, *LossMax);
}


// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

double ICM_Gaussian_LossDistrib_2F::compute_expectedlosstranche(double& tranche_down, 
												  double& tranche_up, 
												  const double& lossunit,
												  //ICM_QMatrix<double>* ShiftMatrix,
												  const int& Tenor,
												  const int& Issuer)
{

	// ---------------------------------
	// ALGORITHM IS COMPLETELY MODIFIED
	// ---------------------------------

	double	Loss2F;	

	double	LossMin;
	double	LossMax;

	LossMin	=	tranche_down;
	LossMax	=	tranche_up;

	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;
	// ----------------------

	// resize arrays
	its_Current_DefProb.clear();
	its_Current_DefProb.resize(its_nbnames);

	// ----------------------------------
	// BARRIERS ARE COMPLETELY COMPUTED
	// ----------------------------------
	// ----------------------------------

	// ----------------------------------
	// INTEGRATOR
	// ----------------------------------

	qIntegratorChoice	TheIntegratorChoice;

	double	x_min	=	-6.0;
	double	x_max	=	6.0;
	int		x_nbsteps	=	21;

	double	y_min	=	-6.0;
	double	y_max	=	6.0;
	int		y_nbsteps	=	21;

	int		x_step	=	0;
	int		y_step	=	0;

	// ----------------------------------

	TheIntegratorChoice	=	itsIntegrationMethod_Common_Factor;
	
	TheIntegrator.SetIntegrationType(TheIntegratorChoice);
	TheIntegrator.SetIntegrationType_2(TheIntegratorChoice);

	if (TheIntegratorChoice == qGAUSS_LEGENDRE)
	{
		TheIntegrator.SetLowBound(x_min);
		TheIntegrator.SetUpBound(x_max);
		TheIntegrator.SetLowBound_2(y_min);
		TheIntegrator.SetUpBound_2(y_max);
	}

	TheIntegrator.SetIntegrationStep(its_IntegrationStep_Common_Factor);
	TheIntegrator.SetIntegrationStep_2(its_IntegrationStep_Sector_Factor);

	// ----------------------------------
	TheParameterVector.Reset();
	TheParameterVector.Append(this);
	TheParameterVector.Append(&TheIntegratorChoice);
	TheParameterVector.Append(&tranche_down);
	TheParameterVector.Append(&tranche_up);
	TheParameterVector.Append(&x_step);
	TheParameterVector.Append(&y_step);

	its_lossunit	=	lossunit;
	// ----------------------------------

	its_Current_DefProb.clear();
	its_Current_DefProb.resize(its_nbnames);


	TheIntegrator.Integrate(x_min, x_max, y_min, y_max, The_Function_To_Integrate, &TheParameterVector, Loss2F);

	// ----------------------------------
	
	return	Loss2F;
}


// ------------------------------------------------------------------------------------------
// Knowing Factor value z, and Sector Value y, Evaluate Default Prob for each name at time T
// ------------------------------------------------------------------------------------------
//
// suppose its_Current_DefProb has been correctly set to its_nbnames
//

void	ICM_Gaussian_LossDistrib_2F	::	ComputeConditionalDefaultProbabilities(double& value_factor, double& value_sector)
{
	int	j;

	// ----------------------------------
	double	pp;
	double	TheRatio;
	// ----------------------------------

	if (itsIntegrationMethod_Common_Factor == qGAUSS_HERMITE)
		value_factor	*=	SQRT2;

	if (itsIntegrationMethod_Sector_Factor == qGAUSS_HERMITE)
		value_sector	*=	SQRT2;

	for (j=0; j<its_nbnames; j++)
	{
		pp	=	its_coeff_a[j];
		pp	-=	its_coeff_b[j] * value_factor;
		pp	-=	its_coeff_c[j] * value_sector;

		TheRatio	=	its_coeff_den[j];

		if (! CHECK_EQUAL(fabs(TheRatio), 0.0))
		{
			//	Normal Cumulative Function	// NAG no error nor warnings
			its_Current_DefProb[j] = NAG_cumul_normal(pp); // / TheRatio);
		}
		else
			if (pp > 0.0)	// + infinite => N(+infinite)
				its_Current_DefProb[j] = 1.0;
			else			// - infinite => N(-infinite)
				its_Current_DefProb[j] = 0.0;
	}
}


// ------------------------------------------------------------------------------------------
// COMPUTE EXPECTED LOSS TRANCHE
// ------------------------------------------------------------------------------------------
//
// Algorithm. cf. 'All your hedges in one basket' 
//

double	ICM_Gaussian_LossDistrib_2F	::	ComputeExpectedLossTranche(const double& LossMin, const double& LossMax)
{
	int		l;
	int		lup, ldown;
	double	loss_level;

	double	exp_loss_tranche_up;

	ldown	=	floor(LossMin / its_lossunit);
	lup		=	floor(LossMax / its_lossunit);
	
	// Resize
	its_ProbCond->ResizeWithCopy(lup+1, 0.);
	its_lossdistrib.resize(lup+1);

	// -------------------------------------------------------
	// Compute of the Loss Distribution for this given LossUp
	// --> gets its_lossdistrib[] and its_taildistrib
	// -------------------------------------------------------
	LossProbabilityDistribution(lup);

	exp_loss_tranche_up	=	0.0;

	// loop, loss units
	for (l=ldown+1; l<=lup; l++)
	{
		loss_level	=	l * its_lossunit;
		exp_loss_tranche_up		+=	its_lossdistrib[l] * (loss_level - LossMin);
	}

	// add the tail
	exp_loss_tranche_up		+=	(LossMax - LossMin) * its_taildistrib;
	
	return	exp_loss_tranche_up;
}


// ------------------------------------------------------------------------------------------
// COMPUTES LOSS PROBABILITY DISTRIBUTION until a given loss level, keep the tail
// ------------------------------------------------------------------------------------------

void ICM_Gaussian_LossDistrib_2F::LossProbabilityDistribution(const int& lup)
{
	int		l;
	double	cumul_distrib;
		
	// START ALGORITHM, NO LOSS
	its_index_loss	=	0;		// global variable
	
	its_lossdistrib[its_index_loss]	=	Compute_Conditional_Distribution_Loss_Zero();
	
	// THEN RECURSION
	for (its_index_loss=1; its_index_loss<=lup; its_index_loss++)
		its_lossdistrib[its_index_loss]	=	Compute_Conditional_Distribution_Losses();

	// GET the TAIL
	cumul_distrib	=	0.0;
	for (l=0; l<=lup; l++)
		cumul_distrib	+=	its_lossdistrib[l];

	its_taildistrib	=	1.0 - cumul_distrib;
}


// ------------------------------------------------------------------------------------------
// START THE COMPUTATION OF LOSS DISTRIBUTION
// ------------------------------------------------------------------------------------------

double	ICM_Gaussian_LossDistrib_2F	::	Compute_Conditional_Distribution_Loss_Zero()
{
	double pk;

	int		k;
	double	CondProbaZeroLoss;
	
	CondProbaZeroLoss	=	1.0;

	int	its_index_z;

	its_index_z	=	its_factor_state * its_IntegrationStep_Common_Factor + its_sector_state;

	// its_With no elements in the basket
	// its_index_loss = 0
	its_ProbCond->SetElt(0, its_index_z, 0, CondProbaZeroLoss);

	for (k=1; k<=its_nbnames; k++) 
	{
		// Conditional default probability for Credit k
		pk	=	its_Current_DefProb[k-1];

		// loss level l = 0.0
		// no credit defaults!

		CondProbaZeroLoss	*=	(1-pk);

		// its_index_loss = 0
		its_ProbCond->SetElt(k, its_index_z, 0, CondProbaZeroLoss);
	}

	return CondProbaZeroLoss;
}


double ICM_Gaussian_LossDistrib_2F	::	Compute_Conditional_Distribution_Losses()
{
	int	k;
	
	// for the current level of loss
	double	loss_Basket_size_k; // proba de loss=L pour un panier de taille k
	
	double	tmp2; // proba de loss=L-1 pour un panier de taille k
	double	pk;     // proba cond de défaut du nom k
	int		ind_lastloss;

	// start its_With Loss = 0.0
	ind_lastloss	=	0;
	int	its_index_z;

	its_index_z	=	its_factor_state * its_IntegrationStep_Common_Factor + its_sector_state;

	its_ProbCond->SetElt(0, its_index_z, its_index_loss, 0.);

	// its_index_Loss: Loss Distribution for a basket of size (0) Knowing factor value its_index_z
	loss_Basket_size_k	=	its_ProbCond->Elt(0, its_index_z, its_index_loss);

	if (IsHomogeneous())
	{
		for (k=1; k<=its_nbnames; k++) 
		{
			// Conditional default probability for Credit k
			pk	=	its_Current_DefProb[k-1];

			// its_index_Loss: Loss Distribution for a basket of size (k-1) Knowing factor value its_index_z
			// tmp1

			// Loss Distribution for a basket of size k Knowing factor value its_index_z
			tmp2	=	loss_Basket_size_k * (1-pk) + pk * its_ProbCond->Elt(k-1, its_index_z, its_index_loss-1);

			its_ProbCond->SetElt(k, its_index_z, its_index_loss, tmp2);

			loss_Basket_size_k = tmp2;
		}
	}
	else
	{
		for (k=1; k<=its_nbnames; k++) 
		{
			// current level of loss - loss rate of names k
			ind_lastloss = its_index_loss - its_int_lossrates[k-1];

			// Conditional default probability for Credit k
			pk	=	its_Current_DefProb[k-1];

			// its_index_Loss: Loss Distribution for a basket of size (k-1) Knowing factor value its_index_z
			// tmp1

			// Loss Distribution for a basket of size k Knowing factor value its_index_z
			tmp2	=	loss_Basket_size_k * (1-pk);
			
			if (ind_lastloss >= 0)
				tmp2	+=	pk * its_ProbCond->Elt(k-1, its_index_z, ind_lastloss);
			// else: the loss of name k is too big

			its_ProbCond->SetElt(k,its_index_z, its_index_loss, tmp2);

			loss_Basket_size_k = tmp2;
		}
	}

	return tmp2;
}


void ICM_Gaussian_LossDistrib_2F :: SetIntegratorType(qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep)
{
	if ((TheStep % 2 == 0) && (TheIntegratorType == qGAUSS_LEGENDRE))
		TheIntegratorType = qGAUSS_HERMITE;
	else if ((TheStep % 2 == 1) && (TheIntegratorType == qGAUSS_HERMITE))
		TheIntegratorType = qGAUSS_LEGENDRE;

}

// --------------------------------------------------------------------------------------
// HEDGES
// --------------------------------------------------------------------------------------
	

//----------------------------------------------//
//					View	
//----------------------------------------------//
void ICM_Gaussian_LossDistrib_2F::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Loss Dsitribution Gaussian 2 Factors ----------------- \n");
		
	ICM_Distribution::View(id, fOut);

	// TO DO...

	if (ficOut == NULL)
	{
		fclose(fOut);
	}
}



void ICM_Gaussian_LossDistrib_2F :: Compute_Conditional_Probabilities_Coeffs()
{
	int		k;
	int		size;
	int		sector_id;
	
	double	den, one_over_den;
	double	beta_value;
	double	lambda_value;
	double	barrier_value;

	double	rho;
	double	gamma;

	// -----------------------------------------------
	size	=	its_nbnames;

	// -----------------------------------------------
	switch (its_CorrelationType)
	{
		case TFCT_FULL:

			for (k=0; k<size; k++)
			{
				// Sector Id of Credit k
				sector_id	=	GetSectorMemberShip(k);

				beta_value		=	its_Betas[sector_id];
				lambda_value	=	its_Lambdas[sector_id];
				barrier_value	=	its_barrier[k];

				if (CHECK_EQUAL(fabs(beta_value), 1.0))
					ICMTHROW(ERR_INVALID_DATA,"Beta values must be strictly between -1.0 and 1.0 - precompute_coeffs method in ICM_Gaussian_LossDistrib_2F class!");
				
				den = sqrt(1.0 - beta_value * beta_value);
				one_over_den = 1.0 / den;
				
				its_coeff_a[k]	=	one_over_den;
				its_coeff_b[k]	=	one_over_den;
				its_coeff_c[k]	=	one_over_den;

				its_coeff_den[k]	=	den;
				
				its_coeff_a[k]	*=	barrier_value;
				its_coeff_b[k]	*=	beta_value * lambda_value;
				its_coeff_c[k]	*=	beta_value * sqrt(1.0 - lambda_value * lambda_value);
			}

			break;

		case TFCT_SAME_INTER_DIFF_INTRA:

			gamma	=	its_Single_inter_sector_correlation;

			for (k=0; k<size; k++)
			{
				// Sector Id of Credit k
				sector_id	=	GetSectorMemberShip(k);

				beta_value		=	its_Betas[sector_id];

				if (CHECK_EQUAL(fabs(beta_value), 0.0))
					ICMTHROW(ERR_INVALID_DATA,"Beta values must not be null - precompute_coeffs method in ICM_Gaussian_LossDistrib_2F class!");
				
				lambda_value	=	sqrt(gamma) / beta_value;
				barrier_value	=	its_barrier[k];

				if (CHECK_EQUAL(fabs(beta_value), 1.0))
					ICMTHROW(ERR_INVALID_DATA,"Beta values must be strictly between -1.0 and 1.0 - precompute_coeffs method in ICM_Gaussian_LossDistrib_2F class!");
				
				den = sqrt(1.0 - beta_value * beta_value);
				one_over_den = 1.0 / den;
				
				its_coeff_a[k]	=	one_over_den;
				its_coeff_b[k]	=	one_over_den;
				its_coeff_c[k]	=	one_over_den;
				
				its_coeff_den[k]	=	den;

				its_coeff_a[k]	*=	barrier_value;
				its_coeff_b[k]	*=	beta_value * lambda_value;
				its_coeff_c[k]	*=	beta_value * sqrt(1.0 - lambda_value * lambda_value);
			}

			break;

		case TFCT_SAME_INTER_SAME_INTRA:

			rho		=	its_Single_intra_sector_correlation;
			gamma	=	its_Single_inter_sector_correlation;

			if (rho <= 0.0)
				ICMTHROW(ERR_INVALID_DATA,"Rho must not be null or negative - precompute_coeffs method in ICM_Gaussian_LossDistrib_2F class!");

			beta_value	=	sqrt(rho);
			if (CHECK_EQUAL(fabs(beta_value), 1.0))
				ICMTHROW(ERR_INVALID_DATA,"Beta values must be < 1.0 - precompute_coeffs method in ICM_Gaussian_LossDistrib_2F class!");

			lambda_value	=	sqrt(gamma) / beta_value;

			den = sqrt(1.0 - beta_value * beta_value);
			one_over_den = 1.0 / den;

			for (k=0; k<size; k++)
			{				
				barrier_value	=	its_barrier[k];

				its_coeff_a[k]	=	one_over_den;
				its_coeff_b[k]	=	one_over_den;
				its_coeff_c[k]	=	one_over_den;
				
				its_coeff_den[k]	=	den;

				its_coeff_a[k]	*=	barrier_value;
				its_coeff_b[k]	*=	beta_value * lambda_value;
				its_coeff_c[k]	*=	beta_value * sqrt(1.0 - lambda_value * lambda_value);
			}

			break;
		}
}


// if they have been computed once, then I just I have to update the barrier level
// as I have got no Beta Term Structure right now

void ICM_Gaussian_LossDistrib_2F :: Update_Conditional_Probabilities_Coeffs()
{
	int		k;
	int		size;
	
	double	den;
	double	barrier_value;

	// -----------------------------------------------
	size	=	its_nbnames;
	// -----------------------------------------------

	for (k=0; k<size; k++)
	{
		barrier_value	=	its_barrier[k];
		
		den	=	its_coeff_den[k];

		its_coeff_a[k]	=	barrier_value / den;
	}

}

