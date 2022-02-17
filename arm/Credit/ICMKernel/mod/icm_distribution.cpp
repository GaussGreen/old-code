#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_distribution.h"
#include "ICMKernel\glob\icm_maths.h"


CtxtDistrib::~CtxtDistrib() 
{
	if (m_ProbaCond_Up) delete m_ProbaCond_Up;
	m_ProbaCond_Up = NULL;
	if (m_ProbaCond_Down) delete m_ProbaCond_Down;
	m_ProbaCond_Down = NULL;
}

void 
CtxtDistrib::Init()
{
	m_IsUsed=false;
	m_IsHomogeneous=false;
	m_IsBegin=false;

	m_nbnames=0;
	m_LossRates.Resize(0);
	m_discretizationstep=40;
	m_CopulaType=qGAUSSIAN;
	m_IntegrationMethod=qGAUSS_HERMITE;
	m_IntegrationStep=40;
	m_LossUnit=0.;
	m_DistribType=qDISTRIB_STD;

	m_TrancheDown=0.;
	m_TrancheUp=0.;
	m_TrancheDown_stepup=0.;
	m_TrancheUp_stepup=0.;

	m_TailProbaLosseUp=0.;
	m_TailProbaLosseDown=0.;

	m_IsLongShort=false;

	m_ProbaCond_Up=NULL;
	m_ProbaCond_Down=NULL;

	//standard distrib parameters
	m_pdef_maturity.Resize(0);
	m_beta_up_maturity.Resize(0);
	m_beta_dw_maturity.Resize(0);

	//collat
	m_collat_pdef_start.Resize(0);

	//term structure review
	m_ts_pdef_up_maturity_t1.Resize(0);
	m_ts_pdef_up_maturity_t2.Resize(0);

	m_ts_pdef_dw_maturity_t1.Resize(0);
	m_ts_pdef_dw_maturity_t2.Resize(0);

	m_ts_t1_corr=0.;
	m_ts_t2_corr=0.;

	m_ts_beta_up_maturity_t1.Resize(0);
	m_ts_beta_up_maturity_t2.Resize(0);

	m_ts_beta_dw_maturity_t1.Resize(0);
	m_ts_beta_dw_maturity_t2.Resize(0);

	//collateral forward + ts

	m_collat_beta_up_start.Resize(0);
	m_collat_beta_dw_start.Resize(0);

	m_ts_collat_pdef_up_start_t1.Resize(0);
	m_ts_collat_pdef_up_start_t2.Resize(0);

	m_ts_collat_pdef_dw_start_t1.Resize(0);
	m_ts_collat_pdef_dw_start_t2.Resize(0);

	m_ts_collat_beta_up_start_t1.clear();
	m_ts_collat_beta_up_start_t2.clear();

	m_ts_collat_beta_dw_start_t1.Resize(0);
	m_ts_collat_beta_dw_start_t2.Resize(0);

	m_ts_collat_start_t1_corr=0.;
	m_ts_collat_start_t2_corr=0.;

	m_losses_reset_beg=0.;
	m_losses_reset_end=0.;

	m_barriers_maturity.Resize(0);
	m_collat_barriers_start.Resize(0);
	m_ts_barriers_up_maturity_t1.Resize(0);
	m_ts_barriers_up_maturity_t2.Resize(0);
	m_ts_barriers_dw_maturity_t1.Resize(0);
	m_ts_barriers_dw_maturity_t2.Resize(0);
	m_ts_collat_barriers_up_start_t1.Resize(0);
	m_ts_collat_barriers_up_start_t2.Resize(0);
	m_ts_collat_barriers_dw_start_t1.Resize(0);
	m_ts_collat_barriers_dw_start_t2.Resize(0);

	//coefs
	m_a_up_maturity.clear();
	m_b_up_maturity.clear();
	m_a_dw_maturity.clear();
	m_b_dw_maturity.clear();
	m_collat_a_up_start.clear();
	m_collat_b_up_start.clear();
	m_collat_a_dw_start.clear();
	m_collat_b_dw_start.clear();
	m_ts_a_up_maturity_t1.clear();
	m_ts_b_up_maturity_t1.clear();
	m_ts_a_up_maturity_t2.clear();
	m_ts_b_up_maturity_t2.clear();
	m_ts_a_dw_maturity_t1.clear();
	m_ts_b_dw_maturity_t1.clear();
	m_ts_a_dw_maturity_t2.clear();
	m_ts_b_dw_maturity_t2.clear();
	m_ts_collat_a_up_start_t1.clear();
	m_ts_collat_b_up_start_t1.clear();
	m_ts_collat_a_dw_start_t1.clear();
	m_ts_collat_b_dw_start_t1.clear();
	m_ts_collat_a_up_start_t2.clear();
	m_ts_collat_b_up_start_t2.clear();
	m_ts_collat_a_dw_start_t2.clear();
	m_ts_collat_b_dw_start_t2.clear();

	m_ts_stepup_pdef_up_maturity_t1.Resize(0);
	m_ts_stepup_pdef_up_maturity_t2.Resize(0);
	m_ts_stepup_pdef_dw_maturity_t1.Resize(0);
	m_ts_stepup_pdef_dw_maturity_t2.Resize(0);
	m_ts_stepup_beta_up_maturity_t1.clear();
	m_ts_stepup_beta_up_maturity_t2.clear();
	m_ts_stepup_barriers_up_maturity_t1.clear();
	m_ts_stepup_barriers_up_maturity_t2.clear();
	m_ts_stepup_beta_dw_maturity_t1.clear();
	m_ts_stepup_beta_dw_maturity_t2.clear();
	m_ts_stepup_barriers_dw_maturity_t1.clear();
	m_ts_stepup_barriers_dw_maturity_t2.clear();
	m_ts_stepup_a_up_maturity_t1.clear();
	m_ts_stepup_b_up_maturity_t1.clear();
	m_ts_stepup_a_up_maturity_t2.clear();
	m_ts_stepup_b_up_maturity_t2.clear();
	m_ts_stepup_a_dw_maturity_t1.clear();
	m_ts_stepup_b_dw_maturity_t1.clear();
	m_ts_stepup_a_dw_maturity_t2.clear();
	m_ts_stepup_b_dw_maturity_t2.clear();
	m_ts_stepup_t1_corr=0.;
	m_ts_stepup_t2_corr=0.;

	m_ProbasLossesUp.clear();
	m_ProbasLossesDown.clear();

	m_vn_contexts.clear();
	m_YearTerm=0.;
	m_Recoveries.clear();
	m_Notionals.clear();
	m_TotalNotionals=0.;
}

CtxtDistrib& 
CtxtDistrib::operator= (const CtxtDistrib& ref) 
{
	if (this!=&ref) 
	{
		this->~CtxtDistrib(); 
		new(this)CtxtDistrib(ref); 
	}
	return *this; 
}

CtxtDistrib::CtxtDistrib(const CtxtDistrib& i) 
{
	Init();

	m_IsUsed=i.m_IsUsed;
	m_IsHomogeneous=i.m_IsHomogeneous;
	m_IsBegin=i.m_IsBegin;

	m_nbnames=i.m_nbnames;
	m_LossUnit=i.m_LossUnit;
	m_DistribType=i.m_DistribType;
	m_LossRates=i.m_LossRates;
	m_discretizationstep=i.m_discretizationstep;
	m_CopulaType=i.m_CopulaType;
	m_IntegrationMethod=i.m_IntegrationMethod;
	m_IntegrationStep=i.m_IntegrationStep;

	m_TrancheDown=i.m_TrancheDown;
	m_TrancheUp=i.m_TrancheUp;
	m_TrancheDown_stepup=i.m_TrancheDown_stepup;
	m_TrancheUp_stepup=i.m_TrancheUp_stepup;

	m_OutLosses=i.m_OutLosses;
	m_TailProbaLosseUp=i.m_TailProbaLosseUp;
	m_TailProbaLosseDown=i.m_TailProbaLosseDown;

	m_IsLongShort = i.m_IsLongShort;

	//matrices stockage probas conditionelles
	m_ProbaCond_Up=i.m_ProbaCond_Up;
	m_ProbaCond_Down=i.m_ProbaCond_Down;

	//standard distrib parameters
	m_pdef_maturity=i.m_pdef_maturity;

	m_beta_up_maturity=i.m_beta_up_maturity;
	m_beta_dw_maturity=i.m_beta_dw_maturity;

	//collat
	m_collat_pdef_start=i.m_collat_pdef_start;

	m_collat_beta_up_start=i.m_collat_beta_up_start;
	m_collat_beta_dw_start=i.m_collat_beta_dw_start;

	//term structure review
	m_ts_pdef_up_maturity_t1=i.m_ts_pdef_up_maturity_t1;
	m_ts_pdef_up_maturity_t2=i.m_ts_pdef_up_maturity_t2;

	m_ts_pdef_dw_maturity_t1=i.m_ts_pdef_dw_maturity_t1;
	m_ts_pdef_dw_maturity_t2=i.m_ts_pdef_dw_maturity_t2;

	m_ts_beta_up_maturity_t1=i.m_ts_beta_up_maturity_t1;
	m_ts_beta_up_maturity_t2=i.m_ts_beta_up_maturity_t2;

	m_ts_beta_dw_maturity_t1=i.m_ts_beta_dw_maturity_t1;
	m_ts_beta_dw_maturity_t2=i.m_ts_beta_dw_maturity_t2;

	m_ts_t1_corr=i.m_ts_t1_corr;
	m_ts_t2_corr=i.m_ts_t2_corr;

	//collateral forward + ts
	m_ts_collat_pdef_up_start_t1=i.m_ts_collat_pdef_up_start_t1;
	m_ts_collat_pdef_up_start_t2=i.m_ts_collat_pdef_up_start_t2;

	m_ts_collat_pdef_dw_start_t1=i.m_ts_collat_pdef_dw_start_t1;
	m_ts_collat_pdef_dw_start_t2=i.m_ts_collat_pdef_dw_start_t2;

	m_ts_collat_beta_up_start_t1=i.m_ts_collat_beta_up_start_t1;
	m_ts_collat_beta_up_start_t2=i.m_ts_collat_beta_up_start_t2;

	m_ts_collat_beta_dw_start_t1=i.m_ts_collat_beta_dw_start_t1;
	m_ts_collat_beta_dw_start_t2=i.m_ts_collat_beta_dw_start_t2;

	m_ts_collat_start_t1_corr=i.m_ts_collat_start_t1_corr;
	m_ts_collat_start_t2_corr=i.m_ts_collat_start_t2_corr;

	//step up
	m_losses_reset_beg=i.m_losses_reset_beg;
	m_losses_reset_end=i.m_losses_reset_end;

	//barriers
	m_barriers_maturity=i.m_barriers_maturity;
	m_collat_barriers_start=i.m_collat_barriers_start;
	m_ts_barriers_up_maturity_t1=i.m_ts_barriers_up_maturity_t1;
	m_ts_barriers_up_maturity_t2=i.m_ts_barriers_up_maturity_t2;
	m_ts_barriers_dw_maturity_t1=i.m_ts_barriers_dw_maturity_t1;
	m_ts_barriers_dw_maturity_t2=i.m_ts_barriers_dw_maturity_t2;
	m_ts_collat_barriers_up_start_t1=i.m_ts_collat_barriers_up_start_t1;
	m_ts_collat_barriers_up_start_t2=i.m_ts_collat_barriers_up_start_t2;
	m_ts_collat_barriers_dw_start_t1=i.m_ts_collat_barriers_dw_start_t1;
	m_ts_collat_barriers_dw_start_t2=i.m_ts_collat_barriers_dw_start_t2;

	//coefs
	m_a_up_maturity=i.m_a_up_maturity;
	m_b_up_maturity=i.m_b_up_maturity;
	m_a_dw_maturity=i.m_a_dw_maturity;
	m_b_dw_maturity=i.m_b_dw_maturity;
	m_collat_a_up_start=i.m_collat_a_up_start;
	m_collat_b_up_start=i.m_collat_b_up_start;
	m_collat_a_dw_start=i.m_collat_a_dw_start;
	m_collat_b_dw_start=i.m_collat_b_dw_start;
	m_ts_a_up_maturity_t1=i.m_ts_a_up_maturity_t1;
	m_ts_b_up_maturity_t1=i.m_ts_b_up_maturity_t1;
	m_ts_a_up_maturity_t2=i.m_ts_a_up_maturity_t2;
	m_ts_b_up_maturity_t2=i.m_ts_b_up_maturity_t2;
	m_ts_a_dw_maturity_t1=i.m_ts_a_dw_maturity_t1;
	m_ts_b_dw_maturity_t1=i.m_ts_b_dw_maturity_t1;
	m_ts_a_dw_maturity_t2=i.m_ts_a_dw_maturity_t2;
	m_ts_b_dw_maturity_t2=i.m_ts_b_dw_maturity_t2;
	m_ts_collat_a_up_start_t1=i.m_ts_collat_a_up_start_t1;
	m_ts_collat_b_up_start_t1=i.m_ts_collat_b_up_start_t1;
	m_ts_collat_a_dw_start_t1=i.m_ts_collat_a_dw_start_t1;
	m_ts_collat_b_dw_start_t1=i.m_ts_collat_b_dw_start_t1;
	m_ts_collat_a_up_start_t2=i.m_ts_collat_a_up_start_t2;
	m_ts_collat_b_up_start_t2=i.m_ts_collat_b_up_start_t2;
	m_ts_collat_a_dw_start_t2=i.m_ts_collat_a_dw_start_t2;
	m_ts_collat_b_dw_start_t2=i.m_ts_collat_b_dw_start_t2;

	//step up term structure
	m_ts_stepup_pdef_up_maturity_t1=i.m_ts_stepup_pdef_up_maturity_t1;
	m_ts_stepup_pdef_up_maturity_t2=i.m_ts_stepup_pdef_up_maturity_t2;
	m_ts_stepup_pdef_dw_maturity_t1=i.m_ts_stepup_pdef_dw_maturity_t1;
	m_ts_stepup_pdef_dw_maturity_t2=i.m_ts_stepup_pdef_dw_maturity_t2;
	m_ts_stepup_beta_up_maturity_t1=i.m_ts_stepup_beta_up_maturity_t1;
	m_ts_stepup_beta_up_maturity_t2=i.m_ts_stepup_beta_up_maturity_t2;
	m_ts_stepup_barriers_up_maturity_t1=i.m_ts_stepup_barriers_up_maturity_t1;
	m_ts_stepup_barriers_up_maturity_t2=i.m_ts_stepup_barriers_up_maturity_t2;
	m_ts_stepup_beta_dw_maturity_t1=i.m_ts_stepup_beta_dw_maturity_t1;
	m_ts_stepup_beta_dw_maturity_t2=i.m_ts_stepup_beta_dw_maturity_t2;
	m_ts_stepup_barriers_dw_maturity_t1=i.m_ts_stepup_barriers_dw_maturity_t1;
	m_ts_stepup_barriers_dw_maturity_t2=i.m_ts_stepup_barriers_dw_maturity_t2;
	m_ts_stepup_a_up_maturity_t1=i.m_ts_stepup_a_up_maturity_t1;
	m_ts_stepup_b_up_maturity_t1=i.m_ts_stepup_b_up_maturity_t1;
	m_ts_stepup_a_up_maturity_t2=i.m_ts_stepup_a_up_maturity_t2;
	m_ts_stepup_b_up_maturity_t2=i.m_ts_stepup_b_up_maturity_t2;
	m_ts_stepup_a_dw_maturity_t1=i.m_ts_stepup_a_dw_maturity_t1;
	m_ts_stepup_b_dw_maturity_t1=i.m_ts_stepup_b_dw_maturity_t1;
	m_ts_stepup_a_dw_maturity_t2=i.m_ts_stepup_a_dw_maturity_t2;
	m_ts_stepup_b_dw_maturity_t2=i.m_ts_stepup_b_dw_maturity_t2;
	m_ts_stepup_t1_corr=i.m_ts_stepup_t1_corr;
	m_ts_stepup_t2_corr=i.m_ts_stepup_t2_corr;

	m_vn_contexts = i.m_vn_contexts;
	m_YearTerm = i.m_YearTerm;

	m_Recoveries = i.m_Recoveries;
	m_Notionals = i.m_Notionals;
	m_TotalNotionals = i.m_TotalNotionals;

}
void 
CtxtDistrib::ComputeAll()
{
	if (m_vn_contexts.size()>0)
	{	for (int i=0;i<m_vn_contexts.size();i++)
		{ m_vn_contexts[i].ComputeAll();}}

	IsHomogeneous();
	ComputeALLBarriers();
	ComputeALLCoefs();
}


bool operator > (const CtxtDistrib& d1, const CtxtDistrib& d2)
{    
    int    i, j;
	double result = true;

	if (
		(d2.m_YearTerm<=d1.m_YearTerm) ||
		(d2.m_LossUnit>d1.m_LossUnit)
		)
	{result=false;}

	if (d1.m_LossRates.size() != d2.m_LossRates.size())
	{result=false;}

	if (result)
	for (i=0;i<d1.m_LossRates.size();i++)
	{
		if (d2.m_LossRates[i]>d1.m_LossRates[i])
		{result = false;break;}

		if (d2.m_pdef_maturity[i]<=d1.m_pdef_maturity[i])
		{result = false;break;}

/*
		if (d2.m_barriers_maturity[i]<=d1.m_barriers_maturity[i])
		{result = false;break;}

		if (d2.m_a_up_maturity[i]<=d1.m_a_up_maturity[i])
		{result = false;break;}

		if (d2.m_b_up_maturity[i]<=d1.m_b_up_maturity[i])
		{result = false;break;}

		if (d2.m_a_dw_maturity[i]<=d1.m_a_dw_maturity[i])
		{result = false;break;}

		if (d2.m_b_dw_maturity[i]<=d1.m_b_dw_maturity[i])
		{result = false;break;}
*/
	}


	return (result);
}

bool ICM_Distribution::itsFwdCollatTScase=false;

void 
ICM_Distribution::Set(const int& nbnames)
{
	SetNbNames(nbnames);

	if (its_ProbCond) delete its_ProbCond;
	its_ProbCond= new ICM_QCubix<double>(nbnames+1,20,1,0.);

	// 17783 if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
	// 17783 its_ProbCond_Perturb= new ICM_QCubix<double>(nbnames+1,20,1,0.);

	// 17783 if (its_lossdistrib_perturb) delete its_lossdistrib_perturb;
	// 17783 its_lossdistrib_perturb = new ICM_QMatrix<double>(1,nbnames);

	its_int_lossrates.resize(nbnames);
	its_dbl_lossrates.resize(nbnames);
	its_lossdistrib.resize(1);
	// 17783 its_taildistrib_perturb.resize(nbnames);

	//Long Short
	/*********************************************************************/

	its_sorted_indices.resize(nbnames);

	/*********************************************************************/

}

// virtual 
ICM_Distribution::~ICM_Distribution()
{
	if (its_ProbCond)
		delete its_ProbCond;
	its_ProbCond = NULL;

	// 17783 if (its_ProbCond_Perturb)
	// 17783 	delete its_ProbCond_Perturb;
	// 17783 its_ProbCond_Perturb = NULL;

	// 17783 if (its_lossdistrib_perturb)
	// 17783 	delete its_lossdistrib_perturb;
	// 17783 its_lossdistrib_perturb = NULL;

	// 17783 if (itsShifts)
	// 17783 	delete itsShifts;
	// 17783 itsShifts = NULL;

	its_unique_beta.clear();
	its_pdef_at_maturity.clear();
	its_pdef_perturb.clear();
	its_barrier.clear();
	its_barrier_perturb.clear();
	its_lossdistrib.clear();
	// 17783 its_taildistrib_perturb.clear();

	//Long Short
	/*********************************************************************/
	its_sorted_indices.clear();
	/*********************************************************************/

	//collateral forward
	its_collat_fwd_unique_beta.clear();
	itsFwdCollatTScase = false;

	its_unique_beta_t2.clear();
	its_pdef_at_maturity_t2.clear();
	its_barrier_t2.clear();
	its_collat_fwd_unique_beta_t2.clear();
	its_pdef_perturb_t2.clear();
	its_barrier_perturb_t2.clear();

	if (its_ProbCond_stepup)
		delete its_ProbCond_stepup;
	its_ProbCond_stepup = NULL;

	its_ind_x_stepup=-1;					
};

//----- Utilities

void 
ICM_Distribution::Init()
{
	its_ctxt=NULL;
	its_ishomogeneous=true;

	its_ProbCond=NULL;
	// 17783 its_ProbCond_Perturb=NULL;
	// 17783 its_lossdistrib_perturb=NULL;
	// 17783 itsShifts = NULL;

	its_unique_beta.clear(); 
	its_pdef_at_maturity.clear();
	its_pdef_perturb.clear();
	its_barrier.clear();
	its_barrier_perturb.clear();
	its_lossdistrib.clear();
	// 17783 its_taildistrib_perturb.clear();
	its_int_lossrates.clear();
	its_dbl_lossrates.clear();

	its_nbnames=0;
	its_lup=0;
	its_taildistrib=0.;
	its_ind_name=0;
	its_lossunit = 0.;

	its_ind_x=-1;
	its_ind_Loss=1;

	its_normPi = 1./sqrt(PI);
	its_2Pi = 2.*PI;
	its_sqrt2 = sqrt(2.);
	its_min_pdef = 0.;

	//LongShort
	/*********************************************************************/
	its_nbnames_short=0;
	its_LSminloss=0; 
	its_sorted_indices.clear();
	/*********************************************************************/

	//collateral forward
	itsFwdCollatTScase = false;
	its_collat_fwd_unique_beta.clear();

	its_unique_beta_t2.clear();
	its_pdef_at_maturity_t2.clear();
	its_barrier_t2.clear();
	itsTSR = false;
	its_t1=0.;
	its_t2=0.;
	its_collat_t1=0.;
	its_collat_t2=0.;
	its_collat_fwd_unique_beta_t2.clear();
	its_pdef_perturb_t2.clear();
	its_barrier_perturb_t2.clear();

	its_ProbCond_stepup=NULL;

	its_ind_x_stepup=-1;

	its_IsUp = true;
}


void ICM_Distribution::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char  fOutName[200];
	int i=0;

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

    fprintf(fOut, "\t\t\t ----------------- Distribution ----------------- \n");
	fprintf(fOut,"ind_Loss : %i\n",its_ind_Loss);
	fprintf(fOut,"ind_x : %i\n",its_ind_x);
	fprintf(fOut,"ind_name : %i\n",its_ind_name);
	fprintf(fOut,"nbnames : %i\n",its_nbnames);
	fprintf(fOut,"lup : %i\n",its_lup);
	fprintf(fOut,"taildistrib : %f\n",its_taildistrib);
	fprintf(fOut,"lossunit : %f\n",its_lossunit);
	fprintf(fOut, "\t\t\t ----------------- Beta at maturity ----------------- \n");
	for (i=0;i<its_nbnames;i++)
	{fprintf(fOut,"%f\n",its_unique_beta[i]);}
	fprintf(fOut, "\t\t\t ----------------- Default Probability at maturity ----------------- \n");
	for (i=0;i<its_nbnames;i++)
	{fprintf(fOut,"%f\n",its_pdef_at_maturity[i]);}
	fprintf(fOut, "\t\t\t ----------------- Barrier at maturity ----------------- \n");
	for (i=0;i<its_nbnames;i++)
	{fprintf(fOut,"barrier : %f\n",its_barrier[i]);}

	fprintf(fOut,"\n\n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



void ICM_Distribution::BitwiseCopy(const ARM_Object* src)
{
    ICM_Distribution* srcdistrib = (ICM_Distribution *) src;

	its_ishomogeneous = srcdistrib->its_ishomogeneous;
	its_ind_Loss = srcdistrib->its_ind_Loss;
	its_ind_x = srcdistrib->its_ind_x;
	its_ind_name = srcdistrib->its_ind_name;
	its_nbnames = srcdistrib->its_nbnames;
	its_lup = srcdistrib->its_lup;
	its_taildistrib = srcdistrib->its_taildistrib;
	its_lossunit = srcdistrib->its_lossunit;

	its_unique_beta = srcdistrib->its_unique_beta;
	its_pdef_at_maturity = srcdistrib->its_pdef_at_maturity;
	its_pdef_perturb = srcdistrib->its_pdef_perturb;
	its_barrier = srcdistrib->its_barrier;
	its_barrier_perturb = srcdistrib->its_barrier_perturb;
	its_lossdistrib = srcdistrib->its_lossdistrib;
	// 17783 its_taildistrib_perturb = srcdistrib->its_taildistrib_perturb;
	its_int_lossrates = srcdistrib->its_int_lossrates;
	its_dbl_lossrates = srcdistrib->its_dbl_lossrates;

	its_normPi = srcdistrib->its_normPi;
	its_2Pi = srcdistrib->its_2Pi;
	its_sqrt2 = srcdistrib->its_sqrt2;

	if (srcdistrib->its_ProbCond)
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond = (ICM_QCubix<double>*) srcdistrib->its_ProbCond->Clone();
	}
	
	// 17783 if (srcdistrib->its_ProbCond_Perturb)
	// 17783 {
	// 17783 	if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
	// 17783 	its_ProbCond_Perturb = (ICM_QCubix<double>*) srcdistrib->its_ProbCond_Perturb->Clone();
	// 17783 }

	// 17783 if (srcdistrib->its_lossdistrib_perturb)
	// 17783 {
	// 17783 	if (its_lossdistrib_perturb) delete its_lossdistrib_perturb;
	// 17783 	its_lossdistrib_perturb = (ICM_QMatrix<double>*) srcdistrib->its_lossdistrib_perturb->Clone();
	// 17783 }

	// 17783 if (srcdistrib->itsShifts)
	// 17783 {
	// 17783 	if (itsShifts) delete itsShifts;
	// 17783 	itsShifts = (ICM_QMatrix<double>*) srcdistrib->itsShifts->Clone();
	// 17783 }

	//LongShort
	/*********************************************************************/

	its_nbnames_short = srcdistrib->its_nbnames_short;
	its_LSminloss= srcdistrib->its_LSminloss;
	its_sorted_indices = srcdistrib->its_sorted_indices;
	
	/*********************************************************************/

	//Forward collateral
	its_collat_fwd_unique_beta = srcdistrib->its_collat_fwd_unique_beta;
}

// -------------
//	Copy Method 
// -------------
void ICM_Distribution::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);
     BitwiseCopy(src);
}

ARM_Object* ICM_Distribution::Clone(void)
{
     ICM_Distribution* theClone = new ICM_Distribution();

     theClone->Copy(this);
 
     return(theClone);
}

void
ICM_Distribution::compute_barrier()
{
	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;

	int nbnames = its_nbnames;
	its_barrier.resize(nbnames);

	if (IsTSR()==false) {
	for (int k=0;k<nbnames;k++)
	{
		//cas standard
		if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
			its_barrier[k] = Limit_case_Minus;
		else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
			its_barrier[k] = Limit_case_Plus;
		else
			its_barrier[k] = NAG_deviates_normal_dist(its_pdef_at_maturity[k]);
	}}
	else
	{	//cas TSR

		if (GetT2())
		{
			for (int k=0;k<nbnames;k++)
			{
			if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
				its_barrier[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
				its_barrier[k] = Limit_case_Plus;
			else
				its_barrier[k] = NAG_deviates_normal_dist(its_pdef_at_maturity[k]);
			}
		}
		else
		{	//cas t>T1
			for (int k=0;k<nbnames;k++)
			{
			if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
				its_barrier[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
				its_barrier[k] = Limit_case_Plus;
			else
				its_barrier[k] = NAG_deviates_normal_dist(its_pdef_at_maturity[k]);
			}
		}

		its_barrier_t2.resize(nbnames);

		for (int k=0;k<nbnames;k++)
		{
			if (fabs(its_pdef_at_maturity_t2[k]) < DB_TOL)
				its_barrier_t2[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity_t2[k]-1.0) < DB_TOL)
				its_barrier_t2[k] = Limit_case_Plus;
			else
				its_barrier_t2[k] = NAG_deviates_normal_dist(its_pdef_at_maturity_t2[k]);
		}
	}

}

void ICM_Distribution::compute_barrier_perturb()
{
	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;

	int nbnames = its_nbnames;
	its_barrier_perturb.resize(nbnames);

	if (IsTSR()==false)
	{
	for (int k=0;k<nbnames;k++)
	{
		if (fabs(its_pdef_perturb[k]) < DB_TOL)
			its_barrier_perturb[k] = Limit_case_Minus;
		else if (fabs(its_pdef_perturb[k]-1.0) < DB_TOL)
			its_barrier_perturb[k] = Limit_case_Plus;
		else
			its_barrier_perturb[k] = NAG_deviates_normal_dist(its_pdef_perturb[k]);
	}}
	else
	{
		if (GetT2())
		{
			for (int k=0;k<nbnames;k++)
			{	
			if (fabs(its_pdef_perturb[k]) < DB_TOL)
				its_barrier_perturb[k] = Limit_case_Minus;
			else if (fabs(its_pdef_perturb[k]-1.0) < DB_TOL)
				its_barrier_perturb[k] = Limit_case_Plus;
			else
				its_barrier_perturb[k] = NAG_deviates_normal_dist(its_pdef_perturb[k]);
			}
		}else
		{
			for (int k=0;k<nbnames;k++)
			{	
			if (fabs(its_pdef_perturb[k]) < DB_TOL)
				its_barrier_perturb[k] = Limit_case_Minus;
			else if (fabs(its_pdef_perturb[k]-1.0) < DB_TOL)
				its_barrier_perturb[k] = Limit_case_Plus;
			else
				its_barrier_perturb[k] = NAG_deviates_normal_dist(its_pdef_perturb[k]);
			}
		}

		its_barrier_perturb_t2.resize(nbnames);

		for (int k=0;k<nbnames;k++)
		{
		if (fabs(its_pdef_perturb_t2[k]) < DB_TOL)
			its_barrier_perturb_t2[k] = Limit_case_Minus;
		else if (fabs(its_pdef_perturb_t2[k]-1.0) < DB_TOL)
			its_barrier_perturb_t2[k] = Limit_case_Plus;
		else
			its_barrier_perturb_t2[k] = NAG_deviates_normal_dist(its_pdef_perturb_t2[k]);
		}
	}

}


// ----------------------------------------------------------------------
// Compute barriers
// --------------f--------------------------------------------------------
void 
CtxtDistrib::ComputeBarriers(int& size,const ARM_Vector& input,ARM_Vector & output)
{
	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;

	output.Resize(size);

	for (int k=0;k<input.size();k++)
	{
		if (fabs(input[k]) < DB_TOL)
			output[k] = Limit_case_Minus;
		else if (fabs(input[k]-1.0) < DB_TOL)
			output[k] = Limit_case_Plus;
		else
			output[k] = NAG_deviates_normal_dist(input[k]);
	}
}

void CtxtDistrib::ComputeALLBarriers()
{
	if (m_IsUsed==false) {return;}


	if (!m_barriers_maturity.size() && m_pdef_maturity.size()) ComputeBarriers(m_nbnames,m_pdef_maturity,m_barriers_maturity);
	if (!m_collat_barriers_start.size() && m_collat_pdef_start.size()) ComputeBarriers(m_nbnames,m_collat_pdef_start,m_collat_barriers_start);
	if (!m_ts_barriers_up_maturity_t1.size() && m_ts_pdef_up_maturity_t1.size()) ComputeBarriers(m_nbnames,m_ts_pdef_up_maturity_t1,m_ts_barriers_up_maturity_t1);
	if (!m_ts_barriers_up_maturity_t2.size() && m_ts_pdef_up_maturity_t2.size()) ComputeBarriers(m_nbnames,m_ts_pdef_up_maturity_t2,m_ts_barriers_up_maturity_t2);
	if (!m_ts_barriers_dw_maturity_t1.size() && m_ts_pdef_dw_maturity_t1.size()) ComputeBarriers(m_nbnames,m_ts_pdef_dw_maturity_t1,m_ts_barriers_dw_maturity_t1);
	if (!m_ts_barriers_dw_maturity_t2.size() && m_ts_pdef_dw_maturity_t2.size()) ComputeBarriers(m_nbnames,m_ts_pdef_dw_maturity_t2,m_ts_barriers_dw_maturity_t2);
	if (!m_ts_collat_barriers_up_start_t1.size() && m_ts_collat_pdef_up_start_t1.size()) ComputeBarriers(m_nbnames,m_ts_collat_pdef_up_start_t1,m_ts_collat_barriers_up_start_t1);
	if (!m_ts_collat_barriers_up_start_t2.size() && m_ts_collat_pdef_up_start_t2.size()) ComputeBarriers(m_nbnames,m_ts_collat_pdef_up_start_t2,m_ts_collat_barriers_up_start_t2);
	if (!m_ts_collat_barriers_dw_start_t1.size() && m_ts_collat_pdef_dw_start_t1.size()) ComputeBarriers(m_nbnames,m_ts_collat_pdef_dw_start_t1,m_ts_collat_barriers_dw_start_t1);
	if (!m_ts_collat_barriers_dw_start_t2.size() && m_ts_collat_pdef_dw_start_t2.size()) ComputeBarriers(m_nbnames,m_ts_collat_pdef_dw_start_t2,m_ts_collat_barriers_dw_start_t2);
	if (!m_ts_stepup_barriers_up_maturity_t1.size() && m_ts_stepup_pdef_up_maturity_t1.size()) ComputeBarriers(m_nbnames,m_ts_stepup_pdef_up_maturity_t1,m_ts_stepup_barriers_up_maturity_t1);
	if (!m_ts_stepup_barriers_up_maturity_t2.size() && m_ts_stepup_pdef_up_maturity_t2.size()) ComputeBarriers(m_nbnames,m_ts_stepup_pdef_up_maturity_t2,m_ts_stepup_barriers_up_maturity_t2);
	if (!m_ts_stepup_barriers_dw_maturity_t1.size() && m_ts_stepup_pdef_dw_maturity_t1.size()) ComputeBarriers(m_nbnames,m_ts_stepup_pdef_dw_maturity_t1,m_ts_stepup_barriers_dw_maturity_t1);
	if (!m_ts_stepup_barriers_dw_maturity_t2.size() && m_ts_stepup_pdef_dw_maturity_t2.size()) ComputeBarriers(m_nbnames,m_ts_stepup_pdef_dw_maturity_t2,m_ts_stepup_barriers_dw_maturity_t2);

}

void CtxtDistrib::ComputeCoefs(const ARM_Vector& barriers,
					  const ARM_Vector& betas,
					  std::vector<double>& coefsA,
					  std::vector<double>& coefsB)
{

	if (barriers.size()==0) return;
	if (betas.size()==0) return;
//	if (coefsA.size()) return;

	double beta_value=0.;
	coefsA.resize(barriers.size());
	coefsB.resize(barriers.size());

	for (int k=0; k<barriers.size(); k++)
	{
		if (fabs(beta_value = betas[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		double den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		coefsA[k]	=	den;
		coefsA[k]	*=	barriers[k];

		coefsB[k]	=	den;
		coefsB[k]	*=	beta_value;
	}
}

void CtxtDistrib::ComputeALLCoefs()
{
	ComputeCoefs(m_barriers_maturity,m_beta_up_maturity,m_a_up_maturity, m_b_up_maturity);
	ComputeCoefs(m_barriers_maturity,m_beta_dw_maturity,m_a_dw_maturity, m_b_dw_maturity);
	ComputeCoefs(m_collat_barriers_start,m_collat_beta_up_start,m_collat_a_up_start, m_collat_b_up_start);
	ComputeCoefs(m_collat_barriers_start,m_collat_beta_dw_start,m_collat_a_dw_start, m_collat_b_dw_start);
	ComputeCoefs(m_ts_barriers_up_maturity_t1,m_ts_beta_up_maturity_t1,m_ts_a_up_maturity_t1, m_ts_b_up_maturity_t1);
	ComputeCoefs(m_ts_barriers_up_maturity_t2,m_ts_beta_up_maturity_t2,m_ts_a_up_maturity_t2, m_ts_b_up_maturity_t2);
	ComputeCoefs(m_ts_barriers_dw_maturity_t1,m_ts_beta_dw_maturity_t1,m_ts_a_dw_maturity_t1, m_ts_b_dw_maturity_t1);
	ComputeCoefs(m_ts_barriers_dw_maturity_t2,m_ts_beta_dw_maturity_t2,m_ts_a_dw_maturity_t2, m_ts_b_dw_maturity_t2);
	ComputeCoefs(m_ts_collat_barriers_up_start_t1,m_ts_collat_beta_up_start_t1,m_ts_collat_a_up_start_t1, m_ts_collat_b_up_start_t1);
	ComputeCoefs(m_ts_collat_barriers_dw_start_t1,m_ts_collat_beta_dw_start_t1,m_ts_collat_a_dw_start_t1, m_ts_collat_b_dw_start_t1);
	ComputeCoefs(m_ts_collat_barriers_up_start_t2,m_ts_collat_beta_up_start_t2,m_ts_collat_a_up_start_t2, m_ts_collat_b_up_start_t2);
	ComputeCoefs(m_ts_collat_barriers_dw_start_t2,m_ts_collat_beta_dw_start_t2,m_ts_collat_a_dw_start_t2, m_ts_collat_b_dw_start_t2);
	ComputeCoefs(m_ts_stepup_barriers_up_maturity_t1,m_ts_stepup_beta_up_maturity_t1,m_ts_stepup_a_up_maturity_t1, m_ts_stepup_b_up_maturity_t1);
	ComputeCoefs(m_ts_stepup_barriers_up_maturity_t2,m_ts_stepup_beta_up_maturity_t2,m_ts_stepup_a_up_maturity_t2, m_ts_stepup_b_up_maturity_t2);
	ComputeCoefs(m_ts_stepup_barriers_dw_maturity_t1,m_ts_stepup_beta_dw_maturity_t1,m_ts_stepup_a_dw_maturity_t1, m_ts_stepup_b_dw_maturity_t1);
	ComputeCoefs(m_ts_stepup_barriers_dw_maturity_t2,m_ts_stepup_beta_dw_maturity_t2,m_ts_stepup_a_dw_maturity_t2, m_ts_stepup_b_dw_maturity_t2);

}

/**
void CptStrikesLosses(double k1,double k2,ICM_Pricer* p,int nblosses, double& Outk1,double& Outk2)
{
	ICM_Pricer_Distrib_Smile* pricer= (ICM_Pricer_Distrib_Smile*)p;
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*)p->GetModel(); 
	ICM_Mez* ftd=(ICM_Mez*)(pricer->GetSecurity()); 
	ICM_Collateral* col = ftd->GetCollateral();

	double sum = col->SumNotionals();
	double percent_k1 = k1/sum;
	double percent_k2 = k2/sum;

	for (int i=0;i<nblosses;i++)
	{
		sum -= (1.-mod->GetRecoveryRate(col->GetIssuersLabels(i)))*
				col->GetIssuersNotional(col->GetIssuersLabels(i));
	}

	Outk1 = percent_k1 * sum;
	Outk2 = percent_k2 * sum;
}
**/ 

void CtxtDistrib::ResetAll()
{
m_barriers_maturity.clear();
m_collat_barriers_start.clear();
m_ts_barriers_up_maturity_t1.clear();
m_ts_barriers_up_maturity_t2.clear();
m_ts_barriers_dw_maturity_t1.clear();
m_ts_barriers_dw_maturity_t2.clear();
m_ts_collat_barriers_up_start_t1.clear();
m_ts_collat_barriers_up_start_t2.clear();
m_ts_collat_barriers_dw_start_t1.clear();
m_ts_collat_barriers_dw_start_t2.clear();
m_ts_stepup_barriers_up_maturity_t1.clear();
m_ts_stepup_barriers_up_maturity_t2.clear();
m_ts_stepup_barriers_dw_maturity_t1.clear();
m_ts_stepup_barriers_dw_maturity_t2.clear();

m_a_up_maturity.clear();
m_b_up_maturity.clear();
m_a_dw_maturity.clear();
m_b_dw_maturity.clear();
m_collat_a_up_start.clear();
m_collat_b_up_start.clear();
m_collat_a_dw_start.clear();
m_collat_b_dw_start.clear();
m_ts_a_up_maturity_t1.clear();
m_ts_b_up_maturity_t1.clear();
m_ts_a_up_maturity_t2.clear();
m_ts_b_up_maturity_t2.clear();
m_ts_a_dw_maturity_t1.clear();
m_ts_b_dw_maturity_t1.clear();
m_ts_a_dw_maturity_t2.clear();
m_ts_b_dw_maturity_t2.clear();
m_ts_collat_a_up_start_t1.clear();
m_ts_collat_b_up_start_t1.clear();
m_ts_collat_a_dw_start_t1.clear();
m_ts_collat_b_dw_start_t1.clear();
m_ts_collat_a_up_start_t2.clear();
m_ts_collat_b_up_start_t2.clear();
m_ts_collat_a_dw_start_t2.clear();
m_ts_collat_b_dw_start_t2.clear();
m_ts_stepup_a_up_maturity_t1.clear();
m_ts_stepup_b_up_maturity_t1.clear();
m_ts_stepup_a_up_maturity_t2.clear();
m_ts_stepup_b_up_maturity_t2.clear();
m_ts_stepup_a_dw_maturity_t1.clear();
m_ts_stepup_b_dw_maturity_t1.clear();
m_ts_stepup_a_dw_maturity_t2.clear();
m_ts_stepup_b_dw_maturity_t2.clear();
}