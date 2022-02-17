#include "firsttoinc.h"
#include "ICMKernel\pricer\icm_VanillaFwdPricer.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\mod\modelmulticurves.h"

void ICM_VanillaFwdPricer::Set(ICM_Pricer *pricer)
{
	if (pricer->GetModel()->GetName()!=ICM_MODELMULTICURVES)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this security or this model ");

	ICM_Pricer_Basket::Set(pricer->GetSecurity(), pricer->GetModel(),GetParameters(),GetAsOfDate());

	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ICM_Security* InitialSec = (ICM_Security*) pricer->GetSecurity();
	ARM_Date StartDate = InitialSec->GetStartDateNA();
	ARM_Date Maturity = InitialSec->GetEndDateNA();

	ICM_Cds* FwdStartSec= (ICM_Cds*) InitialSec->Clone();
	FwdStartSec->CptCashFlowDatesCredit(AsOf,StartDate);
	ICM_Cds* FwdEndSec= (ICM_Cds*) InitialSec->Clone();
	FwdEndSec->CptCashFlowDatesCredit(AsOf,Maturity);

	SetFwdStartSec(FwdStartSec);
	SetFwdEndSec(FwdEndSec);

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) pricer->GetModel();
	ICM_Correlation* corr = model->GetCorrelation();

	ICM_Correlation* FwdStartCorr =(ICM_Correlation*) corr->Clone();
	FwdStartCorr->ResetBetas();
	SetFwdStartCorr(FwdStartCorr);
	ICM_Correlation* FwdEndCorr = (ICM_Correlation*) corr->Clone();
	FwdEndCorr->ResetBetas();
	SetFwdEndCorr(FwdEndCorr);

	itsFwdStartPricer = pricer->CloneOnlyParams(FwdStartSec,model);
	itsFwdEndPricer = pricer->CloneOnlyParams(FwdEndSec,model);

}

// --------------------------------------------------------------------
// Compute Price
// --------------------------------------------------------------------
double ICM_VanillaFwdPricer::ComputePrice(qCMPMETH measure)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) itsFwdStartPricer->GetModel();
	ICM_Correlation* initCorr = model->GetCorrelation();
	ICM_Correlation* CpyFwdStartCorr = NULL;
	ICM_Correlation* CpyFwdEndCorr = NULL;

	//Pricing FeeLeg for FwdStartPricer
	if (itsIsSensiLinkedtoCorr) //dans le cas d'un calcul de sensi à la correlation
	{							//on garde les strikes éq mais on utilise la nelle nappe de correlations
		CpyFwdStartCorr = (ICM_Correlation*)GetFwdStartCorr()->Clone();
		((ICM_Smile_Correlation*)CpyFwdStartCorr)->SetVolCurves(((ICM_Smile_Correlation*)initCorr)->GetVolCurves());
		model->SetCorrelationPtr(CpyFwdStartCorr);
	}
	else if (itsForcedRescaling) {initCorr->ResetBetas();}
	else {model->SetCorrelationPtr(GetFwdStartCorr());}

	double Price1 = itsFwdStartPricer->ComputePrice(measure);

	//Pricing DefLeg for FwdEndPricer
	if (itsIsSensiLinkedtoCorr)
	{
		CpyFwdEndCorr = (ICM_Correlation*)GetFwdEndCorr()->Clone();
		((ICM_Smile_Correlation*)CpyFwdEndCorr)->SetVolCurves(((ICM_Smile_Correlation*)initCorr)->GetVolCurves());
		model->SetCorrelationPtr(CpyFwdEndCorr);
	}
	else if (itsForcedRescaling) {initCorr->ResetBetas();}
	else {model->SetCorrelationPtr(GetFwdEndCorr());}

	double Price2 = itsFwdEndPricer->ComputePrice(measure);

	if (!itsForcedRescaling)
	{	model->SetCorrelationPtr(initCorr);
		model->SetCorrelation(GetFwdEndCorr());	}

	if (CpyFwdStartCorr) delete CpyFwdStartCorr;
	if (CpyFwdEndCorr) delete CpyFwdEndCorr;

	double value = Price2-Price1;
	SetPriceFollowMode(measure,value);

	if (measure==qCMPPRICE)
	{
		double value = 0.;
		value = itsFwdEndPricer->Price(qCMPFEELEGPV)-itsFwdStartPricer->Price(qCMPFEELEGPV);
		SetPriceFollowMode(qCMPFEELEGPV,value);
		value = itsFwdEndPricer->Price(qCMPDEFLEGPV)-itsFwdStartPricer->Price(qCMPDEFLEGPV);
		SetPriceFollowMode(qCMPDEFLEGPV,value);
		// value = itsFwdEndPricer->GetAccruedPV()-itsFwdStartPricer->GetAccruedPV();
		// SetPriceFollowMode(qCMPACCRUED,value);
		value = itsFwdEndPricer->Price(qCMPACCRUED)-itsFwdStartPricer->Price(qCMPACCRUED);
		SetPriceFollowMode(qCMPACCRUED,value);
	}

	return (value);
}

