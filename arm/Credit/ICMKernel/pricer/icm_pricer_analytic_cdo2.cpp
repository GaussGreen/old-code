#include "firsttoinc.h"
#include "ICMKernel\pricer\ICM_Pricer_analytic_cdo2.h"
//#include "ICMKernel\pricer\ICM_Pricer_homogeneous.h"
//#include "ICMKernel\inst\icm_nthtd.h"
#include "ICMKernel\inst\icm_cdo2.h"
#include "ICMKernel\crv\ICM_Constant_Piecewise.h"
#include "ICMKernel\glob\icm_betas_correlation.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\modelmulticurves.h"

// --------------------------------------------------------------------
// Method Set 
// --------------------------------------------------------------------
void 
ICM_Pricer_Analytic_Cdo2::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	ICM_Pricer_Distrib::Set(sec, mod,ICM_Parameters(),asof);
	// SetParameters((ICM_Parameters*)parameters->Clone());
	SetParameters(parameters); 
	CptTranchesModel();
}


// --------------------------------------------------------------------
// Compute Price
// --------------------------------------------------------------------
double ICM_Pricer_Analytic_Cdo2::ComputePrice(qCMPMETH measure)
{
	ARM_Security* Initial_Security = GetSecurity();
	ARM_Model*	Initial_Model = GetModel();

	//Set Intermediate Model & Security
	ICM_Pricer_Distrib::Set(GetIntermediateSecurity(),GetIntermediateModel(),GetParameters(),GetAsOfDate());

	double Price = ICM_Pricer::ComputePrice(measure  );

	//Reset Previous model
	ICM_Pricer_Distrib::Set(Initial_Security,Initial_Model,GetParameters(),GetAsOfDate());

	return (Price);
}

// --------------------------------------------------------------------
// Quick Generation of Intermediate Model
// --------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2::CptTranchesModel()
{
	if (itsIntermediateModel0) return;
	
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	const ICM_DefaultCurve* RefCurve = model->GetDefaultCurves(0);

	// ARM_Date Maturity = cdo2->GetMaturity();
	ARM_Date Maturity = cdo2->GetEndDateNA() ;
	int size = 0;
	int size_adj = 0;
	int nbnames = 0;
	char** Tenors=NULL;
	// char** UnionIssuers=NULL;
	// ARM_Vector* notionals = NULL ;
	int k = 0;
	
	const ICM_Parameters& FlowsMatrix = GetParameters();

	ARM_Vector* BETA = FlowsMatrix.GetColVect("BETA");
	if (!BETA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no BETA vector found");

	double Beta_fix = BETA->Elt(0);


	std::vector<std::string> UnionIssuers; 
	ARM_Vector notionals; 
	// cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals,nbnames) ;
	cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals) ;
	//Tenors = model->GetUnionTenorsForCollateral(UnionIssuers,nbnames,sizeTenors);
	ARM_Vector VectorYF;
	model->GetUnionYFForCollateral(UnionIssuers,VectorYF);
	size = VectorYF.size(); 
	

	// delete notionals;
	// FreePointerTabChar(UnionIssuers,nbnames);

	//Suppression of non usefull Buckets-----------------------
	for (k=0;k<size;k++)
	{
		ARM_Date End = model->GetStartDate() + K_YEAR_LEN*VectorYF[k];

		if (End>Maturity) {size_adj=k+1;Maturity=End;break;}
	}

	//---------------------------------------------------------

	size = size_adj;
	
	int sizePtf = cdo2->GetPortfolio()->GetNbSec();

	// ICM_DefaultCurve** DefCurves = new ICM_DefaultCurve*[sizePtf];
	std::vector<const ICM_DefaultCurve*> DefCurves(sizePtf); 
	ARM_Vector NotionalsUnderlying (sizePtf); 
	std::vector<std::string> TranchesLabels(sizePtf); 
	
	//Calcul des Spread Implicites par maturité
	for (int i=0;i<sizePtf;i++)
	{
		ARM_Vector spread(size,0.);
		ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(i);
		ICM_Mez* Mez = (ICM_Mez*) InitMez->Clone();
		// 14514 NotionalsUnderlying[i] = InitMez->GetInitialNotional();
		NotionalsUnderlying.Elt(i) = InitMez->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0);

		Mez->CptCashFlowDatesCredit(Maturity);

		ICM_Pricer_Distrib  pricer_tmp ; 
		pricer_tmp.Set(Mez, model,ICM_Parameters(),GetAsOfDate());
		double Spread = pricer_tmp.ComputeSpread()/10000.;
		ICM_DistribLoss  distrib = pricer_tmp.getDistribLoss(); 
		// delete pricer_tmp;

		spread[size-1]=Spread;
		
		for (int j=0;j<size-1;j++)
		{
			ARM_Date End = model->GetStartDate() + K_YEAR_LEN*VectorYF[i];
			Mez->CptCashFlowDatesCredit(End);
			ICM_Pricer_Distrib pricer_tmp ; 
			pricer_tmp .Set(Mez, model,ICM_Parameters(),GetAsOfDate());
			pricer_tmp.setDistribLoss(distrib); 
			Spread = pricer_tmp.ComputeSpread()/10000.;
			spread[j]=Spread;
			// delete pricer_tmp;
		}

		delete Mez;

		std::stringstream sstr; sstr<<"Tranche_"<<i; 
		TranchesLabels[i]=sstr.str(); 
		DefCurves[i] = CptDefaultCurve(VectorYF,spread,model,TranchesLabels[i]);

	}
	
	
	ARM_Vector Recovery( sizePtf) ;
	for (k=0;k<sizePtf;k++) {Recovery[k]=0.;}

	ICM_Beta_Correlation* Corr = new ICM_Beta_Correlation(model->GetStartDate(),Beta_fix,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);

	itsIntermediateModel0 = new ICM_ModelMultiCurves(// sizePtf, 
											   DefCurves,
											   model->GetZeroCurve(),
											   Recovery,
											   Corr);

	if (Corr) delete Corr;
	// if (Recovery) delete[] Recovery;
	
	if (itsIntermediateModel)
		delete itsIntermediateModel;
	SetIntermediateModel((ICM_ModelMultiCurves*)itsIntermediateModel0->Clone());

	if (itsIntermediateSecurity0)
		delete itsIntermediateSecurity0;
	itsIntermediateSecurity0 = NULL;

	ICM_Mez* mez = (ICM_Mez*) GetSecurity();
	mez = (ICM_Mez*)mez->Clone();
	mez->SetIssuersInfos(TranchesLabels,NotionalsUnderlying);

	SetIntermediateSecurity0(mez);

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;
	SetIntermediateSecurity(dynamic_cast<ARM_Security*>(mez->Clone()));

	for (k=0;k<sizePtf;k++) delete DefCurves[k];

	// SetExpectedLoss(NULL);
	clearDistribLoss() ; 

}


// --------------------------------------------------------------------------------------------------------------------
// Quick Generation of Default Curve
// --------------------------------------------------------------------------------------------------------------------
ICM_DefaultCurve* 
ICM_Pricer_Analytic_Cdo2::CptDefaultCurve(const ARM_Vector& YFvector,const ARM_Vector& spreads,ICM_ModelMultiCurves* model,const std::string& label)
{
	double recovery = 0.;
	int j=0,size=0;

	size = spreads.size();

	qCDS_ADJ adj = model->GetDefaultCurves(0)->GetCdsAdj();  	
	// ARM_Currency* ccy = model->GetDefaultCurves()[0]->GetCurrency();  	
	std::string ccy = model->GetDefaultCurves(0)->GetCurrency();  	
	bool issummitcurve = model->GetDefaultCurves(0)->GetIsSummitCurve();
	std::string calibrationMethod = model->GetDefaultCurves(0)->GetCalibrationData();
	qDEFCURVE_CALIB_ALGO calibrationAlgo= model->GetDefaultCurves(0)->GetCalibrationAlgo();
	int lag = model->GetDefaultCurve(0)->GetLagStartDate();

	ICM_DefaultCurve* DefCurve = new ICM_Constant_Piecewise(model->GetStartDate(),
												 YFvector,
												 spreads,
												 recovery,
												 model->GetZeroCurve(),
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												 adj,
												 ccy,
												 label,
												 issummitcurve,
												 NULL,false,//2 NULL,
												 K_QUARTERLY,
												 calibrationAlgo,
												 calibrationMethod,lag);			 

	return (DefCurve);
	
}


// --------------------------------------------------------------------
// View Method
// --------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Analytic CDO² Pricer ----------------- \n\n");

	if (itsIntermediateModel)
	{
	for (int i=0; i<itsIntermediateModel->GetNbDefCurves(); i++)
	{
		itsIntermediateModel->GetDefaultCurves(i)->View(id, fOut);
	}
	fprintf(fOut, "\n");
	}

	if (itsIntermediateSecurity) itsIntermediateSecurity->View(id, fOut);
	if (itsIntermediateModel) itsIntermediateModel->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- Terminal CDO Pricer ----------------- \n\n");
	
	ICM_Pricer_Security::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


// --------------------------------------------------------------------
// Quick Generation of Intermediate Model for sensitivity
// --------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2::ResetTranchesModel(const std::string& label)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	const ICM_DefaultCurve* RefCurve = model->GetDefaultCurves(0);

	// ARM_Date Maturity = cdo2->GetMaturity();
	ARM_Date Maturity = cdo2->GetEndDateNA() ;
	int size = 0;
	int size_adj = 0;
	int nbnames = 0;
	// char** UnionIssuers=NULL;
	std::vector<std::string> UnionIssuers; 
	ARM_Vector  notionals ;
	int k = 0;
	
	if (itsIntermediateModel0 == NULL) 
	{
		CptTranchesModel();
		return;
	}

	const ICM_Parameters& FlowsMatrix = GetParameters();

	ARM_Vector* BETA = FlowsMatrix.GetColVect("BETA");
	if (!BETA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no BETA vector found");

	double Beta_fix = BETA->Elt(0);


	// cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals,nbnames) ;
	cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals) ;
	//Tenors = model->GetUnionTenorsForCollateral(UnionIssuers,nbnames,sizeTenors);
	ARM_Vector VectorYF;
	model->GetUnionYFForCollateral(UnionIssuers,VectorYF);
	// delete notionals;
	size = VectorYF.size(); 
	// FreePointerTabChar(UnionIssuers,nbnames);

	//Suppression of non usefull Buckets-----------------------
	for (k=0;k<size;k++)
	{
		ARM_Date End = model->GetStartDate() + K_YEAR_LEN*VectorYF[k];

		if (End>Maturity) {size_adj=k+1;Maturity=End;break;}
	}

	
	//---------------------------------------------------------
	size = size_adj;
	
	int sizePtf = cdo2->GetPortfolio()->GetNbSec();

	// ICM_DefaultCurve** DefCurves = new ICM_DefaultCurve*[sizePtf];
	std::vector<const ICM_DefaultCurve*> DefCurves (sizePtf); 
	// double* NotionalsUnderlying = new double[sizePtf];
	ARM_Vector NotionalsUnderlying(sizePtf); 
	// char** TranchesLabels = new char*[sizePtf];
	std::vector<std::string> TranchesLabels(sizePtf); 

	int* Appartenance = new int[sizePtf];
	memset(Appartenance,'\0',sizeof(int)*sizePtf);

	for (k=0;k<sizePtf;k++)
	{
		ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(k);
		if (InitMez->GetCollateral()->ContainIssuer(label)) {Appartenance[k] = 1;}
	}
	
	//Calcul des Spread Implicites par maturité
	for (int i=0;i<sizePtf;i++)
	{
		std::stringstream sstr; sstr<<"Tranche_"<<i; 
		TranchesLabels[i] = sstr.str() ; // new char[50];
		// sprintf(TranchesLabels[i],"Tranche_%i",i);

		ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(i);
		// 14514 NotionalsUnderlying[i] = InitMez->GetInitialNotional();
		NotionalsUnderlying.Elt(i) = InitMez->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0);

		if (Appartenance[i] == 0)
		{
			ICM_DefaultCurve* item = (ICM_DefaultCurve*) itsIntermediateModel0->GetDefaultCurves(i)->Clone();
			item->SetLabel(TranchesLabels[i]);
			DefCurves[i]=item; 
			continue;
		}

		ARM_Vector spread(size);
		ICM_Mez* Mez = (ICM_Mez*) InitMez->Clone();

		Mez->CptCashFlowDatesCredit(Maturity);

		ICM_Pricer_Distrib  pricer_tmp ; 
		pricer_tmp .Set(Mez, model,ICM_Parameters(),GetAsOfDate());
		double Spread = pricer_tmp.ComputeSpread()/10000.;

		ICM_DistribLoss distribLoss = pricer_tmp.getDistribLoss(); 
		// delete pricer_tmp;

		spread[size-1]=Spread;
		
		for (int j=0;j<size-1;j++)
		{
			ARM_Date End = model->GetStartDate() + K_YEAR_LEN*VectorYF[i];
			Mez->CptCashFlowDatesCredit(End);

			ICM_Pricer_Distrib pricer_tmp;
			pricer_tmp.Set(Mez, model,ICM_Parameters(),GetAsOfDate());
			pricer_tmp.setDistribLoss(distribLoss); 
			Spread = pricer_tmp.ComputeSpread()/10000.;
			spread[j]=Spread;

			// delete pricer_tmp;
		}

		delete Mez;
		DefCurves[i] = CptDefaultCurve(VectorYF,spread,model,TranchesLabels[i]);

	}
	
	if (Appartenance) delete[] Appartenance;

	if (itsIntermediateModel)
		delete itsIntermediateModel;

	ARM_Vector  Recovery ( sizePtf );
	for (k=0;k<sizePtf;k++) {Recovery[k]=0.;}

	ICM_Correlation* Corr = new ICM_Beta_Correlation(model->GetStartDate(),Beta_fix,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);


	itsIntermediateModel = new ICM_ModelMultiCurves(// sizePtf,
											   DefCurves,
											   model->GetZeroCurve(),
											   Recovery,
											   Corr);

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;
	itsIntermediateSecurity = NULL;

	ICM_Mez* mez = (ICM_Mez*) GetSecurity();
	mez = (ICM_Mez*)mez->Clone();
	mez->SetIssuersInfos(TranchesLabels,NotionalsUnderlying);

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;

	SetIntermediateSecurity(mez);

	for (k=0;k<sizePtf;k++) delete DefCurves[k];

	clearDistribLoss(); // SetExpectedLoss(NULL);

}

//	-------------------------------------------------------------------------
// virtual 
void ICM_Pricer_Analytic_Cdo2::BeforePrice(const std::string& label,qSENSITIVITY_TYPE type /** = ICMSPREAD_TYPE**/ ) 
{
	// if (strcmp(label,"NONE") == NULL) 
	if (label=="NONE")
		ResetTranchesModel(label);
	else 
		CptTranchesModel();
} 
//	-------------------------------------------------------------------------

void 
ICM_Pricer_Analytic_Cdo2::Init() 
{
	SetName(ICM_PRICER_ANALYTIC_CDO2);

	itsIntermediateModel = NULL;
	itsIntermediateSecurity = NULL;

	itsIntermediateModel0 = NULL;
	itsIntermediateSecurity0 = NULL;
}
//	-------------------------------------------------------------------------

// virtual 
void ICM_Pricer_Analytic_Cdo2::Reset(void)
{}
//	-------------------------------------------------------------------------

// 
ICM_Pricer_Analytic_Cdo2::~ICM_Pricer_Analytic_Cdo2(void)
{ 
	if (itsIntermediateModel)
		delete itsIntermediateModel;
	itsIntermediateModel = NULL;

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;
	itsIntermediateSecurity = NULL;

	if (itsIntermediateModel0)
		delete itsIntermediateModel0;
	itsIntermediateModel0 = NULL;

	if (itsIntermediateSecurity0)
		delete itsIntermediateSecurity0;
	itsIntermediateSecurity0 = NULL;
}

