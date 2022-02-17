#include "firsttoinc.h"
//#include "ARMKernel\glob\armglob.h"
#include "ICMKernel\pricer\ICM_Pricer_analytic_cdo2_smile.h"
//#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile.h"
//#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\inst\icm_nthtd.h"
#include "ICMKernel\inst\icm_cdo2.h"
#include "ARMKernel\crv\volint.h"
#include "ICMKernel/str/icm_util_str.h"
#include "ICMKernel\util\icm_correlation_fitting.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\modelmulticurves.h"

// --------------------------------------------------------------------
// Method Set 
// --------------------------------------------------------------------
void 
ICM_Pricer_Analytic_Cdo2_Smile::setCorrelTpsDefaut(const ICM_CorrMatrix*src)
{
	if (itsCorrelTpsDefaut) delete itsCorrelTpsDefaut; itsCorrelTpsDefaut=0; 
	if (src) itsCorrelTpsDefaut=dynamic_cast<ICM_CorrMatrix*>(unconst(*src).Clone()); 
}
void ICM_Pricer_Analytic_Cdo2_Smile::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	//la matrice parameters du pricer doit contenir les champs suivant :
    //INTEGRATION_METHOD_UNDERLYING 0   'Underlying method 
	//INTEGRATION_STEP_UNDERLYING 20	'Underlyinq step 
    //FORCED_STRIKE_DOWN_CDO2 3%        'Forced Strike Down for Cdo²
    //FORCED_STRIKE_UP_CDO2 6%          'Forced Strike Up for Cdo²
    //FORCED_CORREL_DOWN_CDO2 20%       'Forced Correlation Down for Cdo²
    //FORCED_CORREL_UP_CDO2 25%         'Forced Correlation Up for Cdo²
    //BETA" >0 or -999                  'Beta for intertranches correlation
    //NUM_SIMUL_CORREL  100000          'NbSimul for intertranches correlation
    //UNDER_RESCAL  1 or 0              'Forced Rescaling for underlyings
	//USE_MF_CDO2  1 or 0               'Use a MF model if necessary (2 or more tranche on the same port

	//la matrice parameters du ptf du Cdo² doit contenir les champs suivant :
	//ITRX_EUR_IG_S2 vecteur de dim = nb Cdos - proportions EUR
	//CDX_NA_IG_S3 vecteur de dim = nb Cdos - proportions USD
	//FORCED_STRIKE_DOWN_UNDERLYING vecteur de dim = nb Cdos - strikes down forcés Cdo1 à n en %
	//FORCED_STRIKE_UP_UNDERLYING vecteur de dim = nb Cdos - strikes up forcés Cdo1 à n en %
	//FORCED_CORREL_DOWN_UNDERLYING vecteur de dim = nb Cdos - correl down forcés Cdo1 à n en %
	//FORCED_CORREL_UP_UNDERLYING vecteur de dim = nb Cdos - correl up forcés Cdo1 à n en %

	

	ICM_Pricer_Distrib_Smile::Set(sec, mod,parameters,asof);
	ARM_Vector* GoTranchesModel = NULL;

	//Integrationstep for underlying Cdos
	if (! parameters.empty())
	{
	ARM_Vector* IntegrationMethod = parameters.GetColVect("INTEGRATION_METHOD_UNDERLYING");
	if (IntegrationMethod) itsIntegrationMethodUnderlying = IntegrationMethod->Elt(0);

	ARM_Vector* IntegrationStep = parameters.GetColVect("INTEGRATION_STEP_UNDERLYING");
	if (IntegrationStep) itsIntegrationStepUnderlying = IntegrationStep->Elt(0);

	ARM_Vector* NumSimulations = parameters.GetColVect("NUM_SIMUL_CORREL");
	if (NumSimulations) itsNumSimulationsForCorrelation = NumSimulations->Elt(0);

	ARM_Vector* ForcedRescalingForUnderlyings = parameters.GetColVect("UNDER_RESCAL");
	if (ForcedRescalingForUnderlyings) itsForcedRescalingForUnderlyings = ForcedRescalingForUnderlyings->Elt(0);

	GoTranchesModel = parameters.GetColVect("CPT_TRANCHE_MODEL");

	ARM_Vector* Vect_CorrelationType = parameters.GetColVect("CORRELATION_TYPE");
	if (Vect_CorrelationType) its_Correlation_Fit_Type = (qCORRELATION_FIT_TYPE) ((int) Vect_CorrelationType->Elt(0));
	}

	if (GoTranchesModel)
	{
	//calcul des correlations individuelles dfes Cdos
	try{ComputeImplicitUnderlyingCorrelations(itsImplicitUnderlyingCorrelations);}
	catch (...)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Pricer Analytic Cdo² Smile : unable to Compute Implicit Underlying Correlations");}

	//calcul de la correlation forcée du CDO²
	try{CptForcedCorrelationForCdo2(parameters);}
 
	catch (...)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Pricer Analytic Cdo² Smile : unable to Compute Forced Correlation For Cdo2");}

	//calcul du model intermédiaire avec rescaling pour chaque underlying Cdo
	CptTranchesModel();
	}

}

// ---------------------------------------------------------------------------------------------------------------------
//	JLA 	adoption : can be nul. 
void 
ICM_Pricer_Analytic_Cdo2_Smile::SetForcedCorrelationForCdo2(ICM_Correlation* correl)
{
	if (itsForcedCorrelationForCdo2) delete itsForcedCorrelationForCdo2; itsForcedCorrelationForCdo2=0 ;
	itsForcedCorrelationForCdo2=correl; 
}
// ---------------------------------------------------------------------------------------------------------------------
//	JLA   adoption
void 
ICM_Pricer_Analytic_Cdo2_Smile::SetImplicitUnderlyingCorrelations(const std::vector<ICM_Smile_Correlation*>&ref)
{
	std::vector<ICM_Smile_Correlation*>::iterator it = itsImplicitUnderlyingCorrelations.begin(); 
	while (it!=itsImplicitUnderlyingCorrelations.end()) 
	{
		delete *it; 
		++it; 
	}
	itsImplicitUnderlyingCorrelations=ref; 
}
// ---------------------------------------------------------------------------------------------------------------------
// calcul de la correlation forcée du CDO² need FORCED_STRIKE_DOWN,FORCED_STRIKE_UP,FORCED_CORREL_DOWN,FORCED_CORREL_UP
// on crée une Correlation à partir des strikes du Cdo² et les corrélations pour ces strikes inputées
// ---------------------------------------------------------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::CptForcedCorrelationForCdo2(const ICM_Parameters& parameters)
{
	// if (parameters == NULL) return;
	if (parameters.empty()) return ;


	ARM_Vector* ForcedStrikeDown = parameters.GetColVect("FORCED_STRIKE_DOWN_CDO2");
	ARM_Vector* ForcedStrikeUp = parameters.GetColVect("FORCED_STRIKE_UP_CDO2");
	ARM_Vector* ForcedCorrelDown = parameters.GetColVect("FORCED_CORREL_DOWN_CDO2");
	ARM_Vector* ForcedCorrelUp = parameters.GetColVect("FORCED_CORREL_UP_CDO2");
	ARM_Vector* Beta = parameters.GetColVect("BETA");
	bool BcptCorrelUpdown = true;
	
	//Use sectorial correlation paramters if necessary
	ARM_Vector* UseMF = parameters.GetColVect("USE_MF_CDO2");
	int DecompWithMF = 0;
	
	ICM_ModelMultiCurves* MMC = (ICM_ModelMultiCurves*) GetModel();
	ICM_Security* security = (ICM_Security*) GetSecurity();
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) security;
	int sizeptf = cdo2->GetPortfolio()->GetNbSec();
	ARM_Date AsOf = MMC->GetStartDate();

	//Sectorial Correlation
	vector<int> sectormember;
	int compteur_secteur = 0;

	double dBeta = 0;
	if (Beta)
	{	if (Beta->Elt(0) != CREDIT_DEFAULT_VALUE) dBeta = Beta->Elt(0); }

	//Test Sector Correlation
	if (UseMF)
		DecompWithMF = UseMF->Elt(0);

	if (DecompWithMF == 1)
	{
		int i=0;
		std::vector<std::string>issuers ; 
		ICM_QMatrix<int>* Appartient = cdo2->GenAppMatrix(issuers);
	
		sectormember.resize(sizeptf);
		for (i=0; i<sizeptf; i++)
			sectormember[i] = -1;
	
		for (i=0; i<sizeptf; i++)
		{
			//Test deja init
			if (sectormember[i] == -1)
			{
				sectormember[i] = compteur_secteur;
				for (int j=i; j<sizeptf; j++)
				{		
					if (Appartient->ColAsStdVector(i) == Appartient->ColAsStdVector(j))
						sectormember[j] = compteur_secteur;
				}
				compteur_secteur++;
			}
		}
		//Test : more than 1 tranche / portfolio
		//Flat correl or Sector Correl
		if (compteur_secteur == sizeptf)
			its_Correlation_Fit_Type = qCORREL_1F_SINGLE_VALUE;
		else
			its_Correlation_Fit_Type = qCORREL_2F_SAME_INTER_SAME_INTRA;

		//Purge
		if (Appartient) delete Appartient;
	}
	
	// Cas correlations & strikes forcés !Pas de Beta & pas de correlation sectorielle
	if ((ForcedStrikeDown) && (ForcedStrikeUp) && (ForcedCorrelDown) && (ForcedCorrelUp) && (!dBeta) && (its_Correlation_Fit_Type == qCORREL_1F_SINGLE_VALUE))
	{
		SetForcedCorrelationForCdo2(0); 

		ARM_Vector* YT = new ARM_Vector(1,5.);
		ARM_Vector* Strikes = new ARM_Vector(2,0.);
		Strikes->Elt(0) = 100.*ForcedStrikeDown->Elt(0);
		Strikes->Elt(1) = 100.*ForcedStrikeUp->Elt(0);
		ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());
		mVolatility->Elt(0,0) = 100.*ForcedCorrelDown->Elt(0);
		mVolatility->Elt(0,1) = 100.*ForcedCorrelUp->Elt(0);

		ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,
												security->GetCurrencyUnit());	


		ARM_VolCurve** VVOL = new ARM_VolCurve*[1];
		VVOL[0]=VOL;

		vector<string> label(1); 
		label[0]="COD2_INDEX" ; 

		ARM_Vector* propo = new ARM_Vector(1,1.);
		ARM_Vector* strikelow = new ARM_Vector(1,ForcedStrikeDown->Elt(0));
		ARM_Vector* strikehigh = new ARM_Vector(1,ForcedStrikeUp->Elt(0));

		SetForcedCorrelationForCdo2 (
			new ICM_Smile_Correlation(AsOf,"CDO2_FORCED_CORR",(const ARM_VolCurve**)VVOL,label,propo,strikelow,strikehigh)
			);

		if (propo) delete propo;
		if (strikelow) delete strikelow;
		if (strikehigh) delete strikehigh;
		delete VOL; VOL = NULL;
		delete[] VVOL;
		if (YT) delete YT;
	}
	//Sectorial Forced correlztion
	else if ((ForcedStrikeDown) && (ForcedStrikeUp) && (ForcedCorrelDown) && (ForcedCorrelUp) && (!dBeta) && (its_Correlation_Fit_Type != qCORREL_1F_SINGLE_VALUE))
	{
		vector<string> TranchesLabels(sizeptf) ;
		for (int i=0;i<sizeptf;i++)
		{
		
			if (cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME"))
				TranchesLabels[i] = cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME")->Getvalue(i,0) ;
			else
			{
				std::stringstream sstr ; sstr<<"NAME_"<<i; 
				sstr>>TranchesLabels[i] ;
			}
			
		}
		qTWO_FACTORS_CORRELATION_TYPE correlType; 
		switch(its_Correlation_Fit_Type)
		{
		case qCORREL_1F_SINGLE_VALUE: correlType=TFCT_FULL; break ;
		case qCORREL_1F_BETAS: correlType=TFCT_SAME_INTER_DIFF_INTRA; break ;
		case qCORREL_2F_SAME_INTER_SAME_INTRA: correlType=TFCT_SAME_INTER_SAME_INTRA; break ;
		default:
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't convert from CorrelFitType to CorrelationSector Type"); 
		} 

		SetForcedCorrelationForCdo2( new ICM_Correlation_Sector(AsOf,"ysaitpaskoimettre",
									correlType,
									TranchesLabels,
									compteur_secteur,
									// sizeptf,
									sectormember,
									ForcedCorrelDown->Elt(0),
									ForcedCorrelUp->Elt(0))
									);
		// FreePointerTabChar(TranchesLabels,sizeptf);
	}
	else if (dBeta)// Si paramétre Beta :: Génération d'une matrice intertranche
	{
		bool IsNTD = false;
		int sizeptf = cdo2->GetPortfolio()->GetNbSec();
		std::vector<std::string>_UnionIssuers ; 
		int i=0,k=0;
		ICM_QMatrix<int>* MatAppartenance = cdo2->GenAppMatrix(_UnionIssuers);
		int nbnames = _UnionIssuers.size(); 
		ICM_CorrMatrix* useless  = NULL;
		ARM_Vector ProbaDefautTranche;
		ARM_Vector VarianceTranche;
		ARM_Vector Recovery(nbnames);
		// vector<double> Notionals;Notionals.resize(nbnames);
		ICM_QMatrix<double>* TupTdown = new ICM_QMatrix<double>(2,sizeptf);
		vector<string> TranchesLabels (sizeptf); 
		ARM_Vector PdefIssuers(nbnames);
		ICM_CorrMatrix* CorrelTpsDefaut = NULL;
		vector<int> Nbdef;Nbdef.resize(sizeptf);
		double bump = 0.05;
		ICM_CorrMatrix* CorrelTpsDefaut_tmp=NULL;
		ICM_CorrMatrix* CorrelTpsDefaut_tmp_up=NULL;
		ICM_CorrMatrix* CorrelTpsDefaut_tmp_down=NULL;
		double RecovCoef = 1.;


		for (i=0;i<sizeptf;i++)
		{
			if (cdo2->GetPortfolio()->GetSecurities()[i]->GetName() == ICM_NTD)
			{
			IsNTD = true;
			ICM_Nthtd* ntd = (ICM_Nthtd*) cdo2->GetPortfolio()->GetSecurities()[i];
			Nbdef[i] = ntd->GetLastNumDefault(); 
			}
			else
			{
			ICM_Mez* mez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurities()[i];
			(*TupTdown)(0,i) =	mez->GetPercentLow(mez->GetFeeLeg()->GetStartDate()); 
			(*TupTdown)(1,i) =	(*TupTdown)(0,i) + mez->GetPercentHight(mez->GetFeeLeg()->GetStartDate()); 
			}
			if (cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME"))
			{
				TranchesLabels[i] =cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME")->Getvalue(i,0) ;
			}
			else
			{
				std::stringstream sstr ; sstr<<"NAME_"<<i  ;
				sstr>>TranchesLabels[i] ;
			}
		}

		//Recovery coef test
		if (cdo2->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = cdo2->GetCollateral()->GetRecovCoef();

		//  sort issuers CDO2
		vector<int> vpermut_UnionIssuers(nbnames);
		vector<string> vUnionIssuers(nbnames);
		for (int kk =0; kk< nbnames; kk++) {
			//vpermut_UnionIssuers[kk] = kk;
			vUnionIssuers[kk] = string(_UnionIssuers[kk]);
		}
		vector<string> v_UnionIssuersSorted= SortVectorString(vUnionIssuers, vpermut_UnionIssuers, K_INCREASING);
		
		ICM_QMatrix<int>* MatApp = cdo2->AppMatrix(v_UnionIssuersSorted);

	#ifdef _DEBUG		
		FILE* pFile =NULL;
		pFile = fopen("C:\\temp\\SortTrancheFilleEtBasket.txt", "w");
		fprintf(pFile, "Appartenance non sortée : \n"); 
		MatAppartenance->View("",pFile); 
		fprintf(pFile, "Appartenance Sortée : \n");
		MatApp->View("",pFile); 
		fprintf(pFile, "Names : \n"); 
		for ( i =0; i<nbnames ; i++) {
			fprintf(pFile, "Phirst security %d : %s ------> NewFirstSecurity : %s\n",i, _UnionIssuers[i].c_str(), v_UnionIssuersSorted[i].c_str());
		}
		for ( i =0; i< nbnames; i++) {
			fprintf(pFile, "permutation %d : %s etait en position %d \n", i, v_UnionIssuersSorted[i].c_str(), vpermut_UnionIssuers[i]); 
		}
		
		if(pFile) fclose(pFile);
	#endif

		ARM_Date Maturity = cdo2->GetEndDateNA() ;
		ARM_Vector Notionals(nbnames) ;
		ARM_Vector collatNotionals ; 
		cdo2->GetCollateral()->GetConstantNotionals(collatNotionals); 
		for (i=0;i<nbnames;i++)
		{
			PdefIssuers[i] = MMC->GetDefaultCurve(_UnionIssuers[vpermut_UnionIssuers[i]])->DefaultProba(Maturity );
			Recovery[i] = MIN(MMC->GetRecoveryRate(_UnionIssuers[vpermut_UnionIssuers[i]])*RecovCoef,1);
			if (cdo2->GetBinaryFlg()) Recovery[vpermut_UnionIssuers[i]] = cdo2->GetBinary();
			//ICM_Ftd* ftd = (ICM_Ftd*) cdo2->GetPortfolio()->GetSecurities()[vpermut_UnionIssuers[i]];
			int issuerRank = cdo2->GetCollateral()->getIssuerPosition(_UnionIssuers[vpermut_UnionIssuers[i]]); 
			Notionals[i] =collatNotionals[issuerRank];
		}

		if (IsNTD)
		{
			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,					
								   sizeptf,
								   MatApp->Getnbrows(),
								   Notionals,
								   Recovery,					
								   dBeta, 
								   Nbdef,			
								   *MatApp,	
								   PdefIssuers,
								   TranchesLabels,
								   useless,	
								   CorrelTpsDefaut,
								   ProbaDefautTranche,
								   VarianceTranche,false );


			if (BcptCorrelUpdown) {
			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,					
								   sizeptf,
								   MatApp->Getnbrows(),
								   Notionals,
								   Recovery,					
								   sqrt(MIN(MAX(dBeta*dBeta+bump,0.),1.)), 
								   Nbdef,			
								   *MatApp,	
								   PdefIssuers,
								   TranchesLabels,
								   useless,	
								   CorrelTpsDefaut_tmp_up,
								   ProbaDefautTranche,
								   VarianceTranche,false);


			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,					
								   sizeptf,
								   MatApp->Getnbrows(),
								   Notionals,
								   Recovery,					
								   sqrt(MIN(MAX(dBeta*dBeta-bump,0.),1.)), 
								   Nbdef,			
								   *MatApp,	
								   PdefIssuers,
								   TranchesLabels,
								   useless,	
								   CorrelTpsDefaut_tmp_down,
								   ProbaDefautTranche,
								   VarianceTranche,false);
			}

		}
		else
		{
			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,
								 sizeptf,
								 MatApp->Getnbrows(),
								 Notionals,
								 Recovery,					
								 dBeta, 
								 TupTdown,			
								 *MatApp,	
								 PdefIssuers,
								 TranchesLabels,
								 useless,	
								 CorrelTpsDefaut,
								 ProbaDefautTranche,
								 VarianceTranche,false);

			if (BcptCorrelUpdown){
			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,
								 sizeptf,
								 MatApp->Getnbrows(),
								 Notionals,
								 Recovery,					
								 sqrt(MIN(MAX(dBeta*dBeta+bump,0.),1.)),
								 TupTdown,			
								 *MatApp,	
								 PdefIssuers,
								 TranchesLabels,
								 useless,	
								 CorrelTpsDefaut_tmp_up,
								 ProbaDefautTranche,
								 VarianceTranche,false);


			MatriceCorrelInterTranche(itsNumSimulationsForCorrelation,
								 sizeptf,
								 MatApp->Getnbrows(),
								 Notionals,
								 Recovery,					
								 sqrt(MIN(MAX(dBeta*dBeta-bump,0.),1.)),
								 TupTdown,			
								 *MatApp,	
								 PdefIssuers,
								 TranchesLabels,
								 useless,	
								 CorrelTpsDefaut_tmp_down,
								 ProbaDefautTranche,
								 VarianceTranche,false);
			}

			// We store a copy of CorrelTpsDefaut, for later view. 
			setCorrelTpsDefaut(CorrelTpsDefaut); 

		}

		ICM_Correlation*	TheCorrelation;
		TheCorrelation	=	NULL;
		TheCorrelation	=	MMC->GetCorrelation();

		double	inter_correl=0,inter_correl_up=0,inter_correl_down=0 ;
		double	intra_correl=0,intra_correl_up=0,intra_correl_down=0 ;
		int n_inter=0;
		int n_intra=0;
		vector<double>	Betas(compteur_secteur);
		vector<double>	Lambdas(compteur_secteur);

		switch (its_Correlation_Fit_Type)
		{
			case	qCORREL_1F_SINGLE_VALUE:

				itsAverageCorrTpsDefault_Down = CorrelTpsDefaut->Average();
				itsAverageCorrTpsDefault_Up = itsAverageCorrTpsDefault_Down;
				if (BcptCorrelUpdown) 
				{
					itsAverageCorrTpsDefault_positive_bump_Down = CorrelTpsDefaut_tmp_up->Average();
					itsAverageCorrTpsDefault_positive_bump_Up = itsAverageCorrTpsDefault_positive_bump_Down;
					itsAverageCorrTpsDefault_negative_bump_Down = CorrelTpsDefaut_tmp_down->Average();
					itsAverageCorrTpsDefault_negative_bump_Up = itsAverageCorrTpsDefault_negative_bump_Down;
				}
				
				SetForcedCorrelationForCdo2( new ICM_Beta_Correlation(MMC->GetStartDate(),NEG_SQRT(itsAverageCorrTpsDefault_Down),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0) );

				break;

			case	qCORREL_1F_BETAS:
				
				break;

			case	qCORREL_2F_SAME_INTER_SAME_INTRA:
				
				for (i=0; i<sizeptf;i++)
				{
					for (int j=i+1; j<sizeptf;j++)
					{
						if (sectormember[i] == sectormember[j])
						{
							intra_correl += CorrelTpsDefaut->GetMatrix().Getvalue(i,j);
							if (BcptCorrelUpdown)
							{
								intra_correl_up += CorrelTpsDefaut_tmp_up->GetMatrix().Getvalue(i,j);
								intra_correl_down += CorrelTpsDefaut_tmp_down->GetMatrix().Getvalue(i,j);
							}
							n_intra++;
						}
						else
						{
							inter_correl += CorrelTpsDefaut->GetMatrix().Getvalue(i,j);
							if (BcptCorrelUpdown)
							{
								inter_correl_up += CorrelTpsDefaut_tmp_up->GetMatrix().Getvalue(i,j);
								inter_correl_down += CorrelTpsDefaut_tmp_down->GetMatrix().Getvalue(i,j);
							}
							n_inter++;
						}

					}
				}

				if (n_intra != 0)
				{
					intra_correl /= n_intra;
					if (BcptCorrelUpdown)
					{
						intra_correl_up /= n_intra;
						intra_correl_down /= n_intra;
					}
				}
				if (n_inter != 0) 
				{
					inter_correl /= n_inter;
					if (BcptCorrelUpdown)
					{
						inter_correl_up /= n_inter;
						inter_correl_down /= n_inter;
					}
				}

				// imposed...
				itsAverageCorrTpsDefault_Down = intra_correl;
				itsAverageCorrTpsDefault_Up = inter_correl; 

				if (BcptCorrelUpdown) 
				{
					itsAverageCorrTpsDefault_positive_bump_Up		=	inter_correl_up;
					itsAverageCorrTpsDefault_positive_bump_Down		=	intra_correl_up;
					itsAverageCorrTpsDefault_negative_bump_Up	=	inter_correl_down;
					itsAverageCorrTpsDefault_negative_bump_Down	=	intra_correl_down;
				}
				SetForcedCorrelationForCdo2( new ICM_Correlation_Sector(
									AsOf,
									"ysaitpakoimettre",
									TFCT_SAME_INTER_SAME_INTRA, // alias qCORREL_2F_SAME_INTER_SAME_INTRA
									TranchesLabels,
									compteur_secteur,
									sectormember,
									intra_correl,
									inter_correl)
									);
				break;

			case	qCORREL_2F_SAME_INTER_DIFF_INTRA:

				ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
				break;

			case	qCORREL_2F_DIFF_INTER_DIFF_INTRA:
				
				ICMTHROW(ERR_INVALID_DATA,"Not implemented  yet!");				
				break;
		}
		itsAverageCorrTpsDefaultFlg = true;

		// MEMORY
		if (CorrelTpsDefaut) delete CorrelTpsDefaut;CorrelTpsDefaut=NULL;
		if (CorrelTpsDefaut_tmp_up) delete CorrelTpsDefaut_tmp_up;CorrelTpsDefaut_tmp_up=NULL;
		if (CorrelTpsDefaut_tmp_down) delete CorrelTpsDefaut_tmp_down;CorrelTpsDefaut_tmp_down=NULL;
		if (TupTdown) delete TupTdown;
		if (MatApp) delete MatApp;
		if (MatAppartenance) delete MatAppartenance; 
	}
	else //On utilise la base correlation standart
	{
		
		// JLA itsForcedCorrelationForCdo2 = (ICM_Correlation*) MMC->GetCorrelation()->Clone();
		SetForcedCorrelationForCdo2 ( dynamic_cast<ICM_Correlation*> ( MMC->GetCorrelation()->Clone() ) );
		//On force le rescaling 
		itsForcedCorrelationForCdo2->ResetBetas();
		ICMTHROW(ERR_INVALID_ARGUMENT,"Pricer Analytic Cdo² Smile : \n unable to construct pricer without parameters (FORCED) or Beta ");

	}
}

// ---------------------------------------------------------------------------------------------------------------------
// calcul des correlations des underlying Cdos
// on crée une Correlation à partir des strikes des Cdos et les corrélations pour ces strikes inputées
// ---------------------------------------------------------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::ComputeImplicitUnderlyingCorrelations(vector<ICM_Smile_Correlation*>& UnderlyingCorrelations,const bool& Reset) 
{

	ARM_Model* model = GetModel();
	ARM_Date AsOf = GetModel()->GetStartDate();
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) GetSecurity();
	int sizePtf = cdo2->GetPortfolio()->GetNbSec();
	int i=0;

	ICM_Parameters* parameters = cdo2->GetPortfolio()->GetParams();

	ARM_Vector* ForcedStrikeDown = parameters?parameters->GetColVect("FORCED_STRIKE_DOWN_UNDERLYING"):0;
	ARM_Vector* ForcedStrikeUp = parameters?parameters->GetColVect("FORCED_STRIKE_UP_UNDERLYING"):0;
	ARM_Vector* ForcedCorrelDown = parameters?parameters->GetColVect("FORCED_CORREL_DOWN_UNDERLYING"):0;
	ARM_Vector* ForcedCorrelUp = parameters?parameters->GetColVect("FORCED_CORREL_UP_UNDERLYING"):0;

	if (UnderlyingCorrelations.size()>0)
	{
		for (int i=0;i<UnderlyingCorrelations.size();i++)
		{
			if (UnderlyingCorrelations[i]) delete UnderlyingCorrelations[i];
			UnderlyingCorrelations[i] = NULL;
		}
	}
	UnderlyingCorrelations.clear();
	UnderlyingCorrelations.resize(sizePtf);

	if ((ForcedStrikeDown) && (ForcedStrikeUp) && (ForcedCorrelDown) && (ForcedCorrelUp) && (!itsForcedRescalingForUnderlyings))
	{

		ARM_Date AsOf = GetModel()->GetStartDate();
		int size = ForcedStrikeDown->GetSize();

		vector<ARM_Currency*> Ccy;Ccy.resize(size);
		vector<double> strikeDown;strikeDown.resize(size);
		vector<double> strikeUp;strikeUp.resize(size);
		vector<double> CorrelDown;CorrelDown.resize(size);
		vector<double> CorrelUp;CorrelUp.resize(size);
		vector<string> IndexName;IndexName.resize(size);
		
		for (i=0; i<size;i++)
		{
		Ccy[i] = new ARM_Currency("EUR");
		strikeDown[i] = ForcedStrikeDown->Elt(i);
		strikeUp[i] = ForcedStrikeUp->Elt(i);
		CorrelDown[i] = ForcedCorrelDown->Elt(i);
		CorrelUp[i] = ForcedCorrelUp->Elt(i);
		char TMP[255]; sprintf(TMP,"FORCED_CORR_UNDER%d",i);
		IndexName[i] = (string)TMP;
		}

		FixedBaseCorrelationMult_Vector(AsOf,Ccy,strikeDown,strikeUp,CorrelDown,CorrelUp,IndexName,UnderlyingCorrelations);

		for (i=0; i<size;i++) if (Ccy[i]) delete	Ccy[i];
	}
	else //Real Base correlations for underlying Cdos
	{
		ICM_Correlation* Correlation = ((ICM_ModelMultiCurves*)model)->GetCorrelation();
		ICM_Smile_Correlation* Smile_correl = (ICM_Smile_Correlation*)Correlation;

		int sizeprop = 0;
		ARM_Vector** Proportions_IndexName = NULL;
		//Récupération des proportions
		vector<string> IndexName; 
		CptProportionsIndexName(sizeprop,Proportions_IndexName,IndexName);	

		for (i=0; i<sizePtf;i++)
		{
			UnderlyingCorrelations[i] = (ICM_Smile_Correlation*) Smile_correl->Clone();
			ARM_Vector V(sizeprop,0.);
			for (int il=0; il<sizeprop;il++) {V.Elt(il)=Proportions_IndexName[il]->Elt(i);

			UnderlyingCorrelations[i]->SetProportions(&V);}

		}

		if (Proportions_IndexName) delete[] Proportions_IndexName;
	}
}

// --------------------------------------------------------------------
// Compute Price
// --------------------------------------------------------------------
double ICM_Pricer_Analytic_Cdo2_Smile::ComputePrice(qCMPMETH measure )
{
	ARM_Security* Initial_Security = GetSecurity();
	ARM_Model*	Initial_Model = GetModel();
	double Price=0.;
	//Set Intermediate Model & Security
	try {
		ICM_Pricer_Distrib_Smile::Set(GetIntermediateSecurity(),GetIntermediateModel(),GetParameters(),GetAsOfDate());
		Price = ICM_Pricer::ComputePrice(measure);
	}catch(...){
		ICM_Pricer_Distrib_Smile::Set(Initial_Security,Initial_Model,GetParameters(),GetAsOfDate());
		ICMTHROW(ERR_INVALID_ARGUMENT,"Pricer Analytic Cdo² Smile : unable to Compute Price");
	}
	//Reset Previous model
	ICM_Pricer_Distrib_Smile::Set(Initial_Security,Initial_Model,GetParameters(),GetAsOfDate());

	return (Price);
}



double ICM_Pricer_Analytic_Cdo2_Smile::DoPrice(qCMPMETH measure){
	double res =0.;
	ARM_Date ExecutionDate = GetModel()->GetStartDate();
	switch(measure) {
	case qCMPAVGCORRDEF:
		if (! GetAverageCorrTpsDefaultFlg() ) ComputeGetAverageCorrTpsDefault();
		res = GetAverageCorrTpsDefault_Down(); 
		break;
	default :
		string pmeas;
		ICM_EnumsCnv::toString( measure, pmeas);
		ICMTHROW( ERR_INVALID_ARGUMENT,
			"Pricer Analytic Cdo² Smile : unable to Compute Pricing Measure :"<< pmeas);
	}
	return res;	

}

void ICM_Pricer_Analytic_Cdo2_Smile::ComputeGetAverageCorrTpsDefault(){
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) GetSecurity();
	ICM_Parameters* parameters = cdo2->GetPortfolio()->GetParams();
	if (parameters) CptForcedCorrelationForCdo2( *parameters)  ; 
	else CptForcedCorrelationForCdo2( ICM_Parameters() )  ; 
}



// -------------------------------------------------------------------------------------------
// this function returns proportions of each index (with index names) for each underlying index
// -------------------------------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::CptProportionsIndexName(int& sizeprop,
															 ARM_Vector**& Proportions_IndexName,
															 vector<string>&IndexName)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;

	//Gestion du rescaling pour les Underlying cdo
	if (cdo2->GetPortfolio()->GetParams())
	if (model->GetCorrelation()->GetName() == ICM_SMILE_CORRMATRIX)
	{
		ICM_Smile_Correlation* corr = (ICM_Smile_Correlation*) model->GetCorrelation();
		IndexName = corr->GetLabels();
		sizeprop  = corr->GetIndexSize();
		ICM_Parameters* parameters = cdo2->GetPortfolio()->GetParams();
		int size_parameters = parameters->GetDblParams()->GetNbCol();

		Proportions_IndexName = new ARM_Vector*[sizeprop];
		int k = 0;
		
		for (int i=0;i<size_parameters;i++)
		{
			int flg = false;
			const char* name = parameters->GetDblParams()->GetColName(i);

			std::string n = name;
			int pos = n.find("FORCED",0);
			if (pos != string::npos) continue;	

			for (int j=0;j<sizeprop;j++)
			{
			if (strcmp(name,IndexName[j].c_str()) == NULL)
			{flg = true; break;}
			}

			if (flg == false)
				continue;

			Proportions_IndexName[k] = parameters->GetDblParams()->GetColVect((char*)name);
			k++;
		}
	}
}



// ---------------------------------------------------------------------------------------------
// Quick Generation of Intermediate Model
// ---------------------------------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::CptTranchesModel()
{
	if (itsIntermediateModel0) return;

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	const ICM_DefaultCurve* RefCurve = model->GetDefaultCurves(0);
	ARM_Date AsOf = RefCurve->GetAsOfDate();

	ARM_Date Maturity = cdo2->GetEndDateNA();
	int sizeTenors = 0;
	int sizeTenors_Adj = 0;
	int nbnames = 0;
	ARM_Vector Tenors_; 
	std::vector<std::string> UnionIssuers; 
	ARM_Vector notionals ;
	int k = 0;
	ARM_Vector Proportions_new;
	
	// const ICM_Parameters& FlowsMatrix = GetParameters();

	ICM_Correlation* Correlation = model->GetCorrelation();
	ICM_Smile_Correlation* Smile_correl = (ICM_Smile_Correlation*)Correlation->Clone();

	ARM_Vector* Proportions = NULL;
	Smile_correl->GetProportions(Proportions);

	cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals) ;
	
	// ARM_Vector _Tenors; 
	model->GetUnionYFForCollateral(UnionIssuers,Tenors_);
	// { Tenors.resize(_Tenors.size()); for(int i=0;i<_Tenors.size();i++) Tenors[i]=_Tenors[i]; }
	sizeTenors = Tenors_.size();


	double Matu = (Maturity-AsOf)/365.;
	for (k=0;k<sizeTenors;k++) 
	{if ((Tenors_[k]>Matu) || CHECK_EQUAL(Tenors_[k],Matu)) {sizeTenors_Adj=k+1;break;}}

	sizeTenors = sizeTenors_Adj;
	ARM_Vector Tenors(sizeTenors_Adj); 
	for (k=0;k<sizeTenors-1;k++) Tenors[k]=Tenors_[k]; 
	// Tenors.Resize(sizeTenors); 
	Tenors[sizeTenors-1] = Matu;
	
	int sizePtf = cdo2->GetPortfolio()->GetNbSec();

	// ICM_DefaultCurve** DefCurves = new ICM_DefaultCurve*[sizePtf];
	std::vector<const ICM_DefaultCurve*> DefCurves (sizePtf); 
	ARM_Vector NotionalsUnderlying (sizePtf); 
	vector<double> TradedCoefs;TradedCoefs.resize(sizePtf);
	std::vector<std::string> TranchesLabels(sizePtf); 
	int i=0,j=0;
	k = 0;

	// ICM_Pricer_Security* pricer_tmp = NULL;
	// double* Recovery = new double[sizePtf];
	ARM_Vector Recovery (sizePtf); 

	//Calcul des Spread Implicites par maturité
	for (i=0;i<sizePtf;i++)
	{
		// TranchesLabels[i] = new char[50];
		if (cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME"))
		{
			TranchesLabels[i]=cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME")->Getvalue(i,0) ;
		}
		else
		{
			std::stringstream sstr ; sstr<<"NAME_"<<i; 
			TranchesLabels[i]=sstr.str(); 
		}

		double dRecovery = 0.;
		model->SetCorrelation(itsImplicitUnderlyingCorrelations[i]);

		ARM_Vector spread (sizeTenors); 
		ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(i);
		ICM_Mez* Mez = (ICM_Mez*) InitMez->Clone();
		NotionalsUnderlying.Elt(i) = InitMez->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0);
		TradedCoefs[i] = InitMez->GetTradedCoef(); 

		Mez->CptCashFlowDatesCredit(cdo2->GetEndDateNA());

		ICM_Pricer_Distrib_Smile pricer_tmp ; 
		pricer_tmp.Set(Mez, model,ICM_Parameters(),GetAsOfDate());
		//On set les parametres pour le pricing des underlyings Cdo² 
		int parameter = 0;
		GetCopulaType(parameter);
		pricer_tmp.SetCopulaType(parameter);
		pricer_tmp.SetIntegrationMethod(itsIntegrationMethodUnderlying);
		GetFreedomDegree(parameter);
		pricer_tmp.SetFreedomDegree(parameter);
		pricer_tmp.SetIntegrationStep1(itsIntegrationStepUnderlying);
		pricer_tmp.SetRescalType(GetRescalType());
		pricer_tmp.SetTermStructurePricing(GetTermStructurePricing());
		//}

		double Spread = pricer_tmp.ComputeSpread()/10000.;

		if (Mez->GetName()==ICM_NTD) 
		{
			dRecovery=pricer_tmp.AvgRecovery(Mez->GetCollateral()->GetIssuersLabels());
		 //round à 5% sur la recovery
		 dRecovery=((double)((int)(dRecovery*20.+1.e-2)))*0.05;
		}

		//l'expected loss est calculée à maturités
		ICM_DistribLoss distribLoss = pricer_tmp.getDistribLoss() ;  
		// delete pricer_tmp;

		spread[sizeTenors-1]=Spread;

		//For each tenor of the default curve we compute the breakeven spread
		for (int j=0;j<sizeTenors-1;j++)
		{
			ARM_Date End = Tenors[j]*365.+ AsOf.GetJulian();
			Mez->CptCashFlowDatesCredit(End);

			ICM_Pricer_Distrib_Smile pricer_tmp ; 
			pricer_tmp.Set(Mez, model,ICM_Parameters(),GetAsOfDate());
			//On set les parametres pour le pricing des underlyings Cdo² 
			int parameter = 0;
			GetCopulaType(parameter);
			pricer_tmp.SetCopulaType(parameter);
			pricer_tmp.SetIntegrationMethod(itsIntegrationMethodUnderlying);
			GetFreedomDegree(parameter);
			pricer_tmp.SetFreedomDegree(parameter);
			pricer_tmp.SetIntegrationStep1(itsIntegrationStepUnderlying);
			pricer_tmp.SetRescalType(GetRescalType());
			pricer_tmp.SetTermStructurePricing(GetTermStructurePricing());
			//}

			//On set l'expected loss calculée à maturité
			pricer_tmp.setDistribLoss(distribLoss); 
			Spread = pricer_tmp.ComputeSpread()/10000.;
			spread[j]=Spread;
			// delete pricer_tmp;
		}

		//Now we can build the default curve for the underlying Cdo N°i
		//On impose un recovery = 0 pour chaque courbe sous-jacente
		string ReferenceCurveName = (string) InitMez->GetCollateral()->GetIssuersLabels(0);
		DefCurves[i] = CptDefaultCurve(Tenors,spread,model,TranchesLabels[i],dRecovery,ReferenceCurveName);
		Recovery[i] = dRecovery;

		((ICM_Smile_Correlation*)itsImplicitUnderlyingCorrelations[i])->SetSlices(((ICM_Smile_Correlation*)model->GetCorrelation())->GetSlices());
		((ICM_Smile_Correlation*)itsImplicitUnderlyingCorrelations[i])->SetAlready_rescal(true);	


		delete Mez;
	}

	//Now we can set the proportion of index for Cdo²
	Proportions_new = ARM_Vector(Proportions->GetSize(),0.);
	Proportions_new.InitElt(Proportions->GetSize()-1,1.);
	Smile_correl->SetProportions(&Proportions_new);

	//On remet la correlation initiale
	model->SetCorrelation(Smile_correl);

	itsIntermediateModel0 = new ICM_ModelMultiCurves(/** sizePtf,
											   DefCurves, **/ 
											   DefCurves,
											   model->GetZeroCurve(),
											   Recovery,
											   Smile_correl,
											   NULL,
											   true,
											   model->GetCpnInfCurve(),
											   model->GetCpnIRCurve());

	if (itsForcedCorrelationForCdo2)
		itsIntermediateModel0->SetCorrelation(itsForcedCorrelationForCdo2);

	// if (Recovery) delete[] Recovery;

	if (itsIntermediateModel)
		delete itsIntermediateModel;
	SetIntermediateModel((ICM_ModelMultiCurves*)itsIntermediateModel0->Clone());

	//On réajuste des notionels avec le traded notional
	for (i=0;i<sizePtf;i++)	NotionalsUnderlying.Elt(i) *= TradedCoefs[i];

	
	ICM_Mez* mez = (ICM_Mez*) GetSecurity();
	

	mez = (ICM_Mez*)mez->Clone();
	mez->SetIssuersInfos(TranchesLabels,NotionalsUnderlying);

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;
	SetIntermediateSecurity(mez);

	for (k=0;k<sizePtf;k++) delete DefCurves[k];
	// delete[] DefCurves;
	TradedCoefs.clear();

	Smile_correl->SetProportions(Proportions);
	clearDistribLoss(); 

	if (Proportions) delete Proportions;
	Proportions = NULL;

	if (Smile_correl) delete Smile_correl;
	Smile_correl = NULL;
}

// --------------------------------------------------------------------
// View Method
// --------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t ------------------------- Analytic CDO² Pricer --------------------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

    fprintf(fOut, "Integration step for underlying : %i\n",itsIntegrationStepUnderlying);
	
	switch (itsIntegrationMethodUnderlying)
	{
	case 0:
		fprintf(fOut, "Integration method for underlying : GAUSS_LEGENDRE\n");
		break;
	case 1:
		fprintf(fOut, "Integration method for underlying : GAUSS_HERMITE\n");
		break;
	case 2:
		fprintf(fOut, "Integration method for underlying : TRAPEZE\n");
		break;
	default:
		fprintf(fOut, "Integration method for underlying : UNKNOWN !\n");
		break;
	}
	fprintf(fOut, "Num Simulation step for underlying : %i\n",itsNumSimulationsForCorrelation);
	fprintf(fOut, "Forced rescaling for underlyings : %i\n",itsForcedRescalingForUnderlyings);

    fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t -------------- Implicit Underlyings Correlations ------------------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

	if ((itsImplicitUnderlyingCorrelations.size()>0) && (itsIntermediateModel))
	{
		for (int i=0; i<itsIntermediateModel->GetNbDefCurves(); i++)
		{
		if (itsIntermediateModel) fprintf(fOut, "Tranche : %s\n ",itsIntermediateModel->GetDefaultCurves(i)->GetLabel().c_str());
		if (itsIntermediateModel) fprintf(fOut, "--------------------------------");

		if (itsImplicitUnderlyingCorrelations.size()>i)
		{	if (itsImplicitUnderlyingCorrelations[i]) itsImplicitUnderlyingCorrelations[i]->View(id, fOut); }
		fprintf(fOut, "\n");
		}
	}

	if (itsCorrelTpsDefaut) 
	{
		fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
		fprintf(fOut, "\t\t\t -------------- Correlation Matrix		      ------------------------ \n");
		fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");
		itsCorrelTpsDefaut->View(id,fOut); 
	}

    fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t -------------- Forced Correlations for Cdo² ------------------------ \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

	if (itsForcedCorrelationForCdo2) itsForcedCorrelationForCdo2->View(id, fOut);
	fprintf(fOut,"\nAverage Correlation :%f\n",itsAverageCorrTpsDefault_Down);
	fprintf(fOut,"Average Correlation Bump + :%f\n",itsAverageCorrTpsDefault_positive_bump_Down);
	fprintf(fOut,"Average Correlation Bump - :%f\n\n",itsAverageCorrTpsDefault_negative_bump_Down);

	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
	fprintf(fOut, "\t\t\t ------------------------- CDO² on deduced CDOs --------------------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

	if (itsIntermediateSecurity) itsIntermediateSecurity->View(id, fOut);

	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
	fprintf(fOut, "\t\t\t ------------------------- Intermediate Model for CDO² -------------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

	if (itsIntermediateModel) itsIntermediateModel->View(id, fOut);

	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
	fprintf(fOut, "\t\t\t ----------- Curves details of Intermediate Model for CDO² ---------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");

	if (itsIntermediateModel)
	{
	for (int i=0; i<itsIntermediateModel->GetNbDefCurves(); i++)
	{
		itsIntermediateModel->GetDefaultCurves(i)->View(id, fOut);
	}
	fprintf(fOut, "\n");
	}

	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t ------------------------- average spread on CDO at Maturity -------- \n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");
	ICM_Cdo2* cdo2 = dynamic_cast<ICM_Cdo2*>(GetSecurity());
	if (cdo2){
		for ( int j=0; j< cdo2->GetPortfolio()->GetNbSec(); j++){
			ICM_Ftd* mez = dynamic_cast<ICM_Ftd*>(cdo2->GetPortfolio()->GetSecurities()[j]);
			int NbIssuers = mez->GetCollateral()->GetNbIssuers();	
			double averageSpread = 0.0;
			for (int k = 0; k< NbIssuers; k++) {
				string name = string(mez->GetCollateral()->GetIssuersLabels(k));	
				ICM_ModelMultiCurves* MMC = (ICM_ModelMultiCurves*) GetModel();
				double spread = MMC->GetDefaultCurve(name)->ImpliedSpreadInterpol(mez->GetEndDateNA());
				averageSpread += spread;	
			}
			averageSpread /= NbIssuers;
			fprintf(fOut, "\t CDO %d, average spread at maturity : %lf\n",j, averageSpread);
		}

	}
	fprintf(fOut, "\n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");
	fprintf(fOut, "\t\t\t ---------------- Homogeneous Strike Security Pricer  --------------- \n\n");
	fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");
	
	ICM_Pricer_Distrib_Smile::View(id, fOut);


	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


// --------------------------------------------------------------------
// Quick Generation of Intermediate Model for sensitivity
// --------------------------------------------------------------------
void ICM_Pricer_Analytic_Cdo2_Smile::ResetTranchesModel(const std::string&  label,qSENSITIVITY_TYPE type)
{

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	const ICM_DefaultCurve* RefCurve = model->GetDefaultCurves(0);
	ARM_Date AsOf = RefCurve->GetAsOfDate();

	ARM_Date Maturity = cdo2->GetEndDateNA();
	int sizeTenors = 0;
	int sizeTenors_Adj = 0;
	int nbnames = 0;
	//JLA ARM_Vector Tenors;
	ARM_Vector Tenors_ ;
	// char** UnionIssuers=NULL;
	std::vector<std::string> UnionIssuers; 
	ARM_Vector  notionals ;
	int k = 0;
	
	if (itsIntermediateModel0 == NULL) 
	{CptTranchesModel();return;}

	ICM_Correlation* Correlation = model->GetCorrelation();

	if (type==ICM_CORREL_BET_CDO_UP_TYPE)
	{	if (itsIntermediateModel) delete itsIntermediateModel;
		itsIntermediateModel = (ICM_ModelMultiCurves*) itsIntermediateModel0->Clone();

		// waiting for the fitting interface
		if (its_Correlation_Fit_Type != qCORREL_1F_SINGLE_VALUE)
		{
			ICM_Correlation_Sector	Corr_plus;

			Corr_plus	=	*((ICM_Correlation_Sector*)(GetForcedCorrelationForCdo2()->Clone()));

			// update parameters

			switch (its_Correlation_Fit_Type)
			{
				case	qCORREL_1F_SINGLE_VALUE:
					// ???

				break;

				case	qCORREL_1F_BETAS:
						// ???

					break;

				case	qCORREL_2F_SAME_INTER_SAME_INTRA:

					Corr_plus.SetSingle_Intra_Sector_Correlation(itsAverageCorrTpsDefault_positive_bump_Down);
					Corr_plus.SetSingle_Inter_Sector_Correlation(itsAverageCorrTpsDefault_positive_bump_Up);
					
					break;

				default:
					break;
			}
			
			itsIntermediateModel->SetCorrelation(&Corr_plus); 
		}
		else
		{
			ICM_Beta_Correlation Corr_plus(model->GetStartDate(),NEG_SQRT(itsAverageCorrTpsDefault_positive_bump_Down),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);

			itsIntermediateModel->SetCorrelation(&Corr_plus); 
		}

		clearDistribLoss(); 
		return;
	}

	if (type==ICM_CORREL_BET_CDO_DW_TYPE)
	{	if (itsIntermediateModel) delete itsIntermediateModel;
		itsIntermediateModel = (ICM_ModelMultiCurves*) itsIntermediateModel0->Clone();

		// waiting for the fitting interface
		if (its_Correlation_Fit_Type != qCORREL_1F_SINGLE_VALUE)
		{
			ICM_Correlation_Sector	Corr_moins;
			Corr_moins	=	*((ICM_Correlation_Sector*)(GetForcedCorrelationForCdo2()->Clone()));

			// update parameters

			switch (its_Correlation_Fit_Type)
			{
				case	qCORREL_1F_SINGLE_VALUE:
					// ???

				break;

				case	qCORREL_1F_BETAS:
						// ???

					break;

				case	qCORREL_2F_SAME_INTER_SAME_INTRA:

					Corr_moins.SetSingle_Intra_Sector_Correlation(itsAverageCorrTpsDefault_negative_bump_Down);
					Corr_moins.SetSingle_Inter_Sector_Correlation(itsAverageCorrTpsDefault_negative_bump_Up);
					
					break;

				default:
					break;
			}
			
			itsIntermediateModel->SetCorrelation(&Corr_moins); 
		}
		else
		{
			ICM_Beta_Correlation Corr_moins(GetAsOfDate(),NEG_SQRT(itsAverageCorrTpsDefault_negative_bump_Down),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);

			itsIntermediateModel->SetCorrelation(&Corr_moins); 
		}
		
		clearDistribLoss(); 
		return;
	}

	ICM_Smile_Correlation* Smile_correl = (ICM_Smile_Correlation*)Correlation->Clone();

	cdo2->GetPortfolio()->GetIssuersDatas(UnionIssuers,notionals) ;
	
	// ARM_Vector _Tenors; 
	model->GetUnionYFForCollateral(UnionIssuers,Tenors_);
	// { Tenors.resize(_Tenors.size()); for(int i=0;i<_Tenors.size();i++) Tenors[i]=_Tenors[i]; }
	sizeTenors = Tenors_.size();


	double Matu = (Maturity-AsOf)/365.;
	for (k=0;k<sizeTenors;k++) 
	{if ((Tenors_[k]>Matu) || CHECK_EQUAL(Tenors_[k],Matu)) {sizeTenors_Adj=k+1;break;}}
	sizeTenors = sizeTenors_Adj;
	ARM_Vector  Tenors(sizeTenors) ;
	for(k=0;k<sizeTenors-1;k++) Tenors[k]=Tenors_[k]; 
	Tenors[sizeTenors-1]=Matu;
	
	int sizePtf = cdo2->GetPortfolio()->GetNbSec();

	// ICM_DefaultCurve** DefCurves = new ICM_DefaultCurve*[sizePtf];
	std::vector<const ICM_DefaultCurve*> DefCurves(sizePtf); 
	vector<string> TranchesLabels(sizePtf);
	//char** TranchesLabels = new char*[sizePtf];

	int* Appartenance = new int[sizePtf];
	memset(Appartenance,'\0',sizeof(int)*sizePtf);
	for (k=0;k<sizePtf;k++)
	{
		if ((type!=ICMSPREAD_TYPE)&&(type!=ICMRECOVERY_TYPE)) {Appartenance[k] = 1;}
			else if (((type=ICMSPREAD_TYPE)||(type!=ICMRECOVERY_TYPE)) && label=="NONE" /** strcmp(label,"NONE") ==0 **/ ) {Appartenance[k] = 1;}
				else 			
					{ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(k);
					if (InitMez->GetCollateral()->ContainIssuer(label)) {Appartenance[k] = 1;}}
	}

	ICM_Pricer_Security* pricer_tmp = NULL;
	ARM_Vector Recovery (sizePtf);
	vector<ICM_Smile_Correlation*> ImplicitUnderlyingCorrelations;

	if ((type==ICMCORRELATION_TYPE)||(type==ICMCORREL_STRIKE_DOWN_TYPE)||(type==ICMCORREL_STRIKE_UP_TYPE))
	{ ComputeImplicitUnderlyingCorrelations(ImplicitUnderlyingCorrelations,false);}	

	//Calcul des Spread Implicites par maturité
	for (int i=0;i<sizePtf;i++)
	{
		double dRecovery = 0.;
		char tmp[50];
		//TranchesLabels[i] = new char[50];

		if (cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME")){
			sprintf(tmp ,(const char*)cdo2->GetPortfolio()->GetParams()->GetColVectStr("NAME")->Getvalue(i,0).c_str(),i);
			TranchesLabels[i] = tmp;
		}else{
			sprintf(tmp,"NAME_%i",i);
			TranchesLabels[i] = tmp;
		}

		string Name = TranchesLabels[i];

		ICM_Mez* InitMez = (ICM_Mez*) cdo2->GetPortfolio()->GetSecurity(i);

		if (Appartenance[i] == 0)
		{
			ICM_DefaultCurve* item = (ICM_DefaultCurve*) itsIntermediateModel0->GetDefaultCurves(i)->Clone(); 
			item->SetLabel(TranchesLabels[i]);
			DefCurves[i]= item; 
			Recovery[i] = itsIntermediateModel0->GetRecoveryRates()[i];
			continue;
		}

		// JLA ARM_Vector spread(sizeTenors);
		ARM_Vector spread(sizeTenors); 
		ICM_Mez* Mez = (ICM_Mez*) InitMez->Clone();

		ICM_Smile_Correlation* Corr = NULL;

		if ((type==ICMCORRELATION_TYPE)||(type==ICMCORREL_STRIKE_DOWN_TYPE)||(type==ICMCORREL_STRIKE_UP_TYPE))
		{Corr = (ICM_Smile_Correlation*)itsImplicitUnderlyingCorrelations[i]->Clone();
		 Corr->SetCorrelations(*(ICM_Smile_Correlation*)ImplicitUnderlyingCorrelations[i]);}
		else
		{Corr = (ICM_Smile_Correlation*)itsImplicitUnderlyingCorrelations[i]->Clone();}

		if (type==ICM_INDX_SPREAD_RESCALING){Corr->ResetBetas();}

		model->SetCorrelation(Corr);
		Mez->CptCashFlowDatesCredit(cdo2->GetEndDateNA());

		ICM_Pricer_Distrib_Smile pricer_tmp ; 
		pricer_tmp .Set(Mez, model,ICM_Parameters(),GetAsOfDate());
		//On set les parametres pour le pricing des underlyings Cdo² 
		int parameter = 0;
		GetCopulaType(parameter);
		pricer_tmp.SetCopulaType(parameter);
		pricer_tmp.SetIntegrationMethod(itsIntegrationMethodUnderlying);
		GetFreedomDegree(parameter);
		pricer_tmp.SetFreedomDegree(parameter);
		pricer_tmp.SetIntegrationStep1(itsIntegrationStepUnderlying);
		pricer_tmp.SetRescalType(GetRescalType());
		pricer_tmp.SetTermStructurePricing(GetTermStructurePricing());
		//}

		double Spread = pricer_tmp.ComputeSpread()/10000.;
		// If Spread is higher than 20 000 bps then we consider the name is in default 
		if (Spread > 2.) Spread = 10.;

		if (Mez->GetName()==ICM_NTD) 
		{	dRecovery=pricer_tmp.AvgRecovery(Mez->GetCollateral()->GetIssuersLabels());
			//round à 5% sur la recovery
			dRecovery=((double)((int)(dRecovery*20.+1.e-2)))*0.05;
		}
 
		ICM_DistribLoss distribLoss ( pricer_tmp.getDistribLoss() ) ;
			// map<ORDER_DOUBLE_WITH_ORDER,double>* EL = new map<ORDER_DOUBLE_WITH_ORDER,double>(*pricer_tmp->GetExpectedLoss());
		// delete pricer_tmp;

		spread[sizeTenors-1]=Spread;
		
		for (int j=0;j<sizeTenors-1;j++)
		{
			ARM_Date End = Tenors[j]*365.+ AsOf.GetJulian();
			Mez->CptCashFlowDatesCredit(End);

			ICM_Pricer_Distrib_Smile pricer_tmp ;
			pricer_tmp.Set(Mez, model,ICM_Parameters(),GetAsOfDate());
			//On set les parametres pour le pricing des underlyings Cdo² 
			int parameter = 0;
			GetCopulaType(parameter);
			pricer_tmp.SetCopulaType(parameter);
			pricer_tmp.SetIntegrationMethod(itsIntegrationMethodUnderlying);
			GetFreedomDegree(parameter);
			pricer_tmp.SetFreedomDegree(parameter);
			pricer_tmp.SetIntegrationStep1(itsIntegrationStepUnderlying);
			pricer_tmp.SetRescalType(GetRescalType());
			pricer_tmp.SetTermStructurePricing(GetTermStructurePricing());
			//}	
 
			pricer_tmp.setDistribLoss(distribLoss); 
			Spread = pricer_tmp.ComputeSpread()/10000.;
	
			// If Spread is higher than 20 000 bps then we consider the name is in default 
			if (Spread > 2.) Spread = 10.;
			spread[j]=Spread;

			// delete pricer_tmp;
		}

		//On impose un recovery = 0 pour chaque courbe sous-jacente
		string ReferenceCurveName = InitMez->GetCollateral()->GetIssuersLabels(0);
		DefCurves[i] = CptDefaultCurve(Tenors,spread,model,Name,dRecovery,ReferenceCurveName);
		Recovery[i] = dRecovery;

		if (Mez) delete Mez;
		if (Corr) delete Corr;

	}
	
	if (ImplicitUnderlyingCorrelations.size()>0)
	{
		for (int i=0;i<ImplicitUnderlyingCorrelations.size();i++)
		{
			if (ImplicitUnderlyingCorrelations[i]) delete ImplicitUnderlyingCorrelations[i];
			ImplicitUnderlyingCorrelations[i] = NULL;
		}
	}
	ImplicitUnderlyingCorrelations.clear();

	if (Appartenance) delete[] Appartenance;

	if (itsIntermediateModel) delete itsIntermediateModel;

	itsIntermediateModel = new ICM_ModelMultiCurves(// sizePtf,
											   DefCurves,
											   model->GetZeroCurve(),
											   Recovery,
											   Smile_correl,
											   NULL,
											   true,
											   model->GetCpnInfCurve(),
											   model->GetCpnIRCurve());

	if (itsForcedCorrelationForCdo2) itsIntermediateModel->SetCorrelation(itsForcedCorrelationForCdo2);

	// if (Recovery) delete[] Recovery;

	for (k=0;k<sizePtf;k++) delete DefCurves[k];
	// delete[] DefCurves;
	//FreePointerTabChar(TranchesLabels,sizePtf);

	clearDistribLoss(); 

	if (Smile_correl) delete Smile_correl;
	Smile_correl = NULL;
}


double ICM_Pricer_Analytic_Cdo2_Smile::ComputeImplicitSpread(const string& name,const string& tenor)
{
	if (itsIntermediateModel0 == NULL) { CptTranchesModel();}

	const ICM_DefaultCurve* defcurve = itsIntermediateModel0->GetDefaultCurve(name);
	if (defcurve==NULL)	return (0.);

	double spread = defcurve->ImpliedSpreadInterpol( tenor);

	return (spread);
}


/*----------------------------------------------------------------------------*
  Computing the Tranche's Duration
*----------------------------------------------------------------------------*/ 
double ICM_Pricer_Analytic_Cdo2_Smile::ComputeDuration(void)
{

	ICM_Cds* cds = (ICM_Cds*) GetIntermediateSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();

	double Duration =0.,bp = 0.01,Fair_Spread=0.,spread=0.;
	double FeeLeg = 0.,Accrued = 0.,cds_NPV =0.,init_NPV=0.;

	double notional = cds->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0)*cds->GetTradedCoef();

	if (cds->GetFeeLeg()->GetLegType() != K_FIXED_LEG)
	{
		ICM_Cds* Fixed_Cds = (ICM_Cds*) cds->Clone();
		Fixed_Cds->GetFeeLeg()->SetCreditLegType(qRunning_Leg); //Force Fix Leg
		Fixed_Cds->SetCreditSpread(bp);							//Force Coupon Rate != 0
		Fixed_Cds->GetFeeLeg()->SetRefSpread(1.);				//Force Ref Spread if exist (CM)
		Fixed_Cds->GetFeeLeg()->GetCreditInfos()->SetFwdType(qCPN_FIXED); //Force Fixed coupon
		SetIntermediateSecurity(Fixed_Cds);									//New Security
		ResetRootPricer();										//Reset intermediate pricing parameters		
		SetInitialPriceFlg(false);								//Reset Initial NPV
		cds_NPV = Price(qCMPPRICE);										//Pricing 
		FeeLeg = Price(qCMPFEELEGPV);								//Feeleg Price
		Accrued = Price(qCMPACCRUED);								//accrued pv
		Duration = 100.*(FeeLeg-Accrued)/(notional*bp);				//Duration	
		SetIntermediateSecurity(cds);										//Reset Initial Cds
		ResetRootPricer();										//Reset intermediate pricing parameters		
		SetInitialPriceFlg(false);								//Reset Initial NPV
		if (Fixed_Cds) delete Fixed_Cds;Fixed_Cds=NULL;			//Free Memory
	}
	else //Fixed Leg
	{
		spread = cds->GetCreditSpread();						//Get Initial Spread
		cds_NPV = Price(qCMPPRICE);										//Pricing 
		FeeLeg = Price(qCMPFEELEGPV);								//Feeleg Price
		Accrued = Price(qCMPACCRUED);								//accrued pv
		Fair_Spread = ComputeSpread(0.);						//comopute breakeven spread

		if(FeeLeg) Duration = (FeeLeg-Accrued)/(notional*(100.*spread)/10000.);
		else Duration = cds_NPV/(notional*(100.*spread - Fair_Spread)/10000.);
	}

	SetDuration(fabs(Duration));
	return fabs(Duration) ; 
}


/*----------------------------------------------------------------------------*
  Compute Breakeven Spread
*----------------------------------------------------------------------------*/ 
double ICM_Pricer_Analytic_Cdo2_Smile::ComputeSpread(const double& MtM )
{
	ICM_Cds* OriginalCds = (ICM_Cds*) GetIntermediateSecurity();

	if (GetSpreadFlg()) return GetSpread(); //Already computed
	if (OriginalCds->GetSecurityType() == qCM_TRANCHE) return ComputeImpliedPartRate(); //CM-Cds spread

	ICM_Cds* cds = (ICM_Cds*) OriginalCds->Clone();
	ICM_Leg* Feeleg = (ICM_Leg*) cds->GetFeeLeg();
	qCredit_Leg_Type feelegtype = Feeleg->GetCreditLegType();

	SetIntermediateSecurity(cds);

	double Result = 0.,NPV = 0.,bp = 0.01,InitialSpread=0.,NewNPV=0.;

	if (Feeleg->GetLegType() == K_FIXED_LEG)
	{
		InitialSpread = Feeleg->GetCreditSpread();
		NPV = Price(qCMPPRICE);

		if  (((GetSecurity()->GetName() == ICM_MEZ ) || 
			  (GetSecurity()->GetName() == ICM_NTD ) ||
			  (GetSecurity()->GetName() == ICM_CDO2 )) 
			  &&  (InitialSpread) 
			  &&  (feelegtype != qStyle_None_Leg))
		{	Result = InitialSpread * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
			Result = fabs(Result*100);	
			SetSpread(Result);
			if (cds) delete cds; cds=NULL;
			SetIntermediateSecurity(OriginalCds);
			return (Result);
		}

		ResetRootPricer();

		Feeleg->SetCreditSpread(bp);
		NewNPV = Price(qCMPPRICE);
		Result = bp * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
		Feeleg->SetCreditSpread(InitialSpread);
	}
	else
	{
		ResetRootPricer();

		Feeleg->SetCreditSpread(bp);
		NewNPV = Price(qCMPPRICE);
		Result = bp * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
		Feeleg->SetCreditSpread(InitialSpread);
	}

	ResetRootPricer();
	SetIntermediateSecurity(OriginalCds);
	if (cds) delete cds; cds=NULL;

	Result = fabs(Result*100.);	
	SetSpread(Result);

	return (Result);
}


// virtual 
void 
ICM_Pricer_Analytic_Cdo2_Smile::BeforePrice(const std::string& label  ,qSENSITIVITY_TYPE type /** = ICMSPREAD_TYPE **/ ) 
{
	if ( (label!="NONE") && ((type == ICMSPREAD_TYPE) || (type == ICMRECOVERY_TYPE))) 
	{ 
		// why ICMSPREAD_TYPE ? was a call with default value parameter
		ResetTranchesModel(label,ICMSPREAD_TYPE); }
	else 
	{ ResetTranchesModel("NONE",type);}
} 

void 
ICM_Pricer_Analytic_Cdo2_Smile::Init() 
{
	SetName(ICM_PRICER_ANALYTIC_CDO2_STRIKE);

	itsIntermediateModel = NULL;
	itsIntermediateSecurity = NULL;

	itsIntermediateModel0 = NULL;

	itsImplicitUnderlyingCorrelations.clear();
	itsForcedCorrelationForCdo2 = NULL;

	itsIntegrationStepUnderlying = INTEGRATION_STEP;
	itsIntegrationMethodUnderlying = qGAUSS_HERMITE;
	itsNumSimulationsForCorrelation = 100000;
	itsForcedRescalingForUnderlyings = 0;
	itsAverageCorrTpsDefault_Down = 0.;	
	itsAverageCorrTpsDefault_Up = 0.;	
	itsAverageCorrTpsDefaultFlg = false;
	itsAverageCorrTpsDefault_positive_bump_Down =0.;
	itsAverageCorrTpsDefault_positive_bump_Up =0.;
	itsAverageCorrTpsDefault_negative_bump_Down =0.;
	itsAverageCorrTpsDefault_negative_bump_Up =0.;

	its_Correlation_Fit_Type	=	qCORREL_1F_SINGLE_VALUE;

	//its_Sectorial_Correlation	=	NULL;
	itsCorrelTpsDefaut=0; 
}

// virtual 
void 
ICM_Pricer_Analytic_Cdo2_Smile::Reset(void)
{}

ICM_Pricer_Analytic_Cdo2_Smile::~ICM_Pricer_Analytic_Cdo2_Smile(void)
{ 
	if (itsCorrelTpsDefaut) delete itsCorrelTpsDefaut ; itsCorrelTpsDefaut=0; 

	if (itsIntermediateModel)
		delete itsIntermediateModel;
	itsIntermediateModel = NULL;

	if (itsIntermediateSecurity)
		delete itsIntermediateSecurity;
	itsIntermediateSecurity = NULL;

	if (itsIntermediateModel0)
		delete itsIntermediateModel0;
	itsIntermediateModel0 = NULL;

	if (itsForcedCorrelationForCdo2)
		delete itsForcedCorrelationForCdo2;
	itsForcedCorrelationForCdo2 = NULL;

	if (itsImplicitUnderlyingCorrelations.size()>0)
	{
		for (int i=0;i<itsImplicitUnderlyingCorrelations.size();i++)
		{
			if (itsImplicitUnderlyingCorrelations[i]) delete itsImplicitUnderlyingCorrelations[i];
			itsImplicitUnderlyingCorrelations[i] = NULL;
		}
	}
	itsImplicitUnderlyingCorrelations.clear();
}

// virtual 
void ICM_Pricer_Analytic_Cdo2_Smile::MarketData(ARM_Security* sec,vector<string>& DMKT)
{
	int nbdmk=0;

//			DMKT.resize(nbdmk+1);
//			ARM_Currency*	ccy		=	sec->GetCurrencyUnit();
//			DMKT[nbdmk]= GetZeroCurveName(ccy->GetCcyName(),GetAsOfDate());
	
//			nbdmk++;
	DMKT.resize(nbdmk+1);
	DMKT[nbdmk]= GetDefaultCurveName("SECTORIAL_CORRELATION", GetAsOfDate());			// HARD CODED!!! as a DEFAULT CURVE, bouhhh


	return;
}