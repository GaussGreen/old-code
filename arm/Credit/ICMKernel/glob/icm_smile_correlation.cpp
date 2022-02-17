#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ARMKernel\crv\volint.h"
#include "ICMKernel\inst\icm_mez.h"

#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"
#include "ICMKernel\glob\icm_index_correlation.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include "ICMKernel\crv\icm_distriblosscalculator.h"
#include <nag.h>
#include "ICMKernel\util\icm_RootFinderND.h"
#include "ICMKernel\glob\icm_betas_correlation.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel/crv/icm_volinterpol.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ARMKernel\crv\zeroflat.h"

#ifndef EPSILON_RESCALING
#define EPSILON_RESCALING 1.e-3
#endif


static void __stdcall objfun_DiffSensisEqUp(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);

static void __stdcall objfun_DiffSensisEqDown(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);

// --------------------------------------------------------------------
_Slice_Correl::~_Slice_Correl() 
	{if (itsVolCurve) delete itsVolCurve;
	//if (itsTermVolCurve) delete itsTermVolCurve;
	if (itsIndex) delete itsIndex;
	itsSmileStrikeLow.clear();
	itsSmileStrikeHigh.clear();
	itsMaturities.clear();
}
// --------------------------------------------------------------------
void _Slice_Correl::Init() 
	{itsVolCurve=NULL;
	//itsTermVolCurve=NULL;
	itsIndex=NULL;
	itsSmileStrikeLow.clear();
	itsSmileStrikeHigh.clear();
	itsMaturities.clear();
	}
// --------------------------------------------------------------------
_Slice_Correl::_Slice_Correl(const _Slice_Correl& input) 
	{Init();
	 if (input.itsVolCurve)
	 { 	
		 itsVolCurve = new ICM_VolInterpol(input.itsVolCurve);
		 //itsTermVolCurve = (ICM_VolInterpol*)(input.itsTermVolCurve->Clone());
	 }
	 if (input.itsIndex) itsIndex = (ICM_Credit_Index*) input.itsIndex->Clone();
	itsProportion = input.itsProportion;
	itsSmileStrikeLow = input.itsSmileStrikeLow;
	itsSmileStrikeHigh = input.itsSmileStrikeHigh;
	itsMaturities = input.itsMaturities;}
// --------------------------------------------------------------------
void _Slice_Correl::SetVolCurve(const ARM_VolCurve* volcurve)
{	if (itsVolCurve){delete itsVolCurve;itsVolCurve=NULL;}
	//  test if volcurve is nulle
	if ( volcurve == NULL ) {
		if (itsVolCurve != NULL) {
			delete itsVolCurve;
			itsVolCurve = NULL;
		}
		itsVolCurve=NULL;
	} else {
		itsVolCurve = dynamic_cast<ARM_VolLInterpol*>(((const ARM_Object*)(volcurve))->Clone());
	}
}
// --------------------------------------------------------------------
void _Slice_Correl::SetIndex(const ICM_Credit_Index* index)
{if (itsIndex){delete itsIndex;itsIndex=NULL;}
	 itsIndex = (ICM_Credit_Index*) ((const ARM_Object*)(index))->Clone();
}
// --------------------------------------------------------------------
// inline 
void ICM_Smile_Correlation::Init(void)
{
	ICM_Correlation::Init(); 
	SetName(ICM_SMILE_CORRMATRIX);
	itsForcedStrikeType = qStrike_NONE;

	its_Callib_PtfBespoke=NULL;
	its_Callib_Index=NULL;
	its_Callib_model=NULL;
	its_Callib_K1_PtfBespoke=0.;
	its_Callib_K2_PtfBespoke=0.;
	its_Callib_NoIndex=0;
	its_Callib_SizeEq_Index=0.;
	its_Callib_ActiveMaturity =0.;
	itsSlices.clear();
	its_Already_rescal = false;
	its_never_rescal = false;
	its_TermStructureRescaling = qNoTermStructure;
	its_PricerType = ICM_PRICER_HOMOGENEOUS_SMILE;
	its_Rescaling_Full = false;
	its_Interpolation_At_Maturity=false;

	its_Callib_ELBespoke=1.;
	its_Callib_ELIndex=1.;
	its_NormalizeByEL=false;
	its_IndexJoinOptim=false;
	itsRescalType = qRescal_Std;
	its_BaseCorrelShift=0.;
}
// --------------------------------------------------------------------
	// ----------------------------
	//	Copy of members data
	// ----------------------------
void ICM_Smile_Correlation::BitwiseCopy(const ARM_Object* src)
{
	ICM_Smile_Correlation* CorrMatrix = (ICM_Smile_Correlation *) src;

	int initial_size = CorrMatrix->GetSize();

	itsSlices = CorrMatrix->itsSlices;
	its_Already_rescal = CorrMatrix->its_Already_rescal;
	its_never_rescal = CorrMatrix->its_never_rescal;
	its_PricerType = CorrMatrix->its_PricerType;
	its_TermStructureRescaling = CorrMatrix->its_TermStructureRescaling;
	its_Rescaling_Full = CorrMatrix->its_Rescaling_Full;
	its_NormalizeByEL = CorrMatrix->its_NormalizeByEL;
	its_IndexJoinOptim = CorrMatrix->its_IndexJoinOptim;
	itsRescalType = CorrMatrix->itsRescalType;
	its_BaseCorrelShift = CorrMatrix->its_BaseCorrelShift;
}
// --------------------------------------------------------------------
// virtual 
void ICM_Smile_Correlation::ResetBetas()
{
	int size = itsSlices.size();
	if (its_never_rescal) return;

	for (int i=0; i<size; i++)
	{
		itsSlices[i].GetSmileStrikeLow().clear();
		itsSlices[i].GetSmileStrikeHigh().clear();
		itsSlices[i].GetMaturities().clear();
	}

	its_Already_rescal = false;
}
// --------------------------------------------------------------------
// virtual 
void ICM_Smile_Correlation::SetProportionsInfos(const std::string& indexname,
						 const double& proportion,
						 const double& forcedstrikelow ,
						 const double& forcedstrikehigh )
{
	int i =0;

	i = GetLabelNo(indexname);
	itsSlices[i].SetProportion(proportion);
	
	if ((itsSlices[i].GetSmileStrikeLow().size()>0) && (forcedstrikelow != CREDIT_DEFAULT_VALUE))
		itsSlices[i].SetSmileStrikeLow(forcedstrikelow);
	else if (itsSlices[i].GetSmileStrikeLow().size()==0)
	{
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(forcedstrikelow);
		itsSlices[i].SetSmileStrikeLow(SmileStrikeLow);
	}


	if ((itsSlices[i].GetSmileStrikeHigh().size()>0) && (forcedstrikehigh != CREDIT_DEFAULT_VALUE))
		itsSlices[i].SetSmileStrikeHigh(forcedstrikehigh);
	else if (itsSlices[i].GetSmileStrikeHigh().size()==0)
	{
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(forcedstrikehigh);
		itsSlices[i].SetSmileStrikeHigh(SmileStrikeHigh);
	}

}
// --------------------------------------------------------------------
// inline 
void ICM_Smile_Correlation::Set(const ARM_Date& AsOf,
				const std::string& name,  
				const ARM_VolCurve** VolCurves, 
				const std::vector<std::string>& labels,
				const ARM_Vector* Proportions,
				const ICM_Credit_Index** IndexVector ) 
{
	ICM_Correlation::Set(AsOf,labels,name,(ARM_IRIndex*)0,(ARM_IRIndex*)0);

	// SetName(ICM_SMILE_CORRMATRIX);

	for (int i=0;i<labels.size();i++)
	{
		_Slice_Correl S;
		S.SetProportion((*Proportions)[i]);
		S.SetVolCurve((VolCurves)[i]);
		S.SetIndex((IndexVector)[i]);
		itsSlices.push_back(S);
	}
}


// --------------------------------------------------------------------
//
// virtual 
vector<ARM_VolCurve*> ICM_Smile_Correlation::GetVolCurves()
	{
		vector<ARM_VolCurve*> tabvols;
		tabvols.resize(itsSlices.size());

		for (int i=0;i<itsSlices.size();i++)
		{tabvols[i] = itsSlices[i].GetVolCurve();}

		return tabvols;
	}		
void ICM_Smile_Correlation::Set(const ARM_Date& AsOf,
				const string& name, 
				const ARM_VolCurve** VolCurves, 
				// char** label,
				const std::vector<std::string>& labels,
				const ARM_Vector* Proportions,
				const ARM_Vector* SmileStrikeLow,
				const ARM_Vector* SmileStrikeHigh)
{
	
	ICM_Correlation::Set(AsOf,labels,name,(ARM_IRIndex*)0,(ARM_IRIndex*)0);

	// SetName(ICM_SMILE_CORRMATRIX);

	for (int i=0;i<labels.size();i++)
	{
		_Slice_Correl S;
		
		std::vector<double> SHigh;
		if (SmileStrikeHigh == NULL)
			SHigh.push_back(-1.0);
		else
			SHigh.push_back((*SmileStrikeHigh)[i]);
		S.SetSmileStrikeHigh(SHigh);
		
		std::vector<double> SLow;
		if (SmileStrikeLow == NULL)
			SLow.push_back(-1.0);
		else
			SLow.push_back((*SmileStrikeLow)[i]);
		S.SetSmileStrikeLow(SLow);
		

		S.SetProportion((*Proportions)[i]);
		S.SetVolCurve((VolCurves)[i]);
		//S.SetTermVolCurve((VolCurves)[i]);
		itsSlices.push_back(S);
	}

	its_never_rescal = true;
}


// --------------------------------------------------------------------
//
//		fullStrikeLow (Up) Structure
//
//							col0				col1
//					row0:	yearterm1			yearterm2	
//		index1		row1:	strike low		strike low
//		index2		row2:	strike low		strike low
//
void ICM_Smile_Correlation::Set(const ARM_Date& AsOf,
								const std::string& name, 
								const ARM_VolCurve** VolCurves, 
								const std::vector<std::string>& labels,
								const ARM_Vector& Proportions,
								const ICM_QMatrix<double>& fullStrikeLow,
								const ICM_QMatrix<double>& fullStrikeUp)
{
	ICM_Correlation::Set(AsOf,labels,name,(ARM_IRIndex*)0,(ARM_IRIndex*)0);

	// SetName(ICM_SMILE_CORRMATRIX);

	for (int i=0;i<labels.size();i++)
	{
		_Slice_Correl S;

		S.SetSmileStrikeLow(fullStrikeLow.RowAsStdVector(i+1)); 
		S.SetSmileStrikeHigh(fullStrikeUp.RowAsStdVector(i+1)); 
		S.SetMaturities(fullStrikeLow.RowAsStdVector(0));  // assuming same maturities for up & low. 
		S.SetProportion(Proportions[i]);
		S.SetVolCurve((VolCurves)[i]);
		itsSlices.push_back(S);
	}

	its_never_rescal = true;
}
// --------------------------------------------------------------------
// Calculation of equivalents strikes
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq(ICM_Pricer* pricer, const qRescalType& rescal)
{
	//Set Pricer Name
	SetPricerType(pricer->GetName());

	itsRescalType = rescal;

	const ICM_Parameters& parameters = pricer->GetParameters();
	if (!parameters.empty())
	{	ARM_Vector* VResc = parameters.GetColVect("TERMS_RESCALING");
		if (VResc) 
		{
			qTERM_STRUCTURE type;
			bool ret;
			ICM_EnumsCnv::cnv(VResc->Elt(0),type , ret);
			SetTermStructureRescaling(type);
		}

		VResc = parameters.GetColVect("FULL_RESCALING");
		if (VResc) 
		{
			if (VResc->Elt(0)==1.) 
			{its_Rescaling_Full = true;}
		}
	}

	if ((its_Already_rescal)||(its_never_rescal)) return;

	its_NormalizeByEL = false;
	its_IndexJoinOptim = false;

	//Select Rescaling method
	switch (rescal)
	{
    case qRescal_ING_Maturity :
		{
			its_Interpolation_At_Maturity = true;
			ComputeStrikes_ING(pricer);
			break;
		}
	case qRescal_Eqty_Digital_Maturity : 
		{
			its_Interpolation_At_Maturity = true;
			ComputeStrikesEq_Equity(pricer);
			break;
		}
	case qRescal_Eqty_Mixed_NormByEL_Maturity :
		{	its_IndexJoinOptim = true;	}
	case qRescal_Eqty_NormByEL_Maturity:
		{
			its_NormalizeByEL = true;
			its_Interpolation_At_Maturity = true;

			if (its_IndexJoinOptim)
				ComputeStrikesEq_Equity_CombinIDX(pricer);
			else if (its_TermStructureRescaling==qTermStructure)
				ComputeStrikesEq_Equity_Term(pricer);
			else
				ComputeStrikesEq_Equity(pricer);
			break;
		}
	case qRescal_Eqty_Mixed_Maturity :
		{	its_Interpolation_At_Maturity = true;
			its_IndexJoinOptim = true;
			ComputeStrikesEq_Equity_CombinIDX(pricer);
			break;
		}
	case qRescal_Eqty_Maturity:
		{	its_Interpolation_At_Maturity = true;}
	case qRescal_Eqty:
		{
			if (its_TermStructureRescaling==qTermStructure)
				ComputeStrikesEq_Equity_Term(pricer);
			else
				ComputeStrikesEq_Equity(pricer);
			break;
		}
	case qRescal_Eloss_Maturity: 
		{	its_Interpolation_At_Maturity = true;}
	case qRescal_Eloss:
		{
			if (its_TermStructureRescaling==qTermStructure)
				ComputeStrikesEq_ELoss_Term(pricer);
			else
				ComputeStrikesEq_ELoss(pricer);
			break;
		}
	case qRescal_Std_Maturity:
		{	its_Interpolation_At_Maturity = true;}
	case qRescal_Std:
	default :
		if (its_TermStructureRescaling==qTermStructure) 
			ComputeStrikesEq_Term(pricer);
		else 
			ComputeStrikesEq_standart(pricer);
	}
}

void ICM_Smile_Correlation::ComputeStrikesEq_standart(ICM_Pricer* pricer)
{
	if ((its_Already_rescal)||(its_never_rescal)) return;

	double tolerance_en_0 = 0.;
	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);
	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = K1 + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());
	int size = itsSlices.size();
	int i = 0;

	//Test : tranche equity
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	//Test : tranche [x - 100]
	bool flag100 = false;
	if CHECK_NULL(K2-1.) flag100=true;

	double FixedSensiPtf = 0.,SensiPtf = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_long_nbnames;
	ICM_Collateral*collat=cdo->GetCollateral(); 
	int bespoke_nbnames = collat->GetNbIssuers();
	// If Short Names in collateral they are removed from rescaling computation 
	if (collat->IsLongShort()) 
		bespoke_long_nbnames = bespoke_nbnames - collat->GetNbIssuersShort();
	else
		bespoke_long_nbnames = bespoke_nbnames;
	
	its_Callib_model= (ARM_Model*) pricer->GetModel()->Clone();

	//Borne Inf de la calibration
		//Mezzanine
 		double _inf = 1.e-4;
		//Equity
		double _infeq = 1.e-4;
	
	//Test si resultat correct
	double test_inf=0., test_sup=0.;

	//Borne Sup de la dichotomie : varie en fonction de la maturité (a modifier si existence plot <5Y)
	//Si <3Y on elargi le range	
	double _sup = 0.,_sup2 = 0.;
	if ((Maturity-AsOf)/365. >= 3.)
		_sup = MIN(K2+0.3,0.99);
	else
		_sup = MIN(K2+0.5,0.99);
	
	//Estimation de l'Average Loss du portefeuille
	double result = 0., RecovCoef = 1.;
	double AVG_LR_CDO = 0.,AVG_LR_TRX = 0.;
	double PorS = 0.;

	if (cdo->GetBinaryFlg()) AVG_LR_CDO = 1 - cdo->GetBinary();
	else 
	{
		if (cdo->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = cdo->GetCollateral()->GetRecovCoef();
		for (i=0; i<bespoke_nbnames; i++) 
		{
			// JLA: we don't care of 'AsOf' since we just want to perform 'long or short' .. 
			if (collat->GetIssuersNotional(i,AsOf) > 0)
				AVG_LR_CDO += MIN(((ICM_ModelMultiCurves*)its_Callib_model)->GetRecoveryRate(cdo->GetCollateral()->GetIssuersLabels(i)) * RecovCoef,1);
		}
		AVG_LR_CDO /= (double) bespoke_long_nbnames; AVG_LR_CDO = 1 - AVG_LR_CDO;
	}

	//Les nominaux du portefeuille sont pris en compte dans l'estimation de Strike Equivalents
	//for (i=0; i<cdo->GetNbIssuers(); i++) *(cdo->GetIssuersNotionals()+i) = 10000000.;

	if (!its_Rescaling_Full)
		{cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}

	cdo->SetSubAmount(0.);

	double inf_attachment=0.;

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		// double val_inf, val_sup,time,result = 0.;
		double time,result = 0.;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;
		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		//Borne Inf
		if (idx_inf!=CREDIT_DEFAULT_VALUE) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_inf]);
		
		//Borne Sup : méthode à modifier si besoin
		if ((isequal==false) && (idx_sup!=CREDIT_DEFAULT_VALUE)) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_sup]);

		//only maturity
		if (its_Interpolation_At_Maturity)
		{
			itsSlices[NoIndex].GetMaturities().clear();
			itsSlices[NoIndex].GetMaturities().push_back(time);
		}


		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{
			result = 0.;
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
			its_Callib_PtfBespoke= cdo;
			its_Callib_K1_PtfBespoke= K1;
			its_Callib_K2_PtfBespoke= K2;
			its_Callib_NoIndex=NoIndex;
			its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
			PorS = cdo->GetPorS();

		
			int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
			for (int il_corr=0;il_corr<szcorrel;il_corr++)
			{inf_attachment = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100.;
			if (!CHECK_NULL(inf_attachment)) break;}

			//Maturité de l'indice de référence
			ARM_Date Maturity_2 = AsOf;
			Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
			its_Callib_Index = GenerateCdoIndex(Start,Maturity_2,NoIndex,PorS);

			//Estimation de la size equivalente (Average Loss de l'indice)
			AVG_LR_TRX = 0.;
			for (i=0; i<((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers(); i++) AVG_LR_TRX += ((ICM_ModelMultiCurves*)its_Callib_model)->GetDefaultCurve(((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetIssuersLabels(i))->GetRecovery();	
			AVG_LR_TRX /= (double) ((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers(); AVG_LR_TRX = 1 - AVG_LR_TRX;

			its_Callib_SizeEq_Index=(AVG_LR_TRX/AVG_LR_CDO)*(K2-K1)*
									(((double)bespoke_long_nbnames)/((double)((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers()));

			_sup2 = _sup - its_Callib_SizeEq_Index;

			try
			{
				//Tranche Equity
				if (flagequity)
				{
					//Test : existence de la solution
					test_inf = DiffSensisEquity(_infeq);
					test_sup = DiffSensisEquity(_sup);

					//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
					if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
						result = _infeq;
					//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
					else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
						result = _sup;
					//Cas Std : minimisation
					else
						//Pas de test car test en 0 -> size = 0 et prix Equity = 0
						result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEquity,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
				}
				//Tranche Mezzanine
				else
				{
					//Test : tranche [x-100] : rescaling de la tranche equity corespondante
					if (flag100)
					{							
						//Test : existence de la solution
						test_inf = DiffSensisEqDown(_infeq);
						test_sup = DiffSensisEqDown(_sup);

						//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
						if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
							result = _inf;
						//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
						else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
							result = _sup;
						//Cas Std : minimisation
						else
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDown,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
					}
					else
					{
						//Test existence de la solution
						test_inf = DiffSensis(_inf);
						test_sup = DiffSensis(_sup2);

						//Cas 1 : ELoss Bespoke > ELoss Traxx -> détachement par défaut à _inf
						if (test_inf < tolerance_en_0)
							result = _inf;
						//Cas 2 : ELoss Bespoke < ELoss Traxx -> On ne peut garder la meme size equivalente
						//typiquement pour les tranches mezzanines x - 100%
						else if (test_inf * test_sup > 0.)
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensis,(*this))).Dichotomy(_inf,_sup,100,1.E-4,1.E-4);
						//Cas std
						else
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensis,(*this))).Dichotomy(_inf,_sup2,100,1.E-4,1.E-4);
					}
				}
			}
			catch (...)
			{
				if (cdo) delete cdo;
				if (its_Callib_model) delete its_Callib_model;
				if (its_Callib_Index) delete its_Callib_Index;
				ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
				itsSlices.clear();
			}

			if (its_Callib_Index) delete its_Callib_Index;
			}

			if ((flagequity) && (result<inf_attachment)) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(inf_attachment);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(inf_attachment);
			}
			else if (flagequity) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result);
			}
			else if (flag100) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(1.);
			}
			else
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(MIN(result+its_Callib_SizeEq_Index,1.0));
			}

		}
	}

	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;

}


// --------------------------------------------------------------------
// Calculation of equivalents strikes
// --------------------------------------------------------------------
ICM_Mez* ICM_Smile_Correlation::GenerateCdoIndex(ARM_Date& Start,
												 ARM_Date& Maturity,
												 const int& NoIndex,
												 const double& PorS)
{
	ICM_Credit_Index* Index = itsSlices[NoIndex].GetIndex();
	int size = Index->GetNbCurves();
	std::vector<std::string> label(size); 
	for(int i=0;i<size;i++) label[i]= Index->GetLabel(i); 

	std::vector<double> IssuersNotionals(size); 
	for ( i=0; i<size; i++) {IssuersNotionals[i] = 10000000.;}
	double FixedRate = 0.;

	ICM_Mez* cdo = new ICM_Mez(Start,
						   Maturity,
						   (ARM_Date*)0,
						   (ARM_Date*)0,
						   FixedRate,
						   K_MATUNADJUSTED, 
						   K_UNADJUSTED, 
						   0.,
						   10000000.,
						   label,
						   IssuersNotionals,
						   // size,
						   K_QUARTERLY,
						   KACTUAL_360,
						   100., 
						   qACCRUED_SETTLED,
						   string(Index->GetCurrencyUnit()->GetCcyName()),//ARM_DEFAULT_COUNTRY, 
						   0.0,
						   K_SHORTSTART,
						   DEFAULT_CREDIT_LAG_INDX,
							DEFAULT_FRQ_DEFLEG, // const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/,
							CREDIT_DEFAULT_VALUE,// const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
							std::string(),// const std::string&  payCalName /* = NULL*/,
							qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
							qStandart_Recovery_Leg,// const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
							INCLUDE_MATURITY// const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);
						   ); 

	if (!its_Rescaling_Full)
	{ cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}

	cdo->SetSubAmount(0.);

	// if (IssuersNotionals) delete[] IssuersNotionals;

	return (cdo);
}

// --------------------------------------------------------------------
// Calculation of Diff sensitivities
// --------------------------------------------------------------------
double ICM_Smile_Correlation::DiffSensis(double Strike_down)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	double DefPv_Index_K1 = 0.,DefPv_bespoke_K1 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;

	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();;
	
	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA(); 
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetDefLeg()->GetMaturity();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;
	
	//On cap le point de détéchament max à 1
	double K2_index = MIN(Strike_down + its_Callib_SizeEq_Index,1);
	double K1_index = Strike_down;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	double Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
	double Beta_index_i_k1 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K1_index);

	Beta_index_i_k2 /= 100.;
	Beta_index_i_k1 /= 100.;

	//Equity tranche Up
	ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(&Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx; pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);
	
	//Bespoke Price (fonction du type de pricer collat fwd ou non)
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Equity Tranche Down
	traxx->SetMezzAmount(K1_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K1_PtfBespoke*SumNot_bespoke);

	ICM_Beta_Correlation Correl_k1(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k1),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(&Correl_k1);
	
	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K1_traxx;pricer_K1_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_traxx.SetFaster(true);pricer_K1_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K1 = pricer_K1_traxx.Price(qCMPDEFLEGPV);
	
	//Bespoke Price (fonction du type de pricer collat fwd ou non)
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result = (DefPv_Index_K2-DefPv_Index_K1)/(SumNot_traxx*(K2_index-K1_index));
    result -= (DefPv_bespoke_K2-DefPv_bespoke_K1)/(SumNot_bespoke*(its_Callib_K2_PtfBespoke-its_Callib_K1_PtfBespoke))/df_bespoke*df_index;

	return (result);
}


// --------------------------------------------------------------------
// Calculation of Diff sensitivities
// --------------------------------------------------------------------
double ICM_Smile_Correlation::DiffSensisEquity(double Strike_up)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
	ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_up;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	//Update de la correl
	if (K2_index <= last_strike)
		//Cas std : on prend la correl système
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
	else
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*last_strike);

	Beta_index_i_k2 /= 100.;
	
	ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(&Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
		switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}

	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result = (DefPv_Index_K2)/(SumNot_traxx*Strike_up);
	result -= (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K2_PtfBespoke)/df_bespoke*df_index;

	return (result);
}


// --------------------------------------------------------------------
// View Mathod
// --------------------------------------------------------------------
void ICM_Smile_Correlation::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

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

	int size = GetSize();

	fprintf(fOut, "\n ======> Parameters of pricing for correlation:\n");
	if (its_never_rescal) fprintf(fOut, "\tRescal : NO\n"); else  fprintf(fOut, "\tRescal : YES\n");
	if (its_Already_rescal) fprintf(fOut, "\tAlready Rescal : YES\n"); else fprintf(fOut, "\tAlready Rescal : NO\n"); 
	if (its_TermStructureRescaling==qTermStructure) fprintf(fOut, "\tTermStructure Rescaling : YES\n"); else fprintf(fOut, "\tTermStructure Rescaling : NO\n"); 
	if (its_Rescaling_Full) fprintf(fOut, "\tRescaling Full : YES\n"); else fprintf(fOut, "\tRescaling Full : NO\n"); 
	if (its_IndexJoinOptim) fprintf(fOut, "\tIndex Join Optim : YES\n"); else fprintf(fOut, "\tIndex Join Optim : NO\n"); 
	if (its_Interpolation_At_Maturity) fprintf(fOut, "\tInterpolation At Maturity : YES\n"); else fprintf(fOut, "\tInterpolation At Maturity : NO\n"); 
	if (its_NormalizeByEL) fprintf(fOut, "\tNormalize by EL : YES\n\n"); else fprintf(fOut, "\tNormalize by EL : NO\n\n"); 

	fprintf(fOut, "\n ======> Correlations by Strike/Currency/Maturities/Index :\n\n");
	fprintf(fOut, "\n Maturity\t\tProportions\t\tSmileStrikeLow\t\tSmileStrikeHigh\n");	

	for (i = 0; i<size; i++)
	{
		for (int m = 0; m<itsSlices[i].GetMaturities().size(); m++)
		{
		fprintf(fOut, "%f\t\t",itsSlices[i].GetMaturities()[m]);	
		fprintf(fOut, "%f\t\t",itsSlices[i].GetProportion());	
		fprintf(fOut, "%f\t\t",itsSlices[i].GetSmileStrikeLow()[m]);	
		fprintf(fOut, "%f\n",itsSlices[i].GetSmileStrikeHigh()[m]);
		}

		if (itsSlices[i].GetMaturities().size()==0)
		{
		fprintf(fOut, "%f\t\t",999.);	
		fprintf(fOut, "%f\t\t",itsSlices[i].GetProportion());	
		if (itsSlices[i].GetSmileStrikeLow().size()>0) fprintf(fOut, "%f\t\t",itsSlices[i].GetSmileStrikeLow()[0]);	else fprintf(fOut, "%f\t\t",-999.);
		if (itsSlices[i].GetSmileStrikeHigh().size()>0) fprintf(fOut, "%f\t\t\n",itsSlices[i].GetSmileStrikeHigh()[0]);	else fprintf(fOut, "%f\t\t\n",-999.);
		}
	}

	if (itsSlices[0].GetVolCurve())
	{
		int k =0;

		for (i = 0; i<size; i++)
		{
			if (itsSlices[i].GetVolCurve()) itsSlices[i].GetVolCurve()->View(id,fOut); 
			if (itsSlices[i].GetIndex())
			{ if (itsSlices[i].GetIndex()) itsSlices[i].GetIndex()->View(id,fOut);} 
			fprintf(fOut, "\n");
		}

	}

	ICM_Correlation::View(id,fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


//Equivalent Strike Down
double ICM_Smile_Correlation::GetEqStrikeDown(const std::string& indexname) 
{
	//if (its_Already_rescal==false) {return CREDIT_DEFAULT_VALUE;}
	int size = itsSlices.size();

	for (int i=0; i<size; i++)
	{
		if (itsSlices[i].GetIndex()==NULL)
			return itsSlices[i].GetSmileStrikeLow()[0];
		// else if (strcmp(indexname,itsSlices[i].GetIndex()->GetLabels()[0])==NULL)
		else if (indexname==itsSlices[i].GetIndex()->GetLabels()[0]) 
			return itsSlices[i].GetSmileStrikeLow()[0];
	}

	return CREDIT_DEFAULT_VALUE;

}

//Equivalent Strike Up
double ICM_Smile_Correlation::GetEqStrikeUp(const std::string& indexname) 
{
	//if (its_Already_rescal==false) {return CREDIT_DEFAULT_VALUE;}
	int size = itsSlices.size();

	for (int i=0; i<size; i++)
	{
		if (itsSlices[i].GetIndex()==NULL)
			return itsSlices[i].GetSmileStrikeLow()[0];
		// if (strcmp(indexname,itsSlices[i].GetIndex()->GetLabels()[0])==NULL)
		if (indexname==itsSlices[i].GetIndex()->GetLabels()[0])
			return itsSlices[i].GetSmileStrikeHigh()[0];
	}

	return CREDIT_DEFAULT_VALUE;
}

void ICM_Smile_Correlation::GetSmileStrikeUp(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes) 
{
	int size = itsSlices.size();

	for (int i=0; i<size; i++)
	{
		if (itsSlices[i].GetIndex()==NULL)
		{vStrikes = itsSlices[i].GetSmileStrikeHigh();
		 vMatu = itsSlices[i].GetMaturities();}
		// else if (strcmp(indexname,itsSlices[i].GetIndex()->GetLabels()[0])==NULL)
		else if (indexname==itsSlices[i].GetIndex()->GetCreditIndexName())
		{vStrikes = itsSlices[i].GetSmileStrikeHigh();
		 vMatu = itsSlices[i].GetMaturities();}
	}

}

void ICM_Smile_Correlation::GetSmileStrikeDown(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes) 
{
	int size = itsSlices.size();

	for (int i=0; i<size; i++)
	{
		if (itsSlices[i].GetIndex()==NULL)
		{vStrikes = itsSlices[i].GetSmileStrikeLow();
		 vMatu = itsSlices[i].GetMaturities(); 
		} 
		// else if (strcmp(indexname,itsSlices[i].GetIndex()->GetLabels()[0])==NULL)
		else if (indexname==itsSlices[i].GetIndex()->GetCreditIndexName())
		{vStrikes = itsSlices[i].GetSmileStrikeLow();
		 vMatu = itsSlices[i].GetMaturities(); 
		}			
	}

}


//Equivalent Correl Strike Down
double ICM_Smile_Correlation::GetCorrelStrikeDown(double maturity) 
{
	//for ING rescaling only
	if (itsRescalType==qRescal_ING_Maturity)
	{
		qCorrel_By_Strike StrikeTypeOld = itsForcedStrikeType;
		itsForcedStrikeType=qStrike_LOW;
		double corr = its_BaseCorrelShift + GetCompositeCorrel("ISSUER","ISSUER",maturity,CREDIT_DEFAULT_VALUE,maturity);
		itsForcedStrikeType = StrikeTypeOld;
		return corr;
	}

	double Correlation=0; 
	if ((its_never_rescal)||(its_TermStructureRescaling==qTermStructure))
	{
		int size_index = itsSlices.size();
		for (int i=0; i<size_index; i++)
		if (itsSlices[i].GetProportion()) 
		{
			double Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(maturity,100.*itsSlices[i].GetSmileStrikeLow()[0])/100.;
			Correlation += itsSlices[i].GetProportion()*Inter_Corr;
		}
	}
	else 
	{
		int size_index = itsSlices.size();
		int size_matu = itsSlices[0].GetMaturities().size();
		double Inter_Corr = 0.,interpol=0.;
		vector<double> VCorrelation;

		for (int i=0; i<size_index; i++)
		if (itsSlices[i].GetProportion()) 
			{VCorrelation.clear();
			for (int j=0; j<size_matu; j++)
				{Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(itsSlices[i].GetMaturities()[j],100.*itsSlices[i].GetSmileStrikeLow()[j])/100.;
				 VCorrelation.push_back(Inter_Corr);}
			// Interpolation temporelle, toujours linéaire
			interpol = VectorInterpol(itsSlices[i].GetMaturities(),VCorrelation,maturity,K_LINEAR);
			Correlation += itsSlices[i].GetProportion()*interpol;}
	}
	return (Correlation);
}

//Equivalent Correl Strike Up
double ICM_Smile_Correlation::GetCorrelStrikeUp(double maturity) 
{
	//for ING rescaling only
	if (itsRescalType==qRescal_ING_Maturity)
	{
		qCorrel_By_Strike StrikeTypeOld = itsForcedStrikeType;
		itsForcedStrikeType=qStrike_UP;
		double corr = its_BaseCorrelShift + GetCompositeCorrel("ISSUER","ISSUER",maturity,CREDIT_DEFAULT_VALUE,maturity);
		itsForcedStrikeType = StrikeTypeOld;
		return corr;
	}

	double Correlation=0; 
	//if (its_Already_rescal==false) {return CREDIT_DEFAULT_VALUE;}
	if ((its_never_rescal)||(its_TermStructureRescaling==qTermStructure))
	{
		int size_index = itsSlices.size();
		for (int i=0; i<size_index; i++)
		if (itsSlices[i].GetProportion()) 
		{
			double Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(maturity,100.*itsSlices[i].GetSmileStrikeHigh()[0])/100.;
			Correlation += itsSlices[i].GetProportion()*Inter_Corr;
		}
	}
	else 
	{
		int size_index = itsSlices.size();
		int size_matu = itsSlices[0].GetMaturities().size();
		double Inter_Corr = 0.,interpol=0.;
		vector<double> VCorrelation;

		for (int i=0; i<size_index; i++)
		if (itsSlices[i].GetProportion()) 
			{VCorrelation.clear();
			for (int j=0; j<size_matu; j++)
				{Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(itsSlices[i].GetMaturities()[j],100.*itsSlices[i].GetSmileStrikeHigh()[j])/100.;
				 VCorrelation.push_back(Inter_Corr);
			}

			// Interpolation temporelle, toujours linéaire
			interpol = VectorInterpol(itsSlices[i].GetMaturities(),VCorrelation,maturity,K_LINEAR);
			Correlation += itsSlices[i].GetProportion()*interpol;
		}
	}
	return (Correlation);
}


void ICM_Smile_Correlation::Find(const std::string& info,int& nthLine,int& nthCol,int& index) 
{
	index = nthLine = nthCol = -1;

	//On a une structure ISSUER_STRIKE_MATU

	int  separate = '|';
	char Tmplabel[255];memset(Tmplabel,'\0',sizeof(Tmplabel));
	char TmpStrike[255];memset(TmpStrike,'\0',sizeof(TmpStrike));
	char TmpMatu[255];memset(TmpMatu,'\0',sizeof(TmpMatu));
	strcpy(Tmplabel,info.c_str());
	char* pdest = strchr(Tmplabel, separate );
	if (pdest == NULL) return;
	*pdest = '\0';
	strcpy(TmpStrike,pdest+1);
	pdest = strchr(TmpStrike, separate );
	if (pdest == NULL) return;
	*pdest = '\0';
	strcpy(TmpMatu,pdest+1);

	index = GetLabelNo(Tmplabel);
	ARM_VolLInterpol* volcurve = (ARM_VolLInterpol*) itsSlices[index].GetVolCurve();

	ARM_Vector* Strikes = volcurve->GetStrikes();
	ARM_CRV_TERMS* ExpiryTerms = &(volcurve->itsYearTermsX);
	ARM_Matrix* Volatilities = volcurve->GetVolatilities();
	int nLines = Volatilities->GetNumLines();
    int nCols  = Volatilities->GetNumCols();	

	double Strike = atof(TmpStrike);

	for (int i=0; i<nCols; i++)
	{
		if (fabs(Strike-Strikes->Elt(i))<DB_TOL)
		{
			nthCol = i;
			for (int j=0; j<nLines; j++)
			{
			if (strcmp(TmpMatu,(*ExpiryTerms)[j])==NULL)
			{
				nthLine = j;
				return;
			}
			}
		}
	}
}

// --------------------------------------------------------------------
// Generate Fixed Base Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* FixedBaseCorrelationMult(const ARM_Date& AsOf,
											vector<ARM_Currency*> Ccy,
											const vector<double>& strikeDown,
											const vector<double>& strikeUp,
											const vector<double>& CorrelDown,
											const vector<double>& CorrelUp,
											const vector<string>& IndexName,
											const vector<double>& Proportion)
{
	int size = strikeDown.size();
	// const ARM_VolCurve** VVOL = new const ARM_VolCurve*[size];
	std::vector<const ARM_VolCurve*> VVOL (size); 
	// char** label = new char*[size];
	std::vector<std::string> labels(size); 
	ARM_Vector* propo = new ARM_Vector(size,0.);
	ARM_Vector* strikelow = new ARM_Vector(size,0.);
	ARM_Vector* strikehigh = new ARM_Vector(size,0.);

	for (int i=0; i<size;i++)
	{
	ARM_Vector* YT = new ARM_Vector(1,5.);
	ARM_Vector* Strikes = new ARM_Vector(2,0.);
	Strikes->Elt(0) = 100.*strikeDown[i];
	Strikes->Elt(1) = 100.*strikeUp[i];
	ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());
	mVolatility->Elt(0,0) = 100.*CorrelDown[i];
	mVolatility->Elt(0,1) = 100.*CorrelUp[i];
	propo->Elt(i)= Proportion[i];
	strikelow->Elt(i) = strikeDown[i];
	strikehigh->Elt(i) = strikeUp[i];

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,Ccy[i]);	


	VVOL[i]=VOL;
	labels[i] = IndexName[i] ;
	if (YT) delete YT;

	}

	string name = "FixedCorr";

// FIXMEFRED: mig.vc8 (28/05/2007 10:28:50):cast
	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
//	delete[] label;

	for (int j=0; i<size;i++){delete VVOL[j]; VVOL[j] = NULL;}
	//delete[] VVOL;

	return (ForcedCorrelation);
}


// --------------------------------------------------------------------
// Generate Fixed Base Correlation
// --------------------------------------------------------------------
void FixedBaseCorrelationMult_Vector(const ARM_Date& AsOf,
											const vector<ARM_Currency*>& Ccy,
											const vector<double>& strikeDown,
											const vector<double>& strikeUp,
											const vector<double>& CorrelDown,
											const vector<double>& CorrelUp,
											const vector<string>& IndexName,
											vector<ICM_Smile_Correlation*>& Correlation)
{
	int size = strikeDown.size();
	Correlation.resize(size);

	for (int i=0; i<size; i++)
		Correlation[i] = FixedBaseCorrelation(AsOf,Ccy[i],strikeDown[i],strikeUp[i],
											CorrelDown[i],CorrelUp[i],IndexName[i]);

}



// --------------------------------------------------------------------
// Equivalent Strikes Computation Using Equity stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_Equity(ICM_Pricer* pricer)
{
	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	// pour la normalisation par l'EL du ptf bespoke
	if (its_NormalizeByEL)
	{its_Callib_ELBespoke=CptPtf_ELoss(pricer,K_ZEROCOUPON);}

	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = K1 + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());
	int size = itsSlices.size();
	int i = 0;
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double FixedSensiPtf = 0.,SensiPtf = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_nbnames = cdo->GetCollateral()->GetNbIssuers();

	its_Callib_model= (ARM_Model*) pricer->GetModel()->Clone();

	//Borne inf de la dichotomie
 	double _inf = 0., _infeq=1.e-4;
	
	//Borne Sup de la dichotomie : varie en fonction de la maturité (a modifier si existence plot <5Y)
	double _sup = 0.,_sup2 = 0.;
	_sup = MIN(K2+0.5,0.99);

	//Test resultat numerique
	double test_inf=0., test_sup = 0.;

	//Estimation de l'Average Loss du portefeuille
	double result = 0.;
	double AVG_LR_CDO = 0.,AVG_LR_TRX = 0.;
	double PorS = 0.;

	if (!its_Rescaling_Full)
	{ cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}
	
	double inf_attachment=0.;

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		// double val_inf, val_sup,time,result = 0.;
		double time,result_up = 0., result_down=0.;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;
		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		//Borne Inf
		if (idx_inf!=CREDIT_DEFAULT_VALUE) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_inf]);

		//Borne Sup
		if ((isequal==false) && (idx_sup!=CREDIT_DEFAULT_VALUE)) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_sup]);

		//only maturity
		if (its_Interpolation_At_Maturity)
		{
			itsSlices[NoIndex].GetMaturities().clear();
			itsSlices[NoIndex].GetMaturities().push_back(time);
		}
		
		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{
			result = 0.;
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
				its_Callib_PtfBespoke= cdo;
				its_Callib_K1_PtfBespoke= K1;
				its_Callib_K2_PtfBespoke= K2;
				its_Callib_NoIndex=NoIndex;
				its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
				PorS = cdo->GetPorS();

				int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
				for (int il_corr=0;il_corr<szcorrel;il_corr++)
				{inf_attachment = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100.;
				if (!CHECK_NULL(inf_attachment)) break;}

				//Maturité de l'indice de référence
				ARM_Date Maturity_2 = AsOf;
				Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
				its_Callib_Index = GenerateCdoIndex(Start,Maturity_2,NoIndex,PorS);

				//calcul de l'El d'indice pour la normalisation
				if (its_NormalizeByEL)
					{its_Callib_ELIndex = CptPtf_ELoss_Index((ICM_ModelMultiCurves *)its_Callib_model,
													itsSlices[NoIndex].GetIndex(),Maturity_2,K_ZEROCOUPON);}
				else
					{its_Callib_ELIndex = 1.;}

				try
				{
					if (flagequity)
					{
						//Test : existence de la solution
						if (itsRescalType==qRescal_Eqty_Digital_Maturity){
							test_inf = DiffSensisEqUp_digit(_infeq);
							test_sup = DiffSensisEqUp_digit(_sup);}
						else {
							test_inf = DiffSensisEqUp(_infeq);
							test_sup = DiffSensisEqUp(_sup);}

						if (fabs(test_sup)<=1.E-4)
							result_up = _sup;
						//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
						else if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
							result_up = _inf;
						//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
						else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
							result_up = _sup;
						//Cas Std : minimisation
						else
						{
							if (itsRescalType==qRescal_Eqty_Digital_Maturity)
							{result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp_digit,(*this))).Dichotomy(_infeq,_sup,100,1.E-6,1.E-4);}
							else
							{result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);}
						}
					}
					else
					{
						//Tranche Up
							_sup = MIN(K2+0.5,0.99);
							//Test : existence de la solution
						if (itsRescalType==qRescal_Eqty_Digital_Maturity){
							test_inf = DiffSensisEqUp_digit(_infeq);
							test_sup = DiffSensisEqUp_digit(_sup);}
						else {
							test_inf = DiffSensisEqUp(_infeq);
							test_sup = DiffSensisEqUp(_sup);}

							if (fabs(test_sup)<=1.E-4)
								result_up = _sup;
							//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
							else if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
								result_up = _inf;
							//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
							else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
								result_up = _sup;
							//Cas Std : minimisation
							else{
								if (itsRescalType==qRescal_Eqty_Digital_Maturity){
									result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp_digit,(*this))).Dichotomy(_infeq,_sup,100,1.E-6,1.E-4);}
								else {
									result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);}
							}
						
						//Tranche Down
							_sup = MIN(K1+0.5,0.99);
							//Test : existence de la solution
							if (itsRescalType==qRescal_Eqty_Digital_Maturity){
							test_inf = DiffSensisEqDown_digit(_infeq);
							test_sup = DiffSensisEqDown_digit(_sup);}
							else{
							test_inf = DiffSensisEqDown(_infeq);
							test_sup = DiffSensisEqDown(_sup);}

							//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
							if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
								result_down = _inf;
							//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
							else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
								result_down = _sup;
							//Cas Std : minimisation
							else{
								if (itsRescalType==qRescal_Eqty_Digital_Maturity){
									result_down = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDown_digit,(*this))).Dichotomy(_infeq,_sup,100,1.E-6,1.E-4);}
								else {
									result_down = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDown,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);}
							}
					}
				}
				catch (...)
				{
					if (cdo) delete cdo;
					if (its_Callib_model) delete its_Callib_model;
					if (its_Callib_Index) delete its_Callib_Index;
					ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
					itsSlices.clear();
				}
				if (its_Callib_Index) delete its_Callib_Index;
				
			}

			if (flagequity) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
			else
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result_down);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
		}
	}

	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}

double ICM_Smile_Correlation::DiffSensisEqDown(double Strike_down)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    
	//Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_down;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//MaJ : caractéristiques de la tranche down
	double init_sub = bespoke->GetSubAmount(bespoke ->GetFeeLeg()->GetStartDate());
	double init_mezz= bespoke->GetMezzAmount(bespoke ->GetFeeLeg()->GetStartDate());	

	bespoke->SetMezzAmount(its_Callib_K1_PtfBespoke*SumNot_bespoke);
	bespoke->SetSubAmount(0.);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	//Update de la correl
	if (K2_index <= last_strike)
		//Cas std : on prend la correl système
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
	else
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*last_strike);

	Beta_index_i_k2 /= 100.;
	
	ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(&Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	///Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex/Strike_down) : 1.);
	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke/its_Callib_K1_PtfBespoke):1.);

//	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex) : 1.);
//	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke):1.);

    //On ramène l'EL ZC à la maturité de la bespoke
    double result1 = (DefPv_Index_K2)/(SumNot_traxx*Strike_down* Callib_ELIndex );
//	if ((its_Callib_ELIndex!=1.) && (result1>0.99))
//		{result1 = 1.;}

	double result2 = (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K1_PtfBespoke* Callib_ELBespoke);
//	if ((its_Callib_ELBespoke!=1.) && (result2>0.99))
//		{result2 = 1.;}

	double result = result1;
	if (its_Callib_ELBespoke!=1.)
	{result -= result2;}
	else
	{result -= result2/df_bespoke*df_index;}

	//Remise a jour carac mezz
	bespoke->SetMezzAmount(init_mezz);
	bespoke->SetSubAmount(init_sub);

	return (result);
}
	
double ICM_Smile_Correlation::DiffSensisEqUp(double Strike_up)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_up;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub = bespoke->GetSubAmount(bespoke ->GetFeeLeg()->GetStartDate());
	if (!CHECK_NULL(bespoke->GetSubAmount(bespoke ->GetFeeLeg()->GetStartDate())))
		bespoke->SetSubAmount(0.);
		
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	//Update de la correl
	if (K2_index <= last_strike)
		//Cas std : on prend la correl système
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
	else
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*last_strike);

	Beta_index_i_k2 /= 100.;
	
	ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(&Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex/Strike_up) : 1.);
	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke/its_Callib_K2_PtfBespoke):1.);

    //On ramène l'EL ZC à la maturité de la bespoke
    double result1 = (DefPv_Index_K2)/(SumNot_traxx*Strike_up*Callib_ELIndex);
	if ((its_Callib_ELIndex!=1.) && (result1>0.99) && (Strike_up>=0.99))
		{result1 = 1.;}

	double result2 = (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K2_PtfBespoke*Callib_ELBespoke);
	if ((its_Callib_ELBespoke!=1.) && (result2>0.99) && (Strike_up>=0.99))
		{result2 = 1.;}

	double result = result1;
	if (its_Callib_ELBespoke!=1.)
	{result -= result2;}
	else
	{result -= result2/df_bespoke*df_index;}

	//Remise a jour sub amount init
	bespoke->SetSubAmount(init_sub);

	return (result);
}


// --------------------------------------------------------------------
// Equivalent Strikes Computation Using ELoss stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_ELoss(ICM_Pricer* pricer)
{
	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	int size = itsSlices.size();
	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double eloss_bespoke = 0.,eloss_indice = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();

	//Eloss du portefeuille : sommes des eloss ind
	its_Callib_model= (ICM_ModelMultiCurves*) pricer->GetModel()->Clone();
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	//Nominal total du portefeuille
	//nominal_total = cdo->GetCollateral()->SumNotionals();
	
	//calcul de l'EL du portefeuille
    eloss_bespoke = CptPtf_ELoss(pricer);

					
	//Strike Bespoke -> NB Eloss bespoke
	double StrikeUpEloss = 0.;
	double StrikeDownEloss = 0.;
	if (!CHECK_NULL(eloss_bespoke))
	{
		StrikeUpEloss = (cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate()) + cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate())) / fabs(eloss_bespoke);
		StrikeDownEloss = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate()) / fabs(eloss_bespoke);
	}

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		double time,result_up = 0., result_down=0.;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;

		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		//Determination des bornes de projection en maturité
		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		//Borne Inf
		if (idx_inf!=CREDIT_DEFAULT_VALUE) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_inf]);

		//Borne Sup
		if ((isequal==false) && (idx_sup!=CREDIT_DEFAULT_VALUE)) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_sup]);

		//only maturity
		if (its_Interpolation_At_Maturity)
		{
			itsSlices[NoIndex].GetMaturities().clear();
			itsSlices[NoIndex].GetMaturities().push_back(time);
		}
		
		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{			
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
				//Maturité de l'indice de référence
				its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];

				ARM_Date Maturity_2 = AsOf;
				Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));

		
				//Estimation de l'eloss sur l'indice (exprimée en %)	
				eloss_indice = CptPtf_ELoss_Index(mod,itsSlices[NoIndex].GetIndex(),Maturity_2);
				
				//Strike Indice -> NB ELoss Indice
				int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
				vector<double> StrikeElossIndex(szcorrel);
				for (int il_corr=0;il_corr<szcorrel;il_corr++)
					StrikeElossIndex[il_corr] = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100./eloss_indice;

				ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes();
				vector<double> StrikeIndex(V->GetSize());
				for (int ij=0;ij<V->GetSize();ij++) {StrikeIndex[ij]=V->Elt(ij);}
				
				//Strike Equivalent	
				try
				{
					if (flagequity)
						result_up = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeUpEloss)/100;
					else
					{
						//Tranche Up
						result_up = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeUpEloss)/100;
						
						//Tranche Down
						result_down = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeDownEloss)/100;
					}
				}
				catch (...)
				{
					if (cdo) delete cdo;
					if (its_Callib_model) delete its_Callib_model;
					if (its_Callib_Index) delete its_Callib_Index;
					ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
					itsSlices.clear();
				}
				if (its_Callib_Index) delete its_Callib_Index;

			}

			if (flagequity) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
			else
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result_down);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
		}
	}
	//Traces
	//fclose(stream);
	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}

//----------------------------------------------------------
//	Rescaling resolution by terms 
//----------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_Term(ICM_Pricer* pricer)
{
	if ((its_Already_rescal)||(its_never_rescal)) return;

	double tolerance_en_0 = 0.;
	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);
	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = K1 + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());
	int size = itsSlices.size();
	int i = 0;

	//Test : tranche equity
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	//Test : tranche [x - 100]
	bool flag100 = false;
	if CHECK_NULL(K2-1.) flag100=true;

	double FixedSensiPtf = 0.,SensiPtf = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_nbnames = cdo->GetCollateral()->GetNbIssuers();

	its_Callib_model= (ARM_Model*) pricer->GetModel()->Clone();

	//Borne Inf de la calibration
		//Mezzanine
 		double _inf = 1.e-4;
		//Equity
		double _infeq = 1.e-4;
	
	//Test si resultat correct
	double test_inf=0., test_sup=0.;

	//Borne Sup de la dichotomie : varie en fonction de la maturité (a modifier si existence plot <5Y)
	//Si <3Y on elargi le range	
	double _sup = 0.,_sup2 = 0.;
	if ((Maturity-AsOf)/365. >= 3.)
		_sup = MIN(K2+0.3,0.99);
	else
		_sup = MIN(K2+0.5,0.99);
	
	//Estimation de l'Average Loss du portefeuille
	double result = 0., RecovCoef = 1.;
	double AVG_LR_CDO = 0.,AVG_LR_TRX = 0.;
	double PorS = 0.;

	if (cdo->GetBinaryFlg()) AVG_LR_CDO = 1 - cdo->GetBinary();
	else 
	{
		if (cdo->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = cdo->GetCollateral()->GetRecovCoef();
		for (i=0; i<bespoke_nbnames; i++) 
			AVG_LR_CDO += MIN(((ICM_ModelMultiCurves*)its_Callib_model)->GetRecoveryRate(cdo->GetCollateral()->GetIssuersLabels(i)) * RecovCoef,1);	
		AVG_LR_CDO /= (double) bespoke_nbnames; AVG_LR_CDO = 1 - AVG_LR_CDO;
	}

	//Les nominaux du portefeuille sont pris en compte dans l'estimation de Strike Equivalents
	if (!its_Rescaling_Full)
		{ cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}

	cdo->SetSubAmount(0.);

	double inf_attachment=0.;

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		double time,result = 0.;
		its_Callib_SizeEq_Index = 0.;
		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		//itsSlices[NoIndex].GetMaturities().push_back(time);
		
		int IndMaturity=0;

		//for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		//{
			result = 0.;
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
			its_Callib_PtfBespoke= cdo;
			its_Callib_K1_PtfBespoke= K1;
			its_Callib_K2_PtfBespoke= K2;
			its_Callib_NoIndex=NoIndex;
			//its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
			its_Callib_ActiveMaturity = time;
			PorS = cdo->GetPorS();

		
			int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
			for (int il_corr=0;il_corr<szcorrel;il_corr++)
			{inf_attachment = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100.;
			if (!CHECK_NULL(inf_attachment)) break;}

			//Maturité de l'indice de référence
			ARM_Date Maturity_2 = AsOf;
			Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
			its_Callib_Index = GenerateCdoIndex(Start,Maturity_2,NoIndex,PorS);

			//Estimation de la size equivalente (Average Loss de l'indice)
			AVG_LR_TRX = 0.;
			for (i=0; i<((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers(); i++) 
				AVG_LR_TRX += ((ICM_ModelMultiCurves*)its_Callib_model)->GetDefaultCurve(((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetIssuersLabels(i))->GetRecovery();	
			AVG_LR_TRX /= (double) ((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers(); AVG_LR_TRX = 1 - AVG_LR_TRX;

			its_Callib_SizeEq_Index=(AVG_LR_TRX/AVG_LR_CDO)*(K2-K1)*
									(((double)bespoke_nbnames)/((double)((ICM_Mez*)its_Callib_Index)->GetCollateral()->GetNbIssuers()));

			_sup2 = _sup - its_Callib_SizeEq_Index;

			try
			{
				//Tranche Equity
				if (flagequity)
				{
					//Test : existence de la solution
					test_inf = DiffSensisEquityTerm(_infeq);
					test_sup = DiffSensisEquityTerm(_sup);

					//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
					if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
						result = _infeq;
					//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
					else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
						result = _sup;
					//Cas Std : minimisation
					else
						//Pas de test car test en 0 -> size = 0 et prix Equity = 0
						result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEquityTerm,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
				}
				//Tranche Mezzanine
				else
				{
					//Test : tranche [x-100] : rescaling de la tranche equity corespondante
					if (flag100)
					{							
						//Test : existence de la solution
						test_inf = DiffSensisEqDownTerm(_infeq);
						test_sup = DiffSensisEqDownTerm(_sup);

						//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
						if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
							result = _inf;
						//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
						else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
							result = _sup;
						//Cas Std : minimisation
						else
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDownTerm,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
					}
					else
					{
						//Test existence de la solution
						test_inf = DiffSensisTerms(_inf);
						test_sup = DiffSensisTerms(_sup2);

						//Cas 1 : ELoss Bespoke > ELoss Traxx -> détachement par défaut à _inf
						if (test_inf < tolerance_en_0)
							result = _inf;
						//Cas 2 : ELoss Bespoke < ELoss Traxx -> On ne peut garder la meme size equivalente
						//typiquement pour les tranches mezzanines x - 100%
						else if (test_inf * test_sup > 0.)
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisTerms,(*this))).Dichotomy(_inf,_sup,100,1.E-4,1.E-4);
						//Cas std
						else
							result = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisTerms,(*this))).Dichotomy(_inf,_sup2,100,1.E-4,1.E-4);
					}
				}
			}
			catch (...)
			{
				if (cdo) delete cdo;
				if (its_Callib_model) delete its_Callib_model;
				if (its_Callib_Index) delete its_Callib_Index;
				ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
				itsSlices.clear();
			}

			if (its_Callib_Index) delete its_Callib_Index;
			//}

			if ((flagequity) && (result<inf_attachment)) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(inf_attachment);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(inf_attachment);
			}
			else if (flagequity) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result);
			}
			else if (flag100) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(1.);
			}
			else
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(MIN(result+its_Callib_SizeEq_Index,1.0));
			}

		}
	}

	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;

}


// --------------------------------------------------------------------
// Calculation of Diff sensitivities Term Structure
// --------------------------------------------------------------------
double ICM_Smile_Correlation::DiffSensisTerms(double Strike_down)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	double DefPv_Index_K1 = 0.,DefPv_bespoke_K1 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;

	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();;
	
	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;
	
	//On cap le point de détéchament max à 1
	double K2_index = MIN(Strike_down + its_Callib_SizeEq_Index,1);
	double K1_index = Strike_down;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	mod->SetCorrelation(this);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->SetAlready_rescal(true);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->KeepOneSlice(its_Callib_NoIndex);

	vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(K2_index);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);

	vector<double> SmileStrikeLow;SmileStrikeLow.push_back(K1_index);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	pricer_K2_traxx.SetTermStructurePricing(qNoTermStructure);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);
	
	//Bespoke Price (fonction du type de pricer collat fwd ou non)
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Equity Tranche Down
	traxx->SetMezzAmount(K1_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K1_PtfBespoke*SumNot_bespoke);

	mod->SetCorrelation(this);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->SetAlready_rescal(true);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->KeepOneSlice(its_Callib_NoIndex);

	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K1_traxx;pricer_K1_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_traxx.SetFaster(true);pricer_K1_traxx.SetStepForDistribLoss(-1);
	pricer_K1_traxx.SetTermStructurePricing(qNoTermStructure);
	DefPv_Index_K1 = pricer_K1_traxx.Price(qCMPDEFLEGPV);
	
	//Bespoke Price (fonction du type de pricer collat fwd ou non)
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			pricer_K1_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			pricer_K1_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K1_bespoke;pricer_K1_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K1_bespoke.SetFaster(true);pricer_K1_bespoke.SetStepForDistribLoss(-1);
			pricer_K1_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K1 = pricer_K1_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(Bespoke_YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result = (DefPv_Index_K2-DefPv_Index_K1)/(SumNot_traxx*(K2_index-K1_index));
    result -= (DefPv_bespoke_K2-DefPv_bespoke_K1)/(SumNot_bespoke*(its_Callib_K2_PtfBespoke-its_Callib_K1_PtfBespoke))/df_bespoke*df_index;

	return (result);
}


// --------------------------------------------------------------------
// Generate Fixed Base Correlation inside existing Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* FixedBaseCorrelationWithExistingBC(ARM_VolCurve* volBC,	
														  const ARM_Date& AsOf,
														  ARM_Currency* Ccy,
														  const double&	strikeDown,
														  const double& strikeUp,
														  const double& CorrelDown,
														  const double& CorrelUp,
														  const string& IndexName,
														  const double& yfmaturity)
{
	int i=0,k=0,l=1;
	int placeYT=volBC->GetExpiryTerms()->GetSize();
	int placeStrike=0;

	//la maturité est-elle déjà dans la liste des yearterms
	for (i=0;i<volBC->GetExpiryTerms()->GetSize();i++)
	{ if CHECK_EQUAL(yfmaturity,volBC->GetExpiryTerms()->Elt(i)) {l=0;break;}
	if (yfmaturity<volBC->GetExpiryTerms()->Elt(i)) {placeYT = i;}}

	ARM_Vector* YT = new ARM_Vector(volBC->GetExpiryTerms()->GetSize()+l,yfmaturity);
	
	for (i=0;i<volBC->GetExpiryTerms()->GetSize();i++)
	{ if (k == placeYT) {k++;}
	  YT->Elt(k)=volBC->GetExpiryTerms()->Elt(i);
	  k++; }

	for (i=0;i<volBC->GetStrikes()->GetSize();i++)
	{	if CHECK_EQUAL(strikeUp,volBC->GetStrikes()->Elt(i))
		{placeStrike=i;}}

	ARM_Vector* Strikes = (ARM_Vector*) volBC->GetStrikes()->Clone();
	ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());

	mVolatility->Elt(placeYT,placeStrike) = 100.*CorrelUp;
	if (placeStrike>0)
	{mVolatility->Elt(placeYT,placeStrike-1) = 100.*CorrelDown;}

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,Ccy);	


	vector<const ARM_VolCurve*> VVOL(1);
	VVOL[0]=VOL;

	std::vector<std::string> labels ;
	labels.push_back(IndexName); 
	// char** label = new char*[1];
	// label[0] = (char*)IndexName.c_str();

	ARM_Vector* propo = new ARM_Vector(1,1.);
	ARM_Vector* strikelow = new ARM_Vector(1,strikeDown);
	ARM_Vector* strikehigh = new ARM_Vector(1,strikeUp);
	string name = "FixedCorr";

	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
	// delete[] label;
	delete VOL; VOL = NULL;
	if (YT) delete YT;

	return (ForcedCorrelation);
}

// --------------------------------------------------------------------
// Generate Fixed Base Correlation inside existing Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* QuickBaseCorrelation(const ARM_Date& AsOf,
										    ARM_Currency* Ccy,
										    const double&	strikeDown,
										    const double& strikeUp,
										    const double& CorrelDown,
										    const double& CorrelUp,
										    const string& IndexName,
										    const double& yfmaturity,
											ARM_Vector* YearTerms,
											ARM_Vector* Strikes,
											ARM_Matrix* BC)
{
	int i=0,k=0,l=1;
	int placeYT=0;
	int placeStrike=0;

	//la maturité est-elle déjà dans la liste des yearterms
	for (i=0;i<YearTerms->GetSize();i++)
	{ if CHECK_EQUAL(yfmaturity,YearTerms->Elt(i)) {placeYT=i;l=0;break;}
	if (yfmaturity<YearTerms->Elt(i)) {placeYT = i;}
	else {placeYT = i+1;}
	}

	ARM_Vector* YT = new ARM_Vector(YearTerms->GetSize()+l,yfmaturity);
	
	for (i=0;i<YearTerms->GetSize();i++)
	{ if (k == placeYT) {k++;}
	  YT->Elt(k)=YearTerms->Elt(i);
	  k++; }

	for (i=0;i<Strikes->GetSize();i++)
	{	if CHECK_EQUAL(strikeUp,Strikes->Elt(i))
		{placeStrike=i;break;}}

	ARM_Vector* GStrikes = (ARM_Vector*) Strikes->Clone();
	ARM_Matrix* mVolatility = new ARM_Matrix(YearTerms->GetSize(),Strikes->GetSize());

	for (i=0;i<YT->GetSize();i++)
	for (k=0;k<Strikes->GetSize();k++)
		{if (i != placeYT)
		{mVolatility->Elt(i,k)=BC->Elt(i,k);}}
	
	mVolatility->Elt(placeYT,placeStrike) = 100.*CorrelUp;
	if (placeStrike>0)
	{mVolatility->Elt(placeYT,placeStrike-1) = 100.*CorrelDown;}

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												GStrikes, mVolatility, 1, K_SMILE_VOL,Ccy);	


	vector<const ARM_VolCurve*> VVOL(1);
	VVOL[0]=VOL;

	std::vector<std::string> labels ; 
	labels.push_back(IndexName); 

	ARM_Vector* propo = new ARM_Vector(1,1.);
	ARM_Vector* strikelow = new ARM_Vector(1,strikeDown);
	ARM_Vector* strikehigh = new ARM_Vector(1,strikeUp);
	string name = "FixedCorr";

	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
	// delete[] label;
	delete VOL; VOL = NULL;
	if (YT) delete YT;

	return (ForcedCorrelation);
}

// --------------------------------------------------------------------
// Generate Fixed Base Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* FixedBaseCorrelation2M(const ARM_Date& AsOf,
											ARM_Currency* Ccy,
											const double& strikeDown,
											const double& strikeUp,
											const double& CorrelDown,
											const double& CorrelUp,
											const double& Maturity,
											const double& PrevCorrelDown,
											const double& PrevCorrelUp,
											const double& PrevMaturity,
											const string& IndexName)
{
	ARM_Vector* YT = new ARM_Vector(2,0.);
	YT->Elt(0)=PrevMaturity;
	YT->Elt(1)=Maturity;

	ARM_Vector* Strikes = new ARM_Vector(2,0.);
	Strikes->Elt(0) = 100.*strikeDown;
	Strikes->Elt(1) = 100.*strikeUp;
	ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());
	mVolatility->Elt(1,0) = 100.*CorrelDown;
	mVolatility->Elt(1,1) = 100.*CorrelUp;

	mVolatility->Elt(0,0) = 100.*PrevCorrelDown;
	mVolatility->Elt(0,1) = 100.*PrevCorrelUp;

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,Ccy);	


	vector<const ARM_VolCurve*> VVOL(1);
	VVOL[0]=VOL;

	std::vector<std::string> labels ; 
	labels.push_back(IndexName); 

	// char** label = new char*[1];
	// label[0] = (char*)IndexName.c_str();

	ARM_Vector* propo = new ARM_Vector(1,1.);
	ARM_Vector* strikelow = new ARM_Vector(1,strikeDown);
	ARM_Vector* strikehigh = new ARM_Vector(1,strikeUp);
	string name = "FixedCorr";

	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
	// delete[] label;
	delete VOL; VOL = NULL;
	if (YT) delete YT;

	return (ForcedCorrelation);
}


void ICM_Smile_Correlation::GetCorrelationTerms(ARM_Vector& vect)
{
	vect = *itsSlices[0].GetVolCurve()->GetExpiryTerms();

	if (itsSlices.size()==1) {return;}

	ARM_Vector* mergedates = NULL;

	for (int j=1; j<itsSlices.size(); j++)
	{
		ARM_Vector* ET = itsSlices[j].GetVolCurve()->GetExpiryTerms();
		MergeDates(&mergedates,ET,&vect);
		vect = *mergedates;
		if (mergedates) delete mergedates;mergedates=NULL;
	}
}


// --------------------------------------------------------------------
// Generate Fixed Base Correlation
// --------------------------------------------------------------------
void UpdateSingleCorrelation(ICM_Smile_Correlation* correl,
							 const double& strike_down,
						     const double& strike_up,
							 const double& yt,
							 const double& value,
							 bool& succeed)
{
	succeed = false;

	vector<ARM_VolCurve*> VVol = correl->GetVolCurves();
	correl->GetSlices()[0].SetSmileStrikeLow(strike_down);
	correl->GetSlices()[0].SetSmileStrikeHigh(strike_up);

	ARM_VolLInterpol* volcurve = (ARM_VolLInterpol*) VVol[0];
	
	ARM_Vector* Strikes = volcurve->GetStrikes();
	ARM_Vector* ExpiryTerms = volcurve->GetExpiryTerms();
	ARM_Matrix* Volatilities = volcurve->GetVolatilities();
	int nLines = Volatilities->GetNumLines();
    int nCols  = Volatilities->GetNumCols();
	
	int nthCol=-1;
	int nthLine=-1;

	for (int i=0; i<nCols; i++)
	{
		if (fabs(strike_up*100.-Strikes->Elt(i))<DB_TOL)
		{
			nthCol = i;
			for (int j=0; j<nLines; j++)
			{
				if (fabs(yt-ExpiryTerms->Elt(j))<DB_TOL)
				{nthLine = j;break;}
			}
		}
	}

	if ((nthCol>=0) && (nthCol>=0))
	{
	Volatilities->Elt(nthLine,nthCol)= value*100.;
	succeed = true;
	}
}


// --------------------------------------------------------------------
// Generate Base Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* FixedBaseCorrelationEmpty(const ARM_Date& AsOf,
												ARM_Currency* Ccy,
												ARM_Vector& strikes,
												ARM_Vector& yearterms,
												ARM_Matrix& matrix)
{
	ARM_Vector* YT = (ARM_Vector*)yearterms.Clone();
	ARM_Vector* Strikes = (ARM_Vector*)strikes.Clone();
	(*Strikes) *= 100.;
	//ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());
	ARM_Matrix* mVolatility = (ARM_Matrix*) matrix.Clone();
	(*mVolatility)*=100.;

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,Ccy);	


	vector<const ARM_VolCurve*> VVOL(1);
	VVOL[0]=VOL;

	std::vector<std::string> labels ; 
	labels.push_back("CORR"); 

	// char** label = new char*[1];
	// label[0] = "CORR";

	ARM_Vector* propo = new ARM_Vector(1,1.);
	ARM_Vector* strikelow = new ARM_Vector(1,0.);
	ARM_Vector* strikehigh = new ARM_Vector(1,0.);
	string name = "FixedCorr";

	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
	// delete[] label;
	delete VOL; VOL = NULL;
	if (YT) delete YT;

	return (ForcedCorrelation);
}


// --------------------------------------------------------------------
// Generate Fixed Base Correlation
// --------------------------------------------------------------------
ICM_Smile_Correlation* FixedBaseCorrelation(const ARM_Date& AsOf,
											ARM_Currency* Ccy,
											const double& strikeDown,
											const double& strikeUp,
											const double& CorrelDown,
											const double& CorrelUp,
											const string& IndexName,
											double yfmaturity)
{
	ARM_Vector* YT = new ARM_Vector(1,yfmaturity);
	ARM_Vector* Strikes = new ARM_Vector(2,0.);
	Strikes->Elt(0) = 100.*strikeDown;
	Strikes->Elt(1) = 100.*strikeUp;
	ARM_Matrix* mVolatility = new ARM_Matrix(YT->GetSize(),Strikes->GetSize());
	mVolatility->Elt(0,0) = 100.*CorrelDown;
	mVolatility->Elt(0,1) = 100.*CorrelUp;

	ARM_VolCurve* VOL = (ARM_VolLInterpol*) new ARM_VolLInterpol( AsOf,(ARM_Vector*)YT->Clone(),
												Strikes, mVolatility, 1, K_SMILE_VOL,Ccy);	


	vector<const ARM_VolCurve*> VVOL(1);
	VVOL[0]=VOL;

	std::vector<std::string> labels ; 
	labels.push_back(IndexName); 

	// char** label = new char*[1];
	// label[0] = (char*)IndexName.c_str();

	ARM_Vector* propo = new ARM_Vector(1,1.);
	ARM_Vector* strikelow = new ARM_Vector(1,strikeDown);
	ARM_Vector* strikehigh = new ARM_Vector(1,strikeUp);
	string name = "FixedCorr";

	ICM_Smile_Correlation* ForcedCorrelation = new ICM_Smile_Correlation(AsOf,name,&(*VVOL.begin()),labels,propo,strikelow,strikehigh);
	//ForcedCorrelation->SetInterpType(K_ICM_STEPUP_RIGHT_MATURITY);	

	if (propo) delete propo;
	if (strikelow) delete strikelow;
	if (strikehigh) delete strikehigh;
	// delete[] label;
	delete VOL; VOL = NULL;
	if (YT) delete YT;

	return (ForcedCorrelation);
}


// --------------------------------------------------------------------
// Calculation of Diff sensitivities
// --------------------------------------------------------------------
double ICM_Smile_Correlation::DiffSensisEquityTerm(double Strike_up)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_up;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	mod->SetCorrelation(this);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->SetAlready_rescal(true);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->KeepOneSlice(its_Callib_NoIndex);

	//Update de la correl
	if (K2_index <= last_strike)
	{
		//Cas std : on prend la correl système
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(K2_index);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
	}
	else
	{
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(last_strike);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
	}

	vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	pricer_K2_traxx.SetTermStructurePricing(qNoTermStructure);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
		switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}

	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result = (DefPv_Index_K2)/(SumNot_traxx*Strike_up);
	result -= (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K2_PtfBespoke)/df_bespoke*df_index;

	return (result);
}

double ICM_Smile_Correlation::DiffSensisEqDownTerm(double Strike_down)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    
	//Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_down;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//MaJ : caractéristiques de la tranche down
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	double init_mezz= bespoke->GetMezzAmount(bespoke->GetFeeLeg()->GetStartDate());	

	bespoke->SetMezzAmount(its_Callib_K1_PtfBespoke*SumNot_bespoke);
	bespoke->SetSubAmount(0.);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	mod->SetCorrelation(this);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->SetAlready_rescal(true);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->KeepOneSlice(its_Callib_NoIndex);

	//Update de la correl
	if (K2_index <= last_strike)
	{
		//Cas std : on prend la correl système
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(K2_index);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);
	}
	else
	{
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(last_strike);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);
	}

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	pricer_K2_traxx.SetTermStructurePricing(qNoTermStructure);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	///Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex/Strike_down) : 1.);
	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke/its_Callib_K1_PtfBespoke):1.);

//	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex) : 1.);
//	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke):1.);

    //On ramène l'EL ZC à la maturité de la bespoke
    double result = (DefPv_Index_K2)/(SumNot_traxx*Strike_down*Callib_ELIndex);
	result -= (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K1_PtfBespoke*Callib_ELBespoke)/df_bespoke*df_index;

	//Remise a jour carac mezz
	bespoke->SetMezzAmount(init_mezz);
	bespoke->SetSubAmount(init_sub);

	return (result);
}


//-----------------------------------------------------------------------------
//Valuation Price using implied spread
//-----------------------------------------------------------------------------
 
double TSBaseCorrelCallib::ValuationCDOIndex(const double& x)
{
	ARM_Date AsOf;
	double NotionalCDO = 0.;
	string ccy("");
	if ( itsCorrel->GetVolCurves().size() != 0)
		ccy = string(itsCorrel->GetVolCurves().operator[](0)->GetCurrency()->GetCcyName());
	

	ICM_Mez* cdo = CdoIndexDefinition(its_Global_StartDate,
									its_Global_Maturity,
									its_Global_Collateral,
									itsK1,
									itsK2,
									x,
									NotionalCDO,
									its_Global_IncMatu,
									its_Global_AdjStart,
									its_Global_CreditLag, ccy ); 

	ICM_Pricer_Advisor Advisor;
	ICM_Smile_Correlation* Correlation = (ICM_Smile_Correlation*)itsMmc->GetCorrelation();

	ICM_Pricer*	pricer = Advisor.GeneratePricer(cdo,itsMmc,itsCname,CREDIT_DEFAULT_VALUE,&itsParams,Correlation->GetAsOfDate());

	Correlation->GetSlices()[0].SetSmileStrikeHigh(itsK2);
	Correlation->GetSlices()[0].SetSmileStrikeLow(itsK1);
	
	its_Global_FeePV=pricer->Price(qCMPFEELEGPV);
	its_Global_DefPV=pricer->Price(qCMPDEFLEGPV);
	its_Global_NPV=pricer->Price(qCMPPRICE);

	double value = its_Global_NPV;

	if (cdo) delete cdo;cdo=NULL;
	if (pricer) delete pricer;pricer=NULL;

	return (value);
}

ICM_QMatrix<double>* TSBaseCorrelCallib::GenerateTSBaseCorrelation(ARM_ZeroCurve* ircurve,
							   ICM_DefaultCurve* defcurveindex,
							   ICM_VolInterpol& basecorrel,
							   int creditlag)
{
	ARM_Date AsOf = ircurve->GetAsOfDate();
	ARM_Vector YFMaturities = *basecorrel.GetExpiryTerms();
	ARM_Vector Maturities = *basecorrel.GetExpiryTerms();
	Maturities *= 365.;
	Maturities += AsOf.GetJulian();

	ARM_Vector Strikes = *basecorrel.GetStrikes();
	Strikes /= 100.;

	ARM_Matrix Volatilities = *basecorrel.GetVolatilities();
	Volatilities /= 100.;
	ARM_Matrix VolatilitiesTS = Volatilities;

	int nLines = Volatilities.GetNumLines();
    int nCols  = Volatilities.GetNumCols();

	ARM_Date StartDate = AsOf +1.;
	double NotionalCDO = 0.;

	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = defcurveindex;

	ICM_Parameters Params;
	ICM_Parameters ParamsTS;

	ARM_CLASS_NAME cname;

	//definition du parametrage
	ARM_Vector PIntegrationStep(1,60);
	ARM_Vector PCopula(1,1.);
	ARM_Vector PIntegrationStep2(1,0.);
	ARM_Vector PFreedomDeg(1,0.);
	ARM_Vector TERMresc(1,1.);
	Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
	Params.Push(&PCopula,"COPULA");
	Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
	Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");

	ParamsTS =  Params ;
	ParamsTS.Push(&TERMresc,"TERMS_RESCALING");

	cname = ICM_PRICER_HOMOGENEOUS_SMILE;

	double _inf=0.,_sup=1.,slope=1.E-5,result=0.,guess=CREDIT_DEFAULT_VALUE;
	double step = 0.05;
	bool autoguess = false;

	ICM_Smile_Correlation*	BaseCorrelCal = FixedBaseCorrelationEmpty(AsOf,
												&ARM_Currency(defcurveindex->GetCurrency().c_str()),
												Strikes,
												YFMaturities,
												Volatilities);

	ICM_ModelMultiCurves mmc(DefaultCurves,ircurve,NULL,BaseCorrelCal);

	Volatilities *= 0.;
	ICM_Smile_Correlation*	BaseCorrelCalTS = FixedBaseCorrelationEmpty(AsOf,
												&ARM_Currency(defcurveindex->GetCurrency().c_str()),
												Strikes,
												YFMaturities,
												Volatilities);

	ICM_ModelMultiCurves mmcTS(DefaultCurves,ircurve,NULL,BaseCorrelCalTS);
	itsCorrel = BaseCorrelCalTS;

	ICM_Pricer* pricer = NULL;
	ICM_Pricer* pricerTS = NULL;

	itsMmc = &mmc;
	itsCname = cname; 
	itsParams =  Params;

	
	ICM_Mez* cdo = NULL;
	its_Global_StartDate = StartDate;

	//JLA char* Collateral[TRX_EUR_NBNAMES];
	//JLA for (int k=0; k<TRX_EUR_NBNAMES; k++)
	//JLA {Collateral[k]=(char*)defcurveindex->GetLabel().c_str();}
	//JLA memcpy(its_Global_Collateral,Collateral,sizeof(TypeCollateral));
	
	its_Global_Collateral.resize(TRX_EUR_NBNAMES) ;
	for (int k=0; k<TRX_EUR_NBNAMES; k++)
		its_Global_Collateral[k]=defcurveindex->GetLabel(); 

	its_Global_NBnames = TRX_EUR_NBNAMES;
	its_Global_IncMatu = INCLUDE_MATURITY;
	its_Global_AdjStart = false;
	its_Global_CreditLag = creditlag;

	ICM_QMatrix<double>* Matrix = new ICM_QMatrix<double>(Maturities.GetSize(),Strikes.GetSize(),0.);

	for (int j=0; j<Maturities.GetSize(); j++)
	for (int i=0; i<Strikes.GetSize();i++)
	{
		its_Global_FeePV=0.;
		its_Global_DefPV=0.;
		its_Global_NPV=0.;

		if (CHECK_NULL(Strikes.Elt(i)) && (i==0))
		{i++;}

		itsK1 = 0.;
		itsK2 = Strikes.Elt(i);

		if (i>0) {itsK1=Strikes.Elt(i-1);}

		if (lt(itsK2,BC_MIN_STRIKE))
			{continue;}
		else if (eq(itsK2,BC_MIN_STRIKE))
			{itsK1 = 0.;}
		
		itsPrice = 0.;

		its_Global_Maturity=(ARM_Date)Maturities[j];
		itsMaturity = YFMaturities.Elt(j);

		result = ValuationCDOIndex(0.01);
		result = its_Global_DefPV/its_Global_FeePV*0.01;

		bool succeed = false;
		UpdateSingleCorrelation(itsCorrel,itsK1,itsK2,itsMaturity,result,succeed);
		(*Matrix)(j,i)=result;

		if (succeed==false)
		{ICMLOG("TSBaseCorrelCallib::GenerateTSBaseCorrelation failed for "<<" i="<<i<<" j="<<j);}
		
		//try again with another guess
		if (result == CREDIT_DEFAULT_VALUE) 
		{	if (autoguess) 
			{guess += step; //increment guess with step
				//no solution available return CREDIT_DEFAULT_VALUE
				if (guess>1.) {autoguess=false;;guess=0.001;continue;}
			}
			else //start autoguess	
			{	guess = MAX(0.0001,guess-0.2);autoguess=true;}
			i--;continue;
		}
		else {autoguess = false;}

		guess=CREDIT_DEFAULT_VALUE;
	}

	if (BaseCorrelCal) delete BaseCorrelCal; BaseCorrelCal=NULL;
	if (itsCorrel) delete itsCorrel; itsCorrel=NULL;
// 	if (ParamsTS) delete ParamsTS;ParamsTS=NULL;

	return Matrix;

}


// --------------------------------------------------------------------
// Equivalent Strikes Computation Using ELoss stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_ELoss_Term(ICM_Pricer* pricer)
{
	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	int size = itsSlices.size();
	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double eloss_bespoke = 0.,eloss_indice = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date TrancheStart = ((ARM_SwapLeg*)cdo->GetFeeLeg())->GetStartDate();
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_nbnames = cdo->GetCollateral()->GetNbIssuers();
	double recovery_ratio = 0., poids = 0., nominal_total =0.;

	//Eloss du portefeuille : sommes des eloss ind
	its_Callib_model= (ICM_ModelMultiCurves*) pricer->GetModel()->Clone();
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	//Nominal total du portefeuille
	nominal_total = cdo->GetCollateral()->SumNotionals(cdo->GetEndDateNA());
	
	//calcul de l'EL du portefeuille
    eloss_bespoke = CptPtf_ELoss(pricer);

				
	double StrikeUpEloss = 0.;
	double StrikeDownEloss = 0.;
	if (!CHECK_NULL(eloss_bespoke))
	{
		StrikeUpEloss = (cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate()) + cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate())) / fabs(eloss_bespoke);
		StrikeDownEloss = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate()) / fabs(eloss_bespoke);
	}

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		double time,result_up = 0., result_down=0.;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;

		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		//Determination des bornes de projection en maturité
		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		if (itsSlices[NoIndex].GetProportion() != 0.)
		{
			//Maturité de l'indice de référence
			//its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
			its_Callib_ActiveMaturity = time;
			ARM_Date Maturity_2 = AsOf;
			Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
				ARM_Currency* Accy = itsSlices[NoIndex].GetIndex()->GetCurrencyUnit();
					
			//Estimation de l'eloss sur l'indice (exprimée en %)	
			eloss_indice = CptPtf_ELoss_Index(mod,itsSlices[NoIndex].GetIndex(),Maturity_2);
	
			//Strike Indice -> NB ELoss Indice
			int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
			vector<double> StrikeElossIndex(szcorrel);
			for (int il_corr=0;il_corr<szcorrel;il_corr++)
				StrikeElossIndex[il_corr] = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100./eloss_indice;

			ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes();
			vector<double> StrikeIndex(V->GetSize());
			for (int ij=0;ij<V->GetSize();ij++) {StrikeIndex[ij]=V->Elt(ij);}
				
			//Strike Equivalent	
			try
			{
				if (flagequity)
					result_up = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeUpEloss)/100;
				else
				{
					//Tranche Up
					result_up = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeUpEloss)/100;
						
					//Tranche Down
					result_down = LinearVectorInterpol(StrikeElossIndex, StrikeIndex, StrikeDownEloss)/100;
				}
			}
			catch (...)
			{
				if (cdo) delete cdo;
				if (its_Callib_model) delete its_Callib_model;
				if (its_Callib_Index) delete its_Callib_Index;
				ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
				itsSlices.clear();
			}
			if (its_Callib_Index) delete its_Callib_Index;

		}
		if (flagequity) 
		{
			itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
			itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
		}
		else
		{
			itsSlices[NoIndex].GetSmileStrikeLow().push_back(result_down);
			itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
		}
	}
	//Traces
	//fclose(stream);
	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}


// --------------------------------------------------------------------
// Equivalent Strikes Computation Using Equity stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_Equity_Term(ICM_Pricer* pricer)
{
	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	if (its_NormalizeByEL)
	{its_Callib_ELBespoke=CptPtf_ELoss(pricer,K_ZEROCOUPON);}

	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = K1 + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());
	int size = itsSlices.size();
	int i = 0;
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double FixedSensiPtf = 0.,SensiPtf = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_nbnames = cdo->GetCollateral()->GetNbIssuers();

	its_Callib_model= (ARM_Model*) pricer->GetModel()->Clone();

	//Borne inf de la dichotomie
 	double _inf = 0., _infeq=1.e-4;
	
	//Borne Sup de la dichotomie : varie en fonction de la maturité (a modifier si existence plot <5Y)
	double _sup = 0.,_sup2 = 0.;
	_sup = MIN(K2+0.5,0.99);

	//Test resultat numerique
	double test_inf=0., test_sup = 0.;

	//Estimation de l'Average Loss du portefeuille
	double result = 0.;
	double AVG_LR_CDO = 0.,AVG_LR_TRX = 0.;
	double PorS = 0.;

	if (!its_Rescaling_Full)
		{cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}
	
	double inf_attachment=0.;

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		// double val_inf, val_sup,time,result = 0.;
		double time,result_up = 0., result_down=0.;
		//int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;
		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		result = 0.;
		if (itsSlices[NoIndex].GetProportion() != 0.)
		{
			its_Callib_PtfBespoke= cdo;
			its_Callib_K1_PtfBespoke= K1;
			its_Callib_K2_PtfBespoke= K2;
			its_Callib_NoIndex=NoIndex;
			//its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
			its_Callib_ActiveMaturity = time;
			PorS = cdo->GetPorS();

			int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
			for (int il_corr=0;il_corr<szcorrel;il_corr++)
			{inf_attachment = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100.;
			if (!CHECK_NULL(inf_attachment)) break;}

			//Maturité de l'indice de référence
			ARM_Date Maturity_2 = AsOf;
			Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
			its_Callib_Index = GenerateCdoIndex(Start,Maturity_2,NoIndex,PorS);

			//calcul de l'El d'indice pour la normalisation
			if (its_NormalizeByEL)
			{its_Callib_ELIndex = CptPtf_ELoss_Index((ICM_ModelMultiCurves *)its_Callib_model
													,itsSlices[NoIndex].GetIndex(),
													Maturity_2,K_ZEROCOUPON);}
			else
			{its_Callib_ELIndex = 1.;}

			try
			{
				if (flagequity)
				{
					//Test : existence de la solution
					test_inf = DiffSensisEqUpTerm(_infeq);
					test_sup = DiffSensisEqUpTerm(_sup);

					if (fabs(test_sup)<=1.E-4)
							result_up = _sup;
					//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
					else if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
						result_up = _inf;
					//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
					else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
						result_up = _sup;
					//Cas Std : minimisation
					else
						result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUpTerm,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
				}
				else
				{
					//Tranche Up
						//Test : existence de la solution
						test_inf = DiffSensisEqUpTerm(_infeq);
						test_sup = DiffSensisEqUpTerm(_sup);

						//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
						if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
							result_up = _inf;
						//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
						else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
							result_up = _sup;
						//Cas Std : minimisation
						else
							result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUpTerm,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
					
						//Tranche Down
						_sup = MIN(K1+0.5,0.99);
						//Test : existence de la solution
						test_inf = DiffSensisEqDownTerm(_infeq);
						test_sup = DiffSensisEqDownTerm(_sup);

						//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
						if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
							result_down = _inf;
						//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
						else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
							result_down = _sup;
						//Cas Std : minimisation
						else
							result_down = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDownTerm,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
				}
			}
			catch (...)
			{
				if (cdo) delete cdo;
				if (its_Callib_model) delete its_Callib_model;
				if (its_Callib_Index) delete its_Callib_Index;
				ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
				itsSlices.clear();
			}
			if (its_Callib_Index) delete its_Callib_Index;
				
			}

			if (flagequity) 
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
			else
			{
				itsSlices[NoIndex].GetSmileStrikeLow().push_back(result_down);
				itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
			}
	}

	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}
	
double ICM_Smile_Correlation::DiffSensisEqUpTerm(double Strike_up)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_up;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	traxx->SetMezzAmount(K2_index*SumNot_traxx);

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	if (!CHECK_NULL(bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate())))
		bespoke->SetSubAmount(0.);
		
	bespoke->SetMezzAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	mod->SetCorrelation(this);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->SetAlready_rescal(true);
	(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->KeepOneSlice(its_Callib_NoIndex);

	//Update de la correl
	if (K2_index <= last_strike)
	{
		//Cas std : on prend la correl système
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(K2_index);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);
	}
	else
	{
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(last_strike);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeHigh(SmileStrikeHigh);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		(dynamic_cast<ICM_Smile_Correlation*>(mod->GetCorrelation()))->itsSlices[0].SetSmileStrikeLow(SmileStrikeLow);
	}

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	pricer_K2_traxx.SetTermStructurePricing(qNoTermStructure);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			pricer_K2_bespoke.SetTermStructurePricing(qNoTermStructure);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex/Strike_up) : 1.);
	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke/its_Callib_K2_PtfBespoke):1.);

//	double Callib_ELIndex =((its_Callib_ELIndex!=1.) ? (its_Callib_ELIndex) : 1.);
//	double Callib_ELBespoke =((its_Callib_ELBespoke!=1.)?(its_Callib_ELBespoke):1.);

    //On ramène l'EL ZC à la maturité de la bespoke
    double result1 = (DefPv_Index_K2)/(SumNot_traxx*Strike_up*Callib_ELIndex);
	if ((its_Callib_ELIndex!=1.) && (result1>0.99) && (Strike_up>=0.99))
		{result1 = 1.;}

	double result2 = (DefPv_bespoke_K2)/(SumNot_bespoke*its_Callib_K2_PtfBespoke*Callib_ELBespoke);
	if ((its_Callib_ELBespoke!=1.) && (result2>0.99) && (Strike_up>=0.99))
		{result2 = 1.;}

	double result = result1;
	if (its_Callib_ELBespoke!=1.)
	{result -= result2;}
	else
	{result -= result2/df_bespoke*df_index;}

	//Remise a jour sub amount init
	bespoke->SetSubAmount(init_sub);

	return (result);
}


// --------------------------------------------------------------------
// Equivalent Strikes Computation Using Equity stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikesEq_Equity_CombinIDX(ICM_Pricer* pricer)
{
	using namespace OptimTools;

	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	// pour la normalisation par l'EL du ptf bespoke
	if (its_NormalizeByEL)
	{its_Callib_ELBespoke=CptPtf_ELoss(pricer,K_ZEROCOUPON);}

	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = K1 + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());
	int size = itsSlices.size();
	int NoIndex=0;
	int i = 0;
	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double FixedSensiPtf = 0.,SensiPtf = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();
	int bespoke_nbnames = cdo->GetCollateral()->GetNbIssuers();

	its_Callib_model= (ARM_Model*) pricer->GetModel()->Clone();

	//Borne inf de la dichotomie
 	double _inf = 0., _infeq=1.e-4;
	
	//Borne Sup de la dichotomie : varie en fonction de la maturité (a modifier si existence plot <5Y)
	double _sup = 0.,_sup2 = 0.;
	_sup = MIN(K2+0.5,0.99);

	//Test resultat numerique
	double test_inf=0., test_sup = 0.;

	//Estimation de l'Average Loss du portefeuille
	double result = 0.;
	double AVG_LR_CDO = 0.,AVG_LR_TRX = 0.;
	double PorS = 0.;

	if (!its_Rescaling_Full)
	{ cdo->CptCashFlowDatesCredit(K_ZEROCOUPON);}
	
	double inf_attachment=0.;

	_context ctxt;
	ctxt.m_Callib_PtfBespoke = cdo;
	ctxt.m_Callib_model = its_Callib_model;
	ctxt.m_Active_correl = this;
	ctxt.m_Slices = &GetSlices();
	ctxt.m_Callib_ELBespoke = its_Callib_ELBespoke;
	ctxt.m_Callib_K1_PtfBespoke = K1;
	ctxt.m_Callib_K2_PtfBespoke = K2;
	ctxt.m_PricerType = pricer->GetName();

	for (NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		// double val_inf, val_sup,time,result = 0.;
		double time;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;
		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		//Borne Inf
		if (idx_inf!=CREDIT_DEFAULT_VALUE) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_inf]);

		//Borne Sup
		if ((isequal==false) && (idx_sup!=CREDIT_DEFAULT_VALUE)) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_sup]);

		//only maturity
		if (its_Interpolation_At_Maturity)
		{
			itsSlices[NoIndex].GetMaturities().clear();
			itsSlices[NoIndex].GetMaturities().push_back(time);
		}

		_context_index ci;

		
		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{
			result = 0.;
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
				its_Callib_PtfBespoke= cdo;
				its_Callib_K1_PtfBespoke= K1;
				its_Callib_K2_PtfBespoke= K2;
				its_Callib_NoIndex=NoIndex;
				its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];
				PorS = cdo->GetPorS();

				int szcorrel = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->GetSize();
				for (int il_corr=0;il_corr<szcorrel;il_corr++)
				{inf_attachment = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetStrikes()->Elt(il_corr)/100.;
				if (!CHECK_NULL(inf_attachment)) break;}

				//Maturité de l'indice de référence
				ARM_Date Maturity_2 = AsOf;
				Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));
				its_Callib_Index = GenerateCdoIndex(Start,Maturity_2,NoIndex,PorS);

				//calcul de l'El d'indice pour la normalisation
				if (its_NormalizeByEL)
					{its_Callib_ELIndex = CptPtf_ELoss_Index((ICM_ModelMultiCurves *)its_Callib_model,
													itsSlices[NoIndex].GetIndex(),Maturity_2,K_ZEROCOUPON);}
				else
					{its_Callib_ELIndex = 1.;}

				ci.m_Callib_Index.push_back((ICM_Mez*)its_Callib_Index);
				ci.m_Callib_ELIndex.push_back(its_Callib_ELIndex);
				ci.m_Callib_ActiveMaturity.push_back(its_Callib_ActiveMaturity);
			}
		}
		
		ctxt.m_Callib_Index_Vector.push_back(ci);
	}


	double result_up = 0., result_down=0.;
	vector<double> bound_inf_X;
	vector<double> bound_sup_X;

//	bound_inf_X.clear();bound_sup_X.clear();
//	for (i=0;i<size;i++){bound_inf_X.push_back(_infeq);bound_sup_X.push_back(_sup);X[i] = (_infeq+_sup)/2.;}


	try
	{
		if (flagequity)
		{
/*			//Test : existence de la solution
			test_inf = DiffSensisEqUp(_infeq);
			test_sup = DiffSensisEqUp(_sup);

			//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
			if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
				result_up = _inf;
			//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
			else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
				result_up = _sup;
			//Cas Std : minimisation
			else
				result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
*/
			vector<double> X;X.resize(size);
			bound_inf_X.clear();bound_sup_X.clear();
			for (i=0;i<size;i++){bound_inf_X.push_back(_infeq);bound_sup_X.push_back(_sup);X[i] = (_infeq+_sup)/2.;}

			result_up = OptimTools::NagMultiOptimisator((void*) &ctxt,objfun_DiffSensisEqUp,X,bound_inf_X,bound_sup_X,1.E-3,100);

			itsSlices[NoIndex].GetSmileStrikeLow().push_back(0.);
			itsSlices[NoIndex].GetSmileStrikeHigh().push_back(result_up);
		}
		else
		{
			//Tranche Up
			//Test : existence de la solution
/*				test_inf = DiffSensisEqUp(_infeq);
				test_sup = DiffSensisEqUp(_sup);

				//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
				if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
					result_up = _inf;
				//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
				else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
					result_up = _sup;
				//Cas Std : minimisation
				else
					result_up = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqUp,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
*/
			vector<double> X1;X1.resize(size);
			bound_inf_X.clear();bound_sup_X.clear();
			for (i=0;i<size;i++){bound_inf_X.push_back(_infeq);bound_sup_X.push_back(_sup);X1[i] = K2;}

			result_up = OptimTools::NagMultiOptimisator((void*) &ctxt,objfun_DiffSensisEqUp,X1,bound_inf_X,bound_sup_X,1.E-3,100);

			for (NoIndex=0;NoIndex<size;NoIndex++)
			{itsSlices[NoIndex].GetSmileStrikeHigh().push_back(X1[NoIndex]);}
			
			//Tranche Down
			_sup = MIN(K1+0.5,0.99);


/*			//Test : existence de la solution
			test_inf = DiffSensisEqDown(_infeq);
			test_sup = DiffSensisEqDown(_sup);

			//Cas 1 : Eloss Bespoke > Eloss Trax  -> detachement par defaut à _infeq
			if ( (test_inf * test_sup > 0.) && (test_inf <0.) )
				result_down = _inf;
			//Cas 2 : Eloss Bespoke < Eloss Trax  -> detachement par defaut à _sup
			else if ( (test_inf * test_sup > 0.) && (test_sup > 0.) )
				result_down = _sup;
			//Cas Std : minimisation
			else
				result_down = RootFinder1D(ff1::mem_call(&ICM_Smile_Correlation::DiffSensisEqDown,(*this))).Dichotomy(_infeq,_sup,100,1.E-4,1.E-4);
*/
			vector<double> X2;X2.resize(size);
			bound_inf_X.clear();bound_sup_X.clear();
			for (i=0;i<size;i++){bound_inf_X.push_back(_infeq);bound_sup_X.push_back(_sup);X2[i] = K1;}

			result_up = OptimTools::NagMultiOptimisator((void*) &ctxt,objfun_DiffSensisEqDown,X2,bound_inf_X,bound_sup_X,1.E-3,100);

			for (NoIndex=0;NoIndex<size;NoIndex++)
			{itsSlices[NoIndex].GetSmileStrikeLow().push_back(X2[NoIndex]);}
		}
	}
	catch (...)
	{
		if (cdo) delete cdo;
		if (its_Callib_model) delete its_Callib_model;
		ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
		itsSlices.clear();
	}
	
	for (NoIndex=0;NoIndex<size;NoIndex++)
		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{
			if (ctxt.m_Callib_Index_Vector[NoIndex].m_Callib_Index[IndMaturity]) 
				delete ctxt.m_Callib_Index_Vector[NoIndex].m_Callib_Index[IndMaturity];
			ctxt.m_Callib_Index_Vector[NoIndex].m_Callib_Index[IndMaturity]=NULL;
		}
			
	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}


// ---------------------------------------------------------------------------------
// Optimisation simultanée des parameters for DiffSensisEqUp
// ---------------------------------------------------------------------------------
static void __stdcall objfun_DiffSensisEqUp(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	int i=0;
	vector<double> V_Strike_up;
	for (i=0;i<n;i++)
	{
		if (x[i]<=1.e-20)
		{	*objf = 10.E20;
			return;	}
		V_Strike_up.push_back(x[i]);
	}

    _context* p = (_context*)comm->p ;

	int nbindex = p->m_Slices->size();

	double result=0.;

	ICM_Smile_Correlation* ActivCorrel = p->m_Active_correl;
	ICM_Smile_Correlation* Fullcorrel = (ICM_Smile_Correlation*) ActivCorrel->Clone();

	for (i=0;i<nbindex;i++)
	{
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(V_Strike_up[i]);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		Fullcorrel->GetSlices()[i].SetSmileStrikeHigh(SmileStrikeHigh);
		Fullcorrel->GetSlices()[i].SetSmileStrikeLow(SmileStrikeLow);
	}

	Fullcorrel->SetAlready_rescal(true);

	// -------------------------------------------------------------------
	// calcul de l'EL pour la bespoke 
	// -------------------------------------------------------------------
	double DefPv_bespoke_K2 = 0.;
	ICM_Mez* bespoke = (ICM_Mez*) p->m_Callib_PtfBespoke;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) p->m_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();

	//Bespoke Mat
	ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
	double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;
	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	if (!CHECK_NULL(bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate())))
			bespoke->SetSubAmount(0.);

	bespoke->SetMezzAmount(p->m_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
	mod->SetCorrelation(Fullcorrel);


	switch (p->m_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
	double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
	double ELBespoke = (DefPv_bespoke_K2)/(SumNot_bespoke*p->m_Callib_K2_PtfBespoke*p->m_Callib_ELBespoke)/df_bespoke;

	//Remise a jour sub amount init
	bespoke->SetSubAmount(init_sub);

	if (initzc) delete initzc;
	initzc = NULL;

	if (Fullcorrel) delete Fullcorrel;
	Fullcorrel=NULL;

	// -------------------------------------------------------------------
	// calcul de l'EL pour les indices
	// -------------------------------------------------------------------

	for (i=0;i<nbindex;i++)
	{
		int matuNO = p->m_Callib_MaturityNo;
		_context_index cti = p->m_Callib_Index_Vector[i];
		double actMatu = cti.m_Callib_ActiveMaturity[matuNO];

		double DefPv_Index_K2 = 0.;

		ICM_Mez* traxx = (ICM_Mez*) cti.m_Callib_Index[matuNO];
		ARM_ZeroCurve* indxzc = mod->GetDefaultCurve((*(p->m_Slices))[i].GetIndex()->GetLabels()[0])->GetZeroCurve();

		//Index Mat
		ARM_Date Maturity = traxx->GetEndDateNA();
		double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
		
		double K2_index = V_Strike_up[i];
		double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

		traxx->SetMezzAmount(K2_index*SumNot_traxx);

		//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
		//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
		//solutions et donc des pb d'estimation numérique.

		double last_strike = 1. ,Beta_index_i_k2 = 0.;

		//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
		if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
			last_strike = 0.22;
		else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
			last_strike = 0.3;

		//Update de la correl
		if (K2_index <= last_strike)
			//Cas std : on prend la correl système
			Beta_index_i_k2 = (*(p->m_Slices))[i].GetVolCurve()->ComputeVolatility(actMatu,100.*K2_index);
		else
			//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
			Beta_index_i_k2 = (*(p->m_Slices))[i].GetVolCurve()->ComputeVolatility(actMatu,100.*last_strike);

		Beta_index_i_k2 /= 100.;
	
		ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
		mod->SetCorrelation(&Correl_k2);

		//Index Price
		mod->SetZeroCurve(indxzc);
		ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
		DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

		double df_index = indxzc->DiscountPrice(YF_Maturity);
	
	    //On ramène l'EL ZC à la maturité de la bespoke
		double intermediate = (DefPv_Index_K2)/(SumNot_traxx*V_Strike_up[i]*cti.m_Callib_ELIndex[matuNO]);
		intermediate -= ELBespoke*df_index;
	
		result+= fabs(intermediate);
		}

	*objf = result;
}

// ---------------------------------------------------------------------------------
// Optimisation simultanée des parameters for DiffSensisEqDown
// ---------------------------------------------------------------------------------
static void __stdcall objfun_DiffSensisEqDown(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	int i=0;
	vector<double> V_Strike_down;
	for (i=0;i<n;i++)
	{
		if (x[i]<=1.e-20)
		{	*objf = 10.E20;
			return;	}
		V_Strike_down.push_back(x[i]);
	}

    _context* p = (_context*)comm->p ;

	int nbindex = p->m_Slices->size();

	double result=0.;

	ICM_Smile_Correlation* ActivCorrel = p->m_Active_correl;
	ICM_Smile_Correlation* Fullcorrel = (ICM_Smile_Correlation*) ActivCorrel->Clone();

	for (i=0;i<nbindex;i++)
	{
		vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(V_Strike_down[i]);
		vector<double> SmileStrikeLow;SmileStrikeLow.push_back(0.);
		Fullcorrel->GetSlices()[i].SetSmileStrikeHigh(SmileStrikeHigh);
		Fullcorrel->GetSlices()[i].SetSmileStrikeLow(SmileStrikeLow);
	}

	Fullcorrel->SetAlready_rescal(true);

	// -------------------------------------------------------------------
	// calcul de l'EL pour la bespoke 
	// -------------------------------------------------------------------
	double DefPv_bespoke_K2 = 0.;
	ICM_Mez* bespoke = (ICM_Mez*) p->m_Callib_PtfBespoke;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) p->m_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();

	//Bespoke Mat
	ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
	double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;
	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 

	bespoke->SetMezzAmount(p->m_Callib_K2_PtfBespoke*SumNot_bespoke);

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	double init_mezz = bespoke->GetMezzAmount(bespoke->GetFeeLeg()->GetStartDate());
	
	bespoke->SetSubAmount(0.);
	bespoke->SetMezzAmount(p->m_Callib_K1_PtfBespoke*SumNot_bespoke);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
	mod->SetCorrelation(Fullcorrel);


	switch (p->m_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke;pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
	double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
	double ELBespoke = (DefPv_bespoke_K2)/(SumNot_bespoke*p->m_Callib_K1_PtfBespoke*p->m_Callib_ELBespoke)/df_bespoke;

	//Remise a jour sub amount init
	bespoke->SetMezzAmount(init_mezz);
	bespoke->SetSubAmount(init_sub);

	if (initzc) delete initzc;
	initzc = NULL;

	if (Fullcorrel) delete Fullcorrel;
	Fullcorrel=NULL;

	// -------------------------------------------------------------------
	// calcul de l'EL pour les indices
	// -------------------------------------------------------------------

	for (i=0;i<nbindex;i++)
	{
		int matuNO = p->m_Callib_MaturityNo;
		_context_index cti = p->m_Callib_Index_Vector[i];
		double actMatu = cti.m_Callib_ActiveMaturity[matuNO];
		double DefPv_Index_K2 = 0.;

		ICM_Mez* traxx = (ICM_Mez*) cti.m_Callib_Index[matuNO];
		ARM_ZeroCurve* indxzc = mod->GetDefaultCurve((*(p->m_Slices))[i].GetIndex()->GetLabels()[0])->GetZeroCurve();

		//Index Mat
		ARM_Date Maturity = traxx->GetEndDateNA();
		double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
		
		double K2_index = V_Strike_down[i];
		double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

		traxx->SetMezzAmount(K2_index*SumNot_traxx);

		//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
		//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
		//solutions et donc des pb d'estimation numérique.

		double last_strike = 1. ,Beta_index_i_k2 = 0.;

		//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
		if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
			last_strike = 0.22;
		else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
			last_strike = 0.3;

		//Update de la correl
		if (K2_index <= last_strike)
			//Cas std : on prend la correl système
			Beta_index_i_k2 = (*(p->m_Slices))[i].GetVolCurve()->ComputeVolatility(actMatu,100.*K2_index);
		else
			//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
			Beta_index_i_k2 = (*(p->m_Slices))[i].GetVolCurve()->ComputeVolatility(actMatu,100.*last_strike);

		Beta_index_i_k2 /= 100.;
	
		ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
		mod->SetCorrelation(&Correl_k2);

		//Index Price
		mod->SetZeroCurve(indxzc);
		ICM_Pricer_Distrib_Smile pricer_K2_traxx;pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),ActivCorrel->GetAsOfDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
		DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

		double df_index = indxzc->DiscountPrice(YF_Maturity);
	
	    //On ramène l'EL ZC à la maturité de la bespoke
		double intermediate = (DefPv_Index_K2)/(SumNot_traxx*V_Strike_down[i]*cti.m_Callib_ELIndex[matuNO]);
		intermediate -= ELBespoke*df_index;
	
		result+= fabs(intermediate);
		}

	*objf = result;
}

// ----------------------------------------------------------------------------------
// Digital Case
// ----------------------------------------------------------------------------------
double ICM_Smile_Correlation::DiffSensisEqUp_digit(double Strike_up)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    //Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_up;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub_traxx = traxx->GetSubAmount(traxx->GetFeeLeg()->GetStartDate());
	traxx->SetSubAmount(K2_index*SumNot_traxx);
	traxx->SetMezzAmount((EPSILON_RESCALING)*SumNot_traxx);

	//Modif : point d'attachement à 0 sur les tranches mezz
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	bespoke->SetSubAmount(its_Callib_K2_PtfBespoke*SumNot_bespoke);
	bespoke->SetMezzAmount((EPSILON_RESCALING)*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.,Beta_index_i_k2e = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	//Update de la correl
	if (K2_index <= last_strike)
	{
		//Cas std : on prend la correl système
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
		Beta_index_i_k2e = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*(K2_index+EPSILON_RESCALING));
	}
	else
	{
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*last_strike);
		Beta_index_i_k2e = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*(last_strike+EPSILON_RESCALING));
	}

	Beta_index_i_k2 /= 100.;
	Beta_index_i_k2e /= 100.;
	
	ICM_Smile_Correlation* Correl_k2 = FixedBaseCorrelation(mod->GetStartDate(),
															  indxzc->GetCurrencyUnit(),
															  K2_index,
															  K2_index+EPSILON_RESCALING,
															  Beta_index_i_k2,
															  Beta_index_i_k2e,
															  "NAME_CORR");

	mod->SetCorrelation(Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx; pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	//Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result1 = (DefPv_Index_K2)/(SumNot_traxx*EPSILON_RESCALING);
	double result2 = (DefPv_bespoke_K2)/(SumNot_bespoke*EPSILON_RESCALING);

	double result = result1 - result2/df_bespoke*df_index;;

	//Remise a jour sub amount init
	bespoke->SetSubAmount(init_sub);
	traxx->SetSubAmount(init_sub_traxx);

	//Free mem
	if (Correl_k2) delete Correl_k2;
	return (result);
}


double ICM_Smile_Correlation::DiffSensisEqDown_digit(double Strike_down)
{
	double DefPv_Index_K2 = 0.,DefPv_bespoke_K2 = 0.;
	int size_prop = itsSlices.size();

	ICM_Mez* bespoke = (ICM_Mez*) its_Callib_PtfBespoke;
	ICM_Mez* traxx = (ICM_Mez*) its_Callib_Index;

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ICM_Correlation* Correl = NULL;
	ARM_ZeroCurve* initzc = (ARM_ZeroCurve*) mod->GetZeroCurve()->Clone();
	ARM_ZeroCurve* indxzc = mod->GetDefaultCurve(itsSlices[its_Callib_NoIndex].GetIndex()->GetLabels()[0])->GetZeroCurve();

	//Index Mat
    ARM_Date Maturity = traxx->GetEndDateNA();
    double YF_Maturity = (Maturity-mod->GetStartDate())/365.;
    
	//Bespoke Mat
    ARM_Date Bespoke_Maturity = bespoke->GetEndDateNA();
    double Bespoke_YF_Maturity = (Bespoke_Maturity-mod->GetStartDate())/365.;

	double K2_index = Strike_down;

	double SumNot_bespoke = bespoke->GetCollateral()->SumNotionals(bespoke->GetEndDateNA()); 
	double SumNot_traxx = traxx->GetCollateral()->SumNotionals(traxx->GetEndDateNA()); 

	//MaJ : caractéristiques de la tranche down
	double init_sub_traxx = traxx->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	double init_mezz_traxx= traxx->GetMezzAmount(bespoke->GetFeeLeg()->GetStartDate());	

	traxx->SetMezzAmount((EPSILON_RESCALING)*SumNot_traxx);
	traxx->SetSubAmount(K2_index*SumNot_traxx);

	//MaJ : caractéristiques de la tranche down
	double init_sub = bespoke->GetSubAmount(bespoke->GetFeeLeg()->GetStartDate());
	double init_mezz= bespoke->GetMezzAmount(bespoke->GetFeeLeg()->GetStartDate());	

	bespoke->SetMezzAmount((EPSILON_RESCALING)*SumNot_bespoke);
	bespoke->SetSubAmount((its_Callib_K1_PtfBespoke)*SumNot_bespoke);

	//Dans le cas des tranches equity : on ne prend pas en compte le niveau de correlation associé au 
	//plots supérieurs aux plots de marché (typiquement le strike 50%) pour eviter d'avoir plusieures
	//solutions et donc des pb d'estimation numérique.

	double last_strike = 1. ,Beta_index_i_k2 = 0.,Beta_index_i_k2e = 0.;

	//Estimation du dernier strike utilisable (en fonction de la ccy de l'index)
	if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"EUR") == 0)
		last_strike = 0.22;
	else if (strcmp(indxzc->GetCurrencyUnit()->GetCcyName(),"USD") == 0)
		last_strike = 0.3;

	//Update de la correl
	if (K2_index <= last_strike)
	{
		//Cas std : on prend la correl système
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*K2_index);
		Beta_index_i_k2e = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*(K2_index+EPSILON_RESCALING));
	}
	else
	{
		//Cas Strike lointain : on fige la correl a partir du dernier strike de marche (EUR :22%, USD : 30%)
		Beta_index_i_k2 = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*last_strike);
		Beta_index_i_k2e = itsSlices[its_Callib_NoIndex].GetVolCurve()->ComputeVolatility(its_Callib_ActiveMaturity,100.*(last_strike+EPSILON_RESCALING));
	}

	Beta_index_i_k2 /= 100.;
	Beta_index_i_k2e /= 100.;
	
	ICM_Smile_Correlation* Correl_k2 = FixedBaseCorrelation(mod->GetStartDate(),
															  indxzc->GetCurrencyUnit(),
															  K2_index,
															  K2_index+EPSILON_RESCALING,
															  Beta_index_i_k2,
															  Beta_index_i_k2e,
															  "NAME_CORR");
	
	//ICM_Beta_Correlation Correl_k2(mod->GetStartDate(),NEG_SQRT(Beta_index_i_k2),"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	mod->SetCorrelation(Correl_k2);

	//Index Price
	mod->SetZeroCurve(indxzc);
	ICM_Pricer_Distrib_Smile pricer_K2_traxx; pricer_K2_traxx.Set(traxx,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_traxx.SetFaster(true);pricer_K2_traxx.SetStepForDistribLoss(-1);
	DefPv_Index_K2 = pricer_K2_traxx.Price(qCMPDEFLEGPV);

	///Bespoke Price
	mod->SetZeroCurve(initzc);
	switch (its_PricerType)
	{
	case ICM_PRICER_HOMOGENEOUS_SMILE:
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
		{
			ICM_Pricer_Distrib_Smile_Collat_Fwd pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	default :
		{
			ICM_Pricer_Distrib_Smile pricer_K2_bespoke; pricer_K2_bespoke.Set(bespoke,mod,ICM_Parameters(),mod->GetStartDate());pricer_K2_bespoke.SetFaster(true);pricer_K2_bespoke.SetStepForDistribLoss(-1);
			DefPv_bespoke_K2 = pricer_K2_bespoke.Price(qCMPDEFLEGPV);
			break;
		}
	}
	
	//Modif CN : prise en compte des diff de discount factors
    double df_bespoke = initzc->DiscountPrice(Bespoke_YF_Maturity);
    double df_index = indxzc->DiscountPrice(YF_Maturity);

	if (initzc) delete initzc;

    //On ramène l'EL ZC à la maturité de la bespoke
    double result1 = (DefPv_Index_K2)/(SumNot_traxx*EPSILON_RESCALING);
	double result2 = (DefPv_bespoke_K2)/(SumNot_bespoke*EPSILON_RESCALING);

	double result = result1;
	result -= result2/df_bespoke*df_index;

	//Remise a jour carac mezz
	bespoke->SetMezzAmount(init_mezz);
	bespoke->SetSubAmount(init_sub);
	traxx->SetMezzAmount(init_mezz_traxx);
	traxx->SetSubAmount(init_sub_traxx);

	//Free mem
	if (Correl_k2) delete Correl_k2;
	return (result);
}

// virtual 
ICM_Correlation* ICM_Smile_Correlation::GenerateShiftCorrel(const std::string& label,
											qSENSITIVITY_TYPE typesensi,
											double epsilon )
{
	int size = itsSlices.size();
	ICM_Smile_Correlation* Correlation = (ICM_Smile_Correlation*) Clone();
	ARM_Date date;

	double eps = epsilon;
	int nthLine,nthCol,NoIndex; 

	ARM_VolCurve* volcurve=Correlation->itsSlices[0].GetVolCurve();
	ARM_Matrix* Volatilities = volcurve->GetVolatilities();
	int nLines_0 = Volatilities->GetNumLines();
	int nCols_0  = Volatilities->GetNumCols();	

	if ((typesensi == ICMCORREL_STRIKE_DOWN_TYPE) && 
		(nLines_0==1) && (nCols_0==2) && //flat correlation case
		(itsSlices.size()==1))
		{nthLine=0;nthCol=0;NoIndex=0;eps=1.;}
	else if ((typesensi == ICMCORREL_STRIKE_UP_TYPE) && 
		(nLines_0==1) && (nCols_0==2) && //flat correlation case
		(itsSlices.size()==1))
		{nthLine=0;nthCol=1;NoIndex=0;eps=1.;}
	else Find(label,nthLine,nthCol,NoIndex); 

	if ((nthLine == -1) || (nthCol == -1) || (NoIndex == -1)) return Correlation;

	(Correlation->itsSlices[NoIndex].GetVolCurve())->BumpVolatility(eps,nthLine+1,nthCol+1);

	return (Correlation);
}


// --------------------------------------------------------------------
// ING Equivalent Strikes Computation Using ELoss stripping
// --------------------------------------------------------------------
void ICM_Smile_Correlation::ComputeStrikes_ING(ICM_Pricer* pricer)
{
	//Rescaling deja effectue
	if ((its_Already_rescal)||(its_never_rescal)) return;

	its_BaseCorrelShift=0.;

	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity()->Clone();
	cdo->GetFeeLeg()->SetCreditLegType(qRunning_Leg);
	cdo->GetDefLeg()->SetCreditLegType(qStandart_Recovery_Leg);

	cdo->SetPorS(1.);		// Cas cdo standart
	cdo->SetTradedCoef(1.); // Cas cdo standart

	int size = itsSlices.size();
	double K1 = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double K2 = cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate()) + K1;

	bool flagequity = false;
	if CHECK_NULL(K1) flagequity=true;

	double eloss_bespoke = 0.,eloss_indice = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	ARM_Date Maturity = cdo->GetDefLeg()->GetMaturity();

	//Eloss du portefeuille : sommes des eloss ind
	its_Callib_model= (ICM_ModelMultiCurves*) pricer->GetModel()->Clone();
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) its_Callib_model;
	ARM_ZeroFlat* zcflat= new ARM_ZeroFlat((ARM_Date)its_Callib_model->GetZeroCurve()->GetAsOfDate(),0.);
	mod->SetZeroCurve(zcflat);

	for (int j=0;j<mod->GetNbDefCurves();j++)
	{  
		ICM_DefaultCurve* dc = (ICM_DefaultCurve*) mod->GetDefaultCurve(j)->Clone();
		dc->SetZeroCurve((ARM_ZeroFlat*)zcflat->Clone());
		mod->SetDefaultCurve(j,dc);
	}

	if (zcflat) delete zcflat;zcflat=NULL;
		
	//Nominal total du portefeuille
	//nominal_total = cdo->GetCollateral()->SumNotionals();
	
	//calcul de l'EL du portefeuille
    eloss_bespoke = CptPtf_ELoss(pricer);

					
	//Strike Bespoke -> NB Eloss bespoke
	double StrikeUpEloss = 0.;
	double StrikeDownEloss = 0.;
	if (!CHECK_NULL(eloss_bespoke))
	{
		StrikeUpEloss = (cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate()) + cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate())) / fabs(eloss_bespoke);
		StrikeDownEloss = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate()) / fabs(eloss_bespoke);
	}

	double EL_Composite=0.;

	for (int NoIndex=0;NoIndex<size;NoIndex++)
	{
		bool isequal=false;
		double time,result_up = 0., result_down=0.;
		int idx_inf,idx_sup ;
		its_Callib_SizeEq_Index = 0.;

		itsSlices[NoIndex].GetSmileStrikeHigh().clear();
		itsSlices[NoIndex].GetSmileStrikeLow().clear();
		itsSlices[NoIndex].GetMaturities().clear();

		//Determination des bornes de projection en maturité
		ARM_Vector* V = ((ARM_VolLInterpol *)itsSlices[NoIndex].GetVolCurve())->GetExpiryTerms();
		vector<double> VX;VX.resize(V->GetSize());
		for (int ij=0;ij<V->GetSize();ij++) {VX[ij]=V->Elt(ij);}
		time=(Maturity-its_Callib_model->GetStartDate())/365.;

		Bornes(VX,time,idx_inf,idx_sup,isequal);

		//Borne Inf
		if (idx_inf!=CREDIT_DEFAULT_VALUE) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_inf]);

		//Borne Sup
		if ((isequal==false) && (idx_sup!=CREDIT_DEFAULT_VALUE)) itsSlices[NoIndex].GetMaturities().push_back(VX[idx_sup]);

		//only maturity
		if (its_Interpolation_At_Maturity)
		{
			itsSlices[NoIndex].GetMaturities().clear();
			itsSlices[NoIndex].GetMaturities().push_back(time);
		}
		
		for (int IndMaturity=0;IndMaturity<itsSlices[NoIndex].GetMaturities().size();IndMaturity++)
		{			
			if (itsSlices[NoIndex].GetProportion() != 0.)
			{
				//Maturité de l'indice de référence
				its_Callib_ActiveMaturity = itsSlices[NoIndex].GetMaturities()[IndMaturity];

				ARM_Date Maturity_2 = AsOf;
				Maturity_2.AddDays((int)(365.* its_Callib_ActiveMaturity+0.5));

		
				//Strike Equivalent	
				try
				{
					//Estimation de l'eloss sur l'indice (exprimée en %)	
					EL_Composite += itsSlices[NoIndex].GetProportion() *
							CptPtf_ELoss_Index(mod,itsSlices[NoIndex].GetIndex(),Maturity_2);
				}
				catch (...)
				{
					if (cdo) delete cdo;
					if (its_Callib_model) delete its_Callib_model;
					if (its_Callib_Index) delete its_Callib_Index;
					ICMTHROW(ERR_INVALID_ARGUMENT,"unable to rescale (maturity,index)=("<<its_Callib_ActiveMaturity<<","<<itsSlices[NoIndex].GetIndex()->GetIndexName()<<")");
					itsSlices.clear();
				}
				if (its_Callib_Index) delete its_Callib_Index;

			}

			itsSlices[NoIndex].GetSmileStrikeLow().push_back(K1);
			itsSlices[NoIndex].GetSmileStrikeHigh().push_back(K2);

		}
	}

	itsForcedStrikeType=qStrike_NONE;
	double Pindex_Idx = GetCompositeCorrel("ISSUER","ISSUER",its_Callib_ActiveMaturity,EL_Composite,its_Callib_ActiveMaturity);
	double Pindex_Besp = GetCompositeCorrel("ISSUER","ISSUER",its_Callib_ActiveMaturity,eloss_bespoke,its_Callib_ActiveMaturity);

	its_BaseCorrelShift = Pindex_Idx - Pindex_Besp;

	//Traces
	//fclose(stream);
	if (cdo) delete cdo;
	if (its_Callib_model) delete its_Callib_model;
	its_Already_rescal = true;
}
// virtual 
double ICM_Smile_Correlation::GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity ,
								  double strike ,
								  double actualYF )
	{

		if (itsRescalType==qRescal_ING_Maturity)
		{return its_BaseCorrelShift + GetCompositeCorrel(issuer1,issuer2,maturity,strike,actualYF);}

		double CheckValue = 0.;
		int k=0;
		int size_index = itsSlices.size();

		vector<int> Vsize_matu(size_index);
		vector<int> VSizeSmileStrikeLow(size_index);
		vector<int> VSizeSmileStrikeUp(size_index);
		bool AllMatuEmpty = true;
		bool AllSmileStrikeEmpty = true;

		for (k=0;k<size_index;k++)
		{
			Vsize_matu[k]= itsSlices[k].GetMaturities().size();
			VSizeSmileStrikeLow[k]= itsSlices[k].GetSmileStrikeLow().size();
			VSizeSmileStrikeUp[k]= itsSlices[k].GetSmileStrikeHigh().size();
		}

		for (k=0;k<size_index;k++)
		{if (Vsize_matu[k]!=0) {AllMatuEmpty=false;break;}}

		for (k=0;k<size_index;k++)
		{if (VSizeSmileStrikeLow[k]!=0) {AllSmileStrikeEmpty=false;break;}}

		//if (its_Already_rescal==false) {return CREDIT_DEFAULT_VALUE;}
		double UsualMaturity = maturity;
		if ((actualYF!=CREDIT_DEFAULT_VALUE)&&(its_TermStructureRescaling==qTermStructure)) 
			{UsualMaturity=actualYF;}

		if	(AllMatuEmpty && (!AllSmileStrikeEmpty))
		{
			for (k=0; k<size_index; k++)
			{
				Vsize_matu[k] = 1;
				itsSlices[k].GetMaturities().push_back(UsualMaturity);
			}
		}
		else if	(AllMatuEmpty) 
		{
			return CREDIT_DEFAULT_VALUE;
		}

		double Inter_Corr = 0.,Correlation = 0.,interpol=0.;
		

		double PrevMaturity = 0.;

		switch (itsForcedStrikeType)
		{
		case qStrike_NONE :
			{
				for (int i=0; i<size_index; i++)
					if (itsSlices[i].GetProportion()) 
						{
							vector<double> VCorrelation(Vsize_matu[i]) ; // VCorrelation.clear();
							for (int j=0; j<Vsize_matu[i]; j++)
							{
								Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(itsSlices[i].GetMaturities()[j],100.*strike)/100.;
								// VCorrelation.push_back(Inter_Corr);
								VCorrelation[j]=Inter_Corr;
							}
							interpol = VectorInterpol(itsSlices[i].GetMaturities(),VCorrelation,UsualMaturity,K_LINEAR);
							Correlation += itsSlices[i].GetProportion()*interpol;
						}
				break;
			}
		case qStrike_UP :
			{for (int i=0; i<size_index; i++)
				if (itsSlices[i].GetProportion()) 
				{
					// VCorrelation.clear();
					vector<double> VCorrelation(Vsize_matu[i]) ; // VCorrelation.clear();
				for (int j=0; j<Vsize_matu[i]; j++)
				{Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(itsSlices[i].GetMaturities()[j],100.*itsSlices[i].GetSmileStrikeHigh()[j])/100.;
				// VCorrelation.push_back(Inter_Corr);
				VCorrelation[j]=Inter_Corr;
				}
				interpol = VectorInterpol(itsSlices[i].GetMaturities(),VCorrelation,UsualMaturity,K_LINEAR);
				Correlation += itsSlices[i].GetProportion()*interpol;
				}
			break;}
		case qStrike_LOW :
			{for (int i=0; i<size_index; i++)
				if (itsSlices[i].GetProportion()) 
				{ // VCorrelation.clear();
					vector<double> VCorrelation(Vsize_matu[i]) ; // VCorrelation.clear();
				for (int j=0; j<Vsize_matu[i]; j++)
				{Inter_Corr = (itsSlices[i].GetVolCurve())->ComputeVolatility(itsSlices[i].GetMaturities()[j],100.*itsSlices[i].GetSmileStrikeLow()[j])/100.;
				 // VCorrelation.push_back(Inter_Corr);
				VCorrelation[j]=Inter_Corr;
				}
				interpol = VectorInterpol(itsSlices[i].GetMaturities(),VCorrelation,UsualMaturity,K_LINEAR);
				Correlation += itsSlices[i].GetProportion()*interpol;
				}
			break;}
		}

		if (AllMatuEmpty)
		{for (int ik=0; ik<size_index; ik++)
		{itsSlices[ik].GetMaturities().clear();}}

		return (Correlation);
	}
//
void ICM_Smile_Correlation::SetCorrelations(ICM_Smile_Correlation& correl)
	{
		if (correl.itsSlices.size()!=itsSlices.size())
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Smile_Correlation:: Incompatible size for SmileCorrelation");
			
		for (int il=0;il<correl.itsSlices.size();il++)
		{itsSlices[il].SetVolCurve(correl.itsSlices[il].GetVolCurve());}
	}
// 
void ICM_Smile_Correlation::SetInterpType(const int& interpType)
	{
		for (int j=0; j<itsSlices.size(); j++)
		{itsSlices[j].GetVolCurve()->SetInterpType(interpType);}
    }
// virtual 
int ICM_Smile_Correlation::GetInterpType() 
{return itsSlices[0].GetVolCurve()->GetInterpType();}

// virtual 
double ICM_Smile_Correlation::GetCompositeCorrel(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity ,
								  double strike ,
								  double actualYF )
	{
		int i=0;
		int size_index = itsSlices.size();
		vector<double> Maturities;
		vector<double> Strikes;

		for (i=0; i<size_index; i++)
		{
			vector<double> Maturities_;
			vector<double> Strikes_;

			vector<double> Maturities_2;
			vector<double> Strikes_2;

			ARM_Vector* strikes = ((ARM_VolLInterpol*)(itsSlices[i].GetVolCurve()))->GetStrikes();
			ARM_Vector* yearterms = ((ARM_VolLInterpol*)(itsSlices[i].GetVolCurve()))->GetExpiryTerms();

			vc2stdv(strikes->GetSize(),strikes->GetElt(),Strikes_);
			vc2stdv(yearterms->GetSize(),yearterms->GetElt(),Maturities_);

			MergeVector(Maturities,Maturities_,Maturities_2);
			MergeVector(Strikes,Strikes_,Strikes_2);
		
			Maturities = Maturities_2;
			Strikes = Strikes_2;
		}

		ARM_Matrix* Mat = new ARM_Matrix(Maturities.size(),Strikes.size(),0.);

		for (int k1=0;k1<Maturities.size();k1++)
			for (int k2=0;k2<Strikes.size();k2++)
			{
				for (int i=0; i<size_index; i++)
				{
					Mat->Elt(k1,k2) += itsSlices[i].GetProportion() *
						(itsSlices[i].GetVolCurve())->ComputeVolatility(Maturities[k1],Strikes[k2]);
				}
			}

// FIXMEFRED: mig.vc8 (28/05/2007 10:30:03):cast
			ARM_VolLInterpol* newvol = (ARM_VolLInterpol*) (itsSlices[0].GetVolCurve())->Clone();
		newvol->SetVolatilities(Mat);
		ARM_Vector* VStrikes = new ARM_Vector(Strikes.size(),&(*Strikes.begin()));
		newvol->SetStrikes(VStrikes);
		ARM_Vector* VMaturities= new ARM_Vector(Maturities.size(),&(*Maturities.begin()));
		newvol->SetExpiryTerms(VMaturities);

		double correl =0.;

		switch (itsForcedStrikeType)
		{
		case qStrike_NONE :
			{
			correl = newvol->ComputeVolatility(maturity,100.*strike)/100.;
			break;
			}
		case qStrike_UP :
			{
			correl = newvol->ComputeVolatility(maturity,100.*itsSlices[0].GetSmileStrikeHigh()[0])/100.;
			break;
			}
		case qStrike_LOW :
			{
			correl = newvol->ComputeVolatility(maturity,100.*itsSlices[0].GetSmileStrikeLow()[0])/100.;
			break;
			}
		}

		if (newvol) delete newvol; 
		newvol = NULL;

		return (correl);
	}