
#include "ARMKernel\glob\firsttoinc.h" 

#include "ICMKernel\crv\icm_defcurve_index_pwc.h"

#include "ICMKernel\pricer\icm_pricer_cds.h"
#include "ICMKernel\util\icm_brentsolver.h"

#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"


// ---------------------------------------------------------------------
ICM_DefcurveIndex::ICM_DefcurveIndex(const ICM_DefcurveIndex&ref) : ICM_Constant_Piecewise(ref)
{
	itsRefSpread=ref.itsRefSpread; 
	itsSltPricer=ref.itsSltPricer; 
}
// ---------------------------------------------------------------------
// virtual 
ARM_Object* 
ICM_DefcurveIndex::Clone()
{
	return new ICM_DefcurveIndex(*this); 
}
// ---------------------------------------------------------------------
// Evaluate NPV of default pricer & default Index for an intensity value
// ---------------------------------------------------------------------
double ICM_DefcurveIndex::Evaluate(const double& x)
{
	int indx = GetCurrent_indice();

	// Variation du lambda
	ResetLambda(indx, x );

	// NPV before variation
	Get_pricer()->ResetPricer();
	double npv = (*Get_pricer()).Price(qCMPPRICE);

	npv -= (*itsSltPricer).Price(qCMPPRICE);
	
	return (npv);
}

// ---------------------------------------------------------------------
// Calibration for cds index curve
// ---------------------------------------------------------------------
void ICM_DefcurveIndex::Calibrate()
{

	switch (itsCalibrationAlgo)
	{
	case qDEFCURVE_DICHO: 
	case qDEFCURVE_NEWTON:
		break; // continue 
	default:
		ICM_DefaultCurve::Calibrate(); 
	}


	if ((GetRecovery() > 1.0) || (GetRecovery() < 0.0))
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Recovery must be comrpised between [0 ; 100%]!") ;
	}
	if (GetCalibrationData()!="STD") 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefcurveIndex: can't support calibrationMethod "<<GetCalibrationData()); 


	//	JLA 
	//
	//		Calibration code is working over itsInterpolXXX vectors
	itsInterpolYF = *itsYearTerms;	
	itsInterpolLambda=*itsLambdas; 
	itsInterpolDates=*itsDates; 
	itsInterpolRates=*itsRates; 
	itsInterpolSP = *itsSurvivalProba; // first element is 1


	//
	//
	ARM_Date AsOf = GetAsOfDate();
	double Nominal = DEFAULT_NOMINAL;
	double rate = GetInterpolRates().Elt(1);
	double Refrate = itsRefSpread[1];
	double result = 1.;


	// Init Lambda
	double Lambda=0.;
	Lambda	=	rate;
	if (1.0-GetRecovery())
		{Lambda	/=	(1.0-GetRecovery());}

	// GetLambdas()->Elt(0)=Lambda;
	SetInterpolLambda(0,Lambda); 
	ICM_DefaultCurveModel* model = NULL;

	model = new ICM_DefaultCurveModel(this,GetZeroCurve(),NULL,false);

	double _inf = 0.;
	double _sup = 0.;

	int	nb_iter=0;

	// En sortie de cette boucle on a calibre tous les lambdas sauf le dernier
	// for (int it=1; it<GetDates()->GetSize(); it++)
	for (int it=1; it<GetInterpolDates().GetSize(); it++)
	{
		// GetCurrent_indice() = it;
		SetCurrent_indice(it); 
	
		if (GetIsNameInDefault()) {ResetLambda(GetCurrent_indice(),1000.);continue;}

		ARM_Date EndDate = GetInterpolDates().Elt(it); 
 		rate = GetInterpolRates().Elt(it); // GetRate(it);
		Refrate = itsRefSpread[it]; // ca marche car mêmes dimensions.. 

		ARM_Date FstCpnEffDate;
		PreviousRegularDate(GetSTDCDS_Frequency(), EndDate,GetSTDCDS_StartDate(),FstCpnEffDate);

 		ICM_Cds cds(GetSTDCDS_StartDate(),
							   EndDate,
							   0,
							   &FstCpnEffDate,
							   GetSTDCDS_ProtectionStartDate(),
							   EndDate,
							   Refrate,
							   GetSTDCDS_Notional(),GetSTDCDS_Notional(),	// premium/default fixe
							   GetSTDCDS_Frequency(),
							   GetSTDCDS_Basis(),
							   GetSTDCDS_AccOnDef(),
							   GetSTDCDS_Currency(),	//GetSTDCDS_Currency() required because of the const
							   GetSTDCDS_Stub(),
							   GetSTDCDS_CreditLag(),
							   -1,
							   GetSTDCDS_Adjusted(),
							   GetSTDCDS_IncludeMaturity(),
							   GetSTDCDS_AdjustedStartDateOnBusinessDay(),
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
							   qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
							   );	


		ICM_Pricer_Cds pricer; pricer.Set(&cds,model,ICM_Parameters(),GetAsOfDate());

		pricer.SetFaster(true);
		Get_pricer() = &pricer;

		//vector of yearterms
		ARM_Vector YT;
		YT.push_back((EndDate-AsOf)/365.);

		ARM_Vector Spreads;
		Spreads.push_back(rate);
	
		//definition for the flat curve
		double Recovery = GetRecovery();
		ICM_Constant_Piecewise FLATDC(AsOf,YT,Spreads,Recovery,GetZeroCurve(),0,0,
				/*qCredit_Special_None_YF*/ qCredit_Default, // useless 
				GetSTDCDS_Currency(),"DUMMY",true,NULL,false,//2NULL
				K_QUARTERLY,GetCalibrationAlgo(),
				GetCalibrationData(), GetLagStartDate());

		ICM_DefaultCurveModel model(dynamic_cast<ICM_DefaultCurve *>(&FLATDC),GetZeroCurve(),NULL,false);
		// ICM_Pricer_Cds sltpricer(&cds,&model);
		ICM_Pricer_Cds sltpricer; sltpricer.Set(&cds,&model,ICM_Parameters(),GetAsOfDate());

		itsSltPricer = &sltpricer;

		_inf = Lambda-1.E-2;
		_sup = Lambda+1.E-2;

		result = RootFinder1D(ff1::mem_call(&ICM_DefcurveIndex::Evaluate,(*this))).ZeroBracket(_inf,_sup);	
		if (!result) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::Calibrate: can't bracket "
			<< GetLabel()
			<< " from inf="<<_inf<<" sup="<<_sup); 
		
		// since we can't check "dichotomy failure" vs "other exception" we catch any exception and consider this a failure
		Lambda=-1; 
		
		try {
			switch(itsCalibrationAlgo)
			{
			case qDEFCURVE_DICHO:
				Lambda = RootFinder1D(ff1::mem_call(&ICM_DefcurveIndex::Evaluate,(*this))).Dichotomy(_inf,_sup,nb_iter,100,1.E-5,1.E-2);
				break; 
			case qDEFCURVE_NEWTON:
				Lambda = RootFinder1D(ff1::mem_call(&ICM_DefcurveIndex::Evaluate,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,0.01,1.E-5,1.E-2);
				break ;
			default:
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't calibrate with "<<itsCalibrationAlgo); 
			}
			
			if (Lambda<0.0) 
				ICMLOG("ICM_DefaultCurve::Calibrate failed: "<< GetLabel()<<" "<<EndDate<<", Lambda="<<Lambda); 
		}
		catch (std::exception&e)
		{
			ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<GetLabel()<<" "<<EndDate<<",Lambda="<<Lambda
				<<"["<<e.what()<<"]"); 
		}
		catch (Exception&e)
		{
			ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<GetLabel()<<" "<<EndDate<<",Lambda="<<Lambda
				<<"["<<e.GetMessage()<<"]"); 
		}
		catch(...) {
			ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<GetLabel()<<" "<<EndDate<<",Lambda="<<Lambda
				<<"[unknown exception]"); 
		}


		ResetLambda(GetCurrent_indice(),Lambda);

	}

	if (model) 
		delete model;
	model = NULL;

	//	At the end of the calibration step :
	//
	*itsLambdas = itsInterpolLambda	; 
	*itsSurvivalProba = itsInterpolSP; 
}

ICM_DefaultCurve* ICM_DefcurveIndex::GenerateShiftCurve(const vector<string>& pTerm, 
																   const ARM_Vector& epsilon ,
																   qSENSITIVITY_TYPE mode ) const
{
   ICM_DefaultCurve* dc = NULL;
   bool check = false;	

   int szinp = epsilon.GetSize();
   int size = GetRates()->GetSize()-1;	

   ARM_Vector Refrates(itsRefSpread.size()-1);
   for (int i=1;i<itsRefSpread.size();i++) {Refrates.Elt(i-1) = itsRefSpread[i];}

   ARM_Vector* ModifiedRates = new ARM_Vector(size,0.);
   memcpy(ModifiedRates->GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* size);
   double recovery = GetRecovery();

   int k=0,l=0;
 
   switch (mode)
   {
   case ICMRECOVERY_TYPE :
	   {
		recovery += epsilon.Elt(0);
		break;
	   }
   case ICMSPRELSHIFT_TYPE : 	
   case ICMSPREAD_TYPE : 		
	   {
		   if (pTerm.size()>0 && pTerm[0]=="ALL") 
		   {
			   if(mode == ICMSPREAD_TYPE) for (k=0; k<size; k++) ModifiedRates->Elt(k) += epsilon.Elt(0)/100.;
			   else if (mode == ICMSPRELSHIFT_TYPE) for (k=0; k<size; k++) ModifiedRates->Elt(k) *=  (1.0 + epsilon.Elt(0)); 
		   }
		   else 
		   {
			   for (l=0; l<szinp; l++)
				for (k=0; k<size; k++)
				{
				 if (pTerm[l] == GetTerm(k))
				 {
					if (mode == ICMSPREAD_TYPE)
						ModifiedRates->Elt(k) = ModifiedRates->Elt(k) + epsilon.Elt(l)/100.;
					else if (mode == ICMSPRELSHIFT_TYPE)
						ModifiedRates->Elt(k) = ModifiedRates->Elt(k) *  (1.0 + epsilon.Elt(l));

					check = true;
					k=ARM_NB_TERMS;
				 }
				}
		   }

		//if (!check) throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, "Invalid Plot");
		break;
	   }		   
   case ICM_DTR_TYPE : 
	   {
		   if (epsilon.Elt(0) != CREDIT_DEFAULT_VALUE) recovery = epsilon.Elt(0);
		   for (k=0; k<size; k++) ModifiedRates->Elt(k) = 20.;
		   break;
	   }
   
   }


   if (GetTerms().empty()) {
		ICMTHROW(ERR_SQUARE_OR_SIZE_PB, "Unable to bump curve : Invalid constructor");
	} else {
		 dc = (ICM_DefaultCurve*) new ICM_DefcurveIndex(GetAsOfDate(),
														GetTerms(),
														ModifiedRates,
														&Refrates,
														recovery,
														GetZeroCurve(),
														GetSTDCDS_Adjusted(),
														GetSTDCDS_AdjustedStartDateOnBusinessDay(),
														GetCdsAdj(),
														GetCurrency(),
														GetLabel(),
														GetIsSummitCurve(),
														//2 NULL,
														GetVolCurve(),
														GetPayFrequency(),
														GetCalibrationAlgo(), //shoule be get
														GetCalibrationData(),
														GetLagStartDate());
	
   }
				   
   if (ModifiedRates)
	delete ModifiedRates;
   ModifiedRates = NULL;

   return(dc);
}

/** JLA Removed 
ICM_DefaultCurve* ICM_DefcurveIndex::xGenerateShiftCurve(double epsilon, 
															 qSENSITIVITY_TYPE mode) const
{
   ICM_DefaultCurve* dc = NULL;
   bool check = false;	

   ARM_Vector Refrates(itsRefSpread.size()-1);
   for (int i=1;i<itsRefSpread.size();i++) {Refrates.Elt(i-1) = itsRefSpread[i];}	

   int size = GetRates()->GetSize()-1;	

   ARM_Vector* ModifiedRates = new ARM_Vector(size,0.);
   memcpy(ModifiedRates->GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* size);
   double recovery = GetRecovery();

   if (mode != ICMRECOVERY_TYPE)
   {
	   if (mode != ICMSPRELSHIFT_TYPE)
	   {
		   for (int l=0; l<size; l++)
		   {
			   ModifiedRates->Elt(l) = ModifiedRates->Elt(l) + epsilon/100.;
		   }
	   }
	   else
		   for (int l=0; l<size; l++)
		   {
			   ModifiedRates->Elt(l) = ModifiedRates->Elt(l) *(1. + epsilon);
		   }
   }
  else
  {
	recovery += epsilon;
  }	

   
  if ( GetTerms().empty()) {
		ICMTHROW(ERR_SQUARE_OR_SIZE_PB, "Unable to bump curve : Invalid constructor");
	} else {		
		dc = (ICM_DefaultCurve*) new ICM_DefcurveIndex(GetAsOfDate(),
														GetTerms(),
														ModifiedRates,
														&Refrates,
														recovery,
														GetZeroCurve(),
														GetSTDCDS_Adjusted(),
														GetSTDCDS_AdjustedStartDateOnBusinessDay(),
														GetCdsAdj(),
														GetCurrency(),
														GetLabel(),
														GetIsSummitCurve(),
														//? NULL,
														GetVolCurve(),
														GetPayFrequency(),
														GetCalibrationAlgo(),
														GetCalibrationData(),
														GetLagStartDate());
   }

   if (ModifiedRates)
	delete ModifiedRates;
   ModifiedRates = NULL;

   return(dc);
}
**/ 