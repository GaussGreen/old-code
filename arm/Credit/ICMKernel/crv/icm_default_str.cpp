
#include "ICMKernel\crv\icm_default_str.h"
#ifndef unix 
#define ITOA _itoa
#else
#define ITOA intToStr
#endif
// ----------------------------
//	Copy of members data
// ----------------------------
// void ICM_DefCurvStr::BitwiseCopy(const ARM_Object* srccds)
// {
//     ICM_DefCurvStr* DefCurvStr = (ICM_DefCurvStr *)srccds;
// 	itsPente = DefCurvStr->itsPente;
// 	itsSpreadDate = DefCurvStr->itsSpreadDate;
// }

// -------------
//	Copy Method 
// -------------
// void ICM_DefCurvStr::Copy(const ARM_Object* srccds)
// {
//      ICM_DefaultCurve::Copy(srccds);
//  
//      BitwiseCopy(srccds);
// }

ICM_DefCurvStr::ICM_DefCurvStr(const ICM_DefCurvStr&ref) : ICM_DefaultCurve(ref)
{
	itsPente=ref.itsPente; 
	itsSpreadDate=ref.itsSpreadDate; 
}
// --------------
//	Clone Method
// --------------
ARM_Object* ICM_DefCurvStr::Clone(void)
{
     // ICM_DefCurvStr* theClone = new ICM_DefCurvStr();
     // theClone->Copy(this);
     // return(theClone);
	return new ICM_DefCurvStr(*this); 
}


// ---------------------
//	Init of members data
// ---------------------

void ICM_DefCurvStr::Init()
{
	itsPente = 0.;
	itsSpreadDate = 0.;
	SetIsSummitCurve(false);
}


// ------------------------------
//  Constructeurs & Set + complet
// ------------------------------
ICM_DefCurvStr::ICM_DefCurvStr(ICM_DefaultCurve* dc)
{
	Init();
	Set(dc->GetRates(),
		dc->GetYearTerms(),
		dc->GetRecovery(),
		dc->GetZeroCurve(),
		dc->GetLabel());
}

ICM_DefCurvStr::ICM_DefCurvStr( const ARM_Vector* Spread,
								const ARM_Vector* YF_plot,
								double Recovery,
								ARM_ZeroCurve* ZCurv,
								const std::string& label)
{
	Init();
	Set(Spread,
		YF_plot,
		Recovery,
		ZCurv,
		label);
}


void ICM_DefCurvStr::Set(const ARM_Vector* SpreadInit,
						 const ARM_Vector* YF_plot,
						 double Recovery,
						 ARM_ZeroCurve* ZCurv,
						 const std::string& label)
{
	int i = 0;
	ARM_Date AsOfDate;
	AsOfDate.Today();
	SetAsOfDate(AsOfDate);
	int sizeSpread	= SpreadInit->GetSize(); 
	int sizePlot	= YF_plot->GetSize();
	
	ARM_Vector Spread = ARM_Vector(sizeSpread + 1,1.);
	memcpy(Spread.GetElt()+1,SpreadInit->GetElt(),sizeof(double)*sizeSpread);
	Spread.Elt(0) = SpreadInit->Elt(0);

	ARM_Vector Plot = ARM_Vector(sizePlot + 1,1.);
	memcpy(Plot.GetElt()+1,YF_plot->GetElt(),sizeof(double)*sizePlot);
	Plot.Elt(0) = 0.;

	SetRates((ARM_Vector*)Spread.Clone());
	SetYearTerms((ARM_Vector*)Plot.Clone());
	SetRecovery(Recovery);

	// Les lambdas sont constants par morceaux et lambda final = lambda(final - 1)
	ARM_Vector Lambda = ARM_Vector(sizeSpread + 1,1.);
	for (i=0;i<sizeSpread;i++)
		Lambda.Elt(i) = GetRate(i+1)/(1-Recovery);
	
	Lambda.Elt(sizeSpread) = Lambda.Elt(sizeSpread-1);
	SetLambdas((ARM_Vector*)Lambda.Clone());

	SetZeroCurve((ARM_ZeroCurve*)ZCurv->Clone());
	SetLabel(label);

	CptSProbaPlot();
	Calibrate();

	// Term
	//char Terms[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	vector<string> Terms(YF_plot->GetSize());
	for(int j = 0; j < YF_plot->GetSize(); j++)
	{	
		char buffer[65];
		ITOA(YF_plot->Elt(j),buffer,10);
		string bufferS = string(buffer);
		Terms[j] = bufferS + string("Y");
	}

	SetTerms(Terms);
}

 
// ------------------------------------------------------
// Constructeur et Set avec un vecteur de date/YF 
// ------------------------------------------------------
ICM_DefCurvStr::ICM_DefCurvStr(const vector<double>& Spread,
							   const vector<double>& YForDates,
							   const double& Recovery,
							   ARM_ZeroCurve* ZCurv,
							   const std::string& Label,
							   const bool& IsDates)
{
	Init();

	Set(Spread,
		YForDates,
		Recovery,
		ZCurv,
		Label,
		IsDates);
}

void ICM_DefCurvStr::Set(const vector<double>& Spread,
						 const vector<double>& YForDates,
						 const double& Recovery,
						 ARM_ZeroCurve* ZCurv,
						 const std::string& Label,
						 const bool& IsDates)
{
	int SizeSpread = Spread.size();
	int SizeDate = YForDates.size();
	int size = SizeSpread;
	if (SizeSpread>SizeDate)
		size = SizeDate;

	// On transforme le std vecteur de spread en ARM_Vector* de Spread
	ARM_Vector* Rates = new ARM_Vector(size,1.);
	for (int j=0;j<size;j++)
		Rates->Elt(j) = Spread[j];

	// on transforme les date Excel en YF s'il le faut
	double Convention = 365.;
	ARM_Date AsOf;
	AsOf.Today();
	ARM_Vector* YF_Plot = new ARM_Vector(size,1.); 
	if (IsDates)
	{
		for (int i=0;i<size;i++)
			YF_Plot->Elt(i) = YFfromDate(AsOf,YForDates[i],Convention);
	}
	else 
		for (int i=0;i<size;i++)
			YF_Plot->Elt(i) = YForDates[i];



	Set(Rates,YF_Plot,Recovery,ZCurv,Label);
}



// ------------------
// Survival Fonction 
// ------------------
inline double ICM_DefCurvStr::SurvivalFunction(const double& yearterm) const
{
	double result = 0.;
	double time = 0.;
	double Lamb = 0.;
	
	if (yearterm < 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "YearTerm :  YearTerm < 0");
	}

	int indice = Find_InfEqual(yearterm,GetYearTerms());
	
	if (indice == CREDIT_DEFAULT_VALUE)
	{
		Lamb = GetLambda(0);
		result = exp(-(yearterm)*Lamb);
	}
	else
	{
		time = GetYearTerm(indice);
		Lamb = GetLambda(indice);
		result = GetSurvivalProba(indice) * exp(-(yearterm-time)*Lamb);
	}
	
	return result;
}
// ------------------------
// Cpt Survival Proba Plot 
// ------------------------
void ICM_DefCurvStr::CptSProbaPlot()
{
	int i = 0;
	double delta_YF = 0.;

	int size = GetYearTerms()->GetSize();
	double Lambda_temp = GetLambda(0);
	double YF_Temp = GetYearTerm(1);

	// On ajoute un element pour préciser qu'en 0 la proba de survie est de 1
	ARM_Vector Proba_Temp = ARM_Vector(size,1.);	
	Proba_Temp.Elt(0) = 1.;
	Proba_Temp.Elt(1) = exp(-YF_Temp * Lambda_temp);
		
	for (i=1;i<size;i++)
	{ 
		Lambda_temp = GetLambda(i-1);
		delta_YF = GetYearTerm(i) - GetYearTerm(i-1);
		Proba_Temp.Elt(i) = Proba_Temp.Elt(i-1) * exp(-delta_YF*Lambda_temp);;
	}

	SetSurvivalProba((ARM_Vector*)Proba_Temp.Clone());
}
// -------------------------------------------------------------
// Fonction Calibrate 
// -------------------------------------------------------------
// virtual 
void ICM_DefCurvStr::Calibrate(/** ICM_DefCurvStr* defcurve**/ )
{
#ifdef _DEBUG
	FILE * pFile = NULL;
	pFile = fopen("C:\\temp\\testCalibrate.txt","w");
#endif
	int i = 0;
	double Lambda_temp = 0.;
	double plot = 0.;
	double Spread = 0.;
	double NPVbefore, NPVafter, Derivee = 0.;
	double Pas = 0.25; 
	double Recov = GetRecovery();
	double tol = 1E-8;
	double test = 0.;
	double pNew = 1E-6;
	int IterMax = 20;
	int NbIter = 0;
	int size = GetYearTerms()->GetSize();
	ARM_Vector Lambda = ARM_Vector(GetLambdas());

#ifdef _DEBUG
	fprintf(pFile,"\n original spread :\n");
	((ARM_Vector*)GetRates())->View("",pFile);
	Lambda.View("",pFile);
	fprintf(pFile,"\n original Lambda :\n");
	Lambda.View("",pFile);
	fprintf(pFile,"\n original Survival proba :\n");
	((ARM_Vector*)GetSurvivalProba())->View("",pFile);
#endif
	
	for (i=0;i<size-1;i++) 
	{
		test = 10.;
		plot   = GetYearTerm(i+1);
		Spread = GetRate(i+1);

		NPVbefore = NPVCds(this, plot, Spread, 0.25);
#ifdef _DEBUG
	fprintf(pFile," \t i : %d\n", i);
	fprintf(pFile," \tPlot : %lf\n", plot);
	fprintf(pFile," \tNPVbefore : %lf\n", NPVbefore);
#endif
		NbIter = 0;

		while ( (fabs(test) > tol) && (NbIter < IterMax))
		{
			if (fabs(test - 10)>tol) {NPVbefore = test;}

			Lambda.Elt(i) += pNew;
			SetLambdas((ARM_Vector*)Lambda.Clone());
			CptSProbaPlot();

			NPVafter = NPVCds(this, plot, Spread, 0.25);
		
			Derivee = (NPVafter - NPVbefore)/pNew;

			Lambda.Elt(i) -= pNew;
			Lambda.Elt(i) -= NPVbefore/Derivee;
#ifdef _DEBUG
	fprintf(pFile," \tnew Survival proba :\n");
	((ARM_Vector*)GetSurvivalProba())->View("",pFile);
	fprintf(pFile,"\n \tNPVafter : %lf\n", NPVafter);
	fprintf(pFile," \tDerivee : %lf\n", Derivee);
	fprintf(pFile," New Lambda :\n");
	Lambda.View("",pFile);
#endif

			SetLambdas((ARM_Vector*)Lambda.Clone());
			CptSProbaPlot();

				
			test = NPVCds(this, plot, Spread, 0.25);
#ifdef _DEBUG
	fprintf(pFile,"\n \tnew Survival proba :\n");
	((ARM_Vector*)GetSurvivalProba())->View("",pFile);
	fprintf(pFile,"\n \t test = NPVCds : %lf\n", test);
#endif
			NbIter ++;
		}

		if (NbIter == IterMax)
		{
#ifdef _DEBUG
	if (pFile) fclose(pFile);
#endif
			
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Maximum number of iterations exceeded");
		}
	}

	Lambda.Elt(size-1) = Lambda.Elt(size-2);
	SetLambdas((ARM_Vector*)Lambda.Clone());

	CptSProbaPlot();

#ifdef _DEBUG
	if (pFile) fclose(pFile);
#endif
	if (itsYearTerms ) itsInterpolYF = *itsYearTerms;	
	if (itsLambdas) itsInterpolLambda=*itsLambdas; 
	if (itsDates) itsInterpolDates=*itsDates; 
	if (itsRates) itsInterpolRates=*itsRates; 
	if (itsSurvivalProba) itsInterpolSP = *itsSurvivalProba; // first element is
}			


// ----------------------------------------------
// Nouvelle Version ConstantPiecewise compliant
// ----------------------------------------------
ICM_DefaultCurve* ICM_DefCurvStr::GenerateShiftCurve(const vector<string>& pTerm, 
																   const ARM_Vector& epsilon ,
																   qSENSITIVITY_TYPE mode ) const 
{
   ICM_DefaultCurve* dc = NULL;
   bool check = false;	

   int szinp = epsilon.GetSize();
   int sizeRate = GetRates()->GetSize()-1;	

   ARM_Vector* ModifiedRates = new ARM_Vector(sizeRate,0.);
   memcpy(ModifiedRates->GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* sizeRate);
   double recovery = GetRecovery();

   int k=0,l=0;
 
   if (mode == ICMRECOVERY_TYPE)
	   recovery += epsilon.Elt(0);

   if (pTerm.size()>0 && pTerm[0]=="ALL") 
   {
		if (mode==ICMSPRELSHIFT_TYPE) for (k=0; k<sizeRate; k++) ModifiedRates->Elt(k) *= (1.0 + epsilon.Elt(0));
		else if (mode==ICMSPREAD_TYPE) for (k=0; k<sizeRate; k++) ModifiedRates->Elt(k) += epsilon.Elt(0)/100.;
   }
   else 
   {
	   for (l=0; l<szinp; l++)
	   {
		  for (k=0; k<sizeRate; k++)
		  {
			  if (pTerm[l] == GetTerm(k) )
			 {
				if (mode != ICMRECOVERY_TYPE)
				{
					if (mode != ICMSPRELSHIFT_TYPE)
						ModifiedRates->Elt(k) += epsilon.Elt(l)/100.;
					else
						ModifiedRates->Elt(k) *= (1.0 + epsilon.Elt(l));
				}
				check = true;
				k=ARM_NB_TERMS;
			 }
		  }
	   }
   }
  
   if (!check) 
	   throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, "Invalid Plot");


   int sizeYearTerm = GetYearTerms()->GetSize() - 1;
   ARM_Vector* PlotYearFraction = new ARM_Vector(sizeYearTerm, 0.);
   memcpy(PlotYearFraction->GetElt(),GetYearTerms()->GetElt()+1 ,sizeof(double)* sizeYearTerm);

  
   dc = (ICM_DefaultCurve*) new ICM_DefCurvStr(ModifiedRates,
											   PlotYearFraction,
											   recovery,
											   GetZeroCurve(),
											   GetLabel());
	
   // Delete
   if (ModifiedRates)
	delete ModifiedRates;
   ModifiedRates = NULL;

   if (PlotYearFraction)
	   delete PlotYearFraction;
   PlotYearFraction = NULL;

   return(dc);
}



/** JLA Removed 
ICM_DefaultCurve* ICM_DefCurvStr::xGenerateShiftCurve(double epsilon, 
															 qSENSITIVITY_TYPE mode) const
{
   ICM_DefaultCurve* dc = NULL;
   bool check = false;	

   int sizeRate = GetRates()->GetSize()-1;	

   ARM_Vector* ModifiedRates = new ARM_Vector(sizeRate,0.);
   memcpy(ModifiedRates->GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* sizeRate);
   double recovery = GetRecovery();

   if (mode != ICMRECOVERY_TYPE)
   {
	   if (mode != ICMSPRELSHIFT_TYPE)
	   {
		   for (int l=0; l<sizeRate; l++)
		      ModifiedRates->Elt(l) += epsilon/100.;
	   }
	   else
		   for (int l=0; l<sizeRate; l++)
		      ModifiedRates->Elt(l) *= (1. + epsilon);
   }
  else
	recovery += epsilon;

   int sizeYearTerm = GetYearTerms()->GetSize() - 1;
   ARM_Vector* PlotYearFraction = new ARM_Vector(sizeYearTerm, 0.);
   memcpy(PlotYearFraction->GetElt(),GetYearTerms()->GetElt()+1 ,sizeof(double)* sizeYearTerm);

  dc = (ICM_DefaultCurve*) new ICM_DefCurvStr(ModifiedRates,
											   PlotYearFraction,
											   recovery,
											   GetZeroCurve(),
											   GetLabel());

   // Delete   
   if (ModifiedRates)
	delete ModifiedRates;
   ModifiedRates = NULL;

   return(dc);

}
**/ 

// ------------------------------------------------------------------
// Fonction convertissant une DefCurveStr en ConstantPiecewise Curve
// ------------------------------------------------------------------
ICM_Constant_Piecewise* ICM_DefCurvStr::ConvertStrToARM(void)
{
	ARM_Date asOf;
	asOf.Today();

	int sizeYearTerm = GetYearTerms()->GetSize() - 1;
	double* YearFractions = new double[sizeYearTerm];
	for (int i = 0; i<sizeYearTerm;i++)
		YearFractions[i] = GetYearTerms()->GetElt()[i+1];

	int sizeIntensity = GetLambdas()->GetSize()-1;

	ARM_Vector Intensity = ARM_Vector(sizeIntensity);
	for (i=0;i<sizeIntensity;i++)
		Intensity.Elt(i) = GetLambda(i);

	double Recovery = GetRecovery();
	ARM_ZeroCurve* zc = (ARM_ZeroCurve*) GetZeroCurve()->Clone();

	ICM_Constant_Piecewise* Curve = new ICM_Constant_Piecewise(asOf,YearFractions,&Intensity,1, Recovery,zc,
														GetSTDCDS_Adjusted(),
														GetSTDCDS_AdjustedStartDateOnBusinessDay(),
														GetCdsAdj(),
														GetCurrency(),GetLabel(),
														NULL,
														GetCalibrationAlgo(),
														GetCalibrationData(),
														GetLagStartDate());
	Curve->SetRates((ARM_Vector*)unconst(*this->GetRates()).Clone());
	Curve->SetIsSummitCurve(false);


	// Delete
	if (YearFractions)
		delete[] YearFractions;
	YearFractions = NULL;

	if (zc)
		delete zc;
	zc = NULL;



	return Curve;
}
// -------------------------------------------------------------
// Computes Forward Spread between 2 dates. Return Fwd en bp
// -------------------------------------------------------------
double ICM_DefCurvStr::FwdSpread(const ARM_Date& t1, const ARM_Date& t2) const 
{
	double Forward = 0.;
	ARM_Date AsOf;
	AsOf.Today();
	double DateInit = (t1.GetJulian() - AsOf.GetJulian())/365.25;
	double Duree = (t2.GetJulian() - t1.GetJulian())/365.25;
	Forward = SpreadForward(this,DateInit,Duree);
	Forward *= 10000.;

	return Forward;
}

/** 
double ICM_DefCurvStr::FwdSpread(const double& t1, const double& t2)
{
	double Forward = 0.;
	ARM_Date AsOf;
	AsOf.Today();
	double DateInit = t1;
	double Duree = t2 - t1;
	Forward = SpreadForward(this,DateInit,Duree);
	Forward *= 10000.;

	return Forward;
}
**/ 
// ------------------------------------------------------------------
// View File
// ------------------------------------------------------------------
void ICM_DefCurvStr::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Discount Curve For Str----------------- \n");


   fprintf(fOut, "\n\nYearTerms\tRate\t\tIntensity\tSurvivalProba\n");

   for (int i = 0; i < GetRates()->GetSize(); i++)
   {	
	fprintf(fOut, "\n");
	fprintf(fOut, "%f\t", GetYearTerm(i));
	fprintf(fOut, "%f\t", GetRate(i));
	fprintf(fOut, "%f\t", GetLambda(i));
	fprintf(fOut, "%f\t", GetSurvivalProba(i));
   }	
    
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


// -------------------------------------------------------------------------------------------------------------------
// virtual 
ICM_DefaultCurve* ICM_DefCurvStr::GenDefCurve(const ARM_Vector& dates,
									  const ARM_Vector& spreads,
									  const double& recovery,
									  const std::string& label,
									  const bool& isdates,
									  ARM_ZeroCurve* ircurve) const 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurvStr::GenDefCurve: NOT IMPLEMENTED"); 
}
