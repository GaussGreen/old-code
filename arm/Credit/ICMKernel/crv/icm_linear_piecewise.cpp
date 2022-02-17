

#include "ICMKernel\crv\ICM_linear_Piecewise.h"

  

//	-----------------------------------------------------------------------------------
ICM_Linear_Piecewise::ICM_Linear_Piecewise(const ICM_Linear_Piecewise&ref) : ICM_DefaultCurve(ref)
{}
//	-----------------------------------------------------------------------------------
ICM_Linear_Piecewise::ICM_Linear_Piecewise(void) 
{ 
	Init();
}

//	-----------------------------------------------------------------------------------
// 	virtual 
double 
ICM_Linear_Piecewise::DefProbInverse(const double& PriceIn) const 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Linear_Piecewise::DefProbInverse: not implemented"); 
}


//	-----------------------------------------------------------------------------------
// virtual     
/**void 
ICM_Linear_Piecewise::Copy(const ARM_Object* srcZc)
{
	ICM_DefaultCurve::Copy(srcZc);
	// BitwiseCopy(srcZc);
}
**/ 
//	-----------------------------------------------------------------------------------
// virtual     
ARM_Object* 
ICM_Linear_Piecewise::Clone(void)
{
	return new ICM_Linear_Piecewise(*this); 
}

//	-----------------------------------------------------------------------------------
// virtual 
void 
ICM_Linear_Piecewise::Calibrate() 
{
	ICM_DefaultCurve::Calibrate(); 
}
//	-----------------------------------------------------------------------------------
// virtual 
// void 
// ICM_Linear_Piecewise::Calibrate_Stress_Test_Guess_Brent() 
// {
// 	ICM_DefaultCurve::Calibrate(); 
// }

//	-----------------------------------------------------------------------------------
void 
ICM_Linear_Piecewise::Init(void)
{
	SetName(ICM_LINEAR_PIECEWISE);
}

//	-----------------------------------------------------------------------------------
// virtual 
double 
ICM_Linear_Piecewise::SurvivalFunction(const double& yearterm) const 
{
	if (GetIsNameInDefault()) return 0.;

	bool equal = false;
	int i = 0;
	// Inf_Equal(yearterm,i,equal);
	// Inf_EqualInterpol(yearterm,i,equal);
	i = position(yearterm); 
	if (eq(itsInterpolYF.Elt(i),yearterm)) 
	// if (equal) 
		// return GetSurvivalProba(i);
		return GetInterpolSP().Elt(i); // GetSurvivalProba(i);

	double SP_prev = GetInterpolSP().Elt(i); // GetSurvivalProba(i);
	double Lambda = GetInterpolLambda().Elt(i); // GetLambda(i);
	double Sproba = 0.;
	int last = GetInterpolDates().GetSize();

	// if ((i<1) || (yearterm > GetYearTerm(last-1)) )
	if ((i<1) || (yearterm > GetInterpolYF().Elt(last-1)) )
	 Sproba = SP_prev*exp(-Lambda* (yearterm - GetInterpolYF().Elt(i)));
	else 
	{
		double YT_prev = GetInterpolYF().Elt(i);
		double YT_next = GetInterpolYF().Elt(i+1);
		double Lambda_prev = GetInterpolLambda().Elt(i);
		double Lambda_next = GetInterpolLambda().Elt(i+1); 
		double Step_time= YT_next - YT_prev;
		double Coef_dir = (Lambda_next - Lambda_prev)/Step_time;
		double Cte		= -(Lambda_next * YT_prev - Lambda_prev * YT_next) / Step_time;
		double Lambda_courant = Coef_dir * yearterm + Cte;
		Sproba = SP_prev*exp(-(Lambda_courant + Lambda_prev) * (yearterm - YT_prev)/2.);
	}
	
	return (Sproba);

}

//	-----------------------------------------------------------------------------------
// virtual 
void 
ICM_Linear_Piecewise::CptTermsSurvivalProba(void)
{
	// int size = GetDates()->GetSize();
	int size = GetInterpolDates().GetSize();
	// GetSurvivalProba()->InitElt(0,1.);
	SetInterpolSP(0,1); 

	double SP_prev = 1.;
	double Lambda  = GetInterpolLambda().Elt(0); // GetLambda(0);
	double Sproba  = 0.;
	double YT_prev = 0.;
	double YT_next = GetInterpolYF().Elt(1); // GetYearTerm(1);
	
	Sproba = exp(-Lambda*YT_next);
	SP_prev = Sproba;
	// GetSurvivalProba()->InitElt(1,Sproba);
	SetInterpolSP(1,Sproba); 

	int i=0;

	for (i=1; i<size-1; i++)
	{
		YT_prev = GetInterpolYF().Elt(i); // GetYearTerms()->Elt(i);
		YT_next = GetInterpolYF().Elt(i+1); // GetYearTerms()->Elt(i+1);

		double Lambda_prev = GetInterpolLambda().Elt(i); // GetLambda(i);
		double Lambda_next = GetInterpolLambda().Elt(i+1); // GetLambda(i+1); 
		double Step_time= YT_next - YT_prev;
		Sproba = SP_prev*exp(-(Lambda_next + Lambda_prev) * (YT_next - YT_prev)/2.);
		SP_prev = Sproba;
	 
		// GetSurvivalProba()->InitElt(i+1,Sproba);
		SetInterpolSP(i+1,Sproba); 
	}	
}


//	-----------------------------------------------------------------------------------
// virtual 
void 
ICM_Linear_Piecewise::ResetLambda(const int& i,const double& lambda)
{
	// int size = GetDates()->GetSize();
	int size = GetInterpolDates().GetSize();

	double SP_prev = 1.;
	double YT = GetInterpolYF().Elt(i); // GetYearTerm(i);
	double YT_prev = 0.;
	double SP = 1.;
	double Lambda_prev = 0.;
	
	if (i == 1) 
	{	
		// GetLambdas()->InitElt(i-1,lambda);			
		// GetLambdas()->InitElt(i,lambda);			
		SetInterpolLambda(i-1,lambda); 
		SetInterpolLambda(i,lambda); 

		//Reset de la survival proba du plot actif
		SP = exp(-lambda *YT);
		// GetSurvivalProba()->InitElt(i,SP);
		SetInterpolSP(i,SP); 
		return;
	}		

	if (i>1)
	{	
		// GetLambdas()->InitElt(i,lambda);			
		SetInterpolLambda(i,lambda); 

		//Reset de la survival proba du plot actif
		SP_prev = GetInterpolSP().Elt(i-1); // GetSurvivalProba()->Elt(i-1);
		YT_prev = GetInterpolYF().Elt(i-1); // GetYearTerm(i-1);
		Lambda_prev = GetInterpolLambda().Elt(i-1); // GetLambda(i-1);
		SP = SP_prev * exp(-(lambda + Lambda_prev) * (YT - YT_prev)/2.);
		// GetSurvivalProba()->InitElt(i,SP);
		SetInterpolSP(i,SP); 
		return;
	}		
}


//	-----------------------------------------------------------------------------------
// virtual 
void 
ICM_Linear_Piecewise::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[40];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

   fprintf(fOut, "\n\n--------------------------------------------------------------------------------\n");
   fprintf(fOut, "---------------- Linear Piece Wise Default Curve -----------------------------\n");
   fprintf(fOut, "--------------------------------------------------------------------------------\n\n");

	ICM_DefaultCurve::View(id, fOut);

   	if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


//	-----------------------------------------------------------------------------------
ICM_Linear_Piecewise::ICM_Linear_Piecewise(const ARM_Date& asOf,
						 const vector<string>& terms,
						 ARM_Vector* rates,
						 double Recovery,
						 ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
						 qCDS_ADJ adj,  	
						 const std::string& ccy,
						 const std::string&  label,
						 int Lag,
						 bool issummitcurve
						 //2 ICM_DefaultCurve* defcurve
						 )
{

	Init();

	Set(asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label, Lag,issummitcurve);
		
}

//	-----------------------------------------------------------------------------------
void 
ICM_Linear_Piecewise::Set (const ARM_Date& asOf,
					   const vector<string>& terms,
                       ARM_Vector* rates,
					   double Recovery,
					   ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
					   qCDS_ADJ adj,  	
                       const std::string& ccy, 
					   const std::string&  label,
					   int Lag,
					   bool issummitcurve
					   //2 ICM_DefaultCurve* defcurv
					   )
{
	if (true)
		ICM_DefaultCurve::Set (asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label,issummitcurve,NULL,K_QUARTERLY,qDEFCURVE_DICHO,"STD", Lag,ICM_Parameters());
	else
		ICM_DefaultCurve::Set (asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label,issummitcurve,NULL,K_QUARTERLY,qDEFCURVE_DICHO,"STD", Lag,ICM_Parameters());
}



//	-----------------------------------------------------------------------------------
ICM_DefaultCurve* 
ICM_Linear_Piecewise::GenerateShiftCurve(const vector<string>& pTerm, 
																   const ARM_Vector& epsilon ,
																   qSENSITIVITY_TYPE mode ) const
{
   ICM_DefaultCurve* dc = NULL;
   bool check = false;	

   int szinp = epsilon.GetSize();
   int size = GetRates()->GetSize()-1;	

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
	   for (l=0; l<szinp; l++)
		for (k=0; k<size; k++)
		{
         if (pTerm[l] == GetTerm(k))
         {
				if (mode == ICMSPREAD_TYPE)
				{
					ModifiedRates->Elt(k) = ModifiedRates->Elt(k) + epsilon.Elt(l)/100.;
				}
				else if (mode == ICMSPRELSHIFT_TYPE)
					ModifiedRates->Elt(k) = ModifiedRates->Elt(k) *  (1.0 + epsilon.Elt(l));
			}
			check = true;
            k=ARM_NB_TERMS;
		}

		//if (!check) throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, "Invalid Plot");
		break;
	   }		   
   case ICM_DTR_TYPE : 
	   {
		   for (k=0; k<size; k++) ModifiedRates->Elt(k) = 20.;
		   break;
	   }
   default :
	   ICMTHROW(ERR_SQUARE_OR_SIZE_PB,"Unknown ShiftCurve Type"); 
	   // throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, "Unknown ShiftCurve Type");
   }

   dc = (ICM_DefaultCurve*) new ICM_Linear_Piecewise(GetAsOfDate(),
														GetTerms(),
														ModifiedRates,
														recovery,
														GetZeroCurve(),
														GetSTDCDS_Adjusted(),
														GetSTDCDS_AdjustedStartDateOnBusinessDay(),
														GetCdsAdj(),
														GetCurrency(),
														GetLabel(),
														GetLagStartDate(),
														GetIsSummitCurve());
				   
   if (ModifiedRates)
	delete ModifiedRates;
   ModifiedRates = NULL;

   return(dc);
}


//	-----------------------------------------------------------------------------------
/** 
ICM_DefaultCurve* 
ICM_Linear_Piecewise::xGenerateShiftCurve(double epsilon, qSENSITIVITY_TYPE mode) const 
{
	ARM_Vector* Vepsilon = new ARM_Vector(1,epsilon);

	ICM_DefaultCurve* dc = GenerateShiftCurve(GetTerms(),Vepsilon,mode );

	if (Vepsilon) delete Vepsilon;

    return(dc);
} **/ 