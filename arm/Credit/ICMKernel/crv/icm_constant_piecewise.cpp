
#include "ICMKernel\crv\ICM_Constant_Piecewise.h"
#include "ICMKernel/util/icm_utils.h"


ICM_Constant_Piecewise::ICM_Constant_Piecewise(const ICM_Constant_Piecewise& ref) : ICM_DefaultCurve(ref) 
{
	// non members.
}

void ICM_Constant_Piecewise::View(char* id, FILE* ficOut)
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
   fprintf(fOut, "---------------- Constant Piece Wise Default Curve -----------------------------\n");
   fprintf(fOut, "--------------------------------------------------------------------------------\n\n");

	ICM_DefaultCurve::View(id, fOut);

   	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

ICM_Constant_Piecewise::ICM_Constant_Piecewise(const ARM_Date& asOf,
												const ARM_Vector& Dates,
												const ARM_Vector& Inputs,
												int Intensity_or_Survivalprobas,
												double Recovery,
												ARM_ZeroCurve* zc,  
												int intRule,
												int adjStartDate,	
												qCDS_ADJ adj ,  
												const std::string& ccy,
												const std::string& label, 
												const ARM_VolCurve* VolCurve,
												// bool	isBrentCalib,
												qDEFCURVE_CALIB_ALGO calibAlgo,
												const std::string& calibData,
												int LagStartDate,
												const ICM_Parameters& params)
{
	Init();

	SetAsOfDate(asOf) ;

	ARM_Vector* Intensity = 0; // JLA new ARM_Vector(Inputs->GetSize(),0.) ;

	if (Intensity_or_Survivalprobas == -1) // Inputs are survival probabilities
		Intensity = CptLambdas(Inputs,Dates);
	else	// Inputs are hazard rates
		Intensity = (ARM_Vector*)unconst(Inputs).Clone();
		
	ICM_DefaultCurve::Set(asOf,Dates,*Intensity, Recovery,zc,intRule,adjStartDate,adj,
					ccy,label,VolCurve, calibAlgo,calibData, LagStartDate,params);

	if(Intensity)
		delete Intensity;
	Intensity = NULL;
}


ICM_Constant_Piecewise::ICM_Constant_Piecewise(const ARM_Date& asOf,
												double* YearFractions,
												ARM_Vector* Inputs,
												int Intensity_or_Survivalprobas,
												double Recovery,
												ARM_ZeroCurve* zc,
												int intRule,
												int adjStartDate,	
												qCDS_ADJ adj ,  
												const std::string& ccy,
												const std::string&  label,
												const ARM_VolCurve* VolCurve,
												qDEFCURVE_CALIB_ALGO	calibAlgo,
												const std::string& calibData,
												int LagStartDate)
{
	Init();

	SetAsOfDate(asOf) ;

	ARM_Vector* Intensity = NULL;

	if (Intensity_or_Survivalprobas == -1) // Inputs are survival probabilities
		Intensity = CptLambdas(*Inputs,YearFractions);
	else	// Inputs are hazard rates
		Intensity = (ARM_Vector*)Inputs->Clone();
	
	ICM_DefaultCurve::Set(asOf,YearFractions,Intensity, Recovery,zc,intRule,
				adjStartDate,adj,ccy,label,VolCurve, calibAlgo,calibData,LagStartDate);

	if(Intensity)
		delete Intensity;
	Intensity = NULL;
}

ICM_Constant_Piecewise::ICM_Constant_Piecewise(const ARM_Date& asOf,
						 const vector<string>& terms,
						 ARM_Vector* rates,
						 double Recovery,
						 ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,												
						 qCDS_ADJ adj,  	
						 const std::string& ccy,
						 const std::string&  label,
						 bool issummitcurve,
						 const ARM_VolCurve* VolCurve,
						 long PayFreq,
						 qDEFCURVE_CALIB_ALGO	calibAlgo,
						 const std::string &calibData,
						 int LagStartDate,
						 const ICM_Parameters& params,
						 ARM_Vector* UpFronts,
					     vector<ARM_ReferenceValue*>& AmortRefValues)
{

	Init();

	SetAsOfDate(asOf) ;
 
	Set(asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label,issummitcurve,
					VolCurve,PayFreq, calibAlgo,calibData, LagStartDate,params,UpFronts,AmortRefValues);
		
}

void ICM_Constant_Piecewise::Set (const ARM_Date& asOf,
								  const std::vector<std::string>&terms,
                       ARM_Vector* rates,
					   double Recovery,
					   ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
					   qCDS_ADJ adj,  	
                       const std::string&  ccy,
					   const std::string&  label,
					   bool issummitcurve,
					   const ARM_VolCurve* VolCurve,
					   long PayFreq,
					   qDEFCURVE_CALIB_ALGO calibAlgo,
					   const std::string& calibData,
					   int LagStartDate,
					   const ICM_Parameters&params,
					   ARM_Vector* UpFronts,
					   vector<ARM_ReferenceValue*>& AmortRefValues)
{
	if (true)
		ICM_DefaultCurve::Set (asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label,issummitcurve,VolCurve,PayFreq, calibAlgo,calibData, LagStartDate,params,UpFronts,AmortRefValues);
	else
		ICM_DefaultCurve::Set (asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,label,issummitcurve,VolCurve,PayFreq, calibAlgo,calibData, LagStartDate,params,UpFronts,AmortRefValues);
}

ICM_Constant_Piecewise::ICM_Constant_Piecewise(const ARM_Date& asOf,
						 const ARM_Vector& YForDates,
						 const ARM_Vector& rates,
						  double Recovery,
						 ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
						qCDS_ADJ adj ,
						const std::string&  ccy,
						 const std::string&  label,
						 bool issummitcurve,
						 const ARM_VolCurve* VolCurve,
						 const bool& isdatesininput,
						 long PayFreq,
						 qDEFCURVE_CALIB_ALGO	calibAlgo,
						 const std::string& calibData,
						 int LagStartDate)
{

	Init();
 
	Set(asOf,YForDates,rates,Recovery,zc,intRule,adjStartDate,
		adj,ccy,label,issummitcurve,VolCurve,isdatesininput, PayFreq, 
		calibAlgo,calibData, LagStartDate);
		
}

void ICM_Constant_Piecewise::Set (const ARM_Date& asOf,
						 const ARM_Vector& YForDates,
						 const ARM_Vector& rates,
						  double Recovery,
						 ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,	
						qCDS_ADJ adj ,
						const std::string&  ccy,
						 const std::string&  label,
						 bool issummitcurve,
						 const ARM_VolCurve* VolCurve,
						 const bool& isdatesininput,
						 long PayFreq,
						 qDEFCURVE_CALIB_ALGO	calibAlgo,
						 const std::string& calibData,
						 int LagStartDate)
{
	if (true)
		ICM_DefaultCurve::Set (asOf,YForDates,rates,Recovery,zc,intRule,
		adjStartDate,adj,ccy,label,issummitcurve,VolCurve,isdatesininput,PayFreq, calibAlgo,calibData, LagStartDate);
	else
		ICM_DefaultCurve::Set (asOf,YForDates,rates,Recovery,zc,intRule,
			adjStartDate,adj, ccy,label,issummitcurve,VolCurve,isdatesininput,PayFreq, calibAlgo,calibData, LagStartDate);
}


void ICM_Constant_Piecewise::Init(void)
{
	SetName(ICM_CST_PIECEWISE);
}

double ICM_Constant_Piecewise::SurvivalFunction(const double& yearterm) const 
{
	if (GetIsNameInDefault()) return 0.;
	bool equal = false;
	int i = 0;
	// Inf_EqualInterpol(yearterm,i,equal);
	i=position(yearterm); 
	if (eq(itsInterpolYF.Elt(i),yearterm))  
	// if (equal) 
		return GetInterpolSP().Elt(i);
	double SP_prev = GetInterpolSP().Elt(i);
	double Lambda = GetInterpolLambda().Elt(i);
	double Sproba = SP_prev*exp(-Lambda*(yearterm - GetInterpolYF().Elt(i)));	
	return (Sproba);
}

void ICM_Constant_Piecewise::CptTermsSurvivalProba(void)
{

	int size = GetInterpolDates().GetSize();
	// GetSurvivalProba()->Elt(0)=1.;
	SetInterpolSP(0,1); 

	double SP_prev = 0.;
	double Lambda  = 0.;
	double Sproba  = 0.;
	double YT	   = 0.;
	double YT_prev = 0.;
	
	for (int i = 1; i<size; i++)
	{
	 /** SP_prev = GetSurvivalProba(i-1);
	 Lambda  = GetLambda(i-1);
	 YT_prev = GetYearTerm(i-1); **/ 
	SP_prev =GetInterpolSP().Elt(i-1); 
	Lambda   =GetInterpolLambda().Elt(i-1); 
	YT_prev= GetInterpolYF().Elt(i-1); 

	 // YT = GetYearTerm(i);
	YT = GetInterpolYF().Elt(i);

	 Sproba = SP_prev*exp(-Lambda*(YT - YT_prev));
	 
	 // GetSurvivalProba()->Elt(i)=Sproba;
	 SetInterpolSP(i,Sproba);
	}
}


void ICM_Constant_Piecewise::ResetLambda(const int& i,const double& lambda)
{
	int size = GetInterpolDates().GetSize();
	// GetLambdas()->Elt(i-1)=lambda;			
	SetInterpolLambda(i-1,lambda); 

	// if (i == (size-1)) GetLambdas()->Elt(i)=lambda;		
	if (i == (size-1)) SetInterpolLambda(i,lambda);		

	//Reset de la survival proba du plot actif
	double SP_prev = GetInterpolSP().Elt(i-1);
	double YT = (GetInterpolYF().Elt(i) - GetInterpolYF().Elt(i-1));
	double SP = SP_prev * exp(-lambda *YT);
	SetInterpolSP(i,SP);
	
}


ICM_DefaultCurve* ICM_Constant_Piecewise::GenerateShiftCurve(const vector<string>&  pTerm, 
																   const ARM_Vector& epsilon ,
																   qSENSITIVITY_TYPE mode ) const 
{

   ICM_DefaultCurve* dc = NULL;
   //bool check = false;	
   // if (! epsilon) ICMTHROW(ERR_INVALID_ARGUMENT,"Unable to bump curve : no epsilon  \n"); 
   int szinp = epsilon.GetSize();
   if (szinp <0) ICMTHROW(ERR_INVALID_ARGUMENT,"Unable to bump curve : no epsilon  \n"); 
			   
   int size = GetRates()->GetSize()-1;	


   ARM_Vector ModifiedRates(size,0.);
   memcpy(ModifiedRates.GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* size);
   double recovery = GetRecovery();

   int k=0,l=0;
 
   if(GetTerms().empty())
	{
		   ARM_Vector dates(GetLambdas()->GetSize()-1); 				   
		   for (int i=0;i<GetLambdas()->GetSize()-1;i++) dates.Elt(i)=GetDates()->Elt(i+1); 
		   switch (mode)
		   {
		   case ICM_DTR_TYPE:
			   {
				   ARM_Vector lambdas(GetLambdas()->GetSize()-1,1000.); 
				   
				   if (epsilon.Elt(0) != CREDIT_DEFAULT_VALUE) recovery = epsilon.Elt(0);
					dc= new ICM_Constant_Piecewise( GetAsOfDate() ,
												dates ,
												lambdas ,
												0 , // this is not -1 
												recovery ,
												GetZeroCurve() ,  
												GetSTDCDS_Adjusted(),
												GetSTDCDS_AdjustedStartDateOnBusinessDay(),
												GetCdsAdj(),
												GetCurrency()  ,
												GetLabel() , 
												GetVolCurve(),
												GetCalibrationAlgo(),
												GetCalibrationData(),
												this->GetLagStartDate(),
												this->GetParameters()
												) ;

					break;
			   }
			case ICM_GREEK_THETA_TYPE :
			   {
				   ARM_Date asOf = GetAsOfDate();
				   ARM_Date  NewAsOf = asOf.NextBusinessDay((int)epsilon.Elt(0)) ;
				   std::auto_ptr<ARM_ZeroCurve> zc(dynamic_cast<ARM_ZeroCurve*>(GetZeroCurve()->Clone()));
				   zc->SetAsOfDate(NewAsOf);

				   ARM_Vector lambdas(GetLambdas()->GetSize()-1,0.); 
				   for (int i=0;i<GetLambdas()->GetSize()-1;i++) lambdas.Elt(i)=GetLambdas()->Elt(i); 
				   dc= new ICM_Constant_Piecewise( NewAsOf ,
													dates ,
													lambdas ,
													0 , // this is not -1 
													recovery ,
													zc.get() ,  
													GetSTDCDS_Adjusted(),
													GetSTDCDS_AdjustedStartDateOnBusinessDay(),
													GetCdsAdj(),
													GetCurrency()  ,
													GetLabel() , 
													// NULL,
													GetVolCurve(),
												GetCalibrationAlgo(),
												GetCalibrationData(),
													GetLagStartDate(),
													GetParameters()
													) ;
			   break;
			   }
		   default:
			   {
				   ICMTHROW(ERR_INVALID_ARGUMENT,"Unable to bump curve : Invalid constructor type "); 
			   }
		   }// end switch
	   } 
		else	
	   {
		   switch (mode)
		   {
		   case ICM_GREEK_THETA_TYPE : {
			 ARM_Date NewAsOf = (ARM_Date(GetAsOfDate())).NextBusinessDay((int)epsilon.Elt(0)) ;
			 std::auto_ptr<ARM_ZeroCurve> zc(dynamic_cast<ARM_ZeroCurve*>(GetZeroCurve()->Clone()));
			 zc->SetAsOfDate(NewAsOf);
			 dc = new ICM_Constant_Piecewise(NewAsOf,
									GetTerms(),
								   &ModifiedRates,
								   recovery,
								   zc.get(),
								   GetSTDCDS_Adjusted(),
								   GetSTDCDS_AdjustedStartDateOnBusinessDay(),
								   GetCdsAdj(),
								   GetCurrency(),
								   GetLabel(),
								   GetIsSummitCurve(),
								   //2 NULL,
								   GetVolCurve(),
								   GetPayFrequency(),
												GetCalibrationAlgo(),
												GetCalibrationData(),
								   GetLagStartDate(),
								   GetParameters()
								   );

			return(dc);
			}
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
					   if (mode==ICMSPREAD_TYPE) for(k=0;k<size;k++) ModifiedRates.Elt(k) = ModifiedRates.Elt(k) + epsilon.Elt(0)/100.;
					   else for(k=0;k<size;k++) ModifiedRates.Elt(k) = ModifiedRates.Elt(k) *  (1.0 + epsilon.Elt(0));
				   }
				   else 
				   {
					   for (l=0; l<szinp; l++)
						for (k=0; k<size; k++)
						{
						 if ( pTerm[l] == GetTerm(k) )
						 {
							if (mode == ICMSPREAD_TYPE)
								ModifiedRates.Elt(k) = ModifiedRates.Elt(k) + epsilon.Elt(l)/100.;
							else if (mode == ICMSPRELSHIFT_TYPE)
								ModifiedRates.Elt(k) = ModifiedRates.Elt(k) *  (1.0 + epsilon.Elt(l));

							//check = true;
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
				   for (k=0; k<size; k++) ModifiedRates.Elt(k) = 20.;
				   break;
			   }
		   } // end Switch

		   dc = new ICM_Constant_Piecewise(GetAsOfDate(),
			   GetTerms(),
			   &ModifiedRates,
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
			   GetCalibrationAlgo(),
			   GetCalibrationData(),
			   GetLagStartDate(),
			   GetParameters()
			   );  
	} // end Else
   return(dc);
}

/** JLA Removed
ICM_DefaultCurve* ICM_Constant_Piecewise::xGenerateShiftCurve(double epsilon, 
															 qSENSITIVITY_TYPE mode) const 
{
   ICM_DefaultCurve* dc = NULL;
   //bool check = false;	

   int size = GetRates()->GetSize()-1;	

   ARM_Vector* ModifiedRates = new ARM_Vector(size,0.);
   memcpy(ModifiedRates->GetElt(),GetRates()->GetElt()+1 ,sizeof(double)* size);
   double recovery = GetRecovery();
   // Here we really want to use the GetYearTerms . 
   ARM_Vector* ModifiedYearsTerms = new ARM_Vector((ARM_Vector*)GetYearTerms(), size, 1, size, 0);


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

  // if( ! GetpTerms()) {
  if( GetTerms().empty()) {
	  if ( GetRate(0) == CREDIT_DEFAULT_VALUE) {
			ICMTHROW(ERR_SQUARE_OR_SIZE_PB,"Unable to bump curve : Invalid constructor"); 
		} else {
			dc = (ICM_DefaultCurve*) new ICM_Constant_Piecewise(GetAsOfDate(),
														*ModifiedYearsTerms,
														*ModifiedRates,
														recovery,
														GetZeroCurve(),
														GetSTDCDS_Adjusted(),
														GetSTDCDS_AdjustedStartDateOnBusinessDay(),
														GetCdsAdj(),
														GetCurrency(),
														GetLabel(),
														GetIsSummitCurve(),
														GetVolCurve(),
														false,
														//2 NULL,
														GetPayFrequency(),
												GetCalibrationAlgo(),
												GetCalibrationData(),
														GetLagStartDate()
														);	
		}	
	} else {
		dc = (ICM_DefaultCurve*) new ICM_Constant_Piecewise(GetAsOfDate(),
														GetTerms(),
														ModifiedRates,
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
												GetCalibrationAlgo(),
												GetCalibrationData(),
														GetLagStartDate()
														);
   }

   if (ModifiedRates)
		delete ModifiedRates;
   ModifiedRates = NULL;
	if (ModifiedYearsTerms)
		delete ModifiedYearsTerms;
	ModifiedYearsTerms = NULL;
   return(dc);
}

**/ 

ICM_DefaultCurve* ICM_Constant_Piecewise::GenDefCurve(const ARM_Vector& dates,
											  const ARM_Vector& spreads,
											  const double& recovery,
											  const string& label,
											  const bool& isdates,
											  ARM_ZeroCurve* ircurve) const 
{

	ARM_ZeroCurve* IRcurve = ((ircurve) ? ircurve : GetZeroCurve());

	ICM_DefaultCurve* DefCurve = new ICM_Constant_Piecewise(GetAsOfDate(),
												 dates,
												 spreads,
												 recovery,
												 IRcurve,
												 GetSTDCDS_Adjusted(),
												 GetSTDCDS_AdjustedStartDateOnBusinessDay(),
												 // GetCdsAdj(),
												 GetCdsAdj(),
												 // qCredit_Default, //isdates?qCredit_Special_None_Date:qCredit_Special_None_YF, useless
												 GetCurrency(),
												 label.c_str(),
												 GetIsSummitCurve(),
												 GetVolCurve(),
												 isdates,
												//2 NULL,
												K_QUARTERLY,
												GetCalibrationAlgo(),
												GetCalibrationData(),
												GetLagStartDate()
												 );
	return (DefCurve);
}

//**************************************************************************************************************/
//
//	Coompute Hazard rate term structure given the survival probabilites and dates.
//
//**************************************************************************************************************/

ARM_Vector* ICM_Constant_Piecewise::CptLambdas(const ARM_Vector& InputSP, const ARM_Vector& InputDates) const
{
	int size = InputSP.GetSize();
	int i = 0 ;
	ARM_Vector* Lambdas = new ARM_Vector(size,0.);
	ARM_Vector InputYF(size,0.);

	for (i = 0 ; i < size ; i++)
	{ 
		InputYF.Elt(i) = ( InputDates.Elt(i) - GetAsOfDate().GetJulian() ) /K_YEAR_LEN;
	}
	if( InputSP.Elt(0) != 1.)
	{
		if (eq(InputSP.Elt(0),0.)) Lambdas->Elt(0) = 1000.; 
		else Lambdas->Elt(0) = - log(InputSP.Elt(0))/InputYF.Elt(0);

		for (i = 1 ; i < size ; i++)
		{ 
			if (eq(InputSP.Elt(i),0.)) Lambdas->Elt(i) = 1000.; 
			else Lambdas->Elt(i) = - log(InputSP.Elt(i)/InputSP.Elt(i-1))/(InputYF.Elt(i)-InputYF.Elt(i-1)); 
		}
	}
	else
	for (i = 0 ; i < size-1 ; i++)
	{ 
		if (eq(InputSP.Elt(i),0.)) Lambdas->Elt(i) = 1000.; 
		else Lambdas->Elt(i) = - log(InputSP.Elt(i+1)/InputSP.Elt(i))/(InputYF.Elt(i+1)-InputYF.Elt(i)); 
	}	
	return (Lambdas);
}

//**************************************************************************************************************/
//
//	Coompute Hazard rate term structure given the survival probabilites and year fractions.
//
//**************************************************************************************************************/

ARM_Vector* ICM_Constant_Piecewise::CptLambdas(const ARM_Vector& InputSP, double* InputYF) const 
{
	int size = InputSP.GetSize();
	ARM_Vector* Lambdas = new ARM_Vector(size,0.);
	int i = 0 ;
	
	if(/** i==0 && **/ InputSP.Elt(0) != 1. && InputYF[0] != 0. )
	{
		if (eq(InputSP.Elt(0),0.)) Lambdas->Elt(0) = 1000.; 
		else Lambdas->Elt(0) = - log(InputSP.Elt(0))/InputYF[0];; 
	
		for (i = 1 ; i < size ; i++)
		{ 
			if (eq(InputSP.Elt(i),0.)) Lambdas->Elt(i) = 1000.; 
			else Lambdas->Elt(i) = - log(InputSP.Elt(i)/InputSP.Elt(i-1))/(InputYF[i]-InputYF[i-1]) ;
		}
	}
	else // Sp(0)==1 OR YF(0)=0
	for (i = 0 ; i < size-1 ; i++)
	{ 
		if (eq(InputSP.Elt(i),0.)) Lambdas->Elt(i) = 1000.; 
		else Lambdas->Elt(i) = - log(InputSP.Elt(i+1)/InputSP.Elt(i))/(InputYF[i+1]-InputYF[i]);
	}		
	return (Lambdas);
}



double ICM_Constant_Piecewise::DefProbInverse(const double& DefaultProba) const 
{
	if ((lt(DefaultProba,0.)) || (leq(1.,DefaultProba)) )
		ICMTHROW(ERR_INVALID_DATA, "Default Proba is " << DefaultProba << " and must lie between 0 and 1 ");

	double result = 0;

	const ARM_Vector&Yt = GetInterpolYF();
	const ARM_Vector& Sp = GetInterpolSP();
	int size = Yt.GetSize();


	// ! Les probabilités de defaut sont supposées croissantes
	for (int i=0; i<Sp.GetSize()-1; i++) 
	{
		if ( eq(DefaultProba, (1.-Sp.Elt(i))) )
			return Yt.Elt(i);
		else if ( (lt((1.-Sp.Elt(i)),DefaultProba)) && (lt(DefaultProba , (1.-Sp.Elt(i+1)))) )
		{   
			// result = log(Sp.Elt(i)/(1.-DefaultProba))/GetLambda(i);
			result = log(Sp.Elt(i)/(1.-DefaultProba))/GetInterpolLambda().Elt(i); // GetLambda(i);
			result = Yt.Elt(i) + result;
			return result;
		}
	}

	return  (Yt.Elt(size-1) + log(Sp.Elt(size-1)/(1.-DefaultProba))/GetInterpolLambda().Elt(size-1)); 		
}

// virtual 
// void ICM_Constant_Piecewise::Copy(const ARM_Object* srcZc)
// {
// 	ICM_DefaultCurve::Copy(srcZc);
// 	// BitwiseCopy(srcZc);
// }

// virtual 
ARM_Object* ICM_Constant_Piecewise::Clone(void)
{
//     ICM_Constant_Piecewise* theClone = new ICM_Constant_Piecewise();
// 	theClone->Copy(this);
// 	return(theClone);
	return new ICM_Constant_Piecewise(*this); 
}

