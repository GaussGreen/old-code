#include "ARMKernel\glob\firsttoinc.h"

#include "ICMKernel\inst\icm_security.h"
#include "ICMKernel\str\icm_util_str.h"
#include "ICMKernel\util\icm_utils.h"


void ICM_Security::Init()
{
	//Accrual basis
	itsAccrualBasis = KACTUAL_360;

	//Flag YearFraction Date AsOf (if already compute)
	// itsFlgYF = false;
	itsLastAsOf = -1.;

	itsIncludeMaturity = false;
	itsAccruedOnDefault = qACCRUED_SETTLED;


	itsSecurityType = qRUNNING;
	itsFwdsType.resize(0);

	itsPaymentFreq = K_QUARTERLY;
	itsPastFixings.clear(); 
	itsFrozenMaturity = CREDIT_DEFAULT_VALUE;
}


void ICM_Security::Set(const ARM_Date& StartDateNA,
					const ARM_Date& EndDateNA,
					const int& nbflows,
					const ARM_ReferenceValue& premiumNotionals,
					const int& AccrualBasis,
					const std::string& Ccy)
{
	itsStartDateNA = StartDateNA;	
	itsEndDateNA = EndDateNA;		
	

	SetNumFlows(nbflows);
	ARM_Security::SetMaturity(EndDateNA); // unadj? 

	itsAccrualBasis = AccrualBasis;
	
	SetCcy(Ccy); 
	SetAmount(&unconst(premiumNotionals),100); 

// 	itsFlgYF = false;
	itsPastFixings.clear(); 
}
 
void ICM_Security::Set(const ARM_Date& StartDate,
					const ARM_Date& EndDate,
					ARM_Vector* AccStartDates,
					ARM_Vector* AccEndDates,
					ARM_Vector* InterestDays,
					ARM_Vector* PayDates,
					ARM_Vector* CouponRates,
					ARM_Vector* Notionals,
					ARM_Vector* NotionalXchange,
					bool includeMaturity,
					const double& InitialNotional,
					const int& AccrualBasis,
					const std::string& Ccy)
{
	int size = AccStartDates->GetSize();

	
	itsStartDateNA = StartDate;
	itsEndDateNA = EndDate;
	

	ARM_Security::SetMaturity((ARM_Date)AccEndDates->Elt(size-1)); // adj. 
	

	itsAccStartDates =*AccStartDates;
	itsAccEndDates =*AccEndDates;
	itsInterestDays= *InterestDays;
	itsYFInterestDays.Resize(0); 
	itsYFAccStartDates.Resize(0);
	itsYFAccEndDates.Resize(0);
	itsPayDates = *PayDates;
	itsYFPayDates.Resize(0);
	itsCouponRates = *CouponRates;
	itsFwdSpreads.Resize(itsCouponRates.GetSize());
	itsAdjFwdSpreads.Resize(itsCouponRates.GetSize());
	itsUnAdjFwdSpreads.Resize(itsCouponRates.GetSize());
	itsFwdPV01s.Resize(itsCouponRates.GetSize());
	itsPartRates.Resize(itsCouponRates.GetSize());
	for (int i=0;i<itsPartRates.GetSize();i++) {itsPartRates.Elt(i)=1.;}
	itsNotionals = *Notionals;
	itsNotionalXchange = *NotionalXchange;
	itsAccrualBasis = AccrualBasis;
	itsStartRiskDates = *AccStartDates;
	itsEndRiskDates = *AccEndDates;

	SetCcy(Ccy) ;

//	itsFlgYF = false;
	itsPeriodDefLegPV.Resize(itsAccStartDates.GetSize());
	itsPeriodFeeLegPV.Resize(itsAccStartDates.GetSize());
	itsPeriodRefDefLegPV.Resize(itsAccStartDates.GetSize());
	itsPeriodRecDefLegPV.Resize(itsAccStartDates.GetSize());
	itsPastFixings.clear(); 
	itsIncludeMaturity=includeMaturity ;
}


void ICM_Security::BitwiseCopy(const ARM_Object* src)
{
    ICM_Security* gen = (ICM_Security*) src;

	itsStartDateNA = gen->itsStartDateNA;
	itsEndDateNA= gen->itsEndDateNA;

	itsAccStartDates = gen->itsAccStartDates;
	itsAccEndDates = gen->itsAccEndDates;
	itsInterestDays = gen->itsInterestDays ;
	itsYFInterestDays=gen->itsYFInterestDays ;
	itsYFAccStartDates = gen->itsYFAccStartDates;
	itsYFAccEndDates = gen->itsYFAccEndDates;
	itsPayDates = gen->itsPayDates;
	itsYFPayDates = gen->itsYFPayDates;
	itsCouponRates = gen->itsCouponRates;
	itsFwdSpreads = gen->itsFwdSpreads;
	itsAdjFwdSpreads = gen->itsAdjFwdSpreads ;
	itsUnAdjFwdSpreads = gen->itsUnAdjFwdSpreads ;
	itsFwdPV01s = gen->itsFwdPV01s ;
	itsNotionals = gen->itsNotionals;
	itsNotionalXchange = gen->itsNotionalXchange;

	itsAccrualBasis = gen->itsAccrualBasis;
	itsStartRiskDates = gen->itsStartRiskDates;
	itsYFStartRiskDates = gen->itsYFStartRiskDates;
	itsEndRiskDates = gen->itsEndRiskDates;
	itsYFEndRiskDates = gen->itsYFEndRiskDates;

	SetCcy(gen->GetCcy()); 

	itsPeriodDefLegPV = gen->itsPeriodDefLegPV;
	itsPeriodFeeLegPV = gen->itsPeriodFeeLegPV;
	itsPeriodRefDefLegPV = gen->itsPeriodRefDefLegPV;
	itsPeriodRecDefLegPV = gen->itsPeriodRecDefLegPV;

	itsLastAsOf = gen->itsLastAsOf;

	itsIncludeMaturity = gen->itsIncludeMaturity;

	
	itsAccruedOnDefault = gen->itsAccruedOnDefault;
	itsIncludeMaturity = gen->itsIncludeMaturity;
	itsSecurityType = gen->itsSecurityType;

	itsFwdStartDates = gen->itsFwdStartDates;
	itsFwdEndDates = gen->itsFwdEndDates;
	itsResetDates = gen->itsResetDates;
	itsYFFwdStartDates = gen->itsYFFwdStartDates  ;
	itsYFFwdEndDates = gen->itsYFFwdEndDates  ;
	itsYFResetDates = gen->itsYFResetDates;

	itsFwdsType = gen->itsFwdsType;

	itsPaymentFreq = gen->itsPaymentFreq;
	itsPastFixings = gen->itsPastFixings ; 

	itsFrozenMaturity = gen->itsFrozenMaturity;
	itsPartRates = gen->itsPartRates;
}

void ICM_Security::Copy(const ARM_Object* src)
{
	ARM_Security::Copy(src);
    BitwiseCopy(src);
}

ARM_Object* ICM_Security::Clone(void)
{
	 ICM_Security* theClone = new ICM_Security();

     theClone->Copy(this);
 
     return(theClone);
}

// JLA . 
ARM_Object* ICM_Security::Clone() const 
{
	return unconst(this)->Clone(); 
}


// *************************************************************
// View Matrix 
// *************************************************************
void ICM_Security::View(char* id, FILE* ficOut)
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

		if (itsAccStartDates.GetSize()>0)
		{

		// Affichage de la matrice
		fprintf(fOut, "\t\t\t ----------------- Periods Informations ----------------- \n\n");

		if (IsFrozenMaturity()) 
		{ 
			char frozend[20];
			GetFrozenMaturity().JulianToStrDate(frozend);
			fprintf(fOut,"Frozen Maturity:%14s\n",frozend);
		}

		char d0[20]; 
		itsStartDateNA.JulianToStrDate(d0) ;fprintf(fOut,"StartDateNA: %s\n",d0); 
		itsEndDateNA.JulianToStrDate(d0) ;fprintf(fOut,"EndDateNA: %s\n",d0); 
		int size = itsAccStartDates.GetSize();
	    char d1[20];
	    char d2[20];
	    char d3[20];
	    char d4[20];
	    char d5[20];
		char d6[20];
		char d7[20];
		char d8[20];
		char d9[20];
		char d10[20];

	    int  i;
		double notional = 0., cfValue = 0.;

	    fprintf(fOut, "%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t\n",
         "Acc Start Dates",
		 "Acc End Dates", 
		 "Payment Dates",
		 "NbInterestDays",
		 "DefaultProba",
		 "DiscountRate",
		 "CpnRates(%)",
		 "FwdTypes",
		 "FwdRates(%)",
		 "PartRate",
		 "Notionals",
		 "NotionalsXChg",
		 "RefDefLegPV",
		 "RecDefLegPV",
		 "DefLegPV",
		 "FeeLegPV",
		 "PeriodPV",
		 "Reset Dates",
		 "Fwd StartDates", 
		 "Fwd End Dates",
		 "Risk StartDates",
		 "Risk EndDates");

/*	    fprintf(fOut, "%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\n\n",
         "Acc Start Dates","Acc End Dates", "Payment Dates","Reset Dates","Fwd StartDates", "Fwd End Dates","Risk StartDates","Risk EndDates",
		 "NbInterestDays","CpnRates(bp)","FwdTypes","FwdRates(bp)","PartRate","Notionals","NotionalsXChg","RefDefLegPV","RecDefLegPV","DefLegPV","FeeLegPV","PeriodPV");
*/
			for (i = 0; i < size; i++)
			{
			((ARM_Date) (itsAccStartDates)[i]).JulianToStrDate(d1);
			((ARM_Date) (itsAccEndDates)[i]).JulianToStrDate(d2);
			((ARM_Date) (itsPayDates)[i]).JulianToStrDate(d3);
			if (itsResetDates.GetSize()>0) ((ARM_Date) (itsResetDates)[i]).JulianToStrDate(d6); else ((ARM_Date) (itsAccStartDates)[i]).JulianToStrDate(d6);
			if (itsFwdStartDates.GetSize()>0) ((ARM_Date) (itsFwdStartDates)[i]).JulianToStrDate(d4); else ((ARM_Date) (itsAccStartDates)[i]).JulianToStrDate(d4);
			if (itsFwdEndDates.GetSize()>0) ((ARM_Date) (itsFwdEndDates)[i]).JulianToStrDate(d5); else ((ARM_Date) (itsAccStartDates)[i]).JulianToStrDate(d5);
			if (itsStartRiskDates.GetSize()>0) ((ARM_Date) (itsStartRiskDates)[i]).JulianToStrDate(d7); else ((ARM_Date) (itsAccStartDates)[i]).JulianToStrDate(d7);
			if (itsEndRiskDates.GetSize()>0) ((ARM_Date) (itsEndRiskDates)[i]).JulianToStrDate(d8); else ((ARM_Date) (itsAccEndDates)[i]).JulianToStrDate(d8);
			
			if (itsDefaultProbability.GetSize()>0) {sprintf(d9,"%10.4lf",itsDefaultProbability.Elt(i));} else {memset(d9,'\0',20);}
			if (itsDiscountRate.GetSize()>0) {sprintf(d10,"%10.4lf",itsDiscountRate.Elt(i));} else {memset(d10,'\0',20);}

			string FwdType = "Fixed";
			if (itsFwdsType.size()>0) FwdType = ConvertFwdType(itsFwdsType[i]);

			fprintf(fOut,"%14s\t%14s\t%14s\t%10.4lf\t%14s\t%14s\t%10.4lf\t%14s\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%14s\t%14s\t%14s\t%14s\t%14s\t\n", 
				d1, 
				d2, 
				d3,
				itsInterestDays.Elt(i),
				d9,
				d10,
				itsCouponRates.Elt(i),
				FwdType.c_str(),
				itsFwdSpreads.Elt(i),
				itsPartRates.Elt(i),
				itsNotionals.Elt(i),
				itsNotionalXchange.Elt(i),
				itsPeriodRefDefLegPV.Elt(i),
				itsPeriodRecDefLegPV.Elt(i),
				itsPeriodDefLegPV.Elt(i),
				itsPeriodFeeLegPV.Elt(i),
				itsPeriodDefLegPV.Elt(i) + itsPeriodFeeLegPV.Elt(i),d6, d4, d5,d7,d8);

			}
		}

		//Pastfixings --------------------------------------------------------------
		fprintf(fOut,"\n\n");
		fprintf(fOut,"PastFixings size: %d\n",itsPastFixings.size())  ;
		if (itsPastFixings.size()!=0) 
		{
			std::map<double,double>::const_iterator it=itsPastFixings.begin() ;
			while (it!=itsPastFixings.end()) 
			{
				char d[20]; 
				ARM_Date(it->first).JulianToStrDate(d); 
				fprintf(fOut,"%s : %10.4lf\n",d,it->second) ;
				++it; 
			}
		}
		fprintf(fOut, "\n");

		//informations -------------------------------------------------------------

		if ( ficOut == NULL )
		{
		fclose(fOut);
		}
	}


// -------------------------------------------------------------------------------------------
//	Retrieves the current period index such that index = max{i in N/ startdate[i]<= AsofDate}
// -------------------------------------------------------------------------------------------
int ICM_Security::PeriodIndex(const ARM_Date& AsofDate)   
{
	
	int index = 0;
	const ARM_Vector& startDates = GetAccStartDates();
	int numflows = startDates.size();
	const ARM_Vector& endDates = GetAccEndDates();

	if (AsofDate.GetJulian() > endDates.Elt(numflows-1))
		return (-1);

	while (startDates.Elt(index) <= AsofDate.GetJulian()) 
	{ index++;
      if (index >= numflows) break;	}

	index--;
	return (index);
}

// -------------------------------------------------------------------------------------------
//	Is in last period
// -------------------------------------------------------------------------------------------
bool ICM_Security::IsLastPeriod(const ARM_Date& AsofDate)   
{
	
	int index = 0;
	const ARM_Vector& startDates = GetAccStartDates();
	const ARM_Vector& endDates = GetAccEndDates();
	int numflows = startDates.size();
	if ((AsofDate.GetJulian() >= startDates.Elt(numflows-1)) && (AsofDate.GetJulian() <= endDates.Elt(numflows-1)))
		return true;

	return (false);
}

// -------------------------------------------------------------------------------------------
//	calculates the daily accrued coupon for a given period (i)
//
//		
// -------------------------------------------------------------------------------------------
double ICM_Security::DailyAccrued(int& index,bool includenot)
{
	double result = 0.,accdc=0.,coupon=0.;
	if (includenot) coupon = FullCoupon(index);	
	else coupon = FullCouponRate(index);	

	accdc= itsInterestDays.Elt(index); 
	result = (coupon/accdc);

	return (result);
}


// ------------------------------------------------------------------------------------
//	Calculates the differential accrued between two yf for t1 & t2 in a same period
// ------------------------------------------------------------------------------------
double ICM_Security::DiffAccrued  (const ARM_Date& AsOf,const double& yf1,const double& yf2,const bool& includenot)
{
	ARM_Date date1 = (ARM_Date)(AsOf.GetJulian()+round(365.*yf1));
	ARM_Date date2 = (ARM_Date)(AsOf.GetJulian()+round(365.*yf2));

	return DiffAccrued(date1,date2,includenot);
}

// ------------------------------------------------------------------------------------
//	Calculates the differential accrued between two dates for t1 & t2 in a same period
// ------------------------------------------------------------------------------------
double ICM_Security::DiffAccrued  (const ARM_Date& t1, const ARM_Date& t2,const bool& includenot)
{
	int index1 = PeriodIndex(t1);
	double result =0.;

	if ((index1>=0) && (t2.GetJulian() <=  GetMaturity().GetJulian()))  // GetMaturity().GetJulian()))
		result = DailyAccrued(index1,includenot) * DaysBetweenDates(itsAccrualBasis,t1,t2);
	else
		result = 0.;

	return (result);
}
 
// ------------------------------------------------------------------------------------
//	Calculates the coupon rate for a given index
// ------------------------------------------------------------------------------------
double ICM_Security::FullCoupon(const int& index) 
{ 
	if (index<0) return 0.;
	double result = itsNotionals.Elt(index)* FullCouponRate(index);
	return result;
}


// ------------------------------------------------------------------------------------
//	Calculates the full coupon for a given period
// ------------------------------------------------------------------------------------
double ICM_Security::FullCouponRate(const int& index)
{
	ARM_Date FwdFixedDate;
	double result = 0.,rate=0.;

	//Cas du ZeroCoupon, on retourne spread*notional
	if (GetPaymentFreq()==0) 
	{
		result = itsYFInterestDays.Elt(0); 
		rate = itsCouponRates.Elt(0);

		if (itsFwdsType.size()>0)
		switch (itsFwdsType[0])
		{
		case (qCPN_FIXED_AND_FWD) :
		rate+=itsFwdSpreads.Elt(0)*itsPartRates.Elt(0);
		break;
		case (qCPN_FWD) :
		rate=itsFwdSpreads.Elt(0)*itsPartRates.Elt(0);
		break;
		case (qCPN_FIXED) :
		default :;
		};

		return result*rate/100.;
	}

	if (index >= 0)
	{	
		// result = CountYears(itsAccrualBasis,itsAccStartDates.Elt(index),Maturity);
		result = itsYFInterestDays.Elt(index); 
		rate = itsCouponRates.Elt(index);

		if (itsFwdsType.size()>0)
		switch (itsFwdsType[index])
		{
		case (qCPN_FIXED_AND_FWD) :
			rate+=itsFwdSpreads.Elt(index)*itsPartRates.Elt(index);
		break;
		case (qCPN_FWD) :
			rate=itsFwdSpreads.Elt(index)*itsPartRates.Elt(index);
		break;
		case (qCPN_FIXED) :
		default :;
		};
		result *= rate/100.;
	}
	else return 0.;

	return (result);
}

// ------------------------------------------------------------------------------------
//	Calculates the full coupon for a given date
// ------------------------------------------------------------------------------------
double ICM_Security::FullCoupon(const ARM_Date& date)
{
	int index = PeriodIndex(date);
	double result = FullCoupon(index);

	return (result);
}

// ------------------------------------------------------------------------------------
//	Returns the Notional at Date
// ------------------------------------------------------------------------------------
double ICM_Security::Notional(const ARM_Date& date)
{
	double result = 0.;
	int index = PeriodIndex(date);

	if (index >= 0)	result = itsNotionals.Elt(index);
	else return (0.);

	return (result);
}

// ------------------------------------------------------------------------------------
//	Returns the Notional at Date
// ------------------------------------------------------------------------------------
double ICM_Security::NotionalXChange(const ARM_Date& date)
{
	double result = 0.;
	int index = PeriodIndex(date);

	if (index >= 0) result = itsNotionalXchange.Elt(index);
	else return (0.);

	return (result);
}

// ------------------------------------------------------------------------------------
//	Calculates the accrued coupon from the last coupon date before Date until Date
// ------------------------------------------------------------------------------------
double ICM_Security::AccruedCoupon(const ARM_Date& date,const bool& includenot)
{
	ARM_Date AsOf = date;
	int gap = GetSettlementGap();
	double result = 0.;
	ARM_Date fsd;
	int numflows = GetAccStartDates().size();

	if (AsOf >= GetAccEndDates().Elt(numflows-1))
		return 0.;

	int index = PeriodIndex(AsOf);

	if (index >= 0)
		fsd = (GetAccStartDates()).Elt(PeriodIndex(AsOf));
	else
		return 0.;

	result = DiffAccrued(fsd, AsOf,includenot); 
 
	return (result);
}

//----------------------------------------------------------------------
// Computes Year Fractions for all vectors
//----------------------------------------------------------------------
void ICM_Security::ComputeYF(const ARM_Date &AsOf)
{

	if (fabs(itsLastAsOf-AsOf.GetJulian())<1.E-8)
		return;

	itsLastAsOf = AsOf.GetJulian();

	if (!(itsAccStartDates.GetSize()>0) &&
		(itsAccEndDates.GetSize()>0) &&
		(itsInterestDays.GetSize()>0) &&
		(itsPayDates.GetSize()>0))
	{
		    throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "Dates not yet computed");
	}
		
	int size = itsAccStartDates.GetSize();
	itsYFInterestDays.Resize(size); 
	itsYFAccStartDates.Resize(size);
	itsYFAccEndDates.Resize(size);
	itsYFPayDates.Resize(size);
	itsYFFwdStartDates.Resize(size);
    itsYFFwdEndDates.Resize(size);
    itsYFResetDates.Resize(size);
	itsYFStartRiskDates.Resize(size);
	itsYFEndRiskDates.Resize(size);
	
	for (int i=0; i<size; i++)
	{
		if (AsOf.GetJulian() <= itsAccStartDates.Elt(i))
			itsYFAccStartDates.Elt(i)= (itsAccStartDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		if (AsOf.GetJulian() <= itsAccEndDates.Elt(i))
			itsYFAccEndDates.Elt(i) =  (itsAccEndDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		itsYFInterestDays.Elt(i)= CountYears(itsAccrualBasis,itsInterestDays.Elt(i),
			itsAccStartDates.Elt(i),itsAccEndDates.Elt(i));

		if (AsOf.GetJulian() <= itsPayDates.Elt(i))
			itsYFPayDates.Elt(i)= (itsPayDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		if (AsOf.GetJulian() <= itsFwdStartDates.Elt(i))
			itsYFFwdStartDates.Elt(i)= (itsFwdStartDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		//if (FwdEndDates.Elt(0) ==0) { // ie maturityDate in Index

			if (AsOf.GetJulian() <= itsFwdEndDates.Elt(i))
			itsYFFwdEndDates.Elt(i)= (itsFwdEndDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		if (AsOf.GetJulian() <= itsResetDates.Elt(i))
			itsYFResetDates.Elt(i)= (itsResetDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		if (AsOf.GetJulian() <= itsStartRiskDates.Elt(i))
			itsYFStartRiskDates.Elt(i) =  (itsStartRiskDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;

		if (AsOf.GetJulian() <= itsEndRiskDates.Elt(i))
			itsYFEndRiskDates.Elt(i) =  (itsEndRiskDates.Elt(i) - AsOf.GetJulian())/K_YEAR_LEN;
	}
//	itsFlgYF = true;
}
 
//	JLA 
//	this method will count the flows whose paymentDate > asofDate
//	and return the position of the first of these flows. 
//	If there are no flows, then first will be meaningless.
unsigned long 
ICM_Security::countPaymentFlows(unsigned long& first,const ARM_Date&asof) const
{
	const ARM_Vector& payDates = GetPayDates(); 

	first=0; 
	unsigned long i =0; 
	unsigned long size=payDates.size() ;
	for(i=0;i<payDates.size();i++)
	{
		if (payDates.Elt(i)<=asof.GetJulian()) { size--; first++; }
		else return size;	// returns at the first element > asofDate
	}
	return size ;
}

//	------------------------------------------------------------------------------------------
//
//		ForwardSpreads are replaced by PastFixings 
//		when it happens that those ones are known 
// 
void 
ICM_Security::updateFwdSpreadsFromPastFixings()
{
	int i; 
	for (i=0;i<itsResetDates.size();i++) 
	{
		std::map<double,double>::const_iterator it = itsPastFixings.find( itsResetDates.Elt(i) ); 
		if (it!=itsPastFixings.end()) 
		{
			// temporary trace.. 
			//ICMLOG("Replacing fwdSpread #"<<i<<"="<<itsFwdSpreads.Elt(i)<<" with past fixing="<<it->second); 
			itsFwdSpreads.Elt(i)=it->second; 
		}
	}
}
void 
ICM_Security::setPastFixings(const std::map<double,double>&fixings)
{
	itsPastFixings=fixings; 
}
//	------------------------------------------------------------------------------------------
//		This will returns the past fixing at a specific reset Date 
//		is it exists. 
bool 
ICM_Security::getPastFixing(const ARM_Date&resetDate,double&value) const
{
	value=0; 
	std::map<double,double>::const_iterator it = itsPastFixings.find( resetDate.GetJulian()); 
	if (it==itsPastFixings.end()) return false; 
	value=it->second; 
	return true; 
}
//	------------------------------------------------------------------------------------------
//		This is a helper that find the next flow occuring after asofDate 
//		and finds the associated reset Date. Such a reset Date can not exist
bool 
ICM_Security::getLastFixingDate(const ARM_Date&asofDate,ARM_Date&resetDate) const
{
	//	this could happen for non-reset schedules (such as fixed ones ?)
	if (itsResetDates.size()==0)
		return false; 
	const ARM_Vector& payDates = GetPayDates(); 
	unsigned long i=0; 
	while ( i<payDates.size()  && payDates.Elt(i)<=asofDate.GetJulian() ) i++ ;
	// check if we are out of schedule ?
	if (i==payDates.size()) 
		return false; 
	//	i is the required position. check that ResetDates is OK : should not happen. 
	if (i>=itsResetDates.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"getLastFixingDate "<<asofDate.toString()<<" : wrong ResetDates array!"); 
	resetDate=ARM_Date(itsResetDates.Elt(i)); 
	return true ; 
}
//	------------------------------------------------------------------------------------------
int 
ICM_Security::PaymentIndex(const ARM_Date& AsofDate) const
{
	const ARM_Vector& payDates = GetPayDates(); 
	unsigned long i=0; 
	while ( i<payDates.size()  && payDates.Elt(i)<=AsofDate.GetJulian() ) i++ ;
	if (i==payDates.size()) 
	{
		// ICMLOG("ICM_Security::PaymentIndex: no flow after "<<AsofDate.toString()); 
		return -1; 
	}
	return i; 
}

// 
ARM_Date 
ICM_Security::GetMaturity() const 
{
	return ARM_Security::GetMaturity(); 
}