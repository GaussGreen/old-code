#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\inst\icm_corridorleg.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\pricer\icm_Pricer_corridor.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ICMKernel\mod\modelmulticurves.h"


void ICM_CorridorLeg::BitwiseCopy(const ARM_Object* srcleg)
{
    ICM_CorridorLeg* leg = (ICM_CorridorLeg *) srcleg;
	if (leg->itsRefValue_kup)
		itsRefValue_kup = (ARM_ReferenceValue*) leg->itsRefValue_kup->Clone();
	else {
		if (itsRefValue_kup) delete itsRefValue_kup; itsRefValue_kup = NULL;
	}
	if (leg->itsRefValue_kdw)
		itsRefValue_kdw = (ARM_ReferenceValue*) leg->itsRefValue_kdw->Clone();
	else {
		if (itsRefValue_kdw) delete itsRefValue_kdw; itsRefValue_kdw = NULL;
	}
	itsLeverageFloatIdx = leg->itsLeverageFloatIdx;
	itsResetfreq = leg->itsResetfreq; 
	if (leg->itsSpreads) 
		itsSpreads = (ARM_ReferenceValue*) leg->itsSpreads->Clone();
	else {
		if (itsSpreads) delete itsSpreads; itsSpreads = NULL;
	}

	itsSwapStart = leg->itsSwapStart;
	itsSwapEnd = leg->itsSwapEnd;
	itsSubStart = leg->itsSubStart;
	itsSubEnd = leg->itsSubEnd;
	itsSubPayDate = leg->itsSubPayDate;
	itsSubReset = leg->itsSubReset;
	itsSubFwdStart = leg->itsSubFwdStart;
	itsSubFwdEnd = leg->itsSubFwdEnd;

	itsResetProbas = leg->itsResetProbas;
	itsVolK1 = leg->itsVolK1;
	itsVolK2 = leg->itsVolK2;
	itsFwdSpreadAdj = leg->itsFwdSpreadAdj;
	itsResetProbas_Float = leg->itsResetProbas_Float;
	itsFwdSpreadAdj_Float = leg->itsFwdSpreadAdj_Float;

}

void ICM_CorridorLeg::Copy(const ARM_Object* srcleg)
{
     ICM_Leg::Copy(srcleg);
 
     BitwiseCopy(srcleg);
}


ARM_Object* ICM_CorridorLeg::Clone(void)
{
     ICM_CorridorLeg* theClone = new ICM_CorridorLeg();

     theClone->Copy(this);
 
     return(theClone);
}

/*----------------------------------------------------------------------------*
    General Constructor of a Corridor
*----------------------------------------------------------------------------*/ 

ICM_CorridorLeg::ICM_CorridorLeg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
			const ARM_Date* FstCpnEffDate,
			ARM_ReferenceValue* Spreads,
			ARM_IRIndex* floatingIdx,
			const double LeverageFloatIdx,
			ICM_Credit_Index* creditIdx,
			ARM_ReferenceValue* RefValue_kup,
			ARM_ReferenceValue* RefValue_kdw,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED */,
			const int& AccruedDayCount /*= KACTUAL_365*/,
			const int& freq /*= K_ANNUAL*/, 
			const int& Resetfreq /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /*= NULL*/)
	{
		Init();

		ICM_Leg::Set(startDate,endDate,refDate,FstCpnEffDate,0.,AccruedOnDefault,
			AccruedDayCount,0.,K_RCV,freq,dayCount,K_COMP_PROP,
            payTiming,intRule,stubRule,discountCcy,payCalName,K_NX_NONE,
			EXCLUDE_MATURITY,K_ADJUSTED,qCorridor_Leg,CREDIT_DEFAULT_VALUE,ISSUER_UNDEFINE);

		SetName(ICM_CORRIDORLEG);
		SetSpreads(Spreads);
		SetRefValue_kup(RefValue_kup); 
		SetRefValue_kdw(RefValue_kdw);
		SetLeverageFloatIdx(LeverageFloatIdx);
		SetResetfreq(Resetfreq);
		SetCreditIndex((ICM_Credit_Index*)creditIdx->Clone());
		SetIRIndex(floatingIdx);
		GetIRIndex()->SetPayFrequency(freq);
		GetIRIndex()->SetIntRule(intRule);
		GetIRIndex()->SetPayTiming(payTiming);
		GetIRIndex()->SetResetFrequency(freq);
		CptCashFlowDatesCredit();
		GenCorridorSchedule();
	}


/*----------------------------------------------------------------------------*
    General Constructor of a Corridor with Notional
*----------------------------------------------------------------------------*/ 

ICM_CorridorLeg::ICM_CorridorLeg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
			const ARM_Date* FstCpnEffDate,
			ARM_ReferenceValue* Spreads,
			const double Notional,
			const long recieveOrPay,
			ARM_IRIndex* floatingIdx,
			const double LeverageFloatIdx,
			ICM_Credit_Index* creditIdx,
			ARM_ReferenceValue* RefValue_kup,
			ARM_ReferenceValue* RefValue_kdw,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED */,
			const int& AccruedDayCount /*= KACTUAL_365*/,
			const int& freq /*= K_ANNUAL*/, 
			const int& Resetfreq /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /*= NULL*/)
	{
		Init();

		ICM_Leg::Set(startDate,endDate,refDate,FstCpnEffDate,0.,AccruedOnDefault,
			AccruedDayCount,0.,K_RCV,freq,dayCount,K_COMP_PROP,
            payTiming,intRule,stubRule,discountCcy,payCalName,K_NX_NONE,
			EXCLUDE_MATURITY,K_ADJUSTED,qCorridor_Leg,CREDIT_DEFAULT_VALUE,ISSUER_UNDEFINE); 
		
		/*ICM_Leg::Set(startDate,  endDate, 	refDate,	FstCpnEffDate,	0., 	Notional,
			NULL,  NULL,	NULL, freq,   dayCount,     payTiming,  intRule, stubRule,    discountCcy,
            payCalName,	qCorridor_Leg,		EXCLUDE_MATURITY,	K_ADJUSTED,	floatingIdx, CREDIT_DEFAULT_VALUE,
			ISSUER_UNDEFINE,	K_NX_NONE,	AccruedOnDefault); */

		ARM_ReferenceValue refvalfix(Notional , 1 /* price */, 0 /* K_CONSTANT */);
		SetAmount(&refvalfix,100.0);
		
		SetInitRcvOrPay(recieveOrPay);
		SetName(ICM_CORRIDORLEG);
		//SetName(ICM_LEG);
		Spreads->operator /=(100); //from spread in BP but like in icm_leg inn %  : /100
		SetSpreads(Spreads);
		SetRefValue_kup(RefValue_kup); // in BP
		SetRefValue_kdw(RefValue_kdw); // in BP
		SetLeverageFloatIdx(LeverageFloatIdx);
		SetResetfreq(Resetfreq);
		SetCreditIndex((ICM_Credit_Index*)creditIdx->Clone());
		SetIRIndex(floatingIdx);
		GetIRIndex()->SetPayFrequency(freq);
		GetIRIndex()->SetIntRule(intRule);
		GetIRIndex()->SetPayTiming(payTiming);
		GetIRIndex()->SetResetFrequency(freq);
		CptCashFlowDatesCredit();
		GenCorridorSchedule();
	}



/* ****************************************************************************************************************** */
/*!	\fn void View()
/* ****************************************************************************************************************** */
void ICM_CorridorLeg::View(char* id, FILE* ficOut)
{
	FILE* fOut;
    char fOutName[200];
	bool status = true;
   
 
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
	fprintf(fOut, "\t\t\t ----------------- Corridor Leg Credit Index  ----------------- \n");
	GetCreditIndex()->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- Leg   ----------------- \n");
	ICM_Leg::View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- Corridor Leg Informations ----------------- \n");

	fprintf(fOut,"Swap Start\tSwap End\tStart Date\tEnd Date\tReset Date\tFwd Start\tFwd End\t\tProba\tVolK1\tVolK2\tFwdSpreadAdj(bp)\tFwdSpreadAdjFloat(bp)\tProba RangeFloat\t \n");
	
	string name = GetSingleName();
	fprintf(fOut, "\tName %s \n",name.c_str());

	string Sname = GetSingleName();
	fprintf(fOut, "\tName %s \n",Sname.c_str());
	for (int i=0; i<itsSwapStart.size();i++) {
		for (int j=0;j<itsSubReset[i].size();j++)
		{	
			char d1[20];	
			char d2[20];
			char d3[20];
			char d4[20];
			char d5[20];
			char d6[20];
			char d0[20];
			
			if (itsSwapStart.size()>0.) ((ARM_Date) (itsSwapStart)[i]).JulianToStrDate(d5); 
			if (itsSwapEnd.size()>0.) ((ARM_Date) (itsSwapEnd)[i]).JulianToStrDate(d6); 
			if (itsSubStart.size()>0.) ((ARM_Date) (itsSubStart)[i][j]).JulianToStrDate(d3); 
			if (itsSubEnd.size()>0.) ((ARM_Date) (itsSubEnd)[i][j]).JulianToStrDate(d4); 
			if (itsSubReset.size()>0.) ((ARM_Date) (itsSubReset)[i][j]).JulianToStrDate(d0); 
			if (itsSubFwdStart.size()>0.) ((ARM_Date) (itsSubFwdStart)[i][j]).JulianToStrDate(d1);
			if (itsSubFwdEnd.size()>0.) ((ARM_Date) (itsSubFwdEnd)[i][j]).JulianToStrDate(d2);

			if (itsSwapStart.size()>0.) fprintf(fOut,"%14s",d5);else  fprintf(fOut,"\t\t");
			if (itsSwapEnd.size()>0.) fprintf(fOut,"%14s",d6);	else  fprintf(fOut,"\t\t");
			if (itsSubStart.size()>0.) fprintf(fOut,"%14s",d3); else  fprintf(fOut,"\t\t");
			if (itsSubEnd.size()>0.) fprintf(fOut,"%14s",d4);	else  fprintf(fOut,"\t\t");
			if (itsSubReset.size()>0.) fprintf(fOut,"%14s",d0);	else  fprintf(fOut,"\t\t");
			if (itsSubFwdStart.size()>0.) fprintf(fOut,"\t%14s",d1); else fprintf(fOut,"\t\t");
			if (itsSubFwdEnd.size()>0.) fprintf(fOut,"\t%14s\t",d2); else fprintf(fOut,"\t\t");

			if ((itsResetProbas.size()>0. ) && itsResetProbas[i].size() >0 ) fprintf(fOut,"\t%1.4lf",itsResetProbas[i][j]); 
	
			if (itsVolK1.size()>0.  && itsVolK1[i].size() >0) fprintf(fOut,"\t%1.4lf",itsVolK1[i][j]); else fprintf(fOut,"\t");
			if (itsVolK2.size()>0. &&  itsVolK2[i].size() >0 ) fprintf(fOut,"\t%1.4lf",itsVolK2[i][j]); else fprintf(fOut,"\t");
			if (itsFwdSpreadAdj.size()>0. && itsFwdSpreadAdj[i].size() >0) fprintf(fOut,"\t%1.2lf",itsFwdSpreadAdj[i][j]); else fprintf(fOut,"\t");
			if (itsFwdSpreadAdj_Float.size()>0. && itsFwdSpreadAdj_Float[i].size() >0 ) fprintf(fOut,"\t\t\t%1.2lf",itsFwdSpreadAdj_Float[i][j]); else fprintf(fOut,"\t\t\t");
			if (itsResetProbas_Float.size()>0. && itsResetProbas_Float[i].size() > 0) fprintf(fOut,"\t\t%1.4lf",itsResetProbas_Float[i][j]); else fprintf(fOut,"\t\t");	
#ifdef _DEBUG
			long dayOfWeek = (ARM_Date((itsSubReset)[i][j])).GetDayOfWeek();
			if (itsSubReset.size()>0.) fprintf(fOut,"\t : %ld",dayOfWeek);
#endif
			fprintf(fOut,"\n");
		}
	}
	fprintf(fOut,"\n");
	
	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

/* ****************************************************************************************************************** */
/*!	GenCorridorSchedule
/* ****************************************************************************************************************** */

void ICM_CorridorLeg::GenCorridorSchedule()
{
	ICM_Security* security = GetCreditInfos();
	int	NbFlows = security->GetAccStartDates().GetSize();
	int Ssubsize = 0;
	ARM_Vector* SchedSubStart=NULL;
	ARM_Vector* SchedSubEnd=NULL;
	ICM_Credit_Index* idxindex = GetCreditIndex();
	ARM_IRIndex* irindex =  GetIRIndex();
	qCDS_ADJ fwdadj = GetCreditIndex()->GetAdjForTenor();
	double tenor = GetCreditIndex()->GetYearTerm();
	int resetgap = idxindex->GetResetGap();
	int resettiming = idxindex->GetResetTiming();
	int paytiming = idxindex->GetPayTiming();
	int paygap = idxindex->GetResetGap();
	int resetfrq = idxindex->GetResetFrequency();
	int fwdrule = idxindex->GetFwdRule();
	int intrule = idxindex->GetIntRule();

	double paydate=0.,reset=0.,fwdstart=0.,fwdend=0.;

	itsSwapStart.Resize(NbFlows);
	itsSwapEnd.Resize(NbFlows);

	itsSubStart.clear();	
	itsSubEnd.clear();		
	itsSubReset.clear();	itsSubReset.resize(NbFlows);
	itsSubFwdStart.clear(); itsSubFwdStart.resize(NbFlows);
	itsSubFwdEnd.clear();	itsSubFwdEnd.resize(NbFlows);
	itsSubPayDate.clear();	

	if (resetfrq != K_WEEKLY || idxindex->GetCM_resetOccur() != 0 ) {
		itsSubStart.resize(NbFlows);
		itsSubEnd.resize(NbFlows);
		itsSubPayDate.resize(NbFlows);
		for (int i=0;i<NbFlows;i++)
		{
			ARM_Date startdate = security->GetAccStartDates().Elt(i);
			ARM_Date enddate = security->GetAccEndDates().Elt(i);
			ARM_Date _paydate = security->GetPayDates().Elt(i);
			ARM_Date resetDate = security->GetResetDates().Elt(i);

			itsSwapStart.Elt(i) = startdate.GetJulian();
			itsSwapEnd.Elt(i) = enddate.GetJulian();

			SchedSubStart= CptStartDates(startdate,enddate,resetfrq,fwdrule,
								K_LONGSTART,intrule,GetCurrencyUnit()->GetCcyName(),/*adjFirstdate*/1);
			SchedSubEnd= CptEndDates(SchedSubStart,enddate,fwdrule,intrule,GetCurrencyUnit()->GetCcyName());

			int subsize = SchedSubStart->GetSize();
			Ssubsize += subsize;
			ARM_Vector vSubStart(subsize);
			ARM_Vector vSubEnd(subsize);
			ARM_Vector vSubReset(subsize);
			ARM_Vector vSubFwdStart(subsize);
			ARM_Vector vSubFwdEnd(subsize);
			ARM_Vector vSubPayDate(subsize);

			for (int j=0;j<subsize;j++)
			{
				double _start = SchedSubStart->Elt(j);
				double _end = SchedSubEnd->Elt(j);

				CptPeriodDates(_start,_end,resetgap,resettiming,paytiming,GetCurrencyUnit(),
						tenor,fwdadj,paygap,_paydate.GetJulian(),idxindex->GetCM_resetWeekDay(),
						idxindex->GetCM_resetOccur(),paydate,reset,fwdstart,fwdend);
				
				vSubStart.Elt(j) = _start;
				vSubEnd.Elt(j) = _end;
				vSubReset.Elt(j) = reset;
				vSubFwdStart.Elt(j) = fwdstart ;
				vSubFwdEnd.Elt(j) = fwdend ;
				vSubPayDate.Elt(j) = paydate;
			}
			if (SchedSubStart) delete SchedSubStart;
			if (SchedSubEnd) delete SchedSubEnd;

			itsSubStart[i] = vSubStart;
			itsSubEnd[i] = vSubEnd;
			itsSubReset[i] = vSubReset;
			itsSubFwdStart[i] = vSubFwdStart;
			itsSubFwdEnd[i] = vSubFwdEnd;
			itsSubPayDate[i] = vSubPayDate;
		}
	}
	else  // cas weekly
	{
		for (int i=0;i<NbFlows;i++)
		{
			ARM_Date startdate = security->GetAccStartDates().Elt(i);
			ARM_Date enddate = security->GetAccEndDates().Elt(i);
			if (resetgap <0) {
// FIXMEFRED: mig.vc8 (28/05/2007 10:30:19):fabs(int) doesnt exist
				startdate.PreviousBusinessDay(abs(resetgap), GetCurrencyUnit()->GetCcyName());
				enddate.PreviousBusinessDay(abs(resetgap), GetCurrencyUnit()->GetCcyName());
			}else {
				startdate.NextBusinessDay(resetgap, GetCurrencyUnit()->GetCcyName());
				enddate.NextBusinessDay(resetgap, GetCurrencyUnit()->GetCcyName());
			}
			itsSwapStart.Elt(i) = startdate.GetJulian();
			itsSwapEnd.Elt(i) = enddate.GetJulian();

			ARM_Vector resetDates;
			ARM_Vector fwdstart;
			ARM_Vector fwdend;
			// find fixing date in this new period.
			CptResetDateWeekly(startdate.GetJulian(), enddate.GetJulian(), 
						idxindex->GetCM_resetWeekDay(), 
						tenor, GetCurrencyUnit(),fwdadj,
						resetDates,  fwdstart, fwdend);		
			itsSubReset[i] = resetDates;
			itsSubFwdStart[i] = fwdstart;
			itsSubFwdEnd[i] = fwdend;
		}

	}
	
}

/* ****************************************************************************************************************** */
/*!	ExtractBounds
/* ****************************************************************************************************************** 

void ICM_CorridorLeg::ExtractBounds(double startdate,double enddate,int& begin,int& end)
{
	begin = end = -1;
	bool first = true;

	for (int i=0; i<itsSwapStart.size();i++)
	{
		if ((itsSwapStart[i]==startdate) && (itsSwapEnd[i]==enddate))
		{
		 if (first)
			{
			 begin = end = i;
			 first = false;
			}
		 else
			{end = i;}
		}
		else if (itsSwapStart[i]>enddate) 
		{break;}
	}

}
*/
void ICM_CorridorLeg::CptCoupons(ARM_Object* object,const ARM_Date& asof) {

	ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves*>(object);
	if (! model) {
		ICMTHROW(ERR_INVALID_ARGUMENT,"Bad model for Corridor leg compute coupons"); 
	}
	// MarketDataMng from ModelMultiCurves :
	auto_ptr<ICM_MktDataMng> pMng(dynamic_cast<ICM_MktDataMng*>(model->GetMktDataMng()->Clone()));
	// ZC curve will always come from the MMC and not in the marketDatatMng
	pMng.get()->insert(model->GetZeroCurve());
	// string formatedDate = asof.toString('.'); // for format in ICM_Pricer	
	ICM_Pricer_Corridor  pricer_tmp; pricer_tmp.Set(this,(ARM_Model*)pMng.get(),ICM_Parameters(),asof );
	pricer_tmp.ComputePrice(qCMPPRICE);
	// if (pricer_tmp) delete pricer_tmp;
}
