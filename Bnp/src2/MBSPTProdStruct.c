#include "string.h"

#include "MBSPTUtil.h"
#include "MBSPTProdStruct.h"

//case sensitive!!
//given code, return corresponding enum
//-1 means unrecognized code

MTGPTACYPROG MBSPT_ProgramLookup( char * code )
{
	if(strcmp( code, "FN" )==0)            return(_FN);
	if(strcmp( code, "FN15" )==0)          return(_FN15);
	if(strcmp( code, "GN" )==0)            return(_GN);
	if(strcmp( code, "GN15" )==0)          return(_GN15);
//
	return(_ILLEGAL);
}//MBSPT_ProgramLookup

char * MBSPT_DealStruct_Construct( MBSPT_DealStruct * deal_struct,
			long DealMaturityDate,
			MTGPTACYPROG AgencyProgram,
			int oterm, int cpnFreq, double cpn,
			int PayDelay
			)
{
	//check
	if(Dateok(DealMaturityDate)) return("Bad DealMaturityDate");
	if(oterm <= 0) return("Bad oterm");
	if( cpnFreq <= 0) return("Bad deal coupon freq");
	if(cpn < 0.0) return( "Bad deal coupon");
	if(PayDelay < 0) return("Bad PayDelay");
	//
	MBSPT_DealStruct_Destruct(deal_struct);
	//
	deal_struct->DealMaturityDate = DealMaturityDate;
	deal_struct->AgencyProgram = AgencyProgram;
	deal_struct->oterm = oterm;
	deal_struct->cpnFreq = cpnFreq;
	deal_struct->cpn = cpn;
	deal_struct->PayDelay = PayDelay;
	return(0);
}//MBSPT_DealStruct_Construct

char * MBSPT_DealStruct_Construct_FromFile( MBSPT_DealStruct * deal_struct,
			char * data_dir,
			long DealMaturityDate,
			int oterm, double cpn,
			int PayDelay
			)
{
	char * err=0;
	FILE * stream;
	char string[MAXSTRBUFSZ];
	MTGPTACYPROG AgcyProg;
	int cpn_freq;
///
	if( !(stream=_OpenMBSFile(data_dir, FilePTDeals)))
		return(nrerror("Can't open PTDeals.dat"));
//
	fgets(string, 80, stream);
	fscanf(stream, "%s\n",string);
	if( (AgcyProg = MBSPT_ProgramLookup(string)) == _ILLEGAL ) return("Bad AgencyProgram");

	fgets(string, 80, stream);
	fscanf(stream, "%d\n", &cpn_freq);	
	if(err=MBSPT_DealStruct_Construct( deal_struct, DealMaturityDate, AgcyProg, oterm, cpn_freq, cpn, PayDelay ))
		return(err);
//
	return(err);
}//MBSPT_DealStruct_Construct_FromFile

void MBSPT_DealStruct_Destruct( MBSPT_DealStruct * deal_struct )
{
	MBSPT_CashFlowStruct_Destruct( &(deal_struct->cashflow_struct) );
}//BSPT_DealStruct_Destruct

char * MBSPT_DealStruct_SetTimes ( MBSPT_DealStruct * deal_struct, long TradeDate, long SettleDate )
{
	//check
	if((Dateok(TradeDate)) || (Dateok(SettleDate)) ||(TradeDate > SettleDate)) return("Bad Trade/Settle dates");
	if(SettleDate > deal_struct->DealMaturityDate)
			return("Bad FwdSettle Date");
	if(month_diff(TradeDate,SettleDate)>1 || month_diff(TradeDate,SettleDate)<0)
		return("Can't handle when Fwd Settle is earlier than trade date or more than 1 mo ahead");
	deal_struct->TradeDate = TradeDate;
	deal_struct->SettleDate = SettleDate;
	return( MBSPT_CashFlowStruct_Construct( &(deal_struct->cashflow_struct),
				TradeDate, deal_struct
			));
}//MBSPT_DealStruct_SetTimes

//this is a helper, not to be used directly
//assume security has wala as of settle
//only the Year and Month of dealMatDate is used, it will always operates as if it is the first of the month, regardless of whether it is a holiday
//duePeriodEndDates starts with the 1st of the month following the first record date on or after settledate
//cachflowDates is the corresponding pay date
//their last are the last pay date and the corresponding 1st of month
char * MBSPT_CashFlowStruct_Construct( MBSPT_CashFlowStruct * cashflow_struct,
		long spotSettleDate, MBSPT_DealStruct * deal_struct)
{
	long yy, mm, dd, cpnDate, /*recDate,*/ dealMatDate;
	int i, n;
	int /*movingWam,*/ cpnInterval, payDelay, oterm, cpnFreq;
	double cpn;
	int payDelayMo, payDay;
	long tmpdate;
	//
	dealMatDate = deal_struct->DealMaturityDate;
	oterm = deal_struct->oterm;
	cpnFreq = deal_struct->cpnFreq;
	cpn = deal_struct->cpn;
	payDelay = deal_struct->PayDelay;
	//check:
	if( (Dateok(spotSettleDate)) || ( spotSettleDate > dealMatDate)) return("Bad settle date");
	//cleanup
	MBSPT_CashFlowStruct_Destruct( cashflow_struct );

	cpnInterval = 12/cpnFreq;

	//recDate should be the last record date strictly before settleDate
	//cpnDate is the 1st of the month following recDate
	n = month_diff(spotSettleDate, dealMatDate);
	Dsplit(spotSettleDate,&yy,&mm,&dd);
	cpnDate=DateYMMDD(yy,mm,1);
	cashflow_struct->numCashFlows = n;


	//fill in others:
	cashflow_struct->DealMaturityDate = dealMatDate;
	cashflow_struct->SettleDate = spotSettleDate;
	cashflow_struct->oterm = oterm;
	cashflow_struct->cpnFreq = cpnFreq;
	cashflow_struct->cpn = cpn;
	cashflow_struct->PayDelay =  payDelay;
	cashflow_struct->NbInitialCashFlowToSkip = month_diff(cashflow_struct->SettleDate, deal_struct->SettleDate);
	//various dates/times
	if( cashflow_struct->numCashFlows > 0)
	{
		cashflow_struct->duePeriodEndDates = mbspt_lvector(0, cashflow_struct->numCashFlows - 1);
		cashflow_struct->CashFlowDates     = mbspt_lvector(0, cashflow_struct->numCashFlows - 1);
		cashflow_struct->duePeriodEndTimes = mbspt_dvector(0, cashflow_struct->numCashFlows - 1);
		cashflow_struct->CashFlowTimes     = mbspt_dvector(0, cashflow_struct->numCashFlows - 1);
	}
	
	payDelayMo = (cashflow_struct->PayDelay/30);
	payDay = cashflow_struct->PayDelay - payDelayMo * 30;
//	Dsplit(cpnDate, &yy, &mm, &dd);
//	cpnDate = DateYMMDD( yy, mm, 1);
	cpnDate = Date_Add(cpnDate, MONTH, 1, 0, 1);//1st duePeriodEndDate fir cashflows owned by buyer
	for( i = 0; i < cashflow_struct->numCashFlows; ++i)
	{
		cashflow_struct->duePeriodEndDates[i] = cpnDate;
		cashflow_struct->CashFlowDates[i] = Date_Add(
				cashflow_struct->duePeriodEndDates[i],
				MONTH, payDelayMo,
				0, 0);
		cashflow_struct->CashFlowDates[i] = Date_Add(
				cashflow_struct->CashFlowDates[i],
				DAY, payDay,
				0,1);
		if(!isBusDay(cashflow_struct->CashFlowDates[i]))
			cashflow_struct->CashFlowDates[i] = Date_Add( cashflow_struct->CashFlowDates[i], DAY, 0, 1,0);
		cashflow_struct->duePeriodEndTimes[i] = ((double) DateDiff( spotSettleDate, cashflow_struct->duePeriodEndDates[i], _DACT ))/365.0;
		cashflow_struct->CashFlowTimes[i]     = ((double) DateDiff( spotSettleDate, cashflow_struct->CashFlowDates[i]   , _DACT ))/365.0;
		cpnDate = Date_Add(cpnDate,MONTH,cpnInterval,0,1);
	}

	//get accrual_factor: 30/360

	n = month_diff(deal_struct->SettleDate, dealMatDate);
	if(n==0)
		cashflow_struct->accrualFactor = 0.0;
	else
	{
		i=0;
		while(cashflow_struct->duePeriodEndDates[i] <= deal_struct->SettleDate) ++i;

		Dsplit(cashflow_struct->duePeriodEndDates[i] , &yy, &mm, &dd);
		tmpdate = DateYMMDD(yy,mm, 1);
		tmpdate = Date_Add(tmpdate,MONTH,-1,0,0);
#ifdef MBSPTCORRECTDAYCT 
	cashflow_struct->accrualFactor = ((double) DateDiff(tmpdate, deal_struct->SettleDate, _30))/360.0;
#else
	cashflow_struct->accrualFactor = ((double) DateDiff(tmpdate, deal_struct->SettleDate, _DACT))/365.0;
#endif//MBSPTCORRECTDAYCT 

	}

	return(0);
}//MBSPT_CashFlowStruct_Construct

void MBSPT_CashFlowStruct_Destruct( MBSPT_CashFlowStruct * cashflow_struct )
{
	if( cashflow_struct->numCashFlows > 0)
	{
		mbspt_free_lvector( cashflow_struct->duePeriodEndDates, 0, cashflow_struct->numCashFlows - 1);
		mbspt_free_dvector( cashflow_struct->duePeriodEndTimes, 0, cashflow_struct->numCashFlows - 1);
		mbspt_free_lvector( cashflow_struct->CashFlowDates,     0, cashflow_struct->numCashFlows - 1);
		mbspt_free_dvector( cashflow_struct->CashFlowTimes,     0, cashflow_struct->numCashFlows - 1);
	}
	cashflow_struct->duePeriodEndDates=0;
	cashflow_struct->duePeriodEndTimes=0;
	cashflow_struct->CashFlowDates=0;
	cashflow_struct->CashFlowTimes=0;
}

//given sch_amort_rate, vector of smm[0..nb-1], coupon (annualized), and df[0..nb_zeros-1][0..nb-1]
//assume io and po's memory have been alloc by user
//add the pv of underlying cashflow of io and pv to the io and pv
//user has to provide smm and sch_amort_rate in a special way at the beginning and end of the backward induction process
//to do it right
//1% coupon <-> 0.01
//1% smm    <-> 0.01 (similar for sch_amort_rate)
//zeros: 0          <-> maturity at first FNM/FH pay day following (i.e., the next 25/15 of a month)
//       nbzeros-1  <-> maturity of (Prepay_Delay + 1)-th FNM/FH pay day following
//       i's        <-> FN/FH paydays in between
//start, end: strat and end indices of space nodes
//IO an dPO are prices with $100 notional
char * MBSPT_AddCashFlowPV( double * io, double * po,
						int start, int end,
						double coupon, 
						double sch_amort_rate, 
						double * smm,
						int nb_zeros, double ** zeros, int NbZerosToSkip )
{
	int i, j;
	double tmp;
	double resid;
	for(i=start; i <= end; ++i)
	{
		resid = (1.0-smm[i])*(1.0-sch_amort_rate);
		//first po
		po[i] = 100.0 * 
					( smm[i] * (1.0 - sch_amort_rate ) * zeros[nb_zeros-1][i] //unscheduled
					+ sch_amort_rate * zeros[0][i]
					)
					+ po[i] * resid; //scheduled pay down
		//then io
		tmp =0.0;
		if(NbZerosToSkip<=0) NbZerosToSkip = 0;
		if( NbZerosToSkip)
		{
			for(j=NbZerosToSkip; j < nb_zeros; ++j ) tmp += zeros[j][i];
		}
		else
			for(j=1; j < nb_zeros; ++j ) tmp += zeros[j][i];

		tmp *= smm[i] * (1.0 - sch_amort_rate );//coupon payments continue until unscheduled paydown occurs
		if( NbZerosToSkip == 0)
			tmp += zeros[0][i];//coupon payment
		tmp *= coupon/12.0;
		io[i] = 100.0 * tmp + resid * io[i];
	};
	////
	return(0);
}//MBSPT_CashFlowPV