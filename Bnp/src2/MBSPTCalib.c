#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "math.h"
#include "MBSPTCalib.h"
#include "swp_h_all.h"
#include "swp_h_df.h"
#include "utdates.h"
//**********************//
//All rates and yields, as input output, will be quoted in %, i.e., 1% <-> 1
//*********************//

double  VolFactor (     double         a,
		double         t)
{
	double
		x;

	if (fabs (a * t) < MBS_SMALL)    //If a=0 or t=0, A(a,t)=1
		return (1.0);

	x = (1.0 - exp (- a * t)) / a / t;

	return (x);

} //VolFactor

//assume completely unitialized data, in particular, pointers are all NULL
//units have to be days first, then months, in each case the terms must be increasing
//if NbMMTerms > 0 but date_units, MMTerms are NULL, still alloc memory but don't copy
char * MM_Data_Construct( MM_Data * mm_data,
						 YBASIS YearBase, DIFFBASIS diffBase, int delay,
						 int NbMMTerms, DATE_UNIT * date_units, int * MMTerms, double * MMYields,
						 long Ref_Date)
{
	int i;
	int hasFoundMo=0;
	int prev = 0;
	//check: have to be days first, then months, in each the terms must be increasing
	if( YearBase == _YACT) return("Can't handle ACT for MM");
	for(i=0; i < NbMMTerms; ++i)
	{
		if(date_units[i] != DAY && date_units[i] != MONTH ) return("Wrong MM date units");
		if(hasFoundMo)
		{
			if(date_units[i] == DAY) return("MM units has to be DAY followed by MONTH");
			if(MMTerms[i] <= prev ) return("MMTerms have to be increasing");
		}
		else
		{
			if(date_units[i] == MONTH)
			{
				hasFoundMo = 1;
			}
			else
			{
				if(MMTerms[i] <= prev) return("MMTerms have to be increasing");
			}
		}
		prev = MMTerms[i];
	}
	//fill:
	mm_data->YearBase=YearBase;
	mm_data->diffBase=diffBase;
	mm_data->delay = delay;
	mm_data->NbMMTerms = NbMMTerms;
	mm_data->Ref_Date = Ref_Date;
	if(NbMMTerms > 0)
	{
		mm_data->date_units = mbspt_DateUnitVector(0, NbMMTerms-1);
		mm_data->MMTerms = mbspt_ivector(0, NbMMTerms-1);
		mm_data->MMYields = mbspt_dvector(0, NbMMTerms-1);
	}
	if(date_units)
	{
		for(i=0; i < NbMMTerms; ++i) mm_data->date_units[i] = date_units[i];
	}
	if(MMTerms)
	{
		for(i=0; i < NbMMTerms; ++i) mm_data->MMTerms[i] = MMTerms[i];
	}
	if(MMYields)
	{
		for(i=0; i < NbMMTerms; ++i) mm_data->MMYields[i] = MMYields[i];
	}

	return 0;
} //MM_Data_Construct

char * MM_Data_Copy( MM_Data * new_data, MM_Data * old_data )
{
	return(
		MM_Data_Construct( new_data, old_data->YearBase, old_data->diffBase, old_data->delay,
							old_data->NbMMTerms, old_data->date_units, old_data->MMTerms, old_data->MMYields,
							old_data->Ref_Date
					));
};

void MM_Data_Destruct( MM_Data * mm_data )
{
	if( mm_data->NbMMTerms > 0)
	{
		mbspt_free_DateUnitVector(mm_data->date_units, 0, mm_data->NbMMTerms-1);
		mbspt_free_ivector(mm_data->MMTerms, 0, mm_data->NbMMTerms-1);
		mbspt_free_dvector(mm_data->MMYields, 0, mm_data->NbMMTerms-1);
	}
	mm_data->date_units = 0;
	mm_data->MMTerms = 0;
	mm_data->MMYields = 0;
}//char * MM_Data_Destruct

//assume completely unitialized data, in particular, pointers are all NULL
//if NbMMTerms > 0 but date_units, MMTerms are NULL, still alloc memory but don't copy
char * Zero_Data_Construct( Zero_Data * zero_data,
						   int frequency, YBASIS YBase, /*CBASIS CvgBase,*/
						   int NbTimePoints, long * dates, double * DayPoints, double * Rates,
						   long Ref_Date )
{
	int i;
	double prev=0.0;
//check
	if( YBase == _YACT) return("Can't handle ACT for zero rates");
	if( dates )
	{
		for(i=0; i < NbTimePoints; ++i)
		{
			if(DayPoints[i] <= prev) return("Zero_Data timepoints must be increasing");
			prev = DayPoints[i];
		}
	}
//fill:
	zero_data->Frequency = frequency;
	zero_data->YearBase = YBase;
	//zero_data->CvgBase = CvgBase;
	zero_data->NbTimePoints=NbTimePoints;
	zero_data->Ref_Date = Ref_Date;
	if(NbTimePoints>0)
	{
		zero_data->dates=mbspt_lvector(0, NbTimePoints - 1);
		zero_data->DayPoints=mbspt_dvector( 0, NbTimePoints - 1);
		zero_data->TimePoints=mbspt_dvector( 0, NbTimePoints - 1);
		zero_data->Rates=mbspt_dvector(0, NbTimePoints - 1);
	}
	if( dates )
	{
		for(i=0; i < NbTimePoints; ++i)
		{
			zero_data->dates[i] = dates[i];
		}
	}
	if( DayPoints )
		for(i=0; i < NbTimePoints; ++i)
		{
			zero_data->DayPoints[i] = DayPoints[i];
			zero_data->TimePoints[i] = DayPoints[i]/(double) zero_data->YearBase;
		}
	if( Rates )
		for(i=0; i < NbTimePoints; ++i)
		{
			zero_data->Rates[i] = Rates[i];
		}
	return(0);
}//Zero_Data_Construct

void Zero_Data_Destruct( Zero_Data * zero_data )
{
	if(zero_data->NbTimePoints > 0)
	{
		mbspt_free_lvector(zero_data->dates, 0, zero_data->NbTimePoints-1);
		mbspt_free_dvector(zero_data->DayPoints, 0, zero_data->NbTimePoints-1);
		mbspt_free_dvector(zero_data->TimePoints, 0, zero_data->NbTimePoints-1);
		mbspt_free_dvector(zero_data->Rates, 0, zero_data->NbTimePoints-1);
	}
	zero_data->dates = 0;
	zero_data->DayPoints = 0;
	zero_data->TimePoints=0;
	zero_data->Rates=0;
}//Zero_Data_Destruct

char * Zero_Data_Copy( Zero_Data * new_data, Zero_Data * old_data )
{
	return( Zero_Data_Construct( new_data,
		old_data->Frequency, old_data->YearBase,
		old_data->NbTimePoints, old_data->dates, old_data->DayPoints, old_data->Rates,
		old_data->Ref_Date
		));
}//Zero_Data_Copy

char * MBSDF_Data_Construct_FromZeroData(  MBSDF_Data * df_data, Zero_Data * zero_data )
{
	char * err=0;
	if(zero_data==0) return("Invalid zero_data");
	MBSDF_Data_Destruct(df_data);
	df_data->market[0] = '\0';
	df_data->use_zero_data = 1;
	if(err=Zero_Data_Copy(&(df_data->zero_data), zero_data )) return(err);
	return(err);
}//MBSDF_Data_Construct_FromZeroData

char * MBSDF_Data_Construct_FromName(  MBSDF_Data * df_data, char * mkt )
{
	char * err=0;
	if( (mkt == 0) || (mkt[0]=='\0')) return("Invalid market name");
	if(strlen(mkt) >= MAXSTRBUFSZ) return("market name too long");

	MBSDF_Data_Destruct(df_data);
	df_data->use_zero_data = 0;
	strcpy(df_data->market, mkt);
	return(err);

}//MBSDF_Data_Construct

void MBSDF_Data_Destruct( MBSDF_Data * df_data )
{
	df_data->market[0] = '\0';
	Zero_Data_Destruct(&(df_data->zero_data));
	df_data->use_zero_data = 0;
}//MBSDF_Data_Destruct

char * MBSDF_Data_Copy( MBSDF_Data * new_struct, MBSDF_Data * old_struct)
{
	char * err=0;
	MBSDF_Data_Destruct(new_struct);
	if(old_struct->use_zero_data)
	{
		if(err=MBSDF_Data_Construct_FromZeroData(new_struct, &(old_struct->zero_data))) return(err);
	}
	else
	{

		strcpy(new_struct->market,old_struct->market);
	}
	return(err);
}//MBSDF_Data_Copy

//all units will be converted into months
char * Yield_Data_Construct( Yield_Data * yield_data,
							YBASIS YearBase, DIFFBASIS diffBase, int Delay, int Frequency,
							int NbTimePoints, DATE_UNIT * date_units, int * yield_terms, double *Yields,
							long Ref_Date)
{
	int i;
	int prev = 0;
	//convert and check
	if(date_units && yield_terms)
	{
		for(i=0; i < NbTimePoints; ++i)
		{
			if(date_units[i] != MONTH && date_units[i] != YEAR) return("Ilegal yield date units");
			if(date_units[i] == YEAR)
			{
				date_units[i] = MONTH;
				yield_terms[i] *= 12;
			}
			if( yield_terms[i] <= prev ) return("Yield terms must be increasing");
			prev = yield_terms[i];
		}
	}
	//fill
	yield_data->YearBase = YearBase;
	yield_data->diffBase = diffBase;
	yield_data->Delay = Delay;
	yield_data->Frequency = Frequency;
	yield_data->NbTimePoints = NbTimePoints;
	yield_data->Ref_Date = Ref_Date;
	if(NbTimePoints > 0)
	{
		yield_data->date_units = mbspt_DateUnitVector(0,NbTimePoints-1);
		yield_data->yield_terms = mbspt_ivector(0,NbTimePoints-1);
		yield_data->Yields = mbspt_dvector(0, NbTimePoints-1);
	}
	if( date_units )
	{
		for(i=0; i < NbTimePoints; ++i) yield_data->date_units[i] = date_units[i];
	}
	if(yield_terms)
	{
		for(i=0; i < NbTimePoints; ++i) yield_data->yield_terms[i] = yield_terms[i];
	}
	if( Yields )
	{
		for(i=0; i < NbTimePoints; ++i) yield_data->Yields[i] = Yields[i];
	}
	return(0);
}//Yield_Data_Construct

void Yield_Data_Destruct( Yield_Data * yield_data )
{
	if(yield_data->NbTimePoints>0)
	{
		mbspt_free_DateUnitVector(yield_data->date_units, 0, yield_data->NbTimePoints-1);
		mbspt_free_ivector(yield_data->yield_terms, 0, yield_data->NbTimePoints-1);
		mbspt_free_dvector(yield_data->Yields, 0, yield_data->NbTimePoints-1);
	}
	yield_data->date_units=0;
	yield_data->yield_terms=0;
	yield_data->Yields=0;
}//Yield_Data_Destruct

char * Yield_Data_Copy( Yield_Data * new_data, Yield_Data * old_data )
{
	return(
		Yield_Data_Construct( new_data,
			old_data->YearBase, old_data->diffBase, old_data->Delay, old_data->Frequency,
			old_data->NbTimePoints, old_data->date_units, old_data->yield_terms, old_data->Yields,
			old_data->Ref_Date
			));
}//Yield_Data_Copy

//using basic info from Yield_Data such as basis and frequency, and knowing the time to maturity
//generate the sets of start dates, end dates, and coverages
// starts, ends, and cgs are memory-alloc here!!
char * yield_date_gen( int *nb, long ** starts, long ** ends, double ** cvgs,
					  Yield_Data * yield_data, DATE_UNIT date_unit, int term )
{
	int factor = (date_unit == MONTH)? 12:1;
	long firstReset, firstEnd;
	int interval;
	int i;
	int offset;

#ifdef MBSPTCORRECTDAYCT
	int j;
	long y,m,d, m1;
	long dates[2];
#endif//MBSPTCORRECTDAYCT

//check
	if( date_unit != MONTH && date_unit != YEAR ) return( "Illegal date_unit" );
//
	*nb = (term * yield_data->Frequency) /factor;
//assign dates
	firstReset = Date_Add( yield_data->Ref_Date, DAY, yield_data->Delay, 1, 0);
	interval = 12 / yield_data->Frequency;

#ifdef MBSPTCORRECTDAYCT
	firstEnd = Date_Add( yield_data->Ref_Date, MONTH, (* nb) * interval, 0, 1);
	if(!isBusDay(firstEnd))
	{
		Dsplit(firstEnd, &y, &m, &d);
		firstEnd = Date_Add(firstEnd, DAY, 1, 1, 0);
		Dsplit(firstEnd,&y,&m1,&d);
		if(m1>m) firstEnd = Date_Add(firstEnd, DAY, -1, 1, 0);
	}

	firstEnd = Date_Add( firstEnd, MONTH, -(* nb) * interval, 0, 1);
	if(!isBusDay(firstEnd))
	{
		Dsplit(firstEnd, &y, &m, &d);
		firstEnd = Date_Add(firstEnd, DAY, 1, 1, 0);
		Dsplit(firstEnd,&y,&m1,&d);
		if(m1>m) firstEnd = Date_Add(firstEnd, DAY, -1, 1, 0);
	}

#else MBSPTCORRECTDAYCT
	firstEnd = firstReset;
#endif

	if(firstEnd > firstReset)
	{
		offset = 0;
		(*nb) += 1;
	}
	else 
	{
		firstEnd = firstReset;
		offset=1;
	}
//alloc memory
	(*starts) = mbspt_lvector(0, *nb-1);
	(*ends)   = mbspt_lvector(0, *nb-1);
	(*cvgs)   = mbspt_dvector(0, *nb-1);
#ifdef MBSPTCORRECTDAYCT
	i=0;
	dates[0] = Date_Add(firstReset, MONTH, i * interval,          0, 1);
	dates[1] = Date_Add(firstEnd,   MONTH, (i+offset) * interval, 0, 1);
	for(j=0; j < 2; ++j)
	{
		if(!isBusDay(dates[j]))
		{
			Dsplit(dates[j], &y, &m, &d);
			dates[j] = Date_Add(dates[j], DAY, 1, 1, 0);
			Dsplit(dates[j], &y, &m1, &d);
			if(m1 > m) dates[j] = Date_Add(dates[j], DAY, -1, 1, 0);
		}
	}
	(*starts)[i] = dates[0];
	(*ends)[i] = dates[1];

	for(i=1; i < *nb; ++i)
	{
		dates[1] = Date_Add(firstEnd,   MONTH, (i+offset) * interval, 0, 1);
		for(j=1; j < 2; ++j)
		{
			if(!isBusDay(dates[j]))
			{
				Dsplit(dates[j], &y, &m, &d);
				dates[j] = Date_Add(dates[j], DAY, 1, 1, 0);
				Dsplit(dates[j], &y, &m1, &d);
				if(m1 > m) dates[j] = Date_Add(dates[j], DAY, -1, 1, 0);
			}
		}
		(*starts)[i] = (*ends)[i-1];
		(*ends)[i] = dates[1];
	}
#else
	for(i=0; i < *nb; ++i)
	{
		(*starts)[i] = Date_Add(firstReset, MONTH, i * interval, 0, 1);
		(*ends)[i] = Date_Add(firstReset, MONTH, (i+1) * interval, 0, 1);
	}
#endif//MBSPTCORRECTDAYCT

//cvgs
	for(i=0; i < *nb; ++i)
	{

#ifdef MBSPTCORRECTDAYCT
		if( yield_data->YearBase == _YACT ) (*cvgs)[i] = 1.0 /(double) yield_data->Frequency;
		else (*cvgs)[i] = ((double) DateDiff((*starts)[i], (*ends)[i], yield_data->diffBase))/(double) yield_data->YearBase;
#else
		(*cvgs)[i] =  1.0 /(double) yield_data->Frequency;
#endif//MBSPTCORRECTDAYCT

	}
	return(0);
}//yield_date_gen

//use df implicit in zero_data and basis/frequency spec in yield_data to calc YTM
//yield 1% <-> 0.01
double calc_yield( Zero_Data * zero_data, Yield_Data * yield_data, DATE_UNIT date_unit, int term )
{

	int nb, i;
	long * starts=0;
	long * ends=0;
	double * cvgs=0;
	double floatingLeg;
	double ann;
//
	yield_date_gen(&nb, &starts, &ends, &cvgs, yield_data, date_unit, term );
	floatingLeg = DF(zero_data, starts[0], 0) - DF(zero_data,ends[nb-1], 0);
	//calc annuity
	ann = 0.0;
	for( i=0; i < nb; ++i )
	{
		ann += cvgs[i] * DF( zero_data, ends[i], 0); 
	}
///clean up
	mbspt_free_lvector(starts, 0, nb-1);
	mbspt_free_lvector(ends,    0, nb-1);
	mbspt_free_dvector(cvgs,   0, nb-1);
//
	return( floatingLeg/ ann * 100.0);
}//calc_yield

//using basic info such as basis to calc the start and maturity dates of MM rate, also cvg
//memory to be alloc by user,
//start is a single date
char * mm_yield_date_gen( long * start, long ** ends, double ** cvgs, MM_Data * mm_data )
{
	int i;

#ifdef MBSPTCORRECTDAYCT
	long y, m, d, m1;
#endif//MBSPTCORRECTDAYCT

	*start = Date_Add(mm_data->Ref_Date, DAY, mm_data->delay, 1,0);
	for(i=0; i < mm_data->NbMMTerms; ++i)
	{

#ifdef MBSPTCORRECTDAYCT
		if( mm_data->date_units[i] == DAY )
			(*ends)[i] = Date_Add( *start, mm_data->date_units[i], mm_data->MMTerms[i], 1, 0);
		else
		{
			(*ends)[i] = Date_Add( *start, mm_data->date_units[i], mm_data->MMTerms[i], 0, 1);
			if( !isBusDay((*ends)[i]))
			{
				Dsplit((*ends)[i], &y, &m, &d );
				(*ends)[i] = Date_Add((*ends)[i], DAY, 1, 1,0);
				Dsplit( (*ends)[i], &y, &m1, &d );
				if( m1 > m) (*ends)[i] = Date_Add((*ends)[i], DAY, -1, 1, 0);
			}
		}
#else
		(*ends)[i] = Date_Add( *start, mm_data->date_units[i], mm_data->MMTerms[i], 0, 1);
#endif//MBSPTCORRECTDAYCT

		(*cvgs)[i] = DateDiff( *start, (*ends)[i], mm_data->diffBase)/(double) mm_data->YearBase;
	}
	return(0);
}//mm_yield_date_gen

//convert MMYields into zero_rates
//zero_data is completely uninitialized, not even memory is alloc
char * MMYields_To_ZeroRates( Zero_Data * zero_data, MM_Data * mm_data, int frequency, YBASIS YearBase )
{
	char * err;
	int i;
	double lambda = (double) frequency;
	long start, *ends;
	double * day_points;
	double *cvgs, df0, df;
//clean up time_points and Rates
	Zero_Data_Destruct(zero_data);
//copy basic info:
	ends = mbspt_lvector(0, mm_data->NbMMTerms-1);
	cvgs = mbspt_dvector(0, mm_data->NbMMTerms-1);
	day_points = mbspt_dvector(0, mm_data->NbMMTerms-1);

	mm_yield_date_gen(&start, &ends, &cvgs, mm_data);
	for(i=0; i < mm_data->NbMMTerms; ++i)
	{
		day_points[i] = (double) DateDiff(mm_data->Ref_Date, ends[i], _DACT );
	}
	if(err=Zero_Data_Construct( zero_data,
								frequency, YearBase,
								mm_data->NbMMTerms, ends, day_points, 0, mm_data->Ref_Date )) return(err);
//dfs:
	df0 = 1.0/(1.0 + mm_data->MMYields[0]/100.0 *
					((double) DateDiff(zero_data->Ref_Date, start, mm_data->diffBase))/((double) mm_data->YearBase)
			);

	for(i=0; i < zero_data->NbTimePoints; ++i)
	{
		df = 1.0/(1.0 + mm_data->MMYields[i]/100.0 * cvgs[i]);//forward df
		zero_data->Rates[i] = DFToZR( df * df0, zero_data->TimePoints[i], zero_data->Frequency );
	}

	//clean up
	mbspt_free_lvector(ends, 0, mm_data->NbMMTerms-1);
	mbspt_free_dvector(cvgs, 0, mm_data->NbMMTerms-1);
	mbspt_free_dvector(day_points, 0, mm_data->NbMMTerms-1);
	//

	return 0;
}//MMYields_To_ZeroRates

//use linear_interp to get zero rate
double zero_rate_from_t (Zero_Data * zero_data, double t, int extrapolate )
{
	return( linear_interp(zero_data->NbTimePoints, zero_data->TimePoints, zero_data->Rates, t, extrapolate ) );
}//zero_rate_from_t

double zero_rate( Zero_Data * zero_data, long date, int extrapolate )
{
	return( zero_rate_from_t( zero_data,
		(double) (DateDiff(zero_data->Ref_Date, date, _DACT)) / (double) zero_data->YearBase,
		extrapolate));

}//zero_rate

double DFToZR( double df, double t, int freq )
{
	if( freq==0)
		return( 100.0 * log(df)/(-t));
	else
		return( (100.0 * (double) freq) * (pow( df, -1.0/t/((double) freq) ) - 1.0));
}//DFToZR

double zrToDF( double rate, double t, int freq)
{
	if(freq==0)
		return( exp( - rate/100.0 * t ));
	else
		return( pow( 1.0 + rate/100.0 / (double) freq, -t * ((double) freq)));
}//zrToDF

double DFFromT( Zero_Data * zero_data, double t, int extrapolate )
{
	double rate = zero_rate_from_t( zero_data, t, extrapolate);
	return(zrToDF( rate, t, zero_data->Frequency ));
}//DFFromT

double DF( Zero_Data * zero_data, long date, int extrapolate )
{
	double t = (double) (DateDiff(zero_data->Ref_Date, date, _DACT)) / (double) zero_data->YearBase;
	return( DFFromT( zero_data, t, extrapolate ));
}//DF

double MBS_swp_f_df( long start, long end, char * mkt)
{
	Date date1, date2;
	long y,m,d;
	Dsplit(start,&y,&m,&d);
	date1 = date(d,m,y);
	Dsplit(end,&y,&m,&d);
	date2 = date(d,m,y);
	return swp_f_df(date1, date2, mkt);
}//MBS_swp_f_df

double MBSDF( MBSDF_Data * df_data, long start, long end )
{
	if(df_data->use_zero_data)
	{
		return DF(&(df_data->zero_data), end, 1)/DF(&(df_data->zero_data), start, 1);
	}
	return MBS_swp_f_df(start, end, df_data->market);
}//DF_MBSBK

double MBSDF_FromT( MBSDF_Data * df_data, long start, double t)
{
	long d1, d2, date1, date2;
	double df1, df2;
	d1 = (long) (t * 365.0);
	d2 = d1 + 1;
	date1 = Date_Add( start, DAY, d1, 0, 0);
	date2 = Date_Add( start, DAY, d2, 0, 0);
	df1 = MBSDF(df_data, start, date1);
	df2 = MBSDF(df_data, start, date2);
	return LogInterp( (t), (((double) d1)/365.0), ( ((double) d2)/365.0), (df1), (df2) );
}//DF_MBSBK_FromT

//yield_data is both input and output
//convert all date_units in yield_data to months if necessary,
//let interp_interval = 12/interp_frequency (months)
//starting at Ref_Date + Delay + yield_data->yield_terms[0] (months), at intervals of interp_interval, and ending at the last yield_term
//if yeilds are provided, just use it, else, use linterp
//assume that last_date - 0-th date is a multiple of interp_interval

char * Yield_Data_Interpolate( Yield_Data * yield_data, int interp_frequency )
{
	int i;
	double * terms, * yields;
	int newNb, oldNb;
	int n;
	int interp_interval;

	interp_interval = 12 / interp_frequency;
	if( 12 != interp_interval * interp_frequency ) return("Illegal interp_frequency");
//copy old data
	oldNb = yield_data->NbTimePoints;
	terms = mbspt_dvector(0, oldNb-1);
	yields = mbspt_dvector(0, oldNb-1);
	for(i=0; i < oldNb; ++i)
	{
		terms[i] = (double) yield_data->yield_terms[i];
		if( yield_data->date_units[i] == YEAR) terms[i] *= 12.0;
		yields[i] = yield_data->Yields[i];
	}
//check
	n = yield_data->yield_terms[oldNb-1];
	if(yield_data->date_units[oldNb-1] == YEAR) n *= 12;
	if( yield_data->date_units[0] == YEAR) n -= yield_data->yield_terms[0] * 12;
	else n -= yield_data->yield_terms[0];
	if( n != (n/interp_interval) * interp_interval ) return("First and Last Dates do not fit with interp_frequency");
//prepare new data
	newNb = 1 + (round(terms[oldNb-1]) - round(terms[0]))/interp_interval;
	Yield_Data_Destruct(yield_data);
	Yield_Data_Construct(yield_data,
						yield_data->YearBase, yield_data->diffBase, yield_data->Delay, yield_data->Frequency,
						newNb, 0, 0, 0,
						yield_data->Ref_Date );
//interpolate
	yield_data->date_units[0] = MONTH;
	yield_data->yield_terms[0] = round(terms[0]);
	yield_data->Yields[0] = yields[0];
	for(i = 1; i < newNb; ++i)
	{
		yield_data->date_units[i] = MONTH;
		yield_data->yield_terms[i] = yield_data->yield_terms[i-1] + interp_interval;
		yield_data->Yields[i] = linterpFlat(oldNb, terms, yields, (double) yield_data->yield_terms[i]);
	}
// clean up
	mbspt_free_dvector(terms, 0, oldNb-1);
	mbspt_free_dvector(yields, 0, oldNb-1);
	return(0);
}//Yield_Data_Interpolate

//output: zrs (array), day_points (array)
//annuity: both input and output
//k = index of zrs (zero rate) and day_points (number of days from yield_data->Ref_Date) to be filled in
//annuity, as input, is the value of an annuity at yield_data->Ref_Date, with frequency etc., as encoded in yield_data,
//					with last payment on day_points[k-1] (from Ref_Date)
//yield = YTM with maturity day_points[k] = Ref_Date + period * 12/frequency (months)
//zrs[k] = the implied zero-rate
//annuity = value of the annuity with last payment on day_points[k]
//assume k >= 1

void MBS_ZR_BootStrap( int k, int period, double * zrs, double * annuity, double * day_points,
					  double yield, Zero_Data * zero_data, Yield_Data * yield_data)
{
	double t;
	double cvg;
	double df;
	int nb;
	long * starts=0;
	long *ends=0;
	double * cvgs=0;
	double df0;
	//update day_points
	yield_date_gen( &nb, (&starts), (&ends), (&cvgs), yield_data, MONTH, period * 12 /yield_data->Frequency );
	day_points[k] = (double) DateDiff(zero_data->Ref_Date, ends[nb-1], _DACT);
	//coverage for the added payments at end
	cvg = cvgs[nb-1];

	//old annuity enables last df to be calc
	df0 = DF( zero_data, starts[0], 0 );
	df = (100.0 * df0 - yield * *annuity ) / ( 100.0 * df0 + yield * cvg );
	t = day_points[k]/((double) zero_data->YearBase);
	zrs[k] = DFToZR(df,t, zero_data->Frequency);
	//update annuity
	*annuity += cvg * df;
	//clean up
	mbspt_free_lvector(starts, 0, nb-1);
	mbspt_free_lvector(ends, 0, nb-1);
	mbspt_free_dvector(cvgs, 0, nb-1);
}//MBS_ZR_BootStrap

//zero_data will be internally initialized, including memory alloc, in this routine
//mm_data and yield_data are inputs
//assume that yield_data starts with terms  1Y, mm_data ends with  < 1Y
//also assume that mm_data has values at terms equal to all the coupon dates of a yield of 1Y maturity
//interp_freq is the frequency of interpolation in zero-rate-space
char * MBS_Strip( Zero_Data * zero_data, MM_Data * mm_data, Yield_Data * yield_data, int interp_freq,
				 int zero_cpd_freq, YBASIS zero_ybase, long zero_Ref_Date)
{
	Zero_Data initial_zeros=_MBS_ALL_ZEROS;
	Yield_Data local_yield_data=_MBS_ALL_ZEROS;
	char * err;
	double * zrs;
	double annuity, t, cvg;
	double * day_points;
	long date;
	int i, max_step;
	int zero_NbTimePoints;
	int offset;
	long * starts, * ends;
	double * cvgs;
	int nb;

//offset is used because we woudl like to switch between the (wrong) case of treating 1Y as MMYield AND as swap
#ifndef MBSPTCORRECTDAYCT
	offset = 0; //if 1Y is considered MM
#else
	offset = 1;
#endif//MBSPTCORRECTDAYCT

	//set up zero_data
	zero_NbTimePoints = yield_data->yield_terms[ yield_data->NbTimePoints - 1 ];
	if( yield_data->date_units[ yield_data->NbTimePoints - 1 ] == MONTH ) zero_NbTimePoints /= 12;
	zero_NbTimePoints = (zero_NbTimePoints-1) * interp_freq + offset;
	zero_NbTimePoints += mm_data->NbMMTerms;//64
	
	Zero_Data_Destruct( zero_data);
	if(err = Zero_Data_Construct( zero_data, zero_cpd_freq, zero_ybase,
		zero_NbTimePoints, 0, 0, 0,
		zero_Ref_Date )) return(err);

	//get zero rates up to 1Y
	MMYields_To_ZeroRates( &initial_zeros, mm_data, zero_data->Frequency, zero_data->YearBase);

	if(err = Yield_Data_Construct( &local_yield_data,
			yield_data->YearBase,  yield_data->diffBase, yield_data->Delay, yield_data->Frequency,
			yield_data->NbTimePoints + 1 - offset, 0, 0, 0,
			yield_data->Ref_Date )) return(err);

	if(offset==0)//if 1Y is considered MM, then we have to back out 1Y Swap Rate
	{
		local_yield_data.Yields[0] = calc_yield(&initial_zeros, &local_yield_data, MONTH, 12 );
		local_yield_data.date_units[0] = MONTH;
		local_yield_data.yield_terms[0] = 12;
	}

	for(i=1-offset; i < local_yield_data.NbTimePoints; ++i)
	{
		local_yield_data.Yields[i] = yield_data->Yields[i-1+offset];
		local_yield_data.date_units[i] = yield_data->date_units[i-1+offset];
		local_yield_data.yield_terms[i] = yield_data->yield_terms[i-1+offset];
	}

	//interpolate yields
	if(err=Yield_Data_Interpolate( &local_yield_data, local_yield_data.Frequency ))
	{
		Zero_Data_Destruct( &initial_zeros );
		Yield_Data_Destruct(&local_yield_data);
		return(err);
	}
	//prepare initial values for bootstrap:
	max_step = local_yield_data.yield_terms[local_yield_data.NbTimePoints-1];
	if( local_yield_data.date_units[local_yield_data.NbTimePoints-1] == YEAR ) max_step *= 12;
	max_step /= 12;
	max_step *= interp_freq;
	max_step += offset;

	zrs = mbspt_dvector(0, max_step + offset);//allow one more for initial point when time = 0
	day_points = mbspt_dvector(0, max_step + offset);
	annuity = 0.0;
	day_points[0] = 0;

	for(i=1; i <= local_yield_data.Frequency - offset; ++i)
	{
		yield_date_gen(&nb, &starts, &ends, &cvgs, &local_yield_data, MONTH, (i * 12 ) / local_yield_data.Frequency);

		date = ends[nb-1];
		day_points[i] = (double) DateDiff(zero_data->Ref_Date, date, _DACT );
		zrs[i] = zero_rate(&initial_zeros, date, 0 );
		cvg = cvgs[nb-1];
		t = ((double) DateDiff( zero_data->Ref_Date, date, _DACT ))/((double) zero_data->YearBase);
		annuity += cvg * zrToDF(zrs[i], t, zero_data->Frequency);
		/////clean up for next round
		mbspt_free_lvector(starts, 0, nb-1);
		mbspt_free_lvector(ends, 0, nb-1);
		mbspt_free_dvector(cvgs, 0, nb-1);
		starts=0;
		ends=0;
		cvgs=0;
	}
	//bootstrap:
	for(i= local_yield_data.Frequency - offset + 1; i <= max_step; ++i)
	{
		MBS_ZR_BootStrap(i, i, zrs, &annuity, day_points,
			local_yield_data.Yields[ i - local_yield_data.Frequency ], zero_data, &local_yield_data);
	}
	zrs[0] = zrs[1];
	//interpolate in zero_rate space and form new zero_data
	for(i=0; i < initial_zeros.NbTimePoints; ++i)
	{
		zero_data->dates[i] = initial_zeros.dates[i];
		zero_data->DayPoints[i] = initial_zeros.DayPoints[i];
		zero_data->TimePoints[i] = initial_zeros.TimePoints[i];
		zero_data->Rates[i] = initial_zeros.Rates[i];
	}
	for(i=0; i < max_step - interp_freq; ++i)
	{

#ifndef MBSPTCORRECTDAYCT
		date = Date_Add(zero_data->Ref_Date, MONTH, 12 + (i+1-offset)*12/interp_freq, 0, 1);
#else
		date = Date_Add(zero_data->Ref_Date, DAY, (int) day_points[local_yield_data.Frequency - offset + 1 +i], 0, 0);
#endif

		zero_data->dates[i+ initial_zeros.NbTimePoints] = date;
		zero_data->DayPoints[i+ initial_zeros.NbTimePoints] = DateDiff(zero_data->Ref_Date, date, _DACT);
		zero_data->TimePoints[i+ initial_zeros.NbTimePoints] = zero_data->DayPoints[i+ initial_zeros.NbTimePoints]/((double) zero_data->YearBase);
		zero_data->Rates[i+ initial_zeros.NbTimePoints] = linterpFlat(max_step+1, day_points, zrs, zero_data->DayPoints[i+ initial_zeros.NbTimePoints] );
	}
	//clean up
	Zero_Data_Destruct( &initial_zeros );
	Yield_Data_Destruct( &local_yield_data );
	mbspt_free_dvector(zrs, 0, max_step + offset);
	mbspt_free_dvector(day_points, 0, max_step + offset);
	return(0);
}//MBS_Strip

//does memory alloc
char * Vol_Surface_Construct( Vol_Surface * vol_surface,
							 int NbSwaptions, double * expiries, double *maturities, double * vols,
							 int Rate_Freq, VOL_TYPE vol_type, long RefDate)
{
	int i;
	//check:
	if(Dateok(RefDate)) return("Bad RefDate for vol_surface");
	if(Rate_Freq <= 0) return("Bad freq for vol_surface construction");
	if( NbSwaptions <= 0) return("Empty vol_surface in construction");
	for(i=0; i < NbSwaptions; ++i)
	{
		if( (expiries[i] < 0.0) || (maturities[i] <= 0.0) || (vols[i] < 0.0))
			return("Bad inputs for vol_surface construction");
	}
	//clean up
	Vol_Surface_Destruct( vol_surface );
	//set simple attributes
	vol_surface->NbSwaptions = NbSwaptions;
	vol_surface->vol_type = vol_type;
	vol_surface->Rate_Freq = Rate_Freq;
	vol_surface->RefDate = RefDate;
	//alloc
	if( NbSwaptions > 0)
	{
		vol_surface->expiries   = mbspt_dvector(0, NbSwaptions-1);
		vol_surface->maturities = mbspt_dvector(0, NbSwaptions-1);
		vol_surface->vols       = mbspt_dvector(0, NbSwaptions-1);
	}
	//set
	if(expiries)
	{
		for(i=0; i < NbSwaptions; ++i)
		{
			vol_surface->expiries[i] = expiries[i];
		}
	}
	if(maturities)
	{
		for(i=0; i < NbSwaptions; ++i)
		{
			vol_surface->maturities[i] = maturities[i];
		}
	}
	if(vols)
	{
		for(i=0; i < NbSwaptions; ++i)
		{
			vol_surface->vols[i] = vols[i];
		}
	}
	//
	return(0);
}//Vol_Surface_Construct

void Vol_Surface_Destruct( Vol_Surface * vol_surface )
{
	if(vol_surface->NbSwaptions > 0)
	{
		mbspt_free_dvector(vol_surface->expiries,   0, vol_surface->NbSwaptions - 1);
		mbspt_free_dvector(vol_surface->maturities, 0, vol_surface->NbSwaptions - 1);
		mbspt_free_dvector(vol_surface->vols,       0, vol_surface->NbSwaptions - 1);
	}

	vol_surface->expiries   = 0;
	vol_surface->maturities = 0;
	vol_surface->vols       = 0;
}//Vol_Surface_Destruct

char * Vol_Surface_Copy( Vol_Surface * new_data, Vol_Surface * old_data )
{
	return( Vol_Surface_Construct( new_data,
				old_data->NbSwaptions, old_data->expiries, old_data->maturities, old_data->vols,
				old_data->Rate_Freq, old_data->vol_type,
				old_data->RefDate
				));
}//Vol_Surface_Copy

//alloc memory: NbCol, resets, starts, ends
//assume NbCol, resets, starts, ends all NULL as input!!
//NbRow <-> single number = number of swaptions, maxNbCol <-> single number
//NnCol <-> an array = number of dates involved for each swaption
//each row of the matrix are a series of times, measured from vol_surface->RefDate, using ACT/365,
//of each swptn
char * vol_surface_date_gen( int * NbRow, int * maxNbCol,
						  int ** NbCol,
						  double *** resets, double *** starts, double *** ends,
						  Vol_Surface * vol_surface )
{
	int NbSwaptions, maxReset;
	int cpnInterval;
	int i,j;
	long end, start;
	int expMon;

	//check
	if(vol_surface->NbSwaptions == 0 ) return("Vol Surface is empty!");
	if( (*NbCol) || (*resets) || (*starts) || (*ends) ) return("Some output pointers pass NOT NULL!");
	//
	NbSwaptions = vol_surface->NbSwaptions;
	(*NbCol) = mbspt_ivector(0, NbSwaptions-1);
	cpnInterval = 12 / vol_surface->Rate_Freq;
	maxReset=0;
	for(i=0; i < NbSwaptions; ++i)
	{
		(*NbCol)[i] = round( 12.0 * vol_surface->maturities[i] )/cpnInterval;
		if(maxReset < (*NbCol)[i]) maxReset = (*NbCol)[i];
	}

	*NbRow = NbSwaptions;
	*maxNbCol = maxReset;
	(*resets) = mbspt_dmatrix(0, (*NbRow) - 1, 0, (*maxNbCol) - 1);
	(*starts) = mbspt_dmatrix(0, (*NbRow) - 1, 0, (*maxNbCol) - 1);
	(*ends)   = mbspt_dmatrix(0, (*NbRow) - 1, 0, (*maxNbCol) - 1);
	for(i=0; i < (*NbRow); ++i)
	{
		expMon = round( 12.0 * vol_surface->expiries[i]);
		for(j=0; j < (*NbCol)[i]; ++j)
		{
			start = Date_Add(vol_surface->RefDate, MONTH, expMon + j* cpnInterval,0,1);
			end   = Date_Add(vol_surface->RefDate, MONTH, expMon + (j+1)* cpnInterval,0,1);
			(*starts)[i][j] = DateDiff( vol_surface->RefDate, start, _DACT)/365.0;
			(*resets)[i][j] = (*starts)[i][j];
			(*ends)[i][j]   = DateDiff( vol_surface->RefDate, end, _DACT)/365.0;
		}
	}
	return(0);
}//vol_surface_date_gen

char * MBS_BKTree_Construct_FromFile( MBS_BKTree * tree,
			char * data_dir,
			//Zero_Data * zero_data,
			MBSDF_Data * df_data,
			int NbSwaptions, double Swaptions_Mat, double* SwaptionsExp, double * vols,
			double beta, double OAS,
			double BaseIndexValue, double RefiSpread,
			long RefDate)
{
	double * swaptions_mat;
	int i;
	char * err=0;
	FILE * stream;
	char string[MAXSTRBUFSZ];
	VOL_TYPE vol_type;
	int vol_index_freq;
	int PPY;
	double DX_COEFF;
	double XLimit_Coeff, NSigma;
	double MtgRate;
///
	//
	tree->cach=0;
	//
	swaptions_mat = mbspt_dvector(0, NbSwaptions-1);
	for( i=0; i < NbSwaptions; ++i)
		swaptions_mat[i] = Swaptions_Mat;

	if( !(stream=_OpenMBSFile(data_dir, FileBKTree)))
		return(nrerror("Can't open BKTree.dat"));

	fgets(string, 80, stream);
	fscanf(stream, "%ld\n", &vol_type);

	fgets(string, 80, stream);
	fscanf(stream, "%ld\n", &vol_index_freq);

	fgets(string, 80, stream);
	fscanf(stream, "%ld\n", &PPY);

	fgets(string, 80, stream);
	fscanf(stream, "%lf\n", &DX_COEFF);	

	fgets(string, 80, stream);
	fscanf(stream, "%lf\n", &XLimit_Coeff);	

	fgets(string, 80, stream);
	fscanf(stream, "%lf\n", &NSigma);
	
	MtgRate = RefiSpread;
	if(RefiSpread < 0.0 ) MtgRate = BaseIndexValue - RefiSpread;
////
	if(err=MBS_BKTree_Construct( tree, df_data,//zero_data,
			NbSwaptions, SwaptionsExp, swaptions_mat, vols, vol_index_freq, vol_type,
			RefDate,
			PPY, DX_COEFF, beta, OAS, XLimit_Coeff, NSigma,
			BaseIndexValue, MtgRate ))
		goto CLEANUP;
////
CLEANUP:
	mbspt_free_dvector(swaptions_mat,0,NbSwaptions-1);
	return(err);
}//MBS_BKTree_Construct_FromFile


//assume cach has not had any mem alloc to it
//NbZerosCached = max number of fwd zeros to be cached;
//==1 if we only cach the (trivial) pv (==1.0)
char * BKTreeCachStruct_Construct( BKTreeCachStruct * cach, MBS_BKTree * tree )
{
	cach->tree = tree;
	cach->cachOn = 0;
	cach->NbXNodes = 0;
	cach->NbTNodes = 0;
	cach->currIdx_zeros = 0;
	cach->zeros=0;
	cach->currIdx_zeros_OAS=0;
	cach->NbFwdZerosWithOAS = 0;
	cach->zeros_OAS = 0;
	cach->currIdx_DF=0;
	cach->DF=0;
	cach->currIdx_DF_OAS=0;
	cach->currIdx_probs=0;
	cach->pu=0;
	cach->p0=0;
	cach->pd=0;
	return(0);
}//BKTreeCachStruct_Construct

//assume cach has not had any mem alloc to it
//NbZerosCached = max number of fwd zeros to be cached;
// NbZerosWithOASCached = max number of zeros with OAS to be cached
//==1 if we only cach the (trivial) pv (==1.0)
char * BKTreeCachStruct_Reset( BKTreeCachStruct * cach, int NbZerosCached, int NbZerosWithOASCached )
{
	int lastIdx;
	MBS_BKTree * tree;
	if(!cach) return(0);
	if( NbZerosWithOASCached > NbZerosCached ) return("Cannot cach more zeros with OAS han those without" );
/////
	tree = cach->tree;
	BKTreeCachStruct_Destruct(cach);
	cach->cachOn = 1;
///
	lastIdx = tree->EndTimeSlices + tree->cach->NbFwdZerosWithOAS;
	cach->NbXNodes = 2 * (2*lastIdx-1);//much more generous than needed
	cach->NbTNodes = NbZerosCached;
	cach->NbFwdZerosWithOAS = NbZerosWithOASCached;
///
	if( cach->NbXNodes > 0)
	{
		cach->DF     = mbspt_dvector(0, cach->NbXNodes-1);
		cach->DF_OAS = mbspt_dvector(0, cach->NbXNodes-1);
		cach->pu     = mbspt_dvector(0, cach->NbXNodes-1);
		cach->p0     = mbspt_dvector(0, cach->NbXNodes-1);
		cach->pd     = mbspt_dvector(0, cach->NbXNodes-1);

	}
	if( (cach->NbTNodes > 0) && (cach->NbXNodes > 0))
	{
		cach->zeros = mbspt_dmatrix( 0, cach->NbTNodes - 1, 0, cach->NbXNodes - 1);
	}
	if( (cach->NbFwdZerosWithOAS > 0) && (cach->NbXNodes > 0))
	{
		cach->zeros_OAS = mbspt_dmatrix( 0, cach->NbFwdZerosWithOAS - 1,
					0, cach->NbXNodes - 1);
	}
///
	cach->currIdx_zeros = cach->currIdx_zeros_OAS
		= cach->currIdx_DF = cach->currIdx_DF_OAS
		= cach->currIdx_probs 
		= lastIdx;//cach->tree->NbTimeSlices;
	//
	return(0);
}//BKTreeCachStruct_Reset

void BKTreeCachStruct_Destruct( BKTreeCachStruct * cach )
{
	if(cach)
	{
		if( cach->NbXNodes > 0)
		{
			mbspt_free_dvector(cach->DF,     0, cach->NbXNodes-1);
			mbspt_free_dvector(cach->DF_OAS, 0, cach->NbXNodes-1);
			mbspt_free_dvector(cach->pu,     0, cach->NbXNodes-1);
			mbspt_free_dvector(cach->p0,     0, cach->NbXNodes-1);
			mbspt_free_dvector(cach->pd,     0, cach->NbXNodes-1);
		}
		if( (cach->NbTNodes > 0) && (cach->NbXNodes > 0))
		{
			mbspt_free_dmatrix( cach->zeros, 0, cach->NbTNodes - 1, 0, cach->NbXNodes - 1);
		}
		if( (cach->NbFwdZerosWithOAS > 0) && (cach->NbXNodes > 0))
		{
			mbspt_free_dmatrix( cach->zeros_OAS,
								0, cach->NbFwdZerosWithOAS - 1,
								0, cach->NbXNodes - 1);
		}
	//
		cach->DF = 0;
		cach->DF_OAS = 0;
		cach->zeros = 0;
		cach->zeros_OAS = 0;
		cach->pu    = 0;
		cach->p0    = 0;
		cach->pd    = 0;
	}

}//BKTreeCachStruct_Destruct

char * MBS_BKTree_ConstructHelper( MBS_BKTree * tree,
							long RefDate, int PPY, double DX_COEFF, double beta, double OAS,
							double XLimit_Coeff, double NSigma,
							double RefiBaseIndex, double RefiRate)
{
	//checks
	if(Dateok(RefDate)) return("Bad RefDate for BKTRee construction");
	if( PPY <= 0) return("Bad PPY");
	if(DX_COEFF <= 0.0) return("Bad DX_COEFF");
	if( beta < 0.0 ) return("Bad beta");
	if( XLimit_Coeff < 0.0 ) return("Bad XLimit_Coeff");
	if(NSigma < 0.0 ) return("Bad NSigma");
	if(RefiBaseIndex < 0.0) return("Bad RefiBaseIndex");

	///
	MBS_BKTree_Destruct(tree);
///
	tree->RefDate = RefDate;
	tree->PPY = PPY;
	tree->DX_COEFF = DX_COEFF;
	tree->beta = beta;
	tree->OAS = OAS;
	tree->XLimit_Coeff = XLimit_Coeff;
	tree->NSigma = NSigma;
	tree->RefiBaseIndex = RefiBaseIndex;
	tree->RefiRate = RefiRate;
//////
	if(tree->RefiRate < 0.0) tree->RefiRate = tree->RefiBaseIndex - tree->RefiRate;
	tree->NodeMax = ((int) sqrt(XLimit_Coeff * ((double) tree->PPY))) * ((int) NSigma);
	if (tree->NodeMax < 10)  tree->NodeMax = 10000;   // 10000 = Big! 10 = Small! 
	
	tree->cach = (BKTreeCachStruct *) calloc(1,sizeof(BKTreeCachStruct));
	BKTreeCachStruct_Construct( tree->cach, tree );
///
	return(0);
}//MBS_BKTree_ConstructHelper

char * MBS_BKTree_Construct( MBS_BKTree * tree,
							 //dfs
							//Zero_Data * zero_data,
							MBSDF_Data * df_data,
							 //vol_info
							int NbSwaptions, double * expiries, double *maturities, double * vols,
							int Rate_Freq, VOL_TYPE vol_type,
							//others
							long RefDate,
							int PPY, double DX_COEFF, double beta, double OAS, 
							double XLimit_Coeff, double NSigma,
							double RefiBaseIndex, double RefiRate
							)
{
	char * err = 0;
	tree->cach = 0;

	if(Dateok(RefDate))
		return("Bad RefDate in MBS_BKTree_Construct");
	if(err=MBS_BKTree_ConstructHelper(tree, RefDate, PPY, DX_COEFF, beta, OAS, XLimit_Coeff, NSigma,
				RefiBaseIndex, RefiRate))
				return(err);
	if(err=MBSDF_Data_Copy(&(tree->df_data), df_data))
		return(err);
	if(err=Vol_Surface_Construct(&(tree->vol_surface),
				NbSwaptions, expiries, maturities, vols,
				Rate_Freq, vol_type,
				RefDate))
			return(err);
	return(err);
}//MBS_BKTree_Construct

void MBS_BKTree_Destruct( MBS_BKTree * tree)
{
	//Zero_Data_Destruct(&(tree->zero_data));
	MBSDF_Data_Destruct(&(tree->df_data));
	Vol_Surface_Destruct(&(tree->vol_surface));
////
	if( tree->NbTimeSlices > 0)
	{
		mbspt_free_dvector( tree->tree_times,   0, tree->NbTimeSlices-1);
		mbspt_free_dvector( tree->DFS,          0, tree->NbTimeSlices-1);
		mbspt_free_dvector( tree->sigmas,       0, tree->NbTimeSlices-1);
		mbspt_free_ivector( tree->Bottom,       0, tree->NbTimeSlices-1);
		mbspt_free_ivector( tree->Top,          0, tree->NbTimeSlices-1);
		mbspt_free_dvector( tree->Drift,        0, tree->NbTimeSlices-1);
	}

	tree->tree_times  = 0;
	tree->DFS         = 0;
	tree->sigmas      = 0;
	tree->Bottom      = 0;
	tree->Top         = 0;
	tree->Drift       = 0;
	BKTreeCachStruct_Destruct( tree->cach );
	free ((BKTreeCachStruct *) tree->cach);
}//MBS_BKTree_Destruct

//Alloc memory for DFS in tree
//fill in DFS using data from Zero_Data

char * MBS_BKTree_FitCrv( MBS_BKTree * tree  )
{
	int i;
	MBSDF_Data * df_data;
	//
	df_data = &(tree->df_data);
	//check
	if((df_data->use_zero_data) && tree->RefDate != df_data->zero_data.Ref_Date) return("Crv and Tree have different RefDates!");
	//alloc memory
	if(tree->NbTimeSlices>0)
	{
		if(!tree->DFS)//if no memory is alloc
			tree->DFS = mbspt_dvector(0, tree->NbTimeSlices-1);
//		if(!tree->FwdZCpn)//if no memory is alloc
//			tree->FwdZCpn = mbspt_dvector(0, tree->NbTimeSlices-1);
	}
	//DFS
	for(i=0; i < tree->NbTimeSlices; ++i)
		tree->DFS[i] = MBSDF_FromT( df_data,
								tree->RefDate,
								tree->tree_times[i]);
	//FwdZCpn
	return(0);

}//MBS_BKTree_FitCrv

char * MBS_BKTree_FitVol( MBS_BKTree * tree)
{
	char * err = 0;
	int i, j, num;
	double prev, mat;
	int MatMon;
	double * IntVol, *FwdZeroYield, *SigmaValues;
	int NbCpn, cpn_interval;
	double ann, atm_strike, df0, df1;
	double x, y, t1, t2;
	int hasMetZero;
	double maxSigma;
	int NbRow=0, NbCol=0, * colNb;
	double ** resets, ** starts, **ends;
	double DF, DT;
	Vol_Surface * vol_surface;
	MBSDF_Data * df_data;
	//
	vol_surface = &(tree->vol_surface);
	df_data = &(tree->df_data);
//
	num = vol_surface->NbSwaptions;
//check it is vol curve and expiries are increasing
	if(tree->RefDate != vol_surface->RefDate) return("Tree and Vol dates differ!");
	if( num==0) return("Vol_Surface is empty!");
	mat = vol_surface->maturities[0];
	prev = vol_surface->expiries[0];
	for(i=1; i < num; ++i)
	{
		if(vol_surface->maturities[i] != mat ) return("Can only fit a vol curve!!");
		if( prev >= vol_surface->expiries[i] ) return("Vol expiries must be increasing!");
		prev = vol_surface->expiries[i];
	}
//other variables to be set for mem alloc:
	MatMon = round( 12.0 * mat );
	cpn_interval = 12 / vol_surface->Rate_Freq;
	NbCpn = MatMon/cpn_interval;
//alloc mem
	IntVol = mbspt_dvector(0, num);
	FwdZeroYield = mbspt_dvector(0, NbCpn - 1 );
	SigmaValues = mbspt_dvector(0,num);
	colNb = 0; //to silence compiler!!
	resets = 0;
	starts = 0;
	ends = 0;
	//alloc and calc swaption times
	err = vol_surface_date_gen( &NbRow, &NbCol,
		&colNb,
		&resets, &starts, &ends, vol_surface );//num == NbRow!!
	if(err) goto CLEANUP1;
//calc integrated vols, using approximations
	IntVol[0] = 0.0;
	for(i=0; i < num; ++i)
	{
		//first calc various forward rates
//		df0 = DFFromT(zero_data, starts[i][0], 0 );
		df0 = MBSDF_FromT(df_data, tree->RefDate, starts[i][0]);
		df1 = df0;
		ann = 0.0;
		for(j = 0; j < NbCpn; ++j)
		{
//			df1 = DFFromT(zero_data, ends[i][j], 0);
			df1 = MBSDF_FromT(df_data, tree->RefDate, ends[i][j]);
			ann += df1;
			FwdZeroYield[j] = DFToZR( df1/df0, ((double) j+1)/((double) vol_surface->Rate_Freq), vol_surface->Rate_Freq)/100.0;
		}
		ann /= (double) vol_surface->Rate_Freq;
		atm_strike = (df0 - df1)/ann;
		//next, get auxiliary quantities
		x = 0.0;
		y = 0.0;
		t1 = starts[i][0];
		t2 = t1;
		for(j=0; j < NbCpn; ++j)
		{
			t2 = ends[i][j];
			x += atm_strike
				* pow( 1.0 + atm_strike / ((double) vol_surface->Rate_Freq), (double) -(j+2))
				* (j+1)
				* atm_strike;
			y += atm_strike
				* pow( 1.0 + FwdZeroYield[j] / ((double) vol_surface->Rate_Freq), (double) -(j+2))
				* (j+1)
				* FwdZeroYield[j]
				* VolFactor( tree->beta, t2 - t1 );
		}
		x += pow( 1.0 + atm_strike / (double) vol_surface->Rate_Freq, (double) -(NbCpn+1) )
				* (double) NbCpn
				* atm_strike;
		y += pow( 1.0 + FwdZeroYield[NbCpn-1] / (double) vol_surface->Rate_Freq, (double) -(NbCpn+1) )
				* (double) NbCpn
				* FwdZeroYield[NbCpn-1]
				* VolFactor( tree->beta, t2 - t1 );
		IntVol[i+1] = pow( vol_surface->vols[i]/100.0 * x / y, 2.0);
	}
	//from IntVol back out the piecewisely constant sigma
	SigmaValues[0] = 0.0;
	t1 = 0.0;
	hasMetZero = 0;
	maxSigma = 0.0;
	for(i=1; i <= num; i++)
	{
		t2 = starts[i-1][0];
		SigmaValues[i] = (t2 * IntVol[i]
						- t1 * IntVol[i-1] * exp( - 2.0 * tree->beta * (t2 - t1))
						) / (t2 - t1)
						/ VolFactor( 2.0 * tree->beta, t2 - t1);
		if(SigmaValues[i] < MBS_SMALL)
		{
			hasMetZero = 0;
		};
		if(SigmaValues[i] >= MBS_SMALL && hasMetZero)
		{
			err = "Can't handle vol curve somewhere zero and somewher not";
			goto CLEANUP;
		}
		if(SigmaValues[i] > 1.0) 
		{
			err = "Some vol > 1.0";
			goto CLEANUP;
		}
		SigmaValues[i] = sqrt(SigmaValues[i]);
		if( maxSigma < SigmaValues[i] ) maxSigma = SigmaValues[i];
		t1 = t2;
	}
//copy sigmas
	tree->sigmas = mbspt_dvector(0, tree->NbTimeSlices-1);
	i = 1;
	j = 0;
	while( (i < num) && ( j < tree->NbTimeSlices ))
	{
		if(tree->tree_times[j] > starts[i-1][0] - MBS_SMALL)
		{
			++i;
		}
		tree->sigmas[j] = SigmaValues[i];
		++j;
	}
	for( i = j; i < tree->NbTimeSlices; ++i)
	{
		tree->sigmas[i] = SigmaValues[num];
	}
//adjust sigmas as what we have is NOT instantaneous FwdRate but of finite maturity
	if( tree->beta > MBS_SMALL)
	{
		for(i=0; i < tree->NbTimeSlices - 1; ++i)
		{
			DT = tree->tree_times[i+1] - tree->tree_times[i];
			DF = tree->DFS[i]/tree->DFS[i+1];
			tree->sigmas[i] *= DF * log( DF ) / (DF - 1.0 )
				*(1.0
					- exp( -tree->beta * DT )
				  )
				/ ( tree->beta * DT );
		}
	}
	else
	{
		for(i=0; i < tree->NbTimeSlices - 1; ++i)
		{
			DT = tree->tree_times[i+1] - tree->tree_times[i];
			DF = tree->DFS[i]/tree->DFS[i+1];
			tree->sigmas[i] *= DF * log( DF ) / (DF - 1.0 );
		}
	}
	tree->sigmas[i] = tree->sigmas[i-1];
//while we are there, set DX:
	tree->DX = maxSigma * sqrt( tree->DX_COEFF / ((double) tree->PPY) );
//set isZV
	tree->isZV = 1;
	for(i=0; i < vol_surface->NbSwaptions; ++i)
	{
		if( vol_surface->vols[i] >= MBS_SMALL )
		{
			tree->isZV = 0;
			break;
		}
	}
//cleanup
CLEANUP:
	mbspt_free_dvector( IntVol,       0, num);
	mbspt_free_dvector( FwdZeroYield, 0, NbCpn-1);
	mbspt_free_dvector(SigmaValues, 0, num);
	mbspt_free_ivector(colNb, 0, NbRow-1);
	mbspt_free_dmatrix( resets, 0, NbRow-1, 0, NbCol-1);
	mbspt_free_dmatrix( starts, 0, NbRow-1, 0, NbCol-1);
	mbspt_free_dmatrix( ends,   0, NbRow-1, 0, NbCol-1);

	return(err);

CLEANUP1:
	mbspt_free_dvector( IntVol,       0, num);
	mbspt_free_dvector( FwdZeroYield, 0, NbCpn-1);
	mbspt_free_dvector(SigmaValues, 0, num);
	return(err);
}//MBS_BKTree_FitVol


//int  MBS_BKTree_CalcBranching( int i, int NbPer, MBS_BKTree * tree )
//{
//	if( NbPer + 1 <= tree->NodeMax ) return(0);
//	else return( (i == tree->Bottom[NbPer] ) - ( i == tree->Top[NbPer] ));
//}//MBS_BKTree_CalcBranching


char * MBS_BKTree_CalcVolAndDT( double * vol, double * DT, int NbPer, MBS_BKTree * tree )
{
	*DT = tree->tree_times[NbPer+1] - tree->tree_times[NbPer];
	*vol = tree->sigmas[NbPer] * sqrt(*DT);
	return(0);
}//MBS_BKTree_CalcVolAndDT


//pu, p0, pd are shallow copied
char * MBS_BKTree_GetProbs( double ** pu, double ** p0, double ** pd, MBS_BKTree * tree, int NbPer)
{
	double probs[3];
	if(!(tree->cach->cachOn)) return("Cach not on for probs!");
	if(!(tree->cach->currIdx_probs == NbPer))//if we have not cached, calc
	{
		MBS_BKTree_GetProb(probs, 0, tree, NbPer, tree->Top[NbPer]);//to trigger cache
	}
	(*pu) = tree->cach->pu;
	(*p0) = tree->cach->p0;
	(*pd) = tree->cach->pd;
	return(0);
}//MBS_BKTree_GetProbs

//user responsiible for preparing a length 3 array probs (pu, p0, pd)
//if VolAndDT == NULL< infer (Vols, DT) from tree
//else, use the input (Vols, DT)
char * MBS_BKTree_GetProb( double * probs, double * VolAndDT, MBS_BKTree * tree, int NbPer, int i )
{

	double Vol, DT, kappa, alpha;
	int k;
	int start, end, j;
	double * pu=0, *p0=0, *pd=0;
	if( (!(tree->cach->cachOn)) || (!(tree->cach->currIdx_probs == NbPer)))//if we have not cached, calc
	{
		if(VolAndDT)
		{
			Vol = VolAndDT[0];
			DT  = VolAndDT[1];
		}
		else
		{
			MBS_BKTree_CalcVolAndDT(&Vol, &DT, NbPer, tree);	
		}

		if(tree->cach->cachOn)
		{
			start = tree->Bottom[NbPer];
			end   = tree->Top[NbPer];
			tree->cach->currIdx_probs = NbPer;
			pu = tree->cach->pu;
			p0 = tree->cach->p0;
			pd = tree->cach->pd;
		}
		else
		{
			start = end = i;
			pu = probs + 2 - i;
			p0 = probs + 1 - i;
			pd = probs     - i;			
		}
///
		for(j=start; j <= end; ++j)
		{
			if( tree->isZV)
			{
				pu[j]=pd[j]=0.0;
			}
			else
			{
				k = MBS_BKTree_CalcBranching(j, NbPer, tree);
				kappa = (double) k;
				alpha = tree->Drift[NbPer] - tree->beta * DT * ((double) (j-NbPer)) - kappa;
		
				pu[j] = 0.5 * ( 
								pow (Vol / tree->DX, 2.0)
							+   alpha * ( alpha + 1.0 )
							);
				pd[j] = 0.5 * ( 
								pow (Vol / tree->DX, 2.0)
							+   alpha * ( alpha - 1.0 )
							);
			}
			p0[j] = 1.0 - pu[j] - pd[j];
		}
	}
	else
	{
		pu = tree->cach->pu;
		p0 = tree->cach->p0;
		pd = tree->cach->pd;
	}
//copy results
	probs[2] = pu[i];
	probs[1] = p0[i];
	probs[0] = pd[i];
//
	return(0);
}//MBS_BKTree_GetProb

////////////////////////////////////////////////
//Mem Alloc for tree->Bottom, Top, and Drift
//geometry of the lattice: as NbPer increases, (bottom, top) increases as (0, 2 * NbPer), with center the middle
//		until NbPer = some predefined NodeMax, which may be lowered in the course of program
//		thereaftet (bottom, top) = (NbPer - NodeMax, NbPer + NodeMax), again the middle being the center node
//Before the critical point NbPer = NodeMax is reached, always normal bramching, which means:
//		(NbPer, i) --> (NbPer+1,i), (NbPer+1, i+1), (NbPer+1, i+2), [k being the time
//		After the critical point, special branching at top and bottom
//			Top:   i --> i-1, i, i+1
//			Botom: i -->    i + 1, i+2, i+3
//At each point before the critical point is reached, we test for the up probability for top, and down for bottom
//		whenever either of them is close to 1.0, we will decrease NodeMax so that the critical point is immediately reached
///////////////////////////////////////////////
char * MBS_BKTree_SolveDrift( MBS_BKTree * tree )
{

	char * err=0;
	double
		Vol,                                                            // Volatility of the short term rate over the period
		FwdRate1,                                                       // 1 period rate in one period
		x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c;                               // Doubles used in the intermediate calculations of the drift 
	double
		*currDiscount,  //at NbPer it is 1-Perd DF between (tree_times[NbPer], tree_times[NbPer+1])
		*currStatePr,       // State Prices (i.e. Discounted Expected Value of $1 received at this node) 
		*NextDiscount,  //at NbPer it is 1-Perd DF between (tree_times[NbPer+1], tree_times[NbPer])
		*nxtStatePr;      // State Prices in one period 
	int
		Bottom,                                                         // Bottom node index in the tree 
		Top,                                                            // Top node index in the tree 
		CutFlag,                                                        // =TRUE if tree is already cut, =FALSE otherwise 
		i,                                                              // Node index (i.e. vertical) 
		k,                                                              // If k=0 the tree is not cut; if k=1 the tree is cut and branches up or down 
		NbPer;                                                          // Time step index (i.e. horizontal) 
	char string[MAXSTRBUFSZ];
	int NbNodes;//NbTimeNodes
	int NbSpaceNodes;
	int NodeMax, tmpNodeMax;
	double DT;
	int currCenter, nxtCenter;
	double kappa;
	double probs[3];
	double VolAndDT[2];
	double TriggerCutFlag;
////mem alloc:
	NbNodes = tree->NbTimeSlices;
	tree->Top    = mbspt_ivector(0, NbNodes - 1);//add 2 more in case we want to price at 2 time slices from now to make derivative calc
	tree->Bottom = mbspt_ivector(0, NbNodes - 1);
	tree->Drift  = mbspt_dvector(0, NbNodes - 1);
////
	NbSpaceNodes = tree->NbSpaceNodes;
	currDiscount  = mbspt_dvector (0, NbSpaceNodes-1);//we need drift for up to NbPer=NbNodes-3, so need space for up to NbPer=NbNodes-1
	currStatePr   = mbspt_dvector (0, NbSpaceNodes-1);
	NextDiscount = mbspt_dvector (0, NbSpaceNodes-1);
	nxtStatePr  = mbspt_dvector (0, NbSpaceNodes-1);
	NodeMax = tree->NodeMax;
///
	CutFlag     = 0;
	currStatePr[0]  = 1.0;
	currDiscount[0] = tree->DFS[1];
	tree->Drift[ NbNodes-1] = 0.0;
	tree->Drift[NbNodes-2]=0.0;
	//we can calc Drift Only up to the Period (times[NbNodes - 3], times[NbNodes - 2]), i.e., Drift[NbNodes-3];
	//so the last 2 are bogus
	//The last drift will let us calc transition prob from times[NbNodes - 3] to times[NbNodes - 2]; this in turn,
	//allows prices on times[NbNodes - 2] to be backward inducted, but remember times[NbNodes-1] is only an artificial device to allow transitional prob t be calc., s we have what we need
	for( NbPer = 0; NbPer < NbNodes - 2; ++NbPer)
	{
		Bottom = NbPer - NodeMax;
		if( Bottom < 0) Bottom = 0;
		Top = NbPer + NodeMax;
		if(Top > 2*NbPer) Top = 2 * NbPer;
		tree->Bottom[NbPer] = Bottom;
		tree->Top[NbPer]    = Top;
		if(tree->isZV)
		{
			tree->Bottom[NbPer] = tree->Top[NbPer] = (Top + Bottom)/2;
			continue;
		}
		currCenter = NbPer;
		nxtCenter  = NbPer + 1;
		MBS_BKTree_CalcVolAndDT( &Vol, &DT, NbPer, tree ); // Volatility over the period = sigma * sqrt(time step) 

		//get next 1-Prd DF between [NbPer+1, NbPer+2]
		for (i = Bottom; i <= Top + 2; i++) //if next time slice stays pre-critical it is [Bottom, Top+2]; if post-critical [Bottom+1, Top+1]
		{
			FwdRate1 = (tree->DFS[NbPer+1]/tree->DFS[NbPer+2] - 1.0)
						* exp (tree->DX * (i - nxtCenter));
			NextDiscount[i] = 1. / (1. + FwdRate1);
		}  // for i 

		x1 = x2 = x3 = y1 = y2 = y3 = z1 = z2 = z3 = 0.0;// Intermediate calculations for the drift 
		a = b = c = 0.0;

		//if the model is r = exp(x), dx = (a - b*x) dt + sigma dW
		//then drift means a(t) * DT/DX
		//certain combination of currStatePrc, prevStatePrc, and curr probabilities  equal to DF
		//but the probabilities are expressible in terms of drift
		//==> yielding a quadraric equation of drift:
		for (i = Bottom; i <= Top; i++)
		{
			k = MBS_BKTree_CalcBranching( i, NbPer, tree );
			kappa = (double) k;
			// If the number of nodes in the next period is greater than the maximum   
			// number, the tree is cut i.e. the top node (i=Top) branches down         
			// (i.e. to Top+1, Top, Top-1 instead of Top+2, Top+1, Top) and the bottom 
			// node branches up (i.e. to bottom+1, bottom+2, bottom+3).                
			// The shift in branching is given by k.             

			x1 = currStatePr[i] * currDiscount[i] * NextDiscount[i+2+k];
			x2 = currStatePr[i] * currDiscount[i] * NextDiscount[i+2+k] * ((double) (i-NbPer));
			x3 = currStatePr[i] * currDiscount[i] * NextDiscount[i+2+k] * ((double) ((i-NbPer) * (i-NbPer)));
			y1 = currStatePr[i] * currDiscount[i] * NextDiscount[i+k];
			y2 = currStatePr[i] * currDiscount[i] * NextDiscount[i+k] * ((double) (i-NbPer));
			y3 = currStatePr[i] * currDiscount[i] * NextDiscount[i+k] * ((double) ((i-NbPer) * (i-NbPer)));
			z1 = currStatePr[i] * currDiscount[i] * NextDiscount[i+1+k];
			z2 = currStatePr[i] * currDiscount[i] * NextDiscount[i+1+k] * ((double) (i-NbPer));
			z3 = currStatePr[i] * currDiscount[i] * NextDiscount[i+1+k] * ((double) ((i-NbPer) * (i-NbPer)));


			c += 0.5 * (
						pow (Vol / tree->DX, 2.) * ( x1 + y1 - 2.0 * z1 )
					+   tree->beta * DT * ( y2 - x2 )
					+   pow( tree->beta * DT, 2.0 ) * ( x3 + y3 - 2.0 * z3 )
					+   kappa * ( (kappa - 1.0) * x1 + (kappa + 1.0) * y1 - (2.0 * kappa ) * z1 )
					+   2.0 * kappa * tree->beta * DT * ( x2 + y2 - 2.0 * z2 )
					)
				+ z1;

			b += 0.5 * ( x1 - y1 )
				+ tree->beta * DT * ( 2.0 * z2 - y2 - x2 )
				+ kappa * ( 2.0 * z1 - x1 - y1 );//= -2 * kappa * a

			a += 0.5 * ( x1 + y1 - 2.0 * z1 );

		}  // for i 

		c -= tree->DFS[NbPer+2];

		//coeffs for the quadraric equation is formed
		//now check its solvability
		if (b*b - 4.*a*c < 0. || fabs(a) < 1E-14)
		{
			sprintf(string, "Quadratic equation for Drift[%d] has no root", NbPer);
			err = string;
			goto CLEANUP;
		}

		tree->Drift[NbPer] = (-b - sqrt (b*b - 4.*a*c)) / 2.0 / a;

		///update StatePrc etc.
		//first StatePrc
		VolAndDT[0] = Vol;
		VolAndDT[1] = DT;
		for (i = Bottom; i <= Top; i ++)                                // Calculate state prices in one period 
		{
			k = MBS_BKTree_CalcBranching( i, NbPer, tree); //(NbPer + 1 > NodeMax) * ((i == Bottom) - (i == Top));
			MBS_BKTree_GetProb( probs, VolAndDT, tree, NbPer, i );
			nxtStatePr[i+2+k] += probs[2] * currStatePr[i] * currDiscount[i];       // There are 3 branches originating from node i and they contribute to 3 state 
			nxtStatePr[i+1+k] += probs[1] * currStatePr[i] * currDiscount[i];       // prices in one period.                                                       
			nxtStatePr[i+k]   += probs[0] * currStatePr[i] * currDiscount[i];
		}  // for i 

		for (i = Bottom; i <= Top + 4; i++)  // Permute the values so that they can be used for the next iteration of the loop,
			//up to Top+4 because we want to be sure that that nxtStatePrc, etc, 2 steps from now, are all set to zero.
		{                                    // Again we don't use all the nodes if the tree is cut.                            
			currDiscount[i]  = NextDiscount[i];
			currStatePr[i]   = nxtStatePr[i];
			NextDiscount[i]  = nxtStatePr[i] = 0.0;
		}  // for i 

		//check probability
		if (!CutFlag)  // If tree has not been cut yet we check that the probabilities at the top and 
		{      // bottom of the tree are not too big or too small
			TriggerCutFlag = 0;

			tmpNodeMax = tree->NodeMax;
			tree->NodeMax = NbPer + 1;//to force prob calculation using normal branching

			MBS_BKTree_GetProb( probs, VolAndDT, tree, NbPer, 0 );
			if ((probs[2] > 0.999) || (probs[1] < .001) || (probs[0] < .001))  // Bottom of the tree (i=0): pu might be close to 1 and pd close to 0 
			{
				TriggerCutFlag = 1; // We cut the tree at the next period 
			}  // if 

			MBS_BKTree_GetProb( probs, VolAndDT, tree, NbPer, 2 * NbPer );
			if ((probs[2] < .001) || (probs[1] < .001) || (probs[0] > .999))  // Top of the tree (i=2*NbPer): pd might be close to 1 and pu close to 0 
			{
				TriggerCutFlag = 1;
			}  // if 

			tree->NodeMax = tmpNodeMax;

			if(TriggerCutFlag)
			{
				CutFlag = TRUE;
				tree->NodeMax = NbPer + 1;
				NodeMax = tree->NodeMax;
			}
		}  // if 
	}// We need Drift[NbNodes-2] => We need FwdVol[NbNodes-2] and Discount[NbNodes].    
////
CLEANUP:
	mbspt_free_dvector (  nxtStatePr,         0,      NbSpaceNodes-1);
	mbspt_free_dvector (  NextDiscount,       0,      NbSpaceNodes-1);
	mbspt_free_dvector (  currStatePr,        0,      NbSpaceNodes-1);
	mbspt_free_dvector (  currDiscount,       0,      NbSpaceNodes-1);
	return(mbs_err_handler(err,0,0));
}//MBS_BKTree_SolveDrift

//calcs the 1-Prd Disc Factor, with or without OAS for time slice NbPer to NbPer + 1
//returned in the array staring with address *values and size *sz, which may be bigger than needed
//use tree->Bottom and Top to know which are not junk
// values always shallow copy!!
char * MBS_BKTree_Calc1PrdDF( int * sz, double ** values, int useOAS, int NbPer, MBS_BKTree * tree )
{
	double rateDT, multiplier, spreadDT;
	int i;
//handle only cach on case
	if(!((tree->cach)->cachOn)) return("Cach is not on!");
//set sz, *values
//if values in cach, no work needs to be done
//else, reset appropriate index and will proceed to calc
	*sz = (tree->cach)->NbXNodes;
	if(useOAS)
	{
		(*values) = (tree->cach)->DF_OAS;
		if(  (tree->cach)->currIdx_DF_OAS == NbPer ) return(0);
	}
	else
	{
		(*values) = (tree->cach)->DF;
		if( (tree->cach)->currIdx_DF == NbPer ) return(0);
	}

	(tree->cach)->currIdx_DF_OAS = NbPer;
	(tree->cach)->currIdx_DF = NbPer;

//#ifdef MBSPTCONSISTENTOAS
	spreadDT = tree->OAS / 10000.0 * (tree->tree_times[NbPer+1] - tree->tree_times[NbPer]);
//#else
//	spreadDT = tree->OAS / 10000.0 * (tree->tree_times[NbPer+1] - tree->tree_times[NbPer]);
//#endif//MBSPTCONSISTENTOAS

	//
	rateDT = (tree->DFS[NbPer]/tree->DFS[NbPer+1] - 1.0) * exp(tree->DX * ((double) (tree->Bottom[NbPer] - NbPer)));
	multiplier = exp( tree->DX );

	for(i=tree->Bottom[NbPer]; i <= tree->Top[NbPer]; ++i)
	{

//#ifdef MBSPTCONSISTENTOAS
//	#ifdef MBSEXPOAS
		(tree->cach->DF_OAS)[i] = exp(-spreadDT)/(1.0+rateDT);
//	#else
//		(tree->cach->DF_OAS)[i] = 1.0/(1.0+spreadDT) / ( 1.0 + rateDT);
//	#endif//MBSEXPOAS
//#else
//		(tree->cach->DF_OAS)[i] = 1.0 / ( 1.0 + rateDT + spreadDT );
//#endif//MBSPTCONSISTENTOAS

		(tree->cach->DF)[i]     = 1.0 / ( 1.0 + rateDT );
		rateDT *= multiplier;
	}

	return(0);
}//MBS_BKTree_Calc1PrdDF

char * MBS_BKTree_BackwardInduct( double * values, //both in and out
								 int useOAS,
								 int outTimeIdx, MBS_BKTree * tree )
{
	int i, k, sz;
	double * dfs=0, tmp, *pu, *p0, *pd;
	char * err=0;

//	MBS_BKTree_CalcVolAndDT( VolAndDT, VolAndDT + 1, outTimeIdx, tree );
	if(err=MBS_BKTree_Calc1PrdDF( &sz, &dfs, useOAS, outTimeIdx, tree )) return(err);
	MBS_BKTree_GetProbs(&pu, &p0, &pd, tree, outTimeIdx);
	tmp = values[tree->Top[outTimeIdx]-1];
	for(i=tree->Bottom[outTimeIdx]; i < tree->Top[outTimeIdx]; ++i)
	{
		k = MBS_BKTree_CalcBranching(i, outTimeIdx, tree );
		values[i] = dfs[i] * (values[i+2+k] * pu[i] + values[i+1+k] * p0[i] + values[i+k] * pd[i]);//that the middle branch of normal branching is 1 higher than the current node means there is no trouncing of values
	}
//special handling for top, for if there is special branching, one node, i.e., Top-1, has been trounced upon!
	//i = tree->Top[outTimeIdx];
	k = MBS_BKTree_CalcBranching(i, outTimeIdx, tree );
	if(k==0) tmp = values[i+k];
	values[i] = dfs[i] * (values[i+2+k] * pu[i] + values[i+1+k] * p0[i] + tmp * pd[i]);
///
	return(err);
}//MBS_BKTree_BackwardInduct


//zero_prices is shallow copied
//if useOAS:
//(*zero_prices)[i][j] means the price at time <-> NbPer, with maturity times[NbPer + i], at space node j
//i=0, ... *lastFwdIndex
//else
//same but zero_prices have OAS
char * MBS_BKTree_ZeroPrices( int * lastFwdIndex, double *** zero_prices, int useOAS,
							 int NbPer, MBS_BKTree * tree )
{
	char * err;
	double ** zeros=0;
	int lastIndex, local_idx;
	double OAS_Factor = 1.0, OAS_Multiple;
	double OAS_scaled = tree->OAS / 10000.0;
	int i , j;
	//
	if(err=MBS_BKTree_ZeroPricesHelper(&lastIndex, &zeros, NbPer, tree)) return(err);
	local_idx = lastIndex;
	if( local_idx > tree->cach->NbFwdZerosWithOAS - 1 ) local_idx = tree->cach->NbFwdZerosWithOAS - 1;
//
	if(useOAS && (tree->cach->currIdx_zeros_OAS != NbPer)) //if not yet cached
	{
		tree->cach->currIdx_zeros_OAS = NbPer;
		for(i=0; i <= local_idx; ++i)
		{
			for(j=tree->Bottom[NbPer]; j <= tree->Top[NbPer]; ++j)
			{
				tree->cach->zeros_OAS[i][j] = tree->cach->zeros[i][j] * OAS_Factor;
			}


			if( i < local_idx )
			{
				OAS_Multiple = OAS_scaled
									* (tree->tree_times[NbPer + i + 1] 
										- tree->tree_times[NbPer + i]
									   );
			}
			else
				OAS_Multiple = 0.0;

			OAS_Factor *= exp(-OAS_Multiple);


		}
	}
//prepare to return
	if(useOAS)
	{
		*lastFwdIndex = local_idx;
		*zero_prices  = tree->cach->zeros_OAS;
	}
	else
	{
		*lastFwdIndex = lastIndex;
		*zero_prices = tree->cach->zeros;
	}
//
	return(0);
}//MBS_BKTree_ZeroPrices

//zero_prices is shallow copied
//(*zero_prices)[i][j] means the price at time <-> NbPer, with maturity times[NbPer + i], at space node j
//i=0, ... *lastFwdIndex
//get zeros without OAS
char * MBS_BKTree_ZeroPricesHelper( int * lastFwdIndex, double *** zero_prices,
							 int NbPer, MBS_BKTree * tree )
{
	int i;
	double * d_ptr;
	int prev_NbZeros;

	if(!tree->cach->cachOn) return( "Cach not on for zero_prices");
	if( NbPer >= tree->EndTimeSlices) return( "NbPer exceeds bounds");
	if( NbPer < tree->EndTimeSlices-1 && tree->cach->currIdx_zeros != NbPer + 1 && tree->cach->currIdx_zeros != NbPer)
		return("Zero prices must be done in a backward induction fashion");
	///
	*lastFwdIndex = tree->EndTimeSlices - 1 - NbPer;
	prev_NbZeros = (*lastFwdIndex);
	if(*lastFwdIndex > tree->cach->NbTNodes - 1) *lastFwdIndex = tree->cach->NbTNodes - 1;
	if(prev_NbZeros > tree->cach->NbTNodes ) prev_NbZeros = tree->cach->NbTNodes;

////if cached
	if(tree->cach->currIdx_zeros == NbPer) 
	{
		(*zero_prices) = tree->cach->zeros;
		return(0);
	}
//if not cached, calc
	tree->cach->currIdx_zeros = NbPer;
	//last time slice
	if( prev_NbZeros == 0)//tree->EndTimeSlices - 1 + tree->cach->NbFwdZerosWithOAS - 1)
	{
		set_vec_to_val( (tree->cach->zeros)[0], 1.0, 0, 2 * NbPer );
		*zero_prices = tree->cach->zeros;
		return(0);
	}
	for(i=0; i < prev_NbZeros; ++i)
	{
		MBS_BKTree_BackwardInduct( tree->cach->zeros[i], 0, NbPer, tree );
	}

//#ifndef MBSPTCORRECTINITDF
//	if(NbPer==0) {(*zero_prices) = tree->cach->zeros; return(0);}
//#endif//#ifndef MBSPTCORRECTINITDF

	//now permute the zeros by permuting the pointers
	d_ptr = tree->cach->zeros[ prev_NbZeros - 1];
	for(i = prev_NbZeros - 1; i > 0; --i)
	{
		tree->cach->zeros[i] = tree->cach->zeros[i-1];
	}
	//if we have not reached the limit yet
	if( prev_NbZeros < tree->cach->NbTNodes )
	{
		tree->cach->zeros[0] = tree->cach->zeros[prev_NbZeros]; //assign an unoccupied array, to be filled with 1
		tree->cach->zeros[prev_NbZeros] = d_ptr;
	}
	else
	{
		tree->cach->zeros[0] = d_ptr;
	}
	//
	set_vec_to_val( (tree->cach->zeros)[0], 1.0, 0, 2 * NbPer );
	(*zero_prices) = tree->cach->zeros;
//
	return(0);
}//MBS_BKTree_ZeroPricesHelper

//yields is pre-alloc by user
//works only for monthly spaced time slices
//if term exceeds the number of cached zeros available, automatically, shrink term to whatever is available
char * MBS_BKTree_Yield_Approx( double * yields,
							  int freq, double term,
							  int NbPer, MBS_BKTree * tree)
{
	double ** zeros=0;
	int lastFwdIndex=0;
	char * err;
	double floating, ann;
	int i, j, start, end, interval;
	////
	if(tree->PPY != 11) return("Can approx yields only for monthly lattice!");
////
	start = (tree->PPY + 1)/freq;//first coupon
	interval = start;//in months

//#ifndef MBSPTCORRECTINITDF
//	if(NbPer==0) start--;
//#endif//#ifndef MBSPTCORRECTINITDF
	 
	end   = start + (round(term * ((double) freq)) - 1) * interval; //last coupon


	if(err = MBS_BKTree_ZeroPrices( &lastFwdIndex, &zeros, 0, 
								NbPer, tree )) return(err);

	if( lastFwdIndex < end ) end = lastFwdIndex;

	end = (end /(freq * interval) ) * (freq*interval);


	for(i=tree->Bottom[NbPer]; i  <= tree->Top[NbPer]; ++i)
	{
		floating = 1.0 - zeros[end][i];
		ann = 0.0;
		for(j= start; j <= end; j += interval )
		{
			ann += zeros[j][i];
		}
		ann /= ((double) freq);
		yields[ i] = 100.0 * floating / ann;
	}
///
	return(0);
}//MBS_BKTree_Yield_Approx

//assume times has been set
char * MBS_BKTree_Calibrate(
	MBS_BKTree * tree )
{
	char * err=0;
		////fit curve
	if( err = MBS_BKTree_FitCrv( tree))
		return(err);
	//fit vol:
	if( err = MBS_BKTree_FitVol( tree))
		return(err);
	//drift:
	if( err = MBS_BKTree_SolveDrift( tree ))
		return(err);
	//
	return(err);
}//MBS_BKTree_Calibrate

