#ifndef MBSPTCALIB_H
#define MBSPTCALIB_H

#include "MBSPTUtil.h"

double  VolFactor (     double         a,
		double         t);

typedef struct
{
	YBASIS YearBase;	//Money market basis: should be 360 <-> term_data.MMB
	DIFFBASIS diffBase; // ACT
	int delay;
	int NbMMTerms; //Number of MM instruments
	DATE_UNIT * date_units;
	int * MMTerms; // amount of units
	double * MMYields; // 1.0 for 1%
	long Ref_Date; //day when pv = 1
} MM_Data;

char * MM_Data_Construct( MM_Data * mm_data,
						 YBASIS YearBase, DIFFBASIS diffBase, int delay,
						 int NbMMTerms, DATE_UNIT * date_units, int * MMTerms, double * MMYields,
						 long Ref_Date);
void MM_Data_Destruct( MM_Data * mm_data );

char * MM_Data_Copy( MM_Data * new_data, MM_Data * old_data );

typedef struct
{
	int Frequency; //compounding frequency = 0 <=> continuous compounding
	YBASIS YearBase; // 365
	//CBASIS CvgBase; // ACT
	int NbTimePoints; 
	long * dates;
	double * DayPoints; //in days relative to RefDate
	double * TimePoints; //DayPoints /YearBase
	double * Rates; // 0.01 for 1%
	long Ref_Date;
} Zero_Data; //Zero_Data enables us to construct a func df: Days --> R, df(t) = df spanning from Ref_Date to (Ref_Date + t)

char * Zero_Data_Construct( Zero_Data * zero_data,
						   int frequency, YBASIS YearBase,
						   int NbTimePoints, long * dates, double * DayPoints, double * Rates,
						   long Ref_Date );
void Zero_Data_Destruct( Zero_Data * zero_data );

char * Zero_Data_Copy( Zero_Data * new_data, Zero_Data * old_data );

typedef struct
{
	char market[MAXSTRBUFSZ];
	int use_zero_data;//==0 <-> use market
	Zero_Data zero_data;
} MBSDF_Data;

char * MBSDF_Data_Construct_FromZeroData( MBSDF_Data * df_data, Zero_Data * zero_data );

char * MBSDF_Data_Construct_FromName( MBSDF_Data * df_data, char * mkt );

void MBSDF_Data_Destruct( MBSDF_Data * df_data );

char * MBSDF_Data_Copy( MBSDF_Data * new_data, MBSDF_Data * old_data);

typedef struct
{
  YBASIS YearBase; //360
  DIFFBASIS diffBase; //_30
  int Delay; //2 bd
  int Frequency;
  int NbTimePoints;
  DATE_UNIT * date_units;
  int * yield_terms;
  double * Yields;
  long Ref_Date;
} Yield_Data; //Yield_Data enables construction of function Y: Days --> Y(t) = yield to maturity from RefDate to (t + Ref_Date)

char * Yield_Data_Construct( Yield_Data * yield_data,
							YBASIS YearBase, DIFFBASIS diffBase, int Delay, int Frequency,
							int NbTimePoints, DATE_UNIT * date_units, int * yield_terms, double *Yields,
							long Ref_Date);
void Yield_Data_Destruct( Yield_Data * yield_data );

char * Yield_Data_Copy( Yield_Data * new_data, Yield_Data * old_data );

char * yield_date_gen( int *nb, long ** starts, long ** ends, double ** cvgs,
					  Yield_Data * yield_data, DATE_UNIT date_unit, int term );

double calc_yield( Zero_Data * zero_data, Yield_Data * yield_data, DATE_UNIT date_unit, int term );

char * mm_yield_date_gen( long * start, long ** ends, double ** cvgs, MM_Data * mm_data );

char * MMYields_To_ZeroRates( Zero_Data * zero_data, MM_Data * mm_data, int frequency, YBASIS YearBase );//convert MMYields into zero_rates

double zero_rate_from_t (Zero_Data * zero_data, double t, int extrapolate );
double zero_rate( Zero_Data * zero_data, long date, int extrapolate);
double DFToZR( double df, double t, int freq );
double zrToDF( double rate, double t, int freq);
double DFFromT( Zero_Data * zero_data, double t, int extrapolate );
double DF( Zero_Data * zero_data, long date, int extrapolate);

double MBS_swp_f_df( long start, long end, char * mkt);

double MBSDF( MBSDF_Data * df_data, long start, long end );
double MBSDF_FromT( MBSDF_Data * df_data, long start, double t);

char * Yield_Data_Interpolate( Yield_Data * yield_data, int interp_frequency );//intrep_frequency in months
void MBS_ZR_BootStrap( int k, int period, double * zrs, double * annuity, double * day_points,
					  double yield, Zero_Data * zero_data, Yield_Data * yield_data);

char * MBS_Strip( Zero_Data * zero_data, MM_Data * mm_data, Yield_Data * yield_data, int interp_freq,
				 int zero_cpd_freq, YBASIS zero_ybase, long zero_Ref_Date );
typedef enum VOL_TYPE
{
	LOGNORMAL = 0,
	NORMAL
} VOL_TYPE;

typedef struct
{
	int NbSwaptions;
	double * expiries; //2. for a 2 x 10 swaption
	double * maturities; //underlyings 10. for a 2 x 10
	double * vols; // 1 for 1%
	int Rate_Freq;
	VOL_TYPE vol_type;
	long RefDate;
} Vol_Surface;

char * Vol_Surface_Construct( Vol_Surface * vol_surface,
							 int NbSwaptions, double * expiries, double *maturities, double * vols,
							 int Rate_Freq, VOL_TYPE vol_type, long RefDate);
void Vol_Surface_Destruct( Vol_Surface * vol_surface );

char * Vol_Surface_Copy( Vol_Surface * new_data, Vol_Surface * old_data );

char * vol_surface_date_gen( int * NbRow, int * maxNbCol,
						  int ** NbCol,
						  double *** resets, double *** starts, double *** ends,
						  Vol_Surface * vol_surface );

//ZeroCachStruct:
//for caching various zero coupon prices for efficiency's sake
//if the tree has T time-slice and S space-nodes
//then zeros[i][j] has meaning only for j >= currIndex and it means the value at time slice currIndex and space node i
//of a zero-coupon maturiing at time slice j

struct tagMBS_BKTree;
typedef struct tagMBS_BKTree MBS_BKTree;

struct tagBKTreeCachStruct;
typedef struct tagBKTreeCachStruct BKTreeCachStruct;


struct tagBKTreeCachStruct
{
	int cachOn;
	int NbXNodes; // number of space nodes = 2 * tree->NbTimeSlices - 1
	//zeros:
	int NbTNodes; // = max number of (fwd) zeros to be cached
	int currIdx_zeros;
	//int NbZerosCached;
	double ** zeros;
	//zeros with OAS
	int currIdx_zeros_OAS;
	int NbFwdZerosWithOAS; // = max number of (fwd) zeros with OAS to be cached <= NbTNodes
	double ** zeros_OAS;
	//zeros[i,j], for i=0, .. min( (NbTNodes - currIndex_zeros + 1), NbZerosCached ) - 1,
	// j in the suitable continguous subinterval of [0, NbXNodes-1],
	//gives the value at currIdx_zeros-th time slice j-th space node,
	//of a zero coupon bond maturing at (currIdx_zeros + i)-th time slice
	//
	//1Prd DF's with and without OAS
	int currIdx_DF;
	double * DF;
	int currIdx_DF_OAS;
	double * DF_OAS;
	int currIdx_probs;
	double * pu;
	double * p0;
	double * pd;
	//
	MBS_BKTree * tree;
};

char * BKTreeCachStruct_Construct( BKTreeCachStruct * cach, MBS_BKTree * tree);

char * BKTreeCachStruct_Reset( BKTreeCachStruct * cach, int NbZerosCached, int NbZerosWithOASCached);

void BKTreeCachStruct_Destruct( BKTreeCachStruct * cach );

struct tagMBS_BKTree
{
	long RefDate;
	int PPY;//parts per year, 12 for roughly month spacing; for cashflow times interval dT, N = number of time slice intervals = round(PPY*dT); ASSUME dT is close to integer multiple of 1/PPY
	double DX_COEFF;
	double beta;
	double OAS; // 25 for 25bp
	double XLimit_Coeff; //right now 20.0
	double NSigma;//right now 4
	double RefiBaseIndex;
	double RefiRate;
	//Zero_Data zero_data;
	MBSDF_Data df_data;
	Vol_Surface vol_surface;
///////////
	int NodeMax; //we have only initial;it will be updated iteratively in drift calculatin
	double DX; //size of jump in log space
	int NbTimeSlices;
	int NbSpaceNodes;
	int EndTimeSlices; // = NbTimeSlices - 1, we start backward induction for zeros at tree_times[EndTimeSlices-1]
	double * tree_times;//always starts from 0.0 <-> RefDate
	double * DFS; //DFS[i] = values at tree_times[0] of $1 payoff at tree_times[i]; i.e., DFS[0] = 1.0
	double * sigmas;//sigmas[i] = constant sigma value prevailing over [tree_times[i], tree_times[i+1]); if i = NbTimeSlices, treat tree_Times[i+1] = infinity, this last one is most likely not used
	int isZV;//zero-vol
	int * Bottom;
	int * Top;
	double * Drift;
	BKTreeCachStruct * cach;
};

char * MBS_BKTree_ConstructHelper( MBS_BKTree * tree,
							long RefDate, int PPY, double DX_COEFF, double beta, double OAS, 
							double XLimit_Coeff, double NSigma,
							double RefiBaseIndex, double RefiRate);

void MBS_BKTree_Destruct( MBS_BKTree * tree);

char * MBS_BKTree_Construct_FromFile( MBS_BKTree * tree,
			char * data_dir,
			//Zero_Data * zero_data,
			MBSDF_Data * df_data,
			int NbSwaptions, double Swaptions_Mat, double* SwaptionsExp, double * vols,
			double beta, double OAS,
			double BaseIndexValue, double RefiSpread,
			long RefDate);

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
							);

char * MBS_BKTree_FitCrv( MBS_BKTree * tree );

char * MBS_BKTree_FitVol( MBS_BKTree * tree );

//**********MBS_BKTree_CalcBranching***************
//given time slice NbPer, and space node i (the lowest one is 0), and the tree,
//clac the brahcing adjustment needed:
//middle branch is always from (i, NbPer) --> (i+1+k, NbPer + 1), normal branching means k = 0
//int  MBS_BKTree_CalcBranching( int i, int NbPer, MBS_BKTree * tree );
#define MBS_BKTree_CalcBranching( i, NbPer, tree )							\
 ( (NbPer) + 1 <= (tree)->NodeMax ) ? 0 :									\
	( ((i) == (tree)->Bottom[(NbPer)]) - ((i) == (tree)->Top[(NbPer)]) );
//**********END MACRO************

char * MBS_BKTree_CalcVolAndDT( double * vol, double * DT, int NbPer, MBS_BKTree * tree );
char * MBS_BKTree_GetProbs( double ** pu, double ** p0, double ** pd, MBS_BKTree * tree, int NbPer);
char * MBS_BKTree_GetProb( double * probs, double * VolAndDT, MBS_BKTree * tree, int NbPer, int i );
char * MBS_BKTree_SolveDrift( MBS_BKTree * tree );
char * MBS_BKTree_Calc1PrdDF( int * sz, double ** values, int useOAS, int NbPer, MBS_BKTree * tree );
char * MBS_BKTree_BackwardInduct( double * values, //both in and out
								 int useOAS,
								 int outTimeIdx, MBS_BKTree * tree );

char * MBS_BKTree_ZeroPrices( int * lastFwdIndex, double *** zero_prices, int useOAS,
							 int NbPer, MBS_BKTree * tree );

char * MBS_BKTree_ZeroPricesHelper( int * lastFwdIndex, double *** zero_prices,
							 int NbPer, MBS_BKTree * tree );

char * MBS_BKTree_Yield_Approx( double * yields,
							  int freq, double term,
							  int NbPer, MBS_BKTree * tree);
					
char * MBS_BKTree_Calibrate(
	MBS_BKTree * tree );

#endif//MBSPTCALIB_H