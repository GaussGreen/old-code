#ifndef MBSPTCALLER_H
#define MBSPTCALLER_H

#include "MBSPTUtil.h"
#include "MBSPTCalib.h"
#include "MBSPTProdStruct.h"
#include "MBSPrepay.h"
#include "MBSPPFuncs.h"

int MBS_BKTree_create_times( 
		MBS_BKTree * tree,
		MBSPT_DealStruct * deal_struct, 
		MBS_Prepay_Engine * prepay_engine
	);

char * MBS_PricesGetSMM( double ** smm,
			int NbPer,
			MBSPT_CashFlowStruct * cashflowstruct,
			MBS_Prepay_Engine * prepay_engine,
			MBS_BKTree * tree
);

char * MBS_PricesBackwardInduct(
			double * io, double * po,
			int NbPer,
			double coupon,//may be replaced later
			double * smm,
			MBSPT_CashFlowStruct * cashflowstruct,
			MBS_Prepay_Engine * prepay_engine,
			MBS_BKTree * tree
			);

char * mbs_iopo_initialize_prices(double * io, double * po, int start, int end);

char * mbs_iopo_initialize_zeros( MBS_BKTree * tree, int UntilNbPer );

char * MBS_IOPO_Price_Adjust( double * io, double * po,
			MBSPT_DealStruct * deal_struct,
			MBS_Prepay_Engine * prepay_engine,
			MBS_BKTree * tree
							 );

char * MBSPT_IOPO_Price(
	    double * io, double * po,
		MBSPT_DealStruct * deal_struct,		
		MBS_Prepay_Engine * prepay_engine,
		MBS_BKTree * tree
						);

char * MBSPT_MBSDF_Prepare_Helper( MBSDF_Data * df_data,
						  int NbRates, double * rates, 
						int NbMMTerms, DATE_UNIT * mm_date_units, int * MMTerms,
						DATE_UNIT * yield_date_units, int * yield_terms,
						  //mm:
						  YBASIS mm_ybase, DIFFBASIS mm_diffBase, int mm_delay,
						  //swap
						  YBASIS yield_ybase, DIFFBASIS yield_diffBase, int yield_freq, int yield_delay,
						  //zero spec
						  int zero_freq, int interp_freq,
						  long RefDate
						  );

char * MBSPT_MBSDF_Prepare_FromRatesAndFile( MBSDF_Data * df_data,
						char * data_dir,
						int NbRates, double * rates,
						int NbMMTerms, DATE_UNIT * mm_date_units, int * MMTerms,
						DATE_UNIT * yield_date_units, int * yield_terms,
						long RefDate
						  );

char * MBSPT_MBSDF_Prepare_FromRatesAndFile_Special( MBSDF_Data * df_data,
						char * data_dir,
						 double* tsycurve, double* swpSprd,
						long RefDate, long spot
						  );

char * MBSPT_Structures_Prepare_FromFile(
			MBS_BKTree * tree, 	MBSPT_DealStruct * deal_struct, MBS_Prepay_Engine * prepay_engine,
			char * data_dir, 
			long valuation, long FwdSettle,
			MBSDF_Data * df_data,
			double refisprd,
			double Rate10Y, //needed in MBS_Prepay_Engine_SetParams for converting ratio form of refi-incent to diff form
			MTGPTACYPROG AgcyProg,
			int NbSwaptions, double* volcurve, long* volexp, double meanrev,
			long swapmat, long* dealmat, double* coupon, double inoas, long paydelays, 
			double* pptweaks);

char * MBSPT_Pricer( double outputs[3],
					MBS_BKTree * tree, 	MBSPT_DealStruct * deal_struct, MBS_Prepay_Engine * prepay_engine);

char * MBSPT_OAS_Solver(
			double * out,
			double mktPrice,
			char * data_dir, 
			long valuation, long fwdSettle,
			char * mkt,
			double refisprd,
			double Rate10Y, //needed in MBS_Prepay_Engine_SetParams for converting ratio form of refi-incent to diff form
			MTGPTACYPROG AgcyProg,
			int NbSwaptions, double* volcurve, long* volexp, double meanrev,
			long swapmat, long* dealmat, double* coupon, double guessed_oas, long paydelays, 
			double* pptweaks,
			double tolerance,
			int max_trial,
			double oas_inc//in bp
		);

char * MBS_Prepay_GetIndexValue(//int * numRow, int * numCol,
								double ** indices,
								int endIndex,
								MBS_BKTree * tree, MBS_Prepay_Engine * prepay_engine
								);

void MBSPT_Structures_Cleanup(
							  MBSDF_Data * df_data,
							  MBS_BKTree * tree, 
							 MBSPT_DealStruct * deal_struct,
							 MBS_Prepay_Engine * prepay_engine );

char * MBS_Prepay_Projection_ParamsFromFile(
	char * data_dir,
	MTGPTACYPROG AgencyProg, long oterm,
	double gwac,
	double Speedup, double rational,
	double RefiRateSpread, double Rate10Y,
	long YYYYMMDD, long spotYYYYMMDD,
	long  wam, int subgroup, int num_future_rates, double * futureRates,
	int useFixedPrepay,
	int * numProjections, double rates[MAXOTERM], double cpr[MAXOTERM], double principal_payments[MAXOTERM],
	int IndexMat
	);

char * MBS_Prepay_Projection_ParamsExplicit(
	char * data_dir,
	MTGPTACYPROG AgencyProg, long oterm,
	double gwac,
	double Speedup, double rational,
	double FebruaryEffect,
	double * IndexLevel, double * AmortLevel,
	double * SeasoningLevel, double * SeasoiningMat,
	int SFSNb, double * Amorts,
	double Seasonality[12],
	int Delay, int NbFixedPrepay, double * FixedPrepay,
	int nuniform,
	double RefiRateSpread, double Rate10Y,
	long TradeYYYYMMDD, long spotYYYYMMDD,
	long  wam, int subgroup, int num_future_rates, double * futureRates,
	int useFixePrepay,
	int * numProjections, double rates[MAXOTERM], double cpr[MAXOTERM], double principal_payments[MAXOTERM],
	int IndexMat, int EndIndexMat, int IndexFreq
		);

#endif//MBSPTCALLER_H