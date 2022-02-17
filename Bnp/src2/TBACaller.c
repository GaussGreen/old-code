#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "TBAUtil.h"
#include "TBAProdStruct.h"
#include "TBACalib.h"
#include "MBSPrepay.h"
#include "MBSPTProdStruct.h"
#include "TBACaller.h"

char
	Benchmarks[14][5] = {"O/N","1mo"," 2mo"," 3mo"," 6mo"," 1yr"," 2yr"," 3yr"," 4yr"," 5yr"," 7yr","10yr","20yr","30yr"};   /* Benchmark instrument labels */

/*****  mbsf  ***************************************************************
*
*       Main routine.
*		Application Entry Point
*/
char * mbsdf ( char * data_dir, long* valuation, long* sptdays, double *refisprd, double* tsycurve, double* swpSprd, double* volcurve, long* volexp, double *meanrev,
			 long *choicen, long *swapmat, long* dealmat, double* coupon, double *inoas, long *paydelays, 
			 double* pptweaks, double* inprices, double *outputs)
{
HISTORY     history;              /* these structures are defined in TBAProdStruct.h */
MBSDENSITY     density;
TRIGGER     trigger;
DEAL        deal;
PREPAY      prepay;

TERM_DATA       term_data;        /* Global structure of term structure data */
TREE_DATA       tree_data;        /* Global structure of tree data */
DEAL_DATA       deal_data;        /* Global structure of deal data */
PREPAY_DATA     prepay_data;      /* Global structure of prepayment data */

#ifdef TIEOUT

	MBSDF_Data df_data=_MBS_ALL_ZEROS;
	MBS_BKTree tree=_MBS_ALL_ZEROS;
	MBSPT_DealStruct deal_struct=_MBS_ALL_ZEROS;
	MBS_Prepay_Engine prepay_engine=_MBS_ALL_ZEROS;
	//double io, po;
#endif//#TIEOUT


int    i,j;
double multiplier[10], price[20],price0;
double  SFS;            
char * err =0;
		

	term_data.Today = *valuation;
	term_data.Spotdays = *sptdays;
	deal_data.Refirate = *refisprd;

	for (i=0; i<6; i++) term_data.MMYield[i] = tsycurve[i];
	for (i=0; i<8; i++) {
		term_data.SwapYield[i] = tsycurve[6+i]+swpSprd[i]/100.0;
		term_data.YieldS[i] = term_data.SwapYield[i];

	}

	term_data.Beta = *meanrev;
	for (i=0; i<NBSWAPTION; i++) {
		term_data.SwaptionVol[i] = volcurve[i];
		term_data.SwaptionExp[i] = volexp[i];
	}
	deal_data.casenumber = *choicen;
	term_data.IndexMat = *swapmat;
	deal_data.Coupon = coupon[0];
	deal_data.Gwac = coupon[1];
	term_data.OAS = *inoas;
	deal_data.pay_delay = *paydelays;
	deal_data.Exer[1] = dealmat[0];
	deal_data.Term = dealmat[1];
	deal_data.Speedup = pptweaks[0];
	deal_data.rational = pptweaks[1];
	deal_data.Inpio = inprices[0];
	deal_data.Inppo = inprices[1];
	deal_data.collatoral = inprices[2];
		//Manager (&term_data,
		//	&tree_data,
		//	&deal_data,
		//	&prepay_data); 
	for (i=0; i < NTRIALS; i++) {
//		printf ("Calibrating interest rate tree ...\n");

		


		if( err = Manager (
			data_dir,
			&term_data,
			&tree_data,
			&deal_data,
			&prepay_data,
			&history,
			&density)) return( err );
	

		
			/* I/O manager */
		if (i==0) multiplier[0] = deal_data.implvol;            
		if (i==1) {
			if (price[0]>price0) multiplier[1] = 2*deal_data.implvol;     
			/* if need to find implied vol*/
			else multiplier[1] = 0.7*deal_data.implvol;
		}

		for (j=0; j<NBSWAPTION; j++) term_data.SwaptionVol[j] *= multiplier[i];
		

	
		if (deal_data.casenumber == 10) price0 = deal_data.collatoral;      /* given collatoral */
		if (deal_data.casenumber == 11) price0 = deal_data.Inpio;       /* given io */
		if (deal_data.casenumber == 12) price0 = deal_data.Inppo;       /* given po */

		if( err = Input_Check (   &term_data,                                             
		/* Check valilidity of inputs */
				&deal_data,
				&prepay_data,
				&tree_data,
				&history,
				&density)) return( err );
#ifndef TIEOUT
		Zero_6m (&term_data); 
		Time (  term_data.ValueDate,                                            
		
		// Set up the different time arrays and allocate memory for the tree 
			
			tree_data.Ppy,
			deal_data.Maturity,
			deal_data.NbExer,
			deal_data.Exer,
			prepay_data.NbFwd,
			prepay_data.FwdMat,
			deal_data.Strike,
			deal_data.Freq,
			'E',                //* Bogus value for the 'E'uropean or 'A'merican argument
			&tree_data);

		if( err =Build_Tree (    &term_data,                                             
				&tree_data,
				&prepay_data)) return(nrerror(err));

		if(err=impp (    &term_data,                                             /* Main routine for SFS calculation */
				&tree_data,
				&deal_data,
				&history,                   
				/* these structures are defined in dtype.h */
				&density,
				&trigger,
				&deal,
				&prepay,
				&price[i],
				&SFS)) return(nrerror(err));
#else
	if(err=MBSPT_MBSDF_Prepare_FromRatesAndFile_Special((&df_data),
		data_dir,
		tsycurve, swpSprd,
		*valuation, *sptdays)) goto CLEANUP;

	if(err=MBSPT_Structures_Prepare_FromFile( &tree, &deal_struct, &prepay_engine,
			data_dir, *valuation, *sptdays,
			&df_data,
			*refisprd, tsycurve[11]+swpSprd[5]/100.0,
			_FN,
			5, volcurve, volexp,
			*meanrev, *swapmat,
			dealmat, coupon, *inoas, *paydelays, pptweaks))
		goto CLEANUP;


	if(err=copyTermStruct(&tree, &deal_struct, &tree_data))
		goto CLEANUP;


	if(err=impp (    &term_data,                                             /* Main routine for SFS calculation */
				&tree_data,
				&deal_data,
				&history,                   
				/* these structures are defined in dtype.h */
				&density,
				&trigger,
				&deal,
				&prepay,
				&price[i],
				&SFS)) goto CLEANUP;

#endif//#ifndef TIEOUT
	

//////////////////////////

		if (deal_data.Hedge != 'N')

		if(err=Tweak ( SFS,
				&term_data,
				&tree_data,
				&deal_data,
				&prepay_data,
				&history,
				&density,
				&trigger,
				&deal,
				&prepay,
				&price[i])) return(nrerror(err));

		for (j=0; j < 4; j++) outputs[j] = deal_data.output[j];

		if (deal_data.casenumber < 10 || deal_data.casenumber >=13 ) break;             /* if no need for implied vol, get out */
		
		if (fabs(price[i]- price0) < ACURA) {   
//			printf("vol multiplier is %lf \n", multiplier[i]);
//			printf("Price is %lf \n", price[i]);
			break;          /* if good enough get out */
		}



		if (i>=1) multiplier[i+1]= multiplier[i-1]+(price0-price[i-1])*
			(multiplier[i]-multiplier[i-1])/(price[i]-price[i-1]); /* linear interpolation */


			
		for (j=0; j<NBSWAPTION; j++) term_data.SwaptionVol[j] /= multiplier[i];	
		
		Tree_Free (     tree_data.NbNodes, &tree_data);                                     /* Free memory allocated for the tree */



} /* end of i loop */

#ifdef TIEOUT
CLEANUP:
		MBSDF_Data_Destruct(&df_data);
		MBS_BKTree_Destruct(&tree);
		MBS_Prepay_Engine_Cleanup(&prepay_engine);
		MBSPT_DealStruct_Destruct(&deal_struct);
		return(err);
#else//#TIEOUT
	return (0);
#endif

}//mbsdf

/*****  impp ******
/*
* main pricing routine
*/
char *  impp(
		TERM_DATA    *term_data, 
        TREE_DATA    *tree_data,
	    DEAL_DATA    *deal_data,
	    HISTORY      *history,					
        MBSDENSITY      *density,
        TRIGGER      *trigger,
        DEAL	       *deal,
	    PREPAY       *prepay,
	    double       *price,                       /* output price, maybe io,po, or collatoral */
		double       *HedgePrice)	//another output neededd if Hedge != 'N'
{
FILE 	*ofpios;
int 	Top, Bottom;				/* number of urates, as input from results.dat */

double	*io, *po,spd[NTRIALS],rtn[NTRIALS],oas[NTRIALS];
double IO_spd, IO_rtn, PO_spd, PO_rtn;                /* gradients of io. po to speedup and speepness */
double
		* weights,
		*tmpSeasoningMat,
		**tmpAmortLevel,
		**Discount,                                                      /* One period discount factor at each node at period NbPer */
		**DiscountS,                                                     /* One period discount factor including Option Adjusted spread */
		**pu,                                                            /* Array of probabilities in the up state in the interest rate tree */
		**pd,                                                            /* Array of probabilities in the down state in the interest rate tree */
		**p0,                                                            /* Array of probabilities in the flat state in the interest rate tree */
		**ParYield,                                                      /* Par Yield in the lattice */
		**Amort,                                                         /* Array of main amortization factor in the interest rate lattice */
		**IO, **PO,                                                         /* SFS prices in the lattice */
		**Zero,                                                         /* Set of zeros maturing used to calculate par yield */
		***ZeroDelay;							/* Set of zeros at each node with maturity upto deal_data->Delay */
	int
		CpPaid,                                                         /* Is the current period a coupon payment date */
		CouponNb,                                                       /* Coupon Number: = deal_data->Term at the end of the tree */
		F,                                                              /* Frequency of underlying payment as a integer */
		IndexF,                                                         /* Frequency of index payment as a integer */
		TotNbZero,                                                      /* Total number of zeros needed to calculate the par yield index */
		TotNbZero1,	                                                /* Totalnumber of zeros actually needed to calculate the par yield index */
		NbZero,                                                         /* Current number of zeros being priced */
		TotPer,                                                         /* Total number of period in the rate tree including the N additional periods */
		Lockout,                                                        /* Lockout date of the SFS expressed in number of nodes */
		ExEnd,                                                          /* Maturity of the SFS expressed in number of nodes */
		NbPer,                                                          /* Current time period */
		N,                                                              /* Number of additional periods at the beginning of the tree */
		i, j, k,q,qm;
	long
		Date,                                                           /* Current coupon date (YYYYMMDD) */
		Month;                                                          /* Current month of the year (1 to 12) */
	long	today, matday, wam, effmat;
	int 	casenumber = deal_data->casenumber;				/* we have 10 cases to deal with */
	int mo;
	char * err;
	int EndIndexMat;

	MBS_Prepay_Engine prepay_engine={0};

	err = 0; //to prevent compler generating warning
	N = 0;   //FIX                                                               /* We add two artificial periods at the beginning of the tree so that at time 0 */
										/* the tree has 5 nodes. With these 5 nodes we can calculate derivatives.       */
	F        = Conv_Freq (deal_data->Freq);
	CouponNb = deal_data->Term;
	IndexF   = Conv_Freq (term_data->IndexFreq);

		TotPer = tree_data->NbNodes + N;
	if (tree_data->ExEnd>180&& deal_data->Always == 'N') TotPer = tree_data->ExEnd+F + N;                        /* sometimes we need to only build trees 1 yr more than maturity */	
	Lockout = tree_data->ExerPosition[0] + N;
	ExEnd  = tree_data->ExEnd + N;

	TotNbZero = term_data->IndexMat * F;                                    /* If we have a 10yr index and monthly mortgage cash flows we need 10*12 zeros */
										/* Note that the program cannot handle the case where the index has a higher   */
	
	Discount  = mbs_dmatrix (0,1,0, 2 * TotPer);                                    /* coupon frequency than the underlying mortgage.                              */
	DiscountS = mbs_dmatrix (0, 1, 0, 2*TotPer);
	pu        = mbs_dmatrix (0, 1,0, 2 * TotPer);
	pd        = mbs_dmatrix (0,  1, 0, 2*TotPer);
	p0        = mbs_dmatrix (0,  1, 0, 2*TotPer);
	Zero      = mbs_dmatrix (0, TotNbZero, 0, 2 * TotPer);                      /* Note that we allocate one more zero than needed */
	ParYield  = mbs_dmatrix (0,1, 0,2 * ExEnd);      
	
	Amort     = mbs_dmatrix (0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);                                     /* we allocate only that much memory.                               */
	IO        = mbs_dmatrix (0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);        /* We to price SFSNb * 2 structures */
	PO	  = mbs_dmatrix(0, 2*deal_data->SFSNb -1, 0, 2*ExEnd);
	io	  = mbs_dvector (0, NTRIALS);
	po        = mbs_dvector (0, NTRIALS);

	ZeroDelay = (double ***) malloc ((unsigned) (TotPer+1)*sizeof(double*));
	if(!ZeroDelay) merror_register("Can't malloc (double ***) ZeroDelay" );
	else
		for (i=0; i<=1; i++) ZeroDelay[i] = mbs_dmatrix (0, deal_data->Delay, 0, 2*TotPer);                   								/* allocate memory to ***ZeroDelay */	

	today = term_data->Today;
	matday = deal_data->Exer[1];
	today /= 100;
	matday /=100;
	wam = (matday/100 - today/100)*12 + (-today+matday)             /* WAM */
		-(matday/100 - today/100)*100;                           /* WAM minus lockout */
	effmat = wam - deal_data->NbFixedPrepay;

for (q=0; q<NTRIALS; q++) {                             /* do iteration if we have to find implied xxx */
	if(!TweakSCurve(deal_data, today))
	{
		MBS_Prepay_Engine_Cleanup(&prepay_engine);
		mbs_free_dmatrix (IO, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
		mbs_free_dmatrix (PO, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
		mbs_free_dmatrix (Amort, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
		mbs_free_dmatrix (ParYield, 0, 1, 0,2 * ExEnd);
		mbs_free_dmatrix (Zero, 0, TotNbZero, 0, 2 * TotPer);
		mbs_free_dmatrix (p0, 0,  1, 0, 2*TotPer);
		mbs_free_dmatrix (pd, 0,  1, 0, 2*TotPer);
		mbs_free_dmatrix (pu, 0,  1, 0, 2*TotPer);
		mbs_free_dmatrix (DiscountS, 0,  1, 0, 2*TotPer);
		mbs_free_dmatrix (Discount, 0,  1, 0, 2*TotPer);
		for (i=0; i<= 1; i++) mbs_free_dmatrix (ZeroDelay[i], 0, deal_data->Delay, 0, 2*TotPer);
		mbs_free_dvector (io, 0, NTRIALS);
		mbs_free_dvector (po, 0, NTRIALS);
		return( nrerror("deal_data.AmortLevel * Seasonality must be between 0 and 1!"));
	}
	Date   = deal_data->Maturity;                                     /* This is the actual maturity date plus 1yr if actual maturity if 15yr or more AND tree is 1yr beyond it,
otherwise it is ten yr more from final maturity */
	NbZero =  1;

	CouponNb = deal_data->Term;

	if( q > 0 )
		MBS_Prepay_Engine_Cleanup(&prepay_engine);

/*************
	if( err = MBS_Prepay_Engine_Params_Reader( PPFDATADIR, &prepay_engine,
			'F', deal_data->Term,
			 deal_data->Gwac, 
			 deal_data->Speedup, deal_data->rational,
			 FEBEFFECT)) return err;
	if( err = MBS_Prepay_Engine_SetParams( &prepay_engine, deal_data->Refirate, term_data->SwapYield[5])) return err;
************/

//************
				tmpAmortLevel = mbs_dmatrix(0, MAXNUMP-1, 0, NBPOINTS-1);
				for(i=0; i < MAXNUMP; ++i)
				{
					for(j=0; j < NBPOINTS; ++j)
						tmpAmortLevel[i][j] = deal_data->AmortLevel[i][j];
				}
				tmpSeasoningMat = mbs_dvector(0,NUMRAMPS-1);
				for(i=0; i < NUMRAMPS; ++i)
					tmpSeasoningMat[i] = (double) deal_data->SeasoningMat[i];
				EndIndexMat = term_data->IndexMat;
				if( deal_data->Term > 180 && deal_data->Always == 'N') EndIndexMat = 1;
				MBS_Prepay_Engine_Initialize( &prepay_engine,
						MAXNUMP, NBPOINTS, deal_data->IndexLevel, tmpAmortLevel,
						deal_data->Seasonality, FEBEFFECT,
						NUMRAMPS, deal_data->SeasoningLevel, tmpSeasoningMat,
						deal_data->Speedup, deal_data->rational,
						2 * deal_data->SFSNb, deal_data->TriggerRate, deal_data->Amort,
						deal_data->Delay,
						deal_data->NbFixedPrepay, deal_data->InitialSMM,
						_FN, deal_data->Term, deal_data->Gwac,
						 history, density,
						 deal_data->Nbweights,
						 term_data->IndexMat, EndIndexMat, Conv_Freq(term_data->IndexFreq));
				mbs_free_dmatrix(tmpAmortLevel, 0, MAXNUMP-1, 0, NBPOINTS-1);
				mbs_free_dvector(tmpSeasoningMat,0, NUMRAMPS-1);

//**********/


	NbPer  = TotPer;
	do {
		CpPaid = (int) floor (tree_data->Accrued[NbPer-N]);             /* CpPaid is TRUE when a cash flow is paid i.e. Accrued = 1. */
		Bottom = (int) max(0., NbPer - tree_data->NodeMax);
		Top = (int) min(2*NbPer, NbPer+tree_data->NodeMax);


		TotNbZero1 = TotNbZero;
		if ( (deal_data->Term-CouponNb)<(TotNbZero-F) && deal_data->Always == 'N') {
			TotNbZero1 = deal_data->Term-CouponNb +F;  /* first 15 years, 10yr index */
		} /* end if */

		if (q == 0) {
			Lattice (       Discount[0],                    // chnage NbPer to 0                  /* Set up the lattice at the current period */
					DiscountS[0],
					pu[0],
					pd[0],
					p0[0],
					NbPer,
					N,
					term_data,
					tree_data);
			Zero_Calc (     Zero,                                /* Discounted expected value of the set of zeros maturing on each coupon date */
					&NbZero,
					TotNbZero1,
					CpPaid,
					Discount[0],
					pu[0],
					pd[0],
					p0[0],
					NbPer,
					TotPer,
					tree_data);
	
			//for (i=Bottom; i<=Top; i++) {
				
			//	if (CouponNb < Oasramp) 
			//		oasfactor = (double) CouponNb/Oasramp;
			//	else oasfactor = 1.0;

			//	if (term_data->TStree == term_data->TSdiscount) 
			//		 DiscountS[0][i] =
			//	               1./(1./Discount[0][i]+term_data->OAS*oasfactor*tree_data->Length[NbPer-N]/10000.); 
			//	else 
			//		DiscountS[0][i] = 1./(1.+(1./Discount[0][i]-1.)*tree_data->FwdRateS[NbPer-N]
			//			*tree_data->Delta[NbPer-N]                      /* tree rates are different from discount rates */
			//			/tree_data->FwdRate[NbPer-N]+term_data->OAS*oasfactor*tree_data->Length[NbPer-N]/10000.);
			//}

			for (i=0; i <= deal_data->Delay; i++) {
				for (j=Bottom; j<= Top; j++) ZeroDelay[0][i][j] = Zero[i][j]/(1. + term_data->OAS*(double)i*30./365./10000.);                    							/* used for lag in prepayment, so oas is included, approximately */
			}


		} /*end if q==0 */
		
		else if (casenumber >= 13 && casenumber <= 15) {				// if do oas analysis, we must redo discountS */
			for (i=Bottom; i<=Top; i++) {
				if (term_data->TStree == term_data->TSdiscount) 

//#ifdef USECONSISTENTOAS
//					 DiscountS[0][i] /=
//				               1.0/(1.0 +term_data->OAS*tree_data->Length[NbPer-N]/10000.0); 
//#else
					 DiscountS[0][i] =
				               1./(1./Discount[0][i]+term_data->OAS*tree_data->Length[NbPer-N]/10000.); 
//#endif

				else 
					DiscountS[0][i] = 1./(1.+(1./Discount[0][i]-1.)*tree_data->FwdRateS[NbPer-N]
						*tree_data->Delta[NbPer-N]                      /* tree rates are different from discount rates */
						/tree_data->FwdRate[NbPer-N]+term_data->OAS*tree_data->Length[NbPer-N]/10000.);
			}
		} //end else if

		if (NbPer <= ExEnd)                                             /* If we sit before the maturity we do the SFS calculations */
		{
			Month = Date / 10000;                                   /* Current date's year (note the integer division) */
			Month = (Date - Month * 10000) / 100;                   /* Current date's month (see Dsplit in date_sub.c) */
			mo = Month - 1;
			if(mo <=0) mo += 12;
		if (CpPaid)
			{
				if (q==0) {
							Par_Yield (   ParYield[0],               /* Value of the par yield in the lattice. We need to calculate it */
							Zero,                           /* on coupon payment dates only. compute only once*/
							(TotNbZero1/F)*F,
							F,
							IndexF,
							tree_data->NodeMax,
							NbPer);
		 		}


				Main_Prepay (
						Amort,
						ParYield[0],
						CouponNb,
						mo,
						NbPer,
						N,
						&prepay_engine,
						tree_data);

			}  /* if */


			for (i = 0; i < 2 * deal_data->SFSNb; i++) {
				IOPO_Price (  IO[i],                         /* Calculate the IO and PO prices in the lattice */
					PO[i],
					ZeroDelay[0],
					DiscountS[0],              /* Note that the discount factor used for the IO includes the OAS */
					Discount[0],
					pu[0],
					pd[0],
					p0[0],
					Amort[i],
					deal_data->Coupon,
					CouponNb,
					F,
					CpPaid,
					NbPer,
					Lockout,
					ExEnd,
					tree_data,
					deal_data);

			}

 			CouponNb -= CpPaid;                                     /* Decrement CouponNb if one coupon has been paid at this period */
		}

		Date = Nxtmth ( Date,                                           /* We reevaluate the current date on a coupon payment by going backward */
				(long) (- CpPaid * 12 / F),                     /* one coupon period at a time.                                         */
				(long) 1);

		NbPer--;

	} while (NbPer >= N);                                                   /* We stop at NbPer=N as it corresponds to time = today */

	for (i = 0; i < 2 * deal_data->SFSNb; i++)
	{
		for (j = deal_data->NbFixedPrepay; j >= 1; j--)                 /* Calculate the fixed amortization for the next few coupons */
		{                                                               /* Note that we have to go backward to do the calculations correctly */
			IO[i][N] *= (1. - deal_data->FixedPrepay[j-1]);


			for (k = 1; k <= j; k++)
			{// We should not have amortized the coupons
			 // falling between today and the i coupon. Discounting should include OAS
				IO[i][N] += deal_data->Coupon / F * deal_data->FixedPrepay[j-1] * Zero[k][N]
					/ (1.0 + term_data->OAS * ( tree_data->Length[0] + ((double) (k-1)) * 30.0 / 365.0 ) / 10000.0 );
			}
														/* We correct it here.                      */
			PO[i][N] = PO[i][N]*(1. - deal_data->FixedPrepay[j-1])+100.*deal_data->FixedPrepay[j-1]*Zero[j][N]
				   	/ (1.0 + term_data->OAS * ( tree_data->Length[0] + ((double) (k-1)) * 30.0 / 365.0 ) / 10000.0 );
		


		}  /* for j */
		if ( deal_data->invfloater == 'Y' )                   /* inverse floater */
			IO[i][N] -= max(0, deal_data->Coupon/F-  
				term_data->MMYield[0]/F) * (tree_data->Accrued[0] - DiscountS[N][N]);
		else if ( deal_data->invfloater == 'C' )                   /* cap  */
			IO[i][N] -= min(deal_data->Coupon/F,  
				term_data->MMYield[0]/F) * (tree_data->Accrued[0] - DiscountS[N][N]); /* last one is coupon*/
	        else	IO[i][N] -= deal_data->Coupon / F * tree_data->Accrued[0];     /* Subtract today's accrued interests */

	}  /* for i */
	 
	deal->ptime = term_data->Today; 				/* today's date */
	deal->otime = Nxtmth(deal_data->Exer[1], -deal_data->Term, 1);
	deal->term = deal_data->Term;					/*coupon */
	deal->Gwac = deal_data->Gwac;  					/* gross coupon */
	deal->Delay = deal_data->Delay;					/* delay */
	prepay->Speedup = deal_data->Speedup;				/*speedup factor */
	prepay->rational = deal_data->rational;				/*rationality*/
	trigger->ngroup = deal_data->SFSNb;				/*number of sinking fund swaps */
	trigger->slow = deal_data->Amort[0];				/*additional smm of slow group */
    trigger->fast = deal_data->Amort[deal_data->SFSNb];		/* additional smm of fast group*/


	i = trigger->ngroup - 1;
	for (i ; i >= 0; i--) {
		trigger->trigrs[trigger->ngroup-1-i] = deal_data->TriggerRate[i];//+deal_data->Refirate-term_data->SwapYield[5]; //trigger rates are refinance rates
		trigger->slowio[trigger->ngroup-1-i] = IO[i][N];				/* slow io price */
		trigger->slowpo[trigger->ngroup-1-i] = PO[i][N];  			/* slow po price */
		trigger->fastio[trigger->ngroup-1-i] = IO[i+deal_data->SFSNb][N], 	/* fast io price */
		trigger->fastpo[trigger->ngroup-1-i] = PO[i+deal_data->SFSNb][N]; 	/* fast po price */			
	}
               		        
	if (casenumber==1) qm = q - (q/3)*3;						/* 3 is enough for 2d linear interpolation */
	else qm =q-(q/2)*2;									/* 2 is for 1d linear interpolation */

	if (deal_data->SFSNb > 1) {						/* for the case we use new prepayment function */

		weights = mbs_dvector(0,prepay_engine.SFSNb-1);
		if(err=(get_Trigger_Distribution(&prepay_engine, deal->ptime, wam, weights))) return(err);

		io[qm] = po[qm] = 0.0;

		for(i=0; i < deal_data->SFSNb; ++i)
		{
			io[qm] += IO[i][N] * weights[i] + IO[i+deal_data->SFSNb][N] * weights[ i+deal_data->SFSNb];
			po[qm] += PO[i][N] * weights[i] + PO[i+deal_data->SFSNb][N] * weights[ i+deal_data->SFSNb];
		}
		mbs_free_dvector( weights, 0, prepay_engine.SFSNb-1);
	} 
	else {	io[qm] = IO[0][N];			/* for the case where we don't use new prepayment function */
		po[qm] = PO[0][N];
	}
	
	rtn[qm] = deal_data->rational;					/* set rtn spd to new values */
	spd[qm] = deal_data->Speedup;					/* print out current results */
	if (casenumber <10) {
//		printf("IO = %7.4lf\n", io[qm]);
//		printf("PO = %7.4lf\n", po[qm]);
//		printf("Collateral = %7.4lf\n", io[qm]+po[qm]);
//		printf("rationality = %7.4f\n speedup = %7.4f\n", rtn[qm], spd[qm]);
	}


	if (casenumber == 0) break; 				/* if we only want to calculate io,po and coll, it's done. */
	
	if (casenumber == 1) {						/* calibrate both speedup and rationality */
		
		if (fabs(io[qm]-deal_data->Inpio)<=ACURA && fabs(po[qm]-deal_data->Inppo)<=ACURA) {
			rs(io,po,rtn,spd,deal_data->Inpio,deal_data->Inppo, &(deal_data->rational),
						&(deal_data->Speedup), &(IO_spd),&(IO_rtn),&(PO_spd),&(PO_rtn)); /* find the gradients */
			deal_data->Speedup = spd[qm];
			deal_data->rational= rtn[qm]; /*set back to correct value */
		        break; /* if good enough, get out early */
		}
		
		if (q==0) { 
		 	deal_data->rational = rtn[0]*(1.0 + (100./deal_data->Coupon)*((io[0]+po[0])/(deal_data->Inppo+deal_data->Inpio)-1.)); /* guess, coupon number roughly equals to dColl/drtn */
			deal_data->Speedup = spd[0]*(deal_data->Inppo/deal_data->Inpio)/(po[0]/io[0]);
		
		}
		
		if (q==1) {
			deal_data->rational =rtn[1]*(1.0 + (100./deal_data->Coupon)*((io[1]+po[1])/(deal_data->Inppo+deal_data->Inpio)-1.));  
			deal_data->Speedup = spd[1] *(deal_data->Inppo/deal_data->Inpio)/(po[1]/io[1]);			/* guess too */
		}
		
		if (q>=2) rs(io,po,rtn,spd,deal_data->Inpio,deal_data->Inppo,
						&(deal_data->rational),&(deal_data->Speedup),
						&(IO_spd), &(IO_rtn), &(PO_spd), &(PO_rtn));  	    /* 2d linear interpolation when q>=2*/
	}
	if (casenumber == 2) {						/* given io+po and rational, solve for speed */
		
		if (fabs(io[qm]+po[qm]-deal_data->collatoral) <= ACURA/3. ) break;	    /* if good enuogh, get out early */
		
		if (q ==0) deal_data->Speedup = spd[0]*(1.0+50.*((po[0]+io[0])/(deal_data->collatoral)-1.));   /* guess work */
		
		if (q >=1) deal_data->Speedup = spd[0]+(deal_data->collatoral
						- io[0]-po[0])*(spd[1]-spd[0])/(io[1]+po[1]-io[0]-po[0]);   /* linear interpolation */
	}
	if (casenumber == 3) {						/* given Io+po and speed, solve for rational */
		
		if (fabs(io[qm]+po[qm]-deal_data->collatoral) <= ACURA/3. ) break;         /* if good enough, exit now */
		
		if (q ==0) deal_data->rational = rtn[0]*(1.0 + 25.*((io[0]+po[0])
						/(deal_data->collatoral)-1.));		    /* guess */	
		
		if (q >=1) deal_data->rational = rtn[0]+(deal_data->collatoral 
						- io[0]-po[0])*(rtn[1]-rtn[0])/(io[1]+po[1]-io[0]-po[0]);   /* linear interpolation */
	}
	if (casenumber == 4) {						/* given io figure out speedup */
		
		if (fabs(io[qm]-deal_data->Inpio) <= ACURA/3. ) break;   			    /* if good enuogh, get out early */
		
		if (q ==0) deal_data->Speedup = spd[0]*io[0]/deal_data->Inpio;   /* guess work */
		
		if (q >=1) deal_data->Speedup = spd[0]+(deal_data->Inpio - io[0])
						*(spd[1]-spd[0])/(io[1]-io[0]);                             /* linear interpolation */
	}

	if (casenumber == 5) {						/* given io figure out rational */
		
		if (fabs(io[qm]-deal_data->Inpio) <= ACURA/3. ) break;	   			    /* if good enuogh, get out early */
		
		if (q ==0) deal_data->rational =    rtn[0]*(io[0]/deal_data->Inpio);   /*guess work */
		
		if (q >=1) deal_data->rational = rtn[0]+(deal_data->Inpio - io[0])
						*(rtn[1]-rtn[0])/(io[1]-io[0]);                             /* linear interpolation */
	}


	if (casenumber == 6) {						/* given po figure out speedup */
		
		if (fabs(po[qm]-deal_data->Inppo) <= ACURA/3. ) break;	   			    /* if good enuogh, get out early */
		
		if (q ==0) deal_data->Speedup = spd[0]*deal_data->Inppo/po[0];   /* guess work */
		
		if (q >=1) deal_data->Speedup = spd[0]+(deal_data->Inppo - po[0])
						*(spd[1]-spd[0])/(po[1]-po[0]);                             /* linear interpolation */
	}

	if (casenumber == 7) {						/* given po figure out rational */
		
		if (fabs(po[qm]-deal_data->Inppo) <= ACURA/3. ) break;	   			    /* if good enuogh, get out early */
		
		if (q ==0) deal_data->rational =    rtn[0]*deal_data->Inppo/po[0]; /* guess */
		
		if (q >=1) deal_data->rational = rtn[0]+(deal_data->Inppo - po[0])
						*(rtn[1]-rtn[0])/(po[1]-po[0]);                             /* linear interpolation */
	}

	if (casenumber == 8) {						/* given po/io figure out speedup */
		
		if (fabs((po[qm]/io[qm])/deal_data->stripratio-1.0) 
						<= ACURA/30.)  break;	            /* if good enuogh, get out early */
		
		if (q ==0) deal_data->Speedup = spd[0]*deal_data->stripratio/(po[0]/io[0]);   /* guess work */
		
		if (q >=1) deal_data->Speedup = spd[0]+(deal_data->stripratio - po[0]/io[0])
						*(spd[1]-spd[0])/(po[1]/io[1]-po[0]/io[0]);                  /* linear interpolation */
	}

	if (casenumber == 9) {						/* given po/io figure out rational */
		
		if (fabs((po[qm]/io[qm])/deal_data->stripratio-1.0) 
						<= ACURA/30. ) break;	            /* if good enuogh, get out early */
		
		if (q ==0) deal_data->rational = rtn[0]*deal_data->stripratio/(po[0]/io[0]);   /* guess work */
		
		if (q >=1) deal_data->rational = rtn[0]+(deal_data->stripratio - po[0]/io[0])
						*(rtn[1]-rtn[0])/(po[1]/io[1]-po[0]/io[0]);                  /* linear interpolation */
	
	}

	if (casenumber == 10) {				/* implied vol for collatoral */
		*price = io[0]+po[0];
		break;
	}

	if (casenumber == 11) {				/* implied vol for io */
		*price = io[0];	
		break;
	}

	if (casenumber == 12) {				/* implied vol for po */
		*price = po[0];
		break;
	}

	if (casenumber == 13) {                         /* oas analysis */
		if (fabs(io[qm]+po[qm]- deal_data->collatoral) < ACURA/3.0) {
//		 	printf("oas is %lf.\n", term_data->OAS);
//			printf("IO is %lf, PO is %lf.\n", io[qm], po[qm]);
			break; 			/* if good enough get out */
		}	
		
		oas[qm] = term_data->OAS;		
		
		if (q==0) {									/*first try */
			if ((io[0]+po[0])> deal_data->collatoral)  term_data->OAS += 20.;      /* if calculated price too large, increase oas */
				else term_data->OAS -= 20.;       /* if calculated price too small, decrease oas */
		}
		
		if (q>=1) term_data->OAS = oas[0]+(deal_data->collatoral-io[0]-po[0])
			*(oas[1]-oas[0])/(io[1]+po[1]-io[0]-po[0]);				 /* linear interpolation */			
	}  /* end case 13 */

	if (casenumber == 14) {                         /*  io oas analysis */
		if (fabs(io[qm]-deal_data->Inpio) < ACURA/3.0) {
//		 	printf("IO oas is %lf.\n", term_data->OAS);
//			printf("IO is %lf, PO is %lf, collatoral is %lf.\n", io[qm], po[qm], io[qm]+po[qm]);
			break; 			/* if good enough get out */
		}	
		
		oas[qm] = term_data->OAS;		
		
		if (q==0) {									/*first try */
			if (io[0]> deal_data->Inpio)  term_data->OAS += 20.;      /* if calculated price too large, increase oas */
				else term_data->OAS -= 20.;       /* if calculated price too small, decrease oas */
		}
		
		if (q>=1) term_data->OAS = oas[0]+(deal_data->Inpio-io[0])
			*(oas[1]-oas[0])/(io[1]-io[0]);				 /* linear interpolation */			
	} /* end case 14 */

	if (casenumber == 15) {                         /* po  oas analysis */
		if (fabs(po[qm]- deal_data->Inppo) < ACURA/3.0) {
//		 	printf("oas is %lf.\n", term_data->OAS);
//			printf("IO is %lf, PO is %lf, collatoral is %lf.\n", io[qm], po[qm], io[qm]+po[qm]);
			break; 			/* if good enough get out */
		}	
		
		oas[qm] = term_data->OAS;		
		
		if (q==0) {									/*first try */
			if (po[0]> deal_data->Inppo)  term_data->OAS += 20.;      /* if calculated price too large, increase oas */
				else term_data->OAS -= 20.;       /* if calculated price too small, decrease oas */
		}
		
		if (q>=1) term_data->OAS = oas[0]+(deal_data->Inppo-po[0])
			*(oas[1]-oas[0])/(po[1]-po[0]);				 /* linear interpolation */			
	} /* end case 15 */

}/* end of q loop */

	if (casenumber == 1) {
		
//		printf("the gradient of IO with respect to speed up factor is %lf, to steepness factor is %lf. \n",
//			IO_spd, IO_rtn);
//		printf("the gradient of PO with respect to speed up factor is %lf, to steepness factor is %lf. \n",
//			PO_spd, PO_rtn);
//		printf("the gradient of collatoral with respect to speed up factor is %lf \n",
//			IO_spd+PO_spd);
//		printf("the gradient of collatoral with respect to steepness factor is %lf\n",
//			IO_rtn+PO_rtn);
	}


	if (casenumber != 0 && casenumber < 10) {							     /* print out new prepayment */
		ofpios = 0; //fopen("final.out", "w");

		for (i = 0; i < MAXNUMP; i++) {
			for (j = 0; j < NBPOINTS; j++) { 
				if (deal_data->AlorCPR =='C')	deal_data->AmortLevel[i][j] = 100. -100* 
           			 pow (1. - deal_data->AmortLevel[i][j] ,  (double) F);    /* Convert the monthly CPR to annual CPR */
				else deal_data->AmortLevel[i][j]= smmal1(deal_data->AmortLevel[i][j], effmat)+deal_data->NbFixedPrepay/12.;   /* if average life is needed, convert SMM to it */
			}
		} /* for i */
	 

		for (j=0; j<NBPOINTS; j++) {    /* write output */
			for (i=0; i<MAXNUMP; i++) fprintf(ofpios, "  %lf", deal_data->AmortLevel[i][j]);
			fprintf(ofpios, "\n");        /* print out new prepayment */
		}
		fprintf(ofpios, "\n\n");
		for (j=0; j<NBPOINTS; j++) {  /* write old ppf */
//			for (i=0; i<MAXNUMP; i++) fprintf(ofpios, "   %lf", deal_data->AmortLevel1[i][j]);
			fprintf(ofpios, "\n");      /* print out old prepayment */
		}
		fclose(ofpios);
	} /* end if casenumber */
	MBS_Prepay_Engine_Cleanup(&prepay_engine);
	mbs_free_dmatrix (IO, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
	mbs_free_dmatrix (PO, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
	mbs_free_dmatrix (Amort, 0, 2 * deal_data->SFSNb - 1, 0, 2 * ExEnd);
	mbs_free_dmatrix (ParYield, 0, 1, 0,2 * ExEnd);
	mbs_free_dmatrix (Zero, 0, TotNbZero, 0, 2 * TotPer);
	mbs_free_dmatrix (p0, 0,  1, 0, 2*TotPer);
	mbs_free_dmatrix (pd, 0,  1, 0, 2*TotPer);
	mbs_free_dmatrix (pu, 0,  1, 0, 2*TotPer);
	mbs_free_dmatrix (DiscountS, 0,  1, 0, 2*TotPer);
	mbs_free_dmatrix (Discount, 0,  1, 0, 2*TotPer);
	for (i=0; i<= 1; i++) mbs_free_dmatrix (ZeroDelay[i], 0, deal_data->Delay, 0, 2*TotPer);

	deal_data->output[0] = io[qm];
    deal_data->output[1] = po[qm];
	deal_data->output[2] = io[qm]+po[qm];
    deal_data->output[3] = term_data->OAS;

	mbs_free_dvector (io, 0, NTRIALS);
	mbs_free_dvector (po, 0, NTRIALS);


	if (deal_data->Hedge == 'N')	*HedgePrice = 0.;
	if (deal_data->Hedge == 'Y')    *HedgePrice = io[qm]+po[qm];
	if (deal_data->Hedge == 'I')	*HedgePrice = io[qm];
	if (deal_data->Hedge == 'P')    *HedgePrice = po[qm];
	return(0);
}/* impp*/ 

/*****  Tweak  **************************************************************/
/*
*       Tweak an individual deal.
*/
char * Tweak (    double          Option,                                         /* Original option price */
		TERM_DATA       *term_data,                                     /* Structure of term structure data */
		TREE_DATA       *tree_data,                                     /* Structure of tree data */
		DEAL_DATA       *deal_data,                                     /* Structure of deal data */
		PREPAY_DATA     *prepay_data,                                   /* Structure of prepayment data */
		HISTORY         *history,
		MBSDENSITY         *density,
		TRIGGER         *trigger,
		DEAL            *deal,
		PREPAY          *prepay,
		double          *price)
{
	int
		i;
	long
		Date,                                                           /* Maturity of the current money market instrument */
		IndexMM[]   = {0, 1, 2, 3, 6, 12},                              /* Maturity in months of the money market instruments */
		IndexSwap[] = {4, 6, 8, 10, 14, 20, 40, 60};                    /* Maturities in half years of the benchmark swaps */
	double
		F,                                                              /* 'A'nnual or 'S'emi-Annual frequency of the yield curve expressed as a double */
		RateTweak[14],                                                  /* Interest rate tweak of the option */
		Duration[14],                                                   /* Duration of benchmark instruments */
		Time;                                                           /* Number of days / Money market basis (used in tweaking money market rates */
	char * err;

	F = Conv_AoS (term_data->AoS);

	for (i=0; i < 14; i++)                                                  /* Reset values to zero */
		RateTweak[i] = Duration[i] = 0.;


//	printf ("Tweaking %4s rate ...\n", Benchmarks[0]);

	Time = 1. / (double) term_data->MMB;                                    /* Time = days / 360 or days / 365 */

	Duration[0] = 1. / (double) term_data->MMB / pow (1. + term_data->MMYield[0] / 100. / (double) term_data->MMB, 2.);     /* PVBP of O/N instrument */

	term_data->MMYield[0] += .01;                                           /* Tweak the O/N rate by one basis point */

	Zero_6m (term_data);                                                    /* Generate zero rates with tweaked value of rate. */

	if( err=Build_Tree (    term_data,                                              /* Re-build the two dimension tree */
			tree_data,
			prepay_data)) return(nrerror(err));

	if(err=impp (   term_data,                                      /* Re-calculate option price */
				tree_data,
				deal_data,
				history,
				density,
				trigger,
				deal,
				prepay,
				price,
				&(RateTweak[0]))) return(nrerror(err));

	term_data->MMYield[0] -= .01;                                           /* Reset the O/N rate to its original value */

	RateTweak[0] = (RateTweak[0] - Option) / Duration[0] * 10000.;


	i = 1;                                                                  /* Tweaking of the money market rates */
	do
	{
//		printf ("Tweaking %4s rate ...\n", Benchmarks[i]);

		Date = Nxtmth ( term_data->ValueDate,
				IndexMM[i],
				(long) 1);

		Time = Daysact(term_data->ValueDate, Date) / (double) term_data->MMB;   /* Time = days / 360 or days / 365 */

		Duration[i] = Time / pow (1. + term_data->MMYield[i] / 100. * Time, 2.);        /* PVBP of a money market instrument */

		term_data->MMYield[i] += .01;                                   /* Tweak the money market rate by one basis point */

		Zero_6m (term_data);                                            /* Generate zero rates with tweaked value of rate. */

		if(err=Build_Tree (    term_data,                                      /* Re-build the two dimension tree */
				tree_data,
				prepay_data)) return(nrerror(err));

		if( err=impp (   term_data,                                      /* Re-calculate option price */
					tree_data,
					deal_data,
					history,
					density,
					trigger,
					deal,
					prepay,
					price,
					&(RateTweak[i]))) return(nrerror(err));

		term_data->MMYield[i] -= .01;                                   /* Reset the money market rate to its original value */

		RateTweak[i] = (RateTweak[i] - Option) / Duration[i] * 10000.;  /* The Notional of the benchmark instrument we have to hold is determined by the */
										/* sensitivity of the option and the duration of the benchmark instrument.       */
		i++;
										/* We stop if the last benchmark instrument has a maturity greater than the */
	} while ((i < 6) && (Daysact (Date, deal_data->Maturity) > 0));         /* maturity of the underlying as the underlying is not sensitive to the     */
										/* following benchmark instruments.                                         */

	i = 0;                                                                  /* Tweaking of the swap rates */
	while ((i < 8) && (Daysact (Date, deal_data->Maturity) > 0))
	{
//		printf ("Tweaking %4s rate ...\n", Benchmarks[i+6]);

		term_data->SwapYield[i] += .01;                                 /* Tweak the current swap yield by one basis point */
		term_data->YieldS[i] += .01;

		if ( i == 6) deal_data->Refirate += .01;

		Zero_6m (term_data);                                            /* Generate zero rates with tweaked value of rate. */

		if(err=Build_Tree (    term_data,                                      /* Re-build the two dimension tree */
				tree_data,
				prepay_data)) return(nrerror(err));

		if(err=impp (   term_data,                                      /* Re-calculate option price */
					tree_data,
					deal_data,
					history,
					density,
					trigger,
					deal,
					prepay,
					price,
					&(RateTweak[i+6]))) return(nrerror(err));

		Duration[i+6] = 100. / term_data->SwapYield[i]
				* (1. - 1. / pow (1. + term_data->Zero[IndexSwap[i]+3] / F, term_data->Period[IndexSwap[i]+3] * F / 365.));

		term_data->SwapYield[i] -= .01;                                 /* Reset swap yield to its original value */
		term_data->YieldS[i] -= .01;
		if ( i == 6) deal_data->Refirate -= .01;

		RateTweak[i+6] = (RateTweak[i+6] - Option) / Duration[i+6] * 10000.;

		Date = Nxtmth(  term_data->ValueDate,
				IndexSwap[i] * 6,
				(long) 1);

		i++;
	}

//	PrintTweak (    Option,                                                 /* Print tweaks into an Ascii file */
//			RateTweak,
//			Duration,
//			prepay_data);
	return(0);

}  /* Tweak */



/*****  PrintTweak  *********************************************************/
/*
*      Print the tweaking results in a ASCII file.
*/
void PrintTweak (double         Option,                                         /* Option price */
		 double         *RateTweak,                                     /* Interest rate tweak of the option */
		 double         *Duration )                                     /* Duration of benchmark instruments */
{
	int
		i;
	FILE
		*stream;


	stream = 0; //fopen ("TWEAK.prn", "w");

	fprintf (stream,"***************************  Hedge Results  ***********************************\n\n");

	fprintf (stream,"          PVBP       Option\n\n");
	fprintf (stream,"Price           %10.3lf\n\n",Option);

	for(i = 0; i < 14; i++)
		fprintf (stream, "%4s %10.4lf %10.3lf\n", Benchmarks[i], Duration[i], RateTweak[i]);


	fclose  (stream);

	return;

}  /* PrintTweak */

/*****  Main_Prepay  ********************************************************/
/*
*       Prepayment function (SFS specific). It includes schedu-
*       led amortization & seasonality (time dependent) as well as seasoning
*       & index amortization (index dependent).
*		It is insured that it lies in te range [0,1]
*/
void    Main_Prepay (   double   ** Amort,                                 /* Array of smms of various speed groups in the interest rate lattice */
			double          *ParYield,                              /* Par Yield in the lattice at the current period */
			int             CouponNb,                               /* Coupon Number: = deal_data->Term at the end of the tree */
			long            Month,                                  /* Current month of the year (1 to 12) */
			int             NbPer,									/* Current time period */
			int             N,                           
            MBS_Prepay_Engine *	prepay_engine,
			TREE_DATA       *tree_data                             /* Structure of tree data */
			)
{
	double 
		*PAmort;
	double *smm; //smm at current period with trigger
	int
		Bottom,                                                         /* Bottom node index in the interest rate tree */
		Top,                                                            /* Top node index in the interest rate tree */
		i, j;

	PAmort = mbs_dvector(0, NBPOINTS); 
	smm = mbs_dvector(0,prepay_engine->SFSNb - 1 );
/* assign NBPOINTS spaces to PAmort */	

	Bottom  = (int) max (0., NbPer - tree_data->NodeMax);                    /* If NbPer>NodeMax the interest rate tree is cut at the top and the bottom */
	Top     = (int) min (2*NbPer, NbPer + tree_data->NodeMax);


/* Use the right prepayment function at each period */

	for (i = Bottom; i <= Top; i ++)
	{

		MBS_Prepay_Engine_SMM( smm, prepay_engine,
			CouponNb/*paymentNb*/, ParYield[i],
			Month, NbPer - N );
		for(j=0; j < prepay_engine->SFSNb; ++j)
			Amort[j][i] = smm[j];
	}  /* for i */

	mbs_free_dvector(PAmort,0,NBPOINTS);
	mbs_free_dvector(smm,0,prepay_engine->SFSNb-1);
	return;

}  // Main_Prepay


//if it is plain vanilla PT and input S curve is quoted in CPR
//this is not needed, the work is done in MBS_Prepay_Engine_SetParams
int TweakSCurve( DEAL_DATA *deal_data, long today )//1 means succcessful
{
	long matday, wam, effmat;
	int i, j;
	double maxSeasonality, F;
	matday = deal_data->Exer[1];
	today /= 100;
	matday /=100;
	wam = (matday/100 - today/100)*12 + (-today+matday)             //wam
		-(matday/100 - today/100)*100;                           //WAM minus lockout
	effmat = wam - deal_data->NbFixedPrepay;
	F = Conv_Freq (deal_data->Freq);

	if (deal_data->AlorCPR == 'C'){ //if input is CPR
		for (i = 0; i < MAXNUMP; i++) {
			for (j = 0; j < NBPOINTS; j++)  
				deal_data->AmortLevel[i][j] = 1. - 
           			 pow (1. - deal_data->AmortLevel1[i][j] / 100., 1. / (double) F);    /* Convert the annual CPR to SMM */
			} //for i
		for (i=0; i < MAXNUMP; i++) 
			spdbase(deal_data->AmortLevel[i],deal_data->Speedup,deal_data->rational,NBPOINTS); //speed up and steepen
	} // end if 
	else if  (deal_data->Schedule == 'N') {  // if input average life and no scheduled amortization, suitable for busted PACs
		for (i=0; i < MAXNUMP; i++) {
			for (j=0; j < NBPOINTS; j++) 
				deal_data->AmortLevel[i][j] = alsmm1(12.*deal_data->AmortLevel1[i][j]-deal_data->NbFixedPrepay, effmat);  //turn average life to SMM
		}	
		for (i=0; i<MAXNUMP; i++) 
				spdbase(deal_data->AmortLevel[i], deal_data->Speedup, deal_data->rational, NBPOINTS); /*speed up and steepen */
			
	  	for (j=0; j < 2* deal_data->SFSNb; j++) deal_data->Amort[j] = alsmm1(12.*deal_data->Amort[j]-deal_data->NbFixedPrepay,effmat);	
	}
	else if (deal_data->Balloon == 'Y')  { /*if avg life and sch amort and balloon, we left out non-balloon,if there is such a thing, set balloon schedule to original term*/
		for (i=0; i<MAXNUMP; i++) {
			for (j=0; j < NBPOINTS; j++) 
				deal_data->AmortLevel[i][j] = alsmm2(12.*deal_data->AmortLevel1[i][j]-deal_data->NbFixedPrepay,deal_data->Gwac, deal_data->BalloonSch,effmat, deal_data->Term - effmat);
		}
		for (i=0; i<MAXNUMP; i++)
				spdbase(deal_data->AmortLevel[i], deal_data->Speedup, deal_data->rational, NBPOINTS); /* speed up and steepen */
		for (j=0; j<2*deal_data->SFSNb; j++) deal_data->Amort[j] = alsmm2(12. *deal_data->Amort[j]-deal_data->NbFixedPrepay, deal_data->Gwac, deal_data->BalloonSch, effmat, deal_data->Term-effmat); 
	}
//check:
	maxSeasonality = deal_data->Seasonality[0];
	for(i=1; i < 12; ++i)
		maxSeasonality = max( maxSeasonality, deal_data->Seasonality[i] );
	for (i=0; i<MAXNUMP; i++)
		for (j=0; j < NBPOINTS; j++) 
		{
			double level = deal_data->AmortLevel[i][j] * maxSeasonality;
			if(level > 1.0 || level < 0.0 ) return 0;
		}
	return 1;
}//TweakSCurve
