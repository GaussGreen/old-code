/*********************************************

        SWP_H_BOND_TRIN_COMPUTE.H

authors: A.P.      , O.VE.      , A.B.      , HdL
version: 24/10/94

Last modified: Dec 2      , 1994 by K L Chau
Reason: use new SrtMkt structure
**********************************************/

#define MAXI 10000

#define EUROPEAN 0
#define AMERICAN 1
#define MID_ATLANTIC 2
#define SHOUT_LOCK 3
#define SHOUT_RESET 4
#define BOTH_SHOUT 5
#define YOU_CHOOSE 6
#define INSTALMENT 7
#define EXTINGUISH 8
#define EURO_DIGIT 9
#define AMER_DIGIT 10
#define PAY_DIGIT 11

#include "srt_h_ts.h"

typedef struct model_struct /* idem forex tree	*/
{
  int step_num;   /* Number of steps in the tree */
  int size_tsv;   /* Size of vol_vector 					*/
  int mat_vol_nb; /* NB associated to opt mat in vol range	*/
  double **tsv;   /* Volatility structure : dates + vol 			*/
  TermStruct *ts;
  double **date_cumvol; /* Cumulative_volatility structure : dates + cum_vol
                         */
  double cum_vol_opt;   /* Cum_vol for the whole option maturity 		*/
  String repo;          /* Yield curve : REPO 					*/
  String swap;          /* Yield Curve : SWAP					*/
  double df_swap_all;   /* Discount factor for the entire option maturity	*/
  double fwd_opt_mat;   /* fwd bond price at option maturity	*/
  double delta_lns;     /* discretisation step applied to the log of the spot
                         */
  double h;             /* Time step used in the tree				*/
  double diff_up;       /* Sup-Smid						*/
  double diff_do;       /* Smid-Sdown						*/
  int hit_level;        /* extremum reached 					*/
} Model;

typedef struct tree_step_struct /* diff from forex tree	*/
{
  double expected_spot; /* Expected value of the underlying */
  double df_swap;       /* Discount factor from step to end of option maturity
                         */
  double cum_vol;       /* Cumulative vol from step to end of option maturity	*/
  double df_swap_end; /* Discount factor from step to end of tree = sub_period*/
  double df_swap_h;   /* Discount factor from step to step+1		*/
  double drift_h;     /* Drift for exp. spot from step to step+1		*/
  double vol_h;       /* Cumulative vol from step to step + 1			*/
  double vol;         /* Vol from step to option maturity (used in B&S)	*/
  double vol_inst;    /* Vol from step to step+1 (used in B&S)	*/
  double time_end;    /* Time from step to end of tree			*/
  double time_start;  /* Time from start to step				*/
  double p_u;         /* probability to go up					*/
  double p_d;         /* Probability to go middle				*/
  double p_m;         /* Probability to go down				*/
  int d_index;        /* */
  int flag_date;      /* */
  double cpn;         /* FOR INSTALMENT*/
  double clean_strike;  /* strike is constant in yield -- a price */
  double clean_barrier; /* barrier is constant in yield -- a price */
  int yield_hit_level;  /* hit level will depend on yield */
} Step;

typedef struct node_struct {
  double intrinsic;    /* (S-K) for a call      , (K-S) for a put	*/
  double asset;        /* Underlying value			*/
  double disc_premium; /* Discounted expected value of the option  */
  double forward;      /* From NODE to OPTION MATURITY		*/
  int position;
  double p_u; /*  for extinguish			*/
  double p_m; /*  for extinguish			*/
  double p_d; /*  for extinguish			*/
} Node;

typedef struct deal_struct /* idem forex tree			*/
{
  int type;
  int call_put;
  double spot;
  Date today;
  Date sub_period; /* end of the tree (for partial options) */
  Date opt_mat;
  double strike;
  double coupon;       /* FOR THE BOND 			*/
  double first_coupon; /* FOR THE BOND 			*/
  SwapDP p;
  int d_set;
  int ret_type;           /* Premium (0) or Greeks ( 1 or 2 ) 	*/
  int size_date_cpn;      /* FOR INSTALMENT & MID			*/
  Date *date;             /* FOR INSTALMENT & MID			*/
  double *cpn;            /* FOR INSTALMENT			*/
  double barrier;         /* used for extinguish...		*/
  int do_up;              /* down or up extinguish...	  	*/
  double corr_fx;         /* QUANTO correlationbetween domestic and fx */
  int quanto;             /* QUANTO flag 0 no 1 yes */
  double yield_strike;    /* strike is a yield      */
  double yield_barrier;   /* barrier is a also a yield*/
  int yield_strike_flag;  /* flag for strikes yield only 0 no 1 yes */
  int yield_barrier_flag; /* flag for barrier yields  0 no 1 yes */
} Deal;

/************* PLEASE REMEMBER WE ALREADY HAVE THE FOLLOWING STRUCUTRE *******
typedef struct {
        String 		bond_name;
        SwapDP		p;
        double 		coupon;
        int		settle_days;
        double 		clean;
        double 		yield;
        char		swap_name[SRTBUFSZ]	;
        char		repo_name[SRTBUFSZ]	;
        TermStruct	*ts;
        SrtBndAuxParam 	aux_par;
} SrtBndDesc;
******************************************************************************/

typedef struct func_struct /* idem forex tree	*/
{
  double (*PAY_OFF)(Deal, Model, Step, Node);  /* in the middle of the tree
                                                */
  double (*END_TREE)(Deal, Model, Step, Node); /* at the end of tree  */
} Func;

double bond_trinomial_fct(Deal contract, Model model, Model model_fx,
                          double *vega);
