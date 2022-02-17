#ifndef SRT_H_GRFNPARAMS_H
#define SRT_H_GRFNPARAMS_H

#include "grf_h_types.h"

/* ------------------------------------------------------------------------
        STRUCT:		SrtGrfnParam

        PURPOSE:	Store all information required for a Grfn pricing

        CONTENTS:

        ------------------------------------------------------------------------
        These can be user defined from the spreadsheet and are important
   ------------------------------------------------------------------------
                long 		num_MCarlo_paths;
                                Number of paths used in Monte Carlo
                double 		max_time_per_slice;
                                Maximum time in years authorised between two
time steps . long 		min_nodes_per_path; Minimum number of time steps
to create for discretisation . SRT_Boolean 	force_mc; Force the use of Monte
Carlo for a "neutral" tableau . SRT_Boolean     Jumping; Force the use if
possible of the Jumping Numeraire method (MC Core). SrtMCSamType
sample_type; Monte Carlo random number generation method (ABS  , Sobol  , ...) .
                SrtMCDfSchType 	difference_scheme;
                                Monte Carlo discretisation scheme (Euler or
Milshtein) . SRT_Boolean 	mc_renormalize; Force the Cheyette type Monte
Carlo control variate on dfs .

- ------------------------------------------------------------------------
        These can be user defined from the spreadsheet and can be useful
   ------------------------------------------------------------------------ .
                long        rand_seed;
                                Initial seed for the random number generator.
                long 		width_phi;
                                Number of discrete phi used in a tree
discretisation of Cheyette . SRT_Boolean		imm_exit; Prevents the
calculation of the grfn function . SRT_Boolean		skip_evaluation; Skip
the evaluation of the grfn deal: initialisation is done . SRT_Boolean
insert_vol_dates; Force the time discretisation to incorporate the volatility
                                dates .

                SRT_Boolean		end_of_day_flg;
                                If true  , today's line is conisdered as being
                                past and payments occurring today will not be
                                taken into account
                                and fixings will have of be used
                SRT_Boolean		end_of_day_fixings;
                                If true  , today's line is conisdered as being
                                past and fixings will be required for events
                                occuring today.
   ------------------------------------------------------------------------
        These can be user defined from the spreadsheet and are for debug
   ------------------------------------------------------------------------ .
                SRT_Boolean 	prin_steps;
                                To print steps information in a file .
                SRT_Boolean 	prin_deal;
                                To print deal information in a file .
                SRT_Boolean 	prin_ok;
                                To print a cell from a grfn tableau in a file .

   ------------------------------------------------------------------------
        These can be user defined from the spreadsheet and are USELESS :
        they have no effect at all  , but are there for backward compatibility
   ------------------------------------------------------------------------ .
                SRT_Boolean		spline ;
                                Only for backward compatibility : not used any
longer . SRT_Boolean 	dbg; Only for backward compatibility : not used any
longer . long 		max_nodes_per_path; Only for backward compatibility :
not used any longer . SRT_Boolean 	insert_tenor_dates; Only for backward
compatibility : not used any longer .

   ------------------------------------------------------------------------
        These parameters are important  , but are set up by the program
   ------------------------------------------------------------------------ .
                long 		step_num;
                                Number of steps in a tree discretisation: set
automatically . SrtMdlType		imp_mdl_type; New field added for
backward compatibility : the model used . SrtImpType 	imp_type; This will be
set to either IMPMONTECARLO or IMPBACKWARDLATTICE . int node_dim; Dimension of a
tree node according to the length of the grfn tableau .

        -----------------------------------------------------------------------
*/

typedef struct {
  long num_MCarlo_paths;
  double max_time_per_slice;
  long min_nodes_per_path;
  long min_pde_num_mesh;
  SRT_Boolean force_mc;
  SRT_Boolean jumping;
  SrtMCSamType sample_type;
  SrtMCDfSchType difference_scheme;
  SRT_Boolean mc_renormalize;

  long rand_seed;
  long width_phi;
  SRT_Boolean imm_exit;
  SRT_Boolean skip_evaluation;
  SRT_Boolean insert_vol_dates;

  SRT_Boolean end_of_day_flg;
  SRT_Boolean end_of_day_fixings;
  SRT_Boolean end_of_day_payment;

  SRT_Boolean prin_steps;
  SRT_Boolean prin_deal;
  SRT_Boolean prin_ok;

  SRT_Boolean spline;
  SRT_Boolean dbg;
  long max_nodes_per_path;
  SRT_Boolean insert_tenor_dates;

  long step_num;
  SrtMdlType imp_mdl_type;
  SrtImpType imp_type;
  int node_dim;

  SrtClsdFrmType closed_form_type;

  SRT_Boolean calib;
  SRT_Boolean first_time;
  SRT_Boolean exfrontierfirst_time; /*parameters needed for the exercise
                                       frontier optimization*/
  SRT_Boolean exfrontierlast_time;  /*""*/
  int exfrontiercounter;            /*""*/
  SrtMeasureType measure;

  Forback lsm;
  SRT_Boolean minim;
  long colmininf;
  long colminsup;
  long colmincible;
  double minmaxtime;
  long minfreedom;
  int exfrontier;
  int exfrontiernumpoints;
  double gammaexfrontier;
  SrtReceiverType recpay;
  int MaxNumDisc;
  int MinNumDisc;

} SrtGrfnParam;

/* ------------------------------------------------------------------------
        STRUCT:		SrtMdl  , SrtMdlPtr

        PURPOSE:	Store external information required for a Grfn pricing
   ------------------------------------------------------------------------ */

typedef struct {

  SrtFncWrapper declare_fixing_wrapper; /* declare cashflow that has been fixed
                                   but has not taken place to caller */
  SrtFncWrapper historical; /* get historical market data from caller */

} SrtMdl, *SrtMdlPtr;

/* ------------------------------------------------------------------------
        Functions to initialise or fill in the structure...
   ------------------------------------------------------------------------ */

Err srt_f_set_GrfnParams(int numParams, String *paramStrings,
                         String *valueStrings, SrtGrfnParam *param);

Err srt_f_set_default_GrfnParams(SrtGrfnParam *param);

#endif
