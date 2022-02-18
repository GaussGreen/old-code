#ifndef SRT_H_ALL_H
#define SRT_H_ALL_H

/* PRINT TO DEBUG: NOT USED... */
#define PR
#define PR2
#define PS
#define PRS
/*
#include "math.h"
*/
#include "ctype.h"
#include "float.h"
#include "limits.h"
#include "setjmp.h"
#include "signal.h"
#include "stdarg.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"

/* SRT UTILITY LIBRARY */
#include "utallhdr.h"

/* SRT NUMERICAL ROUTINES LIBRARY */
#include "num_h_allhdr.h"

/* SRT SWAPS LIBRARY */
#include "srt_h_ts.h"
#include "swp_h_all.h"

/* These things are more or less utilities as well... */
#include "srt_h_types.h"

/* MARKET LIST */
/* #include "srt_h_fwdobj.h" */
#include "srt_h_bnd_aux.h"
#include "srt_h_und_fct.h"
#include "srt_h_und_list.h"
#include "srt_h_und_struct.h"

/* Pre-initialized products routines */
#include "srt_h_product_list.h"

/* SRT LIST AND TERM STRUCT */
#include "srt_h_ts.h"
#include "srt_h_ts_eq.h"
#include "srt_h_ts_fx.h"
#include "srt_h_ts_irm.h"
#include "srt_h_ts_irm_2f.h"
/* CORRELATION LIST */
#include "srt_h_correlation_list.h"

/* INTEREST RATE MODELS */
#include "srt_h_alltrefct.h"
#include "srt_h_autocal_future.h"
#include "srt_h_betaetaM.h"
#include "srt_h_bnd_aux.h"
#include "srt_h_grfn.h"
#include "srt_h_grfn_event.h"
#include "srt_h_grfn_undinfo.h"
#include "srt_h_grfnparams.h"
#include "srt_h_iostruct.h"
#include "srt_h_irmfct.h"
#include "srt_h_lgm2dtre.h"
#include "srt_h_lgm2dtreefct.h"
#include "srt_h_lgmclsdfrm.h"
#include "srt_h_lgmnmint.h"
#include "srt_h_make_time_step.h"
#include "srt_h_mc_core.h"
#include "srt_h_mc_genrand.h"
#include "srt_h_mdl_types.h"
#include "srt_h_mdlinistp.h"
#include "srt_h_readrequest.h"
#include "srt_h_sample.h"
#include "srt_h_step_list.h"
#include "srt_h_tminfstruct.h"
#include "srt_h_trebdt.h"
#include "srt_h_tree_copy.h"
#include "srt_h_tree_core.h"
#include "srt_h_trestruct.h"
#include "srt_h_twoundtre.h"
#include "srt_h_und_struct.h"
#include "srt_h_vegatreche2dr.h"
#include "srt_h_vegatrelgm.h"
#include "srt_h_vegatrelog.h"
#include "srt_h_ytat.h"

/*CALIBRATION*/
#include "srt_h_calib.h"
#include "srt_h_calibparams.h"
#include "srt_h_calibutils.h"
#include "srt_h_twofaccorrel.h"

/* FX OPTIONS */
#include "srt_h_quanto_fct.h"

/* RESET CAPS */
#include "srt_h_fwdvol.h"
#include "srt_h_fwdvolfnc.h"
#include "srt_h_lgmresetcap.h"

/*	Chey Beta PDE	*/
#include "srt_h_cheybeta_bck_pde.h"
#include "srt_h_cheybeta_fwd_pde.h"
#include "srt_h_cheybeta_fwd_pricing.h"
#include "srt_h_cheybeta_grid.h"
#include "srt_h_cheybeta_mdl_def.h"

#endif
