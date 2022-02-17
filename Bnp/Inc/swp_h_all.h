#ifndef SWP_H_ALL_H
#define SWP_H_ALL_H

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

/* Some utilities: memory allocation  , errors  ,.. */
#include "utallhdr.h"

/* For external functions */
#include "swp_h_external_fct.h"

/* These things are more or less utilities as well... */
#include "swp_h_gen_err_mesg.h"
#include "swp_h_string_interp.h"

/* Currency Parameters */
#include "swp_h_ccy_param.h"

/* For Dates and Schedules generations */
#include "swp_h_swapdp.h"

/* Curves Objects */
#include "swp_h_ycinstr.h"
#include "swp_h_zc_types.h"
/* #include "srt_h_fwdobj.h" */
#include "swp_h_cmt_param.h"
#include "swp_h_cmt_types.h"

/* The Curves and the List */
#include "swp_h_curve_struct.h"

/* Uses SwapDateDir */
#include "swp_h_swap_compute.h"
#include "swp_h_swap_generic.h"
#include "swp_h_swap_pricing.h"

#include "swp_h_datelist.h"
#include "swp_h_irrrng.h"
#include "swp_h_swap_simple.h"
#include "swp_h_swap_utils.h"

/* ZERO CURVE */

/* RATES THINGS */
#include "srt_h_rtdesc.h"
#include "srt_h_rtdescevl.h"

/* Df things */
#include "swp_h_df.h"

/* Spread things */
#include "swp_h_spread.h"

/* SWAP TOOLS */
#include "swp_h_capfloor.h"
#include "swp_h_df.h"
#include "swp_h_gen_deal.h"
#include "swp_h_swaption.h"

/* BOND THINGS */
#include "swp_h_bndcrv.h"
#include "swp_h_bond_compute.h"
#include "swp_h_bond_cox_compute.h"
#include "swp_h_bond_trin_compute.h"
#include "swp_h_bond_vol.h"

/* CMT CMS */
#include "swp_h_cmt_capfloor.h"
#include "swp_h_cmt_fwdtreas.h"
#include "swp_h_cmt_init.h"
#include "swp_h_cmt_swap.h"
#include "swp_h_cmt_types.h"

/* CURVE LIST */
#include "swp_h_curve_fct.h"
#include "swp_h_curve_list.h"

/* ========================================================================== */
#endif
