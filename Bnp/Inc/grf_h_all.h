#ifndef GRF_H_ALL_H
#define GRF_H_ALL_H

/*<%%STA>-----------------------------------------------------------------
  DESCRIPTION     :Grfn include file with all function declarations that are
not considered public, and which do not operate on the COMLL structure. (Those
are in grf_h_lang.h; public functions are declared in grf_h_public.h or one of
the files it includes).
This file #includes grf_h_lex.h;

<%%END>---------------------------------------------------------------------*/

#include "grf_h_lang.h"
#include "grf_h_lex.h"
#include "grf_h_public.h"
#include "grf_h_types.h"
#include "num_h_allhdr.h"
#include "srt_h_correlation_list.h"
#include "srt_h_grfn.h"
#include "srt_h_grfn_undinfo.h"
#include "srt_h_grfnparams.h"
#include "srt_h_sample.h"
#include "srt_h_step_list.h"
#include "srt_h_ts.h"
#include "srt_h_types.h"
#include "srt_h_und_struct.h"
#include "swp_h_all.h"
#include "utallhdr.h"

/*	Added Toto 29Nov1999	*/
#include "grf_h_mdlcomm.h"

/*============grfn internal structures ==========================*/

/** grf_f_types.c **/
Err grfn_check_tableau_and_dates(long* nrows, long ncols, GrfnCell** sprdsht, Date* eventdates);

void grfn_free_GrfnEvent(GrfnEvent* event);
void grfn_free_inGrfnDeal(GrfnDeal* gd);
Err  grfn_copy_to_GrfnDeal(
     GrfnDeal*  gd,
     long       nrows,
     long       ncols,
     GrfnCell** sprdsht,
     Date*      eventdates,
     long       numgrng,
     GrfnRng*   grng,
     long       auxwidth,
     long*      auxlen,
     double**   aux);

GrfnCell** GrfnCellmatrix(long nrow, long ncol, long slen);

void grfn_free_GrfnCellmatrix(GrfnCell** m, long nrow, long ncol);

Err grfn_extend_tableau(GrfnCell*** tab, long* nrow, long ncol);

/*======================================*/

Err srt_f_set_GrfnCell(
    int tabRows, int tabCols, char*** tabStrings, int** tabMask, GrfnCell*** tableau);

void free_the_world(SrtStpPtr sptr, GrfnDeal* gd, SrtUndInfo* und_info);

Err attach_params_to_GrfnEvent(GrfnDeal* gd, SrtStpPtr sptr, SrtUndPtr und, int und_index);

Err get_deal_und_info(GrfnDeal* gd, SrtUndInfo* und_info);

int first_unkn_date_index(GrfnDeal* gd, SRT_Boolean end_of_day_flag);

/*=========grfn language processing =============================*/

/* grf_f_lang.c */
Err grfn_interp_comll_name(String name, GrfnSymbol* gs);
Err grfn_interp_name(COMLL_PTR c, String name, GrfnDeal* gd);
Err grfn_interp_func(String name, COMLL_PTR args, GrfnDeal* gd);
Err grfn_interp_cell_func(String name, COMLL_PTR args, GrfnDeal* gd);
Err grfn_check_ranges(GrfnDeal* gd);

EXTERN SRT_Boolean grfn_lang_cur_amcell;

/* grf_f_evlcllfct.c */
void grf_f_evlcllfct(
    double** cells, int c1, int r1, int r2, int ind, COMLLType type, double* answer);

void grf_f_evlcllfct_interp(
    double*   x_values,
    double**  cells,
    int       c1,
    int       r1,
    int       r2,
    double    x,
    COMLLType type,
    double*   answer);

/* grf_f_lang_finance.c */
Err grfn_interp_fin_func(String name, GrfnSymbol* gs, GrfnDeal* gd, COMLL_PTR lastarg);

/* grf_f_underlying.c */
Err  grfn_list_deal_underlyings(char** list, GrfnDeal* gd);
void grfn_free_und_name(void);
int  grfn_store_und_name(String x, String y);

int    get_grf_und_index(String c);
String get_grf_und_name(int index);

/* grf_f_lang_pvrng.c */
Err grfn_pvrange_setup(String name, GrfnDeal* gd, COMLL_PTR top, int index);
Err grfn_dirtyprice_setup(String name, GrfnDeal* gd, COMLL_PTR top, int index);

/* grfn_event_create.c */
Err grfn_set_df_dates_in_event(GrfnEvent* nev, GrfnDeal* gd, int index);

/* grf_f_get_event.c */
Err grfn_create_future_event(
    GrfnDeal* gd, GrfnEvent** event, double step_date, double next_step_date);

Err grfn_create_historical_event(
    GrfnDeal* gd, GrfnEvent** event, double step_date, double next_step_date);

Err grfn_attach_past_info_to_future(GrfnDeal* gd, GrfnEvent** today_event);

/*
Err grfn_attach_events(GrfnDeal *gd, SrtStpPtr sptr, SRT_Boolean flag);
*/

/*======evaluate================================================*/

/* grf_f_eval_event.c */
Err grfn_eval_event(
    GrfnEvent*      event,    /* Grfn Event as attached to a step   */
    SrtSample*      sample,   /* State variables                    */
    GrfnDeal*       gd,       /* GRFN Workspace                     */
    GrfnLatticeVec  lattice,  /* Discounted Vector in tree          */
    EvalEventDfsFct evaldf,   /* Function for calculation of dfs    */
    SrtUndInfo*     und_info, /* Underlying Information             */
    double*         cashflow);        /* CashFlow calculated via the comll  */

Err grfn_eval_event_ammc(
    int             type_eval, /* 0: fwd, 1: bwd, 2: nostatus */
    GrfnEvent*      event,
    SrtSample*      sample,
    GrfnDeal*       gd,
    GrfnLatticeVec  lattice,
    double*         fwd,
    double*         cur,
    EvalEventDfsFct evaldf,
    SrtUndInfo*     und_info,
    double*         cashflow);

/*
Err grfn_eval_event_ammc2
    (    int               type_eval,
                GrfnEvent       *event ,
                SrtSample       *sample,
        GrfnDeal        *gd,
        GrfnLatticeVec   lattice,
                double           *fwd,
                double           *cur,
        EvalEventDfsFct       evaldf,
        SrtUndInfo      *und_info,
                double			*cashflow );
*/

/*===== read to and from a flat file ========*/

Err grf_f_readgrfdl(FILE* in, GrfnDeal* gd);
Err grf_f_pringrfdl(FILE* out, GrfnDeal* gd);
Err grf_f_prgrdl(
    String     filename,
    int        numeventdates,
    Date*      eventdates,
    long       nrows,
    long       ncols,
    GrfnCell** sprdsht,
    long       numgrng,
    GrfnRng*   grng,
    long       auxwidth,
    long*      auxlen,
    double**   aux);

#include "grf_h_fltfl.h"

/*===============================================================*/

#endif
