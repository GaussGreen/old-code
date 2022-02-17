/*
        grf_f_mdlcomm.cxx
        ===============

        Author: Toto 26Nov1999

        Objective: define communications between any model and GRFN
*/

#include "grf_h_all.h"
#include "srt_h_all.h"

#define MAXSTRING 256

/*	Copy a spreadsheet cell	*/
static void copy_grfn_cell(GrfnCell *to, GrfnCell *from) {
  to->type = from->type;
  GrfnCopyStatus(to->status, from->status);
  to->str_alloced = from->str_alloced;
  to->dval = from->dval;
  to->sval = NULL;

  if (from->type != GRFNBCELL && from->sval != NULL &&
      from->str_alloced == SRT_YES) {
    to->sval = (char *)calloc(MAXSTRING, sizeof(char));
    strcpy(to->sval, from->sval);
  }
}

/*	Insert a blank event date in the spreadsheet table	*/
static void insert_blank_evt(long date, int insert_before, long *tab_rw,
                             long tab_cl, long **evt_dts, GrfnCell ***tab) {
  int i, j;
  GrfnCell **new_tab;
  long *new_evt_dts;

  new_evt_dts = (long *)calloc(*tab_rw + 1, sizeof(long));
  new_tab = GrfnCellmatrix(*tab_rw + 1, tab_cl, 0);

  if (insert_before < *tab_rw) {
    for (i = *tab_rw - 1; i >= insert_before; i--) {
      new_evt_dts[i + 1] = (*evt_dts)[i];
      for (j = 0; j < tab_cl; j++) {
        copy_grfn_cell(&(new_tab[i + 1][j]), &((*tab)[i][j]));
      }
    }
  }

  if (insert_before > 0) {
    for (i = 0; i < insert_before; i++) {
      new_evt_dts[i] = (*evt_dts)[i];
      for (j = 0; j < tab_cl; j++) {
        copy_grfn_cell(&(new_tab[i][j]), &((*tab)[i][j]));
      }
    }
  }

  new_evt_dts[insert_before] = date;
  for (j = 0; j < tab_cl; j++) {
    new_tab[insert_before][j].type = GRFNDCELL;
    new_tab[insert_before][j].dval = 0.0;
    new_tab[insert_before][j].sval = NULL;
    new_tab[insert_before][j].str_alloced = SRT_NO;
  }

  (*tab_rw)++;
  free(*evt_dts);
  grfn_free_GrfnCellmatrix(*tab, *tab_rw - 1, tab_cl);
  *evt_dts = new_evt_dts;
  *tab = new_tab;
}

/*	The function to pass to grfn_eval_event to calculate
        deterministic DF for historical valuation	*/
static void calc_det_df(GrfnEvent *event, SrtSample *sample,
                        SrtUndInfo *und_info) {
  SrtUndPtr und;
  long ref;
  char *yc_name;
  int i, j;

  /*	Loop through the underlyings in UndInfo one by one
          and get their yeild curves	*/
  for (i = 0; i < und_info->no_of_underlyings; i++) {
    /*	Get the Underlying	*/
    und = lookup_und(und_info->und_data[i].und_name);

    /*	Get its reference date	*/
    ref = get_today_from_underlying(und);

    /*	Get the curve name	*/
    get_underlying_discname(und, &yc_name);

    /*	Compute discount factors	*/
    for (j = 0; j < event->dflen[i]; j++) {
      event->df[i][j] = swp_f_df(ref, event->dfd[i][j], yc_name);
    }
  }
}

/*
        Init function
        =============

        Inputs:	GRFN tableau
                        Domestic underlying name
                        GRFN param

        Output:	A pointer on a FIRSTAllMkts structure        ,
                        fully allocated        ,
                        with requests on market parameters (df maturities)

        Initialise        , check and process GRFN tableau
        Value historical events
        Choose method
*/

Err FIRSTInitMktStruct(
    /*	GRFN Tableau	*/
    int num_evt_dts, long *evt_dts, long tab_rw, long tab_cl, char ***tab_str,
    int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	GRFN choice: backward (-1) or forward (+1)	*/
    int *back_or_for,
    /*	Output	*/
    GRFNCOMMSTRUCT comm) {

  GrfnDeal *gd;
  SrtUndInfo *und_info;
  Err err;
  SrtUndPtr und;
  long ref, first_stop;
  int loc_tab_rw, loc_tab_cl;
  long *loc_evt_dts;
  GrfnCell **tab;
  GrfnEvent **event;
  int dte_idx, dte_idx2, first_stop_idx;
  char *dom_crv_nme;
  double cash_flow;
  int i;
  int front_evt;

  /*	Diverse initialisations	*/

  /*	Allocate global structures	*/
  gd = (GrfnDeal *)malloc(sizeof(GrfnDeal));
  und_info = (SrtUndInfo *)malloc(sizeof(SrtUndInfo));

  /*	Default err to NULL	*/
  err = NULL;

  /*	Check dimensionality	*/
  if (num_evt_dts != tab_rw) {
    free(gd);
    free(und_info);

    return serror("Number of event dates (%d) is inconsistent "
                  "with number of rows (%d) in GRFN table",
                  num_evt_dts, tab_rw);
  }

  /*	Get underlying from market list	*/

  /*	Get the domestic underlying	*/
  if (!(und = lookup_und(dom_nme))) {
    free(gd);
    free(und_info);
    return serror("Could not find % s underlying in market list", dom_nme);
  }
  /*	Get the reference date	*/
  ref = get_today_from_underlying(und);
  first_stop = ref + (prm->end_of_day_flg == SRT_YES ? 1 : 0);

  /*	Stage 1 - Process tableau in the GrfnCell structure		*/
  /*	===================================================		*/

  /*	Copy event dates to a local array	*/
  loc_tab_rw = num_evt_dts;
  loc_tab_cl = tab_cl;
  loc_evt_dts = (long *)calloc(num_evt_dts, sizeof(long));
  memcpy(loc_evt_dts, evt_dts, num_evt_dts * sizeof(long));

  /*	Copy table to GrfnCell structure	*/
  if (err = srt_f_set_GrfnCell(loc_tab_rw, loc_tab_cl, tab_str, mask, &tab)) {
    free(gd);
    free(und_info);
    free(loc_evt_dts);
    return err;
  }

  /*	Remove empty rows        , and check that dates are increasing */
  if (err = grfn_check_tableau_and_dates(&loc_tab_rw, loc_tab_cl, tab,
                                         loc_evt_dts)) {
    free(gd);
    free(und_info);
    free(loc_evt_dts);
    grfn_free_GrfnCellmatrix(tab, tab_rw, loc_tab_cl);
    return err;
  }

  /*	Insert first stop	*/
  if (loc_evt_dts[0] > first_stop) {
    /*	Insert front	*/
    /*
                    insert_blank_evt (	first_stop        ,
                                                            0        ,
                                                            &loc_tab_rw        ,
                                                            loc_tab_cl        ,
                                                            &loc_evt_dts ,
                                                            &tab);
    */
    first_stop_idx = 0;
  } else if (loc_evt_dts[loc_tab_rw - 1] < first_stop) {
    /*	Insert back	*/
    insert_blank_evt(first_stop, loc_tab_rw, &loc_tab_rw, loc_tab_cl,
                     &loc_evt_dts, &tab);
    first_stop_idx = loc_tab_rw - 1;
  } else if (loc_evt_dts[loc_tab_rw - 1] == first_stop) {
    /*	Insert blank event after today */
    insert_blank_evt(first_stop + 1, loc_tab_rw, &loc_tab_rw, loc_tab_cl,
                     &loc_evt_dts, &tab);
    first_stop_idx = loc_tab_rw - 2;
  } else {
    dte_idx = 0;
    while (loc_evt_dts[dte_idx] < first_stop) {
      dte_idx++;
    }

    /*	Today is not there: insert in the right place	*/
    /*
                    if (loc_evt_dts[dte_idx] != first_stop)
                    {
                            insert_blank_evt (	first_stop        ,
                                                                    dte_idx ,
                                                                    &loc_tab_rw
             , loc_tab_cl        , &loc_evt_dts        , &tab);
                    }
    */
    first_stop_idx = dte_idx;
  }

  /*	Stage 2 - Process tableau in the GrfnDeal structure		*/
  /*	===================================================		*/

  /*	Make copies of all the inputs and place them inside the GrfnDeal */
  if (err = grfn_copy_to_GrfnDeal(gd, loc_tab_rw, loc_tab_cl, tab, loc_evt_dts,
                                  0, NULL, aux_width, aux_len, aux)) {
    free(gd);
    free(und_info);
    free(loc_evt_dts);
    grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
    return err;
  }

  /*	Initialise the payments and their dates		*/
  gd->pay_dates = NULL;
  gd->event_pay_dates = NULL;
  gd->pay_amounts = NULL;
  gd->numpay = 0;

  /*	Set today	*/
  gd->today = ref;

  /*	Store the domestic underlying name and curve */
  strcpy((char *)gd->domestic_und, dom_nme);
  dom_crv_nme = get_discname_from_underlying(und);
  strcpy((char *)gd->disc_curve, dom_crv_nme);

  /*	Find index of first event date on or after today	*/
  gd->first_unkn_index = first_stop_idx;

  /*	Set up the end-of-day flags (for payments/exercise... and fixings) */
  gd->end_of_day_index = (prm->end_of_day_flg == SRT_YES ? 1 : 0);
  gd->end_of_day_fixing = (prm->end_of_day_fixings == SRT_YES ? 1 : 0);
  gd->end_of_day_payment = (prm->end_of_day_payment == SRT_YES ? 1 : 0);

  /*	Sets today's and previous dates in the GrfnDeal */
  if (gd->first_unkn_index == 0) {
    /*	No dates in the past	*/
    gd->nowdt = gd->today;
    gd->prvdt = gd->nowdt;
  } else {
    /*	Dates in the past */
    gd->nowdt = gd->event_dates[0];
    gd->prvdt = gd->nowdt;
  }

  /*	Stage 3 - Parse tableau        , pre-calculate schedules and get
                  information on required discount factors
          ========================================================
   */

  /*	Parse the deal and get info concerning underlyings in the table	*/
  if (err = get_deal_und_info(gd, und_info)) {
    free(und_info);
    free(loc_evt_dts);
    grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
    grfn_free_inGrfnDeal(gd);
    free(gd);
    return err;
  }

  /* Transforms the tableau strings into a Grfn command representation */

  event = (GrfnEvent **)calloc(loc_tab_rw, sizeof(GrfnEvent *));
  for (i = 0; i < loc_tab_rw; i++) {
    event[i] = NULL;
  }

  /*	Start at the top left hand corner of tableau
          (minus one row for potential AM event )*/
  gd->I = -1;
  gd->J = 0;

  /*	Parse historical events	*/
  gd->is_history = SRT_YES;

  /* Set the Historical Fixing function in the GrfnDeal */
  swp_f_GetFixingFunc(&(gd->fixing_fct));

  /*	While the next event date is before the first future (unknown) date */
  for (dte_idx = 0; dte_idx < first_stop_idx; dte_idx++) {
    /*	Creates the historical event        ,
            taking fixings        , cvg into account
            and evaluating what can be done	*/
    if (err = grfn_create_historical_event(gd, &(event[dte_idx]),
                                           loc_evt_dts[dte_idx],
                                           loc_evt_dts[dte_idx + 1])) {
      free(und_info);
      free(loc_evt_dts);
      grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
      grfn_free_inGrfnDeal(gd);
      free(gd);
      for (i = 0; i < loc_tab_rw; i++)
        grfn_free_GrfnEvent(event[i]);
      free(event);
      return err;
    }
  }

  /*	Parse non-historical events	*/

  gd->is_history = SRT_NO;
  /*	Loop through events continues */
  for (dte_idx = first_stop_idx; dte_idx < loc_tab_rw; dte_idx++) {
    /*	Create Grfn event */
    if (err = grfn_create_future_event(
            gd, &(event[dte_idx]), loc_evt_dts[dte_idx],
            dte_idx < loc_tab_rw - 1 ? loc_evt_dts[dte_idx + 1] : 0)) {
      free(und_info);
      free(loc_evt_dts);
      grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
      grfn_free_inGrfnDeal(gd);
      free(gd);
      for (i = 0; i < loc_tab_rw; i++)
        grfn_free_GrfnEvent(event[i]);
      free(event);
      return err;
    }
  }

  /*	Cut useless events off the end by finding last real event	*/
  dte_idx = loc_tab_rw - 1;

  while (dte_idx > 0 && event[dte_idx] == NULL) {
    dte_idx--;
    loc_tab_rw--;
  }

  if (prm->lsm == SRT_FORBACK) {
    prm->imp_type = IMPBACKWARDLATTICE;
  } else {
    /*	Set the Grfn direction	*/
    if (GrfnIsStatus(gd->sumstatus, GRFNCSFUTUREREF) &&
        GrfnIsStatus(gd->sumstatus, GRFNCSPASTREF)) {
      err = serror(GRERR_METHOD_CONFLICT);
      free(und_info);
      free(loc_evt_dts);
      grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
      grfn_free_inGrfnDeal(gd);
      free(gd);
      for (i = 0; i < loc_tab_rw; i++)
        grfn_free_GrfnEvent(event[i]);
      free(event);
      return err;
    }

    if (GrfnIsStatus(gd->sumstatus, GRFNCSPASTREF)) {
      gd->method = GRFNMONTECARLO;
    } else if (GrfnIsStatus(gd->sumstatus, GRFNCSFUTUREREF)) {
      gd->method = GRFNBACKWARDLATTICE;
    } else {
      gd->method = GRFNNOMETHOD;
    }

    if (gd->method == GRFNMONTECARLO ||
        (prm->force_mc == SRT_YES && gd->method != GRFNBACKWARDLATTICE)) {
      prm->imp_type = IMPMONTECARLO;
    } else {
      prm->imp_type = IMPBACKWARDLATTICE;
    }
  }

  *back_or_for = (prm->imp_type == IMPMONTECARLO ? 1 : -1);

  /*	Stage 4 - Value historical events
          =================================	*/

  /*	Set the Historical Fixing function in the GrfnDeal */
  swp_f_GetFixingFunc(&(gd->fixing_fct));

  /*	Sets the "is history" flag to true for evaluation	*/
  gd->is_history = SRT_YES;

  /*	Loop on historical events and evaluate	*/
  for (dte_idx = 0; dte_idx < first_stop_idx; dte_idx++) {
    /*	Evaluate the current event	*/
    if (err = grfn_eval_event(event[dte_idx],
                              NULL,        /*	No state variables	*/
                              gd, NULL,    /*	No PV so far */
                              calc_det_df, /*	Just get DFS */
                              und_info, &cash_flow)) {
      free(und_info);
      free(loc_evt_dts);
      grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
      grfn_free_inGrfnDeal(gd);
      free(gd);
      for (i = 0; i < loc_tab_rw; i++)
        grfn_free_GrfnEvent(event[i]);
      free(event);
      return err;
    }
  }

  /*	Makes sure all needed past variable information (va        , vb...)
          is attached to today's step */

  /*	If there is a front event (am event before first event date)	*/
  if (first_stop_idx > 0 && loc_evt_dts[first_stop_idx] > first_stop &&
      GrfnIsStatus(gd->rowstatus[first_stop_idx - 1], GRFNCSAMERICAN)) {
    front_evt = 1;
  } else {
    front_evt = 0;
  }

  if (err = grfn_attach_past_info_to_future(
          gd, &(event[first_stop_idx - front_evt]))) {
    free(und_info);
    free(loc_evt_dts);
    grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
    grfn_free_inGrfnDeal(gd);
    free(gd);
    for (i = 0; i < loc_tab_rw; i++)
      grfn_free_GrfnEvent(event[i]);
    free(event);
    return err;
  }

  /*	Sets the "is history" flag to false for evaluation	*/
  gd->is_history = SRT_NO;

  /*	Stage 5 - Fill structure to be returned
          =======================================	*/

  comm->gd = gd;
  comm->und_info = und_info;

  comm->num_und = und_info->no_of_underlyings;
  comm->und_data = (Und_Data *)(und_info->und_data);
  for (i = 0; i < comm->num_und; i++) {
    comm->und_ptr[i] = lookup_und(comm->und_data[i].und_name);
  }

  comm->num_cols = prm->node_dim = gd->sswidth;

  comm->num_evt = loc_tab_rw - first_stop_idx + front_evt;
  comm->dts = (long *)calloc(comm->num_evt, sizeof(long));
  comm->tms = (double *)calloc(comm->num_evt, sizeof(double));
  comm->evt = (FIRSTMktAtT *)calloc(comm->num_evt, sizeof(FIRSTMktAtT));
  comm->am = (int *)calloc(comm->num_evt, sizeof(int));

  if (front_evt) {
    comm->evt[0].evt = event[first_stop_idx - 1];
    comm->am[0] = 1;
    comm->dts[0] = first_stop;
    comm->tms[0] = (double)(first_stop - ref) / DAYS_IN_YEAR;
  }

  for (dte_idx = 0, dte_idx2 = first_stop_idx; dte_idx2 < loc_tab_rw;
       dte_idx++, dte_idx2++) {
    comm->evt[dte_idx + front_evt].evt = event[dte_idx2];
    comm->am[dte_idx + front_evt] =
        (GrfnIsStatus(gd->rowstatus[dte_idx2], GRFNCSAMERICAN) ? 1 : 0);

    comm->dts[dte_idx + front_evt] = loc_evt_dts[dte_idx2];
    comm->tms[dte_idx + front_evt] =
        (double)(loc_evt_dts[dte_idx2] - ref) / DAYS_IN_YEAR;
  }

  /*	Stage 6 - Free & return
          =======================	*/

  free(loc_evt_dts);
  grfn_free_GrfnCellmatrix(tab, loc_tab_rw, loc_tab_cl);
  for (i = 0; i < first_stop_idx - front_evt; i++)
    grfn_free_GrfnEvent(event[i]);
  free(event);

  return NULL;
}

/*
        Functions used to get the global deal information
        =================================================
*/

/*	1.	Get underlyings to be used to value the deal	*/

Err FIRSTGetUndFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Number of underlyings needed to value the deal	*/
    int *num_und,
    /*	Pointer is allocated inside        ,
            and must be freed by call to FIRSTFreeUndFromDeal	*/
    SrtUndPtr **und_ptr) {
  int i;

  *num_und = comm->und_info->no_of_underlyings;
  *und_ptr = (SrtUndPtr *)calloc(*num_und, sizeof(SrtUndPtr));

  for (i = 0; i < *num_und; i++) {
    (*und_ptr)[i] = lookup_und(comm->und_info->und_data[i].und_name);
  }

  return NULL;
}

/*	2.	Free the previous function result	*/

Err FIRSTFreeUndFromDeal(int num_und, SrtUndPtr **und_ptr) {
  if (*und_ptr)
    free(*und_ptr);
  return NULL;
}

/*	3.	Get number of columns	*/

Err FIRSTGetNumColFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Information	*/
    int *num_col) {
  *num_col = comm->num_cols;
  return NULL;
}

/*	4.	Get pv of past (to be added to last column)	*/

Err FIRSTGetPvOfPastFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	PV to be added to last column before simulation (MC) or after (Tree)
     */
    double *pv_of_past) {
  *pv_of_past = comm->gd->pv_of_past;
  return NULL;
}

/*	5.	Get the maximum number of dfs required	*/

Err FIRSTGetMaxNumDfFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Max num df	*/
    int *max_num_df) {
  int i, j;
  int m = 0;

  for (i = 0; i < comm->num_evt; i++) {
    if (comm->evt[i].evt) {
      for (j = 0; j < comm->und_info->no_of_underlyings; j++) {
        m = IMAX(m, comm->evt[i].evt->dflen[j]);
      }
    }
  }

  *max_num_df = m;
  return NULL;
}

/*	6.	Get event dates	*/

Err FIRSTGetEvtDatesFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Number of event dates	*/
    int *num_evt,

    /*	Pointers are allocated inside        ,
            and must be freed by call to FIRSTFreeEvtDatesFromDeal	*/

    /*	Event dates and times	*/
    long **evt_dts, double **evt_tms) {
  *num_evt = comm->num_evt;

  *evt_dts = (long *)calloc(*num_evt, sizeof(long));
  *evt_tms = (double *)calloc(*num_evt, sizeof(double));

  if (!*evt_dts || !*evt_tms) {
    return "Memory allocation error in MdlComm -FIRSTGetEvtDatesFromDeal";
  }

  memcpy(*evt_dts, comm->dts, *num_evt * sizeof(long));
  memcpy(*evt_tms, comm->tms, *num_evt * sizeof(double));

  return NULL;
}

/*	7.	Free the previous function result	*/

Err FIRSTFreeEvtDatesFromDeal(int num_evt, long **evt_dts, double **evt_tms) {
  if (*evt_dts)
    free(*evt_dts);
  if (*evt_tms)
    free(*evt_tms);

  return NULL;
}

/*	8.	Get events and df maturities required to calculate them	*/

Err FIRSTGetEventInfoFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	As output from FIRSTGetEvtDatesFromDeal	*/
    int num_evt,
    /*	As output from FIRSTGetUndFromDeal	*/
    int num_und,
    /*	As output from FIRSTGetMaxNumDfFromDeal	*/
    int max_num_df,

    /*	All pointers are allocated inside        ,
            and must be freed by call to FIRSTFreeEventInfoFromDeal	*/

    /*	Events themselves	*/
    FIRSTMktAtT **evt,

    /*	Wether American	*/
    int **am,

    /*	Information relative to df required for event evaluation	*/

    /*	num_df_mat[i][j] = number of df required for event i from underlying j
     */
    int ***num_df_mat,
    /*	df_mat_dts[i][j][k] and df_mat_tms[i][j][k]
            = maturity of the df number k required from underlying j at event i
     */
    long ****df_mat_dts, double ****df_mat_tms)

{
  int i, j, k;
  int m = 0;

  *num_df_mat = imatrix(0, num_evt - 1, 0, num_und - 1);
  *df_mat_dts = l3tensor(0, num_evt - 1, 0, num_und - 1, 0, max_num_df - 1);
  *df_mat_tms = f3tensor(0, num_evt - 1, 0, num_und - 1, 0, max_num_df - 1);

  if (!*num_df_mat || !*df_mat_dts || !*df_mat_tms) {
    return "Memory allocation error in MdlComm -FIRSTGetEventInfoFromDeal";
  }

  for (i = 0; i < comm->num_evt; i++) {
    if (comm->evt[i].evt) {
      for (j = 0; j < comm->und_info->no_of_underlyings; j++) {
        (*num_df_mat)[i][j] = comm->evt[i].evt->dflen[j];

        for (k = 0; k < (*num_df_mat)[i][j]; k++) {
          (*df_mat_dts)[i][j][k] = comm->evt[i].evt->dfd[j][k];
          (*df_mat_tms)[i][j][k] = comm->evt[i].evt->dft[j][k];
        }
      }
    } else {
      for (j = 0; j < comm->und_info->no_of_underlyings; j++) {
        (*num_df_mat)[i][j] = 0;
      }
    }
  }

  *evt = comm->evt;
  *am = comm->am;
  return NULL;
}

/*	9.	Free the previous function result	*/

Err FIRSTFreeEventInfoFromDeal(int num_evt, int num_und, int max_num_df,
                               FIRSTMktAtT **evt, int **am, int ***num_df_mat,
                               long ****df_mat_dts, double ****df_mat_tms) {
  if (*num_df_mat)
    free_imatrix(*num_df_mat, 0, num_evt - 1, 0, num_und - 1);
  if (*df_mat_dts)
    free_l3tensor(*df_mat_dts, 0, num_evt - 1, 0, num_und - 1, 0,
                  max_num_df - 1);
  if (*df_mat_tms)
    free_f3tensor(*df_mat_tms, 0, num_evt - 1, 0, num_und - 1, 0,
                  max_num_df - 1);

  return NULL;
}

/*
        Functions used to get the local deal information
        ================================================
*/

/*	1.	Get event and df maturities required at one given time	*/

Err FIRSTGetLocalEventInfoFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Current time	*/
    double time,
    /*	Return NULL if no event at this date        , or pointer on the event
     */
    FIRSTMktAtT **evt,

    /*	POINTERS num_df_mat        , df_mat_dts and df_mat_tms
            ARE NOT ALLOCATED INSIDE
            they must be allocated prior to function call	*/

    /*	num_df_mat[i] = number of df required from underlying i	*/
    int *num_df_mat,
    /*	df_mat_dts[i][j] and df_mat_tms[i][j]
            = maturity of the df number j required from underlying i	*/
    long **df_mat_dts, double **df_mat_tms) {
  int i, j, k, found_time;

  if (time > comm->tms[comm->num_evt - 1] + ONE_MINUTE) {
    *evt = NULL;
    return NULL;
  }

  if (time < comm->tms[0] - ONE_MINUTE) {
    *evt = NULL;
    return NULL;
  }

  i = 0;
  while (comm->tms[i] < time - ONE_MINUTE) {
    i++;
  }

  if (comm->tms[i] > time + ONE_MINUTE) {
    if (comm->am[i - 1]) {
      found_time = i - 1;
    } else {
      *evt = NULL;
      return NULL;
    }
  } else {
    found_time = i;
  }

  *evt = &(comm->evt[found_time]);

  if ((*evt) && (*evt)->evt) {
    for (j = 0; j < comm->und_info->no_of_underlyings; j++) {
      num_df_mat[j] = (*evt)->evt->dflen[j];

      for (k = 0; k < num_df_mat[j]; k++) {
        df_mat_dts[j][k] = (*evt)->evt->dfd[j][k];
        df_mat_tms[j][k] = (*evt)->evt->dft[j][k];
      }
    }
  } else {
    for (j = 0; j < comm->und_info->no_of_underlyings; j++) {
      num_df_mat[j] = 0;
    }
  }

  return NULL;
}

/*	2.	Set DF values prior to valuation	*/

Err FIRSTSetDFValue(
    /*	Pointer on the event to be valued        ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Index of the underlying	*/
    int und_idx,
    /*	Index of the df	*/
    int df_idx,
    /*	Value of df	*/
    double df_val) {
  evt->evt->df[und_idx][df_idx] = df_val;
  return NULL;
}

/*	3.	Set state variable values prior to valuation	*/

Err FIRSTSetSVValue(
    /*	Pointer on the event to be valued        ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Index of the underlying	*/
    int und_idx,
    /*	Type of the statevar (SPOT        , R        , PHI        , ...)
     */
    int sv_type,
    /*	Value of statevar	*/
    double sv_val) {
  evt->smp.und[und_idx].sv[sv_type] = sv_val;
  return NULL;
}

/*
        Functions used to eval event
        ============================

        DFs and stevars are supposed to have been set accordingly by caller
        before calling this function
*/

Err FIRSTEvalEvent(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Pointer on the event to be valued        ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Number of columns	*/
    int num_col,
    /*	Information needed by Forward/Backward models only
            Please ask Eric Fournie what these are about...
            If you don't want to use these        , please set them to default
     */
    int type_eval, /*	Default is 2	*/
    double *fwd,   /*	Default is NULL	*/
    double *cur,   /*	Default is NULL	*/
    /*	Vector of PVs so far of columns        , will be updated by function
            tree only	*/
    double *col_pvs,
    /*	Cash flow
            MC only	*/
    double *cash_flow)

{
  if (evt && evt->evt) {
    return grfn_eval_event_ammc(type_eval, evt->evt, &(evt->smp), comm->gd,
                                col_pvs, fwd, cur, NULL, comm->und_info,
                                cash_flow);
  } else {
    return NULL;
  }
}

// IMPLEMENTATION OF PVR
Err FIRSTEvalEventCredit(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Pointer on the event to be valued        ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Number of columns	*/
    int num_col,
    /*	Information needed by Forward/Backward models only
            Please ask Eric Fournie what these are about...
            If you don't want to use these        , please set them to default
     */
    int type_eval, /*	Default is 2	*/
    double *fwd,   /*	Default is NULL	*/
    double *cur,   /*	Default is NULL	*/
    /*	Vector of PVs so far of columns        , will be updated by function
            tree only	*/
    double *col_pvNR, double *col_pvR, double *col_pvSource, int *sourceStatus,
    /*	Cash flow
            MC only	*/
    double *cash_flow)

{
  if (evt && evt->evt) {
    return grfn_eval_event_credit(type_eval, evt->evt, &(evt->smp), comm->gd,
                                  col_pvNR, col_pvR, col_pvSource, sourceStatus,
                                  fwd, cur, NULL, comm->und_info, cash_flow);
  } else {
    return NULL;
  }
}

#if 0
Err FIRSTEvalEvent2(
/*	As output from FIRSTInitMktStruct	*/
GRFNCOMMSTRUCT			comm        ,
/*	Pointer on the event to be valued        , 
	must be a valid non-empty event	*/
FIRSTMktAtT*			evt        ,
/*	Number of columns	*/
int						num_col        ,
/*	Information needed by Forward/Backward models only
	Please ask Eric Fournie what these are about...
	If you don't want to use these        , please set them to default	*/
int						type_eval        ,	/*	Default is 2	*/
double*					fwd        ,		/*	Default is NULL	*/
double*					cur        ,		/*	Default is NULL	*/
/*	Vector of PVs so far of columns        , will be updated by function
	tree only	*/
double*					col_pvs        ,
/*	Cash flow
	MC only	*/
double*					cash_flow)

{
	if (evt && evt->evt)
	{
		return grfn_eval_event_ammc2 (	type_eval        ,
										evt->evt        ,
										&(evt->smp)        ,
										comm->gd        ,
										col_pvs        ,
										fwd        ,
										cur        ,
										NULL        ,
										comm->und_info        ,
										cash_flow);
	}
	else
	{
		return NULL;
	}
}

#endif

/*
        Functions used to free everything after valuation (or on error)
        ===============================================================
*/

Err FIRSTFreeMktStruct(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm) {
  int i;

  grfn_free_inGrfnDeal(comm->gd);
  free(comm->gd);

  free(comm->und_info);

  free(comm->dts);
  free(comm->tms);

  free(comm->am);

  for (i = 0; i < comm->num_evt; i++) {
    grfn_free_GrfnEvent(comm->evt[i].evt);
  }
  free(comm->evt);

  return NULL;
}

/*
        A simple example
        ================
*/

/*	1.	Function to free underlying	*/

Err SrtFreeDetermUnd(void *undPtr) {
  /*	Nothing to free inside in this case	*/
  free(undPtr);
  return NULL;
}

/*	2.	Initialise deterministic model	*/

Err SrtInitDetermUnd(char *undName, char *ycName) {
  Err err;
  SrtUndListPtr und_list;
  SrtUndPtr und;
  SrtIrDesc *ir_und;
  long today;
  SrtCurvePtr crv;
  char *ccy;

  /*	Get the Yield Curve */
  if ((crv = lookup_curve(ycName)) == NULL) {
    return serror("Fatal: (SrtInitDetermUnd) Cannot find yield curve");
  }
  if (!ISCURVETYPE(crv, YIELD_CURVE)) {
    return serror("Fatal: (SrtInitDetermUnd) Must define a yield curve object");
  }

  /*	Extract today and currency from this curve */
  today = get_today_from_curve(crv);
  ccy = get_curve_ccy(crv);

  /*	Get the underlying list and check if not empty */
  und_list = get_underlying_list();
  if (und_list == NULL) {
    return serror("No Underlying list defined: call SrtInit before");
  }

  /*	Makes the Underlying Name UpperCase */
  strupper(undName);

  /*	Create space for SrtUndDesc */
  und = (SrtUndPtr)malloc(sizeof(SrtUndDesc));

  /*	Fill general info	*/
  strcpy(und->underl_name, undName);
  strupper(und->underl_name);
  strip_white_space(und->underl_name);
  strcpy(und->underl_lbl, "IR_UND");
  und->underl_ccy = ccy;
  und->underl_type = INTEREST_RATE_UND;

  /*	Fill IR info	*/
  ir_und = (SrtIrDesc *)malloc(sizeof(SrtIrDesc));
  strcpy(ir_und->mdl_lbl, "DETERMINISTIC");
  ir_und->mdl_type = DETERMINISTIC;
  ir_und->mdl_dim = ONE_FAC;
  ir_und->ts = NULL;
  strcpy(ir_und->yc_name, ycName);
  ir_und->spec = NULL;

  und->spec_desc = ir_und;

  /*	Insert (overwrite & update ticker) the object
          in the underlying list (with the same name) */
  err = srt_f_lstins(und_list, undName, 0.0, OBJ_PTR_UND, und, SrtFreeDetermUnd,
                     &(und->underl_ticker));

  if (err) {
    return serror("Error in initialising %s object", und->underl_lbl);
  }

  /*	Return a success message */
  return NULL;
}

/*	3.	Value GRFN deal in deterministic model	*/

Err FirstValueDealInDeterministicIrModel(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	Premium	*/
    double *price) {
  Err err;

  long ref;

  int back_or_for;
  FIRSTAllMkts comm;

  int num_und;
  SrtUndPtr *und_ptr;
  char **yc_name;

  int num_col;

  double pv_of_past;

  int max_num_df;

  int num_evt;
  long *evt_dts;
  double *evt_tms;

  FIRSTMktAtT *evt;
  int *am;
  int **num_df_mat;
  long ***df_mat_dts;
  double ***df_mat_tms;

  double *pvs;

  int i, j, k;

  double cf;

  smessage("This is an example        , "
           "all underlyings are supposed to be deterministic "
           "and interest rates only");

  /*	Init and check deal        , value history	*/
  if (err = FIRSTInitMktStruct(init_num_evt_dts, init_evt_dts, tab_rw, tab_cl,
                               tab_str, mask, aux_width, aux_len, aux, dom_nme,
                               prm, &back_or_for, &comm)) {
    return err;
  }

  /*	Get relevant underlyings	*/
  if (err = FIRSTGetUndFromDeal(&comm, &num_und, &und_ptr)) {
    FIRSTFreeMktStruct(&comm);
    return err;
  }

  /*	Get today's date	*/
  ref = get_today_from_underlying(und_ptr[0]);

  /*	Get curve names	*/
  yc_name = (char **)calloc(num_und, sizeof(char *));
  for (j = 0; j < num_und; j++) {
    get_underlying_discname(und_ptr[j], &(yc_name[j]));
  }

  /*	Get number of columns	*/
  if (err = FIRSTGetNumColFromDeal(&comm, &num_col)) {
    FIRSTFreeMktStruct(&comm);
    FIRSTFreeUndFromDeal(num_und, &und_ptr);
    free(yc_name);
    return err;
  }

  /*	Get pv of cash-flows fixed in the past        , but not paid yet
   */
  if (err = FIRSTGetPvOfPastFromDeal(&comm, &pv_of_past)) {
    FIRSTFreeMktStruct(&comm);
    FIRSTFreeUndFromDeal(num_und, &und_ptr);
    free(yc_name);
    return err;
  }

  /*	Get maximum number of discount factors required (for allocation)
   */
  if (err = FIRSTGetMaxNumDfFromDeal(&comm, &max_num_df)) {
    FIRSTFreeMktStruct(&comm);
    FIRSTFreeUndFromDeal(num_und, &und_ptr);
    free(yc_name);
    return err;
  }

  /*	Get actual event dates (as opposed to initial ones        , today is
     added and empty events are removed)	*/
  if (err = FIRSTGetEvtDatesFromDeal(&comm, &num_evt, &evt_dts, &evt_tms)) {
    FIRSTFreeMktStruct(&comm);
    FIRSTFreeUndFromDeal(num_und, &und_ptr);
    free(yc_name);
    return err;
  }

  /*	Get events and maturities of discount factors required	*/
  if (err = FIRSTGetEventInfoFromDeal(&comm, num_evt, num_und, max_num_df, &evt,
                                      &am, &num_df_mat, &df_mat_dts,
                                      &df_mat_tms)) {
    FIRSTFreeMktStruct(&comm);
    FIRSTFreeUndFromDeal(num_und, &und_ptr);
    FIRSTFreeEvtDatesFromDeal(num_evt, &evt_dts, &evt_tms);
    free(yc_name);
    return err;
  }

  pvs = (double *)calloc(num_col, sizeof(double));

  /*	If forward: add pv of past	*/
  if (back_or_for == 1) {
    (*price) += pv_of_past;
    pvs[num_col - 1] = pv_of_past;
  }

  /*	Main loop	*/
  for (i = (back_or_for == 1 ? 0 : num_evt - 1);
       (back_or_for == 1 ? (i < num_evt) : (i > -1)); i += back_or_for) {
    /*	Set values of dfs to fwd dfs	*/
    for (j = 0; j < num_und; j++) {
      for (k = 0; k < num_df_mat[i][j]; k++) {
        if (err = FIRSTSetDFValue(
                evt + i, j, k,
                swp_f_df(evt_dts[i], df_mat_dts[i][j][k], yc_name[j]))) {
          FIRSTFreeMktStruct(&comm);
          FIRSTFreeUndFromDeal(num_und, &und_ptr);
          FIRSTFreeEvtDatesFromDeal(num_evt, &evt_dts, &evt_tms);
          FIRSTFreeEventInfoFromDeal(num_evt, num_und, max_num_df, &evt, &am,
                                     &num_df_mat, &df_mat_dts, &df_mat_tms);
          free(yc_name);
          free(pvs);
          return err;
        }
      }
    }

    /*	Value event	*/
    if (err =
            FIRSTEvalEvent(&comm, evt + i, num_col, 2, NULL, NULL, pvs, &cf)) {
      FIRSTFreeMktStruct(&comm);
      FIRSTFreeUndFromDeal(num_und, &und_ptr);
      FIRSTFreeEvtDatesFromDeal(num_evt, &evt_dts, &evt_tms);
      FIRSTFreeEventInfoFromDeal(num_evt, num_und, max_num_df, &evt, &am,
                                 &num_df_mat, &df_mat_dts, &df_mat_tms);
      free(yc_name);
      free(pvs);
      return err;
    }

    /*	If fwd: increment price	*/
    if (back_or_for == 1) {
      (*price) += cf * swp_f_df(ref, evt_dts[i], yc_name[0]);
    } else {
      if (i > 0) {
        /*	Discount	*/
        for (j = 0; j < num_col; j++) {
          pvs[j] *= swp_f_df(evt_dts[i - 1], evt_dts[i], yc_name[0]);
        }
      }
    }
  }

  /*	If backward: price is pv	*/
  if (back_or_for == -1) {
    *price =
        pv_of_past + pvs[num_col - 1] * swp_f_df(ref, evt_dts[0], yc_name[0]);
  }

  FIRSTFreeMktStruct(&comm);
  FIRSTFreeUndFromDeal(num_und, &und_ptr);
  FIRSTFreeEvtDatesFromDeal(num_evt, &evt_dts, &evt_tms);
  FIRSTFreeEventInfoFromDeal(num_evt, num_und, max_num_df, &evt, &am,
                             &num_df_mat, &df_mat_dts, &df_mat_tms);
  free(yc_name);
  free(pvs);

  return NULL;
}