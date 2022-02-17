#include "grf_h_all.h"

typedef struct {
  Date date;
  int comll_index;
  int df_index;
} Data;

/* Compare two Data        , first by date        , then by COMLL index */

static int compare(const void *ptr_a, const void *ptr_b) {
  Data *a = (Data *)(ptr_a);
  Data *b = (Data *)(ptr_b);

  if (a->date < b->date)
    return -1;
  if (a->date > b->date)
    return 1;
  if (a->comll_index < b->comll_index)
    return -1;
  if (a->comll_index > b->comll_index)
    return 1;

  return 0;
}

/******************************************************************************\

    FUNCTION        :   grfn_set_df_dates_in_event

    DESCRIPTION     :   Each part of the comll in the new_event->icom field
                        is assumed to have an array of dates in its dfind field
                        This function sorts all these dates into one large
                        array of dates        , (in years) which goes into the
dft field of event.

                                                Df are gathered only underlying
by underlying (there is a loop in grf_f_get_event that calls this function und
by und). It is therefore important to make sure the COMLL we are working on is
for the underlying we consider.

                        The dates in the dfind & dfindfloat fields are replaced
                                                with indices into the big dft
array.

    NOTE            :   [American events]

                        The start day of a swap (for example) may have been
                        moved forward by "n" days + amdt; "n" is assumed
                        to have been taken care of before this function was
                        called. "amdt" does not mean that all the future
                        payments of the swap occur at time amdt in the day.

                        It means only that the times until then are reduced
                        by amdt*years_in_day.  The payment days do NOT change.

\******************************************************************************/

Err grfn_set_df_dates_in_event(GrfnEvent *new_event, GrfnDeal *gd,
                               int und_index)

{
  /*..Local variables..*/
  int i = 0, j = 0, k = 0, k_float = 0;
  int ndf, sz, maxl;
  int tsz;
  int comll_index;
  double *dft = NULL, *df = NULL, amdt;
  Date *dfd = NULL;
  Data *data;
  Date today, idate;
  COMLL_PTR comll, temp;

  /*..Initialise local variables..*/
  ndf = 0;
  sz = 0;
  tsz = 0;
  maxl = 0;
  amdt = 0.0;
  comll = new_event->icom;
  temp = comll;
  today = gd->today;
  idate = gd->event_dates[gd->I];

  /*..Main Body..*/
  if (gd->american == SRT_YES) {
    amdt = gd->amdt;
    idate += gd->amdays;
  }

  /* Determine the number of df's in the full COMLL that concern this underlying
   */
  while (temp != NULL) {
    sz++;

    /* Increment df size only if und_index of und in COMLL == input und index */
    tsz += (temp->dfindlen) * ((temp->ivec[0] == und_index) ? 1 : 0);
    /* Same for potential floating leg */
    tsz += (temp->dfindfloatlen) * ((temp->ivec[0] == und_index) ? 1 : 0);

    temp = temp->next;
  }

  /* Check if there is gathering work to be done : df's might be equired for
   * this und */
  if (tsz > 0) {
    /* Memory allocation for the vector into which df's will be transfered */
    data = srt_calloc(tsz, sizeof(Data));
    if (!data)
      return serror("Allocation failure (1) in grfn_init_event");

    /* We want to put all the dates        , and the index of the comlls to
   which they belong        , in one array of Data (loop on COMLL for the
   transfer) */
    k = 0;
    temp = comll;
    comll_index = 0;
    while (temp != NULL) {
      /* Transfer df_dates only if COMLL und index == input und index */
      if (temp->ivec[0] == und_index) {
        /* Loop on all df's for the fixed leg (no spreads) */
        for (j = 0; j < temp->dfindlen; j++) {
          /* Add into data a new point (#k) with this df date        , and the
             comll_index (there is a cheat: even index for fixed leg: 2*index)
           */
          data[k].date = temp->dfind[j];
          data[k].cxxomll_index = 2 * comll_index;
          k++;
        }
        /* Loop on the df's for the floating leg (for spreads...) */
        for (j = 0; j < temp->dfindfloatlen; j++) {
          /* Add into data a new point (#k) with this df date        , and the
             comll_index
                  (there is a cheat: odd index for floating leg: 2*index + 1) */
          data[k].date = temp->dfindfloat[j];
          data[k].cxxomll_index = 2 * comll_index + 1;
          k++;
        }
      }

      temp = temp->next;
      comll_index++;

    } /* END while(temp != NULL) loop on all COMLL's */

    /* Sort array of dates (by dates first and then by comll_index if needed) */
    qsort(data, tsz, sizeof(Data), compare);

    /* Set df index in data s.t. no dates are counted twice (ndf is number of
     * unique dates) */
    data[0].df_index = 0;
    ndf = 1;
    for (i = 1; i < tsz; i++) {
      if (data[i].date == data[i - 1].date)
        data[i].df_index = data[i - 1].df_index;
      else {
        data[i].df_index = data[i - 1].df_index + 1;
        ndf++;
      }
    }

    /* Change dates in the original COMLL (dfind & dfindfloat) to indices in big
     * df array */
    for (temp = comll, comll_index = 0; temp != NULL;
         temp = temp->next, comll_index++) {
      k = 0;
      k_float = 0;

      for (j = 0; j < tsz; j++) {
        /* Reconciliate fixed leg with even comll_index in data */
        if ((data[j].cxxomll_index == 2 * comll_index) && (temp->dfindlen > 0))
          temp->dfind[k++] = data[j].df_index;
        else
            /* Reconciliate floating leg with odd comll_index in data */
            if ((data[j].cxxomll_index == 2 * comll_index + 1) &&
                (temp->dfindfloatlen > 0))
          temp->dfindfloat[k_float++] = data[j].df_index;
      }

    } /* END for loop */

    /* Move dates from data to array of times */
    dfd = srt_calloc(ndf, sizeof(Date));
    dft = srt_calloc(ndf, sizeof(double));
    df = srt_calloc(ndf, sizeof(double));

    if (dfd == NULL || dft == NULL || df == NULL)
      return serror("Allocation failure (3) in grfn_init_event");

    k = 1;

    /*??? FIX --> What if data[0].date == idate */

    dft[0] = (double)(data[0].date - idate - amdt) * YEARS_IN_DAY;
    dfd[0] = (Date)data[0].date;

    for (i = 1; i < tsz; i++) {
      if (data[i].df_index > data[i - 1].df_index) {
        dft[k] = (double)(data[i].date - idate - amdt) * YEARS_IN_DAY;
        dfd[k++] = (Date)data[i].date;
      }
    }

    srt_free(data);

  } /* end if tsz > 0 */

  /* We have done the dirty work now assign the results to new->event */

  new_event->d = idate;
  new_event->t = (double)(idate - today + amdt) * YEARS_IN_DAY;

  new_event->dfd[und_index] = dfd;
  new_event->dft[und_index] = dft;
  new_event->df[und_index] = df;
  new_event->dflen[und_index] = ndf;

  /*  Removed a chunk of code to srt_f_grfn.cxx */

  return NULL;

} /* END grfn_set_und_df_dates_in_event(..) */
