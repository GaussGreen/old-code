
/* =============================================================================

   FILENAME:    swp_f_datelist.cxx

   PURPOSE:     Generation of a swap schedule

   =============================================================================
 */

#include <swp_h_all.h"

SrtDateList new_DateList(int len) {
  SrtDateList x;

  memset(&x, 0, sizeof(SrtDateList));
  x.date = (long *)(srt_calloc(len, sizeof(long)));
  x.len = len;

  return x;
}

/* -----------------------------------------------------------------------------
  PURPOSE:	From THEORETICAL (i.e. non holidays adjusted) start and
                        end dates (as payment -or start of accrual period-
  dates)        , generates a list of EFFECTIVE payment (or start of accrual
  period) dates        , that are business days adjusted according to
  SrtBusDayConv for a swap with periodicity compd

  AMENDMENT: E.Auld 5 Jan 94 --- added bus_day_start.
                        previously        , swap->start was used        , which
  caused an error when swap->start was not a business day.
  -----------------------------------------------------------------------------
*/

SrtDateList SwapDP_to_DateList(SrtSwapDP *swap, SrtBusDayConv conv) {
  int i, k, length;
  Date cur_date, count_date, bus_day_start;
  long *d;
  DateList list;

  if (swap->direction == FWD) {
    length = 3 + swap->nfp;
  } else {
    length = (year(swap->end) - year(swap->start) + 1) * swap->compd + 1;
  }

  list.date = (long *)srt_calloc(length, sizeof(Date));
  d = (long *)srt_calloc(length, sizeof(Date));

  list.len = 0;
  list.dir = swap->direction;

  switch (swap->direction) {
  case FWD:
    if (swap->start == swap->first_full_fixing) {
      list.type = NOTBROKEN;
      cur_date = swap->start;
      count_date = swap->start;
      list.date[0] = cur_date;
      cur_date = add_unit(count_date, (int)(12 / swap->compd), SRT_MONTH,
                          NO_BUSDAY_CONVENTION);
      k = 1;
    } else {
      /* Sets the first date in the list (broken coupon) */
      /* before going for full swap */
      list.type = BROKEN;
      list.date[0] = swap->start;
      cur_date = swap->first_full_fixing;
      count_date = swap->first_full_fixing;
      k = 0;
    }

    for (i = 1; i <= swap->nfp; i++) {
      /* The effective period start (for accrual coupon..) is the modified
       * theoretical date */
      list.date[i] = cur_date;
      /* The next theoretical payment date (== period start for accrual) is non
       * modified */
      cur_date = add_unit(count_date, (int)(12 / swap->compd) * (i + k),
                          SRT_MONTH, NO_BUSDAY_CONVENTION);
    }

    /* Sets swap->end to the theoretical end date of the swap */
    swap->end = list.date[i - 1];
    list.len = i;
    break;

  case BKWD:
    cur_date = swap->end;
    count_date = swap->end;
    list.type = NOTBROKEN;

    for (i = 0; cur_date >= swap->start; i++) {
      /* The theroretical period start (for accrual coupon..) */
      d[i] = cur_date;
      /* The previous theoretical payment date */
      cur_date = add_unit(count_date, -(int)(12 / swap->compd) * (i + 1),
                          SRT_MONTH, NO_BUSDAY_CONVENTION);
    }

    length = i;

    /* Real business days modified start date of the swap */
    bus_day_start = bus_date_method(swap->start, conv);

    /* Here        , cur_date is before swap->start: check the next one to see
            if broken: has to be done with bus_day_start if for instance ,
            swap->start = SATURDAY        , and cur_date + period = MONDAY:
            the swap is then not considered as broken (effect. start on MONDAY)
     */
    if (add_unit(cur_date, (int)(12 / swap->compd), SRT_MONTH,
                 NO_BUSDAY_CONVENTION) > bus_day_start) {
      d[i] = bus_day_start;
      length += 1;
      list.type = BROKEN;
    }

    /* Checks that applying business day rule will not set d[length-2] =
     * d[length-1] */
    if (bus_date_method(d[length - 2], conv) ==
        bus_date_method(d[length - 1], conv)) {
      length--;
    }

    /* Inverts the list */
    for (i = 0; i < length; i++) {
      list.date[i] = d[length - 1 - i];
    }
    list.len = length;

    break;

  } /* END switch */

  /* Gets the theoretical previous payment date from start or first full coupon
   */
  if (list.type == NOTBROKEN) {
    list.prev =
        add_unit(list.date[list.len - 1], -(12 / swap->compd) * list.len,
                 SRT_MONTH, NO_BUSDAY_CONVENTION);
  } else {
    list.prev =
        add_unit(list.date[list.len - 1], -(12 / swap->compd) * (list.len - 1),
                 SRT_MONTH, NO_BUSDAY_CONVENTION);
  }
  /* Sets the associated bus day modified date */
  list.prev = bus_date_method(list.prev, conv);

  /* Now stores the real swap dates        , as the modified business date */
  for (i = 0; i < list.len; i++) {
    list.date[i] = bus_date_method(list.date[i], conv);
  }

  /* Free memory */
  srt_free(d);

  return (list);
}

/***** return the part of dl that is after the date after
        this is a NEW date list (not just a pointer to a half
        of dl) ***/
SrtDateList DateListafter(SrtDateList dl, Date after) {
  DateList ndl;
  int i = 0;

  while (i < dl.len && dl.date[i] <= after)
    i++;
  ndl.len = dl.len - i;
  ndl.date = (long *)srt_calloc(ndl.len, sizeof(Date));
  memcpy(ndl.date, dl.date + i, ndl.len * sizeof(Date));
  return ndl;
}

/* Frees the dates in the date list */
Err swp_f_free_in_DateList(SrtDateList datelist) {
  if (datelist.date)
    srt_free(datelist.date);

  return NULL;
}
