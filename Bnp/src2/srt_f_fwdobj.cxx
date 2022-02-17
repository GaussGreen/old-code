/* =======================================================================================

   FILENAME :       srt_f_fwdobj.cxx

   PURPOSE:         structures and functions to store a forward curve for use in
   Grfn A ForwardCurve is stored as a SrtDvdObj which is a SrtListHdr

   =======================================================================================
 */

#include "math.h"
#include "srt_h_fwdobj.h"
#include "utallhdr.h"

/* ----------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------
   FUNCTIONS TO CREATE        , INSERT AN ELEMENT IN OR DELETE A SrtDvdObj ==
   SrtList
   ------------------------------------------------------------------------------------
 */

/* Creates a SrtDvdObj with only one first element containing today */
Err srt_f_dvdobj_create(Date today, char *dvd_name, SrtDvdObj **dvdobj) {
  Err err = NULL;
  char list_name[SRTBUFSZ];
  double dvd_today;
  long ticker;

  sprintf(list_name, "%SrtDvdObj", dvd_name);

  err = srt_f_lstcreate(dvdobj, list_name);
  if (err)
    return serror("Error in creating the %s SrtDvdObj", dvd_name);

  dvd_today = 0.0;
  err = srt_f_lstins(*dvdobj, "_TODAY", (double)today, OBJ_DOUBLE,
                     (void *)(&dvd_today), NULL, &ticker);
  if (err)
    return serror("Error in inserting today in the %s SrtDvdObj", dvd_name);

  return err;

} /* Err srt_f_dvdobj_create(...) */

/* ----------------------------------------------------------------------------------
 */

/* Delete a SrtDvdObject and everything it contains */
Err srt_f_dvdobj_delete(SrtDvdObj **dvdobj) {
  /* Delete the SrtDvdObj = SrtList */
  if (srt_f_lstfree(*dvdobj, SRT_YES) != NULL)
    return serror("Error in srt_f_dvdobj_delete");

  /* Free the memory allocated for the object */
  srt_free(*dvdobj);

  return NULL;

} /* END Err srt_f_dvdobj_delete(...) */

/* ----------------------------------------------------------------------------------
 */

/* Insert a point (date        , value) in the dividend object */
Err srt_f_dvdobj_insertpoint(SrtDvdObj *dvdobj, Ddate date, double dvd) {
  Err err = NULL;
  long ticker;

  /* Insert a forward point == object in the SrtList = SrtDvdObj */
  err = srt_f_lstins(dvdobj, "DVD", (double)date, OBJ_DOUBLE, (void *)(&dvd),
                     NULL, &ticker);
  if (err)
    return err;

  /* Return a success message */
  return NULL;

} /* END Err srt_f_fwdobj_insertpoint(...) */

/* ----------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------
   STATIC FUNCTIONS TO FIND THE PREVIOUS OR NEXT ELEMENT IN THE SrtDvdObj with a
   ----------------------------------------------------------------------------------
 */

/* Find the next element in the list with the right name */
static SrtListAtom *next_element_with_name(SrtListAtom *lst, char *name,
                                           int skip) {
  /* Check that the element exists else nothing to do */
  if (!lst)
    return NULL;

  /* If the skip flag is set to "NO" and points to the right element: return it
   */
  if (!skip && !strcmp(lst->element->name, name))
    return lst;

  /* Go from next to next until the right name is found or end of list */
  lst = lst->next;
  while (lst && strcmp(lst->element->name, name) != 0)
    lst = lst->next;

  /* Return the element (even if NULL) */
  return lst;

} /* END SrtListAtom *next_element_with_name (... ) */

/* ----------------------------------------------------------------------------------
 */

/* Find the previous element in the list with the right name */
static SrtListAtom *previous_element_with_name(SrtListAtom *lst, char *name,
                                               int skip) {
  /* Check that the element exists else nothing to do */
  if (!lst)
    return NULL;

  /* If the skip flag is set to "NO" and points to the right element: return it
   */
  if (!skip && !strcmp(lst->element->name, name))
    return lst;

  /* Go from prev to prev until the right name is found or end of list */
  lst = lst->previous;
  while (lst && strcmp(lst->element->name, name) != 0)
    lst = lst->previous;

  /* Return the element (even if NULL) */
  return lst;

} /* END SrtListAtom* previous_element_with_name(...) */

/* ----------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------
   FUINCTIONS TO INTERPOLATE OR EXTRAPOLATE INFORMATION FROM A SrtDvdObj
   ----------------------------------------------------------------------------------
 */

/* Extracts today from the element labelled "_TODAY" in the SrtDvdObj */
Err srt_f_dvdobj_extracttoday(SrtDvdObj *dvdobj, Ddate *today) {
  SrtListAtom *top;

  if (!dvdobj)
    return serror("No Dividend Curve defined in srt_f_dvdobj_extracttoday");

  /* Go to the First element in the list that is labelled "_TODAY" (not
   * "DVD"...) */
  top = previous_element_with_name(dvdobj->tail, "_TODAY", 0);

  /* Today is the Key of the element */
  *today = (double)(top->element->key);

  /* Return a success message */
  return NULL;

} /* END Err srt_f_fwdobj_extracttoday(...) */

/* ----------------------------------------------------------------------------------
 */

/* Dividend to start_date computed from SrtDvdObj */
Err srt_f_dvdobj_dvd(SrtDvdObj *dvdobj, double start_date, double *dvd) {
  double time;
  double cum_dvd;
  SrtListAtom *top, *bot, *tmplst;
  double date_today;
  Err err = NULL;

  *dvd = 0.0;

  /* Check if the dividend obj is defined */
  if (!dvdobj)
    return serror("no dvd curve defined in srt_f_dvdobj_dvd");
  /* Get date today */
  err = srt_f_dvdobj_extracttoday(dvdobj, &date_today);
  if (err)
    return err;

  /* Go to the Last element in the list that is labelled "DVD" (not "_TODAY"...)
   */
  bot = previous_element_with_name(dvdobj->tail, "DVD", 0);

  /* Go to the First element in the list that is labelled "DVD" (not
   * "_TODAY"...) */
  top = next_element_with_name(dvdobj->head, "DVD", 0);
  if (!top)
    return serror("%s: dvd empty", dvdobj->name);

  /* case where the start date is before the first dividend date */
  if (top->element->key >= (double)(start_date)) {
    if (top->element->key == date_today)
      time = 0.0;
    else
      time = (start_date - date_today) * YEARS_IN_DAY;

    *dvd = exp(-time * top->element->val.dval);
    return NULL;
  }

  /* case where there is at least one dividend date before the start date */
  cum_dvd = 0.0;
  while (top && (top->element->key <= (double)start_date)) {
    if (previous_element_with_name(top, "DVD", 1)) {
      tmplst = previous_element_with_name(top, "DVD", 1);
      time = (top->element->key - tmplst->element->key) * YEARS_IN_DAY;
    } else
      time = (top->element->key - date_today) * YEARS_IN_DAY;

    cum_dvd += top->element->val.dval * time;
    top = next_element_with_name(top, "DVD", 1);
  }

  if (top != NULL) {
    tmplst = previous_element_with_name(top, "DVD", 1);
    time = (start_date - tmplst->element->key) * YEARS_IN_DAY;
    cum_dvd += top->element->val.dval * time;
  }

  else if ((top == NULL) && ((double)start_date > bot->element->key)) {
    time = (start_date - bot->element->key) * YEARS_IN_DAY;
    cum_dvd += bot->element->val.dval * time;
  }

  *dvd = exp(-cum_dvd);
  return NULL;

} /* END Err srt_f_fwdobj_fwd(...) */

/* Computes the implied df (inverse of forward ratio) from start_date to
 * end_date */
Err srt_f_dvdobj_df(SrtDvdObj *dvdobj, double start_date, double end_date,
                    double *df) {
  double start_dvd, end_dvd;
  Err err = NULL;

  /* Compute the forward for the start date */
  err = srt_f_dvdobj_dvd(dvdobj, start_date, &start_dvd);
  if (err)
    return err;

  /* Compute the forward for the end date */
  err = srt_f_dvdobj_dvd(dvdobj, end_date, &end_dvd);
  if (err)
    return err;

  /* The implied discount factor is the ratio of forwards */
  *df = start_dvd / end_dvd;

  /* Return a success message */
  return NULL;

} /* END Err srt_f_dvdobj_df(...) */

/* ----------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------
         INITIALISATION OF A SrtDvdObj FROM RAW DATA (DATES + VALUES )
   ----------------------------------------------------------------------------------
 */

/* Allocate memory and build a SrtDvdObj from a vector of forwards: date / value
 */
Err srt_f_dvdobj_init(Date today,
                      double **dividends, /* dividend curve [date][value] */
                      int num_cols, int num_rows, char *dvd_name,
                      SrtDvdObj **dvdobj) {
  Err err;
  int i;

  /* Check the size of the inputs: num cols */
  if (num_cols != 2)
    return serror("need two columns for the dividend curve: date - value");

  /* Builds an empty SrtDvdObj (it is a double linked list) */
  err = srt_f_dvdobj_create(today, dvd_name, dvdobj);
  if (err) {
    srt_free(*dvdobj);
    return err;
  }

  /* Insert all the points of the forward curve one by one */
  for (i = 0; i < num_rows; i++) {
    err = srt_f_dvdobj_insertpoint(*dvdobj, dividends[0][i], dividends[1][i]);
    if (err) {
      srt_f_dvdobj_delete(dvdobj);
      return err;
    }
  }

  /* Return a success message */
  return NULL;

} /* END Err srt_f_dvdobj_init(...) */

/* ----------------------------------------------------------------------------
 */
