#include "math.h"
#include "srt_h_repo_obj.h"
#include "utallhdr.h"

static SrtListAtom *next_element_with_name(SrtListAtom *lst, char *name,
                                           int skip) {
  if (!lst)
    return NULL;

  if (!skip && !strcmp(lst->element->name, name))
    return lst;

  lst = lst->next;
  while (lst && strcmp(lst->element->name, name) != 0)
    lst = lst->next;

  return lst;
}

static SrtListAtom *previous_element_with_name(SrtListAtom *lst, char *name,
                                               int skip) {
  if (!lst)
    return NULL;

  if (!skip && !strcmp(lst->element->name, name))
    return lst;

  lst = lst->previous;
  while (lst && strcmp(lst->element->name, name) != 0)
    lst = lst->previous;

  return lst;
}

Err srt_f_repo_obj_create(Date today, char *repo_crv_name,
                          SrtRepoObj **repo_obj) {
  Err err = NULL;
  char list_name[SRTBUFSZ];
  double repo_today;
  long ticker;

  sprintf(list_name, "%s_SrtRepoObject", repo_crv_name);

  err = srt_f_lstcreate(repo_obj, list_name);
  if (err)
    return serror("error in creating the %s SrtRepoObject", repo_crv_name);

  /* insert a fake point cotaining ony today in the SrtList = SrtRepoObj */
  repo_today = 0.0;
  err = srt_f_lstins(*repo_obj, "_TODAY", (double)today, OBJ_DOUBLE,
                     (void *)(&repo_today), NULL, &ticker);
  if (err)
    return err;

  if (err)
    return serror("Error in inserting today in the %s SrtRepoObj",
                  repo_crv_name);

  return err;
}
Err srt_f_repo_obj_delete(SrtRepoObj **repo_obj) {
  if (srt_f_lstfree(*repo_obj, SRT_YES) != NULL)
    return serror("error in srt_f_lstfree");

  srt_free(*repo_obj);

  return NULL;
}

Err srt_f_repo_obj_insert_point(SrtRepoObj *repo_obj, Ddate date, double repo) {
  Err err = NULL;
  long ticker;

  err = srt_f_lstins(repo_obj, "REPO", (double)date, OBJ_DOUBLE,
                     (void *)(&repo), NULL, &ticker);
  if (err)
    return err;

  return NULL;
}

Err srt_f_repo_obj_init(Date today, double **repos, int num_cols, int num_rows,
                        char *repo_crv_name, SrtRepoObj **repo_obj) {
  Err err = NULL;
  int i;

  if (num_cols != 2)
    return serror("need two columns for the repo curve");

  err = srt_f_repo_obj_create(today, repo_crv_name, repo_obj);
  if (err) {
    srt_free(*repo_obj);
    return err;
  }

  for (i = 0; i < num_rows; i++) {
    err = srt_f_repo_obj_insert_point(*repo_obj, repos[0][i], repos[1][i]);
    if (err) {
      srt_f_repo_obj_delete(repo_obj);
      return err;
    }
  }

  return NULL;
}

Err srt_f_repo_obj_extracttoday(SrtRepoObj *repo_obj, Ddate *today) {
  SrtListAtom *top;

  if (!repo_obj)
    return serror("no repo curve defined in srt_f_repo_obj_extracttoday");

  top = previous_element_with_name(repo_obj->tail, "_TODAY", 0);

  *today = (double)(top->element->key);

  return NULL;

} /* END Err srt_f_fwdobj_extracttoday(...) */

Err srt_f_repo_obj_repo(SrtRepoObj *repo_obj, double start_date, double *repo) {
  double time;
  double cum_repo;
  SrtListAtom *top, *bot, *tmplst;
  double date_today;
  Err err = NULL;

  *repo = 0.0;

  /* Check if the repo obj is defined */
  if (!repo_obj)
    return serror("no repo curve defined in srt_f_repo_obj_repo");
  /* Get date today */
  err = srt_f_repo_obj_extracttoday(repo_obj, &date_today);
  if (err)
    return err;

  /* Go to the Last element in the list that is labelled "DVD" (not "_TODAY"...)
   */
  bot = previous_element_with_name(repo_obj->tail, "REPO", 0);

  /* Go to the First element in the list that is labelled "DVD" (not
   * "_TODAY"...) */
  top = next_element_with_name(repo_obj->head, "REPO", 0);
  if (!top)
    return serror("%s: repo empty", repo_obj->name);

  if (top->element->key >= (double)(start_date)) {
    if (top->element->key == date_today)
      time = 0.0;
    else
      time = (start_date - date_today) * YEARS_IN_DAY;

    *repo = exp(-time * top->element->val.dval);
    return NULL;
  }

  cum_repo = 0.0;
  while ((top) && (top->element->key <= (double)(start_date))) {
    if (previous_element_with_name(top, "REPO", 1)) {
      tmplst = previous_element_with_name(top, "REPO", 1);
      time = (top->element->key - tmplst->element->key) * YEARS_IN_DAY;
    }

    else
      time = (top->element->key - date_today) * YEARS_IN_DAY;

    cum_repo += top->element->val.dval * time;

    top = next_element_with_name(top, "REPO", 1);
  }

  if (top != NULL) {
    tmplst = previous_element_with_name(top, "REPO", 1);
    time = (start_date - tmplst->element->key) * YEARS_IN_DAY;
    cum_repo += top->element->val.dval * time;
  } else if ((top == NULL) && ((double)start_date > bot->element->key)) {
    time = (start_date - bot->element->key) * YEARS_IN_DAY;
    cum_repo += bot->element->val.dval * time;
  }

  *repo = exp(-cum_repo);
  return NULL;

} /* END Err srt_f_repo_obj_repo(...) */
