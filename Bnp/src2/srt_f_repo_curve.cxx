
#include "srt_h_repo_obj.h"
#include "swp_h_all.h"
#include "utallhdr.h"

Err srt_f_init_EQRepoCurve(Date today, char *ccy, double **repo_curve,
                           long ncols, long nrows, String repo_crv_name) {
  SrtRepoObj *repo_obj;
  Err err = NULL;
  SrtCurveList *curve_list;

  curve_list = get_curve_list();
  if (curve_list == NULL)
    return serror("no curve list defined");

  err = srt_f_repo_obj_init(today, repo_curve, ncols, nrows, repo_crv_name,
                            &repo_obj);
  if (err)
    return err;

  err = swp_f_addcurvetolist(curve_list, repo_crv_name, "REPO_CURVE", NULL, ccy,
                             today, 0, NULL, (void *)(repo_obj));

  if (err)
    return serror("can not initialise the repo curve");

  return NULL;
}