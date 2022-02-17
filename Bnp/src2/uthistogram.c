/* -------------------------------------------------------------------------------------

   FILE NAME : 	  UTHISTOGRAM.c

   PURPOSE:       functions to store data for and create histograms
                  (using the histograms structures)

   -------------------------------------------------------------------------------------
 */
#include "uthistogram.h"
#include "math.h"
#include "utallhdr.h"
#include "utlist.h"

EXTERND SrtHistList *_srt_global_histlistptr = NULL;

/*
static SRT_Boolean  _srt_allocated_histlst = SRT_NO;
*/
static long _srt_hist_path_num = 1;

/** FIVE FUNCTIONS TO BE USED FROM THE OUTSIDE WORLD **/

/* Function to be called when loading the executable file */
Err create_histogram_list(String histo_list_name) {
  Err err;

  if (err = srt_f_crehistlist(&_srt_global_histlistptr, histo_list_name))
    return err;

  return NULL;
}

/* Function to be called when exiting the executable file */
Err destroy_histogram_list() {
  Err err = NULL;

  /* OVE: 10-may-1996
  remove this _srt_allocated_histlst
          if (_srt_allocated_histlst==SRT_YES)
  */
  err = srt_f_frehistlist(&_srt_global_histlistptr);
  /*
  _srt_allocated_histlst=SRT_NO;
  */

  return err;
}

Err srt_f_ininewhist(long num_path) {
  Err err = NULL;

  _srt_hist_path_num = num_path;
  /*
  if (_srt_allocated_histlst==SRT_NO)
          if (err = srt_f_crehistlist(&_srt_global_histlistptr))
                  return(err);

  _srt_allocated_histlst=SRT_YES;
  */

  return err;
}

int srt_f_hislsthas(String name) {
  int has;

  has = srt_f_lsthas(*_srt_global_histlistptr, name, 0);

  return has;
}

Err srt_f_sendhistval(String name, long path, double value) {
  SrtHistData *shd;
  SRT_Boolean bool;
  Err err = NULL;

  strupper(name);
  strip_white_space(name);

  bool = srt_f_lsthas(*_srt_global_histlistptr, name, 0);
  if (bool == SRT_NO)
    shd = srt_f_crehistdata(name, _srt_hist_path_num);
  else {
    shd = srt_f_gethistdata(name);
    if (shd->num_path != _srt_hist_path_num) {
      if (shd->data)
        srt_free(shd->data);
      shd->data = srt_calloc(_srt_hist_path_num, sizeof(double));
      shd->num_path = _srt_hist_path_num;
    }
  }
  if (path > _srt_hist_path_num)
    return serror("Path number is too big !");

  shd->data[path] = value;

  return err;
}

/* ========================================================================= */

/** Allocates space for the referenced SrtHistData and the data inside
        and adds this SrtHistData as a SrtObject into _srt_global_histlistptr
        under this new name  **/
SrtHistData *srt_f_crehistdata(String name, long npath) {
  SrtHistData *data;
  Err err;
  data = srt_malloc(sizeof(SrtHistData));
  if (!data)
    return NULL;
  data->data = srt_calloc(npath, sizeof(double));
  data->num_path = npath;
  if (!data->data) {
    srt_free(data);
    return NULL;
  }

  err = srt_f_addhistdata(name, data);
  if (err)
    return NULL;
  return data;
}

/** Free the space used by SrtHistData and the data inside**/
Err srt_f_frehistdata(SrtHistData *data) {
  if (!data)
    return NULL;
  if (data->data)
    srt_free(data->data);
  srt_free(data);
  return NULL;
}

/** Free the space used by SrtHistData and the data inside (for the list) **/
Err srt_f_frehistval(void *histdata) {
  Err err = NULL;
  SrtHistData *data = (SrtHistData *)histdata;

  err = srt_f_frehistdata(data);
  return err;
}

/** Allocates space for a SrtHistStr
        and returns his address as a *SrtHistStr **/
Err srt_f_crehiststr(SrtHistStr **histo) {
  if (*histo)
    return NULL;
  *histo = srt_malloc(sizeof(SrtHistStr));
  if (!(*histo))
    return serror("Couldn't create SrtHistStr");
  return NULL;
}

/** Free  contents and address of the *SrtHistStr **/
Err srt_f_frehiststr(SrtHistStr *histo) {
  if (!histo)
    return NULL;
  if (histo->abs)
    free_dvector(histo->abs, -1, histo->num_seg);
  histo->abs = NULL;
  if (histo->dens)
    free_ivector(histo->dens, -1, histo->num_seg);
  histo->dens = NULL;
  if (histo->name)
    histo->name[0] = 0;
  histo->num_seg = 0;
  srt_free(histo);
  histo = NULL;
  return NULL;
}

/** Allocates space for a *SrtHistList
        and returns this address as a **SrtHistList  **/
Err srt_f_crehistlist(SrtHistList **hist_list, String hist_list_name) {
  Err err;

  err = srt_f_lstcreate(hist_list, hist_list_name);
  if (err)
    return serror("Couldn't create SrtHistList");
  return NULL;
}

/** Free contents (SrtHistData) and addresses of a SrtHistList	**/
Err srt_f_frehistlist(SrtHistListPtr *hist_list) {
  /*
          SrtHistData *data;
          SrtHistAtom *atom;
  */
  Err err;

  err = srt_f_lstfree(*hist_list, SRT_YES);
  srt_free(*hist_list);
  if (err)
    return (err);
  return NULL;
}

/** Puts the input SrtHistData  into the _srt_global_histlistptr list
        and reference it as name**/
Err srt_f_addhistdata(String name, SrtHistData *data) {
  Err err;
  long ticker;

  err = srt_f_lstins(_srt_global_histlistptr, name, 0.0, OBJ_PTR_HIST,
                     (void *)data, &srt_f_frehistval, &ticker);

  return NULL;
}

/** Picks up in _srt_global_histlistptr the data_base necessary to create the
        histogram referenced with name**/
SrtHistData *srt_f_gethistdata(String name) {
  Err err;
  SrtObject *obj;

  err = srt_f_lstgetobj(*_srt_global_histlistptr, name, 0, &obj);
  if (err)
    return NULL;
  return (SrtHistData *)obj->val.pval;
}

/** Creates the effective histogram referenced as name
                (with num_seg segments in the histogram)
        and puts it into _srt_global_histogram **/
Err srt_f_makehistogram(String name_x, long num_seg_x, String name_y,
                        long num_seg_y,
                        /* OUTPUT */
                        double **histo_x, long *x_size, double **histo_y,
                        long *y_size, double ***histo_values) /* [x][y] */
{
  long i;
  long n_x;
  long n_y;
  double scale_x, min_x, max_x;
  double scale_y, min_y, max_y;
  SrtHistData *hist_data_y = NULL;
  SrtHistData *hist_data_x = NULL;
  SRT_Boolean use_y = SRT_NO;
  SrtHistoType type_x;
  SrtHistoType type_y;
  double dens;

  if (!_srt_global_histlistptr)
    return serror("Histogram list not available ");

  if (!name_x)
    return serror("Need at least one name to make histograms");

  strupper(name_x);
  strip_white_space(name_x);

  if (name_y) {
    strip_white_space(name_y);
    strupper(name_y);
  }

  hist_data_x = srt_f_gethistdata(name_x);
  if (!hist_data_x)
    return serror("Couldn't find SrtHistData for %s", name_x);

  /* If really make an histogram */
  if (num_seg_x) {
    type_x = HISTOGRAM;
  } else
  /* If just output the data */
  {
    type_x = SCATTERED;
  }

  /* Histogram only on x */
  if (!name_y[0]) {
    use_y = SRT_NO;
  } else
  /* Cross histogram */
  {
    use_y = SRT_YES;
    hist_data_y = srt_f_gethistdata(name_y);
    if (!hist_data_y)
      return serror("Couldn't find SrtHistData for %s", name_y);
    /* If really make an histogram */
    if (num_seg_y) {
      type_y = HISTOGRAM;
    } else
    /* If just output the data */
    {
      type_y = SCATTERED;
    }
    /* Check for inconsistancy */
    if ((type_y != type_x)) {
      return serror("Inconsistancy in Histogram type");
    }
    if (hist_data_x->num_path != hist_data_y->num_path) {
      return serror("Number of paths do not match: cannot create histogram");
    }
  }

  if (type_x == HISTOGRAM) {
    *x_size = num_seg_x;
    if (use_y == SRT_NO)
      *y_size = 1;
    else
      *y_size = num_seg_y;

    dens = 1.0 / hist_data_x->num_path;

    /* Gets max and min to make sure that there is more than one value for x (or
     * y) */
    max_x = hist_data_x->data[0];
    min_x = max_x;
    for (i = 0; i < hist_data_x->num_path; i++) {
      if (hist_data_x->data[i] > max_x)
        max_x = hist_data_x->data[i];
      else if (hist_data_x->data[i] < min_x)
        min_x = hist_data_x->data[i];
    }
    scale_x = (max_x - min_x) / (num_seg_x);
    if (scale_x == 0.0)
      *x_size = 1;

    if (use_y == SRT_YES) {
      dens /= hist_data_y->num_path;
      max_y = hist_data_y->data[0];
      min_y = max_y;
      for (i = 0; i < hist_data_y->num_path; i++) {
        if (hist_data_y->data[i] > max_y)
          max_y = hist_data_y->data[i];
        else if (hist_data_y->data[i] < min_y)
          min_y = hist_data_y->data[i];
      }
      scale_y = (max_y - min_y) / (num_seg_y);
      if (scale_y == 0.0)
        *y_size = 1;
    }

    /* Memory allocation */
    (*histo_x) = dvector(0, *x_size - 1);
    (*histo_y) = dvector(0, *y_size - 1);
    (*histo_values) = dmatrix(0, *x_size - 1, 0, *y_size - 1);

    if (!(*histo_x) || !(*histo_y) || !(*histo_values)) {
      if (*histo_x)
        free_dvector((*histo_x), 0, *x_size - 1);
      if (*histo_y)
        free_dvector((*histo_y), 0, *y_size - 1);
      if (*histo_values)
        free_dmatrix((*histo_values), 0, *x_size - 1, 0, *y_size - 1);
      return serror("Memory Allocation Failure (1) in srt_f_makehistogram");
    }

    /* Sets all values to 0 */
    memset(&(*histo_values)[0][0], 0, (*x_size) * (*y_size) * sizeof(double));

    /* Sets the points into the right segment */
    for (i = 0; i < hist_data_x->num_path; i++) {

      if (scale_x != 0.0)
        n_x = (long)floor((hist_data_x->data[i] - min_x) / scale_x);
      else
        n_x = 0;
      if (use_y == SRT_YES) {
        if (scale_y != 0.0)
          n_y = (long)floor((hist_data_y->data[i] - min_y) / scale_y);
        else
          n_y = 0;

      } else {
        n_y = 0;
      }
      /* Safety procedure for when n_x corresponds to the max */
      if (n_x >= *x_size)
        n_x = *x_size - 1;
      if (n_y >= *y_size)
        n_y = *y_size - 1;

      (*histo_values)[n_x][n_y] += dens;
    }

    /* Sets the Histo_x and histo_y values */
    for (i = 0; i < *x_size; i++) {
      (*histo_x)[i] = min_x + (double)(i + 0.5) * scale_x;
    }
    if (use_y == SRT_YES) {
      for (i = 0; i < *y_size; i++)
        (*histo_y)[i] = min_y + (double)(i + 0.5) * scale_y;
    } else {
      (*histo_y)[0] = 0.0;
    }

  } /* END if (type_x == HISTOGRAM) */
  else if (type_x == SCATTERED) {
    *x_size = hist_data_x->num_path;
    if (use_y == SRT_NO)
      *y_size = 1; /* returns only one column */
    else
      *y_size = 2; /* returns two columns */

    /* Memory allocation */
    (*histo_x) = dvector(0, *x_size - 1);
    (*histo_y) = dvector(0, *y_size - 1);
    (*histo_values) = dmatrix(0, *x_size - 1, 0, *y_size - 1);

    if (!(*histo_x) || !(*histo_y) || !(*histo_values)) {
      if (*histo_x)
        free_dvector((*histo_x), 0, *x_size - 1);
      if (*histo_y)
        free_dvector((*histo_y), 0, *y_size - 1);
      if (*histo_values)
        free_dmatrix((*histo_values), 0, *x_size - 1, 0, *y_size - 1);
      return serror("Memory Allocation Failure (2) in srt_f_makehistogram");
    }

    /* Sets all values to 0 */
    memset(&(*histo_values)[0][0], 0, (*x_size) * (*y_size) * sizeof(double));

    /* Copy all the values coming from data and sets histo_x to number of path
     */
    for (i = 0; i < hist_data_x->num_path; i++) {
      (*histo_values)[i][0] = hist_data_x->data[i];
      (*histo_x)[i] = i;
    }
    (*histo_y)[0] = 0.0;

    if (use_y == SRT_YES) {
      (*histo_y)[1] = 0.0;
      for (i = 0; i < hist_data_y->num_path; i++)
        (*histo_values)[i][1] = hist_data_y->data[i];
    }
  }

  return NULL;
}
