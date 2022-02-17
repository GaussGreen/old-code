/* -------------------------------------------------------------------------

        FILE NAME	: srt_f_correlation_list.cxx

        PURPOSE		: use a SrtList to store a TermStructure of correlation
   matrixes gives a range of functions that allow to use this list properly
   ------------------------------------------------------------------------- */

#include "math.h"
#include "srt_h_all.h"

#define SWAP(X, Y)                                                             \
  {                                                                            \
    void *_VPTR;                                                               \
    _VPTR = X;                                                                 \
    X = Y;                                                                     \
    Y = _VPTR;                                                                 \
  }

/*  ----------------------------------------------------------------------
        Two useful static functions defined at the bottom of the file
        ----------------------------------------------------------------------
 */
static Err sort_underlyings(int ndates, int ncorr, double ***correl, int n_und,
                            String **und_names, String **sort_und_names);

static Err sort_names(int n_names, String **names);

/* -----------------------------------------------------------------------
   Creates (allocates space for) any  SrtCorrLstPtr and gives it a name...
   ----------------------------------------------------------------------- */
Err srt_f_corrlstcreate(SrtCorrLstPtr *cls, String list_name) {
  Err err = NULL;
  if (err = srt_f_lstcreate(cls, list_name))
    return err;

  return NULL;
}

/* -----------------------------------------------------------------------
   Deletes (frees the space for) any SrtCorrLstPtr
   ----------------------------------------------------------------------- */
Err srt_f_corrlstdelete(SrtCorrLstPtr *cls) {
  Err err;

  if (!*cls)
    return NULL;

  /* Frees everything: even the SrtCorrVal and their contents */
  err = srt_f_lstfree(*cls, SRT_YES);
  if (err)
    return serror("Error in freeing correlation list");

  srt_free(*cls);

  *cls = NULL;

  return NULL;
}

/* -----------------------------------------------------------------------
   Frees a SrtCorrLstVal (the pval of an element of the SrtCorrLst list
   ----------------------------------------------------------------------- */
Err srt_f_corrvalfree(void *corrvalptr) {
  Err err = NULL;
  int i;
  SrtCorrLstVal *corval = (SrtCorrLstVal *)corrvalptr;

  if (!corval)
    return NULL;

  if (corval->correl)
    free_dmatrix(corval->correl, 0, corval->nund - 1, 0, corval->nund - 1);
  corval->correl = NULL;

  if (corval->coeff)
    free_dmatrix(corval->coeff, 0, corval->nund - 1, 0, corval->nund - 1);
  corval->coeff = NULL;

  if (corval->und_names) {
    for (i = 0; i < corval->nund; i++) {
      if (corval->und_names[i])
        free(corval->und_names[i]);
      corval->und_names[i] = NULL;
    }
    free_svector(corval->und_names, 0, corval->nund - 1);
  }
  corval->und_names = NULL;

  srt_free(corval);
  corval = NULL;

  return err;
}

/* -----------------------------------------------------------------------
   Initialises (fills in) the SrtCorrLstPtr the_corr_list according to
   data taken from spreadsheet...
                und_names is a two column array: [0..1][0..ncorr-1]
                correl is a full matrix [0..ndates-1][0..ncorr-1]        , which
   is supposed to store a full term structure of correlation matrixes
   ----------------------------------------------------------------------- */

Err srt_f_init_Corr_TermStruct(int ndates, int ncorr, double **correl,
                               double *dates, String **und_names,
                               SrtCorrLstPtr *the_corr_list) {
  SrtCorrLstVal *corrval;
  Err err = NULL;
  int cur_date_ind, date_index, corr_index;
  int n_und;
  double d_n_und;
  double today;
  SrtUndPtr und;
  String *sort_und_name;
  int i, j;
  long ticker;

  /* Clear the contents of the correlation list */
  /* To clear it we have to destroy the correlation matrix and recreate it */

  if (the_corr_list) {
    err = srt_f_corrlstdelete(the_corr_list);
    if (err)
      return err;
  }

  err = srt_f_corrlstcreate(the_corr_list, "Init Corr Matrix");
  if (err)
    return err;

  /* If n_und = n        , ncorr = n(n-1)/2 */
  d_n_und = 0.5 * (sqrt(8.0 * (double)ncorr + 1.0) + 1.0);
  n_und = (int)(d_n_und);
  if ((double)n_und != d_n_und)
    return serror("Wrong number of correlations input in the matrix");

  /* Memory allocation for sort_und_names */
  sort_und_name = svector(0, n_und - 1);

  /* Sorts underlyings and modifies the correlation matrix accordingly */
  if (err = sort_underlyings(ndates, ncorr, &correl, n_und, und_names,
                             &sort_und_name))
    return err;

  /* DO NOT DO THIS FOR FX STOCH RATES : corr is needed before UND initilisation
                  Makes sure underlyings are initialised
          for (i=0; i<n_und;i++)
          {
                  und = lookup_und(sort_und_name[i]);
                  if (!und)
                          return serror("Underlying %s not found in
     init_Corr_TermStruct !"        , sort_und_name[i]);
          }
  */

  /* Gets today */

  /*	TERRIBLE CODE !!!!!!!: we need one underlying to be initalised
          to be able to init the corr matrix... */

  today = -1;
  i = 0;
  while (today < 0 && i < n_und) {
    und = lookup_und(sort_und_name[i]);
    today = (double)get_today_from_underlying(und);
    i++;
  }

  if (today < 0 && ndates > 1) {
    /* we have a problem only if there is a Term Structure */
    return serror("Cannot initialise the Corr Term Structure when no "
                  "underlying has been initialised yet");
  }

  /* Make sure there is more than today */
  cur_date_ind = 0;
  while (dates[cur_date_ind] < today) {
    cur_date_ind++;
    if (cur_date_ind >= ndates)
      return serror("init_Corr_TermStruct err: no data after today !!");
  }

  /* Loop on all dates after today */
  for (date_index = cur_date_ind; date_index < ndates; date_index++) {
    /* Creation of a new SrtCorrLstVal */
    corrval = (SrtCorrLstVal *)srt_malloc(sizeof(SrtCorrLstVal));
    if (!corrval)
      return serror("Allocation failure (1) in srt_f_init_Corr_TermStruct");

    /* Transfer minimal information */
    corrval->date = dates[date_index];
    corrval->time = (dates[date_index] - today) * YEARS_IN_DAY;
    corrval->ncorrel = ncorr;
    corrval->nund = n_und;

    /* Allocate space for the underlying names */
    corrval->und_names = svector(0, n_und - 1);
    if (!corrval->und_names)
      return serror("Allocation failure (2) in srt_f_init_Corr_TermStruct");

    /* Copy underlying names */
    for (i = 0; i < n_und; i++) {
      corrval->und_names[i] =
          (char *)malloc((strlen(sort_und_name[i]) + 1) * sizeof(char));
      if (!corrval->und_names[i])
        return serror("Allocation failure (3) in srt_f_init_Corr_TermStruct");
      strncpy(corrval->und_names[i], sort_und_name[i],
              strlen(sort_und_name[i]));
      corrval->und_names[i][strlen(sort_und_name[i])] = '\0';
    }

    /* Allocate space for the correlation matrix */
    corrval->correl = dmatrix(0, n_und - 1, 0, n_und - 1);
    if (!corrval->correl)
      return serror("Allocation failure (4) in srt_f_init_Corr_TermStruct");

    corrval->coeff = dmatrix(0, n_und - 1, 0, n_und - 1);
    if (!corrval->coeff)
      return serror("Allocation failure (5) in srt_f_init_Corr_TermStruct");

    /* Sets diagonal values to 1.0 */
    for (i = 0; i < n_und; i++)
      corrval->correl[i][i] = 1.0;

    /* Transfer remaining correlation information (using corr_index        ,
                    knowing that underlying names have been sorted accordingly)
     */
    corr_index = 0;
    for (i = 0; i < n_und - 1; i++) /* Do not need the last one */
    {
      for (j = i + 1; j < n_und; j++) {
        corrval->correl[i][j] = correl[corr_index][date_index];
        corrval->correl[j][i] = corrval->correl[i][j];
        corr_index++;
      }
    }

    /* Checks whether the correlation matrix is a real correlation matrix */
    if (1) {
      err = compute_eigen_from_correl(corrval->correl, n_und, corrval->coeff);

    } else {
      err = compute_coeff_from_correl(corrval->correl, n_und, corrval->coeff);
    }
    if (err) {
      sprintf(err, "%s at date %.0f", err, corrval->date);
      return err;
    }

    /* The pval value of the object is set to point to corrval */
    srt_f_lstins(*the_corr_list, "CorrAtom", corrval->date,
                 OBJ_PTR_CorrelStruct, (void *)corrval, &srt_f_corrvalfree,
                 &ticker);

  } /* END loop on date_index = cur_date_ind; date_index< ndates */

  /* Free the svector sort_und_names */
  free_svector(sort_und_name, 0, n_und - 1);

  /* Return a success message: err == NULL */
  return err;

} /* END of srt_f_init_Corr_TermStruct(...) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   Initialises and fills in a local SrtCorrLst        , extracting information
   from the SrtCorrLstPtr big_corr_list according to
   the underlyings used in the deal
                und_names is a vector with the underlying names: [0..n_und]
   Please note that the local correlation matrix will be defined according
   to the underlying order as input in und_names
   ----------------------------------------------------------------------- */

Err srt_f_make_deal_corrlist(String *und_names, int n_und, String list_name,
                             SrtCorrLstPtr big_corr_list,
                             SrtCorrLstPtr *deal_corr_list) {
  Err err = NULL;
  SrtLst *ls, *ls1;
  SrtCorrLstVal *corrval, *corrval1;
  long ticker;
  int i, j;
  int *und_index;

  /* Checks that the correlation list is properly created */
  if (!big_corr_list)
    return serror("Global correlation list not properly initialised in "
                  "make_deal_corrlist");

  /* Creates (allocate space for) the local SrtCorrLst used for the deal  */
  if (err = srt_f_corrlstcreate(deal_corr_list, list_name))
    return err;
  ls1 = (*deal_corr_list)->head;

  /* Checks that the main list is properly filled in */
  ls = big_corr_list->head;
  if (!ls->element)
    return serror("Global correlation list empty in make_deal_corrlist: cannot "
                  "compute correlation");
  corrval = (SrtCorrLstVal *)ls->element->val.pval;

  /* Memory allocation for indexes that will be used to refer to the
   * underlyings*/
  und_index = ivector(0, n_und);

  /* Goes through the whole und name list and extracts the indexes */
  for (i = 0; i < n_und; i++) {
    j = 0;
    while (strcmp(corrval->und_names[j], und_names[i]) != 0) {
      j++;
      if (j >= corrval->nund) {
        free_ivector(und_index, 0, n_und);
        return serror("Could not find und %s in local deal corr", und_names[i]);
      }
    }
    und_index[i] = j;
  } /* END loop on number of underlyings */

  /* Duplicates the correlation information needed        , passing through the
   * list*/
  while (ls != NULL) {
    corrval = (SrtCorrLstVal *)ls->element->val.pval;

    /* Creation of a new SrtCorrLstVal */
    corrval1 = (SrtCorrLstVal *)srt_malloc(sizeof(SrtCorrLstVal));
    if (!corrval1) {
      free_ivector(und_index, 0, n_und);
      return serror("Allocation failure (1) in srt_f_make_deal_corrlist");
    }

    /* Transfer minimal information */
    corrval1->date = corrval->date;
    corrval1->time = corrval->time;
    corrval1->nund = n_und;
    corrval1->ncorrel = (int)(0.5 * n_und * (n_und - 1));

    /* Memory allocation */
    corrval1->und_names = svector(0, n_und - 1);
    if (!corrval1->und_names) {
      free_ivector(und_index, 0, n_und);
      return serror("Allocation failure (2) in srt_f_init_Corr_TermStruct");
    }
    corrval1->correl = dmatrix(0, n_und - 1, 0, n_und - 1);
    if (!corrval1->correl) {
      free_ivector(und_index, 0, n_und);
      return serror("Allocation failure (3) in srt_f_init_Corr_TermStruct");
    }
    corrval1->coeff = dmatrix(0, n_und - 1, 0, n_und - 1);
    if (!corrval1->coeff) {
      free_ivector(und_index, 0, n_und);
      return serror("Allocation failure (4) in srt_f_init_Corr_TermStruct");
    }

    /* Copy underlying names */
    for (i = 0; i < n_und; i++) {
      corrval1->und_names[i] =
          (char *)malloc((strlen(und_names[i]) + 1) * sizeof(char));
      if (!corrval->und_names[i]) {
        free_ivector(und_index, 0, n_und);
        return serror("Allocation failure (5) in srt_f_init_Corr_TermStruct");
      }
      strncpy(corrval1->und_names[i], und_names[i], strlen(und_names[i]));
      corrval1->und_names[i][strlen(und_names[i])] = '\0';
    }

    /* Copies underlying names and fills in correlation matrix */
    /*	memcpy(corrval1->und_names        , und_names        ,
     * n_und*sizeof(char*));*/

    for (i = 0; i < n_und; i++) {
      for (j = i; j < n_und; j++) {
        corrval1->correl[i][j] = corrval->correl[und_index[i]][und_index[j]];
        corrval1->correl[j][i] = corrval1->correl[i][j];
      }
    }

    /* Computes the linear coefficients used to correlate the numbers */
    if (1) {
      err = compute_eigen_from_correl(corrval1->correl, n_und, corrval1->coeff);
    } else {
      err = compute_coeff_from_correl(corrval1->correl, n_und, corrval1->coeff);
    }
    if (err) {
      free_ivector(und_index, 0, n_und);
      return err;
    }

    /* The pval value of the object is set to point to corrval1 */
    srt_f_lstins(*deal_corr_list, "CorrAtom", corrval1->date,
                 OBJ_PTR_CorrelStruct, corrval1, &srt_f_corrvalfree, &ticker);

    /* Goes to next element in the main list */
    ls = ls->next;

  } /* END of while (ls!= NULL) */

  /* Frees memory allocated for indexes */
  free_ivector(und_index, 0, n_und);

  /* Return a success string (err = NULL) */
  return err;

} /* END of srt_f_make_deal_corrlist(...) */

/* ---------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   Gets the correlation for two underlyings at a certain time (in yrs
   from tdy)        , according to the initialisation of the SrtCorrLstPtr
   correlation list There is no interpolation: the correlation is supposed to be
   constant between two dates
   ----------------------------------------------------------------------- */

Err srt_f_get_corr_from_CorrList(SrtCorrLstPtr cls, String und_name1,
                                 String und_name2, double time, double *value) {
  SrtLst *ls;
  Err err = NULL;
  int i, j;
  SrtCorrLstVal *corrval;

  /* Check if cls has been initialised */
  if (!cls)
    return serror("CorrList not initialised in srt_f_get_corr_from_CorrList");

  /* Moves along the cls until finds a corrval->time > time */
  ls = cls->head;
  while ((ls != NULL) &&
         (((SrtCorrLstVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = cls->tail;

  /* Selects the SrtCorrLstVal element associated to the date (time)  */
  corrval = (SrtCorrLstVal *)ls->element->val.pval;
  if (!corrval->correl)
    return serror("Correlation TermStructure void for time : %f", time);

  /* Finds the underlying number in the corrval->und_names[] for first und */
  i = 0;
  while (strcmp(und_name1, corrval->und_names[i]) != 0) {
    i++;
    if (i >= corrval->nund)
      return serror("Could not find und in correlation matrix: %s", und_name1);
  }

  /* Finds the underlying number in the corrval->und_names[] for second und */
  j = 0;
  while (strcmp(und_name2, corrval->und_names[j]) != 0) {
    j++;
    if (j >= corrval->nund)
      return serror("Could not find und in correlation matrix: %s", und_name2);
  }

  /* Picks up the corresponding value in the corelation matrix */
  *value = corrval->correl[i][j];

  /* Returns a success message: err = NULL */
  return err;

} /* END of srt_f_get_corr_from_CorrList(...) */

/* ====================================================================== */

/* ----------------------------------------------------------------------
   Sorts the underlying names by alphabetical order        , and modifies the
   associated correlation matrix accordingly.
   Space has to be allocated before for *sort_und_names        ,
   The function fills in sort_und_names[...] and modifies **correl
   ---------------------------------------------------------------------- */

static Err sort_underlyings(int ndates, int ncorr, double ***correl, int n_und,
                            String **und_names, String **sort_und_names) {
  Err err = NULL;
  int i, j, k;
  double rubbish;
  char rubbishstr[256];

  /* First sorts so that UND_NAME[...][0] < UND_NAME[...][1] */
  for (i = 0; i < ncorr; i++) {
    if (strcmp(und_names[i][0], und_names[i][1]) > 0) {
      /* OVE: SWAP is dangerous with NR allocation procedures when freed */
      /* Let us hope that all strings have been allocated bigger than necessary
       */
      strncpy(rubbishstr, und_names[i][1], strlen(und_names[i][1]) + 1);
      strncpy(und_names[i][1], und_names[i][0], strlen(und_names[i][0]) + 1);
      strncpy(und_names[i][0], rubbishstr, strlen(rubbishstr) + 1);
    }
    /* Return an error message if same underlying name repeated in matrix */
    if (strcmp(und_names[i][0], und_names[i][1]) == 0)
      return serror("Duplicate underlying name in correlation matrix: %s",
                    und_names[i][0]);
  }

  /* Start loop on all und_names[...] */
  for (i = 0; i < ncorr; i++) {
    for (j = i + 1; j < ncorr; j++) {
      /* Sort first name first (when copying        , copy the '\0' as well) */
      if (strcmp(und_names[i][0], und_names[j][0]) > 0) {
        /*Swaps both names */
        strncpy(rubbishstr, und_names[i][0], strlen(und_names[i][0]) + 1);

        strncpy(und_names[i][0], und_names[j][0], strlen(und_names[j][0]) + 1);

        strncpy(und_names[j][0], rubbishstr, strlen(rubbishstr) + 1);

        strncpy(rubbishstr, und_names[i][1], strlen(und_names[i][1]) + 1);

        strncpy(und_names[i][1], und_names[j][1], strlen(und_names[j][1]) + 1);

        strncpy(und_names[j][1], rubbishstr, strlen(rubbishstr) + 1);

        /* Swaps the entire correlation line */
        for (k = 0; k < ndates; k++) {

          rubbish = (*correl)[i][k];
          (*correl)[i][k] = (*correl)[j][k];
          (*correl)[j][k] = rubbish;
        }
      } /* END if (und_names[i][0] > und_names[j][0]) */
      else
          /* If first name are equal        , then sort or check second one*/
          if (strcmp(und_names[i][0], und_names[j][0]) == 0) {
        if (strcmp(und_names[i][1], und_names[j][1]) > 0) {
          /*Swaps both names (when copying        , copy the '\0' as well)  */
          strncpy(rubbishstr, und_names[i][1], strlen(und_names[i][1]) + 1);
          strncpy(und_names[i][1], und_names[j][1],
                  strlen(und_names[j][1]) + 1);
          strncpy(und_names[j][1], rubbishstr, strlen(rubbishstr) + 1);

          /* Swaps the entire correlation line */
          for (k = 0; k < ndates; k++) {
            rubbish = (*correl)[i][k];
            (*correl)[i][k] = (*correl)[j][k];
            (*correl)[j][k] = rubbish;
          }
        } else if (strcmp(und_names[i][1], und_names[j][1]) == 0)
          return serror("Duplicate correlation couple in matrix :%s-%s",
                        und_names[i][0], und_names[i][1]);
      } /* END if first names are equal */
    }   /* END for (j=i+1; j<ncorr:j++) loop*/
  }     /* END for (i=0; i<ncorr:i++) loop */

  /* Copy the names according to the alphabetical order of the underlyings ,
     assuming that the names are well defined (i.e. UND1 UND2; UND1 UND3; ...*/
  (*sort_und_names)[0] = und_names[0][0];
  j = 0;
  for (i = 0; i < ncorr; i++) {
    if (strcmp(und_names[i][1], (*sort_und_names)[j]) > 0) {
      j++;
      if (j < n_und) {
        (*sort_und_names)[j] = und_names[i][1];
      } else
        return serror("Bad (number of) names in correlation matrix");
    }
  } /* END of loop on second name */

  /* Check that the names are just what they should be */
  for (i = 0; i < ncorr; i++) {
    for (j = 0; j < 2; j++) {
      k = 0;
      while (strcmp((*sort_und_names)[k], und_names[i][j]) != 0) {
        k++;
        if (k >= n_und)
          return serror("Bad (number of) names in correlation matrix");
      }
    }
  } /* END of check that the names are properly defined */

  return err;
} /* END of sort_underlyings */

/* ====================================================================== */

/* ----------------------------------------------------------------------
   Sorts names by alphabetical order
   Modifies names[0..n_names-1]
   ---------------------------------------------------------------------- */

static Err sort_names(int n_names, String **names) {
  int i, j;
  Err err = NULL;
  char rubbishstr[256];

  for (i = 0; i < n_names; i++) {
    for (j = i + 1; j < n_names; j++) {
      if (strcmp((*names)[i], (*names)[j]) > 0) {
        strncpy(rubbishstr, (*names)[i], strlen((*names)[i]));
        strncpy((*names)[i], (*names)[j], strlen((*names)[j]));
        strncpy((*names)[j], rubbishstr, strlen(rubbishstr));
      }
    }
  }
  return err;
}
