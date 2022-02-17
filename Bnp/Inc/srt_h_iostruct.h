/* SRT_H_IOSTRUCT.H */

#ifndef IOSTRUCT_H
#define IOSTRUCT_H

/* -------------------------------------------------------------------------- */

typedef SrtList SrtIOStruct;

/* SrtIOStruct is a double linked list.
   Each element is contains a request for a computation (price  , tau shift  ,
   sigma shift  ,...) and the result of the computation.
   IO <=> Input/Output */

/* -------------------------------------------------------------------------- */

typedef enum SRT_IO_TYPE {
  IO_PREMIUM,
  IO_RATE_SHIFT,
  IO_SIGMA_SHIFT,
  IO_TAU_SHIFT,
  IO_MAX,
  IO_EXPSQPAY,
  IO_STDEV,
  IO_COLPVS
} SRT_IO_TYPE;

/* -------------------------------------------------------------------------- */

typedef enum SRT_SHIFT_TYPE {
  SH_NONE,
  SH_ABSOLUTE, /* the shift value is absolute */
  SH_RELATIVE  /* the shift value is a % of the value shifted */

} SRT_SHIFT_TYPE;

/* -------------------------------------------------------------------------- */

typedef struct SrtIOVal {
  char id_name[32];
  int type;
  Date bucket_start; /* for tau and sigma shifts only */
  Date bucket_end;
  double shift_value;          /* for rate  , sigma or tau shifts */
  int shift_type;              /* either relative or absolute */
  SRT_Boolean done;            /* wether the job has been done or not */
  char computation_origin[32]; /* name of the function that computed 'result */
  long lval;
  double dval; /* value computed */
  void *pval;  /* for the vectors */
} SrtIOVal;

/* ---------------------- type2string function ------------------------------ */

Err srt_f_IOtype2string(int type, char name[32]);

/* -> convert the type number to a string */

/* -------------------------- name encoding function ------------------------ */

Err srt_f_IOencodename(char id_name[32], int type, int bucket_no,
                       double shift_value, int shift_type, char name[32]);

/* -> build the name of an element using its field values.
        Each element must be reached thanks to its name.
        Each name must be different.				*/

/* ------------------------ Show function ----------------------------------- */

Err srt_f_IOvalshow(SrtIOVal *ioval, FILE *file_out);

/* -> print out the contents of SrtIOVal '*ioval into file 'file_out */

/* ----------------------- Free function for IOVal -------------------------- */

Err srt_f_IOvalfree(void *ioval);

/* -> free SrtIOVal '*ioval */

/* ---------------------- Free function for IOSTruct ------------------------ */

Err srt_f_IOstructfree(SrtIOStruct **l);

/* -> free all the list including all the objects pointed
This function is strictly equivalent srt_lstfree(.  ,SRT_YES) */

/* ---------------------------- Conversion function ------------------------- */

SrtObjVal *srt_f_IOval2obj(SrtIOVal *ioval);

/* -> convert a SrtIOVal object into a SrtObjVal */

/* ------------------------- Generic set function --------------------------- */

Err srt_f_IOstructset(SrtIOStruct *l, char id_name[32], int type, int bucket_no,
                      Date bucket_start, Date bucket_end, double shift_value,
                      int shift_type, SRT_Boolean done, char comp_origin[32],
                      long lvalue, double dvalue, void *pvalue);

/* -> put SrtIOVal '*ioval into the list '*l */

/* ------------------------- Generic get function --------------------------- */

Err srt_f_IOstructget(SrtIOStruct l, char id_name[32], int type, int bucket_no,
                      double shift_value, int shift_type, SrtIOVal **ioval_pp);

/* -------------------------------------------------------------------------- */
/* -------------------------- Premium functions ----------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetpremiumval(SrtIOStruct l, double *premium);

/* -> return as '*premium the premium from the list 'l */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetpremium(SrtIOStruct l, SrtIOVal **ioval);

/* -> return as '*ioval the SrtIOVal corresponding to the premium
        in list 'l */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetpremium(SrtIOStruct *l, SRT_Boolean done, double result,
                             char comp_origin[32]);

/* -> set the list cell corresponding to "PREMIUM" with the parameters */

/* -------------------------------------------------------------------------- */
/* ------------------- Standard Deviation functions ------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetstdevval(SrtIOStruct l, double *stdev);

/* -> return as '*stdev the standard deviation from the list 'l */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetstdev(SrtIOStruct l, SrtIOVal **ioval);

/* -> return as '*ioval the SrtIOVal corresponding to the standard deviation
        in list 'l */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetstdev(SrtIOStruct *l, SRT_Boolean done, double result,
                           char comp_origin[32]);

/* -> set the list cell corresponding to "STDEV" with the parameters */

Err srt_f_IOstructsetcolpvs(SrtIOStruct *l, SRT_Boolean done, double *pv_cols,
                            long num_cols, char comp_origin[32]);

Err srt_f_IOstructgetcolpvs(SrtIOStruct l, double **pvs, long *num_pvs);

Err srt_f_IOstructsetexfrontier(SrtIOStruct *l, SRT_Boolean done,
                                double *exboundary, long exfrontierlength,
                                char comp_origin[32]);

Err srt_f_IOstructgetexfrontier(SrtIOStruct l, double **exboundary,
                                long *exfrontierlength);

/* -------------------------------------------------------------------------- */
/* ------------------------ Rate shift functions ---------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetrshiftval(SrtIOStruct l, char id_name[32],
                               double shift_value, int shift_type,
                               double *value);

/* -> return as '*value the value corresponding to a rate shift
        of value 'r_shift and of type shift_type  */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetrshift(SrtIOStruct l, char id_name[32], double shift_value,
                            int shift_type, SrtIOVal **ioval_pp);

/* -> return as '*ioval the SrtIOVal corresponding to a rate shift
        of value 'r_shift and of type shift_type  */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetrshift(SrtIOStruct *l, char id_name[32],
                            double shift_value, int shift_type,
                            SRT_Boolean done, double value,
                            char comp_origin[32], void *ptr);

/* -> put in the list '*l a cell with the corresponding parameters */

/* -------------------------------------------------------------------------- */
/* ---------------------- Sigma shift functions ----------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetsigmashiftval(SrtIOStruct l, int bucket_no,
                                   double shift_value, int shift_type,
                                   double *value);

/* -> return as '*ioval the SrtIOVal corresponding to a tau shift
        of value 'sigma_shift and of type shift_type  */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetsigmashift(SrtIOStruct l, int bucket_no,
                                double shift_value, int shift_type,
                                SrtIOVal **ioval_pp);

/* -> return as '*value from the list 'l the value corresponding to a
        sigma shift of type 'shift_type and shift value 'sigma_shift  ,
        for the bucket of value 'bucket_value */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetsigmashift(SrtIOStruct *l, int bucket_no,
                                Date bucket_start, Date bucket_end,
                                double shift_value, int shift_type,
                                SRT_Boolean done, double value,
                                char comp_origin[32]);

/* -> put in the list '*l a cell corresponding to the prameters */

/* -------------------------------------------------------------------------- */
/* ----------------------- Tau shift functions ------------------------------ */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgettaushiftval(SrtIOStruct l, int bucket_no,
                                 double shift_value, int shift_type,
                                 double *value);

/* -> return as '*value from the list 'l the value corresponding to a tau shift
        of type 'shift_type and shift value 'tau_shift  ,
        for the bucket of value 'bucket_value */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgettaushift(SrtIOStruct l, int bucket_no, double shift_value,
                              int shift_type, SrtIOVal **ioval_pp);

/* -> return as '*ioval the SrtIOVal corresponding to a tau shift
        of value 'tau_shift and of type shift_type  */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsettaushift(SrtIOStruct *l, int bucket_no, Date bucket_start,
                              Date bucket_end, double shift_value,
                              int shift_type, SRT_Boolean done, double value,
                              char comp_origin[32]);

/* -> put in the list '*l a cell corresponding to the parameters */

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructcreate(SrtIOStruct **l, char *name);

/* -> create a list (memory allocation and so on) with 1 element
        which represents the price */

#endif
