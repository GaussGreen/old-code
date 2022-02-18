#ifndef GRF_H_FLTFL
#define GRF_H_FLTFL

#include <stdio.h>

#include "grf_h_fltfldefs.h"
#include "grf_h_pubtypes.h"
#include "utdates.h"
#include "uterror.h"
#include "uttypes.h"

/*==========================================================*/

Err grf_f_readauxrng(FILE* in, long* auxwidthptr, long** auxlenptr, double*** auxptr);

Err grf_f_readcells(FILE* in, Date** eventdates, GrfnCell*** m, long* nr, long* nc);

Err grf_f_readgrfrng(FILE* in, long* numgrng, GrfnRng** grng);

/*==========================================================*/

Err grf_f_pringrfrng(FILE* out, long numgrng, GrfnRng* grng);

Err grf_f_prinauxrng(FILE* out, long auxwidth, long* auxlen, double** aux);

Err grf_f_princells(FILE* out, Date* eventdates, GrfnCell** m, long nr, long nc);

/*==========================================================*/
/*==========================================================*/

/*** grf_f_dbg_fct.c ***/

void print_sprdsht(long numeventdates, Date* eventdates, GrfnCell** sprdsht, FILE* f);

/*==========================================================*/
/*==========================================================*/

#endif
