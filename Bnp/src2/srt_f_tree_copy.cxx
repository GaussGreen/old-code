/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_TREE_COPY.C                                     */
/*                                                                            */
/*      PURPOSE:        Functions needed to compute greeks                    */
/*                                                                            */
/*      AUTHORS:        Arnaud Sahuguet                                       */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           February 1995                                         */
/*                                                                            */
/*      REASON:         add greek functions                                   */
/*                                                                            */
/******************************************************************************/

#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */

Err srt_f_tminfdup(SrtIRMTmInf *tminf, SrtIRMTmInf **tminf_dup) {
  Err err = 0;
  if (tminf == NULL) {
    (*tminf_dup) = NULL;
    return err;
  }
  (*tminf_dup) = (SrtIRMTmInf *)srt_malloc(sizeof(SrtIRMTmInf));
  (*tminf_dup)->df = tminf->df;
  (*tminf_dup)->ev.onef.sig = tminf->ev.onef.sig;
  (*tminf_dup)->ev.onef.sig2 = tminf->ev.onef.sig2;
  (*tminf_dup)->ev.onef.lambda = tminf->ev.onef.lambda;
  (*tminf_dup)->ev.onef.tau = tminf->ev.onef.tau;
  (*tminf_dup)->ev.onef.F = tminf->ev.onef.F;
  (*tminf_dup)->rf.onef.G = tminf->rf.onef.G;
  (*tminf_dup)->rf.onef.H = tminf->rf.onef.H;
  (*tminf_dup)->ev.onef.Psi = tminf->ev.onef.Psi;
  (*tminf_dup)->rf.onef.stdev_x = tminf->rf.onef.stdev_x;
  (*tminf_dup)->fwd_sam = tminf->fwd_sam;
  (*tminf_dup)->yp = tminf->yp;
  return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_trinfdup(SrtCheTreInf *trinf, SrtCheTreInf **trinf_dup) {
  Err err = 0;
  if (trinf == NULL) {
    (*trinf_dup) = NULL;
    return err;
  }
  (*trinf_dup) = (SrtCheTreInf *)srt_malloc(sizeof(SrtCheTreInf));
  (*trinf_dup)->max_r_index = trinf->max_r_index;
  (*trinf_dup)->min_r_index = trinf->min_r_index;
  (*trinf_dup)->max_phi_index = trinf->max_phi_index;
  (*trinf_dup)->min_phi_index = trinf->min_phi_index;
  (*trinf_dup)->logrmin = trinf->logrmin;
  (*trinf_dup)->rmin = trinf->rmin;
  (*trinf_dup)->u = trinf->u;
  return err;
}

/* -------------------------------------------------------------------------- */

SrtStpPtr srt_f_copystp(SrtStpPtr input_dest, SrtStpPtr input_source) {
  SrtStpPtr new_stp;
  if (input_source == NULL)
    return input_dest;
  else {
    new_stp = (SrtStpPtr)srt_malloc(sizeof(SrtStp));
    new_stp->index = input_source->index;
    new_stp->date = input_source->date;
    new_stp->time = input_source->time;
    new_stp->ddate = input_source->ddate;
    new_stp->delta_t = input_source->delta_t;
    new_stp->sqrt_delta_t = input_source->sqrt_delta_t;
    new_stp->tminf = srt_calloc(1, sizeof(void *));
    srt_f_tminfdup((SrtIRMTmInf *)(input_source->tminf[0]),
                   (SrtIRMTmInf **)&(new_stp->tminf[0]));
    srt_f_trinfdup((SrtCheTreInf *)(input_source->trinf),
                   (SrtCheTreInf **)&(new_stp->trinf));
    new_stp->e = input_source->e;
    new_stp->prev = input_dest;
    new_stp->next = NULL;

    if (input_dest != NULL)
      input_dest->next = new_stp;

    input_dest = new_stp;

    return srt_f_copystp(input_dest, input_source->next);
  }
}

/* -------------------------------------------------------------------------- */

Err srt_f_stpdup(SrtStpPtr stp, SrtStpPtr *stp_dup) {
  Err err = 0;
  SrtStpPtr top;
  top = gototop(stp);
  (*stp_dup) = srt_f_copystp(NULL, top);
  return err;
}

/* -------------------------------------------------------------------------- */
/*
void my_print_Step(SrtStpPtr s        , char file_name[]        , int index)
{
        SrtCheTreInf *trinf;
        SrtCheTmInf *tminf;
        FILE* f;
        char name[32]="";
        char index_n[2]="";

        strncat(name        ,file_name        ,4);
        name[4] = '\0';
        strncat(name        ,".dta"        ,strlen(".dta"));
        sprintf(index_n        ,"%d"        ,index);
        strncat(name        ,index_n        ,strlen(index_n));

        f = fopen(name        ,"w");
        s = gototop(s);
        fprintf(f        ,"steps:\n");
        while(s!=NULL)
        {
                tminf=s->tminf;
                trinf=s->trinf;

fprintf(f        ,"index %2d\tdate %.4lf\t time %.4lf\t delta_t%.4lf\n"        ,
                s->index        ,
                s->ddate        ,
                s->time        ,
                s->delta_t);

if ( tminf != NULL ) fprintf(f        ,"\tdf %.4lf\n\tsig %.4lf\n\tlambda
%.4lf\n\ttau
%.4lf\n\tF %.4lf\n\tG %.4lf\n\tH %.4lf\n\tPsi %.4lf\n\tStd %.4lf\n\trate
%.8lf\n\tPhi %.4lf\n\typ %.4lf\n\t\n"        , tminf->df        , tminf->sig ,
                                tminf->lambda        ,
                                tminf->tau        ,
                                tminf->F        ,
                                tminf->G        ,
                                tminf->H        ,
                                tminf->Psi        ,
                                tminf->stdev_r        ,
                                100* tminf->fwd_sam.short_rate        ,
                                tminf->fwd_sam.phi        ,
                                tminf->yp	);
else fprintf(f        ,"\t tminf = NULL \n");

if (trinf != NULL ) fprintf(f        ,"\t\trmin %.8lf\n\t\tu
%.4lf\n\t\tmax_r_index %d\n\t\tmin_r_index %d\n"        , trinf->rmin        ,
trinf->u        , trinf->max_r_index        , trinf->min_r_index	); else
fprintf(f        ,"\t\t trinf = NULL \n");

                s = s->next;
        }
        fclose(f);
}

*/
