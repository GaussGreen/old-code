#ifndef NAGH02
#define NAGH02
#ifdef __cplusplus
extern "C"
{
#endif

/* <nagh02.h>
 *
 * Copyright 1997 Numerical Algorithms Group
 *
 * Include file for NAG C Library h02 Chapter
 *
 * Mark 5, 1997.
 */

#include <stdio.h>       /* For FILE type (also NULL)                */
#include <nag_h02mesg.h> /* Message codes and extern of message list */

#define ALL_NODES -999

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL h02bbc(Integer n, Integer m, double a[], Integer tda,
		   double bl[], double bu[], Boolean intvar[], double cvec[], 
		   double h[], Integer tdh,  
		   void (*qphess)(Integer n, Integer jthcol, double h[],
				  Integer tdh, double x[], double hx[], 
				  Nag_Comm *comm),
		   double x[], double *obj, Nag_H02_Opt *options, 
		   Nag_Comm *comm, NagError *fail);

extern void h02bbx(Boolean root, Integer num_nodes, Nag_EndState endstate, 
		   double obj, Integer n, Integer m, Integer istate[], 
		   double ax[], double bl[], double bu[], double lamda[], 
		   double x[], char **names, Nag_H02_Opt *opt, 
		   Nag_FileSt *stream, Nag_Comm *comm);

extern void h02bby(Integer nsolved, Integer nparent, Integer depth, 
		   Integer ix, double obj, double bl[], double bu[], double x[], 
		   double x_branch, Nag_NodeStatus status, Nag_H02_Opt *opt, 
		   Nag_FileSt *stream, Nag_Comm *comm);

extern Nag_BB_Fail h02bbz(Nag_MIP_Problem *problem, Nag_MIP_Solution *best, 
			  double bl_best[], double bu_best[],
			  Nag_H02_Opt *options, Nag_Comm *comm, Nag_FileSt *stream, 
			  NagError *ovflow, Nag_EndState *endstate);

extern NAG_DLL_EXPIMP void NAG_CALL h02buc(char *mps_file, Boolean minimising, Integer *n, Integer *m,
		   double **a, double **bl, double **bu, 
		   Boolean **intvar, double **cvec, double **x, 
		   Nag_H02_Opt *options, NagError *fail);

extern NAG_DLL_EXPIMP void NAG_CALL h02bvc(double **a, double **bl, double **bu, 
		   Boolean **intvar, double **cvec, double **x);

extern Boolean h02bux(Boolean allocate, Boolean alloc_hess, Nag_MIP_Problem *prob, 
		      int *alloc_status);

extern void h02buy(int alloc_status, Nag_MIP_Problem *prob, 
		   char *caller, NagError *fail);

extern void h02buz(FILE *fp_mps, Nag_MPS_Opt *opt, char *caller, 
		   Nag_MIP_Problem *prob,
		   Nag_FileSt *stream, NagError *fail);

extern Boolean h02xxy(char str[], int nc, int *field_code);

extern Boolean h02xxx(char str[], int nc, int *field_code);

extern NAG_DLL_EXPIMP void NAG_CALL h02xxc(Nag_H02_Opt *opt);

extern void h02xya(int field_code, Nag_H02_Opt *options, char buf[]);

extern NAG_DLL_EXPIMP void NAG_CALL h02xyc(const char *name, const char *opt_file, Nag_H02_Opt *opt,
		   Boolean print, const char *outfile, NagError *fail);

extern void h02xyx(int field_code, Nag_H02_Opt *options,
		   Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);

extern void h02xyy(int field_code, Nag_H02_Opt *options,
		   Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);

extern void h02xyz(Boolean (*valid_field)(char *str, int nc, int *field_code),
		   int *state, int *i, char line[], Integer *linenum,
		   int *field_code, FILE *fp, Nag_H02_Opt *opt, char str[],
		   Nag_Opt_Found *found);

extern NAG_DLL_EXPIMP void NAG_CALL h02xzc(Nag_H02_Opt *opt, char *name, NagError *fail);

extern void h02zxx(Nag_MIP_Problem *prob);
extern void h02zxy(Nag_MIP_Problem *prob, char *name, NagError *fail);
#else
extern void h02bbc();
extern void h02bbx();
extern void h02bby();
extern Nag_BB_Fail h02bbz();
extern void h02buc();
extern void h02bvc();
extern Boolean h02bux();
extern void h02buy();
extern void h02buz();
extern void h02xxc();
extern Boolean h02xxy();
extern Boolean h02xxx();
extern void h02xya();
extern void h02xyc();
extern void h02xyx();
extern void h02xyy();
extern void h02xyz();
extern void h02xzc();
extern void h02yxc();
extern void h02zxx();
extern void h02zzy();
#endif

#ifdef __cplusplus
}
#endif
#endif /* ndef NAGH02 */
