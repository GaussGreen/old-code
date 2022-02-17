/*=====================================================================*/
/*  lstypes.h:   Contains:                                             */
/*   1. typedefs for function pointers: gcomp(function evaluation)     */
/*                                      parsh(user derivatives)        */
/*            passed by user to grgsub                                 */
/*                                                                     */
/*                                                                     */
/*  NOTE:  gcomp pointer is required, parsh pointer is optional        */
/*         and is not checked or used unless user makes call to        */
/*         lsgrg_setparameter() to set the derivative computation mode */
/*         to analytic deriatives. In general, it is a good thing      */
/*         set parsh pointer to NULL before call to grgsub             */
/*                                                                     */
/*   ----  Protptypes for the functions are supplied below -------     */
/*   ----  Please see user documentation for argument specifics --     */
/*                                                                     */
/*  P_GCOMP is type for pointer to function which evaluates the        */
/*          problem functions                                          */
/*    void gcomp(double g[], double x[]);                              */
/*  where x is the array of variable values, g the array of function   */
/*  values (computed by gcomp)                                         */
/*                                                                     */
/*  P_PARSH is type for pointer to function which returns the          */
/*     derivatives of the problem functions                            */
/* void lsgrg_user_parsh(double x[],long n,double paij[],long iprow[], */
/*                     long ipcol[], long *nnonZeros);                 */
/*                                                                     */
/*                                                                     */
/*   P_GCOMPX is pointer to gcomp function which allows individual     */
/*    function access.  This will not be implemented in the initial    */
/*    version of lsgrgc3.                                              */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*=====================================================================*/

/*
#ifdef USER_TERMINATION_ENABLED
typedef long (*P_GCOMP)  (double *, double *);
typedef long (*P_GCOMPX) (double *, double *, long*, long[]);
#else
typedef void (*P_GCOMP)  (double *, double *);
typedef void (*P_GCOMPX) (double *, double *, long*, long[]);
#endif
typedef void (*P_PARSH)  (double[],long, double[],long[],long[], long*);
*/

#ifdef USER_TERMINATION_ENABLED
	typedef long (LSGRGSolver::*P_GCOMP)  (double *, double *);
	typedef long (LSGRGSolver::*P_GCOMPX) (double *, double *, long*, long[]);
#else
	typedef void (LSGRGSolver::*P_GCOMP)  (double *, double *);
	typedef void (LSGRGSolver::*P_GCOMPX) (double *, double *, long*, long[]);
#endif
	typedef void (LSGRGSolver::*P_PARSH)  (double[],long, double[],long[],long[], long*);

