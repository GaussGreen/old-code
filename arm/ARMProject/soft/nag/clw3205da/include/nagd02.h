#ifndef NAGD02
#define NAGD02
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagd02.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library d02 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2145 (Feb 1998).
 */

/* Note :
 * d02qfg is called only by d02qfv.
 * A subsection of the original d02qfv was placed in
 * this external function to allow compilation on PC's.
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02cjc(Integer neq, 
            NAG_D02CJC_FUN fcn,
            double NAG_HUGE *t, double NAG_HUGE y[], double tend,
            double tol, Nag_ErrorControl err_c,
            NAG_D02CJC_OUTFUN output,
            NAG_D02CJC_GFUN g,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d02cjc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02ejc(Integer neq, NAG_D02EJC_FUN fcn,
	    NAG_D02EJC_PFUN pederv,
            double NAG_HUGE *t, double NAG_HUGE y[], double tend, 
            double tol, Nag_ErrorControl err_c,
	    NAG_D02EJC_OUTFUN output, NAG_D02EJC_GFUN g,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d02ejc();
#endif

#ifdef NAG_PROTO
extern void d02eju(double NAG_HUGE *xout, NAG_D02EJC_OUTFUN output,
            double NAG_HUGE chk[], Integer n, 
            double NAG_HUGE *ysav, double NAG_HUGE acor[], double NAG_HUGE *x, Integer NAG_HUGE *nqu, double NAG_HUGE *hu,
            double NAG_HUGE *h, Integer NAG_HUGE *iflag, double NAG_HUGE *xlast, double NAG_HUGE *xnew, 
            double dir, Nag_User NAG_HUGE *comm, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02eju();
#endif

#ifdef NAG_PROTO
extern void d02ejv(Integer n,double NAG_HUGE *x, double NAG_HUGE y[], double NAG_HUGE *hu, double NAG_HUGE *xlast, 
            double NAG_HUGE *h, double NAG_HUGE *rnqu, double NAG_HUGE ysav[], double NAG_HUGE acor[],
            double NAG_HUGE chk[], Integer NAG_HUGE *imon, Integer NAG_HUGE *nstps, double NAG_HUGE *glast,
	    NAG_D02EJC_GFUN g,
            Integer NAG_HUGE *iflag, double NAG_HUGE d[], Integer NAG_HUGE *ifin, double dir, 
            double NAG_HUGE *root, NAG_D02EJC_OUTFUN output, 
            Integer NAG_HUGE *path, double NAG_HUGE *xout, 
            Nag_User NAG_HUGE *comm, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02ejv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02gac(Integer neq, 
            NAG_D02GAC_FUN fcn,
            double a, double b, double NAG_HUGE u[], Integer NAG_HUGE v[],
            Integer mnp, Integer NAG_HUGE *np, double NAG_HUGE x[], double NAG_HUGE y[],
            double tol, Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d02gac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02gbc(Integer neq, 
            NAG_D02GBC_FUN fcnf,
            NAG_D02GBC_GFUN fcng,
            double a, double b, double NAG_HUGE c[], double NAG_HUGE d[], double NAG_HUGE gam[],
            Integer mnp, Integer NAG_HUGE *np, double NAG_HUGE x[], double NAG_HUGE y[], double tol,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d02gbc();
#endif

#ifdef NAG_PROTO
extern void d02gbz(double NAG_HUGE e[], double NAG_HUGE f[], Integer n, double NAG_HUGE c[],
            
            double NAG_HUGE d[], Integer NAG_HUGE * ind, double NAG_HUGE w[]);
#else
extern void d02gbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02nmc(Integer neq, Integer neqmax, double NAG_HUGE *t, double tout, double NAG_HUGE y[],
            double NAG_HUGE ydoti[], double NAG_HUGE rwork[], double NAG_HUGE rtol[], 
            double NAG_HUGE a_tol[], Integer itol, Integer NAG_HUGE inform[], double NAG_HUGE ysave[], 
            Integer ny2dim, double NAG_HUGE wkjac[], Integer nwkjac, Integer NAG_HUGE jacpvt[],
            Integer njcpvt, Integer NAG_HUGE *imon, Integer NAG_HUGE *inln, Integer NAG_HUGE *ires, 
            Integer NAG_HUGE *irevcm, Integer itask, Integer jtrace, Integer NAG_HUGE *ifail,
            Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nmc();
#endif

#ifdef NAG_PROTO
extern void d02nmq(double NAG_HUGE wm[], Integer NAG_HUGE iwm[], double NAG_HUGE x[], Integer n, Integer NAG_HUGE *ier,
            Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nmq();
#endif

#ifdef NAG_PROTO
extern void d02nmt(double NAG_HUGE a[], Integer lda, Integer n, double NAG_HUGE rdae[], 
            double h, double el0);
#else
extern void d02nmt();
#endif

#ifdef NAG_PROTO
extern void d02nmu(Integer neq, double NAG_HUGE y[], double NAG_HUGE *yh, Integer nyh, double NAG_HUGE ewt[], 
            double NAG_HUGE rtem[], double NAG_HUGE savr[], double NAG_HUGE *ydot, double NAG_HUGE wm[], 
            Integer NAG_HUGE *iwm, Integer NAG_HUGE *ifj, double h, double el0, double NAG_HUGE *tn, 
            Integer NAG_HUGE *ifunc, double NAG_HUGE rdae[], double NAG_HUGE rworkx[], Integer NAG_HUGE *irevcm,
            Integer NAG_HUGE inform[], Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nmu();
#endif

#ifdef NAG_PROTO
extern void d02nmv(Integer n, double NAG_HUGE y[], double NAG_HUGE ydoti[], double NAG_HUGE yh[], Integer nyh,
            double NAG_HUGE savr[], double NAG_HUGE acor[], double NAG_HUGE ewt[], Integer NAG_HUGE *ifunc, 
            Integer NAG_HUGE *inln, double NAG_HUGE *h, double NAG_HUGE *el0, double NAG_HUGE rdae[],
            Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nmv();
#endif

#ifdef NAG_PROTO
extern void d02nmw(Integer meth, double NAG_HUGE elco[], double NAG_HUGE tesco[]);
#else
extern void d02nmw();
#endif

#ifdef NAG_PROTO
extern void d02nmx(Integer neq, double NAG_HUGE y[], double NAG_HUGE *yh, Integer nyh, double NAG_HUGE ewt[],
            double NAG_HUGE ydot[], double NAG_HUGE savr[], double NAG_HUGE acor[], Integer NAG_HUGE *inln, 
            Integer NAG_HUGE *istep, double NAG_HUGE *el0, double NAG_HUGE *h, double NAG_HUGE *tn, double NAG_HUGE *hmin,
            double NAG_HUGE *hmxi, double NAG_HUGE rdae[], double NAG_HUGE rworkx[], Integer NAG_HUGE *irevcm,
            Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nmx();
#endif

#ifdef NAG_PROTO
extern void d02nnm(Integer inln);
#else
extern void d02nnm();
#endif

#ifdef NAG_PROTO
extern void d02nnn(double NAG_HUGE x[], Integer n, Integer imsg, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nnn();
#endif

#ifdef NAG_PROTO
extern void d02nnw(Integer neq, double NAG_HUGE *t, double tout, double NAG_HUGE *h0, double NAG_HUGE y[],
            double NAG_HUGE ydoti[], double NAG_HUGE ewt[], double NAG_HUGE rtol[], double NAG_HUGE a_tol[], 
            Integer itol, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nnw();
#endif

#ifdef NAG_PROTO
extern void d02nnx(Integer NAG_HUGE *n, Integer itol, double NAG_HUGE rtol[], double NAG_HUGE a_tol[], 
            double NAG_HUGE ycur[], double NAG_HUGE ewt[], Integer NAG_HUGE *iewset, 
            Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02nnx();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02nsc(Integer neq, Integer neqmax, char NAG_HUGE *jceval, Integer nwkjac, 
            double NAG_HUGE rwork[], Integer NAG_HUGE *ifail);
#else
extern void d02nsc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02nvc(Integer neqmax, Integer ny2dim, Integer maxord, char NAG_HUGE *method, 
            Boolean petzld, double NAG_HUGE const_[], double tcrit, double hmin, 
            double hmax, double h0, Integer maxstp, Integer mxhnil, 
            char NAG_HUGE *norm, double NAG_HUGE rwork[], Integer NAG_HUGE *ifail);
#else
extern void d02nvc();
#endif

#ifdef NAG_PROTO
extern void d02nvz(Integer neqmax, Integer ny2dim, Integer maxord, Integer NAG_HUGE *meth,
            Boolean petzld, double NAG_HUGE const_[], double tcrit, double hmin,
            double hmax, double h0, Integer maxstp, Integer mxhnil, 
            char NAG_HUGE *norm, double NAG_HUGE rwork[], Integer NAG_HUGE *ierr, Boolean NAG_HUGE *report);
#else
extern void d02nvz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pcc(Integer neq, 
            NAG_D02PCC_FUN f,
            double twant, double NAG_HUGE *tgot,
            double NAG_HUGE ygot[], double NAG_HUGE ypgot[], double NAG_HUGE ymax[],
            Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm, NagError NAG_HUGE * fail);
#else
extern void d02pcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pdc(Integer neq,
            NAG_D02PCC_FUN f,
            double NAG_HUGE *tnow, double NAG_HUGE ynow[],
            double NAG_HUGE ypnow[], Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm,
            NagError NAG_HUGE *fail);
#else
extern void d02pdc();
#endif

#ifdef NAG_PROTO
extern void d02pdm(Boolean ask, char NAG_HUGE *srname, Integer NAG_HUGE * state, Nag_ODE_RK NAG_HUGE *opt);
#else
extern void d02pdm();
#endif

#ifdef NAG_PROTO
extern void d02pdp(Integer NAG_HUGE *ier, char NAG_HUGE *srname, Integer NAG_HUGE *ifail, Nag_ODE_RK NAG_HUGE *opt);
#else
extern void d02pdp();
#endif

#ifdef NAG_PROTO
extern double d02pdq(double NAG_HUGE u[], double NAG_HUGE v[], double NAG_HUGE wt[], Integer neq);
#else
extern double d02pdq();
#endif

#ifdef NAG_PROTO
extern void d02pdr(double v1v1, double v0v1, double v0v0,
            double NAG_HUGE *rold, double NAG_HUGE *rho, double NAG_HUGE root1[],
            double NAG_HUGE root2[], Boolean NAG_HUGE * rootre);
#else
extern void d02pdr();
#endif

#ifdef NAG_PROTO
extern void d02pds(double alpha, double beta, double NAG_HUGE r1[],
            double NAG_HUGE r2[]);
#else
extern void d02pds();
#endif

#ifdef NAG_PROTO
extern void d02pdt(Integer neq, double NAG_HUGE v[], double havg, double x,
            double NAG_HUGE y[], NAG_D02PCC_FUN f,
            double NAG_HUGE fxy[], double NAG_HUGE wt[],
            double scale, double vdotv, double NAG_HUGE z[],
            double NAG_HUGE *zdotz, double NAG_HUGE vtemp[], Nag_ad02pd NAG_HUGE *ad02pd_1,
            Nag_bd02pd NAG_HUGE *bd02pd_1, Nag_User NAG_HUGE *comm);
#else
extern void d02pdt();
#endif

#ifdef NAG_PROTO
extern void d02pdu(Integer neq,
            NAG_D02PCC_FUN f,
            double x, double NAG_HUGE y[],
            double hnow, double havg, double xend,
            Integer maxfcn, double NAG_HUGE wt[], double NAG_HUGE fxy[],
            double NAG_HUGE v0[], Boolean NAG_HUGE * unsure, Boolean NAG_HUGE * stif, double NAG_HUGE v1[],
            double NAG_HUGE v2[], double NAG_HUGE v3[], double NAG_HUGE vtemp[], Nag_ODE_RK NAG_HUGE *opt,
            Nag_User NAG_HUGE *comm);
#else
extern void d02pdu();
#endif

#ifdef NAG_PROTO
extern void d02pdv(Integer neq, 
            NAG_D02PCC_FUN f,
            double havg, Integer NAG_HUGE * jflstp,
            Boolean toomch, Integer maxfcn, double NAG_HUGE work[],
            Integer NAG_HUGE * ier, Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm,
            char NAG_HUGE buf[]);
#else
extern void d02pdv();
#endif

#ifdef NAG_PROTO
extern void d02pdw(
            NAG_D02PCC_FUN f,
            Integer neq, double NAG_HUGE y[], double tol,
            double NAG_HUGE weight[], double NAG_HUGE zy[], double NAG_HUGE zyp[],
            double NAG_HUGE zerror[], double NAG_HUGE zynew[], double NAG_HUGE zerres[],
            double NAG_HUGE zstage[], Integer NAG_HUGE * ier, Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm);
#else
extern void d02pdw();
#endif

#ifdef NAG_PROTO
extern void d02pdx(Integer neq, double NAG_HUGE y[], double NAG_HUGE yp[], double h,
            double NAG_HUGE ynew[], double NAG_HUGE stages[], double NAG_HUGE thres[],
            double NAG_HUGE *err, Boolean main, double NAG_HUGE weight[], Nag_dd02pd NAG_HUGE *dd02pd_1);
#else
extern void d02pdx();
#endif

#ifdef NAG_PROTO
extern void d02pdy(double tnow, double NAG_HUGE y[], double NAG_HUGE yp[],
            double tstg, double NAG_HUGE ystg[], double NAG_HUGE ypstg[],
            double NAG_HUGE *htry, double NAG_HUGE weight[], Boolean NAG_HUGE * cutbak, 
            Nag_ODE_RK NAG_HUGE *opt);
#else
extern void d02pdy();
#endif

#ifdef NAG_PROTO
extern void d02pdz(
            NAG_D02PCC_FUN f,
            Integer neq, double tnow, double NAG_HUGE y[],
            double NAG_HUGE yp[], double NAG_HUGE stages[], double tol,
            double NAG_HUGE *htry, double NAG_HUGE weight[], double NAG_HUGE ynew[],
            double NAG_HUGE errest[], double NAG_HUGE *err, Boolean main,
            double hmin, double NAG_HUGE thres[], Boolean NAG_HUGE * phase2,
            Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm);
#else
extern void d02pdz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02ppc(Nag_ODE_RK NAG_HUGE *opt);
#else
extern void d02ppc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pvc(Integer neq, double tstart, double NAG_HUGE ystart[],
            double tend, double tol, double NAG_HUGE thres[],
            Nag_RK_method method, Nag_RK_task task, Nag_ErrorAssess errass,
            double hstart, Nag_ODE_RK NAG_HUGE *opt, NagError NAG_HUGE * fail);
#else
extern void d02pvc();
#endif

#ifdef NAG_PROTO
extern void d02pvx(double NAG_HUGE *mcheps, double NAG_HUGE *dwarf);
#else
extern void d02pvx();
#endif

#ifdef NAG_PROTO
extern void d02pvy(Nag_gd02pd NAG_HUGE *gd02pd_1);
#else
extern void d02pvy();
#endif

#ifdef NAG_PROTO
extern void d02pvz(Integer method, Integer NAG_HUGE * vecstg, Boolean NAG_HUGE * reqstg,
            Integer NAG_HUGE * lintpl, Nag_dd02pd NAG_HUGE *dd02pd_1, Nag_ed02pd NAG_HUGE *ed02pd_1,
            Nag_gd02pd NAG_HUGE *gd02pd_1);
#else
extern void d02pvz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pwc(double tend_new, Nag_ODE_RK NAG_HUGE *opt, NagError NAG_HUGE *fail);
#else
extern void d02pwc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pxc(Integer neq, double twant, Nag_SolDeriv reqest, Integer nwant,
            double NAG_HUGE ywant[], double NAG_HUGE ypwant[],
            NAG_D02PDC_FUN f,
            Nag_ODE_RK NAG_HUGE *opt, Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d02pxc();
#endif

#ifdef NAG_PROTO
extern void d02pxy(
            NAG_D02PDC_FUN f,
            Integer neq, Integer nwant, double NAG_HUGE y[],
            double NAG_HUGE yp[], double NAG_HUGE yold[], double NAG_HUGE ypold[],
            double NAG_HUGE stages[], Boolean calstg, double NAG_HUGE xstage[],
            double NAG_HUGE ytemp[],
            double NAG_HUGE p[], Nag_ad02pd NAG_HUGE *ad02pd_1, Nag_bd02pd NAG_HUGE *bd02pd_1,
            Nag_dd02pd NAG_HUGE *dd02pd_1, Nag_User NAG_HUGE *comm);
#else
extern void d02pxy();
#endif

#ifdef NAG_PROTO
extern void d02pxz(double NAG_HUGE y[], double NAG_HUGE yp[], double NAG_HUGE p[],
            double twant, char reqest, Integer nwant,
            double NAG_HUGE ywant[], double NAG_HUGE ypwant[],
            Nag_bd02pd NAG_HUGE *bd02pd_1, Nag_dd02pd NAG_HUGE *dd02pd_1);
#else
extern void d02pxz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pyc(Integer NAG_HUGE *totfcn, Integer NAG_HUGE * stpcst, double NAG_HUGE *waste,
            Integer NAG_HUGE *stpsok, double NAG_HUGE *hnext, Nag_ODE_RK NAG_HUGE *opt);
#else
extern void d02pyc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02pzc(Integer neq, double NAG_HUGE rmserr[], double NAG_HUGE *errmax, double NAG_HUGE *terrmx,
            Nag_ODE_RK NAG_HUGE *opt, NagError NAG_HUGE *fail);
#else
extern void d02pzc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02qfc(Integer neqf, 
            NAG_D02QFC_FUN fcn,
            double NAG_HUGE *t, double NAG_HUGE y[], double tout,
            NAG_D02QFC_GFUN g,
            Nag_User NAG_HUGE *comm, Nag_ODE_Adams NAG_HUGE *opt, NagError NAG_HUGE *fail);
#else
extern void d02qfc();
#endif

#ifdef NAG_PROTO
extern void d02qfg(Boolean NAG_HUGE *ret,
            int NAG_HUGE *irevcm, double NAG_HUGE *twant, Integer NAG_HUGE *kwant, double gwant,
            Integer neq, double NAG_HUGE *t, double NAG_HUGE y[], double NAG_HUGE ypout[],
            double x, double NAG_HUGE yy[], double NAG_HUGE r[], double NAG_HUGE r2d[], 
            void (NAG_HUGE *rintrp)(double t, double NAG_HUGE yint[], double NAG_HUGE ypint[],
                           Integer neq, double x, double NAG_HUGE yy[], double NAG_HUGE p[],
                           double NAG_HUGE phi[], Nag_ad02 NAG_HUGE *ad02qf),
            Integer neqg, double NAG_HUGE gold[],
            double NAG_HUGE gnew[], double NAG_HUGE tgv[], double NAG_HUGE gv[], double NAG_HUGE tkt[],
            double NAG_HUGE tlbmr[], double NAG_HUGE trbmr[], double NAG_HUGE proot[],
            double NAG_HUGE rootd[], Integer NAG_HUGE mmreq[], Integer NAG_HUGE indxg[],
            Integer NAG_HUGE igsc[], Integer NAG_HUGE needgk[],
            Integer NAG_HUGE *kroot, Integer NAG_HUGE *inroot,
            Nag_ad02 NAG_HUGE *ad02qf, Nag_bd02 NAG_HUGE *bd02qf,
            Nag_dd02 NAG_HUGE *dd02qf, Nag_ed02 NAG_HUGE *ed02qf,
            Nag_fd02 NAG_HUGE *fd02qf, Nag_gd02 NAG_HUGE *gd02qf);
#else
extern void d02qfg();
#endif

#ifdef NAG_PROTO
extern double d02qfn(double NAG_HUGE v[], Integer ncomp);
#else
extern double d02qfn();
#endif

#ifdef NAG_PROTO
extern void d02qfp(int NAG_HUGE *irevcm, Integer neq, double NAG_HUGE *twant,
            double b, double NAG_HUGE y[], double NAG_HUGE yprime[],
            double NAG_HUGE etol[], Integer morder, double small1,
            double big, double NAG_HUGE spy[], double NAG_HUGE pv[], double NAG_HUGE yp[],
            double NAG_HUGE sf[], double NAG_HUGE *h, Nag_ODE_Adams NAG_HUGE *opt);
#else
extern void d02qfp();
#endif

#ifdef NAG_PROTO
extern void d02qfq(int NAG_HUGE *irevcm, double NAG_HUGE *twant, Integer neqn,
            double NAG_HUGE y[], double NAG_HUGE *x, double NAG_HUGE *h, double NAG_HUGE *eps,
            double NAG_HUGE wt[], Boolean NAG_HUGE *start, double NAG_HUGE *hold,
            Integer NAG_HUGE *k, Integer NAG_HUGE *kold, Boolean NAG_HUGE *crash, double NAG_HUGE phi[],
            double NAG_HUGE p[], double NAG_HUGE yp[], double NAG_HUGE psi[], double NAG_HUGE alpha[],
            double NAG_HUGE beta[], double NAG_HUGE sig[], double NAG_HUGE v[], double NAG_HUGE w[],
            double NAG_HUGE g[], Boolean NAG_HUGE *phase1, Integer NAG_HUGE *ns, Boolean NAG_HUGE *nornd,
            Integer NAG_HUGE *ksteps, double twou, double fouru,
            double NAG_HUGE *xold, Integer NAG_HUGE *kprev, Integer NAG_HUGE *ivc, Integer NAG_HUGE iv[],
            Integer NAG_HUGE *kgi, double NAG_HUGE gi[], Integer NAG_HUGE *nsucc, Integer NAG_HUGE *nfail,
            Nag_ODE_Adams NAG_HUGE *opt);
#else
extern void d02qfq();
#endif

#ifdef NAG_PROTO
extern void d02qfr(double x, double NAG_HUGE y[], double xout,
            double NAG_HUGE yout[], double NAG_HUGE ypout[], Integer nint,
            Integer neqn, Integer kold, double NAG_HUGE phi[], Integer ivc,
            Integer NAG_HUGE iv[], Integer kgi, double NAG_HUGE gi[], double NAG_HUGE alpha[],
            double NAG_HUGE og[], double NAG_HUGE ow[], double ox,
            double NAG_HUGE oy[]);
#else
extern void d02qfr();
#endif

#ifdef NAG_PROTO
extern void d02qfs(double t, double NAG_HUGE yint[], double NAG_HUGE ypint[],
            Integer neq, double x, double NAG_HUGE yy[], double NAG_HUGE p[],
            double NAG_HUGE phi[], Nag_ad02 NAG_HUGE *ad02qf_1);
#else
extern void d02qfs();
#endif

#ifdef NAG_PROTO
extern void d02qft(double NAG_HUGE *a, double NAG_HUGE *fa, double NAG_HUGE *b, double NAG_HUGE *fb,
            double rez, double aez, Integer NAG_HUGE *iflag, Nag_dd02 NAG_HUGE *dd02qf_1);
#else
extern void d02qft();
#endif

#ifdef NAG_PROTO
extern void d02qfu(int NAG_HUGE *irevcm, Integer NAG_HUGE *k, double gwant,
            Integer neqg, Integer NAG_HUGE *kroot, Integer NAG_HUGE *inroot,
            double NAG_HUGE gold[], double NAG_HUGE proot[], double NAG_HUGE rootd[],
            double NAG_HUGE gp[], Integer NAG_HUGE igsc[], double tout,
            Nag_ODE_Adams NAG_HUGE *opt, Nag_bd02 NAG_HUGE *bd02qf_1);
#else
extern void d02qfu();
#endif

#ifdef NAG_PROTO
extern void d02qfv(int NAG_HUGE *irevcm, double NAG_HUGE *twant, Integer NAG_HUGE *kwant, double gwant,
            Integer neq, double NAG_HUGE *t, double NAG_HUGE y[], double NAG_HUGE ypout[],
            double tout, double x, double NAG_HUGE yy[], double NAG_HUGE r[], double NAG_HUGE r2d[],
            void (NAG_HUGE *rintrp)(double t, double NAG_HUGE yint[], double NAG_HUGE ypint[],
                           Integer neq, double x, double NAG_HUGE yy[], double NAG_HUGE p[],
                           double NAG_HUGE phi[], Nag_ad02 NAG_HUGE *ad02qf),
            Integer neqg, double NAG_HUGE gold[], double NAG_HUGE gnew[],
            double NAG_HUGE gp[], double NAG_HUGE tgv[], double NAG_HUGE gv[], double NAG_HUGE tkt[],
            double NAG_HUGE tlbmr[], double NAG_HUGE trbmr[], double NAG_HUGE proot[],
            double NAG_HUGE rootd[], Integer NAG_HUGE mmreq[], Integer NAG_HUGE indxg[],
            Integer NAG_HUGE igsc[], Integer NAG_HUGE needgk[], Integer NAG_HUGE *kroot,
            Integer NAG_HUGE *inroot, Integer NAG_HUGE *izflag, Nag_ODE_Adams NAG_HUGE *opt,
            Nag_ad02 NAG_HUGE *ad02qf, Nag_bd02 NAG_HUGE *bd02qf);
#else
extern void d02qfv();
#endif

#ifdef NAG_PROTO
extern void d02qfw(Integer neqg, Integer NAG_HUGE *kroot, Integer NAG_HUGE *inroot,
            double NAG_HUGE tkt[], double NAG_HUGE gold[], double NAG_HUGE proot[],
            double NAG_HUGE rootd[], double NAG_HUGE gp[], Integer NAG_HUGE needgk[],
            Integer NAG_HUGE igsc[], double t, Nag_bd02 NAG_HUGE *bd02qf_1);
#else
extern void d02qfw();
#endif

#ifdef NAG_PROTO
extern void d02qfx(int NAG_HUGE *irevcm, double NAG_HUGE *twant, Integer NAG_HUGE *kwant,
            double gwant, Integer neq, double NAG_HUGE *t,
            double NAG_HUGE y[], double tout, double NAG_HUGE rtol[],
            double NAG_HUGE a_tol[], Integer NAG_HUGE *idid, double NAG_HUGE ypout[],
            double NAG_HUGE yp[], double NAG_HUGE yy[], double NAG_HUGE wt[], double NAG_HUGE p[],
            double NAG_HUGE phi[], double NAG_HUGE gold[], double NAG_HUGE gnew[],
            double NAG_HUGE tgv[], double NAG_HUGE gv[], double NAG_HUGE gp[], double NAG_HUGE tkt[],
            double NAG_HUGE tlbmr[], double NAG_HUGE trbmr[], double NAG_HUGE proot[],
            double NAG_HUGE rootd[], double tstop, double NAG_HUGE *h,
            double NAG_HUGE *eps, double NAG_HUGE *x, double hmax,
            Integer maxnum, Integer NAG_HUGE *nsucc, Integer NAG_HUGE *nfail, Integer NAG_HUGE indxg[],
            Integer NAG_HUGE igsc[], Integer NAG_HUGE mmreq[], Integer NAG_HUGE needgk[], Integer NAG_HUGE *kroot,
            Integer NAG_HUGE *inroot, Integer neqg, Integer NAG_HUGE *badcmp, Nag_ODE_Adams NAG_HUGE *opt,
            Nag_ad02 NAG_HUGE *ad02qf_1, Nag_bd02 NAG_HUGE *bd02qf_1);
#else
extern void d02qfx();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02qwc(Nag_Start NAG_HUGE *state, Integer neqf, Boolean vectol, double NAG_HUGE a_tol[],
            double NAG_HUGE rtol[], Boolean onestp, Boolean crit, double tcrit,
            double hmax, Integer maxstp, Integer neqg, Boolean NAG_HUGE *alterg,
            Boolean sophst, Nag_ODE_Adams NAG_HUGE *opt, NagError NAG_HUGE *fail);
#else
extern void d02qwc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02qyc(Nag_ODE_Adams NAG_HUGE *opt);
#else
extern void d02qyc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02qzc(Integer neqf, double twant, Integer nwant,
            double NAG_HUGE ywant[], double NAG_HUGE ypwant[], Nag_ODE_Adams NAG_HUGE *opt,
            NagError NAG_HUGE *fail);
#else
extern void d02qzc();
#endif

#ifdef NAG_PROTO
extern void d02qzz(double t, double NAG_HUGE yint[], double NAG_HUGE ypint[],
            Integer nint, Integer neq, double x, double NAG_HUGE yy[],
            double NAG_HUGE p[], double NAG_HUGE phi[], Nag_ad02 NAG_HUGE *ad02qf_1);
#else
extern void d02qzz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02rac(Integer neq, double NAG_HUGE *deleps,
            NAG_D02RAC_FUN fcn,
            Integer numbeg, Integer nummix,
            NAG_D02RAC_GFUN g,
            Nag_MeshSet init, Integer mnp, Integer NAG_HUGE * np, double NAG_HUGE x[],
            double NAG_HUGE y[], 
            double tol, double NAG_HUGE abt[],
            NAG_D02RAC_JFUN jacobf,
            NAG_D02RAC_JGFUN jacobg,            
            NAG_D02RAC_JEPSFUN jaceps,
            NAG_D02RAC_JGEPSFUN jacgep,
            Nag_User NAG_HUGE *comm,
            NagError NAG_HUGE * fail);
#else
extern void d02rac();
#endif

#ifdef NAG_PROTO
extern void d02rar (Integer m, Integer n, 
		    Integer NAG_HUGE irn[], Integer nirn, Integer NAG_HUGE ip[], Integer nip, 
		    double NAG_HUGE h[], double NAG_HUGE x[], double NAG_HUGE f[],
		    double NAG_HUGE hmax[], double umin, double uaim,
		    double umax, double NAG_HUGE *eps, double NAG_HUGE *eps1,
		    Boolean NAG_HUGE * adjust,
		    double NAG_HUGE a[], 
		    Integer NAG_HUGE ig[], Integer nig,
		    double NAG_HUGE w[], 
		    Integer nw, 
		    double NAG_HUGE z[], 
		    Integer NAG_HUGE * ind,
		    Integer NAG_HUGE * ifail,
		    Nag_ad02ra *ad02ra);
#else
extern void d02rar();
#endif

#ifdef NAG_PROTO
extern void d02ras(double NAG_HUGE a[], double NAG_HUGE c[], double NAG_HUGE del[], double NAG_HUGE y[],
            Integer m, Integer NAG_HUGE * n, Integer p, Integer r, Integer NAG_HUGE ir[],
            Integer NAG_HUGE ic[], double NAG_HUGE u[], Integer NAG_HUGE * mtnmax, Integer nmax,
            double NAG_HUGE x[], Boolean NAG_HUGE * sing);
#else
extern void d02ras();
#endif

#ifdef NAG_PROTO
extern void d02rat(double NAG_HUGE a[], double NAG_HUGE c[], double NAG_HUGE del[], Integer m,
            Integer NAG_HUGE * n, Integer p, Integer r, Integer NAG_HUGE ir[],
            Integer NAG_HUGE ic[], Boolean NAG_HUGE * sing, Integer NAG_HUGE ica[], double NAG_HUGE aux[],
            Integer NAG_HUGE * mtnmax, Integer nmax, Integer NAG_HUGE * mmax2);
#else
extern void d02rat();
#endif

#ifdef NAG_PROTO
extern void d02rau(Integer neq, Integer NAG_HUGE * np, Integer p, Integer r,
		   double NAG_HUGE x[], double NAG_HUGE y[], Integer iy,
		   
		   NAG_D02RAZ_fcn fcn,
		   NAG_D02RAZ_g g,
		   NAG_D02RAZ_fcnep fcnep,
		   NAG_D02RAZ_jacobg jacobg,
		   NAG_D02RAZ_jacobe jacobe,            
		   double NAG_HUGE b20[], 
		   NAG_D02RAZ_fcna fcna,
		   
		   double NAG_HUGE a1[], double NAG_HUGE b1[],
		   double NAG_HUGE a[], double NAG_HUGE c[],
		   double NAG_HUGE del[],
		   Boolean NAG_HUGE * casi, Boolean NAG_HUGE * sing, 
		   Integer NAG_HUGE ir[], Integer NAG_HUGE ic[],
		   double NAG_HUGE uu[], double NAG_HUGE res[],
		   Integer NAG_HUGE * lin,
		   Integer NAG_HUGE * mtnmax, Integer nmax, Integer NAG_HUGE * mmax2, 
		   double NAG_HUGE hx[],
		   double NAG_HUGE gradf[], double NAG_HUGE aux[], 
		   Integer NAG_HUGE ica[],
		   double NAG_HUGE xau[],
		   double NAG_HUGE *epsnu, 
		   double NAG_HUGE hmax[],
		   Integer NAG_HUGE ig[], 
		   Integer nig,
		   double NAG_HUGE h[],
		   Integer NAG_HUGE irn[],
		   Integer NAG_HUGE ip1[],
		   Integer nip, 
		   double NAG_HUGE f[], double NAG_HUGE alpha[],
		   Integer NAG_HUGE * iflag, Integer lp,
		   Nag_User NAG_HUGE *comm,
		   Nag_ad02ra *ad02ra);
#else
extern void d02rau();
#endif

#ifdef NAG_PROTO
extern void d02rav(Integer neq, Integer NAG_HUGE * np, Integer p, Integer r,
            double NAG_HUGE alpha[], double NAG_HUGE a1[], double NAG_HUGE b1[], double NAG_HUGE x[],
            double NAG_HUGE y[], Integer iy, double NAG_HUGE a2[], double NAG_HUGE c2[],
            double NAG_HUGE del[], 

            NAG_D02RAZ_fcn fcn,
            NAG_D02RAZ_g g,
            NAG_D02RAZ_fcnep fcnep,
            NAG_D02RAZ_fcna fcna,
            NAG_D02RAZ_fcnb fcnb,
            NAG_D02RAZ_jacobe jacobe,
            NAG_D02RAZ_jacobg jacobg,

            double NAG_HUGE a10[], double NAG_HUGE b10[],
            double NAG_HUGE gam[], double NAG_HUGE a20[], double NAG_HUGE b20[],
            Integer NAG_HUGE * jerror, double NAG_HUGE *eps, Integer NAG_HUGE ir[], Integer NAG_HUGE ic[],
            double NAG_HUGE uu[], double NAG_HUGE res[], Integer NAG_HUGE * mtnmax, Integer nmax,
            Integer NAG_HUGE * mmax2, double NAG_HUGE f[], double NAG_HUGE hx[], double NAG_HUGE sk[],
            double NAG_HUGE gradf[], double NAG_HUGE aux[], Integer NAG_HUGE ica[],
            double NAG_HUGE xau[], Integer lp, Integer mp, Integer NAG_HUGE * lin,
            double NAG_HUGE *epsnu, Integer NAG_HUGE * nu, Integer NAG_HUGE * inwt, Boolean NAG_HUGE * casi,
            double NAG_HUGE hmax[], Integer NAG_HUGE ig[], Integer nig, double NAG_HUGE h[],
            Integer NAG_HUGE irn[], Integer NAG_HUGE ip[], Integer nip, Nag_User NAG_HUGE *comm, Nag_ad02ra *ad02ra);
#else
extern void d02rav();
#endif

#ifdef NAG_PROTO
extern void d02raw(Integer NAG_HUGE * i0, Integer NAG_HUGE * n, Integer NAG_HUGE * np, double NAG_HUGE c[],
            double NAG_HUGE bb[], double NAG_HUGE x[], Integer nmax,
            
            double xbar, double NAG_HUGE alf[]);
#else
extern void d02raw();
#endif

#ifdef NAG_PROTO
extern void d02rax(Integer k, Integer p, Integer q, Integer n, Integer m,
            double NAG_HUGE a[], double NAG_HUGE x[], Integer nmax, double NAG_HUGE y[],
            double NAG_HUGE s[], Integer NAG_HUGE * mtnmax, double NAG_HUGE alf[],
            double NAG_HUGE c[], Integer NAG_HUGE * ierror);
#else
extern void d02rax();
#endif

#ifdef NAG_PROTO
extern void d02ray(Integer neq, Integer nmax, Integer NAG_HUGE * n, Integer p, Integer r,
            Integer NAG_HUGE * mtnmax, Integer NAG_HUGE * mmax2,
            double NAG_HUGE *a, double NAG_HUGE *b, double tol,
            double NAG_HUGE x[], double NAG_HUGE y[],
            Integer iy,
            double NAG_HUGE abt[],
            NAG_D02RAZ_fcn fcn,
            NAG_D02RAZ_g g,
            NAG_D02RAZ_fcnep fcnep,
            NAG_D02RAZ_fcna fcna,
            NAG_D02RAZ_fcnb fcnb,
            NAG_D02RAZ_jacobe jacobe,
            NAG_D02RAZ_jacobg jacobg,
            NAG_D02RAZ_jaceps jaceps,
            NAG_D02RAZ_jacgep jacgep,
            double NAG_HUGE a10[], double NAG_HUGE b10[], double NAG_HUGE gam[],
            double NAG_HUGE a20[], double NAG_HUGE b20[], double NAG_HUGE alpha[],
            double NAG_HUGE a1[], double NAG_HUGE b1[], double NAG_HUGE ej[], double NAG_HUGE a2[],
            double NAG_HUGE c2[], double NAG_HUGE del[], double NAG_HUGE uu[], double NAG_HUGE res[],
            double NAG_HUGE f[], double NAG_HUGE hx[], double NAG_HUGE sk[], double NAG_HUGE gradf[],
            double NAG_HUGE aux[], double NAG_HUGE xau[], 
            Integer NAG_HUGE ic[], Integer NAG_HUGE ir[], Integer NAG_HUGE iqj[], Integer NAG_HUGE ica[],
            double NAG_HUGE *deleps, Integer lp,
            Integer mp, Integer init, Integer NAG_HUGE * lin,
            double NAG_HUGE h2[], double NAG_HUGE hmax[], 
            Integer NAG_HUGE ig[], 
            Integer nig,
            Integer NAG_HUGE ip2[],
            Integer nip,
            Integer NAG_HUGE irn[], Integer NAG_HUGE * iflag,
		   Nag_User NAG_HUGE *comm,
		   Nag_ad02ra *ad02ra);
#else
extern void d02ray();
#endif

#ifdef NAG_PROTO
extern void d02raz(Integer neq, Integer nmax, Integer NAG_HUGE * np, Integer numbeg,
            Integer nummix, double NAG_HUGE *a, double NAG_HUGE *b, double tol,
            double NAG_HUGE x[], double NAG_HUGE y[], Integer iy, double NAG_HUGE abt[],
            NAG_D02RAZ_fcn fcn,
            NAG_D02RAZ_g g,
            NAG_D02RAZ_fcnep fcnep,
            NAG_D02RAZ_fcna fcna,
            NAG_D02RAZ_fcnb fcnb,
            NAG_D02RAZ_jacobe jacobe,
            NAG_D02RAZ_jacobg jacobg,
            NAG_D02RAZ_jaceps jaceps,
            NAG_D02RAZ_jacgep jacgep,
            double NAG_HUGE a1[], double NAG_HUGE b1[],
            double NAG_HUGE c1[], double NAG_HUGE d1[], double NAG_HUGE gam[], double NAG_HUGE a2[],
            double NAG_HUGE b2[], double NAG_HUGE work[], Integer NAG_HUGE * lwork, Integer NAG_HUGE * iwork,
            Integer liwork, double NAG_HUGE *deleps, Integer lp, Integer mp,
            Integer init, Integer NAG_HUGE * lin, Integer NAG_HUGE * iflag, Nag_User NAG_HUGE *comm, Nag_ad02ra *ad02ra);
#else
extern void d02raz();
#endif

#ifdef NAG_PROTO
extern void d02xjy(double t, Integer k, double NAG_HUGE *yh, Integer nyh, double NAG_HUGE dky[], 
            Integer NAG_HUGE *iflag, Integer neq, double NAG_HUGE *h, double NAG_HUGE *tn, double NAG_HUGE *hu,
            Integer NAG_HUGE *nq, double NAG_HUGE *odcode, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02xjy();
#endif

#ifdef NAG_PROTO
extern void d02xjz(double t, Integer k, double NAG_HUGE *yh, Integer nyh, double NAG_HUGE dky[],
            Integer NAG_HUGE *iflag, Integer neq, double NAG_HUGE *h, double NAG_HUGE *tn, double NAG_HUGE *hu,
            Integer NAG_HUGE *nq, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02xjz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d02xkc(double xsol, double NAG_HUGE sol[], Integer m, double NAG_HUGE *w, Integer neqmax,
            Integer iw, double NAG_HUGE w2[], Integer neq, double x, Integer nq,
            double hu, double h, Integer NAG_HUGE *ifail, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02xkc();
#endif

#ifdef NAG_PROTO
extern void d02xkz(double t, Integer k, double NAG_HUGE *yh, Integer nyh, double NAG_HUGE dky[], 
            Integer NAG_HUGE *iflag, Integer neq, double NAG_HUGE acor[], double h, double tn,
            double hu, Integer nq, Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern void d02xkz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL d02zac(Integer n, double NAG_HUGE v[], double NAG_HUGE w[], Integer NAG_HUGE *ifail,
              Nag_ODE_BDF NAG_HUGE *intern_comm);
#else
extern double d02zac();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGD02 */
