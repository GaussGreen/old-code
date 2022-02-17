#ifndef NAGG13
#define NAGG13
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg13.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g13 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2162 (Feb 1998).
 */

/* Message codes and extern of message list */
#include <nag_g13mesg.h>

#ifdef NAG_PROTO
extern void g13aac(double NAG_HUGE x[], Integer nx, Integer nd, Integer nds, Integer ns,
            double NAG_HUGE xd[], Integer NAG_HUGE *nxd, NagError NAG_HUGE *fail);
#else
extern void g13aac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13abc(double NAG_HUGE x[], Integer nx, Integer nk, double NAG_HUGE *mean, double NAG_HUGE *var,
            double NAG_HUGE r[], double NAG_HUGE *stat, NagError NAG_HUGE *fail);
#else
extern void g13abc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13acc(double NAG_HUGE r[], Integer nk, Integer nl, double NAG_HUGE p[], double NAG_HUGE v[],
            double NAG_HUGE ar[], Integer NAG_HUGE *nvl, NagError NAG_HUGE *fail);
#else
extern void g13acc();
#endif

#ifdef NAG_PROTO
extern void g13aeu(Integer id, double NAG_HUGE *ex, double NAG_HUGE *alpha, double NAG_HUGE *a, Integer na,
            double NAG_HUGE *w, double NAG_HUGE *beta, double NAG_HUGE *b, double NAG_HUGE *phi,
            double NAG_HUGE *theta, double NAG_HUGE *sphi, double NAG_HUGE *stheta, Integer np,
            Integer nq, Integer nps, Integer nqs, Integer ns, Integer npd);
#else
extern void g13aeu();
#endif

#ifdef NAG_PROTO
extern void g13aex(double NAG_HUGE *zb, Integer nzb, double eps, double NAG_HUGE *pg, Integer NAG_HUGE *kc);
#else
extern void g13aex();
#endif

#ifdef NAG_PROTO
extern void g13aey(Integer NAG_HUGE *mpqs, double NAG_HUGE *pa, Integer kwph,
            Integer npar, double NAG_HUGE *wb, double ef, Integer NAG_HUGE *mc,
            Integer NAG_HUGE *ierr);
#else
extern void g13aey();
#endif

#ifdef NAG_PROTO
extern void g13aez(double NAG_HUGE *par, Integer npar, Integer NAG_HUGE *mpqs, double NAG_HUGE *wa,
            Integer kwph);
#else
extern void g13aez();
#endif

#ifdef NAG_PROTO
extern void g13ahx(double NAG_HUGE *aex, Integer nst, Integer nd, Integer nds,
            Integer ns, Integer kdr, double NAG_HUGE *uv, double NAG_HUGE *dv);
#else
extern void g13ahx();
#endif

#ifdef NAG_PROTO
extern void g13ahy(double NAG_HUGE *st, Integer nst, Integer np, Integer nd,
            Integer nq, Integer nps, Integer nds, Integer nqs,
            Integer ns, double NAG_HUGE *aex, double NAG_HUGE *aal, double NAG_HUGE *aexr);
#else
extern void g13ahy();
#endif

#ifdef NAG_PROTO
extern void g13ajy(Integer NAG_HUGE *mr, Integer NAG_HUGE *np, Integer NAG_HUGE *nd, Integer NAG_HUGE *nq,
            Integer NAG_HUGE *nps, Integer NAG_HUGE *nds, Integer NAG_HUGE *nqs, Integer NAG_HUGE *ns,
            Integer NAG_HUGE *npd, Integer NAG_HUGE *ndd, Integer NAG_HUGE *nqd, Integer NAG_HUGE *mpqs,
            Integer NAG_HUGE *npar);
#else
extern void g13ajy();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13bec(Nag_ArimaOrder NAG_HUGE *arimav, Integer nseries, Nag_TransfOrder NAG_HUGE *transfv,
            double NAG_HUGE para[], Integer npara, Integer nxxy, double NAG_HUGE xxy[],
            Integer tdxxy, double NAG_HUGE sd[], double NAG_HUGE *rss, double NAG_HUGE *objf,
            double NAG_HUGE *df, Nag_G13_Opt NAG_HUGE *options, NagError NAG_HUGE *fail);
#else
extern void g13bec();
#endif

#ifdef NAG_PROTO
extern void g13bej(double NAG_HUGE *xxy, Integer tdxxy, Nag_G13_Opt NAG_HUGE *opt, Integer n,
            Integer nxsp, double NAG_HUGE *pxs, Integer ipxs,
            double NAG_HUGE *wds, Integer iwds, Integer nfr, Integer NAG_HUGE *mt,
            Integer tdmt, double NAG_HUGE *w, double NAG_HUGE *sttf, Integer isttf,
            Integer NAG_HUGE *nsttf, Integer NAG_HUGE *mr, double c, Integer kzef,
            double NAG_HUGE *bf, Integer nbf, Integer kss, Integer NAG_HUGE *ierr);
#else
extern void g13bej();
#endif

#ifdef NAG_PROTO
extern void g13bek(double NAG_HUGE *a, double NAG_HUGE *b, double NAG_HUGE *apb, Integer na, Integer nb,
            Integer NAG_HUGE *napb);
#else
extern void g13bek();
#endif

#ifdef NAG_PROTO
extern void g13bel(double NAG_HUGE *x, Integer n, double NAG_HUGE *px, double NAG_HUGE *wd, Integer nnb,
            Integer nnp, Integer nnq, Integer nnr, Integer npx, Integer NAG_HUGE *nex,
            double NAG_HUGE *ex, double NAG_HUGE *ez);
#else
extern void g13bel();
#endif

#ifdef NAG_PROTO
extern void g13bem(Integer n4, Integer NAG_HUGE *mis, Integer NAG_HUGE *mrn, Integer NAG_HUGE *mord,
            Integer NAG_HUGE *mtyp, Integer NAG_HUGE *mser, Integer nxs);
#else
extern void g13bem();
#endif

#ifdef NAG_PROTO
extern void g13ben(Integer NAG_HUGE *mqab, Integer NAG_HUGE *mqsg, Integer ja, Integer jb,
            Integer nxsp, Integer krn, Integer NAG_HUGE *msn, Integer NAG_HUGE *mspa,
            Integer NAG_HUGE *mspb, Integer NAG_HUGE *k, Integer kspa, Integer kspb,
            Integer NAG_HUGE *mis, Integer NAG_HUGE *mrn);
#else
extern void g13ben();
#endif

#ifdef NAG_PROTO
extern void g13bep(Integer NAG_HUGE *mpab, Integer NAG_HUGE *mpsg, Integer krn, Integer NAG_HUGE *msn,
            Integer NAG_HUGE *mspa, Integer NAG_HUGE *mspb, Integer NAG_HUGE *k, Integer kspa,
            Integer kspb, Integer kspc, Integer NAG_HUGE *mis, Integer NAG_HUGE *mrn);
#else
extern void g13bep();
#endif

#ifdef NAG_PROTO
extern void g13beq(double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE *r, Integer npd, double NAG_HUGE *v);
#else
extern void g13beq();
#endif

#ifdef NAG_PROTO
extern void g13ber(double NAG_HUGE *beta, Integer NAG_HUGE *mop, Integer NAG_HUGE *msn, Integer NAG_HUGE *mqab,
            Integer NAG_HUGE *mqsg, Integer nxsp, double eps, double NAG_HUGE *wd,
            Integer NAG_HUGE *ksfs);
#else
extern void g13ber();
#endif

#ifdef NAG_PROTO
extern void g13bes(Integer NAG_HUGE *mpab, Integer nms, Integer nmsq, Integer numw,
            Integer NAG_HUGE *msn, Integer NAG_HUGE *mspa, Integer NAG_HUGE *mspb, double NAG_HUGE *par,
            Integer npar, Integer kef, Integer nas, Integer nbs,
            Integer nbvd, double NAG_HUGE *h, Integer tdh, double NAG_HUGE *wq, Integer tdwq,
            double NAG_HUGE *wz, double NAG_HUGE *a, double NAG_HUGE *b, double NAG_HUGE *gca, double NAG_HUGE *delta,
            Integer np, Integer nq, Integer nps, Integer ns, Integer npd,
            double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE *r, double s, Integer ndfva,
            Integer kfl, double NAG_HUGE *g, Integer ngh, double NAG_HUGE *phi, double NAG_HUGE *sphi);
#else
extern void g13bes();
#endif

#ifdef NAG_PROTO
extern void g13bet(double NAG_HUGE *a, Integer tda, Integer n, double NAG_HUGE *d);
#else
extern void g13bet();
#endif

#ifdef NAG_PROTO
extern void g13beu(Integer lac, Integer lbc, Integer NAG_HUGE *msn, Integer NAG_HUGE *maspa,
            Integer NAG_HUGE *maspb, Integer nmsu, double NAG_HUGE *a, double NAG_HUGE *b,
            Integer nbvd, double NAG_HUGE *s, double NAG_HUGE *g, double NAG_HUGE *h,
            Integer tdh, Integer NAG_HUGE *mpab);
#else
extern void g13beu();
#endif

#ifdef NAG_PROTO
extern void g13bev(double NAG_HUGE *xxy, Integer tdxxy, Integer n, Integer NAG_HUGE *mr,
            double NAG_HUGE *para, double c, Integer kfc, Integer nxsp,
            double NAG_HUGE *pxs, Integer ipxs, double NAG_HUGE *bf, double NAG_HUGE *wds,
            Integer iwds, Integer nbv, Integer kef, Integer NAG_HUGE *mt,
            Integer tdmt, Integer NAG_HUGE *mpab, Integer NAG_HUGE *mqab, double NAG_HUGE *a,
            Integer ida, double NAG_HUGE *b, Integer idb, double NAG_HUGE *w,
            Integer idw, Integer kab, Integer nbf, Integer lbc);
#else
extern void g13bev();
#endif

#ifdef NAG_PROTO
extern void g13bew(double NAG_HUGE *beta, Integer nbeta, double NAG_HUGE *para, Integer npar,
            double NAG_HUGE *bf, double NAG_HUGE *pxs, Integer ipxs,
            double NAG_HUGE *wds, Integer iwds, double NAG_HUGE *c, Integer kef, Integer NAG_HUGE *mop,
            Integer NAG_HUGE *mis, Integer NAG_HUGE *mrn, Integer kfbs, Integer klbs);
#else
extern void g13bew();
#endif

#ifdef NAG_PROTO
extern void g13bex(Integer NAG_HUGE *mt, Integer tdmt, Integer kcol, Integer nxsp,
            Integer NAG_HUGE *nnb, Integer NAG_HUGE *nnp, Integer NAG_HUGE *nnq, Integer NAG_HUGE *nnr,
            Integer NAG_HUGE *nwd, Integer NAG_HUGE *ngw, Integer NAG_HUGE *npx);
#else
extern void g13bex();
#endif

#ifdef NAG_PROTO
extern void g13bey(Integer NAG_HUGE *mpqs, Integer nxsp, Integer nbf, Integer kfc,
            Integer kef, Integer NAG_HUGE *mt, Integer tdmt, Integer NAG_HUGE *mpab,
            Integer NAG_HUGE *mqab, Integer tdmqab, Integer NAG_HUGE *kcp, Integer NAG_HUGE *mpsg,
            Integer NAG_HUGE *mqsg, Integer tdmqsg, Integer np, Integer nps,
            Integer ns, Integer NAG_HUGE *msn, Integer NAG_HUGE *mspa, Integer NAG_HUGE *mspb,
            Integer NAG_HUGE *mop, Integer NAG_HUGE *mis, Integer NAG_HUGE *mrn);
#else
extern void g13bey();
#endif

#ifdef NAG_PROTO
extern void g13bez(Integer nxs, Integer NAG_HUGE *mt, Integer tdmt, Integer np,
            Integer nq, Integer nps, Integer nqs, Integer ns,
            Integer ndd, Integer nqd, Integer kef, Integer NAG_HUGE *nbss,
            Integer NAG_HUGE *nbv, Integer NAG_HUGE *nbf, Integer nxsp);
#else
extern void g13bez();
#endif

#ifdef NAG_PROTO
extern void g13bfr(double NAG_HUGE *ex, double NAG_HUGE *alpha, double NAG_HUGE *a, double NAG_HUGE *pa,
            Integer na, double NAG_HUGE *w, double NAG_HUGE *beta, double NAG_HUGE *b, double NAG_HUGE *wa,
            Integer nrmp, Integer kdq, Integer np, Integer nq,
            Integer nps, Integer nqs, Integer ns, Integer npd,
            Integer NAG_HUGE *mpqs, Integer nas, Integer nbs);
#else
extern void g13bfr();
#endif

#ifdef NAG_PROTO
extern void g13bfs(double NAG_HUGE *f, double NAG_HUGE *r, double NAG_HUGE *alpha, Integer npd,
            double NAG_HUGE *phi, Integer np, double NAG_HUGE *sphi, Integer nps,
            double NAG_HUGE *dva, Integer iddv, Integer ns,  double v);
#else
extern void g13bfs();
#endif

#ifdef NAG_PROTO
extern void g13bft(Integer NAG_HUGE *mpab, Integer nms, Integer num, Integer numw,
            Integer NAG_HUGE *msn, Integer NAG_HUGE *mspa, Integer NAG_HUGE *mspb, Integer npar,
            Integer nas, Integer nbs, Integer nbvd,  double NAG_HUGE *h,
            Integer tdh, double NAG_HUGE *wq, Integer tdwq, double NAG_HUGE *wz,
            double NAG_HUGE *a, double NAG_HUGE *b, double NAG_HUGE *gc, NagError NAG_HUGE *fail);
#else
extern void g13bft();
#endif

#ifdef NAG_PROTO
extern void g13bfu(double NAG_HUGE *f, double NAG_HUGE *r, double NAG_HUGE *alpha, Integer npd, Integer ks,
            double v);
#else
extern void g13bfu();
#endif

#ifdef NAG_PROTO
extern void g13bfv(double NAG_HUGE *phi, Integer np, double NAG_HUGE *sphi, Integer nps,
            double NAG_HUGE *f, double NAG_HUGE *alpha, Integer npd, Integer ns,
            double NAG_HUGE *delta);
#else
extern void g13bfv();
#endif

#ifdef NAG_PROTO
extern void g13bfw(double NAG_HUGE *h, Integer tdh, Integer i, Integer j,
            Integer ksgn, double NAG_HUGE *f, Integer NAG_HUGE *kfl,
            Integer NAG_HUGE *msn, Integer NAG_HUGE *ix, Integer NAG_HUGE *jx);
#else
extern void g13bfw();
#endif

#ifdef NAG_PROTO
extern void g13bfx(Integer itc, double s, double d, double NAG_HUGE *para,
            Integer npara, double NAG_HUGE *sd, double df,
            Integer npe, Integer NAG_HUGE *mtyp, Integer NAG_HUGE *mser, Nag_Comm NAG_HUGE *comm,
            Nag_FileSt NAG_HUGE *stream, Nag_G13_Opt NAG_HUGE *opt);
#else
extern void g13bfx();
#endif

#ifdef NAG_PROTO
extern void g13bfy(Integer NAG_HUGE *mr, Integer npara, double NAG_HUGE *para, Integer npe,
            Integer nxy, double NAG_HUGE *xxy, Integer tdxxy,
            Integer nxsp, Integer NAG_HUGE *mt, Integer tdmt, Integer kef,
            Integer nit, double NAG_HUGE *zsp, Integer NAG_HUGE *itc, double NAG_HUGE *sd,
            double df, double NAG_HUGE *cm,
            Integer tdcm, double NAG_HUGE *s, double NAG_HUGE *d, Integer kzef, double NAG_HUGE *res,
            Integer NAG_HUGE *mwa, double NAG_HUGE *beta, double NAG_HUGE *betac,
            Integer ngh, double NAG_HUGE *g, double NAG_HUGE *gc, Integer ih, double NAG_HUGE *h,
            Integer tdh, double NAG_HUGE *hc, Integer tdhc, double NAG_HUGE *delt,
            Integer iwds, double NAG_HUGE *wds, double NAG_HUGE *bf,
            Integer ipxs, double NAG_HUGE *pxs,
            double NAG_HUGE *c, double NAG_HUGE *sttf,
            Integer isttf, Integer NAG_HUGE *nsttf, Integer NAG_HUGE *ierr,
            Nag_Comm NAG_HUGE *comm, Nag_FileSt NAG_HUGE *stream,
            Nag_G13_Opt NAG_HUGE *opt);
#else
extern void g13bfy();
#endif

#ifdef NAG_PROTO
extern void g13bhx(double NAG_HUGE *sttf, Integer kpst, Integer nwd, double NAG_HUGE *xc,
            Integer lxc, double NAG_HUGE *zc, Integer lzc, double NAG_HUGE *para,
            Integer kppa, double NAG_HUGE *xn, double NAG_HUGE *zn, Integer nnv,
            Integer kustf);
#else
extern void g13bhx();
#endif

#ifdef NAG_PROTO
extern void g13bhy(double NAG_HUGE *st, Integer nst, Integer np, Integer nd, Integer nq,
            Integer nps, Integer nds, Integer nqs, Integer ns,
            double NAG_HUGE *phi, double NAG_HUGE *theta, double NAG_HUGE *sphi, double NAG_HUGE *stheta,
            double NAG_HUGE *c, double rms, Integer nfv, double NAG_HUGE *fva,
            double NAG_HUGE *fsd, Integer kfva, double NAG_HUGE *aex, double NAG_HUGE *aal,
            double NAG_HUGE *aexr);
#else
extern void g13bhy();
#endif

#ifdef NAG_PROTO
extern void g13bhz(Integer NAG_HUGE *mr, double NAG_HUGE *para, Integer npara, double NAG_HUGE *sttf,
            Integer NAG_HUGE *nsttf, Integer NAG_HUGE *mt, Integer tdmt, Integer nxsp,
            Integer NAG_HUGE *mrx, Integer tdmrx, double NAG_HUGE *parx, Integer tdparx,
            double NAG_HUGE *xxyn, Integer tdxxyn, Nag_G13_Opt NAG_HUGE *opt,
            double NAG_HUGE *rmsxy, Integer nfv, Integer kzef, double NAG_HUGE *fva,
            double NAG_HUGE *fsd, double NAG_HUGE *psi, double NAG_HUGE *aex, double NAG_HUGE *aal, double NAG_HUGE *zn,
            Integer ips, double NAG_HUGE *wa, Integer iwa, Integer npxy);
#else
extern void g13bhz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13bjc(Nag_ArimaOrder NAG_HUGE *arimav, Integer nseries, Nag_TransfOrder NAG_HUGE *transfv,
            double NAG_HUGE para[], Integer npara, Integer nev, Integer nfv,
            double NAG_HUGE xxy[], Integer tdxxy, double NAG_HUGE rmsxy[],
            Integer NAG_HUGE mrx[], Integer tdmrx, double NAG_HUGE parx[], Integer ldparx,
            Integer tdparx, double NAG_HUGE fva[], double NAG_HUGE fsd[], Nag_G13_Opt NAG_HUGE *options,
            NagError NAG_HUGE *fail);
#else
extern void g13bjc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13bxc(Nag_G13_Opt NAG_HUGE *options);
#else
extern void g13bxc();
#endif

#ifdef NAG_PROTO
extern void g13bxz(Nag_ArimaOrder NAG_HUGE *arimav, Nag_TransfOrder NAG_HUGE *transfv, Integer nseries,
            Integer NAG_HUGE mr[], Integer NAG_HUGE mt[]);
#else
extern void g13bxz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13byc(Integer nseries, Nag_TransfOrder NAG_HUGE *transfv, NagError NAG_HUGE *fail);
#else
extern void g13byc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13bzc(Nag_TransfOrder NAG_HUGE *transfv);
#else
extern void g13bzc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cac(Integer nx, Integer mtx,  double px,  Nag_LagWindow lag_win,
             Integer mw, Integer ic, Integer nc,  double NAG_HUGE c[],
             Integer kc, Integer l, Integer lg, Integer nxg,
             double NAG_HUGE xg[],  Integer NAG_HUGE *ng,  double NAG_HUGE stats[],
             NagError NAG_HUGE *fail);
#else
extern void g13cac();
#endif

#ifdef NAG_PROTO
extern double g13caw(Integer iw,  double s3, double pi);
#else
extern double g13caw();
#endif

#ifdef NAG_PROTO
extern void g13cax(double NAG_HUGE xg[],  Integer nx, Integer m,  double pi);
#else
extern void g13cax();
#endif

#ifdef NAG_PROTO
extern void g13cay(double NAG_HUGE xg[],  Integer nx, Integer mtx);
#else
extern void g13cay();
#endif

#ifdef NAG_PROTO
extern void g13caz(Integer nx, Integer mtx,  double px,  Integer iw,
             Integer mw, Integer ic, Integer nc,  double NAG_HUGE c[],
             Integer kc, Integer l, Integer lg, Integer nxg,
             double NAG_HUGE xg[],  Integer NAG_HUGE *ng,  double NAG_HUGE stats[],
             Integer NAG_HUGE *ierror);
#else
extern void g13caz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cbc(Integer nx, NagMeanOrTrend mt_correction,  double px,  Integer mw,
             double pw,  Integer l, Integer kc, Nag_LoggedSpectra lg_spect,
             double NAG_HUGE x[], double NAG_HUGE *NAG_HUGE *g,  Integer NAG_HUGE *ng,  double NAG_HUGE stats[],
             NagError NAG_HUGE *fail);
#else
extern void g13cbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13ccc(Integer nxy, Integer mtxy,  double pxy,  Integer iw,
             Integer mw, Integer is, Integer ic, Integer nc,
             double NAG_HUGE cxy[], double NAG_HUGE cyx[],  Integer kc, Integer l,
             Integer nxyg,  double NAG_HUGE xg[], double NAG_HUGE yg[],  Integer NAG_HUGE *ng,
             NagError NAG_HUGE *fail);
#else
extern void g13ccc();
#endif

#ifdef NAG_PROTO
extern void g13ccy(double NAG_HUGE x[],  Integer k, Integer m,  double NAG_HUGE *s1,
             double NAG_HUGE *s2);
#else
extern void g13ccy();
#endif

#ifdef NAG_PROTO
extern void g13ccz(Integer nxy, Integer mtxy,  double pxy,  Integer iw,
            Integer mw, Integer is, Integer ic, Integer nc,
            double NAG_HUGE cxy[], double NAG_HUGE cyx[],  Integer kc, Integer l,
            Integer nxyg,  double NAG_HUGE xg[], double NAG_HUGE yg[],  Integer NAG_HUGE *ng,
            Integer NAG_HUGE *ierror);
#else
extern void g13ccz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cdc(Integer nxy, NagMeanOrTrend mt_correction, double pxy,  Integer mw,
             Integer is,  double pw,  Integer l, Integer kc,
             double NAG_HUGE x[], double NAG_HUGE y[],  Complex NAG_HUGE *NAG_HUGE *g, Integer NAG_HUGE *ng, NagError NAG_HUGE *fail);
#else
extern void g13cdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cec(double NAG_HUGE xg[], double NAG_HUGE yg[], Complex NAG_HUGE xyg[],  Integer ng,
            double NAG_HUGE stats[], double NAG_HUGE ca[], double NAG_HUGE calw[], double NAG_HUGE caup[],
            double NAG_HUGE *t, double NAG_HUGE sc[], double NAG_HUGE sclw[], double NAG_HUGE scup[],
            NagError NAG_HUGE *fail);
#else
extern void g13cec();
#endif

#ifdef NAG_PROTO
extern void g13cez(double NAG_HUGE xg[], double NAG_HUGE yg[], double NAG_HUGE xyrg[],
             double NAG_HUGE xyig[],  Integer ng,  double NAG_HUGE stats[],
             double NAG_HUGE ca[], double NAG_HUGE calw[], double NAG_HUGE caup[],
             double NAG_HUGE *t, double NAG_HUGE sc[], double NAG_HUGE sclw[],
             double NAG_HUGE scup[],  Integer NAG_HUGE *ierror);
#else
extern void g13cez();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cfc(double NAG_HUGE xg[], double NAG_HUGE yg[], Complex NAG_HUGE xyg[], Integer ng,
            double NAG_HUGE stats[], double NAG_HUGE gn[], double NAG_HUGE gnlw[], double NAG_HUGE gnup[],
            double NAG_HUGE ph[], double NAG_HUGE phlw[], double NAG_HUGE phup[],
            NagError NAG_HUGE *fail);
#else
extern void g13cfc();
#endif

#ifdef NAG_PROTO
extern void g13cfz(double NAG_HUGE xg[], double NAG_HUGE yg[], double NAG_HUGE xyrg[],
             double NAG_HUGE xyig[],  Integer ng,  double NAG_HUGE stats[],
             double NAG_HUGE gn[], double NAG_HUGE gnlw[], double NAG_HUGE gnup[],
             double NAG_HUGE ph[], double NAG_HUGE phlw[], double NAG_HUGE phup[],
             Integer NAG_HUGE *ierror);
#else
extern void g13cfz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13cgc(double NAG_HUGE xg[], double NAG_HUGE yg[], Complex NAG_HUGE xyg[], Integer ng,
            double NAG_HUGE stats[],  Integer l, Integer n,  double NAG_HUGE er[],
            double NAG_HUGE *erlw, double NAG_HUGE *erup, double NAG_HUGE rf[], double NAG_HUGE *rfse,
            NagError NAG_HUGE *fail);
#else
extern void g13cgc();
#endif

#ifdef NAG_PROTO
extern void g13cgz(double NAG_HUGE xg[], double NAG_HUGE yg[], double NAG_HUGE xyrg[],
             double NAG_HUGE xyig[],  Integer ng,  double NAG_HUGE stats[],  Integer l,
             Integer n,  double NAG_HUGE er[], double NAG_HUGE *erlw, double NAG_HUGE *erup,
             double NAG_HUGE rf[], double NAG_HUGE *rfse,  Integer NAG_HUGE *ierror);
#else
extern void g13cgz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13eac(Integer n, Integer m, Integer p, double NAG_HUGE s[], Integer tds,
            double NAG_HUGE a[], Integer tda, double NAG_HUGE b[], Integer tdb,
            double NAG_HUGE q[], Integer tdq, double NAG_HUGE c[], Integer tdc,
            double NAG_HUGE r[], Integer tdr, double NAG_HUGE k[], Integer tdk, 
            double NAG_HUGE h[], Integer tdh, double tol, NagError NAG_HUGE *fail);
#else
extern void g13eac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13ebc(Integer n, Integer m, Integer p, double NAG_HUGE s[], Integer tds,
            double NAG_HUGE a[], Integer tda, double NAG_HUGE b[], Integer tdb, double NAG_HUGE q[],
            Integer tdq, double NAG_HUGE c[], Integer tdc, double NAG_HUGE r[], Integer tdr,
            double NAG_HUGE k[], Integer tdk, double NAG_HUGE h[], Integer tdh,
            double tol, NagError NAG_HUGE *fail);
#else
extern void g13ebc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13ecc(Integer n, Integer m, Integer p,
            Nag_ab_input inp_ab,double NAG_HUGE t[], Integer tdt,
            double NAG_HUGE ainv[], Integer tda, double NAG_HUGE b[], Integer tdb, 
            double NAG_HUGE rinv[], Integer tdr, double NAG_HUGE c[], Integer tdc,
            double NAG_HUGE qinv[], Integer tdq, double NAG_HUGE x[], double NAG_HUGE rinvy[], 
            double NAG_HUGE z[], double tol,
            NagError NAG_HUGE *fail);
#else
extern void g13ecc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13edc(Integer n, Integer m, Integer p, 
            double NAG_HUGE t[], Integer tdt, double NAG_HUGE ainv[], Integer tda,
            double NAG_HUGE ainvb[], Integer tdai, double NAG_HUGE rinv[], Integer tdr,
            double NAG_HUGE c[], Integer tdc, double NAG_HUGE qinv[], Integer tdq, double NAG_HUGE x[],
            double NAG_HUGE rinvy[], double NAG_HUGE z[], double tol,  NagError NAG_HUGE *fail);
#else
extern void g13edc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13ewc(Integer n, Integer p,  Nag_ObserverForm reduceto, double NAG_HUGE a[],
             Integer tda, double NAG_HUGE c[], Integer tdc, double NAG_HUGE u[], Integer tdu,
             NagError NAG_HUGE *fail);
#else
extern void g13ewc ();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13exc(Integer n, Integer m, Nag_ControllerForm reduceto, double NAG_HUGE a[],
            Integer tda, double NAG_HUGE b[], Integer tdb, double NAG_HUGE u[], Integer tdu, 
            NagError NAG_HUGE *fail);
#else
extern void g13exc();
#endif

#ifdef NAG_PROTO
extern void g13exy(MatrixType uplo, Integer n, double NAG_HUGE t[],
            Integer tdt, double NAG_HUGE *rcond, Integer NAG_HUGE *ifail);
     
#else
extern void g13exy();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g13xzc(Nag_G13_Opt NAG_HUGE *options);
#else
extern void g13xzc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG13 */
