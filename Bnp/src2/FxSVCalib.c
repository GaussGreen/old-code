// FxSVCalib.c : FXSV underlying initialization and calibration

#include "srt_h_all.h"
#include "FxSVCalib.h"
#include "math.h"
#include "Fx3FCalib.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nagd01.h"
#include "nagd02.h"
#include "nagf02.h"
#include "nage04.h"
#include "intde2.h"
#include "opfnctns.h"
#include "AffineModelGen.h"

// Memory allocation and copy construction:

Err InitFxSVUnd( SFxSVUndDesc *und, char *dom_und, char *for_und, double spot,
				long *dates, int ndates, double *beta, double *alpha, double *gamma, double ***rho )
{
	SrtUndPtr	domundptr = NULL;
	int			i, j, k;

	if (!und) return "Und must be allocated before calling InitFxSVUnd";
	domundptr = lookup_und(dom_und);
	if (!domundptr) return "Domestic underlying not found in InitFxSVUnd";
	und->today = get_today_from_underlying(domundptr);

	strcpy(und->dom_und, dom_und);
	strupper(und->dom_und);
	rem_tick_string(und->dom_und, und->dom_und);
	strcpy(und->for_und, for_und);
	strupper(und->for_und);
	rem_tick_string(und->for_und, und->for_und);
	und->spot = spot;
	und->ntimes = ndates;
	und->dates = (long *) calloc(ndates, sizeof(long));
	und->times = (double *) calloc(ndates, sizeof(double));
	und->beta = (double *) calloc(ndates, sizeof(double));
	und->alpha = (double *) calloc(ndates, sizeof(double));
	und->gamma = (double *) calloc(ndates, sizeof(double));
	und->rho = f3tensor(0, ndates-1, 0, 3, 0, 3);

	if (!und->dates || !und->times || !und->beta || !und->alpha || !und->gamma || !und->rho)
		return "Memory failure in InitFxSVUnd";

	memcpy(und->dates, dates, ndates * sizeof(long));
	for (i=0; i < ndates; i++) und->times[i] = (dates[i] - und->today) * YEARS_IN_DAY;

	memcpy(und->beta, beta, ndates * sizeof(double));
	memcpy(und->alpha, alpha, ndates * sizeof(double));
	memcpy(und->gamma, gamma, ndates * sizeof(double));

	for (i=0; i < ndates; i++) for (j=0; j < 4; j++) for (k=0; k < 4; k++) und->rho[i][j][k] = rho[i][j][k];

	return NULL;
}

// Destructor:
char *FreeFxSVUnd(void *ptr)
{
	SrtUndPtr		undptr = (SrtUndPtr) ptr;
	SFxSVUndDesc	*und = NULL;

	if (!ptr) return NULL;
	und = (SFxSVUndDesc *) undptr->spec_desc;
	if (!und) goto FREE_RETURN;

	free(und->dates);
	free(und->times);
	free(und->beta);
	free(und->alpha);
	free(und->gamma);
	if (und->rho) free_f3tensor(und->rho, 0, und->ntimes-1, 0, 3, 0, 3);

FREE_RETURN:
	free(und);
	free(ptr);
	return NULL;
}

Err FxSVFullFromUnd( SFxSVUndFull *o, SrtUndPtr und, long last_date )
{
	Err				err = NULL;
	double			tau_d, tau_f, T_last;

	o->pmdl = (SFxSVUndDesc *) und->spec_desc;

	if (!o->pmdl || strcmp(und->underl_lbl, "FXSV_UND"))
		return serror("Underlying %s improperly initialised or of incorrect type", und->underl_name);

	T_last = (last_date - o->pmdl->today) * YEARS_IN_DAY;

	err = Get_LGM_TermStructure2(o->pmdl->dom_und, &o->sigtms_d, &o->sig_d, &o->nsig_d, &tau_d);
	if (err) return err;
	o->lam_d = 1.0 / tau_d;

	err = Get_LGM_TermStructure2(o->pmdl->for_und, &o->sigtms_f, &o->sig_f, &o->nsig_f, &tau_f);
	if (err) return err;
	o->lam_f = 1.0 / tau_f;

	o->ntimes = o->nsig_d + o->nsig_f + o->pmdl->ntimes + 2;
	o->times = calloc(o->ntimes, sizeof(double));
	if (!o->times) return serror("Memory failure");

	memcpy(o->times, o->sigtms_d, o->nsig_d * sizeof(double));
	memcpy(o->times + o->nsig_d, o->sigtms_f, o->nsig_f * sizeof(double));
	memcpy(o->times + o->nsig_d + o->nsig_f, o->pmdl->times, o->pmdl->ntimes * sizeof(double));
	o->times[o->ntimes-2] = 0.0;
	o->times[o->ntimes-1] = T_last;

	num_f_sort_vector(o->ntimes, o->times);
	num_f_unique_vector(&o->ntimes, o->times);
	while (o->times[o->ntimes-1] > T_last + 1e-5) o->ntimes--;

	return NULL;
}

Err FxSVFullFree( SFxSVUndFull *o )
{
	free(o->sig_d);  free(o->sigtms_d);
	free(o->sig_f);  free(o->sigtms_f);
	free(o->times);
	return NULL;
}

// System of ODE - right hand part evaluation function:

static void NAG_CALL EvalDerivatives(Integer neq, double t, double y[], double yp[], Nag_User *comm)
{
	SFxSVComm_ODE	*p = (SFxSVComm_ODE *) comm->p;
	double			Gamma_d, Gamma_f;
	double			c1, c2, c3, c4, c5, c6, c7, c8;
	double			u_re2, u_im2, u_re_u_im;
	double			*yg[4], *ypg[4];
	int				i, offs;

	for (i=0, offs = 6; i < 4; i++) if (p->want[i])
	{
		yg[i] = y + offs;
		ypg[i] = yp + offs;
		offs += 6;
	}

	Gamma_d = -p->sig_d * ( 1.0 - exp(-p->lam_d*(p->T_pay - t)) ) / p->lam_d;
	Gamma_f = -p->sig_f * ( 1.0 - exp(-p->lam_f*(p->T_pay - t)) ) / p->lam_f;
	u_re2 = p->u_re * p->u_re;
	u_im2 = p->u_im * p->u_im;
	u_re_u_im = 2.0 * p->u_re * p->u_im;

	// dy1/dt, dy2/dt
	c1 = Gamma_d * Gamma_d / 2.0 + Gamma_f * Gamma_f / 2.0 - p->rho[0][1] * Gamma_d * Gamma_f;
	c2 = p->rho[0][3] * p->alpha * Gamma_d + p->gamma;
	c3 = p->alpha * p->alpha;
	c4 = p->alpha * (p->rho[1][3] * Gamma_f - p->rho[0][3] * Gamma_d);

	yp[0] = c1 * (u_re2 - u_im2 - p->u_im) - c2 * y[2] - c3 * y[4]
		+ c4 * (p->u_re * y[3] + p->u_im * y[2]) - c3 * (y[2]*y[2] - y[3]*y[3]) / 2.0;

	yp[1] = c1 * (p->u_re + u_re_u_im) - c2 * y[3] - c3 * y[5]
		- c4 * (p->u_re * y[2] - p->u_im * y[3]) - c3 * y[2] * y[3];

	for (i=0; i < 4; i++) if (p->want[i])
	{
		ypg[i][0] = - c2 * yg[i][2] - c3 * yg[i][4]
			+ c4 * (p->u_re * yg[i][3] + p->u_im * yg[i][2]) - c3 * (yg[i][2]*y[2] - yg[i][3]*y[3]);
		ypg[i][1] = - c2 * yg[i][3] - c3 * yg[i][5]
			- c4 * (p->u_re * yg[i][2] - p->u_im * yg[i][3]) - c3 * (y[2] * yg[i][3] + yg[i][2] * y[3]);
	}

	if (p->is_cur_idx && p->want[1])
	{
		c6 = p->rho[0][3] * Gamma_d;
		c7 = (p->rho[1][3] * Gamma_f - p->rho[0][3] * Gamma_d);

		ypg[1][0] += -c6 * y[2] - 2.0 * p->alpha * y[4] + c7 * (p->u_re * y[3] + p->u_im * y[2])
			- p->alpha * (y[2]*y[2] - y[3]*y[3]);

		ypg[1][1] += -c6 * y[3] - 2.0 * p->alpha * y[5] - c7 * (p->u_re * y[2] - p->u_im * y[3])
			- 2.0 * p->alpha * y[2] * y[3];
	}
	if (p->is_cur_idx && p->want[2])
	{
		ypg[2][0] += -y[2];
		ypg[2][1] += -y[3];
	}

	// dy3/dt, dy4/dt
	c1 = p->beta * (p->rho[1][2] * Gamma_f - p->rho[0][2] * Gamma_d);
	c2 *= 2.0;
	c3 *= 2.0;
	c4 *= 2.0;
	c5 = p->alpha * p->rho[2][3] * p->beta;

	yp[2] = c1 * (u_re2 - u_im2 - p->u_im) + p->gamma * y[2] - c2 * y[4]
		+ c5 * (p->u_re * y[3] + p->u_im * y[2]) + c4 * (p->u_re * y[5] + p->u_im * y[4])
		- c3 * (y[2] * y[4] - y[3] * y[5]);

	yp[3] = c1 * (p->u_re + u_re_u_im) + p->gamma * y[3] - c2 * y[5]
		- c5 * (p->u_re * y[2] - p->u_im * y[3]) - c4 * (p->u_re * y[4] - p->u_im * y[5])
		- c3 * (y[2] * y[5] + y[3] * y[4]);

	for (i=0; i < 4; i++) if (p->want[i])
	{
		ypg[i][2] = p->gamma * yg[i][2] - c2 * yg[i][4]
		+ c5 * (p->u_re * yg[i][3] + p->u_im * yg[i][2]) + c4 * (p->u_re * yg[i][5] + p->u_im * yg[i][4])
		- c3 * (y[2] * yg[i][4] - y[3] * yg[i][5] + yg[i][2] * y[4] - yg[i][3] * y[5]);

		ypg[i][3] = p->gamma * yg[i][3] - c2 * yg[i][5]
		- c5 * (p->u_re * yg[i][2] - p->u_im * yg[i][3]) - c4 * (p->u_re * yg[i][4] - p->u_im * yg[i][5])
		- c3 * (y[2] * yg[i][5] + y[3] * yg[i][4] + yg[i][2] * y[5] + yg[i][3] * y[4]);
	}

	if (p->is_cur_idx && p->want[0])
	{
		c6 = (p->rho[1][2] * Gamma_f - p->rho[0][2] * Gamma_d);
		c7 = p->alpha * p->rho[2][3];

		ypg[0][2] += c6 * (u_re2 - u_im2 - p->u_im) + c7 * (p->u_re * y[3] + p->u_im * y[2]);
		ypg[0][3] += c6 * (p->u_re + u_re_u_im) - c7 * (p->u_re * y[2] - p->u_im * y[3]);
	}
	if (p->is_cur_idx && p->want[1])
	{
		c6 = 2.0 * p->rho[0][3] * Gamma_d;
		c7 = 2.0 * (p->rho[1][3] * Gamma_f - p->rho[0][3] * Gamma_d);
		c8 = p->rho[2][3] * p->beta;

		ypg[1][2] += -c6 * y[4] + c8 * (p->u_re * y[3] + p->u_im * y[2]) + c7 * (p->u_re * y[5] + p->u_im * y[4])
			- 4.0 * p->alpha * (y[2] * y[4] - y[3] * y[5]);

		ypg[1][3] += -c6 * y[5] - c8 * (p->u_re * y[2] - p->u_im * y[3]) - c7 * (p->u_re * y[4] - p->u_im * y[5])
			- 4.0 * p->alpha * (y[2] * y[5] + y[3] * y[4]);
	}
	if (p->is_cur_idx && p->want[2])
	{
		ypg[2][2] += y[2] - 2.0 * y[4];
		ypg[2][3] += y[3] - 2.0 * y[5];
	}
	if (p->is_cur_idx && p->want[3])
	{
		c6 = p->alpha * p->beta;

		ypg[3][2] += c6 * (p->u_re * y[3] + p->u_im * y[2]);
		ypg[3][3] += -c6 * (p->u_re * y[2] - p->u_im * y[3]);
	}

	// dy5/dt, dy6/dt
	c1 = p->beta * p->beta / 2.0;
	c2 = 2.0 * p->gamma;
	c3 *= 2.0;
	c5 *= 2.0;

	yp[4] = c1 * (u_re2 - u_im2 - p->u_im) + c2 * y[4] + c5 * (p->u_re * y[5] + p->u_im * y[4])
		- c3 * (y[4] * y[4] - y[5] * y[5]) / 2.0;

	yp[5] = c1 * (p->u_re + u_re_u_im) + c2 * y[5] - c5 * (p->u_re * y[4] - p->u_im * y[5])
		- c3 * y[4] * y[5];

	for (i=0; i < 4; i++) if (p->want[i])
	{
		ypg[i][4] = c2 * yg[i][4] + c5 * (p->u_re * yg[i][5] + p->u_im * yg[i][4])
		- c3 * (y[4] * yg[i][4] - y[5] * yg[i][5]);

		ypg[i][5] = c2 * yg[i][5] - c5 * (p->u_re * yg[i][4] - p->u_im * yg[i][5])
		- c3 * (y[4] * yg[i][5] + yg[i][4] * y[5]);
	}

	if (p->is_cur_idx && p->want[0])
	{
		c6 = 2.0 * p->alpha * p->rho[2][3];

		ypg[0][4] += p->beta * (u_re2 - u_im2 - p->u_im) + c6 * (p->u_re * y[5] + p->u_im * y[4]);
		ypg[0][5] += p->beta * (p->u_re + u_re_u_im) - c6 * (p->u_re * y[4] - p->u_im * y[5]);
	}
	if (p->is_cur_idx && p->want[1])
	{
		c6 = 2.0 * p->rho[2][3] * p->beta;

		ypg[1][4] += c6 * (p->u_re * y[5] + p->u_im * y[4]) - 4.0 * p->alpha * (y[4] * y[4] - y[5] * y[5]);
		ypg[1][5] += -c6 * (p->u_re * y[4] - p->u_im * y[5]) - 8.0 * p->alpha * y[4] * y[5];
	}
	if (p->is_cur_idx && p->want[2])
	{
		ypg[2][4] += 2.0 * y[4];
		ypg[2][5] += 2.0 * y[5];
	}
	if (p->is_cur_idx && p->want[3])
	{
		c6 = 2.0 * p->alpha * p->beta;

		ypg[3][4] += c6 * (p->u_re * y[5] + p->u_im * y[4]);
		ypg[3][5] += -c6 * (p->u_re * y[4] - p->u_im * y[5]);
	}
}

static Err FxSVCalcPhi(SFxSVComm_InvFT *q, double u_re, double u_im, double *h_re, double *h_im)
{
	Err				err = NULL;
	int				idx_d, idx_f, idx_x, i, j, k, neq;
	SFxSVUndFull	*o = q->o;
	SFxSVComm_ODE	comm_RK;
	Nag_User		comm_Nag;
	NagError		fail;
	Nag_ODE_RK		opt;
	double			y[30], yp[30], ymax[30], RKthres[30], tgot, norm;
	double			*yg[4], *ypg[4], yg_re, yg_im;
	const double	RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7;		// adjustable

	memset(&fail, 0, sizeof(NagError));
	memset(&opt, 0, sizeof(Nag_ODE_RK));

	for (k=0, neq = 6; k < 4; k++) if (q->want[k])
	{
		yg[k] = y + neq;
		ypg[k] = yp + neq;
		neq += 6;
	}

	// Preinitialize comm structure for NAG

	comm_Nag.p = &comm_RK;
	comm_RK.T_pay = q->T_pay;
	comm_RK.u_re = u_re;
	comm_RK.u_im = u_im;
	comm_RK.lam_d = o->lam_d;
	comm_RK.lam_f = o->lam_f;
	memcpy(comm_RK.want, q->want, 4 * sizeof(int));

	memset(y, 0, neq * sizeof(double));	// Final y values are all zeros.
	for (j=0; j < neq; j++) RKthres[j] = thres_min;

	idx_d = o->nsig_d-1;
	idx_f = o->nsig_f-1;
	idx_x = o->pmdl->ntimes-1;

	// Proceed backwards integrating the system of ODE using Runge-Kutta

	for (i=o->ntimes-2; i >= 0; i--)
	{
		while (idx_d > 0 && o->sigtms_d[idx_d-1] > o->times[i] + 1e-5) idx_d--;
		while (idx_f > 0 && o->sigtms_f[idx_f-1] > o->times[i] + 1e-5) idx_f--;
		while (idx_x > 0 && o->pmdl->times[idx_x-1] > o->times[i] + 1e-5) idx_x--;

		// Fill in local constant coefficients in comm_Nag
		comm_RK.sig_d = o->sig_d[idx_d];
		comm_RK.sig_f = o->sig_f[idx_f];
		comm_RK.beta = o->pmdl->beta[idx_x];
		comm_RK.alpha = o->pmdl->alpha[idx_x];
		comm_RK.gamma = o->pmdl->gamma[idx_x];
		comm_RK.rho = o->pmdl->rho[idx_x];
		comm_RK.is_cur_idx = (idx_x >= q->idx_from && idx_x <= q->idx_to);

		// Calculate solution at time i

		nag_ode_ivp_rk_setup(neq, o->times[i+1], y, o->times[i], RKtol, RKthres,
			Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0, &opt, &fail);
		if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }

		nag_ode_ivp_rk_range(neq, EvalDerivatives, o->times[i], &tgot, y, yp, ymax, &opt, &comm_Nag, &fail);
		if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }

		nag_ode_ivp_rk_free(&opt);

		for (j=0; j < neq; j++)
		{
			RKthres[j] = thres_coef * ymax[j];
			if (RKthres[j] < thres_min) RKthres[j] = thres_min;
		}
	}
	norm = exp(y[0] + y[2] + y[4]);
	h_re[0] = norm * cos(y[1] + y[3] + y[5]);
	h_im[0] = norm * sin(y[1] + y[3] + y[5]);

	for (k=0; k < 4; k++) if (q->want[k])
	{
		yg_re = yg[k][0] + yg[k][2] + yg[k][4];
		yg_im = yg[k][1] + yg[k][3] + yg[k][5];
		h_re[k+1] = h_re[0] * yg_re - h_im[0] * yg_im;
		h_im[k+1] = h_re[0] * yg_im + h_im[0] * yg_re;
	}

FREE_RETURN:
	nag_ode_ivp_rk_free(&opt);

	return err;
}

// Complex Fourier transform of the density of log(FFXt/FFX0) h(u):

Err FxSVDensityFT(char *fxundname, long fix_date, long pay_date,
				  double u_re, double u_im, double *h_re, double *h_im, int *want, int idx)
{
	Err				err = NULL;
	SrtUndPtr		und;
	SFxSVUndFull	und_full;
	SFxSVComm_InvFT	comm_InvFT;
//	double			h_re_test, h_im_test;

	memset(&und_full, 0, sizeof(SFxSVUndFull));
	memset(&comm_InvFT, 0, sizeof(SFxSVComm_InvFT));

	und = lookup_und(fxundname);
	if (!und) return serror("Cannot find underlying %s", fxundname);

	err = FxSVFullFromUnd(&und_full, und, fix_date);
	if (err) goto FREE_RETURN;

	comm_InvFT.o = &und_full;
	comm_InvFT.T_pay = (pay_date - und_full.pmdl->today) * YEARS_IN_DAY;
	memcpy(comm_InvFT.want, want, 4 * sizeof(int));
	comm_InvFT.idx_from = comm_InvFT.idx_to = idx;

//	err = FxSVCalcPhi(&comm_InvFT, u_re, u_im, h_re, h_im);
//	if (err) goto FREE_RETURN;

	err = AffineCalcPhi(&comm_InvFT, u_re, u_im, h_re, h_im);
	if (err) goto FREE_RETURN;


FREE_RETURN:
	FxSVFullFree(&und_full);

	return err;
}

static double TailIntegrand(double v, void *comm)
{
	SFxSVComm_InvFT		*q = (SFxSVComm_InvFT *) comm;
	double				dd1, dd2, vk = v * q->k[q->want_k];

	dd1 = 1.0 + v * v;
	dd2 = dd1 * v;
	return cos(vk) / dd1 + sin(vk) / dd2;
}

static double IntegrateTail(SFxSVComm_InvFT *q)
{
	const double	tiny = 1.0e-307, halfpi = 1.5707963267949, tol = 1e-6;		// adjustable
	const int		lenaw = 8000;
    double			aw[8000], z, int_err, freq = fabs(q->k[q->want_k]);

	if (freq < 1e-16) return halfpi - atan(q->v_max);

	intdeoini(lenaw, tiny, tol, aw);
	intdeo(TailIntegrand, q->v_max, freq, aw, &z, &int_err, q);
	if (int_err < 0.0) { smessage("Tail integration failed");  return log(-1.0); }

	return z;
}

Err FxSVInitCache(SFxSVCache *cache, int maxpts, int nk)
{
	cache->maxpts = maxpts;
	cache->v = (double *) calloc(maxpts, sizeof(double));
	cache->fn = f3tensor(0, maxpts-1, 0, nk-1, 0, 4);
	if (!cache->v || !cache->fn) return serror("Memory failure");
	memset(cache->v, 0, maxpts * sizeof(double));
	return NULL;
}

Err FxSVFreeCache(SFxSVCache *cache, int nk)
{
	if (cache->next) FxSVFreeCache(cache->next, nk);
	free(cache->next);
	free(cache->v);
	if (cache->fn) free_f3tensor(cache->fn, 0, cache->maxpts-1, 0, nk-1, 0, 4);
	return NULL;
}

// Calculation of cos(vk)Re[Zeta(v)] + sin(vk)Im[Zeta(v)]

static double NAG_CALL Z_func(double v, Nag_User *comm)
{
	Err				err = NULL;
	SFxSVComm_InvFT	*q = (SFxSVComm_InvFT *) comm->p;
	double			phi_re[5], phi_im[5], zeta_re[5], zeta_im[5];
	double			dd1, dd2;
	int				i, j;

	// Retrieve integrand value from cache or calculate it:

	if (fabs(v - q->c_head->v[0]) < 1e-16)
	{
		q->c_cur = q->c_head;
		q->c_cur->cur_pt = 0;	// reset the counter if restarted from 0
	}
	else if (fabs(v - q->c_cur->v[q->c_cur->cur_pt]) > 1e-16)		// not yet calculated or no match
	{
		if (q->c_cur->v[q->c_cur->cur_pt] != 0.0)					// no match (branch point)
		{
			while (q->c_cur->next && fabs(v - q->c_cur->next->v[0]) > 1e-16)
				q->c_cur = q->c_cur->next;							// search in existing branches

			if (!q->c_cur->next)			// if not found - create a new branch
			{
				q->c_cur->next = (SFxSVCache *) malloc(sizeof(SFxSVCache));
				if (!q->c_cur->next) { smessage("Memory failure");  return log(-1.0); }
				memset(q->c_cur->next, 0, sizeof(SFxSVCache));
				err = FxSVInitCache(q->c_cur->next, q->c_cur->maxpts, q->nk);
				if (err) { smessage(err);  return log(-1.0); }
			}
			q->c_cur = q->c_cur->next;
			q->c_cur->cur_pt = 0;
		}

		if (q->c_cur->v[q->c_cur->cur_pt] == 0.0)	// point not yet calculated -> calculate it
		{
//			err = FxSVCalcPhi(q, v, -1.0, phi_re, phi_im);
			err = AffineCalcPhi(q, v, -1.0, phi_re, phi_im);

			if (err) { smessage(err);  return log(-1.0); }		// return NaN if error
			q->count++;		// only increase the function call counter if not retrieving from cache
			q->c_cur->v[q->c_cur->cur_pt] = v;

			// calculate zeta for all functions:
			dd1 = 1.0 + v * v;
			dd2 = dd1 * v;

			for (j=0; j < 5; j++) if (j==0 || q->want[j-1])
			{
				zeta_re[j] = ((j==0) - phi_re[j]) / dd1 + phi_im[j] / dd2;
				zeta_im[j] = -phi_im[j] / dd1 + ((j==0) - phi_re[j]) / dd2;
			}

			// calculate integrand for all strikes and all functions:
			for (i=0; i < q->nk; i++)
			{
				dd1 = cos(v * q->k[i]);
				dd2 = sin(v * q->k[i]);

				for (j=0; j < 5; j++) if (j==0 || q->want[j-1])
					q->c_cur->fn[q->c_cur->cur_pt][i][j] = dd1 * zeta_re[j] + dd2 * zeta_im[j];
			}
		}
	}

	if (++q->c_cur->cur_pt >= q->c_cur->maxpts)
	{ smessage("maxpts exceeded in Z_func. Contact FIRST");  return log(-1.0); }

	return q->c_cur->fn[q->c_cur->cur_pt-1][q->want_k][q->want_fn];
}

static Err FxSVDoIntegration(SFxSVComm_InvFT *q, int idx_k, int idx_fn, double *res)
{
	Err				err = NULL;
	const double	epsabs = 1e-7, epsrel = 1e-4;		// adjustable
	double			int_err;
	Nag_User		comm_Nag;
	NagError		fail;
	Nag_QuadProgress qp;

	memset(&fail, 0, sizeof(NagError));
	memset(&qp, 0, sizeof(Nag_QuadProgress));
	comm_Nag.p = q;
	q->want_k = idx_k;
	q->want_fn = idx_fn;

	nag_1d_quad_gen_1( Z_func, 0.0, q->v_max, epsabs, epsrel, 200,
		res, &int_err, &qp, &comm_Nag, &fail );
		if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }

	if (idx_fn == 0) *res += IntegrateTail(q);

FREE_RETURN:
	NAG_FREE(qp.sub_int_beg_pts);
	NAG_FREE(qp.sub_int_end_pts);
	NAG_FREE(qp.sub_int_result);
	NAG_FREE(qp.sub_int_error);
	return err;
}

// Function called by NAG while calculating moments of Xt:

static void NAG_CALL EvalDerivatives2(Integer neq, double t, double y[], double yp[], Nag_User *comm)
{
	SFxSVComm_ODE	*p = (SFxSVComm_ODE *) comm->p;
	double			Gamma_d, Gamma_f;
	double			c1, c2, c3, c4, c5;

	Gamma_d = -p->sig_d * ( 1.0 - exp(-p->lam_d*(p->T_pay - t)) ) / p->lam_d;
	Gamma_f = -p->sig_f * ( 1.0 - exp(-p->lam_f*(p->T_pay - t)) ) / p->lam_f;

	// dy8/dt, dy13/dt
	c1 = Gamma_d * Gamma_d / 2.0 + Gamma_f * Gamma_f / 2.0 - p->rho[0][1] * Gamma_d * Gamma_f;
	c2 = p->rho[0][3] * p->alpha * Gamma_d + p->gamma;
	c3 = p->alpha * p->alpha;
	c4 = 2.0 * p->alpha * (p->rho[1][3] * Gamma_f - p->rho[0][3] * Gamma_d);

	yp[0] = c1 - c2 * y[1] - c3 * y[2];
	yp[3] = 2.0 * c1 - c2 * y[4] - c3 * y[5] + c4 * y[1] + c3 * y[1]*y[1];

	// dy10/dt, dy15/dt
	c1 = p->beta * (p->rho[1][2] * Gamma_f - p->rho[0][2] * Gamma_d);
	c2 *= 2.0;
	c3 *= 4.0;
	c4 *= 2.0;
	c5 = 2.0 * p->alpha * p->rho[2][3] * p->beta;

	yp[1] = c1 + p->gamma * y[1] - c2 * y[2];
	yp[4] = 2.0 * c1 + p->gamma * y[4] - c2 * y[5] + c5 * y[1] + c4 * y[2] + c3 * y[1]*y[2];

	// dy12/dt, dy17/dt
	c1 = p->beta * p->beta / 2.0;
	c2 = 2.0 * p->gamma;
	c5 *= 2.0;

	yp[2] = c1 + c2 * y[2];
	yp[5] = 2.0 * c1 + c2 * y[5] + c5 * y[2] + c3 * y[2]*y[2];
}

// Calculate mean and std of log(FFX_T/FFX_0)

static Err FxSVCalcMoments(SFxSVComm_InvFT *q, double *mean, double *std)
{
	Err				err = NULL;
	int				idx_d, idx_f, idx_x, i, j;
	SFxSVUndFull	*o = q->o;
	SFxSVComm_ODE	comm_RK;
	Nag_User		comm_Nag;
	NagError		fail;
	Nag_ODE_RK		opt;
	double			y[6], yp[6], ymax[6], RKthres[6], tgot;
	const double	RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7;		// adjustable

	memset(&fail, 0, sizeof(NagError));
	memset(&opt, 0, sizeof(Nag_ODE_RK));

	// Preinitialize comm structure for NAG

	comm_Nag.p = &comm_RK;
	comm_RK.T_pay = q->T_pay;
	comm_RK.u_re = comm_RK.u_im = 0.0;
	comm_RK.lam_d = o->lam_d;
	comm_RK.lam_f = o->lam_f;

	memset(y, 0, 6 * sizeof(double));	// Final y values are all zeros.
	for (j=0; j < 6; j++) RKthres[j] = thres_min;

	idx_d = o->nsig_d-1;
	idx_f = o->nsig_f-1;
	idx_x = o->pmdl->ntimes-1;

	// Proceed backwards integrating the system of 6 ODE using Runge-Kutta

	for (i=o->ntimes-2; i >= 0; i--)
	{
		while (idx_d > 0 && o->sigtms_d[idx_d-1] > o->times[i] + 1e-5) idx_d--;
		while (idx_f > 0 && o->sigtms_f[idx_f-1] > o->times[i] + 1e-5) idx_f--;
		while (idx_x > 0 && o->pmdl->times[idx_x-1] > o->times[i] + 1e-5) idx_x--;

		// Fill in local constant coefficients in comm_Nag
		comm_RK.sig_d = o->sig_d[idx_d];
		comm_RK.sig_f = o->sig_f[idx_f];
		comm_RK.beta = o->pmdl->beta[idx_x];
		comm_RK.alpha = o->pmdl->alpha[idx_x];
		comm_RK.gamma = o->pmdl->gamma[idx_x];
		comm_RK.rho = o->pmdl->rho[idx_x];

		// Calculate solution at time i

		nag_ode_ivp_rk_setup(6, o->times[i+1], y, o->times[i], RKtol, RKthres,
			Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0, &opt, &fail);
		if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }

		nag_ode_ivp_rk_range(6, EvalDerivatives2, o->times[i], &tgot, y, yp, ymax, &opt, &comm_Nag, &fail);
		if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }

		nag_ode_ivp_rk_free(&opt);

		for (j=0; j < 6; j++)
		{
			RKthres[j] = thres_coef * ymax[j];
			if (RKthres[j] < thres_min) RKthres[j] = thres_min;
		}
	}
	*mean = y[0] + y[1] + y[2];
	*std = sqrt( -(y[3] + y[4] + y[5]) );

FREE_RETURN:
	nag_ode_ivp_rk_free(&opt);

	return err;
}

Err FxSVOptions(char *fxundname, long fix_date, long pay_date, int nK,
				double *K, char **rec_pay_str, int *want, int idx, double **res)
{
	Err				err = NULL;
	SrtUndPtr		und;
	SFxSVUndFull	und_full;
	SFxSVComm_InvFT	comm_InvFT;
	SFxSVCache		cache;
	const double	pi = 3.14159265358979, nstd = 20.0;		// adjustable
	const int		maxpts = 8000;
	SrtReceiverType	rec_pay;
	long			spot_date;
	char			*yc_dom, *yc_for;
	double			ffx0, df_dom, df_for, z, std;//, z_test, std_test;
	int				i, j;

	memset(&und_full, 0, sizeof(SFxSVUndFull));
	memset(&comm_InvFT, 0, sizeof(SFxSVComm_InvFT));
	memset(&cache, 0, sizeof(SFxSVCache));

	und = lookup_und(fxundname);
	if (!und) return serror("Cannot find underlying %s", fxundname);

	err = FxSVFullFromUnd(&und_full, und, fix_date);
	if (err) goto FREE_RETURN;

	comm_InvFT.o = &und_full;
	comm_InvFT.T_pay = (pay_date - und_full.pmdl->today) * YEARS_IN_DAY;

	memset(want, 0, 4 * sizeof(int));		// affine
	
	memcpy(comm_InvFT.want, want, 4 * sizeof(int));
	comm_InvFT.idx_from = comm_InvFT.idx_to = idx;

	// Calculate upper limit of integration
//	err = FxSVCalcMoments(&comm_InvFT, &z, &std);
//	if (err) goto FREE_RETURN;

	err = AffineCalcMoments(&comm_InvFT, &z, &std);
	if (err) goto FREE_RETURN;

	comm_InvFT.v_max = nstd / std;

	comm_InvFT.nk = nK;
	comm_InvFT.k = (double *) calloc(nK, sizeof(double));
	if (!comm_InvFT.k) { err = serror("Memory failure");  goto FREE_RETURN; }

	err = FxSVInitCache(&cache, maxpts, nK);
	if (err) goto FREE_RETURN;

	comm_InvFT.c_head = comm_InvFT.c_cur = &cache;

	spot_date = add_unit(und_full.pmdl->today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc_dom = get_discname_from_underlying( lookup_und(und_full.pmdl->dom_und) );
	yc_for = get_discname_from_underlying( lookup_und(und_full.pmdl->for_und) );

	df_dom = 1.0; //swp_f_df(spot_date, pay_date, yc_dom);		affine
	df_for = 1.0; //swp_f_df(spot_date, pay_date, yc_for);		affine

	ffx0 = und_full.pmdl->spot * df_for / df_dom;

	for (i=0; i < nK; i++) comm_InvFT.k[i] = log(K[i] / ffx0);

	// Integrate for all strikes and all functions:

	for (i=0; i < nK; i++)
	{
		err = interp_rec_pay(rec_pay_str[i], &rec_pay);
		if (err) goto FREE_RETURN;

		for (j=0; j < 5; j++) if (j==0 || want[j-1])
		{
			err = FxSVDoIntegration(&comm_InvFT, i, j, &z);
			if (err) goto FREE_RETURN;

			res[i][j] = z / pi * ffx0;
			if (j==0 && rec_pay == SRT_PUT && K[i] > ffx0) res[i][j] += K[i] - ffx0;
			if (j==0 && rec_pay == SRT_CALL && K[i] < ffx0) res[i][j] += ffx0 - K[i];

//			res[i][j] *= swp_f_df(und_full.pmdl->today, pay_date, yc_dom);		affine
		}
	}

FREE_RETURN:
	FxSVFullFree(&und_full);
	free(comm_InvFT.k);
	FxSVFreeCache(&cache, nK);

	return err;
}

static Err EvalF(SFxSVComm_Calib *r, int idx_from, int idx_to, double *f)
{
	Err				err = NULL;
	SFxSVComm_InvFT	comm_InvFT;
	SFxSVCache		cache;
	int				i;
	double			z, std, *mkt_val;
	const double	nstd = 20.0;			// adjustable
	const int		maxpts = 8000;

	memset(&comm_InvFT, 0, sizeof(SFxSVComm_InvFT));
	memset(&cache, 0, sizeof(SFxSVCache));

	comm_InvFT.o = r->o;
	comm_InvFT.T_pay = r->T_pay[idx_to];
	comm_InvFT.idx_from = idx_from;
	comm_InvFT.idx_to = idx_to;
	comm_InvFT.want[0] = (idx_from <= idx_to);		// calibrating beta

	// Calculate upper limit of integration:
	err = FxSVCalcMoments(&comm_InvFT, &z, &std);
	if (err) goto FREE_RETURN;

	comm_InvFT.v_max = nstd / std;

	if (idx_from <= idx_to)			// calibrating beta
	{
		comm_InvFT.nk = 1;
		comm_InvFT.k = r->k[idx_to];
		mkt_val = r->mkt_val[idx_to];
	}
	else							// calibrating smile params
	{
		comm_InvFT.nk = r->nk[idx_to] - 1;
		comm_InvFT.k = r->k[idx_to] + 1;
		mkt_val = r->mkt_val[idx_to] + 1;
	}

	// Initialize integrands cache:
	err = FxSVInitCache(&cache, maxpts, comm_InvFT.nk);
	if (err) goto FREE_RETURN;

	comm_InvFT.c_head = comm_InvFT.c_cur = &cache;

	// Evaluate the residuals at the given point:
	for (i=0; i < comm_InvFT.nk; i++)
	{
		err = FxSVDoIntegration(&comm_InvFT, i, 0, &z);
		if (err) goto FREE_RETURN;

		f[i] = z - mkt_val[i];
	}
	if (idx_from <= idx_to)			// calibrating beta
	{
		err = FxSVDoIntegration(&comm_InvFT, 0, 1, &f[1]);
		if (err) goto FREE_RETURN;
	}

FREE_RETURN:
	FxSVFreeCache(&cache, comm_InvFT.nk);
	return err;
}

static Err CalibrateBeta(SFxSVComm_Calib *r)
{
	Err				err = NULL;
	clock_t			t1, t2;
	SFxSVUndDesc	*pmdl = r->o->pmdl;
	int				i, j, k;
	double			x, f[2];
	const double	tol = 1e-4;		// adjustable

	// Calibrate beta from r->idx_from to r->idx_to

	r->o->ntimes = r->ntimes_last;
	for (i=r->idx_from; i <= r->idx_to; i = j+1)
	{
		for (j=i; !r->calib[j][0]; j++);	// Find next beta calibration date

		x = pmdl->beta[i > 0 ? i-1 : 0];	// First guess
		for (k=i; k <= j; k++) pmdl->beta[k] = x;

		// Set last time in r->o to pmdl->times[j]:
		while ( fabs(r->o->times[r->o->ntimes-1] - pmdl->times[j]) > 1e-5 ) r->o->ntimes++;
		t1 = clock();

		// Launch Newton:
		do
		{
			err = EvalF(r, i, j, f);
			if (err) return err;

			if ( fabs(f[1]) > 1e-16 ) pmdl->beta[i] -= f[0] / f[1];
			if ( fabs(f[1]) <= 1e-16 || pmdl->beta[i] <= r->bl[0] || pmdl->beta[i] >= r->bu[0] )
			{
				smessage("Volatility calibration failed at exercise date %d", j+1);
				for (k=i; k <= j; k++) pmdl->beta[k] = x;
				break;
			}
			for (k=i+1; k <= j; k++) pmdl->beta[k] = pmdl->beta[i];

		} while ( fabs(f[0]) > tol );

		t2 = clock();
		smessage ("Volatility calibration at step %d, time in sec: %.2f", j+1, (double) (t2 - t1) / CLOCKS_PER_SEC);
	}
	return NULL;
}

void NAG_CALL LsqFun(Integer m, Integer n, double x[], double fvec[], Nag_Comm *comm)
{
	Err				err = NULL;
	SFxSVComm_Calib *r = (SFxSVComm_Calib *) comm->p;
	SFxSVUndDesc	*pmdl = r->o->pmdl;
	int				i;
	double			xx[4];

	smessage("Smile parameters calibration at step %d, pass %d...", r->idx_to+1, comm->nf);

	for (i=1, n=0; i < 4; i++) if (r->calib[r->idx_to][i])		// check bounds
	{
		xx[i] = x[n++];
		if (xx[i] <= r->bl[i]) comm->flag = -2 * i;
		if (xx[i] >= r->bu[i]) comm->flag = -2 * i + 1;
		if (comm->flag < 0) return;
	}

	for (i=r->idx_from; i <= r->idx_to; i++)	// modify the underlying according to x
	{
		if (r->calib[r->idx_to][1]) pmdl->alpha[i] = xx[1];
		if (r->calib[r->idx_to][2]) pmdl->gamma[i] = xx[2];
		if (r->calib[r->idx_to][3]) pmdl->rho[i][2][3] = pmdl->rho[i][3][2] = xx[3];
	}

	err = CalibrateBeta(r);
	if (err) { smessage(err);  comm->flag = -99;  return; }

	// r->o->ntimes is set correctly at r->idx_to after CalibrateBeta

	err = EvalF(r, r->idx_to+1, r->idx_to, fvec);
	if (err) { smessage(err);  comm->flag = -199;  return; }
}

static Err CalcCorrBounds(double **rho, double *bl, double *bu, double tol)
{
	double			a[16], r[4];
	double			b[2][2];
	double			c;
	NagError		fail;
	int				i, j, k;

	memset(&fail, 0, sizeof(NagError));
	for (i=0; i < 4; i++) for (j=0; j <= i; j++) a[4*i + j] = rho[i][j];

	nag_real_symm_eigenvalues(4, a, 4, r, &fail);
	if (fail.code != NE_NOERROR) return serror(fail.message);
	if (r[0] <= 0.0) return serror("Input correlation matrix is not positive definite");

	b[0][0] = -1.0;			b[0][1] = rho[2][3];
	b[1][0] = rho[2][3];	b[1][1] = 1.0;
	
	for (k=0; k < 2; k++)
	{
		for (i=0; i < 4; i++) for (j=0; j <= i; j++) a[4*i + j] = rho[i][j];
		a[14] = b[k][k];

		nag_real_symm_eigenvalues(4, a, 4, r, &fail);
		if (fail.code != NE_NOERROR) return serror(fail.message);
		if (r[0] > 0) b[k][!k] = b[k][k];

		while (b[k][1] - b[k][0] > tol)
		{
			c = (b[k][0] + b[k][1]) / 2.0;
			for (i=0; i < 4; i++) for (j=0; j <= i; j++) a[4*i + j] = rho[i][j];
			a[14] = c;

			nag_real_symm_eigenvalues(4, a, 4, r, &fail);
			if (fail.code != NE_NOERROR) return serror(fail.message);

			b[k][(r[0] > 0) ^ k] = c;
		}
	}
	*bl = b[0][1];
	*bu = b[1][0];
	return NULL;
}

#define MAXNK 50

Err FxSVCalibrate(SrtUndPtr			und,		// Starting points and params not calibrated must be initialized
				  int				**calib,	// Flags calibrate { beta, alpha, gamma, rho } (per ex date)
				  int				*nK,		// Number of strikes per exercise date
				  double			**K,		// Strikes per exercise date
				  double			**vol)		// Vols per exercise date per strike
{
	Err				err = NULL;
	clock_t			t1, t2;
	SFxSVUndDesc	*pmdl = (SFxSVUndDesc *) und->spec_desc;
	SFxSVUndFull	und_full;
	SFxSVComm_Calib	comm_Calib;
	Nag_Comm		comm;
	NagError		fail;
	Nag_E04_Opt		options;
	int				i, j, m, n;
	double			rhobl, rhobu;
	double			x[3], fvec[MAXNK], fjac[MAXNK][3], x_save[3], fsumsq;
	long			spot_date, pay_date;
	char			*yc_dom, *yc_for;
	double			ffx0, df_dom, df_for;
	const double	pi = 3.14159265358979, tol = 1e-4;		// adjustable

	memset(&und_full, 0, sizeof(SFxSVUndFull));
	memset(&comm_Calib, 0, sizeof(SFxSVComm_Calib));
	memset(&fail, 0, sizeof(NagError));
	memset(&comm, 0, sizeof(Nag_Comm));
	memset(&options, 0, sizeof(Nag_E04_Opt));

	nag_opt_init(&options);
	options.list = 0;
	options.print_level = Nag_NoPrint;
//	options.deriv_check = 0;
	options.optim_tol = tol;

	err = FxSVFullFromUnd(&und_full, und, pmdl->dates[pmdl->ntimes-1]);
	if (err) goto FREE_RETURN;

	comm_Calib.k = dmatrix(0, pmdl->ntimes-1, 0, MAXNK-1);
	comm_Calib.mkt_val = dmatrix(0, pmdl->ntimes-1, 0, MAXNK-1);
	comm_Calib.T_pay = (double *) calloc(pmdl->ntimes, sizeof(double));

	if (!comm_Calib.k || !comm_Calib.mkt_val || !comm_Calib.T_pay)
	{ err = serror("Memory failure");  goto FREE_RETURN; }

	spot_date = add_unit(pmdl->today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	yc_dom = get_discname_from_underlying( lookup_und(pmdl->dom_und) );
	yc_for = get_discname_from_underlying( lookup_und(pmdl->for_und) );

	comm.p = &comm_Calib;
	comm_Calib.o = &und_full;
	comm_Calib.calib = calib;
	comm_Calib.nk = nK;

	// Initialize strikes and mkt prices:
	for (i=0; i < pmdl->ntimes; i++)
	{
		pay_date = add_unit(pmdl->dates[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);
		df_dom = swp_f_df(spot_date, pay_date, yc_dom);
		df_for = swp_f_df(spot_date, pay_date, yc_for);
		ffx0 = pmdl->spot * df_for / df_dom;

		if (nK[i] > MAXNK) { err = serror("Maximum %d strikes allowed per ex date", MAXNK); goto FREE_RETURN; }

		comm_Calib.T_pay[i] = (pay_date - pmdl->today) * YEARS_IN_DAY;

		for (j=0; j < nK[i]; j++)
		{
			comm_Calib.k[i][j] = log(K[i][j] / ffx0);
			comm_Calib.mkt_val[i][j] = srt_f_optblksch( ffx0, K[i][j], vol[i][j], pmdl->times[i], 1.0,
				(K[i][j] < ffx0 ? SRT_PUT : SRT_CALL), PREMIUM );

			comm_Calib.mkt_val[i][j] *= pi / ffx0;
		}
	}

	// Set constant bounds:
	comm_Calib.bl[0] = 0.0;
	comm_Calib.bu[0] = 1.0;
	comm_Calib.bl[1] = 0.0;
	comm_Calib.bu[1] = 20.0;
	comm_Calib.bl[2] = -1.0;
	comm_Calib.bu[2] = 10.0;

	und_full.ntimes = 1;

	for (i=0; i < pmdl->ntimes; i = j+1)
	{
		comm_Calib.ntimes_last = und_full.ntimes;

		// Search for the next smile calibration date:
		for (j=i; j < pmdl->ntimes && !calib[j][1] && !calib[j][2] && !calib[j][3]; j++);
		if (j == pmdl->ntimes)			// no more smile calibration dates found
		{								// just calibrate beta then
			for (j--; j >= i && !calib[j][0]; j--);
			if (j >= i)
			{
				comm_Calib.idx_from = i;
				comm_Calib.idx_to = j;

				err = CalibrateBeta(&comm_Calib);
				if (err) goto FREE_RETURN;
			}
			if (++j > 0) for (; j < pmdl->ntimes; j++) pmdl->beta[j] = pmdl->beta[j-1];
			break;						// nothing more to calibrate
		}
		else							// step j is a smile calibration date
		{
			if (!calib[j][0])
			{
				err = serror("Error: Smile calibration date %d is not a volatility calibration date", j+1);
				goto FREE_RETURN;
			}

			comm_Calib.idx_from = i;
			comm_Calib.idx_to = j;

			t1 = clock();

			// Initialize first guesses:
			n = 0;
			if (calib[j][1]) x[n++] = pmdl->alpha[i==0 ? 0 : i-1];
			if (calib[j][2]) x[n++] = pmdl->gamma[i==0 ? 0 : i-1];
			if (calib[j][3])
			{
				comm_Calib.bl[3] = -1.0;
				comm_Calib.bu[3] = 1.0;

				for (m=i; m <= j; m++)
				{
					err = CalcCorrBounds(pmdl->rho[m], &rhobl, &rhobu, 0.01);
					if (err) goto FREE_RETURN;

					if (rhobl > comm_Calib.bl[3]) comm_Calib.bl[3] = rhobl;
					if (rhobu < comm_Calib.bu[3]) comm_Calib.bu[3] = rhobu;
				}
				if (comm_Calib.bu[3] - comm_Calib.bl[3] < 0.01)
				{
					calib[j][3] = 0;
					smessage("Rho calibration not feasible at step %d", j+1);
				}
				else
				{
					x[n] = pmdl->rho[i==0 ? 0 : i-1][2][3];
					if (x[n] <= comm_Calib.bl[3]) x[n] = comm_Calib.bl[3] + 0.001;
					if (x[n] >= comm_Calib.bu[3]) x[n] = comm_Calib.bu[3] - 0.001;
					n++;
				}
			}

			// Save the starting point:
			memcpy(x_save, x, 3 * sizeof(double));

			if (n > 0)		// if there's smth to calibrate
			{
				if (nK[j]-1 < n) { err = serror("Number of options must be greater or equal to the number of "
					"calibrated parameters at exercise date %d", j+1);  goto FREE_RETURN; }

				// Launch non-linear least squares minimization:

				nag_opt_lsq_no_deriv(nK[j]-1, n, LsqFun, x, &fsumsq, fvec, &fjac[0][0], 3, &options, &comm, &fail);

				if (fail.code == NE_USER_STOP)
				{
					switch (fail.errnum)
					{
						case -1: smessage("Alpha has hit the upper bound at step %d", j+1);  break;
						case -2: smessage("Alpha has hit the lower bound at step %d", j+1);  break;
						case -3: smessage("Gamma has hit the upper bound at step %d", j+1);  break;
						case -4: smessage("Gamma has hit the lower bound at step %d", j+1);  break;
						case -5: smessage("Rho has hit the upper bound at step %d", j+1);  break;
						case -6: smessage("Rho has hit the lower bound at step %d", j+1);  break;
						default: err = serror("Calibration failed at date %d", j+1);  goto FREE_RETURN;
					}
					smessage("Smile calibration skipped at step %d", j+1);
					memcpy(x, x_save, 3 * sizeof(double));
				}
				else if (fail.code != NE_NOERROR) { err = serror(fail.message);  goto FREE_RETURN; }
			}

			// Modify the underlying with calibrated values:
			for (m=i; m <= j; m++)
			{
				n = 0;
				if (calib[j][1]) pmdl->alpha[m] = x[n++];
				if (calib[j][2]) pmdl->gamma[m] = x[n++];
				if (calib[j][3]) pmdl->rho[m][2][3] = pmdl->rho[m][3][2] = x[n++];
			}

			t2 = clock();
			smessage ("Smile calibration at step %d, time in sec: %.2f", j+1, (double) (t2 - t1) / CLOCKS_PER_SEC);

			if (n == 0 || fail.code == NE_USER_STOP)
			{
				calib[j][1] = calib[j][2] = calib[j][3] = 0;	// That means restart the calibration from i
				j = i-1;										// skipping the smile calibration at step j
				und_full.ntimes = comm_Calib.ntimes_last;
			}
		}			// if step j is a smile calibration date
	}			// for (i=0; i < pmdl->ntimes; i = j+1)

FREE_RETURN:
	FxSVFullFree(&und_full);
	if(comm_Calib.k) free_dmatrix(comm_Calib.k, 0, pmdl->ntimes-1, 0, MAXNK-1);
	if(comm_Calib.mkt_val) free_dmatrix(comm_Calib.mkt_val, 0, pmdl->ntimes-1, 0, MAXNK-1);
	free(comm_Calib.T_pay);
	nag_opt_free(&options, "", &fail);

	return err;
}
