/*--------------------------------------------------------------
	FILE: CheyBetaPDE.h
	PURPOSE: Cheyette beta PDE (forward & backward)
	AUTHOR: Dimitri Mayevski
	DATE: 22/11/2002
  --------------------------------------------------------------*/

#ifndef __CHEYBETAPDE_H__
#define __CHEYBETAPDE_H__

#include "Product.h"
#include "srt_h_cheybeta_new.h"

typedef struct _SCheyBetaPDEGrid
{
	int nx, nphi, ninst;
	double *phi, **ad, ***pv;	/* ad for forward, pv for backward */
} SCheyBetaPDEGrid;

typedef struct _SCheyBetaPDE {
	SProductDesc	*g;
	SCheyBeta		*pmdl;
	SCheyBetaPDEGrid grids[3];
	int				dir;	/* direction: grids[dir] -> grids[!dir] */
	int				nt, nx, nphi, last, savetime, ix0, isfwd;
	double			*times, *x, *df, *fr, *gamT, *phi_min, *phi_max, *phi_mid, *cumxvar;
	double			cutcoef;
	CNPDE_TEMP		*tmppde;
} SCheyBetaPDE;

Err CheyBetaPDE_Init(SCheyBetaPDE *pde,
		SProductDesc *g, int nt, int nx, int nphi, double cutcoef, int isfwd, SCheyBeta *pmdl);
Err CheyBetaPDE_Free(SCheyBetaPDE *pde);

/* Used for calibration */
Err CheyBetaPDE_ExpandFwd(SCheyBetaPDE *pde, double T);
Err CheyBetaPDE_ResetFwd(SCheyBetaPDE *pde, double T);
Err CheyBetaPDE_PriceSwaption(SCheyBetaPDE *pde, int iswp, double *result);

/* Used for pricing */
Err CheyBetaPDE_PriceBkwd(SCheyBetaPDE *pde, double *pv);

#endif  /* #ifndef __CHEYBETAPDE_H__ */