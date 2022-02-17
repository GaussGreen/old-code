/*--------------------------------------------------------------
	FILE: CheyBetaMC.h
	PURPOSE: Cheyette beta Monte-Carlo
	AUTHOR: Dimitri Mayevski
	DATE: 26/02/2003
  --------------------------------------------------------------*/

#ifndef __CHEYBETAMC_H__
#define __CHEYBETAMC_H__

#include "Product.h"
#include "srt_h_cheybeta_new.h"

Err CheyBetaMC_Price(SProductDesc *g, int nt, long npaths, double cutcoef,
					 SCheyBeta *pmdl, int rgmode, double *pv, double *stddev);


#endif  /* #ifndef __CHEYBETAMC_H__ */