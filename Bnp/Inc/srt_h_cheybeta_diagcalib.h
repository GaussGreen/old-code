/*--------------------------------------------------------------
	FILE: srt_h_cheybeta_diagcalib.h
	PURPOSE: Cheyette beta diag calib (new)
	AUTHOR: Dimitri Mayevski
	DATE: 2/10/2002
  --------------------------------------------------------------*/

#ifndef __SRT_H_CHEYBETA_DIAGCALIB_H__
#define __SRT_H_CHEYBETA_DIAGCALIB_H__

/*	If end_date != 0 tenors are ignored, calibrating to diagonal swaptions using
	freq and basis. If end_date = 0 there are two cases: if tenors are present
	calibrate to corresponding swaptions using freq and basis, otherwise use refrate
	defaults and calibrate to corresponding caplets */

Err cheybeta_diagcalib_new(
	SCheyBeta		*pmdl,				/*	beta, lambda, ycname and today must be
											initialized inside. Sigmas are returned.
											To be freed by caller */
	char			*vc_name,			/*	Name of the market vol curve */
	double			vol_shift,
	int				shift_type,			/*	0:	Additive
											1:	Multiplicative */
	int				nex,
	long			*ex_dates,			/*	Exercise dates */
	long			end_date,			/*	End date for diagonal */
	char			**tenors,			/*	Instruments' tenors */
	char			*freq,
	char			*basis,
	char			*refrate,			/*	Name of the reference rate */
	double			*strikes,
	int				strike_type,		/*	0: ATM
											1: CASH
											2: SWAP
											3: STD */
	double			max_std,
	int				nt,					/*  Parameters of the		*/
	int				nx,					/*  CheyBeta forward PDE	*/
	int				nphi,				/*  pricing routine			*/
	double			cutcoef,
	Err				(*get_cash_vol)(				/*	Function to get cash vol (for LGM calib) */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power)
);

#endif  /* #ifndef __SRT_H_CHEYBETA_DIAGCALIB_H__ */