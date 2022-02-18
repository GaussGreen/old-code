/* ==========================================================================
   FILE_NAME:	opsabrgenericcalib.c

   PURPOSE:		MODIFIED IMPLEMENTATION OF SABR Calibration to SABR or Market Smile

   DATE:		03/03/04

   AUTHOR:		P.T.
   ========================================================================== */

#include "opfnctns.h"
#include <OPSABRCALIB.H>
#include <OPSABRGENERIC.H>
#include <OPSABRGENERICCALIB.H>

#include "math.h"

#include "utallhdr.h"

// Added for the NAG optimizer in opsabrgencalibpro

#include "nag.h"
#include "nag_stdlib.h"
#include "nage04.h"
#include "nagx02.h"

static Opt_params params;

/////////////////////////////////////////////////////////////////////////////////////////////
//
// op_sabrgen(), Returns BSlog(@K)
//
//	Returns vol(K) Given one ATM Vol and some a,b,c, params
//	input SigmaBeta ATMLog, ATmNorm
//	output Log or Norm
//	generic version, specify the kind of loc vol you wanna use
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err srt_f_optsabrgenvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           a,
    double           b,
    double           c,
    double           Rho,
    SrtDiffusionType input,
    SrtDiffusionType output,
    SrtDiffusionType volLoc,
    double*          vol)
{
    double res;
    double SigmaBeta;
    Err    err = NULL;

    FuncVolLocType vol_loc;
    GetLocVolFromDiffusionType(volLoc, &vol_loc);

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Fatal: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Fatal: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Fatal: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    // range check
    if (vol_parameters_check_from_SrtDiffusionType(a, b, c, 30, volLoc) == 0)
    {
        err = "Error: local volatility parameters out of range";
        goto FREE_RETURN;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta
        SigmaBeta =
            op_sabrgen_calib(Forward, Forward, Maturity, Sigma, Alpha, a, b, c, Rho, vol_loc);
    else  // Case input = Beta
        SigmaBeta = Sigma;

    if (output == SRT_BETAVOL)
        res = SigmaBeta;
    else
    {  // Convert to LOG
        res = op_sabrgen(Forward, Strike, Maturity, SigmaBeta, Alpha, a, b, c, Rho, vol_loc);
        if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
        {
            err = srt_f_optsarbvol(
                Forward, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
            if (err)
                goto FREE_RETURN;
        }
    }

    *vol = res;

FREE_RETURN:

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
// op_sabrgen_pricing(),
//
// Compute the SABRGEN Price from the parameters.
// Input:	ATMLOG vol
// Output:	Price
// Generic Version
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err op_sabrgen_pricing(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha,
    double a,
    double b,
    double c,
    double Rho,
    /*Loc Vol type*/
    SrtDiffusionType volLoc,
    /*	Vol input	*/
    double atmvol,
    /*	Premium	*/
    double* prem)
{
    double Betavol;
    // assign the proper loc vol
    FuncVolLocType vol_loc;
    GetLocVolFromDiffusionType(volLoc, &vol_loc);

    Betavol = op_sabrgen_calib(forward, forward, maturity, atmvol, Alpha, a, b, c, Rho, vol_loc);

    *prem = op_sabrgen(forward, strike, maturity, Betavol, Alpha, a, b, c, Rho, vol_loc);

    *prem = srt_f_optblksch(forward, strike, *prem, maturity, 1, SRT_CALL, SRT_PREMIUM);

    return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	op_sabrgen_model_sens()
//
//	Sensitivity to model parameters
//	All parameters checks are supposed to have been done before call
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err op_sabrgen_model_sens(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha,
    double a,
    double b,
    double c,
    double Rho,
    /*Loc Vol type*/
    SrtDiffusionType volLoc,
    /*	Vol input	*/
    double input_vol,
    /*	Options	*/
    /*	Size of bumps	*/
    double Alpha_shift,
    double a_shift,
    double b_shift,
    double c_shift,
    double Rho_shift,
    /*	Factor that multiplies results	*/
    double Alpha_fac,
    double a_fac,
    double b_fac,
    double c_fac,
    double Rho_fac,
    /*	Answers	*/
    double* Alpha_sens,
    double* a_sens,
    double* b_sens,
    double* c_sens,
    double* Rho_sens)
{
    double frozen_vol;
    double prem, Alpha_prem, a_prem, b_prem, c_prem, Rho_prem;
    Err    err;

    frozen_vol = input_vol;

    /*	Calculate premium before shift	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b,
            c,
            Rho,
            volLoc,
            frozen_vol,
            &prem))
    {
        return err;
    }

    /*	Calculate premium after shift of Alpha	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha + Alpha_shift,
            a,
            b,
            c,
            Rho,
            volLoc,
            frozen_vol,
            &Alpha_prem))
    {
        return err;
    }

    /*	Calculate premium after shift of a	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a + a_shift,
            b,
            c,
            Rho,
            volLoc,
            frozen_vol,
            &a_prem))
    {
        return err;
    }

    /*	Calculate premium after shift of b	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b + b_shift,
            c,
            Rho,
            volLoc,
            frozen_vol,
            &b_prem))
    {
        return err;
    }

    /*	Calculate premium after shift of c	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b,
            c + c_shift,
            Rho,
            volLoc,
            frozen_vol,
            &c_prem))
    {
        return err;
    }

    /*	Calculate premium after shift of Rho	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b,
            c,
            Rho + Rho_shift,
            volLoc,
            frozen_vol,
            &Rho_prem))
    {
        return err;
    }

    /*	Calculate sensitivites	*/
    *Alpha_sens = Alpha_fac * (Alpha_prem - prem) / Alpha_shift;
    *a_sens     = a_fac * (a_prem - prem) / a_shift;
    *b_sens     = b_fac * (b_prem - prem) / b_shift;
    *c_sens     = c_fac * (c_prem - prem) / c_shift;
    *Rho_sens   = Rho_fac * (Rho_prem - prem) / Rho_shift;

    return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	op_sabrgen_vega_volga()
//
//	Vega/Volga
//	All parameters checks are supposed to have been done before call
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err op_sabrgen_vega_volga(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha,
    double a,
    double b,
    double c,
    double Rho,
    /*Loc Vol type*/
    SrtDiffusionType volLoc,
    /*	Vol input	*/
    double input_vol,
    /*	Options	*/
    /*	Size of bumps	*/
    double vega_shift,
    double volga_shift,
    /*	Multiplicative or additive	*/
    /*	0: Additive, 1: Multiplicative	*/
    int vega_mult,
    int volga_mult,
    /*	Calculate vega or volga	*/
    /*	0: Vega, 1: Volga	*/
    int vega_volga,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double* sens)
{
    double bump_vol;
    double prem, prem2, vega1, vega2;
    double act_vega_shift, act_volga_shift;
    Err    err;

    bump_vol = input_vol;

    /*	Calculate premium before shift	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b,
            c,
            Rho,
            volLoc,
            bump_vol,
            &prem))
    {
        return err;
    }

    /*	Calculate actual shifts	*/

    if (vega_mult == 1)
    {
        act_vega_shift = vega_shift * bump_vol;
    }
    else
    {
        act_vega_shift = vega_shift;
    }

    if (volga_mult == 1)
    {
        act_volga_shift = volga_shift * bump_vol;
    }
    else
    {
        act_volga_shift = volga_shift;
    }

    /*	Calculate premium after shift	*/
    if (err = op_sabrgen_pricing(
            forward,
            strike,
            maturity,
            disc,
            call_put,
            Alpha,
            a,
            b,
            c,
            Rho,
            volLoc,
            bump_vol + act_vega_shift,
            &prem2))
    {
        return err;
    }

    /*	Calculate vega before shift	*/
    vega1 = (prem2 - prem) / vega_shift;
    *sens = fac * vega1;

    /*	If volga is to be calculated, calculate delta after shift	*/
    if (vega_volga == 1)
    {
        if (err = op_sabrgen_vega_volga(
                forward,
                strike,
                maturity,
                disc,
                call_put,
                Alpha,
                a,
                b,
                c,
                Rho,
                volLoc,
                bump_vol + act_volga_shift,
                vega_shift,
                volga_shift,
                vega_mult,
                volga_mult,
                0,
                1.0,
                &vega2))
        {
            return err;
        }

        /*	Calculate volga	*/
        *sens = fac * (vega2 - vega1) / volga_shift;
    }

    return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	SabrGenPriceGradient()
//
//	Used by the optimizer
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err SabrGenPriceGradient(
    double  strike,
    double  paramSABR[],
    double* price,
    double* gradient,
    int     nbr_paramSABR)  // useless ????????????????????
{
    Err              err = NULL;
    double           Alpha, a, b, c, Rho, forward, maturity, ATMvol, Betavol;
    SrtDiffusionType volLoc;
    FuncVolLocType   vol_loc;

    int out_range = 0, i;

    Alpha    = paramSABR[1];
    a        = paramSABR[2];
    b        = paramSABR[3];
    c        = paramSABR[4];
    Rho      = paramSABR[5];
    volLoc   = (SrtDiffusionType)(paramSABR[6]);  // need to cast it seeems to bo OK
    forward  = paramSABR[7];
    maturity = paramSABR[8];
    ATMvol   = paramSABR[9];

    // assign the proper loc vol
    GetLocVolFromDiffusionType(volLoc, &vol_loc);

    Betavol = op_sabrgen_calib(forward, strike, maturity, ATMvol, Alpha, a, b, c, Rho, vol_loc);

    if ((Alpha < 0) || (Rho < -1) || (Rho > 1))
        out_range = 1;

    if (out_range == 1)
    {
        *price = -1000000;
        for (i = 1; i <= 6; i++)
        {
            gradient[i] = 1000000;
        }
    }
    else
    {
        err = op_sabrgen_pricing(
            forward, strike, maturity, 1, SRT_CALL, Alpha, a, b, c, Rho, volLoc, ATMvol, price);

        err = op_sabrgen_model_sens(
            forward,
            strike,
            maturity,
            1,
            SRT_CALL,
            Alpha,
            a,
            b,
            c,
            Rho,
            volLoc,
            ATMvol,
            1e-4,  // think about hose guys
            1e-4,  // think about hose guys
            1e-4,  // think about hose guys
            1e-4,  // think about hose guys
            1e-4,  // think about hose guys
            1,
            1,
            1,
            1,
            1,
            &(gradient[1]),
            &(gradient[2]),
            &(gradient[3]),
            &(gradient[4]),
            &(gradient[5]));

        gradient[6] = 0; /* No need to calculate it Forward		*/
        gradient[7] = 0; /* No need to calculate it Maturity		*/
        gradient[8] = 0; /* No need to calculate it LocVoltype	*/

        err = op_sabrgen_vega_volga(
            forward,
            strike,
            maturity,
            1,
            SRT_CALL,
            Alpha,
            a,
            b,
            c,
            Rho,
            volLoc,
            ATMvol,
            1e-5,
            1e-4,
            1,
            1,
            0,
            1,
            &(gradient[9]));
    }

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	init_SABRGEN_parameters()
//
//	Used by the optimizer
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err init_SABRGEN_parameters(
    double  forward,
    double* strikes,
    double* market_vols,
    int     nbr_strikes,
    double  ATMVol,
    double* Alpha,
    int     freeze_Alpha,
    double* a,
    int     freeze_a,
    double* b,
    int     freeze_b,
    double* c,
    int     freeze_c,
    double* Rho,
    int     freeze_Rho)
{
    Err err = NULL;

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	opsabrgencalib()
//
//	Calibration routine calibrate parameters to market smile
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err opsabrgencalib(
    double           forward,
    double           maturity,
    int              nbr_strikes,
    double*          strikes,      // ACHTUNG: The Vector starts at 1
    double*          market_vols,  // ACHTUNG: The Vector starts at 1
    double*          ATMVol,
    double*          Alpha,
    int              freeze_Alpha,  // if 0, Alpha is not calibrated
    double*          a,
    int              freeze_a,  // If 0, a is not calibrated
    double*          b,
    int              freeze_b,  // If 0, b is not calibrated
    double*          c,
    int              freeze_c,  // If 0, c is not calibrated
    double*          Rho,
    int              freeze_Rho,  // if 0, Rho is not calibrated
    SrtDiffusionType VolLoc,      // Fake SABR parameter, flag for the loc vol ...
    double*          fitting_error)
{
    Err     err                = NULL;
    double *weights_of_strikes = NULL,
           *paramSABR          = NULL,  // Parameters of SabrPriceGradient
               *market_prices  = NULL;
    long *use_paramSABR        = NULL,  // 0 for frozen parameter, 1 for parameter to optimize
        nbr_paramSABR;
    int nbr_iter = 200;  // Maximum number of iterations in Levenberg-Marquardt
    int freeze_ATMVol;
    int i;

    // Memory allocation
    nbr_paramSABR      = 9;  // Alpha,a,b,c,Rho,LocVoltype,Forward,Maturity and ATMVol
    weights_of_strikes = dvector(1, nbr_strikes);
    paramSABR          = dvector(1, nbr_paramSABR);
    market_prices      = dvector(1, nbr_strikes);
    use_paramSABR      = lngvector(1, nbr_paramSABR);

    // Computes the ATM vol by interpolation
    err = compute_ATMVol(forward, strikes, market_vols, nbr_strikes, &freeze_ATMVol, ATMVol);

    // Computes a first guess for the BVMSABR parameters
    err = init_SABRGEN_parameters(
        forward,
        strikes,
        market_vols,
        nbr_strikes,
        *ATMVol,
        Alpha,
        freeze_Alpha,
        a,
        freeze_a,
        b,
        freeze_b,
        c,
        freeze_c,
        Rho,
        freeze_Rho);

    // Assignment of weights to strikes for fitting
    err = compute_weights(strikes, nbr_strikes, forward, *ATMVol, maturity, weights_of_strikes);

    // Initializes SABR parameters
    paramSABR[1] = *Alpha;
    paramSABR[2] = *a;
    paramSABR[3] = *b;
    paramSABR[4] = *c;
    paramSABR[5] = *Rho;
    paramSABR[6] = VolLoc;
    paramSABR[7] = forward;
    paramSABR[8] = maturity;
    paramSABR[9] = *ATMVol;

    // Chooses the parameters to be calibrated
    use_paramSABR[1] = freeze_Alpha;
    use_paramSABR[2] = freeze_a;
    use_paramSABR[3] = freeze_b;
    use_paramSABR[4] = freeze_c;
    use_paramSABR[5] = freeze_Rho;
    use_paramSABR[6] = 0;  // Not calibrating the loc vol type
    use_paramSABR[7] = 0;  // Not calibrating the forward
    use_paramSABR[8] = 0;  // Not calibrating the maturity
    use_paramSABR[9] = freeze_ATMVol;

    // Compute the market prices
    for (i = 1; i <= nbr_strikes; i++)
    {
        market_prices[i] =
            srt_f_optblksch(forward, strikes[i], market_vols[i], maturity, 1, SRT_CALL, PREMIUM);
    }

    // Call Levenberg-Marcquardt

    err = levenberg_marquardt_select(
        strikes,             // From [1] to [nbr_strikes]
        market_prices,       // From [1] to [nbr_strikes]
        weights_of_strikes,  // From [1] to [nbr_strikes]
        nbr_strikes,
        paramSABR,      // From [1] to [nparam]
        use_paramSABR,  // From [1] to [nparam]
        nbr_paramSABR,
        nbr_iter,
        SabrGenPriceGradient,
        fitting_error);

    // Fills the calibrated parameters
    *Alpha  = paramSABR[1];
    *a      = paramSABR[2];
    *b      = paramSABR[3];
    *c      = paramSABR[4];
    *Rho    = paramSABR[5];
    *ATMVol = paramSABR[9];

    // Free memory

    free_dvector(weights_of_strikes, 1, nbr_strikes);
    free_dvector(market_prices, 1, nbr_strikes);
    free_dvector(paramSABR, 1, nbr_paramSABR);
    free_lngvector(use_paramSABR, 1, nbr_paramSABR);

    return err;
}

static void NAG_CALL objfun(Integer n, double x[], double* f, double g[], Nag_Comm* comm);

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	opsabrgencalibpro()  //with NAG because you worth it
//
//	Calibration routine calibrate parameters to market smile
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err opsabrgencalibpro(
    double           forward,
    double           maturity,
    int              nbr_strikes,
    double*          strikes,
    double*          market_vols,
    double*          InitGuess,
    double*          lbounds,
    double*          ubounds,
    SrtDiffusionType VolLoc,
    double*          fitting_error)
{
    Err             err      = NULL;
    int             nb_param = 6;
    int             freeze_ATMVol;
    double          g[6];
    FuncVolLocType  vol_loc;
    double          objf, ATMVol;
    static NagError fail;
    Nag_BoundType   bound;  // NAG: says there are constraints
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer

    Opt_params params;
    e04xxc(&options);
    bound               = Nag_Bounds;
    fail.print          = FALSE;
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    options.optim_tol   = 0.00001;

    GetLocVolFromDiffusionType(VolLoc, &vol_loc);

    // Computes the ATM vol by interpolation
    err = compute_ATMVol(forward, strikes, market_vols, nbr_strikes, &freeze_ATMVol, &ATMVol);

    if (freeze_ATMVol == 0)
    {
        lbounds[0] = ATMVol;
        ubounds[0] = ATMVol;
    }

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = nbr_strikes;
    params.Vec_marketvols = market_vols;
    params.Vec_strikes    = strikes;
    params.vol_loc        = vol_loc;

    nagcomm.p = (void*)&params;

    e04jbc(
        nb_param, objfun, bound, lbounds, ubounds, InitGuess, &objf, g, &options, &nagcomm, &fail);

    *fitting_error = sqrt(objf / nbr_strikes);

    return err;
}

static void NAG_CALL objfun(Integer n, double x[], double* objf, double g[], Nag_Comm* comm)
{
    double         res;
    Opt_params*    pparams        = (Opt_params*)(comm->p);
    int            nb_strikes     = pparams->nb_strikes;
    double         forward        = pparams->forward;
    double         maturity       = pparams->maturity;
    double*        Vec_strikes    = pparams->Vec_strikes;
    double*        Vec_marketvols = pparams->Vec_marketvols;
    double         sigma, sigmabeta;
    int            i;
    FuncVolLocType vol_loc = pparams->vol_loc;

    res = 0.0;

    sigmabeta =
        op_sabrgen_calib(forward, forward, maturity, x[0], x[1], x[2], x[3], x[4], x[5], vol_loc);

    for (i = 0; i < nb_strikes; i++)
    {
        sigma = op_sabrgen(
            forward,
            Vec_strikes[i + 1],
            maturity,
            sigmabeta,
            x[1],
            x[2],
            x[3],
            x[4],
            x[5],
            vol_loc);

        res += (sigma - Vec_marketvols[i + 1]) * (sigma - Vec_marketvols[i + 1]);
    }
    *objf = res;
}
