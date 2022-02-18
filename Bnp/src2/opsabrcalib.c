/* Calibrates the SABR parameters to a market smile */

#include "opfnctns.h"
#include <OPSABRCALIB.H>

#include "math.h"

#include "utallhdr.h"

/* Computes the value and the gradient of the SABR price at a given strike */

Err SabrPriceGradient(
    double strike, double paramSABR[], double* price, double* gradient, int nbr_paramSABR)
{
    Err    err = NULL;
    double alpha, beta, rho, forward, maturity, ATMvol;
    int    out_range = 0, i;

    alpha    = paramSABR[1];
    beta     = paramSABR[2];
    rho      = paramSABR[3];
    forward  = paramSABR[4];
    maturity = paramSABR[5];
    ATMvol   = paramSABR[6];

    if ((alpha < 0) || (beta < 0) || (beta > 1) || (rho < -1) || (rho > 1))
        out_range = 1;

    /*

    *vol=SabrLognormalVol(
                                                    ATMvol,
                                                    alpha,
                                                    beta,
                                                    rho,
                                                    forward,
                                                    strike,
                                                    maturity,
                                                    SABR_ATM_LOG); */

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
        err = op_sabr_pricing(
            forward, strike, maturity, 1, SRT_CALL, alpha, beta, rho, ATMvol, SABR_ATM_LOG, price);

        err = op_sabr_model_sens(
            forward,
            strike,
            maturity,
            1,
            SRT_CALL,
            alpha,
            beta,
            rho,
            ATMvol,
            SABR_ATM_LOG,
            1e-4,
            1e-4,
            1e-4,
            SABR_ATM_LOG,
            1,
            1,
            1,
            &(gradient[1]),
            &(gradient[2]),
            &(gradient[3]));

        /*
        gradient[1] = SabrDsigmaDalpha(
                                                        ATMvol,
                                                        alpha,
                                                        beta,
                                                        rho,
                                                        forward,
                                                        strike,
                                                        maturity,
                                                        SABR_ATM_LOG);







        gradient[2] = SabrDsigmaDbeta(
                                                        ATMvol,
                                                        alpha,
                                                        beta,
                                                        rho,
                                                        forward,
                                                        strike,
                                                        maturity,
                                                        SABR_ATM_LOG);


        gradient[3] = SabrDsigmaDrho(
                                                        ATMvol,
                                                        alpha,
                                                        beta,
                                                        rho,
                                                        forward,
                                                        strike,
                                                        maturity,
                                                        SABR_ATM_LOG); */

        gradient[4] = 0; /* No need to calculate it */
        gradient[5] = 0; /* No need to calculate it */
                         /*	gradient[6]= SabrDsigmaDvolatm(
                                                                                 ATMvol,
                                                                                 alpha,
                                                                                 beta,
                                                                                 rho,
                                                                                 forward,
                                                                                 strike,
                                                                                 maturity); */
        err = op_sabr_vega_volga(
            forward,
            strike,
            maturity,
            1,
            SRT_CALL,
            alpha,
            beta,
            rho,
            ATMvol,
            SABR_ATM_LOG,
            1e-5,
            1e-4,
            1,
            1,
            SABR_ATM_LOG,
            0,
            1,
            &(gradient[6]));
    }

    return err;
}

Err compute_weights(
    double* strikes,
    int     nbr_strikes,
    double  forward,
    double  ATMVol,
    double  maturity,
    double* weights_of_strikes)
{
    Err    err         = NULL;
    double nbr_max_std = 2; /* After this limit, strikes are unreliable */
    double maxweight   = 100;
    double std, moneyness;
    double coeff = 0.5; /* Puts more weights near the money */
    int    i;

    std = forward * ATMVol * sqrt(maturity);

    for (i = 1; i <= nbr_strikes; i++)
    {
        moneyness = fabs((strikes[i] - forward) / std);
        if (moneyness > nbr_max_std)
        {
            weights_of_strikes[i] = maxweight;
        }
        else
        {
            weights_of_strikes[i] = 1.0 + coeff * moneyness * moneyness;
        };
    }

    return err;
}

/* Interpolates the ATM vol from a market smile */
Err compute_ATMVol(
    double  forward,
    double* strikes,
    double* market_vols,
    int     nbr_strikes,
    int*    freeze_ATMVol,
    double* ATMVol)
{
    Err err = NULL;
    int i;
    i = 1;

    if (strikes[1] >= forward)
    {
        smessage("Need in the money strikes");
        return err;
    }
    else
    {
        while ((i < nbr_strikes) && (strikes[i] < forward))
        {
            i++;
        }
        if (i == nbr_strikes)
        {
            smessage("Need out of the money strikes");
            return err;
        }
        else if (strikes[i] == forward) /* if the ATMVol is given, we don't calibrate it */
        {
            *ATMVol        = market_vols[i];
            *freeze_ATMVol = 0;
            return err;
        }
        else
        /* Interpolation of the ATMVol. For the moment, linear interpolation  To be changed */
        {
            *ATMVol = (forward - strikes[i - 1]) / (strikes[i] - strikes[i - 1]) * market_vols[i] +
                      (strikes[i] - forward) / (strikes[i] - strikes[i - 1]) * market_vols[i - 1];
            *freeze_ATMVol = 1; /* We'll calibrate the ATMVol */
            return err;
        }
    }
}

Err init_SABR_parameters(
    double  forward,
    double* strikes,
    double* market_vols,
    int     nbr_strikes,
    double  ATMVol,
    double* alpha,
    int     freeze_alpha,
    double* beta,
    int     freeze_beta,
    double* rho,
    int     freeze_rho)
{
    Err err = NULL;

    return err;
}

Err opsabrcalib(
    double  forward,
    double  maturity,
    int     nbr_strikes,
    double* strikes,     /*ACHTUNG: The Vector starts at 1*/
    double* market_vols, /*ACHTUNG: The Vector starts at 1*/
    double* ATMVol,
    double* alpha,
    int     freeze_alpha, /* if 0, alpha is not calibrated */
    double* beta,
    int     freeze_beta, /* If 0, beta is not calibrated */
    double* rho,
    int     freeze_rho, /* if 0, rho is not calibrated */
    double* fitting_error)
{
    Err     err                = NULL;
    double *weights_of_strikes = NULL, *paramSABR = NULL, /* Parameters of SabrPriceGradient */
                                           *market_prices = NULL;
    long *use_paramSABR = NULL, /* 0 for frozen parameter, 1 for parameter to optimize */
        nbr_paramSABR;
    int nbr_iter = 25; /* Maximum number of iterations in Levenberg-Marquardt */
    int freeze_ATMVol;
    int i;

    /* Memory allocation */
    nbr_paramSABR      = 6; /* alpha,beta,rho,Forward,Maturity and ATMVol */
    weights_of_strikes = dvector(1, nbr_strikes);
    paramSABR          = dvector(1, nbr_paramSABR);
    market_prices      = dvector(1, nbr_strikes);
    use_paramSABR      = lngvector(1, nbr_paramSABR);

    /* Computes the ATM vol by interpolation */
    err = compute_ATMVol(forward, strikes, market_vols, nbr_strikes, &freeze_ATMVol, ATMVol);

    /* Computes a first guess for the SABR parameters */
    err = init_SABR_parameters(
        forward,
        strikes,
        market_vols,
        nbr_strikes,
        *ATMVol,
        alpha,
        freeze_alpha,
        beta,
        freeze_beta,
        rho,
        freeze_rho);

    /* Assignment of weights to strikes for fitting */
    err = compute_weights(strikes, nbr_strikes, forward, *ATMVol, maturity, weights_of_strikes);

    /* Initializes SABR parameters */
    paramSABR[1] = *alpha;
    paramSABR[2] = *beta;
    paramSABR[3] = *rho;
    paramSABR[4] = forward;
    paramSABR[5] = maturity;
    paramSABR[6] = *ATMVol;

    /* Chooses the parameters to be calibrated */
    use_paramSABR[1] = freeze_alpha;
    use_paramSABR[2] = freeze_beta;
    use_paramSABR[3] = freeze_rho;
    use_paramSABR[4] = 0; /* Not calibrating the forward */
    use_paramSABR[5] = 0; /* Not calibrating the maturity */
    use_paramSABR[6] = freeze_ATMVol;

    /* Compute the market prices */
    for (i = 1; i <= nbr_strikes; i++)
    {
        market_prices[i] =
            srt_f_optblksch(forward, strikes[i], market_vols[i], maturity, 1, SRT_CALL, PREMIUM);
    }

    /* Call Levenberg-Marcquardt */

    err = levenberg_marquardt_select(
        strikes,            /* From [1] to [nbr_strikes] */
        market_prices,      /* From [1] to [nbr_strikes] */
        weights_of_strikes, /* From [1] to [nbr_strikes] */
        nbr_strikes,
        paramSABR,     /* From [1] to [nparam] */
        use_paramSABR, /* From [1] to [nparam] */
        nbr_paramSABR,
        nbr_iter,
        SabrPriceGradient,
        fitting_error);

    /* Fills the calibrated parameters */
    *alpha  = paramSABR[1];
    *beta   = paramSABR[2];
    *rho    = paramSABR[3];
    *ATMVol = paramSABR[6];

    /* Free memory */

    free_dvector(weights_of_strikes, 1, nbr_strikes);
    free_dvector(market_prices, 1, nbr_strikes);
    free_dvector(paramSABR, 1, nbr_paramSABR);
    free_lngvector(use_paramSABR, 1, nbr_paramSABR);

    return err;
}
