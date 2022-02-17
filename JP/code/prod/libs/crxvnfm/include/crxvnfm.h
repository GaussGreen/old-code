/******************************************************************************
 * VNFM approximation implementation for credit spreads
 * Converts spot vols to Black-Scholes CDS swaption vols
 *****************************************************************************/

#ifndef _crxvnfm_h
#define _crxvnfm_h

#ifdef  __cplusplus
extern "C" {
#endif

/******************************************************************************
 * CrxCreditVNFMBootstrapSpotVols
 * Bootstraps a set of benchmark CDS and corresponding Black-Scholes vols
 * to an equivalent set of spot vols using the credit VNFM approximation.
 * The vol must be an implied BS vol - if your vols are multi-q, then you
 * should convert them before calling this function.
 * Note that the annuities and par spreads calculated are approximate - 
 * they do not include accrual recovery on default - so do not use them other
 * than as approximations.
 *****************************************************************************/
int CrxCreditVNFMBootstrapSpotVols (
    TDate valueDate,
    /* CREDIT VOL PARAMETERS */
    int           numVolCR,      /**<(I) Num. of input credit vol points     */
    const TDate*  volDateCR,     /**<(I) Credit vol dates                    */
    double*       volCR,         /**<(I/O) Input BS vols, output spot if cal */
    int*          volUsedCR,     /**<(I/O) Input/output 1 if should/did cal  */
    int           skipFlagCR,    /**<(I) If TRUE, skip & go on if fail BM cal*/
    TDateInterval cdsFreq,       /**<(I) Underlying BM CDS fee frequency     */
    long          cdsDCC,        /**<(I) Underlying BM CDS fee day-count     */
    double        recovery,      /**<(I) Recovery Rate                       */
    const TDate*  cdsStartDates, /**<(I) Underlying BM CDS start dates       */
    const TDate*  cdsEndDates,   /**<(I) Underlying BM CDS maturity dates    */
    double        betaCR,        /**<(I) Credit mean-reversion               */
    double        backboneCR,    /**<(I) CR backbone: 0=normal, 1=lognormal  */
    double        factorWeightCR,/**<(I) CR alpha factor weight              */
    const TCurve* crCurve,       /**<(I) Credit clean spread curve           */
    double        leftCR2Q,      /**<(I) Left Credit 2Q smile (1=lognormal)  */
    double        rightCR2Q,     /**<(I) Right Credit 2Q smile (1=lognormal) */
    double        CR2QF,         /**<(I) Forward shift 'F' parameter for 2Q  */
    double        corrIRCR,      /**<(I) IR/CR spot-vol correlation          */
    int           numVolIR,      /**<(I) Num. of IR spot vol points          */
    const TDate*  volDateIR,     /**<(I) IR vol dates                        */
    const double* volIR,         /**<(I) IR spot vols                        */
    double        betaIR,        /**<(I) IR mean-reversion                   */
    double        backboneIR,    /**<(I) IR backbone: 0=normal, 1=lognormal  */
    double        factorWeightIR,/**<(I) IR alpha factor weight              */
    const TCurve* irCurve,       /**<(I) IR zero curve                       */
    double*       parSpreadAry,  /**<(O)Output BM par spreads, where cal     */
    double*       annuityAry     /**<(O)Output BM annuity value, where cal   */
    );

/******************************************************************************
 * CrxCreditVNFMBlackScholesVol
 * Given credit and interest rate spot vols, mean reversions and
 * correlations, returns the Black-Scholes implied volatility of a given
 * swaption.
 *****************************************************************************/
int CrxCreditVNFMBlackScholesVol (
    TDate valueDate,
    TDateInterval cdsFreq,       /**<(I) Underlying BM CDS fee frequency     */
    long          cdsDCC,        /**<(I) Underlying BM CDS fee day-count     */
    double        recovery,      /**<(I) Recovery Rate                       */
    TDate         optionExpiryDate, /**<(I) Expiry of option                 */
    TDate         cdsStartDate,  /**<(I) Underlying BM CDS start dates       */
    TDate         cdsEndDate,    /**<(I) Underlying BM CDS maturity dates    */
    /* CREDIT VOL, CURVE, AND SMILE PARAMETERS */
    int           numVolCR,      /**<(I) Num. of input credit vol points     */
    const TDate*  volDateCR,     /**<(I) Credit spot vol dates               */
    double*       volCR,         /**<(I) Credit spot vols                    */
    double        betaCR,        /**<(I) Credit mean-reversion               */
    double        backboneCR,    /**<(I) CR backbone: 0=normal, 1=lognormal  */
    double        factorWeightCR,/**<(I) CR alpha factor weight              */
    const TCurve* crCurve,       /**<(I) Credit clean spread curve           */
    double        leftCR2Q,      /**<(I) Left Credit 2Q smile (1=lognormal)  */
    double        rightCR2Q,     /**<(I) Right Credit 2Q smile (1=lognormal) */
    double        CR2QF,         /**<(I) Forward shift 'F' parameter for 2Q  */
    /* IR VOL AND CURVE PARAMETERS */
    int           numVolIR,      /**<(I) Num. of IR spot vol points          */
    const TDate*  volDateIR,     /**<(I) IR vol dates                        */
    const double* volIR,         /**<(I) IR spot vols                        */
    double        betaIR,        /**<(I) IR mean-reversion                   */
    double        backboneIR,    /**<(I) IR backbone: 0=normal, 1=lognormal  */
    double        factorWeightIR,/**<(I) IR alpha factor weight              */
    const TCurve* irCurve,       /**<(I) IR zero curve                       */
    /* IR/CR CORRELATION */
    double        corrIRCR,      /**<(I) IR/CR spot-vol correlation          */
    double*       bsVol          /**<(O)Output Black-Scholes volatility      */
    );

#ifdef __cplusplus
}
#endif

#endif