/***************************************************************
 * Module:      cranalyst
 * File:        indexcmcds.h
 * Description: composite fixing adjustment for index CMCDS
 ***************************************************************/
#ifndef _indexcmcds_H
#define _indexcmcds_H

#ifdef __cplusplus
extern "C" 
//{
#endif

#include "crcrv.h"

/**
 * Composite fixing adjustment for cmcds index.
 *
 * Return SUCCESS/FAILURE.
 */

int CrxCompositeCMCDSCrv(TDate         today,            /**<(I) Trade date                           */
                         TDate         valDate,          /**<(I) Spot-date                            */
                         TDateInterval cmcdsTenor,       /**<(I) CMCDS tenor, e.g. 3Y                 */
                         const TCurve  *discZC,          /**<(I) Interest rate zero curve             */
                         const TCurve  *sprdZC,          /**<(I) Index clean curve                    */
                         int           fwdAdj,           /**<(I) Y/N flag for fwd adjustment          */
                         int           matAdj,           /**<(I) Y/N flag for maturity adjustment     */
                         int           nbFwds,           /**<(I) Number of fwds                       */
                         TDate         *fwdSttDates,     /**<(I) Start dates of fwds to be adjusted   */
                         TDate         *fwdEndDates,     /**<(I) Fwd end dates                        */
                         double        recovery,         /**<(I) Recovery rate                        */
                         TDateInterval cdsFreq,          /**<(I) CDS fee payment frequency            */
                         TDayCountConv cdsDCC,           /**<(I) CDS fee daycount                     */
                         KStubLocation feeStub,          /**<(I) CDS fee stub type                    */
                         KAccrualConv  PayAccFlag,       /**<(I) Pay accrued on default or not        */
                         KProtPayConv  protPayType,      /**<(I) Pay on default, or at maturity       */
                         TDateInterval frequency,        /**<(I) Protection leg integration freq.     */
                         TDateInterval delay,            /**<(I) Payment delay                        */
                         TCurve        **newSprdZC);     /**<(O) New Clean spread zero curve          */


/*************************************************************
 ** Generate a clean sprd curve from fowards
 **/
int CrxFwdsToCrv(TDate         today,            /**<(I) Trade date                           */
                 TDate         valDate,          /**<(I) Spot-date                            */
                 const TCurve  *discZC,          /**<(I) Interest rate zero curve             */
                 int           nbFwds,           /**<(I) Number of fwds                       */
                 TDate         *fwdSttDates,     /**<(I) Start dates of fwds to be adjusted   */
                 TDate         *fwdEndDates,     /**<(I) Fwd end dates                        */
                 double        *targetFwdRates,  /**<(I) Target fwd rates                     */
                 double        recovery,         /**<(I) Recovery rate                        */
                 TDateInterval cdsFreq,          /**<(I) CDS fee payment frequency            */
                 TDayCountConv cdsDCC,           /**<(I) CDS fee daycount                     */
                 KStubLocation dateStub,         /**<(I) CDS fee stub type                    */
                 KAccrualConv  payAccFlag,       /**<(I) Pay accrued on default or not        */
                 KProtPayConv  protPayType,      /**<(I) Pay on default, or at maturity       */
                 TDateInterval frequency,        /**<(I) Protection leg integration freq.     */
                 TDateInterval delay,            /**<(I) Payment delay                        */
                 int           nPoints,          /**<(I) number of spot dates                 */
                 TDate         *spotDates,       /**<(I) spot dates                           */
                 double        *inGuessRates,    /**<(I) guesss rates, size of nPoints        */
                 double        *actFwdRates,     /**<(IO) actual fwd rates, size of nbFwds    */
                 TCurve        **newSprdZC);     /**<(O) New Clean spread zero curve          */

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

                         
#endif    /* _indexcmcds_H */
