/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CCY_PARAM                                       */
/*                                                                            */
/*      PURPOSE:                                                              */
/*                                                                            */
/*      AUTHORS:        Eric AULD                     			      */
/*                                                                            */
/*      DATE:           16th September 1993                                   */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    Contains default settings for each currency           */
/*                      currently supported                                   */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*              Tables        , files        , and global variables accessed. */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           18th March        , 1994 */
/*                                                                            */
/*      VERSION:        <not applicable>                                      */
/*                                                                            */
/*      REASON:         Restructuring for re-use                              */
/*                                                                            */
/*      REQUEST NO:     <not applicable>                                      */
/*                                                                            */
/*      DESCRIPTION:    <not applicable>                                      */
/*                                                                            */
/******************************************************************************/

#include "swp_h_all.h"

static Err IMM_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end);

static Err AUD_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end);

static Err GBP_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end);

/* !!! this is wrong        , should be 3rd Wed - 2BD !!! */
static Err PIBOR_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                     Date *fut_start, Date *fut_end);
/* ------------------------------------------------------------------------------------------
 */
/* Dynamic ccy data - if present this takes precedence over static data defined
 * below         */
/* ------------------------------------------------------------------------------------------
 */
typedef struct {
  SrtCcyParam *this;
  void *next;
} CcyParamListItem;

static CcyParamListItem *srt_dyn_ccy_params = 0;

/* ------------------------------------------------------------------------------------------
 */
/*
        {
                DEM        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                FRF        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_365        ,
                SRT_ANNUAL        ,
                LIN_R         ,
                1        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                ITL        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R         ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                ESP        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                BEF        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                XEU        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                NLG        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                ESB        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                PTE        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                IEP        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                SRT_SEMIANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                GRD        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_365        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                ATS        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                FIM        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_ACT_360        ,
                BASIS_30_360        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
        {
                LUF        ,
                SRT_YC_CONV_SRT        ,
                MODIFIED_SUCCEEDING        ,
                MODIFIED_SUCCEEDING        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                BASIS_ACT_365        ,
                SRT_ANNUAL        ,
                LIN_R        ,
                2        ,IMM_fut        ,
                SRT_NO        ,SRT_NO        ,SRT_NO        ,SRT_NO ,SRT_NO } ,
*/

static SrtCcyParam default_ccy_params[] = {
    {AED, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {ARS, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {AUD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 1, AUD_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {BGN, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {BHD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {BRL, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {BYR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {CAD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {CHF, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {CLP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {CNY, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {COO, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_ACT, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {COP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {CZK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {DKK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {EEK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {EGP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {EUR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {GBP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 0, GBP_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {HKD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_QUARTERLY,
     LIN_R, 1, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {HRK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {HUF, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {IDR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {ILS, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_ANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {INO, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {INR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {JPB, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {JPY, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {KRO, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_QUARTERLY,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {KRW, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_QUARTERLY,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {KWD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {KZT, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {LBP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {LTL, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {LVL, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {MAD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {MXN, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_MONTHLY,
     LIN_R, 1, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {MYR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {NOK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {NZD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {OMR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {PEN, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {PHP, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 1, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {PKR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {PLN, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_ACT, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {QAR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {ROL, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {RUB, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {SAR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_360, SRT_QUARTERLY,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {SEK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {SGD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {SIT, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {SKK, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {THB, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {THO, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {TND, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {TRL, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_365, SRT_ANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {TWD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {TWO, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_QUARTERLY,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {UAH, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {USD, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {// Stanley Mrose:	11.02.2003
     UTF, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_ACT, // Nick Beckmann
     BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_SEMIANNUAL, LIN_R, 2,
     IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {VEB, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, SRT_SEMIANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {XAU, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_360, BASIS_ACT_360, BASIS_ACT_360, BASIS_30_360, SRT_ANNUAL,
     LIN_R, 2, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {ZAR, SRT_YC_CONV_SRT, MODIFIED_SUCCEEDING, MODIFIED_SUCCEEDING,
     BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, BASIS_ACT_365, SRT_QUARTERLY,
     LIN_R, 0, IMM_fut, SRT_NO, SRT_NO, SRT_NO, SRT_NO, SRT_NO},
    {LASTCCYCODE, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, 0, 0, 0, 0, 0},

};

/* --------------------------------------------------------------------------------------
 */

/* ======================================================================================
        note:
        this implementation does
        NOT always allocate any ccy_params.
        rather        , when get_ccy_param is called        ,
        a pointer to the appropriate static.
        ccy_param in memory is returned.

        Hence the free_function is sometimes a dummy function.
  ========================================================================================
*/

Err swp_f_get_CcyParam_from_CcyStr(char *CcyStr, SrtCcyParam **ccy_param) {
  int i = 0;
  Err err;
  CcyCode ccy_code;
  CcyParamListItem *item = srt_dyn_ccy_params;

  err = interp_ccy_string(CcyStr, &ccy_code);
  if (err)
    return err;

  /* Check dynamic definitions first since ccy data may be overloaded */

  while (item) {
    if (item->this->currency == ccy_code) {
      *ccy_param = item->this;
      return NULL;
    }
    item = item->next;
  }

  /* Not found in dynamic definitions        , check the hardcoded ones */
  while (default_ccy_params[i].cxxurrency != LASTCCYCODE) {
    if (ccy_code == default_ccy_params[i].cxxurrency) {
      *ccy_param = &(default_ccy_params[i]);
      return NULL;
    }
    i++;
  }

  return serror("Currency %s not recognized by SORT", CcyStr);
}
/* -------------------------------------------------------------------------------------
 */

SrtCcyParam *new_CcyParam() {
  SrtCcyParam *cps;
  cps = (SrtCcyParam *)srt_calloc(1, sizeof(SrtCcyParam));
  if (!cps)
    return NULL;
  cps->customized = SRT_YES;
  return cps;
}

/* -------------------------------------------------------------------------------------
 */

int free_CcyParam(SrtCcyParam *cps) {
  if (cps->customized)
    srt_free(cps);
  return 0;
}

/* -------------------------------------------------------------------------------------
 */

static Err ccy_set_create_dyn_params(SrtCcyParam *p) {

  CcyParamListItem *item, *temp, *prev;

  if (srt_dyn_ccy_params) {
    temp = srt_dyn_ccy_params;
    do {
      if (p->currency == temp->this->currency) {
        srt_free(temp->this);
        temp->this = p;
        return 0;
      } else {
        prev = temp;
        temp = temp->next;
      }
    } while (temp);

    {
      item = (CcyParamListItem *)srt_calloc(1, sizeof(CcyParamListItem));
      if (!item)
        return "Out of memory";
      item->this = p;
      item->next = 0;
      prev->next = item;
    }
  } else {
    item = (CcyParamListItem *)srt_calloc(1, sizeof(CcyParamListItem));
    if (!item)
      return "Out of memory";
    item->this = p;
    item->next = 0;
    srt_dyn_ccy_params = item;
  }
  return 0;
}

/* -------------------------------------------------------------------------------------
 */

void ccy_clear_dynamic_defs(void) {
  CcyParamListItem *temp = srt_dyn_ccy_params, *prev;
  while (temp) {
    srt_free(temp->this);
    prev = temp;
    temp = temp->next;
    srt_free(prev);
  }
  srt_dyn_ccy_params = 0;
  delete_dynamic_ccy_codes();
}

/* -------------------------------------------------------------------------------------
 */

Err ccy_add_definition(
    String ccy_code, SrtYCConvTyp curve_conv, BusDayConv cash_bus_day_conv,
    BusDayConv swap_bus_day_conv, BasisCode cash_basis_code,
    BasisCode m1fut_basis_code, BasisCode m3fut_basis_code,
    BasisCode swap_basis_code, SrtCompounding compd, InterpMethod interp_method,
    int spot_lag, SRT_Boolean ilv_insert_flg, SRT_Boolean ilv_no_overwrite_flg,
    SRT_Boolean toy_insert_flg, SRT_Boolean toy_no_overwrite_flg) {
  SrtCcyParam *p = (SrtCcyParam *)srt_calloc(1, sizeof(SrtCcyParam));
  Err err;
  CcyCode ccy;
  if (!p)
    return "Out of memory";

  /* Temporary adjustment - do not override existing ccy definitions */
  /*******************************************************************/
  if (!interp_default_ccy_string(ccy_code, &ccy)) {
    srt_free(p);
    return 0;
  }
  /*******************************************************************/

  err = ccy_get_or_create_code(ccy_code, &ccy);
  if (err)
    return err;
  p->currency = ccy;
  p->curve_conv = curve_conv;
  p->cash_bus_day_conv = cash_bus_day_conv;
  p->swap_bus_day_conv = swap_bus_day_conv;
  p->cash_basis_code = cash_basis_code;
  p->m1fut_basis_code = m1fut_basis_code;
  p->m3fut_basis_code = m3fut_basis_code;
  p->swap_basis_code = swap_basis_code;
  p->compd = compd;
  p->interp_method = interp_method;
  p->spot_lag = spot_lag;
  p->ilv_insert_flg = ilv_insert_flg;
  p->ilv_no_overwrite_flg = ilv_no_overwrite_flg;
  p->toy_insert_flg = toy_insert_flg;
  p->toy_no_overwrite_flg = toy_no_overwrite_flg;
  p->customized = SRT_FALSE;

  switch (p->currency) {
  case GBP:
    p->fut_func = GBP_fut;
    break;
  case AUD:
    p->fut_func = AUD_fut;
    break;
  default:
    p->fut_func = IMM_fut;
    break;
  }

  return ccy_set_create_dyn_params(p);
}

/* -------------------------------------------------------------------------------------
 */

Err ccy_string_get_param(SrtCcyParam *ccy_param, String param_name,
                         String *param_val) {

  Err err;
  Message val;

  if (!strcmp(param_name, "CASHBUSDAYCONV")) {
    val = ccy_param->cash_bus_day_conv;
    if (err = (translate_bus_day_conv(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "SWAPBUSDAYCONV")) {
    val = ccy_param->swap_bus_day_conv;
    if (err = (translate_bus_day_conv(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "CASHBASIS")) {
    val = ccy_param->cash_basis_code;
    if (err = (translate_basis(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "SWAPBASIS")) {
    val = ccy_param->swap_basis_code;
    if (err = (translate_basis(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "M1FUTBASIS")) {
    val = ccy_param->m1fut_basis_code;
    if (err = (translate_basis(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "M3FUTBASIS")) {
    val = ccy_param->m3fut_basis_code;
    if (err = (translate_basis(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "COMPOUNDING")) {
    val = ccy_param->compd;
    if (err = (translate_compounding(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "INTERPMETHOD")) {
    val = ccy_param->interp_method;
    if (err = (translate_interp_method(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "INTERLEAVE")) {
    val = ccy_param->ilv_insert_flg;
    if (err = (translate_insert_flg(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "TURNOFYEAR")) {
    val = ccy_param->toy_insert_flg;
    if (err = (translate_insert_flg(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "INTERLEAVEOVERWRITE")) {
    val = ccy_param->ilv_no_overwrite_flg;
    if (err = (translate_no_overwrite_flg(param_val, val)))
      return err;
    return 0;
  }
  if (!strcmp(param_name, "TURNOFYEAROVERWRITE")) {
    val = ccy_param->toy_no_overwrite_flg;
    if (err = (translate_no_overwrite_flg(param_val, val)))
      return err;
    return 0;
  }

  return serror("Don't know param_name %s", param_name);
}

/* -------------------------------------------------------------------------------------
 */

Err ccy_string_set_param(SrtCcyParam *ccy_param, String param_name,
                         String param_val) {
  Err err;
  Message val;

  if (!strcmp(param_name, "CASHBUSDAYCONV")) {
    if (err = (interp_bus_day_conv(param_val, &val)))
      return err;
    ccy_param->cash_bus_day_conv = val;
    return 0;
  }
  if (!strcmp(param_name, "SWAPBUSDAYCONV")) {
    if (err = (interp_bus_day_conv(param_val, &val)))
      return err;
    ccy_param->swap_bus_day_conv = val;
    return 0;
  }
  if (!strcmp(param_name, "CASHBASIS")) {
    if (err = (interp_basis(param_val, &val)))
      return err;
    ccy_param->cash_basis_code = val;
    return 0;
  }
  if (!strcmp(param_name, "SWAPBASIS")) {
    if (err = (interp_basis(param_val, &val)))
      return err;
    ccy_param->swap_basis_code = val;
    return 0;
  }
  if (!strcmp(param_name, "M1FUTBASIS")) {
    if (err = (interp_basis(param_val, &val)))
      return err;
    ccy_param->m1fut_basis_code = val;
    return 0;
  }
  if (!strcmp(param_name, "M3FUTBASIS")) {
    if (err = (interp_basis(param_val, &val)))
      return err;
    ccy_param->m3fut_basis_code = val;
    return 0;
  }
  if (!strcmp(param_name, "COMPOUNDING")) {
    if (err = (interp_compounding(param_val, &val)))
      return err;
    ccy_param->compd = val;
    return 0;
  }
  if (!strcmp(param_name, "INTERPMETHOD")) {
    if (err = (interp_interp_method(param_val, &val)))
      return err;
    ccy_param->interp_method = val;
    return 0;
  }
  if (!strcmp(param_name, "INTERLEAVE")) {
    if (err = (interp_insert_flg(param_val, &val)))
      return err;
    ccy_param->ilv_insert_flg = val;
    return 0;
  }
  if (!strcmp(param_name, "TURNOFYEAR")) {
    if (err = (interp_insert_flg(param_val, &val)))
      return err;
    ccy_param->toy_insert_flg = val;
    return 0;
  }
  if (!strcmp(param_name, "INTERLEAVEOVERWRITE")) {
    if (err = (interp_no_overwrite_flg(param_val, &val)))
      return err;
    ccy_param->ilv_no_overwrite_flg = val;
    return 0;
  }
  if (!strcmp(param_name, "TURNOFYEAROVERWRITE")) {
    if (err = (interp_no_overwrite_flg(param_val, &val)))
      return err;
    ccy_param->toy_no_overwrite_flg = val;
    return 0;
  }

  return serror("Don't know param name: %s", param_name);
}

/* -------------------------------------------------------------------------------------
 */

/**********************************************************************
        future date functions
***********************************************************************/

/*
        Futures:


        futures written like
        "sep94m" monthly
        "sep94q" quarterly

        *nm set to 1 montly        , 3 quarterly

        assumes "70"is 1970 but "69"is 2069.

*/

/* last trading day is 2bd before 3rd Wednesday */
static Err IMM_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end) {
  Date d1;

  /* find third wednesday */
  d1 = third_wednesday(sm, sy);
  *fut_last_trading = d1 - 2; /* FIX WHEN WE HAVE HOLIDAYS */
  *fut_start = d1 - 2 + spot;

  /* monthly */
  if (nm == 1) {
    sm += 1;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = third_wednesday(sm, sy);
    *fut_end = d1 - 2 + spot;
  }
  /* quarterly */
  else if (nm == 3) {
    sm += 3;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = third_wednesday(sm, sy);
    *fut_end = d1 - 2 + spot;
  }
  return 0;
}

/* -------------------------------------------------------------------------------------
 */

/* last trading day is 1bd before 2nd Friday */
static Err AUD_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end) {
  Date d1;

  /* find second friday */
  d1 = second_friday(sm, sy);
  *fut_last_trading = d1 - 1; /* FIX WHEN WE HAVE HOLIDAYS */
  *fut_start = d1 - 1 + spot;

  /* monthly */
  if (nm == 1) {
    sm += 1;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = second_friday(sm, sy);
    *fut_end = d1 - 1 + spot;
  }
  /* quarterly */
  else if (nm == 3) {
    sm += 3;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = second_friday(sm, sy);
    *fut_end = d1 - 1 + spot;
  }
  return 0;
}

/* -------------------------------------------------------------------------------------
 */

/**** GBP futures trade until the third wednesday; spot is today ***/
static Err GBP_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                   Date *fut_start, Date *fut_end) {
  Date d1;
  d1 = third_wednesday(sm, sy);
  *fut_last_trading = d1;
  *fut_start = d1;

  /* monthly */
  if (nm == 1) {
    sm += 1;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = third_wednesday(sm, sy);
    *fut_end = d1;
  }
  /* quarterly */
  else if (nm == 3) {
    sm += 3;
    if (sm > 12) {
      sy += 1;
      sm -= 12;
    }
    d1 = third_wednesday(sm, sy);
    *fut_end = d1;
  }
  return 0;
}

/* -------------------------------------------------------------------------------------
 */

/* last trading day is 2bd before 11th Thursday of quarter        ,spot is 1 day
 */
static Err PIBOR_fut(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                     Date *fut_start, Date *fut_end) {

  /* FIX for holidays */

  if (nm != 3)
    return serror("Only understand QUARTERLY Pibor futures.");

  *fut_last_trading = eleventh_thursday(sm - 2, sy) - 2;
  *fut_start = *fut_last_trading + spot;
  *fut_end = eleventh_thursday((sm + 1 < 13 ? sm + 1 : 1),
                               (sm + 1 < 13 ? sy : sy + 1)) -
             2 + spot;

  return 0;
}

/* -------------------------------------------------------------------------------------
 */
