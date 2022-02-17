/* ===================================================================================

   FILENAME:        utcurrency.h

   PURPOSE:         the types and structures needed to deal with a currency

   ==================================================================================
 */

#ifndef UTCURRENCY_H
#define UTCURRENCY_H

/* ---------------------------------------------------------------------------
                          CURRENCY CODE
   ---------------------------------------------------------------------------
 */

/*
typedef enum 	CcyCodeEnum_ {
        GBP      ,
        DEM      , DM=DEM      ,
        FRF      , FFR=FRF      ,
        ITL      ,
        USD      ,
        JPY      ,
        ESP      ,
        BEF      ,
        LUF      ,
        CAD      ,
        CHF      ,
        XEU      , ECU = XEU      ,
        SEK      ,
        DKK      ,
        NLG      ,
        ATS      ,
        AUD      ,
        FIM      ,
        HKD      ,
        HUF      ,
        RUB      ,
        NZD      ,
        ESB      ,
        PTE      ,
        EUR      ,
        IDR      ,
        IEP      ,
        MYR      ,
        SGD      ,
        THB      ,
        TWD      ,
        BHD      ,
        GRD      ,
        KWD      ,
        NOK      ,
        OMR      ,
        QAR      ,
        SAR      ,
        ZAR      ,
        CNY      ,
        CZK      ,
        INR      ,
        JPB      ,
        KRW      ,
        MXN      ,
        PHP      ,
        TRL      ,
        PLN      ,
        ARS      ,
        BRL      ,
        CLP      ,
        COP      ,
        VEB      ,
        KRO      ,
        DEFAULT_CCY      ,
        LASTCCYCODE
}CcyCodeEnum;
*/

typedef enum CcyCodeEnum_ {
  AED,
  ARS,
  AUD,
  BGN,
  BHD,
  BRL,
  BYR,
  CAD,
  CHF,
  CLP,
  CNY,
  COO,
  COP,
  CZK,
  DKK,
  EEK,
  EGP,
  EUR,
  GBP,
  HKD,
  HRK,
  HUF,
  IDR,
  ILS,
  INO, // P McCallum: 23.10.2003
  INR,
  JPB,
  JPY,
  KRO,
  KRW,
  KWD,
  KZT,
  LBP,
  LTL,
  LVL,
  MAD,
  MXN,
  MYR,
  NOK,
  NZD,
  OMR,
  PEN,
  PHP,
  PKR,
  PLN,
  QAR,
  ROL,
  RUB,
  SAR,
  SEK,
  SGD,
  SIT,
  SKK,
  THB,
  THO,
  TND,
  TRL,
  TWD,
  TWO,
  UAH,
  USD,
  UTF, // Stanley Mrose: 12.02.2003
  VEB,
  XAU,
  ZAR,
  DEFAULT_CCY,
  LASTCCYCODE
} CcyCodeEnum;

typedef long CcyCode;

/* ---------------------------------------------------------------------------
                                CURRENCIES TRANSLATION
   ---------------------------------------------------------------------------
 */

Err ccy_get_or_create_code(String cs, CcyCode *ccy);

void delete_dynamic_ccy_codes(void);

Err interp_ccy_string(String cs, CcyCode *ccy);

int interp_default_ccy_string(String cs, CcyCode *ccy);

/* Copy the output string - since if ccy dynamic it can be deleted */
Err translate_ccy(String *cs, CcyCode ccy);

int translate_default_ccy(String *cs, CcyCode ccy);

#endif
