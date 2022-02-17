#include "crtdbg.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "MBSPPFuncs.h"
#include "MBSPTUtil.h"
#include "MBSPrepay.h"
#include "TBAUtil.h"

static long noleap[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static long leap[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static long cumdays[12] = {0,   31,  59,  90,  120, 151,
                           181, 212, 243, 273, 304, 334};

/*****  Manager  ************************************************************/
/*
 *       Manage input & output  , reading from ascii files into data structures.
 */
char *Manager(char *data_dir,
              TERM_DATA *term_data,     /* Structure of term structure data  */
              TREE_DATA *tree_data,     /* Structure of tree data */
              DEAL_DATA *deal_data,     /* Structure of deal structure data  */
              PREPAY_DATA *prepay_data, /* Structure of prepayment data */
              HISTORY *history,         // historical index values
              MBSDENSITY *density) // density function for prepay propensity
{
  long Date;  /* Working date (YYYYMMDD) */
  int F, wam, /* Frequency of underlying payment as a integer */
      i, j, k;
  double ratio1, ratio2;
  char string[MAXBUFF];
  FILE *stream;
  char *err;

  stream = _OpenMBSFile(data_dir, FileTerm);
  if (!stream)
    return (nrerror("Can't open Term.dat")); /* Open the yield curve data file
                                                (see termodel.h for a
                                                description of the inputs) */

  fgets(string, 80, stream); /* Read the comment line */
  fgets(string, 80, stream);
  // fscanf (stream  , "%ld \n"  , &(term_data->Today));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%ld \n"  , &(term_data->Spotdays));
  fgets(string, 80, stream);
  if (fscanf(stream, "%d \n", &(term_data->MMB)) != 1)
    return (nrerror("Can't read term_data.MMB"));
  fgets(string, 80, stream);
  if (fscanf(stream, "%s \n", &(term_data->AoS)) != 1)
    return (nrerror("Can't read term_data.AoS"));
  fgets(string, 80, stream);
  if (fscanf(stream, "%s \n", term_data->YearBase) != 1)
    return (nrerror("Can't read term_data.YearBase"));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf(stream  , "%lf\n"  , &(deal_data->Refirate)); /* today's refinance
  // rate */
  fgets(string, 80, stream);
  if (fscanf(stream, "%ld\n", &(term_data->TStree)) != 1)
    return (nrerror(
        "Can't read term_data.TStree")); /* treasury tree or swap tree */
  fgets(string, 80, stream);
  if (fscanf(stream, "%ld\n", &(term_data->TSdiscount)) != 1)
    return (
        nrerror("Can't read term_data.TSdiscount")); /* discount using treasury
                                                        or swap rates?*/
  fgets(string, 80, stream);
  for (i = 0; i < 6; i++)
    fgets(string, 80, stream);
  //        fscanf (stream  , "%lf \n"  , &(term_data->MMYield[i]));
  fgets(string, 80, stream);
  for (i = 0; i < 8; i++) {
    fgets(string, 80, stream);
    //        fscanf (stream  , "%lf \n"  , &(term_data->SwapYield[i])); /*
    //        Treasury yield for this maturity */ term_data->YieldS[i] =
    //        term_data->SwapYield[i];
  } /* end for */
  fgets(string, 80, stream);
  for (i = 0; i < 8; i++) {
    fgets(string, 80, stream);
    //   fscanf (stream  , "%lf \n"  , &(spread));                           /*
    //   Swap spread for this maturity */ if (term_data->TSdiscount ==1)
    //                   term_data->YieldS[i] += spread/100.; /* if Swap
    //                   discount  , add spread */
    //   if (term_data->TStree == 1) term_data->SwapYield[i] += spread / 100.;
    //   /* if Swap tree  , add spread to Yield in tree building */
  } /* for i */

  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%lf \n"  , &(term_data->Beta));
  fgets(string, 80, stream);
  for (i = 0; i < NBSWAPTION; i++)
    fgets(string, 80, stream);
  //        fscanf (stream  , "%lf %d \n"  ,
  //                &(term_data->SwaptionVol[i])  ,
  //                &(term_data->SwaptionExp[i]));

  stream = _OpenMBSFile(data_dir, FileDeal);
  if (!stream)
    return (nrerror(
        "Can't open Deal.dat")); /* open the deal data file (see struct.h for a
                                    description of the inputs) */
  for (i = 0; i < 15; i++)
    fgets(string, 80, stream); /* read comments and casenumber */
  fgets(string, 80, stream);
  // fscanf(stream  , "%ld\n"  , &(deal_data->casenumber));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf(stream  , "%c\n"  , &(deal_data->Scenario));                 /* for
  // calculating duration  , we parallel shift yield curve. */
  fgets(string, 80, stream);
  if (fscanf(stream, "%lf\n", &(deal_data->Shift)) != 1)
    return (nrerror("Can't read deal_data.Shift"));
  if (deal_data->Scenario == 'Y') {
    for (i = 0; i < 6; i++) {
      term_data->MMYield[i] += deal_data->Shift;
    }
    for (i = 0; i < 8; i++) {
      term_data->SwapYield[i] += deal_data->Shift;
      term_data->YieldS[i] += deal_data->Shift;
    }
  }

  if (deal_data->Refirate < 0.0)
    deal_data->Refirate = term_data->SwapYield[5] - deal_data->Refirate;
  /* negative refirate means negative spread over 10yr treasury yield */

  fgets(string, 80, stream);
  if (fscanf(stream, "%d \n", &(tree_data->Ppy)) != 1)
    return (nrerror("Can't read tree_data.Ppy"));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%d \n"  , &(term_data->IndexMat));
  fgets(string, 80, stream);
  if (fscanf(stream, "%c \n", &(term_data->IndexFreq)) != 1)
    return (nrerror("Can't read term_data.IndexFreq"));

  deal_data->NbExer =
      2; /* 2 "exercise dates": Lockout date and maturity of the SFS */
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  //  fscanf (stream  , "%ld \n"  , &(deal_data->Exer[1])); /* The maturity of
  //  the deal is stored as the second exercise date */
  deal_data->Strike[1] = 0.; /* Bogus value. Will never be used */
  deal_data->Maturity =
      Nxtmth(deal_data->Exer[1], /* The final maturity is the maturity of the
                                    deal plus the index maturity */
             (long)term_data->IndexMat * 12, (long)1);

  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%lf \n"  , &(deal_data->Coupon));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%lf \n"  , &(deal_data->Gwac)); /*gross coupon*/
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%ld \n"  , &(deal_data->pay_delay));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%d \n"  , &(deal_data->Term));
  tree_data->NodeMax =
      (int)sqrt(20. * (double)tree_data->Ppy) *
      Nsigma; /* Maximun number of nodes in spatial direction. */
  fgets(string, 80, stream);
  if (fscanf(stream, "%c \n", &(deal_data->Freq)) != 1)
    return (nrerror("Can't read deal_data.Freq"));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%lf \n"  , &(term_data->OAS));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->Ppform)) != 1)
    return (nrerror(
        "Can't read deal_data.Ppform")); /* ratio or absolute shift ppf */
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->AlorCPR)) != 1)
    return (
        nrerror("Can't read deal_data.AlorCPR")); /* avg life or annual CPR */

  fgets(string, 80, stream); /* Read the prepayment function */
  fgets(string, 80, stream);
  /*      fgets  (string  , 80  , stream); */

  if (deal_data->AlorCPR == 'A') {
    /* if avg life  , there is only one column*/
    for (i = 0; i < NBPOINTS; i++) {
      if (fscanf(stream, "%lf %lf\n", &(deal_data->IndexLevel[i]),
                 &(deal_data->AmortLevel[0][i])) != 2)
        return (nrerror("Can't read deal_data.IndexLevel or AmortLevel"));
      deal_data->AmortLevel1[0][i] = deal_data->AmortLevel[0][i];
      for (j = 1; j < MAXNUMP; j++) {
        deal_data->AmortLevel[j][i] = deal_data->AmortLevel[0][i];
        deal_data->AmortLevel1[j][i] = deal_data->AmortLevel[j][i];
      }
    }
  } else { /* if CPR  , there are MAXNUMP columns */
    for (i = 0; i < NBPOINTS; i++) {
      if (fscanf(stream, "%lf", &(deal_data->IndexLevel[i])) != 1)
        return (nrerror("Can't read deal_data.IndexLevel"));

      for (j = 0; j < MAXNUMP; j++) {
        if (fscanf(stream, "%lf", &(deal_data->AmortLevel[j][i])) != 1)
          return (nrerror("Can't read deal_data.AmortLevel"));
        deal_data->AmortLevel1[j][i] = deal_data->AmortLevel[j][i];
      }
      fscanf(stream, "\n");
    }
  }

  if (deal_data->Ppform == 'A' && deal_data->Scenario == 'Y') {
    for (i = 0; i < NBPOINTS; i++)
      deal_data->IndexLevel[i] -= deal_data->Shift;
  } /* if there is a shift of y.c.  , we should shift scenario */
  if (deal_data->Ppform == 'R') {  // if ratio  , we change from refi/gwac ratio
    for (i = 0; i < NBPOINTS; i++) // to shift in 10yr yield
      deal_data->IndexLevel[i] =
          transformIncentive(deal_data->IndexLevel[i], deal_data->Gwac,
                             deal_data->Refirate, term_data->SwapYield[5]);
  } else
    for (i = 0; i < NBPOINTS; i++)
      deal_data->IndexLevel[i] +=
          term_data->SwapYield[5]; /* turn relative index into absolute value */

  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf (stream  , "%lf\n"  , &(deal_data->Speedup)); /* Read speed up
  // factor of PPF */
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  //  fscanf (stream  , "%lf\n"  , &(deal_data->rational)); /* Read rational
  //  speed up */
  fgets(string, 80, stream);
  if (fscanf(stream, "%c \n", &(deal_data->Hedge)) != 1)
    return (nrerror("Can't read deal_data.Hedge")); /* Compute hedge or not */
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->Schedule)) != 1)
    return (nrerror(
        "Can't read deal_data.Schedule")); /* scheduled amortization or not */
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->Balloon)) != 1)
    return (nrerror("Can't read dal_data.Balloon")); /* Balloon */
  fgets(string, 80, stream);
  if (fscanf(stream, "%ld\n", &(deal_data->BalloonSch)) != 1)
    return (nrerror("Can't read deal_data.BalloonSch"));
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->Always)) != 1)
    return (nrerror("Can't read deal_data.BalloonSch"));
  if (deal_data->Always == 'N' &&
      ((deal_data->Exer[1] - term_data->Today) / 10000 > 15.))
    deal_data->Maturity = Nxtmth(deal_data->Exer[1], (long)12, (long)1);
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf(stream  , "%lf\n"  , &(deal_data->Inpio));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf(stream  , "%lf\n"  , &(deal_data->Inppo));
  fgets(string, 80, stream);
  fgets(string, 80, stream);
  // fscanf(stream  , "%lf\n"  , &(deal_data->collatoral));
  fgets(string, 80, stream);
  if (fscanf(stream, "%lf\n", &(deal_data->stripratio)) != 1)
    return (nrerror("Can't read deal_data.stripio"));
  fgets(string, 80, stream);
  if (fscanf(stream, "%lf\n", &(deal_data->implvol)) != 1)
    return (nrerror("Can't read deal_data.impvol"));
  fgets(string, 80, stream);
  if (fscanf(stream, "%c\n", &(deal_data->invfloater)) != 1)
    return (nrerror("Can't read deal_data.invfloater"));

  stream = _OpenMBSFile(data_dir, FilePrepay);
  if (!stream)
    return (nrerror("Can't open Prepay.dat"));
  /* open the prepayment data file (see termodel.h for a description of the
   * inputs) */

  prepay_data->Spot = 0.; /* This is the one factor version of the model: we
                             don't use the second dimension */
  prepay_data->SpotVol = 0.;
  prepay_data->Beta = 0.;

  prepay_data->NbFwd = 1; /* We don't use any forward curve for the prepayment
                             so we put bogus values here */
  prepay_data->FwdMat[0] = deal_data->Exer[1];
  prepay_data->FwdPrepay[0] = 0.;

  fgets(string, 80, stream);
  if (fscanf(stream, "%d \n", &(deal_data->SFSNb)) != 1)
    return (nrerror("Can't read deal_data.SFSNb"));

  fgets(string, 80, stream); /* First group of SFS */
  fgets(string, 80, stream);
  ratio1 = 0.96 * deal_data->Gwac / deal_data->Refirate;
  ratio2 = 0.20 * ratio1 / 0.95;
  j = 0;
  k = 0;
  if (fscanf(stream, "%lf \n", &(deal_data->Amort[0])) != 1)
    return (nrerror("Can't read deal_data.Amort"));
  deal_data->Amort[0] /= 100.;
  for (i = 0; i < deal_data->SFSNb; i++) {
    if ((0.96 * deal_data->Gwac + REFI_THRESHOLD - deal_data->Refirate) <=
        0.01) {
      deal_data->TriggerRate[i] =
          ratio1 + /* use linear spacing*/
          (ratio2 - ratio1) * (double)i / (double)(2 * deal_data->SFSNb - i);
      deal_data->TriggerRate[i] *= deal_data->Refirate;
    } else { /* otherwise cubic spacing */
      ratio1 = 0.96 * deal_data->Gwac + REFI_THRESHOLD - deal_data->Refirate;
      ratio2 = pow((term_data->SwapYield[5] - 1.0) / ratio1, 1. / 3.) + 1.0;
      deal_data->TriggerRate[i] =
          deal_data->Refirate - 1.0 +
          ratio1 * pow(1.0 - (double)i / (double)deal_data->SFSNb * ratio2, 3);
    }
    deal_data->TriggerRate[i] -= deal_data->Refirate;
    deal_data->TriggerRate[i] +=
        term_data->SwapYield[5]; /* change from refirate/gwac to 10yr yield */
    deal_data->Amort[i] = deal_data->Amort[0];
  } /* for i */

  fgets(string, 80, stream); /* Second group of SFS */
  fgets(string, 80, stream);

  if (fscanf(stream, " %lf \n", &(deal_data->Amort[deal_data->SFSNb])) != 1)
    return (nrerror("Can't read deal_data.Amort[deal_data.SFSNb]"));
  deal_data->Amort[deal_data->SFSNb] /= 100.;
  for (i = deal_data->SFSNb; i < 2 * deal_data->SFSNb; i++) {
    deal_data->TriggerRate[i] = deal_data->TriggerRate[i - deal_data->SFSNb];
    deal_data->Amort[i] = deal_data->Amort[deal_data->SFSNb];
  } /* for i */

  fgets(string, 80, stream);
  if (fscanf(stream, " %ld \n", &(deal_data->Nbweights)) != 1)
    return (nrerror("Can't read deal_data.Nbweights"));

  fgets(string, 80, stream);
  for (i = 0; i < NUMRAMPS; i++) {
    if (fscanf(stream, "%lf %lf \n", &(deal_data->SeasoningLevel[i]),
               &(deal_data->SeasoningMat[i])) != 2)
      return (nrerror(
          "Can't read either deal_data.SeasoningLevel or SeasoningMat"));
    deal_data->SeasoningLevel[i] =
        transformIncentive(deal_data->SeasoningLevel[i], deal_data->Gwac,
                           deal_data->Refirate, term_data->SwapYield[5]);
  }
  fgets(string, 80, stream);
  for (i = 0; i < 12; i++)
    if (fscanf(stream, "%lf \n", &(deal_data->Seasonality[i])) != 1)
      return (nrerror("Can't read deal_data.Seasonality"));

  fgets(string, 80, stream);
  if (fscanf(stream, "%d \n", &(deal_data->Delay)) != 1)
    return (nrerror("Can't read deal_data.Delay"));

  fgets(string, 80, stream);
  if (fscanf(stream, "%d \n", &(deal_data->NbFixedPrepay)) != 1)
    return (nrerror("Can't read deal_data.NbFixedPrepay"));
  fgets(string, 80, stream);
  for (i = 0; i < deal_data->NbFixedPrepay; i++) {
    if (fscanf(stream, "%lf \n", &(deal_data->FixedPrepay[i])) != 1)
      return (nrerror("Can't read deal_data.FixedPrepay"));
    deal_data->FixedPrepay[i] /= 100.;
    // convert from cpr to smm
    if (deal_data->FixedPrepay[i] > 1.0)
      return (nrerror("Some FixedPrepay > 1.0!"));
    deal_data->FixedPrepay[i] =
        1.0 - pow(1.0 - deal_data->FixedPrepay[i], 1.0 / 12.0);
    deal_data->InitialSMM[i] = deal_data->FixedPrepay[i];

  } /* for i */
  // convert to principal paydown rate
  wam = month_diff(term_data->Today, deal_data->Exer[1]);

  //#ifndef MBSPTFIXEDPREPAYISCPR
  //	if( err = get_initial_principal_payments( deal_data->NbFixedPrepay  ,
  //deal_data->FixedPrepay  , wam  , deal_data->Gwac  , deal_data->NbFixedPrepay
  //, deal_data->InitialSMM )) return(err); #endif

  stream = _OpenMBSFile(TBADATADIR, FileHolidays);
  if (!stream)
    return (nrerror("Can't open HOLIDAY.dat"));

  term_data->numHolidays = 0;
  while (fgets(string, 80, stream)) {
    string[8] = '\0';
    term_data->Holidays[term_data->numHolidays++] = atol(string);
  }

  F = Conv_Freq(deal_data->Freq);

  term_data->ValueDate = Nxtbusday(
      term_data
          ->Today, /* Value date is today shifted by (Spotdays) business days */
      term_data->Spotdays, term_data->Holidays, term_data->numHolidays);

  Date = deal_data->Maturity;
  do
    Date = Nxtmth(Date, /* Calculate the coupon date falling before today by
                           going backward from */
                  (long)(-12 / F), /* maturity one coupon at a time. */
                  (long)1);
  while (Date > term_data->ValueDate);

  i = max(0, deal_data->NbFixedPrepay - deal_data->Delay) +
      1; /* Calculate the lockout date  , that is the date on which we do the
            last    */
  deal_data->Exer[0] = Nxtmth(Date, /* reading of the index. It is located
                                       NbFixedPrepay - Delay + 1 coupons   */
                              (long)i * 12 / F, /* after today. */
                              (long)1);
  deal_data->Exer[0] = min(deal_data->Exer[0], deal_data->Exer[1]);
  deal_data->Strike[0] = 0.; /* Bogus value. Will never be used */

  // fclose (stream);

  // if (deal_data->SFSNb > 1)//for the case we use new prepayment function
  {
    //	if( err=manager1( data_dir  , history  ,density ))
    //return(nrerror(err));/* read input */
    if (err = history_reader(data_dir, history))
      return (err);
    if (err = density_reader(data_dir, density))
      return (err);
    adjust_density(density, deal_data->Speedup * deal_data->rational,
                   deal_data->Gwac, deal_data->Refirate,
                   term_data->SwapYield[5], deal_data->Nbweights);
    adjust_history(history, deal_data->Term, deal_data->Refirate,
                   term_data->SwapYield[5]);
  }
  return (0);

} /* Manager */

/*****  Input_Check
 * ************************************************************/
/*
 *       Check that the input read in Manager() are correct.
 */
char *Input_Check(TERM_DATA *term_data, /* Structure of term structure data  */
                  DEAL_DATA *deal_data, /* Structure of deal structure data  */
                  PREPAY_DATA *prepay_data, /* Structure of prepayment data */
                  TREE_DATA *tree_data,     /* Structure of tree data */
                  HISTORY *history,         // historical index values
                  MBSDENSITY *density) // density function for prepay propensity
{
  int i;
  long mo, dd, yy, mm, ValueDate; /* Value date */
  char str[100];
  double prev;
  char *err;

  if (deal_data->casenumber < 0 || deal_data->casenumber > 13)
    return (nrerror("Invalid casenumber!"));

  if (Dateok(term_data->Today))
    return (nrerror("Incorrect format for today's date!"));

  ValueDate = Nxtbusday(
      term_data
          ->Today, /* Value date is today shifted by (Spotdays) business days */
      term_data->Spotdays, term_data->Holidays, term_data->numHolidays);

  if ((term_data->MMB != 360) && (term_data->MMB != 365))
    return (nrerror("Money Market basis should be 360 or 365!"));

  if ((term_data->AoS != 'A') && (term_data->AoS != 'S'))
    return (
        nrerror("Specify 'A' for Annual yield curve or 'S' for semi-annual!"));

  if ((!strcmp(term_data->YearBase, "ACT")) &&
      (!strcmp(term_data->YearBase, "365")))
    return (
        nrerror("Specify 'ACT' for Actual or '365' for 365 Fixed year basis!"));

  for (i = 0; i < 6; i++)
    if ((term_data->MMYield[i] < .0001) || (term_data->MMYield[i] > 100.))
      return (nrerror("6 Money Market yields are needed!"));

  for (i = 0; i < 8; i++)
    if ((term_data->SwapYield[i] < .0001) || (term_data->SwapYield[i] > 100.))
      return (nrerror("8 swap yields are needed!"));

  if (Dateok(deal_data->Maturity))
    return (nrerror("Incorrect format for maturity's date!"));

  if (deal_data->Maturity <= ValueDate)
    return (nrerror("Maturity falls before value date!"));

  if ((term_data->IndexFreq != 'A') && (term_data->IndexFreq != 'S'))
    return (nrerror("Specify a frequency for the Index: 'A' or 'S'!"));

  if ((deal_data->Freq != 'A') && (deal_data->Freq != 'S') &&
      (deal_data->Freq != 'Q') && (deal_data->Freq != 'M'))
    return (nrerror("Specify a frequency for the underlying cash flows: 'A'  , "
                    "'S'  , 'Q' or 'M'!"));

  if (term_data->IndexMat + term_data->SwaptionExp[NBSWAPTION - 1] > MAXMAT)
    return (nrerror("Not enough zeros to calibrate spot volatility!"));

  if (tree_data->Ppy <= 0)
    return (nrerror("tree_data.Ppy has to be positive integer"));
  if (deal_data->Ppform != 'R' && deal_data->Ppform != 'A')
    return (nrerror("deal_data.Ppform must be either R or A!"));
  if (deal_data->AlorCPR != 'A' && deal_data->AlorCPR != 'C')
    return (nrerror("deal_data.AlorCPR must be either A or C!"));

  prev = 0.0;
  for (i = 0; i < NBPOINTS; ++i) {
    if (deal_data->IndexLevel[i] <= prev)
      return (nrerror("deal_data.IndexLevel must be increasing!"));
    prev = deal_data->IndexLevel[i];
  }

  if (deal_data->Hedge != 'Y' && deal_data->Hedge != 'N' &&
      deal_data->Hedge != 'I' && deal_data->Hedge != 'P')
    return (nrerror("deal_data.Hedge has to be Y or N!"));
  if (deal_data->Schedule != 'Y' && deal_data->Schedule != 'N')
    return (nrerror("deal_data.Schedule has to be Y or N!"));
  if (deal_data->Balloon != 'Y' && deal_data->Balloon != 'N')
    return (nrerror("deal_data.Balloon has to be Y or N!"));
  if (deal_data->Always != 'Y' && deal_data->Always != 'N')
    return (nrerror("deal_data.Always has to be Y or N!"));
  if (deal_data->invfloater != 'Y' && deal_data->invfloater != 'N')
    return (nrerror("deal_data.invfloater has to be Y or N!"));
  if (deal_data->BalloonSch < 0)
    return (nrerror("deal_data.BalloonSch must be non-negative!"));
  if (deal_data->stripratio < 0.0)
    return (nrerror("deal_data.stripratio must be non-negative!"));
  if (deal_data->implvol < 0.0)
    return (nrerror("deal_data.impvol must be non-negative!"));
  //
  if (deal_data->SFSNb <= 0)
    return (nrerror("deal_data.SFSNb must be positive integer!"));
  if (deal_data->Amort[0] < 0.0 || deal_data->Amort[0] > 1.0)
    return (nrerror("deal_data.Amort[0] must be between 0 and 1!"));
  if (deal_data->Amort[deal_data->SFSNb] < 0.0 ||
      deal_data->Amort[deal_data->SFSNb] > 1.0)
    return (nrerror("deal_data.Amort[SFSNb] must be between 0 and 1!"));
  if (deal_data->Nbweights <= 0)
    return (nrerror("deal_data.Nbweights must be positive integer!"));

  prev = 0.0;
  for (i = 0; i < 3; ++i) {
    if (deal_data->SeasoningLevel[i] < prev)
      return (nrerror("deal_data.SeasoningLevel must be increasing!"));
    prev = deal_data->SeasoningLevel[i];
    if (deal_data->SeasoningMat[i] <= 0)
      return (nrerror("deal_data.SeasoningMat must be positive integers!"));
  }
  for (i = 0; i < 12; ++i)
    if (deal_data->Seasonality[i] < 0.0)
      return (nrerror("deal_data.Seasonality must be non-negative!"));

  if (deal_data->Delay < 0)
    return (nrerror("deal_data.Delay must be non-negative integer!"));
  if (deal_data->NbFixedPrepay < 0)
    return (nrerror("deal_data.NbFixedPrepay must be non-negtaive integer!"));

  //#ifdef CORRECTDELAY
  if (deal_data->Delay && deal_data->NbFixedPrepay < deal_data->Delay)
    return (nrerror("NbFixedPrepay < Delay!"));
  //#endif

  for (i = 0; i < deal_data->NbFixedPrepay; ++i)
    if (deal_data->FixedPrepay[i] < 0.0 || deal_data->FixedPrepay[i] > 1.0)
      return (nrerror("deal_data.FixedPrepay must be between 0 and 1!"));

  if (density->nirates <= 0)
    return (nrerror("density.nirates must be positive integer!"));
  if (err = MBSDENSITY_Check(density))
    return (err);

  if (err = HISTORY_Check(history))
    return (err);

  mo = Nxtmth(term_data->Today, -1, 1);
  Dsplit(mo, &yy, &mm, &dd);
  mo = DateYMMDD(yy, mm, 1);

  if (YYYYMM2FirstOfMon(history->dates[history->Number - 1]) < mo)
    return (nrerror("Incomplete history"));

  //#ifdef USESRTDATES
  //	for(i=1; i < term_data->numHolidays; ++i)
  //		if(Dateok(toExcelDate(term_data->Holidays[i])) !=0)
  //		{
  //			sprintf( str  , "Invalid holiday->date[%d]=%ld"  , i  ,
  //term_data->Holidays[i] ); 			return(nrerror(str));
  //		}
  //#else
  //	for(i=0; i < term_data->numHolidays; ++i)
  //		if(Dateok(term_data->Holidays[i]) !=0)
  //		{
  //			sprintf( str  , "Invalid holiday->date[%d]=%ld"  , i  ,
  //term_data->Holidays[i] ); 			return(nrerror(str));
  //		}
  //#endif

  return (0);

} /* Input_Check */

/*****  Conv_Freq  ******
 *       Convert the character Freq ('A'nnual  , 'S'emi-Annual  , 'Q'uarterly or
 *       'M'onthly) into an integer (1  , 2  , 4 or 12).
 */
int Conv_Freq(char Freq) {
  switch (Freq) {
  case 'A':
    return (1);
  case 'S':
    return (2);
  case 'Q':
    return (4);
  case 'M':
    return (12);
  default:
    return (0);
  } /* switch */
} /* Conv_Freq */

/*****  Conv_AoS  ***********************************************************/
/*
 *       Convert the character AoS into a double (1. or 2.).
 */
double Conv_AoS(char AoS) {
  if (AoS == 'A')
    return (1.);
  else
    return (2.);
} /* Conv_AoS */

/*******************************************************
 *Dates Related Utilities:
 ********************************************************/

/***Holiday***
 *This routine checks whether the date passed in is a holiday
 *or not according to the HOLIDAYS file. It returns (0) if
 *it is not and (-1) if it is. The format is YYYYMMDD.
 */
int Holiday(long date_i, long *holidays, int numHolidays) {

  int xflag = ERR;
  if (isHoliday(date_i))
    return (ERR);
  else
    return (OKAY);
} /*Holiday*/

long Nxtbusday(long date, int len, long *holidays, int num) {
  return (MBSNxtbusday(date, len));
}

long Daysact(long date1_i, long date2_i) {
  return (DateDiff(date1_i, date2_i, _DACT));
}

/* pksort3
 *  Sort routine of three arrays
 */
void pksort3(int n, double *arr, int *brr, int *crr) {
  int i, j, b, c;
  double a;

  for (j = 1; j < n; j++) {
    a = arr[j];
    b = brr[j];
    c = crr[j];
    i = j - 1;
    while (i >= 0 && arr[i] > a) {
      arr[i + 1] = arr[i];
      brr[i + 1] = brr[i];
      crr[i + 1] = crr[i];
      i--;
    }
    arr[i + 1] = a;
    brr[i + 1] = b;
    crr[i + 1] = c;
  }
} /*pksort3*/

/**error handling***/

/***memory (de-)allocation**/

//***************merror_register*************************
// record error and exit
////////////////////////////////////////////////
void merror_register(char *error_text) {
  mbs_err_handler(error_text, 1, 1);
  /*********
          FILE *outError;
          char string[100];
          string[0]='\0';
          strcat(string  ,OUTPUTDIR);
          strcat(string  ,"ErrorMessages.txt");
          outError=fopen(string  , "wt");

          if(outError)
          {
                  fprintf (outError  ,"run-time error...\n");
                  fprintf (outError  ,"%s\n"  ,error_text);
          }
          fclose(outError);
          exit (1);
          *************/
}

/*****  mbs_ivector
 * ************************************************************/
/*
 *       Allocation of memory for an array of integers starting at nl  , ending
 *       at nh.
 */
int *mbs_ivector(int nl, int nh) {
  return (mbspt_ivector(nl, nh));
  /************
  int
          *v;


  v = (int *) calloc ((unsigned) (nh-nl+1)  , sizeof (int));
  if (!v)
          merror_register("Can't calloc ivector");

  return (v-nl);
***********/
} /* mbs_ivector */

/*****  mbs_dvector
 * ************************************************************/
/*
 *       Allocation of memory for an array of doubles starting at nl  , ending
 *       at nh.
 */
double *mbs_dvector(int nl, int nh) {
  return (mbspt_dvector(nl, nh));
  /********
  double
          *v;


  v = (double *) calloc ((unsigned) (nh-nl+1)  , sizeof (double));
  if (!v)
          merror_register("Can't calloc dvector");

  return (v-nl);
********/
} /* mbs_dvector */

/*****  mbs_dmatrix
 * ************************************************************/
/*
 *       Allocation of memory for a matrix of doubles (nrl  , nrh) x (ncl  ,
 * nch).
 */
double **mbs_dmatrix(int nrl, int nrh, int ncl, int nch) {
  return (mbspt_dmatrix(nrl, nrh, ncl, nch));
  /********
  int
          i;
  double
          **m;


  m = (double **) calloc ((unsigned) (nrh-nrl+1)  , sizeof (double*));
  if (!m)
          merror_register("Can't calloc (double **) for dmatrix");
  m -= nrl;

  for (i=nrl; i<=nrh; i++)
  {
          m[i] = (double *) calloc ((unsigned) (nch-ncl+1)  , sizeof (double));
          if (!m[i])
                  merror_register("Can't calloc *(double *) in dmatrix");
          else
                  m[i] -= ncl;
  }  // for i

  return (m);
**********/
} /* mbs_dmatrix */

/*****  mbs_free_ivector
 * *******************************************************/
/*
 *       Free array of int v (nl x nh).
 */
void mbs_free_ivector(int *v, int nl, int nh) {
  mbspt_free_ivector(v, nl, nh);
  /*******
  free ((char *) (v+nl));
*****/
} /* mbs_free_ivector */

/*****  mbs_free_dvector
 * *******************************************************/
/*
 *       Free array of doubles v (nl x nh).
 */
void mbs_free_dvector(double *v, int nl, int nh) {
  mbspt_free_dvector(v, nl, nh);
  /*********
  free ((char *) (v+nl));
  *********/

} /* mbs_free_dvector */

/*****  mbs_free_dmatrix
 * *******************************************************/
/*
 *       Free matrix of doubles m (nrl  , nrh) x (ncl  , nch).
 */
void mbs_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch) {
  mbspt_free_dmatrix(m, nrl, nrh, ncl, nch);
  /***********
  int
          i;


  for (i=nrh; i>=nrl; i--)
          free ((char *) (m[i] + ncl));

  free ((char *) (m+nrl));
********/
} /* mbs_free_dmatrix */

/*****  Tree_Free  **********************************************************/
/*
 *       Free the arrays constituting the tree.
 */
void Tree_Free(int NbNodes, TREE_DATA *tree_data) {
  mbs_free_ivector(tree_data->ExerFlag, 0, NbNodes);
  mbs_free_dvector(tree_data->FwdStrike, 0, NbNodes);
  mbs_free_dvector(tree_data->Accrued, 0, NbNodes);
  mbs_free_dvector(tree_data->Length, 0, NbNodes);
  mbs_free_dvector(tree_data->PrepayVol, 0, NbNodes);
  mbs_free_dvector(tree_data->FwdPrepay, 0, NbNodes);
  mbs_free_dvector(tree_data->Drift, 0, NbNodes);
  mbs_free_dvector(tree_data->FwdVol, 0, NbNodes);
  mbs_free_dvector(tree_data->FwdRate, 0, NbNodes + 1);
  mbs_free_dvector(tree_data->Delta, 0, NbNodes + 1);
  mbs_free_dvector(tree_data->FwdRateS, 0, NbNodes + 1);
  mbs_free_dvector(tree_data->ZeroCoupon, 0, NbNodes + 1);
  mbs_free_dvector(tree_data->ZeroCouponS, 0, NbNodes + 1);

  return;

} /* Tree_Free */

/*misc*/

/* Given IO and PO prices(I and P)  , figure out the required speedup and
 * rationality (ss and rr) */
/* inputs: io[3]  , po[3]  , rtn[3]  , spd[3]  , I  , P */
/* which are 3 io  , po prices  , at 3 rtn  ,spd (rationality and speed) values
 * that form a plane */
/* we do interpolation to find rr and ss at I  , P */
void rs(double *io, double *po, double *rtn, double *spd, double I, double P,
        double *rr, double *ss, double *IO_spd, double *IO_rtn, double *PO_spd,
        double *PO_rtn) {
  double dd, a1, b1, a2, b2, dr1, dr2, ds1, ds2, dio1, dio2, dpo1, dpo2, dio,
      dpo, temp;

  dio = I - io[0]; /* difference of io  ,po  ,rtn  ,spd  , */
  dio1 = io[1] - io[0];
  dio2 = io[2] - io[0];
  dpo = P - po[0];
  dpo1 = po[1] - po[0];
  dpo2 = po[2] - po[0];
  dr1 = rtn[1] - rtn[0];
  dr2 = rtn[2] - rtn[0];
  ds1 = spd[1] - spd[0];
  ds2 = spd[2] - spd[0];

  dd = dr2 * ds1 - dr1 * ds2;          /* determinant */
  a1 = (dr2 * dio1 - dr1 * dio2) / dd; /* calculate gradients*/
  b1 = (ds1 * dio2 - ds2 * dio1) / dd;
  a2 = (dr2 * dpo1 - dr1 * dpo2) / dd;
  b2 = (ds1 * dpo2 - ds2 * dpo1) / dd;

  dd = b2 * a1 - b1 * a2; /* linear interpolation */
  *rr = rtn[0] - (a2 * dio - a1 * dpo) / dd;
  *ss = spd[0] + (b2 * dio - b1 * dpo) / dd;
  *IO_spd = a1;
  *IO_rtn = b1;
  *PO_spd = a2;
  *PO_rtn = b2; /* output gradients */

  if (*rr < 0.0) { /* if negative  , return min */
    temp = (rtn[1] < rtn[2]) ? rtn[1] : rtn[2];
    *rr = (temp < rtn[0]) ? temp : rtn[0];
  }

  if (*ss < 0.0) {
    temp = (spd[1] < spd[2]) ? spd[1] : spd[2];
    *ss = (temp < spd[0]) ? temp : spd[0];
  }
} /*rs*/

/*  convert smm   with */
/* no scheduled amortization into  average life              */
double smmal1(double smm,
              int N) /* smm is SMM  , N is remaining term in months */
{

  double al = 0.; /* al is average life in years */

  al = (1.0 - pow((1.0 - smm), (double)N)) / smm;

  return al / 12.;
} /*smmal1*/
