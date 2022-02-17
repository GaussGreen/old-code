/* ===================================================================================
   FILENAME:      BreakEven.cxx

   PURPOSE:       Computes Break Even with sliding vol assumption
   ===================================================================================
 */

#pragma warning(disable : 4786) // Disable long name warnings

#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"
//#include "swp_h_amortswaption.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_vol.h"

Err BreakEvenFromMkt(char *yc_name, char *vc_name, long spot_date, int spot_lag,
                     char *swap_freq, char *swap_basis, char *ref_rate_name,
                     int NMaturity, char **maturity_tenor, int NUnderlying,
                     char **underlying_tenor, int NbOfBDayPerYear,
                     double **BEMatrix) {
  Err err = NULL;
  int i, j;

  double power;

  double swap1, vol1;
  long start1, start1_lag, end1;

  double swap2, vol2;
  long start2, end2;

  double BE2;

  long today;

  today = add_unit(spot_date, -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

  for (j = 0; j < NUnderlying; ++j) {
    err = add_tenor(spot_date, maturity_tenor[0], MODIFIED_SUCCEEDING, &start1);
    if (err) {
      return err;
    }

    err = add_tenor(start1, underlying_tenor[j], NO_BUSDAY_CONVENTION, &end1);
    if (err) {
      return err;
    }

    err = swp_f_ForwardRate(start1, end1, swap_freq, swap_basis, yc_name,
                            ref_rate_name, &swap1);
    if (err) {
      return err;
    }

    err = swp_f_vol(vc_name, start1, end1, swap1, &vol1, &power);
    if (err) {
      return err;
    }

    BEMatrix[0][j] = vol1 / sqrt(NbOfBDayPerYear);
  }

  for (i = 1; i < NMaturity; ++i) {
    for (j = 0; j < NUnderlying; ++j) {
      err =
          add_tenor(spot_date, maturity_tenor[i], MODIFIED_SUCCEEDING, &start1);
      if (err) {
        return err;
      }

      err = add_tenor(start1, underlying_tenor[j], NO_BUSDAY_CONVENTION, &end1);
      if (err) {
        return err;
      }

      err = swp_f_ForwardRate(start1, end1, swap_freq, swap_basis, yc_name,
                              ref_rate_name, &swap1);
      if (err) {
        return err;
      }

      err = swp_f_vol(vc_name, start1, end1, swap1, &vol1, &power);
      if (err) {
        return err;
      }

      err = add_tenor(spot_date, maturity_tenor[i - 1], MODIFIED_SUCCEEDING,
                      &start2);
      if (err) {
        return err;
      }

      err = add_tenor(start2, underlying_tenor[j], NO_BUSDAY_CONVENTION, &end2);
      if (err) {
        return err;
      }

      err = swp_f_ForwardRate(start2, end2, swap_freq, swap_basis, yc_name,
                              ref_rate_name, &swap2);
      if (err) {
        return err;
      }

      err = swp_f_vol(vc_name, start2, end2, swap2, &vol2, &power);
      if (err) {
        return err;
      }

      start1_lag = add_unit(start1, -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

      BE2 = 1 + (start1_lag - today) * (log(vol1 * vol1) - log(vol2 * vol2)) /
                    (start1 - start2);

      if (BE2 > 0) {
        BEMatrix[i][j] = vol1 * sqrt(BE2) / sqrt(NbOfBDayPerYear);
      } else {
        smessage("Negative Break Even for maturity %d and underlying %d ",
                 i + 1, j + 1);
        BEMatrix[i][j] = 0;
      }
    }
  }

  return err;
}

Err BreakEvenFromVol(long spot_date, int spot_lag, int NMaturity,
                     char **maturity_tenor, int NUnderlying,
                     char **underlying_tenor, int NbOfBDayPerYear,
                     double **VolMatrix, double **BEMatrix) {
  Err err = NULL;
  int i, j;

  double vol1;
  long start1, start1_lag;

  double vol2;
  long start2;

  double BE2;

  long today;

  today = add_unit(spot_date, -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

  for (j = 0; j < NUnderlying; ++j) {
    BEMatrix[0][j] = VolMatrix[0][j] / sqrt(NbOfBDayPerYear);
  }

  for (i = 1; i < NMaturity; ++i) {
    for (j = 0; j < NUnderlying; ++j) {
      err =
          add_tenor(spot_date, maturity_tenor[i], MODIFIED_SUCCEEDING, &start1);
      if (err) {
        return err;
      }

      err = add_tenor(spot_date, maturity_tenor[i - 1], MODIFIED_SUCCEEDING,
                      &start2);
      if (err) {
        return err;
      }

      vol1 = VolMatrix[i][j];
      vol2 = VolMatrix[i - 1][j];

      start1_lag = add_unit(start1, -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

      BE2 = 1 + (start1_lag - today) * (log(vol1 * vol1) - log(vol2 * vol2)) /
                    (start1 - start2);

      if (BE2 > 0) {
        BEMatrix[i][j] = vol1 * sqrt(BE2) / sqrt(NbOfBDayPerYear);
      } else {
        smessage("Negative Break Even for maturity %d and underlying %d ",
                 i + 1, j + 1);
        BEMatrix[i][j] = 0;
      }
    }
  }

  return err;
}
