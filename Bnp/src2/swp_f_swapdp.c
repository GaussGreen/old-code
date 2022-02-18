/* =========================================================================

   FILENAME: swp_f_swapdp.c

   FUNCTION: swp_f_initSwapDP

   PURPOSE:  to initialise a SwapDP from start, end_nfp, compString, basisStr
   ========================================================================= */

#include "math.h"
#include "swp_h_all.h"

Err swp_f_setSwapDP(long start, long end_nfp, SrtCompounding comp, BasisCode basis, SwapDP* sdp)
{
    Err         err = NULL;
    SrtDateList datelist;

    if (err = srt_f_gen_test_date(start))
    {
        return err;
    }

    sdp->start             = start;
    sdp->first_full_fixing = sdp->start;
    sdp->nfp               = end_nfp;
    sdp->end               = end_nfp;

    if ((end_nfp < 30317) && (end_nfp > 100))
    {
        return serror("End date %d too long ago", sdp->end);
    }

    if ((0 < end_nfp) && (end_nfp < 100))
    {
        sdp->direction = FWD;
    }
    else
    {
        sdp->direction = BKWD;
    }

    if (sdp->direction == BKWD && sdp->start > sdp->end)
    {
        return serror("Start %d must be before end: %d", sdp->start, sdp->end);
    }

    sdp->compd      = comp;
    sdp->basis_code = basis;

    /* Now, sets a proper theoretical EndDate if direction = FWD */
    if (sdp->end < 100)
    {
        datelist = SwapDP_to_DateList(sdp, NO_BUSDAY_CONVENTION);
        sdp->end = datelist.date[datelist.len - 1];
        swp_f_free_in_DateList(datelist);
    }

    return NULL;

} /* END Err swp_f_setSwapDP(...) */

/* ------------------------------------------------------------------------------ */

Err swp_f_initSwapDP(long start, long end_nfp, String compStr, String basisStr, SwapDP* sdp)
{
    Err err = NULL;

    SrtCompounding comp;
    BasisCode      basis;

    /* Interpret compounding */
    if (err = interp_compounding(compStr, &comp))
    {
        return err;
    }

    /* Interpret basis */
    if (err = interp_basis(basisStr, &basis))
    {
        return err;
    }

    /* Sets the elements in the structure */
    err = swp_f_setSwapDP(start, end_nfp, comp, basis, sdp);

    return err;

} /* END Err swp_f_initSwapDP(...) */
