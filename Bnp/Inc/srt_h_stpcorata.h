/* -------------------------------------------------------------------------

        AUTHOR		: O Van Eyseren
        DATE		: Feb 23 1995
        FILE NAME	: srt_h_stpcorata.h
        PURPOSE		: function to attache the right correlation/coefficients
                                  matrix to the stps used (coeff for MC  ,
   correl for the tree)
   ------------------------------------------------------------------------- */
#ifndef STPCORATA_H
#define STPCORATA_H

Err srt_f_attach_correl_to_stp(SrtStpPtr stp, SrtCorrLst *cls);

#endif
