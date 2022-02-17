#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mint      lv_ipage1[2] = {78, 60};

/*Translated by FOR_C++, v0.1, on 06/08/90 at 14:18:53 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/08/90 at 14:18:51
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    Purpose:    Set or retrieve page width and length for printing.

    Usage:      CALL imsl_page (IOPT, IPAGE)

    Arguments:
       IOPT   - Page attribute option.  (Input)
                IOPT   Description of Attribute
                -1, 1  Page width.
                -2, 2  Page length.
                Negative values of IOPT indicate the setting IPAGE is
                input.  Positive values of IOPT indicate the setting
                IPAGE is output.
       IPAGE  - Value of page attribute.  (Input, if IOPT is negative;
                output, if IOPT is positive.)
                IOPT    Description of Attribute   Settings for IPAGE
                -1, 1   Page width (in characters) 10, 11, ...
                -2, 2   Page length (in lines)     10, 11, ...
  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_page(Mint iopt, Mint *ipage)
#else
void imsl_page(iopt, ipage)
	Mint             iopt, *ipage;
#endif
{

	imsl_e1psh("imsl_page"); 
        if (iopt < -2 || iopt > 2 || iopt == 0) 
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_PAGE_OPTION);
	if (iopt < 0) {
		if (iopt == -1) {
			if (*ipage < 10) {
				imsl_e1sti(1, *ipage);
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_BAD_PAGE_WIDTH);
			}
		} else if (iopt == -2) {
			if (*ipage < 10) {
				imsl_e1sti(1, *ipage);
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_BAD_PAGE_LENGTH);
			}
		}
	}
	if (!imsl_n1rty(0)){               /* Retrieve page characteristic. */
		if (iopt > 0) {
			*ipage = lv_ipage1[iopt - 1];
		} else {                /* Set page characteristic. */
			lv_ipage1[-iopt - 1] = *ipage;
		}
	}
	imsl_e1pop("imsl_page");
        return;
}

