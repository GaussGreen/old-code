#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define fmax(a,b) ((a>b)?a:b)

#ifdef ANSI
        static void l_w3trl(Mchar*, Mint*, Mint*);
#else
        static void l_w3trl();
#endif

/* -----------------------------------------------------------------------
    Purpose:    Write a title messsage.

    Usage:      CALL WRTRL (TITLE, LINLEN)

    Arguments:
       TITLE  - Character string specifying the title. (Input)
                TITLE = ' ' suppresses printing of the title.
       LINLEN - Length of a output line.  (Input)

    Remarks:
    1. Use '%/' within titles to cause a new line to occur.
       Long titles is automatically wrapped.

    2. Output is written to the unit specified by subroutine UMACH.

    3. Maximum length of a long title is 256 characters.
  -----------------------------------------------------------------------
 */
void imsl_wrtrl(title, linlen)
	Mchar           *title;
	Mint             linlen;
{
	Mchar            *ipt1, *ipt;
	Mint             done;
	Mint             _l0, ilen, ltitle;

	imsl_e1psh("WRTRL");
	/*
	 * If title is a blank title, print nothing.
	 */
	if (strlen(title) == 0) goto L_9000;

	if (linlen <= 0) {
		    /* The length for an output line must be greater */
		    /* than zero. */
		imsl_ermes(IMSL_TERMINAL, IMSL_ZERO_LINE);
		goto L_9000;
	}
	ltitle = strlen(title);
	if (ltitle > 255) {
		    /* The title message length is too long.  Maximum */
		    /* length of a title is 255 characters.	      */
		imsl_ermes(IMSL_TERMINAL, IMSL_LONG_TITLE);
		goto L_9000;
	}
	                                    /* Print multiple lines. */
	done = 0;
	ipt = title;
	                                    /* dountil */
        do {
            ipt1 = strchr(ipt,(Mint)'\n');
	    if (ipt1 == title) {
                    _l0 = -1;
		    l_w3trl(ipt1, &linlen, &_l0);
		    ipt += 1;
		    if (ltitle==1) done = 1;
	    } else if (ipt1 == 0 && ipt >= title+ltitle) {done = 1;
	    } else if (ipt1 == 0) {
		    ilen = ltitle - (Mint)(ipt-title);
		    l_w3trl(ipt, &linlen, &ilen);
		    done = 1;
	    } else {
		    ilen = (Mint) (ipt1-ipt) - 1;
		    l_w3trl(ipt, &linlen, &ilen);
		    ipt = ipt1 + 1;
	    }
        } while(!done);

L_9000: imsl_e1pop("WRTRL");
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Write a title message.

    Usage:      CALL W3TRL (TITLE, LINLEN, ILEN)

    Arguments:  See WRTRL.
  -----------------------------------------------------------------------
 */
static void l_w3trl(title, linlen, ilen)
	Mchar           *title;
	Mint            *linlen, *ilen;
{
	Mint             _l0, il, iln, ipgopt;
        FILE             *nout;

	imsl_e1psh("W3TRL");
	/*
	 * Get the number of line printed in the current page_c.
	 */
	if (*ilen == -1) {
		imsl_write_options(3, &ipgopt);
                imsl_umach(2,&nout);
		fprintf(nout, "\n");
                _l0 = ipgopt+1;
		imsl_write_options(-3, &_l0);
		                    /* Print multiple lines if needed. */
	} else if (*ilen >= 1) {
		il = 0;
		                    /* dowhile (il.LE.ilen) */
                do {
		    if (!(il < *ilen)) break;
		    iln = *ilen - il;
		    if (iln <= *linlen) {
			imsl_write_line(iln, title + il);
			il += iln;
		    } else {
			imsl_write_line(*linlen, title+il);
			il += *linlen;
		    }
                } while (il <= *ilen);
	}
	imsl_e1pop("W3TRL");
	return;
}


