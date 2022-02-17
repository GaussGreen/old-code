#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*-----------------------------------------------------------------------
    Purpose:    Set or retrieve an option for printing a matrix.

    Usage:      CALL imsl_write_options (IOPT, ISET, ISCOPE)

    Arguments:
       IOPT   - Indicator of option type.  (Input)
                IOPT    Description of Option Type
                -1, 1   Horizontal centering or left justification of
                        matrix to be printed
                -2, 2   Method for printing large matrices
                -3, 3   Paging
                -4, 4   Method for printing NaN (not a number)
                -5, 5   Title option
                -6, 6   Default format for real and complex numbers
                  0     Reset all the current settings saved in internal
                        variables back to their last setting made with an
                        invocation of imsl_write_options with ISCOPE = 1.  (This
                        option is used internally by routines printing a
                        matrix, and is not useful otherwise.)
                If IOPT is negative, ISET and ISCOPE are input and
                are saved in internal variables.  If IOPT is positive,
                ISET and ISCOPE are output and receive the values of the
                internal variables.  If IOPT = 0, ISET and ISCOPE are
                not referenced.
       ISET   - Setting for option selected by IOPT.  (Input, if IOPT is
                negative; output, if IOPT is positive; not referenced if
                IOPT = 0)
                IOPT    ISET  Meaning
                -1, 1     0   Matrix is left justified.
                          1   Matrix is centered horizontally on page.
                -2, 2     0   A complete row is printed before the next
                              row is printed.  Wrapping is used if
                              necessary.
                          m   Here, m is a positive integer.  Let n1 be
                              the maximum number of columns beginning
                              with column 1 that fit across the page (as
                              determined by the widths of the printing
                              formats).  First, columns 1 through n1 are
                              printed for rows 1 through m.  Let n2 be
                              the maximum number of columns beginning
                              with column n1 + 1 that fit across the
                              page.  Second, columns n1 + 1 through n1 +
                              n2 are printed for rows 1 through m.  This
                              continues until the last columns are
                              printed for rows 1 through m.  Printing
                              continues in this fashion for the next m
                              rows, etc.
                 -3, 3   -2   Printing begins on the next line and no
                              paging occurs.
                         -1   Paging is on.  Every invocation of a WR***
                              routine begins on a new page and paging
                              occurs within each invocation as is needed.
                          0   Paging is on.  The first invocation of a
                              WR*** routine begins on a new page.
                              Subsequent paging occurs as is needed.
                              With this option, every invocation of a
                              WR*** routine ends with a call to imsl_write_options to
                              reset this option to k, a positive integer
                              giving the number of lines printed on the
                              current page.
                          k   Here, k is a positive integer.  Paging is
                              on, and k lines have been printed on the
                              current page.  If k is less than the page
                              length IPAGE (see IMSL routine imsl_page), then
                              IPAGE - k lines are printed before a new
                              page instruction is issued.  If k is
                              greater than or equal to IPAGE, then the
                              first invocation of a WR*** routine begins
                              on a new page.  In any case, subsequent
                              paging occurs as is needed.  With this
                              option, every invocation of a WR*** routine
                              ends with a call to imsl_write_options to reset the
                              value of k.
                 -4, 4    0   NaN is printed as '..........'
                          1   NaN is printed as '          '
                 -5, 5    0   Title appears only on first page.
                          1   Title appears on the first page and all
                              continuation pages.
                 -6, 6    0   Format is '(W10.4)'.  See Remark 2.
                          1   Format is '(W12.6)'.  See Remark 2.
                          2   Format is '(1PE12.5)'.
       ISCOPE - Indicator of the scope of the option.  (Input if
                IOPT is nonzero; not referenced if IOPT = 0)
                ISCOPE  Action
                  0     Setting is temporarily active for the next
                        invocation of a WR*** matrix printing routine.
                  1     Setting is active until it is changed by
                        another invocation of imsl_write_options.

    Remarks:
    1. This subprogram can be invoked repeatedly before using a WR***
       routine to print a matrix.   The matrix printing routines retrieve
       these settings to determine the printing options.  It is not
       necessary to call imsl_write_options if a default value of a printing option
       is desired.  The defaults are as follows.
          IOPT  Default Value for ISET
           1       0
           2    1000
           3      -2
           4       0
           5       0
           6       0

    2. The W format is a special format that can be used to automatically
       select a pretty FORTRAN format--either E, F, or I.  The format is
       specified as 'Wn.d'.  Here, n is the field width, and d is the
       number of significant digits generally printed.
  -----------------------------------------------------------------------
 */

void imsl_write_options(iopt, isetx)
	Mint             iopt, *isetx;
{
	Mint             i;

	imsl_e1psh("imsl_write_options");  
        if (iopt < -6 || iopt > 6)
            imsl_ermes(IMSL_TERMINAL, IMSL_ILLEGAL_WRITE_OPTION);
        
	if (iopt == -1) {
            if (*isetx <0 || *isetx > 1) {
                imsl_e1stl(1,"IMSL_SET_CENTERING");
                imsl_e1sti(1,0);
                imsl_e1sti(2,1);
                imsl_e1sti(3,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPTION_VALUE);
            }
	} else if (iopt == -2) {
            if (*isetx <0) {
                imsl_e1stl(1,"IMSL_SET_ROW_WRAP");
                imsl_e1sti(1,0);
                imsl_e1sti(2,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPT_VALUE);
            }
	} else if (iopt == -3) {
            if (*isetx < -2) {
                imsl_e1stl(1,"IMSL_SET_PAGING");
                imsl_e1sti(1,-2);
                imsl_e1sti(2,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPT_VALUE);
            }
	} else if (iopt == -4) {
            if (*isetx <0 || *isetx > 1) {
                imsl_e1stl(1,"IMSL_SET_NAN_CHAR");
                imsl_e1sti(1,0);
                imsl_e1sti(2,1);
                imsl_e1sti(3,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPTION_VALUE);
            }
	} else if (iopt == -5) {
            if (*isetx <0 || *isetx > 1) {
                imsl_e1stl(1,"IMSL_SET_TITLE_PAGE");
                imsl_e1sti(1,0);
                imsl_e1sti(2,1);
                imsl_e1sti(3,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPTION_VALUE);
            }
	} else if (iopt == -6) {
            if (*isetx <0 || *isetx > 2) {
                imsl_e1stl(1,"IMSL_SET_FORMAT");
                imsl_e1sti(1,0);
                imsl_e1sti(2,2);
                imsl_e1sti(3,*isetx);
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_OPTION_VALUE);
            }
        }
	if (!imsl_n1rty(0))  imsl_w1opt(iopt,isetx);
	imsl_e1pop("imsl_write_options");
	return;
}

	
void imsl_w1opt(iopt, isetx)
	Mint             iopt, *isetx;
{
	Mint i;
	static Mint      isety[6] = {0, 1000, -2, 0, 0, 0};

	if (iopt > 0) {         /* Retrieve printing option */
		*isetx = isety[iopt - 1];
	} else if (iopt < 0) {  /* Set printing option. */
		isety[-iopt - 1] = *isetx;
        } else {
                isety[0] = 0;
                isety[1] = 1000;
                isety[2] = -2;
                isety[3] = 0;
                isety[4] = 0;
                isety[5] = 0;
        }
	return;
}
