/*Translated by FOR_C++, v0.1, on 06/12/90 at 12:58:53 */ /*FOR_C++
Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */

#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef IMSL_MACHINE_NT
#include <windows.h>
#endif

#define ITRUE 1
#define IFALSE 0
#define NSPACE 2
#define	NDENTR	3

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)

#define MAT(I_,J_)	(mat+(I_)*(ldmat)+(J_))


Mint            PROTO(imsl_c1tci_f, (Mchar *, Mint, Mint *));
static Mchar*   PROTO(l_i1tri,(Mint,Mint,Mint,Mint,Mchar*,Mint));
static void     PROTO(l_w3irl_f,(Mchar*,Mint*,Mint,Mint,Mchar*,Mchar**,Mchar**,
	                Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,
                        Mint,Mint,Mint));
static void     PROTO(l_w10rl,(Mint,Mint,Mchar**,Mchar**,Mint,Mint,Mint,Mint,
                        Mint,Mint*,Mint*));
static void     PROTO(l_w11rl,(Mint,Mint,Mint,Mint,Mint,Mint*,Mint*,Mint*,
                        Mint*,Mint*));
static void     PROTO(l_w4rrl,(Mchar*,Mint,Mint,Mint,Mint*,Mint*,Mint*,Mint*,
                        Mint*));
static VA_LIST_HACK  PROTO(l_i_write_matrix,(Mchar*,Mint,Mint,Mint[],va_list));

Mchar		*imsl_output_string;
Mint		imsl_output_string_length;
int		imsl_return_string, imsl_write_to_console;


#ifdef ANSI
void imsl_i_write_matrix(Mchar *title, Mint nra, Mint nca, Mint a[], ...)
#else
void imsl_i_write_matrix(title, nra, nca, a, va_alist)
    Mchar        *title;
    Mint         nra, nca, a[];
    va_dcl
#endif
{
    va_list     argptr;

    VA_START (argptr, a);
    imsl_e1psh("imsl_i_write_matrix");

    IMSL_CALL(l_i_write_matrix(title, nra, nca, a, argptr));

    imsl_e1pop("imsl_i_write_matrix");
    va_end (argptr);
}

#ifdef ANSI
static VA_LIST_HACK l_i_write_matrix(Mchar *title, Mint nra, Mint nca, Mint a[],
                             va_list argptr)
#else
static VA_LIST_HACK l_i_write_matrix(title, nra, nca, a, argptr)
    Mchar        *title;
    Mint         nra, nca, *a;
    va_list     argptr;
#endif

{
    Mchar        local[12];
    Mint         arg_number  = 4, iclab= -1, irlab= -1;
    Mint         itring = 0, itri = 0, irl = 0, icl = 0, code, ner=1;
    Mint         cda = ((nca==0)?1:nca);
    Mchar        **rlabel=0;
    Mchar        **clabel=0;
    Mchar        *fmt=0;
    Mint         i = 1, m, j, i1, nrow, mij, ndig, ipgopt;
    Mint         transpose = 1;
    Mchar        **user_output_string=0;

    imsl_return_string = 0;
    imsl_write_to_console = 0;
    while (i) {
        arg_number++;
        switch(code=va_arg(argptr, Mint)) {
            case IMSL_PRINT_ALL:
		if (!itri) {
                   itri = 1;
                   itring = 0;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_PRINT_OPTIONS);
                   goto END_OPTIONS;
                }
            case IMSL_PRINT_UPPER:
		if (!itri) {
                   itri = 1;
                   itring = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_PRINT_OPTIONS);
                   goto END_OPTIONS;
                }
            case IMSL_PRINT_LOWER:
		if (!itri) {
                   itri = 1;
                   itring = -1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_PRINT_OPTIONS);
                   goto END_OPTIONS;
                }
            case IMSL_PRINT_UPPER_NO_DIAG:
		if (!itri) {
                   itri = 1;
                   itring = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_PRINT_OPTIONS);
                   goto END_OPTIONS;
                }
            case IMSL_PRINT_LOWER_NO_DIAG:
		if (!itri) {
                   itri = 1;
                   itring = -2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_PRINT_OPTIONS);
                   goto END_OPTIONS;
                }
            case IMSL_A_COL_DIM:
		arg_number++;
                cda = va_arg(argptr, Mint);
                break;
            case IMSL_ROW_LABELS:
		if (!irl) {
                   irl = 1;
                   arg_number++;
                   rlabel = va_arg(argptr, Mchar**);
                   irlab = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_NO_ROW_LABELS:
		if (!irl) {
                   irl = 1;
                   irlab = 0;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_ROW_NUMBER_ZERO:
		if (!irl) {
                   irl = 1;
                   irlab = 3;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_ROW_NUMBER:
		if (!irl) {
                   irl = 1;
                   irlab = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_COL_LABELS:
		if (!icl) {
                   icl = 1;
                   arg_number++;
                   clabel = va_arg(argptr, Mchar**);
                   iclab = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COL_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_NO_COL_LABELS:
		if (!icl) {
                   icl = 1;
                   iclab = 0;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COL_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_COL_NUMBER_ZERO:
		if (!icl) {
                   icl = 1;
                   iclab = 3;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COL_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_COL_NUMBER:
		if (!icl) {
                   icl = 1;
                   iclab = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_COL_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_WRITE_FORMAT:
		arg_number++;
                fmt = va_arg(argptr, Mchar*);
                imsl_null_pointer("fmt", -1, (void *)fmt);
                break;
	    case IMSL_TRANSPOSE:
                transpose = 0;
                break;
#ifdef IMSL_MACHINE_NT
	    case IMSL_RETURN_STRING:
		imsl_return_string = 1;
		imsl_output_string = NULL;
		imsl_output_string_length = 0;
		user_output_string = va_arg(argptr, Mchar**);
		arg_number++;
		break;
	    case IMSL_WRITE_TO_CONSOLE:
		imsl_write_to_console = 1;
		break;
#endif
	    case 0:
		goto END_OPTIONS;
	    default:
                                    /* Argument number %(I2) is an unknown */
                                    /* optional argument %(I1). */
	    	imsl_e1sti (1,code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    	break;
	}
    }
    END_OPTIONS:
#ifdef IMSL_MACHINE_NT
    if (imsl_return_string && imsl_write_to_console) {
	imsl_ermes (IMSL_TERMINAL, IMSL_RETURN_STRING_ONLY);
	goto RETURN;
    } else if (imsl_write_to_console) {
	HANDLE	console_handle = GetStdHandle(STD_OUTPUT_HANDLE);

	/* If no console is available, create a default one. */
	if (!console_handle)
	    imsl_create_console();
    }
#endif
				/* Check FMT for correct form */

    imsl_c1iarg(nra, "nra", 0, -1, &ner);
    imsl_c1iarg(cda, "a_col_dim", 1, -1, &ner);
    if (cda >= 1 && nca != 0) {
	imsl_c12ile(nca, "nca", cda, "a_col_dim", &ner);
    } else {
	ner += 1;
    }
    imsl_c1iarg(nca, "nca", 0, -1, &ner);
    imsl_null_pointer("title", -1, (void *)title);
    if (nra !=0 && nca !=0) imsl_null_pointer("a", -1, (void *)a);
    if (irlab == 2) {
        imsl_null_pointer("rlabel", -1, (void *)rlabel);
        m = (transpose?nra:nca);
        if (rlabel) for (i=0;i<m;i++) {
            imsl_null_pointer("rlabel", i, (void *)rlabel[i]);
            if (rlabel[i]==0)  goto end1;
        }
    }
end1:
    if (iclab == 2) {
        imsl_null_pointer("clabel", -1, (void *)clabel);
        m = (transpose?nca:nra);
        if (clabel) for (i=0;i<m;i++) {
            imsl_null_pointer("clabel", i, (void *)clabel[i]);
            if (clabel[i]==0)  goto end2;
        }
    }
end2:
    if (imsl_n1rty(0) > 0)
	goto RETURN;

    if (!fmt) {
	fmt = local;            /*determine optimal integer format */
	m = 1;
        if (transpose) {
            for (j=1; j<=nra; j++) {
    		    imsl_w12rl(itring, 1, nca, j, &i1, &nrow);
		    for (i = i1; i<=(i1 + nrow - 1); i++) {
			mij = abs(*(a+(j-1)*cda+(i-1)));
			if (mij > m)   m = mij;
                    }
            }
	} else {
            for (j=1; j <=nca; j++) {
    		    imsl_w12rl(itring, 1, nra, j, &i1, &nrow);
		    for (i = i1; i<=(i1 + nrow - 1); i++) {
			mij = abs(*(a+(j-1)+(i-1)*cda));
			if (mij > m)   m = mij;
                    }
            }
	}
	ndig = log10(m + .01) + 2;
        if (ndig > 9) {sprintf(fmt+1,"%2dd",ndig);
        } else sprintf(fmt+1,"%1dd", ndig);
        fmt[0]='%';
    }

    imsl_w1opt(3, &ipgopt);
    if ((ipgopt > 0) || (ipgopt == -2))   imsl_write_line(1, " ");

    if (irlab == -1) {
    	if (nra==1) irlab = 0;
    	else irlab = 1;
    }
    if (iclab == -1) {
    	if (nca==1) iclab = 0;
    	else iclab = 1;
    }

    if (nra == 0 || nca == 0)
        imsl_i_wrirl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                irlab, iclab);
    else if (transpose) {
        imsl_i_m1ran(nra, cda, a, a);
        imsl_i_wrirl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                irlab, iclab);
        imsl_i_m1ran (cda, nra, a, a);
    } else {
        imsl_i_wrirl(title, nca, nra, a, cda, itring, fmt, rlabel, clabel,
                irlab, iclab);
    }
    if (imsl_return_string)
    {
       if (user_output_string != 0)
          *user_output_string = imsl_output_string;
    }

RETURN:
    return (argptr);
}

#ifdef ANSI
void imsl_null_pointer(Mchar *description, Mint i, void *ptr)
#else
void imsl_null_pointer (description, i, ptr)
Mchar     *description;
Mint      i;
#ifdef COMPUTER_PMXUX
Mchar      *ptr;
#else
void       *ptr;
#endif
#endif
{
    if (!ptr) {
        imsl_e1stl(1, description);
        imsl_e1sti(1, i);
        if (i<0)  imsl_ermes(IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
        else imsl_ermes(IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINT);
    }
    return;
}
/* Structured by FOR_STRUCT, v0.2, on 06/12/90 at 12:58:46
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  WRIRL (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 5, 1985

    Purpose:    Print an integer rectangular matrix with a given format
                and labels.

    Usage:      CALL WRIRL (TITLE, NRMAT, NCMAT, MAT, LDMAT, ITRING, FMT,
                            RLABEL, CLABEL)

    Arguments:
       TITLE  - Character string specifying the title.  (Input)
                TITLE = ' ' suppresses printing of the title.
       NRMAT  - Number of rows.  (Input)
       NCMAT  - Number of columns.  (Input)
       MAT    - NRMAT by NCMAT matrix to be printed.  (Input)
       LDMAT  - Leading dimension of MAT exactly as specified in the
                dimension statement in the calling program.  (Input)
       ITRING - Triangle option.  (Input)
                ITRING   Action
                   0     Full matrix is printed.
                   1     Upper triangle of MAT is printed.
                   2     Upper triangle of MAT excluding the diagonal of
                         MAT is printed.
                  -1     Lower triangle of MAT is printed.
                  -2     Lower triangle of MAT excluding the diagonal of
                         MAT is printed.
       FMT    - Character string containing formats.  (Input)
                For example, FMT = '(I10)' specifies this I format for
                the entire matrix.  FMT = '(2I10, 3I5)' specifies
                an I10 format for columns 1 and 2 and an I5 format for
                columns 3, 4 and 5.  If the end of FMT is encountered and
                if some columns of the matrix remain, format control
                continues with the first format in FMT.  FMT may only
                contain the I edit descriptor, e.g., the X edit
                descriptor is not allowed.  FMT must contain exactly one
                set of parentheses.
       RLABEL - CHARACTER*(*) vector of labels for rows of MAT.  (Input)
                If rows are to be numbered consecutively 1, 2, ...,
                NRMAT, use RLABEL(1) = 'NUMBER'.  If no row labels are
                desired, use RLABEL(1) = 'NONE'.  Otherwise RLABEL is a
                vector of length NRMAT containing the labels.
       CLABEL - CHARACTER*(*) vector of labels for columns of MAT.
                (Input)
                If columns are to be numbered consecutively 1, 2, ...,
                NCMAT, use CLABEL(1) = 'NUMBER'.  If no column labels are
                desired, use CLABEL(1) = 'NONE'.  Otherwise, CLABEL is a
                vector of length NCMAT + 1 containing the column
                headings.  CLABEL(1) is the heading for the row labels,
                and for J = 1, 2, ..., NCMAT, CLABEL(1+J) is the column
                heading for the Jth column.

    Remarks:
    1. The output appears in the following form:
                           TITLE
       CLABEL(1)  CLABEL(2)  CLABEL(3)  CLABEL(4)
       RLABEL(1)      xxxxx      xxxxx      xxxxx
       RLABEL(2)      xxxxx      xxxxx      xxxxx

    2. Use '%/' within titles or labels to cause a new line to occur.
       Long titles or labels are automatically wrapped.

    3. A page width of 78 characters is used.  Page width and page length
       can be reset by invoking IMSL routine imsl_page.

    4. Printing options, such as centering and paging, can be selected
       by invoking IMSL routine imsl_write_options.

    5. Output is written to the unit specified by subroutine UMACH.

    Keywords:   Utilities; Output

    GAMS:       N1

    Chapters:   STAT/LIBRARY Utilities
                MATH/LIBRARY Utilities

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

void imsl_i_wrirl(title, nrmat, ncmat, mat, ldmat, itring, fmt, rlabel,
                    clabel, irlab, iclab)
	Mchar           *title;
	Mint             nrmat, ncmat, *mat, ldmat, itring, irlab, iclab;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
{
	Mchar            *tit;
	Mint             n;
	Mint             icentr, ihotzr, ipagel, ipagew,
                        ipgopt, ititle, jcol1, jcol2, jrow1, jrow2,
                        large, ldata, lrowlb, ltitle, maxfw, maxrlw, maxwl,
                        nrowlb, nstag, ntitle, wfmt, ido;


	imsl_write_format(nrmat, ncmat, ldmat, itring, fmt, &maxfw, "WdiouxX", "h", &n,
                &wfmt);

	if (imsl_n1rty(0)>3)  goto L_9000;

        if (imsl_write_initialize(&ipagew, &ipagel, &icentr, &large, &ipgopt,
                &ititle, rlabel, clabel, &iclab, &irlab, ncmat, &maxrlw, title,
                &ltitle, &ntitle, 1, maxfw, &lrowlb, &nrowlb, &ihotzr, nrmat,
                &maxwl)) goto L_9000;

        ido = 0;
L_10:   imsl_write_controller(&ido, large, iclab, irlab, clabel, rlabel, maxrlw,
                ihotzr, title, &tit, &nstag, lrowlb, &ldata, &jrow1, &jrow2,
                &jcol1, &jcol2, nrmat, ncmat, ipagew, ipgopt,
                ntitle, &ipagel, nrowlb, itring, ititle, fmt, &ltitle, 1);

        if (jrow1 > 0) {
	    l_w3irl_f(tit, mat, ldmat, itring, fmt, rlabel, clabel, jrow1,
                    jrow2, jcol1, jcol2, lrowlb, ldata, nstag, ltitle, irlab,
                    iclab, maxrlw, ihotzr, ipagew, icentr);
            if (ido > 0)  goto L_10;
        }

L_9000:	return;
}

/*-----------------------------------------------------------------------
    Revised:    October 16, 1985

    Purpose:    Get length of longest row label and the number of lines
                in the label requiring the most lines.

    Usage:      CALL W10RL (IRLAB, ICLAB, RLABEL, CLABEL, MAXRLW,
                            IHOTZR, NDENTR, IROW1, IROW2, LROWLB, NROWLB)

    Arguments:
       IRLAB  - Row label option.  (Input)
                IRLAB   Meaning
                  0     No row labels.
                  1     Row labels are 1, 2, 3, ...
                  2     Row labels are in RLABEL.
       ICLAB  - Column label option.  (Input)
                ICLAB   Meaning
                  0     No column labels.
                  1     Column labels are 1, 2, 3, ...
                  2     Column labels are in CLABEL.
       RLABEL - See WRRRL.  (Input)
       CLABEL - See WRRRL.  (Input)
       MAXRLW - Maximum permitted width for a row label.  (Input)
       IHOTZR - Length of hot zone for row label.  (Input)
       NDENTR - Indentation for continuation lines of a row label.
                (Input)
       IROW1  - First row label to be considered.  (Input)
       IROW2  - Last row label to be considered.  (Input)
       LROWLB - Length of the longest row label.  (Output)
       NROWLB - Number of lines in the row label requiring the most
                lines.  (Output)

  -----------------------------------------------------------------------
 */
static void  l_w10rl(irlab, iclab, rlabel, clabel, maxrlw, ihotzr, ndentr,
		    irow1, irow2, lrowlb, nrowlb)
	Mint             irlab, iclab;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mint             maxrlw, ihotzr, ndentr, irow1, irow2, *lrowlb,
			*nrowlb;
{
	Mint             done;
	Mint             ibeg, iend, irow, length, nline;


	if (iclab <= 1 || iclab == 3) {*lrowlb = 0;
	} else {
		nline = 0;
		l_w4rrl(clabel[0], maxrlw, 0, ihotzr, &nline, lrowlb, &ibeg,
			&iend, &done);
	}
	if (irlab == 0) {*nrowlb = 0;
	} else if (irlab == 1) {
		length = log10(irow2 + .01) + F_ONE;
		*lrowlb = max(*lrowlb, length);
		*nrowlb = 1;
	} else if (irlab == 3) {
		length = log10(irow2 - .99) + F_ONE;
		*lrowlb = max(*lrowlb, length);
		*nrowlb = 1;
	} else {
		*nrowlb = 0;
		for (irow = irow1; irow <= irow2; irow++) {
			nline = 0;
			l_w4rrl(rlabel[irow-1], maxrlw, ndentr, ihotzr,
				&nline, &length, &ibeg, &iend, &done);
			*lrowlb = max(*lrowlb, length);
			*nrowlb = max(*nrowlb, nline);
		}
	}
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Check if lv_part of matrix to be printed on this page is
                empty.

    Usage:      CALL W11RL (ITRING, IROW1, IROW2, ICOL1, ICOL2,
                            EMPTY, JROW1, JROW2, JCOL1, JCOL2)

    Arguments:
       ITRING - See WRRRL.  (Input)
       IROW1  - First row of matrix to be printed.  (Input)
       IROW2  - Last row of matrix to be printed.  (Input)
       ICOL1  - First column of matrix to be printed.  (Input)
       ICOL2  - Last column of matrix to be printed.  (Input)
       EMPTY  - LOGICAL variable.  (Output)
                EMPTY = .TRUE. means this lv_part of the matrix to be
                printed is empty, so printing for this page can be
                skipped.
       JROW1  - First nonempty row of matrix to be printed.  (Output if
                EMPTY = .FALSE.)
       JROW2  - Last nonempty row of matrix to be printed.  (Output if
                EMPTY = .FALSE.)
       JCOL1  - First nonempty column of matrix to be printed.  (Output
                if EMPTY = .FALSE.)
       JCOL2  - Last nonempty column of matrix to be printed.  (Output if
                EMPTY = .FALSE.)
  -----------------------------------------------------------------------
 */
static void l_w11rl(itring, irow1, irow2, icol1, icol2, empty,
	   jrow1, jrow2, jcol1, jcol2)
	Mint             itring, irow1, irow2, icol1, icol2;
	Mint             *empty;
	Mint             *jrow1, *jrow2, *jcol1, *jcol2;
{

	*empty = IFALSE;
	if (!itring) {
		*jrow1 = irow1;
		*jrow2 = irow2;
		*jcol1 = icol1;
		*jcol2 = icol2;
	} else if (itring == 1) {
		if (irow1 > icol2) {
			*empty = ITRUE;
		} else {
			*jrow1 = irow1;
			*jrow2 = min(irow2, icol2);
			*jcol1 = max(irow1, icol1);
			*jcol2 = icol2;
		}
	} else if (itring == 2) {
		if (irow1 >= icol2) {
			*empty = ITRUE;
		} else {
			*jrow1 = irow1;
			*jrow2 = min(irow2, icol2 - 1);
			*jcol1 = max(irow1 + 1, icol1);
			*jcol2 = icol2;
		}
	} else if (itring == -1) {
		if (irow2 < icol1) {
			*empty = ITRUE;
		} else {
			*jrow1 = max(irow1, icol1);
			*jrow2 = irow2;
			*jcol1 = icol1;
			*jcol2 = min(irow2, icol2);
		}
	} else if (itring == -2) {
		if (irow2 <= icol1) {
			*empty = ITRUE;
		} else {
			*jrow1 = max(irow1, icol1 + 1);
			*jrow2 = irow2;
			*jcol1 = icol1;
			*jcol2 = min(irow2 - 1, icol2);
		}
	}
	return;
}


/*-----------------------------------------------------------------------
    Purpose:    Write an integer rectangular matrix with a specified
                format and labels.

    Usage:      CALL W3IRL (TITLE, MAT, LDMAT, ITRING, FMT,
                            RLABEL, CLABEL, IROW1, IROW2, ICOL1,
                            ICOL2, LROWLB, LDATA, NSTAG, LTITLE, IRLAB,
                            ICLAB, MAXRLW, IHOTZR, IPAGEW,
                            ICENTR)

    Arguments:
       TITLE  - See WRIRL.
       MAT    - See WRIRL.
       LDMAT  - See WRIRL.
       ITRING - See WRIRL.
       FMT    - See WRIRL.
       RLABEL - See WRIRL.
       CLABEL - See WRIRL.
       IROW1  - Starting row to be printed.  (Input)
       IROW2  - Last row to be printed.  (Input)
       ICOL1  - Starting column to be printed.  (Input)
       ICOL2  - Last row to be printed.  (Input)
       LROWLB - Length of longest row label to be printed.  (Input)
       LDATA  - Width of longest printed row of data.  (Input)
       NSTAG  - Number of spaces indented for the most indented
                (staggered) row.  (Input)
       LTITLE - Length of the longest line in the title.  (Input)
       IRLAB  - Row label option.  (Input)
                IRLAB   Meaning
                  0     No row labels.
                  1     Row labels are 1, 2, 3, ... .
                  2     Row labels are in RLABEL.
       ICLAB  - Column label option.  (Input)
                ICLAB   Meaning
                  0     No column labels.
                  1     Column labels are 1, 2, 3, ... .
                  2     Column labels are in CLABEL.
       MAXRLW - Maximum permitted row label width.  (Input)
       IHOTZR - Length of hot zone for row labels.  (Input)
                The hot zone is the area where line breaks can occur.
       IPAGEW - Page width.  (Input)
       ICENTR - Centering option.  (Input)
                ICENTR   Action
                  0      No centering
                  1      Center
  -----------------------------------------------------------------------
 */

static void l_w3irl_f(title, mat, ldmat, itring, fmt, rlabel, clabel,
	irow1, irow2, icol1, icol2, lrowlb, ldata, nstag, ltitle, irlab,
	iclab, maxrlw, ihotzr, ipagew, icentr)

	Mchar           *title;
	Mint            *mat, ldmat, itring;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mint             irow1, irow2, icol1, icol2, lrowlb, ldata, nstag,
                        ltitle, irlab, iclab, maxrlw, ihotzr, ipagew,
                        icentr;
{
	Mchar            *fmti, line[256], blanks[256];
	Mint             i, done, icol, ipad, irow, istag, lline, lout, middle,
                         nline, nlinei, nwidth;

				/* Get width of output, LOUT. */
	for (i=0;i<256;i++) blanks[i]=' ';blanks[255] = '\0';

	lout = max(ltitle, ldata);
	if (icentr == 1) middle = ipagew / 2;
	else middle = lout / 2;

        strcpy(line,blanks);

	nline = 0;
	imsl_write_title (title, ipagew, &nline, line, middle);
	strcpy(line,blanks);

        imsl_write_labels (iclab, clabel, line, blanks, &nline, lrowlb, maxrlw, ihotzr,
                    middle, icol1, icol2, ipagew, nstag, fmt, ldata, 0);
	strcpy(line,blanks);

                                                	/* Print MAT. */
	for (irow = irow1; irow <= irow2; irow++) {
		nline = 1;
		imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr, NDENTR,
			nline, &lline, line, &done);
		if (lline > 0)  lline += NSPACE;
		istag = 0;
		for (icol = icol1; icol <= icol2; icol++) {
			                        /*Get format width*/
			imsl_w5rrl_f(icol, iclab, clabel, fmt, &nwidth, &ipad,
				&fmti, &nlinei, 1);
		                	/* Get padding before the number. */
			ipad -= nwidth;
				       /* Check to see if current line fits */
			if (lline + ipad + nwidth > ipagew - nstag + istag) {
				                        /* Line won't fit*/
				lline = ldata;
				imsl_c1nter(middle, &lline, line);
				imsl_write_line(lline, line);
				strcpy(line,blanks);
				nline += 1;
				imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw,
                                        ihotzr, NDENTR, nline, &lline, line,
					&done);
				if (lline > 0)  lline += NSPACE;
				istag += 1;
				if (istag == 5)  istag = 0;
				lline += istag;
			}
   			          /* Append entry MAT(IROW,ICOL) to LINE. */
			lline += ipad;
                        strncpy(line+lline, l_i1tri(itring,irow,icol,
                                 *MAT(icol-1,irow-1), fmti, nwidth), nwidth);

			lline += nwidth + NSPACE;
		}
		                    /*Print last line for IROW-th row of A.*/
		lline = ldata;
		imsl_c1nter(middle, &lline, line);
		imsl_write_line(lline, line);
		strcpy(line,blanks);
L_100:
		if (!done) {
			/* Number of lines for IROW-th row label exceeds that
			 * for the IROW-th row of data. Print the remainder
			 * of the row label.*/
			nline += 1;
			imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr,
				NDENTR, nline, &lline, line, &done);
			lline = ldata;
			imsl_c1nter(middle, &lline, line);
			imsl_write_line(lline, line);
			strcpy(line,blanks);
			goto L_100;
		}
	}
}
/*
  -----------------------------------------------------------------------
    Purpose:    Convert an integer matrix element to a character string
                according to a supplied format.  If the matrix element is
                not in the triangle to be printed, blanks are returned.

    Usage:      WRITE (NOUT,1) I1TRI(ITRI, I, J, NUM, FMT)

    Arguments:
       ITRI   - Triangle option.  (Input)
                  ITRI    Action
                    0     Full matrix is printed.
                    1     Upper triangle of matrix is printed.
                    2     Upper triangle of matrix excluding the
                          diagonal is printed.
                   -1     Lower triangle of matrix is printed.
                   -2     Lower triangle of matrix excluding the
                          diagonal is printed.
       I      - Row number of matrix element.  (Input)
       J      - Column number of matrix element.  (Input)
       NUM    - The matrix element to be written.  (Input)
       FMT    - The format to be used in writing NUM.  (Input)
       I1TRI  - A CHARACTER*40 string which contains blanks if the matrix
                element is not to be printed.  (Output)
  ----------------------------------------------------------------------
 */
static Mchar *l_i1tri(itri, i, j, num, fmt, nwidth)
	Mint             itri, i, j, num, nwidth;
	Mchar           *fmt;
{
	static Mchar     i1tri_v[41];
        Mchar           *stars = "*****************************************";
	Mint             iwrite, il;
        Mchar            tmpfmt[12], *iend;

	iwrite = 0;
	if (itri == 0) {iwrite = 1;
	} else if (itri == 1 && i <= j) {iwrite = 1;
	} else if (itri == 2 && i < j) {iwrite = 1;
	} else if (itri == -1 && i >= j) {iwrite = 1;
	} else if (itri == -2 && i > j) iwrite = 1;

	if (!iwrite) {strcpy(i1tri_v, "                                        ");
	} else {
            iend = strchr(fmt+1,(Mint)'%');
            if (iend) {
                il = (Mint)(iend-fmt);
                il = (il>11) ? 11:il;
                strncpy (tmpfmt, fmt, il);
                tmpfmt[il] = '\0';
            } else {
                strncpy (tmpfmt, fmt, 11);
                tmpfmt[11]='\0';
            }
            sprintf(i1tri_v,tmpfmt,num);
            if (strlen(i1tri_v) > nwidth)  strncpy(i1tri_v,stars,nwidth);
        }
	return (i1tri_v);
}

void imsl_write_title (title, ipagew, nline, line, middle)
	Mchar    *title, *line;
	Mint     ipagew, *nline, middle;
{
	Mint ibeg, iend, done=0, lline;

	if (strlen(title)) {
					  /* Print title. */
		while (!done) {
			*nline += 1;
			l_w4rrl(title, ipagew, 0, 10, nline, &lline, &ibeg,
				&iend, &done);
			if (lline > 0) {
					/* Write NLINE-th line of title. */
				strncpy(line, &title[ibeg-1],iend-ibeg+1);
				imsl_c1nter(middle, &lline, line);
				imsl_write_line(lline, line);
			} else if (lline == 0) {
			    		   /* Write a blank line. */
				imsl_write_line(1, " ");
			}
		}
		                           /* Title printing completed. */
	}
	return;
}

void imsl_write_labels (iclab, clabel, line, blanks, nline, lrowlb, maxrlw, ihotzr,
                    middle, icol1, icol2, ipagew, nstag, fmt, ldata, comp)
        Mint     iclab, *nline, lrowlb, maxrlw, ihotzr, middle, icol1,
                icol2, ipagew, nstag, ldata, comp;
        Mchar    **clabel, *blanks, *fmt, *line;

{
        Mint     nlinei, nline0, done0, dones, lline, istag, ipad,
                maxclw, icol11, lfield, nwidth[2], icol, ihotzc, lclab, done,
                ibeg, iend;
        Mchar    temp[71], *fmti[2];

	if (iclab == 0) return;      /* Skip printing of column headings. */
	if (iclab == 2) {
		if (strlen(*clabel))  goto L_30;
		for (icol = icol1; icol <= icol2; icol++) {
			if (strlen(clabel[icol]))  goto L_30;
		}
		/* All column headings are blank. Skip printing of column
		 * headings.  */
		return;
	}
L_30:   nline0 = 1;  	                         /* Print column headings. */
	strcpy(line,blanks);
	if (iclab == 1 || iclab == 3) {
		lline = lrowlb;
		done0 = ITRUE;
	} else {
		imsl_w8rrl(1, iclab, clabel, lrowlb, maxrlw, ihotzr, 0, nline0,
			&lline, line, &done0);
	}
	if (lline > 0)  lline += NSPACE;
	icol = icol1 - 1;
	istag = -1;
L_40:   *nline = 0;        /* Begin new subset of column labels. */
	istag += 1;
	if (istag == 5)  istag = 0;
	icol11 = icol;
	         /* Begin next line for the current subset of column labels.*/
L_50:	lline += istag;
	dones = ITRUE;
	*nline += 1;
	icol = icol11;
	if (icol2 == 0)  goto L_70;
L_60:   icol += 1;              /* Get format and widths for ICOL-th column.*/
        if (!comp) {
	    imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth, &lfield, fmti,
                    &nlinei, 1);
        } else {
            imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth, &lfield, fmti,
                    &nlinei, 2);
        }

				/* Check to see if current line fits */
	if (lline + lfield > ipagew - nstag + istag) {
		/* Line won't fit so print old line, and begin new line.*/
		lline = ldata;
		imsl_c1nter(middle, &lline, line);
		imsl_write_line(lline, line);
		nline0 += 1;
		if (iclab == 1 || iclab == 3 || done0) {
			lline = lrowlb;
			strcpy(line,blanks);
		} else {
			imsl_w8rrl(1, iclab, clabel, lrowlb, maxrlw, ihotzr, 0,
				nline0, &lline, line, &done0);
		}
		if (lline > 0)
			lline += NSPACE;
		if (dones) {   /* Begin new subset of column labels. */
			icol -= 1;
			goto L_40;
		} else goto L_50;  /* Start next line for current subset. */
	}
	maxclw = min(max(nwidth[0] + nwidth[0] / 2, 15), 40);
	ihotzc = maxclw / 3;
	if (iclab == 1) {
		if (*nline == 1) {
                        lclab = log10(icol + .01) + 1;
			strncpy(temp,blanks,lfield-lclab);
                        sprintf(&temp[max(lfield,lclab)-lclab],"%d",icol);
                        lclab = max(lfield,lclab);
		} else {
			lclab = 0;
		}
		done = ITRUE;
	} else if (iclab == 3) {
		if (*nline == 1) {
                        lclab = log10(icol - .01) + 1;
                        if (lclab <= 0) lclab = 1;
			strncpy(temp,blanks,lfield-lclab);
                        sprintf(&temp[max(lfield,lclab)-lclab],"%d",icol-1);
                        lclab = max(lfield,lclab);
		} else {
			lclab = 0;
		}
		done = ITRUE;
	} else {
		l_w4rrl(clabel[icol], maxclw, 0, ihotzc, nline, &lclab,
			&ibeg, &iend, &done);
	}
	dones = dones && done;
				/* Get padding before the column label. */
	ipad = lfield - lclab;
	lline += ipad;
			/* Append the NLINE-th line for CLABEL(ICOL) to LINE. */
	if (iclab == 1 || iclab == 3) {
		if (lclab > 0)  strncpy(line+lline,temp,lclab);
	} else {
		if (lclab > 0)
			strncpy(&line[lline],&clabel[icol][ibeg-1],iend-ibeg+1);
	}
	lline += lclab + NSPACE;
L_70:
	if (icol == icol2) {
		lline = ldata;
		imsl_c1nter(middle, &lline, line);
		imsl_write_line(lline, line);
		if (!dones || !done0) {
 			                  /* Begin a new line. */
			nline0 += 1;
			strcpy(line,blanks);
			if (iclab == 1 || iclab == 3 || done0) {
				lline = lrowlb;
			} else {
				imsl_w8rrl(1, iclab, clabel, lrowlb, maxrlw,
                                        ihotzr, 0, nline0, &lline, line,
					&done0);
			}
			if (lline > 0)  lline += NSPACE;
			goto L_50;
		}
	} else {
		goto L_60;
	}
        return;
}

/*-----------------------------------------------------------------------
    Purpose:    Get length, beginning and ending position for the lv_part of
                a label or title to be printed on NLINE-th line.

    Usage:      CALL W4RRL (LABEL, MAXWL, INDENT, IHOTZ, NLINE, LENGTH,
                            IBEG, IEND, DONE)

    Arguments:
       LABEL  - Character string containing the label or title.  (Input)
       MAXWL  - Maximum width allowed for lv_part of label on line.  (Input)
       INDENT - Indentation used for continuation lines.  (Input)
       IHOTZ  - Length of hot zone where break in label may occur.
                (Input)
       NLINE  - For positive NLINE, the number of the line where lv_part of
                label appears; NLINE = 0 means output the number of lines
                required by the label.  (Input if NLINE is positive;
                input/output if NLINE = 0 on input)
       LENGTH - Length of NLINE-th lv_part of label.  (Output)
                For NLINE = 0, LENGTH is the length of the longest lv_part
                of the label.  (For continuation lines of a label, INDENT
                is included in LENGTH.)  LENGTH = 0 if the label is a
                null or blank line.
                For positive NLINE, LENGTH = 0 if the lv_part of the label
                to be printed is a null or blank line.  Otherwise,
                LENGTH is the length of the lv_part of the label to be
                printed on the NLINE-th line.  (For NLINE greater than
                1, indent is included in LENGTH.)
       IBEG   - Beginning position.  (Output if NLINE is positive and
                LENGTH is positive)
       IEND   - Ending position.  (Output if NLINE is positive and
                LENGTH is positive)
                LABEL(IBEG:IEND) is the lv_part of the label which is to
                appear on line number NLINE.
       DONE   - LOGICAL variable.  (Output)
                DONE = .TRUE. if information for the last lv_part of the
                label has been retrieved.
  -----------------------------------------------------------------------
 */
static void l_w4rrl(label, maxwl, indent, ihotz, nline, length, ibeg, iend,
                        done)
	Mchar           *label;
	Mint             maxwl, indent, ihotz, *nline, *length, *ibeg, *iend;
	Mint            *done;
{
	Mint             i, ibeg1, iblank, iend1, ihot, il, iline, j,
                        len1, llab, maxht, ndent;
	Mchar            *ctmp, *itmp;

        llab=strlen(label);
        ctmp = (Mchar *)malloc(llab);
	*done = 0;
	*length = 0;
	iline = 0;
	if (llab == 0) {
		*done = ITRUE;
		goto L_9000;
	}
	ibeg1 = 1;
	ndent = 0;
	maxht = maxwl - ihotz;
L_10:
	iline += 1;
	ihot = -ndent + ibeg1 + maxht;
	if (llab < ibeg1)  iline -= 1;
	for (i = ibeg1; i <= llab; i++) {
		if (i == llab) {             /* End of label. */
			*done = ITRUE;
			if (iline == *nline || *nline == 0) {
				len1 = llab - ibeg1 + 1 + ndent;
				if (*nline == 0) {
					*length = max(*length, len1);
				} else {
					*length = len1;
					*ibeg = ibeg1;
					*iend = llab;
				}
			}
		} else if (label[i-1]=='\n') {
					/* End of ILINE-th line. */
			if (i == llab)  *done = ITRUE;
			if (iline == *nline || *nline == 0) {
				iend1 = i-1;
				if (iend1 < ibeg1) {
						/* Null line encountered. */
					len1 = 0;
/*				} else if (label[ibeg1-1]==' ') {
						*//* Blank line encountered. *//*
					len1 = 0;*/
				} else {
					len1 = iend1 - ibeg1 + 1 + ndent;
				}
				if (*nline == 0) {
					*length = max(*length, len1);
				} else {
					*length = len1;
					if (*length > 0) {
						*ibeg = ibeg1;
						*iend = iend1;
					}
					goto L_9000;
				}
			}
			ibeg1 = i + 1;
			ndent = indent;
			goto L_10;
		} else if (i >= ihot) {
					/* In hot zone.  */
			itmp = strchr(&label[i], (Mint)'\n');
                        il = ((itmp)?(Mint)(itmp - &label[i]) + 1:0);
			if (il > ihotz + 1) il = 0;

			if (il >= 1) {
						/* New line in hot zone. */
				if (i + il == llab)  *done = ITRUE;
				if (iline == *nline || *nline == 0) {
					len1 = i + il - ibeg1 + ndent - 1;
					if (*nline == 0) {
						*length = max(*length, len1);
					} else {
						*ibeg = ibeg1;
						*iend = i + il - 1;
						*length = len1;
						goto L_9000;
					}
				}
				ibeg1 = i + il + 1;
				ndent = indent;
				goto L_10;
			} else if (llab < ihot + ihotz) {
						/* End of label in hot zone. */
				*done = ITRUE;
				if (iline == *nline || *nline == 0) {
					len1 = llab - ibeg1 + 1 + ndent;
					if (*nline == 0) {
						*length = max(*length, len1);
					} else {
						*length = len1;
						*ibeg = ibeg1;
						*iend = llab;
					}
				}
				goto L_9000;
			} else {
				itmp = strchr(&label[i-1], (Mint)' ');
                                iblank = ((itmp)?(Mint)(itmp-&label[i-1])+1:0);
				if (iblank > ihotz + 1)  iblank = 0;
				if (iblank >= 1) {
						/* Blank in hot zone. */
					if (iline == *nline || *nline == 0) {
						iend1 = i + iblank - 2;
						len1 = ndent + iend1 - ibeg1 + 1;
						if (*nline == 0) {
							*length = max(*length, len1);
						} else {
							*length = len1;
							*ibeg = ibeg1;
							*iend = iend1;
							goto L_9000;
						}
					}
					ndent = indent;
					ibeg1 = i + iblank;
					goto L_10;
				} else {
						/*Look in hot zone for any
						non-alphanumeric character.*/
					for (j = i;j<=(i + ihotz - 1); j++) {
						if (!isalnum((Mint)label[j-1]))
                                                        break;
					}
					if (iline == *nline || *nline == 0) {
						iend1 = j - 1;
						len1 = ndent + iend1 - ibeg1 + 1;
						if (*nline == 0) {
							*length = max(*length, len1);
						} else {
							*ibeg = ibeg1;
							*iend = iend1;
							*length = ndent + *iend - *ibeg + 1;
							goto L_9000;
						}
					}
					ibeg1 = j;
					ndent = indent;
					goto L_10;
				}
			}
		}
	}
L_9000:
	if (*nline == 0)  *nline = iline;
        free(ctmp);
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Determine width of field for column from the format
                and column label and the number of lines required by
                the column label.

    Usage:      CALL W5RRL (ICOL, ICLAB, CLABEL, FMT, NWIDTH, LFIELD,
                            FMTI, NLINES, MODE)

    Arguments:
       ICOL   - Column number.  (Input)
       ICLAB  - Column label option.  (Input)
                ICLAB   Meaning
                  0     No column labels.
                  1     Column labels are 1, 2, 3, ...
                  2     Column labels are in CLABEL.
       CLABEL - CHARACTER(*)*(*) containing the column labels.  (Input)
       FMT    - CHARACTER*20 string containing the format.  (Input)
       NWIDTH - Width of field specified in FMT.  (Output)
       LFIELD - Length of field required by FMT and LABEL.  (Output)
       FMTI   - CHARACTER*20 string with format for column
                ICOL.  (Output)
       NLINES - Number of lines required by the column label.  (Output)
       MODE   - COMPLEX = 2, OTHER = 1.  (Input)
  -----------------------------------------------------------------------
 */
void imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth, lfield, fmti, nlines, mode)
	Mint             icol, iclab, mode;
	Mchar           **clabel;
	Mchar           *fmt;
	Mint            *nwidth, *lfield;
	Mchar           **fmti;
	Mint            *nlines;
{
	Mint             done, nwid, div;
	Mint             ibeg, iend, ihotz, length, maxwl;
	Mchar            aa='.', bb='W';


	if (mode == 1) {
	    *fmti = imsl_w7rrl(icol, fmt);
	    imsl_w6rrl(*fmti, 1, &aa, &bb, nwidth);
            if (*nwidth==0)  *nwidth=10;
            div = 40;
            nwid = *nwidth;
        }else {
            fmti[0] = imsl_w7rrl(2*icol-1, fmt);
	    imsl_w6rrl(fmti[0], 1, &aa, &bb, nwidth);
            if (nwidth[0]==0)  nwidth[0]=10;
            fmti[1] = imsl_w7rrl(2*icol, fmt);
	    imsl_w6rrl(fmti[1], 1, &aa, &bb, nwidth+1);
            if (nwidth[1]==0)  nwidth[1]=10;
            div = 83;
            nwid = nwidth[0] + nwidth[1] + 3;
        }

	if (iclab == 0) {
		*lfield = nwid;
		*nlines = 0;
	} else if (iclab == 1) {
		*lfield = log10(icol + .01) + F_ONE;
		*lfield = max(*lfield, nwid);
		*nlines = 1;
	} else if (iclab == 3) {
		*lfield = log10(icol - .99) + F_ONE;
		*lfield = max(*lfield, nwid);
		*nlines = 1;
	} else {
		*lfield = nwid;
		maxwl = min(max(nwid + nwid / 2, 15), div);
		ihotz = maxwl / 3;
		*nlines = 0;
		l_w4rrl(clabel[icol], maxwl, 0, ihotz, nlines, &length,
			&ibeg, &iend, &done);
		*lfield = max(*lfield, length);
	}
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Get integer between two delimiters in a character string.

    Usage:      CALL W6RRL (STRING, IDEL, DEL1, DEL2, I)

    Arguments:
       STRING - Character string.  (Input)
                STRING must contain DEL1, DEL2 and no blank between
                DEL1 and DEL2.  Also only numbers may appear between
                the delimiters, and there must be at least one digit.
       IDEL   - Delimeter option.  (Input)
                IDEL1  Action
                  0    Delimiters are input in DEL1 and DEL2.
                  1    Width is retrieved from a FORTRAN format.
       DEL1   - CHARACTER*1 string containing first delimiter.  (Input
                if IDEL = 0)
       DEL2   - CHARACTER*1 string containing second delimiter.  (Input
                if IDEL = 0)
       I      - Integer between the two delimiters.  (Output)
  -----------------------------------------------------------------------
*/
void imsl_w6rrl(strng, idel, del1, del2, i)
	Mchar           *strng;
	Mint            idel;
	Mchar           *del1, *del2;
	Mint            *i;
{
	Mint             lenint;
        Mchar            *i1;

	if (!idel) {
		i1 = strchr(strng, (Mint)(*del1));
		if (i1) {
                        i1 += 1;
                        lenint = (Mint)(strchr(i1,(Mint)(*del2))-i1);
                } else lenint = 0;
	} else {
		i1 = strng+1;
                while (strchr(" #+-0",(Mint)(*i1))) i1++;
	        lenint = 0;
                while (isdigit((Mint)(*(i1+lenint)))) lenint++;
	}
	if (lenint) imsl_c1tci_f(i1, lenint, i);
	else *i = 0;
        return;
}
/*-----------------------------------------------------------------------
    Usage:      CALL W7RRL (ICOL, FMT)

    Arguments:
       ICOL   - Number of the column format to be obtained.  (Input)
       FMT    - User specified format list.  (Input)
  -----------------------------------------------------------------------
 */
Mchar *imsl_w7rrl(icol, fmt)
	Mint             icol;
	Mchar           *fmt;
{
	Mchar            *new;
	static Mint      onefmt, lscol, n;
        static Mchar     *now;
	Mint             i;

	if (icol == 1) {
                now = strchr(fmt,(Mint)'%');
                new = now;
                n = 1;
                while ((new = strchr(new+1,(Mint)'%')))  n++;
                onefmt = (n>1?0:1);
        }
	if (onefmt==1) return(now);

	if (icol<=lscol || icol==1) {
            lscol = n*((icol-1)/n);
	    new = fmt+strlen(fmt)-1;
	} else new = now;
        for (i=lscol+1; i<=icol; i++) {
	    if (!(new=strchr(new+1,(Mint)'%')))  new = strchr(fmt,(Mint)'%');
        }
	now = new+1;
        lscol = icol;
        return(new);
}

/* -----------------------------------------------------------------------
    Purpose:    Construct beginning of line with row label or leading
                blanks.

    Usage:      CALL W8RRL (IROW, IRLAB, RLABEL, LROWLB, MAXRLW,
                            IHOTZR, INDENT, NLINE, LLINE, LINE, DONE)

    Arguments:
       IROW   - Row number of matrix to be printed.  (Input)
       IRLAB  - Row label option.  (Input)
                IRLAB   Meaning
                  0     No row labels.
                  1     Row labels are 1, 2, 3, ...
                  2     Row labels are in RLABEL.
       RLABEL - CHARACTER*(*) vector of row labels.  (Input if
                IRLAB = 2)
       LROWLB - Length of row label field.  (Input)
       MAXRLW - Maximum row label length permitted.  (Input)
       IHOTZR - Hot zone for row labels.  (Input)
       INDENT - Indentation for continuation lines for row labels.
                (Input)
                INDENT must equal 1, 2, or 3.
       NLINE  - Number of the line for the lv_part of the row label.
                (Input)
                For example, NLINE = 2 means this is the second line,
                i.e., the first continuation line.  NLINE must be
                positive.
       LLINE  - Length of LINE output.  (Output)
       LINE   - Character string of length LLINE with beginning of
                the NLINE-th line for the IROW-th of the matrix to
                be printed.  (Output)
       DONE   - LOGICAL variable.  (Output)
                DONE = .TRUE. means the row label has been finished.
  -----------------------------------------------------------------------
 */
void imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr, indent, nline,
        lline, line, done)
	Mint             irow, irlab;
	Mchar           **rlabel;
	Mint             lrowlb, maxrlw, ihotzr, indent, nline, *lline;
	Mchar           *line;
	Mint            *done;
{
	Mchar            temp[40];
	Mint             i, ibeg, iend, length, ndig;

	if (irlab == 0) {
		*done = ITRUE;
	} else if (irlab == 1) {
		ndig = log10(irow + .01) + F_ONE;

                for (i=0;i<lrowlb-ndig;i++)  temp[i]=' ';
                sprintf(&temp[max(lrowlb,ndig)-ndig],"%d",irow);
		ndig = max(ndig,lrowlb);
                temp[ndig]='\0';
		l_w4rrl(temp, maxrlw, indent, ihotzr, &nline, &length, &ibeg,
			&iend, done);
		if (length >= 1) {
			if (nline == 1) {
				strncpy(line, &temp[ibeg-1], iend-ibeg+1);
			} else {
				strncpy(&line[indent],&temp[ibeg-1],
                                        iend-ibeg+1);
			}
		}
	} else if (irlab == 3) {
		ndig = log10(irow - .99) + F_ONE;
                if (ndig <= 0)  ndig = 1;

                for (i=0;i<lrowlb-ndig;i++)  temp[i]=' ';
                sprintf(&temp[max(lrowlb,ndig)-ndig],"%d",irow-1);
		ndig = max(ndig,lrowlb);
                temp[ndig]='\0';
		l_w4rrl(temp, maxrlw, indent, ihotzr, &nline, &length, &ibeg,
			&iend, done);
		if (length >= 1) {
			if (nline == 1) {
				strncpy(line, &temp[ibeg-1], iend-ibeg+1);
			} else {
				strncpy(&line[indent],&temp[ibeg-1],
                                        iend-ibeg+1);
			}
		}
	} else {
		l_w4rrl(rlabel[irow-1], maxrlw, indent, ihotzr, &nline,
			&length, &ibeg, &iend, done);
		if (length >= 1) {
			if (nline == 1) {
				strncpy(line,&rlabel[irow-1][ibeg-1],
                                        iend-ibeg+1);
			} else {
                                strncpy(&line[indent],&rlabel[irow-1][ibeg-1],
                                        iend-ibeg+1);
			}
		}
	}
	*lline = lrowlb;
	return;
}

/*-----------------------------------------------------------------------
    Usage:      CALL W9IRL (NRMAT, NCMAT, LDMAT, ITRING, FMT, MAXFW)

    Arguments:
       NRMAT  - See WRIRL.  (Input)
       NCMAT  - See WRIRL.  (Input)
       LDMAT  - See WRIRL.  (Input)
       ITRING - See WRIRL.  (Input)
       FMT    - See WRIRL.  (Input)
       MAXFW  - Maximum field width of FMT.  (Output)
  -----------------------------------------------------------------------
 */
/* NRMAT, NCMAT, and IDMAT are not used, but leave the calling sequence
   intact. */
void imsl_write_format(nrmat, ncmat, ldmat, itring, fmt, maxfw, formats,
             prefixes, n, wfmt)
	Mint             nrmat, ncmat, ldmat, itring;
	Mchar           *fmt, *formats, *prefixes;
	Mint            *maxfw, *n, *wfmt;
{
	Mint            lfmt, iwarn;
	Mchar           *isrch, *i;

	*maxfw = -1;
	lfmt = strlen(fmt);
	if (lfmt <= 2) {
			/* The specified format is an invalid C format */
			/* string: fmt = %(L1). */
		imsl_e1stl(1, fmt);
		imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_INVALID);
		return;
	}
	i = &fmt[0];
	iwarn = 0;
	*n = 0;
        *wfmt = 0;
	while (i<&fmt[lfmt]) {
 		isrch = strchr(i, (Mint)'%');
                if (isrch) {
    		    if (isrch != i){
	        	if (!iwarn){
			/* A character which is not a valid conversion */
                        /* specification character for this data type */
                        /* has been encountered in a format string. */
                        /* Such characters are ignored.  fmt = %(L1). */
                            imsl_e1stl(1,fmt);
			    imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_INVALID_CONV);
			}
		    }
		    i = imsl_write_conversion(isrch,maxfw,formats,prefixes, wfmt);
		    *n += 1;
		    if (!i) return;
		} else {
                        if (*n == 0) {
			    /* A format with no conversion specification */
                            /* has been encountered.  The format is */
                            /* FMT = %(L1). */
                            imsl_e1stl(1,fmt);
                            imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_NO_CONV);
                        } else {
                            imsl_e1stl(1,fmt);
			    imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_INVALID_CONV);
                        }
                        i = &fmt[lfmt];
                }
	}
}

Mchar *imsl_write_conversion(conv,maxfw,formats,prefixes, wfmt)
	Mchar	*conv, *formats, *prefixes;
	Mint     *maxfw, *wfmt;
{
	Mchar    *i=conv+1, *j, *null=0;
	Mint     nwidth, nsig;

	while (strchr("+-#0 ",(Mint)*i) && *i) i++;

        if (*i == '*') {
		 /* A C conversion specification with a '*'	    */
                 /* for a field width is not allowed.  fmt=%(L1).   */
		 imsl_e1stl(1,conv);
		 imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_STAR);
		 return(null);
	}
	j = i;
	while (isdigit((Mint)(*i))) i++;
	if (j == i){
		*maxfw = max(*maxfw,10);
                nwidth = 0;
	}else{
		imsl_c1tci_f(j,(Mint)(i-j),&nwidth);
		*maxfw = max(*maxfw,nwidth);
                if (nwidth>40)  {
		    /* A field width greater than 40 has been	    */
                    /* encountered.  All field widths must be 40    */
                    /* or less. fmt = %(L1)			    */
                    imsl_e1stl(1,conv);
                    imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_WIDE_FIELD);
                    return(null);
                }
	}
	if (*i == '.') {
		i++;
                j=i;
		while (isdigit((Mint)(*i))) i++;
		imsl_c1tci_f(j,(Mint)(i-j),&nsig);
                if (nsig>nwidth)  {
		    /* The number of significant digits to be printed	*/
                    /* is greater than the field width.  This		*/
                    /* conversions specification is not allowed		*/
                    /* here. fmt = %(L1)				*/
                    imsl_e1stl(1,conv);
                    imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_TOO_NARROW);
                    return(null);
                }
	}
                /* Skip over the single allowed prefix */

	if (strlen(prefixes) && *i) if (strchr(prefixes,(Mint)*i)) i++;

	if (strchr(formats, (Mint)*i) && *i) {
                if (*i == 'W') *wfmt += 1;
		i++;
		return(i);
	}else{
			/* An invalid conversion specification has been */
			/* encountered, fmt = %(L1). */
		imsl_e1stl(1,conv);
		imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_BAD_CONV);
		return(null);
	}
}

Mint imsl_c1tci_f(chrstr, slen, num)
	Mchar           *chrstr;
	Mint             slen, *num;
{
	Mint             imach5=INT_MAX;
	Mint             sgn=1;
        Mchar           *i;

	*num = 0;
	i = chrstr;
	while (*i==' ' && i < chrstr+slen) i++;
	if (i==chrstr+slen) return(0);
	while (i < chrstr+slen) {
		switch (*i)
		{
			case '-':
				sgn = -1;
				break;
			case '+':
				break;
			default:
				if (!isdigit(*i)) {
					return((Mint)*i);
				} else if (*num > (imach5-(Mint)(*i)) / 10) {
				       return(-2);
				}else *num = *num*10+(Mint)(*i)-(Mint)('0');
		}
		i++;
        }
        *num = *num*sgn;
	return(0);
}

/*-----------------------------------------------------------------------
    Purpose:    Write a line to output device.

    Usage:      CALL W1LIN (LLINE, LINE)

    Arguments:
       LLINE  - Length of line.  (Input)
                LLINE = -1 means begin a new page.
       LINE   - Character string containing the line to be printed.
                (Input if LLINE is positive)
  -----------------------------------------------------------------------
 */
void imsl_write_line(lline, line)
	Mint             lline;
	Mchar           *line;
{
	Mint             ipagel, ipgopt;
	FILE            *nout = NULL;

	if (lline == 0) return;
	imsl_w1opt(3, &ipgopt);
	if (imsl_return_string) {
		Mchar	*tmp_ptr, *tmp_string;
		int	tmp_string_malloced = 0;
		int	add_to_line = 0;

		if (lline == -1) {
			tmp_string = "\f";
			lline = 0;
			add_to_line++;
			if (ipgopt >= 0)  {ipgopt=0;imsl_w1opt(-3, &ipgopt);}
		} else if (lline >= 1) {
			tmp_ptr = tmp_string = (Mchar*)malloc((lline+3)*sizeof(Mchar));
			tmp_string_malloced = 1;
			if (ipgopt >= 0) {
				imsl_page(2, &ipagel);
				if (ipgopt >= ipagel) {
					strcpy(tmp_string, "\f");
					add_to_line++;
					tmp_ptr++;
					ipgopt = 0;
				}
				ipgopt+=1;
				imsl_w1opt(-3, &ipgopt);
			}
			sprintf(tmp_ptr, "%*.*s\n", lline, lline, line);
			add_to_line++;
		}
		imsl_output_string_length += lline + add_to_line;
		if (imsl_output_string == NULL)
			imsl_output_string = (Mchar*)malloc(imsl_output_string_length*sizeof(*imsl_output_string));
		else
			imsl_output_string = (Mchar*)realloc(imsl_output_string, imsl_output_string_length*sizeof(*imsl_output_string));
		strcpy(imsl_output_string+imsl_output_string_length-lline-1, tmp_string);
		if (tmp_string_malloced) free(tmp_string);
#if 0
	} else if (imsl_write_to_console) {
		HANDLE	console_handle = GetStdHandle(STD_OUTPUT_HANDLE);

		if (lline == -1) {
			WriteConsole(console_handle, "\f", 1, NULL, NULL);
			if (ipgopt >= 0)  {ipgopt=0;imsl_w1opt(-3, &ipgopt);}
		} else if (lline >= 1) {
			Mchar	*tmp_ptr, *tmp_string;

			tmp_ptr = tmp_string = (Mchar*)malloc((lline+3)*sizeof(Mchar));
			if (ipgopt >= 0) {
				imsl_page(2, &ipagel);
				if (ipgopt >= ipagel) {
					sprintf(tmp_string, "\f");
					lline++;
					tmp_ptr++;
					ipgopt = 0;
				}
				ipgopt+=1;
				imsl_w1opt(-3, &ipgopt);
			}
			sprintf(tmp_ptr, "%*.*s\n", lline, lline, line);
			WriteConsole(console_handle, tmp_string, lline+1, NULL, NULL);
			free(tmp_string);
		}
#endif
	} else {
		if (lline == -1) {
			imsl_umach(2,&nout);
			fprintf(nout,"\f");
			if (ipgopt >= 0)  {ipgopt=0;imsl_w1opt(-3, &ipgopt);}
		} else if (lline >= 1) {
			if (ipgopt >= 0) {
				imsl_page(2, &ipagel);
				if (ipgopt >= ipagel) {
					imsl_umach(2,&nout);
					fprintf(nout, "\f");
					ipgopt = 0;
				}
				ipgopt+=1;
				imsl_w1opt(-3, &ipgopt);
			}
			imsl_umach(2,&nout);
			fprintf(nout,"%*.*s\n",lline,lline,line);
		}

	}
}


/* -----------------------------------------------------------------------
    Purpose:    Center a character string.

    Usage:      CALL C1NTER (MIDDLE, LLINE, LINE)

    Arguments:
       MIDDLE - Middle point on page for line.  (Input)
       LLINE  - Length of LINE.  (Input/Output)
       LINE   - Character string containing the line.  (Input/Output)
  -----------------------------------------------------------------------
 */
void    imsl_c1nter(middle, lline, line)
	Mint             middle, *lline;
	Mchar           *line;
{
	Mint             i, left;

	if (*lline != 0)  {
		left = middle - *lline / 2;
		if (left > 0) {
			for (i = *lline;i>=1;i--) line[left+i-1]=line[i-1];
			for (i=0; i<left;i++) line[i]=' ';
			*lline += left;
		}
	}
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Get the row numbers for a column to be printed.

    Usage:      CALL W12RL (ITRING, IROW1, IROW2, JCOL, IR1, NROW)

    Arguments:
       ITRING - See WRRRL.  (Input)
       IROW1  - First row of matrix to be printed.  (Input)
       IROW2  - Last row of matrix to be printed.  (Input)
       JCOL   - Column of matrix to be printed.  (Input)
       IR1    - First row of matrix containing a number to be printed.
                (Output)
       NROW   - Number of elements of column JCOL to be printed.
                (Output)
  -----------------------------------------------------------------------
 */
void imsl_w12rl(itring, irow1, irow2, jcol, ir1, nrow)
	Mint             itring, irow1, irow2, jcol, *ir1, *nrow;
{
	Mint             ir2;

	if (itring == 0) {
		*ir1 = irow1;
		ir2 = irow2;
	} else if (itring == 1) {
		*ir1 = irow1;
		ir2 = min(jcol, irow2);
	} else if (itring == 2) {
		*ir1 = irow1;
		ir2 = min(jcol - 1, irow2);
	} else if (itring == -1) {
		*ir1 = max(jcol, irow1);
		ir2 = irow2;
	} else if (itring == -2) {
		*ir1 = max(jcol + 1, irow2);
		ir2 = irow2;
	}
	*nrow = ir2 - *ir1 + 1;
	if (*nrow < 0)  *nrow = 0;
	return;
}


Mint    imsl_write_initialize(ipagew, ipagel, icentr, large, ipgopt, ititle, rlabel,
                       clabel, iclab, irlab, ncmat, maxrlw, title, ltitle, ntitle,
                       mode, maxfw, lrowlb, nrowlb, ihotzr, nrmat, maxwl)
Mint     *ipagew, *ipagel, *icentr, *ipgopt, *ititle, *iclab, *irlab, ncmat, *maxrlw,
        *ltitle, *ntitle, mode, maxfw, *lrowlb, *nrowlb, *large, *ihotzr, nrmat, *maxwl;
Mchar    **rlabel, **clabel, *title;
{
        Mint     ibeg, iend, done, minpw;

	imsl_page(1, ipagew);             /* Get page width */
	imsl_page(2, ipagel);	        /* and length. */
	imsl_w1opt(1, icentr);     /* Get centering option. */
	imsl_w1opt(2, large);      /* Option for printing large matrices.*/
	imsl_w1opt(3, ipgopt);     /* Paging option. */
	if (*ipgopt > 0 && *ipgopt >= *ipagel)  *ipgopt = 0;

        imsl_w1opt(5, ititle);  /* Title-on-each-page option. */

	if (ncmat == 1) {*maxrlw = *ipagew*2/3;
	} else *maxrlw = *ipagew/4;

	*ihotzr = *maxrlw/3;
	                                    /*Get width of title, LTITLE, and number of
                                            lines for title, NTITLE.*/
	if (!strlen(title)) {
		*ntitle = 0;
		*ltitle = 0;
	} else {
		*ntitle = 0;
		l_w4rrl(title, *ipagew, 0, 10, ntitle, ltitle, &ibeg,
			&iend, &done);
	}
	                                    /*Get the width of the widest
                                            row label, LROWLB. Get the number of
	                                    lines in the row label requiring the
                                            most lines, NROWLB.*/

        if (mode == 1){*maxwl = min(max(maxfw+maxfw/2,15),40);
        } else *maxwl=min(max(maxfw+maxfw/2,15),83);

        l_w10rl(*irlab, *iclab, rlabel, clabel, *maxrlw, *ihotzr, NDENTR, 1,
                nrmat, lrowlb, nrowlb);

	/* Minimum page width is determined from the longest row label, the
	 * longest permitted column label, the format specification, spacing
	 * between columns, and the maximum stagger for long rows.*/

	minpw = *lrowlb + NSPACE + *maxwl + 4;
	if (*ipagew < minpw) {
		imsl_e1sti(1, *ipagew);
		imsl_e1sti(2, minpw);
		if (maxfw >= 12) {
			/* A page width of %(i1) was retrieved via IMSL	 */
			/* routine imsl_page.  The format, fmt, contains */
			/* a field width of %(i3).  Either the largest   */
			/* field width in FMT must be decreased or the   */
			/* page width increased to at least %(i2) by a   */
			/* call to imsl_page.				 */
			imsl_e1sti(3, maxfw);
			imsl_ermes(IMSL_TERMINAL, IMSL_FORMAT_TOO_WIDE);
		} else {
			/* A page width of %(i1) was retrieved via IMSL  */
			/* routine imsl_page.  The page width must be    */
			/* increased to at least %(i2) by a call to      */
			/* imsl_page.					 */
			imsl_ermes(IMSL_TERMINAL, IMSL_NARROW_PAGE);
		}
		return(1);
	}
	if (*ipagew > 255) {
			/* A page width of %(i1) was retrieved via IMSL */
			/* routine imsl_page.  The page width must be	*/
			/* reset to no greater than 255 by a call to	*/
			/* imsl_page.					*/
		imsl_e1sti(1, *ipagew);
		imsl_ermes(IMSL_TERMINAL, IMSL_WIDE_PAGE);
		return(1);
	}
        return (0);
}

void  imsl_write_controller(ido, large, iclab, irlab, clabel, rlabel, maxrlw, ihotzr,
                title, tit, nstag, lrowlb, ldata, jrow1, jrow2,
                jcol1, jcol2, nrmat, ncmat, ipagew, ipgopt,
                ntitle, ipagel, nrowlb, itring, ititle, fmt, ltitle, mode)

Mint         *ido, large, iclab, irlab, maxrlw, ihotzr, *nstag, lrowlb, *ldata,
            *jrow1, *jrow2, *jcol1, *jcol2, nrmat, ncmat, mode,
            ipagew, ipgopt, ntitle, *ipagel, nrowlb, itring, ititle, *ltitle;
Mchar        **clabel, **rlabel, *title, **tit, *fmt;

{
        static Mint nclab0, ibeg, iend, done,  nclab, nclab1, ndata, length,
                   nrowl1, icol1, icol2, irow1, irow2, istag, icol,
		   nwidth[2], ptitle, irow, empty, len1, nlrem,
                   nline, nlines, nrlab, irow11, left, lfield;
        static  Mchar *nul="";
        Mchar    *fmti[2];

        switch(*ido){
            case 0: goto L_5;
            case 70: goto L_70;
            case 120: goto L_120;
            case 200: goto L_200;
        }
L_5:    if (large == 0) {
		/* Wrapping rows is requested. Get the number of lines
		 * required for printing a row of MAT, NDATA. Get the width
		 * of the longest line, LDATA. Get the maximum indentation
		 * for staggered lines, NSTAG. Get the number of lines
		 * required for column labels.*/
		nclab0 = 0;
		if (iclab == 2) l_w4rrl(*clabel, maxrlw, 0, ihotzr, &nclab0,
					&length, &ibeg, &iend, &done);
		*nstag = 0;
L_10:
		istag = 0;
		nclab = 0;
		nclab1 = 0;
		ndata = 1;
		if (lrowlb == 0) {length = -NSPACE;
		} else length = lrowlb;

		*ldata = lrowlb;
		if (nrmat == 0 && iclab == 0)  goto L_30;
		for (icol = 1; icol <= ncmat; icol++) {
			imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth, &lfield,
                                fmti, &nline, mode);
			len1 = length + NSPACE + lfield;
			if (len1 > ipagew - *nstag + istag) {
				ndata += 1;
				istag += 1;
				if (istag == 5)
					istag = 0;
				if (istag > *nstag) {
					*nstag = istag;
					goto L_10;
				}
				*ldata = max(*ldata, length);
				nclab += nclab1;
				length = lrowlb + istag + NSPACE + lfield;
			} else {
				length = len1;
				nclab1 = max(nline, nclab1);
			}
		}
L_30:           *ldata = max(*ldata, length);
		nclab = max(nclab0, nclab + nclab1);
		if (ipgopt == -2) {
			    /*No paging is requested. Print entire matrix,
			     wrapping rows, and no paging.*/
                        *ido = 0;
                        *tit = title;
                        *jrow1 = 1;
                        *jrow2 = nrmat;
                        *jcol1 = 1;
                        *jcol2 = ncmat;
                        return;
		} else {                /* Paging is requested. */
			irow1 = 1;
			ptitle = ITRUE;
			if (ipgopt <= 0)   imsl_write_line(-1, " ");
			if (ipgopt == -1) {nlines = ntitle + 1;
			} else 	nlines = ipgopt + ntitle + 1;
			nlrem = *ipagel - nlines - nclab;
		                        /* Get number of lines required
                                        for first row label.*/
        		if (irlab <= 1 || irlab == 3) {nrowl1 = 1;
			} else {
				nrowl1 = 0;
				l_w4rrl(*rlabel, maxrlw, NDENTR, ihotzr,
                                        &nrowl1, &length, &ibeg, &iend,
					&done);
			}
			/* If current page length will not hold title,
			 * column labels and the longest row of MAT, multiply
			 * pagelength by an integer, so that page breaks
			 * will occur periodically.*/
			*ipagel *= ((ntitle+nclab+
                                    max(nrowlb, ndata))/(*ipagel))+1;
			if (*ipagel-nlines-nclab<max(nrowl1, ndata)) {
				    /*Not enough room on first page to print
				    first row of A. Begin a new page. */
				imsl_write_line(-1, " ");
				nlines = ntitle + 1;
			}
			 /* Find the number of rows of the data matrix that
			 * can be printed on the remaining lines of the
			 * page.*/
        L_40:		nlines += nclab;
			nlrem = *ipagel - nlines;
			if (irlab <= 1 || irlab == 3) {
				irow2 = min(irow1+nlrem/ndata-1,nrmat);
				nlines += (irow2-irow1+1)*ndata;
			} else {
				irow2 = 0;
				for (irow = irow1; irow <= nrmat; irow++) {
					nrlab = 0;
					l_w4rrl(rlabel[irow-1], maxrlw,
                                                NDENTR, ihotzr, &nrlab,
                                                &length, &ibeg, &iend,
						&done);
					nlrem -= max(nrlab, ndata);
					if (nlrem < 0)
						goto L_60;
					nlines += max(nrlab, ndata);
					irow2 = irow;
				}
			}
	L_60:           l_w11rl(itring, irow1, irow2, 1, ncmat, &empty,
                                jrow1, jrow2, jcol1, jcol2);
			if (empty)   {
                            ido = 0;
                            *jrow1 = -1;
                            return;
                        }
						/* Print rows JROW1 to JROW2. */
			if (ptitle) {*tit = title;
			} else  *tit = nul;
                        *ido = 70;
                        return;
	L_70:           if (irow2 < nrmat) {
				irow1 = irow2 + 1;
				if (ititle == 0) {
					ptitle = IFALSE;
					if (*ltitle < *ldata)
						*ltitle = 1;
					nlines = 1;
				} else 	nlines = ntitle + 1;
						/* Begin a new page. */
				imsl_write_line(-1, " ");
				goto L_40;
			}
                        *ido = 0;
                        *jrow1 = -1;
                        return;
		}
	} else {		/* No wrapping is requested. */
		if (ipgopt == -2) {
			irow1 = 1;          /* No paging is requested. */
			icol1 = 1;
			*nstag = 0;
			ptitle = ITRUE;
	L_80:           irow2 = min(irow1 + large - 1, nrmat);
	L_90:           if (lrowlb == 0) {length = -NSPACE;
			} else length = lrowlb;
			icol2 = 0;
			*ldata = lrowlb;
			if (nrmat == 0 && iclab == 0) goto L_110;
			for (icol = icol1; icol <= ncmat; icol++) {
				imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth,
                                        &lfield, fmti, &nline, mode);
				length += NSPACE + lfield;
				if (length > ipagew)
					goto L_110;
				*ldata = length;
				icol2 = icol;
			}
	L_110:		l_w11rl(itring, irow1, irow2, icol1, icol2,
                                &empty, jrow1, jrow2, jcol1, jcol2);
			if (empty)  {
                            ido = 0;
                            *jrow1 = -1;
                            return;
                        }
			if (ptitle) {
                            *tit = title;
	                    ptitle = IFALSE;
			    if (*ltitle < *ldata) *ltitle = 1;
                        } else *tit=nul;
                        *ido = 120;
                        return;
	L_120:          if (nrmat == 0 && iclab == 0) {
			} else if (icol2 < ncmat) {
				icol1 = icol2 + 1;
				imsl_write_line(1, " ");
				goto L_90;
			} else if (irow2 < nrmat) {
				irow1 = irow2 + 1;
				icol1 = 1;
				imsl_write_line(1, " ");
				goto L_80;
			}
		} else {		/* Paging is requested. */
			left = large;
			ndata = 1;
			irow1 = 1;
			irow11 = 1;
			icol1 = 1;
			ptitle = ITRUE;
			if (ipgopt <= 0) imsl_write_line(-1, " ");
			if (ipgopt == -1) {nlines = ntitle + 1;
			} else nlines = ipgopt + ntitle + 1;

	L_130:          if (lrowlb == 0) {length = -NSPACE;
			} else length = lrowlb;
			*ldata = lrowlb;
			icol2 = 0;
			/* Get column ICOL2 so that columns ICOL1 through
			 * ICOL2 can fit on page without wrapping.*/
			if (nrmat == 0 && iclab == 0) goto L_150;
			for (icol = icol1; icol <= ncmat; icol++) {
				imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth,
                                        &lfield, fmti, &nline, mode);
				length += NSPACE + lfield;
				if (length > ipagew) goto L_150;
				icol2 = icol;
				*ldata = length;
			}
			/* Get maximum number of lines required by column
			 * labels. */
	L_150:		if (iclab == 0) {nclab = 0;
			} else if (iclab == 1 || iclab == 3) {nclab = 1;
			} else {
                        	nclab = 0;
				l_w4rrl(*clabel, maxrlw, 0, ihotzr, &nclab,
					&length, &ibeg, &iend, &done);
				for (icol = icol1; icol <= icol2; icol++) {
					imsl_w5rrl_f(icol, iclab, clabel, fmt,
                                                nwidth, &lfield, fmti,
                                                &nline, mode);
					nclab = max(nclab, nline);
				}
			}
			/* Get number of lines required for first row label.*/
	L_170:		if (irlab <= 1 || irlab == 3) { nrowl1 = 1;
			} else {
				nrowl1 = 0;
				l_w4rrl(rlabel[irow1-1], maxrlw, NDENTR,
                                        ihotzr, &nrowl1, &length, &ibeg,
					&iend, &done);
			}
			/* If current page length will not hold title,
			 * column labels and the row IROW1 of MAT, multiply
			 * pagelength by an integer, so that page breaks
			 * will occur periodically.*/
			*nstag = 0;
			*ipagel *= ((ntitle+nclab+max(nrowl1,ndata))/(*ipagel))+1;
			if (*ipagel - nlines - nclab < max(nrowl1, ndata)) {
				/* Not enough room on page to print row
				   IROW1 of A. Begin a new page. */
				nlines = 1;
				if (ptitle)  nlines += ntitle + 1;
				imsl_write_line(-1, " ");
			}
			/* Get the number of rows that can fit on the
			 * remaining lv_part of the page. */
			nlines += nclab;
			nlrem = *ipagel - nlines;
			if (irlab <= 1 || irlab == 3) {
				irow2 = min(irow1+left-1, min(irow1+nlrem-1,
					     nrmat));
				nlines += irow2 - irow1 + 1;
			} else {
				irow2 = 0;
				for (irow = irow1; irow <= min(irow1+left - 1,
							   nrmat); irow++) {
					nrlab = 0;
					l_w4rrl(rlabel[irow-1], maxrlw,
                                                NDENTR, ihotzr, &nrlab,
                                                &length, &ibeg, &iend,
						&done);
					nlrem -= max(nrlab, ndata);
					if (nlrem < 0)
						goto L_190;
					nlines += max(nrlab, ndata);
					irow2 = irow;
				}
			}
	L_190:          l_w11rl(itring, irow1, irow2, icol1, icol2,
                                &empty, jrow1, jrow2, jcol1, jcol2);
			if (empty)  {
                            ido = 0;
                            *jrow1 = -1;
                            return;
                        }
			if (ptitle) {*tit=title;
                        }else *tit=nul;
                        *ido = 200;
                        return;
	L_200:          if (nrmat == 0 && iclab == 0) {
			} else if (irow2 < irow1 + left - 1 && irow2 != nrmat) {
				left -= irow2 - irow1 + 1;
				irow1 = irow2 + 1;
				/* Begin new page and continue printing
				 * columns ICOL1 through ICOL2.*/
				imsl_write_line(-1, " ");
				if (ititle == 0) {
					nlines = 1;
					ptitle = IFALSE;
					if (*ltitle < *ldata)  *ltitle = 1;
				} else {
					ptitle = ITRUE;
					nlines = ntitle + 1;
				}
				goto L_170;
			} else if (icol2 < ncmat) {
					/* Begin next set of columns. */
				irow1 = irow11;
				icol1 = icol2 + 1;
				left = large;
				if (nlines == *ipagel) {
					imsl_write_line(-1, " ");
					if (ititle == 0) {
						nlines = 1;
						ptitle = IFALSE;
						if (*ltitle < *ldata)
							*ltitle = 1;
					} else {
						nlines = ntitle + 1;
						ptitle = ITRUE;
					}
				} else {
					imsl_write_line(1, " ");
					nlines += 1;
					ptitle = IFALSE;
					if (*ltitle < *ldata)
						*ltitle = 1;
				}
				goto L_130;
			} else if (irow2 < nrmat) {
						/* Begin new set of rows. */
				irow1 = irow2 + 1;
				irow11 = irow1;
				left = large;
				icol1 = 1;
				if (nlines == *ipagel) {
					imsl_write_line(-1, " ");
					if (ititle == 0) {
						nlines = 1;
						ptitle = IFALSE;
						if (*ltitle < *ldata)
							*ltitle = 1;
					} else {
						nlines = ntitle + 1;
						ptitle = ITRUE;
					}
				} else {
					imsl_write_line(1, " ");
					nlines += 1;
					ptitle = IFALSE;
					if (*ltitle < *ldata)  *ltitle = 1;
				}
				goto L_130;
			}
		}
	}
        *ido = 0;
        *jrow1 = -1;
}

#undef ITRUE
#undef IFALSE
#undef NSPACE
#undef NDENTR

