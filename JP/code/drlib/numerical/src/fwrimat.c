/*Translated by FOR_C++, v0.1, on 06/14/90 at 14:41:55 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef IMSL_MACHINE_NT
/* MAR 2002/5/3 - removed this as code doesn't compile properly when 
   DOUBLE is defined */
/*#include <windows.h>*/
#endif

#define NSPACE 2
#define	NDENTR	3
#define	BIG	1.0e7
#define	EPS	1.0e-3

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)


Mint            PROTO(imsl_c1tci_f,(Mchar *, Mint, Mint *));
static void     PROTO(l_i1mt,(Mint,Mfloat*,Mint,Mint*));
static void     PROTO(l_f3trr_f,(Mint,Mint,Mfloat*,Mint,Mint,Mchar*,Mchar*));
static Mchar*   PROTO(l_fmtrr_f,(Mint,Mint,Mfloat*,Mint,Mint,Mint,Mint));
static Mchar*   PROTO(l_w1ri,(Mint,Mint,Mint,Mfloat*,Mchar*,Mint));
static void     PROTO(l_w13rl_f,(Mint,Mint,Mint,Mint,Mint,Mfloat*,Mint,Mchar*,
                        Mchar**));
static void     PROTO(l_w3rrl_f,(Mchar*,Mfloat*,Mint,Mint,Mchar*,Mchar**,
		Mchar**,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,
		Mint,Mint,Mint,Mint,Mchar*,Mchar**));
static VA_LIST_HACK  PROTO(l_f_write_matrix,(Mchar*,Mint,Mint,Mfloat*,va_list));
static void     PROTO(l_i1max,(Mint, Mfloat*, Mint, Mint*));
static void     PROTO(l_w2rrl,(Mchar*,Mint,Mint,Mfloat[],Mint,Mint,Mchar*,
                        Mchar**,Mchar**,Mint,Mint,Mchar**));

extern Mchar	*imsl_output_string;
extern Mint	imsl_output_string_length;
extern int	imsl_return_string, imsl_write_to_console;


#ifdef ANSI
void imsl_f_write_matrix(Mchar *title, Mint nra, Mint nca, Mfloat *a, ...)
#else
void imsl_f_write_matrix(title, nra, nca, a, va_alist)
    Mchar 	*title;
    Mint	nra, nca;
    Mfloat	*a;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START (argptr, a);
    E1PSH("imsl_f_write_matrix", "imsl_d_write_matrix");

    IMSL_CALL(l_f_write_matrix(title, nra, nca, a, argptr));

    E1POP("imsl_f_write_matrix", "imsl_d_write_matrix");
    va_end (argptr);
}

#ifdef ANSI
static VA_LIST_HACK l_f_write_matrix(Mchar *title, Mint nra, Mint nca, Mfloat *a,
                            va_list argptr)
#else
static VA_LIST_HACK l_f_write_matrix(title, nra, nca, a, argptr)
    Mchar        *title;
    Mint         nra, nca;
    Mfloat        *a;
    va_list     argptr;
#endif

{
    Mint	 arg_number = 4, ner = 1, irlab= -1, iclab= -1, m;
    Mint         itring = 0, itri=0, irl=0, icl=0, code;
    Mint         cda = (nca==0?1:nca);
    Mchar        **rlabel=0;
    Mchar        **clabel=0;
    Mchar        *fmt=0;
    Mint         i = 1, ifmtop, ipgopt;
    Mint         transpose = 1;
    Mchar        local[11];
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
                   irlab = 2;
                   irl = 1;
                   arg_number++;
                   rlabel = va_arg(argptr, Mchar**);
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_NO_ROW_LABELS:
		if (!irl) {
                   irlab = 0;
                   irl = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_ROW_NUMBER_ZERO:
		if (!irl) {
                   irlab = 3;
                   irl = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_ROW_NUMBER:
		if (!irl) {
                   irlab = 1;
                   irl = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_ROW_LABEL_OPT);
                   goto END_OPTIONS;
                }
            case IMSL_COL_LABELS:
		if (!icl) {
                   icl = 1;
                   iclab = 2;
                   arg_number++;
                   clabel = va_arg(argptr, Mchar**);
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
#if 0
        /* MAR: 2002/5/3 - hashed out because of (a) compilation problems and
           (b) we don't want a console appearing ever */
	HANDLE	console_handle = GetStdHandle(STD_OUTPUT_HANDLE);

	/* If no console is available, create a default one. */
	if (!console_handle)
	    imsl_create_console();
#endif
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
	fmt = local;            /* recover format */
        imsl_w1opt(6, &ifmtop);
        if (ifmtop == 1) strcpy(fmt, "%12.6W");
        else if (ifmtop == 2) strcpy(fmt, "%12.5e");
        else strcpy(fmt, "%10.4W");
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
        imsl_f_wrrrl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                     irlab, iclab);
    else if (transpose) {
        imsl_f_m1ran(nra, cda, a, a);
        imsl_f_wrrrl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                irlab, iclab);
        imsl_f_m1ran (cda, nra, a, a);
    } else {
        imsl_f_wrrrl(title, nca, nra, a, cda, itring, fmt, rlabel, clabel,
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



/* Structured by FOR_STRUCT, v0.2, on 06/14/90 at 14:41:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  W2RRL/DW2RRL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 5, 1985

    Purpose:    Print a real rectangular matrix with a given format and
                labels.

    Usage:      CALL W2RRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, CHWK)

    Arguments:  See WRRRL.

    Chapters:   STAT/LIBRARY Utilities
                MATH/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

static void l_w2rrl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel,
             irlab, iclab, chwk)
	Mchar           *title;
	Mint             nra, nca;
	Mfloat           a[];
	Mint             lda, itring, irlab, iclab;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mchar           **chwk;
{
	Mchar            fmta[11], *itmp, *tit;
	Mint             wfmt;
	Mint             icentr, ihotzr, ipagel, ipagew, ipgopt,
                         ititle, jcol1, jcol2, jrow1, jrow2, large,
                         ldata, lenint, lrowlb, ltitle, maxfw, maxrlw,
                         maxwl, nrowlb, nsig, nstag, ntitle, n, nwidth, nochwk;
        Mint             ido;


	/* Check for terminal errors and return WFMT = .TRUE. if there is a
	 * single W format and the maximum field width, MAWFW.*/

	imsl_write_format(nra, nca, lda, itring, fmt, &maxfw, "WeEfgGdiouxX", "", &n, &wfmt);

	if (imsl_n1rty(0) > 0) goto L_9000;

        if (imsl_write_initialize(&ipagew, &ipagel, &icentr, &large, &ipgopt,
                &ititle, rlabel, clabel, &iclab, &irlab, nca, &maxrlw, title,
                &ltitle, &ntitle, 1, maxfw, &lrowlb, &nrowlb, &ihotzr, nra,
                &maxwl)) goto L_9000;

	if (wfmt==1 && n == 1) {
					/* Get one format for A. */
		itmp = strchr(fmt, (Mint)'%')+1;
	        lenint = 0;
                while (isdigit((Mint)(*(itmp+lenint)))) lenint++;
                if (lenint) {imsl_c1tci_f(itmp,lenint,&nwidth);
                } else nwidth = 10;
		itmp += lenint;
                if (*itmp == '.') {
                    itmp += 1;
                    lenint = 0;
                    while (isdigit((Mint)(*(itmp+lenint)))) lenint++;
                    if (lenint) {imsl_c1tci_f(itmp,lenint,&nsig);
                    } else nsig = 4;
                } else nsig = 4;
		strcpy(fmta, l_fmtrr_f(nra, nca, a, lda, itring, nwidth,
                        nsig));
	}
	nochwk = 1;
	if (wfmt>0 && n > 1) nochwk = 0;
        wfmt = (((wfmt==1)&&n==1)?1:0);

        ido = 0;

L_10:   imsl_write_controller(&ido, large, iclab, irlab, clabel, rlabel, maxrlw,
                ihotzr, title, &tit, &nstag, lrowlb, &ldata, &jrow1, &jrow2,
                &jcol1, &jcol2, nra, nca, ipagew, ipgopt,
                ntitle, &ipagel, nrowlb, itring, ititle, fmt, &ltitle, 1);

       if (jrow1 > 0) {
	    if (nochwk==0) l_w13rl_f(itring, jrow1, jrow2, jcol1, jcol2, a, lda,
                                   fmt, chwk);
	    l_w3rrl_f(tit, a, lda, itring, fmt, rlabel, clabel, wfmt,
		    jrow1, jrow2, jcol1, jcol2, lrowlb, ldata, nstag, ltitle,
		    irlab, iclab, maxrlw, ihotzr, ipagew, icentr, fmta, chwk);
        }
        if (ido > 0)  goto L_10;

	/*
	 * Reset options back to their global state.
	 */
L_9000: return;
}
/*-----------------------------------------------------------------------
    IMSL Name:  WRRRL/DWRRRL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 5, 1985

    Purpose:    Print a real rectangular matrix with a given format and
                labels.

    Usage:      CALL WRRRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL)

    Arguments:
       TITLE  - Character string specifying the title.  (Input)
                TITLE = ' ' suppresses printing of the title.
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - NRA by NCA matrix to be printed.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement in the calling program.  (Input)
       ITRING - Triangle option.  (Input)
                ITRING   Action
                   0     Full matrix is printed.
                   1     Upper triangle of A is printed.
                   2     Upper triangle of A excluding the diagonal of A
                         is printed.
                  -1     Lower triangle of A is printed.
                  -2     Lower triangle of A excluding the diagonal of A
                         is printed.
       FMT    - Character string containing formats.  (Input)
                For example, FMT = '(F10.3)' specifies this F format for
                the entire matrix.  FMT = '(2E10.3, 3F10.3)' specifies
                an E format for columns 1 and 2 and an F format for
                columns 3, 4 and 5.  If the end of FMT is encountered and
                if some columns of the matrix remain, format control
                continues with the first format in FMT.  Even though the
                matrix A is real, an I format can be used to print the
                integer lv_part of matrix elements of A.   The most useful
                format is a special format, called "W format," that can
                be used to specify pretty formats automatically.  Set
                FMT = '(W10.4)' if you want a single D, E, F, or I
                format selected automatically with field width 10 and
                with 4 significant digits.  See Remark 4 for a general
                description of the W format.  FMT may contain only D, E,
                F, G, I, or W edit descriptors, e.g., the X descriptor is
                not allowed.  FMT must contain exactly one set of
                parentheses.
       RLABEL - CHARACTER*(*) vector of labels for rows of A.  (Input)
                If rows are to be numbered consecutively 1, 2, ..., NRA,
                use RLABEL(1) = 'NUMBER'.  If no row labels are desired,
                use RLABEL(1) = 'NONE'.  Otherwise RLABEL is a vector of
                length NRA containing the labels.
       CLABEL - CHARACTER*(*) vector of labels for columns of A.  (Input)
                If columns are to be numbered consecutively 1, 2, ...,
                NCA, use CLABEL(1)  = 'NUMBER'.  If no column labels are
                desired, use CLABEL(1) = 'NONE'.  Otherwise, CLABEL is a
                vector of length NCA + 1 containing the column headings.
                CLABEL(1) is the heading for the row labels, and for
                J = 1, 2, ..., NCA, CLABEL(1+J) is the column heading
                for the J-th column.

    Remarks:
    1. Automatic workspace is used only if W formats are used and IMSL
       routine imsl_write_options has previously been invoked with IOPT = -2 and
       ISET = 0.  In this case, workspace usage is
                WRRRL    10*NCA character units, or
                DWRRRL   10*NCA character units.
       Workspace may be explicitly provided, if desired, by use of
       W2RRL/DW2RRL.  The reference is
                CALL W2RRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, CHWK)
       The additional argument is
       CHWK   - CHARACTER*10 work vector of length NCA.  This workspace
                is referenced only if W formats are used and IMSL routine
                imsl_write_options has previously been invoked with IOPT = -2 and
                ISET = 0.  Otherwise, CHWK is not referenced and can be a
                CHARACTER*10 vector of length one.

    2. The output appears in the following form:
                           TITLE
       CLABEL(1)  CLABEL(2)  CLABEL(3)  CLABEL(4)
       RLABEL(1)      xxxxx      xxxxx      xxxxx
       RLABEL(2)      xxxxx      xxxxx      xxxxx

    3. Use '%/' within titles or labels to create a new line.  Long
       titles or labels are automatically wrapped.

    4. For printing numbers whose magnitudes are unknown, the G format
       in FORTRAN is useful; however, the decimal points will generally
       not be aligned when printing a column of numbers.  The W format is
       a special format used by this routine to select a D, E, F, or I
       format so that the decimal points will be aligned.  The W format
       is specified as '(Wn.d)'.  Here n is the field width and d is the
       number of significant digits generally printed.  Valid values for
       n are 3, 4, ..., 40.  Valid values for d are 1, 2, ..., n-2.
       If FMT specifies one format and that format is a W format, all
       elements of the matrix A are examined to determine one FORTRAN
       format for printing.  If FMT specifies more than one format,
       FORTRAN formats are generated separately from each W format.

    5. A page width of 78 characters is used.  Page width and page length
       can be reset by invoking IMSL routine PGOPT.

    6. Horizontal centering, method for printing large matrices, paging,
       method for printing NaN (not a number), and printing a title on
       each page can be selected by invoking IMSL routine imsl_write_options.

    7. Output is written to the unit specified by IMSL routine UMACH.
  -----------------------------------------------------------------------
 */
void imsl_f_wrrrl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel,
             irlab, iclab)
	Mchar           *title;
	Mint             nra, nca;
	Mfloat           a[];
	Mint             lda, itring, irlab, iclab;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
{
	Mchar            **chwk = NULL;
	Mint             wfmt;
	Mint             i, maxfw, n;


	/*
	 * Check for terminal errors and return WFMT = .TRUE. if there is a
	 * single W format and the maximum field width, MAWFW.
	 */
	imsl_write_format(nra, nca, lda, itring, fmt, &maxfw, "WeEfgGdiouxX",
                          "", &n, &wfmt);

	if (imsl_n1rty(0) != 0)  goto L_9000;

	if (wfmt>0&&n>1) {
		chwk = (Mchar **)malloc(nca*sizeof(Mchar*));
		chwk[0] = (Mchar *)malloc(11*nca*sizeof(Mchar));
		for (i=1;i<nca;i++)  chwk[i] = chwk[0] + i*11;
		l_w2rrl(title, nra, nca, a, lda, itring, fmt, rlabel,
			     clabel, irlab, iclab, chwk);
		free(chwk[0]);
		free(chwk);
		goto L_9000;
	}
	l_w2rrl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel,
                irlab, iclab, chwk);
L_9000: return;
}

/*  -----------------------------------------------------------------------
    Purpose:    Store formats and convert W formats to FORTRAN formats.

    Usage:      CALL W13RL (ITRING, IROW1, IROW2, ICOL1, ICOL2, A, LDA,
                            FMT, CHWK)

    Arguments:
       ITRING - See WRRRL.  (Input)
       IROW1  - First row of matrix to be printed.  (Input)
       IROW2  - Last row of matrix to be printed.  (Input)
       ICOL1  - First column of matrix to be printed.  (Input)
       ICOL2  - Last column of matrix to be printed.  (Input)
       A      - Matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in
                the dimension statement of the calling program.  (Input)
       FMT    - Character string containing format.  (Input)
                See WRRRL.
       CHWK   - CHARACTER*10 vector of length ICOL2-ICOL1 + 1 containing
                some output formats for columns ICOL1, ... , ICOL2 of A.
                (Output)
                A FORTRAN format is output in CHWK(I) only if the format
                in FMT for column ICOL1+I-1 is a W format.
  -----------------------------------------------------------------------
 */
static void l_w13rl_f(itring, irow1, irow2, icol1, icol2, a, lda, fmt, chwk)
	Mint             itring, irow1, irow2, icol1, icol2;
	Mfloat           a[];
Mint             lda;
Mchar           *fmt;
Mchar           *chwk[];
{
	Mchar            *fmti, *W, *N, aa='.', bb='W';
	Mint             i, icol, ir1, nrow, nsig, nwidth;

	i = 0;
	for (icol = icol1; icol <= icol2; icol++) {
		i += 1;
		fmti = imsl_w7rrl(icol, fmt);
		imsl_w6rrl(fmti, 1, &aa, &bb, &nwidth);
                N = strchr(fmti+1,(Mint)'%');
                W = strchr(fmti,(Mint)'W');
                if (!N && W) N = W+1;
		if (W && W < N) {
			            /* I-th format is W. Compute actual FORTRAN
                                       format to be used.*/
			imsl_w6rrl(fmti, 0, &aa, &bb, &nsig);
			imsl_w12rl(itring, irow1, irow2, icol, &ir1, &nrow);
			if (nrow != 0) {
			    strcpy(chwk[i-1], l_fmtrr_f(nrow, 1,
                              &a[(icol-1)*lda+ir1-1], nrow, 0, nwidth, nsig));
                        }
		}
	}
	return;
}

/*  -----------------------------------------------------------------------
    Purpose:    Write a real rectangular matrix with a specified format
                and labels.

    Usage:      CALL W3RRL (TITLE, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, WFMT, IROW1, IROW2, ICOL1, ICOL2,
                            LROWLB, LDATA, NSTAG, LTITLE, IRLAB, ICLAB,
                            MAXRLW, IHOTZR, IPAGEW, ICENTR,
                            FMTA, NOCHWK, CHWK)

    Arguments:
       TITLE  - See WRRRL.
       A      - See WRRRL.
       LDA    - See WRRRL.
       ITRING - See WRRRL.
       FMT    - See WRRRL.
       RLABEL - See WRRRL.
       CLABEL - See WRRRL.
       WFMT   - LOGICAL variable.  (Input)
                WFMT = .TRUE. means a single W format is specified.
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
                  1     Row labels are 1, 2, 3, ...
                  2     Row labels are in RLABEL.
       ICLAB  - Column label option.  (Input)
                ICLAB   Meaning
                  0     No column labels.
                  1     Column labels are 1, 2, 3, ...
                  2     Column labels are in CLABEL.
       MAXRLW - Maximum permitted row label width.  (Input)
       IHOTZR - Length of hot zone for row labels.  (Input)
                The hot zone is the area where line breaks can occur.
       IPAGEW - Page width.  (Input)
       ICENTR - Centering option.  (Input)
                ICENTR   Action
                  0      No centering
                  1      Center
       FMTA   - Character string of length 10 containing a single FORTRAN
                format.  (Input if WFMT is .TRUE.)
       NOCHWK - Indicator of storage of formats.  (Input)
                NOCHWK
                  0     Formats are input in CHWK
                  1     Formats are not input in CHWK
       CHWK   - CHARACTER*10 vector of length ICOL2-ICOL1+1 containing
                formats.  (Input if NOCHWK = 0)
                If NOCHWK = 1, CHWK is not referenced.
  -----------------------------------------------------------------------
 */
static void l_w3rrl_f(title,a,lda,itring,fmt,rlabel,clabel,wfmt,irow1,
      irow2, icol1, icol2, lrowlb, ldata, nstag, ltitle, irlab, iclab,
      maxrlw, ihotzr, ipagew, icentr, fmta, chwk)
	Mchar           *title;
	Mfloat           a[];
	Mint             lda, itring;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mint            wfmt;
	Mint             irow1, irow2, icol1, icol2, lrowlb, ldata, nstag,
			ltitle, irlab, iclab, maxrlw, ihotzr, ipagew, icentr;
	Mchar           *fmta;
	Mchar           *chwk[];
{
	Mchar            blanks[256], *fmti, line[256], *N, *W;
	Mint             i, done, icol, ipad, irow, istag, lfield, lline, lout, middle,
                         nline, nlinei, nwidth;


	for(i=0;i<256;i++) blanks[i]=' ';blanks[255] = '\0';
				/* Get width of output, LOUT. */
	lout = max(ltitle, ldata);
	if (icentr == 1) middle = ipagew / 2;
	else middle = lout / 2;

	strcpy(line,blanks);

	nline = 0;                                  /* Print title. */
	imsl_write_title(title, ipagew, &nline, line, middle);
	strcpy(line,blanks);

        imsl_write_labels (iclab, clabel, line, blanks, &nline, lrowlb, maxrlw,
                      ihotzr, middle, icol1, icol2, ipagew, nstag,
                      fmt, ldata, 0);
        strcpy(line,blanks);
	                                                    /* Print A. */
	for (irow = irow1; irow <= irow2; irow++) {
		nline = 1;
		imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr, NDENTR,
			nline, &lline, line, &done);
		if (lline > 0)  lline += NSPACE;
		istag = 0;
		for (icol = icol1; icol <= icol2; icol++) {
			/* Get format*/
			imsl_w5rrl_f(icol, iclab, clabel, fmt, &nwidth, &lfield,
                                &fmti, &nlinei, 1);
			                /* Get padding before the number. */
			ipad = lfield - nwidth;
			                /* Check to see if current line fits*/
			if (lline + lfield > ipagew - nstag + istag) {
				/* Line won't fit so print old line, and
				 * begin new line.*/
				lline = ldata;
				imsl_c1nter(middle, &lline, line);
				imsl_write_line(lline, line);
                                strcpy(line,blanks);
				nline += 1;
				imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw,
                                        ihotzr, NDENTR, nline, &lline, line, &done);
				if (lline > 0)  lline += NSPACE;
				istag += 1;
				if (istag == 5)  istag = 0;
				lline += istag;
			}
			if (wfmt) fmti = fmta;
                        N = strchr(fmti+1,(Mint)'%');
                        W = strchr(fmti,(Mint)'W');
                        if (!N && W) N = W+1;
		        if (W && W < N) fmti = chwk[icol-icol1];
			                /* Append entry A(IROW,ICOL) to LINE. */
			lline += ipad;
			strncpy(&line[lline], l_w1ri(itring, irow, icol,
				a+(icol-1)*lda+irow-1,fmti,nwidth),nwidth);
			lline += nwidth + NSPACE;
		}
		/** Print last line for IROW-th row of A.*/
		lline = ldata;
		imsl_c1nter(middle, &lline, line);
		imsl_write_line(lline, line);
                strcpy(line,blanks);
L_110:
		if (!done) {
			/*
			 * Number of lines for IROW-th row label exceeds that
			 * for the IROW-th row of data. Print the remainder
			 * of the row label.*/
			nline += 1;
			imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr,
				NDENTR, nline, &lline, line, &done);
			lline = ldata;
			imsl_c1nter(middle, &lline, line);
			imsl_write_line(lline, line);
                        strcpy(line,blanks);
			goto L_110;
		}
	}
	return;
}

/*-----------------------------------------------------------------------
    Purpose:    Convert a matrix element to a character string according
                to a supplied format.  Dots are returned if not-a-number
                is input.  If the matrix element is not in the triangle
                to be printed, blanks are returned.

    Usage:      WRITE (NOUT,1) W1RI(ITRING,I,J,X,FMT)

    Arguments:
       ITRING - Triangle option.  (Input)
                 ITRING   Action
                    0     Full matrix is printed.
                    1     Upper triangle of A is printed.
                    2     Upper triangle of A excluding the diagonal of A
                          is printed.
                   -1     Lower triangle of A is printed.
                   -2     Lower triangle of A excluding the diagonal of A
                          is printed.
       I      - Row number of matrix element.  (Input)
       J      - Column number of matrix element.  (Input)
       X      - The matrix element to be written.  (Input)
       FMT    - The format to be used in writing X.  (Input)
       W1RI   - A CHARACTER*40 string which contains blanks if the matrix
                element is not to be printed.  Otherwise, if X is
                not-a-number W1RI contains dots, and if X is not
                not-a-number W1RI contains a representation of X in the
                form specified by FMT.
                (Output)
  ----------------------------------------------------------------------
 */
static Mchar *l_w1ri(itring, i, j, x, fmt, nwidth)
	Mint             itring, i, j, nwidth;
	Mfloat           *x;
	Mchar           *fmt;
{
	static Mchar     w1ri_v[41];
	Mint             iwrite;


	iwrite = 0;
	if (itring == 0) iwrite = 1;
	else if (itring == 1) {
		if (i <= j)  iwrite = 1;
	} else if (itring == 2) {
		if (i < j) iwrite = 1;
	} else if (itring == -1) {
		if (i >= j) iwrite = 1;
	} else if (itring == -2) {
		if (i > j) iwrite = 1;
	}
	if (iwrite == 0) {
		strcpy(w1ri_v, "                                        ");
	} else {
		strcpy(w1ri_v, imsl_w1iss(x, fmt, nwidth));
	}
	return (w1ri_v);
}

/*  -----------------------------------------------------------------------
    Purpose:    Return a pretty FORTRAN format for a real rectangular
                matrix.

    Usage:      FMTRR(NRA, NCA, A, LDA, ITRING, NWIDTH, NSIG)

    Arguments:
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - NRA by NCA matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ITRING - Triangle option.  (Input)
                ITRING  Action
                   0    Full matrix is examined.
                   1    Upper triangle of A is examined.
                   2    Upper triangle of A excluding the diagonal of A
                        is examined.
                  -1    Lower triangle of A is examined.
                  -2    Lower triangle of A excluding the diagonal of A
                        is examined.
       NWIDTH - Width of FORTRAN format field.  (Input)
                NWIDTH = 10 or greater is suggested for pretty formats.
                NWIDTH must equal 3, 4, ..., or 40.
       NSIG   - Number of significant digits generally printed.  (Input)
                NSIG = NWIDTH-6 is suggested for pretty formats.  NSIG
                must not exceed NWIDTH-2.
       FMTRR  - A CHARACTER string of length 10 that contains the FORTRAN
                format.  (Output)
                FMTRR = '(1PE  .  )' or '(0PF  .  )' or '(I  )' where
                appropriate integers appear right justified in the
                blanks.

    Remarks:
    1. The format returned in FMTRR is determined by the first rule
       listed below that is satisfied.
         a.  If NWIDTH is less than 7, an F format is returned.
         b.  If all elements of A equal NaN (not a number) or NRA or NCA
             is than or equal to zero, an E (or D) format is returned.
         c.  If all elements of A not equal to NaN are integers and can
             be printed with an I format with field width NWIDTH, an I
             format is returned.
         d.  Let AMAX be the value in A with largest absolute value.
             If 0 < ABS(AMAX) < 0.001 or ABS(AMAX) is greater than or
             equal to 1.E7, an E (or D) format is returned.  The E format
             will contain MIN(NSIG,ISIG) significant digits for AMAX.
             ISIG = INT(0.05-LOG10(AMACH(4)) for FMTRR and
             ISIG = INT(0.05-LOG10(DMACH(4)) for DFMTRR.  Here the use of
             ISIG is intended to prevent printing of garbage digits.
         e.  If MIN(NSIG,ISIG) significant digits can be printed using
             an F format with field width NWIDTH, an F format is
             returned.
         f.  An E format is returned containing MIN(NSIG,ISIG).

    2. The F format returned suppresses the printing of trailing
       zeroes after the decimal point.
  -----------------------------------------------------------------------
 */


static Mchar *l_fmtrr_f(nra, nca, a, lda, itring, nwidth, nsig)
	Mint             nra, nca;
	Mfloat           a[];
Mint             lda, itring, nwidth, nsig;
{
	Mchar            temp[11], temp1[11], width[3];
	static Mchar     fmtrr_v[11];
	Mint             i, iflag, indx, ir1, j, nrow;
	Mfloat           aij, amax, itmp;


	sprintf(width, "%d", nwidth);
	for (j = 1; j <= nca; j++) {
		imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		for (i = ir1; i <= (ir1 + nrow - 1); i++) {
			aij = a[(j - 1) * lda + i - 1];
			if (!imsl_ifnan(aij))  goto L_30;
		}
	}
	/* Either all elements of A are not a number or NRA or NCA is less
	 * than one. Return an arbitrary format. */
	aij = imsl_amach(6);
	strcpy(fmtrr_v, imsl_fmtx(&aij, nwidth, nsig));
	goto L_9000;
	                        /* Determine if I format is appropriate.*/
L_30:
	amax = F_ZERO;
	for (j = 1; j <= nca; j++) {
		imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		if (nrow > 0)  {
		    l_i1mt(nrow, &a[(j - 1) * lda + ir1 - 1], nwidth, &iflag);
		    if (iflag == 0) break;
                }
	}
	if (iflag == 1) {
		                        /* Return an I format. */
		strcpy(fmtrr_v, "%");
                strcat(fmtrr_v, width);
                strcat(fmtrr_v, "d");
		goto L_9000;
	} else {
		                        /* Find largest absolute value. */
		amax = F_ZERO;
		for (j = 1; j <= nca; j++) {
			imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
			if (nrow > 0)  {
			    l_i1max(nrow, &a[(j-1)*lda+ir1-1], 1, &indx);
			    if (indx > 0){
				    itmp = a[(j-1)*lda+ir1+indx-2];
				    if (itmp < 0) itmp = -itmp;
				    amax = max(itmp,amax);
			    }
	                }
		}
		strcpy(temp, imsl_fmtx(&amax, nwidth, nsig));
		if (temp[strlen(temp)-1] == 'f') {
			/* If f format, see if there are trailing zeroes that
			 * can be stripped by modifing the f format.*/
			strcpy(temp1, temp);
			l_f3trr_f(nra, nca, a, lda, itring, temp1, temp);
		}
		strcpy(fmtrr_v, temp);
	}
L_9000: return (fmtrr_v);
}

/*  -----------------------------------------------------------------------
    Purpose:    Change an F format so that trailing zeroes are not
                printed after the decimal point for a matrix.

    Usage:      CALL F3TRR (NRA, NCA, A, LDA, ITRING, FFMTI, FFMTO)

    Arguments:
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - NRA by NCA matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ITRING - Triangle option.  (Input)
                ITRING  Action
                   0    Full matrix is used.
                   1    Upper triangle of A is used.
                   2    Upper triangle of A excluding the diagonal of A
                        is used.
                  -1    Lower triangle of A is used.
                  -2    Lower triangle of A excluding the diagonal of A
                        is used.
       FFMTI  - Character string of length 10 containing the F format.
                (Input)
                FFMTI = '( PF  .  )' where the blanks are filled with
                integers that are right justified.
       FFMTO  - Character string of length 10 containing the changed F
                format.  (Output)
                FFMTO = '( PF  .  )' where the blanks are filled with
                integers that are right justified.
  -----------------------------------------------------------------------
 */
static void l_f3trr_f(nra, nca, a, lda, itring, ffmti, ffmto)
	Mint             nra, nca;
	Mfloat           a[];
Mint             lda, itring;
Mchar           *ffmti;
Mchar           *ffmto;
{
	Mint             i, intx, ir1, j, ndec, ndec1, ndec2, nrow;
	Mfloat           rntx, aij, pwr;
        Mchar            tmp[4], *itmp;


	strcpy(ffmto, ffmti);
	itmp = strchr(ffmto,(Mint)'.');
        if (isdigit((Mint)itmp[1])) ndec = (Mint)itmp[1] - (Mint)'0';
	if (isdigit((Mint)itmp[2])) ndec = 10*ndec + (Mint)itmp[2] - (Mint)'0';
	if (ndec != 0) {
	    ndec1 = 0;
            pwr = pow((double)F_TEN,(double)ndec);
	    for (j = 1; j <= nca; j++) {
		    imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		    for (i = ir1; i <= (ir1 + nrow - 1); i++) {
			    ndec2 = ndec;
                            aij = a[(j-1)*lda+i-1];
			    if (!imsl_ifnan(aij)) {
			        rntx = aij*pwr +((a[(j-1)*lda+i-1]>F_ZERO)?F_HALF:-F_HALF);
			        intx = fmod(rntx, F_TEN);
			        if (intx != 0)  return;
	                        while (intx == 0) {
			            ndec2 -= 1;
			            if (ndec2 == 0)  break;
			            rntx /= F_TEN;
			            intx = fmod(rntx, F_TEN);
			        }
			        ndec1 = max(ndec1, ndec2);
	                    }
		    }
	    }
            sprintf(tmp,"%df",ndec1);
            strcpy(itmp+1,tmp);
        }
	return;
}


/*  -----------------------------------------------------------------------
    Purpose:    Return a pretty FORTRAN format for a scalar.

    Usage:      FMTX(X, NWIDTH, NSIG)

    Arguments:
       X      - Number.  (Input)
       NWIDTH - Width of FORMAT field.  (Input)
                NWIDTH = 10 or greater is suggested for pretty formats.
                NWIDTH must equal 3, 4, ..., or 40.
       NSIG   - Number of significant digits generally printed.  (Input)
                NSIG may not exceed NWIDTH-2.
       FMTX   - CHARACTER string of length 10 that contains the FORMAT.
                (Output)
                FMTX = '(1PE  .  )' or '(0PF  .  )' where the blanks are
                filled in with appropriate integers.  See Remark.

    Remark:
       The format returned in FMTX is determined by the first rule listed
       below that is satisfied.
         a. If NWIDTH is less than 7, an F format is returned in FMTX.
         b. If X is NaN (not a number), an E (or D) format is returned in
            FMTX.
         c. If 0 < ABS(X) < 0.001 or ABS(X) is greater than or equal to
            1.E7 an E (or D) format is returned.  The E format will
            contain MIN(NWIDTH-6,ISIG,NSIG) significant digits for X.
            ISIG = INT(0.05-LOG10(AMACH(4)) for FMTX and
            ISIG = INT(0.05-LOG10(DMACH(4)) for DFMTX.  Here the use of
            ISIG is intended to prevent printing of garbage digits.
         d. If MIN(ISIG,NSIG) significant digits can be printed using an
            F format with field width NWIDTH, an F format is returned.
         e. An E format is returned containing MIN(NWIDTH-6,ISIG,NSIG)
            significant figures.
  -----------------------------------------------------------------------
 */

Mchar *imsl_fmtx(x, nwidth, nsig)
	Mfloat           *x;
	Mint             nwidth, nsig;
{
	Mchar            drdec[11], width[11];
	static Mchar     fmtx_v[11];
	Mint             ie, isig, itemp1, itemp2, nldec, nrdec, nsig1;
	Mfloat           absx, xlog10;

	sprintf(width, "%d", nwidth);
	ie = 0;
	/* Set number of significant digits on the basis of AMACH(4)*/
	isig = .05 - log10(imsl_amach(4));
	nsig1 = min(isig, nsig);
	if (imsl_ifnan(*x)) {
		if (nwidth < 7) {
		                	/* Return an arbitrary F format. */
                        strcpy(fmtx_v,"%");
                        strcat(fmtx_v,width);
                        strcat(fmtx_v,"f");
			goto L_9000;
		} else {
			                /* Set E format indicator. */
			ie = 1;
			goto L_20;
		}
	}
	absx = fabs(*x);
	if (nwidth >= 7) {
		/*
		 * When numbers are very large or small use an E FORMAT.
		 */
		if ((absx > F_ZERO) && (absx < EPS))  ie = 1;
		if (absx >= BIG)  ie = 1;
		if (ie == 1)  goto L_20;
	}
	if (ie == 0) {
		/* F FORMAT */
		if (absx == F_ZERO) {
			        /*Set number of digits to right of decimal*/
			nrdec = nsig1;
			goto L_20;
		} else {
			if (absx >= F_ONE) {
				itemp1 = (Mint) (log10(absx));
				nldec = itemp1 + 1;
				nrdec = nsig1 - nldec;
				if (nrdec > 0) {
					/*
					 * Check to see when number is
					 * printed with NRDEC digits to the
					 * right if rounding produces another
					 * digit at the left.
					 */
					itemp2 = (Mint) (log10(absx + F_HALF *
                                                   pow((double)F_TEN, (double)(-nrdec))));
					if (itemp2 > itemp1) {
						nrdec -= 1;
						nldec += 1;
					}
				}
			} else {
				nldec = 1;
				xlog10 = -log10(absx);
				itemp1 = (Mint) (xlog10);
				nrdec = nsig1 + itemp1;
				if ((Mfloat) (itemp1) == xlog10) {
				        	/* NUMBER IS A POWER OF 10 */
					nrdec -= 1;
				} else {
					/*
					 * Check to see when number is
					 * printed with NRDEC digits to the
					 * right if rounding produces another
					 * digit at the left.
					 */
					xlog10 = log10(absx + F_HALF * pow((double)F_TEN, (double)(-nrdec)));
					if (xlog10 > F_ZERO) {
						/* When rounded the number is
						 * 1.0*/
						nrdec -= 1;
					} else {
						itemp2 = (Mint) (-xlog10);
						if (itemp2 < itemp1)
							nrdec -= 1;
					}
				}
			}
			if (nrdec < 0)  nrdec = 0;
	L_10:
			if (nldec + 2 + nrdec > nwidth) {
				/*
				 * Number does not fit in the specified F
				 * format
				 */
				if (nwidth >= 7) {
					ie = 1;
				} else {
					if (nrdec > 0) {
						nrdec -= 1;
					} else {
						nldec -= 1;
					}
					goto L_10;
				}
			}
		}
	}
L_20:
	if (ie == 1) {
		nrdec = min(nsig1 - 1, nwidth - 7);
		sprintf(drdec, "%d", nrdec);
                strcpy(fmtx_v, "%");
                strcat(fmtx_v, width);
                strcat(fmtx_v, ".");
                strcat(fmtx_v, drdec);
                strcat(fmtx_v, "e");
	} else {
		sprintf(drdec, "%d", nrdec);
                strcpy(fmtx_v, "%");
                strcat(fmtx_v, width);
                strcat(fmtx_v, ".");
                strcat(fmtx_v, drdec);
                strcat(fmtx_v, "f");
	}
L_9000:
	return (fmtx_v);
}				/* end of function */

/*  -----------------------------------------------------------------------
    Purpose:    Indicate if an I format for printing is appropriate.

    Usage:      CALL I1MT (N, X, NWIDTH, IFLAG)

    Arguments:
       N      - Length of X.  (Input)
       X      - Vector of length N.  (Input)
       NWIDTH - Field width.  (Input)
       IFLAG  - Indicator for I format.  (Output)
                IFLAG = 1 means an I format is appropriate because (1)
                all elements of X (not equal to NaN) have no fractional
                parts and (2) the field width NWIDTH is large enough for
                printing the elements of X.  Otherwise, IFLAG = 0.
  -----------------------------------------------------------------------
 */
static void l_i1mt(n, x, nwidth, iflag)
	Mint             n;
	Mfloat           x[];
Mint             nwidth, *iflag;
{
	Mint             i;
	Mfloat           xi, xint, xint1, xint2;


	*iflag = 0;
	xint1 = INT_MIN;
	xint2 = INT_MIN;
	if (nwidth < log10(imsl_amach(2))) {
		xint1 = -pow((double)F_TEN,(double)(nwidth-1));
		xint2 = -F_TEN * xint1;
	}
	for (i = 1; i <= n; i++) {
		xi = x[i - 1];
		if (!imsl_ifnan(xi)) {
		    if (fabs(xi) > INT_MAX) goto L_9000;
		    xint = (Mfloat) ((Mint) (xi));
		    if (xint != xi) goto L_9000;
		    if (xint <= xint1) goto L_9000;
		    if (xint >= xint2) goto L_9000;
		}
	}
	*iflag = 1;
L_9000: return;
}

/*  -----------------------------------------------------------------------
    Purpose:    Find the smallest index of the component of a real vector
                x (possibly containing NaN) having maximum absolute
                value.

    Usage:      CALL I1MAX (N, X, INCX, INDX)

    Arguments:
       N      - Length of vector x.  (Input)
                The vector x refers to a specific subvector obtained from
                the vector X.  See the description of the argument INCX.
       X      - Vector of length N*INCX.  (Input)
       INCX   - Increment between the referenced elements of X.  (Input)
                X(1+((I-1)*INCX)) (which is defined to be x(I)) is
                referenced, for I = 1, 2, ... , N.  INCX must be
                positive.
       INDX   - The smallest index of the component of x having maximum
                absolute value.  (Output)
                Components of x equal to NaN are excluded from
                computations.  If N = 0 or all elements of x equal NaN,
                INDX = 0.
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_i1max(Mint n, Mfloat *x, Mint incx, Mint *indx)
#else
static void l_i1max(n, x, incx, indx)
    Mint             n;
    Mfloat           x[];
    Mint             incx, *indx;
#endif
{
	Mint             i, j;
	Mfloat           absxi, xi, xmax;


	*indx = 0;
	j = 1;
	for (i = 1; i <= n; i++) {
		xi = x[j-1];
		if (!imsl_ifnan(xi)) {
			absxi = fabs(xi);
			if (*indx == 0) {
				xmax = absxi;
				*indx = i;
			} else {
				if (absxi > xmax) {
					xmax = absxi;
					*indx = i;
				}
			}
		}
		j += incx;
	}
	return;
}
/*  -----------------------------------------------------------------------
    Purpose:    Convert a number to a character string according to a
                supplied format.  Dots or blanks are returned if NaN
                (not a number) is used.

    Usage:      W1ISS(X, FMT)

    Arguments:
       X      - The number to be written.  (Input)
       FMT    - The format to be used in writing X.  (Input)
       W1ISS  - A CHARACTER*40 string which contains a character
                representation of X in the form specified by FMT.
                (Output)
                If X equals NaN,
                W1ISS = '........................................'.
                (IMSL subroutine imsl_write_options can be invoked to change the
                character representation of NaN.)
                If X equals positive machine infinity,
                W1ISS = '++++++++++++++++++++++++++++++++++++++++'.
                If X equals negative machine infinity,
                W1ISS = '----------------------------------------'.
  -----------------------------------------------------------------------
 */
Mchar *imsl_w1iss(x, fmt, nwidth)
	Mfloat          *x;
	Mchar           *fmt;
        Mint            nwidth;
{
	static Mchar     w1iss_v[41], fmt1[10], *iend, *ibeg;
        Mchar            *stars="****************************************";
	Mint             inan, ix;
        Mfloat           y;

                                    /* Separate out the format */
        iend = strchr(fmt+1,(Mint)'%');
        if (iend) {
                strncpy(fmt1,fmt,(Mint)(iend-fmt));
                fmt1[(Mint)(iend-fmt)]='\0';
        } else strcpy(fmt1,fmt);
	if (imsl_ifnan(*x)) {
		                    /* Retrieve option for printing NaN */
		imsl_w1opt(4, &inan);
		if (inan == 0) {
			strcpy(w1iss_v, "........................................");
		} else {
			strcpy(w1iss_v, "                                        ");
		}
	} else if (*x == imsl_amach(7)) {
		strcpy(w1iss_v, "++++++++++++++++++++++++++++++++++++++++");
	} else if (*x == imsl_amach(8)) {
		strcpy(w1iss_v, "----------------------------------------");
	} else {
            iend = fmt1+strcspn(fmt1,"diouxX");
            ibeg = strchr(fmt1+1,(Mint)'%');
            if (*iend!='\0' && (!ibeg || iend < ibeg)) {
                y = INT_MAX;
                if (*x <= y && *x >= -y) {
	            ix = *x;
	            sprintf(w1iss_v,fmt1,ix);
                    if (strlen(w1iss_v) > nwidth) strncpy(w1iss_v, stars,
                                                            nwidth);
                    return(w1iss_v);
		} else strncpy(w1iss_v, stars, nwidth);
	    } else {
                sprintf(w1iss_v,fmt1,*x);
                if (strlen(w1iss_v) > nwidth)  strncpy(w1iss_v, stars,
                                                    nwidth);
                }
        }
	return (w1iss_v);
}
#undef NSPACE
#undef NDENTR
#undef BIG
#undef EPS
