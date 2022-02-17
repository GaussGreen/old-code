 /*Translated by FOR_C++, v0.1, on 06/19/90 at 15:57:48 */
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

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)

#define A(I_,J_)	(a+(I_)*(lda)+(J_))


Mint            PROTO(imsl_c1tci_f,(Mchar *, Mint, Mint *));
static void     PROTO(l_i1cmx,(Mint,Mf_complex[],Mint,Mint,Mint*));
static void     PROTO(l_c1fmt,(Mint,Mf_complex*,Mint,Mint,Mint*));
static Mchar*   PROTO(l_c1tri,(Mint,Mint,Mint,Mf_complex[],Mint*,Mchar**));
static Mchar*   PROTO(l_f1tcr_f,(Mint,Mint,Mf_complex[],Mint,Mint,Mint,Mint,
                                    Mint));
static void     PROTO(l_f3tcr_f,(Mint,Mint,Mf_complex[],Mint,Mint,Mint,Mchar*,
                                    Mchar*));
static void     PROTO(l_w4crl_f,(Mint,Mint,Mint,Mint,Mint,Mf_complex[],Mint,
                                    Mchar*,Mchar**));
static void     PROTO(l_w3crl_f ,(Mchar*,Mf_complex[],Mint,Mint,Mchar*,Mchar**,
                                    Mchar**,Mint,Mint,Mint,Mint,Mint,Mint,
                                    Mint,Mint,Mint,Mint,Mint,Mint,Mint,Mint,
                                    Mint,Mchar*,Mchar**));
static VA_LIST_HACK  PROTO(l_c_write_matrix,(Mchar*,Mint,Mint,Mf_complex*,va_list));
static void     PROTO(l_w2crl,(Mchar*,Mint,Mint,Mf_complex[],Mint,Mint,Mchar*,
                        Mchar**,Mchar**,Mint,Mint,Mchar**));

extern Mchar	*imsl_output_string;
extern Mint	imsl_output_string_length;
extern int	imsl_return_string, imsl_write_to_console;


#ifdef ANSI
void imsl_c_write_matrix(Mchar *title, Mint nra, Mint nca, Mf_complex a[], ...)
#else
void imsl_c_write_matrix(title, nra, nca, a, va_alist)
    Mchar        *title;
    Mint         nra, nca;
    Mf_complex    *a;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START (argptr, a);
    E1PSH("imsl_c_write_matrix", "imsl_z_write_matrix");

    IMSL_CALL(l_c_write_matrix(title, nra, nca, a, argptr));

    E1POP("imsl_c_write_matrix", "imsl_z_write_matrix");
    va_end (argptr);
}

#ifdef ANSI
static VA_LIST_HACK l_c_write_matrix(Mchar *title, Mint nra, Mint nca,
                                Mf_complex *a, va_list argptr)
#else
static VA_LIST_HACK l_c_write_matrix(title, nra, nca, a, argptr)
    Mchar        *title;
    Mint         nra, nca;
    Mf_complex    *a;
    va_list     argptr;
#endif

{
    Mint         arg_number  = 4, ner = 1, irlab = -1, iclab = -1, m;
    Mint         itring = 0, itri = 0, irl = 0, icl = 0, code;
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
                   irl = 1;
                   irlab = 2;
                   arg_number++;
                   rlabel = va_arg(argptr, Mchar**);
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
            case IMSL_ROW_NUMBER:
		if (!irl) {
                   irl = 1;
                   irlab = 1;
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
            case IMSL_COL_NUMBER:
		if (!icl) {
                   icl = 1;
                   iclab = 1;
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
    imsl_c1iarg(nra, "nra", 0, -1, &ner);
    if (cda != nca) imsl_c1iarg(cda, "a_col_dim", 1, -1, &ner);
    else ner += 1;
    if (cda >= 1 && nca != 0) {
	imsl_c12ile(nca, "nca", cda, "a_col_dim", &ner);
    } else {
	ner += 1;
    }
    imsl_c1iarg(nca, "nca", 0, -1, &ner);
    imsl_null_pointer("title", -1, (void *)title);
    if (nra != 0 && nca != 0) imsl_null_pointer("a", -1, (void *)a);
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
    if (iclab==-1) {
    	if (nca==1) iclab = 0;
    	else iclab = 1;
    }
    if (nra == 0 || nca == 0)
        imsl_c_wrcrl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                     irlab, iclab);
    else if (transpose) {
        imsl_c_m1ran (nra, cda, a, a);
        imsl_c_wrcrl(title, nra, nca, a, nra, itring, fmt, rlabel, clabel,
                irlab, iclab);
        imsl_c_m1ran (cda, nra, a, a);
    } else {
        imsl_c_wrcrl(title, nca, nra, a, cda, itring, fmt, rlabel, clabel,
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



/* Structured by FOR_STRUCT, v0.2, on 06/19/90 at 15:57:42
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  W2CRL/DW2CRL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 5, 1985

    Purpose:    Print a complex rectangular matrix with a given format
                and labels.

    Usage:      CALL W2CRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, CHWK)

    Arguments:  See WRCRL.

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

static void l_w2crl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel,
                irlab, iclab, chwk)
	Mchar           *title;
	Mint             nra, nca;
	Mf_complex       *a;
	Mint             lda, itring, irlab, iclab;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mchar           **chwk;
{

	Mchar            fmta[11], *itmp, *tit;
	Mint             wfmt, n;
	Mint             icentr, ihotzr, ipagel, ipagew, ipgopt,
                         ititle, jcol1, jcol2, jrow1, jrow2, large,
                        ldata, lenint, lrowlb, ltitle, maxfw, maxrlw,
                        maxwl, nrowlb, nsig, nstag, ntitle, nwidth[2],
                        nochwk, ido;



	/* Check for terminal errors and return WFMT = .TRUE. if there is a
	 * single W format and the maximum field width, MAWFW.*/

        imsl_write_format(nra, nca, lda, itring, fmt, &maxfw,
                "diouxXWeEfgG", "", &n, &wfmt);
	if (imsl_n1rty(0) != 0) goto L_9000;

        if (imsl_write_initialize(&ipagew, &ipagel, &icentr, &large, &ipgopt,
                &ititle, rlabel, clabel, &iclab, &irlab, nca, &maxrlw, title,
                &ltitle, &ntitle, 2, maxfw, &lrowlb, &nrowlb, &ihotzr, nra,
                &maxwl)) goto L_9000;
	if (wfmt==1 && n==1) {
		                                  /* Get one format for A. */
		itmp = strchr(fmt, (Mint)'%')+1;
                lenint = 0;
                while (isdigit((Mint)(*(itmp+lenint)))) lenint++;
                if (lenint) {imsl_c1tci_f(itmp,lenint,&nwidth[0]);
                } else nwidth[0] = 10;
                itmp += lenint;
                if (*itmp == '.') {
                    itmp += 1;
                    lenint = 0;
                    while (isdigit((Mint)(*(itmp+lenint)))) lenint++;
                    if (lenint)  {imsl_c1tci_f(itmp,lenint,&nsig);
                    } else nsig = 4;
                } else nsig = 4;
		nwidth[1] = nwidth[0];
		strcpy(fmta, l_f1tcr_f(nra, nca, a, lda, itring, 2,
                       nwidth[0], nsig));
	}
	nochwk = 1;
	if (wfmt>0 && n>1) nochwk = 0;
        wfmt = (((wfmt==1)&&n==1)?1:0);

        ido = 0;
L_10:   imsl_write_controller(&ido, large, iclab, irlab, clabel, rlabel, maxrlw,
                ihotzr, title, &tit, &nstag, lrowlb, &ldata, &jrow1, &jrow2,
                &jcol1, &jcol2, nra, nca, ipagew, ipgopt, ntitle, &ipagel,
                nrowlb, itring, ititle, fmt, &ltitle, 2);

       if (jrow1 > 0) {
	    if (nochwk == 0)
	       l_w4crl_f(itring, jrow1, jrow2, jcol1, jcol2, a, lda,
			fmt, chwk);
	       l_w3crl_f(tit, a, lda, itring, fmt, rlabel, clabel,
		    wfmt, jrow1, jrow2, jcol1, jcol2, lrowlb, ldata, nstag,
		    ltitle, irlab, iclab, maxrlw, ihotzr, ipagew,
		    icentr, fmta, chwk);
        }
        if (ido > 0)  goto L_10;



L_9000: return;
}

/*-----------------------------------------------------------------------
    Purpose:    Print a complex rectangular matrix with a given format
                and labels.

    Usage:      CALL WRCRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL)

    Arguments:
       TITLE  - Character string specifying the title.  (Input)
                TITLE = ' ' suppresses printing of the title.
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - Complex NRA by NCA matrix to be printed.  (Input)
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
                Because a complex number consists of two parts (a real
                and an imaginary lv_part), two edit descriptors are used for
                printing a single complex number.  FMT = '(E10.3, F10.3)'
                specifies an E format for the real lv_part and an F format
                for the imaginary lv_part.  FMT = '(F10.3)' uses an F format
                for both the real and imaginary parts.  If the end of FMT
                is encountered and if all columns of the matrix have not
                been printed, format control continues with the first
                format in FMT.  Even though the matrix A is complex, an I
                format can be used to print the integer parts of the real
                and imaginary components of each complex number.  The
                most useful format is a special format, called
                "W format," that can be used to specify pretty formats
                automatically.  Set FMT = '(W10.4)' if you want a single
                D, E, F, or I format selected automatically with field
                width 10 and with 4 significant digits.  See Remark 4 for
                a general description of the W format.  FMT may contain
                only D, E, F, G, I, or W edit descriptors, e.g., the X
                descriptor is not allowed.  FMT must contain exactly one
                set of parentheses.
       RLABEL - CHARACTER*(*) vector of labels for rows of A.  (Input)
                If rows are to be numbered consecutively 1, 2, ..., NRA,
                use RLABEL(1) = 'NUMBER'.  If no row labels are desired,
                use RLABEL(1) = 'NONE'.  Otherwise RLABEL is a vector of
                length NRA containing the labels.
       CLABEL - CHARACTER*(*) vector of labels for columns of A.  (Input)
                If columns are to be numbered consecutively 1, 2, ...,
                NCA, use CLABEL(1) = 'NUMBER'.  If no column labels are
                desired, use CLABEL(1) = 'NONE'.  Otherwise, CLABEL is a
                vector of length NCA + 1 containing the column headings.
                CLABEL(1) is the heading for the row labels, and for
                J = 1, 2, ..., NCA, CLABEL(1+J) is the column heading
                for the J-th column.

    Remarks:
    1. Automatic workspace is used only if W formats are used and IMSL
       routine imsl_write_options has previously been invoked with IOPT = -2 and
       ISET = 0.  In this case, workspace usage is
                WRCRL    20*NCA character units, or
                DWRCRL   20*NCA character units.
       Workspace may be explicitly provided, if desired, by use of
       W2CRL/DW2CRL.  The reference is
                CALL W2CRL (TITLE, NRA, NCA, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, CHWK)
       The additional argument is
       CHWK   - CHARACTER*10 work vector of length 2*NCA.  This workspace
                is referenced only if W formats are used and IMSL routine
                imsl_write_options has been previously invoked with IOPT = -1 and
                ISET = 0.  Otherwise, CKWK is not referenced and can be a
                CHARACTER*10 vector of length one.

    2. The output appears in the following form:
                                 TITLE
       CLABEL(1)      CLABEL(2)     CLABEL(3)     CLABEL(4)
       RLABEL(1)  (xxxxx,xxxxx) (xxxxx,xxxxx) (xxxxx,xxxxx)
       RLABEL(2)  (xxxxx,xxxxx) (xxxxx,xxxxx) (xxxxx,xxxxx)

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

void imsl_c_wrcrl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel, irlab,
                iclab)
	Mchar           *title;
	Mint             nra, nca, irlab, iclab;
	Mf_complex        a[];
Mint             lda, itring;
Mchar           *fmt;
Mchar           **rlabel;
Mchar           **clabel;
{
        Mchar            **chwk = NULL;
	Mint             wfmt;
	Mint             i, maxfw, n;


	/* Check for terminal errors and return WFMT = .TRUE. if there is a
	 * single W format and the maximum field width, MAWFW.*/
        imsl_write_format(nra, nca, lda, itring, fmt, &maxfw,
                            "WeEfgGdiouxX", "", &n, &wfmt);
	if (imsl_n1rty(0) != 0)  goto L_9000;
	if (wfmt>0 && n>1) {
                chwk = (Mchar **)malloc(2*nca*sizeof(Mchar*));
		chwk[0] = (Mchar *)malloc(2*nca*11*sizeof(Mchar));
                for (i=1;i<2*nca;i++) chwk[i] = *chwk + i*11;
        	l_w2crl(title, nra, nca, a, lda, itring, fmt, rlabel,
                                clabel, irlab, iclab, chwk);
                free(*chwk);
                free(chwk);
                goto L_9000;
	}
	l_w2crl(title, nra, nca, a, lda, itring, fmt, rlabel, clabel,
                irlab, iclab, chwk);

L_9000: return;
}

/*
    Purpose:    Write a complex rectangular matrix with a specified
                format and labels.

    Usage:      CALL W3CRL (TITLE, A, LDA, ITRING, FMT, RLABEL,
                            CLABEL, WFMT, IROW1, IROW2, ICOL1, ICOL2,
                            LROWLB, LDATA, NSTAG, LTITLE, IRLAB, ICLAB,
                            MAXRLW, IHOTZR, IPAGEW, ICENTR,
                            FMTA, NOCHWK, CHWK)

    Arguments:
       TITLE  - See WRCRL.
       A      - See WRCRL.
       LDA    - See WRCRL.
       ITRING - See WRCRL.
       FMT    - See WRCRL.
       RLABEL - See WRCRL.
       CLABEL - See WRCRL.
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
       NDENTR - Indentation for continuation lines for row labels.
                (Input)
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
       CHWK   - CHARACTER*10 vector of length 2*(ICOL2-ICOL1+1)
                containing FORTRAN formats.  (Input if NOCHWK = 0)
                If NOCHWK = 1, CHWK is not referenced.
*/

static  void l_w3crl_f(title, a, lda, itring, fmt, rlabel, clabel, wfmt, irow1,
	irow2, icol1, icol2, lrowlb, ldata, nstag, ltitle, irlab, iclab,
	maxrlw, ihotzr, ipagew, icentr, fmta, chwk)
	Mchar           *title;
	Mf_complex       *a;
	Mint             lda, itring;
	Mchar           *fmt;
	Mchar           **rlabel;
	Mchar           **clabel;
	Mint             wfmt;
	Mint             irow1, irow2, icol1, icol2, lrowlb, ldata, nstag,
                        ltitle, irlab, iclab, maxrlw, ihotzr, ipagew, icentr;
	Mchar           *fmta;
	Mchar           **chwk;
{
	Mchar            blanks[256], *fmti[2], line[256], *N, *W;
	Mint             i, done, icol, ipad, irow, istag, lfield, lline, lout, middle,
                         nline, nlinei, nwidth[2];

        for (i=0;i<256;i++) blanks[i]=' '; blanks[255]='\0';
        strcpy(line,blanks);
                                	/* Get width of output, LOUT. */
	lout = max(ltitle, ldata);
	middle  = ((icentr == 1)?ipagew/2:lout/2);
                                        /* Print title */
        nline = 0;
        imsl_write_title(title, ipagew, &nline, line, middle);
        strcpy(line,blanks);

        imsl_write_labels (iclab, clabel, line, blanks, &nline, lrowlb, maxrlw, ihotzr,
                      middle, icol1, icol2, ipagew, nstag, fmt, ldata, 1);
        strcpy(line,blanks);
							/* Print A. */
	for (irow = irow1; irow <= irow2; irow++) {
		nline = 1;
		imsl_w8rrl(irow, irlab, rlabel, lrowlb, maxrlw, ihotzr, NDENTR,
                        nline, &lline, line, &done);
		if (lline > 0)  lline += NSPACE;
		istag = 0;
		for (icol = icol1; icol <= icol2; icol++) {
			                                /* Get format */
			imsl_w5rrl_f(icol, iclab, clabel, fmt, nwidth, &lfield,
                               fmti, &nlinei, 2);
			                            /* Get padding before the number. */
			ipad = lfield - (nwidth[0] + nwidth[1] + 3);
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
                        if (wfmt) fmti[0] = fmta;
                        N = strchr(fmti[0]+1,(Mint)'%');
                        W = strchr(fmti[0],(Mint)'W');
                        if (!N && W) N = W+1;
                        if (W && W < N) fmti[0] = chwk[(icol-icol1)*2];
                        if (wfmt) fmti[1] = fmta;
                        N = strchr(fmti[1]+1,(Mint)'%');
                        W = strchr(fmti[1],(Mint)'W');
                        if (!N && W) N = W+1;
                        if (W && W < N) fmti[1] = chwk[(icol-icol1)*2+1];
			            /* Append entry A(IROW,ICOL) to LINE. */
	                lline += ipad;
			strncpy(&line[lline], l_c1tri(itring, irow, icol,
			    A(icol-1, irow-1), nwidth, fmti),
                                nwidth[0]+nwidth[1]+3);
			lline += nwidth[0] + nwidth[1] + 3 + NSPACE;
		}
		                            /* Print last line for IROW-th row of A.*/
		lline = ldata;
		imsl_c1nter(middle, &lline, line);
		imsl_write_line(lline, line);
                strcpy(line,blanks);
L_110:          if (!done) {
			/* Number of lines for IROW-th row label exceeds that
			 * for the IROW-th row of data. Print the remainder
			 * of the row label. */
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
/*  -----------------------------------------------------------------------
    Purpose:    Convert a matrix element to a character string according
                to a supplied format.  Dots are returned if not-a-number
                is input.  If the matrix element is not in the triangle
                to be printed, blanks are returned.

    Usage:      C1TRI(ITRING, I, J, X, NWID, FMT)

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
       X      - The COMPLEX element to be written.  (Input)
       NWID   - Vector of length 2 containing the width of each
                format specification in FMT. (Input)
                NWID(1) and NWID(2) must be greater than 1 and less than
                or equal to 40.
       FMT    - CHARACTER*(*) vector of length two two be used in writing
                X.  (Input)
                FMT(1) is used to write the real lv_part and FMT(2) is used
                to write the imaginaray lv_part.
       C1TRI  - A CHARACTER*83 string which contains blanks if the
                MATRIX element is not to be printed.  Otherwise, if X is
                not-a-number C1TRI contains dots, and if X is not
                not-a-number C1TRI contains a representation of X in the
                form specified by FMT.
                (Output)
  -----------------------------------------------------------------------
 */

static Mchar *l_c1tri(itring, i, j, x, nwid, fmt)
	Mint             itring, i, j;
	Mf_complex        *x;
	Mint             nwid[2];
        Mchar           *fmt[2];
{
	Mchar            wmiss1[41], wmiss2[41];
	static Mchar     c1tri_v[84];
	Mint             iwrite;
	Mfloat           xi, xr;

	iwrite = 0;
	if (itring == 0) {iwrite = 1;
	} else if (itring == 1) {if (i <= j) iwrite = 1;
	} else if (itring == 2) {if (i < j) iwrite = 1;
	} else if (itring == -1) {if (i >= j) iwrite = 1;
	} else if (itring == -2) if (i > j) iwrite = 1;

	if (iwrite == 0) strcpy(c1tri_v, "                                                                                   ");
	else {
		xi = x->im;
		xr = x->re;
		strcpy(wmiss1, imsl_w1iss(&xr, fmt[0], nwid[0]));
		strcpy(wmiss2, imsl_w1iss(&xi, fmt[1], nwid[1]));
                c1tri_v[0]='(';
		strncpy(c1tri_v+1, wmiss1, nwid[0]);
                strncpy(c1tri_v+1+nwid[0], ",", 1);
                strncpy(c1tri_v+2+nwid[0], wmiss2, nwid[1]);
                strncpy(c1tri_v+2+nwid[0]+nwid[1], ")", 1);
	}
	return (c1tri_v);
}

/*  -----------------------------------------------------------------------
    Purpose:    Store formats and convert W formats to FORTRAN formats.

    Usage:      CALL W4CRL (ITRING, IROW1, IROW2, ICOL1, ICOL2, A, LDA,
                            FMT, CHWK)

    Arguments:
       ITRING - See WRCRL.  (Input)
       IROW1  - First row of matrix to be printed.  (Input)
       IROW2  - Last row of matrix to be printed.  (Input)
       ICOL1  - First column of matrix to be printed.  (Input)
       ICOL2  - Last column of matrix to be printed.  (Input)
       A      - Matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in
                the dimension statement of the calling program.  (Input)
       FMT    - Character string containing format.  (Input)
                See WRCRL.
       CHWK   - CHARACTER*10 vector of length 2*(ICOL2-ICOL1 + 1)
                containing some output formats for real and imaginary
                parts of elements of columns ICOL1, ... , ICOL2 of A.
                (Output)
                A FORTRAN format is output in CHWK(I) only if the
                (ICOL1+I-1)-th format in FMT is a W format.
  -----------------------------------------------------------------------
 */


static void l_w4crl_f(itring, irow1, irow2, icol1, icol2, a, lda, fmt, chwk)
	Mint             itring, irow1, irow2, icol1, icol2;
	Mf_complex       *a;
	Mint             lda;
	Mchar           *fmt;
	Mchar           **chwk;
{
	Mchar            *fmti, *N, *W, aa='.', bb='W';
	Mint             i, icol, ir1, nrow, nsig, nwidth;


	i = 0;
	for (icol = icol1; icol <= icol2; icol++) {
		                                        /* REAL lv_part. */
		i += 1;
		fmti = imsl_w7rrl(2*icol-1, fmt);
		imsl_w6rrl(fmti, 1, &aa, &bb, &nwidth);
		if (nwidth == 0)  nwidth = 10;
		imsl_w12rl(itring, irow1, irow2, icol, &ir1, &nrow);
                N = strchr(fmti+1,(Mint)'%');
                W = strchr(fmti,(Mint)'W');
                if (!N && W) N = W+1;
                if (W && W < N) {
			            /* I-th format is W. Compute actual FORTRAN format to
			             * be used.*/

			imsl_w6rrl(fmti, 0, &aa, &bb, &nsig);
                        if (nsig == 0)  nsig = 4;
			if (nrow <= 0) {
				i += 1;
				goto L_10;
			}
			strcpy(chwk[i-1], l_f1tcr_f(nrow, 1, A(icol-1,
				ir1-1), nrow, 0, 0, nwidth, nsig));
		}
                                		/* Imaginary lv_part. */
		i += 1;
		fmti = imsl_w7rrl(2*icol, fmt);
		imsl_w6rrl(fmti, 1, &aa, &bb, &nwidth);
                if (nwidth == 0)  nwidth = 10;
                N = strchr(fmti+1,(Mint)'%');
                W = strchr(fmti,(Mint)'W');
                if (!N && W) N = W+1;
                if (W && W < N) {
			            /* I-th format is W. Compute actual FORTRAN format to
			             * be used.*/
			imsl_w6rrl(fmti, 0, &aa, &bb, &nsig);
                        if (nsig == 0)  nsig = 4;
			strcpy(chwk[i-1], l_f1tcr_f(nrow,1,A(icol-1,ir1-1),
                               nrow, 0, 1, nwidth, nsig));
		}
L_10: ;
	}
	return;
}

/*  -----------------------------------------------------------------------
    Purpose:    Return a pretty FORTRAN format for a complex rectangular
                matrix.

    Usage:      F1TCR(NRA, NCA, A, LDA, ITRING, IPART, NWIDTH, NSIG)

    Arguments:
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - NRA by NCA complex matrix.  (Input)
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
       IPART  - Real/Imaginary lv_part option.
                IPART   Meaning
                   0    Only the real lv_part of A is examined.
                   1    Only the imaginary lv_part of A is examined.
                   2    Both the real and imaginary lv_part of A is
                        examined.
       NWIDTH - Width of FORTRAN format field.  (Input)
                NWIDTH = 10 or greater is suggested for pretty formats.
                NWIDTH must equal 3, 4, ..., or 40.
       NSIG   - Number of significant digits generally printed.  (Input)
                NSIG = NWIDTH-6 is suggested for pretty formats.  NSIG
                must not exceed NWIDTH-2.
       F1TCR  - A CHARACTER string of length 10 that contains the FORTRAN
                format.  (Output)
                F1TCR = '(1PE  .  )' or '(0PF  .  )' or '(I  )' where
                appropriate integers appear right justified in the
                blanks.

    Remarks:
    1. The format returned in FMTCR is determined by the first rule
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
             ISIG = INT(0.05-LOG10(AMACH(4)) for F1TCR and
             ISIG = INT(0.05-LOG10(DMACH(4)) for DF1TCR.  Here the use of
             ISIG is intended to prevent printing of garbage digits.
         e.  If MIN(NSIG,ISIG) significant digits can be printed using
             an F format with field width NWIDTH, an F format is
             returned.
         f.  An E format is returned containing MIN(NSIG,ISIG).

    2. The F format returned suppresses the printing of trailing
       zeroes after the decimal point.
  -----------------------------------------------------------------------
 */

static Mchar *l_f1tcr_f(nra, nca, a, lda, itring, ipart, nwidth, nsig)
	Mint             nra, nca;
	Mf_complex       *a;
	Mint             lda, itring, ipart, nwidth, nsig;
{
	Mchar            temp[11], temp1[11], width[3];
	static Mchar     f1tcr_v[11];
	Mint             i, iflag, indx, ipart1, ir1, j, nrow;
	Mfloat           aij, amax;

	sprintf(width, "%d", nwidth);
	for (j = 1; j <= nca; j++) {
		imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		for (i = ir1; i <= (ir1 + nrow - 1); i++) {
			if (ipart == 0) {
				aij = (*A(j - 1, i - 1)).re;
				if (!imsl_ifnan(aij))  goto L_30;
			} else if (ipart == 1) {
				aij = (*A(j - 1, i - 1)).im;
				if (!imsl_ifnan(aij))  goto L_30;
			} else if (ipart == 2) {
				aij = (*A(j - 1, i - 1)).re;
				if (!imsl_ifnan(aij))  goto L_30;
				aij = (*A(j - 1, i - 1)).im;
				if (!imsl_ifnan(aij))  goto L_30;
			}
		}
	}
	                /*Either all elements of A are not a number or NRA or NCA is less
			 * than one. Return an arbitrary format.*/
	aij = imsl_amach(6);
	strcpy(f1tcr_v, imsl_fmtx(&aij, nwidth, nsig));
	goto L_9000;
	                                    /*Determine if I format is appropriate.*/
L_30:
	amax = F_ZERO;
	if (ipart == 0) {
		ipart1 = 0;
	} else if (ipart == 1) {
		ipart1 = 1;
	} else {
		ipart1 = 0;
	}
	for (j = 1; j <= nca; j++) {
		imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		if (nrow <= 0)  goto L_40;
		l_c1fmt(nrow, A(j - 1, ir1 - 1), nwidth, ipart1, &iflag);
		if (iflag == 0)  goto L_50;
		if (ipart == 2) {
			l_c1fmt(nrow, A(j-1,ir1-1), nwidth, 1, &iflag);
			if (iflag == 0)  goto L_50;
		}
L_40:
		;
	}
L_50:
	if (iflag == 1) {
		                                            /* Return an I format. */
		strcpy(f1tcr_v, "%");
                strcat(f1tcr_v, width);
                strcat(f1tcr_v, "d");
		goto L_9000;
	} else {
		                                     /* Find largest absolute value. */
		amax = F_ZERO;
		for (j = 1; j <= nca; j++) {
			imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
			if (nrow <= 0)  goto L_60;
			l_i1cmx(nrow, A(j - 1, ir1 - 1), 1, ipart1, &indx);
			if (indx > 0) {
				if (ipart1 == 0) {
					amax = max(fabs(((*A(j-1,ir1+indx-2)).re)),
						    amax);
				} else {
					amax = max(fabs(((*A(j-1,ir1+indx-2)).im)),
						    amax);
				}
			}
			if (ipart == 2) {
				l_i1cmx(nrow, A(j-1, ir1-1), 1, 1, &indx);
				if (indx > 0) {
					amax = max(fabs(((*A(j-1,ir1+indx-2)).im)),
						    amax);
				}
			}
	L_60:
			;
		}
		strcpy(temp, imsl_fmtx(&amax, nwidth, nsig));
		if (temp[strlen(temp)-1]=='f') {
			/* If F format, see if there are trailing zeroes that
			 * can be stripped by modifing the F format.*/
			strcpy(temp1, temp);
			l_f3tcr_f(nra, nca, a, lda, itring, ipart,
				temp1, temp);
		}
		strcpy(f1tcr_v, temp);
	}
L_9000:
	return (f1tcr_v);
}				/* end of function */
/*  -----------------------------------------------------------------------
    Purpose:    Change an F format so that trailing zeroes are not
                printed after the decimal point for a complex matrix.

    Usage:      CALL F3TCR (NRA, NCA, A, LDA, ITRING, IPART, FFMTI,
                            FFMTO)

    Arguments:
       NRA    - Number of rows.  (Input)
       NCA    - Number of columns.  (Input)
       A      - NRA by NCA complex matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ITRING - Triangle option.  (Input)
                ITRING  Action
                  0     Full matrix is used.
                  1     Upper triangle of A is used.
                  2     Upper triangle of A excluding the diagonal of A
                        is used.
                 -1     Lower triangle of A is used.
                 -2     Lower triangle of A excluding the diagonal of A
                        is used.
       IPART  - Real/Imaginary lv_part option.  (Input)
                IPART   Action
                  0     Only the real lv_part of A is examined.
                  1     Only the imaginary lv_part of A is examined.
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

static void l_f3tcr_f(nra, nca, a, lda, itring, ipart, ffmti, ffmto)
	Mint             nra, nca;
	Mf_complex       *a;
	Mint             lda, itring, ipart;
	Mchar           *ffmti;
	Mchar           *ffmto;
{
	Mint             i, intx, ipart1, ir1, j, ndec, ndec1, ndec2, nrow;
	Mfloat           aij, rntx, pwr;
        Mchar            *itmp, temp[4];

	strcpy(ffmto, ffmti);
        itmp = strchr(ffmto,(Mint)'.');
        if (isdigit((Mint)itmp[1])) ndec = (Mint)itmp[1] - (Mint)'0';
        if (isdigit((Mint)itmp[2])) ndec = 10*ndec + (Mint)itmp[2] - (Mint)'0';
        if (ndec != 0) {
            pwr = pow((double)F_TEN, (double)ndec);
	    ndec1 = 0;
	    if (ipart == 0) {
		ipart1 = 0;
	    } else if (ipart == 1) {
		ipart1 = 1;
	    } else {
		ipart1 = 0;
	    }
	    for (j = 1; j <= nca; j++) {
		    imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
		    for (i = ir1; i <= (ir1 + nrow - 1); i++) {
			    ndec2 = ndec;
			    if (ipart1 == 0) {
				    aij = (*A(j - 1, i - 1)).re;
			    } else {
				    aij = (*A(j - 1, i - 1)).im;
			    }
			    if (!imsl_ifnan(aij)) {
			        rntx = aij*pwr + ((aij>=F_ZERO)?F_HALF:-F_HALF);
                                intx = fmod(rntx,F_TEN);
			        if (intx != 0) return;
	                        while (intx == 0) {
			            ndec2 -= 1;
			            if (ndec2 == 0) break;
			            rntx /= F_TEN;
                                    intx = fmod(rntx,F_TEN);
			        }
			        ndec1 = max(ndec1, ndec2);
	                    }
		    }
	    }
	    if (ipart == 2) {
		    for (j = 1; j <= nca; j++) {
			    imsl_w12rl(itring, 1, nra, j, &ir1, &nrow);
			    for (i = ir1; i <= (ir1 + nrow - 1); i++) {
				    ndec2 = ndec;
				    aij = (*A(j - 1, i - 1)).im;
				    if (imsl_ifnan(aij)) {
                                        rntx = aij*pwr + ((aij>=F_ZERO)?F_HALF:-F_HALF);
                                        intx = fmod(rntx,F_TEN);
				        if (intx != 0) return;
	                                while (intx == 0) {
				            ndec2 -= 1;
				            if (ndec2 == 0)break;
				            rntx /= F_TEN;
                                            intx = fmod(rntx,F_TEN);
				        }
				        ndec1 = max(ndec1, ndec2);
                                    }
			    }
		    }
	    }
            sprintf(temp,"%df",ndec1);
            strcpy(itmp+1,temp);
        }
	return;
}
/*  -----------------------------------------------------------------------
    Purpose:    Find the smallest index of the component of a complex
                vector x (possibly containing NaN) having maximum
                absolute value.

    Usage:      CALL I1CMX (N, X, INCX, IOPT, INDX)

    Arguments:
       N      - Length of vector x.  (Input)
                The vector x refers to a specific subvector obtained from
                the vector X.  See the description of the argument INCX.
       X      - COMPLEX vector of length N*INCX.  (Input)
       INCX   - Increment between the referenced elements of X.  (Input)
                X(1+((I-1)*INCX)) (which is defined to be x(I)) is
                referenced, for I = 1, 2, ... , N.  INCX must be
                positive.
       IOPT   - Real/imaginary option.  (Input)
                IOPT  Action
                  0   Only the real lv_part of x is examined.
                  1   Only the imaginary lv_part of x is examined.
       INDX   - The smallest index of the component of x having maximum
                absolute value.  (Output)
                Components of x equal to NaN are excluded from
                computations.  If N = 0 or all elements of x equal NaN,
                INDX = 0.
  -----------------------------------------------------------------------
 */

static void l_i1cmx(n, x, incx, iopt, indx)
	Mint             n;
	Mf_complex        x[];
Mint             incx, iopt, *indx;
{
	Mint             i, j;
	Mfloat           absxi, xi, xmax;


	*indx = 0;
	j = 1;
	for (i = 1; i <= n; i++) {
		if (iopt == 0) {
			xi = (x[j - 1]).re;
		} else {
			xi = (x[j - 1]).im;
		}
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
    Purpose:    Indicate if an I format for printing is appropriate.

    Usage:      CALL C1FMT (N, X, NWIDTH, IOPT, IFLAG)

    Arguments:
       N      - Length of X.  (Input)
       X      - COMPLEX vector of length N.  (Input)
       NWIDTH - Field width.  (Input)
       IOPT   - Real/imaginary option.  (Input)
                IOPT  Action
                  0   Only the real lv_part of X is examined.
                  1   Only the imaginary lv_part of X is examined.
       IFLAG  - Indicator for I format.  (Output)
                IFLAG = 1 means an I format is appropriate because (1)
                all elements of X (not equal to NaN) have no fractional
                parts and (2) the field width NWIDTH is large enough for
                printing the elements of X.  Otherwise, IFLAG = 0.
  -----------------------------------------------------------------------
 */

static void l_c1fmt(n, x, nwidth, iopt, iflag)
	Mint             n;
	Mf_complex        x[];
        Mint             nwidth, iopt, *iflag;
{
	Mint             i;
	Mfloat           xi, xint, xint1, xint2;

	*iflag = 0;
	xint1 = INT_MIN;
	xint2 = INT_MAX;
	if (nwidth < log10(imsl_amach(2))) {
		xint1 = -pow((double)F_TEN, (double)(nwidth-1));
		xint2 = -F_TEN * xint1;
	}
	for (i = 1; i <= n; i++) {
		if (iopt == 0) {
			xi = x[i - 1].re;
		} else {
			xi = x[i - 1].im;
		}
		if (!imsl_ifnan(xi))  {
		    if (fabs(xi) > INT_MAX)  goto L_9000;
		    xint = (Mint) (xi);
		    if (xint != xi)  goto L_9000;
		    if (xint <= xint1)  goto L_9000;
		    if (xint >= xint2)  goto L_9000;
                }
	}
	*iflag = 1;
L_9000:
	return;
}
