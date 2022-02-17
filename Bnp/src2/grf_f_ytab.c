
/*  A Bison parser  , made from GrfnLang.y
 by  GNU Bison version 1.25
  */

#define YYBISON 1 /* Identify Bison output.  */

#define alloca

#define QNAME 258
#define NAME 259
#define NUMBER 260
#define ASSIGN 261
#define LE 262
#define GE 263
#define AND 264
#define OR 265
#define IF 266
#define EQ 267
#define DO 268
#define IFFIX 269
#define UMINUS 270

#line 1 "GrfnLang.y"

#include "grf_h_all.h"

#line 4 "GrfnLang.y"
typedef union {
  double dval;
  char sval[300];
  COMLL_PTR cval;
} YYSTYPE;
#include "stdio.h"

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif

#define YYFINAL 63
#define YYFLAG 32768
#define YYNTBASE 28

#define YYTRANSLATE(x) ((unsigned)(x) <= 270 ? yytranslate[x] : 31)

static const char yytranslate[] = {
    0,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    23, 24, 19, 17, 27, 18, 2,  20, 2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    15, 2,  16, 2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 25, 2, 26, 21, 2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 2, 2, 2, 2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2,  2, 1, 2, 3, 4,
    5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 22};

#if YYDEBUG != 0
static const short yyprhs[] = {0,  0,  2,  6,  10, 14, 18, 22, 26,
                               30, 34, 38, 42, 46, 50, 54, 58, 61,
                               65, 70, 75, 84, 91, 93, 95, 97, 99};

static const short yyrhs[] = {
    29, 0,  29, 13, 29, 0,  29, 12, 29, 0,  29, 9,  29, 0,  29, 10, 29, 0,
    29, 7,  29, 0,  29, 8,  29, 0,  29, 15, 29, 0,  29, 16, 29, 0,  29, 17,
    29, 0,  29, 18, 29, 0,  29, 19, 29, 0,  29, 20, 29, 0,  29, 21, 29, 0,
    29, 6,  29, 0,  18, 29, 0,  23, 29, 24, 0,  4,  23, 30, 24, 0,  4,  25,
    30, 26, 0,  11, 23, 29, 27, 29, 27, 29, 24, 0,  14, 23, 29, 27, 29, 24,
    0,  4,  0,  3,  0,  5,  0,  29, 0,  30, 27, 29, 0};

#endif

#if YYDEBUG != 0
static const short yyrline[] = {0,   27,  34,  41,  55,  69,  83,  97,  111,
                                125, 139, 153, 167, 181, 195, 209, 216, 230,
                                231, 238, 245, 262, 275, 284, 291, 301, 305};
#endif

#if YYDEBUG != 0 || defined(YYERROR_VERBOSE)

static const char *const yytname[] = {
    "$",      "error", "$undefined.", "QNAME", "NAME",      "NUMBER",
    "ASSIGN", "LE",    "GE",          "AND",   "OR",        "IF",
    "EQ",     "DO",    "IFFIX",       "'<'",   "'>'",       "'+'",
    "'-'",    "'*'",   "'/'",         "'^'",   "UMINUS",    "'('",
    "')'",    "'['",   "']'",         "'  ,'", "statement", "expression",
    "args",   NULL};
#endif

static const short yyr1[] = {0,  28, 29, 29, 29, 29, 29, 29, 29,
                             29, 29, 29, 29, 29, 29, 29, 29, 29,
                             29, 29, 29, 29, 29, 29, 29, 30, 30};

static const short yyr2[] = {0, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                             3, 3, 2, 3, 4, 4, 8, 6, 1, 1, 1, 1, 3};

static const short yydefact[] = {
    0,  23, 22, 24, 0, 0, 0,  0, 1, 0, 0,  0, 0,  16, 0,  0,
    0,  0,  0,  0,  0, 0, 0,  0, 0, 0, 0,  0, 0,  25, 0,  0,
    0,  0,  17, 15, 6, 7, 4,  5, 3, 2, 8,  9, 10, 11, 12, 13,
    14, 18, 0,  19, 0, 0, 26, 0, 0, 0, 21, 0, 20, 0,  0,  0};

static const short yydefgoto[] = {61, 29, 30};

static const short yypact[] = {
    87, -32768, 31,     -32768, -9,  6,      87,  87,     162,   87,  87,
    87, 87,     -32768, 105,    87,  87,     87,  87,     87,    87,  87,
    87, 87,     87,     87,     87,  87,     87,  162,    -14,   -18, 24,
    52, -32768, 178,    -16,    -16, 208,    193, 223,    178,   -16, -16,
    27, 27,     14,     14,     14,  -32768, 87,  -32768, 87,    87,  162,
    68, 124,    87,     -32768, 143, -32768, 38,  49,     -32768};

static const short yypgoto[] = {-32768, 0, 45};

#define YYLAST 244

static const short yytable[] = {
    8,  24, 25, 26, 27, 28, 13, 14, 51, 50, 49, 32, 33, 50, 11, 35, 36, 37, 38,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 12, 15, 16, 17, 18, 19, 28, 20, 21,
    62, 22, 23, 24, 25, 26, 27, 28, 26, 27, 28, 63, 54, 52, 55, 56, 9,  31, 10,
    59, 15, 16, 17, 18, 19, 0,  20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 15, 16,
    17, 18, 19, 53, 20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 1,  2,  3,  0,  0,
    57, 0,  0,  4,  0,  0,  5,  0,  0,  0,  6,  0,  0,  0,  0,  7,  15, 16, 17,
    18, 19, 0,  20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 0,  0,  34, 15, 16, 17,
    18, 19, 0,  20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 0,  0,  58, 15, 16, 17,
    18, 19, 0,  20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 0,  0,  60, 15, 16, 17,
    18, 19, 0,  20, 21, 0,  22, 23, 24, 25, 26, 27, 28, 15, 16, 17, 18, 19, 0,
    20, 0,  0,  22, 23, 24, 25, 26, 27, 28, 16, 17, 18, 0,  0,  20, 0,  0,  22,
    23, 24, 25, 26, 27, 28, 16, 17, 0,  0,  0,  20, 0,  0,  22, 23, 24, 25, 26,
    27, 28, 16, 17, 0,  0,  0,  0,  0,  0,  22, 23, 24, 25, 26, 27, 28};

static const short yycheck[] = {
    0,  17, 18, 19, 20, 21, 6,  7,  26, 27, 24, 11, 12, 27, 23, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 23, 6,  7,  8,  9,  10, 21, 12, 13,
    0,  15, 16, 17, 18, 19, 20, 21, 19, 20, 21, 0,  50, 27, 52, 53, 23, 10, 25,
    57, 6,  7,  8,  9,  10, -1, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, 6,  7,
    8,  9,  10, 27, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, 3,  4,  5,  -1, -1,
    27, -1, -1, 11, -1, -1, 14, -1, -1, -1, 18, -1, -1, -1, -1, 23, 6,  7,  8,
    9,  10, -1, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, -1, -1, 24, 6,  7,  8,
    9,  10, -1, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, -1, -1, 24, 6,  7,  8,
    9,  10, -1, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, -1, -1, 24, 6,  7,  8,
    9,  10, -1, 12, 13, -1, 15, 16, 17, 18, 19, 20, 21, 6,  7,  8,  9,  10, -1,
    12, -1, -1, 15, 16, 17, 18, 19, 20, 21, 7,  8,  9,  -1, -1, 12, -1, -1, 15,
    16, 17, 18, 19, 20, 21, 7,  8,  -1, -1, -1, 12, -1, -1, 15, 16, 17, 18, 19,
    20, 21, 7,  8,  -1, -1, -1, -1, -1, -1, 15, 16, 17, 18, 19, 20, 21};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */

/* Skeleton output parser for bison  ,
   Copyright (C) 1984  , 1989  , 1990 Free Software Foundation  , Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2  , or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful  ,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not  , write to the Free Software
   Foundation  , Inc.  , 675 Mass Ave  , Cambridge  , MA 02139  , USA.  */

/* As a special exception  , when this file is copied by Bison into a
   Bison output file  , you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

#ifndef alloca
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined(__STDC__) && defined(sparc)) || defined(__sparc__) ||            \
    defined(__sparc) || defined(__sgi)
#include "alloca.h"
#else /* not sparc */
#if defined(MSDOS) && !defined(__TURBOC__)
#include "malloc.h"
#else /* not MSDOS  , or __TURBOC__ */
#if defined(_AIX)
#include "malloc.h"
#pragma alloca
#else /* not MSDOS  , __TURBOC__  , or _AIX */
#ifdef __hpux
#ifdef __cplusplus
extern "C" {
void *alloca(unsigned int);
};
#else  /* not __cplusplus */
void *alloca();
#endif /* not __cplusplus */
#endif /* __hpux */
#endif /* not _AIX */
#endif /* not MSDOS  , or __TURBOC__ */
#endif /* not sparc.  */
#endif /* not GNU C.  */
#endif /* alloca not defined.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions  , each action
   as one case of the switch.  */

#define yyerrok (yyerrstatus = 0)
#define yyclearin (yychar = YYEMPTY)
#define YYEMPTY -2
#define YYEOF 0
#define YYACCEPT return (0)
#define YYABORT return (1)
#define YYERROR goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR  , for GCC.
   Once GCC version 2 has supplanted version 1  , this can go.  */
#define YYFAIL goto yyerrlab
#define YYRECOVERING() (!!yyerrstatus)
#define YYBACKUP(token, value)                                                 \
  do                                                                           \
    if (yychar == YYEMPTY && yylen == 1) {                                     \
      yychar = (token), yylval = (value);                                      \
      yychar1 = YYTRANSLATE(yychar);                                           \
      YYPOPSTACK;                                                              \
      goto yybackup;                                                           \
    } else {                                                                   \
      yyerror("syntax error: cannot back up");                                 \
      YYERROR;                                                                 \
    }                                                                          \
  while (0)

#define YYTERROR 1
#define YYERRCODE 256

#ifndef YYPURE
#define YYLEX yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX yylex(&yylval, &yylloc, YYLEX_PARAM)
#else
#define YYLEX yylex(&yylval, &yylloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX yylex(&yylval, YYLEX_PARAM)
#else
#define YYLEX yylex(&yylval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant  , generate the variables here */

#ifndef YYPURE

int yychar;     /*  the lookahead symbol		*/
YYSTYPE yylval; /*  the semantic value of the		*/
                /*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc; /*  location data for the lookahead	*/
                /*  symbol				*/
#endif

int yynerrs; /*  number of parse errors so far       */
#endif       /* not YYPURE */

#if YYDEBUG != 0
int yydebug; /*  nonzero means print parse trace	*/
/* Since this is uninitialized  , it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
int yyparse(void);
#endif

#if __GNUC__ > 1 /* GNU C and GNU C++ define this.  */
#define __yy_memcpy(TO, FROM, COUNT) __builtin_memcpy(TO, FROM, COUNT)
#else /* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void __yy_memcpy(to, from, count) char *to;
char *from;
int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void __yy_memcpy(char *to, char *from, int count) {
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else  /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

Err yyparse(COMLL_PTR *ret_ptr, GrfnDeal *grfndeal) {
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus; /*  number of tokens to shift before error messages enabled
                    */
  int yychar1 =
      0; /*  lookahead token as an internal (translated) token number */

  short yyssa[YYINITDEPTH];   /*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH]; /*  the semantic value stack		*/

  short *yyss = yyssa;   /*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa; /*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE yylsa[YYINITDEPTH]; /*  the location stack			*/
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;

#define YYPOPSTACK (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval; /*  the variable used to return		*/
  /*  semantic values from the action	*/
  /*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state  , which is found in  yystate  .  */
/* In all cases  , when you get here  , the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1) {
    /* Give user a chance to reallocate the stack */
    /* Use copies of these so that the &'s don't force the real ones into
     * memory. */
    YYSTYPE *yyvs1 = yyvs;
    short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
    YYLTYPE *yyls1 = yyls;
#endif

    /* Get the current used size of the three stacks  , in elements.  */
    int size = yyssp - yyss + 1;

#ifdef yyoverflow
    /* Each stack pointer address is followed by the size of
       the data in use in that stack  , in bytes.  */
#ifdef YYLSP_NEEDED
    /* This used to be a conditional around just the two extra args  ,
       but that might be undefined if yyoverflow is a macro.  */
    yyoverflow("parser stack overflow", &yyss1, size * sizeof(*yyssp), &yyvs1,
               size * sizeof(*yyvsp), &yyls1, size * sizeof(*yylsp),
               &yystacksize);
#else
    yyoverflow("parser stack overflow", &yyss1, size * sizeof(*yyssp), &yyvs1,
               size * sizeof(*yyvsp), &yystacksize);
#endif

    yyss = yyss1;
    yyvs = yyvs1;
#ifdef YYLSP_NEEDED
    yyls = yyls1;
#endif
#else /* no yyoverflow */
    /* Extend the stack our own way.  */
    if (yystacksize >= YYMAXDEPTH) {
      yyerror("parser stack overflow");
      return 2;
    }
    yystacksize *= 2;
    if (yystacksize > YYMAXDEPTH)
      yystacksize = YYMAXDEPTH;
    yyss = (short *)alloca(yystacksize * sizeof(*yyssp));
    __yy_memcpy((char *)yyss, (char *)yyss1, size * sizeof(*yyssp));
    yyvs = (YYSTYPE *)alloca(yystacksize * sizeof(*yyvsp));
    __yy_memcpy((char *)yyvs, (char *)yyvs1, size * sizeof(*yyvsp));
#ifdef YYLSP_NEEDED
    yyls = (YYLTYPE *)alloca(yystacksize * sizeof(*yylsp));
    __yy_memcpy((char *)yyls, (char *)yyls1, size * sizeof(*yylsp));
#endif
#endif /* no yyoverflow */

    yyssp = yyss + size - 1;
    yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
    yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
    if (yydebug)
      fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

    if (yyssp >= yyss + yystacksize - 1)
      YYABORT;
  }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  goto yybackup;
yybackup:

  /* Do appropriate processing given the current state.  */
  /* Read a lookahead token if we need one and don't already have one.  */
  /* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY) {
#if YYDEBUG != 0
    if (yydebug)
      fprintf(stderr, "Reading a token: ");
#endif
    yychar = YYLEX;
  }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0) /* This means end of input. */
  {
    yychar1 = 0;
    yychar = YYEOF; /* Don't call YYLEX any more */

#if YYDEBUG != 0
    if (yydebug)
      fprintf(stderr, "Now at end of input.\n");
#endif
  } else {
    yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
    if (yydebug) {
      fprintf(stderr, "Next token is %d (%s", yychar, yytname[yychar1]);
      /* Give the individual parser a way to print the precise meaning
         of a token  , for further debugging info.  */
#ifdef YYPRINT
      YYPRINT(stderr, yychar, yylval);
#endif
      fprintf(stderr, ")\n");
    }
#endif
  }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce  , -yyn is rule number.
     Positive => shift  , yyn is new state.
       New state is final state => don't bother to shift  ,
       just return success.
     0  , or most negative number => error.  */

  if (yyn < 0) {
    if (yyn == YYFLAG)
      goto yyerrlab;
    yyn = -yyn;
    goto yyreduce;
  } else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

    /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s)  , ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three  , turn off error status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1 - yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug) {
    int i;

    fprintf(stderr, "Reducing via rule %d (line %d)  , ", yyn, yyrline[yyn]);

    /* Print the symbols being reduced  , and their result.  */
    for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
      fprintf(stderr, "%s ", yytname[yyrhs[i]]);
    fprintf(stderr, " -> %s\n", yytname[yyr1[yyn]]);
  }
#endif

  switch (yyn) {

  case 1:
#line 27 "GrfnLang.y"
  {
    yyval.cval = yyvsp[0].cval;
    return NULL;
    ;
    break;
  }
  case 2:
#line 34 "GrfnLang.y"
  {
    COMLL_PTR c;
    c = comll_gotobot(yyvsp[-2].cval);
    c = comll_insert_after(c);
    c->type = COMLL_POP;
    yyval.cval = comll_join(yyvsp[-2].cval, yyvsp[0].cval);
    ;
    break;
  }
  case 3:
#line 41 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(DTOL(yyvsp[-2].cval->dval) == DTOL(yyvsp[0].cval->dval));
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_EQ;
    };
    break;
  }
  case 4:
#line 55 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(DTOL(yyvsp[-2].cval->dval) && DTOL(yyvsp[0].cval->dval));
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_AND;
    };
    break;
  }
  case 5:
#line 69 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(DTOL(yyvsp[-2].cval->dval) || DTOL(yyvsp[0].cval->dval));
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_OR;
    };
    break;
  }
  case 6:
#line 83 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval <= yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_LE;
    };
    break;
  }
  case 7:
#line 97 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval >= yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_GE;
    };
    break;
  }
  case 8:
#line 111 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval < yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_LT;
    };
    break;
  }
  case 9:
#line 125 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval > yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_B_GT;
    };
    break;
  }
  case 10:
#line 139 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval + yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_A_PLUS;
    };
    break;
  }
  case 11:
#line 153 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval - yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_A_MINUS;
    };
    break;
  }
  case 12:
#line 167 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval * yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_A_TIMES;
    };
    break;
  }
  case 13:
#line 181 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval =
          (double)(yyvsp[-2].cval->dval / yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_A_DIVIDE;
    };
    break;
  }
  case 14:
#line 195 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[-2].cval, COMLL_REAL) &&
        comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[-2].cval->dval = pow(yyvsp[-2].cval->dval, yyvsp[0].cval->dval);
      comll_free_list(yyvsp[0].cval);
      yyval.cval = yyvsp[-2].cval;
    } else {
      c = comll_gotobot(yyvsp[-2].cval);
      yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
      c = comll_insert_after(c);
      c->type = COMLL_A_POW;
    };
    break;
  }
  case 15:
#line 209 "GrfnLang.y"
  {
    if (!comll_atom(yyvsp[-2].cval, COMLL_VAR))
      yyerror(GRERR_BAD_ASSIGNMENT);
    yyvsp[-2].cval->type = COMLL_ASSIGN;
    yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);

    ;
    break;
  }
  case 16:
#line 216 "GrfnLang.y"
  {
    COMLL_PTR c;
    if (comll_atom(yyvsp[0].cval, COMLL_REAL)) {
      yyvsp[0].cval->dval *= -1.0;
    } else {
      c = comll_gotobot(yyvsp[0].cval);
      c = comll_insert_after(c);
      c->type = COMLL_REAL;
      c->dval = -1.0;
      c = comll_insert_after(c);
      c->type = COMLL_A_TIMES;
    }
    yyval.cval = yyvsp[0].cval;
    ;
    break;
  }
  case 17:
#line 230 "GrfnLang.y"
  {
    yyval.cval = yyvsp[-1].cval;
    ;
    break;
  }
  case 18:
#line 231 "GrfnLang.y"
  {
    Err err;
    if (err = grfn_interp_func(yyvsp[-3].sval, yyvsp[-1].cval, grfndeal)) {
      yyerror(err);
    }
    yyval.cval = yyvsp[-1].cval;
    ;
    break;
  }
  case 19:
#line 238 "GrfnLang.y"
  {
    Err err;
    if (err = grfn_interp_cell_func(yyvsp[-3].sval, yyvsp[-1].cval, grfndeal)) {
      yyerror(err);
    }
    yyval.cval = yyvsp[-1].cval;
    ;
    break;
  }
  case 20:
#line 245 "GrfnLang.y"
  {
    COMLL_PTR c, d;
    c = comll_gotobot(yyvsp[-5].cval);
    c = comll_insert_after(c);
    d = comll_gotobot(yyvsp[-3].cval);
    d = comll_insert_after(d);
    c->type = COMLL_BRANCH;
    c->next = yyvsp[-3].cval;
    yyvsp[-3].cval->prev = c;
    d->type = COMLL_JUMP;
    d->next = yyvsp[-1].cval;
    yyvsp[-1].cval->prev = d;
    c->jmp = d;
    c = comll_gotobot(yyvsp[-1].cval);
    d->jmp = c;
    yyval.cval = yyvsp[-5].cval;
    ;
    break;
  }
  case 21:
#line 262 "GrfnLang.y"
  {
    /* Added by C.Godart on 2002/11/7 */
    if (comll_atom(yyvsp[-3].cval, COMLL_REAL))
      yyerror(
          "first argument of IFFIX should be real val") if (grfndeal
                                                                ->is_history) {
        yyval.cval = yyvsp[-3].cval;
        comll_free_list(yyvsp[-1].cval);
      }
    else {
      yyval.cval = yyvsp[-1].cval;
      comll_free_list(yyvsp[-3].cval);
    };
    break;
  }
  case 22:
#line 275 "GrfnLang.y"
  {
    COMLL_PTR c;
    Err err;
    c = comll_insert_after(NULL);
    if (err = grfn_interp_name(c, yyvsp[0].sval, grfndeal)) {
      yyerror(err);
    }
    yyval.cval = c;
    ;
    break;
  }
  case 23:
#line 284 "GrfnLang.y"
  {
    COMLL_PTR c;
    c = comll_insert_after(NULL);
    c->type = COMLL_STRING;
    strncpy(c->sval, yyvsp[0].sval, COMLLBUFSZ);
    yyval.cval = c;
    ;
    break;
  }
  case 24:
#line 291 "GrfnLang.y"
  {
    COMLL_PTR c;
    c = comll_insert_after(NULL);
    c->type = COMLL_REAL;
    c->dval = yyvsp[0].dval;
    yyval.cval = c;
    ;
    break;
  }
  case 25:
#line 301 "GrfnLang.y"
  {
    yyval.cval = yyvsp[0].cval;
    yyval.cval->nargs = 1;
    ;
    break;
  }
  case 26:
#line 305 "GrfnLang.y"
  {
    yyval.cval = comll_join(yyvsp[0].cval, yyvsp[-2].cval);
    yyval.cval->nargs = yyvsp[-2].cval->nargs + 1;
    ;
    break;
  }
  }
  /* the action file gets copied in in place of this dollarsign */

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug) {
    short *ssp1 = yyss - 1;
    fprintf(stderr, "state stack now");
    while (ssp1 != yyssp)
      fprintf(stderr, " %d", *++ssp1);
    fprintf(stderr, "\n");
  }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0) {
    yylsp->first_line = yylloc.first_line;
    yylsp->first_column = yylloc.first_column;
    yylsp->last_line = (yylsp - 1)->last_line;
    yylsp->last_column = (yylsp - 1)->last_column;
    yylsp->text = 0;
  } else {
    yylsp->last_line = (yylsp + yylen - 1)->last_line;
    yylsp->last_column = (yylsp + yylen - 1)->last_column;
  }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to  ,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab: /* here on detecting error */

  if (!yyerrstatus)
  /* If not already recovering from an error  , report this error.  */
  {
    ++yynerrs;

#ifdef YYERROR_VERBOSE
    yyn = yypact[yystate];

    if (yyn > YYFLAG && yyn < YYLAST) {
      int size = 0;
      char *msg;
      int x, count;

      count = 0;
      /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
      for (x = (yyn < 0 ? -yyn : 0); x < (sizeof(yytname) / sizeof(char *));
           x++)
        if (yycheck[x + yyn] == x)
          size += strlen(yytname[x]) + 15, count++;
      msg = (char *)malloc(size + 15);
      if (msg != 0) {
        strcpy(msg, "parse error");

        if (count < 5) {
          count = 0;
          for (x = (yyn < 0 ? -yyn : 0); x < (sizeof(yytname) / sizeof(char *));
               x++)
            if (yycheck[x + yyn] == x) {
              strcat(msg, count == 0 ? "  , expecting `" : " or `");
              strcat(msg, yytname[x]);
              strcat(msg, "'");
              count++;
            }
        }
        yyerror(msg);
        free(msg);
      } else
        yyerror("parse error; also virtual memory exceeded");
    } else
#endif /* YYERROR_VERBOSE */
      yyerror("parse error");
  }

  goto yyerrlab1;
yyerrlab1: /* here on error raised explicitly by an action */

  if (yyerrstatus == 3) {
    /* if just tried and failed to reuse lookahead token after an error  ,
     * discard it.  */

    /* return failure if at end of input */
    if (yychar == YYEOF)
      YYABORT;

#if YYDEBUG != 0
    if (yydebug)
      fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

    yychar = YYEMPTY;
  }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3; /* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault
    : /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token  , ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop: /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss)
    YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug) {
    short *ssp1 = yyss - 1;
    fprintf(stderr, "Error: state stack now");
    while (ssp1 != yyssp)
      fprintf(stderr, " %d", *++ssp1);
    fprintf(stderr, "\n");
  }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0) {
    if (yyn == YYFLAG)
      goto yyerrpop;
    yyn = -yyn;
    goto yyreduce;
  } else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token  , ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;
}
#line 311 "GrfnLang.y"
