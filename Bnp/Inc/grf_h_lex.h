/*-----------------------------------------------------------------
  DESCRIPTION     :  Grfn include file with all the yacc and lex functions as
well.
---------------------------------------------------------------------*/

#ifndef GRF_H_LEX_H
#define GRF_H_LEX_H

/* ========================== LEX AND YACC ============================== */

/**** LEX USES THIS STRING AS ITS INPUT ***/
EXTERN String grfn_input_string;

/*** TELL LEX TO QUIT PARSING WHEN IT REACHES THE END OF GRFN_INPUT_STRING **/
#define yywrap() 1

#define GRFN_YYPARSE_DECL Err yyparse(COMLL_PTR *ret_ptr, GrfnDeal *grfndeal)
GRFN_YYPARSE_DECL;

Err grfn_yyerror(const String fmt, GrfnDeal *gd);

/*** WHAT LEX AND YACC CALL IF THERE IS AN ERROR **/
#define yyerror(fmt)                                                           \
  { return grfn_yyerror(fmt, grfndeal); }

/*--- declarations of lexer ------*/
int yylex(void);

#endif
