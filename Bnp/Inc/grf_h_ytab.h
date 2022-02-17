typedef union{
	double dval;
	char sval[300];
	COMLL_PTR cval;
} YYSTYPE;
#define	QNAME	258
#define	NAME	259
#define	NUMBER	260
#define	ASSIGN	261
#define	LE	262
#define	GE	263
#define	AND	264
#define	OR	265
#define	IF	266
#define	EQ	267
#define	DO	268
#define	IFFIX	269
#define	UMINUS	270


extern YYSTYPE yylval;
