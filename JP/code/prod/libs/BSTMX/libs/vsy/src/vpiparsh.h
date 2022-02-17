/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	vpiparsh.h
 * Function:	
 * Author:	C. Daher - David Liu
 ***************************************************************/
#ifndef	_vpiparsh_H
#define	_vpiparsh_H
#include "drlstd.h"	/* Definitions and error handling */

/*--------------------------------------------------------------
 *
 */

/* init/clean the lexer */
extern	int	vpiparsLexInit(const char *s);
extern	int	vpiparsLexClean();

/* parser */
extern	int	vpiparsparse();



/*--------------------------------------------------------------
 *
 */


/*#define INSPARNSYMS	32*/	/* maximum number of symbols */
#define INSPARMAXS	256	/* maximum number of char in symbol */
#define	MAXSTRLEN	256	/* maximum string length */
#define NUMSYMMAX	256	/* maximum number of symbols */
#define	MAXCF		360		/* maximum number of cf in input */

/*#define	__DEBUG__*/

typedef	enum	{
	SYMINT,
	SYMDOUBLE,
	SYMTDATE,
	SYMLONG,
	SYMSTRING,
	SYMVPBASE,		/* KVPBase data (stream) */
	SYMDATEINT,		/* KDateInterval data */
	SYMRATE,		/* KRate data (rate definition) */
	SYMCPLXRATE,	/* KCplxRate data (complex rate)*/
	SYMINDICATOR,	/* KIndicator data */
	SYMNIL
} TSymType;




/*--------------------------------------------------------------
 * Definition of a symbol.
 */

typedef	struct {
	char		name[INSPARMAXS];	/* symbol */
	int		setFlag;		/* T=symbol is set */
	TSymType	type;			/* symbol type */

	int		ival;
	double		dval;			/* double value */
	long		lval;			/* long value */
	TDate 		dtval;			/* date value */
	TDateInterval	itval;			/* date interval value */
	TDayCount	dcval;			/* day count value */
	char		strval[INSPARMAXS];	/* string value */

	void		*ptr;			/* symbol value */
} TvpiparsSymbol;


/* Symbol table */
typedef	struct {
	int		numSym;			/* current # of symbols */
	int		numSymMax;		/* maximum # of symbols */
	TvpiparsSymbol	*table;			/* array of symbols */
} TvpiparsSymbolTable;




/*
 * Global Variables
 */
extern	TvpiparsSymbolTable	*vpiparsSymbTableG;












TvpiparsSymbolTable*	_vpiparsSymbolTableNew(
				int numSymMax);
extern	int		_vpiparsSymbolTableClear(
				TvpiparsSymbolTable *symbtab);
extern	int		_vpiparsSymbolTableLookup(
				TvpiparsSymbolTable *symbtab,
				char *s,
				TvpiparsSymbol **sym);
extern	int		_vpiparsSymbolFree(TvpiparsSymbol *symb);

extern	int		_vpiparsSymbolTablePrint(TvpiparsSymbolTable *symbtab);
extern	int		_vpiparsSymbolPrint(TvpiparsSymbol *symb);


extern	int		_vpiparsSymbolSetValue(TvpiparsSymbol *symb,
					TSymType type,
					void *value);
extern	int		_vpiparsSymbolGetValue(TvpiparsSymbol *symb,
					TSymType type,
					void *value);



extern	int		_vpiparsSetDiscName(
					char *discZcName);


/*
 * Create values for symbols.
 */



/* V-Components */

extern	int	_vpiparsSymbolSetKDateInterval(
	TvpiparsSymbol	*symb,          /* (O) */
	TDateInterval	interval,       /* (I) */
	int		dateAdjType,    /* (I) */
	char		*holidayFile,   /* (I) */
	long		badDayConv);    /* (I) */



extern	int	_vpiparsSymbolSetKRate_Fixed(
	TvpiparsSymbol *symb,		/* (O) */
	double fixedRate);

extern	int	_vpiparsSymbolSetKRate_Float1(
	TvpiparsSymbol *symb,		/* (O) */
	TDateInterval mat,		/* (I) */
	TDateInterval freq,		/* (I) */
	long dayCc,			/* (I) */
	char *curveName);		/* (I) */

extern	int	_vpiparsSymbolSetKRate_Float2(
	TvpiparsSymbol *symb,		/* (O) */
	TDateInterval mat,		/* (I) */
	TDateInterval freq,		/* (I) */
	long dayCc,			/* (I) */
	TvpiparsSymbol *spotOffSym,     /* (I) spot offset interval */
	char *curveName);		/* (I) */

extern int _vpiparsSymbolSetKRate_Cplx(
    TvpiparsSymbol *symb,       /* (O)               */
    int            numIns,      /* (I) # of instr    */
    TvpiparsSymbol **insSym,    /* (I) rate index    */
    char           *formula);   /* (I) index formula */

extern int _vpiparsSymbolSetIndicator(
    TvpiparsSymbol *symb,       /* (O)                  */
    TvpiparsSymbol *underSym,   /* (I) rate index       */
    double          lbarrier,   /* (I) lower barrier    */
    double          hbarrier,   /* (I) higher barrier   */
    char            ioWindow,   /* (I) Inside or outside*/
    char            smooth);    /* (I) Smoothing method */

/*----------------------------------------------------------------------
 */

/*
 * Utility structure used to set up schedules of arbitrary dates,
 * strikes, barries, notionals, etc. in the parser.
 * It consists of a table
 * where each column contains elements of the same type.
 */
typedef	void*	OVpiList;



/*
 * Initializes a (empty) table "oListNew" with a given number of columns
 * and no rows. The number and type of columns is given
 * by the string "types" wholse length gives the number of columns
 * and content the type of each column.
 * Available types are 'D' for dates and 'd' for doubles.
 * For example, to create an empty table with 5 columns of
 * types (date,date,date,double,double), pass the string "DDDdd".
 * Returns SUCCESS/FAILURE.
 */
int	vpiListInit(
	OVpiList *oListNew,		/* (O) new list */
	char *types);			/* (I) Ex: "DDd" for date,date,double */


/*
 * Copies a table "oListOld" into a new "oListNew" and
 * adds a new row of elements.
 * The table "oListOld" must already have been
 * initialized by vpiListInit. The character string "types" has the
 * same meaning as in vpiListInit, and is there to check the consistency
 * of the passed argument set with the table.
 * WARNING: "oListOld" is destroyed on exit.
 * For example, to add a row to a table of type "DDDdd", use
 * <PRE>
 *     OVpiList	list1, list2; 
 *     vpiListInit(&list1, "DDDdd");
 *     vpiListAdd(&list2, &list1, "DDDdd",
 *              date1, date2, date1, val1, val2);
 * </PRE>
 * Returns SUCCESS/FAILURE.
 */
int	vpiListAdd(
	OVpiList *oListNew,		/* (O) new list */
	OVpiList *oListOld,		/* (I) old list */
	char *types,			/* (I) Ex: "DDd" for date,date,double */
	...);				/* (I) all arguments */

/*
 * Frees all memory allocated for OVpiList (garbage collection).
 */
int	vpiListCleanup();



/*----------------------------------------------------------------------
 * V-Streams
 */

extern int _vpiparsSymbolSetRoot(
        TvpiparsSymbol *insSymb);               /* (I) instr symb */


extern	int	_vpiparsSymbolSetKVPWBundle(
	TvpiparsSymbol *symb,		/* (O) */
	int numIns,			/* (I) # of instr */
	TvpiparsSymbol **insSymb,	/* (I) array of instr symb */
	double *insWei);		/* (I) array of instr weight */

/*
 * Cash flow parsing routines
 */

extern	int	_vpiparsSymbolSetKVPCashFlows(
	TvpiparsSymbol *symb,		/* (O) */
	OVpiList *ocfSchedule,		/* (I) date and amount schedule */
	char *discZcName);		/* (I) curve name */



/*
 * Option parsing routines
 */

extern	int	_vpiparsSymbolSetKVPOption1(
	TvpiparsSymbol *symb,		/* (O) */
	char* type,			/* (I) call/put */
	TBoolean american,		/* (I) Amer/Euro */
	OVpiList *ooptSchedule,		/* (I) reset,acc, etc. scehdule  */
	int nDays,			/* (I) # of notifcation days */
	TvpiparsSymbol *underSym,	/* (I) underlying bundle */
	char *discZcName);		/* (I) discount curve name */

extern	int	_vpiparsSymbolSetKVPOption2(
	TvpiparsSymbol *symb,		/* (O) */
	char* type,			/* (I) call/put */
	TBoolean american,		/* (I) Amer/Euro */
	OVpiList *ooptSchedule,		/* (I) reset,acc, etc. scehdule  */
	int nDays,			/* (I) # of notifcation days */
	TvpiparsSymbol *underSym,	/* (I) underlying bundle */
	char *discZcName);		/* (I) discount curve name */

extern	int	_vpiparsSymbolSetKVPOptionSimple(
	TvpiparsSymbol *symb,		/* (O) */
	char* type,			/* (I) call/put */
	TBoolean american,		/* (I) Amer/Euro */
	TDate startDate,		/* (I) This date is not included */
	TDate matDate,			/* (I) */
	TDateInterval freq,		/* (I) */
	int nDays,			/* (I) # of notifcation days */
	double strike,			/* (I) strike */
	TvpiparsSymbol *underSym,	/* (I) underlying bundle */
	char *discZcName);	        /* (I) discount curve */

/*
 * KIO parsing routines
 */

extern 	int	_vpiparsSymbolSetKVPKnockIO1(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	char ioWindow,                  /* (I) Knock In/Out window */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */

	OVpiList *oknockSchedule,	/* (I) obs, settle, barrLo, 
					 *     barrHi, rebate scehdule  */

	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIO2(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	char ioWindow,                  /* (I) Knock In/Out window */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */

	OVpiList *oknockSchedule,	/* (I) obs, settle, barrLo, 
					 *     barrHi, rebate scehdule  */

	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIOSimple(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	char ioWindow,                  /* (I) Knock In/Out window */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */
	TDate startDate,                /* (I) This date is not included */
	TDate matDate,                  /* (I) */
	TDateInterval freq,             /* (I) */
	int nDays,                      /* (I) # of notifcation days */
 
	double barrierLo,               /* (I) Low barrier */
	double barrierHi,               /* (I) High barrier */
	double rebate,                  /* (I) Rebate */
	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIO1New(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	char ioWindow,                  /* (I) Knock In/Out window */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */

	OVpiList *oknockSchedule,	/* (I) obs, settle, barrLo, 
					 *     barrHi, rebate scehdule  */

	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIO2New(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	char ioWindow,                  /* (I) Knock In/Out window */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */

	OVpiList *oknockSchedule,	/* (I) obs, settle, barrLo, 
					 *     barrHi, rebate scehdule  */

	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIOSimpleNew(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
	char smooth,              	/* (I) Node smoothing method */
	TDate startDate,                /* (I) This date is not included */
	TDate matDate,                  /* (I) */
	TDateInterval freq,             /* (I) */
	int nDays,                      /* (I) # of notifcation days */
 
	double rebate,                  /* (I) Rebate */
	TvpiparsSymbol *underSym,       /* (I) underlying bundle */
	char *discZcName);

extern 	int	_vpiparsSymbolSetKVPKnockIO2Idx(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out              */
	char ioWindow1,                 /* (I) Knock In/Out window 1     */
	TvpiparsSymbol *rateSym1,       /* (I) Knock In/Out rate index 1 */
    char ioWindow2,                 /* (I) Knock In/Out window 2     */
	TvpiparsSymbol *rateSym2,       /* (I) Knock In/Out rate index 2 */
	char smooth,              	    /* (I) Node smoothing method     */
	OVpiList *oknockSchedule,	    /* (I) obs, settle, barrLo, 
					                 *     barrHi, rebate schedule   */
	TvpiparsSymbol *underSym,       /* (I) underlying bundle         */
	char *discZcName);


/*
 * Floating leg parsing routines:
 * the structure TFlpList is used  in the parser generated C file
 * to construct arbitrary schedules.
 */

extern	int	_vpiparsSymbolSetKVPFloatLegSimple(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	TDate startDate,		/* (I) start date */
	TDate matDate,			/* (I) maturity date */
	TDateInterval freq,		/* (I) frequency */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg1(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg2(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg3(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg4(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. schedule 
					 *     including formula string */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *discZcName);	        /* (I) discount curve */

/* 
 * Float Leg with a set of multiple index rates 
 */

extern	int	_vpiparsSymbolSetKVPFloatLegNewSimple(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	TDate startDate,		/* (I) start date */
	TDate matDate,			/* (I) maturity date */
	TDateInterval freq,		/* (I) frequency */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLegNew1(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLegNew2(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLegNew3(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLegNew4(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. schedule 
					 *     including formula string */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLegNew5(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym,	/* (I) rate index */
	char *discZcName);	        /* (I) discount curve */

/* 
 * Float Leg w/ 2 set of index rates 
 */
extern	int	_vpiparsSymbolSetKVPFloatLeg2IdxSimple(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	TDate startDate,		/* (I) start date */
	TDate matDate,			/* (I) end date */
	TDateInterval freq,		/* (I) frequency */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym1,	/* (I) rate index 1 */
	TvpiparsSymbol *rateSym,	/* (I) rate index 2 */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg2Idx(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, etc. scehdule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym1,	/* (I) rate index 1 */
	TvpiparsSymbol *rateSym2,	/* (I) rate index 2 */
	char *formula,			/* (I) payment formula */
	char *discZcName);	        /* (I) discount curve */

extern	int	_vpiparsSymbolSetKVPFloatLeg2Idx2(
	TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, formula schedule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym1,	/* (I) rate index 1 */
	TvpiparsSymbol *rateSym2,	/* (I) rate index 2 */
	char *discZcName);	        /* (I) discount curve */

/*
 * Float leg with 3 index rates
 */
extern	int	_vpiparsSymbolSetKVPFloatLeg3Idx(
        TvpiparsSymbol *symb,		/* (O) new symbol */
	OVpiList *oflpSchedule,		/* (I) reset,acc, formula schedule  */
	long dayCc,			/* (I) day count convention */
	long stubConv,			/* (I) stub convention */
	TvpiparsSymbol *rateSym1,	/* (I) rate index 1 */
	TvpiparsSymbol *rateSym2,	/* (I) rate index 2 */
	TvpiparsSymbol *rateSym3,	/* (I) rate index 3 */
	char *discZcName);	        /* (I) discount curve */


/*
 * Protection leg 
 */
extern int     _vpiparsSymbolSetKVPProtLegSimple(
        TvpiparsSymbol *symb,        /* (O) */
        TDate startDate,             /* (I) start date  */
        TDate endDate,               /* (I) end date */
        double notional,             /* (I) notional */
        char   *recovery,             /* (I) recovery */
        long   defConv,              /* (I) stub convention */
        char *discZcName);           /* (I) CDS curve */

/*
 * Default protection simple schedule 
 */
extern	int	_vpiparsSymbolSetKVPDefExposureSimple(
	TvpiparsSymbol *symb,		/* (O) */
	TDate startDate,		/* (I) Start date */
	TDate endDate,			/* (I) End date */
	TDateInterval freq,		/* (I) freq     */
        double   recovery,              /* (I) recovery rate       */
	TvpiparsSymbol *underSym,	/* (I) underlying bundle */
	char *discZcName);	        /* (I) discount curve */


/*
 * Default protection general schedule 
 */
extern	int	_vpiparsSymbolSetKVPDefExposureGeneral(
	TvpiparsSymbol *symb,		/* (O) */
	OVpiList *oDefProtSchedule,	/* (I) start, end, settle scehdule  */
        double   recovery,              /* (I) recovery rate       */
	TvpiparsSymbol *underSym,	/* (I) underlying */
	char *discZcName);		/* (I) discount curve name */

/*
 * Default protection general schedule with rebate
 */
extern	int	_vpiparsSymbolSetKVPDefKnockIn(
	TvpiparsSymbol *symb,		/* (O) */
	OVpiList *oDefProtSchedule,	/* (I) start, end, settle, rebate */
	TvpiparsSymbol *underSym,	/* (I) underlying */
	char *discZcName);		/* (I) discount curve name */


/*
 * Reset bank
 */

extern	int	_vpiparsAddRateResetRank(
	TvpiparsSymbol *rateSym,	/* (I) underlying bundle */
	TDate resetDate,		/* (I) */
	double value);			/* (I) */



#endif	/* _vpiparsh_H */


