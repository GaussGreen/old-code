#include "edginc/config.hpp"/* A Bison parser, made by GNU Bison 1.875b.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     DEXP = 258,
     DEXPARRAY = 259,
     DTEXP = 260,
     BEXP = 261,
     BEXPARRAY = 262,
     IEXP = 263,
     IEXPARRAY = 264,
     SCHEDEXP = 265,
     TABFUNCEXP = 266,
     MAX = 267,
     MIN = 268,
     IFUNC = 269,
     IARRAYFUNC = 270,
     BFUNC = 271,
     BARRAYFUNC = 272,
     DFUNC = 273,
     DARRAYFUNC = 274,
     DTFUNC = 275,
     SUBARRAY = 276,
     DVAR = 277,
     DVARARRAY = 278,
     IVAR = 279,
     IVARARRAY = 280,
     BVAR = 281,
     BVARARRAY = 282,
     DTVAR = 283,
     EQUALS = 284,
     NEQUALS = 285,
     UNDEFINED_SYMBOL = 286,
     INTERNAL_ERROR = 287,
     OR = 288,
     AND = 289,
     NEG = 290
   };
#endif
#define DEXP 258
#define DEXPARRAY 259
#define DTEXP 260
#define BEXP 261
#define BEXPARRAY 262
#define IEXP 263
#define IEXPARRAY 264
#define SCHEDEXP 265
#define TABFUNCEXP 266
#define MAX 267
#define MIN 268
#define IFUNC 269
#define IARRAYFUNC 270
#define BFUNC 271
#define BARRAYFUNC 272
#define DFUNC 273
#define DARRAYFUNC 274
#define DTFUNC 275
#define SUBARRAY 276
#define DVAR 277
#define DVARARRAY 278
#define IVAR 279
#define IVARARRAY 280
#define BVAR 281
#define BVARARRAY 282
#define DTVAR 283
#define EQUALS 284
#define NEQUALS 285
#define UNDEFINED_SYMBOL 286
#define INTERNAL_ERROR 287
#define OR 288
#define AND 289
#define NEG 290




/* Copy the first part of user declarations.  */
#line 52 "FRParser.ypp"

// note: perl script adds include of config.hpp
#define EDR_FRPARSER_SOURCE_FILE // allows us to see full header file
#include "edginc/FRParseException.hpp"
#include "edginc/FRParser.hpp"

DRLIB_BEGIN_NAMESPACE
class FRParser;

#define YYERROR_VERBOSE
void  yyerror (const char *s);
static int yylex(void* lvalp, void* parser);
int yyparse(void* param);

#define YYPARSE_PARAM param
#define YYLEX_PARAM param
#define MY_PARAM (*((FRParser::Params*) param))
#define TRY try{
#define CATCH } catch (exception& e) { \
    /* save the exception */ \
    (*((FRParser::Params*) param)).ourError = ModelException(e); \
    YYERROR; /* revert back to C style error handling */ \
}
// macro to ease memory management and error trapping
#define CHECK_PTR(__p) { \
    MY_PARAM.theCtrl->store(__p); \
    if ((__p) == 0){ YYERROR;} \
}

/* stores pointer passed in, checks it's not null and terminates
   earlier TRY statement */
#define CHK_CATCH(__p) CHECK_PTR(__p) CATCH

// parameters to the parser
class FRParser::Params{
public:
    const char*   thePos; // where we are in the string
    const char*   prevPos; // where we were in the string last time
    FRController* theCtrl;
    FRIfaces::IRValue* finalVal;// what the expression evualuates to
    ModelException     ourError;  /* any errors that occur in our code - ie
                                     yylex or code per rule */

    Params(const char* expression, FRController* theCtrl):
        thePos(expression), prevPos(thePos), theCtrl(theCtrl), finalVal(0){}
};
static string errorMsg; // from yyerror() = can't pass Params to it

/** creates FlexRule::IRValue which represents the value
    of this expression at the specific timepoint. The getRValue()
    uses this method and then saves it in the FRController. Return
    null if the save has already been done */
FRIfaces::IRValue* FRParser::createRValue(
    int           index,
    FRController* frCtrl) const{
    if (index != frCtrl->getIndex()){
        throw ModelException("FRParser::createRValue", "Internal error: "
                             "invoking parser for expression at some other "
                             "timepoint");
    }
    return parseRValue(expression.c_str(), frCtrl);
}

/** Parses supplied string. Throws FRParseException for parse type errors.
    Object returned needs to be freed */
FRIfaces::IRValue* FRParser::parseRValue(
    const char*   expression,
    FRController* frCtrl){
    static const string routine("FRParser::createRValue");
    Params params(expression, frCtrl);
    if (yyparse((void*)&params) != 0){
        string err1("Failed parsing: "+string(expression));
        string err2("                "+
                    string(params.prevPos-expression, ' ')+"^");
        FRParseException* cause = 
            dynamic_cast<FRParseException*>(params.ourError.getCause());
        // 
        if ((params.ourError.empty() && !errorMsg.empty()) || cause){
            string what;
            if (cause){
                what = string(cause->what());
                if (!errorMsg.empty()){
                    what += '\n';
                }
            }
            string m(what+errorMsg+"\n"+err1+"\n"+err2);
            errorMsg.erase(); // tidy up
            throw FRParseException(m);
        }
        errorMsg.erase(); // tidy up
        params.ourError.addMsg(err1);
        params.ourError.addMsg(err2);
        params.ourError.addMsg(routine+": Failed");
        throw params.ourError;
    }
    return params.finalVal;
}
    
/** returns the id for this RValueBase */
const string& FRParser::getID() const{
    return expression;
}

FRParser::FRParser(): RValueBase(TYPE) {}

IObject* FRParser::defaultConstructor(){
    return new FRParser();
}

// these have to be outside of the class - c++ is just so rubbish at times
CClassConstSP const FRParser::TYPE = CClass::registerClassLoadMethod(
    "FRParser", typeid(FRParser), load);

/*
  Arrays: want arrays of doubles, ints and bools.
  Syntax: myArray[time idx][array index] is a double/int etc
  myArray[time idx] is an array
  myArray is equivalent to my[idx]
  Rules: array+array, array+scalar, array*array, array*scalar, 
         array-array, array-scalar, array/scalar, 
         MAX(array, array), MAX(array, scalar) repeat for MIN
  Functions: need a sort which returns an int array (operating on an db array)
             need a sort which also takes in this int array, plus start, end
             need a sum/rainbow which takes in  "  "  "      "    "      "
  Should we be able to create a bool array containing 'picked' etc? Not obvious
  that we do.

  Have dVarArray. Have '[]' rule which is?
  If another dVarArray then won't be able to dereference it.
  What if it's a dExpArray? Sounds ok but it is making the assumption that the
  only way we get 2D expressions is from a variable.
  What happens is tokenizer returns different type for plain variables and for
  variables followed by '['. So we return DEXPARRAY for the former and 
  DVARARRAY for the latter.


Need to:
1. Ensure tokeniser supports them
2. Do rules/code for array components, array+array and array+scalar.
(i) Need to update RValUnion and all corresponding code. Done.
Check code compiles.
(ii) Consequently need array version of IndexOperator. (And so need
array version of invokeGetValue). More importantly need array versions of
unary and binary operator.
(iii) Need to review what getValue() returns for arrays. Alter TGetValue
      prototype so that it is 
      typedef double (TGetValue)(void* structure, int index);
      But what about the length of the array? - this is needed! Could put it
      in RT structure as first element (so common throughout arrays?)
      Avoid worrying about performance for now (eg perhaps could have optional
      method that returns double* - but implementations would have to switch)
      Need to sort out impact of this change (mainly FRUtil I think)

3. Allow them (db array) to be created (using an array of variable names). Test
4. Add this support for ints. Test
5. Write sort functions. Test.
6. Add other rules. Test
 */


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 216 "FRParser.ypp"
typedef union YYSTYPE {
    FRIfaces::IRValueDouble* dExp;  // for r-values evaluating to doubles
    FR::RDoubleArray*        dValArray; /* for on the fly arrays of r 
                                           double vals */
    FRIfaces::IRValueDoubleArray* dExpArray; /* for general db arrays */
    FRIfaces::IRValueBool*   bExp;  // for r-values evaluating to bools
    FR::RBoolArray*          bValArray; /* for on the fly arrays of r 
                                           bool vals */
    FRIfaces::IRValueBoolArray* bExpArray; /* for general db arrays */
    FRIfaces::IRValueInt*    iExp;  // for r-values evaluating to ints
    FR::RIntArray*           iValArray; /* for on the fly arrays of r 
                                           int vals */
    FRIfaces::IRValueIntArray* iExpArray; /* for general int arrays */
    FRIfaces::IRValueDate*   dtExp; // for r-values evaluating to dates
    FRIfaces::IRValueSchedule* schedExp; // for schedules
    FRIfaces::IRValueTabulatedFunc* tabFuncExp; // 'tabulated function'
    const FRFunction*        func;  // functions (various signatures)
    FRInternalFunction*      funcArgs;  // functions (various signatures)
    FRIfaces::ILValueExpression* var; // variable before '[]' operator
} YYSTYPE;
/* Line 191 of yacc.c.  */
#line 329 "FRParser.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 341 "FRParser.cpp"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  97
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   4081

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  54
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  39
/* YYNRULES -- Number of rules. */
#define YYNRULES  389
/* YYNRULES -- Number of states. */
#define YYNSTATES  755

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   290

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    45,     2,     2,     2,     2,     2,     2,
      48,    49,    41,    39,    51,    38,    40,    42,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    50,     2,
      36,     2,    37,    33,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    46,     2,    47,    44,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    52,     2,    53,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    34,    35,
      43
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     5,     7,     9,    11,    13,    15,    17,
      19,    21,    23,    28,    33,    36,    40,    44,    48,    52,
      56,    60,    64,    68,    72,    76,    80,    84,    88,    92,
      96,   100,   104,   110,   112,   114,   119,   124,   128,   132,
     136,   140,   143,   147,   151,   157,   164,   171,   173,   177,
     182,   187,   189,   194,   199,   203,   207,   211,   215,   219,
     223,   227,   231,   235,   239,   243,   247,   251,   255,   259,
     262,   266,   273,   280,   287,   294,   301,   308,   314,   320,
     326,   328,   332,   337,   342,   344,   349,   353,   357,   361,
     368,   375,   381,   383,   385,   387,   389,   391,   393,   395,
     397,   399,   401,   403,   406,   410,   414,   417,   420,   424,
     427,   430,   434,   437,   439,   441,   446,   450,   454,   458,
     462,   466,   470,   474,   478,   482,   486,   490,   494,   498,
     502,   505,   507,   513,   519,   525,   531,   540,   544,   546,
     548,   553,   557,   561,   565,   569,   573,   577,   581,   585,
     589,   592,   594,   600,   606,   612,   618,   627,   631,   633,
     635,   638,   643,   645,   649,   653,   657,   661,   665,   669,
     673,   677,   681,   685,   689,   693,   697,   701,   705,   709,
     715,   721,   727,   733,   742,   746,   749,   753,   757,   761,
     765,   769,   773,   777,   781,   785,   791,   797,   803,   809,
     813,   817,   821,   825,   829,   833,   837,   841,   845,   851,
     857,   863,   869,   872,   876,   880,   884,   888,   892,   896,
     900,   904,   908,   914,   920,   926,   932,   936,   940,   944,
     948,   952,   956,   960,   964,   968,   974,   980,   986,   992,
     995,   999,  1003,  1007,  1011,  1015,  1019,  1023,  1027,  1031,
    1037,  1043,  1049,  1055,  1059,  1063,  1067,  1071,  1075,  1079,
    1083,  1087,  1091,  1097,  1103,  1109,  1115,  1118,  1122,  1126,
    1130,  1134,  1138,  1142,  1146,  1150,  1154,  1160,  1166,  1172,
    1178,  1184,  1188,  1192,  1196,  1200,  1204,  1208,  1212,  1216,
    1220,  1226,  1232,  1238,  1244,  1250,  1256,  1262,  1265,  1269,
    1273,  1277,  1281,  1285,  1289,  1293,  1297,  1301,  1307,  1313,
    1319,  1325,  1331,  1337,  1343,  1347,  1351,  1355,  1359,  1363,
    1367,  1371,  1375,  1379,  1385,  1391,  1397,  1403,  1409,  1415,
    1421,  1424,  1428,  1432,  1436,  1440,  1444,  1448,  1452,  1456,
    1460,  1466,  1472,  1478,  1484,  1490,  1496,  1502,  1506,  1510,
    1514,  1518,  1522,  1526,  1530,  1534,  1538,  1544,  1550,  1556,
    1562,  1568,  1574,  1580,  1583,  1587,  1591,  1595,  1599,  1603,
    1607,  1611,  1615,  1619,  1625,  1631,  1637,  1643,  1647,  1651,
    1655,  1659,  1663,  1667,  1671,  1675,  1679,  1685,  1691,  1697
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      55,     0,    -1,    56,    -1,    57,    -1,    58,    -1,    59,
      -1,    76,    -1,    77,    -1,    78,    -1,    31,    -1,    32,
      -1,     6,    -1,    64,    46,    57,    47,    -1,    78,    46,
      57,    47,    -1,    45,    56,    -1,    56,    29,    56,    -1,
      56,    30,    56,    -1,    48,    56,    49,    -1,    57,    29,
      57,    -1,    57,    30,    57,    -1,    58,    29,    58,    -1,
      58,    30,    58,    -1,    59,    29,    59,    -1,    59,    30,
      59,    -1,    57,    36,    57,    -1,    57,    37,    57,    -1,
      58,    36,    58,    -1,    58,    37,    58,    -1,    59,    36,
      59,    -1,    59,    37,    59,    -1,    56,    35,    56,    -1,
      56,    34,    56,    -1,    56,    33,    56,    50,    56,    -1,
      80,    -1,     8,    -1,    63,    46,    57,    47,    -1,    77,
      46,    57,    47,    -1,    57,    39,    57,    -1,    57,    38,
      57,    -1,    57,    41,    57,    -1,    57,    42,    57,    -1,
      38,    57,    -1,    48,    57,    49,    -1,    59,    38,    59,
      -1,    56,    33,    57,    50,    57,    -1,    12,    48,    57,
      51,    57,    49,    -1,    13,    48,    57,    51,    57,    49,
      -1,    82,    -1,    77,    40,    77,    -1,    12,    48,    77,
      49,    -1,    13,    48,    77,    49,    -1,     3,    -1,    62,
      46,    57,    47,    -1,    76,    46,    57,    47,    -1,    58,
      39,    58,    -1,    58,    39,    57,    -1,    57,    39,    58,
      -1,    58,    38,    58,    -1,    58,    38,    57,    -1,    57,
      38,    58,    -1,    58,    41,    58,    -1,    58,    41,    57,
      -1,    57,    41,    58,    -1,    58,    42,    58,    -1,    58,
      42,    57,    -1,    57,    42,    58,    -1,    58,    44,    58,
      -1,    58,    44,    57,    -1,    57,    44,    58,    -1,    38,
      58,    -1,    48,    58,    49,    -1,    12,    48,    58,    51,
      58,    49,    -1,    12,    48,    57,    51,    58,    49,    -1,
      12,    48,    58,    51,    57,    49,    -1,    13,    48,    58,
      51,    58,    49,    -1,    13,    48,    57,    51,    58,    49,
      -1,    13,    48,    58,    51,    57,    49,    -1,    56,    33,
      58,    50,    58,    -1,    56,    33,    57,    50,    58,    -1,
      56,    33,    58,    50,    57,    -1,    84,    -1,    76,    40,
      76,    -1,    12,    48,    76,    49,    -1,    13,    48,    76,
      49,    -1,     5,    -1,    65,    46,    57,    47,    -1,    59,
      39,    57,    -1,    59,    38,    57,    -1,    48,    59,    49,
      -1,    12,    48,    59,    51,    59,    49,    -1,    13,    48,
      59,    51,    59,    49,    -1,    56,    33,    59,    50,    59,
      -1,    92,    -1,    10,    -1,    11,    -1,    22,    -1,    24,
      -1,    26,    -1,    28,    -1,    23,    -1,    25,    -1,    27,
      -1,    52,    -1,    69,    58,    -1,    70,    51,    58,    -1,
      70,    51,    57,    -1,    70,    53,    -1,    69,    57,    -1,
      72,    51,    57,    -1,    72,    53,    -1,    69,    56,    -1,
      74,    51,    56,    -1,    74,    53,    -1,    71,    -1,     4,
      -1,    66,    46,    57,    47,    -1,    76,    39,    58,    -1,
      76,    39,    57,    -1,    76,    39,    76,    -1,    76,    38,
      58,    -1,    76,    38,    57,    -1,    76,    38,    76,    -1,
      76,    41,    58,    -1,    76,    41,    57,    -1,    58,    41,
      76,    -1,    57,    41,    76,    -1,    76,    41,    76,    -1,
      76,    42,    58,    -1,    76,    42,    57,    -1,    76,    42,
      76,    -1,    38,    76,    -1,    86,    -1,    78,    33,    76,
      50,    76,    -1,    78,    33,    76,    50,    58,    -1,    78,
      33,    58,    50,    76,    -1,    56,    33,    76,    50,    76,
      -1,    21,    48,    76,    51,    57,    51,    57,    49,    -1,
      48,    76,    49,    -1,    73,    -1,     9,    -1,    67,    46,
      57,    47,    -1,    77,    39,    57,    -1,    77,    39,    77,
      -1,    77,    38,    57,    -1,    77,    38,    77,    -1,    77,
      41,    57,    -1,    57,    41,    77,    -1,    77,    41,    77,
      -1,    77,    42,    57,    -1,    77,    42,    77,    -1,    38,
      77,    -1,    88,    -1,    78,    33,    77,    50,    77,    -1,
      78,    33,    77,    50,    57,    -1,    78,    33,    57,    50,
      77,    -1,    56,    33,    77,    50,    77,    -1,    21,    48,
      77,    51,    57,    51,    57,    49,    -1,    48,    77,    49,
      -1,    75,    -1,     7,    -1,    45,    78,    -1,    68,    46,
      57,    47,    -1,    90,    -1,    76,    37,    76,    -1,    76,
      36,    76,    -1,    76,    37,    58,    -1,    76,    36,    58,
      -1,    77,    37,    77,    -1,    77,    36,    77,    -1,    77,
      37,    57,    -1,    77,    36,    57,    -1,    78,    34,    78,
      -1,    78,    35,    78,    -1,    78,    29,    78,    -1,    78,
      30,    78,    -1,    77,    29,    77,    -1,    77,    30,    77,
      -1,    76,    29,    76,    -1,    76,    30,    76,    -1,    78,
      33,    78,    50,    78,    -1,    78,    33,    78,    50,    56,
      -1,    78,    33,    56,    50,    78,    -1,    56,    33,    78,
      50,    78,    -1,    21,    48,    78,    51,    57,    51,    57,
      49,    -1,    48,    78,    49,    -1,    16,    48,    -1,    79,
      56,    51,    -1,    79,    78,    51,    -1,    79,    57,    51,
      -1,    79,    77,    51,    -1,    79,    58,    51,    -1,    79,
      76,    51,    -1,    79,    59,    51,    -1,    79,    60,    51,
      -1,    79,    61,    51,    -1,    79,    64,    46,    47,    51,
      -1,    79,    63,    46,    47,    51,    -1,    79,    62,    46,
      47,    51,    -1,    79,    65,    46,    47,    51,    -1,    79,
      56,    49,    -1,    79,    78,    49,    -1,    79,    57,    49,
      -1,    79,    77,    49,    -1,    79,    58,    49,    -1,    79,
      76,    49,    -1,    79,    59,    49,    -1,    79,    60,    49,
      -1,    79,    61,    49,    -1,    79,    64,    46,    47,    49,
      -1,    79,    63,    46,    47,    49,    -1,    79,    62,    46,
      47,    49,    -1,    79,    65,    46,    47,    49,    -1,    14,
      48,    -1,    81,    56,    51,    -1,    81,    78,    51,    -1,
      81,    57,    51,    -1,    81,    77,    51,    -1,    81,    58,
      51,    -1,    81,    76,    51,    -1,    81,    59,    51,    -1,
      81,    60,    51,    -1,    81,    61,    51,    -1,    81,    64,
      46,    47,    51,    -1,    81,    63,    46,    47,    51,    -1,
      81,    62,    46,    47,    51,    -1,    81,    65,    46,    47,
      51,    -1,    81,    56,    49,    -1,    81,    78,    49,    -1,
      81,    57,    49,    -1,    81,    77,    49,    -1,    81,    58,
      49,    -1,    81,    76,    49,    -1,    81,    59,    49,    -1,
      81,    60,    49,    -1,    81,    61,    49,    -1,    81,    64,
      46,    47,    49,    -1,    81,    63,    46,    47,    49,    -1,
      81,    62,    46,    47,    49,    -1,    81,    65,    46,    47,
      49,    -1,    18,    48,    -1,    83,    56,    51,    -1,    83,
      78,    51,    -1,    83,    57,    51,    -1,    83,    77,    51,
      -1,    83,    58,    51,    -1,    83,    76,    51,    -1,    83,
      59,    51,    -1,    83,    60,    51,    -1,    83,    61,    51,
      -1,    83,    64,    46,    47,    51,    -1,    83,    63,    46,
      47,    51,    -1,    83,    62,    46,    47,    51,    -1,    83,
      65,    46,    47,    51,    -1,    83,    56,    49,    -1,    83,
      78,    49,    -1,    83,    57,    49,    -1,    83,    77,    49,
      -1,    83,    58,    49,    -1,    83,    76,    49,    -1,    83,
      59,    49,    -1,    83,    60,    49,    -1,    83,    61,    49,
      -1,    83,    64,    46,    47,    49,    -1,    83,    63,    46,
      47,    49,    -1,    83,    62,    46,    47,    49,    -1,    83,
      65,    46,    47,    49,    -1,    19,    48,    -1,    85,    56,
      51,    -1,    85,    78,    51,    -1,    85,    57,    51,    -1,
      85,    77,    51,    -1,    85,    58,    51,    -1,    85,    76,
      51,    -1,    85,    59,    51,    -1,    85,    60,    51,    -1,
      85,    61,    51,    -1,    85,    64,    46,    47,    51,    -1,
      85,    63,    46,    47,    51,    -1,    85,    62,    46,    47,
      51,    -1,    85,    66,    46,    47,    51,    -1,    85,    65,
      46,    47,    51,    -1,    85,    56,    49,    -1,    85,    78,
      49,    -1,    85,    57,    49,    -1,    85,    77,    49,    -1,
      85,    58,    49,    -1,    85,    76,    49,    -1,    85,    59,
      49,    -1,    85,    60,    49,    -1,    85,    61,    49,    -1,
      85,    64,    46,    47,    49,    -1,    85,    68,    46,    47,
      49,    -1,    85,    63,    46,    47,    49,    -1,    85,    67,
      46,    47,    49,    -1,    85,    62,    46,    47,    49,    -1,
      85,    66,    46,    47,    49,    -1,    85,    65,    46,    47,
      49,    -1,    15,    48,    -1,    87,    56,    51,    -1,    87,
      78,    51,    -1,    87,    57,    51,    -1,    87,    77,    51,
      -1,    87,    58,    51,    -1,    87,    76,    51,    -1,    87,
      59,    51,    -1,    87,    60,    51,    -1,    87,    61,    51,
      -1,    87,    64,    46,    47,    51,    -1,    87,    68,    46,
      47,    51,    -1,    87,    63,    46,    47,    51,    -1,    87,
      67,    46,    47,    51,    -1,    87,    62,    46,    47,    51,
      -1,    87,    66,    46,    47,    51,    -1,    87,    65,    46,
      47,    51,    -1,    87,    56,    49,    -1,    87,    78,    49,
      -1,    87,    57,    49,    -1,    87,    77,    49,    -1,    87,
      58,    49,    -1,    87,    76,    49,    -1,    87,    59,    49,
      -1,    87,    60,    49,    -1,    87,    61,    49,    -1,    87,
      64,    46,    47,    49,    -1,    87,    68,    46,    47,    49,
      -1,    87,    63,    46,    47,    49,    -1,    87,    67,    46,
      47,    49,    -1,    87,    62,    46,    47,    49,    -1,    87,
      66,    46,    47,    49,    -1,    87,    65,    46,    47,    49,
      -1,    17,    48,    -1,    89,    56,    51,    -1,    89,    78,
      51,    -1,    89,    57,    51,    -1,    89,    77,    51,    -1,
      89,    58,    51,    -1,    89,    76,    51,    -1,    89,    59,
      51,    -1,    89,    60,    51,    -1,    89,    61,    51,    -1,
      89,    64,    46,    47,    51,    -1,    89,    68,    46,    47,
      51,    -1,    89,    63,    46,    47,    51,    -1,    89,    67,
      46,    47,    51,    -1,    89,    62,    46,    47,    51,    -1,
      89,    66,    46,    47,    51,    -1,    89,    65,    46,    47,
      51,    -1,    89,    56,    49,    -1,    89,    78,    49,    -1,
      89,    57,    49,    -1,    89,    77,    49,    -1,    89,    58,
      49,    -1,    89,    76,    49,    -1,    89,    59,    49,    -1,
      89,    60,    49,    -1,    89,    61,    49,    -1,    89,    64,
      46,    47,    49,    -1,    89,    68,    46,    47,    49,    -1,
      89,    63,    46,    47,    49,    -1,    89,    67,    46,    47,
      49,    -1,    89,    62,    46,    47,    49,    -1,    89,    66,
      46,    47,    49,    -1,    89,    65,    46,    47,    49,    -1,
      20,    48,    -1,    91,    56,    51,    -1,    91,    78,    51,
      -1,    91,    57,    51,    -1,    91,    77,    51,    -1,    91,
      58,    51,    -1,    91,    76,    51,    -1,    91,    59,    51,
      -1,    91,    60,    51,    -1,    91,    61,    51,    -1,    91,
      64,    46,    47,    51,    -1,    91,    63,    46,    47,    51,
      -1,    91,    62,    46,    47,    51,    -1,    91,    65,    46,
      47,    51,    -1,    91,    56,    49,    -1,    91,    78,    49,
      -1,    91,    57,    49,    -1,    91,    77,    49,    -1,    91,
      58,    49,    -1,    91,    76,    49,    -1,    91,    59,    49,
      -1,    91,    60,    49,    -1,    91,    61,    49,    -1,    91,
      64,    46,    47,    49,    -1,    91,    63,    46,    47,    49,
      -1,    91,    62,    46,    47,    49,    -1,    91,    65,    46,
      47,    49,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   299,   299,   300,   301,   302,   303,   304,   305,   306,
     307,   310,   311,   314,   319,   323,   325,   327,   328,   330,
     332,   335,   338,   341,   344,   346,   348,   351,   354,   357,
     360,   362,   364,   366,   371,   372,   376,   381,   383,   385,
     387,   389,   393,   394,   401,   403,   405,   407,   409,   418,
     425,   435,   436,   440,   445,   447,   450,   454,   456,   459,
     463,   465,   468,   472,   474,   477,   481,   483,   485,   488,
     492,   494,   496,   498,   501,   503,   505,   508,   510,   513,
     516,   518,   527,   534,   544,   545,   548,   551,   554,   555,
     559,   563,   567,   572,   575,   578,   581,   584,   587,   590,
     592,   594,   605,   609,   612,   614,   621,   625,   628,   632,
     636,   639,   643,   648,   649,   650,   654,   660,   666,   672,
     678,   684,   690,   696,   702,   708,   714,   720,   726,   732,
     738,   743,   745,   750,   756,   762,   767,   772,   777,   778,
     779,   783,   789,   795,   801,   807,   813,   819,   825,   831,
     837,   842,   844,   849,   855,   861,   866,   871,   876,   877,
     878,   883,   887,   889,   895,   901,   907,   913,   919,   925,
     931,   937,   940,   943,   949,   955,   961,   967,   973,   979,
     984,   990,   996,  1001,  1006,  1011,  1017,  1018,  1019,  1020,
    1021,  1022,  1023,  1024,  1025,  1026,  1027,  1028,  1029,  1037,
    1038,  1039,  1040,  1041,  1042,  1043,  1044,  1045,  1046,  1047,
    1048,  1049,  1058,  1064,  1065,  1066,  1067,  1068,  1069,  1070,
    1071,  1072,  1073,  1074,  1075,  1076,  1084,  1085,  1086,  1087,
    1088,  1089,  1090,  1091,  1092,  1093,  1094,  1095,  1096,  1105,
    1111,  1112,  1113,  1114,  1115,  1116,  1117,  1118,  1119,  1120,
    1121,  1122,  1123,  1131,  1132,  1133,  1134,  1135,  1136,  1137,
    1138,  1139,  1140,  1141,  1142,  1143,  1152,  1158,  1159,  1160,
    1161,  1162,  1163,  1164,  1165,  1166,  1167,  1168,  1169,  1170,
    1172,  1180,  1181,  1182,  1183,  1184,  1185,  1186,  1187,  1188,
    1189,  1190,  1192,  1193,  1195,  1196,  1198,  1207,  1213,  1214,
    1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,  1223,  1225,
    1226,  1228,  1229,  1231,  1239,  1240,  1241,  1242,  1243,  1244,
    1245,  1246,  1247,  1248,  1249,  1251,  1252,  1254,  1255,  1257,
    1266,  1272,  1273,  1274,  1275,  1276,  1277,  1278,  1279,  1280,
    1281,  1282,  1284,  1285,  1287,  1288,  1290,  1298,  1299,  1300,
    1301,  1302,  1303,  1304,  1305,  1306,  1307,  1308,  1310,  1311,
    1313,  1314,  1316,  1325,  1330,  1331,  1332,  1333,  1334,  1335,
    1336,  1337,  1338,  1339,  1340,  1341,  1342,  1350,  1351,  1352,
    1353,  1354,  1355,  1356,  1357,  1358,  1359,  1360,  1361,  1362
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "DEXP", "DEXPARRAY", "DTEXP", "BEXP", 
  "BEXPARRAY", "IEXP", "IEXPARRAY", "SCHEDEXP", "TABFUNCEXP", "MAX", 
  "MIN", "IFUNC", "IARRAYFUNC", "BFUNC", "BARRAYFUNC", "DFUNC", 
  "DARRAYFUNC", "DTFUNC", "SUBARRAY", "DVAR", "DVARARRAY", "IVAR", 
  "IVARARRAY", "BVAR", "BVARARRAY", "DTVAR", "EQUALS", "NEQUALS", 
  "UNDEFINED_SYMBOL", "INTERNAL_ERROR", "'?'", "OR", "AND", "'<'", "'>'", 
  "'-'", "'+'", "'.'", "'*'", "'/'", "NEG", "'^'", "'!'", "'['", "']'", 
  "'('", "')'", "':'", "','", "'{'", "'}'", "$accept", "line", "bExp", 
  "iExp", "dExp", "dtExp", "schedExp", "tabFuncExp", "dVar", "iVar", 
  "bVar", "dtVar", "dVarArray", "iVarArray", "bVarArray", "expArrayBegin", 
  "dValArrayMid", "dValArray", "iValArrayMid", "iValArray", 
  "bValArrayMid", "bValArray", "dExpArray", "iExpArray", "bExpArray", 
  "bFuncArgsMid", "bFuncArgs", "iFuncArgsMid", "iFuncArgs", 
  "dFuncArgsMid", "dFuncArgs", "dArrayFuncArgsMid", "dArrayFuncArgs", 
  "iArrayFuncArgsMid", "iArrayFuncArgs", "bArrayFuncArgsMid", 
  "bArrayFuncArgs", "dtFuncArgsMid", "dtFuncArgs", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,    63,   288,   289,    60,    62,    45,    43,
      46,    42,    47,   290,    94,    33,    91,    93,    40,    41,
      58,    44,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    54,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    57,    57,    57,    57,    57,    57,
      57,    57,    57,    57,    57,    57,    57,    57,    57,    57,
      57,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    70,    70,    71,    72,    72,    73,
      74,    74,    75,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    77,    77,
      77,    77,    77,    77,    77,    77,    77,    77,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    80,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    80,    81,    81,    81,    81,    81,    81,    81,    81,
      81,    81,    81,    81,    81,    81,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    83,
      83,    83,    83,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    86,    86,    86,    86,    86,    86,    86,    86,    86,
      86,    86,    86,    86,    86,    86,    86,    87,    87,    87,
      87,    87,    87,    87,    87,    87,    87,    87,    87,    87,
      87,    87,    87,    87,    88,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    88,    88,    88,    88,    88,
      89,    89,    89,    89,    89,    89,    89,    89,    89,    89,
      89,    89,    89,    89,    89,    89,    89,    90,    90,    90,
      90,    90,    90,    90,    90,    90,    90,    90,    90,    90,
      90,    90,    90,    91,    91,    91,    91,    91,    91,    91,
      91,    91,    91,    91,    91,    91,    91,    92,    92,    92,
      92,    92,    92,    92,    92,    92,    92,    92,    92,    92
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     4,     4,     2,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     5,     1,     1,     4,     4,     3,     3,     3,
       3,     2,     3,     3,     5,     6,     6,     1,     3,     4,
       4,     1,     4,     4,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       3,     6,     6,     6,     6,     6,     6,     5,     5,     5,
       1,     3,     4,     4,     1,     4,     3,     3,     3,     6,
       6,     5,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     2,     3,     3,     2,     2,     3,     2,
       2,     3,     2,     1,     1,     4,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       2,     1,     5,     5,     5,     5,     8,     3,     1,     1,
       4,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       2,     1,     5,     5,     5,     5,     8,     3,     1,     1,
       2,     4,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     5,
       5,     5,     5,     8,     3,     2,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     5,     5,     5,     5,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     5,     5,
       5,     5,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     5,     5,     5,     5,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     5,     5,     5,     5,     2,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     5,
       5,     5,     5,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     5,     5,     5,     5,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     5,     5,     5,     5,
       5,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       5,     5,     5,     5,     5,     5,     5,     2,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     5,     5,     5,
       5,     5,     5,     5,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     5,     5,     5,     5,     5,     5,     5,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       5,     5,     5,     5,     5,     5,     5,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     5,     5,     5,     5,
       5,     5,     5,     2,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     5,     5,     5,     5,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     5,     5,     5,     5
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned short yydefact[] =
{
       0,    51,   114,    84,    11,   159,    34,   139,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    95,    99,
      96,   100,    97,   101,    98,     9,    10,     0,     0,     0,
     102,     0,     2,     3,     4,     5,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   113,     0,   138,     0,   158,
       6,     7,     8,     0,    33,     0,    47,     0,    80,     0,
     131,     0,   151,     0,   162,     0,    92,     0,     0,   212,
     297,   185,   330,   239,   266,   363,     0,     0,    41,    69,
       0,   130,   150,     0,    14,     0,     0,     0,     0,   160,
       0,     0,     0,     0,     0,     0,     0,     1,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   110,   107,   103,     0,   106,     0,
     109,     0,   112,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    93,
      94,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    17,    42,    70,    88,   137,   157,
     184,    15,    16,     0,     0,     0,     0,     0,     0,     0,
      31,    30,    18,    19,    24,    25,    38,    59,    37,    56,
      39,    62,   125,   146,    40,    65,    68,    20,    21,    26,
      27,    58,    57,    55,    54,    61,    60,   124,    64,    63,
      67,    66,    22,    23,    28,    29,    87,    43,    86,     0,
       0,     0,     0,     0,     0,     0,   105,   104,   108,   111,
     177,   178,   166,   164,   165,   163,   120,   119,   121,   117,
     116,   118,    81,   123,   122,   126,   128,   127,   129,     0,
     175,   176,   170,   168,   169,   167,   143,   144,   141,   142,
      48,   145,   147,   148,   149,     0,   173,   174,     0,     0,
       0,     0,     0,     0,   171,   172,     0,   199,   186,   201,
     188,   203,   190,   205,   192,   206,   193,   207,   194,     0,
       0,     0,     0,   204,   191,   202,   189,   200,   187,   226,
     213,   228,   215,   230,   217,   232,   219,   233,   220,   234,
     221,     0,     0,     0,     0,   231,   218,   229,   216,   227,
     214,   253,   240,   255,   242,   257,   244,   259,   246,   260,
     247,   261,   248,     0,     0,     0,     0,   258,   245,   256,
     243,   254,   241,   281,   267,   283,   269,   285,   271,   287,
     273,   288,   274,   289,   275,     0,     0,     0,     0,     0,
       0,     0,   286,   272,   284,   270,   282,   268,   314,   298,
     316,   300,   318,   302,   320,   304,   321,   305,   322,   306,
       0,     0,     0,     0,     0,     0,     0,   319,   303,   317,
     301,   315,   299,   347,   331,   349,   333,   351,   335,   353,
     337,   354,   338,   355,   339,     0,     0,     0,     0,     0,
       0,     0,   352,   336,   350,   334,   348,   332,   377,   364,
     379,   366,   381,   368,   383,   370,   384,   371,   385,   372,
       0,     0,     0,     0,   382,   369,   380,   367,   378,   365,
       0,     0,     0,    82,    49,     0,     0,     0,    83,    50,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      52,    35,    12,    85,   115,   140,   161,    53,    36,     0,
       0,     0,     0,     0,     0,    13,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    32,    44,    78,    79,
      77,    91,   135,   155,   182,   181,   154,   134,   133,   132,
     153,   152,   180,   179,   210,   197,   209,   196,   208,   195,
     211,   198,   237,   224,   236,   223,   235,   222,   238,   225,
     264,   251,   263,   250,   262,   249,   265,   252,   294,   278,
     292,   277,   290,   276,   296,   280,   295,   279,   293,   291,
     327,   311,   325,   309,   323,   307,   329,   313,   328,   312,
     326,   310,   324,   308,   360,   344,   358,   342,   356,   340,
     362,   346,   361,   345,   359,   343,   357,   341,   388,   375,
     387,   374,   386,   373,   389,   376,    45,    72,    73,    71,
      89,    46,    75,    76,    74,    90,     0,     0,     0,     0,
       0,     0,   136,   156,   183
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    31,    77,    85,    86,    80,   175,   176,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    87,    88,    83,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -47
static const short yypact[] =
{
     779,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -42,   -29,
     -10,    -6,    54,    64,    94,    95,   101,   105,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,  2777,  2777,  2777,
     -47,    15,   111,  4021,  4037,   453,     1,    44,    55,    65,
     108,   113,   125,  2777,   -46,   -47,    13,   -47,    21,   -47,
    3993,  4007,   384,   845,   -47,   845,   -47,   845,   -47,   845,
     -47,   845,   -47,   845,   -47,   845,   -47,  2777,  2777,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,  2777,   111,   128,   132,
     453,   162,   166,   384,   -47,  4021,  4037,  3993,  4007,   177,
     140,  3582,  3598,  3790,   732,  3512,  3801,   -47,  2777,  2777,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
    2777,  2777,  2777,  2777,   111,  4021,  4037,  2777,   -47,  2777,
     -47,  2777,   -47,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,   -47,
     -47,   258,  2812,  2828,    85,    36,   131,   210,   213,   226,
     239,   -28,    -5,   303,   550,  2844,  2860,   458,   153,   179,
     243,   244,   260,   285,    40,   180,   641,  3286,  2876,  2892,
     485,   222,   229,   289,   307,   314,   315,   228,   347,  3169,
    3293,  2908,  2924,  1049,   232,   245,   336,   344,   345,   348,
     353,   358,   361,   399,   422,  3176,  3316,  2940,  2956,  3194,
     248,   299,   362,   369,   370,   377,   378,   397,   411,   576,
     600,  3205,  3323,  2972,  2988,  3223,   329,   330,   420,   423,
     428,   439,   447,   452,   466,   665,   714,  3234,  3346,  3004,
    3020,  3252,   346,   395,   474,   491,   492,   501,  1065,  2794,
    3263,  3054,  3070,  3332,  3526,  3540,  3086,  3102,  3348,  3554,
    3568,   503,  3036,  3359,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,    63,  3422,  3438,  3486,   796,  3377,   846,
     -13,    22,  1081,  1081,  1081,  1081,   183,   219,   183,   219,
     128,   132,   162,   166,   128,   132,   132,  1170,  1170,  1170,
    1170,   183,   219,   183,   219,   128,   132,   162,   128,   132,
     128,   132,   320,   320,   320,   320,   183,   -47,   183,  3833,
    3849,  3865,  3881,  3897,  3913,  3929,  4021,  4037,  4021,   111,
     288,   288,  1170,   288,  1170,   288,   183,   219,     8,   183,
     219,     8,    42,   128,   132,   162,   128,   132,   162,  3945,
     699,   699,  1081,   699,  1081,   699,   183,   165,   183,   165,
      53,   128,   166,   128,   166,  3961,   177,   177,   392,  3454,
    3470,  3392,  3407,  3497,    10,   -26,  3977,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   895,
     945,   995,  1045,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,  1127,  1177,  1227,  1277,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,  1327,  1377,  1427,  1477,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,  1527,  1577,  1627,  1677,  1727,
    1777,  1827,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
    1877,  1927,  1977,  2027,  2077,  2127,  2177,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,  2227,  2277,  2327,  2377,  2427,
    2477,  2527,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
    2577,  2627,  2677,  2727,   -47,   -47,   -47,   -47,   -47,   -47,
    2777,  2777,  2777,   -47,   -47,  2777,  2777,  2777,   -47,   -47,
    2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,  2777,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,  2777,
    2777,  2777,  2777,  2777,  2777,   -47,   398,   416,   421,   435,
     459,   462,   467,   468,   497,   506,   507,   518,   527,   538,
     558,   559,   570,   450,   510,   582,   583,   594,   606,   610,
     617,   618,   640,   647,   661,   664,   674,   675,   679,   682,
     708,   726,   763,  3614,  3630,  3646,  3662,  3815,  3678,  3694,
    3710,  3726,  3819,  3118,  3134,  3150,   111,  4021,  4037,  4021,
    4037,   453,  3993,  4007,   384,   384,  4007,  3993,  4037,  3993,
    4021,  4007,   111,   384,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,  2777,  2777,  2777,  3742,
    3758,  3774,   -47,   -47,   -47
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
     -47,   -47,   135,     0,   255,   597,   122,   256,   187,   871,
     921,   971,   -33,    -1,    12,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   535,   186,   138,   -47,   -47,   -47,   -47,   -47,
     -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47,   -47
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned short yytable[] =
{
      33,   143,   144,   163,   164,   137,    67,   138,   145,   146,
     147,   148,   149,   150,   151,    97,    98,    99,   152,    68,
     168,   413,   102,   414,   153,   154,   220,    78,   236,    91,
     252,   155,   156,   157,   158,   159,   160,   161,    69,   163,
     164,   162,    70,   135,   415,   167,   416,   127,   149,   150,
     151,    98,    99,   172,   152,   185,   168,   198,   221,   211,
     237,   227,   253,   243,   139,   259,   140,   271,   276,   143,
     144,   222,   141,   238,   142,   254,   145,   146,   147,   148,
     149,   150,   151,   150,   151,   405,   152,   406,   152,   435,
     128,   436,    98,    99,   160,   161,   100,   101,   102,   162,
     294,   129,    71,   302,   303,   304,   305,   306,   308,   310,
     314,   130,    72,   573,   121,   122,   321,   323,   325,   328,
     330,   123,   124,   125,   126,   336,   338,   339,   340,   341,
     342,   343,   344,   345,   403,    32,   404,   346,    52,   348,
      98,    99,    73,    74,   100,   101,   102,   356,   359,    75,
     363,   366,   369,    76,   131,   372,   374,   376,   378,   132,
     381,   383,   385,    84,    90,   389,    89,    96,   396,    98,
      99,   133,   111,   100,   101,   102,   120,   188,   134,   201,
     407,   214,   408,   230,     0,   246,    51,   262,   171,   284,
     184,   183,   197,   196,   210,   209,   226,   225,   242,   241,
     258,   257,   427,   270,   428,   159,   160,   161,   152,   153,
     154,   162,   162,    82,   283,    95,   155,   156,   157,   158,
     159,   160,   161,   168,   109,   110,   162,   111,   429,   437,
     430,   438,     0,   291,   292,   293,   300,   301,   299,   182,
     177,   195,   190,   208,   203,   224,   216,   240,   232,   256,
     248,   269,   264,   275,   280,    34,   409,   143,   144,   410,
     118,   119,   282,   120,   145,   146,   147,   148,   149,   150,
     151,   449,   411,   450,   152,     0,   349,   457,   451,   458,
     452,   471,    79,   472,    92,   412,   298,    98,    99,   431,
     432,   100,   101,   102,   473,   313,   474,   496,   136,   497,
     388,   386,   387,   393,   394,   395,   433,   397,   173,   398,
     186,   189,   199,   202,   212,   215,   228,   231,   244,   247,
     260,   263,   272,   277,   145,   146,   147,   148,   149,   150,
     151,   434,   163,   164,   152,   453,   165,   166,   167,   370,
     371,   373,   375,   377,   379,   380,   382,   384,   498,   168,
     499,   392,   417,   454,   418,   295,   123,   124,   125,   126,
     455,   456,   307,   309,   311,   315,   316,   317,   318,   319,
     320,   322,   324,   326,   329,   331,   153,   154,   521,   523,
     522,   524,   475,   155,   156,   157,   158,   159,   160,   161,
     476,   477,   347,   162,   478,   546,   459,   547,   460,   479,
     352,   354,   357,   360,   480,   364,   367,   481,   500,   339,
     340,   341,   342,   163,   164,   501,   502,   165,   166,   167,
     390,    98,    99,   503,   504,   100,   101,   102,   143,   144,
     168,   339,   340,   341,   342,   145,   146,   147,   148,   149,
     150,   151,   589,   505,   548,   152,   549,   664,   482,   665,
     483,   153,   154,   339,   340,   341,   342,   506,   155,   156,
     157,   158,   159,   160,   161,   666,   525,   667,   162,   526,
     668,   484,   669,   485,   527,   339,   340,   341,   342,   343,
     344,   345,   121,   122,   670,   528,   671,   121,   122,   123,
     124,   125,   126,   529,   123,   124,   125,   126,   530,   698,
     339,   340,   341,   342,   343,   344,   345,   425,   672,   426,
     673,   674,   531,   675,   121,   122,   676,   678,   677,   679,
     550,   123,   124,   125,   126,   339,   340,   341,   342,   343,
     344,   345,   143,   144,   447,    50,   448,   551,   552,   145,
     146,   147,   148,   149,   150,   151,   680,   553,   681,   152,
     339,   340,   341,   342,   570,   682,   684,   683,   685,   699,
     633,   635,    81,     0,    94,   638,   640,   686,     0,   687,
     643,   644,   645,     0,   647,   649,   688,     0,   689,    98,
      99,     0,     0,   100,   101,   102,     0,   690,   181,   691,
     194,     0,   207,   660,   223,     0,   239,    35,   255,   419,
     268,   420,   274,   279,     0,   143,   144,   692,   694,   693,
     695,   281,   145,   146,   147,   148,   149,   150,   151,   696,
       0,   697,   152,     0,     0,   507,    93,   508,     0,   153,
     154,   700,   702,   701,   703,   297,   155,   156,   157,   158,
     159,   160,   161,   704,   312,   705,   162,     0,     0,   509,
     174,   510,   187,   327,   200,   706,   213,   707,   229,   708,
     245,   709,   261,     0,   273,   278,   710,   712,   711,   713,
     163,   164,     0,     0,   165,   166,   167,     0,   350,   351,
     353,   355,   358,   361,   362,   365,   368,   168,     0,   714,
     439,   715,   440,     0,   143,   144,   716,   296,   717,     0,
     391,   145,   146,   147,   148,   149,   150,   151,   646,     0,
     718,   152,   719,   720,   532,   721,   533,   654,   332,   333,
     334,   335,   337,   722,   724,   723,   725,   655,   726,   662,
     727,   728,   663,   729,     0,   155,   156,   157,   158,   159,
     160,   161,     0,   153,   154,   162,   749,   750,   751,     0,
     155,   156,   157,   158,   159,   160,   161,   730,     0,   731,
     162,   143,   144,   534,   653,   535,     0,     0,   145,   146,
     147,   148,   149,   150,   151,   732,   656,   733,   152,   661,
       0,   288,     1,     2,     3,     4,     5,     6,     7,     0,
       0,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,     0,     0,
      25,    26,   734,     0,   735,   634,   636,    27,     0,     0,
     639,   641,     0,     0,    28,   143,   144,    29,     0,   648,
     650,    30,   145,   146,   147,   148,   149,   150,   151,     0,
       0,     0,   152,     0,     0,     0,   577,   658,     1,     2,
       3,     4,     5,     6,     7,   169,   170,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,     0,   163,   164,     0,     0,   165,
     166,   167,     0,    27,     0,     0,     0,     0,     0,     0,
      28,     0,   168,    29,     0,     0,   579,    30,     1,     2,
       3,     4,     5,     6,     7,     0,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,   178,     0,   191,     0,   204,     0,
     217,     0,   233,    27,   249,     0,   265,     0,     0,     0,
      28,     0,   596,    29,     0,     0,     0,    30,     1,     2,
       3,     4,     5,     6,     7,     0,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,   179,     0,   192,     0,   205,     0,
     218,     0,   234,    27,   250,     0,   266,     0,     0,     0,
      28,     0,   597,    29,     0,     0,     0,    30,     1,     2,
       3,     4,     5,     6,     7,     0,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,   180,     0,   193,     0,   206,     0,
     219,     0,   235,    27,   251,     0,   267,     0,     0,     0,
      28,     0,   598,    29,     0,     0,     0,    30,     1,     2,
       3,     4,     5,     6,     7,     0,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,     0,     0,     0,     0,   121,   122,
       0,     0,     0,    27,     0,   123,   124,   125,   126,     0,
      28,     0,   599,    29,   143,   144,     0,    30,   469,     0,
     470,   145,   146,   147,   148,   149,   150,   151,     0,     0,
       0,   152,   652,     0,   554,     0,   555,   105,   106,   107,
     108,     0,   109,   110,     0,   111,   657,   659,     0,     0,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,   637,
       0,     0,     0,     0,   642,    27,     0,     0,     0,     0,
       0,     0,    28,   651,   600,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,   114,   115,   116,   117,
       0,   118,   119,     0,   120,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   601,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   602,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   603,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   604,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   605,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   606,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   607,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   608,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   609,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   610,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   611,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   612,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   613,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   614,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   615,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   616,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   617,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   618,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   619,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   620,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   621,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   622,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   623,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   624,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   625,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   626,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   627,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   628,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   629,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   630,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   631,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,     0,   632,    29,     0,     0,     0,    30,
       1,     2,     3,     4,     5,     6,     7,     0,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    27,     0,     0,     0,     0,
       0,     0,    28,   153,   154,    29,     0,     0,     0,    30,
     155,   156,   157,   158,   159,   160,   161,     0,     0,     0,
     162,   103,   104,   556,     0,   557,     0,     0,   105,   106,
     107,   108,     0,   109,   110,     0,   111,   112,   113,     0,
       0,   399,     0,   400,   114,   115,   116,   117,     0,   118,
     119,     0,   120,   103,   104,     0,     0,   401,     0,   402,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   112,
     113,     0,     0,   421,     0,   422,   114,   115,   116,   117,
       0,   118,   119,     0,   120,   103,   104,     0,     0,   423,
       0,   424,   105,   106,   107,   108,     0,   109,   110,     0,
     111,   112,   113,     0,     0,   443,     0,   444,   114,   115,
     116,   117,     0,   118,   119,     0,   120,   103,   104,     0,
       0,   445,     0,   446,   105,   106,   107,   108,     0,   109,
     110,     0,   111,   112,   113,     0,     0,   465,     0,   466,
     114,   115,   116,   117,     0,   118,   119,     0,   120,   103,
     104,     0,     0,   467,     0,   468,   105,   106,   107,   108,
       0,   109,   110,     0,   111,   112,   113,     0,     0,   490,
       0,   491,   114,   115,   116,   117,     0,   118,   119,     0,
     120,   103,   104,     0,     0,   492,     0,   493,   105,   106,
     107,   108,     0,   109,   110,     0,   111,   112,   113,     0,
       0,   515,     0,   516,   114,   115,   116,   117,     0,   118,
     119,     0,   120,   103,   104,     0,     0,   517,     0,   518,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   112,
     113,     0,     0,   540,     0,   541,   114,   115,   116,   117,
       0,   118,   119,     0,   120,   153,   154,     0,     0,   542,
       0,   543,   155,   156,   157,   158,   159,   160,   161,     0,
       0,     0,   162,   103,   104,     0,     0,   571,     0,     0,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   112,
     113,     0,     0,     0,     0,   560,   114,   115,   116,   117,
       0,   118,   119,     0,   120,   103,   104,     0,     0,     0,
       0,   561,   105,   106,   107,   108,     0,   109,   110,     0,
     111,   112,   113,     0,     0,     0,     0,   565,   114,   115,
     116,   117,     0,   118,   119,     0,   120,   103,   104,     0,
       0,     0,     0,   566,   105,   106,   107,   108,     0,   109,
     110,     0,   111,   103,   104,     0,     0,     0,     0,   746,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   103,
     104,     0,     0,     0,     0,   747,   105,   106,   107,   108,
       0,   109,   110,     0,   111,     0,     0,     0,   163,   164,
       0,   748,   165,   166,   167,   163,   164,     0,     0,   165,
     166,   167,     0,     0,     0,   168,     0,     0,   461,     0,
     462,     0,   168,   121,   122,   486,     0,   487,     0,     0,
     123,   124,   125,   126,   163,   164,     0,     0,   165,   166,
     167,     0,     0,   494,     0,   495,     0,     0,     0,     0,
       0,   168,   121,   122,   511,     0,   512,     0,     0,   123,
     124,   125,   126,   163,   164,     0,     0,   165,   166,   167,
       0,     0,   519,     0,   520,     0,     0,     0,     0,     0,
     168,   121,   122,   536,     0,   537,     0,     0,   123,   124,
     125,   126,   163,   164,     0,     0,   165,   166,   167,     0,
       0,   544,     0,   545,     0,     0,     0,     0,     0,   168,
       0,     0,   558,     0,   559,    98,    99,     0,     0,   100,
     101,   102,    98,    99,     0,     0,   100,   101,   102,     0,
       0,     0,     0,     0,     0,   441,     0,   442,     0,     0,
       0,     0,   463,     0,   464,    98,    99,     0,     0,   100,
     101,   102,    98,    99,     0,     0,   100,   101,   102,     0,
       0,   121,   122,     0,     0,   488,     0,   489,   123,   124,
     125,   126,   513,     0,   514,    98,    99,   121,   122,   100,
     101,   102,     0,   562,   123,   124,   125,   126,   163,   164,
       0,     0,   165,   166,   167,   538,     0,   539,     0,   567,
       0,     0,     0,     0,     0,   168,   153,   154,     0,     0,
     572,     0,     0,   155,   156,   157,   158,   159,   160,   161,
       0,   143,   144,   162,     0,     0,     0,   578,   145,   146,
     147,   148,   149,   150,   151,     0,   153,   154,   152,     0,
       0,     0,   592,   155,   156,   157,   158,   159,   160,   161,
       0,   103,   104,   162,     0,     0,     0,   593,   105,   106,
     107,   108,     0,   109,   110,     0,   111,   112,   113,     0,
       0,     0,   574,     0,   114,   115,   116,   117,     0,   118,
     119,     0,   120,   103,   104,     0,     0,     0,   575,     0,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   112,
     113,     0,     0,     0,   590,     0,   114,   115,   116,   117,
       0,   118,   119,     0,   120,   121,   122,     0,     0,     0,
     591,     0,   123,   124,   125,   126,   163,   164,     0,     0,
     165,   166,   167,     0,     0,     0,   576,     0,     0,     0,
       0,   153,   154,   168,     0,     0,     0,   594,   155,   156,
     157,   158,   159,   160,   161,   143,   144,     0,   162,     0,
       0,   289,   145,   146,   147,   148,   149,   150,   151,   153,
     154,     0,   152,     0,     0,   563,   155,   156,   157,   158,
     159,   160,   161,   143,   144,     0,   162,     0,     0,   564,
     145,   146,   147,   148,   149,   150,   151,   153,   154,     0,
     152,     0,     0,   568,   155,   156,   157,   158,   159,   160,
     161,   103,   104,     0,   162,     0,     0,   569,   105,   106,
     107,   108,     0,   109,   110,     0,   111,   112,   113,     0,
       0,   285,     0,     0,   114,   115,   116,   117,     0,   118,
     119,     0,   120,   103,   104,     0,     0,   286,     0,     0,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   112,
     113,     0,     0,   736,     0,     0,   114,   115,   116,   117,
       0,   118,   119,     0,   120,   103,   104,     0,     0,   737,
       0,     0,   105,   106,   107,   108,     0,   109,   110,     0,
     111,   112,   113,     0,     0,   738,     0,     0,   114,   115,
     116,   117,     0,   118,   119,     0,   120,   103,   104,     0,
       0,   739,     0,     0,   105,   106,   107,   108,     0,   109,
     110,     0,   111,   112,   113,     0,     0,   741,     0,     0,
     114,   115,   116,   117,     0,   118,   119,     0,   120,   103,
     104,     0,     0,   742,     0,     0,   105,   106,   107,   108,
       0,   109,   110,     0,   111,   112,   113,     0,     0,   743,
       0,     0,   114,   115,   116,   117,     0,   118,   119,     0,
     120,   103,   104,     0,     0,   744,     0,     0,   105,   106,
     107,   108,     0,   109,   110,     0,   111,   103,   104,     0,
       0,   752,     0,     0,   105,   106,   107,   108,     0,   109,
     110,     0,   111,   103,   104,     0,     0,   753,     0,     0,
     105,   106,   107,   108,     0,   109,   110,     0,   111,   121,
     122,     0,     0,   754,     0,     0,   123,   124,   125,   126,
     163,   164,     0,     0,   165,   166,   167,     0,     0,   287,
       0,     0,     0,     0,   121,   122,     0,   168,   121,   122,
     290,   123,   124,   125,   126,   123,   124,   125,   126,     0,
       0,     0,   103,   104,   740,     0,     0,     0,   745,   105,
     106,   107,   108,     0,   109,   110,     0,   111,   103,   104,
     580,     0,     0,     0,     0,   105,   106,   107,   108,     0,
     109,   110,     0,   111,   103,   104,   581,     0,     0,     0,
       0,   105,   106,   107,   108,     0,   109,   110,     0,   111,
     103,   104,   582,     0,     0,     0,     0,   105,   106,   107,
     108,     0,   109,   110,     0,   111,   103,   104,   583,     0,
       0,     0,     0,   105,   106,   107,   108,     0,   109,   110,
       0,   111,   103,   104,   584,     0,     0,     0,     0,   105,
     106,   107,   108,     0,   109,   110,     0,   111,   103,   104,
     585,     0,     0,     0,     0,   105,   106,   107,   108,     0,
     109,   110,     0,   111,   103,   104,   586,     0,     0,     0,
       0,   105,   106,   107,   108,     0,   109,   110,     0,   111,
     103,   104,   587,     0,     0,     0,     0,   105,   106,   107,
     108,     0,   109,   110,     0,   111,   103,   104,   588,     0,
       0,     0,     0,   105,   106,   107,   108,     0,   109,   110,
       0,   111,   143,   144,   595,     0,     0,     0,     0,   145,
     146,   147,   148,   149,   150,   151,   153,   154,     0,   152,
       0,     0,     0,   155,   156,   157,   158,   159,   160,   161,
     103,   104,     0,   162,     0,     0,     0,   105,   106,   107,
     108,     0,   109,   110,     0,   111,   112,   113,     0,     0,
       0,     0,     0,   114,   115,   116,   117,     0,   118,   119,
       0,   120
};

static const short yycheck[] =
{
       0,    29,    30,    29,    30,    51,    48,    53,    36,    37,
      38,    39,    40,    41,    42,     0,    29,    30,    46,    48,
      46,    49,    35,    51,    29,    30,    59,    27,    61,    29,
      63,    36,    37,    38,    39,    40,    41,    42,    48,    29,
      30,    46,    48,    43,    49,    35,    51,    46,    40,    41,
      42,    29,    30,    53,    46,    55,    46,    57,    59,    59,
      61,    61,    63,    63,    51,    65,    53,    67,    68,    29,
      30,    59,    51,    61,    53,    63,    36,    37,    38,    39,
      40,    41,    42,    41,    42,    49,    46,    51,    46,    49,
      46,    51,    29,    30,    41,    42,    33,    34,    35,    46,
     100,    46,    48,   103,   104,   105,   106,   107,   108,   109,
     110,    46,    48,    50,    29,    30,   116,   117,   118,   119,
     120,    36,    37,    38,    39,   125,   126,   127,   128,   129,
     130,   131,   132,   133,    49,     0,    51,   137,     0,   139,
      29,    30,    48,    48,    33,    34,    35,   147,   148,    48,
     150,   151,   152,    48,    46,   155,   156,   157,   158,    46,
     160,   161,   162,    28,    29,   165,    28,    29,   168,    29,
      30,    46,    44,    33,    34,    35,    44,    55,    43,    57,
      49,    59,    51,    61,    -1,    63,     0,    65,    53,    49,
      55,    53,    57,    55,    59,    57,    61,    59,    63,    61,
      65,    63,    49,    65,    51,    40,    41,    42,    46,    29,
      30,    46,    46,    27,    76,    29,    36,    37,    38,    39,
      40,    41,    42,    46,    41,    42,    46,    44,    49,    49,
      51,    51,    -1,    98,    99,   100,   101,   102,   100,    53,
      53,    55,    55,    57,    57,    59,    59,    61,    61,    63,
      63,    65,    65,    67,    68,     0,    46,    29,    30,    46,
      41,    42,    76,    44,    36,    37,    38,    39,    40,    41,
      42,    49,    46,    51,    46,    -1,   141,    49,    49,    51,
      51,    49,    27,    51,    29,    46,   100,    29,    30,    46,
      46,    33,    34,    35,    49,   109,    51,    49,    43,    51,
     165,   163,   164,   165,   166,   167,    46,    49,    53,    51,
      55,    55,    57,    57,    59,    59,    61,    61,    63,    63,
      65,    65,    67,    68,    36,    37,    38,    39,    40,    41,
      42,    46,    29,    30,    46,    46,    33,    34,    35,   153,
     154,   155,   156,   157,   158,   159,   160,   161,    49,    46,
      51,   165,    49,    46,    51,   100,    36,    37,    38,    39,
      46,    46,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,    29,    30,    49,    49,
      51,    51,    46,    36,    37,    38,    39,    40,    41,    42,
      46,    46,   137,    46,    46,    49,    49,    51,    51,    46,
     145,   146,   147,   148,    46,   150,   151,    46,    46,   409,
     410,   411,   412,    29,    30,    46,    46,    33,    34,    35,
     165,    29,    30,    46,    46,    33,    34,    35,    29,    30,
      46,   431,   432,   433,   434,    36,    37,    38,    39,    40,
      41,    42,    50,    46,    49,    46,    51,    49,    49,    51,
      51,    29,    30,   453,   454,   455,   456,    46,    36,    37,
      38,    39,    40,    41,    42,    49,    46,    51,    46,    46,
      49,    49,    51,    51,    46,   475,   476,   477,   478,   479,
     480,   481,    29,    30,    49,    46,    51,    29,    30,    36,
      37,    38,    39,    46,    36,    37,    38,    39,    46,    49,
     500,   501,   502,   503,   504,   505,   506,    49,    49,    51,
      51,    49,    46,    51,    29,    30,    49,    49,    51,    51,
      46,    36,    37,    38,    39,   525,   526,   527,   528,   529,
     530,   531,    29,    30,    49,     0,    51,    46,    46,    36,
      37,    38,    39,    40,    41,    42,    49,    46,    51,    46,
     550,   551,   552,   553,    51,    49,    49,    51,    51,    49,
     560,   561,    27,    -1,    29,   565,   566,    49,    -1,    51,
     570,   571,   572,    -1,   574,   575,    49,    -1,    51,    29,
      30,    -1,    -1,    33,    34,    35,    -1,    49,    53,    51,
      55,    -1,    57,   593,    59,    -1,    61,     0,    63,    49,
      65,    51,    67,    68,    -1,    29,    30,    49,    49,    51,
      51,    76,    36,    37,    38,    39,    40,    41,    42,    49,
      -1,    51,    46,    -1,    -1,    49,    29,    51,    -1,    29,
      30,    49,    49,    51,    51,   100,    36,    37,    38,    39,
      40,    41,    42,    49,   109,    51,    46,    -1,    -1,    49,
      53,    51,    55,   118,    57,    49,    59,    51,    61,    49,
      63,    51,    65,    -1,    67,    68,    49,    49,    51,    51,
      29,    30,    -1,    -1,    33,    34,    35,    -1,   143,   144,
     145,   146,   147,   148,   149,   150,   151,    46,    -1,    49,
      49,    51,    51,    -1,    29,    30,    49,   100,    51,    -1,
     165,    36,    37,    38,    39,    40,    41,    42,   573,    -1,
      49,    46,    51,    49,    49,    51,    51,   579,   121,   122,
     123,   124,   125,    49,    49,    51,    51,   589,    49,   594,
      51,    49,   594,    51,    -1,    36,    37,    38,    39,    40,
      41,    42,    -1,    29,    30,    46,   746,   747,   748,    -1,
      36,    37,    38,    39,    40,    41,    42,    49,    -1,    51,
      46,    29,    30,    49,   578,    51,    -1,    -1,    36,    37,
      38,    39,    40,    41,    42,    49,   590,    51,    46,   593,
      -1,    49,     3,     4,     5,     6,     7,     8,     9,    -1,
      -1,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    -1,    -1,
      31,    32,    49,    -1,    51,   560,   561,    38,    -1,    -1,
     565,   566,    -1,    -1,    45,    29,    30,    48,    -1,   574,
     575,    52,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    -1,    46,    -1,    -1,    -1,    50,   592,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    -1,    29,    30,    -1,    -1,    33,
      34,    35,    -1,    38,    -1,    -1,    -1,    -1,    -1,    -1,
      45,    -1,    46,    48,    -1,    -1,    50,    52,     3,     4,
       5,     6,     7,     8,     9,    -1,    -1,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    53,    -1,    55,    -1,    57,    -1,
      59,    -1,    61,    38,    63,    -1,    65,    -1,    -1,    -1,
      45,    -1,    47,    48,    -1,    -1,    -1,    52,     3,     4,
       5,     6,     7,     8,     9,    -1,    -1,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    53,    -1,    55,    -1,    57,    -1,
      59,    -1,    61,    38,    63,    -1,    65,    -1,    -1,    -1,
      45,    -1,    47,    48,    -1,    -1,    -1,    52,     3,     4,
       5,     6,     7,     8,     9,    -1,    -1,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    53,    -1,    55,    -1,    57,    -1,
      59,    -1,    61,    38,    63,    -1,    65,    -1,    -1,    -1,
      45,    -1,    47,    48,    -1,    -1,    -1,    52,     3,     4,
       5,     6,     7,     8,     9,    -1,    -1,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    -1,    -1,    -1,    -1,    29,    30,
      -1,    -1,    -1,    38,    -1,    36,    37,    38,    39,    -1,
      45,    -1,    47,    48,    29,    30,    -1,    52,    49,    -1,
      51,    36,    37,    38,    39,    40,    41,    42,    -1,    -1,
      -1,    46,   577,    -1,    49,    -1,    51,    36,    37,    38,
      39,    -1,    41,    42,    -1,    44,   591,   592,    -1,    -1,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,   562,
      -1,    -1,    -1,    -1,   567,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,   576,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    47,    48,    -1,    -1,    -1,    52,
       3,     4,     5,     6,     7,     8,     9,    -1,    -1,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    29,    30,    48,    -1,    -1,    -1,    52,
      36,    37,    38,    39,    40,    41,    42,    -1,    -1,    -1,
      46,    29,    30,    49,    -1,    51,    -1,    -1,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    51,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    51,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    51,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    49,
      -1,    51,    36,    37,    38,    39,    -1,    41,    42,    -1,
      44,    29,    30,    -1,    -1,    49,    -1,    51,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    51,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    51,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    51,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    49,
      -1,    51,    36,    37,    38,    39,    -1,    41,    42,    -1,
      44,    29,    30,    -1,    -1,    49,    -1,    51,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    51,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    51,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    51,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    49,
      -1,    51,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    -1,    46,    29,    30,    -1,    -1,    51,    -1,    -1,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    -1,    -1,    51,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    -1,
      -1,    51,    36,    37,    38,    39,    -1,    41,    42,    -1,
      44,    29,    30,    -1,    -1,    -1,    -1,    51,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    -1,    -1,    51,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    -1,    -1,    51,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    -1,    -1,    51,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    -1,    -1,    -1,    29,    30,
      -1,    51,    33,    34,    35,    29,    30,    -1,    -1,    33,
      34,    35,    -1,    -1,    -1,    46,    -1,    -1,    49,    -1,
      51,    -1,    46,    29,    30,    49,    -1,    51,    -1,    -1,
      36,    37,    38,    39,    29,    30,    -1,    -1,    33,    34,
      35,    -1,    -1,    49,    -1,    51,    -1,    -1,    -1,    -1,
      -1,    46,    29,    30,    49,    -1,    51,    -1,    -1,    36,
      37,    38,    39,    29,    30,    -1,    -1,    33,    34,    35,
      -1,    -1,    49,    -1,    51,    -1,    -1,    -1,    -1,    -1,
      46,    29,    30,    49,    -1,    51,    -1,    -1,    36,    37,
      38,    39,    29,    30,    -1,    -1,    33,    34,    35,    -1,
      -1,    49,    -1,    51,    -1,    -1,    -1,    -1,    -1,    46,
      -1,    -1,    49,    -1,    51,    29,    30,    -1,    -1,    33,
      34,    35,    29,    30,    -1,    -1,    33,    34,    35,    -1,
      -1,    -1,    -1,    -1,    -1,    49,    -1,    51,    -1,    -1,
      -1,    -1,    49,    -1,    51,    29,    30,    -1,    -1,    33,
      34,    35,    29,    30,    -1,    -1,    33,    34,    35,    -1,
      -1,    29,    30,    -1,    -1,    49,    -1,    51,    36,    37,
      38,    39,    49,    -1,    51,    29,    30,    29,    30,    33,
      34,    35,    -1,    51,    36,    37,    38,    39,    29,    30,
      -1,    -1,    33,    34,    35,    49,    -1,    51,    -1,    51,
      -1,    -1,    -1,    -1,    -1,    46,    29,    30,    -1,    -1,
      51,    -1,    -1,    36,    37,    38,    39,    40,    41,    42,
      -1,    29,    30,    46,    -1,    -1,    -1,    50,    36,    37,
      38,    39,    40,    41,    42,    -1,    29,    30,    46,    -1,
      -1,    -1,    50,    36,    37,    38,    39,    40,    41,    42,
      -1,    29,    30,    46,    -1,    -1,    -1,    50,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    -1,    50,    -1,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    -1,    50,    -1,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    -1,    50,    -1,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    -1,
      50,    -1,    36,    37,    38,    39,    29,    30,    -1,    -1,
      33,    34,    35,    -1,    -1,    -1,    50,    -1,    -1,    -1,
      -1,    29,    30,    46,    -1,    -1,    -1,    50,    36,    37,
      38,    39,    40,    41,    42,    29,    30,    -1,    46,    -1,
      -1,    49,    36,    37,    38,    39,    40,    41,    42,    29,
      30,    -1,    46,    -1,    -1,    49,    36,    37,    38,    39,
      40,    41,    42,    29,    30,    -1,    46,    -1,    -1,    49,
      36,    37,    38,    39,    40,    41,    42,    29,    30,    -1,
      46,    -1,    -1,    49,    36,    37,    38,    39,    40,    41,
      42,    29,    30,    -1,    46,    -1,    -1,    49,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    -1,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    -1,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    -1,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    49,
      -1,    -1,    36,    37,    38,    39,    -1,    41,    42,    -1,
      44,    29,    30,    -1,    -1,    49,    -1,    -1,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    -1,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    -1,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    -1,    36,    37,    38,    39,
      -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,    49,
      -1,    -1,    36,    37,    38,    39,    -1,    41,    42,    -1,
      44,    29,    30,    -1,    -1,    49,    -1,    -1,    36,    37,
      38,    39,    -1,    41,    42,    -1,    44,    29,    30,    -1,
      -1,    49,    -1,    -1,    36,    37,    38,    39,    -1,    41,
      42,    -1,    44,    29,    30,    -1,    -1,    49,    -1,    -1,
      36,    37,    38,    39,    -1,    41,    42,    -1,    44,    29,
      30,    -1,    -1,    49,    -1,    -1,    36,    37,    38,    39,
      29,    30,    -1,    -1,    33,    34,    35,    -1,    -1,    49,
      -1,    -1,    -1,    -1,    29,    30,    -1,    46,    29,    30,
      49,    36,    37,    38,    39,    36,    37,    38,    39,    -1,
      -1,    -1,    29,    30,    49,    -1,    -1,    -1,    49,    36,
      37,    38,    39,    -1,    41,    42,    -1,    44,    29,    30,
      47,    -1,    -1,    -1,    -1,    36,    37,    38,    39,    -1,
      41,    42,    -1,    44,    29,    30,    47,    -1,    -1,    -1,
      -1,    36,    37,    38,    39,    -1,    41,    42,    -1,    44,
      29,    30,    47,    -1,    -1,    -1,    -1,    36,    37,    38,
      39,    -1,    41,    42,    -1,    44,    29,    30,    47,    -1,
      -1,    -1,    -1,    36,    37,    38,    39,    -1,    41,    42,
      -1,    44,    29,    30,    47,    -1,    -1,    -1,    -1,    36,
      37,    38,    39,    -1,    41,    42,    -1,    44,    29,    30,
      47,    -1,    -1,    -1,    -1,    36,    37,    38,    39,    -1,
      41,    42,    -1,    44,    29,    30,    47,    -1,    -1,    -1,
      -1,    36,    37,    38,    39,    -1,    41,    42,    -1,    44,
      29,    30,    47,    -1,    -1,    -1,    -1,    36,    37,    38,
      39,    -1,    41,    42,    -1,    44,    29,    30,    47,    -1,
      -1,    -1,    -1,    36,    37,    38,    39,    -1,    41,    42,
      -1,    44,    29,    30,    47,    -1,    -1,    -1,    -1,    36,
      37,    38,    39,    40,    41,    42,    29,    30,    -1,    46,
      -1,    -1,    -1,    36,    37,    38,    39,    40,    41,    42,
      29,    30,    -1,    46,    -1,    -1,    -1,    36,    37,    38,
      39,    -1,    41,    42,    -1,    44,    29,    30,    -1,    -1,
      -1,    -1,    -1,    36,    37,    38,    39,    -1,    41,    42,
      -1,    44
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,     4,     5,     6,     7,     8,     9,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    31,    32,    38,    45,    48,
      52,    55,    56,    57,    58,    59,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    48,    48,    48,
      48,    48,    48,    48,    48,    48,    48,    56,    57,    58,
      59,    76,    77,    78,    56,    57,    58,    76,    77,    78,
      56,    57,    58,    59,    76,    77,    78,     0,    29,    30,
      33,    34,    35,    29,    30,    36,    37,    38,    39,    41,
      42,    44,    29,    30,    36,    37,    38,    39,    41,    42,
      44,    29,    30,    36,    37,    38,    39,    46,    46,    46,
      46,    46,    46,    46,    56,    57,    58,    51,    53,    51,
      53,    51,    53,    29,    30,    36,    37,    38,    39,    40,
      41,    42,    46,    29,    30,    36,    37,    38,    39,    40,
      41,    42,    46,    29,    30,    33,    34,    35,    46,    10,
      11,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    76,    77,    78,    56,    57,    58,    59,    60,    61,
      62,    63,    64,    65,    76,    77,    78,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    76,    77,    78,
      56,    57,    58,    59,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    76,    77,    78,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    76,
      77,    78,    56,    57,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    76,    77,    78,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    76,    77,
      78,    57,    58,    59,    76,    77,    57,    58,    59,    76,
      77,    76,    77,    78,    49,    49,    49,    49,    49,    49,
      49,    56,    56,    56,    57,    58,    59,    76,    77,    78,
      56,    56,    57,    57,    57,    57,    57,    58,    57,    58,
      57,    58,    76,    77,    57,    58,    58,    58,    58,    58,
      58,    57,    58,    57,    58,    57,    58,    76,    57,    58,
      57,    58,    59,    59,    59,    59,    57,    59,    57,    57,
      57,    57,    57,    57,    57,    57,    57,    58,    57,    56,
      76,    76,    58,    76,    58,    76,    57,    58,    76,    57,
      58,    76,    76,    57,    58,    76,    57,    58,    76,    57,
      77,    77,    57,    77,    57,    77,    57,    77,    57,    77,
      77,    57,    77,    57,    77,    57,    78,    78,    56,    57,
      58,    76,    77,    78,    78,    78,    57,    49,    51,    49,
      51,    49,    51,    49,    51,    49,    51,    49,    51,    46,
      46,    46,    46,    49,    51,    49,    51,    49,    51,    49,
      51,    49,    51,    49,    51,    49,    51,    49,    51,    49,
      51,    46,    46,    46,    46,    49,    51,    49,    51,    49,
      51,    49,    51,    49,    51,    49,    51,    49,    51,    49,
      51,    49,    51,    46,    46,    46,    46,    49,    51,    49,
      51,    49,    51,    49,    51,    49,    51,    49,    51,    49,
      51,    49,    51,    49,    51,    46,    46,    46,    46,    46,
      46,    46,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      46,    46,    46,    46,    46,    46,    46,    49,    51,    49,
      51,    49,    51,    49,    51,    49,    51,    49,    51,    49,
      51,    49,    51,    49,    51,    46,    46,    46,    46,    46,
      46,    46,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      46,    46,    46,    46,    49,    51,    49,    51,    49,    51,
      51,    51,    51,    49,    49,    51,    51,    51,    49,    49,
      51,    51,    51,    50,    50,    50,    50,    50,    50,    50,
      47,    47,    47,    47,    47,    47,    47,    47,    47,    50,
      50,    50,    50,    50,    50,    47,    47,    47,    47,    47,
      47,    47,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    47,    57,    58,    57,    58,    59,    57,    58,
      57,    58,    59,    57,    57,    57,    56,    57,    58,    57,
      58,    59,    76,    77,    78,    78,    77,    76,    58,    76,
      57,    77,    56,    78,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    49,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    51,    49,    51,
      49,    51,    49,    51,    49,    51,    49,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    51,    51,    51,    57,
      57,    57,    49,    49,    49
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrlab1


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)         \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (cinluded).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  /* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 299 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].bExp;;}
    break;

  case 3:
#line 300 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].iExp;;}
    break;

  case 4:
#line 301 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].dExp; ;}
    break;

  case 5:
#line 302 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].dtExp;}
    break;

  case 6:
#line 303 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].dExpArray;;}
    break;

  case 7:
#line 304 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].iExpArray;;}
    break;

  case 8:
#line 305 "FRParser.ypp"
    {MY_PARAM.finalVal = yyvsp[0].bExpArray;;}
    break;

  case 9:
#line 306 "FRParser.ypp"
    {YYABORT ;}
    break;

  case 10:
#line 307 "FRParser.ypp"
    {YYABORT ;}
    break;

  case 11:
#line 310 "FRParser.ypp"
    {yyval.bExp = yyvsp[0].bExp;}
    break;

  case 12:
#line 311 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createIndexFunc<FR::RValueBool, bool, 
             FRIfaces::IRValueBool>(yyvsp[-3].var,yyvsp[-1].iExp,MY_PARAM.theCtrl);} CHK_CATCH(yyval.bExp);}
    break;

  case 13:
#line 314 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createArrayBracketOperator<FR::RValueBool, bool, 
             FRIfaces::IRValueBoolArray::RT, FRIfaces::IRValueBool::RT>(
                 yyvsp[-3].bExpArray->getRT(), yyvsp[-1].iExp->getRT());
    } CHK_CATCH(yyval.bExp);}
    break;

  case 14:
#line 319 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createUnaryFunc<FR::RValueBool, bool,
             FR::RConstBool, NotOp, FR::LValueBool>(yyvsp[0].bExp);
    } CHK_CATCH(yyval.bExp);}
    break;

  case 15:
#line 323 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryBoolFunc<EqOp<bool> >(yyvsp[-2].bExp, yyvsp[0].bExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 16:
#line 325 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createBinaryBoolFunc< NotEqOp<bool> >(yyvsp[-2].bExp, yyvsp[0].bExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 17:
#line 327 "FRParser.ypp"
    {yyval.bExp = yyvsp[-1].bExp;}
    break;

  case 18:
#line 328 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryInt2BoolFunc< EqOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 19:
#line 330 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryInt2BoolFunc< NotEqOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 20:
#line 332 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createBinaryDouble2BoolFunc< EqOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.bExp);}
    break;

  case 21:
#line 335 "FRParser.ypp"
    { 
    TRY{ yyval.bExp=createBinaryDouble2BoolFunc<NotEqOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.bExp);}
    break;

  case 22:
#line 338 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createBinaryDate2BoolFunc<
             EqOp<const DateTime::Date&> >(yyvsp[-2].dtExp, yyvsp[0].dtExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 23:
#line 341 "FRParser.ypp"
    { 
    TRY{ yyval.bExp= createBinaryDate2BoolFunc<
             NotEqOp<const DateTime::Date&> >(yyvsp[-2].dtExp, yyvsp[0].dtExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 24:
#line 344 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryInt2BoolFunc< LTOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 25:
#line 346 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryInt2BoolFunc<GTOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 26:
#line 348 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryDouble2BoolFunc<LTOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.bExp);}
    break;

  case 27:
#line 351 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryDouble2BoolFunc<GTOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.bExp);}
    break;

  case 28:
#line 354 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryDate2BoolFunc<
             LTOp<const DateTime::Date&> >(yyvsp[-2].dtExp, yyvsp[0].dtExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 29:
#line 357 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBinaryDate2BoolFunc<
             GTOp<const DateTime::Date&> >(yyvsp[-2].dtExp, yyvsp[0].dtExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 30:
#line 360 "FRParser.ypp"
    { 
    TRY{ yyval.bExp = createBoolAnd(yyvsp[-2].bExp, yyvsp[0].bExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 31:
#line 362 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createBoolOr(yyvsp[-2].bExp, yyvsp[0].bExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 32:
#line 364 "FRParser.ypp"
    {
    TRY{ yyval.bExp = createBoolCondition(yyvsp[-4].bExp, yyvsp[-2].bExp, yyvsp[0].bExp);} CHK_CATCH(yyval.bExp);}
    break;

  case 33:
#line 366 "FRParser.ypp"
    {
    TRY{ yyval.bExp = FRInternalFunctionBool::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.bExp);}
    break;

  case 34:
#line 371 "FRParser.ypp"
    { yyval.iExp = yyvsp[0].iExp;;}
    break;

  case 35:
#line 372 "FRParser.ypp"
    {
    TRY{ yyval.iExp = createIndexFunc<FR::RValueInt, 
             int, FRIfaces::IRValueInt>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);
    } CHK_CATCH(yyval.iExp);}
    break;

  case 36:
#line 376 "FRParser.ypp"
    {
    TRY{ yyval.iExp = createArrayBracketOperator<FR::RValueInt, int,
             FRIfaces::IRValueIntArray::RT, FRIfaces::IRValueInt::RT>(
                 yyvsp[-3].iExpArray->getRT(), yyvsp[-1].iExp->getRT());
    } CHK_CATCH(yyval.iExp);}
    break;

  case 37:
#line 381 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryIntFunc<AddOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 38:
#line 383 "FRParser.ypp"
    {
    TRY{ yyval.iExp = createBinaryIntFunc<SubOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 39:
#line 385 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryIntFunc<MultOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 40:
#line 387 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryIntFunc<DivOp<int> >(yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 41:
#line 389 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createUnaryFunc<FR::RValueInt, int,
             FR::RConstInt, NegOp<int>,
             FR::LValueInt>(yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 42:
#line 393 "FRParser.ypp"
    { yyval.iExp = yyvsp[-1].iExp;;}
    break;

  case 43:
#line 394 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryFunc<FRIfaces::IRValueInt, 
             FR::RValueInt, const DateTime::Date&,
             const DateTime::Date&, int, FR::RConstInt, 
             SubOp<const DateTime::Date&>,
             FRIfaces::IRValueDate, FRIfaces::IRValueDate >(yyvsp[-2].dtExp, yyvsp[0].dtExp);
    } CHK_CATCH(yyval.iExp);}
    break;

  case 44:
#line 401 "FRParser.ypp"
    {
    TRY{ yyval.iExp = createIntCondition(yyvsp[-4].bExp, yyvsp[-2].iExp, yyvsp[0].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 45:
#line 403 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryIntFunc<MaxOp<int> >(yyvsp[-3].iExp, yyvsp[-1].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 46:
#line 405 "FRParser.ypp"
    { 
    TRY{ yyval.iExp = createBinaryIntFunc<MinOp<int> >(yyvsp[-3].iExp, yyvsp[-1].iExp);} CHK_CATCH(yyval.iExp);}
    break;

  case 47:
#line 407 "FRParser.ypp"
    {
    TRY{ yyval.iExp = FRInternalFunctionInt::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.iExp);}
    break;

  case 48:
#line 409 "FRParser.ypp"
    {
    TRY{ yyval.iExp = new Binary<FR::RValueInt, int, 
             FRIfaces::IRValueIntArray::RT*, 
             FRIfaces::IRValueIntArray::RT*, 
             ArrayMultOp<int, FRIfaces::IRValueIntArray::RT*>,
             FRIfaces::IRValueInt::RT, 
             FRIfaces::IRValueIntArray::RT*,
             FRIfaces::IRValueIntArray::RT*>(yyvsp[-2].iExpArray->getRT(), yyvsp[0].iExpArray->getRT());
    } CHK_CATCH(yyval.iExp);}
    break;

  case 49:
#line 418 "FRParser.ypp"
    {
    TRY{ yyval.iExp = new Unary<FR::RValueInt, int, 
             FRIfaces::IRValueIntArray::RT*,
             MaxMinArrayOp<int, FRIfaces::IRValueIntArray::RT*,
             GTOp<int> >, 
             FRIfaces::IRValueInt::RT, 
             FRIfaces::IRValueIntArray::RT*>(yyvsp[-1].iExpArray->getRT()); } CHK_CATCH(yyval.iExp);}
    break;

  case 50:
#line 425 "FRParser.ypp"
    {
    TRY{ yyval.iExp = new Unary<FR::RValueInt, int, 
             FRIfaces::IRValueIntArray::RT*,
             MaxMinArrayOp<int, FRIfaces::IRValueIntArray::RT*,
             LTOp<int> >, 
             FRIfaces::IRValueInt::RT, 
             FRIfaces::IRValueIntArray::RT*>(yyvsp[-1].iExpArray->getRT()); } CHK_CATCH(yyval.iExp);}
    break;

  case 51:
#line 435 "FRParser.ypp"
    { yyval.dExp = yyvsp[0].dExp;}
    break;

  case 52:
#line 436 "FRParser.ypp"
    {
    TRY{ yyval.dExp = createIndexFunc<FR::RValueDouble, 
             double, FRIfaces::IRValueDouble>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);
    }CHK_CATCH(yyval.dExp);}
    break;

  case 53:
#line 440 "FRParser.ypp"
    {
    TRY{ yyval.dExp = createArrayBracketOperator<FR::RValueDouble, double, 
             FRIfaces::IRValueDoubleArray::RT, FRIfaces::IRValueDouble::RT>(
                 yyvsp[-3].dExpArray->getRT(), yyvsp[-1].iExp->getRT());
    } CHK_CATCH(yyval.dExp);}
    break;

  case 54:
#line 445 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFunc<AddOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 55:
#line 447 "FRParser.ypp"
    {
    TRY{ yyval.dExp = createBinaryDoubleFuncInt2<AddOp<double> >(yyvsp[-2].dExp, yyvsp[0].iExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 56:
#line 450 "FRParser.ypp"
    {
    TRY{ yyval.dExp = createBinaryDoubleFuncInt1<AddOp<double> >(yyvsp[-2].iExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 57:
#line 454 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFunc<SubOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 58:
#line 456 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt2<SubOp<double> >(yyvsp[-2].dExp, yyvsp[0].iExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 59:
#line 459 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt1<SubOp<double> >(yyvsp[-2].iExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 60:
#line 463 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFunc<MultOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 61:
#line 465 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt2<MultOp<double> >(yyvsp[-2].dExp, yyvsp[0].iExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 62:
#line 468 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt1<MultOp<double> >(yyvsp[-2].iExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 63:
#line 472 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFunc<DivOp<double> >(yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 64:
#line 474 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt2<DivOp<double> >(yyvsp[-2].dExp, yyvsp[0].iExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 65:
#line 477 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt1<DivOp<double> >(yyvsp[-2].iExp, yyvsp[0].dExp);
    } CHK_CATCH(yyval.dExp);}
    break;

  case 66:
#line 481 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFunc<PowOp>(yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 67:
#line 483 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt2<PowOp>(yyvsp[-2].dExp, yyvsp[0].iExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 68:
#line 485 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createBinaryDoubleFuncInt1<PowOp>(yyvsp[-2].iExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 69:
#line 488 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createUnaryFunc< FR::RValueDouble, double, 
             FR::RConstDouble, NegOp<double>, 
             FR::LValueDouble >(yyvsp[0].dExp); } CHK_CATCH(yyval.dExp);}
    break;

  case 70:
#line 492 "FRParser.ypp"
    { yyval.dExp = yyvsp[-1].dExp;;}
    break;

  case 71:
#line 494 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFunc<MaxOp<double> >(yyvsp[-3].dExp, yyvsp[-1].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 72:
#line 496 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFuncInt1<MaxOp<double> >(yyvsp[-3].iExp, yyvsp[-1].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 73:
#line 498 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFuncInt2<MaxOp<double> >(yyvsp[-3].dExp, yyvsp[-1].iExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 74:
#line 501 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFunc<MinOp<double> >(yyvsp[-3].dExp, yyvsp[-1].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 75:
#line 503 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFuncInt1<MinOp<double> >(yyvsp[-3].iExp, yyvsp[-1].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 76:
#line 505 "FRParser.ypp"
    {
    TRY{ yyval.dExp=createBinaryDoubleFuncInt2<MinOp<double> >(yyvsp[-3].dExp, yyvsp[-1].iExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 77:
#line 508 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = createDoubleCondition(yyvsp[-4].bExp, yyvsp[-2].dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 78:
#line 510 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = new FRIntToDouble(yyvsp[-2].iExp); CHECK_PTR(yyval.dExp);
    yyval.dExp = createDoubleCondition(yyvsp[-4].bExp, yyval.dExp, yyvsp[0].dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 79:
#line 513 "FRParser.ypp"
    { 
    TRY{ yyval.dExp = new FRIntToDouble(yyvsp[0].iExp); CHECK_PTR(yyval.dExp);
    yyval.dExp = createDoubleCondition(yyvsp[-4].bExp, yyvsp[-2].dExp, yyval.dExp);} CHK_CATCH(yyval.dExp);}
    break;

  case 80:
#line 516 "FRParser.ypp"
    {
    TRY{ yyval.dExp = FRInternalFunctionDouble::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.dExp);}
    break;

  case 81:
#line 518 "FRParser.ypp"
    {
    TRY{ yyval.dExp = new Binary<FR::RValueDouble, double, 
             FRIfaces::IRValueDoubleArray::RT*, 
             FRIfaces::IRValueDoubleArray::RT*, 
             ArrayMultOp<double, FRIfaces::IRValueDoubleArray::RT*>,
             FRIfaces::IRValueDouble::RT, 
             FRIfaces::IRValueDoubleArray::RT*,
             FRIfaces::IRValueDoubleArray::RT*>(yyvsp[-2].dExpArray->getRT(), yyvsp[0].dExpArray->getRT());
    } CHK_CATCH(yyval.dExp);}
    break;

  case 82:
#line 527 "FRParser.ypp"
    {
    TRY{ yyval.dExp = new Unary<FR::RValueDouble, double, 
             FRIfaces::IRValueDoubleArray::RT*,
             MaxMinArrayOp<double, FRIfaces::IRValueDoubleArray::RT*,
             GTOp<double> >, 
             FRIfaces::IRValueDouble::RT, 
             FRIfaces::IRValueDoubleArray::RT*>(yyvsp[-1].dExpArray->getRT()); } CHK_CATCH(yyval.dExp);}
    break;

  case 83:
#line 534 "FRParser.ypp"
    {
    TRY{ yyval.dExp = new Unary<FR::RValueDouble, double, 
             FRIfaces::IRValueDoubleArray::RT*,
             MaxMinArrayOp<double, FRIfaces::IRValueDoubleArray::RT*,
             LTOp<double> >, 
             FRIfaces::IRValueDouble::RT, 
             FRIfaces::IRValueDoubleArray::RT*>(yyvsp[-1].dExpArray->getRT()); } CHK_CATCH(yyval.dExp);}
    break;

  case 84:
#line 544 "FRParser.ypp"
    { yyval.dtExp = yyvsp[0].dtExp;}
    break;

  case 85:
#line 545 "FRParser.ypp"
    {
    TRY{ yyval.dtExp = createIndexFunc<FR::RValueDate, const DateTime::Date&,
             FRIfaces::IRValueDate>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);} CHK_CATCH(yyval.dtExp);}
    break;

  case 86:
#line 548 "FRParser.ypp"
    {
    TRY{ yyval.dtExp = createBinaryDateFuncInt2<AddOp<const DateTime::Date&, int>,
             const DateTime::Date&, int>(yyvsp[-2].dtExp, yyvsp[0].iExp);} CHK_CATCH(yyval.dtExp);}
    break;

  case 87:
#line 551 "FRParser.ypp"
    {
    TRY{ yyval.dtExp = createBinaryDateFuncInt2<SubOp<const DateTime::Date&, int>,
             const DateTime::Date&, int>(yyvsp[-2].dtExp, yyvsp[0].iExp);} CHK_CATCH(yyval.dtExp);}
    break;

  case 88:
#line 554 "FRParser.ypp"
    { yyval.dtExp = yyvsp[-1].dtExp;}
    break;

  case 89:
#line 555 "FRParser.ypp"
    {
    TRY{ yyval.dtExp=createBinaryDateFunc< MaxOp<const DateTime::Date&>,
             const DateTime::Date&, const DateTime::Date&>(yyvsp[-3].dtExp, yyvsp[-1].dtExp);
    } CHK_CATCH(yyval.dtExp);}
    break;

  case 90:
#line 559 "FRParser.ypp"
    {
    TRY{ yyval.dtExp=createBinaryDateFunc< MinOp<const DateTime::Date&>, 
             const DateTime::Date&, const DateTime::Date&>(yyvsp[-3].dtExp, yyvsp[-1].dtExp);
    } CHK_CATCH(yyval.dtExp);}
    break;

  case 91:
#line 563 "FRParser.ypp"
    { 
    TRY{ yyval.dtExp = createConditionFunc< FRIfaces::IRValueDate, FR::RValueDate,
             const DateTime::Date&, FRIfaces::IRValueDate>(yyvsp[-4].bExp, yyvsp[-2].dtExp, yyvsp[0].dtExp);
    } CHK_CATCH(yyval.dtExp);}
    break;

  case 92:
#line 567 "FRParser.ypp"
    { 
    TRY{ yyval.dtExp = FRInternalFunctionDate::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.dtExp);}
    break;

  case 93:
#line 572 "FRParser.ypp"
    {yyval.schedExp = yyvsp[0].schedExp;}
    break;

  case 94:
#line 575 "FRParser.ypp"
    {yyval.tabFuncExp = yyvsp[0].tabFuncExp;}
    break;

  case 95:
#line 578 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 96:
#line 581 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 97:
#line 584 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 98:
#line 587 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 99:
#line 590 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 100:
#line 592 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 101:
#line 594 "FRParser.ypp"
    {yyval.var = yyvsp[0].var;}
    break;

  case 103:
#line 609 "FRParser.ypp"
    {
    TRY{ yyval.dValArray = new FR::RDoubleArray(); yyval.dValArray->addArg(yyvsp[0].dExp);} CHK_CATCH(yyval.dValArray);
;}
    break;

  case 104:
#line 612 "FRParser.ypp"
    { 
    TRY{ yyval.dValArray = yyvsp[-2].dValArray; yyvsp[-2].dValArray->addArg(yyvsp[0].dExp);} CATCH;}
    break;

  case 105:
#line 614 "FRParser.ypp"
    { 
    TRY{ // need to cast int to double
        FRIfaces::IRValueDouble* rValue = new FRIntToDouble(yyvsp[0].iExp);
        CHECK_PTR(rValue);
        yyval.dValArray = yyvsp[-2].dValArray; yyvsp[-2].dValArray->addArg(rValue);} CATCH;}
    break;

  case 106:
#line 621 "FRParser.ypp"
    {yyval.dValArray=yyvsp[-1].dValArray;}
    break;

  case 107:
#line 625 "FRParser.ypp"
    {
    TRY{ yyval.iValArray = new FR::RIntArray(); yyval.iValArray->addArg(yyvsp[0].iExp);} CHK_CATCH(yyval.iValArray);
;}
    break;

  case 108:
#line 628 "FRParser.ypp"
    { 
    TRY{ yyval.iValArray = yyvsp[-2].iValArray; yyvsp[-2].iValArray->addArg(yyvsp[0].iExp);} CATCH;}
    break;

  case 109:
#line 632 "FRParser.ypp"
    {yyval.iValArray=yyvsp[-1].iValArray;}
    break;

  case 110:
#line 636 "FRParser.ypp"
    {
    TRY{ yyval.bValArray = new FR::RBoolArray(); yyval.bValArray->addArg(yyvsp[0].bExp);} CHK_CATCH(yyval.bValArray);
;}
    break;

  case 111:
#line 639 "FRParser.ypp"
    { 
    TRY{ yyval.bValArray = yyvsp[-2].bValArray; yyvsp[-2].bValArray->addArg(yyvsp[0].bExp);} CATCH;}
    break;

  case 112:
#line 643 "FRParser.ypp"
    {yyval.bValArray=yyvsp[-1].bValArray;}
    break;

  case 113:
#line 648 "FRParser.ypp"
    {yyval.dExpArray=yyvsp[0].dValArray;}
    break;

  case 114:
#line 649 "FRParser.ypp"
    {yyval.dExpArray=yyvsp[0].dExpArray;}
    break;

  case 115:
#line 650 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = createArrayIndexFunc<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 116:
#line 654 "FRParser.ypp"
    { // add a double to an array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             AddOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 117:
#line 660 "FRParser.ypp"
    { // add an int to a double array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueInt,
             AddOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 118:
#line 666 "FRParser.ypp"
    { // add an array to an array
    TRY{ yyval.dExpArray = new ArrayBinary<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             AddOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 119:
#line 672 "FRParser.ypp"
    { // subtract a double from an array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             SubOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 120:
#line 678 "FRParser.ypp"
    { // subtract an int from a double array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueInt,
             SubOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 121:
#line 684 "FRParser.ypp"
    { // subtract an array from an array
    TRY{ yyval.dExpArray = new ArrayBinary<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             SubOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 122:
#line 690 "FRParser.ypp"
    { // multiply an array by a scalar
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             MultOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 123:
#line 696 "FRParser.ypp"
    { // multiply an array by an int scalar
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueInt,
             MultOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 124:
#line 702 "FRParser.ypp"
    { // multiply a scalar by an array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             MultOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[0].dExpArray, yyvsp[-2].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 125:
#line 708 "FRParser.ypp"
    { // multiply a scalar int by an array
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueInt,
             MultOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[0].dExpArray, yyvsp[-2].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 126:
#line 714 "FRParser.ypp"
    { // multiple element by element an array by an array
    TRY{ yyval.dExpArray = new ArrayBinary<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             MultOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 127:
#line 720 "FRParser.ypp"
    { // divide an array by a scalar
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             DivOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 128:
#line 726 "FRParser.ypp"
    { // divide an array by an int scalar
    TRY{ yyval.dExpArray = new ArrayBinaryOnScalar<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueInt,
             DivOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 129:
#line 732 "FRParser.ypp"
    { // divide element by element an array by an array
    TRY{ yyval.dExpArray = new ArrayBinary<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             DivOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 130:
#line 738 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new ArrayUnary<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray,
             NegOp<double>, FRIfaces::IRValueDoubleArray::RT>(yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 131:
#line 743 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = FRInternalFunctionDoubleArray::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.dExpArray);}
    break;

  case 132:
#line 745 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new ArrayCondition<FR::RValueDoubleArray, 
             FRIfaces::IRValueDoubleArray::RT, double,
             FRIfaces::IRValueDoubleArray>(yyvsp[-4].bExpArray, yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 133:
#line 750 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new ArrayConditionUsingScalar<FR::RValueDoubleArray, 
             FRIfaces::IRValueDoubleArray::RT, double,
             FRIfaces::IRValueDoubleArray, FRIfaces::IRValueDouble,
             false>(yyvsp[-4].bExpArray, yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 134:
#line 756 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new ArrayConditionUsingScalar<FR::RValueDoubleArray, 
             FRIfaces::IRValueDoubleArray::RT, double,
             FRIfaces::IRValueDoubleArray, FRIfaces::IRValueDouble,
             true>(yyvsp[-4].bExpArray, yyvsp[0].dExpArray, yyvsp[-2].dExp); // reverse parameters and negate condition
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 135:
#line 762 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new ArrayConditionFromBool<FR::RValueDoubleArray, 
             FRIfaces::IRValueDoubleArray::RT, double,
             FRIfaces::IRValueDoubleArray>(yyvsp[-4].bExp, yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 136:
#line 767 "FRParser.ypp"
    {
    TRY{ yyval.dExpArray = new SubArray<FR::RValueDoubleArray, double,
             FRIfaces::IRValueDoubleArray, 
             FRIfaces::IRValueDoubleArray::RT>(yyvsp[-5].dExpArray, yyvsp[-3].iExp, yyvsp[-1].iExp);
    }CHK_CATCH(yyval.dExpArray);}
    break;

  case 137:
#line 772 "FRParser.ypp"
    { yyval.dExpArray = yyvsp[-1].dExpArray;;}
    break;

  case 138:
#line 777 "FRParser.ypp"
    {yyval.iExpArray=yyvsp[0].iValArray;}
    break;

  case 139:
#line 778 "FRParser.ypp"
    {yyval.iExpArray=yyvsp[0].iExpArray;}
    break;

  case 140:
#line 779 "FRParser.ypp"
    {
    TRY{ yyval.iExpArray = createArrayIndexFunc<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 141:
#line 783 "FRParser.ypp"
    { // add an int to an array
    TRY{ yyval.iExpArray = new ArrayBinaryOnScalar<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             AddOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 142:
#line 789 "FRParser.ypp"
    { // add an array to an array
    TRY{ yyval.iExpArray = new ArrayBinary<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             AddOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 143:
#line 795 "FRParser.ypp"
    { // subtract a int from an array
    TRY{ yyval.iExpArray = new ArrayBinaryOnScalar<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             SubOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 144:
#line 801 "FRParser.ypp"
    { // subtract an array from an array
    TRY{ yyval.iExpArray = new ArrayBinary<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             SubOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 145:
#line 807 "FRParser.ypp"
    { // multiply an array by a scalar
    TRY{ yyval.iExpArray = new ArrayBinaryOnScalar<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             MultOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 146:
#line 813 "FRParser.ypp"
    { // multiply a scalar by an array
    TRY{ yyval.iExpArray = new ArrayBinaryOnScalar<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             MultOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[0].iExpArray, yyvsp[-2].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 147:
#line 819 "FRParser.ypp"
    { // multiple element by element an array by an array
    TRY{ yyval.iExpArray = new ArrayBinary<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             MultOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 148:
#line 825 "FRParser.ypp"
    { // divide an array by a scalar
    TRY{ yyval.iExpArray = new ArrayBinaryOnScalar<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             DivOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 149:
#line 831 "FRParser.ypp"
    { // divide element by element an array by an array
    TRY{ yyval.iExpArray = new ArrayBinary<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             DivOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 150:
#line 837 "FRParser.ypp"
    {
    TRY{ yyval.iExpArray = new ArrayUnary<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray,
             NegOp<int>, FRIfaces::IRValueIntArray::RT>(yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 151:
#line 842 "FRParser.ypp"
    {
    TRY{ yyval.iExpArray = FRInternalFunctionIntArray::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.iExpArray);}
    break;

  case 152:
#line 844 "FRParser.ypp"
    {
TRY{ yyval.iExpArray = new ArrayCondition<FR::RValueIntArray, 
         FRIfaces::IRValueIntArray::RT, int,
         FRIfaces::IRValueIntArray>(yyvsp[-4].bExpArray, yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 153:
#line 849 "FRParser.ypp"
    {
TRY{ yyval.iExpArray = new ArrayConditionUsingScalar<FR::RValueIntArray, 
         FRIfaces::IRValueIntArray::RT, int,
         FRIfaces::IRValueIntArray, FRIfaces::IRValueInt,
         false>(yyvsp[-4].bExpArray, yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 154:
#line 855 "FRParser.ypp"
    {
TRY{ yyval.iExpArray = new ArrayConditionUsingScalar<FR::RValueIntArray, 
         FRIfaces::IRValueIntArray::RT, int,
         FRIfaces::IRValueIntArray, FRIfaces::IRValueInt,
         true>(yyvsp[-4].bExpArray, yyvsp[0].iExpArray, yyvsp[-2].iExp); // reverse parameters and negate condition
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 155:
#line 861 "FRParser.ypp"
    {
    TRY{ yyval.iExpArray = new ArrayConditionFromBool<FR::RValueIntArray, 
             FRIfaces::IRValueIntArray::RT, int,
             FRIfaces::IRValueIntArray>(yyvsp[-4].bExp, yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 156:
#line 866 "FRParser.ypp"
    {
    TRY{ yyval.iExpArray = new SubArray<FR::RValueIntArray, int,
             FRIfaces::IRValueIntArray, 
             FRIfaces::IRValueIntArray::RT>(yyvsp[-5].iExpArray, yyvsp[-3].iExp, yyvsp[-1].iExp);
    }CHK_CATCH(yyval.iExpArray);}
    break;

  case 157:
#line 871 "FRParser.ypp"
    { yyval.iExpArray = yyvsp[-1].iExpArray;;}
    break;

  case 158:
#line 876 "FRParser.ypp"
    {yyval.bExpArray=yyvsp[0].bValArray;}
    break;

  case 159:
#line 877 "FRParser.ypp"
    {yyval.bExpArray=yyvsp[0].bExpArray;}
    break;

  case 160:
#line 878 "FRParser.ypp"
    { 
    TRY{ yyval.bExpArray = new ArrayUnary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueBoolArray,
             NotOp, FRIfaces::IRValueBoolArray::RT>(yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 161:
#line 883 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = createArrayIndexFunc<FR::RValueBoolArray, bool,
             FRIfaces::IRValueBoolArray>(yyvsp[-3].var, yyvsp[-1].iExp, MY_PARAM.theCtrl);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 162:
#line 887 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = FRInternalFunctionBoolArray::checkArgs(yyvsp[0].funcArgs);} CHK_CATCH(yyval.bExpArray);}
    break;

  case 163:
#line 889 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             GTOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 164:
#line 895 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             LTOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 165:
#line 901 "FRParser.ypp"
    { // compare array with scalar, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinaryOnScalar<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             GTOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 166:
#line 907 "FRParser.ypp"
    { // compare array with scalar, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinaryOnScalar<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDouble,
             LTOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 167:
#line 913 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             GTOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 168:
#line 919 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             LTOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 169:
#line 925 "FRParser.ypp"
    { // compare array with scalar, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinaryOnScalar<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             GTOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 170:
#line 931 "FRParser.ypp"
    { // compare array with scalar, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinaryOnScalar<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueInt,
             LTOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 171:
#line 937 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayOr(yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 172:
#line 940 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayAnd(yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 173:
#line 943 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueBoolArray,
             FRIfaces::IRValueBoolArray,
             EqOp<bool>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 174:
#line 949 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueBoolArray,
             FRIfaces::IRValueBoolArray,
             NotEqOp<bool>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 175:
#line 955 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             EqOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 176:
#line 961 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueIntArray,
             FRIfaces::IRValueIntArray,
             NotEqOp<int>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].iExpArray, yyvsp[0].iExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 177:
#line 967 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             EqOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 178:
#line 973 "FRParser.ypp"
    { // compare two arrays, return a bool array
    TRY{ yyval.bExpArray = new ArrayBinary<FR::RValueBoolArray, bool,
             FRIfaces::IRValueDoubleArray,
             FRIfaces::IRValueDoubleArray,
             NotEqOp<double>, FRIfaces::IRValueBoolArray::RT>(yyvsp[-2].dExpArray, yyvsp[0].dExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 179:
#line 979 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = new ArrayCondition<FR::RValueBoolArray, 
             FRIfaces::IRValueBoolArray::RT, bool,
             FRIfaces::IRValueBoolArray>(yyvsp[-4].bExpArray, yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 180:
#line 984 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = new ArrayConditionUsingScalar<FR::RValueBoolArray, 
             FRIfaces::IRValueBoolArray::RT, bool,
             FRIfaces::IRValueBoolArray, FRIfaces::IRValueBool,
             false>(yyvsp[-4].bExpArray, yyvsp[-2].bExpArray, yyvsp[0].bExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 181:
#line 990 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = new ArrayConditionUsingScalar<FR::RValueBoolArray, 
             FRIfaces::IRValueBoolArray::RT, bool,
             FRIfaces::IRValueBoolArray, FRIfaces::IRValueBool,
             true>(yyvsp[-4].bExpArray, yyvsp[0].bExpArray, yyvsp[-2].bExp); // reverse parameters and negate condition
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 182:
#line 996 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = new ArrayConditionFromBool<FR::RValueBoolArray, 
             FRIfaces::IRValueBoolArray::RT, bool,
             FRIfaces::IRValueBoolArray>(yyvsp[-4].bExp, yyvsp[-2].bExpArray, yyvsp[0].bExpArray);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 183:
#line 1001 "FRParser.ypp"
    {
    TRY{ yyval.bExpArray = new SubArray<FR::RValueBoolArray, bool,
             FRIfaces::IRValueBoolArray, 
             FRIfaces::IRValueBoolArray::RT>(yyvsp[-5].bExpArray, yyvsp[-3].iExp, yyvsp[-1].iExp);
    }CHK_CATCH(yyval.bExpArray);}
    break;

  case 184:
#line 1006 "FRParser.ypp"
    { yyval.bExpArray = yyvsp[-1].bExpArray;;}
    break;

  case 185:
#line 1011 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionBool(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 186:
#line 1017 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 187:
#line 1018 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 188:
#line 1019 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 189:
#line 1020 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 190:
#line 1021 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 191:
#line 1022 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 192:
#line 1023 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 193:
#line 1024 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 194:
#line 1025 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 195:
#line 1026 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 196:
#line 1027 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 197:
#line 1028 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 198:
#line 1029 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 199:
#line 1037 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 200:
#line 1038 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 201:
#line 1039 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 202:
#line 1040 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 203:
#line 1041 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 204:
#line 1042 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 205:
#line 1043 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 206:
#line 1044 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 207:
#line 1045 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 208:
#line 1046 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 209:
#line 1047 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 210:
#line 1048 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 211:
#line 1049 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 212:
#line 1058 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionInt(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 213:
#line 1064 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 214:
#line 1065 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 215:
#line 1066 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 216:
#line 1067 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 217:
#line 1068 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 218:
#line 1069 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 219:
#line 1070 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 220:
#line 1071 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 221:
#line 1072 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 222:
#line 1073 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 223:
#line 1074 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 224:
#line 1075 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 225:
#line 1076 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 226:
#line 1084 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 227:
#line 1085 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 228:
#line 1086 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 229:
#line 1087 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 230:
#line 1088 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 231:
#line 1089 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 232:
#line 1090 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 233:
#line 1091 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 234:
#line 1092 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 235:
#line 1093 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 236:
#line 1094 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 237:
#line 1095 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 238:
#line 1096 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 239:
#line 1105 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionDouble(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 240:
#line 1111 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 241:
#line 1112 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 242:
#line 1113 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 243:
#line 1114 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 244:
#line 1115 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 245:
#line 1116 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 246:
#line 1117 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 247:
#line 1118 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 248:
#line 1119 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 249:
#line 1120 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 250:
#line 1121 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 251:
#line 1122 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 252:
#line 1123 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 253:
#line 1131 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 254:
#line 1132 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 255:
#line 1133 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 256:
#line 1134 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 257:
#line 1135 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 258:
#line 1136 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 259:
#line 1137 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 260:
#line 1138 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 261:
#line 1139 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 262:
#line 1140 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 263:
#line 1141 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 264:
#line 1142 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 265:
#line 1143 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 266:
#line 1152 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionDoubleArray(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 267:
#line 1158 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 268:
#line 1159 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 269:
#line 1160 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 270:
#line 1161 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 271:
#line 1162 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 272:
#line 1163 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 273:
#line 1164 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 274:
#line 1165 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 275:
#line 1166 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 276:
#line 1167 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 277:
#line 1168 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 278:
#line 1169 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 279:
#line 1170 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 280:
#line 1172 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 281:
#line 1180 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 282:
#line 1181 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 283:
#line 1182 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 284:
#line 1183 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 285:
#line 1184 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 286:
#line 1185 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 287:
#line 1186 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 288:
#line 1187 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 289:
#line 1188 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 290:
#line 1189 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 291:
#line 1190 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 292:
#line 1192 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 293:
#line 1193 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 294:
#line 1195 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 295:
#line 1196 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 296:
#line 1198 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 297:
#line 1207 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionIntArray(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 298:
#line 1213 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 299:
#line 1214 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 300:
#line 1215 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 301:
#line 1216 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 302:
#line 1217 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 303:
#line 1218 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 304:
#line 1219 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 305:
#line 1220 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 306:
#line 1221 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 307:
#line 1222 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 308:
#line 1223 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 309:
#line 1225 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 310:
#line 1226 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 311:
#line 1228 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 312:
#line 1229 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 313:
#line 1231 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 314:
#line 1239 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 315:
#line 1240 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 316:
#line 1241 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 317:
#line 1242 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 318:
#line 1243 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 319:
#line 1244 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 320:
#line 1245 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 321:
#line 1246 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 322:
#line 1247 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 323:
#line 1248 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 324:
#line 1249 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 325:
#line 1251 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 326:
#line 1252 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 327:
#line 1254 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 328:
#line 1255 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 329:
#line 1257 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 330:
#line 1266 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionBoolArray(MY_PARAM.theCtrl, yyvsp[-1].func);
    } CHK_CATCH(yyval.funcArgs);}
    break;

  case 331:
#line 1272 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 332:
#line 1273 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 333:
#line 1274 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 334:
#line 1275 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 335:
#line 1276 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 336:
#line 1277 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 337:
#line 1278 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 338:
#line 1279 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 339:
#line 1280 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 340:
#line 1281 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 341:
#line 1282 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 342:
#line 1284 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 343:
#line 1285 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 344:
#line 1287 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 345:
#line 1288 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 346:
#line 1290 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 347:
#line 1298 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 348:
#line 1299 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 349:
#line 1300 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 350:
#line 1301 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 351:
#line 1302 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 352:
#line 1303 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 353:
#line 1304 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 354:
#line 1305 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 355:
#line 1306 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 356:
#line 1307 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 357:
#line 1308 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 358:
#line 1310 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 359:
#line 1311 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 360:
#line 1313 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 361:
#line 1314 "FRParser.ypp"
    { 
    TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 362:
#line 1316 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 363:
#line 1325 "FRParser.ypp"
    {
    TRY{ yyval.funcArgs = new FRInternalFunctionDate(MY_PARAM.theCtrl, yyvsp[-1].func);} CHK_CATCH(yyval.funcArgs);}
    break;

  case 364:
#line 1330 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH;}
    break;

  case 365:
#line 1331 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 366:
#line 1332 "FRParser.ypp"
    { TRY{yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH;}
    break;

  case 367:
#line 1333 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 368:
#line 1334 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH;}
    break;

  case 369:
#line 1335 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 370:
#line 1336 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH;}
    break;

  case 371:
#line 1337 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH;}
    break;

  case 372:
#line 1338 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH;}
    break;

  case 373:
#line 1339 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 374:
#line 1340 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 375:
#line 1341 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 376:
#line 1342 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 377:
#line 1350 "FRParser.ypp"
    {TRY{ yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExp);} CATCH ;}
    break;

  case 378:
#line 1351 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].bExpArray);} CATCH;}
    break;

  case 379:
#line 1352 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExp);} CATCH ;}
    break;

  case 380:
#line 1353 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].iExpArray);} CATCH;}
    break;

  case 381:
#line 1354 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExp);} CATCH ;}
    break;

  case 382:
#line 1355 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dExpArray);} CATCH;}
    break;

  case 383:
#line 1356 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].dtExp);} CATCH ;}
    break;

  case 384:
#line 1357 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].schedExp);} CATCH ;}
    break;

  case 385:
#line 1358 "FRParser.ypp"
    {TRY{  yyval.funcArgs = yyvsp[-2].funcArgs; yyval.funcArgs->addArg(yyvsp[-1].tabFuncExp);} CATCH ;}
    break;

  case 386:
#line 1359 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 387:
#line 1360 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 388:
#line 1361 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;

  case 389:
#line 1362 "FRParser.ypp"
    { TRY{ yyval.funcArgs = yyvsp[-4].funcArgs; yyval.funcArgs->addArg(yyvsp[-3].var);} CATCH;}
    break;


    }

/* Line 999 of yacc.c.  */
#line 5010 "FRParser.cpp"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* Return failure if at end of input.  */
      if (yychar == YYEOF)
        {
	  /* Pop the error token.  */
          YYPOPSTACK;
	  /* Pop the rest of the stack.  */
	  while (yyss < yyssp)
	    {
	      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
	      yydestruct (yystos[*yyssp], yyvsp);
	      YYPOPSTACK;
	    }
	  YYABORT;
        }

      YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
      yydestruct (yytoken, &yylval);
      yychar = YYEMPTY;

    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*----------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action.  |
`----------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      yyvsp--;
      yystate = *--yyssp;

      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 1370 "FRParser.ypp"


/** Invoked when FRParser Class is 'loaded' */
void FRParser::load(CClassSP& clazz){
    REGISTER(FRParser, clazz);
    SUPERCLASS(FR::RValueBase);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(expression, "expression to parse");
    clazz->setPublic(); // make visible to EAS/spreadsheet
    Addin::registerConstructor("FR_EXPRESSION",
                               Addin::FLEX_PAYOFF,
                               "Creates a FLEX_SPI Expression",
                               TYPE);
    // register our [grammar specific] functions
    FRFunction::registerParserSpecificFunction(
        MAX, "MAX", "Computes max of two doubles/ints", 2);
    FRFunction::registerParserSpecificFunction(
        MIN, "MIN", "Computes min of two doubles/ints", 2);
    FRFunction::registerParserSpecificFunction(
        SUBARRAY, "SUBARRAY", "Returns a subset of the supplied array "
        "starting at specified index and of specified length", 3);
}

/** Returns bison type of expression and populates [static] yylval value
    with correct type */
static int expTypeToParserType(
    FRIfaces::VarType    varType,
    FRIfaces::IRValue*   rValue,
    YYSTYPE&             yylval){
    switch (varType){
    case FRIfaces::doubleType:
        yylval.dExp = &dynamic_cast<FRIfaces::IRValueDouble&>(*rValue);
        return DEXP;
    case FRIfaces::doubleArrayType:
        yylval.dExpArray = 
            &dynamic_cast<FRIfaces::IRValueDoubleArray&>(*rValue);
        return DEXPARRAY;
    case FRIfaces::intType:
        yylval.iExp = &dynamic_cast<FRIfaces::IRValueInt&>(*rValue);
        return IEXP;
    case FRIfaces::intArrayType:
        yylval.iExpArray = 
            &dynamic_cast<FRIfaces::IRValueIntArray&>(*rValue);
        return IEXPARRAY;
    case FRIfaces::boolType:
        yylval.bExp = &dynamic_cast<FRIfaces::IRValueBool&>(*rValue);
        return BEXP;
    case FRIfaces::boolArrayType:
        yylval.bExpArray = 
            &dynamic_cast<FRIfaces::IRValueBoolArray&>(*rValue);
        return BEXPARRAY;
    case FRIfaces::dateType:
        yylval.dtExp = &dynamic_cast<FRIfaces::IRValueDate&>(*rValue);
        return DTEXP;
    case FRIfaces::scheduleType:
        yylval.schedExp = &dynamic_cast<FRIfaces::IRValueSchedule&>(*rValue);
        return SCHEDEXP;
    case FRIfaces::tabulatedFuncType:
        yylval.tabFuncExp = 
            &dynamic_cast<FRIfaces::IRValueTabulatedFunc&>(*rValue);
        return TABFUNCEXP;
    default:
        throw ModelException("expTypeToParserType", "Unrecognised "
                             "variable type");
    }
}

/** Maps FRIfaces::VarType to parser's variable type
    (ie DVAR IVAR BVAR DTVAR SCHEDVAR TABFUNCVAR) */
static int varTypeToParserVarType(
    FRIfaces::VarType varType){
    switch (varType){
    case FRIfaces::doubleType:
        return DVAR;
    case FRIfaces::doubleArrayType:
        return DVARARRAY;
    case FRIfaces::intType:
        return IVAR;
    case FRIfaces::intArrayType:
        return IVARARRAY;
    case FRIfaces::boolType:
        return BVAR;
    case FRIfaces::boolArrayType:
        return BVARARRAY;
    case FRIfaces::dateType:
        return DTVAR;
    case FRIfaces::scheduleType:
        //return SCHEDVAR;
    case FRIfaces::tabulatedFuncType:
        //return TABFUNCVAR;
        throw ModelException("varTypeToParserVarType", "Unsupported "
                             "variable type");
    default:
        throw ModelException("varTypeToParserVarType", "Unrecognised "
                             "variable type");
    }
}

/** Maps FRIfaces::VarType to parser's function type
    (ie DFUNC DTFUNC) */
static int funcTypeToParserFuncType(
    FRIfaces::VarType varType){
    switch (varType){
    case FRIfaces::doubleType:
        return DFUNC;
    case FRIfaces::doubleArrayType:
        return DARRAYFUNC;
    case FRIfaces::dateType:
        return DTFUNC;
    case FRIfaces::intType:
        return IFUNC;
    case FRIfaces::intArrayType:
        return IARRAYFUNC;
    case FRIfaces::boolType:
        return BFUNC;
    case FRIfaces::boolArrayType:
        return BARRAYFUNC;
        // have to add these in possibly ...
    case FRIfaces::scheduleType:
    case FRIfaces::tabulatedFuncType:
    default:
        throw ModelException("funcTypeToParserFuncType", "Unrecognised "
                             "function type");
    }
}


/* This is the string parsing function - it just chops the string into
   different symbols/values and returns them to the grammatical parser */
static int yylex(void* lvalpParam, void* parserParam){
    YYSTYPE*  lvalp = (YYSTYPE*)lvalpParam;
    FRParser::Params& myParams = *(FRParser::Params*)parserParam;
    try{
        myParams.prevPos = myParams.thePos; // save previous location
        const char*& pos = myParams.thePos;//we will modify the vale of thePos
        YYSTYPE& yylval = *lvalp; // for ease (looks more familiar)
        /* Ignore whitespace, get first nonwhite character.  */
        char   c;
        while ((c = *pos) == ' ' || c == '\t'){
            pos++;
        }
    
        /* Char starts a number => parse the number.         */
        if (c == '.' || isdigit (c)) {
            /* see if it's an int */
            const char* tmpPos = pos;
            while (isdigit(*tmpPos)){
                tmpPos++;
            }
            if (*tmpPos != '.'){ 
                // must be an int
                char* endPos;
                int iVal = strtol(pos, &endPos, 10 /* base 10 */);
                pos = endPos;
                // now convert to RValueInt
                yylval.iExp = new FR::RConstInt(iVal);
                myParams.theCtrl->store(yylval.iExp); // memory management
                return IEXP;
            } else if (isdigit(c) || isdigit(*(pos+1))){
                // so either of form .55, or 55. or 55.55
                char* endPos;
                double dVal = strtod(pos, &endPos);
                pos = endPos;
                // now convert to RValueInt
                yylval.dExp = new FR::RConstDouble(dVal);
                myParams.theCtrl->store(yylval.dExp); // memory management
                return DEXP;
            } else {
                // using '.' as an operator
                pos++; // move onto next character
                return c;
            }
        }
        /* Char starts an identifier => read the name.       */
        if (isalpha(c) || c == '_'){
            const char* varBegin = pos;
            do{
                pos++; /* Get another character. */
                c = *pos;
            } while (c != '\0' && (isalnum(c) || c == '_'));
            const char* varEnd = pos;
            // construct string holding variable name
            string varName(varBegin, varEnd - varBegin);
            // test for true/false
            bool isTrue = varName == "true";
            if (isTrue || varName == "false"){
                yylval.bExp = new FR::RConstBool(isTrue);
                myParams.theCtrl->store(yylval.bExp); // memory management
                return BEXP;
            } else {
                // what is the variable/function followed by ? eg '(','[' etc
                while ((*varEnd == ' ' || *varEnd == '\t') && *varEnd != '\0'){
                    varEnd++;
                }
                if (*varEnd == '('){
                    // see if it's a function
                    const FRFunction* func = FRFunction::find(varName);
                    if (!func){
                        // wrap our FRParseException inside a ModelException
                        // (since ourError is not a pointer)
                        myParams.ourError = ModelException(
                            FRParseException("Undefined function '"+
                                             varName+"'"));
                        return UNDEFINED_SYMBOL;
                    }
                    yylval.func = func;
                    int parserType = yylval.func->getGrammarType();
                    if (yylval.func->isGeneric()){
                        // map return value
                        parserType = funcTypeToParserFuncType(
                            (FRIfaces::VarType)parserType);
                    }
                    return parserType;
                }
                // look up variable name using controller
                FRIfaces::ILValueExpression* var = 
                    myParams.theCtrl->getLValueExpression(varName);
                if (!var){
                    // wrap our FRParseException inside a ModelException
                    // (since ourError is not a pointer)
                    myParams.ourError = ModelException(
                        FRParseException("Undefined variable '"+varName+"'"));
                    return UNDEFINED_SYMBOL;
                }
                FRIfaces::VarType varType =
                    var->getType(); // one of doubleType/intType etc
                // a variable followed by a '[' means the value at a
                // specific simulation date unless the variable is
                // simulation date independent
                if (*varEnd == '[' && !var->isSDI()){
                    // return a new type which is the LValueExpression
                    // - however need to determine type of variable
                    yylval.var = var;
                    int parserType = varTypeToParserVarType(varType);
                    return parserType;
                } else {
                    FRIfaces::IRValue* rValue = 
                        var->getRValue(myParams.theCtrl);
                    // then map varType into DEXP/IEXP etc
                    int parserType = expTypeToParserType(varType, rValue, 
                                                         yylval);
                    return parserType;
                }
            }
        }
        pos++;

        /* A very manual way of checking for other tokens - we don't actually
           need this at the moment as we have no '=' operator so in the grammar
           rules we could just look for '=' '=' but I've left it here as
           reminder of what to do if we did have a '=' operator */
        if (c == '=' && *pos == '='){
            pos++;
            return EQUALS;
        }
        if (c == '!' && *pos == '='){
            pos++;
            return NEQUALS;
        }
        if (c == '&' && *pos == '&'){
            pos++;
            return AND;
        }
        if (c == '|' && *pos == '|'){
            pos++;
            return OR;
        }
        /* Any other character is a token by itself.        */
        return c;
    } catch (exception& e){
        myParams.ourError = ModelException(e);
        return INTERNAL_ERROR;
    }
}

void yyerror (const char *s)  /* Called by yyparse on error */{
    static const string undefinedString("$undefined");
    try{
        errorMsg = string(s); // record error and then fail
        // improve error message by replacing '$undefined'
        unsigned int undefined = errorMsg.find(undefinedString);
        if (undefined != string::npos){
            // if found then swap
            errorMsg.replace(undefined, undefinedString.size(), "symbol");
        }
        // improve error message by replacing '$' with 'end of expression'
        unsigned int dollarPos = errorMsg.find('$');
        if (dollarPos != string::npos){
            // if found then swap
            errorMsg.replace(dollarPos, 1, "end of expression");
        }
        if (!errorMsg.empty()){
            errorMsg[0] = toupper(errorMsg[0]); // capitalize first letter
        }
    } catch (exception& e){
        errorMsg = e.what(); // best we can do
    }
}

// external symbol to force this file to be linked in
bool loadFRParser() {
    return (FRParser::TYPE != 0);
}


DRLIB_END_NAMESPACE


