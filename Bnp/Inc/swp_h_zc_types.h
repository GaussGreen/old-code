/* ===========================================================================

         FILENAME:     swp_h_zc_types.h

     PURPOSE:      Some structures used in stripping

   =========================================================================== */

#ifndef SWP_H_ZC_TYPES_H
#define SWP_H_ZC_TYPES_H

enum generate_messages
{
    GENERATE_SWAP,
    GENERATE_BOND,
    GENERATE_CAP,
    COMPUTE_DISC_FACTOR,
    COMPUTE_ZERO_RATE,
    COMPUTE_FRA,
    COMPUTE_FWD_DISC_FACTOR,
    COMPUTE_FWD_ZERO_RATE
};

/***
        structure underlying YC_obj, ZC_obj, Arg, ...
***/

typedef struct
{
    int      Dwidth;
    int*     Dlengths;
    double** Dfields;

    int    Iwidth;
    int*   Ilengths;
    long** Ifields;
} Field_List;

Err Field_List_set_Ilength(Field_List* fl, int l);

Err Field_List_set_Dlength(Field_List* fl, int l);

Err Field_List_set_Dwidth(Field_List* fl, int w);

Err Field_List_set_Iwidth(Field_List* fl, int w);

void free_Field_List(Field_List* fl);

Field_List* new_Field_List(void);

/*************** for field lists *******************/

#define field_Dget(fptr, FIELDNAME, index) (fptr->Dfields[FIELDNAME][index])
#define field_Iget(fptr, FIELDNAME, index) (fptr->Ifields[FIELDNAME][index])

#define field_Dset(fptr, FIELDNAME, index, val) \
    {                                           \
        fptr->Dfields[FIELDNAME][index] = val;  \
    }

#define field_Iset(fptr, FIELDNAME, index, val) \
    {                                           \
        fptr->Ifields[FIELDNAME][index] = val;  \
    }

#define field_Dq(fl, FIELDNAME)                            \
    {                                                      \
        (0 <= FIELDNAME && FIELDNAME < fl->Dwidth ? 1 : 0) \
    }
#define field_Iq(fl, FIELDNAME)                            \
    {                                                      \
        (0 <= FIELDNAME && FIELDNAME < fl->Iwidth ? 1 : 0) \
    }

#define field_Dlength(fl, FIELDNAME) fl->Dlengths[FIELDNAME]
#define field_Ilength(fl, FIELDNAME) fl->Ilengths[FIELDNAME]

/* field list based types */

/******** Zc_obj --- hold zero curve ****/

enum ZC_DFields
{
    ZC_DATE,
    ZC_TIME,
    ZC_ZERO_RATE,
    ZC_DF,
    ZC_DWIDTH
};
enum ZC_IFields
{
    ZC_IWIDTH
};

typedef struct
{
    double      ZC_TODAY;
    double      ZC_SPOT;
    CcyCode     ZC_CURRENCY;
    int         ZC_INTERP_METHOD;
    int         ZC_NUMRATE;
    Field_List* data;
} ZC_Obj;

#define ZC_field_Dget(zcptr, FIELDNAME, index) field_Dget(zcptr->data, FIELDNAME, index)

#define ZC_field_Iget(zcptr, FIELDNAME, index) field_Iget(zcptr->data, FIELDNAME, index)

#define ZC_field_Dateget(zcptr, FIELDNAME, index) DTOL(field_Dget(zcptr->data, FIELDNAME, index))

#define ZC_field_Dateset(zcptr, FIELDNAME, index, val)          \
    {                                                           \
        field_Dset(zcptr->data, FIELDNAME, index, (double)val); \
    }

#define ZC_Dateget(zcptr, FIELDNAME) DTOL((zcptr->FIELDNAME))

#define ZC_Dateset(zcptr, FIELDNAME, val) \
    {                                     \
        (zcptr->FIELDNAME) = (double)val; \
    }

#define ZC_field_Dset(zcptr, FIELDNAME, index, val)     \
    {                                                   \
        field_Dset(zcptr->data, FIELDNAME, index, val); \
    }

#define ZC_field_Iset(zcptr, FIELDNAME, index, val)     \
    {                                                   \
        field_Dset(zcptr->data, FIELDNAME, index, val); \
    }

#define ZC_field_Dlength(zcptr, FIELDNAME) field_Dlength(zcptr->data, FIELDNAME)

#define ZC_field_Ilength(zcptr, FIELDNAME) field_Ilength(zcptr->data, FIELDNAME)

#define ZC_Dget(zcptr, FIELDNAME) (zcptr->FIELDNAME)

#define ZC_Iget(zcptr, FIELDNAME) (zcptr->FIELDNAME)

#define ZC_Dset(zcptr, FIELDNAME, val) \
    {                                  \
        (zcptr->FIELDNAME) = (val);    \
    }

#define ZC_Iset(zcptr, FIELDNAME, val) \
    {                                  \
        (zcptr->FIELDNAME) = (val);    \
    }

#define ZC_field_set_Dlength(zcptr, len) Field_List_set_Dlength(zcptr->data, len)

#define ZC_field_set_Ilength(zcptr, len) Field_List_set_Ilength(zcptr->data, len)

/******** YC_obj --- hold swap curve, futures rates & dates ****/

enum YC_DFields
{
    YC_SWAPEND_DATE,
    YC_SWAP_RATE,
    YC_SWAP_FIRST_FULL_FIXING,

    YC_CASH_DATE,
    YC_CASH_RATE,

    YC_M1_FUTSTARTDATE,
    YC_M1_FUTENDDATE,
    YC_M1_FUTSETTLEDATE,
    YC_M1_FUTPRICE,

    YC_M3_FUTSTARTDATE,
    YC_M3_FUTENDDATE,
    YC_M3_FUTSETTLEDATE,
    YC_M3_FUTPRICE,

    YC_TOY_RATE,
    YC_TOY_STARTDATE,
    YC_TOY_ENDDATE,
    YC_DWIDTH
};

enum YC_IFields
{
    YC_SWAP_NUM_FULL_PERIOD,
    YC_SWAP_COMPD,
    YC_SWAP_BASIS_CODE,
    YC_CASH_BASIS_CODE,
    YC_FUT1M_BASIS_CODE,
    YC_FUT3M_BASIS_CODE,
    YC_IWIDTH
};

typedef struct
{
    CcyCode YC_CURRENCY;
    Ddate   YC_TODAY;
    Ddate   YC_SPOT;
    int     YC_OVERNIGHT;
    Ddate   YC_OVERNIGHT_DATE;
    double  YC_OVERNIGHT_RATE;
    int     YC_TOMNEXT;
    Ddate   YC_TOMNEXT_DATE;
    double  YC_TOMNEXT_RATE;

    int         YC_NUMCASH;
    int         YC_NUMSWAP;
    int         YC_NUMFUT1M;
    int         YC_NUMFUT3M;
    int         YC_NUMTOY;
    Field_List* data;
} YC_Obj;

#define YC_field_Dget(ycptr, FIELDNAME, index) field_Dget(ycptr->data, FIELDNAME, index)
#define YC_field_Iget(ycptr, FIELDNAME, index) field_Iget(ycptr->data, FIELDNAME, index)

#define YC_field_Dset(ycptr, FIELDNAME, index, val)     \
    {                                                   \
        field_Dset(ycptr->data, FIELDNAME, index, val); \
    }
#define YC_field_Iset(ycptr, FIELDNAME, index, val)     \
    {                                                   \
        field_Iset(ycptr->data, FIELDNAME, index, val); \
    }

#define YC_field_Dateget(ycptr, FIELDNAME, index) DTOL(field_Dget(ycptr->data, FIELDNAME, index))
#define YC_field_Dateset(ycptr, FIELDNAME, index, val)          \
    {                                                           \
        field_Dset(ycptr->data, FIELDNAME, index, (double)val); \
    }

#define YC_Dateget(ycptr, FIELDNAME) DTOL((ycptr->FIELDNAME))
#define YC_Dateset(ycptr, FIELDNAME, val) \
    {                                     \
        (ycptr->FIELDNAME) = (double)val; \
    }

#define YC_field_Dlength(ycptr, FIELDNAME) field_Dlength(ycptr->data, FIELDNAME)
#define YC_field_Ilength(ycptr, FIELDNAME) field_Ilength(ycptr->data, FIELDNAME)

#define YC_Dget(ycptr, FIELDNAME) (ycptr->FIELDNAME)
#define YC_Iget(ycptr, FIELDNAME) (ycptr->FIELDNAME)

#define YC_Dset(ycptr, FIELDNAME, val) \
    {                                  \
        (ycptr->FIELDNAME) = (val);    \
    }
#define YC_Iset(ycptr, FIELDNAME, val) \
    {                                  \
        (ycptr->FIELDNAME) = (val);    \
    }

#define YC_field_set_Dlength(ycptr, len) Field_List_set_Dlength(ycptr->data, len)
#define YC_field_set_Ilength(ycptr, len) Field_List_set_Ilength(ycptr->data, len)

/****************** 	Leg --- one leg of a swap ***********/

enum Leg_DFields
{
    LEG_DATE,
    LEG_TIME_IN_YEARS,
    LEG_FIXING_DATE,
    LEG_CVG_START_LAG,
    LEG_CVG_END_LAG,
    LEG_PAYMENT,
    LEG_COVERAGE,
    LEG_DISC_FACTOR,
    LEG_STRIKE,
    LEG_VOLATILITY,
    LEG_SWAP_VOLATILITY,
    LEG_SWAP_FORWARD,
    LEG_FORWARD,
    LEG_CM_RATE,
    LEG_NOTIONAL,
    LEG_DWIDTH
};
enum Leg_IFields
{
    LEG_IWIDTH
};

typedef struct
{
    double         LEG_NOTIONAL;
    Ddate          LEG_TODAY;
    Ddate          LEG_START;
    Ddate          LEG_FIRST_FULL_FIXING;
    Ddate          LEG_END;
    double         LEG_STRIKE;
    double         LEG_SPREAD;
    Ddate          LEG_VALUE_DATE;
    double         LEG_INITIAL_NOT;
    double         LEG_FINAL_NOT;
    double         LEG_VALUE;
    double         LEG_FIRST_FIXING_PAYMENT;
    int            LEG_FIRST_FIXING_TYPE;
    int            LEG_TYPE;
    int            LEG_REC_PAY;
    int            LEG_NUM_PAYMENT;
    int            LEG_INDEX_START;
    int            LEG_INDEX_END;
    BasisCode      LEG_BASIS_CODE;
    SrtCompounding LEG_COMPOUNDING;
    int            LEG_NUM_FULL_PERIOD;
    Field_List*    data;
} Leg_Obj;

#define Leg_field_get(legptr, FIELDNAME, index) Leg_field_Dget(legptr, FIELDNAME, index)

#define Leg_field_set(legptr, FIELDNAME, index, val) Leg_field_Dset(legptr, FIELDNAME, index, val)

#define Leg_field_Dget(legptr, FIELDNAME, index) field_Dget(legptr->data, FIELDNAME, index)
#define Leg_field_Iget(legptr, FIELDNAME, index) field_Iget(legptr->data, FIELDNAME, index)

#define Leg_field_Dset(legptr, FIELDNAME, index, val)    \
    {                                                    \
        field_Dset(legptr->data, FIELDNAME, index, val); \
    }
#define Leg_field_Iset(legptr, FIELDNAME, index, val)    \
    {                                                    \
        field_Iset(legptr->data, FIELDNAME, index, val); \
    }

#define Leg_field_Dlength(legptr, FIELDNAME) field_Dlength(legptr->data, FIELDNAME)
#define Leg_field_Ilength(legptr, FIELDNAME) field_Ilength(legptr->data, FIELDNAME)

#define Leg_get(legptr, FIELDNAME) (legptr->FIELDNAME)
#define Leg_Dget(legptr, FIELDNAME) (legptr->FIELDNAME)
#define Leg_Iget(legptr, FIELDNAME) (legptr->FIELDNAME)

#define Leg_set(legptr, FIELDNAME, val) \
    {                                   \
        (legptr->FIELDNAME) = (val);    \
    }
#define Leg_Dset(legptr, FIELDNAME, val) \
    {                                    \
        (legptr->FIELDNAME) = (val);     \
    }
#define Leg_Iset(legptr, FIELDNAME, val) \
    {                                    \
        (legptr->FIELDNAME) = (val);     \
    }

#define Leg_field_set_Dlength(legptr, len) Field_List_set_Dlength(legptr->data, len)
#define Leg_field_set_Ilength(legptr, len) Field_List_set_Ilength(legptr->data, len)

#define Leg_field_Dateget(legptr, FIELDNAME, index) DTOL(field_Dget(legptr->data, FIELDNAME, index))
#define Leg_field_Dateset(legptr, FIELDNAME, index, val)         \
    {                                                            \
        field_Dset(legptr->data, FIELDNAME, index, (double)val); \
    }
#define Leg_Dateget(legptr, FIELDNAME) DTOL((legptr->FIELDNAME))
#define Leg_Dateset(legptr, FIELDNAME, val) \
    {                                       \
        (legptr->FIELDNAME) = (double)val;  \
    }

/*********** Arg --- for passing information ***************/
/******************** initialise *******************************/

typedef double (*Access_Func)();

typedef struct
{
    double         ARG_PV;
    double         ARG_NOTIONAL;
    Ddate          ARG_TODAY;
    Ddate          ARG_SPOT;
    Ddate          ARG_START;
    Ddate          ARG_END;
    Ddate          ARG_FIRST_FULL_FIXING;
    double         ARG_STRIKE;
    double         ARG_SPREAD;
    double         ARG_BOND_STRIKE;
    double         ARG_VOLATILITY;
    double         ARG_INITIAL_NOT;
    double         ARG_FINAL_NOT;
    Ddate          ARG_VALUE_DATE;
    double         ARG_FIRST_FIXING_PAYMENT;
    int            ARG_TYPE;
    int            ARG_REC_PAY;
    int            ARG_INDEX_START;
    int            ARG_INDEX_END;
    BasisCode      ARG_BASIS_CODE;
    SrtCompounding ARG_COMPD;
    int            ARG_NUM_DATES;
    int            ARG_NUM_FULL_PERIOD;
    int            ARG_INFO_TYPE;
    int            ARG_BID_OFF;
    int            ARG_FIRST_FIXING_TYPE;
    SwapDateDir    ARG_DATE_DIR;
    Access_Func    ARG_DATE;
    Access_Func    ARG_SET_FUNC1;
    Access_Func    ARG_SET_FUNC2;
    Access_Func    ARG_GET_FUNC1;
    Access_Func    ARG_GET_FUNC2;
    void*          ARG_TARGET;
    void*          ARG_RESULT_TARGET;
} Arg_Obj;

#define Arg_get(aptr, FIELDNAME) (aptr->FIELDNAME)
#define Arg_Dget(aptr, FIELDNAME) (aptr->FIELDNAME)
#define Arg_Iget(aptr, FIELDNAME) (aptr->FIELDNAME)
#define Arg_Dateget(aptr, FIELDNAME) DTOL(aptr->FIELDNAME)

#define Arg_set(aptr, FIELDNAME, val) (aptr->FIELDNAME = (val))
#define Arg_Dset(aptr, FIELDNAME, val) (aptr->FIELDNAME = (double)(val))
#define Arg_Iset(aptr, FIELDNAME, val) (aptr->FIELDNAME = (int)(val))
#define Arg_Dateset(aptr, FIELDNAME, val) (aptr->FIELDNAME = (double)(val))

#define Arg_set_func(aptr, FIELDNAME, pointer_func) \
    {                                               \
        (aptr)->FIELDNAME = pointer_func;           \
    }

#define Arg_field_get(aptr, mess, index) (*(aptr->mess))(aptr->ARG_TARGET, index)
#define Arg_field_set(aptr, mess, index, val)              \
    {                                                      \
        if (aptr->mess)                                    \
            (*(aptr->mess))(aptr->ARG_TARGET, index, val); \
    }

/****************** Swap ******************************/

typedef struct
{
    int      SWAP_TYPE;
    int      SWAP_REC_PAY;
    int      SWAP_NUM_LEG;
    int      SWAP_INFO_TYPE;
    int      SWAP_BID_OFF;
    double   SWAP_NOTIONAL;
    double   SWAP_VALUE;
    Ddate    SWAP_VALUE_DATE;
    double   SWAP_INITIAL_NOT;
    double   SWAP_FINAL_NOT;
    double   SWAP_VOLATILITY;
    double   SWAP_STRIKE;
    Leg_Obj* SWAP_LEG[10];      /* fix this later !!!! */
    double*  SWAP_VOLATILITIES; /* fix this later !!!! */
} Swap_Obj;

#define Swap_get(aptr, FIELDNAME) (aptr->FIELDNAME)
#define Swap_Dget(aptr, FIELDNAME) (aptr->FIELDNAME)
#define Swap_Iget(aptr, FIELDNAME) (aptr->FIELDNAME)

#define Swap_set(aptr, FIELDNAME, val) (aptr->FIELDNAME = (val))
#define Swap_Dset(aptr, FIELDNAME, val) (aptr->FIELDNAME = (val))
#define Swap_Iset(aptr, FIELDNAME, val) (aptr->FIELDNAME = (val))

#define Swap_field_get(sptr, FIELDNAME, index) (sptr->FIELDNAME[index])
#define Swap_field_set(sptr, FIELDNAME, index, val) \
    {                                               \
        sptr->FIELDNAME[index] = val;               \
    }

/*********** 	allocation	*******************/

ZC_Obj*   new_ZC_Obj();
YC_Obj*   new_YC_Obj();
Leg_Obj*  new_Leg();
Arg_Obj*  new_Arg();
Swap_Obj* new_Swap();

void free_YC_Obj(YC_Obj* yco);
void free_ZC_Obj(ZC_Obj* zco);
void free_Leg(Leg_Obj* leg);
void free_Arg(Arg_Obj* ga);
void free_Swap(Swap_Obj* s);

void leg_allocate(Leg_Obj** leg, int pay_size);
Err  init_Arg(Arg_Obj** ga, Message m);

#endif
