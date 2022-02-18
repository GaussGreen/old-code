/* ==========================================================================
   FILE_NAME:	srt_f_rtdesc.c

   PURPOSE:     function to deal with rates descriptions...
   ========================================================================== */

#include "math.h"
#include "srt_h_all.h"
#include "stdarg.h"
#include "utallhdr.h"

#define MAXNFP 250

/*------------for accessing interest rate model information------------*/
void srt_f_rtsetyp(SrtRtFnc rt, void* ptr)
{
    rt->yp = ptr;
}
void* srt_f_rtgetyp(SrtRtFnc rt)
{
    return rt->yp;
}
void srt_f_rtsettm(SrtRtFnc rt, double* tm)
{
    rt->t = tm;
}
double* srt_f_rtgettm(SrtRtFnc rt)
{
    return rt->t;
}
void srt_f_rtsetdf(SrtRtFnc rt, double* df)
{
    rt->df = df;
}
double* srt_f_rtgetdf(SrtRtFnc rt)
{
    return rt->df;
}
void srt_f_rtsetdt(SrtRtFnc rt, Ddate* dt)
{
    rt->date = dt;
}
Ddate* srt_f_rtgetdt(SrtRtFnc rt)
{
    return rt->date;
}
void srt_f_rtsetcvg(SrtRtFnc rt, Ddate* cvg)
{
    rt->cvg = cvg;
}
Ddate* srt_f_rtgetcvg(SrtRtFnc rt)
{
    return rt->cvg;
}
int srt_f_rtgetlen(SrtRtFnc rt)
{
    return rt->len;
}

/*-----------------------for translating between strings and enums-----*/
static char* srtrtfncnames[] = {"SWAP", "LVL", "IMMSWP", "IMMLVL", "DF", "FRA"};

static Err srt_f_rttrn(SrtRtFncType type, char** name)
{
    if ((int)type < 0 || (int)type >= LASTSRTRTFNCTYPE)
    {
        return serror("Bad SrtRtFncType");
    }

    *name = srtrtfncnames[(int)type];

    return NULL;
}

static Err srt_f_rtintp(SrtRtFncType* type, char* name)
{
    int i;

    for (i = 0; i < (int)LASTSRTRTFNCTYPE; i++)
    {
        if (!strcmp(name, srtrtfncnames[i]))
        {
            *type = (SrtRtFncType)i;
            return NULL;
        }
    }

    return serror("Don't know %s\n", name);
}

/*-------------------------functions operating on  those structures--*/

/*
        allocate a SrtRtFnc
*/
static SrtRtFnc srt_f_rtalloc(void)
{
    SrtRtFnc rt;

    rt = srt_calloc(1, sizeof(struct _srt_rtfnc_desc));

    return rt;
}

/*
        free a SrtRtFnc
*/
void srt_f_rtfre(SrtRtFnc rt)
{
    if (!rt)
        return;
    if (rt->date)
        srt_free(rt->date);
    if (rt->cvg)
        srt_free(rt->cvg);
    if (rt->yp)
        srt_free(rt->yp);
    if (rt->t)
        srt_free(rt->t);
    if (rt->df)
        srt_free(rt->df);
    srt_free(rt);
}

/*-----------------input validation-----------------------------*/
/*
        Error check the arguments describing a swap rate, then populate a
        SrtSwapDP
*/
static Err srt_f_sefb_to_SwapDP(
    Date start, Date end_or_nfp, char* freq, char* basis, SrtSwapDP* sdp)
{
    Err err;

    err = srt_f_gen_test_date(start);
    if (err)
        return err;
    sdp->start = start;

    if (end_or_nfp > start)
    {
        err = srt_f_gen_test_date(end_or_nfp);
        if (err)
            return err;
        sdp->end       = end_or_nfp;
        sdp->direction = BKWD;
    }
    else if (0 < (long)end_or_nfp && (long)end_or_nfp < 361)
    {
        sdp->nfp               = end_or_nfp;
        sdp->direction         = FWD;
        sdp->first_full_fixing = sdp->start;
    }
    else
        return serror("bad end_or_nfp srt_f_rtdesc");

    err = interp_compounding(freq, &sdp->compd);
    if (err)
        return err;

    err = interp_basis(basis, &sdp->basis_code);
    if (err)
        return err;

    return NULL;
}

/*------------------------evaluation------------------------*/

/*
        evaluate rate function assuming discount factors have been set
*/
Err srt_f_rtevl(SrtRtFnc rt, double short_rate, double* answer)
{
    int    i;
    double x = 0;

    switch (rt->type)
    {
    case SRTRTFNCSWP:
    case SRTRTFNCIMMSWP:
    case SRTRTFNCFRA:
    {
        for (i = 1; i < rt->len; i++)
        {
            x += rt->df[i] * rt->cvg[i];
        }
        /*
                if rates are incredibly high, x will by zero; in this case
                replace the swap rate with the short_rate, to avoid numerical
                overflow
        */
        if (x < 1.0e-9)
        {
            x = short_rate;
        }
        else
        {
            x = (rt->df[0] - rt->df[rt->len - 1]) / x;
        }
        break;
    }
    case SRTRTFNCLVL:
    case SRTRTFNCIMMLVL:
    {
        for (i = 1; i < rt->len; i++)
        {
            x += rt->df[i] * rt->cvg[i];
        }
        break;
    }
    case SRTRTFNCDF:
    {
        x = rt->df[0];
        /*
                test for very high rates --->overflow
        */
        if (x < 1.0e-9)
        {
            x = 0.0;
        }
        else
        {
            x = rt->df[1] / x;
        }
        break;
    }
    default:
    {
        return serror("srt_f_rtevl bad input");
    }
    }
    if (answer)
    {
        *answer = x;
    }
    return NULL;
}

/*--------------functions to create descriptors of rate functions-------*/

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :srt_f_rtmkimmswp
  AUTHOR          :E.Auld
  DATE            : 25 Apr 1995
  DESCRIPTION     :given input like the swap tools function
    @immswp(...), creates a SrtRtFnc describing an IMM swap.  Note that
    the only difference in an IMM swap is the dates; note also the implicit
    assumption in all this that IMM swaps trade at par ==> this is the
    assumption made in swap tools when it says that the rates for IMM futures
    (which are 3 month LIBOR) are also the rates between the future set dates
    (when actually this is for a period that might differ several days from
    the period over which 3 month LIBOR applies.)
    Syntax: start is only used if startfut is NULL
            end is only used if endfut is NULL
    if end is used and is not a number of full periods, then it must be an
    IMM date. (Or there will be an error). Num periods is interpreted as
    the number of periods after the first fixing date >= start.
    FREQUENCY: if frequency is less than monthly, dates must fall on IMM
    quarterly months, e.g. {m,j,s,d}.


  MODIFIES        :*rt
  CALL            :

<%%END>---------------------------------------------------------------------*/

static SrtErr intpfutstr(char* futstr, SrtMonth* m, int* y, int* nm)
{
    char   buf[4];
    SrtErr err;
    int    len = strlen(futstr);
    if (len != 5 && len != 6)
        return serror("bad fut str %s", futstr);
    strupper(futstr);
    strncpy(buf, futstr, 3);
    buf[3] = '\0';
    err    = interp_month(buf, m);
    if (err)
        return err;
    *y = atoi(futstr + 3);
    if (*y > 80)
        *y = *y + 1900;
    else
        *y = *y + 2000;
    if (len == 5 || futstr[5] == 'Q')
    {
        if (*m != SRT_MAR && *m != SRT_JUN && *m != SRT_SEP && *m != SRT_DEC)
        {
            return serror("Quarterly futs must be in {m,j,s,d}: %s", futstr);
        }
        *nm = 3;
    }
    else if (futstr[5] == 'M')
        *nm = 1;
    else
        return serror("bad fut str %s", futstr);
    return NULL;
}

SrtErr srt_f_rtmkimmswp(
    Date      start,
    char*     startfutstr,
    Date      end_or_nfp,
    char*     endfutstr,
    char*     freqstr,
    char*     basisstr,
    SrtRtFnc* rt)
{
    double *       cvg, *t;
    Ddate*         swpdate;
    int            i, len, nm = -1, nm2 = -1, nfp = -1, y, y2, flg, diffm;
    Date           end = ((Date)-1), first_non_stub = ((Date)-1);
    SrtErr         err;
    SrtCompounding freq;
    SrtBasisCode   basis;
    SrtMonth       m, m2;

    *rt = NULL;

    err = interp_compounding(freqstr, &freq);
    if (err)
        return err;

    err = interp_basis(basisstr, &basis);
    if (err)
        return err;

    /*
     * check if start is a future date or a date
     */

    if (startfutstr)
    {
        err = intpfutstr(startfutstr, &m, &y, &nm);
        if (err)
            return err;
        start          = third_wednesday(m, y);
        first_non_stub = start;
    }
    else
    {
        err = srt_f_gen_test_date(start);
        if (err)
            return err;
    }

    /*
     * check if end is a future date or a date or an nfp
     */

    if (endfutstr)
    {
        err = intpfutstr(endfutstr, &m, &y, &nm2);
        if (err)
            return err;
        end = third_wednesday(m, y);
        if (nm > 0 && nm != nm2)
        {
            return serror("futs must be same kind %s, %s for IMM swap", startfutstr, endfutstr);
        }
    }
    else if (end_or_nfp > MAXNFP)
    {
        err = srt_f_gen_test_date(end_or_nfp);
        if (err)
            return err;
        end = end_or_nfp;
    }
    else if (end_or_nfp > 0 && end_or_nfp <= MAXNFP)
    {
        nfp = end_or_nfp;
    }
    else
    {
        return serror("bad end date %d for immswp", end_or_nfp);
    }
    /*
     * deal with case where we know start and nfp
     */
    if (nfp > 0)
    {
        m = month(start);
        y = year(start) + 1900;
        if (nm < 0 || freq < 12)
            nm = 3;
        if (nm == 3)
        {
            m              = m + 3 - m % 3 - 3 * (m % 3 == 0);
            first_non_stub = third_wednesday(m, y);
            if (first_non_stub < start)
            {
                m += 3;
                y += (m > 12);
                m -= 12 * (m > 12);
                first_non_stub = third_wednesday(m, y);
            }
        }
        else
        {
            first_non_stub = third_wednesday(m, y);
            if (first_non_stub < start)
            {
                m += 1;
                y += (m > 12);
                m -= 12 * (m > 12);
                first_non_stub = third_wednesday(m, y);
            }
        }
        m2 = (m + (12 / freq) * nfp) % 12;
        if (m2 == 0)
            m2 = 12;
        y2  = y + nfp / freq + (m2 < m);
        end = third_wednesday(m2, y2);
    }

    /*
     * deal with case where we know start and end
     */

    /*
     * check that if freq is less than monthly, end date falls on {m,j,s,d}
     */
    m2 = month(end);
    y2 = year(end) + 1900;
    if ((int)freq < 12)
    {
        if (m2 != SRT_MAR && m2 != SRT_JUN && m2 != SRT_SEP && m2 != SRT_DEC)
        {
            return serror("%s IMM swaps must end on {m,j,s,d}", freqstr);
        }
    }
    if (end != third_wednesday(m2, y2))
    {
        return serror("IMM swaps must end exactly on future dates.");
    }

    /*
     * get first_non_stub
     */
    if (nfp < 0)
    {
        m     = month(start);
        y     = year(start) + 1900;
        diffm = m2 - m + 12 * (y2 - y);
        nfp   = diffm / (12 / freq);
        m     = (m2 - nfp * (12 / freq)) % 12;
        if (m % 12 == 0)
            m = 12;
        while (m < 0)
            m += 12;
        y              = y2 - nfp / freq - (m2 < m);
        first_non_stub = third_wednesday(m, y);
        if (first_non_stub < start)
        {
            m += (12 / freq);
            y += (m > 12);
            m -= 12 * (m > 12);
            first_non_stub = third_wednesday(m, y);
            nfp--;
        }
    }

    if (first_non_stub > end || start >= end || start > first_non_stub || nfp < 0)
    {
        return serror("imm_swp bad dates");
    }

    /*
     * (at last) generate the dates
     */
    len     = nfp + 2 - (start == first_non_stub);
    *rt     = srt_f_rtalloc();
    t       = srt_calloc(len, sizeof(double));
    cvg     = srt_calloc(len, sizeof(double));
    swpdate = srt_calloc(len, sizeof(Ddate));

    if (!(*rt) || !t || !cvg || !swpdate)
        return serror("imm_swp alloc failed");

    (*rt)->type = SRTRTFNCSWP;
    (*rt)->len  = len;
    (*rt)->cvg  = cvg;
    (*rt)->date = swpdate;

    flg          = (first_non_stub > start);
    swpdate[0]   = start;
    swpdate[flg] = first_non_stub;
    m            = month(first_non_stub);
    y            = 1900 + year(first_non_stub);
    for (i = 0; i < nfp; i++)
    {
        m += 12 / freq;
        y += (m > 12);
        m -= 12 * (m > 12);
        swpdate[flg + i + 1] = third_wednesday(m, y);
        cvg[flg + i + 1]     = coverage(DTOL(swpdate[flg + i]), DTOL(swpdate[flg + i + 1]), basis);
    }
    if (flg)
    {
        cvg[flg] = coverage(DTOL(swpdate[flg - 1]), DTOL(swpdate[flg]), basis);
    }
    return NULL;
}

/*
        create a descriptor of a swap rate
*/
static Err srt_f_rtmkswp(Date start, Date end_or_nfp, char* freq, char* basis, SrtRtFnc* ret)
{
    SrtRtFnc    rt;
    SrtDateList dl;
    SwapDP      sdp;
    int         i;
    Err         err;

    err = srt_f_sefb_to_SwapDP(start, end_or_nfp, freq, basis, &sdp);
    if (err)
        return err;

    rt = srt_f_rtalloc();
    if (!rt)
        return serror("srt_f_rtalloc failed");
    dl       = SwapDP_to_DateList(&sdp, MODIFIED_SUCCEEDING);
    rt->date = srt_calloc(dl.len, sizeof(Ddate));
    rt->cvg  = srt_calloc(dl.len, sizeof(double));
    rt->df   = srt_calloc(dl.len, sizeof(double));
    rt->len  = dl.len;
    if (!rt->date || !rt->cvg || !rt->df)
    {
        return serror("alloc failure in srt_f_rtdesc");
    }
    for (i = 0; i < dl.len; i++)
    {
        rt->date[i] = dl.date[i];
    }
    for (i = 1; i < dl.len; i++)
    {
        rt->cvg[i] = coverage(dl.date[i - 1], dl.date[i], sdp.basis_code);
    }
    if (dl.date)
        srt_free(dl.date);
    rt->type = SRTRTFNCSWP;
    *ret     = rt;
    return NULL;
}

/*
        create a descriptor of an fra
*/
static Err srt_f_rtmkfra(Date start, Date end_or_nfp, char* basis, SrtRtFnc* ret)
{
    Err          err;
    SrtRtFnc     rt;
    SrtBasisCode basis_code;

    err = srt_f_gen_test_date(start);
    if (err)
        return err;

    if (end_or_nfp > start)
    {
        err = srt_f_gen_test_date(end_or_nfp);
        if (err)
            return err;
    }
    else if (0 < (long)end_or_nfp && (long)end_or_nfp < 361)
    {
        end_or_nfp = add_unit(start, (int)end_or_nfp, SRT_MONTH, MODIFIED_SUCCEEDING);
    }
    else
        return serror("bad end_or_nfp srt_f_rtdesc");

    err = interp_basis(basis, &basis_code);
    if (err)
        return err;

    rt = srt_f_rtalloc();
    if (!rt)
        return serror("failed to alloc");

    rt->date = srt_calloc(2, sizeof(Ddate));
    rt->cvg  = srt_calloc(2, sizeof(double));
    rt->df   = srt_calloc(2, sizeof(double));
    rt->len  = 2;

    if (!rt->date || !rt->cvg || !rt->df)
    {
        return serror("alloc failure in srt_f_rtdesc");
    }

    rt->date[0] = start;
    rt->date[1] = end_or_nfp;
    rt->cvg[1]  = coverage(start, end_or_nfp, basis_code);

    rt->type = SRTRTFNCFRA;
    *ret     = rt;
    return NULL;
}

SrtErr srt_f_rtmk(
    char*     rtname,
    Date      start,
    char*     stastr,
    Date      end_or_nfp,
    char*     endstr,
    char*     freq,
    char*     basis,
    SrtRtFnc* ret)
{
    SrtErr       err;
    SrtRtFncType type;

    err = srt_f_rtintp(&type, rtname);
    if (err)
        return err;

    switch (type)
    {
    case SRTRTFNCIMMSWP:
    case SRTRTFNCIMMLVL:
    {
        err = srt_f_rtmkimmswp(start, stastr, end_or_nfp, endstr, freq, basis, ret);
        break;
    }
    case SRTRTFNCSWP:
    case SRTRTFNCLVL:
    {
        err = srt_f_rtmkswp(start, end_or_nfp, freq, basis, ret);
        break;
    }
    case SRTRTFNCFRA:
    {
        err = srt_f_rtmkfra(start, end_or_nfp, basis, ret);
        break;
    }
    case SRTRTFNCDF:
    {
        err = srt_f_rtmkfra(start, end_or_nfp, "A5", ret);
        break;
    }
    default:
    {
        err = serror("problem in srt_f_rtmk");
        break;
    }
    }

    if (*ret)
        (*ret)->type = type;
    return err;
}

/*--------------debugging/testing---------------------------------*/
/* dd/mm/yy */
/* 01234567 */
static char* dtstr(Date dt)
{
    static char buf[14];
    SrtMonth    m;
    int         d, y;
    y = year(dt);
    if (y > 99)
        y -= 100;
    m = month(dt);
    d = day(dt);
    sprintf(buf, "%2d/%2d/%d", d, m, y);
    if (d < 10)
        buf[0] = '0';
    if (m < 10)
        buf[3] = '0';
    if (y < 10)
        buf[6] = '0';
    if (y == 10)
        buf[7] = '0';
    buf[8] = '\0';
    return buf;
}

void srt_f_rtprint(FILE* out, SrtRtFnc rt)
{
    int   i;
    char* name = "bad rt type";
    if (!out)
    {
        out = stdout;
    }
    if (!rt)
    {
        fprintf(out, "NULL SrtRtFnc\n");
        return;
    }

    srt_f_rttrn(rt->type, &name);
    fprintf(out, "type: %s\n", name);

    for (i = 0; i < rt->len; i++)
    {
        if (rt->date)
        {
            fprintf(out, "%s %.0lf  ", dtstr(DTOL(rt->date[i])), rt->date[i]);
        }
        if (rt->cvg)
        {
            fprintf(out, "%.6lf  ", rt->cvg[i]);
        }
        if (rt->df)
        {
            fprintf(out, "%.6lf  ", rt->df[i]);
        }
        if (rt->t)
        {
            fprintf(out, "%.6lf  ", rt->t[i]);
        }
        fprintf(out, "\n");
    }
    if (rt->yp)
    {
        fprintf(out, "YTatt_params set\n\n");
    }
    else
    {
        fprintf(out, "YTatt_params not set\n\n");
    }
}

/*main()
{
#define HANDLE {if(err){printf("Err: %s\n",(char*)err); exit(0);}}

        SrtRtFnc rt;
        Date s;
        Err err;

        s = date(11,11,1994);
        err = srt_f_rtmk("SWAP",s,NULL,10,NULL,"SA","A5",&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);

        err = srt_f_rtmk("LVL",s,NULL,10,NULL,"A","BB",&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);

        err = srt_f_rtmk("FRA",s,NULL,6,NULL,NULL,"MM",&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);

        err = srt_f_rtmk("DF",s,NULL,12,NULL,NULL,NULL,&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);

        err = srt_f_rtmk("IMMSWP",s,NULL,10,NULL,"SA","A5",&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);

        err = srt_f_rtmk("IMMSWP",0,"DEC95Q",10,NULL,"SA","A5",&rt);
        HANDLE;
        srt_f_rtprint(NULL,rt);
        srt_f_rtfre(rt);
}
*/
