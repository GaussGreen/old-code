#include "grf_h_all.h"

/* functions to read grfn deals, other information from an ascii file */

#define BUFLEN 1000

#define GETLNG(LINE, T)                               \
    {                                                 \
        String TOK2, TOK = strtok(LINE, GRFSEPSTR);   \
        if (!TOK)                                     \
            return serror("failed to read long");     \
        else                                          \
        {                                             \
            T = strtol(TOK, &TOK2, 10);               \
            if (strlen(TOK2) != 0)                    \
                return serror("failed to read long"); \
        }                                             \
    }

#define GETDBL(LINE, T)                                 \
    {                                                   \
        String TOK2, TOK = strtok(LINE, GRFSEPSTR);     \
        if (!TOK)                                       \
            return serror("failed to read double");     \
        else                                            \
        {                                               \
            T = strtod(TOK, &TOK2);                     \
            if (strlen(TOK2) != 0)                      \
                return serror("failed to read double"); \
        }                                               \
    }

#define GETSTR(LINE, T)                          \
    {                                            \
        String TOK = strtok(LINE, GRFSEPSTR);    \
        if (!TOK)                                \
            return serror("failed to read str"); \
        else                                     \
        {                                        \
            T = TOK;                             \
        }                                        \
    }

#define GETBLN(LINE)                                  \
    {                                                 \
        String TOK = strtok(LINE, GRFSEPSTR);         \
        if (!TOK)                                     \
            return serror("failed to read bln");      \
        if (strlen(TOK) != 1 || TOK[0] != GRFBLNCHAR) \
            return serror("failed to read bln");      \
    }

/* returns 0 if EOF, otherwise reads in next line from *in into buf */
static long getline(String line, long buflen, FILE* in)
{
    if (fgets(line, buflen, in) == NULL)
        return 0;
    return strlen(line);
}

/* returns 0 if EOF, else return next line that does not start with
        comment char, with white space stripped, and capitalized */
static long getgoodline(String line, long buflen, FILE* in, String commentstr)
{
    long len = 0;

    while (1)
    {
        if (!getline(line, buflen, in))
            return 0;
        strupper(line);
        strip_white_space(line);
        len = strlen(line);
        if (len > 0 && strpbrk(line, commentstr) != line)
            break;
    }
    return len;
}

static Err find_line(String inputbuf, long buflen, String to_find, FILE* in)
{
    long len = getgoodline(inputbuf, buflen, in, GRFCOMSTR);
    if (!len)
        return serror("hit EOF looking for %s", to_find);
    if (strcmp(to_find, inputbuf))
        return serror("couldn't find %s", to_find);
    return NULL;
}

static Err read_event_row(String line, long rowlen)
{
    long   i;
    Date   d;
    char   c[2];
    String tok;
    String tok2;
    c[0] = GRFSEPCHAR;
    c[1] = '0';

    tok = strtok(line, c);
    if (!tok)
        return serror("failed to read event row");
    d = (Date)strtol(tok, &tok2, 10);
    if (strlen(tok2) != 0)
        return serror("evdate bad format");

    /* declare date to mad */

    for (i = 0; i < rowlen; i++)
    {
        tok = strtok(NULL, c);
        if (!tok)
            return serror("failed on event %d in row", i);
        if (strlen(tok) == 1 && tok[0] == GRFBLNCHAR)
        /* declare blank cell to mad */
        {
        }
        else
        /* declare tok (non blank cell) to mad */
        {
        }
    }
    return NULL;
}

#define GETGRFCELL(T)                                 \
    {                                                 \
        String TOK = strtok(NULL, GRFSEPSTR);         \
        if (!TOK)                                     \
            return serror("failed to read str");      \
        if (strlen(TOK) == 1 && TOK[0] == GRFBLNCHAR) \
        {                                             \
            T.type = GRFNBCELL;                       \
        }                                             \
        else                                          \
        {                                             \
            T.sval        = new_string(TOK);          \
            T.type        = GRFNSCELL;                \
            T.str_alloced = SRT_YES;                  \
        }                                             \
    }

Err grf_f_readcells(FILE* in, Date** eventdates, GrfnCell*** m, long* nr, long* nc)
{
    long        nrow, ncol, i, j, len;
    Err         err;
    static char inputbuf[BUFLEN];
    GrfnCell**  gc;
    Date*       dates;

    err = find_line(inputbuf, BUFLEN, GRFBEGEV, in);
    if (err)
        return err;

    len = getgoodline(inputbuf, BUFLEN, in, GRFCOMSTR);
    if (!len)
        return serror("premature EOF");

    GETLNG(inputbuf, nrow);
    GETLNG(NULL, ncol);

    if (nrow < 0 || ncol < 0 || nrow >= GRFMAXROW || ncol >= GRFMAXCOL)
        return serror("bad event dimensions");

    gc    = GrfnCellmatrix(nrow, ncol, 0);
    dates = srt_calloc(nrow, sizeof(Date));

    for (i = 0; i < nrow; i++)
    {
        len = getgoodline(inputbuf, BUFLEN, in, GRFCOMSTR);
        if (!len)
            return serror("premature EOF");
        GETLNG(inputbuf, dates[i]);
        for (j = 0; j < ncol; j++)
        {
            GETGRFCELL(gc[i][j]);
        }
    }

    err = find_line(inputbuf, BUFLEN, GRFENDEV, in);
    if (err)
        return err;

    *eventdates = dates;
    *m          = gc;
    *nr         = nrow;
    *nc         = ncol;
    return NULL;
}

Err grf_f_readauxrng(FILE* in, long* auxwidthptr, long** auxlenptr, double*** auxptr)
{
    long        i, j, auxwidth, *auxlen;
    double**    aux;
    Err         err;
    static char inputbuf[BUFLEN];
    long        nrow = 0, len, toalloc = 0;
    String      buf;

    *auxptr      = NULL;
    *auxlenptr   = NULL;
    *auxwidthptr = 0;

    err = find_line(inputbuf, BUFLEN, GRFBEGAUX, in);
    if (err)
        return err;

    len = getgoodline(inputbuf, BUFLEN, in, GRFCOMSTR);
    if (!len)
        return serror("premature EOF");
    if (!strcmp(inputbuf, GRFENDAUX))
        return NULL;

    GETLNG(inputbuf, auxwidth);
    if (auxwidth < 1 || auxwidth > 100)
        return serror("bad auxdimensions");

    aux          = (double**)srt_calloc(auxwidth, sizeof(double*));
    auxlen       = (long*)srt_calloc(auxwidth, sizeof(long));
    *auxptr      = aux;
    *auxlenptr   = auxlen;
    *auxwidthptr = auxwidth;

    len = getgoodline(inputbuf, BUFLEN, in, GRFCOMSTR);
    if (!len)
        return serror("premature EOF");
    GETLNG(inputbuf, auxlen[(i = 0)]);
    if (auxlen[i] < 1 || auxlen[i] > 500)
        return serror("bad auxdimensions");
    toalloc += auxlen[i];
    nrow = IMAX(nrow, auxlen[i]);
    for (i = 1; i < auxwidth; i++)
    {
        GETLNG(NULL, auxlen[i]);
        if (auxlen[i] < 1 || auxlen[i] > 500)
            return serror("bad auxdimensions");
        toalloc += auxlen[i];
        nrow = IMAX(nrow, auxlen[i]);
    }

    aux[0] = (double*)srt_calloc(toalloc, sizeof(double));
    for (i = 1; i < auxwidth; i++)
    {
        aux[i] = aux[i - 1] + auxlen[i - 1];
    }

    for (i = 0; i < nrow; i++)
    {
        len = getgoodline(inputbuf, BUFLEN, in, GRFCOMSTR);
        if (!len)
            return serror("premature EOF reading aux");
        buf = inputbuf;
        for (j = 0; j < auxwidth; j++)
        {
            if (i < auxlen[j])
            {
                GETDBL(buf, aux[j][i]);
                buf = NULL;
            }
            else
            {
                GETBLN(buf);
                buf = NULL;
            }
        }
    }

    err = find_line(inputbuf, BUFLEN, GRFENDAUX, in);
    if (err)
        return err;

    return NULL;
}

Err grf_f_readgrfrng(FILE* in, long* numgrng, GrfnRng** grng)
{
    Err         err;
    static char inputbuf[BUFLEN];

    *numgrng = 0;
    *grng    = NULL;

    err = find_line(inputbuf, BUFLEN, GRFBEGGRN, in);
    if (err)
        return err;
    err = find_line(inputbuf, BUFLEN, GRFENDGRN, in);
    if (err)
        return err;

    return NULL;
}

Err grf_f_readgrfdl(FILE* in, GrfnDeal* gd)
{
    Date*      eventdates;
    long       nrow;
    long       ncol;
    GrfnCell** sprdsht;
    long       numgrng = 0;
    GrfnRng*   grng;
    long       auxwidth;
    long*      auxlen;
    double**   aux;
    Err        err;

    err = grf_f_readauxrng(in, &auxwidth, &auxlen, &aux);
    if (err)
        return err;
    err = grf_f_readcells(in, &eventdates, &sprdsht, &nrow, &ncol);
    if (err)
        return err;
    err = grf_f_readgrfrng(in, &numgrng, &grng);
    if (err)
        return err;

    grfn_copy_to_GrfnDeal(
        gd, nrow, ncol, sprdsht, eventdates, numgrng, grng, auxwidth, auxlen, aux);

    if (aux[0])
        srt_free(aux[0]);
    if (aux)
        srt_free(aux);
    if (auxlen)
        srt_free(auxlen);
    if (sprdsht)
    {
        grfn_free_GrfnCellmatrix(sprdsht, nrow, ncol);
    }
    if (eventdates)
        srt_free(eventdates);
    if (grng)
        srt_free(grng); /*FIX*/

    return NULL;
}
