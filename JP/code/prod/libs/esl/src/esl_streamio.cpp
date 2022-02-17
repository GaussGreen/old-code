#if defined(_MSC_VER)
#pragma once
#pragma warning(disable:4786)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include "esl_log.h"
#include "esl_date.h"
#include "esl_types.h"
#include "esl_stdinput.h"
#include "irx/irxflow.h"
#include "irx/zerocurve.h"

#include <sys/types.h>
#include <sys/stat.h>

void EslHeaderInput(std::string const& key, std::istream& is)
{
    static char const*  routine = "EslHeaderInput ";

    std::string token("<" + key + ">");
    char        buf[1024];

    // title or preamble
    is.getline(buf, sizeof(buf));

    // check if has module separator
    if (!is || strstr(buf, token.c_str()) == 0)
        throw EslException(routine) << "expected '" << token << "' not found";
}


void EslFooterInput(std::string const& key, std::istream& is)
{
    static char const*  routine = "EslFooterInput ";

    std::string token("</" + key + ">");
    char        buf[1024];

    bool rv = false;

    // read until the end of the section
    while (is) {
        is.getline(buf, sizeof(buf));
        if (strstr(buf, token.c_str())) {
            rv = true;
            break;
        }
    }
    if (!rv)
        throw EslException(routine) << "expected '" << token << "' not found";
}



static
void  TermInputOne(T_CURVE*           crv,  ///< structure of zero curve data
                   std::string const& tag,  ///< tag
                   std::istream&      is)   ///< input stream              
{
    static char const*  routine = "TermInputOne ";

    char    buf[1024];

    // title or preamble
    EslHeaderInput(tag, is);

    // value date
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%ld", &crv->Today) != 1) 
        throw EslException(routine) << "Could not read value date";

    crv->Today = IRDateFromYMDDate(crv->Today);

    // Today's date and spot days are not available
    crv->ValueDate = crv->Today;
    crv->SpotDays  = 0;
    
    // money market basis
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &crv->MMB) != 1) 
        throw EslException(routine) << "Could not read money market basis";
        
    // yield curve frequency
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%c", &crv->SwapFreq) != 1) 
        throw EslException(routine) << "Could not read benchmark swap frequency";
        
    // swap DCC
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%s", crv->SwapDCC) != 1) 
        throw EslException(routine) << "Could not read swap DCC";
        
    // number of zeros
    int cnt;
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &cnt) != 1) 
        throw EslException(routine) << "Could not read number of zeros";
        
    if (cnt > MAXNBDATE)
         throw EslException(routine) << "number of rates exceeds " << MAXNBDATE;

    // zero dates and rates
    is.getline(buf, sizeof(buf));

    long    dates[MAXNBDATE];
    double  zeros[MAXNBDATE];

    for (int i = 0; i < cnt; ++i)
    {
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%ld %lf", &dates[i], &zeros[i]) != 2) 
            throw EslException(routine) << "Could not read date and rate in line # " << i+1;
        
        dates[i] = IRDateFromYMDDate(dates[i]);
        zeros[i] /= 100.;
    }
        

    /* We require at least 100 years of zero curve for swaption vol bootstrapping */
    IRDate LastDate = Nxtmth(crv->ValueDate, 1200L, 1L);

    /* If the zero curve does not extend up to last date we add an extra point */
    if (dates[cnt-1] < LastDate)
    {
        cnt += 1;                       
        dates[cnt-1] = LastDate;
        zeros[cnt-1] = zeros[cnt-2];    /* Flat zero curve */
    }
                
    /* copy local curve object to the output */
#ifdef ESL_NEW_CURVE
    if (irxZeroCurveConstructFromRates(
                crv,
                crv->ValueDate,
                cnt,
                dates, 
                zeros,
                IRX_ANNUAL_RATE,
                IRX_ACT_365F) != SUCCESS)
        throw EslException(routine) << "irxZeroCurveConstructFromRates failed";
#else
    for (int k = 0; k < cnt; ++k)
    {
        crv->ZeroDate[k] = dates[k];
        crv->Zero[k]     = zeros[k];
    }
    crv->NbZero = cnt;
#endif

    if (Term_Check_W (crv) == FAILURE)                    
        throw EslException(routine) << "Curve check failed";

    EslFooterInput(tag, is);
}


//-----------------------------------------------------------------------------
// Read all three zero curves at once - stream version. Expects all three 
// curves in the stream in the right order - zero, discount and risky
//-----------------------------------------------------------------------------
void  EslTermInput(T_CURVE       crvs[3],  ///< array of zero curves
                   std::istream& is)       ///< input stream
{
    TermInputOne(&crvs[0], "zero",     is);
    TermInputOne(&crvs[1], "disczero", is);
    TermInputOne(&crvs[2], "riskzero", is);
}


//-----------------------------------------------------------------------------
// Read all three zero curves at once - file version.
//-----------------------------------------------------------------------------
void  EslTermInput(T_CURVE       crvs[3])  ///< array of zero curves
{
    std::ifstream z1("zero.dat");
    TermInputOne(&crvs[0], "zero", z1);

    std::ifstream z2("disczero.dat");
    TermInputOne(&crvs[1], "disczero", z2);

    struct stat sstr;
    char const* file = "disczero.dat";

    if (stat("riskzero.dat", &sstr) == 0)
       file = "riskzero.dat";

    std::ifstream z3(file);
    TermInputOne(&crvs[2], "riskzero", z3);
}







//-----------------------------------------------------------------------------
// Read base vol data - stream version.
//-----------------------------------------------------------------------------
static
void BaseVolInput(MKTVOL_DATA*  mktvol,  /**< Volatility data          */
                  std::istream& is)      /**< input stream             */
{
    static char const*  routine = "BaseVolInput ";

    char    buf[1024];

    // title or preamble
    EslHeaderInput("basevol", is);

    // base vol frequency
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%c", &mktvol->Freq) != 1) 
        throw EslException(routine) << "Could not read base vol frequency";

    // number of base vols
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &mktvol->NbVol) != 1) 
        throw EslException(routine) << "Could not read number of base vols";

    if (mktvol->NbVol > MAXNBDATE)
        throw EslException(routine) << "number of base vols exceeds" << MAXNBDATE;

    // base vol dates and rates
    is.getline(buf, sizeof(buf));

    int i;
    for (i = 0; i < mktvol->NbVol; i++)
    {
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%ld %lf", &mktvol->VolDate[i], &mktvol->Vol[i]) != 2)
            throw EslException(routine) << "Could not read base vol date and rate # " << i+1;

        mktvol->VolDate[i] = IRDateFromYMDDate(mktvol->VolDate[i]);
        mktvol->Vol[i] /= 100.;
    }


    /* Eliminate dates falling before base date */
    int j = 0;
    while (mktvol->VolDate[j] <= mktvol->BaseDate)
        j++;

    mktvol->NbVol -= j;

    for (i = 0; i < mktvol->NbVol; i++)
    {
        mktvol->VolDate[i] = mktvol->VolDate[i+j];
        mktvol->Vol[i]     = mktvol->Vol[i+j];
    }

    EslFooterInput("basevol", is);
}


//-----------------------------------------------------------------------------
// Read base vol data and throw it away - needed for single stream input
//-----------------------------------------------------------------------------
static
void BaseVolSkip(std::istream& is)      ///< input stream
{
    static char const*  routine = "BaseVolSkip ";

    char    buf[1024];

    // title or preamble
    EslHeaderInput("basevol", is);

    // base vol frequency
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    char freq;
    if (!is || sscanf(buf, "%c", &freq) != 1) 
        throw EslException(routine) << "Could not read base vol frequency";

    // number of base vols
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    int nbvol;
    if (!is || sscanf(buf, "%d", &nbvol) != 1) 
        throw EslException(routine) << "Could not read number of base vols";

    if (nbvol > MAXNBDATE)
        throw EslException(routine) << "number of base vols exceeds" << MAXNBDATE;

    // base vol dates and rates
    is.getline(buf, sizeof(buf));

    int i;
    for (i = 0; i < nbvol; i++)
    {
        is.getline(buf, sizeof(buf));
        long    date;
        double  vol;
        if (!is || sscanf(buf, "%ld %lf", &date, &vol) != 2)
            throw EslException(routine) << "Could not read base vol date and rate # " << i+1;
    }

    EslFooterInput("basevol", is);
}



//-----------------------------------------------------------------------------
// Read base vol data - file version.
//-----------------------------------------------------------------------------
/*
static
void BaseVolInput(MKTVOL_DATA*  mktvol)  ///< Volatility data 
{
    std::ifstream is("basevol.dat");
    BaseVolInput(mktvol, is);
}
*/



typedef std::vector<std::vector<double> >   matrix_t;

//-----------------------------------------------------------------------------
// Read swap vol data - stream version.
//-----------------------------------------------------------------------------
static
void SwapVolInput(std::vector<IRDate>&   Expiry,
                  std::vector<IRDate>&   FwdMat,
                  matrix_t&             Matrix,
                  std::vector<char>&    DoMoY,
                  std::istream&         is)
{
    static char const*  routine = "SwapVolInput ";

    char    buf[1024];

    // title or preamble
    EslHeaderInput("swapvol", is);

    // number of expiries
    int rows;
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &rows) != 1) 
        throw EslException(routine) << "Could not read number of rows";

    if (rows > MAXNBDATE)
        throw EslException(routine) << "number of expiries exceeds " << MAXNBDATE;

    // number of columns
    int cols;
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &cols) != 1) 
        throw EslException(routine) << "Could not read number of columns";

    if (cols > MAXNBDATE)
        throw EslException(routine) << "number of forward maturities exceeds " << MAXNBDATE;

    // resize the storage
    int i, j;
    Expiry.resize(rows);
    FwdMat.resize(cols);
    Matrix.resize(rows);
    for(i=0; i<rows; ++i)
        Matrix[i].resize(cols);
    DoMoY.resize(rows);

    // swaption matrix
    is.getline(buf, sizeof(buf));
    for (j = 0; j < cols; j++)
    {
        is >> FwdMat[j];
        if (!is)
            throw EslException(routine) << "Could not read forward maturity # " << j+1;
        
        FwdMat[j] *= 12;    /* Convert to months */
    }
  
  
    for (i = 0; i < rows; i++)
    {
        char token[64];

        is >> token;
        if (!is)
            throw EslException(routine) << "Could not read expiry # " << i+1;

        Expiry[i] = atol(token);

        if (strchr(token, 'D') != NULL)
        {       
            DoMoY[i] = 'D';
        }
        else if (strchr(token, 'Y') != NULL)
        {
            DoMoY[i] = 'Y';
        }
        else
        {
            DoMoY[i] = 'M';
        }   
    
        for (j = 0; j < cols; j++)
        {
            is >> token;
            if (!is)
                throw EslException(routine) << "Could not read volatility # " << i+1 << ',' << j+1;

            Matrix[i][j] = atof(token) / 100.;
        }
    }
    // read to the end of line
    is.getline(buf, sizeof(buf));

    EslFooterInput("swapvol", is);
}



//-----------------------------------------------------------------------------
// Read swap vol data and throw it away - needed for single stream input
//-----------------------------------------------------------------------------
static
void SwapVolSkip(std::istream& is)
{
    static char const*  routine = "SwapVolSkip ";

    char    buf[1024];

    // title or preamble
    EslHeaderInput("swapvol", is);

    // number of expiries
    int rows;
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &rows) != 1) 
        throw EslException(routine) << "Could not read number of rows";

    if (rows > MAXNBDATE)
        throw EslException(routine) << "number of expiries exceeds " << MAXNBDATE;

    // number of columns
    int cols;
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &cols) != 1) 
        throw EslException(routine) << "Could not read number of columns";

    if (cols > MAXNBDATE)
        throw EslException(routine) << "number of forward maturities exceeds " << MAXNBDATE;

    // swaption matrix
    int i, j;
    is.getline(buf, sizeof(buf));
    for (j = 0; j < cols; j++)
    {
        long fwdmat;

        is >> fwdmat;
        if (!is)
            throw EslException(routine) << "Could not read forward maturity # " << j+1;
    }
  
    for (i = 0; i < rows; i++)
    {
        char token[64];

        is >> token;
        if (!is)
            throw EslException(routine) << "Could not read expiry # " << i+1;

        for (j = 0; j < cols; j++)
        {
            is >> token;
            if (!is)
                throw EslException(routine) << "Could not read volatility # " << i+1 << ',' << j+1;
        }
    }
    // read to the end of line
    is.getline(buf, sizeof(buf));

    EslFooterInput("swapvol", is);
}




//-----------------------------------------------------------------------------
// Read swap vol data - file version.
//-----------------------------------------------------------------------------
/*
static
void SwapVolInput(std::vector<IRDate>&   Expiry,
                  std::vector<IRDate>&   FwdMat,
                  matrix_t&             Matrix,
                  std::vector<char>&    DoMoY)
{
    std::ifstream is("swapvol.dat");
    SwapVolInput(Expiry, FwdMat, Matrix, DoMoY, is);
}
*/


static
void FixVolInput(MKTVOL_DATA*      mktvol,  /**< volatility data      */
                 char const*       Index,   /**< index to calibrate   */
                 T_CURVE const*    t_curve, /**< term structure data  */
                 std::istream&     is)      /**< input stream         */
{
    static char const*  routine = "FixVolInput ";

    mktvol->CalibFlag = TRUE;

    long IdxMat = atoi(Index) * 12;

    /* Read conventions from t_curve */
    mktvol->BaseDate = t_curve->Today;
    mktvol->Freq     = t_curve->SwapFreq;

    if (!strcmp(t_curve->SwapDCC, "360"))
    {
        mktvol->DCC = '0';
    }
    else if (!strcmp(t_curve->SwapDCC, "365"))
    {
        mktvol->DCC = '5';
    }
    else
    {
        mktvol->DCC = '3';
    }


    /* Read full swaption matrix */
    std::vector<IRDate>      Expiry;
    std::vector<IRDate>      FwdMat;
    matrix_t                VolMatrix;
    std::vector<char>       DoMoY;

    SwapVolInput(Expiry,
                 FwdMat,
                 VolMatrix,
                 DoMoY,
                 is);

    /* Process the final maturity index */
    size_t i;
    for (i = 0; i < Expiry.size(); i++)
    {
        long Mat = IdxMat - Expiry[i];

        /* We got out of the swaption matrix */
        if (Mat < FwdMat[0])
            break;

        size_t j = 0; 
        while (j < FwdMat.size()-1 && Mat >= FwdMat[j])
            j++;

        /* Use higher end of bracket: we don't interpolate to avoid stub */
        /* This includes the case FwdMat[i] > FwdMat[NbCol-1] so that    */
        /* we use a flat volatility after the last forward maturity.     */
        if (2 * Mat >= FwdMat[j-1] + FwdMat[j])
        {                                                                       
            Mat = FwdMat[j];                                          
            mktvol->Vol[i] = VolMatrix[i][j];
        }
        else                            
        {
            Mat = FwdMat[j-1];
            mktvol->Vol[i] = VolMatrix[i][j-1];
        }

        if( DoMoY[i] == 'M')
            mktvol->VolDate[i] = Nxtmth (mktvol->BaseDate, Expiry[i], 1L);
        else if (DoMoY[i] == 'Y')
            mktvol->VolDate[i] = Nxtmth (mktvol->BaseDate, 12 * Expiry[i], 1L);
        else if (DoMoY[i] == 'D')
            mktvol->VolDate[i] = Nxtday (mktvol->BaseDate, Expiry[i]);

        mktvol->SwapSt [i] = mktvol->VolDate[i];
        mktvol->SwapMat[i] = Nxtmth (mktvol->VolDate[i], Mat, 1L);
        mktvol->VolUsed[i] = TRUE;
    }
        
    if (i == 0)
        throw EslException(routine) << "nyFix calibration falls outside swaption matrix";

    mktvol->NbVol = i;
}






static
void CMSVolInput(MKTVOL_DATA*      mktvol,  /**< volatility data      */
                 char const*       Index,   /**< index to calibrate   */
                 T_CURVE const*    t_curve, /**< term structure data  */
                 std::istream&     is)      /**< input stream         */
{
    static char const*  routine = "CMSVolInput ";

    long IdxMat = atoi(Index) * 12;

    /* Read conventions from t_curve */
    mktvol->BaseDate = t_curve->Today;
    mktvol->Freq     = t_curve->SwapFreq;

    if (!strcmp(t_curve->SwapDCC, "360"))
    {
        mktvol->DCC = '0';
    }
    else if (!strcmp(t_curve->SwapDCC, "365"))
    {
        mktvol->DCC = '5';
    }
    else
    {
        mktvol->DCC = '3';
    }


    /* Read full swaption matrix */
    std::vector<IRDate>      Expiry;
    std::vector<IRDate>      FwdMat;
    matrix_t                VolMatrix;
    std::vector<char>       DoMoY;

    SwapVolInput(Expiry,
                 FwdMat,
                 VolMatrix,
                 DoMoY,
                 is);

    /* Find required Cms column */
    size_t j = 0; 
    while (j < FwdMat.size()-1 && IdxMat > FwdMat[j])
        j++;

    if (IdxMat != FwdMat[j])
        throw EslException(routine) << "Cms index is not in swaption matrix";

    mktvol->NbVol = Expiry.size();

    for (int i = 0; i < mktvol->NbVol; i++)
    {                                                                       
        /* 
         *  VolDate and SwapSt are identical: i.e. we don't calibrate 
         *  options on forward starting swaps (e.g. mid curve options) 
         */
        if( DoMoY[i] == 'M')
            mktvol->VolDate[i] = Nxtmth (mktvol->BaseDate, Expiry[i], 1L);
        else if (DoMoY[i] == 'Y')
            mktvol->VolDate[i] = Nxtmth (mktvol->BaseDate, 12 * Expiry[i], 1L);
        else if (DoMoY[i] == 'D')
            mktvol->VolDate[i] = Nxtday (mktvol->BaseDate, Expiry[i]);

        mktvol->SwapSt[i]  = mktvol->VolDate[i];
        mktvol->SwapMat[i] = Nxtmth (mktvol->VolDate[i], IdxMat, 1L);
        mktvol->Vol[i]     = VolMatrix[i][j];
        mktvol->VolUsed[i] = TRUE;
    }
}




static
void BaseVolInput(MKTVOL_DATA*      mktvol,  /**< volatility data      */
                  char const*       Index,   /**< index to calibrate   */
                  T_CURVE const*    t_curve, /**< term structure data  */
                  std::istream&     is)      /**< input stream         */
{
    static char const*  routine = "BaseVolInput ";

    mktvol->CalibFlag = TRUE;

    long IdxMat = atoi(Index);

    /* Read conventions from t_curve */
    mktvol->BaseDate = t_curve->Today;

    if (t_curve->MMB == 365)
    {
        mktvol->DCC = '5';
    }
    else
    {
        mktvol->DCC = '0';
    }

    /* Read vol curve */
    BaseVolInput(mktvol, is);

    if (12 / Conv_Freq (mktvol->Freq) != IdxMat)
        throw EslException(routine) << "Base vol curve frequency different from calibration";


    for (int i = 0; i < mktvol->NbVol; i++)
    {
        mktvol->SwapSt [i] = mktvol->VolDate[i];
        mktvol->SwapMat[i] = Nxtmth(mktvol->VolDate[i], IdxMat, 1L);
        mktvol->VolUsed[i] = TRUE;
    }
}



//-----------------------------------------------------------------------------
// Read vol data - stream version.
//-----------------------------------------------------------------------------
void EslMktVolInput(MKTVOL_DATA*      mktvol,  ///< volatility data
                    char const*       Index,   ///< index to calibrate
                    T_CURVE const*    t_curve, ///< term structure data
                    std::istream&     is)      ///< input stream
{
    static char const*  routine = "MktVolInput ";
    /* No calibration case */
    if (strstr (Index, "nil") != NULL)
    {
        BaseVolSkip(is);
        SwapVolSkip(is);

        mktvol->CalibFlag = FALSE;

        /* Initialize unused variables */
        mktvol->BaseDate = 0;
        mktvol->NbVol = 0;
        mktvol->Freq = 'z';
        mktvol->DCC = 'z';
        mktvol->SkipFlag = FALSE;
        return;
    }

    mktvol->CalibFlag = TRUE;

    /* Search for * character in calibration index name */
    if (strchr(Index, '*') == NULL)
        mktvol->SkipFlag = FALSE;  /* vol points skipping not allowed */
    else
        mktvol->SkipFlag = TRUE;   /* vol points skipping allowed      */


    if (strstr(Index, "yCms"))  // CMS indices
    {
        BaseVolSkip(is);
        CMSVolInput(mktvol, Index, t_curve, is);
    }
    else
    if (strstr (Index, "yFix")) // Final maturity indices
    {
        BaseVolSkip(is);
        FixVolInput(mktvol, Index, t_curve, is);
    }
    else
    if (strchr(Index, 'm'))     // Base vol indices
    {
        BaseVolInput(mktvol, Index, t_curve, is);
        SwapVolSkip(is);
    }
    else
        throw EslException(routine) << "Incorrect calibration index";

    if (MktVol_Check_W (mktvol) != SUCCESS)
        throw EslException(routine) << "Volatility matrix check failed";
}



//-----------------------------------------------------------------------------
// Read vol data - file version.
//-----------------------------------------------------------------------------
void EslMktVolInput(MKTVOL_DATA*      mktvol,  ///< volatility data
                    char const*       Index,   ///< index to calibrate
                    T_CURVE const*    t_curve) ///< term structure data
{
    static char const*  routine = "MktVolInput ";

    /* No calibration case */
    if (strstr (Index, "nil") != NULL)
    {
        mktvol->CalibFlag = FALSE;

        /* Initialize unused variables */
        mktvol->BaseDate = 0;
        mktvol->NbVol = 0;
        mktvol->Freq = 'z';
        mktvol->DCC = 'z';
        mktvol->SkipFlag = FALSE;
        return;
    }

    mktvol->CalibFlag = TRUE;

    /* Search for * character in calibration index name */
    if (strchr(Index, '*') == NULL)
        mktvol->SkipFlag = FALSE;  /* vol points skipping not allowed */
    else
        mktvol->SkipFlag = TRUE;   /* vol points skipping allowed      */


    if (strstr(Index, "yCms"))  // CMS indices
    {
        std::ifstream is("swapvol.dat");
        CMSVolInput(mktvol, Index, t_curve, is);
    }
    else
    if (strstr (Index, "yFix")) // Final maturity indices
    {
        std::ifstream is("swapvol.dat");
        FixVolInput(mktvol, Index, t_curve, is);
    }
    else
    if (strchr(Index, 'm'))     // Base vol indices
    {
        std::ifstream is("basevol.dat");
        BaseVolInput(mktvol, Index, t_curve, is);
    }
    else
        throw EslException(routine) << "Incorrect calibration index";

    if (MktVol_Check_W (mktvol) != SUCCESS)
        throw EslException(routine) << "Volatility matrix check failed";
}


