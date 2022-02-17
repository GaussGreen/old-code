#include <iostream>
#include <fstream>
#include <math.h>
#include "esl_log.h"
#include "fix123head.h"

#include <sys/types.h>
#include <sys/stat.h>


void Fix3_ParamInput (   
             MKTVOL_DATA*    mktvol_data,                /* Volatility data               */
             FIX3_TREE_DATA* tree_data,                  /* Tree data structure           */
             int             NbFactor,                   /* Number of factors             */
             char            OverWriteString[6][MAXBUFF],/* Overwrite strings             */
             std::istream&   is)                         /* input stream                  */
{

    static char const* routine = "Fix3_ParamInput ";
    int     i; 

    if (NbFactor < 1 || NbFactor > 3)
        throw EslException(routine) << "Nb of factors must be 1, 2 or 3";

    char    buf[1024];

    EslHeaderInput("model", is);
    is.getline(buf, sizeof(buf));

    if (NbFactor == 1)
    {
        // one factor mean reversion
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[0]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        // overwrite it
        if (strstr(OverWriteString[3], "nil") == NULL)
        {
            if (sscanf(OverWriteString[3], "%lf", &mktvol_data->Beta[0]) != 1) 
                throw EslException(routine) << "Could not read mean reversion override";
        }


        // one factor weight
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[0]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        // overwrite it
        if (strstr(OverWriteString[2], "nil") == NULL)
        {
            if (sscanf(OverWriteString[2], "%lf", &mktvol_data->Alpha[0]) != 1) 
                throw EslException(routine) << "Could not read factor weight override";
        }


        // one factor ppy
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%d", &tree_data->Ppy) != 1) 
            throw EslException(routine) << "Could not read ppy";

        // overwrite it
        if (strstr(OverWriteString[0], "nil") == NULL)
        {
            if (sscanf(OverWriteString[0], "%d", &tree_data->Ppy) != 1) 
                throw EslException(routine) << "Could not read ppy override";
        }

        /* Skip two and three factor parameters in the file */
        for (i = 0; i < 12; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is)
                throw EslException(routine) << "Could not find two factor parameters";
        }

        for (i = 0; i < 20; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is) 
                throw EslException(routine) << "Could not find three factor parameters";
        }

        mktvol_data->Alpha[1] = -999.;
        mktvol_data->Alpha[2] = -999.;
        mktvol_data->Beta[1]  = -999.;
        mktvol_data->Beta[2]  = -999.;
        mktvol_data->Rho[0]   = -999.;
        mktvol_data->Rho[1]   = -999.;
        mktvol_data->Rho[2]   = -999.;

    }
    else if (NbFactor == 2)
    {
        /* Skip one factor parameters in the file */
        for (i = 0; i < 5; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is) 
                throw EslException(routine) << "Could not find one factor parameters";
        }

        // two factor mean reversion
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[0]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[1]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        // overwrite it
        if (strstr(OverWriteString[3], "nil") == NULL)
        {
            if (sscanf(OverWriteString[3], "%lf %lf", 
                        &mktvol_data->Beta[0],
                        &mktvol_data->Beta[1]) != 2) 
                throw EslException(routine) << "Could not read mean reversion override";
        }


        // two factor weight
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[0]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[1]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        // overwrite it
        if (strstr(OverWriteString[2], "nil") == NULL)
        {
            if (sscanf(OverWriteString[2], "%lf %lf", 
                        &mktvol_data->Alpha[0],
                        &mktvol_data->Alpha[1]) != 2) 
                throw EslException(routine) << "Could not read factor weight override";
        }


        // two factor correlation
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Rho[0]) != 1) 
            throw EslException(routine) << "Could not read factor correlation";

        // overwrite it
        if (strstr(OverWriteString[4], "nil") == NULL)
        {
            if (sscanf(OverWriteString[4], "%lf", 
                        &mktvol_data->Rho[0]) != 1) 
                throw EslException(routine) << "Could not read factor correlation override";
        }


        // two factor ppy
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%d", &tree_data->Ppy) != 1) 
            throw EslException(routine) << "Could not read ppy";

        // overwrite it
        if (strstr(OverWriteString[0], "nil") == NULL)
        {
            if (sscanf(OverWriteString[0], "%d", &tree_data->Ppy) != 1) 
                throw EslException(routine) << "Could not read ppy override";
        }


        for (i = 0; i < 20; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is) 
                throw EslException(routine) << "Could not find three factor parameters";
        }

        mktvol_data->Alpha[2] = -999.;
        mktvol_data->Beta[2]  = -999.;
        mktvol_data->Rho[1]   = -999.;
        mktvol_data->Rho[2]   = -999.;

    }
    else if (NbFactor == 3)
    {
        /* Skip one and two factor parameters in the file */
        for (i = 0; i < 5; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is) 
                throw EslException(routine) << "Could not find one factor parameters";
        }

        for (i = 0; i < 12; i++)
        {
            is.getline(buf, sizeof(buf));
            if (!is) 
                throw EslException(routine) << "Could not find two factor parameters";
        }

        // three factor mean reversion
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[0]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[1]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Beta[2]) != 1) 
            throw EslException(routine) << "Could not read mean reversion";

        // overwrite it
        if (strstr(OverWriteString[3], "nil") == NULL)
        {
            if (sscanf(OverWriteString[3], "%lf %lf %lf", 
                        &mktvol_data->Beta[0],
                        &mktvol_data->Beta[1],
                        &mktvol_data->Beta[2]) != 3) 
                throw EslException(routine) << "Could not read mean reversion override";
        }


        // three factor weight
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[0]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[1]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Alpha[2]) != 1) 
            throw EslException(routine) << "Could not read factor weight";

        // overwrite it
        if (strstr(OverWriteString[2], "nil") == NULL)
        {
            if (sscanf(OverWriteString[2], "%lf %lf %lf", 
                        &mktvol_data->Alpha[0],
                        &mktvol_data->Alpha[1],
                        &mktvol_data->Alpha[2]) != 3) 
                throw EslException(routine) << "Could not read factor weight override";
        }


        // three factor correlation
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Rho[0]) != 1) 
            throw EslException(routine) << "Could not read factor correlation";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Rho[1]) != 1) 
            throw EslException(routine) << "Could not read factor correlation";

        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%lf", &mktvol_data->Rho[2]) != 1) 
            throw EslException(routine) << "Could not read factor correlation";

        // overwrite it
        if (strstr(OverWriteString[4], "nil") == NULL)
        {
            if (sscanf(OverWriteString[4], "%lf %lf %lf", 
                        &mktvol_data->Rho[0],
                        &mktvol_data->Rho[1],
                        &mktvol_data->Rho[2]) != 3) 
                throw EslException(routine) << "Could not read factor correlation override";
        }


        // three factor ppy
        is.getline(buf, sizeof(buf));
        is.getline(buf, sizeof(buf));
        if (!is || sscanf(buf, "%d", &tree_data->Ppy) != 1) 
            throw EslException(routine) << "Could not read ppy";

        // overwrite it
        if (strstr(OverWriteString[0], "nil") == NULL)
        {
            if (sscanf(OverWriteString[0], "%d", &tree_data->Ppy) != 1) 
                throw EslException(routine) << "Could not read ppy override";
        }
    }

    // QLeft
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%lf", &mktvol_data->QLeft) != 1) 
        throw EslException(routine) << "Could not read QLeft";

    // QRight
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%lf", &mktvol_data->QRight) != 1) 
        throw EslException(routine) << "Could not read QRight";

    // forward shift
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%lf", &mktvol_data->FwdShift) != 1) 
        throw EslException(routine) << "Could not read forward shift";

    // number of calibration iterations
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%d", &mktvol_data->CetNbIter) != 1) 
        throw EslException(routine) << "Could not read number of calibration iterations";

    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
    mktvol_data->QRight = 1. - mktvol_data->QRight;

    // overwrite it
    if (strstr(OverWriteString[1], "nil") == NULL)
    {
        if (strstr(OverWriteString[1], "N") != NULL)
        {
            mktvol_data->QLeft    = 0.;
            mktvol_data->QRight   = 0.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }       
        else if (strstr(OverWriteString[1], "L") != NULL)
        {
            mktvol_data->QLeft    = 1.;
            mktvol_data->QRight   = 1.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }
        else
        {
            if (sscanf(OverWriteString[1], "%lf %lf %lf %d", 
                        &mktvol_data->QLeft,
                        &mktvol_data->QRight,
                        &mktvol_data->FwdShift,
                        &mktvol_data->CetNbIter) != 4)
                throw EslException(routine) << "Could not read Q overrides";
            
            mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
            mktvol_data->QRight = 1. - mktvol_data->QRight;
        }
    }
    
    is.getline(buf, sizeof(buf));
    is.getline(buf, sizeof(buf));
    if (!is || sscanf(buf, "%lf", &mktvol_data->Bbq) != 1) 
        throw EslException(routine) << "Could not read backbone";

    // overwrite it
    if (strstr(OverWriteString[5], "nil") == NULL)
    {
        if (sscanf(OverWriteString[5], "%lf", &mktvol_data->Bbq) != 1)
            throw EslException(routine) << "Could not read backbone override";
    }
    mktvol_data->Bbq  = 1. - mktvol_data->Bbq;

    EslFooterInput("model", is);

    if (Fix3_Param_Check(NbFactor, mktvol_data, tree_data) != SUCCESS)              
        throw EslException(routine) << "Model parameters check failed";
}


void Fix3_ParamInput (   
             MKTVOL_DATA*    mktvol_data,                /* Volatility data               */
             FIX3_TREE_DATA* tree_data,                  /* Tree data structure           */
             int             NbFactor,                   /* Number of factors             */
             char            OverWriteString[6][MAXBUFF])/* Overwrite strings             */
{
    std::ifstream is("modelParameters.dat");
    Fix3_ParamInput(mktvol_data, tree_data, NbFactor, OverWriteString, is);
}


static void copyStream(std::ostream& os, std::string const& file)
{
    char buf[1024];

    std::ifstream is(file.c_str());

    while(is) {
        is.getline(buf, sizeof(buf));
        if (buf[0]) // skipe empty lines
            os << buf << '\n';
    }
}

// concatenate input files into one stream
void Fix3_ConcatenateFix3Input(char const* output, char const* deal)
{
    static char const*  routine = "Fix3_ConcatenateFix3Input ";
    std::ofstream os(output);

    if (!os)
        throw EslException(routine) << "Could not open output file";

    // copy deal file
    os << "<deal>\n";
    copyStream(os, deal);
    os << "</deal>\n";

    // copy zero curve
    os << "<zero>\n";
    copyStream(os, "zero.dat");
    os << "</zero>\n";

    // copy discount curve
    os << "<disczero>\n";
    copyStream(os, "disczero.dat");
    os << "</disczero>\n";

    // copy risky curve
    os << "<riskzero>\n";
    struct stat sstr;
    if (stat("riskzero.dat", &sstr) == 0)
        copyStream(os, "riskzero.dat");
    else
        copyStream(os, "disczero.dat");
    os << "</riskzero>\n";

    // copy basevol
    os << "<basevol>\n";
    copyStream(os, "basevol.dat");
    os << "</basevol>\n";

    // copy swapvol
    os << "<swapvol>\n";
    copyStream(os, "swapvol.dat");
    os << "</swapvol>\n";

    // copy model parameters
    os << "<model>\n";
    copyStream(os, "modelParameters.dat");
    os << "</model>\n";
}


