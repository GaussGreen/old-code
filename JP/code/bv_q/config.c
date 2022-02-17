/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: config.c
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/

#include "config.h"
#include "genaddin.h"
#include "matrixo.h"

#define  CATEGORY_ID  1

CATEGORY CATEGORIES[] = 
{
    { 
        "q3bv",
        "q3bv"
    } 
};

/*
***************************************************************************
** The ADDINS array describes each add-in function to be provided for this
** add-in.
**
** For each add-in we provide the following:
** 1. Name of the function (from the user's point of view).
** 2. Name of the wrapper function (a C-function called by the interface
**    layer).
** 3. Header file of the wrapper function.
** 4. Category number of the function.
** 5. Number of inputs.
** 6. Number of outputs.
** 7. Array of inputs + outputs.
** 8. Sundry flags.
** 9. Description of function.
**
** The array of inputs and outputs is an array of the PARAM data structure.
**
** Within the PARAM structure, there are the following values:
** 1. Name of the input or output field.
** 2. Data type of the field.
** 3. Array type of the field (default = "SCALAR").
** 4. Sundry flags.
** 5. Description of the field.
**
** Note that the description fields are currently being discarded by the
** code generation function.
**
** The ADDINS do not need to be in any particular order - the code
** generation process will sort them into alphabetical order. However it
** might help you when maintaining the file to keep them in order.
**
** The ADDIN and PARAM structures are defined in the config.h file in the
** iface/include directory.
***************************************************************************
*/

ADDIN ADDINS[] =
{
    {
        "PRICER",
        "Q3MQBivarPricerL",
        "mqbvl.h",
        CATEGORY_ID,
        26,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 27 */
        },
        0,
        NULL,
        ALIB_SC_0
    },
    {
        "BUS_PRICER",
        "Q3MQBivarBusPricerL",
        "mqbvl.h",
        CATEGORY_ID,
        27,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"HolFile",     "CHAR_BLOCK","ARRAY",        0, ""}, /* 26 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 27 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 28 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SM_PRICER",
        "Q3SmileBivarPricerL",
        "mqbvl.h",
        CATEGORY_ID,
        26,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 27 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SM_BUS_PRICER",
        "Q3SmileBivarBusPricerL",
        "mqbvl.h",
        CATEGORY_ID,
        27,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"HolFile",     "CHAR_BLOCK","ARRAY",        0, ""}, /* 26 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 27 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 28 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "MQ_MIDCURVE_PRICER",
        "Q3MQBivarMCMidCurvePricerL",
        "mqbvmcl.h",
        CATEGORY_ID,
        28,
        1,
        {   
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats",     "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats",     "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq",        "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC",         "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"MktFreq1",    "LONG",      "ARRAY:1",      0, ""}, /*  7 */
            {"MktDCC1",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  8 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  9 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 11 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 12 */
            {"MktDCC2",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 13 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 14 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 15 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 16 */
            {"FixedTyp",    "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 17 */
            {"stubTyp",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 18 */
            {"stubAtEnd",   "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 19 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 20 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:20", 0, ""}, /* 21 */
            {"VnfmPrms",    "DOUBLE",    "ARRAY",        0, ""}, /* 22 */
            {"Cross Corr",  "DOUBLE",    "ARRAY:1",      0, ""}, /* 23 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 25 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"Strike",      "DOUBLE",    "ARRAY:1",      0, ""}, /* 27 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 28 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 29 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SM_MIDCURVE_PRICER",
        "Q3SmileBivarMCMidCurvePricerL",
        "mqbvmcl.h",
        CATEGORY_ID,
        28,
        1,
        {
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats",     "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats",     "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq",        "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC",         "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"MktFreq1",    "LONG",      "ARRAY:1",      0, ""}, /*  7 */
            {"MktDCC1",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  8 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  9 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 11 */
            {"MktFreq2",    "LONG",      "ARRAY:1",      0, ""}, /* 12 */
            {"MktDCC2",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 13 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 14 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 15 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 16 */
            {"FixedTyp",    "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 17 */
            {"stubTyp",     "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 18 */
            {"stubAtEnd",   "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 19 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 20 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:20", 0, ""}, /* 21 */
            {"VnfmPrms",    "DOUBLE",    "ARRAY",        0, ""}, /* 22 */
            {"Corr",        "DOUBLE",    "ARRAY:1",      0, ""}, /* 23 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 25 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"Strike",      "DOUBLE",    "ARRAY:1",      0, ""}, /* 27 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 28 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 29 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "PRICER_SIMPLE",
        "Q3MQBivarPricerSimpleL",
        "mqbvl.h",
        CATEGORY_ID,
        13,
        1,
        {       
            {"Expiry",      "DOUBLE",    "ARRAY:1",      0, ""}, /*  1 */
            {"RateType1",   "LONG",      "ARRAY:1",      0, ""}, /*  2 */
            {"FwdRate1",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  3 */
            {"Smile1",      "DOUBLE",    "ARRAY",        0, ""}, /*  4 */
            {"SigATM1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  5 */
            {"RateType2",   "LONG",      "ARRAY:1",      0, ""}, /*  6 */
            {"FwdRate2",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  7 */
            {"Smile2",      "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"SigATM2",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"Corr",        "DOUBLE",    "ARRAY:1",      0, ""}, /* 10 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 11 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 12 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 14 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SM_PRICER_SIMPLE",
        "Q3SmileBivarPricerSimpleL",
        "mqbvl.h",
        CATEGORY_ID,
        13,
        1,
        {       
            {"Expiry",      "DOUBLE",    "ARRAY:1",      0, ""}, /*  1 */
            {"RateType1",   "LONG",      "ARRAY:1",      0, ""}, /*  2 */
            {"FwdRate1",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  3 */
            {"Smile1",      "DOUBLE",    "ARRAY",        0, ""}, /*  4 */
            {"SigATM1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  5 */
            {"RateType2",   "LONG",      "ARRAY:1",      0, ""}, /*  6 */
            {"FwdRate2",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  7 */
            {"Smile2",      "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"SigATM2",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"Corr",        "DOUBLE",    "ARRAY:1",      0, ""}, /* 10 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 11 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 12 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 14 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SIMPLE_CORR_CALIB",
        "Q3SmileBivarSimpleCorrCalibL",
        "mqbvl.h",
        CATEGORY_ID,
        13,
        1,
        {       
            {"Premium",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  1 */
            {"Expiry",      "DOUBLE",    "ARRAY:1",      0, ""}, /*  2 */
            {"FwdRate1",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  3 */
            {"Smile1",      "DOUBLE",    "ARRAY",        0, ""}, /*  4 */
            {"SigATM1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  5 */
            {"FwdRate2",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  6 */
            {"Smile2",      "DOUBLE",    "ARRAY",        0, ""}, /*  7 */
            {"SigATM2",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  8 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /*  9 */
            {"Strike",      "DOUBLE",    "ARRAY:1",      0, ""}, /* 10 */
            {"PayPrms",     "DOUBLE",    "ARRAY:2",      0, ""}, /* 11 */
            {"InitialGuess","DOUBLE",    "ARRAY:1",      0, ""}, /* 12 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 14 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

#ifdef DO_NOT_RELEASE
    {                                                                       
        "MC_PRICER",
        "Q3MQBivarMCPricerL",
        "mqbvmcl.h",
        CATEGORY_ID,
        26,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 27 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {                                                                       
        "MC_PRICER",
        "Q3SmileBivarMCPricerL",
        "mqbvmcl.h",
        CATEGORY_ID,
        26,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"",            "DOUBLE",    "ARRAY:3"            }  /* 27 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {                                                                       
        "FA_MARGINAL_DIST",
        "Q3MQFAMarginalDistL",
        "marginalDist.h",
        CATEGORY_ID,
        26,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 22 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 23 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 25 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 26 */
            {"outputs",     "MATRIX",    "SCALAR"            }  /* 27 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {                                                                       
        "MARGINAL_DIST",
        "Q3MQMarginalDistL",
        "marginalDist.h",
        CATEGORY_ID,
        13,
        1,
        {       
            {"Expiry",      "DOUBLE",    "ARRAY:1",      0, ""}, /*  1 */
            {"RateType1",   "LONG",      "ARRAY:1",      0, ""}, /*  2 */
            {"FwdRate1",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  3 */
            {"Smile1",      "DOUBLE",    "ARRAY",        0, ""}, /*  4 */
            {"SigATM1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  5 */
            {"RateType2",   "LONG",      "ARRAY:1",      0, ""}, /*  6 */
            {"FwdRate2",    "DOUBLE",    "ARRAY:1",      0, ""}, /*  7 */
            {"Smile2",      "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"SigATM2",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"Corr",        "DOUBLE",    "ARRAY:1",      0, ""}, /* 10 */
            {"OptTyp",      "LONG",      "ARRAY:1",      0, ""}, /* 11 */
            {"PayPrms",     "DOUBLE",    "ARRAY",        0, ""}, /* 12 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"outputs",     "MATRIX",    "SCALAR"            }
        },
        0,
        NULL,
        ALIB_SC_0
    },

#endif


    {
        "CORR",
        "Q3MQBivarCorrL",
        "mqbvl.h",
        CATEGORY_ID,
        24,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 22 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 23 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"",            "DOUBLE",    "ARRAY:1"            }  /* 25 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "SM_CORR",
        "Q3SmileBivarCorrL",
        "mqbvl.h",
        CATEGORY_ID,
        24,
        1,
        {       
            {"Today",       "DATE",      "ARRAY:1",      0, ""}, /*  1 */
            {"ValDate",     "DATE",      "ARRAY:1",      0, ""}, /*  2 */
            {"IdxDats1",    "DATE",      "ARRAY",        0, ""}, /*  3 */
            {"IdxRats1",    "DOUBLE",    "LINK_SIZE:3",  0, ""}, /*  4 */
            {"Freq1",       "LONG",      "ARRAY:1",      0, ""}, /*  5 */
            {"DCC1",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /*  6 */
            {"rateRstDts1", "DATE",      "ARRAY:3",      0, ""}, /*  7 */
            {"Sml1",        "DOUBLE",    "ARRAY",        0, ""}, /*  8 */
            {"AtmVol1",     "DOUBLE",    "ARRAY:1",      0, ""}, /*  9 */
            {"VnfmPrms1",   "DOUBLE",    "ARRAY",        0, ""}, /* 10 */
            {"IdxDats2",    "DATE",      "ARRAY",        0, ""}, /* 11 */
            {"IdxRats2",    "DOUBLE",    "LINK_SIZE:11", 0, ""}, /* 12 */
            {"Freq2",       "LONG",      "ARRAY:1",      0, ""}, /* 13 */
            {"DCC2",        "CHAR_BLOCK","ARRAY:1",      0, ""}, /* 14 */
            {"rateRstDts2", "DATE",      "ARRAY:3",      0, ""}, /* 15 */
            {"Sml2",        "DOUBLE",    "ARRAY",        0, ""}, /* 16 */
            {"AtmVol2",     "DOUBLE",    "ARRAY:1",      0, ""}, /* 17 */
            {"VnfmPrms2",   "DOUBLE",    "ARRAY",        0, ""}, /* 18 */
            {"DscDats",     "DATE",      "ARRAY",        0, ""}, /* 19 */
            {"DscRats",     "DOUBLE",    "LINK_SIZE:19", 0, ""}, /* 20 */
            {"Corr",        "DOUBLE",    "ARRAY:2",      0, ""}, /* 21 */
            {"PayDate",     "DATE",      "ARRAY:1",      0, ""}, /* 22 */
            {"SetlTyp",     "LONG",      "ARRAY:1",      0, ""}, /* 23 */
            {"Trace",       "LONG",      "ARRAY:1",      0, ""}, /* 24 */
            {"",            "DOUBLE",    "ARRAY:1"            }  /* 25 */
        },
        0,
        NULL,
        ALIB_SC_0
    },

    {
        "VERSION",             /* Addin function name     */
        "Q3BivarVersionL",     /* C Wrapper function name */
        "mqbvl.h",             /* Header file             */
        CATEGORY_ID,           /* Category number         */
        0,                     /* Number of Inputs        */
        1,                     /* number of outputs       */
        {               
            { "",       "CHAR_BLOCK",   "ARRAY:1"}
        },
        0,
        NULL
    },


    {
        "ERR_LOG",       /* Addin function name     */
        "GtoErrMsgSetL", /* C Wrapper function name */
        "cerrorl.h",     /* Header file             */
        CATEGORY_ID,     /* Category number         */
        1,               /* Number of Inputs        */
        0,               /* number of outputs       */
        {               
            { "Flag",    "LONG",    "ARRAY",   0,    "Input parameter 1."}
        },
        0,
        NULL
    },

    {
        "ERR_LOG_FILE_NAME",  /* Function name     */
        "GtoErrMsgFileNameL", /* C-function name   */
        "cerrorl.h",          /* Header file       */
        CATEGORY_ID,          /* Category number   */
        2,                    /* Number of inputs  */
        0,                    /* Number of outputs */
        {
            {"FileName", "CHAR_BLOCK", "SCALAR", 0,
             "Name of alternative error log file."},
            {"Append", "LONG", "SCALAR", 0, 
             "Should the file be appended to if it already exists."}
        },
        0,
        NULL,
        ALIB_SC_4
    },
        
    {
        "LOGGING",        /* Addin function name     */
        "GtoLoggingSetL", /* C Wrapper function name */
        "cerrsupl.h",     /* Header file             */
        CATEGORY_ID,      /* Category number         */
        1,                /* Number of Inputs        */
        0,                /* number of outputs       */
        {               
            { "Flag",   "LONG",         "ARRAY",   0,   "Input parameter 1."}
        },
        0,
        NULL
    }

};

size_t NUM_CATEGORIES = sizeof(CATEGORIES)/sizeof(CATEGORY);
size_t NUM_ADDINS = sizeof(ADDINS)/sizeof(ADDIN);
    

int main()
{
    int status;

    OLD_ADDIN *oldStyleAddinPtr=NULL;

#ifdef USE_OLD_STYLE_ADDIN
    OLD_ADDIN oldAddin;
    CATEGORY  oldAddinCategory;

    oldAddinCategory.name   = "JPM Analytics";
    oldAddinCategory.xlName = "JPM Analytics";

    oldAddin.name = "TEST";
    oldAddin.cat  = &oldAddinCategory;
#endif

    status = GtoGenerateAddinCode2 
        (".",              /* Generated code directory */
         "../src",         /* Source code directory */
         "Q3BV",           /* Function prefix */
         "Q3BV",           /* Library name - not used much */
         NULL,             /* Remote service name - not used here */
         TARGET,           /* Flags denoting what is to be generated */
         oldStyleAddinPtr, /* Definition of old-style add-in */
         NUM_ADDINS,       /* Number of add-in functions */
         ADDINS,           /* Definition of add-in functions */
         NUM_CATEGORIES,   /* Number of categories */
         CATEGORIES,       /* Definition of categories */
         "GtoMatrixRegister", /* Object registration function - not used here */
         FALSE,            /* Add hook functions - not used here */
         "");
    return status;
}

