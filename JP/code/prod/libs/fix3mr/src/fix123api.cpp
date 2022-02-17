// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/24/2003 Afshin Bayrooti
//
// $Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/fix123API.cpp,v 1.1 2004/04/19 15:52:05 markss Exp $
//

#pragma hdrstop

#ifdef _MSC_VER
#pragma warning(disable:4786)
#pragma warning(disable:4100)
#pragma warning(disable:4512)
#pragma warning(disable:4018)
#pragma warning(disable:4127)
#endif


#include "fix123API.h"
#include "fix123head.h"

#include <sstream>
#include <math.h>

#ifndef ERROR
#define  ERROR 1E-7  // Error tolerance used by fix3 - redefined here for the meantime
#endif

namespace IR {

// put here temporarily until a more logical place for these definitions is found
char Date::noLeapYear[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
char Date::leapYear[13] = {0,31,29,31,30,31,30,31,31,30,31,30,31};

// FIX3Error class static member declarations
string FIX3Error::_errorMessage;
bool FIX3Error::_errorDetected;


//
//   T _ C U R V E _ A P I
//

void ToPublicObject(const _T_CURVE* instancePtr, T_CURVE_API* publicInstancePtr)
{
    // Copy the data from the internal object to the public represenation
    publicInstancePtr->TodayDate = instancePtr->Today;
    publicInstancePtr->ValueDate = instancePtr->ValueDate;
    FrequencyEnumAsChar::toPublic(instancePtr->SwapFreq, publicInstancePtr->SwapFreq);
    YearBasisEnumAsString::toPublic(instancePtr->SwapDCC, publicInstancePtr->SwapDCC);
    publicInstancePtr->MMB = instancePtr->MMB;
    DynamicCArrayToVector(instancePtr->ZeroDate, instancePtr->NbZero, publicInstancePtr->ZeroDate);
    DynamicCArrayToVector(instancePtr->Zero, instancePtr->NbZero, publicInstancePtr->Zero);
}

void FromPublicObject(_T_CURVE* instancePtr, const T_CURVE_API* publicInstancePtr)
{
    // Validate the data

    if (publicInstancePtr->ZeroDate.size() != publicInstancePtr->Zero.size()) 
    {
        throw Error(std::string("Number of dates '") 
            << publicInstancePtr->ZeroDate.size()
            << "' does not equal number of rates '"
            << publicInstancePtr->Zero.size());
    }

    instancePtr->Today = publicInstancePtr->TodayDate;
    instancePtr->ValueDate = publicInstancePtr->ValueDate;
    FrequencyEnumAsChar::fromPublic(instancePtr->SwapFreq, publicInstancePtr->SwapFreq);
    YearBasisEnumAsString::fromPublic(instancePtr->SwapDCC, publicInstancePtr->SwapDCC);
    

    instancePtr->MMB = publicInstancePtr->MMB;
    VectorToCArray(publicInstancePtr->ZeroDate, instancePtr->ZeroDate);
    VectorToCArray(publicInstancePtr->Zero, instancePtr->Zero);

    // Fill in the rest of the structure
    // ?????????? Add code to Date class to calculate numdays between two dates
    instancePtr->SpotDays = 0;  // in the meantime set it to zero

    instancePtr->NbZero = publicInstancePtr->ZeroDate.size();
}


//
//   T R E E _ D A T A _  A P I
//

void ToPublicObject(const _TREE_DATA* instancePtr, TREE_DATA_API* publicInstancePtr)
{
    // Copy the data from the internal object to the public represenation

    publicInstancePtr->CvDiff = instancePtr->CvDiff;
    publicInstancePtr->CvDisc = instancePtr->CvDisc;
    publicInstancePtr->NbFactor = instancePtr->NbFactor;
    publicInstancePtr->NbSigmaMax = instancePtr->NbSigmaMax;
}

void FromPublicObject(_TREE_DATA* instancePtr, const TREE_DATA_API* publicInstancePtr)
{
    // Initialize the internal data representation
    Tree_Init(instancePtr);

    if ((publicInstancePtr->CvDiff < 0) || (publicInstancePtr->CvDiff > 2))
    {
        throw Error("Curve to Diffuse must equal 1, 2 or 3");
    }
    instancePtr->CvDiff = publicInstancePtr->CvDiff;

    if ((publicInstancePtr->CvDisc < 0) || (publicInstancePtr->CvDisc > 2))
    {
        throw Error("Curve to Discount must equal 1, 2 or 3");
    }
    instancePtr->CvDisc = publicInstancePtr->CvDisc;

    if ((publicInstancePtr->NbFactor != 1) &&
        (publicInstancePtr->NbFactor != 2) &&
        (publicInstancePtr->NbFactor != 3))
    {
        throw Error("Number of factors in tree must equal 1, 2 or 3");
    }
    instancePtr->NbFactor = publicInstancePtr->NbFactor;
    
    if (publicInstancePtr->NbSigmaMax < 3)
    {
        throw Error("Tree can not be cut at less than three standard deviations");
    }
    instancePtr->NbSigmaMax = publicInstancePtr->NbSigmaMax;
    
}


//
//   M K T V O L _ D A T A _ A P I
//

void ToPublicObject(const _MKTVOL_DATA* instancePtr, MKTVOL_DATA_API* publicInstancePtr)
{
    publicInstancePtr->BaseDate = instancePtr->BaseDate;
    publicInstancePtr->NbVol = instancePtr->NbVol;
    DynamicCArrayToVector(instancePtr->VolDate, instancePtr->NbVol, publicInstancePtr->VolDate);
    DynamicCArrayToVector(instancePtr->Vol, instancePtr->NbVol, publicInstancePtr->Vol);
    DynamicCArrayToVector(instancePtr->VolUsed, instancePtr->NbVol, publicInstancePtr->VolUsed);

    publicInstancePtr->Freq = instancePtr->Freq;
    publicInstancePtr->DCC = instancePtr->DCC;
    DynamicCArrayToVector(instancePtr->SwapSt, instancePtr->NbVol, publicInstancePtr->SwapSt);
    DynamicCArrayToVector(instancePtr->SwapMat, instancePtr->NbVol, publicInstancePtr->SwapMat);

    publicInstancePtr->CalibFlag = static_cast<bool>(instancePtr->CalibFlag);
    publicInstancePtr->SkipFlag = static_cast<bool>(instancePtr->SkipFlag);
}

void FromPublicObject(_MKTVOL_DATA* instancePtr, const MKTVOL_DATA_API* publicInstancePtr)
{   
    MktVol_Init(instancePtr);

    instancePtr->BaseDate = publicInstancePtr->BaseDate;
    int NbVol = publicInstancePtr->NbVol;
    instancePtr->NbVol = NbVol;
    CheckArraySize("Vols", NbVol, "Vol Dates", publicInstancePtr->VolDate.size());
    VectorToCArray(publicInstancePtr->VolDate, instancePtr->VolDate);
    CheckArraySize("Vols", NbVol, "Vols", publicInstancePtr->Vol.size());
    VectorToCArray(publicInstancePtr->Vol, instancePtr->Vol);

    CheckArraySize("VolUsed", NbVol, "Vols", publicInstancePtr->VolUsed.size());
    VectorToCArray(publicInstancePtr->VolUsed, instancePtr->VolUsed);

    instancePtr->Freq = publicInstancePtr->Freq;
    instancePtr->DCC = publicInstancePtr->DCC;
    CheckArraySize("Vols", NbVol, "Swap Start Dates", publicInstancePtr->SwapSt.size());
    VectorToCArray(publicInstancePtr->SwapSt, instancePtr->SwapSt);
    CheckArraySize("Vols", NbVol, "Swap Maturity Dates", publicInstancePtr->SwapMat.size());
    VectorToCArray(publicInstancePtr->SwapMat, instancePtr->SwapMat);

    instancePtr->CalibFlag = static_cast<int>(publicInstancePtr->CalibFlag);
    instancePtr->SkipFlag = static_cast<int>(publicInstancePtr->SkipFlag);

}


//
//   O P T _ O U T _ D A T A _ A P I
//

void ToPublicObject(const _OPT_OUT_DATA* instancePtr, OPT_OUT_DATA_API* publicInstancePtr)
{
    publicInstancePtr->Option = instancePtr->Option;
    CArrayToVector(instancePtr->Price, 10, publicInstancePtr->Price);
}

void FromPublicObject(_OPT_OUT_DATA* instancePtr, const OPT_OUT_DATA_API* publicInstancePtr)
{
    instancePtr->Option = publicInstancePtr->Option;
    VectorToCArray(publicInstancePtr->Price, instancePtr->Price);
}



void PopulateMktVolData(
    const MARKET_DATA &market_data,
    const string& calibrationIndex,
    MKTVOL_DATA &mktvol_data)
{
    MktVol_Init(&mktvol_data);

    // process market Volatilities - function throws errors as no heap memory allocation has
    // occurred yet
    CalculateMktVol(
		calibrationIndex, 
		market_data.t_curve[0], 
		market_data.baseVol_data, 
		market_data.swapVol_data, 
		mktvol_data);

}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


// Market Volatility Calculations

void CalculateMktVol(
    const string &calibrationIndex,
    const T_CURVE &t_curve,
    const BaseVolData &baseVol_data,
    const SwapVolData &swapVol_data,
    MKTVOL_DATA &mktvol_data)
{

    char    *StrIdx = NULL;     
    long    IdxMat;             // Maturity of the index (final or forward) 
    long    Mat;
    int     NbRows=0;           // Swaption volatility matrix 
    int     NbCol=0;
    int     i, j;
    int     status = FAILURE;   // Error status = FAILURE initially         
    char    Index[100];

    strcpy(Index, calibrationIndex.c_str());

    // No calibration case 
    if (strstr (Index, "nil") != NULL)
    {
        mktvol_data.CalibFlag = FALSE;

        // Initialize unused variables 
        mktvol_data.BaseDate = 0;
        mktvol_data.NbVol = 0;
        mktvol_data.Freq = 'z';
        mktvol_data.DCC = 'z';
        mktvol_data.SkipFlag = FALSE;

        return;
    }

    mktvol_data.CalibFlag = TRUE;


    // Search for * character in calibr ation index name 
    StrIdx = strchr(Index, '*');

    if (StrIdx == NULL)
    {
        mktvol_data.SkipFlag = FALSE;  // vol points skipping not allowed 
    }
    else
    {
        mktvol_data.SkipFlag = TRUE;   // vol points skipping allowed      
        *StrIdx = '\0';                 // terminated index name before '*' 
    }


    // Cms indices 
    StrIdx = strstr (Index, "yCms");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index) * 12;


        // Read conventions from t_curve 
        mktvol_data.BaseDate = t_curve.Today;
        mktvol_data.Freq     = t_curve.SwapFreq;

        if (!strcmp(t_curve.SwapDCC, "360"))
        {
            mktvol_data.DCC = '0';
        }
        else if (!strcmp(t_curve.SwapDCC, "365"))
        {
            mktvol_data.DCC = '5';
        }
        else
        {
            mktvol_data.DCC = '3';
        }

        // convert maturities from years to months
        vector<long> FwdSwapMaturity(swapVol_data.FwdSwapMaturity);
        NbRows = swapVol_data.SwaptionExpiry.size();
        NbCol = FwdSwapMaturity.size();
        // this should maybe be in a to/from public function
        for (i = 0; i < NbCol; i++)
            FwdSwapMaturity[i] *= 12;  // convert to months
        Matrix<double> VolMatrix(swapVol_data.VolMatrix);

        // Find required Cms column 
        j = 0; 
        while ((j < NbCol-1) && (IdxMat > FwdSwapMaturity[j]))
            j++;

        if (IdxMat != FwdSwapMaturity[j])
            throw Error("Cms index is not in swaption matrix!");
            
        mktvol_data.NbVol = NbRows;

        for (i = 0; i < mktvol_data.NbVol; i++)
        {                                                                       
            
            // VolDate and SwapSt are identical: i.e. we don't calibrate 
            // options on forward starting swaps (e.g. mid curve options) 
             
            mktvol_data.VolDate[i] = Nxtmth (mktvol_data.BaseDate, swapVol_data.SwaptionExpiry[i], 1L);
            mktvol_data.SwapSt[i]  = mktvol_data.VolDate[i];
            mktvol_data.SwapMat[i] = Nxtmth (mktvol_data.VolDate[i], IdxMat, 1L);
            mktvol_data.Vol[i]     = VolMatrix.at(i, j);
            mktvol_data.VolUsed[i] = TRUE;
        }

        MktVol_Check_W (&mktvol_data);

        return;
    }


    // Final maturity indices 
    StrIdx = strstr (Index, "yFix");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index) * 12;


        // Read conventions from t_curve 
        mktvol_data.BaseDate = t_curve.Today;
        mktvol_data.Freq     = t_curve.SwapFreq;

        if (!strcmp(t_curve.SwapDCC, "360"))
        {
            mktvol_data.DCC = '0';
        }
        else if (!strcmp(t_curve.SwapDCC, "365"))
        {
            mktvol_data.DCC = '5';
        }
        else
        {
            mktvol_data.DCC = '3';
        }


        // convert maturities from years to months
        vector<long> FwdSwapMaturity(swapVol_data.FwdSwapMaturity);
        NbRows = swapVol_data.SwaptionExpiry.size();
        NbCol = FwdSwapMaturity.size();
        for (i = 0; i < NbCol; i++)
            FwdSwapMaturity[i] *= 12;  // convert to months
        // convert swapvols from percent
        Matrix<double> VolMatrix(swapVol_data.VolMatrix);

        // Process the final maturity index 
        for (i = 0; i < NbRows; i++)
        {
            Mat = IdxMat - swapVol_data.SwaptionExpiry[i];

            // We got out of the swaption matrix 
            if (Mat < FwdSwapMaturity[0])
                break;

            j = 0; 
            while ((j < NbCol-1) && (Mat >= FwdSwapMaturity[j]))
                j++;

            // Use higher end of bracket: we don't interpolate to avoid stub 
            // This includes the case FwdMat[i] > FwdMat[NbCol-1] so that    
            // we use a flat volatility after the last forward maturity.     
            if (2 * Mat >= FwdSwapMaturity[j-1] + FwdSwapMaturity[j])
            {                                                                       
                Mat = FwdSwapMaturity[j];
                mktvol_data.Vol[i] = (VolMatrix).at(i,j);
            }
            else                            
            {
                Mat = FwdSwapMaturity[j-1];
                mktvol_data.Vol[i] = VolMatrix.at(i,j-1);
            }

            mktvol_data.VolDate[i] = Nxtmth (mktvol_data.BaseDate, swapVol_data.SwaptionExpiry[i], 1L);
            mktvol_data.SwapSt[i]  = mktvol_data.VolDate[i];
            mktvol_data.SwapMat[i] = Nxtmth (mktvol_data.VolDate[i], Mat, 1L);
            mktvol_data.VolUsed[i] = TRUE;

        }  // for i 
        
        if (i == 0)
            throw Error("nyFix calibration falls outside swaption matrix");

        mktvol_data.NbVol = i;

        MktVol_Check_W (&mktvol_data);
   
        return;
    }


    // Base vol indices 
    StrIdx = strchr(Index, 'm');

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index);


        // Read conventions from t_curve 
        mktvol_data.BaseDate = t_curve.Today;

        if (t_curve.MMB == 365)
        {
            mktvol_data.DCC = '5';
        }
        else
        {
            mktvol_data.DCC = '0';
        }

        //
        // Base Vol Processing
        //

        // check number of input vol dates equals number of input vols
        if (baseVol_data.VolDate.size() != baseVol_data.Vol.size())
            throw Error(string("Number of base vol dates (") << 
                baseVol_data.VolDate.size() << 
                ") must equal number of base vol values (" <<
                baseVol_data.Vol.size() << ")");
        mktvol_data.NbVol = baseVol_data.VolDate.size();
        mktvol_data.Freq = baseVol_data.Freq;

        int maxVols = CArraySize(mktvol_data.Vol);

        if (mktvol_data.NbVol > maxVols)
            throw Error(string("Number of input base vols (") <<
                mktvol_data.NbVol <<
                ") is greater than internal static array size (" <<
                maxVols << ")");

        VectorToCArray (baseVol_data.VolDate, mktvol_data.VolDate);
        VectorToCArray (baseVol_data.Vol, mktvol_data.Vol);

        // Eliminate dates falling before base date
        j = 0;
        while (mktvol_data.VolDate[j] <= mktvol_data.BaseDate)
            j++;

        mktvol_data.NbVol -= j;

        for (i = 0; i < mktvol_data.NbVol; i++)
        {
            mktvol_data.VolDate[i] = mktvol_data.VolDate[i+j];
            mktvol_data.Vol[i] = mktvol_data.Vol[i+j];
        }

        if (12 / Conv_Freq (mktvol_data.Freq) != IdxMat)
            throw Error("Base vol curve frequency different from calibration index");

        for (i = 0; i < mktvol_data.NbVol; i++)
        {
            mktvol_data.SwapSt[i]  = mktvol_data.VolDate[i];
            mktvol_data.SwapMat[i] = Nxtmth (mktvol_data.VolDate[i], IdxMat, 1L);
            mktvol_data.VolUsed[i] = TRUE;
        }

        MktVol_Check_W(&mktvol_data);
    }


    if (StrIdx == NULL)
        throw Error("Incorrect calibration index!");

    return;
}


void Param_Input (
    MKTVOL_DATA   *mktvol_data,   
    TREE_DATA     *tree_data,     
    const MODEL_DATA& model_data,
    const MARKET_DATA &market_data)
{

    int     i; 
    int     readerror;          // Reading error status
    double  norm;

    int NbFactor = tree_data->NbFactor;

    // reuse as much as possible of existing logic in function.  Could do with a rewrite eventually
    // when added into fix3 to support the new interface
    char OverWriteString[6][250];
    strcpy(OverWriteString[0], model_data.treePPY.c_str());
    strcpy(OverWriteString[1], model_data.modelType.c_str());
    strcpy(OverWriteString[2], model_data.factorWeights.c_str());
    strcpy(OverWriteString[3], model_data.meanReversion.c_str());
    strcpy(OverWriteString[4], model_data.correlations.c_str());
    strcpy(OverWriteString[5], model_data.backbone.c_str());

//    if (!market_data.useModelParameters)
    if (market_data.modelParameters_data.IsNull())
    {
        // If not using model parameters file, use model data entered with the instrument

        readerror = sscanf (OverWriteString[0],
                            "%d \n", 
                            &(tree_data->Ppy));

        if (readerror != 1)
        {        
            throw Error("Could not find Periods per year overwrite required! (Param_Input)");
        }

    
        if (strstr(OverWriteString[1], "N"))
        {
            mktvol_data->QLeft    = 0.;
            mktvol_data->QRight   = 0.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }
        else if (strstr(OverWriteString[1], "L"))
        {
            mktvol_data->QLeft    = 1.;
            mktvol_data->QRight   = 1.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }
        else
        {
            readerror = sscanf (OverWriteString[1],
                            "%lf %lf %lf %d\n", 
                            &(mktvol_data->QLeft),
                            &(mktvol_data->QRight),
                            &(mktvol_data->FwdShift),
                            &(mktvol_data->CetNbIter));

            if (readerror != 4)
            {      
                throw Error("Could not find Q overwrite required! (Param_Input)");
            }

            // At input level q=0 means log-normal whereas
            // internally q=0 means normal: we switch here
            mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
            mktvol_data->QRight = 1. - mktvol_data->QRight;
        }

        if (NbFactor == 1)
        {
            readerror = sscanf (OverWriteString[2], 
                                "%lf \n", 
                                &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                throw Error("Could not find factor weight overwrite required! (Param_Input)");
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                throw Error("Could not find mean reversion overwrite required! (Param_Input)");
            }

            // Fill in the unused parameters with N/A values

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
            readerror = sscanf (OverWriteString[2], 
                                "%lf \t%lf \n", 
                                &(mktvol_data->Alpha[0]),
                                &(mktvol_data->Alpha[1]));

            if (readerror != 2)
            {
                throw Error("Could not find file %s: factor weight overwrite required! (Param_Input)");
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]));

            if (readerror != 2)
            {        
                throw Error("Could not find file %s: mean reversion overwrite required! (Param_Input)");
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                throw Error("Could not find file %s: correlation overwrite required! (Param_Input)");
            }

            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;
        }
        else if (NbFactor == 3)
        {
            readerror = sscanf (OverWriteString[2], 
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Alpha[0]),
                                &(mktvol_data->Alpha[1]),
                                &(mktvol_data->Alpha[2]));

            if (readerror != 3)
            {        
                throw Error("Could not find file %s: factor weight overwrite required! (Param_Input)");
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]),
                                &(mktvol_data->Beta[2]));

            if (readerror != 3)
            {        
                throw Error("Could not find file %s: mean reversion overwrite required! (Param_Input)");
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Rho[0]),
                                &(mktvol_data->Rho[1]),
                                &(mktvol_data->Rho[2]));

            if (readerror != 3)
            {        
                throw Error("Could not find file %s: correlation overwrite required! (Param_Input)");
            }
        }  // if NbFactor == 3

        readerror = sscanf (OverWriteString[5],
                            "%lf \n", 
                            &(mktvol_data->Bbq));

        if (readerror != 1)
        {      
            throw Error("Could not find file %s: Bbq overwrite required! (Param_Input)");
        }

        // At input level q=0 means log-normal whereas
        // internally q=0 means normal: we switch here
        mktvol_data->Bbq  = 1. - mktvol_data->Bbq;
    
    }
    else  // use modelParameters Data 
    {        
        
         //  Parameter file exists: use the overwrite strings if they are not 
         // "nil", otherwise use the values in the parameter file.
         //
        if (NbFactor == 1)
        {
            mktvol_data->Beta[0] = market_data.modelParameters_data->oneFactorMR;

            // if overwrite exists, use it
            if (strstr (OverWriteString[3], "nil") == NULL)  // Need to use strstr() here, strcmp won't do 
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \n", 
                                    &(mktvol_data->Beta[0]));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for mean reversion! (Param_Input)");
                }
            }

            mktvol_data->Alpha[0] = market_data.modelParameters_data->oneFactorVol;

            // use overwrite if it exists
            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \n", 
                                    &(mktvol_data->Alpha[0]));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for factor weight! (Param_Input)");
                }
            }

            tree_data->Ppy = market_data.modelParameters_data->oneFactorPPY;

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for Ppy! (Param_Input)");
                }
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

            mktvol_data->Beta[0] = market_data.modelParameters_data->twoFactorMR[0];
            mktvol_data->Beta[1] = market_data.modelParameters_data->twoFactorMR[1];
            
            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]));

                if (readerror != 2)
                {        
                    throw Error("Could not read overwrite string for mean reversion! (Param_Input)");
                }
            }

            mktvol_data->Alpha[0] = market_data.modelParameters_data->twoFactorVol[0];
            mktvol_data->Alpha[1] = market_data.modelParameters_data->twoFactorVol[1];

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]));

                if (readerror != 2)
                {        
                    throw Error("Could not read overwrite string for factor weight! (Param_Input)");
                }
            }

            mktvol_data->Rho[0] = market_data.modelParameters_data->twoFactorCorrelation;

            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \n", 
                                    &(mktvol_data->Rho[0]));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for correlation! (Param_Input)");
                }
            }

            tree_data->Ppy = market_data.modelParameters_data->twoFactorPPY;

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for Ppy! (Param_Input)");
                }
            }

            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;

        }
        else if (NbFactor == 3)
        {

            mktvol_data->Beta[0] = market_data.modelParameters_data->threeFactorMR[0];
            mktvol_data->Beta[1] = market_data.modelParameters_data->threeFactorMR[1];
            mktvol_data->Beta[2] = market_data.modelParameters_data->threeFactorMR[2];

            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]),
                                    &(mktvol_data->Beta[2]));

                if (readerror != 3)
                {        
                    throw Error("Could not read overwrite string for mean reversion! (Param_Input)");
                }
            }

            mktvol_data->Alpha[0] = market_data.modelParameters_data->threeFactorVol[0];
            mktvol_data->Alpha[1] = market_data.modelParameters_data->threeFactorVol[1];
            mktvol_data->Alpha[2] = market_data.modelParameters_data->threeFactorVol[2];
            
            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]),
                                    &(mktvol_data->Alpha[2]));

                if (readerror != 3)
                {        
                    throw Error("Could not read overwrite string for factor weight! (Param_Input)");
                }
            }

            mktvol_data->Rho[0] = market_data.modelParameters_data->threeFactorCorrelation12;
            mktvol_data->Rho[1] = market_data.modelParameters_data->threeFactorCorrelation13;
            mktvol_data->Rho[2] = market_data.modelParameters_data->threeFactorCorrelation23;

            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Rho[0]),
                                    &(mktvol_data->Rho[1]),
                                    &(mktvol_data->Rho[2]));

                if (readerror != 3)
                {        
                    throw Error("Could not read overwrite string for correlation! (Param_Input)");
                }
            }

            tree_data->Ppy = market_data.modelParameters_data->threeFactorPPY;

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    throw Error("Could not read overwrite string for Ppy! (Param_Input)");
                }
            }

        }  // if NbFactor ==

        mktvol_data->QLeft = market_data.modelParameters_data->qLeft;
        mktvol_data->QRight = market_data.modelParameters_data->qRight;
        mktvol_data->FwdShift = market_data.modelParameters_data->FwdShift;
        mktvol_data->CetNbIter = market_data.modelParameters_data->CetNbIter;

        mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
        mktvol_data->QRight = 1. - mktvol_data->QRight;

        if (strstr (OverWriteString[1], "nil") == NULL)
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

                readerror = sscanf (OverWriteString[1],
                                    "%lf %lf %lf %d\n", 
                                    &(mktvol_data->QLeft),
                                    &(mktvol_data->QRight),
                                    &(mktvol_data->FwdShift),
                                    &(mktvol_data->CetNbIter));

                if (readerror != 4)
                {      
                    throw Error("Could not read overwrite string for Qweight! (Param_Input)");
                }       
                
                mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                mktvol_data->QRight = 1. - mktvol_data->QRight;
            }
        }  // if then else
        
        if (strstr (OverWriteString[5], "nil") != NULL)
        {
            mktvol_data->Bbq = market_data.modelParameters_data->backBone;

// ?????????????? if backone is optional parameter, must check if it is supplied by
// the user as it is required if the instrument model data sets this to null.            
    
            mktvol_data->Bbq = 1. - mktvol_data->Bbq;

        }
        else
        {
            readerror = sscanf (OverWriteString[5],
                                "%lf \n", 
                                &(mktvol_data->Bbq));

            if (readerror != 1)
            {      
                throw Error("Could not read overwrite string for Backbone! (Param_Input)");
            }       
            
            mktvol_data->Bbq  = 1. - mktvol_data->Bbq;

        }  // if then else

    }  // if then else
    
    // Set total vol constants
    norm = 0.;
    for (i = 0; i < NbFactor; i++) norm += mktvol_data->Alpha[i] * mktvol_data->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        throw Error("Total alpha is too small !");
    }       
    if (IS_EQUAL(mktvol_data->Bbq,1))
    {
        mktvol_data->VolNorm = 0.;
        mktvol_data->VolLogn = norm;
    }
    else
    if (IS_EQUAL(mktvol_data->Bbq,0))
    {
        mktvol_data->VolNorm = norm;
        mktvol_data->VolLogn = 0.;
    }
    else
    {
        throw Error("Bbq parameter must be either 0 or 1 !");
    }       

    // Check validity of input
    if (Param_Check(NbFactor, mktvol_data, tree_data) == FAILURE)              
    {
        throw Error(FIX3Error::GetError());
    }
        

} 

// formulate error message and throw error
void CheckArraySize( string name, int size, string vectorName, int vectorSize )
{
    if (vectorSize != size)
    {
        std::ostringstream oss;
        oss << "Number of " << vectorName << " entered (" << vectorSize <<
            ") must equal number of " << name << " (" << size << ")";

        throw Error(oss.str());
    }
}




} // IR
