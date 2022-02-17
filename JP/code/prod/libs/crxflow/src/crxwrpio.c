/*********************************************************************************
 * CRXWRPIO.C
 * read IR and credit info from files
 *
 ********************************************************************************/

#include "crxwrpio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <alib/zr2simp.h>     /* GtoZerosToSimplePoint  */
#include <alib/zr2coup.h>     /* GtoZerosToCouponsPoint */
#include <alib/gtomat.h>      /* GtoArray2DNew          */

#include "crcrv.h"
#include "crxmacros.h"

static char irZCFile0[]  = "ir_curve0";
static char irZCFile1[]  = "ir_curve1";
static char irZCFile2[]  = "ir_curve2";
static char irInfoFile[] = "ir_info";
static char irMParFile[] = "ir_mpar";
static char irVolFile[]  = "ir_voldiag";
static char irVolFile0[]  = "ir_voldiag0";

static char crZCFile[]   = "cr_zcurve";
static char crVolFile[]  = "cr_voldiag";
static char crMParFile[] = "cr_mpar";
static char crInfoFile[] = "cr_info";

static char bsZCFile[]   = "bs_zcurve";
static char bsVolFile[]  = "bs_voldiag0";
static char bsMParFile[] = "bs_mpar";
static char bsInfoFile[] = "bs_info";

static char numFile[]    = "numerics";

/*********************************************************************************
 * Conver Date and Curve to Alib
 ********************************************************************************/
int CrxCRInputFillAlibData(
    CR_INPUT *crInput);                   /* (I/O) IR_INPUT                     */

/*********************************************************************************
 * Conver Date and Curve to Alib
 ********************************************************************************/
int CrxIRInputFillAlibData(
    IR_INPUT *irInput);                   /* (I/O) IR_INPUT                     */

/*!
 *   Conver Date and Curve to Alib
 *   -        crInput:                       (I/O) IR_INPUT                     
 */
int CrxCRInputFillAlibData(
    CR_INPUT *crInput)                    /* (I/O) IR_INPUT                     */
{
    static char routine[] = "CRInputFillAlibData";
    int         status    = FAILURE;
    int         i;
    
    /* convert dates */
    if ((crInput->AlibBaseDate = CrxDrDate2TDate(crInput->BaseDate)) < 0)
    {
	DR_Error("%s: CR base date = %8d", routine, crInput->BaseDate);
	goto RETURN;
    }
    
    
    for(i = 0;i < crInput->NbVol;i++)
    {
        if ((crInput->AlibSwapSt[i]  = CrxDrDate2TDate(crInput->SwapSt[i]))  < 0 ||
            (crInput->AlibSwapMat[i] = CrxDrDate2TDate(crInput->SwapMat[i])) < 0 ||
            (crInput->AlibVolDate[i] = CrxDrDate2TDate(crInput->VolDate[i])) < 0)
    	{
            DR_Error("%s Swap [%d] start date = %8d, mat date = %8d, Vol date = %8d",
                     routine, i, 
                     crInput->SwapSt[i],
                     crInput->SwapMat[i],
                     crInput->VolDate[i]);
            goto RETURN;
        }
    }

    /* convert curves. */
    if(CrxDrCurve2TCurve(&crInput->AlibZeroCurve,
                         &crInput->ZeroCurve) == FAILURE)
    {
        DR_Error("Conversion of irInput ZCurve %d failed.\n",i);
        goto RETURN;
    }
    
    status = SUCCESS;
 RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/*!
 * Conver Date and Curve to Alib
 * -          irInput                        (I/O) IR_INPUT                     
 */
int CrxIRInputFillAlibData(
    IR_INPUT *irInput)                    /* (I/O) IR_INPUT                     */
{
    static char routine[] = "IRInputFillAlibData";
    int         status    = FAILURE;
    int         i;
    
    /* convert dates */
    if ((irInput->AlibBaseDate = CrxDrDate2TDate(irInput->BaseDate)) < 0)
    {
	DR_Error("%s: IR base date = %8d", routine, irInput->BaseDate);
	goto RETURN;
    }

    if(irInput->CorrSwapSt)
    {
        if ((irInput->AlibCorrSwapSt = CrxDrDate2TDate(irInput->CorrSwapSt)) < 0)
        {
	    DR_Error("%s: IR swap start date = %8d", routine, irInput->CorrSwapSt);
	    goto RETURN;
        }
    }
    else 
        irInput->AlibCorrSwapSt = 0;

    if (irInput->CorrSwapMat)
    {
        if ((irInput->AlibCorrSwapMat = CrxDrDate2TDate(irInput->CorrSwapMat)) < 0)
        {
	    DR_Error("%s: IR swap mat date = %8d", routine, irInput->CorrSwapMat);
	    goto RETURN;
        }
    }
    else
         irInput->AlibCorrSwapMat = 0;

    for(i = 0;i < irInput->NbVol;i++){
        if ((irInput->AlibSwapSt[i]  = CrxDrDate2TDate(irInput->SwapSt[i]))  < 0 ||
            (irInput->AlibSwapMat[i] = CrxDrDate2TDate(irInput->SwapMat[i])) < 0 ||
            (irInput->AlibVolDate[i] = CrxDrDate2TDate(irInput->VolDate[i])) < 0)
    	{
            DR_Error("%s Swap [%d] start date = %8d, mat date = %8d, Vol date = %8d",
                     routine, i, 
                     irInput->SwapSt[i],
                     irInput->SwapMat[i],
                     irInput->VolDate[i]);
            goto RETURN;
        }
    }

    /* convert curves, for curve 0, 1, and 2. */
    for(i = 0;i < 3;i++){
        if(CrxDrCurve2TCurve(&irInput->AlibZeroCurve[i],
                             &irInput->ZeroCurve[i]) == FAILURE)
        {
            DR_Error("Conversion of irInput ZCurve %d failed.\n",i);
            goto RETURN;
        }
    }

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
    return status;
}

/*!
 *   Conver Date and Curve to Alib
 *   -        spInput:                       (I/O) SP_INPUT                     
 */
int CrxSPInputFillAlibData(
    SP_INPUT *spInput)                    /* (I/O) IR_INPUT                     */
{
    static char routine[] = "SPInputFillAlibData";
    int         status    = FAILURE;
    int         i;
    
    /* convert dates */
    if ((spInput->AlibBaseDate = CrxDrDate2TDate(spInput->BaseDate)) < 0)
    {
	DR_Error("%s: SP base date = %8d", routine, spInput->BaseDate);
	goto RETURN;
    }
    
    
    for(i = 0;i < spInput->NbVol;i++)
    {
        if ((spInput->AlibSwapSt[i]  = CrxDrDate2TDate(spInput->SwapSt[i]))  < 0 ||
            (spInput->AlibSwapMat[i] = CrxDrDate2TDate(spInput->SwapMat[i])) < 0 ||
            (spInput->AlibVolDate[i] = CrxDrDate2TDate(spInput->VolDate[i])) < 0)
    	{
            DR_Error("%s Swap [%d] start date = %8d, mat date = %8d, Vol date = %8d",
                     routine, i, 
                     spInput->SwapSt[i],
                     spInput->SwapMat[i],
                     spInput->VolDate[i]);
            goto RETURN;
        }
    }

    /* convert curves. */
    if(CrxDrCurve2TCurve(&spInput->AlibZeroCurve,
                         &spInput->BasisZCurve) == FAILURE)
    {
        DR_Error("Conversion of irInput ZCurve %d failed.\n",i);
        goto RETURN;
    }
    
    status = SUCCESS;
 RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/* free MAW Input */
void CrxMAWFreeInput(
    CRX_INPUT     *crxInput)              /* (I)  CRX Input                     */
{
    int         i, j;

    if(crxInput != NULL)
    {
        /* free IRInput */
        for(i = 0;i < crxInput->nbIRInput;i++)
        {
            if(&(crxInput->irInput[i]) != NULL)
            {
                for(j = 0;j < 3;j++)
                {
                    GtoFreeTCurve((crxInput->irInput[i]).AlibZeroCurve[j]);
                }
/*
                free(&(crxInput->irInput[i]));
*/
            }
        }

        if(crxInput->irInput != NULL) free(crxInput->irInput);
        
        /* free CRInput */
        for(i = 0;i < crxInput->nbCRInput;i++)
        {
            if(&(crxInput->crInput[i]) != NULL)
            {
                
                GtoFreeTCurve((crxInput->crInput[i]).AlibZeroCurve);
/*
                free(&(crxInput->crInput[i]));
*/
            }
        }

        if(crxInput->crInput != NULL) free(crxInput->crInput);


        /* free correlation matrix */
        GtoArray2DFree((void **)crxInput->corrInput);
        
    }
}

/**
 * Read input from wrapper files
 * - Input
 *   -# pathDir:                                  path to input files
 * - Output
 *   -# Today:                                    today
 *   -# crxInput                                  CRX Input
 */
int CrxMAWReadInput(
    long          *Today,                 /* (O)  today                         */
    CRX_INPUT     *crxInput,              /* (O)  CRX Input                     */
    char          *pathDir)               /* (I)  path to input files           */
    
{
    static      char  routine[] = "CrxMAWReadInput";
    int         status    = FAILURE;

    char        fnameZC0[256];            /* file name for zcurve0           */
    char        fnameZC1[256];            /* file name for zcurve1           */
    char        fnameZC2[256];            /* file name for zcurve2           */
    char        fnameInfo[256];           /* file name for zcurve2           */
    char        fnameVolDiag[256];        /* file name for zcurve2           */
    char        fnameModlPar[256];        /* file name for zcurve2           */
    
    int         i, j;

    char        ***CorrOWS_MAW = NULL;   /* set to 'nil'     */
    int         NbAssets_MAW=0;          /* Total # of MAW assets */

    /* init crx input ----------------------------------------------------------*/
    crxInput->nbIRInput = 0;
    crxInput->irInput   = NULL;

    crxInput->nbCRInput = 0;
    crxInput->crInput   = NULL;

    crxInput->corrInput = NULL;

    
    /* read summary.dat and allocate memory to input structures */
    if (summaryInputEnvData("summary.dat", 
                            Today,
                            &crxInput->nbIRInput,
                            &(crxInput->irInput),
                            NULL,          
                            NULL,          
                            NULL,
                            NULL,
                            NULL,          
                            NULL,          
                            &crxInput->nbCRInput,          
                            &(crxInput->crInput)) == FAILURE) goto RETURN;

    
    /* read IR input -----------------------------------------------------------*/
    for(i = 0;i < crxInput->nbIRInput;i++){
        if(pathDir != NULL){
            sprintf(fnameZC0,    "%s/%s_%d.dat",pathDir,irZCFile0,i);
            sprintf(fnameZC1,    "%s/%s_%d.dat",pathDir,irZCFile1,i);
            sprintf(fnameZC2,    "%s/%s_%d.dat",pathDir,irZCFile2,i);
        
            sprintf(fnameInfo,   "%s/%s_%d.dat",pathDir,irInfoFile,i);
            sprintf(fnameVolDiag,"%s/%s_%d.dat",pathDir,irVolFile,i);
            sprintf(fnameModlPar,"%s/%s_%d.dat",pathDir,irMParFile,i);
        } else {
            sprintf(fnameZC0,    "%s_%d.dat",irZCFile0,i);
            sprintf(fnameZC1,    "%s_%d.dat",irZCFile1,i);
            sprintf(fnameZC2,    "%s_%d.dat",irZCFile2,i);
        
            sprintf(fnameInfo,   "%s_%d.dat",irInfoFile,i);
            sprintf(fnameVolDiag,"%s_%d.dat",irVolFile,i);
            sprintf(fnameModlPar,"%s_%d.dat",irMParFile,i);
        }
      
        if(irInputEnvData(&(crxInput->irInput[i]),
                          fnameZC0,
                          fnameZC1,
                          fnameZC2,
                          fnameInfo,
                          fnameVolDiag,
                          fnameModlPar,
                          "Y",
                          "nil","nil","nil",
                          "nil","nil","nil") == FAILURE)
        {
            DR_Error("Fail to process env data for IR.");
            goto RETURN;
        }

        /* fill out Alib TDate and TCurve */
        if(CrxIRInputFillAlibData(&(crxInput->irInput[i]))){
            DR_Error("Fail to convert data for IR to Alib format.");
            goto RETURN;
        }

    }
    
    
    
    /* read CR input -----------------------------------------------------------*/
    for(i = 0;i < crxInput->nbCRInput;i++){
        if(pathDir != NULL){
            sprintf(fnameZC0,    "%s/%s_%d.dat",pathDir,crZCFile,i);
            sprintf(fnameInfo,   "%s/%s_%d.dat",pathDir,crInfoFile,i);
            sprintf(fnameVolDiag,"%s/%s_%d.dat",pathDir,crVolFile,i);
            sprintf(fnameModlPar,"%s/%s_%d.dat",pathDir,crMParFile,i);
        } else {
            sprintf(fnameZC0,    "%s_%d.dat",crZCFile,i);
            sprintf(fnameInfo,   "%s_%d.dat",crInfoFile,i);
            sprintf(fnameVolDiag,"%s_%d.dat",crVolFile,i);
            sprintf(fnameModlPar,"%s_%d.dat",crMParFile,i);
        }
    

        if(crInputEnvData(&(crxInput->crInput[i]),
                          fnameZC0,
                          fnameInfo,
                          fnameVolDiag,
                          fnameModlPar,
                          "Y",
                          "1",
                          "nil", //alpha
                          "nil", //beta
                          "nil", //rho
                          "nil", //backbone
                          "nil") == FAILURE)
        {
            DR_Error("Fail to process env data for CR.");
            goto RETURN;
        }
        
        /* fill out Alib TDate and TCurve */
        if(CrxCRInputFillAlibData(&(crxInput->crInput[i]))){
            DR_Error("Fail to convert data for IR to Alib format.");
            goto RETURN;
        }

    }


    /* read IR-CR correlation --------------------------------------------------*/
    NbAssets_MAW = crxInput->nbIRInput + crxInput->nbCRInput;
    if (Alloc_CorrOWS(NbAssets_MAW, &(CorrOWS_MAW)) == FAILURE) goto RETURN;

    for (i=0; i<NbAssets_MAW; i++)
    {
        for (j=0; j<NbAssets_MAW; j++)
        {
            strcpy (CorrOWS_MAW[i][j], "nil");
        }
    }

    /* Allocate memory and initialize corr_MAW matrix */
    crxInput->corrInput = (double **)GtoArray2DNew((int)NbAssets_MAW, 
                                                   (int)NbAssets_MAW,
                                                   sizeof(double));
    if (crxInput->corrInput == NULL) goto RETURN;

    if (corrInputEnvData(crxInput->corrInput,
                         NbAssets_MAW,
                         "correlation.dat",
                         CorrOWS_MAW) == FAILURE) goto RETURN;


    
    status = SUCCESS;

 RETURN:
    
    Free_CorrOWS(NbAssets_MAW, CorrOWS_MAW);

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/* free basis MAW Input */
void BSMAWFreeInput(
    BS_INPUT     *bsInput)              /* (I)  CRX Input                     */
{
    int         i, j;

    if(bsInput != NULL)
    {
        /* free IRInput */
        for(i = 0;i < bsInput->nbIRInput;i++)
        {
            if(&(bsInput->irInput[i]) != NULL)
            {
                for(j = 0;j < 3;j++)
                {
                    GtoFreeTCurve((bsInput->irInput[i]).AlibZeroCurve[j]);
                }
            }
        }

        if(bsInput->irInput != NULL) free(bsInput->irInput);
        
        /* free CRInput */
        for(i = 0;i < bsInput->nbSPInput;i++)
        {
            if(&(bsInput->spInput[i]) != NULL)
            {
                
                GtoFreeTCurve((bsInput->spInput[i]).AlibZeroCurve);
            }
        }

        if(bsInput->spInput != NULL) free(bsInput->spInput);


        /* free correlation matrix */
        GtoArray2DFree((void **)bsInput->corrInput);
        
    }
}

/**
 * Read input from wrapper files
 * - Input
 *   -# pathDir:                                  path to input files
 * - Output
 *   -# Today:                                    today
 *   -# bsInput                                   basis Input
 */
int BSMAWReadInput(
    long          *Today,                 /* (O)  today                       */
    BS_INPUT      *bsInput,               /* (O)  basis Input                 */
    char          *pathDir)               /* (I)  path to input files         */
    
{
    static      char  routine[] = "BSMAWReadInput";
    int         status    = FAILURE;

    char        fnameZC0[256];            /* file name for zcurve0           */
    char        fnameZC1[256];            /* file name for zcurve1           */
    char        fnameZC2[256];            /* file name for zcurve2           */
    char        fnameInfo[256];           /* file name for info              */
    char        fnameVolDiag[256];        /* file name for voldiag           */
    char        fnameModlPar[256];        /* file name for model par         */
    char        fnameNum[256];            /* file name for numerics          */
  
    int         i, j;

    char        ***CorrOWS_MAW = NULL;   /* set to 'nil'     */
    int         NbAssets_MAW=0;          /* Total # of MAW assets */

    FILE        *stream = NULL;

    /* init crx input ----------------------------------------------------------*/
    bsInput->nbIRInput = 0;
    bsInput->irInput   = NULL;

    bsInput->nbSPInput = 0;
    bsInput->spInput   = NULL;

    bsInput->corrInput = NULL;

    
    /* read summary.dat and allocate memory to input structures */
    if (summaryInputEnvData("summary.dat", 
                            Today,
                            &bsInput->nbIRInput,
                            &(bsInput->irInput),
                            &bsInput->nbSPInput,
                            &(bsInput->spInput),
                            NULL,          
                            NULL,          
                            NULL,
                            NULL,
                            NULL,          
                            NULL) == FAILURE) goto RETURN;

    
    /* read IR input ---------------------------------------------------------*/
    for(i = 0;i < bsInput->nbIRInput;i++){
        if(pathDir != NULL){
            sprintf(fnameZC0,    "%s/%s_%d.dat",pathDir,irZCFile0,i);
            sprintf(fnameZC1,    "%s/%s_%d.dat",pathDir,irZCFile1,i);
            sprintf(fnameZC2,    "%s/%s_%d.dat",pathDir,irZCFile2,i);
        
            sprintf(fnameInfo,   "%s/%s_%d.dat",pathDir,irInfoFile,i);
            sprintf(fnameVolDiag,"%s/%s_%d.dat",pathDir,irVolFile0,i);
            sprintf(fnameModlPar,"%s/%s_%d.dat",pathDir,irMParFile,i);
        } else {
            sprintf(fnameZC0,    "%s_%d.dat",irZCFile0,i);
            sprintf(fnameZC1,    "%s_%d.dat",irZCFile1,i);
            sprintf(fnameZC2,    "%s_%d.dat",irZCFile2,i);
        
            sprintf(fnameInfo,   "%s_%d.dat",irInfoFile,i);
            sprintf(fnameVolDiag,"%s_%d.dat",irVolFile0,i);
            sprintf(fnameModlPar,"%s_%d.dat",irMParFile,i);
        }
      
        if(irInputEnvData(&(bsInput->irInput[i]),
                          fnameZC0,
                          fnameZC1,
                          fnameZC2,
                          fnameInfo,
                          fnameVolDiag,
                          fnameModlPar,
                          "Y",
                          "nil","nil","nil",
                          "nil","nil","nil") == FAILURE)
        {
            DR_Error("Fail to process env data for IR.");
            goto RETURN;
        }

        /* fill out Alib TDate and TCurve */
        if(CrxIRInputFillAlibData(&(bsInput->irInput[i]))){
            DR_Error("Fail to convert data for IR to Alib format.");
            goto RETURN;
        }

    }
    
    
    
    /* read SP input ---------------------------------------------------------*/
    for(i = 0;i < bsInput->nbSPInput;i++){
        if(pathDir != NULL){
            sprintf(fnameZC0,    "%s/%s_%d.dat",pathDir,bsZCFile,i);
            sprintf(fnameInfo,   "%s/%s_%d.dat",pathDir,bsInfoFile,i);
            sprintf(fnameVolDiag,"%s/%s_%d.dat",pathDir,bsVolFile,i);
            sprintf(fnameModlPar,"%s/%s_%d.dat",pathDir,bsMParFile,i);
        } else {
            sprintf(fnameZC0,    "%s_%d.dat",bsZCFile,i);
            sprintf(fnameInfo,   "%s_%d.dat",bsInfoFile,i);
            sprintf(fnameVolDiag,"%s_%d.dat",bsVolFile,i);
            sprintf(fnameModlPar,"%s_%d.dat",bsMParFile,i);
        }
    

        if(spInputEnvData(&(bsInput->spInput[i]),
                          fnameZC0,
                          fnameInfo,
                          fnameVolDiag,
                          fnameModlPar,
                          "Y",
                          "1",
                          "nil", 
                          "nil", 
                          "nil", 
                          "nil", 
                          "nil") == FAILURE)
        {
            DR_Error("Fail to process env data for CR.");
            goto RETURN;
        }
        
        /* fill out Alib TDate and TCurve */
        if(CrxSPInputFillAlibData(&(bsInput->spInput[i]))){
            DR_Error("Fail to convert data for IR to Alib format.");
            goto RETURN;
        }

    }


    /* read IR-SP correlation ------------------------------------------------*/
    NbAssets_MAW = bsInput->nbIRInput + bsInput->nbSPInput;
    if (bsInput->nbSPInput > 0)
    {
        if (Alloc_CorrOWS(NbAssets_MAW, &(CorrOWS_MAW)) == FAILURE) 
            goto RETURN;
        
        for (i=0; i<NbAssets_MAW; i++)
        {
            for (j=0; j<NbAssets_MAW; j++)
            {
                strcpy (CorrOWS_MAW[i][j], "nil");
            }
        }

        /* Allocate memory and initialize corr_MAW matrix */
        bsInput->corrInput = (double **)GtoArray2DNew(
                                (int)NbAssets_MAW, 
                                (int)NbAssets_MAW,
                                sizeof(double));
    
        if (bsInput->corrInput == NULL) goto RETURN;

        if (corrInputEnvData(bsInput->corrInput,
                         NbAssets_MAW,
                         "correlation.dat",
                         CorrOWS_MAW) == FAILURE) goto RETURN;
    }

    /* read PPY, NbCut, NbCET ------------------------------------------------*/
    if(pathDir != NULL)
        sprintf(fnameNum,    "%s/%s.dat",pathDir,numFile);
    else
        sprintf(fnameNum,    "%s.dat",numFile);

    stream = fopen(fnameNum, "r");
    if (stream == NULL)
    {
        bsInput->PPY      = 4;
        bsInput->nbCutoff = 5;
        bsInput->nbCetIter= 0;
    }
    else
    {
        int     readerror;          /* Reading error status */
        char    ErrorMsg[MAXBUFF];

        /* Read PPY */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "PPY", 
                                "BSMAWReadInput", 
                                fnameNum) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%d \n", &(bsInput->PPY));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read PPY in %s! (BSMAWReadInput)", 
                     fnameNum);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* Read Cutoff stdev*/
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Cutoff stdev", 
                                "BSMAWReadInput", 
                                fnameNum) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%lf \n", &(bsInput->nbCutoff));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read cutoff stdev in %s! (BSMAWReadInput)", 
                     fnameNum);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* Read nb CET iter */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "NB CET iter", 
                                "BSMAWReadInput", 
                                fnameNum) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%d \n", &(bsInput->nbCetIter));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read nb CET iter in %s! (BSMAWReadInput)", 
                     fnameNum);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }
    

    status = SUCCESS;

 RETURN:
    if (stream != NULL)
        fclose(stream);
    
    Free_CorrOWS(NbAssets_MAW, CorrOWS_MAW);

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/*********************************************************************************
 *       This function does:
 *       1) Mode 0: Find and skip mode
 *           Skip everything until it finds line start with '#' (ignore spaces)
 *           Flag an error message if it can't find
 *       2) Mode 1: Skip mode
 *           Skip all the empty lines and spaces
 *           Read the first non-empty character and check it is '#' 
 *           Flag an error message if it isn't.
 *
 *       The routine returns with the file pointer at the start of the next line
 *
 ********************************************************************************/
int     FindAndSkipComLine (int      Mode,
                            FILE    *stream,
                            char     *SectionLabel,
                            char     *Routine,
                            char     *FileName)
{
    char    ErrorMsg[MAXBUFF];
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i;
    int     length;

    if (stream == NULL) goto RETURN;

    if ((Mode != 0) && (Mode != 1))
    {
        DR_Error("FindAndSkipComLine: Mode has to be 0 or 1!");
        goto RETURN;
    }

    while (fgets (string, MAXBUFF, stream) != NULL)
    {
        length = strlen(string);

        for (i = 0; i < length; i++)
        {
            if (string[i] == '#')
            {
                status = SUCCESS;
                goto RETURN;
            }
            else if ((string[i] != ' ' ) && 
                     (string[i] != '\n') && 
                     (string[i] != '\t')) 
            {
                if (Mode == 0) break;
                else goto RETURN;
            }
        }
    }  

    RETURN:


    if (status == FAILURE)
    {
        sprintf (ErrorMsg, "%s: %s section line expected and not found in file"
                    " %s!", Routine, SectionLabel, FileName);
        DR_Error(ErrorMsg);
    }
        
    return (status);

}  /* FindAndSkipComLine */


/*********************************************************************************
 *      Functioning the same as FindAndSkipComLine() with added feature:
 *      3) Mode 2: same as Mode 0 but suppressing error message
 *      4) check optionalSeq number '0','1', ... '9',
 *      if '1' is passed in, only line begin with #1 will return SUCCESS
 ********************************************************************************/
int     FindAndSkipComLineOptional (int      Mode,
                                    FILE    *stream,
                                    char     *SectionLabel,
                                    char     *Routine,
                                    char     *FileName,
                                    char     optionalSeq)
{
    char    ErrorMsg[MAXBUFF];
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i;
    int     length;

    if (stream == NULL) goto RETURN;

    if ((Mode != 0) && (Mode != 1) && (Mode != 2))
    {
        DR_Error("FindAndSkipComLineOptional: Mode has to be 0, 1, 2!");
        goto RETURN;
    }

    while (fgets (string, MAXBUFF, stream) != NULL)
    {
        length = strlen(string);

        for (i = 0; i < length; i++)
        {
            if (string[i] == '#')
            {
                if (i+1 < length && string[i+1] == optionalSeq)
                {
                    status = SUCCESS;
                }
                goto RETURN;
            }
            else if ((string[i] != ' ' ) && 
                     (string[i] != '\n') && 
                     (string[i] != '\t')) 
            {
                if (Mode == 0 || Mode == 2) break;
                else goto RETURN;
            }
        }
    }  

    RETURN:


    if (status == FAILURE && Mode != 2)
    {
        sprintf (ErrorMsg, "%s: %s section line expected and not found in file"
                    " %s!", Routine, SectionLabel, FileName);
        DR_Error(ErrorMsg);
    }
        
    return (status);

}  /* FindAndSkipComLineOptional */

/*********************************************************************************
 *       Read term structure input for DR Wrapper.
 *
 ********************************************************************************/
int  Term_Input_W (T_CURVE   *t_curve,    /* (O) Structure of zero curve data   */
                   char      *FileName)   /* (I) File name including extension  */
{
    int     i;
    int     readerror;          /* Reading error status */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;
    char    CurveDCC[MAXBUFF];
   


    /* Basic Check */
    if ((t_curve  == NULL) ||
        (FileName == NULL)) goto RETURN;

    /* Open the yield curve data file (see termodel.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, stream, "start date", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", &(t_curve->ValueDate));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read value date in %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Today's date and spot days are not available */
    t_curve->Today = t_curve->ValueDate;
    t_curve->SpotDays  = 0;
    
    if (FindAndSkipComLine (1, stream, "dcc for zero rates", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%s \n", CurveDCC);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read dcc for zero rates in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (strstr(CurveDCC, "30/360") != NULL)
    {
        t_curve->CurveDCC = '3'; 
    }
    else if (strstr(CurveDCC, "ACT/365F") != NULL)
    {
        t_curve->CurveDCC = '5';
    }
    else if (strstr(CurveDCC, "ACT/360") != NULL)
    {
        t_curve->CurveDCC = '0';
    }
    else
    {
        sprintf (ErrorMsg, "dcc in file %s has to be '30/360' or"
                " 'ACT/365F' or 'ACT/360'! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

        
    if (FindAndSkipComLine (1, stream, "compounding frequency", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(t_curve->CurveFreq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read compounding frequency in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (1, stream, "number of zeros", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(t_curve->NbZero));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of zeros in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if ((t_curve->NbZero > MAX_ZRATE) ||
        (t_curve->NbZero < 0))
    {        
        sprintf (ErrorMsg, "Nb of rates in file %s must be <= %d and >= 0! (Term_Input_W)", FileName, MAX_ZRATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if (FindAndSkipComLine (1, stream, "zero dates and rates", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < t_curve->NbZero; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \n",
                            &(t_curve->ZeroDate[i]),
                            &(t_curve->Zero[i]));
                        
        t_curve->Zero[i] /= 100.;
                
        if (readerror != 2)
        {        
            sprintf (ErrorMsg, "Could not read zero date and rate #%d in file %s! (Term_Input_W)", i + 1, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    t_curve->InterpType = SRM3_LINEAR_INTERP;

    status = SUCCESS;
        
    RETURN:

    if (status == FAILURE)
    {
        DR_Error("Term_Input_W: failed!");
    }

    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

}  /* Term_Input_W */

/*********************************************************************************
 *       Check validity of DR Wrapper term structure inputs.
 *
 ********************************************************************************/
int     Term_Check_W (T_CURVE  *t_curve)  /* (I) Structure of zero curve data   */
{

    int     i;
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];

    if (t_curve == NULL) goto RETURN;

    if (Dateok(t_curve->ValueDate))
    {
        sprintf (ErrorMsg, "Incorrect format for value date! (Term_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((t_curve->CurveDCC != '0') && 
        (t_curve->CurveDCC != '3') &&
        (t_curve->CurveDCC != '5'))
    {
        sprintf (ErrorMsg, "DCC should be '0' or '3' or '5'! (Term_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }

    if (t_curve->CurveDCC != '5')
    {
        sprintf (ErrorMsg, "only 'ACT/365F' is supported for DCC! (Term_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }


    if (   (t_curve->CurveFreq != 'A')
        && (t_curve->CurveFreq != 'S')
    	&& (t_curve->CurveFreq != 'Q')
    	&& (t_curve->CurveFreq != 'M'))
    {
        sprintf (ErrorMsg, "Incorrect compounding frequency (A,S,Q or M)! (Term_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }

    if (t_curve->CurveFreq != 'A')
    {
        sprintf (ErrorMsg, "Only curves of annual frequency is supported! (Term_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }

    if ((t_curve->NbZero > MAX_ZRATE) ||
        (t_curve->NbZero < 0))
    {        
        sprintf (ErrorMsg, "Nb of rates must be <= %d and >= 0! (Term_Input_W)", MAX_ZRATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    for (i = 0; i < t_curve->NbZero; i++)
    {
        if (Dateok(t_curve->ZeroDate[i]))
        {
            sprintf (ErrorMsg, "Incorrect format for zero date! (Term_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    for (i = 1; i < t_curve->NbZero; i++)
        if (t_curve->ZeroDate[i] <= t_curve->ZeroDate[i-1])
        {
            sprintf (ErrorMsg, "Zero dates must be entered in ascending order! (Term_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

    if (t_curve->NbZero != 0)
    {
        if (t_curve->ZeroDate[0] <= t_curve->ValueDate)
        {
            sprintf (ErrorMsg, "Some zero dates are <= Value date! (Term_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if(t_curve->InterpType!=SRM3_LINEAR_INTERP && t_curve->InterpType!=SRM3_FLATFWD_INTERP)
    {
            sprintf (ErrorMsg, "Interpolation Type is not supported! (Term_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;        
    }


    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Term_Check_W: failed!");
    }


    return (status);

}  /* Term_Check_W */

/*********************************************************************************
 *       Read volatility input
 *
 ********************************************************************************/
int     VolDiag_Input_W (
            long    *BaseDate,            /* (O) Volatility data                */
            int     *VolUnit,
            int     *NbVol,
            long    *VolDate,
            long    *SwapSt,
            long    *SwapMat,
            double  *Vol,
            char    *VolType,
            char    *FileName)            /* (I) File name including extension  */
{
    int     i;
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    if ((BaseDate == NULL) ||
        (VolUnit  == NULL) ||
        (NbVol    == NULL) ||
        (VolDate  == NULL) ||
        (SwapSt   == NULL) ||
        (SwapMat  == NULL) ||
        (Vol      == NULL) ||
        (VolType  == NULL) ||
        (FileName == NULL)) goto RETURN;

    /* Open and read the volatility data file */
	stream = fopen (FileName, "r");


    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (VolDiag_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
	
	
    if (FindAndSkipComLine (1, stream, "start date", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }


    readerror = fscanf (stream, "%ld \n", BaseDate);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read start date in file %s! (VolDiag_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (FindAndSkipComLine (1, stream, "volatility unit", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }


    readerror = fscanf (stream, "%d \n", VolUnit);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read vol unit in file %s! (VolDiag_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, stream, "Nb of vol points", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }


	readerror = fscanf (stream, "%d \n", NbVol);

    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of vols in file %s! (VolDiag_Input_W)", FileName);
		DR_Error(ErrorMsg);
		goto RETURN;
    }


    if ((*NbVol > MAX_VOL) ||
        (*NbVol < 0))
    {        
        sprintf (ErrorMsg, "Nb of vols in file %s must be >= 0 and <= %d! (VolDiag_Input_W)", FileName, MAX_VOL);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if (FindAndSkipComLine (1, stream, "vol dates and rates", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

	
	for (i = 0; i < *NbVol; i++)
    {
	  readerror = fscanf (stream, "%ld %ld %ld %lf %c\n", 
						  &(VolDate[i]),
						  &(SwapSt[i]),
						  &(SwapMat[i]),
						  &(Vol[i]),
						  &(VolType[i]));
		
        if (readerror != 5)
        {        
            sprintf (ErrorMsg, "Could not read vol dates/rate/type #%d in file %s! (VolDiag_Input_W)", i+1, FileName);
			DR_Error(ErrorMsg);
			goto RETURN;
        }

        Vol[i] /= 100.;
    }

    status = SUCCESS;


 RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status == FAILURE)
    {
        DR_Error("VolDiag_Input_W: failed!");
    }
    

    return (status);

}  /* VolDiag_Input_W */

/*********************************************************************************
 *       Read volatility input
 *
 ********************************************************************************/
int     VolDiagCR_Input_W (
    long    *BaseDate,                    /* (O) Volatility data                */
    int     *VolUnit,
    int     *NbVol,
    long    *VolDate,
    long    *SwapSt,
    long    *SwapMat,
    double  *Vol,
    char    *FileName)                    /* (I) File name including extension  */
{
    int     i;
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    if ((BaseDate == NULL) ||
        (VolUnit  == NULL) ||
        (NbVol    == NULL) ||
        (VolDate  == NULL) ||
        (SwapSt   == NULL) ||
        (SwapMat  == NULL) ||
        (Vol      == NULL) ||
        (FileName == NULL)) goto RETURN;

    /* Open and read the volatility data file */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (VolDiagCR_Input_W)", FileName);
        DR_Error(ErrorMsg);

        sprintf (ErrorMsg, "Using one-point volatility default");
        DR_Error(ErrorMsg);

		NbVol = 0;
		status = SUCCESS;

        goto RETURN;
    }
        
    if (FindAndSkipComLine (1, stream, "start date", "VolDiagCR_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", BaseDate);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read start date in file %s! (VolDiagCR_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (FindAndSkipComLine (1, stream, "volatility unit", "VolDiagCR_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", VolUnit);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read vol unit in file %s! (VolDiagCR_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, stream, "Nb of vol points", "VolDiagCR_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbVol);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of vols in file %s! (VolDiagCR_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if ((*NbVol > MAX_VOL) ||
        (*NbVol < 0))
    {        
        sprintf (ErrorMsg, "Nb of vols in file %s must be >= 0 and <= %d! (VolDiagCR_Input_W)", FileName, MAX_VOL);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, stream, "vol dates and rates", "VolDiagCR_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < *NbVol; i++)
    {
        readerror = fscanf (stream, "%ld %ld %ld %lf \n", 
                            &(VolDate[i]),
                            &(SwapSt[i]),
                            &(SwapMat[i]),
                            &(Vol[i])
                            );

        if (readerror != 4)
        {        
            sprintf (ErrorMsg, "Could not read vol dates/rate/type #%d in file %s! (VolDiagCR_Input_W)", i+1, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        Vol[i] /= 100.;
    }

    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status == FAILURE)
    {
        DR_Error("VolDiagCR_Input_W: failed!");
    }
    

    return (status);

}  /* VolDiagCR_Input_W */




/*********************************************************************************
 *       Check validity of volatility inputs.
 *
 ********************************************************************************/
int     VolDiag_Check_W (long    BaseDate, /* Volatility data                   */
                         int     VolUnit,
                         int     NbVol,
                         long    *VolDate,
                         long    *SwapSt,
                         long    *SwapMat,
                         double  *Vol,
                         char    *VolType)    
{
    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];

    /* Basic Check */
    if ((VolDate == NULL) ||
        (SwapSt  == NULL) ||
        (SwapMat == NULL) ||
        (Vol     == NULL) ||
        (VolType == NULL)) goto RETURN;


    if (Dateok(BaseDate))
    {
        sprintf (ErrorMsg, "Incorrect format for start date! (VolDiag_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((VolUnit != 0) &&
        (VolUnit != 1))
    {
        sprintf (ErrorMsg, "Volatility unit must be 0 or 1! (VolDiag_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (VolUnit != 0)
    {
        sprintf (ErrorMsg, "Only 0 is supported for volatility unit! (VolDiag_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((NbVol > MAX_VOL) ||
        (NbVol < 0))
    {        
        sprintf (ErrorMsg, "Nb of vols must be >= 0 and <= %d! (VolDiag_Check_W)", MAX_VOL);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    for (i=0; i<NbVol; i++)
    {
        if ((Dateok(VolDate[i]))||
            (Dateok(SwapSt [i]))||
            (Dateok(SwapMat[i]))) 
        {
            sprintf (ErrorMsg, "incorrect date format #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (SwapMat[i] <= SwapSt[i])
        {
            sprintf (ErrorMsg, "rate end <= rate start #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (VolDate[i] >  SwapSt[i])
        {
            sprintf (ErrorMsg, "expiration date > rate start #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (Vol[i] < 0.0)
        {
            sprintf (ErrorMsg, "vol %d < 0! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if ((VolType[i] != 'C') &&
            (VolType[i] != 'S'))
        {
            sprintf (ErrorMsg, "Volatility type #%d has to be 'C' or 'S' ! "
                               "(VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }


    for (i = 1; i < NbVol; i++)
    {
        if (VolDate[i] <= VolDate[i-1])
        {
            sprintf (ErrorMsg, "Volatilty dates must be entered in ascending order! (VolDiag_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (NbVol != 0)
    {
        if (BaseDate >= VolDate[0])
        {
            sprintf (ErrorMsg, "Base date must be before vol date! (VolDiag_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    status = SUCCESS;
        
    RETURN:

    if (status == FAILURE)
    {
        DR_Error("VolDiag_Check_W: failed!");
    }
        
    return (status);

}  /* VolDiag_Check_W */

/*********************************************************************************
 *       Check validity of credit volatility inputs.
 *
 ********************************************************************************/
int     VolDiagCR_Check_W (
    long    BaseDate,                     /* Volatility data                    */
    int     VolUnit,
    int     NbVol,
    long    *VolDate,
    long    *SwapSt,
    long    *SwapMat,
    double  *Vol)    
{
    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];

    /* Basic Check */
    if ((VolDate == NULL) ||
        (SwapSt  == NULL) ||
        (SwapMat == NULL) ||
        (Vol     == NULL)) goto RETURN;


    if (Dateok(BaseDate))
    {
        sprintf (ErrorMsg, "Incorrect format for start date! (VolDiagCR_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((VolUnit != 0) &&
        (VolUnit != 1))
    {
        sprintf (ErrorMsg, "Volatility unit must be 0 or 1! (VolDiagCR_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (VolUnit != 0)
    {
        sprintf (ErrorMsg, "Only 0 is supported for volatility unit! (VolDiagCR_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((NbVol > MAX_VOL) ||
        (NbVol < 0))
    {        
        sprintf (ErrorMsg, "Nb of vols must be >= 0 and <= %d! (VolDiagCR_Check_W)", MAX_VOL);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    for (i=0; i<NbVol; i++)
    {
        if ((Dateok(VolDate[i]))||
            (Dateok(SwapSt [i]))||
            (Dateok(SwapMat[i]))) 
        {
            sprintf (ErrorMsg, "incorrect date format #%d! (VolDiagCR_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (SwapMat[i] <= SwapSt[i])
        {
            sprintf (ErrorMsg, "rate end <= rate start #%d! (VolDiagCR_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (VolDate[i] >  SwapSt[i])
        {
            sprintf (ErrorMsg, "expiration date > rate start #%d! (VolDiagCR_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (Vol[i] < 0.0)
        {
            sprintf (ErrorMsg, "vol %d < 0! (VolDiagCR_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }


    for (i = 1; i < NbVol; i++)
    {
        if (VolDate[i] <= VolDate[i-1])
        {
            sprintf (ErrorMsg, "Volatilty dates must be entered in ascending order! (VolDiagCR_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (NbVol != 0)
    {
        if (BaseDate >= VolDate[0])
        {
            sprintf (ErrorMsg, "Base date must be before vol date! (VolDiag_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    status = SUCCESS;
        
    RETURN:

    if (status == FAILURE)
    {
        DR_Error("VolDiagCR_Check_W: failed!");
    }
        
    return (status);

}  /* VolDiagCR_Check_W */


/*********************************************************************************
 *  	Read model parameters
 *
 *       Alpha, Beta and Rho must have enough memory for MaxNbFactor before 
 *       calling this function
 *
 ********************************************************************************/
int     Param_Input (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName,                 /* (I) File name including extension  */
    int        MAWFlg)
{
    int status;

    if (MAWFlg == TRUE)
        status = Param_Input_MAW(
                        NbFactor,
                        Alpha,
                        Beta,
                        Rho,
                        BackboneCoeff,
                        QLeft,
                        QRight,
                        FwdShift,
                        MaxNbFactor,
                        FileName);
    else
        status = Param_Input_Classic(
                        NbFactor,
                        Alpha,
                        Beta,
                        Rho,
                        BackboneCoeff,
                        QLeft,
                        QRight,
                        FwdShift,
                        MaxNbFactor,
                        FileName);

    if (status == FAILURE)
    {
        DR_Error("Param_Input: failed!");
    }


    return (status);

}

int     Param_Input_Classic (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName)                 /* (I) File name including extension  */
{
    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;
    int     NbRhoNeeded, MaxNbRho;

    /* Basic Checks */
    if ((NbFactor       == NULL) ||
        (BackboneCoeff  == NULL) ||
        (QLeft          == NULL) ||
        (QRight         == NULL) ||
        (FwdShift       == NULL) ||
        (FileName       == NULL)) goto RETURN;

    if ((MaxNbFactor < 0L) ||
        (MaxNbFactor > 3L)) goto RETURN;

    if (MaxNbFactor != 0L)
    {
        if ((Alpha          == NULL) ||
            (Beta           == NULL)) goto RETURN;
    }

    if (MaxNbFactor > 1L)
    {
        if (Rho == NULL) goto RETURN;
    }


    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Param_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if (FindAndSkipComLine (1, stream, "Number of factors", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

	readerror = fscanf (stream, "%d \n", NbFactor);

    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read Number of factors! (Param_Input)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((*NbFactor < 0) ||
        (*NbFactor > MaxNbFactor))
    {
        DR_Error("Nb of factors out of range! (Param_Input)");
        goto RETURN;
    }

    /* assign dummy numbers into useless array spaces */
    for (i= *NbFactor; i<MaxNbFactor; i++)
    {
        Alpha[i]        = -99999.;
        Beta[i]         = -99999.;
    }

    NbRhoNeeded = (*NbFactor) * ((*NbFactor) - 1) / 2;
    MaxNbRho    = MaxNbFactor * (MaxNbFactor - 1) /2;

    for (i= NbRhoNeeded; i<MaxNbRho; i++)
    {
        Rho[i]          = -99999.;
    }

    if (*NbFactor == 0)
    {
        *BackboneCoeff   = -99999.;
        *QLeft           = -99999.;
        *QRight          = -99999.;
        *FwdShift        = -99999.;
        
        status = SUCCESS;
        goto RETURN;
    }


    if (FindAndSkipComLine (1, stream, "mean reversion", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i=0; i<*NbFactor; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Beta[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read mean reversion in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (FindAndSkipComLine (1, stream, "factor weights", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i=0; i<*NbFactor; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Alpha[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read factor weights in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (FindAndSkipComLine (1, stream, "factor correlations", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i=0; i < NbRhoNeeded; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Rho[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read factor correlations in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (FindAndSkipComLine (1, stream, "backbone coefficient", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf \n", BackboneCoeff);

    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read backbone coefficient in file %s! (Param_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, stream, "Smile parameters", "Param_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf %lf %lf \n", 
                        &(*QLeft),
                        &(*QRight),
                        &(*FwdShift));

    if (readerror != 3)
    {      
        sprintf (ErrorMsg, "Could not read smile parameters in file %s! (Param_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    *QLeft  = 1. - *QLeft;
    *QRight = 1. - *QRight;


    status = SUCCESS;

    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }


    if (status == FAILURE)
    {
        DR_Error("Param_Input: failed!");
    }


    return (status);

}  /* Param_Input_Classic */

/*********************************************************************************
 *  	Read model parameters
 *
 *       Alpha, Beta and Rho must have enough memory for MaxNbFactor before 
 *       calling this function
 *
 ********************************************************************************/
int     Param_Input_MAW (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName)                 /* (I) File name including extension  */
{
    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;
    int     NbRhoNeeded, MaxNbRho;
    int     NbTDInp;
    long    TDDate;


    /* Basic Checks */
    if ((NbFactor       == NULL) ||
        (BackboneCoeff  == NULL) ||
        (QLeft          == NULL) ||
        (QRight         == NULL) ||
        (FwdShift       == NULL) ||
        (FileName       == NULL)) goto RETURN;

    if ((MaxNbFactor < 0L) ||
        (MaxNbFactor > 3L)) goto RETURN;

    if (MaxNbFactor != 0L)
    {
        if ((Alpha          == NULL) ||
            (Beta           == NULL)) goto RETURN;
    }

    if (MaxNbFactor > 1L)
    {
        if (Rho == NULL) goto RETURN;
    }


    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, 
                "Could not open file %s! (Param_Input_MAW)", 
                FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if (FindAndSkipComLine (1, 
                            stream, 
                            "Title line", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine( 1, 
                                stream, 
                                "Section: VNFM parameters", 
                                "Param_Input_MAW", 
                                FileName) == FAILURE)
    {        
        goto RETURN;
    }
    if (FindAndSkipComLine (1, 
                            stream, 
                            "Number of factors", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	readerror = fscanf (stream, "%d \n", NbFactor);

    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read Number of factors! (Param_Input)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((*NbFactor < 0) ||
        (*NbFactor > MaxNbFactor))
    {
        DR_Error("Nb of factors out of range! (Param_Input_MAW)");
        goto RETURN;
    }

    /* assign dummy numbers into useless array spaces */
    for (i= *NbFactor; i<MaxNbFactor; i++)
    {
        Alpha[i]        = -99999.;
        Beta[i]         = -99999.;
    }

    NbRhoNeeded = (*NbFactor) * ((*NbFactor) - 1) / 2;
    MaxNbRho    = MaxNbFactor * (MaxNbFactor - 1) /2;

    for (i= NbRhoNeeded; i<MaxNbRho; i++)
    {
        Rho[i]          = -99999.;
    }

    if (*NbFactor == 0)
    {
        *BackboneCoeff   = -99999.;
        *QLeft           = -99999.;
        *QRight          = -99999.;
        *FwdShift        = -99999.;
        
        status = SUCCESS;
        goto RETURN;
    }


    if (FindAndSkipComLine (1, 
                            stream, 
                            "Number of benchmark dates", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d ", &NbTDInp);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                "Could not read number of time dependent input dates %s! "
                "(Param_Input_MAW)", 
                FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* if classic version of the tree one and only one time dep input date */
    if (NbTDInp != 1)
    {
        sprintf(ErrorMsg, 
                "Only one time dependent input date should be specified for "
                "classic version in file %s! "
                "(Param_Input_MAW)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    if (FindAndSkipComLine (1, 
                            stream, 
                            "Term structure, MR, weights, correlation", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%ld ", &TDDate);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                "Could not read number of time dependent dates %s! "
                "(Param_Input_MAW)", 
                FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    for (i=0; i<*NbFactor; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Beta[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read mean reversion in file %s! "
                    "(Param_Input_MAW)", 
                    FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    for (i=0; i<*NbFactor; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Alpha[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read factor weights in file %s! "
                    "(Param_Input_MAW)", 
                    FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }


    for (i=0; i < NbRhoNeeded; i++)
    {
        readerror = fscanf (stream, "%lf ", &(Rho[i]));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read factor correlations in file %s! "
                    "(Param_Input_MAW)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (FindAndSkipSectionLine( 0, 
                                stream, 
                                "End Section", 
                                "Param_Input_MAW", 
                                FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine( 1, 
                                stream, 
                                "Section: Backbone", 
                                "Param_Input_MAW", 
                                FileName) == FAILURE)
    {        
        goto RETURN;
    }


    if (FindAndSkipComLine (1, 
                            stream, 
                            "backbone coefficient", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf \n", BackboneCoeff);

    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                "Could not read backbone coefficient in file %s! "
                "(Param_Input_MAW)", 
                FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipSectionLine( 0, 
                                stream, 
                                "End Section", 
                                "Param_Input_MAW", 
                                FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine( 1, 
                                stream, 
                                "Section: Smile", 
                                "Param_Input_MAW", 
                                FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipComLine (1, 
                            stream, 
                            "Nb of Smile parameters", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld ", &NbTDInp);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                "Could not read number of time dependent dates %s! "
                "(Param_Input_MAW)", 
                FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (1, 
                            stream, 
                            "Smile parameters", 
                            "Param_Input_MAW", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld %lf %lf %lf \n", 
                        &TDDate,
                        &(*QLeft),
                        &(*QRight),
                        &(*FwdShift));

    if (readerror != 4)
    {      
        sprintf (ErrorMsg, "Could not read smile parameters in file %s! (Param_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    *QLeft  = 1. - *QLeft;
    *QRight = 1. - *QRight;


    status = SUCCESS;

    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }


    if (status == FAILURE)
    {
        DR_Error("Param_Input: failed!");
    }


    return (status);

}  /* MAW_Param_Input */

/*********************************************************************************
 *  	Check model parameters
 *
 ********************************************************************************/
int     Param_Check (long       NbFactor,   
                     double     *Alpha,      
                     double     *Beta,      
                     double     *Rho,       
                     double     BackboneCoeff,
                     double     FwdShift)
{
    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     NbRho;

    /* Basic checks */
    if ((NbFactor < 0L) ||
        (NbFactor > 3L)) goto RETURN;

    if (NbFactor != 0L)
    {
        if ((Alpha          == NULL) ||
            (Beta           == NULL)) goto RETURN;
    }

    if (NbFactor > 1L)
    {
        if (Rho == NULL) goto RETURN;
    }


    /* Check Alpha, Beta and Rho */
    if (NbFactor != 0)
    {
        for (i = 0; i < NbFactor; i++)
        {
            if (Alpha[i] < 0.0)
            {
                DR_Error("Weight out of range! (Param_Check)");
                goto RETURN;
            }

            if (Beta[i] < 0.0)
            {
                DR_Error("Beta out of range! (Param_Check)");
                goto RETURN;
            }
        }

        NbRho = (NbFactor - 1) * NbFactor / 2;

        for (i=0; i<NbRho; i++)
        {
            if ((Rho[i] < -1.0) || (Rho[i] > 1.0))
            {
                DR_Error("Correlation out of range! (Param_Check)");
                goto RETURN;
            }
        }

        if (IS_EQUAL(FwdShift, -1.0))
        {
            DR_Error("Fwd shift = -1.0. (Param_Check)");
            goto RETURN;
        }
            
        if ((BackboneCoeff > 1.0) ||
            (BackboneCoeff < 0.0))
        {
            DR_Error ("Backbone Coeff out of range! (Param_Check)");
            goto RETURN;
        }
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error("Param_Check: failed!");
    }
        
    return (status);

}  /* Param_Check */

/*********************************************************************************
 *       This function does:
 *       1) Mode 0: Find and skip mode
 *           Skip everything until it finds line start with '###' (ignore spaces)
 *           Flag an error message if it can't find
 *       2) Mode 1: Skip mode
 *           Skip all the empty lines and spaces
 *           Read the first non-empty character and check it is '#' followed by '##' 
 *           Flag an error message if it isn't.
 *       3) Mode 2: Find and Skip mode with no error message printed
 *           Same as Mode 0, but no error message is printed if failed to find
 *           next section line.
 *           Note the status is still FAILURE in this case
 *
 *       The routine returns with the file pointer at the start of the next line
 *
 ********************************************************************************/
int     FindAndSkipSectionLine (int      Mode,
                                FILE    *stream,
                                char     *SectionLabel,
                                char     *Routine,
                                char     *FileName)
{
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i;
    int     length;

    if (stream == NULL) goto RETURN;

    if ((Mode != 0) && (Mode != 1) && (Mode != 2))
    {
        DR_Error("FindAndSkipSectionLine: Mode has to be 0 or 1 or 2!");
        goto RETURN;
    }

    while (fgets (string, MAXBUFF, stream) != NULL)
    {
        length = strlen(string);

        for (i = 0; i < length; i++)
        {
            if (string[i] == '#')
            {
                if ((i > length - 3) || 
                    (string[i+1] != '#') || 
                    (string[i+2] != '#')) 
                {
                    if (Mode == 0) break;
                    else goto RETURN;
                }
                else
                {
                    status = SUCCESS;
                    goto RETURN;
                }

            }
            else if ((string[i] != ' ' ) && 
                     (string[i] != '\n') && 
                     (string[i] != '\t')) 
            {
                if ((Mode == 0) || (Mode == 2)) break;
                else goto RETURN;
            }
        }
    }  

    RETURN:


    if (status == FAILURE)
    {
        if (Mode != 2)
        {
            DR_Error ("%s: %s section line expected and not found in file"
                      " %s!", Routine, SectionLabel, FileName);
        }
    }
        
    return (status);

}  /* FindAndSkipSectionLine */



/*********************************************************************************
 *      Read summary.dat
 *
 *      Check all inputs  
 *
 ********************************************************************************/
int    summaryInputEnvData(
    char     *FNameSummary,               /* (I) Summary filename               */
    long     *Today,                      /* (O)                                */
    int      *NbIrInp,                    /* (O)                                */
    IR_INPUT **IrInp,                     /* (O)                                */
    long     *NbSpInp,                    /* (O)                                */
    SP_INPUT **SpInp,                     /* (O)                                */
    long     *NbFxInp,                    /* (O)                                */
    FX_INPUT **FxInp,                     /* (O)                                */
    long     *NbEqInp,                    /* (O)                                */
    EQ_INPUT **EqInp,                     /* (O)                                */
    int      *NbCrInp,                    /* (O)                                */
    CR_INPUT **CrInp)                     /* (O)                                */
{
    int    status = FAILURE;

    int    readerror;
    char   ErrorMsg[MAXBUFF];

    long   i;

    FILE  *stream = NULL;

    long NbIrInpL = 0,
         NbSpInpL = 0,
         NbFxInpL = 0,
         NbEqInpL = 0,
         NbCrInpL = 0;

    long TodayL;

    IR_INPUT *IrInpL = NULL;
    SP_INPUT *SpInpL = NULL;
    FX_INPUT *FxInpL = NULL;
    EQ_INPUT *EqInpL = NULL;
    CR_INPUT *CrInpL = NULL;

    /* basic check */
    if ((FNameSummary == NULL) ||
        (Today        == NULL)) goto RETURN;

    stream = fopen (FNameSummary, "r");
    if (stream == NULL)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: Cannot open file %s!", FNameSummary);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((NbIrInp == NULL) || (IrInp == NULL)) goto RETURN;
    if ((NbSpInp == NULL) != (SpInp == NULL)) goto RETURN;
    if ((NbFxInp == NULL) != (FxInp == NULL)) goto RETURN;
    if ((NbEqInp == NULL) != (EqInp == NULL)) goto RETURN;
    if ((NbCrInp == NULL) != (CrInp == NULL)) goto RETURN;

    /* --------------------- */
    /*  ENVIRONMENT SECTION  */
    /* --------------------- */

    if (FindAndSkipSectionLine(1, /* skip mode */
                               stream, 
                               "Env Section start",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;

    /* read Today */
    if (FindAndSkipComLine (1, stream, 
                           "Today",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(TodayL));
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read today's date.");
        goto RETURN;
    }

    /* check Today */
    if (Dateok(TodayL))
    {
        DR_Error("summaryInputEnvData: incorrect today's date format!");
        goto RETURN;
    }

    if (FindAndSkipSectionLine( 0,  /* skip mode */
                                stream, 
                                "Env Section end",
                                "summaryInputEnvData", 
                                FNameSummary
                                ) == FAILURE) goto RETURN;

    /* ------------ */
    /*  IR SECTION  */
    /* ------------ */

    if (FindAndSkipSectionLine(1,  /* skip mode */
                               stream, 
                               "IR Section start",
                               "summaryInputEnvData", 
                               FNameSummary
                               ) == FAILURE) goto RETURN;

    /* read Nb IR */
    if (FindAndSkipComLine (1, stream, 
                           "Nb IR start",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &NbIrInpL);
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read Nb IR.");
        goto RETURN;
    }

    if (NbIrInpL < 1)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: there must be at least one IR!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* allocate memory for IR */
    IrInpL = calloc(NbIrInpL, sizeof(IR_INPUT));
    if (IrInpL == NULL) goto RETURN;
   
    /* skip to next section */
    if (FindAndSkipSectionLine(0,  /* Find and skip mode */
                               stream, 
                               "IR Section end",
                               "summaryInputEnvData", 
                               FNameSummary
                               ) == FAILURE) goto RETURN;

    /* ------------ */
    /*  SP SECTION  */
    /* ------------ */

    if (FindAndSkipSectionLine(1,  /* Find and skip mode */
                            stream, 
                            "SP Section start",
                            "summaryInputEnvData", 
                            FNameSummary) == FAILURE) goto RETURN;

    /* read Nb SP */
    if (FindAndSkipComLine (1, stream, 
                           "Nb SP",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(NbSpInpL));
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read Nb SP.");
        goto RETURN;
    }

    if (NbSpInpL < 0)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: negative nb of SP!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (NbSpInpL > 0)
    {
        char OverWriteString[MAXBUFF];
        char *getserror;
        /* allocate memory for SP */
        SpInpL = calloc(NbSpInpL, sizeof(SP_INPUT));
        if (SpInpL == NULL) goto RETURN;

        /* read SP currency denomination */
        if (FindAndSkipComLine (1, stream, 
                               "SP currency denomination",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;

        for (i = 0; i < NbSpInpL; i++)
        {
            getserror = fgets(OverWriteString, MAXBUFF, stream);
            readerror = sscanf (OverWriteString, "%ld \t%ld \t %ld\n", 
                                &(SpInpL[i].BaseIR_id),
                                &(SpInpL[i].LiborCrv),
                                &(SpInpL[i].DiscCrv));
            if ( readerror != 3)
            {
                readerror = sscanf (
                                 OverWriteString, "%ld \t%ld\n", 
                                &(SpInpL[i].BaseIR_id),
                                &(SpInpL[i].LiborCrv));
                if (readerror != 2)
                {
                    DR_Error("summaryInputEnvData: Cannot read Base IR id "
                             "/ reference curve.");
                    goto RETURN;
                }
                SpInpL[i].DiscCrv = SpInpL[i].LiborCrv;
            }
        }
    }
    
    if (FindAndSkipSectionLine(0, stream, 
                               "SP Section end",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;


    /* ------------ */
    /*  FX SECTION  */
    /* ------------ */

    if (FindAndSkipSectionLine(1, stream, 
                        "FX Section start",
                        "summaryInputEnvData", 
                        FNameSummary) == FAILURE) goto RETURN;

    /* read Nb FX */
    if (FindAndSkipComLine (1, stream, 
                           "Nb FX",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(NbFxInpL));
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read Nb FX.");
        goto RETURN;
    }

    if (NbFxInpL < 0)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: negative nb of FX!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (NbFxInpL > 0)
    {
        /* allocate memory for FX */
        FxInpL = calloc(NbFxInpL, sizeof(FX_INPUT));
        if (FxInpL == NULL) goto RETURN;

        /* read foreign IR id */
        if (FindAndSkipComLine (1, stream, 
                               "Foreign IR id",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;

        for (i=0; i<NbFxInpL; i++)
        {
            readerror = fscanf (stream, "%ld \n", 
                                &(FxInpL[i].ForeignIR_id));
            if (readerror != 1)
            {
                DR_Error("summaryInputEnvData: Cannot read foreign IR id.");
                goto RETURN;
            }
        }
    }

    if (FindAndSkipSectionLine(0, stream, 
                              "FX Section end",
                              "summaryInputEnvData", 
                              FNameSummary) == FAILURE) goto RETURN;


    /* ------------ */
    /*  EQ SECTION  */
    /* ------------ */

    if (FindAndSkipSectionLine(1, stream, 
                        "EQ Section start",
                        "summaryInputEnvData", 
                        FNameSummary) == FAILURE) goto RETURN;

    /* read Nb EQ */
    if (FindAndSkipComLine (1, stream, 
                           "Nb EQ",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(NbEqInpL));
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read Nb EQ.");
        goto RETURN;
    }

    if (NbEqInpL < 0)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: negative nb of EQ!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (NbEqInpL > 0)
    {
        /* allocate memory for EQ */
        EqInpL = calloc(NbEqInpL, sizeof(EQ_INPUT));
        if (EqInpL == NULL) goto RETURN;

        /* read Equity base ir id */
        if (FindAndSkipComLine (1, stream, 
                               "Equity BaseIR id",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;

        for (i=0; i<NbEqInpL; i++)
        {
            readerror = fscanf (stream, "%ld \n", 
                                &(EqInpL[i].BaseIR_id));
            if (readerror != 1)
            {
                DR_Error("summaryInputEnvData: Cannot read Equity BaseIR id.");
                goto RETURN;
            }
        }
    }

    if (FindAndSkipSectionLine(0, stream, 
                               "EQ Section end",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;


    /* ------------ */
    /*  CR SECTION  */
    /* ------------ */

    if (FindAndSkipSectionLine(1, stream, 
                        "CR Section start",
                        "summaryInputEnvData", 
                        FNameSummary) == FAILURE) goto RETURN;

    /* read Nb CR */
    if (FindAndSkipComLine (1, stream, 
                           "Nb CR",
                           "summaryInputEnvData", 
                           FNameSummary) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(NbCrInpL));
    if (readerror != 1)
    {
        DR_Error("summaryInputEnvData: Cannot read Nb CR.");
        goto RETURN;
    }

    if (NbCrInpL < 0)
    {
        sprintf (ErrorMsg, "summaryInputEnvData: negative nb of CR!");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (NbCrInpL > 0)
    {
        /* allocate memory for CR */
        CrInpL = calloc(NbCrInpL, sizeof(CR_INPUT));
        if (CrInpL == NULL) goto RETURN;

        /* read base ir id for credit */
        if (FindAndSkipComLine (1, stream, 
                               "Credit BaseIR id",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;

        for (i=0; i<NbCrInpL; i++)
        {
            readerror = fscanf (stream, "%ld \n", 
                                &(CrInpL[i].BaseIR_id));
            if (readerror != 1)
            {
                DR_Error("summaryInputEnvData: Cannot read Credit BaseIR id.");
                goto RETURN;
            }
        }
    }
    
    if (FindAndSkipSectionLine(0, stream, 
                               "CR Section end",
                               "summaryInputEnvData", 
                               FNameSummary) == FAILURE) goto RETURN;



    /* Checks */
    if ((NbSpInpL > 0) && (NbSpInp == NULL)) goto RETURN;
    if ((NbFxInpL > 0) && (NbFxInp == NULL)) goto RETURN;
    if ((NbEqInpL > 0) && (NbEqInp == NULL)) goto RETURN;
    if ((NbCrInpL > 0) && (NbCrInp == NULL)) goto RETURN;

    /* IR outputs */
    *NbIrInp    = NbIrInpL;
    *IrInp      = IrInpL;

    /* SP outputs */
    if (NbSpInp != NULL)
    {
        *NbSpInp    = NbSpInpL;
        *SpInp      = SpInpL;
    }

    /* FX outputs */
    if (NbFxInp != NULL)
    {
        *NbFxInp    = NbFxInpL;
        *FxInp      = FxInpL;
    }

    /* EQ outputs */
    if (NbEqInp != NULL)
    {
        *NbEqInp    = NbEqInpL;
        *EqInp      = EqInpL;
    }

    /* CR outputs */
    if (NbCrInp != NULL)
    {
        *NbCrInp    = NbCrInpL;
        *CrInp      = CrInpL;
    }

    *Today = TodayL;

    status = SUCCESS;

RETURN:

    if (stream != NULL) fclose(stream);

    if (status == FAILURE)
    {
        if (IrInpL != NULL) free(IrInpL);
        if (FxInpL != NULL) free(FxInpL);
        if (SpInpL != NULL) free(SpInpL);
        if (EqInpL != NULL) free(EqInpL);
        if (CrInpL != NULL) free(CrInpL);

        DR_Error("summaryInputEnvData: failed!");
    }

    return (status);

} /* summaryInputEnvData */

/*----------------------IRINPUT-------------------------------------------------*/
/*********************************************************************************
 *    irInput
 *    Read the following files
 *      1) the 3 zero curves
 *      2) IR info file
 *      3) IR volatility file
 *      4) IR parameter file
 *
 ********************************************************************************/
int  irInputEnvData (
                IR_INPUT  *irInput,       /* (O) ir Input Structure             */
                char      *FNameZC0,      /* (I) file name for zero curve 0     */
                char      *FNameZC1,      /* (I) file name for zero curve 1     */
                char      *FNameZC2,      /* (I) file name for zero curve 2     */
                char      *FNameInfo,     /* (I) IR info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS)      /* (I) 2q string                      */
{
    int    status       = FAILURE;
    FILE   *stream      = NULL;
    int    readerror;
    int    mawFlag      = FALSE;
    char   ErrorMsg[MAXBUFF];

    /* basic check */
    if ((irInput      == NULL) ||
        (FNameZC0     == NULL) ||
        (FNameZC1     == NULL) ||
        (FNameZC2     == NULL) ||
        (FNameInfo    == NULL) ||
        (FNameVoldiag == NULL) ||
        (FNameModlPar == NULL) ||
        (CalibFlag    == NULL) ||
        (NbFactorOWS  == NULL) ||
        (AlphaOWS     == NULL) ||
        (BetaOWS      == NULL) ||
        (RhoOWS       == NULL) ||
        (BackboneOWS  == NULL) ||
        (SmileOWS     == NULL)) goto RETURN;

    if (strcmp(CalibFlag, "Y"   ) && 
        strcmp(CalibFlag, "Y*"  ) &&
	strcmp(CalibFlag, "N"   ) )
    {
        DR_Error ("Calibration flag must be 'Y' or 'Y*' or 'N'!(irInputEnvData)");
        goto RETURN;
    }

    stream = fopen("model.dat","r");
    if (stream != NULL)
    {
        mawFlag = TRUE;
        fclose(stream);
    }

    /**********************/
    /*  Read Zero Curves  */
    /**********************/

    /* Read index curve */
    if (Term_Input_W (&(irInput->ZeroCurve[0]), FNameZC0) == FAILURE)
    {
        goto RETURN;
    }

    /* Read discount curve */
    if (Term_Input_W (&(irInput->ZeroCurve[1]), FNameZC1) == FAILURE)
    {
        goto RETURN;
    }

    /* Read riskzero curve, if failed use discount */
    {
    FILE   *streamRisk  = NULL;
    streamRisk = fopen (FNameZC2, "r");
    if (streamRisk == NULL)
    {
        if (Term_Input_W (&(irInput->ZeroCurve[2]), FNameZC1) == FAILURE)
        {
            goto RETURN;
        }
    }
    else
    {
        fclose(streamRisk);
        if (Term_Input_W (&(irInput->ZeroCurve[2]), FNameZC2) == FAILURE)
        {
            goto RETURN;
        }
    }
    }

    /********************/
    /*  Read IR info    */
    /********************/
    {
        char    SwapDCC[MAXBUFF];
        int     MMBinput;

        /* Open and read the IR info file */
        stream = fopen (FNameInfo, "r");

        if (stream == NULL)
        {
            sprintf (ErrorMsg, "Could not open file %s! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Money Market Basis */
        if (FindAndSkipComLine (1, stream, "money market basis", "irInputEnvData", FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%d \n", &(MMBinput));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read money market basis in %s! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (MMBinput == 360) irInput->MMB = '0';
        else if (MMBinput == 365) irInput->MMB = '5';
        else
        {
            sprintf (ErrorMsg, "MMB in %s has to be 360 or 365! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Swap Day Count Convention */
        if (FindAndSkipComLine (1, stream, "swap day count convention", "irInputEnvData", FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s \n", SwapDCC);
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read swap dcc in file %s! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                
        
        if (strstr(SwapDCC, "30/360") != NULL)
        {
            irInput->SwapDCC = '3'; 
        }
        else if (strstr(SwapDCC, "ACT/365F") != NULL)
        {
            irInput->SwapDCC = '5';
        }
        else if (strstr(SwapDCC, "ACT/360") != NULL)
        {
            irInput->SwapDCC = '0';
        }
        else
        {
            sprintf (ErrorMsg, "swap dcc in file %s has to be '30/360' or"
                    " 'ACT/365F' or 'ACT/360'! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read swap frequency */
        if (FindAndSkipComLine (1, stream, "swap frequency", "irInputEnvData", FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%c \n", &(irInput->SwapFreq));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read swap frequency in %s! (irInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }

    /**********************/
    /*  Read Vol diagonal */
    /**********************/
    if (strchr(CalibFlag, 'Y') != NULL)
    {   
	  irInput->NbVol = 3;
        if (VolDiag_Input_W (
                        &(irInput->BaseDate),
                        &(irInput->VolUnit),
                        &(irInput->NbVol),
                        irInput->VolDate,
                        irInput->SwapSt,
                        irInput->SwapMat,
                        irInput->Vol,
                        irInput->VolType,
                        FNameVoldiag) == FAILURE) goto RETURN;

        /* Calib is off either CalibFlag from deal = 'N' or NbVol from env = 0 */
        if (irInput->NbVol < 0)
        {
            sprintf (ErrorMsg, "Nb Vol has to be >= 0! (irInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        else if (irInput->NbVol == 0)
        {
            irInput->CalibFlag = FALSE;
            irInput->SkipFlag = FALSE;
            irInput->VolUnit = -99999;
            irInput->NbVol = -99999;
            irInput->Freq = 'z';
            irInput->DCC = 'z';
        }
        else 
        {
            irInput->CalibFlag = TRUE;

            /* set Freq and DCC for vols */
            /* only need to check the first 'Vol Type' entry because 
               in irInputCheck we will check that all Vol Types are equal */
            if (irInput->VolType[0] == 'C')
            {
                irInput->Freq = 'I'; 
                irInput->DCC = irInput->MMB;
            }
            else 
            {
                irInput->Freq = irInput->SwapFreq;
                irInput->DCC = irInput->SwapDCC;
            }

            /* set Skip Flag */
            if (strchr(CalibFlag, '*') == NULL)
            {
                irInput->SkipFlag = FALSE;  /* vol points skipping not allowed */
            }
            else
            {
                irInput->SkipFlag = TRUE;   /* vol points skipping allowed      */
            }

        }
    }
    else
    {
        irInput->CalibFlag = FALSE;
        irInput->SkipFlag = FALSE;
        irInput->BaseDate = -99999;
        irInput->VolUnit = -99999;
        irInput->NbVol = -99999;
        irInput->Freq = 'z';
        irInput->DCC = 'z';
    }


    /*************************/
    /*  Read parameter file  */
    /*************************/

    {
        int     isFileNeeded; /* TRUE = we need some data from text file */
        long    NbFactorInput = 0L;

        /* Check if 'NbFactorOWS' is not given, all related overwrites are also not given */
        if (strstr(NbFactorOWS, "nil") != NULL)
        {
            if ((strstr(AlphaOWS, "nil") == NULL) ||
                (strstr(BetaOWS,  "nil") == NULL) ||
                (strstr(RhoOWS,   "nil") == NULL))
            {
                DR_Error ("Need nb of assets overwrite in order to overwrite "
                          "factor parameters (irInputEnvData)");
                goto RETURN;
            }
        }

        /* if not all overwrites are given, read IR param file */
        if (strstr (NbFactorOWS, "nil") == NULL) 
        {

            readerror = sscanf (NbFactorOWS, "%ld \n", &(NbFactorInput));
            if (readerror != 1)
            {        
                DR_Error("irInputEnvData: cannot read nb of factor.");
                goto RETURN;
            }

            if (NbFactorInput == 1L)
            {
                isFileNeeded = ((strstr(AlphaOWS,      "nil") != NULL) ||
                                (strstr(BetaOWS,       "nil") != NULL) ||
                                (strstr(BackboneOWS,   "nil") != NULL) ||
                                (strstr(SmileOWS,      "nil") != NULL));
            }
            else if (NbFactorInput > 1L)
            {
                isFileNeeded = ((strstr(AlphaOWS,      "nil") != NULL) ||
                                (strstr(BetaOWS,       "nil") != NULL) ||
                                (strstr(RhoOWS,        "nil") != NULL) ||
                                (strstr(BackboneOWS,   "nil") != NULL) ||
                                (strstr(SmileOWS,      "nil") != NULL));
            }
            else 
            {
                DR_Error("irInputEnvData: nb of factor overwrite < 1.");
                goto RETURN;
            }
        }
        else isFileNeeded = TRUE;

        if (isFileNeeded)
        {
            if (Param_Input (&(irInput->NbFactor),
                             irInput->Alpha,
                             irInput->Beta,
                             irInput->Rho,
                             &(irInput->BackboneCoeff),
                             &(irInput->QLeft),
                             &(irInput->QRight),
                             &(irInput->FwdShift),
                             3,
                             FNameModlPar,
                             mawFlag) == FAILURE) goto RETURN;

        }

        /* overwrite nb of factors */
        if (strstr (NbFactorOWS, "nil") == NULL) 
        {
            if (isFileNeeded)
            {
                if (irInput->NbFactor != NbFactorInput)
                {
                    DR_Error("inconsistent nb of factor inputs! (irInputEnvData)");
                    goto RETURN;
                }
            }
            else
            {
                irInput->NbFactor = NbFactorInput;
            }
        }
    
        if ((irInput->NbFactor != 1) &&
            (irInput->NbFactor != 2) &&
            (irInput->NbFactor != 3))
        {
            sprintf (ErrorMsg, "Number of factor must be 1, 2 or 3! (irInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* overwrite alpha, beta and rho */
        if (irInput->NbFactor == 1)
        {
            if (strstr (AlphaOWS, "nil") == NULL) 
            {
                readerror = sscanf (AlphaOWS, 
                                    "%lf \n", 
                                    &(irInput->Alpha[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (strstr (BetaOWS, "nil") == NULL) 
            {
                readerror = sscanf (BetaOWS, 
                                    "%lf \n", 
                                    &(irInput->Beta[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }


            irInput->Alpha[1]        = -99999.;
            irInput->Alpha[2]        = -99999.;
            irInput->Beta[1]         = -99999.;
            irInput->Beta[2]         = -99999.;
            irInput->Rho[0]          = -99999.;
            irInput->Rho[1]          = -99999.;
            irInput->Rho[2]          = -99999.;
        }
        else if (irInput->NbFactor == 2)
        {
            if (strstr (AlphaOWS, "nil") == NULL) 
            {
                readerror = sscanf (AlphaOWS, 
                                    "%lf \t%lf \n", 
                                    &(irInput->Alpha[0]),
                                    &(irInput->Alpha[1]));

                if (readerror != 2)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (strstr (BetaOWS, "nil") == NULL) 
            {
                readerror = sscanf (BetaOWS, 
                                    "%lf \t%lf \n", 
                                    &(irInput->Beta[0]),
                                    &(irInput->Beta[1]));

                if (readerror != 2)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (strstr (RhoOWS, "nil") == NULL) 
            {
                readerror = sscanf (RhoOWS, 
                                    "%lf \n", 
                                    &(irInput->Rho[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor correlation! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            irInput->Alpha[2]        = -99999.;
            irInput->Beta[2]         = -99999.;
            irInput->Rho[1]          = -99999.;
            irInput->Rho[2]          = -99999.;
        }
        else if (irInput->NbFactor == 3)
        {
            if (strstr (AlphaOWS, "nil") == NULL) 
            {
                readerror = sscanf (AlphaOWS, 
                                    "%lf \t%lf \t%lf \n", 
                                    &(irInput->Alpha[0]),
                                    &(irInput->Alpha[1]),
                                    &(irInput->Alpha[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (strstr (BetaOWS, "nil") == NULL) 
            {
                readerror = sscanf (BetaOWS, 
                                    "%lf \t%lf \t%lf \n", 
                                    &(irInput->Beta[0]),
                                    &(irInput->Beta[1]),
                                    &(irInput->Beta[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (strstr (RhoOWS, "nil") == NULL) 
            {
                readerror = sscanf (RhoOWS, 
                                    "%lf \t%lf \t%lf \n", 
                                    &(irInput->Rho[0]),
                                    &(irInput->Rho[1]),
                                    &(irInput->Rho[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor correlation! (irInputEnvData)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }
        }

        /* overwrite backbone coefficient */
        if (strstr (BackboneOWS, "nil") == NULL) 
        {
            readerror = sscanf (BackboneOWS, "%lf \n", &(irInput->BackboneCoeff));
            if (readerror != 1)
            {        
                DR_Error("irInputEnvData: cannot read backbone coefficient.");
                goto RETURN;
            }
        }

        /* overwrite smile parameters */
        if (strstr (SmileOWS, "nil") == NULL) 
        {
            readerror = sscanf (SmileOWS, 
                                "%lf \t%lf \t%lf \n", 
                                &(irInput->QLeft),
                                &(irInput->QRight),
                                &(irInput->FwdShift));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not read overwrite string for 2q! (irInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            irInput->QLeft = 1.0 - irInput->QLeft;
            irInput->QRight = 1.0 - irInput->QRight;
        }
    }

    status = SUCCESS;

RETURN:

    if (stream     != NULL) fclose(stream);
    
    if (status == FAILURE)
    {
        DR_Error("irInputEnvData: failed!");
    }
    return (status);

} /* irInputEnvData */


/*****  irInputCheck  ****************************************************/
/*
 *    Checks the integrity of the irInput structure.
 *    Must be called before Build_Sim
 */
int  irInputCheck (IR_INPUT  *irInput)
{
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];
    long    i;

    if (irInput == NULL) goto RETURN;

    /* Check zero curves */
    for (i=0; i<3; i++)
    {
        if (Term_Check_W(&(irInput->ZeroCurve[i])) == FAILURE) goto RETURN;
    }

    /* check Calib Flag */
    if ((irInput->CalibFlag != TRUE) &&
        (irInput->CalibFlag != FALSE))
    {
        DR_Error("invalid Calib Flag! (irInputCheck)");
        goto RETURN;
    }

    /* check benchmark yield curve conventions */
    if ((irInput->MMB != '0') &&
        (irInput->MMB != '5'))
    {
        sprintf (ErrorMsg, "MM basis has to be '0' or '5'! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if ((irInput->SwapDCC != '0') &&
        (irInput->SwapDCC != '3') &&
        (irInput->SwapDCC != '5'))
    {
        sprintf (ErrorMsg, "swap dcc has to be '0' or"
                " '3' or '5'! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (   (irInput->SwapFreq != 'A')
        && (irInput->SwapFreq != 'S')
    	&& (irInput->SwapFreq != 'Q')
    	&& (irInput->SwapFreq != 'M'))
    {
        sprintf (ErrorMsg, "swap freq has to be 'A' or"
                " 'S' or 'Q' or 'M'! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Check vol inputs */
    if (irInput->CalibFlag == TRUE)
    {

         if (VolDiag_Check_W (irInput->BaseDate,        /* Volatility data */
                              irInput->VolUnit,
                              irInput->NbVol,
                              irInput->VolDate,
                              irInput->SwapSt,
                              irInput->SwapMat,
                              irInput->Vol,
                              irInput->VolType) == FAILURE) goto RETURN;

         /* NbVol = 0 is not allowed for IR */
         if (irInput->NbVol == 0)
         {
            sprintf (ErrorMsg, "Nb vol can't be zero! (irInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
         }
         
         /* Check VolType[i] are all 'C' or all 'S' */
         /* Not RETURN inside VolDiag_Check_W because other assets may
            allow combination of 'C' and 'S' */
         if ((irInput->VolType[0] != 'C') &&
             (irInput->VolType[0] != 'S'))
         {
            sprintf (ErrorMsg, "vol type has to be 'C' or 'S'! (irInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
         }

         for (i=1; i<irInput->NbVol; i++)
         {
             if (irInput->VolType[i] != irInput->VolType[0])
             {
                sprintf (ErrorMsg, "vol types have to be equal! (irInputCheck)");
                DR_Error(ErrorMsg);
                goto RETURN;
             }
         }
     
        if ((irInput->Freq != 'A') &&
            (irInput->Freq != 'S') &&
            (irInput->Freq != 'Q') &&
            (irInput->Freq != 'M') &&
            (irInput->Freq != 'I')) 
        {
            sprintf (ErrorMsg, "vol freq has to be 'A' or"
                    " 'S' or 'Q' or 'M' or 'I'! (irInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* irInput->Freq = 'I' iff irInput->VolType is 'C' */
        if ((irInput->Freq == 'I') !=
            (irInput->VolType[0] == 'C'))
        {
            sprintf (ErrorMsg, "Freq='I' iff VolType='C'! (irInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

             
        if ((irInput->DCC != '0') &&
            (irInput->DCC != '3') &&
            (irInput->DCC != '5')) 
        {
            sprintf (ErrorMsg, "vol dcc has to be '0' or"
                    " '3' or '5'! (irInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }
        
        
    /* Check model parameters */

    if (Param_Check(irInput->NbFactor,
                    irInput->Alpha,
                    irInput->Beta,
                    irInput->Rho,
                    irInput->BackboneCoeff,
                    irInput->FwdShift) == FAILURE) goto RETURN;

    /* NbFactor cannot be zero for IR */
    if (irInput->NbFactor == 0)
    {
        sprintf (ErrorMsg, "Nb of factor can't be zero! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* backbone has to be zero */
    if (IS_NOT_ZERO(irInput->BackboneCoeff))
    {
        sprintf (ErrorMsg, "Only zero backbone is supported! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Check swap defn for correlations */
    if ((Dateok(irInput->CorrSwapSt )) ||
        (Dateok(irInput->CorrSwapMat)))
    {
        sprintf (ErrorMsg, "invalid Corr Swap date format! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (irInput->CorrSwapSt >= irInput->CorrSwapMat)
    {
        sprintf (ErrorMsg, "Corr Swap maturity date before start date! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((irInput->CorrSwapDCC != '0') &&
        (irInput->CorrSwapDCC != '3') &&
        (irInput->CorrSwapDCC != '5'))
    {
        sprintf (ErrorMsg, "Corr Swap dcc must be '0' or '3' or '5'! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if ((irInput->CorrSwapFreq != 'A') &&
        (irInput->CorrSwapFreq != 'S') &&
        (irInput->CorrSwapFreq != 'Q') &&
        (irInput->CorrSwapFreq != 'M')) 
    {
        sprintf (ErrorMsg, "Corr Swap freq has to be 'A' or"
                " 'S' or 'Q' or 'M'! (irInputCheck)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Check Others */
    if ((irInput->CrvToDiff < 0L) || 
        (irInput->CrvToDiff > 2L))
    {
        DR_Error("invalid curve to diffuse choice! (irInputCheck)");
        goto RETURN;
    }

    if ((irInput->CrvToDisc < 0L) || 
        (irInput->CrvToDisc > 2L))
    {
        DR_Error("invalid curve to discount choice! (irInputCheck)");
        goto RETURN;
    }

    if ((irInput->SkipFlag != TRUE) &&
        (irInput->SkipFlag != FALSE))
    {
        DR_Error("invalid Skip Flag! (irInputCheck)");
        goto RETURN;
    }


    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("irInputCheck: failed!");
    }
    return (status);

} /* irInputCheck */

/*----------------------CRINPUT-------------------------------------------------*/
/*****  crInputInit  ******************************************************/
/*
 *   Initialise the Init_CR_Input instance. Must be called before using the
 *   structure. 
 *   This function is introduced to initialize the CalibFlag to FALSE
 */
int  crInputInit(CR_INPUT  *cr_input)
{
    if (cr_input == NULL) return FAILURE;
    memset(cr_input,0,sizeof(*cr_input));  
    return SUCCESS;
} /* crInputInit */

/*********************************************************************************
 *    crInput
 *    Read the following files
 *      1) the credit zero curves
 *      2) CR info file
 *      3) CR volatility file
 *      4) CR parameter file
 * 
 *    The following inputs should be hardcoded when calling this function
 *      CalibFlag       = "N"
 *      NbFactorOWS     = "1"
 *      RhoOWS          = "nil"
 *      BackboneOWS     = "0" (to avoid reading file if all other OWS are given)
 *
 ********************************************************************************/
int  crInputEnvData (
                CR_INPUT  *crInput,       /* (O) cr Input Structure             */
                char      *FNameZC,       /* (I) file name for zero curve 0     */
                char      *FNameInfo,     /* (I) CR info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS)      /* (I) 2q string                      */
{
    int    status     = FAILURE;
    int    readerror;
    char   ErrorMsg[MAXBUFF], Buffer[MAXBUFF];
    FILE    *stream = NULL;
    char SwapDCC [MAXBUFF];

    /* basic check */
    if ((crInput      == NULL) ||
        (FNameZC      == NULL) ||
        (FNameInfo    == NULL) ||
        (FNameVoldiag == NULL) ||
        (FNameModlPar == NULL) ||
        (CalibFlag    == NULL) ||
        (NbFactorOWS  == NULL) ||
        (AlphaOWS     == NULL) ||
        (BetaOWS      == NULL) ||
        (RhoOWS       == NULL) ||
        (BackboneOWS  == NULL) ||
        (SmileOWS     == NULL)) goto RETURN;

    if (strcmp(CalibFlag, "Y" ) && 
        strcmp(CalibFlag, "Y*") &&
        strcmp(CalibFlag, "Y*S" ) &&
		strcmp(CalibFlag, "YS"  ) &&
        strcmp(CalibFlag, "N" ))
    {
        DR_Error ("Calibration flag must be 'Y' or 'Y*' or 'N'!(crInputEnvData)");
        goto RETURN;
    }
    
    /**********************/
    /*  Read Zero Curves  */
    /**********************/

    /* Read index curve */
    if (Term_Input_W (&(crInput->ZeroCurve), FNameZC) == FAILURE)
    {
        goto RETURN;
    }


    /********************/
    /*  Read CR info    */
    /********************/
    {
        /* Open and read the CR info file */
        stream = fopen (FNameInfo, "r");

        if (stream == NULL)
        {
            sprintf (ErrorMsg, "Could not open file %s! (crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Standard Recovery Rate */
        if (FindAndSkipComLine (1, stream, "standard recovery rate", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%lf \n", &(crInput->RecoveryRate));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read recovery rate in %s! (crInputEnvData)", 
                FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read recovery rate convention */
        if (FindAndSkipComLine (1, stream, "recovery rate convention", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%c \n", &(crInput->RecoveryRateConv));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read recovery rate convention in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }           
        
        /* read SPN */
        if (FindAndSkipComLine (1, stream, "SPN", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%ld \n", &(crInput->SPN));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read SPN in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                

        /* read counter party name */
        if (FindAndSkipComLine (1, stream, "counterparty name", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        if( fgets (crInput->CtpyName, MAXBUFF, stream) == NULL )
        {
           sprintf (ErrorMsg, "Could not read counterparty name in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                

    
        /* read currency */
        if (FindAndSkipComLine (1, stream, "currency", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s \n", crInput->Currency);
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read currency in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                


        /* read seniority */
        if (FindAndSkipComLine (1, stream, "senirority", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        if( fgets (crInput->Seniority, MAXBUFF, stream) == NULL )
        {
            sprintf (ErrorMsg, "Could not read seniority in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                


        /* read curve id */
        if (FindAndSkipComLine (1, stream, "curve id", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%ld \n", &(crInput->CurveId));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read curve id in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                

    
        /* read default date */
        if (FindAndSkipComLine (1, stream, "default date", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s\n", Buffer);
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read default date in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                
        if ( strcmp(Buffer, "N") == 0 )
            crInput->DefDate = 0;
        else
        {
            
            readerror = sscanf(Buffer, "%ld \n", &(crInput->DefDate));
            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read default date in file %s! "
                    "(crInputEnvData)", FNameInfo);
                DR_Error(ErrorMsg);
                goto RETURN;
            }                                
        }

        /* read default correlation tweak amount */
        if (FindAndSkipComLine (1, stream, "default corr twk", "crInputEnvData", 
            FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, "%lf \n", &crInput->CorrTwk);
        if (readerror != 1)
        {        
         	sprintf (ErrorMsg, "Could not read default correlation twk amount in file %s! "
                "(crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }             
    }
        

    if (strchr(CalibFlag, 'Y') != NULL)
	{
	    /* reading dcc and cds frequency */
        /* these values are needed for calibration */
        /* currently, CalibFlag is hardcoded in the call 
            of this function in each product manager, where
            this flag can be Y iff DCC, freq and vols
            are present in the environment files */ 

        /* read Swap Day Count Convention */
        if (FindAndSkipComLine (1, stream, "swap day count convention", "crInputEnvData", FNameInfo) == FAILURE)
        {        
            goto RETURN;			
        }
        else
		{

			readerror = fscanf (stream, "%s \n", SwapDCC);

			if (readerror != 1)
			{        
				sprintf (ErrorMsg, "Could not read cds dcc in file %s! (crInputEnvData)", FNameInfo);
				DR_Error(ErrorMsg);
				goto RETURN;
			}                                
        
			if (strstr(SwapDCC, "30/360") != NULL)
			{
				crInput->DCC = '3'; 
			}
			else if (strstr(SwapDCC, "ACT/365F") != NULL)
			{
				crInput->DCC = '5';
			}
			else if (strstr(SwapDCC, "ACT/360") != NULL)
			{
				crInput->DCC = '0';
			}
			else
			{
				sprintf (ErrorMsg, "cds dcc in file %s has to be '30/360' or"
						" 'ACT/365F' or 'ACT/360'! (crInputEnvData)", FNameInfo);
				DR_Error(ErrorMsg);
				goto RETURN;
			}

			/* read swap frequency */
			if (FindAndSkipComLine (1, stream, "cds frequency", "crInputEnvData", FNameInfo) == FAILURE)
			{        
				goto RETURN;
			}

			readerror = fscanf (stream, "%c \n", &(crInput->Freq));
			if (readerror != 1)
			{        
				sprintf (ErrorMsg, "Could not read cds frequency in %s! (crInputEnvData)", FNameInfo);
				DR_Error(ErrorMsg);
				goto RETURN;
			}
		}
    }
    
    /* OptionalSeq 0, additional recovery parameters for recovery correlation model */
    if (FindAndSkipComLineOptional(2, stream, "recovery dispersion", "crInputEnvData", 
                FNameInfo, '0') != FAILURE)
    {        
        readerror = fscanf (stream, "%lf %lf \n", 
                        &(crInput->RecoveryDispersion),
                        &(crInput->RecoveryBeta));
        if (readerror != 2)
        {        
            sprintf (ErrorMsg, "Could not read recovery additional parameters in %s! (crInputEnvData)", 
                FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    // optionalSeq 1
    if (FindAndSkipComLineOptional(2, stream, "IrCr CorrAdjFlag ", "crInputEnvData", FNameInfo, '1') != FAILURE)
    {
        readerror = fscanf (stream, "%d \n", &(crInput->CorrAdjFlag));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read IrCr CorrAdjFlag parameters in %s! (crInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    /**********************/
    /*  Read Vol diagonal */
    /**********************/

    /* currently, CalibFlag is hardcoded in the call 
            of this function in each product manager, where
            this flag can be Y iff DCC, freq and vols
            are present in the environment files */ 


    if (strchr(CalibFlag, 'Y') != NULL)
    {   
        if (VolDiagCR_Input_W (
                        &(crInput->BaseDate),
                        &(crInput->VolUnit),
                        &(crInput->NbVol),
                        crInput->VolDate,
                        crInput->SwapSt,
                        crInput->SwapMat,
                        crInput->Vol,
                        FNameVoldiag) == FAILURE) goto RETURN;

        /* Calib is off either CalibFlag from deal = 'N' or NbVol from env = 0 */
        if (crInput->NbVol < 0)
        {
            sprintf (ErrorMsg, "Nb Vol has to be >= 0! (crInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        else if (crInput->NbVol == 0)
        {
            crInput->CalibFlag = FALSE;
            crInput->SkipFlag = FALSE;
            crInput->VolUnit = -99999;
            crInput->NbVol = -99999;
            crInput->Freq = 'z';
            crInput->DCC = 'z';
        }
        else 
        {
            crInput->CalibFlag = TRUE;

            /* set Skip Flag */
            if (strchr(CalibFlag, '*') == NULL)
            {
                crInput->SkipFlag = FALSE;  /* vol points skipping not allowed */
            }
            else
            {
                crInput->SkipFlag = TRUE;   /* vol points skipping allowed      */
            }


        }
    }
    else
    {
        crInput->CalibFlag = FALSE;
        crInput->SkipFlag = FALSE;
        crInput->BaseDate = -99999;
        crInput->VolUnit = -99999;
        crInput->NbVol = -99999;
        crInput->Freq = 'z';
        crInput->DCC = 'z';
    }

    /*************************/
    /*  Read parameter file  */
    /*************************/

    {
        int     isFileNeeded;       /* TRUE = we need some data from text file */
        int     NbFactor = 1;
		
        if (strstr(NbFactorOWS, "1") == NULL) 
        {
            sprintf (ErrorMsg, "Overwrite for nb factor has to be '1'! (crInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (strstr (RhoOWS, "nil") == NULL)
        {
            sprintf (ErrorMsg, "Credit can only support 1 factor. No overwrite "
                "for rho is needed! (crInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
         
/*
        if (strstr(BackboneOWS, "0") == NULL)
        {
            sprintf (ErrorMsg, "Backbone OWS has to be '0' for credit (crInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
*/

        /* if not all overwrites are given, read CR param file */
        isFileNeeded = ((strstr(AlphaOWS,      "nil") != NULL) ||
                        (strstr(BetaOWS,       "nil") != NULL) ||
                        (strstr(SmileOWS,      "nil") != NULL));

        if (isFileNeeded)
        {

            if (Param_Input (&(crInput->NbFactor),
                             &(crInput->Alpha),
                             &(crInput->Beta),
                             NULL,
                             &(crInput->BackboneCoeff),
                             &(crInput->QLeft),
                             &(crInput->QRight),
                             &(crInput->FwdShift),
                             1,
                             FNameModlPar,
                             FALSE) == FAILURE) goto RETURN;
        }

        /* overwrite nb of factors */
        if (NbFactor != 1)
        {
            sprintf (ErrorMsg, "Number of factor must be 1! (crInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* overwrite alpha, beta and rho */
        if (strstr (AlphaOWS, "nil") == NULL) 
        {
            readerror = sscanf (AlphaOWS, 
                                "%lf \n", 
                                &(crInput->Alpha));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read overwrite string for spot vol! (crInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }

        if (strstr (BetaOWS, "nil") == NULL) 
        {
            readerror = sscanf (BetaOWS, 
                                "%lf \n", 
                                &(crInput->Beta));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (crInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }


        /* overwrite backbone coefficient */
        if (strstr (BackboneOWS, "nil") == NULL) 
        {
            readerror = sscanf (BackboneOWS, "%lf \n", &(crInput->BackboneCoeff));
            if (readerror != 1)
            {        
                DR_Error("crInputEnvData: cannot read backbone coefficient.");
                goto RETURN;
            }
        }

        /* overwrite smile parameters */
        if (strstr (SmileOWS, "nil") == NULL) 
        {
            readerror = sscanf (SmileOWS, 
                                "%lf \t%lf \t%lf \n", 
                                &(crInput->QLeft),
                                &(crInput->QRight),
                                &(crInput->FwdShift));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not read overwrite string for 2q! (crInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            crInput->QLeft = 1.0 - crInput->QLeft;
            crInput->QRight = 1.0 - crInput->QRight;
        }
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL) fclose(stream);

    if (status == FAILURE)
    {
        DR_Error("crInputEnvData: failed!");
    }
    return (status);

} /* crInputEnvData */


/*****  crInputCheck  ****************************************************/
/*
 *    Checks the integrity of the crInput structure.
 *    Must be called before Build_Sim
 */
int  crInputCheck (CR_INPUT  *crInput)
{
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];

    if (crInput == NULL) goto RETURN;

    /* Check zero curves */

    if (Term_Check_W(&(crInput->ZeroCurve)) == FAILURE) goto RETURN;

        
    /* Check model parameters */

    if (Param_Check(1L,
                    &(crInput->Alpha),
                    &(crInput->Beta),
                    NULL,
                    crInput->BackboneCoeff,
                    crInput->FwdShift) == FAILURE) goto RETURN;

    /* Backbone has to be zero */
    if (IS_NOT_ZERO(crInput->BackboneCoeff))
    {
        DR_Error("Only zero backbone is supported! (crInputCheck)");
        goto RETURN;
    }

    /* Check credit info */
    if (crInput->RecoveryRate < 0.0)
    {
        DR_Error("crInputCheck: credit recovery rate < 0.0!");
        goto RETURN;
    }

    if ((crInput->RecoveryRateConv != 'I') &&
        (crInput->RecoveryRateConv != 'E'))
    {
        DR_Error("crInputCheck: credit recovery rate convention must be 'I' or 'E'.");
        goto RETURN;
    }

    if( strlen( crInput->Currency ) != 3 )
    {
        DR_Error("crInputCheck: credit currency name invalid -- length != 3.");
        goto RETURN;
    }

    if( crInput->DefDate && Dateok( crInput->DefDate ) )
    {
        DR_Error("crInputCheck: credit default date invalid format");
        goto RETURN;
    }

    if (ABS(crInput->CorrTwk) > 1.0)
    {
        DR_Error("crInputCheck: default corr twk is > 1 or < -1!");
        goto RETURN;
    }

    /* check vol inputs */
    if(crInput->CalibFlag == TRUE)
    {
                 if (VolDiagCR_Check_W (crInput->BaseDate,        /* Volatility data */
										crInput->VolUnit,
										crInput->NbVol,
										crInput->VolDate,
										crInput->SwapSt,
										crInput->SwapMat,
										crInput->Vol) == FAILURE) goto RETURN;
				 
        /* NbVol = 0 is not allowed for CR */
        if (crInput->NbVol == 0)
        {
            sprintf (ErrorMsg, "Nb vol can't be zero! (crInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
     
        if ((crInput->Freq != 'A') &&
            (crInput->Freq != 'S') &&
            (crInput->Freq != 'Q') &&
            (crInput->Freq != 'M')) 
        {
            sprintf (ErrorMsg, "vol freq has to be 'A' or"
                    " 'S' or 'Q' or 'M'! (crInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
             
        if ((crInput->DCC != '0') &&
            (crInput->DCC != '3') &&
            (crInput->DCC != '5')) 
        {
            sprintf (ErrorMsg, "vol dcc has to be '0' or"
                    " '3' or '5'! (crInputCheck)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("crInputCheck: failed!");
    }
    return (status);

} /* crInputCheck */


/*********************************************************************************
 *    spInput
 *    Read the following files
 *      1) the basis zero curves
 *      2) SP info file
 *      3) SP volatility file
 *      4) SP parameter file
 * 
 *    The following inputs should be hardcoded when calling this function
 *      CalibFlag       = "N"
 *      NbFactorOWS     = "1"
 *      RhoOWS          = "nil"
 *      BackboneOWS     = "0" (to avoid reading file if all other OWS are given)
 *
 ********************************************************************************/
int  spInputEnvData (
                SP_INPUT  *spInput,       /* (O) sp Input Structure             */
                char      *FNameZC,       /* (I) file name for zero curve 0     */
                char      *FNameInfo,     /* (I) SP info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS)      /* (I) 2q string                      */
{
    int    status     = FAILURE;
    int    readerror;
    char   ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    /* basic check */
    if ((spInput      == NULL) ||
        (FNameZC      == NULL) ||
        (FNameInfo    == NULL) ||
        (FNameVoldiag == NULL) ||
        (FNameModlPar == NULL) ||
        (CalibFlag    == NULL) ||
        (NbFactorOWS  == NULL) ||
        (AlphaOWS     == NULL) ||
        (BetaOWS      == NULL) ||
        (RhoOWS       == NULL) ||
        (BackboneOWS  == NULL) ||
        (SmileOWS     == NULL)) goto RETURN;

    if (strcmp(CalibFlag, "Y" ) && 
        strcmp(CalibFlag, "Y*") &&
        strcmp(CalibFlag, "Y*S" ) &&
		strcmp(CalibFlag, "YS"  ) &&
        strcmp(CalibFlag, "N" ))
    {
        DR_Error ("Calibration flag must be 'Y' or 'Y*' or 'N'!(spInputEnvData)");
        goto RETURN;
    }
    
    /**********************/
    /*  Read Zero Curves  */
    /**********************/

    /* Read index curve */
    if (Term_Input_W (&(spInput->BasisZCurve), FNameZC) == FAILURE)
    {
        goto RETURN;
    }


    /********************/
    /*  Read SP info    */
    /********************/
    {
        int     basisType;
        int     tenor;
        char    SwapDCC[MAXBUFF];


        /* Open and read the SP info file */
        stream = fopen (FNameInfo, "r");

        if (stream == NULL)
        {
            sprintf (ErrorMsg, 
                     "Could not open file %s! (spInputEnvData)", 
                     FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Basis Type */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Basis type", 
                                "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%d \n", &(basisType));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                     "Could not read basis type in %s! (spInputEnvData)", 
                     FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if ( basisType != 1 &&
             basisType != 0) 
        {
            sprintf (ErrorMsg, 
                     "Basis type in %s has to be 0 or 1! (spInputEnvData)", 
                     FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        spInput->BasisType = (basisType==0)? 'S' : 'P' ;

        /* read basis rate tenor */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Basis rate tenor", "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, "%d \n", &(tenor));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                     "Could not read basis rate tenor in %s! (spInputEnvData)", 
                     FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Basis day count Convention */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "LIBOR day count convention", "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s \n", SwapDCC);
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read swap dcc in file %s! (spInputEnvData)", 
                    FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                
        
        if (strstr(SwapDCC, "30/360") != NULL)
        {
            spInput->BasisDCC = '3'; 
        }
        else if (strstr(SwapDCC, "ACT/365F") != NULL ||
                 strstr(SwapDCC, "ACT/365")  != NULL )
        {
            spInput->BasisDCC = '5';
        }
        else if (strstr(SwapDCC, "ACT/360") != NULL)
        {
            spInput->BasisDCC = '0';
        }
        else
        {
            sprintf (ErrorMsg, "Basis dcc in file %s has to be '30/360' or"
                    " 'ACT/365F' or 'ACT/360'! (spInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read Basis frequency */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Basis frequency", 
                                "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%c \n", &(spInput->BasisFreq));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read swap frequency in %s! (spInputEnvData)", 
                    FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read reference rate tenor */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Reference rate tenor", "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, "%d \n", &(tenor));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                     "Could not read basis rate tenor in %s! (spInputEnvData)", 
                     FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read reference rate day count Convention */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Reference rate day count convention", 
                                "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s \n", SwapDCC);
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read swap dcc in file %s! (spInputEnvData)", 
                    FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }                                
        
        if (strstr(SwapDCC, "30/360") != NULL)
        {
            spInput->LiborDCC = '3'; 
        }
        else if (strstr(SwapDCC, "ACT/365F") != NULL ||
                 strstr(SwapDCC, "ACT/365")  != NULL)
        {
            spInput->LiborDCC = '5';
        }
        else if (strstr(SwapDCC, "ACT/360") != NULL)
        {
            spInput->LiborDCC = '0';
        }
        else
        {
            sprintf (ErrorMsg, "LIBOR dcc in file %s has to be '30/360' or"
                    " 'ACT/365F' or 'ACT/360'! (spInputEnvData)", FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read reference rate frequency */
        if (FindAndSkipComLine (1, 
                                stream, 
                                "Reference rate frequency", 
                                "spInputEnvData", 
                                FNameInfo) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%c \n", &(spInput->LiborFreq));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, 
                    "Could not read swap frequency in %s! (spInputEnvData)", 
                    FNameInfo);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }
        
    /**********************/
    /*  Read Vol diagonal */
    /**********************/

    /* currently, CalibFlag is hardcoded in the call 
            of this function in each product manager, where
            this flag can be Y iff DCC, freq and vols
            are present in the environment files */ 


    if (strchr(CalibFlag, 'Y') != NULL)
    {           
        if (VolDiag_Input_W (
                        &(spInput->BaseDate),
                        &(spInput->VolUnit),
                        &(spInput->NbVol),
                        spInput->VolDate,
                        spInput->SwapSt,
                        spInput->SwapMat,
                        spInput->Vol,
                        spInput->VolType,
                        FNameVoldiag) == FAILURE) goto RETURN;

        /* Calib is off either CalibFlag from deal = 'N' or NbVol from env = 0 */
        if (spInput->NbVol < 0)
        {
            sprintf (ErrorMsg, 
                    "Nb spread Vol has to be >= 0! (spInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        else if (spInput->NbVol == 0)
        {
            spInput->VolUnit = -99999;
            spInput->NbVol = -99999;
        }
    }
    else
    {
        spInput->BaseDate = -99999;
        spInput->VolUnit = -99999;
        spInput->NbVol = -99999;
    }

    /*************************/
    /*  Read parameter file  */
    /*************************/

    {
        int     isFileNeeded;       /* TRUE = we need some data from text file */
        int     NbFactor = 1;
		
        if (strstr(NbFactorOWS, "1") == NULL) 
        {
            sprintf (ErrorMsg, 
                     "Overwrite for nb factor has to be '1'! (spInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (strstr (RhoOWS, "nil") == NULL)
        {
            sprintf (ErrorMsg, 
                    "Basis spread can only support 1 factor. No overwrite "
                    "for rho is needed! (spInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
         

        /* if not all overwrites are given, read CR param file */
        isFileNeeded = ((strstr(AlphaOWS,      "nil") != NULL) ||
                        (strstr(BetaOWS,       "nil") != NULL) ||
                        (strstr(SmileOWS,      "nil") != NULL));

        if (isFileNeeded)
        {

            if (Param_Input (&(spInput->NbFactor),
                             &(spInput->Alpha),
                             &(spInput->Beta),
                             NULL,
                             &(spInput->BackboneCoeff),
                             &(spInput->QLeft),
                             &(spInput->QRight),
                             &(spInput->FwdShift),
                             1,
                             FNameModlPar,
                             FALSE) == FAILURE) goto RETURN;
        }

        /* overwrite nb of factors */
        if (NbFactor != 1)
        {
            sprintf (ErrorMsg, "Number of factor must be 1! (spInputEnvData)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* overwrite alpha, beta and rho */
        if (strstr (AlphaOWS, "nil") == NULL) 
        {
            readerror = sscanf (AlphaOWS, 
                                "%lf \n", 
                                &(spInput->Alpha));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, 
                        "Could not read overwrite string for spot vol! "
                        "(crInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }

        if (strstr (BetaOWS, "nil") == NULL) 
        {
            readerror = sscanf (BetaOWS, 
                                "%lf \n", 
                                &(spInput->Beta));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, 
                        "Could not read overwrite string for mean reversion! "
                        "(spInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }


        /* overwrite backbone coefficient */
        if (strstr (BackboneOWS, "nil") == NULL) 
        {
            readerror = sscanf (BackboneOWS, 
                                "%lf \n", 
                                &(spInput->BackboneCoeff));
            if (readerror != 1)
            {        
                DR_Error("crInputEnvData: cannot read backbone coefficient.");
                goto RETURN;
            }
        }

        /* overwrite smile parameters */
        if (strstr (SmileOWS, "nil") == NULL) 
        {
            readerror = sscanf (SmileOWS, 
                                "%lf \t%lf \t%lf \n", 
                                &(spInput->QLeft),
                                &(spInput->QRight),
                                &(spInput->FwdShift));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not read overwrite string for 2q! (crInputEnvData)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            spInput->QLeft = 1.0 - spInput->QLeft;
            spInput->QRight = 1.0 - spInput->QRight;
        }
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL) fclose(stream);

    if (status == FAILURE)
    {
        DR_Error("spInputEnvData: failed!");
    }
    return (status);

} /* spInputEnvData */


#ifndef OPEN_FILE
#define OPEN_FILE(filename, routine)                                \
{                                                                   \
    fp = fopen((filename),"w");                                     \
    if (fp == NULL)                                                 \
    {                                                               \
        DR_Error("%s: Unable to open file %s", (routine), (filename));\
        goto RETURN;                                                  \
    }                                                               \
}                                                                   
#endif



/*****  DR_ATOF  ****************************************************/
/*
*      Convert a string into double if the string starts with a number
*      otherwise, return error
*/
static int DR_ATOF(char    *InputString,    /* (I) */
                   double  *OutputDouble)   /* (O) */
{
    int     status = FAILURE;
    long    readerror;

    readerror = sscanf (InputString, "%lf \n", OutputDouble);
    if (readerror != 1) goto RETURN;

    status = SUCCESS;

RETURN:
    return(status);
}


/*****  corrInputEnvData  ****************************************************/
/*
 *    Read corr matrix from correlation.dat.
 *    Memory should have been allocated already.
 *    Replaced with OverWriteString when given
 *    Check inputs
 *
 */
int corrInputEnvData(double  **corr,               /* (O) corr matrix */
                     long     NbAssets,            /* (I) nb of assets    */
                     char     *FName,              /* (I) Corr file name  */
                     char   ***CorrOWS)            /* (I) Corr overwrite string
*/
{
    int       status = FAILURE;
    int       readerror;
    char      ErrorMsg[MAXBUFF];

    FILE     *stream = NULL;

    long    i, j;

    long    NbAssetsL = 0L;
    int    isFileNeeded = FALSE;


    /* basic check */
    if ((corr           == NULL) ||
        (FName          == NULL) ||
        (CorrOWS        == NULL)) goto RETURN;

    if (NbAssets < 1L) goto RETURN;

    /* no need to read file if all overwrites are given */
    for (i=0; i<NbAssets; i++)
    {
        for (j=0; j<NbAssets; j++)
        {
            if (strstr(CorrOWS[i][j], "nil") != NULL)
            {
                isFileNeeded = TRUE;
                break;
            }
        }
        if (isFileNeeded) break;
    }

    if (isFileNeeded)
    {
        stream = fopen (FName, "r");
        if (stream == NULL)
        {
            sprintf (ErrorMsg, "corrInputEnvData: Cannot open file %s!",
                     FName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* read NbAssets and check */
        if (FindAndSkipComLine (1, stream,
                               "Nb assets",
                               "corrInputEnvData",
                               FName) == FAILURE) goto RETURN;

        readerror = fscanf (stream, "%ld \n", &(NbAssetsL));
        if (readerror != 1)
        {
            DR_Error("corrInputEnvData: Cannot read Nb assets");
            goto RETURN;
        }

        if (NbAssetsL != NbAssets)
        {
            DR_Error("corrInputEnvData: Nb assets inconsistent!");
            goto RETURN;
        }

        /* read corr matrix */
        if (FindAndSkipComLine (1, stream,
                               "Corr Matrix",
                               "corrInputEnvData",
                               FName) == FAILURE) goto RETURN;

        for (i=0; i<NbAssets; i++)
        {
            for (j=0; j<NbAssets; j++)
            {
                readerror = fscanf(stream, "%lf ", &(corr[i][j]));
                if (readerror != 1)
                {
                    DR_Error("corrInputEnvData: Cannot read corr data.");
                    goto RETURN;
                }
            }
        }
    }

    /* use overwrite strings if given */
    for (i=0; i<NbAssets; i++)
    {
        for (j=0; j<NbAssets; j++)
        {
            if (strstr(CorrOWS[i][j], "nil") == NULL)
            {
                if (DR_ATOF(CorrOWS[i][j], &(corr[i][j]))
                    == FAILURE)
                {
                    DR_Error("corrInputEnvData: Cannot read corr overwrite.");
                    goto RETURN;
                }
            }
        }
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL) fclose(stream);

    if (status == FAILURE)
    {
        DR_Error("corrInputEnvData: failed!");
    }

    return (status);

} /* corrInputEnvData */



/*****  Alloc_CorrOWS *************************************************
 *
 *      allocate memory for correlation overwrite matrix
 */

int Alloc_CorrOWS   (long       NbAssets,   /* (I) */
                     char   ****CorrOWS)    /* (O) */
{
    int status = FAILURE;
    int i, j;

    char ***CorrOWSL = NULL;;

    CorrOWSL = calloc(NbAssets, sizeof(char **));
    if (CorrOWSL == NULL) goto RETURN;

    for(i=0; i<NbAssets; i++)
    {
        CorrOWSL[i] = calloc(NbAssets, sizeof(char *));
        if (CorrOWSL[i] == NULL) goto RETURN;
    }

    for (i=0; i<NbAssets; i++)
        for (j=0; j<NbAssets; j++)
        {
            CorrOWSL[i][j] = calloc(MAXBUFF, sizeof(char));
            if (CorrOWSL[i][j] == NULL) goto RETURN;
        }

    *CorrOWS = CorrOWSL;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        Free_CorrOWS(NbAssets, CorrOWSL);
    }
    return (status);

}/* Alloc_CorrOWS */



/*****  Free_CorrOWS *************************************************
 *
 *      free memory for correlation overwrite matrix
 */

int Free_CorrOWS   (long   NbAssets,       /* (I) */
                    char   ***CorrOWS)
{
    int status = FAILURE;
    int i, j;

    if (CorrOWS == NULL) return(SUCCESS);
    if (NbAssets <= 0L) goto RETURN;

    for (i=0; i<NbAssets; i++)
    {
        if (CorrOWS[i] != NULL)
        {
            for (j=0; j<NbAssets; j++)
            {
                if (CorrOWS[i][j] != NULL) free(CorrOWS[i][j]);
            }
            free(CorrOWS[i]);
        }
    }
    free(CorrOWS);

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Free_CorrOWS: failed.");
    }
    return(status);

}/* Free_CorrOWS */





/*****************************  PrintTrancheEnv  ****************************/
/*                                                                          */
/* Print market information in files in the directory specified.            */
/* This function dumps the input files for tranches and cdo2 trades,        */
/* which is useful to debug Kapital streams using the wrapper.              */
/* This function should not be used for other trades as dummy hardcoded     */
/* inputs are printed when we don't need them for tranches and cdo2 trades. */
/*                                                                          */
/****************************************************************************/
int PrintTrancheEnv(const IR_INPUT   *irInp,  /* (I) IR_INPUT structure */
					long	      	  nbCR,   /* (I) number of credit names */
					const CR_INPUT   *crInp,  /* (I) array[nbCR] of CR_INPUT struct  */
                    const char       *dir)    /* (I) dir where the files are printed */
{
    static char routine[] = "PrintTrancheEnv";
    int   status = FAILURE;
    long  i,j;
    char
        CR_FILES_HEAD[4][MAXBUFF] = {
                    "cr_zcurve_", 
                    "cr_info_",    
                    "cr_voldiag_",
                    "cr_mpar_"},
        IR_FILES_HEAD[MAXBUFF] = "ir_curve",
        FILEBUF[MAXBUFF];

                    
    FILE *fp = NULL;

    if (dir==NULL) goto RETURN;

    /* ir_info_0 */
    sprintf(FILEBUF, "%sir_info_0.dat",dir);
    OPEN_FILE(FILEBUF, routine);

    fprintf(fp,"# Money Market basis (360 or 365)\n");
    fprintf(fp,"360\n");
    fprintf(fp,"# Swap Day Count Convention \n");
    fprintf(fp,"ACT/360\n");
    fprintf(fp,"# Swap Frequency\n");
    fprintf(fp,"A\n");
    DR_FREE(fp,fclose);

    /* ir_mpar_0 */
    sprintf(FILEBUF, "%sir_mpar_0.dat",dir);
    OPEN_FILE(FILEBUF, routine);

    fprintf(fp,"# Number of Factors\n");
    fprintf(fp,"1\n");
    fprintf(fp,"# Factor mean reversion coefficients\n");
    fprintf(fp,"0.0050\n");
    fprintf(fp,"# Factor weights\n");
    fprintf(fp,"0.0964"); 
    fprintf(fp,"# Factor correlations\n");
    fprintf(fp,"\n"); 
    fprintf(fp,"# Backbone coefficient\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# Smile parameters\n");
    fprintf(fp,"0 0 0 \n");
    DR_FREE(fp,fclose);

    /* ir_voldiag_0 */
    sprintf(FILEBUF, "%sir_voldiag_0.dat",dir);
    OPEN_FILE(FILEBUF, routine);

    fprintf(fp,"# Start date\n");
    fprintf(fp,"%ld\n",irInp->ZeroCurve[0].ValueDate);
    fprintf(fp,"# Units\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# No of entries\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# Option exp date, swap start date, swap mat date, volatility (percentage), Cap/Swaption\n");
    DR_FREE(fp,fclose);

    /* ir_curvej_0 */
    for (j=0;j<3;++j)
    {
        sprintf(FILEBUF, "%s%s%ld_0.dat",dir, IR_FILES_HEAD, j);
        OPEN_FILE(FILEBUF, routine);

        fprintf(fp,"# Start date\n");
        fprintf(fp,"%ld\n",irInp->ZeroCurve[j].ValueDate);
        fprintf(fp,"# Day Count Convention\n");
        fprintf(fp,"ACT/360\n");
        fprintf(fp,"# Compounding Frequency\n");
        fprintf(fp,"A\n");
        fprintf(fp,"# No of entries\n");
        fprintf(fp,"%d\n", irInp->ZeroCurve[j].NbZero);
        fprintf(fp,"#zero maturity (yyyymmdd)  and rates (annual compounding)\n");
        for (i=0; i<irInp->ZeroCurve[j].NbZero;++i)
        {
            fprintf(fp,"%ld\t%lf\n", irInp->ZeroCurve[j].ZeroDate[i], irInp->ZeroCurve[j].Zero[i]*100);
        }
        DR_FREE(fp,fclose);

    }

    for (j=0;j<nbCR; ++j)
    {
        /* cr_zcurve_j */
        sprintf(FILEBUF, "%s%s%ld.dat", dir, CR_FILES_HEAD[0], j);
        OPEN_FILE(FILEBUF, routine);

        fprintf(fp,"# Start date\n");
        fprintf(fp,"%ld\n",crInp[j].ZeroCurve.ValueDate);
        fprintf(fp,"# Day Count Convension\n");
        fprintf(fp,"ACT/365F\n");
        fprintf(fp,"# Compounding Frequency\n");
        fprintf(fp,"A\n");
        fprintf(fp,"# No of entries\n");
        fprintf(fp,"%d\n",crInp[j].ZeroCurve.NbZero);
        fprintf(fp,"#zero maturity (yyyymmdd)  and rates (annual compounding)\n");
        for (i=0;i<crInp[j].ZeroCurve.NbZero;++i)
        {
            fprintf(fp,"%ld\t%lf\n", crInp[j].ZeroCurve.ZeroDate[i], crInp[j].ZeroCurve.Zero[i]*100);
        }
        DR_FREE(fp,fclose);

        
        /* cr_info_j */
        sprintf(FILEBUF, "%s%s%ld.dat", dir, CR_FILES_HEAD[1], j);
        OPEN_FILE(FILEBUF, routine);

        fprintf(fp,"# Recovery rate (e.g. 0.1 is 10%%)\n");
        fprintf(fp,"%lf\n", crInp[j].RecoveryRate);
        fprintf(fp,"# Recovery rate convention (I = inc accrued interest E = exclude accrued interest)\n");
        fprintf(fp,"I\n");
        fprintf(fp,"# SPN\n");
        fprintf(fp,"%ld\n",crInp[j].SPN);
        fprintf(fp,"# Counterparty name\n");
        fprintf(fp,"%s\n",crInp[j].CtpyName);
        fprintf(fp,"# Currency\n");
        fprintf(fp,"EUR\n");
        fprintf(fp,"# Seniority\n");
        fprintf(fp,"Subordinated LT2\n");
        fprintf(fp,"# Curve ID\n");
        fprintf(fp,"%ld\n",crInp[j].CurveId);
        fprintf(fp,"# Defaulted\n");
        fprintf(fp,"%ld\n",crInp[j].DefDate);
        fprintf(fp,"# Default Correlation Amount\n");
        fprintf(fp,"0\n");
        fprintf(fp,"# CDS Day Count Convention\n");
        fprintf(fp,"ACT/360\n");
        fprintf(fp,"# Swap Frequency\n");
        fprintf(fp,"Q\n");
        fprintf(fp,"# Recovery dispersion and rec beta\n");
        fprintf(fp,"%lf\t%lf\n", crInp[j].RecoveryDispersion, crInp[j].RecoveryBeta);
        DR_FREE(fp,fclose);



        /* cr_voldiag_j */
        sprintf(FILEBUF, "%s%s%ld.dat", dir, CR_FILES_HEAD[2], j);
        OPEN_FILE(FILEBUF, routine);

        fprintf(fp,"# Start date\n");
        fprintf(fp,"%ld\n",crInp[j].ZeroCurve.ValueDate);
        fprintf(fp,"# Units\n");
        fprintf(fp,"0\n");
        fprintf(fp,"# No of entries\n");
        fprintf(fp,"0\n");
        fprintf(fp,"# Option exp date, swap start date, swap mat date, volatility (percentage)\n");
        DR_FREE(fp,fclose);

        /* cr_mpar_j */
        sprintf(FILEBUF, "%s%s%ld.dat", dir, CR_FILES_HEAD[3], j);
        OPEN_FILE(FILEBUF, routine);

        fprintf(fp,"# Number of Factors\n");
        fprintf(fp,"1\n");
        fprintf(fp,"# Factor mean reversion coefficients\n");
        fprintf(fp,".0484 \n");
        fprintf(fp,"# Factor weights\n");
        fprintf(fp,".1414 \n");
        fprintf(fp,"# Factor correlations\n");
        fprintf(fp,"# Backbone coefficient\n");
        fprintf(fp,"0\n");
        fprintf(fp,"# Smile parameters\n");
        fprintf(fp,"0 0 0 \n");
        DR_FREE(fp,fclose);
    }


    sprintf(FILEBUF, "%ssummary.dat",dir);
    OPEN_FILE(FILEBUF, routine);

    fprintf(fp,"### Environment Section\n");
    fprintf(fp,"# today's date\n");
    fprintf(fp,"%ld\n",irInp->ZeroCurve[0].ValueDate);
    fprintf(fp,"### End Section\n");
    fprintf(fp,"### IR Section\n");
    fprintf(fp,"# Nb IR\n");
    fprintf(fp,"1\n");
    fprintf(fp,"### End Section\n");
    fprintf(fp,"### Basis Section\n");
    fprintf(fp,"# Nb Basis\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# Currency denomination [0..NbIr-1], reference curve [0, 1, 2 ]\n");
    fprintf(fp,"### End Section\n");
    fprintf(fp,"### FX Section\n");
    fprintf(fp,"# Nb FX\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# Foreign IR identifier (Domestic IR id is always 0)\n");
    fprintf(fp,"### End Section\n");
    fprintf(fp,"### Equity Section\n");
    fprintf(fp,"# Nb Equity\n");
    fprintf(fp,"0\n");
    fprintf(fp,"# Currency denomination\n");
    fprintf(fp,"### End Section\n");
    fprintf(fp,"### Credit Section\n");
    fprintf(fp,"# Nb Credit\n");
    fprintf(fp,"%ld\n",nbCR);
    fprintf(fp,"# Currency denomination [0..NbIr-1]\n");
    for (i=0; i<nbCR; ++i)
    {
        fprintf(fp,"0\n");
    }
    fprintf(fp,"### End Section\n");
    DR_FREE(fp,fclose);

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        DR_FREE(fp, fclose);
    }    

    return (status);

}  /* PrintTrancheEnv */







/**--------------------------------------------------------------
 * Reads the market data from a DR Wrapper type 2
 * and creates a CrxTDrWrapperData.
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * Convenience routine that calls 
 * "CrxTDrWrapperDataGetFull" with CRX_DRW_TYPE2_2CURVES.
 */


static	char	todayFnam[] = "today.dat";
static	char	discZcFnam[] = "disczero.dat";
static	char	riskZcFnam[] = "riskzero.dat";
static	char	zcFnam[] = "zero.dat";
static	char	swvolFnam[] = "swapvol.dat";
static	char	mparamFnam[] = "modelParameters.dat";




/*f--------------------------------------------------------------
 * Reads the market data from a DR wrapper
 * and creates a "CrxTDrWrapperData".
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * The "options" argument dteremines the type of the wrapper
 * and the version. Possible values are:\\
 * CRX_DRW_TYPE2_2CURVES: Reads a type 2 wrapper
 * with two zero curves (disczero.dat and zero.dat)
 * CRX_DRW_TYPE2_3CURVES: Reads a type 2 wrapper
 * with three curves (disczero.dat, zero.dat and riskzero.dat).
 * DRI_DRW_BASIS: Reads a basis type wrapper
 * with four curves (disczero.dat, zero.dat, riskzero.dat, basiszero.dat),
 * plus basis spot vol curve.
 * Return NULL TCurve pointer if file not available.
 */

int CrxTDrWrapperDataGetFull(
	char *pathdir,			/* (I) tmp direcory name (or NULL) */
	long options,			/* (I) see description */
	CrxTDrWrapperData **drWrapper)	/* (O) dr wrapper data */
{
static	char	routine[] = "CrxTDrWrapperDataGet";
	int	status = FAILURE;
	char	fnam[256];

    FILE   *fp = NULL;
    char    buf1[256], *buf = buf1;

	/* */
	if ((*drWrapper = NEW(CrxTDrWrapperData)) == NULL) {
		goto RETURN;
	}

	(*drWrapper)->fDiscZcCurve = NULL;
	(*drWrapper)->fZcCurve = NULL;
	(*drWrapper)->fRiskZcCurve = NULL;
	(*drWrapper)->fBasisZcCurve = NULL;
	(*drWrapper)->fBvCurve = NULL;
	(*drWrapper)->fBSVolCurve = NULL;
	(*drWrapper)->fCmsSwMat = NULL;

	/* Swaption vol Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, swvolFnam);
	} else {
		strcpy(fnam, swvolFnam);
	}

#ifdef _SKIP

	if (CrxTSwaptionMatrix2DFileRead(&(*drWrapper)->fCmsSwMat, fnam,
		TSWAPTION_MATRIX_FMT_LON) != SUCCESS) goto RETURN;
#endif

	/* Disc Zero Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, discZcFnam);
	} else {
		strcpy(fnam, discZcFnam);
	}
	if (CrxTCurveFileRead(&((*drWrapper)->fDiscZcCurve), fnam) != SUCCESS)
		goto RETURN;


	/* Zero Curve (index) */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, zcFnam);
	} else {
		strcpy(fnam, zcFnam);
	}
	if (CrxTCurveFileRead(&((*drWrapper)->fZcCurve), fnam) != SUCCESS)
		goto RETURN;

	/* Zero Curve (risky) */
	if (options >= CRX_DRW_TYPE2_3CURVES) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, riskZcFnam);
	    } else {
		strcpy(fnam, riskZcFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fRiskZcCurve = NULL;
	    else if (CrxTCurveFileRead(&((*drWrapper)->fRiskZcCurve), fnam) 
                    != SUCCESS)
		goto RETURN;
	} else {
	    (*drWrapper)->fRiskZcCurve = NULL;
	}


    /* Read additional parameters from zero.dat - day count conv's
	 * and swap frequency */
    if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    DR_Error("%s: can't open `%s'.\n",
		      routine, fnam);
	    goto RETURN;
	}


#ifdef _SKIP
    READ_DATA(DRL_STRING_T, &buf,	        "start date");
    READ_DATA(DRL_LONG_T, &((*drWrapper)->fMMDenom),  "mm basis");

    READ_DATA(DRL_STRING_T, &buf,	        "swap frequency");
    if (CrxStrLongValueScan(buf, "swap frequency", 
			    &((*drWrapper)->fCmsSwMat->swapPayFreq),
		"A", 1L,
		"S", 2L,
		NULL) != SUCCESS) goto RETURN;

    READ_DATA(DRL_STRING_T, &buf,	        "swap basis");
    if (CrxStrLongValueScan(buf, "swap basis", &((*drWrapper)->fSwDcc),
		"ACT", (TDayCount) GTO_B30_360,
		"360", (TDayCount) GTO_ACT_360,
		"365", (TDayCount) GTO_ACT_365F,
		NULL) != SUCCESS) goto RETURN;

        fclose(fp); fp = NULL;
        

	/* Base vol Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, bvFnam);
	} else {
		strcpy(fnam, bvFnam);
	}
	if (CrxTCurveLondonBaseVolFileRead(&((*drWrapper)->fBvCurve), fnam)
			!= SUCCESS) goto RETURN;

	/* Basis vol Curve */
	if (options == DRI_DRW_BASIS) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, bsvolFnam);
	    } else {
		strcpy(fnam, bsvolFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fBSVolCurve = NULL;
	    else if (CrxTCurveLondonBaseVolFileRead(
			&((*drWrapper)->fBSVolCurve), 
			fnam) != SUCCESS) 
		goto RETURN;
	} else {
	    (*drWrapper)->fBSVolCurve = NULL;
	}

#endif

	/*
	 * Set today date to value date (Kapital)
	 */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, todayFnam);
	} else {
		strcpy(fnam, todayFnam);
	}
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    DR_Error("%s: can't open `%s'.\n",
		      routine, fnam);
	    goto RETURN;
	}
	READ_DATA(DRL_STRING_T, &buf,	        "today");
	if (CrxTDateScanYMD(buf, &((*drWrapper)->fToday)) != SUCCESS) {
		goto RETURN;
	}
        fclose(fp); fp = NULL;


	/*
	 * Read model parameters
	 */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, mparamFnam);
	} else {
		strcpy(fnam, mparamFnam);
	}
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    DR_Error("%s: can't open `%s'.\n",
		      routine, fnam);
	    goto RETURN;
	}

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f1Beta),
		"oneFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f1Weight),
		"oneFactorVolatility1");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f1Ppy),
		"oneFactorPPY");

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Beta1),
		"twoFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Beta2),
		"twoFactorMeanReversion2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Weight1),
		"twoFactorVolatility1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Weight2),
		"twoFactorVolatility2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Corr12),
		"twoFactorCorrelation1and2");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f2Ppy),
		"twoFactorPPY");

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta1),
		"threeFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta2),
		"threeFactorMeanReversion2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta3),
		"threeFactorMeanReversion3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight1),
		"threeFactorVolatility1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight2),
		"threeFactorVolatility2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight3),
		"threeFactorVolatility3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr12),
		"threeFactorCorrelation1and2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr13),
		"threeFactorCorrelation1and3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr23),
		"threeFactorCorrelation2and3");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f3Ppy),
		"threeFactorPPY");


        fclose(fp); fp = NULL;


	/* OK */
	status = SUCCESS;
RETURN:

        if (fp != NULL) fclose(fp);
	if (status != SUCCESS)
		DR_Error("%s: failed.\n", routine);
	return(status);
}



/*f--------------------------------------------------------------
 * Frees a DRW allocated by "CrxTDrWrapperDataGet".
 */

int CrxTDrWrapperDataFree(CrxTDrWrapperData *drWrapper)
{
	if (drWrapper != NULL) {
		GtoFreeTCurve(drWrapper->fDiscZcCurve);
		GtoFreeTCurve(drWrapper->fZcCurve);
		GtoFreeTCurve(drWrapper->fRiskZcCurve);
		GtoFreeTCurve(drWrapper->fBasisZcCurve);
#ifdef _SKIP
		GtoFreeTCurve(drWrapper->fBvCurve);
		GtoFreeTCurve(drWrapper->fBSVolCurve);
		CrxTSwaptionMatrix2DFree(drWrapper->fCmsSwMat);
#endif

		FREE(drWrapper);
	}
	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * Creates the output price file for a DRW type II.
 */

int CrxTDrWrapperDataPutPrice(double value)
{
static  char    routine[] = "CrxTDrWrapperDataPutPrice";
        FILE    *fp = NULL;
static  char    fnam[] = "price";

        if ((fp = fopen(fnam, "w")) == NULL) {
                DR_Error("%s: can't open `%s'.\n",
                         routine, fnam);
                return(FAILURE);
        }
        fprintf(fp, "%.8f", value);
        fclose(fp);
        return(SUCCESS);
}





int CrxTEqStatGet(
    char 	     *pathdir,		/* (I) tmp direcory name (or NULL) */
    char 	     *eqStaFnam,	/* (I) file name (equity.sta)  */
    long	     busDayConv,    /* (I) Business Day Conv    */
    char	     *holidayFile,  /* (I) Holiday File       */
    TEqStatData  **eqStatData)  /* (O) equity static data */
{
    static	char	routine[] = "CrxTEqStatGet";
    int	status = FAILURE;

    char   fnam[256];
    FILE   *fp = NULL;

    long   divPaySettleDays = 0;  

    char   buf[256];
    char   *label = "Dividend";

    long   numDivs;
    TDate  *exDivDates = NULL;
    TDate  *payDivDates = NULL;
    double *divAmount = NULL;
    long   *divTypes = NULL;
    char   *divStrs = NULL;

    long     	numSettle, settleDays;
    char 	    settleType;
    TDateList   *lastTradingDates = NULL;
    TDateList   *settlementDates = NULL;
  
    TDate       stockDate = 0;
    TDateList   *stmPeriodDates = NULL;
    long        numStmPeriods = 1;
    long	    stmPeriod[1];

    long i;
    
    TDividendList     *divList = NULL;
    TEquitySettlement *stm = NULL;


    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	    sprintf(fnam, "%s/%s", pathdir, eqStaFnam);
    } else {
	    strcpy(fnam, eqStaFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	    DR_Error("%s: can't open `%s'.\n",
	         	  routine, fnam);
	    goto RETURN;
    }
    

    /* 
     * Skip all before divident line start
     */
    if(CrxFAdvanceToToken(fp,
			  label) == FAILURE)
    {
	    DR_Error("%s: Cannot find string token %s.\n", routine, label);
	    goto RETURN;
    }

    /* 
     * The file pointer stops after the token string.  
     * Need to skip rest of the line 
     */
    if(fgets(buf, sizeof(buf), fp) == NULL) 
    {
	    DR_Error("%s: EOF encountered before divident list starts.\n",
                     routine);
	    goto RETURN;
    }

    /* 
     * Read the divident list
     */
    READ_DATA(DRL_LONG_T, &numDivs,    "number of dividents");

    if (numDivs <= 0)     /* No dividend. Create a dummy one */ 
    {
	numDivs = 1;
	if((exDivDates = NEW_ARRAY(TDate, numDivs))  == NULL ||
	   (divAmount  = NEW_ARRAY(double, numDivs)) == NULL ||
	   (divStrs    = NEW_ARRAY(char, numDivs))   == NULL)
		goto RETURN;

	exDivDates[0] = 100;
	divAmount[0]  = 0.0;
	divStrs[0]    = 'D';
    }
    else
    {
    	if(CrxLilVectArrayFpReadV(fp, 
			      	  numDivs,
			      	  DRL_TDATE_T,   (void*) &exDivDates,
			      	  DRL_DOUBLE_T,  (void*) &divAmount,
			      	  DRL_CHAR_T,    (void*) &divStrs,
			      	  DRL_NULL_T) == FAILURE)
    	{  
        	DR_Error("%s: Cannot read divident lists.\n", routine);
        	goto RETURN;
    	}
    }

    /*
     * Read settlement date
     */
    READ_DATA(DRL_LONG_T, &numSettle,    "number of settlement points");

    if (numSettle <= 0)    /* Rolling settlement */
    {
	    settleType = 'R';
	    READ_DATA(DRL_LONG_T, &settleDays,    "settlement days");

    	stmPeriodDates = GtoNewEmptyDateList(0);
    	stmPeriod[0] = settleDays;

	    divPaySettleDays = settleDays;
    }
    else	 /* Fixed settlement */
    {
	settleType = 'F';

	if ((lastTradingDates = GtoNewEmptyDateList(0)) == NULL ||
	    (settlementDates  = GtoNewEmptyDateList(0)) == NULL)
		goto RETURN;

	lastTradingDates->fNumItems = numSettle;
	settlementDates->fNumItems = numSettle;

	if(CrxLilVectArrayFpReadV(fp,
                               numSettle,
                               DRL_TDATE_T, (void*) &(lastTradingDates->fArray),
                               DRL_TDATE_T, (void*) &(settlementDates->fArray),
                               DRL_NULL_T) == FAILURE)
	{
	    DR_Error("%s: Cannot read fixed settlement lists.\n", routine);
            goto RETURN;
	}
    }


    /* 
     * create divident pay dates from ex-divident dates and settlement dates
     */
    if((payDivDates = NEW_ARRAY(TDate, numDivs)) == NULL ||
       (divTypes    = NEW_ARRAY(long, numDivs)) == NULL)
    {
	DR_Error("%s: Error allocating memory for payDivDates or divTypes.\n",
		  routine);
	goto RETURN;
    }

    /* 
     * Generate pay dates and dividend types 
     * For continuous dividend, it requires ex-div dates 
     * and payment dates are equal
     */
    for (i=0; i<=numDivs-1; i++){

	switch (toupper(divStrs[i]))
	{
	case 'D':
		divTypes[i] = GtoDIVIDENDRATE_AMOUNT;
		if(GtoDateFromBusDaysOffset(exDivDates[i],
                                    	    divPaySettleDays,
                                    	    holidayFile,
                                    	    &payDivDates[i]) == FAILURE)
	    		goto RETURN;
		break;
	case 'Y':
		divTypes[i] = GtoDIVIDENDRATE_PERCENT;
		if(GtoDateFromBusDaysOffset(exDivDates[i],
                                    	    divPaySettleDays,
                                    	    holidayFile,
                                    	    &payDivDates[i]) == FAILURE)
	    		goto RETURN;
		break;
	case 'C':
		divTypes[i] = GtoDIVIDENDRATE_CONTINUOUS;
		payDivDates[i] = exDivDates[i];
		break;
	default:
		DR_Error("%s: Unknown divident type %c for the %dth divident."
			  "\n",
                  	  routine, divStrs[i], i+1);
        	goto RETURN;
	}
    }


    /* Construct TDividendList and TEquitySettment
     */
    if((divList = GtoNewTDividendList(
                      numDivs,
			          exDivDates,
			          payDivDates,
				      divTypes,
				      divAmount)) == NULL ||
       (stm = GtoEqSettleNew( stockDate,
			      stmPeriodDates,
    			  numStmPeriods,
    			  stmPeriod, 
			      lastTradingDates,
			      settlementDates,
			      holidayFile,
			      busDayConv,
			      holidayFile,
			      busDayConv,
			      GTO_EQ_STM_DIRN_INC)) == NULL)	
	goto RETURN;
			 

    /* 
     * Construct eqStatData
     */

    if((*eqStatData = CrxNewTEqStatData(divList,
				      settleType,
				      stm)) == NULL) 
	goto RETURN;


    status = SUCCESS;

  RETURN:

    FREE(exDivDates);
    FREE(payDivDates);
    FREE(divAmount);
    FREE(divStrs);
    FREE(divTypes);
    GtoFreeDateList(stmPeriodDates);
    GtoFreeDateList(lastTradingDates);
    GtoFreeDateList(settlementDates);
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
	DR_Error("%s: failed.\n", routine);

    return(status);
}


TEqStatData* CrxNewTEqStatData(
    TDividendList     *divList,
	char		      settleType,
    TEquitySettlement *stm)
{
    static      char    routine[] = "CrxNewTEqStatData";
    int status = FAILURE;

    TEqStatData *eqStatData = NULL;
    if ((eqStatData = NEW(TEqStatData)) == NULL)
	goto RETURN;

    eqStatData->divList = divList;
    eqStatData->settleType = settleType;
    eqStatData->stm = stm;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        CrxFreeTEqStatData(eqStatData);
        DR_Error("%s: Failed\n", routine);
        return NULL;
    }

    return eqStatData;
}
  

void CrxFreeTEqStatData(TEqStatData *eqStatData)
{
    if (eqStatData == NULL)
	return;

    GtoFreeTDividendList(eqStatData->divList);
    GtoEqSettleFree(eqStatData->stm);

    FREE(eqStatData);
}



/*f---------------------------------------------------------------------
 * Read data from equity.dyn in London format
 * Basis of vol curve is set to 4L (quarterly)
 *
 * Returns SUCCESS/FAILURE.
 */

int
CrxEqDynDataGet(   
    char    *pathdir,	       	/* (I) tmp direcory name (or NULL) */
    char 	*eqDynFnam,		/* (I) file name (equity.dyn)  */
    double  *spotValue,            	/* (O)  */
    double  *correlation,          	/* (O)  */
    TCurve **indxVolCurve)         	/* (O) index volatility curve */
{
    static	char	routine[] = "CrxEqDynDataGet";
    int	status = FAILURE;

    char	fnam[256];
    FILE *fp = NULL;

    TDate baseDate;

    long numVolPts;
    long basis = 4L;
    long dcc = GTO_ACT_365F;

    TCurve *vc = *indxVolCurve = NULL;

    long i;

    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	sprintf(fnam, "%s/%s", pathdir, eqDynFnam);
    } else {
	strcpy(fnam, eqDynFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	DR_Error("%s: can't open `%s'.\n",
		     routine, fnam);
	goto RETURN;
    }
    
    /* Read constants
     */
    READ_DATA(DRL_TDATE_T,  &baseDate,	        "base date");
    READ_DATA(DRL_DOUBLE_T, spotValue,	        "spot value");
    READ_DATA(DRL_DOUBLE_T, correlation,	"correlation");
    READ_DATA(DRL_LONG_T,   &numVolPts,	        "num vol points");

    /* Build vol curve
     */
    if((vc = GtoNewTCurve(baseDate, numVolPts, basis, dcc)) == NULL)
	goto RETURN;

    for(i=0; i<numVolPts; i++)
    {
	if(CrxFScanVType(fp,DRL_TDATE_T, 
			 (void*) &(vc->fArray[i].fDate)) == FAILURE ||
	   CrxFScanVType(fp,DRL_PERCENT_T, 
			 (void*) &(vc->fArray[i].fRate)) == FAILURE)
	{  
	    DR_Error("%s: Cannot read volatility curve data.\n", routine);
	    goto RETURN;
	}
    }
    
    *indxVolCurve = vc;
    
    status = SUCCESS;
    
  RETURN:
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
    {
	    GtoFreeTCurve(vc);
	    DR_Error("%s: failed.\n", routine);
    }

    return(status);
}




/*f---------------------------------------------------------------------
 * Given a swap zero curve, build a funding curve for a given index.
 * - read funding spreads from equity DR wrapper environment
 * - recover MM and par swap rates from swap curve
 * - add spreads to swap rates, and generate the corresponding zero curve
 *
 * Returns SUCCESS/FAILURE.
 */

int
CrxTFundingCurveGet(   
    char *pathdir,			    /* (I) tmp direcory name (or NULL) */
    CrxTDrWrapperData *drWrapper,	/* (I) dr wrapper data */
    TCurve        **indxZC)     /* (O) index zero curve with funding spread */
{
    static	char	routine[] = "CrxTFundingCurveGet";
    int	status = FAILURE;

    static	char	eqStaFnam[] = "equity.sta";
    char	fnam[256];
    FILE *fp = NULL;

    TCurve *zc = drWrapper->fZcCurve; 
    TDate valueDate = zc->fBaseDate;
    int mmDenom = (int)(drWrapper->fMMDenom);
    long mmDayCount = mmDenom==365 ? GTO_ACT_365F : GTO_ACT_360;
    int swapFreq = (int)(drWrapper->fCmsSwMat->swapPayFreq);

    TDateInterval payInterval;
    long swapDcc = (long)(drWrapper->fSwDcc);

    long numRates = 0;
    double *rates = NULL;
    double *prices = NULL;
    TDate *dates = NULL;
    char *names = NULL;
    TDateInterval *rateMats = NULL;


    double temp;

    long i;
    
    *indxZC = NULL;

    /* Open equity.sta 
     */
    if (pathdir != NULL) {
	sprintf(fnam, "%s/%s", pathdir, eqStaFnam);
    } else {
	strcpy(fnam, eqStaFnam);
    }

    if ((fp = fopen(fnam, "r")) == NULL) 
    {
	DR_Error("%s: can't open `%s'.\n",
		      routine, fnam);
	goto RETURN;
    }
    
    if (CrxFScanVType(fp, DRL_LONG_T, &numRates) != SUCCESS)
    {
	DR_Error("%s: can't read number of benchmarks.\n", routine);
	goto RETURN;
    }

    if((dates=NEW_ARRAY(TDate, numRates)) == NULL ||
       (prices=NEW_ARRAY(double, numRates)) == NULL ||
       (names=NEW_ARRAY(char, numRates)) == NULL)
    {
	DR_Error("%s: Error allocating memory.\n", routine);
	goto RETURN;
    }

     if(CrxLilVectArrayFpReadV(fp, 
			       numRates,
			       DRL_TDATEINTERVAL_T,  (void*) &rateMats,
			       DRL_PERCENT_T,  (void*) &rates,
			       DRL_NULL_T) == FAILURE)
    {  
        DR_Error("%s: Cannot read reset dates.\n", routine);
        goto RETURN;
    }

    /* Fill in arrays for zc generating and compute the funding
     * rates with spreads. Rates <=1y are MM, >1y par rates
     */
    for(i=0; i<numRates; i++)
    {
	if( GtoDtFwdAny(valueDate,
			rateMats+i,
			dates+i) == FAILURE)
	    goto RETURN;

	if(dates[i]-valueDate < 367) /* less than a year = MM */
	{
	    names[i] = 'M';
	    if(GtoZerosToSimplePoint(zc,
				     INTERP_METHOD,
				     valueDate,
				     dates[i],
				     mmDayCount,
				     &temp) == FAILURE)
	    goto RETURN;
	}
	else
	{
	    names[i] = 'S';
	    if(GtoFreq2TDateInterval(swapFreq,
				     &payInterval) == FAILURE ||
	       GtoZerosToCouponsPoint(zc,
				      INTERP_METHOD,
				      valueDate,
				      &payInterval,
				      dates[i],
				      swapDcc,
				      GTO_STUB_NONE,
				      FALSE,
				      &temp) == FAILURE)
	    goto RETURN;
	}

	rates[i] += temp;
	prices[i] = 1.;
    }
   
    /* Generate funding curve
     */
    if(((*indxZC) = GtoNPiZC(valueDate,
			     mmDenom,
			     swapFreq,
			     swapDcc,
			     rates,
			     dates,
			     prices,
			     names,
			     GtoSETINSTR,
			     numRates,
			     NULL,
			     NULL,
			     0, /* num futures */
			     0, /* num FRAs */
			     NULL,
			     0.,
			     0,
			     NULL,
			     0,
			     GtoLINEARINTERP, /* for coupons */
			     INTERP_METHOD, /* for zero rates */
			     0,
			     0)) /* no bus day adjust */
       == NULL) goto RETURN;
			 
    status = SUCCESS;

  RETURN:

    FREE(rates);
    FREE(dates);
    FREE(prices);
    FREE(rateMats);
    FREE(names);
    
    if (fp != NULL) fclose(fp);

    if (status != SUCCESS)
	DR_Error("%s: failed.\n", routine);

    return(status);
}




void    dlinterp (long, double *, long, long, double, double);

static int DrLondonForwardPrice(
        double   spotPrice,       /* (I) spot price                          */
        int      NbDiv,           /* (I) Number of dividend / forward prices */
        double   *dividends,      /* (I) Dividends / Forward prices          */
        TDate    *divDates,       /* (I) ex-Dividend dates                   */
        TDate    *payDivDates,    /* (I) Dividend pay dates                  */
        long     *divTypes,       /* (I) Dividend / forward type             */
        int      NbSettle,        /* (I) Nb of settlement dates(<0 if rolling*/
        TDate    *LastTrading,    /* (I) Last trading date of settl period   */
        TDate    *SettleDate,     /* (I) Corresponding settlement date       */
        char     settleType,      /* (I) Settlement type:'F'ixed or 'R'olling*/
	    long	 busDayConv,	  /* (I) Business day conv  		     */
	    char	 *holidayFile,    /* (I) "NONE" for weekends only
                                   *     "No_Weekends" for no adjustments    */
	    TCurve	 *zcCoF,	      /* (I) Funding Zero Curve	             */ 
        TDate    *fwdDates,       /* (I) Forward dates to compute fwd price  */
        int      NbFwd,           /* (I) Total number of time points         */
        double   *fwdPrices);     /* (O) Fwd price at eachtime forward date  */


/**---------------------------------------------------------------
 * Compute forward price given the equity static data and funding curve.
 * Currently two underlying routines are available: 
 * 1) the ALIB GtoEqForwardPrice function
 * 2) based on London code with proper modification for settlement and dividend
 *    dates adjustments.   
 */
int
CrxForwardPriceGen(
	double		spotPrice,	    /* (I) spot asset price   */
	TEqStatData	*eqStatData,	/* (I) Equity static data */
	TCurve		*zcCoF,		    /* (I) Funding zero curve */
	long	    busDayConv,     /* (I) Business Day Conv    */
	char	    *holidayFile,   /* (I) Holiday File       */
	long		numFwd,		    /* (I) Number of fwd points */
	TDate		*fwdDate,	    /* (I) Forward dates to compute fwdPrice
				        	     *     fwdDate[0] must be the spotDate
					             *     where spotPrice is known  */
	double		*fwdPrice)	    /* (O) Forward price      */
{
	static  char routine[] = "CrxForwardPriceGen"; 
	int     status = FAILURE; 

	TDividendList     *divList   = eqStatData->divList; 

	TEqStmPrivate     *stm       = (TEqStmPrivate *)eqStatData->stm;  
	char		  settleType = eqStatData->settleType; 
	int	i;

	long	numSettle; 
	long	numDiv = divList->fNumItems;

	TDate 	*exDivDates  = NULL;
	TDate 	*payDivDates = NULL;
	double	*divAmounts  = NULL;
	long	*divTypes    = NULL;

	if (settleType == 'R')
		numSettle = stm->stmPeriods[0];
	else
		numSettle = stm->lastTradingDates->fNumItems;

	if ((exDivDates  = NEW_ARRAY(TDate,  numDiv)) == NULL ||
	    (payDivDates = NEW_ARRAY(TDate,  numDiv)) == NULL ||
	    (divAmounts  = NEW_ARRAY(double, numDiv)) == NULL ||
	    (divTypes    = NEW_ARRAY(long,   numDiv)) == NULL )
		goto done;

	for (i=0; i<=divList->fNumItems-1; i++){
		exDivDates[i] = divList->fArray[i].exDivDate;
		payDivDates[i] = divList->fArray[i].payDivDate;
		divAmounts[i] = divList->fArray[i].divAmount;
		divTypes[i]   = divList->fArray[i].divType;
	}

	if (DrLondonForwardPrice(spotPrice,
				 numDiv,
				 divAmounts,
				 exDivDates, 
				 payDivDates, 
				 divTypes, 
				 numSettle,
				 stm->lastTradingDates->fArray,
				 stm->settlementDates->fArray,
				 settleType,
				 busDayConv,
				 holidayFile,
				 zcCoF,			 	
				 fwdDate,
				 numFwd,
				 fwdPrice) == FAILURE)
		goto done;


        status = SUCCESS;

done:
	if (status != SUCCESS)
	    DR_Error("%s: failed.\n", routine);

	FREE(exDivDates);
	FREE(payDivDates);
	FREE(divAmounts);
	FREE(divTypes);

    return (status);

}



/*****  DrLondonForwardPrice  **********************************************/
/*
*       Compute forward price at each forward date given the funding
*       curve and dividend curve.  The first point of forward dates
*	    fwdDate[0] should always be the spot date at which spot price 
*	    is given.
*       interpolation to cope with commodities and equity.
*/
static
int DrLondonForwardPrice(
        double   spotPrice,       /* (I) spot price                          */
        int      NbDiv,           /* (I) Number of dividend / forward prices */
        double   *dividends,      /* (I) Dividends / Forward prices          */
        TDate    *divDates,       /* (I) ex-Dividend dates                   */
        TDate    *payDivDates,    /* (I) Dividend pay dates                  */
        long     *divTypes,       /* (I) Dividend / forward type             */
        int      NbSettle,        /* (I) Nb of settlement dates(<0 if rolling*/
        TDate    *LastTrading,    /* (I) Last trading date of settl period   */
        TDate    *SettleDate,     /* (I) Corresponding settlement date       */
        char     settleType,      /* (I) Settlement type:'F'ixed or 'R'olling*/
	    long	 busDayConv,	  /* (I) Business day conv  		     */
	    char	 *holidayFile,    /* (I) "NONE" for weekends only
                                   *     "No_Weekends" for no adjustments    */
	    TCurve	 *zcCoF,	      /* (I) Funding Zero Curve	             */ 
        TDate    *fwdDates,       /* (I) Forward dates to compute fwd price  */
        int      NbFwd,           /* (I) Total number of time points         */
        double   *fwdPrices)      /* (O) Fwd price at eachtime forward date  */
{
#define __DEBUG__
#undef __DEBUG__
	static  char routine[] = "DrLondonForwardPrice";
	int	status = FAILURE; 

        double  y, dt; 	    

        int     i1,        /* Node index of the previous volatility point */
                i2,        /* Node index of the current volatility point  */
                i3,        /* Node index of the next volatility point     */
                i, j, k;  

	int 	numMerge;  /* Total number of combined dates from divDates
			      and fwdDates  */

	TDate	*adjFwdDates = NULL;
	double 	*fullDateZero = NULL;
	double 	*fullSettleDateZero = NULL;
	double 	*fullPrices = NULL;
	TDate 	*fullDates = NULL;
	TDate 	*fullSettleDates = NULL;

	TDate	spotDate;
	double	spotZero=1.0;

	double 	divDateZero, payDateZero;

	/* 
	 * Adjust the dividend for the delay of payment.
	 * If a dividend is paid after the ex-date, its value at 
 	 * ex-date is the PV over the period pay-date to ex-date.
	 * each discrete div is adjusted by this amount here.
	 */
	for (i=0; i <= NbDiv-1; i++)
	{
	    if (payDivDates[i] > divDates[i])
	    {
	      	if ((GtoDiscountDate(divDates[i],
                                     zcCoF,
                                     INTERP_METHOD,
                                     &divDateZero) == FAILURE) ||
	    	    (GtoDiscountDate(payDivDates[i],
                                     zcCoF,
                                     INTERP_METHOD,
                                     &payDateZero) == FAILURE))
                     goto done;
	
		dividends[i] *= payDateZero/divDateZero;
	    }
 	    else
	    if (payDivDates[i] < divDates[i])
	    {
		DR_Error("%s: Dividend pay date (%s) before ex-date (%s).\n",
			  routine, GtoFormatDate(payDivDates[i]),
			  GtoFormatDate(divDates[i]));
		goto done;
	    }
	}

	/* holiday adjustment for the input forward dates */

	if ((adjFwdDates = NEW_ARRAY(TDate,  NbFwd)) == NULL)
	{
	    DR_Error("%s: Error allocating memory for adjFwdDates.\n",
                      routine);
            goto done;
        }

	for (i = 0; i <= NbFwd-1; i++)
	{
	    if (GtoBusinessDay( fwdDates[i],
				busDayConv,
				holidayFile,
				&adjFwdDates[i]) == FAILURE)
		goto done;
	}

	spotDate = adjFwdDates[0];

	/*
	 *  Merge the adjFwdDates and divDates into one single array
	 *  fullDates and sort out in ascending order and remove
	 *  duplicate dates if exist.
	 */

	if ((fullDates = NEW_ARRAY(TDate,  NbFwd + NbDiv)) == NULL)
	{
	    DR_Error("%s: Error allocating memory for fullDates.\n",
                      routine);
            goto done;
        }
	
	j = 0;

	/* only the divDates after spotDate will be included */
	for (i = 0; i <= NbDiv-1; i++)
	{
	    if (divDates[i] >= spotDate)
	    {
	    	fullDates[j] = divDates[i];
		j++;
	    }
	}

	for (i = 0; i<=NbFwd-1; i++, j++)
	    fullDates[j] = adjFwdDates[i];
	
	numMerge = j;

#ifdef  __DEBUG__
	DR_Error("%s: Merged forward dates before sorting\n", routine);
	for (i=0; i<=numMerge-1; i++)
	    DR_Error("\t%d\t%s\n", i+1, GtoFormatDate(fullDates[i]));
#endif	


	if (CrxVTypeVectSort(fullDates,
			     &numMerge,
			     DRL_TDATE_T,
			     TRUE) == FAILURE)
	    goto done;

#ifdef  __DEBUG__
	DR_Error("%s: Merged forward dates after sorting\n", routine);
	for (i=0; i<=numMerge-1; i++)
	    DR_Error("\t%d\t%s\n", i+1, GtoFormatDate(fullDates[i]));
#endif	

        /* 
         *  Calculate settlement date for each forward date.
         */

	if ((fullSettleDates  = NEW_ARRAY(TDate,  numMerge)) == NULL )
	{
            DR_Error("%s: Error allocating memory for fullSettleDates.\n",
                      routine);
            goto done;
        }

        if (settleType == 'R')     /* Rolling settlement */
	{
	    for (i = 0; i <= numMerge-1; i++)                                  
            {        
		if(GtoDateFromBusDaysOffset(fullDates[i],
                                	    (long)NbSettle,
                                	    holidayFile,
                               		    &fullSettleDates[i]) == FAILURE)
            		goto done;
            }  
	}
        else    		   /* Fixed settlement */
        {
	    i1 = i2 = 0;
	    for (k = 0; k < NbSettle; k++)
            {        
		/* 
		 * All forward dates before current LastTrading[k]
		 * will settle on corresponding SettleDate[k].
		 */
		while ((fullDates[i2] <= LastTrading[k]) && (i2 <= numMerge-1))	
                    i2++; 

                for (i = i1; i < i2; i++)
                    fullSettleDates[i] = SettleDate[k];
                        
                i1 = i2;
            }  
                
	    /*
	     * If there are some forward date left after the last 
             * settlement date, just assume it settles immediately.
	     */
            for (i = i1; i <= numMerge-1; i++)        
                fullSettleDates[i] = fullDates[i];
        }  


	/* 
	 * Compute the zero coupon price at each forward date
	 * and forward settlement date. 
	 */

	if ((fullDateZero       = NEW_ARRAY(double, numMerge)) == NULL ||
	    (fullSettleDateZero = NEW_ARRAY(double, numMerge)) == NULL )
	{
	    DR_Error("%s: Error allocating memory for fullDateZero.\n",
                      routine);
            goto done;
        }

	for (i=0; i <= numMerge-1; i++)
	{
	    if ((GtoDiscountDate(fullDates[i],
                                 zcCoF,
                                 INTERP_METHOD,
                                 &fullDateZero[i]) == FAILURE) ||
	        (GtoDiscountDate(fullSettleDates[i],
                                 zcCoF,
                                 INTERP_METHOD,
                                 &fullSettleDateZero[i]) == FAILURE))
                goto done;
		
	    /* discount between fwd date and spotDate
	     * This allows to handle the case where spot price
 	     * is known at forward starting date
	     * or the settlement date of spot trade is different
 	     * from that of value date of zero curves, 
	     * e.g. bond that follows T+1 settlement convention.
	     */
	    if (i==0) spotZero = fullDateZero[0];
	    fullDateZero[i] /= spotZero;
	    fullSettleDateZero[i] /= spotZero;
	}


	/* forward prices at each date of the full list */

	if ((fullPrices  = NEW_ARRAY(double, numMerge)) == NULL) 
	{
	    DR_Error("%s: Error allocating memory for fullPrices.\n",
                      routine);
            goto done;
        }

        /* 
         *  Adjust quoted spot level for settlemnt delay and pay dividend.
	 *  Get zero corresponding to spot settlement date
         */
        fullPrices[0] = spotPrice * fullSettleDateZero[0];

#ifdef SPOTDIVIDEND
        /* Pay today's dividend if discrete */ 
        for (k = 0; k < NbDiv; k++)         
        {        
            if (fullDates[0] == divDates[k])
            {
		switch (divTypes[k]){
		    case GtoDIVIDENDRATE_AMOUNT:    /* Fixed $ dividend */
		    {
		    	fullPrices[0] -= dividends[k];
                    	break;
                    }
                    case GtoDIVIDENDRATE_PERCENT:   /* Fixed dividend yield */
                    {
		    	/* 
		     	 * Dividend yield dividend[k] is applied to the
		     	 * quoted stock price. 
		     	 */
                    	fullPrices[0] -= spotPrice * dividends[k];
                    	break;                                          
                    }
                    default:
		    {
			break;
	            }
                }  /* switch */
            } 
        } 
#endif

	
        /* 
         *	Calculate forward curve.
         */
	
        i1 = i2 = 0;
                          
        for (k = 0; k <= NbDiv-1; k++)
        {        
	   /* 
	    * If current type is a continuous dividend then 
	    * the following are as well 
	    */
           if (divTypes[k] == GtoDIVIDENDRATE_CONTINUOUS)    
		break;                                                  
        
           /* 
	    *  All forward dates before the current dividend date
	    */        
	   while ((fullDates[i2] < divDates[k]) && (i2 < numMerge-1))    
		i2++;     /* It exists as dividend dates are critical dates */

           if (i2 == i1)
           	continue;

           switch (divTypes[k])
           {
           case GtoDIVIDENDRATE_AMOUNT:          /* Fixed $ dividend */
           {
		/*
		 * Forward price is flat between two dividend dates
		 */
                for (i = i1 + 1; i <= i2; i++)
                    fullPrices[i] = fullPrices[i1];      

                /* 
		 * If fullDates[i2] = divDates[k] 
		 * Subtract PV of $ dividend from previous ex-dividend price
		 */
                if (fullDates[i2] == divDates[k])               
                    fullPrices[i2] -= dividends[k] * fullDateZero[i2]; 

                break;
           }
           case GtoDIVIDENDRATE_PERCENT:         /* Fixed dividend yield */
           {
                for (i = i1 + 1; i <= i2; i++)
                    fullPrices[i] = fullPrices[i1];
                        
                if (fullDates[i2] == divDates[k])
                    fullPrices[i2] -= (fullPrices[i2]/fullSettleDateZero[i2])  
				      * dividends[k] * fullDateZero[i2];	

                break;
           }
           default:          /* Linear interpolation for commodity */
           {
                for (i = i1 + 1; i <= i2; i++)
                {        
                    dlinterp ( fullDates[i], 
                               &y, 
                               fullDates[i1], 
                               fullDates[i2], 
                               fullPrices[i1], 
                               dividends[k]);
                                                        
		    /*
		     * Forward prices have to be discounted
		     */
                    fullPrices[i] = y * fullDateZero[i];   
                }  
           }

           }  /* switch */
        
           i1 = i2;

        }  /* for k */


        if (k < NbDiv)      /* There are some continuous dividend left */
        {
           while ((fullDates[i2] < divDates[k]) && (i2 < numMerge-1))
                i2++;

           /* Flat forward from last discrete dividend date */
           for (i = i1 + 1; i <= i2; i++)        
                fullPrices[i] = fullPrices[i1];
        }  

        i3 = 0;

        for (; k < NbDiv; k++)
        {
           if (k == NbDiv - 1) /* Only point left: we fill dates until last */
                i3 = numMerge - 1;
           else                /* Some continuous dividends left afterwards */
           {
                while ((fullDates[i3] < divDates[k+1]) && (i3 < numMerge-1))
                     i3++;    
           }  
                        
	   /* Continuous dividends are not input in % */
           for (i = i2 + 1; i <= i3; i++) {
		if (GtoDayCountFraction(fullDates[i-1],
					fullDates[i],
					GTO_ACT_365F,
					&dt) == FAILURE)
			goto done;
                fullPrices[i] = fullPrices[i-1] / pow (1 + dividends[k], dt);   
	   }

           i2 = i3;
        }  

        /* 
	 * Flat forward stock after the last discrete dividend 
	 * or forward price 
	 */
        for (i = i2 + 1; i <= numMerge-1; i++)     
            fullPrices[i] = fullPrices[i1];
                                   
        /* 
         * Forward value forward curve.
         */

        for (i = 0; i <= numMerge-1; i++)
            fullPrices[i] /= fullSettleDateZero[i];     

	/*
	 * Select the forward prices at desired input forward dates
	 */
	for (k = 0; k <= NbFwd-1; k++)
	for (i = 0; i <= numMerge-1; i++)
	{ 
	    if (adjFwdDates[k] == fullDates[i])
	    {
		fwdPrices[k] = fullPrices[i];
		break;
	    }
	}	

        status = SUCCESS;

done:
	if (status != SUCCESS)
	    DR_Error("%s: failed.\n", routine);

	FREE(fullDates);
	FREE(fullSettleDates);
        FREE(fullPrices);
        FREE(fullDateZero);
        FREE(fullSettleDateZero);
	FREE(adjFwdDates);

        return (status);

#undef __DEBUG__

} 



/*****  dlinterp  ***********************************************************/
/*
*       Linear interpolation between 2 Dates.
*/
void    dlinterp (  long d,
                    double *y,
                    long d0,
                    long d1,
                    double y0,
                    double y1)
{
	long	a, b;
        double  l0, l1;

	GtoDaysDiff(d, d1, GTO_ACT_365, &a);
	GtoDaysDiff(d0, d1, GTO_ACT_365, &b);
        l0 = ((double) a)/((double) b);

	GtoDaysDiff(d0, d, GTO_ACT_365, &a); 
	GtoDaysDiff(d0, d1, GTO_ACT_365, &b);
        l1 = ((double) a)/((double) b);

        *y=l0*y0+l1*y1;

        return;

}  /* dlinterp */




/* REQUIRE should be used to check input conditions. */
/* REQUIRE assumes existence of 'done' label for error cases. */
/* REQUIRE assumes that the function is named 'routine'. */
#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    GtoErrMsg("%s: Required condition (%s) fails!\n",routine,#cond);\
    goto done;\
}} while(0)

/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveDRWRead(FILE *in_fp)
{
    static char routine[] = "CrxBondRepoCurveDRWRead";

    CrxTBondRepoCurve* repoCurve = NULL;
    
    TDate      spotSettleDate;
    long       dcc;
    int        numDates;
    TDate*     dates = NULL;
    double*    rates = NULL;

    /* read data from file */
    if (CrxReadDate (in_fp, &spotSettleDate)) goto done;
    if (CrxReadGenericData (in_fp, "DayCountConv", 
                            (CrxTStringConverter)CrxDayCountConvFromString,
                            (void*)&dcc)) goto done;
    if (CrxReadInt (in_fp, &numDates)) goto done;
    REQUIRE (numDates > 0);
    if ((dates = NEW_ARRAY(TDate, numDates)) == NULL) goto done;
    if ((rates = NEW_ARRAY(double, numDates)) == NULL) goto done;
    if (CrxReadArrays (in_fp, (size_t)numDates, "DF",
                       dates,
                       rates)) goto done;

    repoCurve = CrxBondRepoCurveMake(spotSettleDate,
                                     dcc,
                                     numDates,
                                     dates,
                                     rates);

  done:

    FREE(dates);
    FREE(rates);
    if (repoCurve == NULL)
        GtoErrMsg("%s: Failed!\n", routine);

    return repoCurve;
}

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondRepoCurve
***************************************************************************
*/
void CrxBondRepoCurveDRWWrite(FILE *out_fp, CrxTBondRepoCurve *p)
{
    int i;

    fputs("# spotSettleDate\n", out_fp);
    fprintf(out_fp, "%ld\n", CrxDateToYYYYMMDD(p->spotSettleDate));

    fputs("# dcc\n", out_fp);
    fprintf(out_fp, "%s\n", CrxDayCountConvToString(p->dcc));

    fputs("# numDates\n", out_fp);
    fprintf(out_fp, "%-4d\n", p->numDates);

    fputs("# dates rates\n", out_fp);
    for (i = 0; i < p->numDates; ++i)
    {
        fprintf(out_fp, "%ld %-10.18g\n", CrxDateToYYYYMMDD(p->dates[i]),
                                          p->rates[i]);
    }
}




/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceDRWRead(FILE *in_fp)
{
    static char routine[] = "CrxBondPriceDRWRead";

    CrxTBondPrice* bondPrice = NULL;
    
    TDate      settleDate;
    double     cleanPrice;

    /* read data from file */
    if (CrxReadDate (in_fp, &settleDate)) goto done;
    if (CrxReadDouble (in_fp, &cleanPrice)) goto done;
    
    bondPrice = CrxBondPriceMake(settleDate,
                                 cleanPrice);

  done:

    if (bondPrice == NULL)
        GtoErrMsg("%s: Failed!\n", routine);

    return bondPrice;
}

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondPrice
***************************************************************************
*/
void CrxBondPriceDRWWrite(FILE *out_fp, CrxTBondPrice *p)
{
    fputs("# settleDate\n", out_fp);
    fprintf(out_fp, "%ld\n", CrxDateToYYYYMMDD(p->settleDate));

    fputs("# cleanPrice\n", out_fp);
    fprintf(out_fp, "%.18g\n", p->cleanPrice);
}


/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveDRWRead(FILE *in_fp)
{
    static char routine[] = "CrxBondPriceVolCurveDRWRead";

    CrxTBondPriceVolCurve* volCurve = NULL;
    
    TDate      today;
    int        numDates;
    TDate*     dates = NULL;
    double*    vols = NULL;

    /* read data from file */
    if (CrxReadDate (in_fp, &today)) goto done;
    if (CrxReadInt (in_fp, &numDates)) goto done;
    REQUIRE (numDates > 0);
    if ((dates = NEW_ARRAY(TDate, numDates)) == NULL) goto done;
    if ((vols = NEW_ARRAY(double, numDates)) == NULL) goto done;
    if (CrxReadArrays (in_fp, (size_t)numDates, "DF",
                       dates,
                       vols)) goto done;

    volCurve = CrxBondPriceVolCurveMake(today,
                                        numDates,
                                        dates,
                                        vols);

  done:

    FREE(dates);
    FREE(vols);
    if (volCurve == NULL)
        GtoErrMsg("%s: Failed!\n", routine);

    return volCurve;
}

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondPriceVolCurve
***************************************************************************
*/
void CrxBondPriceVolCurveDRWWrite(FILE *out_fp, CrxTBondPriceVolCurve *p)
{
    int i;

    fputs("# today\n", out_fp);
    fprintf(out_fp, "%ld\n", CrxDateToYYYYMMDD(p->today));

    fputs("# numDates\n", out_fp);
    fprintf(out_fp, "%-4d\n", p->numDates);

    fputs("# dates vols\n", out_fp);
    for (i = 0; i < p->numDates; ++i)
    {
        fprintf(out_fp, "%ld %-10.18g\n", CrxDateToYYYYMMDD(p->dates[i]),
                                          p->vols[i]);
    }
}


/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveDRWRead(FILE *in_fp)
{
    static char routine[] = "CrxBondSpreadVolCurveDRWRead";

    CrxTBondSpreadVolCurve* volCurve = NULL;
    
    TDate      today;
    int        numDates;
    TDate*     dates = NULL;
    double*    vols = NULL;

    /* read data from file */
    if (CrxReadDate (in_fp, &today)) goto done;
    if (CrxReadInt (in_fp, &numDates)) goto done;
    REQUIRE (numDates > 0);
    if ((dates = NEW_ARRAY(TDate, numDates)) == NULL) goto done;
    if ((vols = NEW_ARRAY(double, numDates)) == NULL) goto done;
    if (CrxReadArrays (in_fp, (size_t)numDates, "DF",
                       dates,
                       vols)) goto done;

    volCurve = CrxBondSpreadVolCurveMake(today,
                                         numDates,
                                         dates,
                                         vols);

  done:

    FREE(dates);
    FREE(vols);
    if (volCurve == NULL)
        GtoErrMsg("%s: Failed!\n", routine);

    return volCurve;
}

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondSpreadVolCurve
***************************************************************************
*/
void CrxBondSpreadVolCurveDRWWrite(FILE *out_fp, CrxTBondSpreadVolCurve *p)
{
    int i;

    fputs("# today\n", out_fp);
    fprintf(out_fp, "%ld\n", CrxDateToYYYYMMDD(p->today));

    fputs("# numDates\n", out_fp);
    fprintf(out_fp, "%-4d\n", p->numDates);

    fputs("# dates vols\n", out_fp);
    for (i = 0; i < p->numDates; ++i)
    {
        fprintf(out_fp, "%ld %-10.18g\n", CrxDateToYYYYMMDD(p->dates[i]),
                p->vols[i]);
    }
}

#define YYYYMMDD(x) x>0?CrxTDate2DrDate(x):0

/*f
***************************************************************************
** Writes to DR wrapper file for TCurve
***************************************************************************
*/
void CrxWriteTCurve (FILE* fp, TCurve *tc)
{
    fprintf (fp, "# No information for money market and swap conventions\n");
    fprintf (fp, "# Base date\n%ld\n",
             YYYYMMDD(tc->fBaseDate));
    fprintf (fp, "# Money market basis (unknown)\n360\n");
    fprintf (fp, "# Annual or semi-annual curve (unknown)\nS\n");
    fprintf (fp, "# Year basis for swaps (unknown)\nACT\n");

    if (IS_EQUAL(tc->fBasis, GTO_ANNUAL_BASIS) && 
        tc->fDayCountConv == GTO_ACT_365F)
    {
        int i;

        fprintf (fp, "# No of entries\n%d\n",
                 tc->fNumItems);
        fprintf (fp, "#zero maturity yyyymmdd rates (ACT/365F annual)\n");
        /* no conversions necessary */
        for (i = 0; i < tc->fNumItems; ++i)
        {
            fprintf (fp, "%ld %.12g\n", 
                     YYYYMMDD(tc->fArray[i].fDate),
                     1e2 * tc->fArray[i].fRate);
        }
    }
    else
    {
        /* need to convert all the rates to ACT/365F annually compounded */
        int     numDates;
        TDate  *dates = NULL;
        double *rates = NULL;

        if (GtoCurveDatesAndRates (tc,
                                   GTO_FLAT_FORWARDS,
                                   NULL,
                                   0,
                                   NULL,
                                   NULL,
                                   GTO_ANNUAL_BASIS,
                                   GTO_ACT_365F,
                                   &numDates,
                                   &dates,
                                   &rates) == SUCCESS)
        {
            int i;
            fprintf (fp, "# No of entries\n%d\n", numDates);
            fprintf (fp, "#zero maturity yyyymmdd rates (ACT/365F annual)\n");
            for (i = 0; i < tc->fNumItems; ++i)
            {
                fprintf (fp, "%ld %.12f\n", YYYYMMDD(dates[i]), 1e2*rates[i]);
            }
            FREE (dates);
            FREE (rates);
        }
        else
        {
            fprintf (fp, "# error converting rates to annual, ACT/365F\n");
        }
    }
}













