/****************************************************************
 * Module:      VNFM
 * Submodule:   Spread sheet object definition for VnfmData
 * Function:
 * Author:      David Liu, 7/6/1999
 *****************************************************************/
#include <string.h>
#include "macros.h"    /* FREE, FREE_ARRAY */
#include "cerror.h"    /* GtoErrMsg */ 
#include "classreg.h"  /* GtoRegisterClasses */

#include "drlmatrixo.h"	/* Matrix object */

#define _vnfm_SOURCE
#include "vnfmdata.h"


/*-------------------------------------------------------------------
 *                CLASS REGISTRATION RECORD 
 *-------------------------------------------------------------------
 */ 

TObjectType VNFMDATA_OBJECT_TYPE =
{ "VNFMDATA",                 		/* class short name */ 
  NULL,                  	      	/* no base class */
  (TObjectFreeFunc*)   VnfmDataDelete, 	/* destructor */
  NULL,                  		/* no base class to get */ 
  (TObjectGetFunc*)    VnfmDataGetField,/* get method */ 
  NULL,                  		/* no coerceTo function  */ 
  NULL,                  		/* no coerceFrom */  
  NULL,                  		/* no put function */ 
  (TObjectEncodeFunc*) VnfmDataEncode,  /* encode function */
  (TObjectDecodeFunc*) VnfmDataDecode,  /* decode function */
  (TObjectCopyFunc *)  VnfmDataMakeCopy /* copy function */
};


GTO_EXPORT(int)         VnfmObjectRegister ( 
        TObjectManager * om) 
{ 
    static char         routine[]="VnfmObjectRegister"; 
    int                 status=FAILURE; 
    
    if (om IS NULL) 
    {
        GtoErrMsg ("%s: input object manager is NULL\n",routine); 
        goto     done; 
    } 
    
    if (GtoRegisterClasses(om) ISNT SUCCESS ||
        GtoRegisterObjectType (om, &MAT_OBJECT_TYPE) ISNT SUCCESS ||
        GtoRegisterObjectType (om, &VNFMDATA_OBJECT_TYPE) ISNT SUCCESS)
    {
        goto done; /* failure */
    }


    status=SUCCESS; 
    
done: 
    if (status IS FAILURE) 
        GtoErrMsg ("%s: Failed to register class.\n", 
                   routine); 
    return (status) ;
}



GTO_EXPORT(VnfmData*)  VnfmDataCreate (
        TDateL *refDateL,       /* 'D' (I) reference date */
        TDateL *zcDatesL,       /* 'D' (I) array of zero coupon dates */     
        FloatL *zcRatesL,       /* 'F' (I) array of zero coupon rates */    
 
        FloatL *backboneqL,     /* 'F' (I) back bone Q */
        FloatL *betaL,          /* 'F' (I) array of mr coeff */
        FloatL *alphaL,         /* 'F' (I) array of weight coeff */
        TDateL *dateL,          /* 'D' (I) array of dates */
        FloatL *sigmaL,         /* 'F' (I) volatility arrays */
        FloatL *rhoL)           /* 'F' (I) correlation arrays */
{ 
	static char routine[]="VnfmDataCreate";
	int   status = FAILURE;

	VnfmData  *vnfmData = NULL; 
	TCurve    *zcCurve  = NULL;

	/* get zero curve */
        if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
                != SUCCESS) goto done;
#ifndef NO_LOGGING
        GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,
                DRL_TCURVE_FMT_PERCENT));
#endif
 
        /* get model parameters */
        if (VnfmWrapRead(&vnfmData,
                        backboneqL, betaL, alphaL,
                        dateL, sigmaL, rhoL,
                        zcCurve)
                != SUCCESS) goto done;
 
        /* compute coeffs */
        if (VnfmComputeCoeff(vnfmData) != SUCCESS)
                goto done;
 
#ifndef NO_LOGGING
        GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog);\
        VnfmPrintCoeff(vnfmData, vnfmFpLog));
#endif

	/* made it through OK */
	status = SUCCESS;

done:

	if (status IS FAILURE)
	GtoErrMsg("%s: Failed.\n", routine);

	return (vnfmData); 

} 




GTO_EXPORT(int) VnfmDataEncode (
     TEncBuffer *eb, char *objName, VnfmData *that)
{
    return SUCCESS;
}



GTO_EXPORT(VnfmData*) VnfmDataDecode (
    TDecBuffer *db, char *objName)
{
    return NULL;
}



GTO_EXPORT(void) VnfmDataPrint (VnfmData* that)
{
    static char         routine[]="VnfmDataPrint"; 

    VnfmPrintCoeff(that, vnfmFpLog);

}




GTO_EXPORT(VnfmData *) VnfmDataMakeCopy (
        VnfmData  *that)
{ 
    static char routine[]="VnfmDataMakeCopy";
    int   status = FAILURE;

    int	i, j, m,
	rIdx, jIdx,
	nDim = that->fNf;

    int 	nDates = that->fNDates,
		nZDates = that->fZcCurve->fNumItems;

    VnfmData  *copy = NULL; 

    /* Allocate a VnfmData structure */
    copy = VnfmNew(nDates, nDim, nZDates+1);
    if (copy == NULL) goto done;

    copy->fBackBoneQ = that->fBackBoneQ;

    for (j=0; j<=nDim-1; j++){
	copy->fBeta[j] = that->fBeta[j];
	copy->fAlpha[j] = that->fAlpha[j];
    }
	
    for(m=0; m<=nDates-1; m++) {
	copy->fTime[m] = that->fTime[m];
	copy->fDt[m] = that->fDt[m];
	copy->fSigma[m] = that->fSigma[m];
	copy->fDate[m] = that->fDate[m];

	for (i=0; i<=nDim-1; i++){
	    copy->fSigma[i][m] = that->fSigma[i][m];

	    for (j=i+1; j<=nDim-1; j++) {
	    	rIdx = RHOIDX(i, j);
		jIdx = JIDX(i, j);
	    	copy->fRho[rIdx][m] = that->fRho[rIdx][m];
	    	copy->fRho[rIdx][m] = that->fRho[rIdx][m];
		copy->fJ[jIdx][m] = that->fJ[jIdx][m];
	    }
	}
    }

    copy->fZcCurve = GtoCopyCurve(that->fZcCurve);
    if (copy->fZcCurve == NULL) goto done;

    for (m=1; m<=nZDates-1;m++) {
	copy->fZTime[m] = that->fZTime[m];
	copy->fRate[m] = that->fRate[m];
	copy->fZero[m] = that->fZero[m];
    }

    if (copy IS (VnfmData*)NULL)
        GtoErrMsg("%s: Failed.\n", routine);

    status = SUCCESS;

done:
    if (status IS FAILURE)
	GtoErrMsg("%s: Failed.\n", routine);

    return (copy); 

} 



GTO_EXPORT(void)  VnfmDataDelete (
        VnfmData     *that)
{     
    VnfmFree(that);
} 




GTO_EXPORT(TVar *) VnfmDataGetField (
        VnfmData       *that, 
        char           *fieldName)
{ 
    static char         routine[]="VnfmDataGetField"; 
    int   status = FAILURE;

    int	    numField = 11;

    char    **fieldStr = NULL; 
    TVar     *retval = NULL; 

    if((fieldStr = NEW_ARRAY(char*, numField)) == NULL)
	goto done;

    fieldStr[0] = "numdates";
    fieldStr[1] = "dates";
    fieldStr[2] = "zcdates";
    fieldStr[3] = "zcrates";
    fieldStr[4] = "backbone";
    fieldStr[5] = "numfactor";
    fieldStr[6] = "beta";
    fieldStr[7] = "alpha";
    fieldStr[8] = "sigma";
    fieldStr[9] = "rhoij";
    fieldStr[10]= "jintegij";

    
    /**
     ** If the field name is an empty string, then return all the field-names. 
     **/ 
    if (fieldName[0] IS '\0') {
	retval = GtoNewVarStringArray (numField, fieldStr);
        if (retval IS NULL) 
            goto done; 
    }
    else 
    { 
        if (strcmp (fieldName, fieldStr[0]) == 0)
        { 
            retval = GtoNewVarLong (that->fNDates); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp (fieldName, fieldStr[1]) == 0)
        { 
            retval = GtoNewVarDateArray (that->fNDates, that->fDate); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[2]) == 0) 
        {
	    int i;
	    TDate *dates;
	    int	numZDates = that->fZcCurve->fNumItems;

            dates = NEW_ARRAY(TDate, numZDates); 
	    for(i=0; i<=numZDates-1; i++)
		dates[i] = that->fZcCurve->fArray[i].fDate;

            retval = GtoNewVarDateArray (numZDates, dates); 

	    FREE(dates);

            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[3]) == 0) 
        {
	    int i;
	    double *rates;
	    int	numZDates = that->fZcCurve->fNumItems;

            rates = NEW_ARRAY(double, numZDates); 
	    for(i=0; i<=numZDates-1; i++)
		rates[i] = that->fZcCurve->fArray[i].fRate;

            retval = GtoNewVarDoubleArray (numZDates, rates); 

	    FREE(rates);

            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[4]) == 0) 
        {
            retval = GtoNewVarDouble ((1e0 - that->fBackBoneQ) * 0.5e0); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[5]) == 0) 
        {
            retval = GtoNewVarLong (that->fNf); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[6]) == 0) 
        {
            retval = GtoNewVarDoubleArray (that->fNf, that->fBeta); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[7]) == 0) 
        {
            retval = GtoNewVarDoubleArray (that->fNf, that->fAlpha); 
            if (retval IS NULL) 
                goto done; 
        } 
        else if (strcmp(fieldName, fieldStr[8]) == 0) 
        {
            retval = GtoNewVarDoubleArray (that->fNDates, that->fSigma[0]); 
            if (retval IS NULL) 
                goto done; 
	}
        else if (strncmp(fieldName, "rho", 3) == 0) 
        {
	    double *data = NULL;
	    int i, j, m,
	        rIdx;
	    int nDim = NDIM;

            /**
             ** Parse the column number from the 
             ** input field-name
             **/
            char    iNum[2]; 
            char    *jNum = fieldName + 4; 

	    iNum[0] = fieldName[3];
	    iNum[1] = '\0';

            if (iNum IS NULL || jNum IS NULL) 
                goto done;	 /* failed */  

            i = atoi (iNum) - 1;
            j = atoi (jNum) - 1;
	    rIdx = RHOIDX(i,j);	

	    if (i > (that->fNf - 1) || j > (that->fNf - 1) ||
		i<0 || j<0 || i==j){
		status = FAILURE;
                GtoErrMsg("%s: invalid input factor numbers %c%c.\n"
			  "To get correlation between factor 1 and 2, for "
			  "example, use rho12.\n"
			  "i != j, and both factors should not "
			  "exceed the total factor number %d\n", 
			  routine, iNum[0], jNum[0], that->fNf);
		goto done;
	    }
                        
            data = NEW_ARRAY(double, that->fNDates); 
            if (data IS NULL) 
                goto done; /* failed */  

	    for(m=0; m<=that->fNDates-1; m++) {
	        data[m] = that->fRho[rIdx][m];	
	    }

            retval = GtoNewVarDoubleArray (that->fNDates, data); 

	    FREE(data);

            if (retval IS NULL) 
                goto done; 
	}
        else if (strncmp(fieldName, "jinteg", 6) == 0) 
        {
	    double *data = NULL;
	    int i, j, m,
		jIdx;
	    int nDim = NDIM;

            /**
             ** Parse the column number from the 
             ** input field-name
             **/
            char    iNum[2]; 
            char    *jNum = fieldName + 7; 

	    iNum[0] = fieldName[6];
	    iNum[1] = '\0';

            if (iNum IS NULL || jNum IS NULL) 
                goto done;	 /* failed */  

            i = atoi (iNum) - 1;
            j = atoi (jNum) - 1;
	    jIdx = JIDX(i,j);	

	    if (i > (that->fNf - 1) || j > (that->fNf - 1) ||
		i<0 || j<0 || i==j){
		status = FAILURE;
                GtoErrMsg("%s: invalid input factor numbers %c%c.\n"
			  "To get J integral between factor 1 and 2, for "
			  "example, use jinteg12.\n"
			  "i != j, and both factors should not "
			  "exceed the total factor number %d\n", 
			  routine, iNum[0], jNum[0], that->fNf);
		goto done;
	    }
                        
            data = NEW_ARRAY(double, that->fNDates); 
            if (data IS NULL) 
                goto done; /* failed */  

	    for(m=0; m<=that->fNDates-1; m++) {
		data[m] = that->fJ[jIdx][m];	
	    }

            retval = GtoNewVarDoubleArray (that->fNDates, data); 

	    FREE(data);

            if (retval IS NULL) 
                goto done; 
	}
    } 
    status = SUCCESS;

    done:
         if (status IS FAILURE)
            GtoErrMsg("%s: Failed.\n", routine);

	 FREE(fieldStr);
	 
         return retval;

} 

