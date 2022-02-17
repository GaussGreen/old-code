/****************************************************************************/
/*      Utilities to maintain a claim bank structure.                       */
/****************************************************************************/
/*      ClaimBank.c                                                         */
/****************************************************************************/

/*
$Header: 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"


/*****  Hyb4_CbkInit  ****************************************************/
/*
 *      Initialise CBK contents
 */

void    Hyb4_CbkInit(CLAIM_BANK *CBK)      /* (I/O) the target claim bank   */
{
    Esl_CbkInit(CBK);
}

/*****  Hyb4_CbkAlloc  ****************************************************/
/*
 *      Allocate memory and initialise the claim bank
 *      Returns SUCCESS or FAILURE
 */

int     Hyb4_CbkAlloc
            (CLAIM_BANK *CBK,          /* (I/O) the target claim bank   */
             int         TotNbSlices,  /* (I) Required size of the bank */
             HYB4_TREE_DATA  *tree_data,    /* (I) Tree data structure       */
             int         dimension)    /* (I) Dim of desired slice      */
{
    int     status = FAILURE;   /* Error status */
    int     i;

    /* check if claim bank is initialised */
    if ((CBK == NULL) || (tree_data == NULL)) goto FREE_MEM_AND_RETURN;

    if ((CBK->Slices   != NULL) ||
        (CBK->EvDates  != NULL) ||
        (CBK->ErDates  != NULL) ||
        (CBK->auxSlice != NULL))
    {
        DR_Error("Hyb4_CbkAlloc: claim bank is not empty!\n");
        return (status); /* return here to avoid freeing the non-empty bank */
    }


    if (TotNbSlices <= 0) goto FREE_MEM_AND_RETURN;
    CBK->TotNbSlices = TotNbSlices;

    /* allocate memory for EvalDates and EarliestUseDates */
    CBK->EvDates = (long *)DR_Array(LONG, 0, TotNbSlices-1);
    CBK->ErDates = (long *)DR_Array(LONG, 0, TotNbSlices-1);

    if ((CBK->EvDates == NULL) || (CBK->ErDates == NULL))
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* allocate memory for the aux slice */
    CBK->auxSlice = Hyb4_Alloc_Slice(tree_data, dimension);
    if (CBK->auxSlice == NULL) goto FREE_MEM_AND_RETURN;
  
    /* allocate memory for the slices */
    CBK->Slices = (TSLICE *)DR_Array(DOUBLE_PTR, 0, TotNbSlices-1);

    if (CBK->Slices == NULL) goto FREE_MEM_AND_RETURN;

    for (i=0; i < TotNbSlices; i++)
    {
        (CBK->Slices)[i] = Hyb4_Alloc_Slice(tree_data, dimension);
        if ((CBK->Slices)[i] == NULL) goto FREE_MEM_AND_RETURN;
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    if (status != SUCCESS) 
    {
        DR_Error("Hyb4_CbkAlloc: Unable to allocate memory!\n");
        Hyb4_CbkFree(CBK, tree_data, dimension);
    }
    return (status);

} /* Hyb4_CbkAlloc */


/*****  Hyb4_CbkFree  ****************************************************/
/*
 *      Allocate memory and initialise the claim bank
 *      Returns SUCCESS or FAILURE
 */

int     Hyb4_CbkFree
            (CLAIM_BANK *CBK,           /* (I/O) the target claim bank   */
             HYB4_TREE_DATA  *tree_data,
             int         dimension)     /* (I) Dim of desired slice      */
{
    int     status = FAILURE;   /* Error status */
    int     i;

    if ((CBK == NULL) || (tree_data == NULL)) goto RETURN;

    if (Free_DR_Array(CBK->EvDates,LONG,
                      0,CBK->TotNbSlices-1) == FAILURE) goto RETURN;

    if (Free_DR_Array(CBK->ErDates,LONG,
                      0,CBK->TotNbSlices-1) == FAILURE) goto RETURN;

    Hyb4_Free_Slice(CBK->auxSlice,tree_data, dimension);

    if (CBK->Slices != NULL)
    {
        for (i=0;i<CBK->TotNbSlices;i++)
        {
            Hyb4_Free_Slice(CBK->Slices[i],tree_data, dimension);
        }
        if (Free_DR_Array(CBK->Slices,DOUBLE_PTR,
                          0,CBK->TotNbSlices-1) == FAILURE) goto RETURN;
    }

    /* set counters to zero and pointer to NULL */
    Hyb4_CbkInit(CBK);

    status = SUCCESS;

RETURN:

    if (status != SUCCESS) DR_Error("Hyb4_CbkFree: Failed!\n");
    return (status);

} /* Hyb4_CbkFree */


/*****  Hyb4_CbkSizeFromDL  ******************************************************/
/*
 *      Given a full list of payoff evaluation dates and
 *      a full list of earliest usage dates, it calculates the optimal size
 *      for a claim bank using these dates
 *      NOTE: contents of the datelists WILL BE CHANGED
 *      Returns the optimal size 
 */

int     Hyb4_CbkSizeFromDL
            (long   *EvDL,              /* (I) payoff eval date list      */
             int     NbEvDates,         /* (I) size of the eval date list */
             long   *ErDL,              /* (I) earliest usage date list   */
             int     NbErDates)         /* (I) size of the usage list     */
{
    char    ErrorMsg[MAXBUFF];
    int     EvIdx = NbEvDates - 1;
    int     ErIdx  = NbErDates - 1;
    int     BankSize = 0;
    int     MaxSize = 0;
    int     i;

    /* basic checks */
    if (NbEvDates != NbErDates) return(-999);
    if (NbEvDates == 0) return(0);

    /* error if nb dates > 0 and array == NULL */
    if ((EvDL == NULL) || (ErDL == NULL)) return(-999);

    /* ensure ErDates <= EvDates */
    for (i=0; i<NbEvDates; i++)
    {
        if (ErDL[i] > EvDL[i])
        {
            sprintf (ErrorMsg, 
                     "Hyb4_CbkSizeFromDL: Earliest usage date #%d (%ld) > "
                     "evaluation date (%ld)!\n",
                     i+1, YMDDateFromIRDate(ErDL[i]), YMDDateFromIRDate(EvDL[i]));
            DR_Error (ErrorMsg);
            return(-999);
        }
    }

    /* sort the date lists separately */
    if (SortDateList(NbEvDates, EvDL, NULL) == FAILURE) return(-999);
    if (SortDateList(NbEvDates, ErDL, NULL) == FAILURE) return(-999);
    
    /* to calculate max bank size, only need to exhaust EvDL */
    while (EvIdx>=0)
    {
        if (EvDL[EvIdx] >= ErDL[ErIdx])
        {
            BankSize++;
            if (MaxSize<BankSize) MaxSize = BankSize;
            EvIdx--;
        }
        else
        {
            BankSize--;
            ErIdx--;
        }
    } /* while */

    return(MaxSize);

} /* Hyb4_CbkSizeFromDL */


/*****  Hyb4_CbkCalcSize  ******************************************/
/*
 *      Given a critical date list, pick out one event type and the 
 *      selected supp date for the earliest usage date, and 
 *      calculates the optimal size for a claim bank using these dates
 *      Returns the optimal size 
 */

int     Hyb4_CbkCalcSize
            (int         NbCritDate,          /* CritDate list size */
             CRIT_DATE  *CritDate,            /* critical datelist  */
             int         EventType,           /* eval event type    */
             int         ErDateIdx)           /* suppdate idx       */
{
    int     MaxSize = -999;
    long   *EvDL = NULL;
    long   *ErDL  = NULL;
    int     NbEvDates = 0;
    int     NbErDates = 0;
    int     i;

    if (CritDate == NULL) goto FREE_MEM_AND_RETURN;

    /* extract eval date and earliest use lists */
    for (i=0; i<NbCritDate; i++)
    {
        if (CritDate[i].Type == EventType)
        {
            if (AddDateToList(&NbEvDates,
                              &EvDL, 
                              CritDate[i].CritDate
                             ) == FAILURE) goto FREE_MEM_AND_RETURN;

            if (AddDateToList(&NbErDates,
                              &ErDL, 
                              CritDate[i].SuppDate[ErDateIdx]
                             ) == FAILURE) goto FREE_MEM_AND_RETURN;
        }
    } /* for i */

    /* now calculate the max size */
    MaxSize = Hyb4_CbkSizeFromDL(EvDL,
                            NbEvDates,
                            ErDL,
                            NbErDates);

FREE_MEM_AND_RETURN:

    Free_DR_Array(EvDL, LONG,0,NbEvDates-1);
    Free_DR_Array(ErDL, LONG,0,NbErDates-1);

    return(MaxSize);

} /* Hyb4_CbkCalcSize */


/*****  Hyb4_CbkDev  *****************************************************/
/*
 *      Hyb4_Dev the active elements of the bank
 *      Returns SUCCESS or FAILURE
 */

int     Hyb4_CbkDev
            (CLAIM_BANK *CBK,             /* (I/O) payoff eval date list   */
             long        CurrentDate,     /* (I) Current Date              */
             int         t,               /* (I) Current time point        */
             int         T,               /* (I) Last time point           */
             int         DCurve,          /* (I) Discount curve            */
             int         DMode,           /* (I) Nb of dims of dev required */
             HYB4_DEV_DATA    *dev_data,       /* (I) Hyb4_Dev data structure        */
             HYB4_TREE_DATA   *tree_data)      /* (I) Tree data structure       */
{
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i,j;
    int     minJ;               /* offset for the lowest active slice */
    int     activeOfsLoBnd = 0; /* any active slice offset must be >= this */
    TSLICE  tmpSlice = NULL;    /* vars for swapping */
    long    tmpEvDate;
    long    tmpErDate;

    if ((CBK       == NULL) || 
        (dev_data  == NULL) || 
        (tree_data == NULL)) goto RETURN;


    if (CBK->TotNbSlices == 0)    return(SUCCESS);
    if (CBK->NbActiveSlices == 0) return(SUCCESS);
    if (CBK->Locked) goto RETURN;

    /* compact the bank if necessary */
    if (CurrentDate < CBK->MaxErDate)
    {
        CBK->MaxErDate = -1;  /* reset max EarliestUseDate */

        /* compact the bank and reduce NbActiveSlices */
        for (i=0; i<CBK->NbActiveSlices; i++)
        {
            if (CurrentDate < CBK->ErDates[i])
            {
                /* i is now the lowest obsolete slice offset */

                /* find the lowest active slice offset j (> i) */
                minJ = -1;
                for (j = MAX(i+1,activeOfsLoBnd); 
                     j < CBK->NbActiveSlices; 
                     j++)
                {
                    if (CurrentDate >= CBK->ErDates[j])
                    {
                        minJ = j;
                        break;
                    }
                }

                /* no j found implies all j >= i are obsolete */
                /* thus NbActiveSlice = i                     */
                if (minJ == -1)
                {
                    CBK->NbActiveSlices = i;
                    break;
                }
                else
                {
                    /* swap the contents of minJ and i */
                    tmpSlice  = CBK->Slices[i];
                    tmpEvDate = CBK->EvDates[i];
                    tmpErDate = CBK->ErDates[i];

                    CBK->Slices[i]  = CBK->Slices[minJ];
                    CBK->EvDates[i] = CBK->EvDates[minJ];
                    CBK->ErDates[i] = CBK->ErDates[minJ];

                    CBK->Slices[minJ]  = tmpSlice;
                    CBK->EvDates[minJ] = tmpEvDate;
                    CBK->ErDates[minJ] = tmpErDate;

                    /* next active slice offset must be at least minJ + 1 */
                    activeOfsLoBnd = minJ + 1;
                }
            } /* if i is obsolete */

            /* here, i is always an active slice offset */
            /* update the MaxErDate                   */
            CBK->MaxErDate = MAX(CBK->MaxErDate, CBK->ErDates[i]);

        } /* for i */

    } /* (CurrentDate < CBK->MaxErDate) */


    /* now dev the active slices */
    for (i=0; i<CBK->NbActiveSlices; i++)
    {            
        if (Hyb4_Dev(CBK->Slices[i],
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* for i */

    status = SUCCESS;

RETURN:

    if (status == FAILURE) DR_Error("Hyb4_CbkDev: Failed!\n");
    return (status);

} /* Hyb4_CbkDev */

             
/*****  Hyb4_CbkPopSlice  **********************************************/
/*
 *      Release one slice for external use
 *      The bank is locked for all activities except for Hyb4_CbkPushSlice.
 *      Caller is expected to free the slice unless he/she/it pushes it
 *      back to the bank via Hyb4_CbkPushSlice
 *      Returns NULL if failed, or pointer to free slice
 */

TSLICE  Hyb4_CbkPopSlice(CLAIM_BANK *CBK)
{
    TSLICE tmpSlice = NULL;

    if (CBK == NULL) goto RETURN;
    if (CBK->Locked) goto RETURN;
    if (CBK->TotNbSlices == 0) goto RETURN;
    if (CBK->NbActiveSlices == CBK->TotNbSlices) goto RETURN;

    CBK->Locked = TRUE;
    tmpSlice = CBK->Slices[CBK->TotNbSlices-1];

    /* set the slice pointer to NULL to maintain safe FREE */
    CBK->Slices[CBK->TotNbSlices-1] = NULL;

RETURN:

    /* return the top slice of the bank                    */
    return(tmpSlice);

} /* Hyb4_CbkPopSlice */


/*****  Hyb4_CbkPushSlice  **********************************************/
/*
 *      Push a slice onto the claim bank
 *      Returns SUCCESS or FAILURE
 */

int     Hyb4_CbkPushSlice
            (CLAIM_BANK *CBK,
             TSLICE      Slice,
             long        EvDate,       
             long        ErDate)
{
    int     status = FAILURE;   /* Error status */
    int     i;

    if ((CBK == NULL) || (Slice == NULL)) goto RETURN;

    if (!CBK->Locked) goto RETURN;
    if (CBK->Slices[CBK->TotNbSlices-1] != NULL) goto RETURN;

    if (ErDate > EvDate) goto RETURN;

    /* make sure the EvDate is in strictly decreasing order */
    if (CBK->NbActiveSlices > 0)
    {
        if (CBK->EvDates[0]<=EvDate)
        {
            DR_Error("Hyb4_CbkPushSlice: "
                     "New slice is evaluated "
                     "later than elements in the bank!\n"); 
            goto RETURN;
        }
    }

    /* shift the lowest obsolete slice to the top */
    if (CBK->NbActiveSlices < CBK->TotNbSlices-1)
    {
        CBK->Slices[CBK->TotNbSlices-1] = CBK->Slices[CBK->NbActiveSlices];
    }

    /* nudge up the elements one by one */
    for (i=CBK->NbActiveSlices-1;i>=0;i--)
    {
        CBK->Slices[i+1] = CBK->Slices[i];
        CBK->EvDates[i+1] = CBK->EvDates[i];
        CBK->ErDates[i+1] = CBK->ErDates[i];
    }

    /* add the element to 0 offset */
    CBK->Slices[0]  = Slice;
    CBK->EvDates[0] = EvDate;
    CBK->ErDates[0] = ErDate;
    (CBK->NbActiveSlices)++;

    /* update the internal parameters */
    CBK->MaxErDate = MAX(CBK->MaxErDate, ErDate);

    CBK->Locked = FALSE;

    status = SUCCESS;

RETURN:
    
    if (status == FAILURE) DR_Error("Hyb4_CbkPushSlice: Failed!\n"); 
    return (status);

} /* Hyb4_CbkPushSlice */

/*****  Hyb4_CbkReadSlice  ***************************************************/
/*
 *      Given an EvalDate and a claimbank,
 *      Returns a pointer to the slice with the desired EvalDate
 *      Returns NULL if not found
 *      Note: the returned slice is for READ ONLY, do not overwrite or free.
 */

TSLICE   Hyb4_CbkReadSlice(CLAIM_BANK *CBK,
                      long        EvDate)
{
    int ofs;

    ofs = Hyb4_CbkGetOffset(CBK, EvDate, CbkEXACT);

    if (ofs < 0)
        return(NULL);
    else
        return(CBK->Slices[ofs]);
}



/*****  Hyb4_CbkGetOffset  ***************************************************/
/*
 *      Given an EvalDate and a claimbank,
 *      Returns the offset of the desired slice in the bank
 *      if mode = 
 *          CbkEXACT:  returns -999 if no matching EvalDate is found
 *
 *          CbkLOWER:  returns the nearest offset that is LOWER or equal;
 *                     or -999 if all dates > EvalDate
 *          CbkHIGHER: returns the nearest offset that is HIGHER or equal;
 *                     or -999 if all dates < EvalDate
 */

int     Hyb4_CbkGetOffset(CLAIM_BANK *CBK,
                     long        EvDate,
                     int         mode)
{
    int  offset = -999;

    if (CBK == NULL) goto RETURN;

    offset = GetDLOffset(CBK->NbActiveSlices,
                         CBK->EvDates,
                         EvDate,
                         mode);
RETURN:

    return (offset);

} /* CbkGetOfs */












