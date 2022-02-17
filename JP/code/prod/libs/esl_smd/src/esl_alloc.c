/****************************************************************************/
/*        Memory allocation.	                                            */
/****************************************************************************/
/*        ALLOC.c                                                           */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "esl_alloc.h"
#include "esl_error.h"


/*****  MktVol_Init  **********************************************************/
/**
*       Initialize MKTVOL_DATA structure.
*/
void    MktVol_Init (
		MKTVOL_DATA    *mktvol_data /** Market Vol */
		)
{
	int    i;

	for (i=0; i<3; i++)
	{
        mktvol_data->Beta[i]  = 0e0;
        mktvol_data->Alpha[i] = 0e0;
        mktvol_data->Rho[i]   = 0e0;
	}

    mktvol_data->FilterSpotVolFlag = FALSE;
    mktvol_data->SmoothingFlag     = 'N';
    mktvol_data->TraceFlag         = 'Y';

    return;

}  /* Dev_Init */



/*****  Opt_Out_Data_Init  **************************************************/
/**
 *		Initialize the OPT_OUT_DATA structure
 */
void     Opt_Out_Data_Init(
		OPT_OUT_DATA  *ood 
		)
{
	int  i,j;

	if (ood == NULL) return;

	ood->Option = 0.0;

	for (i=0; i<10; i++) ood->Price[i] = 0.0;

	for (i=0; i<ESL_NB_EVENT; i++) 
        {
            ood->prob_calc[i] = FALSE;

	    for (j=0; j<OPT_OUT_DATA_SIZE; j++)
            {
                 ood->prob_out_data[i].EventDate[j] = 0;
                 ood->prob_out_data[i].EventProb[j] = 0.;
            }
        }

	return;

} /* Opt_Out_Data_Init */

/*****  DR_Array  ***********************************************************/
/**
*       Allocation of memory for an array of arbitrary type. Returns void*.
*/
void    *DR_Array ( int     type    /** (I) Type         */ 
                   ,int     nl      /** (I) Lower bound  */
                   ,int     nh      /** (I) Higher bound */
		)
{
    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v;

            v = (CRIT_DATE *) calloc ((unsigned) (nh-nl+1), sizeof (CRIT_DATE));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            /* Shift pointer and cast to void* */

            return ((void *) (v-nl));

        }         
        case INT:
        {
            int *v;

            v = (int *) calloc ((unsigned) (nh-nl+1), sizeof (int));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_PTR:
        {
            int **v;

            v = (int **) calloc ((unsigned) (nh-nl+1), sizeof (int *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_D_PTR:
        {
            int ***v;

            v = (int ***) calloc ((unsigned) (nh-nl+1), sizeof (int **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG:
        {
            long *v;

            v = (long *) calloc ((unsigned) (nh-nl+1), sizeof (long));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_PTR:
        {
            long **v;

            v = (long **) calloc ((unsigned) (nh-nl+1), sizeof (long *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_D_PTR:
        {
            long ***v;

            v = (long ***) calloc ((unsigned) (nh-nl+1), sizeof (long **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE:
        {
            double *v;

            v = (double *) calloc ((unsigned) (nh-nl+1), sizeof (double));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_PTR:
        {
            double **v;

            v = (double **) calloc ((unsigned) (nh-nl+1), sizeof (double *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_D_PTR:
        {
            double ***v;

            v = (double ***) calloc ((unsigned) (nh-nl+1), sizeof (double **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }
        case CHAR:
        {
            char *v;

            v = (char *) calloc ((unsigned) (nh-nl+1), sizeof (char));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        } 
        case CHAR_PTR:
        {
            char **v;

            v = (char **) calloc ((unsigned) (nh-nl+1), sizeof (char *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in DR_Array!");
            return (NULL);
        }                          
    }  /* switch */

}  /* DR_Array */

/*****  DR_Matrix  **********************************************************/
/**
 *      Allocation of memory for a matrix of arbitrary type. Returns void *.
 */
void  *DR_Matrix (int     type    /** (I) Type                */
                 ,int     nrl     /** (I) Lower row index     */
                 ,int     nrh     /** (I) Higher row index    */
                 ,int     ncl     /** (I) Lower column index  */
                 ,int     nch     /** (I) Higher column index */
		)
{
    int i;

    switch (type)
    {
        case INT:
        {
            int **m;

            m = (int **) DR_Array (INT_PTR, nrl, nrh);

            if (m == NULL)
            {
                DR_Error ("DR_Matrix: allocation failure #1!");
                return (NULL);
            }

            for (i = nrl; i <= nrh; i++)
            {
                m[i] = (int *) DR_Array (INT, ncl, nch);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Matrix: allocation failure #2!");
                    return (NULL);
                }
            }

            return ((void *) m);

        }                          
        case LONG:
        {
            long **m;

            m = (long **) DR_Array (LONG_PTR, nrl, nrh);

            if (m == NULL)
            {
                DR_Error ("DR_Matrix: allocation failure #1!");
                return (NULL);
            }

            for (i = nrl; i <= nrh; i++)
            {
                m[i] = (long *) DR_Array (LONG, ncl, nch);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Matrix: allocation failure #2!");
                    return (NULL);
                }
            }

            return ((void *) m);

        }                          
        case DOUBLE:
        {
            double **m;

            m = (double **) DR_Array (DOUBLE_PTR, nrl, nrh);

            if (m == NULL)
            {
                DR_Error ("DR_Matrix: allocation failure #1!");
                return (NULL);
            }

            for (i = nrl; i <= nrh; i++)
            {
                m[i] = (double *) DR_Array (DOUBLE, ncl, nch);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Matrix: allocation failure #2!");
                    return (NULL);
                }
            }

            return ((void *) m);

        }                          
        default:
        {
            DR_Error ("Program bug: unsupported type in DR_Matrix!");
            return (NULL);
        }                          
    }  /* switch */

} /* DR_Matrix */


/*****  DR_Cube  ****************************************************************/
/**
 *       Allocation of memory for a 3D array of arbitrary type. Returns void *.
 */
void    *DR_Cube (  int     type    /** (I) Type                */
                   ,int     nl      /** (I) Lower bound         */
                   ,int     nh      /** (I) Higher bound        */
                   ,int     nrl     /** (I) Lower row index     */
                   ,int     nrh     /** (I) Higher row index    */
                   ,int     ncl     /** (I) Lower column index  */
                   ,int     nch     /** (I) Higher column index */
		)
{
    int i;
    int j;

    switch (type)
    {
        case INT:
        {
            int ***m;

            m = (int ***) DR_Array (INT_D_PTR, nl, nh);

            if (m == NULL)
            {
                DR_Error ("DR_Cube: allocation failure #1!");
                return (NULL);
            }

            for (i = nl; i <= nh; i++)
            {
                m[i] = (int **) DR_Array (INT_PTR, nrl, nrh);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Cube: allocation failure #2!");
                    return (NULL);
                }

                for (j = nrl; j <= nrh; j++)
                {
                    m[i][j] = (int *) DR_Array (INT, ncl, nch);

                    if (m[i][j] == NULL)
                    {
                        DR_Error ("DR_Cube: allocation failure #3!");
                        return (NULL);
                    }
                }  /* for j */
            }  /* for i */

            return ((void *) m);

        }                          
        case LONG:
        {
            long ***m;

            m = (long ***) DR_Array (LONG_D_PTR, nl, nh);

            if (m == NULL)
            {
                DR_Error ("DR_Cube: allocation failure #1!");
                return (NULL);
            }

            for (i = nl; i <= nh; i++)
            {
                m[i] = (long **) DR_Array (LONG_PTR, nrl, nrh);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Cube: allocation failure #2!");
                    return (NULL);
                }

                for (j = nrl; j <= nrh; j++)
                {
                    m[i][j] = (long *) DR_Array (LONG, ncl, nch);

                    if (m[i][j] == NULL)
                    {
                        DR_Error ("DR_Cube: allocation failure #3!");
                        return (NULL);
                    }
                }  /* for j */
            }  /* for i */

            return ((void *) m);

        }                          
        case DOUBLE:
        {
            double ***m;

            m = (double ***) DR_Array (DOUBLE_D_PTR, nl, nh);

            if (m == NULL)
            {
                DR_Error ("DR_Cube: allocation failure #1!");
                return (NULL);
            }

            for (i = nl; i <= nh; i++)
            {
                m[i] = (double **) DR_Array (DOUBLE_PTR, nrl, nrh);

                if (m[i] == NULL)
                {
                    DR_Error ("DR_Cube: allocation failure #2!");
                    return (NULL);
                }

                for (j = nrl; j <= nrh; j++)
                {
                    m[i][j] = (double *) DR_Array (DOUBLE, ncl, nch);

                    if (m[i][j] == NULL)
                    {
                        DR_Error ("DR_Cube: allocation failure #3!");
                        return (NULL);
                    }
                }  /* for j */
            }  /* for i */

            return ((void *) m);

        }                          
        default:
        {
            DR_Error ("Program bug: unsupported type in DR_Cube!");
            return (NULL);
        }                          
    }  /* switch */

}  /* DR_Cube */


/*****  Free_DR_Array  ******************************************************/
/**
*       Free DR array.
*/
int     Free_DR_Array ( void   *Array   /** (I) Array        */
                       ,int     type    /** (I) Type         */
                       ,int     nl      /** (I) Lower bound  */
                       ,int     nh      /** (I) Higher bound */
		)
{
    /* To avoid warning message */
    nh += 0;

    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v = (CRIT_DATE *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
        case INT:
        {
            int *v = (int *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_PTR:
        {
            int **v = (int **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_D_PTR:
        {
            int ***v = (int ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG:
        {
            long *v = (long *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_PTR:
        {
            long **v = (long **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_D_PTR:
        {
            long ***v = (long ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double *v = (double *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_PTR:
        {
            double **v = (double **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_D_PTR:
        {
            double ***v = (double ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }
        case CHAR:
        {
            char *v = (char *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case CHAR_PTR:
        {
            char **v = (char **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Array!");
            return (FAILURE);
        }                          
    }  /* switch */

    Array = NULL;

    return (SUCCESS);

}  /* Free_DR_Array */



/*****  Free_DR_Matrix  *****************************************************/
/**
*       Free DR matrix.
*/
int     Free_DR_Matrix (void   *Matrix      /** (I) Matrix              */
                       ,int     type        /** (I) Type                */
                       ,int     nrl         /** (I) Lower row index     */
                       ,int     nrh         /** (I) Higher row index    */
                       ,int     ncl         /** (I) Lower column index  */
                       ,int     nch         /** (I) Higher column index */
		)
{
    int i;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int **m = (int **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case LONG:
        {
            long **m = (long **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double **m = (double **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Matrix!");
            return (FAILURE);
        }                          
    }  /* switch */

    Matrix = NULL;

    return (SUCCESS);

}  /* Free_DR_Matrix */



/*****  Free_DR_Cube  *******************************************************/
/**
*       Free 3D DR array.
*/
int     Free_DR_Cube (  void   *Cube    /** (I) Cube                */
                       ,int     type    /** (I) Type                */
                       ,int     nl      /** (I) Lower bound         */
                       ,int     nh      /** (I) Higher bound        */
                       ,int     nrl     /** (I) Lower row index     */
                       ,int     nrh     /** (I) Higher row index    */
                       ,int     ncl     /** (I) Lower column index  */
                       ,int     nch     /** (I) Higher column index */
		)
{
    int i;
    int j;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int ***c = (int ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case LONG:
        {
            long ***c = (long ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case DOUBLE:
        {
            double ***c = (double ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                        
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Cube!");
            return (FAILURE);
        }                          
    }  /* switch */

    Cube = NULL;

    return (SUCCESS);

}  /* Free_DR_Cube */

