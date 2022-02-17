#ifndef ESL_ALLOC_DOT_H
#define ESL_ALLOC_DOT_H

/** NOTE: This file should be only included through 'esl_alloc.c'
 */


#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"

/*****  MktVol_Init  **********************************************************/
/**
*       Initialize MKTVOL_DATA structure.
*/
void    MktVol_Init (
		MKTVOL_DATA    *mktvol_data /** Market Vol */
		);

/*****  Opt_Out_Data_Init  **************************************************/
/**
 *		Initialize the OPT_OUT_DATA structure
 */
void     Opt_Out_Data_Init(
		OPT_OUT_DATA  *ood 
		);

/*****  DR_Array  ***********************************************************/
/**
*       Allocation of memory for an array of arbitrary type. Returns void*.
*/
void    *DR_Array ( int     type    /** (I) Type         */ 
                   ,int     nl      /** (I) Lower bound  */
                   ,int     nh      /** (I) Higher bound */
		);

/*****  DR_Matrix  **********************************************************/
/**
 *      Allocation of memory for a matrix of arbitrary type. Returns void *.
 */
void  *DR_Matrix (int     type    /** (I) Type                */
                 ,int     nrl     /** (I) Lower row index     */
                 ,int     nrh     /** (I) Higher row index    */
                 ,int     ncl     /** (I) Lower column index  */
                 ,int     nch     /** (I) Higher column index */
		);

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
		);

/*****  Free_DR_Array  ******************************************************/
/**
*       Free DR array.
*/
int     Free_DR_Array ( void   *Array   /** (I) Array        */
                       ,int     type    /** (I) Type         */
                       ,int     nl      /** (I) Lower bound  */
                       ,int     nh      /** (I) Higher bound */
		);

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
		);

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
		);


#ifdef  __cplusplus
}
#endif


#endif


