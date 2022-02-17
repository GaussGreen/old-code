#ifndef CHFINTFC_C_H
#define CHFINTFC_C_H

/****************************************************************************/
/*      C interface for CHFInterface                                 .      */
/****************************************************************************/
/*      CHFIntfc_c.h                                                        */
/****************************************************************************/

#ifndef SUCCESS
#define SUCCESS 0
#endif

#ifndef FAILURE
#define FAILURE -1
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
* CHF_create_handle             
*
* Allocate and initialize a handle to a CHFInterface object.
*
*     today :               (I)
*     poolType :            (I)
*     wac (e.g. 0.05)       (I)
*     originalTerm (months) (I)    term at pool origination
*     originalBalance :     (I)    balance at pool origination, e.g. 1.0
*     paramaterDirectory:   (I)
*
* Returns NULL on failure, or valid handle on success.
*/
	void* CHF_create_handle(
		long today,
		int poolType,
		double wac,					
		double originalTerm,		
		double originalBalance,		
		char* parameterDir);					
									

/*****************************************************************************
* CHF_Cpr             
*
* Calculate CPR for some point in the time-line
*
*     handle :              (I)
*     path (for debug log)  (I)
*     sim time (months) :   (I)
*     wala :                (I)
*     remainingBalance :    (I)
*       blah blah
*     cpr:                  (O)
*
 */
	int CHF_Cpr(
		void*,	
        int p,
		int t,
        double wala,
		double remainingBalance,
        int bufSz,
		double* mtg,         // idx 0 = current
	    double* swap2Y, 
		double* swap10Y,
		double* cpr					
		);

/*****************************************************************************
* CHF_cleanup             
*
* Destroy CHF handle
*
*     handle :              (I)
*
*/
	int CHF_destroy_handle(
		void*);

	
/*****************************************************************************
* CHF_hist_mtg_rate             
*
* Get a historical mortgage rate
*
*       handle :            (I)
*       month               (I)
*       rate                (O)
*/
	int CHF_hist_rate(
        void* handle,
        int month, 
        double* rate);
	
/*****************************************************************************
* CHF_get_error             
*
* 
*
*       
*       
*       
*/

    void CHF_get_error( 
        char*, 
        int sz);


#ifdef __cplusplus
}
#endif


#endif /* CHFINTFC_C_H */