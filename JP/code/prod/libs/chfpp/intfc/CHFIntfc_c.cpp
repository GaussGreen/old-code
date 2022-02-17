/****************************************************************************/
/*      C interface for CHFInterface                                 .      */
/****************************************************************************/
/*      CHFIntfc_c.cpp                                                      */
/****************************************************************************/

#include "CHFIntfc_c.h"
#include "CHFIntfc.h"

#include <iostream>
#include <sstream>

extern ostringstream CHFErrorstream; 

/* implement C style functions for SRM3 */

void* CHF_create_handle(
    long today,
	int poolType,
	double wac,				
	double originalTerm,		
	double originalBalance,		
	char* paramaterDirectory)			
{
	CHFInterface *i = NULL;

    CHFInterface::PoolType pt = CHFInterface::UNKNOWNPOOLTYPE;

    /* Agency enum from srm3; can't include this b/c of conflicts
    PP_FNMA,
    PP_FHLMC,
    PP_GOLD,
    PP_GNMAI,
    PP_GNMAII,
    PP_WHOLE,*/
    
    switch(poolType) {

    case 0:
        if (originalTerm <= 180. )
            pt = CHFInterface::FNMA15;
        else 
            pt = CHFInterface::FNMA30;
        break;

    case 3:  // GNMAI
    case 4:  // GNMAII
        pt = CHFInterface::GNMA30;
        break;

    default:
        CHFErrorstream << "CHF_create_handle failed; poolType not supported : " 
            << (int)poolType << std::endl;
        return (void*) NULL;

    }
   
    try {
        i = new CHFInterface( 
            today,
            pt,
            wac,
            originalTerm,
            originalBalance,
            paramaterDirectory);

		i->initialize();
    }

    catch (...) {
        CHFErrorstream << "CHF_create_handle failed " << std::endl;
        if (i) {
            delete i;
            i = NULL;
        }
        
    }
    
	return (void*) i;
}

int CHF_Cpr(
	void* handle,
    int p,
    int t,
	double wala, 
    double remainingBalance,
    int bufSz,
    double* mtg,         // idx 0 = current
    double* swap2, 
    double* swap10 ,    
	double* cpr					
	) 
{
	
	CHFInterface *i = (CHFInterface*) handle;

    *cpr = 0.;

    try {
        /*if (p==0)*/
        *cpr = i->calculateCPR(
            p,
            t,
            wala,
            remainingBalance,
            bufSz,
            mtg,
            swap2,
            swap10
            );
            
    }
    catch (...) {
        CHFErrorstream << "CHF_Cpr failed" 
            << std::endl;
        return FAILURE;
    }

	return SUCCESS;
		
}

int CHF_destroy_handle(void* handle) {

	CHFInterface* i = (CHFInterface*)handle;

	delete i;
	
	return SUCCESS;
}

int CHF_hist_rate(void* handle, int i, double* rate) {

	*rate = 0.05;

    return 1;

//	CHFInterface* intfc = (CHFInterface*)handle;
//
//	// add exception handling
//	reurn intfc->getHistoricalRate(i)


}

void CHF_get_error(char* s, int sz) {
    int i;
    string e = CHFErrorstream.str();
    for (i=0; i<sz; s[i]=(char)NULL, i++);
    for (i=0; i<sz && i<e.size(); i++)
        s[i] = e[i];

    CHFErrorstream.str("");
}