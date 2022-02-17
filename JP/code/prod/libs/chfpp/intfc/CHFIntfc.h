#ifndef CHFINTFC_H
#define CHFINTFC_H

/****************************************************************************/
/*      CHFInterface class: interface to CHF prepayment model        .      */
/****************************************************************************/
/*      CHFIntfc.h                                                          */
/****************************************************************************/

#include <stdafx.h> // Prepay.h requires that this is included first
#include "prepay.h"
#include <string>
using namespace std;


#define CHFINTFCDDG 0

class CHFInterface {

public:

	enum PoolType {
        UNKNOWNPOOLTYPE,
		FNMA30,
		GNMA30,
		FNMA15
		// need to add other types
	};

	CHFInterface(
        long today,
		PoolType poolType, 
		double wac, 
		double originalTerm, 
		double originalBalance,    
		string dataDir);

	~CHFInterface();
	
	double calculateCPR(
        int p,
        int t,
		double wala, 
		double remainingBalance, 
        int bufSz,
        double* mtg,         // idx 0 = current
	    double* swap2Y, 
		double* swap10Y   
		);

	void initialize();

    // get historical mortgage rate; placeholder
	double getHistoricalRate(int i);
	
private:

	CHFInterface();

	CForecast* m_cForecast;

    long m_today;
	PoolType m_poolType;	
	double m_wac; 
	double m_originalTerm;          // term at pool origination
	double m_originalBalance;       // balance at origination, e.g. 1.0
	string m_CHFDataDirectory;      // location of CHF data files
};


//class CHFInterfaceError {
//
//public:
//	string m_errorMessage;
//
//	CHFInterfaceError(string msg) { m_errorMessage = msg; }
//
//};

// For use by SRM3


#endif /* CHFINTFC_H */


