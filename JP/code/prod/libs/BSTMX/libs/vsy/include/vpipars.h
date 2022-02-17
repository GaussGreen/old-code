/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpipar.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpipar_H
#define	_vpipar_H

#include "vpbase.h"
#include "krstbank.h"	// Past rate resets information (KResetBank)




extern	void	VPInsParEval(
	const char *inp,			    // (I) 
	KVector(SharedPointer<KVPAtom>) *vpArrayO,  // (O) array of instruments 
	SharedPointer<KVPInstr> 	*vpRootO,    // (O) root instruments
	KMap(String,double) 	&vpConstTable,	    // (O) Constants table
	KResetBank &vpResetBank,		    // (O) Resets table
	String &vpDiscName);			    // (O) discount curve name


extern	SharedPointer<KVPInstr>	VPInsParEvalFile(
	const char* fnam,		// (I) file name
	KResetBank &vpResetBank,	// (O) Resets table
	String &vpDiscName);		// (O) discount curve name







#endif	/* _vpipar_H */

