/*----------------------------------------------------------------------------*/
/*  CORBA Object Interface													*/ 
/*----------------------------------------------------------------------------*/


#ifndef ARM_CORBA_REQUEST_CALL_I_H
#define ARM_CORBA_REQUEST_CALL_I_H


#include "req.h"


class ARM_CORBA_REQUEST_CALL_i : public virtual ARM_CorbaRequestModule::ARM_CORBA_REQUEST_CALLBOAImpl
{
	public:

	ARM_CORBA_REQUEST_CALL_i(void) {};
   	~ARM_CORBA_REQUEST_CALL_i(void) {};


	virtual CORBA::Long Send(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST& inRequest,
				 ARM_CorbaRequestModule::ARM_CORBA_REQUEST*& rRequest,
				 CORBA::Environment &IT_env = CORBA::IT_chooseDefaultEnv());
};


#ifdef ARM_CORBA_REQUEST_CALL_i_cpp

#include <stdio.h>

#endif	// ARM_CORBA_REQUEST_CALL_i_cpp

#endif	// ARM_CORBA_REQUEST_CALL_I_H'


// EOF %M%
