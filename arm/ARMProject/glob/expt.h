/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: expt.h,v $
 * Revision 1.9  2004/06/25 14:57:00  ranger
 * added copy constructor
 *
 * Revision 1.8  2004/03/26 07:37:13  ebenhamou
 * added macro to avoid typing line and file..
 *
 * Revision 1.7  2004/03/15 11:14:35  ebenhamou
 * added better truncation error
 *
 * Revision 1.6  2003/11/06 14:48:26  jpriaudel
 * Ajout d'un code d'erreur
 *
 * Revision 1.5  2003/10/14 15:06:59  ebenhamou
 * handle leak for string constructor
 *
 * Revision 1.4  2003/10/02 17:51:02  ebenhamou
 * added accessor to message
 *
 * Revision 1.3  2003/09/26 07:44:03  ebenhamou
 * added log
 *
 *
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : expt.h                                                       */
/*                                                                            */
/* DESCRIPTION : Exceptions classes                                           */
/*                                                                            */
/* DATE        : Tue Jun 11 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _EXPT_H
#define _EXPT_H





#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <string>
/// string support
#ifndef std
	using std::string;
#endif
#include <sstream>


#include "armglob.h"


extern "C"
{
#include "retcode.h"
}

#define RET_CODE	long


/*----------------------------------------------------------------------------*/
/* Error codes by type of error                                               */
/*----------------------------------------------------------------------------*/

#define ERR_RANGE_ARRAY_OVER_UNDER		100
#define ERR_CONDITION_NOT_MEET			101
#define ERR_UNIMP_METHOD_CALL			102
#define ERR_YEAR_TERMS_PAIR				103
#define ERR_YEAR_TERMS_DISCRETE_TIME	104
#define ERR_TIME_INDEX					105
#define ERR_STATE_INDEX					106
#define ERR_MEMORY_ALLOCATION			107
#define ERR_INDEX_OUT_OF_RANGE			108
#define ERR_INVALID_COUPON_FREQ			109 
#define ERR_PAYOFF_CALCULATION_PB		110
#define ERR_BACK_PROPAGATION_PB			111
#define ERR_PRICING_PB			        112
#define ERR_OBJECT_NULL                 113
#define ERR_PERSISTENT_NULL             114
#define ERR_DATE_NOT_VALID				115
#define ERR_CONTAINER_SIZE_NOT_VALID	116
#define ERR_CONTAINER_EMPTY				117
#define ERR_CONTAINER_INDEX				118
#define ERR_INVALID_VECTOR_SIZE			119
#define ERR_INVALID_MODEL   			120 
#define ERR_BAD_OBJECT_CLASS            121
#define ERR_SENSITIVITY_CALCULATION		122
#define ERR_PB_OPEN_FILE				123
#define ERR_SHIFT_MODEL					124
#define ERR_OBJECT_UNK					125
#define ERR_CALC_IMPL_VOL_NIMP			126
#define ERR_PB_COUNTRY_VACATIONS		127
#define ERR_INVALID_OPERATION			128
#define ERR_INVALID_INPUT				129


/*--- Math. error codes ----*/

#define ERR_CALC_DET_PB					200
#define ERR_SQUARE_OR_SIZE_PB			201
#define ERR_MATRIX_INVER_PB				202
#define ERR_MATRIX_LIN_SOLVE_PB			203
#define ERR_TOL_PB						204
#define ERR_FUNCTION_EVAL_OVFL			205
#define ERR_INITIAL_VALUE_PB			206
#define ERR_MAX_ITER_NUM_EXD			207
#define ERR_INVALID_ARGUMENT			208
#define ERR_INVALID_DATA				209
#define ERR_PROPAGATION_PB				210
#define ERR_NEGATIVE_YIELD				211
#define ERR_NEWTON_ROOT_PB				212
#define ERR_BOND_CALC_PB				213
#define ERR_SWAP_CALC_PB				215
#define ERR_FOREX_CALC_PB				216
#define ERR_NUMERICAL_CALL				217
#define ERR_CAPFLOOR_CALC_PB		    218
#define ERR_SWAPLEG_PB	        	    219
#define ERR_CCY_NOT_EQUAL	      	    220



#define ERR_DATA_ACCESS             300



/*----------------------------------------------------------------------------*/







class Exception
{
    private:

	    long		errLine;
	    char*       errFunc;
	    char*		errFile;
	    
	    RET_CODE	errCode;
	    char*		message;

        /// name of the file where it puts full message in case of truncation
	    string itsLongMessageLog;
    
        void Init(void);

    public:
	
       //  Exception(void);

	    Exception(const Exception& e);
	
	    Exception(long line, char* file, RET_CODE code, char* mesg,
		          char* func = NULL);
	
	    /// string version of the Exception object!
	    Exception(long line, char* file, RET_CODE code, const string& mesg,
		          const string& func = "");
	
	    /// accessor only to the text and not the full text
	    // char* GetMessage() const { return message; }
		std::string GetMessage() const { if (!message) return ""; return message; }
	    virtual ~Exception(void);
	    
	    long GetErrCode(void)
	    {
		    return(errCode);
	    }
	    
	    //	JLA. We don't want to use this, since we don't want 
		//	to have long messages to crash ARM
		//
		// void GetErrorMessage(char* msg)
	    // {
		//     sprintf(msg, " ?----> ERROR : %ld : %s, FILE : %s, LINE : %ld", 
		// 	        errCode, message, errFile, errLine);
		//     
		//     char buf[ARM_EXCEPTION_MSG_MAX_SIZE];
		//     
		//     if (errFunc)
		//     {
		// 	    sprintf(buf, ", FUNC : %s", errFunc);
		// 	    
		// 	    strcat(msg, buf);
		//     }  
	    // }

		string	GetErrorString(void)
		{
			char	vErrCode[16];

			string	vErrorMsg(" ?----> ERROR : ");
			sprintf(vErrCode, "%i", errCode);
			vErrorMsg += string(vErrCode) + " : " + message;
			vErrorMsg += string(", FILE : ") + errFile;
			sprintf(vErrCode, "%i", errLine);
			vErrorMsg += string(", LINE : ") + vErrCode;

			return	vErrorMsg;
		}
	
	    ARM_RET_CODE GetRetCode(void)
	    {
		    // char         buf[ARM_EXCEPTION_MSG_MAX_SIZE];
		    ARM_RET_CODE retCode;
		    MEMSET(&retCode, 0, sizeof(ARM_RET_CODE));
		    retCode.code = this->errCode;
		    // GetErrorMessage(buf);
			std::string tmp = GetErrorString(); 
		    strncpy(retCode.msg, tmp.c_str(),ARM_EXCEPTION_MSG_MAX_SIZE);
			retCode.msg[ARM_EXCEPTION_MSG_MAX_SIZE-1]=0; 
		    return(retCode); 
	    }
	    
	    virtual	void DebugPrint(void);

	    /// to set a message that is not too big
	    void SetGoodSizeMessage(const string& initialMessage );
};



class MathException : public Exception
{
	private:

	public:
			// MathException(void){}

			MathException(long line, char* file, RET_CODE code, 
                          char* mesg, char* func = NULL) :
				          Exception(line, file, code, mesg, func)
			{
			}

			void DebugPrint(void)
			{
				Exception::DebugPrint();
			}
};


/// Macro to avoid typing lengthous message!
#define ARM_THROW( errorMessageNb, errorMessage ) \
	throw Exception(__LINE__, __FILE__, errorMessageNb, errorMessage );

class	ArmLogger
{ 
public:	static void logDebugger(const std::string& file,long line,const std::string& msg); 
} ;

#ifndef ARMTHROW
	#define ARMTHROW(errARMTHROW,msgARMTHROW) \
		{ \
		std::stringstream _sstr; _sstr<<msgARMTHROW<<std::endl; \
			ArmLogger::logDebugger(__FILE__,__LINE__,_sstr.str()); \
			throw Exception(__LINE__, __FILE__, errARMTHROW, _sstr.str()); \
		}
#endif // ARMTHROW

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/