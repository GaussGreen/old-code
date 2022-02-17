#ifndef SWP_H_GEN_ERR_MESG_H
#define SWP_H_GEN_ERR_MESG_H

/*define global error message constants*/

#define EM_NOERROR NULL
#define EM_BADRNG "Ranges invalid - must be equal  , no single cells"
#define EM_BLANK "An input has been left blank"
#define EM_STRING "An input  , that should not be a string  , is a string"
#define EM_RESVAL "Input can only take on restricted positive discrete values"
#define EM_POSNUM "Argument must be positive or 0"
#define EM_STRINP "Argument must be a string"
#define EM_DATEPB "Problem with date inputs"
#define EM_ERRSTR "ERR or string value in range"
#define EM_BADFILE "Problem with file open or read"
#define EM_YIELD_CURVE "No INTEREST RATE or SPREAD should be NEGATIVE"
#define EM_YC_BLANK "No INTEREST RATE or SPREAD should be left BLANK"
#define EM_YC_STRING "No INTEREST RATE or SPREAD should be a STRING"
#define EM_BASIS "Basis inputs should 1  ,2 or 3"
#define EM_BASIS_BLANK "No BASIS input should be left BLANK"
#define EM_BASIS_STRING "No BASIS input should be a STRING"
#define EM_SWAP "No SWAP input should be NEGATIVE"
#define EM_SWAP_DATES "In SWAP input  , problems with entered dates"
#define EM_SWAP_123 "SWAP Basis must be 1 2 or 3"
#define EM_SWAP_BLANK "No SWAP input should be left BLANK"
#define EM_SWAP_STRING "No SWAP input should be a STRING"
#define EM_CORR "Correlation must be between -1 and 1"
#define EM_BID_OFF "Bid must be < Offer"
#define EM_TENOR "A TENOR has been entered incorrectly"
#define EM_NOSTRING "An input  , that should be a string  , is not"
#define EM_SENS "An error in the sensitivity calculation"
#define EM_TENOR_DATE1 "A shorter TENOR should be added in the list"
#define EM_TENOR_DATE2 "A longer TENOR should be added in the list"
#define EM_NOT_INTEGER "Input must be an integer."
#define EM_TOO_SMALL "Input too small."
#define EM_TOO_BIG "Input too large."
#define EM_TODAY "Dates cannot be before today."
#define EM_POS_RATE "Rates and futures prices must be positive."
#define EM_UNKN_YC "Unknown yc name"

#endif
