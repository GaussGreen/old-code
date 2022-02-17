/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: helpertext.cpp,v $
 * Revision 1.1  2003/10/08 16:43:31  ebenhamou
 * Initial revision
 *
 */


/*! \file helpertext.cpp
 *
 *  \brief this is more a configuration file
 *  we made lots of efforts to make this as simple as possible..
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gphelp/helpertext.h"

CC_BEGIN_NAMESPACE( ARM )

const string AdditionalHelpText[] = 
{
	"",
	"",
	"	 ####   #####     ##    #    #  #####   ",
	"	#    #  #    #   #  #   ##   #  #    #  ",
	"	#       #    #  #    #  # #  #  #    #  ",
	"	#  ###  #####   ######  #  # #  #    #  ",
	"	#    #  #   #   #    #  #   ##  #    #	",
	"	 ####   #    #  #    #  #    #  #####	",
	"",
	"",
	"	#####   #####      #    #    #	",
	"	#    #  #    #     #     #  #	",
	"	#    #  #    #     #      ##	",
	"	#####   #####      #      ##	",
	"	#       #   #      #     #  #	",
	"	#       #    #     #    #    #	",
	"",
	"",
	"	#    #  ######  #       #####           ######     #    #       ######	",
	"	#    #  #       #       #    #          #          #    #       #		",
	"	######  #####   #       #    #          #####      #    #       #####	",
	"	#    #  #       #       #####           #          #    #       #		",
	"	#    #  #       #       #               #          #    #       #		",
	"	#    #  ######  ######  #               #          #    ######  ######	",
	"",
	"",
	"",
	" VERSION 1.0:",
	" LAST REVISION of the helper file: March 2004",
	" Eric Benhamou",
	"",
	" Copyright (c) CDC IXIS CM July 2003 Paris",
	"",
	" HotLine Email List : LD-M-GenPricerHelp",
	" This email list should be used to report:",
	"		-bugs with level of importance 1 being not very serious 3 being a serious crash or problem!",
	"		-suggestion and comments on the generic pricer",
	"		-questions",
	"		-comments on the documentation",
	"",
	"",
	"",
	" GENERAL INTRODUCTION:"
	"	A deal description is done by a table with:",
	"		- the upper row with different (case insensitive) non empty strings representing column headers",
	"		- the leftmost column containing dates sorted in increasing order",
	"		- the rightmost column containing cashflow description",
	"",
	" The leftmost column represents event dates.. usually reset dates. Cashflows in the current row",
	" have to be in the future. For past cashflow, you are forced to specify another event date",
	" on an upper row!",
	"",
	"	A cell can refer to other cell(s) by calling them with their column header[ reference ]",
	" where the reference is given by either an integer or by i(or I) or i +/- integer",
	" by convention, i is the current row, hence calling the current row cell with column header",
	" 'Rate' is given by Rate[i].",
	"",
	"	By choice, cell can refer only in the future using the function PV!",
	" reference to the past can be done freely! A reference to a row upward indicates ",
	" that the deal is forward looking while a reference with PV to a row downward indi",
	" cates that the deal is backward looking",
	"",	
	"	A function is called by using the function name followed by an opening parenthesis '(' followed ",
	" by the list of arguments seperated by a comma ',' followed by a closing parenthesis ')'",
	" Arguments can be optional. If not specified, the default argument is taken",
	"",
	"	A date can either be an excel date or coded with the format d(d)mmmyyy(yy) with",
	" d(d) being the day, mmm the month, yy(yy) the year",
	" if the year is specified as yy, by convention it is 20yy!",
	" month being one of the following:",
	" Jan, Feb, Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec",
	" Beware that there is no separator between d(d)mmyy(yy) to avoid confusion with other symbols!",
	" Beware that excel does not make the difference between dates and double",
	" Hence some double may be taken as a date. In this particular case, you can just",
	" type the double * 1.0 and this will enforce that the double is recognised as a double"
	"",
	"	Last but not least, the deal description is case insensitive to avoid confusion."
	" ",
	"	PS: For advanced payoff programming in the world of Grand Prix, you should",
	" take advantage of the pre-compiled shared node technology. To do so,",
	" be aware that shared nodes are only the ones that has the same parent node's date.",
	" This comes from the inplace tecbnology on Grand Prix functor!",
	" This means that in order to share a node with different parent dates,",
	" you will need to point to an intermediate cell. If you do not know what you",
	" do or what this means, please do not touch this and ask for assistance.",
	" ",
};

const size_t AdditionalHelpTextSize = sizeof(AdditionalHelpText)/sizeof(AdditionalHelpText[0]);


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

