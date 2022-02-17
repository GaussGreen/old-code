/*************************************************************************** 
*  Name:        ddlmaker/include/tools.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Vadim Zeitlin                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: tools.h,v 1.2 2005/07/06 09:28:30 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __TOOLS_H__
#define __TOOLS_H__

// output message to stdout
void Output(const char *fmt, ...);

// output the message to stderr and terminates the application
//
// FIXME: error handling is fairly trivial: any error kills the program
void Error(const char *fmt, ...);

#endif
