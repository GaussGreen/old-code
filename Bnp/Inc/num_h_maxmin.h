/*******************************************************************************
$Workfile: num_h_maxmin.h $

SOURCE SAFE HISTORY

$Log: /Sort/Development/SORT/inc/num_h_maxmin.h $
 *
 * 3     8/04/99 5:19p D_desbiez
 * Add some parentheses to lock priority order
 *
 * 2     2/07/99 14:44 M_morrell
 * Fixed SourceSafe history stuff
*******************************************************************************/
#if !defined num_h_maxmin_H
#define num_h_maxmin_H

static const char* num_h_maxminHId =
    "{ID}\t$Logfile: /Sort/Development/SORT/inc/num_h_maxmin.h $\t$Revision: 3 $";
/* $NoKeywords: $ */



#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

#endif /* if num_h_maxmin_H not defined */
