/* =============================================================================

   FILENAME:    swp_h_datelist.h

   PURPOSE:     Geeration of a swap schedule

   =============================================================================
 */
#ifndef SWP_H_DATELIST_H
#define SWP_H_DATELIST_H

typedef enum { BROKEN, NOTBROKEN, LASTDATELISTTYPE } DateListType;

typedef struct {
  long *date;        /*** array of dates from start to end **/
  int len;           /** length of *date 	***/
  SwapDateDir dir;   /** FWD or BKWD 	**/
  DateListType type; /** BROKEN or NOTBROKEN **/
  Date prev;         /** date before 2nd date  , if BROKEN **/
} DateList, SrtDateList;

SrtDateList new_DateList(int len);
SrtDateList SwapDP_to_DateList(SrtSwapDP *swap, SrtBusDayConv conv);
SrtDateList DateListafter(SrtDateList dl, Date after);

Err swp_f_free_in_DateList(SrtDateList datelist);

#endif
