/* ======================================================
   FILENAME:  num_h_stat.h
   
   PURPOSE:   a few functions to do basic statistics 
   ====================================================== */

#ifndef NUM_H_STAT_H
#define NUM_H_STAT_H

/*
 * returns stdev of data1, or of its returns or differences depending:
 * if flg=0 returns stdev of log of returns
 * if flg=1 returns stdev of first differences
 * if flg=2 returns stdev
 * if aveptr is not NULL it is set equal to 
 * the average/avg difference/avg return. Almost no error checking
 * is performed.
 * If flg=0 or flg=1 then data1 is overwritten by first differences/log returns
 * for index 0 to len-2
 */
double srt_f_stdev(double *data1, int len, int flg, double *aveptr);

/*
 * computes correlation between data1 and data2.  These are first transformed
 * by calling srt_f_stdev.  It is assumed that the lengths are the same;
 * e.g. this will break if you choose flg1 = 0 or 1 and flg2 = 2. (for example)
 */
double srt_f_corr(double *data1,double *data2,int len, int flg1, int flg2);

#endif
