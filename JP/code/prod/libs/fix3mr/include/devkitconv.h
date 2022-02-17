/****************************************************************************/
/*      Do all conversions necessary to use the dev-kit                     */
/*      conversion, ...                                                     */
/****************************************************************************/
/*      DEVKITCONV.h                                                        */
/****************************************************************************/

/*
$Header:
*/

#ifndef DEVKITCONV_H
#define DEVKITCONV_H



 /* total nb bytes in char block is this +1 */
#ifndef GTO_MAX_STR_LEN
#define GTO_MAX_STR_LEN    127 
#endif

 /* nb bytes for array count: 1 => max 255 elements in char block array */
#ifndef GTO_NUM_BYTES_FOR_STR_CNT
#define GTO_NUM_BYTES_FOR_STR_CNT    1
#endif

/* positioning of first element of each char block in array */
#ifndef GTO_STR_POS
#define GTO_STR_POS(idx)                          \
            (                                     \
                (idx-1)*(GTO_MAX_STR_LEN+1)       \
                      + GTO_NUM_BYTES_FOR_STR_CNT \
            )
#endif

/* total nb bytes in char block array  */
#ifndef DR_CHAR_BLOCK_ARRAY_SIZE
#define DR_CHAR_BLOCK_ARRAY_SIZE(NbBlocks)    \
            (                                 \
                GTO_STR_POS(NbBlocks +1) -1   \
            )
#endif


 /* alib date structure */
typedef long int TDATE;


/* types for Dev Kit convert usage */
typedef enum
{
    DEVKIT_INT,
    DEVKIT_LONG,
    DEVKIT_DOUBLE,
    DEVKIT_CHAR,
    DEVKIT_TDATE
} DEVKIT_TYPE;


int TDate2DRDate(TDATE    date,         /* (I) Alib date format */
                 long    *DrDate);      /* (O) YYYYMMDD         */

int DRDate2TDate(long     DrDate,       /* (I) YYYYMMDD  */
                 TDATE   *odate);       /* (O) Alib date */

int IsPtrNull(DEVKIT_TYPE    type,      /* (I) Type    */
              void          *ptr);      /* (I) Pointer */


#endif  /* DEVKITCONV_H */
