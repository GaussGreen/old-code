/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/exoticoption/asianoption.h
// Purpose:     Names of elements and attributes used in XML for asian option
// Created:     March 29, 2005
// RCS-ID:      $Id: asianoption.h,v 1.5 2006/05/27 15:17:51 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/exoticoption/asianoption.h
    @brief Contains the names of the elements used in the XML description of
           Asian option.
 */

#ifndef _ITO33_XML_FINANCE_EXOTICOPTION_ASIANOPTION_H_
#define _ITO33_XML_FINANCE_EXOTICOPTION_ASIANOPTION_H_


/**
   @name Tag name macros
*/
//@{

#define XML_TAG_ASIAN_OPTION_ROOT                 "asian_option"
#define XML_TAG_ASIAN_OPTION_CURRENT_AVERAGE      "current_average"
#define XML_TAG_ASIAN_OPTION_AVG_START_DATE       "average_start_date"
#define XML_TAG_ASIAN_OPTION_AVG_END_DATE         "average_end_date"   
#define XML_TAG_ASIAN_OPTION_NB_SAMPLING_AVERAGES "number_of_sampling_averages"
#define XML_TAG_ASIAN_OPTION_NB_SAMPLES_USED      "number_of_samples_used"
//@}

#endif // _ITO33_XML_FINANCE_EXOTICOPTION_ASIANOPTION_H_
