//
// file   : ito33/tests/convergence_parameter_value.h
// author : ZHANG Yunzhi
// RCS-ID:  $Id: convergence_parameter_value.h,v 1.1 2004/09/17 18:02:15 zhang Exp $
//

#ifndef _ITO33_TEST_CONVERGENCE_PARAMETER_VALUE_H_
#define _ITO33_TEST_CONVERGENCE_PARAMETER_VALUE_H_

#include "ito33/xml/write.h"

namespace ito33
{
  namespace XML 
  {
    /// allow to display spot and a value next to each other
    struct ConvergenceParameterValue
    {
      ConvergenceParameterValue() {}

      void init(std::string priceTag,
                std::string diffTag,
                std::string ratingTag)
      {
        m_priceTag = priceTag;
        m_diffTag = diffTag;
        m_ratingTag = ratingTag;
      }

      std::string m_priceTag;
      std::string m_diffTag;
      std::string m_ratingTag;

      std::vector<double> m_prices;
      std::vector<double> m_diffs;
      std::vector<double> m_ratings;
    };
  
  template <> inline Tag MakeNamedTag(const char *name, 
                     const ConvergenceParameterValue &parametervalue, 
                     Tag& parent)
    {
      Tag tag(name, parent);
     
      size_t n;
      for ( n = 0; n < parametervalue.m_prices.size(); n++)
      {
        if(n == 0)
          tag.Element("value")
                .Attr("price",(parametervalue.m_prices[n]));
        else if (n == 1)
          tag.Element("value")
                .Attr("price",(parametervalue.m_prices[n]))
                  .Attr("diff",(parametervalue.m_diffs[n]));
        else
          tag.Element("value")
                .Attr("price",(parametervalue.m_prices[n]))
                  .Attr("diff",(parametervalue.m_diffs[n]))
                    .Attr("rating", (parametervalue.m_ratings[n]));
      }

    return tag;
    }
  } //end XML
}

#endif // #define _ITO33_TEST_CONVERGENCE_PARAMETER_VALUE_H_
