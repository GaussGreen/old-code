//
// file   : ito33/tests/convergence_parameter_value.h
// author : ZHANG Yunzhi
// RCS-ID:  $Id: comparison_parameter_value.h,v 1.2 2004/11/23 14:37:43 zhang Exp $
//

#ifndef _ITO33_TEST_COMPARISON_PARAMETER_VALUE_H_
#define _ITO33_TEST_COMPARISON_PARAMETER_VALUE_H_

#include "ito33/xml/write.h"

namespace ito33
{
  namespace XML 
  {
    /// allow to display spot and a value next to each other
    struct ComparisonParameterValue
    {
      ComparisonParameterValue() {}

      void Init(std::string nameTag,
                std::string firstTag,
                std::string secondTag,
                std::string differenceTag)
      {
        m_nameTag = nameTag;
        m_firstTag = firstTag;
        m_secondTag = secondTag;
        m_differenceTag = differenceTag;
      }

      std::string m_nameTag;
      std::string m_firstTag;
      std::string m_secondTag;
      std::string m_differenceTag;

      std::vector<std::string> m_names;
      std::vector<double> m_firsts;
      std::vector<double> m_seconds;
      std::vector<double> m_diffs;
    };
  
  template <> inline Tag MakeNamedTag(const char *name, 
                     const ComparisonParameterValue &parametervalue, 
                     Tag& parent)
    {
      Tag tag(name, parent);

      tag.Attr("name", parametervalue.m_nameTag.c_str())
         .Attr("first", parametervalue.m_firstTag.c_str())
         .Attr("second", parametervalue.m_secondTag.c_str())
         .Attr("difference", parametervalue.m_differenceTag.c_str());
     
      size_t n;
      for ( n = 0; n < parametervalue.m_names.size(); n++)
      {
          tag.Element("result")
                .Attr("name", (parametervalue.m_names[n]))
                .Attr("first",(parametervalue.m_firsts[n]))
                .Attr("second", (parametervalue.m_seconds[n]))
                .Attr("difference", 
                        (
                          fabs(parametervalue.m_diffs[n]) < 1.e-7 ||
                          ( 
                            parametervalue.m_firsts[n] != 0 && 
                            fabs(parametervalue.m_diffs[n] 
                                  / parametervalue.m_firsts[n]) < 1.e-7
                          )
                          ?
                          0
                          :
                          parametervalue.m_diffs[n]
                        )
                     );
      }

    return tag;
    }
  } //end XML
}

#endif // #define _ITO33_TEST_COMPARISON_PARAMETER_VALUE_H_
