
#ifndef _ITO33_REGRESSION_OUTPUTS_SET_H_
#define _ITO33_REGRESSION_OUTPUTS_SET_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <string>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/common.h"

#include "ito33/xml/write.h"
#include "ito33/xml/read.h"

#include "ito33/xml/finance/bondlike/bondlikeoutput.h"
#include "ito33/xml/finance/modeloutput.h"

#include "ito33/tests/error_quality.h"

namespace ito33
{


template <class T>
class RegressionOutputSet
{
public:
  typedef T TypeOutput;

  RegressionOutputSet() : m_nBasicOutputChange(0), m_nBasicOutputAddition(0) {}

  bool Read(const std::string& strFile)
  {
    static xml::init s_xmlInit;

    std::ifstream inFile(strFile.c_str(), std::ios::in);
    if (!inFile)
      return false;
    inFile.close();

    xml::tree_parser *parser = new xml::tree_parser(strFile.c_str());

    const xml::node& nodeRoot = parser->get_document().get_root_node();

    xml::node::const_iterator pNode;

    // re-initialization
    m_pNames.clear();
    m_pOutputs.clear();
    m_nBasicOutputAddition = 0;
    m_nBasicOutputChange = 0;

    for(pNode = nodeRoot.begin(); pNode != nodeRoot.end(); pNode++)
    {
      if( *pNode->get_name() == '_' )
      {
        TypeOutput output;     // this 
        if( ito33::XML::Restore(*pNode, output) )
        {
          m_pNames.push_back(pNode->get_name());
          m_pOutputs.push_back(output);
        }
      }
    }
    
    delete parser;

    return true;
  }

  void Write( const std::string& strFile,
              const std::string &strRoot,
              const std::string &strAttrName,
              const std::string &strAtrValue)
  {
    std::ofstream ofs(strFile.c_str());
    XML::RootTag tagRoot(strRoot.c_str(), ofs);
    tagRoot.precision(10);
    tagRoot.Attr(strAttrName.c_str(), strAtrValue.c_str());

    size_t n;
    for(n = 0; n < m_pNames.size(); n++)
    {
      XML::Tag tagOutput(m_pNames[n].c_str(), tagRoot);
      m_pOutputs[n].Dump(tagOutput);
    }
  }

  ErrorQuality CheckModelOutput(const std::string& strName,
                          const TypeOutput& output,
                          shared_ptr<TypeOutput>& pOutputOld)
  {
    // initialisation of output
    ErrorQuality resultQuality = ErrorQuality_pass;
    pOutputOld.reset();

    size_t n;
    for(n = 0;
        n < m_pNames.size() && strName != m_pNames[n];
        n++)
      ;

    if(n < m_pNames.size()) //find same test
    {
      pOutputOld = make_ptr( new TypeOutput(m_pOutputs[n]) );
      resultQuality = ComparePrice(m_pOutputs[n].GetPrice(), output.GetPrice())
                         .strictLevel;
      if( resultQuality == ErrorQuality_fail )
      {
        m_nBasicOutputChange++;

        m_pOutputs[n] = output;
      }
    }
    else
    {
      m_nBasicOutputAddition++;

      m_pNames.push_back(strName);
      m_pOutputs.push_back(output);
    }

    return resultQuality;
  }
  
  bool Changed() { return m_nBasicOutputChange > 0; }

  bool Added() { return m_nBasicOutputAddition > 0; }

  size_t GetNbBasicOutputChange()
  {
    return m_nBasicOutputChange;
  }

private:

  std::vector<std::string>     m_pNames;
  std::vector<TypeOutput> m_pOutputs;

  size_t m_nBasicOutputChange;
  size_t m_nBasicOutputAddition;
};


}

#endif // #define _ITO33_REGRESSION_OUTPUTS_SET_H_
