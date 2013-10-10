////////////////////////////////////////////////////////////////////////////////
//  Copied from
//  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile
////////////////////////////////////////////////////////////////////////////////

#ifndef JSONFILTER_H_
#define JSONFILTER_H_

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"

bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
   // if the jsonVec is empty, then no JSON file was provided so all
   // events should pass
   if (jsonVec.empty())
   {
      return true;
   }
   bool (* funcPtr) (edm::LuminosityBlockRange const &,
                     edm::LuminosityBlockID const &) = &edm::contains;
   edm::LuminosityBlockID lumiID (event.id().run(), 
                                  event.id().luminosityBlock());
   std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
      std::find_if (jsonVec.begin(), jsonVec.end(),
                    boost::bind(funcPtr, _1, lumiID) );
   return jsonVec.end() != iter;

}

#endif  // JSONFILTER_H_

