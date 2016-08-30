#include "json/json.h"
#include "json/json-forwards.h"

#include "SeqLib/MiniRules2.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const std::vector<std::string> m_propertyNames = {"reference", "mapQuality", "isReverseStrand"};

bool ParseFilterObject(const string& filterName, const Json::Value& filterObject) {
  
  // filter object parsing variables
  Json::Value null(Json::nullValue);
  Json::Value propertyValue;
    
  // store results
  map<string, string> propertyTokens;
    
  // iterate over known properties
  vector<string>::const_iterator propertyNameIter = m_propertyNames.begin();
  vector<string>::const_iterator propertyNameEnd  = m_propertyNames.end();
  for ( ; propertyNameIter != propertyNameEnd; ++propertyNameIter ) {
    const string& propertyName = (*propertyNameIter);

    // if property defined in filter, add to token list
    propertyValue = filterObject.get(propertyName, null);
    if ( propertyValue != null ) {
      std::cerr << "[" << propertyName << "] -- " << propertyValue.asString() << std::endl;
      propertyTokens.insert( make_pair(propertyName, propertyValue.asString()) );
    }
  }
  
  // add this filter to engin
  //m_filterEngine.addFilter(filterName);
  
  // add token list to this filter
  //return AddPropertyTokensToFilter(filterName, propertyTokens);
  return true;
}

const string GetScriptContents(string script) {
  
  // open script for reading
  FILE* inFile = fopen(script.c_str(), "rb");
  if ( !inFile ) {
    cerr << "bamtools filter ERROR: could not open script: "
	 << script << " for reading" << endl;
    return string();
  }
    
  // read in entire script contents  
  char buffer[1024];
  ostringstream docStream("");
  while ( true ) {
        
    // peek ahead, make sure there is data available
    char ch = fgetc(inFile);
    ungetc(ch, inFile);
    if( feof(inFile) )
      break;
        
    // read next block of data
    if ( fgets(buffer, 1024, inFile) == 0 ) {
      cerr << "bamtools filter ERROR: could not read script contents" << endl;
      return string();
    }
        
    docStream << buffer;
  }
    
  // close script file
  fclose(inFile);
    
  // import buffer contents to document, return
  return docStream.str();
}

int main(int argc, char** argv) {

  std::string script = "/xchip/gistic/Jeremiah/GIT/SeqLib/test.json";
  const string document = GetScriptContents(script);

  // set up JsonCPP reader and attempt to parse script
  Json::Value root;
  Json::Reader reader;
  if ( !reader.parse(document, root) ) {
    // use built-in error reporting mechanism to alert user what was wrong with the script
    cerr  << "bamtools filter ERROR: failed to parse script - see error message(s) below" << endl;
    return false;     
  }

  // see if root object contains multiple filters
  const Json::Value filters = root["filters"];

  cerr << " root size " << root.size() << std::endl;
  cerr << " filters size " << filters.size() << std::endl;

  // iterate over any filters found
  int filterIndex = 0;
  Json::Value::const_iterator filtersIter = filters.begin();
  Json::Value::const_iterator filtersEnd  = filters.end();
  for ( ; filtersIter != filtersEnd; ++filtersIter, ++filterIndex ) {
    Json::Value filter = (*filtersIter);
            
    // convert filter index to string
    string filterName;
      
    // if id tag supplied
    const Json::Value id = filter["id"];
    if ( !id.isNull() ) 
      filterName = id.asString();


    // use array index 
    else {
      stringstream convert;
      convert << filterIndex;
      filterName = convert.str();
    }

    cerr << "filter " << filterName << std::endl;            
    // create & parse filter 
    bool success = true;
    success &= ParseFilterObject(filterName, filter);
  }

  /// make mini rules
  cerr << " ... maing mini rules " << endl;
  
  SeqLib::MiniRulesCollection mrc(script);


}
