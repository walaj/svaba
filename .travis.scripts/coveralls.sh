#!/bin/bash

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
#if [ "$CXX" == "g++" ];
#then
  ##sudo pip install cpp-coveralls
  cd test
  ./configure --with-boost=${BOOST_ROOT} ###LDFLAGS=--coverage CXXFLAGS=--coverage
  make
  ./variant_test
  #cpp-coveralls -r ../ -e examples -e doxy -e R -e rtdocs --verbose -t ${COVERALLS_TOKEN}
#fi
