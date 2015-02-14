#ifndef ATUPLE_H
#define ATUPLE_H

// define a structure for holding tumor and normal counts for
// Discovar contigs
struct ATuple {
  size_t tum = 0;
  size_t norm= 0;

  // constructor from strings
  ATuple(string ttum, string nnorm) { 
    tum = static_cast<size_t>(stoi(ttum));
    norm = static_cast<size_t>(stoi(nnorm));
  }

  // constructor from ints
  ATuple(int ttum, int nnorm) { 
    tum = static_cast<size_t>(ttum);
    norm = static_cast<size_t>(nnorm);
  }


};

#endif
