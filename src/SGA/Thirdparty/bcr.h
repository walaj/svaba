/* The MIT License

   Copyright (c) 2012 Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef BCR_H
#define BCR_H

#include <stdint.h>

struct bcr_s;
typedef struct bcr_s bcr_t;

struct bcritr_s;
typedef struct bcritr_s bcritr_t;

extern int bcr_verbose;

#ifdef __cplusplus
extern "C" {
#endif

	bcr_t *bcr_init(int is_threaded, const char *tmpfn);
	void bcr_destroy(bcr_t *b);
	void bcr_append(bcr_t *b, int len, uint8_t *seq);
	void bcr_build(bcr_t *b);

    // Returns the sequence index of the i-th lexicographically ordered read.
    // Only valid after bcr_build has been called.
    // Added by jts
    int64_t bcr_getLexicographicIndex(bcr_t* b, int64_t i);

	bcritr_t *bcr_itr_init(const bcr_t *b);
	const uint8_t *bcr_itr_next(bcritr_t *itr, int *l);

#ifdef __cplusplus
}
#endif

#endif
