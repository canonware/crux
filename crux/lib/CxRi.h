#ifndef CxRi_h
#define CxRi_h

#include "Cx.h"

typedef struct {
    // Array of integers that is used for randomly iterating over all integers
    // from [0..nints).  The following algorithm is used for the random
    // iteration:
    //
    //   1) Initialize ind <-- 0.
    //
    //   2) Choose a random slot r in [ind,nints).
    //
    //   3) If arr[r] == 0, set arr[r] <-- r + 1.
    //
    //   4) If arr[ind] == 0, set arr[ind] <-- ind + 1.
    //
    //   5) Swap arr[ind] <--> arr[r].  arr[ind] contains (choice+1).
    //
    //   6) Increment ind.
    //
    //   7) If ind < nints, go to step 2.
    //
    // arr can be re-initialized by iterating over the first ind elements and
    // clearing the elements associated with the values stored, then clearing
    // the first ind elements.
    //
    // By using 0 as "invalid", it is possible to allocate zeroed memory, which
    // is typically much faster for large allocations.
    uint32_t *arr;
    // Total number of elements in arr.  This is not necessarily the number of
    // elements that are actually being used.
    uint32_t arrLen;
    // Number of integers to randomly iterate over.  nints <= arrLen.
    uint32_t nints;
    // Number of integers that have been chosen from arr.
    uint32_t ind;
} CxtRi;

// Constructor.
void
CxRiNew(CxtRi *aRi);

// Destructor.
void
CxRiDelete(CxtRi *aRi);

// Initialize the internals of aRi such that [0..aNints) will be randomly
// sampled without replacement.  If aRi was previously in the process of
// sampling, re-initialize the internals.
bool
CxRiInit(CxtRi *aRi, uint32_t aNints);

// Get the number of integers being sampled.
uint32_t
CxRiNintsGet(CxtRi *aRi);

// Get the number of times sampling has occurred since CxRiInit() was called.
uint32_t
CxRiIndGet(CxtRi *aRi);

// Sample without replacement.
uint32_t
CxRiRandomGet(CxtRi *aRi);

#endif // CxRi_h
