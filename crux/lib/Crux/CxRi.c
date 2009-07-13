#include "CxRi.h"
#include "../SFMT/SFMT.h"

void
CxRiNew(CxtRi *aRi, sfmt_t *aPrng) {
    CxmCheckPtr(aPrng);

    aRi->prng = aPrng;
    aRi->arr = NULL;
    aRi->arrLen = 0;
    aRi->nints = 0;
    aRi->ind = 0;
}

void
CxRiDelete(CxtRi *aRi) {
    if (aRi->arr != NULL) {
	free(aRi->arr);
    }
}

bool
CxRiInit(CxtRi *aRi, uint32_t aNints) {
    if (aRi->arr == NULL) {
	// Allocate an array with one slot per integer.
	aRi->arrLen = aNints;
	if (aNints > 0) {
	    aRi->arr = (uint32_t *) calloc(aNints, sizeof(uint32_t));
	    if (aRi->arr == NULL) {
		return true;
	    }
	}

	// ind was already initialized to 0, so there is no need to do so here.
    } else if (aRi->arrLen < aNints) {
	// arr isn't big enough.

	// Deallocate old array.
	free(aRi->arr);

	// Allocate new array.
	aRi->arrLen = aNints;
	aRi->arr = (uint32_t *) calloc(aNints, sizeof(uint32_t));
	if (aRi->arr == NULL) {
	    return true;
	}

	// Reset ind.
	aRi->ind = 0;
    } else {
	uint32_t i;

	// Undo damage to arr.
	for (i = 0; i < aRi->ind; i++) {
	    CxmAssert(i < aRi->arrLen);
	    CxmAssert(aRi->arr[i] <= aRi->arrLen);
	    if (aRi->arr[i] - 1 >= aRi->ind) {
		aRi->arr[aRi->arr[i] - 1] = 0;
	    }
	    aRi->arr[i] = 0;
	}

	// Reset ind.
	aRi->ind = 0;
    }

    // Finally, set nints.
    aRi->nints = aNints;

    // Make sure arr is properly initialized.
    CxmAssert(aRi->ind == 0);
#ifdef CxmDebug
    if (aRi->arr != NULL) {
	uint32_t i;
	bool error = false;

	for (i = 0; i < aRi->nints; i++) {
	    if (aRi->arr[i] != 0) {
		error = true;
		fprintf(stderr, "%s:%d:%s(): Cell %u is %u\n",
			__FILE__, __LINE__, __func__, i, aRi->arr[i]);
	    }
	}
	CxmAssert(error == false);
    }
#endif

    return false;
}

uint32_t
CxRiNintsGet(CxtRi *aRi) {
    return aRi->nints;
}

uint32_t
CxRiIndGet(CxtRi *aRi) {
    return aRi->ind;
}

uint32_t
CxRiRandomGet(CxtRi *aRi) {
    uint32_t ret, r;

    CxmAssert(aRi->nints > 0);

    // This algorithm is explained in CxRi.h.

    // If all integers have been iterated over, re-initialize.
    if (aRi->ind == aRi->nints) {
	CxRiInit(aRi, aRi->nints);
    }

    // 2)
    r = aRi->ind + gen_rand64_range(aRi->prng, aRi->nints - aRi->ind);
    CxmAssert(r >= aRi->ind);
    CxmAssert(r < aRi->nints);

    // 3)
    if (aRi->arr[r] == 0) {
	aRi->arr[r] = r + 1;
    }

    // 4)
    CxmAssert(aRi->ind < aRi->nints);
    if (aRi->arr[aRi->ind] == 0) {
	aRi->arr[aRi->ind] = aRi->ind + 1;
    }

    // 5)
    ret = aRi->arr[r];
    aRi->arr[r] = aRi->arr[aRi->ind];
    aRi->arr[aRi->ind] = ret;

    // 6)
    aRi->ind++;

    CxmAssert(ret != 0);
    CxmAssert(ret <= aRi->nints);
#ifdef CxmDebug
    {
	uint32_t i;

	for (i = 0; i < aRi->ind - 1; i++) {
	    CxmAssert(aRi->arr[i] != ret);
	}
    }
#endif

    // Decrement, in order to deal with the fact that all numbers in the array
    // are one greater than actual (allows 0 to mean "invalid").
    ret--;

    return ret;
}
