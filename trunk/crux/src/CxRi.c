/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

void
CxRiNew(CxtRi *aRi)
{
    aRi->arr = NULL;
    aRi->arrLen = 0;
    aRi->nints = 0;
    aRi->ind = 0;
}

void
CxRiDelete(CxtRi *aRi)
{
    if (aRi->arr != NULL)
    {
	CxmFree(aRi->arr);
    }
}

void
CxRiInit(CxtRi *aRi, uint32_t aNints)
{
    if (aRi->arr == NULL)
    {
	/* Allocate an array with one slot per integer. */
	aRi->arrLen = aNints;
	if (aNints > 0)
	{
	    aRi->arr = (uint32_t *) CxmCalloc(aNints, sizeof(uint32_t));
	}

	/* ind was already initialized to 0, so there is no need to do so
	 * here. */
    }
    else if (aRi->arrLen < aNints)
    {
	/* arr isn't big enough. */

	/* Deallocate old array. */
	CxmFree(aRi->arr);

	/* Allocate new array. */
	aRi->arrLen = aNints;
	aRi->arr = (uint32_t *) CxmCalloc(aNints, sizeof(uint32_t));

	/* Reset ind. */
	aRi->ind = 0;
    }
    else
    {
	uint32_t i;

	/* Undo damage to arr. */
	for (i = 0; i < aRi->ind; i++)
	{
	    CxmAssert(i < aRi->arrLen);
	    CxmAssert(aRi->arr[i] <= aRi->arrLen);
	    if (aRi->arr[i] - 1 >= aRi->ind)
	    {
		aRi->arr[aRi->arr[i] - 1] = 0;
	    }
	    aRi->arr[i] = 0;
	}

	/* Reset ind. */
	aRi->ind = 0;
    }

    /* Finally, set nints. */
    aRi->nints = aNints;

    /* Make sure arr is properly initialized. */
    CxmAssert(aRi->ind == 0);
#ifdef CxmDebug
    if (aRi->arr != NULL)
    {
	uint32_t i;

	for (i = 0; i < aRi->nints; i++)
	{
	    if (aRi->arr[i] != 0)
	    {
		fprintf(stderr, "%s:%d:%s(): Cell %u is %u\n",
			__FILE__, __LINE__, __FUNCTION__, i, aRi->arr[i]);
	    }
	}
    }
#endif
}

uint32_t
CxRiNintsGet(CxtRi *aRi)
{
    return aRi->nints;
}

uint32_t
CxRiIndGet(CxtRi *aRi)
{
    return aRi->ind;
}

uint32_t
CxRiRandomGet(CxtRi *aRi, CxtMt *aMt)
{
    uint32_t retval, r;

    CxmAssert(aRi->nints > 0);

    /* This algorithm is explained in ri.h. */

    /* If all integers have been iterated over, re-initialize. */
    if (aRi->ind == aRi->nints)
    {
	CxRiInit(aRi, aRi->nints);
    }

    /* 2) */
    r = aRi->ind + CxMtUint32RangeGet(aMt, aRi->nints - aRi->ind);
    CxmAssert(r >= aRi->ind);
    CxmAssert(r < aRi->nints);

    /* 3) */
    if (aRi->arr[r] == 0)
    {
	aRi->arr[r] = r + 1;
    }

    /* 4) */
    CxmAssert(aRi->ind < aRi->nints);
    if (aRi->arr[aRi->ind] == 0)
    {
	aRi->arr[aRi->ind] = aRi->ind + 1;
    }

    /* 5) */
    retval = aRi->arr[r];
    aRi->arr[r] = aRi->arr[aRi->ind];
    aRi->arr[aRi->ind] = retval;

    /* 6) */
    aRi->ind++;

    CxmAssert(retval != 0);
    CxmAssert(retval <= aRi->nints);
#ifdef CxmDebug
    {
	uint32_t i;

	for (i = 0; i < aRi->ind - 1; i++)
	{
	    CxmAssert(aRi->arr[i] != retval);
	}
    }
#endif

    /* Decrement, in order to deal with the fact that all numbers in the array
     * are one greater than actual (allows 0 to mean "invalid"). */
    retval--;

    return retval;
}
