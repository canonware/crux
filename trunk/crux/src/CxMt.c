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

/* This implementation of the Mersenne Twister random number generator is based
 * on the reference implementation made available by the originators.  The
 * original copyright and license follow:
 *
 ******************************************************************************
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.                          
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. The names of its contributors may not be used to endorse or promote 
 *      products derived from this software without specific prior written 
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

/* Constant vector a. */
#define CxmMtMatrixA 0x9908b0dfUL

/* Most significant w-r bits. */
#define CxmMtUpperMask 0x80000000UL

/* Least significant r bits. */
#define CxmMtLowerMask 0x7fffffffUL

/* Initializes mt[CxmMtN] with a seed. */
CxmpInline void
CxpMtInitGenRand(CxtMt *aMt, uint32_t aS)
{
    uint32_t *mt = aMt->mt;
    uint32_t mti;

    mt[0]= aS & 0xffffffffUL;
    for (mti=1; mti<CxmMtN; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }

    aMt->mti = mti;
}

/* generates a random number on [0,0xffffffff]-interval */
CxmpInline uint32_t
CxpMtUint32Get(CxtMt *aMt)
{
    uint32_t *mt = aMt->mt;
    uint32_t mti = aMt->mti;
    uint32_t y;
    static uint32_t mag01[2]={0x0UL, CxmMtMatrixA};
    /* mag01[x] = x * CxmMtMatrixA  for x=0,1 */

    CxmCheckPtr(aMt);
    CxmAssert(aMt->magic == CxmMtMagic);
    CxmAssert(aMt->mti != CxmMtN + 1);

    if (mti >= CxmMtN) { /* generate CxmMtN words at one time */
        int32_t kk;

        if (mti == CxmMtN+1)   /* if CxpMtInitGenRand() has not been called, */
            CxpMtInitGenRand(aMt, 5489UL); /* a default initial seed is used */

        for (kk=0;kk<CxmMtN-CxmMtM;kk++) {
            y = (mt[kk]&CxmMtUpperMask)|(mt[kk+1]&CxmMtLowerMask);
            mt[kk] = mt[kk+CxmMtM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<CxmMtN-1;kk++) {
            y = (mt[kk]&CxmMtUpperMask)|(mt[kk+1]&CxmMtLowerMask);
            mt[kk] = mt[kk+(CxmMtM-CxmMtN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[CxmMtN-1]&CxmMtUpperMask)|(mt[0]&CxmMtLowerMask);
        mt[CxmMtN-1] = mt[CxmMtM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    aMt->mti = mti;
    return y;
}

CxmpInline int64_t
CxpMtSint64Get(CxtMt *aMt)
{
    int64_t retval;

    retval = (int64_t) CxpMtUint32Get(aMt);
    retval &= 0xfffffffeU; /* Clear least significant bit. */
    retval <<= 31;
    retval |= ((int64_t) CxpMtUint32Get(aMt));

    return retval;
}

void
CxMtNew(CxtMt *aMt)
{
    CxmCheckPtr(aMt);

    memset(aMt, 0, sizeof(CxtMt));

    aMt->mti = CxmMtN + 1;

#ifdef CxmDebug
    aMt->magic = CxmMtMagic;
#endif
}

void
CxMtDelete(CxtMt *aMt)
{
    CxmCheckPtr(aMt);
    CxmAssert(aMt->magic == CxmMtMagic);

#ifdef CxmDebug
    memset(aMt, 0x5a, sizeof(CxtMt));
#endif
}

void
CxMtUint32Seed(CxtMt *aMt, uint32_t aSeed)
{
    CxpMtInitGenRand(aMt, aSeed);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void
CxMtStringSeed(CxtMt *aMt, char *aSeed, uint32_t aSeedLen)
{
    uint32_t *mt = aMt->mt;
    int32_t i, j, k;
    uint32_t *init_key = (uint32_t *) aSeed;
    uint32_t key_length = aSeedLen / sizeof(uint32_t);

    CxmCheckPtr(aMt);
    CxmAssert(aMt->magic == CxmMtMagic);
    CxmCheckPtr(aSeed);
    CxmAssert(aSeedLen >= sizeof(uint32_t));

    CxpMtInitGenRand(aMt, 19650218UL);
    i=1; j=0;
    k = (CxmMtN>key_length ? CxmMtN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=CxmMtN) { mt[0] = mt[CxmMtN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=CxmMtN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=CxmMtN) { mt[0] = mt[CxmMtN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

void
CxMtStateGet(CxtMt *aMt, int32_t *rMtI,
	     uint32_t **rMt, uint32_t *rMtLen)
{
    CxmCheckPtr(aMt);
    CxmAssert(aMt->magic == CxmMtMagic);

    *rMtI = aMt->mti;
    *rMt = aMt->mt;
    *rMtLen = CxmMtN;
}

void
CxMtStateSet(CxtMt *aMt, int32_t aMtI,
	     const uint32_t *aMtA, uint32_t aMtLen)
{
    CxmCheckPtr(aMt);
    CxmAssert(aMt->magic == CxmMtMagic);
    CxmAssert(aMtLen == CxmMtN);

    aMt->mti = aMtI;
    memcpy(aMt->mt, aMtA, sizeof(uint32_t) * aMtLen);
}

int64_t
CxMtSint64Get(CxtMt *aMt)
{
    return CxpMtSint64Get(aMt);
}

int64_t
CxMtSint64RangeGet(CxtMt *aMt, int64_t aRange)
{
    int64_t retval, above;

    above = 0x7fffffffffffffffLL - (0x7fffffffffffffffLL % aRange);
    while (1)
    {
	retval = CxpMtSint64Get(aMt);
	if (retval < above)
	{
	    retval %= aRange;
	    break;
	}
    }

    return retval;
}

/* generates a random number on [0,0x7fffffff]-interval */
int32_t
CxMtSint32Get(CxtMt *aMt)
{
    return (int32_t) (CxpMtUint32Get(aMt)>>1);
}

int32_t
CxMtSint32RangeGet(CxtMt *aMt, int32_t aRange)
{
    int32_t retval, above;

    above = 0x7fffffffL - (0x7fffffffL % aRange);
    while (1)
    {
	retval = (int32_t) (CxpMtUint32Get(aMt)>>1);
	if (retval < above)
	{
	    retval %= aRange;
	    break;
	}
    }

    return retval;
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t
CxMtUint32Get(CxtMt *aMt)
{
    return CxpMtUint32Get(aMt);
}

uint32_t
CxMtUint32RangeGet(CxtMt *aMt, uint32_t aRange)
{
    uint32_t retval, above;

    above = 0xffffffffLU - (0xffffffffLU % aRange);
    while (1)
    {
	retval = CxpMtUint32Get(aMt);
	if (retval < above)
	{
	    retval %= aRange;
	    break;
	}
    }

    return retval;
}

/* generates a random number on [0,1]-real-interval */
double
CxMtReal1Get(CxtMt *aMt)
{
    return CxpMtUint32Get(aMt)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double
CxMtReal2Get(CxtMt *aMt)
{
    return CxpMtUint32Get(aMt)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double
CxMtReal3Get(CxtMt *aMt)
{
    return (((double)CxpMtUint32Get(aMt)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double
CxMtRes53Get(CxtMt *aMt)
{ 
    uint32_t a=CxpMtUint32Get(aMt)>>5, b=CxpMtUint32Get(aMt)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
