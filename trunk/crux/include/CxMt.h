//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

// This implementation of the Mersenne Twister random number generator is based
// on the reference implementation made available by the originators.  The
// original copyright and license follow:
//
//==============================================================================
//
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.                          
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//==============================================================================

typedef struct CxsMt CxtMt;

// Period parameters.
#define CxmMtN 624
#define CxmMtM 397

// Mersenne Twister random number generator state.
struct CxsMt
{
#ifdef CxmDebug
    uint32_t magic;
#define CxmMtMagic 0x89326918
#endif

    // The array for the state vector.
    uint32_t mt[CxmMtN];
    int32_t mti;
};

// Mt.

// Constructor.
void
CxMtNew(CxtMt *aMt);

// Destructor.
void
CxMtDelete(CxtMt *aMt);

// Seed the random number generator.
void
CxMtUint32Seed(CxtMt *aMt, uint32_t aSeed);

// Seed the random number generator.  Up to 2496 bytes are used.
void
CxMtStringSeed(CxtMt *aMt, char *aSeed, uint32_t aSeedLen);

// Get the current internal state.
void
CxMtStateGet(CxtMt *aMt, int32_t *rMtI,
	     uint32_t **rMt, uint32_t *rMtLen);

// Set the internal state.
void
CxMtStateSet(CxtMt *aMt, int32_t aMtI,
	     const uint32_t *aMta, uint32_t aMtLen);

// Get a random integer in the range [0,2^63).
int64_t
CxMtSint64Get(CxtMt *aMt);

// Get a random integer in the range [0,aRange).
int64_t
CxMtSint64RangeGet(CxtMt *aMt, int64_t aRange);

// Get a random integer in the range [0,0x7fffffff].
int32_t
CxMtSint32Get(CxtMt *aMt);

// Get a random integer in the range [0,aRange).
int32_t
CxMtSint32RangeGet(CxtMt *aMt, int32_t aRange);

// Get a random integer in the range [0,0xffffffff].
uint32_t
CxMtUint32Get(CxtMt *aMt);

// Get a random integer in the range [0,aRange).
uint32_t
CxMtUint32RangeGet(CxtMt *aMt, uint32_t aRange);

// Get a random real in the range [0,1].
double
CxMtReal1Get(CxtMt *aMt);

// Get a random real in the range [0,1).
double
CxMtReal2Get(CxtMt *aMt);

// Get a random real in the range (0,1).
double
CxMtReal3Get(CxtMt *aMt);

// Get a random real in the range [0,1) with 53-bit resolution.
double
CxMtRes53Get(CxtMt *aMt);
