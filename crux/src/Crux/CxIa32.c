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

#include "Crux/_cruxmodule.h"

#ifdef CxmCpuIa32

bool CxgIa32UseSse2;

// Check for the cpuid instruction.  If the ID bit (21) in EFLAGS is writable
// then cpuid is supported.
static bool
CxpIa32HasCpuid(void)
{
    bool rVal;
    int before, after;

    asm volatile (
	// Push EFLAGS onto the stack.
	"pushf;"
	// Pop the value of EFLAGS into eax.
	"popl %%eax;"
	// Copy to ecx.
	"movl %%eax, %%ecx;"
	// Toggle the ID bit in eax.
	"xorl $0x00200000, %%eax;"
	// Push eax onto the stack.
	"pushl %%eax;"
	// Pop into EFLAGS.
	"popf;"
	// Push EFLAGS onto the stack.
	"pushf;"
	// Pop the value of EFLAGS into eax.
	"popl %%eax;"
	: "=c" (before), "=a" (after)
	:
	: "cc" // Clobbers condition code register.
	);

    rVal = (before != after) ? true : false;

    return rVal;
}

// Set eax, call the cpuid instruction, and return e[abcd]x.
static void
CxpIa32Cpuid(unsigned aEax, unsigned *rAbcd)
{
    asm volatile (
	"movl %%ebx, %%esi;" // Preserve rAbcd.
	"cpuid;"
	"xchgl %%ebx, %%esi;" // Restore rAbcd.
	: "=a" (rAbcd[0]), "=S" (rAbcd[1]), "=c" (rAbcd[2]), "=d" (rAbcd[3])
	: "0" (aEax)
	);
}

void
CxIa32CpuInit(void)
{
    unsigned abcd[4];

    // If the cpuid instruction is supported, and the SSE2 feature flag is set,
    // enable the use of SSE2.
    if (CxpIa32HasCpuid())
    {
	CxpIa32Cpuid(1, abcd);

	// Mask everything but bit 26 (SSE2 feature flag) of edx.
	if ((abcd[3] & 0x04000000U) != 0)
	{
	    CxgIa32UseSse2 = true;
	}
	else
	{
	    CxgIa32UseSse2 = false;
	}
    }
    else
    {
	// cpuid is missing.  This is an old CPU.
	CxgIa32UseSse2 = false;
    }
}

#endif // CxmCpuIa32
