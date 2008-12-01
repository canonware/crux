#include "CxAmd64.h"

#ifdef CxmCpuAmd64

bool CxgAmd64UseSse2;

// Set eax, call the cpuid instruction, and return e[abcd]x.
static void
CxpAmd64Cpuid(unsigned aEax, unsigned *rAbcd)
{
    __asm__ volatile (
	"movl %%ebx, %%esi;" // Preserve rAbcd.
	"cpuid;"
	"xchgl %%ebx, %%esi;" // Restore rAbcd.
	: "=a" (rAbcd[0]), "=S" (rAbcd[1]), "=c" (rAbcd[2]), "=d" (rAbcd[3])
	: "0" (aEax)
	);
}

void
CxAmd64CpuInit(void)
{
    unsigned abcd[4];

    // If the cpuid instruction is supported, and the SSE2 feature flag is set,
    // enable the use of SSE2.
    CxpAmd64Cpuid(1, abcd);

    // Mask everything but bit 26 (SSE2 feature flag) of edx.
    if ((abcd[3] & 0x04000000) != 0)
    {
	CxgAmd64UseSse2 = true;
    }
    else
    {
	CxgAmd64UseSse2 = false;
    }
}

#endif // CxmCpuAmd64
