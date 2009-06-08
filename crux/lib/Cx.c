#define Cx_c
#include "Cx.h"
#ifdef CxmHaveMallopt
#include <malloc.h>
#endif

unsigned CxNcpus = 1;

static pthread_once_t CxpThreadedOnce = PTHREAD_ONCE_INIT;

void
CxInit(void) {
#ifdef CxmCpuInit
    CxmCpuInit();
#endif
#ifdef CxmHaveMallopt
    // Disable ptmalloc's sliding mmap() threshold, which can cause unbounded
    // fragmentation.
    mallopt(M_MMAP_THRESHOLD, 128*1024);
#endif
}

#if (defined(CxmOsBSD) || defined(CxmOsLinux) || defined(CxmOsSolaris))
CxmpInline unsigned
CxpNcpus(void)
{
    return sysconf(_SC_NPROCESSORS_ONLN);
}
#elif (defined(CxmOsDarwin))
#include <mach/mach_init.h>
#include <mach/mach_host.h>

CxmpInline unsigned
CxpNcpus(void)
{
    kern_return_t error;
    natural_t n;
    processor_info_array_t pinfo;
    mach_msg_type_number_t pinfocnt;

    error = host_processor_info(mach_host_self(), PROCESSOR_BASIC_INFO, &n,
      &pinfo, &pinfocnt);
    if (error != KERN_SUCCESS)
	return (1); /* Error. */
    else
	return (n);
}
#else
CxmpInline unsigned
CxpNcpus(void)
{
    // We lack a way to determine the number of CPUs on this platform, so
    // assume 1 CPU.
    return (1);
}
#endif

static void
CxpThreaded(void) {
    CxNcpus = CxpNcpus();
}

void
CxThreaded(void) {
    pthread_once(&CxpThreadedOnce, CxpThreaded);
}
