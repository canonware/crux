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

/* Pseudo-opaque type. */
typedef struct CxsXep CxtXep;

typedef uint32_t CxtXepv;

#define CxmXepvNone 0
#define CxmXepvCode 1

typedef enum
{
    CxeXepsTry,
    CxeXepsCatch
} CxtXeps;

struct CxsXep
{
    volatile CxmQr(CxtXep) link;
    volatile CxtXepv value;
    volatile bool isHandled;
    volatile bool isLinked;
    volatile CxtXeps state;
    volatile const char *filename;
    volatile uint32_t lineNum;
    jmp_buf context;
};

void
CxXepInit(void);

void
CxXepShutdown(void);

#define CxmXepBegin()							\
    {									\
	CxtXep _xep

#define CxmXepTry							\
	CxpXepLink(&_xep);						\
	switch (setjmp(_xep.context))					\
	{								\
	    case CxmXepvNone:						\
	    case CxmXepvCode:

#define CxmXepCatch(aValue)						\
		break;							\
	    case (aValue):

#define CxmXepMcatch(aValue)						\
	    case (aValue):

#define CxmXepAcatch							\
		break;							\
	    default:							\
		if (_xep.state != CxeXepsCatch)			\
		{							\
		    break;						\
		}

#define CxmXepEnd()							\
	}								\
	CxpXepUnlink(&_xep);						\
    }

#define CxmXepValue() (_xep.value)

#define CxmXepFilename() (_xep.filename)

#define CxmXepLineNum() (_xep.lineNum)

void
CxXepThrowE(CxtXepv aValue, volatile const char *aFilename,
	    uint32_t aLineNum);

#define CxmXepThrow(aValue) CxXepThrowE((aValue), __FILE__, __LINE__)

#define CxmXepRetry() CxpXepRetry(&_xep)

#define CxmXepHandled() CxpXepHandled(&_xep)

/* Private, but visible here so that the cpp macros above don't cause
 * compilation warnings. */
void
CxpXepRetry(CxtXep *aXep);

void
CxpXepHandled(CxtXep *aXep);

void
CxpXepLink(CxtXep *aXep);

void
CxpXepUnlink(CxtXep *aXep);
