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

#define CW_XEPV_NONE 0
#define CW_XEPV_CODE 1

typedef enum
{
    CW_XEPS_TRY,
    CW_XEPS_CATCH
} cw_xeps_t;

struct CxsXep
{
    volatile CxmQr(CxtXep) link;
    volatile CxtXepv value;
    volatile bool is_handled;
    volatile bool is_linked;
    volatile cw_xeps_t state;
    volatile const char *filename;
    volatile uint32_t line_num;
    jmp_buf context;
};

void
CxXepInit(void);

void
CxXepShutdown(void);

#define CxmXepBegin()							\
    {									\
	CxtXep _xep

#define xep_try								\
	CxpXepLink(&_xep);						\
	switch (setjmp(_xep.context))					\
	{								\
	    case CW_XEPV_NONE:						\
	    case CW_XEPV_CODE:

#define CxmXepCatch(a_value)						\
		break;							\
	    case (a_value):

#define CxmXepMcatch(a_value)						\
	    case (a_value):

#define CxmXepAcatch							\
		break;							\
	    default:							\
		if (_xep.state != CW_XEPS_CATCH)			\
		{							\
		    break;						\
		}

#define CxmXepEnd()							\
	}								\
	CxpXepUnlink(&_xep);						\
    }

#define CxmXepValue() (_xep.value)

#define CxmXepFilename() (_xep.filename)

#define CxmXepLineNum() (_xep.line_num)

void
CxXepThrowE(CxtXepv a_value, volatile const char *a_filename,
	    uint32_t a_line_num);

#define CxmXepThrow(a_value) CxXepThrowE((a_value), __FILE__, __LINE__)

#define CxmXepRetry() CxpXepRetry(&_xep)

#define CxmXepHandled() CxpXepHandled(&_xep)

/* Private, but visible here so that the cpp macros above don't cause
 * compilation warnings. */
void
CxpXepRetry(CxtXep *a_xep);

void
CxpXepHandled(CxtXep *a_xep);

void
CxpXepLink(CxtXep *a_xep);

void
CxpXepUnlink(CxtXep *a_xep);
