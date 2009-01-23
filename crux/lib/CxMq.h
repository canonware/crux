#ifndef CxMq_h
#define CxMq_h

#include "Cx.h"

// Pseudo-opaque type.
typedef struct CxsMq CxtMq;

struct CxsMq {
    unsigned msgCount;
    unsigned msgSize;
    unsigned msgsVecCount;
    unsigned msgsBeg;
    unsigned msgsEnd;
    union {
	uint8_t *one;
	uint16_t *two;
	uint32_t *four;
	uint64_t *eight;
	void *x; // Don't care.
    } msgs;

    pthread_mutex_t mtx;
    pthread_cond_t cnd;

    bool getStop;
    bool putStop;
};

bool
CxMqNew(CxtMq *aMq, unsigned aMsgSize, unsigned aMinMsgs);

void
CxMqDelete(CxtMq *aMq);

unsigned
CxMqCount(CxtMq *aMq);

bool
CxMqTryGet(CxtMq *aMq, ...);

bool
CxMqTimedGet(CxtMq *aMq, const struct timespec *aTimeout, ...);

bool
CxMqGet(CxtMq *aMq, ...);

bool
CxMqPut(CxtMq *aMq, ...);

bool
CxMqGetStart(CxtMq *aMq);

bool
CxMqGetStop(CxtMq *aMq);

bool
CxMqPutStart(CxtMq *aMq);

bool
CxMqPutStop(CxtMq *aMq);

#endif // CxMq_h
