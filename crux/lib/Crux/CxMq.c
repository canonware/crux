#include "CxMq.h"

#include <stdarg.h>
#include <sys/time.h>
#include <time.h>

bool
CxMqNew(CxtMq *aMq, unsigned aMsgSize, unsigned aMinMsgs) {
    CxmCheckPtr(aMq);

    aMq->msgCount = 0;

    switch (aMsgSize) {
	case 1: {
	    aMq->msgSize = 1;
	    break;
	} case 2: {
	    aMq->msgSize = 2;
	    break;
	} case 4: {
	    aMq->msgSize = 4;
	    break;
	} case 8: {
	    aMq->msgSize = 8;
	    break;
	} default: {
	    CxmNotReached();
	}
    }

    aMq->msgsVecCount = aMinMsgs;
    aMq->msgsBeg = 0;
    aMq->msgsEnd = 0;

    aMq->msgs.x = malloc(aMq->msgsVecCount * aMq->msgSize);
    if (aMq->msgs.x == NULL) {
	return true;
    }

    if (pthread_mutex_init(&aMq->mtx, NULL)) {
	return true;
    }
    if (pthread_cond_init(&aMq->cnd, NULL)) {
	return true;
    }

    aMq->getStop = false;
    aMq->putStop = false;

    return false;
}

void
CxMqDelete(CxtMq *aMq) {
    CxmCheckPtr(aMq);

    pthread_mutex_destroy(&aMq->mtx);
    pthread_cond_destroy(&aMq->cnd);

    free(aMq->msgs.x);
}

unsigned
CxMqCount(CxtMq *aMq) {
    unsigned ret;

    CxmCheckPtr(aMq);

    pthread_mutex_lock(&aMq->mtx);
    ret = aMq->msgCount;
    pthread_mutex_unlock(&aMq->mtx);

    return ret;
}

bool
CxMqTryGet(CxtMq *aMq, ...) {
    bool ret;
    union {
	uint8_t *one;
	uint16_t *two;
	uint32_t *four;
	uint64_t *eight;
	void *x; // Don't care.
    } rMsg;
    va_list ap;

    CxmCheckPtr(aMq);

    va_start(ap, aMq);
    rMsg.x = (void *) va_arg(ap, void *);
    va_end(ap);

    pthread_mutex_lock(&aMq->mtx);

    if (aMq->getStop) {
	ret = true;
	goto RETURN;
    }
    if (aMq->msgCount > 0) {
	switch (aMq->msgSize) {
	    case 1: {
		*rMsg.one = aMq->msgs.one[aMq->msgsBeg];
		break;
	    } case 2: {
		*rMsg.two = aMq->msgs.two[aMq->msgsBeg];
		break;
	    } case 4: {
		*rMsg.four = aMq->msgs.four[aMq->msgsBeg];
		break;
	    } case 8: {
		*rMsg.eight = aMq->msgs.eight[aMq->msgsBeg];
		break;
	    } default: {
		CxmNotReached();
	    }
	}
	aMq->msgCount--;
	aMq->msgsBeg = (aMq->msgsBeg + 1) % aMq->msgsVecCount;
    } else {
	ret = true;
	goto RETURN;
    }

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

CxmpInline bool
CxpMqTimedWait(pthread_cond_t *aCnd, pthread_mutex_t *aMtx,
  const struct timespec *aTimeout) {
    bool ret;
    int error;
    struct timeval now;
    struct timespec timeout;
    struct timezone tz;

    CxmCheckPtr(aCnd);
    CxmCheckPtr(aMtx);
    CxmCheckPtr(aTimeout);

    // Set timeout.
    memset(&tz, 0, sizeof(struct timezone));
    gettimeofday(&now, &tz);
    timeout.tv_nsec = now.tv_usec * 1000 + aTimeout->tv_nsec;
    // Carry if nanoseconds overflowed.
    timeout.tv_sec = (now.tv_sec + aTimeout->tv_sec
		      + (timeout.tv_nsec / 1000000000));
    // Chop off the number of nanoseconds to be less than one second.
    timeout.tv_nsec %= 1000000000;

    error = pthread_cond_timedwait(aCnd, aMtx, &timeout);
    if (error == 0) {
	ret = false;
    } else {
	ret = true;
    }

    return ret;
}

bool
CxMqTimedGet(CxtMq *aMq, const struct timespec *aTimeout, ...) {
    bool ret, timedOut;
    union {
	uint8_t *one;
	uint16_t *two;
	uint32_t *four;
	uint64_t *eight;
	void *x; // Don't care.
    } rMsg;
    va_list ap;

    CxmCheckPtr(aMq);
    CxmCheckPtr(aTimeout);

    va_start(ap, aTimeout);
    rMsg.x = (void *) va_arg(ap, void *);
    va_end(ap);

    pthread_mutex_lock(&aMq->mtx);

    if (aMq->getStop) {
	ret = true;
	goto RETURN;
    }
    // A spurious wakeup will cause the timeout interval to start over.  This
    // isn't a big deal as long as spurious wakeups don't occur continuously,
    // since the timeout period is merely a lower bound on how long to wait.
    for (timedOut = false; aMq->msgCount == 0 && timedOut == false;) {
	timedOut = CxpMqTimedWait(&aMq->cnd, &aMq->mtx, aTimeout);
	if (aMq->getStop) {
	    ret = true;
	    goto RETURN;
	}
    }
    if (aMq->msgCount > 0) {
	switch (aMq->msgSize) {
	    case 1: {
		*rMsg.one = aMq->msgs.one[aMq->msgsBeg];
		break;
	    } case 2: {
		*rMsg.two = aMq->msgs.two[aMq->msgsBeg];
		break;
	    } case 4: {
		*rMsg.four = aMq->msgs.four[aMq->msgsBeg];
		break;
	    } case 8: {
		*rMsg.eight = aMq->msgs.eight[aMq->msgsBeg];
		break;
	    } default: {
		CxmNotReached();
	    }
	}
	aMq->msgCount--;
	aMq->msgsBeg = (aMq->msgsBeg + 1) % aMq->msgsVecCount;
    } else {
	ret = true;
	goto RETURN;
    }

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqGet(CxtMq *aMq, ...) {
    bool ret;
    union {
	uint8_t *one;
	uint16_t *two;
	uint32_t *four;
	uint64_t *eight;
	void *x; // Don't care.
    } rMsg;
    va_list ap;

    CxmCheckPtr(aMq);

    va_start(ap, aMq);
    rMsg.x = (void *) va_arg(ap, void *);
    va_end(ap);

    pthread_mutex_lock(&aMq->mtx);

    if (aMq->getStop) {
	ret = true;
	goto RETURN;
    }
    while (aMq->msgCount == 0) {
	pthread_cond_wait(&aMq->cnd, &aMq->mtx);
	if (aMq->getStop) {
	    ret = true;
	    goto RETURN;
	}
    }

    switch (aMq->msgSize) {
	case 1: {
	    *rMsg.one = aMq->msgs.one[aMq->msgsBeg];
	    break;
	} case 2: {
	    *rMsg.two = aMq->msgs.two[aMq->msgsBeg];
	    break;
	} case 4: {
	    *rMsg.four = aMq->msgs.four[aMq->msgsBeg];
	    break;
	} case 8: {
	    *rMsg.eight = aMq->msgs.eight[aMq->msgsBeg];
	    break;
	} default: {
	    CxmNotReached();
	}
    }
    aMq->msgCount--;
    aMq->msgsBeg = (aMq->msgsBeg + 1) % aMq->msgsVecCount;

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqPut(CxtMq *aMq, ...) {
    bool ret;
    union {
	uint32_t four;
	uint64_t eight;
    } aMsg;
    va_list ap;

    CxmCheckPtr(aMq);

    va_start(ap, aMq);
    switch (aMq->msgSize) {
	case 1: // Compiler-promoted to 32 bits.
	case 2: // Compiler-promoted to 32 bits.
	case 4: {
	    aMsg.four = va_arg(ap, uint32_t);
	    break;
	} case 8: {
	    aMsg.eight = va_arg(ap, uint64_t);
	    break;
	} default: {
	    CxmNotReached();
	}
    }
    va_end(ap);

    pthread_mutex_lock(&aMq->mtx);

    if (aMq->msgCount == 0) {
	pthread_cond_broadcast(&aMq->cnd);
    }
    if (aMq->putStop) {
	ret = true;
	goto RETURN;
    } else {
	if (aMq->msgCount >= aMq->msgsVecCount) {
	    union {
		uint8_t *one;
		uint16_t *two;
		uint32_t *four;
		uint64_t *eight;
		void *x; // Don't care.
	    } tMsgs;
	    unsigned i, offset;

	    // Array overflow.  Double the array and copy the messages.
	    tMsgs.x = malloc(aMq->msgsVecCount * 2 * aMq->msgSize);
	    if (tMsgs.x == NULL) {
		ret = true;
		goto RETURN;
	    }

	    switch (aMq->msgSize) {
		case 1: {
		    for (i = 0, offset = aMq->msgsBeg;
			 i < aMq->msgCount;
			 i++, offset = (offset + 1) % aMq->msgsVecCount) {
			tMsgs.one[i] = aMq->msgs.one[offset];
		    }
		    break;
		} case 2: {
		    for (i = 0, offset = aMq->msgsBeg;
			 i < aMq->msgCount;
			 i++, offset = (offset + 1) % aMq->msgsVecCount) {
			tMsgs.two[i] = aMq->msgs.two[offset];
		    }
		    break;
		} case 4: {
		    for (i = 0, offset = aMq->msgsBeg;
			 i < aMq->msgCount;
			 i++, offset = (offset + 1) % aMq->msgsVecCount) {
			tMsgs.four[i] = aMq->msgs.four[offset];
		    }
		    break;
		} case 8: {
		    for (i = 0, offset = aMq->msgsBeg;
			 i < aMq->msgCount;
			 i++, offset = (offset + 1) % aMq->msgsVecCount) {
			tMsgs.eight[i] = aMq->msgs.eight[offset];
		    }
		    break;
		} default: {
		    CxmNotReached();
		}
	    }
	    free(aMq->msgs.x);
	    aMq->msgs.x = tMsgs.x;
	    aMq->msgsBeg = 0;
	    aMq->msgsEnd = aMq->msgCount;
	    aMq->msgsVecCount *= 2;
	}

	switch (aMq->msgSize) {
	    case 1: {
		// The compiler promotes 8 bit args to 32 bits.
		aMq->msgs.one[aMq->msgsEnd] = aMsg.four;
		break;
	    } case 2: {
		// The compiler promotes 16 bit args to 32 bits.
		aMq->msgs.two[aMq->msgsEnd] = aMsg.four;
		break;
	    } case 4: {
		aMq->msgs.four[aMq->msgsEnd] = aMsg.four;
		break;
	    } case 8: {
		aMq->msgs.eight[aMq->msgsEnd] = aMsg.eight;
		break;
	    } default: {
		CxmNotReached();
	    }
	}
	aMq->msgCount++;
	aMq->msgsEnd = (aMq->msgsEnd + 1) % aMq->msgsVecCount;
    }

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqGetStart(CxtMq *aMq) {
    bool ret;

    CxmCheckPtr(aMq);
    pthread_mutex_lock(&aMq->mtx);

    if (aMq->getStop == false) {
	ret = true;
	goto RETURN;
    }
    aMq->getStop = false;

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqGetStop(CxtMq *aMq) {
    bool ret;

    CxmCheckPtr(aMq);
    pthread_mutex_lock(&aMq->mtx);

    if (aMq->getStop) {
	ret = true;
	goto RETURN;
    }
    pthread_cond_broadcast(&aMq->cnd);
    aMq->getStop = true;

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqPutStart(CxtMq *aMq) {
    bool ret;

    CxmCheckPtr(aMq);
    pthread_mutex_lock(&aMq->mtx);

    if (aMq->putStop == false) {
	ret = true;
	goto RETURN;
    }
    aMq->putStop = false;

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}

bool
CxMqPutStop(CxtMq *aMq) {
    bool ret;

    CxmCheckPtr(aMq);
    pthread_mutex_lock(&aMq->mtx);

    if (aMq->putStop) {
	ret = true;
	goto RETURN;
    }
    aMq->putStop = true;

    ret = false;
    RETURN:
    pthread_mutex_unlock(&aMq->mtx);
    return ret;
}
