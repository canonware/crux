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

/* List definitions. */
#define CxmQliHead							\
struct									\
{									\
    uint32_t CxQlihFirst;						\
}

#define CxmQliHeadInitializer(aHead) {UINT_MAX}

#define CxmQliElm	CxmQri

/* List functions. */
#define CxmQliNew(aHead)						\
    do									\
    {									\
	(aHead)->CxQlihFirst = UINT_MAX;				\
    } while (0)

#define CxmQliElmNew(aArr, aInd, aField)				\
    CxmQriNew((aArr), (aInd), aField)

#define CxmQliFirst(aHead) ((aHead)->CxQlihFirst)

#define CxmQliLast(aHead, aArr, aField)					\
    ((CxmQliFirst(aHead) != UINT_MAX)					\
      ? CxmQriPrev((aArr), CxmQliFirst(aHead), aField) : UINT_MAX)

#define CxmQliNext(aHead, aArr, aInd, aField)				\
    ((CxmQliLast(aHead, aArr, aField) != (aInd))			\
      ? CxmQriNext(aArr, aInd, aField) : UINT_MAX)

#define CxmQliPrev(aHead, aArr, aInd, aField)				\
    ((CxmQliFirst(aHead) != (aInd))					\
			     ? CxmQriPrev(aArr, aInd, aField)		\
			     : UINT_MAX)

#define CxmQliBeforeInsert(aHead, aArr, aElmind, aInd, aField)		\
    do									\
    {									\
	CxmQriBeforeInsert(aArr, aElmind, aInd, aField);		\
	if (CxmQliFirst(aHead) == (aElmind))				\
	{								\
	    CxmQliFirst(aHead) = (aInd);				\
	}								\
    } while (0)

#define CxmQliAfterInsert(aArr, aElmind, aInd, aField)			\
    CxmQriAfterInsert(aArr, aElmind, aInd, aField)

#define CxmQliHeadInsert(aHead, aArr, aInd, aField)			\
    do									\
    {									\
	if (CxmQliFirst(aHead) != UINT_MAX)				\
	{								\
	    CxmQriBeforeInsert(aArr, CxmQliFirst(aHead), (aInd),	\
			      aField);					\
	}								\
	CxmQliFirst(aHead) = (aInd);					\
    } while (0)

#define CxmQliTailInsert(aHead, aArr, aInd, aField)			\
    do									\
    {									\
	if (CxmQliFirst(aHead) != UINT_MAX)				\
	{								\
	    CxmQriBeforeInsert(aArr, CxmQliFirst(aHead), (aInd),	\
			      aField);					\
	}								\
	CxmQliFirst(aHead) = CxmQriNext(aArr, aInd, aField);		\
    } while (0)

#define CxmQliRemove(aHead, aArr, aInd, aField)				\
    do									\
    {									\
	if (CxmQliFirst(aHead) == (aInd))				\
	{								\
	    CxmQliFirst(aHead) = CxmQriNext(aArr, CxmQliFirst(aHead),	\
					 aField);			\
	}								\
	if (CxmQliFirst(aHead) != (aInd))				\
	{								\
	    CxmQriRemove(aArr, aInd, aField);				\
	}								\
	else								\
	{								\
	    CxmQliFirst(aHead) = UINT_MAX;				\
	}								\
    } while (0)

#define CxmQliHeadRemove(aHead, aArr, aField)				\
    do									\
    {									\
	uint32_t t = CxmQliFirst(aHead);				\
	CxmQliRemove(aHead, aArr, t, aField);				\
    } while (0)

#define CxmQliTailRemove(aHead, aArr, aField)				\
    do									\
    {									\
	uint32_t t = CxmQliLast(aHead, aArr, aField);			\
	CxmQliRemove((aHead), aArr, t, aField);				\
    } while (0)

#define CxmQliForeach(aVar, aHead, aArr, aField)			\
    CxmQriForeach((aVar), aArr, CxmQliFirst(aHead), aField)

#define CxmQliReverseForeach(aVar, aHead, aArr, aField)			\
    CxmQriReverseForeach((aVar), aArr, CxmQliFirst(aHead), aField)
