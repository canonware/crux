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

/* Ring definitions. */
#define CxmQri								\
struct									\
{									\
    uint32_t CxQrieNext;						\
    uint32_t CxQriePrev;						\
}

/* Ring functions. */
#define CxmQriNew(aArr, aInd, aField)					\
    do									\
    {									\
	(aArr)[aInd].aField.CxQrieNext = (aInd);			\
	(aArr)[aInd].aField.CxQriePrev = (aInd);			\
    } while (0)

#define CxmQriNext(aArr, aInd, aField)					\
    ((aArr)[aInd].aField.CxQrieNext)

#define CxmQriPrev(aArr, aInd, aField)					\
    ((aArr)[aInd].aField.CxQriePrev)

#define CxmQriBeforeInsert(aArr, aElmind, aInd, aField)			\
    do									\
    {									\
	(aArr)[aInd].aField.CxQriePrev					\
	  = (aArr)[aElmind].aField.CxQriePrev;				\
	(aArr)[aInd].aField.CxQrieNext = (aElmind);			\
	(aArr)[(aArr)[aInd].aField.CxQriePrev].aField.CxQrieNext	\
	  = (aInd);							\
	(aArr)[aElmind].aField.CxQriePrev = (aInd);			\
    } while (0)

#define CxmQriAfterInsert(aArr, aElmind, aInd, aField)			\
    do									\
    {									\
	(aArr)[aInd].aField.CxQrieNext					\
	  = (aArr)[aElmind].aField.CxQrieNext;				\
	(aArr)[aInd].aField.CxQriePrev = (aElmind);			\
	(aArr)[(aArr)[aInd].aField.CxQrieNext].aField.CxQriePrev	\
	  = (aInd);							\
	(aArr)[aElmind].aField.CxQrieNext = (aInd);			\
    } while (0)

#define CxmQriMeld(aArr, aIndA, aIndB, aField)				\
    do									\
    {									\
	uint32_t t;							\
	(aArr)[(aArr)[aIndA].aField.CxQriePrev].aField.CxQrieNext	\
	  = (aIndB);							\
	(aArr)[(aArr)[aIndB].aField.CxQriePrev].aField.CxQrieNext	\
	  = (aIndA);							\
	t = (aArr)[aIndA].aField.CxQriePrev;				\
	(aArr)[aIndA].aField.CxQriePrev					\
	  = (aArr)[aIndB].aField.CxQriePrev;				\
	(aArr)[aIndB].aField.CxQriePrev = t;				\
    } while (0)

/* CxmQriMeld() and CxmQriSplit() are functionally equivalent, so there's no
 * need to have two copies of the code. */
#define CxmQriSplit(aArr, aIndA, aIndB, aField)				\
    CxmQriMeld((aArr), (aIndA), (aIndB), aField)

#define CxmQriRemove(aArr, aInd, aField)				\
    do									\
    {									\
	(aArr)[(aArr)[aInd].aField.CxQriePrev].aField.CxQrieNext	\
	  = (aArr)[aInd].aField.CxQrieNext;				\
	(aArr)[(aArr)[aInd].aField.CxQrieNext].aField.CxQriePrev	\
	  = (aArr)[aInd].aField.CxQriePrev;				\
	(aArr)[aInd].aField.CxQrieNext = (aInd);			\
	(aArr)[aInd].aField.CxQriePrev = (aInd);			\
    } while (0)

#define CxmQriForeach(var, aArr, aInd, aField)				\
    for ((var) = (aInd);						\
	 (var) != UINT_MAX;						\
	 (var) = (((aArr)[var].aField.CxQrieNext != (aInd))		\
	 ? (aArr)[var].aField.CxQrieNext : UINT_MAX))

#define CxmQriOthersForeach(var, aArr, aInd, aField)			\
    for ((var) = ((aInd) != UINT_MAX)					\
		  ? (aArr)[aInd].aField.CxQrieNext : UINT_MAX;		\
	 (var) != (aInd);						\
	 (var) = (aArr)[var].aField.CxQrieNext)

#define CxmQriReverseForeach(var, aArr, aInd, aField)			\
    for ((var) = ((aInd) != UINT_MAX)					\
	       ? CxmQriPrev(aArr, aInd, aField) : UINT_MAX;		\
	 (var) != UINT_MAX;						\
	 (var) = (((var) != (aInd))					\
	 ? (aArr)[var].aField.CxQriePrev : UINT_MAX))

#define CxmQriOthersReverseForeach(var, aArr, aInd, aField)		\
    for ((var) = ((aInd) != UINT_MAX)					\
		  ? (aArr)[aInd].aField.CxQriePrev : UINT_MAX;		\
	 (var) != (aInd);						\
	 (var) = (aArr)[var].aField.CxQriePrev)
