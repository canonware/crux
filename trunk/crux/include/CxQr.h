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
#define CxmQr(aType)							\
struct									\
{									\
    aType *CxQreNext;							\
    aType *CxQrePrev;							\
}

/* Ring functions. */
#define CxmQrNew(aQr, aField)						\
    do									\
    {									\
	(aQr)->aField.CxQreNext = (aQr);				\
	(aQr)->aField.CxQrePrev = (aQr);				\
    } while (0)

#define CxmQrNext(aQr, aField) ((aQr)->aField.CxQreNext)

#define CxmQrPrev(aQr, aField) ((aQr)->aField.CxQrePrev)

#define CxmQrBeforeInsert(aQrElm, aQr, aField)				\
    do									\
    {									\
	(aQr)->aField.CxQrePrev = (aQrElm)->aField.CxQrePrev;		\
	(aQr)->aField.CxQreNext = (aQrElm);				\
	(aQr)->aField.CxQrePrev->aField.CxQreNext = (aQr);		\
	(aQrElm)->aField.CxQrePrev = (aQr);				\
    } while (0)

#define CxmQrAfterInsert(aQrElm, aQr, aField)				\
    do									\
    {									\
	(aQr)->aField.CxQreNext = (aQrElm)->aField.CxQreNext;		\
	(aQr)->aField.CxQrePrev = (aQrElm);				\
	(aQr)->aField.CxQreNext->aField.CxQrePrev = (aQr);		\
	(aQrElm)->aField.CxQreNext = (aQr);				\
    } while (0)

#define CxmQrMeld(aQrA, aQrB, aType, aField)				\
    do									\
    {									\
	aType *t;							\
	(aQrA)->aField.CxQrePrev->aField.CxQreNext = (aQrB);		\
	(aQrB)->aField.CxQrePrev->aField.CxQreNext = (aQrA);		\
	t = (aQrA)->aField.CxQrePrev;					\
	(aQrA)->aField.CxQrePrev = (aQrB)->aField.CxQrePrev;		\
	(aQrB)->aField.CxQrePrev = t;					\
    } while (0)

/* CxmQrMeld() and CxmQrSplit() are functionally equivalent, so there's no need
 * to have two copies of the code. */
#define CxmQrSplit(aQrA, aQrB, aType, aField)				\
    CxmQrMeld((aQrA), (aQrB), aType, aField)

#define CxmQrRemove(aQr, aField)					\
    do									\
    {									\
	(aQr)->aField.CxQrePrev->aField.CxQreNext			\
	  = (aQr)->aField.CxQreNext;					\
	(aQr)->aField.CxQreNext->aField.CxQrePrev			\
	  = (aQr)->aField.CxQrePrev;					\
	(aQr)->aField.CxQreNext = (aQr);				\
	(aQr)->aField.CxQrePrev = (aQr);				\
    } while (0)

#define CxmQrForeach(aVar, aQr, aField)					\
    for ((aVar) = (aQr);						\
	 (aVar) != NULL;						\
	 (aVar) = (((aVar)->aField.CxQreNext != (aQr))			\
	 ? (aVar)->aField.CxQreNext : NULL))

#define CxmQrReverseForeach(aVar, aQr, aField)				\
    for ((aVar) = ((aQr) != NULL) ? CxmQrPrev(aQr, aField) : NULL;	\
	 (aVar) != NULL;						\
	 (aVar) = (((aVar) != (aQr))					\
	 ? (aVar)->aField.CxQrePrev : NULL))
