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
#define CxmQr(a_type)							\
struct									\
{									\
    a_type *CxmQre_next;							\
    a_type *CxmQre_prev;							\
}

/* Ring functions. */
#define CxmQrNew(a_CxmQr, a_field)						\
    do									\
    {									\
	(a_CxmQr)->a_field.CxmQre_next = (a_CxmQr);				\
	(a_CxmQr)->a_field.CxmQre_prev = (a_CxmQr);				\
    } while (0)

#define CxmQrNext(a_CxmQr, a_field) ((a_CxmQr)->a_field.CxmQre_next)

#define CxmQrPrev(a_CxmQr, a_field) ((a_CxmQr)->a_field.CxmQre_prev)

#define CxmQrBeforeInsert(a_CxmQrelm, a_CxmQr, a_field)			\
    do									\
    {									\
	(a_CxmQr)->a_field.CxmQre_prev = (a_CxmQrelm)->a_field.CxmQre_prev;		\
	(a_CxmQr)->a_field.CxmQre_next = (a_CxmQrelm);				\
	(a_CxmQr)->a_field.CxmQre_prev->a_field.CxmQre_next = (a_CxmQr);		\
	(a_CxmQrelm)->a_field.CxmQre_prev = (a_CxmQr);				\
    } while (0)

#define CxmQrAfterInsert(a_CxmQrelm, a_CxmQr, a_field)				\
    do									\
    {									\
	(a_CxmQr)->a_field.CxmQre_next = (a_CxmQrelm)->a_field.CxmQre_next;		\
	(a_CxmQr)->a_field.CxmQre_prev = (a_CxmQrelm);				\
	(a_CxmQr)->a_field.CxmQre_next->a_field.CxmQre_prev = (a_CxmQr);		\
	(a_CxmQrelm)->a_field.CxmQre_next = (a_CxmQr);				\
    } while (0)

#define CxmQrMeld(a_CxmQr_a, a_CxmQr_b, a_type, a_field)			\
    do									\
    {									\
	a_type *t;							\
	(a_CxmQr_a)->a_field.CxmQre_prev->a_field.CxmQre_next = (a_CxmQr_b);	\
	(a_CxmQr_b)->a_field.CxmQre_prev->a_field.CxmQre_next = (a_CxmQr_a);	\
	t = (a_CxmQr_a)->a_field.CxmQre_prev;					\
	(a_CxmQr_a)->a_field.CxmQre_prev = (a_CxmQr_b)->a_field.CxmQre_prev;	\
	(a_CxmQr_b)->a_field.CxmQre_prev = t;					\
    } while (0)

/* CxmQrMeld() and CxmQrSplit() are functionally equivalent, so there's no need to
 * have two copies of the code. */
#define CxmQrSplit(a_CxmQr_a, a_CxmQr_b, a_type, a_field)			\
    CxmQrMeld((a_CxmQr_a), (a_CxmQr_b), a_type, a_field)

#define CxmQrRemove(a_CxmQr, a_field)					\
    do									\
    {									\
	(a_CxmQr)->a_field.CxmQre_prev->a_field.CxmQre_next			\
	  = (a_CxmQr)->a_field.CxmQre_next;					\
	(a_CxmQr)->a_field.CxmQre_next->a_field.CxmQre_prev			\
	  = (a_CxmQr)->a_field.CxmQre_prev;					\
	(a_CxmQr)->a_field.CxmQre_next = (a_CxmQr);				\
	(a_CxmQr)->a_field.CxmQre_prev = (a_CxmQr);				\
    } while (0)

#define CxmQrForeach(var, a_CxmQr, a_field)					\
    for ((var) = (a_CxmQr);						\
	 (var) != NULL;							\
	 (var) = (((var)->a_field.CxmQre_next != (a_CxmQr))			\
	 ? (var)->a_field.CxmQre_next : NULL))

#define CxmQrReverseForeach(var, a_CxmQr, a_field)				\
    for ((var) = ((a_CxmQr) != NULL) ? CxmQrPrev(a_CxmQr, a_field) : NULL;	\
	 (var) != NULL;							\
	 (var) = (((var) != (a_CxmQr))					\
	 ? (var)->a_field.CxmQre_prev : NULL))
