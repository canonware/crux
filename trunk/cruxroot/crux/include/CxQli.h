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
    uint32_t qlih_first;						\
}

#define CxmQliHeadInitializer(a_head) {UINT_MAX}

#define CxmQliElm	CxmQri

/* List functions. */
#define CxmQliNew(a_head)							\
    do									\
    {									\
	(a_head)->qlih_first = UINT_MAX;				\
    } while (0)

#define CxmQliElm_new(a_arr, a_ind, a_field)				\
    CxmQriNew((a_arr), (a_ind), a_field)

#define CxmQliFirst(a_head) ((a_head)->qlih_first)

#define CxmQliLast(a_head, a_arr, a_field)				\
    ((CxmQliFirst(a_head) != UINT_MAX)					\
      ? CxmQriPrev((a_arr), CxmQliFirst(a_head), a_field) : UINT_MAX)

#define CxmQliNext(a_head, a_arr, a_ind, a_field)				\
    ((CxmQliLast(a_head, a_arr, a_field) != (a_ind))			\
      ? CxmQriNext(a_arr, a_ind, a_field) : UINT_MAX)

#define CxmQliPrev(a_head, a_arr, a_ind, a_field)				\
    ((CxmQliFirst(a_head) != (a_ind)) ? CxmQriPrev(a_arr, a_ind, a_field)	\
				    : UINT_MAX)

#define CxmQliBeforeInsert(a_head, a_arr, a_elmind, a_ind, a_field)	\
    do									\
    {									\
	CxmQriBeforeInsert(a_arr, a_elmind, a_ind, a_field);		\
	if (CxmQliFirst(a_head) == (a_elmind))				\
	{								\
	    CxmQliFirst(a_head) = (a_ind);				\
	}								\
    } while (0)

#define CxmQliAfterInsert(a_arr, a_elmind, a_ind, a_field)		\
    CxmQriAfterInsert(a_arr, a_elmind, a_ind, a_field)

#define CxmQliHead_insert(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (CxmQliFirst(a_head) != UINT_MAX)				\
	{								\
	    CxmQriBeforeInsert(a_arr, CxmQliFirst(a_head), (a_ind),	\
			      a_field);					\
	}								\
	CxmQliFirst(a_head) = (a_ind);					\
    } while (0)

#define CxmQliTailInsert(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (CxmQliFirst(a_head) != UINT_MAX)				\
	{								\
	    CxmQriBeforeInsert(a_arr, CxmQliFirst(a_head), (a_ind),	\
			      a_field);					\
	}								\
	CxmQliFirst(a_head) = CxmQriNext(a_arr, a_ind, a_field);		\
    } while (0)

#define CxmQliRemove(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (CxmQliFirst(a_head) == (a_ind))				\
	{								\
	    CxmQliFirst(a_head) = CxmQriNext(a_arr, CxmQliFirst(a_head),	\
					 a_field);			\
	}								\
	if (CxmQliFirst(a_head) != (a_ind))				\
	{								\
	    CxmQriRemove(a_arr, a_ind, a_field);				\
	}								\
	else								\
	{								\
	    CxmQliFirst(a_head) = UINT_MAX;				\
	}								\
    } while (0)

#define CxmQliHead_remove(a_head, a_arr, a_field)				\
    do									\
    {									\
	uint32_t t = CxmQliFirst(a_head);					\
	CxmQliRemove(a_head, a_arr, t, a_field);				\
    } while (0)

#define CxmQliTailRemove(a_head, a_arr, a_field)				\
    do									\
    {									\
	uint32_t t = CxmQliLast(a_head, a_arr, a_field);			\
	CxmQliRemove((a_head), a_arr, t, a_field);			\
    } while (0)

#define CxmQliForeach(a_var, a_head, a_arr, a_field)			\
    CxmQriForeach((a_var), a_arr, CxmQliFirst(a_head), a_field)

#define CxmQliReverseForeach(a_var, a_head, a_arr, a_field)		\
    CxmQriReverseForeach((a_var), a_arr, CxmQliFirst(a_head), a_field)
