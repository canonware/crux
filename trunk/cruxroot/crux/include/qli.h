/* -*- mode: c ; c-file-style: "canonware-c-style" -*-
 ******************************************************************************
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
#define qli_head							\
struct									\
{									\
    uint32_t qlih_first;						\
}

#define qli_head_initializer(a_head) {UINT_MAX}

#define qli_elm	qri

/* List functions. */
#define qli_new(a_head)							\
    do									\
    {									\
	(a_head)->qlih_first = UINT_MAX;				\
    } while (0)

#define qli_elm_new(a_arr, a_ind, a_field)				\
    qri_new((a_arr), (a_ind), a_field)

#define qli_first(a_head) ((a_head)->qlih_first)

#define qli_last(a_head, a_arr, a_field)				\
    ((qli_first(a_head) != UINT_MAX)					\
      ? qri_prev((a_arr), qli_first(a_head), a_field) : UINT_MAX)

#define qli_next(a_head, a_arr, a_ind, a_field)				\
    ((qli_last(a_head, a_arr, a_field) != (a_ind))			\
      ? qri_next(a_arr, a_ind, a_field) : UINT_MAX)

#define qli_prev(a_head, a_arr, a_ind, a_field)				\
    ((qli_first(a_head) != (a_ind)) ? qri_prev(a_arr, a_ind, a_field)	\
				    : UINT_MAX)

#define qli_before_insert(a_head, a_arr, a_elmind, a_ind, a_field)	\
    do									\
    {									\
	qri_before_insert(a_arr, a_elmind, a_ind, a_field);		\
	if (qli_first(a_head) == (a_elmind))				\
	{								\
	    qli_first(a_head) = (a_ind);				\
	}								\
    } while (0)

#define qli_after_insert(a_arr, a_elmind, a_ind, a_field)		\
    qri_after_insert(a_arr, a_elmind, a_ind, a_field)

#define qli_head_insert(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (qli_first(a_head) != UINT_MAX)				\
	{								\
	    qri_before_insert(a_arr, qli_first(a_head), (a_ind),	\
			      a_field);					\
	}								\
	qli_first(a_head) = (a_ind);					\
    } while (0)

#define qli_tail_insert(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (qli_first(a_head) != UINT_MAX)				\
	{								\
	    qri_before_insert(a_arr, qli_first(a_head), (a_ind),	\
			      a_field);					\
	}								\
	qli_first(a_head) = qri_next(a_arr, a_ind, a_field);		\
    } while (0)

#define qli_remove(a_head, a_arr, a_ind, a_field)			\
    do									\
    {									\
	if (qli_first(a_head) == (a_ind))				\
	{								\
	    qli_first(a_head) = qri_next(a_arr, qli_first(a_head),	\
					 a_field);			\
	}								\
	if (qli_first(a_head) != (a_ind))				\
	{								\
	    qri_remove(a_arr, a_ind, a_field);				\
	}								\
	else								\
	{								\
	    qli_first(a_head) = UINT_MAX;				\
	}								\
    } while (0)

#define qli_head_remove(a_head, a_arr, a_field)				\
    do									\
    {									\
	uint32_t t = qli_first(a_head);					\
	qli_remove(a_head, a_arr, t, a_field);				\
    } while (0)

#define qli_tail_remove(a_head, a_arr, a_field)				\
    do									\
    {									\
	uint32_t t = qli_last(a_head, a_arr, a_field);			\
	qli_remove((a_head), a_arr, t, a_field);			\
    } while (0)

#define qli_foreach(a_var, a_head, a_arr, a_field)			\
    qri_foreach((a_var), a_arr, qli_first(a_head), a_field)

#define qli_reverse_foreach(a_var, a_head, a_arr, a_field)		\
    qri_reverse_foreach((a_var), a_arr, qli_first(a_head), a_field)
