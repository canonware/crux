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
#define qri								\
struct									\
{									\
    uint32_t qrie_next;							\
    uint32_t qrie_prev;							\
}

/* Ring functions. */
#define qri_new(a_arr, a_ind, a_field)					\
    do									\
    {									\
	(a_arr)[a_ind].a_field.qrie_next = (a_ind);			\
	(a_arr)[a_ind].a_field.qrie_prev = (a_ind);			\
    } while (0)

#define qri_next(a_arr, a_ind, a_field)					\
    ((a_arr)[a_ind].a_field.qrie_next)

#define qri_prev(a_arr, a_ind, a_field)					\
    ((a_arr)[a_ind].a_field.qrie_prev)

#define qri_before_insert(a_arr, a_elmind, a_ind, a_field)		\
    do									\
    {									\
	(a_arr)[a_ind].a_field.qrie_prev				\
	  = (a_arr)[a_elmind].a_field.qrie_prev;			\
	(a_arr)[a_ind].a_field.qrie_next = (a_elmind);			\
	(a_arr)[(a_arr)[a_ind].a_field.qrie_prev].a_field.qrie_next	\
	  = (a_ind);							\
	(a_arr)[a_elmind].a_field.qrie_prev = (a_ind);			\
    } while (0)

#define qri_after_insert(a_arr, a_elmind, a_ind, a_field)		\
    do									\
    {									\
	(a_arr)[a_ind].a_field.qrie_next				\
	  = (a_arr)[a_elmind].a_field.qrie_next;			\
	(a_arr)[a_ind].a_field.qrie_prev = (a_elmind);			\
	(a_arr)[(a_arr)[a_ind].a_field.qrie_next].a_field.qrie_prev	\
	  = (a_ind);							\
	(a_arr)[a_elmind].a_field.qrie_next = (a_ind);			\
    } while (0)

#define qri_meld(a_arr, a_ind_a, a_ind_b, a_field)			\
    do									\
    {									\
	uint32_t t;							\
	(a_arr)[(a_arr)[a_ind_a].a_field.qrie_prev].a_field.qrie_next	\
	  = (a_ind_b);							\
	(a_arr)[(a_arr)[a_ind_b].a_field.qrie_prev].a_field.qrie_next	\
	  = (a_ind_a);							\
	t = (a_arr)[a_ind_a].a_field.qrie_prev;				\
	(a_arr)[a_ind_a].a_field.qrie_prev				\
	  = (a_arr)[a_ind_b].a_field.qrie_prev;				\
	(a_arr)[a_ind_b].a_field.qrie_prev = t;				\
    } while (0)

/* qri_meld() and qri_split() are functionally equivalent, so there's no need to
 * have two copies of the code. */
#define qri_split(a_arr, a_ind_a, a_ind_b, a_field)			\
    qri_meld((a_arr), (a_ind_a), (a_ind_b), a_field)

#define qri_remove(a_arr, a_ind, a_field)				\
    do									\
    {									\
	(a_arr)[(a_arr)[a_ind].a_field.qrie_prev].a_field.qrie_next	\
	  = (a_arr)[a_ind].a_field.qrie_next;				\
	(a_arr)[(a_arr)[a_ind].a_field.qrie_next].a_field.qrie_prev	\
	  = (a_arr)[a_ind].a_field.qrie_prev;				\
	(a_arr)[a_ind].a_field.qrie_next = (a_ind);			\
	(a_arr)[a_ind].a_field.qrie_prev = (a_ind);			\
    } while (0)

#define qri_foreach(var, a_arr, a_ind, a_field)				\
    for ((var) = (a_ind);						\
	 (var) != UINT_MAX;						\
	 (var) = (((a_arr)[var].a_field.qrie_next != (a_ind))		\
	 ? (a_arr)[var].a_field.qrie_next : UINT_MAX))

#define qri_others_foreach(var, a_arr, a_ind, a_field)			\
    for ((var) = ((a_ind) != UINT_MAX)					\
		  ? (a_arr)[a_ind].a_field.qrie_next : UINT_MAX;	\
	 (var) != (a_ind);						\
	 (var) = (a_arr)[var].a_field.qrie_next)

#define qri_reverse_foreach(var, a_arr, a_ind, a_field)			\
    for ((var) = ((a_ind) != UINT_MAX)					\
	       ? qri_prev(a_arr, a_ind, a_field) : UINT_MAX;		\
	 (var) != UINT_MAX;						\
	 (var) = (((var) != (a_ind))					\
	 ? (a_arr)[var].a_field.qrie_prev : UINT_MAX))

#define qri_others_reverse_foreach(var, a_arr, a_ind, a_field)		\
    for ((var) = ((a_ind) != UINT_MAX)					\
		  ? (a_arr)[a_ind].a_field.qrie_prev : UINT_MAX;	\
	 (var) != (a_ind);						\
	 (var) = (a_arr)[var].a_field.qrie_prev)
