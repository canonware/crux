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
    uint32_t CxQrieNext;							\
    uint32_t CxQriePrev;							\
}

/* Ring functions. */
#define CxmQriNew(a_arr, a_ind, a_field)					\
    do									\
    {									\
	(a_arr)[a_ind].a_field.CxQrieNext = (a_ind);			\
	(a_arr)[a_ind].a_field.CxQriePrev = (a_ind);			\
    } while (0)

#define CxmQriNext(a_arr, a_ind, a_field)					\
    ((a_arr)[a_ind].a_field.CxQrieNext)

#define CxmQriPrev(a_arr, a_ind, a_field)					\
    ((a_arr)[a_ind].a_field.CxQriePrev)

#define CxmQriBeforeInsert(a_arr, a_elmind, a_ind, a_field)		\
    do									\
    {									\
	(a_arr)[a_ind].a_field.CxQriePrev				\
	  = (a_arr)[a_elmind].a_field.CxQriePrev;			\
	(a_arr)[a_ind].a_field.CxQrieNext = (a_elmind);			\
	(a_arr)[(a_arr)[a_ind].a_field.CxQriePrev].a_field.CxQrieNext	\
	  = (a_ind);							\
	(a_arr)[a_elmind].a_field.CxQriePrev = (a_ind);			\
    } while (0)

#define CxmQriAfterInsert(a_arr, a_elmind, a_ind, a_field)		\
    do									\
    {									\
	(a_arr)[a_ind].a_field.CxQrieNext				\
	  = (a_arr)[a_elmind].a_field.CxQrieNext;			\
	(a_arr)[a_ind].a_field.CxQriePrev = (a_elmind);			\
	(a_arr)[(a_arr)[a_ind].a_field.CxQrieNext].a_field.CxQriePrev	\
	  = (a_ind);							\
	(a_arr)[a_elmind].a_field.CxQrieNext = (a_ind);			\
    } while (0)

#define CxmQriMeld(a_arr, a_ind_a, a_ind_b, a_field)			\
    do									\
    {									\
	uint32_t t;							\
	(a_arr)[(a_arr)[a_ind_a].a_field.CxQriePrev].a_field.CxQrieNext	\
	  = (a_ind_b);							\
	(a_arr)[(a_arr)[a_ind_b].a_field.CxQriePrev].a_field.CxQrieNext	\
	  = (a_ind_a);							\
	t = (a_arr)[a_ind_a].a_field.CxQriePrev;				\
	(a_arr)[a_ind_a].a_field.CxQriePrev				\
	  = (a_arr)[a_ind_b].a_field.CxQriePrev;				\
	(a_arr)[a_ind_b].a_field.CxQriePrev = t;				\
    } while (0)

/* CxmQriMeld() and CxmQriSplit() are functionally equivalent, so there's no need to
 * have two copies of the code. */
#define CxmQriSplit(a_arr, a_ind_a, a_ind_b, a_field)			\
    CxmQriMeld((a_arr), (a_ind_a), (a_ind_b), a_field)

#define CxmQriRemove(a_arr, a_ind, a_field)				\
    do									\
    {									\
	(a_arr)[(a_arr)[a_ind].a_field.CxQriePrev].a_field.CxQrieNext	\
	  = (a_arr)[a_ind].a_field.CxQrieNext;				\
	(a_arr)[(a_arr)[a_ind].a_field.CxQrieNext].a_field.CxQriePrev	\
	  = (a_arr)[a_ind].a_field.CxQriePrev;				\
	(a_arr)[a_ind].a_field.CxQrieNext = (a_ind);			\
	(a_arr)[a_ind].a_field.CxQriePrev = (a_ind);			\
    } while (0)

#define CxmQriForeach(var, a_arr, a_ind, a_field)				\
    for ((var) = (a_ind);						\
	 (var) != UINT_MAX;						\
	 (var) = (((a_arr)[var].a_field.CxQrieNext != (a_ind))		\
	 ? (a_arr)[var].a_field.CxQrieNext : UINT_MAX))

#define CxmQriOthersForeach(var, a_arr, a_ind, a_field)			\
    for ((var) = ((a_ind) != UINT_MAX)					\
		  ? (a_arr)[a_ind].a_field.CxQrieNext : UINT_MAX;	\
	 (var) != (a_ind);						\
	 (var) = (a_arr)[var].a_field.CxQrieNext)

#define CxmQriReverseForeach(var, a_arr, a_ind, a_field)			\
    for ((var) = ((a_ind) != UINT_MAX)					\
	       ? CxmQriPrev(a_arr, a_ind, a_field) : UINT_MAX;		\
	 (var) != UINT_MAX;						\
	 (var) = (((var) != (a_ind))					\
	 ? (a_arr)[var].a_field.CxQriePrev : UINT_MAX))

#define CxmQriOthersReverseForeach(var, a_arr, a_ind, a_field)		\
    for ((var) = ((a_ind) != UINT_MAX)					\
		  ? (a_arr)[a_ind].a_field.CxQriePrev : UINT_MAX;	\
	 (var) != (a_ind);						\
	 (var) = (a_arr)[var].a_field.CxQriePrev)
