/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 * Public interface for the bhp (binomial heap) class, and its helper, the bhpi
 * class.
 *
 ******************************************************************************/

/* Typedef to allow easy function pointer passing. */
typedef int32_t bhp_prio_comp_t (const void *, const void *);

/* Pseudo-opaque types. */
typedef struct cw_bhp_s cw_bhp_t;
typedef struct cw_bhpi_s cw_bhpi_t;

struct cw_bhp_s
{
#ifdef CW_DBG
    uint32_t magic;
#define CW_BHP_MAGIC 0xbf917ca1
#endif
    bool is_malloced;
    cw_bhpi_t *head;
    uint64_t num_nodes;
    bhp_prio_comp_t *priority_compare;
};

struct cw_bhpi_s
{
#ifdef CW_DBG
    uint32_t magic;
#define CW_BHPI_MAGIC 0xbf90ec15
#endif
    bool is_malloced;

    const void *priority;
    const void *data;
    struct cw_bhpi_s *parent;
    struct cw_bhpi_s *child;
    struct cw_bhpi_s *sibling;
    uint32_t degree;
};

cw_bhp_t *
bhp_new(cw_bhp_t *a_bhp, bhp_prio_comp_t *a_prio_comp);

void
bhp_delete(cw_bhp_t *a_bhp);

void
bhp_dump(cw_bhp_t *a_bhp);

void
bhp_insert(cw_bhp_t *a_bhp, const void *a_priority, const void *a_data,
	   cw_bhpi_t *a_bhpi);

bool
bhp_min_find(cw_bhp_t *a_bhp, void **r_priority, void **r_data);

bool
bhp_min_del(cw_bhp_t *a_bhp, void **r_priority, void **r_data,
	    cw_bhpi_t **r_bhpi);

uint64_t
bhp_size_get(cw_bhp_t *a_bhp);

void
bhp_union(cw_bhp_t *a_a, cw_bhp_t *a_b);

int32_t
bhp_uint32_priority_compare(const void *a_a, const void *a_b);

int32_t
bhp_int32_priority_compare(const void *a_a, const void *a_b);

int32_t
bhp_uint64_priority_compare(const void *a_a, const void *a_b);
