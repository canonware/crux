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

typedef struct
{
    PyObject_HEAD
    cw_tr_t *tr;
} TreeObject;

typedef struct
{
    PyObject_HEAD
    TreeObject *tree;
    cw_tr_node_t node;
} NodeObject;

typedef struct
{
    PyObject_HEAD
    TreeObject *tree;
    cw_tr_edge_t edge;
} EdgeObject;

void
crux_tree_init(void);
void
crux_node_init(void);
void
crux_edge_init(void);

TreeObject *
tree_new(void);
PyObject *
tree_ntaxa_get(TreeObject *self);
PyObject *
tree_nedges_get(TreeObject *self);
PyObject *
tree_base_get(TreeObject *self);
void
tree_base_set_cargs(TreeObject *self, NodeObject *a_node);
PyObject *
tree_base_set(TreeObject *self, PyObject *args);

NodeObject *
node_new(TreeObject *a_tree);
PyObject *
node_tree(NodeObject *self);
PyObject *
node_taxon_num_get(NodeObject *self);
void
node_taxon_num_set_cargs(NodeObject *self, uint32_t a_taxon_num);
PyObject *
node_taxon_num_set(NodeObject *self, PyObject *args);
PyObject *
node_edge(NodeObject *self);
PyObject *
node_degree(NodeObject *self);

EdgeObject *
edge_new(TreeObject *a_tree);
PyObject *
edge_tree(EdgeObject *self);
PyObject *
edge_node_cargs(EdgeObject *self, uint32_t a_ind);
PyObject *
edge_node(EdgeObject *self, PyObject *args);
void
edge_next_cargs(EdgeObject *self, uint32_t a_ind, EdgeObject **r_edge,
		uint32_t *r_next_end);
PyObject *
edge_next(EdgeObject *self, PyObject *args);
void
edge_prev_cargs(EdgeObject *self, uint32_t a_ind, EdgeObject **r_edge,
		uint32_t *r_prev_end);
PyObject *
edge_prev(EdgeObject *self, PyObject *args);
PyObject *
edge_length_get(EdgeObject *self);
void
edge_length_set_cargs(EdgeObject *self, double a_length);
PyObject *
edge_length_set(EdgeObject *self, PyObject *args);
void
edge_attach_cargs(EdgeObject *self, NodeObject *a_node_a, NodeObject *a_node_b);
PyObject *
edge_attach(EdgeObject *self, PyObject *args);
PyObject *
edge_detach(EdgeObject *self);
