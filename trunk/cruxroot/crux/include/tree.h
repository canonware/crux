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

PyObject *
tree_ntaxa_get(TreeObject *self);
PyObject *
tree_nedges_get(TreeObject *self);
PyObject *
tree_base_get(TreeObject *self);
PyObject *
tree_base_set(TreeObject *self, PyObject *args);

PyObject *
node_tree(NodeObject *self);
PyObject *
node_taxon_num_get(NodeObject *self);
PyObject *
node_taxon_num_set(NodeObject *self, PyObject *args);
PyObject *
node_edge(NodeObject *self);
PyObject *
node_degree(NodeObject *self);

PyObject *
edge_tree(EdgeObject *self);
PyObject *
edge_node(EdgeObject *self, PyObject *args);
PyObject *
edge_next(EdgeObject *self, PyObject *args);
PyObject *
edge_prev(EdgeObject *self, PyObject *args);
PyObject *
edge_length_get(EdgeObject *self);
PyObject *
edge_length_set(EdgeObject *self, PyObject *args);
PyObject *
edge_attach(EdgeObject *self, PyObject *args);
PyObject *
edge_detach(EdgeObject *self);
