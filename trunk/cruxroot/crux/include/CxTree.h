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
    CxtTr *tr;
} CxtCxtTreeObject;

typedef struct
{
    PyObject_HEAD
    CxtCxtTreeObject *tree;
    CxtTrNode node;
} CxtCxtNodeObject;

typedef struct
{
    PyObject_HEAD
    CxtCxtTreeObject *tree;
    CxtTrEdge edge;
} CxtCxtEdgeObject;

void
CxTreeInit(void);
void
CxNodeInit(void);
void
CxEdgeInit(void);

CxtCxtTreeObject *
CxTreeNew(void);
PyObject *
CxTreeNtaxaGet(CxtCxtTreeObject *self);
PyObject *
CxTreeNedgesCget(CxtCxtTreeObject *self);
PyObject *
CxTreeBaseGet(CxtCxtTreeObject *self);
void
CxTreeBaseSetCargs(CxtCxtTreeObject *self, CxtCxtNodeObject *a_node);
PyObject *
CxTreeBaseSet(CxtCxtTreeObject *self, PyObject *args);

CxtCxtNodeObject *
CxNodeNew(CxtCxtTreeObject *a_tree);
PyObject *
CxNodeTree(CxtCxtNodeObject *self);
PyObject *
CxNodeTaxonNumGet(CxtCxtNodeObject *self);
void
CxNodeTaxonNumSetCargs(CxtCxtNodeObject *self, uint32_t a_taxon_num);
PyObject *
CxNodeTaxonNumSet(CxtCxtNodeObject *self, PyObject *args);
PyObject *
CxNodeEdge(CxtCxtNodeObject *self);
PyObject *
CxNodeDegree(CxtCxtNodeObject *self);

CxtCxtEdgeObject *
CxEdgeNew(CxtCxtTreeObject *a_tree);
PyObject *
CxEdgeTree(CxtCxtEdgeObject *self);
PyObject *
CxEdgeNodeCargs(CxtCxtEdgeObject *self, uint32_t a_ind);
PyObject *
CxEdgeNode(CxtCxtEdgeObject *self, PyObject *args);
void
CxEdgeNextCargs(CxtCxtEdgeObject *self, uint32_t a_ind, CxtCxtEdgeObject **r_edge,
		uint32_t *r_next_end);
PyObject *
CxEdgeNext(CxtCxtEdgeObject *self, PyObject *args);
void
CxEdgePrevCargs(CxtCxtEdgeObject *self, uint32_t a_ind, CxtCxtEdgeObject **r_edge,
		uint32_t *r_prev_end);
PyObject *
CxEdgePrev(CxtCxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeLengthGet(CxtCxtEdgeObject *self);
void
CxEdgeLengthSetCargs(CxtCxtEdgeObject *self, double a_length);
PyObject *
CxEdgeLengthSet(CxtCxtEdgeObject *self, PyObject *args);
void
CxEdgeAttachCargs(CxtCxtEdgeObject *self, CxtCxtNodeObject *a_node_a, CxtCxtNodeObject *a_node_b);
PyObject *
CxEdgeAttach(CxtCxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeDetach(CxtCxtEdgeObject *self);
