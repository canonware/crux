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
} CxtTreeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
    bool valid;
} CxtNodeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrEdge edge;
    bool valid;
} CxtEdgeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtEdgeObject *edge;
    CxtTrRing ring;
} CxtRingObject;

void
CxTreeInit(void);
void
CxNodeInit(void);
void
CxEdgeInit(void);
void
CxRingInit(void);

/* Tree. */
CxtTreeObject *
CxTreeNew(void);
PyObject *
CxTreeNtaxaGet(CxtTreeObject *self);
PyObject *
CxTreeNedgesCget(CxtTreeObject *self);
PyObject *
CxTreeBaseGet(CxtTreeObject *self);
void
CxTreeBaseSet(CxtTreeObject *self, CxtNodeObject *aNode);
PyObject *
CxTreeBaseSetPargs(CxtTreeObject *self, PyObject *args);

/* Node. */
CxtNodeObject *
CxNodeNew(CxtTreeObject *a_tree);
PyObject *
CxNodeTree(CxtNodeObject *self);
PyObject *
CxNodeTaxonNumGet(CxtNodeObject *self);
void
CxNodeTaxonNumSet(CxtNodeObject *self, uint32_t a_taxon_num);
PyObject *
CxNodeTaxonNumSetPargs(CxtNodeObject *self, PyObject *args);
PyObject *
CxNodeRing(CxtNodeObject *self);
PyObject *
CxNodeDegree(CxtNodeObject *self);

/* Edge. */
CxtEdgeObject *
CxEdgeNew(CxtTreeObject *a_tree);
PyObject *
CxEdgeTree(CxtEdgeObject *self);
PyObject *
CxEdgeRingsGet(CxtEdgeObject *self);
PyObject *
CxEdgeLengthGet(CxtEdgeObject *self);
void
CxEdgeLengthSet(CxtEdgeObject *self, double aLength);
PyObject *
CxEdgeLengthSetPargs(CxtEdgeObject *self, PyObject *args);
void
CxEdgeAttach(CxtEdgeObject *self, CxtNodeObject *aNodeA,
	     CxtNodeObject *aNodeB);
PyObject *
CxEdgeAttachPargs(CxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeDetach(CxtEdgeObject *self);

/* Ring. */
CxtRingObject *
CxRingNew(CxtEdgeObject *aEdge, uint32_t aEnd);
PyObject *
CxRingTree(CxtRingObject *self);
PyObject *
CxRingNode(CxtRingObject *self);
PyObject *
CxRingEdge(CxtRingObject *self);
PyObject *
CxRingOther(CxtRingObject *self);
PyObject *
CxRingNext(CxtRingObject *self);
PyObject *
CxRingPrev(CxtRingObject *self);
