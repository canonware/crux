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
} CxtNodeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrEdge edge;
} CxtEdgeObject;

void
CxTreeInit(void);
void
CxNodeInit(void);
void
CxEdgeInit(void);

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
CxNodeEdge(CxtNodeObject *self);
PyObject *
CxNodeDegree(CxtNodeObject *self);

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *a_tree);
PyObject *
CxEdgeTree(CxtEdgeObject *self);
PyObject *
CxEdgeNode(CxtEdgeObject *self, uint32_t aInd);
PyObject *
CxEdgeNodePargs(CxtEdgeObject *self, PyObject *args);
void
CxEdgeNext(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
	   uint32_t *rNextEnd);
PyObject *
CxEdgeNextPargs(CxtEdgeObject *self, PyObject *args);
void
CxEdgePrev(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
	   uint32_t *rPrevEnd);
PyObject *
CxEdgePrevPargs(CxtEdgeObject *self, PyObject *args);
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
