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

typedef struct CxsTreeObject CxtTreeObject;
typedef struct CxsNodeObject CxtNodeObject;
typedef struct CxsEdgeObject CxtEdgeObject;
typedef struct CxsRingObject CxtRingObject;

struct CxsTreeObject
{
    PyObject_HEAD
    CxtTr *tr;
    bool GcCleared;

#define CxmTreeObjectAuxCount 1
#define CxmTreeObjectAuxRf 0
    void *aux[CxmTreeObjectAuxCount];
};

struct CxsNodeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
    bool GcCleared;

#define CxmNodeObjectAuxCount 1
    void *aux[CxmNodeObjectAuxCount];
};

struct CxsEdgeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtRingObject *ringA;
    CxtRingObject *ringB;
    CxtTrEdge edge;
    bool GcDetached:1;
    bool GcCleared:1;

#define CxmEdgeObjectAuxCount 2
#define CxmEdgeObjectAuxMp 0
#define CxmEdgeObjectAuxRf 1
    void *aux[CxmEdgeObjectAuxCount];
};

struct CxsRingObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtEdgeObject *edge;
    CxtTrRing ring;
    bool GcDetached:1;
    bool GcCleared:1;

#define CxmRingObjectAuxCount 1
#define CxmRingObjectAuxMp 0
    void *aux[CxmRingObjectAuxCount];
};

extern PyObject *CxgTreeException;
extern PyObject *CxgTreeValueError;
extern PyObject *CxgTreeTypeError;

extern PyObject *CxgEdgeException;
extern PyObject *CxgEdgeValueError;

extern PyObject *CxgRingException;
extern PyObject *CxgRingValueError;

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

#if (!defined(CxmUseInlines))
void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd);
void
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmTreeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmTreeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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

#if (!defined(CxmUseInlines))
void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd);
void
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmNodeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmNodeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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

#if (!defined(CxmUseInlines))
void *
CxEdgeAuxGet(CxtEdgeObject *self, unsigned aInd);
void
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxEdgeAuxGet(CxtEdgeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmEdgeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmEdgeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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

#if (!defined(CxmUseInlines))
void *
CxRingAuxGet(CxtRingObject *self, unsigned aInd);
void
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxRingAuxGet(CxtRingObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmRingObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmRingObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif
