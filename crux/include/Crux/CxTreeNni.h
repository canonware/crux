//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

bool
CxTreeNniBEdgeSetsGet(CxtTreeObject *self, CxtEdgeObject *aEdge,
		      CxtEdgeObject ***rSetA, unsigned *rNSetA,
		      CxtRingObject **rRingA,
		      CxtEdgeObject ***rSetB, unsigned *rNSetB,
		      CxtRingObject **rRingB);

bool
CxTreeNni(CxtTreeObject *self, CxtEdgeObject *aBisect,
	  CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB);
PyObject *
CxTreeNniPargs(CxtTreeObject *self, PyObject *args);

bool
CxTreeNniNEdgesGet(CxtTreeObject *self, unsigned *rNEdges);

bool
CxTreeNniEdgeGet(CxtTreeObject *self, unsigned aEdge, CxtEdgeObject **rEdge);

bool
CxTreeNniEdgeOffset(CxtTreeObject *self, unsigned aEdge, unsigned *rOffset);

bool
CxTreeNniNNeighborsGet(CxtTreeObject *self, unsigned *rNNeighbors);
PyObject *
CxTreeNniNNeighborsGetPargs(CxtTreeObject *self);

bool
CxTreeNniNeighborGet(CxtTreeObject *self, unsigned aNeighbor,
		     CxtEdgeObject **rBisect, CxtEdgeObject **rReconnectA,
		     CxtEdgeObject **rReconnectB);
PyObject *
CxTreeNniNeighborGetPargs(CxtTreeObject *self, PyObject *args);
