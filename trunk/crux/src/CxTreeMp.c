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

#include "../include/_cruxmodule.h"

typedef struct CxsTreeMpData CxtTreeMpData;
typedef struct CxsTreeMpPs CxtTreeMpPs;

struct CxsTreeMpData
{
    // XXX Use these.
    //PyObject *cTMatrix;
    //uint64_t cTMatrixSeq;

    unsigned edgeAuxInd;
    unsigned ringAuxInd;

    bool eliminateUninformative;
    char **taxa;
    unsigned nTaxa;
    unsigned nChars;
    bool *charsMask;
};

// Character (in the systematics sense of the word).
typedef char CxtTreeMpC;

// Partial parsimony score information.
struct CxsTreeMpPs
{
    // Parent which most recently used this node's partial score when caching
    // its results.  Both children must still point to the parent in order for
    // the cached results to be valid.
    CxtTreeMpPs *parent;

    // Sum of the subtree scores, and this node's score, given particular
    // children.  In order for this to be useful, both childrens' parent
    // pointers must still point to this node.
    unsigned subtreesScore;
    unsigned nodeScore;

    // chars points to an array of Fitch parsimony state sets.  Each element in
    // the array contains a bitmap representation of a subset of {ACGT} in the 4
    // least significant bits.  T is the least significant bit.  1 means that a
    // nucleotide is in the set.
    //
    // There are nchars character state sets.
    //
    // achars is the actual allocation, which is padded in order to
    // be able to guarantee that chars is 16 byte-aligned.
    CxtTreeMpC *chars;
    unsigned nChars;
    CxtTreeMpC *aChars;
};

static bool
CxpTreeMpDataGet(CxtTreeObject *self, CxtTreeMpData **rData);

//==============================================================================
//
// CxTreeMpPs.
//

CxmpInline CxtTreeMpPs *
CxpTreeMpPsNew(void)
{
    CxtTreeMpPs *rVal;

    rVal = (CxtTreeMpPs *) malloc(sizeof(CxtTreeMpPs));
    if (rVal == NULL)
    {
	PyErr_NoMemory();
	goto RETURN;
    }

    rVal->parent = NULL;
    rVal->chars = NULL;

    RETURN:
    return rVal;
}

CxmpInline void
CxpTreeMpPsDelete(CxtTreeMpPs *aPs)
{
    if (aPs->chars != NULL)
    {
	free(aPs->chars);
    }

    free(aPs);
}

#if (0) // Unused (so far).
CxmpInline CxtTreeMpC
CxpTreePsCharGet(CxtTreeMpPs *aPs, unsigned aOffset)
{
    CxtTreeMpC rVal;

    CxmCheckPtr(aPs->chars);
    CxmAssert(aOffset < aPs->nchars);

    rVal = aPs->chars[aOffset >> 1];
    rVall >>= ((aOffset & 1) * 4);

    return rVal;
}
#endif

CxmpInline void
CxpTreePsCharSet(CxtTreeMpPs *aPs, CxtTreeMpC aChar, unsigned aOffset)
{
    CxmCheckPtr(aPs->chars);
    CxmAssert((aChar & 0xfU) == aChar);
    CxmAssert(aOffset < aPs->nChars + ((32 - (aPs->nChars & 0x1fU)) & 0x1fU));

    if ((aOffset & 1) == 0)
    {
	aPs->chars[aOffset >> 1]
	    = (aChar << 4)
	    | (aPs->chars[aOffset >> 1] & 0xfU);
    }
    else
    {
	aPs->chars[aOffset >> 1]
	    = (aPs->chars[aOffset >> 1] & 0xf0U)
	    | aChar;
    }
}

static bool
CxpTreePsPrepare(CxtTreeMpData *aData, CxtTreeMpPs *aPs)
{
    bool rVal;

    // Clean up old character vector if it isn't the right size for
    // aData->nchars characters.
    if (aPs->chars != NULL && aPs->nChars != aData->nChars)
    {
	free(aPs->aChars);
    }

    // Allocate character vector if necessary.
    if (aPs->chars == NULL)
    {
	if (aData->nChars != 0)
	{
	    unsigned nPad, i;

	    // Calculate the number of pad bytes to append, such that the total
	    // number of bytes is a multiple of 16 (total number of taxonomical
	    // characters is a multiple of 32).
	    nPad = (32 - (aData->nChars & 0x1fU)) & 0x1fU;
	    CxmAssert(((aData->nChars + nPad) & 0x1fU) == 0);

	    // Tack on 8 bytes; all modern systems provide at least 8 byte
	    // alignment.
	    aPs->aChars = (CxtTreeMpC *) malloc(sizeof(CxtTreeMpC)
						* (((aData->nChars + nPad)
						    >> 1))
						+ 8);
	    if (aPs->aChars == NULL)
	    {
		PyErr_NoMemory();
		rVal = true;
		goto RETURN;
	    }

	    // Make sure that aPs->chars is 16 byte-aligned.  Assume that
	    // aPs->aChars is at least 8 byte-aligned.
	    aPs->chars = &aPs->aChars[((unsigned) aPs->aChars) & 0xfU];

	    aPs->nChars = aData->nChars + nPad;

	    // Set all pad characters to {ACGT}.  This allows the pad characters
	    // to be calculated along with the actual characters, without
	    // affecting the score.
	    for (i = aData->nChars; i < aData->nChars + nPad; i++)
	    {
		CxpTreePsCharSet(aPs, 0xfU, i);
	    }
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

//
// End CxTreeMpPs.
//
//==============================================================================

static void
CxpTreeMpCleanupFinal(CxtTreeObject *aTree, void *aData, unsigned aInd)
{
    CxtTreeMpData *data = (CxtTreeMpData *) aData;

    if (data->charsMask != NULL)
    {
	free(data->charsMask);
    }

    if (data->taxa != NULL)
    {
	free(data->taxa);
    }

    free(data);
}

static bool
CxpTreeMpInitEdge(CxtEdgeObject *aEdge, unsigned aInd)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtTreeMpPs *ps;

    if (CxpTreeMpDataGet(CxEdgeTree(aEdge), &data))
    {
	rVal = true;
	goto RETURN;
    }

    ps = CxEdgeAuxGet(aEdge, aInd);
    if (ps == NULL)
    {
	ps = CxpTreeMpPsNew();
	if (ps == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	if (CxEdgeAuxSet(aEdge, aInd, ps))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    if (CxpTreePsPrepare(data, ps))
    {
	rVal = true;
	goto RETURN;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static void
CxpTreeMpCleanupEdge(CxtEdgeObject *aEdge, void *aData, unsigned aInd)
{
    CxtTreeMpPs *ps = (CxtTreeMpPs *) aData;

    CxpTreeMpPsDelete(ps);
    CxEdgeAuxSet(aEdge, aInd, NULL);
}

static bool
CxpTreeMpInitRing(CxtRingObject *aRing, unsigned aInd)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtTreeMpPs *ps;
    uint32_t taxonNum;

    if (CxpTreeMpDataGet(CxRingTree(aRing), &data))
    {
	rVal = true;
	goto RETURN;
    }

    ps = CxRingAuxGet(aRing, aInd);
    if (ps == NULL)
    {
	ps = CxpTreeMpPsNew();
	if (ps == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	if (CxRingAuxSet(aRing, aInd, ps))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    if (CxpTreePsPrepare(data, ps))
    {
	rVal = true;
	goto RETURN;
    }

    // If this is a leaf node, initialize the character state sets and scores.
    taxonNum = CxNodeTaxonNumGet(CxRingNode(aRing));
    if (taxonNum != CxmTrNodeNone)
    {
	unsigned i, j;
	char *chars;

	ps->subtreesScore = 0;
	ps->nodeScore = 0;

	CxmAssert(taxonNum < data->nTaxa);
	chars = data->taxa[taxonNum];
	for (i = j = 0; i < data->nChars; i++)
	{
	    // Ignore uninformative characters.
	    if (data->charsMask[i] == false)
	    {
		continue;
	    }

	    switch (chars[i])
	    {
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		    // Treat gaps as uncertainty.  This isn't the only way to do
		    // things, and may need to be made configurable.
		case '-':
		{
		    CxpTreePsCharSet(ps, 0xf, j);
		    break;
		}
		case 'V':
		case 'v':
		{
		    CxpTreePsCharSet(ps, 0xe, j);
		    break;
		}
		case 'H':
		case 'h':
		{
		    CxpTreePsCharSet(ps, 0xd, j);
		    break;
		}
		case 'M':
		case 'm':
		{
		    CxpTreePsCharSet(ps, 0xc, j);
		    break;
		}
		case 'D':
		case 'd':
		{
		    CxpTreePsCharSet(ps, 0xb, j);
		    break;
		}
		case 'R':
		case 'r':
		{
		    CxpTreePsCharSet(ps, 0xa, j);
		    break;
		}
		case 'W':
		case 'w':
		{
		    CxpTreePsCharSet(ps, 0x9, j);
		    break;
		}
		case 'A':
		case 'a':
		{
		    CxpTreePsCharSet(ps, 0x8, j);
		    break;
		}
		case 'B':
		case 'b':
		{
		    CxpTreePsCharSet(ps, 0x7, j);
		    break;
		}
		case 'S':
		case 's':
		{
		    CxpTreePsCharSet(ps, 0x6, j);
		    break;
		}
		case 'Y':
		case 'y':
		{
		    CxpTreePsCharSet(ps, 0x5, j);
		    break;
		}
		case 'C':
		case 'c':
		{
		    CxpTreePsCharSet(ps, 0x4, j);
		    break;
		}
		case 'K':
		case 'k':
		{
		    CxpTreePsCharSet(ps, 0x3, j);
		    break;
		}
		case 'G':
		case 'g':
		{
		    CxpTreePsCharSet(ps, 0x2, j);
		    break;
		}
		case 'T':
		case 't':
		{
		    CxpTreePsCharSet(ps, 0x1, j);
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	    j++;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static void
CxpTreeMpCleanupRing(CxtRingObject *aRing, void *aData, unsigned aInd)
{
    CxtTreeMpPs *ps = (CxtTreeMpPs *) aData;

    CxpTreeMpPsDelete(ps);
    CxRingAuxSet(aRing, aInd, NULL);
}

static bool
CxpTreeMpDataGet(CxtTreeObject *self, CxtTreeMpData **rData)
{
    bool rVal;
    unsigned treeAuxInd;
    CxtTreeMpData *data;

    // Get aux indices for MP data.
    if (CxTreeAuxSearch(self, "MP", &treeAuxInd))
    {
	// No aux registration.
	data = (CxtTreeMpData *) malloc(sizeof(CxtTreeMpData));
	if (data == NULL)
	{
	    PyErr_NoMemory();
	    rVal = true;
	    goto RETURN;
	}

	// Create MP-specific aux mappings for edges and rings.
	if (CxTreeEdgeAuxRegister(self, "MP", NULL,
				  CxpTreeMpInitEdge, NULL,
				  CxpTreeMpCleanupEdge,
				  &data->edgeAuxInd)
	    || CxTreeRingAuxRegister(self, "MP", NULL,
				     CxpTreeMpInitRing, NULL,
				     CxpTreeMpCleanupRing,
				     &data->ringAuxInd)
	    || CxTreeAuxRegister(self, "MP", (void *) data,
				 CxpTreeMpCleanupFinal, NULL, &treeAuxInd))
	{
	    free(data);
	    rVal = true;
	    goto RETURN;
	}

	// Initialize data.
	data->eliminateUninformative = false;
	data->taxa = NULL;
	data->nTaxa = 0;
	data->nChars = 0;
	data->charsMask = NULL;
    }
    else
    {
	data = (CxtTreeMpData *) CxTreeAuxData(self, treeAuxInd);
    }

    *rData = data;
    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpTreeMpPrepareRecurse(CxtTreeObject *self, CxtTreeMpData *aData,
			CxtRingObject *aRing)
{
    bool rVal;
    CxtEdgeObject *edge;
    CxtRingObject *ring;

    // Prepare aRing.
    if (CxpTreeMpInitRing(aRing, aData->ringAuxInd))
    {
	rVal = true;
	goto RETURN;
    }

    // Recurse into subtrees.
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	// Prepare edge before recursing.
	edge = CxRingEdge(ring);
	if (CxpTreeMpInitEdge(edge, aData->edgeAuxInd))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Prepare ring.
	if (CxpTreeMpInitRing(ring, aData->ringAuxInd))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Recurse.
	CxpTreeMpPrepareRecurse(self, aData, ring);
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpTreeMpPrepare(CxtTreeObject *self, bool aUninformativeEliminate,
		 char **aTaxa, uint32_t aNTaxa, uint32_t aNChars)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtNodeObject *node;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    if ((node = CxTreeBaseGet(self)) != NULL)
    {
	CxtRingObject *ringFirst, *ring;
	CxtEdgeObject *edge;
	uint32_t i, ninformative;

	// Set fields in data.
	data->eliminateUninformative = aUninformativeEliminate;
	data->taxa = aTaxa;
	data->nTaxa = aNTaxa;
	data->nChars = aNChars;

	data->charsMask = (bool *) malloc(sizeof(bool) * aNChars);
	if (data->charsMask == NULL)
	{
	    PyErr_NoMemory();
	    rVal = true;
	    goto RETURN;
	}

	if (aUninformativeEliminate)
	{
	    uint32_t codes[15];
	    uint32_t j, k, x, y;

	    // Preprocess the character data.  Eliminate uninformative
	    // characters, but keep track of their contribution to the parsimony
	    // score, were they to be left in.
	    ninformative = 0;
	    for (i = 0; i < aNChars; i++)
	    {
		for (k = 0; k < 15; k++)
		{
		    codes[k] = 0;
		}

		for (j = 0; j < aNTaxa; j++)
		{
		    switch (aTaxa[j][i])
		    {
			case 'N':
			case 'n':
			case 'X':
			case 'x':
			case '-':
			{
			    // Treat gaps as uncertainty.  This isn't the only
			    // way to do things, and may need to be made
			    // configurable.
			    break;
			}
			case 'V':
			case 'v':
			{
			    codes[14]++;
			    break;
			}
			case 'H':
			case 'h':
			{
			    codes[13]++;
			    break;
			}
			case 'M':
			case 'm':
			{
			    codes[12]++;
			    break;
			}
			case 'D':
			case 'd':
			{
			    codes[11]++;
			    break;
			}
			case 'R':
			case 'r':
			{
			    codes[10]++;
			    break;
			}
			case 'W':
			case 'w':
			{
			    codes[9]++;
			    break;
			}
			case 'A':
			case 'a':
			{
			    codes[8]++;
			    break;
			}
			case 'B':
			case 'b':
			{
			    codes[7]++;
			    break;
			}
			case 'S':
			case 's':
			{
			    codes[6]++;
			    break;
			}
			case 'Y':
			case 'y':
			{
			    codes[5]++;
			    break;
			}
			case 'C':
			case 'c':
			{
			    codes[4]++;
			    break;
			}
			case 'K':
			case 'k':
			{
			    codes[3]++;
			    break;
			}
			case 'G':
			case 'g':
			{
			    codes[2]++;
			    break;
			}
			case 'T':
			case 't':
			{
			    codes[1]++;
			    break;
			}
			default:
			{
			    CxmNotReached();
			}
		    }
		}

		// Count the number of states in which two or more taxa exist.
		data->charsMask[i] = false;
		for (x = 1; x < 15; x++)
		{
		    for (y = 1; y < 15; y++)
		    {
			if ((x & y) == 0 && codes[x] >= 2 && codes[y] >= 2)
			{
			    if (data->charsMask[i] == false)
			    {
				data->charsMask[i] = true;
				ninformative++;
			    }
			}
		    }
		}
	    }
	}
	else
	{
	    // Use all characters, regardless of whether they are informative.
	    ninformative = aNChars;

	    for (i = 0; i < aNChars; i++)
	    {
		data->charsMask[i] = true;
	    }
	}

	// Prepare the tree.
	ringFirst = CxNodeRing(node);
	if (ringFirst != NULL)
	{
	    ring = ringFirst;
	    do
	    {
		// Prepare edge before recursing.
		edge = CxRingEdge(ring);
		if (CxpTreeMpInitEdge(edge, data->edgeAuxInd))
		{
		    rVal = true;
		    goto RETURN;
		}

		// Prepare ring.
		if (CxpTreeMpInitRing(ring, data->ringAuxInd))
		{
		    rVal = true;
		    goto RETURN;
		}

		// Recurse.
		if (CxpTreeMpPrepareRecurse(self, data, CxRingOther(ring)))
		{
		    rVal = true;
		    goto RETURN;
		}

		ring = CxRingNext(ring);
	    } while (ring != ringFirst);
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeMpPrepare(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *taxa, *tobj;
    uint32_t elim, ntaxa, i, j;
    char **tarr
#ifdef CxmCcSilence
	= NULL
#endif
	;

    if (PyArg_ParseTuple(args, "O!i", &PyList_Type, &taxa, &elim) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    ntaxa = PyList_Size(taxa);

    // Make sure that all taxa have the same number of characters.
    if (ntaxa > 0)
    {
	uint32_t nchars;

	tobj = PyList_GetItem(taxa, 0);
	// Don't worry about raising ValueError error here, since the for
	// loop below will do so.
	if (PyString_Check(tobj))
	{
	    nchars = PyString_Size(tobj);
	}

	// Create an array of string pointers.
	tarr = (char **) malloc(sizeof(char *) * ntaxa);
	if (tarr == NULL)
	{
	    rVal = PyErr_NoMemory();
	    goto RETURN;
	}

	for (i = 0; i < ntaxa; i++)
	{
	    tobj = PyList_GetItem(taxa, i);
	    if (PyString_Check(tobj) == 0 || PyString_Size(tobj) != nchars)
	    {
		free(tarr);
		CxError(CxgTreeValueError,
			"Character string expected");
		rVal = NULL;
		goto RETURN;
	    }
	    tarr[i] = PyString_AsString(tobj);

	    // Make sure characters are valid codes.
	    for (j = 0; j < nchars; j++)
	    {
		switch (tarr[i][j])
		{
		    case 'N':
		    case 'n':
		    case 'X':
		    case 'x':
		    case 'V':
		    case 'v':
		    case 'H':
		    case 'h':
		    case 'M':
		    case 'm':
		    case 'D':
		    case 'd':
		    case 'R':
		    case 'r':
		    case 'W':
		    case 'w':
		    case 'A':
		    case 'a':
		    case 'B':
		    case 'b':
		    case 'S':
		    case 's':
		    case 'Y':
		    case 'y':
		    case 'C':
		    case 'c':
		    case 'K':
		    case 'k':
		    case 'G':
		    case 'g':
		    case 'T':
		    case 't':
		    case '-':
		    {
			break;
		    }
		    default:
		    {
			free(tarr);
			CxError(CxgTreeValueError,
				"Invalid character '%c'", tarr[i][j]);
			rVal = NULL;
			goto RETURN;
		    }
		}
	    }
	}

	// XXX Recurse through the tree and make sure that the taxa are
	// numbered correctly.

	// Do preparation.  tarr is managed internally from here on, so it
	// shouldn't be freed here unless there is an error.
	if (CxpTreeMpPrepare(self, elim, tarr, ntaxa, nchars))
	{
	    free(tarr);
	    rVal = NULL;
	    goto RETURN;
	}
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

static void
CxpTreeMpFinishRecurse(CxtTreeMpData *aData, CxtRingObject *aRing)
{
    CxtRingObject *ring;
    CxtEdgeObject *edge;

    // Clean up aRing.
    CxpTreeMpCleanupRing(aRing, CxRingAuxGet(aRing, aData->ringAuxInd),
			 aData->ringAuxInd);

    // Recurse into subtrees.
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	// Clean up edge before recursing.
	edge = CxRingEdge(ring);
	CxpTreeMpCleanupEdge(edge, CxEdgeAuxGet(edge, aData->edgeAuxInd),
			     aData->edgeAuxInd);

	// Clean up ring.
	CxpTreeMpCleanupRing(ring, CxRingAuxGet(ring, aData->ringAuxInd),
			     aData->ringAuxInd);

	// Recurse.
	CxpTreeMpFinishRecurse(aData, CxRingOther(ring));
    }
}

PyObject *
CxTreeMpFinish(CxtTreeObject *self)
{
    PyObject *rVal;
    CxtNodeObject *base;

    base = CxTreeBaseGet(self);
    if (base != NULL)
    {
	CxtTreeMpData *data;
	CxtRingObject *ringFirst, *ring;
	CxtEdgeObject *edge;

	if (CxpTreeMpDataGet(self, &data))
	{
	    rVal = NULL;
	    goto RETURN;
	}

	// Clean up the tree.
	ringFirst = CxNodeRing(base);
	if (ringFirst != NULL)
	{
	    ring = ringFirst;
	    do
	    {
		// Clean up edge before recursing.
		edge = CxRingEdge(ring);
		CxpTreeMpCleanupEdge(edge, CxEdgeAuxGet(edge, data->edgeAuxInd),
				     data->edgeAuxInd);

		// Clean up ring.
		CxpTreeMpCleanupRing(ring, CxRingAuxGet(ring, data->ringAuxInd),
				     data->ringAuxInd);

		// Recurse.
		CxpTreeMpFinishRecurse(data, CxRingOther(ring));

		ring = CxRingNext(ring);
	    } while (ring != ringFirst);
	}
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

// XXX Move up.
CxmpInline void
CxpTreeMpCacheInvalidate(CxtTreeMpPs *aPs)
{
    CxmCheckPtr(aPs);

    // Reset this node's parent pointer, to keep the old parent from using an
    // invalid cached value.
    aPs->parent = NULL;
}

// Calculate the partial score for aP, using aA and aB as children.  However, do
// some extra bookkeeping in order to be able to cache the results, and later
// recognize that precisely the same calculation was cached.
CxmpInline void
CxpTreeMpCachePscore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
    CxmError("XXX Not implemented");
}

// Unconditionally calculate the final score of a tree, using aA and aB as
// children.
CxmpInline unsigned
CxpTreeMpFscore(CxtTreeMpPs *aA, CxtTreeMpPs *aB, unsigned aMaxScore)
{
    CxmError("XXX Not implemented");
}

static CxtTreeMpPs *
CxpTreeMpScoreRecurse(CxtTreeMpData *aData, CxtRingObject *aRing,
		      CxtEdgeObject *aBisect)
{
    CxtTreeMpPs *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned degree;
    bool adjacent;
    CxtRingObject *ring;

    // Get the degree of the node.  Don't count the bisection edge (only an
    // issue if this node is adjacent to the bisection).
    degree = 1;
    adjacent = false;
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	if (CxRingEdge(ring) != aBisect)
	{
	    degree++;
	}
	else
	{
	    adjacent = true;
	}
    }

    switch (degree)
    {
	case 1:
	{
	    // Leaf node.  Do nothing.
	    rVal = (CxtTreeMpPs *) CxRingAuxGet(aRing, aData->ringAuxInd);
	    break;
	}
	case 2:
	{
	    // This is a trifurcating node that is adjacent to the bisection.
	    // Return the child node's ps, since this node's ps is
	    // irrelevant.
	    CxmAssert(adjacent);

	    // Clear the cache for the view that is being bypassed.  This is
	    // critical to correctness of the caching machinery, since each view
	    // should never be claimed as the parent of more than two other
	    // views.
	    CxpTreeMpCacheInvalidate((CxtTreeMpPs *)
				     CxRingAuxGet(aRing, aData->ringAuxInd));

	    // Get the ring element that connects to the other portion of the
	    // subtree on this side of the bisection.
	    for (ring = CxRingNext(aRing);
		 ring != aRing;
		 ring = CxRingNext(ring))
	    {
		if (CxRingEdge(ring) != aBisect)
		{
		    rVal = CxpTreeMpScoreRecurse(aData, CxRingOther(ring),
						 aBisect);
		    break;
		}
	    }
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		CxtTreeMpPs *psA, *psB;

		// This is a normal trifurcating node.  This is the common case,
		// and is handled separately from the code below for performance
		// reasons.

		// Recursively calculate partial scores for the subtrees.
		ring = CxRingNext(ring);
		psA = CxpTreeMpScoreRecurse(aData, ring, aBisect);

		ring = CxRingNext(ring);
		psB = CxpTreeMpScoreRecurse(aData, ring, aBisect);

		// Calculate the partial score for this node.
		rVal = (CxtTreeMpPs *) CxRingAuxGet(aRing, aData->ringAuxInd);
		CxpTreeMpCachePscore(rVal, psA, psB);
	    }
	    // Fall through if this node is adjacent to the bisection.
	}
	default:
	{
	    // This is a multifurcating node.
	    CxmError("XXX Not implemented");
	}
    }

    return rVal;
}

PyObject *
CxTreeMp(CxtTreeObject *self)
{
    PyObject *rVal;
    unsigned score;
    CxtNodeObject *base;
    CxtRingObject *ring;

    base = CxTreeBaseGet(self);
    if (base != NULL && (ring = CxNodeRing(base)) != NULL)
    {
	CxtTreeMpData *data;
	CxtEdgeObject *edge;
	CxtRingObject *ringA, *ringB;
	CxtTreeMpPs *psA, *psB;

	if (CxpTreeMpDataGet(self, &data))
	{
	    rVal = NULL;
	    goto RETURN;
	}

	edge = CxRingEdge(ring);
	CxEdgeRingsGet(edge, &ringA, &ringB);

	// Calculate partial scores for the subtrees on each end of edge.
	psA = CxpTreeMpScoreRecurse(data, ringA, NULL);
	psB = CxpTreeMpScoreRecurse(data, ringB, NULL);

	// Calculate the final score.
	score = CxpTreeMpFscore(psA, psB, UINT_MAX);
    }
    else
    {
	score = 0;
    }

    rVal = Py_BuildValue("i", score);
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrBestNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBestNeighborsMp(self->tr, maxHold);

	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrBetterNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBetterNeighborsMp(self->tr, maxHold);

	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrAllNeighborsMp(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrAllNeighborsMp(self->tr);

	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    return rVal;
}

PyObject *
CxTreeNheldGet(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	rVal = Py_BuildValue("i", CxTrNheldGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    return rVal;
}

PyObject *
CxTreeHeldGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int held;
    uint32_t neighbor, score, bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "i", &held) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (held >= CxTrNheldGet(self->tr))
    {
	Py_INCREF(PyExc_ValueError);
	rVal = PyExc_ValueError;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrHeldGet(self->tr, held, &neighbor, &score);
//	CxTrTbrNeighborGet(self->tr, neighbor,
//			    &bisect, &reconnectA, &reconnectB);

	rVal = Py_BuildValue("(i(iii))", score,
			       bisect, reconnectA, reconnectB);
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}
