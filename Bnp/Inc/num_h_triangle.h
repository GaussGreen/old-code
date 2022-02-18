/*****************************************************************************/
/*                                                                           */
/*  Include file for programs that call Triangle.                            */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#ifndef NUM_H_TRIANGLE_H
#define NUM_H_TRIANGLE_H

/*****************************************************************************/
/* structure pour stocker la triangulation */

typedef struct
{
    double* pointlist;               /* In / out */
    double* pointattributelist;      /* In / out */
    int*    pointmarkerlist;         /* In / out */
    int     numberofpoints;          /* In / out */
    int     numberofpointattributes; /* In / out */

    int*    trianglelist;               /* In / out */
    double* triangleattributelist;      /* In / out */
    double* triangleadoubleist;         /* In only */
    int*    neighborlist;               /* Out only */
    int     numberoftriangles;          /* In / out */
    int     numberofcorners;            /* In / out */
    int     numberoftriangleattributes; /* In / out */

    int* segmentlist;       /* In / out */
    int* segmentmarkerlist; /* In / out */
    int  numberofsegments;  /* In / out */

    double* holelist;      /* In / pointer to array copied out */
    int     numberofholes; /* In / copied out */

    double* regionlist;      /* In / pointer to array copied out */
    int     numberofregions; /* In / copied out */

    int*    edgelist;       /* Out only */
    int*    edgemarkerlist; /* Not used with Voronoi diagram; out only */
    double* normlist;       /* Used only with Voronoi diagram; out only */
    int     numberofedges;  /* Out only */

} SrtTriangulationIO, *SrtTriangulationIOPtr;

/*-----------------------------------------------------------------------------------------------------------------*/

/* THE FUNCTION THAT DOES THE FULL TRIANGULATION */
void triangulate_srt(char*, SrtTriangulationIO*);

/*-----------------------------------------------------------------------------------------------------------------*/

/* Labels that signify the result of TriPoint location.  The result of a        */
/*   search indicates that the TriPoint falls in the interior of a TriTriangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */

typedef enum
{
    INTRIANGLE,
    ONEDGE,
    ONVERTEX,
    OUTSIDE
} TriLocateResultType;

/*****************************************************************************/
/*                                                                           */
/*  The basic mesh data structures                                           */
/*                                                                           */
/*  There are three:  points, triangles, and shell edges (abbreviated        */
/*  `TriShell').  These three data structures, linked by pointers, comprise    */
/*  the mesh.  A TriPoint simply represents a TriPoint in space and its properties.*/
/*  A TriTriangle is a TriTriangle.  A shell edge is a special data structure used */
/*  to represent impenetrable segments in the mesh (including the outer      */
/*  boundary, boundaries of holes, and internal boundaries separating two    */
/*  triangulated regions).  Shell edges represent boundaries defined by the  */
/*  user that triangles may not lie across.                                  */
/*                                                                           */
/*  A TriTriangle consists of a list of three vertices, a list of three         */
/*  adjoining triangles, a list of three adjoining shell edges (when shell   */
/*  edges are used), an arbitrary number of optional user-defined floating-  */
/*  TriPoint attributes, and an optional area constraint.  The latter is an     */
/*  upper bound on the permissible area of each TriTriangle in a region, used   */
/*  for mesh refinement.                                                     */
/*                                                                           */
/*  For a TriTriangle on a boundary of the mesh, some or all of the neighboring */
/*  triangles may not be present.  For a TriTriangle in the interior of the     */
/*  mesh, often no neighboring shell edges are present.  Such absent         */
/*  triangles and shell edges are never represented by NULL pointers; they   */
/*  are represented by two special records:  `dummytri', the TriTriangle that   */
/*  fills "outer space", and `dummysh', the omnipresent shell edge.          */
/*  `dummytri' and `dummysh' are used for several reasons; for instance,     */
/*  they can be dereferenced and their contents examined without causing the */
/*  memory protection exception that would occur if NULL were dereferenced.  */
/*                                                                           */
/*  However, it is important to understand that a TriTriangle includes other    */
/*  information as well.  The pointers to adjoining vertices, triangles, and */
/*  shell edges are ordered in a way that indicates their geometric relation */
/*  to each other.  Furthermore, each of these pointers contains orientation */
/*  information.  Each pointer to an adjoining TriTriangle indicates which face */
/*  of that TriTriangle is contacted.  Similarly, each pointer to an adjoining  */
/*  shell edge indicates which side of that shell edge is contacted, and how */
/*  the shell edge is oriented relative to the TriTriangle.                     */
/*                                                                           */
/*  Shell edges are found abutting edges of triangles; either sandwiched     */
/*  between two triangles, or resting against one TriTriangle on an exterior    */
/*  boundary or hole boundary.                                               */
/*                                                                           */
/*  A shell edge consists of a list of two vertices, a list of two           */
/*  adjoining shell edges, and a list of two adjoining triangles.  One of    */
/*  the two adjoining triangles may not be present (though there should      */
/*  always be one), and neighboring shell edges might not be present.        */
/*  Shell edges also store a user-defined integer "boundary marker".         */
/*  Typically, this integer is used to indicate what sort of boundary        */
/*  conditions are to be applied at that location in a finite element        */
/*  simulation.                                                              */
/*                                                                           */
/*  Like triangles, shell edges maintain information about the relative      */
/*  orientation of neighboring objects.                                      */
/*                                                                           */
/*  Points are relatively simple.  A TriPoint is a list of floating TriPoint       */
/*  numbers, starting with the x, and y coordinates, followed by an          */
/*  arbitrary number of optional user-defined floating-TriPoint attributes,     */
/*  followed by an integer boundary marker.  During the segment insertion    */
/*  phase, there is also a pointer from each TriPoint to a TriTriangle that may    */
/*  contain it.  Each pointer is not always correct, but when one is, it     */
/*  speeds up segment insertion.  These pointers are assigned values once    */
/*  at the beginning of the segment insertion phase, and are not used or     */
/*  updated at any other time.  Edge swapping during segment insertion will  */
/*  render some of them incorrect.  Hence, don't rely upon them for          */
/*  anything.  For the most part, points do not have any information about   */
/*  what triangles or shell edges they are linked to.                        */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Handles                                                                  */
/*                                                                           */
/*  The oriented TriTriangle (`TriOrientedTriangle') and oriented shell edge (`edge') data  */
/*  structures defined below do not themselves store any part of the mesh.   */
/*  The mesh itself is made of `TriTriangle's, `TriShell's, and `TriPoint's.         */
/*                                                                           */
/*  Oriented triangles and oriented shell edges will usually be referred to  */
/*  as "handles".  A handle is essentially a pointer into the mesh; it       */
/*  allows you to "hold" one particular part of the mesh.  Handles are used  */
/*  to specify the regions in which one is traversing and modifying the mesh.*/
/*  A single `TriTriangle' may be held by many handles, or none at all.  (The   */
/*  latter case is not a memory leak, because the TriTriangle is still          */
/*  connected to other triangles in the mesh.)                               */
/*                                                                           */
/*  A `TriOrientedTriangle' is a handle that holds a TriTriangle.  It holds a specific side */
/*  of the TriTriangle.  An `edge' is a handle that holds a shell edge.  It     */
/*  holds either the left or right side of the edge.                         */
/*                                                                           */
/*  Navigation about the mesh is accomplished through a set of mesh          */
/*  manipulation primitives, further below.  Many of these primitives take   */
/*  a handle and produce a new handle that holds the mesh near the first     */
/*  handle.  Other primitives take two handles and glue the corresponding    */
/*  parts of the mesh together.  The exact position of the handles is        */
/*  important.  For instance, when two triangles are glued together by the   */
/*  bond() primitive, they are glued by the sides on which the handles lie.  */
/*                                                                           */
/*  Because points have no information about which triangles they are        */
/*  attached to, I commonly represent a TriPoint by use of a handle whose       */
/*  origin is the TriPoint.  A single handle can simultaneously represent a     */
/*  TriTriangle, an edge, and a TriPoint.                                          */
/*                                                                           */
/*****************************************************************************/

/* The TriTriangle data structure.  Each TriTriangle contains three pointers to    */
/*   adjoining triangles, plus three pointers to vertex points, plus three   */
/*   pointers to shell edges (defined below; these pointers are usually      */
/*   `dummysh').  It may or may not also contain user-defined attributes     */
/*   and/or a floating-TriPoint "area constraint".  It may also contain extra   */
/*   pointers for nodes, when the user asks for high-order elements.         */
/*   Because the size and structure of a `TriTriangle' is not decided until     */
/*   runtime, I haven't simply defined the type `TriTriangle' to be a struct.   */

typedef double** TriTriangle; /* doublely:  typedef TriTriangle *TriTriangle   */

/* An oriented TriTriangle:  includes a pointer to a TriTriangle and orientation.  */
/*   The orientation denotes an edge of the TriTriangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge is always        */
/*   directed to TriPoint counterclockwise about the corresponding TriTriangle.    */

typedef struct
{
    TriTriangle* tri;
    int          orient; /* Ranges from 0 to 2. */
} TriOrientedTriangle;

/* The shell data structure.  Each shell edge contains two pointers to       */
/*   adjoining shell edges, plus two pointers to vertex points, plus two     */
/*   pointers to adjoining triangles, plus one shell marker.                 */

typedef double** TriShell; /* doublely:  typedef TriShell *TriShell   */

/* An oriented shell edge:  includes a pointer to a shell edge and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */

typedef struct
{
    TriShell* sh;
    int       shorient; /* Ranges from 0 to 1. */
} TriOrientedShell;

/* The TriPoint data structure.  Each TriPoint is actually an array of doubles.      */
/*   The number of doubles is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a TriTriangle, is appended after the    */
/*   doubles.                                                                  */

typedef double* TriPoint;

/*****************************************************************************/
/*                                                                           */
/*  Mesh manipulation primitives.  Each TriTriangle contains three pointers to  */
/*  other triangles, with orientations.  Each pointer points not to the      */
/*  first byte of a TriTriangle, but to one of the first three bytes of a       */
/*  TriTriangle.  It is necessary to extract both the TriTriangle itself and the   */
/*  orientation.  To save memory, I keep both pieces of information in one   */
/*  pointer.  To make this possible, I assume that all triangles are aligned */
/*  to four-byte boundaries.  The `decode' routine below decodes a pointer,  */
/*  extracting an orientation (in the range 0 to 2) and a pointer to the     */
/*  beginning of a TriTriangle.  The `encode' routine compresses a pointer to a */
/*  TriTriangle and an orientation into a single pointer.  My assumptions that  */
/*  triangles are four-byte-aligned and that the `unsigned long' type is     */
/*  long enough to hold a pointer are two of the few kludges in this program.*/
/*                                                                           */
/*  Shell edges are manipulated similarly.  A pointer to a shell edge        */
/*  carries both an address and an orientation in the range 0 to 1.          */
/*                                                                           */
/*  The other primitives take an oriented TriTriangle or oriented shell edge,   */
/*  and return an oriented TriTriangle or oriented shell edge or TriPoint; or they */
/*  change the connections in the data structure.                            */
/*                                                                           */
/*****************************************************************************/

/********* Mesh manipulation primitives begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/********* Primitives for triangles                                  *********/
/*                                                                           */

/* decode() converts a pointer to an oriented TriTriangle.  The orientation is  */
/*   extracted from the two least significant bits of the pointer.           */

#define decode(ptr, TriOrientedTriangle)                                            \
    (TriOrientedTriangle).orient = (int)((unsigned long)(ptr) & (unsigned long)3l); \
    (TriOrientedTriangle).tri =                                                     \
        (TriTriangle*)((unsigned long)(ptr) ^ (unsigned long)(TriOrientedTriangle).orient)

/* encode() compresses an oriented TriTriangle into a single pointer.  It       */
/*   relies on the assumption that all triangles are aligned to four-byte    */
/*   boundaries, so the two least significant bits of (TriOrientedTriangle).tri are zero.*/

#define encode(TriOrientedTriangle) \
    (TriTriangle)(                  \
        (unsigned long)(TriOrientedTriangle).tri | (unsigned long)(TriOrientedTriangle).orient)

/* sym() finds the abutting TriTriangle, on the same edge.  Note that the       */
/*   edge direction is necessarily reversed, because TriTriangle/edge handles   */
/*   are always directed counterclockwise around the TriTriangle.               */

#define sym(TriOrientedTriangle1, TriOrientedTriangle2)              \
    ptr = (TriOrientedTriangle1).tri[(TriOrientedTriangle1).orient]; \
    decode(ptr, TriOrientedTriangle2);

#define symself(TriOrientedTriangle)                               \
    ptr = (TriOrientedTriangle).tri[(TriOrientedTriangle).orient]; \
    decode(ptr, TriOrientedTriangle);

/* These primitives determine or set the origin, destination, or apex of a   */
/* TriTriangle.                                                                 */

#define org(TriOrientedTriangle, pointptr) \
    pointptr = (TriPoint)(TriOrientedTriangle).tri[plus1mod3[(TriOrientedTriangle).orient] + 3]

#define dest(TriOrientedTriangle, pointptr) \
    pointptr = (TriPoint)(TriOrientedTriangle).tri[minus1mod3[(TriOrientedTriangle).orient] + 3]

#define apex(TriOrientedTriangle, pointptr) \
    pointptr = (TriPoint)(TriOrientedTriangle).tri[(TriOrientedTriangle).orient + 3]

/**                                                                         **/
/********* Mesh manipulation primitives end here                     *********/
/**                                                                         **/

/*  locate()   Find a TriTriangle or edge containing a given TriPoint.             */

TriLocateResultType locate(TriPoint searchpoint, TriOrientedTriangle* searchtri);

void triangledeinit(void);

TriTriangle* getdummytri();

#endif

/*-------------- End of file
 * -----------------------------------------------------------------------------------*/
