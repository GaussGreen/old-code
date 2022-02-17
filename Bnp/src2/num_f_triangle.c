/*****************************************************************************/
/*                                                                           */
/*  A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.      */
/*  (TriTriangle.c)                                                             */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

/* For single precision (which will save some memory and reduce paging),     */
/*   define the symbol SINGLE by using the -DSINGLE compiler switch or by    */
/*   writing "#define SINGLE" below.                                         */
/*                                                                           */
/* For double precision (which will allow you to refine meshes to a smaller  */
/*   edge length), leave SINGLE undefined.                                   */
/*                                                                           */
/* Double precision uses more memory, but improves the resolution of the     */
/*   meshes you can generate with Triangle.  It also reduces the likelihood  */
/*   of a floating exception due to overflow.  Finally, it is much faster    */
/*   than single precision on 64-bit architectures like the DEC Alpha.  I    */
/*   recommend double precision unless you want to generate a mesh for which */
/*   you do not have enough memory.                                          */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

/* If yours is not a Unix system, define the NO_TIMER compiler switch to     */
/*   remove the Unix-specific timing code.                                   */

#define NO_TIMER 

/* To insert lots of self-checks for internal errors, define the SELF_CHECK  */
/*   symbol.  This will slow down the program significantly.  It is best to  */
/*   define the symbol using the -DSELF_CHECK compiler switch, but you could */
/*   write "#define SELF_CHECK" below.  If you are modifying this code, I    */
/*   recommend you turn self-checks on.                                      */

/* #define SELF_CHECK */


/* On some machines, the exact arithmetic routines might be defeated by the  */
/*   use of internal extended precision floating-TriPoint registers.  Sometimes */
/*   this problem can be fixed by defining certain values to be volatile,    */
/*   thus forcing them to be stored to memory and rounded off.  This isn't   */
/*   a great solution, though, as it slows Triangle down.                    */
/*                                                                           */
/* To try this out, write "#define INEXACT volatile" below.  Normally,       */
/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */

#define INEXACT /* Nothing */
/* #define INEXACT volatile */

/* Maximum number of characters in a file name (including the null).         */

#define FILENAMESIZE 512

/* Maximum number of characters in a line read from a file (including the    */
/*   null).                                                                  */

#define INPUTLINESIZE 512

/* For efficiency, a variety of data structures are allocated in bulk.  The  */
/*   following constants determine how many of each structure is allocated   */
/*   at once.                                                                */

#define TRIPERBLOCK 4092           /* Number of triangles allocated at once. */
#define SHELLEPERBLOCK 508       /* Number of shell edges allocated at once. */
#define POINTPERBLOCK 4092            /* Number of points allocated at once. */
#define VIRUSPERBLOCK 1020   /* Number of virus triangles allocated at once. */
/* Number of encroached segments allocated at once. */
#define BADSEGMENTPERBLOCK 252
/* Number of skinny triangles allocated at once. */
#define BADTRIPERBLOCK 4092
/* Number of splay tree nodes allocated at once. */
#define SPLAYNODEPERBLOCK 508

/* The TriPoint marker DEADPOINT is an arbitrary number chosen large enough to  */
/*   (hopefully) not conflict with user boundary markers.  Make sure that it */
/*   is small enough to fit into your machine's integer size.                */

#define DEADPOINT -1073741824

/* The next line is used to outsmart some very stupid compilers.  If your    */
/*   compiler is smarter, feel free to replace the "int" with "void".        */
/*   Not that it matters.                                                    */

#define VOID int

/* Two constants for algorithms based on random sampling.  Both constants    */
/*   have been chosen empirically to optimize their respective algorithms.   */

/* Used for the TriPoint location scheme of Mucke, Saias, and Zhu, to decide    */
/*   how large a random sample of triangles to inspect.                      */
#define SAMPLEFACTOR 11
/* Used in Fortune's sweepline Delaunay algorithm to determine what fraction */
/*   of boundary edges should be maintained in the splay tree for TriPoint      */
/*   location on the front.                                                  */
#define SAMPLERATE 10

/* A number that speaks for itself, every kissable digit.                    */

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

/* Another fave.                                                             */

#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732

/* And here's one for those of you who are intimidated by math.              */

#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333

#include "stdio.h"
#include "string.h"
#include "math.h"
#ifndef NO_TIMER
#include "time.h"
#endif /* NO_TIMER */

#include "num_h_triangle.h"

/* The following obscenity seems to be necessary to ensure that this program */
/* will port to Dec Alphas running OSF/1, because their stdio.h file commits */
/* the unpardonable sin of including stdlib.h.  Hence, malloc(), free(), and */
/* exit() may or may not already be defined at this TriPoint.  I declare these  */
/* functions explicitly because some non-ANSI C compilers lack stdlib.h.     */

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
extern void exit();
extern double strtod();
extern long strtol();
#endif /* _STDLIB_H_ */

/* A few forward declarations.                                               */

void poolrestart();

/* Labels that signify whether a record consists primarily of pointers or of */
/*   floating-TriPoint words.  Used to make decisions about data alignment.     */

enum wordtype {POINTER, FLOATINGPOINT};


/* Labels that signify the result of site insertion.  The result indicates   */
/*   that the TriPoint was inserted with complete success, was inserted but     */
/*   encroaches on a segment, was not inserted because it lies on a segment, */
/*   or was not inserted because another TriPoint occupies the same location.   */

enum insertsiteresult {SUCCESSFULPOINT, ENCROACHINGPOINT, VIOLATINGPOINT,
                       DUPLICATEPOINT};

/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction TriTriangle, along the left edge of the direction TriTriangle,  */
/*   or along the right edge of the direction TriTriangle.                      */

enum finddirectionresult {WITHIN, LEFTCOLLINEAR, RIGHTCOLLINEAR};

/* Labels that signify the result of the circumcenter computation routine.   */
/*   The return value indicates which edge of the TriTriangle is shortest.      */

enum circumcenterresult {OPPOSITEORG, OPPOSITEDEST, OPPOSITEAPEX};


/* A queue used to store encroached segments.  Each segment's vertices are   */
/*   stored so that one can check whether a segment is still the same.       */

struct badsegment {
  TriOrientedShell encsegment;                          /* An encroached segment. */
  TriPoint segorg, segdest;                                /* The two vertices. */
  struct badsegment *nextsegment;     /* Pointer to next encroached segment. */
};

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the TriTriangle.  Each TriTriangle's vertices are    */
/*   stored so that one can check whether a TriTriangle is still the same.      */

struct badface {
  TriOrientedTriangle badfacetri;                              /* A bad TriTriangle. */
  REAL key;                             /* cos^2 of smallest (apical) angle. */
  TriPoint faceorg, facedest, faceapex;                  /* The three vertices. */
  struct badface *nextface;                 /* Pointer to next bad TriTriangle. */
};

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not TriPoint directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `TriPoint' or to a TriTriangle (in encoded format, so that */
/*   an orientation is included).  In the latter case, the origin of the     */
/*   oriented TriTriangle is the apex of a "circle event" of the sweepline      */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */

struct event {
  REAL xkey, ykey;                              /* Coordinates of the event. */
  VOID *eventptr;       /* Can be a TriPoint or the location of a circle event. */
  int heapposition;              /* Marks this event's position in the heap. */
};

/* A node in the splay tree.  Each node holds an oriented ghost TriTriangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a TriTriangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a TriTriangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */

struct splaynode {
  TriOrientedTriangle keyedge;                  /* Lprev of an edge on the front. */
  TriPoint keydest;            /* Used to verify that splay node is still live. */
  struct splaynode *lchild, *rchild;              /* Children in splay tree. */
};

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* itemwordtype is set to POINTER or FLOATINGPOINT, and is used to suggest   */
/*   what sort of word the record is primarily made up of.  alignbytes       */
/*   determines how new records should be aligned in memory.  itembytes and  */
/*   itemwords are the length of a record in bytes (after rounding up) and   */
/*   words.  itemsperblock is the number of items allocated at once in a     */
/*   single block.  items is the number of currently allocated items.        */
/*   maxitems is the maximum number of items that have been allocated at     */
/*   once; it is the current number of items plus the number of records kept */
/*   on deaditemstack.                                                       */

struct memorypool {
  VOID **firstblock, **nowblock;
  VOID *nextitem;
  VOID *deaditemstack;
  VOID **pathblock;
  VOID *pathitem;
  enum wordtype itemwordtype;
  int alignbytes;
  int itembytes, itemwords;
  int itemsperblock;
  long items, maxitems;
  int unallocateditems;
  int pathitemsleft;
};

/* Variables used to allocate memory for triangles, shell edges, points,     */
/*   viri (triangles being eaten), bad (encroached) segments, bad (skinny    */
/*   or too large) triangles, and splay tree nodes.                          */

struct memorypool triangles;
struct memorypool shelles;
struct memorypool points;
struct memorypool viri;
struct memorypool badsegments;
struct memorypool badtriangles;
struct memorypool splaynodes;

/* Variables that maintain the bad TriTriangle queues.  The tails are pointers  */
/*   to the pointers that have to be filled in to enqueue an item.           */

struct badface *queuefront[64];
struct badface **queuetail[64];

static REAL xmin, xmax, ymin, ymax;                              /* x and y bounds. */
static REAL xminextreme;        /* Nonexistent x value used as a flag in sweepline. */
static int inpoints;                                     /* Number of input points. */
static int inelements;                                /* Number of input triangles. */
static int insegments;                                 /* Number of input segments. */
static int holes;                                         /* Number of input holes. */
static int regions;                                     /* Number of input regions. */
static long edges;                                       /* Number of output edges. */
static int mesh_dim;                                  /* Dimension (ought to be 2). */
static int nextras;                              /* Number of attributes per TriPoint. */
static int eextras;                           /* Number of attributes per TriTriangle. */
static long hullsize;                            /* Number of edges of convex hull. */
static int triwords;                                   /* Total words per TriTriangle. */
static int shwords;                                  /* Total words per shell edge. */
static int pointmarkindex;             /* Index to find boundary marker of a TriPoint. */
static int point2triindex;         /* Index to find a TriTriangle adjacent to a TriPoint. */
static int highorderindex;    /* Index to find extra nodes for high-order elements. */
static int elemattribindex;              /* Index to find attributes of a TriTriangle. */
static int areaboundindex;               /* Index to find area bound of a TriTriangle. */
static int checksegments;           /* Are there segments in the triangulation yet? */
static int readnodefile;                             /* Has a .node file been read? */
static long samples;                /* Number of random samples for TriPoint location. */
static unsigned long randomseed;                     /* Current random number seed. */

static REAL splitter;       /* Used to split REAL factors for exact multiplication. */
static REAL epsilon;                             /* Floating-TriPoint machine epsilon. */
static REAL resulterrbound;
static REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
static REAL iccerrboundA, iccerrboundB, iccerrboundC;

static long incirclecount;                   /* Number of incircle tests performed. */
static long counterclockcount;       /* Number of counterclockwise tests performed. */
static long hyperbolacount;        /* Number of right-of-hyperbola tests performed. */
static long circumcentercount;    /* Number of circumcenter calculations performed. */
static long circletopcount;         /* Number of circle top calculations performed. */

/* Switches for the triangulator.                                            */
/*   poly: -p switch.  refine: -r switch.                                    */
/*   quality: -q switch.                                                     */
/*     minangle: minimum angle bound, specified after -q switch.             */
/*     goodangle: cosine squared of minangle.                                */
/*   vararea: -a switch without number.                                      */
/*   fixedarea: -a switch with number.                                       */
/*     maxarea: maximum area bound, specified after -a switch.               */
/*   regionattrib: -A switch.  convex: -c switch.                            */
/*   firstnumber: inverse of -z switch.  All items are numbered starting     */
/*     from firstnumber.                                                     */
/*   edgesout: -e switch.  voronoi: -v switch.                               */
/*   neighbors: -n switch.  geomview: -g switch.                             */
/*   nobound: -B switch.  nopolywritten: -P switch.                          */
/*   nonodewritten: -N switch.  noelewritten: -E switch.                     */
/*   noiterationnum: -I switch.  noholes: -O switch.                         */
/*   noexact: -X switch.                                                     */
/*   order: element order, specified after -o switch.                        */
/*   nobisect: count of how often -Y switch is selected.                     */
/*   steiner: maximum number of Steiner points, specified after -S switch.   */
/*     steinerleft: number of Steiner points not yet used.                   */
/*   incremental: -i switch.  sweepline: -F switch.                          */
/*   dwyer: inverse of -l switch.                                            */
/*   splitseg: -s switch.                                                    */
/*   docheck: -C switch.                                                     */
/*   quiet: -Q switch.  verbose: count of how often -V switch is selected.   */
/*   useshelles: -p, -r, -q, or -c switch; determines whether shell edges    */
/*     are used at all.                                                      */
/*                                                                           */
/* Read the instructions to find out the meaning of these switches.          */

static int poly, refine, quality, vararea, fixedarea, regionattrib, convex;
static int firstnumber;
static int edgesout, voronoi, neighbors, geomview;
static int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
static int noholes, noexact;
static int incremental, sweepline, dwyer;
static int splitseg;
static int docheck;
static int quiet, verbose;
static int useshelles;
static int order;
static int nobisect;
static int steiner, steinerleft;
static REAL minangle, goodangle;
static REAL maxarea;


/* Triangular bounding box points.                                           */

static TriPoint infpoint1, infpoint2, infpoint3;

/* Keep base address so we can free() it later. */
static TriTriangle *dummytribase;      

/* Pointer to the omnipresent shell edge.  Referenced by any TriTriangle or     */
/*   shell edge that isn't really connected to a shell edge at that          */
/*   location.                                                               */
							  
static TriShell *dummysh;
static TriShell *dummyshbase;         /* Keep base address so we can free() it later. */

/* Pointer to a recently visited TriTriangle.  Improves TriPoint location if       */
/*   proximate points are inserted sequentially.                             */

static TriOrientedTriangle recenttri;


/* Pointer to the `TriTriangle' that occupies all of "outer space".             */

static TriTriangle *dummytri;



/********* Mesh manipulation primitives begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/********* Primitives for triangles                                  *********/
/*                                                                           */
/* The following edge manipulation primitives are all described by Guibas    */
/*   and Stolfi.  However, they use an edge-based data structure, whereas I  */
/*   am using a TriTriangle-based data structure.                               */

/* lnext() finds the next edge (counterclockwise) of a TriTriangle.             */

#define lnext(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  (TriOrientedTriangle2).tri = (TriOrientedTriangle1).tri;                                            \
  (TriOrientedTriangle2).orient = plus1mod3[(TriOrientedTriangle1).orient]

#define lnextself(TriOrientedTriangle)                                                    \
  (TriOrientedTriangle).orient = plus1mod3[(TriOrientedTriangle).orient]

/* lprev() finds the previous edge (clockwise) of a TriTriangle.                */

#define lprev(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  (TriOrientedTriangle2).tri = (TriOrientedTriangle1).tri;                                            \
  (TriOrientedTriangle2).orient = minus1mod3[(TriOrientedTriangle1).orient]

#define lprevself(TriOrientedTriangle)                                                    \
  (TriOrientedTriangle).orient = minus1mod3[(TriOrientedTriangle).orient]

/* onext() spins counterclockwise around a TriPoint; that is, it finds the next */
/*   edge with the same origin in the counterclockwise direction.  This edge */
/*   will be part of a different TriTriangle.                                   */

#define onext(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  lprev(TriOrientedTriangle1, TriOrientedTriangle2);                                                  \
  symself(TriOrientedTriangle2);

#define onextself(TriOrientedTriangle)                                                    \
  lprevself(TriOrientedTriangle);                                                         \
  symself(TriOrientedTriangle);

/* oprev() spins clockwise around a TriPoint; that is, it finds the next edge   */
/*   with the same origin in the clockwise direction.  This edge will be     */
/*   part of a different TriTriangle.                                           */

#define oprev(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  sym(TriOrientedTriangle1, TriOrientedTriangle2);                                                    \
  lnextself(TriOrientedTriangle2);

#define oprevself(TriOrientedTriangle)                                                    \
  symself(TriOrientedTriangle);                                                           \
  lnextself(TriOrientedTriangle);

/* dnext() spins counterclockwise around a TriPoint; that is, it finds the next */
/*   edge with the same destination in the counterclockwise direction.  This */
/*   edge will be part of a different TriTriangle.                              */

#define dnext(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  sym(TriOrientedTriangle1, TriOrientedTriangle2);                                                    \
  lprevself(TriOrientedTriangle2);

#define dnextself(TriOrientedTriangle)                                                    \
  symself(TriOrientedTriangle);                                                           \
  lprevself(TriOrientedTriangle);

/* dprev() spins clockwise around a TriPoint; that is, it finds the next edge   */
/*   with the same destination in the clockwise direction.  This edge will   */
/*   be part of a different TriTriangle.                                        */

#define dprev(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  lnext(TriOrientedTriangle1, TriOrientedTriangle2);                                                  \
  symself(TriOrientedTriangle2);

#define dprevself(TriOrientedTriangle)                                                    \
  lnextself(TriOrientedTriangle);                                                         \
  symself(TriOrientedTriangle);

/* rnext() moves one edge counterclockwise about the adjacent TriTriangle.      */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rnext(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  sym(TriOrientedTriangle1, TriOrientedTriangle2);                                                    \
  lnextself(TriOrientedTriangle2);                                                        \
  symself(TriOrientedTriangle2);

#define rnextself(TriOrientedTriangle)                                                    \
  symself(TriOrientedTriangle);                                                           \
  lnextself(TriOrientedTriangle);                                                         \
  symself(TriOrientedTriangle);

/* rnext() moves one edge clockwise about the adjacent TriTriangle.             */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rprev(TriOrientedTriangle1, TriOrientedTriangle2)                                             \
  sym(TriOrientedTriangle1, TriOrientedTriangle2);                                                    \
  lprevself(TriOrientedTriangle2);                                                        \
  symself(TriOrientedTriangle2);

#define rprevself(TriOrientedTriangle)                                                    \
  symself(TriOrientedTriangle);                                                           \
  lprevself(TriOrientedTriangle);                                                         \
  symself(TriOrientedTriangle);

#define setorg(TriOrientedTriangle, pointptr)                                             \
  (TriOrientedTriangle).tri[plus1mod3[(TriOrientedTriangle).orient] + 3] = (TriTriangle) pointptr

#define setdest(TriOrientedTriangle, pointptr)                                            \
  (TriOrientedTriangle).tri[minus1mod3[(TriOrientedTriangle).orient] + 3] = (TriTriangle) pointptr

#define setapex(TriOrientedTriangle, pointptr)                                            \
  (TriOrientedTriangle).tri[(TriOrientedTriangle).orient + 3] = (TriTriangle) pointptr

#define setvertices2null(TriOrientedTriangle)                                             \
  (TriOrientedTriangle).tri[3] = (TriTriangle) NULL;                                         \
  (TriOrientedTriangle).tri[4] = (TriTriangle) NULL;                                         \
  (TriOrientedTriangle).tri[5] = (TriTriangle) NULL;

/* Bond two triangles together.                                              */

#define bond(TriOrientedTriangle1, TriOrientedTriangle2)                                              \
  (TriOrientedTriangle1).tri[(TriOrientedTriangle1).orient] = encode(TriOrientedTriangle2);                       \
  (TriOrientedTriangle2).tri[(TriOrientedTriangle2).orient] = encode(TriOrientedTriangle1)

/* Dissolve a bond (from one side).  Note that the other TriTriangle will still */
/*   think it's connected to this TriTriangle.  Usually, however, the other     */
/*   TriTriangle is being deleted entirely, or bonded to another TriTriangle, so   */
/*   it doesn't matter.                                                      */

#define dissolve(TriOrientedTriangle)                                                     \
  (TriOrientedTriangle).tri[(TriOrientedTriangle).orient] = (TriTriangle) dummytri

/* Copy a TriTriangle/edge handle.                                              */

#define TriOrientedTrianglecopy(TriOrientedTriangle1, TriOrientedTriangle2)                                       \
  (TriOrientedTriangle2).tri = (TriOrientedTriangle1).tri;                                            \
  (TriOrientedTriangle2).orient = (TriOrientedTriangle1).orient

/* Test for equality of TriTriangle/edge handles.                               */

#define TriOrientedTriangleequal(TriOrientedTriangle1, TriOrientedTriangle2)                                      \
  (((TriOrientedTriangle1).tri == (TriOrientedTriangle2).tri) &&                                      \
   ((TriOrientedTriangle1).orient == (TriOrientedTriangle2).orient))

/* Primitives to infect or cure a TriTriangle with the virus.  These rely on    */
/*   the assumption that all shell edges are aligned to four-byte boundaries.*/

#define infect(TriOrientedTriangle)                                                       \
  (TriOrientedTriangle).tri[6] = (TriTriangle)                                               \
                     ((unsigned long) (TriOrientedTriangle).tri[6] | (unsigned long) 2l)

#define uninfect(TriOrientedTriangle)                                                     \
  (TriOrientedTriangle).tri[6] = (TriTriangle)                                               \
                     ((unsigned long) (TriOrientedTriangle).tri[6] & ~ (unsigned long) 2l)

/* Test a TriTriangle for viral infection.                                      */

#define infected(TriOrientedTriangle)                                                     \
  (((unsigned long) (TriOrientedTriangle).tri[6] & (unsigned long) 2l) != 0)

/* Check or set a TriTriangle's attributes.                                     */

#define elemattribute(TriOrientedTriangle, attnum)                                        \
  ((REAL *) (TriOrientedTriangle).tri)[elemattribindex + (attnum)]

#define setelemattribute(TriOrientedTriangle, attnum, value)                              \
  ((REAL *) (TriOrientedTriangle).tri)[elemattribindex + (attnum)] = value

/* Check or set a TriTriangle's maximum area bound.                             */

#define areabound(TriOrientedTriangle)  ((REAL *) (TriOrientedTriangle).tri)[areaboundindex]

#define setareabound(TriOrientedTriangle, value)                                          \
  ((REAL *) (TriOrientedTriangle).tri)[areaboundindex] = value

/********* Primitives for shell edges                                *********/
/*                                                                           */
/*                                                                           */

/* sdecode() converts a pointer to an oriented shell edge.  The orientation  */
/*   is extracted from the least significant bit of the pointer.  The two    */
/*   least significant bits (one for orientation, one for viral infection)   */
/*   are masked out to produce the real pointer.                             */

#define sdecode(sptr, edge)                                                   \
  (edge).shorient = (int) ((unsigned long) (sptr) & (unsigned long) 1l);      \
  (edge).sh = (TriShell *)                                                      \
              ((unsigned long) (sptr) & ~ (unsigned long) 3l)

/* sencode() compresses an oriented shell edge into a single pointer.  It    */
/*   relies on the assumption that all shell edges are aligned to two-byte   */
/*   boundaries, so the least significant bit of (edge).sh is zero.          */

#define sencode(edge)                                                         \
  (TriShell) ((unsigned long) (edge).sh | (unsigned long) (edge).shorient)

/* ssym() toggles the orientation of a shell edge.                           */

#define ssym(edge1, edge2)                                                    \
  (edge2).sh = (edge1).sh;                                                    \
  (edge2).shorient = 1 - (edge1).shorient

#define ssymself(edge)                                                        \
  (edge).shorient = 1 - (edge).shorient

/* spivot() finds the other shell edge (from the same segment) that shares   */
/*   the same origin.                                                        */

#define spivot(edge1, edge2)                                                  \
  sptr = (edge1).sh[(edge1).shorient];                                        \
  sdecode(sptr, edge2)

#define spivotself(edge)                                                      \
  sptr = (edge).sh[(edge).shorient];                                          \
  sdecode(sptr, edge)

/* snext() finds the next shell edge (from the same segment) in sequence;    */
/*   one whose origin is the input shell edge's destination.                 */

#define snext(edge1, edge2)                                                   \
  sptr = (edge1).sh[1 - (edge1).shorient];                                    \
  sdecode(sptr, edge2)

#define snextself(edge)                                                       \
  sptr = (edge).sh[1 - (edge).shorient];                                      \
  sdecode(sptr, edge)

/* These primitives determine or set the origin or destination of a shell    */
/*   edge.                                                                   */

#define sorg(edge, pointptr)                                                  \
  pointptr = (TriPoint) (edge).sh[2 + (edge).shorient]

#define sdest(edge, pointptr)                                                 \
  pointptr = (TriPoint) (edge).sh[3 - (edge).shorient]

#define setsorg(edge, pointptr)                                               \
  (edge).sh[2 + (edge).shorient] = (TriShell) pointptr

#define setsdest(edge, pointptr)                                              \
  (edge).sh[3 - (edge).shorient] = (TriShell) pointptr

/* These primitives read or set a shell marker.  Shell markers are used to   */
/*   hold user boundary information.                                         */

#define mark(edge)  (* (int *) ((edge).sh + 6))

#define setmark(edge, value)                                                  \
  * (int *) ((edge).sh + 6) = value

/* Bond two shell edges together.                                            */

#define sbond(edge1, edge2)                                                   \
  (edge1).sh[(edge1).shorient] = sencode(edge2);                              \
  (edge2).sh[(edge2).shorient] = sencode(edge1)

/* Dissolve a shell edge bond (from one side).  Note that the other shell    */
/*   edge will still think it's connected to this shell edge.                */

#define sdissolve(edge)                                                       \
  (edge).sh[(edge).shorient] = (TriShell) dummysh

/* Copy a shell edge.                                                        */

#define shellecopy(edge1, edge2)                                              \
  (edge2).sh = (edge1).sh;                                                    \
  (edge2).shorient = (edge1).shorient

/* Test for equality of shell edges.                                         */

#define shelleequal(edge1, edge2)                                             \
  (((edge1).sh == (edge2).sh) &&                                              \
   ((edge1).shorient == (edge2).shorient))

/********* Primitives for interacting triangles and shell edges      *********/
/*                                                                           */
/*                                                                           */

/* tspivot() finds a shell edge abutting a TriTriangle.                         */

#define tspivot(TriOrientedTriangle, edge)                                                \
  sptr = (TriShell) (TriOrientedTriangle).tri[6 + (TriOrientedTriangle).orient];                        \
  sdecode(sptr, edge)

/* stpivot() finds a TriTriangle abutting a shell edge.  It requires that the   */
/*   variable `ptr' of type `TriTriangle' be defined.                           */

#define stpivot(edge, TriOrientedTriangle)                                                \
  ptr = (TriTriangle) (edge).sh[4 + (edge).shorient];                            \
  decode(ptr, TriOrientedTriangle)

/* Bond a TriTriangle to a shell edge.                                          */

#define tsbond(TriOrientedTriangle, edge)                                                 \
  (TriOrientedTriangle).tri[6 + (TriOrientedTriangle).orient] = (TriTriangle) sencode(edge);             \
  (edge).sh[4 + (edge).shorient] = (TriShell) encode(TriOrientedTriangle)

/* Dissolve a bond (from the TriTriangle side).                                 */

#define tsdissolve(TriOrientedTriangle)                                                   \
  (TriOrientedTriangle).tri[6 + (TriOrientedTriangle).orient] = (TriTriangle) dummysh

/* Dissolve a bond (from the shell edge side).                               */

#define stdissolve(edge)                                                      \
  (edge).sh[4 + (edge).shorient] = (TriShell) dummytri

/********* Primitives for points                                     *********/
/*                                                                           */
/*                                                                           */

#define pointmark(pt)  ((int *) (pt))[pointmarkindex]

#define setpointmark(pt, value)                                               \
  ((int *) (pt))[pointmarkindex] = value

#define point2tri(pt)  ((TriTriangle *) (pt))[point2triindex]

#define setpoint2tri(pt, value)                                               \
  ((TriTriangle *) (pt))[point2triindex] = value


/**                                                                         **/
/********* Mesh manipulation primitives end here                     *********/
/**                                                                         **/


/********* User interaction routines begin here                      *********/
/**                                                                         **/
/**                                                                         **/

/* Fast lookup arrays to speed some of the mesh manipulation primitives.     */

static int plus1mod3[3] = {1, 2, 0};
static int minus1mod3[3] = {2, 0, 1};

/*****************************************************************************/
/*                                                                           */
/*  internalerror()   Ask the user to send me the defective product.  Exit.  */
/*                                                                           */
/*****************************************************************************/

void internalerror()
{
  printf("  Please report this bug to jrs@cs.cmu.edu\n");
  printf("  Include the message above, your input data set, and the exact\n");
  printf("    command line you used to run Triangle.\n");
  exit(1);
}

/*****************************************************************************/
/*                                                                           */
/*  parsecommandline()   Read the command line, identify switches, and set   */
/*                       up options and file names.                          */
/*                                                                           */
/*  The effects of this routine are felt entirely through global variables.  */
/*                                                                           */
/*****************************************************************************/

void parsecommandline(int argc, char **argv)
{
#define STARTINDEX 0
  int i, j;

  poly = refine = quality = vararea = fixedarea = regionattrib = convex = 0;
  firstnumber = 1;
  edgesout = voronoi = neighbors = geomview = 0;
  nobound = nopolywritten = nonodewritten = noelewritten = noiterationnum = 0;
  noholes = noexact = 0;
  incremental = sweepline = 0;
  dwyer = 1;
  splitseg = 0;
  docheck = 0;
  nobisect = 0;
  steiner = -1;
  order = 1;
  minangle = 0.0;
  maxarea = -1.0;
  quiet = verbose = 0;

  for (i = STARTINDEX; i < argc; i++) {
      for (j = STARTINDEX; argv[i][j] != '\0'; j++) {
        if (argv[i][j] == 'p') {
          poly = 1;
	}
        if (argv[i][j] == 'A') {
          regionattrib = 1;
        }
        if (argv[i][j] == 'c') {
          convex = 1;
        }
        if (argv[i][j] == 'z') {
          firstnumber = 0;
        }
        if (argv[i][j] == 'e') {
          edgesout = 1;
	}
        if (argv[i][j] == 'v') {
          voronoi = 1;
	}
        if (argv[i][j] == 'n') {
          neighbors = 1;
	}
        if (argv[i][j] == 'g') {
          geomview = 1;
	}
        if (argv[i][j] == 'B') {
          nobound = 1;
	}
        if (argv[i][j] == 'P') {
          nopolywritten = 1;
	}
        if (argv[i][j] == 'N') {
          nonodewritten = 1;
	}
        if (argv[i][j] == 'E') {
          noelewritten = 1;
	}
        if (argv[i][j] == 'O') {
          noholes = 1;
	}
        if (argv[i][j] == 'X') {
          noexact = 1;
	}
        if (argv[i][j] == 'o') {
          if (argv[i][j + 1] == '2') {
            j++;
            order = 2;
          }
	}
        if (argv[i][j] == 'l') {
          dwyer = 0;
        }
        if (argv[i][j] == 'Q') {
          quiet = 1;
        }
        if (argv[i][j] == 'V') {
          verbose++;
        }
      }
  }
  steinerleft = steiner;
  useshelles = poly || refine || quality || convex;
  goodangle = cos(minangle * PI / 180.0);
  goodangle *= goodangle;
  if (refine && noiterationnum) {
    printf(
      "Error:  You cannot use the -I switch when refining a triangulation.\n");
    exit(1);
  }
  /* Be careful not to allocate space for element area constraints that */
  /*   will never be assigned any value (other than the default -1.0).  */
  if (!refine && !poly) {
    vararea = 0;
  }
  /* Be careful not to add an extra attribute to each element unless the */
  /*   input supports it (PSLG in, but not refining a preexisting mesh). */
  if (refine || !poly) {
    regionattrib = 0;
  }

}

/**                                                                         **/
/**                                                                         **/
/********* User interaction routines begin here                      *********/

/********* Debugging routines begin here                             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  printtriangle()   Print out the details of a TriTriangle/edge handle.       */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about a      */
/*  TriTriangle/edge handle in digestible form.  It's also used when the        */
/*  highest level of verbosity (`-VVV') is specified.                        */
/*                                                                           */
/*****************************************************************************/

void printtriangle(TriOrientedTriangle *t)
{
  TriOrientedTriangle printtri;
  TriOrientedShell printsh;
  TriPoint printpoint;

  printf("TriTriangle x%lx with orientation %d:\n", (unsigned long) t->tri,
         t->orient);
  decode(t->tri[0], printtri);
  if (printtri.tri == dummytri) {
    printf("    [0] = Outer space\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[1], printtri);
  if (printtri.tri == dummytri) {
    printf("    [1] = Outer space\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[2], printtri);
  if (printtri.tri == dummytri) {
    printf("    [2] = Outer space\n");
  } else {
    printf("    [2] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  org(*t, printpoint);
  if (printpoint == (TriPoint) NULL)
    printf("    Origin[%d] = NULL\n", (t->orient + 1) % 3 + 3);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 1) % 3 + 3, (unsigned long) printpoint,
           printpoint[0], printpoint[1]);
  dest(*t, printpoint);
  if (printpoint == (TriPoint) NULL)
    printf("    Dest  [%d] = NULL\n", (t->orient + 2) % 3 + 3);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 2) % 3 + 3, (unsigned long) printpoint,
           printpoint[0], printpoint[1]);
  apex(*t, printpoint);
  if (printpoint == (TriPoint) NULL)
    printf("    Apex  [%d] = NULL\n", t->orient + 3);
  else
    printf("    Apex  [%d] = x%lx  (%.12g, %.12g)\n",
           t->orient + 3, (unsigned long) printpoint,
           printpoint[0], printpoint[1]);
  if (useshelles) {
    sdecode(t->tri[6], printsh);
    if (printsh.sh != dummysh) {
      printf("    [6] = x%lx  %d\n", (unsigned long) printsh.sh,
             printsh.shorient);
    }
    sdecode(t->tri[7], printsh);
    if (printsh.sh != dummysh) {
      printf("    [7] = x%lx  %d\n", (unsigned long) printsh.sh,
             printsh.shorient);
    }
    sdecode(t->tri[8], printsh);
    if (printsh.sh != dummysh) {
      printf("    [8] = x%lx  %d\n", (unsigned long) printsh.sh,
             printsh.shorient);
    }
  }
  if (vararea) {
    printf("    Area constraint:  %.4g\n", areabound(*t));
  }
}

/*****************************************************************************/
/*                                                                           */
/*  printshelle()   Print out the details of a shell edge handle.            */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about a      */
/*  shell edge handle in digestible form.  It's also used when the highest   */
/*  level of verbosity (`-VVV') is specified.                                */
/*                                                                           */
/*****************************************************************************/

void printshelle(TriOrientedShell *s)
{
  TriOrientedShell printsh;
  TriOrientedTriangle printtri;
  TriPoint printpoint;

  printf("shell edge x%lx with orientation %d and mark %d:\n",
         (unsigned long) s->sh, s->shorient, mark(*s));
  sdecode(s->sh[0], printsh);
  if (printsh.sh == dummysh) {
    printf("    [0] = No shell\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printsh.sh,
           printsh.shorient);
  }
  sdecode(s->sh[1], printsh);
  if (printsh.sh == dummysh) {
    printf("    [1] = No shell\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printsh.sh,
           printsh.shorient);
  }
  sorg(*s, printpoint);
  if (printpoint == (TriPoint) NULL)
    printf("    Origin[%d] = NULL\n", 2 + s->shorient);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           2 + s->shorient, (unsigned long) printpoint,
           printpoint[0], printpoint[1]);
  sdest(*s, printpoint);
  if (printpoint == (TriPoint) NULL)
    printf("    Dest  [%d] = NULL\n", 3 - s->shorient);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           3 - s->shorient, (unsigned long) printpoint,
           printpoint[0], printpoint[1]);
  decode(s->sh[4], printtri);
  if (printtri.tri == dummytri) {
    printf("    [4] = Outer space\n");
  } else {
    printf("    [4] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(s->sh[5], printtri);
  if (printtri.tri == dummytri) {
    printf("    [5] = Outer space\n");
  } else {
    printf("    [5] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
}

/**                                                                         **/
/**                                                                         **/
/********* Debugging routines end here                               *********/

/********* Memory management routines begin here                     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  poolinit()   Initialize a pool of memory for allocation of items.        */
/*                                                                           */
/*  This routine initializes the machinery for allocating items.  A `pool'   */
/*  is created whose records have size at least `bytecount'.  Items will be  */
/*  allocated in `itemcount'-item blocks.  Each item is assumed to be a      */
/*  collection of words, and either pointers or floating-TriPoint values are    */
/*  assumed to be the "primary" word type.  (The "primary" word type is used */
/*  to determine alignment of items.)  If `alignment' isn't zero, all items  */
/*  will be `alignment'-byte aligned in memory.  `alignment' must be either  */
/*  a multiple or a factor of the primary word size; powers of two are safe. */
/*  `alignment' is normally used to create a few unused bits at the bottom   */
/*  of each item's pointer, in which information may be stored.              */
/*                                                                           */
/*  Don't change this routine unless you understand it.                      */
/*                                                                           */
/*****************************************************************************/

void poolinit(struct memorypool *pool, int bytecount, int itemcount, 
			  enum wordtype wtype, int alignment)
{
  int wordsize;

  /* Initialize values in the pool. */
  pool->itemwordtype = wtype;
  wordsize = (pool->itemwordtype == POINTER) ? sizeof(VOID *) : sizeof(REAL);
  /* Find the proper alignment, which must be at least as large as:   */
  /*   - The parameter `alignment'.                                   */
  /*   - The primary word type, to avoid unaligned accesses.          */
  /*   - sizeof(VOID *), so the stack of dead items can be maintained */
  /*       without unaligned accesses.                                */
  if (alignment > wordsize) {
    pool->alignbytes = alignment;
  } else {
    pool->alignbytes = wordsize;
  }
  if (sizeof(VOID *) > pool->alignbytes) {
    pool->alignbytes = sizeof(VOID *);
  }
  pool->itemwords = ((bytecount + pool->alignbytes - 1) / pool->alignbytes)
                  * (pool->alignbytes / wordsize);
  pool->itembytes = pool->itemwords * wordsize;
  pool->itemsperblock = itemcount;

  /* Allocate a block of items.  Space for `itemsperblock' items and one    */
  /*   pointer (to TriPoint to the next block) are allocated, as well as space */
  /*   to ensure alignment of the items.                                    */
  pool->firstblock = (VOID **) malloc(pool->itemsperblock * pool->itembytes
                                      + sizeof(VOID *) + pool->alignbytes);
  if (pool->firstblock == (VOID **) NULL) {
    printf("Error:  Out of memory.\n");
    exit(1);
  }
  /* Set the next block pointer to NULL. */
  *(pool->firstblock) = (VOID *) NULL;
  poolrestart(pool);
}

/*****************************************************************************/
/*                                                                           */
/*  poolrestart()   Deallocate all items in a pool.                          */
/*                                                                           */
/*  The pool is returned to its starting state, except that no memory is     */
/*  freed to the operating system.  Rather, the previously allocated blocks  */
/*  are ready to be reused.                                                  */
/*                                                                           */
/*****************************************************************************/

void poolrestart(struct memorypool *pool)
{
  unsigned long alignptr;

  pool->items = 0;
  pool->maxitems = 0;

  /* Set the currently active block. */
  pool->nowblock = pool->firstblock;
  /* Find the first item in the pool.  Increment by the size of (VOID *). */
  alignptr = (unsigned long) (pool->nowblock + 1);
  /* Align the item on an `alignbytes'-byte boundary. */
  pool->nextitem = (VOID *)
    (alignptr + (unsigned long) pool->alignbytes
     - (alignptr % (unsigned long) pool->alignbytes));
  /* There are lots of unallocated items left in this block. */
  pool->unallocateditems = pool->itemsperblock;
  /* The stack of deallocated items is empty. */
  pool->deaditemstack = (VOID *) NULL;
}

/*****************************************************************************/
/*                                                                           */
/*  pooldeinit()   Free to the operating system all memory taken by a pool.  */
/*                                                                           */
/*****************************************************************************/

void pooldeinit(struct memorypool *pool)
{
  while (pool->firstblock != (VOID **) NULL) {
    pool->nowblock = (VOID **) *(pool->firstblock);
    free(pool->firstblock);
    pool->firstblock = pool->nowblock;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  poolalloc()   Allocate space for an item.                                */
/*                                                                           */
/*****************************************************************************/

VOID *poolalloc(struct memorypool *pool)
{
  VOID *newitem;
  VOID **newblock;
  unsigned long alignptr;

  /* First check the linked list of dead items.  If the list is not   */
  /*   empty, allocate an item from the list rather than a fresh one. */
  if (pool->deaditemstack != (VOID *) NULL) {
    newitem = pool->deaditemstack;               /* Take first item in list. */
    pool->deaditemstack = * (VOID **) pool->deaditemstack;
  } else {
    /* Check if there are any free items left in the current block. */
    if (pool->unallocateditems == 0) {
      /* Check if another block must be allocated. */
      if (*(pool->nowblock) == (VOID *) NULL) {
        /* Allocate a new block of items, pointed to by the previous block. */
        newblock = (VOID **) malloc(pool->itemsperblock * pool->itembytes
                                    + sizeof(VOID *) + pool->alignbytes);
        if (newblock == (VOID **) NULL) {
          printf("Error:  Out of memory.\n");
          exit(1);
        }
        *(pool->nowblock) = (VOID *) newblock;
        /* The next block pointer is NULL. */
        *newblock = (VOID *) NULL;
      }
      /* Move to the new block. */
      pool->nowblock = (VOID **) *(pool->nowblock);
      /* Find the first item in the block.    */
      /*   Increment by the size of (VOID *). */
      alignptr = (unsigned long) (pool->nowblock + 1);
      /* Align the item on an `alignbytes'-byte boundary. */
      pool->nextitem = (VOID *)
        (alignptr + (unsigned long) pool->alignbytes
         - (alignptr % (unsigned long) pool->alignbytes));
      /* There are lots of unallocated items left in this block. */
      pool->unallocateditems = pool->itemsperblock;
    }
    /* Allocate a new item. */
    newitem = pool->nextitem;
    /* Advance `nextitem' pointer to next free item in block. */
    if (pool->itemwordtype == POINTER) {
      pool->nextitem = (VOID *) ((VOID **) pool->nextitem + pool->itemwords);
    } else {
      pool->nextitem = (VOID *) ((REAL *) pool->nextitem + pool->itemwords);
    }
    pool->unallocateditems--;
    pool->maxitems++;
  }
  pool->items++;
  return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  pooldealloc()   Deallocate space for an item.                            */
/*                                                                           */
/*  The deallocated space is stored in a queue for later reuse.              */
/*                                                                           */
/*****************************************************************************/

void pooldealloc(struct memorypool *pool, VOID *dyingitem)
{
  /* Push freshly killed item onto stack. */
  *((VOID **) dyingitem) = pool->deaditemstack;
  pool->deaditemstack = dyingitem;
  pool->items--;
}

/*****************************************************************************/
/*                                                                           */
/*  traversalinit()   Prepare to traverse the entire list of items.          */
/*                                                                           */
/*  This routine is used in conjunction with traverse().                     */
/*                                                                           */
/*****************************************************************************/

void traversalinit(struct memorypool *pool)
{
  unsigned long alignptr;

  /* Begin the traversal in the first block. */
  pool->pathblock = pool->firstblock;
  /* Find the first item in the block.  Increment by the size of (VOID *). */
  alignptr = (unsigned long) (pool->pathblock + 1);
  /* Align with item on an `alignbytes'-byte boundary. */
  pool->pathitem = (VOID *)
    (alignptr + (unsigned long) pool->alignbytes
     - (alignptr % (unsigned long) pool->alignbytes));
  /* Set the number of items left in the current block. */
  pool->pathitemsleft = pool->itemsperblock;
}

/*****************************************************************************/
/*                                                                           */
/*  traverse()   Find the next item in the list.                             */
/*                                                                           */
/*  This routine is used in conjunction with traversalinit().  Be forewarned */
/*  that this routine successively returns all items in the list, including  */
/*  deallocated ones on the deaditemqueue.  It's up to you to figure out     */
/*  which ones are actually dead.  Why?  I don't want to allocate extra      */
/*  space just to demarcate dead items.  It can usually be done more         */
/*  space-efficiently by a routine that knows something about the structure  */
/*  of the item.                                                             */
/*                                                                           */
/*****************************************************************************/

VOID *traverse(struct memorypool *pool)
{
  VOID *newitem;
  unsigned long alignptr;

  /* Stop upon exhausting the list of items. */
  if (pool->pathitem == pool->nextitem) {
    return (VOID *) NULL;
  }
  /* Check whether any untraversed items remain in the current block. */
  if (pool->pathitemsleft == 0) {
    /* Find the next block. */
    pool->pathblock = (VOID **) *(pool->pathblock);
    /* Find the first item in the block.  Increment by the size of (VOID *). */
    alignptr = (unsigned long) (pool->pathblock + 1);
    /* Align with item on an `alignbytes'-byte boundary. */
    pool->pathitem = (VOID *)
      (alignptr + (unsigned long) pool->alignbytes
       - (alignptr % (unsigned long) pool->alignbytes));
    /* Set the number of items left in the current block. */
    pool->pathitemsleft = pool->itemsperblock;
  }
  newitem = pool->pathitem;
  /* Find the next item in the block. */
  if (pool->itemwordtype == POINTER) {
    pool->pathitem = (VOID *) ((VOID **) pool->pathitem + pool->itemwords);
  } else {
    pool->pathitem = (VOID *) ((REAL *) pool->pathitem + pool->itemwords);
  }
  pool->pathitemsleft--;
  return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  dummyinit()   Initialize the TriTriangle that fills "outer space" and the   */
/*                omnipresent shell edge.                                    */
/*                                                                           */
/*  The TriTriangle that fills "outer space", called `dummytri', is pointed to  */
/*  by every TriTriangle and shell edge on a boundary (be it outer or inner) of */
/*  the triangulation.  Also, `dummytri' points to one of the triangles on   */
/*  the convex hull (until the holes and concavities are carved), making it  */
/*  possible to find a starting TriTriangle for TriPoint location.                 */
/*                                                                           */
/*  The omnipresent shell edge, `dummysh', is pointed to by every TriTriangle   */
/*  or shell edge that doesn't have a full complement of real shell edges    */
/*  to TriPoint to.                                                             */
/*                                                                           */
/*****************************************************************************/

void dummyinit(int trianglewords, int shellewords)
{
  unsigned long alignptr;

  /* `triwords' and `shwords' are used by the mesh manipulation primitives */
  /*   to extract orientations of triangles and shell edges from pointers. */
  triwords = trianglewords;       /* Initialize `triwords' once and for all. */
  shwords = shellewords;           /* Initialize `shwords' once and for all. */

  /* Set up `dummytri', the `TriTriangle' that occupies "outer space". */
  dummytribase = (TriTriangle *) malloc(triwords * sizeof(TriTriangle)
                                     + triangles.alignbytes);
  if (dummytribase == (TriTriangle *) NULL) {
    printf("Error:  Out of memory.\n");
    exit(1);
  }
  /* Align `dummytri' on a `triangles.alignbytes'-byte boundary. */
  alignptr = (unsigned long) dummytribase;
  dummytri = (TriTriangle *)
    (alignptr + (unsigned long) triangles.alignbytes
     - (alignptr % (unsigned long) triangles.alignbytes));
  /* Initialize the three adjoining triangles to be "outer space".  These  */
  /*   will eventually be changed by various bonding operations, but their */
  /*   values don't really matter, as long as they can legally be          */
  /*   dereferenced.                                                       */
  dummytri[0] = (TriTriangle) dummytri;
  dummytri[1] = (TriTriangle) dummytri;
  dummytri[2] = (TriTriangle) dummytri;
  /* Three NULL vertex points. */
  dummytri[3] = (TriTriangle) NULL;
  dummytri[4] = (TriTriangle) NULL;
  dummytri[5] = (TriTriangle) NULL;

  if (useshelles) {
    /* Set up `dummysh', the omnipresent "shell edge" pointed to by any      */
    /*   TriTriangle side or shell edge end that isn't attached to a real shell */
    /*   edge.                                                               */
    dummyshbase = (TriShell *) malloc(shwords * sizeof(TriShell)
                                    + shelles.alignbytes);
    if (dummyshbase == (TriShell *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
    /* Align `dummysh' on a `shelles.alignbytes'-byte boundary. */
    alignptr = (unsigned long) dummyshbase;
    dummysh = (TriShell *)
      (alignptr + (unsigned long) shelles.alignbytes
       - (alignptr % (unsigned long) shelles.alignbytes));
    /* Initialize the two adjoining shell edges to be the omnipresent shell */
    /*   edge.  These will eventually be changed by various bonding         */
    /*   operations, but their values don't really matter, as long as they  */
    /*   can legally be dereferenced.                                       */
    dummysh[0] = (TriShell) dummysh;
    dummysh[1] = (TriShell) dummysh;
    /* Two NULL vertex points. */
    dummysh[2] = (TriShell) NULL;
    dummysh[3] = (TriShell) NULL;
    /* Initialize the two adjoining triangles to be "outer space". */
    dummysh[4] = (TriShell) dummytri;
    dummysh[5] = (TriShell) dummytri;
    /* Set the boundary marker to zero. */
    * (int *) (dummysh + 6) = 0;

    /* Initialize the three adjoining shell edges of `dummytri' to be */
    /*   the omnipresent shell edge.                                  */
    dummytri[6] = (TriTriangle) dummysh;
    dummytri[7] = (TriTriangle) dummysh;
    dummytri[8] = (TriTriangle) dummysh;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  initializepointpool()   Calculate the size of the TriPoint data structure   */
/*                          and initialize its memory pool.                  */
/*                                                                           */
/*  This routine also computes the `pointmarkindex' and `point2triindex'     */
/*  indices used to find values within each TriPoint.                           */
/*                                                                           */
/*****************************************************************************/

void initializepointpool()
{
  int pointsize;

  /* The index within each TriPoint at which the boundary marker is found.  */
  /*   Ensure the TriPoint marker is aligned to a sizeof(int)-byte address. */
  pointmarkindex = ((mesh_dim + nextras) * sizeof(REAL) + sizeof(int) - 1)
                 / sizeof(int);
  pointsize = (pointmarkindex + 1) * sizeof(int);
  if (poly) {
    /* The index within each TriPoint at which a TriTriangle pointer is found.   */
    /*   Ensure the pointer is aligned to a sizeof(TriTriangle)-byte address. */
    point2triindex = (pointsize + sizeof(TriTriangle) - 1) / sizeof(TriTriangle);
    pointsize = (point2triindex + 1) * sizeof(TriTriangle);
  }
  /* Initialize the pool of points. */
  poolinit(&points, pointsize, POINTPERBLOCK,
           (sizeof(REAL) >= sizeof(TriTriangle)) ? FLOATINGPOINT : POINTER, 0);
}

/*****************************************************************************/
/*                                                                           */
/*  initializetrisegpools()   Calculate the sizes of the TriTriangle and shell  */
/*                            edge data structures and initialize their      */
/*                            memory pools.                                  */
/*                                                                           */
/*  This routine also computes the `highorderindex', `elemattribindex', and  */
/*  `areaboundindex' indices used to find values within each TriTriangle.       */
/*                                                                           */
/*****************************************************************************/

void initializetrisegpools()
{
  int trisize;

  /* The index within each TriTriangle at which the extra nodes (above three)  */
  /*   associated with high order elements are found.  There are three      */
  /*   pointers to other triangles, three pointers to corners, and possibly */
  /*   three pointers to shell edges before the extra nodes.                */
  highorderindex = 6 + (useshelles * 3);
  /* The number of bytes occupied by a TriTriangle. */
  trisize = ((order + 1) * (order + 2) / 2 + (highorderindex - 3)) *
            sizeof(TriTriangle);
  /* The index within each TriTriangle at which its attributes are found, */
  /*   where the index is measured in REALs.                           */
  elemattribindex = (trisize + sizeof(REAL) - 1) / sizeof(REAL);
  /* The index within each TriTriangle at which the maximum area constraint  */
  /*   is found, where the index is measured in REALs.  Note that if the  */
  /*   `regionattrib' flag is set, an additional attribute will be added. */
  areaboundindex = elemattribindex + eextras + regionattrib;
  /* If TriTriangle attributes or an area bound are needed, increase the number */
  /*   of bytes occupied by a TriTriangle.                                      */
  if (vararea) {
    trisize = (areaboundindex + 1) * sizeof(REAL);
  } else if (eextras + regionattrib > 0) {
    trisize = areaboundindex * sizeof(REAL);
  }
  /* If a Voronoi diagram or TriTriangle neighbor graph is requested, make    */
  /*   sure there's room to store an integer index in each TriTriangle.  This */
  /*   integer index can occupy the same space as the shell edges or       */
  /*   attributes or area constraint or extra nodes.                       */
  if ((voronoi || neighbors) &&
      (trisize < 6 * sizeof(TriTriangle) + sizeof(int))) {
    trisize = 6 * sizeof(TriTriangle) + sizeof(int);
  }
  /* Having determined the memory size of a TriTriangle, initialize the pool. */
  poolinit(&triangles, trisize, TRIPERBLOCK, POINTER, 4);

  if (useshelles) {
    /* Initialize the pool of shell edges. */
    poolinit(&shelles, 6 * sizeof(TriTriangle) + sizeof(int), SHELLEPERBLOCK,
             POINTER, 4);

    /* Initialize the "outer space" TriTriangle and omnipresent shell edge. */
    dummyinit(triangles.itemwords, shelles.itemwords);
  } else {
    /* Initialize the "outer space" TriTriangle. */
    dummyinit(triangles.itemwords, 0);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  triangledealloc()   Deallocate space for a TriTriangle, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void triangledealloc(TriTriangle *dyingtriangle)
{
  /* Set TriTriangle's vertices to NULL.  This makes it possible to        */
  /*   detect dead triangles when traversing the list of all triangles. */
  dyingtriangle[3] = (TriTriangle) NULL;
  dyingtriangle[4] = (TriTriangle) NULL;
  dyingtriangle[5] = (TriTriangle) NULL;
  pooldealloc(&triangles, (VOID *) dyingtriangle);
}

/*****************************************************************************/
/*                                                                           */
/*  triangletraverse()   Traverse the triangles, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

TriTriangle *triangletraverse()
{
  TriTriangle *newtriangle;

  do {
    newtriangle = (TriTriangle *) traverse(&triangles);
    if (newtriangle == (TriTriangle *) NULL) {
      return (TriTriangle *) NULL;
    }
  } while (newtriangle[3] == (TriTriangle) NULL);            /* Skip dead ones. */
  return newtriangle;
}

/*****************************************************************************/
/*                                                                           */
/*  shelledealloc()   Deallocate space for a shell edge, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void shelledealloc(TriShell *dyingshelle)
{
  /* Set shell edge's vertices to NULL.  This makes it possible to */
  /*   detect dead shells when traversing the list of all shells.  */
  dyingshelle[2] = (TriShell) NULL;
  dyingshelle[3] = (TriShell) NULL;
  pooldealloc(&shelles, (VOID *) dyingshelle);
}

/*****************************************************************************/
/*                                                                           */
/*  shelletraverse()   Traverse the shell edges, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

TriShell *shelletraverse()
{
  TriShell *newshelle;

  do {
    newshelle = (TriShell *) traverse(&shelles);
    if (newshelle == (TriShell *) NULL) {
      return (TriShell *) NULL;
    }
  } while (newshelle[2] == (TriShell) NULL);                /* Skip dead ones. */
  return newshelle;
}

/*****************************************************************************/
/*                                                                           */
/*  pointdealloc()   Deallocate space for a TriPoint, marking it dead.          */
/*                                                                           */
/*****************************************************************************/

void pointdealloc(TriPoint dyingpoint)
{
  /* Mark the TriPoint as dead.  This makes it possible to detect dead points */
  /*   when traversing the list of all points.                             */
  setpointmark(dyingpoint, DEADPOINT);
  pooldealloc(&points, (VOID *) dyingpoint);
}

/*****************************************************************************/
/*                                                                           */
/*  pointtraverse()   Traverse the points, skipping dead ones.               */
/*                                                                           */
/*****************************************************************************/

TriPoint pointtraverse()
{
  TriPoint newpoint;

  do {
    newpoint = (TriPoint) traverse(&points);
    if (newpoint == (TriPoint) NULL) {
      return (TriPoint) NULL;
    }
  } while (pointmark(newpoint) == DEADPOINT);             /* Skip dead ones. */
  return newpoint;
}


/*****************************************************************************/
/*                                                                           */
/*  getpoint()   Get a specific TriPoint, by number, from the list.             */
/*                                                                           */
/*  The first TriPoint is number 'firstnumber'.                                 */
/*                                                                           */
/*  Note that this takes O(n) time (with a small constant, if POINTPERBLOCK  */
/*  is large).  I don't care to take the trouble to make it work in constant */
/*  time.                                                                    */
/*                                                                           */
/*****************************************************************************/

TriPoint getpoint(int number)
{
  VOID **getblock;
  TriPoint foundpoint;
  unsigned long alignptr;
  int current;

  getblock = points.firstblock;
  current = firstnumber;
  /* Find the right block. */
  while (current + points.itemsperblock <= number) {
    getblock = (VOID **) *getblock;
    current += points.itemsperblock;
  }
  /* Now find the right TriPoint. */
  alignptr = (unsigned long) (getblock + 1);
  foundpoint = (TriPoint) (alignptr + (unsigned long) points.alignbytes
                        - (alignptr % (unsigned long) points.alignbytes));
  while (current < number) {
    foundpoint += points.itemwords;
    current++;
  }
  return foundpoint;
}

/*****************************************************************************/
/*                                                                           */
/*  triangledeinit()   Free all remaining allocated memory.                  */
/*                                                                           */
/*****************************************************************************/

void triangledeinit()
{
  pooldeinit(&triangles);
  free(dummytribase);
  if (useshelles) {
    pooldeinit(&shelles);
    free(dummyshbase);
  }
  pooldeinit(&points);
}

/**                                                                         **/
/**                                                                         **/
/********* Memory management routines end here                       *********/

/********* Constructors begin here                                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  maketriangle()   Create a new TriTriangle with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

void maketriangle(TriOrientedTriangle *newTriOrientedTriangle)
{
  int i;

  newTriOrientedTriangle->tri = (TriTriangle *) poolalloc(&triangles);
  /* Initialize the three adjoining triangles to be "outer space". */
  newTriOrientedTriangle->tri[0] = (TriTriangle) dummytri;
  newTriOrientedTriangle->tri[1] = (TriTriangle) dummytri;
  newTriOrientedTriangle->tri[2] = (TriTriangle) dummytri;
  /* Three NULL vertex points. */
  newTriOrientedTriangle->tri[3] = (TriTriangle) NULL;
  newTriOrientedTriangle->tri[4] = (TriTriangle) NULL;
  newTriOrientedTriangle->tri[5] = (TriTriangle) NULL;
  /* Initialize the three adjoining shell edges to be the omnipresent */
  /*   shell edge.                                                    */
  if (useshelles) {
    newTriOrientedTriangle->tri[6] = (TriTriangle) dummysh;
    newTriOrientedTriangle->tri[7] = (TriTriangle) dummysh;
    newTriOrientedTriangle->tri[8] = (TriTriangle) dummysh;
  }
  for (i = 0; i < eextras; i++) {
    setelemattribute(*newTriOrientedTriangle, i, 0.0);
  }
  if (vararea) {
    setareabound(*newTriOrientedTriangle, -1.0);
  }

  newTriOrientedTriangle->orient = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  makeshelle()   Create a new shell edge with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

void makeshelle(TriOrientedShell *newedge)
{
  newedge->sh = (TriShell *) poolalloc(&shelles);
  /* Initialize the two adjoining shell edges to be the omnipresent */
  /*   shell edge.                                                  */
  newedge->sh[0] = (TriShell) dummysh;
  newedge->sh[1] = (TriShell) dummysh;
  /* Two NULL vertex points. */
  newedge->sh[2] = (TriShell) NULL;
  newedge->sh[3] = (TriShell) NULL;
  /* Initialize the two adjoining triangles to be "outer space". */
  newedge->sh[4] = (TriShell) dummytri;
  newedge->sh[5] = (TriShell) dummytri;
  /* Set the boundary marker to zero. */
  setmark(*newedge, 0);

  newedge->shorient = 0;
}

/**                                                                         **/
/**                                                                         **/
/********* Constructors end here                                     *********/

/********* Determinant evaluation routines begin here                *********/
/**                                                                         **/
/**                                                                         **/

/* The adaptive exact arithmetic geometric predicates implemented herein are */
/*   described in detail in my Technical Report CMU-CS-96-140.  The complete */
/*   reference is given in the header.                                       */

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C.              */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-TriPoint arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-TriPoint error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-TriPoint numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-TriPoint arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void exactinit()
{
  REAL half;
  REAL check, lastcheck;
  int every_other;

  every_other = 1;
  half = 0.5;
  epsilon = 1.0;
  splitter = 1.0;
  check = 1.0;
  /* Repeatedly divide `epsilon' by two until it is too small to add to      */
  /*   one without causing roundoff.  (Also check if the sum is equal to     */
  /*   the previous sum, for machines that round up instead of using exact   */
  /*   rounding.  Not that these routines will work on such machines anyway. */
  do {
    lastcheck = check;
    epsilon *= half;
    if (every_other) {
      splitter *= 2.0;
    }
    every_other = !every_other;
    check = 1.0 + epsilon;
  } while ((check != 1.0) && (check != lastcheck));
  splitter += 1.0;
  if (verbose > 1) {
    printf("Floating TriPoint roundoff is of magnitude %.17g\n", epsilon);
    printf("Floating TriPoint splitter is %.17g\n", splitter);
  }
  /* Error bounds for orientation and incircle tests. */
  resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
  ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
  ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
  ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
  iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
  iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
  iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See my Robust Predicates paper for details.             */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h)  /* h cannot be e or f. */
{
  REAL Q;
  INEXACT REAL Qnew;
  INEXACT REAL hh;
  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  int eindex, findex, hindex;
  REAL enow, fnow;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      Fast_Two_Sum(enow, Q, Qnew, hh);
      enow = e[++eindex];
    } else {
      Fast_Two_Sum(fnow, Q, Qnew, hh);
      fnow = f[++findex];
    }
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        Two_Sum(Q, enow, Qnew, hh);
        enow = e[++eindex];
      } else {
        Two_Sum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
      }
      Q = Qnew;
      if (hh != 0.0) {
        h[hindex++] = hh;
      }
    }
  }
  while (eindex < elen) {
    Two_Sum(Q, enow, Qnew, hh);
    enow = e[++eindex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  while (findex < flen) {
    Two_Sum(Q, fnow, Qnew, hh);
    fnow = f[++findex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See my Robust Predicates paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)   
/* e and h cannot be the same. */
{
  INEXACT REAL Q, sum;
  REAL hh;
  INEXACT REAL product1;
  REAL product0;
  int eindex, hindex;
  REAL enow;
  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;

  Split(b, bhi, blo);
  Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
  hindex = 0;
  if (hh != 0) {
    h[hindex++] = hh;
  }
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
    Two_Sum(Q, product0, sum, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
    Fast_Two_Sum(product1, sum, Q, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

REAL estimate(int elen, REAL *e)
{
  REAL Q;
  int eindex;

  Q = e[0];
  for (eindex = 1; eindex < elen; eindex++) {
    Q += e[eindex];
  }
  return Q;
}

/*****************************************************************************/
/*                                                                           */
/*  counterclockwise()   Return a positive value if the points pa, pb, and   */
/*                       pc occur in counterclockwise order; a negative      */
/*                       value if they occur in clockwise order; and zero    */
/*                       if they are collinear.  The result is also a rough  */
/*                       approximation of twice the signed area of the       */
/*                       TriTriangle defined by the three points.               */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are collinear or nearly so.            */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

REAL counterclockwiseadapt(TriPoint pa, TriPoint pb, TriPoint pc, REAL detsum)
{
  INEXACT REAL acx, acy, bcx, bcy;
  REAL acxtail, acytail, bcxtail, bcytail;
  INEXACT REAL detleft, detright;
  REAL detlefttail, detrighttail;
  REAL det, errbound;
  REAL B[4], C1[8], C2[12], D[16];
  INEXACT REAL B3;
  int C1length, C2length, Dlength;
  REAL u[4];
  INEXACT REAL u3;
  INEXACT REAL s1, t1;
  REAL s0, t0;

  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;
  INEXACT REAL _i, _j;
  REAL _0;

  acx = (REAL) (pa[0] - pc[0]);
  bcx = (REAL) (pb[0] - pc[0]);
  acy = (REAL) (pa[1] - pc[1]);
  bcy = (REAL) (pb[1] - pc[1]);

  Two_Product(acx, bcy, detleft, detlefttail);
  Two_Product(acy, bcx, detright, detrighttail);

  Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
               B3, B[2], B[1], B[0]);
  B[3] = B3;

  det = estimate(4, B);
  errbound = ccwerrboundB * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
  Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
  Two_Diff_Tail(pa[1], pc[1], acy, acytail);
  Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det;
  }

  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
  det += (acx * bcytail + bcy * acxtail)
       - (acy * bcxtail + bcx * acytail);
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Product(acxtail, bcy, s1, s0);
  Two_Product(acytail, bcx, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

  Two_Product(acx, bcytail, s1, s0);
  Two_Product(acy, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

  Two_Product(acxtail, bcytail, s1, s0);
  Two_Product(acytail, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

  return(D[Dlength - 1]);
}

REAL counterclockwise(TriPoint pa, TriPoint pb, TriPoint pc)
{
  REAL detleft, detright, det;
  REAL detsum, errbound;

  counterclockcount++;

  detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
  detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
  det = detleft - detright;

  if (noexact) {
    return det;
  }

  if (detleft > 0.0) {
    if (detright <= 0.0) {
      return det;
    } else {
      detsum = detleft + detright;
    }
  } else if (detleft < 0.0) {
    if (detright >= 0.0) {
      return det;
    } else {
      detsum = -detleft - detright;
    }
  } else {
    return det;
  }

  errbound = ccwerrboundA * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  return counterclockwiseadapt(pa, pb, pc, detsum);
}

/*****************************************************************************/
/*                                                                           */
/*  incircle()   Return a positive value if the TriPoint pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are cocircular or nearly so.           */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

REAL incircleadapt(TriPoint pa, TriPoint pb, TriPoint pc, TriPoint pd, REAL permanent)
{
  INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
  REAL det, errbound;

  INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  REAL bc[4], ca[4], ab[4];
  INEXACT REAL bc3, ca3, ab3;
  REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
  int axbclen, axxbclen, aybclen, ayybclen, alen;
  REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
  int bxcalen, bxxcalen, bycalen, byycalen, blen;
  REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
  int cxablen, cxxablen, cyablen, cyyablen, clen;
  REAL abdet[64];
  int ablen;
  REAL fin1[1152], fin2[1152];
  REAL *finnow, *finother, *finswap;
  int finlength;

  REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
  INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
  REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
  REAL aa[4], bb[4], cc[4];
  INEXACT REAL aa3, bb3, cc3;
  INEXACT REAL ti1, tj1;
  REAL ti0, tj0;
  REAL u[4], v[4];
  INEXACT REAL u3, v3;
  REAL temp8[8], temp16a[16], temp16b[16], temp16c[16];
  REAL temp32a[32], temp32b[32], temp48[48], temp64[64];
  int temp8len, temp16alen, temp16blen, temp16clen;
  int temp32alen, temp32blen, temp48len, temp64len;
  REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8];
  int axtbblen, axtcclen, aytbblen, aytcclen;
  REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
  int bxtaalen, bxtcclen, bytaalen, bytcclen;
  REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
  int cxtaalen, cxtbblen, cytaalen, cytbblen;
  REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
  int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
  REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
  int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
  REAL axtbctt[8], aytbctt[8], bxtcatt[8];
  REAL bytcatt[8], cxtabtt[8], cytabtt[8];
  int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
  REAL abt[8], bct[8], cat[8];
  int abtlen, bctlen, catlen;
  REAL abtt[4], bctt[4], catt[4];
  int abttlen, bcttlen, cattlen;
  INEXACT REAL abtt3, bctt3, catt3;
  REAL negate;

  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;
  INEXACT REAL _i, _j;
  REAL _0;

  adx = (REAL) (pa[0] - pd[0]);
  bdx = (REAL) (pb[0] - pd[0]);
  cdx = (REAL) (pc[0] - pd[0]);
  ady = (REAL) (pa[1] - pd[1]);
  bdy = (REAL) (pb[1] - pd[1]);
  cdy = (REAL) (pc[1] - pd[1]);

  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
  axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
  aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
  ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
  alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

  Two_Product(cdx, ady, cdxady1, cdxady0);
  Two_Product(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
  bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
  bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
  byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
  blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

  Two_Product(adx, bdy, adxbdy1, adxbdy0);
  Two_Product(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
  cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
  cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
  cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
  clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = iccerrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
  Two_Diff_Tail(pa[1], pd[1], ady, adytail);
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
    return det;
  }

  errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
  det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
                                     - (bdy * cdxtail + cdx * bdytail))
          + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
       + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
                                     - (cdy * adxtail + adx * cdytail))
          + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
       + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
                                     - (ady * bdxtail + bdx * adytail))
          + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if ((bdxtail != 0.0) || (bdytail != 0.0)
      || (cdxtail != 0.0) || (cdytail != 0.0)) {
    Square(adx, adxadx1, adxadx0);
    Square(ady, adyady1, adyady0);
    Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
    aa[3] = aa3;
  }
  if ((cdxtail != 0.0) || (cdytail != 0.0)
      || (adxtail != 0.0) || (adytail != 0.0)) {
    Square(bdx, bdxbdx1, bdxbdx0);
    Square(bdy, bdybdy1, bdybdy0);
    Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
    bb[3] = bb3;
  }
  if ((adxtail != 0.0) || (adytail != 0.0)
      || (bdxtail != 0.0) || (bdytail != 0.0)) {
    Square(cdx, cdxcdx1, cdxcdx0);
    Square(cdy, cdycdy1, cdycdy0);
    Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
    cc[3] = cc3;
  }

  if (adxtail != 0.0) {
    axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
    temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx,
                                          temp16a);

    axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
    temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

    axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
    temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (adytail != 0.0) {
    aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
    temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady,
                                          temp16a);

    aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
    temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

    aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
    temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdxtail != 0.0) {
    bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
    temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,
                                          temp16a);

    bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
    temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

    bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
    temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdytail != 0.0) {
    bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
    temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy,
                                          temp16a);

    bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
    temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

    bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
    temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdxtail != 0.0) {
    cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
    temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx,
                                          temp16a);

    cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
    temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

    cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
    temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdytail != 0.0) {
    cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
    temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy,
                                          temp16a);

    cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
    temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

    cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
    temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if ((adxtail != 0.0) || (adytail != 0.0)) {
    if ((bdxtail != 0.0) || (bdytail != 0.0)
        || (cdxtail != 0.0) || (cdytail != 0.0)) {
      Two_Product(bdxtail, cdy, ti1, ti0);
      Two_Product(bdx, cdytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -bdy;
      Two_Product(cdxtail, negate, ti1, ti0);
      negate = -bdytail;
      Two_Product(cdx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

      Two_Product(bdxtail, cdytail, ti1, ti0);
      Two_Product(cdxtail, bdytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
      bctt[3] = bctt3;
      bcttlen = 4;
    } else {
      bct[0] = 0.0;
      bctlen = 1;
      bctt[0] = 0.0;
      bcttlen = 1;
    }

    if (adxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
      axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (cdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail,
                                            temp32a);
      axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
      temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (adytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
      aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail,
                                            temp32a);
      aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
      temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }
  if ((bdxtail != 0.0) || (bdytail != 0.0)) {
    if ((cdxtail != 0.0) || (cdytail != 0.0)
        || (adxtail != 0.0) || (adytail != 0.0)) {
      Two_Product(cdxtail, ady, ti1, ti0);
      Two_Product(cdx, adytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -cdy;
      Two_Product(adxtail, negate, ti1, ti0);
      negate = -cdytail;
      Two_Product(adx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

      Two_Product(cdxtail, adytail, ti1, ti0);
      Two_Product(adxtail, cdytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
      catt[3] = catt3;
      cattlen = 4;
    } else {
      cat[0] = 0.0;
      catlen = 1;
      catt[0] = 0.0;
      cattlen = 1;
    }

    if (bdxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
      bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (adytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail,
                                            temp32a);
      bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
      temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
      bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail,
                                            temp32a);
      bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
      temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }
  if ((cdxtail != 0.0) || (cdytail != 0.0)) {
    if ((adxtail != 0.0) || (adytail != 0.0)
        || (bdxtail != 0.0) || (bdytail != 0.0)) {
      Two_Product(adxtail, bdy, ti1, ti0);
      Two_Product(adx, bdytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -ady;
      Two_Product(bdxtail, negate, ti1, ti0);
      negate = -adytail;
      Two_Product(bdx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

      Two_Product(adxtail, bdytail, ti1, ti0);
      Two_Product(bdxtail, adytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
      abtt[3] = abtt3;
      abttlen = 4;
    } else {
      abt[0] = 0.0;
      abtlen = 1;
      abtt[0] = 0.0;
      abttlen = 1;
    }

    if (cdxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
      cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (bdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail,
                                            temp32a);
      cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
      temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
      cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail,
                                            temp32a);
      cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
      temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }

  return finnow[finlength - 1];
}

REAL incircle(TriPoint pa, TriPoint pb, TriPoint pc, TriPoint pd)
{
  REAL adx, bdx, cdx, ady, bdy, cdy;
  REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  REAL alift, blift, clift;
  REAL det;
  REAL permanent, errbound;

  incirclecount++;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;
  alift = adx * adx + ady * ady;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;
  blift = bdx * bdx + bdy * bdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;
  clift = cdx * cdx + cdy * cdy;

  det = alift * (bdxcdy - cdxbdy)
      + blift * (cdxady - adxcdy)
      + clift * (adxbdy - bdxady);

  if (noexact) {
    return det;
  }

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift
            + (Absolute(cdxady) + Absolute(adxcdy)) * blift
            + (Absolute(adxbdy) + Absolute(bdxady)) * clift;
  errbound = iccerrboundA * permanent;
  if ((det > errbound) || (-det > errbound)) {
    return det;
  }

  return incircleadapt(pa, pb, pc, pd, permanent);
}

/**                                                                         **/
/**                                                                         **/
/********* Determinant evaluation routines end here                  *********/

/*****************************************************************************/
/*                                                                           */
/*  triangleinit()   Initialize some variables.                              */
/*                                                                           */
/*****************************************************************************/

void triangleinit()
{
  points.maxitems = triangles.maxitems = shelles.maxitems = viri.maxitems =
    badsegments.maxitems = badtriangles.maxitems = splaynodes.maxitems = 0l;
  points.itembytes = triangles.itembytes = shelles.itembytes = viri.itembytes =
    badsegments.itembytes = badtriangles.itembytes = splaynodes.itembytes = 0;
  recenttri.tri = (TriTriangle *) NULL;    /* No TriTriangle has been visited yet. */
  samples = 1;            /* Point location should take at least one sample. */
  checksegments = 0;      /* There are no segments in the triangulation yet. */
  incirclecount = counterclockcount = hyperbolacount = 0;
  circumcentercount = circletopcount = 0;
  randomseed = 1;

  exactinit();                     /* Initialize exact arithmetic constants. */
}

/*****************************************************************************/
/*                                                                           */
/*  randomnation()   Generate a random number between 0 and `choices' - 1.   */
/*                                                                           */
/*  This is a simple linear congruential random number generator.  Hence, it */
/*  is a bad random number generator, but good enough for most randomized    */
/*  geometric algorithms.                                                    */
/*                                                                           */
/*****************************************************************************/

unsigned long randomnation(unsigned int choices)
{
  randomseed = (randomseed * 1366l + 150889l) % 714025l;
  return randomseed / (714025l / choices + 1);
}


/********* Point location routines begin here                        *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  makepointmap()   Construct a mapping from points to triangles to improve  */
/*                  the speed of TriPoint location for segment insertion.       */
/*                                                                           */
/*  Traverses all the triangles, and provides each corner of each TriTriangle   */
/*  with a pointer to that TriTriangle.  Of course, pointers will be            */
/*  overwritten by other pointers because (almost) each TriPoint is a corner    */
/*  of several triangles, but in the end every TriPoint will TriPoint to some      */
/*  TriTriangle that contains it.                                               */
/*                                                                           */
/*****************************************************************************/

void makepointmap()
{
  TriOrientedTriangle triangleloop;
  TriPoint triorg;

  if (verbose) {
    printf("    Constructing mapping from points to triangles.\n");
  }
  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  while (triangleloop.tri != (TriTriangle *) NULL) {
    /* Check all three points of the TriTriangle. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      org(triangleloop, triorg);
      setpoint2tri(triorg, encode(triangleloop));
    }
    triangleloop.tri = triangletraverse();
  }
}

/*****************************************************************************/
/*                                                                           */
/*  preciselocate()   Find a TriTriangle or edge containing a given TriPoint.      */
/*                                                                           */
/*  Begins its search from `searchtri'.  It is important that `searchtri'    */
/*  be a handle with the property that `searchpoint' is strictly to the left */
/*  of the edge denoted by `searchtri', or is collinear with that edge and   */
/*  does not intersect that edge.  (In particular, `searchpoint' should not  */
/*  be the origin or destination of that edge.)                              */
/*                                                                           */
/*  These conditions are imposed because preciselocate() is normally used in */
/*  one of two situations:                                                   */
/*                                                                           */
/*  (1)  To try to find the location to insert a new TriPoint.  Normally, we    */
/*       know an edge that the TriPoint is strictly to the left of.  In the     */
/*       incremental Delaunay algorithm, that edge is a bounding box edge.   */
/*       In Ruppert's Delaunay refinement algorithm for quality meshing,     */
/*       that edge is the shortest edge of the TriTriangle whose circumcenter   */
/*       is being inserted.                                                  */
/*                                                                           */
/*  (2)  To try to find an existing TriPoint.  In this case, any edge on the    */
/*       convex hull is a good starting edge.  The possibility that the      */
/*       vertex one seeks is an endpoint of the starting edge must be        */
/*       screened out before preciselocate() is called.                      */
/*                                                                           */
/*  On completion, `searchtri' is a TriTriangle that contains `searchpoint'.    */
/*                                                                           */
/*  This implementation differs from that given by Guibas and Stolfi.  It    */
/*  walks from TriTriangle to TriTriangle, crossing an edge only if `searchpoint'  */
/*  is on the other side of the line containing that edge.  After entering   */
/*  a TriTriangle, there are two edges by which one can leave that TriTriangle.    */
/*  If both edges are valid (`searchpoint' is on the other side of both      */
/*  edges), one of the two is chosen by drawing a line perpendicular to      */
/*  the entry edge (whose endpoints are `forg' and `fdest') passing through  */
/*  `fapex'.  Depending on which side of this perpendicular `searchpoint'    */
/*  falls on, an exit edge is chosen.                                        */
/*                                                                           */
/*  This implementation is empirically faster than the Guibas and Stolfi     */
/*  TriPoint location routine (which I originally used), which tends to spiral  */
/*  in toward its target.                                                    */
/*                                                                           */
/*  Returns ONVERTEX if the TriPoint lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the TriPoint lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the TriPoint lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the TriPoint lies strictly within a TriTriangle.         */
/*  `searchtri' is a handle on the TriTriangle that contains the TriPoint.         */
/*                                                                           */
/*  Returns OUTSIDE if the TriPoint lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the TriPoint is to the right of.  This might      */
/*  occur when the circumcenter of a TriTriangle falls just slightly outside    */
/*  the mesh due to floating-TriPoint roundoff error.  It also occurs when      */
/*  seeking a hole or region TriPoint that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*  However, it can still be used to find the circumcenter of a TriTriangle, as */
/*  long as the search is begun from the TriTriangle in question.               */
/*                                                                           */
/*****************************************************************************/

TriLocateResultType preciselocate(TriPoint searchpoint, TriOrientedTriangle *searchtri)
{
  TriOrientedTriangle backtracktri;
  TriPoint forg, fdest, fapex;
  TriPoint swappoint;
  REAL orgorient, destorient;
  int moveleft;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (verbose > 2) {
    printf("  Searching for TriPoint (%.12g, %.12g).\n",
           searchpoint[0], searchpoint[1]);
  }
  /* Where are we? */
  org(*searchtri, forg);
  dest(*searchtri, fdest);
  apex(*searchtri, fapex);
  while (1) {
    if (verbose > 2) {
      printf("    At (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             forg[0], forg[1], fdest[0], fdest[1], fapex[0], fapex[1]);
    }
    /* Check whether the apex is the TriPoint we seek. */
    if ((fapex[0] == searchpoint[0]) && (fapex[1] == searchpoint[1])) {
      lprevself(*searchtri);
      return ONVERTEX;
    }
    /* Does the TriPoint lie on the other side of the line defined by the */
    /*   TriTriangle edge opposite the TriTriangle's destination?            */
    destorient = counterclockwise(forg, fapex, searchpoint);
    /* Does the TriPoint lie on the other side of the line defined by the */
    /*   TriTriangle edge opposite the TriTriangle's origin?                 */
    orgorient = counterclockwise(fapex, fdest, searchpoint);
    if (destorient > 0.0) {
      if (orgorient > 0.0) {
        /* Move left if the inner product of (fapex - searchpoint) and  */
        /*   (fdest - forg) is positive.  This is equivalent to drawing */
        /*   a line perpendicular to the line (forg, fdest) passing     */
        /*   through `fapex', and determining which side of this line   */
        /*   `searchpoint' falls on.                                    */
        moveleft = (fapex[0] - searchpoint[0]) * (fdest[0] - forg[0]) +
                   (fapex[1] - searchpoint[1]) * (fdest[1] - forg[1]) > 0.0;
      } else {
        moveleft = 1;
      }
    } else {
      if (orgorient > 0.0) {
        moveleft = 0;
      } else {
        /* The TriPoint we seek must be on the boundary of or inside this */
        /*   TriTriangle.                                                 */
        if (destorient == 0.0) {
          lprevself(*searchtri);
          return ONEDGE;
        }
        if (orgorient == 0.0) {
          lnextself(*searchtri);
          return ONEDGE;
        }
        return INTRIANGLE;
      }
    }

    /* Move to another TriTriangle.  Leave a trace `backtracktri' in case */
    /*   floating-TriPoint roundoff or some such bogey causes us to walk  */
    /*   off a boundary of the triangulation.  We can just bounce off  */
    /*   the boundary as if it were an elastic band.                   */
    if (moveleft) {
      lprev(*searchtri, backtracktri);
      fdest = fapex;
    } else {
      lnext(*searchtri, backtracktri);
      forg = fapex;
    }
    sym(backtracktri, *searchtri);

    /* Check for walking off the edge. */
    if (searchtri->tri == dummytri) {
      /* Turn around. */
      TriOrientedTrianglecopy(backtracktri, *searchtri);
      swappoint = forg;
      forg = fdest;
      fdest = swappoint;
      apex(*searchtri, fapex);
      /* Check if the TriPoint really is beyond the triangulation boundary. */
      destorient = counterclockwise(forg, fapex, searchpoint);
      orgorient = counterclockwise(fapex, fdest, searchpoint);
      if ((orgorient < 0.0) && (destorient < 0.0)) {
        return OUTSIDE;
      }
    } else {
      apex(*searchtri, fapex);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  locate()   Find a TriTriangle or edge containing a given TriPoint.             */
/*                                                                           */
/*  Searching begins from one of:  the input `searchtri', a recently         */
/*  encountered TriTriangle `recenttri', or from a TriTriangle chosen from a       */
/*  random sample.  The choice is made by determining which TriTriangle's       */
/*  origin is closest to the TriPoint we are searcing for.  Normally,           */
/*  `searchtri' should be a handle on the convex hull of the triangulation.  */
/*                                                                           */
/*  Details on the random sampling method can be found in the Mucke, Saias,  */
/*  and Zhu paper cited in the header of this code.                          */
/*                                                                           */
/*  On completion, `searchtri' is a TriTriangle that contains `searchpoint'.    */
/*                                                                           */
/*  Returns ONVERTEX if the TriPoint lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the TriPoint lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the TriPoint lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the TriPoint lies strictly within a TriTriangle.         */
/*  `searchtri' is a handle on the TriTriangle that contains the TriPoint.         */
/*                                                                           */
/*  Returns OUTSIDE if the TriPoint lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the TriPoint is to the right of.  This might      */
/*  occur when the circumcenter of a TriTriangle falls just slightly outside    */
/*  the mesh due to floating-TriPoint roundoff error.  It also occurs when      */
/*  seeking a hole or region TriPoint that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*                                                                           */
/*****************************************************************************/

TriLocateResultType locate(TriPoint searchpoint, TriOrientedTriangle *searchtri)
{
  VOID **sampleblock;
  TriTriangle *firsttri;
  TriOrientedTriangle sampletri;
  TriPoint torg, tdest;
  unsigned long alignptr;
  REAL searchdist, dist;
  REAL ahead;
  long sampleblocks, samplesperblock, samplenum;
  long triblocks;
  long i, j;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (verbose > 2) {
    printf("  Randomly sampling for a TriTriangle near TriPoint (%.12g, %.12g).\n",
           searchpoint[0], searchpoint[1]);
  }
  /* Record the distance from the suggested starting TriTriangle to the */
  /*   TriPoint we seek.                                                */
  org(*searchtri, torg);
  searchdist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0])
             + (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
  if (verbose > 2) {
    printf("    Boundary TriTriangle has origin (%.12g, %.12g).\n",
           torg[0], torg[1]);
  }

  /* If a recently encountered TriTriangle has been recorded and has not been */
  /*   deallocated, test it as a good starting TriPoint.                      */
  if (recenttri.tri != (TriTriangle *) NULL) {
    if (recenttri.tri[3] != (TriTriangle) NULL) {
      org(recenttri, torg);
      if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1])) {
        TriOrientedTrianglecopy(recenttri, *searchtri);
        return ONVERTEX;
      }
      dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0])
           + (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
      if (dist < searchdist) {
        TriOrientedTrianglecopy(recenttri, *searchtri);
        searchdist = dist;
        if (verbose > 2) {
          printf("    Choosing recent TriTriangle with origin (%.12g, %.12g).\n",
                 torg[0], torg[1]);
        }
      }
    }
  }

  /* The number of random samples taken is proportional to the cube root of */
  /*   the number of triangles in the mesh.  The next bit of code assumes   */
  /*   that the number of triangles increases monotonically.                */
  while (SAMPLEFACTOR * samples * samples * samples < triangles.items) {
    samples++;
  }
  triblocks = (triangles.maxitems + TRIPERBLOCK - 1) / TRIPERBLOCK;
  samplesperblock = 1 + (samples / triblocks);
  sampleblocks = samples / samplesperblock;
  sampleblock = triangles.firstblock;
  sampletri.orient = 0;
  for (i = 0; i < sampleblocks; i++) {
    alignptr = (unsigned long) (sampleblock + 1);
    firsttri = (TriTriangle *) (alignptr + (unsigned long) triangles.alignbytes
                          - (alignptr % (unsigned long) triangles.alignbytes));
    for (j = 0; j < samplesperblock; j++) {
      if (i == triblocks - 1) {
        samplenum = randomnation((int)
                                 (triangles.maxitems - (i * TRIPERBLOCK)));
      } else {
        samplenum = randomnation(TRIPERBLOCK);
      }
      sampletri.tri = (TriTriangle *)
                      (firsttri + (samplenum * triangles.itemwords));
      if (sampletri.tri[3] != (TriTriangle) NULL) {
        org(sampletri, torg);
        dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0])
             + (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
        if (dist < searchdist) {
          TriOrientedTrianglecopy(sampletri, *searchtri);
          searchdist = dist;
          if (verbose > 2) {
            printf("    Choosing TriTriangle with origin (%.12g, %.12g).\n",
                   torg[0], torg[1]);
          }
        }
      }
    }
    sampleblock = (VOID **) *sampleblock;
  }
  /* Where are we? */
  org(*searchtri, torg);
  dest(*searchtri, tdest);
  /* Check the starting TriTriangle's vertices. */
  if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1])) {
    return ONVERTEX;
  }
  if ((tdest[0] == searchpoint[0]) && (tdest[1] == searchpoint[1])) {
    lnextself(*searchtri);
    return ONVERTEX;
  }
  /* Orient `searchtri' to fit the preconditions of calling preciselocate(). */
  ahead = counterclockwise(torg, tdest, searchpoint);
  if (ahead < 0.0) {
    /* Turn around so that `searchpoint' is to the left of the */
    /*   edge specified by `searchtri'.                        */
    symself(*searchtri);
  } else if (ahead == 0.0) {
    /* Check if `searchpoint' is between `torg' and `tdest'. */
    if (((torg[0] < searchpoint[0]) == (searchpoint[0] < tdest[0]))
        && ((torg[1] < searchpoint[1]) == (searchpoint[1] < tdest[1]))) {
      return ONEDGE;
    }
  }
  return preciselocate(searchpoint, searchtri);
}

/**                                                                         **/
/**                                                                         **/
/********* Point location routines end here                          *********/

/********* Mesh transformation routines begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  insertshelle()   Create a new shell edge and insert it between two       */
/*                   triangles.                                              */
/*                                                                           */
/*  The new shell edge is inserted at the edge described by the handle       */
/*  `tri'.  Its vertices are properly initialized.  The marker `shellemark'  */
/*  is applied to the shell edge and, if appropriate, its vertices.          */
/*                                                                           */
/*****************************************************************************/

void insertshelle(TriOrientedTriangle *tri, int shellemark)
/* tri : Edge at which to insert the new shell edge. */
/* shellemark : Marker for the new shell edge. */
{
  TriOrientedTriangle oppotri;
  TriOrientedShell newshelle;
  TriPoint triorg, tridest;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  /* Mark points if possible. */
  org(*tri, triorg);
  dest(*tri, tridest);
  if (pointmark(triorg) == 0) {
    setpointmark(triorg, shellemark);
  }
  if (pointmark(tridest) == 0) {
    setpointmark(tridest, shellemark);
  }
  /* Check if there's already a shell edge here. */
  tspivot(*tri, newshelle);
  if (newshelle.sh == dummysh) {
    /* Make new shell edge and initialize its vertices. */
    makeshelle(&newshelle);
    setsorg(newshelle, tridest);
    setsdest(newshelle, triorg);
    /* Bond new shell edge to the two triangles it is sandwiched between. */
    /*   Note that the facing TriTriangle `oppotri' might be equal to        */
    /*   `dummytri' (outer space), but the new shell edge is bonded to it */
    /*   all the same.                                                    */
    tsbond(*tri, newshelle);
    sym(*tri, oppotri);
    ssymself(newshelle);
    tsbond(oppotri, newshelle);
    setmark(newshelle, shellemark);
    if (verbose > 2) {
      printf("  Inserting new ");
      printshelle(&newshelle);
    }
  } else {
    if (mark(newshelle) == 0) {
      setmark(newshelle, shellemark);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  Terminology                                                              */
/*                                                                           */
/*  A "local transformation" replaces a small set of triangles with another  */
/*  set of triangles.  This may or may not involve inserting or deleting a   */
/*  TriPoint.                                                                   */
/*                                                                           */
/*  The term "casing" is used to describe the set of triangles that are      */
/*  attached to the triangles being transformed, but are not transformed     */
/*  themselves.  Think of the casing as a fixed hollow structure inside      */
/*  which all the action happens.  A "casing" is only defined relative to    */
/*  a single transformation; each occurrence of a transformation will        */
/*  involve a different casing.                                              */
/*                                                                           */
/*  A "shell" is similar to a "casing".  The term "shell" describes the set  */
/*  of shell edges (if any) that are attached to the triangles being         */
/*  transformed.  However, I sometimes use "shell" to refer to a single      */
/*  shell edge, so don't get confused.                                       */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  flip()   Transform two triangles to two different triangles by flipping  */
/*           an edge within a quadrilateral.                                 */
/*                                                                           */
/*  Imagine the original triangles, abc and bad, oriented so that the        */
/*  shared edge ab lies in a horizontal plane, with the TriPoint b on the left  */
/*  and the TriPoint a on the right.  The TriPoint c lies below the edge, and the  */
/*  TriPoint d lies above the edge.  The `flipedge' handle holds the edge ab    */
/*  of TriTriangle abc, and is directed left, from vertex a to vertex b.        */
/*                                                                           */
/*  The triangles abc and bad are deleted and replaced by the triangles cdb  */
/*  and dca.  The triangles that represent abc and bad are NOT deallocated;  */
/*  they are reused for dca and cdb, respectively.  Hence, any handles that  */
/*  may have held the original triangles are still valid, although not       */
/*  directed as they were before.                                            */
/*                                                                           */
/*  Upon completion of this routine, the `flipedge' handle holds the edge    */
/*  dc of TriTriangle dca, and is directed down, from vertex d to vertex c.     */
/*  (Hence, the two triangles have rotated counterclockwise.)                */
/*                                                                           */
/*  WARNING:  This transformation is geometrically valid only if the         */
/*  quadrilateral adbc is convex.  Furthermore, this transformation is       */
/*  valid only if there is not a shell edge between the triangles abc and    */
/*  bad.  This routine does not check either of these preconditions, and     */
/*  it is the responsibility of the calling routine to ensure that they are  */
/*  met.  If they are not, the streets shall be filled with wailing and      */
/*  gnashing of teeth.                                                       */
/*                                                                           */
/*****************************************************************************/

void flip(TriOrientedTriangle *flipedge)
/* flipedge : Handle for the TriTriangle abc. */
{
  TriOrientedTriangle botleft, botright;
  TriOrientedTriangle topleft, topright;
  TriOrientedTriangle top;
  TriOrientedTriangle botlcasing, botrcasing;
  TriOrientedTriangle toplcasing, toprcasing;
  TriOrientedShell botlshelle, botrshelle;
  TriOrientedShell toplshelle, toprshelle;
  TriPoint leftpoint, rightpoint, botpoint;
  TriPoint farpoint;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  /* Identify the vertices of the quadrilateral. */
  org(*flipedge, rightpoint);
  dest(*flipedge, leftpoint);
  apex(*flipedge, botpoint);
  sym(*flipedge, top);
#ifdef SELF_CHECK
  if (top.tri == dummytri) {
    printf("Internal error in flip():  Attempt to flip on boundary.\n");
    lnextself(*flipedge);
    return;
  }
  if (checksegments) {
    tspivot(*flipedge, toplshelle);
    if (toplshelle.sh != dummysh) {
      printf("Internal error in flip():  Attempt to flip a segment.\n");
      lnextself(*flipedge);
      return;
    }
  }
#endif /* SELF_CHECK */
  apex(top, farpoint);

  /* Identify the casing of the quadrilateral. */
  lprev(top, topleft);
  sym(topleft, toplcasing);
  lnext(top, topright);
  sym(topright, toprcasing);
  lnext(*flipedge, botleft);
  sym(botleft, botlcasing);
  lprev(*flipedge, botright);
  sym(botright, botrcasing);
  /* Rotate the quadrilateral one-quarter turn counterclockwise. */
  bond(topleft, botlcasing);
  bond(botleft, botrcasing);
  bond(botright, toprcasing);
  bond(topright, toplcasing);

  if (checksegments) {
    /* Check for shell edges and rebond them to the quadrilateral. */
    tspivot(topleft, toplshelle);
    tspivot(botleft, botlshelle);
    tspivot(botright, botrshelle);
    tspivot(topright, toprshelle);
    if (toplshelle.sh == dummysh) {
      tsdissolve(topright);
    } else {
      tsbond(topright, toplshelle);
    }
    if (botlshelle.sh == dummysh) {
      tsdissolve(topleft);
    } else {
      tsbond(topleft, botlshelle);
    }
    if (botrshelle.sh == dummysh) {
      tsdissolve(botleft);
    } else {
      tsbond(botleft, botrshelle);
    }
    if (toprshelle.sh == dummysh) {
      tsdissolve(botright);
    } else {
      tsbond(botright, toprshelle);
    }
  }

  /* New TriPoint assignments for the rotated quadrilateral. */
  setorg(*flipedge, farpoint);
  setdest(*flipedge, botpoint);
  setapex(*flipedge, rightpoint);
  setorg(top, botpoint);
  setdest(top, farpoint);
  setapex(top, leftpoint);
  if (verbose > 2) {
    printf("  Edge flip results in left ");
    lnextself(topleft);
    printtriangle(&topleft);
    printf("  and right ");
    printtriangle(flipedge);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  insertsite()   Insert a vertex into a Delaunay triangulation,            */
/*                 performing flips as necessary to maintain the Delaunay    */
/*                 property.                                                 */
/*                                                                           */
/*  The TriPoint `insertpoint' is located.  If `searchtri.tri' is not NULL,     */
/*  the search for the containing TriTriangle begins from `searchtri'.  If      */
/*  `searchtri.tri' is NULL, a full TriPoint location procedure is called.      */
/*  If `insertpoint' is found inside a TriTriangle, the TriTriangle is split into  */
/*  three; if `insertpoint' lies on an edge, the edge is split in two,       */
/*  thereby splitting the two adjacent triangles into four.  Edge flips are  */
/*  used to restore the Delaunay property.  If `insertpoint' lies on an      */
/*  existing vertex, no action is taken, and the value DUPLICATEPOINT is     */
/*  returned.  On return, `searchtri' is set to a handle whose origin is the */
/*  existing vertex.                                                         */
/*                                                                           */
/*  Normally, the parameter `splitedge' is set to NULL, implying that no     */
/*  segment should be split.  In this case, if `insertpoint' is found to     */
/*  lie on a segment, no action is taken, and the value VIOLATINGPOINT is    */
/*  returned.  On return, `searchtri' is set to a handle whose primary edge  */
/*  is the violated segment.                                                 */
/*                                                                           */
/*  If the calling routine wishes to split a segment by inserting a TriPoint in */
/*  it, the parameter `splitedge' should be that segment.  In this case,     */
/*  `searchtri' MUST be the TriTriangle handle reached by pivoting from that    */
/*  segment; no TriPoint location is done.                                      */
/*                                                                           */
/*  `segmentflaws' and `triflaws' are flags that indicate whether or not     */
/*  there should be checks for the creation of encroached segments or bad    */
/*  quality faces.  If a newly inserted TriPoint encroaches upon segments,      */
/*  these segments are added to the list of segments to be split if          */
/*  `segmentflaws' is set.  If bad triangles are created, these are added    */
/*  to the queue if `triflaws' is set.                                       */
/*                                                                           */
/*  If a duplicate TriPoint or violated segment does not prevent the TriPoint      */
/*  from being inserted, the return value will be ENCROACHINGPOINT if the    */
/*  TriPoint encroaches upon a segment (and checking is enabled), or            */
/*  SUCCESSFULPOINT otherwise.  In either case, `searchtri' is set to a      */
/*  handle whose origin is the newly inserted vertex.                        */
/*                                                                           */
/*  insertsite() does not use flip() for reasons of speed; some              */
/*  information can be reused from edge flip to edge flip, like the          */
/*  locations of shell edges.                                                */
/*                                                                           */
/*****************************************************************************/

enum insertsiteresult insertsite(TriPoint insertpoint, TriOrientedTriangle *searchtri, TriOrientedShell *splitedge,
                                 int segmentflaws, int triflaws)
{
  TriOrientedTriangle horiz;
  TriOrientedTriangle top;
  TriOrientedTriangle botleft, botright;
  TriOrientedTriangle topleft, topright;
  TriOrientedTriangle newbotleft, newbotright;
  TriOrientedTriangle newtopright;
  TriOrientedTriangle botlcasing, botrcasing;
  TriOrientedTriangle toplcasing, toprcasing;
  TriOrientedTriangle testtri;
  TriOrientedShell botlshelle, botrshelle;
  TriOrientedShell toplshelle, toprshelle;
  TriOrientedShell brokenshelle;
  TriOrientedShell checkshelle;
  TriOrientedShell rightedge;
  TriOrientedShell newedge;
  TriOrientedShell *encroached;
  TriPoint first;
  TriPoint leftpoint, rightpoint, botpoint, toppoint, farpoint;
  REAL attrib;
  REAL area;
  enum insertsiteresult success;
  TriLocateResultType intersect;
  int doflip;
  int mirrorflag;
  int i;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;         /* Temporary variable used by spivot() and tspivot(). */

  if (verbose > 1) {
    printf("  Inserting (%.12g, %.12g).\n", insertpoint[0], insertpoint[1]);
  }
  if (splitedge == (TriOrientedShell *) NULL) {
    /* Find the location of the TriPoint to be inserted.  Check if a good */
    /*   starting TriTriangle has already been provided by the caller.    */
    if (searchtri->tri == (TriTriangle *) NULL) {
      /* Find a boundary TriTriangle. */
      horiz.tri = dummytri;
      horiz.orient = 0;
      symself(horiz);
      /* Search for a TriTriangle containing `insertpoint'. */
      intersect = locate(insertpoint, &horiz);
    } else {
      /* Start searching from the TriTriangle provided by the caller. */
      TriOrientedTrianglecopy(*searchtri, horiz);
      intersect = preciselocate(insertpoint, &horiz);
    }
  } else {
    /* The calling routine provides the edge in which the TriPoint is inserted. */
    TriOrientedTrianglecopy(*searchtri, horiz);
    intersect = ONEDGE;
  }
  if (intersect == ONVERTEX) {
    /* There's already a vertex there.  Return in `searchtri' a TriTriangle */
    /*   whose origin is the existing vertex.                            */
    TriOrientedTrianglecopy(horiz, *searchtri);
    TriOrientedTrianglecopy(horiz, recenttri);
    return DUPLICATEPOINT;
  }
  if ((intersect == ONEDGE) || (intersect == OUTSIDE)) {
    /* The vertex falls on an edge or boundary. */
    if (checksegments && (splitedge == (TriOrientedShell *) NULL)) {
      /* Check whether the vertex falls on a shell edge. */
      tspivot(horiz, brokenshelle);
      if (brokenshelle.sh != dummysh) {
        /* The vertex falls on a shell edge. */
        if (segmentflaws) {
          if (nobisect == 0) {
            /* Add the shell edge to the list of encroached segments. */
            encroached = (TriOrientedShell *) poolalloc(&badsegments);
            shellecopy(brokenshelle, *encroached);
          } else if ((nobisect == 1) && (intersect == ONEDGE)) {
            /* This segment may be split only if it is an internal boundary. */
            sym(horiz, testtri);
            if (testtri.tri != dummytri) {
              /* Add the shell edge to the list of encroached segments. */
              encroached = (TriOrientedShell *) poolalloc(&badsegments);
              shellecopy(brokenshelle, *encroached);
            }
          }
        }
        /* Return a handle whose primary edge contains the TriPoint, */
        /*   which has not been inserted.                         */
        TriOrientedTrianglecopy(horiz, *searchtri);
        TriOrientedTrianglecopy(horiz, recenttri);
        return VIOLATINGPOINT;
      }
    }
    /* Insert the TriPoint on an edge, dividing one TriTriangle into two (if */
    /*   the edge lies on a boundary) or two triangles into four.      */
    lprev(horiz, botright);
    sym(botright, botrcasing);
    sym(horiz, topright);
    /* Is there a second TriTriangle?  (Or does this edge lie on a boundary?) */
    mirrorflag = topright.tri != dummytri;
    if (mirrorflag) {
      lnextself(topright);
      sym(topright, toprcasing);
      maketriangle(&newtopright);
    } else {
      /* Splitting the boundary edge increases the number of boundary edges. */
      hullsize++;
    }
    maketriangle(&newbotright);

    /* Set the vertices of changed and new triangles. */
    org(horiz, rightpoint);
    dest(horiz, leftpoint);
    apex(horiz, botpoint);
    setorg(newbotright, botpoint);
    setdest(newbotright, rightpoint);
    setapex(newbotright, insertpoint);
    setorg(horiz, insertpoint);
    for (i = 0; i < eextras; i++) {
      /* Set the element attributes of a new TriTriangle. */
      setelemattribute(newbotright, i, elemattribute(botright, i));
    }
    if (vararea) {
      /* Set the area constraint of a new TriTriangle. */
      setareabound(newbotright, areabound(botright));
    }
    if (mirrorflag) {
      dest(topright, toppoint);
      setorg(newtopright, rightpoint);
      setdest(newtopright, toppoint);
      setapex(newtopright, insertpoint);
      setorg(topright, insertpoint);
      for (i = 0; i < eextras; i++) {
        /* Set the element attributes of another new TriTriangle. */
        setelemattribute(newtopright, i, elemattribute(topright, i));
      }
      if (vararea) {
        /* Set the area constraint of another new TriTriangle. */
        setareabound(newtopright, areabound(topright));
      }
    }

    /* There may be shell edges that need to be bonded */
    /*   to the new TriTriangle(s).                       */
    if (checksegments) {
      tspivot(botright, botrshelle);
      if (botrshelle.sh != dummysh) {
        tsdissolve(botright);
        tsbond(newbotright, botrshelle);
      }
      if (mirrorflag) {
        tspivot(topright, toprshelle);
        if (toprshelle.sh != dummysh) {
          tsdissolve(topright);
          tsbond(newtopright, toprshelle);
        }
      }
    }

    /* Bond the new TriTriangle(s) to the surrounding triangles. */
    bond(newbotright, botrcasing);
    lprevself(newbotright);
    bond(newbotright, botright);
    lprevself(newbotright);
    if (mirrorflag) {
      bond(newtopright, toprcasing);
      lnextself(newtopright);
      bond(newtopright, topright);
      lnextself(newtopright);
      bond(newtopright, newbotright);
    }

    if (splitedge != (TriOrientedShell *) NULL) {
      /* Split the shell edge into two. */
      setsdest(*splitedge, insertpoint);
      ssymself(*splitedge);
      spivot(*splitedge, rightedge);
      insertshelle(&newbotright, mark(*splitedge));
      tspivot(newbotright, newedge);
      sbond(*splitedge, newedge);
      ssymself(newedge);
      sbond(newedge, rightedge);
      ssymself(*splitedge);
    }

#ifdef SELF_CHECK
    if (counterclockwise(rightpoint, leftpoint, botpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle prior to edge TriPoint insertion (bottom).\n");
    }
    if (mirrorflag) {
      if (counterclockwise(leftpoint, rightpoint, toppoint) < 0.0) {
        printf("Internal error in insertsite():\n");
        printf("  Clockwise TriTriangle prior to edge TriPoint insertion (top).\n");
      }
      if (counterclockwise(rightpoint, toppoint, insertpoint) < 0.0) {
        printf("Internal error in insertsite():\n");
        printf("  Clockwise TriTriangle after edge TriPoint insertion (top right).\n"
               );
      }
      if (counterclockwise(toppoint, leftpoint, insertpoint) < 0.0) {
        printf("Internal error in insertsite():\n");
        printf("  Clockwise TriTriangle after edge TriPoint insertion (top left).\n"
               );
      }
    }
    if (counterclockwise(leftpoint, botpoint, insertpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle after edge TriPoint insertion (bottom left).\n"
             );
    }
    if (counterclockwise(botpoint, rightpoint, insertpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf(
        "  Clockwise TriTriangle after edge TriPoint insertion (bottom right).\n");
    }
#endif /* SELF_CHECK */
    if (verbose > 2) {
      printf("  Updating bottom left ");
      printtriangle(&botright);
      if (mirrorflag) {
        printf("  Updating top left ");
        printtriangle(&topright);
        printf("  Creating top right ");
        printtriangle(&newtopright);
      }
      printf("  Creating bottom right ");
      printtriangle(&newbotright);
    }

    /* Position `horiz' on the first edge to check for */
    /*   the Delaunay property.                        */
    lnextself(horiz);
  } else {
    /* Insert the TriPoint in a TriTriangle, splitting it into three. */
    lnext(horiz, botleft);
    lprev(horiz, botright);
    sym(botleft, botlcasing);
    sym(botright, botrcasing);
    maketriangle(&newbotleft);
    maketriangle(&newbotright);

    /* Set the vertices of changed and new triangles. */
    org(horiz, rightpoint);
    dest(horiz, leftpoint);
    apex(horiz, botpoint);
    setorg(newbotleft, leftpoint);
    setdest(newbotleft, botpoint);
    setapex(newbotleft, insertpoint);
    setorg(newbotright, botpoint);
    setdest(newbotright, rightpoint);
    setapex(newbotright, insertpoint);
    setapex(horiz, insertpoint);
    for (i = 0; i < eextras; i++) {
      /* Set the element attributes of the new triangles. */
      attrib = elemattribute(horiz, i);
      setelemattribute(newbotleft, i, attrib);
      setelemattribute(newbotright, i, attrib);
    }
    if (vararea) {
      /* Set the area constraint of the new triangles. */
      area = areabound(horiz);
      setareabound(newbotleft, area);
      setareabound(newbotright, area);
    }

    /* There may be shell edges that need to be bonded */
    /*   to the new triangles.                         */
    if (checksegments) {
      tspivot(botleft, botlshelle);
      if (botlshelle.sh != dummysh) {
        tsdissolve(botleft);
        tsbond(newbotleft, botlshelle);
      }
      tspivot(botright, botrshelle);
      if (botrshelle.sh != dummysh) {
        tsdissolve(botright);
        tsbond(newbotright, botrshelle);
      }
    }

    /* Bond the new triangles to the surrounding triangles. */
    bond(newbotleft, botlcasing);
    bond(newbotright, botrcasing);
    lnextself(newbotleft);
    lprevself(newbotright);
    bond(newbotleft, newbotright);
    lnextself(newbotleft);
    bond(botleft, newbotleft);
    lprevself(newbotright);
    bond(botright, newbotright);

#ifdef SELF_CHECK
    if (counterclockwise(rightpoint, leftpoint, botpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle prior to TriPoint insertion.\n");
    }
    if (counterclockwise(rightpoint, leftpoint, insertpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle after TriPoint insertion (top).\n");
    }
    if (counterclockwise(leftpoint, botpoint, insertpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle after TriPoint insertion (left).\n");
    }
    if (counterclockwise(botpoint, rightpoint, insertpoint) < 0.0) {
      printf("Internal error in insertsite():\n");
      printf("  Clockwise TriTriangle after TriPoint insertion (right).\n");
    }
#endif /* SELF_CHECK */
    if (verbose > 2) {
      printf("  Updating top ");
      printtriangle(&horiz);
      printf("  Creating left ");
      printtriangle(&newbotleft);
      printf("  Creating right ");
      printtriangle(&newbotright);
    }
  }

  /* The insertion is successful by default, unless an encroached */
  /*   edge is found.                                             */
  success = SUCCESSFULPOINT;
  /* Circle around the newly inserted vertex, checking each edge opposite */
  /*   it for the Delaunay property.  Non-Delaunay edges are flipped.     */
  /*   `horiz' is always the edge being checked.  `first' marks where to  */
  /*   stop circling.                                                     */
  org(horiz, first);
  rightpoint = first;
  dest(horiz, leftpoint);
  /* Circle until finished. */
  while (1) {
    /* By default, the edge will be flipped. */
    doflip = 1;
    if (checksegments) {
      /* Check for a segment, which cannot be flipped. */
      tspivot(horiz, checkshelle);
      if (checkshelle.sh != dummysh) {
        /* The edge is a segment and cannot be flipped. */
        doflip = 0;
      }
    }
    if (doflip) {
      /* Check if the edge is a boundary edge. */
      sym(horiz, top);
      if (top.tri == dummytri) {
        /* The edge is a boundary edge and cannot be flipped. */
        doflip = 0;
      } else {
        /* Find the TriPoint on the other side of the edge. */
        apex(top, farpoint);
        /* In the incremental Delaunay triangulation algorithm, any of    */
        /*   `leftpoint', `rightpoint', and `farpoint' could be vertices  */
        /*   of the triangular bounding box.  These vertices must be      */
        /*   treated as if they are infinitely distant, even though their */
        /*   "coordinates" are not.                                       */
        if ((leftpoint == infpoint1) || (leftpoint == infpoint2)
                   || (leftpoint == infpoint3)) {
          /* `leftpoint' is infinitely distant.  Check the convexity of */
          /*   the boundary of the triangulation.  'farpoint' might be  */
          /*   infinite as well, but trust me, this same condition      */
          /*   should be applied.                                       */
          doflip = counterclockwise(insertpoint, rightpoint, farpoint) > 0.0;
        } else if ((rightpoint == infpoint1) || (rightpoint == infpoint2)
                   || (rightpoint == infpoint3)) {
          /* `rightpoint' is infinitely distant.  Check the convexity of */
          /*   the boundary of the triangulation.  'farpoint' might be  */
          /*   infinite as well, but trust me, this same condition      */
          /*   should be applied.                                       */
          doflip = counterclockwise(farpoint, leftpoint, insertpoint) > 0.0;
        } else if ((farpoint == infpoint1) || (farpoint == infpoint2)
            || (farpoint == infpoint3)) {
          /* `farpoint' is infinitely distant and cannot be inside */
          /*   the circumcircle of the TriTriangle `horiz'.           */
          doflip = 0;
        } else {
          /* Test whether the edge is locally Delaunay. */
          doflip = incircle(leftpoint, insertpoint, rightpoint, farpoint)
                   > 0.0;
        }
        if (doflip) {
          /* We made it!  Flip the edge `horiz' by rotating its containing */
          /*   quadrilateral (the two triangles adjacent to `horiz').      */
          /* Identify the casing of the quadrilateral. */
          lprev(top, topleft);
          sym(topleft, toplcasing);
          lnext(top, topright);
          sym(topright, toprcasing);
          lnext(horiz, botleft);
          sym(botleft, botlcasing);
          lprev(horiz, botright);
          sym(botright, botrcasing);
          /* Rotate the quadrilateral one-quarter turn counterclockwise. */
          bond(topleft, botlcasing);
          bond(botleft, botrcasing);
          bond(botright, toprcasing);
          bond(topright, toplcasing);
          if (checksegments) {
            /* Check for shell edges and rebond them to the quadrilateral. */
            tspivot(topleft, toplshelle);
            tspivot(botleft, botlshelle);
            tspivot(botright, botrshelle);
            tspivot(topright, toprshelle);
            if (toplshelle.sh == dummysh) {
              tsdissolve(topright);
            } else {
              tsbond(topright, toplshelle);
            }
            if (botlshelle.sh == dummysh) {
              tsdissolve(topleft);
            } else {
              tsbond(topleft, botlshelle);
            }
            if (botrshelle.sh == dummysh) {
              tsdissolve(botleft);
            } else {
              tsbond(botleft, botrshelle);
            }
            if (toprshelle.sh == dummysh) {
              tsdissolve(botright);
            } else {
              tsbond(botright, toprshelle);
            }
          }
          /* New TriPoint assignments for the rotated quadrilateral. */
          setorg(horiz, farpoint);
          setdest(horiz, insertpoint);
          setapex(horiz, rightpoint);
          setorg(top, insertpoint);
          setdest(top, farpoint);
          setapex(top, leftpoint);
          for (i = 0; i < eextras; i++) {
            /* Take the average of the two triangles' attributes. */
            attrib = 0.5 * (elemattribute(top, i) + elemattribute(horiz, i));
            setelemattribute(top, i, attrib);
            setelemattribute(horiz, i, attrib);
          }
          if (vararea) {
            if ((areabound(top) <= 0.0) || (areabound(horiz) <= 0.0)) {
              area = -1.0;
            } else {
              /* Take the average of the two triangles' area constraints.    */
              /*   This prevents small area constraints from migrating a     */
              /*   long, long way from their original location due to flips. */
              area = 0.5 * (areabound(top) + areabound(horiz));
            }
            setareabound(top, area);
            setareabound(horiz, area);
          }
#ifdef SELF_CHECK
          if (insertpoint != (TriPoint) NULL) {
            if (counterclockwise(leftpoint, insertpoint, rightpoint) < 0.0) {
              printf("Internal error in insertsite():\n");
              printf("  Clockwise TriTriangle prior to edge flip (bottom).\n");
            }
            /* The following test has been removed because constrainededge() */
            /*   sometimes generates inverted triangles that insertsite()    */
            /*   removes.                                                    */
/*
            if (counterclockwise(rightpoint, farpoint, leftpoint) < 0.0) {
              printf("Internal error in insertsite():\n");
              printf("  Clockwise TriTriangle prior to edge flip (top).\n");
            }
*/
            if (counterclockwise(farpoint, leftpoint, insertpoint) < 0.0) {
              printf("Internal error in insertsite():\n");
              printf("  Clockwise TriTriangle after edge flip (left).\n");
            }
            if (counterclockwise(insertpoint, rightpoint, farpoint) < 0.0) {
              printf("Internal error in insertsite():\n");
              printf("  Clockwise TriTriangle after edge flip (right).\n");
            }
          }
#endif /* SELF_CHECK */
          if (verbose > 2) {
            printf("  Edge flip results in left ");
            lnextself(topleft);
            printtriangle(&topleft);
            printf("  and right ");
            printtriangle(&horiz);
          }
          /* On the next iterations, consider the two edges that were  */
          /*   exposed (this is, are now visible to the newly inserted */
          /*   TriPoint) by the edge flip.                                */
          lprevself(horiz);
          leftpoint = farpoint;
        }
      }
    }
    if (!doflip) {
      /* The handle `horiz' is accepted as locally Delaunay. */
      /* Look for the next edge around the newly inserted TriPoint. */
      lnextself(horiz);
      sym(horiz, testtri);
      /* Check for finishing a complete revolution about the new TriPoint, or */
      /*   falling off the edge of the triangulation.  The latter will     */
      /*   happen when a TriPoint is inserted at a boundary.                  */
      if ((leftpoint == first) || (testtri.tri == dummytri)) {
        /* We're done.  Return a TriTriangle whose origin is the new TriPoint. */
        lnext(horiz, *searchtri);
        lnext(horiz, recenttri);
        return success;
      }
      /* Finish finding the next edge around the newly inserted TriPoint. */
      lnext(testtri, horiz);
      rightpoint = leftpoint;
      dest(horiz, leftpoint);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  triangulatepolygon()   Find the Delaunay triangulation of a polygon that */
/*                         has a certain "nice" shape.  This includes the    */
/*                         polygons that result from deletion of a TriPoint or  */
/*                         insertion of a segment.                           */
/*                                                                           */
/*  This is a conceptually difficult routine.  The starting assumption is    */
/*  that we have a polygon with n sides.  n - 1 of these sides are currently */
/*  represented as edges in the mesh.  One side, called the "base", need not */
/*  be.                                                                      */
/*                                                                           */
/*  Inside the polygon is a structure I call a "fan", consisting of n - 1    */
/*  triangles that share a common origin.  For each of these triangles, the  */
/*  edge opposite the origin is one of the sides of the polygon.  The        */
/*  primary edge of each TriTriangle is the edge directed from the origin to    */
/*  the destination; note that this is not the same edge that is a side of   */
/*  the polygon.  `firstedge' is the primary edge of the first TriTriangle.     */
/*  From there, the triangles follow in counterclockwise order about the     */
/*  polygon, until `lastedge', the primary edge of the last TriTriangle.        */
/*  `firstedge' and `lastedge' are probably connected to other triangles     */
/*  beyond the extremes of the fan, but their identity is not important, as  */
/*  long as the fan remains connected to them.                               */
/*                                                                           */
/*  Imagine the polygon oriented so that its base is at the bottom.  This    */
/*  puts `firstedge' on the far right, and `lastedge' on the far left.       */
/*  The right vertex of the base is the destination of `firstedge', and the  */
/*  left vertex of the base is the apex of `lastedge'.                       */
/*                                                                           */
/*  The challenge now is to find the right sequence of edge flips to         */
/*  transform the fan into a Delaunay triangulation of the polygon.  Each    */
/*  edge flip effectively removes one TriTriangle from the fan, committing it   */
/*  to the polygon.  The resulting polygon has one fewer edge.  If `doflip'  */
/*  is set, the final flip will be performed, resulting in a fan of one      */
/*  (useless?) TriTriangle.  If `doflip' is not set, the final flip is not      */
/*  performed, resulting in a fan of two triangles, and an unfinished        */
/*  triangular polygon that is not yet filled out with a single TriTriangle.    */
/*  On completion of the routine, `lastedge' is the last remaining TriTriangle, */
/*  or the leftmost of the last two.                                         */
/*                                                                           */
/*  Although the flips are performed in the order described above, the       */
/*  decisions about what flips to perform are made in precisely the reverse  */
/*  order.  The recursive triangulatepolygon() procedure makes a decision,   */
/*  uses up to two recursive calls to triangulate the "subproblems"          */
/*  (polygons with fewer edges), and then performs an edge flip.             */
/*                                                                           */
/*  The "decision" it makes is which vertex of the polygon should be         */
/*  connected to the base.  This decision is made by testing every possible  */
/*  vertex.  Once the best vertex is found, the two edges that connect this  */
/*  vertex to the base become the bases for two smaller polygons.  These     */
/*  are triangulated recursively.  Unfortunately, this approach can take     */
/*  O(n^2) time not only in the worst case, but in many common cases.  It's  */
/*  rarely a big deal for TriPoint deletion, where n is rarely larger than ten, */
/*  but it could be a big deal for segment insertion, especially if there's  */
/*  a lot of long segments that each cut many triangles.  I ought to code    */
/*  a faster algorithm some time.                                            */
/*                                                                           */
/*  The `edgecount' parameter is the number of sides of the polygon,         */
/*  including its base.  `triflaws' is a flag that determines whether the    */
/*  new triangles should be tested for quality, and enqueued if they are     */
/*  bad.                                                                     */
/*                                                                           */
/*****************************************************************************/

void triangulatepolygon(TriOrientedTriangle *firstedge, TriOrientedTriangle *lastedge, 
						int edgecount, int doflip, int triflaws)
{
  TriOrientedTriangle testtri;
  TriOrientedTriangle besttri;
  TriOrientedTriangle tempedge;
  TriPoint leftbasepoint, rightbasepoint;
  TriPoint testpoint;
  TriPoint bestpoint;
  int bestnumber;
  int i;
  TriTriangle ptr;   /* Temporary variable used by sym(), onext(), and oprev(). */

  /* Identify the base vertices. */
  apex(*lastedge, leftbasepoint);
  dest(*firstedge, rightbasepoint);
  if (verbose > 2) {
    printf("  Triangulating interior polygon at edge\n");
    printf("    (%.12g, %.12g) (%.12g, %.12g)\n", leftbasepoint[0],
           leftbasepoint[1], rightbasepoint[0], rightbasepoint[1]);
  }
  /* Find the best vertex to connect the base to. */
  onext(*firstedge, besttri);
  dest(besttri, bestpoint);
  TriOrientedTrianglecopy(besttri, testtri);
  bestnumber = 1;
  for (i = 2; i <= edgecount - 2; i++) {
    onextself(testtri);
    dest(testtri, testpoint);
    /* Is this a better vertex? */
    if (incircle(leftbasepoint, rightbasepoint, bestpoint, testpoint) > 0.0) {
      TriOrientedTrianglecopy(testtri, besttri);
      bestpoint = testpoint;
      bestnumber = i;
    }
  }
  if (verbose > 2) {
    printf("    Connecting edge to (%.12g, %.12g)\n", bestpoint[0],
           bestpoint[1]);
  }
  if (bestnumber > 1) {
    /* Recursively triangulate the smaller polygon on the right. */
    oprev(besttri, tempedge);
    triangulatepolygon(firstedge, &tempedge, bestnumber + 1, 1, triflaws);
  }
  if (bestnumber < edgecount - 2) {
    /* Recursively triangulate the smaller polygon on the left. */
    sym(besttri, tempedge);
    triangulatepolygon(&besttri, lastedge, edgecount - bestnumber, 1,
                       triflaws);
    /* Find `besttri' again; it may have been lost to edge flips. */
    sym(tempedge, besttri);
  }
  if (doflip) {
    /* Do one final edge flip. */
    flip(&besttri);
  }
  /* Return the base TriTriangle. */
  TriOrientedTrianglecopy(besttri, *lastedge);
}


/**                                                                         **/
/**                                                                         **/
/********* Mesh transformation routines end here                     *********/

/********* Divide-and-conquer Delaunay triangulation begins here     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  The divide-and-conquer bounding box                                      */
/*                                                                           */
/*  I originally implemented the divide-and-conquer and incremental Delaunay */
/*  triangulations using the edge-based data structure presented by Guibas   */
/*  and Stolfi.  Switching to a TriTriangle-based data structure doubled the    */
/*  speed.  However, I had to think of a few extra tricks to maintain the    */
/*  elegance of the original algorithms.                                     */
/*                                                                           */
/*  The "bounding box" used by my variant of the divide-and-conquer          */
/*  algorithm uses one TriTriangle for each edge of the convex hull of the      */
/*  triangulation.  These bounding triangles all share a common apical       */
/*  vertex, which is represented by NULL and which represents nothing.       */
/*  The bounding triangles are linked in a circular fan about this NULL      */
/*  vertex, and the edges on the convex hull of the triangulation appear     */
/*  opposite the NULL vertex.  You might find it easiest to imagine that     */
/*  the NULL vertex is a TriPoint in 3D space behind the center of the          */
/*  triangulation, and that the bounding triangles form a sort of cone.      */
/*                                                                           */
/*  This bounding box makes it easy to represent degenerate cases.  For      */
/*  instance, the triangulation of two vertices is a single edge.  This edge */
/*  is represented by two bounding box triangles, one on each "side" of the  */
/*  edge.  These triangles are also linked together in a fan about the NULL  */
/*  vertex.                                                                  */
/*                                                                           */
/*  The bounding box also makes it easy to traverse the convex hull, as the  */
/*  divide-and-conquer algorithm needs to do.                                */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  pointsort()   Sort an array of points by x-coordinate, using the         */
/*                y-coordinate as a secondary key.                           */
/*                                                                           */
/*  Uses quicksort.  Randomized O(n log n) time.  No, I did not make any of  */
/*  the usual quicksort mistakes.                                            */
/*                                                                           */
/*****************************************************************************/

void pointsort(TriPoint *sortarray, int arraysize)
{
  int left, right;
  int pivot;
  REAL pivotx, pivoty;
  TriPoint temp;

  if (arraysize == 2) {
    /* Recursive base case. */
    if ((sortarray[0][0] > sortarray[1][0]) ||
        ((sortarray[0][0] == sortarray[1][0]) &&
         (sortarray[0][1] > sortarray[1][1]))) {
      temp = sortarray[1];
      sortarray[1] = sortarray[0];
      sortarray[0] = temp;
    }
    return;
  }
  /* Choose a random pivot to split the array. */
  pivot = (int) randomnation(arraysize);
  pivotx = sortarray[pivot][0];
  pivoty = sortarray[pivot][1];
  /* Split the array. */
  left = -1;
  right = arraysize;
  while (left < right) {
    /* Search for a TriPoint whose x-coordinate is too large for the left. */
    do {
      left++;
    } while ((left <= right) && ((sortarray[left][0] < pivotx) ||
                                 ((sortarray[left][0] == pivotx) &&
                                  (sortarray[left][1] < pivoty))));
    /* Search for a TriPoint whose x-coordinate is too small for the right. */
    do {
      right--;
    } while ((left <= right) && ((sortarray[right][0] > pivotx) ||
                                 ((sortarray[right][0] == pivotx) &&
                                  (sortarray[right][1] > pivoty))));
    if (left < right) {
      /* Swap the left and right points. */
      temp = sortarray[left];
      sortarray[left] = sortarray[right];
      sortarray[right] = temp;
    }
  }
  if (left > 1) {
    /* Recursively sort the left subset. */
    pointsort(sortarray, left);
  }
  if (right < arraysize - 2) {
    /* Recursively sort the right subset. */
    pointsort(&sortarray[right + 1], arraysize - right - 1);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  pointmedian()   An order statistic algorithm, almost.  Shuffles an array */
/*                  of points so that the first `median' points occur        */
/*                  lexicographically before the remaining points.           */
/*                                                                           */
/*  Uses the x-coordinate as the primary key if axis == 0; the y-coordinate  */
/*  if axis == 1.  Very similar to the pointsort() procedure, but runs in    */
/*  randomized linear time.                                                  */
/*                                                                           */
/*****************************************************************************/

void pointmedian(TriPoint *sortarray, int arraysize, int median, int axis)
{
  int left, right;
  int pivot;
  REAL pivot1, pivot2;
  TriPoint temp;

  if (arraysize == 2) {
    /* Recursive base case. */
    if ((sortarray[0][axis] > sortarray[1][axis]) ||
        ((sortarray[0][axis] == sortarray[1][axis]) &&
         (sortarray[0][1 - axis] > sortarray[1][1 - axis]))) {
      temp = sortarray[1];
      sortarray[1] = sortarray[0];
      sortarray[0] = temp;
    }
    return;
  }
  /* Choose a random pivot to split the array. */
  pivot = (int) randomnation(arraysize);
  pivot1 = sortarray[pivot][axis];
  pivot2 = sortarray[pivot][1 - axis];
  /* Split the array. */
  left = -1;
  right = arraysize;
  while (left < right) {
    /* Search for a TriPoint whose x-coordinate is too large for the left. */
    do {
      left++;
    } while ((left <= right) && ((sortarray[left][axis] < pivot1) ||
                                 ((sortarray[left][axis] == pivot1) &&
                                  (sortarray[left][1 - axis] < pivot2))));
    /* Search for a TriPoint whose x-coordinate is too small for the right. */
    do {
      right--;
    } while ((left <= right) && ((sortarray[right][axis] > pivot1) ||
                                 ((sortarray[right][axis] == pivot1) &&
                                  (sortarray[right][1 - axis] > pivot2))));
    if (left < right) {
      /* Swap the left and right points. */
      temp = sortarray[left];
      sortarray[left] = sortarray[right];
      sortarray[right] = temp;
    }
  }
  /* Unlike in pointsort(), at most one of the following */
  /*   conditionals is true.                             */
  if (left > median) {
    /* Recursively shuffle the left subset. */
    pointmedian(sortarray, left, median, axis);
  }
  if (right < median - 1) {
    /* Recursively shuffle the right subset. */
    pointmedian(&sortarray[right + 1], arraysize - right - 1,
                median - right - 1, axis);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  alternateaxes()   Sorts the points as appropriate for the divide-and-    */
/*                    conquer algorithm with alternating cuts.               */
/*                                                                           */
/*  Partitions by x-coordinate if axis == 0; by y-coordinate if axis == 1.   */
/*  For the base case, subsets containing only two or three points are       */
/*  always sorted by x-coordinate.                                           */
/*                                                                           */
/*****************************************************************************/

void alternateaxes(TriPoint *sortarray, int arraysize, int axis)
{
  int divider;

  divider = arraysize >> 1;
  if (arraysize <= 3) {
    /* Recursive base case:  subsets of two or three points will be      */
    /*   handled specially, and should always be sorted by x-coordinate. */
    axis = 0;
  }
  /* Partition with a horizontal or vertical cut. */
  pointmedian(sortarray, arraysize, divider, axis);
  /* Recursively partition the subsets with a cross cut. */
  if (arraysize - divider >= 2) {
    if (divider >= 2) {
      alternateaxes(sortarray, divider, 1 - axis);
    }
    alternateaxes(&sortarray[divider], arraysize - divider, 1 - axis);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  mergehulls()   Merge two adjacent Delaunay triangulations into a         */
/*                 single Delaunay triangulation.                            */
/*                                                                           */
/*  This is similar to the algorithm given by Guibas and Stolfi, but uses    */
/*  a TriTriangle-based, rather than edge-based, data structure.                */
/*                                                                           */
/*  The algorithm walks up the gap between the two triangulations, knitting  */
/*  them together.  As they are merged, some of their bounding triangles     */
/*  are converted into real triangles of the triangulation.  The procedure   */
/*  pulls each hull's bounding triangles apart, then knits them together     */
/*  like the teeth of two gears.  The Delaunay property determines, at each  */
/*  step, whether the next "tooth" is a bounding TriTriangle of the left hull   */
/*  or the right.  When a bounding TriTriangle becomes real, its apex is        */
/*  changed from NULL to a real TriPoint.                                       */
/*                                                                           */
/*  Only two new triangles need to be allocated.  These become new bounding  */
/*  triangles at the top and bottom of the seam.  They are used to connect   */
/*  the remaining bounding triangles (those that have not been converted     */
/*  into real triangles) into a single fan.                                  */
/*                                                                           */
/*  On entry, `farleft' and `innerleft' are bounding triangles of the left   */
/*  triangulation.  The origin of `farleft' is the leftmost vertex, and      */
/*  the destination of `innerleft' is the rightmost vertex of the            */
/*  triangulation.  Similarly, `innerright' and `farright' are bounding      */
/*  triangles of the right triangulation.  The origin of `innerright' and    */
/*  destination of `farright' are the leftmost and rightmost vertices.       */
/*                                                                           */
/*  On completion, the origin of `farleft' is the leftmost vertex of the     */
/*  merged triangulation, and the destination of `farright' is the rightmost */
/*  vertex.                                                                  */
/*                                                                           */
/*****************************************************************************/

void mergehulls(TriOrientedTriangle *farleft, TriOrientedTriangle *innerleft, 
				TriOrientedTriangle *innerright, TriOrientedTriangle *farright, int axis)
{
  TriOrientedTriangle leftcand, rightcand;
  TriOrientedTriangle baseedge;
  TriOrientedTriangle nextedge;
  TriOrientedTriangle sidecasing, topcasing, outercasing;
  TriOrientedTriangle checkedge;
  TriPoint innerleftdest;
  TriPoint innerrightorg;
  TriPoint innerleftapex, innerrightapex;
  TriPoint farleftpt, farrightpt;
  TriPoint farleftapex, farrightapex;
  TriPoint lowerleft, lowerright;
  TriPoint upperleft, upperright;
  TriPoint nextapex;
  TriPoint checkvertex;
  int changemade;
  int badedge;
  int leftfinished, rightfinished;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  dest(*innerleft, innerleftdest);
  apex(*innerleft, innerleftapex);
  org(*innerright, innerrightorg);
  apex(*innerright, innerrightapex);
  /* Special treatment for horizontal cuts. */
  if (dwyer && (axis == 1)) {
    org(*farleft, farleftpt);
    apex(*farleft, farleftapex);
    dest(*farright, farrightpt);
    apex(*farright, farrightapex);
    /* The pointers to the extremal points are shifted to TriPoint to the */
    /*   topmost and bottommost TriPoint of each hull, rather than the    */
    /*   leftmost and rightmost points.                                */
    while (farleftapex[1] < farleftpt[1]) {
      lnextself(*farleft);
      symself(*farleft);
      farleftpt = farleftapex;
      apex(*farleft, farleftapex);
    }
    sym(*innerleft, checkedge);
    apex(checkedge, checkvertex);
    while (checkvertex[1] > innerleftdest[1]) {
      lnext(checkedge, *innerleft);
      innerleftapex = innerleftdest;
      innerleftdest = checkvertex;
      sym(*innerleft, checkedge);
      apex(checkedge, checkvertex);
    }
    while (innerrightapex[1] < innerrightorg[1]) {
      lnextself(*innerright);
      symself(*innerright);
      innerrightorg = innerrightapex;
      apex(*innerright, innerrightapex);
    }
    sym(*farright, checkedge);
    apex(checkedge, checkvertex);
    while (checkvertex[1] > farrightpt[1]) {
      lnext(checkedge, *farright);
      farrightapex = farrightpt;
      farrightpt = checkvertex;
      sym(*farright, checkedge);
      apex(checkedge, checkvertex);
    }
  }
  /* Find a line tangent to and below both hulls. */
  do {
    changemade = 0;
    /* Make innerleftdest the "bottommost" TriPoint of the left hull. */
    if (counterclockwise(innerleftdest, innerleftapex, innerrightorg) > 0.0) {
      lprevself(*innerleft);
      symself(*innerleft);
      innerleftdest = innerleftapex;
      apex(*innerleft, innerleftapex);
      changemade = 1;
    }
    /* Make innerrightorg the "bottommost" TriPoint of the right hull. */
    if (counterclockwise(innerrightapex, innerrightorg, innerleftdest) > 0.0) {
      lnextself(*innerright);
      symself(*innerright);
      innerrightorg = innerrightapex;
      apex(*innerright, innerrightapex);
      changemade = 1;
    }
  } while (changemade);
  /* Find the two candidates to be the next "gear tooth". */
  sym(*innerleft, leftcand);
  sym(*innerright, rightcand);
  /* Create the bottom new bounding TriTriangle. */
  maketriangle(&baseedge);
  /* Connect it to the bounding boxes of the left and right triangulations. */
  bond(baseedge, *innerleft);
  lnextself(baseedge);
  bond(baseedge, *innerright);
  lnextself(baseedge);
  setorg(baseedge, innerrightorg);
  setdest(baseedge, innerleftdest);
  /* Apex is intentionally left NULL. */
  if (verbose > 2) {
    printf("  Creating base bounding ");
    printtriangle(&baseedge);
  }
  /* Fix the extreme triangles if necessary. */
  org(*farleft, farleftpt);
  if (innerleftdest == farleftpt) {
    lnext(baseedge, *farleft);
  }
  dest(*farright, farrightpt);
  if (innerrightorg == farrightpt) {
    lprev(baseedge, *farright);
  }
  /* The vertices of the current knitting edge. */
  lowerleft = innerleftdest;
  lowerright = innerrightorg;
  /* The candidate vertices for knitting. */
  apex(leftcand, upperleft);
  apex(rightcand, upperright);
  /* Walk up the gap between the two triangulations, knitting them together. */
  while (1) {
    /* Have we reached the top?  (This isn't quite the right question,       */
    /*   because even though the left triangulation might seem finished now, */
    /*   moving up on the right triangulation might reveal a new TriPoint of    */
    /*   the left triangulation.  And vice-versa.)                           */
    leftfinished = counterclockwise(upperleft, lowerleft, lowerright) <= 0.0;
    rightfinished = counterclockwise(upperright, lowerleft, lowerright) <= 0.0;
    if (leftfinished && rightfinished) {
      /* Create the top new bounding TriTriangle. */
      maketriangle(&nextedge);
      setorg(nextedge, lowerleft);
      setdest(nextedge, lowerright);
      /* Apex is intentionally left NULL. */
      /* Connect it to the bounding boxes of the two triangulations. */
      bond(nextedge, baseedge);
      lnextself(nextedge);
      bond(nextedge, rightcand);
      lnextself(nextedge);
      bond(nextedge, leftcand);
      if (verbose > 2) {
        printf("  Creating top bounding ");
        printtriangle(&baseedge);
      }
      /* Special treatment for horizontal cuts. */
      if (dwyer && (axis == 1)) {
        org(*farleft, farleftpt);
        apex(*farleft, farleftapex);
        dest(*farright, farrightpt);
        apex(*farright, farrightapex);
        sym(*farleft, checkedge);
        apex(checkedge, checkvertex);
        /* The pointers to the extremal points are restored to the leftmost */
        /*   and rightmost points (rather than topmost and bottommost).     */
        while (checkvertex[0] < farleftpt[0]) {
          lprev(checkedge, *farleft);
          farleftapex = farleftpt;
          farleftpt = checkvertex;
          sym(*farleft, checkedge);
          apex(checkedge, checkvertex);
        }
        while (farrightapex[0] > farrightpt[0]) {
          lprevself(*farright);
          symself(*farright);
          farrightpt = farrightapex;
          apex(*farright, farrightapex);
        }
      }
      return;
    }
    /* Consider eliminating edges from the left triangulation. */
    if (!leftfinished) {
      /* What vertex would be exposed if an edge were deleted? */
      lprev(leftcand, nextedge);
      symself(nextedge);
      apex(nextedge, nextapex);
      /* If nextapex is NULL, then no vertex would be exposed; the */
      /*   triangulation would have been eaten right through.      */
      if (nextapex != (TriPoint) NULL) {
        /* Check whether the edge is Delaunay. */
        badedge = incircle(lowerleft, lowerright, upperleft, nextapex) > 0.0;
        while (badedge) {
          /* Eliminate the edge with an edge flip.  As a result, the    */
          /*   left triangulation will have one more boundary TriTriangle. */
          lnextself(nextedge);
          sym(nextedge, topcasing);
          lnextself(nextedge);
          sym(nextedge, sidecasing);
          bond(nextedge, topcasing);
          bond(leftcand, sidecasing);
          lnextself(leftcand);
          sym(leftcand, outercasing);
          lprevself(nextedge);
          bond(nextedge, outercasing);
          /* Correct the vertices to reflect the edge flip. */
          setorg(leftcand, lowerleft);
          setdest(leftcand, NULL);
          setapex(leftcand, nextapex);
          setorg(nextedge, NULL);
          setdest(nextedge, upperleft);
          setapex(nextedge, nextapex);
          /* Consider the newly exposed vertex. */
          upperleft = nextapex;
          /* What vertex would be exposed if another edge were deleted? */
          TriOrientedTrianglecopy(sidecasing, nextedge);
          apex(nextedge, nextapex);
          if (nextapex != (TriPoint) NULL) {
            /* Check whether the edge is Delaunay. */
            badedge = incircle(lowerleft, lowerright, upperleft, nextapex)
                      > 0.0;
          } else {
            /* Avoid eating right through the triangulation. */
            badedge = 0;
          }
        }
      }
    }
    /* Consider eliminating edges from the right triangulation. */
    if (!rightfinished) {
      /* What vertex would be exposed if an edge were deleted? */
      lnext(rightcand, nextedge);
      symself(nextedge);
      apex(nextedge, nextapex);
      /* If nextapex is NULL, then no vertex would be exposed; the */
      /*   triangulation would have been eaten right through.      */
      if (nextapex != (TriPoint) NULL) {
        /* Check whether the edge is Delaunay. */
        badedge = incircle(lowerleft, lowerright, upperright, nextapex) > 0.0;
        while (badedge) {
          /* Eliminate the edge with an edge flip.  As a result, the     */
          /*   right triangulation will have one more boundary TriTriangle. */
          lprevself(nextedge);
          sym(nextedge, topcasing);
          lprevself(nextedge);
          sym(nextedge, sidecasing);
          bond(nextedge, topcasing);
          bond(rightcand, sidecasing);
          lprevself(rightcand);
          sym(rightcand, outercasing);
          lnextself(nextedge);
          bond(nextedge, outercasing);
          /* Correct the vertices to reflect the edge flip. */
          setorg(rightcand, NULL);
          setdest(rightcand, lowerright);
          setapex(rightcand, nextapex);
          setorg(nextedge, upperright);
          setdest(nextedge, NULL);
          setapex(nextedge, nextapex);
          /* Consider the newly exposed vertex. */
          upperright = nextapex;
          /* What vertex would be exposed if another edge were deleted? */
          TriOrientedTrianglecopy(sidecasing, nextedge);
          apex(nextedge, nextapex);
          if (nextapex != (TriPoint) NULL) {
            /* Check whether the edge is Delaunay. */
            badedge = incircle(lowerleft, lowerright, upperright, nextapex)
                      > 0.0;
          } else {
            /* Avoid eating right through the triangulation. */
            badedge = 0;
          }
        }
      }
    }
    if (leftfinished || (!rightfinished &&
           (incircle(upperleft, lowerleft, lowerright, upperright) > 0.0))) {
      /* Knit the triangulations, adding an edge from `lowerleft' */
      /*   to `upperright'.                                       */
      bond(baseedge, rightcand);
      lprev(rightcand, baseedge);
      setdest(baseedge, lowerleft);
      lowerright = upperright;
      sym(baseedge, rightcand);
      apex(rightcand, upperright);
    } else {
      /* Knit the triangulations, adding an edge from `upperleft' */
      /*   to `lowerright'.                                       */
      bond(baseedge, leftcand);
      lnext(leftcand, baseedge);
      setorg(baseedge, lowerright);
      lowerleft = upperleft;
      sym(baseedge, leftcand);
      apex(leftcand, upperleft);
    }
    if (verbose > 2) {
      printf("  Connecting ");
      printtriangle(&baseedge);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  divconqrecurse()   Recursively form a Delaunay triangulation by the      */
/*                     divide-and-conquer method.                            */
/*                                                                           */
/*  Recursively breaks down the problem into smaller pieces, which are       */
/*  knitted together by mergehulls().  The base cases (problems of two or    */
/*  three points) are handled specially here.                                */
/*                                                                           */
/*  On completion, `farleft' and `farright' are bounding triangles such that */
/*  the origin of `farleft' is the leftmost vertex (breaking ties by         */
/*  choosing the highest leftmost vertex), and the destination of            */
/*  `farright' is the rightmost vertex (breaking ties by choosing the        */
/*  lowest rightmost vertex).                                                */
/*                                                                           */
/*****************************************************************************/

void divconqrecurse(TriPoint *sortarray, int vertices, int axis, TriOrientedTriangle *farleft, TriOrientedTriangle *farright)
{
  TriOrientedTriangle midtri, tri1, tri2, tri3;
  TriOrientedTriangle innerleft, innerright;
  REAL area;
  int divider;

  if (verbose > 2) {
    printf("  Triangulating %d points.\n", vertices);
  }
  if (vertices == 2) {
    /* The triangulation of two vertices is an edge.  An edge is */
    /*   represented by two bounding triangles.                  */
    maketriangle(farleft);
    setorg(*farleft, sortarray[0]);
    setdest(*farleft, sortarray[1]);
    /* The apex is intentionally left NULL. */
    maketriangle(farright);
    setorg(*farright, sortarray[1]);
    setdest(*farright, sortarray[0]);
    /* The apex is intentionally left NULL. */
    bond(*farleft, *farright);
    lprevself(*farleft);
    lnextself(*farright);
    bond(*farleft, *farright);
    lprevself(*farleft);
    lnextself(*farright);
    bond(*farleft, *farright);
    if (verbose > 2) {
      printf("  Creating ");
      printtriangle(farleft);
      printf("  Creating ");
      printtriangle(farright);
    }
    /* Ensure that the origin of `farleft' is sortarray[0]. */
    lprev(*farright, *farleft);
    return;
  } else if (vertices == 3) {
    /* The triangulation of three vertices is either a TriTriangle (with */
    /*   three bounding triangles) or two edges (with four bounding   */
    /*   triangles).  In either case, four triangles are created.     */
    maketriangle(&midtri);
    maketriangle(&tri1);
    maketriangle(&tri2);
    maketriangle(&tri3);
    area = counterclockwise(sortarray[0], sortarray[1], sortarray[2]);
    if (area == 0.0) {
      /* Three collinear points; the triangulation is two edges. */
      setorg(midtri, sortarray[0]);
      setdest(midtri, sortarray[1]);
      setorg(tri1, sortarray[1]);
      setdest(tri1, sortarray[0]);
      setorg(tri2, sortarray[2]);
      setdest(tri2, sortarray[1]);
      setorg(tri3, sortarray[1]);
      setdest(tri3, sortarray[2]);
      /* All apices are intentionally left NULL. */
      bond(midtri, tri1);
      bond(tri2, tri3);
      lnextself(midtri);
      lprevself(tri1);
      lnextself(tri2);
      lprevself(tri3);
      bond(midtri, tri3);
      bond(tri1, tri2);
      lnextself(midtri);
      lprevself(tri1);
      lnextself(tri2);
      lprevself(tri3);
      bond(midtri, tri1);
      bond(tri2, tri3);
      /* Ensure that the origin of `farleft' is sortarray[0]. */
      TriOrientedTrianglecopy(tri1, *farleft);
      /* Ensure that the destination of `farright' is sortarray[2]. */
      TriOrientedTrianglecopy(tri2, *farright);
    } else {
      /* The three points are not collinear; the triangulation is one */
      /*   TriTriangle, namely `midtri'.                                 */
      setorg(midtri, sortarray[0]);
      setdest(tri1, sortarray[0]);
      setorg(tri3, sortarray[0]);
      /* Apices of tri1, tri2, and tri3 are left NULL. */
      if (area > 0.0) {
        /* The vertices are in counterclockwise order. */
        setdest(midtri, sortarray[1]);
        setorg(tri1, sortarray[1]);
        setdest(tri2, sortarray[1]);
        setapex(midtri, sortarray[2]);
        setorg(tri2, sortarray[2]);
        setdest(tri3, sortarray[2]);
      } else {
        /* The vertices are in clockwise order. */
        setdest(midtri, sortarray[2]);
        setorg(tri1, sortarray[2]);
        setdest(tri2, sortarray[2]);
        setapex(midtri, sortarray[1]);
        setorg(tri2, sortarray[1]);
        setdest(tri3, sortarray[1]);
      }
      /* The topology does not depend on how the vertices are ordered. */
      bond(midtri, tri1);
      lnextself(midtri);
      bond(midtri, tri2);
      lnextself(midtri);
      bond(midtri, tri3);
      lprevself(tri1);
      lnextself(tri2);
      bond(tri1, tri2);
      lprevself(tri1);
      lprevself(tri3);
      bond(tri1, tri3);
      lnextself(tri2);
      lprevself(tri3);
      bond(tri2, tri3);
      /* Ensure that the origin of `farleft' is sortarray[0]. */
      TriOrientedTrianglecopy(tri1, *farleft);
      /* Ensure that the destination of `farright' is sortarray[2]. */
      if (area > 0.0) {
        TriOrientedTrianglecopy(tri2, *farright);
      } else {
        lnext(*farleft, *farright);
      }
    }
    if (verbose > 2) {
      printf("  Creating ");
      printtriangle(&midtri);
      printf("  Creating ");
      printtriangle(&tri1);
      printf("  Creating ");
      printtriangle(&tri2);
      printf("  Creating ");
      printtriangle(&tri3);
    }
    return;
  } else {
    /* Split the vertices in half. */
    divider = vertices >> 1;
    /* Recursively triangulate each half. */
    divconqrecurse(sortarray, divider, 1 - axis, farleft, &innerleft);
    divconqrecurse(&sortarray[divider], vertices - divider, 1 - axis,
                   &innerright, farright);
    if (verbose > 1) {
      printf("  Joining triangulations with %d and %d vertices.\n", divider,
             vertices - divider);
    }
    /* Merge the two triangulations into one. */
    mergehulls(farleft, &innerleft, &innerright, farright, axis);
  }
}

long removeghosts(TriOrientedTriangle *startghost)
{
  TriOrientedTriangle searchedge;
  TriOrientedTriangle dissolveedge;
  TriOrientedTriangle deadtri;
  TriPoint markorg;
  long hullsize;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (verbose) {
    printf("  Removing ghost triangles.\n");
  }
  /* Find an edge on the convex hull to start TriPoint location from. */
  lprev(*startghost, searchedge);
  symself(searchedge);
  dummytri[0] = encode(searchedge);
  /* Remove the bounding box and count the convex hull edges. */
  TriOrientedTrianglecopy(*startghost, dissolveedge);
  hullsize = 0;
  do {
    hullsize++;
    lnext(dissolveedge, deadtri);
    lprevself(dissolveedge);
    symself(dissolveedge);
    /* If no PSLG is involved, set the boundary markers of all the points */
    /*   on the convex hull.  If a PSLG is used, this step is done later. */
    if (!poly) {
      /* Watch out for the case where all the input points are collinear. */
      if (dissolveedge.tri != dummytri) {
        org(dissolveedge, markorg);
        if (pointmark(markorg) == 0) {
          setpointmark(markorg, 1);
        }
      }
    }
    /* Remove a bounding TriTriangle from a convex hull TriTriangle. */
    dissolve(dissolveedge);
    /* Find the next bounding TriTriangle. */
    sym(deadtri, dissolveedge);
    /* Delete the bounding TriTriangle. */
    triangledealloc(deadtri.tri);
  } while (!TriOrientedTriangleequal(dissolveedge, *startghost));
  return hullsize;
}

/*****************************************************************************/
/*                                                                           */
/*  divconqdelaunay()   Form a Delaunay triangulation by the divide-and-     */
/*                      conquer method.                                      */
/*                                                                           */
/*  Sorts the points, calls a recursive procedure to triangulate them, and   */
/*  removes the bounding box, setting boundary markers as appropriate.       */
/*                                                                           */
/*****************************************************************************/

long divconqdelaunay()
{
  TriPoint *sortarray;
  TriOrientedTriangle hullleft, hullright;
  int divider;
  int i, j;

  /* Allocate an array of pointers to points for sorting. */
  sortarray = (TriPoint *) malloc(inpoints * sizeof(TriPoint));
  if (sortarray == (TriPoint *) NULL) {
    printf("Error:  Out of memory.\n");
    exit(1);
  }
  traversalinit(&points);
  for (i = 0; i < inpoints; i++) {
    sortarray[i] = pointtraverse();
  }
  if (verbose) {
    printf("  Sorting points.\n");
  }
  /* Sort the points. */
  pointsort(sortarray, inpoints);
  /* Discard duplicate points, which can really mess up the algorithm. */
  i = 0;
  for (j = 1; j < inpoints; j++) {
    if ((sortarray[i][0] == sortarray[j][0])
        && (sortarray[i][1] == sortarray[j][1])) {
      if (!quiet) {
        printf(
"Warning:  A duplicate TriPoint at (%.12g, %.12g) appeared and was ignored.\n",
               sortarray[j][0], sortarray[j][1]);
      }
/*  Commented out - would eliminate TriPoint from output .node file, but causes
    a failure if some segment has this TriPoint as an endpoint.
      setpointmark(sortarray[j], DEADPOINT);
*/
    } else {
      i++;
      sortarray[i] = sortarray[j];
    }
  }
  i++;
  if (dwyer) {
    /* Re-sort the array of points to accommodate alternating cuts. */
    divider = i >> 1;
    if (i - divider >= 2) {
      if (divider >= 2) {
        alternateaxes(sortarray, divider, 1);
      }
      alternateaxes(&sortarray[divider], i - divider, 1);
    }
  }
  if (verbose) {
    printf("  Forming triangulation.\n");
  }
  /* Form the Delaunay triangulation. */
  divconqrecurse(sortarray, i, 0, &hullleft, &hullright);
  free(sortarray);

  return removeghosts(&hullleft);
}

/**                                                                         **/
/**                                                                         **/
/********* Divide-and-conquer Delaunay triangulation ends here       *********/

/*****************************************************************************/
/*                                                                           */
/*  boundingbox()   Form an "infinite" bounding TriTriangle to insert points    */
/*                  into.                                                    */
/*                                                                           */
/*  The points at "infinity" are assigned finite coordinates, which are used */
/*  by the TriPoint location routines, but (mostly) ignored by the Delaunay     */
/*  edge flip routines.                                                      */
/*                                                                           */
/*****************************************************************************/

void boundingbox()
{
  TriOrientedTriangle inftri;          /* Handle for the triangular bounding box. */
  REAL width;

  if (verbose) {
    printf("  Creating triangular bounding box.\n");
  }
  /* Find the width (or height, whichever is larger) of the triangulation. */
  width = xmax - xmin;
  if (ymax - ymin > width) {
    width = ymax - ymin;
  }
  if (width == 0.0) {
    width = 1.0;
  }
  /* Create the vertices of the bounding box. */
  infpoint1 = (TriPoint) malloc(points.itembytes);
  infpoint2 = (TriPoint) malloc(points.itembytes);
  infpoint3 = (TriPoint) malloc(points.itembytes);
  if ((infpoint1 == (TriPoint) NULL) || (infpoint2 == (TriPoint) NULL)
      || (infpoint3 == (TriPoint) NULL)) {
    printf("Error:  Out of memory.\n");
    exit(1);
  }
  infpoint1[0] = xmin - 50.0 * width;
  infpoint1[1] = ymin - 40.0 * width;
  infpoint2[0] = xmax + 50.0 * width;
  infpoint2[1] = ymin - 40.0 * width;
  infpoint3[0] = 0.5 * (xmin + xmax);
  infpoint3[1] = ymax + 60.0 * width;

  /* Create the bounding box. */
  maketriangle(&inftri);
  setorg(inftri, infpoint1);
  setdest(inftri, infpoint2);
  setapex(inftri, infpoint3);
  /* Link dummytri to the bounding box so we can always find an */
  /*   edge to begin searching (TriPoint location) from.           */
  dummytri[0] = (TriTriangle) inftri.tri;
  if (verbose > 2) {
    printf("  Creating ");
    printtriangle(&inftri);
  }
}


/*****************************************************************************/
/*                                                                           */
/*  removebox()   Remove the "infinite" bounding TriTriangle, setting boundary  */
/*                markers as appropriate.                                    */
/*                                                                           */
/*  The triangular bounding box has three boundary triangles (one for each   */
/*  side of the bounding box), and a bunch of triangles fanning out from     */
/*  the three bounding box vertices (one TriTriangle for each edge of the       */
/*  convex hull of the inner mesh).  This routine removes these triangles.   */
/*                                                                           */
/*****************************************************************************/

long removebox()
{
  TriOrientedTriangle deadtri;
  TriOrientedTriangle searchedge;
  TriOrientedTriangle checkedge;
  TriOrientedTriangle nextedge, finaledge, dissolveedge;
  TriPoint markorg;
  long hullsize;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (verbose) {
    printf("  Removing triangular bounding box.\n");
  }
  /* Find a boundary TriTriangle. */
  nextedge.tri = dummytri;
  nextedge.orient = 0;
  symself(nextedge);
  /* Mark a place to stop. */
  lprev(nextedge, finaledge);
  lnextself(nextedge);
  symself(nextedge);
  /* Find a TriTriangle (on the boundary of the TriPoint set) that isn't */
  /*   a bounding box TriTriangle.                                    */
  lprev(nextedge, searchedge);
  symself(searchedge);
  /* Check whether nextedge is another boundary TriTriangle */
  /*   adjacent to the first one.                        */
  lnext(nextedge, checkedge);
  symself(checkedge);
  if (checkedge.tri == dummytri) {
    /* Go on to the next TriTriangle.  There are only three boundary   */
    /*   triangles, and this next TriTriangle cannot be the third one, */
    /*   so it's safe to stop here.                                 */
    lprevself(searchedge);
    symself(searchedge);
  }
  /* Find a new boundary edge to search from, as the current search */
  /*   edge lies on a bounding box TriTriangle and will be deleted.    */
  dummytri[0] = encode(searchedge);
  hullsize = -2l;
  while (!TriOrientedTriangleequal(nextedge, finaledge)) {
    hullsize++;
    lprev(nextedge, dissolveedge);
    symself(dissolveedge);
    /* If not using a PSLG, the vertices should be marked now. */
    /*   (If using a PSLG, markhull() will do the job.)        */
    if (!poly) {
      /* Be careful!  One must check for the case where all the input   */
      /*   points are collinear, and thus all the triangles are part of */
      /*   the bounding box.  Otherwise, the setpointmark() call below  */
      /*   will cause a bad pointer reference.                          */
      if (dissolveedge.tri != dummytri) {
        org(dissolveedge, markorg);
        if (pointmark(markorg) == 0) {
          setpointmark(markorg, 1);
        }
      }
    }
    /* Disconnect the bounding box TriTriangle from the mesh TriTriangle. */
    dissolve(dissolveedge);
    lnext(nextedge, deadtri);
    sym(deadtri, nextedge);
    /* Get rid of the bounding box TriTriangle. */
    triangledealloc(deadtri.tri);
    /* Do we need to turn the corner? */
    if (nextedge.tri == dummytri) {
      /* Turn the corner. */
      TriOrientedTrianglecopy(dissolveedge, nextedge);
    }
  }
  triangledealloc(finaledge.tri);

  free(infpoint1);                  /* Deallocate the bounding box vertices. */
  free(infpoint2);
  free(infpoint3);

  return hullsize;
}

/********* General mesh construction routines begin here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  delaunay()   Form a Delaunay triangulation.                              */
/*                                                                           */
/*****************************************************************************/

long delaunay()
{
  eextras = 0;
  initializetrisegpools();

  if (!quiet) {
    printf(
      "Constructing Delaunay triangulation by divide-and-conquer method.\n");
  }
  return divconqdelaunay();
}


/********* Segment (shell edge) insertion begins here                *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  finddirection()   Find the first TriTriangle on the path from one TriPoint     */
/*                    to another.                                            */
/*                                                                           */
/*  Finds the TriTriangle that intersects a line segment drawn from the         */
/*  origin of `searchtri' to the TriPoint `endpoint', and returns the result    */
/*  in `searchtri'.  The origin of `searchtri' does not change, even though  */
/*  the TriTriangle returned may differ from the one passed in.  This routine   */
/*  is used to find the direction to move in to get from one TriPoint to        */
/*  another.                                                                 */
/*                                                                           */
/*  The return value notes whether the destination or apex of the found      */
/*  TriTriangle is collinear with the two points in question.                   */
/*                                                                           */
/*****************************************************************************/

enum finddirectionresult finddirection(TriOrientedTriangle *searchtri, TriPoint endpoint)
{
  TriOrientedTriangle checktri;
  TriPoint startpoint;
  TriPoint leftpoint, rightpoint;
  REAL leftccw, rightccw;
  int leftflag, rightflag;
  TriTriangle ptr;           /* Temporary variable used by onext() and oprev(). */

  org(*searchtri, startpoint);
  dest(*searchtri, rightpoint);
  apex(*searchtri, leftpoint);
  /* Is `endpoint' to the left? */
  leftccw = counterclockwise(endpoint, startpoint, leftpoint);
  leftflag = leftccw > 0.0;
  /* Is `endpoint' to the right? */
  rightccw = counterclockwise(startpoint, endpoint, rightpoint);
  rightflag = rightccw > 0.0;
  if (leftflag && rightflag) {
    /* `searchtri' faces directly away from `endpoint'.  We could go */
    /*   left or right.  Ask whether it's a TriTriangle or a boundary   */
    /*   on the left.                                                */
    onext(*searchtri, checktri);
    if (checktri.tri == dummytri) {
      leftflag = 0;
    } else {
      rightflag = 0;
    }
  }
  while (leftflag) {
    /* Turn left until satisfied. */
    onextself(*searchtri);
    if (searchtri->tri == dummytri) {
      printf("Internal error in finddirection():  Unable to find a\n");
      printf("  TriTriangle leading from (%.12g, %.12g) to", startpoint[0],
             startpoint[1]);
      printf("  (%.12g, %.12g).\n", endpoint[0], endpoint[1]);
      internalerror();
    }
    apex(*searchtri, leftpoint);
    rightccw = leftccw;
    leftccw = counterclockwise(endpoint, startpoint, leftpoint);
    leftflag = leftccw > 0.0;
  }
  while (rightflag) {
    /* Turn right until satisfied. */
    oprevself(*searchtri);
    if (searchtri->tri == dummytri) {
      printf("Internal error in finddirection():  Unable to find a\n");
      printf("  TriTriangle leading from (%.12g, %.12g) to", startpoint[0],
             startpoint[1]);
      printf("  (%.12g, %.12g).\n", endpoint[0], endpoint[1]);
      internalerror();
    }
    dest(*searchtri, rightpoint);
    leftccw = rightccw;
    rightccw = counterclockwise(startpoint, endpoint, rightpoint);
    rightflag = rightccw > 0.0;
  }
  if (leftccw == 0.0) {
    return LEFTCOLLINEAR;
  } else if (rightccw == 0.0) {
    return RIGHTCOLLINEAR;
  } else {
    return WITHIN;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  segmentintersection()   Find the intersection of an existing segment     */
/*                          and a segment that is being inserted.  Insert    */
/*                          a TriPoint at the intersection, splitting an        */
/*                          existing shell edge.                             */
/*                                                                           */
/*  The segment being inserted connects the apex of splittri to endpoint2.   */
/*  splitshelle is the shell edge being split, and MUST be opposite          */
/*  splittri.  Hence, the edge being split connects the origin and           */
/*  destination of splittri.                                                 */
/*                                                                           */
/*  On completion, splittri is a handle having the newly inserted            */
/*  intersection TriPoint as its origin, and endpoint1 as its destination.      */
/*                                                                           */
/*****************************************************************************/

void segmentintersection(TriOrientedTriangle *splittri, TriOrientedShell *splitshelle, TriPoint endpoint2)
{
  TriPoint endpoint1;
  TriPoint torg, tdest;
  TriPoint leftpoint, rightpoint;
  TriPoint newpoint;
  enum insertsiteresult success;
  enum finddirectionresult collinear;
  REAL ex, ey;
  REAL tx, ty;
  REAL etx, ety;
  REAL split, denom;
  int i;
  TriTriangle ptr;                       /* Temporary variable used by onext(). */

  /* Find the other three segment endpoints. */
  apex(*splittri, endpoint1);
  org(*splittri, torg);
  dest(*splittri, tdest);
  /* Segment intersection formulae; see the Antonio reference. */
  tx = tdest[0] - torg[0];
  ty = tdest[1] - torg[1];
  ex = endpoint2[0] - endpoint1[0];
  ey = endpoint2[1] - endpoint1[1];
  etx = torg[0] - endpoint2[0];
  ety = torg[1] - endpoint2[1];
  denom = ty * ex - tx * ey;
  if (denom == 0.0) {
    printf("Internal error in segmentintersection():");
    printf("  Attempt to find intersection of parallel segments.\n");
    internalerror();
  }
  split = (ey * etx - ex * ety) / denom;
  /* Create the new TriPoint. */
  newpoint = (TriPoint) poolalloc(&points);
  /* Interpolate its coordinate and attributes. */
  for (i = 0; i < 2 + nextras; i++) {
    newpoint[i] = torg[i] + split * (tdest[i] - torg[i]);
  }
  setpointmark(newpoint, mark(*splitshelle));
  if (verbose > 1) {
    printf(
    "  Splitting edge (%.12g, %.12g) (%.12g, %.12g) at (%.12g, %.12g).\n",
           torg[0], torg[1], tdest[0], tdest[1], newpoint[0], newpoint[1]);
  }
  /* Insert the intersection TriPoint.  This should always succeed. */
  success = insertsite(newpoint, splittri, splitshelle, 0, 0);
  if (success != SUCCESSFULPOINT) {
    printf("Internal error in segmentintersection():\n");
    printf("  Failure to split a segment.\n");
    internalerror();
  }
  if (steinerleft > 0) {
    steinerleft--;
  }
  /* Inserting the TriPoint may have caused edge flips.  We wish to rediscover */
  /*   the edge connecting endpoint1 to the new intersection TriPoint.         */
  collinear = finddirection(splittri, endpoint1);
  dest(*splittri, rightpoint);
  apex(*splittri, leftpoint);
  if ((leftpoint[0] == endpoint1[0]) && (leftpoint[1] == endpoint1[1])) {
    onextself(*splittri);
  } else if ((rightpoint[0] != endpoint1[0]) ||
             (rightpoint[1] != endpoint1[1])) {
    printf("Internal error in segmentintersection():\n");
    printf("  Topological inconsistency after splitting a segment.\n");
    internalerror();
  }
  /* `splittri' should have destination endpoint1. */
}

/*****************************************************************************/
/*                                                                           */
/*  scoutsegment()   Scout the first TriTriangle on the path from one endpoint  */
/*                   to another, and check for completion (reaching the      */
/*                   second endpoint), a collinear TriPoint, and the            */
/*                   intersection of two segments.                           */
/*                                                                           */
/*  Returns one if the entire segment is successfully inserted, and zero if  */
/*  the job must be finished by conformingedge() or constrainededge().       */
/*                                                                           */
/*  If the first TriTriangle on the path has the second endpoint as its         */
/*  destination or apex, a shell edge is inserted and the job is done.       */
/*                                                                           */
/*  If the first TriTriangle on the path has a destination or apex that lies on */
/*  the segment, a shell edge is inserted connecting the first endpoint to   */
/*  the collinear TriPoint, and the search is continued from the collinear      */
/*  TriPoint.                                                                   */
/*                                                                           */
/*  If the first TriTriangle on the path has a shell edge opposite its origin,  */
/*  then there is a segment that intersects the segment being inserted.      */
/*  Their intersection TriPoint is inserted, splitting the shell edge.          */
/*                                                                           */
/*  Otherwise, return zero.                                                  */
/*                                                                           */
/*****************************************************************************/

int scoutsegment(TriOrientedTriangle *searchtri, TriPoint endpoint2, int newmark)
{
  TriOrientedTriangle crosstri;
  TriOrientedShell crossedge;
  TriPoint leftpoint, rightpoint;
  TriPoint endpoint1;
  enum finddirectionresult collinear;
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  collinear = finddirection(searchtri, endpoint2);
  dest(*searchtri, rightpoint);
  apex(*searchtri, leftpoint);
  if (((leftpoint[0] == endpoint2[0]) && (leftpoint[1] == endpoint2[1])) ||
      ((rightpoint[0] == endpoint2[0]) && (rightpoint[1] == endpoint2[1]))) {
    /* The segment is already an edge in the mesh. */
    if ((leftpoint[0] == endpoint2[0]) && (leftpoint[1] == endpoint2[1])) {
      lprevself(*searchtri);
    }
    /* Insert a shell edge, if there isn't already one there. */
    insertshelle(searchtri, newmark);
    return 1;
  } else if (collinear == LEFTCOLLINEAR) {
    /* We've collided with a TriPoint between the segment's endpoints. */
    /* Make the collinear TriPoint be the TriTriangle's origin. */
    lprevself(*searchtri);
    insertshelle(searchtri, newmark);
    /* Insert the remainder of the segment. */
    return scoutsegment(searchtri, endpoint2, newmark);
  } else if (collinear == RIGHTCOLLINEAR) {
    /* We've collided with a TriPoint between the segment's endpoints. */
    insertshelle(searchtri, newmark);
    /* Make the collinear TriPoint be the TriTriangle's origin. */
    lnextself(*searchtri);
    /* Insert the remainder of the segment. */
    return scoutsegment(searchtri, endpoint2, newmark);
  } else {
    lnext(*searchtri, crosstri);
    tspivot(crosstri, crossedge);
    /* Check for a crossing segment. */
    if (crossedge.sh == dummysh) {
      return 0;
    } else {
      org(*searchtri, endpoint1);
      /* Insert a TriPoint at the intersection. */
      segmentintersection(&crosstri, &crossedge, endpoint2);
      TriOrientedTrianglecopy(crosstri, *searchtri);
      insertshelle(searchtri, newmark);
      /* Insert the remainder of the segment. */
      return scoutsegment(searchtri, endpoint2, newmark);
    }
  }
}


/*****************************************************************************/
/*                                                                           */
/*  delaunayfixup()   Enforce the Delaunay condition at an edge, fanning out */
/*                    recursively from an existing TriPoint.  Pay special       */
/*                    attention to stacking inverted triangles.              */
/*                                                                           */
/*  This is a support routine for inserting segments into a constrained      */
/*  Delaunay triangulation.                                                  */
/*                                                                           */
/*  The origin of fixuptri is treated as if it has just been inserted, and   */
/*  the local Delaunay condition needs to be enforced.  It is only enforced  */
/*  in one sector, however, that being the angular range defined by          */
/*  fixuptri.                                                                */
/*                                                                           */
/*  This routine also needs to make decisions regarding the "stacking" of    */
/*  triangles.  (Read the description of constrainededge() below before      */
/*  reading on here, so you understand the algorithm.)  If the position of   */
/*  the new TriPoint (the origin of fixuptri) indicates that the vertex before  */
/*  it on the polygon is a reflex vertex, then "stack" the TriTriangle by       */
/*  doing nothing.  (fixuptri is an inverted TriTriangle, which is how stacked  */
/*  triangles are identified.)                                               */
/*                                                                           */
/*  Otherwise, check whether the vertex before that was a reflex vertex.     */
/*  If so, perform an edge flip, thereby eliminating an inverted TriTriangle    */
/*  (popping it off the stack).  The edge flip may result in the creation    */
/*  of a new inverted TriTriangle, depending on whether or not the new vertex   */
/*  is visible to the vertex three edges behind on the polygon.              */
/*                                                                           */
/*  If neither of the two vertices behind the new vertex are reflex          */
/*  vertices, fixuptri and fartri, the TriTriangle opposite it, are not         */
/*  inverted; hence, ensure that the edge between them is locally Delaunay.  */
/*                                                                           */
/*  `leftside' indicates whether or not fixuptri is to the left of the       */
/*  segment being inserted.  (Imagine that the segment is pointing up from   */
/*  endpoint1 to endpoint2.)                                                 */
/*                                                                           */
/*****************************************************************************/

void delaunayfixup(TriOrientedTriangle *fixuptri, int leftside)
{
  TriOrientedTriangle neartri;
  TriOrientedTriangle fartri;
  TriOrientedShell faredge;
  TriPoint nearpoint, leftpoint, rightpoint, farpoint;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  lnext(*fixuptri, neartri);
  sym(neartri, fartri);
  /* Check if the edge opposite the origin of fixuptri can be flipped. */
  if (fartri.tri == dummytri) {
    return;
  }
  tspivot(neartri, faredge);
  if (faredge.sh != dummysh) {
    return;
  }
  /* Find all the relevant vertices. */
  apex(neartri, nearpoint);
  org(neartri, leftpoint);
  dest(neartri, rightpoint);
  apex(fartri, farpoint);
  /* Check whether the previous polygon vertex is a reflex vertex. */
  if (leftside) {
    if (counterclockwise(nearpoint, leftpoint, farpoint) <= 0.0) {
      /* leftpoint is a reflex vertex too.  Nothing can */
      /*   be done until a convex section is found.     */
      return;
    }
  } else {
    if (counterclockwise(farpoint, rightpoint, nearpoint) <= 0.0) {
      /* rightpoint is a reflex vertex too.  Nothing can */
      /*   be done until a convex section is found.      */
      return;
    }
  }
  if (counterclockwise(rightpoint, leftpoint, farpoint) > 0.0) {
    /* fartri is not an inverted TriTriangle, and farpoint is not a reflex */
    /*   vertex.  As there are no reflex vertices, fixuptri isn't an    */
    /*   inverted TriTriangle, either.  Hence, test the edge between the   */
    /*   triangles to ensure it is locally Delaunay.                    */
    if (incircle(leftpoint, farpoint, rightpoint, nearpoint) <= 0.0) {
      return;
    }
    /* Not locally Delaunay; go on to an edge flip. */
  }        /* else fartri is inverted; remove it from the stack by flipping. */
  flip(&neartri);
  lprevself(*fixuptri);    /* Restore the origin of fixuptri after the flip. */
  /* Recursively process the two triangles that result from the flip. */
  delaunayfixup(fixuptri, leftside);
  delaunayfixup(&fartri, leftside);
}

/*****************************************************************************/
/*                                                                           */
/*  constrainededge()   Force a segment into a constrained Delaunay          */
/*                      triangulation by deleting the triangles it           */
/*                      intersects, and triangulating the polygons that      */
/*                      form on each side of it.                             */
/*                                                                           */
/*  Generates a single edge connecting `endpoint1' to `endpoint2'.  The      */
/*  TriTriangle `starttri' has `endpoint1' as its origin.  `newmark' is the     */
/*  boundary marker of the segment.                                          */
/*                                                                           */
/*  To insert a segment, every TriTriangle whose interior intersects the        */
/*  segment is deleted.  The union of these deleted triangles is a polygon   */
/*  (which is not necessarily monotone, but is close enough), which is       */
/*  divided into two polygons by the new segment.  This routine's task is    */
/*  to generate the Delaunay triangulation of these two polygons.            */
/*                                                                           */
/*  You might think of this routine's behavior as a two-step process.  The   */
/*  first step is to walk from endpoint1 to endpoint2, flipping each edge    */
/*  encountered.  This step creates a fan of edges connected to endpoint1,   */
/*  including the desired edge to endpoint2.  The second step enforces the   */
/*  Delaunay condition on each side of the segment in an incremental manner: */
/*  proceeding along the polygon from endpoint1 to endpoint2 (this is done   */
/*  independently on each side of the segment), each vertex is "enforced"    */
/*  as if it had just been inserted, but affecting only the previous         */
/*  vertices.  The result is the same as if the vertices had been inserted   */
/*  in the order they appear on the polygon, so the result is Delaunay.      */
/*                                                                           */
/*  In truth, constrainededge() interleaves these two steps.  The procedure  */
/*  walks from endpoint1 to endpoint2, and each time an edge is encountered  */
/*  and flipped, the newly exposed vertex (at the far end of the flipped     */
/*  edge) is "enforced" upon the previously flipped edges, usually affecting */
/*  only one side of the polygon (depending upon which side of the segment   */
/*  the vertex falls on).                                                    */
/*                                                                           */
/*  The algorithm is complicated by the need to handle polygons that are not */
/*  convex.  Although the polygon is not necessarily monotone, it can be     */
/*  triangulated in a manner similar to the stack-based algorithms for       */
/*  monotone polygons.  For each reflex vertex (local concavity) of the      */
/*  polygon, there will be an inverted TriTriangle formed by one of the edge    */
/*  flips.  (An inverted TriTriangle is one with negative area - that is, its   */
/*  vertices are arranged in clockwise order - and is best thought of as a   */
/*  wrinkle in the fabric of the mesh.)  Each inverted TriTriangle can be       */
/*  thought of as a reflex vertex pushed on the stack, waiting to be fixed   */
/*  later.                                                                   */
/*                                                                           */
/*  A reflex vertex is popped from the stack when a vertex is inserted that  */
/*  is visible to the reflex vertex.  (However, if the vertex behind the     */
/*  reflex vertex is not visible to the reflex vertex, a new inverted        */
/*  TriTriangle will take its place on the stack.)  These details are handled   */
/*  by the delaunayfixup() routine above.                                    */
/*                                                                           */
/*****************************************************************************/

void constrainededge(TriOrientedTriangle *starttri, TriPoint endpoint2, int newmark)
{
  TriOrientedTriangle fixuptri, fixuptri2;
  TriOrientedShell fixupedge;
  TriPoint endpoint1;
  TriPoint farpoint;
  REAL area;
  int collision;
  int done;
  TriTriangle ptr;             /* Temporary variable used by sym() and oprev(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  org(*starttri, endpoint1);
  lnext(*starttri, fixuptri);
  flip(&fixuptri);
  /* `collision' indicates whether we have found a TriPoint directly */
  /*   between endpoint1 and endpoint2.                           */
  collision = 0;
  done = 0;
  do {
    org(fixuptri, farpoint);
    /* `farpoint' is the extreme TriPoint of the polygon we are "digging" */
    /*   to get from endpoint1 to endpoint2.                           */
    if ((farpoint[0] == endpoint2[0]) && (farpoint[1] == endpoint2[1])) {
      oprev(fixuptri, fixuptri2);
      /* Enforce the Delaunay condition around endpoint2. */
      delaunayfixup(&fixuptri, 0);
      delaunayfixup(&fixuptri2, 1);
      done = 1;
    } else {
      /* Check whether farpoint is to the left or right of the segment */
      /*   being inserted, to decide which edge of fixuptri to dig     */
      /*   through next.                                               */
      area = counterclockwise(endpoint1, endpoint2, farpoint);
      if (area == 0.0) {
        /* We've collided with a TriPoint between endpoint1 and endpoint2. */
        collision = 1;
        oprev(fixuptri, fixuptri2);
        /* Enforce the Delaunay condition around farpoint. */
        delaunayfixup(&fixuptri, 0);
        delaunayfixup(&fixuptri2, 1);
        done = 1;
      } else {
        if (area > 0.0) {         /* farpoint is to the left of the segment. */
          oprev(fixuptri, fixuptri2);
          /* Enforce the Delaunay condition around farpoint, on the */
          /*   left side of the segment only.                       */
          delaunayfixup(&fixuptri2, 1);
          /* Flip the edge that crosses the segment.  After the edge is */
          /*   flipped, one of its endpoints is the fan vertex, and the */
          /*   destination of fixuptri is the fan vertex.               */
          lprevself(fixuptri);
        } else {                 /* farpoint is to the right of the segment. */
          delaunayfixup(&fixuptri, 0);
          /* Flip the edge that crosses the segment.  After the edge is */
          /*   flipped, one of its endpoints is the fan vertex, and the */
          /*   destination of fixuptri is the fan vertex.               */
          oprevself(fixuptri);
        }
        /* Check for two intersecting segments. */
        tspivot(fixuptri, fixupedge);
        if (fixupedge.sh == dummysh) {
          flip(&fixuptri);   /* May create an inverted TriTriangle on the left. */
        } else {
          /* We've collided with a segment between endpoint1 and endpoint2. */
          collision = 1;
          /* Insert a TriPoint at the intersection. */
          segmentintersection(&fixuptri, &fixupedge, endpoint2);
          done = 1;
        }
      }
    }
  } while (!done);
  /* Insert a shell edge to make the segment permanent. */
  insertshelle(&fixuptri, newmark);
  /* If there was a collision with an interceding vertex, install another */
  /*   segment connecting that vertex with endpoint2.                     */
  if (collision) {
    /* Insert the remainder of the segment. */
    if (!scoutsegment(&fixuptri, endpoint2, newmark)) {
      constrainededge(&fixuptri, endpoint2, newmark);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  insertsegment()   Insert a PSLG segment into a triangulation.            */
/*                                                                           */
/*****************************************************************************/

void insertsegment(TriPoint endpoint1, TriPoint endpoint2, int newmark)
{
  TriOrientedTriangle searchtri1, searchtri2;
  TriTriangle encodedtri;
  TriPoint checkpoint;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (verbose > 1) {
    printf("  Connecting (%.12g, %.12g) to (%.12g, %.12g).\n",
           endpoint1[0], endpoint1[1], endpoint2[0], endpoint2[1]);
  }

  /* Find a TriTriangle whose origin is the segment's first endpoint. */
  checkpoint = (TriPoint) NULL;
  encodedtri = point2tri(endpoint1);
  if (encodedtri != (TriTriangle) NULL) {
    decode(encodedtri, searchtri1);
    org(searchtri1, checkpoint);
  }
  if (checkpoint != endpoint1) {
    /* Find a boundary TriTriangle to search from. */
    searchtri1.tri = dummytri;
    searchtri1.orient = 0;
    symself(searchtri1);
    /* Search for the segment's first endpoint by TriPoint location. */
    if (locate(endpoint1, &searchtri1) != ONVERTEX) {
      printf(
        "Internal error in insertsegment():  Unable to locate PSLG TriPoint\n");
      printf("  (%.12g, %.12g) in triangulation.\n",
             endpoint1[0], endpoint1[1]);
      internalerror();
    }
  }
  /* Remember this TriTriangle to improve subsequent TriPoint location. */
  TriOrientedTrianglecopy(searchtri1, recenttri);
  /* Scout the beginnings of a path from the first endpoint */
  /*   toward the second.                                   */
  if (scoutsegment(&searchtri1, endpoint2, newmark)) {
    /* The segment was easily inserted. */
    return;
  }
  /* The first endpoint may have changed if a collision with an intervening */
  /*   vertex on the segment occurred.                                      */
  org(searchtri1, endpoint1);

  /* Find a TriTriangle whose origin is the segment's second endpoint. */
  checkpoint = (TriPoint) NULL;
  encodedtri = point2tri(endpoint2);
  if (encodedtri != (TriTriangle) NULL) {
    decode(encodedtri, searchtri2);
    org(searchtri2, checkpoint);
  }
  if (checkpoint != endpoint2) {
    /* Find a boundary TriTriangle to search from. */
    searchtri2.tri = dummytri;
    searchtri2.orient = 0;
    symself(searchtri2);
    /* Search for the segment's second endpoint by TriPoint location. */
    if (locate(endpoint2, &searchtri2) != ONVERTEX) {
      printf(
        "Internal error in insertsegment():  Unable to locate PSLG TriPoint\n");
      printf("  (%.12g, %.12g) in triangulation.\n",
             endpoint2[0], endpoint2[1]);
      internalerror();
    }
  }
  /* Remember this TriTriangle to improve subsequent TriPoint location. */
  TriOrientedTrianglecopy(searchtri2, recenttri);
  /* Scout the beginnings of a path from the second endpoint */
  /*   toward the first.                                     */
  if (scoutsegment(&searchtri2, endpoint1, newmark)) {
    /* The segment was easily inserted. */
    return;
  }
  /* The second endpoint may have changed if a collision with an intervening */
  /*   vertex on the segment occurred.                                       */
  org(searchtri2, endpoint2);

    /* Insert the segment directly into the triangulation. */
    constrainededge(&searchtri1, endpoint2, newmark);
}

/*****************************************************************************/
/*                                                                           */
/*  markhull()   Cover the convex hull of a triangulation with shell edges.  */
/*                                                                           */
/*****************************************************************************/

void markhull()
{
  TriOrientedTriangle hulltri;
  TriOrientedTriangle nexttri;
  TriOrientedTriangle starttri;
  TriTriangle ptr;             /* Temporary variable used by sym() and oprev(). */

  /* Find a TriTriangle handle on the hull. */
  hulltri.tri = dummytri;
  hulltri.orient = 0;
  symself(hulltri);
  /* Remember where we started so we know when to stop. */
  TriOrientedTrianglecopy(hulltri, starttri);
  /* Go once counterclockwise around the convex hull. */
  do {
    /* Create a shell edge if there isn't already one here. */
    insertshelle(&hulltri, 1);
    /* To find the next hull edge, go clockwise around the next vertex. */
    lnextself(hulltri);
    oprev(hulltri, nexttri);
    while (nexttri.tri != dummytri) {
      TriOrientedTrianglecopy(nexttri, hulltri);
      oprev(hulltri, nexttri);
    }
  } while (!TriOrientedTriangleequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  formskeleton()   Create the shell edges of a triangulation, including    */
/*                   PSLG edges and edges on the convex hull.                */
/*                                                                           */
/*  The PSLG edges are read from a .poly file.  The return value is the      */
/*  number of segments in the file.                                          */
/*                                                                           */
/*****************************************************************************/


int formskeleton(int *segmentlist, int *segmentmarkerlist, int numberofsegments)
{
  char polyfilename[6];
  int index;
  TriPoint endpoint1, endpoint2;
  int segments;
  int segmentmarkers;
  int end1, end2;
  int boundmarker;
  int i;

  if (poly) {
    if (!quiet) {
      printf("Inserting segments into Delaunay triangulation.\n");
    }
    strcpy(polyfilename, "input");
    segments = numberofsegments;
    segmentmarkers = segmentmarkerlist != (int *) NULL;
    index = 0;
    /* If segments are to be inserted, compute a mapping */
    /*   from points to triangles.                       */
    if (segments > 0) {
      if (verbose) {
        printf("  Inserting PSLG segments.\n");
      }
      makepointmap();
    }

    boundmarker = 0;
    /* Read and insert the segments. */
    for (i = 1; i <= segments; i++) {
      end1 = segmentlist[index++];
      end2 = segmentlist[index++];
      if (segmentmarkers) {
        boundmarker = segmentmarkerlist[i - 1];
      }
      if ((end1 < firstnumber) || (end1 >= firstnumber + inpoints)) {
        if (!quiet) {
          printf("Warning:  Invalid first endpoint of segment %d in %s.\n", i,
                 polyfilename);
        }
      } else if ((end2 < firstnumber) || (end2 >= firstnumber + inpoints)) {
        if (!quiet) {
          printf("Warning:  Invalid second endpoint of segment %d in %s.\n", i,
                 polyfilename);
        }
      } else {
        endpoint1 = getpoint(end1);
        endpoint2 = getpoint(end2);
        if ((endpoint1[0] == endpoint2[0]) && (endpoint1[1] == endpoint2[1])) {
          if (!quiet) {
            printf("Warning:  Endpoints of segment %d are coincident in %s.\n",
                   i, polyfilename);
          }
        } else {
          insertsegment(endpoint1, endpoint2, boundmarker);
        }
      }
    }
  } else {
    segments = 0;
  }
  if (convex || !poly) {
    /* Enclose the convex hull with shell edges. */
    if (verbose) {
      printf("  Enclosing convex hull with segments.\n");
    }
    markhull();
  }
  return segments;
}

/**                                                                         **/
/**                                                                         **/
/********* Segment (shell edge) insertion ends here                  *********/

/********* Carving out holes and concavities begins here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  infecthull()   Virally infect all of the triangles of the convex hull    */
/*                 that are not protected by shell edges.  Where there are   */
/*                 shell edges, set boundary markers as appropriate.         */
/*                                                                           */
/*****************************************************************************/

void infecthull()
{
  TriOrientedTriangle hulltri;
  TriOrientedTriangle nexttri;
  TriOrientedTriangle starttri;
  TriOrientedShell hulledge;
  TriTriangle **deadtri;
  TriPoint horg, hdest;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  if (verbose) {
    printf("  Marking concavities (external triangles) for elimination.\n");
  }
  /* Find a TriTriangle handle on the hull. */
  hulltri.tri = dummytri;
  hulltri.orient = 0;
  symself(hulltri);
  /* Remember where we started so we know when to stop. */
  TriOrientedTrianglecopy(hulltri, starttri);
  /* Go once counterclockwise around the convex hull. */
  do {
    /* Ignore triangles that are already infected. */
    if (!infected(hulltri)) {
      /* Is the TriTriangle protected by a shell edge? */
      tspivot(hulltri, hulledge);
      if (hulledge.sh == dummysh) {
        /* The TriTriangle is not protected; infect it. */
        infect(hulltri);
        deadtri = (TriTriangle **) poolalloc(&viri);
        *deadtri = hulltri.tri;
      } else {
        /* The TriTriangle is protected; set boundary markers if appropriate. */
        if (mark(hulledge) == 0) {
          setmark(hulledge, 1);
          org(hulltri, horg);
          dest(hulltri, hdest);
          if (pointmark(horg) == 0) {
            setpointmark(horg, 1);
          }
          if (pointmark(hdest) == 0) {
            setpointmark(hdest, 1);
          }
        }
      }
    }
    /* To find the next hull edge, go clockwise around the next vertex. */
    lnextself(hulltri);
    oprev(hulltri, nexttri);
    while (nexttri.tri != dummytri) {
      TriOrientedTrianglecopy(nexttri, hulltri);
      oprev(hulltri, nexttri);
    }
  } while (!TriOrientedTriangleequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  plague()   Spread the virus from all infected triangles to any neighbors */
/*             not protected by shell edges.  Delete all infected triangles. */
/*                                                                           */
/*  This is the procedure that actually creates holes and concavities.       */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase identifies all   */
/*  the triangles that will die, and marks them as infected.  They are       */
/*  marked to ensure that each TriTriangle is added to the virus pool only      */
/*  once, so the procedure will terminate.                                   */
/*                                                                           */
/*  The second phase actually eliminates the infected triangles.  It also    */
/*  eliminates orphaned points.                                              */
/*                                                                           */
/*****************************************************************************/

void plague()
{
  TriOrientedTriangle testtri;
  TriOrientedTriangle neighbor;
  TriTriangle **virusloop;
  TriTriangle **deadtri;
  TriOrientedShell neighborshelle;
  TriPoint testpoint;
  TriPoint norg, ndest;
  TriPoint deadorg, deaddest, deadapex;
  int killorg;
  TriTriangle ptr;             /* Temporary variable used by sym() and onext(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  if (verbose) {
    printf("  Marking neighbors of marked triangles.\n");
  }
  /* Loop through all the infected triangles, spreading the virus to */
  /*   their neighbors, then to their neighbors' neighbors.          */
  traversalinit(&viri);
  virusloop = (TriTriangle **) traverse(&viri);
  while (virusloop != (TriTriangle **) NULL) {
    testtri.tri = *virusloop;
    /* A TriTriangle is marked as infected by messing with one of its shell */
    /*   edges, setting it to an illegal value.  Hence, we have to       */
    /*   temporarily uninfect this TriTriangle so that we can examine its   */
    /*   adjacent shell edges.                                           */
    uninfect(testtri);
    if (verbose > 2) {
      /* Assign the TriTriangle an orientation for convenience in */
      /*   checking its points.                                */
      testtri.orient = 0;
      org(testtri, deadorg);
      dest(testtri, deaddest);
      apex(testtri, deadapex);
      printf("    Checking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             deadorg[0], deadorg[1], deaddest[0], deaddest[1],
             deadapex[0], deadapex[1]);
    }
    /* Check each of the TriTriangle's three neighbors. */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      /* Find the neighbor. */
      sym(testtri, neighbor);
      /* Check for a shell between the TriTriangle and its neighbor. */
      tspivot(testtri, neighborshelle);
      /* Check if the neighbor is nonexistent or already infected. */
      if ((neighbor.tri == dummytri) || infected(neighbor)) {
        if (neighborshelle.sh != dummysh) {
          /* There is a shell edge separating the TriTriangle from its */
          /*   neighbor, but both triangles are dying, so the shell */
          /*   edge dies too.                                       */
          shelledealloc(neighborshelle.sh);
          if (neighbor.tri != dummytri) {
            /* Make sure the shell edge doesn't get deallocated again */
            /*   later when the infected neighbor is visited.         */
            uninfect(neighbor);
            tsdissolve(neighbor);
            infect(neighbor);
          }
        }
      } else {                   /* The neighbor exists and is not infected. */
        if (neighborshelle.sh == dummysh) {
          /* There is no shell edge protecting the neighbor, so */
          /*   the neighbor becomes infected.                   */
          if (verbose > 2) {
            org(neighbor, deadorg);
            dest(neighbor, deaddest);
            apex(neighbor, deadapex);
            printf(
              "    Marking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
                   deadorg[0], deadorg[1], deaddest[0], deaddest[1],
                   deadapex[0], deadapex[1]);
          }
          infect(neighbor);
          /* Ensure that the neighbor's neighbors will be infected. */
          deadtri = (TriTriangle **) poolalloc(&viri);
          *deadtri = neighbor.tri;
        } else {               /* The neighbor is protected by a shell edge. */
          /* Remove this TriTriangle from the shell edge. */
          stdissolve(neighborshelle);
          /* The shell edge becomes a boundary.  Set markers accordingly. */
          if (mark(neighborshelle) == 0) {
            setmark(neighborshelle, 1);
          }
          org(neighbor, norg);
          dest(neighbor, ndest);
          if (pointmark(norg) == 0) {
            setpointmark(norg, 1);
          }
          if (pointmark(ndest) == 0) {
            setpointmark(ndest, 1);
          }
        }
      }
    }
    /* Remark the TriTriangle as infected, so it doesn't get added to the */
    /*   virus pool again.                                             */
    infect(testtri);
    virusloop = (TriTriangle **) traverse(&viri);
  }

  if (verbose) {
    printf("  Deleting marked triangles.\n");
  }
  traversalinit(&viri);
  virusloop = (TriTriangle **) traverse(&viri);
  while (virusloop != (TriTriangle **) NULL) {
    testtri.tri = *virusloop;

    /* Check each of the three corners of the TriTriangle for elimination. */
    /*   This is done by walking around each TriPoint, checking if it is   */
    /*   still connected to at least one live TriTriangle.                 */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      org(testtri, testpoint);
      /* Check if the TriPoint has already been tested. */
      if (testpoint != (TriPoint) NULL) {
        killorg = 1;
        /* Mark the corner of the TriTriangle as having been tested. */
        setorg(testtri, NULL);
        /* Walk counterclockwise about the TriPoint. */
        onext(testtri, neighbor);
        /* Stop upon reaching a boundary or the starting TriTriangle. */
        while ((neighbor.tri != dummytri)
               && (!TriOrientedTriangleequal(neighbor, testtri))) {
          if (infected(neighbor)) {
            /* Mark the corner of this TriTriangle as having been tested. */
            setorg(neighbor, NULL);
          } else {
            /* A live TriTriangle.  The TriPoint survives. */
            killorg = 0;
          }
          /* Walk counterclockwise about the TriPoint. */
          onextself(neighbor);
        }
        /* If we reached a boundary, we must walk clockwise as well. */
        if (neighbor.tri == dummytri) {
          /* Walk clockwise about the TriPoint. */
          oprev(testtri, neighbor);
          /* Stop upon reaching a boundary. */
          while (neighbor.tri != dummytri) {
            if (infected(neighbor)) {
            /* Mark the corner of this TriTriangle as having been tested. */
              setorg(neighbor, NULL);
            } else {
              /* A live TriTriangle.  The TriPoint survives. */
              killorg = 0;
            }
            /* Walk clockwise about the TriPoint. */
            oprevself(neighbor);
          }
        }
        if (killorg) {
          if (verbose > 1) {
            printf("    Deleting TriPoint (%.12g, %.12g)\n",
                   testpoint[0], testpoint[1]);
          }
          pointdealloc(testpoint);
        }
      }
    }

    /* Record changes in the number of boundary edges, and disconnect */
    /*   dead triangles from their neighbors.                         */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      sym(testtri, neighbor);
      if (neighbor.tri == dummytri) {
        /* There is no neighboring TriTriangle on this edge, so this edge    */
        /*   is a boundary edge.  This TriTriangle is being deleted, so this */
        /*   boundary edge is deleted.                                    */
        hullsize--;
      } else {
        /* Disconnect the TriTriangle from its neighbor. */
        dissolve(neighbor);
        /* There is a neighboring TriTriangle on this edge, so this edge */
        /*   becomes a boundary edge when this TriTriangle is deleted.   */
        hullsize++;
      }
    }
    /* Return the dead TriTriangle to the pool of triangles. */
    triangledealloc(testtri.tri);
    virusloop = (TriTriangle **) traverse(&viri);
  }
  /* Empty the virus pool. */
  poolrestart(&viri);
}

/*****************************************************************************/
/*                                                                           */
/*  regionplague()   Spread regional attributes and/or area constraints      */
/*                   (from a .poly file) throughout the mesh.                */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase spreads an       */
/*  attribute and/or an area constraint through a (segment-bounded) region.  */
/*  The triangles are marked to ensure that each TriTriangle is added to the    */
/*  virus pool only once, so the procedure will terminate.                   */
/*                                                                           */
/*  The second phase uninfects all infected triangles, returning them to     */
/*  normal.                                                                  */
/*                                                                           */
/*****************************************************************************/

void regionplague(REAL attribute, REAL area)
{
  TriOrientedTriangle testtri;
  TriOrientedTriangle neighbor;
  TriTriangle **virusloop;
  TriTriangle **regiontri;
  TriOrientedShell neighborshelle;
  TriPoint regionorg, regiondest, regionapex;
  TriTriangle ptr;             /* Temporary variable used by sym() and onext(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  if (verbose > 1) {
    printf("  Marking neighbors of marked triangles.\n");
  }
  /* Loop through all the infected triangles, spreading the attribute      */
  /*   and/or area constraint to their neighbors, then to their neighbors' */
  /*   neighbors.                                                          */
  traversalinit(&viri);
  virusloop = (TriTriangle **) traverse(&viri);
  while (virusloop != (TriTriangle **) NULL) {
    testtri.tri = *virusloop;
    /* A TriTriangle is marked as infected by messing with one of its shell */
    /*   edges, setting it to an illegal value.  Hence, we have to       */
    /*   temporarily uninfect this TriTriangle so that we can examine its   */
    /*   adjacent shell edges.                                           */
    uninfect(testtri);
    if (regionattrib) {
      /* Set an attribute. */
      setelemattribute(testtri, eextras, attribute);
    }
    if (vararea) {
      /* Set an area constraint. */
      setareabound(testtri, area);
    }
    if (verbose > 2) {
      /* Assign the TriTriangle an orientation for convenience in */
      /*   checking its points.                                */
      testtri.orient = 0;
      org(testtri, regionorg);
      dest(testtri, regiondest);
      apex(testtri, regionapex);
      printf("    Checking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             regionorg[0], regionorg[1], regiondest[0], regiondest[1],
             regionapex[0], regionapex[1]);
    }
    /* Check each of the TriTriangle's three neighbors. */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      /* Find the neighbor. */
      sym(testtri, neighbor);
      /* Check for a shell between the TriTriangle and its neighbor. */
      tspivot(testtri, neighborshelle);
      /* Make sure the neighbor exists, is not already infected, and */
      /*   isn't protected by a shell edge.                          */
      if ((neighbor.tri != dummytri) && !infected(neighbor)
          && (neighborshelle.sh == dummysh)) {
        if (verbose > 2) {
          org(neighbor, regionorg);
          dest(neighbor, regiondest);
          apex(neighbor, regionapex);
          printf("    Marking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
                 regionorg[0], regionorg[1], regiondest[0], regiondest[1],
                 regionapex[0], regionapex[1]);
        }
        /* Infect the neighbor. */
        infect(neighbor);
        /* Ensure that the neighbor's neighbors will be infected. */
        regiontri = (TriTriangle **) poolalloc(&viri);
        *regiontri = neighbor.tri;
      }
    }
    /* Remark the TriTriangle as infected, so it doesn't get added to the */
    /*   virus pool again.                                             */
    infect(testtri);
    virusloop = (TriTriangle **) traverse(&viri);
  }

  /* Uninfect all triangles. */
  if (verbose > 1) {
    printf("  Unmarking marked triangles.\n");
  }
  traversalinit(&viri);
  virusloop = (TriTriangle **) traverse(&viri);
  while (virusloop != (TriTriangle **) NULL) {
    testtri.tri = *virusloop;
    uninfect(testtri);
    virusloop = (TriTriangle **) traverse(&viri);
  }
  /* Empty the virus pool. */
  poolrestart(&viri);
}

/*****************************************************************************/
/*                                                                           */
/*  carveholes()   Find the holes and infect them.  Find the area            */
/*                 constraints and infect them.  Infect the convex hull.     */
/*                 Spread the infection and kill triangles.  Spread the      */
/*                 area constraints.                                         */
/*                                                                           */
/*  This routine mainly calls other routines to carry out all these          */
/*  functions.                                                               */
/*                                                                           */
/*****************************************************************************/

void carveholes(REAL *holelist, int holes, REAL *regionlist, int regions)
{
  TriOrientedTriangle searchtri;
  TriOrientedTriangle triangleloop;
  TriOrientedTriangle *regiontris;
  TriTriangle **holetri;
  TriTriangle **regiontri;
  TriPoint searchorg, searchdest;
  TriLocateResultType intersect;
  int i;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (!(quiet || (noholes && convex))) {
    printf("Removing unwanted triangles.\n");
    if (verbose && (holes > 0)) {
      printf("  Marking holes for elimination.\n");
    }
  }

  if (regions > 0) {
    /* Allocate storage for the triangles in which region points fall. */
    regiontris = (TriOrientedTriangle *) malloc(regions * sizeof(TriOrientedTriangle));
    if (regiontris == (TriOrientedTriangle *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }

  if (((holes > 0) && !noholes) || !convex || (regions > 0)) {
    /* Initialize a pool of viri to be used for holes, concavities, */
    /*   regional attributes, and/or regional area constraints.     */
    poolinit(&viri, sizeof(TriTriangle *), VIRUSPERBLOCK, POINTER, 0);
  }

  if (!convex) {
    /* Mark as infected any unprotected triangles on the boundary. */
    /*   This is one way by which concavities are created.         */
    infecthull();
  }

  if ((holes > 0) && !noholes) {
    /* Infect each TriTriangle in which a hole lies. */
    for (i = 0; i < 2 * holes; i += 2) {
      /* Ignore holes that aren't within the bounds of the mesh. */
      if ((holelist[i] >= xmin) && (holelist[i] <= xmax)
          && (holelist[i + 1] >= ymin) && (holelist[i + 1] <= ymax)) {
        /* Start searching from some TriTriangle on the outer boundary. */
        searchtri.tri = dummytri;
        searchtri.orient = 0;
        symself(searchtri);
        /* Ensure that the hole is to the left of this boundary edge; */
        /*   otherwise, locate() will falsely report that the hole    */
        /*   falls within the starting TriTriangle.                      */
        org(searchtri, searchorg);
        dest(searchtri, searchdest);
        if (counterclockwise(searchorg, searchdest, &holelist[i]) > 0.0) {
          /* Find a TriTriangle that contains the hole. */
          intersect = locate(&holelist[i], &searchtri);
          if ((intersect != OUTSIDE) && (!infected(searchtri))) {
            /* Infect the TriTriangle.  This is done by marking the TriTriangle */
            /*   as infect and including the TriTriangle in the virus pool.  */
            infect(searchtri);
            holetri = (TriTriangle **) poolalloc(&viri);
            *holetri = searchtri.tri;
          }
        }
      }
    }
  }

  /* Now, we have to find all the regions BEFORE we carve the holes, because */
  /*   locate() won't work when the triangulation is no longer convex.       */
  /*   (Incidentally, this is the reason why regional attributes and area    */
  /*   constraints can't be used when refining a preexisting mesh, which     */
  /*   might not be convex; they can only be used with a freshly             */
  /*   triangulated PSLG.)                                                   */
  if (regions > 0) {
    /* Find the starting TriTriangle for each region. */
    for (i = 0; i < regions; i++) {
      regiontris[i].tri = dummytri;
      /* Ignore region points that aren't within the bounds of the mesh. */
      if ((regionlist[4 * i] >= xmin) && (regionlist[4 * i] <= xmax) &&
          (regionlist[4 * i + 1] >= ymin) && (regionlist[4 * i + 1] <= ymax)) {
        /* Start searching from some TriTriangle on the outer boundary. */
        searchtri.tri = dummytri;
        searchtri.orient = 0;
        symself(searchtri);
        /* Ensure that the region TriPoint is to the left of this boundary */
        /*   edge; otherwise, locate() will falsely report that the     */
        /*   region TriPoint falls within the starting TriTriangle.           */
        org(searchtri, searchorg);
        dest(searchtri, searchdest);
        if (counterclockwise(searchorg, searchdest, &regionlist[4 * i]) >
            0.0) {
          /* Find a TriTriangle that contains the region TriPoint. */
          intersect = locate(&regionlist[4 * i], &searchtri);
          if ((intersect != OUTSIDE) && (!infected(searchtri))) {
            /* Record the TriTriangle for processing after the */
            /*   holes have been carved.                    */
            TriOrientedTrianglecopy(searchtri, regiontris[i]);
          }
        }
      }
    }
  }

  if (viri.items > 0) {
    /* Carve the holes and concavities. */
    plague();
  }
  /* The virus pool should be empty now. */

  if (regions > 0) {
    if (!quiet) {
      if (regionattrib) {
        if (vararea) {
          printf("Spreading regional attributes and area constraints.\n");
        } else {
          printf("Spreading regional attributes.\n");
        }
      } else { 
        printf("Spreading regional area constraints.\n");
      }
    }
    if (regionattrib && !refine) {
      /* Assign every TriTriangle a regional attribute of zero. */
      traversalinit(&triangles);
      triangleloop.orient = 0;
      triangleloop.tri = triangletraverse();
      while (triangleloop.tri != (TriTriangle *) NULL) {
        setelemattribute(triangleloop, eextras, 0.0);
        triangleloop.tri = triangletraverse();
      }
    }
    for (i = 0; i < regions; i++) {
      if (regiontris[i].tri != dummytri) {
        /* Make sure the TriTriangle under consideration still exists. */
        /*   It may have been eaten by the virus.                   */
        if (regiontris[i].tri[3] != (TriTriangle) NULL) {
          /* Put one TriTriangle in the virus pool. */
          infect(regiontris[i]);
          regiontri = (TriTriangle **) poolalloc(&viri);
          *regiontri = regiontris[i].tri;
          /* Apply one region's attribute and/or area constraint. */
          regionplague(regionlist[4 * i + 2], regionlist[4 * i + 3]);
          /* The virus pool should be empty now. */
        }
      }
    }
    if (regionattrib && !refine) {
      /* Note the fact that each TriTriangle has an additional attribute. */
      eextras++;
    }
  }

  /* Free up memory. */
  if (((holes > 0) && !noholes) || !convex || (regions > 0)) {
    pooldeinit(&viri);
  }
  if (regions > 0) {
    free(regiontris);
  }
}

/**                                                                         **/
/**                                                                         **/
/********* Carving out holes and concavities ends here               *********/

/********* Mesh quality maintenance begins here                      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  findcircumcenter()   Find the circumcenter of a TriTriangle.                */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  coordinates.  The xi-eta coordinate system is defined in terms of the    */
/*  TriTriangle:  the origin of the TriTriangle is the origin of the coordinate    */
/*  system; the destination of the TriTriangle is one unit along the xi axis;   */
/*  and the apex of the TriTriangle is one unit along the eta axis.             */
/*                                                                           */
/*  The return value indicates which edge of the TriTriangle is shortest.       */
/*                                                                           */
/*****************************************************************************/

enum circumcenterresult findcircumcenter(TriPoint torg, TriPoint tdest, TriPoint tapex, 
										 TriPoint circumcenter, REAL *xi, REAL *eta)
{
  REAL xdo, ydo, xao, yao, xad, yad;
  REAL dodist, aodist, addist;
  REAL denominator;
  REAL dx, dy;

  circumcentercount++;

  /* Compute the circumcenter of the TriTriangle. */
  xdo = tdest[0] - torg[0];
  ydo = tdest[1] - torg[1];
  xao = tapex[0] - torg[0];
  yao = tapex[1] - torg[1];
  dodist = xdo * xdo + ydo * ydo;
  aodist = xao * xao + yao * yao;
  if (noexact) {
    denominator = 0.5 / (xdo * yao - xao * ydo);
  } else {
    /* Use the counterclockwise() routine to ensure a positive (and */
    /*   reasonably accurate) result, avoiding any possibility of   */
    /*   division by zero.                                          */
    denominator = 0.5 / counterclockwise(tdest, tapex, torg);
    /* Don't count the above as an orientation test. */
    counterclockcount--;
  }
  circumcenter[0] = torg[0] - (ydo * aodist - yao * dodist) * denominator;  
  circumcenter[1] = torg[1] + (xdo * aodist - xao * dodist) * denominator;  

  /* To interpolate TriPoint attributes for the new TriPoint inserted at  */
  /*   the circumcenter, define a coordinate system with a xi-axis, */
  /*   directed from the TriTriangle's origin to its destination, and  */
  /*   an eta-axis, directed from its origin to its apex.           */
  /*   Calculate the xi and eta coordinates of the circumcenter.    */
  dx = circumcenter[0] - torg[0];
  dy = circumcenter[1] - torg[1];
  *xi = (dx * yao - xao * dy) * (2.0 * denominator);
  *eta = (xdo * dy - dx * ydo) * (2.0 * denominator);

  xad = tapex[0] - tdest[0];
  yad = tapex[1] - tdest[1];
  addist = xad * xad + yad * yad;
  if ((addist < dodist) && (addist < aodist)) {
    return OPPOSITEORG;
  } else if (dodist < aodist) {
    return OPPOSITEAPEX;
  } else {
    return OPPOSITEDEST;
  }
}


/*****************************************************************************/
/*                                                                           */
/*  highorder()   Create extra nodes for quadratic subparametric elements.   */
/*                                                                           */
/*****************************************************************************/

void highorder()
{
  TriOrientedTriangle triangleloop, trisym;
  TriOrientedShell checkmark;
  TriPoint newpoint;
  TriPoint torg, tdest;
  int i;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  if (!quiet) {
    printf("Adding vertices for second-order triangles.\n");
  }
  /* The following line ensures that dead items in the pool of nodes    */
  /*   cannot be allocated for the extra nodes associated with high     */
  /*   order elements.  This ensures that the primary nodes (at the     */
  /*   corners of elements) will occur earlier in the output files, and */
  /*   have lower indices, than the extra nodes.                        */
  points.deaditemstack = (VOID *) NULL;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each TriTriangle.  If there isn't another TriTriangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent TriTriangle, operate on the edge only if the current TriTriangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (TriTriangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == dummytri)) {
        org(triangleloop, torg);
        dest(triangleloop, tdest);
        /* Create a new node in the middle of the edge.  Interpolate */
        /*   its attributes.                                         */
        newpoint = (TriPoint) poolalloc(&points);
        for (i = 0; i < 2 + nextras; i++) {
          newpoint[i] = 0.5 * (torg[i] + tdest[i]);
        }
        /* Set the new node's marker to zero or one, depending on */
        /*   whether it lies on a boundary.                       */
        setpointmark(newpoint, trisym.tri == dummytri);
        if (useshelles) {
          tspivot(triangleloop, checkmark);
          /* If this edge is a segment, transfer the marker to the new node. */
          if (checkmark.sh != dummysh) {
            setpointmark(newpoint, mark(checkmark));
          }
        }
        if (verbose > 1) {
          printf("  Creating (%.12g, %.12g).\n", newpoint[0], newpoint[1]);
        }
        /* Record the new node in the (one or two) adjacent elements. */
        triangleloop.tri[highorderindex + triangleloop.orient] =
                (TriTriangle) newpoint;
        if (trisym.tri != dummytri) {
          trisym.tri[highorderindex + trisym.orient] = (TriTriangle) newpoint;
        }
      }
    }
    triangleloop.tri = triangletraverse();
  }
}

/********* File I/O routines begin here                              *********/
/**                                                                         **/
/**                                                                         **/


/*****************************************************************************/
/*                                                                           */
/*  transfernodes()   Read the points from memory.                           */
/*                                                                           */
/*****************************************************************************/

void transfernodes(REAL *pointlist, REAL *pointattriblist, int *pointmarkerlist, 
				   int numberofpoints, int numberofpointattribs)
{
  TriPoint pointloop;
  REAL x, y;
  int i, j;
  int coordindex;
  int attribindex;

  inpoints = numberofpoints;
  mesh_dim = 2;
  nextras = numberofpointattribs;
  readnodefile = 0;
  if (inpoints < 3) {
    printf("Error:  Input must have at least three input points.\n");
    exit(1);
  }

  initializepointpool();

  /* Read the points. */
  coordindex = 0;
  attribindex = 0;
  for (i = 0; i < inpoints; i++) {
    pointloop = (TriPoint) poolalloc(&points);
    /* Read the TriPoint coordinates. */
    x = pointloop[0] = pointlist[coordindex++];
    y = pointloop[1] = pointlist[coordindex++];
    /* Read the TriPoint attributes. */
    for (j = 0; j < numberofpointattribs; j++) {
      pointloop[2 + j] = pointattriblist[attribindex++];
    }
    if (pointmarkerlist != (int *) NULL) {
      /* Read a TriPoint marker. */
      setpointmark(pointloop, pointmarkerlist[i]);
    } else {
      /* If no markers are specified, they default to zero. */
      setpointmark(pointloop, 0);
    }
    x = pointloop[0];
    y = pointloop[1];
    /* Determine the smallest and largest x and y coordinates. */
    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
    }
  }

  /* Nonexistent x value used as a flag to mark circle events in sweepline */
  /*   Delaunay algorithm.                                                 */
  xminextreme = 10 * xmin - 9 * xmax;
}



/*****************************************************************************/
/*                                                                           */
/*  writenodes()   Number the points and write them to a .node file.         */
/*                                                                           */
/*  To save memory, the TriPoint numbers are written over the shell markers     */
/*  after the points are written to a file.                                  */
/*                                                                           */
/*****************************************************************************/


void writenodes(REAL **pointlist, REAL **pointattriblist, int **pointmarkerlist)
{
  REAL *plist;
  REAL *palist;
  int *pmlist;
  int coordindex;
  int attribindex;
  TriPoint pointloop;
  int pointnumber;
  int i;

  if (!quiet) {
    printf("Writing points.\n");
  }
  /* Allocate memory for output points if necessary. */
  if (*pointlist == (REAL *) NULL) {
    *pointlist = (REAL *) malloc(points.items * 2 * sizeof(REAL));
    if (*pointlist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for output TriPoint attributes if necessary. */
  if ((nextras > 0) && (*pointattriblist == (REAL *) NULL)) {
    *pointattriblist = (REAL *) malloc(points.items * nextras * sizeof(REAL));
    if (*pointattriblist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for output TriPoint markers if necessary. */
  if (!nobound && (*pointmarkerlist == (int *) NULL)) {
    *pointmarkerlist = (int *) malloc(points.items * sizeof(int));
    if (*pointmarkerlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  plist = *pointlist;
  palist = *pointattriblist;
  pmlist = *pointmarkerlist;
  coordindex = 0;
  attribindex = 0;

  traversalinit(&points);
  pointloop = pointtraverse();
  pointnumber = firstnumber;
  while (pointloop != (TriPoint) NULL) {
    /* X and y coordinates. */
    plist[coordindex++] = pointloop[0];
    plist[coordindex++] = pointloop[1];
    /* Point attributes. */
    for (i = 0; i < nextras; i++) {
      palist[attribindex++] = pointloop[2 + i];
    }
    if (!nobound) {
      /* Copy the boundary marker. */
      pmlist[pointnumber - firstnumber] = pointmark(pointloop);
    }

    setpointmark(pointloop, pointnumber);
    pointloop = pointtraverse();
    pointnumber++;
  }

}

/*****************************************************************************/
/*                                                                           */
/*  numbernodes()   Number the points.                                       */
/*                                                                           */
/*  Each TriPoint is assigned a marker equal to its number.                     */
/*                                                                           */
/*  Used when writenodes() is not called because no .node file is written.   */
/*                                                                           */
/*****************************************************************************/

void numbernodes()
{
  TriPoint pointloop;
  int pointnumber;

  traversalinit(&points);
  pointloop = pointtraverse();
  pointnumber = firstnumber;
  while (pointloop != (TriPoint) NULL) {
    setpointmark(pointloop, pointnumber);
    pointloop = pointtraverse();
    pointnumber++;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  writeelements()   Write the triangles to an .ele file.                   */
/*                                                                           */
/*****************************************************************************/

void writeelements(int **trianglelist, REAL **triangleattriblist)
{
  int *tlist;
  REAL *talist;
  int pointindex;
  int attribindex;
  TriOrientedTriangle triangleloop;
  TriPoint p1, p2, p3;
  TriPoint mid1, mid2, mid3;
  int elementnumber;
  int i;

  if (!quiet) {
    printf("Writing triangles.\n");
  }
  /* Allocate memory for output triangles if necessary. */
  if (*trianglelist == (int *) NULL) {
    *trianglelist = (int *) malloc(triangles.items *
                               ((order + 1) * (order + 2) / 2) * sizeof(int));
    if (*trianglelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for output TriTriangle attributes if necessary. */
  if ((eextras > 0) && (*triangleattriblist == (REAL *) NULL)) {
    *triangleattriblist = (REAL *) malloc(triangles.items * eextras *
                                          sizeof(REAL));
    if (*triangleattriblist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  tlist = *trianglelist;
  talist = *triangleattriblist;
  pointindex = 0;
  attribindex = 0;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  triangleloop.orient = 0;
  elementnumber = firstnumber;
  while (triangleloop.tri != (TriTriangle *) NULL) {
    org(triangleloop, p1);
    dest(triangleloop, p2);
    apex(triangleloop, p3);
    if (order == 1) {
      tlist[pointindex++] = pointmark(p1);
      tlist[pointindex++] = pointmark(p2);
      tlist[pointindex++] = pointmark(p3);
    } else {
      mid1 = (TriPoint) triangleloop.tri[highorderindex + 1];
      mid2 = (TriPoint) triangleloop.tri[highorderindex + 2];
      mid3 = (TriPoint) triangleloop.tri[highorderindex];
      tlist[pointindex++] = pointmark(p1);
      tlist[pointindex++] = pointmark(p2);
      tlist[pointindex++] = pointmark(p3);
      tlist[pointindex++] = pointmark(mid1);
      tlist[pointindex++] = pointmark(mid2);
      tlist[pointindex++] = pointmark(mid3);
    }

    for (i = 0; i < eextras; i++) {
      talist[attribindex++] = elemattribute(triangleloop, i);
    }

    triangleloop.tri = triangletraverse();
    elementnumber++;
  }

}

/*****************************************************************************/
/*                                                                           */
/*  writepoly()   Write the segments and holes to a .poly file.              */
/*                                                                           */
/*****************************************************************************/

void writepoly(int **segmentlist, int **segmentmarkerlist)
{
  int *slist;
  int *smlist;
  int index;
  TriOrientedShell shelleloop;
  TriPoint endpoint1, endpoint2;
  int shellenumber;

  if (!quiet) {
    printf("Writing segments.\n");
  }
  /* Allocate memory for output segments if necessary. */
  if (*segmentlist == (int *) NULL) {
    *segmentlist = (int *) malloc(shelles.items * 2 * sizeof(int));
    if (*segmentlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for output segment markers if necessary. */
  if (!nobound && (*segmentmarkerlist == (int *) NULL)) {
    *segmentmarkerlist = (int *) malloc(shelles.items * sizeof(int));
    if (*segmentmarkerlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  slist = *segmentlist;
  smlist = *segmentmarkerlist;
  index = 0;

  traversalinit(&shelles);
  shelleloop.sh = shelletraverse();
  shelleloop.shorient = 0;
  shellenumber = firstnumber;
  while (shelleloop.sh != (TriShell *) NULL) {
    sorg(shelleloop, endpoint1);
    sdest(shelleloop, endpoint2);
    /* Copy indices of the segment's two endpoints. */
    slist[index++] = pointmark(endpoint1);
    slist[index++] = pointmark(endpoint2);
    if (!nobound) {
      /* Copy the boundary marker. */
      smlist[shellenumber - firstnumber] = mark(shelleloop);
    }

    shelleloop.sh = shelletraverse();
    shellenumber++;
  }

}

/*****************************************************************************/
/*                                                                           */
/*  writeedges()   Write the edges to a .edge file.                          */
/*                                                                           */
/*****************************************************************************/

void writeedges(int **edgelist, int **edgemarkerlist)
{
  int *elist;
  int *emlist;
  int index;
  TriOrientedTriangle triangleloop, trisym;
  TriOrientedShell checkmark;
  TriPoint p1, p2;
  int edgenumber;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */
  TriShell sptr;                      /* Temporary variable used by tspivot(). */

  if (!quiet) {
    printf("Writing edges.\n");
  }
  /* Allocate memory for edges if necessary. */
  if (*edgelist == (int *) NULL) {
    *edgelist = (int *) malloc(edges * 2 * sizeof(int));
    if (*edgelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for edge markers if necessary. */
  if (!nobound && (*edgemarkerlist == (int *) NULL)) {
    *edgemarkerlist = (int *) malloc(edges * sizeof(int));
    if (*edgemarkerlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  elist = *edgelist;
  emlist = *edgemarkerlist;
  index = 0;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  edgenumber = firstnumber;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each TriTriangle.  If there isn't another TriTriangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent TriTriangle, operate on the edge only if the current TriTriangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (TriTriangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == dummytri)) {
        org(triangleloop, p1);
        dest(triangleloop, p2);
        elist[index++] = pointmark(p1);
        elist[index++] = pointmark(p2);
        if (nobound) {
        } else {
          /* Edge number, indices of two endpoints, and a boundary marker. */
          /*   If there's no shell edge, the boundary marker is zero.      */
          if (useshelles) {
            tspivot(triangleloop, checkmark);
            if (checkmark.sh == dummysh) {
              emlist[edgenumber - firstnumber] = 0;
            } else {
              emlist[edgenumber - firstnumber] = mark(checkmark);
            }
          } else {
            emlist[edgenumber - firstnumber] = trisym.tri == dummytri;
          }
        }
        edgenumber++;
      }
    }
    triangleloop.tri = triangletraverse();
  }

}

/*****************************************************************************/
/*                                                                           */
/*  writevoronoi()   Write the Voronoi diagram to a .v.node and .v.edge      */
/*                   file.                                                   */
/*                                                                           */
/*  The Voronoi diagram is the geometric dual of the Delaunay triangulation. */
/*  Hence, the Voronoi vertices are listed by traversing the Delaunay        */
/*  triangles, and the Voronoi edges are listed by traversing the Delaunay   */
/*  edges.                                                                   */
/*                                                                           */
/*  WARNING:  In order to assign numbers to the Voronoi vertices, this       */
/*  procedure messes up the shell edges or the extra nodes of every          */
/*  element.  Hence, you should call this procedure last.                    */
/*                                                                           */
/*****************************************************************************/

void writevoronoi(REAL **vpointlist, REAL **vpointattriblist, int **vpointmarkerlist, int **vedgelist,
                  int **vedgemarkerlist, REAL **vnormlist)
{
  REAL *plist;
  REAL *palist;
  int *elist;
  REAL *normlist;
  int coordindex;
  int attribindex;
  TriOrientedTriangle triangleloop, trisym;
  TriPoint torg, tdest, tapex;
  REAL circumcenter[2];
  REAL xi, eta;
  int vnodenumber, vedgenumber;
  int p1, p2;
  int i;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (!quiet) {
    printf("Writing Voronoi vertices.\n");
  }
  /* Allocate memory for Voronoi vertices if necessary. */
  if (*vpointlist == (REAL *) NULL) {
    *vpointlist = (REAL *) malloc(triangles.items * 2 * sizeof(REAL));
    if (*vpointlist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  /* Allocate memory for Voronoi vertex attributes if necessary. */
  if (*vpointattriblist == (REAL *) NULL) {
    *vpointattriblist = (REAL *) malloc(triangles.items * nextras *
                                        sizeof(REAL));
    if (*vpointattriblist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  *vpointmarkerlist = (int *) NULL;
  plist = *vpointlist;
  palist = *vpointattriblist;
  coordindex = 0;
  attribindex = 0;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  triangleloop.orient = 0;
  vnodenumber = firstnumber;
  while (triangleloop.tri != (TriTriangle *) NULL) {
    org(triangleloop, torg);
    dest(triangleloop, tdest);
    apex(triangleloop, tapex);
    findcircumcenter(torg, tdest, tapex, circumcenter, &xi, &eta);
    /* X and y coordinates. */
    plist[coordindex++] = circumcenter[0];
    plist[coordindex++] = circumcenter[1];
    for (i = 2; i < 2 + nextras; i++) {
      /* Interpolate the TriPoint attributes at the circumcenter. */
      palist[attribindex++] = torg[i] + xi * (tdest[i] - torg[i])
                                     + eta * (tapex[i] - torg[i]);
    }

    * (int *) (triangleloop.tri + 6) = vnodenumber;
    triangleloop.tri = triangletraverse();
    vnodenumber++;
  }


  if (!quiet) {
    printf("Writing Voronoi edges.\n");
  }
  /* Allocate memory for output Voronoi edges if necessary. */
  if (*vedgelist == (int *) NULL) {
    *vedgelist = (int *) malloc(edges * 2 * sizeof(int));
    if (*vedgelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  *vedgemarkerlist = (int *) NULL;
  /* Allocate memory for output Voronoi norms if necessary. */
  if (*vnormlist == (REAL *) NULL) {
    *vnormlist = (REAL *) malloc(edges * 2 * sizeof(REAL));
    if (*vnormlist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  elist = *vedgelist;
  normlist = *vnormlist;
  coordindex = 0;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  vedgenumber = firstnumber;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each TriTriangle.  If there isn't another TriTriangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent TriTriangle, operate on the edge only if the current TriTriangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (TriTriangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == dummytri)) {
        /* Find the number of this TriTriangle (and Voronoi vertex). */
        p1 = * (int *) (triangleloop.tri + 6);
        if (trisym.tri == dummytri) {
          org(triangleloop, torg);
          dest(triangleloop, tdest);
          /* Copy an infinite ray.  Index of one endpoint, and -1. */
          elist[coordindex] = p1;
          normlist[coordindex++] = tdest[1] - torg[1];
          elist[coordindex] = -1;
          normlist[coordindex++] = torg[0] - tdest[0];
        } else {
          /* Find the number of the adjacent TriTriangle (and Voronoi vertex). */
          p2 = * (int *) (trisym.tri + 6);
          /* Finite edge.  Write indices of two endpoints. */
          elist[coordindex] = p1;
          normlist[coordindex++] = 0.0;
          elist[coordindex] = p2;
          normlist[coordindex++] = 0.0;
        }
        vedgenumber++;
      }
    }
    triangleloop.tri = triangletraverse();
  }

}


void writeneighbors(int **neighborlist)
{
  int *nlist;
  int index;
  TriOrientedTriangle triangleloop, trisym;
  int elementnumber;
  int neighbor1, neighbor2, neighbor3;
  TriTriangle ptr;                         /* Temporary variable used by sym(). */

  if (!quiet) {
    printf("Writing neighbors.\n");
  }
  /* Allocate memory for neighbors if necessary. */
  if (*neighborlist == (int *) NULL) {
    *neighborlist = (int *) malloc(triangles.items * 3 * sizeof(int));
    if (*neighborlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  nlist = *neighborlist;
  index = 0;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  triangleloop.orient = 0;
  elementnumber = firstnumber;
  while (triangleloop.tri != (TriTriangle *) NULL) {
    * (int *) (triangleloop.tri + 6) = elementnumber;
    triangleloop.tri = triangletraverse();
    elementnumber++;
  }
  * (int *) (dummytri + 6) = -1;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  elementnumber = firstnumber;
  while (triangleloop.tri != (TriTriangle *) NULL) {
    triangleloop.orient = 1;
    sym(triangleloop, trisym);
    neighbor1 = * (int *) (trisym.tri + 6);
    triangleloop.orient = 2;
    sym(triangleloop, trisym);
    neighbor2 = * (int *) (trisym.tri + 6);
    triangleloop.orient = 0;
    sym(triangleloop, trisym);
    neighbor3 = * (int *) (trisym.tri + 6);
    nlist[index++] = neighbor1;
    nlist[index++] = neighbor2;
    nlist[index++] = neighbor3;

    triangleloop.tri = triangletraverse();
    elementnumber++;
  }

}

/**                                                                         **/
/**                                                                         **/
/********* File I/O routines end here                                *********/

/*****************************************************************************/
/*                                                                           */
/*  quality_statistics()   Print statistics about the quality of the mesh.   */
/*                                                                           */
/*****************************************************************************/

void quality_statistics()
{
  TriOrientedTriangle triangleloop;
  TriPoint p[3];
  REAL cossquaretable[8];
  REAL ratiotable[16];
  REAL dx[3], dy[3];
  REAL edgelength[3];
  REAL dotproduct;
  REAL cossquare;
  REAL triarea;
  REAL shortest, longest;
  REAL trilongest2;
  REAL smallestarea, biggestarea;
  REAL triminaltitude2;
  REAL minaltitude;
  REAL triaspect2;
  REAL worstaspect;
  REAL smallestangle, biggestangle;
  REAL radconst, degconst;
  int angletable[18];
  int aspecttable[16];
  int aspectindex;
  int tendegree;
  int acutebiggest;
  int i, ii, j, k;

  printf("Mesh quality statistics:\n\n");
  radconst = PI / 18.0;
  degconst = 180.0 / PI;
  for (i = 0; i < 8; i++) {
    cossquaretable[i] = cos(radconst * (REAL) (i + 1));
    cossquaretable[i] = cossquaretable[i] * cossquaretable[i];
  }
  for (i = 0; i < 18; i++) {
    angletable[i] = 0;
  }

  ratiotable[0]  =      1.5;      ratiotable[1]  =     2.0;
  ratiotable[2]  =      2.5;      ratiotable[3]  =     3.0;
  ratiotable[4]  =      4.0;      ratiotable[5]  =     6.0;
  ratiotable[6]  =     10.0;      ratiotable[7]  =    15.0;
  ratiotable[8]  =     25.0;      ratiotable[9]  =    50.0;
  ratiotable[10] =    100.0;      ratiotable[11] =   300.0;
  ratiotable[12] =   1000.0;      ratiotable[13] = 10000.0;
  ratiotable[14] = 100000.0;      ratiotable[15] =     0.0;
  for (i = 0; i < 16; i++) {
    aspecttable[i] = 0;
  }

  worstaspect = 0.0;
  minaltitude = xmax - xmin + ymax - ymin;
  minaltitude = minaltitude * minaltitude;
  shortest = minaltitude;
  longest = 0.0;
  smallestarea = minaltitude;
  biggestarea = 0.0;
  worstaspect = 0.0;
  smallestangle = 0.0;
  biggestangle = 2.0;
  acutebiggest = 1;

  traversalinit(&triangles);
  triangleloop.tri = triangletraverse();
  triangleloop.orient = 0;
  while (triangleloop.tri != (TriTriangle *) NULL) {
    org(triangleloop, p[0]);
    dest(triangleloop, p[1]);
    apex(triangleloop, p[2]);
    trilongest2 = 0.0;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dx[i] = p[j][0] - p[k][0];
      dy[i] = p[j][1] - p[k][1];
      edgelength[i] = dx[i] * dx[i] + dy[i] * dy[i];
      if (edgelength[i] > trilongest2) {
        trilongest2 = edgelength[i];
      }
      if (edgelength[i] > longest) {
        longest = edgelength[i];
      }
      if (edgelength[i] < shortest) {
        shortest = edgelength[i];
      }
    }

    triarea = counterclockwise(p[0], p[1], p[2]);
    if (triarea < smallestarea) {
      smallestarea = triarea;
    }
    if (triarea > biggestarea) {
      biggestarea = triarea;
    }
    triminaltitude2 = triarea * triarea / trilongest2;
    if (triminaltitude2 < minaltitude) {
      minaltitude = triminaltitude2;
    }
    triaspect2 = trilongest2 / triminaltitude2;
    if (triaspect2 > worstaspect) {
      worstaspect = triaspect2;
    }
    aspectindex = 0;
    while ((triaspect2 > ratiotable[aspectindex] * ratiotable[aspectindex])
           && (aspectindex < 15)) {
      aspectindex++;
    }
    aspecttable[aspectindex]++;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dotproduct = dx[j] * dx[k] + dy[j] * dy[k];
      cossquare = dotproduct * dotproduct / (edgelength[j] * edgelength[k]);
      tendegree = 8;
      for (ii = 7; ii >= 0; ii--) {
        if (cossquare > cossquaretable[ii]) {
          tendegree = ii;
        }
      }
      if (dotproduct <= 0.0) {
        angletable[tendegree]++;
        if (cossquare > smallestangle) {
          smallestangle = cossquare;
        }
        if (acutebiggest && (cossquare < biggestangle)) {
          biggestangle = cossquare;
        }
      } else {
        angletable[17 - tendegree]++;
        if (acutebiggest || (cossquare > biggestangle)) {
          biggestangle = cossquare;
          acutebiggest = 0;
        }
      }
    }
    triangleloop.tri = triangletraverse();
  }

  shortest = sqrt(shortest);
  longest = sqrt(longest);
  minaltitude = sqrt(minaltitude);
  worstaspect = sqrt(worstaspect);
  smallestarea *= 2.0;
  biggestarea *= 2.0;
  if (smallestangle >= 1.0) {
    smallestangle = 0.0;
  } else {
    smallestangle = degconst * acos(sqrt(smallestangle));
  }
  if (biggestangle >= 1.0) {
    biggestangle = 180.0;
  } else {
    if (acutebiggest) {
      biggestangle = degconst * acos(sqrt(biggestangle));
    } else {
      biggestangle = 180.0 - degconst * acos(sqrt(biggestangle));
    }
  }

  printf("  Smallest area: %16.5g   |  Largest area: %16.5g\n",
         smallestarea, biggestarea);
  printf("  Shortest edge: %16.5g   |  Longest edge: %16.5g\n",
         shortest, longest);
  printf("  Shortest altitude: %12.5g   |  Largest aspect ratio: %8.5g\n\n",
         minaltitude, worstaspect);
  printf("  Aspect ratio histogram:\n");
  printf("  1.1547 - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
         ratiotable[0], aspecttable[0], ratiotable[7], ratiotable[8],
         aspecttable[8]);
  for (i = 1; i < 7; i++) {
    printf("  %6.6g - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
           ratiotable[i - 1], ratiotable[i], aspecttable[i],
           ratiotable[i + 7], ratiotable[i + 8], aspecttable[i + 8]);
  }
  printf("  %6.6g - %-6.6g    :  %8d    | %6.6g -            :  %8d\n",
         ratiotable[6], ratiotable[7], aspecttable[7], ratiotable[14],
         aspecttable[15]);
  printf(
"  (Triangle aspect ratio is longest edge divided by shortest altitude)\n\n");
  printf("  Smallest angle: %15.5g   |  Largest angle: %15.5g\n\n",
         smallestangle, biggestangle);
  printf("  Angle histogram:\n");
  for (i = 0; i < 9; i++) {
    printf("    %3d - %3d degrees:  %8d    |    %3d - %3d degrees:  %8d\n",
           i * 10, i * 10 + 10, angletable[i],
           i * 10 + 90, i * 10 + 100, angletable[i + 9]);
  }
  printf("\n");
}

/*****************************************************************************/
/*                                                                           */
/*  statistics()   Print all sorts of cool facts.                            */
/*                                                                           */
/*****************************************************************************/

void statistics()
{
  printf("\nStatistics:\n\n");
  printf("  Input points: %d\n", inpoints);
  if (refine) {
    printf("  Input triangles: %d\n", inelements);
  }
  if (poly) {
    printf("  Input segments: %d\n", insegments);
    if (!refine) {
      printf("  Input holes: %d\n", holes);
    }
  }

  printf("\n  Mesh points: %ld\n", points.items);
  printf("  Mesh triangles: %ld\n", triangles.items);
  printf("  Mesh edges: %ld\n", edges);
  if (poly || refine) {
    printf("  Mesh boundary edges: %ld\n", hullsize);
    printf("  Mesh segments: %ld\n\n", shelles.items);
  } else {
    printf("  Mesh convex hull edges: %ld\n\n", hullsize);
  }
  if (verbose) {
    quality_statistics();
    printf("Memory allocation statistics:\n\n");
    printf("  Maximum number of points: %ld\n", points.maxitems);
    printf("  Maximum number of triangles: %ld\n", triangles.maxitems);
    if (shelles.maxitems > 0) {
      printf("  Maximum number of segments: %ld\n", shelles.maxitems);
    }
    if (viri.maxitems > 0) {
      printf("  Maximum number of viri: %ld\n", viri.maxitems);
    }
    if (badsegments.maxitems > 0) {
      printf("  Maximum number of encroached segments: %ld\n",
             badsegments.maxitems);
    }
    if (badtriangles.maxitems > 0) {
      printf("  Maximum number of bad triangles: %ld\n",
             badtriangles.maxitems);
    }
    if (splaynodes.maxitems > 0) {
      printf("  Maximum number of splay tree nodes: %ld\n",
             splaynodes.maxitems);
    }
    printf("  Approximate heap memory use (bytes): %ld\n\n",
           points.maxitems * points.itembytes
           + triangles.maxitems * triangles.itembytes
           + shelles.maxitems * shelles.itembytes
           + viri.maxitems * viri.itembytes
           + badsegments.maxitems * badsegments.itembytes
           + badtriangles.maxitems * badtriangles.itembytes
           + splaynodes.maxitems * splaynodes.itembytes);

    printf("Algorithmic statistics:\n\n");
    printf("  Number of incircle tests: %ld\n", incirclecount);
    printf("  Number of orientation tests: %ld\n", counterclockcount);
    if (hyperbolacount > 0) {
      printf("  Number of right-of-hyperbola tests: %ld\n",
             hyperbolacount);
    }
    if (circumcentercount > 0) {
      printf("  Number of circumcenter computations: %ld\n",
             circumcentercount);
    }
    if (circletopcount > 0) {
      printf("  Number of circle top computations: %ld\n",
             circletopcount);
    }
    printf("\n");
  }
}

/*****************************************************************************/
/*                                                                           */
/*  main() or triangulate()   Gosh, do everything.                           */
/*                                                                           */
/*  The sequence is roughly as follows.  Many of these steps can be skipped, */
/*  depending on the command line switches.                                  */
/*                                                                           */
/*  - Initialize constants and parse the command line.                       */
/*  - Read the points from a file and either                                 */
/*    - triangulate them (no -r), or                                         */
/*    - read an old mesh from files and reconstruct it (-r).                 */
/*  - Insert the PSLG segments (-p), and possibly segments on the convex     */
/*      hull (-c).                                                           */
/*  - Read the holes (-p), regional attributes (-pA), and regional area      */
/*      constraints (-pa).  Carve the holes and concavities, and spread the  */
/*      regional attributes and area constraints.                            */
/*  - Enforce the constraints on minimum angle (-q) and maximum area (-a).   */
/*      Also enforce the conforming Delaunay property (-q and -a).           */
/*  - Compute the number of edges in the resulting mesh.                     */
/*  - Promote the mesh's linear triangles to higher order elements (-o).     */
/*  - Write the output files and print the statistics.                       */
/*  - Check the consistency and Delaunay property of the mesh (-C).          */
/*                                                                           */
/*****************************************************************************/

void triangulate_srt(char *triswitches, SrtTriangulationIO *in)
{

	  triangleinit();

	  parsecommandline(1, &triswitches);

	  transfernodes(in->pointlist, in->pointattributelist, in->pointmarkerlist,
					in->numberofpoints, in->numberofpointattributes);



	  hullsize = delaunay();                          /* Triangulate the points. */


}



TriTriangle *getdummytri()
{
	return dummytri;
}
