#ifdef __cplusplus
#include <iostream>
#include <iterator>
#include <list>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <assert.h>

extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>

typedef struct {
    double x, y, z;
} Vertex;

typedef struct {
    Vertex *p1, *p2, *p3;
} Face;

/*
 * the data type that is temporarily(internally) used for working with
 * the centre point of voxels as Point3I (i,j,k) in reference to a
 * bounding box; such voxels will be redefined(embedded in R3) as
 * \Point3d (3 double) as {x,y,z}
 */
typedef struct {
    unsigned int x, y, z;
} Voxel;

typedef struct {
    Vertex start, end;
} LineSegment;

typedef struct {
    char r, g, b;
} Color;

int readObj(char *fileName, Vertex **_Vertices, int *_numV, Face **_Faces, int *_numF);

int voxelizeMesh(int num_threads, FILE *outFile, Vertex *vertices, unsigned int nVertices, Face **mesh, unsigned int nTriangles, Color *mc, int mid, Vertex *vSize, int co);

#ifdef __cplusplus
}

int splitMesh(int num_threads, FILE *outFile, Vertex *vertices, unsigned int nVertices, Face *mesh, unsigned int nTriangles, Color *mc, int mid, Vertex *vSize, int co);

#endif
