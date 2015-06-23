#include "voxelGen.h"

/*
 * Vertices are unique
 * Can we use Vertex pointer to check uniqness?
 */
typedef std::map<Vertex*, std::list<Face*> > Index;

int splitMesh(int num_threads, FILE *outFile, Vertex *vertices, unsigned int nVertices, Face *mesh, unsigned int nTriangles, Color *mc, int mid, Vertex *vSize, int co) {

    /*
     * Make the index
     * every vertex is linked to the list of faces where it is used
     */
    Index index;
    for (unsigned int currentFace = 0; currentFace < nTriangles; ++currentFace) {

        /* these should insert a vertex only when there isn't one with the same coordinates in the set */
        index[mesh[currentFace].p1].push_back(&mesh[currentFace]);
        index[mesh[currentFace].p2].push_back(&mesh[currentFace]);
        index[mesh[currentFace].p3].push_back(&mesh[currentFace]);
    }

    /*
     * Generate the connected components
     * faces must be unique!!!
     */
    std::set<Face *> facesUsed;
    for (unsigned int currentFace = 0; currentFace < nTriangles; ++currentFace) {

        /* Already processed, so continue */
        if (facesUsed.count(&mesh[currentFace])) continue;

        /* the faces that we will pass to voxelizeMesh */
        std::vector<Face *> facesInConnectedComponent;
        std::stack<Face *> facesToCheck;

        /* Initialise stack with an unused face */
        facesToCheck.push(&mesh[currentFace]);
        facesUsed.insert(&mesh[currentFace]);
        facesInConnectedComponent.push_back(&mesh[currentFace]);
        while (!facesToCheck.empty()) {

            /* Check face on top of the stack */
            if (!facesToCheck.top())
                assert(0);
            facesInConnectedComponent.push_back(facesToCheck.top());
            Vertex* vertexList[] = {(*facesToCheck.top()).p1, (*facesToCheck.top()).p2, (*facesToCheck.top()).p3};
            for (int i = 0; i < 3; i++) {
                for (std::list<Face *>::iterator currentFaceInConnectedComponent = index[vertexList[i]].begin();
                        currentFaceInConnectedComponent != index[vertexList[i]].end();
                        ++currentFaceInConnectedComponent) {

                    /* Add the incident faces that are not used yet to the stack of faces to check */
                    if (facesUsed.count(*currentFaceInConnectedComponent) == 0) {
                        facesToCheck.push(*currentFaceInConnectedComponent);
                        facesUsed.insert(*currentFaceInConnectedComponent);
                        if (!(*currentFaceInConnectedComponent))
                            assert(0);
                        facesInConnectedComponent.push_back(*currentFaceInConnectedComponent);
                    }
                }

            }

            facesToCheck.pop();
        }

        /* Voxelise this face */
        //if (voxelizeMesh(num_threads, file_name, vertices, nVertices, facesInConnectedComponent[0], facesInConnectedComponent.size(), mc, mid, vSize, co))
        Face** res = &facesInConnectedComponent[0];
        if (voxelizeMesh(num_threads, outFile, vertices, nVertices, res, facesInConnectedComponent.size(), mc, mid, vSize, co))
            return -1;
    }

    return 0;
}
