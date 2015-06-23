#include "voxelGen.h"

/*
 * here we can possibly ask for parameters from the user.
 * The two parameters are the voxel size(x,y,z) and the connectivity level desired (CO)
 */
int main(int argc, char * argv[]) {
    Vertex *vertices;
    Face *mesh;
    int nVertices, nTriangles;
#ifdef PERF_LOG
    struct timeval stop, start;
    unsigned long long t;
#endif
    Color mc;
    mc.r = 255;
    mc.g = 255;
    mc.b = 255;
    int mid;
    Vertex vSize;
    vSize.x = 1.0;
    vSize.y = 1.0;
    vSize.z = 1.0;
    int co = 26;

    if (argc != 5) {
        fprintf(stderr, "Number of arguments is incorrect. ./read_obj <with_split_mesh: 0 | 1> <num_threads> <input_file_absolute_path> <output_file_absolute_path>\n");
        return -1;
    }
    if (readObj(argv[3], &vertices, &nVertices, &mesh, &nTriangles))
        goto out;

    printf("The file contained %d vertexes and %d faces!!!\n", nVertices, nTriangles);

#ifdef PERF_LOG
    gettimeofday(&start, NULL);
#endif

    /* File stream for output (to do) */
    FILE *outFile;

    if (!(outFile = fopen(argv[4], "w"))) {
        fprintf(stderr, "ERROR: the path for file to store the voxels is invalid %s\n", argv[4]);
        return -1;
    }

    /*voxelize without splitting the Mesh*/
    if (atoi(argv[1]) == 0) {
        Face **mesh_ptr = (Face **) malloc(sizeof(Face*)*nTriangles);
        for (int i = 0; i < nTriangles; i++)
            mesh_ptr[i] = &mesh[i];
        if (voxelizeMesh(atoi(argv[2]), outFile, vertices, nVertices, mesh_ptr, nTriangles, &mc, mid, &vSize, co)) {
            free(mesh_ptr);
            goto out;
        }
        free(mesh_ptr);
    }
    /*voxelize by splitting the Mesh*/
    else {
        if (splitMesh(atoi(argv[2]), outFile, vertices, nVertices, mesh, nTriangles, &mc, mid, &vSize, co))
            goto out;
    }


#ifdef PERF_LOG
    gettimeofday(&stop, NULL);
    t = 1000 * (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000;
    printf("took %llu ms\n", t);
#endif

out:
    fclose(outFile);
    if (vertices)
        free(vertices);
    if (mesh)
        free(mesh);

    return 0;
}
