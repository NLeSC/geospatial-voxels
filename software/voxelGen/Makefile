all: optimized

optimized:
	cc -O3 -fopenmp -c voxelize.c readObj.c
	g++ -O3 -c splitMesh.cpp voxelGen.cpp
	cc -O3 voxelize.o readObj.o splitMesh.o voxelGen.o -fopenmp -lm -lstdc++ -o voxelGen
	rm voxelize.o readObj.o splitMesh.o voxelGen.o

debug:
	cc -g -fopenmp -c voxelize.c readObj.c
	g++ -g -c splitMesh.cpp voxelGen.cpp
	cc -g voxelize.o readObj.o splitMesh.o voxelGen.o -fopenmp -lm -lstdc++ -o voxelGen
	rm voxelize.o readObj.o splitMesh.o voxelGen.o

clean:
	rm voxelGen
