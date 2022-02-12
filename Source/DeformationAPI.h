#pragma once

extern void* mData;

////////////////////////////////////////////////////////////////////////////////

extern "C" int __declspec(dllexport) __stdcall Initialize(
	void* vertexPositions,
	int numVertices,
	int maxNumVertices,
	void* tetrahedraIndices,
	int numTetrahedra,
	int maxNumTetrahedra,
	float lambda,
	float psi,
	float mu,
	float phi,
	float toughness,
	float density);

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Destroy();

////////////////////////////////////////////////////////////////////////////////
//
extern "C" int __declspec(dllexport) __stdcall Deform(
	void* vertexPositions,
	void* vertexVelocities,
	bool useVelocities,
	void* vertexForces,
	bool useForces,
	void* tetrahedraIndices,
	int* numVertices,
	int* numTetrahedra,
	float timestep,
	int numSteps);

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Delaunay3D(
	float* inoutVertices,
	int* inoutNumVertices,
	int maxNumVertices,
	int* outTetrahedraIndices,
	int maxNumTetrahedra,
	int* outNumTetrahedra,
	bool* success
);

////////////////////////////////////////////////////////////////////////////////
//
// cubeMinMax should be {minX, minY, minZ, maxX, maxY, maxZ}
// planeNormalOrigin should be {pX, pY, pZ, nX, nY, nZ}
// outVertices should be preallocated to the maximum number of vertices, which I think is 10
// outTetrahedraIndices should be preallocated to 4 * maximum number of tetrahedra, which I think is 3?
//
extern "C" void __declspec(dllexport) __stdcall TetrahedralizeCubeIntersection(
	float* cubeMinMax, // x, y, z, x, y, z
	float* planeNormalOrigin,
	float* outVertices,
	int* outNumVertices,
	int* outTetrahedraIndices,
	int* outNumTetrahedra,
	bool positivePlane
);

////////////////////////////////////////////////////////////////////////////////