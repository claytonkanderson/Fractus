//
//  TetraGroup.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "Object.hpp"
#include "Tetrahedron.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtx/component_wise.hpp>
#include <chrono>

using namespace std::chrono;

////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    float x_width, y_width, z_width, height, elastic, poisson, toughness;
	glm::ivec3 div;
	glm::vec3 center_pos, angular_vel, com_vel;
    float theta, phi;
	glm::mat4 model;
} TetraGroupInits;

////////////////////////////////////////////////////////////////////////////////
//
// This class owns a vector of tetrahedron and the vertices that they reference.
// All simulated objects in this codebase will be represented by TetraGroups.
// 
class TetraGroup: public Object
{
public:
    
    // CREATORS
	TetraGroup();

	~TetraGroup();

	void Initialize();

	void CreateGrid(glm::vec3 bottomLeftPos, glm::vec3 topRightPos, glm::ivec3 resolution);
    void AddCube(glm::vec3 bottomLeftPos, glm::vec3 topRightPos, int mask);
    
    // MANIPULATORS
    void AddTetrahedron(Tetrahedron * tet) {Tetrahedra.push_back(tet);}
    void UpdateVertices(float deltaT, int numSteps);
    void ZeroForces();
    void ApplyGravity();
    void GroundResponse();
    void ComputeDeformationForces();
    void ComputeFracture();
    void ComputeSeparation();
    void Init(const TetraGroupInits &inits);
    void SetInitialConditions(const glm::mat4 &model, const glm::vec3 &comVel, const glm::vec3 &angularVel);
    void SetConstants(float elasticMod, float poisson, float toughness);
    void Reset() {Init(inits);}
    void SetDrawingColors();
    void UpdateMasses();
    void UpdateBetaMatrices();
    void ConnectivityAddendum();
    bool CheckFractureFaces(const std::vector<Tetrahedron *> &tets, Vertex * vert);
    Vertex * DuplicateVertex(Vertex * oldVert);
    void VertexSplit(Vertex * oldVert, Vertex * newVert, glm::vec3 &planeNormal, glm::vec3 &planePoint);
    int TetrahedraIdentifyConnectivity(Tetrahedron * tet1, Tetrahedron * tet2);
    Vertex * TetrahedraMutualVertex(Tetrahedron *tet0, Tetrahedron *tet1);
    std::pair<Vertex *,Vertex *> TetrahedraMutualVertices(Tetrahedron *tet0, Tetrahedron *tet1);
    Triangle * TetrahedraMutualEdge(Tetrahedron *tet0, Tetrahedron *tet1);
    void FractureTetrahedra(Tetrahedron * tet, Vertex * oldVertex, Vertex * newVertex, const glm::vec3 &planeNormal, bool planeSnap);
    void SplitConnectedComponents(const std::vector<Tetrahedron *> &tets);
    void SplitConnectedComponents(Vertex * suspectVert);
    void SplitTetrahedraByVertex(Tetrahedron * tet, Vertex * oldVertex, Vertex * newVertex);
    void SplitTetrahedraByVertex(Tetrahedron * loneTet, Tetrahedron * otherTet);
    void SplitTetrahedraByEdge(Tetrahedron * tet, Vertex * oldV0, Vertex * oldV1, Vertex * newV0, Vertex * newV1);
    void SplitTetrahedraByEdge(Tetrahedron * loneTet, Tetrahedron * otherTet);
    void SplitTetrahedraByPlane(Tetrahedron * tet, const glm::vec3 &planePoint, const glm::vec3 &planeNormal);
    bool IsIsolated(Tetrahedron *);
    void IsolateTetrahedra(Tetrahedron * tet);
    void ArbitrarySplit(Vertex * fractureVert, glm::vec3 fractureNormal);
    std::vector<Vertex *> GetVertexNeighbors(Vertex * vert);
    glm::vec3 FindBestFracturePlane(Vertex * fracturingVertex, glm::vec3 eigenDir);
    bool ValidFractureVert(Vertex * fracturingVertex);
    bool PlaneIntersectTetrahedra(Tetrahedron * tet, const glm::vec3 &planePoint, const glm::vec3 &planeNormal);
    
    // ACCESSORS
    Tetrahedron * GetTetrahedron(int index) {return Tetrahedra[index];}
    Vertex * GetVertex(int index) {return &Vertices[index];}
    Vertex * GetVertex(glm::vec3 pos, float epsilon);
	unsigned int GetVertexIndex(Vertex * vertex) const;
    std::vector<duration<double>> GetTimers() {return timers;}
	const std::vector<Tetrahedron *> & GetTetrahedra() const { return Tetrahedra; }
    
	size_t GetNumVerts() const { return Vertices.size(); }
	static constexpr size_t GetMaxNumVerts() { return 8000; }

    bool meshView = false;
    bool snapToPlane = true;
    int addendumFrames = 0;
    float Toughness = 4000.0f;
    std::vector<Vertex> Vertices;
    
    // Initial conditions
    TetraGroupInits inits;
    
    // GPU & Drawing
    GLuint shader;
    GLuint VAO, vertexBuffer, normalBuffer;
    GLuint elementbuffer;
    void SetupGPU();
    void UpdateGPU();
    void SetShader(GLuint shader) { this->shader = shader; }
    void Draw(glm::mat4 &PV);
    
    Material mat;
    
private:
    std::vector<Tetrahedron *> Tetrahedra;
    std::vector<duration<double>> timers;
    duration<double> timeSpan;
};
