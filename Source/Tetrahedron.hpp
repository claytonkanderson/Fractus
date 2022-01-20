//
//  Tetrahedron.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "Vertex.hpp"
#include "Triangle.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/quaternion.hpp>
#include <GL/glew.h>
#include <array>

class Tetrahedron
{
public:
    Tetrahedron() {}
    Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3);
  
    void UpdateTriangles();
    
    void SetConstants(float elastic, float poisson)
    {   this->elasticit_modulus = elastic; this->poisson_ratio = poisson;
        this->lambda = elastic * poisson / ((1 + poisson)*(1-2*poisson));
        this->mu = elastic / (2*(1 + poisson));
    }
    void SetConstants(float lambda, float mu, float phi, float psi)
    { this->lambda = lambda; this->mu = mu; this->phi = phi; this->psi = psi; }
    
    void ComputeDeformationForces();
    void ComputeFractureForces();
    Vertex * GetVertex(int i) {return Vertices[i];}
    
    void ComputeMMat(glm::vec3 eigenVector, glm::mat3 &outputMat);
    int NetConnectivity();
    
    float ComputeSignedDistance(glm::vec3 normal, glm::vec3 planePoint);
	bool VertexInTriangle(int triIndex, Vertex * v);
	bool VertexInTetrahedron(Vertex * v);
    
    bool TriangleInTetrahedron(Triangle * tri);
    void ReplaceVertex(Vertex * oldVertex, Vertex * newVertex);
    void UpdateMasses();
    void UpdateBeta();
    
    bool drawFlag = true;
    
	float GetVolume() const { return volume; }

	static constexpr int NumTriangles = 4;
	static constexpr int NumVertices = 4;

	std::array<Vertex *, NumVertices> Vertices;
	std::array<Triangle, NumTriangles> Triangles;

private:
	void ComputeVolume();

    float volume;
    float mass;
    float density = 1000.0f;
    float phi, psi, lambda, mu;
    float elasticit_modulus;
    float poisson_ratio;
    
    glm::mat4 beta; // Inverse of material coordinate transformation
	glm::mat4x3 Pmat; // World coordinate positions
	glm::mat4x3 Vmat; // World coordinate velocities
	glm::mat3 totalStress;
};
