//
//  Vertex.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "glm/glm.hpp"

class Tetrahedron;

////////////////////////////////////////////////////////////////////////////////

class Vertex 
{
public:

	Vertex();

	void Set(const glm::vec3 &p, const glm::vec3 &v, const glm::vec3 &f);
	void ApplyForce(glm::vec3 &force);
	void ZeroForce();

    void Update(float deltaT);
    
    const glm::vec3 & getPos() const { return Position; }
	const glm::vec3 & getVel() const { return Velocity; }
    const glm::vec3 & getForce() const { return Force; }
	const glm::vec3 & GetMaterialCoordinate() const { return MaterialCoordinate; }
    float getMass() const { return mass; }
    
    void setMass(float mass) {this->mass = mass; invMass = 1.0f / mass;}
    void setPos(glm::vec3 pos) {Position = pos;}
    void setVel(glm::vec3 vel) {Velocity = vel;}
	void SetMaterialCoordinate(const glm::vec3 & materialCoordinate) { MaterialCoordinate = materialCoordinate; }

	void SetCompressiveForce(glm::vec3 force);
	void SetTensileForce(glm::vec3 force);
	bool AddConnectedTetrahedron(Tetrahedron * tet);
    
	void RemoveConnection(Tetrahedron * tet);
	bool AddConnection(Tetrahedron * tet);
    
    glm::vec3 DrawForce;
    glm::vec3 Force;
    glm::vec3 DeformationForce;
    std::vector<glm::vec3> CompressiveForces;
    std::vector<glm::vec3> TensileForces;
    std::vector<Tetrahedron *> tetrahedra;
    
    int numForces;
    unsigned int identifier;
    int numConnections = 0;
    
    int maxConnectivity = 40;
    
private:
    glm::vec3 Position;
    glm::vec3 Normal;
    glm::vec3 MaterialCoordinate;
    glm::vec3 Velocity;
    
    float mass;
    float invMass;
};