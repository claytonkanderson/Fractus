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

class Vertex {
public:
    Vertex() {
        mass = 1.0f;
        invMass = 1.0f;
        CompressiveForces.resize(maxConnectivity);
        TensileForces.resize(maxConnectivity);
        tetrahedra.resize(maxConnectivity);
        
        glm::vec3 zero = glm::vec3(0);
        
        for (int i = 0; i < maxConnectivity; i++)
        {
            CompressiveForces[i] = zero;
            TensileForces[i] = zero;
            tetrahedra[i] = nullptr;
        }
        
        numForces = 0;
        
    }
    void Set(const glm::vec3 &p, const glm::vec3 &v, const glm::vec3 &f)
    {Position=p; Velocity = v; Force = f;}
    void ApplyForce(glm::vec3 &force) {Force += force;}
    void ZeroForce() {
        Force = glm::vec3(0);
        numForces = 0;
        
        for (int i = 0; i < maxConnectivity; i++)
        {
            CompressiveForces[i] = glm::vec3(0);
            TensileForces[i] = glm::vec3(0);
        }
    }
    void Update(float deltaT);
    
    glm::vec3 getPos() {return Position;}
    glm::vec3 getVel() {return Velocity;}
    glm::vec3 getForce() {return Force;}
    float getMass() {return mass;}
    
    void setMass(float mass) {this->mass = mass; invMass = 1.0f / mass;}
    void setPos(glm::vec3 pos) {Position = pos;}
    void setVel(glm::vec3 vel) {Velocity = vel;}
    
    void SetCompressiveForce(glm::vec3 force)
    {
        if (numForces < maxConnectivity)
        {
            CompressiveForces[numForces] = force;
        }
    }
    void SetTensileForce(glm::vec3 force)
    {
        if (numForces < maxConnectivity)
        {
            TensileForces[numForces] = force;
            numConnections++;
        }
    }
    bool setTetrahedron(Tetrahedron * tet) {
        if (numConnections < maxConnectivity)
        {
            tetrahedra[numConnections] = tet;
            numConnections++;
            return true;
        }
        return false;
    }
    
    void RemoveConnection(Tetrahedron * tet)
    {
        for (int i = 0; i < tetrahedra.size(); i++)
        {
            if (tetrahedra[i] == tet)
            {
                tetrahedra.erase(tetrahedra.begin() + i);
                numConnections--;
            }
        }
    }
    
    bool AddConnection(Tetrahedron * tet)
    {
        if (numConnections < maxConnectivity)
        {
            tetrahedra[numConnections] = tet;
            numConnections++;
            return true;
        }
        return false;
    }
    
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
    glm::vec3 TexCoord;
    glm::vec3 Velocity;
    
    float mass;
    float invMass;
};