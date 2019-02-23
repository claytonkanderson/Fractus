//
//  Triangle.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "Vertex.hpp"
#include <glm/gtx/norm.hpp>

class Triangle {
public:
    Triangle() {};
    void Init(Vertex *v0,Vertex *v1,Vertex *v2);
    
//    bool Intersect(const Ray &ray, Intersection &hit) const;
    glm::vec3 Normal;
    void UpdateNormal() {this->Normal =
        glm::normalize(glm::cross(Vtx[1]->getPos() - Vtx[0]->getPos(), Vtx[2]->getPos() - Vtx[0]->getPos()));}
	glm::vec3 GetPosition(int index) {return Vtx[index]->getPos();}
    float GetArea() {
        return 0.5f*glm::l2Norm(glm::cross(Vtx[2]->getPos() - Vtx[0]->getPos(), Vtx[1]->getPos() - Vtx[0]->getPos()));}
    void ApplyForce(glm::vec3 &force){
		glm::vec3 tempForce = force / 3.0f;
        Vtx[0]->ApplyForce(tempForce);
        Vtx[1]->ApplyForce(tempForce);
        Vtx[2]->ApplyForce(tempForce);
    }
    
    Vertex * GetVertex(int index) {return Vtx[index];}
    void ReplaceVertex(Vertex * oldVertex, Vertex * newVertex);
    bool VertexInTriangle(Vertex * v)
    {
        bool success = false;
        if (v == Vtx[0]) success = true;
        if (v == Vtx[1]) success = true;
        if (v == Vtx[2]) success = true;
        return success;
    }
    
private:
    Vertex *Vtx[3];
    float ep = 0.0001f;
};