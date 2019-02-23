//
//  Triangle.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Triangle.hpp"

void Triangle::Init(Vertex *v0,Vertex *v1,Vertex *v2)
{
    Vtx[0]=v0; Vtx[1]=v1; Vtx[2]=v2;
    Normal = glm::normalize(glm::cross(Vtx[1]->getPos() - Vtx[0]->getPos(), Vtx[2]->getPos() - Vtx[0]->getPos()));
}

void Triangle::ReplaceVertex(Vertex * oldVertex, Vertex * newVertex)
{
    if (Vtx[0] == oldVertex) Vtx[0] = newVertex;
    if (Vtx[1] == oldVertex) Vtx[1] = newVertex;
    if (Vtx[2] == oldVertex) Vtx[2] = newVertex;
}
