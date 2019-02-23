//
//  TetraGroupConnectivity.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 6/14/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "TetraGroupConnectivity.hpp"

using namespace glm;

void TetraGroup::SplitConnectedComponents(Vertex * suspectVert)
{
    //vector<Tetrahedron *> neighbors;
    //neighbors.resize(suspectVert->numConnections);
    //Graph G;
    //
    //for (int i = 0 ; i < suspectVert->numConnections; i++)
    //    neighbors[i] = suspectVert->tetrahedra[i];
    //
    //
    //for (int i = 0; i < neighbors.size(); i++)
    //{
    //    for (int j = i; j < neighbors.size(); j++)
    //    {
    //        if (i == j) add_edge(i, j, G);
    //        else if (TetrahedraIdentifyConnectivity(neighbors[i],neighbors[j]) == 3)
    //            add_edge(i, j, G);
    //    }
    //}
    //
    //vector<int> component(num_vertices(G));
    //int num = connected_components(G, &component[0]);
    //
    //if (num == 1) return;
    //
    //vector<vector<Tetrahedron *>> connectedComponents;
    //connectedComponents.resize(num);
    //for (int i = 0; i < component.size(); i++)
    //    connectedComponents[component[i]].push_back(neighbors[i]);
    //
    //
    //// Want all of each section of connected components to have the same verts after breaks
    //// newVert one to one with k-th connected component
    //vector<Vertex *> newVerts;
    //newVerts.resize(num - 1);
    //
    //for (int i = 0; i < newVerts.size(); i++)
    //    newVerts[i] = DuplicateVertex(suspectVert);
    //
    //for (int i = 1; i < connectedComponents.size(); i++)
    //{
    //    for (int j = 0; j < connectedComponents[i].size(); j++)
    //    {
    //        newVerts[i-1]->AddConnection(connectedComponents[i][j]);
    //        suspectVert->RemoveConnection(connectedComponents[i][j]);
    //        connectedComponents[i][j]->ReplaceVertex(suspectVert, newVerts[i-1]);
    //    }
    //}
}


void TetraGroup::ArbitrarySplit(Vertex * fractureVert, vec3 fractureNormal)
{
    vector<Vertex *> edgeVertsTop;
    vector<Vertex *> edgeVertsBot;
    
    Vertex * bottomFractureVert = DuplicateVertex(fractureVert);
    
    vector<Vertex *> fractureVertNN = GetVertexNeighbors(fractureVert);

    float planeDist = dot(fractureNormal, fractureVert->getPos());
    
    // Create new edge nodes along each edge that is intersected by fracture plane
    for (int i = 0; i < fractureVertNN.size(); i++)
    {
        vec3 ray_dir = fractureVertNN[i]->getPos() - fractureVert->getPos();
        float denom = dot(ray_dir, fractureNormal);
        if (fabs(denom) < 1e-6) continue;
        vec3 ray_origin = fractureVert->getPos();
        float t = (planeDist - dot(ray_origin, fractureNormal)) / denom;
        if (t < 1e-6 || t > 1) continue;
        
        vec3 intersectionPt = ray_origin + t * ray_dir;
        
        // Create two new verts at edge locations
        Vertex newV0;
        newV0.identifier = numVerts;
        numVerts++;
        newV0.setPos(intersectionPt);
        newV0.setMass(fractureVert->getMass()); // Will be updated later
        newV0.setVel(fractureVert->getVel()); // Not sure what this should be
        if (Vertices.size() > MAXVERTS) assert(0);
        Vertices.push_back(newV0);
        
        edgeVertsTop.push_back(&Vertices.back());
        Vertex * newV1 = DuplicateVertex(&Vertices.back());
        edgeVertsBot.push_back(newV1);
    }
    
    // Update all immediate tetrahedra neighbors of fractureVert
    // Question : How to know which edge verts to use?
    
    vector<Tetrahedron *> neighbors;
    neighbors.resize(fractureVert->numConnections);
    for(int i = 0; i < neighbors.size(); i++) neighbors[i] = fractureVert->tetrahedra[i];
    
    // Case 0 : tetrahedron is not intersected by fracture plane
    // Case 1 : tetrahedron is intersected by plane along two edges
    // Case 2 : tetrahedron is intersected by plane along one edge (plane goes through two verts)
    // Probably not going to implement Case 2 right now
    
    for (int i = 0; i < neighbors.size(); i++)
    {
        if (!PlaneIntersectTetrahedra(neighbors[i],fractureVert->getPos(), fractureNormal)) continue;
        
        std::cout << "Annoyed" << std::endl;
    }

}

vector<Vertex *> TetraGroup::GetVertexNeighbors(Vertex * vert)
{
    vector<Vertex *> neighbors;
    neighbors.reserve(4*vert->numConnections);
    
    for (int i = 0; i < vert->numConnections; i++)
        for (int j = 0 ; j < 4; j++)
            neighbors.push_back(vert->tetrahedra[i]->GetVertex(j));
    
    sort(neighbors.begin(), neighbors.end());
    auto it = unique(neighbors.begin(), neighbors.end());
    neighbors.resize( std::distance(neighbors.begin(),it) );
    it = find(neighbors.begin(), neighbors.end(), vert);
    neighbors.erase(it);
    
    return neighbors;
}
