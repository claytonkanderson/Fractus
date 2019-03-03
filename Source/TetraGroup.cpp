//
//  TetraGroup.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "TetraGroup.hpp"
#include <algorithm>
#include <array>

#include "dsyevh3.h"

using namespace glm;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ComputeDeformationForces()
{
    for (unsigned int i = 0; i < Tetrahedra.size(); i++)
        Tetrahedra[i]->ComputeDeformationForces();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ApplyGravity()
{
    glm::vec3 gravity(0,-9.8,0);
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        glm::vec3 specificGravity = gravity * Vertices[i].getMass();
        Vertices[i].ApplyForce(specificGravity);
	}
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::UpdateVertices(float deltaT, int numSteps)
{
    high_resolution_clock::time_point t1, t2;
    
    for (int j = 0; j < numSteps; j++)
    {
        //        for (int i = 0; i < Vertices.size(); i++)
        //            SHOWVEC(Vertices[i].Force);
        
        //        SHOWVEC(Vertices[0].getPos());
        t1 = high_resolution_clock::now();
        // Sets Vertex.Force
        // Force on vertex due to deformation of surrounding Tetrahedra
        // Also sets Stress Tensors for Tetrahedra
        ComputeDeformationForces();
        t2 = high_resolution_clock::now();
        timeSpan = duration_cast<duration<double>>(t2-t1);
        timers[0] += timeSpan;
        
		for (const auto & vert : Vertices)
		{
			if (isnan(vert.Force.x))
				std::cout << "Nan force" << std::endl;
		}

        t1 = high_resolution_clock::now();
        ComputeFracture();
        t2 = high_resolution_clock::now();
        timeSpan = duration_cast<duration<double>>(t2-t1);
        timers[1] += timeSpan;
        
        ApplyGravity();
        
        t1 = high_resolution_clock::now();
        for (unsigned int i = 0; i < Vertices.size(); i++)
        {
            Vertices[i].Update(deltaT);
        }
        t2 = high_resolution_clock::now();
        timeSpan = duration_cast<duration<double>>(t2-t1);
        timers[3] += timeSpan;
        
        
        t1 = high_resolution_clock::now();
        ComputeSeparation();
        t2 = high_resolution_clock::now();
        timeSpan = duration_cast<duration<double>>(t2-t1);
        timers[2] += timeSpan;
        
        //        SetDrawingColors();
        
        for (unsigned int i = 0; i < Vertices.size(); i++)
            Vertices[i].ZeroForce();
        
        UpdateMasses();
        // Check mass conservation
        //        float totalMass = 0;
        //        for (int i =0 ; i<Vertices.size(); i++)
        //            totalMass += Vertices[i].getMass();
        //
        //        std::cout << "Total mass: " << totalMass << std::endl;
        
        GroundResponse();
        
        t1 = high_resolution_clock::now();
        for (unsigned int i = 0; i < Tetrahedra.size(); i++)
            Tetrahedra[i]->UpdateTriangles();
        t2 = high_resolution_clock::now();
        timeSpan = duration_cast<duration<double>>(t2-t1);
        timers[4] += timeSpan;
        
    }
    
    t1 = high_resolution_clock::now();
    UpdateGPU();
    t2 = high_resolution_clock::now();
    timeSpan = duration_cast<duration<double>>(t2-t1);
    timers[5] += timeSpan;
    
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::GroundResponse()
{
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        if (Vertices[i].getPos().y < 0)
        {
            glm::vec3 deltaPos = Vertices[i].getPos();
            deltaPos.y = - deltaPos.y;
            Vertices[i].setPos(deltaPos);
            float elasticity = 0.9f;
            float friction = 0.1f;
            glm::vec3 velocity = Vertices[i].getVel();
            Vertices[i].setVel(glm::vec3((1-friction)*velocity.x, -elasticity * velocity.y, (1-friction)*velocity.z));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::Initialize()
{
	UpdateBetaMatrices();
	UpdateMasses();
	SetupGPU();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::CreateGrid(vec3 bottomLeftPos, vec3 topRightPos, ivec3 resolution)
{
    float deltaX = (topRightPos.x - bottomLeftPos.x) / resolution.x;
    float deltaY = (topRightPos.y - bottomLeftPos.y) / resolution.y;
    float deltaZ = (topRightPos.z - bottomLeftPos.z) / resolution.z;
    
    for (float x = bottomLeftPos.x; x <= topRightPos.x; x+=deltaX)
    {
        for (float y = bottomLeftPos.y; y <= topRightPos.y; y+=deltaY)
        {
            for (float z = bottomLeftPos.z; z <= topRightPos.z; z+=deltaZ)
            {
                Vertex v;
                v.Set(vec3(x,y,z), vec3(0), vec3(0));
                v.identifier = numVerts;
                Vertices.push_back(v);
                numVerts++;
            }
        }
    }
    
    SHOWVAR(Vertices.size());
    
    int xCount = 0;
    int yCount = 0;
    int zCount = 0;
    int mask = 0;
    
    for (float x = bottomLeftPos.x; x < topRightPos.x; x+=deltaX)
    {
        for (float y = bottomLeftPos.y; y < topRightPos.y; y+=deltaY)
        {
            for (float z = bottomLeftPos.z; z < topRightPos.z; z+=deltaZ)
            {
                mask = 0;
                if (xCount & 0x01) mask = mask | 0x01;
                if (yCount & 0x01) mask = mask | 0x02;
                if (zCount & 0x01) mask = mask | 0x04;
                AddCube(vec3(x,y,z), vec3(x+deltaX, y+deltaY, z+deltaZ), mask);
                zCount++;
            }
            yCount++;
        }
        xCount++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SetDrawingColors()
{
    float maxVal = 0;
    
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        if (length(abs(Vertices[i].Force)) > maxVal) maxVal = length(abs(Vertices[i].Force));
    }
    
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        Vertices[i].DrawForce = Vertices[i].Force / (maxVal/1.5f);
        //        Vertices[i].DrawForce = Vertices[i].Force / 100000.0f;
    }
}

////////////////////////////////////////////////////////////////////////////////

Vertex * TetraGroup::GetVertex(vec3 pos, float epsilon)
{
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        if (l2Norm(Vertices[i].getPos() - pos) < epsilon)
            return &Vertices[i];
    }
    
    return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::AddCube(vec3 bottomLeftPos, vec3 topRightPos, int mask)
{
    vec3 side = topRightPos - bottomLeftPos;
    
    vec3 s000 = bottomLeftPos;
    vec3 s100 = bottomLeftPos + vec3(side.x,0,0);
    vec3 s010 = bottomLeftPos + vec3(0,side.y,0);
    vec3 s001 = bottomLeftPos + vec3(0,0,side.z);
    vec3 s011 = bottomLeftPos + vec3(0,side.y,side.z);
    vec3 s101 = bottomLeftPos + vec3(side.x,0,side.z);
    vec3 s110 = bottomLeftPos + vec3(side.x,side.y,0);
    vec3 s111 = topRightPos;
    
    std::vector<vec3> positions;
    positions.push_back(s000);
    positions.push_back(s100);
    positions.push_back(s010);
    positions.push_back(s001);
    
    positions.push_back(s011);
    positions.push_back(s101);
    positions.push_back(s110);
    positions.push_back(s111);
    
#define swap_pos(i, j) \
{ \
vec3 temp;\
temp = positions[i];    \
positions[i] = positions[j]; \
positions[j] = temp; \
} \

    vec3 temp;
    
    // Flip X
    if (mask & 0x01)
    {
        swap_pos(2,6);
        swap_pos(4,7);
        swap_pos(0,1);
        swap_pos(3,5);
    }
    
    
    // Flip Y
    if (mask & 0x02)
    {
        swap_pos(2,0);
        swap_pos(6,1);
        swap_pos(4,3);
        swap_pos(7,5);
    }
    
    // Flip Z
    if (mask & 0x04)
    {
        swap_pos(2,4);
        swap_pos(6,7);
        swap_pos(3,0);
        swap_pos(5,1);
    }
    
    Vertex * pt;
    glm::vec3 zero(0,0,0);
    
    std::vector<Vertex *> tempVertices;
    tempVertices.resize(8);
    
    for (int i = 0; i < positions.size(); i++)
    {
        if (GetVertex(positions[i],0.001f) == NULL)
        {
            Vertex v;
            v.Set(positions[i], vec3(0), vec3(0));
            v.identifier = numVerts;
            Vertices.push_back(v);
            numVerts++;
        }
    }
    
    for (int i = 0; i < positions.size(); i++)
    {
        pt = GetVertex(positions[i], 0.001f);
        if (pt == NULL) 
			assert(0);
        else
			tempVertices[i] = pt;
    }
    
    Tetrahedron *tet1, *tet2, *tet3, *tet4, *tet5;
    
    tet1 = new Tetrahedron(tempVertices[0],tempVertices[4],tempVertices[2],tempVertices[6]);
    tet2 = new Tetrahedron(tempVertices[0],tempVertices[5],tempVertices[6],tempVertices[1]);
    tet3 = new Tetrahedron(tempVertices[7],tempVertices[5],tempVertices[4],tempVertices[6]);
    tet4 = new Tetrahedron(tempVertices[0],tempVertices[4],tempVertices[5],tempVertices[3]);
    tet5 = new Tetrahedron(tempVertices[0],tempVertices[5],tempVertices[4],tempVertices[6]);
    
    //tempVertices[0]->AddConnectedTetrahedron(tet1); tempVertices[2]->AddConnectedTetrahedron(tet1);
    //tempVertices[4]->AddConnectedTetrahedron(tet1); tempVertices[6]->AddConnectedTetrahedron(tet1);
    //
    //tempVertices[0]->AddConnectedTetrahedron(tet2); tempVertices[5]->AddConnectedTetrahedron(tet2);
    //tempVertices[6]->AddConnectedTetrahedron(tet2); tempVertices[1]->AddConnectedTetrahedron(tet2);
    //
    //tempVertices[7]->AddConnectedTetrahedron(tet3); tempVertices[5]->AddConnectedTetrahedron(tet3);
    //tempVertices[4]->AddConnectedTetrahedron(tet3); tempVertices[6]->AddConnectedTetrahedron(tet3);
    //
    //tempVertices[0]->AddConnectedTetrahedron(tet4); tempVertices[4]->AddConnectedTetrahedron(tet4);
    //tempVertices[5]->AddConnectedTetrahedron(tet4); tempVertices[3]->AddConnectedTetrahedron(tet4);
    //
    //tempVertices[0]->AddConnectedTetrahedron(tet5); tempVertices[5]->AddConnectedTetrahedron(tet5);
    //tempVertices[4]->AddConnectedTetrahedron(tet5); tempVertices[6]->AddConnectedTetrahedron(tet5);
    
    float lambda = 100;
    float mu = 100;
    float psi = 300;
    float phi = 300;
    
    tet1->SetConstants(lambda, mu, psi, phi);
    tet2->SetConstants(lambda, mu, psi, phi);
    tet3->SetConstants(lambda, mu, psi, phi);
    tet4->SetConstants(lambda, mu, psi, phi);
    tet5->SetConstants(lambda, mu, psi, phi);
    
    Tetrahedra.push_back(tet1);
    Tetrahedra.push_back(tet2);
    Tetrahedra.push_back(tet3);
    Tetrahedra.push_back(tet4);
    Tetrahedra.push_back(tet5);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SetConstants(float elasticMod, float poisson, float toughness)
{
    Toughness = toughness;
    for (unsigned int i = 0; i < Tetrahedra.size(); i++)
        Tetrahedra[i]->SetConstants(elasticMod, poisson);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ComputeFracture()
{
    for(unsigned int i = 0; i < Tetrahedra.size(); i++)
        Tetrahedra[i]->ComputeFractureForces();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SetInitialConditions(const mat4 &model, const vec3 &comVel, const vec3 &angularVel)
{
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        Vertices[i].setVel(comVel + cross(angularVel, Vertices[i].getPos()));
        Vertices[i].setPos(vec3(model * vec4(Vertices[i].getPos(),1)));
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::UpdateBetaMatrices()
{
    for (int i = 0; i < Tetrahedra.size(); i++)
        Tetrahedra[i]->UpdateBeta();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::Init(const TetraGroupInits &init)
{
    for (int i = 0; i < Tetrahedra.size(); i++)
        delete(Tetrahedra[i]);
    
    glDeleteBuffers(1,&vertexBuffer);
    glDeleteBuffers(1,&normalBuffer);
    glDeleteVertexArrays(1,&VAO);
    vector<Tetrahedron *>().swap(Tetrahedra);
    vector<Vertex>().swap(Vertices);
    
    Vertices.reserve(MAXVERTS);
    inits = init;
    
    numVerts = 0;
    CreateGrid(vec3(-0.5*inits.x_width,-0.5*inits.y_width,-0.5*inits.z_width),
               vec3(0.5*inits.x_width,0.5*inits.y_width,0.5*inits.z_width), inits.div);
    
    
    for (int i = 0; i < Vertices.size(); i++)
        assert(Vertices[i].numConnections != 0);
    
    SetInitialConditions(inits.model, inits.com_vel, inits.angular_vel);
    SetConstants(inits.elastic, inits.poisson, inits.toughness);
    
	Initialize();

    //    SHOWVEC(Vertices[0].getVel());
    //    SHOWVAR(inits.elastic);
    //    SHOWVAR(inits.toughness);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ComputeSeparation()
{
    mat3 separationTensor;
    vec3 compressiveForce;
    vec3 tensileForce;
    
    mat3 tempMat;
    
    bool wasFracture = false;
    
    if (numVerts < 0.99f*MAXVERTS)
    {
        unsigned int numOldVerts = numVerts;
        for(unsigned int i = 0; i < numOldVerts; i++)
        {
            assert(Vertices[i].numConnections != 0);
            
            separationTensor = mat3(0);
            compressiveForce = vec3(0);
            tensileForce = vec3(0);
            
            // Can move all of this to ComputeFracture
            // by adding compressiveForce accumulator, tensileForce accumulator
            // and separationTensor accumulator to vertex class
            for (int j = 0; j < Vertices[i].numForces; j++)
            {
                compressiveForce += Vertices[i].CompressiveForces[j];
                tensileForce += Vertices[i].TensileForces[j];
            }
            
            // ComputeMMat part of Tetrahedra class, nothing special about Tetrahedra[0]
            Tetrahedra[0]->ComputeMMat(compressiveForce, tempMat);
            separationTensor += -tempMat;
            
            Tetrahedra[0]->ComputeMMat(tensileForce, tempMat);
            separationTensor += tempMat;
            
            for (int j = 0; j < Vertices[i].numForces; j++)
            {
                Tetrahedra[0]->ComputeMMat(Vertices[i].CompressiveForces[j], tempMat);
                separationTensor += tempMat;
                
                Tetrahedra[0]->ComputeMMat(Vertices[i].TensileForces[j], tempMat);
                separationTensor += -tempMat;
                
            }
            
            separationTensor *= 0.5f;
            
            // End move to Compute Fracture
            
            double separationTensor2[3][3];
            
            for (int k = 0; k < 3; k++)
                for (int j = 0; j < 3; j++)
                    separationTensor2[k][j] = separationTensor[k][j];
            
            double eigenVectors[3][3];
            double eVs[3];
            dsyevh3(separationTensor2, eigenVectors, eVs);
            
            float largestEigenvalue = eVs[0];
            vec3 principalEigenVector = vec3(eigenVectors[0][0], eigenVectors[0][1], eigenVectors[0][2]);
            if (largestEigenvalue < eVs[1]) { largestEigenvalue = eVs[1];
                principalEigenVector = vec3(eigenVectors[1][0], eigenVectors[1][1], eigenVectors[1][2]);
            }
            if (largestEigenvalue < eVs[2]) { largestEigenvalue = eVs[2];
                principalEigenVector = vec3(eigenVectors[2][0], eigenVectors[2][1], eigenVectors[2][2]);
            }
            
            
            if (largestEigenvalue < 0) continue;
            if (largestEigenvalue > Toughness){
                wasFracture = true;
//                std::cout << "largestEigenvalue: " << largestEigenvalue << std::endl;
                Vertex * fracturingVertex = &Vertices[i];
                assert(fracturingVertex->numConnections != 0);
                
                if (!ValidFractureVert(fracturingVertex)) continue;
                vec3 planeNormal = FindBestFracturePlane(fracturingVertex, principalEigenVector);
                vec3 planePoint = fracturingVertex->getPos();
                vector<Tetrahedron *> top;
                vector<Tetrahedron *> bottom;
                
                for (int j = 0; j < fracturingVertex->numConnections; j++)
                {
//                    if (PlaneIntersectTetrahedra(fracturingVertex->tetrahedra[j], planePoint, planeNormal))
//                    {
                    
                        if (fracturingVertex->tetrahedra[j]->ComputeSignedDistance(planeNormal, planePoint) >= 0)
                            top.push_back(fracturingVertex->tetrahedra[j]);
                        else
                            bottom.push_back(fracturingVertex->tetrahedra[j]);
                    }
                
                
                if (top.size() == 0) {
                    std::cout << "Top has size 0" << std::endl;
                    continue;
                }
                
                if (bottom.size() == 0) {
                    std::cout << "Bottom has size 0" << std::endl;
                    continue;
                }
                
//                for (int j = 0; j < Vertices.size(); j++)
//                {
//                    assert(Vertices[j].numConnections != 0);
//                }
                
                
                // The problem is if none of the tet's in bot
                // have a shared face to facture on,
                // this vert won't get any connections
                Vertex * bottomVert = nullptr;
                if (bottom.size() != 0)
                {
                    if (CheckFractureFaces(bottom, fracturingVertex))
                    {
                        bottomVert = DuplicateVertex(fracturingVertex);
                        
                        // Fracture Bottom using bottomVert
                        for (int j = 0; j < bottom.size(); j++)
                        {
                            //                    FractureTetrahedra(bottom[j], fracturingVertex, bottomVert, planeNormal, true);
                            assert(bottomVert != nullptr);
                            bottomVert->AddConnection(bottom[j]);
                            fracturingVertex->RemoveConnection(bottom[j]);
                            bottom[j]->ReplaceVertex(fracturingVertex, bottomVert);
                            
                        }
                        
                    }
                }
                
//                for (int l = 0; l < Vertices.size(); l++)
//                    assert(Vertices[l].numConnections != 0);
                
                for (int j = 0; j < top.size(); j++)
                {
                    if (IsIsolated(top[j])) IsolateTetrahedra(top[j]);
                }
                
                
//                for (int l = 0; l < Vertices.size(); l++)
//                    assert(Vertices[l].numConnections != 0);
                
                for (int j = 0; j < bottom.size(); j++)
                {
                    if (IsIsolated(bottom[j])) IsolateTetrahedra(bottom[j]);
                }
                
//                for (int l = 0; l < Vertices.size(); l++)
//                    assert(Vertices[l].numConnections != 0);
            }
        }
    }
    
    if (addendumFrames == 45) {
        ConnectivityAddendum();
        addendumFrames = 0; }
    addendumFrames++;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::FractureTetrahedra(Tetrahedron * tet, Vertex * oldVertex, Vertex * newVertex, const vec3 &planeNormal, bool planeSnap)
{
    // If snapping to a plane, only fracture on that triangle if there exists another tetrahedra
    // that shares the triangle
    
    if (planeSnap)
    {
        float maxFaceDot = -1;
        int triIndex = -1;
        
        // Find Fracture Face
        for (int i = 0; i < 4; i++)
        {
            Triangle * tri = &tet->Triangles[i];
            if (tri->VertexInTriangle(oldVertex))
            {
                int triCounter = 0;
                // Check that another tetrahedra shares this face
                for (int j = 0; j < oldVertex->numConnections; j++)
                    if (oldVertex->tetrahedra[j]->TriangleInTetrahedron(tri)) triCounter++;
                
                // Only allow fracturing on shared faces
                if (triCounter == 2)
                {
                    float tempDot = dot(planeNormal, tri->Normal);
                    tempDot *= tempDot;
                    if(tempDot > maxFaceDot)
                    {
                        triIndex = i;
                        maxFaceDot = tempDot;
                    }
                }
            }
        }
        
        if (triIndex < 0) return;
        assert(maxFaceDot > 0);
        
        Triangle * FractureFace = &tet->Triangles[triIndex];
        
        vec3 normal = FractureFace->Normal;
        vec3 planePoint = FractureFace->GetPosition(0);
        
        Vertex * oldV0 = oldVertex;
        Vertex * oldV1 = FractureFace->GetVertex(1);
        Vertex * oldV2 = FractureFace->GetVertex(2);
        
        if (oldV1 == oldV0) oldV1 = FractureFace->GetVertex(0);
        else if (oldV2 == oldV0) oldV2 = FractureFace->GetVertex(0);
        
        assert(oldV0 == oldVertex);
        assert(oldV0 != oldV1);
        assert(oldV0 != oldV2);
        assert(oldV1 != oldV2);
        
        Vertex * newV0p = newVertex;
        Vertex * newV1p = DuplicateVertex(oldV1);
        Vertex * newV2p = DuplicateVertex(oldV2);
        
        // Also need to update vertex connectivity
        if (oldV0 != newV0p) VertexSplit(oldV0, newV0p, normal, planePoint);
        VertexSplit(oldV1, newV1p, normal, planePoint);
        VertexSplit(oldV2, newV2p, normal, planePoint);
        
        //        for (int i = 0; i < Vertices.size(); i++)
        //            assert(Vertices[i].numConnections != 0);
        
        //        for (int i = 0; i < oldV0->numConnections; i++)
        //        {
        //            assert(oldV0->tetrahedra[i]->VertexInTetrahedron(oldV0));
        //            assert(!(oldV0->tetrahedra[i]->VertexInTetrahedron(newV0p)));
        //        }
        
        assert(oldV0->numConnections != 0);
        assert(oldV1->numConnections != 0);
        assert(oldV2->numConnections != 0);
        assert(newV0p->numConnections != 0);
        assert(newV1p->numConnections != 0);
        assert(newV2p->numConnections != 0);
        
        for (int i = 0; i < oldV1->numConnections; i++)
        {
            assert(oldV1->tetrahedra[i]->VertexInTetrahedron(oldV1));
            assert(!(oldV1->tetrahedra[i]->VertexInTetrahedron(newV1p)));
        }
        
        for (int i = 0; i < oldV2->numConnections; i++)
        {
            assert(oldV2->tetrahedra[i]->VertexInTetrahedron(oldV2));
            assert(!(oldV2->tetrahedra[i]->VertexInTetrahedron(newV2p)));
        }
        
        //        for (int i = 0; i < newV0p->numConnections; i++)
        //        {
        //            assert(newV0p->tetrahedra[i]->VertexInTetrahedron(newV0p));
        //            assert(!(newV0p->tetrahedra[i]->VertexInTetrahedron(oldV0)));
        //        }
        
        for (int i = 0; i < newV1p->numConnections; i++)
        {
            assert(newV1p->tetrahedra[i]->VertexInTetrahedron(newV1p));
            assert(!(newV1p->tetrahedra[i]->VertexInTetrahedron(oldV1)));
        }
        
        for (int i = 0; i < newV2p->numConnections; i++)
        {
            assert(newV2p->tetrahedra[i]->VertexInTetrahedron(newV2p));
            assert(!(newV2p->tetrahedra[i]->VertexInTetrahedron(oldV2)));
        }
    }
    else
    {
        // Create vector of neighbor tetrahedra that will need updating (are also intersected by plane)
        // (Not including itself)
        
        //        vector<Tetrahedron *> tets;
        //        tets.reserve(4*(FractureVertex->maxConnectivity));
        //
        //        for (int i = 0; i < 4; i++)
        //        {
        //            for (int j = 0; j < tet->GetVertex(i)->numConnections; j++)
        //            {
        //                tets.push_back(tet->GetVertex(i)->tetrahedra[j]);
        //            }
        //        }
        //
        //        // Remove duplicates & tet itself
        //        sort(tets.begin(), tets.end());
        //        auto it = unique(tets.begin(), tets.end());
        //        tets.resize( std::distance(tets.begin(),it) );
        //        it = find(tets.begin(), tets.end(), tet);
        //        tets.erase(it);
        //
        //        // Remove those that are not intersected by the fracture plane
        //        for (auto it2 = tets.begin(); it2 != tets.end(); it2++)
        //        {
        //            if (!PlaneIntersectTetrahedra(*it2, FractureVertex->getPos(), planeNormal))
        //                tets.erase(it2);
        //        }
        //
        //        // Remaining tets will have to be updated after tet has been fractured
        //
        //        // Calculate new vertex locations
        //        vector<vec3> vertexPositions;
        //
        //        for (int i = 0; i < 4; i++)
        //            if (tet->GetVertex(i) != FractureVertex)
        //                vertexPositions.push_back(tet->GetVertex(i)->getPos());
        //
        //        vec3 edge0 = vertexPositions[0] - vertexPositions[1];
        //        vec3 edge1 = vertexPositions[1] - vertexPositions[2];
        //        vec3 edge2 = vertexPositions[2] - vertexPositions[0];
        //
        //        // Plane should intersect two of these
        //
        //
        //        // Create new vertices, tetrahedra
        //        // Ensure Mesh consistency
        //
        
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::IsolateTetrahedra(Tetrahedron * tet)
{
    vector<Tetrahedron *> neighbors;
    
    for (int i = 0 ; i < 4; i ++)
        for (int j = 0; j < tet->GetVertex(i)->numConnections; j++)
            neighbors.push_back(tet->GetVertex(i)->tetrahedra[j]);
    
    std::sort(neighbors.begin(), neighbors.end());
    auto it = std::unique(neighbors.begin(), neighbors.end());
    neighbors.resize( std::distance(neighbors.begin(),it) );
    it = find(neighbors.begin(), neighbors.end(), tet);
    neighbors.erase(it);
    
    for (int i = 0; i < neighbors.size(); i++)
    {
        int con = TetrahedraIdentifyConnectivity(tet, neighbors[i]);
        if (con == 1)
        {
            SplitTetrahedraByVertex(tet, neighbors[i]);
        } else if (con == 2)
        {
            SplitTetrahedraByEdge(tet, neighbors[i]);
        } else if (con == 3) { assert(0); }
    }
    
}

////////////////////////////////////////////////////////////////////////////////

bool TetraGroup::IsIsolated(Tetrahedron * tet)
{
    vector<Tetrahedron *> neighbors;
    
    for (int i = 0 ; i < 4; i ++)
        for (int j = 0; j < tet->GetVertex(i)->numConnections; j++)
            neighbors.push_back(tet->GetVertex(i)->tetrahedra[j]);
    
    sort(neighbors.begin(), neighbors.end());
    auto it = unique(neighbors.begin(), neighbors.end());
    neighbors.resize( std::distance(neighbors.begin(),it) );
    it = find(neighbors.begin(), neighbors.end(), tet);
    neighbors.erase(it);
    
    for (int i = 0; i < neighbors.size(); i++)
    {
        if (TetrahedraIdentifyConnectivity(tet, neighbors[i]) == 3) return false;
    }
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SplitTetrahedraByVertex(Tetrahedron * tet, Vertex * oldVertex, Vertex * newVertex)
{
    assert(tet->VertexInTetrahedron(oldVertex));
    
    vec3 normal;
    vec3 planePoint = oldVertex->getPos();
    
    if (tet->Triangles[0].VertexInTriangle(oldVertex)) normal = tet->Triangles[0].Normal;
    else normal = tet->Triangles[1].Normal;
    
    VertexSplit(oldVertex, newVertex, normal, planePoint);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SplitTetrahedraByVertex(Tetrahedron * loneTet, Tetrahedron * otherTet)
{
    assert(TetrahedraIdentifyConnectivity(loneTet, otherTet) == 1);
    Vertex * oldVert = TetrahedraMutualVertex(loneTet, otherTet);
    Vertex * newVert = DuplicateVertex(oldVert);
    loneTet->ReplaceVertex(oldVert, newVert);
    oldVert->RemoveConnection(loneTet);
    newVert->AddConnection(loneTet);
    
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SplitTetrahedraByEdge(Tetrahedron * loneTet, Tetrahedron * otherTet)
{
    assert(TetrahedraIdentifyConnectivity(loneTet, otherTet) == 2);
    
    auto verts = TetrahedraMutualVertices(loneTet, otherTet);
    Vertex * oldV0 = verts.first;
    Vertex * oldV1 = verts.second;
    
    Vertex * newV0 = DuplicateVertex(oldV0);
    Vertex * newV1 = DuplicateVertex(oldV1);
    
    loneTet->ReplaceVertex(oldV0, newV0);
    loneTet->ReplaceVertex(oldV1, newV1);
    
    oldV0->RemoveConnection(loneTet);
    oldV1->RemoveConnection(loneTet);
    newV0->AddConnection(loneTet);
    newV1->AddConnection(loneTet);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SplitTetrahedraByEdge(Tetrahedron * tet, Vertex * oldV0, Vertex * oldV1, Vertex * newV0, Vertex * newV1)
{
    assert(tet->VertexInTetrahedron(oldV0) && tet->VertexInTetrahedron(oldV1));
    vec3 normal;
    vec3 planePoint = oldV0->getPos();
    for (int i = 0; i < 4; i++)
    {
        if (tet->Triangles[i].VertexInTriangle(oldV0) && tet->Triangles[i].VertexInTriangle(oldV1))
            normal = tet->Triangles[i].Normal;
    }
    VertexSplit(oldV0, newV0, normal, planePoint);
    VertexSplit(oldV1, newV1, normal, planePoint);
}

////////////////////////////////////////////////////////////////////////////////

vec3 TetraGroup::FindBestFracturePlane(Vertex * fracturingVertex, vec3 eigenDir)
{
    vec3 planeNormal = vec3(0,0,0);
    float maxFaceDot = -1;
    
    
    for (int i = 0; i < fracturingVertex->numConnections; i++)
    {
        Tetrahedron * tet = fracturingVertex->tetrahedra[i];
        for (int i = 0; i < 4; i++)
        {
            Triangle * tri = &tet->Triangles[i];
            if (tri->VertexInTriangle(fracturingVertex))
            {
                int triCounter = 0;
                // Check that another tetrahedra shares this face
                for (int j = 0; j < fracturingVertex->numConnections; j++)
                    if (fracturingVertex->tetrahedra[j]->TriangleInTetrahedron(tri)) triCounter++;
                
                // Only allow fracturing on shared faces
                if (triCounter == 2)
                {
                    float tempDot = dot(eigenDir, tri->Normal);
                    tempDot *= tempDot;
                    if(tempDot > maxFaceDot)
                    {
                        planeNormal = tri->Normal;
                        maxFaceDot = tempDot;
                    }
                }
            }
        }
    }
    // Tried to fracture where there were no faces to fracture with
    assert(maxFaceDot > 0);
    return planeNormal;
}

////////////////////////////////////////////////////////////////////////////////

bool TetraGroup::ValidFractureVert(Vertex * fracturingVertex)
{
    // At least one tet should share a face to fracture upon for this vertex to be able to fracture
    std::vector<Tetrahedron *> neighbors; 
	neighbors.resize(fracturingVertex->numConnections);
    for (int i = 0; i < fracturingVertex->numConnections; i++) 
		neighbors[i] = fracturingVertex->tetrahedra[i];
    
    return CheckFractureFaces(neighbors, fracturingVertex);
}

////////////////////////////////////////////////////////////////////////////////

bool TetraGroup::CheckFractureFaces(const std::vector<Tetrahedron *> &tets, Vertex * vert)
{
    for(int k = 0; k < tets.size(); k++)
    {
        for (int i = 0; i < 4; i++)
        {
            Triangle * tri = &tets[k]->Triangles[i];
            if (tri->VertexInTriangle(vert))
            {
                int triCounter = 0;
                // Check that another tetrahedra shares this face
                for (int j = 0; j < vert->numConnections; j++)
                    if (vert->tetrahedra[j]->TriangleInTetrahedron(tri)) triCounter++;
                
                // Only allow fracturing on shared faces
                if (triCounter == 2) return true;
            }
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::VertexSplit(Vertex * oldVert, Vertex * newVert, vec3 &planeNormal, vec3 &planePoint)
{
    for (int k = 0; k < oldVert->numConnections; k++)
    {
        Tetrahedron * tempTet = oldVert->tetrahedra[k];
        
        // If greater than 0, keep old vertices
        // So if less than 0, replace vertices
        if (tempTet->ComputeSignedDistance(planeNormal, planePoint) < 0)
        {
            oldVert->RemoveConnection(tempTet);
            newVert->AddConnection(tempTet);
            tempTet->ReplaceVertex(oldVert, newVert); // oldV0 should not be in tempTet anymore
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetraGroup::PlaneIntersectTetrahedra(Tetrahedron * tet, const vec3 &planePoint, const vec3 &planeNormal)
{
	std::array<float, 4> d;
    
    // Positive if point is pointed to by the normal
    for (int i = 0; i < 4; i++)
        d[i] = -dot(planePoint - tet->GetVertex(i)->getPos(), planeNormal);
    
    // If any pair has different sign, then tetrahedra is intersected by plane
    // Probably a terrible way to do this check
    if (d[0]*d[1] < 0) return true;
    if (d[0]*d[2] < 0) return true;
    if (d[0]*d[3] < 0) return true;
    if (d[1]*d[2] < 0) return true;
    if (d[1]*d[3] < 0) return true;
    if (d[2]*d[3] < 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::UpdateMasses()
{
    for (int i = 0; i < Vertices.size(); i++)
        Vertices[i].setMass(0);
    
    for(int i = 0; i < Tetrahedra.size(); i++)
        Tetrahedra[i]->UpdateMasses();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ConnectivityAddendum()
{
    for (int i = 0; i < Tetrahedra.size(); i++)
    {
        if (IsIsolated(Tetrahedra[i])) IsolateTetrahedra(Tetrahedra[i]);
    }
    
    
    for (unsigned int i = 0; i < Vertices.size(); i++)
    {
        SplitConnectedComponents(&Vertices[i]);
    }
//
//    for (unsigned int i = 0; i < Tetrahedra.size(); i++)
//    {
//        Vertex * v0 = Tetrahedra[i]->GetVertex(0);
//        Vertex * v1 = Tetrahedra[i]->GetVertex(1);
//        Vertex * v2 = Tetrahedra[i]->GetVertex(2);
//        Vertex * v3 = Tetrahedra[i]->GetVertex(3);
//        
//        auto comp = [](Vertex * a, Vertex * b) {return a->numConnections < b->numConnections; };
//        
//        std::vector<Vertex *> verts;
//        verts.push_back(v0);
//        verts.push_back(v1);
//        verts.push_back(v2);
//        verts.push_back(v3);
//        std::make_heap(verts.begin(), verts.end(), comp);
//        
//        //        std::cout << "front: " << verts.front()->numConnections << std::endl;
//        //        std::cout << "back: " << verts.back()->numConnections << std::endl;
//        
//        int shouldBreak = 0;
//        
//        if (v0->numConnections == 1) shouldBreak++;
//        if (v1->numConnections == 1) shouldBreak++;
//        if (v2->numConnections == 1) shouldBreak++;
//        if (v3->numConnections == 1) shouldBreak++;
//        //            SHOWVAR(shouldBreak);
//        // Connected by edge
//        if (shouldBreak == 2)
//        {
//            Vertex * dangler1 = verts.front();
//            std::pop_heap(verts.begin(), verts.end(), comp);
//            Vertex * dangler2 = verts.front();
//            
//            Vertex * newVert1 = DuplicateVertex(dangler1);
//            newVert1->AddConnection(Tetrahedra[i]);
//            
//            Vertex * newVert2 = DuplicateVertex(dangler2);
//            newVert2->AddConnection(Tetrahedra[i]);
//            
//            dangler1->RemoveConnection(Tetrahedra[i]);
//            dangler2->RemoveConnection(Tetrahedra[i]);
//            
//            Tetrahedra[i]->ReplaceVertex(dangler1, newVert1);
//            Tetrahedra[i]->ReplaceVertex(dangler2, newVert2);
//        }
//        
//        // Hanging on by one vertex
//        else if (shouldBreak == 3)
//        {
//            Vertex * dangler;
//            if (v0->numConnections != 1) dangler = v0;
//            else if (v1->numConnections != 1) dangler = v1;
//            else if (v2->numConnections != 1) dangler = v2;
//            else if (v3->numConnections != 1) dangler = v3;
//            
//            Vertex * newVert = DuplicateVertex(dangler);
//            newVert->AddConnection(Tetrahedra[i]);
//            dangler->RemoveConnection(Tetrahedra[i]);
//            Tetrahedra[i]->ReplaceVertex(dangler, newVert);
//        }
//    }
    //
    //    Tetrahedron * tet0, *tet1;
    //
    //    for (int i = 0; i < Tetrahedra.size(); i++)
    //    {
    //        tet0 = Tetrahedra[i];
    //        for (int j = i + 1; j < Tetrahedra.size(); j++)
    //        {
    //            tet1 = Tetrahedra[j];
    //
    //            int connectivityNum = TetrahedraIdentifyConnectivity(tet0, tet1);
    //
    //            switch (connectivityNum) {
    //                    // Not connected, no updating needed
    //                case 0:
    //                    break;
    //
    //                    // Connected by a vertex, should find and separate them
    //                case 1:
    //                {
    //                    std::cout << "Case 1" << std::endl;
    //
    //                    Vertex * loneVert = TetrahedraMutualVertex(tet0, tet1);
    //
    //                    vec3 normal, point;
    //
    //                    if (tet0->VertexInTriangle(0, loneVert)) {
    //                        normal = tet0->Triangles[0].Normal;
    //                        point = loneVert->getPos();
    //                    } else {
    //                        normal = tet0->Triangles[1].Normal;
    //                        point = loneVert->getPos();
    //                    }
    //
    //                    Vertex * newV0p = DuplicateVertex(loneVert);
    //                    VertexSplit(loneVert, newV0p, normal, point);
    //
    //                    break;
    //                }
    //
    //                    // Connected by an edge, should find and separate them
    //                case 2:
    //                {
    ////                    std::cout << "Case 2" << std::endl;
    ////                    Vertex * oldV0, * oldV1;
    ////                    pair<Vertex *, Vertex *> VertexPair = TetrahedraMutualVertices(tet0, tet1);
    ////                    oldV0 = VertexPair.first;
    ////                    oldV1 = VertexPair.second;
    ////
    ////                    Triangle * tri = TetrahedraMutualEdge(tet0, tet1);
    ////                    vec3 normal = tri->Normal;
    ////                    vec3 point0 = oldV0->getPos();
    ////                    vec3 point1 = oldV1->getPos();
    ////
    ////                    Vertex * newV0 = DuplicateVertex(oldV0);
    ////                    Vertex * newV1 = DuplicateVertex(oldV1);
    ////
    ////                    VertexSplit(oldV0, newV0, normal, point0);
    ////                    VertexSplit(oldV1, newV1, normal, point1); // Not sure this necessary or correct
    //
    //                    break;
    //                }
    //                    // Connected by a face, no updating needed
    //                case 3:
    //                    break;
    //
    //                case 4:
    //                    assert(0);
    //                    break;
    //
    //                default:
    //                    assert(0);
    //                    break;
    //            }
    //        }
    //    }
    //
}

////////////////////////////////////////////////////////////////////////////////

pair<Vertex *,Vertex *> TetraGroup::TetrahedraMutualVertices(Tetrahedron *tet0, Tetrahedron *tet1)
{
	Vertex * firstDuplicate = nullptr;
	Vertex * secondDuplicate = nullptr;

	for (auto v0 : tet0->Vertices)
	{
		for (auto v1 : tet1->Vertices)
		{
			if (!firstDuplicate && v0 == v1)
				firstDuplicate = v0;
			else if (v0 == v1)
			{
				secondDuplicate = v0;
				return { firstDuplicate, secondDuplicate };
			}
		}
	}
    
	return {nullptr, nullptr};
}

////////////////////////////////////////////////////////////////////////////////

Vertex * TetraGroup::DuplicateVertex(Vertex * oldVert)
{
    Vertex newV0;
    newV0.identifier = numVerts;
    numVerts++;
    newV0.setPos(oldVert->getPos());
    newV0.setMass(oldVert->getMass());
    newV0.setVel(oldVert->getVel());
    if (Vertices.size() > MAXVERTS) assert(0);
    Vertices.push_back(newV0);
    return &Vertices.back();
}

////////////////////////////////////////////////////////////////////////////////

Vertex * TetraGroup::TetrahedraMutualVertex(Tetrahedron *tet0, Tetrahedron *tet1)
{
	for (auto v0 : tet0->Vertices)
	{
		for (auto v1 : tet1->Vertices)
		{
			if (v0 == v1)
				return v0;
		}
	}

	return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

//Triangle * TetraGroup::TetrahedraMutualEdge(Tetrahedron *tet0, Tetrahedron *tet1)
//{
//	auto pair = TetrahedraMutualVertices(tet0, tet1);
//
//	if (!pair.first || !pair.second)
//		return nullptr;
//
//    vector<Vertex *> v0s = tet0->Vertices;
//    vector<Vertex *> v1s = tet1->Vertices;
//    vector<Vertex *> setIntersection; setIntersection.resize(3);
//    sort(v0s.begin(), v0s.end());
//    sort(v1s.begin(), v1s.end());
//    
//    auto it = set_intersection(v0s.begin(), v0s.end(), v1s.begin(), v1s.end(), setIntersection.begin());
//    setIntersection.resize(it-setIntersection.begin());
//    
//    if (setIntersection.size() != 2) assert(0);
//    
//    Vertex * v0 = setIntersection[0];
//    Vertex * v1 = setIntersection[1];
//    
//    Triangle * tri = nullptr;
//    
//    if (tet0->VertexInTriangle(0, v0) && tet0->VertexInTriangle(0, v1))
//    {
//        tri = &tet0->Triangles[0];
//        return tri;
//    }
//    
//    if (tet0->VertexInTriangle(1, v0) && tet0->VertexInTriangle(1, v1))
//    {
//        tri = &tet0->Triangles[1];
//        return tri;
//    }
//    
//    if (tet0->VertexInTriangle(2, v0) && tet0->VertexInTriangle(2, v1))
//    {
//        tri = &tet0->Triangles[2];
//        return tri;
//    }
//    
//    // Two triangles have the edge we want, code should never reach here
//    
//    assert(0);
//    return tri;
//}

////////////////////////////////////////////////////////////////////////////////

int TetraGroup::TetrahedraIdentifyConnectivity(Tetrahedron * tet1, Tetrahedron * tet2)
{
    // If two Tetrahedra share no verts, they aren't connected
    // If two Tetrahedra share one vert, they are connected by one vert
    // If two Tetrahedra share one vert, they are connected by an edge
    // If two Tetrahedra share one vert, they are connected by a face
    
    int count = 0;
    Vertex * v00 = tet1->GetVertex(0);
    Vertex * v01 = tet1->GetVertex(1);
    Vertex * v02 = tet1->GetVertex(2);
    Vertex * v03 = tet1->GetVertex(3);
    
    Vertex * v10 = tet2->GetVertex(0);
    Vertex * v11 = tet2->GetVertex(1);
    Vertex * v12 = tet2->GetVertex(2);
    Vertex * v13 = tet2->GetVertex(3);
    
    if (v00 == v10 || v00 == v11 || v00 == v12 || v00 == v13) count++;
    if (v01 == v10 || v01 == v11 || v01 == v12 || v01 == v13) count++;
    if (v02 == v10 || v02 == v11 || v02 == v12 || v02 == v13) count++;
    if (v03 == v10 || v03 == v11 || v03 == v12 || v03 == v13) count++;
    
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::ZeroForces()
{
    for(unsigned int i = 0; i < Vertices.size(); i++)
        Vertices[i].ZeroForce();
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::SetupGPU()
{
    int numTriangles = 4* (int)Tetrahedra.size();
    std::vector<vec3> positions(3* numTriangles);
    std::vector<vec3> normals(3* numTriangles);
    
    int local_count = 0;
    
    for (int j = 0; j < Tetrahedra.size(); j++)
    {
        Tetrahedron * tet = Tetrahedra[j];
        for (int i = 0; i < 3*4; i+=3)
        {
            positions[local_count] = tet->Triangles[i/3].GetPosition(0);
            positions[local_count+1] = tet->Triangles[i/3].GetPosition(1);
            positions[local_count+2] = tet->Triangles[i/3].GetPosition(2);
            normals[local_count] = tet->Triangles[i/3].Normal;
            normals[local_count+1] = tet->Triangles[i/3].Normal;
            normals[local_count+2] = tet->Triangles[i/3].Normal;
            local_count += 3;
            
        }
    }
    
    mat.ambient = vec3(0.5,0.5,0.5);
    mat.diffuse = vec3(0.5,0.5,0.2);
    mat.specular = vec3(0.5,0.5,0.2);
    mat.shininess = 1.0f;
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &vertexBuffer);
    glGenBuffers(1, &normalBuffer);
    
    glBindVertexArray(VAO);
    
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(glm::vec3), &positions[0], GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STREAM_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::UpdateGPU()
{
    int NumTriangles = 4* (int)Tetrahedra.size();
    std::vector<vec3> positions(3*NumTriangles);
    std::vector<vec3> normals(3*NumTriangles);
    
    int local_count = 0;
    
    for (int j = 0; j < Tetrahedra.size(); j++)
    {
        Tetrahedron * tet = Tetrahedra[j];
        for (int i = 0; i < 3*4; i+=3)
        {
            positions[local_count] = tet->Triangles[i/3].GetPosition(0);
            positions[local_count+1] = tet->Triangles[i/3].GetPosition(1);
            positions[local_count+2] = tet->Triangles[i/3].GetPosition(2);
            normals[local_count] = tet->Triangles[i/3].Normal;
            normals[local_count+1] = tet->Triangles[i/3].Normal;
            normals[local_count+2] = tet->Triangles[i/3].Normal;
            local_count += 3;
            
        }
    }
    
    glUniform3f(glGetUniformLocation(shader, "material.ambient"), mat.ambient.x, mat.ambient.y, mat.ambient.z);
    glUniform3f(glGetUniformLocation(shader, "material.specular"), mat.specular.x, mat.specular.y, mat.specular.z);
    glUniform3f(glGetUniformLocation(shader, "material.diffuse"), mat.diffuse.x, mat.diffuse.y, mat.diffuse.z);
    glUniform1f(glGetUniformLocation(shader, "material.shininess"), mat.shininess);
    
    glBindVertexArray(VAO);
    
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(glm::vec3), &positions[0], GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STREAM_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    vector<vec3>().swap(positions);
    vector<vec3>().swap(normals);
}

////////////////////////////////////////////////////////////////////////////////

void TetraGroup::Draw(mat4 &PV)
{
    if (meshView) glPolygonMode( GL_FRONT_AND_BACK, GL_LINE ); // Wireframe mode
    else glPolygonMode( GL_FRONT_AND_BACK, GL_FILL ); // Normal Shading
    
    int NumTriangles = 4* (int)Tetrahedra.size();
    // Update view matrices
    mat4 M = mat4(1.0f);
    glUseProgram(shader);
    glUniformMatrix4fv(glGetUniformLocation(shader, "MVP"), 1, GL_FALSE, &(PV*M)[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(shader, "model"), 1, GL_FALSE, &(M)[0][0]);
    
    glBindVertexArray(VAO);
    
    // Position
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    
    // Normals
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    
    glDrawArrays(GL_TRIANGLES, 0, 3*NumTriangles);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
}

////////////////////////////////////////////////////////////////////////////////