//
//  Tetrahedron.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Tetrahedron.hpp"
#include "dsyevq3.h"
#include "dsyevh3.h"
#include <algorithm>
#include <array>

using namespace std;
using namespace glm;

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::UpdateTriangles()
{
    for(int i = 0; i < NumTriangles; i++)
        Triangles[i].UpdateNormal();
    
    vec3 TetCenter = Vertices[0]->getPos() + Vertices[1]->getPos() + Vertices[2]->getPos() + Vertices[3]->getPos();
    TetCenter *= 0.25f;
    
    // Triangle exclusion order
    if (dot(Triangles[0].Normal, TetCenter-Vertices[0]->getPos()) > 0) Triangles[0].Normal *= -1;
    if (dot(Triangles[1].Normal, TetCenter-Vertices[1]->getPos()) > 0) Triangles[1].Normal *= -1;
    if (dot(Triangles[2].Normal, TetCenter-Vertices[0]->getPos()) > 0) Triangles[2].Normal *= -1;
    if (dot(Triangles[3].Normal, TetCenter-Vertices[0]->getPos()) > 0) Triangles[3].Normal *= -1;
}

////////////////////////////////////////////////////////////////////////////////

Tetrahedron::Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3)
{
    Vertices[0] = v0;
	Vertices[1] = v1;
	Vertices[2] = v2;
	Vertices[3] = v3;
    
	v0->AddConnectedTetrahedron(this);
	v1->AddConnectedTetrahedron(this);
	v2->AddConnectedTetrahedron(this);
	v3->AddConnectedTetrahedron(this);

    Triangles[0].Init(Vertices[0], Vertices[3], Vertices[1]);
    Triangles[1].Init(Vertices[2], Vertices[1], Vertices[3]);
    Triangles[2].Init(Vertices[0], Vertices[3], Vertices[2]);
    Triangles[3].Init(Vertices[0], Vertices[2], Vertices[1]);
    
    UpdateBeta();
 
    ComputeVolume();
    
    mass = density * volume;
    
    Vertices[0]->setMass(mass / 4.0f);
    Vertices[1]->setMass(mass / 4.0f);
    Vertices[2]->setMass(mass / 4.0f);
    Vertices[3]->setMass(mass / 4.0f);
    
    UpdateTriangles();
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::UpdateBeta()
{
    glm::mat4 m = glm::mat4(
		vec4(Vertices[0]->GetMaterialCoordinate(), 1), 
		vec4(Vertices[1]->GetMaterialCoordinate(), 1),
		vec4(Vertices[2]->GetMaterialCoordinate(), 1),
		vec4(Vertices[3]->GetMaterialCoordinate(), 1)
	);
    
	beta = glm::inverse(m);
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::ComputeVolume()
{
    volume =  1.0f / 6.0f * fabs(dot(cross(Vertices[1]->GetMaterialCoordinate()-Vertices[0]->GetMaterialCoordinate(),
                                           Vertices[2]->GetMaterialCoordinate()-Vertices[0]->GetMaterialCoordinate()),
                                     Vertices[3]->GetMaterialCoordinate()-Vertices[0]->GetMaterialCoordinate()));
}

////////////////////////////////////////////////////////////////////////////////
// Uncertainities
// Don't know why traction based method doesn't yield forces that sum to zero
// Don't know why I have to multiply force by -1
void Tetrahedron::ComputeDeformationForces()
{
    Pmat = mat4x3(Vertices[0]->getPos(),
                  Vertices[1]->getPos(),
                  Vertices[2]->getPos(),
                  Vertices[3]->getPos()
                  );
    
    
    Vmat = mat4x3(Vertices[0]->getVel(),
                  Vertices[1]->getVel(),
                  Vertices[2]->getVel(),
                  Vertices[3]->getVel()
                  );
    
    mat3x4 PBeta = Pmat * beta;
    
    vec3 dx_u1 = PBeta * vec4(1,0,0,0);
    vec3 dx_u2 = PBeta * vec4(0,1,0,0);
    vec3 dx_u3 = PBeta * vec4(0,0,1,0);
    
	std::array<vec3, 3> dx{ dx_u1, dx_u2, dx_u3 };
    
    mat3x4 VBeta = Vmat * beta;
    vec3 dxd_u1 = VBeta * vec4(1,0,0,0);
    vec3 dxd_u2 = VBeta * vec4(0,1,0,0);
    vec3 dxd_u3 = VBeta * vec4(0,0,1,0);
    
	std::array<vec3, 3 > dxd{ dxd_u1, dxd_u2, dxd_u3 };
	std::array<glm::vec3, 3> dots;

    mat3 strainTensor = mat3(1.0f);
    for (int i = 0; i < 3; i++)
        for (int j = 0 ; j < 3; j++)
            strainTensor[i][j] = dot(dx[i],dx[j]) - (i == j ? 1 : 0);
   
    mat3 rateOfStrainTensor = mat3(1.0f);
    for (int i = 0; i < 3; i++)
        for (int j = 0 ; j < 3; j++)
            rateOfStrainTensor[i][j] = dot(dx[i],dxd[j]) + dot(dxd[i],dx[j]);
    
    mat3 elasticStress = mat3(1.0f);
    float strainTrace = strainTensor[0][0] + strainTensor[1][1] + strainTensor[2][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0 ; j < 3; j++)
            elasticStress[i][j] = lambda * strainTrace * (i == j ? 1 : 0) + 2 * mu * strainTensor[i][j];
    
    mat3 viscousStress = mat3(1.0f);
    float rateTrace = rateOfStrainTensor[0][0] + rateOfStrainTensor[1][1] + rateOfStrainTensor[2][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0 ; j < 3; j++)
            viscousStress[i][j] = phi * rateTrace * (i == j ? 1 : 0) + 2 * psi * rateOfStrainTensor[i][j];
    
    totalStress = mat3(1.0f);
    for (int i = 0; i < 3; i++)
        for (int j = 0 ; j < 3; j++)
            totalStress[i][j] = elasticStress[i][j] + viscousStress[i][j];

    for (int i = 0; i < Vertices.size(); i++)
    {
        vec3 forceOnNode_i = vec3(0,0,0);
        for (int j = 0; j < 4; j++)
        {
            float innerProduct = 0;
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    innerProduct += beta[k][i]*beta[l][j]*totalStress[k][l];
                }
            }
            forceOnNode_i += Pmat[j] * innerProduct;
        }
        
        forceOnNode_i *= -volume * 0.5f;
        
		// Apply force to node
        Vertices[i]->ApplyForce(forceOnNode_i);
    }
    
    //    vec3 totalForce(0,0,0);
    //    for (int i = 0; i < 4; i++)
    //        totalForce += Vertices[i]->getForce();
    
    //    std::cout << "Total force: " << l2Norm(totalForce) << std::endl;
    //    SHOWVEC(totalForce);
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::ComputeFractureForces()
{
    double sigma[3][3];
    double eigenVectors[3][3];
    double eigenValues[3];
    
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            sigma[i][j] = totalStress[j][i];
        }
    }
    
    dsyevh3(sigma, eigenVectors, eigenValues);
    
    mat3 sigmaPlus = mat3(0);
    mat3 sigmaMinus = mat3(0);
    mat3 mMat = mat3(0);
    vec3 vec1;
    
    for (int i = 0; i < 3; i++)
    {
        vec1.x = eigenVectors[0][i];
        vec1.y = eigenVectors[1][i];
        vec1.z = eigenVectors[2][i];
        
        
        ComputeMMat(vec1, mMat);
        
        sigmaPlus += fmax(0.0f,(float)eigenValues[i]) * mMat;
        sigmaMinus += fmin(0.0f,(float)eigenValues[i]) * mMat;
    }
    
    // Need to add a deformation force vector to vertex class instead of just one total force vec
    // Need to store compressive
    vec3 zero(0);
    std::array<vec3, 4> fPlus {zero,zero,zero,zero};
    std::array<vec3, 4> fMinus {zero,zero,zero,zero};
    
    for (int i = 0; i < Vertices.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            float innerProductPlus = 0;
            float innerProductMinus = 0;
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    innerProductPlus += beta[k][i]*beta[l][j]*sigmaPlus[l][k];
                    innerProductMinus += beta[k][i]*beta[l][j]*sigmaMinus[l][k];
                }
            }
            fPlus[i] += Pmat[j] * innerProductPlus;
            fMinus[i] += Pmat[j] * innerProductMinus;
        }
        fPlus[i] *= -volume * 0.5f;
        fMinus[i] *= -volume * 0.5f;
    }
    
    vec3 tempfPlus = vec3(0);
    vec3 tempfMinus = vec3(0);
    
    for (int i = 0; i < Vertices.size(); i++)
    {
        tempfPlus = vec3(fPlus[i].x, fPlus[i].y, fPlus[i].z);
        tempfMinus = vec3(fMinus[i].x, fMinus[i].y, fMinus[i].z);
        Vertices[i]->CompressiveForces[Vertices[i]->numForces] = tempfPlus;
        Vertices[i]->TensileForces[Vertices[i]->numForces] = tempfMinus;
        Vertices[i]->numForces += 1;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::UpdateMasses()
{
    ComputeVolume();
    mass = density * volume;
    Vertices[0]->setMass(Vertices[0]->getMass() + mass / 4.0f);
    Vertices[1]->setMass(Vertices[1]->getMass() + mass / 4.0f);
    Vertices[2]->setMass(Vertices[2]->getMass() + mass / 4.0f);
    Vertices[3]->setMass(Vertices[3]->getMass() + mass / 4.0f);
}

////////////////////////////////////////////////////////////////////////////////

bool Tetrahedron::TriangleInTetrahedron(Triangle * tri)
{
    for (int i = 0; i < 4; i++)
    {
        bool v0 = tri->VertexInTriangle(Triangles[i].GetVertex(0));
        bool v1 = tri->VertexInTriangle(Triangles[i].GetVertex(1));
        bool v2 = tri->VertexInTriangle(Triangles[i].GetVertex(2));
        if (v0 && v1 && v2) 
			return true;
    }
    
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::ComputeMMat(vec3 eigenVector, mat3 &outputMat)
{
    // Need to handle zero-eigen vectors
    for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
            outputMat[j][k] = eigenVector[j]*eigenVector[k];
    
    outputMat /= std::max(l2Norm(eigenVector),1.0f); // Bit of a hack
}

////////////////////////////////////////////////////////////////////////////////

void Tetrahedron::ReplaceVertex(Vertex * oldVertex, Vertex * newVertex)
{
    for (int i = 0; i < Vertices.size(); i++)
        if (Vertices[i] == oldVertex)
            Vertices[i] = newVertex;
    for (int i = 0; i < 4; i++) {
		Triangles[i].ReplaceVertex(oldVertex, newVertex); 
	}
}

////////////////////////////////////////////////////////////////////////////////

int Tetrahedron::NetConnectivity()
{
    int i0 = Vertices[0]->numConnections;
    int i1 = Vertices[1]->numConnections;
    int i2 = Vertices[2]->numConnections;
    int i3 = Vertices[3]->numConnections;
    return i0+i1+i2+i3;
}

////////////////////////////////////////////////////////////////////////////////

float Tetrahedron::ComputeSignedDistance(vec3 normal, vec3 planePoint)
{
    vec3 center = Vertices[0]->getPos() + Vertices[1]->getPos() + Vertices[2]->getPos() + Vertices[3]->getPos();
    center /= 4.0f;
    float val = dot(normal, center-planePoint);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

bool Tetrahedron::VertexInTriangle(int triIndex, Vertex * v)
{
	return Triangles[triIndex].VertexInTriangle(v);
}

////////////////////////////////////////////////////////////////////////////////

bool Tetrahedron::VertexInTetrahedron(Vertex * v)
{
	bool success = false;
	if (Vertices[0] == v) success = true;
	if (Vertices[1] == v) success = true;
	if (Vertices[2] == v) success = true;
	if (Vertices[3] == v) success = true;
	return success;
}
////////////////////////////////////////////////////////////////////////////////

//void Tetrahedron::computeiBody()
//{
//    mat3 iBody = mat3(1.0f);
//
//    float x1, x2, x3, x4;
//    float y1, y2, y3, y4;
//    float z1, z2, z3, z4;
//
//    float a, b, c, a_p, b_p, c_p;
//
//    vec3 pos0 = Vertices[0].Position;
//    vec3 pos1 = Vertices[1].Position;
//    vec3 pos2 = Vertices[2].Position;
//    vec3 pos3 = Vertices[3].Position;
//
//    x1 = pos0.x; x2 = pos1.x; x3 = pos2.x; x4 = pos3.x;
//    y1 = pos0.y; y2 = pos1.y; y3 = pos2.y; y4 = pos3.y;
//    z1 = pos0.z; z2 = pos1.z; z3 = pos2.z; z4 = pos3.z;
//
//    a = y1*y1 + y1*y2 + y2*y2 + y1*y3 + y2*y3 + y3*y3 + y1*y4 +
//    y2*y4 + y3*y4 + y4*y4 + z1*z1 + z1*z2 + z2*z2 + z1*z3 + z2*z3 +
//    z3*z3 + z1*z4 + z2*z4 + z3*z4 + z4*z4;
//
//    a *= 6 * density * volume / 60.0f;
//
//    b = x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3 + x1*x4 + x2*x4 +
//    x3*x4 + x4*x4 +z1*z1 + z1*z2 + z2*z2 + z1*z3 + z2*z3 +
//    z3*z3 + z1*z4 + z2*z4 + z3*z4 + z4*z4;
//
//    b *= 6 * density * volume / 60.0f;
//
//    c = x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3 + x1*x4 + x2*x4 +
//    x3*x4 + x4*x4 + y1*y1 + y1*y2 + y2*y2 + y1*y3 + y2*y3 + y3*y3 +
//    y1*y4 + y2*y4 + y3*y4 + y4*y4;
//
//    c *= 6 * density * volume / 60.0f;
//
//    a_p = 2 * y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2 + 2*y2*z2 + y3*z2 +
//    y4*z2 + y1*z3 + y2*z3 + 2*y3*z3 + y4*z3 + y1*z4 + y2*z4 + y3*z4 +
//    2*y4*z4;
//
//    a_p *= 6 * density * volume / 120.0f;
//
//    b_p = 2*x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2 + 2*x2*z2 + x3*z2 +
//    x4*z2 + x1*z3 + x2*z3 + 2*x3*z3 + x4*z3 + x1*z4 + x2*z4 + x3*z4 +
//    2*x4*z4;
//
//    b_p *= 6 * density * volume / 120.0f;
//
//    c_p = 2*x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2 + 2*x2*y2 + x3*y2 +
//    x4*y2 + x1*y3 + x2*y3 + 2*x3*y3 + x4*y3 + x1*y4 + x2*y4 + x3*y4 +
//    2*x4*y4;
//
//    c_p *= 6 * density * volume / 120.0f;
//
//    iBody[0][0] = a;
//    iBody[0][1] = -b_p;
//    iBody[0][2] = -c_p;
//
//    iBody[1][0] = -b_p;
//    iBody[1][1] = b;
//    iBody[1][2] = -a_p;
//
//    iBody[2][0] = -c_p;
//    iBody[2][1] = -a_p;
//    iBody[2][2] = c;
//
//    this->iBody = iBody;
//    this->iBodyInv = inverse(iBody);
//
//    SHOWMAT3(iBody);
//}
