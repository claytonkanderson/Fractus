//
//  Vertex.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/17/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Vertex.hpp"

////////////////////////////////////////////////////////////////////////////////

Vertex::Vertex()
{
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
	ZeroForce();
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::Set(const glm::vec3 & p, const glm::vec3 & v, const glm::vec3 & f)
{
	Position = p; 
	Velocity = v; 
	Force = f;
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::ApplyForce(glm::vec3 & force)
{
	Force += force;
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::ZeroForce()
{
	Force = glm::vec3(0);
	numForces = 0;

	for (int i = 0; i < maxConnectivity; i++)
	{
		CompressiveForces[i] = glm::vec3(0);
		TensileForces[i] = glm::vec3(0);
	}
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::Update(float deltaT)
{
    glm::vec3 a = invMass * Force;
    Velocity += a * deltaT;
    Position += Velocity * deltaT;
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::SetCompressiveForce(glm::vec3 force)
{
	if (numForces < maxConnectivity)
	{
		CompressiveForces[numForces] = force;
	}
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::SetTensileForce(glm::vec3 force)
{
	if (numForces < maxConnectivity)
	{
		TensileForces[numForces] = force;
		numConnections++;
	}
}

////////////////////////////////////////////////////////////////////////////////

bool Vertex::AddConnectedTetrahedron(Tetrahedron * tet)
{
	if (numConnections < maxConnectivity)
	{
		tetrahedra[numConnections] = tet;
		numConnections++;
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////

void Vertex::RemoveConnection(Tetrahedron * tet)
{
	for (size_t i = 0; i < tetrahedra.size(); i++)
	{
		if (tetrahedra[i] == tet)
		{
			tetrahedra.erase(tetrahedra.begin() + i);
			numConnections--;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

bool Vertex::AddConnection(Tetrahedron * tet)
{
	if (numConnections < maxConnectivity)
	{
		tetrahedra[numConnections] = tet;
		numConnections++;
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////