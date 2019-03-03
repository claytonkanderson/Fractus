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