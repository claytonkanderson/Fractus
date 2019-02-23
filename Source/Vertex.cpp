//
//  Vertex.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/17/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Vertex.hpp"

void Vertex::Update(float deltaT)
{
    glm::vec3 a = invMass * Force;
    Velocity += a * deltaT;
    Position += Velocity * deltaT;
}
