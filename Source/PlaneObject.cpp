//
//  PlaneObject.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/17/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "PlaneObject.hpp"

using namespace std;
using namespace glm;

PlaneObject::PlaneObject(float side, vec3 center, vec3 normal)
{
    int NumTriangles = 2;
    int NumVertices = 4;

    Triangles=new Triangle[NumTriangles];
    Vertices = new Vertex[NumVertices];
    
    vec3 p1(1,0,1);
    vec3 p2(1,0,-1);
    vec3 p3(-1,0,1);
    vec3 p4(-1,0,-1);
    
    p1 *= 0.5f * side;
    p2 *= 0.5f * side;
    p3 *= 0.5f * side;
    p4 *= 0.5f * side;
    
    mat4 rotationMatrix(1.0f);
    
    if (dot(normal, vec3(0,1,0)) != dot(normal,normal))
    {
        vec3 axis = normalize(cross(normal, vec3(0,1,0)));
        float angle = acos(normal.y / l2Norm(normal));
        rotationMatrix = rotate(mat4(1.0f), angle, axis);
    }
    
    mat4 translationMatrix = translate(mat4(1.0f), center);
    
    p1 = vec3(rotationMatrix * translationMatrix * vec4(p1, 1));
    p2 = vec3(rotationMatrix * translationMatrix * vec4(p2, 1));
    p3 = vec3(rotationMatrix * translationMatrix * vec4(p3, 1));
    p4 = vec3(rotationMatrix * translationMatrix * vec4(p4, 1));
    
    Vertices[0].Set(p1, vec3(0,0,0), vec3(0,0,0));
    Vertices[1].Set(p2, vec3(0,0,0), vec3(0,0,0));
    Vertices[2].Set(p3, vec3(0,0,0), vec3(0,0,0));
    Vertices[3].Set(p4, vec3(0,0,0), vec3(0,0,0));
    
    // Normal is 1-0 x 2-0
    Triangles[0].Init(&Vertices[0], &Vertices[1], &Vertices[2]);
    Triangles[1].Init(&Vertices[3], &Vertices[2], &Vertices[1]);
    
    SetupGPU();
}

void PlaneObject::SetupGPU()
{
    std::vector<vec3> positions(6, vec3(0,0,0));
    std::vector<vec3> normals(6, vec3(0,0,0));
    
    for (int i = 0; i < 6; i+=3)
    {
        positions[i] = Triangles[i/3].GetPosition(0);
        positions[i+1] = Triangles[i/3].GetPosition(1);
        positions[i+2] = Triangles[i/3].GetPosition(2);
        
    }
    
    for (int i = 0; i < 6; i+=3)
    {
        normals[i] = Triangles[i/3].Normal;
        normals[i+1] = Triangles[i/3].Normal;
        normals[i+2] = Triangles[i/3].Normal;
    }
    
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(glm::vec3), &positions[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    vector<vec3>().swap(positions);
    vector<vec3>().swap(normals);
}

void PlaneObject::Draw(mat4 &PV)
{
    mat4 M(1.0f);
    // Update view matrices
    glUseProgram(shader);
    glUniformMatrix4fv(glGetUniformLocation(shader, "MVP"), 1, GL_FALSE, &(PV*M)[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(shader, "model"), 1, GL_FALSE, &(M)[0][0]);
    
    glBindVertexArray(VAO);
    
    // Position
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    // Normals
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
