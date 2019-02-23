//
//  PlaneObject.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/17/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "Object.hpp"
#include "Triangle.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtx/norm.hpp>
#include <GL/glew.h>

class PlaneObject: public Object
{
public:
    void Draw(glm::mat4 &PV);
    PlaneObject(float size, glm::vec3 center, glm::vec3 normal);
    void SetShader(GLuint shader) {this->shader = shader;}
    ~PlaneObject() {}
private:
    Triangle * Triangles;
    Vertex * Vertices;
    GLuint VAO, vertexBuffer, normalBuffer;
    GLuint shader;
    void SetupGPU();
};
