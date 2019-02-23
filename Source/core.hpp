//
//  core.h
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#define GLM_ENABLE_EXPERIMENTAL

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <glm/mat4x4.hpp>

#define SHOWVAR(v) std::cout << #v << ": " << v << std::endl
#define SHOWVEC(v) std::cout << #v << ": " << v.x << " " << v.y << " " << v.z << std::endl
#define SHOWVEC4(v) std::cout << #v << ": " << v.x << " " << v.y << " " << v.z << " " << v.w << std::endl
#define SHOWMAT(m) std::cout << #m << ": " << std::endl  \
<< m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << m[0][3] << std::endl \
<< m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << m[1][3] << std::endl \
<< m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << m[2][3] << std::endl \
<< m[3][0] << " " << m[3][1] << " " << m[3][2] << " " << m[3][3] << std::endl

#define SHOWMAT3(m) std::cout << #m << ": " << std::endl  \
<< m[0][0] << " " << m[0][1] << " " << m[0][2] << std::endl \
<< m[1][0] << " " << m[1][1] << " " << m[1][2] << std::endl \
<< m[2][0] << " " << m[2][1] << " " << m[2][2] << std::endl \

#define SHOWMAT4x3(m) std::cout << #m << ": " << std::endl  \
<< m[0][0] << " " << m[0][1] << " " << m[0][2] << std::endl \
<< m[1][0] << " " << m[1][1] << " " << m[1][2] << std::endl \
<< m[2][0] << " " << m[2][1] << " " << m[2][2] << std::endl \
<< m[3][0] << " " << m[3][1] << " " << m[3][2] << std::endl \

#define SHOWMAT3x4(m) std::cout << #m << ": " << std::endl  \
<< m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << m[0][3] << std::endl \
<< m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << m[1][3] << std::endl \
<< m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << m[2][3] << std::endl \

#define MAX(x,y,z) std::max(std::max(x,y),z)
#define M_PI 3.1415926535897

typedef struct {
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
} Material;

