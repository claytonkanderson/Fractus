//
//  demos.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 6/4/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "TetraGroup.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>

using namespace glm;

void init(int demo, float &x_width, float &y_width, float &z_width, float &height, float &elastic, float &poisson, float &toughness,
           ivec3 &div, vec3 &center_pos, vec3 &angular_vel, vec3 &com_vel, mat4 & model);
void init(int demo, TetraGroupInits &inits);

