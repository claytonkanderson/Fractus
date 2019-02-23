//
//  demos.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 6/4/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "demos.hpp"

void init(int demo, TetraGroupInits &inits)
{
    switch (demo) {
            // Half broken fry
        case 0:
            inits.x_width = 40.0f;
            inits.y_width = 1.0f;
            inits.z_width = 1.0f;
            inits.height = 18.0f;
            inits.elastic = 1.00e+09;
            inits.poisson = 0.47f;
            inits.toughness = 2.17e+07;
            inits.div = ivec3(32,2,2);
            inits.center_pos = vec3(0,inits.height,0);
            inits.angular_vel = vec3(0, 1, 4);
            inits.com_vel = vec3(0, -70, 0);
            inits.theta = -0.5f; // Double check how the angles are done on the rest of them
            inits.phi = 0.0f;
            inits.model = rotate(mat4(1.0f),inits.theta,glm::vec3(0,0,1));
            inits.model = rotate(mat4(1.0f),inits.phi,glm::vec3(0,1,0)) * inits.model;
            inits.model = translate(mat4(1.0f), inits.center_pos) * inits.model;
            break;
        
            // Plate broken by rotation
        case 1:
            inits.x_width = 20.0f;
            inits.y_width = 1.0f;
            inits.z_width = 20.0f;
            inits.height = 10.0f;
            inits.elastic = 2e+08;
            inits.poisson = 0.45f;
            inits.toughness = 2e+7;
            inits.div = ivec3(8,4,8);
            inits.center_pos = vec3(0,inits.height,0);
            inits.angular_vel = vec3(0, 15.45, 0);
            inits.com_vel = vec3(0, 0, 0);
            inits.theta = 0.0f;
            inits.phi = 0.0f;
            inits.model = rotate(mat4(1.0f),inits.theta,glm::vec3(0,0,1));
            inits.model = rotate(mat4(1.0f),inits.phi,glm::vec3(0,1,0)) * inits.model;
            inits.model = translate(mat4(1.0f), inits.center_pos) * inits.model;
            break;
            
            // Cube semi-shattered
        case 2:
            inits.x_width = 10.0f;
            inits.y_width = 10.0f;
            inits.z_width = 10.0f;
            inits.height = 10.0f;
            inits.elastic = 4.52e+8;
            inits.poisson = 0.45f;
            inits.toughness = 1.01e+8;
            inits.div = ivec3(4,4,4);
            inits.center_pos = vec3(0,inits.height,0);
            inits.angular_vel = vec3(0, 1, 0);
            inits.com_vel = vec3(0, -45, 0);
            inits.model = rotate(mat4(1.0f),inits.theta,glm::vec3(0,0,1));
            inits.model = rotate(mat4(1.0f),inits.phi,glm::vec3(0,1,0)) * inits.model;
            inits.model = translate(mat4(1.0f), inits.center_pos) * inits.model;
            break;
            
            // Single layer slab
        case 3:
            inits.x_width = 10.0f;
            inits.y_width = 1.0f;
            inits.z_width = 30.0f;
            inits.height = 10.0f;
            inits.elastic = 5.24e+08;
            inits.poisson = 0.45f;
            inits.toughness = 7.1e+7;
            inits.div = ivec3(8,2,8);
            inits.center_pos = vec3(0,inits.height,0);
            inits.angular_vel = vec3(5, 0, 0);
            inits.com_vel = vec3(0, -75, 0);
            inits.theta = 0.0f;
            inits.phi = 1.0f;
            inits.model = rotate(mat4(1.0f),inits.theta,glm::vec3(0,0,1));
            inits.model = rotate(mat4(1.0f),inits.phi,glm::vec3(0,1,0)) * inits.model;
            inits.model = translate(mat4(1.0f), inits.center_pos) * inits.model;
            break;
            
        case 4:
            inits.x_width = 4.0f;
            inits.y_width = 2.0f;
            inits.z_width = 2.0f;
            inits.height = 7.0f;
            inits.elastic = 2e+07;
            inits.poisson = 0.45f;
            inits.toughness = 2.4e+6;
            inits.div = ivec3(2,2,1);
            inits.center_pos = vec3(0,inits.height,0);
            inits.angular_vel = vec3(0, 0, 0);
            inits.com_vel = vec3(0, -25, 0);
            inits.theta = 0.0f;
            inits.phi = 0.0f;
            inits.model = rotate(mat4(1.0f),inits.theta,glm::vec3(0,0,1));
            inits.model = rotate(mat4(1.0f),inits.phi,glm::vec3(0,1,0)) * inits.model;
            inits.model = translate(mat4(1.0f), inits.center_pos) * inits.model;
            break;
            
        default:
            break;
    }
}
