//
//  Camera.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Camera.hpp"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

void Camera::translate(glm::vec3 vec)
{
    center = glm::vec3(glm::translate(glm::mat4(1.0f), vec) * glm::vec4(center,1.0f));
    updatePosition();
    look_at = glm::vec3(glm::translate(glm::mat4(1.0f), vec) * glm::vec4(look_at,1.0f));
}


// Once Center has been updated, call this to update x, y, z from theta & phi
// theta measured from the z-axis to point (rotation about y-axis)
// phi measured from y-axis to point (rotation about z-axis)
void Camera::updatePosition()
{
    float x = radius * sinf(phi) * sinf(theta) + center.x;
    float y = radius * cosf(phi) + center.y;
    float z = radius * sinf(phi) * cosf(theta) + center.z;
    
    pos = glm::vec3 (x, y, z);
}


void Camera::rotate(glm::vec3 rotAxis, float deg)
{
    if (rotAxis != glm::vec3(0,0,0))
    {
        pos = glm::vec3((glm::rotate(glm::mat4(1.0f), deg / 180.0f * glm::pi<float>(), rotAxis))*glm::vec4(pos,1.0f));
        look_at = glm::vec3((glm::rotate(glm::mat4(1.0f), deg / 180.0f * glm::pi<float>(), rotAxis))*glm::vec4(look_at,0));
        up = glm::vec3((glm::rotate(glm::mat4(1.0f), deg / 180.0f * glm::pi<float>(), rotAxis))*glm::vec4(up,0));
    }
    
}

void Camera::rotateView(float xDir, float yDir)
{
    theta += 0.005f * xDir;
    theta = fmodf(theta, 2*M_PI);
    phi += 0.005f * yDir;
    
    const float phi1 = phi;
    const float min1 = 0.0001f;
    const float max1 = M_PI;
    phi = glm::clamp(phi1, min1, max1);
    
    updatePosition();
}

void Camera::setAngles(float theta, float phi)
{
    this->theta = fmodf(theta, 2*M_PI);
    //    this->theta = theta;
    
    const float phi1 = phi;
    const float min1 = 0.0001f;
    const float max1 = M_PI;
    this->phi = glm::clamp(phi1, min1, max1);
    updatePosition();
}
