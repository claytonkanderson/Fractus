//
//  Camera.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include <glm/vec3.hpp>

class Camera
{
public:
    Camera() {
        radius = 25.0f;
        theta = 0;
		phi = (float)(0.5*M_PI - 0.1);
        center = glm::vec3(0,0,0);
        
        updatePosition();
        
        look_at = center;                           // d  | This is where the camera looks at
        up = glm::vec3 (0.0f, 1.0f, 0.0f);			// up | What orientation "up" is
        std::cout << "Camera starting at ";
        SHOWVEC(pos);
        std::cout << "Camera looking at ";
        SHOWVEC(look_at);
    }
    glm::vec3 getPos() {return pos;}
    glm::vec3 getLookAt() {return look_at;}
    glm::vec3 getUp() {return up;}
    float getSpeed() {return this->speed;}
    void translate(glm::vec3 vec);
    void setAngles(float theta, float phi);
    void rotateView(float xDir, float yDir);
    void rotate(glm::vec3 rotAxis, float deg);
    
private:
    void updatePosition();
    float speed = 20.0f;
    float theta;
    float phi;
    float radius;
    glm::vec3 center;
    glm::vec3 pos;
    glm::vec3 look_at;
    glm::vec3 up;
};

