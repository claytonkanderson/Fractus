//
//  Object.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include <glm/mat4x4.hpp>

class Object
{
public:
    virtual ~Object() {};
    virtual void Draw(glm::mat4 &MVP) = 0;
};
