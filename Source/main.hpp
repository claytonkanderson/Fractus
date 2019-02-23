//
//  main.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "Scene.hpp"
#include "Window.hpp"
#include "TetraGroup.hpp"
#include "PlaneObject.hpp"
#include "demos.hpp"

#include <GLFW/glfw3.h>

GLFWwindow* window;
void cameraControls();
void TetraGroupControls(TetraGroup &group, bool &pause);
