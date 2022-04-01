#pragma once

#include <glm/vec3.hpp>
#include <array>

class TetraGroup;

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
    void ImplicitUpdate(TetraGroup& group, float timestep);
}

