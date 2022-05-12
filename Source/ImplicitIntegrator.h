#pragma once
#include "core.pb.h"

#include <glm/vec3.hpp>
#include <array>

class TetraGroup;

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
    bool ImplicitUpdate(TetraGroup& group, float timestep, bool saveFrame, IronGames::SimulationFrame* frame, float deltaVThreshold);
}

