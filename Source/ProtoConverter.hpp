#pragma once

#include "core.pb.h"

#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>

////////////////////////////////////////////////////////////////////////////////

namespace ProtoConverter
{
    IronGames::Vector3 Convert(const glm::dvec3& vec);
    IronGames::Matrix3 Convert(const glm::dmat3& mat);
}

////////////////////////////////////////////////////////////////////////////////