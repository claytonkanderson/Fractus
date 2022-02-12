#include "ProtoConverter.hpp"

////////////////////////////////////////////////////////////////////////////////

IronGames::Vector3 ProtoConverter::Convert(const glm::dvec3& vec)
{
    IronGames::Vector3 v;
    v.set_x(vec.x);
    v.set_y(vec.y);
    v.set_z(vec.z);
    return v;
}

////////////////////////////////////////////////////////////////////////////////

IronGames::Matrix3 ProtoConverter::Convert(const glm::dmat3& mat)
{
    IronGames::Matrix3 m;
    *m.mutable_col0() = Convert(mat[0]);
    *m.mutable_col1() = Convert(mat[1]);
    *m.mutable_col2() = Convert(mat[2]);
    return m;
}

////////////////////////////////////////////////////////////////////////////////