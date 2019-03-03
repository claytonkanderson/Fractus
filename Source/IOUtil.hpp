//
//  Vertex.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include "TetraGroup.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace IOUtil
{
	void LoadTetrahedronObj(TetraGroup & outGroup, const char * filename);
}

////////////////////////////////////////////////////////////////////////////////
