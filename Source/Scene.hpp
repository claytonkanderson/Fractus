//
//  Scene.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "Object.hpp"
#include "Light.hpp"
#include "core.hpp"
#include "Camera.hpp"

class Scene {
public:
    Scene()										{SkyColor.Set(0.2f,0.2f,0.5f);}
    void AddObject(Object &obj)					{Objects.push_back(&obj);}
    void AddLight(Light &lgt)					{Lights.push_back(&lgt);}
    void SetSkyColor(const Color sky)			{SkyColor=sky;}
    
    void Draw(glm::mat4 &PV){
        for (unsigned int i = 0 ; i < Objects.size(); i++)
            Objects[i]->Draw(PV);
    }
    
    Object &GetObj(int i)						{return *Objects[i];}
    int GetNumLights()							{return (int)Lights.size();}
    Light &GetLight(int i)						{return *Lights[i];}
    Color GetSkyColor()							{return SkyColor;}

private:
    std::vector<Object*> Objects;
    std::vector<Light*> Lights;
    Color SkyColor;

};
