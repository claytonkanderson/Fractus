#version 330 core

out vec4 color;
in vec3 Normal;
in vec3 FragPos;

void main()
{
//        color = vec4(1.0,0.0,0.0, 1.0);
    
        color = vec4(abs(Normal), 1.0);
//    color = vec4(Normal, 1.0);
}
