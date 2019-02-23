#version 330 core

uniform mat4 MVP;
uniform mat4 model;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 vertNormal;
//layout (location = 2) in vec2 texcoord;

out vec3 FragPos;
out vec3 Normal;
//out vec2 Texcoord;

void main() {
    // Pass some variables to the fragment shader
    
    gl_Position = MVP * vec4(position, 1.0f);
    
    Normal = mat3(transpose(inverse(model))) * vertNormal; // Tutorial Code
    FragPos = vec3(model * vec4(position, 1.0f)); // Tutorial Code
//    Texcoord = texcoord;
    
}
