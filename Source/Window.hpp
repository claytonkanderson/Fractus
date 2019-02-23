//
//  Window.hpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#pragma once

#include "core.hpp"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Camera.hpp"
#include "Shader.hpp"

class Window
{
    
public:
    static Camera * cam;
    static glm::mat4 P;
    static glm::mat4 V;
    static int width;
    static int height;
    static GLuint triangleShader;
    static GLuint regularShader;
    static glm::vec2 oldMousePos;
    static bool mousePressed0;
    static bool mousePressed1;
    
    static void shaderInit();
    static void clean_up();
    static void initialize_objects(); // Loads Shaders
    static GLFWwindow* create_window(int width, int height);
    static void resize_callback(GLFWwindow* window, int width, int height);
    static void idle_callback(GLFWwindow*);
    static void display_callback(GLFWwindow*);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
    static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    static void char_callback(GLFWwindow* window, unsigned codepoint);
    static glm::mat4 getPV() {return P*V;};
    static GLuint * getTriangleShader() {return &triangleShader;}
    static GLuint * getRegularShader() {return &regularShader;}
};

struct DirLight {
    int on;
    glm::vec3 dir;
    
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};

struct PointLight {
    int on;
    glm::vec3 pos;
    
    float constant;
    float linear;
    float quadratic;
    
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};

struct SpotLight {
    int on;
    glm::vec3 pos;
    glm::vec3 dir;
    //    float cutOff;
    //    float outerCutOff;
    float cutoff;
    float spot_exponent;
    
    float constant;
    float linear;
    float quadratic;
    
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};