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

#include <functional>
#include <unordered_map>

////////////////////////////////////////////////////////////////////////////////

class Window
{
    
public:
	Window(std::shared_ptr<Camera> cam, int width, int height);

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
    void initialize_objects(); // Loads Shaders
    static GLFWwindow* create_window(int width, int height);
    static void resize_callback(GLFWwindow* window, int width, int height);
    static void display_callback(GLFWwindow*);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
    static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    static glm::mat4 getPV() {return P*V;};
    static GLuint * getTriangleShader() {return &triangleShader;}
    static GLuint * getRegularShader() {return &regularShader;}

	void register_key_callback(int key, int action, std::function<void()> callback);

	float SecondsSinceLastUpdate = 0.0;

public:
	void setup_callbacks();
	void setup_camera_controls();

	GLFWwindow* glfwWindow = nullptr;
	std::shared_ptr<Camera> camera;
	// Map from key to callback
	static std::unordered_map<int, std::function<void()>> PressKeyCallbacks;
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