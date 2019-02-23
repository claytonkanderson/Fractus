//
//  Window.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "Window.hpp"

int Window::width;
int Window::height;
glm::mat4 Window::P;
glm::mat4 Window::V;
Camera * Window::cam;
GLuint Window::triangleShader;
GLuint Window::regularShader;
glm::vec2 Window::oldMousePos(0.0f, 0.0f);
bool Window::mousePressed0;
bool Window::mousePressed1;

PointLight pointLight;
SpotLight spotLight;
DirLight dirLight;
DirLight dirLight2;

void Window::initialize_objects()
{
    triangleShader = LoadShaders("../../Shaders/TriangleShader.vert", "../../Shaders/TriangleShader.frag");
    regularShader = LoadShaders("../../Shaders/shader.vert", "../../Shaders/shader.frag");
    
    cam = new Camera();

    
    V = glm::lookAt(cam->getPos(), cam->getLookAt(), cam->getUp());
    
    dirLight.on = 1;
    dirLight2.on = 1;
    glUseProgram(regularShader);
    glUniform3f(glGetUniformLocation(regularShader, "viewPos"), cam->getPos().x, cam->getPos().y, cam->getPos().z);
    glUniform1i(glGetUniformLocation(regularShader, "dirLight.on"), dirLight.on);
    glUniform1i(glGetUniformLocation(regularShader, "dirLight2.on"), dirLight2.on);
    
    shaderInit();
}

void Window::clean_up()
{
    glDeleteProgram(triangleShader);
    glDeleteProgram(regularShader);
}

GLFWwindow* Window::create_window(int width, int height)
{
    // Initialize GLFW
    if (!glfwInit())
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return NULL;
    }
    
    // 4x antialiasing
    glfwWindowHint(GLFW_SAMPLES, 4);
    
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    // Create the GLFW window
    GLFWwindow* window = glfwCreateWindow(width, height, "Fracturing Project", NULL, NULL);
    // Check if the window could not be created
    if (!window)
    {
        fprintf(stderr, "Failed to open GLFW window.\n");
        glfwTerminate();
        return NULL;
    }
    
    // Make the context of the window
    glfwMakeContextCurrent(window);
    
    // Set swap interval to 1
    glfwSwapInterval(1);
    // Get the width and height of the framebuffer to properly resize the window
    glfwGetFramebufferSize(window, &width, &height);
    // Call the resize callback to make sure things get drawn immediately
    
    Window::resize_callback(window, width, height);
    
    return window;
}

void Window::resize_callback(GLFWwindow* window, int width, int height)
{
    // Frame buffer is 2x size of window for mac retina graphics
    Window::width = 2*width;
    Window::height = 2*height;
    // Set the viewport size
    glViewport(0, 0, width, height);
    if (height > 0)
    {
        P = glm::perspective(45.0f, (float)width / (float)height, 0.1f, 800.0f);
    }
}

void Window::idle_callback(GLFWwindow* window)
{
    
    
}

void Window::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        if (key == GLFW_KEY_ESCAPE)
        {
            // Close the window. This causes the program to also terminate.
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
        
    }
}

void Window::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (action == 1 && button == 0)
        mousePressed0 = true;
    else if (action == 1 && button == 1)
        mousePressed1 = true;
    else {
        mousePressed0 = false;
        mousePressed1 = false;
    }
}

void Window::scroll_callback(GLFWwindow *window, double xoffset, double yoffset)
{

}

void Window::shaderInit()
{
    pointLight.pos = glm::vec3(0.0f, 0.0f, 1.1f);
    pointLight.constant = 0.0f;
    pointLight.linear = 0.0f;
    pointLight.quadratic = 0.1f;
    
    spotLight.pos = glm::vec3(0,0,5.0f);
    spotLight.dir = glm::vec3(0,0,-1.0f);
    spotLight.constant = 0.0f;
    spotLight.linear = 0.0f;
    spotLight.quadratic = 0.1f;
    spotLight.spot_exponent = 5.f;
    spotLight.cutoff = glm::radians(22.5f);
    
    dirLight.dir = glm::vec3(0.0f,-1.0f,1.0f);
    dirLight.ambient =  2.0f*glm::vec3(.3f,.3f,.3f);
    dirLight.diffuse =  2.0f*glm::vec3(.6f,.6f,.6f);
    dirLight.specular = 0.7f*glm::vec3(.85f,.85f,.85f);
    
    dirLight2.dir = glm::vec3(0.0f,-1.0f,0.0f);
    dirLight2.ambient = glm::vec3(.3f,.3f,.3f);
    dirLight2.diffuse = glm::vec3(.6f,.6f,.6f);
    dirLight2.specular = glm::vec3(.85f,.85f,.85f);
    
    //    glClearColor(0.13f, 0.54f, 0.54f, 0.13f); // Decent like tealish green
//    glClearColor(0.0f, 206.0f/255.f, 209.0f/255.f, 0.5f);
    glm::vec3 pointLightColors[] = {
        10.0f*glm::vec3(0.75f, 0.75f, 0.75f)
    };
    
    
    
    // ***************************************************
    // Directional light 1
    glUseProgram(regularShader);
    glUniform1i(glGetUniformLocation(regularShader, "dirLight.on"), dirLight.on);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight.direction"), dirLight.dir.x,dirLight.dir.y,dirLight.dir.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight.ambient"), dirLight.ambient.x, dirLight.ambient.y, dirLight.ambient.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight.diffuse"), dirLight.diffuse.x,dirLight.diffuse.y,dirLight.diffuse.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight.specular"), dirLight.specular.x,dirLight.specular.y,dirLight.specular.z);
    
    // ***************************************************
    // Directional light 2
    glUniform1i(glGetUniformLocation(regularShader, "dirLight2.on"), dirLight2.on);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight2.direction"), dirLight2.dir.x,dirLight2.dir.y,dirLight2.dir.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight2.ambient"), dirLight2.ambient.x, dirLight2.ambient.y, dirLight2.ambient.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight2.diffuse"), dirLight2.diffuse.x,dirLight2.diffuse.y,dirLight2.diffuse.z);
    glUniform3f(glGetUniformLocation(regularShader, "dirLight2.specular"), dirLight2.specular.x,dirLight2.specular.y,dirLight2.specular.z);
    
    // ***************************************************
    // Point light 1
    glUniform1i(glGetUniformLocation(regularShader, "pointLights[0].on"), 0);
    glUniform3f(glGetUniformLocation(regularShader, "pointLights[0].position"), pointLight.pos.x, pointLight.pos.y, pointLight.pos.z);
    glUniform3f(glGetUniformLocation(regularShader, "pointLights[0].ambient"), 0.2f, 0.2f, 0.2f);
    glUniform3f(glGetUniformLocation(regularShader, "pointLights[0].diffuse"), 0.5f, 0.5f, 0.5f);
    glUniform3f(glGetUniformLocation(regularShader, "pointLights[0].specular"), pointLightColors[0].x,  pointLightColors[0].y,  pointLightColors[0].z);
    glUniform1f(glGetUniformLocation(regularShader, "pointLights[0].constant"), pointLight.constant);
    glUniform1f(glGetUniformLocation(regularShader, "pointLights[0].linear"), pointLight.linear);
    glUniform1f(glGetUniformLocation(regularShader, "pointLights[0].quadratic"), pointLight.quadratic);
    
    // ***************************************************
    // SpotLight
    glUniform1i(glGetUniformLocation(regularShader, "spotLight.on"), 0);
    glUniform3f(glGetUniformLocation(regularShader, "spotLight.position"), spotLight.pos.x, spotLight.pos.y, spotLight.pos.z);
    glUniform3f(glGetUniformLocation(regularShader, "spotLight.direction"), spotLight.dir.x, spotLight.dir.y, spotLight.dir.z);
    glUniform3f(glGetUniformLocation(regularShader, "spotLight.ambient"), 0.2f, 0.2f, 0.2f);
    glUniform3f(glGetUniformLocation(regularShader, "spotLight.diffuse"), 0.5f, 0.5f, 0.5f);
    glUniform3f(glGetUniformLocation(regularShader, "spotLight.specular"), 0.75f, 0.75f, 0.75f);
    glUniform1f(glGetUniformLocation(regularShader, "spotLight.constant"), spotLight.constant);
    glUniform1f(glGetUniformLocation(regularShader, "spotLight.linear"), spotLight.linear);
    glUniform1f(glGetUniformLocation(regularShader, "spotLight.quadratic"), spotLight.quadratic);
    glUniform1f(glGetUniformLocation(regularShader, "spotLight.cutoff"), spotLight.cutoff);
    glUniform1f(glGetUniformLocation(regularShader, "spotLight.spot_exponent"), spotLight.spot_exponent);
}

void Window::cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    glm::vec2 newPos = glm::vec2(xpos, ypos);
    glm::vec2 oldPos = glm::vec2(oldMousePos.x, oldMousePos.y);
    
    glm::vec2 direction = oldPos - newPos;
    
    
    if (mousePressed0)
    {
        glm::vec3 right = glm::normalize(glm::cross(cam->getUp(),cam->getLookAt() - cam->getPos()));
    }
    
    if (mousePressed1)
    {
        cam->rotateView(direction.x, direction.y);
        V = glm::lookAt(cam->getPos(), cam->getLookAt(), cam->getUp());
        glUseProgram(regularShader);
        glUniform3f(glGetUniformLocation(regularShader, "viewPos"), cam->getPos().x, cam->getPos().y, cam->getPos().z);
    }
    
    oldMousePos.x = xpos;
    oldMousePos.y = ypos;
}

void Window::char_callback(GLFWwindow* window, unsigned codepoint)
{

}


void Window::display_callback(GLFWwindow* window)
{
    // Clear the color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
}
