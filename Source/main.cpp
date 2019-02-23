//
//  main.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "main.hpp"

#include <GL/glew.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

void error_callback(int error, const char* description)
{
    // Print error
    fputs(description, stderr);
}

void setup_callbacks()
{
    // Set the error callback
    glfwSetErrorCallback(error_callback);
    // Set the key callback
    glfwSetKeyCallback(window, Window::key_callback);
    // Set the window resize callback
    glfwSetFramebufferSizeCallback(window, Window::resize_callback);
    glfwSetMouseButtonCallback(window, Window::mouse_button_callback);
    glfwSetCursorPosCallback(window, Window::cursor_position_callback);
    glfwSetCharCallback(window, Window::char_callback);
    //    glfwSetScrollCallback(window, Window::GLFWscrollfun);
    glfwSetScrollCallback(window, Window::scroll_callback);
    //    glfwSetMouseWheelCallback(window, Window::GLFWscrollfun);
}

void setup_glew()
{
    glewExperimental = GL_TRUE;
    // Initialize GLEW
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
        glfwTerminate();
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
}

void setup_opengl_settings()
{
    // Setup GLEW
    setup_glew();
    // Enable depth buffering
    glEnable(GL_DEPTH_TEST);
    // Related to shaders and z value comparisons for the depth buffer
    glDepthFunc(GL_LEQUAL);
    // Set polygon drawing mode to fill front and back of each polygon
    // You can also use the paramter of GL_LINE instead of GL_FILL to see wireframes
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    // Disable backface culling to render both sides of polygons
    glDisable(GL_CULL_FACE);
    // Set clear color
    glClearColor(0.2f, 0.2f, 0.5f, 1.0f);
}

void print_versions()
{
    // Get info of GPU and supported OpenGL version
    printf("Renderer: %s\n", glGetString(GL_RENDERER));
    printf("OpenGL version supported %s\n", glGetString(GL_VERSION));
    
    //If the shading language symbol is defined
#ifdef GL_SHADING_LANGUAGE_VERSION
    std::printf("Supported GLSL version is %s.\n", (char *)glGetString(GL_SHADING_LANGUAGE_VERSION));
#endif
}

int main(int argc, const char * argv[]) {
    
    // Create the GLFW window
    window = Window::create_window((int)(1.5*640), (int)(1.5*480));
    // Print OpenGL and GLSL versions
    print_versions();
    // Setup callbacks
    setup_callbacks();
    // Setup OpenGL settings, including lighting, materials, etc.
    setup_opengl_settings();
    // Initialize objects/pointers for rendering
    
    // Resolves faces clipping ontop of one another
    //    glEnable(GL_CULL_FACE);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);
    
    
    // TIMER1 START
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    high_resolution_clock::time_point t2,t3,t4;
    
    Window::initialize_objects();

    Scene scn;
    
    TetraGroup group;
    
    ///////////////////////////////////
    
    float timestep = 0.0001f;
    int physPerDraw = 15;
    bool pause = false;
    
    int demoNum = 3;
    group.meshView = false;
    TetraGroupInits inits;
    
    init(demoNum, inits);
    
    group.Init(inits); // Creates Grid
    group.SetShader(*Window::getRegularShader());
    scn.AddObject(group);
    
    // Create ground
    PlaneObject * plane = new PlaneObject(100.0f, vec3(0,0,0), vec3(0,1,0));
    plane->SetShader(*Window::getRegularShader());
    scn.AddObject(*plane);
    
    // Want to make sure that each vertex is in each of it's neighboring tetrahedra
    
    while (!glfwWindowShouldClose(window))
    {
        // Main render display callback. Rendering of objects is done here.
        Window::display_callback(window);

        glm::mat4 PV = Window::getPV();
        if (!pause)
        {
            group.UpdateVertices(timestep, physPerDraw);
        }
        
        scn.Draw(PV);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        TetraGroupControls(group, pause);
        cameraControls();
        Window::idle_callback(window);
    }
    
    Window::clean_up();
    // Destroy the window
    glfwDestroyWindow(window);
    // Terminate GLFW
    glfwTerminate();
    
    vector<duration<double>> tetraTimers = group.GetTimers();
    
    t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "Program ran for: " << time_span.count() << " seconds." << endl;
    cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
    cout << "Percent spent in ComputeDeformationForces: " << tetraTimers[0].count() / time_span.count() * 100.0<< endl;
    cout << "Percent spent in ComputeFracture: " << tetraTimers[1].count() / time_span.count() * 100.0<< endl;
    cout << "Percent spent in ComputeSeparation: " << tetraTimers[2].count() / time_span.count() * 100.0<< endl;
    cout << "Percent spent in Vertex Update: " << tetraTimers[3].count() / time_span.count() * 100.0<< endl;
    cout << "Percent spent in Triangle Update: " << tetraTimers[4].count() / time_span.count() * 100.0<< endl;
    cout << "Percent spent in GPU Update: " << tetraTimers[5].count() / time_span.count() * 100.0<< endl;
    cout << "Total time account for: " << (tetraTimers[0].count() + tetraTimers[1].count() + tetraTimers[2].count() + tetraTimers[3].count() +
    tetraTimers[4].count() + tetraTimers[5].count())/(float)(time_span.count()) * 100.0f << endl;
    cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
    
    exit(EXIT_SUCCESS);
    
    return 0;
}

void TetraGroupControls(TetraGroup &group, bool &pause)
{
    int shift = 1;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)) shift = -1;
    else shift = 1;
    
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
        group.Reset();
    
    // Elastic
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        group.inits.elastic *= 1 + shift*0.01f;
        SHOWVAR(group.inits.elastic);
    }
    
    // Poisson Ratio
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        group.inits.poisson += (float)0.002*shift;
        SHOWVAR(group.inits.poisson);
    }
    
    // Toughness
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        group.inits.toughness *= 1 + shift*0.01f;
        SHOWVAR(group.inits.toughness);
    }
    
    // Initial Height
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        group.inits.height += shift;
        group.inits.center_pos.y = group.inits.height;
        group.inits.model = rotate(mat4(1.0f),group.inits.theta,glm::vec3(0,0,1));
        group.inits.model = rotate(mat4(1.0f),group.inits.phi,glm::vec3(0,0,1)) * group.inits.model;
        group.inits.model = translate(mat4(1.0f), group.inits.center_pos) * group.inits.model;
        SHOWVAR(group.inits.height);
    }
    
    // Down velocity
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        group.inits.com_vel = vec3(group.inits.com_vel.x, group.inits.com_vel.y + shift, group.inits.com_vel.z);
        SHOWVEC(group.inits.com_vel);
    }
    
    // Angular Vel about X-axis
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        group.inits.angular_vel = vec3(group.inits.angular_vel.x + 0.5f*shift, group.inits.angular_vel.y, group.inits.angular_vel.z);
        SHOWVEC(group.inits.angular_vel);
    }
    
    // Angular Vel about Z-axis
    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    {
        group.inits.angular_vel = vec3(group.inits.angular_vel.x, group.inits.angular_vel.y, group.inits.angular_vel.z + 0.5f*shift);
        SHOWVEC(group.inits.angular_vel);
    }
    
    // Angle to X-axis
    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    {
        group.inits.theta += 0.1f * shift;
        group.inits.model = rotate(mat4(1.0f),group.inits.theta,glm::vec3(0,0,1));
        group.inits.model = rotate(mat4(1.0f),group.inits.phi,glm::vec3(0,1,0)) * group.inits.model;
        group.inits.model = translate(mat4(1.0f), group.inits.center_pos) * group.inits.model;
        SHOWVAR(group.inits.theta);
    }
    
    // Angle to Z-axis,
    if (glfwGetKey(window, GLFW_KEY_9) == GLFW_PRESS)
    {
        group.inits.phi += 0.1f * shift;
        group.inits.model = rotate(mat4(1.0f),group.inits.theta,glm::vec3(0,0,1));
        group.inits.model = rotate(mat4(1.0f),group.inits.phi,glm::vec3(0,1,0)) * group.inits.model;
        group.inits.model = translate(mat4(1.0f), group.inits.center_pos) * group.inits.model;
        SHOWVAR(group.inits.phi);
    }
    
    // Select Demo
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) { init(0, group.inits); }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) { init(1, group.inits); }
    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) { init(2, group.inits); }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) { init(3, group.inits); }
    
    
    // Toggle Wireframe or Full mode
    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS)
    {
        if (shift < 0) group.meshView = false;
        else group.meshView = true;
    }
    
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
    {
        if (shift < 0) pause = false;
        else pause = true;
    }
    
}

void cameraControls()
{
    static double lastTime = glfwGetTime();
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);
    
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        glm::vec3 right = glm::normalize(glm::cross(glm::normalize((Window::cam->getLookAt() - Window::cam->getPos())), Window::cam->getUp()));
        Window::cam->translate(-1.0f*right * Window::cam->getSpeed() * deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        Window::cam->translate(-1.0f*glm::normalize((Window::cam->getLookAt() - Window::cam->getPos())) * deltaTime * Window::cam->getSpeed());
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        glm::vec3 right = glm::normalize(glm::cross(glm::normalize((Window::cam->getLookAt() - Window::cam->getPos())), Window::cam->getUp()));
        Window::cam->translate(right * Window::cam->getSpeed() * deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        Window::cam->translate(glm::normalize((Window::cam->getLookAt() - Window::cam->getPos())) *deltaTime * Window::cam->getSpeed());
    }
    
    glUniform3f(glGetUniformLocation(Window::regularShader, "viewPos"), Window::cam->getPos().x, Window::cam->getPos().y, Window::cam->getPos().z);
    Window::V = glm::lookAt(Window::cam->getPos(), Window::cam->getLookAt(), Window::cam->getUp());
    
    lastTime = currentTime;
}
