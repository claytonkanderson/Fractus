//
//  main.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/16/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "core.hpp"
#include "Scene.hpp"
#include "Window.hpp"
#include "TetraGroup.hpp"
#include "PlaneObject.hpp"
#include "IOUtil.hpp"
#include "demos.hpp"

#include <GLFW/glfw3.h>
#include <GL/glew.h>
#include <chrono>

GLFWwindow* window;
void cameraControls(Window & window);
void TetraGroupControls(TetraGroup &group, GLFWwindow* window, bool pause);

using namespace std;
using namespace std::chrono;

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
    
	auto camera = std::make_shared<Camera>();

    // Create the GLFW window
    auto window = Window(camera, (int)(1.5*640), (int)(1.5*480));
    
    // TIMER1 START
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    high_resolution_clock::time_point t2,t3,t4;
    
    Scene scene;
    TetraGroup group;
    
    float timestep = 0.01667 / 2000.0;
    int physPerDraw = 25;
    bool pause = false;
    
    int demoNum = 5;
    group.meshView = false;
    TetraGroupInits inits;
    
    //init(demoNum, inits);
    
    //group.Init(inits); // Creates Grid

	// For next time : add material coordinate property to vertices,
	// verify usage of material / world positions.
	// Verify locations of UpdateBeta / UpdateMass (as well as their impls)

	IOUtil::LoadTetrahedronObj(group, "BowlTetrahedra.obj");

	for (auto & tet : group.GetTetrahedra())
		tet->SetConstants(0, 5.29e7, 0, 198);
	group.Toughness = 106;

	group.SetInitialConditions(
		glm::translate(glm::mat4(), glm::vec3(0, 1.5, 0)),
		glm::vec3(0, 0, 0), 
		glm::vec3()
	);

    group.SetShader(*Window::getRegularShader());
	scene.AddObject(group);
    
    // Create ground
    PlaneObject * plane = new PlaneObject(100.0f, vec3(0,0,0), vec3(0,1,0));
    plane->SetShader(*Window::getRegularShader());
	scene.AddObject(*plane);
    
    // Want to make sure that each vertex is in each of it's neighboring tetrahedra
    
    while (!glfwWindowShouldClose(window.glfwWindow))
    {
        // Main render display callback. Rendering of objects is done here.
        Window::display_callback(window.glfwWindow);

        glm::mat4 PV = Window::getPV();
        if (!pause)
            group.UpdateVertices(timestep, physPerDraw);
        
		scene.Draw(PV);
        
        glfwSwapBuffers(window.glfwWindow);
        glfwPollEvents();
        
        TetraGroupControls(group, window.glfwWindow, pause);
        cameraControls(window);
    }
    
    Window::clean_up();
    // Destroy the window
    glfwDestroyWindow(window.glfwWindow);
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

////////////////////////////////////////////////////////////////////////////////

void TetraGroupControls(TetraGroup &group, GLFWwindow* window, bool pause)
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

////////////////////////////////////////////////////////////////////////////////

void cameraControls(Window & window)
{
    static double lastTime = glfwGetTime();
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);
    
	window.SecondsSinceLastUpdate = deltaTime;
    
	auto & camera = window.camera;

    glUniform3f(glGetUniformLocation(Window::regularShader, "viewPos"), camera->getPos().x, camera->getPos().y, camera->getPos().z);
    Window::V = glm::lookAt(camera->getPos(), camera->getLookAt(), camera->getUp());
    
    lastTime = currentTime;
}

////////////////////////////////////////////////////////////////////////////////