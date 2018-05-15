
/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *  Homework 4 for 17-18 : Discrete Grains
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */


#include "glfem.h"


int main(void)
{
    int    n = 50;
    double radius = 0.05;
    double mass = 0.1;
    double radiusIn = 0.4;
    double radiusOut = 2.0;
    double vext = 5;
    double dt = 0.05;
    double tEnd = 5.0;
    double tol = 1e-6;
    double t = 0;
    double iterMax = 100;
    double gamma = 0.6;
    double mu = 1e-3; //1.8e-3;
    
    int i;
    
    femPoissonProblem* theProblemu = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");
    femPoissonProblem* theProblemv = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");
    
    femGrains* theGrains = femGrainsCreateSimple(n, radius, mass, radiusIn, radiusOut, gamma);
    femPoissonSolveu(theProblemu, radiusIn, radiusOut, vext, mu, theGrains);
    femPoissonSolvev(theProblemv, radiusIn, radiusOut, vext, mu, theGrains);
    
    
    double* norme = malloc(sizeof(double)*theProblemu->mesh->nNode);
    
    
    GLFWwindow* window = glfemInit("Simulation multi-echelle de milieux granulaires immerges.");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;
    
    do {
        int i,w,h;
        double currentTime = glfwGetTime();
        
        
        
        for (i = 0; i < theProblemu->mesh->nNode; i++) {
            norme[i] = pow(pow(theProblemu->system->B[i], 2) + pow(theProblemv->system->B[i], 2), 0.5);
        }
        
        
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(radiusOut,w,h);
        //glfemPlotField(theProblemu->mesh, theProblemu->system->B);
        //glfemPlotField(theProblemu->mesh, theProblemv->system->B);
        glfemPlotField(theProblemu->mesh, norme);
        for (i=0 ;i < theGrains->n; i++) {
            glColor3f(1.0, 1.0, 1.0);
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);
        char theMessage[256];
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);
            
            femGrainsUpdate(theGrains, dt, tol, iterMax, theProblemu, theProblemv);
            femFullSystemInit(theProblemu->system);
            femFullSystemInit(theProblemv->system);
            femPoissonSolveu(theProblemu, radiusIn, radiusOut, vext, mu, theGrains);
            femPoissonSolvev(theProblemv, radiusIn, radiusOut, vext, mu, theGrains);
            t += dt; }
        
        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
            if (glfwGetKey(window,'R') == GLFW_PRESS)
                theRunningMode = 1;
            if (glfwGetKey(window,'S') == GLFW_PRESS)
                theRunningMode = 0; }
        
    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) != 1);
    
    
    free(norme);
    glfwTerminate();
    femGrainsFree(theGrains);
    femPoissonFree(theProblemu);
    femPoissonFree(theProblemv);
    exit(EXIT_SUCCESS);
    
}

