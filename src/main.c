
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
    double mass = 0.05;
    double radiusIn = 0.4;
    double radiusOut = 2.0;
    double dt = 0.015;
    double tEnd = 5.0;
    double tol = 1e-6;
    double t = 0;
    double iterMax = 100;
    double gamma = 0.6;

    femPoissonProblem* theProblemu = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");
    femPoissonProblem* theProblemv = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");
    
    femGrains* theGrains = femGrainsCreateSimple(n, radius, mass, radiusIn, radiusOut, gamma);
    femPoissonSolve(theProblemu, theProblemv, theGrains);
    
    
    double* B = malloc(sizeof(double)*theProblemu->mesh->nNode);
    double* B1;
    double* B2;
    int option=1;
    GLFWwindow* window = glfemInit("Projet meca");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;
    
    do {
        int i,w,h;
        double currentTime = glfwGetTime();
        
       
               glfwGetFramebufferSize(window,&w,&h);
            glfemReshapeWindows(radiusOut,w,h);
        
        if(option==1){
            
        for (i = 0; i < theProblemu->mesh->nNode; i++) {
            B1 = theProblemu->system->B;
            B2 = theProblemv->system->B;
            B[i] = sqrt( (B1[i]*B1[i]) + (B2[i]*B2[i]) );
        }
             glfemPlotField(theProblemu->mesh, B);
        }
      //glfwGetFramebufferSize(window,&w,&h);
       // glfemReshapeWindows(radiusOut,w,h);
        //glfemPlotField(theProblemu->mesh, theProblemu->system->B);
        //glfemPlotField(theProblemu->mesh, theProblemv->system->B);
        //glfemPlotField(theProblemu->mesh, norme);
        for (i=0 ;i < theGrains->n; i++) {
            if(option==2){
                if(sqrt(theGrains->vx[i]*theGrains->vx[i]+theGrains->vy[i]*theGrains->vy[i])<1.3){
                     glColor3f(1.0, 1.0, 1.0);
                glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);
            }
            else if (sqrt(theGrains->vx[i]*theGrains->vx[i]+theGrains->vy[i]*theGrains->vy[i])<1.6){
                glColor3f(255, 165, 0.0);
                glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);
            }
            else if (sqrt(theGrains->vx[i]*theGrains->vx[i]+theGrains->vy[i]*theGrains->vy[i])<1.9){
                glColor3f(1.0, 0.0, 0.25f);
                glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);}
        else {
            glColor3f(0.0, 0.0, 0.0);
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);
        }
            }
            else {
                glColor3f(0.0, 0.0, 0.0);
                glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);
            }
        }
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
            femPoissonSolve(theProblemu, theProblemv, theGrains);

            t += dt; }
        if(glfwGetKey(window,'V')==GLFW_PRESS) option=1;
        if(glfwGetKey(window,'S')==GLFW_PRESS) option=0;
        if(glfwGetKey(window,'D')==GLFW_PRESS) option=2;
        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
            if (glfwGetKey(window,'R') == GLFW_PRESS)
                theRunningMode = 1;
            if (glfwGetKey(window,'S') == GLFW_PRESS)
                theRunningMode = 0; }
        
    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) != 1);
    
    
    free(B);
    glfwTerminate();
    femGrainsFree(theGrains);
    femPoissonFree(theProblemu);
    femPoissonFree(theProblemv);
    exit(EXIT_SUCCESS);
    
}

