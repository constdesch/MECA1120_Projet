
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
    
    int    n         = 15;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;
    double dt        = 1e-1;
    double tEnd      = 8.0;
    double tol       = 1e-6;
    double t         = 0;
    double iterMax   = 100;
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);
    


	femPoissonProblem* theProblemU = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");
    femPoissonProblem* theProblemV = femPoissonCreate("../data/meca1120-projet-meshMedium.txt");

	femPoissonSolve(theProblemU,0, theGrains);
    femPoissonSolve(theProblemV,1, theGrains);
    
    int i;
    double *B = malloc(sizeof(double) * theProblemU->system->size);
    double *B1 = theProblemU->system->B;
    double *B2 = theProblemV->system->B;
    for (i=0;i<theProblemU->system->size;i++){
        B[i] = sqrt(B1[i]*B1[i] + B2[i]*B2[i]);
    }
    
   
  //  A decommenter pour obtenir l'exemple de la seance d'exercice :-)
  //  femGrains* theGrains = femGrainsCreateTiny(radiusIn,radiusOut);;

    GLFWwindow* window = glfemInit("MECA1120 : Projet final");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1;
    float theVelocityFactor = 0.25;
     
    do {
        int i,w,h;
        double currentTime = glfwGetTime();
        
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(radiusOut,w,h);   


		glfemPlotField(theProblemU->mesh,B);

        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(1,0,0.85f); glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]);
			glColor3f(0, 0, 0); glfemDrawCircle(theGrains->x[i], theGrains->y[i], theGrains->r[i]);
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
  //
  // A decommenter pour pouvoir progresser pas par pas
  //          printf("press CR to compute the next time step >>");
  //          char c= getchar();
  //
            femGrainsUpdate(theGrains,dt,tol,iterMax, theProblemU);
            femPoissonSolve(theProblemU,0, theGrains);
            femPoissonSolve(theProblemV,1, theGrains);
            int j;
            B1 = theProblemU->system->B;
            B2 = theProblemV->system->B;
            for (j=0;j<theProblemU->system->size;j++){
                B[j] = sqrt(B1[j]*B1[j] + B2[j]*B2[j]);
            }
            glfemPlotField(theProblemU->mesh,B);
            t += dt; }
         
        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS) 
        	    theRunningMode = 1; 
          if (glfwGetKey(window,'S') == GLFW_PRESS) 
        	    theRunningMode = 0; }
            
    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	        (!glfwWindowShouldClose(window)));
	   
               
    glfwTerminate();
    free(B);
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
}



