
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0
#define TRUE  1

typedef enum { FEM_TRIANGLE, FEM_QUAD } femElementType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    double *area;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    void(*x2)(double *xsi, double *eta);
    void(*phi2)(double xsi, double eta, double *phi);
    void(*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femFullSystem *system;
    double* area;
} femPoissonProblem;


femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh             *femMeshRead(const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshClean(femMesh *theMesh);
void                 femMeshFree(femMesh *theMesh);
void                 femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemAlloc(femFullSystem* mySystem, int size);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);


typedef struct {
    int n;
    double radiusIn;
    double radiusOut;
    double gravity[2];
    double gamma;
    double *x;
    double *y;
    int *elem; /*tableau d'appartenance a un element*/
    double *vx;
    double *vy;
    double *r;
    double *m;
    double *dvBoundary;
    double *dvContacts;
} femGrains;


femGrains  *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut, double gamma);
femGrains  *femGrainsCreateTiny(double radiusIn, double radiusOut);
void        femGrainsFree(femGrains *myGrains);
void        femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *theProblem,femPoissonProblem* theProblemV);
double      femGrainsContactIterate(femGrains *myGrains, double dt, int iter);


femPoissonProblem   *femPoissonCreate(const char *filename);
void                 femPoissonFree(femPoissonProblem *theProblem);
void femPoissonSolve(femPoissonProblem *theProblemU,femPoissonProblem *theProblemV, femGrains *myGrains);


double      femMin(double *x, int n);
double      femMax(double *x, int n);
void        femError(char *text, int line, char *file);
void        femErrorScan(int test, int line, char *file);
void        femWarning(char *text, int line, char *file);

/* ADDED FUNCTIONS */

int WithinTriangle(double *x, double *y, double xc, double  yc, double Aire);
void ComputeArea (femPoissonProblem *theProblem);
void invefctforme(femPoissonProblem *theProblem, double xx, double yy, double *x, double *y, double *tau);

#endif
