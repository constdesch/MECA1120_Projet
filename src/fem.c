
/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */


#include "fem.h"

/* ===================== ADDED FUNCTIONS ==================== */

/* @area: calcule l'aire du triangle dont les sommets sont A(xa,ya), B(xb,yb) et C(xc,yc)
 * @pre: /
 * @post: retourne l'aire du triangle
 */
double area (double xa, double xb, double xc, double ya, double yb, double yc){
    return 1 /2 * fabs( (xb - xa)*(yc - ya) - (xc - xa)*(yb - ya) );
}

/*
 * @ComputeArea: Rempli le tableau theProblem->area qui contient les aires de chaque triangle du maillage
 *               suivant la même numérotation que le tableau theProblem->elem
 * @pre: theProblem->mesh != NULL et theProblem->mesh->nLocalNode == 3
 * @post: /
 */
void ComputeArea (femPoissonProblem *theProblem){
    femMesh *theMesh = theProblem->mesh;
    int elem;
    int map[3];
    double x[3], y[3];
    for (elem=0; elem<theMesh->nElem;elem++){
        femMeshLocal(theMesh, elem, map, x, y);
        double xa = x[0], xb = x[1], xc = x[2];
        double ya = y[0], yb = y[1], yc = y[2];
        theProblem->area[elem] = 1 / 2 * fabs((xb - xa)*(yc - ya) - (xc - xa)*(yb - ya));
    }
}
/*
 * @WithinTriangle: vérifie si un point P(xc,yc) appartient à un triangle donc les abscices des sommets sont
 *                  contenues dans *x et les ordonnées dans *y suivant la même numérotation. Aire est l'aire
 *                  du triangle.
 * @pre: x et y sont de taille 3, Aire > 0
 * @post: retourne 1 si P appartient au triangle, 0 sinon.
 */
int WithinTriangle(double *x, double *y, double xc, double yc, double Aire) {
    double aire1 = 1 / 2 * area(xc, x[1],x[2],yc,y[1],y[2]);
    double aire2 = 1 / 2 * area(x[0],xc,x[2],y[0],yc,y[2]);
    double aire3 = 1 / 2 * area(x[0],x[1],xc,y[0],y[1],yc);
    double s = aire1 + aire2 + aire3;
    return (Aire == s);
}
/*
 * @invefctforme: donne l'image du point P(xx,yy) évalué sur le triangle parent par les fonctions de forme. Il y a 3
 *                fonctions de forme, tau est de taille 3 et contient l'image de chacune de ces fonctions.
 * @pre: *tau, *x et *y de taille 3, theProblem != NULL,
 * @post: /
 *
 */
void invefctforme(femPoissonProblem *theProblem, double xx, double yy, double *x, double *y, double *tau) {
    double *tab = malloc(sizeof(double) * 2);
    tab[0] = ((xx - x[0])*(y[2] - y[0]) - (yy - y[0])*(x[2] - x[0])) / ((x[1] - x[0])*(y[2] - y[0]) - (y[1] - y[0])*(x[2] - x[0]));
    tab[1] = ((xx - x[0])*(y[0] - y[1]) - (yy - y[0])*(x[0] - x[1])) / ((x[1] - x[0])*(y[2] - y[0]) - (y[1] - y[0])*(x[2] - x[0]));
    femDiscretePhi2(theProblem->space, tab[0], tab[1], tau);
    free(tab);
}

/* ======================== GRAINS & POISSON ======================= */

static const double _gaussQuad4Xsi[4] = { -0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Eta[4] = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000 };
static const double _gaussTri3Xsi[3] = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3] = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3] = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n = 4;
        theRule->xsi = _gaussQuad4Xsi;
        theRule->eta = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight;
    }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n = 3;
        theRule->xsi = _gaussTri3Xsi;
        theRule->eta = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight;
    }
    else Error("Cannot create such an integration rule !");
    return theRule;
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] = 1.0;  eta[0] = 1.0;
    xsi[1] = -1.0;  eta[1] = 1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] = 1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = (1.0 + eta) / 4.0;
    dphidxsi[1] = -(1.0 + eta) / 4.0;
    dphidxsi[2] = -(1.0 - eta) / 4.0;
    dphidxsi[3] = (1.0 - eta) / 4.0;
    dphideta[0] = (1.0 + xsi) / 4.0;
    dphideta[1] = (1.0 - xsi) / 4.0;
    dphideta[2] = -(1.0 - xsi) / 4.0;
    dphideta[3] = -(1.0 + xsi) / 4.0;
    
}

void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0;  eta[0] = 0.0;
    xsi[1] = 1.0;  eta[1] = 0.0;
    xsi[2] = 0.0;  eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] = 1.0;
    dphidxsi[2] = 0.0;
    dphideta[0] = -1.0;
    dphideta[1] = 0.0;
    dphideta[2] = 1.0;
    
}
double* femDiscreteinvetaxsi(double xc, double yc, double x[4], double y[4]) {
    double *tab = malloc(sizeof(double) * 2);
    //((x - x1)*(y2 - y1) - (y - y1)*(x2 - x1)) / ((x3 - x1)*(y2 - y1) - (y3 - y1)*(x2 - x1));
    tab[1] = ((xc - x[0])*(y[0] - y[1]) -(yc - y[0])*(x[1] - x[0])) / ((-x[0] + x[2])*(-y[0] + y[1]) - (y[2] - y[0])*(x[1] - x[0]));
    //theProblem->mesh->xsi = (x - x1 - (((x - x1)*(y2 - y1) - (y - y1)*(x2 - x1)) / ((x3 - x1)*(y2 - y1) - (y3 - y1)*(x2 - x1)))*(x3 - x1)) / (x2 - x1);
    tab[0] = ((xc - x[0]) - ((xc - x[0])*(y[0] - y[1]) - (yc - y[0])*(x[1] - x[0])) / ((-x[0] + x[2])*(-y[0] + y[1]) - (y[2] - y[0])*(x[1] - x[0])))*(x[2] - x[0]) / (x[1] - x[0]);
    return tab;
}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n = 4;
        theSpace->x2 = _q1c0_x;
        theSpace->phi2 = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n = 3;
        theSpace->x2 = _p1c0_x;
        theSpace->phi2 = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
    }
    else Error("Cannot create such a discrete space !");
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi, eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi, eta, phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i, j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace, xsi, eta);
    for (i = 0; i < n; i++) {
        
        femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
        femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
        
        for (j = 0; j < n; j++) {
            printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
            printf(" phi(%d)=%+.1f", j, phi[j]);
            printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
            printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
        }
        printf(" \n");
    }
}



femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));
    
    int i, trash, *elem;
    
    FILE* file = fopen(filename, "r");
    if (file == NULL) Error("No mesh file !");
    
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theMesh->X[i], &theMesh->Y[i]));
    }
    
    char str[256]; if (fgets(str, sizeof(str), file) == NULL) Error("Corrupted mesh file !");
    
    if (!strncmp(str, "Number of triangles", 19)) {
        ErrorScan(sscanf(str, "Number of triangles %d \n", &theMesh->nElem));
        theMesh->elem = malloc(sizeof(int) * 3 * theMesh->nElem);
        theMesh->nLocalNode = 3;
        theMesh->area = malloc(sizeof(double)*theMesh->nElem);
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i * 3]);
            ErrorScan(fscanf(file, "%d : %d %d %d\n", &trash, &elem[0], &elem[1], &elem[2]));
        }
    }
    else if (!strncmp(str, "Number of quads", 15)) {
        printf("%s \n", str);
        ErrorScan(sscanf(str, "Number of quads %d \n", &theMesh->nElem));
        theMesh->elem = malloc(sizeof(int) * 4 * theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i * 4]);
            ErrorScan(fscanf(file, "%d : %d %d %d %d\n", &trash, &elem[0], &elem[1], &elem[2], &elem[3]));
        }
    }
    
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->area);
    free(theMesh->elem);
    free(theMesh);
}

void femMeshClean(femMesh *theMesh)
{
    int i, j, *elem;
    
    
    int *check = malloc(sizeof(int)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i)
        check[i] = 0;
    int *map = malloc(sizeof(int)*theMesh->nNode);
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        for (j = 0; j < theMesh->nLocalNode; ++j) {
            check[elem[j]] = 1;
        }
    }
    int iGlo = 0;
    for (i = 0; i < theMesh->nNode; ++i) {
        if (check[i] != 0) {
            map[i] = iGlo;
            theMesh->X[iGlo] = theMesh->X[i];
            theMesh->Y[iGlo] = theMesh->Y[i];
            iGlo++;
        }
    }
    theMesh->nNode = iGlo;
    
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        for (j = 0; j < theMesh->nLocalNode; ++j) {
            elem[j] = map[elem[j]];
        }
    }
    
    free(check);
    free(map);
}


void femMeshWrite(const femMesh *theMesh, const char *filename)
{
    int i, *elem;
    
    FILE* file = fopen(filename, "w");
    
    
    
    
    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file, "%6d : %14.7e %14.7e \n", i, theMesh->X[i], theMesh->Y[i]);
    }
    
    if (theMesh->nLocalNode == 4) {
        fprintf(file, "Number of quads %d \n", theMesh->nElem);
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i * 4]);
            fprintf(file, "%6d : %6d %6d %6d %6d \n", i, elem[0], elem[1], elem[2], elem[3]);
        }
    }
    else if (theMesh->nLocalNode == 3) {
        fprintf(file, "Number of triangles %d \n", theMesh->nElem);
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i * 3]);
            fprintf(file, "%6d : %6d %6d %6d \n", i, elem[0], elem[1], elem[2]);
        }
    }
    
    fclose(file);
}

femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i, j, n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc];
        }
    }
    
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);
    
    int index = 0;
    int nBoundary = 0;
    
    for (i = 0; i < theEdges->nEdge; i++) {
        if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i], &edges[i + 1]) != 0) {
            edges[index] = edges[i];
            nBoundary++;
        }
        else {
            edges[index] = edges[i];
            edges[index].elem[1] = edges[i + 1].elem[0];
            i = i + 1;
        }
        index++;
    }
    
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;
    for (i = 0; i < theEdges->nEdge; ++i) {
        printf("%6d : %4d %4d : %4d %4d \n", i,
               theEdges->edges[i].node[0], theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0], theEdges->edges[i].elem[1]);
    }
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*)edgeOne)->node;
    int *nodeTwo = ((femEdge*)edgeTwo)->node;
    int  diffMin = fmin(nodeOne[0], nodeOne[1]) - fmin(nodeTwo[0], nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0], nodeOne[1]) - fmax(nodeTwo[0], nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1;
    return  0;
}

femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);
    
    return theSystem;
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;
    double *elem = malloc(sizeof(double) * size * (size + 1));
    mySystem->A = malloc(sizeof(double*) * size);
    mySystem->B = elem;
    mySystem->A[0] = elem + size;
    mySystem->size = size;
    for (i = 1; i < size; i++)
        mySystem->A[i] = mySystem->A[i - 1] + size;
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i, size = mySystem->size;
    for (i = 0; i < size*(size + 1); i++)
        mySystem->B[i] = 0;
}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    
    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;
    
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");
            else               printf(" %+.1e", A[i][j]);
        printf(" :  %+.1e \n", B[i]);
    }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k = 0; k < size; k++) {
        if (fabs(A[k][k]) <= 1e-16) {
            printf("Pivot index %d  ", k);
            printf("Pivot value %e  ", A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        for (i = k + 1; i < size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k + 1; j < size; j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor;
        }
    }
    
    /* Back-substitution */
    
    for (i = size - 1; i >= 0; i--) {
        factor = 0;
        for (j = i + 1; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = (B[i] - factor) / A[i][i];
    }
    
    return(mySystem->B);
}

void  femFullSystemConstrain(femFullSystem *mySystem,
                             int myNode, double myValue)
{
    double  **A, *B;
    int     i, size;
    
    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;
    
    for (i = 0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0;
    }
    
    for (i = 0; i < size; i++)
        A[myNode][i] = 0;
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}


femGrains *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut, double gamma)
{
    int i,nContact = n*(n-1)/2;
    
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = n;
    theGrains->radiusIn = radiusIn;
    theGrains->radiusOut = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -9.81;
    theGrains->gamma = gamma;
    
    theGrains->elem = malloc(n*sizeof(int)); /* Tableau d'appartenance a un element */
    theGrains->x  = malloc(n*sizeof(double));
    theGrains->y  = malloc(n*sizeof(double));
    theGrains->vx = malloc(n*sizeof(double));
    theGrains->vy = malloc(n*sizeof(double));
    theGrains->r  = malloc(n*sizeof(double));
    theGrains->m  = malloc(n*sizeof(double));
    theGrains->dvBoundary = malloc(n * sizeof(double));
    theGrains->dvContacts = malloc(nContact * sizeof(double));
    
    for(i = 0; i < n; i++) {
        theGrains->r[i] = r;
        theGrains->m[i] = m;
        theGrains->x[i] = (i%5) * r * 2.5 - 5 * r + 1e-8;
        theGrains->y[i] = (i/5) * r * 2.5 + 2 * r + radiusIn;
        theGrains->vx[i] = 0.0;
        theGrains->vy[i] = 0.0;
        theGrains->dvBoundary[i] = 0.0; }
    
    for(i = 0; i < nContact; i++)
        theGrains->dvContacts[i] = 0.0;
    
    
    return theGrains;
}

femGrains *femGrainsCreateTiny(double radiusIn, double radiusOut)
{
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = 2;
    theGrains->radiusIn   = radiusIn;
    theGrains->radiusOut  = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -10;
    theGrains->gamma = 0.0;
    
    theGrains->x  = malloc(2 * sizeof(double));
    theGrains->y  = malloc(2 * sizeof(double));
    theGrains->vx = malloc(2 * sizeof(double));
    theGrains->vy = malloc(2 * sizeof(double));
    theGrains->r  = malloc(2 * sizeof(double));
    theGrains->m  = malloc(2 * sizeof(double));
    theGrains->dvBoundary = malloc(2 * sizeof(double));
    theGrains->dvContacts = malloc(sizeof(double));
    
    theGrains->r[0] = 0.1;
    theGrains->r[1] = 0.1;
    theGrains->m[0] = 1.0;
    theGrains->m[1] = 1.0;
    theGrains->x[0] = 0.0;
    theGrains->x[1] = 0.0;
    theGrains->y[0] = -radiusOut + 0.3;
    theGrains->y[1] = -radiusOut + 0.1;
    theGrains->vx[0] = 0.0;
    theGrains->vx[1] = 0.0;
    theGrains->vy[0] = 0.0;
    theGrains->vy[1] = 0.0;
    theGrains->dvBoundary[0] = 0.0;
    theGrains->dvBoundary[1] = 0.0;
    
    theGrains->dvContacts[0] = 0.0;
    
    return theGrains;
}

void femGrainsFree(femGrains *theGrains)
{
    free(theGrains->elem);
    free(theGrains->x);
    free(theGrains->y);
    free(theGrains->vx);
    free(theGrains->vy);
    free(theGrains->r);
    free(theGrains->m);
    free(theGrains->dvBoundary);
    free(theGrains->dvContacts);
    free(theGrains);
}

double femMin(double *x, int n)
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++)
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++)
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femErrorScan(int test, int line, char *file)
{
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}
