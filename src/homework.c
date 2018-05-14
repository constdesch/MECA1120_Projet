#include"fem.h"
#include <math.h>

# ifndef NOPOISSONCREATE

void area (femPoissonProblem *theProblem){
    femMesh *theMesh = theProblem->mesh;
    int elem;
    int map[3];
    double px[3], py[3];
    for (elem=0; elem<theMesh->nElem;elem++){
        femMeshLocal(theMesh, elem, map, px, py);
        double p0x = px[0], p1x = px[1], p2x = px[2], p0y = py[0], p1y = py[1], p2y = py[2];
        theMesh->area[elem] = -p1y * p2x + p0y * (p2x - p1x) + p0x * (p1y - p2y) + p1x * p2y ;
        //printf("area : %f ",theMesh->area[elem]);
    }
    //printf("\n");
}

int withinTriangle(double xc, double yc, double area, double px[3],double py[3]) {
 double p0x = px[0], p1x = px[1], p2x = px[2], p0y = py[0], p1y = py[1], p2y = py[2];
 
 double s = p0y * p2x - p0x * p2y + (p2y - p0y) * xc + (p0x - p2x)*yc;
 double t = p0x * p1y - p0y * p1x + (p0y - p1y) * xc + (p1x - p0x)*yc;
 
 if ( (s<0) != (t<0))
 return 0;
 
 if (area <0.0 ){
 s = -s;
 t = -t;
 area = -area;
 }
 return ( s>0 && t>0 && (s+t) <= area );
 }
/*
int withinTriangle(double xc, double yc, double area, double px[3], double py[3]) {
    double    jac1 = (px[1] - xc)*(py[2] - yc) -(px[2] - xc)*(py[1] - yc);
    double jac2 = (xc - px[0])*(py[2] - py[0]) - (px[2] - px[0])*(yc - py[0]);
    double    jac3 = (px[1] - px[0])*(yc - py[0]) - (xc - px[0])*(py[1] - py[0]);
    return(jac1 > 0 && jac2 > 0 && jac3 > 0);
}*/
void indexoftriangle(femPoissonProblem* theProblem, femGrains* theGrains) {
    int i,j;
    femMesh *theMesh = theProblem->mesh;
    int map[3];
    double px[3], py[3];
    for (i = 0; i < theGrains->n; i++) {
        int ExitFlag = 0;
        for (j = 0; j < theMesh->nElem&&ExitFlag==0; j++) {
            femMeshLocal(theMesh, j, map, px, py);
            if (withinTriangle(theGrains->x[i], theGrains->y[i], theMesh->area[j], px, py) == 1) {
                theGrains->elem[i] = j;
                ExitFlag++;
            }
        }
        printf("%d ", theGrains->elem[i]);
        
    }
    printf("fini\n");
}
femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh = femMeshRead(filename);
    area(theProblem);
    femMeshClean(theProblem->mesh);
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
    }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
    }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}



# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}

# endif
# ifndef NOMESHLOCAL

void femMeshLocal(const femMesh *theMesh, const int iElem, int *map, double *x, double *y)
{
    int j, nLocal = theMesh->nLocalNode;
    
    for (j = 0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal + j];
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
    }
}

# endif
# ifndef NOPOISSONSOLVE
void femPoissonSolve(femPoissonProblem *theProblemU, femPoissonProblem *theProblemV, femGrains* theGrains)
{

    femMesh *theMeshU = theProblemU->mesh;
    femEdges *theEdgesU = theProblemU->edges;
    femFullSystem *theSystemU = theProblemU->system;
    femIntegration *theRuleU = theProblemU->rule;
    femDiscrete *theSpaceU = theProblemU->space;
    
    femMesh *theMeshV = theProblemV->mesh;
    femEdges *theEdgesV = theProblemV->edges;
    femFullSystem *theSystemV = theProblemV->system;
    femIntegration *theRuleV = theProblemV->rule;
    femDiscrete *theSpaceV = theProblemV->space;
    
   
    int e,f;
    for(e = 0; e < theSystemU->size;e++){
        for(f = 0; f < theSystemU->size;f++){
            theProblemU->system->A[e][f] = 0.0;
            theProblemV->system->A[e][f] = 0.0;
        }
    }
    
    
    if (theSpaceU->n > 4) Error("Unexpected discrete space size !");
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4], tau[4], gamma = 0.5, presentU, presentV, pos[2], mu = 1.0;
    int iElem, iInteg, iEdge, i, j, map[4], k, exitflag;
    for (iElem = 0; iElem < theMeshU->nElem; iElem++) {
        femMeshLocal(theMeshU, iElem, map, x, y);
        presentU = 0.0;
        presentV = 0.0;
        exitflag = 0;
        for (k= 0; k < theGrains->n && exitflag==0; k++) {
            if (iElem == theGrains->elem[k]) {
                presentU = theGrains->vx[k];
                presentV = theGrains->vy[k];
                pos[0] = theGrains->x[k];
                pos[1] = theGrains->y[k];
                exitflag++;
            }
        }
        for (iInteg = 0; iInteg < theRuleU->n; iInteg++) {
            double xsi = theRuleU->xsi[iInteg];
            double eta = theRuleU->eta[iInteg];
            if (exitflag != 0) {
                double *tab;
                tab=femDiscreteinvetaxsi(pos[0], pos[1], x, y);
                femDiscretePhi2(theSpaceU, tab[0],tab[1] ,tau);
                free(tab);
            }
            double weight = theRuleU->weight[iInteg];
            femDiscretePhi2(theSpaceU, xsi, eta, phi);
            femDiscreteDphi2(theSpaceU, xsi, eta, dphidxsi, dphideta);
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            double xloc = 0.0;
            double yloc = 0.0;
            for (i = 0; i < theSpaceU->n; i++) {
                xloc   += x[i] * phi[i];
                yloc   += y[i] * phi[i];
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpaceU->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            
            for (i = 0; i < theSpaceU->n; i++) {
                for (j = 0; j < theSpaceU->n; j++) {
                    if (exitflag!=0) {
                        theSystemU->A[map[i]][map[j]] += mu * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight+ gamma*tau[i] * tau[j];
                        theSystemV->A[map[i]][map[j]] += mu * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight+ gamma*tau[i] * tau[j];
                    }
                    else{
                        theSystemU->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight;
                        theSystemV->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight;}
                }
                theSystemU->B[map[i]] += gamma*tau[i] * presentU;
                theSystemV->B[map[i]] += gamma*tau[i] * presentV;
            }
        }
    }
    
    for (iEdge = 0; iEdge < theEdgesU->nEdge; iEdge++) {
        if (theEdgesU->edges[iEdge].elem[1] < 0) {
            for (i = 0; i < 2; i++) {
                double rOut = 2.0;//theGrains->radiusOut;
                double rIn = 0.4;//theGrains->radiusIn;
                int iNode = theEdgesU->edges[iEdge].node[i];
                double xlocU = theMeshU->X[iNode];
                double ylocU = theMeshU->Y[iNode];
                double xlocV = theMeshV->X[iNode];
                double ylocV = theMeshV->Y[iNode];
                double r = sqrt(xlocU*xlocU + ylocU * ylocU);
                double vext = 3.0;
                if (r <(rOut-rIn)) {
                    double vx = 0.0;
                    femFullSystemConstrain(theSystemU, iNode, 0);
                    femFullSystemConstrain(theSystemV, iNode, 0);
                }
                else {
                    double vx =vext *ylocU/rOut ;
                    femFullSystemConstrain(theSystemU, iNode, vx);

                    double vy =- vext *xlocV/rOut;
                    femFullSystemConstrain(theSystemV, iNode, vy);
                }
            }
        }
    }
    
    femFullSystemEliminate(theSystemU);
    femFullSystemEliminate(theSystemV);
}


# endif

# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
    int i, j, iContact;
    int n = myGrains->n;
    
    
    double *x = myGrains->x;
    double *y = myGrains->y;
    double *m = myGrains->m;
    double *r = myGrains->r;
    double *vy = myGrains->vy;
    double *vx = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn = myGrains->radiusIn;
    double rOut = myGrains->radiusOut;
    
    double error = 0.0;
    double gap, rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;
    
    error = 0;
    iContact = 0;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            rx = (x[j] - x[i]);
            ry = (y[j] - y[i]);
            rr = sqrt(rx*rx + ry * ry);
            nx = rx / rr;
            ny = ry / rr;
            if (iter == 0) {
                dv = dvContacts[iContact];
            }
            else {
                vn                    = (vx[i] - vx[j])*nx + (vy[i] - vy[j])*ny;
                gap                   = rr - (r[i] + r[j]);
                dv                    = fmax(0.0, vn + dvContacts[iContact] - gap / dt);
                dv                    = dv - dvContacts[iContact];
                dvContacts[iContact] += dv;
                error                 = fmax(fabs(dv), error);
            }
            vx[i] -= dv * nx * m[j] / (m[i] + m[j]);
            vy[i] -= dv * ny * m[j] / (m[i] + m[j]);
            vx[j] += dv * nx * m[i] / (m[i] + m[j]);
            vy[j] += dv * ny * m[i] / (m[i] + m[j]);
            iContact++;
        }
    }
    
    for (i = 0; i < n; i++) {
        rr = sqrt(x[i] * x[i] + y[i] * y[i]);
        nx = x[i] / rr;
        ny = y[i] / rr;
        if (iter == 0) {
            dv = dvBoundary[i];
        }
        else {
            vn = vx[i] * nx + vy[i] * ny;
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap / dt);
            gap = rr - rIn - r[i];
            dvIn = fmax(0.0, -vn - dvBoundary[i] - gap / dt);
            dv = dvOut - dvIn - dvBoundary[i];
            dvBoundary[i] += dv;
            error = fmax(fabs(dv), error);
        }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny;
    }
    return error;
}

# endif
# ifndef NOUPDATE


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *theProblemU,femPoissonProblem* theProblemV)
{
  
    int i;
    int n = myGrains->n;
    
    double *x = myGrains->x;
    double *y = myGrains->y;
    double *m = myGrains->m;
    double *vy = myGrains->vy;
    double *vx = myGrains->vx;
    double gamma = myGrains->gamma;
    double gx = myGrains->gravity[0];
    double gy = myGrains->gravity[1];
    double xc[3], yc[3],uu[3],vv[3];
    int map[3];
    int iElem;
    int j;
    double u = 0, v = 0;
    //
    // -1- Calcul des nouvelles vitesses des grains sur base de la gravit et de la trainee
    //
    
    for (i = 0; i < n; i++) {
        
        for (iElem = 0; iElem < theProblemU->mesh->nElem; iElem++) {
            femMeshLocal(theProblemU->mesh, iElem, map, xc, yc);
            for (j = 0; j < 3; j++) {
                uu[j] = theProblemU->system->B[map[j]];
                vv[j] = theProblemV->system->B[map[j]];
            }
            if (withinTriangle(x[i], y[i], theProblemU->mesh->area[iElem], xc, yc) == 1) {
                double *tab;
                double tau[3];
                tab = femDiscreteinvetaxsi(x[i], y[i], xc, yc);
                femDiscretePhi2(theProblemU->space, tab[0], tab[1], tau);
                free(tab);
                printf("hey bro [%f %f %f]\n", tau[0], tau[1], tau[2]);
                u = tau[0] * uu[0] + tau[1] * uu[1] + tau[2] * uu[2];
                v = tau[0] * vv[0] + tau[1] * vv[1] + tau[2] * vv[2];
            }
            
            vx[i] += -gamma * dt / m[i] * (vx[i] - u);
            vy[i] += (m[i] * gy - gamma * (vy[i] - v)) * dt / m[i];
        }
    }
    
    
    //
    // -2- Correction des vitesses pour tenir compte des contacts
    //
    
    int iter = 0;
    double error;
    
    do {
        error = femGrainsContactIterate(myGrains, dt, iter);
        iter++;
    } while ((error > tol / dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n", iter - 1, error);
    
    //
    // -3- Calcul des nouvelles positions sans penetrations de points entre eux
    //
    
    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
    }
    indexoftriangle(theProblemU, myGrains);
    indexoftriangle(theProblemV, myGrains);
}
# endif
