#include"fem.h"

# ifndef NOCONTACTITERATE


double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
    int i,j,iContact;
    int n = myGrains->n;
    
    
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;
    
    double error = 0.0;
    double gap,rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;
    
    error = 0;
    iContact = 0;
    for(i = 0; i < n; i++) {
        for(j = i+1; j < n; j++) {
            rx = (x[j]-x[i]);
            ry = (y[j]-y[i]);
            rr = sqrt(rx*rx+ry*ry);
            nx = rx/rr;
            ny = ry/rr;
            if (iter == 0) {
                dv = dvContacts[iContact]; }
            else {
                vn = (vx[i]-vx[j])*nx + (vy[i]-vy[j])*ny ;
                gap = rr - (r[i]+r[j]);
                dv = fmax(0.0, vn + dvContacts[iContact] - gap/dt);
                dv = dv - dvContacts[iContact];
                dvContacts[iContact] += dv;
                error = fmax(fabs(dv),error); }
            vx[i] -= dv * nx * m[j] / ( m[i] + m[j] );
            vy[i] -= dv * ny * m[j] / ( m[i] + m[j] );
            vx[j] += dv * nx * m[i] / ( m[i] + m[j] );
            vy[j] += dv * ny * m[i] / ( m[i] + m[j] );
            iContact++; }}
    
    for(i = 0; i < n; i++) {
        rr = sqrt(x[i]*x[i]+y[i]*y[i]);
        nx = x[i]/rr;
        ny = y[i]/rr;
        if (iter == 0) {
            dv = dvBoundary[i]; }
        else {
            vn = vx[i]*nx + vy[i]*ny ;
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap/dt);
            gap = rr - rIn - r[i];
            dvIn  = fmax(0.0,-vn - dvBoundary[i] - gap/dt);
            dv = dvOut - dvIn - dvBoundary[i];
            dvBoundary[i] += dv;
            error = fmax(fabs(dv),error); }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny; }
    return error;
}

# endif


# ifndef NOUPDATE


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *Problemu, femPoissonProblem *Problemv)
{
    int i, j, iElem;
    int n = myGrains->n;
    
    double *x = myGrains->x;
    double *y = myGrains->y;
    double *m = myGrains->m;
    double *vy = myGrains->vy;
    double *vx = myGrains->vx;
    double gamma = myGrains->gamma;
    double gy = myGrains->gravity[1];
    
    
    double xt[3], yt[3], uLoc[3], vLoc[3],tau[3];
    int mapt[3];
    double uxy = 0.0, vxy = 0.0;
    
    //-1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
    for (i = 0; i < n; i++) {
        for (iElem = 0; iElem < Problemu->mesh->nElem; iElem++) {
            femMeshLocal(Problemu->mesh, iElem, mapt, xt, yt);
            double area = Problemu->area[iElem];
            
            for (j = 0; j < 3; j++) {
                
                uLoc[j] = Problemu->system->B[mapt[j]];
                vLoc[j] = Problemv->system->B[mapt[j]];
            }
            if( WithinTriangle(xt, yt, x[i], y[i],area) ){
                invefctforme(Problemu, x[i], y[i], xt, yt,tau);
            
                uxy = tau[0] * uLoc[0] + tau[1] * uLoc[1] + tau[2] * uLoc[2];
                vxy = tau[0] * vLoc[0] + tau[1] * vLoc[1] + tau[2] * vLoc[2];
            }
            
            vx[i] += - gamma * dt / m[i] * (vx[i] - uxy);
            vy[i] += (m[i] * gy - gamma * (vy[i] - vxy)) * dt / m[i];
            
        }
    }
    
    //-2- Correction des vitesses pour tenir compte des contacts
    
    
    int iter = 0;
    double error;
    
    do {
        error = femGrainsContactIterate(myGrains, dt, iter);
        iter++;
    } while ((error > tol / dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d  error = %14.7e n", iter - 1, error);
    
    
    //-3- Calcul des nouvelles positions sans penetrations de points entre eux
    
    
    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
    }
}

# endif

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh = femMeshRead(filename);
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
    theProblem->area = malloc(sizeof(double)*theProblem->mesh->nElem);
    ComputeArea(theProblem);
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
    free(theProblem->area);
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

void femPoissonSolve(femPoissonProblem *theProblemU,femPoissonProblem *theProblemV, femGrains *myGrains)
{
    femMesh *theMesh = theProblemU->mesh;
    femEdges *theEdges = theProblemU->edges;
    femIntegration *theRule = theProblemU->rule;
    femDiscrete *theSpace = theProblemU->space;
    
    femFullSystem *theSystemU = theProblemU->system;
    femFullSystem *theSystemV = theProblemV->system;
    
    double vext = 5.0;
    double radiusIn = myGrains->radiusIn;
    double radiusOut= myGrains->radiusOut;
    double mu = 1e-1;
    double *xg = myGrains->x;
    double *yg = myGrains->y;
    double *vy = myGrains->vy;
    double *vx = myGrains->vx;
    int nGrains = myGrains->n;
    double gamma = myGrains->gamma;
    
    double x[3], y[3], phi[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3], tau[3];
    int iElem, iInteg, iEdge, i, j, k, map[3];
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femMeshLocal(theMesh, iElem, map, x, y);
        double area = theProblemU->area[iElem];
        
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            
            double jac = fabs((x[0]-x[1]) * (y[0]-y[2]) - (x[0]-x[2]) * (y[0]-y[1]));
            for (i = 0; i < 3; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    theSystemU->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight*mu;
                    theSystemV->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight*mu;
                }
            }
            for (k = 0; k < nGrains; k++) {
                if(WithinTriangle(x,y,xg[k],yg[k],area)){
                    invefctforme(theProblemU, xg[k], yg[k], x, y, tau);
                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            theSystemU->A[map[i]][map[j]] += gamma*tau[i] * tau[j];
                            theSystemV->A[map[i]][map[j]] += gamma*tau[i] * tau[j];
                        }
                        theSystemU->B[map[i]] += gamma * tau[i] * vx[k];
                        theSystemV->B[map[i]] += gamma * tau[i] * vy[k];
                    }
                }
            }
        }
    }
    
    for (iEdge = 0; iEdge < theEdges->nEdge; iEdge++) {
        if (theEdges->edges[iEdge].elem[1] < 0) {
            for (i = 0; i < 2; i++) {
                int iNode = theEdges->edges[iEdge].node[i];
                double xloc = theMesh->X[iNode];
                double yloc = theMesh->Y[iNode];
                double rMid = (radiusIn + radiusOut)/2;
                if ( sqrt(xloc*xloc + yloc * yloc) < rMid){
                    femFullSystemConstrain(theSystemU, iNode, 0.0);
                    femFullSystemConstrain(theSystemV, iNode, 0.0);}
                else{
                    femFullSystemConstrain(theSystemU, iNode, yloc*vext/radiusOut);
                    femFullSystemConstrain(theSystemV, iNode, -xloc*vext/radiusOut);}
            }
        }
    }
    
    femFullSystemEliminate(theSystemU);
    femFullSystemEliminate(theSystemV);
}



# endif






