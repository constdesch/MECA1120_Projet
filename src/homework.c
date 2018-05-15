//Lac Adlane & Roelandts Antoine
#include"fem.h"


double area (double xa, double xb, double xc, double ya, double yb, double yc){
    return 1 /2 * fabs( (xb - xa)*(yc - ya) - (xc - xa)*(yb - ya) );
}

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
int DeltaPointTriangle(double *x, double *y, double xc, double yc, double Aire) {
    double aire1 = 1 / 2 * area(xc, x[1],x[2],yc,y[1],y[2]);
    double aire2 = 1 / 2 * area(x[0],xc,x[2],y[0],yc,y[2]);
    double aire3 = 1 / 2 * area(x[0],x[1],xc,y[0],y[1],yc);
    double s = aire1 + aire2 + aire3;
    return (Aire == s); 
}
/*int DeltaPointTriangle(double xa, double ya, double xb, double yb, double xc, double yc, double xp, double yp)
{
    
    double ABC = 1 / 2 * fabs((xb - xa)*(yc - ya) - (xc - xa)*(yb - ya));
    double PBC = 1 / 2 * fabs((xb - xp)*(yc - yp) - (xc - xp)*(yb - yp));
    double APC = 1 / 2 * fabs((xp - xa)*(yc - ya) - (xc - xa)*(yp - ya));
    double ABP = 1 / 2 * fabs((xb - xa)*(yp - ya) - (xp - xa)*(yb - ya));
    
    double Somme = PBC + APC + ABP;
    
    if (Somme == ABC) return 1;
    
    else return 0;
    
}*/



double invksi(double x, double y, double X1, double Y1, double X2, double Y2, double X3, double Y3)
{
    return ((x - X1)*(Y3 - Y1) - (y - Y1)*(X3 - X1)) / ((X2 - X1)*(Y3 - Y1) - (Y2 - Y1)*(X3 - X1));
    //return (y - Y1 - (Y2 - Y1)*inveta(x, y, X1, Y1, X2, Y2, X3, Y3)) / (Y2 - Y1);
}


double inveta(double x, double y, double X1, double Y1, double X2, double Y2, double X3, double Y3)
{
    return ((x - X1)*(Y1 - Y2) - (y - Y1)*(X1 - X2)) / ((X2 - X1)*(Y3 - Y1) - (Y2 - Y1)*(X3 - X1));
    
    
    //return (y - Y1 - (Y2 - Y1)*invksi(x, y, X1, Y1, X2, Y2, X3, Y3)) / (Y3 - Y1);
}


void invtau(double x, double y, double X1, double Y1, double X2, double Y2, double X3, double Y3, double* invtau)
{
    invtau[0] = 1-invksi(x, y, X1, Y1, X2, Y2, X3,Y3)-inveta(x, y, X1, Y1, X2, Y2, X3, Y3);
    invtau[1] = invksi(x, y, X1, Y1, X2, Y2, X3, Y3);
    invtau[2] = inveta(x, y, X1, Y1, X2, Y2, X3, Y3);
    
}



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
    double uxy, vxy;
    
    //-1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
    for (i = 0; i < n; i++) {
        for (iElem = 0; iElem < Problemu->mesh->nElem; iElem++) {
            femMeshLocal(Problemu->mesh, iElem, mapt, xt, yt);
            double area = Problemu->area[iElem];
            
            for (j = 0; j < 3; j++) {
                
                uLoc[j] = Problemu->system->B[mapt[j]];
                vLoc[j] = Problemv->system->B[mapt[j]];
            }
            //  ‡ priori sert ‡ qued ici DeltaPointTriangle(xt[0], yt[0], xt[1], yt[1], xt[2], yt[2], x[i], y[i])*
            
            invtau(x[i], y[i], xt[0], yt[0], xt[1], yt[1], xt[2], yt[2], tau);
            
            uxy = DeltaPointTriangle(xt, yt, x[i], y[i],area)*(tau[0] * uLoc[0] + tau[1] * uLoc[1] + tau[2] * uLoc[2]);
            vxy = DeltaPointTriangle(xt, yt,x[i], y[i],area)*(tau[0] * vLoc[0] + tau[1] * vLoc[1] + tau[2] * vLoc[2]);
            
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
    double *m = myGrains->m;
    double *vy = myGrains->vy;
    double *vx = myGrains->vx;
    int nGrains = myGrains->n;
    double gamma = myGrains->gamma;
    
    double x[3], y[3], phi[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3], tauetoile[3];
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
            double xloc = 0;
            double yloc = 0;
            for (i = 0; i < 3; i++) {
                xloc += x[i] * phi[i];
                yloc += y[i] * phi[i];
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
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
        }
        
        for (k = 0; k < nGrains; k++) {
            invtau(xg[k], yg[k], x[0], y[0], x[1], y[1], x[2], y[2], tauetoile);
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    theSystemU->A[map[i]][map[j]] += DeltaPointTriangle(x,y,xg[k],yg[k],area)*gamma*tauetoile[i] * tauetoile[j];
                    theSystemV->A[map[i]][map[j]] += DeltaPointTriangle(x,y,xg[k], yg[k],area)*gamma*tauetoile[i] * tauetoile[j];
                }
                theSystemU->B[map[i]] += DeltaPointTriangle(x,y, xg[k], yg[k],area)*gamma * tauetoile[i] * vx[k];
                theSystemV->B[map[i]] += DeltaPointTriangle(x,y, xg[k], yg[k],area)*gamma * tauetoile[i] * vy[k];
                
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






