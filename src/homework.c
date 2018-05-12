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
		
	}
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
void femPoissonSolve(femPoissonProblem *theProblem, int flag, femGrains* theGrains)
{
	femMesh *theMesh = theProblem->mesh;
	femEdges *theEdges = theProblem->edges;
	femFullSystem *theSystem = theProblem->system;
	femIntegration *theRule = theProblem->rule;
	femDiscrete *theSpace = theProblem->space;
    /*indexoftriangle(theProblem, theGrains);
    int z;
    for(z=0; theGrains->n;z++){
        printf("%d ",theGrains->elem[z]);
    }
    printf("\n");*/
    /*int elem;
    printf("area : ");
    for (elem=0; elem<theMesh->nElem;elem++){
        printf("%f ",theMesh->area[elem]);
    }
    printf("\n"); */
    

	if (theSpace->n > 4) Error("Unexpected discrete space size !");
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4], tau[4], gamma = 0.5;
    int iElem, iInteg, iEdge, i, j, map[4];

	for (iElem = 0; iElem < theMesh->nElem; iElem++) {
		femMeshLocal(theMesh, iElem, map, x, y);
		double present = 0.0;
		int k, exitflag = 0;
		for (k= 0; k < theGrains->n && exitflag==0; k++) {
			if (iElem == theGrains->elem[k]) {
				if (flag == 0) {
					present = theGrains->vx[k];
					exitflag++;
				}
				else {
					present = theGrains->vy[k];
					exitflag++;
				}
			}
		}
        /*printf("present : %f\n",present);*/
		for (iInteg = 0; iInteg < theRule->n; iInteg++) {
			double xsi = theRule->xsi[iInteg];
			double eta = theRule->eta[iInteg];
            double xx = x[iInteg];
            double yy = y[iInteg];
			double weight = theRule->weight[iInteg];
			femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscretePhi2(theSpace,  xx,  yy, tau);
			femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
			double dxdxsi = 0.0;
			double dxdeta = 0.0;
			double dydxsi = 0.0;
			double dydeta = 0.0;
			double xloc = 0.0;
			double yloc = 0.0;
			for (i = 0; i < theSpace->n; i++) {
				xloc += x[i] * phi[i];
				yloc += y[i] * phi[i];
				dxdxsi += x[i] * dphidxsi[i];
				dxdeta += x[i] * dphideta[i];
				dydxsi += y[i] * dphidxsi[i];
				dydeta += y[i] * dphideta[i];
			}
			double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
			for (i = 0; i < theSpace->n; i++) {
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
			}
			for (i = 0; i < theSpace->n; i++) {
				for (j = 0; j < theSpace->n; j++) {
						theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight
						+ gamma * tau[i] * tau[j];
				}
				if(present!=0.0)
                theSystem->B[map[i]] +=  gamma * tau[i]*present;
			}
		}
	}

	for (iEdge = 0; iEdge < theEdges->nEdge; iEdge++) {
		if (theEdges->edges[iEdge].elem[1] < 0) {
			for (i = 0; i < 2; i++) {
                double radiusOut = theGrains->radiusOut;
                double radiusIn = theGrains->radiusIn;
				int iNode = theEdges->edges[iEdge].node[i];
				double xloc = theMesh->X[iNode];
				double yloc = theMesh->Y[iNode];
                if (xloc * xloc + yloc*yloc <= radiusIn*radiusIn) femFullSystemConstrain(theSystem, iNode, 0.0);
                else {
                    double vext = 20.0;
                    if (flag == 0){
                        double vx = vext * yloc / radiusOut;
                        femFullSystemConstrain(theSystem, iNode, vx);
                    }
                    if (flag == 1){
                        double vy = -vext * xloc/radiusOut;
                        femFullSystemConstrain(theSystem, iNode, vy);
                    }
                }
			}
		}
	}

	femFullSystemEliminate(theSystem);
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


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *theProblem)
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

	// 
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravit et de la trainee
	//

	for (i = 0; i < n; i++) {
		double fx = m[i] * gx - gamma * vx[i];
		double fy = m[i] * gy - gamma * vy[i];
		vx[i] += fx * dt / m[i];
		vy[i] += fy * dt / m[i];
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
    indexoftriangle(theProblem, myGrains);
}
# endif
