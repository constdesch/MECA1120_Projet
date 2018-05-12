#include"fem.h"

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

void femPoissonSolve(femPoissonProblem *theProblem)
{
	femMesh *theMesh = theProblem->mesh;
	femEdges *theEdges = theProblem->edges;
	femFullSystem *theSystem = theProblem->system;
	femIntegration *theRule = theProblem->rule;
	femDiscrete *theSpace = theProblem->space;

	if (theSpace->n > 4) Error("Unexpected discrete space size !");
	double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
	int iElem, iInteg, iEdge, i, j, map[4];

	for (iElem = 0; iElem < theMesh->nElem; iElem++) {
		femMeshLocal(theMesh, iElem, map, x, y);
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
					theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j]
						+ dphidy[i] * dphidy[j]) * jac * weight;
				}
			}
			for (i = 0; i < theSpace->n; i++) {
				theSystem->B[map[i]] += phi[i] * jac *weight;
			}
		}
	}

	for (iEdge = 0; iEdge < theEdges->nEdge; iEdge++) {
		if (theEdges->edges[iEdge].elem[1] < 0) {
			for (i = 0; i < 2; i++) {
				int iNode = theEdges->edges[iEdge].node[i];
				double xloc = theMesh->X[iNode];
				double yloc = theMesh->Y[iNode];
				femFullSystemConstrain(theSystem, iNode, 0.0);
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
				vn = (vx[i] - vx[j])*nx + (vy[i] - vy[j])*ny;
				gap = rr - (r[i] + r[j]);
				dv = fmax(0.0, vn + dvContacts[iContact] - gap / dt);
				dv = dv - dvContacts[iContact];
				dvContacts[iContact] += dv;
				error = fmax(fabs(dv), error);
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


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
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
}

# endif