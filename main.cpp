/***************************************************************************
					  main.cpp  -  description
						 -------------------
copyright            : (C) 2006 by Joan Anglada
email                : xxxxxx
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <unistd.h>
#include "MTRand.h"     //UTILITZAREM EL MERSENNE TWISTER RNG PQ TÉ UN PERIODE QUE FLIPES
// LA RUTINA RAND() DEL C++ TE UN PERIODE DE NOMÉS  2³¹
//EN CANVI AQUEST TÉ UN PERIODE DE 2^19937-1 (diuen)
//I A MÉS, ÉS (+?)  RÀPID

using namespace std;
//FUNCTION PROTOTYPES:


void initialize();              //LEGEIX parametres en para.in i crea lattice
int randomInt();                //NOMBRES RANDOM -1 o +1
void step();              // Fa un montecarlo step, canvia energia i m.
void plotIt();                  //dibuixa a trvés de la consola la xarxa
void renderScene(void);         //pinta a través de la finestra;
void calcVariables();           //Calcula energia i magnetització inicial
void calc(int a, int b);        //calcula energia d'un lloc
void correlations();

int getS(int x, int y, int dx, int dy);

void writeResults();

MTRand mtrand1;
//GLOBAL VARIABLES:

long int seed;

int trans;
int N, NC;
int steps;
int Jo;
int count = 0;
int latt[1000][1000];
float correl[50];
float correlAv[50];
int nsteps = 0;


float eAv;
float mAv;
float m2;
float m2Av;
float chi;
float energy;
float m;
float cvAv;
float cv;
float energsum2;
double T;
double h;
float H1, H2, Hij, Mij;
//double Z1,Z2;



int main(int argc, char **argv) {
	MTRand mtrand1(5);

	float question = -1;//=-1;
	//char conf="S";
	cout<<"Voleu representació gràfica?(1=S o -1=N)";
	cin>>question;

	if (question > 0) {

		initialize();

		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		glutCreateWindow("ISING");
		glutInitWindowSize(900, 900);
		glutDisplayFunc(renderScene);
		glutIdleFunc(renderScene);
		glEnable(GL_DEPTH_TEST);
		glutMainLoop();

	} else {

		initialize();
		int ncor = 0;

		//cout<<"hi\n";
		calcVariables();
		// cout<<"magnetització inicial="<<m<<"\n";
		nsteps = 0;
		for (int i = 0; i < steps; i++) {
			step();
			nsteps++;

			if (i > trans && (ncor % 2 == 0)) {
				// correlations();
				ncor++;
			}
		}

		eAv = eAv / (float) (steps - trans);
		energy = eAv / N;
		mAv = mAv / (float) (steps - trans);
		m = mAv / N;
		cvAv = cvAv / (float) (steps - trans);
		cvAv = cvAv / (N * N);
		cv = (cvAv - energy * energy) / (T * T);
		m2Av = m2Av / (float) (steps - trans);
		m2Av = m2Av / (N * N);
		chi = (m2Av - m * m) / (T * T);


// for (int l=0;l<12;l++) {correlAv[l]=correlAv[l]/(ncor);}


		//cout<<"temperatura, energia per spin i magnetització per spin:\n";
		cout << T << "     " << energy << "   " << m << "   " << cv << "      " << h << "      " << chi << "\n";
// writeResults();
	}
}

void initialize() {


	//  long int seed;
// OPEN FILE WITH PARAMETERS    
	int k = 0;
	double parameters[10];
	ifstream inFile("para.in");
	if (!inFile) {
		cerr << "Unable to open file datafile.txt";
		exit(1);   // call system to stop
	}

	inFile >> parameters[k];
	while (!inFile.eof()) { // keep reading until end-of-file
		// cout << " " << parameters[k] ;
		k++;
		inFile >> parameters[k]; // sets EOF flag if no value found
		// cout<<parameters[k]
	}
	inFile.close();

//SET PARAMETERS    
	NC = (int) parameters[0];         //LLEGEIX DIMENSIÓ COSTAT XARXA;
	if (NC > 1000) {
		cout << "NO SIGUIS CAPDELLUÇ! ENTRA UNA XARXA DE COSTAT<1000\n";
		exit(1);
	}
	N = NC * NC;
	steps = (int) parameters[1];      // LLEGEIX STEPS TOTALS MC (en modalitat no gràfica)
	trans = (int) parameters[2];      //LONGITUD DEL TRANSITORI;
	h = (double) parameters[3];                //EXTERNAL FIELD
	T = (double) parameters[4];         //TEMPERAATURE
	seed = (int) parameters[5];       //SEED DEL RNG
	Jo = (int) parameters[6];         //CONSTANT ACOBLAMENT J


	// INIT LATTICE
	int g = ((mtrand1.rand() - 0.5) < 0) ? 1 : -1;
	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {

			if (T < 2.269) latt[i][j] = g;
			else latt[i][j] = ((mtrand1.rand() - 0.5) < 0) ? 1 : -1;

		}
	}
	m2 = 0;
	m2Av = 0;
	chi = 0;
	eAv = 0;
	mAv = 0;
	cvAv = 0;
	cv = 0;
	// energ2sum=0;
}

void calc(int a, int b) {
	float sN = 0;
	float sE = 0;
	float sS = 0;
	float sW = 0;
	//TENIM EN COMPTE LES CONDICIONS PERIÒDIQUES DE CONTORN
	if (a == 0) {
		sW = latt[NC - 1][b];
	} else {
		sW = latt[a - 1][b];
	}

	if (a == NC - 1) {
		sE = latt[0][b];
	} else {
		sE = latt[a + 1][b];
	}

	if (b == 0) {
		sN = latt[a][NC - 1];
	} else {
		sN = latt[a][b - 1];
	}

	if (b == NC - 1) {
		sS = latt[a][0];
		//   yes++;
	} else {
		sS = latt[a][b + 1];
	}

	//  Z1=exp(1*(Jo*(sW+sE+sN+sS)+h)/T);
	//Z2=exp(-1*(Jo*(sW+sE+sN+sS)+h)/T);
	Hij = latt[a][b] * (Jo * (sW + sE + sN + sS) + h);
	//H1=1*(Jo*(sW+sE+sN+sS)+h);
	//  H2=-1*(Jo*(sW+sE+sN+sS));

}

void step() {

	int a;
	int b;
	int sign;
	double p1 = 0;
	double p2 = 0;
	float de;
	for (int y = 0; y < N; y++) { //per cada spin

		//CHOOSE SPIN AT SITE a,b
		a = mtrand1.randInt(NC);
		b = mtrand1.randInt(NC);
		//  sign=latt[a][b];
		//CALCULEM PROB +1 ó -1
		calc(a, b);
		//    p1=Z1/(Z1+Z2);
		//  p2=Z2/(Z1+Z2);
		de = 2 * Hij; //CANVI HIPOTETIC EN L'ENERGIA
		if (de < 0 || (exp(-de / T) > mtrand1.rand())) {

			energy = energy + de;
			latt[a][b] = -latt[a][b];
			m = m + 2 * latt[a][b];
			//energ2sum=energ2sum+energy**2
		}

	}

	if (nsteps >= trans) {

		//si HA PASSAT EL TRANSITORI, GUARDEM DADES energia I m.
		eAv = eAv + energy;
		mAv = mAv + m;
		cvAv = cvAv + energy * energy;
		m2Av = m2Av + m * m;
	}

}

void plotIt() {

	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {
			latt[i][j] = 1;//randomInt();

			if (latt[i][j] >= 0) cout << "+" << "  ";
			else cout << "-" << "  ";
			if (j == NC - 1) cout << "\n";
		}
	}
}

void calcVariables() {

	float haux = 0;
	float maux = 0;

	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {
			calc(i, j);
			haux = haux - Hij;
			maux = maux + latt[i][j];

		}
	}

	//energy=-haux;
	m = maux;
	energy = haux / 2;//eAv+haux;
	mAv = 0;//mAv+maux;
	eAv = 0;
}


void renderScene(void) {


	float d;
	float L = 2;
	d = L / NC;
	step();

	nsteps++;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	for (int k = 0; k < NC; k++) {

		for (int j = 0; j < NC; j++) {

			if (latt[k][j] > 0) {
				glColor3f(0.0, 0.0, 0.0);
			} else {
				glColor3f(0.33, 0.7, 0.6);
			}
			glRectf(-1 + k * d, 1 - (j + 1) * d, -1 + (k + 1) * d, 1 - j * d); //DIBUIXA RECTANGLE


		}
	}
	glutSwapBuffers();          //FA LA ANIMACIÓ MÉS SMOOTH

}

void correlations() {
	int s = 0;
	for (int k = 1; k < 12; k++) {
		if (k == 1) { s = 2; }
		else { s = 2 * k; }
		for (int i = 0; i < NC; i++) {
			for (int j = 0; j < NC; j++) {

				correl[k - 1] = correl[k - 1] + latt[i][j] * (getS(i, j, s, 0) + getS(i, j, -s, 0) + getS(i, j, 0, s) +
															  getS(i, j, 0, -s));

			}
		}
		correl[k - 1] = correl[k - 1] / (4 * N);
		correlAv[k - 1] = correlAv[k - 1] + correl[k - 1];
	}
}

int getS(int x, int y, int dx, int dy) {

	// return latt[x+dx][y+dy];

	if ((x + dx) >= NC) { return latt[x + dx - NC][y]; }
	else {
		if ((x + dx) < 0) { return latt[NC + x + dx - 1][y]; }
		else { return latt[x + dx][y + dy]; }
	}

	if ((y + dy) >= NC) { return latt[x][y + dy - NC]; }
	else {
		if ((y + dy) < 0) { return latt[x][NC + y + dy]; }
		else { return latt[x + dx][y + dy]; }
	}

}


void writeResults() {
	int s;
	ofstream results("results.dat");
	ofstream corr("correlation.dat");

	results << "after " << steps << "done\n";
	results << "************************************************************************\n";
	results << "temperatura, energia per spin i magnetització per spin:\n";
	results << T << "     " << energy << "   " << m << "\n";
	results << "CORRELACIONS\n";

	for (int k = 1; k < 12; k++) {
		if (k == 1) { s = 2; }
		else { s = 2 * k; }
		results << "distancia " << s << "  " << correlAv[k - 1] << "\n";
		double lo = log(correlAv[k - 1]);
		corr << s << "    " << correlAv[k - 1] << "      " << lo << "\n";

	}
	results.flush();
	results.close();
	corr.flush();
	corr.close();

}
 
 
