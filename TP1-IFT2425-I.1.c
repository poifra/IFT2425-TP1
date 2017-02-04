//------------------------------------------------------
// module  : Tp-IFT2425-I.1.c
// author  : François Poitras & Guillaume Noël-Martel
// date    :
// version : 1.0
// language: C
// note    :
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <new>
#include <time.h>

//------------------------------------------------
// DEFINITIONS -----------------------------------
//------------------------------------------------
#define CARRE(X) ((X)*(X))
#define CUBE(X)  ((X)*(X)*(X))

float f(float x)
{
	return 0.0;
}
float df(float x, float epsilon)
{
	return (f(x + epsilon) - f(x))/epsilon;
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char** argv)
{
	//---------------------------
	//Algorithme NEWTON
	//---------------------------

	//implementer ici
	float epsilon = 0.000001;
	float x0 = 0.25, x1;
	int maxIter = 10;
	for (int i = 0; i < maxIter; i++)
	{

	}

	//retour sans probleme
	printf("\nFini...\n\n\n");
	return 0;
}