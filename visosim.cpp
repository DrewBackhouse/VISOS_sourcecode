#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <emscripten/bind.h>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <malloc.h>
#include "NeutrinoPropagator.h"
#include "BargerPropagator.h"
#include "EarthDensity.h"
#include "mosc.h"
#include "mosc3.h"

using namespace std;

#ifndef EM_PORT_API
#    if defined(__EMSCRIPTEN__)
#        include <emscripten.h>
#        if defined(__cplusplus)
#            define EM_PORT_API(rettype) extern "C" rettype EMSCRIPTEN_KEEPALIVE
#        else
#            define EM_PORT_API(rettype) rettype EMSCRIPTEN_KEEPALIVE
#        endif
#    else
#        if defined(__cplusplus)
#            define EM_PORT_API(rettype) extern "C" rettype
#        else
#            define EM_PORT_API(rettype) rettype
#        endif
#    endif
#endif

EM_PORT_API(double*) Propagate(double num0, double num1, double num2, double num3, double num4, double num5, double energy, double totaldistance ,int nutype, double density ) 
{
  //vector<double> X;
  //vector<double> Y;
  //vector<double> D;
  const int npoint = 44000;
  double* B = (double*)malloc(npoint * 8 * 3);

  BargerPropagator *bNunue = new BargerPropagator();
  bNunue->UseMassEigenstates(false);

  //set PMNS can't be global, otherwise will be rewritten by other instances of bNunue
  const double non_doubled_angle = true;
  // Set up 6 oscpars
  //(num0)= sin^2(theta12)
  //(num1)= sin^2(theta23)
  //(num2)=sin^2(theta13)
  //(num3)= dm^2_21
  //(num4)= dm^2_32
  //(num5) = dcp
  bNunue->SetMNS( num0, num2, num1, num3, num4, num5, energy, non_doubled_angle, nutype);
  
  const double distPoint = totaldistance/44000;
  double distance = 0;
  
  for(int ii=0; ii<npoint && distance<=totaldistance; ii++){
    // Set PMNS and propagate with given distance
    bNunue->propagateLinear(nutype, distance, density);
    
    //probablity projection.
    const double t_pnue = bNunue->GetProb( nutype, 1);
    const double t_pnumu = bNunue->GetProb( nutype, 2);
    const double t_pnutau = bNunue->GetProb( nutype, 3);
    //const double t_sum = t_pnue + t_pnumu + t_pnutau;


    const double radangle = -135 * (M_PI/180);
    const double cosa = cos(radangle);
    const double sina = sin(radangle);

    const double tmpX =  cosa * t_pnumu - sina * t_pnutau;
    const double tmpY =  sina * t_pnumu + cosa * t_pnutau + 1/sqrt(2);
    const double tmpZ =  t_pnue;

    const double t_X = tmpX;

    const double radangle2 = -atan(sqrt(2));

    const double cosa2 = cos(radangle2);
    const double sina2 = sin(radangle2);
    
    const double t_Y =  cosa2 * tmpY - sina2 * tmpZ;
    //const double t_Z =  sina2 * tmpY + cosa2 * tmpZ;
    //X.push_back(t_X);
    //Y.push_back(t_Y);
    //D.push_back(distance);
    B[ii]=t_X;
    B[ii+44000]=t_Y;
    B[ii+88000]=distance;
    distance += distPoint;
  };
    
  return B;
}


// Function that averages elements of array B

const int npoint = 88000;

double* average(double* B[]){

  const int n = 1000;
  double* A = (double*)malloc(npoint * 8 * 3);

  for (int j=0; j<npoint-n; j++){

    double sum = 0;
    for (int i=j; i<n+j; i++){
      sum = sum + *B[i];
    };
    A[j] = sum/((double)n);
  };
  return A;
}
	
     
EM_PORT_API(void) free_buf(void* buf) {
  free(buf);
}
