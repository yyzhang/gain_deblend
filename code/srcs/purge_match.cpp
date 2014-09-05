//
// purge match code for use with DEBLENDER
//
//  Created by Yuan-Yuan Zhang on 10/23/13.
//
//
#include <algorithm>
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "idl_export.h"
#include <map>
#include <utility>
using namespace std;

#define READONLY 0
#define READWRITE 1
void purge_match(const double*image, IDL_LONG *external_indices, IDL_LONG *external_nele, IDL_LONG *maxima_ind, IDL_LONG *maxima_nele, IDL_LONG *naxes1, IDL_LONG *naxes2, const double *dis){
     long axes1=*naxes1;
     long axes2=*naxes2;
     long nele_all=axes1*axes2; 
     long nex=*external_nele;
     long nma=*maxima_nele;  
     double match_dis=*dis;
     
     for (long i=0;i!=nex;i++){
         double min_dis=100000;
         double dis_temp=10000000;
         long j0=0;

         long ex_x=external_indices[i]%axes1;
         long ex_y=external_indices[i]/axes1;
         for (long j=0;j!=nma;j++){
             if (maxima_ind[j]>-0.5){
                IDL_LONG xind=maxima_ind[j]%axes1;
                IDL_LONG yind=maxima_ind[j]/axes1;
                dis_temp=sqrt((ex_x-xind)*(ex_x-xind)+(ex_y-yind)*(ex_y-yind));
                if (dis_temp<min_dis){ 
                   min_dis=dis_temp;
                   j0=j;
                }
              }
         }
         if (min_dis<match_dis) {maxima_ind[j0]=-1;}
     }

}
int main(int argc, void *argv[])
{
  //cout<<"Greetings from c++ purge code! "<<endl;
  if (argc!=8) return 0;
  purge_match((double*)argv[0], (IDL_LONG *)argv[1], (IDL_LONG *)argv[2], (IDL_LONG *)argv[3], (IDL_LONG *)argv[4], (IDL_LONG *)argv[5], (IDL_LONG *)argv[6], (double  *)argv[7]);
  return 1;
}
