//
// flood shedding code for use with DEBLENDER
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

void flood_shed(const double *image, IDL_LONG *res, IDL_LONG *maxima_ind, IDL_LONG *nele_maxima, double *count_maxima_ind, IDL_LONG *naxes1, IDL_LONG *naxes2){
  IDL_LONG nele=*nele_maxima; 
  IDL_LONG axes1=*naxes1;
  IDL_LONG axes2=*naxes2;
  IDL_LONG nele_all=axes1*axes2; 
  for (long i=0;i<nele_all;i++) res[i]=-1000;
  for (long i=0;i<nele;i++) res[maxima_ind[i]]=i;
  //cout<<axes1<<endl;
  //cout<<axes2<<endl;

  multimap<double, long> res_map;
  for (long i=0;i<nele_all;i++) { if(res[i]<-1) {res_map.insert (pair<double, long>(-image[i], i));}}
  //cout<<"2"<<endl;
  

  multimap<double, long>::iterator iter=res_map.begin();
  while(iter!=res_map.end()){
       long k=(*iter).second;
       long i=k%axes1;
       long j=k/axes1;
       //cout<<k<<endl;

       long imin=(i<1)?i:i-1;
       long imax=(i>=axes1-1)?i:i+1;
       long jmin=(j<1)?j:j-1;
       long jmax=(j>=axes2-1)?j:j+1;
       long ngb_indices[8]={imin+jmin*axes1, imin+j*axes1, imin+jmax*axes1,  i+jmin*axes1, i+jmax*axes1,  imax+jmin*axes1, imax+j*axes1, imax+jmax*axes1};
       IDL_LONG res1=-1000, res2=-1000;
       for (long ngb_i=0;ngb_i !=8;ngb_i++){
           long ngb_ind=ngb_indices[ngb_i];
           if (res[ngb_ind]>-1){res2=res[ngb_ind];}
           if (res1<-1){res1=res2;}
       }
       if (res1==res2 and res1>-1 ){
          res[k]=res1;
          count_maxima_ind[res1]=count_maxima_ind[res1]+1.0;
       }

       res_map.erase(iter);
       iter=res_map.begin();
  }
}  
  

int main(int argc, void *argv[])
{
  //cout<<"Greetings from c++ water shed segmentation code! "<<endl;
  if (argc!=7) return 0;
  flood_shed((double*)argv[0], (IDL_LONG *)argv[1], (IDL_LONG *)argv[2], (IDL_LONG *)argv[3], (double *)argv[4], (IDL_LONG *)argv[5], (IDL_LONG *)argv[6]);
  //cout<<"c++ water shed code finished its job! "<<endl;
  //cout<<endl;
  return 1;
}


