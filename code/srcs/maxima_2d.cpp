//
// maxima_2d code for use with DEBLENDER
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
void maxima_2d(const double*img_beg, IDL_LONG *res_matrix, IDL_LONG *count, IDL_LONG *naxes1, IDL_LONG *naxes2){
   long axes1=*naxes1;
   long axes2=*naxes2;
   long count_area=*count;
   if (count_area==8){         
       for(long i=0;i!=axes1;i++){
         for(long j=0;j!=axes2;j++){
            
            long step=1;
            long ijind=i+j*axes1;
            res_matrix[ijind]=-1;

            long imin=i-step<0?i:i-step;
            long imax=i+step<axes1?i+step:i;
            long jmin=j-step<0?j:j-step;
            long jmax=j+step<axes2?j+step:j;
            
            long ijmin=i+jmin*axes1;
            long ijmax=i+jmax*axes1;
            
            long iminj=imin+j*axes1;
            long iminjmin=imin+jmin*axes1;
            long iminjmax=imin+jmax*axes1;
            
            long imaxj=imax+j*axes1;
            long imaxjmin=imax+jmin*axes1;
            long imaxjmax=imax+jmax*axes1;
            
            if (img_beg[ijmin]<img_beg[ijind] and img_beg[ijmax]<img_beg[ijind] and
                img_beg[iminj]<img_beg[ijind] and img_beg[iminjmin]<img_beg[ijind] and img_beg[iminjmax]<img_beg[ijind] and
                img_beg[imaxj]<img_beg[ijind] and img_beg[imaxjmin]<img_beg[ijind] and img_beg[imaxjmax]<img_beg[ijind]){
                    res_matrix[ijind]=1;
            }
         }
        }
   }
   if (count_area==4){         
       for(long i=0;i!=axes1;i++){
         for(long j=0;j!=axes2;j++){
            
            long step=1;
            long ijind=i+j*axes1;
            res_matrix[ijind]=-1;

            long imin=i-step<0?i:i-step;
            long imax=i+step<axes1?i+step:i;
            long jmin=j-step<0?j:j-step;
            long jmax=j+step<axes2?j+step:j;
            
            long ijmin=i+jmin*axes1;
            long ijmax=i+jmax*axes1;
            
            long iminj=imin+j*axes1;
            long iminjmin=imin+jmin*axes1;
            long iminjmax=imin+jmax*axes1;
            
            long imaxj=imax+j*axes1;
            long imaxjmin=imax+jmin*axes1;
            long imaxjmax=imax+jmax*axes1;
            
            if (count_area==4 and img_beg[ijmin]<img_beg[ijind] and img_beg[ijmax]<img_beg[ijind] and
               	    img_beg[iminj]<img_beg[ijind] and img_beg[imaxj]<img_beg[ijind]){res_matrix[ijind]=1;}
         }
       }
    }
       
}

int main(int argc, void *argv[])
{
  //cout<<"Greetings from c++ maxima_2d code! "<<endl;
  if (argc!=5) return 0;
  maxima_2d((double*)argv[0], (IDL_LONG *)argv[1], (IDL_LONG *)argv[2], (IDL_LONG *)argv[3], (IDL_LONG *)argv[4]);
  //cout<<"c++ maxima_2d code finished its job! "<<endl;
  //cout<<endl;
  return 1;
}
