//
// fmm code for use with DEBLENDER
//
//  Created by Yuan-Yuan Zhang on 10/18/13.
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

void fmm(const double *img_beg, double *inpainted_image, double *inpaintf, IDL_LONG *naxes1, IDL_LONG *naxes2, const double *gradx, const double *grady);
double inpaint_value(const double *inpainted_image, const double *inpaint_f, const double *inpaint_t, long epsilo, long ind, long axes1, long axes2, const double *gradx, const double *grady);
void fmm2(const double *img_beg, double *inpainted_image, double *inpaintf, IDL_LONG *naxes1, IDL_LONG *naxes2, const double *gradx, const double *grady);
double inpaint_value2(const double *inpainted_image, const double *inpaint_f, const double *inpaint_t, long epsilo, long ind, long axes1, long axes2, const double *gradx, const double *grady);
void grad(const double *img_hd, const double *inpaint_f, long ind, long axes1, long axes2, double *gradi);
double solve_t(const double *img_t, const double *inpaint_f, long i1, long j1, long i2, long j2, long axes1, long axes2);

int main(int argc, void *argv[])
{
  //cout<<"Greetings from c++! "<<endl;
  if (argc!=8) return 0;
  const double *img_beg=(double*)argv[0];
  double *inpainted_image=(double *)argv[1];
  double *inpaint_radius=(double *)argv[2];
  IDL_LONG *axes1=(IDL_LONG *)argv[3];
  IDL_LONG *axes2=(IDL_LONG *)argv[4];
  IDL_LONG *flag=(IDL_LONG *)argv[5];
  const double *gradx=(double*)argv[6];
  const double *grady=(double*)argv[7];

  IDL_LONG xc=((*axes1)-(*axes1)%2)/2;
  IDL_LONG yc=((*axes2)-(*axes2)%2)/2;
  if ((*inpaint_radius)>sqrt(xc*yc)){*inpaint_radius=sqrt(xc*yc);}
  IDL_LONG nele=(*axes1)*(*axes2);
  double *inpaint_f=new double[nele];
  double *inpainted_image1=new double[nele];
  double *inpainted_image2=new double[nele];
  for(long i=0;i!=nele;i++){
     IDL_LONG ix=i%(*axes1);
     IDL_LONG iy=(i-i%(*axes1))/(*axes1);
     if (sqrt((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc))<(*inpaint_radius)){inpaint_f[i]=1;}
     else{inpaint_f[i]=0;}
     inpainted_image1[i]=inpainted_image[i];
     inpainted_image2[i]=inpainted_image[i];
  }
  fmm(img_beg, inpainted_image1, inpaint_f, axes1, axes2, gradx, grady);
  for(long i=0;i!=nele;i++){
     IDL_LONG ix=i%(*axes1);
     IDL_LONG iy=(i-i%(*axes1))/(*axes1);
     if (sqrt((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc))<(*inpaint_radius)){inpaint_f[i]=1;}
     else{inpaint_f[i]=0;}
  }
  fmm2(img_beg, inpainted_image2, inpaint_f, axes1, axes2, gradx, grady);
  delete[] inpaint_f;
  
  double sum_b1=0.0;
  double npixels_b1=0.0;
  double aver_b1=0.0;
  double sum_b2=0.0;
  double npixels_b2=0.0;
  double aver_b2=0.0;
  
  for(long i=0;i!=nele;i++){
     IDL_LONG ix=i%(*axes1);
     IDL_LONG iy=(i-i%(*axes1))/(*axes1);
    if (sqrt((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc))>=(*inpaint_radius)-2.0 and sqrt((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc))<(*inpaint_radius)){
       sum_b1=sum_b1+img_beg[i]-inpainted_image1[i];
       npixels_b1=npixels_b1+1;
       sum_b2=sum_b2+img_beg[i]-inpainted_image2[i];
       npixels_b2=npixels_b2+1;
    }
  }


  if (npixels_b1>0 and npixels_b2>0){
     aver_b1=sum_b1/npixels_b1;
     aver_b2=sum_b2/npixels_b2;  
  }
  for(long i=0;i!=nele;i++){
        inpainted_image[i]=inpainted_image1[i];
        IDL_LONG ix=i%(*axes1);
        IDL_LONG iy=(i-i%(*axes1))/(*axes1);
        if (sqrt((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc))<(*inpaint_radius)){
           inpainted_image2[i]=inpainted_image2[i]+aver_b2;
           //;
           if (img_beg[i]<inpainted_image[i]){inpainted_image[i]=img_beg[i];}
           inpainted_image[i]=(inpainted_image[i]*0.7+inpainted_image2[i]*0.3);
        }
   }


  delete[] inpainted_image1;
  delete[] inpainted_image2;

  return 1;
}

void fmm(const double *img_beg, double *inpainted_image, double *inpaint_f, IDL_LONG *naxes1, IDL_LONG *naxes2, const double *gradx, const double *grady){
	//inpaint_f. 1: to be inpainted. 0: area known. -1 narrowband
        //cout<<img_beg[150+2*301]<<endl;
        //initializing inpaint_t
        double *inpaint_t;
        long axes1=(*naxes1);
        long axes2=(*naxes2);
        //cout<<axes1<<"  "<<axes2<<endl;
        IDL_LONG nele=axes1*axes2;
        inpaint_t=new double[nele];
        for (long k=0; k<nele; k++){ 
            inpaint_t[k]=img_beg[k];
            (inpaint_f[k]<0.5)?inpaint_t[k]=0:inpaint_t[k]=1000000;
        }
        //cout<<"Greetings from FMM"<<endl;
	//initializing narrowband and inpaint_f
	multimap<double, long> narrowband;
        for (long k=0; k<nele; k++){
             if (inpaint_f[k]<0.5){
             	long i=k%axes1;
    	     	long j=k/axes1;
   	     	long imin=(i<1)?k:k-1;
             	long imax=(i==axes1-1)?k:k+1;
             	long jmin=(j<1)?k:k-axes1;
             	long jmax=(j==axes2-1)?k:k+axes1;
             	if (inpaint_f[imin]>0.5 or inpaint_f[imax]>0.5 or inpaint_f[jmin]>0.5 or inpaint_f[jmax]>0.5){
                    inpaint_f[k]=-1; //-1 narrowband flag value
                    narrowband.insert (pair<double, long>(inpaint_t[k], k));
                }  
             }
        }

        multimap<double, long>::iterator iter=narrowband.begin();
        while(iter!=narrowband.end()){
             long k=(*iter).second;
             inpaint_f[k]=0; //change flag to known
             long i=k%axes1;
    	     long j=k/axes1;
   	     long imin=(i<1)?k:k-1;
             long imax=(i==axes1-1)?k:k+1;
             long jmin=(j<1)?k:k-axes1;
             long jmax=(j==axes2-1)?k:k+axes1;
             long ngb_indices[4]={imin, imax, jmin, jmax};
             for(int i_ind=0;i_ind<4;i_ind++){
                long ngb_ind=ngb_indices[i_ind];
                if(inpaint_f[ngb_ind] > 0.5){
                  long epsilo2=10;
                  inpainted_image[ngb_ind]=inpaint_value(inpainted_image, inpaint_f, inpaint_t, epsilo2, ngb_ind, axes1, axes2, gradx, grady);//inpaint this pixel
                  inpaint_f[ngb_ind]=-1;//change flag value to NarrowBand
                  long ingb=ngb_ind%axes1, jngb=ngb_ind/axes1;
                  double t1, t2, t3, t4;
                  t1=solve_t(inpaint_t, inpaint_f, ingb-1, jngb, ingb, jngb-1, axes1, axes2);
                  t2=solve_t(inpaint_t, inpaint_f, ingb+1, jngb, ingb, jngb-1, axes1, axes2);
                  t3=solve_t(inpaint_t, inpaint_f, ingb-1, jngb, ingb, jngb-1, axes1, axes2);
                  t4=solve_t(inpaint_t, inpaint_f, ingb+1, jngb, ingb, jngb+1, axes1, axes2);
                  double tarray[4]={t1, t2, t3, t4};
                  inpaint_t[ngb_ind]=(*min_element(tarray, tarray+3));//Update T value */;
                  narrowband.insert (pair<double, long>(inpaint_t[ngb_ind], ngb_ind));
                }
             }
             narrowband.erase(iter);
             iter=narrowband.begin();
             //cout<<"narrow band size: "<<narrowband.size()<<endl;  
        }
        delete [] inpaint_t;
}

double inpaint_value(const double *inpainted_image, const double *inpaint_f, const double *inpaint_t, long epsilo, long ind, long axes1, long axes2, const double *gradx, const double *grady){
       double Iij=-1;
       double Ia=0, s=0;

       long i=ind%axes1, j=ind/axes1;
       for (long k=i-epsilo; k<=i+epsilo;k++){
           for (long l=j-epsilo; l<=j+epsilo;l++){
              if (k>=0 and k<=axes1-1 and l>=0 and l<=axes2-1 and k+l*axes1<axes1*axes2 and (k-i)*(k-i)+(l-j)*(l-j)<=epsilo*epsilo){
                 long kl_ind=k+l*axes1;
                 long rvector[2]={i-k, j-l};
                 double r_len2=(rvector[0])*(rvector[0])+(rvector[1])*(rvector[1]);
                 if (r_len2 >0.1 and inpaint_f[kl_ind] < 0.5 and inpaint_f[kl_ind] > -0.5 ){
                    double dir=fabs( (gradx[ind]*rvector[0]+grady[ind]*rvector[1])/sqrt(r_len2) );
                    double dst=(1.0/r_len2);
                    double lev=1.0/(1.0+fabs(inpaint_t[kl_ind]-inpaint_t[ind]));
                    double w=dst*dir*lev;//lev*dir;// not using dst. Expect lev to behave similarly to dst.
                    //double gradi[2]={0,0};
                    //grad(inpainted_image, inpaint_f, kl_ind, axes1, axes2, gradi);
                    Ia=Ia+w*(inpainted_image[kl_ind]);//+gradi[0]*rvector[0]+gradi[1]*rvector[1]); //not using gradient. plus always makes interpolated gradient larger.
                    s=s+w;
                    //cout<<"Ia: "<<Ia<<" s:"<< s<<"w: "<<w<<endl;
                 }
              }
           }
       }
       if (s>0){Iij=Ia/s;}
       //cout<<"check image inpainted1:"<<Ia<<"  "<<s<<endl;
       return Iij;
}
void grad(const double *img_hd, const double *inpaint_f, long ind, long axes1, long axes2, double *gradi){
     long i=ind%axes1, j=ind/axes1;
     double vecx=0, vecy=0;
     
     long xlo=(i>0)?i-1:i;
     long xhi=(i<axes1-1)?i+1:i;
     long ylo=(j>0)?j-1:j;
     long yhi=(j<axes2-1)?j+1:j;
     //xlo=(inpaint_f[xlo*axes1+j]>0.5)?i:xlo;
     //xhi=(inpaint_f[xhi*axes1+j]>0.5)?i:xhi;
     //ylo=(inpaint_f[i*axes1+ylo]>0.5)?j:ylo;
     //yhi=(inpaint_f[i*axes1+yhi]>0.5)?j:yhi;
     if (inpaint_f[xlo+j*axes1]>0.5 or inpaint_f[xhi+j*axes1]>0.5 or inpaint_f[i+ylo*axes1]>0.5 or inpaint_f[i+yhi*axes1]>0.5){
        xlo=i;
        xhi=i;
        ylo=j;
        yhi=j;
     }

     vecx=img_hd[xhi+j*axes1]-img_hd[xlo+j*axes1];
     vecy=img_hd[i+yhi*axes1]-img_hd[i+ylo*axes1];
     //cout<<"before wweighting:"<<vecx<<"  "<<vecy<<endl;
     vecx=(fabs(xlo-xhi)>0.5)?vecx/double(xhi-xlo):0.0;
     vecy=(fabs(ylo-yhi)>0.5)?vecy/double(yhi-ylo):0.0;
     gradi[0]=vecx;
     gradi[1]=vecy;
     //cout<<"after weighting:"<<vecx<<"  "<<vecy<<endl;
}

double solve_t(const double *img_t, const double *inpaint_f, long i1, long j1, long i2, long j2, long axes1, long axes2){
       double sol=1000000;
       long ind1=i1+j1*axes1, ind2=i2+j2*axes1;
       double f2=1;

       if (i1>=0 and i1<=axes1-1 and j1>=0 and j1<=axes2-1){
           if (inpaint_f[ind1]<0.5){
               sol=f2+img_t[ind1];
               if (i2>=0 and i2<=axes1-1 and j2>=0 and j2<=axes2-1){
                  if (inpaint_f[ind2]<0.5){
                     double r=2.0*f2-(img_t[ind1]-img_t[ind2])*(img_t[ind1]-img_t[ind2]);
                     if (r>=0){
                     	double s=(img_t[ind1]+img_t[ind2]+sqrt(r))/2.0;
                        if (s >=img_t[ind1] and s>=img_t[ind2]) {sol=s;}
                     }
                  }
               }
           }
           else if (i2>=0 and i2<=axes1-1 and j2>=0 and j2<=axes2-1){
                if (inpaint_f[ind2]<0.5){sol=f2+img_t[ind2];}      
           }  
       }
       else if (i2>=0 and i2<=axes1-1 and j2>=0 and j2<=axes2-1){
            if (inpaint_f[ind2]<0.5){sol=f2+img_t[ind2];}      
       }
       
       return sol;
}

void fmm2(const double *img_beg, double *inpainted_image, double *inpaint_f, IDL_LONG *naxes1, IDL_LONG *naxes2, const double *gradx, const double *grady){

        double *inpaint_t;
        long axes1=(*naxes1);
        long axes2=(*naxes2);
        IDL_LONG nele=axes1*axes2;
        inpaint_t=new double[nele];
        for (long k=0; k<nele; k++){inpaint_t[k]=img_beg[k];}

	multimap<double, long> narrowband;
        for (long k=0; k<nele; k++){
             if (inpaint_f[k]<0.5){
             	long i=k%axes1;
    	     	long j=k/axes1;
   	     	long imin=(i<1)?k:k-1;
             	long imax=(i==axes1-1)?k:k+1;
             	long jmin=(j<1)?k:k-axes1;
             	long jmax=(j==axes2-1)?k:k+axes1;
             	if (inpaint_f[imin]>0.5 or inpaint_f[imax]>0.5 or inpaint_f[jmin]>0.5 or inpaint_f[jmax]>0.5){
                    inpaint_f[k]=-1; //-1 narrowband flag value
                    narrowband.insert (pair<double, long>(inpaint_t[k], k));
                }  
             }
        }

        multimap<double, long>::iterator iter=narrowband.begin();
        while(iter!=narrowband.end()){
             long k=(*iter).second;
             inpaint_f[k]=0; //change flag to known
             long i=k%axes1;
    	     long j=k/axes1;
   	     long imin=(i<1)?k:k-1;
             long imax=(i==axes1-1)?k:k+1;
             long jmin=(j<1)?k:k-axes1;
             long jmax=(j==axes2-1)?k:k+axes1;
             long ngb_indices[4]={imin, imax, jmin, jmax};
             for(int i_ind=0;i_ind<4;i_ind++){
                long ngb_ind=ngb_indices[i_ind];
                if(inpaint_f[ngb_ind] > 0.5){
                  long epsilo2=10;
                  inpainted_image[ngb_ind]=inpaint_value2(inpainted_image, inpaint_f, inpaint_t, epsilo2, ngb_ind, axes1, axes2, gradx, grady);//inpaint this pixel
                  inpaint_f[ngb_ind]=-1;//change flag value to NarrowBand
                  narrowband.insert (pair<double, long>(inpaint_t[ngb_ind], ngb_ind));
                }
             }
             narrowband.erase(iter);
             iter=narrowband.begin();
             //cout<<"narrow band size: "<<narrowband.size()<<endl;  
        }
        delete [] inpaint_t;
}

double inpaint_value2(const double *inpainted_image, const double *inpaint_f, const double *inpaint_t, long epsilo, long ind, long axes1, long axes2, const double *gradx, const double *grady){
       double Iij=-1;
       double Ia=0, s=0;

       long i=ind%axes1, j=ind/axes1;
       for (long k=i-epsilo; k<=i+epsilo;k++){
           for (long l=j-epsilo; l<=j+epsilo;l++){
              if (k>=0 and k<=axes1-1 and l>=0 and l<=axes2-1 and k+l*axes1<axes1*axes2 and (k-i)*(k-i)+(l-j)*(l-j)<=epsilo*epsilo){
                 long kl_ind=k+l*axes1;
                 long rvector[2]={i-k, j-l};
                 double r_len2=(rvector[0])*(rvector[0])+(rvector[1])*(rvector[1]);
                 if (r_len2 >0.1 and inpaint_f[kl_ind] < 0.5 and inpaint_f[kl_ind] > -0.5 ){
                    double dir=fabs( (gradx[ind]*rvector[0]+grady[ind]*rvector[1])/sqrt(r_len2) );
                    double dst=(1.0/r_len2);
                    double lev=1.0/(1.0+fabs(inpaint_t[kl_ind]-inpaint_t[ind]));
                    double w=dst*dir*lev;//lev*dir;// not using dst. Expect lev to behave similarly to dst.
                    //double gradi[2]={0,0};
                    //grad(inpainted_image, inpaint_f, kl_ind, axes1, axes2, gradi);
                    Ia=Ia+w*(inpainted_image[kl_ind]);
                    s=s+w;
                 }
              }
           }
       }
       if (s>0){Iij=Ia/s;}
       //cout<<"check image inpainted1:"<<Ia<<"  "<<s<<endl;
       return Iij;
}


