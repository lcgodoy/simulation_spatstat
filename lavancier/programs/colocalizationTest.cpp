/*
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License version 3 as published by the
* Free Software Foundation:
* <http://www.gnu.org/licenses/agpl-3.0.txt>.
*  
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more
* details.
* 
* Authors: Thierry Pecot
*/

#include "CImg.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include "math.h"
#include <sstream>
#include <algorithm>

using namespace cimg_library;
using namespace std;

static double CND(double d)
{
    const double       A1 = 0.31938153;
    const double       A2 = -0.356563782;
    const double       A3 = 1.781477937;
    const double       A4 = -1.821255978;
    const double       A5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;
    
    double
    K = 1.0 / (1.0 + 0.2316419 * fabs(d));
    
    double
    cnd = RSQRT2PI * exp(- 0.5 * d * d) *
    (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
    
    if (d > 0)
    cnd = 1.0 - cnd;
    
    return cnd;
}

template<typename T>
CImgList<T> complex_mul(const CImgList<T> &z1, const CImgList<T> &z2){
    if (!z1.is_sameXYZC(z2)) throw CImgArgumentException("complex_mul inputs have not the same size.");
    CImgList<> dest(z1);
#pragma omp parallel for
    for(unsigned long i = 0; i < dest[0].size(); i++) {
        const T a = z1[0](i), b = z1[1](i), c = z2[0](i), d = z2[1](i);
        dest[0](i) = a * c - b * d;
        dest[1](i) = b * c + a * d;
    }
    return dest;
}

template<typename T>
CImgList<T> complex_conjugate(const CImgList<T> &z){
    CImgList<> dest(2);
    dest[0] = z[0];
    dest[1] = -z[1];
    return dest;
}

//! compute a sum along a dimension of the image
/**
 \param img an input image
 \param axis can be 'x','y','z','c'
 **/
template<typename T>
CImg<T> sum(const CImg<T> & img, const char axis='y'){
    CImg<T> dest;
    if (axis == 'x') {
        dest.assign(1,img._height,img._depth,img._spectrum,0);
        cimg_forXYZC(img,x,y,z,c){ dest(0,y,z,c) += img(x,y,z,c);}
    } else if (axis == 'y'){
        dest.assign(img._width,1,img._depth,img._spectrum,0);
        cimg_forXYZC(img,x,y,z,c){ dest(x,0,z,c) += img(x,y,z,c);}
    } else if (axis == 'z'){
        dest.assign(img._width,img._height,1,img._spectrum,0);
        cimg_forXYZC(img,x,y,z,c){ dest(x,y,0,c) += img(x,y,z,c);}
    } else if (axis == 'c'){
        dest.assign(img._width,img._height,img._depth,1,0);
        cimg_forXYZC(img,x,y,z,c){ dest(x,y,0,c) += img(x,y,z,c);}
    } else
    throw CImgArgumentException("sum(const CImg<>, const char): invalid axis name.");
    return dest;
}

//! compute the mean along a dimension of the image
/**
 \param img an input image
 \param axis can be 'x','y','z','c'
 **/
template<typename T>
CImg<T> mean(const CImg<T> & img, const char axis='y'){
    if (axis == 'x') return sum(img, axis) / img.width();
    else if (axis == 'y') return sum(img, axis) / img.height();
    else if (axis == 'z') return sum(img, axis) / img.depth();
    else if (axis == 'c') return sum(img, axis) / img.spectrum();
    else
    throw CImgArgumentException("mean(const CImg<>, const char): invalid axis name.");
}

//! Compute the covariance matrix of x as a collection of vectors
/**
 \param img a 2D image
 \return ((img-m)*(img-m).get_transpose())/(w.width()-1)
 **/
template<typename T>
CImg<> covariance(const CImg<T> &img){
    if (img.depth()!=1 && img.spectrum()!=1)
    throw CImgArgumentException("Covariance(const CImg<T> &): argument is not a 2D image (%dx%dx%dx%d).",
                                img.width(),img.height(),img.depth(),img.spectrum());
    CImg<> dest = img, m = mean(img,'x');
    cimg_forXY(dest,x,y) dest(x,y) -= m(y);
    dest = dest * dest.get_transpose();
    return dest/(dest.width()-1);
}

//! Shift image for centering the 0 of the Fourier transform
/**
 \param image
 \return the shifted imageusing periodic boundaries
 **/
template<typename T>
CImg<> fftshift(const CImg<T> &img) {
    return img.get_shift(-img.width()/2, -img.height()/2, -img.depth()/2, 0, 2);
}

//! Compute the convolution between two images using the Fast Fourier Transform
/**
 \param img0 : first image
 \param img1 : second image
 \note Computed as \f$ \mathrm{Re}(\mathcal{F}^{-1}(\mathcal{F}(f)\mathcal{F}(g))) \f$
 **/
template<typename T>
CImg<> convolve_fft(const CImg<T> &f, const CImg<T>&g){
    if (!f.is_sameXYZC(g)) throw CImgArgumentException("convolve_fft inputs have not the same size.");
    return fftshift((complex_mul(f.get_FFT(), g.get_FFT()).FFT(true))[0]);
}

//! Compute the correlation between two images using the Fast Fourier Transform
/**
 \param f : first image
 \param g : second image
 \note Computed as \f$  \mathrm{Re}(\mathcal{F}^{-1}(\mathcal{F}(f)\mathcal{F}(g)^{*})) \f$
 **/
template<typename T>
CImg<> correlate_fft(const CImg<T> & f, const CImg<T> &g){
    if (!f.is_sameXYZC(g)) throw CImgArgumentException("correlate_fft inputs have not the same size.");
   f.get_FFT().display();
   return fftshift((complex_mul(f.get_FFT(), complex_conjugate(g.get_FFT())).FFT(true))[0]);
}

! Compute the autocorrelation of the image using the Fast Fourier Transform
**
\param img : input image
\note Computed as \f$ \mathrm{Re}(\mathcal{F}^{-1}(\mathcal{F}(f)\mathcal{F}(f)^{*})) \f$
**/
template<typename T>
CImg<T> autocorrelate_fft(const CImg<T> &img){
   return correlate_fft(img, img);
}

template<typename T>
CImg<T> autocorrelation(const CImg<T> &img,int w,int h,int d){
   
   CImg<T> output=autocorrelate_fft(img)/(w*h*d);
   
   if(d==1){
       output.crop((output.width()-2*w)/2,(output.height()-2*h)/2,
                   output.width()-(output.width()-2*w)/2,output.height()-(output.height()-2*h)/2);
   }
   else{
       output.crop((output.width()-2*w)/2,(output.height()-2*h)/2,(output.depth()-2*d)/2,
                   output.width()-(output.width()-2*w)/2,output.height()-(output.height()-2*h)/2,output.depth()-(output.depth()-2*d)/2);
   }
   return output;
   
}

void loadSequences(const char * i1,const char * i2,const bool tf,int w,int h,int d,int firstFrame,int lastFrame,int nbPlanes,int planeDivider,int maxImages,int channel1,int channel2,CImgList<> & input1,CImgList<> & input2){
   
   if(i2==NULL){
       if(tf){
           for(int t=firstFrame;t<=lastFrame;t++){
               char seq_filename[1024];
               sprintf(seq_filename,i1,t);
               ifstream file(seq_filename);
               CImg<> sequence1(w,h,d,1,0),sequence2(w,h,d,1,0);
               for(int z=0;z<d;z++){
                   for(int y=0;y<h;y++){
                       for(int x=0;x<w;x++){
                           if(planeDivider>-1){
                               if(z%planeDivider==0){
                                   file >> sequence1(x,y,z) >> sequence2(x,y,z);
                               }
                           }
                           else{
                               file >> sequence1(x,y,z) >> sequence2(x,y,z);
                           }
                       }
                   }
               }
               input1.push_back(sequence1);
               input2.push_back(sequence2);
           }
       }
       else{
           CImg<> sequence1(i1);
           if(sequence1.spectrum()>1){
               if(nbPlanes>1){
                   CImg<> inputTmp1(sequence1.width(),sequence1.height(),nbPlanes,1,0),inputTmp2(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                   int planeCpt=0,totalImages=0;
                   if((maxImages>0)&&(maxImages<(sequence1.depth()/nbPlanes))){
                       totalImages = maxImages*nbPlanes;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i++){
                       if(channel1>-1){
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y,i%nbPlanes) = sequence1(x,y,i,channel1);
                               inputTmp2(x,y,i%nbPlanes) = sequence1(x,y,i,channel2);
                           }
                       }
                       else{
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y,i%nbPlanes) = sequence1(x,y,i,0);
                               inputTmp2(x,y,i%nbPlanes) = sequence1(x,y,i,1);
                           }
                       }
                       planeCpt++;
                       if(planeCpt==nbPlanes){
                           inputTmp1.threshold(1);
                           inputTmp2.threshold(1);
                           input1.push_back(inputTmp1);
                           input2.push_back(inputTmp2);
                           planeCpt=0;
                           inputTmp1.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                           inputTmp2.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                       }
                   }
               }
               else{
                   int totalImages=0;
                   if((maxImages>0)&&(maxImages<sequence1.depth())){
                       totalImages = maxImages;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i++){
                       CImg<> inputTmp1(sequence1.width(),sequence1.height(),1,1,0),inputTmp2(sequence1.width(),sequence1.height(),1,1,0);
                       if(channel1>-1){
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y) = sequence1(x,y,i,channel1);
                               inputTmp2(x,y) = sequence1(x,y,i,channel2);
                           }
                       }
                       else{
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y) = sequence1(x,y,i,0);
                               inputTmp2(x,y) = sequence1(x,y,i,1);
                           }
                       }
                       inputTmp1.threshold(1);
                       inputTmp2.threshold(1);
                       input1.push_back(inputTmp1);
                       input2.push_back(inputTmp2);
                   }
               }
           }
           else{
               if(nbPlanes>1){
                   CImg<> inputTmp1(sequence1.width(),sequence1.height(),nbPlanes,1,0),inputTmp2(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                   int planeCpt=0,totalImages=0;
                   if((maxImages>0)&&(maxImages<(sequence1.depth()/nbPlanes))){
                       totalImages = maxImages*nbPlanes;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i++){
                       if(channel1>-1){
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y,i%nbPlanes) = sequence1(x,y,i,channel1);
                               inputTmp2(x,y,i%nbPlanes) = sequence1(x,y,i,channel2);
                           }
                       }
                       else{
                           cimg_forXY(sequence1,x,y){
                               inputTmp1(x,y,i%nbPlanes) = sequence1(x,y,i,0);
                               inputTmp2(x,y,i%nbPlanes) = sequence1(x,y,i,1);
                           }
                       }
                       planeCpt++;
                       if(planeCpt==nbPlanes){
                           inputTmp1.threshold(1);
                           inputTmp2.threshold(1);
                           input1.push_back(inputTmp1);
                           input2.push_back(inputTmp2);
                           planeCpt=0;
                           inputTmp1.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                           inputTmp2.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                       }
                   }
               }
               else{
                   int totalImages=0;
                   if((maxImages>0)&&(maxImages<(sequence1.depth()/2))){
                       totalImages = 2*maxImages;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i+=2){
                       CImg<> inputTmp1(sequence1.width(),sequence1.height(),1,1,0),inputTmp2(sequence1.width(),sequence1.height(),1,1,0);
                       cimg_forXY(sequence1,x,y){
                           inputTmp1(x,y) = sequence1(x,y,i);
                           inputTmp2(x,y) = sequence1(x,y,i+1);
                       }
                       inputTmp1.threshold(1);
                       inputTmp2.threshold(1);
                       input1.push_back(inputTmp1);
                       input2.push_back(inputTmp2);
                   }
               }
           }
       }
   }
   else{
       if(tf){
           for(int t=firstFrame;t<=lastFrame;t++){
               char seq_filename1[1024],seq_filename2[1024];
               sprintf(seq_filename1,i1,t);
               ifstream file_1(seq_filename1);
               sprintf(seq_filename2,i2,t);
               ifstream file_2(seq_filename2);
               if(planeDivider>-1){
                   CImg<> sequence1(w,h,d/planeDivider,1,0),sequence2(w,h,d/planeDivider,1,0);
                   for(int z=0;z<d;z++){
                       if(z%planeDivider==0){
                           for(int y=0;y<h;y++){
                               for(int x=0;x<w;x++){
                                   file_1 >> sequence1(x,y,z/planeDivider);
                                   file_2 >> sequence2(x,y,z/planeDivider);
                               }
                           }
                       }
                       else{
                           int tmp;
                           for(int y=0;y<h;y++){
                               for(int x=0;x<w;x++){
                                   file_1 >> tmp;
                                   file_2 >> tmp;
                               }
                           }
                       }
                   }
                   sequence1.save("image1.tif");sequence2.save("image2.tif");
                   input1.push_back(sequence1);
                   input2.push_back(sequence2);
               }
               else{
                   CImg<> sequence1(w,h,d,1,0),sequence2(w,h,d,1,0);
                   for(int z=0;z<d;z++){
                       for(int y=0;y<h;y++){
                           for(int x=0;x<w;x++){
                               file_1 >> sequence1(x,y,z);
                               file_2 >> sequence2(x,y,z);
                           }
                       }
                   }
                   input1.push_back(sequence1);
                   input2.push_back(sequence2);
               }
           }
       }
       else{
           CImg<> sequence1(i1),sequence2(i2);
           if((sequence1.width()==sequence2.width())&&(sequence1.height()==sequence2.height())&&(sequence1.depth()==sequence2.depth())){
               if(nbPlanes>1){
                   CImg<> inputTmp1(sequence1.width(),sequence1.height(),nbPlanes,1,0),inputTmp2(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                   int planeCpt=0,totalImages=0;
                   if((maxImages>0)&&(maxImages<(sequence1.depth()/nbPlanes))){
                       totalImages = maxImages*nbPlanes;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i++){
                       cimg_forXY(sequence1,x,y){
                           inputTmp1(x,y,i%nbPlanes) = sequence1(x,y,i);
                           inputTmp2(x,y,i%nbPlanes) = sequence2(x,y,i);
                       }
                       planeCpt++;
                       if(planeCpt==nbPlanes){
                           inputTmp1.threshold(1);
                           inputTmp2.threshold(1);
                           input1.push_back(inputTmp1);
                           input2.push_back(inputTmp2);
                           planeCpt=0;
                           inputTmp1.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                           inputTmp2.assign(sequence1.width(),sequence1.height(),nbPlanes,1,0);
                       }
                   }
               }
               else{
                   int totalImages=0;
                   if((maxImages>0)&&(maxImages<sequence1.depth())){
                       totalImages = maxImages;
                   }
                   else{
                       totalImages = sequence1.depth();
                   }
                   for(int i=0;i<totalImages;i++){
                       CImg<> inputTmp1(sequence1.width(),sequence1.height(),1,1,0),inputTmp2(sequence1.width(),sequence1.height(),1,1,0);
                       cimg_forXY(sequence1,x,y){
                           inputTmp1(x,y) = sequence1(x,y,i);
                           inputTmp2(x,y) = sequence2(x,y,i);
                       }
                       inputTmp1.threshold(1);
                       inputTmp2.threshold(1);
                       input1.push_back(inputTmp1);
                       input2.push_back(inputTmp2);
                   }
               }
           }
       }
   }
}

void radiusAndDepthTruncations(CImg<> & autocorrelation1,CImg<> & autocorrelation2,float truncationValue,int w,int h,int d,int & radiusForTruncation,int & depthForTruncation){
   
   if(d==1){
       double max1=autocorrelation1.max(),max2=autocorrelation2.max(),distanceMax1=0.,distanceMax2=0.;
       cimg_forXY(autocorrelation1,x,y){
           if(autocorrelation1(x,y)>(truncationValue*max1)){
               if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.))>distanceMax1){
                   distanceMax1 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
               }
           }
           if(autocorrelation2(x,y)>(truncationValue*max2)){
               if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.))>distanceMax2){
                   distanceMax2 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
               }
           }
       }
       int maxRadius=0;
       if(w>h){maxRadius=h;}
       else{maxRadius=w;}
       if(distanceMax1>distanceMax2){radiusForTruncation = cimg::min(ceil(distanceMax1),maxRadius);}
       else{radiusForTruncation = cimg::min(ceil(distanceMax2),maxRadius);}
   }
   else{
   	//autocorrelation1.display();autocorrelation2.display();
       double max1=autocorrelation1.max(),max2=autocorrelation2.max(),distanceMax1=0.,distanceMax2=0.,distanceMax2D1=0.,distanceMax2D2=0.,distanceMaxZ1=0.,distanceMaxZ2=0.;
       cout << "max1: " << max1 << "   max2: " << max2 << endl;
       //cout << max1 << ' ' << max2 << endl;
       cimg_forXYZ(autocorrelation1,x,y,z){
           if(autocorrelation1(x,y,z)>(truncationValue*max1)){
               if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.))>distanceMax1){
                   distanceMax1 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.));
                   distanceMax2D1 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
                   distanceMaxZ1 = sqrt(pow(double(z-autocorrelation1.depth()/2),2.));
               }
           }
           if(autocorrelation2(x,y,z)>(truncationValue*max2)){
               if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.))>distanceMax2){
                   distanceMax2 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.));
                   distanceMax2D2 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
                   distanceMaxZ2 = sqrt(pow(double(z-autocorrelation1.depth()/2),2.));
               }
           }
       }
       cout << "distanceMax1: " << distanceMax1 << "   distanceMax2: " << distanceMax2 << endl;
       int maxRadius=0,maxDepth=0;
       if(w>h){maxRadius=h;}
       else{maxRadius=w;}
       if(distanceMax1>distanceMax2){radiusForTruncation = cimg::min(ceil(distanceMax1),maxRadius);}
       else{radiusForTruncation = cimg::min(ceil(distanceMax2),maxRadius);}
       if(distanceMax1>distanceMax2){depthForTruncation = cimg::min(ceil(distanceMax1),d);}
       else{depthForTruncation = cimg::min(ceil(distanceMax2),d);}
       cout << "radius for truncation: " << radiusForTruncation << endl;
   }
}

void computeDenominator(CImg<double> & denominator,CImg<> & autocorrelation1,CImg<> & autocorrelation2,int w,int h,int d,int currentFrame,int currentTemporalShift){
   double S1=0.,S2=0.,S3_1=0.,S3_2=0.;
   if(d==1){
       for(int i2=0;i2<h;i2++){
           for(int i1=0;i1<w;i1++){
               double S2_1=0.,S2_2=0.;
               for(int y=i2+1;y<(i2+h+1);y++){
                   for(int x=i1+1;x<(i1+w+1);x++){
                       S2_1 += autocorrelation1(x,y);
                       S2_2 += autocorrelation2(x,y);
                   }
               }
               S2 += S2_1*S2_2;
           }
       }
   }
   else{
       for(int i3=0;i3<d;i3++){
           for(int i2=0;i2<h;i2++){
               for(int i1=0;i1<w;i1++){
                   double S2_1=0.,S2_2=0.;
                   for(int z=i3+1;z<(i3+d+1);z++){
                       for(int y=i2+1;y<(i2+h+1);y++){
                           for(int x=i1+1;x<(i1+w+1);x++){
                               S2_1 += autocorrelation1(x,y,z);
                               S2_2 += autocorrelation2(x,y,z);
                           }
                       }
                   }
                   S2 += S2_1*S2_2;
               }
           }
       }
   }
   
   CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
   if(d==1){
       for(int y=-h;y<=h;y++){
           for(int x=-w;x<=w;x++){
               weightedMatrix(x+w,y+h) = (w-abs(x))*(h-abs(y));
           }
       }
       cimg_forXY(autocorrelation1,x,y){
           S1 += autocorrelation1(x,y)*autocorrelation2(x,y)*weightedMatrix(x,y);
           S3_1 += autocorrelation1(x,y)*weightedMatrix(x,y);
           S3_2 += autocorrelation2(x,y)*weightedMatrix(x,y);
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h),2.);
       denominator(currentFrame,1,currentTemporalShift) = 2.*S2/pow(double(w*h),3.);
       denominator(currentFrame,2,currentTemporalShift) = S3_1*S3_2/pow(double(w*h),4.);
   }
   else{
       for(int z=-d;z<=d;z++){
           for(int y=-h;y<=h;y++){
               for(int x=-w;x<=w;x++){
                   weightedMatrix(x+w,y+h,z+d) = (w-abs(x))*(h-abs(y))*(d-abs(z));
               }
           }
       }
       cimg_forXYZ(autocorrelation1,x,y,z){
           S1 += autocorrelation1(x,y,z)*autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
           S3_1 += autocorrelation1(x,y,z)*weightedMatrix(x,y,z);
           S3_2 += autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h*d),2.);
       denominator(currentFrame,1,currentTemporalShift) = 2.*S2/pow(double(w*h*d),3.);
       denominator(currentFrame,2,currentTemporalShift) = S3_1*S3_2/pow(double(w*h*d),4.);
   }
}

void computeTruncatedDenominator(CImg<double> & denominator,CImg<> & autocorrelation1,CImg<> & autocorrelation2,int w,int h,int d,int currentFrame,int radiusForTruncation,int depthForTruncation,int currentTemporalShift){
   double S1=0.,S2=0.,S3_1=0.,S3_2=0.;
   if(d==1){
       for(int i2=0;i2<h;i2++){
           for(int i1=0;i1<w;i1++){
               double S2_1=0.,S2_2=0.;
               for(int y=max(h-radiusForTruncation,i2+1);y<=min(h+radiusForTruncation,i2+h);y++){
                   for(int x=max(w-radiusForTruncation,i1+1);x<=min(w+radiusForTruncation,i1+w);x++){
                       S2_1 += autocorrelation1(x,y);
                       S2_2 += autocorrelation2(x,y);
                   }
               }
               S2 += S2_1*S2_2;
           }
       }
   }
   else{
       for(int i3=0;i3<d;i3++){
           for(int i2=0;i2<h;i2++){
               for(int i1=0;i1<w;i1++){
                   double S2_1=0.,S2_2=0.;
                   for(int z=max(d-radiusForTruncation,i3+1);z<=min(d+radiusForTruncation,i3+d);z++){
                       for(int y=max(h-radiusForTruncation,i2+1);y<=min(h+radiusForTruncation,i2+h);y++){
                           for(int x=max(w-radiusForTruncation,i1+1);x<=min(w+radiusForTruncation,i1+w);x++){
                               S2_1 += autocorrelation1(x,y,z);
                               S2_2 += autocorrelation2(x,y,z);
                           }
                       }
                   }
                   S2 += S2_1*S2_2;
               }
           }
       }
   }
   CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
   if(d==1){
       for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
           for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
               weightedMatrix(x+w,y+h) = (w-abs(x))*(h-abs(y));
           }
       }
       for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
           for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
               S1 += autocorrelation1(x,y)*autocorrelation2(x,y)*weightedMatrix(x,y);
               S3_1 += autocorrelation1(x,y)*weightedMatrix(x,y);
               S3_2 += autocorrelation2(x,y)*weightedMatrix(x,y);
           }
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h),2.);
       denominator(currentFrame,1,currentTemporalShift) = 2.*S2/pow(double(w*h),3.);
       denominator(currentFrame,2,currentTemporalShift) = S3_1*S3_2/pow(double(w*h),4.);
   }
   else{
       for(int z=-depthForTruncation;z<=depthForTruncation;z++){
           for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
               for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                   weightedMatrix(x+w,y+h,z+d) = (w-abs(x))*(h-abs(y))*(d-abs(z));
               }
           }
       }
       for(int z=-depthForTruncation;z<=depthForTruncation;z++){
           for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
               for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                   S1 += autocorrelation1(x+w,y+h,z+d)*autocorrelation2(x+w,y+h,z+d)*weightedMatrix(x+w,y+h,z+d);
                   S3_1 += autocorrelation1(x+w,y+h,z+d)*weightedMatrix(x+w,y+h,z+d);
                   S3_2 += autocorrelation2(x+w,y+h,z+d)*weightedMatrix(x+w,y+h,z+d);
               }
           }
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h*d),2.);
       denominator(currentFrame,1,currentTemporalShift) = 2.*S2/pow(double(w*h*d),3.);
       denominator(currentFrame,2,currentTemporalShift) = S3_1*S3_2/pow(double(w*h*d),4.);
   }
}

void computeS1(CImg<double> & denominator,CImg<> & autocorrelation1,CImg<> & autocorrelation2,int w,int h,int d,int currentFrame,int radiusForTruncation,int depthForTruncation,int currentTemporalShift){
   double S1=0.;
   
   CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
   if(d==1){
       for(int y=-h;y<=h;y++){
           for(int x=-w;x<=w;x++){
               weightedMatrix(x+w,y+h) = (w-abs(x))*(h-abs(y));
           }
       }
       cimg_forXY(autocorrelation1,x,y){
           S1 += autocorrelation1(x,y)*autocorrelation2(x,y)*weightedMatrix(x,y);
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h),2.);
   }
   else{
       for(int z=-d;z<=d;z++){
           for(int y=-h;y<=h;y++){
               for(int x=-w;x<=w;x++){
                   weightedMatrix(x+w,y+h,z+d) = (w-abs(x))*(h-abs(y))*(d-abs(z));
               }
           }
       }
       cimg_forXYZ(autocorrelation1,x,y,z){
           S1 += autocorrelation1(x,y,z)*autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h*d),2.);
   }
}

void computeTruncatedS1(CImg<double> & denominator,CImg<> & autocorrelation1,CImg<> & autocorrelation2,int w,int h,int d,int currentFrame,int radiusForTruncation,int depthForTruncation,int currentTemporalShift){
   double S1=0.;
   CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
   if(d==1){
       for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
           for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
               weightedMatrix(x+w,y+h) = (w-abs(x))*(h-abs(y));
           }
       }
       for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
           for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
               S1 += autocorrelation1(x+w,y+h)*autocorrelation2(x+w,y+h)*weightedMatrix(x+w,y+h);
           }
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h),2.);
   }
   else{
       for(int z=-depthForTruncation;z<=depthForTruncation;z++){
           for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
               for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                   weightedMatrix(x+w,y+h,z+d) = (w-abs(x))*(h-abs(y))*(d-abs(z));
               }
           }
       }
       
       for(int z=-depthForTruncation;z<=depthForTruncation;z++){
           for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
               for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                   S1 += autocorrelation1(x+w,y+h,z+d)*autocorrelation2(x+w,y+h,z+d)*weightedMatrix(x+w,y+h,z+d);
               }
           }
       }
       denominator(currentFrame,0,currentTemporalShift) = S1/pow(double(w*h*d),2.);
   }
}

void postTreatment(int nbFrames,vector<double> & pval,vector<int> & df,CImg<double> & pvalues){
   int nbRejects=0;
   for(unsigned int t=0;t<nbFrames;t++){
       pval.push_back(pvalues(t,3));
       if(pvalues(t,3)<0.05){nbRejects++;}
   }
   sort(pval.begin(),pval.end());
   
   df.push_back(0);
   for(unsigned int u=1;u<pval.size();u++){
       int nbPvals = df[df.size()-1];
       if(df[df.size()-1]<pval.size()){
           for(unsigned int v=df[df.size()-1];v<pval.size();v++){
               if(pval[v]<(double(u)/double(pval.size()))){
                   nbPvals++;
               }
           }
           df.push_back(nbPvals);
       }
       else{
           df.push_back(pval.size());
       }
   }
}


         main          //////////////////////////////////////////////////////////////////



int main(int argc,char** argv){
   
   TIFFSetWarningHandler(NULL);
   
   const char * inputImage1  = cimg_option("-i1",(char*)NULL,"input image 1");
   const char * inputImage2  = cimg_option("-i2",(char*)NULL,"input image 2 (if input image 1 is not a 2 channel image)");
   const char * inputCellMask  = cimg_option("-cm",(char*)NULL,"cell mask");
   const char * outputFilename  = cimg_option("-o",(char*)NULL,"output file name");
   const bool wholeDenominator = cimg_option("-d",false,"true: computes S1,S2,S3; false: only computes S1");
   const bool saveImagesAsExamples = cimg_option("-s",false,"true: save images for figures; false: process colocalization test");
   const int nbThreads = cimg_option("-nbt",4,"number of cores");
   const int maxImages = cimg_option("-mi",0,"maximum number of images in the sequence");
   const int nbPlanes = cimg_option("-nbp",1,"number of planes for each image");
   const int planeDivider = cimg_option("-pd",-1,"planes to consider if not all");
   //const int dim = cimg_option("-d",128,"image dimension");
   //const bool truncatedCovariance = cimg_option("-t",true,"remove autocovariance coefficients corresponding to noise");
   const float truncationValue = cimg_option("-tv",0.1,"truncation value to remove autocovariance coefficients corresponding to noise");
   const bool textfile = cimg_option("-tf",false,"input files as text files");
   const bool denseEstimation = cimg_option("-de",false,"correlation and pvalue estimation at each point in the image");
   const int maxShift = cimg_option("-ms",0,"temporal shift between sequences");
   const int channel1 = cimg_option("-c1",-1,"first channel to consider");
   const int channel2 = cimg_option("-c2",-1,"first channel to consider");
   const int tfWidth = cimg_option("-tfw",250,"width for images from textfile");
   const int tfHeight = cimg_option("-tfh",250,"height for images from textfile");
   const int tfDepth= cimg_option("-tfd",20,"depth for images from textfile");
   const int tfFirstFrame= cimg_option("-tff",1,"first image index from textfile");
   const int tfLastFrame= cimg_option("-tfl",100,"last image index from textfile");
   const bool randomizationStudy = cimg_option("-r",false,"randomization study");
   const int randomizationWindowWidth = cimg_option("-rww",50,"randomization window width");
   const int randomizationWindowHeight = cimg_option("-rwh",50,"randomization window height");
   const int randomizationWindowDepth = cimg_option("-rwd",5,"randomization window depth");
   const int windowWidthStep = cimg_option("-wws",50,"randomization window width");
   const int windowHeightStep = cimg_option("-whs",50,"randomization window height");
   const int windowDepthStep = cimg_option("-wds",5,"randomization window depth");
   
   cimg::tic();
   
   if(inputImage1!=NULL){
       
       // inputs
       CImgList<> input1,input2;
       CImg<> mask;

       // load sequences
       loadSequences(inputImage1,inputImage2,textfile,tfWidth,tfHeight,tfDepth,tfFirstFrame,tfLastFrame,nbPlanes,planeDivider,maxImages,channel1,channel2,input1,input2);
       
       if((input1.size()>0)&&(input2.size()>0)){

       	if(inputCellMask!=NULL){
       		mask.assign(inputCellMask);
       	}
       	else{
       		mask.assign(input1[0].width(),input1[0].height(),input1[0].depth(),1,1);
       	}

           // temporal analysis
           if(maxShift>0){
               omp_set_num_threads(nbThreads);

               // output
               CImg<double> statisticalTest(input1.size(),2,maxShift+1,1,-1.),correlation(input1.size(),maxShift+1,1,1,-1.),pvalues(input1.size(),maxShift+1,1,1,-1.),denominator(input1.size(),3,maxShift+1,1,-1.);
               string sortedPvaluesFilename,pValueDistributionFunctionFilename,correlationFilename;//denominatorFile,statFile,pvalFile;
               if(outputFilename!=NULL){
                   sortedPvaluesFilename = outputFilename;
                   sortedPvaluesFilename += "_sortedPvalues.txt";
                   pValueDistributionFunctionFilename = outputFilename;
                   pValueDistributionFunctionFilename += "_pValueDistributionFunction.txt";
                   correlationFilename = outputFilename;
                   correlationFilename += "_correlation.txt";
               }
               else{
                   sortedPvaluesFilename = "sortedPvalues.txt";
                   pValueDistributionFunctionFilename = "pValueDistributionFunction.txt";
                   correlationFilename += "correlation.txt";
               }
               ofstream sortedPvalues(sortedPvaluesFilename.c_str()),pValueDistributionFunction(pValueDistributionFunctionFilename.c_str()),correlationOutput(correlationFilename.c_str());
               string statTestFilename;
               if(outputFilename!=NULL){
               	statTestFilename = outputFilename;
               	statTestFilename += "_statTest.txt";
               }
               else{
               	statTestFilename = "statTest.txt";
               }
               ofstream statTest(statTestFilename.c_str());
               
               // new dimensions for FFT
               int powerX=2,powerY=2,powerZ=2;
               while(2*input1[0].width()>pow(2,double(powerX))){powerX++;}
               while(2*input1[0].height()>pow(2,double(powerY))){powerY++;}
               if(input1[0].depth()>1){
                   while(2*input1[0].depth()>pow(2,double(powerZ))){powerZ++;}
               }
               else{
                   powerZ=0;
               }
               
               CImgList<> autocorrelationOverTime1,autocorrelationOverTime2;
               
               for(int t=0;t<input1.size();t++){
                   CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0);
                   cimg_forXYZ(input1[t],x,y,z){
                       input1Tmp(x,y,z) = input1[t](x,y,z);
                       input2Tmp(x,y,z) = input2[t](x,y,z);
                   }
                   
                   // resizing
                   CImg<> currentAutocorrelation1=autocorrelation(input1Tmp,input1[t].width(),input1[t].height(),input1[t].depth()),
                   currentAutocorrelation2=autocorrelation(input2Tmp,input1[t].width(),input1[t].height(),input1[t].depth());
                   autocorrelationOverTime1.push_back(currentAutocorrelation1);
                   autocorrelationOverTime2.push_back(currentAutocorrelation2);
               }
               
               for(int shift=0;shift<=maxShift;shift++){
                   // statistical test estimation over time
#pragma omp parallel for
               	for(int t=0;t<(input1.size()-maxShift);t++){//for(int t=0;t<(input1.size()-shift);t++){//
                       CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0);
                       cimg_forXYZ(input1[t],x,y,z){
                           input1Tmp(x,y,z) = input1[t](x,y,z);
                           input2Tmp(x,y,z) = input2[t+shift](x,y,z);
                       }
                       CImg<> autocorrelation1=autocorrelationOverTime1[t],autocorrelation2=autocorrelationOverTime2[t+shift];
                       
                       // probas
                       CImg<> intersection(input1[t].width(),input1[t].height(),input1[t].depth(),1,0);
                       cimg_forXYZ(intersection,x,y,z){
                           intersection(x,y,z) = input1Tmp(x,y,z)*input2Tmp(x,y,z);
                       }
                       double p1=input1Tmp.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth()),p2=input2Tmp.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth()),p12=intersection.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth());
                       
                       // centering covariances
                       cimg_forXYZ(autocorrelation1,x,y,z){
                           autocorrelation1(x,y,z) -= pow(p1,2.);
                           autocorrelation2(x,y,z) -= pow(p2,2.);
                       }
                       
                       // remove coefficients corresponding to noise
                       int radiusForTruncation=0,depthForTruncation=0;
                       if(truncationValue>0.){
                           radiusAndDepthTruncations(autocorrelation1,autocorrelation2,truncationValue,input1[t].width(),input1[t].height(),input1[t].depth(),radiusForTruncation,depthForTruncation);
                       }
                       /*if(input1[t].depth()==1){
                        double max1=autocorrelation1.max(),max2=autocorrelation2.max(),distanceMax1=0.,distanceMax2=0.;
                        cimg_forXY(autocorrelation1,x,y){
                        if(autocorrelation1(x,y)>(truncationValue*max1)){
                        if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.))>distanceMax1){
                        distanceMax1 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
                        }
                        }
                        if(autocorrelation2(x,y)>(truncationValue*max2)){
                        if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.))>distanceMax2){
                        distanceMax2 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.));
                        }
                        }
                        }
                        int maxRadius=0;
                        if(input1[t].width()>input1[t].height()){maxRadius=input1[t].height();}
                        else{maxRadius=input1[t].width();}
                        if(distanceMax1>distanceMax2){radiusForTruncation = cimg::min(ceil(distanceMax1),maxRadius);}
                        else{radiusForTruncation = cimg::min(ceil(distanceMax2),maxRadius);}
                        }
                        else{
                        double max1=autocorrelation1.max(),max2=autocorrelation2.max(),distanceMax1=0.,distanceMax2=0.;//,distanceMax2D1=0.,distanceMax2D2=0.,distanceMaxZ1=0.,distanceMaxZ2=0.;
                        cimg_forXYZ(autocorrelation1,x,y,z){
                        if(autocorrelation1(x,y,z)>(truncationValue*max1)){
                        if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.))>distanceMax1){
                        distanceMax1 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.));
                        }
                        }
                        if(autocorrelation2(x,y,z)>(truncationValue*max2)){
                        if(sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.))>distanceMax2){
                        distanceMax2 = sqrt(pow(double(x-autocorrelation1.width()/2),2.)+pow(double(y-autocorrelation1.height()/2),2.)+pow(double(z-autocorrelation1.depth()/2),2.));
                        }
                        }
                        }
                        int maxRadius=0,maxDepth=0;
                        if(input1[t].width()>input1[t].height()){maxRadius=input1[t].height();}
                        else{maxRadius=input1[t].width();}
                        if(distanceMax1>distanceMax2){radiusForTruncation = cimg::min(ceil(distanceMax1),maxRadius);}
                        else{radiusForTruncation = cimg::min(ceil(distanceMax2),maxRadius);}
                        if(distanceMax1>distanceMax2){depthForTruncation = cimg::min(ceil(distanceMax1),input1[t].depth());}
                        else{depthForTruncation = cimg::min(ceil(distanceMax2),input1[t].depth());}
                        cout << maxRadius << endl;
                        }
                        }*/
                       
                       if(radiusForTruncation==0){
                           if(wholeDenominator){
                               computeDenominator(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,shift);
                               // outputs
                               // stat test whole denominator
                               statisticalTest(t,0,shift) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                               // stat test S1
                               statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0));
                               // two-sided p-values whole denominator
                               pvalues(t,0,shift) = (1.-CND(-statisticalTest(t,0)));
                               // two-sided p-values S1
                               pvalues(t,1,shift) = (1.-CND(-statisticalTest(t,1)));
                               // one-sided p-values whole denominator
                               pvalues(t,2,shift) = (1.-CND(statisticalTest(t,0)));
                               // one-sided p-values S1
                               pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1)));
                               // normalized correlation
                               correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           }
                           // only computes S1
                           else{
                               // S1 computation
                               computeS1(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,shift);
                               // stat test S1
                               statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0));
                               // two-sided p-values S1
                               pvalues(t,1,shift) = (1.-CND(-statisticalTest(t,1)));
                               // one-sided p-values S1
                               pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1)));
                               // normalized correlation
                               correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           }
                           
                           /*// computes S1, S2 and S3
                            if(wholeDenominator){
                            // denominator computation
                            double S1=0.,S2=0.,S3_1=0.,S3_2=0.;
                            if(input1[t].depth()==1){
                            cimg_forXY(input1[t],i1,i2){
                            double S2_1=0.,S2_2=0.;
                            for(int y=i2+1;y<(i2+input1[t].height()+1);y++){
                            for(int x=i1+1;x<(i1+input1[t].width()+1);x++){
                            S2_1 += autocorrelation1(x,y);
                            S2_2 += autocorrelation2(x,y);
                            }
                            }
                            S2 += S2_1*S2_2;
                            }
                            }
                            else{
                            cimg_forXYZ(input1[t],i1,i2,i3){
                            double S2_1=0.,S2_2=0.;
                            for(int z=i3+1;z<(i3+input1[t].depth()+1);z++){
                            for(int y=i2+1;y<(i2+input1[t].height()+1);y++){
                            for(int x=i1+1;x<(i1+input1[t].width()+1);x++){
                            S2_1 += autocorrelation1(x,y,z);
                            S2_2 += autocorrelation2(x,y,z);
                            }
                            }
                            }
                            S2 += S2_1*S2_2;
                            }
                            }
                            
                            CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                            if(input1[t].depth()==1){
                            for(int y=-input1[t].height();y<=input1[t].height();y++){
                            for(int x=-input1[t].width();x<=input1[t].width();x++){
                            weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                            }
                            }
                            cimg_forXY(autocorrelation1,x,y){
                            S1 += autocorrelation1(x,y)*autocorrelation2(x,y)*weightedMatrix(x,y);
                            S3_1 += autocorrelation1(x,y)*weightedMatrix(x,y);
                            S3_2 += autocorrelation2(x,y)*weightedMatrix(x,y);
                            }
                            denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                            denominator(t,1,shift) = 2.*S2/pow(double(input1[t].width()*input1[t].height()),3.);
                            denominator(t,2,shift) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()),4.);
                            }
                            else{
                            for(int z=-input1[t].depth();z<=input1[t].depth();z++){
                            for(int y=-input1[t].height();y<=input1[t].height();y++){
                            for(int x=-input1[t].width();x<=input1[t].width();x++){
                            weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                            }
                            }
                            }
                            cimg_forXYZ(autocorrelation1,x,y,z){
                            S1 += autocorrelation1(x,y,z)*autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
                            S3_1 += autocorrelation1(x,y,z)*weightedMatrix(x,y,z);
                            S3_2 += autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
                            }
                            denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                            denominator(t,1,shift) = 2.*S2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),3.);
                            denominator(t,2,shift) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),4.);
                            }
                            
                            statisticalTest(t,0,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift)-denominator(t,1,shift)+denominator(t,2,shift));
                            statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                            
                            pvalues(t,0,shift) = 2.*(1.-CND(abs(statisticalTest(t,0,shift))));
                            pvalues(t,1,shift) = 2.*(1.-CND(abs(statisticalTest(t,1,shift))));
                            
                            pvalues(t,2,shift) = (1.-CND(statisticalTest(t,0,shift)));
                            pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1,shift)));
                            
                            correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                            }
                            // only computes S1
                            else{
                            // denominator computation
                            double S1=0.;
                            CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                            if(input1[t].depth()==1){
                            for(int y=-input1[t].height();y<=input1[t].height();y++){
                            for(int x=-input1[t].width();x<=input1[t].width();x++){
                            weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                            }
                            }
                            cimg_forXY(autocorrelation1,x,y){
                            S1 += autocorrelation1(x,y)*autocorrelation2(x,y)*weightedMatrix(x,y);
                            }
                            denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                            }
                            else{
                            for(int z=-input1[t].depth();z<=input1[t].depth();z++){
                            for(int y=-input1[t].height();y<=input1[t].height();y++){
                            for(int x=-input1[t].width();x<=input1[t].width();x++){
                            weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                            }
                            }
                            }
                            cimg_forXYZ(autocorrelation1,x,y,z){
                            S1 += autocorrelation1(x,y,z)*autocorrelation2(x,y,z)*weightedMatrix(x,y,z);
                            }
                            denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                            }
                            
                            statisticalTest(t,0,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                            statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                            
                            pvalues(t,0,shift) = 2.*(1.-CND(abs(statisticalTest(t,0,shift))));
                            pvalues(t,1,shift) = 2.*(1.-CND(abs(statisticalTest(t,1,shift))));
                            
                            pvalues(t,2,shift) = (1.-CND(statisticalTest(t,0,shift)));
                            pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1,shift)));
                            
                            correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                            }*/
                       }
                       else{
                           // computes S1, S2 and S3
                           if(wholeDenominator){
                               // denominator computation
                               computeTruncatedDenominator(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,shift);
                               /*double S1=0.,S2=0.,S3_1=0.,S3_2=0.;
                                if(input1[t].depth()==1){
                                for(int i2=0;i2<input1[t].height();i2++){
                                for(int i1=0;i1<input1[t].width();i1++){
                                double S2_1=0.,S2_2=0.;
                                for(int y=max(input1[t].height()-radiusForTruncation,i2+1);y<=min(input1[t].height()+radiusForTruncation,i2+input1[t].height());y++){
                                for(int x=max(input1[t].width()-radiusForTruncation,i1+1);x<=min(input1[t].width()+radiusForTruncation,i1+input1[t].width());x++){
                                S2_1 += autocorrelation1(x,y);
                                S2_2 += autocorrelation2(x,y);
                                }
                                }
                                S2 += S2_1*S2_2;
                                }
                                }
                                }
                                else{
                                for(int i3=0;i3<input1[t].depth();i3++){
                                for(int i2=0;i2<input1[t].height();i2++){
                                for(int i1=0;i1<input1[t].width();i1++){
                                double S2_1=0.,S2_2=0.;
                                for(int z=max(input1[t].depth()-depthForTruncation,i3+1);z<=min(input1[t].depth()+depthForTruncation,i3+input1[t].depth());z++){
                                for(int y=max(input1[t].height()-radiusForTruncation,i2+1);y<=min(input1[t].height()+radiusForTruncation,i2+input1[t].height());y++){
                                for(int x=max(input1[t].width()-radiusForTruncation,i1+1);x<=min(input1[t].width()+radiusForTruncation,i1+input1[t].width());x++){
                                S2_1 += autocorrelation1(x,y,z);
                                S2_2 += autocorrelation2(x,y,z);
                                }
                                }
                                }
                                S2 += S2_1*S2_2;
                                }
                                }
                                }
                                }
                                CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                                if(input1[t].depth()==1){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                                }
                                }
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                S3_1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                S3_2 += autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                }
                                }
                                denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                                denominator(t,1,shift) = 2.*S2/pow(double(input1[t].width()*input1[t].height()),3.);
                                denominator(t,2,shift) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()),4.);
                                }
                                else{
                                for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                                }
                                }
                                }
                                for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                S3_1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                S3_2 += autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                }
                                }
                                }
                                denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                                denominator(t,1,shift) = 2.*S2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),3.);
                                denominator(t,2,shift) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),4.);
                                }*/
                               
                               statisticalTest(t,0,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift)-denominator(t,1,shift)+denominator(t,2,shift));
                               statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                               
                               pvalues(t,0,shift) = 2.*(1.-CND(abs(statisticalTest(t,0,shift))));
                               pvalues(t,1,shift) = 2.*(1.-CND(abs(statisticalTest(t,1,shift))));
                               
                               pvalues(t,2,shift) = (1.-CND(statisticalTest(t,0,shift)));
                               pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1,shift)));
                               
                               correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           }
                           else{
                               // denominator computation
                               computeTruncatedS1(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,shift);
                               /*double S1=0.;
                                CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                                if(input1[t].depth()==1){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                                }
                                }
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                }
                                }
                                denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                                }
                                else{
                                for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                                }
                                }
                                }
                                for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                }
                                }
                                }
                                denominator(t,0,shift) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                                }*/
                               
                               statisticalTest(t,0,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                               statisticalTest(t,1,shift) = (p12-p1*p2)/sqrt(denominator(t,0,shift));
                               
                               pvalues(t,0,shift) = 2.*(1.-CND(abs(statisticalTest(t,0,shift))));
                               pvalues(t,1,shift) = 2.*(1.-CND(abs(statisticalTest(t,1,shift))));
                               
                               pvalues(t,2,shift) = (1.-CND(statisticalTest(t,0,shift)));
                               pvalues(t,3,shift) = (1.-CND(statisticalTest(t,1,shift)));
                               
                               correlation(t,shift) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           }
                       }
                   }
                   vector<double> pval;
                   vector<int> distributionFunction;
                   int nbRejects=0;
                   for(unsigned int t=0;t<input1.size();t++){
                       pval.push_back(pvalues(t,2,shift));
                       if(pvalues(t,2,shift)<0.05){nbRejects++;}
                   }
                   sort(pval.begin(),pval.end());
                   
                   distributionFunction.push_back(0);
                   for(unsigned int u=1;u<(pval.size()-maxShift);u++){
                       int nbPvals = distributionFunction[distributionFunction.size()-1];
                       if(distributionFunction[distributionFunction.size()-1]<pval.size()){
                           for(unsigned int v=distributionFunction[distributionFunction.size()-1];v<pval.size();v++){
                               if(pval[v]<(double(u)/double(pval.size()))){
                                   nbPvals++;
                               }
                           }
                           distributionFunction.push_back(nbPvals);
                       }
                       else{
                           distributionFunction.push_back(pval.size());
                       }
                   }
                   
                   cout << "Shift: " << shift << "   p-value: min: " << pval[0] << "   lower quartile: " << pval[pval.size()/4] << "   median: " << pval[pval.size()/2] << "   upper quartile: " << pval[3*pval.size()/4] << "   max: " << pval[pval.size()-1] << endl;
                   
                   sortedPvalues << pval[0];
                   pValueDistributionFunction << distributionFunction[0];
                   correlationOutput << correlation(0,shift);
                   for(unsigned int i=1;i<(pval.size()-maxShift);i++){
                       sortedPvalues << '\t' << pval[i];
                       pValueDistributionFunction << '\t' << distributionFunction[i];
                       correlationOutput << '\t' << correlation(i,shift);
                   }
                   sortedPvalues << endl;
                   pValueDistributionFunction << endl;
                   correlationOutput << endl;
               }
               cimg_forXZ(statisticalTest,t,shift){
               	statTest << t << '\t' << shift << '\t' << statisticalTest(t,1,shift) << endl;
               }
           }
           
           // image per image test
           else{
               if(randomizationStudy){
                   omp_set_num_threads(nbThreads);
                   int nbRandomizations = mask.sum()/500;
                   // output files
                   string sortedPvaluesFilename,pValueDistributionFunctionFilename,correlationFilename,colocalizationDisplayHitsFilename,colocalizationDisplayOverlayFilename;
                   if(outputFilename!=NULL){
                       sortedPvaluesFilename = outputFilename;
                       sortedPvaluesFilename += "_sortedPvalues.txt";
                       pValueDistributionFunctionFilename = outputFilename;
                       pValueDistributionFunctionFilename += "_pValueDistributionFunction.txt";
                       correlationFilename = outputFilename;
                       correlationFilename += "_correlation.txt";
                       colocalizationDisplayHitsFilename = outputFilename;
                       colocalizationDisplayHitsFilename += "_displayHits.tif";
                       colocalizationDisplayOverlayFilename = outputFilename;
                       colocalizationDisplayOverlayFilename += "_displayOverlay.tif";
                   }
                   else{
                       sortedPvaluesFilename = "sortedPvalues.txt";
                       pValueDistributionFunctionFilename = "pValueDistributionFunction.txt";
                       correlationFilename = "correlation.txt";
                       colocalizationDisplayHitsFilename = "displayHits.tif";
                       colocalizationDisplayOverlayFilename = "displayOverlay.tif";
                   }
                   ofstream sortedPvalues(sortedPvaluesFilename.c_str()),pValueDistributionFunction(pValueDistributionFunctionFilename.c_str()),correlationOutput(correlationFilename.c_str());//outputDenominator(denominatorFile.c_str()),outputStat(statFile.c_str()),outputPval(pvalFile.c_str());
                   
                   // output
                   CImg<double> statisticalTest(input1.size(),2,nbRandomizations,1,0.),correlation(input1.size(),nbRandomizations,1,1,0.),pvalues(input1.size(),4,nbRandomizations,1,0.),denominator(input1.size(),3,nbRandomizations,1,0.);
                   // new dimensions for FFT
                   int powerX=2,powerY=2,powerZ=2;
                   while(2*randomizationWindowWidth>pow(2,double(powerX))){powerX++;}
                   while(2*randomizationWindowHeight>pow(2,double(powerY))){powerY++;}
                   if(randomizationWindowDepth>1){
                       while(2*randomizationWindowDepth>pow(2,double(powerZ))){powerZ++;}
                   }
                   else{
                       powerZ=0;
                   }
                   
                   // statistical test estimation over time
                   CImgList<> colocalizationDisplayHits,colocalizationDisplayOverlay;
                   const unsigned char blue[] = { 0,0,1 },white[] = { 1,1,1 };
#pragma omp parallel for
                   for(int t=0;t<input1.size();t++){
                       CImg<> currentDisplay1(input1[t].width(),input1[t].height(),input1[t].depth(),3,0),currentDisplay2(input1[t].width(),input1[t].height(),input1[t].depth(),3,0);
                       cimg_forXYZ(currentDisplay1,x,y,z){
                           currentDisplay1(x,y,z,0) = input1[t](x,y,z);
                           currentDisplay1(x,y,z,1) = input2[t](x,y,z);
                       }

                       for(int r=0;r<nbRandomizations;r++){
                           CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0);
                           double rand1,rand2,rand3;
                           int rand1Shift,rand2Shift,rand3Shift;
                           bool out=false;
                           while(!out){
                           	rand1=cimg::rand();rand2=cimg::rand();rand3=cimg::rand();
                           	rand1Shift=int(rand1*(double(input1[t].width())-double(randomizationWindowWidth)));
                           	rand2Shift=int(rand2*(double(input1[t].height())-double(randomizationWindowHeight)));
                           	if(randomizationWindowDepth>1){
                           		rand3Shift=int(rand3*(double(input1[t].depth())-double(randomizationWindowDepth)));
                           		if(mask(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2,rand3Shift+randomizationWindowDepth/2)>0){out=true;}
                           	}
                           	else{
                           		if(mask(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2)>0){out=true;}
                           	}
                           }

                           if(randomizationWindowDepth>1){
                               for(int z=rand3Shift;z<(randomizationWindowDepth+rand3Shift);z++){
                                   for(int y=rand2Shift;y<(randomizationWindowHeight+rand2Shift);y++){
                                       for(int x=rand1Shift;x<(randomizationWindowWidth+rand1Shift);x++){
                                       	if(mask(x,y,z)>0){
                                       		input1Tmp(x-rand1Shift,y-rand2Shift,z-rand3Shift) = input1[t](x,y,z);
                                       		input2Tmp(x-rand1Shift,y-rand2Shift,z-rand3Shift) = input2[t](x,y,z);
                                       	}
                                       }
                                   }
                               }
                           }
                           else{
                               for(int y=rand2Shift;y<(randomizationWindowHeight+rand2Shift);y++){
                                   for(int x=rand1Shift;x<(randomizationWindowWidth+rand1Shift);x++){
                                   	if(mask(x,y)>0){
                                   		input1Tmp(x-rand1Shift,y-rand2Shift) = input1[t](x,y);
                                   		input2Tmp(x-rand1Shift,y-rand2Shift) = input2[t](x,y);
                                   	}
                                   }
                               }
                           }
                           if((input1Tmp.sum()>0)&&(input2Tmp.sum()>0)){
                           	// autocorrelation computation
                           	CImg<> autocorrelation1=autocorrelation(input1Tmp,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth),
                           			autocorrelation2=autocorrelation(input2Tmp,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth);
                           	// probas
                           	CImg<> intersection(randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,1,0);
                           	cimg_forXYZ(intersection,x,y,z){
                           		intersection(x,y,z) = input1Tmp(x,y,z)*input2Tmp(x,y,z);
                           	}
                           	double p1=input1Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth),p2=input2Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth),p12=intersection.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth);
                           	// centering covariances
                           	cimg_forXYZ(autocorrelation1,x,y,z){
                           		autocorrelation1(x,y,z) -= pow(p1,2.);
                           		autocorrelation2(x,y,z) -= pow(p2,2.);
                           	}
                           	// remove coefficients corresponding to noise
                           	int radiusForTruncation=0,depthForTruncation=0;
                           	if(truncationValue>0.){
                           		radiusAndDepthTruncations(autocorrelation1,autocorrelation2,truncationValue,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,radiusForTruncation,depthForTruncation);
                           	}

                           	// with truncation
                           	if(radiusForTruncation==0){
                           		// computes S1, S2 and S3
                           		if(wholeDenominator){
                           			computeDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,t,0);
                           			// outputs
                           			// stat test whole denominator
                           			statisticalTest(t,0,r) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                           			// stat test S1
                           			statisticalTest(t,1,r) = (p12-p1*p2)/sqrt(denominator(t,0));
                           			// two-sided p-values whole denominator
                           			pvalues(t,0,r) = 2.*(1.-CND(abs(statisticalTest(t,0,r))));
                           			// two-sided p-values S1
                           			pvalues(t,1,r) = 2.*(1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1 anti-coloc
                           			pvalues(t,2,r) = (1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1
                           			pvalues(t,3,r) = (1.-CND(statisticalTest(t,1,r)));
                           			// normalized correlation
                           			correlation(t,r) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           		}
                           		// only computes S1
                           		else{
                           			// S1 computation
                           			computeS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,t,radiusForTruncation,depthForTruncation,0);
                           			// stat test S1
                           			statisticalTest(t,1,r) = (p12-p1*p2)/sqrt(denominator(t,0));
                           			// two-sided p-values S1
                           			pvalues(t,1,r) = 2.*(1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1 anti-coloc
                           			pvalues(t,2,r) = (1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1
                           			pvalues(t,3,r) = (1.-CND(statisticalTest(t,1,r)));
                           			// normalized correlation
                           			correlation(t,r) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           		}
                           	}
                           	else{
                           		// computes S1, S2 and S3
                           		if(wholeDenominator){
                           			// denominator computation
                           			computeTruncatedDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,t,radiusForTruncation,depthForTruncation,0);
                           			// stat test whole denominator
                           			statisticalTest(t,0,r) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                           			// stat test S1
                           			statisticalTest(t,1,r) = (p12-p1*p2)/sqrt(denominator(t,0));
                           			// two-sided p-values whole denominator
                           			pvalues(t,0,r) = 2.*(1.-CND(abs(statisticalTest(t,0,r))));
                           			// two-sided p-values S1
                           			pvalues(t,1,r) = 2.*(1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1 anti-coloc
                           			pvalues(t,2,r) = (1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1
                           			pvalues(t,3,r) = (1.-CND(statisticalTest(t,1,r)));
                           			// normalized correlation
                           			correlation(t,r) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           		}
                           		// only compute S1
                           		else{
                           			// S1 computation
                           			computeTruncatedS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,t,radiusForTruncation,depthForTruncation,0);
                           			// stat test S1
                           			statisticalTest(t,1,r) = (p12-p1*p2)/sqrt(denominator(t,0));
                           			// two-sided p-values S1
                           			pvalues(t,1,r) = 2.*(1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1 anti-coloc
                           			pvalues(t,2,r) = (1.-CND(abs(statisticalTest(t,1,r))));
                           			// one-sided p-values S1
                           			pvalues(t,3,r) = (1.-CND(statisticalTest(t,1,r)));
                           			// normalized correlation
                           			correlation(t,r) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                           		}
                           	}
                           	if(correlation(t,r)>0){
                           		if(pvalues(t,3,r)<0.05){
                           			currentDisplay1.draw_circle(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2,4,white,1,1);
                           			currentDisplay2.draw_circle(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2,4,white,1,1);
                           		}
                           	}
                           	else{
                           		if(pvalues(t,2,r)<0.05){
                           			currentDisplay1.draw_circle(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2,4,blue,1,0);
                           			currentDisplay2.draw_circle(rand1Shift+randomizationWindowWidth/2,rand2Shift+randomizationWindowHeight/2,4,blue,1,0);
                           		}
                           	}
                           }
                           else{
                           	r--;
                           }
                       }
                       vector<double> pval;
                       vector<int> distributionFunction;
                       int nbRejects=0;
                       for(unsigned int r=0;r<nbRandomizations;r++){
                           pval.push_back(pvalues(t,3,r));
                           if(pvalues(t,3,r)<0.05){nbRejects++;}
                       }
                       sort(pval.begin(),pval.end());
                       
                       // indicator output
                       cout << "t: " << t << "    Rejections: " << double(nbRejects)/double(nbRandomizations) << "   p-value: min: " << pval[0] << "   lower quartile: " << pval[pval.size()/4] << "   median: " << pval[pval.size()/2] << "   upper quartile: " << pval[3*pval.size()/4] << "   max: " << pval[pval.size()-1] << endl;
                       
                       distributionFunction.push_back(0);
                       for(unsigned int u=1;u<(pval.size()-maxShift);u++){
                           int nbPvals = distributionFunction[distributionFunction.size()-1];
                           if(distributionFunction[distributionFunction.size()-1]<pval.size()){
                               for(unsigned int v=distributionFunction[distributionFunction.size()-1];v<pval.size();v++){
                                   if(pval[v]<(double(u)/double(pval.size()))){
                                       nbPvals++;
                                   }
                               }
                               distributionFunction.push_back(nbPvals);
                           }
                           else{
                               distributionFunction.push_back(pval.size());
                           }
                       }
                       
                       sortedPvalues << pval[0];
                       pValueDistributionFunction << distributionFunction[0];
                       correlationOutput << correlation(t,0);
                       for(unsigned int r=0;r<nbRandomizations;r++){
                           sortedPvalues << '\t' << pval[r];
                           pValueDistributionFunction << '\t' << distributionFunction[r];
                           correlationOutput << '\t' << correlation(t,r);
                       }
                       sortedPvalues << endl;
                       pValueDistributionFunction << endl;
                       correlationOutput << endl;
                       colocalizationDisplayOverlay.push_back(currentDisplay1);
                       colocalizationDisplayHits.push_back(currentDisplay2);
                   }
                   colocalizationDisplayOverlay.save(colocalizationDisplayOverlayFilename.c_str());
                   colocalizationDisplayHits.save(colocalizationDisplayHitsFilename.c_str());
               }
               else{
                   if(denseEstimation){
                       omp_set_num_threads(nbThreads);
                       // output files
                       string pvaluesFilename,correlationFilename,statisticalTestFilename;
                       if(outputFilename!=NULL){
                           pvaluesFilename = outputFilename;
                           pvaluesFilename += "_pvalues.tif";
                           correlationFilename = outputFilename;
                           correlationFilename += "_correlation.tif";
                           statisticalTestFilename = outputFilename;
                           statisticalTestFilename += "_denseField.tif";
                       }
                       else{
                           pvaluesFilename = "pvalues.tif";
                           correlationFilename = "correlation.tif";
                           statisticalTestFilename = "denseField.tif";
                       }
                       
                       // output
                       CImgList<double> statisticalTest,correlation,pvalues,statisticalTestDisplay;
                       // new dimensions for FFT
                       int powerX=2,powerY=2,powerZ=2;
                       while(2*randomizationWindowWidth>pow(2,double(powerX))){powerX++;}
                       while(2*randomizationWindowHeight>pow(2,double(powerY))){powerY++;}
                       if(randomizationWindowDepth>1){
                           while(2*randomizationWindowDepth>pow(2,double(powerZ))){powerZ++;}
                       }
                       else{
                           powerZ=0;
                       }
                       
                       // statistical test estimation over time
                       CImgList<> colocalizationDisplayHits,colocalizationDisplayOverlay;
                       const unsigned char blue[] = { 0,0,1 },white[] = { 1,1,1 };
                       for(int t=0;t<input1.size();t++){

                       	CImg<float> currentDisplay(input1[0].width(),input1[0].height(),input1[0].depth(),1,0);
                           CImg<double> currentStatisticalTest(input1[0].width(),input1[0].height(),input1[0].depth(),2,0.),currentCorrelation(input1[0].width(),input1[0].height(),input1[0].depth(),1,0.),currentPvalues(input1[0].width(),input1[0].height(),input1[0].depth(),1,0.);
						CImg<double> currentStatisticalTest2(input1[0].width(),input1[0].height(),1,1,0.),currentCorrelation2(input1[0].width(),input1[0].height(),1,1,0.),currentPvalues2(input1[0].width(),input1[0].height(),1,1,0.);
                           if(randomizationWindowDepth>1){
                           	for(int z=9;z<10;z++){
                               //for(int z=0;z<(input1[0].depth()-randomizationWindowDepth);z+=windowDepthStep){
#pragma omp parallel for
                                   for(int y=0;y<(input1[0].height()-randomizationWindowHeight);y+=windowHeightStep){
                                       for(int x=0;x<(input1[0].width()-randomizationWindowWidth);x+=windowWidthStep){
                                           CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0);
                                           for(int k=0;k<randomizationWindowDepth;k++){
                                               for(int j=0;j<randomizationWindowHeight;j++){
                                                   for(int i=0;i<randomizationWindowWidth;i++){
                                                   	if((mask(x+i,y+j,z+k)>0)&&(sqrt(pow(i-randomizationWindowWidth/2,2.)+pow(j-randomizationWindowHeight/2,2.))<(randomizationWindowWidth/2+.01))){
                                                   		input1Tmp(i,j,k) = input1[t](x+i,y+j,z+k);
                                                   		input2Tmp(i,j,k) = input2[t](x+i,y+j,z+k);
                                                   	}
                                                   }
                                               }
                                           }

                                           if((input1Tmp.sum()>0)&&(input2Tmp.sum()>0)){
                                           	// autocorrelation computation
                                           	CImg<> autocorrelation1=autocorrelation(input1Tmp,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth),
                                           			autocorrelation2=autocorrelation(input2Tmp,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth);
                                           	// probas
                                           	CImg<> intersection(randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,1,0);
                                           	cimg_forXYZ(intersection,x,y,z){
                                           		intersection(x,y,z) = input1Tmp(x,y,z)*input2Tmp(x,y,z);
                                           	}
                                           	double p1=input1Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth),p2=input2Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth),p12=intersection.sum()/double(randomizationWindowWidth*randomizationWindowHeight*randomizationWindowDepth);
                                           	// centering covariances
                                           	cimg_forXYZ(autocorrelation1,x,y,z){
                                           		autocorrelation1(x,y,z) -= pow(p1,2.);
                                           		autocorrelation2(x,y,z) -= pow(p2,2.);
                                           	}
                                           	// remove coefficients corresponding to noise
                                           	int radiusForTruncation=0,depthForTruncation=0;
                                           	if(truncationValue>0.){
                                           		radiusAndDepthTruncations(autocorrelation1,autocorrelation2,truncationValue,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,radiusForTruncation,depthForTruncation);
                                           	}
                                           	CImg<double> denominator(1,3,1,1,0.);
                                           	// with truncation
                                           	if(radiusForTruncation==0){
                                           		// computes S1, S2 and S3
                                           		if(wholeDenominator){
                                           			computeDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,0,0);
                                           			// outputs
                                           			// stat test whole denominator
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0) = (p12-p1*p2)/sqrt(denominator(0,0)-denominator(0,1)+denominator(0,2));
                                           			// stat test S1
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                           			// two-sided p-values whole denominator
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0))));
                                           			// two-sided p-values S1
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1))));
                                           			// one-sided p-values whole denominator
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0)));
                                           			// one-sided p-values S1
                                           			currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1)));
                                           			// normalized correlation
                                           			currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           		}
                                           		// only computes S1
                                           		else{
                                           			// S1 computation
                                           			computeS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,0,radiusForTruncation,depthForTruncation,0);
                                           			// stat test S1
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                           			// two-sided p-values S1
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1))));
                                           			// one-sided p-values S1
                                           			currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1)));
                                           			// normalized correlation
                                           			currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           		}
                                           	}
                                           	else{
                                           		// computes S1, S2 and S3
                                           		if(wholeDenominator){
                                           			// denominator computation
                                           			computeTruncatedDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,0,radiusForTruncation,depthForTruncation,0);
                                           			// stat test whole denominator
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0) = (p12-p1*p2)/sqrt(denominator(0,0)-denominator(0,1)+denominator(0,2));
                                           			// stat test S1
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                           			// two-sided p-values whole denominator
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0))));
                                           			// two-sided p-values S1
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1))));
                                           			// one-sided p-values whole denominator
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,0)));
                                           			// one-sided p-values S1
                                           			currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1)));
                                           			// normalized correlation
                                           			currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           		}
                                           		// only compute S1
                                           		else{
                                           			// S1 computation
                                           			computeTruncatedS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,randomizationWindowDepth,0,radiusForTruncation,depthForTruncation,0);
                                           			// stat test S1
                                           			currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                           			// two-sided p-values S1
                                           			//currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1))));
                                           			// one-sided p-values S1
                                           			currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2,1)));
                                           			// normalized correlation
                                           			currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           		}
                                           	}
                                           	currentDisplay(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2);
                                           }
                                           else{
                                           	currentDisplay(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,z+randomizationWindowDepth/2) = 1;
                                           }
                                       }
                                   }
                               }
                           }
                           else{
                               for(int y=0;y<(input1[0].height()-randomizationWindowHeight);y+=windowHeightStep){
                                   for(int x=0;x<(input1[0].width()-randomizationWindowWidth);x+=windowWidthStep){
                                       CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),1,1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),1,1,0);
                                       for(int j=0;j<randomizationWindowHeight;j++){
                                           for(int i=0;i<randomizationWindowWidth;i++){
                                               input1Tmp(i,j) = input1[t](x+i,y+j);
                                               input2Tmp(i,j) = input2[t](x+i,y+j);
                                           }
                                       }
                                       // autocorrelation computation
                                       CImg<> autocorrelation1=autocorrelation(input1Tmp,randomizationWindowWidth,randomizationWindowHeight,1),
                                       autocorrelation2=autocorrelation(input2Tmp,randomizationWindowWidth,randomizationWindowHeight,1);
                                       // probas
                                       CImg<> intersection(randomizationWindowWidth,randomizationWindowHeight,1,1,0);
                                       cimg_forXY(intersection,x,y){
                                           intersection(x,y) = input1Tmp(x,y)*input2Tmp(x,y);
                                       }
                                       double p1=input1Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight),p2=input2Tmp.sum()/double(randomizationWindowWidth*randomizationWindowHeight),p12=intersection.sum()/double(randomizationWindowWidth*randomizationWindowHeight);
                                       // centering covariances
                                       cimg_forXY(autocorrelation1,x,y){
                                           autocorrelation1(x,y) -= pow(p1,2.);
                                           autocorrelation2(x,y) -= pow(p2,2.);
                                       }
                                       // remove coefficients corresponding to noise
                                       int radiusForTruncation=0,depthForTruncation=0;
                                       if(truncationValue>0.){
                                           radiusAndDepthTruncations(autocorrelation1,autocorrelation2,truncationValue,randomizationWindowWidth,randomizationWindowHeight,1,radiusForTruncation,depthForTruncation);
                                       }
                                       
                                       CImg<double> denominator(1,3,1,1,0.);
                                       // with truncation
                                       if(radiusForTruncation==0){
                                           // computes S1, S2 and S3
                                           if(wholeDenominator){
                                               computeDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,1,0,0);
                                               // outputs
                                               // stat test whole denominator
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight,0) = (p12-p1*p2)/sqrt(denominator(0,0)-denominator(0,1)+denominator(0,2));
                                               // stat test S1
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                               // two-sided p-values whole denominator
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0))));
                                               // two-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1))));
                                               // one-sided p-values whole denominator
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0)));
                                               // one-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,3) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1)));
                                               // normalized correlation
                                               currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           }
                                           // only computes S1
                                           else{
                                               // S1 computation
                                               computeS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,1,0,radiusForTruncation,depthForTruncation,0);
                                               // stat test S1
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                               // two-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1))));
                                               // one-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,3) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1)));
                                               // normalized correlation
                                               currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           }
                                       }
                                       else{
                                           // computes S1, S2 and S3
                                           if(wholeDenominator){
                                               // denominator computation
                                               computeTruncatedDenominator(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,1,0,radiusForTruncation,depthForTruncation,0);
                                               // stat test whole denominator
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0) = (p12-p1*p2)/sqrt(denominator(0,0)-denominator(0,1)+denominator(0,2));
                                               // stat test S1
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                               // two-sided p-values whole denominator
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0))));
                                               // two-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1))));
                                               // one-sided p-values whole denominator
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,2) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,0)));
                                               // one-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,3) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1)));
                                               // normalized correlation
                                               currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           }
                                           // only compute S1
                                           else{
                                               // S1 computation
                                               computeTruncatedS1(denominator,autocorrelation1,autocorrelation2,randomizationWindowWidth,randomizationWindowHeight,1,0,radiusForTruncation,depthForTruncation,0);
                                               // stat test S1
                                               currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = (p12-p1*p2)/sqrt(denominator(0,0));
                                               // two-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1) = 2.*(1.-CND(abs(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1))));
                                               // one-sided p-values S1
                                               currentPvalues(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,3) = (1.-CND(currentStatisticalTest(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2,1)));
                                               // normalized correlation
                                               currentCorrelation(x+randomizationWindowWidth/2,y+randomizationWindowHeight/2) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                           }
                                       }
                                   }
                               }
                           }
                           /*statisticalTest.push_back(currentStatisticalTest);
                           correlation.push_back(currentStatisticalTest.get_blur(5.,false,true));
                           pvalues.push_back(currentDisplay.get_blur(5.,false,true));//(currentPvalues);*/
                           cimg_forXY(currentStatisticalTest,x,y){
                           	currentStatisticalTest2(x,y) = currentStatisticalTest(x,y,14);
                           }
                           currentStatisticalTest.blur(5.,false,true);
                           currentDisplay.get_blur(5.,false,true);
                           cimg_forXY(currentStatisticalTest,x,y){
                           	currentCorrelation2(x,y) = currentStatisticalTest(x,y,14);
                           	currentPvalues2(x,y) = currentDisplay(x,y,14);
                           }
                           statisticalTest.push_back(currentStatisticalTest2);
                           correlation.push_back(currentCorrelation2);
                           pvalues.push_back(currentPvalues2);
                           //currentDisplay.blur(5.,false,true);
                           //statisticalTestDisplay.push_back(currentDisplay.get_blur(5.,false,true));
                       }
                       pvalues.save(pvaluesFilename.c_str());
                       correlation.save(correlationFilename.c_str());
                       statisticalTest.save(statisticalTestFilename.c_str());
                   }
                   else{
                       if(!saveImagesAsExamples){
                           omp_set_num_threads(nbThreads);
                           
                           //CImg<> input1(dim,dim,1,1,0.),input2(dim,dim,1,1,0.);
                           //cimg_forXY(input1,x,y){
                           //input1(x,y) = round(cimg::rand());
                           //input2(x,y) = round(cimg::rand());
                           //}
                           
                           // output
                           CImg<double> statisticalTest(input1.size(),2,1,1,0.),correlation(input1.size(),1,1,1,0.),pvalues(input1.size(),4,1,1,0.),denominator(input1.size(),3,1,1,0.),truncation(input1.size(),2,1,1,0.);
                           
                           // new dimensions for FFT
                           int powerX=2,powerY=2,powerZ=2;
                           while(2*input1[0].width()>pow(2,double(powerX))){powerX++;}
                           while(2*input1[0].height()>pow(2,double(powerY))){powerY++;}
                           if(input1[0].depth()>1){
                               while(2*input1[0].depth()>pow(2,double(powerZ))){powerZ++;}
                           }
                           else{
                               powerZ=0;
                           }
                           
                           // statistical test estimation over time
#pragma omp parallel for
                           for(int t=0;t<input1.size();t++){
                               CImg<> input1Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0),input2Tmp(int(pow(2,double(powerX))),int(pow(2,double(powerY))),int(pow(2,double(powerZ))),1,0);
                               cimg_forXYZ(input1[t],x,y,z){
                                   input1Tmp(x,y,z) = input1[t](x,y,z);
                                   input2Tmp(x,y,z) = input2[t](x,y,z);
                               }
                               // autocorrelation computation
                               CImg<> autocorrelation1=autocorrelation(input1Tmp,input1[t].width(),input1[t].height(),input1[t].depth()),
                               autocorrelation2=autocorrelation(input2Tmp,input1[t].width(),input1[t].height(),input1[t].depth());

                               // probas
                               CImg<> intersection(input1[t].width(),input1[t].height(),input1[t].depth(),1,0);
                               cimg_forXYZ(intersection,x,y,z){
                                   intersection(x,y,z) = input1Tmp(x,y,z)*input2Tmp(x,y,z);
                               }
                               double p1=input1Tmp.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth()),p2=input2Tmp.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth()),p12=intersection.sum()/double(input1[t].width()*input1[t].height()*input1[t].depth());
                               cout << "p1: " << p1 << "   p2: " << p2 << "   p12: " << p12 << endl;
                               // centering covariances
                               cimg_forXYZ(autocorrelation1,x,y,z){
                                   autocorrelation1(x,y,z) -= pow(p1,2.);
                                   autocorrelation2(x,y,z) -= pow(p2,2.);
                               }
                               // remove coefficients corresponding to noise
                               int radiusForTruncation=0,depthForTruncation=0;
                               if(truncationValue>0.){
                                   radiusAndDepthTruncations(autocorrelation1,autocorrelation2,truncationValue,input1[t].width(),input1[t].height(),input1[t].depth(),radiusForTruncation,depthForTruncation);
                                   truncation(t,0) = radiusForTruncation;
                                   truncation(t,1) = depthForTruncation;
                                   cout << "radius for truncation: " << radiusForTruncation << "   depth for truncation: " << depthForTruncation << endl;
                               }
                               // with truncation
                               if(radiusForTruncation==0){
                                   // computes S1, S2 and S3
                                   if(wholeDenominator){
                                       computeDenominator(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,0);
                                       // outputs
                                       // stat test whole denominator
                                       statisticalTest(t,0) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                                       // stat test S1
                                       statisticalTest(t,1) = (p12-p1*p2)/sqrt(denominator(t,0));
                                       // two-sided p-values whole denominator
                                       pvalues(t,0) = (1.-CND(-statisticalTest(t,0)));
                                       // two-sided p-values S1
                                       pvalues(t,1) = (1.-CND(-statisticalTest(t,1)));
                                       // one-sided p-values whole denominator
                                       pvalues(t,2) = (1.-CND(statisticalTest(t,0)));
                                       // one-sided p-values S1
                                       pvalues(t,3) = (1.-CND(statisticalTest(t,1)));
                                       // normalized correlation
                                       correlation(t) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                   }
                                   // only computes S1
                                   else{
                                       // S1 computation
                                       computeS1(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,0);
                                       // stat test S1
                                       statisticalTest(t,1) = (p12-p1*p2)/sqrt(denominator(t,0));
                                       // two-sided p-values S1
                                       pvalues(t,1) = (1.-CND(-statisticalTest(t,1)));
                                       // one-sided p-values S1
                                       pvalues(t,3) = (1.-CND(statisticalTest(t,1)));
                                       // normalized correlation
                                       correlation(t) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                   }
                               }
                               else{
                                   // computes S1, S2 and S3
                                   if(wholeDenominator){
                                       // denominator computation
                                       computeTruncatedDenominator(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,0);
                                       /*double S1=0.,S2=0.,S3_1=0.,S3_2=0.;
                                        if(input1[t].depth()==1){
                                        for(int i2=0;i2<input1[t].height();i2++){
                                        for(int i1=0;i1<input1[t].width();i1++){
                                        double S2_1=0.,S2_2=0.;
                                        for(int y=max(input1[t].height()-radiusForTruncation,i2+1);y<=min(input1[t].height()+radiusForTruncation,i2+input1[t].height());y++){
                                        for(int x=max(input1[t].width()-radiusForTruncation,i1+1);x<=min(input1[t].width()+radiusForTruncation,i1+input1[t].width());x++){
                                        S2_1 += autocorrelation1(x,y);
                                        S2_2 += autocorrelation2(x,y);
                                        }
                                        }
                                        S2 += S2_1*S2_2;
                                        }
                                        }
                                        }
                                        else{
                                        for(int i3=0;i3<input1[t].depth();i3++){
                                        for(int i2=0;i2<input1[t].height();i2++){
                                        for(int i1=0;i1<input1[t].width();i1++){
                                        double S2_1=0.,S2_2=0.;
                                        for(int z=max(input1[t].depth()-radiusForTruncation,i3+1);z<=min(input1[t].depth()+radiusForTruncation,i3+input1[t].depth());z++){
                                        for(int y=max(input1[t].height()-radiusForTruncation,i2+1);y<=min(input1[t].height()+radiusForTruncation,i2+input1[t].height());y++){
                                        for(int x=max(input1[t].width()-radiusForTruncation,i1+1);x<=min(input1[t].width()+radiusForTruncation,i1+input1[t].width());x++){
                                        S2_1 += autocorrelation1(x,y,z);
                                        S2_2 += autocorrelation2(x,y,z);
                                        }
                                        }
                                        }
                                        S2 += S2_1*S2_2;
                                        }
                                        }
                                        }
                                        }
                                        CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                                        if(input1[t].depth()==1){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                                        }
                                        }
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                        S3_1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                        S3_2 += autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                        }
                                        }
                                        denominator(t,0) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                                        denominator(t,1) = 2.*S2/pow(double(input1[t].width()*input1[t].height()),3.);
                                        denominator(t,2) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()),4.);
                                        }
                                        else{
                                        for(int z=-radiusForTruncation;z<=radiusForTruncation;z++){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                                        }
                                        }
                                        }
                                        for(int z=-radiusForTruncation;z<=radiusForTruncation;z++){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                        S3_1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                        S3_2 += autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                        }
                                        }
                                        }
                                        denominator(t,0) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                                        denominator(t,1) = 2.*S2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),3.);
                                        denominator(t,2) = S3_1*S3_2/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),4.);
                                        }
                                        
                                        statisticalTest(t,0) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                                        statisticalTest(t,1) = (p12-p1*p2)/sqrt(denominator(t,0));
                                        
                                        pvalues(t,0) = 2.*(1.-CND(abs(statisticalTest(t,0))));
                                        pvalues(t,1) = 2.*(1.-CND(abs(statisticalTest(t,1))));
                                        
                                        pvalues(t,2) = (1.-CND(statisticalTest(t,0)));
                                        pvalues(t,3) = (1.-CND(statisticalTest(t,1)));
                                        
                                        correlation(t) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                        }
                                        else{
                                        // denominator computation
                                        double S1=0.;
                                        CImg<> weightedMatrix(autocorrelation1.width(),autocorrelation1.height(),autocorrelation1.depth(),1,0.);
                                        if(input1[t].depth()==1){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        weightedMatrix(x+input1[t].width(),y+input1[t].height()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y));
                                        }
                                        }
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height())*autocorrelation2(x+input1[t].width(),y+input1[t].height())*weightedMatrix(x+input1[t].width(),y+input1[t].height());
                                        }
                                        }
                                        denominator(t,0) = S1/pow(double(input1[t].width()*input1[t].height()),2.);
                                        }
                                        else{
                                        for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth()) = (input1[t].width()-abs(x))*(input1[t].height()-abs(y))*(input1[t].depth()-abs(z));
                                        }
                                        }
                                        }
                                        for(int z=-depthForTruncation;z<=depthForTruncation;z++){
                                        for(int y=-radiusForTruncation;y<=radiusForTruncation;y++){
                                        for(int x=-radiusForTruncation;x<=radiusForTruncation;x++){
                                        S1 += autocorrelation1(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*autocorrelation2(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth())*weightedMatrix(x+input1[t].width(),y+input1[t].height(),z+input1[t].depth());
                                        }
                                        }
                                        }
                                        denominator(t,0) = S1/pow(double(input1[t].width()*input1[t].height()*input1[t].depth()),2.);
                                        }*/
                                       // stat test whole denominator
                                       statisticalTest(t,0) = (p12-p1*p2)/sqrt(denominator(t,0)-denominator(t,1)+denominator(t,2));
                                       // stat test S1
                                       statisticalTest(t,1) = (p12-p1*p2)/sqrt(denominator(t,0));
                                       // two-sided p-values whole denominator
                                       pvalues(t,0) = (1.-CND(-statisticalTest(t,0)));
                                       // two-sided p-values S1
                                       pvalues(t,1) = (1.-CND(-statisticalTest(t,1)));
                                       // one-sided p-values whole denominator
                                       pvalues(t,2) = (1.-CND(statisticalTest(t,0)));
                                       // one-sided p-values S1
                                       pvalues(t,3) = (1.-CND(statisticalTest(t,1)));
                                       // normalized correlation
                                       correlation(t) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                   }
                                   // only compute S1
                                   else{
                                       // S1 computation
                                       computeTruncatedS1(denominator,autocorrelation1,autocorrelation2,input1[t].width(),input1[t].height(),input1[t].depth(),t,radiusForTruncation,depthForTruncation,0);
                                       // stat test S1
                                       statisticalTest(t,1) = (p12-p1*p2)/sqrt(denominator(t,0));
                                       cout << "numerator: " << (p12-p1*p2) << "   sqrt denominator: " << sqrt(denominator(t,0)) << endl;
                                       // two-sided p-values S1
                                       pvalues(t,1) = (1.-CND(-statisticalTest(t,1)));
                                       // one-sided p-values S1
                                       pvalues(t,3) = (1.-CND(statisticalTest(t,1)));
                                       // normalized correlation
                                       correlation(t) = (p12-p1*p2)/sqrt(p1*(1.-p1)*p2*(1.-p2));
                                   }
                               }
                           }
                           
                           // post-treatment
                           //vector<double> s1,s2,s3,denom,stat1,stat2,pval1,pval2,pval3,pval4;
                           vector<double> s,pval;
                           vector<int> distributionFunction;
                           postTreatment(input1.size(),pval,distributionFunction,pvalues);
                           /*int nbRejects1=0,nbRejects2=0,nbRejects3=0,nbRejects4=0;
                            for(unsigned int t=0;t<input1.size();t++){
                            s1.push_back(denominator(t,0));
                            s2.push_back(denominator(t,1));
                            s3.push_back(denominator(t,2));
                            denom.push_back(denominator(t,0)-denominator(t,1)+denominator(t,2));
                            stat1.push_back(statisticalTest(t,0));
                            stat2.push_back(statisticalTest(t,1));
                            pval1.push_back(pvalues(t,0));
                            if(pvalues(t,0)<0.05){nbRejects1++;}
                            pval2.push_back(pvalues(t,1));
                            if(pvalues(t,1)<0.05){nbRejects2++;}
                            pval3.push_back(pvalues(t,2));
                            if(pvalues(t,2)<0.05){nbRejects3++;}
                            pval4.push_back(pvalues(t,3));
                            if(pvalues(t,3)<0.05){nbRejects4++;}
                            }
                            sort(s1.begin(),s1.end());
                            sort(s2.begin(),s2.end());
                            sort(s3.begin(),s3.end());
                            sort(denom.begin(),denom.end());
                            sort(stat1.begin(),stat1.end());
                            sort(stat2.begin(),stat2.end());
                            sort(pval1.begin(),pval1.end());
                            sort(pval2.begin(),pval2.end());
                            sort(pval3.begin(),pval3.end());
                            sort(pval4.begin(),pval4.end());
                            
                            distributionFunction.push_back(0);
                            for(unsigned int u=1;u<pval3.size();u++){
                            int nbPvals = distributionFunction[distributionFunction.size()-1];
                            if(distributionFunction[distributionFunction.size()-1]<pval3.size()){
                            for(unsigned int v=distributionFunction[distributionFunction.size()-1];v<pval3.size();v++){
                            if(pval3[v]<(double(u)/double(pval3.size()))){
                            nbPvals++;
                            }
                            }
                            distributionFunction.push_back(nbPvals);
                            }
                            else{
                            distributionFunction.push_back(pval3.size());
                            }
                            }*/
                           
                           // indicator output
                           cout << "p-value: min: " << pval[0] << "   lower quartile: " << pval[pval.size()/4] << "   median: " << pval[pval.size()/2] << "   upper quartile: " << pval[3*pval.size()/4] << "   max: " << pval[pval.size()-1] << endl;
                           cout << "test statistic: " << statisticalTest(0,1) << endl;
                           // output files
                           string pvaluesFilename,sortedPvaluesFilename,pValueDistributionFunctionFilename,correlationFilename,truncationFilename;//denominatorFile,statFile,pvalFile;
                           if(outputFilename!=NULL){
                               pvaluesFilename = outputFilename;
                               pvaluesFilename += "_pvalues.txt";
                               sortedPvaluesFilename = outputFilename;
                               sortedPvaluesFilename += "_sortedPvalues.txt";
                               pValueDistributionFunctionFilename = outputFilename;
                               pValueDistributionFunctionFilename += "_pValueDistributionFunction.txt";
                               correlationFilename = outputFilename;
                               correlationFilename += "_correlation.txt";
                               truncationFilename = outputFilename;
                               truncationFilename += "_truncatedRadiusAndDepth.txt";
                               //denominatorFile = outputFilename;
                               //denominatorFile += "_denominator.txt";
                               //statFile = outputFilename;
                               //statFile += "_stat.txt";
                               //pvalFile = outputFilename;
                               //pvalFile += "_pvalues.txt";
                           }
                           else{
                               pvaluesFilename = "pvalues.txt";
                               sortedPvaluesFilename = "sortedPvalues.txt";
                               pValueDistributionFunctionFilename = "pValueDistributionFunction.txt";
                               correlationFilename = "correlation.txt";
                               truncationFilename = "truncatedRadiusAndDepth.txt";
                               //denominatorFile = "denominator.txt";
                               //statFile = "stat.txt";
                               //pvalFile = "pvalues.txt";
                           }
                           ofstream pValues(pvaluesFilename.c_str()),sortedPvalues(sortedPvaluesFilename.c_str()),pValueDistributionFunction(pValueDistributionFunctionFilename.c_str()),correlationOutput(correlationFilename.c_str()),truncationOutput(truncationFilename.c_str());//outputDenominator(denominatorFile.c_str()),outputStat(statFile.c_str()),outputPval(pvalFile.c_str());
                           for(unsigned int i=0;i<pval.size();i++){
                               pValues << pvalues(i,3) << endl;
                               sortedPvalues << pval[i] << endl;
                               pValueDistributionFunction << distributionFunction[i] << endl;
                               correlationOutput << correlation(i) << endl;
                               if(input1[0].depth()==1){
                                   truncationOutput << truncation(i,0) << endl;
                               }
                               else{
                                   truncationOutput << truncation(i,0) << '\t' << truncation(i,1) << endl;
                               }
                           }
                           
                       }
                       else{
                           CImg<> input1(inputImage1),input2(inputImage2),output(input1.width(),input1.height(),1,3,0);
                           //input1.threshold(1);input2.threshold(1);
                           cimg_forXY(input1,x,y){
                               output(x,y,0) = input1(x,y,0);
                               output(x,y,1) = input2(x,y,0);
                           }
//                            // output.normalize(0,255).save("output.png").display();
                        }
                    }
                }
            }
        }
        
    }
    cimg::toc();
    
    return(0);
    
}
