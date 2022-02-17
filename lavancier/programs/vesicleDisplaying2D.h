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

#ifndef croixvesicle2D_h
#define croixvesicle2D_h

#include "CImg.h"
#include <vector>
#include <set>
#include <iostream>
#include <map>

using namespace cimg_library;
using namespace std;

typedef pair<unsigned int,unsigned int> coord2D;
typedef vector<coord2D> pair_list2D;
typedef pair<double,double> coord2Ddouble;
typedef vector<coord2Ddouble> pair_list2Ddouble;
typedef vector<double> features2D;
typedef vector<features2D> particleVolume;

class particle2D {
public:
	particleVolume vol;
	int t;
	double xCoord;
	double yCoord;
	double intensity;
	double minIntensity;
	double maxIntensity;

	particle2D() {
		this->t = 0;
		this->xCoord = 0.;
		this->yCoord = 0.;
		this->intensity = 0.;
		this->minIntensity = 0.;
		this->maxIntensity = 0.;
	}
	particle2D(particleVolume volume,int t=0, double x=0., double y=0., double intensity=0., double minIntensity=0., double maxIntensity=0.) {
		this->vol = volume;
		this-> t = t;
		this->xCoord = x;
		this->yCoord = y;
		this->intensity = intensity;
		this->minIntensity = minIntensity;
		this->maxIntensity = maxIntensity;
	}
	void computeMassCenter(){
		this->minIntensity = this->vol[0][3];
		this->maxIntensity = this->vol[0][3];
		for(unsigned int i=0;i<this->vol.size();i++){
			this->xCoord += this->vol[i][0];
			this->yCoord += this->vol[i][1];
			this->intensity += this->vol[i][2];
			if(this->vol[i][2]<this->minIntensity){this->minIntensity=this->vol[i][2];}
			if(this->vol[i][2]>this->maxIntensity){this->maxIntensity=this->vol[i][2];}
		}
		this->xCoord /= double(this->vol.size());
		this->yCoord /= double(this->vol.size());
		this->intensity /= double(this->vol.size());
	}

};

vector< particle2D > getObjectMassCenterCoordinates2D (CImg<int> &detection);

#endif
