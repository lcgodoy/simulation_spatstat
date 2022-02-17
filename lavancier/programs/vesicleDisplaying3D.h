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

#ifndef croixvesicle3D_h
#define croixvesicle3D_h

#include "CImg.h"
#include <vector>
#include <set>
#include <iostream>
#include <map>

using namespace cimg_library;
using namespace std;

typedef vector<int> coord3D;
typedef vector<coord3D> pair_list3D;
typedef vector<double> features3D;
typedef vector<features3D> particleVolume;


class particle3D {
public:
	particleVolume vol;
	int t;
	double xCoord;
	double yCoord;
	double zCoord;
	double intensity;
	double minIntensity;
	double maxIntensity;

	particle3D() {
		this->t = 0;
		this->xCoord = 0.;
		this->yCoord = 0.;
		this->zCoord = 0.;
		this->intensity = 0.;
		this->minIntensity = 0.;
		this->maxIntensity = 0.;
	}
	particle3D(particleVolume volume,int t=0, double x=0., double y=0., double z=0., double intensity=0., double minIntensity=0., double maxIntensity=0.) {
		this->vol = volume;
		this-> t = t;
		this->xCoord = x;
		this->yCoord = y;
		this->zCoord = z;
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
			this->zCoord += this->vol[i][2];
			this->intensity += this->vol[i][3];
			if(this->vol[i][3]<this->minIntensity){this->minIntensity=this->vol[i][3];}
			if(this->vol[i][3]>this->maxIntensity){this->maxIntensity=this->vol[i][3];}
		}
		this->xCoord /= double(this->vol.size());
		this->yCoord /= double(this->vol.size());
		this->zCoord /= double(this->vol.size());
		this->intensity /= double(this->vol.size());
	}

};

vector< particle3D > getObjectMassCenterCoordinates3D (CImg<int> &detection);

#endif
