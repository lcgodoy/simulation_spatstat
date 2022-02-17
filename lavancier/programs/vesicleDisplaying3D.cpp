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

#include "vesicleDisplaying3D.h"

inline unsigned int min (unsigned int val1, unsigned int val2) {return (val1<val2)?val1:val2;}
inline unsigned int max (unsigned int val1, unsigned int val2) {return (val1>val2)?val1:val2;}


// first vector level : frame
// second vector level : Barycentre3D coordinates <x,y>
vector< particle3D > getObjectMassCenterCoordinates3D (CImg<int> &detection) {

	vector< particle3D > result;

	// classify connected components in two pass
	CImg<int> currentDetection = detection.get_label(true);

	CImg<> detection_Stats=currentDetection.get_stats();
	for(int i=1;i<=detection_Stats(1);i++){
		particle3D part;
		part.t=0;
		cimg_forXYZ(currentDetection,x,y,z){
			if(currentDetection(x,y,z)==i){
				features3D feats;
				feats.push_back(double(x));
				feats.push_back(double(y));
				feats.push_back(double(z));
				part.vol.push_back(feats);
			}
		}
		if(part.vol.size()>0){
			part.computeMassCenter();
			result.push_back(part);
		}
	}

	return result;
}

