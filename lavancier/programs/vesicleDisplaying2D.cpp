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

#include "vesicleDisplaying2D.h"

inline unsigned int min (unsigned int val1, unsigned int val2) {return (val1<val2)?val1:val2;}
inline unsigned int max (unsigned int val1, unsigned int val2) {return (val1>val2)?val1:val2;}

void neighbor(int x,int y,CImg<int> & detect,CImg<int> & flag,vector<coord2D> & cell){

	for(int v=-1;v<=1;v++){
		for(int u=-1;u<=1;u++){
			if(((x+u)>=0)&&((x+u)<detect.width())&&((y+v)>=0)&&((y+v)<detect.height())){
				if((flag(x+u,y+v)==0)&&(detect(x+u,y+v)==1)){
					flag(x+u,y+v) = 1;
					coord2D point(x+u,y+v);
					cell.push_back(point);
					neighbor(x+u,y+v,detect,flag,cell);
				}
			}
		}
	}
}

// first vector level : frame
// second vector level : barycentre coordinates <x,y>
vector< particle2D > getObjectMassCenterCoordinates2D (CImg<int> &detection) {

	vector< particle2D > result;

	// classify connected components in two pass
	for (int t = 0; t < detection.depth(); t++) { // for each frame
		CImg<int> detection_t(detection.width(),detection.height(),detection.spectrum(),1,0);
		cimg_forXY(detection,x,y){
			if(detection(x,y,t)>0.){
				detection_t(x,y) = 1;
			}
		}
		detection_t.label(true);
		cimg_forXY(detection,x,y){
			if((detection_t(x,y)>0)&&(detection(x,y,t)==0)){
				detection_t(x,y) = 0;
			}
		}

		CImg<> detection_tStats=detection_t.get_stats();
		for(int i=1;i<=detection_tStats(1);i++){
			particle2D part;
			part.t=t;
			cimg_forXY(detection_t,x,y){
				if(detection_t(x,y)==i){
					features2D feats;
					feats.push_back(double(x));
					feats.push_back(double(y));
					part.vol.push_back(feats);
				}
			}
			if(part.vol.size()>0){
				part.computeMassCenter();
				result.push_back(part);
			}
		}
	}

	return result;

}

