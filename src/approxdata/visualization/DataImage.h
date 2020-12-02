/*
 * ApproxData Library
 * Copyright (c) 2018 Michael Shekelyan <michael.shekelyan@gmail.com>
 */

/*
Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions 
of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef SHEKELYAN_SELIMAGE_H
#define SHEKELYAN_SELIMAGE_H

#define cimg_use_png
#include <approxdata/external/CImg.h>

using namespace cimg_library;

class DataImage{

private:
	
	Point temp1;
	Point temp2;
	
	Box dataBox;

public:

	CImg<float> img;
	
	int w;
	int h;
	
	double xMul;
	double xAdd;
	
	double yMul;
	double yAdd;
	
	int xMax;
	int yMax;
	
	inline Box getDataBox(){
	
		return dataBox;
	}
	
	inline DataImage(int w, int h) : w(w), h(h), img(w,h, 1, 3), temp1(2), temp2(2){
		
		img = 255;
		
		Point datamin(2);
		Point datamax(2);
		
		datamin[0] = 0;
		datamin[1] = 0;
		
		datamax[0] = 1;
		datamax[1] = 1;
		
		xMax = w-1;
		yMax = h-1;
		
		Box b(datamin, datamax);
		
		setDataBox(b);
	}
	
	inline void resize(long width, long height){

		w = width;
		h = height;
	
		img.resize(w,h,-100, -100, 2);

	}
	
	inline void setDataBox(double x1, double y1, double x2, double y2){
	
		Box box(x1, y1, x2, y2);
		
		dataBox = box;
		xMul = 1.0/( box.max[0]-box.min[0] );
		xAdd = -box.min[0];
		
		yMul = 1.0/( box.max[1]-box.min[1] );
		yAdd = -box.min[1];
	}
	
	inline void setDataBox(Box box){
	
		dataBox = box;
		xMul = 1.0/( box.max[0]-box.min[0] );
		xAdd = -box.min[0];
		
		yMul = 1.0/( box.max[1]-box.min[1] );
		yAdd = -box.min[1];
	}

	inline void save(const string& file){
	
		img.save_png(file.c_str());
	}
	
	inline int display(){
	
		img.display();
		
		return 0;
	}
	
	inline int display(const string s){
	
		img.display(s.c_str());
		
		return 0;
	}
	
	inline void toImg(Point& in, Point& out){
	
		out.clear();
		
		const double normX = (in[0]+xAdd)*xMul;
		const double normY = (in[1]+yAdd)*yMul;
		
		out.push_back( normX*xMax );
		out.push_back( yMax-(normY*yMax) );
	}
	
	/*
	inline void drawLine(Point& p1, Point& p2, const float c[3]){
		
		toImg(p1, temp1);
		toImg(p2, temp2);
		
		img.draw_line(temp1[0], temp1[1], temp2[0], temp2[1], c.data());
	}*/
	
	inline void draw(int x, int y, array<float,3>& c, double a){
	
		if (x < 0 || x > xMax)
			return;
			
		if (y < 0 || y > yMax)
			return;
	
		a = UtilMath::makeBetween<double>(0,1, a);
	
		
	
		img(x,y,0) = img(x,y,0)*(1.0-a)+c[0]*a;
		img(x,y,1) = img(x,y,1)*(1.0-a)+c[1]*a;
		img(x,y,2) = img(x,y,2)*(1.0-a)+c[2]*a;
	}
	
	
	inline void drawPoint(Point& p, array<float,3>& c, double a){
		
		
		
		toImg(p, temp1);
		
		if (!dataBox.contains(p))
			return;
		
		draw(temp1[0], temp1[1], c, 50000*a);
	}
	
	
	inline void drawRect(Box& b, array<float,3>& c, double a){
		
		if (dataBox.intersection(b) == SpatialRelation::NO_INTERSECTION)
			return;
		
		
		toImg(b.min, temp1);
		toImg(b.max, temp2);
		
		const int x1 = UtilMath::makeBetween<int>(0, w-1, temp1[0]);
		const int y1 = UtilMath::makeBetween<int>(0, h-1, temp1[1]);
		
		const int x2 = UtilMath::makeBetween<int>(0, w-1, temp2[0]);
		const int y2 = UtilMath::makeBetween<int>(0, h-1, temp2[1]);
		
		for (int x = x1; x <= x2; x++){
			draw(x, y1, c, a);
			draw(x, y2, c, a);
		}
		
		for (int y = y1; y >= y2; y--){
			draw(x1, y, c, a);
			draw(x2, y, c, a);
		}
	}
	
	inline void fillRect(Box& b, array<float,3>& c){
	
		toImg(b.min, temp1);
		toImg(b.max, temp2);
		
		const double pixelVolume = UtilMath::absVal(temp2[0]-temp1[0])*UtilMath::absVal(temp2[1]-temp1[1]);
		
		const double pixelWeight = 50000*b.count/pixelVolume;
		
		const int x1 = UtilMath::makeBetween<int>(0, w-1, temp1[0]);
		const int y1 = UtilMath::makeBetween<int>(0, h-1, temp1[1]);
		
		const int x2 = UtilMath::makeBetween<int>(0, w-1, temp2[0]-1);
		const int y2 = UtilMath::makeBetween<int>(0, h-1, temp2[1]+1);
		
		for (int x = x1; x <= x2; x++){
		
			for (int y = y1; y >= y2; y--){
			
				draw(x, y, c, pixelWeight);
			}
			
		}
		
	}
	
	inline void fillRect(Box& b, array<float,3>& c, double a){
		
		if (dataBox.intersection(b) == SpatialRelation::NO_INTERSECTION)
			return;
		
		toImg(b.min, temp1);
		toImg(b.max, temp2);
		
		const int x1 = UtilMath::makeBetween<int>(0, w-1, temp1[0]);
		const int y1 = UtilMath::makeBetween<int>(0, h-1, temp1[1]);
		
		const int x2 = UtilMath::makeBetween<int>(0, w-1, temp2[0]-1);
		const int y2 = UtilMath::makeBetween<int>(0, h-1, temp2[1]+1);
		
		for (int x = x1; x <= x2; x++){
		
			for (int y = y1; y >= y2; y--){
			
			
				draw(x, y, c, a);
			}
			
		}
		
	}
};

#endif