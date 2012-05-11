#include "DataGrid2D.h"
//#include "CompressionUtils.h"
#include "Common.h"
//#include <fstream>
//#include <iostream>
//using namespace std;
//#include <bzlib.h>
#include <math.h>

#ifndef WIN32
#include <stdlib.h>
#endif

/*
DataGrid2D::DataGrid2D(const string & filename, bool compress)  {
	x = 1;
	y = 1;
	z = -1;
	mode = 1;

	dimension = 2;
	ex = 0;

	loadData(filename, compress);
}
*/

DataGrid2D::DataGrid2D(int xval, int yval, int zval)  { 
	x = xval;
	y = yval;
	z = zval;
	
	origx = xval;
	origy = yval;
	origz = zval;
	dimension = 2;
	
	mode = 0;
	
	data = new row[x];
	
	for (int k =0; k < x; k++)
        {
		data[k].col = new double[y];
                for (int j=0;j<y;j++)
                  data[k].col[j]=0;
//		bzero(data[k].col,sizeof(double)*y);
	}
	
}

DataGrid2D::~DataGrid2D() {
	for (int k =0; k < x; k++) {
		delete [] data[k].col;
	}
	delete [] data;
}

void DataGrid2D::setData(int xval, int yval, int zval, double value) {
	if (data == NULL || xval > x || yval > y) {
		return;
	}

	data[xval].col[yval] = value;
}

double DataGrid2D::getData(int xval, int yval, int zval) {
	if (xval >= x || yval >= y) {
//		cout << "REQUESTING BAD DATA: " << xval << "," << yval <<
//		" and max size is: " << x << "," <<yval << endl;
	}

	return data[xval].col[yval];
}

/*
void DataGrid2D::loadCompressedData(ifstream & input) {

	double * result = CompressionUtils::inflate(this,input);

	int index = 0;

	for (int k =0; k < x; k++) {
		data[k].col = new double[y];
		bzero(data[k].col,sizeof(double)*y);
		for (int j = 0; j < y; j++) {
			double val = result[index];
			data[k].col[j] = val;
			index++;
		}
	}

}
*/

void DataGrid2D::output() {
	for (int k =0; k < x; k++) {
		for (int j = 0; j < y; j++) {
//			cout << k << "," << j << " = " << this->getData(k,j,0) << " ";
		}
//		cout << endl;
	}
	
//	cout << " sized: " << x << " " << y << endl;
}


void DataGrid2D::resize(int newx, int newy, int newz, bool copy) {
	row * newdata = new row[newx];
	for (int k =0; k < newx; k++) {
		newdata[k].col = new double[newy];
                for (int j=0;j<newy;j++)
                   newdata[k].col[j]=0;
//		bzero(newdata[k].col,sizeof(double)*y);
	}

	if (copy) {
		for (int k =0; k < x; k++) {
			for (int j = 0; j < y; j++) {
				double val = this->getData(k,j,0);
				newdata[k].col[j] = val;
			}
		}
	}

        /* delete the old data */
	for (int k =0; k < x; k++)
          delete [] data[k].col;
	delete [] data;

	data = newdata;
	x = newx;
	y = newy;
}

void DataGrid2D::copyTo(DataGrid * newData) {

	for (int k =0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			double val = this->getData(k,j,0);
			newData->setData(k,j,0,val);
		}	
	}
}

/*
void DataGrid2D::loadData(const string & filename, bool compress)
{
	
	ifstream input(filename.c_str(), ios::in | ios::binary);
	
	input.seekg(0);

	int dimension, xsize, ysize;
	
	// if we are decompressin	
	if (!compress) {
		
		int originalx, originaly;
		
		input.read((char*)&dimension, sizeof (int));
		input.read((char*)&xsize, sizeof (int));
		input.read((char*)&ysize, sizeof (int));
		input.read((char*)&originalx, sizeof (int));	
		input.read((char*)&originaly, sizeof (int));	
		
		this->origx =  originalx;
		this->origy =  originaly;
		
		if (verbose) cout << "Data has a current (X,Y): (" << xsize << "," << ysize <<
			") with dimension " << dimension << " but is originally sized (X,Y) : (" << 
			origx << "," << origy << ") " << endl;
	}
	else {
		double xs, ys;
		input.read((char*)&xs, sizeof (double));
		input.read((char*)&ys, sizeof (double));

		xsize = (int) xs;
		ysize = (int) ys;
	
		this->origx = xsize;
		this->origy = ysize;
	
		mode = 0;
	}
	
	if (x <= 0 || y <= 0) {
		cout << "The X and Y values are not correct (" << x << "," << y << "). " <<
		"Perhaps, this is not the correct dimension " << endl;	
	}
	
	x =  xsize;
	y =  ysize;
	
	originalFileSize = x*y*sizeof(double) + sizeof(double)*2;
	
	if (x != y) {
		cout << "You should use a square matrix. Results may not be accurate" << endl;
	}
	
	data = new row[x];
	
	if (!compress) {
		this->loadLevels(input);
		this->loadCompressedData(input);
		return;
	}
	
	for (int k =0; k < x; k++) {
		data[k].col = new double[y];
		bzero(data[k].col,sizeof(double)*y);
	}

	for (int k = 0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			input.read((char*)&data[k].col[j], sizeof (double));
		}
	}

	input.close();
}


void DataGrid2D::loadLevels(ifstream & input) {

	input.read((char*)&mode, sizeof (int));

	int l;
	input.read((char*)&l, sizeof (int));
	(this->levels) = l;

	this->levelStructure = new int[this->levels*2];

	for (int k = 0; k < this->levels*2; k+=2) {
		input.read((char *)&(this->levelStructure[k]),sizeof(int));
		input.read((char *)&(this->levelStructure[k+1]),sizeof(int));
	}
}
*/

void DataGrid2D::stripZeros(double limit) {

	int count = 0;
	for (int k = 0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			double val = this->getData(k,j,0);
			if (fabs(val) < limit) {
				count++;
				val = 0;	
			}
			this->setData(k,j,0,val);
		}		
	}		
}
