#include "DataGrid3D.h"
//#include "CompressionUtils.h"
#include "Common.h"
#include <fstream>
#include <iostream>
using namespace std;
#include <math.h>


DataGrid3D::DataGrid3D(int xval, int yval, int zval)  { 
	x = xval;
	y = yval;
	z = zval;
	
	origx = xval;
	origy = yval;
	origz = zval;
	dimension = 3;
	
	mode = 0;
	
	data = new third[x];
	
	for (int k = 0; k < x; k++) {
		data[k].row = new rows[y];
		for (int j = 0; j < y; j++) {
			data[k].row[j].col = new double[z];
			bzero(data[k].row[j].col,sizeof(double)*z);	
		}
	}	
}

DataGrid3D::~DataGrid3D() {
	for (int k =0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			delete data[k].row[j].col;
		}
		delete data[k].row;
	}
	delete data;
}

void DataGrid3D::setData(int xval, int yval, int zval, double value) {
	if (data == NULL || xval > x || yval > y || zval > z) {
		return;
	}
	
	data[xval].row[yval].col[zval] = value;
}

double DataGrid3D::getData(int xval, int yval, int zval) {
	if (xval >= x || yval >= y || zval >= z) {
		cout << "REQUESTING BAD DATA: " << xval << "," << yval <<
		" and max size is: " << x << "," <<yval << endl;
	}
	
	return data[xval].row[yval].col[zval];
}

void DataGrid3D::output() {
	for (int k =0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			cout << k << "," << j << endl;
			for (int i = 0; i < z; i++) {
				cout <<"\t" << i << " = " << this->getData(k,j,i) << " ";
			}	
			cout << endl;
		}	
	}
	
	cout << " sized: " << x << " " << y << " " << z << endl;
}

void DataGrid3D::resize(int newx, int newy, int newz, bool copy) {
	
//	if (verbose) cout << " have to resize final array to: " << newx <<
//		"," << newy << "," << newz << endl;
		
	third * newdata = new third[newx];
	
	for (int k = 0; k < newx; k++) {
		newdata[k].row = new rows[newy];
		for (int j = 0; j < newy; j++) {
			newdata[k].row[j].col = new double[newz];
			bzero(newdata[k].row[j].col,sizeof(double)*newz);	
		}
	}	
	
	
	if (copy) {
		for (int k =0; k < x; k++) {
			for (int j = 0; j < y; j++) {
				for (int i = 0; i < z; i++) {
					newdata[k].row[j].col[i] = this->getData(k,j,i);
				}
			}	
		}
	}	

	delete data;
	data = newdata;
	
	x = newx;	
	y = newy;
	z = newz;
}

void DataGrid3D::copyTo(DataGrid * newData) {

	for (int k =0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			for (int i = 0; i < z; i++) {
				double val = this->getData(k,j,i);
				newData->setData(k,j,i,val);
			}
		}	
	}
}

void DataGrid3D::stripZeros(double limit) {
	
	int count = 0;
	for (int k = 0; k < x; k++) {
		for (int j = 0; j < y; j++) {
			for (int i = 0; i < z; i++) {
				double val = this->getData(k,j,i);
				if (fabs(val) < limit) {
					count++;
					val = 0;	
				}	
				this->setData(k,j,i,val);
			}
		}		
	}		
//	if (verbose)  cout << "Thresholding removed: " << count << " zeros with limit: " << limit << endl;
}
