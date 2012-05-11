#ifndef DataGrid3D_H_
#define DataGrid3D_H_

//#include <string>
//#include <fstream>
//using namespace std;

#include "DataGrid.h"

struct rows {
	double * col;
};

struct third {
	rows * row;
};

class DataGrid3D : public DataGrid
{
public:
//	DataGrid3D(const string & filename, bool compress);
	DataGrid3D(int xval, int yval, int zval);
	~DataGrid3D();
	double getData(int x, int y, int z);
	void setData(int xval, int yval, int zval, double value);

	void copyTo(DataGrid * newData);

	void output();

	void stripZeros(double limit);
	void resize(int newx, int newy, int newz, bool copy);

protected:
//	void loadData(const string & filename, bool compress);
//	void loadCompressedData(ifstream & input);

	void setNewSize(int xval, int yval);
//	void loadLevels(ifstream & input);

	third * data;

};

#endif /*DataGrid3D_H_*/
