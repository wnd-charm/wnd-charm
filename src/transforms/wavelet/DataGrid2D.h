#ifndef DATAGRID2D_H_
#define DATAGRID2D_H_

//#include <string>
//#include <fstream>
//using namespace std;

#include "DataGrid.h"

struct row {
	double * col;	
};

class DataGrid2D : public DataGrid
{
public:
//	DataGrid2D(const string & filename, bool compress);
	DataGrid2D(int xval, int yval, int zval);
	~DataGrid2D();
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

	row * data;

};

#endif /*DATAGRID2D_H_*/
