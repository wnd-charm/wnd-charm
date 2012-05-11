#ifndef DATAGRID_H_
#define DATAGRID_H_

//#include <string>
//using namespace std;

class DataGrid
{
public:
	DataGrid() {} ;
	virtual ~DataGrid() { };

	int getX() { return (int)x; }
	int getY() { return (int)y; }
	int getZ() { return (int)z; }

	int getOriginalX() { return (int)origx; }
	int getOriginalY() { return (int)origy; }
	int getOriginalZ() { return (int)origz; }

	int getMode() { return mode; }
	void setMode(int m) { mode = m; }

	virtual void setStructure(int * s) { this->levelStructure = s; }
	void setLevels(int l) { this->levels = l; }

//	static int determineDimension(const string &filename);

//	double * decompress(ifstream &input);

	int getDimension() { return dimension; }
	double getOriginalFileSize() { return originalFileSize; }


	virtual void output() = 0;

	virtual double getData(int x, int y, int z) = 0;
	virtual void setData(int xval, int yval, int zval, double value) = 0;

	virtual void copyTo(DataGrid * newData) = 0;

	virtual void stripZeros(double limit) = 0;

	virtual void resize(int newx, int newy, int newz, bool copy) = 0;

	int * levelStructure;
	int levels;

protected:
	virtual void loadData() { };

	int ex;

	int x;
	int y;
	int z;

	int origx;
	int origy;
	int origz;
	
	int periodsize;
	
	int mode;
	int dimension;
	
	double originalFileSize;
};

#endif /*DATAGRID_H_*/
