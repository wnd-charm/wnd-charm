#include "Wavelet.h"
#include "DataGrid3D.h"
#include "DataGrid2D.h"
//#include "DataGrid1D.h"
//#include "DataGrid.h"
#include "wt.h"
#include "Common.h"

//extern bool verbose;
//extern bool parallel;
//extern bool master;
//extern char * slave;


Wavelet::Wavelet()
{
//	remote = new RemoteProcessor();
}

Wavelet::~Wavelet()
{
}

void Wavelet::inverseTransform(DataGrid * data) {

	switch (data->getDimension()) {
		case 1:
//			this->inverseTransform1D(data);
			break;
		case 2:
//			this->inverseTransform2D(data);
			break;
		case 3:
//			this->inverseTransform3D(data);
			break;
	}
}

void Wavelet::transform(DataGrid * data) {

	switch (data->getDimension()) {
		case 1:
			this->transform1D(data);
			break;
		case 2:
			this->transform2D(data);	
			break;
		case 3:
//			this->transform3D(data);	
			break;	
	}
	
	data->stripZeros(limit);
}

void Wavelet::transform1D(DataGrid * data) {

	int size = data->getX();

	double ** arrays = new double*[this->nsteps + 1];
	int * steps = new int[this->nsteps];

	// get the input data into a double array 
	// this is wasted copying because of our dataGrid class
	double * input = new double[size];
	bzero(input,size);
	for (int k = 0; k  < size; k++) {
		input[k] = data->getData(k,0,0); 	
	}

	// for each decomposition level
	for (int step = 0; step < this->nsteps; ++step) {
			
		//int outsize = (size >> 1) + 1;
		//if (size % 2 == 1) { outsize++; }
		int outsize = (int) floor((size + this->dec_len - 1) / 2.);

		//cout << "out size x " << outsize << " size: " << size << endl;
		
		// create approx and detail vectors
		double * approx = new double[outsize];
		double * detail = new double[outsize];
		
		bzero(approx,outsize);
		bzero(detail,outsize);

		// get the approximation data
//		if (parallel && master) {
//			remote->remoteProcess(input, size, approx, detail, outsize, MODE_ZEROPAD, this->mode);
//		}
//		else
                {
  		    d_dec_a(input, size, this, approx, outsize, MODE_ZEROPAD);
		    d_dec_d(input, size, this, detail, outsize, MODE_ZEROPAD);
		}

		size = outsize;
		input = approx;

		// save the detail and then approx. next level the detail
		// will overwrite this approx, which is what we want to happen

		arrays[step] = detail;
		arrays[step+1] = approx;
		steps[step] = outsize;
//		if (verbose) cout << "Level " << step << " is size: " << outsize << endl << endl;
	}
	
	int length = 0;
	for (int k = 0; k < this->nsteps; k++) {
		length += steps[k];
//		if (verbose) cout << "\t" << steps[k] << endl;
	}
	length += steps[this->nsteps-1];

	// save teh structure of the decomposition
	data->levels = this->nsteps;
	data->setStructure(steps);
	
	data->resize(length,0,0,false);
	
	// put the data back into the data grid
	int index = 0;
	for (int k = 0; k < this->nsteps+1; k++) {
		double * array = arrays[k];
		int array_size = steps[k];
		if (k == this->nsteps) { array_size = steps[k-1]; }

		for (int j = 0; j < array_size; j++) {
			data->setData(index,0,0,array[j]); 
			index++;
		}
	}
}

void Wavelet::transform2D(DataGrid * data) {

	int xsize = data->getX();
	int ysize = data->getY();

	int arrayindex = 0, stepindex = 0;
	DataGrid ** arrays = new DataGrid*[3*this->nsteps+1];
	int * steps = new int[this->nsteps*2];
	DataGrid * inputGrid = data;

	int lengthx = 0, lengthy = 0;

	for (int step = 0; step < this->nsteps; ++step)
        {
		int outsizex = (int) floor((xsize + this->dec_len - 1) / 2.);
		int outsizey = (int) floor((ysize + this->dec_len - 1) / 2.);

		DataGrid * iGrid1 = new DataGrid2D(xsize,outsizey,0);
		DataGrid * iGrid2 = new DataGrid2D(xsize,outsizey,0);

		DataGrid * tempGrid1 = new DataGrid2D(outsizex,outsizey,0);
		DataGrid * tempGrid2 = new DataGrid2D(outsizex,outsizey,0);
		DataGrid * tempGrid3 = new DataGrid2D(outsizex,outsizey,0);
		DataGrid * tempGrid4 = new DataGrid2D(outsizex,outsizey,0);

		double * approx = new double[outsizey];
		double * detail = new double[outsizey];

		bzero(approx,outsizey);
		bzero(detail,outsizey);

   	        double * input = new double[ysize];
		// go through each row and get the approx and detail
		for (int j = 0; j < xsize; j++) {
			bzero(input,ysize);

			// copy whole row
			for (int k = 0; k  < ysize; k++) {
				input[k] = inputGrid->getData(j,k,0);
			}

//			if (parallel && master) {
//				remote->remoteProcess(input, ysize, approx, detail, outsizey, MODE_ZEROPAD, this->mode);
//			}
//			else
                        {
				d_dec_a(input, ysize, this, approx, outsizey, MODE_ZEROPAD);
				d_dec_d(input, ysize, this, detail, outsizey, MODE_ZEROPAD);
			}

			for (int index = 0; index < outsizey; index++) {
				iGrid1->setData(j,index,0,approx[index]);
				iGrid2->setData(j,index,0,detail[index]);
			}

			bzero(approx,outsizey);
			bzero(detail,outsizey);

		}
                delete [] input;


/* added here **********************************************************************/

// horizontal cooeficients
/*
delete approx;
delete detail;

		approx = new double[outsizex];
		detail = new double[outsizex];

		bzero(approx,outsizey);
		bzero(detail,outsizey);

   	        input = new double[xsize];
		// go through each row and get the approx and detail
		for (int j = 0; j < ysize; j++) {
			bzero(input,xsize);

			// copy whole row
			for (int k = 0; k  < xsize; k++) {
				input[k] = inputGrid->getData(k,j,0);
			}

//			if (parallel && master) {
//				remote->remoteProcess(input, ysize, approx, detail, outsizey, MODE_ZEROPAD, this->mode);
//			}
//			else
                        {
				d_dec_a(input, xsize, this, approx, outsizex, MODE_ZEROPAD);
				d_dec_d(input, xsize, this, detail, outsizex, MODE_ZEROPAD);
			}

			for (int index = 0; index < outsizex; index++) {
				tempGrid1->setData(index,j,0,approx[index]);
				tempGrid2->setData(index,j,0,detail[index]);
			}

			bzero(approx,outsizex);
			bzero(detail,outsizex);

		}
                delete input;


data->resize(outsizex,outsizey,0,false);
for (int x=0;x<outsizex;x++)
  for (int y=0;y<outsizey;y++)
    data->setData(x,y,0,iGrid2->getData(x,y,0)+tempGrid2->getData(x,y,0));
delete iGrid1;
delete iGrid2;
delete tempGrid1;
delete tempGrid2;
delete tempGrid3;
delete tempGrid4;
delete approx;
delete detail;
return;

added here (end) *********************************************************************
*/



		double * approx1 = new double[outsizex];
		double * detail1 = new double[outsizex];

		bzero(approx1,outsizex);
		bzero(detail1,outsizex);

		double * approx2 = new double[outsizex];
		double * detail2 = new double[outsizex];

		bzero(approx2,outsizex);
		bzero(detail2,outsizex);

		// now go through the columns of the approx and the detail
		// and get their approx and detail (4 vectors now) and
		// store them
  	        double * input1 = new double[xsize];
		double * input2 = new double[xsize];
		for (int j = 0; j < outsizey; j++) {

			bzero(input1,xsize);
			bzero(input2,xsize);

			// copy whole col
			for (int k = 0; k  < xsize; k++) {
				input1[k] = iGrid1->getData(k,j,0);
				input2[k] = iGrid2->getData(k,j,0);
			}

//			if (parallel && master) {
//				remote->remoteProcess(input1, xsize, approx1, detail1, outsizex, MODE_ZEROPAD, this->mode);
//			}
//			else
                        {
				d_dec_a(input1, xsize, this, approx1, outsizex, MODE_ZEROPAD);
				d_dec_d(input1, xsize, this, detail1, outsizex, MODE_ZEROPAD);
			}

//			if (parallel && master) {
//				remote->remoteProcess(input2, xsize, approx2, detail2, outsizex, MODE_ZEROPAD, this->mode);
//			}
//			else
                        {
				d_dec_a(input2, xsize, this, approx2, outsizex, MODE_ZEROPAD);
				d_dec_d(input2, xsize, this, detail2, outsizex, MODE_ZEROPAD);
			}

			for (int index = 0; index < outsizex; index++) {
				tempGrid1->setData(index,j,0,approx1[index]);
				tempGrid2->setData(index,j,0,detail1[index]);
				tempGrid3->setData(index,j,0,approx2[index]);
				tempGrid4->setData(index,j,0,detail2[index]);
			}

		}
                delete [] input1;
                delete [] input2;

		xsize = outsizex;
		ysize = outsizey;
		inputGrid = tempGrid1;

		// save the details, and the approx. the approx
		// gets overwritten in the next level of deconstruction
		arrays[arrayindex] = tempGrid4;
		arrays[arrayindex+1] = tempGrid3;
		arrays[arrayindex+2] = tempGrid2;
		arrays[arrayindex+3] = tempGrid1;
		arrayindex += 3;

		lengthx += outsizex;
		lengthy += outsizey;

//		if (verbose) cout << "Level " << step << ":  X size: " << outsizey <<
//			"and Y size: " << outsizey << endl;

		// store the structure of the decomposition
		steps[stepindex++] = outsizex;
		steps[stepindex++] = outsizey;

                delete iGrid1;
                delete iGrid2;
                delete [] approx;
                delete [] detail;
                delete [] approx1;
                delete [] detail1;
                delete [] approx2;
                delete [] detail2;

	}

	lengthx += xsize;
	lengthy += ysize;

	data->levels = this->nsteps;
	data->setStructure(steps);

	data->resize(lengthx,lengthy,0,false);

	int indexx = 0;
	int indexy = 0;
	int oldx = indexx;
	int oldy = indexy;
	int count = 0;

	// we need some fancy footwork here to store
	// the data properly in a recursive manner and
	// in a way that can be reformed later.
	for (int k = 0; k < this->nsteps; k++) {

		oldx = indexx;
		oldy = indexy;

		DataGrid * tmp1 = arrays[k*3];
		DataGrid * tmp2 = arrays[k*3+1];
		DataGrid * tmp3 = arrays[k*3+2];
		//DataGrid * tmp4 = arrays[k*3+3];

		int lastx = 0, lasty = 0;

		for (int x = 0; x < tmp1->getX(); x++) {
			for (int y = 0; y < tmp1->getY(); y++, indexy++) {
				data->setData(indexx,indexy,0,tmp1->getData(x,y,0));
				count++;
			}
			lasty = indexy;
			indexy = oldy;
			indexx++;
		}
		lastx = indexx;

		indexy = lasty;
		indexx = oldx;

		for (int x = 0; x < tmp2->getX(); x++) {
			for (int y = 0; y < tmp2->getY(); y++,indexy++) {
				data->setData(indexx,indexy,0,tmp2->getData(x,y,0));
				count++;
			}
			indexy = lasty;
			indexx++;
		}

		indexx = lastx;
		indexy = oldy;
		for (int x = 0; x < tmp3->getX(); x++) {
			for (int y = 0; y < tmp3->getY(); y++, indexy++) {
				data->setData(indexx,indexy,0,tmp3->getData(x,y,0));
				count++;
			}
			indexy = oldy;
			indexx++;
		}
		indexx = lastx;
		indexy = lasty;
	}

	oldy = indexy;

	// get last approx vector
	DataGrid * tmp1 = arrays[this->nsteps*3];
	for (int x = 0; x < tmp1->getX(); x++) {
		for (int y = 0; y < tmp1->getY(); y++,indexy++) {
			data->setData(indexx,indexy,0,tmp1->getData(x,y,0));
			count++;
		}
		indexx++;
		indexy = oldy;
	}

        for (int a=0;a<3*this->nsteps+1;a++)
          delete arrays[a];
        delete [] arrays;
        delete [] steps;
}


void Wavelet::transform3D(DataGrid * data) {
	int xsize = data->getX();
	int ysize = data->getY();
	int zsize = data->getZ();

	int arrayindex = 0, stepindex = 0;
	DataGrid ** arrays = new DataGrid*[7*this->nsteps+1];
	int * steps = new int[this->nsteps*3];
	DataGrid * inputGrid = data;

		int lengthx = 0, lengthy = 0, lengthz = 0;

	for (int step = 0; step < this->nsteps; ++step) {

		int outsizex = (int) floor((xsize + this->dec_len - 1) / 2.);
		int outsizey = (int) floor((ysize + this->dec_len - 1) / 2.);
		int outsizez = (int) floor((zsize + this->dec_len - 1) / 2.);

//		if (verbose) cout << " step " << step << "x: " << outsizex << " y: " << outsizey
//			<< " z: " << outsizez << endl;

		DataGrid * l1Grid1 = new DataGrid3D(xsize,ysize,outsizez);
		DataGrid * l1Grid2 = new DataGrid3D(xsize,ysize,outsizez);
		
		DataGrid * l2Grid1 = new DataGrid3D(xsize,outsizey,outsizez);
		DataGrid * l2Grid2 = new DataGrid3D(xsize,outsizey,outsizez);
		DataGrid * l2Grid3 = new DataGrid3D(xsize,outsizey,outsizez);
		DataGrid * l2Grid4 = new DataGrid3D(xsize,outsizey,outsizez);

		DataGrid * l3Grid1 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid2 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid3 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid4 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid5 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid6 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid7 = new DataGrid3D(outsizex,outsizey,outsizez);
		DataGrid * l3Grid8 = new DataGrid3D(outsizex,outsizey,outsizez);
	
		
		double * approx = new double[outsizez];
		double * detail = new double[outsizez];
		
		bzero(approx,outsizez);
		bzero(detail,outsizez);
		
		// go through each row and get the approx and detail
		for (int j = 0; j < xsize; j++) {
			for (int k = 0; k < ysize; k++) {
				double * input = new double[zsize];
				bzero(input,zsize);
				
				// copy whole row
				for (int i = 0; i  < zsize; i++) {
					input[i] = inputGrid->getData(j,k,i); 	
				}
				
//				if (parallel && master) {
//					remote->remoteProcess(input, zsize, approx, detail, outsizez, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input, zsize, this, approx, outsizez, MODE_ZEROPAD);
					d_dec_d(input, zsize, this, detail, outsizez, MODE_ZEROPAD);				
//				}
				
				for (int index = 0; index < outsizez; index++) {
					l1Grid1->setData(j,k,index,approx[index]);	
					l1Grid2->setData(j,k,index,detail[index]);
				}
				
				bzero(approx,outsizez);
				bzero(detail,outsizez);
			}
		}

		
		delete approx;
		delete detail;

		double * approx1 = new double[outsizey];
		double * detail1 = new double[outsizey];
		
		bzero(approx1,outsizey);
		bzero(detail1,outsizey);

		double * approx2 = new double[outsizey];
		double * detail2 = new double[outsizey];
		
		bzero(approx2,outsizey);
		bzero(detail2,outsizey);

		// now go through the columns of the approx and the detail
		// and get their approx and detail (4 vectors now) and
		// store them
		for (int k = 0; k < xsize; k++) {
			
			double * input1 = new double[ysize];
			bzero(input1,ysize);
			
			double * input2 = new double[ysize];
			bzero(input2,ysize);
			
			// copy whole col
			for (int j = 0; j  < outsizez; j++) {
				for (int i = 0; i < ysize; i++) {
					input1[i] = l1Grid1->getData(k,i,j);
					input2[i] = l1Grid2->getData(k,i,j); 
				}	
			
//				if (parallel && master) {
//					remote->remoteProcess(input1, ysize, approx1, detail1, outsizey, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input1, ysize, this, approx1, outsizey, MODE_ZEROPAD);
					d_dec_d(input1, ysize, this, detail1, outsizey, MODE_ZEROPAD);
//				}

//				if (parallel && master) {
//					remote->remoteProcess(input2, ysize, approx2, detail2, outsizey, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input2, ysize, this, approx2, outsizey, MODE_ZEROPAD);
					d_dec_d(input2, ysize, this, detail2, outsizey, MODE_ZEROPAD);
//				}

				for (int index = 0; index < outsizey; index++) {
					l2Grid1->setData(k,index,j,approx1[index]);	
					l2Grid2->setData(k,index,j,detail1[index]);				
					l2Grid3->setData(k,index,j,approx2[index]);	
					l2Grid4->setData(k,index,j,detail2[index]);
				}
			}
		}
		
		delete approx1;
		delete approx2;
		delete detail1;
		delete detail2;
			
		// NOW GO THROUGH AND DECIMATE THE FINAL
		// ROW - THE X ROW	
			   
		approx1 = new double[outsizex];
		detail1 = new double[outsizex];
		approx2 = new double[outsizex];
		detail2 = new double[outsizex];
		double * approx3 = new double[outsizex];
		double * detail3 = new double[outsizex];		
		double * approx4 = new double[outsizex];
		double * detail4 = new double[outsizex];
		
		bzero(approx1,outsizex);
		bzero(detail1,outsizex);
		bzero(approx2,outsizex);
		bzero(detail2,outsizex);
		bzero(approx3,outsizex);
		bzero(detail3,outsizex);
		bzero(approx4,outsizex);
		bzero(detail4,outsizex);
		
		// now go through the columns of the approx and the detail
		// and get their approx and detail (4 vectors now) and
		// store them
		for (int k = 0; k < outsizey; k++) {
			
			double * input1 = new double[xsize];
			double * input2 = new double[xsize];
			double * input3 = new double[xsize];
			double * input4 = new double[xsize];
			bzero(input1,xsize);
			bzero(input2,xsize);
			bzero(input3,xsize);
			bzero(input4,xsize);
			
			// copy whole col
			for (int j = 0; j < outsizez; j++) {
				for (int i = 0; i < xsize; i++) {
					input1[i] = l2Grid1->getData(i,k,j); 	
					input2[i] = l2Grid2->getData(i,k,j); 
					input3[i] = l2Grid3->getData(i,k,j); 
					input4[i] = l2Grid4->getData(i,k,j); 
				}
			
			
//				if (parallel && master) {
//					remote->remoteProcess(input1, xsize, approx1, detail1, outsizex, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input1, xsize, this, approx1, outsizex, MODE_ZEROPAD);
					d_dec_d(input1, xsize, this, detail1, outsizex, MODE_ZEROPAD);
//				}

//				if (parallel && master) {
//					remote->remoteProcess(input2, xsize, approx2, detail2, outsizex, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input2, xsize, this, approx2, outsizex, MODE_ZEROPAD);
					d_dec_d(input2, xsize, this, detail2, outsizex, MODE_ZEROPAD);
//				}

//				if (parallel && master) {
//					remote->remoteProcess(input3, xsize, approx3, detail3, outsizex, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input3, xsize, this, approx3, outsizex, MODE_ZEROPAD);
					d_dec_d(input3, xsize, this, detail3, outsizex, MODE_ZEROPAD);
//				}

//				if (parallel && master) {
//					remote->remoteProcess(input4, xsize, approx4, detail4, outsizex, MODE_ZEROPAD, this->mode);
//				}
//				else {
					d_dec_a(input4, xsize, this, approx4, outsizex, MODE_ZEROPAD);
					d_dec_d(input4, xsize, this, detail4, outsizex, MODE_ZEROPAD);
//				}

				for (int index = 0; index < outsizex; index++) {
					l3Grid1->setData(index,k,j,approx1[index]);
					l3Grid2->setData(index,k,j,detail1[index]);
					l3Grid3->setData(index,k,j,approx2[index]);
					l3Grid4->setData(index,k,j,detail2[index]);
					l3Grid5->setData(index,k,j,approx3[index]);
					l3Grid6->setData(index,k,j,detail3[index]);
					l3Grid7->setData(index,k,j,approx4[index]);
					l3Grid8->setData(index,k,j,detail4[index]);
				}
			}
		}

		delete approx1;
		delete approx2;
		delete approx3;
		delete approx4;
		delete detail1;
		delete detail2;
		delete detail3;
		delete detail4;
		
		xsize = outsizex;
		ysize = outsizey;
		zsize = outsizez;
		inputGrid = l3Grid1;
		
		// save the details, and the approx. the approx
		// gets overwritten in the next level of deconstruction
		arrays[arrayindex] = l3Grid8;
		arrays[arrayindex+1] = l3Grid7;
		arrays[arrayindex+2] = l3Grid6;
		arrays[arrayindex+3] = l3Grid5;
		arrays[arrayindex+4] = l3Grid4;
		arrays[arrayindex+5] = l3Grid3;
		arrays[arrayindex+6] = l3Grid2;
		arrays[arrayindex+7] = l3Grid1;
		arrayindex += 7;
				
		lengthx += outsizex;
		lengthy += outsizey;
		lengthz += outsizez;

//		if (verbose) cout << "Level " << step << ":  X size: " << outsizey <<
//			"and Y size: " << outsizey << " and Z size: " << outsizez << endl;
		
		// store the structure of the decomposition
		steps[stepindex++] = outsizex;
		steps[stepindex++] = outsizey;
		steps[stepindex++] = outsizez;
	}
	
	lengthx += xsize;
	lengthy += ysize;
	lengthz += zsize;
	
	data->levels = this->nsteps;
	data->setStructure(steps);
	
	data->resize(lengthx,lengthy,lengthz,false);
	
	int indexx = 0;
	int indexy = 0;
	int indexz = 0;
	int oldx = indexx;
	int oldy = indexy;
	int oldz = indexz;

	// we need some fancy footwork here to store
	// the data properly in a recursive manner and
	// in a way that can be reformed later.
	for (int k = 0; k < this->nsteps; k++) {
		oldx = indexx;
		oldy = indexy;
		oldz = indexz;
		
		//cout << indexx << " " << indexy << " " << indexz << endl;
				
		int lastx = 0, lasty = 0, lastz = 0;

		for (int g = 0; g < 7; g++) {
			DataGrid * tmp = arrays[k*7+g];

			for (int x = 0; x < tmp->getX(); x++, indexx++) {
				for (int y = 0; y < tmp->getY(); y++, indexy++) {
					for (int z = 0; z < tmp->getZ(); z++, indexz++) {
						data->setData(indexx,indexy,indexz,tmp->getData(x,y,z)); 
//							if (g == 0) { cout << "DDD " << 
//								indexx << " " << indexy << " " << indexz <<  " " << tmp->getData(x,y,z) << endl; }
					}
					if (g == 0) { lastz = indexz; indexz = oldz; }
					else if (g == 1) { indexz = oldz; }
					else if (g == 2) { indexz = lastz; }
					else if (g == 3) { indexz = lastz; } 
					else if (g == 4) { indexz = oldz; } 
					else if (g == 5) { indexz = oldz; } 
					else if (g == 6) { indexz = lastz; } 
				}
				if (g < 4) { lasty = indexy; indexy = oldy; }
				else { indexy = lasty; }
			}

			if (g == 0) { lastx = indexx;  indexz = oldz; }
			else if (g == 1) { indexx = oldx;  indexz = lastz; }
			else if (g == 2) { indexx = lastx; indexz = lastz; }
			else if (g == 3) { indexx = oldx;  indexz = oldz; indexy = lasty; }
			else if (g == 4) { indexx = lastx; indexz = oldz;  }
			else if (g == 5) { indexx = oldx; indexz = lastz; }
			else if (g == 6) { indexx = lastx; indexz = lastz; }			
		}

		indexx = lastx;
		indexy = lasty;
		indexz = lastz;		
	}
	
	oldz = indexz;
	oldy = indexy;
	oldx = indexx;
	
	// get last approx vector
	DataGrid * tmp1 = arrays[this->nsteps*7];
	for (int x = 0; x < tmp1->getX(); x++, indexx++) {
		for (int y = 0; y < tmp1->getY(); y++,indexy++) {
			for (int z = 0; z < tmp1->getZ(); z++, indexz++) {
				data->setData(indexx,indexy,indexz,tmp1->getData(x,y,z)); 
			}
			indexz = oldz;
		}
		indexy = oldy;
	}
}

/*

void Wavelet::inverseTransform1D(DataGrid * data) {


	int * steps = new int[this->nsteps];

	int size = data->getX();

	double * input = new double[size];
	bzero(input,size);
	for (int k = 0; k  < size; k++) {
		input[k] = data->getData(k,0,0);
	}

	int stepsize = data->levelStructure[this->nsteps - 1]; 
	double * a_data = new double[stepsize];
	bzero(a_data,stepsize);
	
	for (int index = 0, k = size-(stepsize); k < size; k++, index++) {
		a_data[index] = input[k]; 
	}
	
	
	for (int step = this->nsteps-1; step >= 0; --step) {
		
		stepsize = data->levelStructure[step]; 
		
		int start = 0;
		for (int j = 0; j < step; j++) {
			start += data->levelStructure[j];
		}

	 	double * d_data = new double[stepsize];
		bzero(d_data,stepsize);
		
//		if (verbose) cout << "Reconstruction level: " << step << " step size : " << stepsize << "start index: " << start << endl;	
		
		for (int index = 0, k = start; index < stepsize; k++, index++) {
			d_data[index] = input[k];  
		}
	
		int reconstruction_len = idwt_buffer_length(stepsize, this->rec_len, MODE_ZEROPAD);

	    double * output = new double[reconstruction_len];
	  
	    d_idwt(a_data, stepsize, d_data, stepsize,this, output, 
	    		reconstruction_len, MODE_ZEROPAD, 1);
	    		
	    a_data = output;
	    size = reconstruction_len;
	}
		
	for (int k = 0; k < size; k++) {
		data->setData(k,0,0,a_data[k]); 	
	}
}

void Wavelet::inverseTransform2D(DataGrid * data) {

	int xsize = data->getX();
	int ysize = data->getY();

	int stepsizey = data->levelStructure[2*this->nsteps - 1];
	int stepsizex = data->levelStructure[2*this->nsteps - 2];
	
	DataGrid * input = data;
	
	DataGrid * a_grid = new DataGrid2D(stepsizex,stepsizey,0);

	if (verbose) cout << "data size : " << xsize << " step x: " << stepsizex << " step y: " << stepsizey << endl;

	int startx = 0, starty = 0;
	for (int j = 0; j < (this->nsteps-1)*2; j+=2) {
		startx += data->levelStructure[j];
		starty += data->levelStructure[j+1];
	}
	
	int indexx, indexy;
	for (int k = startx+stepsizex, indexx = 0; k < startx+2*stepsizex; k++, indexx++) {
		for (int j = starty+stepsizey, indexy = 0; j < starty+2*stepsizey; j++, indexy++) {
			a_grid->setData(indexx,indexy,0,data->getData(k,j,0)); 
		}
	}
	
	for (int step = this->nsteps-1; step >= 0; --step) {
		
 		stepsizex = data->levelStructure[2*step]; 
		stepsizey = data->levelStructure[2*step+1]; 
		
		startx = 0, starty = 0;
		for (int j = 0; j < step*2; j+=2) {
			startx += data->levelStructure[j];
			starty += data->levelStructure[j+1];
		}

		if (verbose) cout << "step: " << step << " step size x: " << stepsizex << " step y: " << stepsizey << endl;

		DataGrid * d_grid1 = new DataGrid2D(stepsizex,stepsizey,0);
	 	for (int k = startx + stepsizex, indexx = 0; k < startx + 2*stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				d_grid1->setData(indexx,indexy,0,input->getData(k,j,0)); 
			}
		}

		DataGrid * d_grid2 = new DataGrid2D(stepsizex,stepsizey,0);
	 	for (int k = startx, indexx = 0; k < startx + stepsizex; k++, indexx++) {
			for (int j = starty + stepsizey, indexy = 0; j < starty + 2*stepsizey; j++, indexy++) {
				d_grid2->setData(indexx,indexy,0,input->getData(k,j,0)); 
			}
		}

		DataGrid * d_grid3 = new DataGrid2D(stepsizex,stepsizey,0);
	 	for (int k = startx, indexx = 0; k < startx+stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				d_grid3->setData(indexx,indexy,0,input->getData(k,j,0)); 
			}
		}

		int reconstruction_lenx = idwt_buffer_length(stepsizex, this->rec_len, MODE_ZEROPAD);
		int reconstruction_leny = idwt_buffer_length(stepsizey, this->rec_len, MODE_ZEROPAD);

		// take cols data, recon
		// take rows data recon,
		// recon both
	   
	   	DataGrid * o_grid1 = new DataGrid2D(reconstruction_lenx,stepsizey,0);
	   	DataGrid * o_grid2 = new DataGrid2D(reconstruction_lenx,stepsizey,0);

   		for (int k = 0; k < stepsizey; k++) {
		    double output1[reconstruction_lenx];
		    double output2[reconstruction_lenx];
		    
		    double a_data[stepsizex];
		    double d_data1[stepsizex];
		    double d_data2[stepsizex];
		    double d_data3[stepsizex];
			for (int i = 0; i < stepsizex; i++) {
				a_data[i] = a_grid->getData(i,k,0);

				d_data1[i] = d_grid1->getData(i,k,0);
				d_data2[i] = d_grid2->getData(i,k,0);
				d_data3[i] = d_grid3->getData(i,k,0);
			}

		    d_idwt(a_data, stepsizex, d_data1, stepsizex,this, output1, 
	    		reconstruction_lenx, MODE_ZEROPAD, 1);

		    d_idwt(d_data2, stepsizex, d_data3, stepsizex,this, output2, 
	    		reconstruction_lenx, MODE_ZEROPAD, 1);
	    		
	    	for (int i = 0; i < reconstruction_lenx; i++) {
				o_grid1->setData(i,k,0,output1[i]);
				o_grid2->setData(i,k,0,output2[i]);
			}
   		}
   		
   		DataGrid * o_grid = new DataGrid2D(reconstruction_lenx,reconstruction_leny,0);
   				
   		for (int k = 0; k < reconstruction_lenx; k++) {
		    double output[reconstruction_leny];
		   		    
		    double a_data[stepsizey];
		    double d_data[stepsizey];
			for (int i = 0; i < stepsizey; i++) {
				a_data[i] = o_grid1->getData(k,i,0);
				d_data[i] = o_grid2->getData(k,i,0);
			}
	
		    d_idwt(a_data, stepsizey, d_data, stepsizey,this, output, 
	    		reconstruction_leny, MODE_ZEROPAD, 1);

	    	for (int i = 0; i < reconstruction_leny; i++) {
				o_grid->setData(k,i,0,output[i]);
			}
   		}		
   		
	    delete a_grid;
	    delete o_grid1;
	    delete o_grid2;	
	    a_grid = o_grid;
	}

	a_grid->copyTo(data);	
}


void Wavelet::inverseTransform3D(DataGrid * data) {
	int xsize = data->getX();
	int ysize = data->getY();
	int zsize = data->getZ();

	int stepsizez = data->levelStructure[3*this->nsteps - 1];
	int stepsizey = data->levelStructure[3*this->nsteps - 2]; 
	int stepsizex = data->levelStructure[3*this->nsteps - 3]; 
	
	DataGrid * input = data;
	
	DataGrid * aaa = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	
	if (verbose) cout << "data size : " << xsize << " step x" << stepsizex 
		<< " step y: " << stepsizey << " step z: " << stepsizez << endl;

	int startx = 0, starty = 0, startz = 0;
	for (int j = 0; j < (this->nsteps-1)*3; j+=3) {
		startx += data->levelStructure[j];
		starty += data->levelStructure[j+1];
		startz += data->levelStructure[j+2];
	}
	
	//if (verbose) cout << "start: " << startx << " " << starty << " "  << startz << endl;			
	
	int indexx =0, indexy=0, indexz=0;
	for (int k = startx+stepsizex, indexx = 0; k < startx+2*stepsizex; k++, indexx++) {
		for (int j = starty+stepsizey, indexy = 0; j < starty+2*stepsizey; j++, indexy++) {
			for (int i = startz+stepsizez, indexz = 0; i < startz+2*stepsizez; i++, indexz++) {
				aaa->setData(indexx,indexy,indexz,data->getData(k,j,i)); 
			}
		}
	}
	
	for (int step = this->nsteps-1; step >= 0; --step) {

 		stepsizex = data->levelStructure[3*step];
		stepsizey = data->levelStructure[3*step+1];
		stepsizez = data->levelStructure[3*step+2];

		if (verbose) cout << "step: " << step << " step size x: " << stepsizex << " step y: " << stepsizey << " z: " << stepsizez << endl;					

		startx = 0, starty = 0, startz = 0;
		for (int j = 0; j < step*3; j+=3) {
			startx += data->levelStructure[j];
			starty += data->levelStructure[j+1];
			startz += data->levelStructure[j+2];
		}

		//if (verbose) cout << "start: " << startx << " " << starty << " "  << startz << endl;			

		DataGrid * aad = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx, indexx = 0; k < startx+stepsizex; k++, indexx++) {
			for (int j = starty+stepsizey, indexy = 0; j < starty+2*stepsizey; j++, indexy++) {
				for (int i = startz+stepsizez, indexz = 0; i < startz + 2*stepsizez; i++, indexz++) {
					aad->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}

						
		DataGrid * ada = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx+stepsizex, indexx = 0; k < startx+2*stepsizex; k++, indexx++) {
			for (int j = starty+stepsizey, indexy = 0; j < starty+2*stepsizey; j++, indexy++) {
				for (int i = startz, indexz = 0; i < startz + stepsizez; i++, indexz++) {
					ada->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}

		DataGrid * add = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx, indexx = 0; k < startx+stepsizex; k++, indexx++) {
			for (int j = starty+stepsizey, indexy = 0; j < starty+2*stepsizey; j++, indexy++) {
				for (int i = startz, indexz = 0; i < startz + stepsizez; i++, indexz++) {
					add->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}
		
		DataGrid * daa = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx+stepsizex, indexx = 0; k < startx+2*stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				for (int i = startz+stepsizez, indexz = 0; i < startz + 2*stepsizez; i++, indexz++) {
					daa->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}
		
		DataGrid * dad = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx, indexx = 0; k < startx+stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				for (int i = startz+stepsizez, indexz = 0; i < startz + 2*stepsizez; i++, indexz++) {
					dad->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}
		
		DataGrid * dda = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx+stepsizex, indexx = 0; k < startx+2*stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				for (int i = startz, indexz = 0; i < startz + stepsizez; i++, indexz++) {
					dda->setData(indexx,indexy,indexz,input->getData(k,j,i));
				} 
			}
		}
		
		DataGrid * ddd = new DataGrid3D(stepsizex,stepsizey,stepsizez);
	 	for (int k = startx, indexx = 0; k < startx+stepsizex; k++, indexx++) {
			for (int j = starty, indexy = 0; j < starty+stepsizey; j++, indexy++) {
				for (int i = startz, indexz = 0; i < startz + stepsizez; i++, indexz++) {
					ddd->setData(indexx,indexy,indexz,input->getData(k,j,i));
				}
			}
		}
	

		int reconstruction_lenx = idwt_buffer_length(stepsizex, this->rec_len, MODE_ZEROPAD);
		int reconstruction_leny = idwt_buffer_length(stepsizey, this->rec_len, MODE_ZEROPAD);
		int reconstruction_lenz = idwt_buffer_length(stepsizez, this->rec_len, MODE_ZEROPAD);

		// take cols data, recon
		// take rows data recon,
		// recon both
	   
	   	DataGrid * aa = new DataGrid3D(reconstruction_lenx,stepsizey,stepsizez);
	   	DataGrid * ad = new DataGrid3D(reconstruction_lenx,stepsizey,stepsizez);
	    DataGrid * da = new DataGrid3D(reconstruction_lenx,stepsizey,stepsizez);
	   	DataGrid * dd = new DataGrid3D(reconstruction_lenx,stepsizey,stepsizez);
	   
   		for (int j = 0; j < stepsizey; j++) {
   			for (int i = 0; i < stepsizez; i++) {

			    double output1[reconstruction_lenx];
			    double output2[reconstruction_lenx];
			    double output3[reconstruction_lenx];
			    double output4[reconstruction_lenx];
			    
			    double aaa_data[stepsizex];
			    double aad_data[stepsizex];
			    double ada_data[stepsizex];
			    double add_data[stepsizex];
			    double daa_data[stepsizex];
			    double dad_data[stepsizex];
			    double dda_data[stepsizex];
			    double ddd_data[stepsizex];
		    	
				for (int k = 0; k < stepsizex; k++) {
					aaa_data[k] = aaa->getData(k,j,i);
					aad_data[k] = aad->getData(k,j,i);
					ada_data[k] = ada->getData(k,j,i);
					add_data[k] = add->getData(k,j,i);
					daa_data[k] = daa->getData(k,j,i);
					dad_data[k] = dad->getData(k,j,i);
					dda_data[k] = dda->getData(k,j,i);
					ddd_data[k] = ddd->getData(k,j,i);
				}		
			

			    d_idwt(aaa_data, stepsizex, aad_data, stepsizex,this, output1, 
		    		reconstruction_lenx, MODE_ZEROPAD, 1);
			    d_idwt(ada_data, stepsizex, add_data, stepsizex,this, output2, 
		    		reconstruction_lenx, MODE_ZEROPAD, 1);
		   		d_idwt(daa_data, stepsizex, dad_data, stepsizex,this, output3, 
		    		reconstruction_lenx, MODE_ZEROPAD, 1);
				d_idwt(dda_data, stepsizex, ddd_data, stepsizex,this, output4,
		    		reconstruction_lenx, MODE_ZEROPAD, 1);
		    		
		    	for (int k = 0; k < reconstruction_lenx; k++) {
					aa->setData(k,j,i,output1[k]);
					ad->setData(k,j,i,output2[k]);
					da->setData(k,j,i,output3[k]);
					dd->setData(k,j,i,output4[k]);
				}
   			}
   		}

   		DataGrid * a = new DataGrid3D(reconstruction_lenx,reconstruction_leny,stepsizez);
   		DataGrid * d = new DataGrid3D(reconstruction_lenx,reconstruction_leny,stepsizez);
   				
   		for (int k = 0; k < reconstruction_lenx; k++) {
   			for (int j = 0; j < stepsizez; j++) {
			    double output1[reconstruction_leny];
			    double output2[reconstruction_leny];
			   		    
			    double aa_data[stepsizey];
			    double ad_data[stepsizey];
			    double da_data[stepsizey];
			    double dd_data[stepsizey];
			    
				for (int i = 0; i < stepsizey; i++) {
					aa_data[i] = aa->getData(k,i,j);
					ad_data[i] = ad->getData(k,i,j);
					da_data[i] = da->getData(k,i,j);
					dd_data[i] = dd->getData(k,i,j);
				}

			    d_idwt(aa_data, stepsizey, ad_data, stepsizey,this, output1,
		    		reconstruction_leny, MODE_ZEROPAD, 1);
		    	d_idwt(da_data, stepsizey, dd_data, stepsizey,this, output2, 
		    		reconstruction_leny, MODE_ZEROPAD, 1);
		    		
		    	for (int i = 0; i < reconstruction_leny; i++) {
					a->setData(k,i,j,output1[i]);
					d->setData(k,i,j,output2[i]);
				}
   			}
   		}	

   		DataGrid * complete = new DataGrid3D(reconstruction_lenx,reconstruction_leny,reconstruction_lenz);
   			
   		for (int k = 0; k < reconstruction_lenx; k++) {
   			for (int j = 0; j < reconstruction_leny; j++) {	
			    double output[reconstruction_lenz];
			   		    
			    double a_data[stepsizez];
			    double d_data[stepsizez];
				for (int i = 0; i < stepsizez; i++) {
					a_data[i] = a->getData(k,j,i);
					d_data[i] = d->getData(k,j,i);
				}
		
			    d_idwt(a_data, stepsizez, d_data, stepsizez,this, output, 
		    		reconstruction_lenz, MODE_ZEROPAD, 1);
		    		
		    	for (int i = 0; i < reconstruction_lenz; i++) {
					complete->setData(k,j,i,output[i]);
				}
	   		}
   		}


	    delete a;  delete d;
	    delete aa; delete ad; delete da; delete dd;
	    delete aaa;	delete aad; delete ada; delete add;
	    delete daa; delete dad; delete dda; delete ddd;
	    aaa = complete;
	}
	aaa->copyTo(data);

}
*/
