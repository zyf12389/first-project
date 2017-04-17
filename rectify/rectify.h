#ifndef RECTIFY_H_
#define RECTIFY_H_
#include<iostream>
#include<fstream>
#include<sstream>
#include<opencv2\opencv.hpp>
using namespace std;
using namespace cv;

void rectify(const Mat& img_T, const Mat& img_W,const string filepath,Mat &roi_T,Mat &roi_W,double &Tx,double &f);
void mystereoRectify(InputArray _cameraMatrix1, InputArray _distCoeffs1,
	InputArray _cameraMatrix2, InputArray _distCoeffs2,
	Size imageSize, InputArray _Rmat, InputArray _Tmat,
	OutputArray _Rmat1, OutputArray _Rmat2,
	OutputArray _Pmat1, OutputArray _Pmat2,
	OutputArray _Qmat, int flags,
	double alpha, Size newImageSize,
	Rect* validPixROI1, Rect* validPixROI2);

#endif