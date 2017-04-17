#ifndef COMPUTEDEPTH_H_
#define COMPUTEDEPTH_H_
#include<opencv2\opencv.hpp>
using namespace cv;

void computeDepth(Mat& ldisp,const double Tx, const double f,const double pixelSize, Mat&depth_W);

void computeDepth_T(Mat& rdisp, const double Tx, const double f, const double pixelSize, Mat &depth_T);

#endif
