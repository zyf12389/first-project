#ifndef FGS_UPSAMPLE_H
#define FGS_UPSAMPLE_H
#include<iostream>
#include<opencv2/opencv.hpp>
#include"..\pp.h"

using namespace cv;

void fgs_upsample(const Mat& rgb,const Mat& rgb_r,const Mat&disparity,const Mat& rdisparity,const int scale_rate,const int maxDis,int multiple,Mat& result,Mat& result_r);

#endif