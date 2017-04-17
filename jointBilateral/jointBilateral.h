#ifndef JOINTBILATERAL_H_
#define  JOINTBILATERAL_H_

#include<iostream>
#include<opencv2\opencv.hpp>
using namespace std;
using namespace cv;

void jointBilateral(const Mat&color, Mat&disparity,const int radius, const int multiple,const int scale_rate, double sigma_c,  double sigma_s, Mat&dst);

#endif