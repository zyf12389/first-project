#pragma once
#include<iostream>
#include<cv.h>
#include<highgui.h>
using namespace std;
using namespace cv;

void costAggregation(const Mat&limg, const Mat& rimg,const int maxDis, Mat* costVol);
void costAggregationRight(const Mat&limg, const Mat&rimg,const int maxDis, Mat* rcostVol);
void cvtMatQX(const Mat& img, unsigned char* qxImg);