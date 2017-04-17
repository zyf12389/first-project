#pragma once
#include<iostream>
#include<cv.h>
#include<highgui.h>
using namespace std;
using namespace cv;
#define RE_COMPUTE_COST
void costAggregationST(const Mat&limg, const Mat& rimg, const int maxDis, Mat* costVol);
void costAggregationRightST(const Mat&limg, const Mat&rimg, const int maxDis, Mat* rcostVol);