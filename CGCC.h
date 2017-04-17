#pragma once
#include<iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

#define CENCUS_WND 13
#define CENCUS_BIT 168
#define LAMBDA_CEN 0.5
#define GRD_BOUND  0.1


#define CG_BORDER_THRES 1.0
#define CG_TAU_1 1.0
#define CG_TAU_2 0.0
#define CG_ALPHA 0.4

void buildcostCG( const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol,Mat* rcostVol);
void buildRightcostCG( const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol );

