#pragma once
#include<iostream>
#include<cv.h>
#include<highgui.h>
using namespace std;
using namespace cv;

#define CENCUS_WND 13
#define CENCUS_BIT 168
#define LAMBDA_CEN 0.5
#define WIN_SIZE 9


void buildcostSAD(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol, Mat* rcostVol);

