#pragma once
#include<iostream>
#include <opencv2/opencv.hpp>
#include<bitset>
using namespace std;
using namespace cv;
#define CENCUS_WND 9
#define CENCUS_BIT 80
void buildCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol);
void buildRightCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol);
