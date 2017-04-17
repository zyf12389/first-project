#pragma once
#include<iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
#define DIFF 0.05

void doFusion(const Mat&lImg, const Mat&rImg, const Mat &ldisp, const Mat& rdisp, const int maxDisc,const int multiple,Mat &resultL,Mat &resultR);
void doFusion2(const Mat&lImg, const Mat&rImg, const Mat &ldisp, const Mat& rdisp, const int maxDisc, const int multiple, Mat &resultL, Mat &resultR);
void improve(const Mat&lImg, Mat& result);