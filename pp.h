#ifndef PP_H
#define PP_H
#include<iostream>
#include<cv.h>
#include<highgui.h>
using namespace std;
using namespace cv;
#define MED_SZ 13
#define SIG_CLR 0.1
#define SIG_DIS 6

void postProcess(Mat &limg,Mat& rimg,Mat &ldisp,Mat &rdisp,const int maxDisp,const int multiple,Mat* costVol,Mat* rcostVol);
void fgs_filter(const Mat& limg, const Mat& rimg, const Mat &ldisp, const Mat &rdisp, const int maxDisp, const int multiple, Mat &result, Mat&result_r);//上采样使用
void onlyDetect(Mat& ldisp, Mat &rdisp, const int multiple);
void upSampling(const Mat &guidance_img, Mat &disp, Mat* costVol,const int maxDisp,const int multiple,const int scale_rate,Mat& result_disp);

#endif