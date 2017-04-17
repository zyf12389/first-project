#include"computeDepth.h"
#include<fstream>
#include<iostream>
using namespace std;
inline void changeZero(Mat img, int m, int n)
{
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			int m1 = m + i;
			int n1 = n + j;
			if (m1>= 0 && m1< img.rows&&n1>= 0 && n1< img.cols)
			{
				if (img.at<uchar>(m1, n1))
				{
					img.at<uchar>(m, n) = img.at<uchar>(m1, n1);
					break;
				}
			}
		}
	}
}
void computeDepth(Mat& ldisp,  const double Tx, const double f, const double pixelSize, Mat&depth_W)
{
	ofstream out("w.txt");
	for (int i = 0; i < ldisp.rows; i++)
	{
		uchar* ldata = (uchar*)ldisp.ptr<uchar>(i);
		for (int j = 0; j < ldisp.cols; j++)
		{
			if (!ldata[j])//视差值存在为0的情况
			{
				cout << "5555555555555555555555555555555555555555555555555555555555\n";
				changeZero(ldisp, i, j);
			}
			depth_W.at<double>(i, j) = (Tx*f) / double(ldata[j]);
			out << depth_W.at<double>(i, j) << ' ';
		}
		out << endl;
	}

	depth_W=depth_W / pixelSize;
	normalize(depth_W, depth_W, 1.0, 0.0, NORM_MINMAX);
}
void computeDepth_T(Mat& rdisp, const double Tx, const double f, const double pixelSize, Mat &depth_T)
{
	ofstream out("t.txt");
	for (int i = 0; i < rdisp.rows; i++)
	{
		uchar* rdata = (uchar*)rdisp.ptr<uchar>(i);
		for (int j = 0; j < rdisp.cols; j++)
		{
			if (!rdata[j])
			{
				cout << "6666666666666666666666666666666666666666666666666666\n";
				changeZero(rdisp, i, j);
			}
			depth_T.at<double>(i, j) = (Tx*f) /double( rdata[j]);
			out << depth_T.at<double>(i, j) << ' ';
		}
		out << endl;
	}
	depth_T= depth_T / pixelSize;
	normalize(depth_T, depth_T, 1.0, 0.0, NORM_MINMAX);
}