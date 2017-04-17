#include"doFusion.h"
#include<fstream>
inline double colorDiff(double* c1, double* c2)
{
	double diff = 0.0;
	for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
	{
		diff += abs(c1[i] - c2[j]);
	}
	return diff/3.0;
}
void improve(const Mat&lImg, Mat& result)
{
	int hei = lImg.rows;
	int wid = lImg.cols;
	int num = 0;
	ofstream out("diff.txt");
	for (int h = 0; h < hei; h++)
	{
		double* ldata = (double*)lImg.ptr<double>(h);
		double* resultdata = (double*)result.ptr<double>(h);
		for (int w = 0; w < wid; w++)
		{
			double *data = ldata + 3 * w;
			double* rdata = resultdata + 3 * w;
			double diff = colorDiff(rdata, data);
			out << diff << ' ';
			if (diff>DIFF)
			{
				num++;
				for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
				{
					rdata[i] = data[j];
				}
			}
		}
		out << endl;
	}
	cout << num << endl;
}
void doFusion(const Mat&lImg, const Mat&rImg, const Mat &ldisp, const Mat& rdisp,const int maxDisc, const int multiple, Mat &resultL,Mat &resultR)
{
	int num = 0;
	int hei = ldisp.rows;
	int wid = ldisp.cols;
	for (int h = 0; h < hei; h++)
	{
		double *rdata = (double*)rImg.ptr<double>(h);
		double *ldata = (double*)lImg.ptr<double>(h);
		uchar*ldispdata = (uchar*)ldisp.ptr<uchar>(h);
		uchar* rdispdata = (uchar*)rdisp.ptr<uchar>(h);
		double* resultdata = (double*)resultL.ptr<double>(h);
		double* resultdata_R = (double*)resultR.ptr<double>(h);
		for (int w = 0; w < wid; w++)
		{
			double *data = rdata + 3 * w;
			double* data_L = ldata + 3 * w;
			int offset = rdispdata[w] / multiple;
			int offset_L = ldispdata[w] / multiple;
			if (w + offset <= wid - 1){
				double* l = resultdata + 3 * (w + offset);
				for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
				{
					l[i] = data[j];
				}
			}
			if (w - offset_L >= 0)
			{
				double* r = resultdata_R + 3 * (w - offset_L);
				for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
				{
					r[i] = data_L[j];
				}
			}
		}
	}
	//for (int h = 0; h < hei; h++)
	//{
	//	double* ldata = (double*)lImg.ptr<double>(h);
	//	double* resultdata = (double*)result.ptr<double>(h);
	//	for (int w = 0; w < wid; w++)
	//	{
	//		double* rdata = resultdata + 3 * w;
	//		double* data = ldata + 3 * w;
	//		if (rdata[0] && rdata[1] && rdata[2])
	//			continue;

	//			num++;
	//			for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
	//			{
	//				rdata[i] = data[j];
	//			}
	//	}
	//}
	//cout << num << endl;
//	improve(lImg,result);
	
}
void doFusion2(const Mat&lImg, const Mat&rImg, const Mat &ldisp, const Mat& rdisp, const int maxDisc, const int multiple, Mat &resultL, Mat &resultR)
{
	int hei = ldisp.rows;
	int wid = ldisp.cols;
#pragma omp parallel for
	for (int h = 0; h < hei; h++)
	{
		double *rdata = (double*)rImg.ptr<double>(h);
		double *ldata = (double*)lImg.ptr<double>(h);
		uchar*ldispdata = (uchar*)ldisp.ptr<uchar>(h);
		uchar* rdispdata = (uchar*)rdisp.ptr<uchar>(h);
		double* resultdata = (double*)resultL.ptr<double>(h);
		double* resultdata_R = (double*)resultR.ptr<double>(h);
		for (int w = 0; w < wid; w++)
		{
			int offset = rdispdata[w] / multiple;
			int offset_L = ldispdata[w] / multiple;
			if (w + offset <= wid - 1)
			{
				double *r = resultdata_R + 3 * w;
				double* data_L = ldata + 3 * (w+offset);
				for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
				{
					r[i] = data_L[j];
				}
			}
			if (w - offset_L >= 0)
			{
				double *l = resultdata + 3 * w;
				double *data_R = rdata + 3 * (w-offset_L);
				for (int i = 0, j = 2; i < 3 && j >= 0; i++, j--)
				{
					l[i] = data_R[j];
				}
			}
		}
	}
}