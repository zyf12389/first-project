#include<iostream>
#include<cv.h>
#include<highgui.h>
using namespace std;
using namespace cv;
#define BORDER_THRES 0.011764

#define TAU_1 0.02745
#define TAU_2 0.0384
#define ALPHA 0.4
double cost(double* lc,double* rc,double* lg,double* rg)
{
	double cdiff=0.0;
	for(int i=0;i<3;i++)
	{
		cdiff+=fabs(lc[i]-rc[i]);
	}
	cdiff/=3;
	double gdiff=fabs(*lg-*rg);
	cdiff=cdiff>TAU_1?TAU_1:cdiff;
	gdiff=gdiff>TAU_2?TAU_2:gdiff;
	return ALPHA*cdiff+(1-ALPHA)*gdiff;
}
double cost(double*lc,double* lg)
{
	double cdiff=0.0;
	for(int i=0;i<3;i++)
	{
		cdiff+=fabs(lc[i]-BORDER_THRES);
	}
	cdiff/=3;
	double gdiff=fabs(lg[0]-BORDER_THRES);
	cdiff=cdiff>TAU_1?TAU_1:cdiff;
	gdiff=gdiff>TAU_2?TAU_2:gdiff;
	return ALPHA*cdiff + (1 - ALPHA)*gdiff;
}
void computeCost(const Mat &limg,const Mat &rimg,const int maxDisp,Mat* costVol)
{
	CV_Assert( limg.type() == CV_64FC3 && rimg.type() == CV_64FC3 );
	int hei=limg.rows;
	int wid=limg.cols;
	Mat lgray,rgray;
	Mat lgrdx,rgrdx;
	Mat temp;
	limg.convertTo(temp,CV_32F);
	cvtColor(temp,lgray,CV_RGB2GRAY);
    rimg.convertTo(temp,CV_32F);
	cvtColor(temp,rgray,CV_RGB2GRAY);
	Sobel(lgray,lgrdx,CV_64F,1,0,1);
	Sobel(rgray,rgrdx,CV_64F,1,0,1);
	lgrdx+=0.5;
	rgrdx+=0.5;
	for(int d=0;d<maxDisp;d++)
	{
		cout<<"-L-C-";
		for(int i=0;i<hei;i++)
		{
			double *lptr=(double*)limg.ptr<double>(i);
			double *rptr=(double*)rimg.ptr<double>(i);
			double *lgptr=(double*)lgrdx.ptr<double>(i);
			double *rgptr=(double*)rgrdx.ptr<double>(i);
			double* c=(double*)costVol[d].ptr<double>(i);
			for(int j=0;j<wid;j++)
			{
				if(j-d>=0)
				{
					double *lc=lptr+3*j;
					double *rc=rptr+3*(j-d);
					double *lg=lgptr+j;
					double *rg=rgptr+j-d;
					c[j]=cost(lc,rc,lg,rg);
				}
				else
				{
				    double *lc=lptr+3*j;
					double *lg=lgptr+j;
					c[j]=cost(lc,lg);
				}
			}
		}
	}
}
void computeCostRight(const Mat &limg, const Mat &rimg, const int maxDisp, Mat* rcostVol)
{
	CV_Assert(limg.type() == CV_64FC3 && rimg.type() == CV_64FC3);
	int hei = limg.rows;
	int wid = limg.cols;
	Mat lgray, rgray;
	Mat lgrdx, rgrdx;
	Mat temp;
	limg.convertTo(temp, CV_32F);
	cvtColor(temp, lgray, CV_RGB2GRAY);
	rimg.convertTo(temp, CV_32F);
	cvtColor(temp, rgray, CV_RGB2GRAY);
	Sobel(lgray, lgrdx, CV_64F, 1, 0, 1);
	Sobel(rgray, rgrdx, CV_64F, 1, 0, 1);
	lgrdx += 0.5;
	rgrdx += 0.5;
	for (int d = 0; d<maxDisp; d++)
	{
		cout << "-R-C-";
		for (int i = 0; i<hei; i++)
		{
			double *lptr = (double*)limg.ptr<double>(i);
			double *rptr = (double*)rimg.ptr<double>(i);
			double *lgptr = (double*)lgrdx.ptr<double>(i);
			double *rgptr = (double*)rgrdx.ptr<double>(i);
			double* c = (double*)rcostVol[d].ptr<double>(i);
			for (int j = 0; j<wid; j++)
			{
				if (j + d <wid)
				{
					double *lc = lptr + 3 * (j+d);
					double *rc = rptr + 3 * j;
					double *lg = lgptr + j+d;
					double *rg = rgptr + j;
					c[j] = cost(lc, rc, lg, rg);
				}
				else
				{
					double *rc = rptr + 3 * j;
					double *rg = rgptr + j;
					c[j] = cost(rc, rg);
				}
			}
		}
	}
}
