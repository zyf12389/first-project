//#include"filter.h"
#include"CANLC\NLCCA.h"
//#include"CAST/STCA.h"
//#include"CenCC.h"
#include"CGCC.h"
//#include"sadcc.h"
#define DOUBLE_MAX 1e10
void match(const Mat& limg,const Mat& rimg,const int maxDis,const int multiple,Mat* costVol,Mat *rcostVol,Mat& ldisp,Mat& rdisp)
{
	cout<<"Begin Match...\n";
//	buildCV(limg, rimg, maxDis, costVol);
	//buildRightCV(limg, rimg, maxDis, rcostVol);
	buildcostCG(limg,rimg,maxDis,costVol,rcostVol);
	//buildcostSAD(limg, rimg, maxDis, costVol, rcostVol);
	//buildRightcostCG(limg, rimg, maxDis, rcostVol);
	//computeCost(limg, rimg, maxDis, costVol);
	//computeCostRight(limg, rimg, maxDis, rcostVol);
	cout<<"after CC\n";
	costAggregation(limg,rimg,maxDis,costVol);
	costAggregationRight(limg,rimg,maxDis,rcostVol);
	//costAggregationST(limg, rimg, maxDis, costVol);
	//costAggregationRightST(limg, rimg, maxDis, rcostVol);
	cout<<"after CA\n";
	cout << "left matching....\n";

	for(int i=0;i<limg.rows;i++)
	{
		uchar* lidstdata=(uchar*)ldisp.ptr<uchar>(i);
		for(int j=0;j<limg.cols;j++)
		{
			double minCost=DOUBLE_MAX;
			int mincol=0;
			for(int d=1;d<maxDis;d++)
			{
				double*cost=(double*)costVol[d].ptr<double>(i);
				if(cost[j]<minCost)
				{
					minCost=cost[j];
					mincol=d;
				}
			}
			lidstdata[j]=mincol*multiple;
		}
	}
	cout << "right matching...\n";
	for (int i = 0; i<rimg.rows; i++)
	{
		uchar* ridstdata = (uchar*)rdisp.ptr<uchar>(i);
		for (int j = 0; j<rimg.cols; j++)
		{
			double minCost = DOUBLE_MAX;
			int mincol = 0;
			for (int d = 1; d<maxDis; d++)
			{
				double*cost = (double*)rcostVol[d].ptr<double>(i);
				if (cost[j]<minCost)
				{
					minCost = cost[j];
					mincol = d;
				}
			}
			ridstdata[j] = mincol * multiple;
		}
	}
}