#include"sadcc.h"
#include<bitset>
double sadvalue(const Mat& src1, const Mat& src2)
{
	Mat diff;
	absdiff(src1,src2,diff);
	Scalar saddiff = sum(diff);
	double avgdiff = (saddiff.val[2] + saddiff.val[1] + saddiff.val[0]) / 3;
	return avgdiff;
}

void buildcostSAD(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol, Mat* rcostVol)
{
	int half_size = WIN_SIZE / 2;
	int hei = lImg.rows;
	int wid = lImg.cols;
	 Mat limg;
	 Mat rimg;
	copyMakeBorder(lImg, limg, half_size, half_size, half_size, half_size, BORDER_REPLICATE);
	copyMakeBorder(rImg, rimg, half_size, half_size, half_size, half_size, BORDER_REPLICATE);
	cout << limg.size() << endl;
	Mat *lcost = new Mat[maxDis];
	Mat *rcost = new Mat[maxDis];
	for (int i = 0; i < maxDis; i++)
	{
		lcost[i] = Mat::zeros(hei+2*half_size, wid+2*half_size,CV_64FC1);
		rcost[i] = Mat::zeros(hei + 2 * half_size, wid + 2 * half_size, CV_64FC1);
	}
	for (int i = 0; i < maxDis; i++)
	{
		copyMakeBorder(costVol[i], lcost[i], half_size, half_size, half_size, half_size, BORDER_REPLICATE);
		copyMakeBorder(rcostVol[i], rcost[i], half_size, half_size, half_size, half_size, BORDER_REPLICATE);
	}
	cout << lcost[3].size() << endl;
	cout << rcost[3].size() << endl;
	//left
	for (int d = 0; d < maxDis; d++)
	{
		for (int h = half_size; h < limg.rows-half_size; h++)
		{
			double* cost = (double*)lcost[d].ptr<double>(h);
			for (int w = half_size; w <limg.cols-half_size; w++)
			{
				Mat leftWin = limg(Range(h - half_size, h + half_size), Range(w - half_size, w + half_size));
				if ((w - half_size - d) >= 0){
					Mat rightWin = rimg(Range(h - half_size, h + half_size), Range(w - half_size - d, w + half_size-d));
					cost[w] = sadvalue(leftWin,rightWin);
				}
				else{
					Mat rightWin = rimg(Range(h - half_size, h + half_size), Range(w - half_size, w + half_size));
					cost[w] = sadvalue(leftWin,rightWin);
				}
			}
		}
	}
	for (int i = 0; i < maxDis; i++)
		costVol[i] = lcost[i](Range(half_size, limg.rows-half_size), Range(half_size, limg.cols-half_size));
	cout << costVol[3].size() << endl;
	//right
	for (int d = 0; d < maxDis; d++)
	{
		for (int h = half_size; h < rimg.rows-half_size; h++)
		{
			double* cost = (double*)rcost[d].ptr<double>(h);
			for (int w = half_size; w < rimg.cols-half_size; w++)
			{
				Mat rightWin = rimg(Range(h - half_size, h + half_size), Range(w - half_size, w + half_size));
				if ((w + half_size + d) <=limg.cols-1 ){
					Mat leftWin = limg(Range(h - half_size, h + half_size), Range(w - half_size + d, w + half_size + d));
					cost[w] = sadvalue(leftWin, rightWin);
				}
				else{
					Mat leftWin = limg(Range(h - half_size, h + half_size), Range(w - half_size, w + half_size));
					cost[w] = sadvalue(leftWin, rightWin);
				}
			}
		}
	}
	for (int i = 0; i < maxDis; i++)
		rcostVol[i] = rcost[i](Range(half_size, rimg.rows-half_size), Range(half_size, rimg.cols-half_size));
	cout << rcostVol[3].size() << endl;
	
}