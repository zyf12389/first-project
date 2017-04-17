#include"fgs_upsample.h"
void fgs_upsample(const Mat& rgb,const Mat& rgb_r,const Mat&disparity_1,const Mat& disparity_2, const int scale_rate,const int maxDis ,int multiple, Mat& result,Mat&result_r)
{
	if (!rgb.data&&!disparity_1.data&&!disparity_2.data)
	{
		cout << "fgs_upsample:can't load image\n";
		system("pause");
	}
	if (0 == multiple)
		multiple = 1;
	int rate = 1 << scale_rate;
	Size size = rgb.size();
	result = Mat::zeros(size, CV_8UC1);
	fgs_filter(rgb,rgb_r,disparity_1,disparity_2, maxDis, multiple,result,result_r);
	result *= rate / multiple;
	result_r *= rate / multiple;
}