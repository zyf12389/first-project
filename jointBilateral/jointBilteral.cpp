#include"jointBilateral.h"

void jointBilateral(const Mat&color, Mat&disparity,const int radius,const int multiple,const int scale_rate, double sigma_c,  double sigma_s, Mat&dst)
{
	if (!color.data || !disparity.data)
	{
		cerr << "can't load image!\n";
		cin.get();
		exit(-1);
	}
	CV_Assert(color.type() == CV_8UC3 && disparity.type() == CV_8UC1);
	if (sigma_c <= 0.0)
		sigma_c = 1.0;
	if (sigma_s <= 0.0)
		sigma_s = 1.0;
	double sigma_c_coeff = -0.5 / (sigma_c*sigma_c);
	double sigma_s_coeff = -0.5 / (sigma_s*sigma_s);
	Mat disparity_Factor;
	disparity /= multiple;
	resize(disparity, disparity_Factor, Size(color.cols, color.rows),0,0,2);
	Mat colorBorder, disparityBorder,dstBorder;
	copyMakeBorder(color, colorBorder, radius, radius, radius, radius, BORDER_REPLICATE);
	copyMakeBorder(disparity_Factor, disparityBorder, radius, radius, radius, radius, BORDER_REPLICATE);
	copyMakeBorder(dst, dstBorder, radius, radius, radius, radius, BORDER_REPLICATE);
	Size size = colorBorder.size();

	for (int i = radius; i < size.height - radius; i++)
	{
		uchar* d_data = (uchar*)disparityBorder.ptr<uchar>(i);
		uchar* c_data = (uchar*)colorBorder.ptr<uchar>(i);
		uchar* dst_data = (uchar*)dstBorder.ptr<uchar>(i);
		for (int j = radius; j < size.width - radius; j++)
		{
			double sumwgt = 0.0, sum = 0.0;
			uchar* c = c_data + 3 * j;
			int b = c[0], g = c[1], r = c[2];
			for (int rx = -radius; rx <=radius; rx++)
			{
				uchar* c_data_2 = (uchar*)colorBorder.ptr<uchar>(i + rx);
				uchar* d_data_2 = (uchar*)disparityBorder.ptr<uchar>(i + rx);
				for (int ry = -radius; ry <=radius; ry++)
				{
					uchar* c2 = c_data_2 + 3 * (j + ry);
					int b2 = c2[0], g2 = c2[1], r2 = c2[2];
					int wx = rx*rx + ry*ry;
					double ws = exp(sigma_s_coeff*wx);
					int colordiff = (abs(b - b2) + abs(g - g2) + abs(r - r2));
					double wc = exp(sigma_c_coeff*colordiff*colordiff);
					sumwgt += (int(d_data_2[j+ry]))*ws*wc;
					sum += (ws*wc);
				}
			}
			sum = 1. / sum;
			int value = cvRound(sumwgt*sum);
			dst_data[j] = uchar(value);
		}
	}
	int factor = 1 << scale_rate;
	dst = dstBorder(Range(radius, size.height - radius), Range(radius, size.width - radius));
	dst = dst*factor;
	cout << dst.size() << endl;
}