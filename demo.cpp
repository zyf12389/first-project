#include<time.h>
#include<omp.h>
#include<fstream>
#include"match.h"
#include"pp.h"
#include"Fusion\doFusion.h"
#include"rectify\rectify.h"
#include"computeDepth\computeDepth.h"
#include"upsample\fgs_upsample.h"
int main()
{
	clock_t start = clock();
	int maxDisc = 16;//最大视差范围
	int multiple = 4;//放大倍数
	for (int i = 1; i <= 1; i++)
	{
		char filename1[256];
		char filename2[256];
		char filepath1[256];
		char filepath2[256];
		if (i < 10)
		{
			sprintf(filename1, "%s%d%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\", i, "_W.jpeg");
			sprintf(filename2, "%s%d%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\", i, "_T.jpeg");
		}
		else
		{
			sprintf(filename1, "%s%d%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\", i, "_W.jpeg");
			sprintf(filename2, "%s%d%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\", i, "_T.jpeg");
		}
		cout << filename1 << endl;
		double pixelSize = 1.4;
		Mat img_T = imread(filename2);
		Mat img_W = imread(filename1);
		string filepath = "dual_cal1.txt";
		Mat limg, rimg;
		double Tx, f;
		rectify(img_T, img_W, filepath, rimg, limg, Tx, f);//立体校正

		if (!limg.data || !rimg.data)
		{
			printf("Error: can not open image\n");
			printf("\nPress any key to continue...\n");
			getchar();
			return -1;
		}
		int hei = limg.rows;
		int wid = limg.cols;
		int scale_rate;
		if (hei / 4 < 500)
		{
			scale_rate = 2;
		}
		else
		{
			scale_rate = 3;
		}


		Mat limg_copy = limg;
		Mat rimg_copy = rimg;

		Mat lImg = limg, rImg = rimg;
		Mat limg1 = limg, rimg1 = rimg;
		Mat lImg1 = limg1, rImg1 = rimg1;
		for (int j = 0; j < scale_rate; j++){
			pyrDown(limg, lImg, Size(limg.cols / 2, limg.rows / 2));
			pyrDown(rimg, rImg, Size(rimg.cols / 2, rimg.rows / 2));
			limg = lImg;
			rimg = rImg;
		}

		for (int j = 0; j < scale_rate; j++){
			resize(limg1, lImg1, Size(limg1.cols / 2, limg1.rows / 2));
			resize(rimg1, rImg1, Size(rimg1.cols / 2, rimg1.rows / 2));
			limg1 = lImg1;
			rimg1 = rImg1;
		}

		if (i < 10)
		{
			string str = "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00" + to_string(i);
			string str1, str2;
			str1 = str + "\\limg_64.jpg";
			str2 = str + "\\rimg_64.jpg";
			imwrite(str1, lImg1);
			imwrite(str2, rImg1);
		}
		else
		{
			string str = "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0" + to_string(i);
			string str1, str2;
			str1 = str + "\\limg_64.jpg";
			str2 = str + "\\rimg_64.jpg";
			imwrite(str1, lImg1);
			imwrite(str2, rImg1);
		}

		cout << lImg.size() << endl;
		cout << scale_rate << endl;


		cvtColor(lImg, lImg, CV_BGR2RGB);
		cvtColor(rImg, rImg, CV_BGR2RGB);
		cvtColor(lImg1, lImg1, CV_BGR2RGB);
		cvtColor(rImg1, rImg1, CV_BGR2RGB);


		lImg.convertTo(lImg, CV_64F, 1 / 255.0f);
		rImg.convertTo(rImg, CV_64F, 1 / 255.0f);
		lImg1.convertTo(lImg1, CV_64F, 1 / 255.0f);
		rImg1.convertTo(rImg1, CV_64F, 1 / 255.0f);

		Mat lP = lImg.clone();
		Mat rP = rImg.clone();
		Mat ldisp = Mat::zeros(lImg.rows, lImg.cols, CV_8UC1);
		Mat rdisp = Mat::zeros(rImg.rows, rImg.cols, CV_8UC1);
		Mat subldisp = Mat::zeros(lImg.rows, lImg.cols, CV_32F);
		Mat subrdisp = Mat::zeros(rImg.rows, rImg.cols, CV_32F);
		Mat result_L = Mat::zeros(lImg.rows, lImg.cols, CV_64FC3);
		Mat result_R = Mat::zeros(rImg.rows, rImg.cols, CV_64FC3);
		Mat* costVol, *rcostVol;
		costVol = new Mat[maxDisc];
		for (int mindex = 0; mindex < maxDisc; mindex++)
			costVol[mindex] = Mat::zeros(lImg.rows, lImg.cols, CV_64FC1);
		rcostVol = new Mat[maxDisc];
		for (int index = 0; index < maxDisc; index++)
			rcostVol[index] = Mat::zeros(rImg.rows, rImg.cols, CV_64FC1);

		match(lP, rP, maxDisc, multiple, costVol, rcostVol, ldisp, rdisp);

		//只进行左右检测
		onlyDetect(ldisp, rdisp, multiple);
		//后处理
		postProcess(lImg, rImg, ldisp, rdisp, maxDisc, multiple, costVol, rcostVol);
		//视点转换
		doFusion2(lImg1, rImg1, ldisp, rdisp, maxDisc, multiple, result_L, result_R);

		if (i < 10)
		{
			string str = "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00" + to_string(i);
			string str1, str2;
			str1 = str + "\\R2L.jpg";
			str2 = str + "\\L2R.jpg";
			imwrite(str1, result_L * 255);
			imwrite(str2, result_R * 255);
		}
		else
		{
			string str = "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0" + to_string(i);
			string str1, str2;
			str1 = str + "\\R2L.jpg";
			str2 = str + "\\L2R.jpg";
			imwrite(str1, result_L * 255);
			imwrite(str2, result_R * 255);
		}

		
		if (i < 10)
		{
			sprintf(filepath1, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\W_disp_16.jpg");
			sprintf(filepath2, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\T_disp_16.jpg");
		}
		else
		{
			sprintf(filepath1, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\W_disp_16.jpg");
			sprintf(filepath2, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\T_disp_16.jpg");
		}
		cout << filepath1 << endl;
		imwrite(filepath1, ldisp);
		imwrite(filepath2, rdisp);
		memset(filepath1, 0, sizeof(char)* 256);
		memset(filepath2, 0, sizeof(char)* 256);
		//	normalize(ldisp, ldisp, 255, 0, NORM_MINMAX);
		//	normalize(rdisp, rdisp, 255, 0, NORM_MINMAX);
		if (i < 10)
		{
			sprintf(filepath1, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\depth_W.jpg");
			sprintf(filepath2, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\00", i, "\\depth_T.jpg");
		}
		else
		{
			sprintf(filepath1, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\depth_W.jpg");
			sprintf(filepath2, "%s%d%s", "C:\\Users\\ww\\Desktop\\双摄项目\\W-T\\wide-tele\\roi3\\0", i, "\\depth_T.jpg");
		}
		Mat depth_W = Mat::zeros(ldisp.rows, ldisp.cols, CV_64FC1);
		Mat depth_T = Mat::zeros(rdisp.rows, rdisp.cols, CV_64FC1);
		computeDepth(ldisp, Tx, f, pixelSize, depth_W);//计算深度，由于校正截取图像中心变化，还没有考虑图像中心的变化
		computeDepth_T(rdisp, Tx, f, pixelSize, depth_T);
		imwrite(filepath1, depth_W * 255);
		imwrite(filepath2, depth_T * 255);

		//同时进行上采样
		Mat upsample_result,upsample_result_r;
		fgs_upsample(limg_copy, rimg_copy,ldisp, rdisp, scale_rate, maxDisc, multiple, upsample_result,upsample_result_r);
		imwrite("upsample_L.jpg", upsample_result);
		imwrite("upsample_R.jpg", upsample_result_r);
		Mat transfer_L = Mat(Size(limg_copy.cols,limg_copy.rows),CV_64FC3);
		Mat transfer_R = Mat(Size(rimg_copy.cols, rimg_copy.rows),CV_64FC3);
		cvtColor(limg_copy, limg_copy, CV_BGR2RGB);
		cvtColor(rimg_copy, rimg_copy, CV_BGR2RGB);
		limg_copy.convertTo(limg_copy, CV_64F, 1 / 255.0f);
		rimg_copy.convertTo(rimg_copy, CV_64F, 1 / 255.0f);
		int rate = 1 << scale_rate;
		doFusion2(limg_copy,rimg_copy,upsample_result,upsample_result_r,maxDisc*rate,1,transfer_L,transfer_R );
		imwrite("R2L.jpg", transfer_L*255);
		imwrite("L2R.jpg", transfer_R*255);

		clock_t end = clock();
		cout << "Time is:" << (end - start) / CLOCKS_PER_SEC << endl;
		waitKey(-1);
		delete[]costVol;
		delete[]rcostVol;
	}
	system("pause");
	return 0;
}