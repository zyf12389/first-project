#include"pp.h"
#include"CANLC\qx_nonlocal_cost_aggregation.h"
#include"CANLC\qx_tree_upsampling.h"
#include"CANLC\NLCCA.h"
#include"disparityWLSFilter.h"
#include<fstream>
#include <omp.h>
#define DOUBLE_MAX 1e10
namespace PP
{
	void lrCheck(Mat& ldisp, Mat& rdisp, int *lvaild, int *rvalid, const int multiple)
	{
		int hei = ldisp.rows;
		int wid = ldisp.cols;
		int imgSize = hei*wid;
		memset(lvaild, 0, imgSize*sizeof(int));
		memset(rvalid, 0, imgSize*sizeof(int));
		int *plvalid = lvaild;
		int *prvalid = rvalid;
		for (int i = 0; i < hei; i++)
		{
			uchar* ldisdata = (uchar*)ldisp.ptr<uchar>(i);
			uchar* rdisdata = (uchar*)rdisp.ptr<uchar>(i);
			for (int j = 0; j < wid; j++)
			{
				int ldis = ldisdata[j] / multiple;
				int rloc = (j - ldis + wid) % wid;
				int rdis = rdisdata[rloc] / multiple;
				if (j - ldis >= 0 && abs(ldis - rdis) <= 1 && ldisdata[j]<255 && ldis != 0)
					*plvalid = 1;
				rdis = rdisdata[j] / multiple;
				int loc = (j + rdis + wid) % wid;
				ldis = ldisdata[loc] / multiple;
				if (j + rdis <= wid - 1 && abs(ldis - rdis) <= 1 && rdisdata[j]<255 && rdis != 0)
					*prvalid = 1;
				plvalid++;
				prvalid++;
			}
		}
	}
	void fillInv(Mat& ldisp, Mat &rdisp, int* lvalid, int *rvalid)
	{
		int hei = ldisp.rows;
		int wid = rdisp.cols;
		int *plvalid = lvalid;
		for (int i = 0; i < hei; i++)
		{
			int *ylvalid = lvalid + i*wid;
			uchar* ldisdata = (uchar*)ldisp.ptr<uchar>(i);
			for (int j = 0; j < wid; j++)
			{
				if (*plvalid == 0)
				{
					int lfirst = j;
					int lfind = 0;
					while (lfirst >= 0)
					{
						if (ylvalid[lfirst])
						{
							lfind = 1;
							break;
						}
						lfirst--;
					}
					int rfind = 0;
					int rfirst = j;
					while (rfirst < wid)
					{
						if (ylvalid[rfirst])
						{
							rfind = 1;
							break;
						}
						rfirst++;
					}
					if (lfind&&rfind)
					{
						ldisdata[j] = ldisdata[lfirst] >= ldisdata[rfirst] ? ldisdata[rfirst] : ldisdata[lfirst];
					}
					else if (lfind){
						ldisdata[j] = ldisdata[lfirst];
					}
					else if (rfind){
						ldisdata[j] = ldisdata[rfirst];
					}
				}
				plvalid++;
			}
		}
		int *prvalid = rvalid;
		for (int i = 0; i < hei; i++)
		{
			int *yrvalid = rvalid + i*wid;
			uchar* rdisdata = (uchar*)rdisp.ptr<uchar>(i);
			for (int j = 0; j < wid; j++)
			{
				if (*prvalid == 0)
				{
					int lfirst = j;
					int lfind = 0;
					while (lfirst >= 0)
					{
						if (yrvalid[lfirst])
						{
							lfind = 1;
							break;
						}
						lfirst--;
					}
					int rfind = 0;
					int rfirst = j;
					while (rfirst < wid)
					{
						if (yrvalid[rfirst])
						{
							rfind = 1;
							break;
						}
						rfirst++;
					}
					if (lfind&&rfind)
					{
						rdisdata[j] = rdisdata[lfirst] >= rdisdata[rfirst] ? rdisdata[rfirst] : rdisdata[lfirst];
					}
					else if (lfind){
						rdisdata[j] = rdisdata[lfirst];
					}
					else if (rfind){
						rdisdata[j] = rdisdata[rfirst];
					}
				}
				prvalid++;
			}
		}
	}
	void wgtMedian(const Mat& lImg, const Mat& rImg, Mat& lDis, Mat& rDis, int* lValid, int* rValid, const int maxDis, const int disSc)
	{
		int hei = lDis.rows;
		int wid = lDis.cols;
		int wndR = MED_SZ / 2;
		double* disHist = new double[maxDis];
		// filter left	
		int* pLValid = lValid;
		for (int y = 0; y < hei; y++) {
			uchar* lDisData = (uchar*)lDis.ptr<uchar>(y);
			float* pL = (float*)lImg.ptr<float>(y);
			for (int x = 0; x < wid; x++) {
				if (*pLValid == 0) {
					// just filter invalid pixels
					memset(disHist, 0, sizeof(double)* maxDis);
					double sumWgt = 0.0f;
					// set disparity histogram by bilateral weight
					for (int wy = -wndR; wy <= wndR; wy++) {
						int qy = (y + wy + hei) % hei;
						// int* qLValid = lValid + qy * wid;
						float* qL = (float*)lImg.ptr<float>(qy);
						uchar* qDisData = (uchar*)lDis.ptr<uchar>(qy);
						for (int wx = -wndR; wx <= wndR; wx++) {
							int qx = (x + wx + wid) % wid;
							// invalid pixel also used
							// if( qLValid[ qx ] && wx != 0 && wy != 0 ) {
							int qDep = qDisData[qx] / disSc;
							if (qDep != 0) {

								double disWgt = wx * wx + wy * wy;
								disWgt = sqrt(disWgt);
								double clrWgt = (pL[3 * x] - qL[3 * qx]) * (pL[3 * x] - qL[3 * qx]) +
									(pL[3 * x + 1] - qL[3 * qx + 1]) * (pL[3 * x + 1] - qL[3 * qx + 1]) +
									(pL[3 * x + 2] - qL[3 * qx + 2]) * (pL[3 * x + 2] - qL[3 * qx + 2]);
								clrWgt = sqrt(clrWgt);
								double biWgt = exp(-disWgt / (SIG_DIS * SIG_DIS) - clrWgt / (SIG_CLR * SIG_CLR));
								disHist[qDep] += biWgt;

								sumWgt += biWgt;
							}
						}
					}
					double halfWgt = sumWgt / 2.0f;
					sumWgt = 0.0f;
					int filterDep = 0;
					for (int d = 0; d < maxDis; d++) {
						sumWgt += disHist[d];
						if (sumWgt >= halfWgt) {
							filterDep = d;
							break;
						}
					}
					// set new disparity
					lDisData[x] = filterDep * disSc;
				}
				pLValid++;
			}
		}
		// filter right depth
		int* pRValid = rValid;
		for (int y = 0; y < hei; y++) {
			uchar* rDisData = (uchar*)rDis.ptr<uchar>(y);
			float* pR = (float*)rImg.ptr<float>(y);
			for (int x = 0; x < wid; x++) {
				if (*pRValid == 0) {
					// just filter invalid pixels
					memset(disHist, 0, sizeof(double)* maxDis);
					double sumWgt = 0.0f;
					// set disparity histogram by bilateral weight
					for (int wy = -wndR; wy <= wndR; wy++) {
						int qy = (y + wy + hei) % hei;
						// int* qRValid = rValid + qy * wid;
						float* qR = (float*)rImg.ptr<float>(qy);
						uchar* qDisData = (uchar*)rDis.ptr<uchar>(qy);
						for (int wx = -wndR; wx <= wndR; wx++) {
							int qx = (x + wx + wid) % wid;
							// if( qRValid[ qx ] && wx != 0 && wy != 0 ) {
							int qDep = qDisData[qx] / disSc;
							if (qDep != 0) {

								double disWgt = wx * wx + wy * wy;
								disWgt = sqrt(disWgt);
								double clrWgt =
									(pR[3 * x] - qR[3 * qx]) * (pR[3 * x] - qR[3 * qx]) +
									(pR[3 * x + 1] - qR[3 * qx + 1]) * (pR[3 * x + 1] - qR[3 * qx + 1]) +
									(pR[3 * x + 2] - qR[3 * qx + 2]) * (pR[3 * x + 2] - qR[3 * qx + 2]);
								clrWgt = sqrt(clrWgt);
								double biWgt = exp(-disWgt / (SIG_DIS * SIG_DIS) - clrWgt / (SIG_CLR * SIG_CLR));
								disHist[qDep] += biWgt;
								sumWgt += biWgt;
							}
							// }

						}
					}
					double halfWgt = sumWgt / 2.0f;
					sumWgt = 0.0f;
					int filterDep = 0;
					for (int d = 0; d < maxDis; d++) {
						sumWgt += disHist[d];
						if (sumWgt >= halfWgt) {
							filterDep = d;
							break;
						}
					}
					// set new disparity
					rDisData[x] = filterDep * disSc;
				}
				pRValid++;
			}
		}

		delete[] disHist;
	}
	void saveChk(Mat &ldisp, Mat& rdisp, const int hei, const int wid, int*lvalid, int*rvalid, bool flag)
	{
		Mat lchk = Mat::zeros(hei, wid, CV_8UC3);
		Mat rchk = Mat::zeros(hei, wid, CV_8UC3);
		int* pLV = lvalid;
		int* pRV = rvalid;
		for (int i = 0; i < hei; i++)
		{
			uchar* lchkdata = (uchar*)lchk.ptr<uchar>(i);
			uchar* rchkdata = (uchar*)rchk.ptr<uchar>(i);
			uchar* ldispdata = (uchar*)ldisp.ptr<uchar>(i);
			uchar* rdispdata = (uchar*)rdisp.ptr<uchar>(i);
			for (int j = 0; j < wid; j++)
			{
				if (*pLV){
					uchar* l = lchkdata + 3 * j;
					l[0] = ldispdata[j];
					l[1] = ldispdata[j];
					l[2] = ldispdata[j];
				}
				else{
					uchar* l = lchkdata + 3 * j;
					l[2] = 255;
				}
				if (*pRV){
					uchar*r = rchkdata + 3 * j;
					r[0] = rdispdata[j];
					r[1] = rdispdata[j];
					r[2] = rdispdata[j];
				}
				else{
					uchar*r = rchkdata + 3 * j;
					r[2] = 255;
				}
				pLV++;
				pRV++;
			}
		}
		if (flag)
		{
			imwrite("l_chk_after_refine.png", lchk);
			imwrite("r_chk_after_refine.png", rchk);
		}
		else
		{
			imwrite("l_chk_before_refine.png", lchk);
			imwrite("r_chk_before_refine.png", rchk);
		}
	}


	//adcensus refinement
	inline double colorDiff(double* c1, double* c2)
	{
		double diff = 0, maxdiff = 0;
		for (int color = 0; color < 3; color++)
		{
			diff = abs(c1[color] - c2[color]);
			maxdiff = maxdiff < diff ? diff : maxdiff;
		}
		return maxdiff;
	}
	void outlierDetect(const Mat &ldisp,const Mat &rdisp, const int maxDis, const int multiple, Mat& ldisparity, Mat& rdisparity)
	{
		int occ = 0;
		int mismatch = 0;
		int rocc = 0;
		int rmismatch = 0;
		Size size = ldisp.size();
		for (int h = 0; h < ldisp.rows; h++)
		{
			uchar* ldata = (uchar*)ldisp.ptr<uchar>(h);
			uchar* rdata = (uchar*)rdisp.ptr<uchar>(h);
			int* ldispptr = (int*)ldisparity.ptr<int>(h);
			int* rdispptr = (int*)rdisparity.ptr<int>(h);
			for (int w = 0; w < ldisp.cols; w++)
			{
				int disp = ldata[w] / multiple;
				int dispr = rdata[w] / multiple;
				//left
				if (abs(disp - rdata[w - disp] / multiple)>1 || w - disp<0)//or w-d<0
				{
					bool occlusion = true;
					for (int d = 0; d < maxDis&&occlusion; d++)
					{
						if (w - d >= 0 && d == rdata[w - d] / multiple){
							occlusion = false;
						}
					}
					if (occlusion)
						occ++;
					else
						mismatch++;
					disp = (occlusion) ? -1 : -5;
				}
				ldispptr[w] = disp*multiple;
				//right
				if (abs(dispr - ldata[w + dispr] / multiple) > 1 || w + dispr>ldisp.cols - 1)
				{
					bool occlusion = true;
					for (int d = 0; d < maxDis&&occlusion; d++)
					{
						if (w + d <= ldisp.cols - 1 && d == ldata[w + d] / multiple){
							occlusion = false;
						}
					}
					if (occlusion)
						rocc++;
					else
						rmismatch++;
					dispr = (occlusion) ? -1 : -5;

				}

				rdispptr[w] = dispr*multiple;
			}
		}
		cout << occ << ' ' << mismatch << endl;
		cout << rocc << ' ' << rmismatch << endl;
	}
	Mat covertdisp(const Mat &disparity,Mat &disp)
	{
		Size size = disparity.size();
		Mat dispmap(size, CV_8U);
		for (int i = 0; i < size.height; i++)
		{
			uchar* uptr = (uchar*)dispmap.ptr<uchar>(i);
			int *iptr = (int*)disparity.ptr<int>(i);
			for (int j = 0; j < size.width; j++)
			{
				uptr[j] = (iptr[j] < 0) ? disp.at<uchar>(i,j) : uchar(iptr[j]);
			}
		}
		return dispmap;
	}
	Mat confidenceMap(const Mat& disparity,const Mat&disp)
	{
		Size size = disparity.size();
		Mat dispmap(size, CV_32F);
		for (int i = 0; i < size.height; i++)
		{
			float* uptr = (float*)dispmap.ptr<float>(i);
			int *iptr = (int*)disparity.ptr<int>(i);
			for (int j = 0; j < size.width; j++)
			{
				uptr[j] = (iptr[j] < 0) ? float(disp.at<uchar>(i,j)) : 255;
			}
		}
		return dispmap;
	}
	void interpolation(Mat& disparity, const Mat& limg, const int multiple)
	{
		Size size = disparity.size();
		Mat disptemp(size, CV_32S);
		int directionsW[] = { 0, 2, 2, 2, 0, -2, -2, -2, 1, 2, 2, 1, -1, -2, -2, -1 };
		int directionsH[] = { 2, 2, 0, -2, -2, -2, 0, 2, 2, 1, -1, -2, -2, -1, 1, 2 };
		for (int h = 0; h < size.height; h++)
		{
			int* disp = (int*)disparity.ptr<int>(h);
			int* temp = (int*)disptemp.ptr<int>(h);
			double* limgptr = (double*)limg.ptr<double>(h);
			for (int w = 0; w < size.width; w++)
			{
				if (disp[w]>0)
					temp[w] = disp[w];
				else
				{
					vector<int>neighbordisp(16, disp[w]);
					vector<double>neighbordiff(16, -1);

					for (int dir = 0; dir < 16; dir++)
					{
						bool inside = true, gotdisp = false;
						int hd = h, wd = w;
						for (int sd = 0; sd < 8 && inside&&!gotdisp; sd++)
						{
							if (sd % 2 == 0)
							{
								hd += directionsH[dir] / 2;
								wd += directionsW[dir] / 2;
							}
							else
							{
								hd += directionsH[dir] - directionsH[dir] / 2;
								wd += directionsW[dir] - directionsW[dir] / 2;
							}
							inside = (hd >= 0 && hd < size.height&&wd >= 0 && wd < size.width);
							if (inside&&disparity.at<int>(hd, wd)>0)
							{
								double* limgptr2 = (double*)limg.ptr<double>(hd);
								neighbordisp[dir] = disparity.at<int>(hd, wd);
								neighbordiff[dir] = colorDiff(limgptr + 3 * w, limgptr2 + 3 * wd);
								gotdisp = true;
							}
						}
					}
					if (disp[w] == -multiple)
					{
						int mindisp = neighbordisp[0];
						for (int dir = 1; dir < 16; dir++)
						{
							if (mindisp>neighbordisp[dir])
								mindisp = neighbordisp[dir];
						}
						temp[w] = mindisp;
					}
					else
					{
						int mindisp = neighbordisp[0];
						double mindiff = neighbordiff[0];
						for (int dir = 1; dir < 16; dir++)
						{
							if (mindiff<0 || (mindiff>neighbordiff[dir] && neighbordiff[dir]>0))
							{
								mindisp = neighbordisp[dir];
								mindiff = neighbordiff[dir];
							}
						}
						temp[w] = mindisp;
					}
				}

			}
		}
		disptemp.copyTo(disparity);
	}
	Mat subpixel(Mat &disparity, Mat* costVol, const int maxDis, int multiple)
	{
		Size size = disparity.size();
		Mat temp(size, CV_32F);
		for (int h = 0; h < size.height; h++)
		{
			uchar* disp = (uchar*)disparity.ptr<uchar>(h);
			float* tempdata = (float*)temp.ptr<float>(h);
			for (int w = 0; w < size.width; w++)
			{
				int d = disp[w] / multiple;
				double interdisp = d;
				if (d>0 && d < maxDis - 1)
				{
					double cost = costVol[d].at<double>(h, w);
					double costplus = costVol[d + 1].at<double>(h, w);
					double costminus = costVol[d - 1].at<double>(h, w);
					double diff = (costplus - costminus) / (2 * (costplus + costminus - 2 * cost));
					if (diff>-1 && diff < 1)
						interdisp -= diff;
				}
				tempdata[w] = interdisp*multiple;
			}
		}
		medianBlur(temp, temp, 3);
		return temp;
	}

	//non-local
	void nonLocalRe(Mat& limg, Mat &rimg, Mat &ldisp, Mat &rdisp, const int maxDisp, const int multiple, Mat* costVol, Mat* rcostVol)
	{
		double sigma = 0.1;
		int hei = ldisp.rows;
		int wid = rdisp.cols;
		int imgSize = hei*wid;
		int *lvalid = new int[imgSize];
		int *rvalid = new int[imgSize];
		Mat *lcost = new Mat[maxDisp];
		Mat *rcost = new Mat[maxDisp];
		for (int i = 0; i < maxDisp; i++)
		{
			lcost[i] = costVol[i];
			rcost[i] = rcostVol[i];
		}
		unsigned char*** left = qx_allocu_3(hei, wid, 3);//allocate memory
		unsigned char*** right = qx_allocu_3(hei, wid, 3);
		cvtMatQX(limg, left[0][0]);
		cvtMatQX(rimg, right[0][0]);
		qx_nonlocal_cost_aggregation m_nlca, m_nlcaR;
		m_nlca.init(hei, wid, maxDisp, sigma);
		m_nlcaR.init(hei, wid, maxDisp, sigma);
		m_nlca.m_left = left;
		m_nlca.m_right = right;
		m_nlcaR.m_left = left;
		m_nlcaR.m_right = right;
		for (int d = 0; d < maxDisp; d++) {
			for (int y = 0; y < hei; y++) {
				double* costl = (double*)lcost[d].ptr<double>(y);
				double* costr = (double*)rcost[d].ptr<double>(y);
				for (int x = 0; x < wid; x++) {
					m_nlca.m_cost_vol[y][x][d] = costl[x];
					m_nlcaR.m_cost_vol_right[y][x][d] = costr[x];
				}
			}
		}
		m_nlca.m_tf.build_tree(left[0][0]);
		m_nlcaR.m_tf_right.build_tree(right[0][0]);
		lrCheck(ldisp, rdisp, lvalid, rvalid, multiple);
		int *plv = lvalid;
		int *rplv = rvalid;
		for (int h = 0; h < hei; h++)
		for (int w = 0; w < wid; w++)
		{
			if (!*plv)
			{
				for (int d = 0; d < maxDisp; d++)
					m_nlca.m_cost_vol[h][w][d] = abs(d - int(ldisp.at<uchar>(h, w)) / multiple);
			}
			/*if (*rplv == 0)
			{
			for (int d = 0; d < maxDisp; d++)
			m_nlcaR.m_cost_vol_right[h][w][d] = abs(d -int(rdisp.at<uchar>(h, w)) / multiple);
			}*/
			plv++;
			//	rplv++;
		}
		m_nlca.m_tf.update_table(sigma / 2);
		m_nlca.m_tf.filter(m_nlca.m_cost_vol[0][0], m_nlca.m_cost_vol_temp[0][0], maxDisp);
		m_nlcaR.m_tf_right.update_table(sigma / 2);
		m_nlcaR.m_tf_right.filter(m_nlcaR.m_cost_vol_right[0][0], m_nlcaR.m_cost_vol_temp[0][0], maxDisp);

		for (int d = 1; d < maxDisp; d++) {
			for (int y = 0; y < hei; y++) {
				double* costl = (double*)lcost[d].ptr<double>(y);
				double* costr = (double*)rcost[d].ptr<double>(y);
				for (int x = 0; x < wid; x++) {
					costl[x] = m_nlca.m_cost_vol[y][x][d];
					//		costr[x] = m_nlcaR.m_cost_vol_right[y][x][d];

				}
			}
		}
		for (int i = 0; i<limg.rows; i++)
		{
			uchar* lidstdata = (uchar*)ldisp.ptr<uchar>(i);
			uchar* rdistdata = (uchar*)rdisp.ptr<uchar>(i);
			for (int j = 0; j<limg.cols; j++)
			{
				double minCostl = DOUBLE_MAX;
				double minCostr = DOUBLE_MAX;
				int mincoll = 0;
				int mincolr = 0;
				for (int d = 1; d<maxDisp; d++)
				{
					double*costl = (double*)lcost[d].ptr<double>(i);
					double*costr = (double*)rcost[d].ptr<double>(i);
					if (costl[j]<minCostl)
					{
						minCostl = costl[j];
						mincoll = d;
					}
					if (costr[j] < minCostr)
					{
						minCostr = costr[j];
						mincolr = d;
					}
				}
				lidstdata[j] = mincoll*multiple;
				//	rdistdata[j] = mincolr*multiple;
			}
		}
		delete[]lvalid;
		delete[]rvalid;
		qx_freeu_3(left);
		left = NULL;//free memory
		qx_freeu_3(right);
		right = NULL;
	}
}
void onlyDetect(Mat& ldisp, Mat &rdisp, const int multiple)
{
	bool flag = false;
	int hei = ldisp.rows;
	int wid = ldisp.cols;
	int *lvalid = new int[hei*wid];
	int* rvalid = new int[hei*wid];
	PP::lrCheck(ldisp, rdisp, lvalid, rvalid, multiple);
	PP::saveChk(ldisp, rdisp, hei, wid, lvalid, rvalid, flag);
	delete[]lvalid;
	delete[]rvalid;
}

void fgs_filter(const Mat& limg,const Mat&rimg,const Mat &ldisp,const Mat &rdisp, const int maxDisp, const int multiple,Mat& result,Mat&result_r)
{
	int wid = limg.cols;
	int hei = limg.rows;
	result = Mat(Size(wid, hei), CV_8UC1);
	Mat ldisparity(Size(ldisp.cols, ldisp.rows), CV_32S);
	Mat rdisparity(Size(rdisp.cols, rdisp.rows), CV_32S);
	Mat llast(Size(wid, hei), CV_8UC1);
	Mat rlast(Size(wid, hei), CV_8UC1);

	PP::outlierDetect(ldisp, rdisp, maxDisp, multiple, ldisparity, rdisparity);
	Mat confidence_map_1, confidence_map_right_1;
	confidence_map_1 = PP::confidenceMap(ldisparity,ldisp);
	confidence_map_right_1 = PP::confidenceMap(rdisparity,rdisp);

	resize(ldisp, llast, Size(wid, hei), 0);
	resize(rdisp, rlast, Size(wid, hei), 0);
	Mat confidence_map = Mat(Size(wid, hei), CV_32F);
	Mat confidence_map_right = Mat(Size(wid, hei), CV_32F);
	resize(confidence_map_1, confidence_map, Size(wid, hei), 0);
	resize(confidence_map_right_1, confidence_map_right, Size(wid, hei), 0);

	Ptr<DisparityWLSFilter> wls_filter = createDisparityWLSFilterGeneric(true);
	wls_filter->setMultiple(multiple);
	wls_filter->setLambda(9000.0);
	wls_filter->setLambda_a(1.0);
	wls_filter->setSigmaColor(1.5);
	wls_filter->setIter(30);
	wls_filter->setConfidenceMap(confidence_map);
	wls_filter->setConfidenceMap_right(confidence_map_right);
	wls_filter->filter(llast, limg, llast, rlast, rlast, Rect(), rimg);
	llast.convertTo(result, CV_8UC1);
	rlast.convertTo(result_r, CV_8UC1);
}

void postProcess(Mat& limg, Mat &rimg, Mat &ldisp, Mat &rdisp, const int maxDisp, const int multiple, Mat* costVol, Mat* rcostVol)
{
	int hei = ldisp.rows;
	int wid = rdisp.cols;
	int imgSize = hei*wid;
	int *lvalid = new int[imgSize];
	int *rvalid = new int[imgSize];
	Mat ltemp, rtemp;
	limg.convertTo(ltemp, CV_32F);
	rimg.convertTo(rtemp, CV_32F);
	cout << "post processing...\n";
	Mat ldisparity(Size(wid, hei), CV_32S);
	Mat rdisparity(Size(wid, hei), CV_32S);
	PP::outlierDetect(ldisp, rdisp, maxDisp, multiple, ldisparity, rdisparity);

	//non-local refinment

//	PP::nonLocalRe(limg, rimg, ldisp, rdisp, maxDisp, multiple, costVol, rcostVol);

	//non-local refinement end

	for (int i = 0; i < 3; i++)
	{
		PP::lrCheck(ldisp, rdisp, lvalid, rvalid, multiple);
		PP::fillInv(ldisp, rdisp, lvalid, rvalid);
		PP::wgtMedian(ltemp, rtemp, ldisp, rdisp, lvalid, rvalid, maxDisp, multiple);
	}
	delete[] lvalid;
	delete[] rvalid;

	PP::outlierDetect(ldisp, rdisp, maxDisp, multiple, ldisparity, rdisparity);
	for (int i = 0; i < 64; i++){
		PP::interpolation(ldisparity, limg, multiple);
		PP::interpolation(rdisparity, rimg, multiple);
	}
	Mat ld(Size(wid, hei), CV_8U);
	Mat rd(Size(wid, hei), CV_8U);
	ld = PP::covertdisp(ldisparity,ldisp);
	rd = PP::covertdisp(rdisparity,rdisp);

	int *lvalid2 = new int[imgSize];
	int *rvalid2 = new int[imgSize];
	for (int i = 0; i < 3; i++)
	{
		PP::lrCheck(ld, rd, lvalid2, rvalid2, multiple);
		PP::fillInv(ld, rd, lvalid2, rvalid2);
		PP::wgtMedian(ltemp, rtemp, ld, rd, lvalid2, rvalid2, maxDisp, multiple);
	}
	ld.copyTo(ldisp);
	rd.copyTo(rdisp);

	//sub-pixel interpolation
	Mat llast(Size(wid, hei), CV_32F);
	Mat rlast(Size(wid, hei), CV_32F);
	PP::subpixel(ld, costVol, maxDisp, multiple).copyTo(llast);
	PP::subpixel(rd, rcostVol, maxDisp, multiple).copyTo(rlast);


	//FGS Filter
	PP::outlierDetect(ldisp, rdisp, maxDisp, multiple, ldisparity, rdisparity);
	Mat confidence_map, confidence_map_right;
	confidence_map = PP::confidenceMap(ldisparity,ldisp);
	confidence_map_right = PP::confidenceMap(rdisparity,rdisp);
	Mat lImg, rImg;
	limg.convertTo(lImg, CV_8UC3, 255);
	rimg.convertTo(rImg, CV_8UC3, 255);
	Ptr<DisparityWLSFilter> wls_filter = createDisparityWLSFilterGeneric(true);
	wls_filter->setMultiple(multiple);
	wls_filter->setLambda(8000.0);
	wls_filter->setLambda_a(1.0);
	wls_filter->setSigmaColor(1.5);
	wls_filter->setIter(20);
	wls_filter->setConfidenceMap(confidence_map);
	wls_filter->setConfidenceMap_right(confidence_map_right);
	wls_filter->filter(llast, lImg, llast, rlast, rlast, Rect(), rImg);
	llast.convertTo(ldisp, CV_8UC1);
	rlast.convertTo(rdisp, CV_8UC1);

	PP::outlierDetect(ldisp, rdisp, maxDisp, multiple, ldisparity, rdisparity);
	bool flag = true;
	PP::lrCheck(ldisp, rdisp, lvalid2, rvalid2, multiple);
	PP::saveChk(ldisp, rdisp, hei, wid, lvalid2, rvalid2, flag);
	delete[]lvalid2;
	delete[]rvalid2;
}