/******************************************************************************************
\Author:	Qingxiong Yang
\Function:	Disparity map upsampling
\Reference:	Qingxiong Yang, Stereo Matching Using Tree Filtering, PAMI 2014.
*******************************************************************************************/
#ifndef QX_TREE_UPSAMPLING_H
#define QX_TREE_UPSAMPLING_H
#include "qx_tree_filter.h"

class qx_tree_upsampling
{
public:
	qx_tree_filter m_tf;
	double***m_cost_vol, ***m_cost_vol_temp;
    qx_tree_upsampling();
    ~qx_tree_upsampling();
    void clean();
	int init(int h,int w,int nr_plane,
		double sigma_range=0.1
		);
	int build_minimum_spanning_tree(unsigned char***guidance_image);
	int disparity_upsampling(double**disparity);//Don't use unsigned char disparity values
private:
	int	m_h,m_w,m_nr_plane; double m_sigma_range;

	unsigned char**m_disparity,**m_disparity_mf;
};
#endif