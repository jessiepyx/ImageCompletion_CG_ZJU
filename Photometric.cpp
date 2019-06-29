#include "Photometric.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"
#include <iostream>
#include <cmath>

#define ATD at<Vec3d>
#define ATU at<uchar>
#define M_DST 0
#define M_SRC 1
#define M_BORDER 2
#define M_BOUNDARY 3

int offset_t[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
// local offset: + offset
#define L_OFFSET(i) (y+offset_t[i][0]),(x+offset_t[i][1])
// mask offset: + patch pos + 1 + offset
#define M_OFFSET(i) (y+offset_y+offset_t[i][0]),(x+offset_x+offset_t[i][1])

// init
Mat Photometric::mask;
Mat Photometric::dst;

// init the mask for correction
// only needs calling once
void Photometric::initMask(Mat image, Mat imageMask, uchar unknown, uchar known)
{
	Mat temp;
	// create dst mat
	dst = Mat(imageMask.size().height, imageMask.size().width, CV_64FC3);
	image.convertTo(temp, CV_64FC3);
	temp.copyTo(dst);
	// create mask, +2 is for border
	// mask = Mat(imageMask.size().height + 2, imageMask.size().width + 2, CV_8U);
	// the same size is okay
	mask = Mat(imageMask.size().height, imageMask.size().width, CV_8U);
	// update mask, treat unknown region as border
	mask.setTo(Scalar(unknown));
	//imageMask.copyTo(mask(roi));
	imageMask.copyTo(mask);
	Mat unknown_roi = mask == unknown;
	Mat known_roi = mask == known;
	// M_BORDER will be useless
	mask.setTo(Scalar(M_BORDER), unknown_roi);
	mask.setTo(Scalar(M_DST), known_roi);
	return;
}

void Photometric::correctE(Mat &patch, int offset_x, int offset_y)
{
	// infos
	int width, height, y, x, i, cnt = 0;
	// need preprocessing
	width = patch.size().width;
	height = patch.size().height;
	Rect patch_mask = Rect(offset_x, offset_y, width, height);
	// src: patch with double type
	// result: the modified patch
	Mat patch_d;
	patch.convertTo(patch_d, CV_64FC3);
	Mat result = Mat(height, width, CV_64FC3);
	patch_d.copyTo(result);
	Mat src = Mat(height, width, CV_64FC3);
	patch_d.copyTo(src);
	Mat bitmap = Mat(height, width, CV_8U);
	for (y = 0; y < height; y++)
	{
		for (x = 0; x < width; x++)
		{
			if (mask.ATU(y + offset_y, x + offset_x) == M_DST)
			{
				result.ATD(y, x) = dst.ATD(y + offset_y, x + offset_x);
				src.ATD(y, x) = dst.ATD(y + offset_y, x + offset_x);
				bitmap.ATU(y, x) = M_DST;
			}
			else if (mask.ATU(y + offset_y, x + offset_x) == M_BORDER)
			{
				mask.ATU(y + offset_y, x + offset_x) = M_SRC;
				bitmap.ATU(y, x) = M_SRC;
			}
			if (x == 0 || y == 0 || x == width - 1 || y == height - 1)
			{
				mask.ATU(y + offset_y, x + offset_x) = M_BOUNDARY;
			}
		}
	}
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b[3], sol[3];
	int total = (height - 2)*(width - 2);
	A = Eigen::SparseMatrix<double>(total, total);
	A.reserve(Eigen::VectorXd::Constant(total, 5));
	for (i = 0; i < 3; i++)
	{
		b[i] = Eigen::VectorXd(total);
		sol[i] = Eigen::VectorXd(total);
	}
	// index
	Mat index = Mat(height, width, CV_32S);
	cnt = 0;
	for (y = 0; y < height; y++)
	{
		for (x = 0; x < width; x++)
		{
			if (mask.ATU(y + offset_y, x + offset_x) == M_DST || mask.ATU(y + offset_y, x + offset_x) == M_SRC)
			{
				index.at<int>(y, x) = cnt;
				cnt++;
			}
		}
	}
	// traverse all f_q in patch
	// we know that the patch is a square
	// may using matrix manipulations if i have enough time
	int ch;
	for (y = 1; y < height - 1; y++)
	{
		for (x = 1; x < width - 1; x++)
		{
			for (ch = 0; ch < 3; ch++)
			{
				double sum_vpq = 0, sum_boundary = 0;
				// neighbors
				double neighbor = 0;
				// traverse neighbors
				for (i = 0; i < 4; i++)
				{
					switch (mask.ATU(M_OFFSET(i)))
					{
					case M_BORDER:
						// border, truncated neighborhood
						break;
					case M_BOUNDARY:
						neighbor += 1.0;
						sum_boundary += src.ATD(L_OFFSET(i))(ch);
						if (bitmap.ATU(y, x) == bitmap.ATU(L_OFFSET(i)))
						{
							sum_vpq += src.ATD(y, x)(ch) - src.ATD(L_OFFSET(i))(ch);
						}
						break;
					case M_SRC:
					case M_DST:
						// in region
						if (ch == 0)
						{
							A.insert(index.at<int>(y, x), index.at<int>(L_OFFSET(i))) = -1.0;
						}
						// gradient
						if (mask.ATU(y + offset_y, x + offset_x) == mask.ATU(M_OFFSET(i)))
						{
							sum_vpq += src.ATD(y, x)(ch) - src.ATD(L_OFFSET(i))(ch);
						}
						neighbor += 1.0;
						break;
					}
				}
				if (ch == 0)
				{
					A.insert(index.at<int>(y, x), index.at<int>(y, x)) = neighbor;
				}
				b[ch](index.at<int>(y, x)) = sum_boundary + sum_vpq;
			}
		}
	}
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "decomposition failed" << std::endl;
		return;
	}
	for (ch = 0; ch < 3; ch++)
	{
		sol[ch] = solver.solve(b[ch]);
		if (solver.info() != Eigen::Success)
		{
			std::cout << "solving failed" << std::endl;
			return;
		}
	}
	for (ch = 0; ch < 3; ch++)
	{
		for (y = 1; y < height - 1; y++)
		{
			for (x = 1; x < width - 1; x++)
			{
				result.ATD(y, x)(ch) = sol[ch](index.at<int>(y, x));
			}
		}
	}
	// update mask
	mask(patch_mask).setTo(Scalar(M_DST));
	// get result
	Mat uresult;
	result.convertTo(uresult, CV_8UC3);
	uresult.copyTo(patch);
	// update dst
	result.copyTo(dst(patch_mask));
	return;
}