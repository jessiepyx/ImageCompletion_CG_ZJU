#ifndef PHOTOMETRIC_H
#define PHOTOMETRIC_H

#include <opencv2/opencv.hpp>

using namespace cv;

class Photometric
{
public:
	Photometric() = default;
	~Photometric() = default;
	static void initMask(Mat image, Mat imageMask, uchar unknown = 0, uchar known = 255);
	static void correctE(Mat &patch, int offset_x, int offset_y);

	// M_BORDER for border, M_DST for base image, M_SRC for src image
	static Mat mask;
	// original image
	static Mat dst;
};

#endif /* PHOTOMETRIC_H */
