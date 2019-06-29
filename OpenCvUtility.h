#ifndef OPENCV_UTILITY_H
#define OPENCV_UTILITY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void LineInterpolation(Point points[2], vector<Point> &PointList);
void GetCurve(const vector<Point>& mouse_points, vector<Point> &PointList);

inline void DrawPoints(vector<Point> PointList, Mat &img, Scalar color, int r)
{
	for (Point &p : PointList)
	{
		circle(img, p, r, color, CV_FILLED);
	}
}

inline Vec3b AlphaBlending(Vec3b pixel1, Vec3b pixel2, double alpha)
{
	Vec3b res;
	for (int i = 0; i < 3; i++)
	{
		res[i] = uchar(pixel1[i] * alpha + pixel2[i] * (1 - alpha));
	}
	return res;
}

#endif /* OPENCV_UTILITY_H */