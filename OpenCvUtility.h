#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void LineInterpolation(Point points[2],vector<Point> &PointList);
void GetCurve(const vector<Point>& mouse_points, vector<Point> &PointList);
void GetMask(const vector<Point>& points, const Mat& mat, Mat1b &mask);
extern Mat3b result_copy;

inline void DrawRect(Point p, int size) {
	rectangle(result_copy, Rect(p.x - size / 2, p.y - size / 2, size, size), CV_RGB(255, 0, 0), 2);
}

inline void DrawPoints(vector<Point> PointList,Mat &img,Scalar color,int r)
{
	for (int i=0;i<PointList.size();i++)
		circle(img,PointList[i],r,color,CV_FILLED);
}