#pragma once

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void LineInterpolation(Point points[2],vector<Point> &PointList);
inline void DrawPoints(vector<Point> PointList,Mat &img,Scalar color,int r)
{
	for (int i=0;i<PointList.size();i++)
		circle(img,PointList[i],r,color,CV_FILLED);
}