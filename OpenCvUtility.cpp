#include "OpenCvUtility.h"

void LineInterpolation(Point points[2],vector<Point> &PointList)
{
	int absx=std::abs(points[0].x-points[1].x);
	int absy=std::abs(points[0].y-points[1].y);
	int x_temp=points[0].x<points[1].x?1:-1;
	int y_temp=points[0].y<points[1].y?1:-1;
	if (PointList.size() == 0) {
		PointList.push_back(points[0]);
	}
	if (points[0].x == points[1].x && points[0].y == points[1].y) {
		return;
	}
	if (absy>absx)
	{
		for (int itery=1;itery<=absy;itery++)
		{
			int y=points[0].y+y_temp*itery;
			int x=(double)(points[1].x-points[0].x)/(points[1].y-points[0].y)*y_temp*itery+points[0].x;
			PointList.push_back(Point(x,y));
		}
	}
	else
	{
		for (int iterx=1;iterx<=absx;iterx++)
		{
			int x=points[0].x+x_temp*iterx;
			int y=(double)(points[1].y-points[0].y)/(points[1].x-points[0].x)*x_temp*iterx+points[0].y;
			PointList.push_back(Point(x,y));
		}
	}
}

void Wang_GetCurve(const vector<Point>& mouse_points, vector<Point> &PointList)
{
	if(PointList.size()!=0)
	{
		PointList.clear();
	}
	for (int i=0;i<mouse_points.size()-1;i++)
	{
		Point points_temp[2]={mouse_points[i],mouse_points[i+1]};
		LineInterpolation(points_temp,PointList);
	}
}

void getMask(const vector<Point>& mouse_points, const Mat& img, Mat1b &mask) {
	vector<Point> mask_points;
	mask = Mat::zeros(img.rows, img.cols, CV_8UC1);
	Wang_GetCurve(mouse_points, mask_points);
	Point points_temp[2] = { mouse_points[mouse_points.size() - 1], mouse_points[0]};
	LineInterpolation(points_temp, mask_points);
	int top = img.rows, bottom = 0;
	int maxx[1010], minx[1010];
	for (int i = 0; i < img.rows; i++)
	{
		maxx[i] = 0;
		minx[i] = img.cols;
	}
	for (int i = 0; i < mask_points.size(); i++) {
		int y = mask_points[i].y;
		int x = mask_points[i].x;
		maxx[y] = max(maxx[y], x);
		minx[y] = min(minx[y], x);
	}
	for (int i = 0; i < mask.rows;i++)
	for (int j = 0; j < mask.cols; j++)
	{
		mask.at<uchar>(i, j) = !(j>=minx[i] && j<=maxx[i]) * 255;
	}
	/*
	for (int i = 0; i < mask_points.size(); i++) {
		int y = mask_points[i].y;
		int x = mask_points[i].x;
		if (y < top) {
			top = y;
		}
		if (y > bottom) {
			bottom = y;
		}
		mask.at<uchar>(y, x) = 255;
	}
	bool inMask;
	for (int i = 0; i < mask.rows; i++) {
		inMask = false;
		const uchar* ptr = mask.ptr<uchar>(i);
		for (int j = 0; j < mask.cols; j++) {
			if (ptr[j] == 255) {
				while (ptr[++j] == 255);
				inMask = !inMask;
			}
			if (i <= top || i >= bottom || inMask == false) {
				mask.at<uchar>(i, j) = 255;
			}
		}
	}*/
}