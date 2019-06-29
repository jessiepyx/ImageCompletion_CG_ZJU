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

void GetCurve(const vector<Point>& mouse_points, vector<Point> &PointList)
{
	if(!PointList.empty())
	{
		PointList.clear();
	}
	for (int i=0;i<mouse_points.size()-1;i++)
	{
		Point points_temp[2]={mouse_points[i],mouse_points[i+1]};
		LineInterpolation(points_temp,PointList);
	}
}
