#include "stdafx.h"
#include "OpenCvUtility.h"
void LineInterpolation(Point points[2],vector<Point> &PointList)
{
	PointList.clear();
	int absx=std::abs(points[0].x-points[1].x);
	int absy=std::abs(points[0].y-points[1].y);
	int x_temp=points[0].x<points[1].x?1:-1;
	int y_temp=points[0].y<points[1].y?1:-1;
	if (absy>absx)
	{
		for (int itery=0;itery<=absy;itery++)
		{
			int y=points[0].y+y_temp*itery;
			int x=(double)(points[1].x-points[0].x)/(points[1].y-points[0].y)*y_temp*itery+points[0].x;
			PointList.push_back(Point(x,y));
		}
	}
	else
	{
		for (int iterx=0;iterx<=absx;iterx++)
		{
			int x=points[0].x+x_temp*iterx;
			int y=(double)(points[1].y-points[0].y)/(points[1].x-points[0].x)*x_temp*iterx+points[0].y;
			PointList.push_back(Point(x,y));
		}
	}
}
void Wang_GetCurve(const vector<Point>& mouse_points,vector<Point> &PointList)
{
	if(PointList.size()!=0)
	{
		PointList.clear();
		PointList.reserve(0);
	}
	vector<Point> temp_list;
	for (int i=0;i<mouse_points.size()-1;i++)
	{
		Point points_temp[2]={mouse_points[i],mouse_points[i+1]};
		LineInterpolation(points_temp,temp_list);
		PointList.insert(PointList.end(),temp_list.begin(),temp_list.end());
		//½«temp_listÖÃÎª0
		temp_list.clear();
		temp_list.reserve(0);
	}
}