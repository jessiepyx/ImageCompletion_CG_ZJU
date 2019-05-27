#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "OpenCvUtility.h"
#include "StructurePropagation.h"

using namespace std;
using namespace cv;

Point points[2]={Point(-1,-1),Point(-1,-1)};
int points_i=0;
Mat3b img;
Mat1b mask;
Mat3b result;
Mat3b result_copy;
vector<Point> PointsList;
Point prev_pt(-1,-1);
vector<Point> mousepoints;
int blocksize=10;
int samplestep=2;
bool iscurve=false;
ofstream file;
void onmouse(int event,int x,int y,int flags,void* parm)
{
	if(!iscurve)
	{
		if (event!=CV_EVENT_LBUTTONDOWN)
			return;
		points[points_i].x=x;
		points[points_i].y=y;
		//cout<<x<<"  "<<y<<"  "<<points[0]<<" "<<points[1]<<endl;
		points_i=(points_i+1)%2;
		if (points[0].x!=-1&&points[1].x!=-1&&points_i==0)
		{
			result.copyTo(result_copy);
			file<<points[0]<<"  "<<points[1]<<endl;
			LineInterpolation(points,PointsList);
			DrawPoints(PointsList,result_copy,Scalar(255,0,255),1);//×ÏÉ«
			circle(result_copy,points[0],3,Scalar(255,0,0),CV_FILLED);//À¶É«
			circle(result_copy,points[1],3,Scalar(255,0,0),CV_FILLED);
//			rectangle(result_copy,RectByCenter(PointsList[0],blocksize),CV_RGB(255,0,0),2);
			imshow("img",result_copy);
		}
	}
	else
	{
		if( event == CV_EVENT_LBUTTONUP)// || !(flags & CV_EVENT_FLAG_LBUTTON) )
			prev_pt = cvPoint(-1,-1);
		if( event == CV_EVENT_LBUTTONDOWN )
		{
			prev_pt = cvPoint(x,y);
//			rectangle(result_copy,RectByCenter(prev_pt,blocksize),CV_RGB(255,0,0),2);
			mousepoints.push_back(prev_pt);
		}
		else if( event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON) )
		{
			CvPoint pt = cvPoint(x,y);
			mousepoints.push_back(pt);
			if( prev_pt.x < 0 )
				prev_pt = pt;
			//cvLine( inpaint_mask, prev_pt, pt, cvScalarAll(255), 5, 8, 0 );
			line( result_copy, prev_pt, pt, cvScalarAll(255), 1, 8, 0 );
			prev_pt = pt;
			imshow( "img", result_copy );
		}
	}
}

int main(int argc, char* argv[])
{
	img = imread("img.jpg", 1);
	mask = imread("mask.bmp", 0);

	//img=imread("curve_test1.png",1);
	//mask=imread("curve_test1.bmp",0);

	threshold(mask,mask,125,255,CV_THRESH_BINARY_INV);
	result.zeros(img.size());
	img.copyTo(result,mask);
	namedWindow("img");
//	namedWindow("mask");
	createTrackbar("BlockSize","img",&blocksize,50);
	createTrackbar("SampleStep","img",&samplestep,20);
	int iscurve_temp=0;
	createTrackbar("iscurve","img",&iscurve_temp, 1);
	setMouseCallback("img",onmouse);

	cout << result.size() << endl;
	cout << mask.size() << endl;

	imshow("img",result);
//	imshow("mask",mask);
	//StructurePropagation SP;
	//SP.SetParm(blocksize,samplestep,iscurve);
	Mat3b Local_Result_Copy(result.size());
	result.copyTo(Local_Result_Copy);
	result.copyTo(result_copy);
	//file.open("test.txt");
	for (;;)
	{
		iscurve=iscurve_temp;
		char c=cvWaitKey(10);
		if (c==27)
			break;
		else if (c=='s')
		{
			file<<blocksize<<"   "<<samplestep<<endl;
			if(iscurve)
//				Wang_GetCurve(mousepoints,PointsList);
			DrawPoints(PointsList,img,CV_RGB(255,0,0),1);
			//SP.SetParm(blocksize,samplestep,iscurve);
			//SP.Run(mask,img,PointsList,Local_Result_Copy);
//			imshow("img",Local_Result_Copy);
		}
		else if (c=='r')
		{
			result.copyTo(result_copy);
			result.copyTo(Local_Result_Copy);
			mousepoints.clear();
			imshow("img",result_copy);
		}
		else if (c=='a')
		{
			imwrite("result.jpg",Local_Result_Copy);
		}
		else if(c=='e')
		{
//			Wang_GetCurve(mousepoints,PointsList);
//			DrawPoints(PointsList,result_copy,CV_RGB(255,0,0),1);
//			imshow("img",result_copy);
		}
	}
	file.close();
	return 0;
}

