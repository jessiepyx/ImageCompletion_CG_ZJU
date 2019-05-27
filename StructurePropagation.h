#pragma  once

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "OpenCvUtility.h"
#include "StructurePropagation.h"

using namespace std;
using namespace cv;

class StructurePropagation
{
public:
	~StructurePropagation(){}
	void Run(const Mat1b &_mask,const Mat& _img,const vector<Point> &points,Mat& result);
	//利用目前类中已经存储的数据继续经行修补
	//void RunAgain();
	void SetParm(int _blocksize,int _samplestep,int _iscurve);
};