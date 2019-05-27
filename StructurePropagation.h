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
	//����Ŀǰ�����Ѿ��洢�����ݼ��������޲�
	//void RunAgain();
	void SetParm(int _blocksize,int _samplestep,int _iscurve);
};