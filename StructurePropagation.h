#pragma  once
#include <Wang_Utility/OpenCvUtility.h>
#include <Wang_Utility/MatrixBasedOnVector.h>
class StructurePropagation
{
public:
	~StructurePropagation(){}
	void Run(const Mat1b &_mask,const Mat& _img,const vector<Point> &points,Mat& result);
	//利用目前类中已经存储的数据继续经行修补
	//void RunAgain();
	void SetParm(int _blocksize,int _samplestep,int _iscurve);
};