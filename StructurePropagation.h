#pragma  once
#include <Wang_Utility/OpenCvUtility.h>
#include <Wang_Utility/MatrixBasedOnVector.h>
class StructurePropagation
{
public:
	~StructurePropagation(){}
	void Run(const Mat1b &_mask,const Mat& _img,const vector<Point> &points,Mat& result);
	//����Ŀǰ�����Ѿ��洢�����ݼ��������޲�
	//void RunAgain();
	void SetParm(int _blocksize,int _samplestep,int _iscurve);
};