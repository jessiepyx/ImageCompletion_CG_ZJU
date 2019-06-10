#pragma  once

#include "OpenCvUtility.h"
#include "PointManager.h"
#include "OpenCvUtility.h"
#include <map>
#include <list>
#include <set>
#include <math.h>

using namespace std;
using namespace cv;

class StructurePropagation
{
public:
	StructurePropagation() {
		ki = 0.75;
		ks = 0.25;
		//ki = 0.25;
		//ks = 0.75;
	}
	~StructurePropagation(){}
	void SetParam(int block_size,int sample_step,int line_or_curve);
	void Run(const Mat &mask, const Mat& img, Mat &smask, vector<vector<Point>> &plist, Mat& result);

//	void TextureCompletion(const Mat1b &_mask, Mat1b &LineMask, const Mat &mat, Mat &result);
//	void TextureCompletion2(Mat1b _mask, Mat1b LineMask, const Mat &mat, Mat &result);

private:
	int block_size;
	int sample_step;
	int line_or_curve;
	double ki;
	double ks;
	PointManager pointManager;

	void getResult(Mat1b mask, int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result);
	void ModifyMask(Mat &LineMask, vector<PointPos>AnchorPoints);
	int *DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	int *BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	double gauss(double x);
	double calcEs(const PointPos &i, const PointPos &xi);
	double calcEi(const Mat &mat, const PointPos &i, const PointPos &xi);
	double calcE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2);
	void calcMij(MyNode &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints);
};