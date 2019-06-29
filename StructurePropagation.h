#ifndef STRUCTURE_PROPAGATION_H
#define STRUCTURE_PROPAGATION_H

#include "OpenCvUtility.h"
#include "PointManager.h"
#include <map>
#include <list>
#include <set>
#include <math.h>

using namespace std;
using namespace cv;

class StructurePropagation
{
public:
	StructurePropagation() = default;
	~StructurePropagation() = default;
	void SetParam(int block_size, int sample_step, int line_or_curve, double ks, double ki);
	void Run(const Mat &mask, const Mat &img, Mat &mask_structure, vector<vector<Point>> &plist, Mat &result);

private:
	int block_size;
	int sample_step;
	int line_or_curve;
	double ks;
	double ki;
	PointManager pointManager;

	int *DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	double computeEs(const PointPos &i, const PointPos &xi);
	double computeEi(const Mat &mat, const PointPos &i, const PointPos &xi);
	double computeE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2);
	int *BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	void computeMij(MyNode &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints);
	void getResult(Mat mask, int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result);
};

#endif /* STRUCTURE_PROPAGATION_H */