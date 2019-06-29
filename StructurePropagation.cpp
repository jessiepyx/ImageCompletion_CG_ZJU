#include "StructurePropagation.h"
#include <algorithm>
#include <queue>

#include "Photometric.h"

/**
 * Set parameters.
 * @param block_size 		side length of square cpatch
 * @param sample_step 		distance between neighboring anchor points
 * @param line_or_curve 	structure is straight line (0) or curve (1)
 * @param ks 				coefficient of Es
 * @param ki 				coefficient of Ei
 */
void StructurePropagation::SetParam(int block_size, int sample_step, int line_or_curve, double ks, double ki)
{
	this->block_size = block_size;
	this->sample_step = sample_step;
	this->line_or_curve = line_or_curve;
	this->ks = ks;
	this->ki = ki;
}

/**
 * Run the structural propagation algorithm.
 * @param mask 				0 for known region, 255 for known region
 * @param img 				0 for unknown region
 * @param mask_structure 	mask after structure propagation
 * @param plist 			list of sample points along structure line/curve
 * @param result 			img after structure propagation
 */
void StructurePropagation::Run(const Mat &mask, const Mat& img, Mat &mask_structure, vector<vector<Point>> &plist, Mat& result)
{
	Mat grayMat = Mat::zeros(img.rows, img.cols, CV_8UC1);
	cvtColor(img, grayMat, CV_BGR2GRAY);

	set<shared_ptr<list<int>>> lineSets;
	pointManager.reset(plist, grayMat, block_size, lineSets);

	int *sampleIndices;
	vector<PointPos> anchorPoints;
	vector<PointPos> samplePoints;
	set<shared_ptr<list<int>>>::iterator itor;
	for (itor = lineSets.begin(); itor != lineSets.end(); itor++)
	{
		pointManager.getSamplePoints(samplePoints, sample_step, **itor);
		if (!samplePoints.empty())
		{
			// compute patch matches
			if ((*itor)->size() == 1)
			{
				pointManager.getAnchorPoints(anchorPoints, **itor);
				sampleIndices = DP(samplePoints, anchorPoints, grayMat);
			}
			else
			{
				pointManager.constructBPMap(**itor);
				sampleIndices = BP(samplePoints, anchorPoints, grayMat);
			}

			// update mask (mark anchored patches as known)
			for (auto p : anchorPoints)
			{
				Point tar = pointManager.getPoint(p);
				for (int j = -block_size / 2; j < block_size / 2; j++)
				{
					for (int k = -block_size / 2; k < block_size / 2; k++)
					{
						int y = j + tar.y;
						int x = k + tar.x;
						if (x >= 0 && y >= 0 && x < mask_structure.cols && y < mask_structure.rows)
						{
							mask_structure.at<uchar>(y, x) = 255;
						}
					}
				}
			}

			getResult(mask, sampleIndices, samplePoints, anchorPoints, result);
		}
	}
}

/**
 * Dynamic programing for single-line structure.
 * @param samplePoints 		known patch centers
 * @param anchorPoints 		unknown patch centers
 * @param mat 				input gray scale image
 * @return					indices
 */
int *StructurePropagation::DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat)
{

	auto *M = (double *)malloc(2 * samplePoints.size() * sizeof(double));
	auto *record = (int *)malloc(samplePoints.size() * anchorPoints.size() * sizeof(int));

	// first anchor point
	for (int xi = 0; xi < samplePoints.size(); xi++)
	{
		M[xi] = ks * computeEs(anchorPoints[0], samplePoints[xi]) + ki * computeEi(mat, anchorPoints[0], samplePoints[xi]);
	}

	// for each anchor point, for each xi, compute M
	int curOffset = 0;
	int preOffset = 0;
	for (int i = 1; i < anchorPoints.size(); i++)
	{
		curOffset = (i % 2) * samplePoints.size();
		preOffset = ((i + 1) % 2) * samplePoints.size();

		for (int xi = 0; xi < samplePoints.size(); xi++)
		{
			// compute E1
			double E1 = ks * computeEs(anchorPoints[i], samplePoints[xi]) + ki * computeEi(mat, anchorPoints[i], samplePoints[xi]);

			// compute E2 and
			double min = INT_MAX;
			int minIndex = 0;
			for (int xj = 0; xj < samplePoints.size(); xj++)
			{
				double tmp = computeE2(mat, anchorPoints[i], anchorPoints[i - 1], samplePoints[xi], samplePoints[xj]) + M[preOffset + xj];
				if (tmp < min)
				{
					min = tmp;
					minIndex = xj;
				}
			}

			// record xi and M
			record[samplePoints.size() * i + xi] = minIndex;
			M[curOffset + xi] = E1 + min;
		}
	}

	// last anchor point
	int *sampleIndices = (int*)malloc(anchorPoints.size() * sizeof(int));
	double min = INT_MAX;
	for (int xi = 0; xi < samplePoints.size(); xi++)
	{
		if (M[curOffset + xi] < min)
		{
			sampleIndices[anchorPoints.size() - 1] = xi;
			min = M[curOffset + xi];
		}
	}

	// trace back
	for (int i = anchorPoints.size() - 2; i >= 0; i--)
	{
		sampleIndices[i] = record[samplePoints.size() * (i + 1) + sampleIndices[i + 1]];
	}

	free(M);
	free(record);

	return sampleIndices;
}

/**
 * Term for structural similarity.
 * @param i 	which anchor
 * @param xi 	which label
 * @return 		Es
 */
double StructurePropagation::computeEs(const PointPos &i, const PointPos &xi)
{
	// get points of curve segment contained in patch
	list<Point*> begin1, begin2;
	list<int> length1, length2;
	pointManager.getPointsinPatch(i, begin1, length1);
	pointManager.getPointsinPatch(xi, begin2, length2);

	int len1 = 0;
	for (auto l : length1)
	{
		len1 += l;
	}

	int len2 = 0;
	for (auto l : length2)
	{
		len2 += l;
	}
	static vector<int> minDist1(len1), minDist2(len2);

	// initialize minimal distance
	for (int i = 0; i < len1; i++)
	{
		minDist1[i] = INT_MAX;
	}
	for (int i = 0; i < len2; i++)
	{
		minDist2[i] = INT_MAX;
	}

	Point pi = pointManager.getPoint(i);
	Point pxi = pointManager.getPoint(xi);
	int offsetx = pxi.x - pi.x;
	int offsety = pxi.y - pi.y;
	list<int>::iterator lenItor1, lenItor2;
	list<Point*>::iterator pointItor1, pointItor2;

	// compute minimal distance
	for (lenItor1 = length1.begin(), pointItor1 = begin1.begin(); lenItor1 != length1.end(); lenItor1++, pointItor1++)
	{
		Point *points1 = *pointItor1;
		for (int i = 0; i < *lenItor1; i++)
		{
			for (lenItor2 = length2.begin(), pointItor2 = begin2.begin(); lenItor2 != length2.end(); lenItor2++, pointItor2++)
			{
				for (int j = 0; j < *lenItor2; j++)
				{
					Point *points2 = *pointItor2;
					int dx = points1[i].x - points2[j].x + offsetx;
					int dy = points1[i].y - points2[j].y + offsety;
					int dist = dx * dx + dy * dy;
					if (dist < minDist1[i])
					{
						minDist1[i] = dist;
					}
					if (dist < minDist2[j])
					{
						minDist2[j] = dist;
					}
				}
			}
		}
	}

	int Es = 0;
	for (auto d : minDist1)
	{
		Es += d;
	}
	for (auto d : minDist2)
	{
		Es += d;
	}
	return (double)Es / minDist1.size();
}

/**
 * Term for completion.
 * @param mat 	input image
 * @param i 	which anchor
 * @param xi 	which label
 * @return 		Ei
 */
double StructurePropagation::computeEi(const Mat &mat, const PointPos &i, const PointPos &xi)
{
	// compute Ei for every boundary patch
	if (pointManager.nearBoundary(i))
	{
		Point pi = pointManager.getPoint(i);
		Point pxi = pointManager.getPoint(xi);
		int offset1 = block_size / 2;
		int offset2 = block_size - offset1;

		int cnt = 0;
		int ssd = 0;
		for (int i = -offset1; i < offset2; i++)
		{
			const auto *ptri = mat.ptr<uchar>(i + pi.y);
			const auto *ptrxi = mat.ptr<uchar>(i + pxi.y);
			for (int j = -offset1; j < offset2; j++)
			{
				if (ptri[j + pi.x] != 0)
				{
					int diff = ptri[j + pi.x] - ptrxi[j + pxi.x];
					ssd += diff * diff;
					cnt++;
				}
			}
		}
		return (double)ssd / cnt;
	}
	else
	{
		return 0;
	}
}

/**
 * Term for coherence.
 * @param mat 	input image
 * @param i1 	which anchor
 * @param i2 	which anchor
 * @param xi1 	which label
 * @param xi2 	which label
 * @return 		E2
 */
double StructurePropagation::computeE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2)
{
	Point p1 = pointManager.getPoint(i1);
	Point p2 = pointManager.getPoint(i2);
	Point px1 = pointManager.getPoint(xi1);
	Point px2 = pointManager.getPoint(xi2);

	int left1, left2, right1, right2;
	int up1, up2, down1, down2;
	if (p1.x > p2.x)
	{
		left1 = 0;
		left2 = p1.x - p2.x;
		right1 = block_size - left2;
		right2 = block_size;
	}
	else
	{
		left2 = 0;
		left1 = p2.x - p1.x;
		right2 = block_size - left1;
		right1 = block_size;
	}

	if (p1.y > p2.y)
	{
		up1 = 0;
		up2 = p1.y - p2.y;
		down1 = block_size - up2;
		down2 = block_size;
	}
	else
	{
		up2 = 0;
		up1 = p2.y - p1.y;
		down2 = block_size - up1;
		down1 = block_size;
	}

	// compute E2 between every pair of neighboring patches
	if (right1 >= 0 && right2 >= 0 && down1 >= 0 && down2 >= 0)
	{
		int cols = right1 - left1;
		int rows = down1 - up1;

		double ssd = 0;
		for (int i = 0; i < rows; i++)
		{
			const auto *ptr1 = mat.ptr<uchar>(i + up1 + px1.y - block_size / 2);
			const auto *ptr2 = mat.ptr<uchar>(i + up2 + px2.y - block_size / 2);
			for (int j = 0; j < cols; j++)
			{
				double diff = ptr1[j + left1 + px1.x - block_size / 2] - ptr2[j + left2 + px2.x - block_size / 2];
				ssd += diff * diff;
			}
		}
		return ssd / (cols * rows);
	}
	else
	{
		return 0;
	}
}

/**
 * Belief propagation for multi-line structure.
 * @param samplePoints 		known patch centers
 * @param anchorPoints 		unknown patch centers
 * @param mat 				input gray scale image
 * @return					indices
 */
int *StructurePropagation::BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat)
{
	//initialization
	int size = pointManager.getPropstackSize();
	anchorPoints.clear();
	anchorPoints.reserve(size);

	list<shared_ptr<MyNode>>::iterator itor;
	list<shared_ptr<MyNode>>::iterator end;
	pointManager.getPropstackItor(itor, end);

	// receive messages from neighbors
	for (; itor != end; itor++)
	{
		shared_ptr<MyNode> n = *itor;
		list<shared_ptr<Edge>> edges;
		n->getEdges(edges);
		// message only for next neighbor (the node that enqueued this node)
		computeMij(*n, n->getEdgeBegin(), mat, samplePoints);
	}

	auto *sampleIndices = (int*)malloc(size * sizeof(int));
	auto *cur = (double*)malloc(samplePoints.size() * sizeof(double));

	list<shared_ptr<MyNode>>::reverse_iterator rev_itor;
	list<shared_ptr<MyNode>>::reverse_iterator rev_end;
	pointManager.getPropstackReverseItor(rev_itor, rev_end);

	// send updated messages back to neighbors
	for (int i = 0; rev_itor != rev_end; rev_itor++, i++)
	{
		shared_ptr<MyNode> n = *rev_itor;
		list<shared_ptr<Edge>>::iterator begin = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		list<shared_ptr<Edge>>::iterator itor = begin;
		anchorPoints.push_back(n->p);

		// messages for all neighbors
		for (itor++; itor != end; itor++)
		{
			computeMij(*n, itor, mat, samplePoints);
		}

		// compute E1 for all possible xi
		int minIndex = 0;
		double min = INT64_MAX;
		for (int xi = 0; xi < samplePoints.size(); xi++)
		{
			cur[xi] = ks * computeEs(n->p, samplePoints[xi]) + ki * computeEi(mat, n->p, samplePoints[xi]);
		}

		// add up all messages sent to this node
		for (itor = begin; itor != end; itor++)
		{
			double **toMptr = (*itor)->getMbyTo(n->id);
			for (int i = 0; i < samplePoints.size(); i++)
			{
				cur[i] += (*toMptr)[i];
			}
		}

		// find the optimal xi
		for (int i = 0; i < samplePoints.size(); i++)
		{
			if (cur[i] < min)
			{
				min = cur[i];
				minIndex = i;
			}
		}
		sampleIndices[i] = minIndex;
	}

	// release resources
	pointManager.getPropstackItor(itor, end);
	for (; itor != end; itor++)
	{
		shared_ptr<MyNode> n = *itor;
		list<shared_ptr<Edge>>::iterator edgeItor = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		for (; edgeItor != end; edgeItor++)
		{
			double **M = (*edgeItor)->getMbyFrom(n->id);
			if (*M != nullptr)
			{
				free(*M);
			}
		}
	}
	free(cur);

	return sampleIndices;
}

/**
 * Messages in belief propagation.
 * @param n					sender struct
 * @param edgeItor 			edge to receiver
 * @param mat 				input image
 * @param samplePoints 		sample points
 */
void StructurePropagation::computeMij(MyNode &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints)
{
	double **Mptr = (*edgeItor)->getMbyFrom(n.id);
	auto end = n.getEdgeEnd();

	if (*Mptr == nullptr)
	{
		*Mptr = (double*)malloc(samplePoints.size() * sizeof(double));
		memset(*Mptr, 0, samplePoints.size() * sizeof(double));

		for (int i = 0; i < samplePoints.size(); i++)
		{
			double E1 = ks * computeEs(n.p, samplePoints[i]) + ki * computeEi(mat, n.p, samplePoints[i]);

			// add up messages sent from Mki (k != j)
			double msg = 0;
			for (auto itor = n.getEdgeBegin(); itor != end; itor++)
			{
				if (itor != edgeItor)
				{
					double **toMptr = (*itor)->getMbyTo(n.id);
					if (*toMptr == nullptr)
					{
						assert(0);
					}
					msg += (*toMptr)[i];
				}
			}

			PointPos tmpPos = pointManager.getPointPos((*edgeItor)->getAnother(n.id));
			for (int j = 0; j < samplePoints.size(); j++)
			{
				// update each item in Mij
				double E2 = computeE2(mat, n.p, tmpPos, samplePoints[i], samplePoints[j]);
				if ((*Mptr)[j] == 0 || E1 + E2 + msg < (*Mptr)[j])
				{
					(*Mptr)[j] = E1 + E2 + msg;
				}
			}
		}
	}
}

void StructurePropagation::getResult(Mat mask, int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat &result)
{
	vector<vector<int> > tmp_mask(mask.rows, vector<int>(mask.cols, 0));
	for (int i = 0; i < mask.rows; i++)
	{
		for (int j = 0; j < mask.cols; j++)
		{
			tmp_mask[i][j] = (mask.at<uchar>(i, j) > 0);
			if (tmp_mask[i][j])
			{
				mask.at<uchar>(i, j) = 255;
			}
		}
	}

	Photometric::initMask(result, mask);

	int offset1 = block_size / 2;
	int offset2 = block_size - offset1;

	// copy all sample patches to corresponding anchor pathces
	for (int i = 0; i < anchorPoints.size(); i++)
	{
		Point src = pointManager.getPoint(samplePoints[sampleIndices[i]]);
		Point tar = pointManager.getPoint(anchorPoints[i]);

		Mat patch = result(Rect(src.x - offset1, src.y - offset1, block_size, block_size)).clone();
		Photometric::correctE(patch, src.x - offset1, src.y - offset1);

		for (int m = -offset1; m < offset2; m++)
		{
			int tary = tar.y + m;
			const Vec3b* srcPtr = result.ptr<Vec3b>(src.y + m);
			for (int n = -offset1; n < offset2; n++)
			{
				Vec3b tmp = result.at<Vec3b>(tar.y + m, tar.x + n);
				if (tmp_mask[tar.y + m][tar.x + n] == 0)
				{
					result.at<Vec3b>(tar.y + m, tar.x + n) = srcPtr[src.x + n];
					tmp_mask[tar.y + m][tar.x + n] = 1;
				}
				else
				{
					result.at<Vec3b>(tar.y + m, tar.x + n) = AlphaBlending(srcPtr[src.x + n], result.at<Vec3b>(tar.y + m, tar.x + n), 0.5);
//                    result.at<Vec3b>(tar.y + m, tar.x + n) = patch.at<Vec3b>(m + offset1, n + offset1);
//                    result.at<Vec3b>(tar.y + m, tar.x + n) = AlphaBlending(patch.at<Vec3b>(m + offset1, n + offset1), result.at<Vec3b>(tar.y + m, tar.x + n), 0.5);
				}
			}
		}
	}
	free(sampleIndices);
}