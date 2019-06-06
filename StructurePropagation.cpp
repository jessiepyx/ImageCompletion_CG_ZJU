#include "StructurePropagation.h"
#include <algorithm>
#include <queue>
#include <time.h>
#define mp make_pair

void StructurePropagation::SetParam(int block_size, int sample_step, int line_or_curve)
{
	this->block_size = block_size;
	this->sample_step = sample_step;
	this->line_or_curve = line_or_curve;
}

void StructurePropagation::Run(const Mat &mask, const Mat& img, Mat &smask, vector<vector<Point>> &plist, Mat& result)
{
	time_t start = time(NULL);
	Mat grayMat = Mat::zeros(img.rows, img.cols, CV_8UC1);
	cvtColor(img, grayMat, CV_BGR2GRAY);

	set<shared_ptr<list<int>>> lineSets;
	pointManager.reset(plist, grayMat, block_size, lineSets);

	int *sampleIndices;
	vector<PointPos> anchorPoints;
	vector<PointPos> samplePoints;
	set<shared_ptr<list<int>>>::iterator itor;
	for (itor = lineSets.begin(); itor != lineSets.end(); itor++) {
		pointManager.getSamplePoints(samplePoints, sample_step, **itor);
		if (samplePoints.size() == 0){
			continue;
		}
		if ((*itor)->size() == 1) {
			puts("DP");
			pointManager.getAnchorPoints(anchorPoints, **itor);
			sampleIndices = DP(samplePoints, anchorPoints, grayMat);
		}
		else {
			puts("BP");
			pointManager.constructBPMap(**itor);
			/*list<shared_ptr<Node>>::iterator begin, end;
			pointManager.getPropstackItor(begin, end);
			shared_ptr<Node> n = *begin;
			list<shared_ptr<Edge>> edges;
			list<shared_ptr<Node>> queue;
			vector<int> visit(Node::totalNum + 1);
			queue.push_back(n);
			while (queue.size()) {
				shared_ptr<Node> tmpn = *queue.begin();
				visit[tmpn->id] = 1;
				tmpn->getEdges(edges);
				list<shared_ptr<Edge>>::iterator itor;
				for (itor = edges.begin(); itor != edges.end(); itor++) {
					ushort nid = (*itor)->getAnother(tmpn->id);
					if (!visit[nid]) {
						queue.push_back(pointManager.getNode(nid));
						Point points[2];
						points[0] = pointManager.getPoint(pointManager.getPointPos(tmpn->id));
						points[1] = pointManager.getPoint(pointManager.getPointPos(nid));
						vector<Point> pointList;
						LineInterpolation(points, pointList);
						DrawPoints(pointList, result, CV_RGB(255, 0, 0), 1);
					}
				}
				queue.pop_front();
			}*/
			sampleIndices = BP(samplePoints, anchorPoints, grayMat);
		}
		time_t end = time(NULL);
		cout << "time consuming:" << end - start << endl;
		ModifyMask(smask, anchorPoints);
		getResult(mask, sampleIndices, samplePoints, anchorPoints, result);
	}
	/*pointManager.getSamplePoints(samplePoints, sampleStep, **lineSets.begin());
	// sampleIndices = DP(samplePoints, anchorPoints, grayMat);
	pointManager.constructBPMap();
	sampleIndices = BP(samplePoints, anchorPoints, grayMat);
	for (int i = 0; i < anchorPoints.size(); i++) {
		cout << sampleIndices[i] << ", ";
	}
	ModifyMask(Linemask, anchorPoints);
	getResult(_mask,sampleIndices, samplePoints, anchorPoints, result);*/
}

void StructurePropagation::ModifyMask(Mat &LineMask, vector<PointPos>AnchorPoints)
{
	int N = LineMask.rows, M = LineMask.cols;
	int offset1 = block_size / 2;
	int offset2 = block_size / 2;
	for (int i = 0; i < AnchorPoints.size(); i++)
	{
		Point tar = pointManager.getPoint(AnchorPoints[i]);
		for (int j = -offset1; j < offset2; j++)
		{
			for (int k = -offset1; k < offset2; k++)
			{
				int nx = tar.y + j;
				int ny = tar.x + k;
				if (nx < 0 || nx >= N || ny < 0 || ny >= M) continue;
				LineMask.at<uchar>(nx, ny) = 2;
			}
		}
	}
}

int **MyAlloc(int N, int M)
{
	int **MyNode;
	MyNode = (int**)malloc(sizeof(int*)* N);
	for (int i = 0; i < N; i++)
	{
		MyNode[i] = (int*)malloc(sizeof(int)* M);
	}
	return MyNode;
}

void MyFree(int **MyNode,int N,int M)
{
	for (int i = 0; i < N; i++)
	{
		free(MyNode[i]);
	}
	free(MyNode);
}

int sqr(int x)
{
	return x * x;
}

int dist(Vec3b V1, Vec3b V2)
{
	return sqr(int(V1[0]) - int(V2[0])) + sqr(int(V1[1]) - int(V2[1])) + sqr(int(V1[2]) - int(V2[2]));
	/*double pr = (V1[0] + V2[0]) * 0.5;
	return sqr(V1[0] - V2[0]) * (2 + (255 - pr) / 256)
	+ sqr(V1[1] - V2[1]) * 4
	+ sqr(V1[2] - V2[2]) * (2 + pr / 256);*/
}

void StructurePropagation::TextureCompletion2(Mat1b _mask, Mat1b LineMask, const Mat &mat, Mat &result)
{
	int N = _mask.rows;
	int M = _mask.cols;
	int knowncount = 0;
	for (int i = 0; i < N;i++)
	for (int j = 0; j < M; j++)
	{
		knowncount += (_mask.at<uchar>(i, j) == 255);
	}
	if (knowncount * 2< N * M)
	{
		for (int i = 0; i < N;i++)
		for (int j = 0; j < M; j++)
			_mask.at<uchar>(i, j) = 255 - _mask.at<uchar>(i, j);
	}
	
	vector<vector<int> >my_mask(N, vector<int>(M, 0)), sum_diff(N, vector<int>(M, 0));
	
	for (int i = 0; i < N;i++)
	for (int j = 0; j < M; j++)
	LineMask.at<uchar>(i, j) = LineMask.at<uchar>(i, j) * 100;
	
	result = mat.clone();
	imshow("mask", _mask);
	imshow("linemask", LineMask);
	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		my_mask[i][j] = (_mask.at<uchar>(i, j) == 255);
		if (my_mask[i][j] == 0 && LineMask.at<uchar>(i, j) > 0)
		{
			my_mask[i][j] = 2;
		}
	}
	int bs = 5;
	int step = 1 * bs;
	auto usable(my_mask);
	int to_fill = 0, filled = 0;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		to_fill += (my_mask[i][j] == 0);
	}
	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		if (my_mask[i][j] == 1)
			continue;
		int k0 = max(0, i - step), k1 = min(N - 1, i + step);
		int l0 = max(0, j - step), l1 = min(M - 1, j + step);
		for (int k = k0; k <= k1; k++)
		for (int l = l0; l <= l1; l++)
			usable[k][l] = 2;
	}
	Mat use = _mask.clone();
	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	if (usable[i][j] == 2)
		use.at<uchar>(i, j) = 255;
	else use.at<uchar>(i, j) = 0;
	//imshow("usable", use);
	int itertime = 0;
	Mat match;
	while (true)
	{
		itertime++;
		int x, y, cnt = -1;
		for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			if (my_mask[i][j] != 0) continue;
			bool edge = false;
			int k0 = max(0, i - 1), k1 = min(N - 1, i + 1);
			int l0 = max(0, j - 1), l1 = min(M - 1, j + 1);
			for (int k = k0; k <= k1;k++)
			for (int l = l0; l <= l1; l++)
				edge |= (my_mask[k][l] == 1);
			if (!edge) continue;
			k0 = max(0, i - bs), k1 = min(N - 1, i + bs);
			l0 = max(0, j - bs), l1 = min(M - 1, j + bs);
			int tmpcnt = 0;
			for (int k = k0; k <= k1; k++)
			for (int l = l0; l <= l1; l++)
				tmpcnt += (my_mask[k][l] == 1);
			if (tmpcnt > cnt)
			{
				cnt = tmpcnt;
				x = i;
				y = j;
			}
		}
		if (cnt == -1) break;

		bool debug = false;
		bool debug2 = false;
		int k0 = min(x, bs), k1 = min(N - 1 - x, bs);
		int l0 = min(y, bs), l1 = min(M - 1 - y, bs);
		int sx, sy, min_diff = INT_MAX;
		for (int i = step; i + step < N; i += step)
		for (int j = step; j + step < M; j += step)
		{
			if (usable[i][j] == 2)continue;
			int tmp_diff = 0;
			for (int k = -k0; k <= k1; k++)
			for (int l = -l0; l <= l1; l++)
			{
				//printf("%d %d %d %d %d %d\n", i + k, j + l, x + k, y + l, N, M);
				if (my_mask[x + k][y + l] != 0)
					tmp_diff += dist(result.at<Vec3b>(i + k, j + l), result.at<Vec3b>(x + k, y + l));
			}
			sum_diff[i][j] = tmp_diff;
			if (min_diff > tmp_diff)
			{
				sx = i;
				sy = j;
				min_diff = tmp_diff;
			}
		}


		if (debug)
		{
			printf("x = %d y = %d\n", x, y);
			printf("sx = %d sy = %d\n", sx, sy);
			printf("mindiff = %d\n", min_diff);
		}
		if (debug2)
		{
			match = result.clone();
		}
		for (int k = -k0; k <= k1; k++)
		for (int l = -l0; l <= l1; l++)
		if (my_mask[x + k][y + l] == 0)
		{
			result.at<Vec3b>(x + k, y + l) = result.at<Vec3b>(sx + k, sy + l);
			my_mask[x + k][y + l] = 1;
			filled++;
			if (debug)
			{
				result.at<Vec3b>(x + k, y + l) = Vec3b(255, 0, 0);
				result.at<Vec3b>(sx + k, sy + l) = Vec3b(0, 255, 0);
			}
			if (debug2)
			{
				match.at<Vec3b>(x + k, y + l) = Vec3b(255, 0, 0);
				match.at<Vec3b>(sx + k, sy + l) = Vec3b(0, 255, 0);
			}
		}
		else
		{
			if (debug)
			{
				printf("(%d,%d,%d) matches (%d,%d,%d)\n", result.at<Vec3b>(x + k, y + l)[0], result.at<Vec3b>(x + k, y + l)[1], result.at<Vec3b>(x + k, y + l)[2], result.at<Vec3b>(sx + k, sy + l)[0], result.at<Vec3b>(sx + k, sy + l)[1], result.at<Vec3b>(sx + k, sy + l)[2]);
			}
		}
		if (debug2)
		{
			imshow("match", match);
		}
		if (debug) return;
		printf("done :%.2lf%%\n", 100.0 * filled / to_fill);
	}
}


void StructurePropagation::TextureCompletion(const Mat1b &_mask, Mat1b &LineMask, const Mat &mat, Mat &result)
{
	const int PatchSize = 5;
	const int step = 2*PatchSize;
	
	int N = mat.rows;
	int M = mat.cols;

	//prework
	int **my_mask, **neighbor;
	int **pixelcnt;
	my_mask = MyAlloc(N, M);
	pixelcnt = MyAlloc(N, M);//count of known pixels around (i,j)
	neighbor = MyAlloc(N, M);
	/*
	my_mask = (int**)malloc(sizeof(int*) * N);
	for (int i = 0; i < N; i++)
	{
	my_mask[i] = (int*)malloc(sizeof(int)* M);
	}*/

	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		int tmp = _mask.at<uchar>(i, j);
		if (tmp == 255) tmp = 1;
		my_mask[i][j] = tmp;
	}

	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		pixelcnt[i][j] = 0;
		neighbor[i][j] = 0;
		for (int k = -PatchSize; k <= PatchSize; k++)
		for (int l = -PatchSize; l <= PatchSize; l++)
		{
			if (i + k < 0 || i + k >= N) continue;
			if (j + l < 0 || j + l >= M) continue;
			pixelcnt[i][j] += (my_mask[i + k][j + l] == 1);
			neighbor[i][j]++;
		}
	}
	/*
	for (int i = 0; i < 20; i++)
	{
	for (int j = 0; j < 20; j++)printf("%d", my_mask[i][j]);
	puts("");
	}*/

	for (int i = 0; i < N;i++)
	for (int j = 0; j < M; j++)
	{
		if (LineMask.at<uchar>(i, j) == 2 && my_mask[i][j] == 0)
		{
			my_mask[i][j] = 2;
		}
	}

	//mask test
	
	/*
	for (int i = 0; i < N; i++)
	for (int j = 0; j < M; j++)
	{
		for (int k = 0; k < 3; k++)
			result.at<Vec3b>(i,j)[k] = 0;
		result.at<Vec3b>(i, j)[my_mask[i][j]] = 255;
	}
	return;*/

	// Get sample patches
	// Sample: [i-PatchSize,i+PatchSize] * [j-PatchSize,j+PatchSize]
	vector<pair<int, int> >Sample;
	for (int i = step; i < N - step; i += step)
	for (int j = step; j < M - step; j += step)
	{
		bool ok = true;
		for (int k = -PatchSize; k <= PatchSize; k++)
		for (int l = -PatchSize; l <= PatchSize; l++)
		{
			ok &= (my_mask[i + k][j + l] == 1);
		}
		if (ok) Sample.emplace_back(i, j);
	}

	printf("size = %d\n", Sample.size());

	//do it!
	priority_queue<pair<int,pair<int,int> > >heap;
	for (int i = 0; i < N; i++)
	{
		//printf("%d\n", i);
		for (int j = 0; j < M; j++)
		if (my_mask[i][j] == 0)
		{
			bool edge;
			edge = (i > 0 && my_mask[i - 1][j] == 1) ||
				(i < N - 1 && my_mask[i + 1][j] == 1) ||
				(j > 0 && my_mask[i][j - 1] == 1) ||
				(j < M - 1 && my_mask[i][j + 1] == 1);

			if (!edge) continue;
			heap.push(mp(pixelcnt[i][j], mp(i, j)));
		}
	}
	puts("done");

	result = mat.clone();

	
	int done = 0;
	while (!heap.empty())
	{
		auto tmp = heap.top();
		printf("ch = %d\n", tmp.first);
		heap.pop();
		int x = tmp.second.first;
		int y = tmp.second.second;
		if (my_mask[x][y] != 0) continue;
		//printf("solve %d %d\n", x, y);
		done++;
		if (done % 100 == 0)
		{
			printf("%d\n", done);
		}
		double minD = 1e9,patchnum = -1;
		int minSpaceD = 1e9;
		const double K1 = 1, K2 = 1;
		for (int i = 0; i < Sample.size(); i++)
		{
			int sx = Sample[i].first, sy = Sample[i].second;
			int D = 0;
			for (int j = -PatchSize; j <= PatchSize;j++)
			for (int k = -PatchSize; k <= PatchSize; k++)
			{
				if (x + j < 0 || x + j >= N) continue;
				if (y + k < 0 || y + k >= M) continue;
				if (my_mask[x + j][y + k] == 0) continue;
				D += dist(result.at<Vec3b>(x + j, y + k), result.at<Vec3b>(sx + j, sy + k));
			}
			if (D < minD)
			{
				minD = D;
			}
		}
		for (int i = 0; i < Sample.size(); i++)
		{
			int sx = Sample[i].first, sy = Sample[i].second;
			int D = 0;
			for (int j = -PatchSize; j <= PatchSize; j++)
			for (int k = -PatchSize; k <= PatchSize; k++)
			{
				if (x + j < 0 || x + j >= N) continue;
				if (y + k < 0 || y + k >= M) continue;
				if (my_mask[x + j][y + k] == 0) continue;
				D += dist(result.at<Vec3b>(x + j, y + k), result.at<Vec3b>(sx + j, sy + k));
			}
			if (D < minD * 1.2)
			{
				int SpaceD = sqr(sx - x) + sqr(sy - y);
				if (SpaceD < minSpaceD)
				{
					minSpaceD = SpaceD;
					patchnum = i;
				}
			}
		}
		int sx = Sample[patchnum].first, sy = Sample[patchnum].second;

		if (done == 30)
		{
			for (int i = -PatchSize; i <= PatchSize;i++)
			for (int j = -PatchSize; j <= PatchSize; j++)
			{
				result.at<Vec3b>(x + i, y + j)[0] = 255;
				result.at<Vec3b>(x + i, y + j)[1] = 0;
				result.at<Vec3b>(x + i, y + j)[2] = 0;

				result.at<Vec3b>(sx + i, sy + j)[0] = 0;
				result.at<Vec3b>(sx + i, sy + j)[1] = 255;
				result.at<Vec3b>(sx + i, sy + j)[2] = 0;
			}
			break;
		}
		//printf("target = %d %d\n", sx, sy);
		for (int j = -PatchSize; j <= PatchSize; j++)
		for (int k = -PatchSize; k <= PatchSize; k++)
		{
			if (x + j < 0 || x + j >= N) continue;
			if (y + k < 0 || y + k >= M) continue;
			if (my_mask[x + j][y + k] != 0) continue;
			my_mask[x + j][y + k] = 3;
			result.at<Vec3b>(x + j, y + k) = result.at<Vec3b>(sx + j, sy + k);
		}
		for (int j = -PatchSize-1; j <= PatchSize+1; j++)
		for (int k = -PatchSize-1; k <= PatchSize+1; k++)
		{
			if (x + j < 0 || x + j >= N) continue;
			if (y + k < 0 || y + k >= M) continue;
			if (my_mask[x + j][y + k] != 0) continue;
			int sx = x + j, sy = y + k;

			bool edge;
			edge = (sx > 0 && my_mask[sx - 1][sy] != 0) ||
				(sx < N - 1 && my_mask[sx + 1][sy] != 0) ||
				(sy > 0 && my_mask[sx][sy - 1] != 0) ||
				(sy < M - 1 && my_mask[sx][sy + 1] != 0);

			if (!edge) continue;

			pixelcnt[sx][sy] = 0;
			for (int l = -PatchSize; l <= PatchSize; l++)
			for (int m = -PatchSize; m <= PatchSize;m++)
			{
				if (sx + l < 0 || sx + l >= N) continue;
				if (sy + m < 0 || sy + m >= M) continue;
				pixelcnt[sx][sy] += (my_mask[sx + l][sy + m] % 2 != 0);
			}
			heap.push(mp(pixelcnt[sx][sy], mp(sx, sy)));
		}
	}

	//free
	MyFree(my_mask, N, M);
	MyFree(pixelcnt, N, M);
	MyFree(neighbor, N, M);
}

double StructurePropagation::gauss(double x)
{
	double sigma = block_size / 2;
	double pi = acos(-1);
	return exp(-x / (2 * sigma*sigma)) / (sqrt(2*pi) * sigma);
}

Vec3b fuse(Vec3b v1, Vec3b v2,double weight)
{
	// printf("%.5lf\n", weight);
	Vec3b res;
	for (int i = 0; i < 3; i++)
	{
		res[i] = int(v1[i] * weight + v2[i] * (1 - weight));
	}
	return res;
}

uchar adjust(uchar v0, double alpha, uchar v1, uchar v2)
{
	int V0 = v0, V1 = v1, V2 = v2;
	double result = V0 + alpha * (v2 - v1);
	if (result > 255 || result < 1) return V0;
	return uchar(result);
}


void StructurePropagation::getResult(Mat1b mask, int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result)
{
	// copy all sample patches to corresponding anchor pathces
	for (int i = 0; i < anchorPoints.size(); i++) {
		cout << sampleIndices[i] << ",";
	}
	cout << "sample number: " << samplePoints.size() << endl;
	cout << "anchor number: " << anchorPoints.size() << endl;
	int offset1 = block_size / 2;
	int offset2 = block_size - offset1;
	int N = mask.rows, M = mask.cols;
	vector<vector<int> > my_mask(N, vector<int>(M, 0));
	for (int i = 0; i < N;i++)
	for (int j = 0; j < M; j++)
	{
		my_mask[i][j] = (mask.at<uchar>(i, j) > 0);
		if (my_mask[i][j]) mask.at<uchar>(i, j) = 255;
	}
	//imshow("mask", mask);
	for (int i = 0; i < anchorPoints.size(); i++) {
		Point src = pointManager.getPoint(samplePoints[sampleIndices[i]]);
		Point tar = pointManager.getPoint(anchorPoints[i]);
		
		/*
		for (int m = -offset1; m < offset2; m++) {
			int tary = tar.y + m;
			const Vec3b* srcPtr = result.ptr<Vec3b>(src.y + m);
			
			if (my_mask[tar.y + m][tar.x - offset1] == 0)
			{
				result.at<Vec3b>(tar.y + m, tar.x - offset1) = srcPtr[src.x - offset1];
				my_mask[tar.y + m][tar.x - offset1] = 1;
			}
			if (my_mask[tar.y + m][tar.x + offset2 - 1] == 0)
			{
				result.at<Vec3b>(tar.y + m, tar.x + offset2 - 1) = srcPtr[src.x +offset2 - 1];
				my_mask[tar.y + m][tar.x + offset2 - 1] = 1;
			}
			int n = -offset1;
			while (n < offset2)
			{
				if (my_mask[tar.y + m][tar.x + n] == 1)
				{
					n++;
					continue;
				}
				//printf("%d %d\n", tar.y + m, tar.x + n);
				int nxt = n;
				while (my_mask[tar.y + m][tar.x + nxt] == 0)
				{
					nxt++;
					//printf("%d %d\n", nxt, offset2);
				}
				for (int ch = 0; ch < 3; ch++)
				{
					int srcdelta = int(result.at<Vec3b>(src.y + m, src.x + nxt)[ch] - int(result.at<Vec3b>(src.y + m, src.x + n - 1)[ch]));
					int tardelta = int(result.at<Vec3b>(tar.y + m, tar.x + nxt)[ch] - int(result.at<Vec3b>(tar.y + m, tar.x + n - 1)[ch]));
					if (srcdelta == 0) srcdelta = 1;
					double alpha = 1.0 * tardelta / srcdelta;
					for (int i = n; i < nxt; i++)
					{
						result.at<Vec3b>(tar.y + m, tar.x + i)[ch] = adjust(result.at<Vec3b>(tar.y + m, tar.x + i - 1)[ch], alpha, result.at<Vec3b>(src.y + m, src.x + i - 1)[ch], result.at<Vec3b>(src.y + m, src.x + i)[ch]);
					}
				}
				n = nxt;
			}
			for (int n = -offset1; n < offset2; n++)
				my_mask[tar.y + m][tar.x + n] = 1;

		}*/
		
		for (int m = -offset1; m < offset2; m++) {
			int tary = tar.y + m;
			const Vec3b* srcPtr = result.ptr<Vec3b>(src.y + m);
			for (int n = -offset1; n < offset2; n++) {
				Vec3b tmp = result.at<Vec3b>(tar.y + m, tar.x + n);
				if (my_mask[tar.y + m][tar.x + n] == 0) {
					result.at<Vec3b>(tar.y + m, tar.x + n) = srcPtr[src.x + n];
					my_mask[tar.y + m][tar.x + n] = 1;
				}
				else {
					result.at<Vec3b>(tar.y + m, tar.x + n) = fuse(srcPtr[src.x + n], result.at<Vec3b>(tar.y + m, tar.x + n), 0.5);
				}
			}
		}

		/*
		for (int m = -offset1; m < offset2; m++) {
			int tary = tar.y + m;
			const Vec3b* srcPtr = result.ptr<Vec3b>(src.y + m);
			for (int n = -offset1; n < offset2; n++) {
				Vec3b tmp = result.at<Vec3b>(tar.y + m, tar.x + n);
				if (tmp[0] == 0 && tmp[1] == 0 && tmp[2] == 0) {
					result.at<Vec3b>(tar.y + m, tar.x + n) = srcPtr[src.x + n];
				}
			}
		}*/
		

		//brightness fix
		/*for (int m = -offset1; m < offset2; m++)
		for (int n=-offset1;n<offset2;n++)
		{
		if (mask.at<uchar>(tar.y + m, tar.x + n) == 255)
		{
		// puts("fuse");
		// printf("%d %d\n", tar.y + m, tar.x + n);
		result.at<Vec3b>(tar.y + m, tar.x + n) = fuse(result.at<Vec3b>(src.y + m, src.x + n), result.at<Vec3b>(tar.y + m, tar.x + n),9*gauss(m*m+n*n));
		}
		else
		{
		result.at<Vec3b>(tar.y + m, tar.x + n) = result.at<Vec3b>(src.y + m, src.x + n);
		}
		mask.at<uchar>(tar.y + m, tar.x + n) = 255;
		}*/
		/*
		if (i == 0) continue;
		Point lasttar = pointManager.getPoint(anchorPoints[i - 1]);
		int ly = lasttar.y, lx = lasttar.x, ty = tar.y, tx = tar.x;
		int sx = src.x, sy = src.y;
		int L = max(lx, tx) - offset1, R = min(lx, tx) + offset2;
		int U = max(ly, ty) - offset1, D = min(ly, ty) + offset2;
		for (int m = U; m < D; m++)
		for (int n = L; n < R; n++)
		{
		double weight, W1, W2, W3, W4;
		if (tx > lx)
		W1 = R - n - 1, W2 = n - L;
		else
		W1 = n - L, W2 = R - n - 1;

		if (ty > ly)
		W3 = D - m - 1, W4 = m - U;
		else
		W3 = m - U, W4 = D - m - 1;
		weight = W1*W3 / (W1*W3 + W2*W4);
		int dy = m - ty, dx = n - tx;
		for (int ch = 0; ch < 3; ch++)
		{
		int ret;
		int v1 = result.at<Vec3b>(m, n)[ch];
		//int tmp1 = sy + dy;
		//int tmp2 = sx + dx;
		int v2 = result.at<Vec3b>(sy + dy, sx + dx)[ch];
		ret = v2 * weight - v1 * (1 - weight);
		result.at<Vec3b>(m, n)[ch] = ret;
		}
		}
		*/
	}
	free(sampleIndices);
}

int *StructurePropagation::BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat) {
	// do initialization
	int size = pointManager.getPropstackSize();
	anchorPoints.clear();
	anchorPoints.reserve(size);

	list<shared_ptr<MyNode>>::iterator itor;
	list<shared_ptr<MyNode>>::iterator end;
	pointManager.getPropstackItor(itor, end);
	// receive message sent from other neighbors
	for (; itor != end; itor++) {
		shared_ptr<MyNode> n = *itor;
		list<shared_ptr<Edge>> edges;
		n->getEdges(edges);
		// calculate message for next neighbor (the node that enqueued this node)
		calcMij(*n, n->getEdgeBegin(), mat, samplePoints);
	}

	// send updated message back to neighbors
	int *sampleIndices = (int*)malloc(size * sizeof(int));
	double *cur = (double*)malloc(samplePoints.size() * sizeof(double));

	list<shared_ptr<MyNode>>::reverse_iterator rev_itor;
	list<shared_ptr<MyNode>>::reverse_iterator rev_end;
	pointManager.getPropstackReverseItor(rev_itor, rev_end);
	for (int i = 0; rev_itor != rev_end; rev_itor++, i++) {
		shared_ptr<MyNode> n = *rev_itor;
		list<shared_ptr<Edge>>::iterator begin = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		list<shared_ptr<Edge>>::iterator itor = begin;

		anchorPoints.push_back(n->p);
		// calculate all the send-back message
		for (itor++; itor != end; itor++) {
			calcMij(*n, itor, mat, samplePoints);
		}
		int minIndex;
		double min = INT64_MAX;
		// calculate E1 for all possible xi
		for (int i = 0; i < samplePoints.size(); i++) {
			cur[i] = ks * calcEs(n->p, samplePoints[i]) + ki * calcEi(mat, n->p, samplePoints[i]);
		}
		// add up all messages sent to this node
		for (itor = begin; itor != end; itor++) {
			double **toMptr = (*itor)->getMbyTo(n->id);
			for (int i = 0; i < samplePoints.size(); i++) {
				cur[i] += (*toMptr)[i];
			}
		}
		// find out the optimal xi
		for (int i = 0; i < samplePoints.size(); i++) {
			if (cur[i] < min) {
				min = cur[i];
				minIndex = i;
			}
		}
		// cout << n->p.lineIndex << " " << n->p.pointIndex << " | ";
		// cout <<  i << " " << minIndex << endl;
		sampleIndices[i] = minIndex;
	}

	// release resources
	pointManager.getPropstackItor(itor, end);
	for (; itor != end; itor++) {
		shared_ptr<MyNode> n = *itor;
		list<shared_ptr<Edge>>::iterator edgeItor = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		for (; edgeItor != end; edgeItor++) {
			double **M = (*edgeItor)->getMbyFrom(n->id);
			if (*M != NULL) {
				free(*M);
			}
		}
	}
	free(cur);

	return sampleIndices;
}

/*int *StructurePropagation::BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat) {
// do initialization
int size = pointManager.getPropstackSize();
anchorPoints.clear();
anchorPoints.reserve(size);

list<shared_ptr<Node>>::iterator itor, begin, end;
pointManager.getPropstackItor(begin, end);
bool changed = true;
// receive message sent from other neighbors
while (changed) {
changed = false;
for (itor = begin; itor != end; itor++) {
shared_ptr<Node> n = *itor;
list<shared_ptr<Edge>> edges;
n->getEdges(edges);
list<shared_ptr<Edge>>::iterator edgeItor;
for (edgeItor = edges.begin(); edgeItor != edges.end(); edgeItor++) {
if (*edgeItor == NULL) {
list<shared_ptr<Edge>>::iterator tmpItor;
int calculable = true;
for (tmpItor = edges.begin(); tmpItor != edges.end(); tmpItor++) {
if (tmpItor != edgeItor) {
calculable = (*tmpItor != NULL);
}
}
if (calculable) {
calcMij(*n, edgeItor, mat, samplePoints);
changed = true;
}
}
}
}
}

for (itor = begin; itor != end; itor++) {
list<shared_ptr<Edge>>::iterator edgeItor;
list<shared_ptr<Edge>> edges;
for (edgeItor = edges.begin(); edgeItor != edges.end(); edgeItor++) {
assert(*(*edgeItor)->getMbyFrom((*itor)->id) != NULL);
}
}

// send updated message back to neighbors
int *sampleIndices = (int*)malloc(size * sizeof(int));
double *cur = (double*)malloc(samplePoints.size() * sizeof(double));

list<shared_ptr<Node>>::iterator rev_itor;
list<shared_ptr<Node>>::iterator rev_end;
pointManager.getPropstackItor(rev_itor, rev_end);
for (int i = 0; rev_itor != rev_end; rev_itor++, i++) {
shared_ptr<Node> n = *rev_itor;
list<shared_ptr<Edge>>::iterator begin = n->getEdgeBegin();
list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
list<shared_ptr<Edge>>::iterator itor = begin;

anchorPoints.push_back(n->p);
int minIndex;
double min = INT_MAX;
// calculate E1 for all possible xi
for (int i = 0; i < samplePoints.size(); i++) {
cur[i] = ks * calcEs(n->p, samplePoints[i]) + ki * calcEi(mat, n->p, samplePoints[i]);
}
// add up all messages sent to this node
for (itor = begin; itor == end; itor++) {
double **toMptr = (*itor)->getMbyTo(n->id);
for (int i = 0; i < samplePoints.size(); i++) {
cur[i] += (*toMptr)[i];
}
}
// find out the optimal xi
for (int i = 0; i < samplePoints.size(); i++) {
if (cur[i] < min) {
min = cur[i];
minIndex = i;
}
}
// cout << n->p.lineIndex << " " << n->p.pointIndex << " | ";
// cout <<  i << " " << minIndex << endl;
sampleIndices[i] = minIndex;
}

// release resources
pointManager.getPropstackItor(itor, end);
for (; itor != end; itor++) {
shared_ptr<Node> n = *itor;
list<shared_ptr<Edge>>::iterator edgeItor = n->getEdgeBegin();
list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
for (; edgeItor != end; edgeItor++) {
double **M = (*edgeItor)->getMbyFrom(n->id);
if (*M != NULL) {
free(*M);
}
}
}
free(cur);

return sampleIndices;
}*/


void StructurePropagation::calcMij(MyNode &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints) {
	double **Mptr = (*edgeItor)->getMbyFrom(n.id);
	list<shared_ptr<Edge>>::iterator end = n.getEdgeEnd();
	if (*Mptr == NULL) {
		*Mptr = (double*)malloc(samplePoints.size() * sizeof(double));
		memset(*Mptr, 0, samplePoints.size() * sizeof(double));
		for (int i = 0; i < samplePoints.size(); i++) {
			// calcalate Ei beforehand
			double E1 = ks * calcEs(n.p, samplePoints[i]) + ki * calcEi(mat, n.p, samplePoints[i]);
			// add up the message sent from Mki (k != j)
			list<shared_ptr<Edge>>::iterator itor = n.getEdgeBegin();
			double msg = 0.0;
			for (; itor != end; itor++) {
				if (itor != edgeItor) {
					double **toMptr = (*itor)->getMbyTo(n.id);
					if (*toMptr == NULL) {
						assert(0);
					}
					msg += (*toMptr)[i];
				}
			}
			PointPos tmpPos = pointManager.getPointPos((*edgeItor)->getAnother(n.id));
			for (int j = 0; j < samplePoints.size(); j++) {
				// try updating tne minimal value of each item in Mij
				double E2 = calcE2(mat, n.p, tmpPos, samplePoints[i], samplePoints[j]);
				if ((*Mptr)[j] == 0 || E1 + E2 + msg < (*Mptr)[j]) {
					(*Mptr)[j] = E1 + E2 + msg;
				}
			}
		}
	}
}

int *StructurePropagation::DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat) {

	double *M = (double *)malloc(2 * samplePoints.size() * sizeof(double));
	int *record = (int *)malloc(samplePoints.size() * anchorPoints.size() * sizeof(int));

	for (int i = 0; i < samplePoints.size(); i++) {
		M[i] = ks * calcEs(anchorPoints[0], samplePoints[i]) +
			ki * calcEi(mat, anchorPoints[0], samplePoints[i]);
	}

	int i, curOffset = 0, preOffset = 0;
	for (i = 1; i < anchorPoints.size(); i++) {
		curOffset = (i % 2) * samplePoints.size();
		preOffset = ((i + 1) % 2) * samplePoints.size();
		// calculate the min value for each xi
		// xi = j
		for (int j = 0; j < samplePoints.size(); j++) {
			double E1 = ks * calcEs(anchorPoints[i], samplePoints[j]) +
				ki * calcEi(mat, anchorPoints[i], samplePoints[j]);
			double min = INT_MAX;
			int minIndex;
			// choose optimal x(i-1)
			// x(i-1) = k
			for (int k = 0; k < samplePoints.size(); k++) {
				double tmp = calcE2(mat, anchorPoints[i], anchorPoints[i - 1], samplePoints[j], samplePoints[k]) + M[preOffset + k];
				if (tmp < min) {
					min = tmp;
					minIndex = k;
				}
			}
			record[samplePoints.size()*i + j] = minIndex;
			M[curOffset + j] = E1 + min;
		}
		/*for (int i = 0; i < samplePoints.size(); i++) {
		cout << M[preOffset + i] << ", ";
		}
		cout << endl;
		for (int j = 0; j < samplePoints.size(); j++) {
		cout << samplePoints.size()*i + j << ", ";
		}
		cout << endl;
		for (int i = 0; i < samplePoints.size(); i++) {
		cout << M[curOffset + i] << ", ";
		}
		cout << endl;*/
	}

	// find out the optimal xi of last anchor point
	int *sampleIndices = (int*)malloc(anchorPoints.size() * sizeof(int));
	double min = INT_MAX;
	for (int j = 0; j < samplePoints.size(); j++) {
		if (M[curOffset + j] < min) {
			sampleIndices[anchorPoints.size() - 1] = j;
			min = M[curOffset + j];
		}
	}

	// trace back
	for (int i = anchorPoints.size() - 2; i >= 0; i--) {
		sampleIndices[i] = record[samplePoints.size()*(i + 1) + sampleIndices[i + 1]];
	}

	free(M);
	free(record);

	return sampleIndices;
}

double StructurePropagation::calcE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2)
{
	int colLeft1, colLeft2, colRight1, colRight2;
	int rowUp1, rowUp2, rowDown1, rowDown2;
	Point p1 = pointManager.getPoint(i1);
	Point p2 = pointManager.getPoint(i2);
	Point px1 = pointManager.getPoint(xi1);
	Point px2 = pointManager.getPoint(xi2);
	// calculate the relative offsets of overlapping area's bounderies
	if (p1.x > p2.x) {
		colLeft1 = 0;
		colLeft2 = p1.x - p2.x;
		colRight1 = block_size - colLeft2;
		colRight2 = block_size;
	}
	else {
		colLeft2 = 0;
		colLeft1 = p2.x - p1.x;
		colRight2 = block_size - colLeft1;
		colRight1 = block_size;
	}

	if (p1.y > p2.y) {
		rowUp1 = 0;
		rowUp2 = p1.y - p2.y;
		rowDown1 = block_size - rowUp2;
		rowDown2 = block_size;
	}
	else {
		rowUp2 = 0;
		rowUp1 = p2.y - p1.y;
		rowDown2 = block_size - rowUp1;
		rowDown1 = block_size;
	}
	if (colRight1 >= 0 && colRight2 >= 0 && rowDown1 >= 0 && rowDown2 >= 0) {
		double ssd = 0.0;
		int cols = colRight1 - colLeft1;
		int rows = rowDown1 - rowUp1;
		// calculate the absolute cooordinates of boundaries
		int xOffset1 = colLeft1 + px1.x - block_size / 2;
		int xOffset2 = colLeft2 + px2.x - block_size / 2;
		int yOffset1 = rowUp1 + px1.y - block_size / 2;
		int yOffset2 = rowUp2 + px2.y - block_size / 2;
		for (int i = 0; i < rows; i++) {
			const uchar *ptr1 = mat.ptr<uchar>(i + yOffset1);
			const uchar *ptr2 = mat.ptr<uchar>(i + yOffset2);
			for (int j = 0; j < cols; j++) {
				double diff = ptr1[j + xOffset1] - ptr2[j + xOffset2];
				ssd += diff*diff;
			}
		}
		// do normlization
		return ssd / (cols * rows);
	}
	else {
		// no overlapping part
		return 0.0;
	}
}

/*double StructurePropagation::calcEs(const PointPos &i, const PointPos &xi) {
vector<Point> points1, points2;
vector<int> minDistance1, minDistance2;
Point pi = pointManager.getPoint(i);
Point pxi = pointManager.getPoint(xi);
int offsetx = pxi.x - pi.x;
int offsety = pxi.y - pi.y;
// get points of curve segment contained in patch
pointManager.getPointsinPatch(i, points1);
pointManager.getPointsinPatch(xi, points2);
minDistance1.resize(points1.size());
minDistance2.resize(points2.size());
// initialize minimal distance
for (int i = 0; i < points1.size(); i++) {
minDistance1[i] = INT_MAX;
}
for (int i = 0; i < points2.size(); i++) {
minDistance2[i] = INT_MAX;
}
// calculate the minimal distances for points in curve segments
for (int i = 0; i < points1.size(); i++) {
for (int j = 0; j < points2.size(); j++) {
int diffx = points1[i].x - points2[j].x + offsetx;
int diffy = points1[i].y - points2[j].y + offsety;
int distance = diffx*diffx + diffy*diffy;
if (distance < minDistance1[i]) {
minDistance1[i] = distance;
}
if (distance < minDistance2[j]) {
minDistance2[j] = distance;
}
}
}
int es1 = 0, es2 = 0;
for (int i = 0; i < minDistance1.size(); i++) {
es1 += minDistance1[i];
}
for (int i = 0; i < minDistance2.size(); i++) {
es2 += minDistance2[i];
}
return ((double)es1 + (double)es2) / minDistance1.size();
// return (double)es1 / minDistance1.size() + (double)es2 / minDistance2.size();
}*/

double StructurePropagation::calcEs(const PointPos &i, const PointPos &xi) {
	static vector<int> minDistance1, minDistance2;
	Point pi = pointManager.getPoint(i);
	Point pxi = pointManager.getPoint(xi);
	int offsetx = pxi.x - pi.x;
	int offsety = pxi.y - pi.y;
	// get points of curve segment contained in patch
	list<Point*> begin1, begin2;
	list<int> length1, length2;
	pointManager.getPointsinPatch(i, begin1, length1);
	pointManager.getPointsinPatch(xi, begin2, length2);

	assert(begin1.size() == length1.size());
	assert(begin2.size() == length2.size());

	list<int>::iterator lenItor;
	int totalLen1 = 0, totalLen2 = 0;
	for (lenItor = length1.begin(); lenItor != length1.end(); lenItor++) {
		totalLen1 += *lenItor;
	}
	for (lenItor = length2.begin(); lenItor != length2.end(); lenItor++) {
		totalLen2 += *lenItor;
	}
	minDistance1.resize(totalLen1);
	minDistance2.resize(totalLen2);


	// initialize minimal distance
	for (int i = 0; i < totalLen1; i++) {
		minDistance1[i] = INT_MAX;
	}
	for (int i = 0; i < totalLen2; i++) {
		minDistance2[i] = INT_MAX;
	}

	list<int>::iterator lenItor1, lenItor2;
	list<Point*>::iterator pointItor1, pointItor2;

	for (lenItor1 = length1.begin(), pointItor1 = begin1.begin(); lenItor1 != length1.end(); lenItor1++, pointItor1++) {
		Point *points1 = *pointItor1;
		for (int i = 0; i < *lenItor1; i++) {
			for (lenItor2 = length2.begin(), pointItor2 = begin2.begin(); lenItor2 != length2.end(); lenItor2++, pointItor2++) {
				for (int j = 0; j < *lenItor2; j++) {
					Point *points2 = *pointItor2;
					int diffx = points1[i].x - points2[j].x + offsetx;
					int diffy = points1[i].y - points2[j].y + offsety;
					int distance = diffx*diffx + diffy*diffy;
					if (distance < minDistance1[i]) {
						minDistance1[i] = distance;
					}
					if (distance < minDistance2[j]) {
						minDistance2[j] = distance;
					}
				}
			}
		}
	}

	int es1 = 0, es2 = 0;
	for (int i = 0; i < minDistance1.size(); i++) {
		es1 += minDistance1[i];
	}
	for (int i = 0; i < minDistance2.size(); i++) {
		es2 += minDistance2[i];
	}
	return ((double)es1 + (double)es2) / minDistance1.size();
}

double StructurePropagation::calcEi(const Mat &mat, const PointPos &i, const PointPos &xi) {
	if (pointManager.nearBoundary(i)) {
		int offset1 = block_size / 2;
		int offset2 = block_size - offset1;
		int ssd = 0;
		int overlappingPixelNum = 0;
		Point pi = pointManager.getPoint(i);
		Point pxi = pointManager.getPoint(xi);
		for (int i = -offset1; i < offset2; i++) {
			const uchar *ptri = mat.ptr<uchar>(i + pi.y);
			const uchar *ptrxi = mat.ptr<uchar>(i + pxi.y);
			for (int j = -offset1; j < offset2; j++) {
				// filter out invalid pixels located in masked area
				if (ptri[j + pi.x] != 0) {
					int diff = ptri[j + pi.x] - ptrxi[j + pxi.x];
					overlappingPixelNum++;
					ssd += diff*diff;
				}
			}
		}
		// do nomalization
		return (double)ssd / overlappingPixelNum;
	}
	else {
		// target patch is contained in mask area totally
		return 0.0;
	}
}