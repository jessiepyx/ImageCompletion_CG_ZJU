#ifndef POINT_MANAGER_H
#define POINT_MANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <unordered_map>
#include <memory>

using namespace std;
using namespace cv;

class Endpoints {
public:
	int trueLineIndex;
	int startIndex;
	int endIndex;
};

class PointPos {
public:
	int lineIndex;
	int pointIndex;

	PointPos(int lineIndex = -1, int pointIndex = -1) : lineIndex(lineIndex), pointIndex(pointIndex) {
	}
};

class Edge;

class MyNode {
private:
	list<shared_ptr<Edge>> edges;
	int edgeNum;
public:
	const PointPos p;
	const ushort id;
	static ushort totalNum;

	MyNode(PointPos p) : p(p), id(++totalNum) {
		edgeNum = 0;
	}
	void push_front(const shared_ptr<Edge> &e) {
		edges.push_front(e);
		edgeNum++;
	}
	void push_back(const shared_ptr<Edge> &e) {
		edges.push_back(e);
		edgeNum++;
	}
	void eraseEdge(shared_ptr<Edge> &e) {
		edges.remove(e);
		edges.push_back(e);
		edgeNum--;
	}
	list<shared_ptr<Edge>>::iterator getEdgeBegin() {
		return edges.begin();
	}
	list<shared_ptr<Edge>>::iterator getEdgeEnd() {
		return edges.end();
	}
	void getEdges(list<shared_ptr<Edge>> &edges) {
		edges = this->edges;
	}
	int getEdgeNum() {
		return edgeNum;
	}
};

class Edge {
private:
	double *Mij;
	double *Mji;
public:
	ushort ni;
	ushort nj;

	Edge(ushort ni, ushort nj) : ni(ni), nj(nj) {
		Mij = Mji = nullptr;
	}

	Edge() {
		Mij = Mji = nullptr;
	}

	inline ushort getAnother(ushort n) {
		if (n == ni) {
			return nj;
		}
		else {
			return ni;
		}
	}

	inline double **getMbyFrom(ushort from) {
		if (from == ni) {
			return &Mij;
		}
		else {
			return &Mji;
		}

	}

	inline double **getMbyTo(ushort to) {
		if (to == ni) {
			return &Mji;
		}
		else  {
			return &Mij;
		}
	}

};

class PointManager {
public:
	PointManager() = default;
	void reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize, set<shared_ptr<list<int>>> &lineSets);
	Point getPoint(PointPos p);
	bool nearBoundary(PointPos p);
	void getPointsinPatch(PointPos p, vector<Point> &ret);
	void getPointsinPatch(const PointPos &p, list<Point*> &begin, list<int> &length);
	void getSamplePoints(vector<PointPos> &samples, int sampleStep, list<int> &line);
	void getSamplePoints(vector<PointPos> &samples, int sampleStep);
	void constructBPMap(list<int> &line);
	void getAnchorPoints(vector<PointPos> &anchors, list<int> &line);
	void getPropstackItor(list<shared_ptr<MyNode>>::iterator &begin, list<shared_ptr<MyNode>>::iterator &end);
	void getPropstackReverseItor(list<shared_ptr<MyNode>>::reverse_iterator &begin, list<shared_ptr<MyNode>>::reverse_iterator &end);
	int getPropstackSize() {
		return propagationStack.size();
	}
	PointPos getPointPos(ushort id) {
		return (*nodes[id])->p;
	}
	shared_ptr<MyNode> getNode(ushort id) {
		return *nodes[id];
	}

private:
	vector<vector<Point>> linePoints; //记录用户绘制的点的信息
	Mat1b mask;
	int blockSize;
	vector<Endpoints> lineEnds; //用于记录经过PointManager再次划分后的线段的首尾信息
	set<PointPos> boundaryPoints; //用于记录所在patch与边界重叠的锚点
	map<int, list<PointPos>> intersectingMap; //用于记录交点，键值为根据交点的真实坐标计算出的hash值
	map<int, list<PointPos>> outIntersectingMap;
	vector<list<shared_ptr<MyNode>>::iterator> nodes; //一张根据Node id查找node的表，记录着Node对象在双向链表中的迭代器
	list<shared_ptr<MyNode>> propagationStack; //记录BP算法中信息传递的顺序

	bool nearBoundary(const Point &p, bool isSample);
	int calcHashValue(int x, int y);
	int addNeighbor(MyNode &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<shared_ptr<MyNode>> &neighbors);
};

#endif /* POINT_MANAGER_H */