#include "PointManager.h"

#define getLineIndex(x) (x-1) >> 24
#define getPointIndex(x) (x-1) << 8 >> 8
#define visit(l, p) ((l << 24) | (p + 1))

ushort MyNode::totalNum = 0;

bool operator<(const PointPos &p1, const PointPos &p2) {
	return (p1.lineIndex == p2.lineIndex) ? p1.pointIndex < p2.pointIndex : p1.lineIndex < p2.lineIndex;
}

bool operator==(const Edge &e1, const Edge &e2) {
	return (e1.ni == e2.ni) & (e1.nj == e2.nj);
}


void PointManager::reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize, set<shared_ptr<list<int>>> &lineSets) {
	// clean up the remnant data
	this->linePoints = linePoints;
	this->mask = mask;
	this->blockSize = blockSize;
	lineEnds.clear();
	boundaryPoints.clear();
	intersectingMap.clear();
	outIntersectingMap.clear();

	// do preparation for the following operation
	Mat visitMat = Mat::zeros(mask.rows, mask.cols, CV_32SC1);
	bool inMask = false;
	Endpoints endpoints;
	for (int j = 0; j < linePoints.size(); j++) {
		int i;
		for (i = 0; i < linePoints[j].size(); i++) {
			int y = linePoints[j][i].y;
			int x = linePoints[j][i].x;
			if (y < 0 || y >= mask.rows || x < 0 || x >= mask.cols) {
				continue;
			}
			else if (mask.at<uchar>(y, x)) {
				// current line get out of the mask area, end this line
				if (inMask == true) {
					if (nearBoundary(linePoints[j][i], true)) {
						boundaryPoints.insert(PointPos(lineEnds.size(), i));
					}
					else {
						endpoints.endIndex = i;
						lineEnds.push_back(endpoints);
						inMask = false;
					}
					endpoints.endIndex = i;
					lineEnds.push_back(endpoints);
					inMask = false;
				}
				int visitRecord = visitMat.at<int>(y, x);
				int lineIndex = getLineIndex(visitRecord);
				int pointIndex = getPointIndex(visitRecord);
				// check if this point is an intersection^M
				if (visitRecord != 0 && lineIndex != lineEnds.size()) {
					// record the position info of overlapping points in the intersectingMap^M
					list<PointPos> *intersectingList = &outIntersectingMap[calcHashValue(x, y)];
					if (intersectingList->size() == 0) {
						intersectingList->push_back(PointPos(lineIndex, pointIndex));
					}
					intersectingList->push_back(PointPos(j, i));
				}
				else {
					// mark the point as visited
					visitMat.at<int>(y, x) = visit(j, i);
				}
			}
			else {
				// current patch crosses the boundary between img and mask
				if (nearBoundary(linePoints[j][i], false)) {
					boundaryPoints.insert(PointPos(lineEnds.size(), i));
				}
				// current line get into mask area, start a new line
				if (inMask == false) {
					endpoints.startIndex = i;
					endpoints.trueLineIndex = j;
					inMask = true;
				}
				int visitRecord = visitMat.at<int>(y, x);
				int lineIndex = getLineIndex(visitRecord);
				int pointIndex = getPointIndex(visitRecord);
				// check if this point is an intersection
				if (visitRecord != 0 && lineIndex != lineEnds.size()) {
					// record the position info of overlapping points in the intersectingMap
					list<PointPos> *intersectingList = &intersectingMap[calcHashValue(x, y)];
					if (intersectingList->size() == 0) {
						intersectingList->push_back(PointPos(lineIndex, pointIndex));
					}
					intersectingList->push_back(PointPos(lineEnds.size(), i));
				}
				else {
					// mark the point as visited
					visitMat.at<int>(y, x) = visit(lineEnds.size(), i);
				}
			}
		}
		// cope with the situation where line ends in mask area
		if (inMask == true) {
			inMask = false;
			endpoints.endIndex = i;
			lineEnds.push_back(endpoints);
		}
	}

	vector<shared_ptr<list<int>>> lineSetRecord(lineEnds.size());
	map<int, list<PointPos>>::iterator mapItor;
	list<PointPos>::iterator listItor;
	for (int i = 0; i < lineEnds.size(); i++) {
		shared_ptr<list<int>> ptr = make_shared<list<int>>();
		ptr->push_back(i);
		lineSetRecord[i] = ptr;
	}
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		shared_ptr<list<int>> ptr = NULL;
		for (listItor = mapItor->second.begin(); listItor != mapItor->second.end(); listItor++) {
			if (lineSetRecord[listItor->lineIndex] != NULL) {
				if (ptr == NULL) {
					ptr = lineSetRecord[listItor->lineIndex];
				}
				else {
					// merge two line sets
					ptr->insert(ptr->end(), lineSetRecord[listItor->lineIndex]->begin(), lineSetRecord[listItor->lineIndex]->end());
				}
			}
		}
		list<int>::iterator itor;
		for (itor = ptr->begin(); itor != ptr->end(); itor++) {
			lineSetRecord[*itor] = ptr;
		}
	}

	for (int i = 0; i < lineSetRecord.size(); i++) {
		lineSets.insert(lineSetRecord[i]);
	}
}

bool PointManager::nearBoundary(const Point &p, bool isSample) {
	int leftBound = MAX(p.x - blockSize / 2, 0);
	int rightBound = MIN(p.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(p.y - blockSize / 2, 0);
	int downBound = MIN(p.y + blockSize - blockSize / 2, mask.rows);
	const uchar *upPtr = mask.ptr<uchar>(upBound);
	const uchar *downPtr = mask.ptr<uchar>(downBound - 1);
	// check if the mask boundary crosses up and down boundary of the patch
	for (int i = leftBound; i < rightBound; i++) {
		if (!upPtr[i] == isSample || !downPtr[i] == isSample) {
			return true;
		}
	}
	// check if the mask boundary crosses left and right boundary of the patch
	for (int i = upBound + 1; i < downBound - 1; i++) {
		if (!mask.at<uchar>(i, leftBound) == isSample || !mask.at<uchar>(i, rightBound - 1) == isSample) {
			return true;
		}
	}
	return false;
}

inline int PointManager::calcHashValue(int x, int y) {
	// calculate the hash value of the visitedMap
	return x + y * mask.cols;
}

inline Point PointManager::getPoint(PointPos p) {
	//get Point Object by PointPos object
	return linePoints[lineEnds[p.lineIndex].trueLineIndex][p.pointIndex];
}

bool PointManager::nearBoundary(PointPos p) {
	// check if the patch crosses the boundary
	return boundaryPoints.count(p);
}

void PointManager::getPointsinPatch(PointPos p, vector<Point> &ret) {
	// get all points of the line segment contained in this patch
	Point center = getPoint(p);
	int leftBound = MAX(center.x - blockSize / 2, 0);
	int rightBound = MIN(center.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(center.y - blockSize / 2, 0);
	int downBound = MIN(center.y + blockSize - blockSize / 2, mask.rows);
	int hashValue = calcHashValue(center.x, center.y);
	list<PointPos> pointPositions;
	// check if the the anchor point is an intersaction
	// if it is, points of sevaral line segments will be returned
	if (intersectingMap.count(hashValue)) {
		pointPositions = intersectingMap[hashValue];
	}
	else {
		pointPositions.push_back(p);
	}
	for (list<PointPos>::iterator p = pointPositions.begin(); p != pointPositions.end(); p++) {
		Endpoints endPoints = lineEnds[p->lineIndex];
		Point *points = &linePoints[endPoints.trueLineIndex][0];
		int beginIndex = p->pointIndex;
		//find the start index of the line segment
		for (int i = p->pointIndex; i >= 0; i--) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				beginIndex = i + 1;
				break;
			}
		}
		// get anchor points forward
		for (int i = beginIndex; i < p->pointIndex; i++) {
			ret.push_back(points[i]);
		}
		// get anchor points backward
		for (int i = p->pointIndex; i < linePoints[endPoints.trueLineIndex].size(); i++) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				break;
			}
			else {
				ret.push_back(points[i]);
			}
		}
	}
	int i = 0;
	i++;
}

void PointManager::getPointsinPatch(const PointPos &p, list<Point*> &begin, list<int> &length) {
	// get all points of the line segment contained in this patch
	Point center = getPoint(p);
	int leftBound = MAX(center.x - blockSize / 2, 0);
	int rightBound = MIN(center.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(center.y - blockSize / 2, 0);
	int downBound = MIN(center.y + blockSize - blockSize / 2, mask.rows);
	int hashValue = calcHashValue(center.x, center.y);
	list<PointPos> pointPositions;
	bool inMask = true;
	// check if the the anchor point is an intersaction
	// if it is, points of sevaral line segments will be returned
	if (intersectingMap.count(hashValue)) {
		pointPositions = intersectingMap[hashValue];
	}
	else if (outIntersectingMap.count(hashValue)) {
		pointPositions = outIntersectingMap[hashValue];
		inMask = false;
	}
	else {
		pointPositions.push_back(p);
	}
	for (list<PointPos>::iterator p = pointPositions.begin(); p != pointPositions.end(); p++) {
		int trueLineIndex = (inMask) ? lineEnds[p->lineIndex].trueLineIndex : p->lineIndex;
		Point *points = &linePoints[trueLineIndex][0];
		int beginIndex = p->pointIndex;
		//find the start index of the line segment
		for (int i = p->pointIndex; i >= 0; i--) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				beginIndex = i + 1;
				break;
			}
		}
		begin.push_back(points + beginIndex);
		// get anchor points backward
		int i;
		for (i = p->pointIndex; i < linePoints[trueLineIndex].size(); i++) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				length.push_back(i - beginIndex);
				break;
			}
		}
		if (i == linePoints[trueLineIndex].size()) {
			length.push_back(linePoints[trueLineIndex].size() - beginIndex);
		}
	}
}

void PointManager::constructBPMap(list<int> &line) {
	map<int, list<PointPos>>::iterator mapItor;
	list<shared_ptr<MyNode>> BFSstack;
	vector<vector<ushort>> pointVisitedMarks(linePoints.size());
	vector<list<shared_ptr<MyNode>>> nodeListBucket(4);
	set<int> intersectionSet;

	nodes.clear();
	propagationStack.clear();
	MyNode::totalNum = 0;

	// initialize the visit map
	pointVisitedMarks.resize(linePoints.size());
	for (int i = 0; i < linePoints.size(); i++) {
		pointVisitedMarks[i].resize(linePoints[i].size());
	}

	int total = 0;
	list<int>::iterator itor;
	for (itor = line.begin(); itor != line.end(); itor++) {
		total += lineEnds[*itor].endIndex - lineEnds[*itor].startIndex;
		intersectionSet.insert(*itor);
	}
	// reserve enough space for the node table
	nodes.reserve(total / blockSize + 1);
	// skip the entry numbered with 0
	nodes.resize(1);

	// enqueue all the intersections
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		if (intersectionSet.count(mapItor->second.begin()->lineIndex) == 0) {
			continue;
		}
		// enqueue the intersection (choose one point's position to represent all)
		BFSstack.push_back(make_shared<MyNode>(*(mapItor->second.begin())));
		// mark all intersecting points as visited
		list<PointPos>::iterator listItor = mapItor->second.begin();
		for (; listItor != mapItor->second.end(); listItor++) {
			pointVisitedMarks[lineEnds[listItor->lineIndex].trueLineIndex][listItor->pointIndex] = MyNode::totalNum;
		}
	}

	// enqueue all the neighbor nodes of intersections
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		if (intersectionSet.count(mapItor->second.begin()->lineIndex) == 0) {
			continue;
		}
		shared_ptr<MyNode> n = *BFSstack.begin();
		list<PointPos>::iterator listItor = mapItor->second.begin();
		int neighborNum = 0;
		for (; listItor != mapItor->second.end(); listItor++) {
			neighborNum += addNeighbor(*n, *listItor, pointVisitedMarks, BFSstack);
		}
		// Enlarge nodeListBucket if necessary
		if (neighborNum > nodeListBucket.size()) {
			nodeListBucket.resize(mapItor->second.size() * 2);
		}
		nodeListBucket[neighborNum - 1].push_front(n);
		nodes.push_back(nodeListBucket[neighborNum - 1].begin());
		BFSstack.pop_front();
	}

	// start propagation
	while (BFSstack.size()) {
		shared_ptr<MyNode> n = *BFSstack.begin();
		int neighborNum = addNeighbor(*n, n->p, pointVisitedMarks, BFSstack);
		nodeListBucket[neighborNum - 1].push_front(n);
		nodes.push_back(nodeListBucket[neighborNum - 1].begin());
		BFSstack.pop_front();
	}

	//generate the sequence for message sending
	while (nodeListBucket[0].size() > 0) {
		shared_ptr<MyNode> n = *nodeListBucket[0].begin();
		if (n->getEdgeNum() != 1) {
			assert(n->getEdgeNum() == 1);
		}
		list<shared_ptr<Edge>>::iterator eItor = n->getEdgeBegin();
		nodes[n->id] = propagationStack.insert(propagationStack.end(), n);
		nodeListBucket[0].pop_front();
		// degrade the adjacent node
		int id = (*eItor)->getAnother(n->id);
		n = *nodes[id];
		int edgeNum = n->getEdgeNum();
		n->eraseEdge(*eItor);
		nodeListBucket[edgeNum - 1].erase(nodes[id]);
		if (edgeNum > 1) {
			nodes[id] = nodeListBucket[edgeNum - 2].insert(nodeListBucket[edgeNum - 2].end(), n);
		}
		else {
			nodes[id] = propagationStack.insert(propagationStack.end(), n);
		}
	}
}

int PointManager::addNeighbor(MyNode &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<shared_ptr<MyNode>> &BFSstack) {
	Endpoints endpoints = lineEnds[pos.lineIndex];
	int lineIndex = endpoints.trueLineIndex;
	int pointIndex = pos.pointIndex;
	int prePointIndex = pointIndex - blockSize / 2;
	int nextPointIndex = pointIndex + blockSize / 2;
	int neighborNum = 0;

	// check the point before current anchor point
	if (prePointIndex >= endpoints.startIndex) {
		int i;
		// try choosing a existed anchor point as its neighbor
		for (i = prePointIndex; i < pointIndex; i++) {
			if (visitedMark[lineIndex][i]) {
				if (nodes.size() > visitedMark[lineIndex][i]) {
					// add an edge between two points
					shared_ptr<Edge> tmpEdge = make_shared<Edge>(n.id, visitedMark[lineIndex][i]);
					n.push_front(tmpEdge);
					(*nodes[visitedMark[lineIndex][i]])->push_back(tmpEdge);
				}
				break;
			}
		}
		// no existed point can be chosen, construct a new anchor point and enqueue it
		if (i == pointIndex) {
			BFSstack.push_back(make_shared<MyNode>(PointPos(pos.lineIndex, prePointIndex)));
			visitedMark[lineIndex][prePointIndex] = MyNode::totalNum;
		}
		neighborNum++;
	}

	// check the point behind current anchor point
	if (nextPointIndex < endpoints.endIndex) {
		int i;
		// try choosing a existed anchor point as its neighbor
		for (i = nextPointIndex; i > pointIndex; i--) {
			if (visitedMark[lineIndex][i]) {
				if (nodes.size() > visitedMark[lineIndex][i]) {
					// add an edge between two points
					shared_ptr<Edge> tmpEdge = make_shared<Edge>(n.id, visitedMark[lineIndex][i]);
					n.push_front(tmpEdge);
					(*nodes[visitedMark[lineIndex][i]])->push_back(tmpEdge);
				}
				break;
			}
		}
		// no existed point can be chosen, construct a new anchor point and enqueue it
		if (i == pointIndex) {
			BFSstack.push_back(make_shared<MyNode>(PointPos(pos.lineIndex, nextPointIndex)));
			visitedMark[lineIndex][nextPointIndex] = MyNode::totalNum;
		}
		neighborNum++;
	}
	return neighborNum;
}

void PointManager::getPropstackItor(list<shared_ptr<MyNode>>::iterator &begin, list<shared_ptr<MyNode>>::iterator &end) {
	begin = propagationStack.begin();
	end = propagationStack.end();
}

void PointManager::getPropstackReverseItor(list<shared_ptr<MyNode>>::reverse_iterator &begin, list<shared_ptr<MyNode>>::reverse_iterator &end) {
	begin = list<shared_ptr<MyNode>>::reverse_iterator(propagationStack.end())++;
	end = list<shared_ptr<MyNode>>::reverse_iterator(propagationStack.begin());
}

void PointManager::getSamplePoints(vector<PointPos> &samples, int sampleStep, list<int> &line) {
	if (lineEnds.size() == 0) {
		return;
	}
	samples.clear();
	Endpoints endpoints = lineEnds[*line.begin()];

	// reserve enough space for samples
	int total = 0;
	int curLine = -1;
	list<int>::iterator itor;
	for (itor = line.begin(); itor != line.end(); itor++) {
		if (lineEnds[*itor].trueLineIndex != curLine) {
			curLine = lineEnds[*itor].trueLineIndex;
			total += linePoints[curLine].size();
		}
		total -= (lineEnds[*itor].endIndex - lineEnds[*itor].startIndex);
	}
	samples.reserve(total / sampleStep);

	// get samples from all line segments outside the mask area
	// sampling step = sampleStep
	itor = line.begin();
	for (int i = 0; i < linePoints.size(); i++) {
		// +blocksize: ensure all samples have complete line segments
		if (endpoints.trueLineIndex != i) {
			continue;
		}
		int beginIndex = blockSize;
		int endIndex;
		while (endpoints.trueLineIndex == i) {
			endIndex = endpoints.startIndex;
			for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
				if (j == beginIndex) {
					int c = 0;
					c++;
				}
				if (!nearBoundary(linePoints[i][j], true)) {
					samples.push_back(PointPos(*itor, j));
				}
			}
			beginIndex = endpoints.endIndex;
			itor++;
			if (itor == line.end()) {
				break;
			}
			endpoints = lineEnds[*itor];
		}
		// -blocksize: ensure all samples have complete line segments
		endIndex = linePoints[i].size() - blockSize;
		for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
			if (!nearBoundary(linePoints[i][j], true)) {
				samples.push_back(PointPos(*(list<int>::reverse_iterator(itor)++), j));
			}
		}
	}
	samples.shrink_to_fit();
}

void PointManager::getSamplePoints(vector<PointPos> &samples, int sampleStep) {
	if (lineEnds.size() == 0) {
		return;
	}
	samples.clear();
	int lineIndex = 0;
	Endpoints endpoints = lineEnds[0];

	// reserve enough space for samples
	int total = 0;
	for (int i = 0; i < linePoints.size(); i++) {
		total += linePoints[i].size();
	}
	for (int i = 0; i < lineEnds.size(); i++) {
		total -= (lineEnds[i].endIndex - lineEnds[i].startIndex);
	}
	samples.reserve(total / sampleStep);

	// get samples from all line segments outside the mask area
	// sampling step = sampleStep
	for (int i = 0; i < linePoints.size(); i++) {
		// +blocksize: ensure all samples have complete line segments
		int beginIndex = blockSize;
		int endIndex;
		while (endpoints.trueLineIndex == i) {
			endIndex = endpoints.startIndex;
			for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
				if (j == beginIndex) {
					int c = 0;
					c++;
				}
				if (!nearBoundary(linePoints[i][j], true)) {
					samples.push_back(PointPos(lineIndex, j));
				}
			}
			beginIndex = endpoints.endIndex;
			++lineIndex;
			if (lineIndex >= lineEnds.size()) {
				break;
			}
			endpoints = lineEnds[lineIndex];
		}
		// -blocksize: ensure all samples have complete line segments
		endIndex = linePoints[i].size() - blockSize;
		for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
			if (!nearBoundary(linePoints[i][j], true)) {
				samples.push_back(PointPos(lineIndex - 1, j));
			}
		}
	}
	samples.shrink_to_fit();
}

void PointManager::getAnchorPoints(vector<PointPos> &anchors, list<int> &line) {
	// only called by DP 
	// BP get anchor points by calling constructBPMap and getPropstackItor sequentially
	anchors.clear();
	Endpoints endpoints = lineEnds[*line.begin()];

	// reserve enough space for anchors
	int total = 0;
	list<int>::iterator itor;
	for (itor = line.begin(); itor != line.end(); itor++) {
		total += (lineEnds[*(itor)].endIndex - lineEnds[*(itor)].startIndex);
	}
	anchors.reserve(total / (blockSize / 2));

	// get anchors from all iline segments inside the mask area
	// sampling step = blockSize / 2
	itor = line.begin();
	for (int i = 0; i < linePoints.size(); i++) {
		while (endpoints.trueLineIndex == i) {
			for (int j = endpoints.startIndex; j < endpoints.endIndex; j += blockSize / 2) {
				anchors.push_back(PointPos(*itor, j));
			}
			++itor;
			if (itor == line.end()) {
				break;
			}
			endpoints = lineEnds[*itor];
		}
	}
	anchors.shrink_to_fit();
}