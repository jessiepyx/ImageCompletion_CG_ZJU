#include <opencv2/opencv.hpp>
#include <iostream>
#include <algorithm>
#include <queue>

using namespace cv;
using namespace std;

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

//全黑色是0，全白色是255
// mask: 二值化的mask图像
// Linemask：暂时理解为结构线
// mat：是之前带有mask的没有进行纹理补全的结果
// result：最后输出的结果
void TextureCompletion2(Mat1b _mask, Mat1b LineMask, const Mat &mat, Mat &result)
{
	int N = _mask.rows;
	int M = _mask.cols;
	int knowncount = 0;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			knowncount += (_mask.at<uchar>(i, j) == 255);
			//统计输入mask中纯白色像素点的个数
		}
	//做了一种优化处理，判断是黑色点多还是白色点多，从而进行后面的操作
	// mask部分是0白色？？
	if (knowncount * 2< N * M)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				_mask.at<uchar>(i, j) = 255 - _mask.at<uchar>(i, j);
	}

	//新建一个my_mask和sum_diff
	vector<vector<int> >my_mask(N, vector<int>(M, 0)), sum_diff(N, vector<int>(M, 0));

	//Linemask扩大这后面白色的255*100变成灰色，黑色依旧是0
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			LineMask.at<uchar>(i, j) = LineMask.at<uchar>(i, j) * 100;

	result = mat.clone();
	imshow("mask", _mask);
	imshow("linemask", LineMask);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			//mymask对应于mask（mask中的黑色遮挡部分mymask为0，mask白色部分mymask为1）
			my_mask[i][j] = (_mask.at<uchar>(i, j) == 255);
			//如果mymask中的一个位置坐标既是遮挡，又是LineMask中的灰色部分，则标注为2	
			if (my_mask[i][j] == 0 && LineMask.at<uchar>(i, j) > 0)
			{
				my_mask[i][j] = 2;
			}
		}
	/*
	my_mask的结构
	1 1 1 1 1 1 1
	1 1 1 1 1 1 1
	1 0 0 0 0 0 1
	1 0 0 0 2 0 1  ---结构线
	1 0 2 2 2 0 1  ---结构线
	1 0 0 0 0 0 1
	1 1 1 1 1 1 1
	*/

	int bs = 5;
	int step = 1 * bs;
	auto usable(my_mask);	//自动生成了一个和mymask相同类型的变量
	int to_fill = 0;	//mymask中未被填充的阴影遮挡的部分（非结构线）
	int filled = 0;		//mymask中未被填充的阴影遮挡的部分（非结构线）
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			to_fill += (my_mask[i][j] == 0);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			//如果my_mask[i][j] == 1说明不需要填充则继续
			if (my_mask[i][j] == 1)
				continue;
			//对于mymask中需要被填充的地方
			//在一个step的矩形邻域内，需要把usable标记为2
			//usable[k][l] == 2说明需要被填充
			//（我的理解是在原来的mask周围扩大了需要补全纹理的范围）
			int k0 = max(0, i - step), k1 = min(N - 1, i + step);
			int l0 = max(0, j - step), l1 = min(M - 1, j + step);
			for (int k = k0; k <= k1; k++)
				for (int l = l0; l <= l1; l++)
					usable[k][l] = 2;
		}
	//按照usable中2的地方生成一个黑白图，其中白色是需要填充的地方值为2
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
						//略过不需要填充的地方以及轮廓线部分
						if (my_mask[i][j] != 0) continue;
						//此时my_mask[i][j]==0
						bool edge = false;
						int k0 = max(0, i - 1), k1 = min(N - 1, i + 1);
						int l0 = max(0, j - 1), l1 = min(M - 1, j + 1);
						//取到像素点的一个小邻域8个像素点，如果这个邻域内的点有一个是1则最后edge==true
						/*
						1 1 1
						1 0 1
						1 1 1
						*/
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								edge |= (my_mask[k][l] == 1);	//或等于 edge = edge | (my_mask==1);
						if (!edge) continue;
						//如果edge==true说明当前像素点是边界点
						//------猜测后面需要对这个像素点进行融合运算！-------
						k0 = max(0, i - bs), k1 = min(N - 1, i + bs);
						l0 = max(0, j - bs), l1 = min(M - 1, j + bs);
						int tmpcnt = 0;
						//此时取到当前像素点周围的一个step大小的矩形邻域
						//tmpcnt计算了这个矩形邻域内不需要填充的像素点的个数
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								tmpcnt += (my_mask[k][l] == 1);
						if (tmpcnt > cnt)
						{
							cnt = tmpcnt;
							x = i;
							y = j;
						}
						//结束for循环的时候xy记录了周围不需要填充像素最多的那个像素点的坐标
					}
				//如果cnt==-1说明所有edge都是false，也就是说所有mymask[i，j]都是1都是不需要填充，跳出while
				if (cnt == -1) break;

				bool debug = false;
				bool debug2 = false;


				//这部分再次遍历全图；比较一个邻域内和整张图片其他邻域内是否有相似的块
				int k0 = min(x, bs), k1 = min(N - 1 - x, bs);
				int l0 = min(y, bs), l1 = min(M - 1 - y, bs);
				int sx, sy, min_diff = INT_MAX;	//最大的int值
				for (int i = step; i + step < N; i += step)
					for (int j = step; j + step < M; j += step)
					{
						//通过usable找到最近的不需要填充的像素点
						//如果==2说明这里的纹理不可用
						if (usable[i][j] == 2)continue;
						int tmp_diff = 0;
						//取到xy周围step的矩形邻域
						for (int k = -k0; k <= k1; k++)
							for (int l = -l0; l <= l1; l++)
							{
								//printf("%d %d %d %d %d %d\n", i + k, j + l, x + k, y + l, N, M);
								//ij表示可以用来比较的不需要填充纹理的坐标点
								//xy表示当前需要被填充的点，由之前的for循环生成
								//[x + k][y + l]表示xy的step邻域内的某点
								//[i + k][j + l]表示ij的step邻域内的某点
								if (my_mask[x + k][y + l] != 0)
									tmp_diff += dist(result.at<Vec3b>(i + k, j + l), result.at<Vec3b>(x + k, y + l));
								//tmp_diff计算了这两个对应点之间，RGB值的差异；显然需要全图搜索找到一个最小的tmpdiff，这说明这两块邻域最像


								//--------------------------这里似乎有文章？？？没有规定这个对应点的范围，没有考虑轮廓线--------------------//



							}
						sum_diff[i][j] = tmp_diff;
						if (min_diff > tmp_diff)
						{
							sx = i;
							sy = j;
							min_diff = tmp_diff;
						}
						//结束循环的时候，得到的是对比xy有最小tmpdiff的点的坐标sx，sy
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
				//用（sx，sy）周围的点的RGB值填充xy周围需要被填充的点
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

vector<vector<int> > getContous(Mat1b LineMask) {

	vector<vector<int> > myMap;
	return myMap;
}


int main()
{
	//四个输入：mask，line，
	
	//读入原图
	Mat3b img = imread("Origin.png");//5,1
	
	//读入二值化的mask图像
	Mat1b mask = Mat::zeros(img.rows, img.cols, CV_8UC1);
	mask = imread("Mask.bmp", 0);
	
	threshold(mask, mask, 125, 255, CV_THRESH_BINARY_INV);
	imshow("img", img);
	waitKey(10);
	imshow("mask", mask);
	waitKey(10);
	//生成带有mask但是没有进行补全的图
	Mat3b result;
	result.zeros(img.size());
	img.copyTo(result, mask);
	imshow("result", result);
	waitKey(10);
	//读入linemask
	Mat1b Linemask = Mat::zeros(img.rows, img.cols, CV_8UC1);
	//最终结果变量
	Mat3b finalResult(img.size());
	img.copyTo(finalResult);

	//TextureCompletion2(mask, Linemask, result, finalResult);

	system("pause");
}


