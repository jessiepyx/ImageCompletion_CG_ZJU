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

//ȫ��ɫ��0��ȫ��ɫ��255
// mask: ��ֵ����maskͼ��
// Linemask����ʱ���Ϊ�ṹ��
// mat����֮ǰ����mask��û�н�������ȫ�Ľ��
// result���������Ľ��
void TextureCompletion2(Mat1b _mask, Mat1b LineMask, const Mat &mat, Mat &result)
{
	int N = _mask.rows;
	int M = _mask.cols;
	int knowncount = 0;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			knowncount += (_mask.at<uchar>(i, j) == 255);
			//ͳ������mask�д���ɫ���ص�ĸ���
		}
	//����һ���Ż������ж��Ǻ�ɫ��໹�ǰ�ɫ��࣬�Ӷ����к���Ĳ���
	// mask������0��ɫ����
	if (knowncount * 2< N * M)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				_mask.at<uchar>(i, j) = 255 - _mask.at<uchar>(i, j);
	}

	//�½�һ��my_mask��sum_diff
	vector<vector<int> >my_mask(N, vector<int>(M, 0)), sum_diff(N, vector<int>(M, 0));

	//Linemask����������ɫ��255*100��ɻ�ɫ����ɫ������0
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			LineMask.at<uchar>(i, j) = LineMask.at<uchar>(i, j) * 100;

	result = mat.clone();
	imshow("mask", _mask);
	imshow("linemask", LineMask);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			//mymask��Ӧ��mask��mask�еĺ�ɫ�ڵ�����mymaskΪ0��mask��ɫ����mymaskΪ1��
			my_mask[i][j] = (_mask.at<uchar>(i, j) == 255);
			//���mymask�е�һ��λ����������ڵ�������LineMask�еĻ�ɫ���֣����עΪ2	
			if (my_mask[i][j] == 0 && LineMask.at<uchar>(i, j) > 0)
			{
				my_mask[i][j] = 2;
			}
		}
	/*
	my_mask�Ľṹ
	1 1 1 1 1 1 1
	1 1 1 1 1 1 1
	1 0 0 0 0 0 1
	1 0 0 0 2 0 1  ---�ṹ��
	1 0 2 2 2 0 1  ---�ṹ��
	1 0 0 0 0 0 1
	1 1 1 1 1 1 1
	*/

	int bs = 5;
	int step = 1 * bs;
	auto usable(my_mask);	//�Զ�������һ����mymask��ͬ���͵ı���
	int to_fill = 0;	//mymask��δ��������Ӱ�ڵ��Ĳ��֣��ǽṹ�ߣ�
	int filled = 0;		//mymask��δ��������Ӱ�ڵ��Ĳ��֣��ǽṹ�ߣ�
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			to_fill += (my_mask[i][j] == 0);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			//���my_mask[i][j] == 1˵������Ҫ��������
			if (my_mask[i][j] == 1)
				continue;
			//����mymask����Ҫ�����ĵط�
			//��һ��step�ľ��������ڣ���Ҫ��usable���Ϊ2
			//usable[k][l] == 2˵����Ҫ�����
			//���ҵ��������ԭ����mask��Χ��������Ҫ��ȫ����ķ�Χ��
			int k0 = max(0, i - step), k1 = min(N - 1, i + step);
			int l0 = max(0, j - step), l1 = min(M - 1, j + step);
			for (int k = k0; k <= k1; k++)
				for (int l = l0; l <= l1; l++)
					usable[k][l] = 2;
		}
	//����usable��2�ĵط�����һ���ڰ�ͼ�����а�ɫ����Ҫ���ĵط�ֵΪ2
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
						//�Թ�����Ҫ���ĵط��Լ������߲���
						if (my_mask[i][j] != 0) continue;
						//��ʱmy_mask[i][j]==0
						bool edge = false;
						int k0 = max(0, i - 1), k1 = min(N - 1, i + 1);
						int l0 = max(0, j - 1), l1 = min(M - 1, j + 1);
						//ȡ�����ص��һ��С����8�����ص㣬�����������ڵĵ���һ����1�����edge==true
						/*
						1 1 1
						1 0 1
						1 1 1
						*/
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								edge |= (my_mask[k][l] == 1);	//����� edge = edge | (my_mask==1);
						if (!edge) continue;
						//���edge==true˵����ǰ���ص��Ǳ߽��
						//------�²������Ҫ��������ص�����ں����㣡-------
						k0 = max(0, i - bs), k1 = min(N - 1, i + bs);
						l0 = max(0, j - bs), l1 = min(M - 1, j + bs);
						int tmpcnt = 0;
						//��ʱȡ����ǰ���ص���Χ��һ��step��С�ľ�������
						//tmpcnt������������������ڲ���Ҫ�������ص�ĸ���
						for (int k = k0; k <= k1; k++)
							for (int l = l0; l <= l1; l++)
								tmpcnt += (my_mask[k][l] == 1);
						if (tmpcnt > cnt)
						{
							cnt = tmpcnt;
							x = i;
							y = j;
						}
						//����forѭ����ʱ��xy��¼����Χ����Ҫ������������Ǹ����ص������
					}
				//���cnt==-1˵������edge����false��Ҳ����˵����mymask[i��j]����1���ǲ���Ҫ��䣬����while
				if (cnt == -1) break;

				bool debug = false;
				bool debug2 = false;


				//�ⲿ���ٴα���ȫͼ���Ƚ�һ�������ں�����ͼƬ�����������Ƿ������ƵĿ�
				int k0 = min(x, bs), k1 = min(N - 1 - x, bs);
				int l0 = min(y, bs), l1 = min(M - 1 - y, bs);
				int sx, sy, min_diff = INT_MAX;	//����intֵ
				for (int i = step; i + step < N; i += step)
					for (int j = step; j + step < M; j += step)
					{
						//ͨ��usable�ҵ�����Ĳ���Ҫ�������ص�
						//���==2˵���������������
						if (usable[i][j] == 2)continue;
						int tmp_diff = 0;
						//ȡ��xy��Χstep�ľ�������
						for (int k = -k0; k <= k1; k++)
							for (int l = -l0; l <= l1; l++)
							{
								//printf("%d %d %d %d %d %d\n", i + k, j + l, x + k, y + l, N, M);
								//ij��ʾ���������ȽϵĲ���Ҫ�������������
								//xy��ʾ��ǰ��Ҫ�����ĵ㣬��֮ǰ��forѭ������
								//[x + k][y + l]��ʾxy��step�����ڵ�ĳ��
								//[i + k][j + l]��ʾij��step�����ڵ�ĳ��
								if (my_mask[x + k][y + l] != 0)
									tmp_diff += dist(result.at<Vec3b>(i + k, j + l), result.at<Vec3b>(x + k, y + l));
								//tmp_diff��������������Ӧ��֮�䣬RGBֵ�Ĳ��죻��Ȼ��Ҫȫͼ�����ҵ�һ����С��tmpdiff����˵����������������


								//--------------------------�����ƺ������£�����û�й涨�����Ӧ��ķ�Χ��û�п���������--------------------//



							}
						sum_diff[i][j] = tmp_diff;
						if (min_diff > tmp_diff)
						{
							sx = i;
							sy = j;
							min_diff = tmp_diff;
						}
						//����ѭ����ʱ�򣬵õ����ǶԱ�xy����Сtmpdiff�ĵ������sx��sy
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
				//�ã�sx��sy����Χ�ĵ��RGBֵ���xy��Χ��Ҫ�����ĵ�
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
	//�ĸ����룺mask��line��
	
	//����ԭͼ
	Mat3b img = imread("Origin.png");//5,1
	
	//�����ֵ����maskͼ��
	Mat1b mask = Mat::zeros(img.rows, img.cols, CV_8UC1);
	mask = imread("Mask.bmp", 0);
	
	threshold(mask, mask, 125, 255, CV_THRESH_BINARY_INV);
	imshow("img", img);
	waitKey(10);
	imshow("mask", mask);
	waitKey(10);
	//���ɴ���mask����û�н��в�ȫ��ͼ
	Mat3b result;
	result.zeros(img.size());
	img.copyTo(result, mask);
	imshow("result", result);
	waitKey(10);
	//����linemask
	Mat1b Linemask = Mat::zeros(img.rows, img.cols, CV_8UC1);
	//���ս������
	Mat3b finalResult(img.size());
	img.copyTo(finalResult);

	//TextureCompletion2(mask, Linemask, result, finalResult);

	system("pause");
}


