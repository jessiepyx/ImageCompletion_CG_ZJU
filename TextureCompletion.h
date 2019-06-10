#ifndef IMAGECOMPLETION_TEXTURECOMPLETION_H
#define IMAGECOMPLETION_TEXTURECOMPLETION_H

#include <opencv2/opencv.hpp>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <string.h>
#include <queue>


using namespace cv;
using namespace std;

void texture(Mat origin, Mat img, Mat mask, Mat &result, Mat Linemask, string listpath);

#endif //IMAGECOMPLETION_TEXTURECOMPLETION_H
