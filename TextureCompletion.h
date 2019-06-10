//
// Created by Alexis Peng on 2019-06-10.
//

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

void texture(Mat3b origin, Mat3b img, Mat1b mask, Mat &result, Mat1b Linemask, string listpath);
#endif //IMAGECOMPLETION_TEXTURECOMPLETION_H
