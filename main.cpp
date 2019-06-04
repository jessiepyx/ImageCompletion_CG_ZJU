#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "StructurePropagation.h"
#include "OpenCvUtility.h"

using namespace std;
using namespace cv;

#define NUM_OF_IMAGES 7
#define USER_DRAW_MASK 0
#define PRE_MADE_MASK 1
#define LINE_STRUCTURE 0
#define CURVE_STRUCTURE 1

Mat img;
Mat mask;
Mat mask_inv;
Mat draw_mask;
Mat show_brush;
Mat img_masked;
Mat draw_structure;
Mat1b mask_structure;
Mat sp_result;
Mat ts_result;
Point pt;
Point prev_pt;
Point points[2] = {Point(-1, -1), Point(-1, -1)};
vector<Point> curvePoints;
vector<vector<Point>> PointsList;
vector<vector<Point>> mousepoints;
StructurePropagation SP;
int brush_size;
int img_current = 0;
int block_size = 10;
int sample_step = 2;
int line_or_curve = LINE_STRUCTURE;
int points_i=0;

void get_input_image();
void get_input_mask(int mask_from);
static void callback_draw_mask(int event, int x, int y, int flags, void* param);
void make_masked_image();
void show_interface();
static void callback_draw_structure(int event, int x, int y, int flags, void* param);

int main(int argc, char* argv[])
{
    get_input_image();
    get_input_mask(PRE_MADE_MASK);
    make_masked_image();
    show_interface();
	return 0;
}

/**
 * Choose one of the original images as input.
 * Press Key '[' to switch to last one, and Key ']' to switch to next one.
 * Press Key 'esc' to confirm the choice.
 */
void get_input_image()
{
    img = imread("img" + to_string(img_current) + ".png", 1);
    imshow("img", img);
    char k = waitKey(0);
    while (k != 27)
    {
        // last image
        if (k == '[')
        {
            img_current = (img_current + NUM_OF_IMAGES - 1) % NUM_OF_IMAGES;
        }
        // next image
        else if (k == ']')
        {
            img_current = (img_current + 1) % NUM_OF_IMAGES;
        }
        img = imread("img" + to_string(img_current) + ".png", 1);
        imshow("img", img);
        k = waitKey(0);
    }
    destroyAllWindows();
}

/**
 * Get a mask for the region of interest from the input image.
 * If user draws the mask,
 *      press Key '[' to get a smaller brush, and Key ']' to get a larger brush;
 *      press Key 'esc' to save the mask, and Key 'r' to reset.
 */
void get_input_mask(int mask_from)
{
    // user draws the mask
    if (mask_from == USER_DRAW_MASK)
    {
        mask = Mat::zeros(img.rows, img.cols, CV_8UC1);
        draw_mask = img.clone();
        show_brush = draw_mask.clone();
        brush_size = 30;
        prev_pt = Point(-1, -1);

        namedWindow("draw mask");
        imshow("draw mask", show_brush);
        setMouseCallback("draw mask", callback_draw_mask);

        char k = waitKey(0);
        while (k != 27)
        {
            // reset
            if (k == 'r')
            {
                mask = Mat::zeros(img.rows, img.cols, CV_8UC1);
                draw_mask = img.clone();
                prev_pt = Point(-1, -1);
            }
            // smaller brush
            else if (k == '[')
            {
                cout << "[" << endl;
                if (brush_size > 1)
                {
                    brush_size--;
                }
            }
            // larger brush
            else if (k == ']')
            {
                if (brush_size < 40)
                {
                    brush_size++;
                }
            }

            show_brush = draw_mask.clone();
            circle(show_brush, pt, brush_size, Scalar(255, 0, 255), -1);
            imshow("draw mask", show_brush);

            k = waitKey(0);
        }
    }
    // load pre-made mask
    else if (mask_from == PRE_MADE_MASK)
    {
        mask = imread("mask" + to_string(img_current) + ".png", 0);
    }
    threshold(mask, mask_inv, 100, 255, CV_THRESH_BINARY_INV);
    destroyAllWindows();
}

/**
 * Mouse callback function for drawing the mask.
 */
static void callback_draw_mask(int event, int x, int y, int flags, void* param)
{
    pt = Point(x, y);
    if ((event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON)) || event == CV_EVENT_LBUTTONDOWN)
    {
        if (prev_pt.x == -1)
        {
            prev_pt = pt;
        }
        line(mask, prev_pt, pt, Scalar(255), 1.5 * brush_size);
        line(draw_mask, prev_pt, pt, Scalar(255, 0, 0), 1.5 * brush_size);
        prev_pt = pt;
    }
    else if (event == CV_EVENT_LBUTTONUP)
    {
        prev_pt = Point(-1, -1);
    }

    show_brush = draw_mask.clone();
    circle(show_brush, pt, brush_size, Scalar(255, 0, 255), -1);
    imshow("draw mask", show_brush);
}

/**
 * Mask the region of interest.
 */
void make_masked_image()
{
	img_masked = Mat::zeros(img.size(), CV_8UC3);
	img.copyTo(img_masked, mask_inv);
	imwrite("img_masked" + to_string(img_current) + ".png", img_masked);
}

/**
 * Show the user interface for drawing structural lines and curves.
 * Press Key 's' to run structure propagation, and Key 't' to run texture synthesis.
 * Press Key 'r' to reset, and Key 'a' to save.
 * Press Key 'e' to show curve points.
 */
void show_interface()
{
	prev_pt = Point(-1, -1);
	sp_result = img_masked.clone();
	draw_structure = img_masked.clone();

	namedWindow("run");
	createTrackbar("Block Size", "run", &block_size, 50);
	createTrackbar("Sample Step", "run", &sample_step, 20);
	createTrackbar("Line or Curve", "run", &line_or_curve, 1);
	imshow("run", draw_structure);
	setMouseCallback("run", callback_draw_structure);

	char k = waitKey(0);
	while (k != 27)
	{
		// structure propagation
		if (k == 's')
		{
            if (line_or_curve == CURVE_STRUCTURE)
            {
                PointsList.resize(mousepoints.size());
                for (int i = 0; i < mousepoints.size(); i++)
                {
                    Wang_GetCurve(mousepoints[i], PointsList[i]);
                }
            }

            for (int i = 0; i < PointsList.size(); i++)
            {
                DrawPoints(PointsList[i], draw_structure, CV_RGB(255, 0, 0), 1);
            }

            mask_structure = Mat::zeros(img.rows, img.cols, CV_8UC1);

            // run structure porpagation
            SP.SetParm(block_size, sample_step, line_or_curve);
            SP.Run(mask_inv, img_masked, mask_structure, PointsList, sp_result);

            draw_structure = sp_result.clone();
            imshow("run", draw_structure);

            PointsList.clear();
            mousepoints.clear();
		}
		// reset
		else if (k == 'r')
        {
		    draw_structure = img_masked.clone();
		    sp_result = img_masked.clone();
            mask_structure = Mat::zeros(img.rows, img.cols, CV_8UC1);
            imshow("run", draw_structure);

            PointsList.clear();
            mousepoints.clear();
        }
		// save
		else if (k == 'a')
        {
            imwrite("sp_result/sp" + to_string(img_current) + ".png", sp_result);
            imwrite("ts_result/ts" + to_string(img_current) + ".png", sp_result);
            imwrite("mask_structure/mask_s" + to_string(img_current) + ".bmp", mask_structure);
        }
		// show curve points
        else if (k == 'e')
        {
            PointsList.resize(mousepoints.size());
            for (int i = 0; i < mousepoints.size(); i++)
            {
                Wang_GetCurve(mousepoints[i], PointsList[i]);
                DrawPoints(PointsList[i], draw_structure, CV_RGB(0, 0, 255), 1);
            }
            imshow("run", draw_structure);
        }
        // texture synthesis
        else if (k == 't')
        {
            SP.TextureCompletion2(mask, mask_structure, sp_result, ts_result);
            imshow("run", ts_result);
        }
        k = waitKey(0);
	}
}

/**
 * Mouse callback function for drawing the structure lines/curves.
 */
static void callback_draw_structure(int event, int x, int y, int flags, void* param)
{
    if (line_or_curve == LINE_STRUCTURE)
    {
        if (event != CV_EVENT_LBUTTONDOWN)
        {
            return;
        }
        points[points_i].x = x;
        points[points_i].y = y;
        points_i = (points_i + 1) % 2;

        if (points[0].x != -1 && points[1].x != -1 && points_i == 0)
        {
            vector<Point> line;
            LineInterpolation(points, line);
            PointsList.push_back(line);

            DrawPoints(line, draw_structure, Scalar(255, 0, 255), 1);
            circle(draw_structure, points[0], 3, Scalar(255,0,0), CV_FILLED);
            circle(draw_structure, points[1], 3, Scalar(255,0,0), CV_FILLED);

            imshow("run", draw_structure);
        }
    }
    else
    {
        if (event == CV_EVENT_LBUTTONUP)
        {
            prev_pt = Point(-1, -1);
            mousepoints.push_back(curvePoints);
            curvePoints = vector<Point>();
        }
        if ( event == CV_EVENT_LBUTTONDOWN )
        {
            prev_pt = Point(x, y);
            curvePoints.push_back(prev_pt);
        }
        else if (event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON))
        {
            pt = Point(x,y);
            curvePoints.push_back(pt);

            if (prev_pt.x < 0)
            {
                prev_pt = pt;
            }
            line(draw_structure, prev_pt, pt, cvScalarAll(255), 1, 8, 0);
            prev_pt = pt;
            imshow("run", draw_structure);
        }
    }
}
