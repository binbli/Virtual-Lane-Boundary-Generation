#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<iostream>
#include<vector>

using namespace cv;
using namespace std;

string binaryWindow = "window", originalWindow = "or", convertedWindow = "co";
pair<int, int> minMax1, minMax2, minMax3;
int channelNum = 0;
vector<int> codes = {-1,
					 COLOR_BGR2GRAY,
					 COLOR_BGR2Lab,
					 COLOR_BGR2YCrCb,
					 COLOR_BGR2HSV,
					 COLOR_BGR2XYZ,
					 COLOR_BGR2Luv,
					 COLOR_BGR2HLS,
//					  COLOR_BGR2BGR565, 
//					  COLOR_BGR2BGR555, 
					 COLOR_BGR2HSV_FULL,
					 COLOR_BGR2HLS_FULL,
					 COLOR_BGR2YUV,
//					  COLOR_BGR2YUV_I420,
//					  COLOR_BGR2YUV_YV12  
};
Mat in, out;
int code = codes[0];

void Detect() {
	Mat grayImage;
	if (code == -1)
		out = in.clone();
	else
		cvtColor(in, out, code);
	if (out.channels() == 1)
		inRange(out, Scalar(minMax1.first), Scalar(minMax1.second), grayImage);
	else
		inRange(out, Scalar(minMax1.first, minMax2.first, minMax3.first),
				Scalar(minMax1.second, minMax2.second, minMax3.second), grayImage);
	imshow(originalWindow, in);
	vector<Mat> channels;
	split(out, channels);
	imshow(convertedWindow, channels[channelNum]);
	imshow(binaryWindow, grayImage);
}

void on_low1_thresh_trackbar(int v, void *ptr) {
	setTrackbarPos("Low 1", binaryWindow, v);
	Detect();
}

void on_high1_thresh_trackbar(int v, void *ptr) {
	setTrackbarPos("High 1", binaryWindow, v);
	Detect();
}

void on_low2_thresh_trackbar(int v, void *ptr) {

	setTrackbarPos("Low 2", binaryWindow, v);
	Detect();
}

void on_high2_thresh_trackbar(int v, void *ptr) {
	setTrackbarPos("High 2", binaryWindow, v);
	Detect();
}

void on_low3_thresh_trackbar(int v, void *ptr) {
	setTrackbarPos("Low 3", binaryWindow, v);
	Detect();
}

void on_high3_thresh_trackbar(int v, void *ptr) {
	setTrackbarPos("High 3", binaryWindow, v);
	Detect();
}

void on_conversion_change(int v, void *ptr) {
	setTrackbarPos("Conversion code", convertedWindow, v);
	code = codes[v];
	Detect();
}

void on_channel_num(int v, void *ptr) {
	setTrackbarPos("Channel number", convertedWindow, v);
	Detect();
}

void CreateTrackbars() {
	createTrackbar("Channel number", convertedWindow, &channelNum, 2, on_channel_num);
	createTrackbar("Conversion code", convertedWindow, &code, codes.size() - 1, on_conversion_change);
	createTrackbar("Low 1", binaryWindow, &minMax1.first, 255, on_low1_thresh_trackbar);
	createTrackbar("High 1", binaryWindow, &minMax1.second, 255, on_high1_thresh_trackbar);
	createTrackbar("Low 2", binaryWindow, &minMax2.first, 255, on_low2_thresh_trackbar);
	createTrackbar("High 2", binaryWindow, &minMax2.second, 255, on_high2_thresh_trackbar);
	createTrackbar("Low 3", binaryWindow, &minMax3.first, 255, on_low3_thresh_trackbar);
	createTrackbar("High 3", binaryWindow, &minMax3.second, 255, on_high3_thresh_trackbar);
}

int main(int argc, char *argv[]) {
	in = imread(argv[1]);
	namedWindow(binaryWindow, WINDOW_AUTOSIZE);
	namedWindow(originalWindow, WINDOW_AUTOSIZE);
	namedWindow(convertedWindow, WINDOW_AUTOSIZE);
	CreateTrackbars();
	Detect();
	waitKey(0);
}
