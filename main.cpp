#include <opencv.hpp>
#include "poissonDisc.h"
#include "delTriangulation.h"

#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


using namespace std;
using namespace cv;

int loadimage(const char* path, int *w, int* h, char* data) {
	

	return 0;
}

int poissonTriangular() {
	//const int width = 800;
	//const int height = 800;
	
	const char* path = "C:\\Users\\a\\Downloads\\src_1643_1.png";
	int width, height, channel;
	unsigned char* data = stbi_load(path, &width, &height, &channel, 4);
	//unsigned char* data
	Mat image(height, width, CV_8UC4, data);
	cvtColor(image, image, COLOR_BGRA2RGBA);
	DelTriangulater* trianguler = new DelTriangulater();
	trianguler->doinit(width, height);
	vector<MPOINTF> points;
	vector<MUInt32> vertexIdx;
	double mindist = 35;
	vector<P> poissonPts = sampling(double(width), double(height), mindist, data);

	for (int i = 0; i < poissonPts.size(); i++)
	{
		Point p(poissonPts[i].first, poissonPts[i].second);
		MPOINTF v = { poissonPts[i].first, poissonPts[i].second };
		points.push_back(v);
		circle(image, p, 3, Scalar(0, 0, 255));
	}

	trianguler->setTriangulePoints(points);
	trianguler->triVertexToIndex(points, vertexIdx);
	int triNum = trianguler->getTriangleNum();

	vector<Point> cvpt(3);
	for (size_t num = 0; num < triNum * 3; num++) {
		MUInt32 idx = vertexIdx[num++];
		cvpt[0] = cv::Point(points[idx].x, points[idx].y);
		idx = vertexIdx[num++];
		cvpt[1] = cv::Point(points[idx].x, points[idx].y);
		idx = vertexIdx[num];
		cvpt[2] = cv::Point(points[idx].x, points[idx].y);
		line(image, cvpt[0], cvpt[1], cv::Scalar(255, 0, 255, 0.0), 1, cv::LINE_AA, 0);
		line(image, cvpt[1], cvpt[2], cv::Scalar(255, 0, 255, 0.0), 1, cv::LINE_AA, 0);
		line(image, cvpt[2], cvpt[0], cv::Scalar(255, 255, 255, 0.0), 1, cv::LINE_AA, 0);
	}

	imshow("poisson disc", image);
	waitKey(0);
	return 0;
}

Point newPoint(int x, int y, int width, int height) {
	return Point{ x * width / 10000 , y * height / 10000 };
}

Point mid(Point a, Point b) {
	return Point{ (a.x + b.x) / 2, (a.y + b.y) / 2 };
}

#define IMAGE_2

int main() {
	vector<Point> points;

#ifdef IMAGE_1
	const char* path = "G:\\test\\image\\1.jpg";
	int width, height, channel;
	unsigned char* data = stbi_load(path, &width, &height, &channel, 4);
	Mat image(height, width, CV_8UC4, data);
	cvtColor(image, image, COLOR_BGRA2RGBA);
	//0
	points.push_back(newPoint(4781, 6000, width, height));
	//1
	points.push_back(newPoint(4937, 10000, width, height));
	//2
	points.push_back(newPoint(0, 0, width, height));
	//5
	points.push_back(newPoint(8750, 9968, width, height));

	points.push_back(mid(points[0], points[1]));
	points.push_back(mid(points[2], points[3]));
#endif // IMAGE_1

#ifdef IMAGE_2
	const char* path = "G:\\test\\image\\2.jpg";
	int width, height, channel;
	unsigned char* data = stbi_load(path, &width, &height, &channel, 4);
	Mat image(height, width, CV_8UC4, data);
	cvtColor(image, image, COLOR_BGRA2RGBA);
	//0
	points.push_back(newPoint(6525, 3875, width, height));
	//1
	points.push_back(newPoint(4406, 7750, width, height));
	//2
	points.push_back(newPoint(1313, 7750, width, height));
	//5
	points.push_back(newPoint(7457, 7750, width, height));

	points.push_back(mid(points[0], points[1]));
	points.push_back(mid(points[2], points[3]));
#endif // IMAGE_2

#ifdef IMAGE_3
	const char* path = "G:\\test\\image\\3.png";
	int width, height, channel;
	unsigned char* data = stbi_load(path, &width, &height, &channel, 4);
	Mat image(height, width, CV_8UC4, data);
	cvtColor(image, image, COLOR_BGRA2RGBA);
	//0
	points.push_back(newPoint(4375, 4156, width, height));
	//1
	points.push_back(newPoint(5333, 7625, width, height));
	//2
	points.push_back(newPoint(1625, 7625, width, height));
	//5
	points.push_back(newPoint(9000, 7625, width, height));

#endif // IMAGE_3

	points.push_back(mid(points[0], points[1]));

	Point hvector = Point((points[3].x - points[2].x), (points[3].y - points[2].y));
	//double hlength = norm(hvector);
	Point vvector = Point(hvector.y, -hvector.x);
	Point vNormal = vvector / norm(vvector);

	double length = norm(points[0] - points[1]);

	Point new2 = points[2] + 0.25 * length * vNormal;
	Point new5 = points[3] + 0.25 * length * vNormal;
	points.push_back(new2);
	points.push_back(new5);

	for (size_t i = 0; i < points.size(); i++)
	{
		cv::circle(image, points[i], 3, Scalar(0, 0, 255));
	}

	imshow("skeleton", image);
	waitKey(0);
	return 0;
}