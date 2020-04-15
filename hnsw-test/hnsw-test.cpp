// hnsw-test.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <hnswlib/hnswlib.h>
#include <stdio.h>

using namespace cv;
using namespace hnswlib;

void rootSift(Mat &descriptors, const float eps = 1e-7)
{
	// Compute sums for L1 Norm
	Mat sums_vec;
	descriptors = abs(descriptors); //otherwise we draw sqrt of negative vals
	reduce(descriptors, sums_vec, 1 /*sum over columns*/, CV_REDUCE_SUM, CV_32FC1);
	for (int row = 0; row < descriptors.rows; row++) {
		int offset = row * descriptors.cols;
		for (int col = 0; col < descriptors.cols; col++) {
			descriptors.at<float>(offset + col) = sqrt(descriptors.at<float>(offset + col) /
				(sums_vec.at<float>(row) + eps) /*L1-Normalize*/);
		}
		// L2 distance
		normalize(descriptors.row(row), descriptors.row(row), 1.0, 0.0, NORM_L2);

	}
	return;
}

int main()
{
	Mat srcImg1 = imread("basketball1.png", 1);
	Mat srcImg2 = imread("basketball2.png", 1);
	if (srcImg1.empty() || srcImg2.empty())
		return -1;
// 	resize(srcImg1, srcImg1, Size(500, 500));
// 	resize(srcImg2, srcImg2, Size(500, 500));

	Ptr<cv::Feature2D> detector1 = xfeatures2d::SIFT::create(128);
	std::vector<KeyPoint> keypoints1, keypoints2;
	Mat descriptors1, descriptors2;
	// 提取图片的特征点及特征点描述
	detector1->detect(srcImg1, keypoints1);
	detector1->compute(srcImg1, keypoints1, descriptors1);
	rootSift(descriptors1);
	//drawKeypoints(srcImg1, keypoints1, srcImg1, Scalar::all(-1), DrawMatchesFlags::DRAW_OVER_OUTIMG);

	Ptr<cv::Feature2D> detector2 = xfeatures2d::SIFT::create(128);
	detector2->detect(srcImg2, keypoints2);
	detector2->compute(srcImg2, keypoints2, descriptors2);
	rootSift(descriptors2);
	//drawKeypoints(srcImg2, keypoints2, srcImg2, Scalar::all(-1), DrawMatchesFlags::DRAW_OVER_OUTIMG);

	// 创建hnsw
	size_t maxCount = descriptors1.rows;
	int M = 32;
	int efConstruction = 1000;// TestFindSelf，140以下有不能查找到自己的情况
	int vecdim = descriptors1.cols;// 特征维度
	hnswlib::InnerProductSpace *L2Space = new hnswlib::InnerProductSpace(vecdim);
	hnswlib::HierarchicalNSW<float> *hnsw = new hnswlib::HierarchicalNSW<float>(L2Space, maxCount, M, efConstruction);
	int count = 0;
	for (size_t i = 0; i < descriptors1.rows; i++)
	{
		float *mass = new float[vecdim];
		for (size_t j = 0; j < vecdim; j++)
		{
			mass[j] = descriptors1.at<float>(i, j);
		}
		hnsw->addPoint((void*)mass, count);
		++count;
		delete[] mass;
	}

	// 搜索
	int nn = 3;
	std::vector<DMatch> matchPoints;
	Mat indices, dists;
	for (int i = 0; i < descriptors2.rows; i++)
	{
		float *mass = new float[vecdim];
		for (int j = 0; j < vecdim; j++)
		{
			mass[j] = descriptors2.at<float>(i, j);
		}
		std::priority_queue<std::pair<float, hnswlib::labeltype>> candidates = hnsw->searchKnn((void*)mass, nn);
		delete[] mass;

		for (int k = 0; k < nn; k++)
		{
			double dAngle = abs(keypoints2.at(i).angle - keypoints1.at(candidates.top().second).angle);
			if (abs(candidates.top().first) < 0.1 && dAngle < 10)
			{
				std::cout << "angle diff: " << dAngle << ", dist: " << candidates.top().first << std::endl;
				DMatch dmatches(candidates.top().second, i, candidates.top().first);
				matchPoints.push_back(dmatches);
			}
			candidates.pop();
		}
	}

	Mat output;
	drawMatches(srcImg1, keypoints1, srcImg2, keypoints2, matchPoints, output);
	namedWindow("result", WINDOW_AUTOSIZE);
	imshow("result", output);
	waitKey(0);

	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
