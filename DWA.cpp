 #include <math.h>
#include <iostream>
#include <utility>
#include <vector>
#include <time.h>
#pragma warning (disable : 4996)
#include <opencv2\opencv.hpp>
#include <stdio.h>     
#include <stdlib.h>
#include <string>
#include "DWA.h"
#include "rs232.h"
#include "laser_filters.h"
#include "LineExtractionFromLaser.h"

using namespace std;
using namespace cv;

#define M_PI          3.1415927F
#define DEG_TO_RAD(deg) ((deg) / 180.0F * M_PI)

#define max_velocity 0.3F//  0.3m/s
#define max_omega 1.346F
#define car_width 0.4F
#define alpha 0.1F
#define beta 2.0F
#define obs_number 10.0F
#define DRAW 
//#define DEBUG



// -----------------------
// dynamic window approach 
// -----------------------
int DWA::cal_DWA(arima_laserscan scans){
	// VMap為畫圖所用
	VMap VelocityMap;
	std::vector<vector<slope_range>> slope_vector;
	//------cluster-------
	cluster_vector clusters = cluster(scans);
	//--------------------------------------------------------------------------------
	// 為了之後判斷laser point有沒有介於兩輪軌跡之間，先計算laser point相對於左右輪座標的radius
	//--------------------------------------------------------------------------------
		for (unsigned i = 0; i < clusters.size(); ++i){
		vector<slope_range> temp;
		for (unsigned j = 0; j < clusters[i].size(); ++j){
			float x = (float)(scans.ranges[clusters[i][j]] * cos((1.5*M_PI) - DEG_TO_RAD(scans.angles[clusters[i][j]])));			
			float y = (float)(scans.ranges[clusters[i][j]] * sin((1.5*M_PI) - DEG_TO_RAD(scans.angles[clusters[i][j]])));
			// 平移半個車身
			float r_left = (((x + car_width) * (x + car_width)) + (y * y)) / (2 * (x + car_width));
			float r_right = (((x - car_width) * (x - car_width)) + (y * y)) / (2 * (x - car_width));
			float r = ((x  * x) + (y * y)) / (2 * x);
			slope_range sr;
			sr.r_left = r_left;
			sr.r_right = r_right;
			sr.r = r;
			sr.y = y;
			temp.push_back(sr);
		}
		slope_vector.push_back(temp);
	}
	//find the first and last point in a cluster
	/*
		//------------------------------------------------------
		//		 y                          y
		//       ^                          ^
		//		 |                          |
		//90---------->x    ===>    180---------->0 x
		//       |                          |
		//       0  laser domain            |  robot facing domain
		//-------------------------------------------------------
	*/
	//-----weight of every curve------
	std::vector<float> weight_vector;
	std::vector<float> clearances;
	//--------------------------------
	std::vector<std::pair<float, float>> tempv;
	// divide omega into 101 intervals including 0 , with Velocity = 0.3m/s
	// from -1.346 ~ 0 ~ 1.346 
	for (int n = 0; n < 101; ++n){
		float w;
		if (n < 50){
			w = -1 * max_omega / 50 * (50 - n);
		}
		else if (n == 50){
			w = 0;
		}
		else{
			w = max_omega / 50 * (n - 50);
		}
		std::pair<float, float> vw(max_velocity, w);
		tempv.push_back(vw);
		float slope, slope_right, slope_left;
		// 軌跡的半徑為速度與角速度之比值
		slope = max_velocity / w; //std::cout <<"slope:"<< slope << " ";
		// 左右輪的軌跡半徑
		slope_right = slope - car_width;
		slope_left = slope + car_width;
		//std::cout << "srv.size:" << slope_range_vector.size() << std::endl;
		//------------------------------------------------------------------------------
		// 如果障礙物出現在左右輪軌跡之間，則計算從原點到障礙物沿著軌跡的弧長當作clearance
		//------------------------------------------------------------------------------
		// clearance = 1000 if no obstacles ahead
		float clearance = 1000;
		float r, y1;
		for (unsigned i = 0; i < slope_vector.size(); ++i){
			for (unsigned j = 1; j < slope_vector[i].size(); ++j){
				r = slope_vector[i][j].r;
				y1 = slope_vector[i][j].y;
				float A = slope_vector[i][j].r_right - slope_right;
				float B = slope_vector[i][j].r_left - slope_left;
				float C = slope_vector[i][j].r_right;
				float D = slope_vector[i][j].r_left;
				//if ((((A > 0 || C < 0) && (B < 0 && D > 0)) || ((B < 0 || D > 0) && (A > 0 && C < 0)) || ((slope == 0) && (D > 0 && C < 0))) && (clearance > r * asin(y1 / r))) //<-這理這樣應該沒錯吧 ? ?
				//{
				//	// iterate through all points,找出最小的clearance
				//	clearance = (float)(r * asin(y1 / r));
				//}

				if (n == 50){
					if (slope_vector[i][j].r_left > 0 && slope_vector[i][j].r_right < 0){
						if (clearance > y1)
							clearance = y1;
					}
				}
				else if (slope > 0){
					if ((slope_vector[i][j].r_right > slope_right || slope_vector[i][j].r_right < 0) && (slope_vector[i][j].r_left < slope_left && slope_vector[i][j].r_left>0)){
						// iterate through all points,找出最小的clearance
						if (clearance > r * asin(y1 / r))
							clearance = r * asin(y1 / r);
					}
				}
				else if (slope < 0){
					if ((slope_vector[i][j].r_left < slope_left || slope_vector[i][j].r_left > 0) && (slope_vector[i][j].r_right > slope_right && slope_vector[i][j].r_right < 0)){
						if (clearance > r * asin(y1 / r))
							clearance = r * asin(y1 / r);
					}
				}
			}
		}
		clearances.push_back(clearance);
		weight_vector.push_back(0);
	}
	VelocityMap.push_back(tempv);
	//-------------------------------
	// Normalize clearance
	//-------------------------------
	// 找出最大最小clearance
	float min_c = 1000, max_c = 0;
	for (int i = 0; i < (int)clearances.size(); i++){
		//fprintf(fptr, "clear %f\n", clearances[i]);
		if (clearances[i] < min_c)
			min_c = clearances[i];
		if (clearances[i] > max_c && clearances[i] != 1000)
			max_c = clearances[i];
	}
	//---------------------------------------------------------------------
	// 1.
	// 當count_obs<=10時,不需要太大的轉彎角度,因此w較大的軌跡不在選擇範圍之內
	// (因為只有前方120度的FOV,太大的轉彎角度永遠不會碰到障礙物,clearances會被當作1.03)
	// 2.
	// 如果軌跡上沒有障礙物,將權重設為1.02;如果障礙物太遠,將權重設為1.02
	// 3.
	// normalize clearances to 0.0~1.0
	//---------------------------------------------------------------------
	for (int j = 0; j < (int)clearances.size(); j++){
		if (count_obs<obs_number && (j <= 30 || j >= 70))
			clearances[j] = 0.0F;
		else if (clearances[j] == 1000)
			clearances[j] = 1.02F;
		else if (clearances[j] > 2.5F)
			clearances[j] = 1.02F;
		//else if (clearances[j] == -1000 )
		//clearances[j] = 0.0F;
		else
			clearances[j] = (clearances[j] - min_c) / (max_c - min_c);
	}
	//-------------------------------------------------------------
	// Target heading
	//-------------------------------------------------------------
	//-------way point direction-------
	int wpd = head_direction;
	//---------------------------------
	float w;//角速度
	float next_direction;//下一秒機器人頭向
	float sita;//頭向與目標方向夾角
	std::vector<float> target;
	for (int i = 0; i <= 100; ++i){
		if (i < 50){
			w = max_omega / 50 * (50 - i);
		}
		else if (i == 50){
			w = 0;
		}
		else{
			w = -1 * max_omega / 50 * (i - 50);
		}
		// 計算經過一秒移動後機器人頭向的角度
		next_direction = robot_head + ((w*1.0F) / M_PI*180.0F);
		if (next_direction >= 360){
			next_direction = next_direction - 360;
		}
		else if (next_direction < 0){
			next_direction = next_direction + 360;
		}
		// 計算機器人頭向與目標方向的夾角
		sita = abs(wpd - next_direction);
		// 夾角介於0~180
		if (sita >= 180)
			sita = 360 - sita;
		// 夾角越小權重越高
		sita = 180 - sita;
		//--------------------------
		target.push_back(sita);		
	}
	//----------------------------------------------------
	// Normalize target
	//----------------------------------------------------
	min_c = 1000;
	max_c = 0;
	for (int i = 0; i < (int)target.size(); i++){
		if (target[i] < min_c)
			min_c = target[i];
		if (target[i] > max_c)
			max_c = target[i];
	}

	for (int j = 0; j < (int)target.size(); j++){
		target[j] = (target[j] - min_c) / (max_c - min_c);
	}		
	int index = 0;// returned index
	float temp_weight = -999;
	//---------------------------------------------- 
	// objective function
	//----------------------------------------------
	if (count_obs <= obs_number){
		std::fprintf(fptr, "target\n");
		std::printf("target\n");
	}
	for (int i = 0; i < (int)weight_vector.size(); ++i){
		float weight_target;
		//-------------------------------------------------
		// 如果障礙物離很近,啟動純避障模式,不把target考慮進來
		//-------------------------------------------------
		if (count_obs > obs_number){
			weight_target = mask[i]*0.01;
		}
		else{
			// 將target乘上0.2-1.0-0.2的mask
			weight_target = target[i] * mask[i];
		}
		weight_vector[i] = weight_vector[i]
			+ weight_target * alpha
			+ clearances[i] * beta;
		weight_vector[i] *= 0.3F;
		//---- find the heaviest weight--------
		if (weight_vector[i] > temp_weight){
			temp_weight = weight_vector[i];
			index = i;
		}
		std::fprintf(fptr, "index: %i h: %f c: %f \n", i, weight_target * alpha, clearances[i] * beta);
	}
	//------manual control---------
	if (head_direction == 0){
		index = 100;
	}
	else if (head_direction == 180){
		index = 0;
	}
	//-----if ( obstacles too close )stop-----
	else if (count_backup > 20)
		Lock = true;
	////----------------------------------------
	//// CLineExtractionFromLaser 
	////----------------------------------------
	//LnExt::CLineExtractionFromLaser lineExtractor;
	//LnExt::ArimaLaserscan laserScans;
	//LnExt::ClusterVector CVcluster;
	//laserScans.angles = scans.angles;
	//laserScans.ranges = scans.ranges;

	//int bestFrom = -1;
	//int bestTo = -1;
	//bool hasValidLineSegment = lineExtractor.LineExtract(laserScans, clusters, bestFrom, bestTo);
	//// ------------------
	//// heading correction
	//// ------------------
	//if (hasValidLineSegment)
	//{
	//	float _beta = 0; // beta is the angle between the global 0-deg direction and the line landmark

	//	// swap "bestFrom" and "bestTo"      
	//	if (laserScans.angles[bestFrom] > laserScans.angles[bestTo])
	//	{
	//		int temp = bestFrom;
	//		bestFrom = bestTo;
	//		bestTo = temp;
	//	}
	//	// compute the including angle between robot's heading and the line
	//	cv::Point2f vecRobotHeading(-1, 0);
	//	cv::Point2f vecFromTo;
	//	vecFromTo.x = laserScans.ranges[bestTo] * cos((float)DEG_TO_RAD(laserScans.angles[bestTo])) - laserScans.ranges[bestFrom] * cos((float)DEG_TO_RAD(laserScans.angles[bestFrom]));
	//	vecFromTo.y = laserScans.ranges[bestTo] * sin((float)DEG_TO_RAD(laserScans.angles[bestTo])) - laserScans.ranges[bestFrom] * sin((float)DEG_TO_RAD(laserScans.angles[bestFrom]));
	//	float distFromTo = sqrt(vecFromTo.x * vecFromTo.x + vecFromTo.y * vecFromTo.y);
	//	float theta = 0;
	//	if (distFromTo != 0) theta = acos(-vecFromTo.x / distFromTo);
	//	// check if the robot is leaving or approaching the line
	//	bool approaching = true; // true if the robot is aproaching or moving parallely to the line
	//	if (vecFromTo.y != 0)
	//	{
	//		// y = mx + b, because m is not infinity
	//		float m = vecFromTo.y / vecFromTo.x;
	//		float xFrom = laserScans.ranges[bestFrom] * cos((float)DEG_TO_RAD(laserScans.angles[bestFrom]));
	//		float yFrom = laserScans.ranges[bestFrom] * sin((float)DEG_TO_RAD(laserScans.angles[bestFrom]));
	//		float b = yFrom - m * xFrom;
	//		if (-b / m > 0) approaching = false;
	//	}
	//	// debug
	//	if (approaching) printf("approaching\n");
	//	else printf("leaving\n");
	//	// robot's heading correction
	//	float robotHeading = (approaching) ? (beta - theta) : (beta + theta);
	//	printf("=====================\n");
	//	printf("Robot's heading is %f\n", robotHeading);
	//	printf("=====================\n");
	//}

#ifdef DRAW
		cv::Mat background(1000, 1000, CV_8UC3, cv::Scalar(255, 255, 255));
		int color = 240 / ((int)clusters.size() + 1);
		for (unsigned i = 0; i < clusters.size(); ++i){
			for (unsigned j = 0; j < clusters[i].size(); ++j){
				int idx = clusters[i][j];
				float x1 = (float)((scans.ranges[idx] * 100 * cos((1.5*M_PI) - DEG_TO_RAD(scans.angles[idx]))) + 500);
				float y1 = (float)((-1)*scans.ranges[idx] * 100 * sin((1.5*M_PI) - DEG_TO_RAD(scans.angles[idx])) + 500); //cout << scans.ranges[i] << " y1 "; cout << endl;
				cv::Point2f pt(x1, y1) ;
				cv::Point2f ptfont(x1+5, y1+5);
				cv::circle(background, pt, 1, cv::Scalar(i*color, 125, 255 - (i*color)), 2, 8, 0);
				if (j == 0 || j == clusters[i].size()-1)
				cv::putText(background,std::to_string(i), ptfont, FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(i*color, 125, 255 - (i*color)), 1, 8, false);

			}
		}
		//----------------------------
		//  CLineExtractionFromLaser
		//----------------------------
		/*if (hasValidLineSegment){
			float Fromx = (float)((scans.ranges[bestFrom] * 100 * cos((1.5*M_PI) - DEG_TO_RAD(scans.angles[bestFrom]))) + 500);
			float Fromy = (float)((-1)*scans.ranges[bestFrom] * 100 * sin((1.5*M_PI) - DEG_TO_RAD(scans.angles[bestFrom])) + 500);
			float Tox = (float)((scans.ranges[bestTo] * 100 * cos((1.5*M_PI) - DEG_TO_RAD(scans.angles[bestTo]))) + 500);
			float Toy = (float)((-1)*scans.ranges[bestTo] * 100 * sin((1.5*M_PI) - DEG_TO_RAD(scans.angles[bestTo])) + 500);
			cv::line(background, cv::Point2f(Fromx, Fromy), cv::Point2f(Tox, Toy), cv::Scalar(0, 255, 0), 3, 8, 0);
		}*/
		//----------------------------
		cv::Point2f pto = cv::Point2f(500, 500);
		for (int m = 0; m < 1; m++){
			for (int n = 0; n < 101; n++){
				float r = VelocityMap[m][n].first / VelocityMap[m][n].second;
				r = r * 100;
				float rout = r + (car_width * 100.0F);
				float rin = r - (car_width * 100.0F);
				if (n == 50 && VelocityMap[m][n].first != 0){
					if (n == index)
						cv::line(background, pto, cv::Point2f(500, 0), cv::Scalar(0, 0, 255), 3, 8, 0);
					else
						cv::line(background, pto, cv::Point2f(500, 0), cv::Scalar(255, 0, 0), 1, 8, 0);
					continue;
				}
				cv::Point2f prev = pto,prevout = pto, previn = pto;
				// draw 600 cycles move curve
				for (int k = 0; k < 600; k++){
					cv::Point2f pt, ptout, ptin;
					//cv::Point2f ptnew;
					if (VelocityMap[m][n].second == 0){
						break;
					}
					else {
						pt = cv::Point2f((float)(pto.x + r - r*cos(VelocityMap[m][n].second*0.1*k)), (float)(-1 * r*sin(VelocityMap[m][n].second*0.1*k) + pto.y));
						ptout = cv::Point2f((float)(pto.x - (car_width * 100) + rout - rout*cos(VelocityMap[m][n].second*0.1*k)), (float)(-1 * rout*sin(VelocityMap[m][n].second*0.1*k) + pto.y));
						ptin = cv::Point2f((float)(pto.x + (car_width * 100) + rin - rin*cos(VelocityMap[m][n].second*0.1*k)), (float)(-1 * rin*sin(VelocityMap[m][n].second*0.1*k) + pto.y));
						
					}
					if (n == index){
						cv::line(background, prev, pt, cv::Scalar(0, 0, 255), 3, 8, 0);
						cv::line(background, prevout, ptout, cv::Scalar(0, 255, 0), 3, 8, 0);
						cv::line(background, previn, ptin, cv::Scalar(255, 0, 0), 3, 8, 0);
					}
					else{
						cv::line(background, prev, pt, cv::Scalar(255, 0, 0), 1, 8, 0);
						//cv::line(background, prevout, ptout, cv::Scalar(200, 105, 0), 1, 8, 0);
						//cv::line(background, previn, ptin, cv::Scalar(0, 105, 200), 1, 8, 0);
					}
					prev = pt;
					prevout = ptout;
					previn = ptin;
				}
			}
		}
		// save every frame
		String image_name = "Saved_image/image.jpg";
		char buffer[10];
		image_name.insert(17, itoa(count_frame, buffer, 10));
		cv::imwrite(image_name, background);
		count_frame++;
		// ----------------
		cv::imshow("Display window", background);
		cv::waitKey(1);
#endif

	return index;
}

// --------------------------------
// cluster function for laser scans
// --------------------------------
DWA::cluster_vector DWA::cluster(arima_laserscan& scans){

	std::vector<std::vector<int> > clusters;
	clusters.clear();
	//-------------------------------------------------
	// 計算每個點與下一個點的距離，如果dist<=0.5m則跳過此點
	// 當dist>0.5時將前面所有跳過的點歸為同一個cluster
	// 當i == scans.angles.size()-1時，將之與i==0比較
	// 因為角度為一圈的關係，有可能第一個與最後一個點為同一cluster
	//------------------------------------------------
	int m = 0;
	for (int i = 0; i <= (int)scans.angles.size() - 1; i++)
	{
		float dist;
		if (i == (int)scans.angles.size() - 1){
			//dist = (float)sqrt(pow(scans.ranges[i] * cos(DEG_TO_RAD(scans.angles[i])) - scans.ranges[0] * cos(DEG_TO_RAD(scans.angles[0])), 2)
			//	+ pow(scans.ranges[i] * sin(DEG_TO_RAD(scans.angles[i])) - scans.ranges[0] * sin(DEG_TO_RAD(scans.angles[0])), 2));
			 m = -1; 
		}
		dist = (float)sqrt(pow(scans.ranges[i] * cos(DEG_TO_RAD(scans.angles[i])) - scans.ranges[m + 1] * cos(DEG_TO_RAD(scans.angles[m + 1])), 2)
					+ pow(scans.ranges[i] * sin(DEG_TO_RAD(scans.angles[i])) - scans.ranges[m + 1] * sin(DEG_TO_RAD(scans.angles[m + 1])), 2));
		m++;
		//else{
		//	dist = (float)sqrt(pow(scans.ranges[i] * cos(DEG_TO_RAD(scans.angles[i])) - scans.ranges[i + 1] * cos(DEG_TO_RAD(scans.angles[i + 1])), 2)
		//		+ pow(scans.ranges[i] * sin(DEG_TO_RAD(scans.angles[i])) - scans.ranges[i + 1] * sin(DEG_TO_RAD(scans.angles[i + 1])), 2));
		//}
		
		std::vector<int> new_cluster;
		if ((dist > 0.5) || (i == (int)scans.angles.size() - 1)) // threshold according to [Using Boosted Features for the Detection of People in 2D Range Data] 
		{
			int j = 0;
			if (clusters.size() != 0)
			{
				int cluster_index = (int)clusters.size() - 1;
				int scan_index = (int)clusters[cluster_index].size() - 1;
				//if (i == (int)scans.angles.size() - 1) scan_index = 0;
				j = clusters[cluster_index][scan_index] + 1;
				if ((cluster_index >= 0) && (i == (int)scans.angles.size() - 1))
				{
					for (int k = 0; k < clusters[0].size(); k++)
					{
						new_cluster.push_back(k);
					}
					clusters.erase(clusters.begin());
				}
			}
			for (; j <= i; j++)
			{
				new_cluster.push_back(j);
			}
			clusters.push_back(new_cluster);
		}
		//if (dist > 0.5) // threshold according to [Using Boosted Features for the Detection of People in 2D Range Data] 
		//{
		//	std::vector<int> new_cluster;
		//	int j = 0;
		//	if (clusters.size() != 0)
		//	{
		//		int cluster_index = (int)clusters.size() - 1;
		//		int scan_index = (int)clusters[cluster_index].size() - 1;
		//		j = clusters[cluster_index][scan_index] + 1;
		//	}
		//	for (; j <= i; j++)
		//	{
		//		new_cluster.push_back(j);
		//	}
		//	clusters.push_back(new_cluster);
		//}
		//else if (i == (int)scans.angles.size() - 1)
		//{
		//	std::vector<int> last_cluster;
		//	int cluster_index = (int)clusters.size() - 1;
		//	int scan_index = 0;
		//	int j = 0;
		//	if (cluster_index >= 0)
		//	{
		//		scan_index = (int)clusters[cluster_index].size() - 1;
		//		j = clusters[cluster_index][scan_index] + 1;

		//		for (int k = 0; k < clusters[0].size(); k++){
		//			last_cluster.push_back(k);
		//		}
		//		clusters.erase(clusters.begin());
		//	}
		//	for (; j <= i; j++)
		//	{
		//		last_cluster.push_back(j);
		//	}
		//	clusters.push_back(last_cluster);
		//}

	}
#ifdef DEBUG
	int count = 0;
	for (unsigned i = 0; i<clusters.size(); ++i){
		for (unsigned j = 0; j<clusters[i].size(); ++j){

			std::cout << clusters[i][j] << " ";
		}
		std::cout << " count " << count; count++;
		std::cout << std::endl;
	}
#endif
	return clusters;
}

void DWA::motorThread(int cport){
	if (!global_laserscan.empty()){
		arima_laserscan scans = global_laserscan.front();
		global_laserscan.pop();
		unsigned char motor[6];
		int index;

		index = cal_DWA(scans);

		std::fprintf(fptr, "%i vw index:%i\n", count_frame - 1, index);
		//---------------------
		//  障礙物太近時急停
		//---------------------
		if (Lock){
			motor[0] = 0x53;
			motor[1] = 0x10;
			motor[2] = 0x02;
			motor[3] = 0x0;
			motor[4] = 0x0;
			motor[5] = 0x45;
			RS232_SendBuf(cport, motor, 6);
			std::cout << "STAHP" << endl;
			Sleep(100);
			/*	motor[1] = 0x13;
				motor[2] = 0x80;
				if (pre_v > 50){
				motor[3] = 0x80;
				pre_v = 101;
				}
				else if (pre_v < 50){
				motor[3] = 0x7F;
				pre_v = -1;
				}
				motor[4] = 0xFF;
				for (int i = 0; i < 1; ++i){
				RS232_SendBuf(cport, motor, 6);
				Sleep(200);
				}*/
			Lock = false;
			return;
		}

		//------------------------
		// manual control
		//-----------------------
		//---temporarily stop-----
		if (head_direction == 360){
			motor[0] = 0x53;
			motor[1] = 0x10;
			motor[2] = 0x02;
			motor[3] = 0x0;
			motor[4] = 0x0;
			motor[5] = 0x45;
			RS232_SendBuf(cport, motor, 6);
			return;
		}
		//---back up--------------
		else if (head_direction == 300){
			motor[0] = 0x53;
			motor[1] = 0x13;
			motor[2] = 0xA0;
			motor[3] = 0x0;
			motor[4] = 0xFF;
			motor[5] = 0x45;
			RS232_SendBuf(cport, motor, 6);
			return;
		}
		if (head_direction == 90)
			index = 50;
			else if (head_direction == 110)
			index = 40;
			else if (head_direction == 70)
			index = 60;
		//------------------------
		//-------------------------------------
		// Give cmd to motor according to index
		//-------------------------------------
		printf("index:%i\n", index);
		motor[0] = 0x53;
		motor[1] = 0x13;
		motor[2] = 0x7F;
		if (index == 50)
			motor[3] = 0;
		else if (index < 50 && index >= 0)
			// 128,130.54,133.08....252.46
			motor[3] = ((2.54)*index + 128);
		else if (index > 50 && index <= 100){
			motor[3] = (127 - (2.54)*(100 - index));
		}
		else{
			motor[2] = 0x0;
			motor[3] = 0x0;
		}
		motor[4] = 0xFF;
		motor[5] = 0x45;
		std::printf("motor: 0x%X\n", motor[3]);
		RS232_SendBuf(cport, motor, 6);
	}
	
	return;
}

void DWA::laserThread(int cport){
	arima_laserscan laserscan;
	// angle_calib = 校正雷射頭安裝時的角度偏差
	int rx_length, angle_calib = 8, L60 = 0, M60 = 0, R60 = 0, count=0;
	float angle = 0, pre_angle = 0, dist = 0, start = 0;
	vector<float> angles, ranges;
	unsigned char buf[6] = { 0 };
	// 計算1.5公尺內的障礙物
	count_obs = 0;
	// 計算0.4公尺內的障礙物
	count_backup = 0;
	std::fprintf(fptr, "%i cycle\n", count_frame);
	while (1){
		rx_length = RS232_PollComport(cport, buf, 1);
		if ((rx_length == 1) && (buf[0] == 0x55)){
			rx_length = RS232_PollComport(cport, &buf[1], 5);
			if ((rx_length == 5) && (buf[5] == 0xff))
			{
				angle = ((float)(buf[1] << 2)) + ((float)(buf[2] >> 6)) + ((float)(buf[2] & 0x3F) / 64.0F) - angle_calib;
				dist = (((float)((buf[3] & 0x3F) << 8)) + (float)buf[4]) / 1000.0F;
				count++;
				if (angle < 0)
					angle += 360;
				//fprintf(fptr, "%f, %f\n", angle, dist);
				// start為一圈scan的起始角度
				if (angles.empty())
					start = angle;
				// check if data gap > 30
				else if ((angle - pre_angle >= 30) ||
					(pre_angle>angle && pre_angle - angle <= 330)){
					return;
				}
				else if (start < angle + 1.0 && start > angle - 1.0)
					break;
				pre_angle = angle;
				//printf("angle = %f, dist = %f\n\n", angle, dist);
				if (dist < 5 && dist > 0 && angle > 120 && angle < 240){
					angles.push_back(angle);
					ranges.push_back(dist);
					/*if (angle > 90 && angle <= 150 && dist <= 0.8)
						L60++;
					else if (angle > 150 && angle <= 210 && dist <= 0.8)
						M60++;
					else if (angle > 210 && angle <= 270 && dist <= 0.8)
						R60++;*/
					if (dist <= 1.0 && angle > 150.0 && angle < 210){
						count_obs++;
						if (dist <= 0.3)
							count_backup++;
					}
				}

			}
		}
	}
	//--5Hz----
	//Sleep(100);
	//---------
	//std::cout << " count " << count << endl;
	if (angles.size() <= 10){
		return;
	}
	//if (L60 >= 40 && M60 >= 40 && R60 >= 40){
	//	trap = true;
	//	Lclose = true;
	//	Rclose = true;
	//}
	//else if (L60 <= 10){
	//	Lclose = false;
	//	
	//}
	//else if (R60 <= 10){
	//	Rclose = false;
	//	
	//}
	vector<float> angles_out, ranges_out;
	Laser_Filters LF;
	LF.median_filter(angles, ranges, angles_out, ranges_out);
	laserscan.angles = angles_out;
	laserscan.ranges = ranges_out;
	global_laserscan.push(laserscan);

	return;
}

void DWA::robotHead(int cport1){
	int rx_length;
	unsigned char buf[20] = { 0 };
	int turn;
		LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds, Frequency;
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartingTime);
		while (1){
			rx_length = RS232_PollComport(cport1, buf, 2);
			if ((rx_length == 2) && (buf[0] == 0x53) && (buf[1] == 0x20))
			{
				rx_length = RS232_PollComport(cport1, &buf[2], 18);
				QueryPerformanceCounter(&EndingTime);
				ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
				ElapsedMicroseconds.QuadPart *= 1000000;
				ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
				//------------------------------------
				// 取回angular velocity
				// 計算機器人目前頭向
				//------------------------------------
				short v;
				v= ((buf[4] <<8) | (buf[5]));
				short w1;
				w1 = ((buf[6] << 8) | (buf[7]));
				w1 *= -1;
				printf("Velocity %d, Omega %d\n",v, w1);
				//printf("short %hd\n", v);
				/*if (buf[6] == 0xFF)
					turn = 1;
				else
					turn = -1;				
				float w = buf[7];
				if (w <= 127){
					w = max_omega / 127.0F * w;
				}
				else{
					w = max_omega / 127.0F * (255.0F - w);
				}
				w *= turn;
				printf("%f zzz\n", w);*/
				robot_head = robot_head + ((w1*((float)ElapsedMicroseconds.QuadPart / 1000000 - time_lapsed)) / M_PI*180.0F);
				if (robot_head > 360)
					robot_head -= 360;
				else if (robot_head < 0)
					robot_head += 360;
				time_lapsed = (float)ElapsedMicroseconds.QuadPart / 1000000;
			}
			if (flag)
				break;
		}
	return;
}
void DWA::twoThread(){

	int bdrate1 = 19200;      
	int	bdrate2 = 230400;       /* 230400 baud */
	int cport1 = 3;
	int cport2 = 5;
	char mode[] = { '8', 'N', '1', 0 };//8-bit data; no parity bit, 1 stop bit
	if (RS232_OpenComport(cport1, bdrate1, mode)){
		std::printf("Can not open comport1\n");
		system("pause");
		return;
	}

	unsigned char start[2],stop[1];//start command
	start[0] = 0xA5;
	start[1] = 0x20;

	if (RS232_OpenComport(cport2, bdrate2, mode))
	{
		std::printf("Can not open comport2\n");
		system("pause");
		return;
	}
	if (RS232_PollComport(cport2, stop, 1) == 0){
		RS232_SendBuf(cport2, start, 2);
	}

	std::thread t(&DWA::robotHead,this,cport1);
	while (1){		
		std::printf("robot head %f\n ", robot_head);
		std::fprintf(fptr, "robot head %f\n ", robot_head);
		std::cout << "start" << endl;
		LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds, Frequency;
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartingTime);
		laserThread(cport2); 
		//Sleep(100);
		if (flag)
			break;
		/*if (trap){
			while (1){
				if (!trap)
					break;
				escapeTrap(cport1, cport2);
			}
			continue;
		}*/
		motorThread(cport1);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		//time_lapsed = (float)ElapsedMicroseconds.QuadPart / 1000000;
		std::cout << (float)ElapsedMicroseconds.QuadPart / 1000000 << " laser time  " << std::endl;
		//std::fprintf(fptr, "time lapsed %f\n ", time_lapsed);
		if (flag)
			break;
		
	}
	t.join();
	return;
}
