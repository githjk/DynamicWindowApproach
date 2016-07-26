#pragma once
#include <vector>
#include <opencv2/imgproc/imgproc.hpp>

namespace LnExt
{



///////////////////////
#define LASER_ANGLE 120 // 120 or 180

#ifndef M_PI
#define M_PI 3.14159265358979323846F
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD(deg) ((deg) / 180.0F * M_PI)
#endif
   

   //////////////////////
   struct ArimaLaserscan
   {
      std::vector<float> ranges; // meter
      std::vector<float> angles; // deg
   };
   typedef std::vector<std::vector<int>> ClusterVector;


   //////////////////////////////
   class CLineExtractionFromLaser
   {
   private:
      float m_minLineLengthM;
      float m_inlierRatio;
      int m_minPtInCluster;
      int m_maxTrialCnt;
      float m_maxDistToLine;
   private:
      void ComputeLineEq(const cv::Point2f sp, const cv::Point2f ep, float& thetaRad, float& rho);
      cv::Point2f CvtLaserScanToXY(const float angleRad, const float rangeM, const float angularOffset = 0, const float scale = 1, const cv::Point2f offset = cv::Point2f(0, 0));
      inline float ComputeDistFromPtToLine(const float thetaRad, const float rho, cv::Point2f pt);
      inline float InnerProduct(const cv::Point2f ptL, const cv::Point2f ptR, const cv::Point2f ptO);
      inline float ComputeDistBetweenTwoPt(const cv::Point2f pt1, const cv::Point2f pt2);
   public:
      CLineExtractionFromLaser();
      ~CLineExtractionFromLaser();
      bool LineExtract(const ArimaLaserscan laserScans, const ClusterVector clusters, int& sIdx, int& eIdx);
   };
}


