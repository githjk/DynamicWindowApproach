#include "LineExtractionFromLaser.h"

using namespace LnExt;

CLineExtractionFromLaser::CLineExtractionFromLaser()
{
   m_minLineLengthM = 1.5; // (warning) TBD
   m_inlierRatio = 0.25F; // inlier / total (warning) TBD
   m_minPtInCluster = (int)(LASER_ANGLE * m_inlierRatio / 2); // 1/2 #inliers
   m_maxTrialCnt = (int)(log(1 - 0.999) / log(1 - m_inlierRatio*m_inlierRatio) + 0.5);
   m_maxDistToLine = 0.1F; // (unit: meter) if the distance between a point and the computed line is smaller than this threshold, consider this point as an inlier
}


CLineExtractionFromLaser::~CLineExtractionFromLaser()
{
}

bool CLineExtractionFromLaser::LineExtract(const ArimaLaserscan laserScans, const ClusterVector clusters, int& sIdx, int& eIdx)
{
   // find largest cluster
   int largestClusterIdx = -1;
   int largestClusterNum = 0;
   for (unsigned int i = 0; i < clusters.size(); i++)
   {
      if (clusters[i].size() > largestClusterNum)
      {
         largestClusterIdx = i;
         largestClusterNum = (int)clusters[i].size();
      }
   }
   // RANSAC
   bool hasValidLineSegment = false;
   int bestFrom = -1;
   int bestTo = -1;
   int bestInlierNum = -1;
   //
   if (largestClusterNum > m_minPtInCluster) // apply RANSAC only in the largest cluster which has more laser scans than "minPtInCluster" 
   {
      for (int t = 0; t < m_maxTrialCnt; t++)
      {
         // select points from the largest cluster
         int from = rand() % clusters[largestClusterIdx].size();
         from = clusters[largestClusterIdx][from];
         int to = rand() % clusters[largestClusterIdx].size();
         to = clusters[largestClusterIdx][to];
         while (from == to)
         {
            to = rand() % clusters[largestClusterIdx].size();
            to = clusters[largestClusterIdx][to];
         }

         // compute an equation for line
         cv::Point2f pt1 = CvtLaserScanToXY(DEG_TO_RAD(laserScans.angles[from]), laserScans.ranges[from]);
         cv::Point2f pt2 = CvtLaserScanToXY(DEG_TO_RAD(laserScans.angles[to]), laserScans.ranges[to]);

         float theta = 0;
         float rho = 0;
         ComputeLineEq(pt1, pt2, theta, rho);

         // compute the distance between laser scans in the largest cluster and the computed line, and find the line segment (sp <---> ep)
         int inlierCnt = 0;
         cv::Point2f sp = pt1;
         cv::Point2f ep = pt2;
         // [20151203_Curtis] only laser scans in the largest cluster can vote
         for (unsigned int i = 0; i < clusters[largestClusterIdx].size(); i++)
         {
            int idx = clusters[largestClusterIdx][i];
            cv::Point2f pt = CvtLaserScanToXY(DEG_TO_RAD(laserScans.angles[idx]), laserScans.ranges[idx]);
            float dist = ComputeDistFromPtToLine(theta, rho, pt);
            if (dist < m_maxDistToLine)
            {
               inlierCnt++;
               bool updateSP = (InnerProduct(pt, ep, sp) < 0);
               bool updateEP = (InnerProduct(pt, sp, ep) < 0);
               if (updateSP)
               {
                  sp = pt;
                  from = idx;
               }
               if (updateEP)
               {
                  ep = pt;
                  to = idx;
               }
               if (updateSP && updateEP)
               {
                  printf("===========================\n");
                  printf("Warning: check what happens\n");
                  printf("===========================\n");
               }
            }
         }

         // check if the found line segment is long enough and has enough inliers
         // besides, try to find as many as inliers instead the longest distance
         float dist = ComputeDistBetweenTwoPt(sp, ep);
         if (inlierCnt > m_minPtInCluster && dist > m_minLineLengthM && inlierCnt > bestInlierNum)
         {
            hasValidLineSegment = true;
            bestFrom = from;
            bestTo = to;
            bestInlierNum = inlierCnt;
         }
      }
   }
 
   //
   sIdx = bestFrom;
   eIdx = bestTo;
   return hasValidLineSegment;
}

// ------------------------------------------------------
// function to compute equation of a line from two points
// ------------------------------------------------------
void CLineExtractionFromLaser::ComputeLineEq(const cv::Point2f sp, const cv::Point2f ep, float& thetaRad, float& rho)
{
   // use "x*cos(thetaRad) + y*sin(thetaRad) = c" as line equation   
   thetaRad = atan2f(sp.x - ep.x, sp.y - ep.y);
   rho = sp.x * cos(thetaRad) + sp.y * sin(thetaRad);
   return;
}

// ----------------------------------------------------------
// function to compute the (x, y) position of a given laser scan
// ----------------------------------------------------------
cv::Point2f CLineExtractionFromLaser::CvtLaserScanToXY(const float angleRad, const float rangeM, const float angularOffset, const float scale, const cv::Point2f offset)
{
   cv::Point2f xyPt;
   // 1). add angularOffset to angle
   float angle = angleRad + angularOffset;
   // 2). multiple scale when compute (x, y)
   xyPt.x = rangeM * cos(angle) * scale;
   xyPt.y = rangeM * sin(angle) * scale;
   // 3). add offset to computed (x, y)
   xyPt.x += offset.x;
   xyPt.y += offset.y;

   return xyPt;
}

inline float CLineExtractionFromLaser::ComputeDistFromPtToLine(const float thetaRad, const float rho, cv::Point2f pt)
{
   // dist(ax + by + c = 0, (x0, y0)) = |ax0 + by0 + c| / sqrt(a * a + b * b)
   return abs(pt.x * cos(thetaRad) + pt.y * sin(thetaRad) - rho);
}

inline float CLineExtractionFromLaser::InnerProduct(const cv::Point2f ptL, const cv::Point2f ptR, const cv::Point2f ptO)
{
   //  ptL     PtR
   //    \    /    : compute the inner product of ptO->ptL and ptO->ptR 
   //      ptO 
   return ((ptL.x - ptO.x) * (ptR.x - ptO.x) + (ptL.y - ptO.y) * (ptR.y - ptO.y));
}

inline float CLineExtractionFromLaser::ComputeDistBetweenTwoPt(const cv::Point2f pt1, const cv::Point2f pt2)
{
   return sqrt(pow(pt1.x - pt2.x, 2) + pow(pt1.y - pt2.y, 2));
}
