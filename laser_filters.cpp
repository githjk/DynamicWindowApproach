#include "laser_filters.h"



Laser_Filters::Laser_Filters() 
{
   ;
}
   
Laser_Filters::~Laser_Filters()
{
   ;
}

bool Laser_Filters::median_filter(std::vector<float>& laser_scans_angle_in, std::vector<float>& laser_scans_dist_in, std::vector<float>& laser_scans_angle_out, std::vector<float>& laser_scans_dist_out)
{
   laser_scans_dist_out = laser_scans_dist_in;
   laser_scans_angle_out = laser_scans_angle_in;
   
   //
   int median_filter_radius = 2;
   for(int i = 0; i < laser_scans_dist_in.size(); i++)
   {
      if((i < median_filter_radius) || (i >= laser_scans_dist_in.size()-median_filter_radius))
      {
         laser_scans_dist_out[i] = laser_scans_dist_in[i];
      }
      else
      {
         float median = find_median(i, laser_scans_dist_in, median_filter_radius);
         if(median == -1)
         {
            laser_scans_dist_out[i] = 0;
            return 0;   
         }
         else
         {
            laser_scans_dist_out[i] = median;
         }
      }
   }
   
   return 1;
}

float Laser_Filters::find_median(int index, std::vector<float>& laser_scans_dist_in, int radius)
{
   int i = 0;
   int j = 0;
   int l = 0;
   int m = radius*2;
  
   // check the boundary of laser_scans_dist_in
   if((index - radius < 0) || (index + radius >= laser_scans_dist_in.size()))
   {
      return -1;
   }
   
   // copy data from laser_scans_dist_in
   float* laser_scans_dist_local = new float[radius*2+1]; 
   for(i = -radius; i <= radius; i++)
   {
      laser_scans_dist_local[i+radius] = laser_scans_dist_in[index+i];
   }
   
   // find the median
   float x;
   while(l < m)
   {
      x = laser_scans_dist_local[radius];
      i = l;
      j = m;
      do
      {
         while(laser_scans_dist_local[i] < x) i++;
         while(x < laser_scans_dist_local[j]) j--;
         if(i<=j)
         {
            // swap
            float temp = laser_scans_dist_local[i];
            laser_scans_dist_local[i] = laser_scans_dist_local[j];
            laser_scans_dist_local[j] = temp;
            i++; 
            j--;
         }
      } while (i <= j);
      if(j < radius) l = i;
      if(radius < i) m = j;
   }  
   
   // debug - make sure the returned variable is median
   /*
   for(i = -radius; i <= radius; i++)
   {
      if(i > 0)
      {
         if(laser_scans_dist_local[radius] > laser_scans_dist_local[radius+i])
         {
            ROS_ERROR("*** [laser_filters] find incorrect median ***");
         }
      }
      else if(i < 0)
      {
         if(laser_scans_dist_local[radius] < laser_scans_dist_local[radius+i])
         {
            ROS_ERROR("*** [laser_filters] find incorrect median ***");
         }      
      }
   }
   */ 
  
   // delete and return
   float ret_val = laser_scans_dist_local[radius];
   delete [] laser_scans_dist_local;
   return ret_val;
}


