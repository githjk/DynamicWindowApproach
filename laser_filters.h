// [20150825_Curtis] remove all member variables so that this class only provides functions for filtering laser scans

#include <vector>

class Laser_Filters
{
public:
   Laser_Filters(); 
   ~Laser_Filters();
   bool median_filter(std::vector<float>& laser_scans_angle_in, std::vector<float>& laser_scans_dist_in, std::vector<float>& laser_scans_angle_out, std::vector<float>& laser_scans_dist_out);
   
private:
   float find_median(int index, std::vector<float>& laser_scans_dist_in, int radius = 2); // find the median between laser_scans_dist_in[index-radius] and laser_scans_dist_in[index+radius]
};

