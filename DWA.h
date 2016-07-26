#include <vector>
#include <queue>
#include <thread>

class DWA{
	typedef std::vector<std::vector<int>> cluster_vector;
	typedef std::vector<std::vector<std::pair<float, float>>> VMap;
public:
	DWA(){
		errno_t err;
		err = fopen_s(&fptr, "Log.txt", "w");
		if (err != 0){
			printf("open file failure...\n");
		}
		head_direction = 90;
		targetMask();
	}
	~DWA(){
		fclose(fptr);
	}
	void twoThread();
	volatile bool flag = false;
	volatile int head_direction ;
private:
	struct arima_laserscan
	{
		std::vector<float> ranges;
		std::vector<float> angles;
	};
	struct slope_range
	{
		float r_right;
		float r_left;
		float r;
		float y;
	};
	void targetMask(){
		for (int i = 0; i < 101; i++){
			if (i < 50)
				mask.push_back((float)i / 50.0F);
			else if (i == 50)
				mask.push_back(1);
			else
				mask.push_back((float)(100 - i) / 50.0F);
		}
	}
	void robotHead(int);
	void motorThread(int);
	void laserThread(int);
	int cal_DWA(arima_laserscan);
	cluster_vector cluster(arima_laserscan&);
	std::queue<arima_laserscan> global_laserscan;
	std::vector<float> mask;
	int count_frame = 0;
	int count_obs, count_backup;
	int pre_v = 50;
	volatile float robot_head = 90.0F;
	float time_lapsed = 0;
	bool Lock = false, trap = false, Lclose = false, Rclose = false;
	FILE* fptr;
};