
#include "DWA.h"
#include "rs232.h"
#include "laser_filters.h"
#include <conio.h>
#include <thread>
#include <Mutex>
#pragma warning (disable : 4996)


void main(int argc, const char** argv){
	DWA *_dwa = new DWA();
	std::thread t(&DWA::twoThread,_dwa);
	//std::thread t = _dwa->exeThread();
	while (1){
		if (kbhit() != 0){
			int c = getch();
			switch (c){
			case 27://ESC
				_dwa->flag = true; break;
			case 100://d
				_dwa->head_direction = 0; break;
			case 101://e
				_dwa->head_direction = 70; break;
			case 119://w
				_dwa->head_direction = 90; break;
			case 113://q
				_dwa->head_direction = 110; break;
			case 97://a
				_dwa->head_direction = 180; break;
			case 115://s for stop
				_dwa->head_direction = 360; break;
			case 120://x for backup
				_dwa->head_direction = 300; break;
			}
		}
		if (_dwa->flag)
			break;
	}
	t.join();
	system("pause");
	return;
}