#ifndef __TIME_RECORD_H__
#define __TIME_RECORD_H__
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif


//timer
void timerStart();
void timerEnd();
double elapsedTime();

//timer1
void timer1Start();
void timer1End();
double elapsedTime1();

void printTime_toStandardOutput();
void printTime_toFile();

typedef struct timeStruct {
	double iniTime;
	double preprocessingTime;
	double copyTimeHostToDevice;
	double kernelTime;
	double copyTimeDeviceToHost;
	double postprocessingTime;
	double totalTime;
} STRUCT_TIME;

extern struct timeval timeStart;
extern struct timeval timeEnd;
extern struct timeval time1Start;
extern struct timeval time1End;

extern STRUCT_TIME strTime;

#ifdef __cplusplus
}
#endif

#endif //__TIME_RECORD_H__
