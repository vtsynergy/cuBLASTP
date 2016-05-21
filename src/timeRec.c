#include "timeRec.h"
#include <stdio.h>
#include <stdlib.h>

struct timeval timeStart;
struct timeval timeEnd;
struct timeval time1Start;
struct timeval time1End;

STRUCT_TIME strTime;

void timerStart()
{
	gettimeofday(&timeStart, NULL);
}

void timerEnd()
{
	gettimeofday(&timeEnd, NULL);
}

void timer1Start()
{
	gettimeofday(&time1Start, NULL);
}

void timer1End()
{
	gettimeofday(&time1End, NULL);
}

//return value is ms
double elapsedTime()
{
	double deltaTime;
	deltaTime = (timeEnd.tv_sec - timeStart.tv_sec) * 1000.0 + 
				(timeEnd.tv_usec - timeStart.tv_usec) / 1000.0;
	
	return deltaTime;
}

double elapsedTime1()
{
	double deltaTime;
	deltaTime = (time1End.tv_sec - time1Start.tv_sec) * 1000.0 + 
				(time1End.tv_usec - time1Start.tv_usec) / 1000.0;
	
	return deltaTime;
}

void printTime_toStandardOutput()
{
	//calculate the total time
	strTime.totalTime = strTime.preprocessingTime +
						strTime.copyTimeHostToDevice +
						strTime.kernelTime +
						strTime.copyTimeDeviceToHost +
						strTime.postprocessingTime;
	printf("Version 21: %.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
			strTime.preprocessingTime,
			strTime.copyTimeHostToDevice,
			strTime.kernelTime,
			strTime.copyTimeDeviceToHost,
			strTime.postprocessingTime,
			strTime.totalTime);

	return;
}

void printTime_toFile()
{
	FILE *pTimeFile;
	pTimeFile = fopen("../runtimeGapped.txt", "at");
	if (pTimeFile == NULL)
	{
		printf("File runtime.txt open error!\n");
		return;
	}
	
	//calculate the total time
	strTime.totalTime = strTime.preprocessingTime +
						strTime.copyTimeHostToDevice +
						strTime.kernelTime +
						strTime.copyTimeDeviceToHost +
						strTime.postprocessingTime;
	fprintf(pTimeFile, "Version 21: %.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
			strTime.preprocessingTime,
			strTime.copyTimeHostToDevice,
			strTime.kernelTime,
			strTime.copyTimeDeviceToHost,
			strTime.postprocessingTime,
			strTime.totalTime);

	fclose(pTimeFile);

	return;
}


