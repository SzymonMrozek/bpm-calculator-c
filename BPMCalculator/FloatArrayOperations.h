//
//  FloatArrayOperations.h
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#ifndef FloatArrayOperations_h
#define FloatArrayOperations_h

float meanValue(float* array, int count);

float* undersample(float * array, int count, int decimator, int * undersampledCount);

void magnitude(float* array, int count);

void normalized(float * array, int count);

float* autocorrelation(float * array, int count);

void filtered(float * array, int count);

float* absoluteMax(float * array, int count);

float* floatMax(float * array, int count);

int * maxIndex(float * array, int count);

float median(float * array, int count);

void copy(float * destination, float * source, int count);
void fillWithZeros(float * destination, int count);
float * copyAndNew(float * vector, int count);
float * newAndInitializedToZero(int count);

void printArray(float * array, int count);

#endif /* FloatArrayOperations_h */
