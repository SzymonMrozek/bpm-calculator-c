//
//  BPMCalculator.c
//  TempoWorker
//
//  Created by Szymon Mrozek on 16.10.2017.
//  Copyright © 2017 Szymon Mrozek. All 12rights reserved.
//

#include "BPMCalculator.h"
#include "WaveletAnalysisContainer.h"
#include "FloatArrayOperations.h"
#include "WaveletTransformer.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <math.h> 

float ** grabOverlappingBlocks(float * samples, int count, float samplingFrequency, int * numberOfResultBuffers, int * resultSingleBufferLength);

int calculateBPM(float * samples, int numberOfSamples, float samplingFrequency) {
    
    int * numberOfBuffers = (int *) malloc(sizeof(int));
    int * numberOfSamplesInBlock = (int *) malloc(sizeof(int));
    
    float ** block = grabOverlappingBlocks(samples, numberOfSamples, samplingFrequency, numberOfBuffers, numberOfSamplesInBlock);
    if (*numberOfSamplesInBlock == 0 || *numberOfBuffers == 0) { return 40; }
    
    int * minimumBufferSize = (int *) malloc(sizeof(int));
    
    // 4 Levels Discrete Wavelet Transform for X of 3s 50% overlapping buffers
    float ** summedSignal = generateWavelet(block, *numberOfBuffers, *numberOfSamplesInBlock, minimumBufferSize);

    float ** correlatedBuffers = (float **) malloc(*numberOfBuffers * sizeof(float *));
    
    for(int i = 0; i < *numberOfBuffers; i ++) {
        correlatedBuffers[i] = autocorrelation(summedSignal[i], *minimumBufferSize);
    }
    
    float * bpm = (float *) malloc((numberOfSamples - 1) * sizeof(float));
    
    for(int i = 0; i < numberOfSamples - 1; i ++) {
        bpm[i] = 60.0 / ((float)i / (samplingFrequency / 16.0));
    }
    
    float * bpmsInRange = (float *) malloc((numberOfSamples - 1) * sizeof(float));
    int * wantedIndexes = (int *) malloc((numberOfSamples - 1) * sizeof(int));
    
    int numberOfElements = 0;
    for(int i = 0; i < numberOfSamples - 1; i ++) {
        float bpmValue = bpm[i];
        if (bpmValue >= 40.0 && bpmValue <= 220.0) {
            wantedIndexes[numberOfElements] = i;
            bpmsInRange[numberOfElements] = bpmValue;
            numberOfElements ++;
        }
    }
    
    bpmsInRange = (float *) realloc(bpmsInRange, numberOfElements * sizeof(float));
    wantedIndexes = (int *) realloc(wantedIndexes, numberOfElements * sizeof(int));
    
    float * estimatedBPMs = (float *) malloc(*numberOfBuffers * sizeof(float));
    float * currentCorrel = (float *) malloc(numberOfElements * sizeof(float));
    
    for(int i = 0; i < *numberOfBuffers; i ++) {
        for( int j = 0; j < numberOfElements; j ++) {
            int readIndex = wantedIndexes[j];
            currentCorrel[j] = correlatedBuffers[i][readIndex];
        }
        int * index = maxIndex(currentCorrel, numberOfElements);
        if( index == NULL) {
            continue;
        }
        estimatedBPMs[i] = bpmsInRange[*index];
        free(index);
    }

    int resulut =  (int)roundf(median(estimatedBPMs, *numberOfBuffers));
    
    free(minimumBufferSize);
    for(int i = 0; i < *numberOfBuffers; i ++) {
        free(correlatedBuffers[i]);
        free(summedSignal[i]);
    }
    free(summedSignal);
    free(correlatedBuffers);
    free(bpm);
    free(bpmsInRange);
    free(wantedIndexes);
    free(estimatedBPMs);
    free(currentCorrel);
    free(numberOfBuffers);
    free(numberOfSamplesInBlock);

    return resulut;
}

float ** grabOverlappingBlocks(float * samples, int count, float samplingFrequency, int * numberOfResultBuffers, int * resultSingleBufferLength) {
    
    // This method split array into 50% overlapping sequences of 3s length
    // Works as matlab `buffer`
    // Assumes that processing block length is even number
    
    float overlapRatio = 0.5; // 50% overlaping
    int threeSecondsBlockLength = 3 * (int)samplingFrequency;
    int overlappedBlockLength = (int)(threeSecondsBlockLength * overlapRatio);
    
    int prototypeNumberOfBuffers = (int) ((double)count / (double)overlappedBlockLength) + 2;
    int realNumberOfBuffers = 0;
    if (count < threeSecondsBlockLength) { return NULL; }

    float ** result = (float **) malloc(prototypeNumberOfBuffers * sizeof(float *));
    for(int i = 0; i < prototypeNumberOfBuffers; i ++) {
        result[i] = (float *) malloc(threeSecondsBlockLength * sizeof(float));
    }
    
    
    // First block contains zeros of count = overlappedBlockLength
    fillWithZeros(result[0], overlappedBlockLength);
    copy(result[0] + overlappedBlockLength, samples, threeSecondsBlockLength - overlappedBlockLength);
    
    int samplesIterator = threeSecondsBlockLength;

    // Blocks in the middle:
    for(; samplesIterator < count - overlappedBlockLength; samplesIterator += (threeSecondsBlockLength - overlappedBlockLength)) {
        realNumberOfBuffers ++;
        copy(result[realNumberOfBuffers], samples + (samplesIterator - overlappedBlockLength) , overlappedBlockLength);
        copy(result[realNumberOfBuffers] + overlappedBlockLength, samples + samplesIterator, (threeSecondsBlockLength - overlappedBlockLength));
    }
    
    // Last block
    realNumberOfBuffers ++;
    copy(result[realNumberOfBuffers], samples + (samplesIterator - overlappedBlockLength), overlappedBlockLength);
    int lastSegmentCount = count - samplesIterator;
    copy(result[realNumberOfBuffers] + overlappedBlockLength, samples + samplesIterator, lastSegmentCount);
    int zerosToFill = threeSecondsBlockLength - (lastSegmentCount + overlappedBlockLength);
    fillWithZeros(result[realNumberOfBuffers] + (overlappedBlockLength + lastSegmentCount), zerosToFill);


    for(int i = realNumberOfBuffers + 1; i < prototypeNumberOfBuffers; i ++) {
        free(result[i]);
    }
    
    *resultSingleBufferLength = threeSecondsBlockLength;
    *numberOfResultBuffers = realNumberOfBuffers + 1;
    
    free(samples);
    
    return result;
}
