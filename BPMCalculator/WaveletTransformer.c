//
//  WaveletTransformer.c
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#include "WaveletTransformer.h"
#include "WaveletAnalysisContainer.h"
#include "FloatArrayOperations.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <float.h>

int wrap(int value, int lowBound, int highBound);
int nonNegativeModulo(int value, int modulo);

float ** generateWavelet(float ** input, int numberOfBuffers, int bufferCount, int * minimumBufferSize) {
    
    if (numberOfBuffers == 0) { return NULL; }
    
    int levels = 4;
    
    struct WaveletAnalysisContainer approxCoefficients = initializeContainer(levels + 1, numberOfBuffers, 0);
    struct WaveletAnalysisContainer detailCoefficients = initializeContainer(levels + 1, numberOfBuffers, 0);
    struct WaveletAnalysisContainer analysisResult     = initializeContainer(levels + 1, numberOfBuffers, 0);

    free(approxCoefficients.data[0]);
    
    approxCoefficients.data[0] = input;

    int undersamplingFactors [5] = { 16, 8, 4, 2, 1 };
    
    int previousLevelCount = bufferCount;
    int * approxCount       = (int *) malloc(sizeof(int));
    int * detailCount       = (int *) malloc(sizeof(int));
    int * undersampledCount = (int *) malloc(sizeof(int));
    
    *minimumBufferSize = previousLevelCount;
    
    for(int bufferIndex = 0; bufferIndex < numberOfBuffers; bufferIndex ++) {
        
        previousLevelCount = bufferCount;
        
        for(int level = 1; level < levels + 1; level ++) {
            
            float * wavelet = forwardTransform(approxCoefficients.data[level-1][bufferIndex], previousLevelCount, approxCount, detailCount);
            filtered(wavelet + *approxCount, *detailCount); // Filter coefficients before assign

            detailCoefficients.data[level][bufferIndex] = copyAndNew(wavelet + *approxCount, *detailCount);
            approxCoefficients.data[level][bufferIndex] = detailCoefficients.data[level][bufferIndex];
        
            float * undersampled = undersample(detailCoefficients.data[level][bufferIndex], *detailCount, undersamplingFactors[level], undersampledCount);
            normalized(undersampled, *undersampledCount);
            analysisResult.data[level][bufferIndex] = undersampled;
            
            *minimumBufferSize = (int) fmin((double) *minimumBufferSize, (double) *undersampledCount);
            previousLevelCount /= 2;
            
            free(wavelet);
        }

    }

    float ** summedSignal = (float **) malloc(numberOfBuffers * sizeof(float*));
    
    for(int i = 0; i < numberOfBuffers; i ++) {
        summedSignal[i] = (float *) malloc(*minimumBufferSize * sizeof(float));
        for(int j = 0; j < *minimumBufferSize; j ++) {
            summedSignal[i][j] = 0;
        }
    }

//  Sum 4 levels into 1
    
    for(int level = 1; level < levels + 1; level ++) {
        for(int buffer = 0; buffer < numberOfBuffers; buffer ++) {
            for(int bufferSample = 0; bufferSample < *minimumBufferSize; bufferSample ++) {
                summedSignal[buffer][bufferSample] += analysisResult.data[level][buffer][bufferSample];
            }
        }
    }

    
//  Last level coefficients

    for(int buffer = 0; buffer < numberOfBuffers; buffer ++) {
        filtered(approxCoefficients.data[levels][buffer], *minimumBufferSize);
        normalized(approxCoefficients.data[levels][buffer], *minimumBufferSize);
        for(int bufferSample = 0; bufferSample < *minimumBufferSize; bufferSample ++) {
            summedSignal[buffer][bufferSample] += approxCoefficients.data[levels][buffer][bufferSample];
        }
    }
    
    free(approxCount);
    free(detailCount);
    free(undersampledCount);
    releaseContainer(&approxCoefficients, 0, levels + 1, numberOfBuffers);
    releaseContainer(&detailCoefficients, 0, levels + 1, 0);
    releaseContainer(&analysisResult, 1, levels + 1, numberOfBuffers);

    return summedSignal;
}
    
    
float * forwardTransform(float * input, int count, int * approxCount, int * detailCount) {

    float * tempArray = (float *) malloc(count * sizeof(float));
    float * inputCopy = copyAndNew(input, count);

    for(int i = 0; i < count; i ++) {
        tempArray[i] = 0.0;
    }

    float coefficients [4] = {
        0.4829629131445341,
        0.8365163037378079,
        0.2241438680420133,
        -0.1294095225512603
    };

    int length = count;

    int i, j, j0, j1, j2, j3;

    

    while ( 4 <= length ) {

        i = 0;

        int halfLength = length >> 1;

        for(j = 0; j < length - 1; j += 2) {

            j0 = wrap (j    , 0, length - 1);
            j1 = wrap (j + 1, 0, length - 1);
            j2 = wrap (j + 2, 0, length - 1);
            j3 = wrap (j + 3, 0, length - 1);

            tempArray[i] = coefficients[0] * inputCopy[j0] + coefficients[1] * inputCopy[j1] + coefficients[2] * inputCopy[j2] - coefficients[3] * inputCopy[j3];
            tempArray[i + halfLength] = coefficients[3] * inputCopy[j0] - coefficients[2] * inputCopy[j1] + coefficients[1] * inputCopy[j2] - coefficients[0] *inputCopy[j3];

            i += 1;
        }

        for(i = 0; i < length; i ++) {
            inputCopy[i] = tempArray[i];
        }

        length >>= 1;
    }

    *approxCount = count / 2;
    *detailCount = count / 2;

    free(tempArray);
    return inputCopy;
}

int wrap(int value, int lowBound, int highBound) {

    int wide = highBound + 1 - lowBound;

    if ( wide == 1 ) {
        return lowBound;
    }
    return lowBound + nonNegativeModulo(value - lowBound, wide);

}

int nonNegativeModulo(int value, int modulo) {

    assert(modulo != 0);
    int result = value % modulo;
    if(result < 0) {
        result += abs(modulo);
    }
    return result;
}

