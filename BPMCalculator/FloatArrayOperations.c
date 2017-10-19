//
//  FloatArrayOperations.c
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#include "FloatArrayOperations.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"

float meanValue(float* array, int count) {
    
    if (count == 0) { return 0.0; }
    float mean = 0.0;
    for(int i = 0; i < count; i ++) {
        mean += array[i] / (float) count;
    }
    return mean;
}

float* undersample(float * array, int count, int decimator, int * undersampledCount) {
    
    float * copy = (float *) malloc(count * sizeof(float));
    
    int j = 0;
    for(int i = 0; i < count; i ++) {
        if ( i % decimator == 0 ) {
            copy[j] = array[i];
            j++;
        }
    }
    
    *undersampledCount = j;
    copy = realloc(copy, j * sizeof(float));
    
    return copy;
}

void magnitude(float* array, int count) {
    for(int i = 0; i < count; i ++) {
        array[i] = fabsf(array[i]);
    }
}

void normalized(float * array, int count) {
 
    float mean = meanValue(array, count);
    
    magnitude(array, count);
    
    for(int i = 0; i < count; i ++) {
        array[i] = array[i] - mean;
    }
}

float* autocorrelation(float * array, int count) {

    float * real = copyAndNew(array, count);
    float * imag = newAndInitializedToZero(count);


    Fft_transform(real, imag, count);
    
    for(int i = 0; i < count; i ++ ) {
        real[i] = real[i] * real[i] + imag[i] * imag[i];
        imag[i] = 0.0;
    }
    
    Fft_inverseTransform(real, imag, count);
    
    free(imag);
    return real;
}

void filtered(float * array, int count) {
    array[0] = 0.0;
    for(int i = 1; i < count; i++) {
        array[i] = 0.01 * array[i] + 0.99 * array[i-1];
    }
}
    
float* absoluteMax(float * array, int count) {
 
    if (count == 0) { return NULL; }
    
    float * copy = copyAndNew(array, count);
    magnitude(copy, count);
    return floatMax(copy, count);
}

float* floatMax(float * array, int count) {
    
    if (count == 0) {
        return NULL;
    }
    
    float maxValue = array[0];
    
    for(int i = 1; i < count; i ++) {
        maxValue = fmax(array[i], maxValue);
    }
    
    float * max = malloc(sizeof(float));
    *max = maxValue;
    return max;
}

int * maxIndex(float * array, int count) {
    
    float * max = floatMax(array, count);
    if( max == NULL) { return NULL; }
    
    for(int i = 0; i < count; i ++) {
        if(array [i] == *max) {
            int * index = malloc(sizeof(int));
            *index = i;
            free(max);
            return index;
        }
    }
    free(max);
    return NULL;
}

int compare(const void * lhs, const void * rhs) {
    return *(float *)lhs > *(float *)rhs;
}


float median(float * array, int count) {
    
    qsort(array, count, sizeof(float), compare);
    
    if( count % 2 == 0) {
        int midIndex = count / 2;
        return ( array[midIndex] + array[midIndex - 1]) / 2.0;
    }
    return array[ count / 2];
    
}

void copy(float * destination, float * source, int count) {
    
    memcpy(destination, source, count * sizeof(float));
}

void fillWithZeros(float * destination, int count) {
    
    for( int i = 0; i < count; i ++) {
        destination[i] = 0.0;
    }
}


float * copyAndNew(float * vector, int count) {
    
    float * result = (float *) malloc(count * sizeof(float));
    memcpy(result, vector, count * sizeof(float));
    return result;
}

float * newAndInitializedToZero(int count) {
    
    float * result = (float *) malloc(count * sizeof(float));
    for( int i = 0; i < count; i ++) {
        result[i] = 0.0;
    }
    return result;
}
