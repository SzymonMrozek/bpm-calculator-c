//
//  WaveletTransformer.h
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#ifndef WaveletTransformer_h
#define WaveletTransformer_h

struct WaveletAnalysisContainer;

float ** generateWavelet(float ** input, int numberOfBuffers, int bufferCount, int * minimumBufferSize);
float * forwardTransform(float * x, int n, int * approxCount, int * detailCount);
float * inverseTransform(float * y, int n);

#endif /* WaveletTransformer_h */
