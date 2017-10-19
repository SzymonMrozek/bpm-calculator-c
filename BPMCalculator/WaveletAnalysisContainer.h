//
//  WaveletAnalysisContainer.h
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#ifndef WaveletAnalysisContainer_h
#define WaveletAnalysisContainer_h

#include <stdio.h>

struct WaveletAnalysisContainer {
    
    float *** data;
    int topLevelLength;
    int middleLevelLength;
    int bottomLevelLength;
    
};

struct WaveletAnalysisContainer initializeContainer(int topLevelLength, int middleLevelLength, int bottomLevelLength);
void releaseContainer (struct WaveletAnalysisContainer *container,int topOffset, int topLevels, int middleLevels);

#endif /* WaveletAnalysisContainer_h */
