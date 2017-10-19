//
//  WaveletAnalysisContainer.c
//  TempoWorker
//
//  Created by Szymon Mrozek on 15.10.2017.
//  Copyright © 2017 Szymon Mrozek. All rights reserved.
//

#include "WaveletAnalysisContainer.h"
#include <stdlib.h>
#include <stdio.h>

struct WaveletAnalysisContainer initializeContainer(int topLevelLength, int middleLevelLength, int bottomLevelLength) {
    
    struct WaveletAnalysisContainer container;
    
    container.topLevelLength = topLevelLength;
    container.middleLevelLength = middleLevelLength;
    container.bottomLevelLength = bottomLevelLength;
    
    container.data = (float ***) malloc(topLevelLength * sizeof(float **));
    
    for(int i = 0; i < topLevelLength; i ++) {
        container.data[i] = (float **) malloc(middleLevelLength * sizeof(float *));
    }
    
    return container;
}

void releaseContainer (struct WaveletAnalysisContainer *container, int topOffset, int topLevels, int middleLevels) {

    for(int top = 0; top < topLevels; top ++) {
        for(int middle = 0; middle < middleLevels; middle ++) {
            if( top >= topOffset) {
                free(container->data[top][middle]);
            }
        }
        free(container->data[top]);
    }
    free(container->data);
}

