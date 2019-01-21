//----------------------------------------------------------------------
// File: main.cpp
// Date: 2017-
// Description: unittest for kernel/utils
// Author: PARC
// Copyright(c): All rights reserved
//----------------------------------------------------------------------

#include <ctime>
#include <random>
#include "gtest/gtest.h"
#include "LogCF.h"

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    LogCF *logcf = LogCF::Instance();
    logcf->mute();

    return RUN_ALL_TESTS();
}