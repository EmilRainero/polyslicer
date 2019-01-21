
#include "Path.h"
#include "gtest/gtest.h"

namespace {
    TEST(PathTest, Join) {
        ASSERT_EQ("hello/world", Path::join("hello", "world"));
        ASSERT_EQ("hello/world", Path::join("hello", "/world"));
        ASSERT_EQ("hello/world", Path::join("hello/", "world"));
        ASSERT_EQ("hello//world", Path::join("hello/", "/world"));
        ASSERT_EQ("hello", Path::join("hello"));
        ASSERT_EQ("hello/world/parc", Path::join("hello", "world", "parc"));
    }
}