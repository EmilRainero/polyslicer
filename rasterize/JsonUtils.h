#include "json/json.h"

#ifndef JSON_UTILS_H
#define JSON_UTILS_H

class JsonUtils {
public:
    static Json::Value parse(const std::string& json, const bool replaceQuoteToDoubleQuote = false);

    static std::string toString(const Json::Value& root, const bool minimal = false);

    static Json::Value readFromFile(const std::string& filename);

    static void writeToFile(const Json::Value& json, const std::string& filename);
};

#endif
