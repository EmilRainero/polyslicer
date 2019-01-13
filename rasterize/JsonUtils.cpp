#include <iostream>
#include <sstream>
#include <fstream>
//#include <boost/filesystem.hpp>
#include <stdexcept>

#include "JsonUtils.h"


Json::Value
JsonUtils::parse(const std::string& json, const bool replaceQuoteToDoubleQuote) {
    std::string jsonTemp(json);
//    if (replaceQuoteToDoubleQuote) {
//        std::replace(jsonTemp.begin(), jsonTemp.end(), '\'', '"');
//
//    }
    Json::CharReaderBuilder rbuilder;
    std::stringstream jsonStream(jsonTemp);
    Json::Value root;
    std::string errors;

    bool success = Json::parseFromStream(rbuilder, jsonStream, &root, &errors);

    if (!success) {
        throw std::runtime_error("JsonCpp:" + errors);
    }

    return root;
}

std::string
JsonUtils::toString(const Json::Value& root, const bool minimal) {
    Json::StreamWriterBuilder wbuilder;
    wbuilder["indentation"] = minimal ? "" : "  ";
    wbuilder["precision"] = 10;
    std::string document = Json::writeString(wbuilder, root);

    return document;
}

Json::Value JsonUtils::readFromFile(const std::string& filename) {

    Json::Value result;

    std::ifstream file1(filename);
    file1 >> result;

    return result;
}

void JsonUtils::writeToFile(const Json::Value& json, const std::string& filename) {
    std::ofstream out(filename);
    out << toString(json);
    out.close();
}
