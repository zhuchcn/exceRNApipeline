#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

// helper functions to process strings
// https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}
inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}
inline std::vector<std::string> split(const std::string& s, std::string delim) {
    std::vector<std::string> vec;
    std::string token;
    size_t start = 0;
    size_t end = s.find(delim);
    while(end != std::string::npos) {
        token = s.substr(start, end - start);
        trim(token);
        vec.push_back(token);
        start = end + delim.length();
        end = s.find(delim, start);
    }
    token = s.substr(start, end - start);
    vec.push_back(token);
    return vec;
}

// Checks if a string s starts with a pattern p.
inline bool startswith(const std::string &s, const std::string &p) {
    if(s.length() < p.length()) return false;
    return s.find(p) == 0;
}

// Checks if a string s ends with a pattern p.
inline bool endswith(const std::string &s, const std::string &p) {
    if(s.length() < p.length()) return false;
    return s.compare(s.size() - p.size(), s.size(), p) == 0;
}

inline void printVectorString(std::vector<std::string>& v){
    std::cout << "[ ";
    for(auto vit = v.begin(); vit < v.end(); ++vit) {
        std::cout << "\"" << *vit << (vit < v.end()-1 ? "\", " : "\" ]\n");
    } 
}

// helper class for generating istream_iterator while reading from a text file.
class Line : public std::string {
public:
    friend std::istream &operator>>(std::istream &is, Line &l){
        std::getline(is, l);
        return is;
    }
};

#endif