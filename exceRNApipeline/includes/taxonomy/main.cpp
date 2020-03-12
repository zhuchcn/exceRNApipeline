#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "utilities.hpp"

using namespace std;

// template <typename IS>
// void catFile(string file) {
//     IS io (file.c_str());
//     if(!io) {
//         cout << "Can't open file" << endl;
//         return;
//     }
//     istream_iterator<Line> end;
//     istream_iterator<Line> it (io);

//     while(it != end) {
//         cout << *(it++) << endl;
//     }
// }

int main() {
    
    string s ("this\tis\ta\tstring");
    string* s2 = &s;
    string delim = "\t";
    vector<string> a = split(s, delim);

    return 0;
}