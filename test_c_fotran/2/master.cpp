#include <iostream>
using namespace std;

extern "C" {
    void square_(int *n, int *out);  
}

void main(){  
   int *input;  
   *input =3;  
   int *output;  
   square_(input,output);  
   cout *output << endl;  // returns 9  
}

