#include<iostream>
#include<fstream>
#include<cstring>
#include<cstdint>
#include "using.h"
#include "../include/Awk.h"
#include "Dict.h"
#include "stats.h"
#include "util.h"
#include "GetOpt.h"
//#include "../include/gzstream.h"
void convert(char* p, char* q){q = p;}
void convert(char*p, string& s){s = string(p);}

template<typename KEY>
void test(int nkeys, const string& test_file){
  // now a big test with every word from file test_file
  
  Dict<KEY, int> d1(nkeys); // set to a small size to test re-sizing
  Array<KEY> all_keys(nkeys);
  Awk awk;
  if(!awk.open(test_file.c_str())){
    cerr<< "Can't open "<< test_file <<endl;
    exit(1);
  }
  int i(0), n(0);
  while(!awk.eof()){
    awk.next();
    for(int j = 1;j <= awk.nf;j++){
      KEY key;
      
      convert(const_cast<char*>(awk[j]),key);
      d1[key]= i + j;
      all_keys[n++] = key;  // save the string for a subsequent check
      cout << "adding "<<key<<" at "<<i+j<<endl;
      if(!d1.has_key(all_keys[n-1])){
        cout<< "missing "<<all_keys[n-1] <<"  at "<<i+j<<endl;
        exit(1);
      }
    }
    i += awk.nf;
  }
  cout << i << " keys were added, we now have " << d1.nkeys()<<endl;
  for(int i = 0;i < n;i++){
    if(d1.has_key(all_keys[i])) continue;
    else{
      cout << "missing "<<all_keys[i]<< " at i = "<<i<<endl;
    }
  }
  cout << "d1 has "<<d1.nkeys()<< " keys\n";
}

int main(int argc, char** argv){
  int nkeys(10000); // blocksize for the Dict Arrays. They will automagically resize as needed.
  string key_type = "string";
  string test_file = "Dict.h";

  GetOpt cl(argc,argv);
  cl.get("nkeys", nkeys);
  cl.get("key_type",key_type);
  cl.get("test_file",test_file);

  if(key_type == "char*") test<char*>(nkeys,test_file);
  else if (key_type == "string") test<string>(nkeys,test_file);
}
