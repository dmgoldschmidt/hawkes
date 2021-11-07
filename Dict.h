#ifndef Dict_H
#define Dict_H
#include "Array.h"
#include "libcrc/include/checksum.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include "stats.h"

using namespace std;

/* Specialized high-perfomance associative array.  Usage:

Dict<char*, int> d;
d["one"] = 1;
d["two"] = 2;
int sum = d["one"] + d["two"];

if(d.has("one"){ (do something)}
for(inti = 0;i < d.nkeys();i++){
(do something with or to d(i))
} 

 */

struct Hash {
  uint64_t h0;
  uint64_t h1;
  uint64_t& operator[](int i){ return i==0? h0 : h1;}
  bool operator==(const Hash& h) const {
    if(h0 != h.h0 || h1 != h.h1) return false;
    return true;
  }
};


template <typename KEY, typename ITEM>
class Dict{
  Array<int> hash_table;
  Array<Hash> key_hash;
  Array<ITEM> item;
  int _nkeys; // no. of distinct keys seen so far  
  LFSR128 my_rand;
  int size; //block size for hash table
  double rehash_at; // max fraction of hash to fill before expanding size and re-hashing

  void rehash(int new_size){ // private function to initialize and/or expand the size of the hash table
    for(int i = 0;i < new_size;i++) hash_table[i] = -1;
    if(new_size < 1.1*_nkeys) throw "Dict internal error";
        for(int i = 0;i <_nkeys;i++){
      //      uint64_t* p = key_hash[i].p;
      my_rand.step(&key_hash[i][0]); // initialize LFSR at the next key hash 
      while(true){
        int j = (int)(my_rand.uniform()*new_size);
        if(hash_table[j] == -1){
          hash_table[j] = i;
          break;
        }
      }
    }
  }

  Hash hash_it(const KEY& key){
    Hash hash;
    unsigned char key_bytes[sizeof(key)+1];
    key_bytes[sizeof(key)] = '1';  // append a constant character for the second crc call
    const unsigned char* p = &key;
    for(int i = 0;i < sizeof(key);i++) key_bytes[i] = p[i]; // copy in the key bytes
    hash[0] = crc_64_ecma(key_bytes,sizeof(key));
    hash[1] = crc_64_ecma(key_bytes,sizeof(key)+1); // 2nd hash appends '1', making a 128-bit hash.
    return hash;
  }

  int find(Hash hash){ // private function to seach the hash table
    if(_nkeys >= rehash_at*size) rehash(size *= 2);   
    my_rand.step(&hash[0]); // re-seed the LFSR with the key hash
    while(true){ // follow the "trail" until either the key is found or an empty slot is found
      uint hti = (uint)(my_rand.uniform()*size);
      if(hash_table[hti] == -1 || key_hash[hash_table[hti]] == hash) return hti;
    }
  }


public:
  Dict(void){}
  Dict(int s = 10000, double r = .1) : size(s), rehash_at(r), key_hash(s), item(s), _nkeys(0){
    if constexpr(std::is_pointer<KEY>::value) throw "Keys cannot be of pointer type\n";
    rehash(size); // initialize the hash table
  }
  ~Dict(void){}
  ITEM& operator[](const KEY& key){
    int i, hti;
    Hash hash = hash_it(key); // hash the key with crc64 (twice!)
    i = hash_table[hti = find(hash)];
    if(i >= 0) return item[i]; // key exists, so return the ref
    else{ // it's new, so add a slot
      key_hash[_nkeys] = hash;  // save the key hash for later comparisons
      hash_table[hti] = _nkeys; // update the hash table
      return item[_nkeys++];
    }
  }
  ITEM& operator()(int i){ // access items sequentially
    if(i < 0 || i >= _nkeys) throw "Bad index value\n";
    return item[i];
  }
  bool has(const KEY& key){
    Hash hash = hash_it(key);
    if(hash_table[find(hash)] >= 0)return true;
    return false;
  }
  int nkeys(void){return _nkeys;}
};

template<typename ITEM>
int Dict<char*, ITEM>::find(Hash hash){
  




#endif
