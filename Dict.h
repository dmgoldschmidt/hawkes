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

if(d.has_key("one"){ (do something)}
for(int i = 0;i < d.nkeys();i++){
   (do something with or to d(i))
 } 

 */

// template<typename KEY>
// char* get_bytes(const KEY& key, int& nbytes);

// template<>
// char* get_bytes(const char* key, int& nbytes);

// template<typename KEY>
// char* get_bytes(const KEY& key, int& nbytes){
//   //  if constexpr(std::is_pointer<KEY>::value) throw "Keys cannot be of pointer type\n";
//   char* p = (char*) &key;
//   nbytes = sizeof(key);
//   return p;
// }

//template<>
const char* get_bytes(const char* key, int& nbytes){
  const char* p = key;
  nbytes = strlen(key);
  return p;
}

const char* get_bytes(const string& key, int& nbytes){
  const char* p = key.c_str();
  nbytes = key.size();
  return p;
}

struct my_Hash {
  uint64_t h0;
  uint64_t h1;
  uint64_t& operator[](int i){ return i==0? h0 : h1;}
  bool operator==(const my_Hash& h) const {
    if(h0 != h.h0 || h1 != h.h1) return false;
    return true;
  }
};
  
template <typename KEY, typename ITEM>
struct Dict{
  Array<int> hash_table; // look up the key hashes, return an index into Array key_hash
  Array<my_Hash> key_hash; /* hashes of all keys seen so far NOTE: we are assuming no collisions
                         * for the 128-bit hash.  So we don't need to store copies of the actual keys.
                         * If this doesn't work out, we can always store the keys themselves.
                         */
  Array<ITEM> item; // storage slots (1 per key) user stores/modifies/retrieves his items here
  int _nkeys; // no. of distinct keys seen so far  
  LFSR128 my_rand; // 128 bit generator (uses the xorshift128+ algorithm)
  int size; //current size for hash table
  double rehash_at; // max fraction of hash table to fill before expanding size

  void rehash(int new_size){ // private function to initialize and/or expand the size of the hash table
    cout << "resizing to "<<new_size<<" with nkeys = "<<_nkeys<<endl;
    size = new_size;
    for(int i = 0;i < new_size;i++) hash_table[i] = -1; // hash table is now empty
    for(int i = 0;i <_nkeys;i++){
      my_rand.reseed(&key_hash[i][0]); // initialize LFSR at the next key hash 
      while(true){ // store the key_hash indices in the longer hash table (it automagically expands)
        int j = (int)(my_rand.uniform()*new_size);
        if(hash_table[j] == -1){
          hash_table[j] = i;
          break;
        }
      }
    }
  }

  my_Hash hash_it(const KEY& key){ 
    my_Hash my_hash;
    const char* p;
    int nbytes;

    p = get_bytes(key,nbytes);
    unsigned char key_bytes[nbytes+1];
    key_bytes[nbytes] = '1';  // append a constant character for the second crc call
    // cout << "nbytes = "<<nbytes<<endl;
    for(int i = 0;i < nbytes;i++){
      key_bytes[i] = p[i]; // copy in the key bytes
    //   cout << format("key_byte[%d] = %d = %c\n",i,key_bytes[i],key_bytes[i]);
    }
    my_hash[0] = crc_64_ecma(key_bytes,nbytes);
    my_hash[1] = crc_64_ecma(key_bytes,nbytes+1); // 2nd hash appends '1', making a 128-bit hash.
    return my_hash;
  }

  int find(const KEY& key, my_Hash& my_hash){ // private function to seach the hash table
    if(_nkeys >= rehash_at*size) rehash(2*size);   
    my_hash = hash_it(key);
    my_rand.reseed(&my_hash[0]); // re-seed the LFSR with the key hash
    for(int i = 0;i < size;i++){ // follow the "trail" until either the key is found or an empty slot is found
      uint hti = (uint)(my_rand.uniform()*size);
      if(hash_table[hti] == -1 || key_hash[hash_table[hti]] == my_hash) return hti;
    }
    throw "Hash table full"; // not found
  }


public:
  Dict(int s = 10000, double r = .9) : size(s), rehash_at(r), key_hash(s), hash_table(s), item(s), _nkeys(0){
    rehash(size); // initialize the hash table
  }
  ~Dict(void){}

  bool has_key(const KEY& key){
    my_Hash my_hash;
    return hash_table[find(key,my_hash)] >= 0;
  }
  //   my_Hash hash = hash_it(key);
  //   if(hash_table[find[hash]] >= 0)return true;
  //   return false;
  // }

  ITEM& operator[](const KEY& key){
    int i;
    uint hti;
    my_Hash my_hash;
    my_hash = hash_it(key);  // hash the key
    my_rand.reseed(&my_hash[0]); // re-seed the LFSR with the key hash
    hti = find(key,my_hash); // returns an index into hash_table
    i = hash_table[hti]; // now we have an index into both key_hash and item
    if(i  >= 0) return item[i]; // key exists, so return the ref
    else{ // it's new, so add a slot
      key_hash[_nkeys] = my_hash;  // save the key hash for later comparisons
      hash_table[hti] = _nkeys; // update the hash table
      return item[_nkeys++];
    }
  }
  ITEM& operator()(int i){ // access items sequentially
    if(i < 0 || i >= _nkeys) throw "Bad index value\n";
    return item[i];
  }

  int nkeys(void){return _nkeys;}
};





#endif
