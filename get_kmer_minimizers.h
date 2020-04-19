////////////////////////////////////////////////////////////////////////////////
// get_kmer_minimizers.h
//   get kmer minimizer header file.
//
//  holding the get kmer_minimizer method which is a basic building block of the algorithm
//
////////////////////////////////////////////////////////////////////////////////
//  author: Alexander Petri

#ifndef GET_KMER_MINIMIZERS_H
#define GET_KMER_MINIMIZERS_H

#include "main.h"
#include "Minimizer.h"
#include "B-tree.hh"
#include "B_tree_node.hh"

using namespace std;
using namespace md;
//using namespace dyn;
typedef typename B_tree<int,int,7,3>::key_t _key_t;
typedef typename B_tree<int,int,7,3>::shifted_key_ptr_t _shifted_key_ptr_t;


/*
 * Print the current state of the forward list to the console (for debugging reasons)
 *
 * @param forward    the forward list to be printed to the console
 */
void print_forward_list(std::forward_list<std::string>& forward){
  ////cout<<"\n";
  for ( auto it = forward.begin(); it != forward.end(); ++it ){
    ////cout<<*it<<"\n";
  }
  ////cout<<"\n";
}

/*!
 * Generate the kmer minimizers of a sequence. Inspired by Kristoffer Sahlins' get_kmer_minimizer, however
 * this method uses a forward list to store the kmers. Does not generate end minimizers!!!
 *
 * @param sequence    the sequence for which minimizers are to be generated
 * @param k_size      the length of the window_kmers
 * @param w_size      the window size (length of the subsequence in which w kmers are present)
 *
 * @return minimizers  the minimizers for the sequence stored in a vector
 *
 * (@param w)         not a param of this function as w can be calculated by w=w_size-k_size+1
 */
std::vector<Minimizer> get_kmer_minimizers(string& sequence, int& k_size, int& w_size){
  int w = w_size - k_size+1;
  string max(k_size, 'Z');
  string curr_min;
  int min_pos=std::numeric_limits<int>::max();
  int min_diff=0;
  curr_min=max;
  std::vector<Minimizer> minimizers;
  std::forward_list<std::string> forward;
  auto beginIt=forward.begin();
  auto otherIt=beginIt;

  // generate the first minimizer of the sequence by finding the minimum kmer out of the first w k_mers
  for(int i=0;i<w;i++){
    //fill the forward list with the kmers
    std::string new_kmer=sequence.substr(i,k_size);
    if(forward.empty()){
      otherIt = forward.insert_after(forward.before_begin(), new_kmer);
    }
    else{
      otherIt = forward.insert_after(otherIt, new_kmer);
    }
  }
  int itnum=0;
  int startpos=0;
  std::string minimum=max;
  //find the minimal k_mer (the minimizer)
  for ( auto it = forward.begin(); it != forward.end(); ++it ){
    if(*it<minimum){
      minimum=*it;
      startpos=itnum;
    }
    itnum++;
  }
  //discard all kmers before the minimizer
  if(startpos!=0){
    for( auto it2 = 0; it2 != startpos; ++it2 ){
      forward.pop_front();
    }
  }
  forward.pop_front();
  //set the otherIt iterator to the last element of the forward list
  for( auto it = forward.begin(); it != forward.end(); ++it ){
    otherIt=it;
  }
  curr_min=minimum;
  //generate Minimizer object and add it to the solution
  Minimizer startmini=Minimizer(startpos,minimum);
  minimizers.push_back(startmini);
  Minimizer mini;
  min_pos=startpos;
  //find the rest of the Minimizers
  for (int i=w;i<sequence.length()-k_size+1;i++){
    string new_kmer=sequence.substr(i,k_size);
    //if the new k_mer is smaller than curr_min: We have found a minimizer! Empty forward list
    if(new_kmer<curr_min){
      curr_min=new_kmer;
      min_pos=i;
      mini.updateMinimizer(min_pos,curr_min);
      minimizers.push_back(mini);
      forward.clear();
    }
   //if we have not found a new minimizer in w consecutive kmers. Find minimum of the set of kmers->new minimizer
    else if(i-min_pos==w){
      otherIt = forward.insert_after(otherIt, new_kmer);
      int itnum=0;
      int pos=0;
      std::string minimum=max;
      for ( auto it = forward.begin(); it != forward.end(); ++it ){
        itnum++;
        if(*it<minimum){
          minimum=*it;
          pos=itnum;
        }
      }
      //only discard all kmers located before the minimizer from the forward list
      if(pos!=0){
        for( auto it2 = 1; it2 != pos; ++it2 ){
          forward.pop_front();
        }
        forward.pop_front();
      }
      for( auto it = forward.begin(); it != forward.end(); ++it ){
        otherIt=it;
      }
      curr_min=minimum;
      min_pos=min_pos+pos;
      mini.updateMinimizer(min_pos,curr_min);
      minimizers.push_back(mini);
    }
    //only add the new kmer to the forward list
    else{
      if(forward.empty()){
        otherIt = forward.insert_after(forward.before_begin(), new_kmer);
      }
      else{
        otherIt = forward.insert_after(otherIt, new_kmer);
      }
    }
  }
  return minimizers;
}
/*!
 * Generate the kmer minimizers of a sequence. Inspired by Kristoffer Sahlins' get_kmer_minimizer, however
 * this method uses a forward list to store the kmers. Does not generate end minimizers!!!
 *
 * @param sequence    the sequence for which minimizers are to be generated
 * @param k_size      the length of the window_kmers
 * @param w_size      the window size (length of the subsequence in which w kmers are present)
 *
 * @return minimizers  the minimizers for the sequence stored in a vector
 *
 * (@param w)         not a param of this function as w can be calculated by w=w_size-k_size+1
 */
std::vector<Minimizer> get_kmer_minimizers_algo(string& sequence, int& k_size, int& w_size,int& posshift){
  int w = w_size - k_size+1;
  string max(k_size, 'Z');
  string curr_min;
  int min_pos=std::numeric_limits<int>::max();
  int min_diff=0;
  curr_min=max;
  std::vector<Minimizer> minimizers;
  std::forward_list<std::string> forward;
  auto beginIt=forward.begin();
  auto otherIt=beginIt;
  int realpos=0;
  // generate the first minimizer of the sequence by finding the minimum kmer out of the first w k_mers
  for(int i=0;i<w;i++){
    //fill the forward list with the kmers
    std::string new_kmer=sequence.substr(i,k_size);
    if(forward.empty()){
      otherIt = forward.insert_after(forward.before_begin(), new_kmer);
    }
    else{
      otherIt = forward.insert_after(otherIt, new_kmer);
    }
  }
  int itnum=0;
  int startpos=0;
  std::string minimum=max;
  //find the minimal k_mer (the minimizer)
  for ( auto it = forward.begin(); it != forward.end(); ++it ){
    if(*it<minimum){
      minimum=*it;
      startpos=itnum;
    }
    itnum++;
  }
  //discard all kmers before the minimizer
  if(startpos!=0){
    for( auto it2 = 0; it2 != startpos; ++it2 ){
      forward.pop_front();
    }
  }
  forward.pop_front();
  //set the otherIt iterator to the last element of the forward list
  for( auto it = forward.begin(); it != forward.end(); ++it ){
    otherIt=it;
  }
  curr_min=minimum;
  //generate Minimizer object and add it to the solution
  realpos=startpos+posshift;
  //Minimizer startmini=Minimizer(realpos,minimum);
  minimizers.push_back(Minimizer(realpos,minimum));
  //minimizerTree->insert(realpos,minimum);
  //Minimizer mini;
  min_pos=startpos;
  //find the rest of the Minimizers
  for (int i=w;i<sequence.length()-k_size+1;i++){
    string new_kmer=sequence.substr(i,k_size);
    //if the new k_mer is smaller than curr_min: We have found a minimizer! Empty forward list
    if(new_kmer<curr_min){
      curr_min=new_kmer;
      min_pos=i;
      realpos=min_pos+posshift;
      //mini.updateMinimizer(realpos,curr_min);
      minimizers.push_back(Minimizer(realpos, curr_min));
      //minimizerTree->insert(realpos,curr_min);
      forward.clear();
    }
   //if we have not found a new minimizer in w consecutive kmers. Find minimum of the set of kmers->new minimizer
    else if(i-min_pos==w){
      otherIt = forward.insert_after(otherIt, new_kmer);
      int itnum=0;
      int pos=0;
      std::string minimum=max;
      for ( auto it = forward.begin(); it != forward.end(); ++it ){
        itnum++;
        if(*it<minimum){
          minimum=*it;
          pos=itnum;
        }
      }
      //only discard all kmers located before the minimizer from the forward list
      if(pos!=0){
        for( auto it2 = 1; it2 != pos; ++it2 ){
          forward.pop_front();
        }
        forward.pop_front();
      }
      for( auto it = forward.begin(); it != forward.end(); ++it ){
        otherIt=it;
      }
      curr_min=minimum;
      min_pos=min_pos+pos;
      realpos=min_pos+posshift;
      //mini.updateMinimizer(realpos,curr_min);

      minimizers.push_back(Minimizer(realpos,curr_min));
      //minimizerTree->insert(realpos,curr_min);
    }
    //only add the new kmer to the forward list
    else{
      if(forward.empty()){
        otherIt = forward.insert_after(forward.before_begin(), new_kmer);
      }
      else{
        otherIt = forward.insert_after(otherIt, new_kmer);
      }
    }
  }
  //cout<<"Printing getkmerminimizers\n";
  for(int i=0;i<minimizers.size();i++){
    minimizers[i].printMinimizer();
  }
  //cout<<"Printing getkmerminimizers done\n";
  return minimizers;
}

#endif
