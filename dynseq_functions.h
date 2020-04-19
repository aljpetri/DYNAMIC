////////////////////////////////////////////////////////////////////////////////
// brute_force.h
//   Algorithm header file.
//
// brute force implementation used as a reference for the dynamic minimizer algorithm.
// This version does use a wt_string imported from the DYNAMIC library to store the sequence
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#ifndef DYNAMIC_SEQUENCE_FUNCTIONS_H
#define DYNAMIC_SEQUENCE_FUNCTIONS_H

#include "main.h"
#include "Variant.h"
#include "B-tree.hh"
#include "include/dynamic.hpp"

/*
* delivers a substring from the dynamic sequence
* @param dynamic_sequence   the dynamic sequence
* @param left               the lower bound for the substring
* @param right              the upper bound for the substring
*
* @return subsequence       the subsequence
*/
std::string dynseq_get_substr(dyn::wt_str dynamic_sequence, int left, int right){
  std::string subsequence="";
  for(int i=left;i<=right;i++){
    subsequence+=dynamic_sequence.at(i);
  }
  return subsequence;
}
/*
* translates the dynamic sequence into a std::string
* @param dynamic_sequence   the dynamic sequence
*
* @return output            the std::string
*/
std::string dynseq_tostring(dyn::wt_str& dynamic_sequence){
  int size=dynamic_sequence.size();
  std::string output="";
  for(int i=0;i<size;i++){
    output+=dynamic_sequence.at(i);
  }
  return output;
}
/*
* updates the dynamic sequence by replacing a substring
* @param dynamic_sequence   the dynamic sequence
* @param left               the lower bound for the elements to be deleted
* @param right              the upper bound for the elements to be deleted
* @param subsequence        the subsequence which is inserted into the dynamic sequence
*/
void dynseq_update_substr(dyn::wt_str& dynamic_sequence, int left, int right,std::string subsequence){
  //cout<<"Size of dynseq:"<<dynamic_sequence.size()<<"\n";
  for(int i=left;i<right;i++){
    //if(i<dynamic_sequence.size()){
      dynamic_sequence.remove(left);
    //  cout<<"Dynseq after: "<<dynseq_tostring(dynamic_sequence)<<"\n";
  //  }
  }
  for(int i=0;i<subsequence.length();i++){
    //cout<<"Adding Element at "<<left+i<<"\n";
    dynamic_sequence.insert(left+i,subsequence.at(i));
    //cout<<"Subsequence at i:"<<subsequence.at(i)<<"\n";
    //cout<<"Dynseq after: "<<dynseq_tostring(dynamic_sequence)<<"\n";
    // use vector<char>in order to not fall out of range
  }
}




#endif
