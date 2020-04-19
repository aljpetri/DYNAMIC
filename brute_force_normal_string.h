////////////////////////////////////////////////////////////////////////////////
// brute_force_normal_string.h
//   Algorithm header file.
//
// brute force implementation used as a reference for the dynamic minimizer algorithm.
// This version does use a std::string to store the sequence
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#ifndef BRUTE_FORCE_NO_H
#define BRUTE_FORCE_NO_H

#include "main.h"
#include "Variant.h"
#include "B-tree.hh"
#include "dynamic_minimizer.h"
#include "DYNAMIC-master/include/dynamic.hpp"


/*!
* Brute force implementation to compute dynamic minimizers.
* @param minimizerTree:     B-tree holding the final minimizers
* @param dynamic_sequence:  the sequence to be altered
* @param variants:          Vector of variants which are applied to the sequence
* @param k_size:            length of the k-mers
* @param w_size:            window size for the minimizer generations
*/
void brute_force_minimizer_computation_normal_string(B_tree<int,std::string,7,3>* minimizerTree, std::string& dynamic_sequence,std::vector<Variant>& variants,int& k_size,int& w_size){
  int shift=0;
  int original_lenth=dynamic_sequence.size();
  //iterate over the variants to be applied to the sequence (This does only update the sequence)
  for(int i=0;i<variants.size();i++){
    int delta=0;
    int v_pos=variants[i].getVariantPosition();
    int v_originalseqlen=variants[i].getVariantOriginalSeqLen();
    std::string v_sequence=variants[i].getVariantSequence();
    std::string newsubsequence="";
    if(v_pos>0){
      int start=v_pos+shift;
      int end=start+v_originalseqlen;
      if(end>dynamic_sequence.size()){
        end=dynamic_sequence.size();
      }

      newsubsequence =dynamic_sequence.substr(0,start)+v_sequence+dynamic_sequence.substr(end);
    }
    else{
      newsubsequence =dynamic_sequence.substr(0,1)+v_sequence+dynamic_sequence.substr(v_originalseqlen+1);

    }
    dynamic_sequence=newsubsequence;
    delta=variants[i].getVariantLength()-v_originalseqlen;
    shift+=delta;

  }

  //std::string finseq=dynseq_tostring(dynamic_sequence);
  cout<<"Final BF-Sequence no dynstring: "<<dynamic_sequence<<"\n";
  //generate the minimizers for the updated sequence
  std::vector<Minimizer> minimizers =get_kmer_minimizers(dynamic_sequence,k_size,w_size);
  //fill the minimizers into the minimizer tree
  fill_minimizer_tree(minimizerTree,minimizers);
}

#endif
