////////////////////////////////////////////////////////////////////////////////
// dynamic_minimizer_no.h
//   Algorithm header file.
//
// Implementation of the dynamic minimizer algorithm
// This version does use a wt_string imported from the DYNAMIC library to store the sequence
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#ifndef DYNAMIC_MINIMIZER_NO_H
#define DYNAMIC_MINIMIZER_NO_H

#include "main.h"
#include "Variant.h"
#include "B-tree.hh"
#include "B_tree_node.hh"
#include "B_tree_operations.h"
#include "dynseq_functions.h"
#include "dynamic_minimizer.h"
#include "DYNAMIC-master/include/dynamic.hpp"

/*!
* Computes the upper (right) bound of the variation-impact-range.
* @param variants:          Vector of variants which are applied to the sequence
* @param variant index:     the index of the variant for which the variation-impact-range is calculated
* @param w_size:            window size for the minimizer generations
* @param k_size:            length of the k-mers
* @param sequence:          the sequence to be altered
*
* @return right:            the position of the upper bound of the variation-impact-range
* @return subseq:           boolean value, which is true if this variation-impact-range overlaps with the subsequent
*                           variation-impact-range and false if not

*/
std::tuple<int,bool> compute_right_bound_no_dynseq(std::vector<Variant>& variants,int& variant_index,int& w_size,int& k_size,std::string& sequence){
  Variant this_variant=variants.at(variant_index);
  bool subseq=false;
  int length=this_variant.getVariantLength();
  int originalseqlen=this_variant.getVariantOriginalSeqLen();
  int this_variant_pos=this_variant.getVariantPosition();
  int this_variant_delta=length - originalseqlen;
  string this_variant_seq=this_variant.getVariantSequence();
  int next_variant_pos=0;
  int right=0;
  if(variant_index+1<variants.size()){ //if this variant is not the last element of variants
    next_variant_pos=variants.at(variant_index+1).getVariantPosition(); //get the position of the next variation
    if(this_variant_pos+length+(2*w_size)-1+2*(k_size-1)>next_variant_pos){ //the next variation is in the variation-range of this variation
      right=this_variant_pos+originalseqlen;//set right to the position after the affected sequence part(last aff position)
      subseq=true;
    }
    else{ //the next variation is not in the variation range of this variant
      right = this_variant_pos + w_size + originalseqlen +(k_size - 2);
    }
  }
  else{//this variant is the last element of variants
    if(this_variant_pos+w_size+originalseqlen+k_size+1>=sequence.size()-1){
      right=sequence.size();
      cout<<"right! "<<right<<"\n";
    }
    else{
        right = this_variant_pos + w_size + originalseqlen + k_size+1;
    }
  }
  return make_tuple(right,subseq);
}
/*!
* Implementation of the dynamic minimizer algorithm.
* @param minimizerTree:     B-tree holding the final minimizers
* @param sequence:  the sequence to be altered
* @param variants:          Vector of variants which are applied to the sequence
* @param k_size:            length of the k-mers
* @param w_size:            window size for the minimizer generations
*/
std::string compute_dynamic_minimizers_no_dynseq(B_tree<int,std::string,7,3>* minimizerTree,std::string sequence,std::vector<Variant>& variants,int& k_size,int& w_size){
std::string dynamic_sequence=sequence;
int previous_shift=0;
int previous_right = 0;
int prevlength = 0;
int prevseqstart = 0;
bool prevseq = false;
bool subseq=false;
int shift=0;
std::vector<int> variants_in_subseq;
std::vector<int> variant_changes;
std::string previous_sequence="";
int var_impact_shift=0;
int appliedshift=0;
//iterate over all variations
  for(int i=0;i<variants.size();i++){
    int offset=0;
    std::tuple<int,int,int> left_infos;
    std::tuple<int,bool> right_infos;
    //apply the shift to the position of the current variation
    variants.at(i).updateVariantPosition(previous_shift);
    Variant this_var = variants.at(i);
    string this_variant_seq=this_var.getVariantSequence();
    int originalseqlen=this_var.getVariantOriginalSeqLen();
    int this_variant_delta= this_var.getVariantLength() - originalseqlen;
    shift=previous_shift+this_variant_delta;
    int variant_index = i;
    //compute the variation-impact-range
    left_infos = compute_left_bound(previous_right,this_var,prevlength,prevseqstart,k_size,w_size,prevseq);
    int left = std::get<0>(left_infos);
    offset=std::get<1>(left_infos);
    int thisstartpos=std::get<2>(left_infos);
    var_impact_shift+=this_variant_delta;
    //std::string whole_sequence=dynseq_tostring(dynamic_sequence);
    right_infos = compute_right_bound_no_dynseq(variants,variant_index,w_size,k_size,dynamic_sequence);
    int right = std::get<0>(right_infos);
    subseq = std::get<1>(right_infos);
    //cout<<"Subsequence: from "<<left<<" to "<<right<<":" "\n";
    //cout<<"Hello World before get_subseq\n";

    //get the substring covering the variation-impact-range
    cout<<"Dynseqsize: "<<dynamic_sequence.size()<<"\n";
    cout<<"Subsequence from "<<left<<" to "<<right<<"\n";
    std::string subsequence =dynamic_sequence.substr(left,right-left);
    //std::string subsequence = dynseq_get_substr(dynamic_sequence,left,right);
    //cout<<"Hello World after get_subseq\n";
    cout<<"Subsequence: "<<subsequence<<"\n";
    cout<<"Variation at position: "<<this_var.getVariantPosition()<<", original "<<originalseqlen<<",  new: "<<this_var.getVariantLength()<<"\n";
    //cout<<"Subsequence: from "<<left<<" to "<<right<<":"<<subsequence <<"\n";
    // apply the changes to the subsequence
    cout<<"give me range: "<<subsequence.size()<<"\n";
    cout<<"Right: "<<right<<"\n";
    cout<<"otherend: "<<offset+originalseqlen+1<<"\n";
    //apply the variation to the subsequence
    std::string newsubsequence="";
    if (offset > 0){
      int subseqend=offset+originalseqlen;
      if(subseqend>subsequence.size()){
      subseqend=subsequence.size();
      }
      newsubsequence =subsequence.substr(0,offset)+this_variant_seq+subsequence.substr(subseqend);

    }
    else{
      newsubsequence=subsequence.substr(0,originalseqlen)+this_variant_seq+subsequence.substr(originalseqlen+1);
    }
    subsequence=newsubsequence;
    //whole_sequence=dynseq_tostring(dynamic_sequence);
    cout<<"Subsequence after: "<<subsequence<<"\n";
    //cout<<"Sequence before update: "<<whole_sequence<<"\n";
    /*if(right+1<whole_sequence.size()){
      dynseq_update_substr(dynamic_sequence,left,right+1,subsequence);
    }
    else{
      dynseq_update_substr(dynamic_sequence,left,right,subsequence);
    }*/
    //update the sequence
    dynamic_sequence.replace(left, right-left,subsequence);
    //dynseq_update_substr(dynamic_sequence,left,right+1,subsequence);

    //whole_sequence=dynseq_tostring(dynamic_sequence);
    cout<<"Sequence: "<<dynamic_sequence<<"\n";
    //cout<<"Size of dynseq after"<<dynamic_sequence.size()<<"\n";
    //variants_in_subseq.push_back(v_pos);
    //variant_changes.push_back(originalseqlen);
    std::string fullsubseq="";
    //if this variation-impact-range overlaps with the previous, generate the full sub sequence for both in order to generate the minimizers for the merged subsequence
    if(prevseq){
      cout<<"prevseq\n";
      cout<<"prevseqstart: "<<prevseqstart<<", sprevseqsize"<<previous_sequence.size()<<", left: "<<left<<"\n";
      int overlap=prevseqstart+previous_sequence.size()-left;
      cout<<"overlap: "<<overlap <<"\n";
      if (overlap<0){
        fullsubseq=previous_sequence+subsequence;
      }
      else{
        fullsubseq=previous_sequence+subsequence.substr(overlap,subsequence.size()-overlap);
      }
    }
    else{
      fullsubseq=subsequence;
    }
    cout<<"Fullsubseq: "<<fullsubseq<<"with size "<<fullsubseq.size()<<"\n";
    //if this is the last variation in the merged impact range;
    if(!subseq){
      //left=left-shift+appliedshift;
      //right=right-shift+appliedshift;
      //cout<<"I'm going to delete the minimizers in the set ["<<thisstartpos<<" ,"<<right<<" ]\n";
      //int suc=minimizerTree->successor(right).key->value;
      //update_minimizerTree(minimizerTree, thisstartpos, right);
      //cout<<"Succ:"<<suc<<"\n";
      //print_minimizerTree(minimizerTree);
      //cout<<"Printing the minimizer Tree done\n";
      //cout<<"Varimpact "<<var_impact_shift<<"\n";
      //cout<<"updating the minimizerTree\n";

      //update the minimizer tree holding the minimizers
      update_minimizerTree(minimizerTree,fullsubseq,thisstartpos,k_size,w_size,var_impact_shift);

      //cout<<"done with applying shifts\n";
      //get_kmer_minimizers_algo(minimizerTree,fullsubseq,k_size,w_size,thisstartpos);
      cout<<"Printing the minimizer Tree\n";
      print_minimizerTree(minimizerTree);
      cout<<"Printing the minimizer Tree done\n";
      appliedshift+=var_impact_shift;
      var_impact_shift=0;
      prevseq=false;
    }
    else{
      previous_sequence=fullsubseq;
      prevseqstart=thisstartpos;
    }
   cout<<"Previous shift: "<<shift<<"\n";
    previous_shift=shift;
    previous_right=right;
    prevlength=this_variant_delta;
    prevseq=subseq;
    cout<<"Full sequence: "<<dynamic_sequence<<"\n";
    //cout<<"Full sequence: "<<dynseq_tostring(dynamic_sequence)<<"\n";
    //return dynamic_sequence;
  }
  cout<<"Algorithm finished!!!\n";
  return dynamic_sequence;
}
#endif
