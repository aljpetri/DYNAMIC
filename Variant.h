////////////////////////////////////////////////////////////////////////////////
// Variant.h
//   variant class header file.
//
//  class to internally represent a variant
//
////////////////////////////////////////////////////////////////////////////////
//  author: Alexander Petri

#ifndef VARIANT_H
#define VARIANT_H

#include "main.h"
/*
* Class to define a variant
* author: Alexander Petri
*
* @param position    the position of the variant
* @param sequence    the sequence of the variants
* @param originalseqlen  the length of the subsequence before the variant was applied
* @param length      the length of the subsequence after the variant was applied
*
*/
class Variant{
private:
  int position;
  string sequence;
  int originalseqlen;
  int length;
public:
  // Constructor
  Variant(int& pos, int& origin, int& len, string& seq){
    position=pos;
    sequence=seq;
    originalseqlen=origin;
    length=len;
  }

  /*
  *
  * updates the position of the variant according to shift
  *
  * @param shift   the amount of bases the variant is shifted by
  *
  */
  void updateVariantPosition(int& shift){
    position+=shift;
  }
  /*
  *returns the position of the variant
  */
  int getVariantPosition(){
    return position;
  }
  /*
  *returns the original length of the variant
  */
  int getVariantOriginalSeqLen(){
    return originalseqlen;
  }
  /*
  *returns the new length of the variant
  */
  int getVariantLength(){
    return length;
  }
  /*
  *returns the new sequence of the variant
  */
  std::string getVariantSequence(){
    return sequence;
  }
  /*
  * prints the variant to the console
  *
  *Output: Variant at position, original: originalsequencelength, new: newlength, sequence newsequence
  */
  void printVariant(){
    cout<<"Variant at "<<position<<" , original: "<<originalseqlen<<", new: "<< length<<", sequence "<<sequence<<"\n";
  }
};

#endif
