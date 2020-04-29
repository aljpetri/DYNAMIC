////////////////////////////////////////////////////////////////////////////////
// generate_random_test_cases.h
//
//   header file including functions, which generate the random sequences and variations.
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#ifndef GEN_RANDOM_TEST_CASES_H
#define GEN_RANDOM_TEST_CASES_H

#include "main.h"
#include "Variant.h"


/*
* Generates a random integer with a value between left and right
*
* @param left     the lower bound for the integer
* @param right    the upper bound for the integer
*
* @return rand_int   the random integer
*/
int generate_random_integer_bounded(int& left,int& right){
  auto e = std::mt19937{};
  //seed engine with system clock
  e.seed(std::chrono::system_clock::now().time_since_epoch().count());
  auto distr = std::uniform_int_distribution<int>{left, right};
  int rand_int=distr(e);
  return rand_int;
}
/*
*Generates a string holding a random sequence of nucleotides (A,T,G,C)
*
* @param length     the length of the to be generated sequence
*
* @return sequence  the random sequence

*/
string generate_random_sequence(int length){
  //std::default_random_engine generator;
  char* temp=new char[length];
  for(int j=0;j<length;j++){
    int left=0,right=3;
    int basenr=generate_random_integer_bounded(left,right);
    char base{};
    switch(basenr){
      case 0:
        base='A';
        break;
      case 1:
        base='C';
        break;
      case 2:
        base='G';
        break;
      case 3:
        base='T';
        break;
      default:
        base='Z';
        break;
    }
    temp[j]=base;
  }
  std::string sequence=temp;
return sequence;
}

/*
* Function, which generates random pseudo vcf entries
*
* @param sequence   the sequence for which the variants are generate_random_sequence
* @param number_of_variations   the number of variants to be generated
*
* @return variants  a vector of variant objects (position,sequence,originalseqlen,length)
*/
vector<Variant> generate_random_variations(string& sequence, int& num_variations){
  std::vector<Variant> variants;
  int prev_variant_end=1;
  int factor=sequence.length()/num_variations;
  int zero=0;
  int seqlenth=sequence.length()-1;
  for(int i=1;i<=num_variations;i++){
      int right=i*factor;
      int position = generate_random_integer_bounded(prev_variant_end,right);
      int original_var_length=0;
      int max=20;
      int possible_length=right-position;
      original_var_length=generate_random_integer_bounded(zero, possible_length);
      int new_var_length = generate_random_integer_bounded(zero,max);
      std::string variant_seq=generate_random_sequence(new_var_length);
      prev_variant_end=position+original_var_length+1;
      Variant this_variant=Variant(position, original_var_length, new_var_length, variant_seq);
      variants.push_back(this_variant);
  }
return variants;
}


#endif

