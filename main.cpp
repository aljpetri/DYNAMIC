////////////////////////////////////////////////////////////////////////////////
// get_kmer_minimizers.h
//   get kmer minimizer header file.
//
//  holding the get kmer_minimizer method which is a basic building block of the algorithm
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#include "Minimizer.h"
#include "main.h"
#include "B-tree.hh"
#include "B_tree_node.hh"
#include "generate_random_test_cases.h"
//#include "Alter_Minimizers_Algo.h"
#include "get_kmer_minimizers.h"
#include "brute_force_normal_string.h"
#include "B_tree_operations.h"
#include "dynamic_minimizer.h"
#include "dynamic_minimizer_no.h"
#include "brute_force.h"
#include "dynseq_functions.h"
#include "include/dynamic.hpp"

using namespace std;
using namespace md;
using namespace dyn;
typedef typename B_tree<int,int,7,3>::key_t _key_t;
typedef typename B_tree<int,int,7,3>::shifted_key_ptr_t _shifted_key_ptr_t;

//typedef typename B_tree<int,int,3,1>::key_t _key_t;
//typedef typename B_tree<int,int,3,1>::shifted_key_ptr_t _shifted_key_ptr_t;

int main(){
  //Predefined sequences for debugging reasons only
  //auto sequence="CCCAACCCGGCGCGGCCAAGAAGCCAGCCAGGCGAAAACAAGGAGGCGAGAGCACCGCGCGAACCGGACGCGGCCCCCCAAAAAAGAAGAAACACGAAGA"s;
  //auto sequence= "TCAACGGTCTCTGAGCGTCAACCTCGTACTTAGAAGGGCGGAACCGCCAGCGTGCCTACTCCAGTCGTCGATTTACATTAACATACGTTCTCAGCTCTAA"s;
  //auto sequence= "ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA"s;
  //auto sequence="231032101233101"s;
  //main parameters needed for the algorithm
  int k = 4;
  int w= 6;
  int seqlen=100;
  int numbervars=5;

  string sequence=generate_random_sequence(seqlen);
  std::string sequence2=sequence;
  //cout<<"Random sequence: "<<sequence<<"\n";
  const uint64_t sigma=4;
  const uint64_t sigma2=0;
  wt_str dynamic_sequence(sigma);
  wt_str dynamic_sequence2(sigma);
  //sstd::string finalsequence="";
  //dynamic_sequence.insert(0,'a');
  //dynamic_sequence.
  std::vector<char> myVector(sequence.begin(), sequence.end());
  for (int i=0;i<myVector.size();i++){
    dynamic_sequence.push_back(myVector[i]);
  }
  for (int i=0;i<myVector.size();i++){
    dynamic_sequence2.push_back(myVector[i]);
  }
  cout<<"Size: "<<dynamic_sequence.size()<<"\n";
  cout<<"aSize: "<<dynamic_sequence.alphabet_size()<<"\n";
  //std::string substr=dynseq_get_substr(dynamic_sequence,0,3);
  //cout<<"Substring: "<<substr<<"\n";
//  cout<<"Dynseq before: "<<dynseq_tostring(dynamic_sequence)<<"\n";

  vector<Variant> variants=generate_random_variations(sequence,numbervars);
  cout<<"Random variations generated\n";
  vector<Variant> variants2=variants;
  vector<Variant> variants3=variants;
  /*vector<Variant> variants;
  int pos=6;
  int origlen=1;
  int len=1;
  std::string vsequence="A";
  Variant this_variant=Variant(pos,origlen,len,vsequence);
  variants.push_back(this_variant);
  pos=15;
  origlen=4;
  len=9;
  vsequence="CGCAGCGAC";
  Variant v_0=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_0);
  pos=34;
  origlen=3;
  len=4;
  vsequence="ACCA";
  Variant v_1=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_1);
  pos=40;
  origlen=3;
  len=6;
  vsequence="AGAGCG";
  Variant v_2=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_2);
  pos=62;
  origlen=2;
  len=9;
  vsequence="GAACACACA";
  Variant v_5=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_5);
  pos=68;
  origlen=1;
  len=9;
  vsequence="GGCGAAAGA";
  Variant v_3=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_3);
  pos=82;
  origlen=2;
  len=12;
  vsequence="CACCGGGCGAGC";
  Variant v_4=Variant(pos,origlen,len,vsequence);
  variants.push_back(v_4);
  /*6: 1 1 A
  15: 4 9 CGCAGCGAC
  34: 3 4 ACCA
  40: 3 6 AGAGCG
  62: 2 9 GAACACACA
  68: 1 9 GGCGAAAGA
  82: 2 12 CACCGGGCGAGC
*/

//vector<Minimizer> minimizers;
/*for(int i=1;i<21;i++){
  int pos=i;
  std::string seq="AAAAA";
  Minimizer m0=Minimizer (pos,seq);
  minimizers.push_back(m0);
}*/
  //cout<<"random variants generated from random sequence \n";
  auto begin = chrono::high_resolution_clock::now();
  vector<Minimizer> minimizers = get_kmer_minimizers(sequence,k,w);
  auto sndtime=std::chrono::system_clock::now();
  auto dur=sndtime-begin;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
  cout<<"Time needed: "<< ms<<"miliseconds\n";

  cout<<"Random variants:\n";
  for(int i=0;i<variants.size();i++){
    variants.at(i).printVariant();
  }
  //for (int i=0;i< variants.size();i++){
  //  variants[i].printVariant();
  //}


  //write the minimizers to the command line
  /*for(int i=0;i<minimizers.size();i++){
    minimizers.at(i).printMinimizer();
  }*/

  //int val=3;
  //cout<<sequence<<"\n";
  //vector<Minimizer> minis=get_kmer_minimizers(sequence,k,w);
  B_tree<int,std::string,7,3>* minimizerTree = new B_tree<int,std::string,7,3>();
  B_tree<int,std::string,7,3>* minimizerTreeAlgo2 = new B_tree<int,std::string,7,3>();
  //generates the B-tree holding the minimizers generated above and fills it
  //B_tree<int,std::pair<std::string,int>,7,3>* minimizerTree = new B_tree<int,std::pair<std::string,int>,7,3>();
  //B_tree<int,std::pair<std::string,int>,7,3>* minimizerTree = new B_tree<int,std::pair<std::string,int>,7,3>();
  //fill_minimizer_tree(minimizerTree,minis);
  //int posshift=0;
  //std::vector<Minimizer> minis=get_kmer_minimizers_algo(sequence,k,w,posshift);
  fill_minimizer_tree(minimizerTree,minimizers);
  fill_minimizer_tree(minimizerTreeAlgo2,minimizers);
  print_minimizerTree(minimizerTree);
  //cout<<"Dynseq before: "<<dynseq_tostring(dynamic_sequence)<<"\n";
  //dynseq_update_substr(dynamic_sequence, 0, 4,"AAAAAAAAAAAAA");
  //cout<<"Dynseq after: "<<dynseq_tostring(dynamic_sequence)<<"\n";
  /*int elem=9;
  int shift=4;
  minimizerTree->search(elem);
  minimizerTree->shift_greater(elem,shift);
  print_minimizerTree(minimizerTree);
*/

  B_tree<int,std::string,7,3>* minimizerTreeBF = new B_tree<int,std::string,7,3>();
  //compute_dynamic_minimizers(minimizerTree,dynamic_sequence,variants,k,w);
  cout<<"Starting normal compute dynamic minimizers\n";
  auto begin3 = chrono::high_resolution_clock::now();

  compute_dynamic_minimizers(minimizerTree,dynamic_sequence2,variants,k,w);
  auto sndtime3=std::chrono::system_clock::now();
  auto dur3=sndtime3-begin3;
  auto msalgo = std::chrono::duration_cast<std::chrono::milliseconds>(dur3).count();
  cout<<"Starting compute dynamic minimizers without dynseq\n";
  auto begin2no = chrono::high_resolution_clock::now();

  sequence=compute_dynamic_minimizers_no_dynseq(minimizerTreeAlgo2,sequence,variants2,k,w);
  auto sndtime2no=std::chrono::system_clock::now();
  auto durno=sndtime2no-begin2no;
  auto msalgono = std::chrono::duration_cast<std::chrono::milliseconds>(durno).count();
  cout<<"Starting brute force\n";
  auto begin2 = chrono::high_resolution_clock::now();
  brute_force_minimizer_computation(minimizerTreeBF,dynamic_sequence,variants3,k,w);
  auto sndtime2=std::chrono::system_clock::now();
  auto dur2=sndtime2-begin2;
  auto msbf = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
  std::vector<Minimizer> newminisbf=minimizer_to_vector(minimizerTreeBF);
  cout<<"Main Hello World!\n";
  //std::vector<Minimizer> newminisalgono=minimizer_to_vector(minimizerTreeAlgo2);
  cout<<"Main Hello World2!\n";
  std::string bf_result=dynseq_tostring(dynamic_sequence);
  cout<<"Main Hello World\n";
  //*print_minimizerTree(minimizerTree);
  std::string algo_result=dynseq_tostring(dynamic_sequence2);
  //for(int i=0;i<newminisbf.size();i++){
  //  Minimizer bfmini=newminisbf[i];
  //  cout<<"Minimizer "<<bfmini.getSequence()<<": "<<bfmini.getPosition()<<"\n";
  //}
  cout<<"Algo: "<<algo_result<<"\n";
  cout<<"Algo no dynseq: "<<sequence<<"\n";
  if(bf_result.compare(algo_result)==0){
    cout<<"The algorithm returned the right sequence!\n";
  }
  if(bf_result.compare(sequence)==0){
    cout<<"The algorithm no dynseq returned the right sequence!\n";
  }
  std::vector<Minimizer> algominis=minimizer_to_vector(minimizerTree);
  cout<<"Main Hello World algominis\n";
  cout<<"Bf-Minimizer      vs        AlgoMinimizer\n";
  bool rightMinis=true;
  bool rightMinisno=true;
  for(int i=0;i<newminisbf.size();i++){
    Minimizer bfmini=newminisbf[i];
    if(i<=algominis.size()){
      Minimizer algomini=algominis[i];
      cout<<"Minimizer "<<bfmini.getSequence()<<": "<<bfmini.getPosition()<<"  vs     "<< algomini.getSequence()<<": "<<algomini.getPosition()<<"   "<<(bfmini.getSequence()==algomini.getSequence() && bfmini.getPosition()==algomini.getPosition())<<"\n";
      if(!(bfmini.getSequence()==algomini.getSequence() && bfmini.getPosition()==algomini.getPosition())){
        rightMinis=false;
        //cout<<"here wrong\n";
      }

    }
  /*if(i<=algominis.size()){
    Minimizer algominino=newminisalgono[i];
    if(!(bfmini.getSequence()==algominino.getSequence() && bfmini.getPosition()==algominino.getPosition())){
      rightMinisno=false;
      //cout<<"here wrong\n";
      }
    }*/
  }
  cout<<"Time needed for algo without dynseq: "<< msalgono<<"miliseconds\n";
  cout<<"Time needed for algo: "<< msalgo<<"miliseconds\n";
  cout<<"Time needed for bf: "<< msbf<<"miliseconds\n";
  //cout<<"Algo: "<<algo_result<<"\n";
  //cout<<"Algo no dynseq: "<<sequence<<"\n";
  if(bf_result.compare(algo_result)==0){

    cout<<"The algorithm returned the right sequence!\n";
  }
  if(bf_result.compare(sequence)==0){

    cout<<"The algorithm no dynseq returned the right sequence!\n";
  }
  if(rightMinis){
    cout<<"The algorithm delivered the right minimizers!\n";
  }
  /*if(rightMinisno){
    cout<<"The algorithm no dynseq delivered the right minimizers!\n";
  }*/
  else{
    cout<<"ERROR\n";
  }
  brute_force_minimizer_computation_normal_string(minimizerTreeBF,sequence2,variants3,k,w);

  //std::vector<Minimizer> newminimethod=minimizer_to_vector(minimizerTree);
  /*int pos=3;
  int shift=2;
  minimizerTree->shift_greater(pos,shift);
  print_minimizerTree(minimizerTree);
  cout<<"End of Tree\n";
  pos=9;
  minimizerTree->shift_greater(pos,shift);
  print_minimizerTree(minimizerTree);
  cout<<"End of Tree\n";
  int left=5;
  int right=23;
  delete_minimizers(minimizerTree,left,right);
  print_minimizerTree(minimizerTree);*/

}
