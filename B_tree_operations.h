////////////////////////////////////////////////////////////////////////////////
// B_tree_operations.h
//   B-tree header file.
//
// additional functions for the access and update of the B-tree
//
////////////////////////////////////////////////////////////////////////////////
// author: Alexander Petri

#ifndef B_TREE_OPERATIONS_H
#define B_TREE_OPERATIONS_H


#include "B-tree.hh"
#include "B_tree_node.hh"

using namespace std;
using namespace md;

/*!
 * Print all elements stored in the B-tree to the command line
 * @param minimizerTree:    the B-tree to be printed
 */
void print_minimizerTree(B_tree<int,std::string,7,3>* minimizerTree){
  int i = 0;
  for(auto elem: *minimizerTree){
    int minikey=elem.first;
    auto elem3 = minimizerTree->search(minikey);
    std::vector<std::string> es2=*elem3.key->satellites;
    std::string sequence=es2[0];
    cout<<"Minimizer at "<<elem.first<<":" <<sequence<<"\n";
    i+= 1;
  }
}
/*!
 * Fill the B-tree with the previously computed minimizers
 * @param minimizerTree:    the B-tree to be filled
 * @param minis:    the set of minimizers to be stored in the B-tree
 */
void fill_minimizer_tree(B_tree<int,std::string,7,3>* minimizerTree,vector<Minimizer> minis){
  for(int i=0;i<minis.size();i++){
    int position=minis[i].getPosition();
    std::string satelliteval =minis[i].getSequence();
    minimizerTree->insert(position,satelliteval);
  }
}

/*!
 * Delete the elements having a key with left<=key<=right.
 * @param minimizerTree:    the B-tree to be altered
 * @param left:  the lower bound of the range
 * @param right:  the upper bound of the range
 */
void delete_minimizers_inefficient(B_tree<int,std::string,7,3>* minimizerTree,int& left,int& right){
  for(int i = left; i <=right; i+=1){
    cout<<"Removing "<<i<<" \n";
    minimizerTree->remove(i);
  }
}

/*!
 * Updating the B-tree by deleting old minimizers and filling the tree with the updated minimizers.
 * @param minimizerTree:    the B-tree to be updated
 * @param fullsubseq:  the substring which the new minimzers are generated for
 * @param thisstartpos:  the index of the substrings starting position in the whole DNA sequences
 * @param k_size: length of the k-mers
 * @param w_size: size of the window
 * @param var_impact_shift: the length by which subsequent minimizers key have to be shifted
 */
void update_minimizerTree(B_tree<int,std::string,7,3>* minimizerTree,std::string& fullsubseq,int& thisstartpos,int& k_size,int& w_size, int& var_impact_shift){
//generate the minimizers for the updated subsequence
  std::vector<Minimizer> newminis=get_kmer_minimizers_algo(fullsubseq, k_size, w_size,thisstartpos);
  //find the positions of the first and last minimizer in the new set
  int start=newminis.front().getPosition();
  cout<<start<<" is the first new minimizer\n";
  int end=newminis.back().getPosition();//end after shift
  //update last minimizers' position to find the right end position for the deletion
  int newend = end-var_impact_shift;
  //find the last minimizer position in the minimizer tree
  int lastminipos=minimizerTree->get_max();

cout<<"Shifting the elements by "<<var_impact_shift<<"\n";
cout<<"Newend "<<newend<<", Lastminipos: "<<lastminipos<<"\n";

//as the B-tree does not deliver stop if after last element set position of last deleted minimizer to maximum in tree
  int suc=0;
  /*if(end>lastminipos){
    end=lastminipos;
    suc=-1;
  }*/
  if(newend>=lastminipos){
    newend=lastminipos;
    suc=0;
  }
  else{
    cout<<"I want to find "<<newend<<" in the tree \n";
    //int suc=minimizerTree->search(newend).key->value;
    //auto succ=minimizerTree->successor(newend).key->value;
    int suc=minimizerTree->successor(newend).key->value;
    //auto suc=minimizerTree->predecessor(succ).key->value;
    cout<<"found suc "<<suc<<"\n";
  }
  cout<<"Hello World\n";
  cout<<"deleting the minimizers in the interval ("<<start<<", "<<newend<<")\n";
  //cout<<"deleting the minimizers in the interval ("<<start<<", "<<newend<<")\n";
  cout<<"last element in minimizerTree "<<lastminipos<<"\n";
  //find the boundaries for the minimizers which are to be deleted
  if(thisstartpos==0){
    start=minimizerTree->get_min();
  }

  //int right=minimizerTree->search(end).key->value;
  //delete all minimizers between left and right
if(!(start>minimizerTree->get_max())){
  delete_minimizers_inefficient(minimizerTree,start,newend);
  //delete_minimizers_iterator(minimizerTree,start,newend);
}
  if(suc>0){
    minimizerTree->shift_greater(suc,var_impact_shift);
  }
  print_minimizerTree(minimizerTree);
  cout<<"New Minimizers to be added:\n";
  for(int i=0;i<newminis.size();i++){
    newminis[i].printMinimizer();
  }
  cout<<"printing new minimizers done\n";
  fill_minimizer_tree(minimizerTree,newminis);
  print_minimizerTree(minimizerTree);
  cout<<"Updating done\n";
}

/*!
 * Copy the elements stored in the B-tree into a vector
 * @param minimizerTree:    the B-tree to be copied
 */
std::vector<Minimizer> minimizer_to_vector(B_tree<int,std::string,7,3>* minimizerTree){
  std::vector<Minimizer> minimizers;
  for(auto elem: *minimizerTree){
    int minikey=elem.first;
    auto elem3 = minimizerTree->search(minikey);
    std::vector<std::string> es2=*elem3.key->satellites;
    std::string sequence=es2[0];
    Minimizer mini=Minimizer(minikey,sequence);
    minimizers.push_back(mini);
  }
  return minimizers;
}

#endif
