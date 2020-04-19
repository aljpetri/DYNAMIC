////////////////////////////////////////////////////////////////////////////////
// Minimizer.h
//   minimizer class header file.
//
//  class to internally represent a minimizer
//
////////////////////////////////////////////////////////////////////////////////
//  author: Alexander Petri

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "main.h"

/*
* Class to internally represent minimizers
*
* @param position    the position of the minimizer
* @param sequence    the sequence of the minimizer
*
*/

class Minimizer{
private:
  int position;
  string sequence;
public:
  //Custom constructor
  Minimizer(int& pos, string& seq){
    position=pos;
    sequence=seq;
  }
  //Default constructor
  Minimizer() = default;
  void alterposition(int& pos){
    position=pos;
  }

  /*
  * prints a minimizer to the console
  * Output: position: sequence
  */
  void printMinimizer(){
    cout<<position<<": "<<sequence<<"\n";
  }


  /*
   *updates the minimizer with the new parameters
   *@param pos: the new position
   *@param seq: the new sequence
   */
  void updateMinimizer(int& pos, string& seq){
    position=pos;
    sequence=seq;
  }
  /*
   *returns the position of the minimizer
   */
  int getPosition(){
    return position;
  }
  /*
   *returns the sequence of the minimizer
   */
  std::string getSequence(){
    return sequence;
  }
};

#endif
