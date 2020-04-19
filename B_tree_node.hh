////////////////////////////////////////////////////////////////////////////////
// B_tree_node.hh
//   B-tree Node header file.
////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019 Simon J. Puglisi and Massimiliano Rossi
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#ifndef MD_B_TREE_NODE__HH
#define MD_B_TREE_NODE__HH


#include <vector>
#include <assert.h>
#include <typeinfo>
#include <type_traits>
#include <limits>

namespace md{

  template< typename T>
  constexpr bool in_range_unsigned(const long long int x)
  {
    return  ((x-std::numeric_limits<T>::min()) <= (std::numeric_limits<T>::max()-std::numeric_limits<T>::min()));
  }

  /*!
   * K is the type of the key values
   * S is the type of the satellites
   * B is the number of the pivots
   * T is the minimum degree of the B-tree
   */
  template< typename K, typename S, size_t B = 63, size_t T = 3>
  class B_tree_node{
  public:
    typedef struct t_key{
      K value;
      std::vector<S>* satellites;

      ~t_key(){
        // if(satellites != nullptr)
        //     delete satellites;
      }

      t_key(){
        value = 0;
        satellites = nullptr;
      }

    } key_t;

    typedef struct t_shifted_key_ptr{
      K shift;
      key_t* key;

      ~t_shifted_key_ptr(){
        // NtD
      }

      t_shifted_key_ptr():
        shift(0),
        key(nullptr)
      {
        // NtD
      }

      t_shifted_key_ptr(key_t* key_, K shift_):
          shift(shift_),
          key(key_)
      {
        // NtD
      }

      void do_shift(K shift_){
        shift += shift_;
      }

      key_t shift_key(){
        key_t tmp = *key;
        tmp.value += shift;
        shift = 0;
        return tmp;
      }

    } shifted_key_ptr_t;



    typedef typename std::conditional<in_range_unsigned<uint8_t>(B),uint8_t,
                              typename std::conditional<in_range_unsigned<uint16_t>(B), uint16_t ,
                                typename std::conditional<in_range_unsigned<uint32_t>(B), uint32_t ,
                                  uint64_t
                                >::type
                              >::type
                            >::type B_t;


    /*!
     * Costructor
     */
    B_tree_node(bool _is_leaf = true);

    /*!
     * Desctructor
     */
    ~B_tree_node();

    /*!
     * Tells if a node is a leaf
     * @return true if the node is a leaf, false otherwise
     */
    bool is_leaf();

    /*!
     * Tells if a node is a full
     * @return true if the node is a full, false otherwise
     */
    bool is_full();

    /*!
     * Finds the element of key value in the subtree rooted in this node.
     * @param  value the key of the element we look for.
     * @return       the element in the tree.
     */
    shifted_key_ptr_t search(K value);
    /*!
     * Finds the element of key value in the subtree rooted in this node and shifts all values having an
     * equal or greater key than key value.
     * ATTENTION: Does not support shifting with a key, which is not in the B-tree
     * @param  value the key of the element we look for.
     * @param  shift: the shift which is applied to the greater elements
     * @return       a pointer to the element in the tree.
     */
    shifted_key_ptr_t shift_greater(K value, K shift);

    /*!
     * Insert the element value with its satellite in the tree
     * @param value     the key of the element that has to be inserted.
     * @param satellite the satellite information attached to the element.
     */
    shifted_key_ptr_t insert(K value, S satellite);

    /*!
     * Remove the element with key value from the tree
     * @param  value the key value of the element to be removed
     * @return       the removed element with its satellite informations.
     */
    key_t remove(K value);

    /*!
     * Split the child at index i into two.
     * @param index the index of the child to be splitted
     */
    void split_child(B_t index);

    /*!
     * Access the i-th child of the node.
     * @param i the index of the child to be returned
     * @return  the i-th child of the node if it exists.
     */
    B_tree_node<K,S,B,T>* get_child(B_t i);

    /*!
     * Access the i-th key of the node.
     * @param i the index of the key to be returned
     * @return  the i-th key of the node if it exists.
     */
    shifted_key_ptr_t get_key(B_t i);

    /*!
     * The number of keys in the node.
     * @return the number of keys in the node.
     */
    B_t get_n_keys();

    /*!
     * Set the child as the i-th child of the node.
     * @param i the index of the child to be inserted into
     */
    void set_child(B_t i, B_tree_node<K,S,B,T>* child);

    /*!
     * Finds the predecessor of the element of key value in the subtree rooted in this node.
     * @param  value the key of the element we look for.
     * @return       a pointer to the element in the tree.
     */
    shifted_key_ptr_t predecessor(K value);

    /*!
     * Finds the successor of the element of key value in the subtree rooted in this node.
     * @param  value the key of the element we look for.
     * @return       a pointer to the element in the tree.
     */
    shifted_key_ptr_t successor(K value);

    friend void swap(B_tree_node<K,S,B,T>& first, B_tree_node<K,S,B,T>& second)
    {
        using std::swap;

        for(B_t i = 0; i < B; ++i){
          swap(first.keys[i], second.keys[i]);
        }
        swap(first.n, second.n);
        swap(first.children, second.children);
    }

    /*!
     * Join two B-trees. It requirest that all the elements of other are greater
     * than the greatest element of this.
     * @param other The other B-tree.
     */
    void join(B_tree_node<K,S,B,T>* other);

    B_tree_node<K,S,B,T>* split(const K &value_, size_t *h_this = nullptr, size_t *h_rhs = nullptr);

    /*!
     * The height h of the tree rooted in this node.
     * Complexity: O(h)
     * @return the height of the tree rooted in this node
     */
    K height();

    /*!
     * The max element of the tree rooted in this node.
     * Complexity: O(log_B(n))
     * @return the value of the max element of the tree rooted in this node
     */
    K get_max();

    /*!
     * The min element of the tree rooted in this node.
     * Complexity: O(log_B(n))
     * @return the value of the min element of the tree rooted in this node
     */
    K get_min();

    /*!
     * Shift the keys in the subtree rooted in this node by shift_
     * @param shift_ the value of the shift
     */
    void shift(K shift_);


    /*!
     * The value of the shift.
     * @return the value of the shift
     */
    K get_shift();

    bool check_integrity(){
      if(is_leaf()){
        for(B_t i = 0; i < n; ++i){
          assert(keys[i].satellites != nullptr);
        }
        for(B_t i = n; i < B; ++i){
          assert(keys[i].satellites == nullptr);
        }
      }else{
        for(B_t i = 0; i < n; ++i){
          assert(keys[i].satellites != nullptr);
          assert(children[i] != nullptr);
        }
        assert(children[n] != nullptr);
        if(n < B) assert(keys[n].satellites == nullptr);
        for(B_t i = n+1; i < B; ++i){
          assert(keys[i].satellites == nullptr);
          assert(children[i] == nullptr);
        }
        for(B_t i = 0; i < n+1; ++i){
          assert(children[i]->check_integrity());
        }
      }
      return true;
    }

  protected:

    void merge_children(B_t i, size_t* h_this = nullptr);

    void fuse_children(B_t i, size_t* h_this = nullptr);
    /*!
     * Shift keys and children pointers to the left from position pos on
     * i.e., pos <- pos + 1
     * @param pos the position that has to be overwritten.
     */
    void shift_left(B_t pos);

    /*!
     * Shift keys and children pointers to the right from position pos on, by offset elements
     * i.e., pos + offset <- pos
     * @param pos the position that has to be overwritten.
     */
    void shift_right(B_t pos, B_t offset = 1);

    /*!
     * Join the lhs subtree on the left spine of this tree, using the max_value
     * as pivot. The max_value must not be in the lhs subtree.
     * @param lhs     the subtree to be joint
     * @param max_value the pivot element
     * @param h_this    the height of the tree
     * @param h_lhs     the height of the subtree to be joint
     */
    void join_left(B_tree_node<K,S,B,T>* lhs, key_t max_value , size_t *h_this = nullptr, size_t *h_lhs = nullptr);

    /*!
     * Join the lhs subtree on the right spine of this tree, using the min_value
     * as pivot. The min_value must not be in the rhs subtree.
     * @param rhs     the subtree to be joint
     * @param min_value the pivot element
     * @param h_this    the height of the tree
     * @param h_rhs     the height of the subtree to be joint
     */
    void join_right(B_tree_node<K,S,B,T>* rhs, key_t min_value , size_t *h_this = nullptr, size_t *h_rhs = nullptr);


  private:
    key_t keys[B];            // The keys of the node.
    B_tree_node<K,S,B,T>** children; // The pointers to the children of the node. (if nullptr the node is a leaf)
    B_t n;                    // The number of keys in the node.
    K _shift;                    // The shift value of the node.

  }; // B_tree_node

  // Ctor
  template< typename K, typename S, size_t B, size_t T>
  B_tree_node<K,S,B,T>::B_tree_node(bool _is_leaf):
    children(nullptr),
    n(0),
    _shift(0)
  {
    if(!_is_leaf){
      children = new B_tree_node<K,S,B,T>*[B+1];
      for(B_t i = 0; i < B+1; ++i)
        children[i] = nullptr;

    }

  }

  // Dtor
  template< typename K, typename S, size_t B, size_t T>
  B_tree_node<K,S,B,T>::~B_tree_node()
  {
    if(children != nullptr){
      for(B_t i = 0; i < n+1; ++i)
        if(children[i] != nullptr)
          delete children[i];

      delete[] children;
    }

    for(B_t i = 0; i < B; ++i){
      if(keys[i].satellites != nullptr){
        assert(i<n);
        delete keys[i].satellites;
      }
    }
  }


  /*!
   * Access the i-th child of the node.
   * @param i the index of the child to be returned
   * @return  the i-th child of the node if it exists.
   */
  template< typename K, typename S, size_t B, size_t T>
  B_tree_node<K,S,B,T>* B_tree_node<K,S,B,T>::get_child(B_t i)
  {
      if(i < n+1 && !is_leaf()){
          return children[i];
      }
      return nullptr;
  }

  /*!
   * Access the i-th key of the node.
   * @param i the index of the key to be returned
   * @return  the i-th key of the node if it exists.
   */

  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::get_key(B_t i)
  {
    if(i < n){
      return shifted_key_ptr_t(&keys[i],_shift);
    }
    return shifted_key_ptr_t( nullptr, _shift );
  }

  /*!
   * The number of keys in the node.
   * @return the number of keys in the node.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::B_t B_tree_node<K,S,B,T>::get_n_keys(){
    return n;
  }

  /*!
   * Set the child as the i-th child of the node.
   * @param i the index of the child to be inserted into
   */
  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::set_child(B_t i, B_tree_node<K,S,B,T>* child)
  {
    if(children == nullptr)
      children = new B_tree_node<K,S,B,T>*[B+1];

    children[i] = child;
  }

  template< typename K, typename S, size_t B, size_t T>
  bool B_tree_node<K,S,B,T>::is_leaf()
  {

    return (children == nullptr);

  }

  template< typename K, typename S, size_t B, size_t T>
  bool B_tree_node<K,S,B,T>::is_full()
  {

    return (n == B);

  }

  template< typename K, typename S, size_t B, size_t T>
  K B_tree_node<K,S,B,T>::height()
  {
    // if(n == 0 && children == nullptr) return 0;
    // else
    if(children == nullptr) return 1;
    else return (children[0]->height() + 1);
  }

  template< typename K, typename S, size_t B, size_t T>
  K B_tree_node<K,S,B,T>::get_max()
  {
    if(children == nullptr) return (keys[n-1].value + _shift);
    else return (children[n]->get_max() + _shift);
  }

  template< typename K, typename S, size_t B, size_t T>
  K B_tree_node<K,S,B,T>::get_min()
  {
    if(children == nullptr) return (keys[0].value  + _shift);
    else return (children[0]->get_min() + _shift);
  }

  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::shift(K shift_)
  {
    _shift += shift_;
  }

  template< typename K, typename S, size_t B, size_t T>
  K B_tree_node<K,S,B,T>::get_shift()
  {
    return _shift;
  }

    /*!
     * Finds the element of key value in the subtree rooted in this node and shifts all values having an
     * equal or greater key than key value.
     * ATTENTION: Does not support shifting with a key, which is not in the B-tree
     * @param  value the key of the element we look for.
     * @param  shift: the shift which is applied to the greater elements
     * @return       a pointer to the element in the tree.
     */
    template< typename K, typename S, size_t B, size_t T>
    typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::shift_greater(K value,K shift)
    {
      cout<<"internal shift: "<<_shift<<"\n";
      value -= _shift;
      cout<<"val"<<value<<"\n";
      size_t l = 0;
      size_t r = n;

      // Find the predecessor and successor in the node by performing a binary search
      while(r > l + 1){
        size_t mid = (l+r)/2;
        if(keys[mid].value <= value){
          l = mid;
        }else{
          r = mid;
        }
      }
      cout<<"left "<<l<<" right "<<r<<"\n";
      //if we have found the key in this node, shift all children and keys which are on the right of the key
      if(keys[l].value == value && keys[l].satellites != nullptr){
        for(B_t s=l;s<=n;s++){
          keys[s].value=keys[s].value+shift;
          auto keyt =get_key(s);
          cout<<"Before error\n";
          if(s>l){
            if(!is_leaf() && children[s] != nullptr){
              auto child=children[s];
              cout<<"shift: "<<child->get_shift()<<"\n";
              child->shift(shift);
              cout<<"shift: "<<child->get_shift()<<"\n";
            }
          }
          cout<<"Shifting completed\n";
        }
        return shifted_key_ptr_t(&keys[l], _shift);
      }

      // If it is greater than the largest element among the pivots the child is the rightmost one
      if(keys[l].value < value){
        l++;
      }
      for(B_t s=l;s<=n;s++){
        keys[s].value=keys[s].value+shift;
        if(s>l){
          if(!(children == nullptr)){
            auto child=children[s];
            cout<<"shift: "<<child->get_shift()<<"\n";
            child->shift(shift);
            cout<<"shift: "<<child->get_shift()<<"\n";
          }
        }
        cout<<"Shift done\n";
      }
      if(!is_leaf() && children[l] != nullptr){
        shifted_key_ptr_t tmp_ptr = children[l]->shift_greater(value,shift);
        tmp_ptr.do_shift(_shift);
        return tmp_ptr;
      }
      return shifted_key_ptr_t( nullptr, _shift);
    }

  /*!
   * Finds the element of key value in the subtree rooted in this node.
   * @param  value the key of the element we look for.
   * @return       a pointer to the element in the tree.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::search(K value)
  {
    // shift the value
    value -= _shift;

    size_t l = 0;
    size_t r = n;

    // Find the predecessor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    if(keys[l].value == value && keys[l].satellites != nullptr) return shifted_key_ptr_t(&keys[l], _shift);

    // If it is greater than the largest element among the pivots the children is the rightmost one
    if(keys[l].value < value) l++;

    if(!is_leaf() && children[l] != nullptr){

      shifted_key_ptr_t tmp_ptr = children[l]->search(value);
      tmp_ptr.do_shift(_shift);

      return tmp_ptr;
    }

    return shifted_key_ptr_t( nullptr, _shift);
  }

  /*!
   * Finds the predecessor of the element of key value in the subtree rooted in this node.
   * @param  value the key of the element we look for.
   * @return       a pointer to the element in the tree.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::predecessor(K value)
  {
    // shift the value
    value -= _shift;

    size_t l = 0;
    size_t r = n;

    // Find the predecessor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    if(is_leaf() && keys[l].value <= value && keys[l].satellites != nullptr) return shifted_key_ptr_t(&keys[l],_shift);

    // If it is greater than the largest element among the pivots the children is the rightmost one
    if(keys[l].value < value) l++;

    shifted_key_ptr_t ans( nullptr, _shift);

    if(!is_leaf() && children[l] != nullptr){
      ans = children[l]->predecessor( value );
      ans.do_shift(_shift);
    }

    if(ans.key == nullptr && l > 0 && keys[l-1].value <= value){
      ans.key = &keys[l-1];
    }

    return ans;
  }

  /*!
   * Finds the successor of the element of key value in the subtree rooted in this node.
   * @param  value the key of the element we look for.
   * @return       a pointer to the element in the tree.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::successor(K value)
  {

    // shift the value
    value -= _shift;

    size_t l = 0;
    size_t r = n;

    // Find the successor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    // If it is greater than the largest element among the pivots the children is the rightmost one
    if(keys[l].value > value) r--;

    if(is_leaf() && r < n && keys[r].value > value && keys[r].satellites != nullptr) return shifted_key_ptr_t(&keys[r], _shift);

    shifted_key_ptr_t ans( nullptr, _shift);

    if(!is_leaf() && children[r] != nullptr){
      ans = children[r]->successor(value);
      ans.do_shift(_shift);
    }

    if(ans.key == nullptr && r < n && keys[r].value > value){
      ans.key = &keys[r];
    }

    return ans;
  }

  /*!
   * Insert the element value with its satellite in the tree
   * @param value     the key of the element that has to be inserted.
   * @param satellite the satellite information attached to the element.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::shifted_key_ptr_t B_tree_node<K,S,B,T>::insert(K value, S satellite)
  {
    // shift the value
    value -= _shift;

    // If I am a leaf, insert the element in the correct position
    if(is_leaf()){
      // If the element is already in the tree, append the satellites to the element.
      shifted_key_ptr_t elem = search(value + _shift);

      if(elem.key != nullptr){
        elem.key->satellites->push_back(satellite);
        return elem;
      }

      // Create the new element
      keys[n].value = value;
      keys[n].satellites = new std::vector<S>();
      keys[n].satellites->push_back(satellite);

      // Bubble the new element in the correct position
      B_t i = 0;

      while(i < n && keys[n-i-1].value > value){
        std::swap(keys[n-i-1],keys[n-i]);
        i++;
      }

      n++;
      return shifted_key_ptr_t( &keys[n-i-1], _shift );
    }

    size_t l = 0;
    size_t r = n;

    // Find the predecessor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    // If we found the element we append the satellite to its list
    if(keys[l].value == value){
      keys[l].satellites->push_back(satellite);
      return shifted_key_ptr_t( &keys[l], _shift);
    }

    if(keys[l].value < value) l++;

    if(children[l] != nullptr){
        // if the child is full, split it
        if(children[l]->is_full()){
          split_child(l);
          // If we found the element we append the satellite to its list
          if(keys[l].value == value){
            keys[l].satellites->push_back(satellite);
            return shifted_key_ptr_t( &keys[l], _shift);
          }

          if(keys[l].value < value) l++;
        }
    }else{
      children[l] = new B_tree_node<K,S,B,T>();
    }

    // Insert the element in the child
    shifted_key_ptr_t ans = children[l]->insert(value,satellite);

    ans.do_shift( _shift );
    return ans;
  }

  /*!
   * Remove the element with key value from the tree
   * @param  value the key value of the element to be removed
   * @return       the removed element with its satellite informations.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree_node<K,S,B,T>::key_t B_tree_node<K,S,B,T>::remove(K value)
  {

    // shift the value
    value -= _shift;

    // Find value in the node
    size_t l = 0;
    size_t r = n;

    // Find the predecessor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    key_t res;
    if(is_leaf()){
      // Case 1. If value is in the node and the node is a leaf

      if(keys[l].value == value){
        // If found remove the element from the node
        n--;
        res = keys[l];
        while(l < n){
          keys[l] = keys[l+1];
          l++;
        }
        keys[n].value = 0;
        keys[n].satellites = nullptr;

        // Remove the memory of the satellites
        assert(res.satellites != nullptr);

      }
    }else if(keys[l].value == value){
      // Case 2. If the value is in the node and the node is an internal node.
      B_tree_node<K,S,B,T>* y = children[l];
      B_tree_node<K,S,B,T>* z = children[l+1];
      if(y->n >= T){
        // Case 2.a) If the child y that precedes value in the node has at least T keys, then
        // find the predecessor of value in the subtree rooted in y.
        shifted_key_ptr_t pred = y->predecessor(value);

        // Recursively delete pred and replace value by pred.
        res = keys[l];

        keys[l] = pred.shift_key();

        y->remove(keys[l].value);

      }else if(z->n >= T){
        // Case 2.b) Simmetrically examine the child z that follows value in the node.
        // If k has at least T keys, then find the successor of value in the subtree rooted in z.
        shifted_key_ptr_t succ = z->successor(value);

        // Recursively delete succ and replace value by pred.
        res = keys[l];

        keys[l] = succ.shift_key();

        z->remove(keys[l].value);
      }else{
        // Case 2.c) Otherwise, if both y and z have less than T - 1 keys.
        // merge y and z.
        merge_children(l);
        res = remove(value);
      }

    }else{
      // Case 3. If value is not in the node. Determine the root of the subtree that must contain value.
      if(keys[l].value < value) l++;
      B_tree_node<K,S,B,T>* c = children[l];
      B_tree_node<K,S,B,T>* lhs = nullptr;
      B_tree_node<K,S,B,T>* rhs = nullptr;
      if(l > 0) lhs = children[l-1];
      if(l < n) rhs = children[l+1];
      // If c has only T - 1 keys
      if(c->n < T){
        // Case 3.a) If c has an immediate sibling with at least T keys.
        if(lhs != nullptr && lhs->n >= T){
          // Move an extra key from this node down into c.

          // // Bubble the element to the right

          c->shift_right(0,1);

          c->keys[0] = keys[l-1];
          c->keys[0].value -= c->_shift; // Shift the key

          // Move a key from the c's sibling
          keys[l-1] = lhs->keys[lhs->n-1];
          keys[l-1].value += lhs->_shift; // Shift the key

          lhs->keys[lhs->n-1].value = 0;
          lhs->keys[lhs->n-1].satellites = nullptr;

          if(!c->is_leaf()){
            c->children[0] = lhs->children[lhs->n];
            // Shift the child
            c->children[0]->_shift+= lhs->_shift;
            c->children[0]->_shift-= c->_shift;
          }

          if(!lhs->is_leaf()) lhs->children[lhs->n] = nullptr;
          lhs->n--;


          res = c->remove(value);
        }else if(rhs != nullptr && rhs->n >= T){
          // Move an extra key from this node down into c.

          c->keys[c->n] = keys[l];
          c->keys[c->n].value -= c->_shift; // Shift the key
          c->n++;
          // Move a key from c's sibling
          keys[l] = rhs->keys[0];
          keys[l].value += rhs->_shift; // Shift the key

          rhs->keys[0].value = 0;
          rhs->keys[0].satellites = nullptr;

          if(!c->is_leaf()){
            c->children[c->n] = rhs->children[0];
            // Shift the child
            c->children[c->n]->_shift+= rhs->_shift;
            c->children[c->n]->_shift-= c->_shift;
          }
          // Bubble the element of rhs to the left

          rhs->shift_left(0);

          res = c->remove(value);
        }else{
          // Case 3.b) If c and both its immediate siblings have T-1 keys.
          if(l > 0){
            merge_children(l-1);
            // res = children[l-1]->remove(value);
          }else{
            merge_children(l);
            // res = children[l]->remove(value);
          }
          value += _shift;
          res = remove(value);
          res.value -= _shift;
        }
      }else{
        res = c->remove(value);
      }


    }
    // shift the result
    res.value += _shift;

    return res;
  }

  /*!
   * Split the child at index i into two.
   * @param index the index of the child to be splitted
   */
  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::split_child(B_t index)
  {

    B_tree_node<K,S,B,T>* lhs = children[index];
    B_tree_node<K,S,B,T>* rhs = new B_tree_node<K,S,B,T>(lhs->is_leaf());

    // Shift rhs
    rhs->_shift = lhs->_shift;

    B_t median_pos = lhs->n/2;

    if(lhs->is_leaf()){

     //  Copy the elements from median + 1 on into rhs
     for(B_t i = median_pos + 1 ; i < lhs->n; ++i){
       std::swap(lhs->keys[i], rhs->keys[i-median_pos-1]);
     }

   }else{

     // Copy the elements from median + 1 on into rhs
     for(B_t i = median_pos + 1 ; i < lhs->n; ++i){
       std::swap(lhs->children[i], rhs->children[i-median_pos-1]);
       std::swap(lhs->keys[i], rhs->keys[i-median_pos-1]);
     }
     std::swap(lhs->children[lhs->n], rhs->children[lhs->n-median_pos-1]);

   }

   // Copy the median into the last key of the root node
   std::swap(lhs->keys[median_pos], keys[n]);
   // shift the median element
   keys[n].value += lhs->_shift;
   children[n+1] = rhs;
   // Bubble the median in the correct place
   B_t i = n;
   while(i > 0 && keys[i-1].value > keys[i].value){

     std::swap(children[i], children[i+1]);
     std::swap(keys[i-1], keys[i]);

     i--;
   }
   n++;
   rhs->n = lhs->n - (median_pos + 1);
   lhs->n = median_pos;

  }

  /*!
   * Shift keys and children pointers to the left from position pos on
   * i.e., pos <- pos + 1
   * @param pos the position that has to be overwritten.
   */
   template< typename K, typename S, size_t B, size_t T>
   void B_tree_node<K,S,B,T>::shift_left(B_t pos)
  {
    // Shift left the elements in the node
    if(is_leaf()){

      for(B_t j = pos; j < n-1;++j ){
        keys[j] = keys[j+1];
      }

    }else{

      for(B_t j = pos; j < n-1;++j ){
        keys[j] = keys[j+1];
        children[j] = children[j+1];
      }
      children[n-1] = children[n];

      children[n] = nullptr;
    }
    keys[n-1].value = 0;
    keys[n-1].satellites = nullptr;
    n--;
  }

  /*!
   * Shift keys and children pointers to the right from position pos on, by offset elements
   * i.e., pos + offset <- pos
   * @param pos the position that has to be overwritten.
   */
   template< typename K, typename S, size_t B, size_t T>
   void B_tree_node<K,S,B,T>::shift_right(B_t pos, B_t offset)
  {
    // Shift left the elements in the node
    if(is_leaf()){

      for(B_t j = 0; j < n-pos;++j ){
        keys[n -1 - j + offset] = keys[n -1 - j];
      }

      for(B_t j = pos; j < pos + offset;++j ){
        keys[j].value = 0;
        keys[j].satellites = nullptr;
      }

    }else{

      children[n + offset] = children[n];
      for(B_t j = 0; j < n-pos; ++j ){
        keys[n -1 - j + offset] = keys[n -1 - j];
        children[n -1 - j + offset] = children[n -1 - j];
      }

      for(B_t j = pos; j < pos + offset;++j ){
        keys[j].value = 0;
        keys[j].satellites = nullptr;
        children[j] = nullptr;
      }

    }
    n = n + offset;
  }

  /*!
   * Merge the children i and i+1 together
   * @param i the position of the key such that the preceding and following children have to be merged
   */
  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::merge_children(B_t i, size_t* h_this)
  {
    assert(i < n);
    assert( (B_t)(children[i]->n + children[i+1]->n + 1) <= B);

    B_tree_node<K,S,B,T>* lhs = children[i];
    B_tree_node<K,S,B,T>* rhs = children[i+1];

    B_t median_pos = lhs->n;
    // Move the i-th key as median of the new node.
    std::swap(lhs->keys[lhs->n], keys[i]);
    // Shift the key in lhs
    lhs->keys[lhs->n].value -= lhs->_shift;

    // Copy the keys of rhs into lhs
    if(lhs->is_leaf()){

     //  Copy the elements from rhs into lhs
     for(B_t i = 0 ; i < rhs->n; ++i){
       std::swap(lhs->keys[median_pos + i + 1], rhs->keys[i]);
       // Shift the key in lhs
       lhs->keys[median_pos + i + 1].value += rhs->_shift;
       lhs->keys[median_pos + i + 1].value -= lhs->_shift;
     }

   }else{

     // Copy the elements from rhs into lhs
     for(B_t i = 0 ; i < rhs->n; ++i){
       std::swap(lhs->children[median_pos + i + 1], rhs->children[i]);
       std::swap(lhs->keys[median_pos + i + 1], rhs->keys[i]);
       // Shift the key in lhs
       lhs->keys[median_pos + i + 1].value += rhs->_shift;
       lhs->keys[median_pos + i + 1].value -= lhs->_shift;
       // Shift the children in lhs
       lhs->children[median_pos + i + 1]->_shift += rhs->_shift;
       lhs->children[median_pos + i + 1]->_shift -= lhs->_shift;
     }
     std::swap(lhs->children[median_pos + rhs->n + 1], rhs->children[rhs->n]);
     // Shift the child in lhs
     lhs->children[median_pos + rhs->n + 1]->_shift += rhs->_shift;
     lhs->children[median_pos + rhs->n + 1]->_shift -= lhs->_shift;

   }
   lhs->n = median_pos + rhs->n + 1;

   // Shift left the elements in the node
   std::swap(children[i],children[i+1]);

   shift_left(i);
   delete rhs;

   // Shrinking the height
   if(n == 0){
     // std::swap(*this,*lhs);
     for(B_t i = 0; i < B; ++i){
       std::swap(keys[i], lhs->keys[i]);
     }
     std::swap(children, lhs->children);
     std::swap(n, lhs->n);
     _shift+= lhs->_shift;

     lhs->children[0] = nullptr;
     delete lhs;

     // Update the height of the subtree
     if(h_this != nullptr) (*h_this)--;
   }

  }

  /*!
   * Fuse the children i and i+1 together
   * @param i the position of the key that has to be fused with following children have to be merged
   */
  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::fuse_children(B_t i, size_t* h_this)
  {
    assert(i < n);

    B_tree_node<K,S,B,T>* lhs = children[i];
    B_tree_node<K,S,B,T>* rhs = children[i+1];

    // If the total number of keys is smaller than or equals to B then merge.
    if( (B_t)(lhs->n + rhs->n + 1) <= B){
      merge_children(i, h_this);
      return;
    }else{
      // I have to balance the keys in the two children

      // Find the median position across the two children
      size_t median_pos = (lhs->n + rhs->n + 1) >> 1;

      // The children are balanced
      if(median_pos == lhs->n || median_pos == rhs->n)
        return;

      if(median_pos < lhs->n){
        // The median key is in the left child

        // I have to make room on the right child for lhs->n - median_pos + 1 elements
        B_t offset = lhs->n - median_pos;

        if(rhs->is_leaf()){

          for(B_t i = 0; i < rhs->n; ++i){
            rhs->keys[rhs->n -1 - i + offset] = rhs->keys[rhs->n -1 - i];
          }

        }else{

          rhs->children[rhs->n + offset] = rhs->children[rhs->n];
          for(B_t i = 0; i < rhs->n; ++i){
            rhs->keys[rhs->n -1 - i + offset] = rhs->keys[rhs->n -1 - i];
            rhs->children[rhs->n -1 - i + offset] = rhs->children[rhs->n -1 - i];
          }

        }

        // Copy the i-th key in the correct position
        rhs->keys[offset - 1] = keys[i];
        // Shift the key in rhs
        rhs->keys[offset - 1].value -= rhs->_shift;

        // Copy the last offset -1 keys from lhs into rhs
        if(rhs->is_leaf()){

          for(B_t i = 0; i < offset-1; ++i){
            rhs->keys[i] = lhs->keys[median_pos + 1 + i];
            // Shift the keys in rhs
            rhs->keys[i].value += lhs->_shift;
            rhs->keys[i].value -= rhs->_shift;

            lhs->keys[median_pos + 1 + i].value = 0;
            lhs->keys[median_pos + 1 + i].satellites = nullptr;
          }

        }else{

          for(B_t i = 0; i < offset - 1; ++i){
            rhs->keys[i] = lhs->keys[median_pos + 1 + i];
            rhs->children[i] = lhs->children[median_pos + 1 + i];

            // Shift the keys in rhs
            rhs->keys[i].value += lhs->_shift;
            rhs->keys[i].value -= rhs->_shift;
            // Shift the children in rhs
            rhs->children[i]->_shift += lhs->_shift;
            rhs->children[i]->_shift -= rhs->_shift;

            lhs->keys[median_pos + 1 + i].value = 0;
            lhs->keys[median_pos + 1 + i].satellites = nullptr;
            lhs->children[median_pos + 1 + i] = nullptr;
          }
          rhs->children[offset - 1] = lhs->children[lhs->n];
          // Shift the child in rhs
          rhs->children[offset -1]->_shift += lhs->_shift;
          rhs->children[offset -1]->_shift -= rhs->_shift;

          lhs->children[lhs->n] = nullptr;

        }

        // Place the median key in the i-th position
        keys[i] = lhs->keys[median_pos];
        // Shift the key
        keys[i].value += lhs->_shift;

        lhs->keys[median_pos].value = 0;
        lhs->keys[median_pos].satellites = nullptr;

        // Update the values of n
        lhs->n = median_pos;
        rhs->n = rhs->n + offset;

      }else{
        // The median key is in the right child

        // I have to make room on the leftt child for lhs->n - median_pos elements
        B_t offset = median_pos - lhs->n -1;

        // place the median_pos in the rhs offset
        median_pos -= (lhs->n + 1);


        // Copy the i-th key in the correct position
        lhs->keys[lhs->n] = keys[i];
        // Shift the key in lhs
        lhs->keys[lhs->n].value -= lhs->_shift;

        // Copy the first offset keys from rhs into lhs
        if(rhs->is_leaf()){

          for(B_t i = 0; i < offset; ++i){
            lhs->keys[lhs->n +1 + i] = rhs->keys[i];
            // Shift the key in lhs
            lhs->keys[lhs->n +1 + i].value += rhs->_shift;
            lhs->keys[lhs->n +1 + i].value -= lhs->_shift;
          }

        }else{

          for(B_t i = 0; i < offset; ++i){
            lhs->keys[lhs->n + 1 + i] = rhs->keys[i];
            lhs->children[lhs->n + 1 + i] = rhs->children[i];
            // Shift the keys in lhs
            lhs->keys[lhs->n +1 + i].value += rhs->_shift;
            lhs->keys[lhs->n +1 + i].value -= lhs->_shift;
            // Shift the children in lhs
            lhs->children[lhs->n +1 + i]->_shift += rhs->_shift;
            lhs->children[lhs->n +1 + i]->_shift -= lhs->_shift;

          }
          lhs->children[lhs->n + 1 + offset] = rhs->children[offset];
          // Shift the child in lhs
          lhs->children[lhs->n +1 + offset]->_shift += rhs->_shift;
          lhs->children[lhs->n +1 + offset]->_shift -= lhs->_shift;

        }

        // Place the median key in the i-th position
        keys[i] = rhs->keys[median_pos];
        // Shift the key
        keys[i].value += rhs->_shift;

        rhs->keys[median_pos].value = 0;
        rhs->keys[median_pos].satellites = nullptr;

        // shift offset + 1 elements of rhs to the left
        if(rhs->is_leaf()){

          for(B_t i = 0; i < (rhs->n - offset - 1); ++i){
            rhs->keys[i] = rhs->keys[i + offset + 1];

            rhs->keys[i + offset + 1].value = 0;
            rhs->keys[i + offset + 1].satellites = nullptr;
          }

        }else{

          for(B_t i = 0; i < (rhs->n - offset - 1); ++i){
            rhs->keys[i] = rhs->keys[i + offset + 1];
            rhs->children[i] = rhs->children[i + offset + 1];

            rhs->keys[i + offset + 1].value = 0;
            rhs->keys[i + offset + 1].satellites = nullptr;
            rhs->children[i + offset + 1] = nullptr;
          }

          rhs->children[rhs->n - offset - 1] = rhs->children[rhs->n];
          rhs->children[rhs->n] = nullptr;

        }

        // Update the values of n
        lhs->n = lhs->n + offset + 1;
        rhs->n = rhs->n - offset - 1;


      }
    }


  }

  /*!
   * Join the lhs subtree on the left spine of this tree, using the max_key
   * as pivot. The max_key must not be in the lhs subtree.
   * @param lhs     the subtree to be joint
   * @param max_key the pivot element
   * @param h_this    the height of the tree
   * @param h_lhs     the height of the subtree to be joint
   */
 template< typename K, typename S, size_t B, size_t T>
 void B_tree_node<K,S,B,T>::join_left(B_tree_node<K,S,B,T>* lhs, key_t max_key , size_t *h_this, size_t *h_lhs)
  {
    B_tree_node<K,S,B,T>* t1 = this;
    B_tree_node<K,S,B,T>* t2 = lhs;

    // Check if the heads are not null
    if(t2 == nullptr){
      return;
    }

    // Get the heights of the two trees.
    K h1, h2;
    if(h_this != nullptr)
      h1 = *h_this;
    else
      h1 = t1->height();

    if(h_lhs != nullptr)
      h2 = *h_lhs;
    else
      h2 = t2->height();

    assert(h2 < h1);

    // Test if the root is full and if so, split it.
    if(t1->is_full()){
      B_tree_node<K,S,B,T>* tmp = new B_tree_node<K,S,B,T>(false);
      // std::swap(*this,tmp);
      for(B_t i = 0; i < B; ++i){
        std::swap(keys[i], tmp->keys[i]);
      }
      std::swap(children, tmp->children);
      std::swap(n, tmp->n);
      // shift+= tmp->shift;

      children[0] = tmp;
      this->split_child(0);

      h1++; (*h_this)++;
    }

    // From this point on we have that h2 < h1.

    // Find the node on the left spine of t1 at height (h1 - h2)
    while( h1 > h2 + 1 ){

      B_tree_node<K,S,B,T>* t1_child = t1->children[0];
      // if the child is full, split it
      if(t1_child->is_full()){
        t1->split_child(0);
      }
      // Shift the tree t2 by -t1->shift;
      t2->_shift -= t1->_shift;
      max_key.value -= t1->_shift;

      t1 = t1_child;
      h1--;
    }

    // t1 has haight h2 + 1

    // Makes room for an extra element in front of t1
    t1->shift_right(0,1);
    // Put the max element of t1 as extra key of t1
    t1->keys[0] = max_key;
    // Shift the key
    t1->keys[0].value -= t1->_shift;

    // Put t2 as the leftmost child of t1
    t1->children[0] = t2;
    // Shift the child
    t2->_shift -= t1->_shift;

    // Fuse t1_child and t2
    t1->fuse_children(0);



  }


  /*!
   * Join the lhs subtree on the right spine of this tree, using the min_key
   * as pivot. The min_key must not be in the rhs subtree.
   * @param rhs     the subtree to be joint
   * @param min_key the pivot element
   * @param h_this    the height of the tree
   * @param h_rhs     the height of the subtree to be joint
   */
 template< typename K, typename S, size_t B, size_t T>
 void B_tree_node<K,S,B,T>::join_right(B_tree_node<K,S,B,T>* rhs, key_t min_key , size_t *h_this, size_t *h_rhs)
  {
    B_tree_node<K,S,B,T>* t1 = this;
    B_tree_node<K,S,B,T>* t2 = rhs;

    // Check if the heads are not null
    if(t2 == nullptr){
      return;
    }

    // Get the heights of the two trees.
    K h1, h2;
    if(h_this != nullptr)
      h1 = *h_this;
    else
      h1 = t1->height();

    if(h_rhs != nullptr)
      h2 = *h_rhs;
    else
      h2 = t2->height();

    assert(h1 >= h2);

    // From this point on we have that h1 >= h2.

    // Test if the root is full and if so, split it.
    if(t1->is_full()){
      B_tree_node<K,S,B,T>* tmp = new B_tree_node<K,S,B,T>(false);
      // std::swap(*this,tmp);
      for(B_t i = 0; i < B; ++i){
        std::swap(keys[i], tmp->keys[i]);
      }
      std::swap(children, tmp->children);
      std::swap(n, tmp->n);

      children[0] = tmp;
      this->split_child(0);

      h1++;(*h_this)++;
    }

    // If t1 and t2 have the same heights
    if (h1 == h2){
      // Create a new node
      B_tree_node<K,S,B,T>* tmp = new B_tree_node<K,S,B,T>(false);
      // std::swap(*this,tmp);
      for(B_t i = 0; i < B; ++i){
        std::swap(keys[i], tmp->keys[i]);
      }
      std::swap(children, tmp->children);
      std::swap(n, tmp->n);
      std::swap(_shift, tmp->_shift);
      // Place the max element of T as the median key in the new node
      t1->keys[0] = min_key;
      // Shift the key
      t1->keys[0].value -= t1->_shift;
      // Place t1 as left child and t2 as right child of tmp
      t1->children[0] = tmp;
      t1->children[1] = t2;
      t1->n = 1;
      (*h_this)++;
      // Fuse t1 and t2
      t1->fuse_children(0, h_this);
    }else{
      // Find the node on the right spine of t1 at height (h1 - h2)
      while( h1 > h2 + 1 ){

        B_tree_node<K,S,B,T>* t1_child = t1->children[t1->n];
        // if the child is full, split it
        if(t1_child->is_full()){
          t1->split_child(t1->n);
        }
        // Shift the tree t2 by -t1->shift;
        t2->_shift -= t1->_shift;
        min_key.value -= t1->_shift;

        t1 = t1->children[t1->n];
        h1--;
      }

      // t1 has height h2 + 1

      // Put the max element of t1 as extra key of t1
      t1->keys[t1->n] = min_key;
      // Shift the key
      t1->keys[t1->n].value -= t1->_shift;

      t1->n ++;

      // Put t2 as the rightmost child of t1
      t1->children[t1->n] = t2;
      // Shift the child
      t2->_shift -= t1->_shift;

      // Fuse t1_child and t2
      t1->fuse_children(t1->n - 1, h_this);

    }




  }


  template< typename K, typename S, size_t B, size_t T>
  void B_tree_node<K,S,B,T>::join(B_tree_node<K,S,B,T>* other)
  {
    B_tree_node<K,S,B,T>* t1 = this;
    B_tree_node<K,S,B,T>* t2 = other;

    // Check if the heads are not null
    if(t2 == nullptr){
      return;
    }

    // Get max and min values
    K t1_max = t1->get_max();
    // K t2_max = t2->get_max();
    // K t2_min = t2->get_min();

    assert(t1_max < t2->get_min());

    // Remove the max element of t1 from t1
    key_t t1_max_key = t1->remove(t1_max);

    // Get the heights of the two trees.
    size_t h1 = t1->height();
    size_t h2 = t2->height();

    if(h1 >= h2){

      t1->join_right(t2,t1_max_key,&h1,&h2);

    }else{
      // From this point on we have that h1 < h2.
      // Swap the pointers of t1 and t2
      // (this is to be consistent that the result of the join is in *this)

      // std::swap(t1,t2);
      for(B_t i = 0; i < B; ++i){
        std::swap(t1->keys[i], t2->keys[i]);
      }
      std::swap(t1->children, t2->children);
      std::swap(t1->n, t2->n);
      std::swap(t1->_shift, t2->_shift);

      std::swap(h1, h2);

      t1->join_left(t2,t1_max_key,&h1,&h2);
    }

  }

  template< typename K, typename S, size_t B, size_t T>
  B_tree_node<K,S,B,T>* B_tree_node<K,S,B,T>::split(const K &value_, size_t *h_this, size_t *h_rhs)
  {
    // Shift the value
    K value = value_ - _shift;

    // If not computed, compute the height of the tree
    size_t tmp_h_this, tmp_h_rhs;
    if(h_this == nullptr){
      tmp_h_this = height();
      h_this = &tmp_h_this;
    }
    if(h_rhs == nullptr){
      h_rhs = &tmp_h_rhs;
    }
    *h_rhs = *h_this;

    // Finds the value in the current node
    size_t l = 0;
    size_t r = n;

    // Find the predecessor and successor of value
    while(r > l + 1){
      size_t mid = (l+r)/2;
      if(keys[mid].value <= value){
        l = mid;
      }else{
        r = mid;
      }
    }

    if(keys[l].value <= value) l++;
    // l contains the index of the child that has to be split.

    // split the current node around the element l
    B_tree_node<K,S,B,T>* rhs = new B_tree_node<K,S,B,T>(is_leaf());
    // Shift rhs
    rhs->_shift = _shift;

    // copy the elements on the right of l into the new node
    if(is_leaf()){

      for(B_t i = l; i < n; ++i){
          rhs->keys[i-l] = keys[i];

          keys[i].value = 0;
          keys[i].satellites = nullptr;
      }
      rhs->n = n-l;
      n = l;

    }else{
      // Disconnect the l-th child from this node
      B_tree_node<K,S,B,T>* lhs_child = children[l];
      children[l] = nullptr;

      // Split the l-th child
      // Get the height of the children subtree
      size_t h_sub_tree = (*h_this) -1;
      size_t h_rhs_sub_tree = (*h_this) -1;
      B_tree_node<K,S,B,T>* rhs_child = lhs_child->split(value, &h_sub_tree, &h_rhs_sub_tree);
      // Shift the rhs_child and the lhs_child
      rhs_child->_shift += _shift;
      lhs_child->_shift += _shift;



      // Popolate the rhs node
      if( l == n ){

        delete rhs;
        rhs = rhs_child;
        // Update the height of rhs
        (*h_rhs) = (h_rhs_sub_tree);

      }else if( h_rhs_sub_tree == ((*h_rhs) -1)){

        // rhs_child and the children of rhs have the same height
        for(B_t i = l; i < n; ++i){
          rhs->keys[i-l] = keys[i];
          rhs->children[i-l+1] = children[i+1];

          keys[i].value = 0;
          keys[i].satellites = nullptr;
          children[i+1]= nullptr;
        }
        rhs->n = n-l;

        // Connect the child
        rhs->children[0] = rhs_child;
        // Shift the rhs_child
        rhs_child->_shift -= _shift;

        rhs->fuse_children(0, h_rhs);

      }else{
        // rhs_child and the children of rhs have different heights

        // Remove the l-th key
        key_t min_rhs = keys[l];
        // Shift the key
        min_rhs.value += _shift;

        keys[l].value = 0;
        keys[l].satellites = nullptr;

        if(l == (B_t)(n-1)){
          // rhs is empty
          delete rhs;
          rhs = children[l+1];
          children[l+1] = nullptr;
          // Shift rhs
          rhs->_shift += _shift;
          // Update the height of rhs
          (*h_rhs) --;

          if((*h_rhs) == h_rhs_sub_tree){
            rhs_child->join_right(rhs,min_rhs,h_rhs,h_rhs);
            rhs = rhs_child;
          }else{
            assert((*h_rhs) > h_rhs_sub_tree);
            rhs->join_left(rhs_child,min_rhs,h_rhs,&h_rhs_sub_tree);
          }

        }else{
          // Popolate the rhs node
          rhs->n = n-l-1;
          rhs->children[0] = children[l+1];
          children[l+1] = nullptr;

          for(B_t i = l+1; i < n; ++i){
            rhs->keys[i-l-1] = keys[i];
            rhs->children[i-l] = children[i+1];

            keys[i].value = 0;
            keys[i].satellites = nullptr;
            children[i+1]= nullptr;
          }

          rhs->join_left(rhs_child,min_rhs,h_rhs,&h_rhs_sub_tree);
        }

      }







      // Popolate the lhs (this) node
      n = l;
      if(l == 0){
        // Swap this with its child
        // std::swap(*this,lhs_child);
        for(B_t i = 0; i < B; ++i){
          std::swap(keys[i], lhs_child->keys[i]);
        }
        std::swap(n, lhs_child->n);
        std::swap(children, lhs_child->children);
        std::swap(_shift, lhs_child->_shift);
        // _shift += lhs_child->_shift;
        lhs_child->children[0] = nullptr;
        delete lhs_child;
        // Update the height of this subtree
        (*h_this) = (h_sub_tree);
      }else{

        // Remove the (l-1)-th key
        key_t max_lhs = keys[l-1];
        // Shift key
        max_lhs.value += _shift;

        keys[l-1].value = 0;
        keys[l-1].satellites = nullptr;

        n = l-1;

        if(l == 1){
          // Swap this with its child
          // std::swap(*this,children[0]);
          B_tree_node<K,S,B,T>* other = children[0];
          for(B_t i = 0; i < B; ++i){
            std::swap(keys[i], other->keys[i]);
          }
          std::swap(n, other->n);
          std::swap(children, other->children);
          _shift += other->_shift;
          other->children[0] = nullptr;
          delete other;
          // Update the height of this subtree
          (*h_this)--;

        }

        join_right(lhs_child,max_lhs,h_this,&h_sub_tree);


      }



    }

    return rhs;
  }




} // md


#endif // MD_B_TREE_NODE__HH
