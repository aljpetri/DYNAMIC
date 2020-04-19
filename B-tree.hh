////////////////////////////////////////////////////////////////////////////////
// B-tree.hh
//   B-tree header file.
//
// representation of B-tree nodes, which is used to store the Minimizers
//
////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019 Simon J. Puglisi and Massimiliano Rossi altered by Alexander Petri
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


#ifndef MD_B_TREE__HH
#define MD_B_TREE__HH

#include <vector>
#include <assert.h>
#include <stack>
#include "B_tree_node.hh"
namespace md{

  /*!
   * T is the type of the values
   * S is the type of the satellites
   * B is the number of the pivots
   */
  template< typename K, typename S, size_t B = 7, size_t T = 3>
  class B_tree{
  public:
    typedef typename std::conditional<in_range_unsigned<uint8_t>(B),uint8_t,
                              typename std::conditional<in_range_unsigned<uint16_t>(B), uint16_t ,
                                typename std::conditional<in_range_unsigned<uint32_t>(B), uint32_t ,
                                  uint64_t
                                >::type
                              >::type
                            >::type B_t;

    typedef typename B_tree_node<K,S,B,T>::key_t key_t;
    typedef typename B_tree_node<K,S,B,T>::shifted_key_ptr_t shifted_key_ptr_t;

    typedef struct tt_key{
      K value;
      std::vector<S> satellites;

      ~tt_key(){
        // NtD
      }

      tt_key():
        value(0),satellites(0)
      {

      }

      tt_key(key_t other):
        value(other.value),satellites(*other.satellites)
      {

      }

    } key_tt;

    B_tree();
    ~B_tree();

    /*!
     * Create a new set and insert a new element.
     * @param value_     the value of the element in the set.
     * @param satellite_ the satellite attached to the element in the set.
     */
    shifted_key_ptr_t make_set(K &value_, S &satellite_);
    /*!
     * Insert a new element in the set.
     * @param value_     the value of the element to be inserted.
     * @param satellite_ the satellite attached to the element to be inserted.
     */
    shifted_key_ptr_t insert(K &value_, S &satellite_);
    /*!
     * Remove the element of value value_ from the set.
     * @param value_ the value of the element to be removed.
     */
    key_tt remove(K &value_);
    /*!
     * Finds if the element y = value_ exists in the set
     * @param  value_ the value of the element y
     * @return        the element y if it exists, nullptr otherwise.
     */
    shifted_key_ptr_t search(const K &value_);
    /*!
     * Shifts all keys being greater or equal than value_ by shift_
     * @param  value_ the value of the element y
     * @param  shift_ the shift to be applied
     * @return        the element y if it exists, nullptr otherwise.
     */
    shifted_key_ptr_t shift_greater(const K &value_, K &shift_);

    shifted_key_ptr_t predecessor(const K &value_);
    shifted_key_ptr_t successor(const K &value_);


    B_tree<K,S,B,T>* shift(K &shift_);
    B_tree<K,S,B,T>* join(B_tree<K,S,B,T>* rhs);
    B_tree<K,S,B,T>* split(const K &value_);
    B_tree<K,S,B,T>* merge(B_tree<K,S,B,T>* rhs);

    // Iterator: return one bitstring at a time
   class iterator{
   public:
     typedef iterator self_type;
     typedef std::pair<K,std::vector<S> > value_type;
     typedef B_tree_node<K,S,B,T>* pointer;
     /*!
      * Create an iterator that start from the first or from the last element.
      * @param current_   the head of the B-tree
      */
     iterator(pointer current_) : _current(current_), _shift(0) {
       if(_current != nullptr){
         stack.push(std::make_pair(nullptr,std::make_pair(0,0)));
         while(!_current->is_leaf() && _current->get_child(0) != nullptr ){
           // Find the leftmost child
           stack.push(std::make_pair(_current,std::make_pair(0,_shift)));
           _shift += _current->get_shift();
           _current=_current->get_child(0);
         }
         _current_index = 0;
         stack.push(std::make_pair(_current,std::make_pair(_current_index,_shift) ));
       }
     }

     self_type operator++() {
       self_type i = *this;
       internal_inc();
       return i;
     }

     self_type operator++(int junk) {
       (void)junk;
       internal_inc();
       return *this;
     }
     const value_type operator->() {
       shifted_key_ptr_t key_ptr = _current->get_key(_current_index);
       K value =  key_ptr.key->value + key_ptr.shift + _shift;
       return std::make_pair((value),*(key_ptr.key->satellites));
     }
     const value_type operator*() {
       shifted_key_ptr_t key_ptr = _current->get_key(_current_index);
       K value =  key_ptr.key->value + key_ptr.shift + _shift;
       return std::make_pair((value),*(key_ptr.key->satellites));
     }
     bool operator==(const self_type& rhs) { return _current == rhs._current; }
     bool operator!=(const self_type& rhs) { return _current != rhs._current; }
   protected:

     void internal_inc()
     {
       _current = stack.top().first;
       _current_index = stack.top().second.first;
       _shift = stack.top().second.second;
       stack.pop();
       bool last = false; // Is the last element of the subtree rooted in current
       if(_current == nullptr) return;

       _current_index++;
       if(_current_index <_current->get_n_keys()){

         stack.push(std::make_pair(_current,std::make_pair(_current_index,_shift)));

       }else{
         last = true;
       }

       if(!_current->is_leaf() && _current->get_child(_current_index)!= nullptr){
         _shift += _current->get_shift();
         _current=_current->get_child(_current_index);

         // Find the leftmost child
         _current_index = 0;
         while(!_current->is_leaf() && _current->get_child(_current_index) != nullptr ){
           stack.push(std::make_pair(_current,std::make_pair(0,_shift)));
           _shift += _current->get_shift();
           _current=_current->get_child(_current_index);
         }
         stack.push(std::make_pair(_current,std::make_pair(_current_index,_shift)));

       }else if(last){
         while(_current != nullptr && _current_index >=_current->get_n_keys()){
           _current = stack.top().first;
           _current_index = stack.top().second.first;
           _shift = stack.top().second.second;
           stack.pop();
         }
         stack.push(std::make_pair(_current,std::make_pair(_current_index,_shift)));

       }

     }

   private:
     std::stack<std::pair<pointer,std::pair<B_t,K>> > stack;
     pointer _current;
     B_t _current_index;
     K _shift;
   };


   iterator begin() {
     return iterator(_head);
   }
   iterator end(){
     return iterator(nullptr);
   }


   // utils to get max and min values
   K get_max(){
       return _head->get_max();
   }

   K get_min(){
       return _head->get_min();
   }

   bool is_empty(){
     return (_head == NULL);
   }

   bool check_integrity(){
     return _head->check_integrity();
   }

  private:
    B_tree_node<K,S,B,T>* _head;
  }; // B_tree

  // Ctor
  template< typename K, typename S, size_t B, size_t T>
  B_tree<K,S,B,T>::B_tree():
    _head(nullptr)
  {
    assert(T<= B && T > 1);
  }

  // Dtor
  template< typename K, typename S, size_t B, size_t T>
  B_tree<K,S,B,T>::~B_tree()
  {
    if(_head != nullptr)
      delete _head;
  }

  /*!
   * Create a new set and insert a new element.
   * @param value_     the value of the element in the set.
   * @param satellite_ the satellite attached to the element in the set.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::make_set(K &value_, S &satellite_)
  {
    _head = new B_tree_node<K,S,B,T>();
    return _head->insert(value_,satellite_);
  }
  /*!
   * Insert a new element in the set.
   * @param value_     the value of the element to be inserted.
   * @param satellite_ the satellite attached to the element to be inserted.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::insert(K &value_, S &satellite_)
  {
    if(_head == nullptr) return make_set(value_,satellite_);

    if(_head->is_full()){
      B_tree_node<K,S,B,T>* tmp = new B_tree_node<K,S,B,T>(false);
      tmp->set_child(0,_head);
      tmp->split_child(0);
      _head = tmp;
    }

    return _head->insert(value_,satellite_);
  }
  /*!
   * Remove the element of value value_ from the set.
   * @param value_ the value of the element to be removed.
   */
  template< typename K, typename S, size_t B, size_t T>
  typename B_tree<K,S,B,T>::key_tt B_tree<K,S,B,T>::remove(K &value_)
  {
    key_t res = _head->remove(value_);
    if(_head->get_n_keys() == 0){
      delete _head;
      _head = nullptr;
    }

    if(res.satellites != nullptr){
      key_tt res_(res);

      delete res.satellites;

      return res_;
    }else{
      key_tt res_;
      return res_;
    }

  }

  /*!
   * Shifts all keys being greater or equal than value_ by shift_
   * @param  value_ the value of the element y
   * @param  shift_ the shift to be applied
   * @return        the element y if it exists, nullptr otherwise.
   */
   template< typename K, typename S, size_t B, size_t T>
   typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::shift_greater(const K &value_,K &shift_)
   {
     if(_head != nullptr)
       return _head->shift_greater(value_,shift_);
     else
       return shifted_key_ptr_t();
   }
  /*!
   * Finds if the element y = value_ exists in the set
   * @param  value_ the value of the element y
   * @return        the element y if it exists, NULL otherwise.
   */
   template< typename K, typename S, size_t B, size_t T>
   typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::search(const K &value_)
   {
     if(_head != nullptr)
       return _head->search(value_);
     else
       return shifted_key_ptr_t();
   }

   template< typename K, typename S, size_t B, size_t T>
   typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::predecessor(const K &value_)
   {
     if(_head != nullptr)
       return _head->predecessor(value_);
     else
       return shifted_key_ptr_t();
   }

   template< typename K, typename S, size_t B, size_t T>
   typename B_tree<K,S,B,T>::shifted_key_ptr_t B_tree<K,S,B,T>::successor(const K &value_)
   {
     if(_head != nullptr)
       return _head->successor(value_);
     else
       return shifted_key_ptr_t();
   }

   template< typename K, typename S, size_t B, size_t T>
   B_tree<K,S,B,T>* B_tree<K,S,B,T>::shift(K &shift_)
   {
     if(_head != nullptr)
       _head->shift(shift_);
     return this;
   }


  template< typename K, typename S, size_t B, size_t T>
  B_tree<K,S,B,T>* B_tree<K,S,B,T>::join(B_tree<K,S,B,T>* rhs)
  {
    if(_head == nullptr){
      _head = rhs->_head;
    }else if(rhs->_head == nullptr){
      // NTD
    }else{
      _head->join(rhs->_head);
    }

    rhs->_head = nullptr;
    delete rhs;
    return this;

  }

  template< typename K, typename S, size_t B, size_t T>
  B_tree<K,S,B,T>* B_tree<K,S,B,T>::split(const K &value_)
  {
    B_tree<K,S,B,T>* rhs = new B_tree<K,S,B,T>();

    if(_head != nullptr){
      rhs->_head = _head->split(value_);
      if(rhs->_head->get_n_keys() == 0){
        delete rhs->_head;
        rhs->_head = nullptr;
      }
      if(_head->get_n_keys() == 0){
        delete _head;
        _head = nullptr;
      }
    }

    return rhs;

  }

  template< typename K, typename S, size_t B, size_t T>
  B_tree<K,S,B,T>* B_tree<K,S,B,T>::merge(B_tree<K,S,B,T>* rhs)
  {
    // ******************************************
    // COMMON PART
    // ******************************************
    if(_head == nullptr){
      _head = rhs->_head;
      rhs->_head = nullptr;
    }
    if(rhs->_head == nullptr){
      delete rhs;
      return this;
    }
    // ******************************************
    // Extract and Insert merge algorithm
    // ******************************************
    //
    // for(auto elem: (*rhs)){
    //   // insert(elem.first,elem.second);
    //   for(auto sat: (elem.second))
    //     insert(elem.first,sat);
    // }
    //
    // // ToDo: check if it is necessary to connect the elements to a new "pointer"
    // // Free space
    // delete rhs;
    // return this;

    // ******************************************
    // Farrach & Thorup merge algorithm
    // ******************************************
    B_tree<K,S,B,T>* A = new B_tree<K,S,B,T>();
    B_tree<K,S,B,T>* D = new B_tree<K,S,B,T>();
    B_tree<K,S,B,T>* C = new B_tree<K,S,B,T>();

    A->_head = _head;
    D->_head = rhs->_head;
    do{
      K min_A = A->_head->get_min();
      K min_D = D->_head->get_min();
      if(min_D < min_A){
        std::swap(A,D);
        std::swap(min_A,min_D);
      }
      B_tree<K,S,B,T>* A_I = A;//new B_tree<K,S,B,T>();
      A_I->_head = A->_head;
      A = A_I->split(min_D);
      // perform the pruning
      B_tree<K,S,B,T>* eq_elem = A_I->split(min_D-1);
      if(eq_elem->_head != nullptr){
        shifted_key_ptr_t min_elem = D->_head->predecessor(min_D);
        for(auto sat: *(eq_elem->_head->get_key(0).key->satellites)){
          min_elem.key->satellites->push_back(sat);
        }
      }
      delete eq_elem;

      C->join(A_I);
    }while(A->_head != nullptr && D->_head != nullptr);

    C->join(A);
    C->join(D);

    _head = C->_head;

    C->_head = nullptr;
    delete C;

    rhs->_head=nullptr;
    delete rhs;
    return this;
  }



} // md


#endif // MD_B_TREE__HH
