#include <iostream>
#include "LinkedList.h"
using namespace std;

LinkedList::LinkedList(){
  head = NULL; // set head to NULL
  tail = NULL;
  current = NULL;
  count = 0;
  idx = 0;
}

long int LinkedList::getCount() {
  return count;
}


Node* LinkedList::getNext() {
  // There has to be at least one node
  if(count > 0) {
    // If getNext(); has alreadby been called
    if(idx > 0) {
      current = current->next;
      idx++;
    }
    // Else, return the first element
    else {
      current = head;
      idx++;
    }
    return current;
  }
  else {
    return NULL;
  }
}

Node* LinkedList::getPrev() {
  // There has to be at least one node
  if(count > 0) {
    // If getNext(); has alreadby been called
    if(idx > 0) {
      current = current->prev;
      idx--;
    }
    // Else, return the first element
    else {
      current = head;
      //idx--;
    }
    return current;
  }
  else {
    return NULL;
  }
}

Node* LinkedList::getHead() {
  return head;
}

// This prepends a new value at the beginning of the list
void LinkedList::prepend(int val){
  //cout << "(" << count << ',' << val << "), ";
  Node *n = new Node();   // create new Node
  n->x = val;             // set value
  n->prev = NULL;
  //  If the list is empty, this is NULL, so the end of the list --> OK

  if(count > 0) {
    n->next = head;         // make the node point to the next node.
    head->prev = head;
  }
  // If we add the first element
  else {
    n->next = NULL;
    tail = n;
  }
  head = n;
  count++;
}

// This prepends a new value at the beginning of the list
void LinkedList::append(int val){
  //cout << "(" << count << ',' << val << "), ";
  Node *n = new Node();   // create new Node
  n->x = val;             // set value
  n->next = NULL;         // make the node point to the next node.
  if(count > 0) {
    n->prev = tail;
    tail->next = n;
  }
  // If we add the first element
  else {
    n->prev = NULL;
    head = n;
  }
  tail = n;
  count++;
}

// returns the first element in the list and deletes the Node.
// caution, no error-checking here!
int LinkedList::popValue(){
  Node *n = head;
  int ret = n->x;

  head = head->next;
  delete n;
  return ret;
}

void LinkedList::printFromHead(){
  Node *temp = head;
  cout << "\n-- Beginning of protein --" << endl;
  while(temp != NULL) {
    cout << temp->x << ", ";
    temp = temp->next;
  }
  cout << "\n-- End of protein --" << endl;
}

void LinkedList::printFromTail(){
  // To check
  Node *temp;
  while(tail != NULL) {
    temp = tail;
    delete tail;
    cout << temp->x << ' ';
    tail = temp->prev;
  }
}
void LinkedList::deleteList(){
  Node *temp;
  while(head != NULL) {
    temp = head->next;
    delete(head);
    count--;
    head = temp;
  }
  // temp should be NULL and deleted when going out of scope
  //delete(temp);
}
