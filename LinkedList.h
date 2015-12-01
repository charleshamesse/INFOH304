struct Node {
  int x;
  Node *next;
  Node *prev;
};

class LinkedList{
private:
  Node *head; // this is the private member variable. It is just a pointer to the first Node
  Node *tail; // this is the private member variable. It is just a pointer to the first Node
public:
  // constructor
  long int count; // Number of nodes
  long int idx; // Current position in LinkedList
  Node *current;

  LinkedList();

  long int getCount();

  Node* getNext();
  Node* getPrev(); // To be corrected
  Node* getHead();

  // This prepends a new value at the beginning of the list
  void prepend(int val);

  // This prepends a new value at the beginning of the list
  void append(int val);

  // returns the first element in the list and deletes the Node.
  // caution, no error-checking here!
  int popValue();

  void printFromHead();
  void printFromTail();
  void deleteList();

};
