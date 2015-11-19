#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "LinkedList.h"
using namespace std;


bool cmpProteins(LinkedList p1, LinkedList p2) {
  Node *current_p1;
  Node *current_p2;
  current_p1 = p1.getHead();
  current_p2 = p2.getHead();
  while(current_p1 != NULL && current_p2 != NULL) {
    if(current_p1->x != current_p2->x) {
      return false;
    }
    current_p1 = p1.getNext();
    current_p2 = p2.getNext();
  }
  if(current_p1 == NULL && current_p2 == NULL)
  return true;
  else
  return false;
}

int main(int argc, char* argv[]) {
  int i = 0;
  int j = 0;
  string line;

  char e = '\0';
  char c;
  unsigned char *b;


  // Init the map
  map<char, int> encoding_table;
  ifstream encoding ("assets/encoding.txt");
  while(getline(encoding,line))
  {
    stringstream   linestream(line);
    string         value;
    bool is_key = true;
    char temp_c;
    while(getline(linestream,value,','))
    {
      if(is_key) {
        //cout << "Key: " << value << ", ";
        temp_c = value.at(0);
      }
      else {
        //cout << "value: " << value << endl;
        // Creating pair
        encoding_table[temp_c] = stoi(value);
      }
      is_key = !is_key;
    }
    //cout << "Created pair:" << encoding_table[temp_c] << endl;

  }
  encoding.close();


  // Read the protein
  LinkedList testProtein;
  ifstream current_protein ("assets/proteins/P00533.fasta");
  if (current_protein.is_open())
  {
    bool reading = false;
    while (current_protein.read((&c), 1))
    {
      if(c == '\n' && !reading) {
        reading = true;
      }
      if(reading) {
        if(c != '\n') {
          testProtein.append(encoding_table[c]);
        }
      }
    }
  }
  current_protein.close();

  cout << "Protein ("<< testProtein.count << " AA):\n";
  testProtein.printFromHead();

  // Init the BLOSUM matrix
  vector< vector<int> > matrix;
  map<char, map<char, int> > maptrix;
  char keys[24];


  ifstream current_matrix("assets/blosum/BLOSUM62.txt");
  if (current_matrix.is_open())
  {

    int k = 0;
    while(getline(current_matrix,line))
    {
      vector<int> row;
      map<char, int> maprow;
      if(line.at(0) != '#') {
        // cout << line << endl;
        // Here we only get lines not commented
        stringstream linestream(line);
        string val;

        int j = 0;
        while(getline(linestream, val, ' ')) {
          // If it is the first line, it contains the keys
          if(k == 0) {
            //cout << "First line" << endl;
            char x;
            stringstream valuestream(val);
            valuestream >> x;

            if (!valuestream.fail()) {
              maptrix[val.at(0)] = maprow;
              keys[j] = val.at(0);
              cout << "Inserting key" << j << ": " << keys[j] << endl;
              cout << "Inserted:" << val.at(0) << endl;

              j++;
            }
          }

          // Here we are at the beginning of the line
          stringstream valuestream(val);
          // First, we need to extract the first char:
          int x;
          valuestream >> x;

          if (!valuestream.fail()) {
            row.push_back(stoi(val));
            cout << "Inserted (" << keys[k-1] << ", " << keys[j] << "):" << val << endl;
            maptrix[keys[k-1]][keys[j]] = stoi(val);
            j++;
          }
        }
        cout << "end of line" << k << endl;
        k++;
      }
      matrix.push_back(row);
    }
  }

  cout << "Map A:G" << maptrix['P']['P'] << endl;
  current_matrix.close();


  // Print the matrix
  /*
  for (vector< vector<int> >::const_iterator j = matrix.begin(); j != matrix.end(); ++j) {
  for (vector<int>::const_iterator i = (*j).begin(); i != (*j).end(); ++i) {
  cout << *i << ' ';
}
cout << endl;
}
*/







// Read the database


cout << "\n\nReading database: \n" << endl;
ifstream myfile ("assets/db/uniprot_sprot.fasta.psq", ios::binary);
if (myfile.is_open())
{
  // Let us create a new linked list
  LinkedList dbProtein;
  while (myfile.read((&c), 1))
  {
    if(c == e) {
      if(i > 0) { // First byte of the file is an empty byte
        //cout << "\ndb:";
        //dbProtein.printFromHead();

        // Compare proteins
        map<char, map<char, int> > score;

        if(cmpProteins(dbProtein, testProtein)) {
          cout << "Match: " << i << endl;
          break;
        }
        else {
          //cout << "Mismatch" << endl;
        }

        dbProtein.deleteList();
      }
      i++;
    }
    else {
      dbProtein.append((uint8_t)c);
    }
    //cout << (int)c << ' ';
  }

  myfile.close();
}

else cout << "\nUnable to open file\n";

return EXIT_SUCCESS;
}
