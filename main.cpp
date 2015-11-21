#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "LinkedList.h"
using namespace std;

/*
 * cmpProteins: returns true if p1 is equal to p2, false if not.
 */
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

/* For now, the main function does the following actions:
 * 1. Read the encoding table
 * 2. Read the input protein
 * 3. Read the BLOSUM matrix
 * 4. Read the database and test for equality
 */
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
  // For each line in the file
  while(getline(encoding,line))
  {
    stringstream   linestream(line);
    string         value;
    bool is_key = true;
    char temp_c;
    // Separate the key from the value (each line is formatted as key, value such as A,10)
    while(getline(linestream,value,','))
    {
      // The following is always true at the beginning of the line
      if(is_key) {
        // We store the key temporarily (e.g. 'A')
        temp_c = value.at(0);
      }
      else {
        // That means we're on the value (e.g. '10')
        encoding_table[temp_c] = stoi(value);
      }
      is_key = !is_key;
    }
    // Here we finish the line with the pair created.
  }
  encoding.close();


  // Read the protein
  LinkedList testProtein;
  ifstream current_protein ("assets/proteins/P00533.fasta");
  if (current_protein.is_open())
  {
    bool reading = false;
    // We read the protein byte per byte
    // The first line of fasta files is to be ignored. For this, we create a boolean value which is false during the first line and becomes true as from the second line.
    while (current_protein.read((&c), 1))
    {
      // If we encounter the first line break
      if(c == '\n' && !reading) {
        reading = true;
      }
      if(reading) {
        // We add each char to the list, except line breaks.
        if(c != '\n') {
          testProtein.append(encoding_table[c]);
        }
      }
    }
  }
  current_protein.close();
  // Print the protein
  cout << "Protein ("<< testProtein.count << " AA):\n";
  testProtein.printFromHead();

  // Init the BLOSUM matrix
  vector< vector<int> > matrix;
  map<char, map<char, int> > maptrix;
  char keys[24];


  ifstream current_matrix("assets/blosum/BLOSUM62.txt");
  if (current_matrix.is_open())
  {

    int k = 0; // line counter
    // For each line in the file
    while(getline(current_matrix,line))
    {
      vector<int> row;
      map<char, int> maprow;
      // If the line is a comment, we do nothing.
      if(line.at(0) != '#') {
        stringstream linestream(line);
        string val;

        int j = 0; // char counter
        // We then separate the current at each white space
        while(getline(linestream, val, ' ')) {
          // If it is the first line, it contains the keys. k is the line counter.
          if(k == 0) {
            // So here's the first line only.
            // The following code tries to fit the read value into a char.
            char x;
            stringstream valuestream(val);
            valuestream >> x;

            // If it did not fail (we're reading a char), we're reading the keys. Thus we store them.
            if (!valuestream.fail()) {
              maptrix[val.at(0)] = maprow;
              keys[j] = val.at(0);
              cout << "Inserting key" << j << ": " << keys[j] << endl;
              cout << "Inserted:" << val.at(0) << endl;
              j++; // char counter
            }
          }

          // Every line
          // The idea is the same as above, with an int.
          stringstream valuestream(val);
          int x;
          valuestream >> x;

          if (!valuestream.fail()) {
            row.push_back(stoi(val));
            cout << "Inserted (" << keys[k-1] << ", " << keys[j] << "):" << val << endl;
            maptrix[keys[k-1]][keys[j]] = stoi(val);
            // End of int
            j++;
          }
        }
        // End of line
        k++;
      }
      matrix.push_back(row);
    }
  }

  // Little test:
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
  // For each byte in the file
  while (myfile.read((&c), 1))
  {
    // If the current byte is the empty byte
    if(c == e) {
      // First byte of the file is an empty byte, so we do not add the protein to the list as it is empty.
      if(i > 0) {
        // We could print the protein using: dbProtein.printFromHead();
        // Compare proteins
        // TO-DO {
        map<char, map<char, int> > score;

        if(cmpProteins(dbProtein, testProtein)) {
          cout << "Match: " << i << endl;
          break;
        }
        else {
          //cout << "Mismatch" << endl;
        }
        //}

        // We delete the list containing the protein and start over.
        dbProtein.deleteList();
      }
      i++;
    }
    // If the byte is not empty, we append it to the current protein.
    else {
      dbProtein.append((uint8_t)c);
    }
  }

  myfile.close();
}

else cout << "\nUnable to open file\n";

return EXIT_SUCCESS;
}
