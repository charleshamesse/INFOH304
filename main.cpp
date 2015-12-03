#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <thread>
#include <pthread.h>
#include "LinkedList.h"
#include "ThreadPool.h"
using namespace std;

/*
* cmpProteins: returns true if p1 is equal to p2, false if not.
* old version using LinkedLists

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
*/
/*
void* execThread(int i) {
cout << "Thread " << i << "<< says Hi!" << endl;
}
*/

// cmpProteins: new version using vectors
int cmpProteins(vector<int> p1, vector<int> p2, int indel, map<int, map<int, int> > *maptrix){
  // prt p1
  // myPrt p2
  int score = 0;
  vector<vector<int> > vecScore;

  for (int i = 0; i<= p2.size(); i++){ // Create an empty score matrix
    vector<int> row;
    for(int j = 0; j<= p1.size(); j++){
      row.push_back(0);
    }
    vecScore.push_back(row);
  }


  for (int i = 1; i<= p2.size(); i++){
    vector<int> row;
    // current_p2 = p2[i-1]
    for(int j = 1; j< p1.size(); j++){
      // current_p1 = p1[j-1]
      int possUP = vecScore[i][j-1] + indel;
      int possLEFT = vecScore[i-1][j] + indel;
      int possUL = vecScore[i-1][j-1] + (*maptrix)[p2[i-1]][p1[j-1]];
      int maxUL = max(possUP, possLEFT);
      int maxZUL = max(0, possUL);
      vecScore[i][j] = max(maxUL, maxZUL); //calculate the maximum score form the possibilities

      if((i==p2.size() || j == p1.size()) && (vecScore[i][j] > score)){ //register the score
        score = vecScore[i][j];
        //maxScoreX = i; maxScoreY = j;
      }

      //cout << ". Score: " << vecScore[i][j] << endl;
    }
  }

  return score;
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
  vector<int> vecTestProtein;
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
          vecTestProtein.push_back(encoding_table[c]);
        }
      }
    }
  }
  current_protein.close();
  // Print the protein
  cout << "Protein ("<< testProtein.count << " AA):\n";

  // Init the BLOSUM matrix
  vector< vector<int> > matrix;
  map<int, map<int, int> > maptrix;
  int keys[24];


  ifstream current_matrix("assets/blosum/BLOSUM62.txt");
  if (current_matrix.is_open())
  {

    int k = 0; // line counter
    // For each line in the file
    while(getline(current_matrix,line))
    {
      vector<int> row;
      map<int, int> maprow;
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
              int key = encoding_table[val.at(0)];
              maptrix[key] = maprow;
              keys[j] = key;
              //maptrix[val.at(0)] = maprow;
              //keys[j] = val.at(0);
              //cout << "Inserting key" << j << ": " << keys[j] << endl;
              //cout << "Inserted:" << key << endl;
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
            //cout << "Inserted (" << keys[k-1] << ", " << keys[j] << "):" << val << endl;
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
  //cout << "Map A:G" << maptrix['P']['P'] << endl;
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


/*

Thread example

*/


/*
int Num_Threads =  thread::hardware_concurrency();
cout << "Number of threads: " << Num_Threads << endl;
pthread_t threads[5];
for(size_t i = 0; i < sizeof threads / sizeof *threads; ++i) {
pthread_create(threads + i, NULL,execThread, i);
}

// Wait for the threads to terminate.
for(size_t i = 0; i < sizeof threads / sizeof *threads; ++i) {
pthread_join(threads[i], NULL);
}
*/

/*
*
* Init the thread pool
*
*/

int nThreads =  std::thread::hardware_concurrency();
std::cout << "Number of threads: " << nThreads << std::endl;
ThreadPool pool(nThreads);
std::vector< std::future<int> > results;
int current_score;




// Read the database
cout << "\n\nReading database: \n" << endl;
ifstream myfile ("assets/db/uniprot_sprot.fasta.psq", ios::binary);
if (myfile.is_open())
{
  // Let us create a new linked list
  LinkedList dbProtein;
  vector<int> vecDbProtein;
  // For each byte in the file
  while (myfile.read((&c), 1) && i < 50)
  {
    // If the current byte is the empty byte
    if(c == e) {
      // First byte of the file is an empty byte, so we do not add the protein to the list as it is empty.
      if(i > 0) {
        // We could print the protein using: dbProtein.printFromHead();
        // Compare proteins
        // TO-DO
        map<char, map<char, int> > score;

        // Give the job to the thread pool
        results.emplace_back(
          pool.enqueue([current_score, vecDbProtein, vecTestProtein, &maptrix] {
            return cmpProteins(vecDbProtein, vecTestProtein, -1, &maptrix);
          })
        );

        /* Without multi-threading:
        int current_score = cmpProteins(vecDbProtein, vecTestProtein, -1, &maptrix);
        cout << "Score " << i << ": " << current_score << endl;
        */

        // We delete the list containing the protein and start over.
        dbProtein.deleteList();
        vecDbProtein.clear();
      }
      i++;
    }
    // If the byte is not empty, we append it to the current protein.
    else {
      dbProtein.append((uint8_t)c);
      vecDbProtein.push_back((uint8_t)c);
    }
  }

  myfile.close();
  // Displaying scores
  int nScore = 1;
  for(auto && result: results) {
    std::cout << "Score " << nScore << ": " << result.get() << endl;
    nScore++;
  }
  std::cout << std::endl;
}

else cout << "\nUnable to open file\n";



return EXIT_SUCCESS;
}
