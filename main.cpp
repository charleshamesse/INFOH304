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

// cmpProteins: version using vectors 
int cmpProteins(vector<int> p1, vector<int> p2, int indel, int openIndel, map<int, map<int, int> > *maptrix){
  int nbInd = 0;
  int score = 0;
  vector<vector<int> > vecScore;
  vector<vector<int> > vecIndel;

  for (int i = 0; i<= p2.size(); i++){ // Create an empty score matrix and an empty indel matrix
    vector<int> row;
    vector<int> row2;
    for(int j = 0; j<= p1.size(); j++){
      row.push_back(0);
      row2.push_back(0);
    }
    vecScore.push_back(row);
    vecIndel.push_back(row2);
  }
  for (int i = 1; i<= p2.size(); i++){
    vector<int> row;
    for(int j = 1; j< p1.size(); j++){
	  int possUP = 0; int possLEFT = 0;
	  if(vecIndel[i][j-1] == 1){
		possUP = vecScore[i][j-1] + indel;
	  }else{
		possUP = vecScore[i][j-1] + indel + openIndel;
	  }
	  if(vecIndel[i-1][j] == 2){
		possLEFT = vecScore[i-1][j] + indel;
	  }else{
		possLEFT = vecScore[i-1][j] + indel + openIndel;
	  }
      int possUL = vecScore[i-1][j-1] + (*maptrix)[p2[i-1]][p1[j-1]];
      int maxUL = max(possUP, possLEFT);
      int maxZUL = max(0, possUL);
      vecScore[i][j] = max(maxUL, maxZUL); //calculate the maximum score from the possibilities
      if (max(maxUL, maxZUL) == possUP && max(maxUL, maxZUL) != 0){
		  vecIndel[i][j] = 1;
	  }
	  else if(max(maxUL, maxZUL) == possLEFT && max(maxUL, maxZUL) != 0){
		  vecIndel[i][j] = 2;
	  }
      if((i==p2.size() || j == p1.size()) && (vecScore[i][j] > score)){ //register the score
        score = vecScore[i][j];
      }
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

  map<char, int> encoding_table;  // Init the map
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
  current_matrix.close();
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
  while (myfile.read((&c), 1))// && i < 5)
  {
    // If the current byte is the empty byte
    if(c == e) {
      // First byte of the file is an empty byte, so we do not add the protein to the list as it is empty.
      if(i > 0) {
        // Compare proteins
        map<char, map<char, int> > score;

        // Give the job to the thread pool
        results.emplace_back(
          pool.enqueue([current_score, vecDbProtein, vecTestProtein, &maptrix] {
            return cmpProteins(vecDbProtein, vecTestProtein, -1, -6, &maptrix);
          })
        );
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
  int nScore = 0;
  for(auto && result: results) {
	  if(nScore%1000 == 0){
		std::cout << "Score " << nScore << ": " << result.get() << endl;
	  }
    nScore++;
  }
  std::cout << std::endl;
}

else cout << "\nUnable to open file\n";

return EXIT_SUCCESS;
}
