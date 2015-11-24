#include <map>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

//this code is just a test where the score matrix is calculate with arbitrary proteins

//calculate the maximum between 3 values

int max3(int a, int b, int c){
    int maxTot = max(a, b);
    maxTot = max(maxTot, c);
return maxTot;
}

//Calculate the score matrix, and print it, return de score

int alignment(vector<char> prt, vector<char> myPrt, int *indel, map<char, map<char, int> > *maptrix){ 
    int score = 0;
    vector<vector<int> > vecScore;
    //int maxScoreX = 0; position of the maximum score
    //int maxScoreY = 0;

    for (int i = 0; i<= myPrt.size(); i++){ //create an empty score matrix
        vector<int> row;
        for(int j = 0; j<= prt.size(); j++){
            row.push_back(0);
        }
        vecScore.push_back(row);
    }

    for (int i = 1; i<= myPrt.size(); i++){
            vector<int> row;
        for(int j = 1; j<= prt.size(); j++){
            int possUP = vecScore[i][j-1] + *indel;
            int possLEFT = vecScore[i-1][j] + *indel;
            int possUL = vecScore[i-1][j-1] + (*maptrix)[myPrt[i-1]][prt[j-1]];
            int maxUL = max(possUP, possLEFT); int maxZUL = max(0, possUL);
            vecScore[i][j] = max(maxUL, maxZUL); //calculate the maximum score form the possibilities
            
            if((i==myPrt.size() || j == prt.size()) && (vecScore[i][j] > score)){ //register the score
				score = vecScore[i][j];
				//maxScoreX = i; maxScoreY = j;
			}
			
            cout << vecScore[i][j] << " ";
        }
        cout << endl;
    }

	//attempt of printing the aligned proteins, not functional yet

    /*
    int i = maxScoreX; int j = maxScoreY;
    vector<char> ansOne;
    vector<char> ansTwo;
    
    if(i<myPrt.size()){ 
        for(int k = myPrt.size()-1; k>i; k--){
            ansOne.insert(ansOne.begin(),myPrt[k]);
            cout << myPrt[k];
        }
    }
    if(j<prt.size()){
        for(int l = myPrt.size()-1; l>j; l--){
            ansTwo.insert(ansTwo.begin(),prt[l]);
            cout << prt[l];
        }
    }
    while (i > 0 && j > 0){

        cout << "(" << i << ", " << j << ")  "; cout << vecScore[i][j] << " ";
        int currentMax = max3(vecScore[i-1][j],vecScore[i][j-1],vecScore[i-1][j-1]);
        if(vecScore[i-1][j-1]==currentMax){
            ansOne.insert(ansOne.begin(),myPrt[i-1]);
            ansTwo.insert(ansTwo.begin(),prt[j-1]);
            cout << "UL:(" << myPrt[i-1] << ", " << prt[j-1] << ")" << endl;
            i--; j--;
        }else if(vecScore[i-1][j]==currentMax){
            ansOne.insert(ansOne.begin(),myPrt[i-1]);
            ansTwo.insert(ansTwo.begin(),'-');
            cout << "Left:(" << myPrt[i-1] << ", " << "-" << ")" << endl;
            i--;
        }else{
            ansOne.insert(ansOne.begin(),'-');
            ansTwo.insert(ansTwo.begin(),prt[j-1]);
            cout << "UP:(" << "-" << ", " << prt[j-1] << ")" << endl;
            j--;
        }
    }
    for(int k = 0; k<ansOne.size(); k++){
        cout << ansOne[k];
    }

    cout << endl;
    for(int l = 0; l<ansTwo.size(); l++){
        cout << ansTwo[l];
    }*/

return score;
}

int main(){

	//creation of the indel and the possible letters

    map<char, map<char, int> > maptrix;
    vector<char> ltrs;
    ltrs.push_back('A'); ltrs.push_back('Z'); ltrs.push_back('E'); ltrs.push_back('R');
    int indel = -6;

	//creation opf the proteins to comparate (pseudo-random)

    vector<vector<char> > chains;
    for(int i = 0; i < 3; i++){
        vector<char> vec;
        cout << i+1 << ": ";
        for(int j = 0; j < 6; j++){
            char vecChr = ltrs[rand() % 4];
            vec.push_back(vecChr);
            cout << vecChr;
        }
        chains.push_back(vec);
        cout << endl;
    }
    
    //creation of our protein (pseudo-random)
    
    vector<char> myChain;
    cout << endl;
    cout << "my chain: ";
    for(int j = 0; j < 6; j++){
        char vecChr = ltrs[rand() % 4];
        myChain.push_back(vecChr);
        cout << vecChr;
    }
    cout << endl; cout << endl; cout << "my \"blosum\" matrix:" << endl; cout << endl;

    //srand(time(NULL)); allow real random numbers
    
    //create of the "blosum" matrix, just used for this code example
    for(int i = 0; i<ltrs.size(); i++){
        char chr = ltrs[i];
        map<char, int> maprow;
        maptrix[chr] = maprow;
        for(int j = 0; j<ltrs.size(); j++){
                if(j>i){
                    maptrix[chr][ltrs[j]] = rand() % 11 - 5;
                }
                else if(i==j){
                    maptrix[chr][ltrs[j]] = rand() % 11 + 5;
                }
                else{
                    maptrix[chr][ltrs[j]] = maptrix[ltrs[j]][chr];
                }
                cout << maptrix[chr][ltrs[j]] << ' ';
        }
        cout << endl;
    }
    cout << endl;

	//vector of the scores of the different sequences (not used here)

    vector<int> score;
    for(int i = 0; i < chains.size(); i++){
        score.push_back(0);
    }

    int finalScore = alignment(chains[0], myChain, &indel, &maptrix);
    cout << finalScore;

return 0;
}
