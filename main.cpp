/*

Copyright (c) 2011 - 2012, Brno University of Technology, Antonínská 548/1,
Brno 601 90, Czech Republic

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, 
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of the Brno University of Technology nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@author David Zouhar <xzouha02@stud.fit.vutbr.cz>

*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <cstdlib>

#include "classifier.h"

#define NUM_STATES	9
#define FEATURE_WIDTH	2
#define NUM_TRAIN_DATA	100
#define NUM_TEST_DATA	50

using namespace std;

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void loadTrainData(double*** data, int featureWidth, int* dataNum, int numTrainSet);

vector<vector<double *> > loadTrainData(string dirWithTrainData, int numTrainSet, int P);

vector<double *> loadTestData(string filePath);


/**
  * MAIN
  */
int main (int argc, char *argv[]) 
{
	string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
	//** TRENOVANI
	for(unsigned int i = 0; i < alphabet.size(); i++) {
		
		if(i != alphabet.find("T")) {
			continue;
		}
		
		cout << alphabet[i] << " -------------------------------------------------------" << endl;
		string characterFolder = "data/";
		characterFolder += alphabet[i];
		
		vector<vector<double *> > trainData = loadTrainData(characterFolder, 100, FEATURE_WIDTH);
		Model m(NUM_STATES, FEATURE_WIDTH);
		m.init_model(trainData[0]);
		
		m.train(trainData);
		
		ofstream file;
		string modelPath = "models/";
		modelPath += alphabet[i];
		modelPath += ".model";
		file.open(modelPath.c_str());
		file << m << endl;
		file.close();
	}
	//*/
	
	
	/** ROZPOZNAVANI
	
	Classifier classifier(FEATURE_WIDTH);
	classifier.addModels("models", alphabet);
	
	int confMatrixCount[alphabet.size()][alphabet.size()];
	for(unsigned int i = 0; i < alphabet.size(); i++) {
		for(unsigned int j = 0; j < alphabet.size(); j++) {
			confMatrixCount[i][j] = 0;
		}
	}
	
	for(unsigned int i = 0; i < alphabet.size(); i++) {
		
		cout << alphabet[i] << " -> ";
		for(int j = 0; j < NUM_TEST_DATA; j++) {
			vector<double *> testData;
			string dataFile = "test/";
			dataFile += alphabet[i];
			dataFile += "/";
			dataFile += alphabet[i];
			dataFile += convertInt(j);
			dataFile += "learn.txt";
		
			testData = loadTestData(dataFile);
		
		
			classifier.insertTestData(testData);
			vector<unsigned int> bestModels = classifier.getResults();
			
			confMatrixCount[i][bestModels[0]]++;
		
			cout << alphabet[bestModels[0]] << " ";
		}
		cout << endl;
		
	}
	
	cout << "CONFUSION MATRIX:" << endl;
	cout << " \t";
	for(unsigned int i = 0; i < alphabet.size(); i++) {
		cout << alphabet[i] << "\t";
	}
	cout << endl;
	for(unsigned int i = 0; i < alphabet.size(); i++) {
		cout << alphabet[i] << "\t";
		for(unsigned int j = 0; j < alphabet.size(); j++) {
			cout << confMatrixCount[i][j]/(double)NUM_TEST_DATA * 100.0 << "\t";
		}
		cout << endl;
	}
	//*/
	
  return 0;
} //end main()


vector<vector<double *> > loadTrainData(string dirWithTrainData, int numTrainSet, int P)
{
	ifstream file;
	vector<vector<double *> > trainData;
	int featureWidth, numberOfData;
	
	for(int i = 0; i < numTrainSet; i++) {
	  
		char c =  dirWithTrainData[dirWithTrainData.size() - 1];
		string pathToDataFile = dirWithTrainData + "/" + c + convertInt(i) + "learn.txt";
		
		
		///cout << pathToDataFile << endl;
		
		
		file.open(pathToDataFile.c_str());
		if(!file.is_open())
		{
			cerr << "File " << pathToDataFile << " cannot open! Data for training cannot be loaded" << endl;
			exit(1);
		}
		
		file >> featureWidth;
		
		/*if(featureWidth != P)
		{
			cerr << "File " << pathToDataFile << " has wrong feature dimensions! Dimension are different than model." << endl;
			exit(1);
		}*/
		
		file >> numberOfData;
		
		vector<double *> oneLetter;
		
		for(int n = 0; n < numberOfData; n++) {
		  
			double *featureVector = new double[featureWidth];
			
			for(int p = 0; p < featureWidth; p++) {
				file >> featureVector[p];
			}
			oneLetter.push_back(featureVector);
		}
		trainData.push_back(oneLetter);
		
		file.close();
	}
	
	return trainData;
}


void loadTrainData(double*** data, int P, int* dataNum, int numTrainSet)
{
	ifstream file;
	int featureWidth;
	
	//cout  << "feature width: " << P << ", data number: " << *dataNum << endl;
	
	string files[5] = {"N0learn.txt", "N1learn.txt", "N2learn.txt", "N3learn.txt", "N4learn.txt"};
	
	
	for(int i = 0; i < numTrainSet; i++) {
		
		file.open(files[i].c_str());
		if(!file.is_open())
		{
			cerr << "File " << files[i] << " cannot open! Data for training cannot be loaded" << endl;
			return;
		}
		
		file >> featureWidth;
		file >> dataNum[i];
		
		if(featureWidth != P)
		{
			cerr << "File " << files[i] << " has wrong feature dimensions! Dimension are different than model." << endl;
			return;
		}
		
		data[i] = new double*[dataNum[i]];
		
		for(int n = 0; n < dataNum[i]; n++) {
			data[i][n] = new double[P];
		}
		
		for(int n = 0; n < dataNum[i]; n++) {
			for(int p = 0; p < P; p++) {
				file >> data[i][n][p];
			}
		}
		
		file.close();
	}
	/*
	for(int i = 0; i < numTrainSet; i++) {
		cout << "__________" << "-" << i << "-" << "__________" << endl;
		cout << "feature width: " << P << ", data number: " << dataNum[i] << endl;
		for(int n = 0; n < dataNum[i]; n++) {
			for(int p = 0; p < P; p++) {
				cout << data[i][n][p] << "\t";
			}
			cout << endl;
		}
		
		cout << "__________" << "-" << i << "-" << "__________" << endl;
	}
	//*/
	
	return;
}

vector<double *> loadTestData(string filePath)
{
	ifstream file;
	vector<double *> testData;
	int P;
	int length;
	
	file.open(filePath.c_str());
	if(!file.is_open()) {
		cerr << "File " << filePath << " cannot open !!!" << endl;
		return testData;
	}
	
	file >> P;
	file >> length;
	
	for(int l = 0; l < length; l++) {
		double *featureVect = new double[P];
		
		for(int p = 0; p < P; p++) {
			file >> featureVect[p];
		}
		testData.push_back(featureVect);
	}
	
	file.close();
	
	return testData;
}

