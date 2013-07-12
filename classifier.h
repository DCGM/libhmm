/**
  * 
  * @author David Zouhar <xzouha02@stud.fit.vutbr.cz>
  */
#ifndef CLASSIFIER
#define CLASSIFIER

#include "hmm.h"

#include <vector>
#include <iostream>

using namespace std;

/**
  * Classifier
  */
class Classifier {
	
	private: 
		vector<Model*> models;
		vector<double* > data;
		int P;
		//double* featureVector;
		
	public:
		Classifier(int dimensionsNum);
		~Classifier();
		
		void addModels(const char *dirWithModels, string alphabet);
		void addModel(const char* filepath);
		void insertTestData(vector<double *> timeData);
		vector<unsigned int> getResults();
};
#endif // CLASSIFIER