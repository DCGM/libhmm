/**
  * 
  * @author David Zouhar <xzouha@stud.fit.vutbr.cz>
  */
#include "classifier.h"

#include "stdlib.h"

#include <fstream>
#include <sstream>

Classifier::Classifier(int dimensionNums) {
	P = dimensionNums;
}

Classifier::~Classifier() {
	for(unsigned int i = 0; i < models.size(); i++) {
		delete models[i];
	}
}

void Classifier::addModels(const char *dirWithModels, string alphabet)
{
	string filepath;
	string strIndex;
	ifstream file;
	
	// search files with models and compute number of models
	for(unsigned int i = 0; i < alphabet.size(); i++)
	{
		stringstream strIndex;
		strIndex << i;
		filepath = string(dirWithModels) + "/" + alphabet[i] + ".model";
		
		models.push_back(new Model(filepath));
	}
}

void Classifier::addModel(const char *filepath)
{
	models.push_back(new Model(filepath));
}

void Classifier::insertTestData(vector<double * > timeData)
{
	data = timeData;
}

vector<unsigned int> Classifier::getResults()
{
	vector<double> results;
	for(unsigned int i = 0; i < models.size(); i++) {
		results.push_back(models[i]->evaluate(data));
	}
	
	double max;
	vector<unsigned int> sortedIndexes;
	
	for(unsigned int i = 0; i < results.size(); i++) {
		max = -INFINITY;
		unsigned int maxIndex;
		for(unsigned int j = 0; j < results.size(); j++) {
			if(results[j] > max) {
				maxIndex = j;
				max = results[j];
			}
		}
		results[maxIndex] = -INFINITY;
		sortedIndexes.push_back(maxIndex);
	}
	
	
	
	return sortedIndexes;
}