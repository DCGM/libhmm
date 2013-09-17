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