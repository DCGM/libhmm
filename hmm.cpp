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
#include <iomanip> //setprecision - for debugging
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <sstream>
#include <cmath>

#include "hmm.h"

using namespace std;




Model::Model(int numberOfStates, int featureWidth)
{
	N = numberOfStates + 2;
	P = featureWidth;
	
	
	// alokation for model
	means = new double*[N];
	for (int i = 0; i < N; ++i) {
		means[i] = new double[P];
	}
  
	vars = new double*[N];
	for (int i = 0; i < N; ++i) {
		vars[i] = new double[P];
	}
	
	trans = new double*[N];
	for (int i = 0; i < N; ++i) {
		trans[i] = new double[N];
	}

	// means and vars for 1. and last state will be zeros, others must be loaded from file
	for (int i = 0; i < P; i++) {
		means[0][i] = means[N-1][i] = vars[0][i] = vars[N-1][i] = 0.0; 
	}
}

/**
  * Loads model from file
  */

Model::Model(const string filename)
{
	int state, i;
	FILE *fp = fopen (filename.c_str(), "r");
	
	if (!fp) {
		cerr << "model \"" << filename << "\" does not exists" << endl;
		exit(1);
	}
	
	fscanf(fp, "%d", &N); 
	N += 2; 
	fscanf(fp, "%d", &P); 
  
	// alokation for model
	means = new double*[N];
	for (int i = 0; i < N; ++i) {
		means[i] = new double[P];
	}
  
	vars = new double*[N];
	for (int i = 0; i < N; ++i) {
		vars[i] = new double[P];
	}
	
	trans = new double*[N];
	for (int i = 0; i < N; ++i) {
		trans[i] = new double[N];
	}

	// means and vars for 1. and last state will be zeros, others must be loaded from file
	for (i = 0; i < P; i++) {
		means[0][i] = means[N-1][i] = vars[0][i] = vars[N-1][i] = 0.0; 
	}
	   
	for (state = 1; state <= N-2; state++) {
		for (i = 0; i < P; i++) {
			fscanf(fp, "%lf", &means[state][i]);
		}
	}
	
	for (state = 1; state <= N-2; state++) {
		for (i=0; i<P; i++) {
			fscanf(fp, "%lf", &vars[state][i]);
		}
	}

	// transition probabilities
	for (state = 0; state < N; state++) {
		for (i = 0; i < N; i++) {
			fscanf(fp, "%lf", &trans[state][i]);
		}
	}

	fclose(fp);
}

/**
  *  dealocates memory
  */ 

Model::~Model()
{
	for (int i = 0; i < N; ++i) {
		delete [] means[i];
	}
	delete [] means;
	
	for (int i = 0; i < N; ++i) {
		delete [] vars[i];
	}
	delete [] vars;
	
	for (int i = 0; i < N; ++i) {
		delete [] trans[i];
	}
	delete [] trans;
}

/**
  * Adding logarithmic probabilities
  * - special formula
  */
double Model::logSum(const double x, const double y) 
{
	if(x == -INFINITY) {
		return y;
	}
	else if(y == -INFINITY) {
		return x;
	}
	else {
		if(x > y) {
			return x + log(1.0 + exp(y-x));
		}
		else {
			return y + log(1.0 + exp(x-y));
		}
	}
}

/**
  * Special multiplie of logarithmic probabilities
  */
double Model::logProduct(const double x, const double y) 
{
	if((x == -INFINITY) || (y == -INFINITY)) {
		return -INFINITY;
	}
	else {
		return x + y;
	}
}

/**
  * Initialize model with one feature vector
  */
void Model:: init_model(double **featureVect,  const int T)
{
	double step = ((double)(N-3)/(double)(T-1));
	int index = 0;
	double sumaForIndex = 0.0;
	double suma[P];
	double number = 0.0;
	
	
	// compute initialize means
	for(int p = 0; p < P; p++)
		suma[p] = 0.0;
	
	for(int t = 0; t < T; t++) {
		if(index == round(sumaForIndex)) {
			for(int p = 0; p < P; p++) {
				suma[p] += featureVect[t][p];
			}
			number++;
		}
		else {
			for(int p = 0; p < P; p++) {;
				means[index+1][p] = suma[p] / number;
				suma[p] = 0.0;
			}
			number = 0.0;
			
			for(int p = 0; p < P; p++) {
				suma[p] += featureVect[t][p];
			}
			number++;
		}
		index = round(sumaForIndex);
		sumaForIndex += step;
	}
	//last state
	for(int p = 0; p < P; p++) {
		means[index+1][p] = suma[p] / number;
	}
	
	//compute initialize variables
	index = 0;
	sumaForIndex = 0.0;
	number = 0.0;
	for(int p = 0; p < P; p++)
		suma[p] = 0.0;
	
	for(int t = 0; t < T; t++) {
		if(index == round(sumaForIndex)) {
			for(int p = 0; p < P; p++) {
				suma[p] += pow(featureVect[t][p] - means[index+1][p], 2.0);
			}
			number++;
		}
		else {
			for(int p = 0; p < P; p++) {
				vars[index+1][p] = sqrt(suma[p] / (number-1));
				if(suma[p] == 0.0) {
					vars[index+1][p] = 1.0;  // variable cannot be zero -> gauss returns nan
				}
				suma[p] = 0.0;
			}
			number = 0.0;
			
			for(int p = 0; p < P; p++) {
				suma[p] += pow(featureVect[t][p] - means[index+2][p], 2.0);
			}
			number++;
		}
		index = round(sumaForIndex);
		sumaForIndex += step;
	}
	
	//last state
	for(int p = 0; p < P; p++) {
		vars[index+1][p] = sqrt(suma[p] / (number-1));
	}
	
	// initialize transition probabilities
	trans[0][1] = 1.0;
	for(int n = 1; n < N-1; n++) {
		trans[n][n] = 0.5;
		trans[n][n+1] = 0.5;
	}
}

void Model::init_model(vector<double *> featureVect)
{
	double step = ((double)(N-3)/(double)(featureVect.size() - 1));
	int index = 0;
	double sumaForIndex = 0.0;
	double suma[P];
	double number = 0.0;
	
	
	// compute initialize means
	for(int p = 0; p < P; p++)
		suma[p] = 0.0;
	
	for(unsigned int t = 0; t < featureVect.size(); t++) {
		if(index == round(sumaForIndex)) {
			for(int p = 0; p < P; p++) {
				suma[p] += featureVect[t][p];
			}
			number++;
		}
		else {
			for(int p = 0; p < P; p++) {
				means[index+1][p] = suma[p] / number;
				suma[p] = 0.0;
			}
			number = 0.0;
			
			for(int p = 0; p < P; p++) {
				suma[p] += featureVect[t][p];
			}
			number++;
		}
		index = round(sumaForIndex);
		sumaForIndex += step;
	}
	//last state
	for(int p = 0; p < P; p++) {
		means[index+1][p] = suma[p] / number;
	}
	
	//compute initialize variables
	index = 0;
	sumaForIndex = 0.0;
	number = 0.0;
	for(int p = 0; p < P; p++)
		suma[p] = 0.0;
	
	for(unsigned int t = 0; t < featureVect.size(); t++) {
		if(index == round(sumaForIndex)) {
			for(int p = 0; p < P; p++) {
				suma[p] += pow(featureVect[t][p] - means[index+1][p], 2.0);
			}
			number++;
		}
		else {
			for(int p = 0; p < P; p++) {
				vars[index+1][p] = sqrt(suma[p] / (number-1));
				if((suma[p] == 0.0) || (number <= 1)) {
					vars[index+1][p] = 1.0;  // variable cannot be zero -> gauss returns nan
				}
				suma[p] = 0.0;
			}
			number = 0.0;
			
			for(int p = 0; p < P; p++) {
				suma[p] += pow(featureVect[t][p] - means[index+2][p], 2.0);
			}
			number++;
		}
		index = round(sumaForIndex);
		sumaForIndex += step;
	}
	
	//last state
	for(int p = 0; p < P; p++) {
		if((number <= 1) || (suma[p] == 0.0)) {
			vars[index+1][p] = 1;
		}
		else {
			vars[index+1][p] = sqrt(suma[p] / (number-1));
		}
	}
	
	// initialize transition probabilities
	for(int n = 0; n < N; n++) {
		for(int nn = 0; nn < N; nn++) {
			trans[n][nn] = 0.0;
		}
	}
	trans[0][1] = 1.0;
	for(int n = 1; n < N-1; n++) {
		trans[n][n] = 0.5;
		trans[n][n+1] = 0.5;
	}
}

/**
  * computes viterbi log likelihoods
  *
  * @return viterbi likelihood for feature vector
  * @param fea     feature vector
  */
double Model::evaluate(vector<double *> fea) 
{
	double *this_round = new double[N];
	double *next_round = new double[N];
	double a, token;
	
	/* initialize the algorithm - inserting token in first all 
	    entry states - means the 1st state */
	for (int i = 0; i < N; i++) {
		this_round[i] = next_round[i] = -INFINITY;
	}
	
	this_round[0] = 0.0;
	
	double likelyhood;
	
	/* run Viterbi - propagate tokens - for all times.  */
	
	for (unsigned int t = 0; t < fea.size(); t++) { // projde vsechy casy feature vektoru
		
		for (int i = 0; i < N; i++) { // projde vsecky stavy v soucasnem case
			
			for (int j = 0; j < N; j++) { // projde vsecky stavy v nasledujicim case
				
				if ( (trans[i][j] != 0.0f) ) { // prechodova pravdepodobnost != 0
					
					likelyhood = 
						// viterbiho pravdepodobnost z predchoziho stavu
						this_round[i]
						+ (
						// vysilaci pravdepodobnost
						LogGauss(fea[t], means[j], vars[j], P)
						// prechodova pravdepodobnost
						+ log(trans[i][j])
						);          
					// prirazeni maximalniho likelyhoodu
					next_round[j] = fmaxf(likelyhood, next_round[j]);
				} // endif (trans[i][j] != 0.0f)
			} // endfor(j: 0-N)
		} // endfor(i: 0-N)
		
		/* token passing for one time 't' is ready, next becomes this_round */
		
		for (int i = 0; i < N; i++) {
			  
			this_round[i] = next_round[i]; 
			next_round[i] = -INFINITY;
		}
	} /* end for all times */
	
	/* finish the Viterbi - see what is connected to exit state and add just 
	    the transition proba */
	next_round[N-1] = -INFINITY;
	
	for (int from = 1; from < N-1; from++) {
		token = this_round[from]; 
		if (token == -INFINITY) {
			continue; /* no token */
		}
		
		a = trans[from][N-1];
		
		if (a == 0.0) {
			continue;  /* no transition */ 
		}
		
		token += log(a);  /* only trans proba, no emission */
		
		if (token > next_round[N-1]) {
			next_round[N-1] = token;
		}
		
	}
	
	return next_round[N-1];
}

/**
  *  Compute alpha probabilities in logarithmic scale
  */
double Model::logForward(vector<double *> featureVect, double **alpha)
{
	double sum;
	double baumWelchProba = 0.0;
	
	// initialize alpha matrix to zeros
	for(int i = 0; i < N; i++) {
		for(unsigned int j = 0; j < featureVect.size(); j++) {
			alpha[i][j] = -INFINITY;
		}
	}
	
	// first step
	for(int j = 1; j < N-1; j++) {
		alpha[j][0] = logProduct(log(trans[0][j]), LogGauss(featureVect[0], means[j], vars[j], P));
	}
	
	// go through time of features
	for(unsigned int t = 1; t < featureVect.size(); t++) {
		// compute alpha for all states int time t (number of states is N - first and last state -> states for help)
		for(int j = 1; j < (N - 1); j++) {
			sum = -INFINITY;
			
			//sum all probabilities from states, from we can reach observed state
			for(int i = 1; i < (N - 1); i++) {
				sum = logSum(sum, logProduct(alpha[i][t-1], log(trans[i][j])) );
			}
			
			// add (multiplie with) emit probability this state for features in time t
			alpha[j][t] = logProduct(sum, LogGauss(featureVect[t], means[j], vars[j], P));
		}
	}
	
	// compute Baum-Welch probability
	sum = -INFINITY;
	for(int i = 1; i < N-1; i++) {
		sum = logSum(sum, logProduct(alpha[i][featureVect.size()-1], log(trans[i][N-1])) );
	}
	baumWelchProba = sum;
	
	/*
	///pomocny vypis
	cout << "----------------log alpha----------------------------------" << endl << "cas:\t  ";
	for(int j = 0; j < T; j++) {
		cout << j << "\t  ";
	}
	cout << endl;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < T; j++) {
			//cout << scientific  << setprecision(3) << "\t" << exp(alpha[i][j]);
			if(alpha[i][j] == -INFINITY) {
				cout << scientific << "\t" << alpha[i][j] << "\t";
			}
			else {
				cout << scientific  << setprecision(3) << "\t" << alpha[i][j];
			}
		}
		cout << endl;
	}
	cout << endl;
	
	
	//*/
	
	//cout << "log alpha baum welch: " << baumWelchProba << endl;
	return baumWelchProba;
}

/**
  *  Compute beta probabilities in logarithmic scale
  */
double Model::logBackward(vector<double *> featureVect, double **beta) 
{
	double sum;
	double baumWelchProba;
	
	// initialize beta matrix to zeros
	for(int i = 0; i < N; i++) {
		for(unsigned int j = 0; j < featureVect.size(); j++) {
			beta[i][j] = -INFINITY;
		}
	}
	
	// first step
	for(int i = 1; i < N-1; i++){
		beta[i][featureVect.size()-1] = log(trans[i][N-1]);
	}
	// go back through time of features
	for(int t = featureVect.size()-2; t >= 0; t--) {
		// compute beta for all states int time t (number of states is N - first and last state -> states for help)
		for(int j = 1; j < (N - 1); j++) {
			sum = -INFINITY;
			
			//sum all probabilities from states, from we can reach observed state
			for(int i = 1; i < (N - 1); i++) {
				sum = logSum(sum, logProduct(log(trans[j][i]), 
							     logProduct(LogGauss(featureVect[t+1], means[i], vars[i], P),
									beta[i][t+1]) 
									) 
					    );
			}
			
			// get probability to state
			beta[j][t] = sum;
		}
	}
	
	sum = -INFINITY;
	for(int i = 1; i < (N - 1); i++) {
		sum = logSum(sum, logProduct(log(trans[0][i]), logProduct(LogGauss(featureVect[0], means[i], vars[i], P), beta[i][0])) );
	}
	baumWelchProba = sum;
	
	/*
	cout << "-----------------beta----------------------------------" << endl << "cas:\t  ";
	for(int j = 0; j < T; j++) {
		cout << j << "\t\t  ";
	}
	cout << endl;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < T; j++) {
			cout << setw(3) << setprecision( 4 ) << "\t" << beta[i][j] ;
		}
		cout << endl;
	}
	cout << endl;
	//*/
	
	// return baum-welch probability
	//cout << "beta baum welch: " << baumWelchProba << endl;
	return baumWelchProba;
}

/**
  * Baum-Welch algorithm for training - FOR ONE FEATURE VECTOR
  * - in logarithmic scale
  */
void Model::logBaum_welch(vector<double *> featureVect, double** newMeans, double** newVars, double ** newTrans)
{
	double baumWelchProba;
	
	double** logAlpha = new double *[N];
	double** logBeta = new double *[N];
	for(int i = 0; i < N; i++) {
		logAlpha[i] = new double[featureVect.size()];
		logBeta[i] = new double[featureVect.size()];
	}
	
	double** logL = new double *[N];
	for(int i = 0; i < N; i++) {
		logL[i] = new double[featureVect.size()];
	}
	
	double** logGamma = new double *[N];
	for(int i = 0; i < N; i++) {
		logGamma[i] = new double[featureVect.size()];
	}
	
	baumWelchProba = logForward(featureVect, logAlpha);
	logBackward(featureVect, logBeta);
	
	
	//compute Ljt - state occupation function
	for(int j = 0; j < N; j++) {
		for(unsigned int t = 0; t < featureVect.size(); t++) {
			logL[j][t] = logProduct(logProduct(logAlpha[j][t], logBeta[j][t]), -baumWelchProba);
		}
	}
	
	
	//compute trans probabilities
	double numerator;
	double denominator;
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			numerator = -INFINITY;
			denominator = -INFINITY;
			for(unsigned int t = 0; t < featureVect.size(); t++) {
				denominator = logSum(denominator, logL[j][t]);
			}
			
			for(unsigned int t = 0; t < featureVect.size()-1; t++) {
				numerator = logSum(numerator, 
						   logProduct(logAlpha[i][t], 
							      logProduct(log(trans[i][j]), 
									 logProduct(LogGauss(featureVect[t+1], means[j], vars[j], P), logBeta[j][t+1])
									 )
							      )
						    );
			}
			newTrans[i][j] = exp(logProduct(numerator, -denominator));
		}
	}
	
	// transition probability from last state
	for(int i = 0; i < N; i++) {
		numerator = -INFINITY;
		denominator = -INFINITY;
		
		numerator = logProduct(logAlpha[i][featureVect.size()-1], logBeta[i][featureVect.size()-1]);
		
		for(unsigned int t = 0; t < featureVect.size(); t++) {
			denominator = logSum(denominator, logL[i][t]);
		}
		newTrans[i][N-1] = exp(logProduct(numerator, -denominator));
	}
	
	//normalize trans proba - sum of proba from one state = 1.0
	double normalizer;
	
	for(int i = 0; i < N; i++) {
		if(i == 0) {
		      newTrans[i][1] = 1.0;
		      continue;
		}
		normalizer = 0.0;
		for(int j = 0; j < N; j++) {
			normalizer += newTrans[i][j];
		}
		for(int j = 0; j < N; j++) {
			if(normalizer == 0.0) {
				newTrans[i][j] = 0.0;
			}
			else {
				newTrans[i][j] /= normalizer;
			}
		}
	}
	
	
	//means and vars
	double sumMeans, sumVars;
	for(int j = 0; j < N; j++) {
		for(int p = 0; p < P; p++) {
			numerator = -INFINITY;
			denominator = -INFINITY;
			for(unsigned int t = 0; t < featureVect.size(); t++) {
				denominator = logSum(denominator, logL[j][t]);
			}
			
			sumMeans = 0.0;
			sumVars = 0.0;
			for(unsigned int t = 0; t < featureVect.size(); t++) {
				sumMeans += exp(logProduct(logL[j][t], -denominator)) * featureVect[t][p];
				sumVars += exp(logProduct(logL[j][t], -denominator)) * pow(featureVect[t][p] - means[j][p], 2.0);
			}
			newMeans[j][p] = sumMeans;
			newVars[j][p] = sqrt(sumVars);
			if(newVars[j][p] < 0.01) {
				newVars[j][p] = 0.1;
			}
		}
	}
}

/**
  * Baum-Welch algorithm for training - FOR MULTIPLIE FEATURE VECTOR
  * - in logarithmic scale
  */
double Model::logMultBaum_welchStep(vector<vector<double *> > featureVect, double** newMeans, double** newVars, double ** newTrans)
{
	double baumWelchProba[featureVect.size()];
	double totalLogLikelihood = 0.0;
	
	double*** logAlpha = new double **[featureVect.size()];
	double*** logBeta = new double **[featureVect.size()];
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		logAlpha[f] = new double *[N];
		logBeta[f] = new double *[N];
		for(int n = 0; n < N; n++) {
			logAlpha[f][n] = new double[featureVect[f].size()];
			logBeta[f][n] = new double[featureVect[f].size()];
		}
	}
	
	double*** logL = new double **[featureVect.size()];
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		logL[f] = new double *[N];
		for(int n = 0; n < N; n++) {
			logL[f][n] = new double[featureVect[f].size()];
		}
	}
	
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		baumWelchProba[f] = logForward(featureVect[f], logAlpha[f]);
		logBackward(featureVect[f], logBeta[f]);
		
		totalLogLikelihood = logProduct(totalLogLikelihood, baumWelchProba[f]);
		cout << f << "-" << baumWelchProba[f] << "\t";
		//cout << f << ". baumWelchProba: " << baumWelchProba[f] << "\tfeatureVect[f].size(): " << featureVect[f].size() << endl;
	}
	cout << endl;
	
	// Compute state occupation function -> Ljt
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		for(int j = 0; j < N; j++) {
			for(unsigned int t = 0; t < featureVect[f].size(); t++) {
				logL[f][j][t] = logProduct(logProduct(logAlpha[f][j][t], logBeta[f][j][t]), -baumWelchProba[f]);
			}
		}
	}
	
	// Compute transition probability
	double numerator;
	double denominator;
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			numerator = -INFINITY;
			denominator = -INFINITY;
			
			for(unsigned int f = 0; f < featureVect.size(); f++) {
				for(unsigned int t = 0; t < featureVect[f].size(); t++) {
					denominator = logSum(denominator, logL[f][j][t]);
				}
				
				for(unsigned int t = 0; t < featureVect[f].size()-1; t++) {
					numerator = logSum(numerator, 
							   logProduct(logAlpha[f][i][t], 
								      logProduct(log(trans[i][j]), 
										 logProduct(LogGauss(featureVect[f][t+1], means[j], vars[j], P), logBeta[f][j][t+1])
										 )
								      )
							    );
				}
			}
			
			newTrans[i][j] = exp(logProduct(numerator, -denominator));
		}
	}
	
	// transition probability from last state
	for(int i = 0; i < N; i++) {
		numerator = -INFINITY;
		denominator = -INFINITY;
		
		for(unsigned int f = 0; f < featureVect.size(); f++) {
			numerator = logSum(numerator, logProduct(logAlpha[f][i][featureVect[f].size()-1], logBeta[f][i][featureVect[f].size()-1]) );
			
			for(unsigned int t = 0; t < featureVect[f].size(); t++) {
				denominator = logSum(denominator, logL[f][i][t]);
			}
		}
		newTrans[i][N-1] = exp(logProduct(numerator, -denominator));
	}
	
	//normalize trans proba
	double normalizer;
	
	for(int i = 0; i < N; i++) {
		if(i == 0) {
		      newTrans[i][1] = 1.0;
		      continue;
		}
		normalizer = 0.0;
		for(int j = 0; j < N; j++) {
			normalizer += newTrans[i][j];
		}
		for(int j = 0; j < N; j++) {
			if(normalizer == 0.0) {
				newTrans[i][j] = 0.0;
			}
			else {
				newTrans[i][j] /= normalizer;
			}
		}
	}
	
	
	// Compute new means and vars
	double sumMeans, sumVars;
	
	for(int j = 1; j < (N - 1); j++) {
		for(int p = 0; p < P; p++) {
			denominator = -INFINITY;
			
			sumMeans = 0.0;
			sumVars = 0.0;
			
			for(unsigned int f = 0; f < featureVect.size(); f++) {
				for(unsigned int t = 0; t < featureVect[f].size(); t++) {
					denominator = logSum(denominator, logL[f][j][t]);
				}
			}
			for(unsigned int f = 0; f < featureVect.size(); f++) {
				for(unsigned int t = 0; t < featureVect[f].size(); t++) {
					sumMeans += exp(logProduct(logL[f][j][t], -denominator)) * featureVect[f][t][p];
					sumVars += exp(logProduct(logL[f][j][t], -denominator)) * pow(featureVect[f][t][p] - means[j][p], 2.0);
				}
			}
			
			newMeans[j][p] = sumMeans;
			newVars[j][p] = sqrt(sumVars);
			if(newVars[j][p] < 0.01) {
				newVars[j][p] = 0.1;
			}
		}
	}
	
	return totalLogLikelihood;
}


void Model::train(vector<vector<double *> > featureVect)
{

	double** newMeans;
	double** newVars;
	double** newTrans;
	newMeans = new double*[N];
	newVars = new double*[N];
	newTrans = new double*[N];
	for (int i = 0; i < N; ++i) {
		newMeans[i] = new double[P];
		newVars[i] = new double[P];
		newTrans[i] = new double[N];
	}
	
	// allocate matrix of alpha values - not use - only for forward function
	double*** alpha = new double**[featureVect.size()];
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		alpha[f] = new double*[N];
		for(int i = 0; i < N; i++) {
			alpha[f][i] = new double[featureVect[f].size()];
		}
	}
	
	double*** beta = new double**[featureVect.size()];
	for(unsigned int f = 0; f < featureVect.size(); f++) {
		beta[f] = new double*[N];
		for(int i = 0; i < N; i++) {
			beta[f][i] = new double[featureVect[f].size()];
		}
	}
	
	/**
	cout << endl << "old trans -------------------" << endl;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			cout << fixed << setprecision(4) << trans[i][j] << "\t"; 
		}
		cout << endl;
	}
	cout << endl << "old means -------------------" << endl;
	for(int p = 0; p < P; p++) {
		for(int n = 0; n < N; n++) {
			cout << fixed << setprecision(4) <<  means[n][p] << "\t"; 
		}
		cout << endl;
	}
	cout << endl << "old vars -------------------" << endl;
	for(int p = 0; p < P; p++) {
		for(int n = 0; n < N; n++) {
			cout << fixed << setprecision(4) << vars[n][p] << "\t"; 
		}
		cout << endl;
	}
	//*/
	
	double oldLogLikelihood = 0.0;
	double newLogLikelihood = 0.0;
	bool isTrained = false;
	double epsilon = featureVect.size() / 10.0;
	///TODO: remove i -> is for cout in testing
	int i = 0.0;
	
	while(!isTrained) {
		i++;
		
		newLogLikelihood = logMultBaum_welchStep(featureVect, newMeans, newVars, newTrans);
		cout << i << ". cycle: " << newLogLikelihood << "\tepsilon: " << abs(newLogLikelihood - oldLogLikelihood) << endl;
		
		if(i > 100) {
			break;
		}
		
		
		if(abs(newLogLikelihood - oldLogLikelihood) < epsilon) {
			isTrained = true;
		}
		oldLogLikelihood = newLogLikelihood;
		
		// copy new values
		for(int n = 0; n < N; n++) {
			for(int nn = 0; nn < N; nn++) {
				trans[n][nn] = newTrans[n][nn];
			}
			for(int p = 0; p < P; p++) {
				means[n][p] = newMeans[n][p];
			}
			for(int p = 0; p < P; p++) {
				vars[n][p] = newVars[n][p];
			}
		}
	}
	
	/**
	
	cout << endl << "new trans -------------------" << endl;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			cout << fixed << setprecision(4) << newTrans[i][j] << "\t"; 
		}
		cout << endl;
	}
	
	cout << endl << "new means -------------------" << endl;
	for(int p = 0; p < P; p++) {
		for(int n = 0; n < N; n++) {
			cout << fixed << setprecision(4) <<  newMeans[n][p] << "\t"; 
		}
		cout << endl;
	}
	
	cout << endl << "new vars -------------------" << endl;
	for(int p = 0; p < P; p++) {
		for(int n = 0; n < N; n++) {
			cout << fixed << setprecision(4) << newVars[n][p] << "\t"; 
		}
		cout << endl;
	}
	//*/
	
}

/**
  * method for printing of model to OutputStream 
  */

std::ostream& Model::toOstream (std::ostream &s) 
{
	std::stringstream ss;
	int state, i;
	
	ss << N-2 << endl;
	ss << P << endl;
	ss << scientific;
	
	//ss << "Means:" << std::endl;
	//ss << "--------------------------------" << std::endl;
	for (state = 1; state <= N-2; state++) {
		for (i = 0; i < P; i++) {
			ss << means[state][i] << " ";
		}
		ss << std::endl;
	}
	//ss << std::endl;
	
	//ss << "Vars:" << std::endl;
	//ss << "--------------------------------" << std::endl;
	for (state = 1; state <= N-2; state++) {
		for (i=0; i<P; i++) {
			ss << vars[state][i] << " ";
		}
		ss << std::endl;
	}
	//ss << std::endl;
	
	// transition probas 
	//ss << "Trans:" << std::endl;
	//ss << "--------------------------------" << std::endl;
	//ss.setf(ios::fixed, ios::floatfield);
	for (state = 0; state < N; state++) {
		for (i = 0; i < N; i++) {
			ss << trans[state][i] << "\t";
		}
		ss << std::endl;
	}
	
	return s << ss.str();
}

///////////////////////////////////////////////////////////////////////////////

/**
  * Function compute probability of N-dimensional observation vector
  *	compute only gaussian function with diagonal covariance matrix
  *		-> it can be compute as product of multiply one-dimensional matrix for every dimension
  */
double gauss(double  *observation, double *means, double *variances, int size)
{
	
	double firstPartNormal;
	double secondPartNormal;
	double result = 1.0;
	
	for(int i = 0; i < size; i++) {
		secondPartNormal = exp(-(pow((observation[i] - means[i]), 2) / (2.0 * pow(variances[i], 2))));
		firstPartNormal = secondPartNormal/(sqrt(2.0 * M_PI));
		result *= firstPartNormal/variances[i];
	}
	if(result == 0.0) {
		return 1e-10;
	}
	return result;

}

/**
  * Function compute logarithmic probability of N-dimensional observation vector
  *	compute only gaussian function with diagonal covariance matrix
  */
double LogGauss(double  *observation, double *means, double *variances, int size)
{
	double cov_det = 1;
	double logPdf = 0;
	int i;
	
	for(i = 0; i < size; i++) {
		cov_det *= variances[i];
		logPdf += SQR(observation[i] - means[i]) / variances[i];
	}
	logPdf += log(cov_det * pow(2 * M_PI, size));
	
	return -0.5 * logPdf;
}