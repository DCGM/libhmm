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


#ifndef HMM_H
#define HMM_H

#include <vector>
#include <string>
#include <iostream>

#define SQR(a) ((a)*(a))

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef INFINITY
#define INFINITY 1.0/0.0
#endif

using namespace std;

double gauss(double  *observation, double *means, double *variances, int size);
double LogGauss(double  *observation, double *means, double *variances, int size);

/**
  * Hidden Markov Model
  */

class Model {
	
	private: 
		int N, P;	// emiting states, feature width
		double** means;
		double** vars;
		double** trans;	// transition matrix
		
		//Methods
		double logSum(const double x, const double y);
		double logProduct(const double x, const double y);
		
		void logBaum_welch(vector<double *> featureVect, double** newMeans, double** newVars, double ** newTrans);
		
		double logForward(vector<double *> featureVect, double **alpha);
		double logBackward(vector<double *> featureVect, double **beta);
		double logMultBaum_welchStep(vector<vector<double *> > featureVect, double** newMeans, double** newVars, double ** newTrans);
		
	public:
		Model(int numberOfStates, int featureWidth);
		Model(string filename);
		~Model();
		
		void init_model(double **featureVect,  const int featureLen);
		void train(double ***featureVect, const int* featureLen, const int featureNum);
		
		double evaluate(vector<double *> fea);
		void init_model(vector<double *>featureVect);
		void train(vector<vector <double *> > featureVect);
		
		std::ostream& toOstream (std::ostream &s);
};

inline std::ostream& operator << (std::ostream &s, Model& m) { return m.toOstream(s); }

class Exception
{
	private:
		string reason;
	public:
		Exception(string d) {
			reason = d;
		}
		void set(string d) { reason = d; }
		string get() { return reason; }
};

#endif // HMM_H
