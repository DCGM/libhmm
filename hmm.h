/**
  * 
  * @author David Zouhar <xzouha02@stud.fit.vutbr.cz>
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
