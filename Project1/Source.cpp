#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <limits.h>
#include <queue>
#include <map>
#include <random>
#include <stdio.h>     
#include <stdlib.h>     
#include <time.h> 
#include <math.h>
#include "Node.h"
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
const int num_words = 4003;
struct com {
	bool operator ()(pair<pair<double, double>, double> a, pair<pair<double, double>, double> b) const {
		return a.second < b.second;
	}
};
int h_plus(int i, VectorXd train_v) {
	if (train_v(i) == 1) {
		return 1;
	}
	else {
		return -1;
	}
}
int h_minus(int i, VectorXd train_v) {
	if (train_v(i) == 1) {
		return -1;
	}
	else {
		return 1;
	}
}
int indicator_h_plus(VectorXd train_v, int i, double label) {
	if (h_plus(i, train_v) != label) {
		return 1;
	}
	else {
		return 0;
	}
}
int indicator_h_minus(VectorXd train_v, int i, double label) {
	if (h_minus(i, train_v) != label) {
		return 1;
	}
	else {
		return 0;
	}
}
pair<double, double> get_err(int n, vector<pair<VectorXd, double>> trainset, vector<double> w_i) {
	vector<double> h_t;
	double res_plus = 0;
	double res_minus = 0;
	for (int i = 0; i < trainset.size(); i++) {
		res_plus += w_i[i] * indicator_h_plus(trainset[i].first, n, trainset[i].second);
		res_minus += w_i[i] * indicator_h_minus(trainset[i].first, n, trainset[i].second);
	}

	return make_pair(res_plus, res_minus);
}
double getsign(double n) {
	if (n >= 0.0) {
		return 1;
	}
	else {
		return -1;
	}
}
double error_h_pos(int word, vector<pair<VectorXd, double>> trainset, vector<double> w_i) {
	double res = 0.0;
	for (int i = 0; i < trainset.size(); i++) {
		res += (w_i[i] * indicator_h_plus(trainset[i].first, word, trainset[i].second));
	}
	return res;
}
double error_h_neg(int word, vector<pair<VectorXd, double>> trainset, vector<double> w_i) {
	double res = 0.0;
	for (int i = 0; i < trainset.size(); i++) {
		res += (w_i[i] * indicator_h_minus(trainset[i].first, word, trainset[i].second));
	}
	return res;
}
double boosting(vector<pair<VectorXd, double>> trainset, vector<pair<VectorXd, double>> testset, int iter) {
	double err = 0.0;
	double index = 0.0;
	double initial_w = 1.0 / 4003.0;
	vector<double> w_i;
	for (int i = 0; i < trainset.size(); i++) {
		w_i.push_back(initial_w);
	}
	//initialize d

	vector<double> D_i;
	for (int i = 0; i < trainset.size(); i++) {
		D_i.push_back(1.0 / double(trainset.size()));
	}
	vector<double> h_t;
	for (int i = 0; i < iter; i++) {
		h_t.push_back(0);
	}
	vector<double> At;
	for (int i = 0; i < iter; i++) {
		At.push_back(0);
	}
	//for t = 1...T
	for (int t = 0; t < iter; t++) {
		//calculate err
		vector<pair<pair<double, double>, double>> error_v;

		////////////////////////////////////////////////////////////////
		//find the best accuracy
		double error = numeric_limits<double>::max();
		int id = 0;
		int indicator = 0;
		cout << num_words << endl;
		
		for (int i = 0; i < num_words; i++) {
			int h_indicator = 0;
			double compare = 0.0;
			double err_neg = error_h_neg(i, trainset, w_i);
			double err_pos = error_h_pos(i, trainset, w_i);
			if (err_neg > err_pos) {
				h_indicator = 1;
				compare = err_pos;
			}
			else {
				h_indicator = -1;
				compare = err_neg;
			}
			if (compare < error) {
				error = compare;
				id = i;
				indicator = h_indicator;
			}
			
		}

		

		h_t[t] = indicator;
		index = id;
		cout << "mei cuo a" << h_t[t]<<" " << id <<" "<< error << endl;
		
		error = 0.1;
		//index = error_v[0].first.first;

		//at
		if (error == 1.0||error ==0.0) {
			At[t] = 0.0;
		}
		else {
			At[t] = (1.0 / 2.0)*log((1.0 - error) / error);
		}
		cout << "at " << t << "th iteration " << At[t] << " and error " << error << endl;


		//Z_t
		double Zt = 0.0;
		vector<double> w_tp1;//wt+1
		
		for (int i = 0; i < trainset.size(); i++) {
			double w_t = 0.0;
			if (h_t[t] == 1) {
				double power = -1.0*(At[t] * trainset[i].second*h_plus(index, trainset[i].first));
				w_t = D_i[i] * exp(power);
			}
			if (h_t[t] == -1) {
				double power = -1.0*(At[t] * trainset[i].second*h_minus(index, trainset[i].first));
				w_t = D_i[i] * exp(power);
			}
			w_i[i] = w_t;
			Zt += w_t;
		}
		for (int i = 0; i < trainset.size(); i++) {
			D_i[i] = w_i[i] / Zt;
		}
		double k = 0;
		for (int i = 0; i < trainset.size(); i++) {
			k += D_i[i];
		}
		
	}
	for (int i = 0; i < testset.size(); i++) {
		VectorXd x = testset[i].first;
		double H = 0.0;
		for (int t = 0; t < iter; t++) {
			//cout << h_t[t] << endl;
			if (h_t[t] == 1) {
				H += At[t] * h_plus(index, testset[i].first);
			}
			if (h_t[t] == -1) {
				H += At[t] * h_minus(index, testset[i].first);
			}
		}
		double sign = getsign(H);
		if (sign != testset[i].second) {
			err += 1.0;
		}
	}
	err = err / testset.size();
	return err;
}

int main()
{
	VectorXd train_v(4003);
	MatrixXd train_M(450, 4004);
	vector<pair<VectorXd, double>> trainset;

	VectorXd test_v(4003);
	MatrixXd test_M(129, 4004);
	vector<pair<VectorXd, double>> testset;

	vector<string> dictionary;
	const char* train;
	const char* test;
	const char* dict;

	MatrixXd transitionM(27, 27);
	train = "pa5train.txt";
	test = "pa5test.txt";
	dict = "pa5dictionary.txt";

	ifstream fp_train(train);
	ifstream fp_test(test);
	ifstream fp_dict(dict);

	//streams
	if (!fp_train)
		cout << "test fails" << endl;
	if (!fp_test)
		cout << "test fails" << endl;
	if (!fp_dict)
		cout << "test fails" << endl;

	double num = 0;
	int row = 0;
	int col = 0;

	while (fp_train >> num) {
		double label = 0;
		if (col < 4004) {
			if (col == 4003) {
				label = num;
			}
			else {
				train_v(col) = num;
			}
			train_M(row, col) = num;
			col++;
		}
		if (col >= 4004) {
			col = 0;
			row++;
			trainset.push_back(make_pair(train_v, label));
		}
	}

	fp_train.close();
	//a1
	num = 0;
	row = 0;
	col = 0;
	while (fp_test >> num) {
		double label = 0;
		if (col < 4004) {
			if (col == 4003) {
				label = num;
			}
			else {
				test_v(col) = num;
			}
			test_M(row, col) = num;
			col++;
		}
		if (col >= 4004) {
			col = 0;
			row++;
			testset.push_back(make_pair(test_v, label));
		}
	}

	fp_test.close();

	string word;
	while (fp_dict >> word) {

		dictionary.push_back(word);
	}
	//cout << dictionary.size() << endl;
	fp_dict.close();

	cout << boosting(trainset, trainset, 4);
	int lala;
	cin >> lala;
}
