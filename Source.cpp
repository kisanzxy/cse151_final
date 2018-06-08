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
struct com {
	bool operator ()(pair<pair<double,double>,double> a, pair<pair<double,double>,double> b) const {
		return a.second < b.second;
	}
};
vector<Node*> tree;
//MatrixXd trainM(56, 819);
int c0 = 0;
int c1 = 0;

double calculate(vector<pair<VectorXd, double>> trainset, vector<pair<VectorXd, double>> testset) {
	//cout << validateset.size() << endl;
	map<double, VectorXd>w_map;
	double validerr = 0.0;
	
	map<VectorXd, double> f_map;
	VectorXd w_t = VectorXd(819);
	double m = 1.0;
	
	for (int i = 0; i < w_t.size(); i++) {
		w_t[i] = 0.0;
	}
	w_map.insert(make_pair(m, w_t));
	//cout << w_t << endl;
	for (unsigned int i = 0; i < trainset.size(); i++) {
		double y_t =trainset[i].second;
		if (y_t != -1.0 && y_t != 1.0) {
			cout << "asdgasg" << endl;
		}
		VectorXd x_t = trainset[i].first;
		
		double wxt = w_map[m].dot(x_t);
		if (y_t*wxt <= 0.0) {
			//validerr+=1.0;
			w_t = w_t + y_t * x_t;
			m += 1.0;
			w_map.insert(make_pair(m, w_t));
		}
		
	}
	for (unsigned int i = 0; i < testset.size();i++) {
		double y_t = testset[i].second;
		VectorXd x_t = testset[i].first;
		double wxt = w_map[m].dot(x_t);
		if (y_t*wxt <= 0.0) {
			validerr += 1.0;
		}
		
	}
	cout <<  validerr  << endl;
	return validerr / testset.size();
}
double calculate_voted(vector<pair<VectorXd, double>> trainset, vector<pair<VectorXd, double>> testset) {
	cout << trainset.size() << " " << testset.size() << " should be the same" << endl;
	map<double, double> c_map;
	map<double, VectorXd>w_map;
	
	double m = 1;
	//double cm = 1.0;
	c_map.insert(make_pair(m, 1.0));
	
	double validerr = 0.0;
	double w1 = 0;
	VectorXd w_t = VectorXd(819);
	for (int i = 0; i < w_t.size(); i++) {
		w_t[i] = 0.0;
	}
	w_map.insert(make_pair(m, w_t));
	//cout << w_t << endl;
	
	for (unsigned int i = 0; i < trainset.size(); i++) {
		double y_t = trainset[i].second;
		if (y_t != -1.0 && y_t != 1.0) {
			cout << "asdgasg" << endl;
		}
		
		VectorXd x_t = trainset[i].first;
		double wxt = w_map[m].dot(x_t);
		if (y_t*wxt <= 0.0) {
			VectorXd w_t = w_map[m] + y_t * x_t;
			
			m = m + 1.0;

			w_map.insert(make_pair(m,w_t));
			c_map.insert(make_pair(m, 1.0));

		}
		else {
			c_map[m] += 1.0;
		}
		
	}
	
	
	for (unsigned int i = 0; i < testset.size();i++) {
		
		double y_t = testset[i].second;
		VectorXd x_t = testset[i].first;
		double res = 0.0;
		
		for (unsigned int j = 1; j < m+1;j++) {
			double sign = -1.0;
			if (w_map[j].dot(x_t) <= 0) {
				sign = -1.0;
			}
			else {
				sign = 1.0;
			}
		//	cout << "debug1" << endl;
			double ite = c_map[j] * sign;
			//cout << "debug1" << endl;
			//double it = w_m[j].second*(signbit(w_m[j].first.dot(x_t)));
			res += ite;
		}
		double check = 1.0;//signbit(res);
		if (res <= 0) {
			check = -1.0;
		}
		else {
			check = 1.0;
		}
		if (check != y_t) {
			validerr += 1.0;
		}
	}
	cout << validerr << endl;
	return validerr / testset.size();
}
double calculate_avr_voted(vector<pair<VectorXd, double>> trainset, vector<pair<VectorXd, double>> testset) {
	cout << trainset.size() << " " << testset.size() << " should be the same" << endl;
	map<double, double> c_map;
	map<double, VectorXd>w_map;

	double m = 1;
	//double cm = 1.0;
	c_map.insert(make_pair(m, 1.0));

	double validerr = 0.0;
	double w1 = 0;
	VectorXd w_t = VectorXd(819);
	for (int i = 0; i < w_t.size(); i++) {
		w_t[i] = 0.0;
	}
	w_map.insert(make_pair(m, w_t));
	//cout << w_t << endl;

	for (unsigned int i = 0; i < trainset.size(); i++) {
		double y_t = trainset[i].second;
		if (y_t != -1.0 && y_t != 1.0) {
			cout << "asdgasg" << endl;
		}

		VectorXd x_t = trainset[i].first;
		double wxt = w_map[m].dot(x_t);
		if (y_t*wxt <= 0.0) {
			VectorXd w_t = w_map[m] + y_t * x_t;

			m = m + 1.0;

			w_map.insert(make_pair(m, w_t));
			c_map.insert(make_pair(m, 1.0));

		}
		else {
			c_map[m] += 1.0;
		}

	}


	for (unsigned int i = 0; i < testset.size(); i++) {

		double y_t = testset[i].second;
		VectorXd x_t = testset[i].first;
		VectorXd res = VectorXd(819);
		for (int i = 0; i < res.size(); i++) {
			res[i] = 0.0;
		}

		for (unsigned int j = 1; j < m + 1; j++) {

			VectorXd ite = c_map[j] * w_map[j];
			res = res + ite;
		}

		double check = res.dot(x_t);//signbit(res);
		double sign = 0.0;
		if (check <= 0.0) {
			sign = -1.0;
		}
		else {
			sign = 1.0;
		}
		if (sign != y_t) {
			validerr += 1.0;
		}
	}
	cout << validerr << endl;
	return validerr / testset.size();
}
int main()
{
//	int trainerr = 0;
	vector<pair<VectorXd, double>> trainset;
	vector<pair<VectorXd, double>> trainset_p1;
	vector<pair<VectorXd, double>> testset_p1;
	VectorXd trainV(819);
	VectorXd testV(819);
	vector<pair<VectorXd, double>> testset;
	vector<pair<VectorXd, double>> validateset;
	VectorXd validateV(819);
	MatrixXd testM(56, 819);

	string testline;
	const char* test;
	const char* train;
	const char* validate;
	const char* feature;
	MatrixXd trainM(3000, 820);
	test = "pa3test.txt";
	train = "pa3train.txt";
	validate = "pa2validation.txt";
	feature = "pa2features.txt";
	ifstream fptest(test);
	ifstream fptrain(train);
	ifstream fpvalidate(validate);
	ifstream fpfeature(feature);
	//testM = readMatrix(test);
	//streams
	if (!fptest)
		cout << "test fails" << endl;
	if (!fptrain)
		cout << "test fails" << endl;
	if (!fpvalidate)
		cout << "test fails" << endl;
	if (!fpfeature)
		cout << "test fails" << endl;

	double num = 0;
	int row = 0;
	int col = 0;
	
	while (fptrain >> num) {
		double label = 0;
		if (col < 820) {
			//trainV(col) = num;
			//trainVt(col) = num;
			if (col == 819)
				label = num;
			else {
				//trainM(row, col) = num;
				trainV(col) = num;
			}
			//trainM(row, col) = num;
			col++;
		}
		if (col >= 820) {
			col = 0;
			row++;
			trainset.push_back(make_pair(trainV, label));
			if (label == 1) {
				trainset_p1.push_back(make_pair(trainV, -1.0));
			}
			if (label == 2) {
				trainset_p1.push_back(make_pair(trainV, 1.0));
			}
			//trainsett.push_back(trainVt);
		}
	}
	cout << trainset.size()<< " "<<trainset_p1.size() << endl;
	fptrain.close();
	num = 0;
	row = 0;
	col = 0;

	while (fptest >> num) {
		double label = 0;
		if (col < 820) {
			//trainV(col) = num;
			//trainVt(col) = num;
			if (col == 819)
				label = num;
			else {
				//trainM(row, col) = num;
				testV(col) = num;
			}
			//trainM(row, col) = num;
			col++;
		}
		if (col >= 820) {
			col = 0;
			row++;
			testset.push_back(make_pair(testV, label));
			if (label == 1) {
				testset_p1.push_back(make_pair(testV, -1.0));
			}
			if (label == 2) {
				testset_p1.push_back(make_pair(testV, 1.0));
			}
			//trainsett.push_back(trainVt);
		}
	}
	//cout << testset_p1.size() << endl;
	fptest.close();/*
	num = 0;
	row = 0;
	col = 0;

	while (fpvalidate >> num) {
		double label = 0;
		if (col < 23) {
			//trainV(col) = num;
			//trainVt(col) = num;
			if (col == 22)
				label = num;
			else {
				//trainM(row, col) = num;
				validateV(col) = num;
			}
			//trainM(row, col) = num;
			col++;
		}
		if (col >= 23) {
			col = 0;
			row++;
			validateset.push_back(make_pair(validateV, label));
			//trainsett.push_back(trainVt);
		}
	}
	cout << validateset.size() << endl;
	fpvalidate.close();
	*/

	double trainerr = calculate(trainset_p1,trainset_p1);
	cout << trainerr << endl;
	
	double voted_trainerr = calculate_voted(trainset_p1, trainset_p1);
	cout << voted_trainerr << endl;
	double voted_avr_trainerr = calculate_avr_voted(trainset_p1, trainset_p1);
	cout << voted_avr_trainerr << endl;
	cout << "finished" << endl;
	int lala;
	cin>>lala;
}
