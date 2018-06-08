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
	bool operator ()(pair<string,double> a, pair<string, double> b) const {
		return a.second < b.second;
	}
};


int c0 = 0;
int c1 = 0;
int times = 0;
//int decimal_compare(double a, double b) {
//	int c1 = a * 10000;
//	int c2 = b * 10000;
//	if (c1 - c2 == 0) {
//		return 1;
//	}
//	else
//		return 0;
//}
bool decimal_compare(double a, double b) {
	int c1 = a * 10000;
	int c2 = b * 10000;
	if (c1 - c2 == 0) {
		return true;
	}
	else
		return false;
}
double get_a(int o, int s, MatrixXd v, int a, map<int, MatrixXd> state_map) {
	double res = 0.0;
	for (int s_prime = 1; s_prime < 81; s_prime++) {
		res += state_map[a](s, s_prime)*v(s_prime);
	}
	return res;
}

double max_a(int o, int s, MatrixXd v, map<int, MatrixXd> state_map) {
	double max = -numeric_limits<double>::max();
	for (int a = 0; a < 4; a++) {
		double term = get_a(o, s, v, a, state_map);
		if (term > max) {
			max = term;
		}
	}
	return max;
}
double get_argmax(int o, int s, MatrixXd v, map<int, MatrixXd> state_map) {
	double max = -numeric_limits<double>::max();
	double a_res = -1.0;
	for (int a = 0; a < 4; a++) {
		double term = get_a(o, s, v, a, state_map);
		if (term > max) {
			max = term;
			a_res = a;
		}
	}
	return a_res;
}
// to be rethought
pair<MatrixXd,MatrixXd> optimal_state_value(vector<double> rewards, map<int, MatrixXd> state_map) {
	MatrixXd res(9, 9);
	MatrixXd direction(9, 9);
	VectorXd v(81);
	VectorXd d(81);
	map<int, VectorXd> v_star;
	for (int i = 0; i < 81; i++) {
		v(i) = 0.0;
	}
	double difference = 1.0;
	int o = 1;
	
	while(difference != 0.0 ){
		difference = 0.0;
		VectorXd v_copy(81);
		v_copy = v;
		for (int s = 0; s < 81; s++) {
			v(s) = rewards[s] + 0.99*max_a(o,s,v_copy,state_map);//
			if (!decimal_compare(v(s), v_copy(s))) {
				difference = 1.0;
			}
		}
		o++;
	}
	for (int i = 0; i < 81; i++) {
		d(i) = get_argmax(o, i, v, state_map);
	}
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 9; j++) {
			res(j, i) = v(9 * i + j);
			direction(j, i) = d(9 * i + j);
		}
	}
	return make_pair(res,direction);
}

int main()
{
	//vector<pair<VectorXd, double>> trainset;
	MatrixXd stateM_a1(81, 81);
	stateM_a1.setZero();
	MatrixXd stateM_a2(81, 81);
	stateM_a2.setZero();
	MatrixXd stateM_a3(81, 81);
	stateM_a3.setZero();
	MatrixXd stateM_a4(81, 81);
	stateM_a4.setZero();
	vector<double> rewards_v;

	
	const char* rewards;
	const char* prob_a4;
	const char* prob_a3;
	const char* prob_a2;
	const char* prob_a1;
	MatrixXd transitionM(27, 27);
	rewards = "rewards.txt";
	prob_a4 = "prob_a4.txt";
	prob_a3 = "prob_a3.txt";
	prob_a2 = "prob_a2.txt";
	prob_a1 = "prob_a1.txt";
	ifstream fp_rewards(rewards);
	ifstream fp_prob_a4(prob_a4);
	ifstream fp_prob_a3(prob_a3);
	ifstream fp_prob_a2(prob_a2);
	ifstream fp_prob_a1(prob_a1);
	//streams
	if (!fp_rewards)
		cout << "test fails" << endl;
	if (!fp_prob_a4)
		cout << "test fails" << endl;
	if (!fp_prob_a3)
		cout << "test fails" << endl;
	if (!fp_prob_a2)
		cout << "test fails" << endl;
	if (!fp_prob_a1)
		cout << "test fails" << endl;
	double num = 0;
	int row = 0;
	int col = 0;
	
	while (fp_rewards >> num) {
		rewards_v.push_back(num);
	}
	for (int i = 0; i < rewards_v.size(); i++) {
		cout << rewards_v[i];
	}
	cout << rewards_v.size()<< endl;
	fp_rewards.close();
	//a1
	num = 0;
	row = 0;
	col = 0;
	while (fp_prob_a1 >> row && fp_prob_a1 >> col && fp_prob_a1 >> num) {
		stateM_a1(row - 1.0, col - 1.0) = num;
	}
	fp_prob_a1.close();
	//checking if parse is correct
	/*for (int i = 0; i < 81; i++) {
		double sum = 0;
		for (int j = 0; j < 81; j++) {
			sum += stateM_a1(i, j);
		}
		cout << sum << endl;
	}*/
	//a2
	num = 0;
	row = 0;
	col = 0;
	while (fp_prob_a2 >> row && fp_prob_a2 >> col && fp_prob_a2 >> num) {
		stateM_a2(row - 1.0, col - 1.0) = num;
	}
	fp_prob_a2.close();

	//a3
	num = 0;
	row = 0;
	col = 0;
	while (fp_prob_a3 >> row && fp_prob_a3 >> col && fp_prob_a3 >> num) {
		stateM_a3(row - 1.0, col - 1.0) = num;
	}
	fp_prob_a3.close();

	//a4
	num = 0;
	row = 0;
	col = 0;
	while (fp_prob_a4 >> row && fp_prob_a4 >> col && fp_prob_a4 >> num) {
		stateM_a4(row - 1.0, col - 1.0) = num;
	}
	fp_prob_a4.close();

	map<int, MatrixXd> state_map;
	state_map.insert(make_pair(0, stateM_a1));
	state_map.insert(make_pair(1, stateM_a2));
	state_map.insert(make_pair(2, stateM_a3));
	state_map.insert(make_pair(3, stateM_a4));

	cout << optimal_state_value(rewards_v, state_map).first << endl;
	cout << optimal_state_value(rewards_v, state_map).second << endl;
	
	int lala;
	cin>>lala;
}
