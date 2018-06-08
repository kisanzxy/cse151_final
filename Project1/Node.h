#pragma once
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

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
class Node
{
public:
	int fake=0;
	double label;
	double m,s;
	vector<pair<VectorXd, double>> info;
	bool visited;
	bool visited_p;
	bool visited_m;
	bool pure;
	Node * left;
	Node * right;
	Node(vector<pair<VectorXd, double>> );
	
	~Node();
	bool isPure();
};

