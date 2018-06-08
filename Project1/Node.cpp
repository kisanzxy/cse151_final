#include "Node.h"



Node::Node(vector<pair<VectorXd, double>> trainset)
{
	label = 0;
	m = 0;
	s = 0;
	visited_p = false;
	visited = false;
	info = trainset;
	left = NULL;
	right = NULL;
}


Node::~Node()
{
	delete(this);
}

bool Node::isPure(){
	//visited = false;
	//bool pure = true;
	double label = info[0].second;
	for (int i = 0; i < info.size(); i++) {
		if (info[i].second != label)
			return false;
	}
	return true;
}