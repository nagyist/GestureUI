#include "usage.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GString.h"
#include <iostream>
#include <cstring>

using namespace GClasses;
using std::string;
using std::vector;
using std::cout;

UsageNode::UsageNode(const char* templ, const char* descr)
{
#ifdef DEBUG_HELPERS
	p0 = NULL; p1 = NULL; p2 = NULL; p3 = NULL;
#endif
	string s;
	while(true)
	{
		while(*templ == ' ')
			templ++;
		if(*templ == '\0')
			break;
		size_t i;
		for(i = 0; templ[i] != ' ' && templ[i] != '\0'; i++)
		{
		}
		s.assign(templ, i);
		templ += i;
		m_parts.push_back(s);
	}
	m_description = descr;
#ifdef DEBUG_HELPERS
	if(m_parts.size() > 0) p0 = m_parts[0].c_str();
	if(m_parts.size() > 1) p1 = m_parts[1].c_str();
	if(m_parts.size() > 2) p2 = m_parts[2].c_str();
	if(m_parts.size() > 3) p3 = m_parts[3].c_str();
#endif
}

UsageNode::~UsageNode()
{
	for(vector<UsageNode*>::iterator it = m_choices.begin(); it != m_choices.end(); it++)
		delete(*it);
}

UsageNode* UsageNode::add(const char* templ, const char* descr)
{
	UsageNode* pNode = new UsageNode(templ, descr);
	m_choices.push_back(pNode);
	return pNode;
}

UsageNode* UsageNode::choice(const char* name)
{
	for(vector<UsageNode*>::iterator it = m_choices.begin(); it != m_choices.end(); it++)
	{
		if(std::strcmp((*it)->tok(), name) == 0)
			return *it;
	}
	return NULL;
}

void UsageNode::sig(string* pS)
{
	vector<string>::iterator it = m_parts.begin();
	pS->append(*it);
	for(it++; it != m_parts.end(); it++)
	{
		pS->append(" ");
		pS->append(*it);
	}
}

void UsageNode::print(int depth, int tabSize, int maxWidth, bool descriptions)
{
	if(!descriptions && m_parts.size() == 1 && m_choices.size() == 0 && m_parts[0].length() > 0 && m_parts[0][0] != '-')
		return;

	// Print the token and args
	for(int i = 0; i < depth; i++)
		for(int j = 0; j < tabSize; j++)
			cout << " ";
	vector<string>::iterator it = m_parts.begin();
	cout << (*it);
	for(it++; it != m_parts.end(); it++)
	{
		cout << " ";
		cout << (*it);
	}
	cout << "\n";

	// Print the description
	if(descriptions)
	{
		GStringChopper sc(m_description.c_str(), 10, maxWidth - depth * tabSize, true);
		while(true)
		{
			const char* szLine = sc.next();
			if(!szLine)
				break;
			for(int i = 0; i <= depth; i++)
				for(int j = 0; j < tabSize; j++)
					cout << " ";
			cout << szLine;
			cout << "\n";
		}
	}

	// Print the children
	for(vector<UsageNode*>::iterator it = m_choices.begin(); it != m_choices.end(); it++)
		(*it)->print(depth + 1, tabSize, maxWidth, descriptions);
}

UsageNode* makeMasterUsageTree()
{
	UsageNode* pRoot = new UsageNode("[app]", "Welcome to the Waffles Wizard. This wizard will help you build a Waffles command.");
	pRoot->choices().push_back(makeLearnUsageTree());
	pRoot->choices().push_back(makePlotUsageTree());
	pRoot->choices().push_back(makeTransformUsageTree());
	pRoot->choices().push_back(makeGenerateUsageTree());
	return pRoot;
};

UsageNode* makeLearnUsageTree()
{
	UsageNode* pRoot = new UsageNode("waffles_learn [command]", "Supervised learning, transduction, cross-validation, etc.");
	UsageNode* pTrain = pRoot->add("train <options> [dataset] <data_opts> [algorithm]", "Trains a supervised learning algorithm. The trained model-file is printed to stdout. (Typically, you will want to pipe this to a file.)");
	{
		UsageNode* pOpts = pTrain->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pTrain->add("[dataset]", "The filename of a dataset in \".arff\" format. The last attribute (or the last n attributes if there are multiple label dims) is assumed to be the label. (If you wish to predict values for another attribute, you should swap columns to move it to the end.)");
		UsageNode* pDO = pTrain->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pPredict = pRoot->add("predict <options> [model-file] [dataset] <data_opts>", "Predict labels for all of the patterns in [dataset]. Results are printed in the form of a \".arff\" file (including both features and predictions) to stdout.");
	{
		UsageNode* pOpts = pPredict->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pPredict->add("[model-file]", "The filename of a trained model. (This is the file to which you saved the output when you trained a supervised learning algorithm.)");
		pPredict->add("[dataset]", "The filename of a dataset in \".arff\" format. (There should already be placeholder labels in this dataset. The placeholder labels will be replaced in the output by the labels that the model predicts.)");
		UsageNode* pDO = pPredict->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pPOP = pRoot->add("predictonepattern <options> [model-file] [data-set] <data_opts> [pattern]", "Predict labels for a single pattern and print the prediction to stdout. Confidence levels are also reported.");
	{
		UsageNode* pOpts = pPOP->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pPOP->add("[model-file]", "The filename of a trained model. (This is the file to which you saved the output when you trained a supervised learning algorithm.)");
		pPOP->add("[data-set]", "The filename of a \".arff\" file from which to obtain meta-data. This can be the training set or the test set. It doesn't matter which, because the data is ignored. Only the meta-data, such as the string names of attribute values, are obtained from this dataset.");
		pPOP->add("[pattern]", "A list of feature values separated by spaces. (A \"?\" may be used for unknown feature values if the model supports using unknown feature values.)");
		UsageNode* pDO = pPOP->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pTest = pRoot->add("test <options> [model-file] [dataset] <data_opts>", "Test a trained model using some test data. Results are printed to stdout for each dimension in the label vector. Predictive accuracy is reported for nominal label dimensions, and mean-squared-error is reported for continuous label dimensions.");
	{
		UsageNode* pOpts = pTest->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pTest->add("[model-file]", "The filename of a trained model. (This is the file to which you saved the output when you trained a supervised learning algorithm.)");
		pTest->add("[dataset]", "The filename of a test dataset in \".arff\" format. (This dataset must have the same number of columns as the dataset with which the model was trained.)");
		UsageNode* pDO = pTest->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pTransduce = pRoot->add("transduce <options> [labeled-set] <data_opts1> [unlabeled-set] <data_opts2> [algorithm]", "Predict labels for [unlabeled-set] based on the examples in [labeled-set]. For most algorithms, this is the same as training on [labeled-set] and then predicting labels for [unlabeled-set]. Some algorithms, however, have no models. These can transduce, even though they cannot be trained. Results are printed to stdout as a \".arff\" file.");
	{
		UsageNode* pOpts = pTransduce->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pTransduce->add("[labeled-set]", "The filename of a dataset in \".arff\" format. The labels in this dataset are used to infer labels for the unlabeled set. The last attribute (or the last n attributes if there are multiple label dims) is assumed to be the label. (If you wish to predict values for another attribute, you should swap columns to move it to the end.)");
		pTransduce->add("[unlabeled-set]", "The filename of a dataset in \".arff\" format. This dataset must have placeholder labels. The placeholder labels will be replaced in the output with the new predicted labels.");
		UsageNode* pDO1 = pTransduce->add("<data_opts1>");
		pDO1->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO1->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		UsageNode* pDO2 = pTransduce->add("<data_opts1>");
		pDO2->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO2->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pTransAcc = pRoot->add("transacc <options> [training-set] <data_opts1> [test-set] <data_opts2> [algorithm]", "Measure the transductive accuracy of [algorithm] with respect to the specified training and test sets. Results are printed to stdout for each dimension in the label vector. Predictive accuracy is reported for nominal label dimensions, and mean-squared-error is reported for continuous label dimensions.");
	{
		UsageNode* pOpts = pTransAcc->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pTransAcc->add("[training-set]", "The filename of a dataset in \".arff\" format. The labels in this dataset are used to infer labels for the unlabeled set. The last attribute (or the last n attributes if there are multiple label dims) is assumed to be the label. (If you wish to predict values for another attribute, you should swap columns to move it to the end.)");
		pTransAcc->add("[test-set]", "The filename of a dataset in \".arff\" format. This dataset must have placeholder labels. The placeholder labels will be replaced in the output with the new predicted labels.");
		UsageNode* pDO1 = pTransAcc->add("<data_opts1>");
		pDO1->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO1->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		UsageNode* pDO2 = pTransAcc->add("<data_opts1>");
		pDO2->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO2->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pSplitTest = pRoot->add("splittest <options> [dataset] <data_opts> [algorithm]", "This shuffles the data, then splits it into two parts, trains with one part, and tests with the other. (This also works with model-free algorithms.) Results are printed to stdout for each dimension in the label vector. Predictive accuracy is reported for nominal label dimensions, and mean-squared-error is reported for continuous label dimensions.");
	{
		UsageNode* pOpts = pSplitTest->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pOpts->add("-trainratio [value]", "Specify the amount of the data (between 0 and 1) to use for training. The rest will be used for testing.");
		pOpts->add("-reps [value]", "Specify the number of repetitions to perform. The default is 1.");
		UsageNode* pDO = pSplitTest->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pCV = pRoot->add("crossvalidate <options> [dataset] <data_opts> [algorithm]", "Perform cross-validation with the specified dataset and algorithm. Results are printed to stdout. (Supports model-free algorithms too.)");
	{
		UsageNode* pOpts = pCV->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pOpts->add("-reps [value]", "Specify the number of repetitions to perform. The default is 5.");
		pOpts->add("-folds [value]", "Specify the number of folds to use. The default is 2.");
		pOpts->add("-succinct", "Just report the average accuracy. Do not report deviation, or results at each fold.");
		pCV->add("[dataset]", "The filename of a dataset in \".arff\" format. The last attribute (or the last n attributes if there are multiple label dims) is assumed to be the label. (If you wish to predict values for another attribute, you should swap columns to move it to the end.)");
		UsageNode* pDO = pCV->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pPR = pRoot->add("precisionrecall <options> [dataset] <data_opts> [algorithm]", "Compute the precision/recall for a dataset and algorithm");
	{
		UsageNode* pOpts = pPR->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pOpts->add("-labeldims [n]", "Specify the number of dimensions in the label (output) vector. The default is 1. (Don't confuse this with the number of class labels. It only takes one dimension to specify a class label, even if there are k possible labels.)");
		pOpts->add("-reps [n]", "Specify the number of reps to perform. More reps means it will take longer, but results will be more accurate. The default is 5.");
		pOpts->add("-samples [n]", "Specify the granularity at which to measure recall. The default is 100.");
		pPR->add("[dataset]", "The filename of a dataset in \".arff\" format. The last attribute (or the last n attributes if there are multiple label dims) is assumed to be the label. (If you wish to predict values for another attribute, you should swap columns to move it to the end.)");
		UsageNode* pDO = pPR->add("<data_opts>");
		pDO->add("-labels [attr_list]", "Specify which attributes to use as labels. (If not specified, the default is to use the last attribute for the label.) [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
		pDO->add("-ignore [attr_list]", "Specify attributes to ignore. [attr_list] is a comma-separated list of zero-indexed attributes. A hyphen may be used to specify a range of values. Example: 0,2-5,7");
	}
	UsageNode* pTS = pRoot->add("trainsparse <options> [sparse-matrix] [algorithm]", "Train the specified algorithm with the sparse matrix. Only incremental learners (such as naivebayes or neuralnet) support this functionality. It will print the trained model-file to stdout.");
	{
		UsageNode* pOpts = pTS->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pOpts->add("-labeldims [n]", "Specify the number of dimensions in the label (output) vector. The default is 1. (Don't confuse this with the number of class labels. It only takes one dimension to specify a class label, even if there are k possible labels.)");
	}
	UsageNode* pPS = pRoot->add("predictsparse <options> [model-file] [sparse-matrix]", "Predict labels for all of the rows in [sparse-matrix]. Label predictions for each row are printed to stdout. (The features are not printed with the predictions.)");
	{
		UsageNode* pOpts = pPS->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pPS->add("[sparse-matrix]", "A dataset in the form of a sparse matrix. The dataset must have placeholder labels, although these labels will be replaced with predicted labels. (The last attribute is assumed to be the class label.)");
		pTest->add("[model-file]", "The filename of a trained model. (This is the file to which you saved the output when you trained a supervised learning algorithm.) Only incremental learning algorithms are supported.");
	}
	UsageNode* pTR = pRoot->add("trainrecurrent <options> [method] [obs-data] [action-data] [context-dims] [algorithm] [algorithm]", "Train a recurrent model of a dynamical system with the specified training [method]. The training data is specified by [obs-data], which specifies the sequence of observations, and [action-data], which specifies the sequence of actions. [context-dims] specifies the number of dimensions in the state-space of the system. The two algorithms specify the two functions of a model of a dynamical system. The first [algorithm] models the transition function. The second [algorithm] models the observation function.");
	{
		UsageNode* pOpts = pTR->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator. (Use this option to ensure that your results are reproduceable.)");
		pOpts->add("-paramdims 2 [wid] [hgt]", "If observations are images, use this option to parameterize the predictions, so only the channel values of each pixel are predicted. (Other values besides 2 dimensions are also supported.)");
		pOpts->add("-state [filename]", "Save the estimated state to the specified file. (Only has effect if moses is used as the training method.)");
		pOpts->add("-validate [interval] 1 [obs] [action]", "Perform validation at [interval]-second intervals with observation data, [obs], and action data, [action]. (Also supports more than 1 validation sequence if desired.)");
		pOpts->add("-out [filename]", "Save the resulting model to the specified file. The default is \"model.twt\".");
		pOpts->add("-noblur", "Do not use blurring. The default is to use blurring. Sometimes blurring improves results. Sometimes not.");
		pOpts->add("-traintime", "Specify how many seconds to train the model. The default is 1 hour.");
		pOpts->add("-isomap", "Use Isomap instead of Breadth-first Unfolding if moses is used as the training method.");

		UsageNode* pMeth = pTR->add("[method]");
		pMeth->add("moses", "Use Temporal-NLDR to estimate state, then build the model using the state estimate.");
		pMeth->add("bptt [depth] [iters-per-grow-sequence]", "Backpropagation Through Time. [depth] specifies the number of instances of the transition function that will appear in the unfolded model. A good value might be 3. [iters-per-grow-sequence] specifies the number of pattern presentations before the sequence is incremented. A good value might be 50000.");
		pMeth->add("evolutionary", "Train with evoluationary optimization");
		pMeth->add("hillclimber", "Train with a hill-climbing algorithm.");
		pMeth->add("annealing [deviation] [decay] [window]", "Train with simulated annealing. Good values might be 2.0 0.5 300");
	}
	return pRoot;
}

UsageNode* makeAlgorithmUsageTree()
{
	UsageNode* pRoot = new UsageNode("[algorithm]", "A supervised learning algorithm, or a transductive algorithm.");
	UsageNode* pAT = pRoot->add("agglomerativetransducer <options>", "A model-free transduction algorithm based on agglomerative clustering. Unlabeled patterns take the label of the cluster into which they are joined. It never joins clusters with different labels. It continues until there are the same number of clusters as labels.");
	{
		UsageNode* pOpts = pAT->add("<options>");
		pOpts->add("-balancingfactor [value]", "The balancing factor is used to prevent large clusters from overpowering small clusters. If [value] is large (like 3.5), it favors joining small clusters together. If [value] is small (like 0.2), then distance between clusters will have more influence.");
	}
	UsageNode* pBag = pRoot->add("bag <contents> end", "A bagging (bootstrap aggregating) ensemble. This is a way to combine the power of many learning algorithms through voting. \"end\" marks the end of the ensemble contents.");
	{
		UsageNode* pContents = pBag->add("<contents>");
		pContents->add("[instance_count] [algorithm]", "Specify the number of instances of a learning algorithm to add to the bagging ensemble.");
	}
	pRoot->add("baseline", "This is one of the simplest of all supervised algorithms. It ignores all features. For nominal labels, it always predicts the most common class in the training set. For continuous labels, it always predicts the mean label in the training set. An effective learning algorithm should never do worse than baseline--hence the name \"baseline\".");
	UsageNode* pBucket = pRoot->add("bucket <contents> end", "This uses cross-validation with the training set to select the best model from a bucket of models. When accuracy is measured across multiple datasets, it will usually do better than the best model in the bucket could do. \"end\" marks the end of the contents of the bucket.");
	{
		UsageNode* pContents = pBucket->add("<contents>");
		pContents->add("[algorithm]", "Add an algorithm to the bucket");
	}
	UsageNode* pCat = pRoot->add("nominaltocat <options> [algorithm]", "This is a filter that wraps another model. It is sort-of the opposite of discretize. It converts each nominal attribute to a categorical distribution by representing each value using the corresponding row of the identity matrix. For example, if a certain nominal attribute has 4 possible values, then a value of 3 would be encoded as the vector 0 0 1 0. When predictions are converted back to nominal values, the mode of the categorical distribution is used as the predicted value. (This is similar to Weka's NominalToBinaryFilter.)");
	{
/*
"_________________               _______________________________",
"|  Data         |               |  NominalToCat               |",
"|Features       |               |                ____________ |",
"|1)size         |--5.3(Real)--->|----5.3(Real)-->|          | |",
"|               |               | /--0.0(Real)-->|          | |",
"|2)color {r,g,b}|--g(Nominal)-->|----1.0(Real)-->|Supervised| |",
"|               |               | \\--0.0(Real)-->| Learner  | |",
"|               |               |                |          | |",
"|Label          |               | /--0.0(Real)<--|          | |",
"|3)class {a,b,c}|<--c(Nominal)--|----0.0(Real)<--|          | |",
"|               |               | \\--1.0(Real)<--|__________| |",
"|_______________|               |_____________________________|",
*/
		UsageNode* pOpts = pCat->add("<options>");
		pOpts->add("-maxvalues [cap]", "If a nominal attribute has more than [cap] possible values, then it will simply enumerate them instead of creating a separate dimension for each possible value. This is somewhat of a safety-net to ensure that the dimensionality doesn't explode. The default is 12.");
		pOpts->add("-notfeatures", "Don't convert features (inputs) to categorical distributions.");
		pOpts->add("-notlabels", "Don't convert the labels (outputs) to categorical distributions.");
	}
	UsageNode* pDT = pRoot->add("decisiontree <options>", "A decision tree.");
	{
		UsageNode* pOpts = pDT->add("<options>");
		pOpts->add("-random [draws]", "Use random divisions (instead of divisions that reduce entropy). Random divisions make the algorithm train faster, and also increase model variance, so it is better suited for ensembles, but random divisions also make the decision tree more vulnerable to problems with irrelevant features. [draws] is typically 1, but if you specify a larger value, it will pick the best out of the specified number of random draws.");
		pOpts->add("-leafthresh [n]", "When building the tree, if the number of samples is <= this value, it will stop trying to divide the data and will create a leaf node. The default value is 1. For noisy data, larger values may be advantageous.");
	}
	UsageNode* pDisc = pRoot->add("discretize <options> [algorithm]", "This is a filter that wraps another model. It discretizes the features before passing them to the model, and then undiscretizes any predicted labels before passing them back to the caller. This is useful with algorithms that only support nominal features, and/or labels, like naivebayes.");
	{
		UsageNode* pOpts = pDisc->add("<options>");
		pOpts->add("-buckets [count]", "Specify the number of buckets to use. The default is the floor of the square root of the number of rows in the training set.");
		pOpts->add("-notfeatures", "Don't discretize the features (inputs).");
		pOpts->add("-notlabels", "Don't discretize the labels (outputs).");
	}
	UsageNode* pGCT = pRoot->add("graphcuttransducer [neighbors]", "This is a model-free transduction algorithm. It uses a min-cut/max-flow graph-cut algorithm to separate each label from all of the others.");
	{
		pGCT->add("[neighbors]", "The number of neighbors to connect with each point in order to form the graph.");
	}
	UsageNode* pKNN = pRoot->add("knn [k] <options>", "The k-Nearest-Neighbor instance-based learning algorithm. It uses Euclidean distance for continuous features and Hamming distance for nominal features.");
	{
		pKNN->add("[k]", "The number of neighbors to use");
		UsageNode* pOpts = pKNN->add("<options>");
		pOpts->add("-equalweight", "Give equal weight to every neighbor. (The default is to use linear weighting for continuous features, and sqared linear weighting for nominal features.");
		pOpts->add("-scalefeatures", "Use a hill-climbing algorithm on the training set to scale the feature dimensions in order to give more accurate results. This increases training time, but also improves accuracy and robustness to irrelevant features.");
	}
	pRoot->add("meanmarginstree", "This is a very simple linear combination tree. (A powerful model can be created using a bagging ensemble of buckets, that each contain one decision tree and one mean margins tree. This combination has been shown to do better than even much larger ensembles of random trees.)");
	UsageNode* pNB = pRoot->add("naivebayes <options>", "The naive Bayes learning algorithm. Supports only discrete features and labels, so it is common to wrap it with a discretizing filter.");
	{
		UsageNode* pOpts = pNB->add("<options>");
		pOpts->add("-ess [value]", "Specifies an equivalent sample size to prevent unsampled values from dominating the joint distribution. Good values typically range between 0 and 1.5.");
	}
	pRoot->add("naiveinstance [neighbors]", "This is an instance learner that assumes each dimension is conditionally independant from other dimensions. It lacks the accuracy of knn in low dimensional feature space, but scales much better to high dimensionality. It only supports continuous values, so it is common to wrap it in a nominaltocat filter.");
	UsageNode* pNT = pRoot->add("neighbortransducer [neighbors] <options>", "This is a model-free transduction algorithm. It is an instance learner that propagates labels where the neighbors are most in agreement. This algorithm does well when classes sample a manifold (such as with text recognition).");
	{
		UsageNode* pOpts = pNT->add("<options>");
		pOpts->add("-friends [intrinsic-dims] [alpha] [beta]", "Use the manifold-friend-finding algorithm instead of the nearest Euclidean neighbors.");
		pOpts->add("-prune", "Prune shortcuts. (Only effective if used with the -friends option.)");
	}
	UsageNode* pNN = pRoot->add("neuralnet <options>", "A single or multi-layer feed-forward neural network of logistic units. It is trained with online backpropagation. Only continuous values are supported, so it is common to wrap it in a nominaltocat filter so it can handle discrete attributes too. It is also common to wrap that in a normalizing filter, to ensure that any continuous inputs are within a reasonable range.");
	{
		UsageNode* pOpts = pNN->add("<options>");
		pOpts->add("-addlayer [size]", "Add a hidden layer with \"size\" logisitic units to the network. You may use this option multiple times to add multiple layers. The first layer added is adjacent to the input features. The last layer added is adjacent to the output labels. If you don't add any hidden layers, the network is just a single layer of sigmoid units.");
		pOpts->add("-learningrate [value]", "Specify a value for the learning rate. The default is 0.1");
		pOpts->add("-momentum [value]", "Specifies a value for the momentum. The default is 0.0");
		pOpts->add("-windowepochs [value]", "Specifies the number of training epochs that are performed before the stopping criteria is tested again. Bigger values will result in a more stable stopping criteria. Smaller values will check the stopping criteria more frequently.");
		pOpts->add("-minwindowimprovement [value]", "Specify the minimum improvement that must occur over the window of epochs for training to continue. [value] specifies the minimum decrease in error as a ratio. For example, if value is 0.02, then training will stop when the mean squared error does not decrease by two percent over the window of epochs. Smaller values will typically result in longer training times.");
		pOpts->add("-dontsquashoutputs", "Don't squash the outputs values with the logistic function. Just report the net value at the output layer. This is often used for regression.");
		pOpts->add("-crossentropy", "Use cross-entropy instead of squared-error for the error signal.");
	}
	UsageNode* pNorm = pRoot->add("normalize <options> [algorithm]", "This is a filter that wraps another model. It scales and translates real-valued data to fall within a specified range. Nominal attributes pass through unaltered.");
	{
		UsageNode* pOpts = pNorm->add("<options>");
		pOpts->add("-range [min] [max]", "Specify the target range. The default is 0 1.");
		pOpts->add("-notfeatures", "Don't normalize the features (inputs).");
		pOpts->add("-notlabels", "Don't normalize the labels (outputs).");
	}
	UsageNode* pPCA = pRoot->add("pca [target-dims] [algorithm]", "This is a filter that wraps another model. It performs principle component analysis to reduce the dimensionality of the features (inputs) before passing them to the wrapped model.");
	{
		pPCA->add("[target-dims]", "The number of target dimensions in the feature vector.");
	}

	return pRoot;
}

UsageNode* makeTransformUsageTree()
{
	UsageNode* pRoot = new UsageNode("waffles_transform [command]", "Transform data, reduce dimensionality, cluster, shuffle rows, swap columns, matrix operations, etc.");
	UsageNode* pAdd = pRoot->add("add [dataset1] [dataset2]", "Adds two matrices together element-wise. Results are printed to stdout.");
	{
		pAdd->add("[dataset1]", "The filename of the first matrix in ARFF format.");
		pAdd->add("[dataset2]", "The filename of the second matrix in ARFF format.");
	}
	UsageNode* pAIC = pRoot->add("addindexcolumn [dataset] <options>", "Add a column that Specify the index of each row.");
	{
		pAIC->add("[dataset]", "The filename of a dataset in ARFF format");
		UsageNode* pOpts = pAIC->add("<options>");
		pOpts->add("-start [value]", "Specify the initial index. (the default is 0).");
		pOpts->add("-increment [value]", "Specify the increment amount. (the default is 1).");
	}
	UsageNode* pAN = pRoot->add("addnoise [dataset] [dev] <options>", "Add Gaussian noise with the specified deviation to all the elements in the dataset. (Assumes that the values are all continuous.)");
	{
		pAN->add("[dataset]", "The filename of a dataset in ARFF format");
		pAN->add("[dev]", "The deviation of the Gaussian noise");
		UsageNode* pOpts = pAIC->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-excludelast [n]", "Do not add noise to the last [n] columns.");
	}
	pRoot->add("agglomerative [dataset] [clusters]", "Performs agglomerative clustering. Outputs the cluster id for each row.");
	UsageNode* pAlign = pRoot->add("align [a] [b]", "Translates and rotates dataset [b] to minimize mean squared difference with dataset [a]. (Uses the Kabsch algorithm.)");
	{
		pAlign->add("[a]", "The filename of a dataset in ARFF format");
		pAlign->add("[b]", "The filename of a dataset in ARFF format");
	}
	UsageNode* pAS = pRoot->add("attributeselector [dataset] <options>", "Make a ranked list of attributes from most to least salient. The ranked list is printed to stdout. Attributes are zero-indexed.");
	{
		pAS->add("[dataset]", "The filename of a dataset in ARFF format");
		UsageNode* pOpts = pAS->add("<options>");
		pOpts->add("-out [n] [filename]", "Save a dataset containing only the [n]-most salient features to [filename].");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-labeldims [n]", "Specify the number of dimensions in the label (output) vector. The default is 1. (Don't confuse this with the number of class labels. It only takes one dimension to specify a class label, even if there are k possible labels.)");
	}
	pRoot->add("autocorrelation [dataset]", "Compute the autocorrelation of the specified time-series data.");
	UsageNode* pBE = pRoot->add("blendembeddings [data-orig] [neighbor-finder] [data-a] [data-b] <options>", "Compute a blended \"average\" embedding from two reduced-dimensionality embeddings of some data.");
	{
		pBE->add("[data-orig]", "The filename of the original high-dimensional data in ARFF format.");
		pBE->add("[data-a]", "The first reduced dimensional embedding of [data-orig]");
		pBE->add("[data-b]", "The second reduced dimensional embedding of [data-orig]");
		UsageNode* pOpts = pBE->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pBFU = pRoot->add("breadthfirstunfolding [dataset] [neighbor-finder] [target_dims] <options>", "A manifold learning algorithm.");
	{
		UsageNode* pOpts = pBFU->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-reps [n]", "The number of times to compute the embedding and blend the results together. The default is 1.");
	}
	UsageNode* pCat = pRoot->add("nominaltocat [dataset] <options>", "Convert the dataset to a categorical distribution of all real values. Columns with only two nominal values are converted to 0 or 1. If there are three or more possible values, a column is created for each value. The column corresponding to the value is set to 1, and the others are set to 0. (This is similar to Weka's NominalToBinaryFilter.)");
	{
		UsageNode* pOpts = pCat->add("<options>");
		pOpts->add("-maxvalues [cap]", "Specify the maximum number of nominal values for which to create new columns. The default is 12.");
	}
	pRoot->add("cholesky [dataset]", "Compute the cholesky decomposition of the specified matrix.");
	UsageNode* pCorr = pRoot->add("correlation [dataset] [attr1] [attr2] <options>", "Compute the linear correlation coefficient of the two specified attributes.");
	{
		pCorr->add("[dataset]", "The filename of a dataset in ARFF format.");
		pCorr->add("[attr1]", "A zero-indexed attribute number.");
		pCorr->add("[attr2]", "A zero-indexed attribute number.");
		UsageNode* pOpts = pCorr->add("<options>");
		pOpts->add("-aboutorigin", "Compute the correlation about the origin. (The default is to compute it about the mean.)");
	}
	UsageNode* pDet = pRoot->add("determinant [dataset]", "Compute the determinant of the specified matrix.");
	{
		pDet->add("[dataset]", "The filename of a dataset in ARFF format.");
	}
	UsageNode* pDisc = pRoot->add("discretize [dataset] <options>", "Discretizes the continuous attributes in the specified dataset.");
	{
		pDisc->add("[dataset]", "The filename of a dataset in ARFF format.");
		UsageNode* pOpts = pDisc->add("<options>");
		pOpts->add("-buckets [count]", "Specify the number of buckets to use. The default is the square root of the number of rows in the dataset.");
		pOpts->add("-colrange [first] [last]", "Specify a range of columns. Only continuous columns in the specified range will be modified. (Columns are zero-indexed.)");
	}
	UsageNode* pDropCols = pRoot->add("dropcolumns [dataset] [column-list]", "Remove one or more columns from a dataset and prints the results to stdout. (The input file is not modified.)");
	{
		pDropCols->add("[column-list]", "A list of zero-indexed columns to drop. Example: 0,2-5,7");
	}
	pRoot->add("dropmissingvalues [dataset]", "Remove all rows that contain missing values.");
	UsageNode* pEx = pRoot->add("export [dataset] <options>", "Print the data as a list of comma separated values without any meta-data.");
	{
		UsageNode* pOpts = pEx->add("<options>");
		pOpts->add("-tab", "Separate with tabs instead of commas.");
		pOpts->add("-space", "Separate with spaces instead of commas.");
	}
	pRoot->add("droprows [dataset] [after-size]", "Removes all rows except for the first [after-size] rows.");
	UsageNode* pIm = pRoot->add("import [dataset] <options>", "Convert a text file of comma separated (or otherwise separated) values to a .arff file. The meta-data is automatically determined. The .arff file is printed to stdout. This makes it easy to operate on structured data from a spreadsheet, database, or pretty-much any other source.");
	{
		UsageNode* pOpts = pIm->add("<options>");
		pOpts->add("-tab", "Data elements are separated with a tab character instead of a comma.");
		pOpts->add("-space", "Data elements are separated with a single space instead of a comma.");
		pOpts->add("-whitespace", "Data elements are separated with an arbitrary amount of whitespace.");
		pOpts->add("-semicolon", "Data elements are separated with semicolons.");
		pOpts->add("-separator [char]", "Data elements are separated with the specified character.");
		pOpts->add("-columnnames", "Use the first row of data for column names.");
	}
	UsageNode* pIsomap = pRoot->add("isomap [dataset] [neighbor-finder] [target_dims] <options>", "Use the Isomap algorithm to reduce dimensionality.");
	{
		UsageNode* pOpts = pIsomap->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pKM = pRoot->add("kmeans [dataset] [clusters]", "Performs k-means clustering. Outputs the cluster id for each row.");
	{
		UsageNode* pOpts = pKM->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	pRoot->add("kmedoids [dataset] [clusters]", "Performs k-medoids clustering. Outputs the cluster id for each row.");
	UsageNode* pLLE = pRoot->add("lle [dataset] [neighbor-finder] [target_dims] <options>", "Use the LLE algorithm to reduce dimensionality.");
	{
		UsageNode* pOpts = pLLE->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pMS = pRoot->add("manifoldsculpting [dataset] [neighbor-finder] [target_dims] <options>", "Use the Manifold Sculpting algorithm to reduce dimensionality.");
	{
		UsageNode* pOpts = pMS->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-continue [dataset]", "Continue refining the specified reduced-dimensional results. (This feature enables Manifold Sculpting to improve upon its own results, or to refine the results from another dimensionality reduction algorithm.)");
		pOpts->add("-scalerate [value]", "Specify the scaling rate. The default is 0.999. A value close to 1 will give better results, but will cause the algorithm to take longer.");
	}
	UsageNode* pMSE = pRoot->add("measuremeansquarederror [dataset1] [dataset2] <options>", "Print the mean squared error between two datasets. (Both datasets must be the same size.)");
	{
		UsageNode* pOpts = pMSE->add("<options>");
		pOpts->add("-fit", "Use a hill-climber to find an affine transformation to make dataset2 fit as closely as possible to dataset1. Report results after each iteration.");
	}
	pRoot->add("mergehoriz [dataset1] [dataset2]", "Merge two (or more) datasets horizontally. All datasets must already have the same number of rows. The resulting dataset will have all the columns of both datasets.");
	pRoot->add("mergevert [dataset1] [dataset2]", "Merge two datasets vertically. Both datasets must already have the same number of columns. The resulting dataset will have all the rows of both datasets.");
	UsageNode* pMDS = pRoot->add("multidimensionalscaling [distance-matrix] [target-dims]", "Perform MDS on the specified [distance-matrix].");
	{
		pMDS->add("[distance-matrix]", "The filename of an arff file that contains the pair-wise distances (or dissimilarities) between every pair of points. It must be a square matrix of real values. Only the upper-triangle of this matrix is actually used. The lower-triangle and diagonal is ignored.");
		UsageNode* pOpts = pMDS->add("<options>");
		pOpts->add("-squareddistances", "The distances in the distance matrix are squared distances, instead of just distances.");
	}
	UsageNode* pMult1 = pRoot->add("multiply [a] [b] <options>", "Matrix multiply [a] x [b]. Both arguments are the filenames of .arff files. Results are printed to stdout.");
	{
		UsageNode* pOpts = pMult1->add("<options>");
		pOpts->add("-transposea", "Transpose [a] before multiplying.");
		pOpts->add("-transposeb", "Transpose [b] before multiplying.");
	}
	pRoot->add("multiplyscalar [dataset1] [scalar]", "Multiply all elements in [dataset1] by the specified scalar. Results are printed to stdout.");
	UsageNode* pNorm = pRoot->add("normalize [dataset] <options>", "Normalize all continuous attributes to fall within the specified range. (Nominal columns are left unchanged.)");
	{
		UsageNode* pOpts = pNorm->add("<options>");
		pOpts->add("-range [min] [max]", "Specify the output min and max values. (The default is 0 1.)");
	}
	pRoot->add("pca [dataset] [target_dims]", "Projects the data into the specified number of dimensions with principle component analysis. (Prints results to stdout. The input file is not modified.)");
	UsageNode* pPC = pRoot->add("principalcomponents [dataset] [k] <options>", "Compute the first k principal components of [dataset] and print the results to stdout.");
	{
		UsageNode* pOpts = pPC->add("<options>");
		pOpts->add("-aboutorigin", "Compute the principal components about the origin. (The default is to compute them about the mean.)");
	}
	pRoot->add("pseudoinverse [dataset]", "Compute the Moore-Penrose pseudo-inverse of the specified matrix of real values.");
	UsageNode* pRMS = pRoot->add("replacemissingvalues [dataset] <options>", "Replace any missing values with another randomly chosen value from the dataset.");
	{
		UsageNode* pOpts = pRMS->add("<options>");
		pOpts->add("-seed", "Specify a seed for the random number generator.");
	}
	pRoot->add("reducedrowechelonform [dataset]", "Convert a matrix to reduced row echelon form. Results are printed to stdout.");
	UsageNode* pShuffle = pRoot->add("shuffle [dataset] <options>", "Shuffle the row order.");
	{
		pShuffle->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pSignif = pRoot->add("significance [dataset] [attr1] [attr2] <options>", "Compute statistical significance values for the two specified attributes.");
	{
		pSignif->add("[dataset]", "The filename of a .arff file.");
		pSignif->add("[attr1]", "A zero-indexed column number.");
		pSignif->add("[attr2]", "A zero-indexed column number.");
		UsageNode* pOpts = pSignif->add("<options>");
		pOpts->add("-tol [value]", "Sets the tolerance value for the Wilcoxon Signed Ranks test. The default value is 0.001.");
	}
	UsageNode* pSortCol = pRoot->add("sortcolumn [dataset] [col] <options>", "Sort the rows in [dataset] such that the values in the specified column are in ascending order and print the results to to stdout. (The input file is not modified.)");
	{
		pSortCol->add("[dataset]", "The filename of a dataset in ARFF format.");
		pSortCol->add("[col]", "The zero-indexed column number in which to sort");
		UsageNode* pOpts = pSortCol->add("<options>");
		pOpts->add("-descending", "Sort in descending order instead of ascending order.");
	}
	UsageNode* pSS = pRoot->add("sparsehuffle [sparse-matrix] <options>", "Shuffles the row order or a sparse matrix.");
	{
		pSS->add("[sparse-matrix]", "The filename of a sparse matrix (not an ARFF file).");
		UsageNode* pOpts = pSS->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pSSplit = pRoot->add("sparsesplit [sparse-matrix] [rows] [filename1] [filename2]", "Splits a sparse matrix into two datasets. Nothing is printed to stdout.");
	{
		pSSplit->add("[sparse-matrix]", "The filename of a sparse matrix.");
		pSSplit->add("[rows]", "The number of rows to put in the first file. The rest go in the second file.");
	}
	UsageNode* pSplit = pRoot->add("split [dataset] [rows] [filename1] [filename2]", "Split a dataset into two datasets. (Nothing is printed to stdout.)");
	{
		pSplit->add("[dataset]", "The filename of a datset in .arff format.");
		pSplit->add("[rows]", "The number of rows to go into the first file. The rest go in the second file.");
	}
	pRoot->add("squareddistance [a] [b]", "Computes the sum and mean squared distance between dataset [a] and [b]. ([a] and [b] are each the names of files in .arff format. They must have the same dimensions.)");
	UsageNode* pSVD = pRoot->add("svd [matrix] <options>", "Compute the singular value decomposition of a matrix.");
	{
		pSVD->add("[matrix]", "A .arff file containing the matrix values.");
		UsageNode* pOpts = pSVD->add("<options>");
		pOpts->add("-ufilename [filename]", "Set the filename to which U will be saved. U is the matrix in which the columns are the eigenvectors of [matrix] times its transpose. The default is u.arff.");
		pOpts->add("-sigmafilename [filename]", "Set the filename to which Sigma will be saved. Sigma is the matrix that contains the singular values on its diagonal. All values in Sigma except the diagonal will be zero. If this option is not specified, the default is to only print the diagonal values (not the whole matrix) to stdout. If this options is specified, nothing is printed to stdout.");
		pOpts->add("-vfilename [filename]", "Set the filename to which V will be saved. V is the matrix in which the row are the eigenvectors of the transpose of [matrix] times [matrix]. The default is v.arff.");
		pOpts->add("-maxiters [n]", "Specify the number of times to iterate before giving up. The default is 100, which should be sufficient for most problems.");
	}
	pRoot->add("swapcolumns [dataset] [col1] [col2]", "Swap two columns in the specified dataset and prints the results to stdout. (Columns are zero-indexed.)");
	UsageNode* pTransition = pRoot->add("transition [action-sequence] [state-sequence] <options>", "Given a sequence of actions and a sequence of states (each in separate datasets), this generates a single dataset to map from action-state pairs to the next state. This would be useful for generating the data to train a transition function.");
	{
		UsageNode* pOpts = pTransition->add("<options>");
		pOpts->add("-delta", "Predict the delta of the state transition instead of the new state.");
	}
	pRoot->add("transpose [dataset]", "Transpose the data such that columns become rows and rows become columns.");

	return pRoot;
}


UsageNode* makeGenerateUsageTree()
{
	UsageNode* pRoot = new UsageNode("waffles_generate [command]", "Generate certain useful datasets");
	UsageNode* pCrane = pRoot->add("crane <options>", "Generate a dataset where each row represents a ray-traced image of a crane with a ball.");
	{
		UsageNode* pOpts = pCrane->add("<options>");
		pOpts->add("-saveimage [filename]", "Save an image showing all the frames.");
		pOpts->add("-ballradius [size]", "Specify the size of the ball. The default is 0.3.");
		pOpts->add("-frames [horiz] [vert]", "Specify the number of frames to render.");
		pOpts->add("-size [wid] [hgt]", "Specify the size of each frame.");
		pOpts->add("-blur [radius]", "Blurs the images. A good starting value might be 5.0.");
		pOpts->add("-gray", "Use a single grayscale value for every pixel instead of three (red, green, blue) channel values.");
	}
	pRoot->add("cube [n]", "returns data evenly distributed on the surface of a unit cube. Each side is sampled with [n]x[n] points. The total number of points in the dataset will be 6*[n]*[n]-12*[n]+8.");
	UsageNode* pDocsToSparse = pRoot->add("docstosparsematrix <options> [folder1] [folder2] ...", "Converts a set of documents to a sparse matrix. [folder1] should contain all of the documents in class1. [folder2] should contain all the documents in class2, and so forth. The words are filtered against a common set of stop words. Also, words less than 4 letters are ignored. Currently, only .txt and .html documents are supported. Other file types are ignored. Each row in the sparse matrix represents one of the documents. If there is more than one folder/class, then an additional column is added for the class label. Subdirectories are not followed. Results are saved to a text file in compressed-column format. A mapping from row number to document filename is printed to stdout.");
	{
		UsageNode* pOpts = pDocsToSparse->add("<options>");
		pOpts->add("-nostem", "Specifies not to stem the words. (The default is to use the Porter stemming algorithm.)");
		pOpts->add("-binary", "Just use the value 1 if the word occurs in a document. The default is to compute the somewhat more meaningful value: a*log(b/c)/d, where a=the number of times the word occurs in this document, b=total number of documents, c=number of documents that contain this word, and d=the max number of times this word occurs in any document.");
		pOpts->add("-outfile [filename]", "Specify the filename for the sparse matrix. The default is \"docs.sparse\".");
		pOpts->add("-vocabfile [filename]", "Save the vocabulary of words to the specified file. The default is to not save the list of words. Note that the words will be stemmed (unless -nostem was specified), so it is normal for many of them to appear misspelled.");
	}
	UsageNode* pES = pRoot->add("entwinedspirals [points] <options>", "Generates points that lie on an entwined spirals manifold.");
	{
		pES->add("[points]", "The number of points with which to sample the manifold.");
		UsageNode* pOpts = pES->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-reduced", "Generate intrinsic values instead of extrinsic values. (This might be useful to empirically measure the accuracy of a manifold learner.)");
	}
	UsageNode* pFishBowl = pRoot->add("fishbowl [n] <option>", "Generate samples on the surface of a fish-bowl manifold.");
	{
		UsageNode* pOpts = pFishBowl->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-opening [size]", "the size of the opening. (0.0 = no opening. 0.25 = default. 1.0 = half of the sphere.)");
	}
	UsageNode* pGRW = pRoot->add("gridrandomwalk [arff-file] [width] [samples] <options>", "Generate a sequence of action-observation pairs by randomly walking around on a grid of observation vectors.");
	{
		pGRW->add("[arff-file]", "The filename of an arff file containing observation vectors arranged in a grid.");
		pGRW->add("[width]", "The width of the grid.");
		pGRW->add("[samples]", "The number of samples to take. In other words, the length of the random walk.");
		UsageNode* pOpts = pGRW->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-start [x] [y]", "Specifies the starting state. The default is to start in the center of the grid.");
		pOpts->add("-obsfile [filename]", "Specify the filename for the observation sequence data. The default is observations.arff.");
		pOpts->add("-actionfile [filename]", "Specify the filename for the actions data. The default is actions.arff.");
	}
	UsageNode* pITON = pRoot->add("imagetranslatedovernoise [png-file] <options>", "Sample a manifold by translating an image over a background of noise.");
	{
		pITON->add("[png-file]", "The filename of a png image.");
		UsageNode* pOpts = pITON->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-reduced", "Generate intrinsic values instead of extrinsic values. (This might be useful to empirically measure the accuracy of a manifold learner.)");
	}
	UsageNode* pNoise = pRoot->add("noise [rows] <options>", "Generate random data by sampling from a distribution.");
	{
		pNoise->add("[rows]", "The number of patterns to generate.");
		UsageNode* pOpts = pNoise->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		UsageNode* pDist = pOpts->add("-dist [distribution]", "Specify the distribution. The default is normal 0 1");
		pDist->add("beta [alpha] [beta]");
		pDist->add("binomial [n] [p]");
		pDist->add("categorical 3 [p0] [p1] [p2]", "A categorical distribution with 3 classes. [p0], [p1], and [p2] specify the probabilities of each of the 3 classes. (This is just an example. Other values besides 3 may be used for the number of classes.)");
		pDist->add("cauchy [median] [scale]");
		pDist->add("chisquare [t]");
		pDist->add("exponential [beta]");
		pDist->add("f [t] [u]");
		pDist->add("gamma [alpha] [beta]");
		pDist->add("gaussian [mean] [deviation]");
		pDist->add("geometric [p]");
		pDist->add("logistic [mu] [s]");
		pDist->add("lognormal [mu] [sigma]");
		pDist->add("normal [mean] [deviation]");
		pDist->add("poisson [mu]");
		pDist->add("softimpulse [s]");
		pDist->add("spherical [dims] [radius]");
		pDist->add("student [t]");
		pDist->add("uniform [a] [b]");
		pDist->add("weibull [gamma]");
	}
	UsageNode* pScalRot = pRoot->add("scalerotate [png-file] <options>", "Generate a dataset where each row represents an image that has been scaled and rotated by various amounts. Thus, these images form an open-cylinder (although somewhat cone-shaped) manifold.");
	{
		UsageNode* pOpts = pScalRot->add("<options>");
		pOpts->add("-saveimage [filename]", "Save a composite image showing all the frames in a grid.");
		pOpts->add("-frames [rotate-frames] [scale-frames]", "Specify the number of frames. The default is 40 15.");
		pOpts->add("-arc [radians]", "Specify the rotation amount. The default is 6.2831853... (2*PI).");
	}
	UsageNode* pSCurve = pRoot->add("scurve [points] <options>", "Generate points that lie on an s-curve manifold.");
	{
		pSCurve->add("[points]", "The number of points with which to sample the manifold");
		UsageNode* pOpts = pSCurve->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-reduced", "Generate intrinsic values instead of extrinsic values. (This might be useful to empirically measure the accuracy of a manifold learner.)");
	}
	UsageNode* pSIR = pRoot->add("selfintersectingribbon [points] <options>", "Generate points that lie on a self-intersecting ribbon manifold.");
	{
		pSIR->add("[points]", "The number of points with which to sample the manifold.");
		UsageNode* pOpts = pSIR->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
	}
	UsageNode* pSR = pRoot->add("swissroll [points] <options>", "Generate points that lie on a swiss roll manifold.");
	{
		pSR->add("[points]", "The number of points with which to sample the manifold.");
		UsageNode* pOpts = pSR->add("<options>");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-reduced", "Generate intrinsic values instead of extrinsic values. (This might be useful to empirically measure the accuracy of a manifold learner.)");
		pOpts->add("-cutoutstar", "Don't sample within a star-shaped region on the manifold.");
	}
	UsageNode* pWI = pRoot->add("windowedimage [png-file] <options>", "Sample a manifold by translating a window over an image. Each pattern represents the windowed portion of the image.");
	{
		pWI->add("[png-file]", "The filename of the png image from which to generate the data.");
		UsageNode* pOpts = pWI->add("<options>");
		pOpts->add("-reduced", "Generate intrinsic values instead of extrinsic values. (This might be useful to empirically measure the accuracy of a manifold learner.)");
		pOpts->add("-stepsizes [horiz] [vert]", "Specify the horizontal and vertical step sizes. (how many pixels to move the window between samples.)");
		pOpts->add("-windowsize [width] [height]", "Specify the size of the window. The default is half the width and height of [png-file].");
	}

	return pRoot;
}

UsageNode* makePlotUsageTree()
{
	UsageNode* pRoot = new UsageNode("waffles_plot [command]", "Visualize data, plot functions, make charts, etc.");
	UsageNode* p3d = pRoot->add("3d [dataset] <options>", "Make a 3d scatter plot. Points are colored with a spectrum according to their order in the dataset.");
	{
		p3d->add("[dataset]", "The filename of an arff file containing the data to plot. It must have exactly 3 continuous attributes.");
		UsageNode* pOpts = p3d->add("<options>");
		pOpts->add("-blast", "Produce a 5-by-5 grid of renderings, each time using a random point of view. It will print the random camera directions that it selects to stdout.");
		pOpts->add("-seed [value]", "Specify a seed for the random number generator.");
		pOpts->add("-size [width] [height]", "Sets the size of the image. The default is 1000 1000.");
		pOpts->add("-pointradius [radius]", "Set the size of the points. The default is 40.0.");
		pOpts->add("-bgcolor [color]", "Set the background color. The default is ffffff.");
		pOpts->add("-cameradistance [dist]", "Set the distance between the camera and the mean of the data. This value is specified as a factor, which is multiplied by the distance between the min and max corners of the data. The default is 1.5. (If the camera is too close to the data, make this value bigger.)");
		pOpts->add("-cameradirection [dx] [dy] [dz]", "Specifies the direction from the camera to the mean of the data. (The camera always looks at the mean.) The default is 0.6 -0.3 -0.8.");
		pOpts->add("-out [filename]", "Specify the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
		pOpts->add("-nolabels", "Don't put axis labels on the bounding box.");
		pOpts->add("-nobox", "Don't draw a bounding box around the plot.");
	}
	UsageNode* pBar = pRoot->add("bar [dataset] <options>", "Make a bar chart.");
	{
		pBar->add("[dataset]", "The filename of an arff file containing the data for the bar chart. The dataset must contain exactly one continuous attribute. Each data row specifies the height of a bar.");
		UsageNode* pOpts = pBar->add("<options>");
		pOpts->add("-log", "Use a logarithmic scale.");
		pOpts->add("-out [filename]", "Specifies the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
	}
	pRoot->add("bigo [dataset]", "Estimate the Big-O runtime of algorithms based on empirical results. Regresses the formula t=a*(n^b+c) to fit the data, where n is the value in attribute 0 (representing the size of the data), and t (representing time) in the other attributes for each algorithm. The values of a, b, and c are reported for each attribute > 0.");
	UsageNode* pEquat = pRoot->add("equation <options> [equations]", "Plot an equation (or multiple equations) in 2D");
	{
		UsageNode* pOpts = pEquat->add("<options>");
		pOpts->add("-out [filename]", "Specify the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
		pOpts->add("-size [width] [height]", "Specify the size of the chart. (The default is 1024 1024.)");
		pOpts->add("-range [xmin] [ymin] [xmax] [ymax]", "Set the range. (The default is: -10 -10 10 10.)");
		pEquat->add("[equations]", "A set of equations separated by semicolons. Since '^' is a special character for many shells, it's usually a good idea to put your equations inside quotation marks. Here are some examples:\n"
		"\"f1(x)=3*x+2\"\n"
		"\"f1(x)=(g(x)+1)/g(x); g(x)=sqrt(x)+pi\"\n"
		"\"h(bob)=bob^2;f1(x)=3+bar(x,5)*h(x)-(x/foo);bar(a,b)=a*b-b;foo=3.2\"\n"
		"Only functions that begin with 'f' followed by a number will be plotted, starting with 'f1', and it will stop when the next number in ascending order is not defined. You may define any number of helper functions or constants with any name you like. Built in constants include: e, and pi. Built in functions include: +, -, *, /, %, ^, abs, acos, acosh, asin, asinh, atan, atanh, ceil, cos, cosh, erf, floor, gamma, lgamma, log, max, min, sin, sinh, sqrt, tan, and tanh. These generally have the same meaning as in C, except '^' means exponent, \"gamma\" is the gamma function, and max and min can support any number (>=1) of parameters. (Some of these functions may not not be available on Windows, but most of them are.) You can override any built in constants or functions with your own variables or functions, so you don't need to worry too much about name collisions. Variables must begin with an alphabet character or an underscore. Multiplication is never implicit, so you must use a '*' character to multiply. Whitespace is ignored.");
	}
	UsageNode* pHist = pRoot->add("histogram [dataset] <options>", "Make a histogram.");
	{
		pHist->add("[dataset]", "The filename of an arff file containing the data for the histogram. The dataset must contain exactly one continuous attribute.");
		UsageNode* pOpts = pHist->add("<options>");
		pOpts->add("-size [width] [height]", "Specify the size of the chart. (The default is 1024 1024.)");
		pOpts->add("-attr [index]", "Specify which attribute is charted. (The default is 0.)");
		pOpts->add("-out [filename]", "Specify the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
		pOpts->add("-range [xmin] [xmax] [ymax]", "Specify the range of the histogram plot");
	}
	UsageNode* pModel = pRoot->add("model [model-file] [dataset] [attr-x] [attr-y]", "Plot the model space of a trained supervised learning algorithm.");
	{
		pModel->add("[model-file]", "The filename of the trained model. (You can use \"waffles_learn train\" to make a model file.)");
		pModel->add("[dataset]", "The filename of a dataset in arff format to be plotted. It can be the training set that was used to train the model, or a test set that it hasn't yet seen.");
		pModel->add("[attr-x]", "The zero-based index of a continuous feature attributes for the horizontal axis.");
		pModel->add("[attr-y]", "The zero-based index of a continuous feature attributes for the vertical axis.");
	}
	UsageNode* pOver = pRoot->add("overview [dataset]", "Generate a matrix of plots of attribute distributions and correlations. This is a useful chart for becoming acquainted with a dataset.");
	{
		pOver->add("[dataset]", "The filename of a dataset in arff format to be charted.");
		UsageNode* pOpts = pOver->add("<options>");
		pOpts->add("-out [filename]", "Specify the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
		pOpts->add("-cellsize [value]", "Change the size of each cell. The default is 100.");
		pOpts->add("-jitter [value]", "Specify how much to jitter the plotted points. The default is 0.03.");
		pOpts->add("-maxattrs [value]", "Specifies the maximum number of attributes to plot. The default is 20.");
	}
	UsageNode* pPDT = pRoot->add("printdecisiontree [model-file] <dataset>", "Print a textual representation of a decision tree to stdout.");
	{
		pPDT->add("[model-file]", "The filename of a trained decision tree model. (You can make one with the command \"waffles_learn train [dataset] decisiontree > [filename]\".)");
		pPDT->add("<dataset>", "An optional filename of the arff file that was used to train the decision tree. The data in this file is ignored, but the meta-data will be used to make the printed model richer.");
	}
	UsageNode* pScat = pRoot->add("scatter [dataset] <options>", "Makes a scatter plot or line graph.");
	{
		pScat->add("[dataset]", "The filename of a dataset in arff format to be plotted. The first attribute specifies the values on the horizontal axis. All other attributes specify the values on the vertical axis for a certain color.");
		UsageNode* pOpts = pScat->add("<options>");
		pOpts->add("-lines", "Draw lines connecting sequential point in the data. (In other words, make a line graph instead of a scatter plot.)");
		pOpts->add("-size [width] [height]", "Specify the size of the chart. (The default is 1024 1024.)");
		pOpts->add("-logx", "Show the horizontal axis on a logarithmic scale");
		pOpts->add("-logy", "Show the vertical axis on a logarithmic scale");
		pOpts->add("-nohorizaxislabels", "Do not draw any labels (or grid lines) for the horizontal axis.");
		pOpts->add("-novertaxislabels", "Do not draw any labels (or grid lines) for the vertical axis.");
		pOpts->add("-textsize [size]", "Sets the label font size. The default is 2.0.");
		pOpts->add("-pointradius [radius]", "Set the size of the point dots. The default is 7.0.");
		pOpts->add("-linethickness [value]", "Specify the line thickness. (The default is 3.0.)");
		pOpts->add("-range [xmin] [ymin] [xmax] [ymax]", "Sets the range. (The default is to determine the range automatically.)");
		pOpts->add("-aspect", "Adjust the range to preserve the aspect ratio. In other words, make sure that both axes visually have the same scale.");
		pOpts->add("-chartcolors [background] [text] [grid]", "Sets colors for the specified areas. (The default is ffffff 000000 808080.)");
		pOpts->add("-linecolors [c1] [c2] [c3] [c4]", "Sets the colors for the first four attributes. The default is 0000a0 a00000 008000 504010 (blue, red, green, brown). (If there are more than four lines, it will just distribute them evenly over the color spectrum.)");
		pOpts->add("-spectrum", "Instead of giving each line a unique color, this will use the color spectrum to indicate the position of each point within the data.");
		pOpts->add("-specmod [cycle]", "Like -spectrum, except it repeats the spectrum with the specified cycle size.");
		pOpts->add("-out [filename]", "Specifies the name of the output file. (The default is plot.png.) It should have the .png extension because other image formats are not yet supported.");
		pOpts->add("-neighbors [neighbor-finder]", "Draw lines connecting each point with its neighbors as determined by the specified neighbor finding algorithm.");
	}
	UsageNode* pStats = pRoot->add("stats [dataset]", "Prints some basic stats about the dataset to stdout.");
	{
		pStats->add("[dataset]", "The filename of an arff file.");
	}
	return pRoot;
}

UsageNode* makeNeighborUsageTree()
{
	UsageNode* pRoot = new UsageNode("[neighbor-finder]", "A neighbor-finding algorithm.");
	UsageNode* pBF = pRoot->add("bruteforce <options> [k]", "The slow way to find the neareast Euclidean-distance neighbors.");
	{
		UsageNode* pOpts = pBF->add("<options>");
		UsageNode* pCC = pOpts->add("-cyclecut [thresh]", "Use CycleCut to break shortcuts and cycles.");
		pCC->add("[thresh]", "The threshold cycle-length for bad cycles.");
		pOpts->add("-normalize", "Normalize distances in local neighborhoods so that all neighborhoood have a uniform amount of total distance.");
		pBF->add("[k]", "The number of neighbors.");
	}
	UsageNode* pKD = pRoot->add("kdtree <options> [k]", "A faster way to find the neareast Euclidean-distance neighbors.");
	{
		UsageNode* pOpts = pKD->add("<options>");
		UsageNode* pCC = pOpts->add("-cyclecut [thresh]", "Use CycleCut to break shortcuts and cycles.");
		pCC->add("[thresh]", "The threshold cycle-length for bad cycles.");
		pKD->add("[k]", "The number of neighbors.");
	}
	UsageNode* pMan = pRoot->add("manifold <options> [k] [t] [alpha] [beta]", "An intelligent neighbor-finder that finds neighborhoods with aligned tangent hyperplanes.");
	{
		UsageNode* pOpts = pMan->add("<options>");
		UsageNode* pCC = pOpts->add("-cyclecut [thresh]", "Use CycleCut to break shortcuts and cycles.");
		pCC->add("[thresh]", "The threshold cycle-length for bad cycles.");
		pMan->add("[k]", "The number of neighbors.");
		pMan->add("[t]", "The number of dimensions in the tangent hyperplanes.");
		pMan->add("[alpha]", "The first parameter. 2.0 might be a good value.");
		pMan->add("[beta]", "The second parameter. 1.0 might be a good value.");
	}
	UsageNode* pSys = pRoot->add("system <options> [action-data] [k]", "A neighbor-finder designed for modeling dynamical systems.");
	{
		UsageNode* pOpts = pSys->add("<options>");
		UsageNode* pCC = pOpts->add("-cyclecut [thresh]", "Use CycleCut to break shortcuts and cycles.");
		pCC->add("[thresh]", "The threshold cycle-length for bad cycles.");
		pSys->add("[action-data] is the filename of an arff file for the sequence of actions given to the system.");
		pSys->add("[k]", "The number of neighbors.");
	}
	return pRoot;
}
