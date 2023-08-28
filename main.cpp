#include "metis.h"
#include "epanet2_2.h"
#include "json.hpp"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <numeric>
#include <random>
#include <cstdlib>

using namespace std;

struct BurstEvent
{
	int nodeIndex;
	string nodeID;
	int startStep;
	int durStep;
	double burstDemand;
	int burstDay;
};

bool operator < (const BurstEvent& a, const BurstEvent& b) {
	if (a.burstDay != b.burstDay) {
		return a.burstDay < b.burstDay;
	} else if (a.startStep != b.startStep) {
		return a.startStep < b.startStep;
	} else if (a.durStep != b.durStep) {
		return a.durStep < b.durStep;
	} else {
		return a.nodeIndex < b.nodeIndex;
	}
}

string inpFile, normalPressureFile, burstPressureFile, normalPressureMeanFile, burstPressureSimularityFile, burstEventFile, burstLabelFile, burstPressureBeforeFile, idIndexMapFile;
double normaldemandNoise, burstdemandNoise, pressureNoise, pressureSensorPercent, pmin, preq, pexp;
int normalDays, burstDays, minBurstDemand, maxBurstDemand, burstEventNum, nodeCount, linkCount, stepsPerDay, maxBurstDuration, minBurstDuration;
long timeStep, duration;
vector<string> nodeIDs;
vector<int> nodeIndexes, burstLabel;
vector<int> burstSelectedNodeIndexes;
vector<vector<double>> normalPressure, normalPressureMeans;
vector<vector<double>> burstPressure;
vector<double> initBaseDemand, burstPressureSimularity;
vector<BurstEvent> burstEvents;
set<int> burstDaySet, reservoirIndex, tankIndex;
map<int, string> nodeIndexToID;
map<string, int> nodeIDToIndex;

void read_JsonConfig() {
	std::ifstream fin("resource/config.json");
	nlohmann::json j;
	fin >> j;
	fin.close();
	inpFile = j["inpFile"];
	burstdemandNoise = j["burstdemandNoise"];
	normaldemandNoise = j["normaldemandNoise"];
	pressureNoise = j["pressureNoise"];
	normalDays = j["normalDays"];
	burstDays = j["burstDays"];
	minBurstDemand = j["minBurstDemand"];
	maxBurstDemand = j["maxBurstDemand"];
	burstEventNum = j["burstEventNum"];
	pressureSensorPercent = j["pressureSensorPercent"];
	pmin = j["pmin"];
	preq = j["preq"];
	pexp = j["pexp"];
	nodeIDs = j["nodeIDs"].get<vector<string>>();
	normalPressureFile = j["normalPressureFile"];
	burstPressureFile = j["burstPressureFile"];
	maxBurstDuration = j["maxBurstDuration"];
	minBurstDuration = j["minBurstDuration"];
	normalPressureMeanFile = j["normalPressureMeanFile"];
	burstPressureSimularityFile = j["burstPressureSimularityFile"];
	burstEventFile = j["burstEventFile"];
	burstLabelFile = j["burstLabelFile"];
	burstPressureBeforeFile = j["burstPressureBeforeFile"];
	idIndexMapFile = j["idIndexMapFile"];
	burstSelectedNodeIndexes = j["burstSelectedNodeIndexes"].get<vector<int>>();
}
//¶ÁÈ¡Êý¾Ý
void read_InpData() {  
	EN_Project ph;
	int errCode;
	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0) {
		cout << "open inp file err, code =" << errCode << ", file name is " << inpFile << endl;
		return; 
	}
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);
	int index;
	for (int i = 0; i < nodeIDs.size(); i++) {
		errCode = EN_getnodeindex(ph, (char*)nodeIDs[i].c_str(), &index);
		if (errCode != 0) {
			cout << "EN_getnodeindex, code=" << errCode << endl;
			return;
		}
		cout << "node id " << nodeIDs[i] << " ======> node index " << index << endl;
		nodeIndexes.push_back(index); 
	}
	double value;
	char id[100];
	int nodeType;
	ofstream out(idIndexMapFile);
	if (!out) {
		cout << "open file failed, fileName = " << idIndexMapFile << endl;
		exit(1);
	}
	for (int i = 1; i <= nodeCount; i++) {  
		EN_getnodevalue(ph, i, EN_BASEDEMAND, &value);
		EN_getnodeid(ph, i, id); 
		EN_getnodetype(ph, i, &nodeType);
		nodeIndexToID[i] = id;
		nodeIDToIndex[id] = i;
		initBaseDemand.push_back(value);
		if (nodeType == EN_RESERVOIR) {
			reservoirIndex.insert(i);
			cout << "node index " << i << " (id "<< id << ")is reservoir" << endl;  
		} else if (nodeType == EN_TANK) {
			tankIndex.insert(i);
			cout << "node index " << i << " (id " << id << ")is tank" << endl;
		}
		out << "node index " << i << " ======> node id " << id << endl;

	}
	out.close();
	EN_gettimeparam(ph, EN_HYDSTEP, &timeStep);
	EN_gettimeparam(ph, EN_DURATION, &duration);
	stepsPerDay = duration / timeStep + 1;
	errCode = EN_deleteproject(ph); 
	if (errCode != 0) {
		cout << "EN_deleteproject, code=" << errCode << endl;
		return;
	}
}

void generate_BurstEvent() { 
	srand(324);
	for (int i = 0; i < burstEventNum * 100000000; i++) {
		BurstEvent be = {
		//	277, 
		//	"",
		//	40, 
		//	5,
		//	rand() % (maxBurstDemand - minBurstDemand + 1) + minBurstDemand,
		//	0, 
		//}; 
		
		    burstSelectedNodeIndexes[rand() % burstSelectedNodeIndexes.size()],  
			"",
			rand() % (stepsPerDay - maxBurstDuration),
			rand() % (maxBurstDuration - minBurstDuration + 1) + minBurstDuration,
			rand() % (maxBurstDemand - minBurstDemand + 1) + minBurstDemand,
			rand() % burstDays,
		}; 
		be.nodeID = nodeIndexToID[be.nodeIndex];
		if (burstDaySet.find(be.burstDay) != burstDaySet.end()) {
			continue;
		}
		if (reservoirIndex.find(be.nodeIndex) != reservoirIndex.end()) {
			continue;
		}
		if (tankIndex.find(be.nodeIndex) != tankIndex.end()) {
			continue;
		}
		burstEvents.push_back(be);
		burstDaySet.insert(be.burstDay);
		if (burstEvents.size() == burstEventNum) {
			break;
		}
	}
	sort(burstEvents.begin(), burstEvents.end());
	ofstream out(burstEventFile);
	if (!out) {
		cout << "open file failed, fileName = " << burstEventFile << endl;
		exit(1);
	}
	for (int i = 0; i < burstEvents.size(); i++) {
		BurstEvent be = burstEvents[i];
		out << be.burstDay << "," << be.startStep << "," << be.durStep << "," << be.nodeID << "," << be.burstDemand << endl;
	}
	out.close();
	cout << "write burst events success£¬ file name " << burstPressureFile << endl;
}

void print_ImportantInformation() {
	cout << "inpFile:\t" << inpFile << endl;
	cout << "normalPressureFile:\t" << normalPressureFile << endl;
	cout << "burstPressureFile:\t" << burstPressureFile << endl;
	cout << "burstdemandNoise:\t" << burstdemandNoise << endl;
	cout << "normaldemandNoise:\t" << normaldemandNoise << endl;
	cout << "pressureNoise:\t" << pressureNoise << endl;
	cout << "pressureSensorPercent:\t" << pressureSensorPercent << endl;
	cout << "pmin:\t" << pmin << endl;
	cout << "preq:\t" << preq << endl;
	cout << "pexp:\t" << pexp << endl;
	cout << "normalDays:\t" << normalDays << endl;
	cout << "burstDays:\t" << burstDays << endl;
	cout << "minBurstDemand:\t" << minBurstDemand << endl;
	cout << "maxBurstDemand:\t" << maxBurstDemand << endl;
	cout << "burstEventNum:\t" << burstEventNum << endl;
	cout << "nodeCount:\t" << nodeCount << endl;
	cout << "linkCount:\t" << linkCount << endl;
	cout << "stepsPerDay:\t" << stepsPerDay << endl;
	cout << "timeStep:\t" << timeStep << endl;
	cout << "duration:\t" << duration << endl;
	cout << "maxBurstDuration:\t" << maxBurstDuration << endl;
	cout << "minBurstDuration:\t" << minBurstDuration << endl;

	cout << "nodeIDs:\t[ ";
	for (int i = 0; i < nodeIDs.size(); i++) {
		cout << nodeIDs[i] << " ";
	}
	cout << "]" << endl; 

	cout << "nodeIndexes:\t[ ";
	for (int i = 0; i < nodeIndexes.size(); i++) {
		cout << nodeIndexes[i] << " ";
	}
	cout << "]" << endl;

	cout << "burstEvents:" << endl;
	for (int i = 0; i < burstEvents.size(); i++) {
		cout << "{ burstDay:" << burstEvents[i].burstDay << ",\t\tstartStep: " << burstEvents[i].startStep << ",\t\tdurStep: " << burstEvents[i].durStep << ",\t\tnodeID: " << burstEvents[i].nodeID << ",\t\tburstDemand: " << burstEvents[i].burstDemand << " }" << endl;
	}
	cout << endl;
}

vector<double> compute_NoisedVector(vector<double> vec, double u, double sigma) {
	random_device rd;
	mt19937 gen(324);
	normal_distribution<double> normal(u, sigma);
	vector<double> result;
	for (int i = 0; i < vec.size(); i++) {
		result.push_back(vec[i] * (1 + normal(gen)));
	}
	return result; 
}

vector<vector<double>> compute_NoisedMatrix(vector<vector<double>> mat, double u, double sigma) {
	random_device rd;
	mt19937 gen(324);
	normal_distribution<double> normal(u, sigma);
	vector<vector<double>> result;
	for (int i = 0; i < mat.size(); i++) {
		vector<double> tmp;
		for (int j = 0; j < mat[i].size(); j++) {
			
			tmp.push_back(mat[i][j] + normal(gen));
		}
		result.push_back(tmp);
	}
	return result;
}

void run_NormalPressure() {
	EN_Project ph;
	int errCode;

	
	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
		cout << "open inp file err, code =" << errCode << ", file name is " << inpFile << endl;
		return;
	}


	errCode = EN_setdemandmodel(ph, EN_PDA, pmin, preq, pexp);
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}
	for (int i = 0; i < normalDays; i++) {
		long t, tstep;
		errCode = EN_openH(ph); 
		if (errCode != 0) {
			cout << "EN_openH, code=" << errCode << endl;
			return;
		}
		vector<double> newBaseDemand = compute_NoisedVector(initBaseDemand, 0, normaldemandNoise);
		for (int k = 1; k <= nodeCount; k++) {
			errCode = EN_setnodevalue(ph, k, EN_BASEDEMAND, newBaseDemand[k - 1]);
			if (errCode != 0) {
				cout << "EN_setnodevalue, code=" << errCode << endl;
				return;
			}
		}
		errCode = EN_initH(ph, EN_NOSAVE); 
		if (errCode != 0) {
			cout << "EN_initH, code=" << errCode << endl;
			return;
		}
		do {
			errCode = EN_runH(ph, &t); 
			if (errCode != 0) {
				cout << "EN_runH, code=" << errCode << endl;
				return;
			}
			if (t % timeStep == 0) {
				vector<double> pressure;
				double value;
				for (int j = 0; j < nodeIndexes.size(); j++) {
					errCode = EN_getnodevalue(ph, nodeIndexes[j], EN_PRESSURE, &value);
					if (errCode != 0) {
						cout << "EN_getnodevalue, code=" << errCode << endl;
						return;
					}
					if (isnan(value)) {
						cout << "perssure value is nan, please rerun" << endl;
						return;
					}
					pressure.push_back(value);
				}
				normalPressure.push_back(pressure);
			}
			errCode = EN_nextH(ph, &tstep);
			if (errCode != 0) {
				cout << "EN_nextH, code=" << errCode << endl;
				return;
			}
		} while (tstep > 0);

		errCode = EN_closeH(ph);
		if (errCode != 0) {
			cout << "EN_closeH, code=" << errCode << endl;
			return;
		}
	}
	errCode = EN_deleteproject(ph);
	if (errCode != 0) {
		cout << "EN_deleteproject, code=" << errCode << endl;
		return;
	}

	normalPressure = compute_NoisedMatrix(normalPressure, 0, pressureNoise); 
	ofstream out(normalPressureFile);
	if (!out) {
		cout << "open file failed, fileName = " << normalPressureFile << endl;
		exit(1);
	}
	for (int i = 0; i < normalPressure.size(); i++) {
		for (int j = 0; j < normalPressure[i].size(); j++) {
			if (j != 0) {
				out << ",";
			}
			out << normalPressure[i][j];
		}
		out << endl;
	}
	out.close();
	cout << "write normal pressure success£¬ file name " << burstPressureFile << endl;
	normalPressure.clear();
	normalPressure.shrink_to_fit();
}

void run_BurstPressure() {
	EN_Project ph;
	int errCode;

	
	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
		cout << "open inp file err, code =" << errCode << ", file name is " << inpFile << endl;
		return;
	}

	
	errCode = EN_setdemandmodel(ph, EN_PDA, pmin, preq, pexp);
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}
	map<int, vector<double>> noisedDemand;
	for (int i = 0; i < burstDays; i++) {
		long t, tstep;
		errCode = EN_openH(ph);
		if (errCode != 0) {
			cout << "EN_openH, code=" << errCode << endl;
			return;
		}
		vector<double> newBaseDemand = compute_NoisedVector(initBaseDemand, 0, burstdemandNoise); 
		if (burstDaySet.find(i) != burstDaySet.end()) {
			noisedDemand[i] = newBaseDemand;
		}
		for (int k = 1; k <= nodeCount; k++) {
			errCode = EN_setnodevalue(ph, k, EN_BASEDEMAND, newBaseDemand[k - 1]);
			if (errCode != 0) {
				cout << "EN_setnodevalue, code=" << errCode << endl;
				return;
			}
		}
		errCode = EN_initH(ph, EN_NOSAVE);
		if (errCode != 0) {
			cout << "EN_initH, code=" << errCode << endl;
			return;
		}
		do {
			errCode = EN_runH(ph, &t); 
			if (errCode != 0) {
				cout << "EN_runH, code=" << errCode << endl;
				return;
			}
			if (t % timeStep == 0) {  
				vector<double> pressure;
				double value;
				for (int j = 0; j < nodeIndexes.size(); j++) {
					errCode = EN_getnodevalue(ph, nodeIndexes[j], EN_PRESSURE, &value);
					if (errCode != 0) {
						cout << "EN_getnodevalue, code=" << errCode << endl;
						return;
					}
					pressure.push_back(value);
				}
				burstPressure.push_back(pressure);
				burstLabel.push_back(0);
			}
			errCode = EN_nextH(ph, &tstep);
			if (errCode != 0) {
				cout << "EN_nextH, code=" << errCode << endl;
				return;
			}
		} while (tstep > 0);

		errCode = EN_closeH(ph);
		if (errCode != 0) {
			cout << "EN_closeH, code=" << errCode << endl;
			return;
		}
	}
	
	for (int i = 0; i < burstEvents.size(); i++) {
		BurstEvent burstEvent = burstEvents[i];

		long t, tstep;
		errCode = EN_openH(ph);
		if (errCode != 0) {
			cout << "EN_openH, code=" << errCode << endl;
			return;
		}
		vector<double> newBaseDemand = noisedDemand[burstEvent.burstDay];
		for (int k = 1; k <= nodeCount; k++) {
			errCode = EN_setnodevalue(ph, k, EN_BASEDEMAND, newBaseDemand[k - 1]);
			if (errCode != 0) {
				cout << "EN_setnodevalue, code=" << errCode << endl;
				return;
			}
		}
		errCode = EN_setnodevalue(ph, burstEvent.nodeIndex, EN_BASEDEMAND, newBaseDemand[burstEvent.nodeIndex - 1] + burstEvent.burstDemand);
		if (errCode != 0) {
			cout << "EN_setnodevalue, code=" << errCode << endl;
			return;
		}
		errCode = EN_initH(ph, EN_NOSAVE);
		if (errCode != 0) {
			cout << "EN_initH, code=" << errCode << endl;
			return;
		}
		do {
			errCode = EN_runH(ph, &t);
			if (errCode != 0) {
				cout << "EN_runH, code=" << errCode << endl;
				return;
			}
			if (t % timeStep == 0 && t / timeStep >= burstEvent.startStep && t / timeStep < burstEvent.startStep + burstEvent.durStep) {
				vector<double> pressure;
				double value;
				int index = burstEvent.burstDay * stepsPerDay + t / timeStep; 
				for (int j = 0; j < nodeIndexes.size(); j++) {
					errCode = EN_getnodevalue(ph, nodeIndexes[j], EN_PRESSURE, &value);
					if (errCode != 0) {
						cout << "EN_getnodevalue, code=" << errCode << endl;
						return;
					}
					burstPressure[index][j] = value;
				}
				burstLabel[index] = 1;
			}
			errCode = EN_nextH(ph, &tstep);
			if (errCode != 0) {
				cout << "EN_nextH, code=" << errCode << endl;
				return;
			}
		} while (tstep > 0);

		errCode = EN_closeH(ph);
		if (errCode != 0) {
			cout << "EN_closeH, code=" << errCode << endl;
			return;
		}
	}
	errCode = EN_deleteproject(ph);
	if (errCode != 0) {
		cout << "EN_deleteproject, code=" << errCode << endl;
		return;
	}
	
	burstPressure = compute_NoisedMatrix(burstPressure, 0, pressureNoise); 
	ofstream out(burstPressureFile);
	if (!out) {
		cout << "open file failed, fileName = " << burstPressureFile << endl;
		exit(1);
	}
	for (int i = 0; i < burstPressure.size(); i++) {
		for (int j = 0; j < burstPressure[i].size(); j++) {
			if (j != 0) {
				out << ",";
			}
			out << burstPressure[i][j];
		}
		out << endl;
	}
	out.close();
	cout << "write burst pressure success£¬ file name " << burstPressureFile << endl;

	ofstream out2(burstLabelFile);
	if (!out2) {
		cout << "open file failed, fileName = " << burstLabelFile << endl;
		exit(1);
	}
	for (int i = 0; i < burstLabel.size(); i++) {
		out2 << burstLabel[i] << endl;
	}
	out2.close();
	cout << "write burst label success£¬ file name " << burstLabelFile << endl;
}

void compute_NormalPressureMean() {
	for (int i = 0; i < stepsPerDay; i++) {
		vector<double> tmpResult;
		for (int j = 0; j < nodeIDs.size(); j++) {
			vector<double> data;
			for (int k = 0; k < normalDays; k++) {
				data.push_back(normalPressure[i * k][j]);
			}
			double mean = accumulate(begin(data), end(data), 0.0) / data.size();
			tmpResult.push_back(mean);
		}
		normalPressureMeans.push_back(tmpResult);
	}
	ofstream out(normalPressureMeanFile);
	if (!out) {
		cout << "open file failed, fileName = " << normalPressureMeanFile << endl;
		exit(1);
	}
	for (int i = 0; i < normalPressureMeans.size(); i++) {
		for (int j = 0; j < normalPressureMeans[i].size(); j++) {
			if (j != 0) {
				out << ",";
			}
			out << normalPressureMeans[i][j];
		}
		out << endl;
	}
	out.close();
	cout << "write normal pressure mean success£¬ file name " << normalPressureMeanFile << endl;
}

double compute_CosineSimilarity(vector<double> v1, vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "compute_CosineSimilarity error" << endl;
		return 0;
	}
	double accum1 = 0.0;
	for_each(begin(v1), end(v1), [&](const double d) {
		accum1 += d * d;
	});
	double accum2 = 0.0;
	for_each(begin(v2), end(v2), [&](const double d) {
		accum2 += d * d;
	});
	double accum3 = 0.0;
	for (int i = 0; i < v1.size(); i++) {
		accum3 += v1[i] * v2[i];
	}
	return accum3 / (sqrt(accum1) * sqrt(accum2));
}

double compute_PearsonSimilarity(vector<double> v1, vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "compute_PearsonSimilarity error" << endl;
		return 0;
	}
	double accum1 = 0.0, accum2 = 0.0, accum3 = 0.0, accum4 = 0.0, accum5 = 0.0;
	for (int i = 0; i < v1.size(); i++) {
		accum1 += v1[i] * v2[i];
		accum2 += v1[i];
		accum3 += v2[i];
		accum4 += v1[i] * v1[i];
		accum5 += v2[i] * v2[i];
	}
	return (v1.size() * accum1 - accum2 * accum3) / sqrt((v1.size() * accum4 - accum2 * accum2) * (v1.size() * accum5 - accum3 * accum3));
}

double compute_EuclideanDistance(vector<double> v1, vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "computeEuclideanDistance error" << endl;
		return 0;
	}
	double accum = 0.0;
	for (int i = 0; i < v1.size(); i++) {
		accum += (v1[i] - v2[i])* (v1[i] - v2[i]);
	}
	return sqrt(accum);
}

double compute_ManhattanDistance(vector<double> v1, vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "compute_ManhattanDistance error" << endl;
		return 0;
	}
	double accum = 0.0;
	for (int i = 0; i < v1.size(); i++) {
		accum += abs(v1[i] - v2[i]);
	}
	return accum;
}

void compute_BurstPressureSimilarity() {
	for (int i = 0; i < burstPressure.size(); i++) {
		double simularity = compute_EuclideanDistance(burstPressure[i], normalPressureMeans[i % stepsPerDay]); 
		burstPressureSimularity.push_back(simularity);
	}
	ofstream out(burstPressureSimularityFile);
	if (!out) {
		cout << "open file failed, fileName = " << burstPressureSimularityFile << endl;
		exit(1);
	}
	for (int i = 0; i < burstPressureSimularity.size(); i++) {
		out << burstPressureSimularity[i] << endl;
	}
	out.close();
	cout << "write burst pressure simularity success£¬ file name " << burstPressureSimularityFile << endl;
}

int main() {
	_setmaxstdio(8096);
	read_JsonConfig();
	read_InpData();
	generate_BurstEvent();
	print_ImportantInformation();

	//run_NormalPressure();
	run_BurstPressure();
} 
