//
//  neural_net.cpp
//  Neural_Net
//
//  Created by Aleksei Petukhov on 16/04/2017.
//  Copyright © 2017 PekkiPo. All rights reserved.
//

using namespace std;

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

class Neuron;
class Neural_Net;

// Define Layer type
typedef  vector<Neuron> Layer;

// Connection from a neuron
struct Connection {
    double weight;
    double deltaWeight;
};

/******* Neuron class ********/
class Neuron {
public:
    Neuron(unsigned numOutputs, unsigned myIndex); // need to tell a Neuron how many neurons in a next layer
    // Get-set
    void setOutputVal(double val) { m_outputVal = val; };
    double getOutputVal(void) const { return m_outputVal; };
    // Gradients
    void calcOutputGradients(double targetVal);
    void calcHiddenGradients(const Layer &nextLayer);
    // Update and feed forward
    void updateInputWeights(Layer &prevLayer);
    void feedForward(const Layer &prevLayer);
private:
    // Parameters
    unsigned m_myIndex;
    double m_outputVal;
    double m_gradient;
    vector<Connection> m_outputWeights;
    // Member functions
    static double transferFunction(double x);
    static double transferFunctionDerivative(double x);
    static double randomWeight(void) { return rand() / double(RAND_MAX);}
    double sunDOW(const Layer &nextLayer) const;
    // Tunable parameters
    static double eta;
    static double alpha;
};

// Adjustable parameters
double Neuron::eta = 0.15;  // 0-1 overall net learning rate
double Neuron::alpha = 0.5; //0-1 momentum. 0 - no momenutm, 0.5 - moderate momentum

// Constructor
Neuron::Neuron(unsigned numOutputs, unsigned myIndex) {
    for (unsigned c = 0; c < numOutputs; ++c) {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight(); // work with just added element
    }
    m_myIndex = myIndex;
}


void Neuron::updateInputWeights(Layer &prevLayer) {
    
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        Neuron &neuron = prevLayer[n];
        double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;
        
        double newDeltaWeight =
        eta // overall learning rate
        * neuron.getOutputVal()
        *m_gradient
        *alpha // momentum
        *oldDeltaWeight;
        
        neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
        neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
        
    }
}

double Neuron::sunDOW(const Layer &nextLayer) const {
    
    double sum = 0.0;
    
    for (unsigned n = 0; n < nextLayer.size() - 1; ++n) {
        sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
    }
    return sum;
}


void Neuron::calcOutputGradients(double targetVal) {
    double delta = targetVal - m_outputVal;
    m_gradient = delta * Neuron::transferFunctionDerivative(m_outputVal);
}

void Neuron::calcHiddenGradients(const Layer &nextLayer) {
    
    double dow = sunDOW(nextLayer); // sum of the derivatives of the next layer
    m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
}

void Neuron::feedForward(const Layer &prevLayer) {
    
    double sum = 0.0;
    
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
    }
    m_outputVal = Neuron::transferFunction(sum);
}


double Neuron::transferFunction(double x) {

    return tanh(x);
    
}

double Neuron::transferFunctionDerivative(double x) {
    
    return 1.0 - x*x;
}



/****** Neural Network class *******/
class Neural_Net {

public:
    Neural_Net(const vector<unsigned> &topology);
    
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultsVals) const;
    
private:
    vector<Layer> m_layers; // number of layers in the network. Defined by the topology
    
    double m_error;
    double m_recentAverageError;
    double m_recentAverageSmoothingFactor;
    
};


Neural_Net::Neural_Net(const vector<unsigned> &topology) {

    unsigned numLayers = topology.size();
    
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {

        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum +  1];

        for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum) {
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
            cout << "Made a neuron" << endl;
            
        }
        m_layers.back().back().setOutputVal(1.0);
    }
}

void Neural_Net::getResults(vector<double> &resultVals) const {
    
    resultVals.clear();
    
    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
        resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
}

void Neural_Net::feedForward(const vector<double> &inputVals) {
    
    assert(inputVals.size() == m_layers[0].size() - 1);
    for(unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);

    }
    
    // Forward propagate
    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
        }
    }
}

void Neural_Net::backProp(const vector<double> &targetVals) {
    
    Layer &outputLayer = m_layers.back();
    m_error = 0.0;
    
    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        double delta = targetVals[n] - outputLayer[n].getOutputVal();
        m_error += delta * delta;
    }
    m_error /= outputLayer.size() - 1; // get average error squared
    m_error = sqrt(m_error);
    
    m_recentAverageError = (m_recentAverageError * m_recentAverageSmoothingFactor + m_error) / (m_recentAverageSmoothingFactor + 1);
    
    // Calculate output layer gradients
    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        outputLayer[n].calcOutputGradients(targetVals[n]);
        
    }
    
    // Calculate gradients on hidden layers
    for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum) {
        Layer &hiddenLayer = m_layers[layerNum]; // & because just reference to the address
        Layer &nextLayer = m_layers[layerNum + 1]; // & - address, * - contents
        
        for (unsigned n = 0; n < hiddenLayer.size(); ++n) {
            hiddenLayer[n].calcHiddenGradients(nextLayer);
        }
        
    }
    
    
    // For all layers from outputs
    for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum) {
        Layer &layer = m_layers[layerNum];
        Layer &prevLayer = m_layers[layerNum - 1];
        
        for (unsigned n = 0; n < layer.size() - 1; ++n) {
            layer[n].updateInputWeights(prevLayer);
        }
    }
    
}

int main() {
    
    vector<unsigned> topology;
    topology.push_back(3); // 3 layers in the net
    topology.push_back(2); // 2 neurons in 1 hidden layer
    topology.push_back(1); // 1 output layer
    Neural_Net net(topology);
    
    vector<double> inputVals; // vector can be variable lenght array
    vector<double> targetVals;
    vector<double> resultVals;
    
    // Training
    net.feedForward(inputVals);
    net.backProp(targetVals);
    
    
    net.getResults(resultVals);
    
}
