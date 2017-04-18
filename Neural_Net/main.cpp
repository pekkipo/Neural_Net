//
//  main.cpp
//  Neural_Net
//
//  Created by Aleksei Petukhov on 16/04/2017.
//  Copyright © 2017 PekkiPo. All rights reserved.
//


/* Structure
 
Neural_Net class interface:
 - feedForward()
 - backProb()
 - getResults()

 
 
 */

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

using namespace std;

class Neuron;
typedef  vector<Neuron> Layer;


// connection from a neuron
struct Connection {
    double weight;
    double deltaWeight;
};

// ******* Neuron class ********
class Neuron {
public:
    Neuron(unsigned numOutputs, unsigned myIndex); // need to tell a Neuron how many neurons in a next layer
    void feedForward(const Layer &prevLayer);
    void setOutputVal(double val) { m_outputVal = val; };
    double getOutputVal(void) const { return m_outputVal; };
    
private:
    static double randomWeight(void) { return rand() / double(RAND_MAX);}
    double m_outputVal;
    vector<Connection> m_outputWeights;
    unsigned m_myIndex;
    static double transferFunction(double x);
    static double transferFunctionDerivative(double x);
    
    
};

Neuron::Neuron(unsigned numOutputs, unsigned myIndex) {
    for (unsigned c = 0; c < numOutputs; ++c) {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight(); // work with just added element
        
    }
    
    m_myIndex = myIndex;
}

void Neuron::feedForward(const Layer &prevLayer) {
    
    double sum = 0.0;
    
    // Sum the previous layer's outputs (which are actually our inputs)
    // Inlcude the bias node from the previous layer
    
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
        // we sum up inputs with their respective connection weights
    }
    
    m_outputVal = Neuron::transferFunction(sum);
    // tr - is an activation function
}

// Transfer function
double Neuron::transferFunction(double x) {
    
    // we ll use a hyperbolic tangent function
    // output range [-1.0...1.0]
    return tanh(x);
    
}

double Neuron::transferFunctionDerivative(double x) {
    
    return 1.0 - x*x;
}

// ****** Neural Net class *******
class Neural_Net {

public:
    Neural_Net(const vector<unsigned> &topology);
    
    void feedForward(const vector<double> &inputVals); // & pass by reference. conts - because
    //this function doesn't change anything in input vals/ We promise not to change it
    void backProp(const vector<double> &targetVals) {};
    void getResults(vector<double> &resultsVals) const {};
    // const correctness here - because member function doesn't modify any information
    
    
private:
    vector<Layer> m_layers; // m_layers[layer_Number][neuron_Num_in_Layer]
    // each layer contain Neurons
    
};


Neural_Net::Neural_Net(const vector<unsigned> &topology) {

    // we need to make a bunch of neurons out of input container
    unsigned numLayers = topology.size();
    
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
    // create a layer and add it to m_layer
        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum +1];
        // may be different (input, output layers have different number of neurons compared to hidden)
        // first condition - if output layer
        
        
        // The layer is made and now one needs to fiil it with neurons and add a bias
        // neuron to the layer
        for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum) {
            // <= to add a bias neuron
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
            // the last element in the container
            cout << "Made a neuron" << endl;
            
        }
        
    }
    
}

void Neural_Net::feedForward(const vector<double> &inputVals) {
    
    assert(inputVals.size() == m_layers[0].size() - 1); // 0th in m_layers is an input layer
    // -1 because we also have additional bias neuron
    // assert. In production do error checking instead
    
    // goes through every imput value
    for(unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);
        // 0 - input layer, i - neuron in the input layer
    }
    
    // Forward propagate
    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        // from one as we skipping the input layer

        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
            // layer number - neuron number
        }
    }
}

int main() {
    
    vector<unsigned> topology;
    topology.push_back(3);
    topology.push_back(2);
    topology.push_back(1);
    Neural_Net net(topology);
    
    vector<double> inputVals; // vector can be variable lenght array
    vector<double> targetVals;
    vector<double> resultVals;
    // Training
    net.feedForward(inputVals);
    net.backProp(targetVals);
    
    
    net.getResults(resultVals);
    
}
