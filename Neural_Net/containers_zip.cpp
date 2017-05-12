//
//  main.cpp
//  Task3
//
//  Created by Aleksei Petukhov on 07/04/2017.
//  Copyright © 2017 PekkiPo. All rights reserved.
//



/*
 
 perform_sorting is the required function
 zip, unzip - complementary functions for convenience
 
 
 */

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;


// Fill vector with pairs of elements from both vectors
template <typename A, typename B>
void zip(const vector<A> &a, const vector<B> &b, vector<pair<A,B>> &zipped)
{
    
    // Check if vectors are of the same size
    
    
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(make_pair(a[i], b[i]));
    }
}

// Break down the zipped vector into two original but sorted containers
template <typename A, typename B>
void unzip(const vector<pair<A, B>> &zipped, vector<A> &a, vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

template<typename A, typename B>
void perform_sorting(vector<A> &container_1, vector<B> &container_2) {
    
    // check vector length's
    if (container_1.size() != container_2.size()) {
        cout << "Vectors have different number of elements" << endl;
        return;
    }
    
    vector<pair<A,B>> zipped_vector;
    
    zip(container_1, container_2, zipped_vector);
    
    
    sort(begin(zipped_vector), end(zipped_vector),
         [&](const pair<A,B> &a, const pair<A,B> &b)
         {
             return a.second < b.second;
         });
    
    unzip(zipped_vector, container_1, container_2);
    
    
    // Display
    for(size_t i=0; i<container_1.size(); i++)
    {
        cout << container_1[i] << " ";
        
    }
    
    cout << endl;
    
    for(size_t i=0; i<container_2.size(); i++)
    {
        cout << container_2[i] << " ";
    }
    
    cout << endl;
    
}



int main(int argc, char* argv[])
{

    
    vector<string> c1 = {"a", "b", "c"};
    vector<int> c2 = {4,7,1};
    
    perform_sorting(c1, c2);

    return 0;
}



