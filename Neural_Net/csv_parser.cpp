//
//  main.cpp
//  Task2_3
//
//  Created by Aleksei Petukhov on 05/04/2017.
//  Copyright © 2017 PekkiPo. All rights reserved.


// I assumed the values in the data to be ints!


#include <istream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class CSVparser {
 
public:
    

    
    void ReadFile(string);
    void WriteToFile(string);
    void CheckValue(string&, int, int);
    string InterpolateValue(int, int);
    
    // Default values for my convenience. Files are in the working directory
    string default_input;// = "data.csv";
    string default_output;// = "output.csv";
    
    // Сtor
    CSVparser(void) : default_input("data.csv"), default_output("output.csv") {}; // member initializer list
    
protected:
    vector<vector<string> > matrix;
    vector<string> row;
    string line;
    string cell;
 };

void CSVparser::CheckValue(string &value, int i, int j) {
    
    if (value == "0") {
        value = InterpolateValue(i, j);
    }
}

string CSVparser::InterpolateValue(int i, int j) {
    
    /* Bilinear interpolation scheme
     y2 --   Q12--------Q22
     |       |     |     |
     y  --   ------P -----
     |       |     |     |
     y1 --   Q11---------Q21
             |     |     |
             x1    x     x2
     */
    
    
    int Q12;
    int Q22;
    int Q11;
    int Q21;
    
    //If bad value is the most left/right or top/bottom - consider lacking surrounding values to be zeros
    
    // First check if this is the first line. Q12 and Q22 might be unavailable
    if (i == 0) {
        // take the value on the last line same column
        Q12 = 0;
        Q22 = 0;
    } else {
        Q12 = atoi(matrix[i-1][j-1].c_str());
        Q22 = atoi(matrix[i-1][j+1].c_str());
    }
    
    // Check if this is the last line. Q11 and Q21 might be unavailable
    if (i == matrix.size()-1) {
        // take the value on the first line same column
        Q11 = 0;
        Q21 = 0;
    } else {
        Q11 = atoi(matrix[i+1][j-1].c_str());
        Q21 = atoi(matrix[i+1][j+1].c_str());
    }
    
    
    // Mirroring method. I think that's an overkill for that task and makes no sense
    // Bad value - most left
    // change it to the most right of that row
    /*
    if (Q12 == 0) {
        Q12 = atoi(matrix[i-1][int(matrix[i].size()) - 1].c_str());
    }
    
    if (Q11 == 0) {
        Q11 = atoi(matrix[i+1][int(matrix[i].size()) - 1].c_str());
    }
    
    // Bad value - most right
    // change it to the most left of that row
    if (Q22 == 0) {
        Q22 = atoi(matrix[i-1][int(matrix[i].size()) - 1].c_str());
    }
    
    if (Q21 == 0) {
        Q21 = atoi(matrix[i+1][int(matrix[i].size()) - 1].c_str());
    }
    */
    
    // Linear interpolation in one dimension
    int f_x_y1 = (((j+1.0) - j) / ((j+1.0) - (j-1.0)))*Q11 + ((j - (j-1.0)) / ((j+1.0) - (j-1.0)))*Q21;
    int f_x_y2 = (((j+1.0) - j) / ((j+1.0) - (j-1.0)))*Q12 + ((j - (j-1.0)) / ((j+1.0) - (j-1.0)))*Q22;
    // Combine
    int result = (((i-1.0) - i) / ((i-1.0) - (i+1.0)))*f_x_y1 + ((i - (i+1.0)) / ((i-1.0) - (i+1.0)))*f_x_y2;
    
    return to_string(result);
    
}

void CSVparser::ReadFile(string filename){
    
    /*
    if (filename == "") {
        filename = default_input; // in case I just want to read the test file in the working folder
    }
    */
    
    ifstream file(filename);
    
    if(!file.is_open()) cout << "ERROR: Cannot open the file. Check the path" << '\n';
    
    while(file)
    {
        getline(file,line);
        stringstream lineStream(line);
        row.clear();
        
        while(getline(lineStream, cell, ','))
            row.push_back(cell);
        
        if(!row.empty())
            matrix.push_back( row );
    }
    
}

void CSVparser::WriteToFile(string filename) {
    
    /*
    if (filename == "") {
        filename = default_output; // in case I just want to read the test file in the working folder
    }
    */
        ofstream output(filename);
        
        for( int i=0; i<int(matrix.size()); i++ )
        {
            for( int j=0; j<int(matrix[i].size()); j++ )
            {
                
                CheckValue(matrix[i][j], i, j);
                cout << matrix[i][j] << " "; // just for displaying
                
                if (j != int(matrix[i].size())-1)
                {
                    output << matrix[i][j] << ',';
                }
                else
                {
                    output << matrix[i][j]; // if the last column, no need for comma
                }

            }
            
            cout << endl;
            output << endl;
            
        }

}


int main(int argc,char *argv[]) {
    
    if (argc != 3) {
    cout << "Not enough input arguments (must be 2: input file and output file" << endl;
    return -1;
    }
    
    string input_file = argv[1];
    string output_file = argv[2];
    
    CSVparser parser;
    parser.ReadFile(input_file);
    parser.WriteToFile(output_file);
    
    return 0;
    
    // Terminal
    // Compile: g++ main.cpp
    // Run: ./a.out file_to_read output_file
}





