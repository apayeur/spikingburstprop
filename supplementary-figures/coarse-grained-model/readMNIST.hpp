//
//  readMNIST.hpp
//  burstprop_mnist
//

#ifndef readMNIST_hpp
#define readMNIST_hpp

//#include "opencv2/core/core.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/highgui/highgui.hpp"
#include <armadillo>
#include <math.h>
#include <iostream>
#include <vector>

int ReverseInt (int i);
void read_Mnist(std::string filename, std::vector<std::vector<double> > &vec);
void read_Mnist(std::string filename, std::vector<arma::mat> &vec);
void read_Mnist_Label(std::string filename, std::vector<double> &vec);
void read_Mnist_Label(std::string filename, arma::colvec &vec);


#endif /* readMNIST_hpp */
