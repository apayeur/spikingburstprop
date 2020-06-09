//
//  main.cpp
//  burstprop_mnist
//
//  Created by Alexandre Payeur on 3/16/18.
//  Copyright © 2018 Alexandre Payeur. All rights reserved.
//

#include <iostream>
#include "readMNIST.hpp"
#include <armadillo>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;


const int nb_input = 28*28;
const int nb_hidden = 200;
const int nb_classes = 10;

const int number_of_training_examples = 60000;
const int number_of_test_images = 10000;

const bool weight_transport = false;

const double T_prediction = 5.; //was 19
const double T_teacher = 0.1*T_prediction; // was 1
const double dt = 0.1;
const int nb_epochs = 25;

const double beta = 1.;
const double lamda = 0.05;
const double tau_u = 0.1;
const double tau_z = 0.1;
const double tau_e = tau_z;
const double tau_W = 100.;  // was 30
const double tau_moving_avg = 5.;
const double e_max = 5.;

const string PATH_TO_MNIST_DATA = "./data/";
const string PATH_TO_MNIST_RESULTS = "./results/";


void initialize_weights(Mat<double> &W1, Mat<double> &W2, Mat<double> &B);
colvec softmax(colvec &x);
colvec clip(colvec &a, double a_min, double a_max);
double loss(colvec &e, colvec &target, double epsilon=1.e-7);
void feedforward(colvec & e_input, colvec & e_hidden, colvec & e_output, colvec & z_hidden, colvec & z_output, Mat<double> &W_onto_hidden, Mat<double> &W_onto_output, colvec &n_hidden, colvec &n_output);
Mat<double> to_categorical(colvec &t);
double compute_test_acc(Mat<double> &Y_test, vector<mat> &X_test, Mat<double> &W1, Mat<double> &W2, colvec &bias1, colvec &bias2);
void compute_training_acc(Mat<double> &Y_train, vector<mat> &X_train, Mat<double> &W1, Mat<double> &W2, colvec &bias1, colvec &bias2);
void construct_heterosyn_matrix(Mat<double> &H, colvec &e_pre, colvec &e_post);
colvec exponential(colvec &x){ return exp(x); }
colvec sigmoid(colvec &x){
    colvec x_clipped = clip(x, -20., 20.);
    return 1.0/(1.0+exp(-x_clipped));
}
colvec h_for_sigmoid(colvec &x){
    return ones<colvec>(x.n_elem) - x;
}
colvec h_for_exp(colvec &x){
    return ones<colvec>(x.n_elem);
}


colvec (*activation_hidden)(colvec &) = &exponential;

int main(int argc, const char * argv[]) {
    //Open files
    ofstream file1(PATH_TO_MNIST_RESULTS + "comp_PSP.dat");
    ofstream file2(PATH_TO_MNIST_RESULTS + "maxWeights.dat");
    ofstream file_e_input(PATH_TO_MNIST_RESULTS + "e_input.dat");
    ofstream file_e_hidden(PATH_TO_MNIST_RESULTS + "e_hidden.dat");
    ofstream file_e_output(PATH_TO_MNIST_RESULTS + "e_output.dat");
    ofstream file_output_error(PATH_TO_MNIST_RESULTS + "output_error.dat");
    ofstream file_hidden_error(PATH_TO_MNIST_RESULTS + "hidden_error.dat");
    ofstream file_teacher(PATH_TO_MNIST_RESULTS + "teacher.dat");
    ofstream file_delta_pred(PATH_TO_MNIST_RESULTS + "delta_prediction.dat");
    ofstream file_testerror(PATH_TO_MNIST_RESULTS + "testerror.dat");
    ofstream file_deltap_squared(PATH_TO_MNIST_RESULTS + "deltap_square.dat");
    ofstream file_p_output(PATH_TO_MNIST_RESULTS + "p_output.dat");
    ofstream file_p_hidden(PATH_TO_MNIST_RESULTS + "p_hidden.dat");
    ofstream file_p_bar_output(PATH_TO_MNIST_RESULTS + "p_bar_output.dat");
    ofstream file_p_bar_hidden(PATH_TO_MNIST_RESULTS + "p_bar_hidden.dat");
    ofstream file_w(PATH_TO_MNIST_RESULTS + "example_weights.dat");

    // Read training data
    string filenameTrainingImages = PATH_TO_MNIST_DATA +"train-images-idx3-ubyte";
    vector<mat> x_train;
    read_Mnist(filenameTrainingImages, x_train);
    string filenameTrainLabels = PATH_TO_MNIST_DATA + "train-labels-idx1-ubyte";
    colvec y_train = zeros<colvec>(number_of_training_examples);
    read_Mnist_Label(filenameTrainLabels, y_train);
    
    // Read test data
    string filenameTestImages = PATH_TO_MNIST_DATA + "t10k-images-idx3-ubyte";
    vector<mat> x_test;
    read_Mnist(filenameTestImages, x_test);
    string filenameTestLabels = PATH_TO_MNIST_DATA + "t10k-labels-idx1-ubyte";
    colvec y_test = zeros<colvec>(number_of_test_images);
    read_Mnist_Label(filenameTestLabels, y_test);
    
    // Convert labels to one-hot representation
    Mat<double> Y_train = to_categorical(y_train);
    Mat<double> Y_test = to_categorical(y_test);

    // Initialize the random generator
    arma_rng::set_seed(1);
        
    // Weights and biases
    Mat<double> W1, W2, B, H1, H2;
    H2.zeros(nb_classes, nb_hidden); // heterosyn matrix
    H1.zeros(nb_hidden, nb_input);
    Mat<double> W2_tmp;
    colvec n1, n2;
        
    // Initialize weights and biases
    n1.zeros(nb_hidden);
    n2.zeros(nb_classes);
    initialize_weights(W1, W2, B);
    
    // Define h(e)
    colvec (*h)(colvec &);
    if (activation_hidden == &exponential) h = &h_for_exp;
    else if (activation_hidden == &sigmoid) h = &h_for_sigmoid;
    
    // Create variables
    colvec e0(nb_input, fill::zeros);
    colvec e1(nb_hidden, fill::zeros);
    colvec e1_cache(nb_hidden, fill::zeros);
    colvec e2(nb_classes, fill::zeros);
    colvec e2_cache(nb_classes, fill::zeros);
    colvec z1(nb_hidden, fill::zeros);
    colvec z2(nb_classes, fill::zeros);
    
    colvec p2(nb_classes, fill::zeros);
    colvec p2_cache(nb_classes, fill::zeros);
    colvec p2_bar(nb_classes, fill::zeros);
    colvec p2_bar_cache(nb_classes, fill::zeros);
    
    colvec e1_bar(nb_hidden, fill::zeros);
    colvec e2_bar(nb_classes, fill::zeros);
    colvec e1_bar_cache(nb_hidden, fill::zeros);
    colvec e2_bar_cache(nb_classes, fill::zeros);

    colvec b1_bar(nb_hidden, fill::zeros);
    colvec b2_bar(nb_classes, fill::zeros);
    colvec b1_bar_cache(nb_hidden, fill::zeros);
    colvec b2_bar_cache(nb_classes, fill::zeros);
    
    colvec p1(nb_hidden, fill::zeros);
    colvec p1_cache(nb_hidden, fill::zeros);
    colvec p1_bar(nb_hidden, fill::zeros);
    colvec p1_bar_cache(nb_hidden, fill::zeros);
    
    colvec prediction_before(nb_classes, fill::zeros);
    colvec prediction_after(nb_classes, fill::zeros);

    colvec small_vector(nb_classes, fill::zeros);
    small_vector.fill(1.e-8);
    
    colvec p_ref(nb_classes, fill::zeros);
    p_ref.fill(0.5);
    
    colvec deltab2(nb_classes, fill::zeros);
    colvec deltab2_cache(nb_classes, fill::zeros);
    colvec deltab1(nb_hidden, fill::zeros);
    colvec deltab1_cache(nb_hidden, fill::zeros);
    colvec b1(nb_hidden, fill::zeros);
    colvec b1_cache(nb_hidden, fill::zeros);
    colvec b2(nb_classes, fill::zeros);
    colvec b2_cache(nb_classes, fill::zeros);
    colvec u(nb_hidden, fill::zeros);
    colvec u_cache(nb_hidden, fill::zeros);
    colvec target(nb_classes, fill::zeros);
    colvec notarget(nb_classes, fill::zeros);
    colvec deltap2(nb_classes, fill::zeros);
    colvec deltap_cache(nb_hidden, fill::zeros);
    colvec hidden_error(nb_hidden, fill::zeros);
    colvec PSP;

    
    // initializations
    int m = number_of_training_examples;
    p2_bar = p_ref;
    p1 = p_ref;
    p1_bar = p_ref;

    
    for (int epoch=0;epoch<nb_epochs;epoch++){
        cout << "Epoch #" << epoch+1 << " ..." <<endl;
        double cost_epoch = 0.;
        
        for (int example(0);example<m;example++){
            
            // FEEDFORWARD PHASE
            // select appropriate example and target in training set
            e0 = vectorise(x_train[example])/255.;
            target = Y_train.col(example);

            // initialize cost for that given example
            double cost = 0.;
            
            for (int t=0;t<int(T_prediction/dt);t++){

                // record a few examples
                /*if (epoch == 0 and (example==65 or example==66 or example==67)){
                    e0.t().raw_print(file_e_input);
                    e1.t().raw_print(file_e_hidden);
                    e2.t().raw_print(file_e_output);
                    notarget.t().raw_print(file_teacher);
                    p2.t().raw_print(file_p_output);
                    p1.t().raw_print(file_p_hidden);
                    p2_bar.t().raw_print(file_p_bar_output);
                    p1_bar.t().raw_print(file_p_bar_hidden);
                    file_w << W1(100, 300) << "\t" << W2(5, 100) << std::endl;
                }*/
                
                // put variables of previous step into "caches"
                e1_cache = e1;
                e2_cache = e2;
                b1_cache = b1;
                b2_cache = b2;
                p1_cache = p1;
                p2_cache = p2;
                e1_bar_cache = e1_bar;
                e2_bar_cache = e2_bar;
                b1_bar_cache = b1_bar;
                b2_bar_cache = b2_bar;
                p1_bar_cache = p1_bar;
                p2_bar_cache = p2_bar;
                
                // construct prediction during interval [example*T, example*T + T_prediction]
                z1 = (1.-dt/tau_z)*z1 + (dt/tau_z)*(W1*e0 + n1);
                e1 = activation_hidden(z1);
                z2 = (1.-dt/tau_z)*z2 + (dt/tau_z)*(W2*e1_cache + n2);
                e2 = softmax(z2);
                //z1 = W1*e0 + n1;
                //e1 = activation_hidden(z1);
                //z2 = W2*e1 + n2;
                //e2 = softmax(z2);
                
                // burst rates during prediction interval
                p2 = p_ref;
                b2 = p2%e2;
                u = u + (dt/tau_u)*(-u + beta*h(e1_cache)%(B*b2_cache));
                //u = beta*h(e1)%(B*b2);
                p1 = sigmoid(u);
                b1 = p1%e1;
                
                // compute moving average
                e2_bar = e2_bar + (dt/tau_moving_avg)*(e2_cache-e2_bar);
                e1_bar = e1_bar + (dt/tau_moving_avg)*(e1_cache-e1_bar);
                b2_bar = b2_bar + (dt/tau_moving_avg)*(b2_cache-b2_bar);
                b1_bar = b1_bar + (dt/tau_moving_avg)*(b1_cache-b1_bar);
                p2_bar = b2_bar/e2_bar;
                p1_bar = b1_bar/e1_bar;
                
                // exponential moving averages (assuming convergence)
                //p2_bar = p2;
                //p1_bar = p1;
                //e1_bar = e1;
                //e2_bar = e2;
                //b1_bar = b1;
                //b2_bar = b2;
                
            }
            
            // BACKWARD/LEARNING PHASE
            // record prediction for that example before learning:
            prediction_before = e2;
            
            if (example % 1000 == 0){
                cout << "*** Example "<<example+1<<":"<<endl;
                cout << "Loss before = " << loss(e2, target) << endl;
                //compute_training_acc(Y_train, x_train, W1, W2, n1, n2);
                compute_test_acc(Y_test, x_test, W1, W2, n1, n2);
            }
            
            for (int t=int(T_prediction/dt);t<int((T_prediction+T_teacher)/dt);t++){
                
                /*
                // record a few examples
                if (epoch == 0 and (example==65 or example==66 or example==67)){
                    e0.t().raw_print(file_e_input);
                    e1.t().raw_print(file_e_hidden);
                    e2.t().raw_print(file_e_output);
                    target.t().raw_print(file_teacher);
                    p2.t().raw_print(file_p_output);
                    p1.t().raw_print(file_p_hidden);
                    p2_bar.t().raw_print(file_p_bar_output);
                    p1_bar.t().raw_print(file_p_bar_hidden);
                    file_w << W1(100, 300) << "\t" << W2(5, 100) << std::endl;
                }*/

                // put variables of previous step into "caches"
                e1_cache = e1;
                e2_cache = e2;
                b1_cache = b1;
                b2_cache = b2;
                p1_cache = p1;
                p2_cache = p2;
                e1_bar_cache = e1_bar;
                e2_bar_cache = e2_bar;
                b1_bar_cache = b1_bar;
                b2_bar_cache = b2_bar;
                p1_bar_cache = p1_bar;
                p2_bar_cache = p2_bar;
                
                // compute output BP
                //////////    ONLY WORK IF P_REL = 0.5     //////////
                //deltab2 = (target - e2_bar);
                deltap2 = p_ref%tanh((target - e2_bar)/(e2_bar + small_vector));
                p2 = p_ref + deltap2;
                //p2 = clip(p2, 0., 1.);
                
                // feedforward propagation
                z1 = (1.-dt/tau_z)*z1 + (dt/tau_z)*(W1*e0 + n1);
                e1 = activation_hidden(z1);
                z2 = (1.-dt/tau_z)*z2 + (dt/tau_z)*(W2*e1_cache + n2);
                e2 = softmax(z2);
                
                // output-layer burst rate
                b2 = p2%e2;

                // hidden-layer variables
                u = u + (dt/tau_u)*(-u + beta*h(e1_cache)%(B*b2_cache));
                p1 = sigmoid(u);
                b1 = p1%e1;
                
                // BP moving average
                e2_bar = e2_bar + (dt/tau_moving_avg)*(e2_cache-e2_bar);
                e1_bar = e1_bar + (dt/tau_moving_avg)*(e1_cache-e1_bar);
                b2_bar = b2_bar + (dt/tau_moving_avg)*(b2_cache-b2_bar);
                b1_bar = b1_bar + (dt/tau_moving_avg)*(b1_cache-b1_bar);
                p2_bar = b2_bar/e2_bar;
                p1_bar = b1_bar/e1_bar;
                
                //file_deltap_squared << sum(deltap%deltap)/nb_hidden << endl;
                
             
                //if (activation_hidden == &sigmoid){
                //     hidden_error = deltap%e1_cache%(1.-e1_cache);
                //}
                //else if (activation_hidden == &exponential){
                //    hidden_error = deltap%e1_cache;
                //}

                // construct heterosynaptic plasticity matrix
                construct_heterosyn_matrix(H1, e0, e1);
                
                // weight update
                n1 = n1 + (dt/tau_W)*((p1_cache-p1_bar_cache)%e1_cache);
                n2 = n2 + (dt/tau_W)*((p2_cache-p2_bar_cache)%e2_cache);
                W1 = W1 + (dt/tau_W)*(((p1_cache-p1_bar_cache)%e1_cache)*e0.t() - lamda*H1);
                W2 = W2 + (dt/tau_W)*(((p2_cache-p2_bar_cache)%e2_cache)*e1_cache.t());

                if (weight_transport) B = W2.t();
                
            }
            // record prediction for that example after learning
            prediction_after = e2;
            (prediction_after - prediction_before).t().raw_print(file_delta_pred);
            
            double local_loss = loss(e2, target);
            if (example % 1000 == 0){
                cout << "Loss after = " << local_loss << endl << endl;
            }
            cost = local_loss;

            //cost /= (T_teacher/dt);
            //if (example % 100 == 0){
            //    cout << "cost for example " << example+1 << " = " << cost << endl;
            //}
            cost_epoch += cost;
        }
        cout << "cost = " << cost_epoch/m << endl;
        double testacc = compute_test_acc(Y_test, x_test, W1, W2, n1, n2);
        file_testerror << epoch + 1 << "\t" << 100.*(1.-testacc) << endl;
        //compute_training_acc(Y_train, x_train, W1, W2, n1, n2);
             
    }
    
    file1.close();
    file2.close();
    file_e_hidden.close();
    file_e_output.close();
    file_output_error.close();
    file_hidden_error.close();
    file_teacher.close();
    file_delta_pred.close();
    file_testerror.close();
    file_deltap_squared.close();
    file_p_output.close();
    file_p_hidden.close();
    file_p_bar_output.close();
    file_p_bar_hidden.close();
    file_w.close();
    
    return 0;
}



void feedforward(colvec & e_input, colvec & e_hidden, colvec & e_output, colvec & z_hidden, colvec & z_output, Mat<double> &W_onto_hidden, Mat<double> &W_onto_output, colvec &n_hidden, colvec &n_output, Mat<double> &Y, colvec &p1, colvec &p1_bar, colvec &p2, colvec &p2_bar, colvec &u){
    
    // put variables of previous step into a "cache"
    colvec hidden_cache = e_hidden;
    colvec output_cache = e_output;
    colvec p1_cache = p1;

    // feedforward propagation
    z_hidden = (1.-dt/tau_z)*z_hidden + (dt/tau_z)*(W_onto_hidden*e_input + n_hidden);
    e_hidden = activation_hidden(z_hidden);
    z_output = (1.-dt/tau_z)*z_output + (dt/tau_z)*(W_onto_output*hidden_cache + n_output);
    e_output = softmax(z_output);
    
    // BP hidden
    u += (dt/tau_u)*(-u + beta*Y*(p2%output_cache));
    p1 = sigmoid(u);
    
    // BP moving average
    p2_bar += (dt/tau_moving_avg)*(p2-p2_bar);
    p1_bar += (dt/tau_moving_avg)*(p1_cache-p1_bar);
    
}

void backward(colvec & e_input, colvec & e_hidden, colvec & e_output, colvec &target, colvec &u, Mat<double> &W_onto_hidden, Mat<double> &W_onto_output, colvec &n_hidden, colvec &n_output){
    //event rates must be those of the previous time step
    /*
    colvec deltab2 = target - e_output;
    colvec PSP = beta*W_onto_output.t()*deltab2;
    u = (1.-dt/tau_b)*u + (dt/tau_b)*PSP;
    colvec deltap = sigmoid(PSP) - 0.5;
    //deltab1 = (1.-dt/tau_b)*deltab1 + (dt/tau_b)*tanh(beta*W2.t()*deltab2_cache)%e1_cache;
    //deltab1 = (1.-dt/tau_b)*deltab1 + (dt/tau_b)*tanh(beta*B.t()*deltab2)%e1%(1.-e1);
    
    
    n_hidden += -gammaWD*n_hidden + (dt/tau_W)*eta*deltap_cache%e_hidden%(1.-e_hidden);
    n_output += -gammaWD*n_output + (dt/tau_W)*eta*deltab2;
    W_onto_hidden += -gammaWD*W_onto_hidden + (dt/tau_W)*eta*(deltap_cache%e_hidden%(1.-e_hidden))*e_input.t();
    W_onto_output += -gammaWD*W_onto_output + (dt/tau_W)*eta*deltab2_cache*e_hidden.t();
     */
}


void initialize_weights(Mat<double> &W1, Mat<double> &W2, Mat<double> &B){
    W1.randn(nb_hidden, nb_input);
    W1 *= sqrt(2.0/(nb_input + nb_hidden));
    W2.randn(nb_classes, nb_hidden);
    W2 *= sqrt(2.0/(nb_classes + nb_hidden));

    if (weight_transport) {
        B = W2.t();
    }
    else {
        B.randn(nb_hidden, nb_classes);
        B *= 0.1*sqrt(2.0/(nb_classes + nb_hidden));
    }
}

void construct_heterosyn_matrix(Mat<double> &H, colvec &e_pre, colvec &e_post){
    H = ((ones<colvec>(e_post.n_elem) - e_max/(e_post+1e-8*ones<colvec>(e_post.n_elem)))%(e_post > e_max*ones<colvec>(e_post.n_elem)))*e_pre.t();
    //for(uword i=0; i<H.n_rows; ++i){
    //    H.col(i) = (ones<colvec>(e_post.n_elem) - e_max/e_post)%(e_post > e_max*ones<colvec>(e_post.n_elem));
    //}
}

colvec softmax(colvec &x){
    double max = x.max();
    colvec onesvec;
    onesvec.ones(x.n_rows);
    colvec e = exp(x-max*onesvec);
    return e/sum(e);
}


double loss(colvec &e, colvec &target, double epsilon){
    colvec output = clip(e, epsilon, 1. - epsilon);
    //uword ind = target.index_max();
    double s = sum(-target % log(output));
    //double s = -log(e(ind));
    return s;
}

colvec clip(colvec &a, double a_min, double a_max){
    uvec indices_a_greater_than_a_max = find(a > a_max);
    uvec indices_a_smaller_than_a_min = find(a < a_min);
    a.elem(indices_a_greater_than_a_max) = a_max*ones<vec>(indices_a_greater_than_a_max.n_elem);
    a.elem(indices_a_smaller_than_a_min) = a_min*ones<vec>(indices_a_smaller_than_a_min.n_elem);
    return a;
}


Mat<double> to_categorical(colvec &t){
    uword m = t.n_elem;
    Mat<double> Y(nb_classes, m, fill::zeros);
    uvec converted_t = conv_to<uvec> ::from(t);
    for (uword i=0;i<m;i++){
        Y(converted_t(i), i) = 1.;
    }
    return Y;
}

double compute_test_acc(Mat<double> &Y_test, vector<mat> &X_test, Mat<double> &W1, Mat<double> &W2, colvec &bias1, colvec &bias2){
    int running_acc(0);
    colvec e0(nb_input, fill::zeros);
    colvec e1(nb_hidden, fill::zeros);
    colvec e1_cache(nb_hidden, fill::zeros);
    colvec e2(nb_classes, fill::zeros);
    colvec e2_cache(nb_classes, fill::zeros);
    colvec z1(nb_hidden, fill::zeros);
    colvec z2(nb_classes, fill::zeros);


    for(int test_index(0);test_index<number_of_test_images;test_index++){
        
        e0 = vectorise(X_test[test_index])/255.;
        z1 = W1*e0 + bias1;
        e1 = activation_hidden(z1);
        z2 = W2*e1 + bias2;
        e2 = softmax(z2);
    
        if (Y_test.col(test_index).index_max() == e2.index_max()) running_acc +=1;
    }
 
    cout << "Test accuracy = " << double(running_acc)/number_of_test_images << endl;
    return double(running_acc)/number_of_test_images;
 }

void compute_training_acc(Mat<double> &Y_train, vector<mat> &X_train, Mat<double> &W1, Mat<double> &W2, colvec &bias1, colvec &bias2){
    int running_acc(0);
    colvec e0(nb_input, fill::zeros);
    colvec e1(nb_hidden, fill::zeros);
    colvec e1_cache(nb_hidden, fill::zeros);
    colvec e2(nb_classes, fill::zeros);
    colvec e2_cache(nb_classes, fill::zeros);
    colvec z1(nb_hidden, fill::zeros);
    colvec z2(nb_classes, fill::zeros);
    
    for(int train_index(0);train_index<number_of_training_examples;train_index++){
        
        e0 = vectorise(X_train[train_index])/255.;
        z1 = W1*e0 + bias1;
        e1 = activation_hidden(z1);
        z2 = W2*e1 + bias2;
        e2 = softmax(z2);
        
        if (Y_train.col(train_index).index_max() == e2.index_max()) running_acc +=1;
    }
    
    cout << "Train accuracy = " << double(running_acc)/number_of_training_examples << endl;
}

