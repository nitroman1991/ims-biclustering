#include <math.h>
#include <omp.h>
#include <getopt.h>
#include <math.h>
#include <time.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <list>
#include <ctime>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

#include "logging.hpp"

#include "arpack++/areig.h"

using namespace std;

typedef double t_ims_real;

char *csv_separator = ",";
string input_filename = "";
bool matlab_input = false;
bool matlab_input_2 = false;

// weight of spatial edges in the incidence matrix
t_ims_real alpha = 0.5;
// edge weight between neighboring points
t_ims_real max_force = 150;
// distance inside which we introduce spatial edges
t_ims_real threshold = 2.8;
// denominator threshold
t_ims_real r = 0.25;
// how many eigenvalues to find
int num_eigens = 10;
//


int help(char *argv[]) {
    cout << "Usage: " << argv[0] << " [--mat] --input=input_file" << endl;
    return 0;
}

int main(int argc, char *argv[]) {
    const struct option longopts[] = {
        {"help",            no_argument,        0, 'h'},
        {"input",           required_argument,  0, 'i'},
        {"mat",             no_argument,        0, 'm'},
        {"mat1",            no_argument,        0, '1'},
        {"mat2",            no_argument,        0, 'n'},
        {"alpha",           required_argument,  0, 'a'},
        {"eigens",          required_argument,  0, 'e'},
        {"r",               required_argument,  0, 'r'},
        {"maxforce",        required_argument,  0, 'f'},
	    {"threshold",       required_argument,  0, 't'},
        {0,0,0,0},
    };

    if (argc < 1) return help(argv);
    timer.reset();
    int index, iarg=0;
    opterr=1;    
    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "fvh", longopts, &index);
        switch (iarg) {
            case 'h':   return help(argv);              break;
            case 'i':   input_filename = optarg;        break;
            case 'a':   alpha = atof(optarg);           break;
            case 'r':   r = atof(optarg);               break;
            case 'f':   max_force = atof(optarg);       break;
            case 'e':   num_eigens = atoi(optarg);      break;
            case 'm':   matlab_input = true;            break;
            case '1':   matlab_input = true;            break;
            case 'n':   matlab_input_2 = true;          break;
	        case 't':   threshold = atof(optarg);       break;
        }
    }

    LOG("threshold = " << threshold);
    // weight of bipartite edges in the incidence matrix
    t_ims_real beta = 1 - alpha;
	
    uint len_spectrum, num_pixels;
    t_ims_real **spectra;
    t_ims_real *maxima, *specsdiag, *pixdiag;
    t_ims_real *xcoord, *ycoord;
    	
	boost::escaped_list_separator<char> lst_separator("\\", csv_separator, "\"");
	string tmp_str = "";
	
	//reading x coordinates 
	ifstream xcoord_file(input_filename + ".x.csv");
	getline(xcoord_file, tmp_str); //in fst string - num of pixels
	num_pixels = atoi(tmp_str.c_str());
	LOG("num_pixels = " << num_pixels << "...\n");
	try{
		xcoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
	}
	catch(bad_alloc& e){
		LOG("Cannot allocate enough memory, aborting...\n");
		return 1;
	}
	getline(xcoord_file, tmp_str); //in tmp_str - string with x coords 
	xcoord_file.close();
	boost::tokenizer<boost::escaped_list_separator<char> > sep1(tmp_str, lst_separator);
	LOG("Starting read x_coord...");
	uint i = 0;
	for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator it = sep1.begin(); it != sep1.end(); ++it){
		xcoord[i] = atof((*it).c_str()); ++i;
    	}
	LOG("End read x_coord...");	

	//reading y coordinates
	ifstream ycoord_file(input_filename + ".y.csv");
	getline(ycoord_file, tmp_str); //in fst string - num of pixels
	num_pixels = atoi(tmp_str.c_str());
	try{
		ycoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
	}
	catch(bad_alloc& e){
		LOG("Cannot allocate enough memory, aborting...\n");
		return 1;
	}
	getline(ycoord_file, tmp_str); //now in tmp_str - string with y coords 
	ycoord_file.close();
	boost::tokenizer<boost::escaped_list_separator<char> > sep2(tmp_str, lst_separator);
	LOG("Starting read y_coord...");    
	i = 0;
	for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator it = sep2.begin(); it != sep2.end(); ++it){
		ycoord[i] = atof((*it).c_str()); ++i;
	}
	LOG("End read y_coord...i2 = " << i);
	
	ifstream spectra_file(input_filename + ".spectra.csv");
	getline(spectra_file, tmp_str); //взяли первую строку, там - число пикселей
	num_pixels = atoi(tmp_str.c_str());
	getline(spectra_file, tmp_str); //взяли вторую строку, там - число mz-значений
	len_spectrum = atoi(tmp_str.c_str());
	LOG("len_spectrum = " << len_spectrum << "...\n");
	try{
		spectra = allocate_2d_with_default<t_ims_real>(num_pixels, len_spectrum, 0);
		specsdiag = allocate_1d_with_default<t_ims_real>(len_spectrum, 0);
		pixdiag = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
		maxima = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
	}
	catch(bad_alloc& e){
		LOG("Cannot allocate enough memory, aborting...\n");
		return 1;
	}
	t_ims_real tmp;
	LOG("Starting read spectra...");	
	for (uint j = 0; j < len_spectrum; ++j) {
		tmp_str.clear();
		getline(spectra_file, tmp_str); //get a next string from the 
		boost::tokenizer<boost::escaped_list_separator<char> > sep3(tmp_str, lst_separator);
		i = 0;		
		for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator it = sep3.begin(); it != sep3.end(); ++it){
			tmp = atof((*it).c_str());
			spectra[i][j] = tmp;
			specsdiag[j] += tmp;
	                pixdiag[i] += tmp;
	                if (tmp > maxima[i]) maxima[i] = tmp;
			++i;
		}
	}
	spectra_file.close();
	LOG("End read spectra...i = " << i);

    

    uint k = (uint)(7*len_spectrum/(t_ims_real)8);

    // find maximal elements
    LOG("Finding maximal elements...");
    uint **U;
	try{
		U = allocate_2d_with_default<uint>(len_spectrum, num_pixels, 0);
    }
	catch(bad_alloc& e){
		LOG("Cannot allocate enough memory, aborting...\n");
		return 1;
	}
	vector<uint> indices(num_pixels);
    for (uint j=0; j<num_pixels; ++j) { indices.push_back(j); }
    for (uint i=0; i<len_spectrum; ++i) {
        sort(begin(indices), end(indices), [&](const uint & j1, const uint & j2) {
            return spectra[i][j1] > spectra[i][j2];
        });
        for (uint j=0; j<k; ++j) {
            U[i][indices[j]] = 1;
        }
        // LOG("\t" << spectra[i][indices[0]] << "\t" << spectra[i][indices[1]] << "\t" << spectra[i][indices[2]]);
    }

    LOG("Filling incidence matrix...");
    t_ims_real **W;
	t_ims_real *on_diag;
	t_ims_real *little;
	try{
		W = allocate_2d_with_default<t_ims_real>(num_pixels, num_pixels, 0);
		on_diag = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
		little = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
	}
	catch(bad_alloc& e){
		LOG("Cannot allocate enough memory, aborting...\n");
		return 1;
	}
    for (uint j1=0; j1<num_pixels; ++j1) {
        for (uint j2=0; j2<j1; ++j2) {
            t_ims_real dist = (xcoord[j1]-xcoord[j2])*(xcoord[j1]-xcoord[j2]) + (ycoord[j1]-ycoord[j2])*(ycoord[j1]-ycoord[j2]);
            if (dist <= threshold * threshold) {
                uint common_pixels = 0;
                for (uint i=0; i<len_spectrum; ++i) {
                    if (U[i][j1] == 1 && U[i][j2] == 1) {
                        ++common_pixels;
                    }
                }
                W[j1][j2] = alpha * max_force * common_pixels / ( (len_spectrum - k) * dist );
                W[j2][j1] = W[j1][j2];
                on_diag[j1] += W[j1][j2];
                on_diag[j2] += W[j1][j2];
                if (common_pixels <= r * (len_spectrum - k)) {
                    little[j1] += W[j1][j2];
                    little[j2] += W[j1][j2];
                }
            }
        }
    }

    LOG("Allocating big matrix...");
    // t_ims_real **bigL = allocate_2d_with_default<t_ims_real>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    // t_ims_real **bigW = allocate_2d_with_default<t_ims_real>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    // arma::mat bigL = arma::zeros<arma::mat>(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // Eigen::MatrixXf bigL(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigL);
    // gsl_matrix *bigW = gsl_matrix_alloc (num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigW);

    int n = num_pixels + len_spectrum;
    int nnz = num_pixels * (num_pixels + len_spectrum) + num_pixels * len_spectrum + len_spectrum;
    int *irow = new int[nnz];
    int *pcol = new int[n+1];
    t_ims_real *A = new t_ims_real[nnz];

    LOG("Filling big matrix in ARPACK's column major format...");

    uint m=0;
    pcol[0] = 0;
    t_ims_real w_inv;
    for (uint j=0; j<n; ++j) { // column by column
        w_inv = 1.0 / ( (j < num_pixels) ? (little[j] + beta * pixdiag[j]) : (beta * specsdiag[j-num_pixels]) ); // will divide by W
        if (j < num_pixels) {
            for (uint i=0; i<num_pixels; ++i) { // top left block W
                irow[m] = i;
                A[m++] = ( (i == j ? (on_diag[j] + beta * pixdiag[j]) : 0) - W[i][j] ) * w_inv;
            }
            for (uint i=0; i<len_spectrum; ++i) { // bottom left block beta*spectra
                irow[m] = num_pixels + i;
                A[m++] = -beta * spectra[j][i] * w_inv;
            }
        } else {
            for (uint i=0; i<num_pixels; ++i) { // top right block beta*transpose(spectra)
                irow[m] = i;
                A[m++] = -beta * spectra[i][j-num_pixels] * w_inv;
            }
            // bottom right block -- diagonal
            irow[m] = j;
            A[m++] = beta * specsdiag[j-num_pixels] * w_inv;
        }
        pcol[j+1] = m;
    }

    LOG("Solving for " << num_eigens << " eigenvalues with ARPACK...");

    int nconv;
    t_ims_real *EigValR = new t_ims_real[num_eigens];
    t_ims_real *EigValI = new t_ims_real[num_eigens];
    t_ims_real *EigVec = new t_ims_real[n * num_eigens];

    nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, num_eigens, "SR");
    LOG("Eigenvalues:");
    for (uint i=0; i<nconv; ++i) {
        LOG("\tlambda[" << (i+1) << "] = " << EigValR[i] << (EigValI[i] >= 0.0 ? '+' : '-') << fabs(EigValI[i]));
    }

    LOG("Printing eigenvalues to " << input_filename << ".val.csv");
    ofstream ofs_val(input_filename + ".val.csv");
    for (uint i=0; i<nconv; ++i) {
        ofs_val << EigValR[i] << ";" << EigValI[i] << "\n";
    }
    ofs_val.close();

    LOG("Printing eigenvectors to " << input_filename << ".vec.csv");
    ofstream ofs_vec(input_filename + ".vec.csv");
    for (uint i=0; i<nconv; ++i) {
        for (uint j=0; j<n; ++j) {
            ofs_vec << (i+1) << ";" << (j+1) << ";" << EigVec[i*n+j] << "\n";
        }
    }
    ofs_vec.close();

    LOG("Printing pixel coords to " << input_filename << ".coords.csv");
    ofstream ofs_coord(input_filename + ".coords.csv");
    for (uint i=0; i<num_pixels; ++i) {
        ofs_coord << xcoord[i] << ";" << ycoord[i] << "\n";
    }
    ofs_coord.close();

    delete maxima;
    delete specsdiag;
    delete pixdiag;
    delete on_diag;
    delete little;
    delete [] EigValR;
    delete [] EigValI;
    delete [] EigVec;
    delete [] irow;
    delete [] pcol;
    delete [] A;
    delete_2d<t_ims_real>(len_spectrum, spectra);
    delete_2d<t_ims_real>(num_pixels, W);

    LOG("All done.");
    return 0;
}
