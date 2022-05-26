#include <cmath>
#include <iostream>

void hard_inv_2d(double* A, double* Ai){
    double det = A[0]*A[3] - A[2]*A[1];
    Ai[0] = A[3]/det;
    Ai[1] = -A[1]/det;
    Ai[2] = -A[2]/det;
    Ai[3] = A[0]/det;
}

void sqmat_t(double* A, double* At, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            At[i+dim*j] = A[j + i*dim];
        }
    }
}

void sqmat_mult(double* A, double* B, double* AB, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            AB[i+dim*j] = 0;
            for (int k = 0; k < dim; k++){
                AB[i+dim*j] += A[i+dim*k] * B[k+dim*j];
            }
        }
    }
}

void eye(double* I, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            if (i == j){
                I[i+dim*j] = 1;
            }
            else{
                I[i+dim*j] = 0;
            }
        }
    }
}

double* affine_to_rot(double* A, int dim, bool hard_inv){
    double* Rot = new double[dim*dim];
    //Getting A'
    double* At = new double[dim*dim];
    sqmat_t(A, At, dim);
    
    //Getting A'*A
    double* AtA = new double[dim*dim];
    sqmat_mult(At, A, AtA, dim);

    //Schur decomposition
    double* Q = new double[dim*dim];
    double* Qt = new double[dim*dim];
    double* U = new double[dim*dim];
    double* U2 = new double[dim*dim];
    //Set U to identity
    eye(U,dim);
    double* R = new double[dim*dim];
    double* B = new double[dim*dim];

    //TODO: Change this for something a bit more sensible
    for (int i = 0; i < 999; i++){
        //Gram-schmidt method
        for (int j = 0; j < dim; j++){
            //Put column in column
            for (int i = 0; i < dim; i++){
                B[i+dim*j] = AtA[i+dim*j];
            }
            //Sum
            double* gschmidt_sum = new double[dim];
            for (int k = 0; k < j; k++){
                double multn = 0;
                double multd = 0;
                //numerator
                for (int i = 0; i < dim; i++){
                    multn += AtA[i+dim*j]*B[i+dim*k];
                }
                //denominator
                for (int i = 0; i < dim; i++){
                    multd += B[i+dim*k]*B[i+dim*k];
                }
                for (int i = 0; i < dim; i++){
                    gschmidt_sum[i] += multn/multd * B[i+dim*k];
                }
            }
            for (int i = 0; i < dim; i++){
                //The negatives are not at the same spot compared to the octave prototype!
                B[i+dim*j] -= gschmidt_sum[i];
            }
        }
        //QR factorization
        //Q
        for (int j = 0; j < dim; j++){
            double norm = 0;
            for (int i = 0; i < dim; i++){
                norm += B[i+dim*j] * B[i+dim*j];
            }
            norm = std::sqrt(norm);
            for (int i = 0; i < dim; i++){
                Q[i+dim*j] = B[i+dim*j]/norm;
            }
        }
        //R
        sqmat_t(Q, Qt, dim);
        sqmat_mult(Qt, AtA, R, dim);
        //Schur decomposition
        sqmat_mult(R, Q, AtA, dim);
        sqmat_mult(U, Q, U2, dim);
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                U[i+dim*j] = U2[i+dim*j];
            }
        }
    }
    //SQRTM*
    double* S = new double[dim*dim];
    double* Ut = new double[dim*dim];
    double* US = new double[dim*dim];
    double* X = new double[dim*dim];
    double sum = 0;
    for (int j = 0; j < dim; j++){
        S[j+dim*j] = sqrt(AtA[j+dim*j]);
        for (int i = j-1; i >= 0; i--){
            for (int k = i+1; k <= j-1; k++){
                sum += S[i+dim*k] * S[k+dim*j];
            }
            S[i+dim*j] = (AtA[i+dim*j] - sum)/(S[i+dim*i] + S[j+dim*j]);
            sum = 0;
        }
    }

    sqmat_t(U, Ut, dim);
    sqmat_mult(U, S, US, dim);
    sqmat_mult(US, Ut, X, dim);

    //Get inverse of X
    if (hard_inv){
        if (dim == 2){
            double* Xi = new double[4];
            hard_inv_2d(X, Xi);
            sqmat_mult(A, Xi, Rot, 2);
        }
        else{
            std::cerr << "Hard inverse implemented for 2 dimensions only!!\n";
            return X;
        }
    }
    else{
        std::cerr << "Only hard inverse implemented for now!\n";
    }

    return Rot;
    /*
    Octave prototype:

    A=[1.00137,0.00535237;-0.000786181,1.00203];
    B = zeros(2,2);
    U = eye(2);
    At = A.'*A;
    T = At;

    report = true;

    for k = 1:999
        %Gschmidt
        B(:,1) = T(:,1);
        B(:,2) = T(:,2) - dot(T(:,2),B(:,1))/dot(B(:,1),B(:,1)) * B(:,1);
        if report
            disp("B is \n");
            disp(B);
            disp("\n");
            report = false;
        endif
        %QR factorization
        Q = [-B(:,1)/norm(B(:,1)),-B(:,2)/norm(B(:,2))];
        R = Q.'*T;
        %Schur decomposition
        T = R*Q;
        U = U*Q;
    endfor
    disp(U);
    disp(T);
    %sqrtm*
    R = zeros(2,2);
    sum = 0;
    for j = 1:2
        R(j,j) = sqrt(T(j,j))
        for i = j-1:-1:1
            for k = i+1:j-1
                sum += R(i,k)*R(k,j)
            endfor
            R(i,j) = (T(i,j) - sum)/(R(i,i) + R(j,j));
            sum = 0;
        endfor
    endfor

    X = U*R*U.'

    Ra = A*inv(X)
    */
}