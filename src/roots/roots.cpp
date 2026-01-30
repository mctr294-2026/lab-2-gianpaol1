#include "roots.hpp"
#include <cmath>

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    const double tolerance = 1e-6;
    const int Max_Iter = 1e6;
    double fa = f(a);
    double fb = f(b);
    if(fa == 0.0){
        *root = a;
        return true;
    }
    if(fb == 0.0){
        *root = b;
        return true;
    }
    if (fa * fb > 0) {
        return false;
    }


    for (int i = 0; i < Max_Iter; i++){
        double c = (a+b)*0.5;
        double fc = f(c);

        if(std::fabs(fc) <= tolerance || std::fabs(b-a) <= tolerance) {
            *root = c;
            return true;
        }
        if(fa*fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        };
    }
    return false;
    };


bool regula_falsi(std::function<double(double)> f, double a, double b, double *root){
    const double tolerance = 1e-6;
    const int Max_Iter = 1e6;

    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) {
        return false; 
    }
    for (int i = 0; i < Max_Iter; ++i) {
        double c = a - ((f(a)*(b-a))/(fb - fa)); // formula given
        double fc = f(c);
        if (std::fabs(fc) < tolerance) { 
            *root = c; 
            return true;
        }
        if (fa * fc < 0) { 
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return false; 
}



bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root){
    const double tolerance = 1e-6;
    const int Max_Iter = 1e6;
    
   

    double x_n = c; // set initial guess
    if (f(a) * f(b) > 0){ 
        return false;
    }
    for(int i = 0; i < Max_Iter; i++){
        double fx_n =f(x_n);
        double gx_n = g(x_n);
        if (std::fabs(gx_n) < tolerance) { // check derivative isnt too small/zero (to avoid division by zero error)
            return false;
        }

        double x_new = x_n - (fx_n/gx_n); 
        double fx = f(x_new);
        if(std::fabs(fx) < tolerance){ 
            *root = x_new;
            return true;
        }

        x_n = x_new; 
    }
    return false; 
}

bool secant(std::function <double(double)> f, double a, double b, double c, double *root){ 
    const double tolerance = 1e-6;
    const int Max_Iter = 1e6;
    
    double x0 = c; // second initial guess
    double x1 = c+0.1; // first initial guess

    if(f(a) * f(b) > 0){
        return false;
    }
    for(int i = 0; i < Max_Iter; i++){
        if (std::fabs((f(x0) - f(x1))) < tolerance){ // avoid division by zero
            return false; 
        }
        double x_new = x0 - f(x0)*((x0 - x1)/(f(x0) - f(x1)));
        if (std::fabs(f(x_new)) < tolerance){
            *root = x_new;
            return true;
        }
        else {
             x1 = x0;
            x0 = x_new;
        }
    }
    return false; 
}

