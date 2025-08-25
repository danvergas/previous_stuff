#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <math.h>

const int n = 13+8; //Размерность задачи
const int m = 2; //Количество параметров пристрелки
double alpha;
double eps = 1e-7;
double dr[m][m]; // матрица частных производных
double k[6][n];

void f(double* x, double* res, double*p, double t, double h);
void load(double*x, double*p);
void next_k(double* x, double*p, double h, double t);
void next_x(double* x);
double error();
void Runge_Kutta(double* x, double*p);
double Runge_Kutta_write(double* x, double*p);
double Runge_Kutta_writeS(double* x, double*p);
void residual(double *p, double *r);
void residualS(double *p, double *r, double *x);
double norm(double *x, int k);
double metric_abs(double *x, double *y, int k);
void gauss(double *b);
double shootingS(double *p);

//Диффур
void f(double* x, double* res, double*p, double t, double h){
    double u = x[1]*(1+alpha*t*t)/(cos(alpha*x[3])+2);
    res[0] = -(alpha*x[1]*x[1]*(1+alpha*t*t)*sin(alpha*x[3]))/(2*pow((cos(alpha*x[3])+2), 2))*h; //p_1
    res[1] = (-x[0])*h; //P_2
    res[2] = u*h; //y
    res[3] = x[2]*h; //x
    // Подсчет функционала
    res[4] = u*u*(cos(alpha*x[3])+2)/(1+alpha*t*t)*h;
    //Для якоби
    double u1 = (alpha*sin(alpha*x[3])/(cos(alpha*x[3])+2))*x[9]-x[7]*(1+alpha*t*t)/(2*(cos(alpha*x[3])+2));
    double u2 = (alpha*sin(alpha*x[3])/(cos(alpha*x[3])+2))*x[10]-x[8]*(1+alpha*t*t)/(2*(cos(alpha*x[3])+2));
    res[5] = (-2*alpha*sin(alpha*x[3])/(1+alpha*t*t)*u1-alpha*alpha*u*u*cos(alpha*x[3])/(1+alpha*t*t)*x[9])*h; //q_1^!
    res[6] = (-2*alpha*sin(alpha*x[3])/(1+alpha*t*t)*u2-alpha*alpha*u*u*cos(alpha*x[3])/(1+alpha*t*t)*x[10])*h; //q_1^2
    res[7] = -x[5]*h; //q_2^1
    res[8] = -x[6]*h; //q_2^2
    res[9] = x[11]*h; //dx_1
    res[10] = x[12]*h; //dx_2
    res[11] = u1*h; //dy_1
    res[12] = u2*h; //dy_2
    // Уравнение в вариациях
    res[13] = x[15]*h;// dx/d lamda_1
    res[14] = x[16]*h;// dx/d lamda_3
    res[15] = (x[13]*alpha*sin(alpha*x[3])*(1+alpha*t*t)*res[1]/pow(cos(alpha*x[3])+2, 2)+x[19]*(1+alpha*t*t)/(cos(alpha*x[3])+2))*h;// dy/d lamda_1
    res[16] = (x[14]*alpha*sin(alpha*x[3])*(1+alpha*t*t)*res[1]/pow(cos(alpha*x[3])+2, 2)+x[20]*(1+alpha*t*t)/(cos(alpha*x[3])+2))*h;// dy/d lamda_3
    u1 = -x[1]*x[1]*(1+alpha*t*t)*(cos(alpha*x[3])*cos(alpha*x[3])+2*alpha*cos(alpha*x[3])+2*alpha*sin(alpha*x[3])*sin(alpha*x[3]))/(2*pow(cos(alpha*x[3])+2, 3));
    u2 = -alpha*x[1]*(1+alpha*t*t)*sin(alpha*x[3])/((cos(alpha*x[3])+2)*(cos(alpha*x[3])+2));
    res[17] = (u1*x[13]+u2*x[19])*h;// dp_1/d lamda_1
    res[18] = (u1*x[14]+u2*x[20])*h;// dp_1/d lamda_3
    res[19] = -x[17]*h;// dp_2/d lamda_1
    res[20] = -x[18]*h;// dp_2/d lamda_3
}

//начальные данные для диффура
void load(double*x, double *p){
    x[0] = p[0];
    x[1] = p[1];
    x[2] = 1;
    x[3] = 0;
    x[4] = 0;
    //Нач данные для якоби
    x[5] = 1; //q_1^1
    x[6] = 0; //q_1^2
    x[7] = 0; //q_2^1
    x[8] = 1; //q_2^2
    x[9] = 0; //dx_1
    x[10] = 0; //dx_2
    x[11] = 0; //dy_1
    x[12] = 0; //dy_2
    x[13] = 0;
    x[14] = 0;
    x[15] = 0;
    x[16] = 0;
    x[17] = 1;
    x[18] = 0;
    x[19] = 0;
    x[20] = 1;
}

void residual(double *p, double *r){
    double x[n];
    load(x, p);
    Runge_Kutta(x, p);
    r[0] = x[0];
    r[1] = x[2];
}

void residualS(double *p, double *r, double *x){
    load(x, p);
    Runge_Kutta(x, p);
    r[0] = x[0];
    r[1] = x[2];
}

void next_k(double* x,  double*p, double h, double t){
    //-------------------------
    //Подсчет K1 (k[0])
    f(x, k[0], p, t, h);
    //-------------------------
    //Подсчет K2 (k[1])
    double tmp[n];
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/2.;
    }
    f(tmp, k[1], p, t+h/2., h);
    //-------------------------
    //Подсчет K3 (k[2])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/4.+k[1][i]/4.;
    }
    f(tmp, k[2], p, t+h/2., h);
    //-------------------------
    //Подсчет K4 (k[3])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]-k[1][i]+2.*k[2][i];
    }
    f(tmp, k[3], p, t+h, h);
    //-------------------------
    //Подсчет K5 (k[4])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*7./27.+k[1][i]*10./27.+k[3][i]/27.;
    }
    f(tmp, k[4], p, t+2.*h/3., h);
    //-------------------------
    //Подсчет K6 (k[5])
    for(int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*0.0448-k[1][i]*0.2+k[2][i]*0.8736+k[3][i]*0.0864-k[4][i]*0.6048;
    }
    f(tmp, k[5], p, t+h/5., h);
    //-------------------------
}

void next_x(double* x){
    for (int i=0; i<n; i++){
        x[i] = x[i]+ k[0][i]/24.+k[3][i]*5./48.+k[4][i]*27./56.+k[5][i]*125./336.;
    }
}

double error(){
    double err = 0;
    double tmp = 0;
    for (int i=0; i<n; i++){
        tmp = -0.125*k[0][i]-k[2][i]*224./336.-k[3][i]*0.0625+k[4][i]*162./336.+k[5][i]*125./336.;
        err+=tmp*tmp;
    }
    err = sqrt(err);
    return err;
}

void Runge_Kutta(double* x, double *p){
    double t = 0;
    double T = 1;
    double x_tmp;
    double I=0;
    double fac;
    double h_new, h;
    h_new = 0.01;
    while (t<T){
        double err = 1;
        while (err > eps){
            h = h_new;
            next_k(x, p, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, p, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
    }
}

std::string to_str(double number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << number;
    std::string str = oss.str();
    
    int lastSignificant = str.find_last_not_of('0');
    if (lastSignificant != std::string::npos) {
        if (str[lastSignificant] == '.') {
            lastSignificant--;
        }
        str = str.substr(0, lastSignificant + 1);
    }
    if (str.empty()) {
        return "0";
    }

    return str;
}


double Runge_Kutta_write(double* x, double *p){
    double t = 0;
    double T = 1;
    double I=0;
    double x_tmp;
    double h_new, h;
    h_new = 0.01;
    std::ofstream fout("output_" +to_str(alpha)+".csv");
    fout << "t,p1,p2,y,x,z,det" << std::endl;
    while (t<T){
        fout << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] <<  "," << x[9]*x[12]-x[10]*x[11] << std::endl;
        double err = 1;
        double fac;
        while (err > eps){
            h = h_new;
            next_k(x, p, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, p, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
    }
    fout << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << "," << x[9]*x[12]-x[10]*x[11] << std::endl;
    fout.close();
    return x[4];
}

double Runge_Kutta_writeS(double* x, double *p){
    double t = 0;
    double T = 1;
    double I=0;
    double x_tmp;
    double h_new, h;
    h_new = 0.01;
    while (t<T){
        double err = 1;
        double fac;
        while (err > eps){
            h = h_new;
            next_k(x, p, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, p, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
    }
    return x[4];
}

double norm(double *x, int k){
    double sum =0;
    for (int i=0; i<k; i++){
        sum += fabs(x[i]);
    }
    return sum;
}

double metric_abs(double *x, double *y, int k){
    double sum =0;
    for (int i=0; i<k; i++){
        sum += fabs(y[i]-x[i]);
    }
    return sum;
}

void gauss(double *b){
    double mx;
    int mxn;
    for (int i=0;i<m;i++){
        mx = fabs(dr[i][i]);
        mxn = i;
        //Находим строку с наибольшим значеием на данной итерации
        for (int j=i+1; j<m; j++){
            if (fabs(dr[j][i]) > mx){
                mx = fabs(dr[j][i]);
                mxn = j;
            }
        }

        for (int j=i; j<m; j++){
            mx = dr[i][j];
            dr[i][j] = dr[mxn][j];
            dr[mxn][j] = mx;
        }
        mx = b[i];
        b[i] = b[mxn];
        b[mxn] = mx;

        mx = dr[i][i];
        dr[i][i] = 1;
        b[i]=b[i]/mx;
        for (int j=i+1; j<m;j++){
            dr[i][j]=dr[i][j]/mx;
        }

        for (int j=0; j<m; j++) {
			if(j != i){
				for (int l=i+1; l<m; l++){
                    dr[j][l] -= dr[j][i]*dr[i][l];
                }
                b[j] -= dr[j][i]*b[i];
				dr[j][i] = 0;
			}
        }
    }
}

double shooting(double *p){
    double d=1e-4;
    double r[m], r_new[m], p_new[m], h[m];
    double r_norm, r_new_norm;
    double gamma = 1.;
    int counter = 0;
    while (true){
        residual(p, r);
        r_norm = norm(r, m);
        if (r_norm < eps){
            return 1;
        }
        for (int i=0; i<m;i++){
            h[i] = r[i];
            for (int j=0; j<m;j++){
                p_new[j]=p[j];
            }
            p_new[i]+=d;
            residual(p_new, r_new);
            for (int j=0; j<m;j++){
                dr[j][i]=(r_new[j]-r[j])/d;
            }
        }

        gauss(h);
        counter++;

        while (true){
            for(int i=0;i<m;i++){
                p_new[i] = p[i] - h[i]*gamma;
            }
            residual(p_new, r_new);
            r_new_norm = norm(r_new, m);
            if (r_new_norm<eps){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                return counter;
            }
            else if (r_new_norm < r_norm){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                break;
            }
            else if (gamma>1e-20){
                gamma = gamma/2;
            }
            else {
                throw std::runtime_error("Gamma less then eps");
            }
        }
        if (counter > 10000){
            throw std::runtime_error("Counter too big");
        }
    }
}

double shootingS(double *p){
    double r[m], r_new[m], p_new[m], h[m], x[n];
    double r_norm, r_new_norm;
    double gamma = 1.;
    int counter = 0;
    while (true){
        residualS(p, r, x);
        r_norm = norm(r, m);
        if (r_norm < eps){
            return 1;
        }
        dr[0][0]=x[17];
        dr[0][1]=x[18];
        dr[1][0]=x[15];
        dr[1][1]=x[16];

        gauss(h);
        counter++;

        while (true){
            for(int i=0;i<m;i++){
                p_new[i] = p[i] - h[i]*gamma;
            }
            residual(p_new, r_new);
            r_new_norm = norm(r_new, m);
            if (r_new_norm<eps){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                return counter;
            }
            else if (r_new_norm < r_norm){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                break;
            }
            else if (gamma>1e-20){
                gamma = gamma/2;
            }
            else {
                throw std::runtime_error("Gamma less then eps");
            }
        }
        if (counter > 10000){
            throw std::runtime_error("Counter too big");
        }
    }
}



int main(){
    int counter;
    double p[m];
    double x[n];
    double func;
    p[0] = 0;
    p[1] = -3;
    alpha = -0.1;
    for (int i=0; i<=110; i++){
        alpha += 0.1;
        try{
            counter = shootingS(p);
        }
        catch (std::runtime_error& e){
            std::cout << "Error alpha=" << alpha <<std::endl;
            std::cout << e.what() << std::endl;
            break;
        }
        //if ((i == 0) || (i == 1) || (i == 10) || (i == 100) || (i == 105) || (i == 106) || (i == 107) || (i == 108)){
            std::cout << "alpha=" << alpha << ":" << std::endl;
            std::cout << "Calculated p = (";
            for(int jm=0;jm<m;jm++){
                std::cout << p[jm];
                if (jm != m-1) std::cout << ", ";;
            }
            std::cout << ")"<< std::endl;
            load(x, p);
            func = Runge_Kutta_writeS(x, p);
            std::cout << "Calculated value of functional: " << func << std::endl;
            std::cout << "-------------------------------------" << std::endl;
        //}
    }

    return 1;
}