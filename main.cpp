#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

#define PI 3.14159265358979323846
#define a 0.0


double f_der(double x, double y){
    return cos(x);
}

double f(double x){
    return sin(x);
}

vector<double> calculate_xi_h (int n){
    vector<double> x_i(n, 0);
    double h = PI / (n - 1);
    x_i[0] = a;
    for (int i = 0; i < n - 1; i++){
        x_i[i + 1] = x_i[i] + h;
    }
    return x_i;
}

vector<double> euler_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    y_i[0] = f(x_i[0]);
    double h = PI / (n - 1);
    for (int i = 0; i < n - 1; i++){
        y_i[i + 1] = y_i[i] + h * f_der(x_i[i], y_i[i]);
    }
    return y_i;
}

vector<double> improved_euler_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    double h = PI / (n - 1);
    y_i[0] = f(x_i[0]);
    for (int i = 0; i < n - 1; i++){
        double k1 = f_der(x_i[i], y_i[i]);
        double k2 = f_der(x_i[i] + h, y_i[i] + h * k1);
        y_i[i + 1] = y_i[i] + h * (k1 + k2) / 2;
    }
    return y_i;
}


vector<double> runge_kutta_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    double h = PI / (n - 1);
    y_i[0] = f(x_i[0]);
    for (int i = 0; i < n - 1; i++){
        double k1 = f_der(x_i[i], y_i[i]);
        double k2 = f_der(x_i[i] + h/2, y_i[i] + h / 2 * k1);
        double k3 = f_der(x_i[i] + h/2, y_i[i] + h / 2 * k2);
        double k4 = f_der(x_i[i] + h, y_i[i] + h * k3);
        y_i[i + 1] = y_i[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    return y_i;
}

vector<double> calculate_x(int n){
    vector<double> x_i(n);
    for (int i = 0; i < n; i++) {
        x_i[i] = (i * PI / (n - 1));
    }
    return x_i;
}

vector<double> exact_solution(int n, vector<double> x_i){
    vector<double> y_i(n);
    for (int i = 0; i < n; i++){
        y_i[i] = f(x_i[i]);
    }
    return y_i;
}
/*
 * task 1 is exact solution
 * task 2 is Euler's solution
 * task 3 is improved Euler's solution
 * task 4 is Runge-Kutta solution
 * task 5 is find local errors of the Euler's solution
 * task 6 is find local errors of the improved Euler's solution
 * task 7 is find local errors of the Runge-Kutta solution
 * task 8 is find global errors of the Euler's solution
 * task 9 is find global errors of the improved Euler's solution
 * task 10 is fidn global errors of the Runge-Kutta solution
 */

int main() {
    int n, n1, n2, task;
    cin >> n >> n1 >> n2 >> task;
    double h = PI / (n - 1);
    vector<double> local_error;
    if (task == 1){
        cout << "xi=\n";
        vector<double> x = calculate_x(n);
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\ny(xi)=\n";
        vector<double> y = exact_solution(n, x);
        for (auto y_i : y) cout << fixed << setprecision(5) << y_i << " ";
    } else if (task == 2){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\nEuler_yi=\n";
        vector<double> y = euler_method(n, x);
        for (auto y_i : y) cout << fixed << setprecision(5) << y_i << " ";
    } else if (task == 3){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\niEuler_yi=\n";
        vector<double> y = improved_euler_method(n, x);
        for (auto y_i : y) cout << fixed << setprecision(5) << y_i << " ";
    } else if (task == 4){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\nRK4_yi=\n";
        vector<double> y = runge_kutta_method(n, x);
        for (auto y_i : y) cout << fixed << setprecision(5) << y_i << " ";
    } else if (task == 5){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        vector<double> euler_values = euler_method(n, x);
        vector<double> exact_values = exact_solution(n, x);
        for (int i = 0; i < n; i++){
            local_error.push_back(abs(exact_values[i] - euler_values[i]));
        }
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\nEuler_LE(xi)=\n";
        for (auto err : local_error) cout << fixed << setprecision(5) << err << " ";
    } else if (task == 6){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        vector<double> improved_euler_values = improved_euler_method(n, x);
        vector<double> exact_values = exact_solution(n, x);
        for (int i = 0; i < n; i++){
            local_error.push_back(abs(exact_values[i] - improved_euler_values[i]));
        }
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\niEuler_LE(xi)=\n";
        for (auto err : local_error) cout << fixed << setprecision(5) << err << " ";
    } else if (task == 7){
        cout << "xi=\n";
        vector<double> x = calculate_xi_h(n);
        vector<double> runge_kutta_values = runge_kutta_method(n, x);
        vector<double> exact_values = exact_solution(n, x);
        for (int i = 0; i < n; i++){
            local_error.push_back(abs(exact_values[i] - runge_kutta_values[i]));
        }
        for (auto x_i : x) cout << fixed << setprecision(5) <<  x_i << " ";
        cout << "\nRK4_LE(xi)=\n";
        for (auto err : local_error) cout << fixed << setprecision(5) << err << " ";
    } else if (task == 8) {
        cout << "ni=\n";
        for (int i = n1; i <= n2; i++) cout << i << " ";
        vector<double> global_errors;

        for (int n_ = n1; n_ <= n2; n_++) {
            vector<double> x_h = calculate_xi_h(n_);
            vector<double> euler_values = euler_method(n_, x_h);
            vector<double> x = calculate_x(n_);
            vector<double> exact_values = exact_solution(n_, x);
            double max_err = 0.0;
            for (int i = 0; i < n_; i++) {
                double local_err = abs(euler_values[i] - exact_values[i]);
                max_err = max(max_err, local_err);
            }
            global_errors.push_back(max_err);
        }
        cout << "\nEuler_GE(n)=\n";
        for (auto i : global_errors) cout << fixed << setprecision(5) << i << " ";
    } else if (task == 9){
        cout << "ni=\n";
        for (int i = n1; i <= n2; i++) cout << i << " ";
        vector<double> global_errors;

        for (int n_ = n1; n_ <= n2; n_++) {
            vector<double> x_h = calculate_xi_h(n_);
            vector<double> improved_euler_values = improved_euler_method(n_, x_h);
            vector<double> x = calculate_x(n_);
            vector<double> exact_values = exact_solution(n_, x);
            double max_err = 0.0;
            for (int i = 0; i < n_; i++) {
                double local_err = abs(improved_euler_values[i] - exact_values[i]);
                max_err = max(max_err, local_err);
            }
            global_errors.push_back(max_err);
        }
        cout << "\niEuler_GE(n)=\n";
        for (auto i : global_errors) cout << fixed << setprecision(5) << i << " ";
    } else if (task == 10){
        cout << "ni=\n";
        for (int i = n1; i <= n2; i++) cout << i << " ";
        vector<double> global_errors;
        for (int n_ = n1; n_ <= n2; n_++) {
            vector<double> x_h = calculate_xi_h(n_);
            vector<double> runge_kutta_values = runge_kutta_method(n_, x_h);
            vector<double> x = calculate_x(n_);
            vector<double> exact_values = exact_solution(n_, x);
            double max_err = 0.0;
            for (int i = 0; i < n_; i++) {
                double local_err = abs(runge_kutta_values[i] - exact_values[i]);
                max_err = max(max_err, local_err);
            }
            global_errors.push_back(max_err);
        }

        cout << "\nRK4_GE(n)=\n";
        for (auto i : global_errors) cout << fixed << setprecision(5) << i << " ";
    }

}

