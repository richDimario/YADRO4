#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

const double PI = 3.14159265358979323846;

class FFT {
public:
    static void fft(vector<complex<double>>& a, bool invert) {
        int n = a.size();
        if (n == 1) return;

        // Разделение на подмассивы в зависимости от кратности длины
        if (n % 2 == 0) {
            // Кратные 2
            fft2(a, invert);
        }
        else if (n % 3 == 0) {
            // Кратные 3
            fft3(a, invert);
        }
        else if (n % 5 == 0) {
            // Кратные 5
            fft5(a, invert);
        }
        else {
            throw runtime_error("FFT длина должна быть кратной 2, 3 или 5");
        }
    }

    // Прямое преобразование Фурье для кратных 2
    static void fft2(vector<complex<double>>& a, bool invert) {
        int n = a.size();
        if (n == 1) return;

        vector<complex<double>> a0(n / 2), a1(n / 2);
        for (int i = 0; i < n / 2; ++i) {
            a0[i] = a[2 * i];
            a1[i] = a[2 * i + 1];
        }

        fft2(a0, invert);
        fft2(a1, invert);

        double angle = 2 * PI / n * (invert ? -1 : 1);
        complex<double> w(1), wn(cos(angle), sin(angle));
        for (int i = 0; i < n / 2; ++i) {
            a[i] = a0[i] + w * a1[i];
            a[i + n / 2] = a0[i] - w * a1[i];
            if (invert) {
                a[i] /= 2;
                a[i + n / 2] /= 2;
            }
            w *= wn;
        }
    }

    // Прямое преобразование Фурье для кратных 3
    static void fft3(vector<complex<double>>& a, bool invert) {
        int n = a.size();
        if (n == 1) return;

        vector<complex<double>> a0(n / 3), a1(n / 3), a2(n / 3);
        for (int i = 0; i < n / 3; ++i) {
            a0[i] = a[3 * i];
            a1[i] = a[3 * i + 1];
            a2[i] = a[3 * i + 2];
        }

        fft3(a0, invert);
        fft3(a1, invert);
        fft3(a2, invert);

        double angle = 2 * PI / n * (invert ? -1 : 1);
        complex<double> w1(1), w2(cos(angle), sin(angle));
        for (int i = 0; i < n / 3; ++i) {
            complex<double> t = w1 * a1[i] + w2 * a2[i];
            a[i] = a0[i] + t;
            a[i + n / 3] = a0[i] - t;
            a[i + 2 * n / 3] = a0[i] - t;
            if (invert) {
                a[i] /= 3;
                a[i + n / 3] /= 3;
                a[i + 2 * n / 3] /= 3;
            }
            w1 *= w2;
        }
    }

    // Прямое преобразование Фурье для кратных 5
    static void fft5(vector<complex<double>>& a, bool invert) {
        int n = a.size();
        if (n == 1) return;

        vector<complex<double>> a0(n / 5), a1(n / 5), a2(n / 5), a3(n / 5), a4(n / 5);
        for (int i = 0; i < n / 5; ++i) {
            a0[i] = a[5 * i];
            a1[i] = a[5 * i + 1];
            a2[i] = a[5 * i + 2];
            a3[i] = a[5 * i + 3];
            a4[i] = a[5 * i + 4];
        }

        fft5(a0, invert);
        fft5(a1, invert);
        fft5(a2, invert);
        fft5(a3, invert);
        fft5(a4, invert);

        double angle = 2 * PI / n * (invert ? -1 : 1);
        complex<double> w1(1), w2(cos(angle), sin(angle));
        for (int i = 0; i < n / 5; ++i) {
            complex<double> t = w1 * a1[i] + w2 * (a2[i] + a3[i] + a4[i]);
            a[i] = a0[i] + t;
            a[i + n / 5] = a0[i] - t;
            a[i + 2 * n / 5] = a0[i] - t;
            a[i + 3 * n / 5] = a0[i] - t;
            a[i + 4 * n / 5] = a0[i] - t;
            if (invert) {
                a[i] /= 5;
                a[i + n / 5] /= 5;
                a[i + 2 * n / 5] /= 5;
                a[i + 3 * n / 5] /= 5;
                a[i + 4 * n / 5] /= 5;
            }
            w1 *= w2;
        }
    }

    // Обратное преобразование Фурье
    static void ifft(vector<complex<double>>& a) {
        fft(a, true);
    }
};

void printVector(const vector<complex<double>>& v) {
    for (const auto& x : v) {
        cout << "(" << real(x) << "," << imag(x) << ") ";
    }
    cout << endl;
}

double calculateError(const vector<complex<double>>& original, const vector<complex<double>>& transformed) {
    double error = 0.0;
    for (size_t i = 0; i < original.size(); ++i) {
        double real_diff = real(original[i]) - real(transformed[i]);
        double imag_diff = imag(original[i]) - imag(transformed[i]);
        error += real_diff * real_diff + imag_diff * imag_diff;
    }
    return sqrt(error);
}

int main() {
    setlocale(LC_ALL, "Russian");
    const int N = 15; // Можно заменить на размер, кратный 2, 3, или 5

    // Генерация случайных комплексных данных
    vector<complex<double>> data(N);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(-1.0, 1.0);

    for (auto& x : data) {
        x = complex<double>(dis(gen), dis(gen));
    }

    cout << "Исходные данные:" << endl;
    printVector(data);

    vector<complex<double>> data_copy = data;

    // Прямое преобразование Фурье
    FFT::fft(data, false);

    cout << "После Прямого Преобразования Фурье:" << endl;
    printVector(data);

    // Обратное преобразование Фурье
    FFT::ifft(data);

    cout << "После Обратного Преобразования Фурье:" << endl;
    printVector(data);

    // Сравнение ошибки
    double error = calculateError(data_copy, data);
    cout << "Ошибка между входными и выходными данными: " << error << endl;

    return 0;
}
