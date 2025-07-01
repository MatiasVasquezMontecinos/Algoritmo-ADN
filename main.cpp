#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

using Vector2fList = vector<Vector2f>;
using MatrixXi = Matrix<int, Dynamic, Dynamic>;
using MatrixXf = Matrix<float, Dynamic, Dynamic>;

mt19937 rng(42);

// Leer CSV y obtener vectores 2D
Vector2fList leerCSV(const string& filename) {
    ifstream file(filename);
    string line;
    Vector2fList data;

    getline(file, line); // Saltar encabezado

    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        int col = 0;
        float dim1 = NAN, dim2 = NAN;
        int chamber_ok = 0, congress_ok = 0;

        while (getline(ss, cell, ',')) {
            if (col == 2 && cell == "House") chamber_ok = 1;
            if (col == 1 && cell == "94") congress_ok = 1;
            if (col == 10) dim1 = stof(cell);
            if (col == 11) dim2 = stof(cell);
            col++;
        }
        if (chamber_ok && congress_ok && !isnan(dim1) && !isnan(dim2)) {
            data.emplace_back(Vector2f(dim1, dim2));
        }
    }
    return data;
}

// Matriz de distancias euclídeas
MatrixXf calcular_matriz_distancias(const Vector2fList& posiciones) {
    int n = posiciones.size();
    MatrixXf D(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            float dist = (posiciones[i] - posiciones[j]).norm();
            D(i, j) = D(j, i) = dist;
        }
    return D;
}

// Generar población inicial
MatrixXi generar_poblacion(int pob_max, int n, int quorum) {
    MatrixXi poblacion = MatrixXi::Zero(pob_max, n);
    for (int i = 0; i < pob_max; ++i) {
        vector<int> idx(n);
        iota(idx.begin(), idx.end(), 0);
        shuffle(idx.begin(), idx.end(), rng);
        for (int j = 0; j < quorum; ++j)
            poblacion(i, idx[j]) = 1;
    }
    return poblacion;
}

// Evaluar población
VectorXf evaluar_poblacion(const MatrixXi& poblacion, const MatrixXf& D) {
    int N = poblacion.rows();
    VectorXf Z = VectorXf::Zero(N);

    for (int i = 0; i < N; ++i) {
        vector<int> idx;
        for (int j = 0; j < poblacion.cols(); ++j)
            if (poblacion(i, j)) idx.push_back(j);
        float suma = 0.0;
        for (size_t a = 0; a < idx.size(); ++a)
            for (size_t b = a + 1; b < idx.size(); ++b)
                suma += D(idx[a], idx[b]);
        Z(i) = suma;
    }
    return Z;
}

// Selección de padres
pair<RowVectorXi, RowVectorXi> seleccionar_padres(const MatrixXi& poblacion) {
    uniform_int_distribution<int> dist(0, poblacion.rows() - 1);
    return {poblacion.row(dist(rng)), poblacion.row(dist(rng))};
}

// Cruzamiento
pair<RowVectorXi, RowVectorXi> cruzar(const RowVectorXi& p1, const RowVectorXi& p2) {
    int punto = uniform_int_distribution<int>(1, p1.size() - 2)(rng);
    RowVectorXi h1(p1.size()), h2(p1.size());
    h1 << p1.head(punto), p2.tail(p1.size() - punto);
    h2 << p2.head(punto), p1.tail(p1.size() - punto);
    return {h1, h2};
}

// Mutación
void mutar(RowVectorXi& individuo, double prob = 0.1) {
    uniform_real_distribution<double> dist(0, 1);
    if (dist(rng) < prob) {
        vector<int> unos, ceros;
        for (int i = 0; i < individuo.size(); ++i)
            (individuo(i) ? unos : ceros).push_back(i);
        if (!unos.empty() && !ceros.empty()) {
            int i = unos[rng() % unos.size()];
            int j = ceros[rng() % ceros.size()];
            swap(individuo(i), individuo(j));
        }
    }
}

// Verificar quórum
bool verificar(const RowVectorXi& ind, int quorum) {
    return ind.sum() == quorum;
}

// Crear nueva generación
MatrixXi nueva_generacion(const MatrixXi& pob, const MatrixXf& D, int quorum, double prob_mut) {
    vector<RowVectorXi> nueva;
    while (nueva.size() < pob.rows()) {
        auto [p1, p2] = seleccionar_padres(pob);
        auto [h1, h2] = cruzar(p1, p2);
        mutar(h1, prob_mut); mutar(h2, prob_mut);
        if (verificar(h1, quorum)) nueva.push_back(h1);
        if (nueva.size() < pob.rows() && verificar(h2, quorum)) nueva.push_back(h2);
    }
    MatrixXi nueva_pob(nueva.size(), pob.cols());
    for (size_t i = 0; i < nueva.size(); ++i)
        nueva_pob.row(i) = nueva[i];
    VectorXf Z = evaluar_poblacion(nueva_pob, D);
    vector<int> orden(nueva_pob.rows());
    iota(orden.begin(), orden.end(), 0);
    sort(orden.begin(), orden.end(), [&Z](int a, int b) { return Z(a) < Z(b); });

    MatrixXi ordenada(nueva_pob.rows(), pob.cols());
    for (size_t i = 0; i < orden.size(); ++i)
        ordenada.row(i) = nueva_pob.row(orden[i]);
    return ordenada;
}

int main() {
    auto start = chrono::steady_clock::now();

    int quorum = 216, pob_max = 38;
    auto datos = leerCSV("H094_members.csv");
    int n = datos.size();

    cout << "🔧 Calculando matriz de distancias..." << endl;
    MatrixXf D = calcular_matriz_distancias(datos);

    cout << "🧬 Generando población inicial..." << endl;
    MatrixXi poblacion = generar_poblacion(pob_max, n, quorum);

    cout << "📊 Evaluando población..." << endl;
    VectorXf fitness = evaluar_poblacion(poblacion, D);
    for (int i = 0; i < fitness.size(); ++i)
        cout << "Individuo " << i << ": Z = " << fitness(i) << endl;

    cout << "🔁 Evolución por generaciones..." << endl;
    float mejor_Z = numeric_limits<float>::max();
    RowVectorXi mejor_ind(n);
    int sin_mejora = 0;

    for (int gen = 0; gen < 10000; ++gen) {
        poblacion = nueva_generacion(poblacion, D, quorum, 0.17);
        VectorXf f = evaluar_poblacion(poblacion, D);
        float fmin = f.minCoeff();
        if (fmin < mejor_Z) {
            mejor_Z = fmin;
            mejor_ind = poblacion.row(min_element(f.data(), f.data() + f.size()) - f.data());
            sin_mejora = 0;
        } else sin_mejora++;

        if (gen % 500 == 0)
            cout << "Gen " << gen << " Mejor Z = " << mejor_Z << endl;
        if (abs(mejor_Z - 9686.93f) < 0.01f || sin_mejora > 1500) break;
    }

    cout << "🎯 Mejor Z encontrado: " << mejor_Z << endl;
    cout << "⏱ Tiempo total: " <<
        chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count() << "s" << endl;

    return 0;
}
