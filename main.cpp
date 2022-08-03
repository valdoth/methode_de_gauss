/***********************************************
       NDIMBIARISOA VALDO TSIARO HASINA
                  L3 MISA
         valdotsiarohasina@gmail.com


  élimination de GAUSS par pivotage partielle
***********************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

class Gauss {
    private:
        float **matrice;
        float *vecteur;
        float *solution;
        int dim;

    public:
        // constructeur
        Gauss(string file);
        // destructeur
        ~Gauss();

        // getters
        float** getMatrice();
        float* getVecteur();
        float* getSolution();
        int getDim();

        // resolution de notre equation
        void triangularisation();
        void resolution();
};

void afficheMatrice(float **matrice, int dim);
void afficheVec(float *vecteur, int dim);
void afficherSolution(float* solution, int dim);
float **initialiseMatrice(int rows, int cols);
float *initialiseVec(int n);



int main() {
    cout << "Méthode élimination de Gauss par pivotage partielle" << endl;

    Gauss matrices("input.txt");
    cout << "La matrice de départ :\n";
    afficheMatrice(matrices.getMatrice(), matrices.getDim());
    cout << "\nLe vecteur :\n";
    afficheVec(matrices.getVecteur(), matrices.getDim());

    cout << "\nLa matrice après triangularisation: " << endl;
    matrices.triangularisation();
    afficheMatrice(matrices.getMatrice(), matrices.getDim());

    cout << "\nLa solution de notre équation linéaire est :\n";
    matrices.resolution();
    afficherSolution(matrices.getSolution(), matrices.getDim());

    return 0;
}




void afficheMatrice(float** matrice, int dim) {
    for (int i=0;i<dim;i++){
        for(int j=0; j<dim; j++)
            cout << setw(10) << matrice[i][j] << " ";
        cout << endl;
    }
}

void afficheVec(float* vecteur, int dim) {
    for(int i=0; i<dim; i++){
        cout << setw(10) << vecteur[i] << endl;
    }
}

void afficherSolution(float* solution, int dim){
    for(int i=0; i<dim; i++){
        cout << "x"<< i+1 << " = " << solution[i] << endl;
    }
}


float **initialiseMatrice(int n) {
    float **matrice(NULL);
    matrice = new float*[n];

    for (int i=0;i<n;i++)
        matrice[i] = new float[n];
    return matrice;
}

// créer un tableau à n dimension
float *initialiseVec(int n) {
    float* vecteur;
    vecteur = new float [n];
    return vecteur;
}


Gauss::Gauss(string file) {
    ifstream data (file);

    // recuperation des données dans notre fichier
    if (data){
        data >> dim;
        matrice = initialiseMatrice(dim);
        vecteur = initialiseVec(dim);
        solution = initialiseVec(dim);

        // recuperer la valeur de la matrice
        for(int i=0; i<dim;i++){
            for(int j=0; j<dim;j++){
                float temp = 0;
                data >> temp;
                matrice[i][j] = temp;
            }
        }

        // recuperation du vecteur
        for(int i=0; i<dim;i++){
            data >> vecteur[i];
        }
    }
    data.close();
}

Gauss::~Gauss() {
    for(int y = 0; y < dim; y++)    //delete colonne
        delete[] matrice[y];
    delete[] matrice;               //supprimé ligne
    delete[] vecteur;                     //supprimer le tableau contenant le second membre
    delete[] solution;              //supprimer le tableau contenant la solution
}

float** Gauss::getMatrice() {
    return matrice;
}
float* Gauss::getVecteur() {
    return vecteur;
}
float* Gauss::getSolution() {
    return solution;
}
int Gauss::getDim() {
    return dim;
}

void Gauss::triangularisation() {
    int i(0),j(0),k(0),lp(0);
    float p(0),*tp(0),t(0);
    for(k=0;k<dim;k++){
        //Recherche d'un plus grand p+ivot
        lp=k;
        p=fabs(matrice[k][k]);
        for(i=k+1;i<dim;i++){
            if(fabs(matrice[i][k])>p){
                p=fabs(matrice[i][k]); lp=i;
            }
        }
        // Permutation des lignes k et lp
        tp=matrice[k];
        matrice[k]=matrice[lp];
        matrice[lp]=tp;
        // Permutation du second membre
        t=vecteur[k];
        vecteur[k]=vecteur[lp];
        vecteur[lp]=t;
        /// Elimination
        for(i=k+1;i<dim;i++){
            for(j=k+1;j<dim;j++){
                matrice[i][j]-=(matrice[i][k]/matrice[k][k]*matrice[k][j]);
            }
            vecteur[i]-=(matrice[i][k]/matrice[k][k]*vecteur[k]);
            matrice[i][k]=0;
        }
    }
}

// resolution de la matrice triangulaire
void Gauss::resolution() {
    for(int i=dim-1;i>=0;i--){
        float sum=0;
        for(int j = i+1;j<dim;j++)
            sum += (matrice[i][j]*solution[j]);
        solution[i] = (vecteur[i]-sum)/(matrice[i][i]);
    }
}
