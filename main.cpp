#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

#define infinito 0x7FFFFFFF

struct Cidade {
    int a;
    int b;
};

struct Celula {
    double val;
    vector <int> caminhos;
    Celula* prox;
};

/
class CidadesVisitas {
private:
    vector <Celula*> ini;
    vector <Celula*> fim;

public:

    CidadesVisitas (int numCidades) {
        ini.resize(numCidades+1);
        fim.resize(numCidades+1);
        for (int x = 0; x < numCidades+1; x++) {
            ini[x] = NULL;
            fim[x] = NULL;
        }
    }


    void inserir (int vertice, double val, vector <int> caminho) {

        Celula *nova = new Celula;


        if (ini[vertice] == NULL) {

            nova->val = val;
            nova->caminhos.resize (caminho.size());
            nova->caminhos = caminho;
            nova->prox = NULL;

            ini[vertice] = nova;
            fim[vertice] = ini[vertice];

        } else {

            nova->val = val;
            nova->caminhos.resize (caminho.size());
            nova->caminhos = caminho;
            nova->prox = NULL;

            fim[vertice]->prox = nova;
            fim[vertice] = nova;

            nova = NULL;
        }
    }


    double getVal (int vertice, vector <int> caminho) {
        double resposta = infinito;
        Celula *temp = new Celula;
        temp = ini[vertice];
        bool encontrado = true;

        if (caminho.size() == 0) {
            return (ini[vertice]->val);
        } else {

            while (temp != NULL) {

                for (int x = 0; x < caminho.size(); x++ ) {
                    if (caminho[x] != temp->caminhos[x]) {
                        encontrado = false;
                    }
                }

                if (encontrado) {
                    resposta = temp->val;
                    break;
                }
                encontrado = true;
                temp = temp->prox;
            }
            return (resposta);
        }
    }
};


class TSPD {
private:

    FILE *arqEntrada;
    FILE *arqSaida;
    int n;
    int controle;
    Cidade *cidades;
    vector <int> auxComb;
    vector <int> vetor;
    vector < vector <int> > combinacaoElementos;
    vector < vector <double> > mAdjacentes;
    vector <CidadesVisitas*>  custo;
public:

    TSPD () {
        char nome [256];
        cout << "Informar nome arquivo entrada : ";
        cin >> nome;
        arqEntrada = fopen (nome, "r");
        arqSaida = fopen ("caminho.txt", "w");
        fscanf (arqEntrada, "%d\n", &n);
        vetor.resize(n-1);
        cidades = new Cidade [n+1];
        cidades[0].a = 0;
        cidades[0].b = 0;
        for (int a = 1; a <= n; a++) {
            fscanf (arqEntrada, "%d %d\n", &cidades[a].a, &cidades[a].b);
        }
        fclose(arqEntrada);
        preencherMAdjacentes();
        for (int a = 0; a < (n-1); a++) {
            vetor[a] = a+2;
        }
    }

    void algoritmoTSPD () {
        time_t inicio, fim;
        inicio = clock();
        double custoAux;
        double custoCaminho;
        custo.resize(n-1);
        vector <int> inicial;
        inicial.resize(1);
        inicial[0] = 1;
        for (int a = 0; a < n-1; a++) {
            custo[a] = new CidadesVisitas (n);
        }


        for (int a = 2; a <= n; a++ ) {
            custo[0]->inserir(a,mAdjacentes[a][1], extrair (inicial, 1));
        }


        for (int s = 1; s < (n-1); s++) {

            for (int a = 2; a <= n; a++) {
                controle = 0;
                int numComb = calcularNumCombinacoes(extrair(vetor,a).size(),s);
                alocarCombinacoes (numComb,s);
                auxComb.resize(s);
                combinacoes (extrair(vetor,a),s,0,auxComb);
                for (int j = 0; j < numComb; j++) {
                    vector <int> vetAux = removerVetor (s,j);
                    custoCaminho = infinito;

                    for (int k = 0; k < s; k++) {
                        custoAux = mAdjacentes [a][vetAux[k]] + custo[s-1]->getVal(vetAux[k], extrair (vetAux, vetAux[k]));

                        if (custoAux < custoCaminho) {
                            custoCaminho = custoAux;
                        }
                    }
                    custo[s]->inserir (a,custoCaminho, vetAux);

                }

            }
        }
        custoCaminho = infinito;
        int cidadeBase = 0;
        for (int a = 2; a <= n; a++) {

            custoAux = mAdjacentes[1][a] + custo[n-2]->getVal(a, extrair (vetor, a));

            if (custoAux < custoCaminho) {
                custoCaminho = custoAux;
                cidadeBase = a;
            }
        }


        fprintf(arqSaida, "NUMERO CIDADES : %d\nCAMINHO\n%d ", n, 1);
        recuperarCaminho (cidadeBase, vetor, (n-2));
        fim = clock();
        fprintf(arqSaida, "\nCUSTO : %0.2f", custoCaminho);
        fprintf(arqSaida,"\nTEMPO : %0.2f",  (( double(fim) - double(inicio) )/CLOCKS_PER_SEC));
        fclose(arqSaida);

    }



    ~TSPD () {}

private:


    void alocarCombinacoes (int linhas, int colunas) {
        combinacaoElementos.resize(linhas);
        for (int a = 0; a < linhas; a++) {
            combinacaoElementos[a].resize(colunas);
        }
    }

    double calcularDistancia (Cidade primeira, Cidade segunda) {
        double distancia = sqrt ( pow ((primeira.a - segunda.a),2) + pow ((primeira.b - segunda.b), 2));
        return (distancia);
    }


    void preencherMAdjacentes () {
        mAdjacentes.resize(n+1);
        for (int a = 0; a < n+1; a++) {
            mAdjacentes[a].resize(n+1);
        }
        for (int a = 0; a <= n; a++) {
            for (int b = 0; b <= n; b++) {
                if (a == 0 || b == 0 || a == b) {
                    mAdjacentes[a][b] = 0.0;
                } else {
                    mAdjacentes[a][b] = calcularDistancia (cidades[a], cidades[b]);
                }
            }
        }
        delete (cidades);
    }


    vector <int> extrair (vector <int> vet, int elemento) {
        vector <int> resposta;
        resposta.resize(vet.size()-1);
        int pos = 0;
        for (int a = 0; a < vet.size(); a++) {
            if (vet[a] != elemento) {
                resposta[pos++] = vet[a];
            }
        }
        return (resposta);
    }


    int calcularNumCombinacoes (int numEle, int tamConj) {
        int combinacaoDividendo = 1;
        int combinacaDivisor = 1;
        int aux = tamConj;
        for (int a = 0; a < tamConj; a++){
            if (numEle != aux) {
                combinacaoDividendo *= numEle--;
            } else if (numEle == aux) {
                numEle--;
                aux--;
            }
        }
        int aux2 = aux;
        for (int a = 0; a < aux2; a++) {
            combinacaDivisor *= aux--;
        }
        return (combinacaoDividendo/combinacaDivisor);
    }


    vector <int> removerVetor (int s, int escolher) {
        vector <int> resp;
        resp.resize(s);
        for (int a = 0; a < s; a++) {
            resp[a] = combinacaoElementos[escolher][a];
        }
        return (resp);
    }



    void combinacoes (vector <int> arr, int len, int posicaoInicial, vector <int> aux) {
        if (len == 0) {
            for (int x = 0; x < aux.size(); x++) {
                combinacaoElementos[controle][x] = aux[x];
            }
            controle++;
            return;
        }
        for (int i = posicaoInicial; i <= arr.size() - len; i++) {
            aux[aux.size() - len] = arr[i];
            combinacoes (arr, len-1, i+1, aux);
        }
    }


    void recuperarCaminho (int caminho, vector <int> vet, int c) {
        fprintf(arqSaida, "%d ", caminho);
        vector <int> aux = extrair(vet,caminho);
        if (aux.size() != 0) {
            int custoAux;
            int custoCaminho = infinito;
            int temp = 0;
            for (int a = 0; a < aux.size(); a++) {
                custoAux = custo[c-1]->getVal (aux[a], extrair(aux,aux[a]));

                if (custoAux < custoCaminho) {
                    custoCaminho = custoAux;
                    temp = aux[a];
                }
            }

            recuperarCaminho (temp, aux, (c-1));
        }
    }


};


int main(int argc, char **argv){
    TSPD *novo = new TSPD();
    novo->algoritmoTSPD();
    delete (novo);
    return(0);
}
