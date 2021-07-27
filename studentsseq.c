/*
Grupo - 01
Alunos:
João Pedro Fidelis Belluzzo - 10716661
Leonardo Cerce Guioto - 10716640
Luis Fernando Costa de Oliveira - 10716532
Rodrigo Augusto Valeretto - 10684792
Thiago Daniel Cagnoni de Pauli - 10716629
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <locale.h>
#include <math.h>

typedef struct{
    int min;
    int max;
    double median; 
    double mean;
    double sd;
}city;

typedef struct{
    int min;
    int max;
    double median; 
    double mean;
    double sd;
    city * cities;
}region;

typedef struct{
    int min;
    int max;
    double median;
    double mean;
    double sd;
    region * regions;
}country;

//Função para calcular mediana de um vetor vec de tamanho size
double median(int * vec, int size){
    if(size%2){
        return (double) vec[size/2];
    }else{
        return ((double)vec[size/2] + (double)vec[size/2-1]) / 2;
    }
}

//Função para calcular a média de um vetor vec de tamanho size
double mean(int * vec, int size){
    double sum = 0;
    for(int i=0; i<size; i++){
        sum += vec[i];
    }
    return sum / size;
}

//Função para calcular desvio padrão de um vetor vec de tamanho size, com a média mean
double sd(int * vec, int size, double mean){
    double sum = 0;
    for(int i=0; i<size; i++){
        sum += pow(vec[i]-mean,2);
    }
    return sqrt(sum / size);
}

//Função de comparação utilizada no quicksort
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main(int argc, char * argv[]){
    FILE * file;
    int R, C, A, seed; //Valores de entrada
    int best_mean_region_index, best_city_index_i, best_city_index_j;
    int *grades;
    double time;
    double best_mean_region, best_mean_city; 
    region * regions;
    country brasil;
    setlocale(LC_ALL, "Portuguese");


    scanf("%d", &R);
    scanf("%d", &C);
    scanf("%d", &A);
    scanf("%d", &seed);

    //Criação dos vetores referentes às regiões
    regions = (region*) malloc(R * sizeof(region));
    for(int i=0; i<R; i++){
        //Criação dos vetores referentes às cidades
        regions[i].cities = (city*) malloc(C * sizeof(city));
    }

    brasil.regions = regions;

    //Inicia a geração das notas com base na seed da entrada
    srand(seed);

    grades = (int*) malloc(R * C * A * sizeof(int));

    for(int i=0; i<R*C*A; i++){
        grades[i] = rand() % 101;
    }

    //Início da contagem do tempo
    time = omp_get_wtime();

    //Realiza a ordenação dos vetores das notas de cada cidade
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            qsort(&(grades[i*(C*A) + j*A]), A, sizeof(int), cmpfunc);
        }
    }

    //Calcula os dados referentes às cidades
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            regions[i].cities[j].min = grades[i*(C*A) + j*A];
            regions[i].cities[j].max = grades[i*(C*A) + j*A + A-1];
            regions[i].cities[j].median = median(&(grades[i*(C*A) + j*A]), A);
            regions[i].cities[j].mean = mean(&(grades[i*(C*A) + j*A]), A);
            regions[i].cities[j].sd = sd(&(grades[i*(C*A) + j*A]), A, regions[i].cities[j].mean);
        }
    }

    //Realiza a ordenação dos vetores das notas de cada região
    for(int i=0; i<R; i++){
        qsort(&(grades[i*(C*A)]), C*A, sizeof(int), cmpfunc);
    }

    //Calcula os dados referentes às regiões
    for(int i=0; i<R; i++){
        double sum = 0;
        regions[i].min = grades[i*C*A];
        regions[i].max = grades[i*C*A + C*A-1];
        regions[i].median = median(&(grades[i*(C*A)]), C*A);
        for(int j=0; j<C; j++){
            sum += regions[i].cities[j].mean;
        }
        regions[i].mean = sum / C;
        regions[i].sd = sd(&(grades[i*(C*A)]), C*A, regions[i].mean);
    }

    //Ordena todas as notas do país
    qsort(grades, R*C*A, sizeof(int), cmpfunc);

    //Calcula os dados correspondentes a todas as notas do país
    brasil.min = grades[0];
    brasil.max = grades[R*C*A-1];
    brasil.median = median(grades, R*C*A);
    double sum = 0;
    for(int i=0; i<R; i++){
        sum += brasil.regions[i].mean;
    }
    brasil.mean = sum / R;
    brasil.sd = sd(grades, R*C*A, brasil.mean);

    //Inicia o cálculo da melhor cidade e melhor região com base na média
    best_mean_region = -1;
    best_mean_city = -1;

    for(int i=0; i<R; i++){
        if(regions[i].mean > best_mean_region){
            best_mean_region = regions[i].mean;
            best_mean_region_index = i;
        }
        for(int j=0; j<C; j++){
            if(regions[i].cities[j].mean > best_mean_city){
                best_mean_city = regions[i].cities[j].mean;
                best_city_index_i = i;
                best_city_index_j = j;
            }
        }
    }

    //Finaliza a contagem do tempo
    time = omp_get_wtime() - time;

    //O código abaixo realiza a impressão dos resultados
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n",
            i, j, regions[i].cities[j].min, regions[i].cities[j].max, regions[i].cities[j].median, regions[i].cities[j].mean,
            regions[i].cities[j].sd);
        }
        printf("\n");
    }

    for(int i=0; i<R; i++){
        printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n",
            i, regions[i].min, regions[i].max, regions[i].median, regions[i].mean, regions[i].sd);
    }

    printf("\nBrasil: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n",
    brasil.min, brasil.max, brasil.median, brasil.mean, brasil.sd);

    printf("\nMelhor região: Região %d\n",best_mean_region_index);

    printf("Melhor cidade: Região %d, Cidade %d\n",best_city_index_i, best_city_index_j);

    printf("\nTempo de resposta sem considerar E/S, em segundos: %.8fs\n",time);

    free(grades);
    for(int i=0; i<R; i++){
        free(regions[i].cities);
    }
    free(regions);
    
    return 0;
}