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
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include <locale.h>

#define T 8

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

int find_max(int * vec, int tam){
    int maximo = INT_MIN;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(max : maximo)
    for(int i=0; i<tam; i++){
        maximo = (maximo < vec[i]) ? vec[i] : maximo;
    }

    return maximo;
}

int find_min(int * vec, int tam){
    int minimo = INT_MAX;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(min : minimo)
    for(int i=0; i<tam; i++){
        minimo = (minimo > vec[i]) ? vec[i] : minimo;
    }

    return minimo;
}

double find_city_mean(int * vec, int tam){
    double mean = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : mean)
    for(int i=0; i<tam; i++){
        mean += vec[i];
    }

    mean /= tam;

    return mean;
}

double find_mean(double * vec, int tam){
    double mean = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : mean)
    for(int i=0; i<tam; i++){
        mean += vec[i];
    }

    mean /= tam;

    return mean;
}

double find_sd(int * vec, int tam, double mean){
    double sd = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : sd)
    for(int i=0; i<tam; i++){
        sd += pow((vec[i] - mean),2);
    }

    sd /= tam;
    sd = sqrt(sd);

    return sd;
}


//Função de comparação utilizada no quicksort
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main(int argc, char * argv[]){

    int provided;
    setlocale(LC_ALL, "Portuguese");

    // Inicia o MPI indicando a utilização de threads com OMP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    // Se não estiver disponível a utilização de múltiplas threads, encerra o processo
    if(provided != MPI_THREAD_MULTIPLE){
        printf("Problema na criação das threads\n");
        return 1;
    }

    int seed; //Valores de entrada
    int best_mean_region_index, best_city_index_i, best_city_index_j;
    int *grades;
    double time;
    double best_mean_region, best_mean_city; 
    region * regions;
    country brasil;

    int root = 0, dst = 0, tag = 0;
    int my_rank, num_proc;
    int RCA[3];

    omp_set_nested(1);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    
    if(my_rank == 0){
        scanf("%d", &(RCA[0]));
        scanf("%d", &(RCA[1]));
        scanf("%d", &(RCA[2]));
        scanf("%d", &seed);
    }

    MPI_Bcast(RCA, 3, MPI_INT, root, MPI_COMM_WORLD);

    grades = (int*) malloc(RCA[0] * RCA[1] * RCA[2] * sizeof(int));

    if(my_rank == 0){
        //Inicia a geração das notas com base na seed da entrada
        srand(seed);

        for(int i=0; i < RCA[0] * RCA[1] * RCA[2]; i++){
            grades[i] = rand() % 101;
        }
    }

    MPI_Bcast(grades, RCA[0] * RCA[1] * RCA[2], MPI_INT, root, MPI_COMM_WORLD);

    
    //Criação dos vetores referentes às regiões
    regions = (region*) malloc(RCA[0] * sizeof(region));

    for(int i=0; i<RCA[0]; i++){
        //Criação dos vetores referentes às cidades
        regions[i].cities = (city*) malloc(RCA[1] * sizeof(city));
    }

    brasil.regions = regions;

    if(my_rank == 1){
        int * max_cidades = (int*) malloc(RCA[0] * RCA[1] * sizeof(int));
        int * max_regioes = (int*) malloc(RCA[0] * sizeof(int));
        int max_brasil;

        #pragma omp parallel for num_threads(T/2) shared(max_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            max_cidades[i] = find_max(&(grades[i * RCA[2]]), RCA[2]);
            //printf("Cidade: %d | Max: %d\n",i,max_cidades[i]);
        }

        #pragma omp parallel for num_threads(T/2) shared(max_regioes)
        for(int i=0; i<RCA[0]; i++){
            max_regioes[i] = find_max(&(max_cidades[i*RCA[1]]), RCA[1]);
            //printf("Região: %d | Max: %d\n",i,max_regioes[i]);
        }

        max_brasil = find_max(max_regioes, RCA[0]);
        //printf("Brasil | Max: %d\n",max_brasil);

    }

    if(my_rank == 2){
        int * min_cidades = (int*) malloc(RCA[0] * RCA[1] * sizeof(int));
        int * min_regioes = (int*) malloc(RCA[0] * sizeof(int));
        int min_brasil;

        #pragma omp parallel for num_threads(T/2) shared(min_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            min_cidades[i] = find_min(&(grades[i * RCA[2]]), RCA[2]);
            //printf("Cidade: %d | Min: %d\n",i,min_cidades[i]);
        }

        #pragma omp parallel for num_threads(T/2) shared(min_regioes)
        for(int i=0; i<RCA[0]; i++){
            min_regioes[i] = find_min(&(min_cidades[i*RCA[1]]), RCA[1]);
            //printf("Região: %d | Min: %d\n",i,min_regioes[i]);
        }

        min_brasil = find_min(min_regioes, RCA[0]);
        //printf("Brasil | Min: %d\n",min_brasil);

    }

    if(my_rank == 3){
        double * mean_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * mean_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double mean_brasil;
        MPI_Request request;
        int dst = 4, tag = 1;

        #pragma omp parallel for num_threads(T/2) shared(mean_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            mean_cidades[i] = find_city_mean(&(grades[i * RCA[2]]), RCA[2]);
            printf("Cidade: %d | Mean: %.2lf\n",i,mean_cidades[i]);
        }

        MPI_Isend(mean_cidades, RCA[0] * RCA[1], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &request);

        #pragma omp parallel for num_threads(T/2) shared(mean_regioes)
        for(int i=0; i<RCA[0]; i++){
            mean_regioes[i] = find_mean(&(mean_cidades[i*RCA[1]]), RCA[1]);
            printf("Região: %d | Mean: %.2lf\n",i,mean_regioes[i]);
        }

        tag = 2;
        MPI_Isend(mean_regioes, RCA[0], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &request);

        mean_brasil = find_mean(mean_regioes, RCA[0]);
        printf("Brasil | Mean: %.2lf\n",mean_brasil);

        tag = 3;
        MPI_Send(&mean_brasil, 1, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);

    }
    if(my_rank == 4){
        double * mean_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * sd_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * mean_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double * sd_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double mean_brasil;
        double sd_brasil;
        int src = 3, tag = 1;
        MPI_Status status;

        MPI_Recv(mean_cidades, RCA[0] * RCA[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        #pragma omp parallel for num_threads(T/2) shared(sd_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            sd_cidades[i] = find_sd(&(grades[i * RCA[2]]), RCA[2], mean_cidades[i]);
            printf("Cidade: %d | SD: %.2lf\n",i,sd_cidades[i]);
        }

        tag = 2;
        MPI_Recv(mean_regioes, RCA[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        #pragma omp parallel for num_threads(T/2) shared(mean_regioes)
        for(int i=0; i<RCA[0]; i++){
            sd_regioes[i] = find_sd(&(grades[i*RCA[1]*RCA[2]]), RCA[1]*RCA[2], mean_regioes[i]);
            printf("Região: %d | SD: %.2lf\n",i,sd_regioes[i]);
        }

        tag = 3;
        MPI_Recv(&mean_brasil, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        sd_brasil = find_sd(grades, RCA[0] * RCA[1] * RCA[2], mean_brasil);
        printf("Brasil | SD: %.2lf\n",sd_brasil);

    }


    MPI_Finalize();
    return 0;
}
