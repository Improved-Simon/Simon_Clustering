/*
  Accompanying code to the paper "Improved Differential and Linear Cryptanalysis on Round-Reduced SIMON"

  Computes a lower boud on the probability of a linear approximation for Simon or Simeck

  Change the variants parameters and window size at line 31~35; then change the corresponding input output difference/mask; 
  Finally, change the corresponding dynamic window.

  Compile with gcc -Wall -Wextra -O3 -march=native -fopenmp -lm
 */

#define _GNU_SOURCE
#define _XOPEN_SOURCE 500
#include "uthash.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <time.h>

// Parameters
//#define SIMECK
#define LWIM
#define SIMON
#define WORDSIZE 16
#define PRECISION 13 // Window size (w)
#define ROUNDS 13

#ifdef LWIM
  #define total_round ROUNDS/2+1
#endif

#ifndef LWIM
  #define total_round ROUNDS+1
#endif


// #define PRINT_DISTRIBUTION
#define PRINT_MAX

// mmap-based memory allocation reduces memory usage for Simon
#define USE_MMAP

#if defined(SIMON) == defined(SIMECK)
#error Please only exactly one of SIMON and SIMECK
#endif

#define ROT(x,n) ( (((uint64_t)(x)<<(n))&((1ULL<<WORDSIZE)-1)) | (uint64_t)(x)>>(WORDSIZE-(n)) )

#define ROTR(x,n) ( (((uint64_t)(x)<<(WORDSIZE-(n)))&((1ULL<<WORDSIZE)-1)) | (uint64_t)(x)>>(n) )

//#define ROT(x,n) ( ((uint64_t)(x)<<((n)&31) ) | ( (uint64_t)(x)>>((32-(n))&31) ) )
#ifdef SIMECK
static const int a = 0;
static const int b = 5;
static const int c = 1;
#endif
#ifdef SIMON
static const int a = 1;
static const int b = 8;
static const int c = 2;
#endif

uint64_t f(uint64_t x) {
  return ROT(x,c)^(ROT(x, a)&ROT(x, b));
}

typedef double distribution[1<<PRECISION][1<<PRECISION];
typedef long long int logVec[1<<PRECISION];

typedef struct {
  uint64_t B;
  uint64_t A[WORDSIZE]; // vector space AX+B
  double proba;
} transition_space;

typedef struct 
{
  long long int id;
  int index;
  UT_hash_handle hh;
} uthash_user;

#ifdef USE_MMAP

// mmap-based allocation can overcommit and leave holes in allocation
// This reduces memory usage for Simon
void *alloc_distribution() {
  void *p = mmap(NULL, sizeof(distribution), PROT_READ | PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
  assert(p != MAP_FAILED);
  // NOHUGEPAGE reduces memory usage but reduces performance
  // madvise(p, sizeof(distribution), MADV_NOHUGEPAGE);
  return p;
}

void free_distribution(void *p) {
  munmap(p, sizeof(distribution));
}

#else

void *alloc_distribution() {
  void *p = calloc(1, sizeof(distribution));
  assert(p);
  return p;
}

void free_distribution(void *p) {
  free(p);
}

#endif

int vector2index(uthash_user *in, long long int vec)
{
  uthash_user *user;
  HASH_FIND_INT( in, &vec, user);
  int index;
  if (user == NULL)
  {
    index = -1;
  }
  else
  {
    index = user->index;
  }
  return index;
}

uint64_t index2vector(logVec in, int index)
{
  uint64_t vector = in[index];
  return vector;
}

long long int *getVecFromWindow( int window[])
{
  long long int *vec_all = calloc(1, sizeof(logVec));
  int *window_shape = calloc(PRECISION, sizeof(int));
  for (int i=0;i<PRECISION;i++)
  {
    window_shape[i] = window[i];
  }
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<1<<PRECISION; i++)
  {
    long long int tmp = i;
    long long int vec_tmp = 0;
    int j = 0;
    while(tmp){
      vec_tmp ^= (tmp&0x1) << window_shape[j];
      tmp = tmp >> 1;
      j += 1;
    };
    vec_all[i] = vec_tmp;   
  }
  return vec_all;
}

uthash_user *getHash_from_logVec( logVec in )
{
  uthash_user *vec, *vec_all = NULL;
//#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<1<<PRECISION; i++)
  {
    vec = (uthash_user*)malloc(sizeof(uthash_user));
    if (vec == NULL)
    {
      exit(-1);
    }    

    vec->id = in[i];
    vec->index = i; 
    HASH_ADD_INT(vec_all, id, vec);  
    
  }
  return vec_all;
}


void round_trans(distribution in, distribution out, transition_space *m, logVec rin, uthash_user *rin_hash, logVec lout, uthash_user *lout_hash) {

  long long int nonzero_total = 0;
  int *nonzero = calloc(1<<PRECISION, sizeof(int));

#pragma omp parallel for schedule(dynamic)
  for (int l=0; l<1<<PRECISION; l++) {
      for (int r=0; r<1<<PRECISION; r++) {
        uint64_t diff_rin = index2vector(rin, r);
	// Generate output space with Gray enumeration
	if (in[r][l]) {
	  uint64_t d = m[l].B;//gamma
    nonzero[l]++;
    for (int i=1;; i++)
    {
      int index_lout = vector2index( lout_hash, d^diff_rin); //get the rout mask
      if( index_lout != -1 )
      {    
        out[l][index_lout] += in[r][l]*m[l].proba;
        int z = __builtin_ctz(i);
        if (z<WORDSIZE && m[l].A[z])
        {
          d ^= m[l].A[z];
        }  
        else
        {
          break; 
        }                   
      }
      else
      { 
        int z = __builtin_ctz(i);
        if (z<WORDSIZE && m[l].A[z])
        {
          d ^= m[l].A[z];
        }  
        else
        {
          break;
        }
                    
      }
    }
	}
      }
  }
#pragma omp critical
  for (int i=0;i<1<<PRECISION;i++)
  {
    nonzero_total += nonzero[i]; 
  }
  free(nonzero);
  double nonzeroDiff = (double)nonzero_total/(double)(1ULL<<(PRECISION*2));
  printf("The non-zero diff in all distribution is %.2lf %%\n", nonzeroDiff * 100);

}

double print_LWIM_prob(distribution in, transition_space *m, logVec lin, uthash_user *lin_hash, logVec rin, uthash_user *rin_hash)
{
  double total_prob = 0;
  long long int match_total = 0;
  long long int nonzero_total = 0;

  double *r_prob = calloc(1<<PRECISION, sizeof(double));
  int *match = calloc(1<<PRECISION, sizeof(int));
  int *nonzero = calloc(1<<PRECISION, sizeof(int));

#pragma omp parallel for schedule(dynamic)
  for (int l=0;l<1<<PRECISION;l++)
  {
    for(int r=0;r<1<<PRECISION;r++)
    {
      uint64_t diff_rin = index2vector(rin, r);
      if(in[r][l])// r is l
      {
        double prob = in[r][l];
        nonzero[l]++;
        //printf("The in[%u][%u] = %2f\n", r, l, prob);
        uint64_t d = m[l].B;
        for(int i=1;;i++)
        {
          int index_lout = vector2index(rin_hash, d^diff_rin);//rin=lout,lin=rout
          if(index_lout != -1 && in[index_lout][l])
          {
            match[l]++;
            r_prob[l] += prob * m[l].proba * in[index_lout][l];
            int z = __builtin_ctz(i);
            if (z<WORDSIZE && m[l].A[z])
            {
              d ^= m[l].A[z];
            }
            else
            {
              break;
            }
          
          }
          else
          {
            int z = __builtin_ctz(i);
            if (z<WORDSIZE && m[l].A[z])
            {
              d ^= m[l].A[z];
            }  
            else
            {
              break;
            }            
          }
          
        }

      }
    }

  }
#pragma omp critical
  for (int i=0;i<1<<PRECISION;i++)
  {
    total_prob += r_prob[i];
    match_total += match[i];
    nonzero_total += nonzero[i]; 
  }
  free(r_prob);
  free(match);
  double nonzeroDiff = (double)nonzero_total/(double)(1ULL<<(PRECISION*2));
  double get_matched = (double)match_total/(double)nonzero_total;
  printf("The non-zero diff in all distribution is %.2lf %%\n", nonzeroDiff * 100);
  printf("The point that matched total %.2lf %%\n", get_matched * 100);
  printf("The proba of the total is %lf!!!\n", log2(total_prob));

  return total_prob;
}


transition_space* init_transition(logVec lin) {
  transition_space* matrix = calloc(1<<PRECISION, sizeof(transition_space));
  for (int vec=0; vec < 1<<PRECISION; vec++) {
    // Get beta from vector
    uint64_t delta = lin[vec];
    // Compute vector space
    uint64_t B = f(0)^f(delta);
    uint64_t X[WORDSIZE];
    for (int i=0; i<WORDSIZE; i++) {
      X[i] = f(1ULL<<i)^f(delta^(1ULL<<i))^B;
    }
    // Gaussian reduce
    int p = (WORDSIZE-1);
    for (int i=0; i<WORDSIZE; i++) {
      // Find pivot
      while (p>=0 && (X[i] & (1ULL<<p)) == 0) {
	for (int j=i+1; j<WORDSIZE-1; j++) {
	  if (X[j] & (1ULL<<p)) {
	    X[i] ^= X[j];
	    goto PIVOT;;
	  }
	}
	p--;
      }
    PIVOT:
      if (p<0)
	break;
      // Reduce
      for (int j=i+1; j<WORDSIZE; j++) {
	if (X[j] & (1ULL<<p)) X[j] ^= X[i];
      }
      if (B & (1ULL<<p)) B ^= X[i];
    }

    matrix[vec].proba = 1;
    matrix[vec].B = B;

    int j = 0;
    for (int i=WORDSIZE-1; i>=0; i--){
      if (X[i])	{
	  matrix[vec].proba /= 2;

    matrix[vec].A[j++] = X[i];
      }
    }
  }
  return matrix;
}

void print_max(distribution d, int round, logVec lin, logVec rin)
{
      // Print max 
      if (round == 1) {
	printf ("Input:\n");
      }
      double max = 0;
#pragma omp parallel for reduction(max: max)
      for (int i=0; i < 1<<PRECISION; i++) {
	for (int j=0; j < 1<<PRECISION; j++) {
	  if (d[i][j]>max) {
	    max = d[i][j];      
	  }
	}
      }
      printf ("Max: %f", log2(max) );
#pragma omp parallel for
      for (int i=0; i < 1<<PRECISION; i++) {
	for (int j=0; j < 1<<PRECISION; j++) {
	  if (d[i][j] > max*0.999999) {
      uint64_t diff_l = index2vector(lin, i);
      uint64_t diff_r = index2vector(rin, j);
#pragma omp critical
	    printf (" (%lu,%lu)", diff_l, diff_r);
	  }
	}
      }
      printf ("\n");
}

void print_distribution(distribution d)
{
  // Print distribution
  for (int i=0; i<257; i++) {
    printf ("[%2i]: ", i);
    for (int j=0; j<257; j++) {
      double m = log2(d[i][j]);
      printf (" %f", m);
    }
    printf ("\n");
  }
}

int main() {
  clock_t start = clock();
  //初始窗口应该包含输入差分


  int start_window[2][PRECISION] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},//left
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}//right
  };

//Dynamic window for 13-round simon32: (0,1)->(1,0).
  int window[total_round][PRECISION] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} , //0 
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} , //1 
    {2, 1, 8, 3, 4, 5, 6, 7, 0, 9, 10, 11, 12, 13, 14, 15} , //2 
    {4, 0, 10, 3, 9, 2, 5, 7, 8, 1, 6, 11, 12, 13, 14, 15} , //3 
    {6, 8, 5, 12, 1, 11, 2, 4, 10, 3, 7, 9, 13, 0, 14, 15} , //4 
    {8, 0, 4, 7, 3, 14, 13, 10, 9, 6, 12, 2, 5, 11, 15, 1} , //5 
    {10, 2, 8, 12, 1, 9, 5, 0, 15, 11, 6, 4, 14, 7, 13, 3} , //6 

    {8, 0, 4, 7, 3, 14, 13, 10, 9, 6, 12, 2, 5, 11, 15, 1} , //5 
    {6, 8, 5, 12, 1, 11, 2, 4, 10, 3, 7, 9, 13, 0, 14, 15} , //4 
    {4, 0, 10, 3, 9, 2, 5, 7, 8, 1, 6, 11, 12, 13, 14, 15} , //3 
    {2, 1, 8, 3, 4, 5, 6, 7, 0, 9, 10, 11, 12, 13, 14, 15} , //2 
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} , //1 
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} , //0 
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} , //0 

  };


  //Dynamic window for 17-round SIMON48: (0,1)->(1,0)
  // int window[total_round][PRECISION] = {
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //0 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
  //   {1, 2, 8, 0, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
  //   {0, 4, 3, 10, 9, 16, 2, 5, 11, 13, 7, 8, 12, 6, 14, 15, 1, 17, 18} , //3 tree
  //   {6, 8, 5, 12, 1, 11, 4, 18, 10, 17, 0, 3, 7, 9, 13, 2, 14, 15, 16} , //4 tree
  //   {0, 4, 8, 3, 14, 10, 7, 13, 9, 20, 6, 2, 19, 12, 16, 5, 18, 1, 11} , //5 tree
  //   {2, 10, 8, 12, 1, 5, 16, 9, 15, 11, 22, 4, 14, 21, 18, 17, 3, 20, 7} , //6 tree
  //   {8, 12, 0, 10, 18, 7, 14, 3, 11, 13, 17, 9, 16, 6, 20, 23, 2, 19, 5} , //7 tree
  //   {14, 1, 13, 5, 8, 16, 20, 12, 9, 19, 11, 18, 15, 4, 2, 21, 22, 17, 7} , //8 tree

  //   {8, 12, 0, 10, 18, 7, 14, 3, 11, 13, 17, 9, 16, 6, 20, 23, 2, 19, 5} , //7 tree
  //   {2, 10, 8, 12, 1, 5, 16, 9, 15, 11, 22, 4, 14, 21, 18, 17, 3, 20, 7} , //6 tree
  //   {0, 4, 8, 3, 14, 10, 7, 13, 9, 20, 6, 2, 19, 12, 16, 5, 18, 1, 11} , //5 tree
  //   {6, 8, 5, 12, 1, 11, 4, 18, 10, 17, 0, 3, 7, 9, 13, 2, 14, 15, 16} , //4 tree
  //   {0, 4, 3, 10, 9, 16, 2, 5, 11, 13, 7, 8, 12, 6, 14, 15, 1, 17, 18} , //3 tree
  //   {1, 2, 8, 0, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //0 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //0 tree
  // };

//Dynamic window for 23-round simon64: (11,40)->(40,11).
  // int window[total_round][PRECISION] = {

  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //add1 tree
  //   {1, 2, 5, 6, 8, 12, 0, 3, 4, 7, 9, 10, 11, 13, 14, 15, 16, 17, 18} , //add2 tree
  //   {0, 2, 3, 4, 6, 7, 8, 9, 10, 13, 14, 16, 20, 1, 5, 11, 12, 15, 17} , //add3 tree

  //   {6, 10, 2, 5, 12, 9, 16, 8, 1, 15, 11, 22, 4, 18, 14, 17, 21, 3, 7} , //0 tree     
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree     
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree     
  //   {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree     
  //   {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree    
  //   {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree    
  //   {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree  
  //   {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree  
  //   {14, 13, 9, 20, 12, 1, 5, 8, 16, 19, 11, 18, 15, 26, 21, 4, 25, 22, 17} , //8 tree

  //   {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree  
  //   {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree  
  //   {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree    
  //   {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree    
  //   {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree     
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree     
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
  //   {6, 10, 2, 5, 12, 9, 16, 8, 1, 15, 11, 22, 4, 18, 14, 17, 21, 3, 7} , //0 tree
    
  //   {0, 2, 3, 4, 6, 7, 8, 9, 10, 13, 14, 16, 20, 1, 5, 11, 12, 15, 17} , //add3 tree
  //   {1, 2, 5, 6, 8, 12, 0, 3, 4, 7, 9, 10, 11, 13, 14, 15, 16, 17, 18} , //add2 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //add1 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //add0 tree
  //   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //add0 tree  
  // };

//Dynamic window for 33-round simon96: (0,1)->(1,0).
// int window[total_round][PRECISION] = {
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //0 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
//     {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree
//     {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree
//     {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree
//     {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree
//     {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree
//     {14, 13, 9, 20, 12, 1, 5, 8, 16, 19, 11, 18, 15, 26, 21, 4, 25, 22, 17} , //8 tree

//     {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree
//     {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree
//     {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree
//     {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree
//     {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //0 tree
// //*********************************************************************************
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
//     {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree
//     {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree
//     {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree
//     {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree
//     {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree
//     {14, 13, 9, 20, 12, 1, 5, 8, 16, 19, 11, 18, 15, 26, 21, 4, 25, 22, 17} , //8 tree

//     {0, 8, 12, 14, 10, 7, 11, 18, 3, 17, 13, 9, 16, 6, 20, 24, 23, 19, 2} , //7 tree
//     {2, 10, 16, 5, 1, 12, 8, 9, 15, 11, 22, 21, 14, 4, 18, 17, 20, 7, 28} , //6 tree
//     {8, 0, 4, 14, 7, 3, 10, 13, 9, 20, 6, 19, 16, 2, 12, 26, 5, 18, 25} , //5 tree
//     {6, 1, 12, 5, 8, 11, 4, 18, 17, 10, 24, 3, 0, 13, 14, 15, 16, 9, 2} , //4 tree
//     {0, 4, 3, 10, 9, 16, 2, 7, 8, 1, 5, 11, 12, 13, 14, 15, 6, 17, 18} , //3 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //2 tree
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} , //1 tree

//     {2, 6, 5, 12, 1, 8, 11, 4, 18, 0, 3, 7, 9, 10, 13, 14, 15, 16, 17} , //0 tree
//     {2, 6, 5, 12, 1, 8, 11, 4, 18, 0, 3, 7, 9, 10, 13, 14, 15, 16, 17} , //0 tree
  // };




  long long int *lin = getVecFromWindow(start_window[0]);

  uthash_user *lin_hash = getHash_from_logVec(lin);
  long long int *rin = getVecFromWindow(start_window[1]);
  uthash_user *rin_hash = getHash_from_logVec(rin);
  transition_space *m = init_transition(lin);


  double (*d)[1<<PRECISION] = alloc_distribution();

  //Change the input difference here
  int lindex = vector2index(lin_hash, 0); //difference on left part
  int rindex = vector2index(rin_hash, 1); //right part

  if(lindex != -1 && rindex != -1)
  {
    d[rindex][lindex] = 1; // Initial state
  }
  else
  {
    printf("Difference is not in the window[0] and window[1], please reset!!\n");
    exit(0);
  }

  printf("lindex = %x, rindex = %d\n", lindex, rindex);

  // Input mask (1,0) -- rotated by w-1 bits
  //d[1<<(PRECISION-1)][0] = 1; // Initial state
  // d[1<<(PRECISION-2)][1<<(PRECISION-1)] = 1;

  for (int round=1; round<total_round; round++) {
#ifdef PRINT_DISTRIBUTION
    {
      print_distribution(d);
    }
#endif
#ifdef PRINT_MAX
    print_max(d, round, lin, rin);
#endif
#ifdef debug_m
    print_debug(m, round);
#endif
    fflush(stdout);
    long long int *lout = getVecFromWindow( window[round] );
    uthash_user *lout_hash = getHash_from_logVec(lout);
    
    // printf("window[%i] = {", round);
    // for(int i=0;i<PRECISION;i++)
    // {
    //   printf("%i,", window[round][i]);
    // }
    // printf("}\n");
    long long int *rout = lin;
    uthash_user *rout_hash = lin_hash;

    double (*tmp)[1<<PRECISION] = alloc_distribution();
    round_trans(d, tmp, m, rin, rin_hash, lout, lout_hash);
    //更新下一轮输入为上一轮输出
    lin = lout;
    lin_hash = lout_hash;
    rin = rout;
    rin_hash = rout_hash;

    free_distribution(d);
    free(m);
    m = NULL;
    d = tmp;
    if(round < total_round)
    {   
      m = init_transition(lin);
    }   
    // Output diff (1,0), change output difference here
    int index_lout = vector2index(lout_hash, 1);//left part
    int index_rout = vector2index(rout_hash, 0);//right part
    if (index_lout != -1 && index_rout != -1)
    {
      printf("round: %2i, proba: %2f\n", round, log2(d[index_rout][index_lout]));
    }
    else
    {
      printf("round: %2i, not in the window\n", round);
    }

    // Output mask (0,1)
    //printf("%2i: %f\n", i, log2(tmp[0][1<<(PRECISION-1)]));
    // printf("%2i: %f\n", i, log2(t[1<<(PRECISION-1)][1<<(PRECISION-2)]));
    fflush(stdout);



    if(round == total_round-1)
    {
    #ifdef PRINT_DISTRIBUTION
        {
          print_distribution(d);
        }
    #endif
    #ifdef PRINT_MAX
        print_max(d, total_round, lin, rin);
    #endif
    }

    #ifdef LWIM
      if(round == total_round-1)
      {
        double proba = print_LWIM_prob(d, m, lin, lin_hash, rin, rin_hash);
      }
    #endif

  }
  clock_t end = clock();
  printf("clock time: %ld tic\n", end - start );
	printf("clock time : %f ms\n", (double)(end - start) / 1000);
	printf("clock time : %f s\n", (double)(end- start) / CLOCKS_PER_SEC);
}
