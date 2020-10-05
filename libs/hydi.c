
#ifndef MAC
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <assert.h>
#include <tgmath.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <pthread.h>
#include <zlib.h>

#include "basic-types.h"
#include "stringutils.h"
#include "fileio.h"
#include "filebuffer.h"
#include "info.h"
#include "mathematics.h"
#include "vstack.h"
#include "manopt.h"
#include "fac.h"

unsigned char mute=0;
typedef unsigned int uint;

//one count structure for each group representing one
//genomic position with bs and oxbs counts respectively 
typedef struct{
  uint16_t samples;
  char *chrom;
  uint64_t pos;
  char strand;
  //flags
  char *f;
  //bs count data
  uint32_t *n;
  uint32_t *k;
  //oxbs count data
  uint32_t *m;
  uint32_t *l;
} hydi_counts_t;

typedef struct{

  double cil1;
  double est1;
  double cir1;
 
  double cil2;
  double est2;
  double cir2;
  
  double cil3;
  double est3;
  double cir3;
  
  long double pval[3];
  long double fdr[3];

} hydi_test_t;

typedef struct {
  uint8_t which;
  hydi_test_t *data;
} hydi_fdrsort_t;

char isOverlapping(long double a, long double b, long double c, long double d) {
  char ovl = (a <= d) && (c <= b) && (a <= b) && (c <= d);
  return(ovl);
}

long double minintervaldist(long double a, long double b, long double c, long double d) {
  
  if(isOverlapping(a,b,c,d)) return 0;
  if(b < c) return c-b;
  if(d < a) return -1*(a-d);

  return 0.0;
}

long double maxintervaldist(long double a, long double b, long double c, long double d) {
  long double left = 0.0;
  long double right = 0.0;
  if(isOverlapping(a,b,c,d)) {
    if(a < c) left = a;
    if(c < a) left = c;
    if(b > d) right =b;
    if(d > b) right =d;
  }
 return(right+left);
}



double pchisq(double x, double df, char lower) {
  if(lower) {
    return gsl_cdf_gamma_P (x, df/2., 2.);
  }
  return gsl_cdf_gamma_Q (x, df/2., 2.); 
}

double pnorm(double x, double sigma, char lower) {
  if(lower) {
    return gsl_cdf_gaussian_P (x, sigma);
  }
  return gsl_cdf_gaussian_Q (x, sigma); 
}



uint64_t hydi_atol(char *str, char* err) {
  uint64_t lnum;
  char *end;

  lnum = strtol(str, &end, 10);  //10 specifies base-10
  if (end == str) {
    *err |= (char)1;
    if(strcmp(str,"NA") && strcmp(str,"na") && strcmp(str,".")) {
      NFO( "Non numeric value '%s' found. Flagging as \"NA\".\n", str);
    }
    return 0;
  }
  return(lnum);
}


char** hydi_header(char *line, uint32_t len, uint32_t* nsamp) {
  uint32_t j;
  stringset_t *set;
  char **names; 

  set = tokensToStringset(NULL, "\t", line, len);   

  if((set->noofstrings-3) % 4 != 0) {
    NFO("%d count columns found in the header. This is not a multiple of 4. Columns missing? Aborting.\n", (set->noofstrings-3));
    exit(-1);
  }

  *nsamp = (set->noofstrings-3)/4; 
  names = ALLOCMEMORY(NULL, NULL, char*, set->noofstrings);

  for(j=0; j < set->noofstrings; j++) {
    names[j] = bl_strdup(set->strings[j].str);
  }

  destructStringset(NULL, set);
  return(names);
}

void hydi_parse(char *line, uint32_t len, hydi_counts_t *tab, 
    uint32_t nsamp, uint64_t nline) {

  uint32_t j,dim;
  stringset_t *set;
  char *str, err;

  set = tokensToStringset(NULL, "\t", line, len);  
  dim = (set->noofstrings-3)/4; 

  if((set->noofstrings-3) % 4 != 0) {
  NFO( "%d count columns found in line %d. This is not a multiple of 4. Columns missing? Aborting.\n", (set->noofstrings-3), nline);
   exit(-1);
  }

  if(dim != nsamp) {
   NFO( "%d samples read in line %d. This does not match the number of samples in the header (n=%d). Aborting.\n", dim, nline, nsamp);
   exit(-1);
  }

  //alloc count matrix
  tab->m = ALLOCMEMORY(NULL, NULL, uint32_t,dim);
  tab->k = ALLOCMEMORY(NULL, NULL, uint32_t,dim);
  tab->n = ALLOCMEMORY(NULL, NULL, uint32_t,dim);
  tab->l = ALLOCMEMORY(NULL, NULL, uint32_t,dim);
  tab->f = ALLOCMEMORY(NULL, NULL, char,dim);

  //chromosome
  str = set->strings[0].str; 
  tab->chrom = bl_strdup(str);
  //position 
  str = set->strings[1].str;
  tab->pos = hydi_atol(str,&err);
  //strand
  str = set->strings[2].str;
  tab->strand = str[0]; 
  //counts by sample
  tab->samples = dim; 
  for(j=0; j < dim; j++){
    
    err = 0;
    
    str = set->strings[(j*4)+3].str;  
    tab->n[j] = hydi_atol(str, &err);
    
    str = set->strings[(j*4)+4].str;  
    tab->k[j] = hydi_atol(str, &err);

    str = set->strings[(j*4)+5].str;  
    tab->m[j] = hydi_atol(str, &err);
    
    str = set->strings[(j*4)+6].str;  
    tab->l[j] = hydi_atol(str, &err);

    //set error flag
    //if(err) fprintf(stderr, "setting error flag for column %d\n", j);

    tab->f[j] = err;
  }

  destructStringset(NULL, set);
}

void hydi_print(hydi_counts_t *tab) {
  uint32_t i;
  printf("%s:%lu (%c)\n", tab->chrom,tab->pos,tab->strand);
  for(i=0; i < tab->samples; i++) {
    printf("excluded sample: %d\n", tab->f[i]);
    printf("%d\t%d\t", tab->n[i], tab->k[i]);
    printf("%d\t%d\n", tab->m[i], tab->l[i]);
  }
}

void hydi_destruct(hydi_counts_t *tab) {
  FREEMEMORY(NULL, tab->chrom);
  FREEMEMORY(NULL, tab->n);
  FREEMEMORY(NULL, tab->k);
  FREEMEMORY(NULL, tab->m);
  FREEMEMORY(NULL, tab->l);
  FREEMEMORY(NULL, tab->f);

  return;
}

void hydi_assert(hydi_counts_t *tab1, hydi_counts_t *tab2, uint64_t line) {
  uint32_t i;

  if(strcmp(tab1->chrom,tab2->chrom)) {
    NFO("Different chromsomes found in line %llu of both input files. Wrong sorting or missing value?\n", line);
    exit(-1);
  }

  if(tab1->pos != tab2->pos){
    NFO("Different genomic positions found in line %llu of both input files. Wrong sorting or missing value? Aborting.\n", line);
    exit(-1);
  }

  if(tab1->strand != tab2->strand){
    NFO("Different genomic strands found in line %llu of both input files. Wrong sorting or missing value? Aborting.\n", line);
    exit(-1);
  }


  for(i=0; i < tab1->samples; i++) {
    if(!tab1->f[i]){
      if(tab1->n[i] < tab1->k[i]) {
        NFO("n < k in sample %d in line %d of group 1. Aborting.\n", i, line);
        exit(-1);
      }
      if(tab1->m[i] < tab1->l[i]) {
        NFO("m < l in sample %d in line %d of group 1. Aborting.\n", i, line);
        exit(-1);
      }
    }
  }
  for(i=0; i < tab2->samples; i++) {
    if(!tab2->f[i]){
      if(tab2->n[i] < tab2->k[i]) {
        NFO("n < k in sample %d in line %d of group 2. Aborting.\n", i, line);
        exit(-1);
      }
      if(tab2->m[i] < tab2->l[i]) {
        NFO("m < l in sample %d in line %d of group 2. Aborting.\n", i, line);
        exit(-1);
      }
    }
  }



  return;
}

long double gmplog(mpz_t x) {
    signed long int ex;
    const long double di = mpz_get_d_2exp(&ex, x);
    return logl(di) + logl(2) * (long double) ex;
}

void gmpfactorial(mpz_t result, int n) {
   mpz_set_ui(result,1);
   for(int i=1; i < n+1; i++)
      mpz_mul_ui(result,result,i);
}

long double binomial(uint k, uint n, long double **F) {
  return (expl(F[n][1]-(F[k][1]+F[n-k][1])));
}

#ifdef MAC
int dbl_order_cmp_decr(void *thunk, const void *a, const void *b){
#else
int dbl_order_cmp_decr(const void *a, const void *b, void *thunk){
#endif
  uint64_t  *first = (uint64_t*) a;
  uint64_t  *secnd = (uint64_t*) b;
  hydi_fdrsort_t *thk = (hydi_fdrsort_t*) thunk; 
  hydi_test_t *dat = thk->data;
  uint8_t which = thk->which;

  //nans and infs are sorted to the top
  char nan1 = isnan(dat[*first].pval[which]) | isinf(dat[*first].pval[which]);
  char nan2 = isnan(dat[*secnd].pval[which]) | isinf(dat[*secnd].pval[which]);
 
  if (nan1 && !nan2) return -1;
  if (nan2 && !nan1) return 1;
  if (nan1 && nan2) return 0;

  if (dat[*first].pval[which] < dat[*secnd].pval[which]) return 1;
  if (dat[*first].pval[which] > dat[*secnd].pval[which]) return -1;

  return 0;
}

void hydi_fdr(hydi_test_t *p, uint64_t n, uint8_t which) {
  uint64_t *order, i, k ;
  const long double one = 1.0;
  long double adj=1, cur; 
  hydi_fdrsort_t data;

  data.which = which;
  data.data = p;

  //init order vector
  order = ALLOCMEMORY(NULL, NULL, uint64_t, n);
  for(i=0; i < n; i++) order[i] = i;

  //sort order vector according too p-values in decreasing manner
#ifdef MAC
  qsort_r(order, n, sizeof(uint64_t), &data, dbl_order_cmp_decr);  
#else
  qsort_r(order, n, sizeof(uint64_t), dbl_order_cmp_decr, &data);  
#endif
  //adjust using BH 
  for(i=0, k=0; i < n; i++) { 
    if(!isnan(p[order[i]].pval[which]) && !isinf(p[order[i]].pval[which])) {
      cur = (long double) (n-k)/((n-k)-(i-k)) * (long double) p[order[i]].pval[which];
      adj = DMIN(adj, cur);
      adj = DMIN(one, adj);
      p[order[i]].fdr[which] = adj; 
#ifdef HYDI_FDR_DEBUG 
      printf("%llu -> %e (adj:%Le (%lld/%lld))\n", order[i], p[order[i]].pval[which], adj, (n-k), (n-k)-(i-k));
#endif 
    } else {
      //skip counter for nan and inf at the top of the ordering
      p[order[i]].fdr[which] = p[order[i]].pval[which];
      k++;
    }  
  } 
 
  FREEMEMORY(NULL,order);

  return;
}

long double dbinoml(uint k, uint n, long double p){  
  if(p ==  .0 && k != 0) return log(0.0);
  if(p ==  .0 && k == 0) return log(1.0);
  if(p == 1.0 && k == n) return log(1.0);
  if(p == 1.0 && k != n) return log(0.0);

  //long double R = F[n][1]-(F[k][1]+F[n-k][1]);
  long double R = fac[n]-(fac[k]+fac[n-k]);
  return(R+logl(p)*k+logl((long double)1.0-p)*(n-k));
}


long double dbinoml_gamma(uint k, uint n, long double p){ 
  if(p ==  .0 && k != 0) return log(0.0);
  if(p ==  .0 && k == 0) return log(1.0);
  if(p == 1.0 && k == n) return log(1.0);
  if(p == 1.0 && k != n) return log(0.0);

  long double R = gsl_sf_lngamma(n+1) - gsl_sf_lngamma(k+1) - gsl_sf_lngamma(n-k+1); 
  return(R+logl(p)*k+logl((long double)1.0-p)*(n-k));
}

long double dbinoml_nofac(uint k, uint n, long double p){
  if(p ==  .0 && k != 0) return log(0.0);
  if(p ==  .0 && k == 0) return log(1.0);
  if(p == 1.0 && k == n) return log(1.0);
  if(p == 1.0 && k != n) return log(0.0);

  return(logl(p)*k+logl((long double)1.0-p)*(n-k));
}

long double dbinoml_wrapper(uint k, uint n, long double p) {
  long double res;
  if(n >= 50000) { //Fm
    res = dbinoml_gamma(k,n,p);
  } else {
    res = dbinoml(k,n,p);
  }
  return(res);
}

long double llbdnl(long double d, uint32_t z1, uint32_t z2, uint32_t n, uint32_t m, long double pi) {
  long double one = 1.0;
  long double t1, t2, t3, t4;

  if(pi == 0) {
    t2 = 0;
  } else {
    t2 = (n-z1)*logl(one-pi);
  }
  if(pi == 1) {
    t1 = 0;
  } else {
    t1 = z1 * logl(pi);
  }
  if(one-pi+d == 0) {
    t4 = 0;
  } else {
    t4 = (m-z2)*logl(one-pi+d);
  }
  if(pi-d == 0) {
    t3 = 0;
  } else {
    t3 = z2*logl(pi-d);
  }
  
  return(t1+t2+t3+t4);
}

long double llbdn_ddl(long double d, uint32_t n, uint32_t x, long double pi) {
  long double xl=x, nl=n, one=1.0;
  return(-xl/(pi-d) + (nl-xl)/(one-pi+d));
}

long double dbdnl(long signed z, uint n, uint m, long double p, long double q) {
  long signed l;
  long double sum = .0, tmp;
  char first = 1;

  if(z >= 0) {
    for(l=0; l <= n; l++) {
#ifdef HYDI_DBNL_DEBUG
      printf("z %ld, l %ld, n %d, m %d\n", z, l, n, m);
#endif
      if(z+l <= n && l <= m) {

#ifdef HYDI_DBNL_DEBUG
        printf("%.52Lf\n",dbinoml(z+l,n,p));
        printf("%.52Lf\n",dbinoml(l,m,q));
#endif
        tmp = dbinoml_wrapper(z+l,n,p) + dbinoml_wrapper(l,m,q);
        if(first) {
          sum = tmp;
          first = 0;
        } else {
          sum = logaddl(sum, tmp);
        }
      }
    }
  } else {
    for(l=0; l <= m; l++) {
#ifdef HYDI_DBNL_DEBUG
      printf("z %ld, l %ld, n %d, m %d\n", z, l, n, m);
#endif

      if(l <= n && l-z <= m && l-z >=0) {
#ifdef HYDI_DBNL_DEBUG
        printf("l:%ld, n:%d, p:%Le, l-z:%ld, m:%d, q:%Le", l, n, p, l-z, m, q);
        printf("%.52Lf\n",dbinoml_gamma(l,n,p));
        printf("%.52Lf\n",dbinoml_gamma(l-z,m,q));
           
        printf("%.52Lf\n",dbinoml_wrapper(l,n,p));
        printf("%.52Lf\n",dbinoml_wrapper(l-z,m,q));

#endif

        tmp = dbinoml_wrapper(l,n,p) + dbinoml_wrapper(l-z,m,q);
        if(first) {
          sum = tmp;
          first = 0;
        } else {
          sum = logaddl(sum, tmp);
        }
      }
    }
  }
  return(sum);
}


long double dbdnl_deriv_pi(uint k, uint l, uint n, uint m, long double pi1, long double d) {
  long double x1 = k, x2=l, n1=n, n2=m;
  return(x1/pi1 - (n1-x1)/(1-pi1) + x2/(pi1-d) - (n2-x2)/(1-pi1+d));
}

long double dbdnl_find_pi(hydi_counts_t *tab,  long double d, char rev, long double eps, uint32_t maxiter) {
  uint i, k=0; 
  long double x, n, grad, dir, l, r, mid;
  long double half = 0.5, one=1.0;

  l= d+eps;
  r = one-eps;

  if(l > one-eps) {
    return(one-eps);
  }

  if(rev) {
    for(x=0, n=0, i=0; i < tab->samples; i++){
      x += tab->l[i];
      n += tab->m[i];
    }
  } else {
    for(x=0, n=0, i=0; i < tab->samples; i++){
      x += tab->k[i];
      n += tab->n[i];
    }
  }

  mid = x/n;
  if(mid < l || mid > r) mid = (l+r)*half;
 

  while(k==0 || (k < maxiter && fabsl(grad) > eps)) {
    x = mid;
    
    if(rev) {
      for(grad=.0, i=0; i < tab->samples; i++) {
        //the reversal here is a bit ugly but better to interpret in the output
        grad += dbdnl_deriv_pi(tab->l[i], tab->k[i], tab->m[i], tab->n[i], x, d);
      }

    } else {
      for(grad=.0, i=0; i < tab->samples; i++) {
        grad += dbdnl_deriv_pi(tab->k[i], tab->l[i], tab->n[i], tab->m[i], x, d);
      }
    }

    dir = (grad > 0) - (grad < 0);
    
    if(dir == 1) {
      l = mid;
      mid = MAX((mid+r)*half,d+eps);
      if(l == mid && mid == d+eps) {
        return(mid);
      }
    } 
    
    if(dir == -1) {
      r = mid;
      mid = MAX((mid+l)*half,d+eps); 
      if(mid == r && mid == d+eps){
        return(mid);
      }
    }

    if(dir == 0) {
      return(x);
    }

    k++;
  }

  return x;
}


long double  dbdnl_optim_pi1(hydi_counts_t *tab, long double d, long double eps, uint32_t maxiter, 
    long double *p1, long double *p2) {
  long double pi1, pi2, sum;
  uint32_t i;

  if(d >= 0) { 
    pi1 = dbdnl_find_pi(tab, d, 0, eps, maxiter); 
    pi2 = MAX(eps,pi1-d);
  } else {
    pi2 = dbdnl_find_pi(tab, -1*d, 1, eps, maxiter); 
    pi1 = MAX(eps,pi2+d);
  }

  for(sum=0, i=0; i < tab->samples; i++) {
    sum += dbinoml_nofac(tab->k[i], tab->n[i], pi1);
    sum += dbinoml_nofac(tab->l[i], tab->m[i], pi2);
  }

  *p1 = pi1;
  *p2 = pi2;
  return sum;
}


long double dbdnl_optim_pi2(hydi_counts_t *tab1, hydi_counts_t *tab2, long double d1, long double d2, 
    long double eps, uint32_t maxiter, long double *p1, long double *p2) {

  long double sum=0;
  long double pi11, pi12, pi21, pi22; 
  sum+= dbdnl_optim_pi1(tab1,d1,eps,maxiter, &pi11, &pi12);
  sum+= dbdnl_optim_pi1(tab2,d2,eps,maxiter, &pi21, &pi22);

  *p1 = pi11;
  *p2 = pi12;
  return sum;
}

long double dbdnl_profile_delta_H0(hydi_counts_t *tab1, hydi_counts_t *tab2, 
    long double eps, uint32_t maxiter, long double *p1, long double *p2){
  
  long double l, r, lt, rt, vl=log(0), vr=log(0);
  long double three = 3.0;
  long double pi11=0, pi12=0, pi21=0, pi22=0;

  l = -1+eps;
  r = 1-eps;
  
  while(fabsl(r-l) >= eps) {
    lt = l + (r-l)/three;
    rt = r - (r-l)/three;
    
    vl = dbdnl_optim_pi2(tab1, tab2, lt, lt, eps, maxiter, &pi11, &pi12);
    vr = dbdnl_optim_pi2(tab1, tab2, rt, rt, eps, maxiter, &pi21, &pi22);

    if(vl < vr) {
      l = lt;
    } else {
      r = rt;
    }
  }
  *p1 = pi11;
  *p2 = pi12;
  return(vl);
}


long double dbdnl_profile_delta_H1(hydi_counts_t *tab, long double eps, uint32_t maxiter,
    long double *p1, long double *p2) {
  long double l, r, lt, rt, vl=log(0), vr=log(0);
  long double three = 3.0;
  long double mone = -1;
  long double one = 1;
  long double pi11=0, pi12=0, pi21, pi22;
  l = mone+eps;
  r = one-eps;
  
  while(fabsl(r-l) >= eps) {
    lt = l + (r-l)/three;
    rt = r - (r-l)/three;
    vl = dbdnl_optim_pi1(tab, lt, eps, maxiter, &pi11, &pi12);
    vr = dbdnl_optim_pi1(tab, rt, eps, maxiter, &pi21, &pi22);
    if(vl < vr) {
      l = lt;
    } else {
      r = rt;
    }
  }

  *p1 = pi11;
  *p2 = pi12;
  return(vl);
}

struct diffscore_deriv_cubic_params{
  long double x0;
  long double n0;
  long double x_;
  long double n_;
  long double delta;
};

double diffscore_deriv_cubic(double x, void *params) {
  long double sum;
  struct diffscore_deriv_cubic_params *p
    = (struct diffscore_deriv_cubic_params *) params;

  sum  = (p->n_) * powl(x,3);
  sum += ((p->n_ + p->n0)*p->delta - (p->x_ + p->n_)) * powl(x,2);
  sum += (p->x_ - (2*p->x0 +p->n_)*p->delta + p->n0*(p->delta*p->delta)) * x;
  sum += (p->delta*(1-p->delta)*p->x0);

  return(sum);
}

void hydi_gsl_handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno) {
  fprintf(stderr, "error handler invoked with %d in %s:%d,  for reason %s\n", gsl_errno, file, line, reason);
}

long double binomprop_diffscore(uint32_t k, uint32_t l, uint32_t n, uint32_t m, long double d) {

  long double mle = ((long double)k)/((long double) n);
  long double scr, var, pi = mle, res=0;
  long double one = 1.0;  
  double x_lo = .0, x_hi= .0;
  int status=0, iter=0, max_iter=100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  struct diffscore_deriv_cubic_params params;
  gsl_error_handler_t* old_handler;
  old_handler = gsl_set_error_handler (&hydi_gsl_handler);

  //find root of cubic function approximating the 
  //derivative of the score function
  params.x0 = k;
  params.n0 = n;
  params.x_ = ((long double)k) + ((long double)l);
  params.n_ = ((long double)n) + ((long double)m);
  params.delta =d;


  double left = diffscore_deriv_cubic(0.0001, &params);
  double right = diffscore_deriv_cubic(0.9999-d, &params);

  //printf("x0: %Lf, n0: %Lf, x_: %Lf, n_: %Lf\n", params.x0, params.n0, params.x_, params.n_);

  if((left < 0 && right >0) || (left >0 && right <0)) { 
    F.function = &diffscore_deriv_cubic;
    F.params = &params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, 0.0001, 0.9999-d);  

    do {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      pi = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 0.0001);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
  } else {
    pi = mle;
  }

  scr = (((long double)l)-((long double)m)*(pi+d))/((pi+d)*(1-(pi+d)));
  var = 1/((pi*(one-pi)/((long double)n))+(pi+d)*(one-(pi+d))/((long double)m));
  res = scr/sqrtl(var);

  if(isnanl(res)) {
    res = 0;
  }

  gsl_set_error_handler(old_handler);
  //printf("mle:%Lf, pi: %Lf, delta:%Lf, score: %Lf, var: %Lf, res:%Lf, k:%d, l:%d, n:%d, m:%d, brend:%d\n", mle, pi, d, scr, var, res, k,l,n,m, brent);
  return res; 
}


long double hydi_cibdnsearch(long double a, long double b, hydi_counts_t *tab, 
    long double p, long double off, int sdir, long double eps, uint32_t maxiter) {
  
  uint32_t k=0;
  long double l, r, mid, grad = 0;
  long double two = 2.0;
  long double pi1, pi2;
  int dir;

  l = a;
  r = b;
  mid = (a+b)/two;

  while(k == 0 || (k < maxiter && fabsl(grad) > eps)) {
    
    grad = off;
    grad +=  dbdnl_optim_pi1(tab, mid, eps, maxiter, &pi1, &pi2); 
        
    dir = (grad > 0) - (grad < 0);

    if(dir == sdir) {
      l = mid;
      mid = (mid+r)/two;
    }

    if(dir == -1*sdir) {
      r = mid;
      mid = (l+mid)/two;
    }

    k += 1;
  }

  return(mid);
}

void hydi_cibdn(hydi_counts_t *tab, long double p, long double d, 
    long double eps, uint32_t maxiter, double chialphahalf, long double *cil, long double *cir) {

  long double mone = -1.0, one=1.0;
  long double a, b, l, r, off = .0;
  uint32_t i;

  a = mone+eps;//(mone+epsd)+p;
  b = one-eps;//p-epsd;

  l = d;
  r = d;

  if(a < b) {
    for(i=0; i < tab->samples; i++) { 
      off +=  llbdnl(d, tab->k[i], tab->l[i], tab->n[i], tab->m[i], p);
    }
    
    if(a < d) {
      l = hydi_cibdnsearch(a, d, tab, p, mone*off+chialphahalf, -1, eps, maxiter);
    } else {
      l = d;
    }
    if(d < b) {
      r = hydi_cibdnsearch(d, b, tab, p, mone*off+chialphahalf, +1, eps, maxiter);
    } else {
      r = d;
    }
    if(fabsl(l) == eps){
      l =0;
    }

    if(fabsl(r) == eps){
      r =0;
    }
  }

  *cir = r;
  *cil = l;

  return;
}

long double hydi_cibdnsearch2(long double a, long double b, hydi_counts_t *tab1, hydi_counts_t *tab2, 
    long double d1, long double d2, long double off, int sdir, long double eps, uint32_t maxiter) {
  
  uint32_t k=0;
  long double l, r, mid, grad = 0;
  long double two = 2.0;
  long double pi1, pi2;
  int dir;

  l = a;
  r = b;
  mid = (a+b)/two;

  while(k == 0 || (k < maxiter && fabsl(grad) > eps)) {
    
    grad = off;
    grad +=  dbdnl_optim_pi2(tab1, tab2, d1-mid, d2-mid, eps, maxiter, &pi1, &pi2); 
        
    dir = (grad > 0) - (grad < 0);

    if(dir == sdir) {
      l = mid;
      mid = (mid+r)/two;
    }

    if(dir == -1*sdir) {
      r = mid;
      mid = (l+mid)/two;
    }

    k += 1;
  }

  return(mid);
}

void hydi_cibdn2(hydi_counts_t *tab1, hydi_counts_t *tab2, long double p1, long double p2, long double d1, long double d2, 
   long double eps, uint32_t maxiter, double chialphahalf, long double *cil, long double *cir) {
  
  long double mone = -1.0, one=1.0, two=2.0, mtwo=-2.0;
  long double a, b, l, r, dd, off = .0;
  uint32_t i;

  a = mone+(MAX(fabsl(d2),fabsl(d1)));//mone+eps;//(mone+epsd)+p;
  b = one-(MAX(fabsl(d2),fabsl(d1))); //one-eps;//p-epsd;
  dd = 0;
  l = dd;
  r = dd;

  if(a < b) {
    for(i=0; i < tab1->samples; i++) {
      off +=  llbdnl(d1, tab1->k[i], tab1->l[i], tab1->n[i], tab1->m[i], p1); 
    }

    for(i=0; i < tab2->samples; i++) {
      off +=  llbdnl(d2, tab2->k[i], tab2->l[i], tab2->n[i], tab2->m[i], p2);
    } 

    if(a < dd) {
      l = hydi_cibdnsearch2(a, dd, tab1, tab2, d1, d2, mone*off+chialphahalf, -1, eps, maxiter);
    } else {
      l = dd;
    }
    if(dd < b) {
    r = hydi_cibdnsearch2(dd, b, tab1, tab2, d1, d2, mone*off+chialphahalf, +1, eps, maxiter);
    } else {
      r = dd;
    }
    if(fabsl(l) == eps){
      l =0;
    }

    if(fabsl(r) == eps){
      r =0;
    }
  }

  *cir = MAX(mtwo,MIN(two,r));
  *cil = MAX(mtwo,MIN(two,l));

  return;
}


void hydi_pool(hydi_counts_t *tab, hydi_counts_t *pool, uint32_t *count) {
  uint32_t k, l, n, m, i, c;
   

  for(k=0, l=0, n=0, m=0, i=0, c=0; i < tab->samples; i++) {
    if(!tab->f[i]) {
      
      k += tab->k[i];
      l += tab->l[i];
      n += tab->n[i];
      m += tab->m[i];

      c++;
    }
  }

  pool->m = ALLOCMEMORY(NULL, NULL, uint32_t,1);
  pool->k = ALLOCMEMORY(NULL, NULL, uint32_t,1);
  pool->n = ALLOCMEMORY(NULL, NULL, uint32_t,1);
  pool->l = ALLOCMEMORY(NULL, NULL, uint32_t,1);
  pool->f = ALLOCMEMORY(NULL, NULL, char,1);

  pool->chrom=bl_strdup(tab->chrom);
  pool->l[0] = l;
  pool->k[0] = k;
  pool->m[0] = m; 
  pool->n[0] = n;
  pool->f[0] = 0;
  pool->pos = tab->pos;
  pool->samples = 1;
  pool->strand = tab->strand; 

  *count = c;

  return;
}

void hydi_lrt(hydi_counts_t *tab1, hydi_counts_t *tab2, 
    uint32_t minrep, long double eps, 
    uint32_t maxiter, double chialphahalf, hydi_test_t *t) {

  long double z1, z2, lambda=0, H00=NAN, H11=NAN, H12=NAN;
  long double pi001, pi002, pi111, pi112, pi121, pi122;
  long double cil, cir;
  uint32_t count1, count2;
  hydi_counts_t pool1, pool2;
  char nestedtest = 0; 
  long double mtwo = -2.0, two=2.0;
 
 
  hydi_pool(tab1, &pool1, &count1);

  if(count1 < minrep) {
    t->pval[0] = NAN;
    pi111 = ((long double)pool1.k[0])/pool1.n[0]; 
    pi112 = ((long double)pool1.l[0])/pool1.m[0];
    t->est1 = (double) (pi111-pi112);
    t->cil1 = t->est1;
    t->cir1 = t->est1;
  } else {
    //score test
    z1 = binomprop_diffscore(pool1.k[0], pool1.l[0], pool1.n[0], pool1.m[0], 0.0);
    t->pval[0] = two*pnorm(fabsl(z1), 1, 0);
    //profile test
    H11 = dbdnl_profile_delta_H1(&pool1, eps, maxiter, &pi111, &pi112); 
    t->est1 = (double)(pi111 - pi112);   
    hydi_cibdn(&pool1, pi111, t->est1, eps, maxiter, chialphahalf, &cil, &cir);
    t->cil1 = (double) cil;
    t->cir1 = (double) cir;   
  }

#ifdef HYDI_LRT_DEBUG
  printf("G1\tpi1: %Lf, pi2: %Lf, cil2:%Lf, cir2:%Lf, z2:%Lf, pval2:%Lf\n",pi111, pi112, t->cil1, t->cir1, z1, t->pval[0]);  
#endif

   hydi_pool(tab2,&pool2, &count2);

  if(count2 < minrep) {
    t->pval[1] = NAN;
    pi121 = ((long double)pool2.k[0])/pool2.n[0]; 
    pi122 = ((long double)pool2.l[0])/pool2.m[0];
    t->est2 = (double)(pi121-pi122);
    t->cil2 = t->est2;
    t->cir2 = t->est2;

  } else {  
    z2 = binomprop_diffscore(pool2.k[0], pool2.l[0], pool2.n[0], pool2.m[0], 0.0);
    t->pval[1] = two*pnorm(fabsl(z2), 1, 0);
    //profile test
    H12 = dbdnl_profile_delta_H1(&pool2, eps, maxiter, &pi121, &pi122);
    t->est2 = (double)(pi121 - pi122);
    hydi_cibdn(&pool2, pi121, t->est2, eps, maxiter, chialphahalf, &cil, &cir);
    t->cil2 = (double) cil;
    t->cir2 = (double) cir;
    
  }

#ifdef HYDI_LRT_DEBUG
  printf("G2\tpi1: %Lf, pi2: %Lf, cil2:%Lf, cir2:%Lf, z2:%Lf, pval2:%Lf\n",pi121, pi122, t->cil2, t->cir2, z2, t->pval[1]);  
#endif


  if(count1 < minrep || count2 < minrep) {
    t->pval[2]=NAN;
    t->est3 = t->est1 - t->est2;
    t->cil3 = t->est3;
    t->cir3 = t->est3;
  } else {

    //only test for differential hydroxymethylation if
    //there is a 
    //- a significant departure from 0 in one of the groups
    //- the departure is positive, i.e. not an overshoot

    if(!nestedtest || (t->pval[0] < 0.05 && t->est1 > .0) || (t->pval[1] < 0.05 && t->est2 >.0) ) {
 
      t->est3 = t->est1 - t->est2;
      H00 = dbdnl_profile_delta_H0(&pool1, &pool2, eps, maxiter, &pi001, &pi002);
           
      if(isnanl(H00) || isnanl(H11) || isnanl(H12) || 
          isinfl(H00) || isinfl(H11) || isinfl(H12) || H00-(H11+H12) > 0) {
        lambda = NAN;
        t->pval[2] = 1.0;
      } else {
        lambda = mtwo * (H00 - (H11+H12)); 
        t->pval[2] = pchisq(lambda, 1.0, 0);
        long double cil3;
        long double cir3;
        hydi_cibdn2(&pool1, &pool2, pi111, pi121, t->est1, t->est2, eps, maxiter, chialphahalf, &cil3, &cir3);
        t->cil3 = (t->est1+cil3)-(t->est2+cir3);
        t->cir3 = (t->est1+cir3)-(t->est2+cil3);
      }
    } else {
      t->pval[2] = NAN;
      t->est3 = t->est1 - t->est2;
      t->cil3 = t->est3;
      t->cir3 = t->est3;

    }
  }

#ifdef HYDI_LRT_DEBUG
  printf("H00:%Lf\n", H00);
  printf("H11:%Lf\n", H11);
  printf("H12:%Lf\n", H12);
  printf("lambda:%Lf\n", lambda);
  printf("p3:%f\n", t->pval[2]);
#endif
      
  hydi_destruct(&pool1);
  hydi_destruct(&pool2);
 
  return;
}


void hydi_estimates_print(long double *p) {
  long double *q;

  q = &p[3];
  printf("p0: %0.56Lf\n p1: %0.56Lf\n p2: %0.56Lf\n",p[0],p[1],p[2]);
  printf("q0: %0.56Lf\n q1: %0.56Lf\n q2: %0.56Lf\n",q[0],q[1],q[2]);
  
  return;
}

void hydi_estimates(hydi_counts_t *tab1, hydi_counts_t *tab2, long double *p) {
  uint32_t i; 
  long double a, b, d, e, f, g, h, c; 
  long double *q;

  
  q = &p[3];

  //getting the ML estimates for p and q
  for(a=.0, b=.0, e=.0, f=.0, i=0; i < tab1->samples; i++) {
    if(!tab1->f[i]) {  
      e  += ((long double)tab1->k[i]);
      a  += ((long double)tab1->n[i]);
      f  += ((long double)tab1->l[i]);
      b  += ((long double)tab1->m[i]); 
      
    }   
  }
  p[1] = e/a;
  q[1] = f/b;


  //getting the ML estimates for p and q
  for(c=.0, d=.0, g=.0, h=.0, i=0; i < tab2->samples; i++) {
    if(!tab2->f[i]) {
      g += ((long double)tab2->k[i]);
      c += ((long double)tab2->n[i]);
      h += ((long double)tab2->l[i]);
      d += ((long double)tab2->m[i]); 
    }
  }
  p[2] = g/c;
  q[2] = h/d;

  p[0] = (e+g)/(a+c);
  q[0] = (f+h)/(b+d);
  
   
  return;
}

/*
 * Filters to facilitate post-processing
 * Overshoot is flagged based on unadjusted
 * p-value test to be more conservative. 
 */

char hydi_overshoot(hydi_test_t *t, double alpha) { 
  char overshoot=0;

  //flag is set if lower bound indicates
  //overshoot and simple p-value
  if (t->cil1 < 0 && t->pval[0] < alpha){
    overshoot += 1;
  }

  //flag is set if lower bound indicates
  //overshoot and simple p-value
  if (t->cil2 < 0 && t->pval[1] < alpha) {
    overshoot += 2;
  }

  //if the overshoot in one of the groups
  //is strong enough to draw the groupwise
  //difference into the negative
  //if (t->cir3 < 0 && overshoot) {
  //  overshoot += 4;
  //}

  return(overshoot);
}

/*
 * Filters to facilitate post-processing
 * Hydroxymeth is flagged using adjusted (!)
 * p-value test to be more conservative
 */

char hydi_hmC(hydi_test_t *t, double alpha) {
  
  char hmC=0;

  //flag set if the lower boundary indicates
  //5hmC and adj. sig.
  if (t->cil1 > 0 && t->pval[0] < alpha){
    hmC += 1;
  }

  //flag set if the lower boundary indicates
  //5hmC and adj. sig.
  if (t->cil2 > 0 && t->pval[1] < alpha) {
    hmC += 2;
  }

  return(hmC);
}


char hydi_isgz(char *fn) {
  uint32_t prefixlen = bl_fileprefixlen(fn);

    if(strncmp(&fn[prefixlen], ".gz", 3) == 0 || 
       strncmp(&fn[prefixlen], ".gzip", 5) == 0 ||
       strncmp(&fn[prefixlen], ".bgz", 4) == 0 || 
       strncmp(&fn[prefixlen], ".bgzip", 6) == 0 ) {
      NFO("guessing gzip from extension of input file '%s'.\n",fn);
      return 1;
    }

    return 0;
}

char hydi_fileopen(circbuffer_t *buf, char *fn) {
  FILE *fp;
  gzFile gzfp;

  if(hydi_isgz(fn)) {
    gzfp = gzopen(fn, "r"); 
    
    if(gzfp == NULL) {
      NFO("gzopen of file '%s' failed: %s. Aborting.\n", fn, strerror(errno));
      exit(EXIT_FAILURE);
    }
    //open buffers
    bl_circBufferInitGz(buf, 100000, gzfp, NULL);

  } else {
    fp = fopen(fn, "r"); 
  
    if(fp == NULL) {
      NFO("fopen of file '%s' failed: %s. Aborting.\n", fn, strerror(errno));
      exit(EXIT_FAILURE);
    }
    //open buffers
    bl_circBufferInit(buf, 100000, fp, NULL);
  }

  return 1;
}

void hydi_fileclose(circbuffer_t *buf) {

  if(buf->gzip) {
    gzclose(buf->gzdev);
  } else {
    fclose(buf->dev);
  }

  return;
}

void hydi_output(FILE *dev, char *chrom, uint64_t pos, char strand,
    hydi_test_t *t, uint32_t overshoot, uint32_t hmC) {

  
  fprintf(dev, "%s\t%ld\t%c\t", chrom, pos, strand);
  fprintf(dev, "%f\t%f\t%f\t",t->cil1, t->est1, t->cir1);
  fprintf(dev, "%Lf\t%Lf\t", t->pval[0], t->fdr[0]);
  fprintf(dev, "%f\t%f\t%f\t",t->cil2, t->est2, t->cir2);
  fprintf(dev, "%Lf\t%Lf\t", t->pval[1], t->fdr[1]);
  fprintf(dev, "%d\t%d\t",overshoot, hmC);
  fprintf(dev, "%f\t%f\t%f\t", t->cil3, t->est3, t->cir3);
  fprintf(dev, "%Lf\t%Lf\t", t->pval[2], t->fdr[2]);
  //courtesy calculation: minimum difference of differences
  if(t->cil3 < 0 && t->cir3 > 0) {
    fprintf(dev, "%f\n", 0.0);
  } else if(fabsl(t->cil3) < fabsl(t->cir3)) {
    fprintf(dev, "%f\n", t->cil3);
  } else {
    fprintf(dev, "%f\n", t->cir3);
  }

  return;
}

void hydi_iter(FILE *dev, char* fn1, char* fn2, uint32_t minrep, 
    long double eps, uint32_t maxiter, double alpha) {
  
  uint32_t len1, len2, nsamp1=0, nsamp2=0;  
  char *str1, *str2, **names1=NULL, **names2=NULL, overshoot,hmC;
  double chialphahalf;
  hydi_test_t test, *p_A1; 
  circbuffer_t buf1, buf2;
  VStack p_A1_; 
  
  chialphahalf = gsl_cdf_chisq_Qinv(alpha, 1);
  chialphahalf *= 0.5;
  fprintf(stderr, "%f\n", chialphahalf);

  bl_vstackInit(&p_A1_, 10000, sizeof(hydi_test_t));
 

  hydi_fileopen(&buf1, fn1);
  hydi_fileopen(&buf2, fn2); 

  uint64_t count =0;
  while ((str1 = bl_circBufferReadLine(&buf1, &len1)) != NULL &&
      (str2 = bl_circBufferReadLine(&buf2, &len2)) != NULL) {
    
    if(count % 50000 == 0) NFO("%ld lines processed\n", count); 

    //here we are skipping the header
    if(count > 0) {
      hydi_counts_t tab1;
      hydi_counts_t tab2;

      hydi_parse(str1, len1, &tab1, nsamp1, count);
      hydi_parse(str2, len2, &tab2, nsamp1, count);
      hydi_assert(&tab1, &tab2, count); 

      hydi_lrt(&tab1, &tab2, minrep, eps, maxiter, chialphahalf, &test);   
      bl_vstackPush(&p_A1_, &test);

      hydi_destruct(&tab1);
      hydi_destruct(&tab2);

    } else {  
      names1 = hydi_header(str1, len1, &nsamp1); 
      names2 = hydi_header(str2, len2, &nsamp2); 
    }
 
    count++;
    assert(strlen(str1) == len1);
    FREEMEMORY(NULL, str1);
    FREEMEMORY(NULL, str2);
  }

  hydi_fileclose(&buf1);
  hydi_fileclose(&buf2);

  bl_circBufferDestruct(&buf1); 
  bl_circBufferDestruct(&buf2); 
  
  p_A1 = ((hydi_test_t*)p_A1_.stackspace);
  hydi_fdr(p_A1, p_A1_.top+1, 0);
  hydi_fdr(p_A1, p_A1_.top+1, 1);
  hydi_fdr(p_A1, p_A1_.top+1, 2);
 

  // output results and estimates
  hydi_fileopen(&buf1, fn1);
  hydi_fileopen(&buf2, fn2); 

  count = 0;

  while ((str1 = bl_circBufferReadLine(&buf1, &len1)) != NULL &&
      (str2 = bl_circBufferReadLine(&buf2, &len2)) != NULL) {

    if(count > 0) {
      hydi_counts_t tab1;
      hydi_counts_t tab2;

      hydi_parse(str1, len1, &tab1, nsamp1, count);
      hydi_parse(str2, len2, &tab2, nsamp1, count);
 
      overshoot = hydi_overshoot(&p_A1[count-1], alpha);
      hmC = hydi_hmC(&p_A1[count-1], alpha);
     
      hydi_output(dev, tab1.chrom, tab1.pos, tab1.strand, &p_A1[count-1], overshoot, hmC);

      hydi_destruct(&tab1);
      hydi_destruct(&tab2);

    }     

    count++;
    assert(strlen(str1) == len1);
    FREEMEMORY(NULL, str1);
    FREEMEMORY(NULL, str2);
  }

  //removing names
  for(uint32_t i=0; i < nsamp1*4+3; i++) {
    FREEMEMORY(NULL, names1[i]);
  }
  for(uint32_t i=0; i < nsamp2*4+3; i++) {
    FREEMEMORY(NULL, names2[i]);
  }

  FREEMEMORY(NULL, names1);
  FREEMEMORY(NULL, names2);  
  bl_vstackDestruct(&p_A1_, NULL);
   
  hydi_fileclose(&buf1);
  hydi_fileclose(&buf2);

  bl_circBufferDestruct(&buf1); 
  bl_circBufferDestruct(&buf2); 

  return;
}

void test() {
    long double one = 1.0;
    long double pi = 1.0;
    long double d = 0.0;
    long double n = 153;
    long double m= 149;
    long double var = 1/((pi*(one-pi)/((long double)n))+(pi+d)*(one-(pi+d))/((long double)m));
    fprintf(stderr, "%Lf\n", var);
    fprintf(stderr, "%Lf\n", var/(-1*var));
    long double res = var/(-1*var);
    fprintf(stderr, "%d", isnanl(res));
}

int main(int argc, char *argv[]) {  
  manopt_arg *unflagged;
  manopt_optionset optset; 
  uint32_t minrep = 3, maxiter=100;
  int32_t factorial = 20000;
  long double eps=0.000001;
  double alpha = 0.05;
  FILE *dev = stdout;
  char *version = "0.0.1", *fn1=NULL, *fn2=NULL, *fn3=NULL;
 
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Testing for the equality of hydroxymethylation from oxbs sequencing\n",
      "hydi is free software under GPL \n   2020 Leibniz Institute on Aging (FLI) ",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de");
  
  manopt_blockseparator(&optset, "INPUT");

  manopt(&optset, REQSTRINGOPT, 0, 'a', "group1", 
      "path/filename of query sequences", "<file>", NULL, &fn1);
 
  manopt(&optset, REQSTRINGOPT, 0, 'b', "group2", 
      "path/filename of query sequences", "<file>", NULL, &fn2);

  /*
  manopt(&optset, REQINTOPT, 0, 'g', "factorial", 
      "generate factorial", NULL, NULL, &factorial); 
  */
 
  manopt_blockseparator(&optset, "OUTPUT");
 
  manopt(&optset, REQSTRINGOPT, 0, 'o', "output", 
      "path/filename of output", "<file>", NULL, &fn3);
 
  manopt_blockseparator(&optset, "CONTROL");
 
 
  manopt(&optset, REQINTOPT, 0, 'm', "maxiter", 
      "maximum iterations in gradient descent", NULL, NULL, &maxiter);

  manopt(&optset, REQDBLOPT, 0, 'p', "alpha", 
      "significance level for flags and confidence intervals", NULL, NULL, &alpha);

  manopt(&optset, REQLDBLOPT, 0, 'e', "epsilon", 
      "precision of numerical approximation", NULL, NULL, &eps);

  /*get unflagged options*/
  unflagged = manopt_getopts(&optset, argc, argv);
   
  if ((!manopt_isset(&optset, 'a', NULL) || !manopt_isset(&optset, 'b', NULL)) &&
      !manopt_isset(&optset, 'g', NULL)){
      manopt_help(&optset, "Input files missing.\n");
    }

  if(unflagged->noofvalues > 1) { 
      manopt_help(&optset, "unknown argument(s)\n");
  }

  if(fn3 == NULL) {
    dev = stdout;
  } else {
    dev = fopen(fn3, "w");
    if(dev == NULL) {
      manopt_help(&optset, "could not open output file. Do you have writing privileges in this directory?\n");
    }
  }

  if(manopt_isset(&optset, 'g', NULL)){   
    mpz_t mfrac;
    mpz_init(mfrac);
    for(uint i = 0; i < factorial; i++) {
       gmpfactorial(mfrac, i);
       printf("%d\t%0.56Lf\n", i, gmplog(mfrac));
       //gmp_printf ("%Zd", mfrac);
    }
    mpz_clear(mfrac);
  } else { 
 
    NFO("reading files %s and %s as input\n", fn1, fn2);
    hydi_iter(dev, fn1, fn2, minrep, eps, maxiter, alpha);
  }

    
  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(NULL, unflagged);
  
  return 0;
}
