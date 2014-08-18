#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "energy.h"
#include "utils.h"
#include "global.h"
#include "constants.h"

long double *V; 
long double *W; 
long double *VBI; 
long double *VM; 
long double **WM; 
long double **WMPrime; 
int *indx; 
long double **PP; 

int alloc_flag = 0;
const float RT = (0.00198721 * 310.15); //* 100.00);
const float RT_ = (0.00198721 * 310.15);

void create_tables(int len) {	
  V = (long double *) malloc(((len+1)*len/2 + 1) * sizeof(long double));
  if (V == NULL) {
    perror("Cannot allocate variable 'V'");
    exit(-1);
  }

  int i;
  WM = (long double **) malloc((len+1)* sizeof(long double *));
  if (WM == NULL) {
    perror("Cannot allocate variable 'WM'");
    exit(-1);
  }   
  for (i = 0; i <= len; i++) {
    WM[i] = (long double *)malloc((len+1)* sizeof(long double));
    if (WM[i] == NULL) {
      perror("Cannot allocate variable 'WM[i]'");
      exit(-1);
    }   
  }

  PP = (long double **) malloc((len+1)* sizeof(long double *));
  if (PP == NULL) {
    perror("Cannot allocate variable 'WM'");
    exit(-1);
  }   
  for (i = 0; i <= len; i++) {
    PP[i] = (long double *)malloc((len+1)* sizeof(long double));
    if (PP[i] == NULL) {
      perror("Cannot allocate variable 'WM[i]'");
      exit(-1);
    }   
  }

  WMPrime = (long double **) malloc((len+1)* sizeof(long double *));
  if (WMPrime == NULL) {
    perror("Cannot allocate variable 'WM'");
    exit(-1);
  }   
  for (i = 0; i <= len; i++) {
    WMPrime[i] = (long double *)malloc((len+1)* sizeof(long double));
    if (WMPrime[i] == NULL) {
      perror("Cannot allocate variable 'WM[i]'");
      exit(-1);
    }   
  }

  VM = (long double *) malloc(((len+1)*len/2 + 1) * sizeof(long double));
  if (VM == NULL) {
    perror("Cannot allocate variable 'V'");
    exit(-1);
  }

  VBI = (long double *) malloc(((len+1)*len/2 + 1) * sizeof(long double));
  if (VBI == NULL) {
    perror("Cannot allocate variable 'V'");
    exit(-1);
  }

  W = (long double *) malloc((len+1) * sizeof(long double));
  if (W == NULL) {
    perror("Cannot allocate variable 'W'");
    exit(-1);
  }

  indx = (int*) malloc((len+1) * sizeof(int));
  if (indx == NULL) {
    perror("Cannot allocate variable 'indx'");
    exit(-1);
  }

  alloc_flag = 1;
	
  init_tables(len);
}


void init_tables(int len) {
  int i, j, LLL;
	
  for (i = 0; i <= len; i++) {
    W[i] = INFINITY; 
    for (j = 0; j <= len; j++) {
      WM[i][j] = INFINITY;
      WMPrime[i][j] = INFINITY;
      PP[i][j] = 0;
    }
  }
	
  LLL = (len)*(len+1)/2 + 1;
  for (i = 0; i < LLL; i++) {
    V[i] = INFINITY;
    VM[i] = INFINITY;
    VBI[i] = INFINITY;
  }

  for (i = 1; i <= (unsigned) len; i++) 
    indx[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */

  return;
}

void free_tables(int len) {
  if (alloc_flag == 1) {
    free(indx);

    int i;
    for (i = 0; i <= len; i++) free(WM[i]);
    free(WM);

    for (i = 0; i <= len; i++) free(WMPrime[i]);
    free(WMPrime);

    free(VM);
    free(VBI);
    free(V);
    free(W);
  }
}


long double Ed3(int i, int j, int k) { return dangle[RNA[i]][RNA[j]][RNA[k]][1];}
long double Ed5(int i, int j, int k) { return dangle[RNA[i]][RNA[j]][RNA[k]][0]; }
long double auPenalty(int i, int j) { return auPen(RNA[i], RNA[j]);}

long double eL1(int i, int j, int ip, int jp) {
  long double energy;
  int size1, size2, size;
  long double loginc; /* SH: Originally unassiged, but needs to be set to 0 so it doesn't throw off later calculations. */
  int lopsided; /* define the asymmetry of an interior loop */

  energy = INFINITY;
  loginc = 0;

  /*SH: These calculations used to incorrectly be within the bulge loop code, moved out here. */
  size1 = ip - i - 1;
  size2 = j - jp - 1;
  size = size1 + size2;

  if (size1 == 0 || size2 == 0) {
    if (size > 30) {
      /* AM: Does not depend upon i and j and ip and jp - Stacking Energies */
      loginc = prelog * log((long double) size / 30.0);
      energy = bulge[30] + eparam[2] + loginc + auPen(RNA[i], RNA[j])
        + auPen(RNA[ip], RNA[jp]);
    } else if (size <= 30 && size != 1) {
      /* Does not depend upon i and j and ip and jp - Stacking Energies  */
      energy = bulge[size] + eparam[2];
      energy += auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]);
    } else if (size == 1) {
      energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[ip], RNA[jp])]
        + bulge[size] + eparam[2];
    }
  } else {
    /* Internal loop */
    lopsided = abs(size1 - size2);

    if (size > 30) {
      loginc = prelog * log((long double) size / 30.0);

      if (!((size1 == 1 || size2 == 1) && gail)) { /* normal internal loop with size > 30*/

        energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + 
          tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[30] + loginc +
          eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
      } else { /* if size is more than 30 and it is a grossely asymmetric internal loop and gail is not zero*/
        energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A,
                                                                                            BASE_A)] + inter[30] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
      }
    }
    else if (size1 == 2 && size2 == 2) { /* 2x2 internal loop */
      energy = iloop22[RNA[i]][RNA[ip]][RNA[j]][RNA[jp]][RNA[i + 1]][RNA[i + 2]][RNA[j - 1]][RNA[j - 2]];
    } else if (size1 == 1 && size2 == 2) {
      energy = iloop21[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]][RNA[j - 2]][RNA[ip]][RNA[jp]];
    } else if (size1 == 2 && size2 == 1) { /* 1x2 internal loop */
      energy = iloop21[RNA[jp]][RNA[ip]][RNA[j - 1]][RNA[i + 2]][RNA[i + 1]][RNA[j]][RNA[i]];
    } else if (size == 2) { /* 1*1 internal loops */
      energy = iloop11[RNA[i]][RNA[i + 1]][RNA[ip]][RNA[j]][RNA[j - 1]][RNA[jp]];
    } else if ((size1 == 2 && size2 == 3) || (size1 == 3 && size2 == 2)) {
      return tstacki23[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]] +
        tstacki23[RNA[jp]][RNA[ip]][RNA[jp + 1]][RNA[ip - 1]];
    }
    else if ((size1 == 1 || size2 == 1) && gail) { /* gail = (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
      energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)]
        + inter[size] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
    } else { /* General Internal loops */
      energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[size]
        + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))])); 
		
    }
  }

  return energy;
}

long double eL(int i, int j, int ip, int jp) {
  long double energy;
  int size1, size2, size;
  long double loginc; /* SH: Originally unassiged, but needs to be set to 0 so it doesn't throw off later calculations. */
  int lopsided; /* define the asymmetry of an interior loop */

  energy = INFINITY;
  loginc = 0;

  /*SH: These calculations used to incorrectly be within the bulge loop code, moved out here. */
  size1 = ip - i - 1;
  size2 = j - jp - 1;
  size = size1 + size2;

  if (size1 == 0 || size2 == 0) {
    if (size > 30) {
      /* AM: Does not depend upon i and j and ip and jp - Stacking Energies */
      loginc = prelog * log((long double) size / 30.0);
      energy = bulge[30] + eparam[2] + loginc + auPen(RNA[i], RNA[j])
        + auPen(RNA[ip], RNA[jp]);
    } else if (size <= 30 && size != 1) {
      /* Does not depend upon i and j and ip and jp - Stacking Energies  */
      energy = bulge[size] + eparam[2];
      energy += auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]);
    } else if (size == 1) {
      energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[ip], RNA[jp])]
        + bulge[size] + eparam[2];
    }
  } else {
    /* Internal loop */
    lopsided = abs(size1 - size2);

    if (size > 30) {
      loginc = prelog * log((long double) size / 30.0);

      if (!((size1 == 1 || size2 == 1) && gail)) { /* normal internal loop with size > 30*/

        energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + 
          tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[30] + loginc +
          eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
      } else { /* if size is more than 30 and it is a grossely asymmetric internal loop and gail is not zero*/
        energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A,
                                                                                            BASE_A)] + inter[30] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
      }
    }
    else if (size1 == 2 && size2 == 2) { /* 2x2 internal loop */
      energy = iloop22[RNA[i]][RNA[ip]][RNA[j]][RNA[jp]][RNA[i + 1]][RNA[i + 2]][RNA[j - 1]][RNA[j - 2]];
    } else if (size1 == 1 && size2 == 2) {
      energy = iloop21[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]][RNA[j - 2]][RNA[ip]][RNA[jp]];
    } else if (size1 == 2 && size2 == 1) { /* 1x2 internal loop */
      energy = iloop21[RNA[jp]][RNA[ip]][RNA[j - 1]][RNA[i + 2]][RNA[i + 1]][RNA[j]][RNA[i]];
    } else if (size == 2) { /* 1*1 internal loops */
      energy = iloop11[RNA[i]][RNA[i + 1]][RNA[ip]][RNA[j]][RNA[j - 1]][RNA[jp]];
    } 
    //else if ((size1 == 2 && size2 == 3) || (size1 == 3 && size2 == 2)) {
    //	return tstacki23[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]] +
    //		tstacki23[RNA[jp]][RNA[ip]][RNA[jp + 1]][RNA[ip - 1]];
    //}
    else if ((size1 == 1 || size2 == 1) && gail) { /* gail = (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
      energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)]
        + inter[size] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
    } else { /* General Internal loops */
      energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[size]
        + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))])); 
		
    }
  }

  return energy;
}

long double eH(int i, int j) {
  /*  Hairpin loop for all the bases between i and j */
  /*  size for size of the loop, energy is the result, loginc is for the extrapolation for loops bigger than 30 */
  int size;
  long double loginc;
  long double energy = INFINITY;
  int key, index, count, kmult;
  long double tlink;

  size = j - i - 1; /*  size is the number of bases in the loop, when the closing pair is excluded */

  /*  look in hairpin, and be careful that there is only 30 values */

  if (size > 30) {
    loginc = prelog * log(((long double) size) / 30.0);
    energy = hairpin[30] + loginc + tstkh[fourBaseIndex(RNA[i], RNA[j],
                                                        RNA[i + 1], RNA[j - 1])] + eparam[4]; /* size penalty + terminal mismatch stacking energy*/
  }

  else if (size <= 30 && size > 4) {
    energy = hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
                                                 RNA[i + 1], RNA[j - 1])] + eparam[4]; /* size penalty + terminal mismatch stacking energy*/
  }

  else if (size == 4) {
    /*  tetraloop */
    key = 0;
    tlink = 0;
    for (index = 0; index < 6; ++index) {
      switch (RNA[i + index]) {
      case BASE_A:
        kmult = 1;
        break;
      case BASE_C:
        kmult = 2;
        break;
      case BASE_G:
        kmult = 3;
        break;
      case BASE_U:
        kmult = 4;
        break;
      default:
        kmult = 0;
        fprintf(stderr, "ERROR: in tetraloop calculation\n");
      }
      key += kmult * (int) pow(10.0, 5 - index);
    }
    /*  if the sequence is in tloop, we use this value */
    for (count = 1; count < numoftloops && tlink == 0; ++count) {
      if (key == tloop[count][0]) {
        tlink = tloop[count][1];
      }
    }
    energy = tlink + hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
                                                         RNA[i + 1], RNA[j - 1])] + eparam[4];
  }

  else if (size == 3) {
    /*  triloop... For the moment, the file triloop.dat is empty */
    /*  else, should have a treatment like the one if size==4 */
    energy = hairpin[size];
    /* AM: Don't include stacking energy terms for triloopls */
    /* + tstkh[RNA[i]][RNA[j]][RNA[i+1]][RNA[j-1]]  */
    /* + eparam[4]; */
    /*  Must be another penalty for terminal AU... Not sure of this */
    energy += auPen(RNA[i], RNA[j]);
  }

  else if (size < 3 && size != 0) {
    /*  no terminal mismatch */
    energy = hairpin[size] + eparam[4];
  } else if (size == 0)
    return INFINITY;

  /*  GGG Bonus => GU closure preceded by GG */
  /*  i-2 = i-1 = i = G, and j = U; i < j */
  if (i > 2) {
    if (RNA[i - 2] == BASE_G && RNA[i - 1] == BASE_G && RNA[i] == BASE_G
        && RNA[j] == BASE_U) {
      energy += gubonus;
      /*  printf ("\n GGG bonus for i %d j %d ", i, j); */
    }
  }

  /*  Poly-C loop => How many C are needed for being a poly-C loop */
  tlink = 1;
  for (index = 1; (index <= size) && (tlink == 1); ++index) {
    if (RNA[i + index] != BASE_C)
      tlink = 0;
  }
  if (tlink == 1) {
    if (size == 3) {
      energy += c3;
    } else {
      energy += cint + size * cslope;
    }
  }

  return energy;
}

long double eS(int i, int j) {
  long double energy;
  /*  not sure about eparam[1], come from MFold.. = 0 */
  energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[i+1], RNA[j-1])] + eparam[1];

  return energy;
}

long double Estackm(int i, int j) {
  return tstackm[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]];
}

long double Estacke(int i, int j) {
  return tstacke[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]];
}
