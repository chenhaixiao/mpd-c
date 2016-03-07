#include "mem.h"
#include "dbg.h"

char *cvector(int nl, int nh)
{
	char *v;

	v = (char *) malloc((unsigned) (nh - nl + 1) * sizeof(char));
	if (!v)
		log_err("allocation failure in cvector()");
	return v - nl;
}

uchar *ucvector(int nl, int nh)
{
	uchar *v;

	v = (uchar *) malloc((unsigned) (nh - nl + 1) * sizeof(uchar));
	if (!v)
		log_err("allocation failure in cvector()");
	return v - nl;
}

int *ivector(int nl, int nh)
{
	int *v;

	v = (int *) malloc((unsigned) (nh - nl + 1) * sizeof(int));
	if (!v)
		log_err("allocation failure in ivector()");
	return v - nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v = (double *) malloc((unsigned) (nh - nl + 1) * sizeof(double));
	if (!v)
		log_err("allocation failure in dvector()");
	return v - nl;
}


int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, **m;

	m = (int **) malloc((unsigned) (nrh - nrl + 1) * sizeof(int *));
	if (!m)
		log_err("allocation failure 1 in imatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (int *) malloc((unsigned) (nch - ncl + 1) * sizeof(int));
		if (!m[i])
			log_err("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int        i;

	for (i = nrh; i >= nrl; i--)
		free((void *) (m[i] + ncl));
	free((void *) (m + nrl));
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int    i;
	double **m;

	m = (double **) malloc((unsigned) (nrh - nrl + 1) * sizeof(double *));
	if (!m)
		log_err("allocation failure 1 in dmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (double *) malloc((unsigned) (nch - ncl + 1) * sizeof(double));
		if (!m[i])
			log_err("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((void *) (m[i] + ncl));
	free((void *) (m + nrl));
}

char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
	int  i;
	char **m;

	m = (char **) malloc((unsigned) (nrh - nrl + 1) * sizeof(char *));
	if (!m)
		log_err("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (char *) malloc((unsigned) (nch - ncl + 1) * sizeof(char));
		if (!m[i])
			log_err("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}


uchar **ucmatrix(int nrl, int nrh, int ncl, int nch)
{
	int   i;
	uchar **m;

	m = (uchar **) malloc((unsigned) (nrh - nrl + 1) * sizeof(uchar *));
	if (!m)
		log_err("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (uchar *) malloc((unsigned) (nch - ncl + 1) * sizeof(uchar));
		if (!m[i])
			log_err("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}



void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((void *) (m[i] + ncl));
	free((void *) (m + nrl));
}


void free_ucmatrix(uchar ** m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((void *) (m[i] + ncl));
	free((void *) (m + nrl));
}

void free_cvector(char *v, int nl, int nh)
{
	free((void *) (v + nl));
}

void free_ucvector(uchar * v, int nl, int nh)
{
	free((void *) (v + nl));
}


void free_ivector(int *v, int nl, int nh)
{
	free((void *) (v + nl));
}

void free_dvector(double *v, int nl, int nh)
{
	free((void *) (v + nl));
}

