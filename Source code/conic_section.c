#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"polynomial_operation.h"
#define d 10e-8
#define f 10e10

double line(double x, double m, double c)
{
	double y;
	y = m * x + c;
	return y;
}

void sublc(int n,double h, double k, double r, double x[], double y[])
{
	if ((x[0] > (h + 2.0*r)) || (x[0] < (h - 2.0*r)) || (y[0] > (k + 2.0*r)) || (y[0] < (k - 2.0*r)))
		printf("Not intersect");
	else
		for (int i = 0; i < n; i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
}

void linecircle(double h, double k, double m, double c, double r)
{
	double a, b, c1,c2;
	double x[2],y[2];
	int n;
	a = m*m + 1.0;
	b = 2.0*m*c - 2.0*m*k - 2.0*h;
	c1 = (c-k)*(c-k) + h*h - r*r;
	c2 = (c - h)*(c - h) + k * k - r * r;
	if (fabs(m) > f)
	{
		Quad(1.0, -2.0*k, c2, &n, y);
		for (int i = 0; i < n; i++)
			x[i] = c;
		sublc(n, h, k, r, x, y);
		exit(0);
	}
	else {
		Quad(a, b, c1, &n, x);
		if (n == 0)
		{
			printf("Not intersect\n");
			exit(1);
		}
		for (int i = 0; i < n; i++)
			y[i] = line(x[i], m, c);
		sublc(n, h, k, r, x, y);
	}
}

void subci(double h, double r1, double r2, double k1, double k2, double x[], double *y,int *n)
{
	*y = (r1*r1 - r2 * r2 - k1 * k1 + k2 * k2) / (2.0*k2 - 2.0*k1);
	double s = r1 * r1 - (*y - k1)*(*y - k1);
	if (s < 0)
	{
		printf("Not intersect");
		*n = 0;
		exit(0);
	}
	else if (s == 0) {
		*n = 1;
		x[0] = h;
	}
	else
	{
		*n = 2;
		x[0] = h + sqrt(s);
		x[1] = h - sqrt(s);
	}
}

void circleintersect(double h1, double k1, double h2, double k2, double r1, double r2)
{
	double m, c,x[2],y[2];
	int n;
	if ((h1 == h2) && (k1 == k2) && (r1 != r2))
	{
		printf("Not intersect");
		exit(0);
	}
	else if ((h1 == h2) && (k1 == k2) && (r1 == r2))
	{
		printf("Same circle");
		exit(0);
	}
	else if ((h1 == h2) || (k1 == k2))
	{
		if (h1 == h2)
		{
			subci(h1, r1, r2, k1, k2, x, &y[0], &n);
			for (int i = 0; i < n; i++)
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[0]);
			exit(0);
		}
		else if (k1 == k2)
		{
			subci(k1, r1, r2, h1, h2, y, &x[0], &n);
			for (int i = 0; i < n; i++)
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[0], y[i]);
			exit(0);
		}
	}
	m = -(h1 - h2) / (k1 - k2);
	c = ((r2*r2 - r1 * r1) + (h1*h1 - h2 * h2) + (k1*k1 - k2 * k2)) / (2 * (h1 - h2));
	linecircle(h1, k1, m, c, r1);
}

void lineline(double m1, double c1, double m2, double c2)
{
	double x, y;
	if ((m1 == m2) && (c1 == c2))
	{
		printf("Same line");
		exit(1);
	}
	else if (m1 == m2)
	{
		printf("No intersect point");
		exit(1);
	}
	x = (c2 - c1) / (m1 - m2);
	y = line(x, m1, c1);
	printf("Intersect point is (%.2f,%.2f)\n", x, y);
}

void para1(double b, double c, double d1, double e, double m1, double c1)
{
	double A, B, C;
	double x[2];
	double y[2];
	int n;
	A = b * m1*m1;
	B = 2.0 * b*m1*c1 + c + d1 * m1;
	C = b * c1*c1 + d1 * c1 + e;
	Quad(A, B, C, &n, x);
	if (n == 0)
	{
		printf("No intersect point");
		exit(0);
	}
	else {
		for (int i = 0; i < n; i++)
		{
			y[i] = line(x[i], m1, c1);
			printf("Intersect point %d is (%.2f,%.2f)\n", i+1, x[i], y[i]);
		}
	}
}

void para2(double a, double c, double d1, double e, double m1, double c1)
{
	double A, B, C;
	double x[2];
	double y[2];
	int n;
	A = a;
	B = c+d1*m1;
	C = d1*c1 + e;
	Quad(A, B, C, &n, x);
	if (n == 0)
	{
		printf("No intersect point");
		exit(0);
	}
	else {
		for (int i = 0; i < n; i++)
		{
			y[i] = line(x[i], m1, c1);
			printf("Intersect point %d is (%.2f,%.2f)\n", i+1, x[i], y[i]);
		}
	}
}

void ellipse_hyper(double a, double b, double c, double d1, double e, double m1, double c1)
{
	double A, B, C;
	double x[2], y[2];
	int n;
	A = a + b * m1*m1;
	B = 2.0 * m1*c1*b + c + d1 * m1;
	C = b * c1*c1 + d1 * c1 + e;
	Quad(A, B, C, &n, x);
	if (n == 0)
	{
		printf("No intersect point");
		exit(0);
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			y[i] = line(x[i], m1, c1);
			printf("Intersect point %d is (%.2f,%.2f)\n", i+1, x[i], y[i]);
		}
	}
}

void lineconic(double a, double b, double c, double d1, double e, double m1, double c1)
{
	if (fabs(a)<d)
		para1(b, c, d1, e, m1, c1);
	else if (fabs(b)<d)
		para2(a, c, d1, e, m1, c1);
	else
		ellipse_hyper(a, b, c, d1, e, m1, c1);
}

void subcir_para(double a, double c, double d1, double e, double h, double k, double r,double x[],double y[],int *n)
{
	if (d1 == 0)
	{
		printf("Error no parabola");
		exit(EXIT_FAILURE);
	}
	double A, B, C, D, E;
	A = a * a;
	B = 2.0*a*c;
	C = d1 * d1 + c * c + 2.0*a*(e + d1 * k);
	D = 2.0*c*(e + d1 * k) - 2.0*h*d1*d1;
	E = d1 * d1*h*h + (e + d1 * k)*(e + d1 * k) - d1 * d1*r*r;
	quartic(A, B, C, D, E, n, x);
	if (*n == 0)
	{
		printf("No intersection");
		exit(0);
	}
	else {
		for (int i = 0; i < *n; i++)
			y[i] = (-a * x[i] * x[i] - c * x[i] - e) / (d1);
	}
}

void subcir_conic(double a, double b, double c, double d1, double e, double k1, double k2, double k3)
{
	double A, B, C, D, E;
	double x[4], y[4];
	int n;
	A = a * k1*k1;
	B = 2.0*k1*k2*a;
	C = a * k2*k2 + 2.0*k3*k1*a + b + c * k1;
	D = 2.0*k3*k2*a + c * k2 + d1;
	E = a * k3*k3 + c * k3 + e;
	quartic(A, B, C, D, E, &n, y);
	if (n == 0)
	{
		printf("No intersection");
		exit(0);
	}
	else {
		for (int i = 0; i < n; i++)
		{
			x[i] = k1 * y[i] * y[i] + k2 * y[i] + k3;
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
		}
	}
}

void subcir_conic2(double a, double b, double c, double d1, double e, double h, double k, double r, int *n, double x[], double y[])
{
	double k1, k2, k3;
	k1 = a - b;
	k2 = -a * k*2.0 - d1;
	k3 = a * h*h + a*k * k - e - a * r*r;
	Quad(k1, k2, k3, n, y);
	if (*n == 0)
	{
		printf("No solution");
		exit(EXIT_FAILURE);
	}
	else if (*n == 1)
	{
		x[0] = h + sqrt(r*r - (y[0] - k)*(y[0] - k));
		x[1] = h - sqrt(r*r - (y[0] - k)*(y[0] - k));
		if(fabs(x[0]-x[1])<d)
			printf("Intersect point 1 is (%.2f,%.2f)\n", x[0], y[0]);
		else
			for	(int i = 0; i < *n; i++)
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[0]);
	}
	else
	{
		double x1[2], x2[2];
		int j1, j2;
		Quad(1.0, -2.0*h, h*h - r * r + (y[0] - k)*(y[0] - k), &j1, x1);
		Quad(1.0, -2.0*h, h*h - r * r + (y[1] - k)*(y[1] - k), &j2, x2);
		for (int i = 0; i < j1 + j2; i++)
		{
			if (i < j1)
				printf("Intersect point %d is (%.2f,%.2f)\n",i+1, x1[i], y[0]);
			else
				printf("Intersect point %d is (%.2f,%.2f)\n",i+1, x2[i-j1], y[1]);
		}
	}
}

void cir_conic(double a, double b, double c, double d1, double e, double h, double k, double r)
{
	double x[4], y[4];
	int n;
	if (fabs(a) < d)
	{
		subcir_para(b, d1, c, e, k, h, r, y, x, &n);
		for(int i=0;i<n;i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
	else if (fabs(b) < d)
	{
		subcir_para(a, c, d1, e, h, k, r, x, y, &n);
		for(int i=0;i<n;i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
	else
	{
		if (fabs(2.0*a*h + c) < d)
		{
			subcir_conic2(a, b, c, d1, e, h, k, r, &n,x, y);
			exit(EXIT_SUCCESS);
		}
		double v1, v2, v3;
		v1 = (a - b) / (2.0*a*h + c);
		v2 = (-2.0*a*k - d1) / (2.0*a*h + c);
		v3 = (a*h*h + k * k - e - a * r*r) / (2.0*a*h + c);
		subcir_conic(a, b, c, d1, e, v1, v2, v3);
	}
}

void subparapara(double b, double c, double d1, double e, double k1, double k2, double k3,double x[],double y[],int *n)
{
	double A, B, C;
	A = b + c * k1;
	B = d1 + c * k2;
	C = e + c * k3;
	Quad(A,B,C, n, y);
	if (*n == 0)
	{
		printf("No intersect");
		exit(0);
	}
	else {
		for (int i = 0; i < *n; i++)
			x[i] = k1 * y[i] * y[i] + k2 * y[i] + k3;
	}
}

void subparapara2(double a, double c, double d1, double e, double k1, double k2, double k3, double x[], double y[], int *n)
{
	double A, B, C, D, E;
	A = a * k1*k1;
	B = 2.0*a*k1*k2;
	C = a * k2*k2 + 2.0*k1*k3*a + c * k1;
	D = 2.0*k2*k3*a + c * k2 + d1;
	E = a * k3*k3 + c * k3 + e;
	quartic(A,B,C,D,E, n, y);
	if (*n == 0)
	{
		printf("No intersect");
		exit(0);
	}
	else {
		for (int i = 0; i < *n; i++)
			x[i] = k1 * y[i] * y[i] + k2 * y[i] + k3;
	}
}

void parapara(double a1, double b1, double c1, double d1, double e1, double a2, double b2, double c2, double d2, double e2)
{
	if ((fabs(a1) < d) && (fabs(a2) < d))
	{
		if (fabs(c1) < d)
		{
			printf("No solution");
			exit(1);
		}
		double k1, k2, k3;
		int n;
		k1 = -b1 / c1;
		k2 = -d1 / c1;
		k3 = -e1 / c1;
		double x[2], y[2];
		subparapara(b2, c2, d2, e2, k1, k2, k3, x, y, &n);
		for(int i=0;i<n;i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
	else if ((fabs(b1) < d) && (fabs(b2) < d))
	{
		if (fabs(d1) < d)
		{
			printf("No solution");
			exit(1);
		}
		double k1, k2, k3;
		int n;
		k1 = -a1 / d1;
		k2 = -c1 / d1;
		k3 = -e1 / d1;
		double x[2], y[2];
		subparapara(a2, d2, c2, e2, k1, k2, k3, y, x, &n);
		for(int i=0;i<n;i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
	else if ((fabs(a1) < d) && (fabs(b2) < d))
	{
		if (fabs(c1) < d)
		{
			printf("No solution");
			exit(1);
		}
		double k1, k2, k3;
		int n;
		k1 = -b1 / c1;
		k2 = -d1 / c1;
		k3 = -e1 / c1;
		double x[4], y[4];
		subparapara2(a2, c2, d2, e2, k1, k2, k3, x, y, &n);
		for (int i = 0; i < n; i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
	else if ((fabs(b1) < d) && (fabs(a2) < d))
	{
		if (fabs(d1) < d)
		{
			printf("No solution");
			exit(1);
		}
		double k1, k2, k3;
		int n;
		k1 = -a1 / d1;
		k2 = -c1 / d1;
		k3 = -e1 / d1;
		double x[4], y[4];
		subparapara2(b2, c2, d2, e2, k1, k2, k3, y, x, &n);
		for (int i = 0; i < n; i++)
			printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
	}
}


void subparaconic(double a,double b, double c, double d1, double e, double k1, double k2, double k3, double x[], double y[], int *n)
{
	double A, B, C, D, E;
	A = a * k1*k1;
	B = 2.0*a*k1*k2;
	C = b + a * k2*k2 + 2.0*k1*k3*a + c * k1;
	D = 2.0*k2*k3*a + c * k2 + d1;
	E = a * k3*k3 + c * k3 + e;
	quartic(A,B,C,D,E, n, y);
	if (*n == 0)
	{
		printf("No intersect");
		exit(0);
	}
	else {
		for (int i = 0; i < *n; i++)
			x[i] = k1 * y[i] * y[i] + k2 * y[i] + k3;
	}
}


void paraconic(double a1, double b1, double c1, double d1, double e1, double a2, double b2, double c2, double d2, double e2)
{
	double k1, k2, k3;
	int n;
	double x[4], y[4];
	if (fabs(a1) < d)
	{
		if (fabs(c1) < d)
		{
			printf("No solution");
			exit(1);
		}
		k1 = -b1 / c1;
		k2 = -d1 / c1;
		k3 = -e1 / c1;
		subparaconic(a2, b2, c2, d2, e2, k1, k2, k3, x, y, &n);
	}
	else if (fabs(b1) < d)
	{
		if (fabs(d1) < d)
		{
			printf("No solution");
			exit(1);
		}
		k1 = -a1 / d1;
		k2 = -c1 / d1;
		k3 = -e1 / d1;
		subparaconic(b2, a2, d2, c2, e2, k1, k2, k3, y, x, &n);
	}
	for (int i = 0; i < n; i++)
		printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
}


void subellipseconic(double a1, double b1, double c1, double d1, double e1, double a2, double b2, double c2, double d2, double e2)
{
	double k1, k2, k3;
	int n;
	double x[4], y[4];
	k1 = a1 * b2 - a2 * b1;
	k2 = a1 * d2 - a2 * d1;
	k3 = a1 * e2 - a2 * e1;
	Quad(k1, k2, k3, &n, y);
	if (n == 0)
	{
		printf("No solution");
		exit(EXIT_FAILURE);
	}
	else if (n == 1)
	{
		int j;
		Quad(a1, c1, b1*y[0] * y[0] + d1 * y[0] + e1, &j, x);
		if (fabs(x[0] - x[1]) < d)
			printf("Intersect point 1 is (%.2f,%.2f)\n", x[0], y[0]);
		else
			for (int i = 0; i < j; i++)
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[0]);
	}
	else
	{
		double x1[2], x2[2];
		int j1, j2;
		Quad(a1, c1, b1*y[0] * y[0] + d1 * y[0] + e1, &j1, x1);
		Quad(a1, c1, b1*y[1] * y[1] + d1 * y[1] + e1, &j2, x2);
		for (int i = 0; i < j1 + j2; i++)
		{
			if (i < j1)
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x1[i], y[0]);
			else
				printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x2[i - j1], y[1]);
		}
	}

}

void ellipse_hyperbola(double a1, double b1, double c1, double d1, double e1, double a2, double b2, double c2, double d2, double e2)
{
	if (fabs(a1*c2 - a2 * c1) < d)
	{
		printf("No solution");
		exit(1);
	}
	double k1, k2, k3;
	int n;
	double x[4], y[4];
	k1 = (a2*b1 - a1 * b2) / (a1*c2 - a2 * c1);
	k2 = (a2*d1 - a1 * d2) / (a1*c2 - a2 * c1);
	k3 = (a2*e1 - a1 * e2) / (a1*c2 - a2 * c1);
	subparaconic(a2, b2, c2, d2, e2, k1, k2, k3, x, y, &n);
	for (int i = 0; i < n; i++)
		printf("Intersect point %d is (%.2f,%.2f)\n", i + 1, x[i], y[i]);
}

