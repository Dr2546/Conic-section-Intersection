#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define m 10e-18
#define m1 10e-8


void Quad(double a, double b, double c, int *n, double x[])
{
	double D = b * b - 4.0*a*c;

	if (D < 0.0) {
		*n = 0;
		return;
	}
	else {
		double d = sqrt(D);
		double A = 2.0*a;
		if (d < 1.0e-14) {
			*n = 1;
			x[0] = -b / A;
		}
		else {			
			*n = 2;
			x[0] = (-b - d) / A;
			x[1] = (-b + d) / A;
		}
	}
}

void cubic(double a, double b, double c, double d,int *n,double x[])
{
	double f = 1.0, f1, x1;
	x1 = 0;
	while (fabs(f) >= m)
	{
		f = a * x1*x1*x1 + b * x1*x1 + c * x1 + d;
		f1 = 3.0*a*x1*x1 + 2.0*b*x1 + c;
		x1 = x1 - f / f1;
	}
	double c1[2];
	c1[0] = a*x1+b;
	c1[1] = -d / x1;
	double y[2];
	Quad(a, c1[0], c1[1], n, y);
	if (0 == n)
		x[0] = x1;
	else if (((fabs(x1 - y[0]) <= m1) || (fabs(x1 - y[1]) <= m1)))
		for (int i = 0; i < *n; i++)
			x[i] = y[i];
	else 
	{
		int i;
		x[0] = x1;
		for (i = 0; i < *n; i++)
			x[i+1] = y[i];
		*n = *n + 1;
	}
}

void sturm1(double b, double c,int *n, double x[])
{
	double p[4];
	int v[2];
	int d;
	int n1;
	int i = 0;
	double dum = (256.0*c*c*c) / (27.0*b*b*b) - b;
	while (i < 2)
	{
		n1 = 0;
		if (i == 0)
		{
			p[0], p[1] = -1.0;
			if (b > m)
				p[2] = 1.0;
			else
				p[2] = -1.0;
		}
		else {
			p[0], p[1] = 1.0;
			if (b > m)
				p[2] = -1.0;
			else
				p[2] = 1.0;
		}
		if (dum > m)
			p[3] = 1.0;
		else
			p[3] = -1.0;
		for (int j = 0; j < 3; j++)
			if (p[j] * p[j + 1] == -1.0)
				n1++;
		v[i] = n1;
		i++;
	}
	d = v[0] - v[1];
	if (abs(d) == 0)
	{
		*n = 0;
		return;
	}
	else
	{
		double f = 1.0, f1, x1;
		x1 = 0;
		while (fabs(f) >= m)
		{
			f =  x1*x1*x1*x1 + b * x1 + c;
			f1 = 4.0*x1*x1*x1+b;
			x1 = x1 - f / f1;
		}
		x[0] = x1;
		double y[3];
		cubic(1.0, x1, x1*x1, b + x1 * x1*x1, n, y);
		for (i = 0; i < *n; i++)
			x[i + 1] = y[i];
		*n = *n + 1;
	}
}

void sturm2(double a,double b, double c,int *n, double x[])
{
	double p[5];
	int v[2];
	int n1, d;
	int i = 0;
	double alpha = (8.0*c) / a - 2.0*a - (9.0*b*b) / (a*a);
	double beta = -(12.0*b*c) / (a*a) - b;
	double dum = (-1.0*beta) / (alpha);
	double p4 = (a*dum*dum) / 2.0 + (3.0*b*dum) / 4.0 - c;
	while (i < 2)
	{
		n1 = 0;
		p[0] = 1.0;
		if (i == 0)
		{
			p[1] = -1.0;
			if (a > m)
				p[2] = -1.0;
			else 
				p[2] = 1.0;
			if (alpha > m)
				p[3] = -1.0;
			else
				p[3] = 1.0;
			
		}else
		{
			p[1] = 1.0;
			if (a > m)
				p[2] = -1.0;
			else
				p[2] = 1.0;
			if (alpha > m)
				p[3] = 1.0;
			else
				p[3] = -1.0;
		}
		if (p4 > m)
			p[4] = 1.0;
		else
			p[4] = -1.0;
		for (int j = 0; j < 4; j++)
			if (p[j] * p[j + 1] == -1.0)
				n1++;
		v[i] = n1;
		i++;
	}
	d = v[0] - v[1];
	if (abs(d) == 0)
	{
		*n = 0;
		return;
	}
	else
	{
		double f = 1.0, f1, x1;
		x1 = 0;
		while (fabs(f) >= m)
		{
			f = x1 * x1*x1*x1 + a * x1*x1 + b * x1 + c;
			f1 = 4.0*x1*x1*x1 + 2.0*a*x1 + b;
			x1 = x1 - f / f1;
		}
		x[0] = x1;
		double y[3];
		cubic(1.0, x1, a + x1 * x1, -c / x1, n, y);
		for (i = 0; i < *n; i++)
			x[i + 1] = y[i];
		*n = *n + 1;
	}
}

void quartic(double a, double b, double c, double d, double e,int *n,double x[])
{
	int i;
	double y[4];
	double A = c / a - (3.0*b*b) / (8.0*a*a);
	double B = (b*b*b) / (8.0*a*a*a) - (b*c) / (2.0*a*a) + d / a;
	double C = (b*b*c) / (16.0*a*a*a) - (3.0*b*b*b*b) / (256.0*a*a*a*a) - (b*d) / (4.0*a*a) + e / a;
	if ((A == 0) && (B == 0))
	{
		if (C > m)
		{
			*n = 0;
			return;
		}
		else if (C <= m)
		{
			*n = 2;
			y[0] = pow(C, 0.25);
			y[1] = -pow(C, 0.25);
			for (i = 0; i < *n; i++)
				x[i] = y[i] - b / (4.0*a);
		}
	}
	else if (A == 0)
	{
		sturm1(B, C, n,y);
		for (i = 0; i < *n; i++)
			x[i] = y[i] - b / (4.0*a);
	}
	else
	{
		sturm2(A, B, C,n, y);
		for (i = 0; i < *n; i++)
			x[i] = y[i] - b / (4.0*a);
	}

}