#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdio>

double f(double x) 
{
	return 2 * cos(2 * x);
}

double F(double x) 
{
	return sin(2 * x);
}

double AnalyticIntegral(double a, double b) 
{
	return F(b) - F(a);
}

double x(double a, double b, double t) {
	return (a + b) / 2 + (b - a) * t / 2;
}

double NumericIntegralNewtonCotes1(double a, double b) 
{
	double I = 0;
	double h = b - a;
	I += f(x(a, b, 0));
	return h * I;
}

double NumericIntegralNewtonCotes2(double a, double b)
{
	double I = 0;
	double h = b - a;
	I += f(x(a, b, -1)) / 2;
	I += f(x(a, b,  1)) / 2;
	return h * I;
}

double NumericIntegralNewtonCotes3(double a, double b)
{
	double I = 0;
	double h = (b - a) / 2;
	I += 1.0 / 3.0 * f(x(a, b, -1));
	I += 4.0 / 3.0 * f(x(a, b,  0));
	I += 1.0 / 3.0 * f(x(a, b,  1));
	return h * I;
}

double NumericIntegralNewtonCotes4(double a, double b)
{
	double I = 0;
	double h = (b - a) / 3;
	I += 3.0 / 8.0 * f(x(a, b, -1.0));
	I += 9.0 / 8.0 * f(x(a, b, -1.0 / 3.0));
	I += 9.0 / 8.0 * f(x(a, b, +1.0 / 3.0));
	I += 3.0 / 8.0 * f(x(a, b, +1.0));
	return h * I;
}

double NumericIntegralNewtonCotes5(double a, double b)
{
	double I = 0;
	double h = (b - a) / 4;
	I += 14.0 / 45.0 * f(x(a, b, -1.0));
	I += 64.0 / 45.0 * f(x(a, b, -0.5));
	I += 24.0 / 45.0 * f(x(a, b, 0));
	I += 64.0 / 45.0 * f(x(a, b, 0.5));
	I += 14.0 / 45.0 * f(x(a, b, 1.0));
	return h * I;
}

double NumericIntegralGauss(double a, double b)
{
	double I = 0;
	double h = (b - a) / 2;
	I += 5.0 / 9.0 * f(x(a, b, -sqrt(3.0 / 5.0)));
	I += 8.0 / 9.0 * f(x(a, b, 0));
	I += 5.0 / 9.0 * f(x(a, b, sqrt(3.0 / 5.0)));
	return h * I;
}

int main()
{
	double a = 0, b = 3.14 / 4;
	const int K = 100;
	double h = (b - a) / K;
	
	double mesh[K + 1];
	double IA = 0;
	double I1K = 0, I2K = 0, I3K = 0, I5K = 0, IGK = 0;
	for (int i = 0; i < K; i++)
	{
		I1K += NumericIntegralNewtonCotes1(a + i * h, a + (i + 1) * h);
		I2K += NumericIntegralNewtonCotes2(a + i * h, a + (i + 1) * h);
		I3K += NumericIntegralNewtonCotes3(a + i * h, a + (i + 1) * h);
		I5K += NumericIntegralNewtonCotes5(a + i * h, a + (i + 1) * h);
		IGK += NumericIntegralGauss(a + i * h, a + (i + 1) * h);
	}
	double I1K2 = 0, I2K2 = 0, I3K2 = 0, I5K2 = 0, IGK2 = 0;
	for (int i = 0; i < 2 * K; i++)
	{
		I1K2 += NumericIntegralNewtonCotes1(a + i * h / 2, a + (i + 1) * h / 2);
		I2K2 += NumericIntegralNewtonCotes2(a + i * h / 2, a + (i + 1) * h / 2);
		I3K2 += NumericIntegralNewtonCotes3(a + i * h / 2, a + (i + 1) * h / 2);
		I5K2 += NumericIntegralNewtonCotes5(a + i * h / 2, a + (i + 1) * h / 2);
		IGK2 += NumericIntegralGauss(a + i * h / 2, a + (i + 1) * h / 2);
	}
	IA = AnalyticIntegral(a, b);
	
	FILE* ResudialFile;
	ResudialFile = fopen("Resudial.txt","w");
	fprintf(ResudialFile, "|------------------------------------------------|\n");
	fprintf(ResudialFile, "|              |                K                |\n");
	fprintf(ResudialFile, "|--------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|   Nodes conut| %14d | %14d |\n",K, 2 * K);
	fprintf(ResudialFile, "|  Newton-Cotes|----------------|----------------|\n");
	fprintf(ResudialFile, "| %12d | %14.8e | %14.8e |\n",1, abs(I1K - IA) / abs(IA), abs(I1K2 - IA) / abs(IA));
	fprintf(ResudialFile, "| %12d | %14.8e | %14.8e |\n",2, abs(I2K - IA) / abs(IA), abs(I2K2 - IA) / abs(IA));
	fprintf(ResudialFile, "| %12d | %14.8e | %14.8e |\n",3, abs(I3K - IA) / abs(IA), abs(I3K2 - IA) / abs(IA));
	fprintf(ResudialFile, "| %12d | %14.8e | %14.8e |\n",5, abs(I5K - IA) / abs(IA), abs(I5K2 - IA) / abs(IA));
	fprintf(ResudialFile, "|         Gauss|----------------|----------------|\n");
	fprintf(ResudialFile, "| %12d | %14.8e | %14.8e |\n",3, abs(IGK - IA) / abs(IA),abs(IGK2 - IA) / abs(IA));
	fclose(ResudialFile);
	
	return 0;
}