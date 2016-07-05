#ifndef _VECTORSPACE_
#define _VECTORSPACE_

#include<vector>
#include<complex>
#include <cmath>

class vectors//矢量运算类
{
private://只能在本类中访问

	double dist, dot;
	int i, j, k, temp;

protected: //可以被子类访问
	//static const int d; //dimension of vector space

public://可以被对象访问

	void vec_cross(double cross[3], double Ain[3], double Bin[3])//实现叉乘Ain cross Bin的成员函数
	{
		for (i = 0; i < 3; i++)
		{
			j = i + 1;
			if (j > 2) { j = j - 3; }
			k = j + 1;
			if (k > 2) { k = k - 3; }
			cross[i] = Ain[j] * Bin[k];
		}
	}
	//重载
	void vec_cross(double cross[3], const double Ain[3], const double Bin[3])//实现叉乘Ain cross Bin的成员函数
	{
		for (i = 0; i < 3; i++)
		{
			j = i + 1;
			if (j > 2) { j = j - 3; }
			k = j + 1;
			if (k > 2) { k = k - 3; }
			cross[i] = Ain[j] * Bin[k];
		}
	}

	double vec_dot(double Ain[3], double Bin[3])//实现点乘Ain dot Bin的成员函数
	{
		dot = 0;
		for (j = 0; j < 3; j++)
		{
			dot += Ain[j] * Bin[j];
		}
		return dot;
	}

	void vec_plus(double plus[3], double Ain[3], double Bin[3])//求出Ain+Bin
	{
		for (j = 0; j < 3; j++)
		{
			plus[j] = Ain[j] + Bin[j];
		}
	}

	void vec_minus(double minus[3], double Ain[3], double Bin[3])//求出Ain-Bin
	{
		for (j = 0; j < 3; j++)
		{
			minus[j] = Ain[j] - Bin[j];
		}
	}
	//重载
	void vec_minus(double minus[3], const double Ain[3], const double Bin[3])//求出Ain-Bin
	{
		for (j = 0; j < 3; j++)
		{
			minus[j] = Ain[j] - Bin[j];
		}
	}

	double vec_distance(double Ain[3], double Bin[3])//求出AB模
	{
		double v_minus[3];
		dist = 0;
		for (j = 0; j < 3; j++)
		{
			v_minus[j] = Ain[j] - Bin[j];
			dist += pow(v_minus[j], 2);
		}
		dist = sqrt(dist);
		return dist;
	}
	//重载
	double vec_distance(const double Ain[3], const double Bin[3])//求出AB模
	{
		double v_minus[3];
		dist = 0;
		for (j = 0; j < 3; j++)
		{
			v_minus[j] = Ain[j] - Bin[j];
			dist += pow(v_minus[j], 2);
		}
		dist = sqrt(dist);
		return dist;
	}
	//求出矢量长度
	double vec_scale(double Ain[3])
	{
		double scale = 0;
		for (j = 0; j < 3; j++)
		{
			scale += pow(Ain[j], 2);
		}
		scale = sqrt(scale);
		return scale;
	}
};

//int vectors::d = 3; //dimension of vector space

#endif