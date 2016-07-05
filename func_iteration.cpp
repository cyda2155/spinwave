#include"headfile_othercpp.h"
#include <omp.h>

void func_iteration(const int Fe_N_Atoms, const int O_N_Atoms, const int N_Atom_total, const int atoms_in_supercell, const double *d_bond, const double lamda_so, double *kkk, const double *Electric_field, const Atom *unit_cell, const Atom *supercell, double *Eigenvalue_out_temp, complex<double> **eigenvectors_out)
{
	complex<double> *Hamiltonian = new complex<double>[N_Atom_total * 2 * N_Atom_total * 2];
	MKL_Complex16  *HAMILTONIAN_MKL = new MKL_Complex16[N_Atom_total * 2 * N_Atom_total * 2];    //2表示自旋,带自旋的2*2块的大矩阵

	double *Eigenvalue = new double[N_Atom_total * 2];
	for (int i = 0; i < N_Atom_total * 2; i++)
	{
		Eigenvalue[i] = 0;
	}


	//考虑迭代自洽的n
	double n[N_Atom_total * 2];

	complex<double> **eigenvectors = new complex<double> *[N_Atom_total * 2]; //只考虑k_band本征矢的各分量
	for (int i = 0; i < N_Atom_total * 2; i++) //行数
	{
		eigenvectors[i] = new complex<double>[N_Atom_total * 2];//new列数
	}

	for (int i = 0; i < N_Atom_total; i++){ n[i] = 1.0; }//各能带迭代前初始化
	for (int i = N_Atom_total; i < N_Atom_total * 2; i++){ n[i] = 0.0; }

	double n0[N_Atom_total * 2], n1[N_Atom_total * 2]; //存储前后给出的电荷密度
	double NDIFFRATIO = 10;
	int count = 1;
	int reset = 0;
	do
	{
		double NDIFF = 0;
		double NTOTAL = 0;
		//每次迭代前必须重置哈密顿量，因为内含+=
		for (int i = 0; i < N_Atom_total * 2; i++)
		{
			for (int j = 0; j < N_Atom_total * 2; j++)
			{
				Hamiltonian[i *  N_Atom_total * 2 + j].real() = 0.0;
				Hamiltonian[i *  N_Atom_total * 2 + j].imag() = 0.0;
			}
		}
		//临时存储同一位点两次循环的n差值的最大值
		double nd_temp = 0;

		//将上一循环n赋给n0
		for (int i = 0; i < N_Atom_total * 2; i++){ n0[i] = n[i]; }

		func_Hamiltonian(N_Atom_total, atoms_in_supercell, lamda_so, Hamiltonian, d_bond, kkk, Electric_field, n, unit_cell, supercell);

		//哈密顿量转换到MKL库对应的矩阵类型
		for (int i = 0; i < N_Atom_total * 2; i++)
		{
			for (int j = 0; j < N_Atom_total * 2; j++)
			{
				HAMILTONIAN_MKL[i * N_Atom_total * 2 + j].real = Hamiltonian[i * N_Atom_total * 2 + j].real();
				HAMILTONIAN_MKL[i * N_Atom_total * 2 + j].imag = Hamiltonian[i * N_Atom_total * 2 + j].imag();
			}
		}
		//对角化，注意HAMILTONIAN_MKL是一维的
		heev_cpp(HAMILTONIAN_MKL, Eigenvalue, N_Atom_total * 2);
		//得到本征矢HAMILTONIAN_MKL和本征值Eigenvalue

		//取出本征矢
		for (int j = 0; j < N_Atom_total * 2; j++)
		{
			for (int i = 0; i < N_Atom_total * 2; i++)
			{
				eigenvectors[i][j].real() = HAMILTONIAN_MKL[i * N_Atom_total * 2 + j].real;
				eigenvectors[i][j].imag() = HAMILTONIAN_MKL[i * N_Atom_total * 2 + j].imag;
			}
		}

		//计算新的电子密度
		for (int i = 0; i < N_Atom_total * 2; i++)
		{
			n[i] = 0;
		}
		for (int i = 0; i < N_Atom_total; i++)
		{
			for (int j = 0; j < Fe_N_Atoms / 2 + O_N_Atoms; j++) //spin up/down各自算，O对应每种都有一个在O site， Fe只有一种在一个Fe site
			{
				n[i] += pow(eigenvectors[i][j].real(), 2) + pow(eigenvectors[i][j].imag(), 2);
				n[i + N_Atom_total] += pow(eigenvectors[i + N_Atom_total][j].real(), 2) + pow(eigenvectors[i + N_Atom_total][j].imag(), 2);
			}
		}

		//比较前后两次计算的电子密度
		for (int i = 0; i < N_Atom_total * 2; i++)
		{
			//将本次循环得到的n赋给n1
			n1[i] = n[i];
			NTOTAL += n1[i];
			NDIFF += n1[i] - n0[i];
		}
		NDIFFRATIO = NDIFF / NTOTAL;
		count += 1;
		if (count > N_loop_times)
		{
			if (count < N_loop_terminate)
			{
				if (NDIFFRATIO < nTor){ break; }
				else if (NDIFFRATIO < nTor_terminate){ break; }
				else { continue; }
			}
			else
			{
				cout << "Too many iterations, n diff ratio=" << NDIFFRATIO << ", loop times=" << count << endl;
				break;
			}
		}
	} while (NDIFFRATIO > nTor);//达到精度后对n的自洽结束

	/*ofstream H("Hamiltonian.dat", ios::app);
	H << "k_distance=" << kkk << endl;
	H << "H:" << endl;*/
	/*cout << "kpoint:" << "(" << ii << ", " << jj << ")" << "  ";
	//cout << "density difference=" << nd_maximum << "  ";
	cout << endl;*///多线程同时cout可能导致冲突
	for (int j = 0; j < N_Atom_total * 2; j++)
	{
		for (int i = 0; i < N_Atom_total * 2; i++)
		{
			eigenvectors_out[i][j] = eigenvectors[i][j];
			//H << Hamiltonian[j * N_Atom_total * 2 + i].real() << "+(" << Hamiltonian[j * N_Atom_total * 2 + i].imag() << "i)" << "  ";
		}
		//H << endl;
		Eigenvalue_out_temp[j] = Eigenvalue[j];
	}
	/*delete[] Hamiltonian;
	delete[] HAMILTONIAN_MKL;
	delete[] Eigenvalue;
	delete[] eigenvectors;*/
}
