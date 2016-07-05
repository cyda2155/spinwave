#include"headfile_main.h"
#include <omp.h>

//取出原子根据CONTCAR文件
string Fe_a = "Fe_a";
string Fe_d = "Fe_d";
string O = "O";
int main()
{
	//IO stream setup
	cout << setiosflags(ios::fixed) << setiosflags(ios::left) << setprecision(10); //设置输出format

	ofstream info("info");
	info << setiosflags(ios::fixed) << setiosflags(ios::left) << setprecision(10);

	//利用读取文件类接收信息
	Read_Input_CONTCAR const contcar;
	Read_Input_KPOINT const kpoint;

	//const变量作用域和位置有关，在开头显式定义作用域为本文件，这里在main内
	const int N_Atom_total = contcar.N_Atom_total;
	const int N_Species = contcar.N_Species;
	const int Kpoints = kpoint.Kpoints;
	const int Nk_total = kpoint.Nk_total;
	const int Nk_max = kpoint.Nk_max;

	cout << "Input succeed!" << endl;

	//variables
	vectors vec;

	//建立元胞
	//用Atom类给出元胞中每个原子的原子种类，位置，轨道信息，unit_cell成员数为原子个数
	//double bondlength[2];
	Atom *unit_cell_temp = new Atom[N_Atom_total]; //Atom类在Atom.h中。赋给元胞每个原子

	int atom_count = 0;
	for (int i = 0; i < N_Species; i++) //按原子种类顺序取出
	{
		for (int j = 0; j < contcar.N_Atom[i]; j++)
		{
			unit_cell_temp[atom_count].species = contcar.Atom_type[i]; //成员带有原子种类信息
			unit_cell_temp[atom_count].atom_label = atom_count; //label即矩阵元编号
			for (int k = 0; k < contcar.N; k++)
			{
				unit_cell_temp[atom_count].position[k] = contcar.Atom_Position_C[atom_count][k]; //元胞原子位置传递
			}
			atom_count++;
		}
	}

	//给出unit_cell每个原子on site energy
	for (int i = 0; i< N_Atom_total; i++)
	{
		//Fe on site energy
		if (unit_cell_temp[i].species == Fe_a)
		{
			unit_cell_temp[i].On_Site = e0;
		}
		if (unit_cell_temp[i].species == Fe_d)
		{
			unit_cell_temp[i].On_Site = e1;
		}
		//O on site energy
		if (unit_cell_temp[i].species == O)
		{
			unit_cell_temp[i].On_Site = e2;
		}
	}

	const Atom *unit_cell = unit_cell_temp;

	//用于存储Fe-O\Fe-Fe两种键长，前面是a位，后面是d位
	double d_bond_temp[4];
	func_d_bond(N_Atom_total, d_bond_temp, unit_cell);
	const double *d_bond = d_bond_temp;

	//建立超胞
	int adjacent_cell_max = 3;  //考虑1个元胞距离内的原子，构成3*3*3元胞的超胞
	int cell_abs_max = adjacent_cell_max / 2;
	const int atoms_in_supercell = pow(adjacent_cell_max, 3) * N_Atom_total; //超胞总原子数
	info << "Number of atoms in the supercell=" << atoms_in_supercell << endl;

	Atom *supercell_temp = new Atom[atoms_in_supercell]; //supercell中的原子仍用atom类进行继承

	atom_count = 0;
	//i,j,zz为三个单位基矢整数倍
	for (int i = -cell_abs_max; i <= cell_abs_max; i++)
	{
		for (int j = -cell_abs_max; j <= cell_abs_max; j++)
		{
			for (int zz = -cell_abs_max; zz <= cell_abs_max; zz++)
			{
				for (int m = 0; m < N_Atom_total; m++)
				{
					for (int k = 0; k < 3; k++)//xyz三个方向
					{
						supercell_temp[atom_count].position[k] = unit_cell[m].position[k] + contcar.lattice[0].basis[k] * i + contcar.lattice[1].basis[k] * j + contcar.lattice[2].basis[k] * zz;
					}
					supercell_temp[atom_count].species = unit_cell[m].species; //这两句考虑了平移不变性
					supercell_temp[atom_count].atom_label = unit_cell[m].atom_label;
					atom_count++;
				}
			}
		}
	}

	const Atom *supercell = supercell_temp;

	//计算哈密顿矩阵元
	const double lamda_so = m_electron * e_charge / pow(h_bar, 2) * pow(lamda_material, 2); //系数为国际单位

	//能带总数
	const int N_orbits = 1; //轨道数
	const int N_spin_components = 2; //自旋取向
	const int N_k_band = N_Atom_total * N_orbits * N_spin_components; //总的能带数
	const int N_spin_average = N_Atom_total * N_orbits; //合计自旋向上向下的带的自旋,假设SOC不明显

	//为了计算占据对应的n用于迭代，首先计算Fe原子数，基态只有Fe是单占据的，参与填充
	int Fe_a_label;
	int Fe_d_label;
	int O_label;
	for (int i = 0; i < N_Species; i++)
	{
		if (contcar.Atom_type[i] == Fe_a)
			Fe_a_label = i;
		if (contcar.Atom_type[i] == Fe_d)
			Fe_d_label = i;
		if (contcar.Atom_type[i] == O)
			O_label = i;
	}
	const int Fe_a_N_Atoms = contcar.N_Atom[Fe_a_label];
	const int Fe_d_N_Atoms = contcar.N_Atom[Fe_d_label];
	const int Fe_N_Atoms = Fe_a_N_Atoms + Fe_d_N_Atoms;
	const int O_N_Atoms = contcar.N_Atom[O_label];
	
	//电场
	const int E_dim = 3;//电场维数这里手动设定
	double *Electric_field = new double[E_dim];
	for (int i = 0; i < E_dim; i++)
	{
		Electric_field[i] = 0;
	}
	//取点设定
	const int circle = 10; //圆环个数
	const int point = pow((circle + 1), 2); //计算的点数
	const double dR = 12.4; //原胞长度

	//自旋维度
	const int spin_dim = 3;
	//输出原胞自旋
	ofstream cell_s("CELLSPIN.dat");
	cell_s << setiosflags(ios::fixed) << setiosflags(ios::left) << setprecision(10);
	cell_s << "x  " << "y  " << "S_x  " << "S_y  " << "S_z  " << endl;

	//对圆环圈数循环
	double r,x,y,theta,Et,En;
	for (int l_circle = 0; l_circle < circle + 1; l_circle++)
	{
		r = dR * l_circle;
		//计算Et
		Et = 2 * V0 * r / ((1 + epsilon) / r0 + (1 - epsilon) / (2 * d - r0)) / pow((pow(d, 2) + pow(r, 2)), 1.5);
		En = 2 * V0 * d / ((1 + epsilon) / r0 + (1 - epsilon) / (2 * d - r0)) / pow((pow(d, 2) + pow(r, 2)), 1.5);
		for (int i_point = 0; i_point < 2 * l_circle + 1; i_point++)
		{
			if (l_circle == 0)
			{
				theta = 0;
			}
			else
			{
				theta = PI / (2 * l_circle) * i_point;
			}
			x = r * cos(theta);
			y = r * sin(theta);
			Electric_field[0] = Et * cos(theta);
			Electric_field[1] = Et * sin(theta);
			Electric_field[2] = En;
			//为便于并行，将电场在内循环内定义为常量，const不影响变量的作用域
			const double *E_field = Electric_field;

			//计算各S.必须对ci整体计算，这里用于存sum(Ne)(c(dagger)_nks*c_nks')项
			complex<double> **Spin_x = new complex<double> *[Nk_total];
			complex<double> **Spin_y = new complex<double> *[Nk_total];
			complex<double> **Spin_z = new complex<double> *[Nk_total];
			for (int k = 0; k < Nk_total; k++) //行数
			{
				//new列数
				Spin_x[k] = new complex<double>[N_Atom_total];
				Spin_y[k] = new complex<double>[N_Atom_total];
				Spin_z[k] = new complex<double>[N_Atom_total];
			}
			for (int k = 0; k < Nk_total; k++)
			{
				for (int i = 0; i < N_Atom_total; i++)
				{
					Spin_x[k][i] = 0;
					Spin_y[k][i] = 0;
					Spin_z[k][i] = 0;
				}
			}

			//openmp并行针对的循环
			//采用openmp并行
#pragma omp parallel
			for (int ii_index = 1; ii_index < Kpoints; ii_index++)
			{
				//循环Kpoints-1次，分别给出K0->K1和K1->K2....的kk网格
#pragma omp for
				for (int jj_index = 0; jj_index < (kpoint.Nk_Kpoint[ii_index - 1] + 1); jj_index++)
				{
					const int ii = ii_index;
					const int jj = jj_index;
					const int k_index = (ii - 1)* Nk_max + jj;
					//kk网络的三个坐标，注意循环是重复赋值的
					double *kk = new double[kpoint.N]; //= { 0.0, 0.0, 0.0 };
					for (int m = 0; m < kpoint.N; m++)
					{
						kk[m] = (kpoint.K_path[ii][m] - kpoint.K_path[ii - 1][m]) / kpoint.Nk_Kpoint[ii - 1] * jj + kpoint.K_path[ii - 1][m];
					}
					//将kk坐标转为真实的倒易空间坐标kkk，lattice_Rec为基矢，同前面实空间原子位置坐标的处理
					double *kkk = new double[kpoint.N];
					for (int m = 0; m < kpoint.N; m++)
					{
						kkk[m] = kk[0] * contcar.lattice_Rec[0].basis[m] + kk[1] * contcar.lattice_Rec[1].basis[m] + kk[2] * contcar.lattice_Rec[2].basis[m];
					}
					//为中间变量申请内存（注意需要定义为循环内的局部指针变量）
					//前面指标代表一个本征矢内的分量（各原子轨道系数），后面指标代表本征矢编号(根据heev_cpp测试结果)
					complex<double> **eigenvectors_out = new complex<double> *[N_k_band]; //存储迭代结果
					for (int i = 0; i < N_k_band; i++) //行数
					{
						eigenvectors_out[i] = new complex<double>[N_k_band];//new列数
					}

					for (int j = 0; j < N_k_band; j++)
					{
						for (int i = 0; i < N_k_band; i++)
						{
							eigenvectors_out[i][j] = 0;
						}
					}
					//输入元胞，k点，电场，迭代求能量本征方程
					double *Eigenvalue_out = new double[N_k_band];//与本征态同时计算的本征值
					func_iteration(Fe_N_Atoms, O_N_Atoms, N_Atom_total, atoms_in_supercell, d_bond, lamda_so, kkk, E_field, unit_cell, supercell, Eigenvalue_out, eigenvectors_out);

					//输出计算的网格位置
					cout << "kpoint:" << "(" << ii << ", " << jj << ")" << "  ";
					cout << endl;

					for (int i = 0; i < N_Atom_total; i++)
					{
						for (int j = 0; j < Fe_N_Atoms / 2 + O_N_Atoms; j++) //最低能态求和得到ck，各个格点i独立
						{
							Spin_x[k_index][i] += std::conj(eigenvectors_out[i][j]) * eigenvectors_out[i + N_Atom_total][j] + std::conj(eigenvectors_out[i + N_Atom_total][j]) * eigenvectors_out[i][j];
							Spin_y[k_index][i] += -I_unit * (std::conj(eigenvectors_out[i][j]) * eigenvectors_out[i + N_Atom_total][j] - std::conj(eigenvectors_out[i + N_Atom_total][j]) * eigenvectors_out[i][j]);
							Spin_z[k_index][i] += std::conj(eigenvectors_out[i][j]) * eigenvectors_out[i][j] - std::conj(eigenvectors_out[i + N_Atom_total][j]) * eigenvectors_out[i + N_Atom_total][j];
						}
					}
										
					//注意循环内new的新指针必须delete否则内存泄露
					delete[] kk;
					delete[] kkk;

					delete[] Eigenvalue_out;
					delete[] eigenvectors_out;
				}
			}
			//omp并行结束

			//原胞自旋
			complex<double> *cellspin = new complex<double>[spin_dim];
			for (int i = 0; i < spin_dim; i++)
			{
				cellspin[i] = 0;
			}
			//格点自旋值
			complex<double> *citespin_x = new complex<double>[N_Atom_total];
			complex<double> *citespin_y = new complex<double>[N_Atom_total];
			complex<double> *citespin_z = new complex<double>[N_Atom_total];
			for (int i = 0; i < N_Atom_total; i++)
			{
				citespin_x[i] = 0;
				citespin_y[i] = 0;
				citespin_z[i] = 0;
			}

			//对每个格点计算非均匀k网格下的自旋平均（给出格点自旋值）
			for (int i = 0; i < N_Atom_total; i++)
			{
				double dk_total = 0;
				for (int ii = 0; ii < Kpoints - 1; ii++)
				{
					complex<double> citespin_x_temp = 0;
					complex<double> citespin_y_temp = 0;
					complex<double> citespin_z_temp = 0;
					double dk = vec.vec_distance(kpoint.K_path[ii + 1], kpoint.K_path[ii]);
					//先对内层k，即两个特殊K点之间的网格求和
					for (int jj = 0; jj < (kpoint.Nk_Kpoint[ii] + 1); jj++)
					{
						int k_index = ii* Nk_max + jj;
						citespin_x_temp += Spin_x[k_index][i];
						citespin_y_temp += Spin_y[k_index][i];
						citespin_z_temp += Spin_z[k_index][i];
					}
					//除Nk_ii:区间内归一化，乘以权重对ii求和
					citespin_x[i] += citespin_x_temp / (double)(kpoint.Nk_Kpoint[ii] + 1) * dk;
					citespin_y[i] += citespin_y_temp / (double)(kpoint.Nk_Kpoint[ii] + 1) * dk;
					citespin_z[i] += citespin_z_temp / (double)(kpoint.Nk_Kpoint[ii] + 1) * dk;
					//总权重
					dk_total += dk;
				}
				citespin_x[i] = citespin_x[i] / dk_total;
				citespin_y[i] = citespin_y[i] / dk_total;
				citespin_z[i] = citespin_z[i] / dk_total;
			}

			//计算总的块自旋
			for (int i = 0; i < N_Atom_total; i++)
			{
				cellspin[0] += citespin_x[i];
				cellspin[1] += citespin_y[i];
				cellspin[2] += citespin_z[i];
			}
			//输出
			cell_s << x << "  " << y << "  ";
			for (int i = 0; i < spin_dim; i++)
			{
				cell_s << 1.0/2.0 * cellspin[i].real() << "  ";
			}
			cell_s << endl;
			//delete myloc
			//指向const的指针（const *）不用delete，保护被指地址的内容不被更改
			//const指针只能指向初始化的地址

			delete[] Spin_x;
			delete[] Spin_y;
			delete[] Spin_z;
			delete[] citespin_x;
			delete[] citespin_y;
			delete[] citespin_z;
			delete[] cellspin;
		}
		cell_s << endl;
	}

	delete[] unit_cell_temp;
	delete[] supercell_temp;
	delete[] Electric_field;
	//delete[] d_bond_temp;
	cout << "Calculation succeed!" << endl;
	return 0;
}
