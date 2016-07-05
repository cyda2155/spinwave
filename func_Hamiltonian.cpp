#include"headfile_othercpp.h"

//注意！！！每次调用前Hamiltonian必须清零
void func_Hamiltonian(const int N_Atom_total, const int atoms_in_supercell, const double lamda_so, complex<double> *Hamiltonian, const double *d_bond, double *kkk, const double *Electric_field, double *n, const Atom *unit_cell, const Atom *supercell)
{
	//给出参数
	//the hopping constant，tij	
	double Matrix = 0.0;//

	//SOC, c(i)->c(i)exp(I_unit*lamda_so*(x cross E)dot sigma), 两者相减，分成3个坐标在t项出现，作为小量按e指数展开
	complex<double> t_lamda_so[3]; //能量单位为eV

	double distance;//Rij

	vectors vec;

	bool flag00, flag01;
	bool flag10, flag11, flag20, flag21;
	bool flag1, flag2, flag;

	for (int atom_in_cell = 0; atom_in_cell < N_Atom_total; atom_in_cell++)//元胞原子循环
	{
		for (int atom_adjacent = 0; atom_adjacent < atoms_in_supercell; atom_adjacent++)//超胞相邻原子循环（实际上这里包含了和atom_in_cell同位原子的情形）
		{
			double *vec_out = new double[3];
			double *vec_out1 = new double[3];
			vec.vec_minus(vec_out, supercell[atom_adjacent].position, unit_cell[atom_in_cell].position);
			complex<double> eikd = exp(I_unit * vec.vec_dot(kkk, vec_out));

			distance = vec.vec_distance(supercell[atom_adjacent].position, unit_cell[atom_in_cell].position);  //给出Rij

			//We can initial the diagonal matrix element
			if (abs(distance) < Tor_onsite)//同位原子，直接对应写到Hamiltonian
			{
				if (unit_cell[atom_in_cell].species == Fe_a || unit_cell[atom_in_cell].species == Fe_d)
				{
					//Hubbard U只加到Fe上
					//上对角块
					Hamiltonian[unit_cell[atom_in_cell].atom_label * N_Atom_total * 2 + unit_cell[atom_in_cell].atom_label] += unit_cell[atom_in_cell].On_Site + U * (n[atom_in_cell + N_Atom_total] - 0.5);
					//下对角块
					Hamiltonian[(unit_cell[atom_in_cell].atom_label + N_Atom_total) * N_Atom_total * 2 + (unit_cell[atom_in_cell].atom_label + N_Atom_total)] += unit_cell[atom_in_cell].On_Site + U * (n[atom_in_cell] - 0.5);
				}
				if (unit_cell[atom_in_cell].species == O)
				{
					//上对角块
					Hamiltonian[unit_cell[atom_in_cell].atom_label * N_Atom_total * 2 + unit_cell[atom_in_cell].atom_label] += unit_cell[atom_in_cell].On_Site;
					//下对角块
					Hamiltonian[(unit_cell[atom_in_cell].atom_label + N_Atom_total) * N_Atom_total * 2 + (unit_cell[atom_in_cell].atom_label + N_Atom_total)] += unit_cell[atom_in_cell].On_Site;
				}
			}

			else
			{
				//Fe-O hopping
				//根据SAGA给定的ad,dd位距离
				flag00 = distance > abs(d_bond[0] - Tor_bond) && distance < abs(d_bond[2] - Tor_bond);
				flag01 = distance > abs(d_bond[1] - Tor_bond) && distance < abs(d_bond[2] - Tor_bond);
				flag10 = unit_cell[atom_in_cell].species == Fe_a && supercell[atom_adjacent].species == O;
				flag11 = unit_cell[atom_in_cell].species == Fe_d && supercell[atom_adjacent].species == O;
				flag20 = unit_cell[atom_in_cell].species == O && supercell[atom_adjacent].species == Fe_a;
				flag21 = unit_cell[atom_in_cell].species == O && supercell[atom_adjacent].species == Fe_d;

				flag1 = (flag10 && flag00) || (flag20 && flag00);
				flag2 = (flag11 && flag01) || (flag21 && flag01);
				flag = flag1 || flag2;

				if (flag == true)
				{
					Matrix = -t;
					//计入SOC，直接对应写到Hamiltonian。注意单带模型中SOC只对跃迁产生影响。而磁性只和跃迁相关
					vec.vec_cross(vec_out1, vec_out, Electric_field);

					//x分量，产生在2*2HAMILTON非对角块，在虚部
					//y分量产生在非对角块实部
					//z分量产生在对角元虚部
					for (int m = 0; m < 3; m++)
					{
						t_lamda_so[m] = I_unit * lamda_so * vec_out1[m];
					}

					//上对角块
					Hamiltonian[unit_cell[atom_in_cell].atom_label * N_Atom_total * 2 + supercell[atom_adjacent].atom_label]
						+= (Matrix * exp(t_lamda_so[2])) * eikd; //i轨道-j轨道，最后按TB展开
					//非对角块为0：不考虑SOC,跃迁仅产生在sigma相同的位点间
					//上非对角块
					Hamiltonian[atom_in_cell * N_Atom_total * 2 + (supercell[atom_adjacent].atom_label + N_Atom_total)] += (Matrix * exp(t_lamda_so[0] - I_unit * t_lamda_so[1])) * eikd;
					//下非对焦块
					Hamiltonian[(atom_in_cell + N_Atom_total) * N_Atom_total * 2 + supercell[atom_adjacent].atom_label] += (Matrix * exp(t_lamda_so[0] + I_unit * t_lamda_so[1])) * eikd;
					//下对角块
					Hamiltonian[(unit_cell[atom_in_cell].atom_label + N_Atom_total) * N_Atom_total * 2 + (supercell[atom_adjacent].atom_label + N_Atom_total)]
						+= (Matrix * exp(-t_lamda_so[2])) * eikd;
				}
			}
			delete[] vec_out;
			delete[] vec_out1;
		}
	}
}
