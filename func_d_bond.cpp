#include"headfile_othercpp.h"

void func_d_bond(const int N_Atom_total, double *d_bond, const Atom *unit_cell)
{
	double distance;//Rij
	vectors vec;
	string atom_label[4];

	for (int i = 0; i < 4; i++)
	{
		d_bond[i] = 13.0; //单位A
		atom_label[i] = "";
	}

	bool position10, position11, position1;
	bool position20, position21, position2;

	//首先根据元胞原子排布计算键长
	for (int atom_in_cell = 0; atom_in_cell < N_Atom_total; atom_in_cell++)//元胞原子循环
	{
		for (int atom_adjacent = 0; atom_adjacent < N_Atom_total; atom_adjacent++)
		{
			if (atom_adjacent != atom_in_cell)
			{
				//求a位Fe-O键长
				position11 = unit_cell[atom_in_cell].species == Fe_a && unit_cell[atom_adjacent].species == O;

				position1 = position11;
				if (position1 == true)
				{
					double *vec_out = new double[3];
					distance = vec.vec_distance(unit_cell[atom_adjacent].position, unit_cell[atom_in_cell].position);
					if (distance < d_bond[0])
					{
						d_bond[0] = distance;
					}
					atom_label[0] = unit_cell[atom_in_cell].species + "-" + unit_cell[atom_adjacent].species;
					delete[] vec_out;
				}
				//求a位Fe - d位Fe键长，参考值3.46
				position11 = unit_cell[atom_in_cell].species == Fe_a && unit_cell[atom_adjacent].species == Fe_d;

				position1 = position11;
				if (position1 == true)
				{
					double *vec_out = new double[3];
					distance = vec.vec_distance(unit_cell[atom_adjacent].position, unit_cell[atom_in_cell].position);
					if (distance < d_bond[2])
					{
						d_bond[2] = distance;
					}
					atom_label[2] = unit_cell[atom_in_cell].species + "-" + unit_cell[atom_adjacent].species;
					delete[] vec_out;
				}

				//求d位Fe-O键长
				position21 = unit_cell[atom_in_cell].species == Fe_d && unit_cell[atom_adjacent].species == O;

				position2 = position21;
				if (position2 == true)
				{
					double *vec_out = new double[3];
					distance = vec.vec_distance(unit_cell[atom_adjacent].position, unit_cell[atom_in_cell].position);
					if (distance < d_bond[1])
					{
						d_bond[1] = distance;
					}
					atom_label[1] = unit_cell[atom_in_cell].species + "-" + unit_cell[atom_adjacent].species;
					delete[] vec_out;
				}

				//求d位Fe-d位Fe键长，参考值3.79
				position21 = unit_cell[atom_in_cell].species == Fe_d && unit_cell[atom_adjacent].species == Fe_d;

				position2 = position21;
				if (position2 == true)
				{
					double *vec_out = new double[3];
					distance = vec.vec_distance(unit_cell[atom_adjacent].position, unit_cell[atom_in_cell].position);
					if (distance < d_bond[3])
					{
						d_bond[3] = distance;
					}
					atom_label[3] = unit_cell[atom_in_cell].species + "-" + unit_cell[atom_adjacent].species;
					delete[] vec_out;
				}
			}
		}
	}
	ofstream d("bond_distance");
	d << setiosflags(ios::fixed) << setiosflags(ios::left) << setprecision(10);
	d << "d_bond_Fe-O = " << d_bond[0] << "(" << atom_label[0] << ")" << ", " << d_bond[1] << "(" << atom_label[1] << ")" << endl;
	d << "d_bond_Fe-Fe = " << d_bond[2] << "(" << atom_label[2] << ")" << ", " << d_bond[3] << "(" << atom_label[3] << ")" << endl;
}
