#include"Pub_head.h"
#include"PhysicsConstant.h"
#include"class_Read_Input.h"
#include"Pub_head.h"
#include"class_vectors.h"

//得到原子种类用于对Read_Input_CONTCAR类初始化
const int Read_Input_CONTCAR::N = 3;

void Read_Input_CONTCAR::Init0()
{
	int m_Max = 1;
	N_Species = 0;
	N_Atom_total = 0;
	N_Atom = new int[m_Max];
	Atom_type = new string[m_Max];
	Atom_Position = new double *[m_Max]; //Direct,首先new行
	Atom_Position_C = new double *[m_Max];//Cartestian
	for (int i = 0; i < m_Max; i++)
	{
		Atom_Position[i] = new double[N];//其次new列
		Atom_Position_C[i] = new double[N];
	}
	cerr << "Input error, program exit" << endl;
}

void Read_Input_CONTCAR::Init()
{
	vectors vec;
	if (get_N_Species() > 0)
	{
		N_Species = get_N_Species();
		N_Atom = new int[N_Species];
		Atom_type = new string[N_Species];
	}
	else
	{
		cout << "N_Species should be positive integers:" << "  ";
		Init0();
		Free();
		system("pause");
		exit(0);
	}

	ofstream info("info", ios::app);
	//开始读取POSCAR文件，这里不同类型顺序读取
	ifstream fin("CONTCAR");  //fin为读入流关键字，不需要声明。CONTCAR代表弛豫后结果

	if (!fin.is_open())
	{
		cout << "Initializing N_Atom_total" << "  ";
		Init0();
		cerr << "CONTCAR is not found in present directory, please check!" << endl;
		Free();
		system("pause");
		exit(0);
	}

	//POSCAR文件第一行为VESTA中的体系名称      
	fin >> system_name;
	info << system_name << endl;
	//来自POSCAR的第二行缩放系数
	fin >> factor;
	//info << factor << endl;  
	info << "Scaling factor=" << factor << endl; //output
	//3-5行为晶格基矢，...采用简化元胞减少计算
	info << "lattice:" << endl;
	for (int i = 0; i< N; i++)//基矢编号
	{
		for (int j = 0; j< N; j++)//xyz坐标
		{
			fin >> lattice[i].basis[j];
			info << lattice[i].basis[j] << "   ";
		}
		info << endl;
	}

	//calculate the lattice constant
	for (int i = 0; i < N; i++)
	{
		lattice[i].constant = vec.vec_scale(lattice[i].basis);
	}

	//第6行原子名称
	//info << "Atom types:" << endl;
	for (int i = 0; i < N_Species; i++)
	{
		fin >> Atom_type[i];
		//info<<Atom_type[i]<<"  ";
	}
	//info << endl;
	//第7行元胞中各原子数
	int N_Atom_total_temp = 0;
	for (int i = 0; i < N_Species; i++)
	{
		fin >> N_Atom[i]; //某种原子的原子数
		//info<<N_Atom[i]<<"  ";
		N_Atom_total_temp += N_Atom[i]; //得到元胞原子总数
	}
	N_Atom_total = N_Atom_total_temp;
	info << "Total atoms=" << N_Atom_total << endl;
	//已知原子数后对原子位置矩阵初始化
	Atom_Position = new double *[N_Atom_total]; //Direct,首先new行
	Atom_Position_C = new double *[N_Atom_total];//Cartestian
	for (int i = 0; i < N_Atom_total; i++)
	{
		Atom_Position[i] = new double[N];//其次new列
		Atom_Position_C[i] = new double[N];
	}

	//第8行坐标类型，分单位1的Direct，和实际晶格长度下的Cartesian
	fin >> is_direct;
	//info<<is_direct<<endl;

	//初始化Direct坐标
	for (int i = 0; i < N_Atom_total; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Atom_Position[i][j] = 0.0;
			Atom_Position_C[i][j] = 0.0;
		}
	}
	//读取元胞所有原子Direct坐标
	for (int i = 0; i<N_Atom_total; i++)
	{
		for (int j = 0; j < N; j++) //以三个基矢为1的坐标
		{
			fin >> Atom_Position[i][j];
		}
	}
	fin.close();
	//以上读取POSCAR完毕
	//Transfer Atom position from Direct to Cartesian 
	for (int i = 0; i<N_Atom_total; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)//转换到xyz坐标
			{
				Atom_Position_C[i][k] += Atom_Position[i][j] * lattice[j].basis[k] * factor;
			}

		}
	}
	get_lattice_Rec();
}

void Read_Input_CONTCAR::Free()
{
	delete[] N_Atom;
	delete[] Atom_type;
	/*
	for (int i = 0; i < N; i++)
	{
	delete[] Atom_Position[i];//首先delete列
	}
	for (int i = 0; i < N; i++)
	{
	delete[] Atom_Position_C[i];
	}*/
	delete[] Atom_Position;//其次delete行
	delete[] Atom_Position_C;
}

//默认构造函数
Read_Input_CONTCAR::Read_Input_CONTCAR()
{
	Init();
}
/*
//构造函数（如果考虑多个文件，应该以文件名为变量写到构造函数中，初始化采用Init0）
Read_Input_CONTCAR::Read_Input_CONTCAR(string get_file_number)
{
if (get_file_number == NULL)
Init0();
else
{}
}
*/
//常数初始化在Init0中的话必须同时终止，否则会重复初始化
//拷贝构造函数 
Read_Input_CONTCAR::Read_Input_CONTCAR(const Read_Input_CONTCAR& struc)
{
	N_Species = struc.N_Species; /*复制常规成员*/
	//m_Max = struc.m_Max;
	factor = struc.factor;
	is_direct = struc.is_direct;
	system_name = struc.system_name;
	const int N_Atom_total = struc.N_Atom_total;

	for (int i = 0; i < N; i++) /*类定义的成员*/
	{
		lattice[i] = struc.lattice[i];
		lattice_Rec[i] = struc.lattice_Rec[i];
	}

	N_Atom = new int[N_Species];  /*复制指针指向的内容*/
	Atom_type = new string[N_Species];
	Atom_Position = new double *[N_Atom_total]; //Direct
	Atom_Position_C = new double *[N_Atom_total];//Cartestian
	for (int i = 0; i < N_Atom_total; i++)
	{
		Atom_Position[i] = new double[N];
		Atom_Position_C[i] = new double[N];
	}
	memcpy(N_Atom, struc.N_Atom, N_Species*sizeof(int));
	memcpy(Atom_type, struc.Atom_type, N_Species*sizeof(string));
	memcpy(Atom_Position, struc.Atom_Position, N_Species*N*sizeof(double));
	memcpy(Atom_Position_C, struc.Atom_Position_C, N_Species*N*sizeof(double));
}

//析构函数
Read_Input_CONTCAR::~Read_Input_CONTCAR()
{
	Free();
}

//public function
int Read_Input_CONTCAR::get_N_Species()
{
	int N_Spec = 1;
	string im;
	int start;
	ifstream con("CONTCAR");
	if (!con.is_open())
	{
		cout << "Initializing N_Species" << "  ";
		Init0();
		cerr << "CONTCAR is not found in present directory, please check!" << endl;
		Free();
		system("pause");
		exit(0);
	}
	for (int i = 1; i <= 6; i++)
	{
		getline(con, im); //此处循环为了使con历经到6
	}
	con.close();

	int len = im.size() - 1;
	if (len>1)
	{
		for (int i = 0; i<len; ++i)
		{
			if (im[i] != ' ')
			{
				start = i;
				break;
			}
		}

		for (int i = start; i<len; ++i)
		{
			if (im[i] == ' '&&im[i + 1] != ' ')
			{
				N_Spec++; //给出原子种类
				//cout << N_Spec << endl;
			}
		}
	}
	return N_Spec;
}

void Read_Input_CONTCAR::get_lattice_Rec()
{
	ofstream info("info", ios::app);
	//Calculate reciprocal lattice
	//倒格矢分母计算
	double omega;
	vectors vec;
	double *vec_out = new double[N];
	vec.vec_cross(vec_out, lattice[1].basis, lattice[2].basis);
	omega = vec.vec_dot(lattice[0].basis, vec_out);
	//倒格矢计算
	info << "Reciprocal lattice:" << endl;
	for (int i = 0; i < N; i++)
	{
		int j = i + 1;
		if (j > 2) { j = j - N; }
		int k = j + 1;
		if (k > 2) { k = k - N; }
		vec.vec_cross(vec_out, lattice[j].basis, lattice[k].basis);
		for (int m = 0; m < N; m++)
		{
			lattice_Rec[i].basis[m] = 2.0 * PI * vec_out[m] / omega;
			info << lattice_Rec[i].basis[m] << "  ";
		}
		info << endl;
	}
	delete[] vec_out;
}


