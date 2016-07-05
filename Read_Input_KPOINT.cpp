#include"Pub_head.h"
#include"class_Read_Input.h"
#include"class_vectors.h"

//得到原子种类用于对Read_Input_KPOINT类初始化
const int Read_Input_KPOINT::N = 3;

void Read_Input_KPOINT::Init0()
{
	int m_Max = 1;
	Kpoints = 0;
	Nk_total = 0;
	Nk_max = 0;
	Nk_Kpoint = new int[m_Max];
	K_path = new double *[m_Max]; //首先new行
	for (int i = 0; i < m_Max; i++)
	{
		K_path[i] = new double[N];//其次new列
	}
	cerr << "Input error, program exit" << endl;
}

void Read_Input_KPOINT::Init()
{
	//read kpath
	ifstream kin("KPOINT"); //记录K点的文件
	ofstream info("info", ios::app);
	if (!kin.is_open())
	{
		cout << "Initializing Kpoints:" << "  ";
		Init0();
		cerr << "KPOINT is not found in present directory, please check!" << endl;
		Free();
		system("pause");
		exit(0);
	}
	//K点数

	kin >> Kpoints;
	info << Kpoints << endl;

	//K网格
	int Nk_max_temp = 0;
	if (Kpoints > 0)
	{
		Nk_Kpoint = new int[Kpoints - 1];
	}
	else
	{
		cout << "Kpoints should be positive integers:" << "  ";
		Init0();
		Free();
		system("pause");
		exit(0);
	}
	for (int i = 0; i < (Kpoints - 1); i++)
	{
		kin >> Nk_Kpoint[i];
		if (Nk_Kpoint[i] > Nk_max_temp)
		{
			Nk_max_temp = Nk_Kpoint[i];
		}
		info << "division of k net=" << Nk_Kpoint[i] << "   ";
	}
	Nk_max_temp += 1; //这里考虑到实际循环变量范围
	Nk_max = Nk_max_temp;
	info << endl;

	//各K点位置
	K_path = new double *[Kpoints]; //首先new行
	for (int i = 0; i < Kpoints; i++)
	{
		K_path[i] = new double[N];//其次new列
	}
	info << "K position" << endl;
	for (int i = 0; i < Kpoints; i++)
	{
		for (int j = 0; j < N; j++)
		{
			kin >> K_path[i][j];
			info << K_path[i][j] << "   ";
		}
		info << endl;
	}
	kin.close();
	Nk_total = (Kpoints - 1) * Nk_max;
}

void Read_Input_KPOINT::Free()
{
	delete[] Nk_Kpoint;
	for (int i = 0; i < N; i++)
	{
		delete[] K_path[i];//首先delete列
	}
	delete[] K_path;//其次delete行
}

//默认构造函数
Read_Input_KPOINT::Read_Input_KPOINT()
{
	Init();
}

//拷贝构造函数 
Read_Input_KPOINT::Read_Input_KPOINT(const Read_Input_KPOINT& struc)
{

	Kpoints = struc.Kpoints; /*复制常规成员*/
	Nk_total = struc.Nk_total;
	Nk_max = struc.Nk_max;

	Nk_Kpoint = new int[Kpoints - 1];
	K_path = new double *[Kpoints]; //首先new行
	for (int i = 0; i < Kpoints; i++)
	{
		K_path[i] = new double[N];//其次new列
	}

	memcpy(Nk_Kpoint, struc.Nk_Kpoint, (Kpoints - 1)*sizeof(int));
	memcpy(K_path, struc.K_path, Kpoints*N*sizeof(double));
}

//析构函数
Read_Input_KPOINT::~Read_Input_KPOINT()
{
	Free();
}

