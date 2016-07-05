#include"Pub_head.h"
//class DArray;
class lattices
{
public:
	double basis[3];//xyz方向
	double constant;
};

class Read_Input_CONTCAR
{
private:
	string is_direct; //type of coordinates
	string system_name; //POSCAR给定的体系名称
	double **Atom_Position;
private:
	void Init();//初始化并读取POSCAR信息
	void Free();//删除动态内存
	void Init0();
	//string keyword;
public:
	static const int N; //维数
	int N_Atom_total; //原胞内总原子数
	int N_Species;//元素种类	
	int *N_Atom;//各种元素的原子数	
	string *Atom_type; //原子类型名称
	lattices lattice[3]; //3表示有三个元胞基矢
	lattices lattice_Rec[3];//存储倒格矢
	double factor;//缩放系数
	double **Atom_Position_C;
public:
	Read_Input_CONTCAR();//默认构造函数
	//Read_Input_CONTCAR(int get_N_Species, double dValue);//构造函数
	Read_Input_CONTCAR(const Read_Input_CONTCAR& struc); // 拷贝构造函数
	~Read_Input_CONTCAR();    // 析构函数 

	int  get_N_Species();//逐个读取用于初始化
	void get_lattice_Rec();//求倒格矢
	//Read_Input_CONTCAR CONTCAR();//读取POSCAR文件返回本类
};

class Read_Input_KPOINT
{
private:

private:
	void Init();//初始化并读取POSCAR信息
	void Free();//删除动态内存
	void Init0();
public:
	static const int N; //维数
	int Kpoints; //高对称K点数
	int Nk_total; //k网格总点数
	int Nk_max; //K附近最密的k点数，用于输出方便
	int *Nk_Kpoint;//K附近k点数（代表网格大小）
	double **K_path;//各高对称K点位置
public:
	Read_Input_KPOINT();//默认构造函数
	Read_Input_KPOINT(const Read_Input_KPOINT& struc); // 拷贝构造函数
	~Read_Input_KPOINT();    // 析构函数
};
