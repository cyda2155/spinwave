//所有extern变量在主程序赋值,对于const类型初始化语句为const int xx = xx;
//只能在文件函数体外使用显式表达式初始化extern const int xx = const.，在其他文件调用时引用该头
/*
extern const int N_Species;
extern const int N_Atom_total;
extern const int Kpoints;
extern const int Nk_total;
extern const int Nk_max;

extern const int atoms_in_supercell;
extern const double lamda_so;
*/
#ifndef _EXTERNAL_H_
#define _EXTERNAL_H_

extern string Fe_a;
extern string Fe_d;
extern string O;

#endif

