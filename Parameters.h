#ifndef _PARAMETERS_H_
#define	_PARAMETERS_H_

#include <complex>
#include <string>

const std::complex<double> I_unit(0.0, 1.0); //虚数单位i
const double Tor_bond = 0.2; //键长范围容许差异
const double Tor_onsite = 0.05; //视作同位
const double nTor = 1e-7; //电子密度精度
const int N_loop_times = 2000; //迭代容许范围
const double nTor_terminate = 1e-5; //跳出循环所要求的最低精度要求
const int N_loop_terminate = 5000; //强制跳出循环的迭代次数

//material specified
const double E0 = 0; //费米能
//on-site energy of Fe and O, e1= [8t^4 / (40K / 11605(K/eV))]^(1/3) - U 根据SAGA OF YIG和PRL关于铁磁系数的推导（只考虑最近邻Jad=40K）
//const double e1 = 1.833, e2 = 0;
const double e0 = -1.5, e1 = -0.44, e2 = -2.5;
const double t = 0.8;
const double U = 8.0; //Hubbard U
const double lamda_material = 1.13e-10; //Compton wavelength in materials, in SI unit, 注意系数部分用国际单位制，Ex可以都用Angstrom

//electric field
const double V0 = 3;
const double r0 = 100;
const double epsilon = 15.3;
const double d = 120;
#endif	/* _PARAMETERS_H_ */
