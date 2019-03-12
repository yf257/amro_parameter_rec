#pragma once
class calculateAMRO
{
private:
	int veloXsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);//derivative of e[kx,ky,kz]
	int veloYsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
	int veloZsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

	int veloX(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);//derivative of eigen[{e[kx,ky,kz],delta},{delta,e[kx+Pi/a,ky+Pi/a,kz]}]
	int veloY(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);
	int veloZ(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);

	int func(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);

	int veloXsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);//derivative of e[kx,ky,kz]
	int veloYsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
	int veloZsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

	int veloXD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);//derivative of eigen[{e[kx,ky,kz],delta},{delta,e[kx+Pi/a,ky+Pi/a,kz]}]
	int veloYD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);
	int veloZD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);

	int funcD(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);


	int fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
	int fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
	int fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out);

	//int taufun(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out);
	int taufun(double *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *one);
	int taufundos(Ipp64f *params, Ipp64f minDos, Ipp64f maxDos, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones);
	Ipp64f *thetas;
	Ipp64f *phis;
	double *params;
	FindFermi *Fermi;
	double *_data;
	Ipp64f *tempfunc;

	Ipp64f *circum;
	int _dataLeng;
	Ipp64f *condout;
	Ipp64f tau;
	Ipp64f final;//time final?
	long steps;//number of time steps?
	Ipp64f h; //delta time
	Ipp64f field45; // 45 tesla in appropriate units
	//FindFermi Fermi("name", params);
	int nPoints;
	Ipp64f *output; //stores evolution of orbit around Fermi surface
	Ipp64f *times; //time steps
	Ipp64f *starts;
	//std::clock_t startT;
	Ipp64f duration;
	
	Ipp64f *field;
	Ipp64f *vzStorage;
	Ipp64f *vz0Storage;
	Ipp64f *DOS;
	Ipp64f *ones;//for inverting
	//ippsSet_64f(1, ones, nPoints);
	Ipp64f *taus;//phi dependent taus

	Ipp64f *vx;
	Ipp64f *vy;
	Ipp64f *vz;

	Ipp64f *argx;
	Ipp64f *argy;
	Ipp64f *argz;

	Ipp64f *tempx;
	Ipp64f *tempy;
	Ipp64f *tempz;

	Ipp64f *k1x;
	Ipp64f *k1y;
	Ipp64f *k1z;
	Ipp64f *k2x;
	Ipp64f *k2y;
	Ipp64f *k2z;
	Ipp64f *k3x;
	Ipp64f *k3y;
	Ipp64f *k3z;
	Ipp64f *k4x;
	Ipp64f *k4y;
	Ipp64f *k4z;//what are these? why do we need k1-4?
	Ipp64f *tempdif;
	Ipp64f *exptau;
	Ipp64f minDos;
	Ipp64f maxDos;
	Ipp64f total;

	Ipp64f *temp ;
	Ipp64f *upper_bound ;
	
	Ipp64f *lower_bound ;
	
	Ipp64f *over_shift;
	
	Ipp64f *less_shift;
	
	Ipp64f *zeros;//for inverting
	

	Ipp64f sum_resdual;
	Ipp64f Resdual;
	int updatetheta(Ipp64f *theta, int Length);
	int updatedata(Ipp64f *data, int Length);
	int updatephi(Ipp64f *phi, int Length);
public:

	int updatepara(double *param, int Length);
	calculateAMRO( double * data, double * param, Ipp64f * theta, int cdev, int gridN, int dataLeng,Ipp64f * phi);
	int printPar();
	Ipp64f returnvalue(double *param);
	int writefile(double * param);
	virtual ~calculateAMRO();
};

