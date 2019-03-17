#pragma once
class FindFermi
{
private:
	
	int func(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);
	int funcf2(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *functemp);
	int funcf1(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *functemp);

	int funcd(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);
	int funcfd2(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *functemp);
	int funcfd1(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *functemp);

	//int nPoints;
	Ipp64f *subMaxR;
	Ipp64f *temp1;
	Ipp64f *temp2;
	Ipp64f *argCos;
	Ipp64f *argSin;
	Ipp64f *argkz;
	Ipp64f *zeros;
	Ipp64f *kz;
	Ipp64f *funcval;
	Ipp64f *tfunc1;
	Ipp64f *tfunc2;
	Ipp64f *tfunc;
	Ipp64f *tfunc3;
	int *lengArr;
	Ipp64f *kx;
	Ipp64f *ky;
	Ipp64f *tempx1;
	Ipp64f *tempy1;
	Ipp64f *tempx2;
	Ipp64f *tempy2;
	Ipp64f *tempx3;
	Ipp64f *tempy3;
	Ipp64f *tempx4;
	Ipp64f *tempy4;
	Ipp64f *tempz;
	Ipp64f *circumf;
	//Ipp64f *startpoint ;
	Ipp64f tol;

	int fineN;
	int NumLeng;
	int cdevs;
	int nfinepoint;
	//Ipp64f finDis;
	int finP;
	Ipp64f params[9];
public:
	



	int nPoints;
	FindFermi(Ipp64f * param, int cdev,int Ngrid);
	int UpdatePar(double *param);
	int PrintPar();
	int ReturnStart(Ipp64f *startpoint);
	int ReturnNumPoint();
	int ReturnCircum(Ipp64f *circum);
	//int interfunc(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int *Laylength, int Laynum, Ipp64f distance, Ipp64f *temp, Ipp64f *temp2, Ipp64f *outx, Ipp64f *outy, Ipp64f*outz);
	int interfunc(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int *Laylength, int Laynum, int finePoint, Ipp64f *temp, Ipp64f *temp2, Ipp64f *outx, Ipp64f *outy, Ipp64f*outz, Ipp64f*outcircum);

	//int interfunc(Ipp64f *theta, Ipp64f *kx, Ipp64f *ky, int Laylength, int clength, int flength, Ipp64f *temp, Ipp64f *temp2, Ipp64f *out);
	virtual ~FindFermi();
};

