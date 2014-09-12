/*
	Written by CL Frank in 2014.6
*/
#ifndef _FINITE_ELEMENT_SOLVER_H_
#define _FINITE_ELEMENT_SOLVER_H_
#include <vector>
#include <set>
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Vec.h"
using std::vector;
using std::pair;
using std::set;

typedef Eigen::Matrix<double,2,6> Matrix2_6d;
typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Triplet<double> Tr;
typedef Eigen::SparseMatrix<double> SparseMatrixType;
#define BIGNUMERIC 1e16


class static_solver
{
public:
	//Young Modulus 
	double		            ex;
	//Possion Ratio
/*	double		           prxy;*/
	//constraint store the index of constraint nodes.
	vector<int>			    constraint; 
	//ext_f is the exterior force applied to nodes.
	vector<pair<int,vec> >  ext_f;
	//du is the displacement vector of nodes need to solve.
	Eigen::VectorXd		    du;
	//nodes of the truss structure
	vector<vec>			    nd;
	//radius array of strut unit
	vector<double>		    link_rad;
	//strut unit,contain the index of two end nodes.
	vector< pair<int,int> > un;
	//stress of strut unit
	vector<double>			stress;
	//将单元刚度矩阵从局部坐标系转换到全局坐标系的转置矩阵,动态指针数组
	//需要析构释放内存
	vector<Matrix2_6d*>	   Transto;
	//每一根杆的长度
	vector<double>		   link_length;

	~static_solver()
	{
		for (vector<Matrix2_6d*>::size_type i=0; i<Transto.size(); i++)
		{
			delete Transto[i];
			Transto[i] = NULL;
		}
	}
	void    readTrussData(const char* filename);
	void    save_results(const char* filename);
	void    Solve();
private:
	bool	SolveforDisplacement();
	void	SolveforStress();
	void	GetStiffnessMat(const vec& nd1,const vec& nd2,double rad,Matrix6d& K,int id);
	void	InsertTriplet(vector<Tr>& tpl,int r,int c,const Matrix6d& K);
	void	SetBigNum(SparseMatrixType& K,int i);
};

#endif