#include "FiniteElemSolver.h"

void static_solver::readTrussData(const char* filename)
{
	FILE *file;
	if((file=fopen(filename,"r"))==NULL)
	{
		printf("file not exists,please check the path of the file");
		return;
	}
	char str[80];
	while (fscanf(file,"%s",str) == 1)
	{
		if (strncmp(str,"elasticModulus",14) == 0)
		{
			float tmpex;
			fscanf(file,"%f",&tmpex);
			ex = tmpex;
		}
		else if (strncmp(str,"nodes",5) == 0)
		{
			int node_count;
			fscanf(file,"%d",&node_count);
			for (int i=0; i<node_count; i++)
			{
				float x,y,z;
				fscanf(file,"%f %f %f",&x,&y,&z);
				nd.push_back(vec(x,y,z));
			}
		}
		else if (strncmp(str,"units",5) == 0)
		{
			int unit_count;
			fscanf(file,"%d",&unit_count);
			un.resize(unit_count);
			link_rad.resize(unit_count);
			for (int i=0; i<unit_count; i++)
			{
				int point1,point2; float rad;
				fscanf(file,"%d %d %f",&point1,&point2,&rad);
				un[i].first = point1;
				un[i].second = point2;
				link_rad[i] = rad;
			}
		}
		else if (strncmp(str,"force",5) == 0)
		{
			int force_count;
			fscanf(file,"%d",&force_count);
			for (int i=0; i<force_count; i++)
			{
				int f_index; float x,y,z;
				fscanf(file,"%d %f %f %f",&f_index,&x,&y,&z);
				ext_f.push_back(pair<int,vec>(f_index,vec(x,y,z)));
			}
		}
		else if (strncmp(str,"constraint",10) == 0)
		{
			int cons_count;
			fscanf(file,"%d",&cons_count);
			for (int i=0; i<cons_count; i++)
			{
				int cons_index;
				fscanf(file,"%d",&cons_index);
				constraint.push_back(cons_index);
			}
		}
		else
		{
			printf("wrong file\n");
		}
	}
	fclose(file);
}

void static_solver::Solve()
{
	if (SolveforDisplacement())
	{
		printf("success\n");
		SolveforStress();
	}
	else
	{
		printf("there are some errors in the truss structure");
		return;
	}
}
bool static_solver::SolveforDisplacement()
{
	int unit_count = un.size();
	int nd_count = nd.size();
	int dimf_count = 3 * nd_count;
	vector<Tr> tripletList;
	tripletList.reserve(36*unit_count);
	Transto.resize(unit_count);
	link_length.resize(unit_count);

	for (int i=0; i<unit_count; i++)
	{
		//p is first node of the strut , q is the second 
		int p = un[i].first, q = un[i].second;
		Matrix6d uk;
		GetStiffnessMat(nd[p],nd[q],link_rad[i],uk,i);
		InsertTriplet(tripletList,p,q,uk);
	}

	SparseMatrixType Kw(3*nd_count,3*nd_count);
	Kw.setFromTriplets(tripletList.begin(),tripletList.end());

	Eigen::VectorXd Fr(dimf_count);
	Fr.setZero();

	//set external force 
	for (size_t i=0; i<ext_f.size(); i++)
	{
		for (int j=0; j<3; j++)
			Fr(3*ext_f[i].first+j) = ext_f[i].second[j];
	}

	//set constraint 
	for (size_t i=0; i< constraint.size(); i++)
	{
		for (int j=0; j<3; j++)
			Fr(3*constraint[i]+j) = 0;

		SetBigNum(Kw,constraint[i]);
	}

	
	Eigen::SimplicialCholesky<SparseMatrixType> K_Solver;
	K_Solver.compute(Kw);
  	if (K_Solver.info() != Eigen::Success)
	{
		return false;
	}         
	
	du = K_Solver.solve(Fr);
	return true;
}

void static_solver::SetBigNum(SparseMatrixType& K,int i)
{
	i *= 3;
	for (int j=0; j<3; j++)
	{
		for (int k=0; k<3; k++)
		{
			K.coeffRef(i+j,i+k) = BIGNUMERIC * K.coeff(i+j,i+k);
		}
	}
}

//     r      c
//   * * *  * * *
//r  * * *  * * *
//   * * *  * * * 
//
//   * * *  * * * 
//c  * * *  * * *
//   * * *  * * * 
void static_solver::InsertTriplet(vector<Tr>& tpl,int r,int c,const Matrix6d& K)
{
	r *= 3; c *= 3;

	for (int i=0,k=0; i<3; i++,k++)
		for (int j=0,l=0; j<3; j++,l++)
			tpl.push_back(Tr(r+k,r+l,K(i,j)));

	for (int i=0,k=0; i<3; i++,k++)
		for (int j=3,l=0; j<6; j++,l++)
			tpl.push_back(Tr(r+k,c+l,K(i,j)));

	for (int i=3,k=0; i<6; i++,k++)
		for (int j=0,l=0; j<3; j++,l++)
			tpl.push_back(Tr(c+k,r+l,K(i,j)));

	for (int i=3,k=0; i<6; i++,k++)
		for (int j=3,l=0; j<6; j++,l++)
			tpl.push_back(Tr(c+k,c+l,K(i,j)));
}

inline void static_solver::GetStiffnessMat(const vec& nd1,const vec& nd2,
										  double rad, Matrix6d& K,int id)
{
	assert(rad > 0);
	const double pi = 3.14159265359;
	double ld = len(nd2 - nd1);
	assert(ld > 0);
	double area = pi * rad * rad;
	double coefficient = ex * area / ld;

	//cos(x,x~) = (x2 - x1) / l
	double md = (nd2[0] - nd1[0]) / ld;
	//cos(x,y~) = (y2 - y1) / l
	double nd = (nd2[1] - nd1[1]) / ld;
	//cos(x,z~) = (z2 - z1) / l
	double rd = (nd2[2] - nd1[2]) / ld;

	Eigen::Matrix<double,2,6>* Te = new Eigen::Matrix<double,2,6>;
	*Te << md, nd, rd, 0, 0, 0,
		   0, 0, 0, md, nd, rd;

	Eigen::Matrix2d Ke;
	Ke << 1, -1, -1, 1;

	K = coefficient * ((*Te).transpose() * Ke * (*Te));

	Transto[id] = Te;
	link_length[id] = ld;
}

void static_solver::SolveforStress()
{
	assert(!un.empty());
	if (un.empty())
		return;

	stress.clear();
	int un_count = un.size();
	stress.resize(un_count);
	Eigen::Matrix<double,1,2> unit;
	unit << -1.0,1.0;
	double coeffs = 0.0;
	int nd1,nd2;

	for (int i=0; i<un_count; i++)
	{
		coeffs = ex / link_length[i];
		nd1 = 3 * un[i].first; 
		nd2 = 3 * un[i].second;
		Eigen::Matrix<double,6,1> dq;
		dq << du[nd1],du[nd1+1],du[nd1+2],
			  du[nd2],du[nd2+1],du[nd2+2];

		stress[i] = coeffs * unit * (*Transto[i]) * dq;
	}
}

void static_solver::save_results(const char* filename)
{
	FILE *fp;
	char triname[256];

	strcpy(triname,filename);

	if ((fp = fopen(triname,"w")) == NULL)
	{
		printf("error");
		return;
	}

	fprintf(fp,"results data\n");
	fprintf(fp,"Nodes Displacement\n");
	int nd_count = du.size() / 3;
	for (int i=0; i<nd_count; i++)
	{
		int n_index = i*3;
		fprintf(fp,"Node%d:%.3f,%.3f,%.3f\n",i,du[n_index],du[n_index+1],du[n_index+2]);
	}

	fprintf(fp,"Element Unit Stress\n");
	for (size_t i=0; i<stress.size(); i++)
	{
		fprintf(fp,"Unit%d:%.3f\n",i,static_cast<float>(stress[i]));
	}
	fclose(fp);
}
