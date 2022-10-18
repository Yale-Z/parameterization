#ifndef TUTTE_H_
#define TUTTE_H_
#include<Eigen/LU>
#include<vector>
#include<Eigen/Sparse>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

};

typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> Mesh;

inline int findBoundary(Mesh& mesh_, std::vector<std::pair<int,Eigen::Vector2d>>& bc_) {
	std::pair<int, Eigen::Vector2d> p;
	for (auto he_it = mesh_.halfedges_begin(); he_it != mesh_.halfedges_end(); ++he_it) {
		if (mesh_.is_boundary(*he_it)) {
			p.first = (*he_it).from().idx();
			bc_.emplace_back(p);
			 auto vh = (*he_it).to();
			while (vh.idx()!=bc_[0].first){
				p.first = vh.idx();
				bc_.emplace_back(p);
				for (auto voh_it = mesh_.voh_begin(vh); voh_it != mesh_.voh_end(vh); ++voh_it){
					if (mesh_.is_boundary(*voh_it)) {
						vh = (*voh_it).to();
						break;
					}
				}
			}
			break;
		}
	}

	std::cout  <<"\n";
	return 0;
}
inline int fixBoundaryCircle(std::vector<std::pair<int, Eigen::Vector2d>>& bc_, double area_) {
	int n = bc_.size();
	const double pi = 3.141592653589793;
	const double r = sqrt(area_ / pi);
	const double unit = 2 * pi / n;
	for (int i = 0; i < n; ++i) {
		const double x = cos(unit * i);
		const double y = sin(unit * i);
		bc_[i].second(0) = r * x ;
		bc_[i].second(1) = r * y ;
		//bc_[i].second(0) = abs(x) > 1e-8 ? r * x : 0.;
		//bc_[i].second(1) = abs(y) > 1e-8 ? r * y : 0.;
	}
	return 0;
}
inline double caculateArea(Mesh& mesh_) {
	double sum = 0.;
	for (auto f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) {
		auto fv_it = mesh_.fv_iter(*f_it);
		const Mesh::Point& a = mesh_.point(fv_it); ++fv_it;
		const Mesh::Point& b = mesh_.point(fv_it); ++fv_it;
		const Mesh::Point& c = mesh_.point(fv_it);
		sum += ((b - a) % (c - a)).norm() * 0.5f;
	}
	return sum;
}
inline int getFaces(Mesh& mesh_, Eigen::Matrix3Xi& F_) {
	const size_t n = mesh_.n_faces();
	F_.resize(3, n);
	for (auto f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) {
		int idx = (*f_it).idx();
		auto fv_it = mesh_.fv_iter(*f_it);
		int i = 0;
		F_.col(idx)(i) = (*fv_it).idx(); ++i;++fv_it;
		F_.col(idx)(i) = (*fv_it).idx(); ++i;++fv_it;
		F_.col(idx)(i) = (*fv_it).idx();
	}
	return 0;
}
inline int getVertices(Mesh& mesh_, Eigen::Matrix3Xd& V_) {
	const size_t n = mesh_.n_vertices();
	V_.resize(3, n);
	for (auto v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
		int idx = (*v_it).idx();
		V_.col(idx) = Eigen::Map<Eigen::Vector3d>(mesh_.point(*v_it).data());
	}
	return 0;
}
inline int cacCot(Eigen::Matrix3Xd& V_, Eigen::Matrix3Xi& F_, Eigen::Matrix3Xd& cot_) {
	cot_.resize(3, F_.cols());
	for (size_t j = 0; j < F_.cols(); j++) {
		const Eigen::Vector3i& v_ids = F_.col(j);
		for (int v = 0; v < 3; v++) {
			int v_id[3];
			v_id[0] = v_ids(v);
			v_id[1] = v_ids((v + 1) % 3);
			v_id[2] = v_ids((v + 2) % 3);
			const Eigen::Vector3d& v1 = V_.col(v_id[0]);//point A of triangle ABC
			const Eigen::Vector3d& v2 = V_.col(v_id[1]);//point B of triangle ABC
			const Eigen::Vector3d& v3 = V_.col(v_id[2]);//point C of triangle ABC
			const double angle = std::acos(std::max(-1.0,
				std::min(1.0, (double)((v2 - v1).normalized().dot((v3 - v1).normalized())))));
			float cot = 1.0 / std::tan(angle);
			cot_(v, j) = cot;//cot of A
		}
	}
	return 0;
}
inline int getWeight(Mesh& mesh_, Eigen::SparseMatrix<double>& laplace_) {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	Eigen::Matrix3Xd cot;
	getFaces(mesh_, F);
	getVertices(mesh_, V);
	cacCot(V, F, cot);
	const size_t n = mesh_.n_vertices();
	laplace_.resize(n, n);
	std::vector<Eigen::Triplet<double>> triplet_laplace;//triplet laplace
	triplet_laplace.reserve(F.cols() * 9);
	for (auto f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) {
		int fid = (*f_it).idx();
		const Eigen::Vector3i& vids = F.col(fid);
		const Eigen::Vector3d& cots = cot.col(fid);
		for (auto fv_it = mesh_.fv_begin(f_it); fv_it != mesh_.fv_end(f_it); ++fv_it) {
			if (mesh_.is_boundary(fv_it))	continue;
			auto v_it = fv_it;
			int a = (*v_it).idx(); v_it++;
			int b = (*v_it).idx(); v_it++;
			int c = (*v_it).idx(); int ai, bi, ci;
			for (int i = 0; i < 3; ++i) {
				if (vids[i] == a) ai = i;
				else if (vids[i] == b) bi = i;
				else if (vids[i] == c) ci = i;
			}
			triplet_laplace.emplace_back(Eigen::Triplet<double>(a, a, cots(bi) + cots(ci)));
			triplet_laplace.emplace_back(Eigen::Triplet<double>(a, b, -cots(ci)));
			triplet_laplace.emplace_back(Eigen::Triplet<double>(a, c, -cots(bi)));
		}
	}
	laplace_.setFromTriplets(triplet_laplace.begin(), triplet_laplace.end());
	return 0;
}
inline int tutte(Mesh& mesh_) {
	std::vector<std::pair<int, Eigen::Vector2d>> bc;
	findBoundary(mesh_,bc);
	Eigen::Matrix3Xi F;
	const double area = caculateArea(mesh_);
	fixBoundaryCircle(bc, area);
	Eigen::SparseMatrix<double> laplace;
	getWeight(mesh_, laplace);
	int n = mesh_.n_vertices();
	Eigen::Matrix2Xd b(2, n);
	b.setZero();
	for (auto& elem : bc) {
		int idx = elem.first;
		b.col(idx) = elem.second;
		laplace.coeffRef(idx, idx) = 1;
	}
	
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
	solver.analyzePattern(laplace);
	solver.factorize(laplace);
	Eigen::MatrixX2d xt = solver.solve(b.transpose());
	Eigen::Matrix2Xd x = xt.transpose();

	//std::cout << laplace << std::endl;
	for (auto v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) {
		int idx = (*v_it).idx();
		mesh_.set_point((*v_it), Mesh::Point(x.col(idx)(0), x.col(idx)(1),0.));
	}
	return 0;
};
#endif // !TUTTE_H_

