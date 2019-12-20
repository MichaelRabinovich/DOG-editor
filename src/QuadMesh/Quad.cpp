#include "Quad.h"

#include <igl/boundary_loop.h>
#include <igl/is_border_vertex.h>
#include <igl/vertex_triangle_adjacency.h>

#include <set>

using namespace std;

void quad_topology(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr, 
                    QuadTopology& quadTop) {
  Eigen::MatrixXi F = Fsqr_to_F(Fsqr); // yuck!
  get_stars_and_bnd_vertices_nbds(V,F,Fsqr,quadTop.stars,quadTop.vi_to_star,quadTop.bnd4,quadTop.bnd3,quadTop.bnd2);
  igl::boundary_loop(F, quadTop.bnd_loop);
  quadTop.is_bnd_v = igl::is_border_vertex(V,F);
  quad_edges(Fsqr,quadTop.E);
  vertex_quad_adjacency(Fsqr, quadTop.VF);
  quad_adjacenecy_list(Fsqr,quadTop.A);
  quadTop.v_n = V.rows();
}


 std::vector<bool> is_border_vertex(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr) {
  Eigen::MatrixXi F = Fsqr_to_F(Fsqr); // yuck!
  return igl::is_border_vertex(V,F);

}
int num_vertex_sharing(int f1, int f2, const Eigen::MatrixXi& F) {
  std::vector<int> F1_V; std::vector<int> F2_V;
  for (int i = 0; i < 4; i++) {
    F1_V.push_back(F(f1,i));
    F2_V.push_back(F(f2,i));
  }
  std::sort(F1_V.begin(), F1_V.end());
  std::sort(F2_V.begin(), F2_V.end());
  std::set<int> intersect;
  std::set_intersection(F1_V.begin(),F1_V.end(),F2_V.begin(),F2_V.end(),
    std::inserter(intersect,intersect.begin()));
  return intersect.size();
}

int num_face_nbd(std::vector<int>& VFi, int f, const Eigen::MatrixXi& F) {
  int cnt = 0;
  for (int i = 0; i < VFi.size(); i++) {
    if (f != VFi[i]) {
      if (num_vertex_sharing(f,VFi[i],F) == 2) {
        cnt++;
      }
    }
  }
  return cnt;
}
void vertex_quad_adjacency(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& VF) {
  int n = F.maxCoeff()+1;
  VF.clear();
  VF.resize(n);

  for(int fi=0; fi<F.rows(); ++fi) {
    for(int i = 0; i < F.cols(); ++i) {
      VF[F(fi,i)].push_back(fi);
    }
  }
  // Order them in cyclic order (only needed if there are 3 or 4 quads around a vertex)
  for (int i = 0; i< n; i++) {
    if (VF[i].size() < 3) continue;
    //if (VF[i].size() == 3) {
    {
      int f1 = VF[i][0], f2 = VF[i][1], f3 = VF[i][2];
      if (num_face_nbd(VF[i],f1,F) >= 2) {
        VF[i][0] = f2;
        VF[i][1] = f1;
      }
      
      if (VF[i].size() == 4) {

        f1 = VF[i][0], f2 = VF[i][1], f3 = VF[i][2];
        if (num_face_nbd(VF[i],f1,F) >= 2) {
          std::swap(VF[i][0],VF[i][2]);
        }
      }

      // now we know that the first one is good (this is always true for 4/2 faces around a vertex)
    }
    // check if the next face is in good order
    int f1 = VF[i][0], f2 = VF[i][1], f3 = VF[i][2];
    if (num_vertex_sharing(f1,f2,F) < 2) {
      VF[i][1] = f3;
      VF[i][2] = f2;
      f1 = VF[i][0], f2 = VF[i][1], f3 = VF[i][2];
      if (num_vertex_sharing(f1,f2,F) < 2) {
        VF[i][1] = VF[i][3];
        VF[i][3] = f1;
      }
    }
    // It can still be bad if there are 4 faces
    if (VF[i].size() == 4) {
      f1 = VF[i][0], f2 = VF[i][1], f3 = VF[i][2]; int f4 = VF[i][3];
      if (num_vertex_sharing(f2,f3,F) < 2) {
        VF[i][2] = f4;
        VF[i][3] = f3;
      }
    }
  }
}

int num_face_sharing(int v1, int v2, std::vector<std::vector<int>>& VF) {
  std::sort(VF[v1].begin(),VF[v1].end());
  std::sort(VF[v2].begin(),VF[v2].end());

  std::set<int> intersect;
  std::set_intersection(VF[v1].begin(),VF[v1].end(),VF[v2].begin(),VF[v2].end(),
    std::inserter(intersect,intersect.begin()));
  return intersect.size();
}

int find_bnd_vertex(std::vector<int>& nbd, std::vector<std::vector<int>>& VF) {
  for (int i = 0; i < nbd.size(); i++) {
    int num_sharing_v = 0;
    for (int j = 0; j < nbd.size(); j++) {
      if ((i!=j) && (num_face_sharing(nbd[i],nbd[j],VF) >= 1 )) {
        num_sharing_v++;
      }
    }
    if (num_sharing_v < 2) return i;
  }
  return -1;
}

void quad_adjacenecy_list(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& A) {
  A.clear(); 
  A.resize(F.maxCoeff()+1);
  
  // Loop over faces
  for(int i = 0;i<F.rows();i++) {
    for(int j = 0;j<F.cols();j++) {
      // Get indices of edge: s --> d
      int s = F(i,j);
      for (int k = 0; k < F.cols(); k++) {
        if ( (k-j)%2 == 1) {
          int d = F(i,k);
          A.at(s).push_back(d);
          A.at(d).push_back(s);
        }
      }
    }
  }
  // Remove duplicates
  for(int i=0; i<(int)A.size();++i) {
    std::sort(A[i].begin(), A[i].end());
    A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
  }
  // (Locally) orient neighbours
  std::vector<std::vector<int>> VF,VFi; igl::vertex_triangle_adjacency(F.maxCoeff()+1, F, VF, VFi);
  for(int i=0; i<(int)A.size();++i) {
    // If the first vertex and the second vertex are "opposite vertices" (i.e they share no face, change the order)
    if (A[i].size() >= 3) {
      //if (i == 1) {for (int idx = 0; idx < A[i].size(); idx++) cout <<"A[i][idx] = " << A[i][idx] << endl;}
      int bnd_idx = find_bnd_vertex(A[i],VF);
      // if there is no bnd vertex, any vertex could be the first
      if (bnd_idx != -1) {
        std::swap(A[i][bnd_idx],A[i][0]);
      }
      // swap second and third if needed
      if (num_face_sharing(A[i][0],A[i][1],VF) < 1) {
        std::swap(A[i][1],A[i][2]);
        if (num_face_sharing(A[i][0],A[i][1],VF) < 1) {
          std::swap(A[i][1],A[i][3]); // can only happen with 5 nbd after cutting
        }
      }
      // If there are 4 vertices we know that we now have 0,1 but can have 3,2 instead of 2,3
      // If there are 3 vertices we should be ok
      if ( (A[i].size () >= 4) && (num_face_sharing(A[i][1],A[i][2],VF) < 1))  {
        std::swap(A[i][2],A[i][3]);  
      }
    }
  }
}

std::vector<int> get_edge_faces(const QuadTopology& quadTop, Edge& e) {
  Eigen::VectorXi e1(2); e1 << e.v1,e.v2;
  return get_edge_faces(quadTop, e1);
}

std::vector<int> get_edge_faces(const QuadTopology& quadTop, const Eigen::RowVectorXi& e) {
  std::vector<int> v1_nbd = quadTop.VF[e(0)], v2_nbd = quadTop.VF[e(1)];
  std::sort(v1_nbd.begin(),v1_nbd.end()); std::sort(v2_nbd.begin(), v2_nbd.end()); // set_intersection needs this to be sorted!
  std::vector<int> edge_faces;
  std::set_intersection(v1_nbd.begin(), v1_nbd.end(), v2_nbd.begin(), v2_nbd.end(), 
                        std::inserter(edge_faces, edge_faces.begin()));

  return edge_faces;
}

void quad_edges(const Eigen::MatrixXi& F, Eigen::MatrixXi& E) {
  // build adjacency matrix
  //Eigen::MatrixXi A;
  std::vector<std::vector<int> > A;
  quad_adjacenecy_list(F,A);
  int nnz = 0; 
  for (int i = 0; i < A.size(); i++) nnz += A[i].size();
  assert(nnz%2 == 0);
  // Resize to fit edges
  E.resize(nnz/2,2); int e_i = 0;
  for (int i =0; i < A.size(); i++) {
    for (int j = 0; j < A[i].size(); j++) {
      if (i < A[i][j]) {
        E(e_i,0) = i;
        E(e_i,1) = A[i][j];
        e_i++;
      }
    }
  }
}

void get_edges_from_face(const Eigen::MatrixXi& F, int sqr_f, std::vector<Edge>& quad_edges) {
  std::set<Edge> tri1_edges,tri2_edges;

  tri1_edges.insert(Edge(F(2*sqr_f,0),F(2*sqr_f,1)));
  tri1_edges.insert(Edge(F(2*sqr_f,0),F(2*sqr_f,2)));
  tri1_edges.insert(Edge(F(2*sqr_f,1),F(2*sqr_f,2)));

  tri2_edges.insert(Edge(F(2*sqr_f+1,0),F(2*sqr_f+1,1)));
  tri2_edges.insert(Edge(F(2*sqr_f+1,0),F(2*sqr_f+1,2)));
  tri2_edges.insert(Edge(F(2*sqr_f+1,1),F(2*sqr_f+1,2)));

  //for (auto edge: tri1_edges ) {cout << "Tri1 edge = " << edge.v1 << "," << edge.v2 << endl;}
  //for (auto edge: tri2_edges ) {cout << "Tri2 edge = " << edge.v1 << "," << edge.v2 << endl;}

  std::set<Edge> all_edges;
  std::set_union(tri1_edges.begin(), tri1_edges.end(), tri2_edges.begin(), tri2_edges.end(), 
                        std::inserter(all_edges, all_edges.begin()));

  //for (auto edge: all_edges ) {cout << "all_edges edge = " << edge.v1 << "," << edge.v2 << endl;}

  std::set<Edge> non_quad_edge;
  std::set_intersection(tri1_edges.begin(), tri1_edges.end(), tri2_edges.begin(), tri2_edges.end(), 
                        std::inserter(non_quad_edge, non_quad_edge.begin()));

  //for (auto edge: non_quad_edge ) {cout << "non_quad_edge edge = " << edge.v1 << "," << edge.v2 << endl;}

  // quad edges are the union of edges "minus" the non-quad edge, which is the shared edge between the triangle faces
  std::set_difference(all_edges.begin(), all_edges.end(), non_quad_edge.begin(), non_quad_edge.end(), 
                        std::inserter(quad_edges, quad_edges.begin()));

  if (quad_edges.size() !=4) {
    exit(1); //I am unforigiving
  }
}


void get_stars_and_bnd_vertices_nbds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& Fsqr, Eigen::VectorXi& stars, std::vector<int>& vi_to_star,
    Eigen::VectorXi& bnd4, Eigen::VectorXi& bnd3, Eigen::VectorXi& bnd2) {

  // get quad adjacency list
  std::vector<std::vector<int> > A;
  quad_adjacenecy_list(Fsqr,A);
  std::vector<bool> is_border = igl::is_border_vertex(V,F);
  vi_to_star.resize(V.rows()); for (int i = 0; i < vi_to_star.size(); i++) vi_to_star[i] = -1; // -1 if there's no star..
  int stars_cnt = 0;

  std::vector<int> stars_v,bnd4_v,bnd3_v,bnd2_v;
  for (int i = 0; i < A.size(); i++) {
    switch (A[i].size()) {
      case 4: {
        if (!is_border[i]) {
          stars_v.push_back(i);
          //cout << "star at index = " << i << endl;
          for (int j = 0; j < A[i].size(); j++) {stars_v.push_back(A[i][j]);}
          vi_to_star[i] = stars_cnt*5; // every star contain 5 indices (the center and nbd)
          stars_cnt++;
        } else {
          bnd4_v.push_back(i);  
          //cout << "bnd 4 at index = " << i << endl;
          for (int j = 0; j < A[i].size(); j++) {bnd4_v.push_back(A[i][j]);}
        }
        break;
      } case 3: {
        bnd3_v.push_back(i);
        //if (i == 1) {for (int idx = 0; idx < A[i].size(); idx++) cout <<"A[i][idx] = " << A[i][idx] << endl;}
        //cout << "bnd 3 at index = " << i << endl;
        for (int j = 0; j < A[i].size(); j++) {bnd3_v.push_back(A[i][j]);}
        //if (i==1) {for (int idx = 0; idx < A[i].size()+1; idx++) cout <<"bnd3_v[idx] = " << bnd3_v[idx] << endl;}
        //if (i == 1) {for (int idx = 0; idx < A[i].size()+1; idx++) cout <<"bnd3_v[idx] = " << A[i][idx] << endl;}
        break;
      } case 2: {
        bnd2_v.push_back(i);
        for (int j = 0; j < A[i].size(); j++) {bnd2_v.push_back(A[i][j]);}
        //cout << "bnd 2 at index = " <<  i << endl;
        break;
      case 5: {
        // add 1 star for first vertices
        bnd4_v.push_back(i);
        for (int j = 0; j < A[i].size()-1; j++) {bnd4_v.push_back(A[i][j]);}
          
        bnd4_v.push_back(i);
        // add 2nd star for other vertices
        for (int j = 1; j < A[i].size(); j++) {bnd4_v.push_back(A[i][j]);}
        
        //cout << "Warning, got 5 nbd, did you cut anything?" << endl;
        break;
      }
      } default: {
        //cerr << "Error: #Nbd should be between 2-5, but it is " << A[i].size() << " for i = " << i << endl;
        //for (int j = 0; j < A[i].size(); j++) cout << "At idx " << j << " it is "<< A[i][j] << endl;
      }

    }
  }
  // flatten vectors to eigen vectors
  ivec_to_eigen(stars_v,stars);
  ivec_to_eigen(bnd4_v,bnd4);
  ivec_to_eigen(bnd3_v,bnd3);
  ivec_to_eigen(bnd2_v,bnd2);
}

void ivec_to_eigen(const std::vector<int>& vec, Eigen::VectorXi& ei_v) {
  ei_v.resize(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    ei_v(i) = vec[i];
  }
}

void mat2_to_vec(const Eigen::MatrixXd& mat, Eigen::VectorXd& vec) {
  // x1,x2,...., y1,y2,...
  //Eigen::MatrixXd tmp = mat.transpose();
  //Eigen::MatrixXd tmp = mat;
  //tmp.resize(mat.cols()*mat.rows(),1);
  //vec = tmp;
  vec.resize(mat.cols()*mat.rows());
  for (int i = 0; i < mat.cols(); i++) {
    for (int j = 0; j < mat.rows(); j++) {
        vec(i*mat.rows()+j) = mat(j,i);
    }
  }
}

void vec_to_mat2(const Eigen::VectorXd& vec, Eigen::MatrixXd& mat) {
  int v_n = vec.rows()/3;
  mat.resize(vec.rows()/3, 3);
  for (int i = 0; i < mat.rows(); i++) {
    //mat.row(i) << vec(i*3), vec(i*3+1),vec(i*3+2);
    mat.row(i) << vec(i), vec(i + v_n),vec(i+2*v_n);
  }
}


Eigen::MatrixXi F_to_Fsqr(const Eigen::MatrixXi& F) {
  Eigen::MatrixXi Fsqr(F.rows()/2,4);
  for (int i = 0; i < F.rows(); i+=2) {
    Fsqr.row(i/2) << F(i+1,0),F(i+1,1),F(i+1,2),F(i,2);
    // which means   s_t,s_t1,s1_t1,s1_t (i.e. in index notation F,F_1,F_12 F_2)
  }
  return Fsqr;
}

Eigen::MatrixXi Fsqr_to_F(const Eigen::MatrixXi& Fsqr) {
  Eigen::MatrixXi Ft(2*Fsqr.rows(),3);
  for (int i = 0; i < Fsqr.rows(); i++) {
    Ft.row(2*i) << Fsqr(i,0),Fsqr(i,2),Fsqr(i,3);
    Ft.row(2*i+1) << Fsqr(i,0),Fsqr(i,1),Fsqr(i,2);
  }
  Eigen::MatrixXi F;Eigen::VectorXi dummy;
  F = Ft;//igl::bfs_orient(Ft,F,dummy);
  return F;
}

void get_planar_square_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int square_h, int square_w) {
  get_planar_square_mesh_V(V,square_h,square_w);
  get_planar_square_mesh_F(F,square_h,square_w);
}


void get_planar_square_mesh_V(Eigen::MatrixXd& V, int square_h, int square_w) {
  const int s = square_h; const int t = square_w;
  V.resize(s*t,3); V.row(0) << 0,0,0;

  Eigen::RowVectorXd xStep(3); xStep << 1,0,0;
  for (int j = 1; j < t; j++) V.row(j) = V.row(j-1) + xStep;
  Eigen::RowVectorXd yStep(3); yStep << 0,1,0;
  for (int i = 1; i < s; i++) {
    V.row(i*t) = V.row((i-1)*t) + yStep;
    for (int j = 1; j < t; j++) {
      V.row(i*t+j) = V.row(i*t+j-1) + xStep;
    }
  }
}

void get_planar_square_mesh_F(Eigen::MatrixXi& F, int square_h, int square_w) {
  const int s = square_h; const int t = square_w;
  F.resize(2*(s-1)*(t-1),3); int f_idx = 0;
  for (int i = 0; i < s-1; i++) {
    for (int j = 0; j < t-1; j++) {
      int s_t = i*t + (j);
      int s1_t = (i+1)*t + (j);
      int s_t1 = i*t + (j+1);
      int s1_t1 = (i+1)*t + (j+1);

      F.row(f_idx++) << s_t,s1_t1,s1_t;
      F.row(f_idx++) << s_t,s_t1,s1_t1;
    }
  }
}
