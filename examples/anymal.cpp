#include <stdio.h> 
#include <ctime>
#include <algorithm>

#include "fix64_scalar.h" 
#include "tiny_double_utils.h"
#include "tiny_matrix3x3.h"
#include "tiny_matrix_x.h"
#include "tiny_mb_constraint_solver.h"
#include "tiny_multi_body.h"
#include "tiny_pose.h"
#include "tiny_quaternion.h"
#include "tiny_urdf_structures.h"
#include "tiny_raycast.h"  
#include "tiny_rigid_body.h"
#include "tiny_urdf_to_multi_body.h"
#include "tiny_vector3.h"
#include "tiny_world.h"
#include "examples/motion_import.h"
#include "examples/tiny_urdf_parser.h"
#include "tiny_metrics.h"
// #include "examples/pybullet_urdf_import.h"
// #include "examples/pybullet_visualizer_api.h"



// typedef PyBulletVisualizerAPI VisualizerAPI;
typedef double Scalar; 
typedef DoubleUtils Utils;
 
std::string to_string(const Scalar &x) { 
  return std::to_string(Utils::getDouble(x));
} 
    
template <typename TinyScalar, typename TinyConstants>
struct UrdfToMultiBody2 { 
  typedef ::TinyUrdfStructures<TinyScalar, TinyConstants> TinyUrdfStructures;

  void convert(TinyUrdfStructures *urdf_structures,
               TinyWorld<TinyScalar, TinyConstants> *world,
               TinyMultiBody<TinyScalar, TinyConstants> *mb) {
 
    char search_path[4096];
    sprintf(search_path, "./data/");
    world->meshcat_viz.m_path_prefix = search_path;
    std::string texture_path = "checker_purple.jpg";
    world->meshcat_viz.convert_visuals(*urdf_structures, texture_path);
 
                   
    TinyUrdfToMultiBody<TinyScalar, TinyConstants>::convert_to_multi_body(
        *urdf_structures, *world, *mb);
    mb->initialize();
  }
}; 

using std::vector;
template <typename TinyScalar, typename TinyConstants>
void do_mix(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &state,
  vector<double> &q, vector<double> &qd) {
  int curr_q = 0, curr_qd = 0;
  for (auto *mb : world.m_multi_bodies) {
    for (int i = 0; i < mb->dof(); ++i)
      state.push_back(TinyConstants::scalar_from_double(q[curr_q++]));
    for (int i = 0; i < mb->dof_qd(); ++i)
      state.push_back(TinyConstants::scalar_from_double(qd[curr_qd++]));
  }
} 
 
   
template <typename TinyScalar, typename TinyConstants>
vector<TinyScalar> get_joints(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &q, vector<TinyScalar> &tans) {
  int curr_q = 0;
  for (auto *mb : world.m_multi_bodies) {
    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->forward_kinematics(tmp_v);

    curr_q += mb->dof();
    auto root = mb->get_world_transform(-1).m_translation;
    tans.push_back(root.x());
    tans.push_back(root.y());
    tans.push_back(root.z());
    for (int i = 0; i < mb->m_links.size(); ++i) {
      auto tr = mb->get_world_transform(i);
      auto pos = tr.m_translation;
      tans.push_back(pos.x());
      tans.push_back(pos.y());
      tans.push_back(pos.z());
    }
  }
  return tans;
}

template <typename TinyScalar, typename TinyConstants>
void adj_get_joints(
  TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &Rtans,
  vector<TinyScalar> &q, vector<double> &dldq) {
  int curr_q = 0; 
  int i_tan = 0;
  for (auto *mb : world.m_multi_bodies) {
    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->adj_constraint_set_zero();
    mb->adj_f_kinematics(tmp_v);

    TinyVector3<TinyScalar, TinyConstants> Rroot; 
    Rroot.set_zero();
    Rroot.m_x = Rtans[i_tan++];
    Rroot.m_y = Rtans[i_tan++];
    Rroot.m_z = Rtans[i_tan++];
    mb->adj_get_world_transform(-1).m_translation += Rroot;

    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyVector3<TinyScalar, TinyConstants> Rpos;
      Rpos.set_zero();
      Rpos.m_x = Rtans[i_tan++];
      Rpos.m_y = Rtans[i_tan++];
      Rpos.m_z = Rtans[i_tan++];
      mb->adj_get_world_transform(i).m_translation += Rpos;
    }

    TinyVectorX<TinyScalar, TinyConstants> Rq(tmp_v.size()), Rqd(0), Rqdd(0);
    Rq.set_zero();
    mb->adj_fk(Rq, Rqd, Rqdd, tmp_v);

    for (int i = 0; i < Rq.m_size; i++)
      dldq.push_back(TinyConstants::getDouble(Rq[i]));
    mb->adj_constraint_set_zero();
 
    curr_q += mb->dof();
  }
}

//input: dq: [sum_dof]
//output: tans: [(num_links + num_mb)*3]
template <typename TinyScalar, typename TinyConstants>
vector<double> forward_get_joints(vector<double> &dq, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));
  get_joints(world, q, tans);
  vector<double> ans;
  for (auto ele : tans)
    ans.push_back(TinyConstants::getDouble(ele));
  return ans;
}
   
//input: dq: [sum_dof], dldans: [(num_links + num_mb)*3]
//output: dldq: [sum_dof]
template <typename TinyScalar, typename TinyConstants>
vector<double> backward_get_joints(vector<double> &dq, vector<double> &dldans, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));

  vector<double> dldq;
  vector<TinyScalar> Rtans;
  for (auto ele : dldans)
    Rtans.push_back(TinyConstants::scalar_from_double(ele));
  adj_get_joints(world, Rtans, q, dldq);
  return dldq;  
}   
 

template <typename TinyScalar, typename TinyConstants>
vector<TinyScalar> get_com(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &q, vector<TinyScalar> &tans) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  int curr_q = 0;
  for (auto *mb : world.m_multi_bodies) {
    TinyVector3 ans(TinyConstants::zero(),TinyConstants::zero(),TinyConstants::zero());
    TinyScalar totm = TinyConstants::zero();
    mb->forward_kinematics(vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof()));
    curr_q += mb->dof();
    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyVector3 tmp = mb->get_world_com(i);
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      ans = ans + tmp * mass;
      totm = totm + mass;
    }
    ans = ans * (1./ totm);
    tans.push_back(ans[0]);
    tans.push_back(ans[1]);
    tans.push_back(ans[2]);
  }
  return tans;
}
 
//input: dq: [sum_dof]
//output: tans: [3*num_mb]
template <typename TinyScalar, typename TinyConstants>
vector<double> forward_get_com(vector<double> &dq, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));
  get_com(world, q, tans);
  vector<double> ans;
  for (auto ele : tans)
    ans.push_back(TinyConstants::getDouble(ele));
  return ans;
}

 

template <typename TinyScalar, typename TinyConstants>
void adj_get_com(
  TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &Rtans,
  vector<TinyScalar> &q, vector<double> &dldq) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;

  for (int ii = 0; ii < Rtans.size(); ii++) {
  }
  int curr_q = 0;
  int curr_mb_i = 0;
  for (auto *mb : world.m_multi_bodies) { 
 
    TinyVector3 Rtans_vec(Rtans[3*curr_mb_i], Rtans[3*curr_mb_i+1], Rtans[3*curr_mb_i+2]);
    curr_mb_i++; 

    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->adj_constraint_set_zero();
    mb->adj_f_kinematics(tmp_v);
    TinyScalar totm = TinyConstants::zero();
    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      totm = totm + mass;
    }

    for (int i = mb->m_links.size()-1; i >= 0; --i) {
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      TinyVector3 Rv = Rtans_vec * mass * (1./totm);
      mb->adj_get_world_com(i, Rv); 
    }

    TinyVectorX<TinyScalar, TinyConstants> Rq(tmp_v.size()), Rqd(0), Rqdd(0);
    Rq.set_zero();
    mb->adj_fk(Rq, Rqd, Rqdd, tmp_v);
 
    for (int i = 0; i < Rq.m_size; i++)
      dldq.push_back(TinyConstants::getDouble(Rq[i]));
    mb->adj_constraint_set_zero();

    curr_q += mb->dof();
  }

} 

//input: dq: [sum_dof], dldans: [3*num_mb]
//output: dldq: [sum_dof]
template <typename TinyScalar, typename TinyConstants>
vector<double> backward_get_com(vector<double> &dq, vector<double> &dldans, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq) {
    q.push_back(TinyConstants::scalar_from_double(ele));
  }

  vector<double> dldq;
  vector<TinyScalar> Rtans;
  for (auto ele : dldans) {
    Rtans.push_back(TinyConstants::scalar_from_double(ele));
  }
  adj_get_com(world, Rtans, q, dldq);

  return dldq;
}

template <typename TinyScalar, typename TinyConstants>
vector<vector<double> > forward_step(vector<double> &q, vector<double> &qd, 
  vector<double> &tau, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> state, taut;
  do_mix(world, state, q, qd);
  for (auto ele : tau)
    taut.push_back(TinyConstants::scalar_from_double(ele));
  TinyVectorX<TinyScalar, TinyConstants> taux(taut);
  world.adj_step_torch(state, taux);  
  vector<double> ansq, ansqd;
  for (auto *mb : world.m_multi_bodies) {
    auto state = mb->adjm_this_state;
    int dof1 = mb->dof(), dof2 = mb->dof_qd(); 
    for (int i = 0; i < dof1; ++i)
      ansq.push_back(TinyConstants::getDouble(state[i]));
    for (int i = 0; i < dof2; ++i)
      ansqd.push_back(TinyConstants::getDouble(state[i+dof1])); 
  } 
  return {ansq, ansqd};
}      
  
template <typename TinyScalar, typename TinyConstants>
vector<vector<double> > backward_step(vector<double> &q, vector<double> &qd, vector<double> &tau, 
  vector<double> &dldq, vector<double> &dldqd,  
  TinyWorld<TinyScalar, TinyConstants> &world) { 
  static int n_step = 0;  

  vector<TinyScalar> Rt, taut, state;
  do_mix(world, Rt, dldq, dldqd);  
  do_mix(world, state, q, qd);
  for (auto ele : tau)
    taut.push_back(TinyConstants::scalar_from_double(ele)); 
  TinyVectorX<TinyScalar, TinyConstants> R(Rt), taux(taut);     
                  
  std::vector<TinyVectorX<TinyScalar, TinyConstants>> ans, our_ans;
  TinyVectorX<TinyScalar, TinyConstants> tR(Rt), ttaux(taut), tstate(state.size());   
            
  ans = world.adj_back_step_man(R, taux, state); 
 
  auto ans_state = ans.back();
  vector<double> ansq, ansqd, anstau;  
  int curr_state = 0, curr_body = 0;
  for (auto *mb : world.m_multi_bodies) {
    int dof1 = mb->dof(), dof2 = mb->dof_qd(), dof3 = mb->dof_u();;
    for (int i = 0; i < dof1; ++i)
      ansq.push_back(TinyConstants::getDouble(ans_state[curr_state++]));
    for (int i = 0; i < dof2; ++i)
      ansqd.push_back(TinyConstants::getDouble(ans_state[curr_state++]));
    for (int i = 0; i < dof3; ++i)
      anstau.push_back(TinyConstants::getDouble(ans[curr_body][i]));
    curr_body++;
  }
   
  // if (n_step==0)
  //   exit(0);
  n_step++;
      

  return {ansq, ansqd, anstau};
}


// #include <ctime>
// #include <fstream>
// #include <iostream>
// #include <vector>
// #include <streambuf>
// #include <string>
// #include <thread>


// #include "fix64_scalar.h" 
// #include "tiny_double_utils.h"
// #include "tiny_matrix3x3.h"
// #include "tiny_matrix_x.h"
// #include "tiny_mb_constraint_solver.h"
// #include "tiny_multi_body.h"
// #include "tiny_pose.h"
// #include "tiny_quaternion.h"
// #include "tiny_urdf_structures.h"
// #include "tiny_raycast.h"  
// #include "tiny_rigid_body.h"
// #include "tiny_urdf_to_multi_body.h"
// #include "tiny_vector3.h"
// #include "tiny_world.h"
// #include "examples/motion_import.h"
// #include "examples/tiny_urdf_parser.h"
// #include "tiny_metrics.h"



// typedef TinyAlgebra<Scalar, Utils> MyAlgebra;
// typedef tds::MultiBodyContactPoint<MyAlgebra> MultiBodyContactPoint;

// // MyAlgebra::Vector3 gravity( MyAlgebra::zero(),
// //                             MyAlgebra::zero(),
// //                             MyAlgebra::fraction(-980, 100));

//  auto ts = MyAlgebra::fraction(1, 100);

typedef double Scalar; 
typedef DoubleUtils Utils;

int main(){
    std::cout << "Hello World\n"; 
    std::string anymal_path = "/home/sgiles/anymal_ws/src/anymal_c_simple_description/urdf/anymal2.urdf";
    std::string plane_path = "/home/sgiles/src/tiny-differentiable-simulator/data/plane_implicit.urdf";
    auto g = Utils::fraction(-981, 100);
    TinyWorld<Scalar, Utils> world(g, false);
    // world = pd.TinyWorld(gravity_z=g, do_vis=False)
    double dt = 0.01;
    world.dt = Utils::scalar_from_double(dt);
    // world.dt = pd.Utils.scalar_from_double(dt)
    TinyUrdfParser<Scalar, Utils> parser;
    // parser = pd.TinyUrdfParser()

    UrdfToMultiBody2<Scalar, Utils> convert_tool;
    // convert_tool = pd.UrdfToMultiBody2()

    auto plane_mb = world.create_multi_body();
    // plane_mb = world.create_multi_body()
    plane_mb->m_isFloating = false;
    TinyUrdfStructures plane_urdf_data = parser.load_urdf(plane_path);
// plane_urdf_data = parser.load_urdf(str(plane_urdf_path))
    convert_tool.convert(&plane_urdf_data, &world, plane_mb);
// convert_tool.convert2(plane_urdf_data, world, plane_mb)

    auto mb = world.create_multi_body();
    mb->m_isFloating = true;
    TinyUrdfStructures urdf_data = parser.load_urdf(anymal_path);
// urdf_data = parser.load_urdf(str(anymal_urdf_path))
    convert_tool.convert(&urdf_data, &world, mb);
// convert_tool.convert2(urdf_data, world, mb)
    int n_step = 0;
    int dof_u = 12;

    world.adj_initialize(TinyVector3<Scalar, Utils>(0,0,g), n_step, dof_u);
    // vector<Scalar> q(19, 0.0), qd(18, 0.0), tau(12, 0.0);
    // q[3] = 1.0;
    // q[6] = 0.60;
    // auto fwd = forward_step(q, qd, tau, world);
    // for (auto x: fwd){
    //     for (auto y: x){
    //         std::cout << y << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    vector<Scalar> dldq(19, 0.0), dldqd(18, 0.0);


    #define N_RUNS 10000

    // q[3] = 1.0;
    // q[6] = 0.60;
    // auto fwd = forward_step(q, qd, tau, world);   
    // dldq[6] = 1.0;
    // auto ans = backward_step(q, qd, tau, dldq, dldqd, world);
    // dldq[6] = 0.0;

    std::clock_t c_start = std::clock();
    vector<Scalar> q(19, 0.0), qd(18, 0.0), tau(12, 0.0);
    q[3] = 1.0;
    q[6] = 0.55; // 0.65
    auto fwd = forward_step(q, qd, tau, world);   
    for (int i = 0; i<N_RUNS; i++){
      dldqd[i % 18] = 1.0;
      auto ans = backward_step(q, qd, tau, dldq, dldqd, world);
      dldqd[i % 18] = 0.0;
    }

    std::clock_t c_end = std::clock();


    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC / N_RUNS;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
    std::cout << "CPU time used: " << c_end-c_start << " clocks\n";
    std::cout << std::endl;
    // for (int i=0; i<q.size(); i++){
    //     ans[0][i] = (ans[0][i] - dldq[i]) / dt;
    // }
    // for (int i=0; i<qd.size(); i++){
    //     ans[1][i] = (ans[1][i] - dldqd[i]) / dt;
    // }
    // for (auto x: ans){
    //     for (auto y: x){
    //         std::cout << y << ", ";
    //     }
    //     std::cout << std::endl;
    // }


    // std::cout << ans[0][0];

    // std::vector<MultiBody<MyAlgebra> *> bodies;
    // bodies.emplace_back(&plane_mb);
    // bodies.emplace_back(&mb);
    // auto dispatcher = world.get_collision_dispatcher();
    // MultiBodyConstraintSolver<MyAlgebra> mb_solver;
    // mb_solver.pgs_iterations_ = 10;

    // // std::cout << "Initialized" << std::endl;
    // // std::cout << "Enter to run.";
    // // std::cin.get();

    // mb.q_[6].set_real(0.55);
    // mb.tau_[0].set_dual(1.);

    // // Assign q and qd to state

    // auto qd_out = mb.qd_;

    // std::clock_t c_start = std::clock();

    // forward_kinematics(mb, mb.q(), mb.qd());
    // // Assign input to tau


    // forward_dynamics(mb, world.get_gravity());
    // integrate_euler_qdd(mb, ts);

    // auto contacts = world.compute_contacts_multi_body(bodies, &dispatcher);
    // auto collisions = mb_solver.resolve_collision2(contacts[0], ts);
    // // qdd_out = [(mb.qd_ - qd_out) / ts for i in range(mb.qd.size())]
    // auto qdd_out = (mb.qd_ - qd_out) / ts;
    
    // std::clock_t c_end = std::clock();

    // // for (int i=0; i<18; i++){
    // //     auto x = qd_out[i];
    // //     std::cout << x.real() << '\t' << x.dual() << std::endl;
    // // }    
    // // for (int i=0; i<18; i++){
    // //     auto x = qdd_out[i];
    // //     std::cout << x.real() << '\t' << x.dual() << std::endl;
    // // }

    // integrate_euler(mb, ts);


    // double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    // std::cout << "collisions: " << collisions.size() << "\n";
    // std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
    // std::cout << "CPU time used: " << c_end-c_start << " clocks\n";
    // return 0;
}