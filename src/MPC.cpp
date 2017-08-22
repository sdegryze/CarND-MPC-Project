#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the number of time points of the reference horizon and the duration of 1 timepoint
size_t N = 40;
double dt = 0.025;
// Set the reference speed
double ref_v = 70;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // `fg` is a vector of the cost constraints
    // `vars` is a vector of variable values (6 state vars & 2 actuators vars for each of N timesteps)
    
    // Reference State Cost
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
    
    // Remember here to sum the cost for all time steps N
    for (int t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start + t], 2);
      // The next line is important to make the car actually drive, with this term we force
      // solutions that have a speed close to the reference speed
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    
    // The next part is analogous to regularization: we are penalizing large values of the
    // actuators, to achieve a more smoother ride
    // I found that this was especially important for the steering angle delta and less so
    // for the throttle a
    
    // Note that there are only N-1 actuators, since the state at t=0 is the given (current) state.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += 70000 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + t], 2);
    }
    
    // Finally, penalize abrupt changes in actuation
    for (int t = 0; t < N - 2; t++) {
      fg[0] += 10000 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
    
    
    // Initial constraints - this sets the first variables to the (known) initial state (e.g.,
    // current position of the car)
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    
    // In the rest of the constraints, we implement the motion model
    // Note that we are starting from 1, since 0 indicates the current state
    for (int t = 1; t < N; t++) {
      // Note that we are not using the delta and a value at (relative) time 1
      // Hence, we have 1 iteration less values for the actuators than for the state
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];
      
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];
      
      // Evaluate the 3rd degree polynomial of the waypoints (i.e., the reference trajectory) at x0
      // this term is needed to calculate the cross-track error for future timepoints
      AD<double> fx0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      // Evaluate the 1st derivative of the 3rd degree polynomial at x0
      AD<double> fpx0 = coeffs[1] + 2. * coeffs[2] * x0 + 3. * coeffs[3] * x0 * x0;
      // Calculate the angle of the tangent to the reference trajectory at x0
      AD<double> psi_des0 = CppAD::atan(fpx0);
      
      
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.
      // Implementation of the motion model and the cross-track error
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 / Lf * delta0 * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - (fx0 - y0 - v0 * CppAD::sin(epsi0) * dt);
      fg[1 + epsi_start + t] = epsi1 - (psi0 - psi_des0 + v0 / Lf * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

MPC_OUTPUT MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  
  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // See above for why the actuators are multiplied by (N - 1): the state in
  // the last timestep cannot be updated anymore with any actuations in that
  // time step anyhow
  size_t n_vars = state.size() * N + 2 * (N - 1);
  // Set the number of constraints (+ 1 for the cost function)
  size_t n_constraints = state.size() * N + 1;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Now set the initial values for the independent variables
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // Set lower and upper limits for the independent variables.
  // by setting all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  
  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  

  // Set lower and upper limits for the constraints
  // Should be both 0 besides initial state
  // Note that there is a difference of 1 with the index in fg. So, constraint with index 2
  // corresponds to fg with index 1
  // This is because the cost function (to be minimized) is put in fg with index 0 and the
  // cost function does not require constraints, obviously
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;
  
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  MPC_OUTPUT output;
  
  vector<double> xpts(N);
  vector<double> ypts(N);
  for (int i = 0; i < N; i++) {
    xpts[i] = solution.x[x_start + i];
    ypts[i] = solution.x[y_start + i];
  }

  // Note that we want the values of delta and a for the first timestep
  // For the other variables, we don't, since that would be simply the initial state
  output.vars = {solution.x[x_start + 1],   solution.x[y_start + 1],
                 solution.x[psi_start + 1], solution.x[v_start + 1],
                 solution.x[cte_start + 1], solution.x[epsi_start + 1],
                 solution.x[delta_start],   solution.x[a_start]};
  
  // sending these points as output is needed to show the MPC predicted trajectory
  output.xpts = xpts;
  output.ypts = ypts;

  return output;
}
