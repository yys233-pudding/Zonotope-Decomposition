#ifndef ZONOOPT_SOLVER_DATA_STRUCUTURES_HPP_
#define ZONOOPT_SOLVER_DATA_STRUCUTURES_HPP_

/**
 * @file SolverDataStructures.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Optimization settings and solution data structures for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <Eigen/Dense>
#include <sstream>

namespace ZonoOpt
{

/**
 * @brief Settings for optimization routines in ZonoOpt library.
 * 
 */
struct OptSettings
{
    // general settings

    /// display optimization progress
    bool verbose = false; 

    /// print every verbose_interval iterations
    int verbosity_interval = 100;

    /// max time for optimization
    double t_max = std::numeric_limits<double>::max(); 

    // ADMM settings

    /// max admm iterations
    int k_max_admm = 5000; 

    /// admm penalty parameter, higher prioritizes feasibility during iterations, lower prioritizes optimality
    zono_float rho = static_cast<zono_float>(10.0); 

    /// dual convergence tolerance 
    zono_float eps_dual = static_cast<zono_float>(1e-2);

    /// primal convergence tolerance
    zono_float eps_prim = static_cast<zono_float>(1e-3);

    /// check infeasibility every k_inf_check iterations
    int k_inf_check = 10; 

    /// use infinity norm for convergence check (if false, scaled 2-norm is used)
    bool inf_norm_conv = true; 

    /// flag to use interval contractor for constraint tightening / implication
    bool use_interval_contractor = true; 

    /// number of interval contractor iterations
    int contractor_iter = 1; 

    // mixed integer settings

    /// 0: best first, 1: best dive
    int search_mode = 0;

    /// flag to perform solution polishing
    bool polish = true;

    /// dual residual convergence tolerance during branch and bound and search
    zono_float eps_dual_search = static_cast<zono_float>(1e-1); 

    /// primal residual convergence tolerance during branch and bound and search
    zono_float eps_prim_search = static_cast<zono_float>(1e-2);

    /// relative convergence tolerance
    zono_float eps_r = static_cast<zono_float>(1e-2); 

    /// absolute convergence tolerance
    zono_float eps_a = static_cast<zono_float>(1e-1); 

    /// max number of branch-and-bound iterations
    int k_max_bnb = static_cast<int>(1e5);

    /// max threads for branch and bound
    int n_threads_bnb = 4;

    /// terminate if more than this many nodes are in branch and bound queue
    int max_nodes = static_cast<int>(1e5);

    /// when applying interval contractor in branch and bound, this is how deep to search the constraint tree for affected variables
    int contractor_tree_search_depth = 10;

    // validity check
    /**
     * @brief Checks whether settings struct is valid
     * 
     * @return validity boolean
     */
    bool settings_valid() const
    {
        const bool general_valid = t_max > 0 && verbosity_interval > 0;
        const bool admm_valid = (rho > 0 && k_max_admm > 0 &&
                eps_dual >= 0 && eps_prim >= 0 && k_inf_check >= 0 && contractor_iter > 0);
        const bool mi_valid = (eps_r >= 0 && eps_a >= 0 && eps_dual_search > 0 && eps_prim_search > 0 &&
                k_max_bnb > 0 && n_threads_bnb >= 0 && max_nodes > 0 && contractor_tree_search_depth > 0 &&
               (search_mode == 0 || search_mode == 1));

        return (general_valid && admm_valid && mi_valid);
    }

    /**
     * @brief displays settings as string
     *
     * @return string
     */
    std::string print() const
    {
        std::stringstream ss;
        ss << "OptSettings structure: " << std::endl;
        ss << "  verbose: " << (verbose ? "true" : "false") << std::endl;
        ss << "  verbosity_interval: " << verbosity_interval << std::endl;
        ss << "  t_max: " << t_max << std::endl;
        ss << "  k_max_admm: " << k_max_admm << std::endl;
        ss << "  rho: " << rho << std::endl;
        ss << "  eps_dual: " << eps_dual << std::endl;
        ss << "  eps_prim: " << eps_prim << std::endl;
        ss << "  k_inf_check: " << k_inf_check << std::endl;
        ss << "  inf_norm_conv: " << (inf_norm_conv ? "true" : "false") << std::endl;
        ss << "  use_interval_contractor: " << (use_interval_contractor ? "true" : "false") << std::endl;
        ss << "  contractor_iter: " << contractor_iter << std::endl;
        ss << "  search_mode: " << search_mode << std::endl;
        ss << "  polish: " << polish << std::endl;
        ss << "  eps_dual_search: " << eps_dual_search << std::endl;
        ss << "  eps_prim_search: " << eps_prim_search << std::endl;
        ss << "  eps_r: " << eps_r << std::endl;
        ss << "  eps_a: " << eps_a << std::endl;
        ss << "  k_max_bnb: " << k_max_bnb << std::endl;
        ss << "  n_threads_bnb: " << n_threads_bnb << std::endl;
        ss << "  max_nodes: " << max_nodes << std::endl;
        ss << "  contractor_tree_search_depth: " << contractor_tree_search_depth;
        return ss.str();
    }
};

/**
 * @brief Solution data structure for optimization routines in ZonoOpt library.
 * 
 */
struct OptSolution
{
    // general

    /// solution vector
    Eigen::Vector<zono_float, -1> z;

    /// objective
    zono_float J = -std::numeric_limits<zono_float>::infinity();

    /// time to compute solution
    double run_time = 0.0;

    /// time to factorize matrices and run interval contractors
    double startup_time = 0.0;

    /// number of iterations
    int iter = 0;

    /// true if optimization has converged
    bool converged = false;

    /// true if optimization problem is provably infeasible
    bool infeasible = false;

    // admm-specific

    /// ADMM primal variable, approximately equal to z when converged
    Eigen::Vector<zono_float, -1> x;

    /// ADMM dual variable
    Eigen::Vector<zono_float, -1> u;

    /// primal residual, corresponds to feasibility
    zono_float primal_residual = std::numeric_limits<zono_float>::infinity();

    /// dual residual, corresponds to optimality
    zono_float dual_residual = std::numeric_limits<zono_float>::infinity();

    /**
     * @brief displays solution as string
     *
     * @return string
     */
    std::string print() const
    {
        std::stringstream ss;
        ss << "OptSolution structure:" << std::endl;
        ss << "  z: vector of length " << z.size() << std::endl;
        ss << "  J: " << J << std::endl;
        ss << "  run_time: " << run_time << std::endl;
        ss << "  startup_time: " << startup_time << std::endl;
        ss << "  iter: " << iter << std::endl;
        ss << "  converged: " << (converged ? "true" : "false") << std::endl;
        ss << "  infeasible: " << (infeasible ? "true" : "false") << std::endl;
        ss << "  x: vector of length " << x.size() << std::endl;
        ss << "  u: vector of length " << u.size() << std::endl;
        ss << "  primal_residual: " << primal_residual << std::endl;
        ss << "  dual_residual: " << dual_residual << std::endl;
        return ss.str();
    }
};


} // end namespace ZonoOpt

#endif