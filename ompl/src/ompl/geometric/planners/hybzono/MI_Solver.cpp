#include "ZonoOpt.hpp"

namespace ZonoOpt::detail
{
    MI_Solver::MI_Solver(const MI_data& data) : data(data), node_queue(comp)
    {
        // check settings validity
        if (!this->data.admm_data->settings.settings_valid())
        {
            throw std::invalid_argument("MI_solver setup: invalid settings.");
        }

        // check number of threads is valid
        if (this->data.admm_data->settings.n_threads_bnb > static_cast<int>(std::thread::hardware_concurrency()) - 1)
        {
            std::stringstream ss;
            ss << "MI_solver setup: number of threads for branch and bound (" << this->data.admm_data->settings.
                n_threads_bnb
                << ") + convergence monitoring (1) exceeds available threads (" << std::thread::hardware_concurrency()
                << ").";
            throw std::invalid_argument(ss.str());
        }
    }

    OptSolution MI_Solver::solve()
    {
        this->multi_sol = false;
        auto sol = solver_core();
        return std::get<OptSolution>(sol);
    }

    std::pair<std::vector<OptSolution>, OptSolution> MI_Solver::multi_solve(const int max_sols)
    {
        this->multi_sol = true;
        auto sol = solver_core(max_sols);
        return std::get<std::pair<std::vector<OptSolution>, OptSolution>>(sol);
    }

    std::unique_ptr<Node, MI_Solver::NodeDeleter> MI_Solver::make_node(const std::shared_ptr<ADMM_data>& data)
    {
        void* mem = pool.allocate(sizeof(Node), alignof(Node));
        Node* node = new(mem) Node(data);
        return {node, NodeDeleter(&pool)};
    }

    std::unique_ptr<Node, MI_Solver::NodeDeleter> MI_Solver::clone_node(const std::unique_ptr<Node, NodeDeleter>& other)
    {
        void* mem = pool.allocate(sizeof(Node), alignof(Node));
        Node* node = new(mem) Node(*other);
        return {node, NodeDeleter(&pool)};
    }

    std::variant<OptSolution, std::pair<std::vector<OptSolution>, OptSolution>> MI_Solver::solver_core(int max_sols)
    {
        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        double run_time;

        // set flags
        this->done = false;
        this->converged = false;

        // verbosity
        std::stringstream ss;
        if (this->data.admm_data->settings.verbose)
        {
            if (this->multi_sol)
            {
                ss << "Finding up to " << max_sols << " solutions to MIQP with " << this->data.admm_data->n_x <<
                    " variables and "
                    << this->data.admm_data->n_cons << " constraints using " << this->data.admm_data->settings.
                    n_threads_bnb << " branch-and-bound threads.";
            }
            else
            {
                ss << "Solving MIQP problem with " << this->data.admm_data->n_x << " variables and "
                    << this->data.admm_data->n_cons << " constraints using " << this->data.admm_data->settings.
                    n_threads_bnb << " branch-and-bound threads.";
            }
            print_str(ss);
        }

        auto return_infeasible_solution = [this, &start
            ]() -> std::variant<OptSolution, std::pair<std::vector<OptSolution>, OptSolution>>
        {
            OptSolution infeasible_solution;
            infeasible_solution.infeasible = true;
            infeasible_solution.run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<
                std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());
            infeasible_solution.startup_time = infeasible_solution.run_time;
            infeasible_solution.iter = 0;
            infeasible_solution.converged = false;
            infeasible_solution.J = std::numeric_limits<zono_float>::infinity();
            infeasible_solution.z = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.x = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.u = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.primal_residual = std::numeric_limits<zono_float>::infinity();
            infeasible_solution.dual_residual = std::numeric_limits<zono_float>::infinity();

            if (this->multi_sol)
            {
                return std::make_pair(this->solutions.get(), infeasible_solution);
            }
            else
            {
                return infeasible_solution;
            }
        };

        // add root node
        this->bnb_data.reset(this->data.admm_data->clone()); // init
        this->bnb_data->settings.verbose = false;
        this->bnb_data->settings.eps_dual = this->data.admm_data->settings.eps_dual_search;
        this->bnb_data->settings.eps_prim = this->data.admm_data->settings.eps_prim_search;
        std::unique_ptr<Node, NodeDeleter> root = this->make_node(this->bnb_data);
        if (!root->run_contractor()) // check for infeasibility during setup via interval contractor
        {
            return return_infeasible_solution();
        }

        // log startup time
        double startup_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        if (this->data.admm_data->settings.verbose)
        {
            ss << "Startup time = " << startup_time << " sec";
            print_str(ss);
        }

        // solve root relaxation
        root->solve();
        if (root->solution.infeasible)
        {
            return return_infeasible_solution();
        }

        // start threads
        std::vector<std::thread> bnb_threads;
        for (int i = 0; i < this->data.admm_data->settings.n_threads_bnb; i++)
        {
            bnb_threads.emplace_back([this]() { worker_loop(); });
        }

        // push root to node queue
        this->push_node(std::move(root));

        // loop and check for exit conditions
        int print_iter = 0;
        if (this->data.admm_data->settings.verbose)
        {
            ss << std::endl << std::setw(10) << "Iter" << std::setw(10) << "Queue" << std::setw(10) <<
                "Threads" << std::setw(10) << "Time [s]" << std::setw(10) << "J_min" << std::setw(10) <<
                "J_max" << std::setw(10) << "Gap [%]" << std::setw(10) << "Feasible" << std::setw(10);
            print_str(ss);
        }

        while (!this->done)
        {
            // check for timeout
            run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());
            if (run_time > this->data.admm_data->settings.t_max)
            {
                this->done = true;
            }

            // check for max nodes
            int queue_size;
            {
                std::lock_guard<std::mutex> lock(pq_mtx);
                queue_size = static_cast<int>(this->node_queue.size());
            }
            if (queue_size > this->data.admm_data->settings.max_nodes)
            {
                this->done = true;
            }

            // check for max iterations
            if (this->iter >= this->data.admm_data->settings.k_max_bnb)
            {
                this->done = true;
            }

            // check for convergence

            // get lower bound / check if there are no nodes remaining
            zono_float J_min = -std::numeric_limits<zono_float>::infinity();
            {
                std::pair<zono_float, bool> J_min_threads_pair;
                std::lock_guard<std::mutex> lock(pq_mtx);
                J_min_threads_pair = this->J_threads.get_min(); // lower bound from active threads

                if (this->node_queue.empty())
                {
                    if (!J_min_threads_pair.second)
                    {
                        this->done = true; // no nodes remaining
                        this->converged = true;
                    }
                    else
                    {
                        J_min = J_min_threads_pair.first;
                    }
                }
                else
                {
                    J_min = std::min(this->node_queue.top()->solution.J, J_min_threads_pair.first);
                }
            }

            // check for convergence based on lower and upper bounds
            zono_float gap_percent = 0;
            if (!this->multi_sol)
            {
                const zono_float gap = std::abs(this->J_max - J_min);
                gap_percent = std::abs(this->J_max - J_min) / std::abs(this->J_max);
                if ((gap_percent < this->data.admm_data->settings.eps_r) || (gap < this->data.admm_data->settings.
                    eps_a))
                {
                    this->done = true;
                    this->converged = true;
                }
            }
            else // check based on number of solutions
            {
                if (this->solutions.size() >= static_cast<size_t>(max_sols))
                {
                    this->done = true;
                    this->converged = true;
                }
            }

            // verbosity
            if (this->data.admm_data->settings.verbose && (this->iter >= print_iter))
            {
                size_t n_threads = this->J_threads.size();
                ss << std::setw(10) << this->iter << std::setw(10) << queue_size << std::setw(10)
                    << n_threads << std::setw(10) << run_time << std::setw(10)
                    << J_min << std::setw(10) << this->J_max << std::setw(10)
                    << gap_percent * 100.0f << std::setw(10)
                    << (this->feasible ? "true" : "false") << std::endl;
                print_str(ss);
                print_iter += this->data.admm_data->settings.verbosity_interval;
            }
        }

        // clean up
        pq_cv_bnb.notify_all(); // notify all threads to stop waiting
        {
            std::lock_guard<std::mutex> lock(pq_mtx);
            this->node_queue.clear();
        }
        this->J_threads.clear();

        for (auto& thread : bnb_threads)
        {
            if (thread.joinable()) thread.join();
        }

        this->node_queue.clear(); // nodes need to be freed before pool goes out of scope to avoid race condition

        // assemble solution
        OptSolution solution;
        solution.z = this->z.get();
        solution.J = this->J_max;
        solution.run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        solution.startup_time = startup_time;
        solution.iter = this->iter;
        solution.converged = this->converged;
        solution.infeasible = !this->feasible;

        solution.x = this->x.get();
        solution.u = this->u.get();
        solution.primal_residual = this->primal_residual;
        solution.dual_residual = this->dual_residual;

        if (this->multi_sol)
        {
            if (this->data.admm_data->settings.verbose)
            {
                ss << "Found " << this->solutions.size() << " solutions to MIQP problem.";
                print_str(ss);
            }

            return std::make_pair(solutions.get(), solution);
        }
        else
        {
            // verbosity
            if (this->data.admm_data->settings.verbose)
            {
                ss << "Converged = " << (this->converged ? "true" : "false") << ", Feasible = " <<
                    (this->feasible ? "true" : "false") << ", Iterations = " << this->iter << ", Solve time = " <<
                    solution.run_time << " sec, Objective = " << solution.J << ", Average QP iterations = " <<
                    static_cast<double>(this->qp_iter) / static_cast<double>(this->iter)
                    << ", Average solve time = " << this->total_run_time.get() / static_cast<double>(this->iter) <<
                    " sec, Average startup time = "
                    << this->total_startup_time.get() / static_cast<double>(this->iter) << " sec" << std::endl;
                print_str(ss);
            }

            return solution;
        }
    }

    void MI_Solver::solve_and_branch(const std::unique_ptr<Node, NodeDeleter>& node)
    {
        // objective prior to solving
        const zono_float J_min_prior = node->solution.J;

        // solve node
        node->solve(&this->done);
        if (this->done) return;

        // cleanup function
        auto cleanup = [&, this]()
        {
            // remove J from J_threads vector
            this->J_threads.remove(J_min_prior);

            // increment nodes evaluated and logging info
            ++this->iter;
            this->qp_iter += node->solution.iter;
            this->total_run_time += node->solution.run_time;
            this->total_startup_time += node->solution.startup_time;
        };

        // return if infeasible, not converged, or no optimal solution exists in branch
        if (!(node->solution.infeasible || !node->solution.converged || (node->solution.J > this->J_max && !this->
            multi_sol)))
        {
            // check if node is integer feasible
            if (is_integer_feasible(node->solution.z.segment(this->data.idx_b.first, this->data.idx_b.second)))
            {
                // rerun with refined tolerance
                if (this->data.admm_data->settings.polish)
                {
                    node->update_convergence_tolerances(this->data.admm_data->settings.eps_prim,
                                                        this->data.admm_data->settings.eps_dual);
                    node->warmstart(node->solution.z, node->solution.u);
                    node->solve();
                }

                // make sure still integer feasible after refining
                if (!is_integer_feasible(node->solution.z.segment(this->data.idx_b.first, this->data.idx_b.second)))
                {
                    branch_most_frac(node);
                    cleanup(); // cleanup function
                    return;
                }

                // make sure return conditions are not met after refining
                if (node->solution.infeasible || !node->solution.converged || (node->solution.J > this->J_max && !this->
                    multi_sol))
                {
                    cleanup(); // cleanup function
                    return;
                }
                if (this->multi_sol) // store solution if doing multisol
                {
                    std::function<bool(const OptSolution&, const OptSolution&)> compare_eq = [this
                        ](const OptSolution& a, const OptSolution& b)
                    {
                        return this->check_bin_equal(a, b);
                    };
                    if (!this->solutions.contains(node->solution, compare_eq)) // new solution
                        this->solutions.push_back(node->solution);
                }
                if (node->solution.J < this->J_max - zono_eps) // check if node is better than current best
                {
                    // update incumbent
                    this->J_max = node->solution.J;
                    this->x.set(node->solution.x);
                    this->z.set(node->solution.z);
                    this->u.set(node->solution.u);
                    this->primal_residual = node->solution.primal_residual;
                    this->dual_residual = node->solution.dual_residual;
                    this->feasible = true;

                    // prune
                    if (!this->multi_sol) this->prune(node->solution.J);
                }
            }
            else
            {
                branch_most_frac(node);
            }
        }

        cleanup(); // cleanup function
    }

    bool MI_Solver::is_integer_feasible(const Eigen::Ref<const Eigen::Vector<zono_float, -1>> xb) const
    {
        const zono_float low = this->data.zero_one_form ? zero : -one;
        constexpr zono_float high = 1;

        for (int i = 0; i < xb.size(); i++)
        {
            if ((std::abs(xb(i) - high) > zono_eps) && (std::abs(xb(i) - low) > zono_eps))
                return false;
        }
        return true;
    }

    void MI_Solver::branch_most_frac(const std::unique_ptr<Node, NodeDeleter>& node)
    {
        // must be at least 1 binary variable
        if (this->data.idx_b.second <= 0)
            return;

        const zono_float low = this->data.zero_one_form ? zero : -one;
        constexpr zono_float high = 1;

        // round and find most fractional variable
        const Eigen::Array<zono_float, -1, 1> xb = node->solution.z.segment(
            this->data.idx_b.first, this->data.idx_b.second).array();
        Eigen::Array<zono_float, -1, 1> l(xb.size());
        Eigen::Array<zono_float, -1, 1> u(xb.size());
        l.setConstant(low);
        u.setConstant(high);
        const Eigen::Array<zono_float, -1, 1> d_l = (xb - l).abs();
        const Eigen::Array<zono_float, -1, 1> d_u = (xb - u).abs();
        const Eigen::Array<zono_float, -1, 1> d = d_l.min(d_u); // distance to rounded value
        int idx_most_frac;
        d.maxCoeff(&idx_most_frac); // index of most fractional variable
        idx_most_frac = this->data.idx_b.first + idx_most_frac; // convert to original index

        // branch on most fractional variable
        std::unique_ptr<Node, NodeDeleter> left = this->clone_node(node);
        std::unique_ptr<Node, NodeDeleter> right = this->clone_node(node);

        // branches
        const bool left_inf = !left->fix_bound(idx_most_frac, low);
        const bool right_inf = !right->fix_bound(idx_most_frac, high);

        // warm start
        left->warmstart(node->solution.z, node->solution.u);
        right->warmstart(node->solution.z, node->solution.u);

        switch (this->data.admm_data->settings.search_mode)
        {
        case (0):
            {
                // best first: push both nodes to queue
                if (!left_inf) this->push_node(std::move(left));
                if (!right_inf) this->push_node(std::move(right));
                break;
            }
        case (1):
            {
                // best dive: push worse nodes to queue, solve better node
                if (left_inf && right_inf) // both branches infeasible
                {
                    return; // nothing to do
                }
                else if (left_inf)
                {
                    this->solve_and_branch(right);
                }
                else if (right_inf)
                {
                    this->solve_and_branch(left);
                }
                else // both branches feasible
                {
                    if (left->get_box().width() > right->get_box().width()) // left is worse
                    {
                        this->push_node(std::move(left));
                        this->solve_and_branch(right);
                    }
                    else // right is worse
                    {
                        this->push_node(std::move(right));
                        this->solve_and_branch(left);
                    }
                }
                break;
            }
        default:
            {
                std::stringstream ss;
                ss << "MI_ADMM_solver: unknown search mode " << this->data.admm_data->settings.search_mode;
                throw std::runtime_error(ss.str());
            }
        }
    }

    void MI_Solver::worker_loop()
    {
        while (!this->done)
        {
            std::unique_ptr<Node, NodeDeleter> node = make_node(this->bnb_data);
            {
                std::unique_lock<std::mutex> lock(pq_mtx);
                pq_cv_bnb.wait(lock, [this]() { return this->done || !this->node_queue.empty(); });
                if (this->done) return;
                node = this->node_queue.pop_top();
                this->J_threads.add(node->solution.J);
                // add J to J_threads vector, need to do this before releasing lock
            }
            solve_and_branch(node);
        }
    }

    void MI_Solver::push_node(std::unique_ptr<Node, NodeDeleter>&& node)
    {
        std::unique_lock<std::mutex> lock(pq_mtx);
        this->node_queue.push(std::move(node));
        pq_cv_bnb.notify_one();
    }

    void MI_Solver::prune(const zono_float J_max)
    {
        // create node with J = J_max
        const std::unique_ptr<Node, NodeDeleter> n = this->make_node(this->bnb_data);
        n->solution.J = J_max;
        {
            std::lock_guard<std::mutex> lock(pq_mtx);
            this->node_queue.prune(n);
        }
    }

    bool MI_Solver::check_bin_equal(const OptSolution& sol1, const OptSolution& sol2) const
    {
        return (sol1.z.segment(this->data.idx_b.first, this->data.idx_b.second) - sol2.z.segment(
                   this->data.idx_b.first, this->data.idx_b.second))
               .cwiseAbs().maxCoeff() < zono_eps;
    }
}
