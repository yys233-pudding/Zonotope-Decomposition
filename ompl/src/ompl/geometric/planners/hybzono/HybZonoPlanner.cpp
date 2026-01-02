#include "HybZono.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/tools/config/SelfConfig.h"
#include "ompl/base/spaces/RealVectorStateSpace.h"
#include <Eigen/Dense>

namespace ompl
{
    namespace geometric
    {
        HybZono::HybZono(const base::SpaceInformationPtr &si) : base::Planner(si, "HybZono")
        {
            specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
            specs_.directed = true;
        }

        void HybZono::clear()
        {
            Planner::clear();
        }

        void HybZono::setup()
        {
            Planner::setup();
            // Setup MPC
            int dim = si_->getStateSpace()->getDimension();
            if (dim != 2) {
                OMPL_WARN("HybZono planner currently only supports 2D spaces");
            }

            // Dynamics: simple integrator x' = u, discretized as x_{k+1} = x_k + u_k * dt
            Eigen::MatrixXd A(2,2);
            A << 1, 0, 0, 1;
            Eigen::MatrixXd B(2,1);
            B << dt_, 0;

            mpcSolver_.set_dyn_matrices(A.sparseView(), B.sparseView());

            // Costs
            Eigen::MatrixXd Q(2,2);
            Q << 1, 0, 0, 1;
            Eigen::MatrixXd R(1,1);
            R << 0.1;
            Eigen::MatrixXd P(2,2);
            P << 1, 0, 0, 1;
            mpcSolver_.set_stage_cost(Q.sparseView(), R.sparseView());
            mpcSolver_.set_terminal_cost(P.sparseView());

            // Constraints: x in [-10,10], u in [-1,1]
            Eigen::MatrixXd Ax_ineq(4,2);
            Ax_ineq << 1, 0, -1, 0, 0, 1, 0, -1;
            Eigen::VectorXd bx_ineq(4);
            bx_ineq << 10, 10, 10, 10;
            mpcSolver_.set_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);

            Eigen::MatrixXd Au_ineq(2,1);
            Au_ineq << 1, -1;
            Eigen::VectorXd bu_ineq(2);
            bu_ineq << 1, 1;
            mpcSolver_.set_input_inequality_constraints(Au_ineq.sparseView(), bu_ineq);

            // Terminal constraints same as state
            mpcSolver_.set_terminal_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);

            // Settings
            HybZonoMPC::MPC_Settings settings;
            settings.n_horizon = static_cast<int>(horizon_ / dt_);
            mpcSolver_.set_mpc_settings(settings);

            // Build controller
            mpcSolver_.build_controller();
        }

        base::PlannerStatus HybZono::solve(const base::PlannerTerminationCondition &ptc)
        {
            auto *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());
            if (!goal)
            {
                OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
                return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
            }

            // Get start state
            if (pdef_->getStartStateCount() == 0)
            {
                OMPL_ERROR("%s: No start state specified", getName().c_str());
                return base::PlannerStatus::INVALID_START;
            }
            const base::State *start_state = pdef_->getStartState(0);

            // Set initial condition
            Eigen::VectorXd x0(2);
            x0 << start_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
                  start_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[1];

            // Sample goal
            base::State *goal_state = si_->allocState();
            goal->sampleGoal(goal_state);
            Eigen::VectorXd goal_vec(2);
            goal_vec << goal_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
                        goal_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[1];
            si_->freeState(goal_state);

            // Set reference to goal
            int n_steps = static_cast<int>(horizon_ / dt_);
            std::vector<Eigen::VectorXd> x_ref(n_steps, goal_vec);

            // Control
            auto [u, success] = mpcSolver_.control(x0, x_ref);
            if (!success) {
                OMPL_WARN("%s: MPC control failed", getName().c_str());
                return base::PlannerStatus::TIMEOUT;
            }

            // Get trajectories
            auto [x_vec, u_vec] = mpcSolver_.get_trajectory();

            // Build path
            auto path = std::make_shared<PathGeometric>(si_);
            path->append(start_state);
            for (const auto& x : x_vec) {
                auto state = si_->allocState();
                state->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = x(0);
                state->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = x(1);
                path->append(state);
            }

            // Add goal
            goal_state = si_->allocState();
            goal_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = goal_vec(0);
            goal_state->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = goal_vec(1);
            path->append(goal_state);
            si_->freeState(goal_state);

            pdef_->addSolutionPath(path);
            return base::PlannerStatus::EXACT_SOLUTION;
        }

        void HybZono::getPlannerData(base::PlannerData &data) const
        {
            Planner::getPlannerData(data);
        }

        void HybZono::setMPCParameters(double horizon, double dt)
        {
            horizon_ = horizon;
            dt_ = dt;
        }
    }
}