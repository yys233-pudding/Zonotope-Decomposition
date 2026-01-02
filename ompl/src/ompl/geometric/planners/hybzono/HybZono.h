/* Author: Generated for integration */

#ifndef OMPL_GEOMETRIC_PLANNERS_HYBZONO_HYBZONO_
#define OMPL_GEOMETRIC_PLANNERS_HYBZONO_HYBZONO_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "MPC_QP.hpp"

namespace ompl
{
    namespace geometric
    {
        class HybZono : public base::Planner
        {
        public:
            HybZono(const base::SpaceInformationPtr &si);

            ~HybZono() override = default;

            void clear() override;

            void setup() override;

            base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc) override;

            void getPlannerData(base::PlannerData &data) const override;

            // Additional methods specific to HybZono
            void setMPCParameters(double horizon, double dt);

        protected:
            double horizon_ = 10.0;
            double dt_ = 0.1;
            HybZonoMPC::MPC_QP mpcSolver_;
        };
    }
}

#endif