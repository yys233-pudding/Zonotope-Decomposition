#ifndef ZONOOPT_EMPTYSET_HPP_
#define ZONOOPT_EMPTYSET_HPP_

/**
 * @file EmptySet.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Empty Set class for ZonoOpt library.
 * @version 1.0
 * @date 2025-09-11
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "ConZono.hpp"
#include "Zono.hpp"

namespace ZonoOpt {

/**
 * @brief Empty Set class.
 *
 * Used to facilitate set operations with trivial solutions when one of the sets is an empty set.
 */
class EmptySet final : public ConZono
{
public:

    /**
     * @brief Default constructor for EmptySet class
     *
     */
    EmptySet() = default;

    /**
     * @brief EmptySet constructor
     *
     * @param n dimension
     */
    explicit EmptySet(int n);

    HybZono* clone() const override;

    std::string print() const override;

    void constraint_reduction() override { /* do nothing */ }

    std::unique_ptr<Zono> to_zono_approx() const override { throw std::runtime_error("to_zono_approx: EmptySet"); }

protected:
    Eigen::Vector<zono_float, -1> do_optimize_over(
            const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, zono_float,
            const OptSettings&, OptSolution* solution) const override;

    Eigen::Vector<zono_float, -1> do_project_point(const Eigen::Vector<zono_float, -1>&, const OptSettings&, OptSolution* solution) const override;

    zono_float do_support(const Eigen::Vector<zono_float, -1>&, const OptSettings&, OptSolution* solution) override;

    bool do_contains_point(const Eigen::Vector<zono_float, -1>&, const OptSettings&, OptSolution*) const override;

    Box do_bounding_box(const OptSettings&, OptSolution*) override;

    bool do_is_empty(const OptSettings&, OptSolution*) const override;

    std::unique_ptr<HybZono> do_complement(zono_float delta_m, bool, const OptSettings&, OptSolution*, int, int) override;
};

}

#endif