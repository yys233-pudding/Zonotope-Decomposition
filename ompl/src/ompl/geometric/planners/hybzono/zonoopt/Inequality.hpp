#ifndef INEQUALITY_HPP_
#define INEQUALITY_HPP_

/**
 * @file Inequality.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Class definitions for zero-one inequalities.
 * @version 1.0
 * @date 2025-06-06
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <vector>

namespace ZonoOpt
{

/**
 * @brief Structure containing term in 0-1 inequality.
 *
 */
struct IneqTerm
{
    /// index of variable
    int idx; 

    /// coefficient of variable
    zono_float coeff; 
};

/**
 * @brief Enumeration to select inequality direction / use equality.
 *
 */
enum IneqType
{
    LESS_OR_EQUAL = -1, // less than or equal to
    GREATER_OR_EQUAL = 1, // greater than or equal to
    EQUAL = 0, // equal to
    LESS = -2, // strictly less than
    GREATER = 2 // strictly greater than
};

/**
 * @brief Inequality class
 *
 */
class Inequality {
public:
    /**
     * @brief Constructs inequality, must specify number of dimensions.
     *
     * @param n_dims number of dimensions in inequality
     */
    explicit Inequality(const int n_dims) : n_dims(n_dims) {}

    /**
     * @brief Constructs inequality with all quantities specified
     *
     * @param n_dims number of dimensions in inequality
     * @param terms vector of terms in inequality
     * @param rhs right-hand side of the inequality
     * @param type inequality type (e.g., less than or equal, greater than or equal, or equal)
     */
    Inequality(const int n_dims, const std::vector<IneqTerm>& terms, const zono_float rhs, const IneqType type) :
        n_dims(n_dims), terms(terms), rhs(rhs), type(type) {}

    /**
     * @brief Adds a term to the inequality
     * 
     * @param idx variable index
     * @param coeff coefficient
     */
    void add_term(const int idx, const zono_float coeff)
    {
        IneqTerm term{};
        term.idx = idx;
        term.coeff = coeff;
        this->terms.push_back(term);
    }

    /**
     * @brief Set the right hand side of the inequality
     * 
     * @param rhs right hand side value
     */
    void set_rhs(const zono_float rhs)
    {
        this->rhs = rhs;
    }

    /**
     * @brief Sets the direction of the inequality or sets it to be an equality
     * 
     * @param type inequality type (e.g., less than or equal, greater than or equal, or equal)
     */
    void set_ineq_type(const IneqType type)
    {
        this->type = type;
    }

    /**
     * @brief Get reference to terms in the inequality.
     * @return reference to terms member
     *
     */
    const std::vector<IneqTerm>& get_terms() const
    {
        return this->terms;
    }

    /**
     * @brief Get right hand side of the inequality.
     * @return right hand side value (rhs member)
     */
    zono_float get_rhs() const
    {
        return this->rhs;
    }

    /**
     * @brief Get inequality type / direction.
     * @return inequality type (type member)
     */
    IneqType get_ineq_type() const
    {
        return this->type;
    }

    /**
     * @brief Get number of dimensions for inequality.
     * @return number of dimensions (n_dims member)
     */
    int get_n_dims() const
    {
        return n_dims;
    }

private:
    /// number of dimensions
    int n_dims;

    /// terms
    std::vector<IneqTerm> terms;

    /// right-hand side of the inequality
    zono_float rhs = 0;

    /// type of the inequality (less than or equal, greater than or equal, or equal)
    IneqType type = EQUAL;
};

}

#endif