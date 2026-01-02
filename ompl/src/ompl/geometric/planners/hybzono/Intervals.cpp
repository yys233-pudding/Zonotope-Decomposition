#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    Box::Box(const size_t size)
    {
        x_lb.resize(static_cast<Eigen::Index>(size));
        x_ub.resize(static_cast<Eigen::Index>(size));
    }

    Box::Box(const std::vector<Interval>& vals)
    {
        x_lb.resize(static_cast<Eigen::Index>(vals.size()));
        x_ub.resize(static_cast<Eigen::Index>(vals.size()));

        for (Eigen::Index i=0; i<static_cast<Eigen::Index>(vals.size()); i++)
        {
            this->x_lb(i) = vals[i].y_min();
            this->x_ub(i) = vals[i].y_max();
        }
    }

    Box::Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub)
    {
        if (x_lb.size() != x_ub.size())
            throw std::invalid_argument("x_l and x_u must have the same size");
        this->x_lb = x_lb;
        this->x_ub = x_ub;
    }

    Box& Box::operator=(const Box& other)
    {
        if (this != &other)
        {
            this->x_lb = other.x_lb;
            this->x_ub = other.x_ub;
        }
        return *this;
    }

    Box::Box(const Box& other)
    {
        this->x_lb = other.x_lb;
        this->x_ub = other.x_ub;
    }

    IntervalView Box::operator[](const size_t i)
    {
        if (i >= static_cast<size_t>(x_lb.size()))
            throw std::out_of_range("Index out of range");
        return IntervalView(&x_lb(static_cast<Eigen::Index>(i)), &x_ub(static_cast<Eigen::Index>(i)));
    }

    Interval Box::operator[](const size_t i) const
    {
        if (i >= static_cast<size_t>(x_lb.size()))
            throw std::out_of_range("Index out of range");
        return Interval(x_lb(static_cast<Eigen::Index>(i)), x_ub(static_cast<Eigen::Index>(i)));
    }

    size_t Box::size() const
    {
        return x_lb.size();
    }

    void Box::project(Eigen::Ref<Eigen::Vector<zono_float, -1>> x) const
    {
        if (x.size() != x_lb.size())
            throw std::invalid_argument("x must have the same size as the Box");
        x = x.cwiseMax(x_lb).cwiseMin(x_ub);
    }

    Box* Box::clone() const
    {
        return new Box(*this);
    }

    zono_float Box::width() const
    {
        zono_float w = 0;
        for (Eigen::Index i=0; i<x_lb.size(); i++)
        {
            w += x_ub(i) - x_lb(i);
        }
        return w;
    }

    Eigen::Vector<zono_float, -1> Box::center() const
    {
        Eigen::Vector<zono_float, -1> c (this->size());
        for (Eigen::Index i=0; i<static_cast<Eigen::Index>(this->size()); i++)
        {
            c(i) = (*this)[i].center();
        }
        return c;
    }

    Box Box::operator+(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box addition: inconsistent dimensions");
        Box out =  *this;
        for (size_t i=0; i<this->size(); ++i)
        {
            out[i] = (*this)[i] + other[i];
        }
        return out;
    }

    Box Box::operator-(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box subtraction: inconsistent dimensions");
        Box out =  *this;
        for (size_t i=0; i<this->size(); ++i)
        {
            out[i] = (*this)[i] - other[i];
        }
        return out;
    }

    Box Box::operator*(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box multiplication: inconsistent dimensions");
        Box out =  *this;
        for (size_t i=0; i<this->size(); ++i)
        {
            out[i] = (*this)[i] * other[i];
        }
        return out;
    }

    Box Box::operator*(zono_float alpha) const
    {
        Box out = *this;
        for (size_t i=0; i<this->size(); ++i)
        {
            out[i] = (*this)[i]*alpha;
        }
        return out;
    }

    Box Box::operator/(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box division: inconsistent dimensions");
        Box out = *this;
        for (size_t i=0; i<this->size(); ++i)
        {
            out[i] = (*this)[i] / other[i];
        }
        return out;
    }

    bool Box::contract(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const Eigen::Vector<zono_float, -1>& b, const int iter)
    {
        if (iter <= 0)
            throw std::invalid_argument("iter must be positive");

        // contract over all constraints
        std::set<int> constraints;
        for (int i=0; i<A.rows(); i++)
        {
            constraints.insert(i);
        }

        // run contractor
        return contract_helper(A, b, iter, constraints);
    }

    bool Box::contract_subset(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm, const Eigen::Vector<zono_float, -1>& b, int iter,
                             const Eigen::SparseMatrix<zono_float>& A, const std::set<int>& inds, int tree_search_depth)
    {
        if (iter <= 0)
            throw std::invalid_argument("iter must be positive");
        if (A_rm.rows() != A.rows() || A_rm.cols() != A.cols())
            throw std::invalid_argument("A_rm must equal A");

        // get affected constraints
        std::set<int> all_constraints, all_vars;

        get_vars_cons(A, A_rm, all_constraints, all_vars, inds, 0, tree_search_depth);

        // run contractor
        return contract_helper(A_rm, b, iter, all_constraints);
    }

    Box Box::linear_map(const Eigen::Matrix<zono_float, -1, -1>& A) const
    {
        // input handling
        if (A.cols() != x_lb.size())
            throw std::invalid_argument("Matrix A must have the same number of columns as the size of the Box");

        // declare
        Box y(A.rows());

        // linear map
        for (int i=0; i<A.rows(); i++)
        {
            y[i] = Interval(0, 0);
            for (int j=0; j<A.cols(); j++)
            {
                y[i].add_assign(y[i], ((*this)[j]*A(i, j)).as_view());
            }
        }
        return y;
    }

    Box Box::linear_map(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const
    {
        // input handling
        if (A.cols() != x_lb.size())
            throw std::invalid_argument("Matrix A must have the same number of columns as the size of the Box");

        // declare
        Box y (A.rows());

        // linear map
        for (int i=0; i<A.rows(); i++)
        {
            y[i] = Interval(0, 0);
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it)
            {
                y[i].add_assign(y[i], ((*this)[it.col()]*it.value()).as_view());
            }
        }
        return y;
    }

    Interval Box::dot(const Eigen::Vector<zono_float, -1>& x) const
    {
        // input handling
        if (x.size() != x_lb.size())
            throw std::invalid_argument("Vector x must have the same size as the Box");

        // declare
        Interval y(0, 0);

        // linear map
        for (int i=0; i<this->x_lb.size(); i++)
            y.add_assign(y, ((*this)[i]*x(i)));
        return y;
    }

    void Box::permute(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P)
    {
        this->x_lb = P*this->x_lb;
        this->x_ub = P*this->x_ub;
    }

    std::string Box::print() const
    {
        std::stringstream ss;
        ss << "Box: " << std::endl;
        for (Eigen::Index i=0; i<x_lb.size(); i++)
        {
            ss << "  " << (*this)[i] << std::endl;
        }
        return ss.str();
    }

    std::ostream& operator<<(std::ostream& os, const Box& box)
    {
        os << box.print();
        return os;
    }

    bool Box::contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const Eigen::Vector<zono_float, -1>& b, const int iter,
                             const std::set<int>& constraints)
    {
        for (int i=0; i<iter; i++)
        {
            // loop through constraints
            for (const int k : constraints)
            {
                // forward propagate
                Interval y(0, 0);
                std::vector<int> cols; // keeping track of columns
                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
                {
                    y.add_assign(y, (*this)[it.col()].to_interval()*it.value());
                    cols.push_back(static_cast<int>(it.col()));
                }

                // check validity
                if (!y.contains(b(k)))
                    return false;
                else
                    y = Interval(b(k), b(k));

                // backward propagate
                for (const int col : cols)
                {
                    Interval x = y; // init
                    zono_float a_col=one;
                    for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
                    {
                        if (it.col() != col)
                        {
                            x.add_assign(x, (*this)[it.col()].to_interval()*(-it.value()));
                        }
                        else
                        {
                            a_col = it.value();
                        }
                    }

                    // update interval
                    if (std::abs(a_col) > zono_eps) // divide by zero protection
                        (*this)[col].intersect_assign((*this)[col], (x * (one/a_col)).as_view());
                }
            }
        }

        return true; // constraints valid
    }

    void Box::get_vars_cons(const Eigen::SparseMatrix<zono_float>& A, const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                           std::set<int>& constraints, std::set<int>& vars, const std::set<int>& new_vars, int depth, int max_depth)
    {
        // immediately copy over constraints and vars
        vars.insert(new_vars.begin(), new_vars.end());

        // find new constraints
        std::set<int> new_constraints;
        for (const int i : new_vars)
        {
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(A, i); it; ++it)
            {
                if (!constraints.count(static_cast<int>(it.row())))
                {
                    new_constraints.insert(static_cast<int>(it.row()));
                }
            }
        }

        // find new vars
        std::set<int> new_new_vars;
        for (const int i : new_constraints)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A_rm, i); it; ++it)
            {
                if (!vars.count(static_cast<int>(it.col())))
                {
                    new_new_vars.insert(static_cast<int>(it.col()));
                }
            }
        }

        // add new constraints to set
        constraints.insert(new_constraints.begin(), new_constraints.end());

        // recurse if able
        depth++;
        if (depth < max_depth && new_new_vars.empty())
            get_vars_cons(A, A_rm, constraints, vars, new_new_vars, depth, max_depth);
    }

}