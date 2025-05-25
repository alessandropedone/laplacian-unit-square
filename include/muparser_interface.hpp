#ifndef MUPARSERX_INTERFACE_HPP
#define MUPARSERX_INTERFACE_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <mpParser.h>

namespace muparser
{
    using matrix_type = std::vector<std::vector<double>>;
    using vector_type = std::vector<double>;
    using string_type = std::string;

    /**
     * @brief muParserX interface for mathematical expression parsing and evaluation
     * Adapted from muParserInterface inside pacs-examples repository (https://github.com/pacs-course/pacs-examples.git)
     */
    class muParserXInterface
    {

    public:
        //! Default constructor
        //!
        //! mup::pckALL_NON_COMPLEX|mup::pckMATRIX means that I do not want the module
        //! for complex numbers but I want to treat arrays and matrices in muparserX
        //! expressions
        //! N is not a template parameter and is set at run time in order to achive more flexibility
        muParserXInterface(const unsigned N = 1)
            : My_e(),
              M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX), M_value{N, 0.0}, N(N)
        {
            M_parser.DefineVar("x", mup::Variable(&M_value));
        }
        //! Constructor that takes a string containing muParserX expression
        muParserXInterface(const string_type expression, const unsigned N = 1) : muParserXInterface(N)
        {
            try
            {
                My_e = expression;
                M_parser.SetExpr(My_e.c_str());
            }
            catch (mup::ParserError &error)
            {
                std::cerr << "Muparsex error with code:" << error.GetCode()
                          << std::endl;
                std::cerr << "While processing expression: " << error.GetExpr()
                          << std::endl;
                std::cerr << "Error Message: " << error.GetMsg() << std::endl;
                throw error;
            }
        }

        /*!
         * The copy constructor
         *
         * MuparserX has a particular design, which obliges to define a special copy
         * constructor The reson is that a muparser engine stores the address of the
         * variables. So a normal copy would do a shallow copy, which is NOT what you
         * want. Moreover, because of a poor design, you may loose the expression.
         * That's why I keep a copy in the class as a string and a redefine in in the
         * muparser engine.
         *
         * @param mpi the muParserXIterface to be copied
         */
        muParserXInterface(muParserXInterface const &mpi)
            : My_e(mpi.My_e),
              M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX), M_value{mpi.N, 0.0}, N(mpi.N)
        {
            M_parser.DefineVar("x", mup::Variable(&M_value));
            M_parser.SetExpr(My_e.c_str());
        }

        /*!
         * The copy assignment operator
         *
         * MuparserX has a particular design, which obliges to define a special copy
         * assignement
         * @param mpi the muParserXInterface to be copied
         * The copy constructor
         */
        muParserXInterface
        operator=(muParserXInterface const &mpi)
        {
            if (this != &mpi)
            {
                this->My_e = mpi.My_e;
                this->M_parser.ClearVar(); // clear the variables!
                this->M_value = mpi.M_value;
                this->N = mpi.N;
                M_parser.DefineVar("x", mup::Variable(&M_value));
                M_parser.SetExpr(My_e.c_str());
            }
            return *this;
        }

        //! Sets the muparserX expression.
        /*!
         * Beware, the input variables are indicated by x[].
         * example of a valid expression: sin(x[0])+x[1]*x[2]
         * @par e The expression
         */
        void
        set_expression(const string_type &e)
        {
            My_e = e;
            M_parser.SetExpr(e.c_str());
        }

        //! Function call operator to evaluate the expression
        /*!
         * Takes a vector of input variables, evaluates the expression, and returns the result.
         * Supports both scalar and matrix outputs.
         *
         * @param x Vector of input variable values.
         * @return Vector of evaluated results.
         */
        matrix_type operator()(const vector_type &x) const
        {
            // Assign input values to the parser's variable storage
            for (unsigned i = 0; i < N; ++i)
            {
                M_value.At(i) = x[i];
            }

            mup::Value val;
            matrix_type res;
            try
            {
                // Evaluate the parsed expression
                val = M_parser.Eval();

                if (val.IsScalar())
                {
                    // If the result is a scalar, return a single-element vector
                    res = matrix_type(1, vector_type(1));
                    res[0][0] = val.GetFloat();
                }
                else
                {
                    // If the result is a matrix, extract its values
                    unsigned rows = val.GetRows();
                    unsigned cols = val.GetCols();
                    res = matrix_type(rows, vector_type(cols));

                    for (unsigned i = 0; i < rows; ++i)
                    {
                        for (unsigned j = 0; j < cols; ++j)
                        {
                            res[i][j] = val.At(i, j).GetFloat();
                        }
                    }
                }
            }
            catch (mup::ParserError &error)
            {
                // Handle parsing errors
                std::cerr << "Muparsex error with code:" << error.GetCode() << std::endl;
                std::cerr << "While processing expression: " << error.GetExpr() << std::endl;
                std::cerr << "Error Message: " << error.GetMsg() << std::endl;
                throw error;
            }
            return res;
        }

    protected:
        /// @brief A copy of the muparserX expression, used for the copy operations
        string_type My_e;
        /// @brief The muparseX engine
        mup::ParserX M_parser;
        /// @brief The muparserX value used to set the variables in the engine
        mutable mup::Value M_value;
        /// @brief The number of variables in the expression
        mutable unsigned N;
    };

    /**
     * \brief A muParserX interface with a vector output
     *
     * This class is a specialisation of the muParserXInterface class, with a vector output.
     */
    class muParserXVectorInterface : public muParserXInterface
    {
    public:
        /*!
         * Constructor of the muParserXVectorInterface class
         *
         * This is a thin wrapper around the muParserXInterface constructor.
         * It is used to create a muParserXInterface that returns a vector.
         *
         * @param expression The expression to be parsed
         * @param N The size of the vector returned by the expression
         */
        muParserXVectorInterface(const string_type expression, const unsigned N = 1) : muParserXInterface(expression, N) {}

        /*!
         * Since the expression may evaluate to a matrix, the result is
         * taken to be the first column of the matrix.
         *
         * @param x Vector of input variable values.
         * @return Vector of evaluated results.
         */
        vector_type operator()(const vector_type &x) const
        {
            matrix_type mat = muParserXInterface::operator()(x);
            vector_type result;
            result.reserve(mat.size());
            for (const auto& row : mat) {
                if (!row.empty()) {
                    result.push_back(row[0]);
                }
            }
            return result;
        };
    };

    /**
     * \brief A muParserX interface with a scalar output
     *
     * This class is a specialisation of the muParserXVectorInterface class, with a scalar output.
     *
     * It is used to create a muParserXInterface that returns a scalar.
     */
    class muParserXScalarInterface : public muParserXVectorInterface
    {
    public:
        /*!
         * Constructor of the muParserXScalarInterface class
         *
         * This is a thin wrapper around the muParserXVectorInterface constructor.
         * It is used to create a muParserXInterface that returns a scalar.
         *
         * @param expression The expression to be parsed
         * @param N The size of the vector returned by the expression
         */
        muParserXScalarInterface(const string_type expression, const unsigned N = 1) : muParserXVectorInterface(expression, N) {}

        /*!
         * Evaluate the expression and return the first element of the result.
         *
         * The expression is evaluated using the muParserXInterface::operator() method.
         * The result of the expression is taken to be the first element of the
         * returned vector.
         *
         * @param x Vector of input variable values.
         * @return The first element of the result of the expression.
         */
        double operator()(const vector_type &x) const
        {
            vector_type vec = muParserXVectorInterface::operator()({x});
            return vec.size() == 0 ? 0.0 : vec[0];
        }
    };
}
#endif // MUPARSERX_INTERFACE_HPP