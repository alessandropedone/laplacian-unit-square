#ifndef MUPARSERX_INTERFACE_HPP
#define MUPARSERX_INTERFACE_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <mpParser.h>

namespace muparser
{
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
            // Assign input values to the parser's variable storage
            for (unsigned i = 0; i < N; ++i)
            {
                M_value.At(i) = x[i];
            }
            mup::Value val;
            double res;
            try
            {
                // Evaluate the parsed expression
                val = M_parser.Eval();
                if (val.IsScalar())
                {
                    // If the result is a scalar, return it directly
                    res = val.GetFloat();
                }
                else
                {
                    std::cerr << "Muparsex error: expression did not evaluate to a scalar." << std::endl;
                    throw mup::ParserError("Expected scalar result.");
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

    private:
        /// @brief A copy of the muparserX expression, used for the copy operations
        string_type My_e;
        /// @brief The muparseX engine
        mup::ParserX M_parser;
        /// @brief The muparserX value used to set the variables in the engine
        mutable mup::Value M_value;
        /// @brief The number of variables in the expression
        mutable unsigned N;
    }; // class muParserXInterface    
} // namespace muparser
#endif // MUPARSERX_INTERFACE_HPP