#ifndef FUNCTION_HPP
#define FUNCTION_HPP

namespace Swarm
{
    /**
     * \brief Abstract base class representing a generic objective function.
     * \tparam T The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2).
     * \see Point
     */
    template <typename T = float, int dim = 2>
    class ObjectiveFunction
    {
    public:
        virtual T evaluate(const Point<T, dim> &p) const = 0;
        virtual ~ObjectiveFunction() = default;
    };

    namespace ObjectiveFunctions
    {
        /**
         * \brief Class representing the Sphere objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The Sphere function is defined as \f[
         * f(x) = \sum_{i=1}^{dim} x_i^2
         * \f] and it is minimized in position \f((0, 0, \dots, 0)\f) with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class SphereFunction : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                sum = p.squareNorm2();
                return sum;
            }
        };

        /**
         * \brief Class representing the Ellipsoid objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The Ellipsoid function is defined as \f[
         * f(x) = \sum_{i=1}^{dim} i \cdot x_i^2
         * \f] and it is minimized in position \f((0, 0, \dots, 0)\f) with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class EllipsoidFunction : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                for (int i = 0; i < dim; ++i)
                {
                    sum += i * p[i] * p[i];
                }
                return sum;
            }
        };

        /**
         * \brief Class representing the Quintic objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The Quintic function is defined as \f[
         * f(x) = \sum_{i=1}^{dim} |x_i^5 - 3x_i^4 + 4x_i^3 + 2x_i^2 - 10x_i - 4|
         * \f] and it is minimized in \f[ x \in \mathbb{X}^d : \mathbb{X} = \{-1, -0.402627\dots, 2\} \f] with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class QuinticFunction : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                sum = (p.pow(5) - p.pow(4) * static_cast<T>(3) + p.pow(3) * static_cast<T>(4) + p.pow(2) * static_cast<T>(2) - p * static_cast<T>(10) - Point<T, dim>(static_cast<T>(4))).abs().norm1();
                return sum;
            }
        };

        /**
         * \brief Class representing the DropWave objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The DropWave function is defined as \f[
         * f(x) = 1 - \frac{1 + \cos(12 \sqrt{\sum_{i=1}^{dim} x_i^2})}{0.5 \cdot \sum_{i=1}^{dim} x_i^2 + 2}
         * \f] and it is minimized in position \f((0, 0, \dots, 0)\f) with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class DropwaveFunction : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                T part_sum = 0.0f;
                part_sum = p.squareNorm2();
                sum = 1 - (1 + cos(12 * sqrt(part_sum))) / (0.5 * part_sum + 2);
                return sum;
            }
        };

        /**
         * \brief Class representing the Alpine1 objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The Alpine1 function is defined as \f[
         * f(x) = \sum_{i=1}^{dim} |x_i \cdot \sin(x_i) + 0.1 \cdot x_i|
         * \f] and it is minimized in position \f((0, 0, \dots, 0)\f) with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class Alpine1Function : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                sum = ((p.sin()) * p + p * 0.1f).abs().norm1();
                return sum;
            }
        };

        /**
         * \brief Class representing the Ackley objective function.
         * \tparam T The data type for the coordinates (default is float).
         * \tparam dim The dimensionality of the points (default is 2).
         * \see ObjectiveFunction
         *
         * The Ackley function is defined as \f[
         * f(x) = -20 \cdot \exp\left(-0.2 \cdot \sqrt{\frac{1}{dim} \sum_{i=1}^{dim} x_i^2}\right) - \exp\left(\frac{1}{dim} \sum_{i=1}^{dim} \cos(2 \pi x_i)\right) + 20 + e
         * \f] and it is minimized in position \f((0, 0, \dots, 0)\f) with value \f(0\f).
         */
        template <typename T = float, int dim = 2>
        class AckleyFunction : public ObjectiveFunction<T, dim>
        {
        public:
            T evaluate(const Point<T, dim> &p) const override
            {
                T sum = 0.0f;
                T part_sum_1 = 0.0f;
                T part_sum_2 = 0.0f;
                part_sum_1 = p.squareNorm2();
                part_sum_2 = (p * static_cast<T>(2 * M_PI)).cos().norm1();
                sum = -20 * exp(-0.2 * sqrt(part_sum_1 / dim)) - exp(part_sum_2 / dim) + 20 + exp(1);
                return sum;
            }
        };
    };
}

#endif