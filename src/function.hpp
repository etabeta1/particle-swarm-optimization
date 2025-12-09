#ifndef FUNCTION_HPP
#define FUNCTION_HPP

namespace Swarm
{
    // Class declarations

    /*
        Definition of a generic objective function.
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - evaluate: given a point, returns the function value at that point
    */
    template <typename T = float, int dim = 2>
    class Function
    {
    public:
        virtual float evaluate(const Point<T, dim>& p) const = 0;
        virtual ~Function() = default;
    };

    /*
        Definition of the sphere function.
        This function is  continuous, convex, unimodal, differentiable,
        separable, highly symmetric, and rotationally invariant. 
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - evaluate: given a point, returns the sum of squares of its coordinates
    */
    template <typename T = float, int dim = 2>
    class SphereFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += p[i] * p[i]; 
            }
            return sum;
        }
    };

    /*
        Definition of the ellipsoid function.
        This function is continuous, convex, differentiable, separable, and unimodal
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - evaluate: given a point, returns the weighted sum of squares of its coordinates (weighted sphere function).

    */
    template <typename T = float, int dim = 2>
    class EllipsoidFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += i*p[i] * p[i]; 
            }
            return sum;
        }
    };

    /*
    Definition of the quintic function.
    The function is continuos and has two distinct global minima with
     f(x*) = 0 at x* = {−1, −1, . . . , −1} or x* = {2, 2, . . . , 2}
    - `T` is the type used to store the coordinates (defaults to `float`)
    - `dim` is the number of dimensions for the vector (defaults to 2)

    - evaluate: given a point, returns the sum of the modulus of some polynomials of fifth grade in the coordinates of the point.
    */
    template <typename T = float, int dim = 2>
    class QuinticFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += abs(pow(p[i],5) - 3*pow(p[i],4) + 4*pow(p[i],3) + 2*pow(p[i],2) -10*p[i] -4); 
            }
            return sum;
        }
    };

    /*
    Definition of the dropwave function.
    The function is continuos and has a global minima in x* = {0, 0, . . . , 0} with
     f(x*) = -1
    - `T` is the type used to store the coordinates (defaults to `float`)
    - `dim` is the number of dimensions for the vector (defaults to 2)

    - evaluate: given a point, returns the value of the DropWave function in that point.
    */   

    template <typename T = float, int dim = 2>
    class DropwaveFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;
            float part_sum = 0.0f;
            for (int i = 0; i < dim; ++i) {
                part_sum = p[i]*p[i];
            }
            sum = 1 - (1 + cos(12 * sqrt(part_sum))) / (0.5 * part_sum + 2);
            return sum;
        }
    };


    /*
    Definition of the alpine1 function.
    The function is continuos and has a global minima in x* = {0,   0, . . . , 0} with
     f(x*) = 0
    - `T` is the type used to store the coordinates (defaults to `float`)
    - `dim` is the number of dimensions for the vector (defaults to 2)

    - evaluate: given a point, returns the sum of the absolute values of the coordinates multiplied by the sine of the coordinates plus 0.1 times the coordinates.

    */
    template <typename T = float, int dim = 2>
    class Alpine1Function : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;

            for (int i = 0; i < dim; ++i) {
                sum += abs(p[i] * sin(p[i]) + 0.1 * p[i]);
            }
            return sum;
        }
    };

    /*
    Definition of the ackley function.
    The function is continuos and has a global minima in x* = {0,   0, . . . , 0} with
     f(x*) = 0

    - `T` is the type used to store the coordinates (defaults to `float`)
    - `dim` is the number of dimensions for the vector (defaults to 2)

    - evaluate: given a point, returns the value of the Ackley function at that point.
    */
    template <typename T = float, int dim = 2>
    class AckleyFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;
            float part_sum_1 = 0.0f;
            float part_sum_2 = 0.0f;
            for (int i = 0; i < dim; ++i) {
               part_sum_1 += p[i]*p[i];
                part_sum_2 += cos(2 * M_PI * p[i]); 
            }
            sum = -20 * exp(-0.2 * sqrt(part_sum_1 / dim)) - exp(part_sum_2 / dim) + 20 + exp(1);
            return sum;
        }
    };
}

#endif