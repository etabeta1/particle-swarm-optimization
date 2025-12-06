#ifndef FUNCTION_HPP
#define FUNCTION_HPP

namespace Swarm
{
    // Class declarations
    template <typename T = float, int dim = 2>
    class Function
    {
    public:
        virtual float evaluate(const Point<T, dim>& p) const = 0;
        virtual ~Function() = default;
    };

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